mod adapters;
mod io;
mod kmers;
mod metrics;


use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::{Arc, Mutex};

use needletail::parse_fastx_file;
use needletail::Sequence;
use rayon::prelude::*;

use crate::types::{
    format_number, FileStats, LogLevel, ProcessConfig, ProcessingStatus, SharedState,
    DUP_SAMPLE_SIZE, FLUSH_INTERVAL, PARALLEL_BATCH,
};

use self::io::{open_writer, write_fastq_record};
use self::kmers::count_kmers_parallel;
use self::metrics::{fingerprint, BatchAccum};

// ---------------------------------------------------------------------------
// Owned record — safe to send across threads
// ---------------------------------------------------------------------------

struct OwnedRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn process_files(paths: Vec<String>, state: Arc<Mutex<SharedState>>, config: ProcessConfig) {
    let n = paths.len();
    for (idx, path) in paths.into_iter().enumerate() {
        {
            let mut s = state.lock().unwrap();
            s.current_file_idx = idx;
            let short = Path::new(&path)
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or(&path)
                .to_string();
            s.log(LogLevel::Info, format!("Starting file {}/{}: {}", idx + 1, n, short));
        }
        process_single_file(path, Arc::clone(&state), &config);
        if matches!(state.lock().unwrap().status, ProcessingStatus::Error(_)) {
            return;
        }
    }

    let mut s = state.lock().unwrap();
    s.status = ProcessingStatus::Completed;
    s.current = None;
    let total: u64 = s.completed_files.iter().map(|f| f.read_count).sum();
    s.log(
        LogLevel::Success,
        format!("All {} file(s) complete — {} total reads", n, format_number(total)),
    );
}

// ---------------------------------------------------------------------------
// Per-file processing — parallel batch model
// ---------------------------------------------------------------------------

fn process_single_file(file_path: String, state: Arc<Mutex<SharedState>>, config: &ProcessConfig) {
    let file_size = std::fs::metadata(&file_path).map(|m| m.len()).unwrap_or(0);

    let trim_path = config.trim_output.then(|| {
        let stem = Path::new(&file_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output")
            .trim_end_matches(".fastq")
            .to_string();
        format!("{}/{}_trimmed.fastq.gz", config.output_dir, stem)
    });

    {
        let mut s = state.lock().unwrap();
        let mut stats = FileStats::new(file_path.clone(), file_size);
        stats.trim_output_path = trim_path.clone();
        s.current = Some(stats);
        s.log(
            LogLevel::Info,
            format!(
                "Opening {} ({:.2} GB) [{} threads]{}",
                file_path,
                file_size as f64 / 1_073_741_824.0,
                rayon::current_num_threads(),
                trim_path.as_deref()
                    .map(|p| format!("  trim→ {}", p))
                    .unwrap_or_default(),
            ),
        );
    }

    let mut reader = match parse_fastx_file(&file_path) {
        Ok(r) => r,
        Err(e) => {
            let msg = format!("Cannot open {}: {}", file_path, e);
            let mut s = state.lock().unwrap();
            s.log(LogLevel::Error, msg.clone());
            s.status = ProcessingStatus::Error(msg);
            return;
        }
    };

    let mut trim_writer = if let Some(ref tp) = trim_path {
        match open_writer(tp) {
            Ok(w) => Some(w),
            Err(e) => {
                let msg = format!("Cannot open trim output {}: {}", tp, e);
                let mut s = state.lock().unwrap();
                s.log(LogLevel::Error, msg.clone());
                s.status = ProcessingStatus::Error(msg);
                return;
            }
        }
    } else {
        None
    };

    // --- Running totals (merged from parallel batches) ---
    let mut read_count = 0u64;
    let mut total_bases = 0u64;
    let mut gc_count = 0u64;
    let mut quality_sum = 0u64;
    let mut quality_bases = 0u64;
    let mut q20_bases = 0u64;
    let mut q30_bases = 0u64;
    let mut adapter_hits = 0u64;
    let mut trimmed_reads = 0u64;
    let mut trimmed_bases_removed = 0u64;
    let mut min_length = u64::MAX;
    let mut max_length = 0u64;
    let mut bytes = 0u64;
    let mut length_histogram: HashMap<u64, u64> = HashMap::with_capacity(512);
    let mut quality_by_pos = vec![(0u64, 0u64); crate::types::MAX_QUAL_POSITION];
    let mut per_tile: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut kmer_total: HashMap<[u8; 4], u64> = HashMap::with_capacity(256);
    let mut base_composition = vec![[0u64; 5]; crate::types::MAX_QUAL_POSITION];
    let mut quality_distribution = vec![0u64; 43];
    let mut dup_seen: HashSet<u64> = HashSet::with_capacity(DUP_SAMPLE_SIZE);
    let mut dup_sampled = 0u64;
    let mut flush_counter = 0u64;
    let mut total_flushed = 0u64;
    let mut error_occurred = false;

    // --- Batch reading loop ---
    loop {
        // Read PARALLEL_BATCH records into owned memory
        let mut batch: Vec<OwnedRecord> = Vec::with_capacity(PARALLEL_BATCH);
        let mut batch_bytes = 0u64;

        while batch.len() < PARALLEL_BATCH {
            match reader.next() {
                None => break,
                Some(Err(e)) => {
                    let mut s = state.lock().unwrap();
                    if config.strict {
                        let msg = format!("Bad record: {}", e);
                        s.log(LogLevel::Error, msg.clone());
                        s.status = ProcessingStatus::Error(msg);
                        error_occurred = true;
                        break;
                    }
                    s.log(LogLevel::Warning, format!("Skipping bad record: {}", e));
                }
                Some(Ok(record)) => {
                    let id = record.id().to_vec();
                    let seq = record.normalize(false).to_vec();
                    let qual = record.qual().map(|q| q.to_vec());
                    batch_bytes += seq.len() as u64 + id.len() as u64 + 10;
                    if let Some(ref q) = qual { batch_bytes += q.len() as u64; }
                    batch.push(OwnedRecord { id, seq, qual });
                }
            }
        }

        if batch.is_empty() || error_occurred {
            break;
        }

        bytes += batch_bytes;

        // --- Parallel stats computation ---
        let cfg = config;
        let acc: BatchAccum = batch
            .par_iter()
            .fold(BatchAccum::default, |mut acc, rec| {
                // Quality trimming
                let (seq, qual): (&[u8], Option<&[u8]>) =
                    if cfg.quality_trim_threshold > 0 {
                        if let Some(ref q) = rec.qual {
                            let cut = adapters::quality_trim_3p(q, cfg.quality_trim_threshold);
                            (&rec.seq[..cut], Some(&q[..cut]))
                        } else {
                            (&rec.seq, rec.qual.as_deref())
                        }
                    } else {
                        (&rec.seq, rec.qual.as_deref())
                    };

                // Adapter detection
                let adapter_pos = adapters::find_adapter_pos_with_custom(seq, &cfg.custom_adapters);
                let (effective_seq, effective_qual, adapter_hit) = match adapter_pos {
                    Some(pos) if cfg.trim_output => {
                        (&seq[..pos], qual.map(|q| &q[..pos]), true)
                    }
                    Some(_) => (seq, qual, true),
                    None => (seq, qual, false),
                };

                // Per-tile
                let tile_info = parse_illumina_tile(&rec.id).and_then(|tile_id| {
                    effective_qual.map(|q| {
                        let phred_sum: u64 = q.iter().map(|&b| b.saturating_sub(33) as u64).sum();
                        (tile_id, phred_sum, q.len() as u64)
                    })
                });

                let fp = fingerprint(effective_seq);
                acc.process(effective_seq, effective_qual, tile_info, fp, adapter_hit);
                acc
            })
            .reduce(BatchAccum::default, BatchAccum::merge);

        // --- Duplication sampling (sequential, only first DUP_SAMPLE_SIZE) ---
        for fp in &acc.fingerprints {
            if dup_sampled >= DUP_SAMPLE_SIZE as u64 { break; }
            dup_seen.insert(*fp);
            dup_sampled += 1;
        }

        // --- K-mer counting for this batch ---
        if !acc.kmer_seqs.is_empty() {
            let kmers = count_kmers_parallel(&acc.kmer_seqs);
            for (k, v) in kmers {
                *kmer_total.entry(k).or_insert(0) += v;
            }
        }

        // --- Trim output (sequential to preserve order) ---
        if config.trim_output {
            for rec in &batch {
                let (seq, qual): (&[u8], Option<&[u8]>) =
                    if config.quality_trim_threshold > 0 {
                        if let Some(ref q) = rec.qual {
                            let cut = adapters::quality_trim_3p(q, config.quality_trim_threshold);
                            (&rec.seq[..cut], Some(&q[..cut]))
                        } else {
                            (&rec.seq, rec.qual.as_deref())
                        }
                    } else {
                        (&rec.seq, rec.qual.as_deref())
                    };

                let adapter_pos = adapters::find_adapter_pos_with_custom(seq, &config.custom_adapters);
                if let Some(pos) = adapter_pos {
                    let trimmed_seq = &seq[..pos];
                    let trimmed_len = trimmed_seq.len() as u64;
                    trimmed_bases_removed += seq.len() as u64 - trimmed_len;
                    if trimmed_len >= config.min_length {
                        trimmed_reads += 1;
                        if let Some(ref mut w) = trim_writer {
                            let tq = qual.map(|q| &q[..pos]);
                            if let Err(e) = write_fastq_record(w, &rec.id, trimmed_seq, tq) {
                                state.lock().unwrap().log(LogLevel::Warning, format!("Trim write error: {}", e));
                            }
                        }
                    }
                } else if let Some(ref mut w) = trim_writer {
                    if let Err(e) = write_fastq_record(w, &rec.id, seq, qual) {
                        state.lock().unwrap().log(LogLevel::Warning, format!("Trim write error: {}", e));
                    }
                }
            }
        }

        // --- Merge batch accumulator into running totals ---
        read_count    += acc.read_count;
        total_bases   += acc.total_bases;
        gc_count      += acc.gc_count;
        quality_sum   += acc.quality_sum;
        quality_bases += acc.quality_bases;
        q20_bases     += acc.q20_bases;
        q30_bases     += acc.q30_bases;
        adapter_hits  += acc.adapter_hits;
        if acc.min_length < min_length { min_length = acc.min_length; }
        if acc.max_length > max_length { max_length = acc.max_length; }
        for (k, v) in acc.length_histogram {
            *length_histogram.entry(k).or_insert(0) += v;
        }
        for i in 0..crate::types::MAX_QUAL_POSITION {
            quality_by_pos[i].0    += acc.quality_by_pos[i].0;
            quality_by_pos[i].1    += acc.quality_by_pos[i].1;
            for j in 0..5 {
                base_composition[i][j] += acc.base_composition[i][j];
            }
        }
        for i in 0..43 {
            quality_distribution[i] += acc.quality_distribution[i];
        }
        for (k, v) in acc.per_tile {
            let e = per_tile.entry(k).or_insert((0, 0));
            e.0 += v.0;
            e.1 += v.1;
        }

        // --- Periodic UI flush ---
        flush_counter += acc.read_count;
        if flush_counter >= FLUSH_INTERVAL {
            flush_counter = 0;
            total_flushed += FLUSH_INTERVAL;

            let mut s = state.lock().unwrap();
            if let Some(ref mut cur) = s.current {
                cur.read_count    = read_count;
                cur.total_bases   = total_bases;
                cur.gc_count      = gc_count;
                cur.quality_sum   = quality_sum;
                cur.quality_bases = quality_bases;
                cur.q20_bases     = q20_bases;
                cur.q30_bases     = q30_bases;
                cur.adapter_hits  = adapter_hits;
                cur.trimmed_reads = trimmed_reads;
                cur.trimmed_bases_removed = trimmed_bases_removed;
                cur.min_length    = min_length;
                cur.max_length    = max_length;
                cur.bytes_processed = bytes.min(file_size);
            }

            if total_flushed % 500_000 == 0 {
                let elapsed = s.elapsed_secs();
                let rps  = if elapsed > 0.0 { read_count as f64 / elapsed } else { 0.0 };
                let mbps = if elapsed > 0.0 { bytes as f64 / 1_048_576.0 / elapsed } else { 0.0 };
                s.log(LogLevel::Info, format!(
                    "{} reads  {:.0} reads/s  {:.1} MB/s",
                    format_number(read_count), rps, mbps,
                ));
            }
        }
    }

    drop(trim_writer);

    let dup_rate = if dup_sampled > 0 {
        (1.0 - dup_seen.len() as f64 / dup_sampled as f64) * 100.0
    } else { 0.0 };

    let mut s = state.lock().unwrap();
    if let Some(ref mut cur) = s.current {
        cur.read_count            = read_count;
        cur.total_bases           = total_bases;
        cur.gc_count              = gc_count;
        cur.quality_sum           = quality_sum;
        cur.quality_bases         = quality_bases;
        cur.q20_bases             = q20_bases;
        cur.q30_bases             = q30_bases;
        cur.adapter_hits          = adapter_hits;
        cur.trimmed_reads         = trimmed_reads;
        cur.trimmed_bases_removed = trimmed_bases_removed;
        cur.min_length            = min_length;
        cur.max_length            = max_length;
        cur.bytes_processed       = file_size;
        cur.length_histogram      = length_histogram;
        cur.quality_by_position   = quality_by_pos;
        cur.kmer_counts           = kmer_total;
        cur.dup_rate_pct          = dup_rate;
        cur.per_tile_quality      = per_tile;
        cur.base_composition      = base_composition;
        cur.quality_distribution  = quality_distribution;
    }

    let elapsed = s.elapsed_secs();
    let mbps = if elapsed > 0.0 { bytes as f64 / 1_048_576.0 / elapsed } else { 0.0 };
    let gc_pct = if total_bases > 0 { gc_count as f64 / total_bases as f64 * 100.0 } else { 0.0 };
    let avg_q  = if quality_bases > 0 { quality_sum as f64 / quality_bases as f64 } else { 0.0 };
    let adapter_pct = if read_count > 0 { adapter_hits as f64 / read_count as f64 * 100.0 } else { 0.0 };
    let trim_note = if config.trim_output {
        format!(" | trimmed {:.1}%", if read_count > 0 { trimmed_reads as f64 / read_count as f64 * 100.0 } else { 0.0 })
    } else { String::new() };
    let tile_count = s.current.as_ref().map(|c| c.per_tile_quality.len()).unwrap_or(0);
    let tile_note = if tile_count > 0 { format!(" | {} tiles", tile_count) } else { String::new() };

    s.log(LogLevel::Success, format!(
        "Done: {} reads | GC {:.1}% | Q{:.1} | {:.2}% adapter | dup ~{:.1}% | {:.1} MB/s{}{}",
        format_number(read_count), gc_pct, avg_q, adapter_pct, dup_rate, mbps, trim_note, tile_note,
    ));

    if let Some(done) = s.current.take() {
        s.completed_files.push(done);
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn parse_illumina_tile(id: &[u8]) -> Option<u32> {
    let s = std::str::from_utf8(id).ok()?;
    let s = s.split_ascii_whitespace().next()?;
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() >= 6 { parts[4].parse().ok() } else { None }
}
