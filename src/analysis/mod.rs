mod adapters;
mod io;
mod kmers;
mod metrics;


use std::collections::HashMap;
use std::path::Path;
use std::sync::{Arc, Mutex};

use hyperloglog::HyperLogLog as Hll;
use needletail::parse_fastx_file;
use needletail::Sequence;
use rayon::prelude::*;

use crate::types::{
    format_number, FileStats, LogLevel, ProcessConfig, ProcessingStatus, SharedState,
    FLUSH_INTERVAL, PARALLEL_BATCH, PHRED_BUCKETS,
};

use self::io::{open_writer, write_fastq_record};
use self::kmers::count_kmers_parallel;
use self::metrics::{fingerprint, BatchAccum};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Apply quality trim + adapter clip to a record, returning effective slices.
#[inline]
fn trim_record<'a>(
    seq: &'a [u8],
    qual: Option<&'a [u8]>,
    config: &ProcessConfig,
) -> (&'a [u8], Option<&'a [u8]>, bool) {
    // 3′ quality trim
    let (seq, qual) = if config.quality_trim_threshold > 0 {
        if let Some(q) = qual {
            let cut = adapters::quality_trim_3p(q, config.quality_trim_threshold);
            (&seq[..cut], Some(&q[..cut]))
        } else {
            (seq, qual)
        }
    } else {
        (seq, qual)
    };

    // Adapter detection
    match adapters::find_adapter_pos_with_custom(seq, &config.custom_adapters) {
        Some(pos) => (&seq[..pos], qual.map(|q| &q[..pos.min(q.len())]), true),
        None => (seq, qual, false),
    }
}

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

    // Paired-end mode: single R1 + R2
    if let Some(ref r2_path) = config.paired_end_r2 {
        assert_eq!(n, 1, "paired-end mode requires exactly one R1 path");
        {
            let mut s = state.lock().unwrap();
            s.log(LogLevel::Info, format!(
                "Paired-end mode: R1={} R2={}",
                paths[0], r2_path
            ));
        }
        process_paired_files(paths[0].clone(), r2_path.clone(), Arc::clone(&state), &config);
        let mut s = state.lock().unwrap();
        if !matches!(s.status, ProcessingStatus::Error(_)) {
            s.status = ProcessingStatus::Completed;
            s.current = None;
            let total: u64 = s.completed_files.iter().map(|f| f.read_count).sum();
            s.log(LogLevel::Success, format!("Paired-end complete — {} read pairs", format_number(total / 2)));
        }
        return;
    }

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
    let mut quality_distribution = vec![0u64; PHRED_BUCKETS];
    // HyperLogLog: ~1% error rate, counts all reads, ~10 KB RAM
    let mut hll = Hll::new(0.01);
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
                let (effective_seq, effective_qual, adapter_hit) =
                    trim_record(&rec.seq, rec.qual.as_deref(), cfg);

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

        // --- Duplication estimation via HyperLogLog (all reads) ---
        for fp in &acc.fingerprints {
            hll.insert(fp);
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
                let (trimmed_seq, trimmed_qual, adapter_hit) =
                    trim_record(&rec.seq, rec.qual.as_deref(), config);
                if adapter_hit {
                    trimmed_bases_removed += rec.seq.len() as u64 - trimmed_seq.len() as u64;
                }
                let keep = trimmed_seq.len() as u64 >= config.min_length;
                if keep {
                    if adapter_hit { trimmed_reads += 1; }
                    if let Some(ref mut w) = trim_writer {
                        if let Err(e) = write_fastq_record(w, &rec.id, trimmed_seq, trimmed_qual) {
                            state.lock().unwrap().log(LogLevel::Warning, format!("Trim write error: {}", e));
                        }
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
        for i in 0..PHRED_BUCKETS {
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

    let dup_rate = if read_count > 0 {
        let unique = hll.len();
        (1.0 - (unique / read_count as f64).min(1.0)) * 100.0
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

// ---------------------------------------------------------------------------
// Paired-end processing — R1 and R2 processed in lockstep
// ---------------------------------------------------------------------------

fn process_paired_files(
    r1_path: String,
    r2_path: String,
    state: Arc<Mutex<SharedState>>,
    config: &ProcessConfig,
) {
    let r1_size = std::fs::metadata(&r1_path).map(|m| m.len()).unwrap_or(0);
    let r2_size = std::fs::metadata(&r2_path).map(|m| m.len()).unwrap_or(0);

    {
        let mut s = state.lock().unwrap();
        let mut stats = FileStats::new(r1_path.clone(), r1_size + r2_size);
        stats.trim_output_path = config.trim_output.then(|| {
            let stem = Path::new(&r1_path)
                .file_stem().and_then(|s| s.to_str())
                .unwrap_or("r1")
                .trim_end_matches(".fastq").to_string();
            format!("{}/{}_trimmed.fastq.gz", config.output_dir, stem)
        });
        s.current = Some(stats);
    }

    let mut r1_reader = match parse_fastx_file(&r1_path) {
        Ok(r) => r,
        Err(e) => {
            let msg = format!("Cannot open R1 {}: {}", r1_path, e);
            let mut s = state.lock().unwrap();
            s.log(LogLevel::Error, msg.clone());
            s.status = ProcessingStatus::Error(msg);
            return;
        }
    };

    let mut r2_reader = match parse_fastx_file(&r2_path) {
        Ok(r) => r,
        Err(e) => {
            let msg = format!("Cannot open R2 {}: {}", r2_path, e);
            let mut s = state.lock().unwrap();
            s.log(LogLevel::Error, msg.clone());
            s.status = ProcessingStatus::Error(msg);
            return;
        }
    };

    // Running totals — same fields as single-file processing
    let mut read_count = 0u64;
    let mut total_bases = 0u64;
    let mut gc_count = 0u64;
    let mut quality_sum = 0u64;
    let mut quality_bases = 0u64;
    let mut q20_bases = 0u64;
    let mut q30_bases = 0u64;
    let mut adapter_hits = 0u64;
    let mut min_length = u64::MAX;
    let mut max_length = 0u64;
    let mut bytes = 0u64;
    let mut length_histogram: HashMap<u64, u64> = HashMap::with_capacity(512);
    let mut quality_by_pos = vec![(0u64, 0u64); crate::types::MAX_QUAL_POSITION];
    let mut base_composition = vec![[0u64; 5]; crate::types::MAX_QUAL_POSITION];
    let mut quality_distribution = vec![0u64; PHRED_BUCKETS];
    let mut per_tile: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut kmer_total: HashMap<[u8; 4], u64> = HashMap::with_capacity(256);
    let mut hll = Hll::new(0.01);
    let mut flush_counter = 0u64;
    let mut total_flushed = 0u64;

    // Read both files in lockstep batch by batch
    loop {
        let mut batch: Vec<(OwnedRecord, OwnedRecord)> = Vec::with_capacity(PARALLEL_BATCH);
        let mut batch_bytes = 0u64;

        loop {
            if batch.len() >= PARALLEL_BATCH { break; }
            let r1 = match r1_reader.next() {
                None => break,
                Some(Err(e)) => {
                    state.lock().unwrap().log(LogLevel::Warning, format!("R1 bad record: {}", e));
                    continue;
                }
                Some(Ok(rec)) => {
                    let id  = rec.id().to_vec();
                    let seq = rec.normalize(false).to_vec();
                    let qual = rec.qual().map(|q| q.to_vec());
                    batch_bytes += seq.len() as u64 + id.len() as u64;
                    OwnedRecord { id, seq, qual }
                }
            };
            let r2 = match r2_reader.next() {
                None => break, // R2 ended early — truncated pair
                Some(Err(e)) => {
                    state.lock().unwrap().log(LogLevel::Warning, format!("R2 bad record: {}", e));
                    continue;
                }
                Some(Ok(rec)) => {
                    let id  = rec.id().to_vec();
                    let seq = rec.normalize(false).to_vec();
                    let qual = rec.qual().map(|q| q.to_vec());
                    batch_bytes += seq.len() as u64 + id.len() as u64;
                    OwnedRecord { id, seq, qual }
                }
            };
            batch.push((r1, r2));
        }

        if batch.is_empty() { break; }
        bytes += batch_bytes;

        // Process R1+R2 together in parallel — each pair is one unit
        let acc: BatchAccum = batch
            .par_iter()
            .fold(BatchAccum::default, |mut acc, (r1, r2)| {
                for rec in [r1, r2] {
                    let (eff_seq, eff_qual, adapter_hit) =
                        trim_record(&rec.seq, rec.qual.as_deref(), config);
                    let tile_info = parse_illumina_tile(&rec.id).and_then(|tile_id| {
                        eff_qual.map(|q| {
                            let phred_sum: u64 = q.iter().map(|&b| b.saturating_sub(33) as u64).sum();
                            (tile_id, phred_sum, q.len() as u64)
                        })
                    });
                    let fp = fingerprint(eff_seq);
                    acc.process(eff_seq, eff_qual, tile_info, fp, adapter_hit);
                }
                acc
            })
            .reduce(BatchAccum::default, BatchAccum::merge);

        // Merge accumulators
        for fp in &acc.fingerprints { hll.insert(fp); }
        if !acc.kmer_seqs.is_empty() {
            let kmers = count_kmers_parallel(&acc.kmer_seqs);
            for (k, v) in kmers { *kmer_total.entry(k).or_insert(0) += v; }
        }

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
        for (k, v) in acc.length_histogram { *length_histogram.entry(k).or_insert(0) += v; }
        for i in 0..crate::types::MAX_QUAL_POSITION {
            quality_by_pos[i].0 += acc.quality_by_pos[i].0;
            quality_by_pos[i].1 += acc.quality_by_pos[i].1;
            for j in 0..5 { base_composition[i][j] += acc.base_composition[i][j]; }
        }
        for i in 0..PHRED_BUCKETS { quality_distribution[i] += acc.quality_distribution[i]; }
        for (k, v) in acc.per_tile { let e = per_tile.entry(k).or_insert((0,0)); e.0+=v.0; e.1+=v.1; }

        flush_counter += acc.read_count;
        if flush_counter >= FLUSH_INTERVAL {
            flush_counter = 0;
            total_flushed += FLUSH_INTERVAL;
            let mut s = state.lock().unwrap();
            if let Some(ref mut cur) = s.current {
                cur.read_count = read_count;
                cur.total_bases = total_bases;
                cur.gc_count = gc_count;
                cur.quality_sum = quality_sum;
                cur.quality_bases = quality_bases;
                cur.q20_bases = q20_bases;
                cur.q30_bases = q30_bases;
                cur.adapter_hits = adapter_hits;
                cur.min_length = min_length;
                cur.max_length = max_length;
                cur.bytes_processed = bytes.min(r1_size + r2_size);
            }
            if total_flushed % 500_000 == 0 {
                let elapsed = s.elapsed_secs();
                let rps = if elapsed > 0.0 { read_count as f64 / elapsed } else { 0.0 };
                s.log(LogLevel::Info, format!("{} reads  {:.0} pairs/s", format_number(read_count), rps / 2.0));
            }
        }
    }

    let dup_rate = if read_count > 0 {
        (1.0 - (hll.len() / read_count as f64).min(1.0)) * 100.0
    } else { 0.0 };

    let mut s = state.lock().unwrap();
    if let Some(ref mut cur) = s.current {
        cur.read_count = read_count;
        cur.total_bases = total_bases;
        cur.gc_count = gc_count;
        cur.quality_sum = quality_sum;
        cur.quality_bases = quality_bases;
        cur.q20_bases = q20_bases;
        cur.q30_bases = q30_bases;
        cur.adapter_hits = adapter_hits;
        cur.min_length = min_length;
        cur.max_length = max_length;
        cur.bytes_processed = r1_size + r2_size;
        cur.length_histogram = length_histogram;
        cur.quality_by_position = quality_by_pos;
        cur.kmer_counts = kmer_total;
        cur.dup_rate_pct = dup_rate;
        cur.per_tile_quality = per_tile;
        cur.base_composition = base_composition;
        cur.quality_distribution = quality_distribution;
    }

    let elapsed = s.elapsed_secs();
    let pairs = read_count / 2;
    let mbps = if elapsed > 0.0 { bytes as f64 / 1_048_576.0 / elapsed } else { 0.0 };
    let gc_pct = if total_bases > 0 { gc_count as f64 / total_bases as f64 * 100.0 } else { 0.0 };
    let avg_q = if quality_bases > 0 { quality_sum as f64 / quality_bases as f64 } else { 0.0 };

    s.log(LogLevel::Success, format!(
        "Done: {} pairs | GC {:.1}% | Q{:.1} | dup ~{:.1}% | {:.1} MB/s",
        format_number(pairs), gc_pct, avg_q, dup_rate, mbps,
    ));

    if let Some(done) = s.current.take() {
        s.completed_files.push(done);
    }
}
