mod adapters;
mod io;
mod kmers;
mod metrics;
mod mmap_reader;


use std::collections::HashMap;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread;

use needletail::parse_fastx_file;
use needletail::Sequence;
use rayon::prelude::*;

use crate::types::{
    format_number, FileStats, LogLevel, OverrepSeq, ProcessConfig, ProcessingStatus,
    QcStatus, SharedState, ADAPTER_MATCH_LEN, ADAPTERS, FLUSH_INTERVAL, OVERREP_SAMPLE,
    PARALLEL_BATCH, PHRED_BUCKETS,
};

use self::io::{open_writer, write_fastq_record};
use self::kmers::count_kmers_parallel;
use self::metrics::{fingerprint, merge_batch_into_totals, BatchAccum};
use self::mmap_reader::RecordRange;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Strip FASTQ/FASTA extensions from a file path to produce a clean output stem.
/// Handles: .fastq, .fq, .fastq.gz, .fq.gz, .fasta, .fa, .fasta.gz, .fa.gz
fn output_stem(path: &str) -> String {
    let p = Path::new(path);
    // First strip .gz if present
    let without_gz = if path.ends_with(".gz") {
        p.file_stem().and_then(|s| s.to_str()).unwrap_or(path)
    } else {
        p.file_name().and_then(|s| s.to_str()).unwrap_or(path)
    };
    // Then strip the biological extension
    for ext in &[".fastq", ".fq", ".fasta", ".fa"] {
        if let Some(s) = without_gz.strip_suffix(ext) {
            return s.to_string();
        }
    }
    without_gz.to_string()
}

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

    // Adapter detection (skip when --no-adapter is set)
    if !config.flags.adapter_detection {
        return (seq, qual, false);
    }
    match adapters::find_adapter_pos_with_custom(seq, &config.custom_adapters) {
        Some(pos) => (&seq[..pos], qual.map(|q| &q[..pos.min(q.len())]), true),
        None => (seq, qual, false),
    }
}

// ---------------------------------------------------------------------------
// Overrepresented sequence + module status helpers
// ---------------------------------------------------------------------------

/// Check a sequence against known adapter prefixes and return the source name.
fn detect_adapter_source(seq: &[u8]) -> String {
    for (name, adapter) in ADAPTERS {
        let n = adapter.len().min(ADAPTER_MATCH_LEN);
        if seq.windows(n).any(|w| w == &adapter[..n]) {
            return name.to_string();
        }
    }
    "No hit".to_string()
}

/// Finalize overrepresented sequence analysis from accumulated counts.
fn finish_overrep(
    map: HashMap<[u8; 50], u64>,
    sampled: usize,
) -> Vec<OverrepSeq> {
    if sampled == 0 { return Vec::new(); }
    let threshold = sampled as f64 * 0.001; // 0.1% of sampled reads
    let mut sorted: Vec<([u8; 50], u64)> = map
        .into_iter()
        .filter(|(_, c)| *c as f64 >= threshold)
        .collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1));
    sorted.truncate(20);
    sorted
        .into_iter()
        .map(|(key, count)| {
            // Trim trailing null bytes from the fixed-size key
            let seq = &key[..key.iter().rposition(|&b| b != 0).map(|p| p + 1).unwrap_or(0)];
            let percentage = count as f64 / sampled as f64 * 100.0;
            let possible_source = detect_adapter_source(seq);
            OverrepSeq {
                sequence: String::from_utf8_lossy(seq).into_owned(),
                count,
                percentage,
                possible_source,
            }
        })
        .collect()
}

/// Compute FastQC-style per-module pass/warn/fail from finished FileStats.
fn compute_module_status(f: &FileStats) -> Vec<(String, QcStatus)> {
    let mut m = Vec::with_capacity(8);

    // Per-base sequence quality — worst position average
    let qual_per_pos = f.avg_qual_per_position();
    let min_pos_q = qual_per_pos.iter().cloned().fold(f64::MAX, f64::min);
    m.push(("Per-base quality".into(), if min_pos_q >= 28.0 { QcStatus::Pass }
        else if min_pos_q >= 20.0 { QcStatus::Warn } else { QcStatus::Fail }));

    // Per-sequence quality (Q30 %)
    let q30 = f.q30_pct();
    m.push(("Per-sequence quality".into(), if q30 >= 80.0 { QcStatus::Pass }
        else if q30 >= 60.0 { QcStatus::Warn } else { QcStatus::Fail }));

    // Per-base sequence content — check A/T balance and C/G balance
    let bcomp = f.base_composition_pct();
    let unbalanced = bcomp.iter().any(|p| {
        (p[0] - p[3]).abs() > 10.0 || (p[1] - p[2]).abs() > 10.0  // |A-T| or |C-G| > 10%
    });
    let very_unbalanced = bcomp.iter().any(|p| {
        (p[0] - p[3]).abs() > 20.0 || (p[1] - p[2]).abs() > 20.0
    });
    m.push(("Per-base sequence content".into(), if very_unbalanced { QcStatus::Fail }
        else if unbalanced { QcStatus::Warn } else { QcStatus::Pass }));

    // GC content — thresholds are deliberately wide: many organisms (bacteria,
    // fungi, AT-rich parasites) fall well outside the mammalian 40–60% range.
    // We only flag extreme values that suggest severe contamination or a
    // sequencing artifact rather than simply a non-model organism.
    let gc = f.gc_content();
    m.push(("GC content".into(), if (20.0..=80.0).contains(&gc) { QcStatus::Pass }
        else if (10.0..=90.0).contains(&gc) { QcStatus::Warn } else { QcStatus::Fail }));

    // N content — max N% at any position
    let max_n_pct = bcomp.iter().map(|p| p[4]).fold(0.0_f64, f64::max);
    m.push(("N content".into(), if max_n_pct < 5.0 { QcStatus::Pass }
        else if max_n_pct < 20.0 { QcStatus::Warn } else { QcStatus::Fail }));

    // Sequence length distribution
    let len_range = f.max_length.saturating_sub(f.effective_min_length());
    m.push(("Sequence length".into(), if len_range == 0 { QcStatus::Pass }
        else if len_range < 50 { QcStatus::Warn } else { QcStatus::Pass }));

    // Sequence duplication level
    m.push(("Duplication level".into(), if f.dup_rate_pct < 20.0 { QcStatus::Pass }
        else if f.dup_rate_pct < 50.0 { QcStatus::Warn } else { QcStatus::Fail }));

    // Overrepresented sequences
    let has_high_overrep = f.overrepresented_sequences.iter().any(|s| s.percentage >= 1.0);
    m.push(("Overrepresented seqs".into(),
        if f.overrepresented_sequences.is_empty() { QcStatus::Pass }
        else if has_high_overrep { QcStatus::Fail } else { QcStatus::Warn }));

    // Adapter content
    m.push(("Adapter content".into(), if f.adapter_pct() < 5.0 { QcStatus::Pass }
        else if f.adapter_pct() < 20.0 { QcStatus::Warn } else { QcStatus::Fail }));

    m
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
        format!("{}/{}_trimmed.fastq.gz", config.output_dir, output_stem(&file_path))
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

    // Use mmap zero-copy reader for plain FASTQ; fall back to needletail for gzip/fasta
    let is_plain_fastq = !file_path.ends_with(".gz")
        && !file_path.ends_with(".fasta")
        && !file_path.ends_with(".fa");

    if is_plain_fastq {
        process_single_file_mmap(file_path, state, config, trim_path, file_size);
        return;
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
    // Deterministic dedup: exact count on first OVERREP_SAMPLE fingerprints
    let mut dup_map: HashMap<u64, u32> = HashMap::with_capacity(OVERREP_SAMPLE);
    let mut dup_sampled = 0usize;
    // Overrepresented sequences: first 50bp of first OVERREP_SAMPLE reads
    let mut overrep_map: HashMap<[u8; 50], u64> = HashMap::with_capacity(4096);
    let mut overrep_sampled = 0usize;
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
        let per_tile_enabled = cfg.flags.per_tile_quality;
        let acc: BatchAccum = batch
            .par_iter()
            .fold(BatchAccum::default, move |mut acc, rec| {
                let (effective_seq, effective_qual, adapter_hit) =
                    trim_record(&rec.seq, rec.qual.as_deref(), cfg);

                // Per-tile (skip when --no-per-tile or --fast)
                let tile_info = if per_tile_enabled {
                    parse_illumina_tile(&rec.id).and_then(|tile_id| {
                        effective_qual.map(|q| {
                            let phred_sum: u64 = q.iter().map(|&b| b.saturating_sub(33) as u64).sum();
                            (tile_id, phred_sum, q.len() as u64)
                        })
                    })
                } else {
                    None
                };

                let fp = fingerprint(effective_seq);
                acc.process(effective_seq, effective_qual, tile_info, fp, adapter_hit);
                acc
            })
            .reduce(BatchAccum::default, BatchAccum::merge);

        // --- Deterministic dedup (skip when --no-duplication or --fast) ---
        if config.flags.duplication_check {
            for &fp in &acc.fingerprints {
                if dup_sampled < OVERREP_SAMPLE {
                    *dup_map.entry(fp).or_insert(0) += 1;
                    dup_sampled += 1;
                }
            }
        }

        // --- Overrepresented sequences + K-mer counting ---
        if !acc.kmer_seqs.is_empty() {
            if config.flags.overrep_sequences {
                for seq in &acc.kmer_seqs {
                    if overrep_sampled < OVERREP_SAMPLE {
                        let mut key = [0u8; 50];
                        let n = seq.len().min(50);
                        key[..n].copy_from_slice(&seq[..n]);
                        *overrep_map.entry(key).or_insert(0) += 1;
                        overrep_sampled += 1;
                    }
                }
            }
            if config.flags.kmer_analysis {
                let kmers = count_kmers_parallel(&acc.kmer_seqs);
                for (k, v) in kmers {
                    *kmer_total.entry(k).or_insert(0) += v;
                }
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
        merge_batch_into_totals(
            &acc,
            &mut read_count, &mut total_bases, &mut gc_count,
            &mut quality_sum, &mut quality_bases, &mut q20_bases, &mut q30_bases,
            &mut adapter_hits, &mut min_length, &mut max_length,
            &mut length_histogram, &mut quality_by_pos,
            &mut base_composition, &mut quality_distribution, &mut per_tile,
        );

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

            if total_flushed.is_multiple_of(500_000) {
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

    // Deterministic dup rate: count reads whose fingerprint appeared >1 time in sample
    let dup_count: u64 = dup_map.values().map(|&c| if c > 1 { (c - 1) as u64 } else { 0 }).sum();
    let dup_rate = if dup_sampled > 0 { dup_count as f64 / dup_sampled as f64 * 100.0 } else { 0.0 };

    let overrep_seqs = finish_overrep(overrep_map, overrep_sampled);

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
        cur.overrepresented_sequences = overrep_seqs;
        cur.module_status = compute_module_status(cur);
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
        "Done: {} reads | GC {:.1}% | Q{:.1} | {:.2}% adapter | dup {:.1}% | {:.1} MB/s{}{}",
        format_number(read_count), gc_pct, avg_q, adapter_pct, dup_rate, mbps, trim_note, tile_note,
    ));

    if let Some(done) = s.current.take() {
        s.completed_files.push(done);
    }
}

// ---------------------------------------------------------------------------
// mmap-based zero-copy processing for plain FASTQ files
// Pipeline: dedicated reader thread fills a bounded channel while the main
// thread drains it with rayon — I/O and compute overlap fully.
// ---------------------------------------------------------------------------

fn process_single_file_mmap(
    file_path: String,
    state: Arc<Mutex<SharedState>>,
    config: &ProcessConfig,
    trim_path: Option<String>,
    file_size: u64,
) {
    let mut mmap_reader = match mmap_reader::MmapFastq::open(&file_path) {
        Ok(r) => r,
        Err(e) => {
            let msg = format!("Cannot open {}: {}", file_path, e);
            let mut s = state.lock().unwrap();
            s.log(LogLevel::Error, msg.clone());
            s.status = ProcessingStatus::Error(msg);
            return;
        }
    };

    let mut trim_writer: Option<Box<dyn std::io::Write>> = if let Some(ref tp) = trim_path {
        match io::open_writer(tp) {
            Ok(w) => Some(w),
            Err(e) => {
                let msg = format!("Cannot open trim output {}: {}", tp, e);
                let mut s = state.lock().unwrap();
                s.log(LogLevel::Error, msg.clone());
                s.status = ProcessingStatus::Error(msg);
                return;
            }
        }
    } else { None };

    // Running totals
    let mut read_count        = 0u64;
    let mut total_bases       = 0u64;
    let mut gc_count          = 0u64;
    let mut quality_sum       = 0u64;
    let mut quality_bases     = 0u64;
    let mut q20_bases         = 0u64;
    let mut q30_bases         = 0u64;
    let mut adapter_hits      = 0u64;
    let mut trimmed_reads     = 0u64;
    let mut trimmed_bases_removed = 0u64;
    let mut min_length        = u64::MAX;
    let mut max_length        = 0u64;
    let mut bytes             = 0u64;
    let mut length_histogram: HashMap<u64, u64> = HashMap::with_capacity(512);
    let mut quality_by_pos    = vec![(0u64, 0u64); crate::types::MAX_QUAL_POSITION];
    let mut base_composition  = vec![[0u64; 5]; crate::types::MAX_QUAL_POSITION];
    let mut quality_distribution = vec![0u64; PHRED_BUCKETS];
    let mut per_tile: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut kmer_total: HashMap<[u8; 4], u64> = HashMap::with_capacity(256);
    let mut dup_map: HashMap<u64, u32> = HashMap::with_capacity(OVERREP_SAMPLE);
    let mut dup_sampled = 0usize;
    let mut overrep_map: HashMap<[u8; 50], u64> = HashMap::with_capacity(4096);
    let mut overrep_sampled = 0usize;
    let mut flush_counter     = 0u64;
    let mut total_flushed     = 0u64;

    // Arc<Mmap> shared between the producer thread (via MmapFastq) and the
    // consumer (main thread).  Both deref to &[u8] for zero-copy access.
    let mmap_arc = mmap_reader.mmap_arc();

    // Bounded channel: 8 batches in flight.  Wider buffer means the producer
    // thread rarely stalls waiting for the consumer — better I/O / compute
    // overlap on NVMe or fast HDDs where seeks are cheap.
    let (tx, rx) = crossbeam_channel::bounded::<(Vec<RecordRange>, u64)>(8);

    // --- Producer thread: reads RecordRange batches, no data copies ---
    let producer = thread::spawn(move || {
        loop {
            let mut batch: Vec<RecordRange> = Vec::with_capacity(PARALLEL_BATCH);
            let mut batch_bytes = 0u64;
            while batch.len() < PARALLEL_BATCH {
                match mmap_reader.next_range() {
                    None => break,
                    Some(range) => {
                        // '@' + id + '\n' + seq + '\n' + "+\n" + qual + '\n' = +6 overhead
                        batch_bytes += range.id_len as u64 + range.seq_len as u64
                            + range.qual_len as u64 + 6;
                        batch.push(range);
                    }
                }
            }
            if batch.is_empty() { break; }
            if tx.send((batch, batch_bytes)).is_err() { break; }
        }
        // tx dropped here signals end-of-stream to the consumer
    });

    // --- Consumer: receive batches, process with rayon while producer reads ahead ---
    // Deref Arc<Mmap> → &[u8] once; safe to share across all rayon workers
    // because Mmap is Sync and outlives all rayon fold calls below.
    let mmap_bytes: &[u8] = &mmap_arc;

    while let Ok((batch, batch_bytes)) = rx.recv() {
        bytes += batch_bytes;

        let per_tile_enabled = config.flags.per_tile_quality;
        let acc: BatchAccum = batch
            .par_iter()
            .fold(BatchAccum::default, move |mut acc, range| {
                let id   = range.id(mmap_bytes);
                let seq  = range.seq(mmap_bytes);
                let qual = range.qual(mmap_bytes);
                let (eff_seq, eff_qual, adapter_hit) =
                    trim_record(seq, Some(qual), config);
                let tile_info = if per_tile_enabled {
                    parse_illumina_tile(id).and_then(|tile_id| {
                        eff_qual.map(|q| {
                            let phred_sum: u64 = q.iter().map(|&b| b.saturating_sub(33) as u64).sum();
                            (tile_id, phred_sum, q.len() as u64)
                        })
                    })
                } else {
                    None
                };
                let fp = fingerprint(eff_seq);
                acc.process(eff_seq, eff_qual, tile_info, fp, adapter_hit);
                acc
            })
            .reduce(BatchAccum::default, BatchAccum::merge);

        if config.flags.duplication_check {
            for &fp in &acc.fingerprints {
                if dup_sampled < OVERREP_SAMPLE {
                    *dup_map.entry(fp).or_insert(0) += 1;
                    dup_sampled += 1;
                }
            }
        }
        if !acc.kmer_seqs.is_empty() {
            if config.flags.overrep_sequences {
                for seq in &acc.kmer_seqs {
                    if overrep_sampled < OVERREP_SAMPLE {
                        let mut key = [0u8; 50];
                        let n = seq.len().min(50);
                        key[..n].copy_from_slice(&seq[..n]);
                        *overrep_map.entry(key).or_insert(0) += 1;
                        overrep_sampled += 1;
                    }
                }
            }
            if config.flags.kmer_analysis {
                let kmers = count_kmers_parallel(&acc.kmer_seqs);
                for (k, v) in kmers { *kmer_total.entry(k).or_insert(0) += v; }
            }
        }

        if config.trim_output {
            for range in &batch {
                let seq  = range.seq(mmap_bytes);
                let id   = range.id(mmap_bytes);
                let qual = range.qual(mmap_bytes);
                let (trimmed_seq, trimmed_qual, adapter_hit) =
                    trim_record(seq, Some(qual), config);
                if adapter_hit {
                    trimmed_bases_removed += seq.len() as u64 - trimmed_seq.len() as u64;
                }
                if trimmed_seq.len() as u64 >= config.min_length {
                    if adapter_hit { trimmed_reads += 1; }
                    if let Some(ref mut w) = trim_writer {
                        if let Err(e) = io::write_fastq_record(w, id, trimmed_seq, trimmed_qual) {
                            state.lock().unwrap().log(LogLevel::Warning, format!("Trim write error: {}", e));
                        }
                    }
                }
            }
        }

        merge_batch_into_totals(
            &acc,
            &mut read_count, &mut total_bases, &mut gc_count,
            &mut quality_sum, &mut quality_bases, &mut q20_bases, &mut q30_bases,
            &mut adapter_hits, &mut min_length, &mut max_length,
            &mut length_histogram, &mut quality_by_pos,
            &mut base_composition, &mut quality_distribution, &mut per_tile,
        );

        flush_counter += acc.read_count;
        if flush_counter >= FLUSH_INTERVAL {
            flush_counter = 0;
            total_flushed += FLUSH_INTERVAL;
            let mut s = state.lock().unwrap();
            if let Some(ref mut cur) = s.current {
                cur.read_count        = read_count;
                cur.total_bases       = total_bases;
                cur.gc_count          = gc_count;
                cur.quality_sum       = quality_sum;
                cur.quality_bases     = quality_bases;
                cur.q20_bases         = q20_bases;
                cur.q30_bases         = q30_bases;
                cur.adapter_hits      = adapter_hits;
                cur.min_length        = min_length;
                cur.max_length        = max_length;
                cur.bytes_processed   = bytes.min(file_size);
            }
            if total_flushed.is_multiple_of(500_000) {
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

    producer.join().expect("reader thread panicked");
    drop(trim_writer);

    let dup_count: u64 = dup_map.values().map(|&c| if c > 1 { (c - 1) as u64 } else { 0 }).sum();
    let dup_rate = if dup_sampled > 0 { dup_count as f64 / dup_sampled as f64 * 100.0 } else { 0.0 };

    let overrep_seqs = finish_overrep(overrep_map, overrep_sampled);

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
        cur.overrepresented_sequences = overrep_seqs;
        cur.module_status = compute_module_status(cur);
    }

    let elapsed = s.elapsed_secs();
    let mbps = if elapsed > 0.0 { bytes as f64 / 1_048_576.0 / elapsed } else { 0.0 };
    let gc_pct = if total_bases > 0 { gc_count as f64 / total_bases as f64 * 100.0 } else { 0.0 };
    let avg_q  = if quality_bases > 0 { quality_sum as f64 / quality_bases as f64 } else { 0.0 };
    let adapter_pct = if read_count > 0 { adapter_hits as f64 / read_count as f64 * 100.0 } else { 0.0 };

    s.log(LogLevel::Success, format!(
        "Done: {} reads | GC {:.1}% | Q{:.1} | {:.2}% adapter | dup {:.1}% | {:.1} MB/s",
        format_number(read_count), gc_pct, avg_q, adapter_pct, dup_rate, mbps,
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
            format!("{}/{}_trimmed.fastq.gz", config.output_dir, output_stem(&r1_path))
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
    let mut dup_map_pe: HashMap<u64, u32> = HashMap::with_capacity(OVERREP_SAMPLE);
    let mut dup_sampled_pe = 0usize;
    let mut overrep_map_pe: HashMap<[u8; 50], u64> = HashMap::with_capacity(4096);
    let mut overrep_sampled_pe = 0usize;
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
        let per_tile_enabled = config.flags.per_tile_quality;
        let acc: BatchAccum = batch
            .par_iter()
            .fold(BatchAccum::default, move |mut acc, (r1, r2)| {
                for rec in [r1, r2] {
                    let (eff_seq, eff_qual, adapter_hit) =
                        trim_record(&rec.seq, rec.qual.as_deref(), config);
                    let tile_info = if per_tile_enabled {
                        parse_illumina_tile(&rec.id).and_then(|tile_id| {
                            eff_qual.map(|q| {
                                let phred_sum: u64 = q.iter().map(|&b| b.saturating_sub(33) as u64).sum();
                                (tile_id, phred_sum, q.len() as u64)
                            })
                        })
                    } else {
                        None
                    };
                    let fp = fingerprint(eff_seq);
                    acc.process(eff_seq, eff_qual, tile_info, fp, adapter_hit);
                }
                acc
            })
            .reduce(BatchAccum::default, BatchAccum::merge);

        // Dedup + overrep
        if config.flags.duplication_check {
            for &fp in &acc.fingerprints {
                if dup_sampled_pe < OVERREP_SAMPLE {
                    *dup_map_pe.entry(fp).or_insert(0) += 1;
                    dup_sampled_pe += 1;
                }
            }
        }
        if !acc.kmer_seqs.is_empty() {
            if config.flags.overrep_sequences {
                for seq in &acc.kmer_seqs {
                    if overrep_sampled_pe < OVERREP_SAMPLE {
                        let mut key = [0u8; 50];
                        let n = seq.len().min(50);
                        key[..n].copy_from_slice(&seq[..n]);
                        *overrep_map_pe.entry(key).or_insert(0) += 1;
                        overrep_sampled_pe += 1;
                    }
                }
            }
            if config.flags.kmer_analysis {
                let kmers = count_kmers_parallel(&acc.kmer_seqs);
                for (k, v) in kmers { *kmer_total.entry(k).or_insert(0) += v; }
            }
        }

        merge_batch_into_totals(
            &acc,
            &mut read_count, &mut total_bases, &mut gc_count,
            &mut quality_sum, &mut quality_bases, &mut q20_bases, &mut q30_bases,
            &mut adapter_hits, &mut min_length, &mut max_length,
            &mut length_histogram, &mut quality_by_pos,
            &mut base_composition, &mut quality_distribution, &mut per_tile,
        );

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
            if total_flushed.is_multiple_of(500_000) {
                let elapsed = s.elapsed_secs();
                let rps = if elapsed > 0.0 { read_count as f64 / elapsed } else { 0.0 };
                s.log(LogLevel::Info, format!("{} reads  {:.0} pairs/s", format_number(read_count), rps / 2.0));
            }
        }
    }

    let dup_pe_count: u64 = dup_map_pe.values().map(|&c| if c > 1 { (c-1) as u64 } else { 0 }).sum();
    let dup_rate = if dup_sampled_pe > 0 { dup_pe_count as f64 / dup_sampled_pe as f64 * 100.0 } else { 0.0 };
    let overrep_seqs_pe = finish_overrep(overrep_map_pe, overrep_sampled_pe);

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
        cur.overrepresented_sequences = overrep_seqs_pe;
        cur.module_status = compute_module_status(cur);
    }

    let elapsed = s.elapsed_secs();
    let pairs = read_count / 2;
    let mbps = if elapsed > 0.0 { bytes as f64 / 1_048_576.0 / elapsed } else { 0.0 };
    let gc_pct = if total_bases > 0 { gc_count as f64 / total_bases as f64 * 100.0 } else { 0.0 };
    let avg_q = if quality_bases > 0 { quality_sum as f64 / quality_bases as f64 } else { 0.0 };

    s.log(LogLevel::Success, format!(
        "Done: {} pairs | GC {:.1}% | Q{:.1} | dup {:.1}% | {:.1} MB/s",
        format_number(pairs), gc_pct, avg_q, dup_rate, mbps,
    ));

    if let Some(done) = s.current.take() {
        s.completed_files.push(done);
    }
}
