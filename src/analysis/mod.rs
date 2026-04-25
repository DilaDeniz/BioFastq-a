mod adapters;
mod io;
mod kmers;
mod metrics;
mod mmap_reader;


use std::collections::HashMap;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread;

use hyperloglog::HyperLogLog;
use needletail::parse_fastx_file;
use needletail::Sequence;
use rayon::prelude::*;

use crate::types::{
    format_number, FileStats, LogLevel, OverrepSeq, ProcessConfig, ProcessingStatus,
    QcStatus, SharedState, ADAPTER_MATCH_LEN, ADAPTERS, FLUSH_INTERVAL,
    MAX_QUAL_VS_LEN_POINTS, OVERREP_SAMPLE, PARALLEL_BATCH, PHRED_BUCKETS,
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

/// Apply all trimming steps in biologically correct order.
/// Returns None if the read is filtered out (fails quality/N/length checks).
/// Returns Some((trimmed_seq, trimmed_qual, adapter_name)) if the read passes.
#[allow(clippy::type_complexity)]
#[inline]
fn trim_record<'a>(
    seq: &'a [u8],
    qual: Option<&'a [u8]>,
    config: &ProcessConfig,
) -> Option<(&'a [u8], Option<&'a [u8]>, Option<&'static str>)> {
    let mut seq = seq;
    let mut qual = qual;

    // 1. Hard global trim (front, then tail)
    let front = config.trim_front_bases as usize;
    let tail  = config.trim_tail_bases as usize;
    if front > 0 && seq.len() > front {
        seq  = &seq[front..];
        qual = qual.map(|q| &q[front.min(q.len())..]);
    } else if front >= seq.len() && front > 0 {
        return None; // trimmed to nothing
    }
    if tail > 0 && seq.len() > tail {
        let cut = seq.len() - tail;
        seq  = &seq[..cut];
        qual = qual.map(|q| &q[..cut.min(q.len())]);
    } else if tail >= seq.len() && tail > 0 {
        return None;
    }

    // 2. Sliding window quality trimming
    if let Some(q) = qual {
        // cut_front: advance start past low-quality 5' bases
        if config.cut_front_window > 0 {
            let start = adapters::sliding_window_cut_front(
                q, config.cut_front_window as usize, config.cut_front_qual);
            if start >= seq.len() { return None; }
            seq  = &seq[start..];
            qual = Some(&q[start..]);
        }
        // cut_right: cut from first bad window to read end
        if config.cut_right_window > 0 {
            let q2 = qual.unwrap_or(q);
            let cut = adapters::sliding_window_cut_right(
                q2, config.cut_right_window as usize, config.cut_right_qual);
            if cut == 0 { return None; }
            seq  = &seq[..cut.min(seq.len())];
            qual = Some(&q2[..cut]);
        }
        // cut_tail: shrink end past low-quality 3' bases
        if config.cut_tail_window > 0 {
            let q3 = qual.unwrap_or(q);
            let cut = adapters::sliding_window_cut_tail(
                q3, config.cut_tail_window as usize, config.cut_tail_qual);
            if cut == 0 { return None; }
            seq  = &seq[..cut.min(seq.len())];
            qual = Some(&q3[..cut]);
        }
    }

    // 3. PolyG trimming (3' end — 2-color Illumina artifact)
    if config.poly_g_min_len > 0 {
        let cut = adapters::trim_poly_base(seq, b'G', config.poly_g_min_len as usize);
        seq  = &seq[..cut];
        qual = qual.map(|q| &q[..cut.min(q.len())]);
    }

    // 4. Existing per-base 3' quality trim
    if config.quality_trim_threshold > 0 {
        if let Some(q) = qual {
            let cut = adapters::quality_trim_3p(q, config.quality_trim_threshold);
            seq  = &seq[..cut.min(seq.len())];
            qual = Some(&q[..cut]);
        }
    }

    // 5. Adapter detection & trim
    let adapter_name = if config.flags.adapter_detection {
        match adapters::find_adapter_with_name(seq, &config.custom_adapters) {
            Some((pos, name)) => {
                seq  = &seq[..pos];
                qual = qual.map(|q| &q[..pos.min(q.len())]);
                Some(name)
            }
            None => None,
        }
    } else {
        None
    };

    // 6. PolyX trimming (any homopolymer at 3' end)
    if config.poly_x_min_len > 0 && !seq.is_empty() {
        // Detect dominant 3' base and trim if it forms a run
        let dominant_base = *seq.last().unwrap();
        // Only trim if it's not already handled by polyG
        if dominant_base != b'G' || config.poly_g_min_len == 0 {
            let cut = adapters::trim_poly_base(seq, dominant_base, config.poly_x_min_len as usize);
            seq  = &seq[..cut];
            qual = qual.map(|q| &q[..cut.min(q.len())]);
        }
    }

    // 7. Per-read filters (after all trimming)
    if !adapters::passes_filters(seq, qual, config.min_avg_quality, config.max_n_bases) {
        return None;
    }
    if seq.len() < config.min_length as usize {
        return None;
    }

    Some((seq, qual, adapter_name))
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
/// Uses an adaptive threshold: max(0.1% of sampled, 99th-percentile of counts).
/// This prevents noisy hits on small datasets and adapts to sequencing depth.
fn finish_overrep(
    map: HashMap<[u8; 50], u64>,
    sampled: usize,
) -> Vec<OverrepSeq> {
    if sampled == 0 { return Vec::new(); }
    // Collect all counts and compute 99th percentile as adaptive floor
    let mut all_counts: Vec<u64> = map.values().cloned().collect();
    all_counts.sort_unstable();
    let p99_idx = (all_counts.len() as f64 * 0.99) as usize;
    let p99_threshold = all_counts.get(p99_idx).cloned().unwrap_or(1);
    // Final threshold: at least 0.1% of reads AND at least the 99th percentile count
    let pct_threshold = (sampled as f64 * 0.001).ceil() as u64;
    let threshold = pct_threshold.max(p99_threshold).max(5);
    let mut sorted: Vec<([u8; 50], u64)> = map
        .into_iter()
        .filter(|(_, c)| *c >= threshold)
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
    let lr = f.long_read_mode;

    // Per-base sequence quality — worst position average
    // Long-read: ONT Q-scores are lower by design (Q8-Q20 typical)
    let qual_per_pos = f.avg_qual_per_position();
    let min_pos_q = qual_per_pos.iter().cloned().fold(f64::MAX, f64::min);
    if lr {
        m.push(("Per-base quality".into(), if min_pos_q >= 10.0 { QcStatus::Pass }
            else if min_pos_q >= 8.0 { QcStatus::Warn } else { QcStatus::Fail }));
    } else {
        m.push(("Per-base quality".into(), if min_pos_q >= 28.0 { QcStatus::Pass }
            else if min_pos_q >= 20.0 { QcStatus::Warn } else { QcStatus::Fail }));
    }

    // Per-sequence quality
    // Long-read mode: use Q10 warn / Q8 fail thresholds; compare mean quality not Q30%
    if lr {
        let avg_q = f.avg_quality();
        m.push(("Per-sequence quality".into(), if avg_q >= 10.0 { QcStatus::Pass }
            else if avg_q >= 8.0 { QcStatus::Warn } else { QcStatus::Fail }));
    } else {
        let q30 = f.q30_pct();
        m.push(("Per-sequence quality".into(), if q30 >= 80.0 { QcStatus::Pass }
            else if q30 >= 60.0 { QcStatus::Warn } else { QcStatus::Fail }));
    }

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

    // Sequence length distribution.
    // Long-read: variable length is expected and normal → always Pass.
    // Short-read: len_range == 0 (all same) → Pass; any variation → Warn.
    let len_range = f.max_length.saturating_sub(f.effective_min_length());
    m.push(("Sequence length".into(),
        if lr || len_range == 0 { QcStatus::Pass } else { QcStatus::Warn }));

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
    let mut reads_filtered = 0u64;
    let mut min_length = u64::MAX;
    let mut max_length = 0u64;
    let mut bytes = 0u64;
    let mut length_histogram: HashMap<u64, u64> = HashMap::with_capacity(512);
    let mut quality_by_pos = vec![(0u64, 0u64); crate::types::MAX_QUAL_POSITION];
    let mut qual_hist_by_pos = vec![[0u64; 43]; crate::types::MAX_QUAL_POSITION];
    let mut per_tile: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut kmer_total: HashMap<[u8; 4], u64> = HashMap::with_capacity(256);
    let mut base_composition = vec![[0u64; 5]; crate::types::MAX_QUAL_POSITION];
    let mut quality_distribution = vec![0u64; PHRED_BUCKETS];
    let mut quality_by_length_bin: HashMap<u32, (u64, u64)> = HashMap::new();
    // HyperLogLog for global duplicate estimation across ALL reads (1% error, ~2KB memory)
    let mut hll: HyperLogLog = HyperLogLog::new(0.01);
    let mut overrep_map: HashMap<[u8; 50], u64> = HashMap::with_capacity(4096);
    let mut overrep_sampled = 0usize;
    let mut trimmed_by_adapter: HashMap<String, u64> = HashMap::new();
    let mut flush_counter = 0u64;
    let mut total_flushed = 0u64;
    let mut error_occurred = false;
    let mut gc_distribution = [0u64; 101];
    let mut dup_sample: HashMap<u64, u32> = HashMap::with_capacity(8192);
    let mut dup_sample_count = 0usize;
    // Long-read / ONT accumulators
    let mut ont_channel_counts: HashMap<u32, u32> = HashMap::new();
    let mut reads_over_time_raw: Vec<(u64, u64)> = Vec::new();
    let mut qual_vs_length: Vec<(u32, f32)> = Vec::new();
    // Auto-detection state
    let mut first1k_lengths: Vec<u64> = Vec::with_capacity(1001);
    let mut ont_header_detected = false;
    let mut long_read_detected = false;
    let mut detection_done = false;

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

        // --- Auto-detect long-read mode from first 1000 reads ---
        if !detection_done && read_count + batch.len() as u64 >= 1000 {
            detection_done = true;
        }
        if !detection_done || (read_count < 1000 && !first1k_lengths.is_empty()) {
            for rec in &batch {
                if first1k_lengths.len() < 1000 {
                    first1k_lengths.push(rec.seq.len() as u64);
                    if !ont_header_detected && is_ont_header(&rec.id) {
                        ont_header_detected = true;
                    }
                }
            }
            if first1k_lengths.len() >= 1000 && !long_read_detected {
                let med = median_length(&first1k_lengths);
                if med > 1000 {
                    long_read_detected = true;
                }
            }
        }

        let is_lr = config.is_long_read || long_read_detected || ont_header_detected;

        // --- Parallel stats computation ---
        let cfg = config;
        let per_tile_enabled = cfg.flags.per_tile_quality;
        let cur_read_count = read_count;
        let cur_qvl_len = qual_vs_length.len();
        let acc: BatchAccum = batch
            .par_iter()
            .enumerate()
            .fold(BatchAccum::default, move |mut acc, (bi, rec)| {
                let Some((effective_seq, effective_qual, adapter_name)) =
                    trim_record(&rec.seq, rec.qual.as_deref(), cfg)
                else {
                    acc.reads_filtered += 1;
                    return acc;
                };
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
                acc.process(effective_seq, effective_qual, tile_info, fp, adapter_name.is_some());
                // Long-read ONT data
                if is_lr {
                    if let Some(ch) = parse_ont_channel(&rec.id) {
                        *acc.ont_channel_counts.entry(ch).or_insert(0) += 1;
                    }
                    if let Some(ts) = parse_ont_start_time(&rec.id) {
                        acc.reads_over_time_raw.push((ts, 1));
                    }
                    // Sample qual-vs-length: sample every Nth read to keep ~2000 points
                    // N is chosen so that with the current batch offset we stay under cap
                    let global_idx = cur_read_count + bi as u64;
                    // sample 1 in every max(1, total_cap_factor) reads
                    // Use a simple modulo: sample if global_idx % sample_every == 0
                    // We don't know total reads; keep at most MAX_QUAL_VS_LEN_POINTS
                    let already = cur_qvl_len + acc.qual_vs_length.len();
                    let sample_n: u64 = if already < MAX_QUAL_VS_LEN_POINTS {
                        // Rough estimate: sample every ~(total_est / max) reads
                        // Use simple approach: every 50th read (gives 2000 for ~100k reads)
                        50
                    } else {
                        u64::MAX // stop sampling
                    };
                    if sample_n < u64::MAX && global_idx.is_multiple_of(sample_n) {
                        if let Some(q) = effective_qual {
                            if !q.is_empty() {
                                let phred_sum: u64 = q.iter().map(|&b| b.saturating_sub(33) as u64).sum();
                                let mean_q = phred_sum as f32 / q.len() as f32;
                                acc.qual_vs_length.push((effective_seq.len() as u32, mean_q));
                            }
                        }
                    }
                }
                acc
            })
            .reduce(BatchAccum::default, BatchAccum::merge);

        // HyperLogLog: insert ALL fingerprints (no sampling cap)
        if config.flags.duplication_check {
            for &fp in &acc.fingerprints {
                hll.insert(&fp);
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
                if let Some((trimmed_seq, trimmed_qual, adapter_name)) =
                    trim_record(&rec.seq, rec.qual.as_deref(), config)
                {
                    if let Some(name) = adapter_name {
                        trimmed_bases_removed += rec.seq.len() as u64 - trimmed_seq.len() as u64;
                        *trimmed_by_adapter.entry(name.to_string()).or_insert(0) += 1;
                        trimmed_reads += 1;
                    }
                    if let Some(ref mut w) = trim_writer {
                        if let Err(e) = write_fastq_record(w, &rec.id, trimmed_seq, trimmed_qual) {
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
            &mut length_histogram, &mut quality_by_pos, &mut qual_hist_by_pos,
            &mut base_composition, &mut quality_distribution, &mut per_tile,
            &mut quality_by_length_bin,
        );
        reads_filtered += acc.reads_filtered;

        // Merge long-read data from batch accum
        for (ch, cnt) in acc.ont_channel_counts {
            *ont_channel_counts.entry(ch).or_insert(0) += cnt;
        }
        // Truncate qual_vs_length to MAX_QUAL_VS_LEN_POINTS
        if qual_vs_length.len() < MAX_QUAL_VS_LEN_POINTS {
            let remaining = MAX_QUAL_VS_LEN_POINTS - qual_vs_length.len();
            let take = acc.qual_vs_length.len().min(remaining);
            qual_vs_length.extend_from_slice(&acc.qual_vs_length[..take]);
        }
        reads_over_time_raw.extend(acc.reads_over_time_raw);

        // GC per-read distribution (reuses already-computed per-batch data)
        for (dst, &src) in gc_distribution.iter_mut().zip(acc.gc_per_read.iter()) {
            *dst += src;
        }
        // Dup level sample (capped at 200k reads for histogram)
        if config.flags.duplication_check && dup_sample_count < OVERREP_SAMPLE {
            for &fp in &acc.fingerprints {
                if dup_sample_count >= OVERREP_SAMPLE { break; }
                *dup_sample.entry(fp).or_insert(0) += 1;
                dup_sample_count += 1;
            }
        }

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

    // Finalize long-read mode
    let is_long_read_final = config.is_long_read || long_read_detected || ont_header_detected;
    let final_reads_over_time = bucket_reads_over_time(reads_over_time_raw);

    // HyperLogLog dup rate: (total - estimated_unique) / total
    let dup_rate = if config.flags.duplication_check && read_count > 0 {
        let unique_est = hll.len() as u64;
        read_count.saturating_sub(unique_est) as f64 / read_count as f64 * 100.0
    } else { 0.0 };

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
        cur.reads_filtered        = reads_filtered;
        cur.min_length            = min_length;
        cur.max_length            = max_length;
        cur.bytes_processed       = file_size;
        cur.length_histogram      = length_histogram;
        cur.quality_by_position   = quality_by_pos;
        cur.qual_hist_by_position = qual_hist_by_pos;
        cur.kmer_counts           = kmer_total;
        cur.dup_rate_pct          = dup_rate;
        cur.per_tile_quality      = per_tile;
        cur.base_composition      = base_composition;
        cur.quality_distribution  = quality_distribution;
        cur.quality_by_length_bin = quality_by_length_bin;
        cur.trimmed_by_adapter    = trimmed_by_adapter;
        cur.overrepresented_sequences = overrep_seqs;
        cur.gc_distribution = gc_distribution.to_vec();
        cur.dup_level_histogram = compute_dup_histogram(&dup_sample);
        cur.long_read_mode = is_long_read_final;
        cur.ont_channel_counts = ont_channel_counts;
        cur.reads_over_time = final_reads_over_time;
        cur.qual_vs_length = qual_vs_length;
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
    let mut reads_filtered    = 0u64;
    let mut min_length        = u64::MAX;
    let mut max_length        = 0u64;
    let mut bytes             = 0u64;
    let mut length_histogram: HashMap<u64, u64> = HashMap::with_capacity(512);
    let mut quality_by_pos    = vec![(0u64, 0u64); crate::types::MAX_QUAL_POSITION];
    let mut qual_hist_by_pos  = vec![[0u64; 43]; crate::types::MAX_QUAL_POSITION];
    let mut base_composition  = vec![[0u64; 5]; crate::types::MAX_QUAL_POSITION];
    let mut quality_distribution = vec![0u64; PHRED_BUCKETS];
    let mut quality_by_length_bin: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut per_tile: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut kmer_total: HashMap<[u8; 4], u64> = HashMap::with_capacity(256);
    let mut hll: HyperLogLog = HyperLogLog::new(0.01);
    let mut overrep_map: HashMap<[u8; 50], u64> = HashMap::with_capacity(4096);
    let mut overrep_sampled = 0usize;
    let mut trimmed_by_adapter: HashMap<String, u64> = HashMap::new();
    let mut flush_counter     = 0u64;
    let mut total_flushed     = 0u64;
    let mut gc_distribution = [0u64; 101];
    let mut dup_sample: HashMap<u64, u32> = HashMap::with_capacity(8192);
    let mut dup_sample_count = 0usize;
    // Long-read / ONT accumulators (mmap path)
    let mut ont_channel_counts_mm: HashMap<u32, u32> = HashMap::new();
    let mut reads_over_time_raw_mm: Vec<(u64, u64)> = Vec::new();
    let mut qual_vs_length_mm: Vec<(u32, f32)> = Vec::new();
    // Auto-detection state (mmap path)
    let mut first1k_lengths_mm: Vec<u64> = Vec::with_capacity(1001);
    let mut ont_header_detected_mm = false;
    let mut long_read_detected_mm = false;
    let mut detection_done_mm = false;

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

        // Auto-detect long-read mode from first 1000 reads (mmap path)
        if !detection_done_mm && read_count + batch.len() as u64 >= 1000 {
            detection_done_mm = true;
        }
        if !detection_done_mm || (read_count < 1000 && !first1k_lengths_mm.is_empty()) {
            for range in &batch {
                if first1k_lengths_mm.len() < 1000 {
                    first1k_lengths_mm.push(range.seq_len as u64);
                    if !ont_header_detected_mm {
                        let id = range.id(mmap_bytes);
                        if is_ont_header(id) {
                            ont_header_detected_mm = true;
                        }
                    }
                }
            }
            if first1k_lengths_mm.len() >= 1000 && !long_read_detected_mm {
                let med = median_length(&first1k_lengths_mm);
                if med > 1000 {
                    long_read_detected_mm = true;
                }
            }
        }

        let is_lr_mm = config.is_long_read || long_read_detected_mm || ont_header_detected_mm;
        let per_tile_enabled = config.flags.per_tile_quality;
        let cur_read_count_mm = read_count;
        let cur_qvl_len_mm = qual_vs_length_mm.len();
        let acc: BatchAccum = batch
            .par_iter()
            .enumerate()
            .fold(BatchAccum::default, move |mut acc, (bi, range)| {
                let id   = range.id(mmap_bytes);
                let seq  = range.seq(mmap_bytes);
                let qual = range.qual(mmap_bytes);
                let Some((eff_seq, eff_qual, adapter_name)) =
                    trim_record(seq, Some(qual), config)
                else {
                    acc.reads_filtered += 1;
                    return acc;
                };
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
                acc.process(eff_seq, eff_qual, tile_info, fp, adapter_name.is_some());
                // Long-read ONT data (mmap path)
                if is_lr_mm {
                    if let Some(ch) = parse_ont_channel(id) {
                        *acc.ont_channel_counts.entry(ch).or_insert(0) += 1;
                    }
                    if let Some(ts) = parse_ont_start_time(id) {
                        acc.reads_over_time_raw.push((ts, 1));
                    }
                    let global_idx = cur_read_count_mm + bi as u64;
                    let already = cur_qvl_len_mm + acc.qual_vs_length.len();
                    let sample_n: u64 = if already < MAX_QUAL_VS_LEN_POINTS { 50 } else { u64::MAX };
                    if sample_n < u64::MAX && global_idx.is_multiple_of(sample_n) {
                        if let Some(q) = eff_qual {
                            if !q.is_empty() {
                                let phred_sum: u64 = q.iter().map(|&b| b.saturating_sub(33) as u64).sum();
                                let mean_q = phred_sum as f32 / q.len() as f32;
                                acc.qual_vs_length.push((eff_seq.len() as u32, mean_q));
                            }
                        }
                    }
                }
                acc
            })
            .reduce(BatchAccum::default, BatchAccum::merge);

        if config.flags.duplication_check {
            for &fp in &acc.fingerprints {
                hll.insert(&fp);
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
                if let Some((trimmed_seq, trimmed_qual, adapter_name)) =
                    trim_record(seq, Some(qual), config)
                {
                    if let Some(name) = adapter_name {
                        trimmed_bases_removed += seq.len() as u64 - trimmed_seq.len() as u64;
                        *trimmed_by_adapter.entry(name.to_string()).or_insert(0) += 1;
                        trimmed_reads += 1;
                    }
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
            &mut length_histogram, &mut quality_by_pos, &mut qual_hist_by_pos,
            &mut base_composition, &mut quality_distribution, &mut per_tile,
            &mut quality_by_length_bin,
        );
        reads_filtered += acc.reads_filtered;

        // Merge long-read data from batch accum (mmap path)
        for (ch, cnt) in acc.ont_channel_counts {
            *ont_channel_counts_mm.entry(ch).or_insert(0) += cnt;
        }
        if qual_vs_length_mm.len() < MAX_QUAL_VS_LEN_POINTS {
            let remaining = MAX_QUAL_VS_LEN_POINTS - qual_vs_length_mm.len();
            let take = acc.qual_vs_length.len().min(remaining);
            qual_vs_length_mm.extend_from_slice(&acc.qual_vs_length[..take]);
        }
        reads_over_time_raw_mm.extend(acc.reads_over_time_raw);

        // GC per-read distribution (reuses already-computed per-batch data)
        for (dst, &src) in gc_distribution.iter_mut().zip(acc.gc_per_read.iter()) {
            *dst += src;
        }
        // Dup level sample (capped at 200k reads for histogram)
        if config.flags.duplication_check && dup_sample_count < OVERREP_SAMPLE {
            for &fp in &acc.fingerprints {
                if dup_sample_count >= OVERREP_SAMPLE { break; }
                *dup_sample.entry(fp).or_insert(0) += 1;
                dup_sample_count += 1;
            }
        }

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

    // Finalize long-read mode (mmap path)
    let is_long_read_final_mm = config.is_long_read || long_read_detected_mm || ont_header_detected_mm;
    let final_reads_over_time_mm = bucket_reads_over_time(reads_over_time_raw_mm);

    let dup_rate = if config.flags.duplication_check && read_count > 0 {
        let unique_est = hll.len() as u64;
        read_count.saturating_sub(unique_est) as f64 / read_count as f64 * 100.0
    } else { 0.0 };

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
        cur.reads_filtered        = reads_filtered;
        cur.min_length            = min_length;
        cur.max_length            = max_length;
        cur.bytes_processed       = file_size;
        cur.length_histogram      = length_histogram;
        cur.quality_by_position   = quality_by_pos;
        cur.qual_hist_by_position = qual_hist_by_pos;
        cur.kmer_counts           = kmer_total;
        cur.dup_rate_pct          = dup_rate;
        cur.per_tile_quality      = per_tile;
        cur.base_composition      = base_composition;
        cur.quality_distribution  = quality_distribution;
        cur.quality_by_length_bin = quality_by_length_bin;
        cur.trimmed_by_adapter    = trimmed_by_adapter;
        cur.overrepresented_sequences = overrep_seqs;
        cur.gc_distribution = gc_distribution.to_vec();
        cur.dup_level_histogram = compute_dup_histogram(&dup_sample);
        cur.long_read_mode = is_long_read_final_mm;
        cur.ont_channel_counts = ont_channel_counts_mm;
        cur.reads_over_time = final_reads_over_time_mm;
        cur.qual_vs_length = qual_vs_length_mm;
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

fn compute_dup_histogram(sample: &HashMap<u64, u32>) -> Vec<f64> {
    if sample.is_empty() { return vec![0.0; 9]; }
    let total = sample.values().map(|&c| c as u64).sum::<u64>() as f64;
    let mut bins = [0u64; 9];
    for &count in sample.values() {
        let idx = match count {
            1 => 0, 2 => 1, 3..=4 => 2, 5..=9 => 3,
            10..=49 => 4, 50..=99 => 5, 100..=499 => 6, 500..=999 => 7, _ => 8,
        };
        bins[idx] += 1;
    }
    bins.iter().map(|&c| c as f64 / total * 100.0).collect()
}

fn parse_illumina_tile(id: &[u8]) -> Option<u32> {
    let s = std::str::from_utf8(id).ok()?;
    let s = s.split_ascii_whitespace().next()?;
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() >= 6 { parts[4].parse().ok() } else { None }
}

// ---------------------------------------------------------------------------
// ONT / long-read helpers
// ---------------------------------------------------------------------------

/// Parse `ch=<N>` from an ONT FASTQ read header.
/// ONT headers look like:
///   @<uuid> runid=<id> read=<N> ch=<channel> start_time=2024-01-15T10:30:00Z ...
/// Matches only when `ch=` is preceded by whitespace or is at the start (space-delimited field).
fn parse_ont_channel(id: &[u8]) -> Option<u32> {
    let s = std::str::from_utf8(id).ok()?;
    // Find " ch=" (space-prefixed) to avoid matching inside other tokens like "arch="
    let search = " ch=";
    let idx = s.find(search).map(|i| i + search.len())
        .or_else(|| s.strip_prefix("ch=").map(|_| 3))?;
    let after = &s[idx..];
    let end = after.find(|c: char| !c.is_ascii_digit()).unwrap_or(after.len());
    if end == 0 { return None; }
    after[..end].parse().ok()
}

/// Parse `start_time=<ISO8601>` from an ONT read header and return Unix seconds.
/// Handles the format: `2024-01-15T10:30:00Z` (UTC, Zulu suffix).
/// Returns None if the field is absent or unparseable.
fn parse_ont_start_time(id: &[u8]) -> Option<u64> {
    let s = std::str::from_utf8(id).ok()?;
    let idx = s.find("start_time=")?;
    let after = &s[idx + 11..];
    // Grab the timestamp token (no whitespace)
    let end = after.find(|c: char| c.is_ascii_whitespace()).unwrap_or(after.len());
    let ts = &after[..end];
    // Expected: YYYY-MM-DDTHH:MM:SS[Z] or YYYY-MM-DDTHH:MM:SS+00:00
    // Must be at least 19 chars: 2024-01-15T10:30:00
    if ts.len() < 19 { return None; }
    let year:  u64 = ts[0..4].parse().ok()?;
    let month: u64 = ts[5..7].parse().ok()?;
    let day:   u64 = ts[8..10].parse().ok()?;
    let hour:  u64 = ts[11..13].parse().ok()?;
    let min:   u64 = ts[14..16].parse().ok()?;
    let sec:   u64 = ts[17..19].parse().ok()?;
    // Validate ranges
    if month == 0 || month > 12 || day == 0 || day > 31 { return None; }
    if hour > 23 || min > 59 || sec > 60 { return None; }
    // Days from epoch (1970-01-01) to year-month-day using Gregorian calendar.
    // We use the algorithm: days since epoch = days_since_epoch(Y,M,D).
    let days = days_since_unix_epoch(year, month, day)?;
    Some(days * 86400 + hour * 3600 + min * 60 + sec)
}

/// Compute number of days from 1970-01-01 to Y-M-D (proleptic Gregorian calendar).
/// Returns None for implausible dates.
fn days_since_unix_epoch(year: u64, month: u64, day: u64) -> Option<u64> {
    if year < 1970 { return None; }
    // Days per month (non-leap)
    const DAYS_IN_MONTH: [u64; 13] = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    let is_leap = |y: u64| (y.is_multiple_of(4) && !y.is_multiple_of(100)) || y.is_multiple_of(400);
    let days_in_feb = |y: u64| if is_leap(y) { 29u64 } else { 28u64 };

    // Days from 1970 to beginning of `year`
    let mut total_days: u64 = 0;
    for y in 1970..year {
        total_days += if is_leap(y) { 366 } else { 365 };
    }
    // Add days for completed months in the current year
    for m in 1..month {
        total_days += if m == 2 { days_in_feb(year) } else { DAYS_IN_MONTH[m as usize] };
    }
    // Add days in the current month (day 1 = 0 extra days)
    total_days += day.checked_sub(1)?;
    Some(total_days)
}

/// Returns true when the header looks like an ONT read (contains ` ch=` AND `start_time=`).
fn is_ont_header(id: &[u8]) -> bool {
    let Ok(s) = std::str::from_utf8(id) else { return false; };
    (s.contains(" ch=") || s.starts_with("ch=")) && s.contains("start_time=")
}

/// Compute median read length from a sample of up to `n` lengths.
/// Returns 0 if the sample is empty.
fn median_length(lengths: &[u64]) -> u64 {
    if lengths.is_empty() { return 0; }
    let mut v = lengths.to_vec();
    v.sort_unstable();
    let mid = v.len() / 2;
    if v.len().is_multiple_of(2) {
        (v[mid - 1] + v[mid]) / 2
    } else {
        v[mid]
    }
}

/// Bucket raw (timestamp_secs, 1) events into 5-minute cumulative intervals.
/// Returns vec of (minutes_since_start, cumulative_reads).
fn bucket_reads_over_time(mut raw: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    if raw.is_empty() { return Vec::new(); }
    raw.sort_unstable_by_key(|&(ts, _)| ts);
    let t0 = raw[0].0;
    const BUCKET_SECS: u64 = 300; // 5 minutes
    let last_ts = raw.last().unwrap().0;
    let n_buckets = ((last_ts - t0) / BUCKET_SECS + 1) as usize;
    let mut buckets = vec![0u64; n_buckets];
    for (ts, _) in &raw {
        let b = ((ts - t0) / BUCKET_SECS) as usize;
        buckets[b.min(n_buckets - 1)] += 1;
    }
    // Convert to cumulative
    let mut result = Vec::with_capacity(n_buckets);
    let mut cumulative = 0u64;
    for (i, &count) in buckets.iter().enumerate() {
        cumulative += count;
        result.push(((i as u64) * 5, cumulative));
    }
    result
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
    let mut reads_filtered = 0u64;
    let mut min_length = u64::MAX;
    let mut max_length = 0u64;
    let mut bytes = 0u64;
    let mut length_histogram: HashMap<u64, u64> = HashMap::with_capacity(512);
    let mut quality_by_pos = vec![(0u64, 0u64); crate::types::MAX_QUAL_POSITION];
    let mut qual_hist_by_pos = vec![[0u64; 43]; crate::types::MAX_QUAL_POSITION];
    let mut base_composition = vec![[0u64; 5]; crate::types::MAX_QUAL_POSITION];
    let mut quality_distribution = vec![0u64; PHRED_BUCKETS];
    let mut quality_by_length_bin: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut per_tile: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut kmer_total: HashMap<[u8; 4], u64> = HashMap::with_capacity(256);
    let mut hll_pe: HyperLogLog = HyperLogLog::new(0.01);
    let mut overrep_map_pe: HashMap<[u8; 50], u64> = HashMap::with_capacity(4096);
    let mut overrep_sampled_pe = 0usize;
    let mut flush_counter = 0u64;
    let mut total_flushed = 0u64;
    let mut gc_distribution = [0u64; 101];
    let mut dup_sample: HashMap<u64, u32> = HashMap::with_capacity(8192);
    let mut dup_sample_count = 0usize;

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

        let per_tile_enabled = config.flags.per_tile_quality;
        let acc: BatchAccum = batch
            .par_iter()
            .fold(BatchAccum::default, move |mut acc, (r1, r2)| {
                for rec in [r1, r2] {
                    let Some((eff_seq, eff_qual, adapter_name)) =
                        trim_record(&rec.seq, rec.qual.as_deref(), config)
                    else {
                        acc.reads_filtered += 1;
                        continue;
                    };
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
                    acc.process(eff_seq, eff_qual, tile_info, fp, adapter_name.is_some());
                }
                acc
            })
            .reduce(BatchAccum::default, BatchAccum::merge);

        if config.flags.duplication_check {
            for &fp in &acc.fingerprints {
                hll_pe.insert(&fp);
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
            &mut length_histogram, &mut quality_by_pos, &mut qual_hist_by_pos,
            &mut base_composition, &mut quality_distribution, &mut per_tile,
            &mut quality_by_length_bin,
        );
        reads_filtered += acc.reads_filtered;

        // GC per-read distribution (reuses already-computed per-batch data)
        for (dst, &src) in gc_distribution.iter_mut().zip(acc.gc_per_read.iter()) {
            *dst += src;
        }
        // Dup level sample (capped at 200k reads for histogram)
        if config.flags.duplication_check && dup_sample_count < OVERREP_SAMPLE {
            for &fp in &acc.fingerprints {
                if dup_sample_count >= OVERREP_SAMPLE { break; }
                *dup_sample.entry(fp).or_insert(0) += 1;
                dup_sample_count += 1;
            }
        }

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

    let dup_rate_pe = if config.flags.duplication_check && read_count > 0 {
        let unique_est = hll_pe.len() as u64;
        read_count.saturating_sub(unique_est) as f64 / read_count as f64 * 100.0
    } else { 0.0 };
    let overrep_seqs_pe = finish_overrep(overrep_map_pe, overrep_sampled_pe);

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
        cur.reads_filtered        = reads_filtered;
        cur.min_length            = min_length;
        cur.max_length            = max_length;
        cur.bytes_processed       = r1_size + r2_size;
        cur.length_histogram      = length_histogram;
        cur.quality_by_position   = quality_by_pos;
        cur.qual_hist_by_position = qual_hist_by_pos;
        cur.kmer_counts           = kmer_total;
        cur.dup_rate_pct          = dup_rate_pe;
        cur.per_tile_quality      = per_tile;
        cur.base_composition      = base_composition;
        cur.quality_distribution  = quality_distribution;
        cur.quality_by_length_bin = quality_by_length_bin;
        cur.overrepresented_sequences = overrep_seqs_pe;
        cur.gc_distribution = gc_distribution.to_vec();
        cur.dup_level_histogram = compute_dup_histogram(&dup_sample);
        cur.module_status = compute_module_status(cur);
    }

    let elapsed = s.elapsed_secs();
    let pairs = read_count / 2;
    let mbps = if elapsed > 0.0 { bytes as f64 / 1_048_576.0 / elapsed } else { 0.0 };
    let gc_pct = if total_bases > 0 { gc_count as f64 / total_bases as f64 * 100.0 } else { 0.0 };
    let avg_q = if quality_bases > 0 { quality_sum as f64 / quality_bases as f64 } else { 0.0 };

    s.log(LogLevel::Success, format!(
        "Done: {} pairs | GC {:.1}% | Q{:.1} | dup {:.1}% | {:.1} MB/s",
        format_number(pairs), gc_pct, avg_q, dup_rate_pe, mbps,
    ));

    if let Some(done) = s.current.take() {
        s.completed_files.push(done);
    }
}
