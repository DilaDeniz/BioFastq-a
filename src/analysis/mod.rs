mod adapters;
mod io;
mod kmers;
mod metrics;

pub use adapters::{find_adapter_pos, find_adapter_pos_with_custom, quality_trim_3p};

use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::{Arc, Mutex};

use needletail::parse_fastx_file;
use needletail::Sequence;

use crate::types::{
    format_number, FileStats, LogLevel, ProcessConfig, ProcessingStatus, SharedState,
    BATCH_SIZE, DUP_SAMPLE_SIZE, FLUSH_INTERVAL, MAX_QUAL_POSITION,
};

use self::io::{open_writer, write_fastq_record};
use self::kmers::count_kmers_parallel;
use self::metrics::{accum_read_stats, fingerprint};

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
    let mut length_histogram: HashMap<u64, u64> = HashMap::new();
    let mut quality_by_pos = vec![(0u64, 0u64); MAX_QUAL_POSITION];
    let mut per_tile: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut kmer_total: HashMap<[u8; 4], u64> = HashMap::new();
    let mut batch_seqs: Vec<Vec<u8>> = Vec::with_capacity(BATCH_SIZE);
    let mut base_composition = vec![[0u64; 5]; MAX_QUAL_POSITION];
    let mut quality_distribution = vec![0u64; 43];
    let mut dup_seen: HashSet<u64> = HashSet::with_capacity(DUP_SAMPLE_SIZE);
    let mut dup_sampled = 0u64;
    let mut total_flushed = 0u64;
    let mut flush_counter = 0u64;

    while let Some(result) = reader.next() {
        let record = match result {
            Ok(r) => r,
            Err(e) => {
                let mut s = state.lock().unwrap();
                if config.strict {
                    let msg = format!("Bad record: {}", e);
                    s.log(LogLevel::Error, msg.clone());
                    s.status = ProcessingStatus::Error(msg);
                    return;
                }
                s.log(LogLevel::Warning, format!("Skipping bad record: {}", e));
                continue;
            }
        };

        let id = record.id();
        let raw_seq = record.normalize(false);
        let raw_qual = record.qual();

        bytes += raw_seq.len() as u64 + id.len() as u64 + 10;
        if let Some(q) = raw_qual { bytes += q.len() as u64; }
        read_count += 1;

        if dup_sampled < DUP_SAMPLE_SIZE as u64 {
            dup_seen.insert(fingerprint(&raw_seq));
            dup_sampled += 1;
        }

        if let (Some(tile), Some(q)) = (parse_illumina_tile(id), raw_qual) {
            let e = per_tile.entry(tile).or_insert((0, 0));
            for &qb in q {
                e.0 += qb.saturating_sub(33) as u64;
                e.1 += 1;
            }
        }

        let (seq, qual): (&[u8], Option<&[u8]>) = if config.quality_trim_threshold > 0 {
            if let Some(q) = raw_qual {
                let cut = adapters::quality_trim_3p(q, config.quality_trim_threshold);
                (&raw_seq[..cut], Some(&q[..cut]))
            } else {
                (&raw_seq, raw_qual)
            }
        } else {
            (&raw_seq, raw_qual)
        };

        let adapter_pos = adapters::find_adapter_pos_with_custom(seq, &config.custom_adapters);

        if let Some(pos) = adapter_pos {
            adapter_hits += 1;
            if config.trim_output {
                let trimmed_seq = &seq[..pos];
                let trimmed_len = trimmed_seq.len() as u64;
                trimmed_bases_removed += seq.len() as u64 - trimmed_len;
                if trimmed_len >= config.min_length {
                    trimmed_reads += 1;
                    if let Some(ref mut w) = trim_writer {
                        let tq = qual.map(|q| &q[..pos]);
                        if let Err(e) = write_fastq_record(w, id, trimmed_seq, tq) {
                            state.lock().unwrap().log(LogLevel::Warning, format!("Trim write error: {}", e));
                        }
                    }
                }
                accum_read_stats(
                    trimmed_seq, qual.map(|q| &q[..pos]),
                    &mut total_bases, &mut gc_count, &mut min_length, &mut max_length,
                    &mut length_histogram, &mut quality_sum, &mut quality_bases,
                    &mut q20_bases, &mut q30_bases, &mut quality_by_pos,
                    &mut base_composition, &mut quality_distribution, &mut batch_seqs,
                );
            } else {
                accum_read_stats(
                    seq, qual,
                    &mut total_bases, &mut gc_count, &mut min_length, &mut max_length,
                    &mut length_histogram, &mut quality_sum, &mut quality_bases,
                    &mut q20_bases, &mut q30_bases, &mut quality_by_pos,
                    &mut base_composition, &mut quality_distribution, &mut batch_seqs,
                );
            }
        } else {
            if config.trim_output {
                if let Some(ref mut w) = trim_writer {
                    if let Err(e) = write_fastq_record(w, id, seq, qual) {
                        state.lock().unwrap().log(LogLevel::Warning, format!("Trim write error: {}", e));
                    }
                }
            }
            accum_read_stats(
                seq, qual,
                &mut total_bases, &mut gc_count, &mut min_length, &mut max_length,
                &mut length_histogram, &mut quality_sum, &mut quality_bases,
                &mut q20_bases, &mut q30_bases, &mut quality_by_pos,
                &mut base_composition, &mut quality_distribution, &mut batch_seqs,
            );
        }

        flush_counter += 1;
        if flush_counter % FLUSH_INTERVAL == 0 {
            if batch_seqs.len() >= BATCH_SIZE {
                let kmers = count_kmers_parallel(&batch_seqs);
                for (k, v) in kmers { *kmer_total.entry(k).or_insert(0) += v; }
                batch_seqs.clear();
            }
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
                cur.trimmed_reads = trimmed_reads;
                cur.trimmed_bases_removed = trimmed_bases_removed;
                cur.min_length = min_length;
                cur.max_length = max_length;
                cur.bytes_processed = bytes.min(file_size);
            }
            total_flushed += FLUSH_INTERVAL;
            if total_flushed % 100_000 == 0 {
                let elapsed = s.elapsed_secs();
                let rps = if elapsed > 0.0 { read_count as f64 / elapsed } else { 0.0 };
                let mbps = if elapsed > 0.0 { bytes as f64 / 1_048_576.0 / elapsed } else { 0.0 };
                s.log(LogLevel::Info, format!(
                    "{} reads  {:.0} reads/s  {:.1} MB/s",
                    format_number(read_count), rps, mbps,
                ));
            }
        }
    }

    if !batch_seqs.is_empty() {
        let kmers = count_kmers_parallel(&batch_seqs);
        for (k, v) in kmers { *kmer_total.entry(k).or_insert(0) += v; }
    }
    drop(trim_writer);

    let dup_rate = if dup_sampled > 0 {
        (1.0 - dup_seen.len() as f64 / dup_sampled as f64) * 100.0
    } else { 0.0 };

    let mut s = state.lock().unwrap();
    if let Some(ref mut cur) = s.current {
        cur.read_count       = read_count;
        cur.total_bases      = total_bases;
        cur.gc_count         = gc_count;
        cur.quality_sum      = quality_sum;
        cur.quality_bases    = quality_bases;
        cur.q20_bases        = q20_bases;
        cur.q30_bases        = q30_bases;
        cur.adapter_hits     = adapter_hits;
        cur.trimmed_reads    = trimmed_reads;
        cur.trimmed_bases_removed = trimmed_bases_removed;
        cur.min_length       = min_length;
        cur.max_length       = max_length;
        cur.bytes_processed  = file_size;
        cur.length_histogram = length_histogram;
        cur.quality_by_position = quality_by_pos;
        cur.kmer_counts      = kmer_total;
        cur.dup_rate_pct     = dup_rate;
        cur.per_tile_quality = per_tile;
        cur.base_composition = base_composition;
        cur.quality_distribution = quality_distribution;
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

fn parse_illumina_tile(id: &[u8]) -> Option<u32> {
    let s = std::str::from_utf8(id).ok()?;
    let s = s.split_ascii_whitespace().next()?;
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() >= 6 { parts[4].parse().ok() } else { None }
}
