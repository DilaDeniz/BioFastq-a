use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use needletail::parse_fastx_file;
use needletail::Sequence;
use rayon::prelude::*;

use crate::types::{
    format_number, FileStats, LogLevel, ProcessingStatus, SharedState, ADAPTER_MATCH_LEN,
    ADAPTERS, BATCH_SIZE, FLUSH_INTERVAL, MAX_QUAL_POSITION,
};

/// Entry point: process a list of files sequentially, updating shared state.
pub fn process_files(paths: Vec<String>, state: Arc<Mutex<SharedState>>) {
    let n = paths.len();
    for (idx, path) in paths.into_iter().enumerate() {
        {
            let mut s = state.lock().unwrap();
            s.current_file_idx = idx;
            let short = std::path::Path::new(&path)
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or(&path)
                .to_string();
            s.log(
                LogLevel::Info,
                format!("Starting file {}/{}: {}", idx + 1, n, short),
            );
        }
        process_single_file(path, Arc::clone(&state));

        // Stop processing further files if an error occurred
        {
            let s = state.lock().unwrap();
            if matches!(s.status, ProcessingStatus::Error(_)) {
                return;
            }
        }
    }

    let mut s = state.lock().unwrap();
    s.status = ProcessingStatus::Completed;
    s.current = None;
    let total: u64 = s.completed_files.iter().map(|f| f.read_count).sum();
    s.log(
        LogLevel::Success,
        format!(
            "All {} file(s) complete — {} total reads",
            n,
            format_number(total)
        ),
    );
}

fn process_single_file(file_path: String, state: Arc<Mutex<SharedState>>) {
    let file_size = std::fs::metadata(&file_path)
        .map(|m| m.len())
        .unwrap_or(0);

    {
        let mut s = state.lock().unwrap();
        s.current = Some(FileStats::new(file_path.clone(), file_size));
        let size_gb = file_size as f64 / 1_073_741_824.0;
        let cpus = rayon::current_num_threads();
        s.log(
            LogLevel::Info,
            format!(
                "Opening {} ({:.2} GB) with {} CPU threads",
                &file_path, size_gb, cpus
            ),
        );
    }

    let reader = match parse_fastx_file(&file_path) {
        Ok(r) => r,
        Err(e) => {
            let mut s = state.lock().unwrap();
            let msg = format!("Cannot open {}: {}", file_path, e);
            s.log(LogLevel::Error, msg.clone());
            s.status = ProcessingStatus::Error(msg);
            return;
        }
    };

    let mut reader = reader;

    // --- Local running totals (SET into shared state on each flush) ---
    let mut read_count: u64 = 0;
    let mut total_bases: u64 = 0;
    let mut gc_count: u64 = 0;
    let mut quality_sum: u64 = 0;
    let mut quality_bases: u64 = 0;
    let mut q20_reads: u64 = 0;
    let mut q30_reads: u64 = 0;
    let mut adapter_hits: u64 = 0;
    let mut min_length: u64 = u64::MAX;
    let mut max_length: u64 = 0;
    let mut bytes: u64 = 0;
    let mut length_histogram: HashMap<u64, u64> = HashMap::new();
    let mut quality_by_pos: Vec<(u64, u64)> = vec![(0, 0); MAX_QUAL_POSITION];
    let mut kmer_total: HashMap<[u8; 4], u64> = HashMap::new();
    let mut batch_seqs: Vec<Vec<u8>> = Vec::with_capacity(BATCH_SIZE);
    let mut total_flushed: u64 = 0;
    let mut flush_counter: u64 = 0;

    while let Some(result) = reader.next() {
        let record = match result {
            Ok(r) => r,
            Err(e) => {
                let mut s = state.lock().unwrap();
                s.log(LogLevel::Warning, format!("Skipping bad record: {}", e));
                continue;
            }
        };

        let seq = record.normalize(false);
        let seq_len = seq.len() as u64;

        read_count += 1;
        total_bases += seq_len;
        bytes += seq_len + record.id().len() as u64 + 10;

        // Length tracking
        if seq_len < min_length {
            min_length = seq_len;
        }
        if seq_len > max_length {
            max_length = seq_len;
        }
        *length_histogram.entry(seq_len).or_insert(0) += 1;

        // GC content (uppercase only; normalize(false) leaves case intact)
        gc_count += seq
            .iter()
            .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
            .count() as u64;

        // Quality scores (FASTQ only)
        if let Some(qual) = record.qual() {
            bytes += qual.len() as u64;
            let mut read_q_sum: u64 = 0;
            let read_q_len = qual.len() as u64;

            for (pos, &q) in qual.iter().enumerate() {
                let phred = q.saturating_sub(33) as u64;
                quality_sum += phred;
                quality_bases += 1;
                read_q_sum += phred;
                if pos < MAX_QUAL_POSITION {
                    quality_by_pos[pos].0 += phred;
                    quality_by_pos[pos].1 += 1;
                }
            }

            if read_q_len > 0 {
                let read_avg = read_q_sum as f64 / read_q_len as f64;
                if read_avg >= 20.0 {
                    q20_reads += 1;
                }
                if read_avg >= 30.0 {
                    q30_reads += 1;
                }
            }
        }

        // Adapter detection via exact prefix match
        if has_adapter(&seq) {
            adapter_hits += 1;
        }

        // Collect sequence for batched parallel k-mer counting
        if seq.len() >= 4 {
            batch_seqs.push(seq.to_vec());
        }

        flush_counter += 1;
        if flush_counter.is_multiple_of(FLUSH_INTERVAL) {
            // Count k-mers for the current batch if it's large enough
            if batch_seqs.len() >= BATCH_SIZE {
                let kmers = count_kmers_parallel(&batch_seqs);
                for (k, v) in kmers {
                    *kmer_total.entry(k).or_insert(0) += v;
                }
                batch_seqs.clear();
            }

            // Push running totals into shared state (SET, not accumulate)
            let mut s = state.lock().unwrap();
            if let Some(ref mut cur) = s.current {
                cur.read_count = read_count;
                cur.total_bases = total_bases;
                cur.gc_count = gc_count;
                cur.quality_sum = quality_sum;
                cur.quality_bases = quality_bases;
                cur.q20_reads = q20_reads;
                cur.q30_reads = q30_reads;
                cur.adapter_hits = adapter_hits;
                cur.min_length = min_length;
                cur.max_length = max_length;
                cur.bytes_processed = bytes.min(file_size);
            }

            total_flushed += FLUSH_INTERVAL;
            if total_flushed.is_multiple_of(100_000) {
                let elapsed = s.elapsed_secs();
                let rps = if elapsed > 0.0 {
                    read_count as f64 / elapsed
                } else {
                    0.0
                };
                s.log(
                    LogLevel::Info,
                    format!(
                        "{} reads processed  ({:.0} reads/sec)",
                        format_number(read_count),
                        rps
                    ),
                );
            }
        }
    }

    // Process remaining batch
    if !batch_seqs.is_empty() {
        let kmers = count_kmers_parallel(&batch_seqs);
        for (k, v) in kmers {
            *kmer_total.entry(k).or_insert(0) += v;
        }
    }

    // Final state update
    let gc_pct = if total_bases > 0 {
        gc_count as f64 / total_bases as f64 * 100.0
    } else {
        0.0
    };
    let avg_q = if quality_bases > 0 {
        quality_sum as f64 / quality_bases as f64
    } else {
        0.0
    };
    let adapter_pct = if read_count > 0 {
        adapter_hits as f64 / read_count as f64 * 100.0
    } else {
        0.0
    };

    {
        let mut s = state.lock().unwrap();
        if let Some(ref mut cur) = s.current {
            cur.read_count = read_count;
            cur.total_bases = total_bases;
            cur.gc_count = gc_count;
            cur.quality_sum = quality_sum;
            cur.quality_bases = quality_bases;
            cur.q20_reads = q20_reads;
            cur.q30_reads = q30_reads;
            cur.adapter_hits = adapter_hits;
            cur.min_length = min_length;
            cur.max_length = max_length;
            cur.bytes_processed = file_size;
            cur.length_histogram = length_histogram;
            cur.quality_by_position = quality_by_pos;
            cur.kmer_counts = kmer_total;
        }

        s.log(
            LogLevel::Success,
            format!(
                "Done: {} reads | GC {:.1}% | Avg Q{:.1} | {:.2}% adapter",
                format_number(read_count),
                gc_pct,
                avg_q,
                adapter_pct
            ),
        );

        // Move current file to completed list
        if let Some(done) = s.current.take() {
            s.completed_files.push(done);
        }
    }
}

/// Returns true if any known adapter sequence is found in the read.
fn has_adapter(seq: &[u8]) -> bool {
    ADAPTERS.iter().any(|(_, adapter)| {
        let search_len = adapter.len().min(ADAPTER_MATCH_LEN);
        let search = &adapter[..search_len];
        seq.windows(search_len).any(|w| w == search)
    })
}

/// Count 4-mers across a batch of sequences using rayon parallel fold/reduce.
fn count_kmers_parallel(sequences: &[Vec<u8>]) -> HashMap<[u8; 4], u64> {
    sequences
        .par_iter()
        .fold(HashMap::new, |mut map: HashMap<[u8; 4], u64>, seq| {
            for window in seq.windows(4) {
                if window
                    .iter()
                    .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
                {
                    let key = [window[0], window[1], window[2], window[3]];
                    *map.entry(key).or_insert(0) += 1;
                }
            }
            map
        })
        .reduce(HashMap::new, |mut a, b| {
            for (k, v) in b {
                *a.entry(k).or_insert(0) += v;
            }
            a
        })
}
