use std::collections::{HashMap, HashSet};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_file;
use needletail::Sequence;
use rayon::prelude::*;

use crate::types::{
    format_number, FileStats, LogLevel, ProcessConfig, ProcessingStatus, SharedState,
    ADAPTER_MATCH_LEN, ADAPTERS, BATCH_SIZE, DUP_SAMPLE_SIZE, FLUSH_INTERVAL, MAX_QUAL_POSITION,
};

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
            s.log(
                LogLevel::Info,
                format!("Starting file {}/{}: {}", idx + 1, n, short),
            );
        }
        process_single_file(path, Arc::clone(&state), &config);

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

// ---------------------------------------------------------------------------
// Per-file processing
// ---------------------------------------------------------------------------

fn process_single_file(file_path: String, state: Arc<Mutex<SharedState>>, config: &ProcessConfig) {
    let file_size = std::fs::metadata(&file_path)
        .map(|m| m.len())
        .unwrap_or(0);

    // Determine trim output path for this file
    let trim_path: Option<String> = if config.trim_output {
        let stem = Path::new(&file_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");
        // Strip .fastq from stem if present (handles .fastq.gz -> stem is .fastq)
        let stem = stem.trim_end_matches(".fastq");
        Some(format!("{}/{}_trimmed.fastq.gz", config.output_dir, stem))
    } else {
        None
    };

    {
        let mut s = state.lock().unwrap();
        let mut stats = FileStats::new(file_path.clone(), file_size);
        stats.trim_output_path = trim_path.clone();
        s.current = Some(stats);
        let size_gb = file_size as f64 / 1_073_741_824.0;
        let cpus = rayon::current_num_threads();
        let trim_msg = trim_path
            .as_deref()
            .map(|p| format!("  trim→ {}", p))
            .unwrap_or_default();
        s.log(
            LogLevel::Info,
            format!(
                "Opening {} ({:.2} GB) [{} threads]{}",
                &file_path, size_gb, cpus, trim_msg
            ),
        );
    }

    // Open reader
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

    // Open trim writer if requested
    let mut trim_writer: Option<Box<dyn Write>> = if let Some(ref tp) = trim_path {
        match open_writer(tp) {
            Ok(w) => Some(w),
            Err(e) => {
                let mut s = state.lock().unwrap();
                let msg = format!("Cannot open trim output {}: {}", tp, e);
                s.log(LogLevel::Error, msg.clone());
                s.status = ProcessingStatus::Error(msg);
                return;
            }
        }
    } else {
        None
    };

    let mut reader = reader;

    // --- Local running totals ---
    let mut read_count: u64 = 0;
    let mut total_bases: u64 = 0;
    let mut gc_count: u64 = 0;
    let mut quality_sum: u64 = 0;
    let mut quality_bases: u64 = 0;
    let mut q20_reads: u64 = 0;
    let mut q30_reads: u64 = 0;
    let mut adapter_hits: u64 = 0;
    let mut trimmed_reads: u64 = 0;
    let mut trimmed_bases_removed: u64 = 0;
    let mut min_length: u64 = u64::MAX;
    let mut max_length: u64 = 0;
    let mut bytes: u64 = 0;
    let mut length_histogram: HashMap<u64, u64> = HashMap::new();
    let mut quality_by_pos: Vec<(u64, u64)> = vec![(0, 0); MAX_QUAL_POSITION];
    let mut per_tile: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut kmer_total: HashMap<[u8; 4], u64> = HashMap::new();
    let mut batch_seqs: Vec<Vec<u8>> = Vec::with_capacity(BATCH_SIZE);
    // Duplication estimation
    let mut dup_seen: HashSet<u64> = HashSet::with_capacity(DUP_SAMPLE_SIZE);
    let mut dup_sampled: u64 = 0;
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

        let id = record.id();
        let seq = record.normalize(false);
        let qual = record.qual();
        let orig_len = seq.len() as u64;

        bytes += orig_len + id.len() as u64 + 10;
        if let Some(q) = qual {
            bytes += q.len() as u64;
        }

        read_count += 1;

        // --- Duplication fingerprint (first DUP_SAMPLE_SIZE reads) ---
        if dup_sampled < DUP_SAMPLE_SIZE as u64 {
            let h = fingerprint(&seq);
            dup_seen.insert(h);
            dup_sampled += 1;
        }

        // --- Per-tile quality (Illumina header parsing) ---
        if let Some(tile) = parse_illumina_tile(id) {
            if let Some(q) = qual {
                let tile_entry = per_tile.entry(tile).or_insert((0, 0));
                for &qb in q {
                    tile_entry.0 += qb.saturating_sub(33) as u64;
                    tile_entry.1 += 1;
                }
            }
        }

        // --- Adapter detection & trimming ---
        let trim_pos = find_adapter_pos(&seq);

        if let Some(pos) = trim_pos {
            adapter_hits += 1;

            if config.trim_output {
                let trimmed_seq = &seq[..pos];
                let trimmed_len = trimmed_seq.len() as u64;
                trimmed_bases_removed += orig_len - trimmed_len;

                if trimmed_len >= config.min_length {
                    trimmed_reads += 1;
                    if let Some(ref mut w) = trim_writer {
                        let trimmed_qual = qual.map(|q| &q[..pos]);
                        if let Err(e) = write_fastq_record(w, id, trimmed_seq, trimmed_qual) {
                            let mut s = state.lock().unwrap();
                            s.log(LogLevel::Warning, format!("Trim write error: {}", e));
                        }
                    }
                }

                // Stats use the TRIMMED length
                let seq_len = trimmed_len;
                total_bases += seq_len;
                gc_count += trimmed_seq
                    .iter()
                    .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
                    .count() as u64;
                if seq_len < min_length { min_length = seq_len; }
                if seq_len > max_length { max_length = seq_len; }
                *length_histogram.entry(seq_len).or_insert(0) += 1;
                if let Some(q) = qual {
                    let q_trimmed = &q[..pos];
                    accum_quality(
                        q_trimmed,
                        &mut quality_sum,
                        &mut quality_bases,
                        &mut q20_reads,
                        &mut q30_reads,
                        &mut quality_by_pos,
                    );
                }
                if seq_len >= 4 {
                    batch_seqs.push(trimmed_seq.to_vec());
                }
            } else {
                // Not trimming — use full read for stats but flag the adapter
                accum_read_stats(
                    &seq,
                    qual,
                    orig_len,
                    &mut total_bases,
                    &mut gc_count,
                    &mut min_length,
                    &mut max_length,
                    &mut length_histogram,
                    &mut quality_sum,
                    &mut quality_bases,
                    &mut q20_reads,
                    &mut q30_reads,
                    &mut quality_by_pos,
                    &mut batch_seqs,
                );
            }
        } else {
            // No adapter found
            if config.trim_output {
                // Write unmodified read to trim output
                if let Some(ref mut w) = trim_writer {
                    if let Err(e) = write_fastq_record(w, id, &seq, qual) {
                        let mut s = state.lock().unwrap();
                        s.log(LogLevel::Warning, format!("Trim write error: {}", e));
                    }
                }
            }
            accum_read_stats(
                &seq,
                qual,
                orig_len,
                &mut total_bases,
                &mut gc_count,
                &mut min_length,
                &mut max_length,
                &mut length_histogram,
                &mut quality_sum,
                &mut quality_bases,
                &mut q20_reads,
                &mut q30_reads,
                &mut quality_by_pos,
                &mut batch_seqs,
            );
        }

        flush_counter += 1;
        if flush_counter.is_multiple_of(FLUSH_INTERVAL) {
            if batch_seqs.len() >= BATCH_SIZE {
                let kmers = count_kmers_parallel(&batch_seqs);
                for (k, v) in kmers {
                    *kmer_total.entry(k).or_insert(0) += v;
                }
                batch_seqs.clear();
            }

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
                cur.trimmed_reads = trimmed_reads;
                cur.trimmed_bases_removed = trimmed_bases_removed;
                cur.min_length = min_length;
                cur.max_length = max_length;
                cur.bytes_processed = bytes.min(file_size);
            }

            total_flushed += FLUSH_INTERVAL;
            if total_flushed.is_multiple_of(100_000) {
                let elapsed = s.elapsed_secs();
                let rps = if elapsed > 0.0 { read_count as f64 / elapsed } else { 0.0 };
                let mbps = if elapsed > 0.0 { bytes as f64 / 1_048_576.0 / elapsed } else { 0.0 };
                s.log(
                    LogLevel::Info,
                    format!(
                        "{} reads  {:.0} reads/s  {:.1} MB/s",
                        format_number(read_count), rps, mbps
                    ),
                );
            }
        }
    }

    // Flush remaining k-mer batch
    if !batch_seqs.is_empty() {
        let kmers = count_kmers_parallel(&batch_seqs);
        for (k, v) in kmers {
            *kmer_total.entry(k).or_insert(0) += v;
        }
    }

    // Flush trim writer
    drop(trim_writer);

    // Compute duplication estimate
    let dup_rate = if dup_sampled > 0 {
        (1.0 - dup_seen.len() as f64 / dup_sampled as f64) * 100.0
    } else {
        0.0
    };

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
    let mbps = {
        // calc at end
        let s = state.lock().unwrap();
        let elapsed = s.elapsed_secs();
        if elapsed > 0.0 { bytes as f64 / 1_048_576.0 / elapsed } else { 0.0 }
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
            cur.trimmed_reads = trimmed_reads;
            cur.trimmed_bases_removed = trimmed_bases_removed;
            cur.min_length = min_length;
            cur.max_length = max_length;
            cur.bytes_processed = file_size;
            cur.length_histogram = length_histogram;
            cur.quality_by_position = quality_by_pos;
            cur.kmer_counts = kmer_total;
            cur.dup_rate_pct = dup_rate;
            cur.per_tile_quality = per_tile;
        }

        let trim_note = if config.trim_output {
            format!(" | trimmed {:.1}%", if read_count > 0 { trimmed_reads as f64 / read_count as f64 * 100.0 } else { 0.0 })
        } else {
            String::new()
        };
        let tile_count = {
            s.current.as_ref().map(|c| c.per_tile_quality.len()).unwrap_or(0)
        };
        let tile_note = if tile_count > 0 {
            format!(" | {} tiles", tile_count)
        } else {
            String::new()
        };

        s.log(
            LogLevel::Success,
            format!(
                "Done: {} reads | GC {:.1}% | Q{:.1} | {:.2}% adapter | dup ~{:.1}% | {:.1} MB/s{}{}",
                format_number(read_count),
                gc_pct,
                avg_q,
                adapter_pct,
                dup_rate,
                mbps,
                trim_note,
                tile_note,
            ),
        );

        if let Some(done) = s.current.take() {
            s.completed_files.push(done);
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Accumulate stats for a single read into local counters.
#[allow(clippy::too_many_arguments)]
fn accum_read_stats(
    seq: &[u8],
    qual: Option<&[u8]>,
    seq_len: u64,
    total_bases: &mut u64,
    gc_count: &mut u64,
    min_length: &mut u64,
    max_length: &mut u64,
    length_histogram: &mut HashMap<u64, u64>,
    quality_sum: &mut u64,
    quality_bases: &mut u64,
    q20_reads: &mut u64,
    q30_reads: &mut u64,
    quality_by_pos: &mut Vec<(u64, u64)>,
    batch_seqs: &mut Vec<Vec<u8>>,
) {
    *total_bases += seq_len;
    if seq_len < *min_length { *min_length = seq_len; }
    if seq_len > *max_length { *max_length = seq_len; }
    *length_histogram.entry(seq_len).or_insert(0) += 1;
    *gc_count += seq
        .iter()
        .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
        .count() as u64;
    if let Some(q) = qual {
        accum_quality(q, quality_sum, quality_bases, q20_reads, q30_reads, quality_by_pos);
    }
    if seq.len() >= 4 {
        batch_seqs.push(seq.to_vec());
    }
}

fn accum_quality(
    qual: &[u8],
    quality_sum: &mut u64,
    quality_bases: &mut u64,
    q20_reads: &mut u64,
    q30_reads: &mut u64,
    quality_by_pos: &mut Vec<(u64, u64)>,
) {
    let mut read_q_sum: u64 = 0;
    let read_q_len = qual.len() as u64;
    for (pos, &q) in qual.iter().enumerate() {
        let phred = q.saturating_sub(33) as u64;
        *quality_sum += phred;
        *quality_bases += 1;
        read_q_sum += phred;
        if pos < quality_by_pos.len() {
            quality_by_pos[pos].0 += phred;
            quality_by_pos[pos].1 += 1;
        }
    }
    if read_q_len > 0 {
        let avg = read_q_sum as f64 / read_q_len as f64;
        if avg >= 20.0 { *q20_reads += 1; }
        if avg >= 30.0 { *q30_reads += 1; }
    }
}

/// Hash the first 50 bytes of a sequence for duplication fingerprinting.
fn fingerprint(seq: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    let len = seq.len().min(50);
    seq[..len].hash(&mut hasher);
    hasher.finish()
}

/// Find the leftmost position where a known adapter begins.
pub fn find_adapter_pos(seq: &[u8]) -> Option<usize> {
    ADAPTERS
        .iter()
        .filter_map(|(_, adapter)| {
            let n = adapter.len().min(ADAPTER_MATCH_LEN);
            let needle = &adapter[..n];
            seq.windows(n).position(|w| w == needle)
        })
        .min()
}

/// Parse Illumina tile ID from a FASTQ header.
/// Expects format: @instrument:run:flowcell:lane:tile:x:y [...]
pub fn parse_illumina_tile(id: &[u8]) -> Option<u32> {
    let s = std::str::from_utf8(id).ok()?;
    // Take only the part before any whitespace
    let s = s.split_ascii_whitespace().next()?;
    let parts: Vec<&str> = s.split(':').collect();
    // Field index 4 (0-based) is the tile number in Illumina 1.8+ format
    if parts.len() >= 6 {
        parts[4].parse().ok()
    } else {
        None
    }
}

/// Open a file for writing; gzip-compress if path ends in `.gz`.
fn open_writer(path: &str) -> io::Result<Box<dyn Write>> {
    let file = std::fs::File::create(path)?;
    if path.ends_with(".gz") {
        Ok(Box::new(BufWriter::new(GzEncoder::new(
            file,
            Compression::default(),
        ))))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

/// Write one FASTQ record to `writer`.
fn write_fastq_record(
    writer: &mut dyn Write,
    id: &[u8],
    seq: &[u8],
    qual: Option<&[u8]>,
) -> io::Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    if let Some(q) = qual {
        writer.write_all(q)?;
    } else {
        // FASTA input — no quality; write placeholder so output is still valid
        let fake: Vec<u8> = vec![b'I'; seq.len()];
        writer.write_all(&fake)?;
    }
    writer.write_all(b"\n")?;
    Ok(())
}

/// Count 4-mers across a batch using rayon parallel fold/reduce.
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
