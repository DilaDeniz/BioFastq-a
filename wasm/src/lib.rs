mod metrics;
mod report;
mod types;

use std::collections::HashMap;
use wasm_bindgen::prelude::*;

use aho_corasick::AhoCorasick;
use metrics::{estimate_dup_rate, fingerprint, kmer_table_to_map, BatchAccum};
use types::{
    FileStats, OverrepSeq, QcStatus, ADAPTER_MATCH_LEN, ADAPTERS, OVERREP_SAMPLE,
};

// ---------------------------------------------------------------------------
// FASTQ / FASTA sequential parser
// ---------------------------------------------------------------------------

/// Single parsed record with borrowed slices into the input buffer.
struct Record<'a> {
    seq: &'a [u8],
    qual: Option<&'a [u8]>,
}

/// Iterate over FASTQ records in `data` (4-line format).
/// Also handles plain FASTA (2-line: >header + sequence) when `data` starts with '>'.
struct FastqIter<'a> {
    data: &'a [u8],
    pos: usize,
    is_fasta: bool,
}

impl<'a> FastqIter<'a> {
    fn new(data: &'a [u8]) -> Self {
        let is_fasta = data.first().copied() == Some(b'>');
        Self { data, pos: 0, is_fasta }
    }

    /// Return the next newline-terminated line as a slice, advancing `pos`.
    fn next_line(&mut self) -> Option<&'a [u8]> {
        if self.pos >= self.data.len() {
            return None;
        }
        let start = self.pos;
        // Use memchr for fast newline search
        let end = memchr::memchr(b'\n', &self.data[self.pos..])
            .map(|i| self.pos + i)
            .unwrap_or(self.data.len());
        self.pos = end + 1;
        // Strip trailing '\r' for Windows-style line endings
        let line = &self.data[start..end];
        if line.last() == Some(&b'\r') {
            Some(&line[..line.len() - 1])
        } else {
            Some(line)
        }
    }
}

impl<'a> Iterator for FastqIter<'a> {
    type Item = Record<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.is_fasta {
                // Skip header line (starts with '>')
                let header = self.next_line()?;
                if header.is_empty() || !header.starts_with(b">") {
                    // Skip blank lines or unexpected content
                    continue;
                }
                let seq = self.next_line()?;
                if seq.is_empty() { continue; }
                return Some(Record { seq, qual: None });
            } else {
                // FASTQ: @id, seq, +, qual
                let id_line = self.next_line()?;
                if id_line.is_empty() { continue; }
                if !id_line.starts_with(b"@") { continue; }
                let seq = self.next_line()?;
                if seq.is_empty() { continue; }
                let _plus = self.next_line()?;
                let qual = self.next_line()?;
                return Some(Record { seq, qual: Some(qual) });
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Adapter detection via Aho-Corasick
// ---------------------------------------------------------------------------

fn build_adapter_automaton() -> AhoCorasick {
    let patterns: Vec<&[u8]> = ADAPTERS.iter().map(|(_, seq)| *seq).collect();
    AhoCorasick::builder()
        .ascii_case_insensitive(false)
        .build(patterns)
        .expect("valid adapter patterns")
}

fn check_adapter(seq: &[u8], ac: &AhoCorasick) -> bool {
    // Quick prefix-match check first (matches original ADAPTER_MATCH_LEN logic)
    for (_, pat) in ADAPTERS {
        let cmp_len = ADAPTER_MATCH_LEN.min(pat.len()).min(seq.len());
        if cmp_len > 0 && seq[..cmp_len] == pat[..cmp_len] {
            return true;
        }
    }
    // Full Aho-Corasick scan for internal matches
    ac.find(seq).is_some()
}

// ---------------------------------------------------------------------------
// Overrepresented sequence detection
// ---------------------------------------------------------------------------

fn detect_overrep(seqs: &[[u8; 50]], total_sampled: u64) -> Vec<OverrepSeq> {
    if seqs.is_empty() || total_sampled == 0 {
        return Vec::new();
    }

    let threshold = (total_sampled as f64 * 0.001) as u64; // 0.1%

    // Count occurrences of each prefix
    let mut counts: HashMap<[u8; 50], u64> = HashMap::with_capacity(1024);
    for key in seqs {
        *counts.entry(*key).or_insert(0) += 1;
    }

    // Filter and sort
    let mut overrep: Vec<_> = counts.into_iter()
        .filter(|(_, cnt)| *cnt >= threshold.max(2))
        .collect();
    overrep.sort_by(|a, b| b.1.cmp(&a.1));
    overrep.truncate(20);

    // Build adapter prefix for source matching
    let adapter_ac = AhoCorasick::builder()
        .build(ADAPTERS.iter().map(|(_, s)| s).copied())
        .unwrap();

    overrep.into_iter().map(|(key, count)| {
        // Trim null bytes from the 50-byte key
        let len = key.iter().position(|&b| b == 0).unwrap_or(50);
        let sequence = String::from_utf8_lossy(&key[..len]).to_string();
        let percentage = count as f64 / total_sampled as f64 * 100.0;
        let possible_source = if adapter_ac.find(key[..len].as_ref()).is_some() {
            "Adapter sequence".to_string()
        } else if sequence.chars().all(|c| c == 'A') {
            "Poly-A".to_string()
        } else if sequence.chars().all(|c| c == 'T') {
            "Poly-T".to_string()
        } else if sequence.chars().all(|c| c == 'G') {
            "Poly-G (NovaSeq signal loss)".to_string()
        } else if sequence.chars().all(|c| c == 'C') {
            "Poly-C".to_string()
        } else {
            "No hit".to_string()
        };
        OverrepSeq { sequence, count, percentage, possible_source }
    }).collect()
}

// ---------------------------------------------------------------------------
// Module status (FastQC-style traffic lights)
// ---------------------------------------------------------------------------

fn compute_module_status(f: &FileStats) -> Vec<(String, QcStatus)> {
    let mut v = Vec::new();

    // Per-base sequence quality
    let mean_q = f.avg_quality();
    v.push(("Per-base quality".to_string(),
        if mean_q >= 28.0 { QcStatus::Pass } else if mean_q >= 20.0 { QcStatus::Warn } else { QcStatus::Fail }));

    // Per-read quality
    let q30 = f.q30_pct();
    v.push(("Per-read quality".to_string(),
        if q30 >= 70.0 { QcStatus::Pass } else if q30 >= 50.0 { QcStatus::Warn } else { QcStatus::Fail }));

    // GC content
    let gc = f.gc_content();
    v.push(("GC content".to_string(),
        if (35.0..=65.0).contains(&gc) { QcStatus::Pass }
        else if (25.0..=75.0).contains(&gc) { QcStatus::Warn }
        else { QcStatus::Fail }));

    // Adapter content
    let ap = f.adapter_pct();
    v.push(("Adapter content".to_string(),
        if ap < 1.0 { QcStatus::Pass } else if ap < 10.0 { QcStatus::Warn } else { QcStatus::Fail }));

    // Sequence duplication
    let dup = f.dup_rate_pct;
    v.push(("Sequence duplication".to_string(),
        if dup < 5.0 { QcStatus::Pass } else if dup < 20.0 { QcStatus::Warn } else { QcStatus::Fail }));

    // Overrepresented sequences
    v.push(("Overrepresented seqs".to_string(),
        if f.overrepresented_sequences.is_empty() { QcStatus::Pass }
        else if f.overrepresented_sequences.len() <= 3 { QcStatus::Warn }
        else { QcStatus::Fail }));

    // N content
    let high_n = f.n_content_per_position().iter().any(|&n| n > 5.0);
    v.push(("N content".to_string(), if high_n { QcStatus::Fail } else { QcStatus::Pass }));

    v
}

// ---------------------------------------------------------------------------
// Long-read auto-detection
// ---------------------------------------------------------------------------

fn is_long_read(acc: &BatchAccum) -> bool {
    if acc.read_count == 0 { return false; }
    let total_reads = acc.read_count;
    // Compute median length from the length histogram
    let half = total_reads / 2;
    let mut cum = 0u64;
    for (len, &cnt) in acc.length_hist.iter().enumerate() {
        cum += cnt;
        if cum >= half {
            return len > 1000;
        }
    }
    false
}

// ---------------------------------------------------------------------------
// Main Wasm entry point
// ---------------------------------------------------------------------------

/// Decompress gzip input (detected via magic bytes `1f 8b`), otherwise pass through unchanged.
fn decompress_if_gzip(data: &[u8]) -> std::borrow::Cow<'_, [u8]> {
    if data.len() >= 2 && data[0] == 0x1f && data[1] == 0x8b {
        use std::io::Read;
        let mut out = Vec::new();
        let _ = flate2::read::GzDecoder::new(data).read_to_end(&mut out);
        std::borrow::Cow::Owned(out)
    } else {
        std::borrow::Cow::Borrowed(data)
    }
}

/// Analyse a FASTQ/FASTA file provided as raw bytes (gzip-compressed or plain).
/// Returns a complete HTML report as a String.
/// `timestamp` is a JS-formatted date string (e.g. from `new Date().toLocaleString()`).
/// `elapsed_ms` is wall-clock time in milliseconds (passed from JS `performance.now()`).
#[wasm_bindgen]
pub fn analyze(data: &[u8], filename: &str, timestamp: &str, elapsed_ms: f64) -> String {
    let compressed_size = data.len() as u64;
    let decompressed = decompress_if_gzip(data);
    let data: &[u8] = decompressed.as_ref();

    let ac = build_adapter_automaton();

    let mut acc = BatchAccum::default();

    for record in FastqIter::new(data) {
        let fp = fingerprint(record.seq);
        let adapter_hit = check_adapter(record.seq, &ac);
        acc.process(record.seq, record.qual, fp, adapter_hit);
    }

    finalize_to_filestats(acc, filename, compressed_size, timestamp, elapsed_ms)
}

/// Version exposed to JavaScript for capability detection.
#[wasm_bindgen]
pub fn version() -> String {
    "2.3.1-wasm".to_string()
}

// ---------------------------------------------------------------------------
// Finalize BatchAccum → FileStats → HTML
// ---------------------------------------------------------------------------

fn finalize_to_filestats(
    acc: BatchAccum,
    filename: &str,
    file_size: u64,
    timestamp: &str,
    elapsed_ms: f64,
) -> String {
    let long_read = is_long_read(&acc);
    let total_sampled = acc.read_count.min(OVERREP_SAMPLE as u64);

    // Duplication estimate from fingerprints
    let (dup_rate_pct, dup_level_histogram) = estimate_dup_rate(&acc.fingerprints);

    // Overrepresented sequences
    let overrep = detect_overrep(&acc.overrep_seqs, total_sampled);

    // Build FileStats
    let mut f = FileStats::new(filename.to_string(), file_size);
    f.read_count = acc.read_count;
    f.total_bases = acc.total_bases;
    f.gc_count = acc.gc_count;
    f.quality_sum = acc.quality_sum;
    f.quality_bases = acc.quality_bases;
    f.q20_bases = acc.q20_bases;
    f.q30_bases = acc.q30_bases;
    f.adapter_hits = acc.adapter_hits;
    f.min_length = if acc.min_length == u64::MAX { 0 } else { acc.min_length };
    f.max_length = acc.max_length;
    f.bytes_processed = file_size;
    f.dup_rate_pct = dup_rate_pct;
    f.dup_level_histogram = dup_level_histogram;
    f.overrepresented_sequences = overrep;
    f.long_read_mode = long_read;
    f.qual_vs_length = acc.qual_vs_length;
    f.per_tile_quality = acc.per_tile;

    // Length histogram: Vec<u64> → HashMap<u64, u64>
    f.length_histogram = acc.length_hist.iter().enumerate()
        .filter(|(_, &v)| v > 0)
        .map(|(i, &v)| (i as u64, v))
        .collect();

    // Quality by position
    let end_pos = acc.active_pos;
    for i in 0..end_pos {
        f.quality_by_position[i] = acc.quality_by_pos[i];
        // Convert u32 per-position histogram to u64
        for j in 0..43 {
            f.qual_hist_by_position[i][j] = acc.qual_hist_by_pos[i][j] as u64;
        }
        f.base_composition[i] = acc.base_composition[i];
    }

    // Quality distribution
    f.quality_distribution = acc.quality_distribution;

    // GC distribution (per-read histogram)
    f.gc_distribution = acc.gc_per_read.to_vec();

    // Quality by length bin → HashMap
    f.quality_by_length_bin = acc.quality_by_len_bin.iter().enumerate()
        .filter(|(_, t)| t.1 > 0)
        .map(|(i, &t)| (i as u32, t))
        .collect();

    // K-mer counts
    f.kmer_counts = kmer_table_to_map(&acc.kmer_table);

    // Module status
    f.module_status = compute_module_status(&f);

    report::generate_html(&f, timestamp, elapsed_ms)
}

// ---------------------------------------------------------------------------
// Panic hook (dev builds only)
// ---------------------------------------------------------------------------
#[cfg(feature = "console_error_panic_hook")]
pub fn set_panic_hook() {
    console_error_panic_hook::set_once();
}
