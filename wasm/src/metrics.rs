use std::collections::HashMap;
use xxhash_rust::xxh3::xxh3_64;
use crate::types::{PHRED_BUCKETS, MAX_QUAL_POSITION, OVERREP_SAMPLE};

pub(crate) static BASE_LUT: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0; t[b'a' as usize] = 0;
    t[b'C' as usize] = 1; t[b'c' as usize] = 1;
    t[b'G' as usize] = 2; t[b'g' as usize] = 2;
    t[b'T' as usize] = 3; t[b't' as usize] = 3;
    t
};

#[inline(always)]
fn base_to_2bit(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

#[inline(always)]
fn kmer4_index(seq: &[u8], pos: usize) -> Option<u8> {
    let a = base_to_2bit(seq[pos])?;
    let b = base_to_2bit(seq[pos + 1])?;
    let c = base_to_2bit(seq[pos + 2])?;
    let d = base_to_2bit(seq[pos + 3])?;
    Some((a << 6) | (b << 4) | (c << 2) | d)
}

/// Count G and C bases using memchr — works in Wasm (Wasm SIMD when available).
#[inline]
pub fn count_gc(seq: &[u8]) -> u64 {
    (memchr::memchr_iter(b'G', seq).count()
        + memchr::memchr_iter(b'C', seq).count()
        + memchr::memchr_iter(b'g', seq).count()
        + memchr::memchr_iter(b'c', seq).count()) as u64
}

pub struct BatchAccum {
    pub read_count: u64,
    pub total_bases: u64,
    pub gc_count: u64,
    pub quality_sum: u64,
    pub quality_bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
    pub adapter_hits: u64,
    pub min_length: u64,
    pub max_length: u64,
    pub length_hist: Vec<u64>,
    pub quality_by_pos: Vec<(u64, u64)>,
    pub qual_hist_by_pos: Vec<[u32; 43]>,
    pub base_composition: Vec<[u64; 5]>,
    pub quality_distribution: Vec<u64>,
    pub kmer_table: [u64; 256],
    pub fingerprints: Vec<u64>,
    pub quality_by_len_bin: [(u64, u64); 21],
    pub gc_per_read: [u64; 101],
    pub overrep_seqs: Vec<[u8; 50]>,
    pub active_pos: usize,
    pub qual_vs_length: Vec<(u32, f32)>,
    pub per_tile: HashMap<u32, (u64, u64)>,
}

impl Default for BatchAccum {
    fn default() -> Self {
        Self {
            read_count: 0,
            total_bases: 0,
            gc_count: 0,
            quality_sum: 0,
            quality_bases: 0,
            q20_bases: 0,
            q30_bases: 0,
            adapter_hits: 0,
            min_length: u64::MAX,
            max_length: 0,
            length_hist: vec![0u64; 2001],
            quality_by_pos: vec![(0, 0); MAX_QUAL_POSITION],
            qual_hist_by_pos: vec![[0u32; 43]; MAX_QUAL_POSITION],
            base_composition: vec![[0u64; 5]; MAX_QUAL_POSITION],
            quality_distribution: vec![0u64; PHRED_BUCKETS],
            kmer_table: [0u64; 256],
            fingerprints: Vec::new(),
            quality_by_len_bin: [(0u64, 0u64); 21],
            gc_per_read: [0u64; 101],
            overrep_seqs: Vec::new(),
            active_pos: 0,
            qual_vs_length: Vec::new(),
            per_tile: HashMap::new(),
        }
    }
}

impl BatchAccum {
    #[inline]
    pub fn process(&mut self, seq: &[u8], qual: Option<&[u8]>, fp: u64, adapter_hit: bool) {
        self.read_count += 1;
        self.fingerprints.push(fp);

        let len = seq.len() as u64;
        self.total_bases += len;
        if len > 0 {
            if len < self.min_length { self.min_length = len; }
            if len > self.max_length { self.max_length = len; }
            self.length_hist[len.min(2000) as usize] += 1;
        }

        if adapter_hit { self.adapter_hits += 1; }

        let end = seq.len().min(MAX_QUAL_POSITION);
        if end > self.active_pos { self.active_pos = end; }
        for pos in 0..end {
            let bi = BASE_LUT[seq[pos] as usize] as usize;
            self.base_composition[pos][bi] += 1;
        }

        let gc_this = count_gc(seq);
        self.gc_count += gc_this;
        if len > 0 {
            let gc_pct = (gc_this * 100 / len).min(100) as usize;
            self.gc_per_read[gc_pct] += 1;
        }

        if seq.len() >= 4 && self.read_count <= OVERREP_SAMPLE as u64 {
            let n = seq.len().min(50);
            let mut key = [0u8; 50];
            key[..n].copy_from_slice(&seq[..n]);
            self.overrep_seqs.push(key);
            let klen = seq.len();
            for pos in 0..klen.saturating_sub(3) {
                if let Some(idx) = kmer4_index(seq, pos) {
                    self.kmer_table[idx as usize] += 1;
                }
            }
        }

        if let Some(q) = qual {
            let n = q.len();

            let phred_sum: u64 = q.iter().map(|&b| b.saturating_sub(33) as u64).sum();
            let q20: u64 = q.iter().filter(|&&b| b >= 53).count() as u64;
            let q30: u64 = q.iter().filter(|&&b| b >= 63).count() as u64;

            self.quality_sum += phred_sum;
            self.quality_bases += n as u64;
            self.q20_bases += q20;
            self.q30_bases += q30;

            let cap = n.min(MAX_QUAL_POSITION);
            if cap > self.active_pos { self.active_pos = cap; }

            for (pos, &b) in q[..cap].iter().enumerate() {
                let phred = b.saturating_sub(33) as u64;
                self.quality_by_pos[pos].0 += phred;
                self.quality_by_pos[pos].1 += 1;
            }
            for (pos, &b) in q[..cap].iter().enumerate() {
                let p = b.saturating_sub(33).min(42) as usize;
                self.qual_hist_by_pos[pos][p] += 1;
            }

            if n > 0 {
                let bucket = (phred_sum as f64 / n as f64) as usize;
                self.quality_distribution[bucket.min(42)] += 1;

                let bin = (n / 100).min(20);
                self.quality_by_len_bin[bin].0 += phred_sum;
                self.quality_by_len_bin[bin].1 += n as u64;

                // Sample qual-vs-length scatter (up to 2000 points)
                if self.qual_vs_length.len() < 2000 {
                    let mean_phred = phred_sum as f32 / n as f32;
                    self.qual_vs_length.push((n as u32, mean_phred));
                }
            }
        }
    }
}

/// Hash the first 50 bytes of a sequence for duplication fingerprinting.
pub fn fingerprint(seq: &[u8]) -> u64 {
    xxh3_64(&seq[..seq.len().min(50)])
}

pub const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Convert kmer count table to HashMap.
pub fn kmer_table_to_map(table: &[u64; 256]) -> std::collections::HashMap<[u8; 4], u64> {
    let mut map = std::collections::HashMap::with_capacity(256);
    for idx in 0u8..=255 {
        let count = table[idx as usize];
        if count == 0 { continue; }
        let a = BASES[((idx >> 6) & 3) as usize];
        let b = BASES[((idx >> 4) & 3) as usize];
        let c = BASES[((idx >> 2) & 3) as usize];
        let d = BASES[(idx & 3) as usize];
        map.insert([a, b, c, d], count);
    }
    map
}

/// Estimate duplication rate from a Vec of fingerprints.
/// Returns (dup_rate_pct, dup_level_histogram_9bins).
pub fn estimate_dup_rate(fps: &[u64]) -> (f64, Vec<f64>) {
    if fps.is_empty() {
        return (0.0, vec![0.0; 9]);
    }
    let mut sorted = fps.to_vec();
    sorted.sort_unstable();

    // Count runs (each run = one unique sequence, run_len = duplication level)
    let mut dup_counts = [0u64; 9]; // bins: 1,2,3,4,5,6,7,8,9+
    let mut unique = 0u64;
    let mut i = 0;
    while i < sorted.len() {
        let mut run_len = 1usize;
        while i + run_len < sorted.len() && sorted[i + run_len] == sorted[i] {
            run_len += 1;
        }
        unique += 1;
        let bin = (run_len - 1).min(8);
        dup_counts[bin] += 1;
        i += run_len;
    }

    let total = fps.len() as f64;
    let dup_rate = ((total - unique as f64) / total * 100.0).max(0.0);

    let hist: Vec<f64> = dup_counts.iter().map(|&c| c as f64 / unique as f64 * 100.0).collect();

    (dup_rate, hist)
}
