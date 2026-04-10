use std::collections::HashMap;
use xxhash_rust::xxh3::xxh3_64;
use crate::types::PHRED_BUCKETS;
use super::mmap_reader::{count_gc_simd, count_n_simd};

use crate::types::MAX_QUAL_POSITION;

// ---------------------------------------------------------------------------
// BatchAccum — thread-local accumulator used by rayon fold/reduce
// ---------------------------------------------------------------------------

/// All stats for a batch of records. Designed to be folded in parallel.
#[derive(Clone)]
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
    pub length_histogram: HashMap<u64, u64>,
    pub quality_by_pos: Vec<(u64, u64)>,
    pub base_composition: Vec<[u64; 5]>,
    pub quality_distribution: Vec<u64>,
    pub per_tile: HashMap<u32, (u64, u64)>,
    pub kmer_seqs: Vec<Vec<u8>>,
    pub fingerprints: Vec<u64>,
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
            length_histogram: HashMap::new(),
            quality_by_pos: vec![(0, 0); MAX_QUAL_POSITION],
            base_composition: vec![[0u64; 5]; MAX_QUAL_POSITION],
            quality_distribution: vec![0u64; PHRED_BUCKETS],
            per_tile: HashMap::new(),
            kmer_seqs: Vec::new(),
            fingerprints: Vec::new(),
        }
    }
}

impl BatchAccum {
    /// Process a single record into this accumulator.
    pub fn process(
        &mut self,
        seq: &[u8],
        qual: Option<&[u8]>,
        tile: Option<(u32, u64, u64)>,
        fp: u64,
        adapter_hit: bool,
    ) {
        self.read_count += 1;
        self.fingerprints.push(fp);

        let len = seq.len() as u64;
        self.total_bases += len;
        if len < self.min_length { self.min_length = len; }
        if len > self.max_length { self.max_length = len; }
        *self.length_histogram.entry(len).or_insert(0) += 1;

        if adapter_hit { self.adapter_hits += 1; }

        // GC count via SIMD (whole sequence, fast)
        self.gc_count += count_gc_simd(seq);

        // Per-base composition (up to MAX_QUAL_POSITION positions)
        let end = seq.len().min(MAX_QUAL_POSITION);
        for pos in 0..end {
            let idx = match seq[pos] {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _           => 4,
            };
            self.base_composition[pos][idx] += 1;
        }
        // N content tracked separately via SIMD for accuracy
        let _ = count_n_simd(seq); // counted in base_composition[..][4] above

        // Quality
        if let Some(q) = qual {
            let mut read_q_sum = 0u64;
            let n = q.len();
            for (pos, &qb) in q.iter().enumerate() {
                let phred = qb.saturating_sub(33) as u64;
                self.quality_sum += phred;
                self.quality_bases += 1;
                read_q_sum += phred;
                if phred >= 20 { self.q20_bases += 1; }
                if phred >= 30 { self.q30_bases += 1; }
                if pos < MAX_QUAL_POSITION {
                    self.quality_by_pos[pos].0 += phred;
                    self.quality_by_pos[pos].1 += 1;
                }
            }
            if n > 0 {
                let bucket = (read_q_sum as f64 / n as f64) as usize;
                self.quality_distribution[bucket.min(42)] += 1;
            }
        }

        // Per-tile
        if let Some((tile_id, phred_sum, count)) = tile {
            let e = self.per_tile.entry(tile_id).or_insert((0, 0));
            e.0 += phred_sum;
            e.1 += count;
        }

        // K-mer input
        if seq.len() >= 4 {
            self.kmer_seqs.push(seq.to_vec());
        }
    }

    /// Merge another accumulator into self (rayon reduce).
    pub fn merge(mut self, other: Self) -> Self {
        self.read_count    += other.read_count;
        self.total_bases   += other.total_bases;
        self.gc_count      += other.gc_count;
        self.quality_sum   += other.quality_sum;
        self.quality_bases += other.quality_bases;
        self.q20_bases     += other.q20_bases;
        self.q30_bases     += other.q30_bases;
        self.adapter_hits  += other.adapter_hits;
        if other.min_length < self.min_length { self.min_length = other.min_length; }
        if other.max_length > self.max_length { self.max_length = other.max_length; }

        for (k, v) in other.length_histogram {
            *self.length_histogram.entry(k).or_insert(0) += v;
        }
        for i in 0..MAX_QUAL_POSITION {
            self.quality_by_pos[i].0    += other.quality_by_pos[i].0;
            self.quality_by_pos[i].1    += other.quality_by_pos[i].1;
            for j in 0..5 {
                self.base_composition[i][j] += other.base_composition[i][j];
            }
        }
        for i in 0..PHRED_BUCKETS {
            self.quality_distribution[i] += other.quality_distribution[i];
        }
        for (k, v) in other.per_tile {
            let e = self.per_tile.entry(k).or_insert((0, 0));
            e.0 += v.0;
            e.1 += v.1;
        }
        self.kmer_seqs.extend(other.kmer_seqs);
        self.fingerprints.extend(other.fingerprints);
        self
    }
}

// ---------------------------------------------------------------------------
// Legacy helpers (still used for sequential trim-mode processing)
// ---------------------------------------------------------------------------

/// Hash the first 50 bytes of a sequence for duplication fingerprinting.
/// Uses xxHash3 — fast, high-quality, low collision rate.
pub fn fingerprint(seq: &[u8]) -> u64 {
    xxh3_64(&seq[..seq.len().min(50)])
}
