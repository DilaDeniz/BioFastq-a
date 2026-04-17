use std::collections::HashMap;
use xxhash_rust::xxh3::xxh3_64;
use crate::types::PHRED_BUCKETS;
use super::mmap_reader::count_gc_simd;
use crate::types::MAX_QUAL_POSITION;

const BASE_LUT: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0; t[b'a' as usize] = 0;
    t[b'C' as usize] = 1; t[b'c' as usize] = 1;
    t[b'G' as usize] = 2; t[b'g' as usize] = 2;
    t[b'T' as usize] = 3; t[b't' as usize] = 3;
    t
};

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
    /// Per-position Phred histogram — enables Q25/Q50/Q75 per position (FastQC boxplots).
    pub qual_hist_by_pos: Vec<[u64; 43]>,
    pub base_composition: Vec<[u64; 5]>,
    pub quality_distribution: Vec<u64>,
    pub per_tile: HashMap<u32, (u64, u64)>,
    pub kmer_seqs: Vec<Vec<u8>>,
    pub fingerprints: Vec<u64>,
    /// Average quality bucketed by read length (bin = len / 100bp).
    pub quality_by_length_bin: HashMap<u32, (u64, u64)>,
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
            qual_hist_by_pos: vec![[0u64; 43]; MAX_QUAL_POSITION],
            base_composition: vec![[0u64; 5]; MAX_QUAL_POSITION],
            quality_distribution: vec![0u64; PHRED_BUCKETS],
            per_tile: HashMap::new(),
            kmer_seqs: Vec::new(),
            fingerprints: Vec::new(),
            quality_by_length_bin: HashMap::new(),
        }
    }
}

impl BatchAccum {
    #[inline]
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
        if len > 0 {
            if len < self.min_length { self.min_length = len; }
            if len > self.max_length { self.max_length = len; }
            *self.length_histogram.entry(len).or_insert(0) += 1;
        }

        if adapter_hit { self.adapter_hits += 1; }

        self.gc_count += count_gc_simd(seq);

        let end = seq.len().min(MAX_QUAL_POSITION);
        for pos in 0..end {
            let idx = BASE_LUT[seq[pos] as usize] as usize;
            self.base_composition[pos][idx] += 1;
        }

        if let Some(q) = qual {
            let n = q.len();

            let phred_sum: u64 = q.iter().map(|&b| b.saturating_sub(33) as u64).sum();
            let q20: u64 = q.iter().filter(|&&b| b >= 53).count() as u64;
            let q30: u64 = q.iter().filter(|&&b| b >= 63).count() as u64;

            self.quality_sum   += phred_sum;
            self.quality_bases += n as u64;
            self.q20_bases     += q20;
            self.q30_bases     += q30;

            let cap = n.min(MAX_QUAL_POSITION);
            for (pos, &b) in q.iter().enumerate().take(cap) {
                let phred = b.saturating_sub(33) as usize;
                self.quality_by_pos[pos].0 += phred as u64;
                self.quality_by_pos[pos].1 += 1;
                self.qual_hist_by_pos[pos][phred.min(42)] += 1;
            }

            if n > 0 {
                let bucket = (phred_sum as f64 / n as f64) as usize;
                self.quality_distribution[bucket.min(42)] += 1;

                // Length-binned quality: bin 0 = 0–99bp, bin 1 = 100–199bp, etc.
                let bin = (n / 100) as u32;
                let e = self.quality_by_length_bin.entry(bin).or_insert((0, 0));
                e.0 += phred_sum;
                e.1 += n as u64;
            }
        }

        if let Some((tile_id, phred_sum, count)) = tile {
            let e = self.per_tile.entry(tile_id).or_insert((0, 0));
            e.0 += phred_sum;
            e.1 += count;
        }

        if seq.len() >= 4 && self.read_count <= crate::types::OVERREP_SAMPLE as u64 {
            self.kmer_seqs.push(seq.to_vec());
        }
    }

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
        for (i, ((qp, qhp), bc)) in self.quality_by_pos.iter_mut()
            .zip(self.qual_hist_by_pos.iter_mut())
            .zip(self.base_composition.iter_mut())
            .enumerate()
        {
            qp.0 += other.quality_by_pos[i].0;
            qp.1 += other.quality_by_pos[i].1;
            for (dst, &src) in qhp.iter_mut().zip(other.qual_hist_by_pos[i].iter()) {
                *dst += src;
            }
            for (dst, &src) in bc.iter_mut().zip(other.base_composition[i].iter()) {
                *dst += src;
            }
        }
        for (dst, &src) in self.quality_distribution.iter_mut().zip(other.quality_distribution.iter()) {
            *dst += src;
        }
        for (k, v) in other.per_tile {
            let e = self.per_tile.entry(k).or_insert((0, 0));
            e.0 += v.0;
            e.1 += v.1;
        }
        for (k, v) in other.quality_by_length_bin {
            let e = self.quality_by_length_bin.entry(k).or_insert((0, 0));
            e.0 += v.0;
            e.1 += v.1;
        }
        self.kmer_seqs.extend(other.kmer_seqs);
        self.fingerprints.extend(other.fingerprints);
        self
    }
}

#[allow(clippy::too_many_arguments)]
pub fn merge_batch_into_totals(
    acc: &BatchAccum,
    read_count:    &mut u64,
    total_bases:   &mut u64,
    gc_count:      &mut u64,
    quality_sum:   &mut u64,
    quality_bases: &mut u64,
    q20_bases:     &mut u64,
    q30_bases:     &mut u64,
    adapter_hits:  &mut u64,
    min_length:    &mut u64,
    max_length:    &mut u64,
    length_histogram:      &mut HashMap<u64, u64>,
    quality_by_pos:        &mut [(u64, u64)],
    qual_hist_by_pos:      &mut [[u64; 43]],
    base_composition:      &mut [[u64; 5]],
    quality_distribution:  &mut [u64],
    per_tile:              &mut HashMap<u32, (u64, u64)>,
    quality_by_length_bin: &mut HashMap<u32, (u64, u64)>,
) {
    *read_count    += acc.read_count;
    *total_bases   += acc.total_bases;
    *gc_count      += acc.gc_count;
    *quality_sum   += acc.quality_sum;
    *quality_bases += acc.quality_bases;
    *q20_bases     += acc.q20_bases;
    *q30_bases     += acc.q30_bases;
    *adapter_hits  += acc.adapter_hits;
    if acc.min_length < *min_length { *min_length = acc.min_length; }
    if acc.max_length > *max_length { *max_length = acc.max_length; }
    for (&k, &v) in &acc.length_histogram {
        *length_histogram.entry(k).or_insert(0) += v;
    }
    for (i, ((qp, qhp), bc)) in quality_by_pos.iter_mut()
        .zip(qual_hist_by_pos.iter_mut())
        .zip(base_composition.iter_mut())
        .enumerate()
    {
        qp.0 += acc.quality_by_pos[i].0;
        qp.1 += acc.quality_by_pos[i].1;
        for (dst, &src) in qhp.iter_mut().zip(acc.qual_hist_by_pos[i].iter()) {
            *dst += src;
        }
        for (dst, &src) in bc.iter_mut().zip(acc.base_composition[i].iter()) {
            *dst += src;
        }
    }
    for (dst, &src) in quality_distribution.iter_mut().zip(acc.quality_distribution.iter()) {
        *dst += src;
    }
    for (&k, &v) in &acc.per_tile {
        let e = per_tile.entry(k).or_insert((0, 0));
        e.0 += v.0;
        e.1 += v.1;
    }
    for (&k, &v) in &acc.quality_by_length_bin {
        let e = quality_by_length_bin.entry(k).or_insert((0, 0));
        e.0 += v.0;
        e.1 += v.1;
    }
}

/// Hash the first 50 bytes of a sequence for duplication fingerprinting.
pub fn fingerprint(seq: &[u8]) -> u64 {
    xxh3_64(&seq[..seq.len().min(50)])
}
