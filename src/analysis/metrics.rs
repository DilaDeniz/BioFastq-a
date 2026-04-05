use std::collections::HashMap;
use std::hash::{DefaultHasher, Hash, Hasher};

use crate::types::MAX_QUAL_POSITION;

/// Accumulate per-read stats into running counters.
#[allow(clippy::too_many_arguments)]
pub fn accum_read_stats(
    seq: &[u8],
    qual: Option<&[u8]>,
    total_bases: &mut u64,
    gc_count: &mut u64,
    min_length: &mut u64,
    max_length: &mut u64,
    length_histogram: &mut HashMap<u64, u64>,
    quality_sum: &mut u64,
    quality_bases: &mut u64,
    q20_bases: &mut u64,
    q30_bases: &mut u64,
    quality_by_pos: &mut Vec<(u64, u64)>,
    base_composition: &mut Vec<[u64; 5]>,
    quality_distribution: &mut Vec<u64>,
    batch_seqs: &mut Vec<Vec<u8>>,
) {
    let len = seq.len() as u64;
    *total_bases += len;
    if len < *min_length { *min_length = len; }
    if len > *max_length { *max_length = len; }
    *length_histogram.entry(len).or_insert(0) += 1;

    for (pos, &base) in seq.iter().enumerate() {
        if pos >= MAX_QUAL_POSITION { break; }
        let idx = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4,
        };
        base_composition[pos][idx] += 1;
        if matches!(base, b'G' | b'g' | b'C' | b'c') {
            *gc_count += 1;
        }
    }

    if let Some(q) = qual {
        accum_quality(q, quality_sum, quality_bases, q20_bases, q30_bases, quality_by_pos, quality_distribution);
    }

    if seq.len() >= 4 {
        batch_seqs.push(seq.to_vec());
    }
}

pub fn accum_quality(
    qual: &[u8],
    quality_sum: &mut u64,
    quality_bases: &mut u64,
    q20_bases: &mut u64,
    q30_bases: &mut u64,
    quality_by_pos: &mut Vec<(u64, u64)>,
    quality_distribution: &mut Vec<u64>,
) {
    let mut read_q_sum = 0u64;
    let n = qual.len() as u64;

    for (pos, &q) in qual.iter().enumerate() {
        let phred = q.saturating_sub(33) as u64;
        *quality_sum += phred;
        *quality_bases += 1;
        read_q_sum += phred;
        if phred >= 20 { *q20_bases += 1; }
        if phred >= 30 { *q30_bases += 1; }
        if pos < quality_by_pos.len() {
            quality_by_pos[pos].0 += phred;
            quality_by_pos[pos].1 += 1;
        }
    }

    if n > 0 {
        let bucket = (read_q_sum as f64 / n as f64) as usize;
        quality_distribution[bucket.min(42)] += 1;
    }
}

/// Hash the first 50 bytes of a sequence for duplication fingerprinting.
pub fn fingerprint(seq: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    seq[..seq.len().min(50)].hash(&mut hasher);
    hasher.finish()
}
