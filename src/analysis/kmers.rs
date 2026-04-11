use std::collections::HashMap;

use rayon::prelude::*;

/// Count 4-mers across a batch using rayon parallel fold/reduce.
pub fn count_kmers_parallel(sequences: &[Vec<u8>]) -> HashMap<[u8; 4], u64> {
    sequences
        .par_iter()
        .fold(HashMap::new, |mut map, seq| {
            for window in seq.windows(4) {
                if window.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                    *map.entry([window[0], window[1], window[2], window[3]])
                        .or_insert(0) += 1;
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
