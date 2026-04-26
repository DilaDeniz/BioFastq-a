use aho_corasick::AhoCorasick;
use std::sync::OnceLock;

use crate::types::{ADAPTER_MATCH_LEN, ADAPTERS};

// Build the Aho-Corasick automaton once at startup from built-in adapter prefixes.
static AC: OnceLock<AhoCorasick> = OnceLock::new();

fn builtin_ac() -> &'static AhoCorasick {
    AC.get_or_init(|| {
        let patterns: Vec<&[u8]> = ADAPTERS
            .iter()
            .map(|(_, a)| &a[..a.len().min(ADAPTER_MATCH_LEN)])
            .collect();
        AhoCorasick::new(patterns).expect("failed to build adapter automaton")
    })
}

/// Find the leftmost adapter hit and return (position, adapter_name).
pub fn find_adapter_with_name(seq: &[u8], custom: &[Vec<u8>]) -> Option<(usize, &'static str)> {
    let builtin_hit = builtin_ac()
        .find(seq)
        .map(|m| (m.start(), ADAPTERS[m.pattern().as_usize()].0));

    let custom_hit: Option<(usize, &'static str)> = custom.iter().filter_map(|adapter| {
        let n = adapter.len().min(ADAPTER_MATCH_LEN);
        if n == 0 { return None; }
        seq.windows(n).position(|w| w == &adapter[..n]).map(|p| (p, "custom"))
    }).min_by_key(|&(p, _)| p);

    match (builtin_hit, custom_hit) {
        (Some(a), Some(b)) => if a.0 <= b.0 { Some(a) } else { Some(b) },
        (a, b) => a.or(b),
    }
}

/// Trim trailing 3′ bases whose Phred score is below `threshold`.
/// Returns the index to cut at (keep `seq[..cut]`).
pub fn quality_trim_3p(qual: &[u8], threshold: u8) -> usize {
    let mut cut = qual.len();
    while cut > 0 && qual[cut - 1].saturating_sub(33) < threshold {
        cut -= 1;
    }
    cut
}

/// Trim a homopolymer run of `base` from the 3' end.
/// Returns cut index (keep seq[..cut]).
/// min_len: minimum run length required to trigger trimming.
pub fn trim_poly_base(seq: &[u8], base: u8, min_len: usize) -> usize {
    let n = seq.len();
    let mut run = 0usize;
    let mut i = n;
    while i > 0 && seq[i - 1] == base {
        run += 1;
        i -= 1;
    }
    if run >= min_len { n - run } else { n }
}

/// Cut-right: slide window 5'→3', cut at first window with mean Phred < threshold.
/// Returns cut index (keep seq[..cut]). Window mean = sum(phred) / window_size.
pub fn sliding_window_cut_right(qual: &[u8], window_size: usize, threshold: u8) -> usize {
    let n = qual.len();
    if window_size == 0 || n < window_size { return n; }
    // Seed the first window sum
    let mut wsum: u32 = qual[..window_size]
        .iter().map(|&b| b.saturating_sub(33) as u32).sum();
    if wsum / (window_size as u32) < (threshold as u32) { return 0; }
    for i in 1..=(n - window_size) {
        wsum -= qual[i - 1].saturating_sub(33) as u32;
        wsum += qual[i + window_size - 1].saturating_sub(33) as u32;
        if wsum / (window_size as u32) < (threshold as u32) { return i; }
    }
    n
}

/// Cut-front: scan 5'→3', advance start while window mean Phred < threshold.
/// Returns start index (keep seq[start..]).
pub fn sliding_window_cut_front(qual: &[u8], window_size: usize, threshold: u8) -> usize {
    let n = qual.len();
    if window_size == 0 || n < window_size { return 0; }
    let mut wsum: u32 = qual[..window_size]
        .iter().map(|&b| b.saturating_sub(33) as u32).sum();
    if wsum / (window_size as u32) >= (threshold as u32) { return 0; }
    for i in 1..=(n - window_size) {
        wsum -= qual[i - 1].saturating_sub(33) as u32;
        wsum += qual[i + window_size - 1].saturating_sub(33) as u32;
        if wsum / (window_size as u32) >= (threshold as u32) { return i; }
    }
    n - window_size + 1 // everything below threshold
}

/// Cut-tail: scan 3'→5', shrink end while window mean Phred < threshold.
/// Returns cut index (keep seq[..cut]).
pub fn sliding_window_cut_tail(qual: &[u8], window_size: usize, threshold: u8) -> usize {
    let n = qual.len();
    if window_size == 0 || n < window_size { return n; }
    let mut cut = n;
    let mut wsum: u32 = qual[cut - window_size..cut]
        .iter().map(|&b| b.saturating_sub(33) as u32).sum();
    if wsum / (window_size as u32) >= (threshold as u32) { return cut; }
    while cut > window_size {
        wsum -= qual[cut - 1].saturating_sub(33) as u32;
        cut -= 1;
        wsum += qual[cut - window_size].saturating_sub(33) as u32;
        if wsum / (window_size as u32) >= (threshold as u32) { return cut; }
    }
    0
}

/// Returns false if the (already trimmed) read fails per-read quality or N filters.
/// Called AFTER all trimming steps, BEFORE counting in stats.
pub fn passes_filters(seq: &[u8], qual: Option<&[u8]>, min_avg_quality: u8, max_n_bases: Option<u32>) -> bool {
    if let Some(max_n) = max_n_bases {
        let n_count = seq.iter().filter(|&&b| b == b'N' || b == b'n').count() as u32;
        if n_count > max_n { return false; }
    }
    if min_avg_quality > 0 {
        if let Some(q) = qual {
            if q.is_empty() { return false; }
            let mean = q.iter().map(|&b| b.saturating_sub(33) as u32).sum::<u32>()
                / q.len() as u32;
            if mean < min_avg_quality as u32 { return false; }
        }
    }
    true
}
