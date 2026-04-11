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

/// Find the leftmost position where any known or custom adapter begins.
pub fn find_adapter_pos_with_custom(seq: &[u8], custom: &[Vec<u8>]) -> Option<usize> {
    // Built-in adapters via Aho-Corasick (single pass over seq)
    let builtin_hit = builtin_ac()
        .find(seq)
        .map(|m| m.start());

    // Custom adapters (rare, simple scan)
    let custom_hit = custom.iter().filter_map(|adapter| {
        let n = adapter.len().min(ADAPTER_MATCH_LEN);
        if n == 0 { return None; }
        seq.windows(n).position(|w| w == &adapter[..n])
    }).min();

    match (builtin_hit, custom_hit) {
        (Some(a), Some(b)) => Some(a.min(b)),
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
