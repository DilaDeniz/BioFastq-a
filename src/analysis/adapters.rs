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
