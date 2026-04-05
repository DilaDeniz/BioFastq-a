use crate::types::{ADAPTER_MATCH_LEN, ADAPTERS};

/// Find the leftmost position where a known adapter begins.
#[allow(dead_code)]
pub fn find_adapter_pos(seq: &[u8]) -> Option<usize> {
    find_adapter_pos_with_custom(seq, &[])
}

/// Find leftmost adapter position, including any custom adapter sequences.
pub fn find_adapter_pos_with_custom(seq: &[u8], custom: &[Vec<u8>]) -> Option<usize> {
    let builtin = ADAPTERS.iter().filter_map(|(_, adapter)| {
        let n = adapter.len().min(ADAPTER_MATCH_LEN);
        seq.windows(n).position(|w| w == &adapter[..n])
    });
    let custom_hits = custom.iter().filter_map(|adapter| {
        let n = adapter.len().min(ADAPTER_MATCH_LEN);
        if n == 0 { return None; }
        seq.windows(n).position(|w| w == &adapter[..n])
    });
    builtin.chain(custom_hits).min()
}

/// Trim trailing bases whose Phred score is below `threshold`.
/// Returns the index to cut at (keep `seq[..cut]`).
pub fn quality_trim_3p(qual: &[u8], threshold: u8) -> usize {
    let mut cut = qual.len();
    while cut > 0 && qual[cut - 1].saturating_sub(33) < threshold {
        cut -= 1;
    }
    cut
}
