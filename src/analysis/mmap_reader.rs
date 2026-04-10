/// Zero-copy FASTQ parser backed by a memory-mapped file.
///
/// Parses the file in-place — no heap allocations for sequence data.
/// Each record is a set of byte slices pointing directly into the mmap.
use memchr::memchr;
use memmap2::Mmap;
use std::fs::File;

pub struct MmapFastq {
    mmap: Mmap,
    pos: usize,
}

/// A single FASTQ record — zero-copy slices into the mmap.
pub struct RawRecord<'a> {
    pub id:   &'a [u8],
    pub seq:  &'a [u8],
    pub qual: &'a [u8],
}

impl MmapFastq {
    pub fn open(path: &str) -> std::io::Result<Self> {
        let file = File::open(path)?;
        // Safety: we treat the mmap as read-only and never mutate it.
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Self { mmap, pos: 0 })
    }

    #[inline]
    pub fn next_record(&mut self) -> Option<RawRecord<'_>> {
        let buf = &self.mmap[self.pos..];
        if buf.is_empty() { return None; }

        // Skip '@'
        let id_end   = memchr(b'\n', buf)?;
        let id       = &buf[1..id_end];  // skip leading '@'

        let seq_start = id_end + 1;
        let seq_buf   = &buf[seq_start..];
        let seq_end   = memchr(b'\n', seq_buf)?;
        let seq       = trim_cr(&seq_buf[..seq_end]);

        let plus_start = seq_start + seq_end + 1;
        let plus_buf   = &buf[plus_start..];
        let plus_end   = memchr(b'\n', plus_buf)?;

        let qual_start = plus_start + plus_end + 1;
        let qual_buf   = &buf[qual_start..];
        let qual_end   = memchr(b'\n', qual_buf).unwrap_or(qual_buf.len());
        let qual       = trim_cr(&qual_buf[..qual_end]);

        self.pos += qual_start + qual_end + 1;

        Some(RawRecord { id, seq, qual })
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize { self.mmap.len() }
}

#[inline(always)]
fn trim_cr(s: &[u8]) -> &[u8] {
    if s.last() == Some(&b'\r') { &s[..s.len()-1] } else { s }
}

// ---------------------------------------------------------------------------
// SIMD base counting — uses memchr's SIMD backend under the hood
// ---------------------------------------------------------------------------

/// Count G+C bases in a sequence slice using SIMD memchr.
#[inline]
pub fn count_gc_simd(seq: &[u8]) -> u64 {
    (memchr::memchr_iter(b'G', seq).count()
        + memchr::memchr_iter(b'C', seq).count()) as u64
}

/// Count N bases.
#[inline]
pub fn count_n_simd(seq: &[u8]) -> u64 {
    memchr::memchr_iter(b'N', seq).count() as u64
}
