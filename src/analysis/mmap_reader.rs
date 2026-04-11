/// Zero-copy FASTQ parser backed by a memory-mapped file.
///
/// Parses the file in-place — no heap allocations for sequence data.
/// Records are represented as byte-offset descriptors (`RecordRange`) that
/// slice directly into the shared `Arc<Mmap>` during parallel processing.
use memchr::memchr;
use memmap2::Mmap;
use std::fs::File;
use std::sync::Arc;

pub struct MmapFastq {
    mmap: Arc<Mmap>,
    pos: usize,
}

/// Compact byte-offset descriptor for one FASTQ record inside the mmap.
///
/// Storing offsets (not slices) lets us collect a whole batch without
/// any `&mut mmap` conflict, then hand `Arc<Mmap>` to rayon workers
/// for fully zero-copy parallel access.
#[derive(Clone)]
pub struct RecordRange {
    pub id_start:   usize,
    pub id_len:     usize,
    pub seq_start:  usize,
    pub seq_len:    usize,
    pub qual_start: usize,
    pub qual_len:   usize,
}

impl RecordRange {
    #[inline]
    pub fn id<'a>(&self, mmap: &'a [u8]) -> &'a [u8] {
        &mmap[self.id_start .. self.id_start + self.id_len]
    }
    #[inline]
    pub fn seq<'a>(&self, mmap: &'a [u8]) -> &'a [u8] {
        &mmap[self.seq_start .. self.seq_start + self.seq_len]
    }
    #[inline]
    pub fn qual<'a>(&self, mmap: &'a [u8]) -> &'a [u8] {
        &mmap[self.qual_start .. self.qual_start + self.qual_len]
    }
}

impl MmapFastq {
    pub fn open(path: &str) -> std::io::Result<Self> {
        let file = File::open(path)?;
        // Safety: we treat the mmap as read-only and never mutate it.
        let mmap = Arc::new(unsafe { Mmap::map(&file)? });
        Ok(Self { mmap, pos: 0 })
    }

    /// Clone the Arc so callers can share the mmap with rayon workers.
    #[inline]
    pub fn mmap_arc(&self) -> Arc<Mmap> { Arc::clone(&self.mmap) }

    /// Return the next record as a `RecordRange` — no borrow of self after return.
    ///
    /// Use this to collect a whole batch into `Vec<RecordRange>` without
    /// conflicting borrows, then process with `Arc<Mmap>` in rayon.
    #[inline]
    pub fn next_range(&mut self) -> Option<RecordRange> {
        let base = self.pos;
        let mmap = &self.mmap[base..];
        if mmap.is_empty() { return None; }

        // Header line: '@' + id + '\n'
        let id_end = memchr(b'\n', mmap)?;
        // id bytes are mmap[1..id_end] (skip '@')
        let id_start = base + 1;
        let id_len   = id_end.saturating_sub(1); // exclude '@'

        // Sequence line
        let seq_start_rel = id_end + 1;
        let seq_buf       = &self.mmap[base + seq_start_rel..];
        let seq_end       = memchr(b'\n', seq_buf)?;
        let seq_len       = trim_len(seq_buf, seq_end);

        // '+' line
        let plus_start_rel = seq_start_rel + seq_end + 1;
        let plus_buf       = &self.mmap[base + plus_start_rel..];
        let plus_end       = memchr(b'\n', plus_buf)?;

        // Quality line
        let qual_start_rel = plus_start_rel + plus_end + 1;
        let qual_buf       = &self.mmap[base + qual_start_rel..];
        let qual_end       = memchr(b'\n', qual_buf).unwrap_or(qual_buf.len());
        let qual_len       = trim_len(qual_buf, qual_end);

        self.pos = base + qual_start_rel + qual_end + 1;

        Some(RecordRange {
            id_start,
            id_len,
            seq_start: base + seq_start_rel,
            seq_len,
            qual_start: base + qual_start_rel,
            qual_len,
        })
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize { self.mmap.len() }
}

/// Return the content length, trimming a trailing '\r' if present.
#[inline(always)]
fn trim_len(buf: &[u8], raw_end: usize) -> usize {
    if raw_end > 0 && buf.get(raw_end - 1) == Some(&b'\r') { raw_end - 1 } else { raw_end }
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
#[allow(dead_code)]
pub fn count_n_simd(seq: &[u8]) -> u64 {
    memchr::memchr_iter(b'N', seq).count() as u64
}
