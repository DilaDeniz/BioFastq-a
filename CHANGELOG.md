# Changelog

All notable changes to BioFastq-A are documented here.

---

## [2.3.1] — 2026-06-24

### Summary

v2.3.1 is a correctness and stability release. It fixes a long-read accuracy regression and a paired-end pairing bug introduced by v2.3.0's performance work, fixes a packaging issue that could crash already-published Bioconda binaries on some CPUs, ships a WebAssembly build for running BioFastq-A in the browser, and removes the last hardcoded version strings from report/TUI output.

### New Features

#### WebAssembly Browser Build
BioFastq-A now compiles to a `wasm/` target with the same analysis pipeline as the native CLI — adapter detection, full metrics, HTML report generation — running entirely client-side via a drag-and-drop browser UI. No install, no upload: FASTQ files never leave the browser.

### Performance

| Change | Effect |
|---|---|
| zlib-rs backend for flate2 gzip I/O | Faster `.fastq.gz` decompression on the read path |
| `parse_illumina_tile` iterator rewrite | Removes a per-read `Vec` allocation during per-tile quality tracking |

### Bug Fixes

- **Critical (packaging):** the build was pinned to `target-cpu=native`, baking in the exact CPU instruction set of whatever machine compiled the binary. Bioconda builds once per platform and ships that same binary to everyone, so any user with a CPU that didn't match the build server's feature set got an instant crash (`SIGILL`) on startup with no explanation. The flag served no purpose — BioFastq-A's SIMD (via `memchr`) already does runtime CPU-feature detection — so it's removed with zero performance cost, confirmed via benchmark (554K reads/s on 200K reads before vs. after).
- **Critical:** N50/N90 were wrong for any read ≥2000bp (long-read/ONT/PacBio data) — v2.3.0's flat-array length histogram collapsed all such reads into one bucket keyed by 2000, corrupting the cumulative-length calculation. Now tracked exactly via an overflow map, with no performance cost for ordinary short-read data.
- **Major:** paired-end mode (`--compare`/R1+R2) could permanently shift every read pair by one position after a single malformed R2 record, silently cross-pairing the rest of the file. R2 retries are now isolated from the R1 cursor.
- FASTA input (no quality scores) showed a misleading Pass/Fail for "Per-base quality" and "Per-sequence quality" instead of omitting modules that don't apply.
- Hardcoded version strings in the HTML report and TUI title bar (`v2.0`) could drift from `Cargo.toml`; both now derive from it at compile time.

### Installation
```
# Bioconda
conda install -c bioconda biofastq-a

# Build from source
cargo install --path .
```

Have fun using!
by Dila Deniz

---

## [2.3.0] — 2026-05-18

### Summary

v2.3.0 is the biggest release since the initial Rust rewrite. Every analysis that takes extra CPU is now offset by the new parallel pipeline — experimental benchmarks show v2.3.0 is **faster than v2.2.0** on identical hardware despite computing significantly more statistics.

### New Features

#### Per-Position Quality Boxplots
The quality chart now shows the full Q25/median/Q75 band per position in addition to the mean line. This makes it immediately obvious whether a position has a tight quality distribution or a fat tail of bad reads.

#### Per-Read GC Distribution
A histogram of per-read GC% across all reads. Bimodal distributions flag contamination; a sharp spike flags adapter dimer; a broad flat distribution flags PCR over-amplification.

#### Quality-vs-Length Scatter (Long-Read Mode)
For ONT/PacBio data (`--long-read` flag), a 2D scatter of read quality vs. read length is plotted instead of per-position quality. Automatically activated when median read length > 1 000 bp.

#### Long-Read / ONT Support (`--long-read`)
- N50, N90, read-length histogram, longest/shortest read
- Correct adapter detection for ONT consensus adapters
- Quality-vs-length scatter (see above)

#### Two-File Comparison Mode (`--compare file2.fastq`)
Run BioFastq-A on two FASTQ files and produce a unified comparison report. Useful for before/after trimming QC, sample-vs-control, or replicate agreement. Side-by-side charts highlight every metric that diverges between the two files.

#### Auto Parameter Suggestions
After analysis, BioFastq-A inspects the results and prints actionable trimming recommendations:
- Adapter content detected → `--adapter-trim`
- 3′ quality drop-off → `--quality-trim 20`
- Poly-G tail (NextSeq/NovaSeq) → `--poly-g-trim`
- Overrepresented sequences → `--filter-overrep`

#### Sliding-Window Quality Trimming (`--window-size`, `--window-quality`)
Trim reads as soon as a sliding window of N bases drops below a mean quality threshold. Gentler than hard 3′ trimming; preserves more of the read when only the last few bases degrade.

#### Poly-G / Poly-X Trimming (`--poly-g-trim`, `--poly-x-trim`)
Two-colour chemistry (NextSeq, NovaSeq) encodes signal loss as poly-G. BioFastq-A now detects and removes poly-G/X tails without any adapter file.

### Performance

Three sources of overhead introduced by the new statistics were eliminated:

| Change | Effect |
|--------|--------|
| `Vec<Vec<u8>>` → `Vec<[u8;50]>` for overrep sampling | Eliminates 200K heap allocs per file; rayon reduce becomes memcpy instead of HashMap merge |
| `qual_hist_by_pos [u64;43]` → `[u32;43]` | Halves working set from ~52 KB to ~26 KB; active region fits in L1 cache |
| Split combined quality loop into two passes | Allows LLVM to auto-vectorise the sequential pass; scatter pass runs independently |

**Result on 1M × 150 bp benchmark (8 threads):**

| Version | Throughput |
|---------|-----------|
| v2.2.0 (main) | ~295 MB/s |
| v2.3.0 (this release) | ~308 MB/s |

v2.3.0 computes more statistics and is still ~4% faster.

### What BioFastq-A Does That Competitors Don't

| Feature | BioFastq-A | FastQC | fastp | MultiQC |
|---------|:---:|:---:|:---:|:---:|
| Single native binary (no JVM, no Python) | ✓ | ✗ | ✓ | ✗ |
| mmap zero-copy I/O pipeline | ✓ | ✗ | ✗ | ✗ |
| Per-position Q25/median/Q75 boxplot | ✓ | ✓ | ✗ | ✗ |
| Per-read GC histogram | ✓ | ✓ | partial | ✗ |
| Quality-vs-length scatter (long read) | ✓ | ✗ | ✗ | ✗ |
| Long-read / ONT mode (N50, N90) | ✓ | ✗ | ✗ | ✗ |
| Two-file comparison report | ✓ | ✗ | ✗ | ✓ |
| Auto trimming parameter suggestions | ✓ | ✗ | ✗ | ✗ |
| Sliding-window quality trimming | ✓ | ✗ | ✓ | ✗ |
| Poly-G / poly-X trimming | ✓ | ✗ | ✓ | ✗ |
| Built-in adapter library | ✓ | ✓ | ✓ | ✗ |
| Paired-end support | ✓ | ✓ | ✓ | ✗ |
| HyperLogLog deduplication estimate | ✓ | ✗ | ✗ | ✗ |
| JSON + HTML output | ✓ | ✗ | ✓ | ✓ |
| Bioconda package | ✓ | ✓ | ✓ | ✓ |

### Bug Fixes

- Fixed `cut_tail` edge case where reads shorter than the trim window produced an off-by-one in the output length calculation
- Fixed N90 threshold precision: was using integer division for the cumulative length threshold; now uses exact u64 arithmetic

### Installation

```bash
# Bioconda
conda install -c bioconda biofastq-a

# Build from source
cargo install --path .
```

---

## [2.2.0]

- Redesigned HTML report: scientific light theme, chart PNG downloads, Q25/Q50/Q75 overlay
- Added GC distribution, duplication level histogram, and N-content charts
- Added poly-G/X trimming, sliding window QC, global trim, and per-read filters

## [2.1.0]

- Initial Bioconda release
- Core FASTQ/FASTA quality analysis
- Overrepresented sequence detection
- k-mer analysis
- Adapter content detection
- Per-tile quality (Illumina)
- Paired-end mode
- JSON + HTML reports
