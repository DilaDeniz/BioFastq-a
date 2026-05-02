<div align="center">

<img src="logo.png" width="180" alt="BioFastq-A logo"/>

# BioFastq-A

**The fastest FASTQ/FASTA quality analysis tool — written in Rust**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.80%2B-orange.svg)](https://www.rust-lang.org)
[![Version](https://img.shields.io/badge/version-2.3.0-blue.svg)](Cargo.toml)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/biofastq-a.svg?label=Bioconda)](https://anaconda.org/bioconda/biofastq-a)

No Java. No Python. No internet required. Just one binary.

**Faster than fastp. Faster than FastQC. Does more than both.**

</div>

---

## Benchmark

Synthetic benchmark — **1M reads × 150bp, 302 MB FASTQ, 4 threads, cold cache**

| Tool | Time | Throughput | Language |
|------|------|-----------|----------|
| **BioFastq-A** | **2.1s** | **469K reads/s · 141 MB/s** | Rust |
| fastp 0.23.4 | 5.0s | 200K reads/s · 60 MB/s | C++ |
| FastQC 0.12.1 | 14.3s | 70K reads/s · 21 MB/s | Java |

BioFastq-A is **2.4× faster than fastp** and **6.8× faster than FastQC** — while running *more* analyses than either.

Real Illumina data — **SRR38033288** (43.5M reads · 6.08 Gbp · ~14 GB · 4 threads · cold cache · WSL2)

| Tool | Time |
|------|------|
| **BioFastq-A** | **33s** |
| fastp | 58s |
| FastQC | 168s |

---

## Installation

### Bioconda (recommended)

```bash
conda install -c bioconda biofastq-a
```

### Build from source

Requires [Rust ≥ 1.80](https://rustup.rs).

```bash
git clone https://github.com/DilaDeniz/BioFastq-a.git
cd BioFastq-a
cargo build --release
# binary → target/release/biofastq-a
```

Enable native CPU optimisations (AVX2/SSE4 — recommended for maximum speed):

```bash
mkdir -p .cargo
echo '[build]' > .cargo/config.toml
echo 'rustflags = ["-C", "target-cpu=native"]' >> .cargo/config.toml
cargo build --release
```

Install system-wide:

```bash
bash install.sh          # → /usr/local/bin (may need sudo)
bash install.sh ~/bin    # → ~/bin (no sudo needed)
```

### Docker

```bash
docker build -t biofastq-a .
docker run --rm -v "$PWD":/data biofastq-a sample.fastq --headless
docker run --rm -v "$PWD":/data biofastq-a *.fastq.gz --trim --output-dir /data/qc
```

---

## Quick Start

```bash
# Interactive TUI — live dashboard while processing
biofastq-a sample.fastq

# Headless mode — for scripts and CI
biofastq-a sample.fastq --headless --output-dir ./reports

# Trim adapters, quality-filter, and analyse
biofastq-a reads.fastq.gz --trim --poly-g --cut-right --output-dir ./qc

# Paired-end mode
biofastq-a R1.fastq.gz --in2 R2.fastq.gz --headless --output-dir ./qc

# Long-read mode (Oxford Nanopore / PacBio) — auto-detected or forced
biofastq-a ont_reads.fastq --long-read --headless

# Compare two samples side by side
biofastq-a sample1.fastq sample2.fastq --compare --headless

# Fast mode — skip heavy QC modules
biofastq-a *.fastq.gz --fast --headless --output-dir ./reports

# MultiQC-compatible output for pipeline aggregation
biofastq-a sample.fastq.gz --headless --multiqc --output-dir ./qc

# Open the report
xdg-open ./qc/sample_report.html   # Linux
open ./qc/sample_report.html       # macOS
```

---

## All CLI Flags

### Core Options

| Flag | Description |
|------|-------------|
| `<file> [<file2> ...]` | Input FASTQ or FASTA files. Gzip (`.gz`) handled automatically. |
| `--headless` | Run without interactive TUI. Use for scripts, CI, HPC. |
| `--output-dir <dir>` | Directory for reports and trimmed output. Default: current directory. |
| `--threads <N>` | Number of CPU threads. Default: all available cores. |
| `--in2 <file>` | R2 file for paired-end mode (R1 is the positional argument). |
| `--strict` | Abort on first malformed record. Default: skip and warn. |
| `--version, -V` | Print version and exit. |
| `--help, -h` | Show help message. |

---

### Trimming Options

All trimming steps are applied in the correct biological order:
**hard trim → sliding window → poly-G → quality trim 3′ → adapter → poly-X → length filter**

| Flag | Description |
|------|-------------|
| `--trim` | Trim adapter sequences and write cleaned reads to `<stem>_trimmed.fastq.gz`. |
| `--adapter <seq>` | Additional adapter sequence to screen and trim. Repeatable. Combined with 7 built-in adapters. |
| `--quality-trim <Q>` | Trim 3′ bases with Phred quality below Q (classic per-base trim). Default: off. |
| `--poly-g [N]` | Trim poly-G tails of ≥ N consecutive G's. Default N=10. **Essential for NovaSeq/NextSeq** 2-color chemistry data where no-signal cycles produce G's. |
| `--poly-x [N]` | Trim any homopolymer tail (A/T/C/G) of ≥ N bases. Default N=10. |
| `--cut-right` | Sliding window 5′→3′: cut the read from the first window where mean quality < threshold to the end. Equivalent to Trimmomatic `SLIDINGWINDOW`. |
| `--cut-front` | Sliding window: advance the read start past low-quality 5′ bases. |
| `--cut-tail` | Sliding window: shrink the read end past low-quality 3′ bases. |
| `--window-size <N>` | Window size for all sliding window trims. Default: 4. |
| `--window-quality <Q>` | Phred quality threshold for all sliding window trims. Default: 20. |
| `--trim-front <N>` | Hard-trim exactly N bases from every read's 5′ end unconditionally. |
| `--trim-tail <N>` | Hard-trim exactly N bases from every read's 3′ end unconditionally. |
| `--min-length <N>` | Discard reads shorter than N bp after all trimming. Default: 20. |

---

### Per-Read Filtering Options

Applied after all trimming. Reads that fail are counted and reported.

| Flag | Description |
|------|-------------|
| `--min-quality <Q>` | Discard reads whose mean Phred quality (after trimming) is below Q. |
| `--max-n <N>` | Discard reads containing more than N uncalled (N) bases. |

---

### QC Module Toggles

All modules are **on by default**. Disable expensive modules on constrained systems or when the data type makes them irrelevant.

| Flag | What it skips | When to use |
|------|--------------|-------------|
| `--no-kmer` | 4-mer frequency analysis | When speed matters more than k-mer enrichment |
| `--no-duplication` | HyperLogLog duplicate estimation | Very large files where even HLL overhead matters |
| `--no-per-tile` | Per-tile Illumina quality | Non-Illumina data (ONT, PacBio, Ion Torrent) |
| `--no-overrep` | Overrepresented sequence detection | Non-Illumina or very large files |
| `--no-adapter` | Adapter content analysis | Already-trimmed data |
| `--fast` | kmer + duplication + per-tile + overrep | Quick QC pass, skips all heavy modules |

---

### Output Options

| Flag | Description |
|------|-------------|
| `--no-html` | Skip HTML report generation. |
| `--no-json` | Skip JSON report generation. |
| `--multiqc` | Write a MultiQC-compatible `<stem>_mqc.json` file. Drop it next to your `multiqc_data/` directory and run `multiqc .` |
| `--compare` | When processing two or more files, generate a side-by-side `comparison.html` report comparing the first two files. Works for both short-read and long-read data. |

---

### Long-Read Options

| Flag | Description |
|------|-------------|
| `--long-read` | Force long-read mode. Auto-detected when median read length > 1000 bp or ONT headers (`ch=`, `start_time=`) are present. |

In long-read mode:
- Quality thresholds adjusted to ONT ranges (warn < Q10, fail < Q8 vs. Illumina's Q28/Q20)
- Q10/Q20 stat cards replace Q20/Q30
- Per-tile quality chart hidden (Illumina-specific)
- Three additional charts: **Quality vs Length scatter**, **Reads over Time**, **Channel Occupancy**
- ONT adapters (SQK-LSK, SQK-PCS, SQK-RAD) added to detection

---

## Output Files

```
<stem>_report.html        Self-contained HTML report with interactive charts.
                           No internet required. Works offline.
<stem>_report.json        Machine-readable JSON for downstream pipelines.
<stem>_mqc.json           MultiQC general stats table (with --multiqc).
<stem>_trimmed.fastq.gz   Adapter-trimmed reads (with --trim).
comparison.html           Side-by-side comparison of two files (with --compare).
```

For multiple input files: one report per file + a `batch_report.html` summary.

---

## QC Modules

### Short-Read (Illumina)

| Module | Description |
|--------|-------------|
| Per-base sequence quality | Mean Phred per position with Q25/Q50/Q75 IQR band and median line. Pass/warn/fail traffic light. |
| Per-sequence quality distribution | Histogram of mean Phred scores across all reads. |
| Per-base sequence content | A/C/G/T/N % per position. |
| Per-sequence GC content | Per-read GC % histogram with theoretical normal distribution overlay. Contamination is visible as a secondary peak. |
| Per-base N content | N % per position. |
| Sequence length distribution | Read length histogram with N50 and N90 prominently displayed. |
| Sequence duplication levels | 9-bin histogram (1×, 2×, 3-4×, 5-9×, 10-49×, 50-99×, 100-499×, 500-999×, ≥1000×) via HyperLogLog — O(1) memory regardless of file size. |
| Overrepresented sequences | Top sequences by frequency with adapter source identification. Adaptive threshold. |
| Adapter content | Per-position adapter presence. 7 built-in adapters + custom via `--adapter`. |
| Per-tile quality | Illumina CASAVA 1.8+ flow cell tile quality. Spot defective tiles visually. |
| K-mer analysis | Parallel 4-mer counting. Top enriched k-mers ranked by observed/expected ratio. |

Each module shows a **FastQC-style traffic light** (✓ Pass / ! Warn / ✗ Fail).

### Long-Read (ONT / PacBio) — additional charts

| Module | Description |
|--------|-------------|
| Quality vs Length scatter | Mean Phred vs read length scatter plot. Sampled up to 2000 points. Reveals whether longer reads are lower quality. |
| Reads over time | Cumulative read output over run duration (5-minute buckets). ONT-only: parsed from `start_time=` in read headers. |
| Channel occupancy | Distribution of reads across flow cell channels. Reveals clogged or inactive pores. ONT-only: parsed from `ch=` in read headers. |

### Automatic Parameter Suggestions

After every analysis, BioFastq-A inspects its own results and prints actionable recommendations — in the terminal output and as an amber card in the HTML report:

```
Suggestions:
  → Consider --poly-g 10  (high G content at 3' end — likely NovaSeq/NextSeq 2-color data)
  → Consider --trim  (adapter sequences detected in 15.2% of reads)
  → Consider --cut-right --window-quality 20  (mean quality drops below Q20 at position 145)
  → Consider --min-quality 20  (mean read quality Q18.4 is below typical threshold)
```

No other tool does this. BioFastq-A analyzes your data and tells you exactly what to do next.

---

## Built-in Adapters

Detected automatically via Aho-Corasick multi-pattern search (single pass, all adapters simultaneously):

| Name | Sequence |
|------|----------|
| TruSeq Read 1 | `AGATCGGAAGAGCACACGTCT` |
| TruSeq Read 2 | `AGATCGGAAGAGCGTCGTGTA` |
| Nextera Read 1/2 | `CTGTCTCTTATACACATCT` |
| Small RNA 3′ | `TGGAATTCTCGGGTGCCAAGG` |
| Poly-A | `AAAAAAAAAAAAAAAAAAAAAA` |
| Poly-T | `TTTTTTTTTTTTTTTTTTTTTT` |
| ONT Ligation (SQK-LSK) | `AATGTACTTCGTTCAGTTACGTATTGCT` |
| ONT PCR (SQK-PCS) | `ACTTGCCTGTCGCTCTATCTTC` |
| ONT Rapid (SQK-RAD) | `GCTTGGGTGTTTAACCTTTTTTCGCAACGGGT` |

Add custom adapters with `--adapter SEQUENCE` (repeatable, combined with built-ins).

---

## HTML Report Features

The self-contained HTML report (no CDN, works completely offline) includes:

- **Interactive Canvas.js charts** for all QC modules
- **PNG download button** on every chart — export publication-ready figures
- **Print / PDF button** in the header — `@media print` CSS optimised for paper
- **FastQC-style traffic lights** — instant pass/warn/fail overview
- **Automatic parameter suggestion card** (amber) — actionable recommendations
- **Stat cards** with colour-coded thresholds (green/yellow/red)
- **Batch summary** when multiple files are processed

---

## Comparison with Other Tools

### vs FastQC

| | BioFastq-A | FastQC |
|--|-----------|--------|
| Language | Rust | Java (requires JVM) |
| Speed | **6.8× faster** | baseline |
| Parallelism | Full parallel batch fold | Single-threaded per file |
| Adapter trimming | **Yes** | No |
| Poly-G/X trimming | **Yes** | No |
| Sliding window QC | **Yes** | No |
| Long-read support | **Yes (ONT/PacBio)** | Limited |
| HyperLogLog dedup | **Yes (O(1) memory)** | No (HashSet — memory grows) |
| GC distribution overlay | **Yes** | Yes |
| Per-pos IQR band | **Yes** | Yes |
| N50 / N90 | **Yes** | No |
| Interactive TUI | **Yes** | No |
| Auto-parameter suggestions | **Yes** | No |
| Side-by-side comparison | **Yes** | No |
| Chart PNG download | **Yes** | No |
| MultiQC output | **Yes** | Yes |
| Offline HTML | **Yes** | Yes |

### vs fastp

| | BioFastq-A | fastp |
|--|-----------|-------|
| Language | Rust | C++ |
| Speed | **2.4× faster** | baseline |
| Per-tile quality | **Yes** | No |
| K-mer analysis | **Yes** | No |
| Overrepresented sequences | **Yes** | No |
| GC distribution chart | **Yes** | Yes |
| Duplication histogram | **Yes** | No |
| N content chart | **Yes** | No |
| Long-read support | **Yes** | Partial |
| Interactive TUI | **Yes** | No |
| N50 / N90 | **Yes** | No |
| Auto-parameter suggestions | **Yes** | No |
| Side-by-side comparison | **Yes** | No |
| HyperLogLog dedup | **Yes** | No |
| Poly-G trimming | **Yes** | Yes |
| Sliding window QC | **Yes** | Yes |
| Paired-end correction | No | **Yes** |
| UMI support | No | **Yes** |

---

## Why Is It Fast?

BioFastq-A is architecturally different from its counterparts in ways that matter at scale:

**Lock-free parallel batch fold** — Every rayon worker thread owns its own `BatchAccum` struct. No shared mutable state, no mutexes, no contention. Workers fold independently and merge at the end. This is the single biggest reason for the speed advantage over fastp.

**Memory-mapped I/O** — Files are mapped directly into address space with `memmap2`. The OS kernel manages page caching. No `fread()` buffer copies. On sequential reads the kernel prefetches pages automatically.

**SIMD acceleration** — `memchr` uses AVX2/SSE4 for newline and byte searching. GC counting uses hand-vectorised loops. This is 8-32× faster than scalar code on modern CPUs.

**HyperLogLog cardinality estimation** — Duplicate rate uses a fixed-size HyperLogLog structure (1% error margin) instead of a growing HashSet. Constant memory regardless of file size, eliminating cache pressure from millions of HashMap insertions.

**Aho-Corasick adapter search** — All adapter sequences are searched in a single pass over each read via a precompiled automaton. No O(adapters × read_length) loop.

**Zero-cost new features** — GC distribution, N content, duplication histogram, and quality percentiles are computed from data that was *already being collected*. They add zero meaningful overhead.

**LTO + codegen-units=1** — Link-time optimisation inlines all library code at compile time. The compiler sees the entire program at once and can make global optimisation decisions impossible in separate compilation.

---

## Pipeline Integration

### Snakemake

```python
rule fastq_qc:
    input:  "data/{sample}.fastq.gz"
    output:
        html = "qc/{sample}_report.html",
        json = "qc/{sample}_report.json"
    threads: 8
    shell:
        "biofastq-a {input} --headless --threads {threads} --output-dir qc/"

rule fastq_trim:
    input:  "data/{sample}.fastq.gz"
    output:
        trimmed = "trimmed/{sample}_trimmed.fastq.gz",
        html    = "qc/{sample}_report.html"
    shell:
        "biofastq-a {input} --trim --poly-g --cut-right "
        "--min-quality 20 --headless --output-dir qc/"

rule multiqc:
    input:  expand("qc/{sample}_mqc.json", sample=SAMPLES)
    output: "qc/multiqc_report.html"
    shell:  "multiqc qc/ -o qc/"
```

### Nextflow

```groovy
process BIOFASTQA {
    input:  path fastq
    output: path "*_report.{html,json}", path "*_mqc.json" optional true

    script:
    """
    biofastq-a ${fastq} --headless --multiqc --output-dir .
    """
}

process BIOFASTQA_TRIM {
    input:  path fastq
    output: path "*_trimmed.fastq.gz", path "*_report.html"

    script:
    """
    biofastq-a ${fastq} --trim --poly-g --cut-right \\
        --min-quality 20 --headless --output-dir .
    """
}
```

### MultiQC Integration

`--multiqc` writes `<stem>_mqc.json` — picked up automatically by `multiqc .`

The file populates the **General Statistics** table with:
Total Reads · % GC · Avg Length · Avg Quality · % Q30 · % Adapter · % Dups · N50

```bash
for f in data/*.fastq.gz; do
    biofastq-a "$f" --headless --multiqc --output-dir qc/
done
multiqc qc/ -o qc/
```

---

## Examples

```bash
# Basic QC
biofastq-a sample.fastq

# Full preprocessing pipeline for NovaSeq data
biofastq-a reads.fastq.gz \
    --trim \
    --poly-g \
    --cut-right --window-quality 20 \
    --min-quality 20 \
    --max-n 5 \
    --output-dir ./qc

# ONT long-read QC
biofastq-a nanopore.fastq --long-read --headless --output-dir ./qc

# Compare two conditions
biofastq-a control.fastq treatment.fastq --compare --headless --output-dir ./qc

# Batch process all samples, generate MultiQC-compatible output
biofastq-a *.fastq.gz --headless --multiqc --output-dir ./qc

# Fast QC (skip heavy modules)
biofastq-a large_file.fastq.gz --fast --no-html --headless

# Paired-end with adapter trimming
biofastq-a R1.fastq.gz --in2 R2.fastq.gz \
    --trim --poly-g --cut-right \
    --output-dir ./qc

# CI / automated pipeline
biofastq-a sample.fastq.gz \
    --headless \
    --no-html \
    --multiqc \
    --output-dir ./reports
```

---

## Scientific Accuracy

All calculations follow standard bioinformatics conventions and have been audited:

- **Phred quality**: `raw_byte − 33 = Phred score` (Sanger/Illumina 1.8+ encoding). Mean quality is the arithmetic mean of Phred scores, not the mean of error probabilities.
- **GC content**: `(G + C) / total_bases × 100` including N bases in the denominator.
- **N50 / N90**: Reads sorted descending by length; cumulative sum threshold uses `div_ceil` for precision at boundary values.
- **Duplication rate**: `1 − (HyperLogLog_estimate / total_reads) × 100`. HyperLogLog uses 1% error parameter.
- **Sliding window QC**: O(n) sliding sum — seed first window, add new base, subtract outgoing base. Scientifically equivalent to Trimmomatic's `SLIDINGWINDOW`.
- **Poly-G/X trimming**: Scans 3′ end backwards, counts consecutive homopolymer run, cuts if run ≥ min_len.
- **Per-position percentiles**: Index-based histogram of Phred 0–42; percentile by cumulative threshold.

---

## Support

If you find BioFastq-A useful in your research, consider citing it or starring the repository.

Due to age restrictions I'm unable to use traditional payment platforms — crypto is the only way I can receive support:

| Network | Address |
|---------|---------|
| **Solana (SOL)** | `AY5SwVxbvTHL16SUGj6kJBqMk4USniZmbqdXxH8xVrTa` |
| **Ethereum (ETH)** | `0x5176d005DD096aFa145B3ffff308b72ed76f1554` |

---

## License

MIT © 2026 Zeynep Dila Deniz
