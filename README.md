<div align="center">

<img src="assets/logo.png" width="160" alt="BioFastq-A logo"/>

# BioFastq-A

**High-performance FASTQ/FASTA quality analysis — written in Rust**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.80%2B-orange.svg)](https://www.rust-lang.org)
[![Version](https://img.shields.io/badge/version-2.2.0-blue.svg)](Cargo.toml)

No Java. No Python. No internet required.

</div>

---

## Benchmark

Real Illumina data — **SRR38033288** (43.5M reads · 6.08 Gbp · ~14 GB on disk · 4 threads · cold cache · WSL2)

| Tool | Time | What it does |
|---|---|---|
| **BioFastq-A** | **33s** | overrep seqs · k-mers · dup · per-tile · N50/N90 · HTML report |
| fastp | 58s | adapter trimming · QC · k-mers |
| FastQC | 168s | similar analysis depth to BioFastq-A |

---

## Quick start

```bash
# Build (one-time, ~20s)
cargo build --release

# Interactive TUI — live dashboard while processing
./target/release/biofastq-a sample.fastq

# Headless — for scripts and CI
./target/release/biofastq-a sample.fastq --headless --output-dir ./reports

# Trim adapters + analyse
./target/release/biofastq-a reads.fastq --trim --output-dir ./qc

# Open the report
xdg-open ./reports/sample_report.html   # Linux
open ./reports/sample_report.html       # macOS
explorer.exe ./reports/sample_report.html  # WSL
```

---

## Installation

<details>
<summary><b>Build from source (recommended)</b></summary>

Requires [Rust ≥ 1.80](https://rustup.rs).

```bash
git clone https://github.com/DilaDeniz/BioFastq-a.git
cd BioFastq-a
cargo build --release
# binary → target/release/biofastq-a
```

Enable native CPU optimisations (AVX2/SSE4 — recommended):

```bash
mkdir -p .cargo
echo '[build]' > .cargo/config.toml
echo 'rustflags = ["-C", "target-cpu=native"]' >> .cargo/config.toml
cargo build --release
```

Install system-wide:

```bash
bash install.sh          # → /usr/local/bin (may need sudo)
bash install.sh ~/bin    # → ~/bin (no sudo)
```

</details>

<details>
<summary><b>Docker</b></summary>

```bash
docker build -t biofastq-a .

# Run (mount current directory as /data)
docker run --rm -v "$PWD":/data biofastq-a sample.fastq --headless
docker run --rm -v "$PWD":/data biofastq-a *.fastq.gz --trim --output-dir /data/qc
```

</details>

<details>
<summary><b>Homebrew (macOS / Linux)</b></summary>

```bash
brew tap DilaDeniz/biofastq-a
brew install biofastq-a
```

</details>

---

## Usage

```
biofastq-a [OPTIONS] <file> [<file2> ...]

OPTIONS:
  --headless             No TUI — for scripts and CI
  --output-dir <dir>     Where to write reports (default: current directory)
  --trim                 Trim adapters; write <stem>_trimmed.fastq.gz
  --min-length <N>       Drop trimmed reads shorter than N bp (default: 20)
  --adapter <seq>        Additional adapter sequence to screen/trim (repeatable)
  --quality-trim <Q>     Trim 3' bases with Phred quality below Q (default: off)
  --threads <N>          Number of CPU threads (default: all cores)
  --strict               Abort on first malformed record (default: skip & warn)
  --paired-end <R2>      Paired-end mode: provide R2 file path
  --version, -V          Print version
  --help, -h             Show help
```

---

## Features

<details>
<summary><b>QC Modules</b></summary>

| Module | Details |
|---|---|
| Per-base quality | Phred per position up to 500 bp · Q20/Q28/Q30 zone shading |
| Per-sequence quality | Read-level mean Phred distribution |
| Base composition | A/C/G/T/N % per position |
| GC content | Overall + FastQC-style pass/warn/fail |
| N content | N % per position |
| Sequence length | Distribution chart · N50 · N90 |
| Duplication | Fingerprint-hashes first 200k reads · deterministic |
| Overrepresented seqs | Top sequences by frequency · adapter source detection |
| Adapter content | 7 built-in sequences + custom via `--adapter` |
| Per-tile quality | Illumina CASAVA 1.8+ tile IDs · bar chart per tile |
| K-mer analysis | Parallel 4-mer counting · top enriched k-mers |

Each module shows a **FastQC-style traffic light** (Pass / Warn / Fail).

</details>

<details>
<summary><b>Output</b></summary>

```
<stem>_report.html      — self-contained HTML report (offline, no CDN)
<stem>_report.json      — machine-readable JSON for pipelines
<stem>_trimmed.fastq.gz — trimmed reads (only with --trim)
```

For multiple input files: one report per file + `batch_report.html` summary.

</details>

<details>
<summary><b>Adapters detected</b></summary>

| Name | Sequence (prefix matched) |
|---|---|
| TruSeq Read 1 | `AGATCGGAAGAGCACACGTCT` |
| TruSeq Read 2 | `AGATCGGAAGAGCGTCGTGTA` |
| Nextera Read 1/2 | `CTGTCTCTTATACACATCT` |
| Small RNA 3′ | `TGGAATTCTCGGGTGCCAAGG` |
| Poly-A | `AAAAAAAAAAAAAAAAAAAAAA` |
| Poly-T | `TTTTTTTTTTTTTTTTTTTTTT` |

Add custom adapters with `--adapter SEQUENCE` (repeatable).

</details>

---

## How it's fast

- **mmap zero-copy reader** — sequence data never copied to heap
- **RecordRange descriptors** — byte offsets into shared mmap, no allocations in hot path  
- **crossbeam I/O pipeline** — reader thread and rayon workers run in parallel
- **BASE_LUT** — 256-entry lookup table replaces 5-way branch per base
- **AVX2 quality loops** — phred sum, Q20, Q30 as separate vectorised passes (32 bytes/cycle)
- **K-mer sampling** — capped at first 200k reads, not the full file

---

## Comparison

<details>
<summary><b>vs FastQC</b></summary>

| | BioFastq-A | FastQC |
|---|---|---|
| Language | Rust | Java |
| Speed (real data) | **33s** / 6.08 Gbp | 168s / 6.08 Gbp |
| Interactive TUI | **Yes** | No |
| Adapter trimming | **Yes** | No |
| N50 / N90 | **Yes** | No |
| Long-read support | **Yes** | Limited |
| Offline / no deps | **Yes** | Requires JVM |
| HTML report | Yes | Yes |
| Per-tile quality | Yes | Yes |
| Duplication estimate | Yes | Yes |

</details>

<details>
<summary><b>vs fastp</b></summary>

| | BioFastq-A | fastp |
|---|---|---|
| Language | Rust | C++ |
| Speed (real data) | **33s** / 6 Gb | 58s / 6 Gb |
| Interactive TUI | **Yes** | No |
| N50 / N90 | **Yes** | No |
| Per-tile quality | **Yes** | No |
| Overrepresented seqs | **Yes** | No |
| FastQC traffic lights | **Yes** | No |
| Multi-file batch | **Yes** | No |
| Paired-end support | Yes | **Yes (default)** |
| Auto adapter detection | No | **Yes** |
| Poly-G tail trim | No | **Yes** |

</details>

---

## Pipeline integration

<details>
<summary><b>Snakemake</b></summary>

```python
rule fastq_qc:
    input:  "data/{sample}.fastq.gz"
    output:
        html = "qc/{sample}_report.html",
        json = "qc/{sample}_report.json"
    shell:
        "biofastq-a {input} --headless --output-dir qc/"
```

</details>

<details>
<summary><b>Nextflow</b></summary>

```groovy
process BIOFASTQA {
    input:  path fastq
    output: path "*_report.{html,json}"
    script:
    """
    biofastq-a ${fastq} --headless --output-dir .
    """
}
```

</details>

---

## Support

If you find this project useful, consider sending a small tip. Due to age restrictions I'm unable to use traditional payment platforms — crypto is the only way I can receive support. Thank you!

| Network | Address |
|---|---|
| **Solana (SOL)** | `AY5SwVxbvTHL16SUGj6kJBqMk4USniZmbqdXxH8xVrTa` |
| **Ethereum (ETH)** | `0x5176d005DD096aFa145B3ffff308b72ed76f1554` |

---

## License

MIT
