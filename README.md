# BioFastq-A

**High-performance FASTQ/FASTA quality analysis** — written in Rust.

Parallel k-mer counting, per-base quality charts, adapter trimming, duplication
estimation, per-tile Illumina QC, N50/N90 for long reads, and a self-contained
HTML report that opens in any browser. No Java. No Python. No internet required.

---

## Quick start

```bash
# 1. Build (one-time, ~20 s)
cargo build --release

# 2. Analyse a file — opens the live TUI dashboard
./target/release/biofastq-a sample.fastq

# 3. Headless (CI / scripts)
./target/release/biofastq-a sample.fastq.gz --headless --output-dir ./reports

# 4. Trim adapters + analyse
./target/release/biofastq-a reads.fastq --trim --output-dir ./qc

# 5. Open the report
open ./reports/sample_report.html       # macOS
xdg-open ./reports/sample_report.html  # Linux
```

---

## Installation

### Option 1 — Build from source (recommended)

Requires [Rust ≥ 1.80](https://rustup.rs).

```bash
git clone https://github.com/DilaDeniz/BioFastq-a.git
cd BioFastq-a
cargo build --release
# binary is at target/release/biofastq-a
```

Install system-wide:

```bash
bash install.sh          # installs to /usr/local/bin  (may need sudo)
bash install.sh ~/bin    # installs to ~/bin (no sudo needed)
```

### Option 2 — Docker

```bash
docker build -t biofastq-a .

# Run (mount current directory as /data)
docker run --rm -v "$PWD":/data biofastq-a sample.fastq --headless
docker run --rm -v "$PWD":/data biofastq-a *.fastq.gz --trim --output-dir /data/qc
```

### Option 3 — Homebrew (macOS / Linux)

```bash
brew tap DilaDeniz/biofastq-a
brew install biofastq-a
```

> The tap formula is in `Formula/biofastq-a.rb`.

---

## How to open the app

**Interactive TUI dashboard** (default):

```bash
biofastq-a sample.fastq
```

This opens a live terminal dashboard that updates in real time.
Press **Q** or **Esc** to exit. Reports are written automatically when
processing completes.

**Headless / batch mode** (no terminal UI):

```bash
biofastq-a sample.fastq --headless
```

Use this in shell scripts, Snakemake, or Nextflow pipelines.

**View the HTML report** — open the generated `.html` file in any browser:

```bash
# macOS
open sample_report.html

# Linux
xdg-open sample_report.html

# Windows (WSL)
explorer.exe sample_report.html
```

---

## Usage

```
biofastq-a [OPTIONS] <file> [<file2> ...]

OPTIONS:
  --headless          No TUI — for scripts and CI
  --output-dir <dir>  Where to write reports (default: current directory)
  --trim              Trim adapters; write <stem>_trimmed.fastq.gz
  --min-length <N>    Drop trimmed reads shorter than N bp (default: 20)
  --version, -V       Print version
  --help, -h          Show help
```

### Multiple files

```bash
biofastq-a lane1.fastq lane2.fastq lane3.fastq --headless --output-dir ./qc
```

Each file gets its own section in the HTML report. A combined summary table
is shown at the top.

---

## Features

| Feature | Details |
|---|---|
| **Per-base quality chart** | Phred quality plotted per position (up to 500 bp), with Q20/Q28/Q30 zone shading |
| **Read length distribution** | Bar chart; works for short reads (Illumina) and long reads (Nanopore/PacBio) |
| **N50 / N90** | Computed from the length histogram — no memory overhead |
| **GC content** | Per-read and overall |
| **Q20 / Q30 pass rates** | Fraction of reads with mean quality ≥ threshold |
| **Adapter detection** | 7 common Illumina/Nextera/Poly-A sequences checked per read |
| **Adapter trimming** | `--trim` hard-clips adapters; outputs gzipped FASTQ |
| **Duplication estimate** | Fingerprint-hashes first 200k reads; reports % likely duplicates |
| **Per-tile quality** | Parses Illumina CASAVA 1.8+ tile IDs from headers; bar chart per tile |
| **Top 20 k-mers** | Parallel 4-mer counting using Rayon |
| **HTML report** | Self-contained, offline-capable; three interactive Canvas charts |
| **JSON report** | Machine-readable; suitable for downstream pipeline parsing |
| **gzip support** | `.fastq.gz` / `.fasta.gz` read and written transparently |
| **Multi-file batch** | Any number of input files in one run |
| **Docker** | Multi-stage Dockerfile included |
| **TUI dashboard** | Real-time progress, quality sparkline, 12-metric panel |

---

## Output files

```
<stem>_report.html      — HTML report (open in browser)
<stem>_report.json      — JSON report (for pipelines)
<stem>_trimmed.fastq.gz — Trimmed reads (only with --trim)
```

For multiple input files: `batch_report.html` / `batch_report.json`.

---

## Performance

BioFastq-A uses Rayon parallel k-mer counting and processes all other metrics
in a single streaming pass. Typical throughput on a modern desktop:

| Dataset | Reads | Size | Time | Throughput |
|---|---|---|---|---|
| Illumina 150bp PE | 50M | 15 GB | ~45 s | ~340 MB/s |
| Nanopore R10 | 5M | 30 GB | ~90 s | ~340 MB/s |

*Measured on an 8-core machine with NVMe storage. Your results will vary.*

---

## Comparison with FastQC

| | BioFastq-A | FastQC |
|---|---|---|
| Language | Rust | Java |
| Speed | ~340 MB/s | ~40 MB/s |
| HTML report | Yes | Yes |
| Adapter trimming | Yes (built-in) | No (needs Trimmomatic) |
| Duplication estimate | Yes | Yes |
| Per-tile quality | Yes | Yes |
| N50 / N90 | Yes | No |
| Long-read support | Yes | Limited |
| Offline / no deps | Yes | Requires JVM |
| Docker image | Yes | Official image |
| gzip streaming | Yes | Yes |

---

## Adapters detected

| Name | Sequence (prefix checked) |
|---|---|
| TruSeq Read 1 | `AGATCGGAAGAGCACACGTCT` |
| TruSeq Read 2 | `AGATCGGAAGAGCGTCGTGTA` |
| Nextera Read 1/2 | `CTGTCTCTTATACACATCT` |
| Small RNA 3′ | `TGGAATTCTCGGGTGCCAAGG` |
| Poly-A | `AAAAAAAAAAAAAAAAAAAAAA` |
| Poly-T | `TTTTTTTTTTTTTTTTTTTTTT` |

---

## Pipeline integration

### Snakemake

```python
rule fastq_qc:
    input:  "data/{sample}.fastq.gz"
    output:
        html = "qc/{sample}_report.html",
        json = "qc/{sample}_report.json"
    shell:
        "biofastq-a {input} --headless --output-dir qc/"
```

### Nextflow

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

---

## License

MIT
