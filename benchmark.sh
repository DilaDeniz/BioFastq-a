#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------------------------
# BioFastq-A vs FastQC vs fastp benchmark
# Usage: bash benchmark.sh
# ---------------------------------------------------------------------------

WORKDIR="$HOME/benchmark_na12878"
R1="$WORKDIR/ERR194147_1.fastq.gz"
R2="$WORKDIR/ERR194147_2.fastq.gz"
RESULTS="$WORKDIR/results.txt"
REPO_DIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$WORKDIR"
echo "Results will be written to: $RESULTS"
> "$RESULTS"

# ---------------------------------------------------------------------------
# 1. Download data
# ---------------------------------------------------------------------------
echo "=== Downloading NA12878 (ERR194147) ==="
if [ ! -f "$R1" ]; then
    wget -c -P "$WORKDIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
else
    echo "R1 already exists, skipping."
fi
if [ ! -f "$R2" ]; then
    wget -c -P "$WORKDIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz
else
    echo "R2 already exists, skipping."
fi

# ---------------------------------------------------------------------------
# 2. Install tools if missing
# ---------------------------------------------------------------------------
echo "=== Checking tools ==="

if ! command -v fastqc &>/dev/null; then
    echo "FastQC not found, installing via conda..."
    conda install -y -c bioconda fastqc
else
    echo "FastQC: $(fastqc --version 2>&1 | head -1)"
fi

if ! command -v fastp &>/dev/null; then
    echo "fastp not found, installing via conda..."
    conda install -y -c bioconda fastp
else
    echo "fastp: $(fastp --version 2>&1 | head -1)"
fi

# ---------------------------------------------------------------------------
# 3. Build latest BioFastq-A
# ---------------------------------------------------------------------------
echo "=== Building latest BioFastq-A ==="
cd "$REPO_DIR"
git pull origin main
cargo build --release 2>&1 | tail -3
BIOFASTQ="$REPO_DIR/target/release/biofastq-a"

drop_cache() {
    sync
    if [ -w /proc/sys/vm/drop_caches ]; then
        echo 3 > /proc/sys/vm/drop_caches
    else
        # WSL: sudo gerekebilir, yoksa atla
        sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches' 2>/dev/null || true
    fi
    sleep 1
}

# ---------------------------------------------------------------------------
# 4. FastQC
# ---------------------------------------------------------------------------
echo ""
echo "=== Benchmarking FastQC ==="
drop_cache
START=$(date +%s%N)
fastqc --threads 4 --outdir "$WORKDIR" "$R1" "$R2" 2>&1 | tail -3
END=$(date +%s%N)
FASTQC_TIME=$(( (END - START) / 1000000 ))
echo "FastQC: ${FASTQC_TIME}ms" | tee -a "$RESULTS"

# ---------------------------------------------------------------------------
# 5. fastp
# ---------------------------------------------------------------------------
echo ""
echo "=== Benchmarking fastp ==="
drop_cache
START=$(date +%s%N)
fastp \
    --in1 "$R1" --in2 "$R2" \
    --thread 4 \
    --json "$WORKDIR/fastp_report.json" \
    --html "$WORKDIR/fastp_report.html" \
    --disable_adapter_trimming \
    2>&1 | tail -3
END=$(date +%s%N)
FASTP_TIME=$(( (END - START) / 1000000 ))
echo "fastp:  ${FASTP_TIME}ms" | tee -a "$RESULTS"

# ---------------------------------------------------------------------------
# 6. BioFastq-A
# ---------------------------------------------------------------------------
echo ""
echo "=== Benchmarking BioFastq-A ==="
drop_cache
START=$(date +%s%N)
"$BIOFASTQ" "$R1" --in2 "$R2" --headless --output-dir "$WORKDIR" 2>&1 | tail -3
END=$(date +%s%N)
BIOFASTQ_TIME=$(( (END - START) / 1000000 ))
echo "BioFastq-A: ${BIOFASTQ_TIME}ms" | tee -a "$RESULTS"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "=============================="
echo "       BENCHMARK RESULTS      "
echo "=============================="
cat "$RESULTS"
echo ""
echo "Speedup vs FastQC: $(echo "scale=1; $FASTQC_TIME / $BIOFASTQ_TIME" | bc)x"
echo "Speedup vs fastp:  $(echo "scale=1; $FASTP_TIME  / $BIOFASTQ_TIME" | bc)x"
