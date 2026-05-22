# BioFastq-A: A High-Performance FASTQ/FASTA Quality Control and Preprocessing Tool Written in Rust

**Zeynep Dila Deniz**

Correspondence: zeynepdiladeniz2010@gmail.com

---

## Abstract

Quality control of raw sequencing reads remains one of the most routinely executed steps in every next-generation sequencing (NGS) and long-read sequencing workflow, yet the tools most widely used for this task were written over a decade ago in languages that impose significant runtime overhead. Here we present BioFastq-A, a FASTQ and FASTA quality analysis and preprocessing tool implemented in Rust that is 2.4× faster than fastp and 6.8× faster than FastQC on identical hardware while computing a strictly broader set of quality metrics than either. BioFastq-A combines a memory-mapped, lock-free parallel batch pipeline with SIMD-accelerated byte searching and HyperLogLog-based duplication estimation to achieve throughputs exceeding 300 MB/s on commodity hardware. The tool supports short-read Illumina data, long-read Oxford Nanopore Technology (ONT) and PacBio data, paired-end mode, integrated adapter and quality trimming, an interactive terminal dashboard, and self-contained offline HTML reports with publication-ready chart exports. BioFastq-A is available as a single native binary on Bioconda for Linux and macOS and is distributed under the MIT license.

**Keywords:** FASTQ quality control, bioinformatics, Rust, next-generation sequencing, long-read sequencing, Oxford Nanopore, preprocessing, adapter trimming

---

## 1. Introduction

The volume of raw sequencing data produced by modern NGS platforms continues to grow at a pace that outstrips most other data-generating technologies in the life sciences. Short-read instruments such as the Illumina NovaSeq X produce terabytes of FASTQ-formatted reads per run, and long-read platforms such as the Oxford Nanopore PromethION are now capable of sustained multi-hundred-gigabase throughput. In every downstream analysis—variant calling, RNA-seq quantification, metagenomics classification, genome assembly—raw reads must first be evaluated for quality and preprocessed to remove low-quality bases and contaminating adapter sequences before meaningful biological conclusions can be drawn.

FastQC (Andrews, 2010) remains the de facto standard for sequencing quality control after more than fifteen years, but its Java-based implementation is inherently single-threaded per file and requires a Java Virtual Machine (JVM) to be installed and maintained. fastp (Chen et al., 2018) improved on FastQC substantially in both speed and functionality by adopting a multithreaded C++ implementation with integrated trimming, but it lacks several quality metrics—per-tile quality, overrepresented sequence detection, k-mer enrichment analysis, duplication level histograms—that practitioners routinely need. In practice, many pipelines chain both tools: fastp for trimming followed by FastQC for the metrics fastp omits, which doubles the I/O cost and complicates workflow management.

We developed BioFastq-A to address these gaps. Our primary design goal was to provide a single tool that performs all quality assessment and preprocessing steps in one pass over the input data, while being faster than any existing tool at each individual task. A secondary goal was first-class support for long-read data, which introduces fundamentally different quality profiles and analysis requirements compared to Illumina short reads.

---

## 2. Design and Implementation

### 2.1 Core Pipeline Architecture

BioFastq-A is implemented in Rust (version 1.80+) and compiled with link-time optimisation (LTO) and `codegen-units=1`, allowing the compiler to perform whole-program inlining and global register allocation decisions that are impossible in separate-compilation models.

For plain uncompressed FASTQ input, BioFastq-A uses memory-mapped I/O (`memmap2` crate) to map the file directly into the process address space. No intermediate `fread()` buffers are allocated; the OS kernel manages page caching and read-ahead automatically. Record boundaries are located using `memchr` (AVX2/SSE4 accelerated), and read slices are handed to worker threads as zero-copy references into the mapped region.

For compressed input (`.fastq.gz`, `.fasta.gz`) and FASTA files, the tool falls back to the `needletail` library, which provides a well-tested streaming parser with gzip decompression.

### 2.2 Lock-Free Parallel Batch Fold

All quality metric accumulation uses Rayon's parallel fold/reduce pattern. Each worker thread independently owns a `BatchAccum` struct—a collection of arrays covering all quality histograms, base content counters, GC distribution bins, and overrepresented sequence candidates. No shared mutable state exists during processing; consequently, no mutexes or atomic operations appear in the hot path.

At the end of each batch, worker accumulators are merged in a binary-tree reduction. This merge step is O(log N) in the number of threads and involves only plain array additions and `Vec::extend` (memcpy) for the overrepresented sequence buffer—deliberately avoiding HashMap merge operations, which would impose O(n × hash_time) cost per reduction level.

Overrepresented sequence candidates are stored as flat `[u8; 50]` arrays rather than heap-allocated `Vec<u8>` values. This eliminates one heap allocation and one heap deallocation per sampled read—up to 200,000 allocations per file—and ensures that the Rayon reduce step copies 50-byte stack values rather than performing pointer-chasing HashMap insertions.

### 2.3 Cache-Optimised Histograms

Per-position quality histograms are stored as `[u32; 43]` arrays (Phred 0–42) rather than `[u64; 43]`. For a maximum read length of 150 bp this halves the working set from approximately 52 KB to 26 KB, allowing the active portion of the histogram to reside in L1 cache (typically 32 KB) alongside the other hot per-position arrays. The per-position quality sum and count arrays—used for mean quality calculation—are accumulated in a separate sequential pass over each read's quality bytes, allowing LLVM to auto-vectorize the inner loop. The histogram scatter-write pass (which has a data-dependent index and cannot be vectorized) runs separately.

### 2.4 Adapter and K-mer Searching

All adapter sequences are compiled into an Aho-Corasick automaton at startup. Each read is searched in a single pass of length O(read_length), regardless of the number of adapters. This contrasts with a naive implementation that would require O(adapters × read_length) per read. Nine adapters are included by default, covering the full range of Illumina TruSeq, Nextera, and Small RNA library preparation kits as well as Oxford Nanopore SQK-LSK, SQK-PCS, and SQK-RAD adapters. Users may add arbitrary sequences via the `--adapter` flag.

Four-mer frequency analysis uses a 256-entry lookup table for base encoding (eliminating a 5-way conditional branch per base) and accumulates counts into a compact `[u32; 256]` array per thread.

### 2.5 HyperLogLog Deduplication

Duplication rate estimation uses a HyperLogLog cardinality estimator with a 1% standard error parameter rather than a HashSet. This requires approximately 2 KB of memory regardless of file size and scales to files with tens or hundreds of millions of reads without any increase in memory footprint or cache pressure. The per-read hash is computed with xxHash3, a non-cryptographic hash function that is measurably faster than MurmurHash or CityHash on modern CPUs with vector unit support.

---

## 3. Features

### 3.1 Short-Read Quality Assessment

For Illumina short-read data, BioFastq-A computes the following quality metrics in a single pass:

- **Per-position quality statistics**: mean Phred score, Q25, median, and Q75 per position, displayed as a boxplot with an IQR band overlay. FastQC-style pass/warn/fail thresholds are applied.
- **Per-read quality distribution**: histogram of mean Phred scores across all reads.
- **Per-position base composition**: A, C, G, T, and N percentages per position.
- **Per-read GC content distribution**: histogram of GC% per read with a theoretical normal distribution overlay. Bimodal distributions indicate contamination; sharp spikes at extreme values indicate adapter dimer or PhiX spike-in.
- **Per-position N content**: fraction of uncalled bases per position.
- **Read length distribution**: histogram with N50 and N90 prominently displayed.
- **Duplication level histogram**: nine-bin histogram covering 1× through ≥1000×, computed via HyperLogLog.
- **Overrepresented sequences**: top sequences ranked by frequency with adapter source identification and adaptive frequency thresholds.
- **Adapter content**: per-position adapter presence for all built-in and custom adapters.
- **Per-tile quality**: Illumina CASAVA 1.8+ flow cell tile quality map, enabling visual identification of defective tile regions that produce systematically lower quality reads.
- **K-mer enrichment**: parallel 4-mer counting with enrichment score (observed/expected ratio) for the top enriched k-mers.

### 3.2 Long-Read Support

When processing Oxford Nanopore or PacBio data—either detected automatically from header format and median read length or forced with the `--long-read` flag—BioFastq-A adjusts quality thresholds to ONT-appropriate ranges (Q8 fail, Q10 warn versus Illumina's Q20/Q28), replaces per-position quality with a quality-versus-length scatter plot (sampled up to 2,000 points), and adds two ONT-specific charts: reads over time (cumulative throughput in 5-minute buckets, parsed from `start_time=` read header fields) and channel occupancy (read distribution across flow cell pore channels, parsed from `ch=` fields). These charts are essential for identifying run failures, pore clogging, and sequencing kinetics issues that are invisible in Illumina data.

### 3.3 Integrated Trimming

BioFastq-A performs adapter trimming, sliding window quality trimming, poly-G/X tail trimming, hard 5′ and 3′ trimming, and per-read quality and N-content filtering in a single integrated pass. All trimming steps are applied in a biologically motivated order: hard trim → sliding window → poly-G/X → quality trim 3′ → adapter → length filter. Trimmed reads are written to `<stem>_trimmed.fastq.gz`.

Poly-G tail trimming is particularly important for two-colour chemistry instruments (NovaSeq, NextSeq) where empty cycles produce G basecalls rather than N. The tool detects and removes these artifactual tails automatically when `--poly-g` is specified, without requiring any reference adapter sequence.

Sliding window quality trimming (equivalent to Trimmomatic's `SLIDINGWINDOW` mode) uses an O(n) sliding sum rather than recomputing the window mean from scratch at each position.

### 3.4 Automatic Parameter Suggestions

After each analysis, BioFastq-A inspects its own results and emits concrete, actionable preprocessing recommendations—both in the terminal output and as a dedicated card in the HTML report. For example, the tool recommends `--poly-g` when it detects elevated G content at the 3′ end consistent with two-colour chemistry signal loss, `--cut-right --window-quality 20` when the mean per-position quality drops below Q20 in the final positions, and `--min-quality 20` when the mean read quality falls below typical Illumina thresholds. This feature removes the iterative trial-and-error cycle of running FastQC, reading the output, deciding on trimming parameters, and re-running that is standard practice in most labs.

### 3.5 Two-File Comparison

The `--compare` flag generates a side-by-side HTML comparison report when two input files are provided. Diverging metrics are highlighted automatically. This is useful for before/after trimming validation, sample-versus-control comparisons, and replicate reproducibility assessment.

### 3.6 Output Formats

BioFastq-A generates a self-contained HTML report (no CDN dependencies, fully functional offline), a machine-readable JSON report for downstream pipeline integration, and an optional MultiQC-compatible JSON summary file (`--multiqc`). The HTML report includes interactive canvas-based charts with per-chart PNG export buttons for producing publication-ready figures without any post-processing.

---

## 4. Performance Evaluation

### 4.1 Benchmarking Setup

We benchmarked BioFastq-A 2.3.0 against fastp 0.23.4 and FastQC 0.12.1 on two datasets. The synthetic benchmark consists of 1,000,000 reads of 150 bp generated from uniform random ACGT sequences, totalling 302 MB of plain FASTQ. The real-data benchmark uses SRR38033288, a public Illumina dataset of 43.5 million reads (6.08 Gbp, ~14 GB) retrieved from the NCBI Sequence Read Archive. All benchmarks were run on 4 CPU threads with a cold page cache to reflect realistic production conditions.

### 4.2 Results

On the synthetic benchmark, BioFastq-A completed analysis in 2.1 seconds (469,000 reads/s, 141 MB/s), versus 5.0 seconds for fastp (200,000 reads/s) and 14.3 seconds for FastQC (70,000 reads/s). BioFastq-A is thus **2.4× faster than fastp** and **6.8× faster than FastQC** on equivalent hardware.

On the real Illumina dataset (SRR38033288), BioFastq-A completed in 33 seconds versus 58 seconds for fastp and 168 seconds for FastQC—a **1.76× speedup over fastp** and **5.1× speedup over FastQC** on a real-world 14 GB file.

Critically, BioFastq-A computes a strictly broader set of statistics than either competitor: per-tile quality, overrepresented sequences, k-mer enrichment, and a full duplication level histogram are all absent from fastp; adapter trimming and per-read quality filtering are absent from FastQC.

Version 2.3.0 introduced three additional statistics relative to version 2.2.0 (per-position quality boxplots, per-read GC histograms, quality-vs-length scatter for long reads) while simultaneously improving throughput by approximately 4% through the three architectural optimisations described in Section 2.2 and 2.3. This demonstrates that the overhead of additional analyses can be eliminated by attention to data layout and parallelism rather than by feature reduction.

---

## 5. Conclusion

BioFastq-A delivers quality control and preprocessing capabilities that exceed those of FastQC and fastp individually, in a single-pass tool that is faster than either. Its Rust implementation eliminates the JVM dependency of FastQC and consistently outperforms the multithreaded C++ implementation of fastp in wall-clock time. The automatic parameter suggestion feature shortens the iterative QC-trim-QC cycle that adds meaningful time to large studies. Native long-read support makes BioFastq-A suitable for ONT and PacBio workflows without requiring a separate tool. We believe BioFastq-A is a practical replacement for the FastQC+fastp combination that is currently standard in most short-read pipelines, and a useful complement to existing long-read QC tools.

---

## Availability

BioFastq-A is available under the MIT license. Source code: https://github.com/DilaDeniz/BioFastq-a. Installation via Bioconda:

```bash
conda install -c bioconda biofastq-a
```

Precompiled binaries for Linux x86_64, Linux aarch64, macOS x86_64, and macOS arm64 are provided. The tool requires no runtime dependencies beyond the binary itself.

---

## References

Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Chen S., Zhou Y., Chen Y., Gu J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890.

Flajolet P., Fusy É., Gandouet O., Meunier F. (2007). HyperLogLog: the analysis of a near-optimal cardinality estimation algorithm. *DMTCS Proceedings*, AH, 127–146.

Aho A.V., Corasick M.J. (1975). Efficient string matching: an aid to bibliographic search. *Communications of the ACM*, 18(6), 333–340.

Bolger A.M., Lohse M., Usadel B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*, 30(15), 2114–2120.

Ewels P., Magnusson M., Lundin S., Käller M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048.
