use std::collections::HashMap;
use std::time::{Duration, Instant};

/// Maximum read length tracked for per-position quality and base composition.
/// 2000 covers Illumina (≤300bp), most amplicon data, and medium-length ONT reads.
/// Positions beyond this cap are still counted in global stats (GC, avg quality)
/// but are excluded from the per-position charts.
pub const MAX_QUAL_POSITION: usize = 2000;
/// How often (in reads) the UI stats panel is refreshed.  Larger = less lock
/// contention; 100k gives ~1-2 s refresh granularity at 50–100k reads/s.
pub const FLUSH_INTERVAL: u64    = 100_000;
/// Records per rayon fold batch.  50k reduces task-dispatch overhead vs 20k
/// while keeping memory bounded (~10 MB per batch at 200 B/record).
pub const PARALLEL_BATCH: usize  = 50_000;
pub const MAX_LOG_ENTRIES: usize = 1_000;
/// Phred scores range 0–42; index directly into a 43-element array.
pub const PHRED_BUCKETS: usize   = 43;
/// Number of reads sampled for overrepresented sequence and duplication analysis.
pub const OVERREP_SAMPLE: usize  = 200_000;

// Adapter sequences to check (name, sequence prefix to match)
pub const ADAPTERS: &[(&str, &[u8])] = &[
    ("TruSeq Read 1",      b"AGATCGGAAGAGCACACGTCT"),
    ("TruSeq Read 2",      b"AGATCGGAAGAGCGTCGTGTA"),
    ("Nextera Read 1",     b"CTGTCTCTTATACACATCT"),
    ("Nextera Read 2",     b"CTGTCTCTTATACACATCT"),
    ("Small RNA 3-prime",  b"TGGAATTCTCGGGTGCCAAGG"),
    ("Poly-A",             b"AAAAAAAAAAAAAAAAAAAAAA"),
    ("Poly-T",             b"TTTTTTTTTTTTTTTTTTTTTT"),
    // Oxford Nanopore Technologies adapters
    ("ONT Ligation (SQK-LSK)", b"AATGTACTTCGTTCAGTTACGTATTGCT"),
    ("ONT PCR (SQK-PCS)",      b"ACTTGCCTGTCGCTCTATCTTC"),
    ("ONT Rapid (SQK-RAD)",    b"GCTTGGGTGTTTAACCTTTTTTCGCAACGGGT"),
];

// Minimum bp of adapter prefix that must match exactly
pub const ADAPTER_MATCH_LEN: usize = 11;

#[derive(Clone, PartialEq)]
pub enum ProcessingStatus {
    Running,
    Completed,
    Error(String),
}

#[derive(Clone, PartialEq, Debug)]
pub enum LogLevel {
    Info,
    Success,
    Warning,
    Error,
}

#[derive(Clone)]
pub struct LogEntry {
    pub timestamp: Duration,
    pub level: LogLevel,
    pub message: String,
}

// ---------------------------------------------------------------------------
// QC module pass/warn/fail status — mirrors FastQC's per-module traffic lights
// ---------------------------------------------------------------------------

#[derive(Clone, PartialEq, Debug)]
pub enum QcStatus { Pass, Warn, Fail }

impl QcStatus {
    pub fn css_class(&self) -> &'static str {
        match self { QcStatus::Pass => "pass", QcStatus::Warn => "warn", QcStatus::Fail => "fail" }
    }
    pub fn icon(&self) -> &'static str {
        match self { QcStatus::Pass => "✓", QcStatus::Warn => "!", QcStatus::Fail => "✗" }
    }
    #[allow(dead_code)]
    pub fn label(&self) -> &'static str {
        match self { QcStatus::Pass => "PASS", QcStatus::Warn => "WARN", QcStatus::Fail => "FAIL" }
    }
}

/// A single overrepresented sequence found in the sample.
#[derive(Clone)]
pub struct OverrepSeq {
    pub sequence: String,
    pub count: u64,
    pub percentage: f64,
    pub possible_source: String,
}

/// Maximum number of qual-vs-length scatter points to keep.
pub const MAX_QUAL_VS_LEN_POINTS: usize = 2000;

/// Controls which optional QC modules and outputs are enabled.
/// All fields default to `true` (everything on).
#[derive(Debug, Clone)]
pub struct FeatureFlags {
    /// Run k-mer frequency analysis (--no-kmer to disable)
    pub kmer_analysis: bool,
    /// Run duplicate read estimation (--no-duplication to disable)
    pub duplication_check: bool,
    /// Parse per-tile Illumina quality (--no-per-tile to disable)
    pub per_tile_quality: bool,
    /// Detect overrepresented sequences (--no-overrep to disable)
    pub overrep_sequences: bool,
    /// Run adapter content analysis (--no-adapter to disable)
    pub adapter_detection: bool,
    /// Write HTML report (--no-html to disable)
    pub html_report: bool,
    /// Write JSON report (--no-json to disable)
    pub json_report: bool,
}

impl Default for FeatureFlags {
    fn default() -> Self {
        Self {
            kmer_analysis: true,
            duplication_check: true,
            per_tile_quality: true,
            overrep_sequences: true,
            adapter_detection: true,
            html_report: true,
            json_report: true,
        }
    }
}

/// Configuration passed to the processing engine.
#[derive(Clone)]
pub struct ProcessConfig {
    /// If set, write adapter-trimmed reads to this path (.fastq or .fastq.gz)
    pub trim_output: bool,
    /// Discard trimmed reads shorter than this (bp)
    pub min_length: u64,
    /// Directory for trimmed output files
    pub output_dir: String,
    /// Additional custom adapter sequences to screen/trim (on top of built-ins)
    pub custom_adapters: Vec<Vec<u8>>,
    /// If > 0, quality-trim the 3' end: drop trailing bases with Phred < threshold
    pub quality_trim_threshold: u8,
    /// Trim poly-G tails >= this many bp (0 = disabled). Auto-detect NovaSeq/NextSeq: not implemented yet, set manually.
    pub poly_g_min_len: u8,
    /// Trim poly-X tails >= this many bp (0 = disabled).
    pub poly_x_min_len: u8,
    /// Sliding window cut-right: window size (0 = disabled). Cuts from first bad window to read end.
    pub cut_right_window: u8,
    /// Sliding window cut-right: Phred quality threshold (default 20).
    pub cut_right_qual: u8,
    /// Sliding window cut-front: window size (0 = disabled). Trims bad quality from 5' end.
    pub cut_front_window: u8,
    /// Sliding window cut-front: Phred quality threshold (default 20).
    pub cut_front_qual: u8,
    /// Sliding window cut-tail: window size (0 = disabled). Trims bad quality from 3' end.
    pub cut_tail_window: u8,
    /// Sliding window cut-tail: Phred quality threshold (default 20).
    pub cut_tail_qual: u8,
    /// Hard-trim this many bases from the 5' end of every read (0 = disabled).
    pub trim_front_bases: u16,
    /// Hard-trim this many bases from the 3' end of every read (0 = disabled).
    pub trim_tail_bases: u16,
    /// Discard reads whose mean Phred quality (post-trim) is below this value (0 = disabled).
    pub min_avg_quality: u8,
    /// Discard reads with more than this many N bases (post-trim). None = disabled.
    pub max_n_bases: Option<u32>,
    /// Abort on first malformed record instead of skipping it
    pub strict: bool,
    /// If set, treat input as paired-end: this is the R2 file path.
    /// The primary input file is R1.
    pub paired_end_r2: Option<String>,
    /// Feature flags controlling which QC modules and outputs are active
    pub flags: FeatureFlags,
    /// Force long-read mode (ONT/PacBio) regardless of auto-detection.
    pub is_long_read: bool,
    /// Hint that the data is specifically from ONT (vs generic long-read / PacBio).
    /// Reserved for future ONT-specific processing; currently informational.
    #[allow(dead_code)]
    pub is_ont: bool,
    /// When true and trim_output is also true, run a pre-pass collecting raw stats
    /// before trimming so a before/after comparison report can be generated.
    pub compare: bool,
}

impl Default for ProcessConfig {
    fn default() -> Self {
        Self {
            trim_output: false,
            min_length: 20,
            output_dir: ".".to_string(),
            custom_adapters: Vec::new(),
            quality_trim_threshold: 0,
            poly_g_min_len: 0,
            poly_x_min_len: 0,
            cut_right_window: 0,
            cut_right_qual: 20,
            cut_front_window: 0,
            cut_front_qual: 20,
            cut_tail_window: 0,
            cut_tail_qual: 20,
            trim_front_bases: 0,
            trim_tail_bases: 0,
            min_avg_quality: 0,
            max_n_bases: None,
            strict: false,
            paired_end_r2: None,
            flags: FeatureFlags::default(),
            is_long_read: false,
            is_ont: false,
            compare: false,
        }
    }
}


/// Per-file statistics accumulated during analysis
#[derive(Clone)]
pub struct FileStats {
    pub file_path: String,
    pub file_size: u64,
    // Counts
    pub read_count: u64,
    pub total_bases: u64,
    pub gc_count: u64,
    // Quality
    pub quality_sum: u64,
    pub quality_bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
    // k-mer
    pub kmer_counts: HashMap<[u8; 4], u64>,
    // Progress
    pub bytes_processed: u64,
    // Length stats
    pub min_length: u64,
    pub max_length: u64,
    pub length_histogram: HashMap<u64, u64>,
    // Per-position quality: (sum_of_phred, count) for each position
    pub quality_by_position: Vec<(u64, u64)>,
    // Adapter contamination
    pub adapter_hits: u64,
    // Trimming stats (filled when --trim is active)
    pub trimmed_reads: u64,
    pub trimmed_bases_removed: u64,
    pub trim_output_path: Option<String>,
    /// Reads discarded by per-read quality/N filters (post-trim)
    pub reads_filtered: u64,
    // Duplication estimate via HyperLogLog (all reads)
    pub dup_rate_pct: f64,
    // Per-tile quality (Illumina only): tile_id -> (phred_sum, count)
    pub per_tile_quality: HashMap<u32, (u64, u64)>,
    // Per-base composition: [A, C, G, T, N] counts per position
    pub base_composition: Vec<[u64; 5]>,
    // Quality score distribution: count of reads per mean Phred score (index = floor(phred))
    pub quality_distribution: Vec<u64>,
    // Per-position Phred histogram for Q25/Q50/Q75 boxplot (FastQC-style)
    pub qual_hist_by_position: Vec<[u64; 43]>,
    // Average quality per length bin (bin = len/100bp): (phred_sum, count)
    pub quality_by_length_bin: HashMap<u32, (u64, u64)>,
    // Trim counts broken down by adapter name
    pub trimmed_by_adapter: HashMap<String, u64>,
    // Overrepresented sequences (sampled from first OVERREP_SAMPLE reads)
    pub overrepresented_sequences: Vec<OverrepSeq>,
    // Per-module QC pass/warn/fail status (FastQC-style traffic lights)
    pub module_status: Vec<(String, QcStatus)>,
    // Per-read GC content histogram: index = GC% (0..=100), value = read count
    pub gc_distribution: Vec<u64>,
    // Duplication level histogram: percentage of reads at each duplication level (9 bins)
    pub dup_level_histogram: Vec<f64>,
    // Long-read / ONT specific fields
    /// True when this file was processed in long-read mode.
    pub long_read_mode: bool,
    /// ONT channel → read count (only populated for ONT reads).
    pub ont_channel_counts: HashMap<u32, u32>,
    /// Cumulative reads over time: (minutes_since_run_start, cumulative_reads).
    /// Bucketed into 5-minute intervals; only non-empty for ONT reads.
    pub reads_over_time: Vec<(u64, u64)>,
    /// Sampled (read_length, mean_phred_quality) points; max 2000 entries.
    pub qual_vs_length: Vec<(u32, f32)>,
    /// When comparison mode is active (--compare --trim), this holds the RAW
    /// (pre-trim) stats; the main FileStats fields hold the post-trim stats.
    pub comparison_stats: Option<Box<FileStats>>,
}

impl FileStats {
    pub fn new(file_path: String, file_size: u64) -> Self {
        Self {
            file_path,
            file_size,
            read_count: 0,
            total_bases: 0,
            gc_count: 0,
            quality_sum: 0,
            quality_bases: 0,
            q20_bases: 0,
            q30_bases: 0,
            kmer_counts: HashMap::new(),
            bytes_processed: 0,
            min_length: u64::MAX,
            max_length: 0,
            length_histogram: HashMap::new(),
            quality_by_position: vec![(0, 0); MAX_QUAL_POSITION],
            adapter_hits: 0,
            trimmed_reads: 0,
            trimmed_bases_removed: 0,
            trim_output_path: None,
            reads_filtered: 0,
            dup_rate_pct: 0.0,
            per_tile_quality: HashMap::new(),
            base_composition: vec![[0u64; 5]; MAX_QUAL_POSITION],
            quality_distribution: vec![0u64; PHRED_BUCKETS],
            qual_hist_by_position: vec![[0u64; 43]; MAX_QUAL_POSITION],
            quality_by_length_bin: HashMap::new(),
            trimmed_by_adapter: HashMap::new(),
            overrepresented_sequences: Vec::new(),
            module_status: Vec::new(),
            gc_distribution: vec![0u64; 101],
            dup_level_histogram: vec![0.0; 9],
            long_read_mode: false,
            ont_channel_counts: HashMap::new(),
            reads_over_time: Vec::new(),
            qual_vs_length: Vec::new(),
            comparison_stats: None,
        }
    }

    pub fn gc_content(&self) -> f64 {
        if self.total_bases == 0 {
            return 0.0;
        }
        self.gc_count as f64 / self.total_bases as f64 * 100.0
    }

    pub fn avg_quality(&self) -> f64 {
        if self.quality_bases == 0 {
            return 0.0;
        }
        self.quality_sum as f64 / self.quality_bases as f64
    }

    pub fn avg_length(&self) -> f64 {
        if self.read_count == 0 {
            return 0.0;
        }
        self.total_bases as f64 / self.read_count as f64
    }

    pub fn adapter_pct(&self) -> f64 {
        if self.read_count == 0 {
            return 0.0;
        }
        self.adapter_hits as f64 / self.read_count as f64 * 100.0
    }

    pub fn q20_pct(&self) -> f64 {
        if self.quality_bases == 0 { return 0.0; }
        self.q20_bases as f64 / self.quality_bases as f64 * 100.0
    }

    pub fn q30_pct(&self) -> f64 {
        if self.quality_bases == 0 { return 0.0; }
        self.q30_bases as f64 / self.quality_bases as f64 * 100.0
    }

    /// Returns 0 if no reads have been seen yet (min_length is u64::MAX in that case)
    pub fn effective_min_length(&self) -> u64 {
        if self.min_length == u64::MAX {
            0
        } else {
            self.min_length
        }
    }

    /// Compute N50 and N90 from the length histogram.
    /// N50 = shortest read length L such that reads of length >= L cover >= 50% of total bases.
    /// N90 = same for 90%.
    pub fn compute_n50_n90(&self) -> (u64, u64) {
        if self.length_histogram.is_empty() || self.total_bases == 0 {
            return (0, 0);
        }
        let mut sorted: Vec<(u64, u64)> = self
            .length_histogram
            .iter()
            .map(|(&l, &c)| (l, c))
            .collect();
        // Sort descending by length
        sorted.sort_by(|a, b| b.0.cmp(&a.0));

        let n50_thresh = self.total_bases.div_ceil(2);
        let n90_thresh = (self.total_bases * 9).div_ceil(10);

        let mut cumsum = 0u64;
        let mut n50 = 0u64;
        let mut n90 = 0u64;

        for (len, cnt) in &sorted {
            cumsum += len * cnt;
            if n50 == 0 && cumsum >= n50_thresh {
                n50 = *len;
            }
            if n90 == 0 && cumsum >= n90_thresh {
                n90 = *len;
            }
            if n50 != 0 && n90 != 0 {
                break;
            }
        }
        (n50, n90)
    }

    /// Average Phred quality per base position, trimmed to last position with data.
    pub fn avg_qual_per_position(&self) -> Vec<f64> {
        let last = self
            .quality_by_position
            .iter()
            .rposition(|(_, c)| *c > 0)
            .map(|p| p + 1)
            .unwrap_or(0);
        self.quality_by_position[..last]
            .iter()
            .map(|(s, c)| if *c == 0 { 0.0 } else { *s as f64 / *c as f64 })
            .collect()
    }

    /// Read length distribution sorted by length ascending.
    pub fn sorted_length_dist(&self) -> Vec<(u64, u64)> {
        let mut v: Vec<_> = self
            .length_histogram
            .iter()
            .map(|(&l, &c)| (l, c))
            .collect();
        v.sort_by_key(|(l, _)| *l);
        v
    }

    pub fn trimmed_pct(&self) -> f64 {
        if self.read_count == 0 {
            return 0.0;
        }
        self.trimmed_reads as f64 / self.read_count as f64 * 100.0
    }

    /// Per-tile quality sorted by tile ID ascending.
    /// Returns vec of (tile_id, avg_phred).
    pub fn sorted_tile_quality(&self) -> Vec<(u32, f64)> {
        let mut v: Vec<(u32, f64)> = self
            .per_tile_quality
            .iter()
            .map(|(&tile, &(sum, cnt))| {
                (tile, if cnt == 0 { 0.0 } else { sum as f64 / cnt as f64 })
            })
            .collect();
        v.sort_by_key(|(tile, _)| *tile);
        v
    }

    /// Per-base composition trimmed to last position with data.
    /// Returns vec of (A%, C%, G%, T%, N%) per position.
    pub fn base_composition_pct(&self) -> Vec<[f64; 5]> {
        let last = self
            .base_composition
            .iter()
            .rposition(|counts| counts.iter().any(|&c| c > 0))
            .map(|p| p + 1)
            .unwrap_or(0);
        self.base_composition[..last]
            .iter()
            .map(|counts| {
                let total: u64 = counts.iter().sum();
                if total == 0 {
                    [0.0; 5]
                } else {
                    let t = total as f64;
                    [
                        counts[0] as f64 / t * 100.0,
                        counts[1] as f64 / t * 100.0,
                        counts[2] as f64 / t * 100.0,
                        counts[3] as f64 / t * 100.0,
                        counts[4] as f64 / t * 100.0,
                    ]
                }
            })
            .collect()
    }

    /// N content per position as percentage.
    #[allow(dead_code)]
    pub fn n_content_per_position(&self) -> Vec<f64> {
        self.base_composition_pct()
            .iter()
            .map(|pct| pct[4])
            .collect()
    }

    /// Q25, Q50 (median), Q75 per base position — FastQC-style boxplot data.
    /// Returns vec of [q25, q50, q75] in Phred units, trimmed to last covered position.
    pub fn qual_percentiles_per_position(&self) -> Vec<[f64; 3]> {
        let last = self.qual_hist_by_position
            .iter()
            .rposition(|h| h.iter().any(|&c| c > 0))
            .map(|p| p + 1)
            .unwrap_or(0);
        self.qual_hist_by_position[..last]
            .iter()
            .map(|hist| {
                let total: u64 = hist.iter().sum();
                if total == 0 { return [0.0; 3]; }
                let q25_thresh = total / 4;
                let q50_thresh = total / 2;
                let q75_thresh = total * 3 / 4;
                let mut cumsum = 0u64;
                let mut q25 = 0usize;
                let mut q50 = 0usize;
                let mut q75 = 0usize;
                for (score, &count) in hist.iter().enumerate() {
                    cumsum += count;
                    if q25 == 0 && cumsum >= q25_thresh { q25 = score; }
                    if q50 == 0 && cumsum >= q50_thresh { q50 = score; }
                    if q75 == 0 && cumsum >= q75_thresh { q75 = score; }
                    if q75 > 0 { break; }
                }
                [q25 as f64, q50 as f64, q75 as f64]
            })
            .collect()
    }

    #[allow(dead_code)]
    /// Average quality per length bin, sorted by bin index.
    pub fn avg_quality_by_length(&self) -> Vec<(u32, f64)> {
        let mut v: Vec<(u32, f64)> = self.quality_by_length_bin
            .iter()
            .map(|(&bin, &(sum, cnt))| {
                (bin * 100, if cnt == 0 { 0.0 } else { sum as f64 / cnt as f64 })
            })
            .collect();
        v.sort_by_key(|(bp, _)| *bp);
        v
    }

    /// Top N k-mers by frequency.
    pub fn top_kmers(&self, n: usize) -> Vec<(String, u64)> {
        let mut kmers: Vec<_> = self
            .kmer_counts
            .iter()
            .map(|(k, &v)| (String::from_utf8_lossy(k).to_string(), v))
            .collect();
        kmers.sort_by(|a, b| b.1.cmp(&a.1));
        kmers.truncate(n);
        kmers
    }
}

/// Global shared state accessed by both the processing thread and the UI thread.
#[derive(Clone)]
pub struct SharedState {
    /// Files that have finished processing
    pub completed_files: Vec<FileStats>,
    /// Stats for the file currently being processed (None when idle/done)
    pub current: Option<FileStats>,
    /// Index of file currently being processed
    pub current_file_idx: usize,
    /// Total number of files to process
    pub total_files: usize,
    pub start_time: Instant,
    pub status: ProcessingStatus,
    pub log_messages: Vec<LogEntry>,
}

impl SharedState {
    pub fn new(total_files: usize) -> Self {
        Self {
            completed_files: Vec::new(),
            current: None,
            current_file_idx: 0,
            total_files,
            start_time: Instant::now(),
            status: ProcessingStatus::Running,
            log_messages: Vec::new(),
        }
    }

    pub fn log(&mut self, level: LogLevel, message: String) {
        let ts = self.start_time.elapsed();
        self.log_messages.push(LogEntry {
            timestamp: ts,
            level,
            message,
        });
        if self.log_messages.len() > MAX_LOG_ENTRIES {
            // Drop oldest half to amortise repeated drains
            self.log_messages.drain(0..MAX_LOG_ENTRIES / 2);
        }
    }

    pub fn elapsed_secs(&self) -> f64 {
        self.start_time.elapsed().as_secs_f64()
    }

    /// Returns stats for the file currently being shown in the UI:
    /// the in-progress file if available, or the last completed file.
    pub fn active_stats(&self) -> Option<&FileStats> {
        self.current.as_ref().or_else(|| self.completed_files.last())
    }

    /// All files: completed + current (if any)
    pub fn all_files(&self) -> Vec<&FileStats> {
        let mut all: Vec<&FileStats> = self.completed_files.iter().collect();
        if let Some(c) = &self.current {
            all.push(c);
        }
        all
    }
}

// ---------------------------------------------------------------------------
// Shared utilities
// ---------------------------------------------------------------------------

pub fn format_number(n: u64) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let len = bytes.len();
    let mut result = String::with_capacity(len + len / 3);
    for (i, &b) in bytes.iter().enumerate() {
        result.push(b as char);
        let remaining = len - i - 1;
        if remaining > 0 && remaining.is_multiple_of(3) {
            result.push(',');
        }
    }
    result
}

pub fn format_duration(d: Duration) -> String {
    let secs = d.as_secs();
    let h = secs / 3600;
    let m = (secs % 3600) / 60;
    let s = secs % 60;
    format!("{:02}:{:02}:{:02}", h, m, s)
}

pub fn format_bases(n: u64) -> String {
    if n >= 1_000_000_000 {
        format!("{:.2} Gb", n as f64 / 1_000_000_000.0)
    } else if n >= 1_000_000 {
        format!("{:.2} Mb", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.2} Kb", n as f64 / 1_000.0)
    } else {
        format!("{} bp", n)
    }
}
