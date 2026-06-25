use std::collections::HashMap;

pub const MAX_QUAL_POSITION: usize = 2000;
pub const PHRED_BUCKETS: usize = 43;
pub const OVERREP_SAMPLE: usize = 200_000;
pub const ADAPTER_MATCH_LEN: usize = 11;

pub const ADAPTERS: &[(&str, &[u8])] = &[
    ("TruSeq Read 1",          b"AGATCGGAAGAGCACACGTCT"),
    ("TruSeq Read 2",          b"AGATCGGAAGAGCGTCGTGTA"),
    ("Nextera Read 1",         b"CTGTCTCTTATACACATCT"),
    ("Nextera Read 2",         b"CTGTCTCTTATACACATCT"),
    ("Small RNA 3-prime",      b"TGGAATTCTCGGGTGCCAAGG"),
    ("Poly-A",                 b"AAAAAAAAAAAAAAAAAAAAAA"),
    ("Poly-T",                 b"TTTTTTTTTTTTTTTTTTTTTT"),
    ("ONT Ligation (SQK-LSK)", b"AATGTACTTCGTTCAGTTACGTATTGCT"),
    ("ONT PCR (SQK-PCS)",      b"ACTTGCCTGTCGCTCTATCTTC"),
    ("ONT Rapid (SQK-RAD)",    b"GCTTGGGTGTTTAACCTTTTTTCGCAACGGGT"),
];

#[derive(Clone, PartialEq, Debug)]
pub enum QcStatus {
    Pass,
    Warn,
    Fail,
}

impl QcStatus {
    pub fn css_class(&self) -> &'static str {
        match self {
            QcStatus::Pass => "pass",
            QcStatus::Warn => "warn",
            QcStatus::Fail => "fail",
        }
    }
    pub fn icon(&self) -> &'static str {
        match self {
            QcStatus::Pass => "✓",
            QcStatus::Warn => "!",
            QcStatus::Fail => "✗",
        }
    }
}

#[derive(Clone)]
pub struct OverrepSeq {
    pub sequence: String,
    pub count: u64,
    pub percentage: f64,
    pub possible_source: String,
}

#[derive(Clone)]
pub struct FileStats {
    pub file_path: String,
    pub file_size: u64,
    pub read_count: u64,
    pub total_bases: u64,
    pub gc_count: u64,
    pub quality_sum: u64,
    pub quality_bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
    pub kmer_counts: HashMap<[u8; 4], u64>,
    pub bytes_processed: u64,
    pub min_length: u64,
    pub max_length: u64,
    pub length_histogram: HashMap<u64, u64>,
    pub quality_by_position: Vec<(u64, u64)>,
    pub adapter_hits: u64,
    pub trimmed_reads: u64,
    pub trimmed_bases_removed: u64,
    pub trim_output_path: Option<String>,
    pub reads_filtered: u64,
    pub dup_rate_pct: f64,
    pub per_tile_quality: HashMap<u32, (u64, u64)>,
    pub base_composition: Vec<[u64; 5]>,
    pub quality_distribution: Vec<u64>,
    pub qual_hist_by_position: Vec<[u64; 43]>,
    pub quality_by_length_bin: HashMap<u32, (u64, u64)>,
    pub trimmed_by_adapter: HashMap<String, u64>,
    pub overrepresented_sequences: Vec<OverrepSeq>,
    pub module_status: Vec<(String, QcStatus)>,
    pub gc_distribution: Vec<u64>,
    pub dup_level_histogram: Vec<f64>,
    pub long_read_mode: bool,
    pub ont_channel_counts: HashMap<u32, u32>,
    pub reads_over_time: Vec<(u64, u64)>,
    pub qual_vs_length: Vec<(u32, f32)>,
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
        if self.quality_bases == 0 {
            return 0.0;
        }
        self.q20_bases as f64 / self.quality_bases as f64 * 100.0
    }

    pub fn q30_pct(&self) -> f64 {
        if self.quality_bases == 0 {
            return 0.0;
        }
        self.q30_bases as f64 / self.quality_bases as f64 * 100.0
    }

    pub fn effective_min_length(&self) -> u64 {
        if self.min_length == u64::MAX {
            0
        } else {
            self.min_length
        }
    }

    pub fn compute_n50_n90(&self) -> (u64, u64) {
        if self.length_histogram.is_empty() || self.total_bases == 0 {
            return (0, 0);
        }
        let mut sorted: Vec<(u64, u64)> = self
            .length_histogram
            .iter()
            .map(|(&l, &c)| (l, c))
            .collect();
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

    pub fn n_content_per_position(&self) -> Vec<f64> {
        self.base_composition_pct().iter().map(|pct| pct[4]).collect()
    }

    pub fn qual_percentiles_per_position(&self) -> Vec<[f64; 3]> {
        let last = self
            .qual_hist_by_position
            .iter()
            .rposition(|h| h.iter().any(|&c| c > 0))
            .map(|p| p + 1)
            .unwrap_or(0);
        self.qual_hist_by_position[..last]
            .iter()
            .map(|hist| {
                let total: u64 = hist.iter().sum();
                if total == 0 {
                    return [0.0; 3];
                }
                let q25_thresh = total / 4;
                let q50_thresh = total / 2;
                let q75_thresh = total * 3 / 4;
                let mut cumsum = 0u64;
                let mut q25 = 0usize;
                let mut q50 = 0usize;
                let mut q75 = 0usize;
                for (score, &count) in hist.iter().enumerate() {
                    cumsum += count;
                    if q25 == 0 && cumsum >= q25_thresh {
                        q25 = score;
                    }
                    if q50 == 0 && cumsum >= q50_thresh {
                        q50 = score;
                    }
                    if q75 == 0 && cumsum >= q75_thresh {
                        q75 = score;
                    }
                    if q75 > 0 {
                        break;
                    }
                }
                [q25 as f64, q50 as f64, q75 as f64]
            })
            .collect()
    }

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

pub fn format_number(n: u64) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let len = bytes.len();
    let mut result = String::with_capacity(len + len / 3);
    for (i, &b) in bytes.iter().enumerate() {
        result.push(b as char);
        let remaining = len - i - 1;
        if remaining > 0 && remaining % 3 == 0 {
            result.push(',');
        }
    }
    result
}

pub fn format_bases(n: u64) -> String {
    if n >= 1_000_000_000 {
        format!("{:.2} Gbp", n as f64 / 1_000_000_000.0)
    } else if n >= 1_000_000 {
        format!("{:.2} Mbp", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.2} Kbp", n as f64 / 1_000.0)
    } else {
        format!("{} bp", n)
    }
}
