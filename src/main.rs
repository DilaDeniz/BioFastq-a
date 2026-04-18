mod analysis;
mod report;
mod tui;
mod types;

use std::io::stdout;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

use crossterm::event::{self, Event, KeyCode};
use crossterm::execute;
use crossterm::terminal::{disable_raw_mode, LeaveAlternateScreen};

use types::{format_number, FeatureFlags, ProcessConfig, ProcessingStatus, SharedState};

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

fn print_help() {
    eprintln!(
        r#"
BioFastq-A {} — High-performance FASTQ/FASTA quality analysis

USAGE:
  biofastq-a [OPTIONS] <file> [<file2> ...]

ARGUMENTS:
  <file>              FASTQ or FASTA file(s) to analyse.
                      Gzipped (.gz) files are handled automatically.

OPTIONS:
  --headless          Run without interactive TUI (for scripts / CI)
  --output-dir <dir>  Directory for reports and trimmed output
                      (default: current directory)
  --trim              Trim adapter sequences and write cleaned reads to
                      <stem>_trimmed.fastq.gz alongside the reports
  --min-length <N>    Discard trimmed reads shorter than N bp (default: 20)
  --adapter <seq>     Additional adapter sequence to screen/trim (repeatable)
  --quality-trim <Q>  Trim 3' bases with Phred quality below Q (default: off)
  --in2 <file>        R2 file for paired-end mode (R1 is the positional argument)
  --threads <N>       Number of CPU threads (default: all cores)
  --strict            Abort on first malformed record (default: skip and warn)
  --version, -V       Print version and exit
  --help, -h          Show this help message

QC MODULE TOGGLES:
  --no-kmer           Skip k-mer frequency analysis (faster, less memory)
  --no-duplication    Skip duplicate read estimation (faster, less memory)
  --no-per-tile       Skip per-tile quality analysis (non-Illumina data)
  --no-overrep        Skip overrepresented sequence detection
  --no-adapter        Skip adapter content analysis
  --fast              Quick mode: disables kmer, duplication, per-tile, overrep

OUTPUT OPTIONS:
  --no-html           Skip HTML report generation
  --no-json           Skip JSON report generation
  --multiqc           Write a MultiQC-compatible <stem>_mqc.json file
                      (drop it next to multiqc_data/ and run multiqc .)

TRIMMING OPTIONS:
  --poly-g [N]        Trim poly-G tails >= N bp (default N=10). Essential for
                      NextSeq/NovaSeq 2-color chemistry data
  --poly-x [N]        Trim any homopolymer tail >= N bp (default N=10)
  --cut-right         Sliding window 5'->3': cut from first window with mean
                      quality < threshold to read end (like Trimmomatic SLIDINGWINDOW)
  --cut-front         Sliding window: trim low-quality bases from 5' end
  --cut-tail          Sliding window: trim low-quality bases from 3' end
  --window-size N     Window size for sliding window trims (default: 4)
  --window-quality N  Phred quality threshold for sliding window trims (default: 20)
  --trim-front N      Hard-trim N bases from every read's 5' end
  --trim-tail N       Hard-trim N bases from every read's 3' end

FILTERING OPTIONS:
  --min-quality Q     Discard reads with mean Phred quality < Q (post-trim)
  --max-n N           Discard reads with more than N uncalled (N) bases

OUTPUT FILES:
  <stem>_report.html      Self-contained HTML report with interactive charts
  <stem>_report.json      Machine-readable JSON report
  <stem>_trimmed.fastq.gz Adapter-trimmed reads (when --trim is used)

EXAMPLES:
  biofastq-a sample.fastq
  biofastq-a reads.fastq.gz --trim --output-dir ./qc
  biofastq-a *.fastq.gz --headless --output-dir ./reports
  biofastq-a reads.fastq.gz --fast --no-html
  biofastq-a reads.fasta --no-per-tile --no-kmer
"#,
        env!("CARGO_PKG_VERSION")
    );
}

struct CliConfig {
    input_files: Vec<String>,
    headless: bool,
    output_dir: String,
    trim: bool,
    min_length: u64,
    custom_adapters: Vec<Vec<u8>>,
    quality_trim: u8,
    strict: bool,
    threads: Option<usize>,
    paired_r2: Option<String>,
    // Feature flags
    no_kmer: bool,
    no_duplication: bool,
    no_per_tile: bool,
    no_overrep: bool,
    no_adapter: bool,
    no_html: bool,
    no_json: bool,
    fast: bool,
    multiqc: bool,
    // New trimming/filtering options
    poly_g: bool,
    poly_g_min: u8,
    poly_x: bool,
    poly_x_min: u8,
    cut_right: bool,
    cut_front: bool,
    cut_tail: bool,
    window_size: u8,
    window_qual: u8,
    trim_front: u16,
    trim_tail: u16,
    min_quality: u8,
    max_n: Option<u32>,
}

fn parse_args() -> Result<CliConfig, String> {
    let args: Vec<String> = std::env::args().skip(1).collect();

    if args.is_empty() || args.iter().any(|a| a == "--help" || a == "-h") {
        print_help();
        std::process::exit(0);
    }

    if args.iter().any(|a| a == "--version" || a == "-V") {
        println!("biofastq-a {}", env!("CARGO_PKG_VERSION"));
        std::process::exit(0);
    }

    let mut input_files = Vec::new();
    let mut headless = false;
    let mut output_dir = ".".to_string();
    let mut trim = false;
    let mut min_length: u64 = 20;
    let mut custom_adapters: Vec<Vec<u8>> = Vec::new();
    let mut quality_trim: u8 = 0;
    let mut strict = false;
    let mut threads: Option<usize> = None;
    let mut paired_r2: Option<String> = None;
    let mut no_kmer = false;
    let mut no_duplication = false;
    let mut no_per_tile = false;
    let mut no_overrep = false;
    let mut no_adapter = false;
    let mut no_html = false;
    let mut no_json = false;
    let mut fast = false;
    let mut multiqc = false;
    let mut poly_g = false;
    let mut poly_g_min: u8 = 0;
    let mut poly_x = false;
    let mut poly_x_min: u8 = 0;
    let mut cut_right = false;
    let mut cut_front = false;
    let mut cut_tail = false;
    let mut window_size: u8 = 4;
    let mut window_qual: u8 = 20;
    let mut trim_front: u16 = 0;
    let mut trim_tail: u16 = 0;
    let mut min_quality: u8 = 0;
    let mut max_n: Option<u32> = None;
    let mut i = 0;

    while i < args.len() {
        match args[i].as_str() {
            "--headless" => headless = true,
            "--trim" => trim = true,
            "--output-dir" => {
                i += 1;
                if i >= args.len() {
                    return Err("--output-dir requires a value".into());
                }
                output_dir = args[i].clone();
            }
            "--min-length" => {
                i += 1;
                if i >= args.len() {
                    return Err("--min-length requires a value".into());
                }
                min_length = args[i]
                    .parse()
                    .map_err(|_| format!("--min-length must be an integer, got '{}'", args[i]))?;
            }
            "--adapter" => {
                i += 1;
                if i >= args.len() {
                    return Err("--adapter requires a sequence value".into());
                }
                custom_adapters.push(args[i].to_uppercase().into_bytes());
            }
            "--quality-trim" => {
                i += 1;
                if i >= args.len() {
                    return Err("--quality-trim requires a Phred value".into());
                }
                quality_trim = args[i]
                    .parse()
                    .map_err(|_| format!("--quality-trim must be 0-42, got '{}'", args[i]))?;
            }
            "--strict" => strict = true,
            "--poly-g" => {
                poly_g = true;
                if i + 1 < args.len() {
                    if let Ok(n) = args[i + 1].parse::<u8>() {
                        poly_g_min = n; i += 1;
                    }
                }
                if poly_g_min == 0 { poly_g_min = 10; }
            }
            "--poly-x" => {
                poly_x = true;
                if i + 1 < args.len() {
                    if let Ok(n) = args[i + 1].parse::<u8>() {
                        poly_x_min = n; i += 1;
                    }
                }
                if poly_x_min == 0 { poly_x_min = 10; }
            }
            "--cut-right" => cut_right = true,
            "--cut-front" => cut_front = true,
            "--cut-tail"  => cut_tail = true,
            "--window-size" => {
                i += 1;
                if i >= args.len() { return Err("--window-size requires a value".into()); }
                window_size = args[i].parse::<u8>()
                    .map_err(|_| format!("--window-size must be 1-50, got '{}'", args[i]))?;
            }
            "--window-quality" => {
                i += 1;
                if i >= args.len() { return Err("--window-quality requires a value".into()); }
                window_qual = args[i].parse::<u8>()
                    .map_err(|_| format!("--window-quality must be 0-42, got '{}'", args[i]))?;
            }
            "--trim-front" => {
                i += 1;
                if i >= args.len() { return Err("--trim-front requires a value".into()); }
                trim_front = args[i].parse::<u16>()
                    .map_err(|_| format!("--trim-front must be an integer, got '{}'", args[i]))?;
            }
            "--trim-tail" => {
                i += 1;
                if i >= args.len() { return Err("--trim-tail requires a value".into()); }
                trim_tail = args[i].parse::<u16>()
                    .map_err(|_| format!("--trim-tail must be an integer, got '{}'", args[i]))?;
            }
            "--min-quality" => {
                i += 1;
                if i >= args.len() { return Err("--min-quality requires a value".into()); }
                min_quality = args[i].parse::<u8>()
                    .map_err(|_| format!("--min-quality must be 0-42, got '{}'", args[i]))?;
            }
            "--max-n" => {
                i += 1;
                if i >= args.len() { return Err("--max-n requires a value".into()); }
                max_n = Some(args[i].parse::<u32>()
                    .map_err(|_| format!("--max-n must be a non-negative integer, got '{}'", args[i]))?);
            }
            "--no-kmer" => no_kmer = true,
            "--no-duplication" => no_duplication = true,
            "--no-per-tile" => no_per_tile = true,
            "--no-overrep" => no_overrep = true,
            "--no-adapter" => no_adapter = true,
            "--no-html" => no_html = true,
            "--no-json" => no_json = true,
            "--fast" => fast = true,
            "--multiqc" => multiqc = true,
            "--in2" => {
                i += 1;
                if i >= args.len() {
                    return Err("--in2 requires a file path".into());
                }
                paired_r2 = Some(args[i].clone());
            }
            "--threads" => {
                i += 1;
                if i >= args.len() {
                    return Err("--threads requires a value".into());
                }
                let n: usize = args[i]
                    .parse()
                    .map_err(|_| format!("--threads must be a positive integer, got '{}'", args[i]))?;
                if n == 0 {
                    return Err("--threads must be at least 1".into());
                }
                threads = Some(n);
            }
            arg if arg.starts_with("--") => {
                return Err(format!("Unknown option: {}  (use --help for usage)", arg));
            }
            path => input_files.push(path.to_string()),
        }
        i += 1;
    }

    if input_files.is_empty() {
        return Err("No input files specified. Use --help for usage.".into());
    }

    Ok(CliConfig {
        input_files,
        headless,
        output_dir,
        trim,
        min_length,
        custom_adapters,
        quality_trim,
        strict,
        threads,
        paired_r2,
        no_kmer,
        no_duplication,
        no_per_tile,
        no_overrep,
        no_adapter,
        no_html,
        no_json,
        fast,
        multiqc,
        poly_g,
        poly_g_min,
        poly_x,
        poly_x_min,
        cut_right,
        cut_front,
        cut_tail,
        window_size,
        window_qual,
        trim_front,
        trim_tail,
        min_quality,
        max_n,
    })
}

// ---------------------------------------------------------------------------
// Output helpers
// ---------------------------------------------------------------------------

fn print_file_summary(f: &types::FileStats) {
    let (n50, n90) = f.compute_n50_n90();
    let fname = std::path::Path::new(&f.file_path)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(&f.file_path);

    println!(
        "  File:            {}\n  \
         Reads:           {}\n  \
         Total bases:     {}\n  \
         Avg length:      {:.1} bp\n  \
         Min / Max:       {} / {} bp\n  \
         N50 / N90:       {} / {} bp\n  \
         GC content:      {:.2}%\n  \
         Avg quality:     Q{:.1}\n  \
         Q20 bases:       {:.1}%\n  \
         Q30 bases:       {:.1}%\n  \
         Adapter cont:    {:.2}%\n  \
         Dup rate (est):  {:.1}%",
        fname,
        format_number(f.read_count),
        types::format_bases(f.total_bases),
        f.avg_length(),
        f.effective_min_length(), f.max_length,
        n50, n90,
        f.gc_content(),
        f.avg_quality(),
        f.q20_pct(),
        f.q30_pct(),
        f.adapter_pct(),
        f.dup_rate_pct,
    );
    if f.trimmed_reads > 0 {
        println!("  Trimmed reads:   {} ({:.1}%)", format_number(f.trimmed_reads), f.trimmed_pct());
        if let Some(ref tp) = f.trim_output_path {
            println!("  Trim output:     {}", tp);
        }
    }
    if !f.per_tile_quality.is_empty() {
        println!("  Illumina tiles:  {}", f.per_tile_quality.len());
    }
    println!();
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let cfg = match parse_args() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    };

    // Configure rayon thread pool
    if let Some(n) = cfg.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .unwrap_or_else(|e| eprintln!("Warning: could not set thread count: {}", e));
    }

    // Validate input files
    for path in &cfg.input_files {
        match std::fs::metadata(path) {
            Ok(m) if m.is_file() => {}
            Ok(_) => {
                eprintln!("Error: '{}' is not a regular file.", path);
                std::process::exit(1);
            }
            Err(e) => {
                eprintln!("Error: Cannot open '{}': {}", path, e);
                std::process::exit(1);
            }
        }
    }

    // Ensure output directory exists
    if let Err(e) = std::fs::create_dir_all(&cfg.output_dir) {
        eprintln!(
            "Error: Cannot create output directory '{}': {}",
            cfg.output_dir, e
        );
        std::process::exit(1);
    }

    // Validate --in2 is only used with a single R1 input
    if cfg.paired_r2.is_some() && cfg.input_files.len() != 1 {
        eprintln!("Error: --in2 requires exactly one R1 input file.");
        std::process::exit(1);
    }

    let flags = FeatureFlags {
        kmer_analysis:    !(cfg.no_kmer        || cfg.fast),
        duplication_check:!(cfg.no_duplication || cfg.fast),
        per_tile_quality: !(cfg.no_per_tile    || cfg.fast),
        overrep_sequences:!(cfg.no_overrep     || cfg.fast),
        adapter_detection: !cfg.no_adapter,
        html_report:       !cfg.no_html,
        json_report:       !cfg.no_json,
    };

    let mut process_config = ProcessConfig {
        trim_output: cfg.trim,
        min_length: cfg.min_length,
        output_dir: cfg.output_dir.clone(),
        custom_adapters: cfg.custom_adapters.clone(),
        quality_trim_threshold: cfg.quality_trim,
        strict: cfg.strict,
        paired_end_r2: cfg.paired_r2.clone(),
        flags,
        ..ProcessConfig::default()
    };
    process_config.poly_g_min_len = if cfg.poly_g { cfg.poly_g_min } else { 0 };
    process_config.poly_x_min_len = if cfg.poly_x { cfg.poly_x_min } else { 0 };
    process_config.cut_right_window = if cfg.cut_right { cfg.window_size } else { 0 };
    process_config.cut_right_qual   = cfg.window_qual;
    process_config.cut_front_window = if cfg.cut_front { cfg.window_size } else { 0 };
    process_config.cut_front_qual   = cfg.window_qual;
    process_config.cut_tail_window  = if cfg.cut_tail { cfg.window_size } else { 0 };
    process_config.cut_tail_qual    = cfg.window_qual;
    process_config.trim_front_bases = cfg.trim_front;
    process_config.trim_tail_bases  = cfg.trim_tail;
    process_config.min_avg_quality  = cfg.min_quality;
    process_config.max_n_bases      = cfg.max_n;

    let n_files = cfg.input_files.len();
    let state = Arc::new(Mutex::new(SharedState::new(n_files)));

    // Panic hook — always restore terminal
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let _ = disable_raw_mode();
        let _ = execute!(stdout(), LeaveAlternateScreen);
        default_hook(info);
    }));

    // Stash flags before process_config is moved into the thread
    let report_flags = process_config.flags.clone();

    // Spawn processing thread
    let state_worker = Arc::clone(&state);
    let files = cfg.input_files.clone();
    let handle = thread::spawn(move || {
        analysis::process_files(files, state_worker, process_config);
    });

    if cfg.headless {
        let _ = handle.join();
    } else {
        let mut terminal = match tui::setup_terminal() {
            Ok(t) => t,
            Err(e) => {
                eprintln!("Failed to initialise terminal: {}", e);
                std::process::exit(1);
            }
        };

        let tick = Duration::from_millis(80);
        loop {
            let snap = state.lock().unwrap().clone();
            let _ = terminal.draw(|f| tui::render(f, &snap));

            if crossterm::event::poll(tick).unwrap_or(false) {
                if let Ok(Event::Key(key)) = event::read() {
                    if key.code == KeyCode::Char('q') || key.code == KeyCode::Esc {
                        break;
                    }
                }
            }

            if matches!(
                snap.status,
                ProcessingStatus::Completed | ProcessingStatus::Error(_)
            ) {
                thread::sleep(Duration::from_millis(600));
                break;
            }
        }

        tui::restore_terminal(&mut terminal);
        let _ = handle.join();
    }

    // ----- Report export & terminal summary -----
    let snap = state.lock().unwrap().clone();

    match &snap.status {
        ProcessingStatus::Completed => {
            println!("\nAnalysis complete!\n");
            for f in snap.all_files() {
                print_file_summary(f);
            }
            if report_flags.html_report {
                match report::export_html(&snap, &cfg.output_dir, &report_flags) {
                    Ok(path) => println!("HTML report:  {}", path),
                    Err(e) => eprintln!("Warning: HTML report failed: {}", e),
                }
            }
            if report_flags.json_report {
                match report::export_json(&snap, &cfg.output_dir) {
                    Ok(path) => println!("JSON report:  {}", path),
                    Err(e) => eprintln!("Warning: JSON report failed: {}", e),
                }
            }
            if cfg.multiqc {
                match report::export_multiqc(&snap, &cfg.output_dir) {
                    Ok(path) => println!("MultiQC data: {}", path),
                    Err(e) => eprintln!("Warning: MultiQC export failed: {}", e),
                }
            }
            let elapsed = snap.elapsed_secs();
            let total_reads: u64 = snap.all_files().iter().map(|f| f.read_count).sum();
            let total_bytes: u64 = snap.all_files().iter().map(|f| f.file_size).sum();
            let rps = if elapsed > 0.0 { total_reads as f64 / elapsed } else { 0.0 };
            let mbs = if elapsed > 0.0 { total_bytes as f64 / 1_048_576.0 / elapsed } else { 0.0 };
            println!("\nProcessing time:  {:.1}s", elapsed);
            println!("Throughput:       {:.0} reads/s  |  {:.0} MB/s", rps, mbs);
        }
        ProcessingStatus::Error(e) => {
            eprintln!("\nAnalysis failed: {}", e);
            std::process::exit(1);
        }
        ProcessingStatus::Running => {
            println!("\nAnalysis was interrupted.");
        }
    }
}
