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

use types::{format_number, ProcessConfig, ProcessingStatus, SharedState};

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

fn print_help() {
    eprintln!(
        r#"
BioFastq-A v2.1.0 — High-performance FASTQ/FASTA quality analysis

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
  --strict            Abort on first malformed record (default: skip and warn)
  --version, -V       Print version and exit
  --help, -h          Show this help message

OUTPUT:
  <stem>_report.html      Self-contained HTML report with interactive charts
  <stem>_report.json      Machine-readable JSON report
  <stem>_trimmed.fastq.gz Adapter-trimmed reads (when --trim is used)

EXAMPLES:
  biofastq-a sample.fastq
  biofastq-a reads.fastq.gz --trim --output-dir ./qc
  biofastq-a *.fastq.gz --headless --output-dir ./reports
  biofastq-a reads.fasta --headless
"#
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
    })
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

    let process_config = ProcessConfig {
        trim_output: cfg.trim,
        min_length: cfg.min_length,
        output_dir: cfg.output_dir.clone(),
        custom_adapters: cfg.custom_adapters.clone(),
        quality_trim_threshold: cfg.quality_trim,
        strict: cfg.strict,
    };

    let n_files = cfg.input_files.len();
    let state = Arc::new(Mutex::new(SharedState::new(n_files)));

    // Panic hook — always restore terminal
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let _ = disable_raw_mode();
        let _ = execute!(stdout(), LeaveAlternateScreen);
        default_hook(info);
    }));

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
                let (n50, n90) = f.compute_n50_n90();
                let fname = std::path::Path::new(&f.file_path)
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or(&f.file_path);

                println!("  File:            {}", fname);
                println!("  Reads:           {}", format_number(f.read_count));
                println!("  Total bases:     {}", types::format_bases(f.total_bases));
                println!("  Avg length:      {:.1} bp", f.avg_length());
                println!(
                    "  Min / Max:       {} / {} bp",
                    f.effective_min_length(),
                    f.max_length
                );
                println!("  N50:             {} bp", n50);
                println!("  N90:             {} bp", n90);
                println!("  GC content:      {:.2}%", f.gc_content());
                println!("  Avg quality:     Q{:.1}", f.avg_quality());
                println!("  Q20 bases:       {:.1}%", f.q20_pct());
                println!("  Q30 bases:       {:.1}%", f.q30_pct());
                println!("  Adapter cont:    {:.2}%", f.adapter_pct());
                println!("  Dup rate (est):  {:.1}%", f.dup_rate_pct);
                if f.trimmed_reads > 0 {
                    println!(
                        "  Trimmed reads:   {} ({:.1}%)",
                        format_number(f.trimmed_reads),
                        f.trimmed_pct()
                    );
                    if let Some(ref tp) = f.trim_output_path {
                        println!("  Trim output:     {}", tp);
                    }
                }
                if !f.per_tile_quality.is_empty() {
                    println!("  Illumina tiles:  {}", f.per_tile_quality.len());
                }
                println!();
            }

            match report::export_html(&snap, &cfg.output_dir) {
                Ok(path) => println!("HTML report:  {}", path),
                Err(e) => eprintln!("Warning: HTML report failed: {}", e),
            }
            match report::export_json(&snap, &cfg.output_dir) {
                Ok(path) => println!("JSON report:  {}", path),
                Err(e) => eprintln!("Warning: JSON report failed: {}", e),
            }

            let elapsed = snap.elapsed_secs();
            let total_reads: u64 = snap.all_files().iter().map(|f| f.read_count).sum();
            let total_bytes: u64 = snap.all_files().iter().map(|f| f.file_size).sum();
            let reads_per_sec = if elapsed > 0.0 { total_reads as f64 / elapsed } else { 0.0 };
            let mb_per_sec = if elapsed > 0.0 { total_bytes as f64 / 1_048_576.0 / elapsed } else { 0.0 };

            println!("Processing time:  {:.1}s", elapsed);
            println!("Throughput:       {:.0} reads/s  |  {:.0} MB/s", reads_per_sec, mb_per_sec);
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
