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

use types::{format_number, ProcessingStatus, SharedState};

// ---------------------------------------------------------------------------
// CLI help / argument parsing
// ---------------------------------------------------------------------------

fn print_help() {
    eprintln!(
        r#"
BioFastq-A v2.0.0 — High-performance FASTQ/FASTA quality analysis

USAGE:
  biofastq-a [OPTIONS] <file> [<file2> ...]

ARGUMENTS:
  <file>     FASTQ or FASTA file(s) to analyse. Gzipped files (.gz) are
             supported automatically.

OPTIONS:
  --headless         Run without interactive TUI (batch/script mode)
  --output-dir <dir> Directory for report files (default: current directory)
  --help, -h         Show this help message

OUTPUT:
  <stem>_report.html — Self-contained HTML report with interactive charts
  <stem>_report.json — Machine-readable JSON report

EXAMPLES:
  biofastq-a sample.fastq
  biofastq-a *.fastq.gz --headless --output-dir ./reports
  biofastq-a reads.fasta --headless
"#
    );
}

struct Config {
    input_files: Vec<String>,
    headless: bool,
    output_dir: String,
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = std::env::args().skip(1).collect();

    if args.is_empty()
        || args.iter().any(|a| a == "--help" || a == "-h")
    {
        print_help();
        std::process::exit(0);
    }

    let mut input_files = Vec::new();
    let mut headless = false;
    let mut output_dir = ".".to_string();
    let mut i = 0;

    while i < args.len() {
        match args[i].as_str() {
            "--headless" => headless = true,
            "--output-dir" => {
                i += 1;
                if i >= args.len() {
                    return Err("--output-dir requires a value".into());
                }
                output_dir = args[i].clone();
            }
            arg if arg.starts_with("--") => {
                return Err(format!("Unknown option: {}", arg));
            }
            path => input_files.push(path.to_string()),
        }
        i += 1;
    }

    if input_files.is_empty() {
        return Err("No input files specified. Use --help for usage.".into());
    }

    Ok(Config {
        input_files,
        headless,
        output_dir,
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

    // Validate that all input files exist before starting
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

    let n_files = cfg.input_files.len();
    let state = Arc::new(Mutex::new(SharedState::new(n_files)));

    // Panic hook — always restore the terminal
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
        analysis::process_files(files, state_worker);
    });

    if cfg.headless {
        // ----- headless mode -----
        let _ = handle.join();
    } else {
        // ----- interactive TUI mode -----
        let mut terminal = match tui::setup_terminal() {
            Ok(t) => t,
            Err(e) => {
                eprintln!("Failed to set up terminal: {}", e);
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

            // Exit UI automatically when all processing is done
            if matches!(
                snap.status,
                ProcessingStatus::Completed | ProcessingStatus::Error(_)
            ) {
                // Give the user a moment to see the final state
                thread::sleep(Duration::from_millis(500));
                break;
            }
        }

        tui::restore_terminal(&mut terminal);
        let _ = handle.join();
    }

    // ----- Report export -----
    let snap = state.lock().unwrap().clone();

    match &snap.status {
        ProcessingStatus::Completed => {
            println!("\nAnalysis complete!\n");

            // Print summary for each file
            for f in snap.all_files() {
                let (n50, n90) = f.compute_n50_n90();
                let fname = std::path::Path::new(&f.file_path)
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or(&f.file_path);
                println!("  File:          {}", fname);
                println!("  Reads:         {}", format_number(f.read_count));
                println!("  Total bases:   {}", types::format_bases(f.total_bases));
                println!("  Avg length:    {:.1} bp", f.avg_length());
                println!(
                    "  Min/Max:       {}/{} bp",
                    f.effective_min_length(),
                    f.max_length
                );
                println!("  N50:           {} bp", n50);
                println!("  N90:           {} bp", n90);
                println!("  GC content:    {:.2}%", f.gc_content());
                println!("  Avg quality:   Q{:.1}", f.avg_quality());
                println!("  Q20 reads:     {:.1}%", f.q20_pct());
                println!("  Q30 reads:     {:.1}%", f.q30_pct());
                println!("  Adapter cont:  {:.2}%", f.adapter_pct());
                println!();
            }

            // Export HTML report
            match report::export_html(&snap, &cfg.output_dir) {
                Ok(path) => println!("HTML report: {}", path),
                Err(e) => eprintln!("Warning: Could not write HTML report: {}", e),
            }

            // Export JSON report
            match report::export_json(&snap, &cfg.output_dir) {
                Ok(path) => println!("JSON report: {}", path),
                Err(e) => eprintln!("Warning: Could not write JSON report: {}", e),
            }

            println!("\nProcessing time: {:.1}s", snap.elapsed_secs());
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
