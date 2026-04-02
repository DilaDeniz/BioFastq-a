use std::collections::HashMap;
use std::env;
use std::fs;
use std::io::{self, stdout};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use crossterm::event::{self, Event, KeyCode};
use crossterm::execute;
use crossterm::terminal::{
    disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen,
};
use needletail::parse_fastx_file;
use needletail::Sequence;
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Gauge, List, ListItem, Paragraph};
use ratatui::Terminal;
use rayon::prelude::*;
use serde::Serialize;

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

#[derive(Clone, PartialEq)]
enum ProcessingStatus {
    Running,
    Completed,
    Error(String),
}

#[derive(Clone)]
enum LogLevel {
    Info,
    Success,
    Error,
}

#[derive(Clone)]
struct LogEntry {
    timestamp: Duration,
    level: LogLevel,
    message: String,
}

#[derive(Clone)]
struct SharedState {
    read_count: u64,
    total_bases: u64,
    gc_count: u64,
    total_quality_sum: u64,
    total_quality_bases: u64,
    kmer_counts: HashMap<[u8; 4], u64>,
    bytes_processed: u64,
    file_size: u64,
    start_time: Instant,
    status: ProcessingStatus,
    log_messages: Vec<LogEntry>,
}

impl SharedState {
    fn new(file_size: u64) -> Self {
        Self {
            read_count: 0,
            total_bases: 0,
            gc_count: 0,
            total_quality_sum: 0,
            total_quality_bases: 0,
            kmer_counts: HashMap::new(),
            bytes_processed: 0,
            file_size,
            start_time: Instant::now(),
            status: ProcessingStatus::Running,
            log_messages: Vec::new(),
        }
    }

    fn log(&mut self, level: LogLevel, message: String) {
        let elapsed = self.start_time.elapsed();
        self.log_messages.push(LogEntry {
            timestamp: elapsed,
            level,
            message,
        });
        if self.log_messages.len() > 1000 {
            self.log_messages.remove(0);
        }
    }
}

#[derive(Serialize)]
struct AnalysisReport {
    read_count: u64,
    average_sequence_length: f64,
    gc_content_percent: f64,
    average_phred_quality: f64,
    top_5_kmers: Vec<(String, u64)>,
    processing_time_seconds: f64,
}

// ---------------------------------------------------------------------------
// Processing engine
// ---------------------------------------------------------------------------

fn process_file(path: String, state: Arc<Mutex<SharedState>>) {
    {
        let mut s = state.lock().unwrap();
        let size_gb = s.file_size as f64 / 1_073_741_824.0;
        s.log(
            LogLevel::Info,
            format!("Opening file: {} ({:.2} GB)", path, size_gb),
        );
        let cpus = rayon::current_num_threads();
        s.log(
            LogLevel::Info,
            format!("Rayon thread pool: {} threads active", cpus),
        );
    }

    let reader = match parse_fastx_file(&path) {
        Ok(r) => r,
        Err(e) => {
            let mut s = state.lock().unwrap();
            s.status = ProcessingStatus::Error(format!("Parse error: {}", e));
            s.log(LogLevel::Error, format!("Failed to open file: {}", e));
            return;
        }
    };

    let flush_interval: u64 = 10_000;
    let batch_size: usize = 50_000;

    let mut local_read_count: u64 = 0;
    let mut local_total_bases: u64 = 0;
    let mut local_gc_count: u64 = 0;
    let mut local_quality_sum: u64 = 0;
    let mut local_quality_bases: u64 = 0;
    let mut local_bytes: u64 = 0;
    let mut batch: Vec<Vec<u8>> = Vec::with_capacity(batch_size);
    let mut local_kmer_counts: HashMap<[u8; 4], u64> = HashMap::new();
    let mut total_flushed: u64 = 0;

    let mut reader = reader;

    while let Some(result) = reader.next() {
        let record = match result {
            Ok(r) => r,
            Err(e) => {
                let mut s = state.lock().unwrap();
                s.log(LogLevel::Error, format!("Skipping bad record: {}", e));
                continue;
            }
        };

        let seq_bytes = record.normalize(false);

        local_read_count += 1;
        local_total_bases += seq_bytes.len() as u64;

        // Estimate bytes processed (header + seq + qual + format overhead)
        local_bytes += seq_bytes.len() as u64 + record.id().len() as u64 + 10;

        // GC count
        for &b in seq_bytes.iter() {
            if b == b'G' || b == b'C' {
                local_gc_count += 1;
            }
        }

        // Quality scores (FASTQ only)
        if let Some(qual) = record.qual() {
            local_bytes += qual.len() as u64;
            for &q in qual.iter() {
                local_quality_sum += q.saturating_sub(33) as u64;
            }
            local_quality_bases += qual.len() as u64;
        }

        // Collect sequences for batched parallel k-mer counting
        if seq_bytes.len() >= 4 {
            batch.push(seq_bytes.to_vec());
        }

        // Flush stats periodically
        if local_read_count.is_multiple_of(flush_interval) {
            // Process k-mer batch with rayon if big enough
            if batch.len() >= batch_size {
                let batch_kmers = count_kmers_parallel(&batch);
                merge_kmer_counts(&mut local_kmer_counts, &batch_kmers);
                batch.clear();
            }

            let mut s = state.lock().unwrap();
            s.read_count += local_read_count;
            s.total_bases += local_total_bases;
            s.gc_count += local_gc_count;
            s.total_quality_sum += local_quality_sum;
            s.total_quality_bases += local_quality_bases;
            s.bytes_processed += local_bytes;
            merge_kmer_counts(&mut s.kmer_counts, &local_kmer_counts);

            total_flushed += local_read_count;
            if total_flushed.is_multiple_of(100_000) {
                let count = s.read_count;
                s.log(
                    LogLevel::Success,
                    format!("Processed {} reads...", count),
                );
            }

            local_read_count = 0;
            local_total_bases = 0;
            local_gc_count = 0;
            local_quality_sum = 0;
            local_quality_bases = 0;
            local_bytes = 0;
            local_kmer_counts.clear();
        }
    }

    // Final k-mer batch
    if !batch.is_empty() {
        let batch_kmers = count_kmers_parallel(&batch);
        merge_kmer_counts(&mut local_kmer_counts, &batch_kmers);
    }

    // Final flush
    let mut s = state.lock().unwrap();
    s.read_count += local_read_count;
    s.total_bases += local_total_bases;
    s.gc_count += local_gc_count;
    s.total_quality_sum += local_quality_sum;
    s.total_quality_bases += local_quality_bases;
    s.bytes_processed += local_bytes;
    merge_kmer_counts(&mut s.kmer_counts, &local_kmer_counts);
    s.status = ProcessingStatus::Completed;
    let count = s.read_count;
    s.log(
        LogLevel::Success,
        format!("Analysis complete! {} total reads processed.", count),
    );
}

fn count_kmers_parallel(sequences: &[Vec<u8>]) -> HashMap<[u8; 4], u64> {
    sequences
        .par_iter()
        .fold(HashMap::new, |mut map: HashMap<[u8; 4], u64>, seq| {
            for window in seq.windows(4) {
                if window
                    .iter()
                    .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
                {
                    let key: [u8; 4] = [window[0], window[1], window[2], window[3]];
                    *map.entry(key).or_insert(0) += 1;
                }
            }
            map
        })
        .reduce(HashMap::new, |mut a, b| {
            for (k, v) in b {
                *a.entry(k).or_insert(0) += v;
            }
            a
        })
}

fn merge_kmer_counts(target: &mut HashMap<[u8; 4], u64>, source: &HashMap<[u8; 4], u64>) {
    for (k, v) in source {
        *target.entry(*k).or_insert(0) += v;
    }
}

// ---------------------------------------------------------------------------
// JSON export
// ---------------------------------------------------------------------------

fn export_report(state: &SharedState, input_path: &str) -> io::Result<String> {
    let mut top_kmers: Vec<(String, u64)> = state
        .kmer_counts
        .iter()
        .map(|(k, &v)| (String::from_utf8_lossy(k).to_string(), v))
        .collect();
    top_kmers.sort_by(|a, b| b.1.cmp(&a.1));
    top_kmers.truncate(5);

    let report = AnalysisReport {
        read_count: state.read_count,
        average_sequence_length: if state.read_count > 0 {
            state.total_bases as f64 / state.read_count as f64
        } else {
            0.0
        },
        gc_content_percent: if state.total_bases > 0 {
            (state.gc_count as f64 / state.total_bases as f64) * 100.0
        } else {
            0.0
        },
        average_phred_quality: if state.total_quality_bases > 0 {
            state.total_quality_sum as f64 / state.total_quality_bases as f64
        } else {
            0.0
        },
        top_5_kmers: top_kmers,
        processing_time_seconds: state.start_time.elapsed().as_secs_f64(),
    };

    let json = serde_json::to_string_pretty(&report).map_err(|e| {
        io::Error::other(e)
    })?;

    let stem = Path::new(input_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("output");
    let output_path = format!("{}_report.json", stem);
    fs::write(&output_path, &json)?;
    Ok(output_path)
}

// ---------------------------------------------------------------------------
// Terminal setup / teardown
// ---------------------------------------------------------------------------

fn setup_terminal() -> io::Result<Terminal<CrosstermBackend<io::Stdout>>> {
    enable_raw_mode()?;
    let mut stdout = stdout();
    execute!(stdout, EnterAlternateScreen)?;
    let backend = CrosstermBackend::new(stdout);
    Terminal::new(backend)
}

fn restore_terminal(terminal: &mut Terminal<CrosstermBackend<io::Stdout>>) {
    let _ = disable_raw_mode();
    let _ = execute!(terminal.backend_mut(), LeaveAlternateScreen);
    let _ = terminal.show_cursor();
}

// ---------------------------------------------------------------------------
// UI rendering
// ---------------------------------------------------------------------------

fn format_number(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

fn format_duration(d: Duration) -> String {
    let secs = d.as_secs();
    let h = secs / 3600;
    let m = (secs % 3600) / 60;
    let s = secs % 60;
    format!("{:02}:{:02}:{:02}", h, m, s)
}

fn ui(frame: &mut ratatui::Frame, snap: &SharedState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(8),
            Constraint::Length(8),
        ])
        .split(frame.area());

    // -- Top bar --
    let user = env::var("USER").unwrap_or_else(|_| "Researcher".to_string());
    let title_line = Line::from(vec![
        Span::styled(
            " Dilo Science Terminal v1.0 ",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(" | ", Style::default().fg(Color::DarkGray)),
        Span::styled(
            format!("Licensed to: {} ", user),
            Style::default()
                .fg(Color::Magenta)
                .add_modifier(Modifier::BOLD),
        ),
    ]);
    let top_bar = Paragraph::new(title_line)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .border_style(Style::default().fg(Color::Cyan))
                .title(" BIOFASTQ-A ")
                .title_style(Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)),
        )
        .alignment(ratatui::layout::Alignment::Center);
    frame.render_widget(top_bar, chunks[0]);

    // -- Main area: left + right panels --
    let main_chunks = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(50), Constraint::Percentage(50)])
        .split(chunks[1]);

    render_processing_panel(frame, main_chunks[0], snap);
    render_results_panel(frame, main_chunks[1], snap);

    // -- Bottom log --
    render_log_panel(frame, chunks[2], snap);
}

fn render_processing_panel(frame: &mut ratatui::Frame, area: ratatui::layout::Rect, snap: &SharedState) {
    let inner_chunks = Layout::default()
        .direction(Direction::Vertical)
        .margin(1)
        .constraints([
            Constraint::Length(2),
            Constraint::Length(3),
            Constraint::Length(2),
            Constraint::Min(0),
        ])
        .split(area);

    let block = Block::default()
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Cyan))
        .title(" PROCESSING ")
        .title_style(Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD));
    frame.render_widget(block, area);

    // Reads/sec
    let elapsed = snap.start_time.elapsed().as_secs_f64();
    let reads_per_sec = if elapsed > 0.0 {
        snap.read_count as f64 / elapsed
    } else {
        0.0
    };
    let status_text = match &snap.status {
        ProcessingStatus::Running => "RUNNING",
        ProcessingStatus::Completed => "COMPLETED",
        ProcessingStatus::Error(_) => "ERROR",
    };
    let status_color = match &snap.status {
        ProcessingStatus::Running => Color::Yellow,
        ProcessingStatus::Completed => Color::LightGreen,
        ProcessingStatus::Error(_) => Color::Red,
    };

    let speed_line = Line::from(vec![
        Span::styled("Status: ", Style::default().fg(Color::White)),
        Span::styled(
            status_text,
            Style::default().fg(status_color).add_modifier(Modifier::BOLD),
        ),
        Span::raw("  "),
        Span::styled("Speed: ", Style::default().fg(Color::White)),
        Span::styled(
            format!("{} reads/sec", format_number(reads_per_sec as u64)),
            Style::default().fg(Color::LightGreen),
        ),
    ]);
    frame.render_widget(Paragraph::new(speed_line), inner_chunks[0]);

    // Progress bar
    let progress = if snap.file_size > 0 {
        (snap.bytes_processed as f64 / snap.file_size as f64).min(1.0)
    } else {
        0.0
    };
    let pct_text = format!("{:.1}%", progress * 100.0);
    let gauge = Gauge::default()
        .block(Block::default().title(" Progress ").borders(Borders::ALL).border_style(Style::default().fg(Color::DarkGray)))
        .gauge_style(
            Style::default()
                .fg(Color::Cyan)
                .bg(Color::DarkGray),
        )
        .ratio(progress)
        .label(pct_text);
    frame.render_widget(gauge, inner_chunks[1]);

    // Elapsed time
    let time_line = Line::from(vec![
        Span::styled("Elapsed: ", Style::default().fg(Color::White)),
        Span::styled(
            format_duration(snap.start_time.elapsed()),
            Style::default().fg(Color::LightGreen),
        ),
    ]);
    frame.render_widget(Paragraph::new(time_line), inner_chunks[2]);
}

fn render_results_panel(frame: &mut ratatui::Frame, area: ratatui::layout::Rect, snap: &SharedState) {
    let inner_chunks = Layout::default()
        .direction(Direction::Vertical)
        .margin(1)
        .constraints([
            Constraint::Length(2),
            Constraint::Length(2),
            Constraint::Length(2),
            Constraint::Length(2),
            Constraint::Min(0),
        ])
        .split(area);

    let block = Block::default()
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Magenta))
        .title(" ANALYSIS RESULTS ")
        .title_style(
            Style::default()
                .fg(Color::Magenta)
                .add_modifier(Modifier::BOLD),
        );
    frame.render_widget(block, area);

    // Read count
    let read_line = Line::from(vec![
        Span::styled("  Read Count:   ", Style::default().fg(Color::White)),
        Span::styled(
            format_number(snap.read_count),
            Style::default()
                .fg(Color::LightGreen)
                .add_modifier(Modifier::BOLD),
        ),
    ]);
    frame.render_widget(Paragraph::new(read_line), inner_chunks[0]);

    // GC content
    let gc_pct = if snap.total_bases > 0 {
        (snap.gc_count as f64 / snap.total_bases as f64) * 100.0
    } else {
        0.0
    };
    let gc_line = Line::from(vec![
        Span::styled("  GC Content:   ", Style::default().fg(Color::White)),
        Span::styled(
            format!("{:.2}%", gc_pct),
            Style::default()
                .fg(Color::LightGreen)
                .add_modifier(Modifier::BOLD),
        ),
    ]);
    frame.render_widget(Paragraph::new(gc_line), inner_chunks[1]);

    // Avg quality
    let avg_qual = if snap.total_quality_bases > 0 {
        format!(
            "{:.1}",
            snap.total_quality_sum as f64 / snap.total_quality_bases as f64
        )
    } else {
        "N/A".to_string()
    };
    let qual_line = Line::from(vec![
        Span::styled("  Avg Quality:  ", Style::default().fg(Color::White)),
        Span::styled(
            avg_qual,
            Style::default()
                .fg(Color::LightGreen)
                .add_modifier(Modifier::BOLD),
        ),
    ]);
    frame.render_widget(Paragraph::new(qual_line), inner_chunks[2]);

    // Avg length
    let avg_len = if snap.read_count > 0 {
        format!("{:.1} bp", snap.total_bases as f64 / snap.read_count as f64)
    } else {
        "N/A".to_string()
    };
    let len_line = Line::from(vec![
        Span::styled("  Avg Length:   ", Style::default().fg(Color::White)),
        Span::styled(
            avg_len,
            Style::default()
                .fg(Color::LightGreen)
                .add_modifier(Modifier::BOLD),
        ),
    ]);
    frame.render_widget(Paragraph::new(len_line), inner_chunks[3]);
}

fn render_log_panel(frame: &mut ratatui::Frame, area: ratatui::layout::Rect, snap: &SharedState) {
    let block = Block::default()
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::DarkGray))
        .title(" EVENT LOG ")
        .title_style(Style::default().fg(Color::White).add_modifier(Modifier::BOLD));

    let max_visible = area.height.saturating_sub(2) as usize;
    let start = if snap.log_messages.len() > max_visible {
        snap.log_messages.len() - max_visible
    } else {
        0
    };

    let items: Vec<ListItem> = snap.log_messages[start..]
        .iter()
        .map(|entry| {
            let ts = format_duration(entry.timestamp);
            let (prefix, color) = match entry.level {
                LogLevel::Info => ("INFO", Color::White),
                LogLevel::Success => (" OK ", Color::LightGreen),
                LogLevel::Error => ("ERR ", Color::Red),
            };
            let line = Line::from(vec![
                Span::styled(
                    format!(" [{}] ", ts),
                    Style::default().fg(Color::DarkGray),
                ),
                Span::styled(
                    format!("{} ", prefix),
                    Style::default().fg(color).add_modifier(Modifier::BOLD),
                ),
                Span::styled(&entry.message, Style::default().fg(color)),
            ]);
            ListItem::new(line)
        })
        .collect();

    let list = List::new(items).block(block);
    frame.render_widget(list, area);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = env::args().collect();
    let headless = args.iter().any(|a| a == "--headless");
    let file_args: Vec<&String> = args.iter().skip(1).filter(|a| *a != "--headless").collect();

    if file_args.is_empty() {
        eprintln!("Usage: biofastq-a <input.fastq|input.fasta> [--headless]");
        eprintln!("  Analyzes FASTQ/FASTA files with real-time TUI dashboard.");
        std::process::exit(1);
    }

    let input_path = file_args[0];

    // Validate file exists
    let metadata = match fs::metadata(input_path) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("Error: Cannot open '{}': {}", input_path, e);
            std::process::exit(1);
        }
    };

    if !metadata.is_file() {
        eprintln!("Error: '{}' is not a file.", input_path);
        std::process::exit(1);
    }

    let file_size = metadata.len();
    let state = Arc::new(Mutex::new(SharedState::new(file_size)));

    // Set up panic hook to restore terminal
    let default_panic = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let _ = disable_raw_mode();
        let _ = execute!(stdout(), LeaveAlternateScreen);
        default_panic(info);
    }));

    // Spawn processing thread
    let state_clone = Arc::clone(&state);
    let path_clone = input_path.clone();
    let processing_handle = thread::spawn(move || {
        process_file(path_clone, state_clone);
    });

    if headless {
        // Headless mode: just wait for processing to finish
        let _ = processing_handle.join();
    } else {
        // Set up TUI
        let mut terminal = match setup_terminal() {
            Ok(t) => t,
            Err(e) => {
                eprintln!("Error: Failed to initialize terminal: {}", e);
                std::process::exit(1);
            }
        };

        let tick_rate = Duration::from_millis(80);

        // Main UI loop
        loop {
            // Clone state snapshot for rendering
            let snap = {
                let s = state.lock().unwrap();
                s.clone()
            };

            let _ = terminal.draw(|f| ui(f, &snap));

            if crossterm::event::poll(tick_rate).unwrap_or(false) {
                if let Ok(Event::Key(key)) = event::read() {
                    if key.code == KeyCode::Char('q') || key.code == KeyCode::Esc {
                        break;
                    }
                }
            }
        }

        // Restore terminal
        restore_terminal(&mut terminal);

        // Wait for processing thread to finish
        let _ = processing_handle.join();
    }

    // Export report if completed
    let snap = state.lock().unwrap();
    match &snap.status {
        ProcessingStatus::Completed => {
            match export_report(&snap, input_path) {
                Ok(path) => println!("\nReport saved to: {}", path),
                Err(e) => eprintln!("\nFailed to write report: {}", e),
            }
            println!("Total reads: {}", format_number(snap.read_count));
            let gc_pct = if snap.total_bases > 0 {
                (snap.gc_count as f64 / snap.total_bases as f64) * 100.0
            } else {
                0.0
            };
            println!("GC Content: {:.2}%", gc_pct);
        }
        ProcessingStatus::Error(e) => {
            eprintln!("\nAnalysis failed: {}", e);
        }
        ProcessingStatus::Running => {
            println!("\nAnalysis cancelled.");
        }
    }
}
