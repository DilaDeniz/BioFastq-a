use std::io;
use std::io::stdout;

use crossterm::execute;
use crossterm::terminal::{
    disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen,
};
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Alignment, Constraint, Direction, Layout};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Gauge, List, ListItem, Paragraph};
use ratatui::Terminal;

use crate::types::{format_bases, format_duration, format_number, LogLevel, ProcessingStatus, SharedState};

pub type Term = Terminal<CrosstermBackend<io::Stdout>>;

pub fn setup_terminal() -> io::Result<Term> {
    enable_raw_mode()?;
    let mut stdout = stdout();
    execute!(stdout, EnterAlternateScreen)?;
    let backend = CrosstermBackend::new(stdout);
    Terminal::new(backend)
}

pub fn restore_terminal(terminal: &mut Term) {
    let _ = disable_raw_mode();
    let _ = execute!(terminal.backend_mut(), LeaveAlternateScreen);
    let _ = terminal.show_cursor();
}

// ---------------------------------------------------------------------------
// Main UI entry point
// ---------------------------------------------------------------------------

pub fn render(frame: &mut ratatui::Frame, snap: &SharedState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3), // title bar
            Constraint::Min(14),   // main panels
            Constraint::Length(8), // log
        ])
        .split(frame.area());

    render_title_bar(frame, chunks[0], snap);

    let main = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(38), Constraint::Percentage(62)])
        .split(chunks[1]);

    render_progress_panel(frame, main[0], snap);
    render_metrics_panel(frame, main[1], snap);
    render_log_panel(frame, chunks[2], snap);
}

// ---------------------------------------------------------------------------
// Title bar
// ---------------------------------------------------------------------------

fn render_title_bar(frame: &mut ratatui::Frame, area: ratatui::layout::Rect, snap: &SharedState) {
    let status_str = match &snap.status {
        ProcessingStatus::Running => "RUNNING",
        ProcessingStatus::Completed => "COMPLETED",
        ProcessingStatus::Error(_) => "ERROR",
    };
    let status_color = match &snap.status {
        ProcessingStatus::Running => Color::Yellow,
        ProcessingStatus::Completed => Color::LightGreen,
        ProcessingStatus::Error(_) => Color::Red,
    };

    let file_progress = if snap.total_files > 1 {
        format!(
            "  [{}/{}]",
            snap.current_file_idx + 1,
            snap.total_files
        )
    } else {
        String::new()
    };

    let title = Line::from(vec![
        Span::styled(
            " BioFastq-A v2.0 ",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(" | ", Style::default().fg(Color::DarkGray)),
        Span::styled(
            status_str,
            Style::default()
                .fg(status_color)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(
            &file_progress,
            Style::default().fg(Color::DarkGray),
        ),
        Span::styled(
            "  |  press Q to exit",
            Style::default().fg(Color::DarkGray),
        ),
    ]);

    let bar = Paragraph::new(title)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .border_style(Style::default().fg(Color::Cyan)),
        )
        .alignment(Alignment::Center);
    frame.render_widget(bar, area);
}

// ---------------------------------------------------------------------------
// Left: progress + quality mini-chart
// ---------------------------------------------------------------------------

fn render_progress_panel(
    frame: &mut ratatui::Frame,
    area: ratatui::layout::Rect,
    snap: &SharedState,
) {
    let block = Block::default()
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Cyan))
        .title(" PROCESSING ")
        .title_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    frame.render_widget(block, area);

    let inner = Layout::default()
        .direction(Direction::Vertical)
        .margin(1)
        .constraints([
            Constraint::Length(1), // status + speed
            Constraint::Length(1), // blank
            Constraint::Length(3), // progress bar
            Constraint::Length(1), // elapsed
            Constraint::Length(1), // blank
            Constraint::Min(4),    // quality mini-chart
        ])
        .split(area);

    let elapsed = snap.start_time.elapsed();
    let elapsed_f = elapsed.as_secs_f64();

    // --- Status & speed ---
    let (reads, bytes_proc, file_size) = if let Some(s) = snap.active_stats() {
        (s.read_count, s.bytes_processed, s.file_size)
    } else {
        (0, 0, 0)
    };

    let rps = if elapsed_f > 0.0 {
        reads as f64 / elapsed_f
    } else {
        0.0
    };

    let speed_line = Line::from(vec![
        Span::styled("Speed: ", Style::default().fg(Color::White)),
        Span::styled(
            format!("{}/s", format_number(rps as u64)),
            Style::default().fg(Color::LightGreen),
        ),
        Span::styled("  Reads: ", Style::default().fg(Color::White)),
        Span::styled(
            format_number(reads),
            Style::default()
                .fg(Color::LightGreen)
                .add_modifier(Modifier::BOLD),
        ),
    ]);
    frame.render_widget(Paragraph::new(speed_line), inner[0]);

    // --- Progress gauge ---
    let progress = if file_size > 0 {
        (bytes_proc as f64 / file_size as f64).min(1.0)
    } else if snap.status == ProcessingStatus::Completed {
        1.0
    } else {
        0.0
    };

    let gauge = Gauge::default()
        .block(
            Block::default()
                .title(" File Progress ")
                .borders(Borders::ALL)
                .border_style(Style::default().fg(Color::DarkGray)),
        )
        .gauge_style(Style::default().fg(Color::Cyan).bg(Color::DarkGray))
        .ratio(progress)
        .label(format!("{:.1}%", progress * 100.0));
    frame.render_widget(gauge, inner[2]);

    // --- Elapsed ---
    let time_line = Line::from(vec![
        Span::styled("Elapsed: ", Style::default().fg(Color::White)),
        Span::styled(
            format_duration(elapsed),
            Style::default().fg(Color::LightGreen),
        ),
    ]);
    frame.render_widget(Paragraph::new(time_line), inner[3]);

    // --- Quality mini bar chart ---
    if let Some(stats) = snap.active_stats() {
        render_quality_sparkline(frame, inner[5], stats);
    }
}

/// Render per-position quality as a Unicode block sparkline
fn render_quality_sparkline(
    frame: &mut ratatui::Frame,
    area: ratatui::layout::Rect,
    stats: &crate::types::FileStats,
) {
    let qual_vec = stats.avg_qual_per_position();
    if qual_vec.is_empty() {
        return;
    }

    let width = area.width.saturating_sub(2) as usize; // inner width
    // Downsample to fit terminal width
    let step = (qual_vec.len() as f64 / width as f64).max(1.0) as usize;
    let sampled: Vec<f64> = (0..width)
        .map(|i| {
            let idx = (i * step).min(qual_vec.len() - 1);
            qual_vec[idx]
        })
        .collect();

    let blocks = [' ', '▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];

    let spans: Vec<Span> = sampled
        .iter()
        .map(|&q| {
            let level = ((q / 40.0) * 8.0).round().min(8.0) as usize;
            let ch = blocks[level];
            let color = if q >= 30.0 {
                Color::LightGreen
            } else if q >= 20.0 {
                Color::Yellow
            } else {
                Color::Red
            };
            Span::styled(ch.to_string(), Style::default().fg(color))
        })
        .collect();

    let header = Line::from(Span::styled(
        "Quality/Position",
        Style::default()
            .fg(Color::DarkGray)
            .add_modifier(Modifier::BOLD),
    ));
    let bar_line = Line::from(spans);

    let sparkline = Paragraph::new(vec![header, bar_line]).block(
        Block::default()
            .borders(Borders::TOP)
            .border_style(Style::default().fg(Color::DarkGray)),
    );
    frame.render_widget(sparkline, area);
}

// ---------------------------------------------------------------------------
// Right: metrics grid
// ---------------------------------------------------------------------------

fn render_metrics_panel(
    frame: &mut ratatui::Frame,
    area: ratatui::layout::Rect,
    snap: &SharedState,
) {
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

    let Some(stats) = snap.active_stats() else {
        return;
    };

    let (n50, n90) = stats.compute_n50_n90();

    let inner = Layout::default()
        .direction(Direction::Vertical)
        .margin(1)
        .constraints([
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Min(0),
        ])
        .split(area);

    let rows: &[(&str, String, Color)] = &[
        (
            "Read Count    ",
            format_number(stats.read_count),
            Color::LightGreen,
        ),
        (
            "Total Bases   ",
            format_bases(stats.total_bases),
            Color::LightGreen,
        ),
        (
            "Avg Length    ",
            format!("{:.1} bp", stats.avg_length()),
            Color::LightGreen,
        ),
        (
            "Min/Max Length",
            format!(
                "{} / {} bp",
                stats.effective_min_length(),
                stats.max_length
            ),
            Color::Cyan,
        ),
        (
            "N50           ",
            format!("{} bp", n50),
            Color::Cyan,
        ),
        (
            "N90           ",
            format!("{} bp", n90),
            Color::Cyan,
        ),
        (
            "GC Content    ",
            format!("{:.2}%", stats.gc_content()),
            Color::LightGreen,
        ),
        (
            "Avg Quality   ",
            format!("Q{:.1}", stats.avg_quality()),
            quality_color(stats.avg_quality()),
        ),
        (
            "Q20 / Q30     ",
            format!("{:.1}% / {:.1}%", stats.q20_pct(), stats.q30_pct()),
            Color::LightGreen,
        ),
        (
            "Adapter Contam",
            format!("{:.2}%", stats.adapter_pct()),
            adapter_color(stats.adapter_pct()),
        ),
    ];

    for (i, (label, value, color)) in rows.iter().enumerate() {
        if i >= inner.len() {
            break;
        }
        let line = Line::from(vec![
            Span::styled(
                format!("  {} ", label),
                Style::default().fg(Color::White),
            ),
            Span::styled(
                value.clone(),
                Style::default()
                    .fg(*color)
                    .add_modifier(Modifier::BOLD),
            ),
        ]);
        frame.render_widget(Paragraph::new(line), inner[i]);
    }
}

fn quality_color(q: f64) -> Color {
    if q >= 30.0 {
        Color::LightGreen
    } else if q >= 20.0 {
        Color::Yellow
    } else {
        Color::Red
    }
}

fn adapter_color(pct: f64) -> Color {
    if pct < 1.0 {
        Color::LightGreen
    } else if pct < 10.0 {
        Color::Yellow
    } else {
        Color::Red
    }
}

// ---------------------------------------------------------------------------
// Bottom: event log
// ---------------------------------------------------------------------------

fn render_log_panel(
    frame: &mut ratatui::Frame,
    area: ratatui::layout::Rect,
    snap: &SharedState,
) {
    let block = Block::default()
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::DarkGray))
        .title(" EVENT LOG ")
        .title_style(
            Style::default()
                .fg(Color::White)
                .add_modifier(Modifier::BOLD),
        );

    let max_visible = area.height.saturating_sub(2) as usize;
    let start = snap
        .log_messages
        .len()
        .saturating_sub(max_visible);

    let items: Vec<ListItem> = snap.log_messages[start..]
        .iter()
        .map(|entry| {
            let ts = format_duration(entry.timestamp);
            let (prefix, color) = match entry.level {
                LogLevel::Info => ("INFO", Color::White),
                LogLevel::Success => (" OK ", Color::LightGreen),
                LogLevel::Warning => ("WARN", Color::Yellow),
                LogLevel::Error => ("ERR ", Color::Red),
            };
            let line = Line::from(vec![
                Span::styled(
                    format!(" [{}] ", ts),
                    Style::default().fg(Color::DarkGray),
                ),
                Span::styled(
                    format!("{} ", prefix),
                    Style::default()
                        .fg(color)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(entry.message.clone(), Style::default().fg(color)),
            ]);
            ListItem::new(line)
        })
        .collect();

    let list = List::new(items).block(block);
    frame.render_widget(list, area);
}
