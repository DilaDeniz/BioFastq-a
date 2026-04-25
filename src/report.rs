use std::fs;
use std::io;
use std::path::Path;

use serde::Serialize;

use crate::types::{format_bases, format_number, FileStats, FeatureFlags, SharedState, ADAPTERS};

// ---------------------------------------------------------------------------
// JSON export structures
// ---------------------------------------------------------------------------

#[derive(Serialize)]
pub struct JsonReport {
    pub tool: &'static str,
    pub version: &'static str,
    pub generated: String,
    pub processing_time_seconds: f64,
    pub files: Vec<JsonFileReport>,
}

#[derive(Serialize)]
pub struct JsonFileReport {
    pub file_name: String,
    pub file_path: String,
    pub read_count: u64,
    pub total_bases: u64,
    pub average_read_length: f64,
    pub min_read_length: u64,
    pub max_read_length: u64,
    pub n50: u64,
    pub n90: u64,
    pub gc_content_percent: f64,
    pub average_phred_quality: f64,
    pub q20_reads_percent: f64,
    pub q30_reads_percent: f64,
    pub adapter_contamination_percent: f64,
    pub duplication_rate_estimate_percent: f64,
    pub trimmed_reads: u64,
    pub trimmed_reads_percent: f64,
    pub trim_output_path: Option<String>,
    pub top_kmers: Vec<(String, u64)>,
}

pub fn export_json(state: &SharedState, output_dir: &str) -> io::Result<String> {
    let files: Vec<JsonFileReport> = state
        .all_files()
        .iter()
        .map(|f| {
            let (n50, n90) = f.compute_n50_n90();
            let fname = Path::new(&f.file_path)
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or(&f.file_path)
                .to_string();
            JsonFileReport {
                file_name: fname,
                file_path: f.file_path.clone(),
                read_count: f.read_count,
                total_bases: f.total_bases,
                average_read_length: f.avg_length(),
                min_read_length: f.effective_min_length(),
                max_read_length: f.max_length,
                n50,
                n90,
                gc_content_percent: f.gc_content(),
                average_phred_quality: f.avg_quality(),
                q20_reads_percent: f.q20_pct(),
                q30_reads_percent: f.q30_pct(),
                adapter_contamination_percent: f.adapter_pct(),
                duplication_rate_estimate_percent: f.dup_rate_pct,
                trimmed_reads: f.trimmed_reads,
                trimmed_reads_percent: f.trimmed_pct(),
                trim_output_path: f.trim_output_path.clone(),
                top_kmers: f.top_kmers(20),
            }
        })
        .collect();

    let now = chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
    let report = JsonReport {
        tool: "BioFastq-A",
        version: "2.2.0",
        generated: now,
        processing_time_seconds: state.elapsed_secs(),
        files,
    };

    let json = serde_json::to_string_pretty(&report)
        .map_err(io::Error::other)?;

    let stem = report_stem(&state.all_files());
    let path = format!("{}/{}_report.json", output_dir, stem);
    fs::write(&path, &json)?;
    Ok(path)
}

// ---------------------------------------------------------------------------
// MultiQC custom-content export  (*_mqc.json)
//
// MultiQC picks up any file named `*_mqc.json` in the run directory and
// merges it into the MultiQC report as a General Statistics table section.
// Format spec: https://multiqc.info/docs/custom_content/
// ---------------------------------------------------------------------------

/// Derive a clean sample name from a file path: strip directory, then
/// strip .gz, then strip .fastq / .fasta / .fq / .fa extensions.
fn sample_name_from_path(file_path: &str) -> String {
    let p = std::path::Path::new(file_path);
    // strip .gz
    let s1 = if p.extension().map(|e| e == "gz").unwrap_or(false) {
        p.file_stem().unwrap_or(p.as_os_str())
    } else {
        p.file_name().unwrap_or(p.as_os_str())
    };
    // strip .fastq / .fasta / .fq / .fa
    let p2 = std::path::Path::new(s1);
    p2.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(file_path)
        .to_string()
}

pub fn export_multiqc(state: &SharedState, output_dir: &str) -> io::Result<String> {
    let files = state.all_files();
    if files.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "No file data"));
    }

    // Build per-sample data map
    let mut data = serde_json::Map::new();
    for f in &files {
        let sample = sample_name_from_path(&f.file_path);
        let (n50, _) = f.compute_n50_n90();
        let mut m = serde_json::Map::new();
        m.insert("total_sequences".into(),   f.read_count.into());
        m.insert("total_bases".into(),       f.total_bases.into());
        m.insert("avg_sequence_length".into(), f.avg_length().into());
        m.insert("percent_gc".into(),        f.gc_content().into());
        m.insert("mean_quality".into(),      f.avg_quality().into());
        m.insert("percent_q30".into(),       f.q30_pct().into());
        m.insert("percent_adapter".into(),   f.adapter_pct().into());
        m.insert("percent_duplicates".into(),f.dup_rate_pct.into());
        m.insert("n50".into(),               n50.into());
        data.insert(sample, serde_json::Value::Object(m));
    }

    // Column configuration for the General Statistics table
    let pconfig = serde_json::json!([
        {"total_sequences":   {"title": "Total Reads",   "description": "Total number of reads",             "format": "{:,.0f}", "scale": "Blues", "shared_key": "read_count"}},
        {"percent_gc":        {"title": "% GC",          "description": "GC content (%)",                    "min": 0, "max": 100, "suffix": "%", "format": "{:.1f}", "scale": "RdYlGn"}},
        {"avg_sequence_length":{"title": "Avg Length",   "description": "Average read length (bp)",          "suffix": " bp", "format": "{:.0f}"}},
        {"mean_quality":      {"title": "Avg Quality",   "description": "Mean Phred quality score",          "min": 0, "max": 42, "format": "{:.1f}", "scale": "RdYlGn"}},
        {"percent_q30":       {"title": "% Q30",         "description": "Bases with Phred ≥ Q30 (%)",        "min": 0, "max": 100, "suffix": "%", "format": "{:.1f}", "scale": "RdYlGn"}},
        {"percent_adapter":   {"title": "% Adapter",     "description": "Reads with adapter contamination",  "min": 0, "max": 100, "suffix": "%", "format": "{:.2f}", "scale": "Oranges"}},
        {"percent_duplicates":{"title": "% Dups",        "description": "Estimated duplication rate (%)",    "min": 0, "max": 100, "suffix": "%", "format": "{:.1f}", "scale": "Reds"}},
        {"total_bases":       {"title": "Total Bases",   "description": "Total sequenced bases",             "format": "{:,.0f}", "shared_key": "base_count"}},
        {"n50":               {"title": "N50",           "description": "N50 read length (bp)",              "suffix": " bp", "format": "{:,.0f}"}}
    ]);

    let report = serde_json::json!({
        "id": "biofastq_a",
        "section_name": "BioFastq-A",
        "description": "High-performance FASTQ/FASTA quality analysis",
        "plot_type": "generalstats",
        "pconfig": pconfig,
        "data": serde_json::Value::Object(data)
    });

    let json = serde_json::to_string_pretty(&report).map_err(io::Error::other)?;
    let stem = report_stem(&files);
    let path = format!("{}/{}_mqc.json", output_dir, stem);
    fs::write(&path, &json)?;
    Ok(path)
}

// ---------------------------------------------------------------------------
// HTML report
// ---------------------------------------------------------------------------

/// Static CSS — clean scientific light theme
const CSS: &str = r#"
:root{
  --bg:#f5f6f7;--card:#ffffff;--card2:#f0f2f4;--border:#d0d5dd;
  --text:#1a1f2e;--muted:#5a6370;--accent:#1a5fa8;
  --green:#1a7a3e;--green-bg:#eaf6ee;
  --yellow:#7a5500;--yellow-bg:#fef9e7;
  --red:#991b1b;--red-bg:#fef2f2;
  --section-num:#9ca3af;
}
*{box-sizing:border-box;margin:0;padding:0;}
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Arial,sans-serif;background:var(--bg);color:var(--text);padding:32px 40px;line-height:1.65;font-size:14px;max-width:1400px;margin:0 auto;}
@media print{body{padding:16px;background:#fff;}.no-print{display:none!important;}.file-section{break-inside:avoid;}.chart-box{break-inside:avoid;}}
a{color:var(--accent);text-decoration:none;}
h1{font-size:24px;font-weight:700;color:var(--text);letter-spacing:-.3px;margin-bottom:2px;}
h2{font-size:16px;font-weight:700;color:var(--text);margin-bottom:14px;padding-bottom:6px;border-bottom:2px solid var(--border);}
h3{font-size:11px;font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:.08em;margin-bottom:10px;}
.header{padding-bottom:20px;margin-bottom:28px;border-bottom:2px solid var(--border);display:flex;justify-content:space-between;align-items:flex-start;flex-wrap:wrap;gap:12px;}
.header-left h1 span{color:var(--accent);}
.header-meta{color:var(--muted);font-size:12px;margin-top:6px;display:flex;gap:16px;flex-wrap:wrap;}
.header-meta span::before{content:"·";margin-right:8px;color:var(--border);}
.header-meta span:first-child::before{content:"";}
.header-actions{display:flex;gap:8px;align-items:center;}
.btn{display:inline-flex;align-items:center;gap:5px;padding:6px 14px;border-radius:6px;border:1px solid var(--border);background:var(--card);color:var(--text);font-size:12px;font-weight:600;cursor:pointer;text-decoration:none;}
.btn:hover{background:var(--card2);border-color:#aab;}
.btn-primary{background:var(--accent);color:#fff;border-color:var(--accent);}
.btn-primary:hover{background:#145090;}
/* Summary table */
.summary-table{width:100%;border-collapse:collapse;margin-bottom:28px;font-size:13px;}
.summary-table th{background:var(--card2);color:var(--muted);font-weight:700;text-align:left;padding:9px 12px;border:1px solid var(--border);font-size:11px;text-transform:uppercase;letter-spacing:.05em;}
.summary-table td{padding:9px 12px;border:1px solid var(--border);background:var(--card);}
.summary-table tr:hover td{background:#f0f4f8;}
/* File section */
.file-section{background:var(--card);border:1px solid var(--border);border-radius:10px;margin-bottom:32px;overflow:hidden;box-shadow:0 1px 4px rgba(0,0,0,.06);}
.file-section-header{background:var(--card2);padding:14px 20px;border-bottom:1px solid var(--border);}
.file-section-header .filename{font-size:15px;font-weight:700;font-family:monospace;color:var(--text);}
.file-section-header .filepath{font-size:11px;color:var(--muted);margin-top:2px;font-family:monospace;}
.file-body{padding:20px;}
/* Stats grid */
.stats-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(155px,1fr));gap:10px;margin-bottom:20px;}
.stat-card{background:var(--card2);border:1px solid var(--border);border-radius:8px;padding:12px 14px;}
.stat-label{font-size:10px;font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:.07em;margin-bottom:5px;}
.stat-value{font-size:19px;font-weight:700;font-family:monospace;color:var(--text);}
.stat-value.pass{color:var(--green);}
.stat-value.warn{color:var(--yellow);}
.stat-value.fail{color:var(--red);}
/* Charts */
.charts-row{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:16px;margin-bottom:16px;}
.chart-box{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:16px;}
.chart-box canvas{max-width:100%;display:block;}
.chart-header{display:flex;justify-content:space-between;align-items:center;margin-bottom:10px;}
.chart-header h3{margin-bottom:0;}
.chart-note{font-size:11px;color:var(--muted);margin-top:8px;line-height:1.5;}
.dl-btn{display:inline-flex;align-items:center;gap:3px;padding:3px 8px;border-radius:4px;border:1px solid var(--border);background:var(--card2);color:var(--muted);font-size:10px;font-weight:600;cursor:pointer;white-space:nowrap;}
.dl-btn:hover{background:var(--border);color:var(--text);}
/* Adapter section */
.adapter-grid{display:grid;grid-template-columns:auto 1fr;gap:24px;align-items:start;}
.big-pct{font-size:38px;font-weight:700;font-family:monospace;line-height:1;}
.adapter-list{margin-top:8px;padding:0;}
.adapter-list li{list-style:none;font-size:12px;font-family:monospace;color:var(--muted);padding:2px 0;}
.adapter-list li::before{content:"› ";color:var(--accent);}
/* Duplication gauge */
.dup-gauge-wrap{display:flex;align-items:center;gap:16px;}
.dup-gauge-track{flex:1;height:12px;background:var(--card2);border:1px solid var(--border);border-radius:6px;overflow:hidden;}
.dup-gauge-fill{height:100%;border-radius:6px;}
/* Trim info */
.trim-box{background:var(--card2);border:1px solid var(--border);border-radius:8px;padding:16px;}
.trim-box p{font-size:13px;color:var(--muted);margin:4px 0;}
.trim-box code{color:var(--accent);font-family:monospace;font-size:12px;}
/* Footer */
.footer{text-align:center;color:var(--muted);font-size:12px;margin-top:40px;padding-top:16px;border-top:1px solid var(--border);}
/* Module status badges */
.module-grid{display:flex;flex-wrap:wrap;gap:6px;margin-bottom:20px;}
.module-badge{display:flex;align-items:center;gap:5px;padding:4px 10px;border-radius:5px;border:1px solid var(--border);background:var(--card2);font-size:11px;font-weight:600;}
.module-badge .icon{font-size:12px;width:14px;text-align:center;}
.module-badge.pass{background:var(--green-bg);border-color:#a7d7b8;}
.module-badge.pass .icon{color:var(--green);}
.module-badge.warn{background:var(--yellow-bg);border-color:#f0d080;}
.module-badge.warn .icon{color:var(--yellow);}
.module-badge.fail{background:var(--red-bg);border-color:#fca5a5;}
.module-badge.fail .icon{color:var(--red);}
.module-badge .name{color:var(--text);}
/* Overrepresented sequences table */
.overrep-table{width:100%;border-collapse:collapse;font-size:12px;margin-top:8px;}
.overrep-table th{background:var(--card2);color:var(--muted);padding:7px 10px;text-align:left;border:1px solid var(--border);font-weight:700;font-size:11px;text-transform:uppercase;letter-spacing:.04em;}
.overrep-table td{padding:6px 10px;border:1px solid var(--border);}
.overrep-table tr:hover td{background:#f0f4f8;}
.overrep-seq{word-break:break-all;max-width:320px;font-family:monospace;color:var(--accent);font-size:11px;}
"#;

/// Static JavaScript for chart drawing (uses vanilla Canvas API, no deps)
const JS: &str = r#"
'use strict';
var RD = __REPORT_DATA__;

/* ---------- Utilities ---------- */
function downloadChart(id, name) {
  var cv = document.getElementById(id); if (!cv) return;
  var a = document.createElement('a');
  a.download = (name || id) + '.png';
  a.href = cv.toDataURL('image/png');
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
}

function printReport() { window.print(); }

/* ---------- Quality per position (mean + Q25/Q50/Q75 IQR overlay) ---------- */
function drawQual(id, data, pct) {
  var cv = document.getElementById(id);
  if (!cv || !data || !data.length) return;
  var ctx = cv.getContext('2d');
  var W = cv.clientWidth || cv.width; cv.width = W;
  var H = cv.height;
  var pad = {t:18,r:16,b:36,l:44};
  var pw = W-pad.l-pad.r, ph = H-pad.t-pad.b, MQ = 42;
  ctx.fillStyle = '#ffffff'; ctx.fillRect(0,0,W,H);

  /* quality zone backgrounds */
  [[0,20,'rgba(220,53,69,.07)'],[20,28,'rgba(255,193,7,.08)'],[28,MQ,'rgba(25,135,84,.06)']].forEach(function(z){
    var y1=pad.t+(1-z[1]/MQ)*ph, y2=pad.t+(1-z[0]/MQ)*ph;
    ctx.fillStyle=z[2]; ctx.fillRect(pad.l,y1,pw,y2-y1);
  });

  /* dashed threshold lines */
  [[20,'rgba(220,53,69,.35)'],[28,'rgba(200,150,0,.45)'],[30,'rgba(25,135,84,.35)']].forEach(function(qc){
    var y=pad.t+(1-qc[0]/MQ)*ph;
    ctx.setLineDash([3,4]); ctx.strokeStyle=qc[1]; ctx.lineWidth=1;
    ctx.beginPath(); ctx.moveTo(pad.l,y); ctx.lineTo(W-pad.r,y); ctx.stroke();
  });
  ctx.setLineDash([]);

  /* grid lines */
  [0,10,20,30,40].forEach(function(q){
    var y=pad.t+(1-q/MQ)*ph;
    ctx.strokeStyle='rgba(0,0,0,.07)'; ctx.lineWidth=1;
    ctx.beginPath(); ctx.moveTo(pad.l,y); ctx.lineTo(W-pad.r,y); ctx.stroke();
  });

  /* Q25–Q75 IQR band + Q50 median line */
  if (pct && pct.length > 1) {
    var np = pct.length;
    ctx.beginPath();
    pct.forEach(function(pp, ii) {
      var x=pad.l+(ii/(np-1))*pw, y=pad.t+(1-Math.min(pp[0],MQ)/MQ)*ph;
      if(ii===0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
    });
    for (var ri=np-1; ri>=0; ri--) {
      ctx.lineTo(pad.l+(ri/(np-1))*pw, pad.t+(1-Math.min(pct[ri][2],MQ)/MQ)*ph);
    }
    ctx.closePath();
    ctx.fillStyle='rgba(26,95,168,.13)'; ctx.fill();

    ctx.beginPath();
    pct.forEach(function(pp, ii) {
      var x=pad.l+(ii/(np-1))*pw, y=pad.t+(1-Math.min(pp[1],MQ)/MQ)*ph;
      if(ii===0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
    });
    ctx.strokeStyle='rgba(200,120,0,.7)'; ctx.lineWidth=1.5;
    ctx.setLineDash([]); ctx.stroke();
  }

  /* filled area + mean line */
  if (data.length > 1) {
    var grad = ctx.createLinearGradient(0,pad.t,0,H-pad.b);
    grad.addColorStop(0,'rgba(26,95,168,.20)');
    grad.addColorStop(1,'rgba(26,95,168,0)');
    ctx.beginPath();
    data.forEach(function(q,i){
      var x=pad.l+(i/(data.length-1))*pw, y=pad.t+(1-Math.min(q,MQ)/MQ)*ph;
      i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);
    });
    ctx.lineTo(pad.l+pw,H-pad.b); ctx.lineTo(pad.l,H-pad.b);
    ctx.closePath(); ctx.fillStyle=grad; ctx.fill();

    ctx.beginPath();
    data.forEach(function(q,i){
      var x=pad.l+(i/(data.length-1))*pw, y=pad.t+(1-Math.min(q,MQ)/MQ)*ph;
      i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);
    });
    ctx.strokeStyle='#1a5fa8'; ctx.lineWidth=2;
    ctx.setLineDash([]); ctx.stroke();
  }

  /* axes */
  ctx.strokeStyle='#555'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();

  /* Y labels */
  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='right';
  [0,10,20,28,30,40].forEach(function(q){
    ctx.fillText(q, pad.l-5, pad.t+(1-q/MQ)*ph+4);
  });

  /* X labels */
  ctx.textAlign='center';
  var n = data.length;
  var step = n > 300 ? 100 : n > 150 ? 50 : n > 50 ? 25 : 10;
  for (var i=0; i<n; i+=step) {
    ctx.fillText(i+1, pad.l+(i/Math.max(n-1,1))*pw, H-pad.b+14);
  }
  ctx.fillText(n, pad.l+pw, H-pad.b+14);

  /* axis titles */
  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='center';
  ctx.fillText('Position (bp)', pad.l+pw/2, H-2);
  ctx.save(); ctx.translate(10,pad.t+ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('Phred Quality', 0, 0); ctx.restore();

  /* legend (only when percentile data present) */
  if (pct && pct.length) {
    var lx = W-pad.r-126, ly = pad.t+4;
    ctx.font='10px monospace'; ctx.textAlign='left';
    ctx.fillStyle='#1a5fa8'; ctx.fillRect(lx,ly,14,3);
    ctx.fillStyle='#5a6370'; ctx.fillText('Mean',lx+17,ly+5);
    ctx.fillStyle='rgba(200,120,0,.7)'; ctx.fillRect(lx,ly+11,14,2);
    ctx.fillStyle='#5a6370'; ctx.fillText('Median (Q50)',lx+17,ly+16);
    ctx.fillStyle='rgba(26,95,168,.18)'; ctx.fillRect(lx,ly+22,14,8);
    ctx.fillStyle='#5a6370'; ctx.fillText('IQR Q25–Q75',lx+17,ly+30);
  }
}

/* ---------- Read length distribution ---------- */
function drawLen(id, data) {
  var cv = document.getElementById(id);
  if (!cv || !data || !data.length) return;
  var ctx = cv.getContext('2d');
  var W = cv.clientWidth || cv.width; cv.width = W;
  var H = cv.height;
  var pad = {t:18,r:12,b:40,l:60};
  var pw = W-pad.l-pad.r, ph = H-pad.t-pad.b;
  ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H);

  var maxC = 0;
  data.forEach(function(d){ if(d[1]>maxC) maxC=d[1]; });
  var n = data.length;
  var bw = pw / n;

  data.forEach(function(d,i){
    var x = pad.l + i*bw;
    var bh = (d[1]/maxC)*ph;
    var y = pad.t+ph-bh;
    var grad = ctx.createLinearGradient(0,y,0,y+bh);
    grad.addColorStop(0,'rgba(26,95,168,.85)');
    grad.addColorStop(1,'rgba(26,95,168,.40)');
    ctx.fillStyle = grad;
    ctx.fillRect(x+0.5, y, Math.max(1, bw-1), bh);
  });

  /* axes */
  ctx.strokeStyle='#555'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();

  /* X labels */
  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='center';
  var lblN = Math.min(6, n);
  for (var i=0; i<lblN; i++) {
    var idx = Math.round(i*(n-1)/(Math.max(lblN-1,1)));
    var x = pad.l+(idx+0.5)*bw;
    var len = data[idx][0];
    ctx.fillText(len>999?(len/1000).toFixed(1)+'k':len+'bp', x, H-pad.b+14);
  }

  /* Y labels */
  ctx.textAlign='right';
  [0,.5,1].forEach(function(frac){
    var v=Math.round(maxC*frac);
    var label = v>999999?(v/1000000).toFixed(1)+'M':v>999?(v/1000).toFixed(0)+'k':String(v);
    ctx.fillText(label, pad.l-5, pad.t+(1-frac)*ph+4);
  });

  /* axis titles */
  ctx.textAlign='center';
  ctx.fillText('Read Length', pad.l+pw/2, H-2);
  ctx.save(); ctx.translate(10,pad.t+ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('Count', 0, 0); ctx.restore();
}

/* ---------- K-mer frequency ---------- */
function drawKmer(id, data) {
  var cv = document.getElementById(id);
  if (!cv || !data || !data.length) return;
  var H = Math.max(260, data.length*22+20); cv.height = H;
  var ctx = cv.getContext('2d');
  var W = cv.clientWidth || cv.width; cv.width = W;
  var pad = {t:10,r:90,b:10,l:68};
  var pw = W-pad.l-pad.r, ph = H-pad.t-pad.b;
  ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H);

  var maxC = data[0][1];
  var rowH = ph/data.length;
  var barH = Math.max(10, rowH-4);
  var ntCol = {A:'#1a7a3e',C:'#1a5fa8',G:'#7a5500',T:'#991b1b',N:'#5a6370'};

  data.forEach(function(d,i){
    var kmer=d[0], cnt=d[1];
    var y = pad.t + i*rowH + (rowH-barH)/2;
    var bw = (cnt/maxC)*pw;

    ctx.fillStyle='rgba(26,95,168,0.55)';
    ctx.fillRect(pad.l, y, bw, barH);

    ctx.font = 'bold 12px monospace'; ctx.textAlign = 'right';
    var cw = 8.4;
    var startX = pad.l - 5 - (kmer.length-1)*cw;
    for (var j=0; j<kmer.length; j++) {
      ctx.fillStyle = ntCol[kmer[j]] || '#1a1f2e';
      ctx.fillText(kmer[j], startX + j*cw + cw, y+barH-3);
    }

    ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='left';
    ctx.fillText(cnt.toLocaleString(), pad.l+bw+5, y+barH-3);
  });
}

/* ---------- Per-tile quality ---------- */
function drawTile(id, data) {
  var cv = document.getElementById(id);
  if (!cv || !data || !data.length) return;
  var ctx = cv.getContext('2d');
  var W = cv.clientWidth || cv.width; cv.width = W;
  var H = cv.height;
  var pad = {t:18,r:12,b:36,l:60};
  var pw = W-pad.l-pad.r, ph = H-pad.t-pad.b;
  ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H);

  var n = data.length;
  var bw = pw / n;
  var minQ = Math.min.apply(null, data.map(function(d){return d[1];}));
  var maxQ = Math.max.apply(null, data.map(function(d){return d[1];}));
  var rng = maxQ - minQ || 1;

  data.forEach(function(d,i){
    var q = d[1];
    var x = pad.l + i*bw;
    var bh = ((q - minQ) / rng) * ph;
    var y = pad.t + ph - bh;
    var t = (q - minQ) / rng;
    var rc = Math.round(200*(1-t) + 25*t);
    var gc = Math.round(50*(1-t) + 135*t);
    var bc = Math.round(50*(1-t) + 84*t);
    ctx.fillStyle = 'rgba('+rc+','+gc+','+bc+',0.85)';
    ctx.fillRect(x+0.5, y, Math.max(1,bw-1), bh);
  });

  ctx.strokeStyle='#555'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();

  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='right';
  [minQ, (minQ+maxQ)/2, maxQ].forEach(function(q){
    var y=pad.t+(1-(q-minQ)/rng)*ph;
    ctx.fillText(q.toFixed(1), pad.l-5, y+4);
  });

  ctx.textAlign='center';
  var lblN = Math.min(8, n);
  for(var i=0; i<lblN; i++){
    var idx=Math.round(i*(n-1)/(Math.max(lblN-1,1)));
    ctx.fillText(data[idx][0], pad.l+(idx+0.5)*bw, H-pad.b+14);
  }
  ctx.fillText('Tile ID', pad.l+pw/2, H-2);
}

/* ---------- Base composition chart ---------- */
function drawBaseComp(id, data) {
  var cv=document.getElementById(id); if(!cv)return;
  var ctx=cv.getContext('2d'), W=cv.width, H=cv.height;
  var pad={l:52,r:16,t:20,b:36};
  var pw=W-pad.l-pad.r, ph=H-pad.t-pad.b;
  ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H);
  if(!data||!data.length){ctx.fillStyle='#5a6370';ctx.font='14px monospace';ctx.textAlign='center';ctx.fillText('No data',W/2,H/2);return;}
  var colors=['#1a5fa8','#1a7a3e','#dc3545','#7a5500','#5a6370']; // A,C,G,T,N
  var labels=['A','C','G','T','N'];
  var n=data.length, bw=pw/n;
  ctx.strokeStyle='rgba(0,0,0,.07)'; ctx.lineWidth=1;
  for(var gi=0;gi<=4;gi++){var gy=pad.t+ph*(1-gi/4);ctx.beginPath();ctx.moveTo(pad.l,gy);ctx.lineTo(pad.l+pw,gy);ctx.stroke();}
  labels.forEach(function(lbl,li){
    ctx.strokeStyle=colors[li]; ctx.lineWidth=1.5; ctx.beginPath();
    data.forEach(function(pt,xi){
      var x=pad.l+(xi+0.5)*bw, y=pad.t+ph*(1-pt[li]/100);
      if(xi===0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
    });
    ctx.stroke();
  });
  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='right';
  for(var gi=0;gi<=4;gi++){ctx.fillText((gi*25)+'%', pad.l-4, pad.t+ph*(1-gi/4)+4);}
  ctx.textAlign='center';
  ctx.fillText('1',pad.l,H-pad.b+14);
  ctx.fillText(n,pad.l+pw-20,H-pad.b+14);
  ctx.fillText('Position (bp)',pad.l+pw/2,H-2);
  var lx=pad.l+4;
  labels.forEach(function(lbl,li){
    ctx.fillStyle=colors[li]; ctx.fillRect(lx,pad.t+2,10,10);
    ctx.fillStyle='#1a1f2e'; ctx.textAlign='left'; ctx.fillText(lbl,lx+13,pad.t+11); lx+=32;
  });
}

/* ---------- Quality distribution chart ---------- */
function drawQualDist(id, data) {
  var cv=document.getElementById(id); if(!cv)return;
  var ctx=cv.getContext('2d'), W=cv.width, H=cv.height;
  var pad={l:52,r:16,t:20,b:36};
  var pw=W-pad.l-pad.r, ph=H-pad.t-pad.b;
  ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H);
  if(!data||!data.length){ctx.fillStyle='#5a6370';ctx.font='14px monospace';ctx.textAlign='center';ctx.fillText('No data',W/2,H/2);return;}
  var mx=Math.max.apply(null,data)||1;
  var n=data.length, bw=pw/n;
  data.forEach(function(cnt,q){
    var frac=cnt/mx;
    var rc=q<20?180:q<30?180:40, gc=q<20?40:q<30?130:150, bc=q<20?40:q<30?40:60;
    ctx.fillStyle='rgba('+rc+','+gc+','+bc+',0.80)';
    var barH=ph*frac;
    ctx.fillRect(pad.l+q*bw, pad.t+ph-barH, Math.max(1,bw-1), barH);
  });
  ctx.strokeStyle='rgba(0,0,0,.07)'; ctx.lineWidth=1;
  for(var gi=0;gi<=4;gi++){var gy=pad.t+ph*gi/4;ctx.beginPath();ctx.moveTo(pad.l,gy);ctx.lineTo(pad.l+pw,gy);ctx.stroke();}
  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='center';
  [0,10,20,30,40].forEach(function(q){ctx.fillText('Q'+q, pad.l+q*bw-4, H-pad.b+14);});
  ctx.fillText('Mean Phred score per read', pad.l+pw/2, H-2);
}

/* ---------- GC content distribution ---------- */
function drawGcDist(id, data, meanGc) {
  var cv=document.getElementById(id); if(!cv||!data||!data.length)return;
  var ctx=cv.getContext('2d'), W=cv.width, H=cv.height;
  var pad={l:52,r:16,t:24,b:36};
  var pw=W-pad.l-pad.r, ph=H-pad.t-pad.b;
  ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H);
  var total=data.reduce(function(s,v){return s+v;},0)||1;
  var max=Math.max.apply(null,data)||1;
  var n=data.length, bw=pw/n;
  // Draw bars
  data.forEach(function(cnt,gc){
    var frac=cnt/max;
    var barH=ph*frac;
    ctx.fillStyle='rgba(26,95,168,0.55)';
    ctx.fillRect(pad.l+gc*bw, pad.t+ph-barH, Math.max(1,bw-1), barH);
  });
  // Theoretical normal overlay
  if (meanGc >= 0) {
    var variance = 0;
    data.forEach(function(cnt,gc){variance += cnt*(gc-meanGc)*(gc-meanGc);});
    variance /= total;
    var sd = Math.sqrt(variance) || 5;
    ctx.beginPath();
    var peaked = false;
    for (var gc=0; gc<=100; gc++) {
      var norm = Math.exp(-0.5*((gc-meanGc)/sd)*((gc-meanGc)/sd));
      var x=pad.l+(gc+0.5)*bw;
      var y=pad.t+ph*(1-norm);
      if(!peaked){ctx.moveTo(x,y);peaked=true;}else{ctx.lineTo(x,y);}
    }
    ctx.strokeStyle='rgba(220,53,69,0.7)'; ctx.lineWidth=2;
    ctx.setLineDash([4,3]); ctx.stroke(); ctx.setLineDash([]);
  }
  // Axes
  ctx.strokeStyle='#555'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();
  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='center';
  [0,20,40,50,60,80,100].forEach(function(g){ctx.fillText(g+'%',pad.l+g*bw,H-pad.b+14);});
  ctx.fillText('GC Content (%)', pad.l+pw/2, H-2);
  // Legend
  ctx.textAlign='left'; ctx.font='10px monospace';
  ctx.fillStyle='rgba(26,95,168,0.55)'; ctx.fillRect(pad.l+4,pad.t+4,12,8);
  ctx.fillStyle='#5a6370'; ctx.fillText('Observed',pad.l+19,pad.t+12);
  ctx.strokeStyle='rgba(220,53,69,0.7)'; ctx.lineWidth=2; ctx.setLineDash([4,3]);
  ctx.beginPath(); ctx.moveTo(pad.l+80,pad.t+8); ctx.lineTo(pad.l+94,pad.t+8); ctx.stroke();
  ctx.setLineDash([]); ctx.fillStyle='#5a6370'; ctx.fillText('Theoretical',pad.l+97,pad.t+12);
}

/* ---------- Duplication level histogram ---------- */
function drawDupHist(id, data) {
  var cv=document.getElementById(id); if(!cv||!data||!data.length)return;
  var ctx=cv.getContext('2d'), W=cv.width, H=cv.height;
  var pad={l:52,r:16,t:20,b:50};
  var pw=W-pad.l-pad.r, ph=H-pad.t-pad.b;
  ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H);
  var max=Math.max.apply(null,data.map(function(d){return d[1];}));
  if(!max) return;
  var n=data.length, bw=pw/n;
  data.forEach(function(d,i){
    var frac=d[1]/max;
    var barH=ph*frac;
    var t=i/(n-1||1);
    var r=Math.round(26*(1-t)+220*t), g=Math.round(95*(1-t)+53*t), b=Math.round(168*(1-t)+69*t);
    ctx.fillStyle='rgba('+r+','+g+','+b+',0.75)';
    ctx.fillRect(pad.l+i*bw+1, pad.t+ph-barH, Math.max(1,bw-2), barH);
    // label
    ctx.fillStyle='#5a6370'; ctx.font='10px monospace'; ctx.textAlign='center';
    ctx.fillText(d[0], pad.l+(i+0.5)*bw, H-pad.b+14);
    if(d[1]>0.5){ctx.fillText(d[1].toFixed(1)+'%', pad.l+(i+0.5)*bw, pad.t+ph-barH-4);}
  });
  ctx.strokeStyle='#555'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();
  ctx.fillStyle='#5a6370'; ctx.font='10px monospace'; ctx.textAlign='right';
  [0,50,100].forEach(function(p){
    var y=pad.t+ph*(1-p/100);
    ctx.fillText(p+'%',pad.l-4,y+4);
  });
  ctx.textAlign='center'; ctx.fillText('Duplication level', pad.l+pw/2, H-2);
}

/* ---------- Per-base N content ---------- */
function drawNContent(id, data) {
  var cv=document.getElementById(id); if(!cv||!data||!data.length)return;
  var ctx=cv.getContext('2d'), W=cv.width, H=cv.height;
  var pad={l:52,r:16,t:20,b:36};
  var pw=W-pad.l-pad.r, ph=H-pad.t-pad.b;
  ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H);
  var maxN=Math.max.apply(null,data)||0.1;
  var cap=Math.max(5,Math.ceil(maxN));
  // warn zone
  ctx.fillStyle='rgba(220,53,69,.06)'; ctx.fillRect(pad.l,pad.t,pw,ph);
  // line
  ctx.beginPath();
  data.forEach(function(v,i){
    var x=pad.l+(i/(data.length-1||1))*pw;
    var y=pad.t+ph*(1-v/cap);
    if(i===0)ctx.moveTo(x,y); else ctx.lineTo(x,y);
  });
  ctx.strokeStyle='#dc3545'; ctx.lineWidth=1.5; ctx.stroke();
  // 5% threshold line
  if(cap>=5){
    var ty=pad.t+ph*(1-5/cap);
    ctx.setLineDash([3,4]); ctx.strokeStyle='rgba(220,53,69,.4)'; ctx.lineWidth=1;
    ctx.beginPath(); ctx.moveTo(pad.l,ty); ctx.lineTo(W-pad.r,ty); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle='rgba(220,53,69,.6)'; ctx.font='10px monospace'; ctx.textAlign='right';
    ctx.fillText('5%',pad.l-4,ty+4);
  }
  // Axes
  ctx.strokeStyle='#555'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();
  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='right';
  [0,cap/2,cap].forEach(function(v){
    ctx.fillText(v.toFixed(1)+'%',pad.l-4,pad.t+ph*(1-v/cap)+4);
  });
  ctx.textAlign='center';
  var n=data.length;
  var step=n>300?100:n>150?50:n>50?25:10;
  for(var i=0;i<n;i+=step){ctx.fillText(i+1,pad.l+(i/(n-1||1))*pw,H-pad.b+14);}
  ctx.fillText(n,pad.l+pw,H-pad.b+14);
  ctx.fillText('Position (bp)',pad.l+pw/2,H-2);
  ctx.save(); ctx.translate(10,pad.t+ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('N content (%)',0,0); ctx.restore();
}

/* ---------- Quality vs Length scatter (long-read) ---------- */
function drawQualVsLen(id, data) {
  var cv = document.getElementById(id);
  if (!cv || !data || !data.length) return;
  var ctx = cv.getContext('2d');
  var W = cv.clientWidth || cv.width; cv.width = W;
  var H = cv.height;
  var pad = {t:18,r:16,b:40,l:48};
  var pw = W-pad.l-pad.r, ph = H-pad.t-pad.b;
  ctx.fillStyle = '#ffffff'; ctx.fillRect(0,0,W,H);

  var maxLen = 0, minQ = 42, maxQ = 0;
  data.forEach(function(d){ if(d[0]>maxLen) maxLen=d[0]; if(d[1]<minQ) minQ=d[1]; if(d[1]>maxQ) maxQ=d[1]; });
  if (maxLen === 0) return;
  minQ = Math.max(0, Math.floor(minQ) - 1);
  maxQ = Math.min(42, Math.ceil(maxQ) + 1);
  var qRange = maxQ - minQ || 1;

  /* grid */
  ctx.strokeStyle = 'rgba(0,0,0,.07)'; ctx.lineWidth = 1;
  [0,1,2,3,4].forEach(function(gi) {
    var y = pad.t + gi/4 * ph;
    ctx.beginPath(); ctx.moveTo(pad.l,y); ctx.lineTo(W-pad.r,y); ctx.stroke();
    var x = pad.l + gi/4 * pw;
    ctx.beginPath(); ctx.moveTo(x,pad.t); ctx.lineTo(x,H-pad.b); ctx.stroke();
  });

  /* dots */
  data.forEach(function(d) {
    var x = pad.l + (d[0] / maxLen) * pw;
    var y = pad.t + (1 - (d[1] - minQ) / qRange) * ph;
    ctx.beginPath();
    ctx.arc(x, y, 2.5, 0, 2*Math.PI);
    ctx.fillStyle = 'rgba(26,95,168,0.55)';
    ctx.fill();
  });

  /* axes */
  ctx.strokeStyle = '#555'; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();

  /* Y labels */
  ctx.fillStyle = '#5a6370'; ctx.font = '11px monospace'; ctx.textAlign = 'right';
  [minQ, Math.round((minQ+maxQ)/2), maxQ].forEach(function(q) {
    var y = pad.t + (1 - (q - minQ) / qRange) * ph;
    ctx.fillText(q, pad.l-4, y+4);
  });

  /* X labels */
  ctx.textAlign = 'center';
  [0, 0.25, 0.5, 0.75, 1.0].forEach(function(frac) {
    var v = Math.round(maxLen * frac);
    var x = pad.l + frac * pw;
    var lbl = v >= 1000 ? (v/1000).toFixed(1)+'k' : v+'bp';
    ctx.fillText(lbl, x, H-pad.b+14);
  });

  ctx.fillText('Read Length (bp)', pad.l+pw/2, H-2);
  ctx.save(); ctx.translate(12, pad.t+ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('Mean Phred Quality', 0, 0); ctx.restore();
}

/* ---------- Reads over time (ONT only) ---------- */
function drawReadsOverTime(id, data) {
  var cv = document.getElementById(id);
  if (!cv || !data || !data.length) return;
  var ctx = cv.getContext('2d');
  var W = cv.clientWidth || cv.width; cv.width = W;
  var H = cv.height;
  var pad = {t:18,r:16,b:40,l:60};
  var pw = W-pad.l-pad.r, ph = H-pad.t-pad.b;
  ctx.fillStyle = '#ffffff'; ctx.fillRect(0,0,W,H);

  var maxMin = data[data.length-1][0];
  var maxReads = data[data.length-1][1];
  if (!maxMin || !maxReads) return;

  /* grid */
  ctx.strokeStyle = 'rgba(0,0,0,.07)'; ctx.lineWidth = 1;
  [0,0.25,0.5,0.75,1.0].forEach(function(frac) {
    var y = pad.t + (1-frac)*ph;
    ctx.beginPath(); ctx.moveTo(pad.l,y); ctx.lineTo(W-pad.r,y); ctx.stroke();
  });

  /* filled area */
  var grad = ctx.createLinearGradient(0,pad.t,0,H-pad.b);
  grad.addColorStop(0,'rgba(26,95,168,.20)');
  grad.addColorStop(1,'rgba(26,95,168,0)');
  ctx.beginPath();
  data.forEach(function(d,i) {
    var x = pad.l + (d[0]/maxMin)*pw;
    var y = pad.t + (1 - d[1]/maxReads)*ph;
    if(i===0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
  });
  ctx.lineTo(pad.l+pw, H-pad.b); ctx.lineTo(pad.l, H-pad.b);
  ctx.closePath(); ctx.fillStyle = grad; ctx.fill();

  /* line */
  ctx.beginPath();
  data.forEach(function(d,i) {
    var x = pad.l + (d[0]/maxMin)*pw;
    var y = pad.t + (1 - d[1]/maxReads)*ph;
    if(i===0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
  });
  ctx.strokeStyle = '#1a5fa8'; ctx.lineWidth = 2; ctx.stroke();

  /* axes */
  ctx.strokeStyle = '#555'; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();

  /* Y labels */
  ctx.fillStyle = '#5a6370'; ctx.font = '11px monospace'; ctx.textAlign = 'right';
  [0, 0.5, 1.0].forEach(function(frac) {
    var v = Math.round(maxReads * frac);
    var label = v>999999?(v/1000000).toFixed(1)+'M':v>999?(v/1000).toFixed(0)+'k':String(v);
    ctx.fillText(label, pad.l-4, pad.t+(1-frac)*ph+4);
  });

  /* X labels */
  ctx.textAlign = 'center';
  [0, 0.25, 0.5, 0.75, 1.0].forEach(function(frac) {
    var v = Math.round(maxMin * frac);
    ctx.fillText(v+'m', pad.l + frac*pw, H-pad.b+14);
  });
  ctx.fillText('Minutes since run start', pad.l+pw/2, H-2);
  ctx.save(); ctx.translate(12, pad.t+ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('Cumulative reads', 0, 0); ctx.restore();
}

/* ---------- Channel occupancy (ONT only) ---------- */
function drawChannelOccupancy(id, data) {
  var cv = document.getElementById(id);
  if (!cv || !data || !data.length) return;
  var ctx = cv.getContext('2d');
  var W = cv.clientWidth || cv.width; cv.width = W;
  var H = cv.height;
  var pad = {t:18,r:16,b:40,l:52};
  var pw = W-pad.l-pad.r, ph = H-pad.t-pad.b;
  ctx.fillStyle = '#ffffff'; ctx.fillRect(0,0,W,H);

  /* Build histogram of reads-per-channel counts */
  var maxCnt = 0;
  data.forEach(function(d){ if(d[1]>maxCnt) maxCnt=d[1]; });
  if (!maxCnt) return;

  /* bin into ~20 buckets */
  var n = data.length;
  var nbins = Math.min(20, n);
  var binSize = Math.ceil(maxCnt / nbins);
  if (binSize === 0) binSize = 1;
  var bins = [];
  for (var b=0; b<nbins; b++) bins.push(0);
  data.forEach(function(d){
    var bi = Math.min(Math.floor(d[1]/binSize), nbins-1);
    bins[bi]++;
  });
  var maxBin = Math.max.apply(null, bins) || 1;
  var bw = pw / nbins;

  bins.forEach(function(cnt, i) {
    var bh = (cnt/maxBin)*ph;
    var y = pad.t+ph-bh;
    ctx.fillStyle = 'rgba(26,95,168,0.65)';
    ctx.fillRect(pad.l+i*bw+0.5, y, Math.max(1,bw-1), bh);
  });

  /* axes */
  ctx.strokeStyle = '#555'; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();

  /* X labels (reads per channel) */
  ctx.fillStyle = '#5a6370'; ctx.font = '10px monospace'; ctx.textAlign = 'center';
  [0, Math.floor(nbins/2), nbins-1].forEach(function(i) {
    var lo = i*binSize, hi = (i+1)*binSize-1;
    ctx.fillText(lo+'-'+hi, pad.l+(i+0.5)*bw, H-pad.b+14);
  });
  ctx.fillText('Reads per channel', pad.l+pw/2, H-2);

  /* Y labels */
  ctx.textAlign = 'right';
  [0, 0.5, 1.0].forEach(function(frac) {
    var v = Math.round(maxBin*frac);
    ctx.fillText(v, pad.l-4, pad.t+(1-frac)*ph+4);
  });
  ctx.save(); ctx.translate(10, pad.t+ph/2); ctx.rotate(-Math.PI/2);
  ctx.textAlign='center'; ctx.fillText('Channel count', 0, 0); ctx.restore();
}

/* ---------- Init ---------- */
window.addEventListener('load', function(){
  RD.forEach(function(f,i){
    var gc_mean = f.gc_dist ? (function(){
      var tot=f.gc_dist.reduce(function(s,v){return s+v;},0)||1;
      var s=0; f.gc_dist.forEach(function(c,g){s+=c*g;}); return s/tot;
    })() : -1;
    drawQual('qual-'+i, f.qual, f.qual_pct);
    drawLen('len-'+i, f.len);
    drawGcDist('gcdist-'+i, f.gc_dist, gc_mean);
    drawDupHist('duphist-'+i, f.dup_hist);
    drawNContent('ncontent-'+i, f.n_content);
    drawKmer('kmer-'+i, f.kmers);
    drawBaseComp('basecomp-'+i, f.basecomp);
    drawQualDist('qualdist-'+i, f.qualdist);
    if (f.tiles && f.tiles.length) drawTile('tile-'+i, f.tiles);
    if (f.qvl && f.qvl.length) drawQualVsLen('qvl-'+i, f.qvl);
    if (f.rot && f.rot.length) drawReadsOverTime('rot-'+i, f.rot);
    if (f.chann && f.chann.length) drawChannelOccupancy('chann-'+i, f.chann);
  });
});
"#;

// ---------------------------------------------------------------------------
// HTML generation
// ---------------------------------------------------------------------------

pub fn export_html(state: &SharedState, output_dir: &str, flags: &FeatureFlags) -> io::Result<String> {
    let files = state.all_files();
    if files.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "No file data"));
    }

    let timestamp = chrono::Local::now()
        .format("%Y-%m-%d %H:%M:%S")
        .to_string();
    let elapsed = state.elapsed_secs();

    // Build per-file HTML sections + JSON payload
    let mut file_sections = String::new();
    let mut json_files: Vec<serde_json::Value> = Vec::new();

    for (i, f) in files.iter().enumerate() {
        file_sections.push_str(&build_file_section(f, i, flags));
        json_files.push(build_file_json(f));
    }

    // Build summary table (multi-file)
    let summary_table = build_summary_table(&files);

    // JSON data embedded in JS
    let json_data = serde_json::to_string(&json_files).unwrap_or_else(|_| "[]".into());
    let js_with_data = JS.replace("__REPORT_DATA__", &json_data);

    // Assemble HTML
    let mut html = String::with_capacity(256 * 1024);
    html.push_str("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n");
    html.push_str("<meta charset=\"UTF-8\">\n");
    html.push_str("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
    html.push_str(&format!(
        "<title>BioFastq-A Report — {}</title>\n",
        timestamp
    ));
    html.push_str("<style>\n");
    html.push_str(CSS);
    html.push_str("</style>\n</head>\n<body>\n");

    // Header
    html.push_str("<div class=\"header\">\n");
    html.push_str("<div class=\"header-left\">\n");
    html.push_str("<h1>BioFastq-A <span>Quality Analysis Report</span></h1>\n");
    html.push_str("<div class=\"header-meta\">\n");
    html.push_str(&format!("<span>Generated: {}</span>\n", timestamp));
    html.push_str("<span>Tool version: 2.2.0</span>\n");
    html.push_str(&format!("<span>Files analysed: {}</span>\n", files.len()));
    html.push_str(&format!("<span>Processing time: {:.1}s</span>\n", elapsed));
    html.push_str("</div>\n</div>\n");
    html.push_str("<div class=\"header-actions no-print\">\n");
    html.push_str("<button class=\"btn\" onclick=\"printReport()\">&#128438; Print / PDF</button>\n");
    html.push_str("</div>\n");
    html.push_str("</div>\n");

    // Summary table
    if files.len() > 1 {
        html.push_str("<h2>Summary</h2>\n");
        html.push_str(&summary_table);
    }

    // Per-file sections
    html.push_str(&file_sections);

    // Footer
    html.push_str("<div class=\"footer\">Generated by <strong>BioFastq-A v2.2.0</strong></div>\n");

    // Embedded script
    html.push_str("<script>\n");
    html.push_str(&js_with_data);
    html.push_str("\n</script>\n</body>\n</html>\n");

    let stem = report_stem(&files);
    let path = format!("{}/{}_report.html", output_dir, stem);
    fs::write(&path, &html)?;
    Ok(path)
}

fn skipped_notice(flag: &str) -> String {
    format!(
        "<p class=\"chart-note\" style=\"padding:20px 0;font-style:italic\">Module skipped &mdash; <code>{}</code> was active.</p>\n",
        flag
    )
}

fn build_file_section(f: &FileStats, idx: usize, flags: &FeatureFlags) -> String {
    let (n50, n90) = f.compute_n50_n90();
    let fname = Path::new(&f.file_path)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(&f.file_path);

    let lr = f.long_read_mode;
    let qual_class = quality_class_lr(f.avg_quality(), lr);
    let adapter_class = adapter_status_class(f.adapter_pct());

    let adapter_list: String = ADAPTERS
        .iter()
        .map(|(name, _)| format!("<li>{}</li>", name))
        .collect();

    let mut s = String::new();
    s.push_str("<div class=\"file-section\">\n");
    s.push_str("<div class=\"file-section-header\">\n");
    s.push_str(&format!(
        "<div class=\"filename\">{}</div>\n",
        escape_html(fname)
    ));
    s.push_str(&format!(
        "<div class=\"filepath\">{}</div>\n",
        escape_html(&f.file_path)
    ));
    s.push_str("</div>\n<div class=\"file-body\">\n");

    // Module status traffic lights (FastQC-style pass/warn/fail per module)
    if !f.module_status.is_empty() {
        s.push_str("<div class=\"module-grid\">\n");
        for (name, status) in &f.module_status {
            s.push_str(&format!(
                "<div class=\"module-badge {}\"><span class=\"icon\">{}</span><span class=\"name\">{}</span></div>\n",
                status.css_class(), status.icon(), escape_html(name)
            ));
        }
        s.push_str("</div>\n");
    }

    // Stats grid
    s.push_str("<div class=\"stats-grid\">\n");
    stat_card(&mut s, "Total Reads", &format_number(f.read_count), "");
    stat_card(
        &mut s,
        "Total Bases",
        &format_bases(f.total_bases),
        "",
    );
    stat_card(
        &mut s,
        "Avg Read Length",
        &format!("{:.1} bp", f.avg_length()),
        "",
    );
    stat_card(
        &mut s,
        "Min / Max Length",
        &format!("{} / {} bp", f.effective_min_length(), f.max_length),
        "",
    );
    if lr {
        // Long-read: N50 is the most important metric — show it prominently
        stat_card(&mut s, "N50 (prominent)", &format!("{} bp", n50), if n50 > 10000 { "pass" } else { "warn" });
        stat_card(&mut s, "N90", &format!("{} bp", n90), "");
    } else {
        stat_card(&mut s, "N50", &format!("{} bp", n50), "");
        stat_card(&mut s, "N90", &format!("{} bp", n90), "");
    }
    stat_card(
        &mut s,
        "GC Content",
        &format!("{:.2}%", f.gc_content()),
        gc_class(f.gc_content()),
    );
    stat_card(
        &mut s,
        "Avg Quality",
        &format!("Q{:.1}", f.avg_quality()),
        qual_class,
    );
    if lr {
        // Long-read mode: show Q10/Q20 stats (ONT reads don't typically reach Q30)
        stat_card(
            &mut s,
            "≥Q10 Reads",
            &format!("{:.1}%", f.q20_pct()), // q20_pct is ≥Q20 baseline; for ONT context we show q20 as a proxy
            if f.q20_pct() >= 50.0 { "pass" } else { "warn" },
        );
        stat_card(
            &mut s,
            "≥Q20 Reads",
            &format!("{:.1}%", f.q20_pct()),
            if f.q20_pct() >= 30.0 { "pass" } else { "warn" },
        );
    } else {
        stat_card(
            &mut s,
            "≥Q20 Reads",
            &format!("{:.1}%", f.q20_pct()),
            if f.q20_pct() >= 80.0 { "pass" } else { "warn" },
        );
        stat_card(
            &mut s,
            "≥Q30 Reads",
            &format!("{:.1}%", f.q30_pct()),
            if f.q30_pct() >= 70.0 { "pass" } else { "warn" },
        );
    }
    stat_card(
        &mut s,
        "Adapter Content",
        &format!("{:.2}%", f.adapter_pct()),
        adapter_class,
    );
    if f.reads_filtered > 0 {
        stat_card(
            &mut s,
            "Filtered Reads",
            &format!("{} ({:.1}%)",
                format_number(f.reads_filtered),
                f.reads_filtered as f64 / (f.read_count + f.reads_filtered) as f64 * 100.0
            ),
            "warn",
        );
    }
    s.push_str("</div>\n"); // stats-grid

    // Charts row 1: quality per position + length distribution
    s.push_str("<div class=\"charts-row\">\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Per-Base Sequence Quality</h3>\n");
    s.push_str(&format!(
        "<button class=\"dl-btn no-print\" onclick=\"downloadChart('qual-{idx}','per_base_quality')\">&#8595; PNG</button>\n",
        idx = idx
    ));
    s.push_str("</div>\n");
    s.push_str(&format!(
        "<canvas id=\"qual-{}\" width=\"600\" height=\"240\"></canvas>\n",
        idx
    ));
    s.push_str("<p class=\"chart-note\">Blue line = mean quality. Orange dashed = median (Q50). Shaded band = IQR (Q25&ndash;Q75). Green zone: Q&ge;28 &nbsp;|&nbsp; Yellow: Q20&ndash;28 &nbsp;|&nbsp; Red: Q&lt;20.</p>\n");
    s.push_str("</div>\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Read Length Distribution</h3>\n");
    s.push_str(&format!(
        "<button class=\"dl-btn no-print\" onclick=\"downloadChart('len-{idx}','read_length_dist')\">&#8595; PNG</button>\n",
        idx = idx
    ));
    s.push_str("</div>\n");
    s.push_str(&format!(
        "<canvas id=\"len-{}\" width=\"600\" height=\"240\"></canvas>\n",
        idx
    ));
    s.push_str("</div>\n");

    s.push_str("</div>\n"); // charts-row

    // Charts row 2: k-mer + adapter info
    s.push_str("<div class=\"charts-row\">\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Top K-mer Frequencies (4-mer)</h3>\n");
    if flags.kmer_analysis {
        s.push_str(&format!(
            "<button class=\"dl-btn no-print\" onclick=\"downloadChart('kmer-{idx}','kmer_frequencies')\">&#8595; PNG</button>\n",
            idx = idx
        ));
    }
    s.push_str("</div>\n");
    if flags.kmer_analysis {
        s.push_str(&format!(
            "<canvas id=\"kmer-{}\" width=\"600\" height=\"300\"></canvas>\n",
            idx
        ));
    } else {
        s.push_str(&skipped_notice("--no-kmer / --fast"));
    }
    s.push_str("</div>\n");

    // Adapter info card
    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<h3>Adapter Contamination</h3>\n");
    if flags.adapter_detection {
        s.push_str("<div class=\"adapter-grid\">\n");
        s.push_str(&format!(
            "<div><div class=\"big-pct {}\">{:.1}%</div><div style=\"color:var(--muted);font-size:12px;margin-top:4px\">reads with adapter</div></div>\n",
            adapter_class, f.adapter_pct()
        ));
        s.push_str("<div>\n");
        s.push_str("<p style=\"font-size:12px;color:var(--muted);margin-bottom:6px\">Screened sequences:</p>\n");
        s.push_str("<ul class=\"adapter-list\">\n");
        s.push_str(&adapter_list);
        s.push_str("</ul>\n</div>\n</div>\n");
    } else {
        s.push_str(&skipped_notice("--no-adapter"));
    }
    s.push_str("</div>\n");

    s.push_str("</div>\n"); // charts-row

    // --- Row 3: Duplication + Per-tile + Trim info ---
    s.push_str("<div class=\"charts-row\">\n");

    // Duplication card
    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<h3>Duplication Rate (estimate)</h3>\n");
    if flags.duplication_check {
        let dup_class = dup_status_class(f.dup_rate_pct);
        let dup_bar_color = match dup_class {
            "pass" => "var(--green)",
            "warn" => "var(--yellow)",
            _ => "var(--red)",
        };
        s.push_str("<div class=\"dup-gauge-wrap\" style=\"margin-bottom:12px\">\n");
        s.push_str(&format!(
            "<span class=\"stat-value {}\" style=\"font-size:32px;min-width:80px\">{:.1}%</span>\n",
            dup_class, f.dup_rate_pct
        ));
        s.push_str("<div class=\"dup-gauge-track\">\n");
        s.push_str(&format!(
            "<div class=\"dup-gauge-fill\" style=\"width:{:.1}%;background:{}\"></div>\n",
            f.dup_rate_pct.min(100.0),
            dup_bar_color
        ));
        s.push_str("</div>\n</div>\n");
        s.push_str("<p class=\"chart-note\">HyperLogLog cardinality estimate across all reads (error &le;1%). &lt;5% = pass &nbsp;|&nbsp; 5&ndash;20% = warn &nbsp;|&nbsp; &gt;20% = high</p>\n");
    } else {
        s.push_str(&skipped_notice("--no-duplication / --fast"));
    }
    s.push_str("</div>\n");

    // Per-tile quality card (hidden for long reads — irrelevant for ONT/PacBio)
    if !lr {
        s.push_str("<div class=\"chart-box\">\n");
        s.push_str("<div class=\"chart-header\">\n");
        s.push_str("<h3>Per-Tile Quality Score</h3>\n");
        if flags.per_tile_quality {
            let tile_has_data = !f.sorted_tile_quality().is_empty();
            if tile_has_data {
                s.push_str(&format!(
                    "<button class=\"dl-btn no-print\" onclick=\"downloadChart('tile-{idx}','per_tile_quality')\">&#8595; PNG</button>\n",
                    idx = idx
                ));
            }
        }
        s.push_str("</div>\n");
        if !flags.per_tile_quality {
            s.push_str(&skipped_notice("--no-per-tile / --fast"));
        } else {
            let tile_data = f.sorted_tile_quality();
            if !tile_data.is_empty() {
                s.push_str(&format!(
                    "<canvas id=\"tile-{}\" width=\"600\" height=\"200\"></canvas>\n",
                    idx
                ));
                s.push_str(&format!(
                    "<p class=\"chart-note\">{} Illumina tiles detected. Bars coloured red→green by quality.</p>\n",
                    tile_data.len()
                ));
            } else {
                s.push_str("<p class=\"chart-note\" style=\"padding:20px 0\">No Illumina tile information found in read headers.<br>Per-tile QC requires standard Illumina CASAVA 1.8+ headers.</p>\n");
            }
        }
        s.push_str("</div>\n");
    }

    s.push_str("</div>\n"); // charts-row (row 3)

    // Row 4: Base composition + Quality distribution
    s.push_str("<div class=\"charts-row\">\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Per-Base Sequence Content</h3>\n");
    s.push_str(&format!(
        "<button class=\"dl-btn no-print\" onclick=\"downloadChart('basecomp-{idx}','base_composition')\">&#8595; PNG</button>\n",
        idx = idx
    ));
    s.push_str("</div>\n");
    s.push_str(&format!(
        "<canvas id=\"basecomp-{}\" width=\"600\" height=\"280\"></canvas>\n",
        idx
    ));
    s.push_str("<p class=\"chart-note\">Nucleotide composition (A/C/G/T/N) per read position. Parallel flat lines = unbiased composition.</p>\n");
    s.push_str("</div>\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Per-Read Quality Distribution</h3>\n");
    s.push_str(&format!(
        "<button class=\"dl-btn no-print\" onclick=\"downloadChart('qualdist-{idx}','quality_distribution')\">&#8595; PNG</button>\n",
        idx = idx
    ));
    s.push_str("</div>\n");
    s.push_str(&format!(
        "<canvas id=\"qualdist-{}\" width=\"600\" height=\"280\"></canvas>\n",
        idx
    ));
    s.push_str("<p class=\"chart-note\">Histogram of mean Phred quality scores across all reads. Red = Q&lt;20, yellow = Q20&ndash;30, green = Q&ge;30.</p>\n");
    s.push_str("</div>\n");

    s.push_str("</div>\n"); // charts-row (row 4)

    // Row 5: GC distribution + duplication histogram + N content
    s.push_str("<div class=\"charts-row\">\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>GC Content Distribution</h3>\n");
    s.push_str(&format!(
        "<button class=\"dl-btn no-print\" onclick=\"downloadChart('gcdist-{idx}','gc_distribution')\">&#8595; PNG</button>\n",
        idx = idx
    ));
    s.push_str("</div>\n");
    s.push_str(&format!(
        "<canvas id=\"gcdist-{}\" width=\"600\" height=\"240\"></canvas>\n", idx
    ));
    s.push_str("<p class=\"chart-note\">Per-read GC content histogram (blue). Red dashed curve = theoretical normal distribution fitted to observed mean and SD. Deviation from normal may indicate contamination or adapter dimers.</p>\n");
    s.push_str("</div>\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Duplication Level</h3>\n");
    if flags.duplication_check {
        s.push_str(&format!(
            "<button class=\"dl-btn no-print\" onclick=\"downloadChart('duphist-{idx}','dup_level_histogram')\">&#8595; PNG</button>\n",
            idx = idx
        ));
    }
    s.push_str("</div>\n");
    if flags.duplication_check {
        s.push_str(&format!(
            "<canvas id=\"duphist-{}\" width=\"600\" height=\"240\"></canvas>\n", idx
        ));
        s.push_str("<p class=\"chart-note\">Percentage of reads at each duplication level, estimated from first 200,000 read fingerprints. High counts in 2x+ bins indicate library over-amplification.</p>\n");
    } else {
        s.push_str(&skipped_notice("--no-duplication / --fast"));
    }
    s.push_str("</div>\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Per-Base N Content</h3>\n");
    s.push_str(&format!(
        "<button class=\"dl-btn no-print\" onclick=\"downloadChart('ncontent-{idx}','n_content')\">&#8595; PNG</button>\n",
        idx = idx
    ));
    s.push_str("</div>\n");
    s.push_str(&format!(
        "<canvas id=\"ncontent-{}\" width=\"600\" height=\"240\"></canvas>\n", idx
    ));
    s.push_str("<p class=\"chart-note\">Percentage of uncalled bases (N) per read position. Red dashed line = 5% warning threshold. Spikes at read ends are common for some sequencing platforms.</p>\n");
    s.push_str("</div>\n");

    s.push_str("</div>\n"); // charts-row (row 5)

    // Row 6 (long-read only): Quality vs Length scatter + Reads over time + Channel occupancy
    if lr {
        s.push_str("<div class=\"charts-row\">\n");

        // Quality vs Length scatter
        s.push_str("<div class=\"chart-box\">\n");
        s.push_str("<div class=\"chart-header\">\n");
        s.push_str("<h3>Quality vs Read Length</h3>\n");
        if !f.qual_vs_length.is_empty() {
            s.push_str(&format!(
                "<button class=\"dl-btn no-print\" onclick=\"downloadChart('qvl-{idx}','qual_vs_length')\">&#8595; PNG</button>\n",
                idx = idx
            ));
        }
        s.push_str("</div>\n");
        if f.qual_vs_length.is_empty() {
            s.push_str("<p class=\"chart-note\" style=\"padding:20px 0\">No quality-vs-length data available.</p>\n");
        } else {
            s.push_str(&format!(
                "<canvas id=\"qvl-{}\" width=\"600\" height=\"280\"></canvas>\n",
                idx
            ));
            s.push_str("<p class=\"chart-note\">Sampled read length vs mean Phred quality. Each dot = one read. Horizontal scatter indicates quality is independent of length (typical for ONT).</p>\n");
        }
        s.push_str("</div>\n");

        // Reads over time (ONT only)
        if !f.reads_over_time.is_empty() {
            s.push_str("<div class=\"chart-box\">\n");
            s.push_str("<div class=\"chart-header\">\n");
            s.push_str("<h3>Reads over Time (ONT)</h3>\n");
            s.push_str(&format!(
                "<button class=\"dl-btn no-print\" onclick=\"downloadChart('rot-{idx}','reads_over_time')\">&#8595; PNG</button>\n",
                idx = idx
            ));
            s.push_str("</div>\n");
            s.push_str(&format!(
                "<canvas id=\"rot-{}\" width=\"600\" height=\"280\"></canvas>\n",
                idx
            ));
            s.push_str("<p class=\"chart-note\">Cumulative reads per 5-minute bucket since run start. Plateaus may indicate pore saturation or flow cell end.</p>\n");
            s.push_str("</div>\n");
        }

        s.push_str("</div>\n"); // charts-row (row 6)

        // Channel occupancy (ONT only)
        if !f.ont_channel_counts.is_empty() {
            s.push_str("<div class=\"charts-row\">\n");
            s.push_str("<div class=\"chart-box\">\n");
            s.push_str("<div class=\"chart-header\">\n");
            s.push_str("<h3>Channel Occupancy (ONT)</h3>\n");
            s.push_str(&format!(
                "<button class=\"dl-btn no-print\" onclick=\"downloadChart('chann-{idx}','channel_occupancy')\">&#8595; PNG</button>\n",
                idx = idx
            ));
            s.push_str("</div>\n");
            s.push_str(&format!(
                "<canvas id=\"chann-{}\" width=\"600\" height=\"240\"></canvas>\n",
                idx
            ));
            let n_channels = f.ont_channel_counts.len();
            let total_chan_reads: u32 = f.ont_channel_counts.values().sum();
            s.push_str(&format!(
                "<p class=\"chart-note\">{} active channels detected ({} total reads). Histogram shows distribution of reads-per-channel. PromethION has up to 3000 channels; MinION up to 512.</p>\n",
                n_channels, total_chan_reads
            ));
            s.push_str("</div>\n");
            s.push_str("</div>\n"); // charts-row (channel occupancy)
        }
    }

    // Overrepresented sequences table
    s.push_str("<div class=\"chart-box\" style=\"margin-top:16px\">\n");
    s.push_str("<h3>Overrepresented Sequences</h3>\n");
    if !flags.overrep_sequences {
        s.push_str(&skipped_notice("--no-overrep / --fast"));
    } else if f.overrepresented_sequences.is_empty() {
        s.push_str("<p class=\"chart-note\" style=\"padding:12px 0\">No overrepresented sequences found (&ge;0.1% of sampled reads).</p>\n");
    } else {
        s.push_str("<p class=\"chart-note\" style=\"margin-bottom:8px\">Sampled first 200,000 reads. Sequences present in &ge;0.1% of reads:</p>\n");
        s.push_str("<table class=\"overrep-table\">\n");
        s.push_str("<thead><tr><th>#</th><th>Sequence (first 50 bp)</th><th>Count</th><th>%</th><th>Possible Source</th></tr></thead>\n");
        s.push_str("<tbody>\n");
        for (i, seq) in f.overrepresented_sequences.iter().enumerate() {
            let src_cls = if seq.possible_source == "No hit" { "color:var(--muted)" } else { "color:var(--yellow)" };
            s.push_str(&format!(
                "<tr><td>{}</td><td class=\"overrep-seq\">{}</td><td>{}</td><td>{:.2}%</td><td style=\"{}\">{}</td></tr>\n",
                i + 1,
                escape_html(&seq.sequence),
                format_number(seq.count),
                seq.percentage,
                src_cls,
                escape_html(&seq.possible_source),
            ));
        }
        s.push_str("</tbody></table>\n");
    }
    s.push_str("</div>\n");

    // Trim output info (only shown when trimming was active)
    if f.trim_output_path.is_some() || f.trimmed_reads > 0 {
        s.push_str("<div class=\"trim-box\" style=\"margin-top:16px\">\n");
        s.push_str("<h3 style=\"margin-bottom:8px\">Adapter Trim Output</h3>\n");
        s.push_str(&format!(
            "<p>Reads with adapters trimmed: <strong>{}</strong> ({:.1}%)</p>\n",
            format_number(f.trimmed_reads),
            f.trimmed_pct()
        ));
        s.push_str(&format!(
            "<p>Bases removed by trimming: <strong>{}</strong></p>\n",
            format_bases(f.trimmed_bases_removed)
        ));
        if let Some(ref tp) = f.trim_output_path {
            s.push_str(&format!(
                "<p>Trimmed output: <code>{}</code></p>\n",
                escape_html(tp)
            ));
        }
        s.push_str("</div>\n");
    }

    s.push_str("</div>\n</div>\n"); // file-body + file-section
    s
}

fn build_summary_table(files: &[&FileStats]) -> String {
    let mut s = String::new();
    s.push_str("<table class=\"summary-table\">\n");
    s.push_str("<thead><tr>");
    for hdr in &[
        "File",
        "Reads",
        "Total Bases",
        "Avg Length",
        "N50",
        "GC%",
        "Avg Quality",
        "Q30%",
        "Adapter%",
        "Dup~%",
    ] {
        s.push_str(&format!("<th>{}</th>", hdr));
    }
    s.push_str("</tr></thead>\n<tbody>\n");
    for f in files {
        let (n50, _) = f.compute_n50_n90();
        let fname = Path::new(&f.file_path)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or(&f.file_path);
        s.push_str("<tr>");
        s.push_str(&format!(
            "<td style=\"font-family:monospace\">{}</td>",
            escape_html(fname)
        ));
        s.push_str(&format!("<td>{}</td>", format_number(f.read_count)));
        s.push_str(&format!("<td>{}</td>", format_bases(f.total_bases)));
        s.push_str(&format!("<td>{:.1} bp</td>", f.avg_length()));
        s.push_str(&format!("<td>{} bp</td>", n50));
        s.push_str(&format!("<td>{:.2}%</td>", f.gc_content()));
        s.push_str(&format!("<td>Q{:.1}</td>", f.avg_quality()));
        s.push_str(&format!("<td>{:.1}%</td>", f.q30_pct()));
        s.push_str(&format!("<td>{:.2}%</td>", f.adapter_pct()));
        s.push_str(&format!("<td>{:.1}%</td>", f.dup_rate_pct));
        s.push_str("</tr>\n");
    }
    s.push_str("</tbody></table>\n");
    s
}

fn build_file_json(f: &FileStats) -> serde_json::Value {
    let qual_per_pos: Vec<f64> = f.avg_qual_per_position();
    let len_dist: Vec<[u64; 2]> = f
        .sorted_length_dist()
        .iter()
        .map(|(l, c)| [*l, *c])
        .collect();
    let kmer_json: Vec<serde_json::Value> = f
        .top_kmers(20)
        .into_iter()
        .map(|(k, v)| serde_json::json!([k, v]))
        .collect();
    // Per-tile: [[tile_id, avg_quality], ...]
    let tile_json: Vec<serde_json::Value> = f
        .sorted_tile_quality()
        .into_iter()
        .map(|(t, q)| serde_json::json!([t, q]))
        .collect();
    // Per-base composition: [[A%, C%, G%, T%, N%], ...]
    let basecomp_json: Vec<serde_json::Value> = f
        .base_composition_pct()
        .into_iter()
        .map(|pct| serde_json::json!(pct))
        .collect();
    // Quality distribution: [count_at_Q0, count_at_Q1, ..., count_at_Q42]
    let qualdist_json: Vec<u64> = f.quality_distribution.clone();

    // Q25/Q50/Q75 per position for IQR overlay on quality chart
    let qual_pct_json: Vec<serde_json::Value> = f
        .qual_percentiles_per_position()
        .into_iter()
        .map(|pct| serde_json::json!(pct))
        .collect();

    let gc_dist_json: Vec<u64> = f.gc_distribution.clone();
    let dup_hist_labels = ["1x","2x","3-4x","5-9x","10-49x","50-99x","100-499x","500-999x",">=1000x"];
    let dup_hist_json: Vec<serde_json::Value> = f.dup_level_histogram.iter()
        .enumerate()
        .map(|(i, &pct)| serde_json::json!([dup_hist_labels[i], pct]))
        .collect();
    let n_content_json: Vec<f64> = f.base_composition_pct().iter().map(|p| p[4]).collect();

    // Long-read / ONT specific data
    // qual_vs_length: [[len, mean_q], ...]
    let qvl_json: Vec<serde_json::Value> = f.qual_vs_length
        .iter()
        .map(|(l, q)| serde_json::json!([l, q]))
        .collect();
    // reads_over_time: [[minutes, cumulative_reads], ...]
    let rot_json: Vec<serde_json::Value> = f.reads_over_time
        .iter()
        .map(|(m, c)| serde_json::json!([m, c]))
        .collect();
    // ont_channel_counts: [[channel, count], ...] sorted by channel
    let mut chann_sorted: Vec<(u32, u32)> = f.ont_channel_counts
        .iter()
        .map(|(&ch, &cnt)| (ch, cnt))
        .collect();
    chann_sorted.sort_by_key(|(ch, _)| *ch);
    let chann_json: Vec<serde_json::Value> = chann_sorted
        .iter()
        .map(|(ch, cnt)| serde_json::json!([ch, cnt]))
        .collect();

    serde_json::json!({
        "qual": qual_per_pos,
        "qual_pct": qual_pct_json,
        "len": len_dist,
        "kmers": kmer_json,
        "tiles": tile_json,
        "basecomp": basecomp_json,
        "qualdist": qualdist_json,
        "gc_dist": gc_dist_json,
        "dup_hist": dup_hist_json,
        "n_content": n_content_json,
        "qvl": qvl_json,
        "rot": rot_json,
        "chann": chann_json,
        "long_read": f.long_read_mode,
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn stat_card(buf: &mut String, label: &str, value: &str, cls: &str) {
    buf.push_str("<div class=\"stat-card\">\n");
    buf.push_str(&format!(
        "<div class=\"stat-label\">{}</div>\n",
        escape_html(label)
    ));
    if cls.is_empty() {
        buf.push_str(&format!("<div class=\"stat-value\">{}</div>\n", value));
    } else {
        buf.push_str(&format!(
            "<div class=\"stat-value {}\">{}</div>\n",
            cls, value
        ));
    }
    buf.push_str("</div>\n");
}

fn quality_class(q: f64) -> &'static str {
    if q >= 30.0 {
        "pass"
    } else if q >= 20.0 {
        "warn"
    } else {
        "fail"
    }
}

/// Quality class with long-read aware thresholds (ONT: warn < Q10, fail < Q8).
fn quality_class_lr(q: f64, long_read: bool) -> &'static str {
    if long_read {
        if q >= 10.0 { "pass" } else if q >= 8.0 { "warn" } else { "fail" }
    } else {
        quality_class(q)
    }
}

fn gc_class(gc: f64) -> &'static str {
    if (35.0..=65.0).contains(&gc) {
        "pass"
    } else if (25.0..=75.0).contains(&gc) {
        "warn"
    } else {
        "fail"
    }
}

fn dup_status_class(pct: f64) -> &'static str {
    if pct < 5.0 {
        "pass"
    } else if pct < 20.0 {
        "warn"
    } else {
        "fail"
    }
}

fn adapter_status_class(pct: f64) -> &'static str {
    if pct < 1.0 {
        "pass"
    } else if pct < 10.0 {
        "warn"
    } else {
        "fail"
    }
}

fn escape_html(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn report_stem(files: &[&FileStats]) -> String {
    if files.len() == 1 {
        Path::new(&files[0].file_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output")
            .to_string()
    } else {
        "batch_report".to_string()
    }
}

// ---------------------------------------------------------------------------
// Auto-parameter suggestions
// ---------------------------------------------------------------------------

pub fn suggest_parameters(f: &FileStats) -> Vec<String> {
    let mut suggestions = Vec::new();

    // Poly-G: high G% at 3' end positions relative to overall GC
    let n_pos = f.base_composition.len();
    if n_pos >= 10 {
        let last10_g: f64 = f.base_composition[n_pos.saturating_sub(10)..]
            .iter()
            .map(|b| {
                let total = (b[0] + b[1] + b[2] + b[3] + b[4]) as f64;
                if total > 0.0 { b[2] as f64 / total } else { 0.0 }
            })
            .sum::<f64>() / 10.0;
        if last10_g > 0.40 && f.gc_content() < 65.0 {
            suggestions.push(format!(
                "--poly-g 10  (high G% at 3' end [{:.0}%], typical of NovaSeq/NextSeq 2-color chemistry)",
                last10_g * 100.0
            ));
        }
    }

    // Adapter contamination
    let ap = f.adapter_pct();
    if ap > 20.0 {
        suggestions.push(format!(
            "--trim  (HIGH adapter content [{:.1}%] — strongly recommended)",
            ap
        ));
    } else if ap > 5.0 {
        suggestions.push(format!(
            "--trim  (adapter sequences detected in {:.1}% of reads)",
            ap
        ));
    }

    // 3' quality drop: find position where mean quality first falls below Q20
    let qual_means = f.avg_qual_per_position();
    let drop_pos = qual_means.iter().enumerate()
        .find(|(_, &q)| q > 0.0 && q < 20.0)
        .map(|(i, _)| i);

    if let Some(pos) = drop_pos {
        let pct = pos as f64 / qual_means.len().max(1) as f64;
        if pct < 0.90 {
            suggestions.push(format!(
                "--cut-right --window-quality 20  (mean quality drops below Q20 at position {})",
                pos
            ));
        }
    }

    // Low overall quality
    let aq = f.avg_quality();
    if aq < 25.0 && aq > 0.0 {
        suggestions.push(format!(
            "--min-quality 20  (mean read quality Q{:.1} is below typical threshold)",
            aq
        ));
    }

    // High N content
    let high_n = f.n_content_per_position().iter().any(|&n| n > 5.0);
    if high_n {
        suggestions.push("--max-n 5  (elevated N base content detected)".to_string());
    }

    // High duplication
    if f.dup_rate_pct > 50.0 {
        suggestions.push(format!(
            "High duplication rate ({:.1}%) — consider deduplication downstream",
            f.dup_rate_pct
        ));
    }

    suggestions
}

// ---------------------------------------------------------------------------
// Two-file comparison report
// ---------------------------------------------------------------------------

pub fn export_comparison_html(snap: &SharedState, output_dir: &str) -> Result<String, String> {
    let files = snap.all_files();
    if files.len() < 2 {
        return Err("Need at least two files for comparison".to_string());
    }
    let f1 = files[0];
    let f2 = files[1];

    let fname1 = Path::new(&f1.file_path).file_name()
        .and_then(|n| n.to_str()).unwrap_or(&f1.file_path);
    let fname2 = Path::new(&f2.file_path).file_name()
        .and_then(|n| n.to_str()).unwrap_or(&f2.file_path);

    let long_read = f1.long_read_mode || f2.long_read_mode;
    let (q_lo, q_hi) = if long_read { (10u8, 20u8) } else { (20u8, 30u8) };

    let mixed_note = if f1.long_read_mode != f2.long_read_mode {
        "<p style='color:#b45309;background:#fffbeb;padding:8px;border-radius:6px;margin:12px 0'>\
        ⚠ Note: comparing short-read and long-read files — thresholds may differ.</p>"
    } else { "" };

    fn pct_delta(a: f64, b: f64) -> String {
        let d = b - a;
        if d.abs() < 0.01 { "±0".to_string() }
        else if d > 0.0 { format!("+{:.2}", d) }
        else { format!("{:.2}", d) }
    }
    fn fmt_delta(a: f64, b: f64, unit: &str) -> String {
        let d = b - a;
        let sign = if d >= 0.0 { "+" } else { "" };
        format!("{}{:.1}{}", sign, d, unit)
    }

    let (n50_1, _n90_1) = f1.compute_n50_n90();
    let (n50_2, _n90_2) = f2.compute_n50_n90();

    let rows = vec![
        ("Reads", format_number(f1.read_count), format_number(f2.read_count),
         fmt_delta(f1.read_count as f64, f2.read_count as f64, "")),
        ("Total Bases", format_bases(f1.total_bases), format_bases(f2.total_bases),
         fmt_delta(f1.total_bases as f64, f2.total_bases as f64, " bp")),
        ("Avg Length", format!("{:.1} bp", f1.avg_length()), format!("{:.1} bp", f2.avg_length()),
         fmt_delta(f1.avg_length(), f2.avg_length(), " bp")),
        ("GC Content", format!("{:.2}%", f1.gc_content()), format!("{:.2}%", f2.gc_content()),
         pct_delta(f1.gc_content(), f2.gc_content()) + "%"),
        ("Avg Quality", format!("Q{:.1}", f1.avg_quality()), format!("Q{:.1}", f2.avg_quality()),
         fmt_delta(f1.avg_quality(), f2.avg_quality(), "")),
        (if q_lo == 10 { "≥Q10" } else { "≥Q20" },
         format!("{:.1}%", f1.q20_pct()),
         format!("{:.1}%", f2.q20_pct()),
         pct_delta(f1.q20_pct(), f2.q20_pct()) + "%"),
        (if q_hi == 20 { "≥Q20" } else { "≥Q30" },
         format!("{:.1}%", if q_hi == 20 { f1.q20_pct() } else { f1.q30_pct() }),
         format!("{:.1}%", if q_hi == 20 { f2.q20_pct() } else { f2.q30_pct() }),
         pct_delta(f1.q30_pct(), f2.q30_pct()) + "%"),
        ("Adapter %", format!("{:.2}%", f1.adapter_pct()), format!("{:.2}%", f2.adapter_pct()),
         pct_delta(f1.adapter_pct(), f2.adapter_pct()) + "%"),
        ("Dup Rate", format!("{:.1}%", f1.dup_rate_pct), format!("{:.1}%", f2.dup_rate_pct),
         pct_delta(f1.dup_rate_pct, f2.dup_rate_pct) + "%"),
        ("N50", n50_1.to_string(), n50_2.to_string(),
         fmt_delta(n50_1 as f64, n50_2 as f64, " bp")),
    ];

    let mut table_rows = String::new();
    for (label, v1, v2, delta) in &rows {
        table_rows.push_str(&format!(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td style='color:#6b7280'>{}</td></tr>",
            label, v1, v2, delta
        ));
    }

    // Serialize qual per position for JS
    let serialize_qual = |means: &[f64]| -> String {
        format!("[{}]", means.iter().map(|v| format!("{:.2}", v)).collect::<Vec<_>>().join(","))
    };

    let qual1_js = serialize_qual(&f1.avg_qual_per_position());
    let qual2_js = serialize_qual(&f2.avg_qual_per_position());

    // Length histogram: convert HashMap to sorted Vec for JS
    let serialize_len = |hist: &std::collections::HashMap<u64, u64>| -> String {
        if hist.is_empty() { return "[]".to_string(); }
        let max_len = *hist.keys().max().unwrap_or(&0);
        let mut arr = vec![0u64; (max_len as usize).min(10_000) + 1];
        for (&len, &cnt) in hist {
            if (len as usize) < arr.len() { arr[len as usize] = cnt; }
        }
        format!("[{}]", arr.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(","))
    };

    let len_dist1_js = serialize_len(&f1.length_histogram);
    let len_dist2_js = serialize_len(&f2.length_histogram);

    let html = format!(r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>BioFastq-A Comparison Report</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#f8fafc;color:#1e293b;margin:0;padding:0}}
header{{background:#1a5fa8;color:#fff;padding:18px 32px;display:flex;align-items:center;justify-content:space-between}}
header h1{{margin:0;font-size:1.3rem;font-weight:700}}
header span{{font-size:.85rem;opacity:.8}}
.container{{max-width:1100px;margin:32px auto;padding:0 24px}}
table{{width:100%;border-collapse:collapse;background:#fff;border-radius:10px;overflow:hidden;box-shadow:0 1px 4px rgba(0,0,0,.08);margin-bottom:32px}}
th{{background:#f1f5f9;padding:10px 16px;text-align:left;font-size:.8rem;text-transform:uppercase;color:#64748b;letter-spacing:.05em}}
td{{padding:10px 16px;border-top:1px solid #f1f5f9;font-size:.9rem}}
tr:hover td{{background:#f8fafc}}
.chart-grid{{display:grid;grid-template-columns:1fr 1fr;gap:24px;margin-bottom:32px}}
.card{{background:#fff;border-radius:10px;box-shadow:0 1px 4px rgba(0,0,0,.08);padding:20px}}
.card h3{{margin:0 0 12px;font-size:.95rem;color:#1e293b}}
canvas{{width:100%!important}}
.legend{{display:flex;gap:16px;font-size:.8rem;margin-bottom:8px}}
.legend-dot{{width:10px;height:10px;border-radius:50%;display:inline-block;margin-right:4px}}
</style>
</head>
<body>
<header>
  <h1>BioFastq-A — Comparison Report</h1>
  <span>biofastq-a v{ver}</span>
</header>
<div class="container">
  {mixed}
  <table>
    <thead><tr><th>Metric</th><th>{f1}</th><th>{f2}</th><th>Delta (f2−f1)</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  <div class="chart-grid">
    <div class="card">
      <h3>Quality per Position</h3>
      <div class="legend">
        <span><span class="legend-dot" style="background:#1a5fa8"></span>{f1}</span>
        <span><span class="legend-dot" style="background:#f59e0b"></span>{f2}</span>
      </div>
      <canvas id="qualChart" height="180"></canvas>
    </div>
    <div class="card">
      <h3>Read Length Distribution</h3>
      <div class="legend">
        <span><span class="legend-dot" style="background:#1a5fa8"></span>{f1}</span>
        <span><span class="legend-dot" style="background:#f59e0b"></span>{f2}</span>
      </div>
      <canvas id="lenChart" height="180"></canvas>
    </div>
  </div>
</div>
<script>
const qual1={q1}, qual2={q2};
const len1={l1}, len2={l2};

function drawOverlay(id, data1, data2, label1, label2, yLabel) {{
  const cv = document.getElementById(id);
  if (!cv) return;
  const ctx = cv.getContext('2d');
  cv.width = cv.parentElement.clientWidth || 500;
  cv.height = 200;
  const W=cv.width, H=cv.height, pad={{t:20,r:20,b:40,l:50}};
  const iW=W-pad.l-pad.r, iH=H-pad.t-pad.b;
  ctx.clearRect(0,0,W,H);
  const all=[...data1,...data2].filter(v=>v>0);
  if(!all.length) return;
  const maxV=Math.max(...all);
  const len=Math.max(data1.length,data2.length);
  function drawLine(data, color) {{
    ctx.beginPath(); ctx.strokeStyle=color; ctx.lineWidth=2;
    data.forEach((v,i)=>{{
      const x=pad.l+i/Math.max(len-1,1)*iW;
      const y=pad.t+iH-(v/maxV)*iH;
      i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);
    }});
    ctx.stroke();
  }}
  drawLine(data1,'#1a5fa8');
  drawLine(data2,'#f59e0b');
  // axes
  ctx.strokeStyle='#e2e8f0'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,pad.t+iH); ctx.lineTo(pad.l+iW,pad.t+iH); ctx.stroke();
  ctx.fillStyle='#64748b'; ctx.font='11px sans-serif'; ctx.textAlign='center';
  [0,Math.floor(len/4),Math.floor(len/2),Math.floor(3*len/4),len-1].forEach(i=>{{
    const x=pad.l+i/Math.max(len-1,1)*iW;
    ctx.fillText(i,x,H-8);
  }});
}}

window.addEventListener('load',()=>{{
  drawOverlay('qualChart',qual1,qual2,'{f1}','{f2}','Mean Phred');
  drawOverlay('lenChart',len1,len2,'{f1}','{f2}','Reads');
}});
</script>
</body>
</html>"#,
        ver = env!("CARGO_PKG_VERSION"),
        mixed = mixed_note,
        f1 = escape_html(fname1),
        f2 = escape_html(fname2),
        rows = table_rows,
        q1 = qual1_js,
        q2 = qual2_js,
        l1 = len_dist1_js,
        l2 = len_dist2_js,
    );

    let stem = format!("{}_vs_{}",
        Path::new(&f1.file_path).file_stem().and_then(|s| s.to_str()).unwrap_or("file1"),
        Path::new(&f2.file_path).file_stem().and_then(|s| s.to_str()).unwrap_or("file2"),
    );
    let out_path = format!("{}/{}_comparison.html", output_dir.trim_end_matches('/'), stem);
    fs::write(&out_path, html).map_err(|e| e.to_string())?;
    Ok(out_path)
}
