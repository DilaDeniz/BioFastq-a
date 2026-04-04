use std::fs;
use std::io;
use std::path::Path;

use serde::Serialize;
use serde_json;

use crate::types::{format_bases, format_number, FileStats, SharedState, ADAPTERS};

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
        version: "2.1.0",
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
// HTML report
// ---------------------------------------------------------------------------

/// Static CSS embedded in the report (dark GitHub-inspired theme)
const CSS: &str = r#"
:root{--bg:#0d1117;--card:#161b22;--card2:#1c2128;--border:#30363d;--text:#c9d1d9;--muted:#8b949e;--accent:#58a6ff;--green:#3fb950;--yellow:#d29922;--orange:#e3b341;--red:#f85149;--purple:#bc8cff;}
*{box-sizing:border-box;margin:0;padding:0;}
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:var(--bg);color:var(--text);padding:24px;line-height:1.6;font-size:14px;}
a{color:var(--accent);text-decoration:none;}
h1{font-size:26px;font-weight:700;color:var(--accent);margin-bottom:4px;}
h2{font-size:18px;font-weight:600;color:var(--text);margin-bottom:12px;}
h3{font-size:13px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.06em;margin-bottom:10px;}
.header{border-bottom:1px solid var(--border);padding-bottom:20px;margin-bottom:28px;}
.header-meta{color:var(--muted);font-size:13px;margin-top:6px;display:flex;gap:20px;flex-wrap:wrap;}
.badge{display:inline-block;padding:2px 8px;border-radius:12px;font-size:11px;font-weight:600;background:var(--card2);border:1px solid var(--border);}
.badge-green{border-color:var(--green);color:var(--green);}
.badge-yellow{border-color:var(--yellow);color:var(--yellow);}
.badge-red{border-color:var(--red);color:var(--red);}
/* Summary table */
.summary-table{width:100%;border-collapse:collapse;margin-bottom:28px;font-size:13px;}
.summary-table th{background:var(--card2);color:var(--muted);font-weight:600;text-align:left;padding:8px 12px;border:1px solid var(--border);}
.summary-table td{padding:8px 12px;border:1px solid var(--border);background:var(--card);}
.summary-table tr:hover td{background:var(--card2);}
/* File section */
.file-section{border:1px solid var(--border);border-radius:10px;margin-bottom:32px;overflow:hidden;}
.file-section-header{background:var(--card2);padding:16px 20px;border-bottom:1px solid var(--border);}
.file-section-header .filename{font-size:16px;font-weight:700;font-family:monospace;color:var(--text);}
.file-section-header .filepath{font-size:11px;color:var(--muted);margin-top:2px;font-family:monospace;}
.file-body{padding:20px;}
/* Stats grid */
.stats-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(155px,1fr));gap:10px;margin-bottom:20px;}
.stat-card{background:var(--card2);border:1px solid var(--border);border-radius:8px;padding:12px 14px;}
.stat-label{font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.05em;margin-bottom:4px;}
.stat-value{font-size:20px;font-weight:700;font-family:monospace;color:var(--text);}
.stat-value.pass{color:var(--green);}
.stat-value.warn{color:var(--yellow);}
.stat-value.fail{color:var(--red);}
/* Charts */
.charts-row{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:16px;margin-bottom:16px;}
.chart-box{background:var(--card2);border:1px solid var(--border);border-radius:8px;padding:16px;}
.chart-box canvas{max-width:100%;display:block;}
.chart-note{font-size:11px;color:var(--muted);margin-top:6px;}
/* Adapter section */
.adapter-grid{display:grid;grid-template-columns:auto 1fr;gap:24px;align-items:start;}
.big-pct{font-size:40px;font-weight:700;font-family:monospace;line-height:1;}
.adapter-list{margin-top:8px;padding:0;}
.adapter-list li{list-style:none;font-size:12px;font-family:monospace;color:var(--muted);padding:2px 0;}
.adapter-list li::before{content:"▸ ";color:var(--accent);}
/* Duplication gauge */
.dup-gauge-wrap{display:flex;align-items:center;gap:16px;}
.dup-gauge-track{flex:1;height:16px;background:var(--card);border:1px solid var(--border);border-radius:8px;overflow:hidden;}
.dup-gauge-fill{height:100%;border-radius:8px;transition:width .3s;}
/* Per-tile table */
.tile-table{width:100%;border-collapse:collapse;font-size:12px;font-family:monospace;}
.tile-table th{background:var(--card2);color:var(--muted);padding:5px 10px;text-align:left;border:1px solid var(--border);}
.tile-table td{padding:4px 10px;border:1px solid var(--border);}
/* Trim info */
.trim-box{background:var(--card2);border:1px solid var(--border);border-radius:8px;padding:16px;}
.trim-box p{font-size:13px;color:var(--muted);margin:4px 0;}
.trim-box code{color:var(--accent);font-family:monospace;font-size:12px;}
/* Footer */
.footer{text-align:center;color:var(--muted);font-size:12px;margin-top:40px;padding-top:20px;border-top:1px solid var(--border);}
"#;

/// Static JavaScript for chart drawing (uses vanilla Canvas API, no deps)
const JS: &str = r#"
'use strict';
var RD = __REPORT_DATA__;

/* ---------- Quality per position ---------- */
function drawQual(id, data) {
  var c = document.getElementById(id);
  if (!c || !data || !data.length) return;
  var ctx = c.getContext('2d');
  var W = c.clientWidth || c.width; c.width = W;
  var H = c.height;
  var p = {t:18,r:12,b:36,l:44};
  var pw = W-p.l-p.r, ph = H-p.t-p.b, MQ = 42;
  ctx.clearRect(0,0,W,H);

  /* quality zone backgrounds */
  [[0,20,'rgba(248,81,73,.10)'],[20,28,'rgba(211,153,52,.10)'],[28,MQ,'rgba(63,185,80,.07)']].forEach(function(z){
    var y1=p.t+(1-z[1]/MQ)*ph, y2=p.t+(1-z[0]/MQ)*ph;
    ctx.fillStyle=z[2]; ctx.fillRect(p.l,y1,pw,y2-y1);
  });

  /* dashed threshold lines */
  [[20,'rgba(248,81,73,.45)'],[28,'rgba(211,153,52,.45)'],[30,'rgba(63,185,80,.35)']].forEach(function(qc){
    var y=p.t+(1-qc[0]/MQ)*ph;
    ctx.setLineDash([3,4]); ctx.strokeStyle=qc[1]; ctx.lineWidth=1;
    ctx.beginPath(); ctx.moveTo(p.l,y); ctx.lineTo(W-p.r,y); ctx.stroke();
  });
  ctx.setLineDash([]);

  /* grid lines */
  [0,10,20,30,40].forEach(function(q){
    var y=p.t+(1-q/MQ)*ph;
    ctx.strokeStyle='rgba(48,54,61,0.6)'; ctx.lineWidth=1;
    ctx.beginPath(); ctx.moveTo(p.l,y); ctx.lineTo(W-p.r,y); ctx.stroke();
  });

  /* filled area + line */
  if (data.length > 1) {
    var grad = ctx.createLinearGradient(0,p.t,0,H-p.b);
    grad.addColorStop(0,'rgba(88,166,255,.25)');
    grad.addColorStop(1,'rgba(88,166,255,0)');
    ctx.beginPath();
    data.forEach(function(q,i){
      var x=p.l+(i/(data.length-1))*pw, y=p.t+(1-Math.min(q,MQ)/MQ)*ph;
      i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);
    });
    ctx.lineTo(p.l+pw,H-p.b); ctx.lineTo(p.l,H-p.b);
    ctx.closePath(); ctx.fillStyle=grad; ctx.fill();

    ctx.beginPath();
    data.forEach(function(q,i){
      var x=p.l+(i/(data.length-1))*pw, y=p.t+(1-Math.min(q,MQ)/MQ)*ph;
      i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);
    });
    ctx.strokeStyle='#58a6ff'; ctx.lineWidth=2; ctx.stroke();
  }

  /* axes */
  ctx.strokeStyle='#30363d'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(p.l,p.t); ctx.lineTo(p.l,H-p.b); ctx.lineTo(W-p.r,H-p.b); ctx.stroke();

  /* Y labels */
  ctx.fillStyle='#8b949e'; ctx.font='11px monospace'; ctx.textAlign='right';
  [0,10,20,28,30,40].forEach(function(q){
    ctx.fillText(q,p.l-5,p.t+(1-q/MQ)*ph+4);
  });

  /* X labels */
  ctx.textAlign='center';
  var n = data.length;
  var step = n > 300 ? 100 : n > 150 ? 50 : n > 50 ? 25 : 10;
  for (var i=0; i<n; i+=step) {
    ctx.fillText(i+1, p.l+(i/Math.max(n-1,1))*pw, H-p.b+14);
  }
  ctx.fillText(n, p.l+pw, H-p.b+14);

  /* axis titles */
  ctx.fillStyle='#8b949e'; ctx.font='11px monospace'; ctx.textAlign='center';
  ctx.fillText('Position (bp)', p.l+pw/2, H-2);
  ctx.save(); ctx.translate(10,p.t+ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('Phred Quality', 0, 0); ctx.restore();
}

/* ---------- Read length distribution ---------- */
function drawLen(id, data) {
  var c = document.getElementById(id);
  if (!c || !data || !data.length) return;
  var ctx = c.getContext('2d');
  var W = c.clientWidth || c.width; c.width = W;
  var H = c.height;
  var p = {t:18,r:12,b:40,l:60};
  var pw = W-p.l-p.r, ph = H-p.t-p.b;
  ctx.clearRect(0,0,W,H);

  var maxC = 0;
  data.forEach(function(d){ if(d[1]>maxC) maxC=d[1]; });
  var n = data.length;
  var bw = pw / n;

  data.forEach(function(d,i){
    var x = p.l + i*bw;
    var bh = (d[1]/maxC)*ph;
    var y = p.t+ph-bh;
    var grad = ctx.createLinearGradient(0,y,0,y+bh);
    grad.addColorStop(0,'rgba(88,166,255,.85)');
    grad.addColorStop(1,'rgba(88,166,255,.45)');
    ctx.fillStyle = grad;
    ctx.fillRect(x+0.5, y, Math.max(1, bw-1), bh);
  });

  /* axes */
  ctx.strokeStyle='#30363d'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(p.l,p.t); ctx.lineTo(p.l,H-p.b); ctx.lineTo(W-p.r,H-p.b); ctx.stroke();

  /* X labels */
  ctx.fillStyle='#8b949e'; ctx.font='11px monospace'; ctx.textAlign='center';
  var lblN = Math.min(6, n);
  for (var i=0; i<lblN; i++) {
    var idx = Math.round(i*(n-1)/(Math.max(lblN-1,1)));
    var x = p.l+(idx+0.5)*bw;
    var len = data[idx][0];
    ctx.fillText(len>999?(len/1000).toFixed(1)+'k':len+'bp', x, H-p.b+14);
  }

  /* Y labels */
  ctx.textAlign='right';
  [0,.5,1].forEach(function(f){
    var v=Math.round(maxC*f);
    var label = v>999999?(v/1000000).toFixed(1)+'M':v>999?(v/1000).toFixed(0)+'k':v;
    ctx.fillText(label, p.l-5, p.t+(1-f)*ph+4);
  });

  /* axis titles */
  ctx.textAlign='center';
  ctx.fillText('Read Length', p.l+pw/2, H-2);
  ctx.save(); ctx.translate(10,p.t+ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('Count', 0, 0); ctx.restore();
}

/* ---------- K-mer frequency ---------- */
function drawKmer(id, data) {
  var c = document.getElementById(id);
  if (!c || !data || !data.length) return;
  var H = Math.max(260, data.length*22+20); c.height = H;
  var ctx = c.getContext('2d');
  var W = c.clientWidth || c.width; c.width = W;
  var p = {t:10,r:90,b:10,l:68};
  var pw = W-p.l-p.r, ph = H-p.t-p.b;
  ctx.clearRect(0,0,W,H);

  var maxC = data[0][1];
  var rowH = ph/data.length;
  var barH = Math.max(10, rowH-4);
  var ntCol = {A:'#3fb950',C:'#58a6ff',G:'#d29922',T:'#f85149',N:'#8b949e'};

  data.forEach(function(d,i){
    var kmer=d[0], cnt=d[1];
    var y = p.t + i*rowH + (rowH-barH)/2;
    var bw = (cnt/maxC)*pw;

    /* bar */
    ctx.fillStyle='rgba(88,166,255,0.65)';
    ctx.fillRect(p.l, y, bw, barH);

    /* k-mer label (colored by nucleotide) */
    ctx.font = 'bold 12px monospace';
    ctx.textAlign = 'right';
    var cw = 8.4;
    var startX = p.l - 5 - (kmer.length-1)*cw;
    for (var j=0; j<kmer.length; j++) {
      var nt = kmer[j];
      ctx.fillStyle = ntCol[nt] || '#c9d1d9';
      ctx.fillText(nt, startX + j*cw + cw, y+barH-3);
    }

    /* count label */
    ctx.fillStyle='#8b949e'; ctx.font='11px monospace'; ctx.textAlign='left';
    ctx.fillText(cnt.toLocaleString(), p.l+bw+5, y+barH-3);
  });
}

/* ---------- Per-tile quality ---------- */
function drawTile(id, data) {
  var c = document.getElementById(id);
  if (!c || !data || !data.length) return;
  var ctx = c.getContext('2d');
  var W = c.clientWidth || c.width; c.width = W;
  var H = c.height;
  var p = {t:18,r:12,b:36,l:60};
  var pw = W-p.l-p.r, ph = H-p.t-p.b, MQ = 42;
  ctx.clearRect(0,0,W,H);

  var n = data.length;
  var bw = pw / n;
  var minQ = Math.min.apply(null, data.map(function(d){return d[1];}));
  var maxQ = Math.max.apply(null, data.map(function(d){return d[1];}));
  var rng = maxQ - minQ || 1;

  data.forEach(function(d,i){
    var q = d[1];
    var x = p.l + i*bw;
    var bh = ((q - minQ) / rng) * ph;
    var y = p.t + ph - bh;
    /* colour: red→yellow→green */
    var t = (q - minQ) / rng;
    var r = Math.round(248*(1-t) + 63*t);
    var g = Math.round(81*(1-t) + 185*t);
    var b2 = Math.round(73*(1-t) + 80*t);
    ctx.fillStyle = 'rgba('+r+','+g+','+b2+',0.8)';
    ctx.fillRect(x+0.5, y, Math.max(1,bw-1), bh);
  });

  ctx.strokeStyle='#30363d'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(p.l,p.t); ctx.lineTo(p.l,H-p.b); ctx.lineTo(W-p.r,H-p.b); ctx.stroke();

  ctx.fillStyle='#8b949e'; ctx.font='11px monospace'; ctx.textAlign='right';
  [minQ, (minQ+maxQ)/2, maxQ].forEach(function(q){
    var y=p.t+(1-(q-minQ)/rng)*ph;
    ctx.fillText(q.toFixed(1),p.l-5,y+4);
  });

  ctx.textAlign='center';
  var lblN = Math.min(8, n);
  for(var i=0; i<lblN; i++){
    var idx=Math.round(i*(n-1)/(Math.max(lblN-1,1)));
    var x=p.l+(idx+0.5)*bw;
    ctx.fillText(data[idx][0], x, H-p.b+14);
  }
  ctx.fillText('Tile ID', p.l+pw/2, H-2);
}

/* ---------- Base composition chart ---------- */
function drawBaseComp(id, data) {
  // data: array of [A%, C%, G%, T%, N%] per position
  var cv=document.getElementById(id); if(!cv)return;
  var ctx=cv.getContext('2d'), W=cv.width, H=cv.height;
  var p={l:52,r:16,t:20,b:36};
  var pw=W-p.l-p.r, ph=H-p.t-p.b;
  ctx.fillStyle='#0d1117'; ctx.fillRect(0,0,W,H);
  if(!data||!data.length){ctx.fillStyle='#8b949e';ctx.font='14px monospace';ctx.fillText('No data',W/2-30,H/2);return;}
  var colors=['#58a6ff','#3fb950','#f78166','#d29922','#8b949e']; // A,C,G,T,N
  var labels=['A','C','G','T','N'];
  var n=data.length, bw=pw/n;
  // grid
  ctx.strokeStyle='#21262d'; ctx.lineWidth=1;
  for(var g=0;g<=4;g++){var y=p.t+ph*(1-g/4);ctx.beginPath();ctx.moveTo(p.l,y);ctx.lineTo(p.l+pw,y);ctx.stroke();}
  // draw lines
  labels.forEach(function(lbl,li){
    ctx.strokeStyle=colors[li]; ctx.lineWidth=1.5; ctx.beginPath();
    data.forEach(function(pt,xi){
      var x=p.l+(xi+0.5)*bw, y=p.t+ph*(1-pt[li]/100);
      if(xi===0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
    });
    ctx.stroke();
  });
  // axes
  ctx.fillStyle='#8b949e'; ctx.font='11px monospace';
  for(var g=0;g<=4;g++){ctx.fillText((g*25)+'%',2,p.t+ph*(1-g/4)+4);}
  ctx.fillText('1',p.l,H-p.b+14);
  ctx.fillText(n,p.l+pw-20,H-p.b+14);
  // legend
  var lx=p.l+4;
  labels.forEach(function(lbl,li){
    ctx.fillStyle=colors[li]; ctx.fillRect(lx,p.t+2,10,10);
    ctx.fillStyle='#c9d1d9'; ctx.fillText(lbl,lx+13,p.t+11); lx+=32;
  });
}

/* ---------- Quality distribution chart ---------- */
function drawQualDist(id, data) {
  // data: array of 43 counts indexed by mean Phred score
  var cv=document.getElementById(id); if(!cv)return;
  var ctx=cv.getContext('2d'), W=cv.width, H=cv.height;
  var p={l:52,r:16,t:20,b:36};
  var pw=W-p.l-p.r, ph=H-p.t-p.b;
  ctx.fillStyle='#0d1117'; ctx.fillRect(0,0,W,H);
  if(!data||!data.length){ctx.fillStyle='#8b949e';ctx.font='14px monospace';ctx.fillText('No data',W/2-30,H/2);return;}
  var max=Math.max.apply(null,data)||1;
  var n=data.length, bw=pw/n;
  data.forEach(function(cnt,q){
    var pct=cnt/max;
    var r=q<20?200:q<30?210:80, g=q<20?80:q<30?160:200, b=q<20?80:q<30?80:100;
    ctx.fillStyle='rgba('+r+','+g+','+b+',0.85)';
    var barH=ph*pct;
    ctx.fillRect(p.l+q*bw, p.t+ph-barH, bw-1, barH);
  });
  // grid lines
  ctx.strokeStyle='#21262d'; ctx.lineWidth=1;
  for(var g=0;g<=4;g++){var y=p.t+ph*g/4;ctx.beginPath();ctx.moveTo(p.l,y);ctx.lineTo(p.l+pw,y);ctx.stroke();}
  // labels
  ctx.fillStyle='#8b949e'; ctx.font='11px monospace';
  [0,10,20,30,40].forEach(function(q){ctx.fillText('Q'+q,p.l+q*bw-8,H-p.b+14);});
  ctx.fillText('Mean Phred per read',p.l+pw/2-60,H-2);
}

/* ---------- Init ---------- */
window.addEventListener('load', function(){
  RD.forEach(function(f,i){
    drawQual('qual-'+i, f.qual);
    drawLen('len-'+i, f.len);
    drawKmer('kmer-'+i, f.kmers);
    drawBaseComp('basecomp-'+i, f.basecomp);
    drawQualDist('qualdist-'+i, f.qualdist);
    if (f.tiles && f.tiles.length) drawTile('tile-'+i, f.tiles);
  });
});
"#;

// ---------------------------------------------------------------------------
// HTML generation
// ---------------------------------------------------------------------------

pub fn export_html(state: &SharedState, output_dir: &str) -> io::Result<String> {
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
        file_sections.push_str(&build_file_section(f, i));
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
    html.push_str("<h1>BioFastq-A Quality Analysis Report</h1>\n");
    html.push_str("<div class=\"header-meta\">\n");
    html.push_str(&format!("<span>Generated: {}</span>\n", timestamp));
    html.push_str("<span>Tool version: 2.0.0</span>\n");
    html.push_str(&format!(
        "<span>Files analysed: {}</span>\n",
        files.len()
    ));
    html.push_str(&format!(
        "<span>Processing time: {:.1}s</span>\n",
        elapsed
    ));
    html.push_str("</div>\n</div>\n");

    // Summary table
    if files.len() > 1 {
        html.push_str("<h2>Summary</h2>\n");
        html.push_str(&summary_table);
    }

    // Per-file sections
    html.push_str(&file_sections);

    // Footer
    html.push_str("<div class=\"footer\">Generated by <strong>BioFastq-A v2.0.0</strong></div>\n");

    // Embedded script
    html.push_str("<script>\n");
    html.push_str(&js_with_data);
    html.push_str("\n</script>\n</body>\n</html>\n");

    let stem = report_stem(&files);
    let path = format!("{}/{}_report.html", output_dir, stem);
    fs::write(&path, &html)?;
    Ok(path)
}

fn build_file_section(f: &FileStats, idx: usize) -> String {
    let (n50, n90) = f.compute_n50_n90();
    let fname = Path::new(&f.file_path)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(&f.file_path);

    let qual_class = quality_class(f.avg_quality());
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
    stat_card(&mut s, "N50", &format!("{} bp", n50), "");
    stat_card(&mut s, "N90", &format!("{} bp", n90), "");
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
    stat_card(
        &mut s,
        "Adapter Content",
        &format!("{:.2}%", f.adapter_pct()),
        adapter_class,
    );
    s.push_str("</div>\n"); // stats-grid

    // Charts row 1: quality per position + length distribution
    s.push_str("<div class=\"charts-row\">\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<h3>Per-Base Sequence Quality</h3>\n");
    s.push_str(&format!(
        "<canvas id=\"qual-{}\" width=\"600\" height=\"240\"></canvas>\n",
        idx
    ));
    s.push_str("<p class=\"chart-note\">Green zone: Q&ge;28 &nbsp;|&nbsp; Orange: Q20&ndash;28 &nbsp;|&nbsp; Red: Q&lt;20</p>\n");
    s.push_str("</div>\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<h3>Read Length Distribution</h3>\n");
    s.push_str(&format!(
        "<canvas id=\"len-{}\" width=\"600\" height=\"240\"></canvas>\n",
        idx
    ));
    s.push_str("</div>\n");

    s.push_str("</div>\n"); // charts-row

    // Charts row 2: k-mer + adapter info
    s.push_str("<div class=\"charts-row\">\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<h3>Top K-mer Frequencies (4-mer)</h3>\n");
    s.push_str(&format!(
        "<canvas id=\"kmer-{}\" width=\"600\" height=\"300\"></canvas>\n",
        idx
    ));
    s.push_str("</div>\n");

    // Adapter info card
    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<h3>Adapter Contamination</h3>\n");
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
    s.push_str("</div>\n");

    s.push_str("</div>\n"); // charts-row

    // --- Row 3: Duplication + Per-tile + Trim info ---
    s.push_str("<div class=\"charts-row\">\n");

    // Duplication card
    let dup_class = dup_status_class(f.dup_rate_pct);
    let dup_bar_color = match dup_class {
        "pass" => "var(--green)",
        "warn" => "var(--yellow)",
        _ => "var(--red)",
    };
    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<h3>Duplication Rate (estimate)</h3>\n");
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
    s.push_str("<p class=\"chart-note\">Estimated from fingerprint hashing of first 200,000 reads. &lt;5% = pass &nbsp;|&nbsp; 5–20% = warn &nbsp;|&nbsp; &gt;20% = high</p>\n");
    s.push_str("</div>\n");

    // Per-tile quality card (only if Illumina data was detected)
    let tile_data = f.sorted_tile_quality();
    if !tile_data.is_empty() {
        s.push_str("<div class=\"chart-box\">\n");
        s.push_str("<h3>Per-Tile Quality Score</h3>\n");
        s.push_str(&format!(
            "<canvas id=\"tile-{}\" width=\"600\" height=\"200\"></canvas>\n",
            idx
        ));
        s.push_str(&format!(
            "<p class=\"chart-note\">{} Illumina tiles detected. Bars coloured red→green by quality.</p>\n",
            tile_data.len()
        ));
        s.push_str("</div>\n");
    } else {
        // No Illumina tiles — show note
        s.push_str("<div class=\"chart-box\">\n");
        s.push_str("<h3>Per-Tile Quality Score</h3>\n");
        s.push_str("<p class=\"chart-note\" style=\"padding:20px 0\">No Illumina tile information found in read headers.<br>Per-tile QC requires standard Illumina CASAVA 1.8+ headers.</p>\n");
        s.push_str("</div>\n");
    }

    s.push_str("</div>\n"); // charts-row (row 3)

    // Row 4: Base composition + Quality distribution
    s.push_str("<div class=\"charts-row\">\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<h3>Per-Base Sequence Content</h3>\n");
    s.push_str(&format!(
        "<canvas id=\"basecomp-{}\" width=\"600\" height=\"280\"></canvas>\n",
        idx
    ));
    s.push_str("<p class=\"chart-note\">Base composition (A/C/G/T/N) per read position. Uniform lines indicate balanced sequence.</p>\n");
    s.push_str("</div>\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<h3>Per-Read Quality Distribution</h3>\n");
    s.push_str(&format!(
        "<canvas id=\"qualdist-{}\" width=\"600\" height=\"280\"></canvas>\n",
        idx
    ));
    s.push_str("<p class=\"chart-note\">Distribution of mean Phred quality scores across all reads. Red=poor, orange=acceptable, green=good.</p>\n");
    s.push_str("</div>\n");

    s.push_str("</div>\n"); // charts-row (row 4)

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

    serde_json::json!({
        "qual": qual_per_pos,
        "len": len_dist,
        "kmers": kmer_json,
        "tiles": tile_json,
        "basecomp": basecomp_json,
        "qualdist": qualdist_json,
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

fn gc_class(gc: f64) -> &'static str {
    if gc >= 35.0 && gc <= 65.0 {
        "pass"
    } else if gc >= 25.0 && gc <= 75.0 {
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
