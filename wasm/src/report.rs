use crate::types::{format_bases, format_number, FileStats, ADAPTERS};

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
.summary-table{width:100%;border-collapse:collapse;margin-bottom:28px;font-size:13px;}
.summary-table th{background:var(--card2);color:var(--muted);font-weight:700;text-align:left;padding:9px 12px;border:1px solid var(--border);font-size:11px;text-transform:uppercase;letter-spacing:.05em;}
.summary-table td{padding:9px 12px;border:1px solid var(--border);background:var(--card);}
.summary-table tr:hover td{background:#f0f4f8;}
.file-section{background:var(--card);border:1px solid var(--border);border-radius:10px;margin-bottom:32px;overflow:hidden;box-shadow:0 1px 4px rgba(0,0,0,.06);}
.file-section-header{background:var(--card2);padding:14px 20px;border-bottom:1px solid var(--border);}
.file-section-header .filename{font-size:15px;font-weight:700;font-family:monospace;color:var(--text);}
.file-section-header .filepath{font-size:11px;color:var(--muted);margin-top:2px;font-family:monospace;}
.file-body{padding:20px;}
.stats-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(155px,1fr));gap:10px;margin-bottom:20px;}
.stat-card{background:var(--card2);border:1px solid var(--border);border-radius:8px;padding:12px 14px;}
.stat-label{font-size:10px;font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:.07em;margin-bottom:5px;}
.stat-value{font-size:19px;font-weight:700;font-family:monospace;color:var(--text);}
.stat-value.pass{color:var(--green);}
.stat-value.warn{color:var(--yellow);}
.stat-value.fail{color:var(--red);}
.charts-row{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:16px;margin-bottom:16px;}
.chart-box{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:16px;}
.chart-box canvas{max-width:100%;display:block;}
.chart-header{display:flex;justify-content:space-between;align-items:center;margin-bottom:10px;}
.chart-header h3{margin-bottom:0;}
.chart-note{font-size:11px;color:var(--muted);margin-top:8px;line-height:1.5;}
.dl-btn{display:inline-flex;align-items:center;gap:3px;padding:3px 8px;border-radius:4px;border:1px solid var(--border);background:var(--card2);color:var(--muted);font-size:10px;font-weight:600;cursor:pointer;white-space:nowrap;}
.dl-btn:hover{background:var(--border);color:var(--text);}
.adapter-grid{display:grid;grid-template-columns:auto 1fr;gap:24px;align-items:start;}
.big-pct{font-size:38px;font-weight:700;font-family:monospace;line-height:1;}
.adapter-list{margin-top:8px;padding:0;}
.adapter-list li{list-style:none;font-size:12px;font-family:monospace;color:var(--muted);padding:2px 0;}
.adapter-list li::before{content:"› ";color:var(--accent);}
.dup-gauge-wrap{display:flex;align-items:center;gap:16px;}
.dup-gauge-track{flex:1;height:12px;background:var(--card2);border:1px solid var(--border);border-radius:6px;overflow:hidden;}
.dup-gauge-fill{height:100%;border-radius:6px;}
.trim-box{background:var(--card2);border:1px solid var(--border);border-radius:8px;padding:16px;}
.trim-box p{font-size:13px;color:var(--muted);margin:4px 0;}
.trim-box code{color:var(--accent);font-family:monospace;font-size:12px;}
.footer{text-align:center;color:var(--muted);font-size:12px;margin-top:40px;padding-top:16px;border-top:1px solid var(--border);}
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
.overrep-table{width:100%;border-collapse:collapse;font-size:12px;margin-top:8px;}
.overrep-table th{background:var(--card2);color:var(--muted);padding:7px 10px;text-align:left;border:1px solid var(--border);font-weight:700;font-size:11px;text-transform:uppercase;letter-spacing:.04em;}
.overrep-table td{padding:6px 10px;border:1px solid var(--border);}
.overrep-table tr:hover td{background:#f0f4f8;}
.overrep-seq{word-break:break-all;max-width:320px;font-family:monospace;color:var(--accent);font-size:11px;}
.wasm-badge{display:inline-flex;align-items:center;gap:6px;background:#e0edff;border:1px solid #b3d0f5;border-radius:6px;padding:4px 10px;font-size:11px;font-weight:600;color:#1a5fa8;margin-left:12px;}
"#;

const JS: &str = r#"
'use strict';
var RD = __REPORT_DATA__;

function downloadChart(id, name) {
  var cv = document.getElementById(id); if (!cv) return;
  var a = document.createElement('a');
  a.download = (name || id) + '.png';
  a.href = cv.toDataURL('image/png');
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
}

function printReport() { window.print(); }

function drawQual(id, data, pct) {
  var cv = document.getElementById(id);
  if (!cv || !data || !data.length) return;
  var ctx = cv.getContext('2d');
  var W = cv.clientWidth || cv.width; cv.width = W;
  var H = cv.height;
  var pad = {t:18,r:16,b:36,l:44};
  var pw = W-pad.l-pad.r, ph = H-pad.t-pad.b, MQ = 42;
  ctx.fillStyle = '#ffffff'; ctx.fillRect(0,0,W,H);
  [[0,20,'rgba(220,53,69,.07)'],[20,28,'rgba(255,193,7,.08)'],[28,MQ,'rgba(25,135,84,.06)']].forEach(function(z){
    var y1=pad.t+(1-z[1]/MQ)*ph, y2=pad.t+(1-z[0]/MQ)*ph;
    ctx.fillStyle=z[2]; ctx.fillRect(pad.l,y1,pw,y2-y1);
  });
  [[20,'rgba(220,53,69,.35)'],[28,'rgba(200,150,0,.45)'],[30,'rgba(25,135,84,.35)']].forEach(function(qc){
    var y=pad.t+(1-qc[0]/MQ)*ph;
    ctx.setLineDash([3,4]); ctx.strokeStyle=qc[1]; ctx.lineWidth=1;
    ctx.beginPath(); ctx.moveTo(pad.l,y); ctx.lineTo(W-pad.r,y); ctx.stroke();
  });
  ctx.setLineDash([]);
  [0,10,20,30,40].forEach(function(q){
    var y=pad.t+(1-q/MQ)*ph;
    ctx.strokeStyle='rgba(0,0,0,.07)'; ctx.lineWidth=1;
    ctx.beginPath(); ctx.moveTo(pad.l,y); ctx.lineTo(W-pad.r,y); ctx.stroke();
  });
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
  ctx.strokeStyle='#555'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();
  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='right';
  [0,10,20,28,30,40].forEach(function(q){
    ctx.fillText(q, pad.l-5, pad.t+(1-q/MQ)*ph+4);
  });
  ctx.textAlign='center';
  var n = data.length;
  var step = n > 300 ? 100 : n > 150 ? 50 : n > 50 ? 25 : 10;
  for (var i=0; i<n; i+=step) {
    ctx.fillText(i+1, pad.l+(i/Math.max(n-1,1))*pw, H-pad.b+14);
  }
  ctx.fillText(n, pad.l+pw, H-pad.b+14);
  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='center';
  ctx.fillText('Position (bp)', pad.l+pw/2, H-2);
  ctx.save(); ctx.translate(10,pad.t+ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('Phred Quality', 0, 0); ctx.restore();
  if (pct && pct.length) {
    var lx = W-pad.r-126, ly = pad.t+4;
    ctx.font='10px monospace'; ctx.textAlign='left';
    ctx.fillStyle='#1a5fa8'; ctx.fillRect(lx,ly,14,3);
    ctx.fillStyle='#5a6370'; ctx.fillText('Mean',lx+17,ly+5);
    ctx.fillStyle='rgba(200,120,0,.7)'; ctx.fillRect(lx,ly+11,14,2);
    ctx.fillStyle='#5a6370'; ctx.fillText('Median (Q50)',lx+17,ly+16);
    ctx.fillStyle='rgba(26,95,168,.18)'; ctx.fillRect(lx,ly+22,14,8);
    ctx.fillStyle='#5a6370'; ctx.fillText('IQR Q25-Q75',lx+17,ly+30);
  }
}

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
  ctx.strokeStyle='#555'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();
  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='center';
  var lblN = Math.min(6, n);
  for (var i=0; i<lblN; i++) {
    var idx = Math.round(i*(n-1)/(Math.max(lblN-1,1)));
    var x = pad.l+(idx+0.5)*bw;
    var len = data[idx][0];
    ctx.fillText(len>999?(len/1000).toFixed(1)+'k':len+'bp', x, H-pad.b+14);
  }
  ctx.textAlign='right';
  [0,.5,1].forEach(function(frac){
    var v=Math.round(maxC*frac);
    var label = v>999999?(v/1000000).toFixed(1)+'M':v>999?(v/1000).toFixed(0)+'k':String(v);
    ctx.fillText(label, pad.l-5, pad.t+(1-frac)*ph+4);
  });
  ctx.textAlign='center';
  ctx.fillText('Read Length', pad.l+pw/2, H-2);
  ctx.save(); ctx.translate(10,pad.t+ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('Count', 0, 0); ctx.restore();
}

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

function drawBaseComp(id, data) {
  var cv=document.getElementById(id); if(!cv)return;
  var ctx=cv.getContext('2d'), W=cv.width, H=cv.height;
  var pad={l:52,r:16,t:20,b:36};
  var pw=W-pad.l-pad.r, ph=H-pad.t-pad.b;
  ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H);
  if(!data||!data.length){ctx.fillStyle='#5a6370';ctx.font='14px monospace';ctx.textAlign='center';ctx.fillText('No data',W/2,H/2);return;}
  var colors=['#1a5fa8','#1a7a3e','#dc3545','#7a5500','#5a6370'];
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

function drawGcDist(id, data, meanGc) {
  var cv=document.getElementById(id); if(!cv||!data||!data.length)return;
  var ctx=cv.getContext('2d'), W=cv.width, H=cv.height;
  var pad={l:52,r:16,t:24,b:36};
  var pw=W-pad.l-pad.r, ph=H-pad.t-pad.b;
  ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H);
  var total=data.reduce(function(s,v){return s+v;},0)||1;
  var max=Math.max.apply(null,data)||1;
  var n=data.length, bw=pw/n;
  data.forEach(function(cnt,gc){
    var frac=cnt/max;
    var barH=ph*frac;
    ctx.fillStyle='rgba(26,95,168,0.55)';
    ctx.fillRect(pad.l+gc*bw, pad.t+ph-barH, Math.max(1,bw-1), barH);
  });
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
  ctx.strokeStyle='#555'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();
  ctx.fillStyle='#5a6370'; ctx.font='11px monospace'; ctx.textAlign='center';
  [0,20,40,50,60,80,100].forEach(function(g){ctx.fillText(g+'%',pad.l+g*bw,H-pad.b+14);});
  ctx.fillText('GC Content (%)', pad.l+pw/2, H-2);
}

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

function drawNContent(id, data) {
  var cv=document.getElementById(id); if(!cv||!data||!data.length)return;
  var ctx=cv.getContext('2d'), W=cv.width, H=cv.height;
  var pad={l:52,r:16,t:20,b:36};
  var pw=W-pad.l-pad.r, ph=H-pad.t-pad.b;
  ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H);
  var maxN=Math.max.apply(null,data)||0.1;
  var cap=Math.max(5,Math.ceil(maxN));
  ctx.fillStyle='rgba(220,53,69,.06)'; ctx.fillRect(pad.l,pad.t,pw,ph);
  ctx.beginPath();
  data.forEach(function(v,i){
    var x=pad.l+(i/(data.length-1||1))*pw;
    var y=pad.t+ph*(1-v/cap);
    if(i===0)ctx.moveTo(x,y); else ctx.lineTo(x,y);
  });
  ctx.strokeStyle='#dc3545'; ctx.lineWidth=1.5; ctx.stroke();
  if(cap>=5){
    var ty=pad.t+ph*(1-5/cap);
    ctx.setLineDash([3,4]); ctx.strokeStyle='rgba(220,53,69,.4)'; ctx.lineWidth=1;
    ctx.beginPath(); ctx.moveTo(pad.l,ty); ctx.lineTo(W-pad.r,ty); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle='rgba(220,53,69,.6)'; ctx.font='10px monospace'; ctx.textAlign='right';
    ctx.fillText('5%',pad.l-4,ty+4);
  }
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
  ctx.strokeStyle = 'rgba(0,0,0,.07)'; ctx.lineWidth = 1;
  [0,1,2,3,4].forEach(function(gi) {
    var y = pad.t + gi/4 * ph; ctx.beginPath(); ctx.moveTo(pad.l,y); ctx.lineTo(W-pad.r,y); ctx.stroke();
    var x = pad.l + gi/4 * pw; ctx.beginPath(); ctx.moveTo(x,pad.t); ctx.lineTo(x,H-pad.b); ctx.stroke();
  });
  data.forEach(function(d) {
    var x = pad.l + (d[0] / maxLen) * pw;
    var y = pad.t + (1 - (d[1] - minQ) / qRange) * ph;
    ctx.beginPath(); ctx.arc(x, y, 2.5, 0, 2*Math.PI);
    ctx.fillStyle = 'rgba(26,95,168,0.55)'; ctx.fill();
  });
  ctx.strokeStyle = '#555'; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b); ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();
  ctx.fillStyle = '#5a6370'; ctx.font = '11px monospace'; ctx.textAlign = 'right';
  [minQ, Math.round((minQ+maxQ)/2), maxQ].forEach(function(q) {
    var y = pad.t + (1 - (q - minQ) / qRange) * ph; ctx.fillText(q, pad.l-4, y+4);
  });
  ctx.textAlign = 'center';
  [0, 0.25, 0.5, 0.75, 1.0].forEach(function(frac) {
    var v = Math.round(maxLen * frac);
    var lbl = v >= 1000 ? (v/1000).toFixed(1)+'k' : v+'bp';
    ctx.fillText(lbl, pad.l + frac*pw, H-pad.b+14);
  });
  ctx.fillText('Read Length (bp)', pad.l+pw/2, H-2);
  ctx.save(); ctx.translate(12, pad.t+ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('Mean Phred Quality', 0, 0); ctx.restore();
}

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
  });
});
"#;

fn escape_html(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn stat_card(buf: &mut String, label: &str, value: &str, cls: &str) {
    buf.push_str("<div class=\"stat-card\">\n");
    buf.push_str(&format!("<div class=\"stat-label\">{}</div>\n", escape_html(label)));
    if cls.is_empty() {
        buf.push_str(&format!("<div class=\"stat-value\">{}</div>\n", value));
    } else {
        buf.push_str(&format!("<div class=\"stat-value {}\">{}</div>\n", cls, value));
    }
    buf.push_str("</div>\n");
}

fn quality_class(q: f64, long_read: bool) -> &'static str {
    if long_read {
        if q >= 10.0 { "pass" } else if q >= 8.0 { "warn" } else { "fail" }
    } else if q >= 30.0 { "pass" } else if q >= 20.0 { "warn" } else { "fail" }
}

fn gc_class(gc: f64) -> &'static str {
    if (35.0..=65.0).contains(&gc) { "pass" }
    else if (25.0..=75.0).contains(&gc) { "warn" }
    else { "fail" }
}

fn dup_class(pct: f64) -> &'static str {
    if pct < 5.0 { "pass" } else if pct < 20.0 { "warn" } else { "fail" }
}

fn adapter_class(pct: f64) -> &'static str {
    if pct < 1.0 { "pass" } else if pct < 10.0 { "warn" } else { "fail" }
}

fn build_file_json(f: &FileStats) -> serde_json::Value {
    let qual_per_pos: Vec<f64> = f.avg_qual_per_position();
    let len_dist: Vec<[u64; 2]> = f.sorted_length_dist().iter().map(|(l, c)| [*l, *c]).collect();
    let kmer_json: Vec<serde_json::Value> = f.top_kmers(20).into_iter()
        .map(|(k, v)| serde_json::json!([k, v])).collect();
    let tile_json: Vec<serde_json::Value> = {
        let mut tv: Vec<(u32, f64)> = f.per_tile_quality.iter()
            .map(|(&t, &(s, c))| (t, if c == 0 { 0.0 } else { s as f64 / c as f64 })).collect();
        tv.sort_by_key(|(t, _)| *t);
        tv.into_iter().map(|(t, q)| serde_json::json!([t, q])).collect()
    };
    let basecomp_json: Vec<serde_json::Value> = f.base_composition_pct().into_iter()
        .map(|pct| serde_json::json!(pct)).collect();
    let qualdist_json: Vec<u64> = f.quality_distribution.clone();
    let qual_pct_json: Vec<serde_json::Value> = f.qual_percentiles_per_position().into_iter()
        .map(|pct| serde_json::json!(pct)).collect();
    let gc_dist_json: Vec<u64> = f.gc_distribution.clone();
    let dup_hist_labels = ["1x","2x","3-4x","5-9x","10-49x","50-99x","100-499x","500-999x",">=1000x"];
    let dup_hist_json: Vec<serde_json::Value> = f.dup_level_histogram.iter()
        .enumerate().map(|(i, &p)| serde_json::json!([dup_hist_labels[i], p])).collect();
    let n_content_json: Vec<f64> = f.base_composition_pct().iter().map(|p| p[4]).collect();
    let qvl_json: Vec<serde_json::Value> = f.qual_vs_length.iter()
        .map(|(l, q)| serde_json::json!([l, q])).collect();

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
        "long_read": f.long_read_mode,
    })
}

fn build_file_section(f: &FileStats, idx: usize) -> String {
    let (n50, n90) = f.compute_n50_n90();
    let fname = f.file_path.rsplit('/').next().unwrap_or(&f.file_path);
    let lr = f.long_read_mode;
    let qual_cls = quality_class(f.avg_quality(), lr);
    let adp_cls = adapter_class(f.adapter_pct());

    let adapter_list: String = ADAPTERS.iter()
        .map(|(name, _)| format!("<li>{}</li>", name))
        .collect();

    let mut s = String::new();
    s.push_str("<div class=\"file-section\">\n");
    s.push_str("<div class=\"file-section-header\">\n");
    s.push_str(&format!("<div class=\"filename\">{}</div>\n", escape_html(fname)));
    s.push_str(&format!("<div class=\"filepath\">{}</div>\n", escape_html(&f.file_path)));
    s.push_str("</div>\n<div class=\"file-body\">\n");

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

    s.push_str("<div class=\"stats-grid\">\n");
    stat_card(&mut s, "Total Reads", &format_number(f.read_count), "");
    stat_card(&mut s, "Total Bases", &format_bases(f.total_bases), "");
    stat_card(&mut s, "Avg Read Length", &format!("{:.1} bp", f.avg_length()), "");
    stat_card(&mut s, "Min / Max Length",
        &format!("{} / {} bp", f.effective_min_length(), f.max_length), "");
    if lr {
        stat_card(&mut s, "N50", &format!("{} bp", n50), if n50 > 10000 { "pass" } else { "warn" });
        stat_card(&mut s, "N90", &format!("{} bp", n90), "");
    } else {
        stat_card(&mut s, "N50", &format!("{} bp", n50), "");
        stat_card(&mut s, "N90", &format!("{} bp", n90), "");
    }
    stat_card(&mut s, "GC Content", &format!("{:.2}%", f.gc_content()), gc_class(f.gc_content()));
    stat_card(&mut s, "Avg Quality", &format!("Q{:.1}", f.avg_quality()), qual_cls);
    if lr {
        stat_card(&mut s, "≥Q20 Reads", &format!("{:.1}%", f.q20_pct()),
            if f.q20_pct() >= 30.0 { "pass" } else { "warn" });
    } else {
        stat_card(&mut s, "≥Q20 Reads", &format!("{:.1}%", f.q20_pct()),
            if f.q20_pct() >= 80.0 { "pass" } else { "warn" });
        stat_card(&mut s, "≥Q30 Reads", &format!("{:.1}%", f.q30_pct()),
            if f.q30_pct() >= 70.0 { "pass" } else { "warn" });
    }
    stat_card(&mut s, "Adapter Content", &format!("{:.2}%", f.adapter_pct()), adp_cls);
    s.push_str("</div>\n");

    // Row 1: Quality per position + Length distribution
    s.push_str("<div class=\"charts-row\">\n");
    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Per-Base Sequence Quality</h3>\n");
    s.push_str(&format!("<button class=\"dl-btn no-print\" onclick=\"downloadChart('qual-{}','per_base_quality')\">&#8595; PNG</button>\n", idx));
    s.push_str("</div>\n");
    s.push_str(&format!("<canvas id=\"qual-{}\" width=\"600\" height=\"240\"></canvas>\n", idx));
    s.push_str("<p class=\"chart-note\">Blue = mean quality. Orange dashed = median (Q50). Shaded band = IQR (Q25&ndash;Q75). Green: Q&ge;28 | Yellow: Q20&ndash;28 | Red: Q&lt;20.</p>\n");
    s.push_str("</div>\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Read Length Distribution</h3>\n");
    s.push_str(&format!("<button class=\"dl-btn no-print\" onclick=\"downloadChart('len-{}','read_length_dist')\">&#8595; PNG</button>\n", idx));
    s.push_str("</div>\n");
    s.push_str(&format!("<canvas id=\"len-{}\" width=\"600\" height=\"240\"></canvas>\n", idx));
    s.push_str("</div>\n");
    s.push_str("</div>\n");

    // Row 2: K-mer + Adapter
    s.push_str("<div class=\"charts-row\">\n");
    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Top K-mer Frequencies (4-mer)</h3>\n");
    s.push_str(&format!("<button class=\"dl-btn no-print\" onclick=\"downloadChart('kmer-{}','kmer_frequencies')\">&#8595; PNG</button>\n", idx));
    s.push_str("</div>\n");
    s.push_str(&format!("<canvas id=\"kmer-{}\" width=\"600\" height=\"300\"></canvas>\n", idx));
    s.push_str("</div>\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<h3>Adapter Contamination</h3>\n");
    s.push_str("<div class=\"adapter-grid\">\n");
    s.push_str(&format!(
        "<div><div class=\"big-pct {}\">{:.1}%</div><div style=\"color:var(--muted);font-size:12px;margin-top:4px\">reads with adapter</div></div>\n",
        adp_cls, f.adapter_pct()
    ));
    s.push_str("<div>\n");
    s.push_str("<p style=\"font-size:12px;color:var(--muted);margin-bottom:6px\">Screened sequences:</p>\n");
    s.push_str("<ul class=\"adapter-list\">\n");
    s.push_str(&adapter_list);
    s.push_str("</ul>\n</div>\n</div>\n");
    s.push_str("</div>\n");
    s.push_str("</div>\n");

    // Row 3: Duplication + Per-tile
    s.push_str("<div class=\"charts-row\">\n");
    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<h3>Duplication Rate (estimate)</h3>\n");
    let dc = dup_class(f.dup_rate_pct);
    let dup_bar_color = match dc {
        "pass" => "var(--green)",
        "warn" => "var(--yellow)",
        _ => "var(--red)",
    };
    s.push_str("<div class=\"dup-gauge-wrap\" style=\"margin-bottom:12px\">\n");
    s.push_str(&format!(
        "<span class=\"stat-value {}\" style=\"font-size:32px;min-width:80px\">{:.1}%</span>\n",
        dc, f.dup_rate_pct
    ));
    s.push_str("<div class=\"dup-gauge-track\">\n");
    s.push_str(&format!(
        "<div class=\"dup-gauge-fill\" style=\"width:{:.1}%;background:{}\"></div>\n",
        f.dup_rate_pct.min(100.0), dup_bar_color
    ));
    s.push_str("</div>\n</div>\n");
    s.push_str("<p class=\"chart-note\">Sort-dedup estimate from first 200,000 read fingerprints. &lt;5% = pass | 5&ndash;20% = warn | &gt;20% = high</p>\n");
    s.push_str("</div>\n");

    if !lr && !f.per_tile_quality.is_empty() {
        s.push_str("<div class=\"chart-box\">\n");
        s.push_str("<div class=\"chart-header\">\n");
        s.push_str("<h3>Per-Tile Quality Score</h3>\n");
        s.push_str(&format!("<button class=\"dl-btn no-print\" onclick=\"downloadChart('tile-{}','per_tile_quality')\">&#8595; PNG</button>\n", idx));
        s.push_str("</div>\n");
        s.push_str(&format!("<canvas id=\"tile-{}\" width=\"600\" height=\"200\"></canvas>\n", idx));
        s.push_str("</div>\n");
    }
    s.push_str("</div>\n");

    // Row 4: Base composition + Quality distribution
    s.push_str("<div class=\"charts-row\">\n");
    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Per-Base Sequence Content</h3>\n");
    s.push_str(&format!("<button class=\"dl-btn no-print\" onclick=\"downloadChart('basecomp-{}','base_composition')\">&#8595; PNG</button>\n", idx));
    s.push_str("</div>\n");
    s.push_str(&format!("<canvas id=\"basecomp-{}\" width=\"600\" height=\"280\"></canvas>\n", idx));
    s.push_str("<p class=\"chart-note\">Nucleotide composition (A/C/G/T/N) per read position. Parallel flat lines = unbiased.</p>\n");
    s.push_str("</div>\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Per-Read Quality Distribution</h3>\n");
    s.push_str(&format!("<button class=\"dl-btn no-print\" onclick=\"downloadChart('qualdist-{}','quality_distribution')\">&#8595; PNG</button>\n", idx));
    s.push_str("</div>\n");
    s.push_str(&format!("<canvas id=\"qualdist-{}\" width=\"600\" height=\"280\"></canvas>\n", idx));
    s.push_str("<p class=\"chart-note\">Histogram of mean Phred quality scores across all reads. Red = Q&lt;20, yellow = Q20&ndash;30, green = Q&ge;30.</p>\n");
    s.push_str("</div>\n");
    s.push_str("</div>\n");

    // Row 5: GC distribution + Dup histogram + N content
    s.push_str("<div class=\"charts-row\">\n");
    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>GC Content Distribution</h3>\n");
    s.push_str(&format!("<button class=\"dl-btn no-print\" onclick=\"downloadChart('gcdist-{}','gc_distribution')\">&#8595; PNG</button>\n", idx));
    s.push_str("</div>\n");
    s.push_str(&format!("<canvas id=\"gcdist-{}\" width=\"600\" height=\"240\"></canvas>\n", idx));
    s.push_str("<p class=\"chart-note\">Per-read GC content histogram. Red dashed = theoretical normal. Deviation may indicate contamination.</p>\n");
    s.push_str("</div>\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Duplication Level</h3>\n");
    s.push_str(&format!("<button class=\"dl-btn no-print\" onclick=\"downloadChart('duphist-{}','dup_level_histogram')\">&#8595; PNG</button>\n", idx));
    s.push_str("</div>\n");
    s.push_str(&format!("<canvas id=\"duphist-{}\" width=\"600\" height=\"240\"></canvas>\n", idx));
    s.push_str("<p class=\"chart-note\">Percentage of reads at each duplication level. High 2x+ bins indicate over-amplification.</p>\n");
    s.push_str("</div>\n");

    s.push_str("<div class=\"chart-box\">\n");
    s.push_str("<div class=\"chart-header\">\n");
    s.push_str("<h3>Per-Base N Content</h3>\n");
    s.push_str(&format!("<button class=\"dl-btn no-print\" onclick=\"downloadChart('ncontent-{}','n_content')\">&#8595; PNG</button>\n", idx));
    s.push_str("</div>\n");
    s.push_str(&format!("<canvas id=\"ncontent-{}\" width=\"600\" height=\"240\"></canvas>\n", idx));
    s.push_str("<p class=\"chart-note\">Percentage of uncalled (N) bases per position. Red dashed = 5% warning threshold.</p>\n");
    s.push_str("</div>\n");
    s.push_str("</div>\n");

    // Long-read: Quality vs Length scatter
    if lr && !f.qual_vs_length.is_empty() {
        s.push_str("<div class=\"charts-row\">\n");
        s.push_str("<div class=\"chart-box\">\n");
        s.push_str("<div class=\"chart-header\">\n");
        s.push_str("<h3>Quality vs Read Length</h3>\n");
        s.push_str(&format!("<button class=\"dl-btn no-print\" onclick=\"downloadChart('qvl-{}','qual_vs_length')\">&#8595; PNG</button>\n", idx));
        s.push_str("</div>\n");
        s.push_str(&format!("<canvas id=\"qvl-{}\" width=\"600\" height=\"280\"></canvas>\n", idx));
        s.push_str("<p class=\"chart-note\">Sampled read length vs mean Phred quality (up to 2000 reads).</p>\n");
        s.push_str("</div>\n");
        s.push_str("</div>\n");
    }

    // Overrepresented sequences
    s.push_str("<div class=\"chart-box\" style=\"margin-top:16px\">\n");
    s.push_str("<h3>Overrepresented Sequences</h3>\n");
    if f.overrepresented_sequences.is_empty() {
        s.push_str("<p class=\"chart-note\" style=\"padding:12px 0\">No overrepresented sequences found (&ge;0.1% of sampled reads).</p>\n");
    } else {
        s.push_str("<p class=\"chart-note\" style=\"margin-bottom:8px\">Sampled first 200,000 reads. Sequences present in &ge;0.1% of reads:</p>\n");
        s.push_str("<table class=\"overrep-table\">\n");
        s.push_str("<thead><tr><th>#</th><th>Sequence (first 50 bp)</th><th>Count</th><th>%</th><th>Possible Source</th></tr></thead>\n<tbody>\n");
        for (i, seq) in f.overrepresented_sequences.iter().enumerate() {
            let src_cls = if seq.possible_source == "No hit" { "color:var(--muted)" } else { "color:var(--yellow)" };
            s.push_str(&format!(
                "<tr><td>{}</td><td class=\"overrep-seq\">{}</td><td>{}</td><td>{:.2}%</td><td style=\"{}\">{}</td></tr>\n",
                i + 1, escape_html(&seq.sequence), format_number(seq.count),
                seq.percentage, src_cls, escape_html(&seq.possible_source),
            ));
        }
        s.push_str("</tbody></table>\n");
    }
    s.push_str("</div>\n");

    s.push_str("</div>\n</div>\n");
    s
}

/// Generate a complete HTML report for one FASTQ file.
/// `timestamp` is a pre-formatted string (e.g. from JS `new Date().toLocaleString()`).
pub fn generate_html(f: &FileStats, timestamp: &str, elapsed_ms: f64) -> String {
    let file_section = build_file_section(f, 0);
    let json_data = serde_json::to_string(&[build_file_json(f)]).unwrap_or_else(|_| "[]".into());
    let js_with_data = JS.replace("__REPORT_DATA__", &json_data);

    let fname = f.file_path.rsplit('/').next().unwrap_or(&f.file_path);
    let elapsed_s = elapsed_ms / 1000.0;

    let mut html = String::with_capacity(512 * 1024);
    html.push_str("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n");
    html.push_str("<meta charset=\"UTF-8\">\n");
    html.push_str("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
    html.push_str(&format!("<title>BioFastq-A Report — {}</title>\n", escape_html(fname)));
    html.push_str("<style>\n");
    html.push_str(CSS);
    html.push_str("</style>\n</head>\n<body>\n");

    html.push_str("<div class=\"header\">\n");
    html.push_str("<div class=\"header-left\">\n");
    html.push_str("<h1>BioFastq-A <span>Quality Analysis Report</span> <span class=\"wasm-badge\">&#9654; WebAssembly</span></h1>\n");
    html.push_str("<div class=\"header-meta\">\n");
    html.push_str(&format!("<span>Generated: {}</span>\n", escape_html(timestamp)));
    html.push_str("<span>Tool version: 2.3.1 (Wasm)</span>\n");
    html.push_str(&format!("<span>Analysis time: {:.1}s</span>\n", elapsed_s));
    html.push_str("<span>Data never left your browser</span>\n");
    html.push_str("</div>\n</div>\n");
    html.push_str("<div class=\"header-actions no-print\">\n");
    html.push_str("<button class=\"btn\" onclick=\"printReport()\">&#128438; Print / PDF</button>\n");
    html.push_str("<button class=\"btn\" onclick=\"saveReport()\">&#8595; Save HTML</button>\n");
    html.push_str("</div>\n");
    html.push_str("</div>\n");

    html.push_str(&file_section);

    html.push_str("<div class=\"footer\">Generated by <strong>BioFastq-A v2.3.1 (WebAssembly)</strong> &mdash; your sequencing data never left your browser.</div>\n");

    html.push_str("<script>\n");
    html.push_str(&js_with_data);
    html.push_str("\nfunction saveReport(){var a=document.createElement('a');a.download='biofastq_report.html';a.href='data:text/html;charset=utf-8,'+encodeURIComponent(document.documentElement.outerHTML);a.click();}\n");
    html.push_str("</script>\n</body>\n</html>\n");

    html
}
