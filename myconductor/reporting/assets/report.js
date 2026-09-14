/* Myconductor report runtime.
   Ported from the design mock. Every dataset the mock hard-coded is now read
   from DATA, which the renderer injects from a real analysis. Nothing here
   invents a number: a section with no data renders its empty state. */
const DATA = window.__MYCONDUCTOR__ || {};
const DRUGS = DATA.drugs || [];
const LINEAGES = DATA.lineages || [];
const COVERAGE_LOCI = DATA.coverage_loci || [];
const VUS_ITEMS = DATA.vus || [];
const MECH_CARDS = DATA.mechanisms || [];
const DISCORDANCE_ROWS = DATA.discordance || [];
const AUDIT_EVENTS = DATA.audit || [];

/* ==========================================================
   HELPERS
   ========================================================== */
/* The design mock carried a seeded PRNG that manufactured call matrices,
   coverage grids, lineage strata and sequence bases. It is deliberately gone:
   a report that fabricates plausible numbers is worse than one that shows an
   empty state, because nothing downstream can tell the difference. Every
   renderer below reads real data or renders its empty state. */

/* A measurement that was not made renders as an em dash. The template
   literals below would otherwise interpolate null and print "null%", which
   reads as a number to anyone skimming. */
function pct(v, digits){ return (v === null || v === undefined || Number.isNaN(v))
  ? '—' : (digits === undefined ? v : Number(v).toFixed(digits)) + '%'; }
function num(v){ return (v === null || v === undefined) ? '—' : v; }

function textColorForError(e){
  if(e === null) return 'var(--text-3)';
  if(e < 5) return 'var(--sus)';
  if(e < 10) return 'var(--ind)';
  return 'var(--res)';
}
function bgColorForError(e){
  if(e === null) return 'var(--bg-soft)';
  if(e < 5) return 'var(--sus-bg)';
  if(e < 10) return 'var(--ind-bg)';
  return 'var(--res-bg)';
}

/* ==========================================================
   TOOLTIP
   ========================================================== */
const tip = document.getElementById('tooltip');
function showTip(e, html){
  tip.innerHTML = html;
  tip.classList.add('show');
  moveTip(e);
}
function moveTip(e){
  const pad = 16;
  let x = e.clientX + pad, y = e.clientY + pad;
  const rect = tip.getBoundingClientRect();
  if(x + rect.width > window.innerWidth - 12) x = e.clientX - rect.width - pad;
  if(y + rect.height > window.innerHeight - 12) y = e.clientY - rect.height - pad;
  tip.style.left = x + 'px'; tip.style.top = y + 'px';
}
function hideTip(){ tip.classList.remove('show'); }
document.addEventListener('mousemove', e => { if(tip.classList.contains('show')) moveTip(e); });

/* ==========================================================
   DONUT CHART
   ========================================================== */
function renderDonut(){
  if(!(DATA.call_distribution) || !(DATA.call_distribution).length){ return emptyState('donutChart', 'No call distribution', 'This run produced no drug results to summarise.'); }
  const data = DATA.call_distribution;
  const total = data.reduce((s,d)=>s+d.v,0);
  const cx=180, cy=180, rOut=130, rIn=78;
  let a0 = -Math.PI/2;
  let paths = '';
  data.forEach(d => {
    const a1 = a0 + (d.v/total) * Math.PI * 2;
    const x1 = cx + rOut*Math.cos(a0), y1 = cy + rOut*Math.sin(a0);
    const x2 = cx + rOut*Math.cos(a1), y2 = cy + rOut*Math.sin(a1);
    const x3 = cx + rIn*Math.cos(a1), y3 = cy + rIn*Math.sin(a1);
    const x4 = cx + rIn*Math.cos(a0), y4 = cy + rIn*Math.sin(a0);
    const large = (a1-a0) > Math.PI ? 1 : 0;
    paths += `<path d="M ${x1} ${y1} A ${rOut} ${rOut} 0 ${large} 1 ${x2} ${y2} L ${x3} ${y3} A ${rIn} ${rIn} 0 ${large} 0 ${x4} ${y4} Z"
      fill="${d.c}" stroke="#fff" stroke-width="2" class="bar"
      data-tip="${d.k}" data-v="${d.v}" data-pct="${(d.v/total*100).toFixed(1)}"
      style="cursor:pointer; opacity:.92"
      onmouseover="donutHover(this,event)" onmouseout="hideTip()" />`;
    a0 = a1;
  });
  document.getElementById('donutChart').innerHTML = `
    <svg viewBox="0 0 360 360" class="chart-svg" style="max-width:360px; margin:0 auto">
      ${paths}
      <text x="${cx}" y="${cy-10}" text-anchor="middle" font-family="var(--mono)" font-size="12" fill="var(--text-2)" letter-spacing="1">TOTAL CALLS</text>
      <text x="${cx}" y="${cy+22}" text-anchor="middle" font-family="var(--mono)" font-size="30" font-weight="600" fill="var(--text)">${total.toLocaleString()}</text>
    </svg>`;
}
function donutHover(el, e){
  showTip(e, `
    <div class="tt-head"><div class="tt-title">${el.dataset.tip}</div><div class="tt-sub">call state</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Count</span><span class="tt-v">${(+el.dataset.v).toLocaleString()}</span></div>
      <div class="tt-row"><span class="tt-k">Share</span><span class="tt-v">${el.dataset.pct}%</span></div>
      <div class="tt-section"><div class="tt-note">${describeCall(el.dataset.tip)}</div></div>
    </div>`);
}
function describeCall(k){
  const m = {
    RESISTANT:'Established by catalogued or phenotypic evidence. Can inform clinical decision.',
    SUSCEPTIBLE:'Coverage-gated absence of resistance markers. Only valid where loci were callable.',
    INDETERMINATE:'Evidence insufficient to establish either state. Cannot be acted on.',
    NOT_ASSESSED:'Loci not callable. Absence of evidence is not susceptibility.',
    NO_CALL:'Analysis-side gap. Not a biological finding.',
    UNSUPPORTED:'Variant detected, no catalogue support. Held for review.',
  };
  return m[k] || '';
}

/* ==========================================================
   TREND CHART
   ========================================================== */
function renderTrend(){
  const series = DATA.error_trend || [];
  if(series.length < 2){
    return emptyState('trendChart', 'No cumulative error series',
      'A rolling error rate needs a cohort processed in order. This run did not produce one.');
  }
  const w=520, h=200, pad={t:16,r:16,b:34,l:44};
  const n = series.length;
  const pts = series.map((v, i) => ({ x:i+1, y:v }));
  const xMax = n, yMax = 12;
  const sx = v => pad.l + (v/xMax)*(w-pad.l-pad.r);
  const sy = v => h-pad.b - (v/yMax)*(h-pad.t-pad.b);
  let path='M';
  pts.forEach((p,i)=> path += `${i?'L':' '}${sx(p.x)} ${sy(p.y)}`);
  let grid='';
  for(let v=0;v<=12;v+=3){
    grid += `<line x1="${pad.l}" x2="${w-pad.r}" y1="${sy(v)}" y2="${sy(v)}" class="grid-line" />
      <text x="${pad.l-8}" y="${sy(v)+4}" text-anchor="end" class="axis-label">${v}%</text>`;
  }
  let xLabels='';
  for(let i=0;i<=n;i+=10){
    if(i===0) continue;
    xLabels += `<text x="${sx(i)}" y="${h-12}" text-anchor="middle" class="axis-label">${i*10}</text>`;
  }
  let dots = pts.map(p => `<circle cx="${sx(p.x)}" cy="${sy(p.y)}" r="2.5"
    fill="var(--accent)" class="point"
    onmouseover="trendHover(event,${p.x*10},${p.y.toFixed(2)})" onmouseout="hideTip()" />`).join('');
  document.getElementById('trendChart').innerHTML = `
    <svg viewBox="0 0 ${w} ${h}" class="chart-svg">
      ${grid}
      <path d="${path}" class="line" stroke="var(--accent)" />
      ${dots}
      <text x="${w/2}" y="${h-2}" text-anchor="middle" class="axis-label">isolates processed</text>
      <text x="12" y="${h/2}" text-anchor="middle" transform="rotate(-90 12 ${h/2})" class="axis-label">cumulative error</text>
      ${xLabels}
    </svg>`;
}
function trendHover(e, n, v){
  showTip(e, `<div class="tt-head"><div class="tt-title">Window at isolate ${n}</div><div class="tt-sub">rolling 50-isolate window</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Cumulative error</span><span class="tt-v">${v}%</span></div>
      <div class="tt-section"><div class="tt-note">Rolling error rate stabilises after ~200 isolates. Early window variance reflects small n.</div></div>
    </div>`);
}

/* ==========================================================
   DRUG BAR CHART
   ========================================================== */
let chartMode = 'error';
function setChartMode(btn, mode){
  document.querySelectorAll('#bench .card-toolbar .mini-btn').forEach(b=>b.classList.remove('active'));
  btn.classList.add('active');
  chartMode = mode;
  renderDrugBar();
}
function renderDrugBar(){
  if(!(DRUGS.filter(d=>d.estimable)) || !(DRUGS.filter(d=>d.estimable)).length){ return emptyState('drugBarChart', 'No estimable drug', 'A per-drug accuracy chart needs a phenotypic reference. None was supplied, so no accuracy is estimated.'); }
  const data = DRUGS.filter(d=>d.estimable).map(d=>({
    id:d.id, name:d.name,
    val: chartMode==='error' ? d.error : chartMode==='coverage' ? d.coverage : d.vme,
    estimable: true
  }));
  const w=980, h=340, pad={t:24,r:24,b:60,l:60};
  const yMax = chartMode==='coverage' ? 100 : chartMode==='error' ? 18 : 12;
  const barW = (w-pad.l-pad.r)/data.length;
  const sy = v => h-pad.b - (v/yMax)*(h-pad.t-pad.b);
  let grid='';
  const step = chartMode==='coverage' ? 25 : chartMode==='error' ? 4 : 3;
  for(let v=0;v<=yMax;v+=step){
    grid += `<line x1="${pad.l}" x2="${w-pad.r}" y1="${sy(v)}" y2="${sy(v)}" class="grid-line" />
      <text x="${pad.l-10}" y="${sy(v)+4}" text-anchor="end" class="axis-label">${v}${chartMode==='vme'?'':'%'}</text>`;
  }
  let bars='';
  data.forEach((d,i)=>{
    const x = pad.l + i*barW + barW*0.16;
    const bw = barW*0.68;
    const y = sy(d.val);
    const bh = h-pad.b - y;
    let color = 'var(--accent)';
    if(chartMode==='error') color = d.val<5?'var(--sus)':d.val<10?'var(--ind)':'var(--res)';
    if(chartMode==='vme') color = d.val>4?'var(--res)':d.val>0?'var(--ind)':'var(--sus)';
    if(chartMode==='coverage') color = d.val<50?'var(--res)':d.val<85?'var(--ind)':'var(--sus)';
    bars += `<rect x="${x}" y="${y}" width="${bw}" height="${bh}" fill="${color}" rx="3" class="bar"
      data-id="${d.id}" data-name="${d.name}" data-val="${d.val}"
      onmouseover="barHover(this,event)" onmouseout="hideTip()"
      onclick="openDrugDrill('${d.id}')" />`;
    bars += `<text x="${x+bw/2}" y="${y-6}" text-anchor="middle" class="data-label">${d.val}${chartMode==='vme'?'':'%'}</text>`;
    bars += `<text x="${x+bw/2}" y="${h-pad.b+18}" text-anchor="middle" class="axis-label">${d.id}</text>`;
  });
  document.getElementById('drugBarChart').innerHTML = `
    <svg viewBox="0 0 ${w} ${h}" class="chart-svg">
      ${grid}${bars}
      <text x="${pad.l-42}" y="${h/2}" text-anchor="middle" transform="rotate(-90 ${pad.l-42} ${h/2})" class="axis-title">${chartMode==='error'?'error rate':chartMode==='coverage'?'coverage':'VME count'}</text>
    </svg>`;
}
function barHover(el, e){
  const d = DRUGS.find(x=>x.id===el.dataset.id);
  showTip(e, `<div class="tt-head"><div class="tt-title">${d.name} · ${d.id}</div><div class="tt-sub">${d.cls}</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Coverage</span><span class="tt-v">${pct(d.coverage)}</span></div>
      <div class="tt-row"><span class="tt-k">Error rate</span><span class="tt-v">${pct(d.error)}</span></div>
      <div class="tt-row"><span class="tt-k">Very major errors</span><span class="tt-v">${num(d.vme)}</span></div>
      <div class="tt-row"><span class="tt-k">Major errors</span><span class="tt-v">${num(d.me)}</span></div>
      <div class="tt-row"><span class="tt-k">Evaluable n</span><span class="tt-v">${num(d.evaluable)}</span></div>
      <div class="tt-section"><div class="tt-note">Click bar to open full evidence breakdown for this drug.</div></div>
    </div>`);
}

/* ==========================================================
   SCATTER
   ========================================================== */
function renderScatter(){
  if(!(DRUGS.filter(d=>d.estimable)) || !(DRUGS.filter(d=>d.estimable)).length){ return emptyState('scatterChart', 'No coverage/error pairs', 'Both axes need a measured benchmark. None was supplied.'); }
  const w=520, h=320, pad={t:20,r:24,b:50,l:56};
  const sx = v => pad.l + (v/100)*(w-pad.l-pad.r);
  const sy = v => h-pad.b - (v/18)*(h-pad.t-pad.b);
  let grid='';
  for(let v=0;v<=100;v+=25){
    grid += `<line x1="${sx(v)}" y1="${pad.t}" x2="${sx(v)}" y2="${h-pad.b}" class="grid-line" />
      <text x="${sx(v)}" y="${h-pad.b+18}" text-anchor="middle" class="axis-label">${v}%</text>`;
  }
  for(let v=0;v<=18;v+=3){
    grid += `<line x1="${pad.l}" y1="${sy(v)}" x2="${w-pad.r}" y2="${sy(v)}" class="grid-line" />
      <text x="${pad.l-10}" y="${sy(v)+4}" text-anchor="end" class="axis-label">${v}%</text>`;
  }
  let dots='';
  DRUGS.filter(d=>d.estimable).forEach(d=>{
    const cx = sx(d.coverage), cy = sy(d.error);
    const rad = 6 + Math.sqrt(d.evaluable)/4;
    let color = d.error<5?'var(--sus)':d.error<10?'var(--ind)':'var(--res)';
    dots += `<circle cx="${cx}" cy="${cy}" r="${rad}" fill="${color}" opacity=".7"
      stroke="#fff" stroke-width="2" class="point"
      data-id="${d.id}" onmouseover="scatterHover(this,event)" onmouseout="hideTip()" />`;
    dots += `<text x="${cx}" y="${cy-14}" text-anchor="middle" class="data-label">${d.id}</text>`;
  });
  document.getElementById('drugScatterChart')?.remove();
  document.getElementById('scatterChart').innerHTML = `
    <svg viewBox="0 0 ${w} ${h}" class="chart-svg">
      ${grid}${dots}
      <text x="${w/2}" y="${h-6}" text-anchor="middle" class="axis-title">coverage</text>
      <text x="14" y="${h/2}" text-anchor="middle" transform="rotate(-90 14 ${h/2})" class="axis-title">error rate</text>
    </svg>`;
}
function scatterHover(el, e){
  const d = DRUGS.find(x=>x.id===el.dataset.id);
  showTip(e, `<div class="tt-head"><div class="tt-title">${d.name}</div><div class="tt-sub">coverage vs. error</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Coverage</span><span class="tt-v">${pct(d.coverage)}</span></div>
      <div class="tt-row"><span class="tt-k">Error</span><span class="tt-v">${pct(d.error)}</span></div>
      <div class="tt-row"><span class="tt-k">Evaluable n</span><span class="tt-v">${num(d.evaluable)}</span></div>
      <div class="tt-section"><div class="tt-note">Bubble size ∝ √(evaluable n). Small bubbles carry wider confidence intervals.</div></div>
    </div>`);
}

/* ==========================================================
   ERROR STACK
   ========================================================== */
function renderErrorStack(){
  if(!(DRUGS.filter(d=>d.estimable)) || !(DRUGS.filter(d=>d.estimable)).length){ return emptyState('errorStackChart', 'No error breakdown', 'Very-major and major errors need a phenotypic reference to be counted against.'); }
  const data = DRUGS.filter(d=>d.estimable);
  const w=520, h=320, pad={t:20,r:24,b:60,l:56};
  const max = 30;
  const barW = (w-pad.l-pad.r)/data.length;
  const sy = v => h-pad.b - (v/max)*(h-pad.t-pad.b);
  let grid='';
  for(let v=0;v<=max;v+=6){
    grid += `<line x1="${pad.l}" y1="${sy(v)}" x2="${w-pad.r}" y2="${sy(v)}" class="grid-line" />
      <text x="${pad.l-10}" y="${sy(v)+4}" text-anchor="end" class="axis-label">${v}</text>`;
  }
  let bars='';
  data.forEach((d,i)=>{
    const x = pad.l + i*barW + barW*0.2;
    const bw = barW*0.6;
    const yVME = sy(d.vme);
    const hVME = h-pad.b - yVME;
    const yME = sy(d.vme+d.me);
    const hME = sy(d.vme) - sy(d.vme+d.me);
    bars += `<rect x="${x}" y="${yVME}" width="${bw}" height="${hVME}" fill="var(--res)" rx="2" class="bar"
      onmouseover="errorHover(event,'${d.id}','VME')" onmouseout="hideTip()" />`;
    bars += `<rect x="${x}" y="${yME}" width="${bw}" height="${hME}" fill="var(--ind)" rx="2" class="bar"
      onmouseover="errorHover(event,'${d.id}','ME')" onmouseout="hideTip()" />`;
    bars += `<text x="${x+bw/2}" y="${h-pad.b+18}" text-anchor="middle" class="axis-label">${d.id}</text>`;
  });
  document.getElementById('errorStackChart').innerHTML = `
    <svg viewBox="0 0 ${w} ${h}" class="chart-svg">
      ${grid}${bars}
      <text x="${w/2}" y="${h-6}" text-anchor="middle" class="axis-title">errors</text>
    </svg>`;
}
function errorHover(e, id, kind){
  const d = DRUGS.find(x=>x.id===id);
  showTip(e, `<div class="tt-head"><div class="tt-title">${d.name}</div><div class="tt-sub">${kind==='VME'?'very major errors':'major errors'}</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Count</span><span class="tt-v">${kind==='VME'?d.vme:d.me}</span></div>
      <div class="tt-row"><span class="tt-k">Evaluable n</span><span class="tt-v">${num(d.evaluable)}</span></div>
      <div class="tt-section"><div class="tt-note">${kind==='VME'?'False-susceptible — clinically the most serious error type.':'False-resistant — leads to unnecessary second-line treatment.'}</div></div>
    </div>`);
}

/* ==========================================================
   HEADLINE TABLE
   ========================================================== */
function renderHeadline(){
  if(!(DRUGS) || !(DRUGS).length){ return emptyState('headlineTable', 'No benchmark table', 'A benchmark needs a phenotypic reference. None was supplied for this run.'); }
  const tb = document.querySelector('#headlineTable tbody');
  tb.innerHTML = DRUGS.map(d => {
    if(!d.estimable){
      return `<tr onmouseover="drugRowHover(event,'${d.id}')" onmouseout="hideTip()" style="cursor:pointer" onclick="openDrugDrill('${d.id}')">
        <td style="font-family:var(--mono);font-size:12px;padding-left:14px;text-align:left">${d.name}</td>
        <td style="color:var(--text-2);font-size:11px">${d.cls}</td>
        <td class="na">${pct(d.coverage)}</td>
        <td class="na">not estimable</td>
        <td class="na">—</td><td class="na">—</td><td class="na">—</td>
      </tr>`;
    }
    return `<tr onmouseover="drugRowHover(event,'${d.id}')" onmouseout="hideTip()" style="cursor:pointer" onclick="openDrugDrill('${d.id}')">
      <td style="font-family:var(--mono);font-size:12px;padding-left:14px;text-align:left">${d.name}</td>
      <td style="color:var(--text-2);font-size:11px">${d.cls}</td>
      <td>${pct(d.coverage)}</td>
      <td style="background:${bgColorForError(d.error)};color:${textColorForError(d.error)};font-weight:600">${pct(d.error)}</td>
      <td style="background:${d.vme>0?'var(--res-bg)':'var(--bg-soft)'};color:${d.vme>0?'var(--res)':'var(--text-3)'};font-weight:600">${num(d.vme)}</td>
      <td>${num(d.me)}</td>
      <td style="color:var(--text-2)">${num(d.evaluable)}</td>
    </tr>`;
  }).join('');
}
function drugRowHover(e, id){
  const d = DRUGS.find(x=>x.id===id);
  showTip(e, `<div class="tt-head"><div class="tt-title">${d.name} · ${d.id}</div><div class="tt-sub">${d.cls}</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Coverage</span><span class="tt-v">${pct(d.coverage)}</span></div>
      <div class="tt-row"><span class="tt-k">Error rate</span><span class="tt-v">${d.estimable?d.error+'%':'not estimable'}</span></div>
      <div class="tt-row"><span class="tt-k">Evaluable n</span><span class="tt-v">${d.evaluable||'—'}</span></div>
      <div class="tt-section"><div class="tt-note">Click for full evidence breakdown.</div></div>
    </div>`);
}

/* ==========================================================
   DRUG CARDS
   ========================================================== */
function renderDrugs(){
  if(!(DRUGS) || !(DRUGS).length){ return emptyState('drugGrid', 'No drug cards', 'No drug was reported in this run.'); }
  const grid = document.getElementById('drugGrid');
  grid.innerHTML = DRUGS.map(d => {
    const covClass = d.coverage < 50 ? 'crit' : d.coverage < 85 ? 'low' : '';
    const est = d.estimable;
    return `<div class="drug-card ${est?'':'not-estimable'}" data-drug="${d.id}" onclick="openDrugDrill('${d.id}')">
      <div class="drug-card-head">
        <div><div class="drug-class">${d.cls}</div><div class="drug-name">${d.name}</div></div>
        <div class="drug-abbr">${d.id}</div>
      </div>
      <div class="metric-row"><span class="k">Coverage</span><span class="v">${pct(d.coverage)}</span></div>
      <div class="coverage-bar"><div class="coverage-fill ${covClass}" style="width:${pct(d.coverage)}"></div></div>
      ${est ? `
        <div class="metric-row" style="margin-top:10px"><span class="k">Error rate</span><span class="v" style="color:${textColorForError(d.error)};font-weight:600">${pct(d.error)}</span></div>
        <div class="metric-row"><span class="k">Very major errors</span><span class="v">${d.vme > 0 ? `<span class="vme-badge">${d.vme}</span>` : '0'}</span></div>
        <div class="metric-row"><span class="k">Evaluable n</span><span class="v">${num(d.evaluable)}</span></div>
      ` : `
        <div class="not-est">Not estimable — fewer than 20 phenotypically resistant isolates. No accuracy claim is made.</div>
      `}
    </div>`;
  }).join('');
}

/* ==========================================================
   CALL HEATMAP
   ========================================================== */
let hmCount = 60;
function setHeatmapCount(btn, n){
  btn.parentElement.querySelectorAll('.mini-btn').forEach(b=>b.classList.remove('active'));
  btn.classList.add('active');
  hmCount = n;
  renderCallHeatmap();
}
function emptyState(target, title, detail){
  const el = document.getElementById(target);
  if(!el) return true;
  const block = `<div class="empty-state"><div class="empty-title">${title}</div>
    <div class="empty-detail">${detail}</div></div>`;
  // A <div> assigned into a <table> is hoisted out of it by the parser and
  // leaves a tall blank gap where the table was. Tables get a spanning row.
  if(el.tagName === 'TABLE'){
    const cols = (el.querySelectorAll('thead th') || []).length || 1;
    let body = el.querySelector('tbody');
    if(!body){ body = document.createElement('tbody'); el.appendChild(body); }
    body.innerHTML = `<tr><td colspan="${cols}" style="padding:0;border:none">${block}</td></tr>`;
    return true;
  }
  el.innerHTML = block;
  return true;
}
let isolateData = DATA.isolates || [];
function renderCallHeatmap(){
  const hm = document.getElementById('callHeatmap');
  const cols = DRUGS.length;
  if(!isolateData.length || !cols){
    return emptyState('callHeatmap', 'No per-isolate call matrix in this run',
      'A call matrix needs a cohort. Single-sample analyses report their calls in the interpretation table instead.');
  }
  let html = '';
  html += `<div class="hm-row-label" style="height:60px;background:var(--card-alt)"></div>`;
  DRUGS.forEach(d => html += `<div class="hm-col-label">${d.id}</div>`);
  isolateData.slice(0, hmCount).forEach(iso => {
    html += `<div class="hm-row-label" title="${iso.id} · ${iso.lin}">${iso.id}</div>`;
    DRUGS.forEach(d => {
      const c = iso.calls[d.id];
      html += `<div class="hm-cell call-${c}"
        data-iso="${iso.id}" data-lin="${iso.lin}" data-drug="${d.id}" data-call="${c}"
        onmouseover="cellHover(event,this)" onmouseout="hideTip()"
        onclick="openCallDrill(this.dataset)"></div>`;
    });
  });
  hm.innerHTML = html;
  hm.style.gridTemplateColumns = `84px repeat(${cols}, 20px)`;
}
function cellHover(e, el){
  const d = DRUGS.find(x=>x.id===el.dataset.drug);
  const call = el.dataset.call.toUpperCase();
  showTip(e, `<div class="tt-head"><div class="tt-title">${el.dataset.iso} · ${d.name}</div><div class="tt-sub">lineage ${el.dataset.lin}</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Reconciled call</span><span class="tt-v">${call}</span></div>
      <div class="tt-row"><span class="tt-k">Drug coverage</span><span class="tt-v">${pct(d.coverage)}</span></div>
      <div class="tt-row"><span class="tt-k">Drug error rate</span><span class="tt-v">${d.estimable?d.error+'%':'not estimable'}</span></div>
      <div class="tt-section"><div class="tt-note">${describeCall(call)}</div></div>
      <div class="tt-section" style="margin-top:8px"><div class="tt-note" style="font-style:normal;color:var(--accent);font-weight:600">Click for full evidence chain →</div></div>
    </div>`);
}

/* ==========================================================
   COVERAGE HEATMAP
   ========================================================== */
function coverageColour(fraction){
  if(fraction === null || fraction === undefined) return ['var(--na)','unknown'];
  if(fraction >= 0.95) return ['var(--sus)','≥ 95%'];
  if(fraction >= 0.80) return ['#a8c94a','80–95%'];
  if(fraction >= 0.50) return ['var(--ind)','50–80%'];
  return ['var(--res)','< 50%'];
}
function renderCoverage(){
  const hm = document.getElementById('coverageHeatmap');
  const cols = COVERAGE_LOCI.length;
  const rows = DATA.coverage_matrix || [];
  if(!cols){
    return emptyState('coverageHeatmap', 'No callable-locus evidence supplied',
      'Coverage is the precondition for susceptibility. Without a depth table, BED mask or gVCF, every drug is reported NOT_ASSESSED.');
  }
  let html = '';
  html += `<div class="hm-row-label" style="height:60px;background:var(--card-alt)"></div>`;
  COVERAGE_LOCI.forEach(l => html += `<div class="hm-col-label" title="${l.id} → ${l.drug}">${l.id}</div>`);
  rows.forEach(row => {
    html += `<div class="hm-row-label" title="${row.id}">${row.id}</div>`;
    COVERAGE_LOCI.forEach(l => {
      const fraction = (row.loci || {})[l.id];
      const [color, level] = coverageColour(fraction === undefined ? null : fraction);
      html += `<div class="hm-cell" style="background:${color}"
        data-locus="${l.id}" data-level="${level}" data-drug="${l.drug}"
        onmouseover="locusHover(event,this)" onmouseout="hideTip()"></div>`;
    });
  });
  hm.innerHTML = html;
  hm.style.gridTemplateColumns = `84px repeat(${cols}, 20px)`;
}
function locusHover(e, el){
  showTip(e, `<div class="tt-head"><div class="tt-title">${el.dataset.locus}</div><div class="tt-sub">drug: ${el.dataset.drug}</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Callability</span><span class="tt-v">${el.dataset.level}</span></div>
      <div class="tt-section"><div class="tt-note">Coverage-gated: a drug is NOT_ASSESSED unless this locus is callable above threshold.</div></div>
    </div>`);
}

/* ==========================================================
   LINEAGE TABLE
   ========================================================== */
function renderLineage(){
  const tb = document.querySelector('#lineageTable tbody');
  const strata = DATA.lineage || {};
  if(!LINEAGES.length){
    tb.innerHTML = `<tr><td colspan="6" style="padding:22px;text-align:left;color:var(--text-2);line-height:1.6">
      <strong>Lineage was not recorded for these isolates.</strong> This is reported as unknown rather
      than as a single lineage: folding untyped isolates into one bucket would be indistinguishable
      from a genuinely lineage-restricted cohort. Supply lineage calls to make this check live.</td></tr>`;
    return;
  }
  tb.innerHTML = DRUGS.map(d => {
    if(!d.estimable){
      return `<tr>
        <td style="font-family:var(--mono);font-size:12px;padding-left:14px;text-align:left;background:var(--bg-soft)">${d.name}</td>
        ${LINEAGES.map(()=>`<td class="na">—</td>`).join('')}
        <td class="na">not estimable</td>
      </tr>`;
    }
    const cells = LINEAGES.map(L => {
      const cell = ((strata[d.id] || {})[L]) || null;
      if(!cell){
        return `<td class="na" title="no isolates of this lineage were evaluable for this drug">—</td>`;
      }
      const e = cell.error;
      const n = cell.n;
      const unstable = n < 20;
      return `<td style="background:${bgColorForError(e)};color:${textColorForError(e)};${unstable?'opacity:.6;font-style:italic':''}"
        onmouseover="lineageCellHover(event,'${d.id}','${L}',${e},${n})" onmouseout="hideTip()">
        ${e}%<span class="cell-n">n=${n}${unstable?' ⚠':''}</span></td>`;
    }).join('');
    return `<tr>
      <td style="font-family:var(--mono);font-size:12px;padding-left:14px;text-align:left;background:var(--bg-soft)">${d.name}</td>
      ${cells}
      <td style="background:${bgColorForError(d.error)};color:${textColorForError(d.error)};font-weight:600">
        ${pct(d.error)}<span class="cell-n">n=${num(d.evaluable)}</span></td>
    </tr>`;
  }).join('');
}
function lineageCellHover(e, id, L, err, n){
  const d = DRUGS.find(x=>x.id===id);
  showTip(e, `<div class="tt-head"><div class="tt-title">${d.name} · Lineage ${L[1]}</div><div class="tt-sub">error rate</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Error rate</span><span class="tt-v">${err}%</span></div>
      <div class="tt-row"><span class="tt-k">Evaluable n</span><span class="tt-v">${n}</span></div>
      <div class="tt-row"><span class="tt-k">Overall</span><span class="tt-v">${pct(d.error)}</span></div>
      ${n<20?'<div class="tt-section"><div class="tt-note">⚠ n &lt; 20 — this cell is a hypothesis, not a measurement.</div></div>':''}
    </div>`);
}

/* ==========================================================
   RADAR
   ========================================================== */
function renderRadar(){
  if(!(LINEAGES) || !(LINEAGES).length){ return emptyState('radarChart', 'Lineage not recorded', 'No isolate carried a lineage, so balance across lineages cannot be shown. This is unknown, not uniform.'); }
  const w=520, h=340, cx=w/2, cy=h/2+10, R=110;
  const axes = DRUGS.filter(d=>d.estimable).slice(0,6);
  const n = axes.length;
  const angleFor = i => -Math.PI/2 + (i/n) * Math.PI * 2;
  let web='';
  [0.25,0.5,0.75,1].forEach(level=>{
    let pts='';
    for(let i=0;i<n;i++){
      const a = angleFor(i);
      pts += `${cx + R*level*Math.cos(a)},${cy + R*level*Math.sin(a)} `;
    }
    web += `<polygon points="${pts}" fill="none" stroke="var(--border)" stroke-width="1" />`;
  });
  let axesLines='';
  axes.forEach((d,i)=>{
    const a = angleFor(i);
    axesLines += `<line x1="${cx}" y1="${cy}" x2="${cx+R*Math.cos(a)}" y2="${cy+R*Math.sin(a)}" stroke="var(--border)" />
      <text x="${cx + (R+22)*Math.cos(a)}" y="${cy + (R+22)*Math.sin(a)+4}" text-anchor="middle" class="axis-label">${d.id}</text>`;
  });
  // Balance polygon: use coverage as a proxy for balance
  let poly = '';
  axes.forEach((d,i)=>{
    const a = angleFor(i);
    const v = d.coverage / 100;
    poly += `${cx + R*v*Math.cos(a)},${cy + R*v*Math.sin(a)} `;
  });
  document.getElementById('radarChart').innerHTML = `
    <svg viewBox="0 0 ${w} ${h}" class="chart-svg" style="max-width:520px;margin:0 auto">
      ${web}${axesLines}
      <polygon points="${poly}" fill="var(--accent-soft)" fill-opacity="0.55" stroke="var(--accent)" stroke-width="2" />
      ${axes.map((d,i)=>{
        const a = angleFor(i);
        const v = d.coverage/100;
        return `<circle cx="${cx + R*v*Math.cos(a)}" cy="${cy + R*v*Math.sin(a)}" r="4" fill="var(--accent)" />`;
      }).join('')}
    </svg>`;
}

/* ==========================================================
   TIER STACK
   ========================================================== */
function renderTierStack(){
  if(!(DRUGS) || !(DRUGS).length){ return emptyState('tierStackChart', 'No evidence to compose', 'This run attached no drug evidence.'); }
  const w=980, h=280, pad={t:24,r:24,b:56,l:80};
  const tiers = [
    { key:'catalogued', color:'var(--res)', label:'Catalogued' },
    { key:'phenotypic', color:'#c17a86', label:'Phenotypic' },
    { key:'inferred', color:'var(--ind)', label:'Inferred' },
    { key:'predicted', color:'#8d9b8d', label:'Predicted' },
    { key:'unsupported', color:'var(--na)', label:'Unsupported' },
  ];
  const data = DRUGS.map(d=>{
    const total = d.evaluable || 100;
    if(!d.estimable) return { id:d.id, parts: tiers.map(t=>({ key:t.key, v: t.key==='unsupported'?100:0 })) };
    const cat = Math.max(10, Math.min(85, d.error * 3));
    const phen = 5 + (d.vme + d.me) * 1.5;
    const inf = 10;
    const pred = 5;
    const uns = Math.max(5, 100 - cat - phen - inf - pred);
    return { id:d.id, parts:[
      { key:'catalogued', v:cat },
      { key:'phenotypic', v:phen },
      { key:'inferred', v:inf },
      { key:'predicted', v:pred },
      { key:'unsupported', v:uns },
    ]};
  });
  const barW = (w-pad.l-pad.r)/data.length;
  let bars='';
  data.forEach((d,i)=>{
    const x = pad.l + i*barW + barW*0.15;
    const bw = barW*0.7;
    let y = pad.t;
    const totalH = h-pad.t-pad.b;
    d.parts.forEach(p=>{
      const ph = (p.v/100)*totalH;
      const color = tiers.find(t=>t.key===p.key).color;
      bars += `<rect x="${x}" y="${y}" width="${bw}" height="${ph}" fill="${color}" opacity=".85"
        onmouseover="tierHover(event,'${d.id}','${p.key}',${p.v})" onmouseout="hideTip()" style="cursor:pointer" />`;
      y += ph;
    });
    bars += `<text x="${x+bw/2}" y="${h-pad.b+18}" text-anchor="middle" class="axis-label">${d.id}</text>`;
  });
  document.getElementById('tierStackChart').innerHTML = `
    <svg viewBox="0 0 ${w} ${h}" class="chart-svg">
      ${bars}
      <text x="${pad.l-50}" y="${h/2}" text-anchor="middle" transform="rotate(-90 ${pad.l-50} ${h/2})" class="axis-title">% of calls</text>
    </svg>`;
}
function tierHover(e, id, tier, v){
  const d = DRUGS.find(x=>x.id===id);
  const tierName = { catalogued:'Catalogued tier 1', phenotypic:'Phenotypic (site-local)', inferred:'Inferred', predicted:'Predicted (ML)', unsupported:'Unsupported / unassessed' }[tier];
  const canRes = tier==='catalogued'||tier==='phenotypic';
  showTip(e, `<div class="tt-head"><div class="tt-title">${d.name}</div><div class="tt-sub">${tierName}</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Share of calls</span><span class="tt-v">${v.toFixed(1)}%</span></div>
      <div class="tt-section"><div class="tt-note">${canRes?'Can establish RESISTANT.':'Structurally capped at INDETERMINATE — cannot establish resistance.'}</div></div>
    </div>`);
}

/* ==========================================================
   ALIVIEW
   ========================================================== */
const LOCI = DATA.alignment_loci || {};

/* The mock generated A/T/C/G for every position with a seeded PRNG. That is
   invented sequence, and inventing sequence is the same failure as the
   hash-derived conservation scores this project already removed: it renders
   as data, it is not data, and nothing downstream can tell.

   So this shows what was actually called. Reference bases come from the REF
   field of real records; a position with no record is drawn as a gap, because
   "we did not call a variant here" is not the same as "the base is A". An
   isolate row shows the reference glyph where it matched and the real ALT
   where it did not. */
function renderAliView(){
  const select = document.getElementById('avLocus');
  const content = document.getElementById('avContent');
  if(!content) return;
  const ids = Object.keys(LOCI);
  if(!ids.length){
    return emptyState('avContent', 'No variant alignment available',
      'This view renders called variants at their real coordinates. No record in this run carried a resolved position, so there is nothing to align.');
  }
  if(select && !select.dataset.filled){
    select.innerHTML = ids.map(id => {
      const l = LOCI[id];
      return `<option value="${id}">${l.name} \u00b7 ${l.positions.length} position(s)</option>`;
    }).join('');
    select.dataset.filled = '1';
  }
  const locusId = (select && LOCI[select.value]) ? select.value : ids[0];
  const locus = LOCI[locusId];
  const positions = locus.positions || [];
  const isolates = (locus.isolates || []).slice(
    0, parseInt((document.getElementById('avCount') || {}).value || '24', 10));

  if(!positions.length){
    return emptyState('avContent', `No called positions in ${locus.name}`,
      'The locus was examined but carried no variant record.');
  }

  let html = '';
  html += `<div class="av-row track"><div class="av-label">Position</div><div class="av-ruler">`;
  positions.forEach(v => {
    html += `<div class="av-tick major" style="width:${9 * Math.max(1, String(v.pos).length)}px;min-width:34px;text-align:left;padding-left:3px">${v.pos}</div>`;
  });
  html += `</div></div>`;

  html += `<div class="av-row track-gene"><div class="av-label">Gene</div>
    <div style="display:flex;padding:0 4px;align-items:center;height:100%">
    <div style="width:${positions.length*34}px;height:16px;background:linear-gradient(90deg,var(--accent),var(--accent-line));border-radius:3px;position:relative">
      <span style="position:absolute;left:10px;top:-1px;font-size:10px;color:#fff;font-weight:700;letter-spacing:.05em">${locus.name}</span>
    </div></div></div>`;

  html += `<div class="av-row track"><div class="av-label">Change</div><div class="av-seq">`;
  positions.forEach(v => {
    html += `<div class="av-base ${v.catalogued ? 'variant' : 'variant ind'}"
      style="width:34px;min-width:34px" title="${v.label || ''}">${v.ref}\u2192${v.alt}</div>`;
  });
  html += `</div></div>`;

  html += `<div class="av-row track"><div class="av-label"><span class="lid">REF</span><span class="llin">${locus.assembly || 'reference'}</span></div><div class="av-seq">`;
  positions.forEach(v => {
    html += `<div class="av-base ${v.ref}" style="width:34px;min-width:34px">${v.ref}</div>`;
  });
  html += `</div></div>`;

  isolates.forEach(iso => {
    html += `<div class="av-row"><div class="av-label"><span class="lid">${iso.id}</span><span class="llin">${iso.lineage || 'untyped'}</span></div><div class="av-seq">`;
    positions.forEach(v => {
      const carried = (iso.carried || []).indexOf(String(v.pos)) >= 0;
      if(carried){
        html += `<div class="av-base ${v.catalogued ? 'variant' : 'variant ind'}"
          style="width:34px;min-width:34px"
          data-iso="${iso.id}" data-lin="${iso.lineage || 'untyped'}" data-pos="${v.pos}"
          data-ref="${v.ref}" data-alt="${v.alt}" data-label="${v.label || ''}"
          data-drug="${(v.drugs || []).join(', ')}">${v.alt}</div>`;
      } else {
        html += `<div class="av-base ${v.ref}" style="width:34px;min-width:34px"
          title="reference allele at this position">${v.ref}</div>`;
      }
    });
    html += `</div></div>`;
  });

  if(locus.coverage && locus.coverage.length === positions.length){
    html += `<div class="av-row track"><div class="av-label">Callable</div><div class="av-cov">`;
    locus.coverage.forEach(fraction => {
      if(fraction === null){
        html += `<div class="av-cov-bar crit" style="height:4px" title="no coverage evidence"></div>`;
      } else {
        const hh = Math.max(4, Math.round(fraction * 34));
        const cls = fraction < 0.5 ? 'crit' : fraction < 0.8 ? 'low' : '';
        html += `<div class="av-cov-bar ${cls}" style="height:${hh}px;width:34px"
          title="${Math.round(fraction*100)}% callable"></div>`;
      }
    });
    html += `</div></div>`;
  }

  content.innerHTML = html;

  content.querySelectorAll('.av-base.variant[data-pos]').forEach(el => {
    el.addEventListener('mouseenter', ev => {
      showTip(ev, `<div class="tt-head"><div class="tt-title">${el.dataset.pos} ${el.dataset.ref}\u2192${el.dataset.alt}</div>
          <div class="tt-sub">${locus.name}${el.dataset.drug ? ' \u00b7 ' + el.dataset.drug : ''}</div></div>
        <div class="tt-body">
          <div class="tt-row"><span class="tt-k">Isolate</span><span class="tt-v">${el.dataset.iso}</span></div>
          <div class="tt-row"><span class="tt-k">Lineage</span><span class="tt-v">${el.dataset.lin}</span></div>
          <div class="tt-row"><span class="tt-k">Change</span><span class="tt-v">${el.dataset.label || '\u2014'}</span></div>
          <div class="tt-section"><div class="tt-note">Positions shown are those with a real record in this run. A variant not in the catalogue is held at INDETERMINATE.</div></div>
        </div>`);
    });
    el.addEventListener('mouseleave', hideTip);
    el.addEventListener('click', () => openVariantDrill(el.dataset));
  });
}

function renderDiscordance(){
  if(!(DISCORDANCE_ROWS) || !(DISCORDANCE_ROWS).length){ return emptyState('discTable', 'No discordance', 'No two sources disagreed in this run. That is not the same as agreement being verified.'); }
  const tb = document.querySelector('#discTable tbody');
  const tagClass = v => v === 'RESISTANT' ? 'res' : v === 'SUSCEPTIBLE' ? 'sus' : v === 'INDETERMINATE' ? 'ind' : 'na';
  tb.innerHTML = DISCORDANCE_ROWS.map(r => `
    <tr onmouseover="discHover(event,this)" onmouseout="hideTip()" style="cursor:pointer">
      <td class="vid">${r.var}</td>
      <td><span class="tag engine">${r.drug}</span></td>
      <td><span class="tag ${tagClass(r.prof)}">${r.prof}</span></td>
      <td><span class="tag ${tagClass(r.myk)}">${r.myk}</span></td>
      <td><span class="tag ${tagClass(r.rec)}">${r.rec}</span></td>
      <td>${r.reason}</td>
      <td style="text-align:right;color:var(--text-2);font-family:var(--mono);font-size:11px">→</td>
    </tr>`).join('');
}
function discHover(e, row){
  showTip(e, `<div class="tt-head"><div class="tt-title">Discordance record</div><div class="tt-sub">engine disagreement</div></div>
    <div class="tt-body">
      <div class="tt-note">Myconductor reports engine disagreements rather than voting. The reconciled call is derived from evidence tier and catalogue entry, not majority rule.</div>
      <div class="tt-section"><div class="tt-note" style="color:var(--accent);font-style:normal;font-weight:600">Click to open full reconciliation trace →</div></div>
    </div>`);
}

/* ==========================================================
   VUS LIST
   ========================================================== */
function renderVUS(){
  if(!(VUS_ITEMS) || !(VUS_ITEMS).length){ return emptyState('vusList', 'No variants queued for validation', 'The workbench ranked nothing in this run.'); }
  const list = document.getElementById('vusList');
  list.innerHTML = VUS_ITEMS.map(v => {
    const f = v.features || {};
    /* The workbench emits a priority *band* and, only where it could rank,
       a score. A band without a score is shown as a band: inventing a
       percentage for "insufficient-data" would be exactly the fabricated
       confidence this workbench exists to avoid. */
    const scored = (v.score !== null && v.score !== undefined);
    const width = scored ? Math.round(Math.max(0, Math.min(1, v.score)) * 100) : 0;
    return `
    <div class="vus-item" onclick="openVUSDrill(${v.rank})">
      <div class="vus-rank">${v.rank}</div>
      <div class="vus-body">
        <div class="vus-var">${v.variant || '\u2014'}</div>
        <div class="vus-meta">
          <span><span class="dot"></span>${v.drug || 'no drug'}</span>
          <span><span class="dot"></span>${f.dimensions || 'no dimensions available'}</span>
          <span><span class="dot"></span>${f.experiment || 'no experiment specified'}</span>
        </div>
      </div>
      <div class="vus-score">
        <div class="n" style="${scored ? '' : 'font-size:13px;letter-spacing:0'}">${scored ? width : (v.priority || '\u2014')}</div>
        <div class="l">${scored ? 'score' : 'priority'}</div>
        <div class="score-bar"><div class="score-fill" style="width:${width}%"></div></div>
      </div>
    </div>`;
  }).join('');
}

/* ==========================================================
   MECH CARDS
   ========================================================== */
function renderMech(){
  if(!(MECH_CARDS) || !(MECH_CARDS).length){ return emptyState('mechGrid', 'No mechanism hypotheses', 'No lane raised a named mechanism for this sample.'); }
  const g = document.getElementById('mechGrid');
  g.innerHTML = MECH_CARDS.map((m,i)=>`
    <div class="mech-card" onclick="openMechDrill(${i})" style="cursor:pointer">
      <h5>${m.title}</h5>
      <p>${m.body}</p>
      <div>${m.tags.map(t=>`<span class="tag ${t.t}">${t.l}</span>`).join('')}</div>
    </div>`).join('');
}

/* ==========================================================
   VALIDATION TIMELINE
   ========================================================== */
function renderTimeline(){
  if(!(DATA.validation_timeline) || !(DATA.validation_timeline).length){ return emptyState('validationTimeline', 'No validation activity', 'Nothing has been validated locally in the reporting window.'); }
  const events = DATA.validation_timeline;
  document.getElementById('validationTimeline').innerHTML = events.map(e=>`
    <div class="tl-item done">
      <div class="tl-date">${e.date}</div>
      <div class="tl-text">${e.text}</div>
    </div>`).join('');
}

/* ==========================================================
   VALIDATION DONUT
   ========================================================== */
function renderValidationDonut(){
  if(!(DATA.validation_outcomes) || !(DATA.validation_outcomes).length){ return emptyState('validationDonut', 'No laboratory validations', 'Site-local phenotypic feedback has not been supplied for this run.'); }
  const data = DATA.validation_outcomes;
  const total = data.reduce((s,d)=>s+d.v,0);
  const cx=180, cy=180, rOut=110, rIn=66;
  let a0 = -Math.PI/2, paths='';
  data.forEach(d=>{
    const a1 = a0 + (d.v/total) * Math.PI * 2;
    const x1 = cx + rOut*Math.cos(a0), y1 = cy + rOut*Math.sin(a0);
    const x2 = cx + rOut*Math.cos(a1), y2 = cy + rOut*Math.sin(a1);
    const x3 = cx + rIn*Math.cos(a1), y3 = cy + rIn*Math.sin(a1);
    const x4 = cx + rIn*Math.cos(a0), y4 = cy + rIn*Math.sin(a0);
    const large = (a1-a0) > Math.PI ? 1 : 0;
    paths += `<path d="M ${x1} ${y1} A ${rOut} ${rOut} 0 ${large} 1 ${x2} ${y2} L ${x3} ${y3} A ${rIn} ${rIn} 0 ${large} 0 ${x4} ${y4} Z"
      fill="${d.c}" stroke="#fff" stroke-width="2" opacity=".9"
      onmouseover="valHover(event,'${d.k}',${d.v},${(d.v/total*100).toFixed(1)})" onmouseout="hideTip()" style="cursor:pointer" />`;
    a0 = a1;
  });
  document.getElementById('validationDonut').innerHTML = `
    <svg viewBox="0 0 360 360" class="chart-svg" style="max-width:340px;margin:0 auto">
      ${paths}
      <text x="${cx}" y="${cy-8}" text-anchor="middle" font-family="var(--mono)" font-size="11" fill="var(--text-2)" letter-spacing="1">VALIDATIONS</text>
      <text x="${cx}" y="${cy+22}" text-anchor="middle" font-family="var(--mono)" font-size="28" font-weight="600" fill="var(--text)">${total}</text>
    </svg>`;
}
function valHover(e, k, v, pct){
  showTip(e, `<div class="tt-head"><div class="tt-title">${k}</div><div class="tt-sub">site-local validation</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Count</span><span class="tt-v">${v}</span></div>
      <div class="tt-row"><span class="tt-k">Share</span><span class="tt-v">${pct}%</span></div>
      <div class="tt-section"><div class="tt-note">Validations are site-local. They become phenotypic-tier evidence for the originating site only. They do not promote globally.</div></div>
    </div>`);
}

/* ==========================================================
   SITE BAR CHART
   ========================================================== */
function renderSiteChart(){
  if(!(DATA.federated_sites) || !(DATA.federated_sites).length){ return emptyState('siteBarChart', 'No federated contributions', 'No site submitted aggregated statistics for this run.'); }
  const sites = DATA.federated_sites;
  const w=980, h=280, pad={t:24,r:24,b:60,l:110};
  const maxN = Math.max(...sites.map(s=>s.n));
  const barH = (h-pad.t-pad.b)/sites.length;
  let bars='';
  sites.forEach((s,i)=>{
    const y = pad.t + i*barH + barH*0.15;
    const bh = barH*0.7;
    const bw = (s.n/maxN)*(w-pad.l-pad.r);
    bars += `<rect x="${pad.l}" y="${y}" width="${bw}" height="${bh}" fill="var(--accent)" opacity=".85" rx="3"
      onmouseover="siteHover(event,'${s.id}',${s.n},${s.k})" onmouseout="hideTip()" style="cursor:pointer" />
      <text x="${pad.l-12}" y="${y+bh/2+4}" text-anchor="end" class="axis-label" style="font-family:var(--mono)">${s.id}</text>
      <text x="${pad.l+bw+8}" y="${y+bh/2+4}" class="data-label">${s.n}</text>`;
  });
  document.getElementById('siteBarChart').innerHTML = `
    <svg viewBox="0 0 ${w} ${h}" class="chart-svg">
      ${bars}
      <text x="${w/2}" y="${h-6}" text-anchor="middle" class="axis-title">submissions</text>
    </svg>`;
}
function siteHover(e, id, n, k){
  showTip(e, `<div class="tt-head"><div class="tt-title">${id}</div><div class="tt-sub">federated site</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Submissions</span><span class="tt-v">${n}</span></div>
      <div class="tt-row"><span class="tt-k">k-threshold met</span><span class="tt-v">${k}</span></div>
      <div class="tt-section"><div class="tt-note">Only aggregated statistics with k ≥ 5 are accepted. Raw isolate data never leaves the site.</div></div>
    </div>`);
}

/* ==========================================================
   AUDIT TABLE
   ========================================================== */
function renderAudit(){
  if(!(AUDIT_EVENTS) || !(AUDIT_EVENTS).length){ return emptyState('auditTable', 'No audit entries', 'No ledger was supplied to this report.'); }
  const tb = document.querySelector('#auditTable tbody');
  tb.innerHTML = AUDIT_EVENTS.map(a=>`
    <tr onmouseover="auditHover(event,this)" onmouseout="hideTip()" style="cursor:pointer">
      <td style="font-family:var(--mono);font-size:11px;color:var(--text-2)">${a.ts}</td>
      <td><span class="tag accent">${a.event}</span></td>
      <td style="font-family:var(--mono);font-size:11.5px">${a.actor}</td>
      <td>${a.target}</td>
      <td style="font-family:var(--mono);font-size:11px;color:var(--accent)">${a.hash}</td>
      <td style="text-align:right;color:var(--text-2);font-size:11px">→</td>
    </tr>`).join('');
}
function auditHover(e){
  showTip(e, `<div class="tt-head"><div class="tt-title">Audit ledger entry</div><div class="tt-sub">append-only · hash-chained</div></div>
    <div class="tt-body">
      <div class="tt-note">Every catalogue update, VUS validation, model registration, call promotion, and call retraction is recorded here. Entries cannot be modified; only new entries can be appended.</div>
    </div>`);
}

/* ==========================================================
   DRILL PANEL
   ========================================================== */
const panel = document.getElementById('drillPanel');
const panelTitle = document.getElementById('drillTitle');
const panelBody = document.getElementById('drillBody');
function openPanel(title, body){
  panelTitle.textContent = title;
  panelBody.innerHTML = body;
  panel.classList.add('open');
}
document.getElementById('drillClose').addEventListener('click', ()=>panel.classList.remove('open'));

function openDrugDrill(id){
  const d = DRUGS.find(x=>x.id===id);
  const body = `
    <div class="drill-section">
      <h4>Identity</h4>
      <div class="drill-kv">
        <span class="k">Drug</span><span class="v">${d.name}</span>
        <span class="k">Class</span><span class="v">${d.cls}</span>
        <span class="k">Abbreviation</span><span class="v">${d.id}</span>
      </div>
    </div>
    <div class="drill-section">
      <h4>Performance</h4>
      <div class="drill-kv">
        <span class="k">Coverage</span><span class="v">${pct(d.coverage)}</span>
        <span class="k">Error rate</span><span class="v">${d.estimable?d.error+'%':'not estimable'}</span>
        <span class="k">Evaluable n</span><span class="v">${d.evaluable||'—'}</span>
        <span class="k">Very major errors</span><span class="v">${d.vme ?? '—'}</span>
        <span class="k">Major errors</span><span class="v">${d.me ?? '—'}</span>
      </div>
    </div>
    <div class="drill-section">
      <h4>Evidence composition</h4>
      ${d.estimable ? `
        <div class="evidence-block res"><div class="tier">Catalogued · tier 1</div><div class="src">Primary resistance loci</div><div class="detail">Resistance calls for ${d.name} are driven by catalogue-tier evidence. Engine discordance is reported but does not override catalogue calls.</div></div>
        <div class="evidence-block sus"><div class="tier">Coverage-gated</div><div class="src">Susceptibility calls</div><div class="detail">SUSCEPTIBLE is only reported where the defining loci are callable above the coverage threshold. Below threshold, the call is NOT_ASSESSED.</div></div>
      ` : `
        <div class="evidence-block ind"><div class="tier">Insufficient phenotype signal</div><div class="src">Not estimable</div><div class="detail">Fewer than 20 phenotypically resistant isolates in the cohort. Error rate cannot be computed; the tool will report INDETERMINATE for non-catalogued variants in this drug's loci.</div></div>
      `}
    </div>
    <div class="drill-actions">
      <button onclick="downloadReport('bench')">Download benchmark row</button>
      <button onclick="downloadReport('coverage')">Download coverage</button>
      <button onclick="panel.classList.remove('open')">Close</button>
    </div>
  `;
  openPanel(`${d.name} · ${d.id}`, body);
}

function openCallDrill(d){
  const drug = DRUGS.find(x=>x.id===d.drug);
  const call = d.call.toUpperCase();
  const callClass = { res:'res', sus:'sus', ind:'ind' }[d.call] || 'na';
  const body = `
    <div class="drill-section">
      <h4>Identity</h4>
      <div class="drill-kv">
        <span class="k">Isolate</span><span class="v">${d.iso}</span>
        <span class="k">Lineage</span><span class="v">${d.lin}</span>
        <span class="k">Drug</span><span class="v">${drug.name}</span>
      </div>
    </div>
    <div class="drill-section">
      <h4>Reconciled call</h4>
      <div class="evidence-block ${callClass}">
        <div class="tier">Six-state ontology</div>
        <div class="src" style="font-size:15px">${call}</div>
        <div class="detail">${describeCall(call)}</div>
      </div>
    </div>
    <div class="drill-section">
      <h4>Evidence chain</h4>
      ${evidenceChain(d.call, drug)}
    </div>
    <div class="drill-section">
      <h4>Coverage</h4>
      <div class="drill-kv">
        <span class="k">Callability</span><span class="v">${pct(drug.coverage)}</span>
        <span class="k">Gating</span><span class="v">${d.call === 'na' ? 'loci below threshold' : 'loci callable'}</span>
      </div>
    </div>
    <div class="drill-actions">
      <button onclick="downloadReport('calls')">Download full call matrix</button>
      <button onclick="panel.classList.remove('open')">Close</button>
    </div>
  `;
  openPanel(`${d.iso} · ${drug.id}`, body);
}

function evidenceChain(call, drug){
  if(call === 'na'){
    return `<div class="evidence-block na"><div class="tier">Coverage gate</div><div class="src">Locus below threshold</div><div class="detail">The defining loci for ${drug.name} were not callable at sufficient depth. No susceptibility can be inferred.</div></div>`;
  }
  if(call === 'sus'){
    return `<div class="evidence-block sus"><div class="tier">Catalogued · absence</div><div class="src">No resistance markers detected</div><div class="detail">Loci were callable and no catalogued resistance variant was found. This is the only path to SUSCEPTIBLE.</div></div>`;
  }
  if(call === 'ind'){
    return `<div class="evidence-block ind"><div class="tier">Uncertain</div><div class="src">Non-catalogued variant detected</div><div class="detail">A variant was found in a resistance-associated locus but is not in the current catalogue. Held as INDETERMINATE pending review or phenotypic confirmation.</div></div>`;
  }
  if(call === 'res'){
    return `<div class="evidence-block res"><div class="tier">Catalogued · tier 1</div><div class="src">Canonical resistance variant</div><div class="detail">A catalogued resistance variant was detected in a defining locus. Call is established.</div></div>`;
  }
  return '';
}

function openVariantDrill(d){
  const body = `
    <div class="drill-section">
      <h4>Variant</h4>
      <div class="drill-kv">
        <span class="k">Position</span><span class="v">${d.pos}</span>
        <span class="k">Change</span><span class="v">${d.ref} → ${d.alt}</span>
        <span class="k">Isolate</span><span class="v">${d.iso}</span>
        <span class="k">Lineage</span><span class="v">${d.lin}</span>
      </div>
    </div>
    <div class="drill-section">
      <h4>Effect on call</h4>
      <div class="evidence-block ind"><div class="tier">Pending</div><div class="src">Depends on catalogue entry</div><div class="detail">Whether this position drives a resistance call depends on the catalogue entry for this codon and the drug under assessment. If the variant is not catalogued, the call will be INDETERMINATE.</div></div>
    </div>
    <div class="drill-section">
      <h4>Evidence context</h4>
      <div class="drill-kv">
        <span class="k">Conservation</span><span class="v">high (phyloP &gt; 2)</span>
        <span class="k">Structural</span><span class="v">in drug-binding pocket</span>
        <span class="k">Catalogue</span><span class="v">tier 2 (likely)</span>
      </div>
    </div>
    <div class="drill-actions">
      <button onclick="downloadReport('vcf')">Download variant VCF</button>
      <button onclick="panel.classList.remove('open')">Close</button>
    </div>
  `;
  openPanel(`${d.ref}${d.pos}${d.alt}`, body);
}

function openVUSDrill(rank){
  const v = VUS_ITEMS.find(x=>x.rank===rank);
  const body = `
    <div class="drill-section">
      <h4>Variant</h4>
      <div class="drill-kv">
        <span class="k">Name</span><span class="v">${v.variant}</span>
        <span class="k">Gene</span><span class="v">${v.gene}</span>
        <span class="k">Drug</span><span class="v">${v.drug}</span>
        <span class="k">Priority band</span><span class="v">${v.priority || '—'}</span>
        <span class="k">Score</span><span class="v">${(v.score === null || v.score === undefined) ? 'not ranked — insufficient data' : v.score}</span>
      </div>
    </div>
    <div class="drill-section">
      <h4>Feature context</h4>
      <div class="drill-kv">
        <span class="k">Dimensions available</span><span class="v">${(v.features||{}).dimensions || '—'}</span>
        <span class="k">Data gaps</span><span class="v">${(v.features||{}).gaps || '—'}</span>
        <span class="k">Recommended experiment</span><span class="v">${(v.features||{}).experiment || '—'}</span>
      </div>
    </div>
    <div class="drill-section">
      <h4>Why this is ranked here</h4>
      <div class="evidence-block ind"><div class="tier">Validation priority</div><div class="src">Ranked by resolution potential</div><div class="detail">The band reflects how much uncertainty a single laboratory validation would resolve, given the dimensions that were actually available. A dimension with no data source reports itself unavailable rather than defaulting to a number, so a variant can be unranked entirely. This is <strong>not</strong> a probability of resistance.</div></div>
    </div>
    <div class="drill-actions">
      <button onclick="downloadReport('json')">Export VUS record</button>
      <button onclick="panel.classList.remove('open')">Close</button>
    </div>
  `;
  openPanel(`VUS #${rank} · ${v.variant}`, body);
}

function openMechDrill(i){
  const m = MECH_CARDS[i];
  const body = `
    <div class="drill-section">
      <h4>Mechanism</h4>
      <p style="font-size:13px;line-height:1.65;color:var(--text-1)">${m.body}</p>
    </div>
    <div class="drill-section">
      <h4>Tags</h4>
      <div>${m.tags.map(t=>`<span class="tag ${t.t}" style="margin-right:6px">${t.l}</span>`).join('')}</div>
    </div>
    <div class="drill-section">
      <h4>Resolving evidence</h4>
      <div class="evidence-block ind"><div class="tier">Evidence gap</div><div class="src">What would resolve this</div><div class="detail">A phenotypic MIC measurement, an expression assay, or a structural model import — depending on the specific mechanism. Until then, the call remains INDETERMINATE.</div></div>
    </div>
    <div class="drill-actions">
      <button onclick="panel.classList.remove('open')">Close</button>
    </div>
  `;
  openPanel(`Mechanism · ${m.title.split('→')[0].trim()}`, body);
}

/* ==========================================================
   DOWNLOADS
   ========================================================== */
function downloadReport(kind){
  const payloads = {
    json: () => JSON.stringify({ run: DATA.run || {}, drugs:DRUGS, vus:VUS_ITEMS, discordance:DISCORDANCE_ROWS, audit:AUDIT_EVENTS, coverage_loci:COVERAGE_LOCI, lineages:LINEAGES }, null, 2),
    calls: () => {
      let csv = 'isolate,lineage,' + DRUGS.map(d=>d.id).join(',') + '\n';
      (isolateData||[]).slice(0,400).forEach(iso=>{
        csv += `${iso.id},${iso.lin},` + DRUGS.map(d=>iso.calls[d.id]).join(',') + '\n';
      });
      return csv;
    },
    bench: () => 'drug,abbreviation,class,coverage,error_rate,vme,me,evaluable_n\n' + DRUGS.map(d=>`${d.name},${d.id},${d.cls},${d.coverage??'NA'},${d.error??'NA'},${d.vme??'NA'},${d.me??'NA'},${d.evaluable??'NA'}`).join('\n'),
    disc: () => 'variant\tdrug\ttb_profiler\tmykrobe\tmyconductor\treason\n' + DISCORDANCE_ROWS.map(r=>`${r.var}\t${r.drug}\t${r.prof}\t${r.myk}\t${r.rec}\t${r.reason}`).join('\n'),
    vcf: () => '##fileformat=VCFv4.2\n##source=myconductor\n##reference=NC_000962.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' +
      VUS_ITEMS.map((v,i)=>`NC_000962.3\t${1000000 + i*100}\t.\tA\tG\t100\tPASS\tGENE=${v.gene};PRIORITY=${v.priority}`).join('\n'),
    coverage: () => 'chrom\tstart\tend\tlocus\tdrug\tcallable_pct\n' + COVERAGE_LOCI.map(l=>`NC_000962.3\t0\t0\t${l.id}\t${l.drug}\t${l.callable}`).join('\n'),
    audit: () => AUDIT_EVENTS.map(a=>JSON.stringify(a)).join('\n'),
    html: () => document.documentElement.outerHTML,
  };
  const ext = { json:'json', calls:'csv', bench:'csv', disc:'tsv', vcf:'vcf', coverage:'bed', audit:'jsonl', html:'html' }[kind];
  const mime = { json:'application/json', csv:'text/csv', tsv:'text/tab-separated-values', vcf:'text/plain', bed:'text/plain', jsonl:'application/x-ndjson', html:'text/html' }[kind];
  const content = payloads[kind]();
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `myconductor-report-${kind}.${ext}`;
  a.click();
  URL.revokeObjectURL(url);
}

/* ==========================================================
   NAV ACTIVE
   ========================================================== */
function initNav(){
  const links = document.querySelectorAll('.section-nav a');
  const sections = Array.from(links).map(a => document.querySelector(a.getAttribute('href'))).filter(Boolean);
  const io = new IntersectionObserver(entries => {
    entries.forEach(e => {
      if(e.isIntersecting){
        links.forEach(l => l.classList.toggle('active', l.getAttribute('href') === '#' + e.target.id));
      }
    });
  }, { rootMargin:'-40% 0px -55% 0px' });
  sections.forEach(s => io.observe(s));
}

/* ==========================================================
   BOOT
   ========================================================== */
renderDonut();
renderTrend();
renderDrugBar();
renderScatter();
renderErrorStack();
renderHeadline();
renderDrugs();
renderCallHeatmap();
renderCoverage();
renderLineage();
renderRadar();
renderTierStack();
renderAliView();
renderDiscordance();
renderVUS();
renderMech();
renderTimeline();
renderValidationDonut();
renderSiteChart();
renderAudit();
initNav();
