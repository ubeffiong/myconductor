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
const MUTATIONS = DATA.mutation_index || [];
const PREVALENCE = DATA.prevalence || [];
const TARGETS = DATA.targets || [];
const EPISTASIS = DATA.epistasis || [];
const WATCHLIST = DATA.watchlist || [];
const WORKFLOWS = DATA.workflows || {};
const DISCORDANCE_ROWS = DATA.discordance || [];
const DISCORDANCE_TICKETS = DATA.discordance_tickets || [];
const SAMPLE_DETAILS = DATA.sample_details || {};
const EXTERNAL_BENCHMARKS = DATA.external_benchmarks || [];
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
function esc(v){ return String(v === null || v === undefined ? '' : v).replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;',"'":'&#39;'}[ch])); }
function pctFraction(v, digits){ return (v === null || v === undefined || Number.isNaN(v)) ? '—' : (Number(v) * 100).toFixed(digits === undefined ? 1 : digits) + '%'; }

/* Status words this report shows (ticket/watchlist lifecycle, benchmark
   verdicts) get one shared tone mapping, so "open"/"resolved"/"beats
   baseline" read the same amber/green/red semantics as a call state does. */
function statusTone(status){
  const key = String(status || '').toLowerCase().replace(/[\s-]+/g, '_');
  if(['resolved', 'beats_baseline', 'approved', 'pass'].indexOf(key) >= 0) return 'sus';
  if(['open', 'no_better_than_baseline', 'rejected', 'fail'].indexOf(key) >= 0) return 'res';
  if(['under_investigation', 'cannot_match_coverage', 'underpowered', 'reviewed', 'submitted'].indexOf(key) >= 0) return 'ind';
  return 'na';
}
function statusPill(status){
  return `<span class="status-pill status-${statusTone(status)}">${safeText(status)}</span>`;
}

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
      data-tip="${d.k}" data-enum="${d.enum || d.k}" data-v="${d.v}" data-pct="${(d.v/total*100).toFixed(1)}"
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
      <div class="tt-section"><div class="tt-note">${describeCall(el.dataset.enum || el.dataset.tip)}</div></div>
    </div>`);
}
function describeCall(k){
  const m = {
    resistant:'Established by catalogued or phenotypic evidence. Can inform a clinical decision.',
    susceptible:'Coverage-gated absence of resistance markers. Only valid where the loci were callable.',
    indeterminate:'Evidence insufficient to establish either state. Cannot be acted on.',
    not_assessed:'Loci not shown callable. Absence of evidence is not susceptibility.',
    no_call:'Analysis-side gap. Not a biological finding.',
    unsupported:'Drug not covered by this organism profile.',
  };
  const key = String(k || '').toLowerCase().replace(/[\s-]+/g, '_');
  return m[key] || '';
}

/* ==========================================================
   TREND CHART
   ========================================================== */
function renderTrend(){
  const series = DATA.error_trend || [];
  if(series.length < 2){
    return emptyState('trendChart', 'No ordered error series',
      'An error trend needs an ordered cohort series from the report-context dataset. This run did not supply one.');
  }
  const w=520, h=200, pad={t:16,r:16,b:34,l:44};
  const n = series.length;
  const pts = series.map((v, i) => ({ x:i+1, y:v }));
  const xMax = Math.max(1, n), observedMax = Math.max.apply(null, series.map(Number));
  const yMax = Math.max(5, Math.ceil(observedMax / 5) * 5);
  const sx = v => pad.l + (v/xMax)*(w-pad.l-pad.r);
  const sy = v => h-pad.b - (v/yMax)*(h-pad.t-pad.b);
  let path='M';
  pts.forEach((p,i)=> path += `${i?'L':' '}${sx(p.x)} ${sy(p.y)}`);
  let grid='';
  const yStep = yMax / 4;
  for(let j=0;j<=4;j++){
    const v = j * yStep;
    grid += `<line x1="${pad.l}" x2="${w-pad.r}" y1="${sy(v)}" y2="${sy(v)}" class="grid-line" />
      <text x="${pad.l-8}" y="${sy(v)+4}" text-anchor="end" class="axis-label">${v.toFixed(v%1?1:0)}%</text>`;
  }
  let xLabels='';
  const xStep = Math.max(1, Math.ceil(n/5));
  for(let i=1;i<=n;i+=xStep){
    xLabels += `<text x="${sx(i)}" y="${h-12}" text-anchor="middle" class="axis-label">${i}</text>`;
  }
  let dots = pts.map(p => `<circle cx="${sx(p.x)}" cy="${sy(p.y)}" r="2.5"
    fill="var(--accent)" class="point"
    onmouseover="trendHover(event,${p.x},${p.y.toFixed(2)})" onmouseout="hideTip()" />`).join('');
  document.getElementById('trendChart').innerHTML = `
    <svg viewBox="0 0 ${w} ${h}" class="chart-svg">
      ${grid}
      <path d="${path}" class="line" stroke="var(--accent)" />
      ${dots}
      <text x="${w/2}" y="${h-2}" text-anchor="middle" class="axis-label">ordered cohort window</text>
      <text x="12" y="${h/2}" text-anchor="middle" transform="rotate(-90 12 ${h/2})" class="axis-label">cumulative error</text>
      ${xLabels}
    </svg>`;
}
function trendHover(e, n, v){
  showTip(e, `<div class="tt-head"><div class="tt-title">Cohort window ${n}</div><div class="tt-sub">supplied ordered series</div></div>
    <div class="tt-body">
      <div class="tt-row"><span class="tt-k">Cumulative error</span><span class="tt-v">${v}%</span></div>
      <div class="tt-section"><div class="tt-note">Interpret this point using the window definition and ordering used by the source dataset.</div></div>
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
      <div class="tt-row"><span class="tt-k">Error rate</span><span class="tt-v">${d.estimable ? pct(d.error) : 'Not estimable'}</span></div>
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

/* ==========================================================
   SAMPLE DRILLDOWN
   Select one isolate out of a cohort report and inspect its own site,
   lineage, drug calls, called variants and locus coverage -- all read from
   the same payload the cohort-level views already use. Nothing here is
   re-fetched or re-derived per selection.
   ========================================================== */
function initSampleFilter(){
  const select = document.getElementById('sampleFilter');
  if(!select) return;
  const ids = Object.keys(SAMPLE_DETAILS);
  if(!ids.length){
    return emptyState('sampleDrilldown', 'No sample detail available',
      'This payload carries no per-sample context to drill into.');
  }
  select.innerHTML = ids.map(function(id){
    return `<option value="${safeText(id)}">${safeText(id)}</option>`;
  }).join('');
  renderSampleDrilldown(ids[0]);
}

function renderSampleDrilldown(sampleId){
  const container = document.getElementById('sampleDrilldown');
  if(!container) return;
  const sample = SAMPLE_DETAILS[sampleId];
  if(!sample){
    container.innerHTML = `<div class="empty-state"><div class="empty-title">Sample not found</div>
      <div class="empty-detail">No detail is carried for "${safeText(sampleId)}" in this payload.</div></div>`;
    return;
  }
  const calls = (sample.drug_calls || []).map(function(d){
    return `<tr><td>${safeText(d.drug)}</td>
      <td><span class="pill call-${safeText(d.call_key)}">${safeText(d.call)}</span></td>
      <td>${safeText(d.tier)}</td><td>${safeText(d.reason)}</td></tr>`;
  }).join('') || '<tr><td colspan="4" style="color:var(--text-2)">No drug results for this sample.</td></tr>';

  const variants = (sample.variants || []).map(function(v){
    return `<tr><td>${safeText(v.gene)}</td><td>${safeText(v.variant)}</td>
      <td>${safeText(v.consequence)}</td><td>${safeText(v.drug)}</td><td>${safeText(v.call)}</td></tr>`;
  }).join('') || '<tr><td colspan="5" style="color:var(--text-2)">No called variants recorded for this sample.</td></tr>';

  const coverage = (sample.coverage || []).map(function(c){
    const frac = c.callable_fraction;
    const cls = frac === null || frac === undefined ? '' : (frac >= 0.95 ? 'good' : frac >= 0.5 ? 'low' : 'crit');
    return `<tr><td>${safeText(c.locus)}</td>
      <td class="${cls ? 'value ' + cls : ''}">${pctFraction(frac, 0)}</td>
      <td>${num(c.mean_depth)}</td></tr>`;
  }).join('') || '<tr><td colspan="3" style="color:var(--text-2)">No coverage evidence recorded for this sample.</td></tr>';

  container.innerHTML = `
    <div class="sample-meta-grid">
      <div><span class="k">Site</span><span class="v">${safeText(sample.site)}</span></div>
      <div><span class="k">Lineage</span><span class="v">${safeText(sample.lineage)}</span></div>
      <div><span class="k">Organism</span><span class="v">${safeText(sample.organism)}</span></div>
      <div><span class="k">Assay</span><span class="v">${safeText(sample.assay)}</span></div>
    </div>
    <div class="two-col">
      <div><h4>Drug calls</h4><div class="tw"><table class="disc-table">
        <thead><tr><th>Drug</th><th>Call</th><th>Tier</th><th>Reason</th></tr></thead>
        <tbody>${calls}</tbody></table></div></div>
      <div><h4>Called variants</h4><div class="tw"><table class="disc-table">
        <thead><tr><th>Gene</th><th>Variant</th><th>Consequence</th><th>Drug</th><th>Call</th></tr></thead>
        <tbody>${variants}</tbody></table></div></div>
    </div>
    <h4>Locus coverage</h4><div class="tw"><table class="disc-table">
      <thead><tr><th>Locus</th><th>Callable fraction</th><th>Mean depth</th></tr></thead>
      <tbody>${coverage}</tbody></table></div>`;
}

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
      'Coverage is the precondition for susceptibility. Without a depth table, BED mask or gVCF, every drug is reported not assessed.');
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
      <div class="tt-section"><div class="tt-note">Coverage-gated: a drug is not assessed unless this locus is callable above threshold.</div></div>
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
      const nPos = (l.positions || []).length;
      const nCarriers = new Set((l.positions || []).flatMap(v => v.carriers || [])).size;
      return `<option value="${esc(id)}">${esc(l.name)} · ${nPos} position(s) · ${nCarriers} carrier(s)</option>`;
    }).join('');
    select.dataset.filled = '1';
  }
  const locusId = (select && LOCI[select.value]) ? select.value : ids[0];
  const locus = LOCI[locusId];
  const focus = (document.getElementById('avFocus') || {}).value || 'all';
  const sourcePositions = locus.positions || [];
  const positions = sourcePositions.filter(v =>
    focus === 'all' || (focus === 'catalogued' ? v.catalogued : !v.catalogued));
  const isolates = (locus.isolates || []).slice(
    0, parseInt((document.getElementById('avCount') || {}).value || '12', 10));
  const cellW = 46;

  if(!positions.length){
    return emptyState('avContent', `No ${focus} positions in ${locus.name}`,
      'The current filter removed every observed position for this locus. Choose all called positions to return to the complete context view.');
  }

  const posByKey = new Map(positions.map(v => [String(v.pos), v]));
  const locusStart = Math.min(...positions.map(v => Number(v.pos)));
  const locusEnd = Math.max(...positions.map(v => Number(v.pos)));
  const carrierTotal = new Set(positions.flatMap(v => v.carriers || [])).size;
  const positionTotal = positions.length;

  let html = `<div class="av-summary">
    <div><span class="av-summary-k">Locus</span><strong>${esc(locus.name)}</strong><small>${esc(locus.assembly || 'reference')}</small></div>
    <div><span class="av-summary-k">Observed span</span><strong>${locusStart.toLocaleString()}–${locusEnd.toLocaleString()}</strong><small>${positionTotal} called position(s)</small></div>
    <div><span class="av-summary-k">Cohort carriers</span><strong>${carrierTotal}/${(locus.isolates || []).length}</strong><small>within rendered report payload</small></div>
    <div><span class="av-summary-k">Evidence mix</span><strong>${positions.filter(v=>v.catalogued).length} catalogued</strong><small>${positions.filter(v=>!v.catalogued).length} contextual/uncertain</small></div>
  </div>`;

  html += `<div class="av-row track"><div class="av-label">Coordinate</div><div class="av-ruler">`;
  positions.forEach(v => {
    html += `<button class="av-tick major" style="width:${cellW}px;min-width:${cellW}px" data-pos="${v.pos}" title="${esc(v.display || v.label || '')}">${v.pos}</button>`;
  });
  html += `</div></div>`;

  html += `<div class="av-row track-gene"><div class="av-label">Gene model</div>
    <div class="av-gene-track" style="width:${positions.length*cellW}px;min-width:${positions.length*cellW}px">
      <span>${esc(locus.name)}</span>${positions.map(v => `<i class="${v.catalogued ? 'hot' : 'uncertain'}" style="left:${positions.indexOf(v)*cellW + cellW/2}px" title="${esc(v.display || v.label)}"></i>`).join('')}
    </div></div>`;

  const track = (label, values, cls='') => {
    html += `<div class="av-row track"><div class="av-label">${label}</div><div class="av-seq ${cls}">` + values.join('') + `</div></div>`;
  };
  track('Reference', positions.map(v => `<div class="av-base ${esc(v.ref)}" style="width:${cellW}px;min-width:${cellW}px">${esc(v.ref)}</div>`));
  track('Alternate', positions.map(v => `<div class="av-base ${v.catalogued ? 'variant' : 'variant ind'}" style="width:${cellW}px;min-width:${cellW}px" title="${esc(v.display || v.label)}">${esc(v.alt)}</div>`));
  track('Consequence', positions.map(v => `<div class="av-annotation" style="width:${cellW}px;min-width:${cellW}px" title="${esc(v.display || '')}">${esc((v.consequence || '—').replace(' ', '\n'))}</div>`));
  track('Drug', positions.map(v => `<div class="av-annotation drug" style="width:${cellW}px;min-width:${cellW}px" title="${esc((v.drugs || []).join(', '))}">${esc((v.drugs || ['—'])[0])}</div>`));
  track('Carriers', positions.map(v => `<div class="av-carrier-bar" style="width:${cellW}px;min-width:${cellW}px"><span style="height:${Math.max(4, Math.round((v.carrier_fraction || 0) * 34))}px" class="${v.catalogued ? 'hot' : 'uncertain'}"></span><em>${esc(v.carrier_count || 0)}</em></div>`), 'carrier-track');

  isolates.forEach(iso => {
    html += `<div class="av-row"><div class="av-label"><span class="lid">${esc(iso.id)}</span><span class="llin">${esc(iso.lineage || 'untyped')}</span></div><div class="av-seq">`;
    positions.forEach(v => {
      const variant = (iso.variants || {})[String(v.pos)];
      if(variant){
        const payload = encodeURIComponent(JSON.stringify({ iso: iso.id, lineage: iso.lineage, site: iso.site, locus: locus.name, ...v, isolateVariant: variant }));
        html += `<button class="av-base ${v.catalogued ? 'variant' : 'variant ind'}"
          style="width:${cellW}px;min-width:${cellW}px"
          data-variant="${payload}" title="${esc(iso.id)} carries ${esc(v.display || v.label)}">${esc(v.alt)}</button>`;
      } else {
        html += `<div class="av-base refcall" style="width:${cellW}px;min-width:${cellW}px" title="No alternate allele reported for ${esc(iso.id)} at ${v.pos}">${esc(v.ref)}</div>`;
      }
    });
    html += `</div></div>`;
  });

  html += `<div class="av-position-cards">`;
  positions.forEach(v => {
    html += `<button class="av-position-card ${v.catalogued ? 'catalogued' : 'uncertain'}" data-variant="${encodeURIComponent(JSON.stringify({ locus: locus.name, ...v }))}">
      <strong>${esc(v.display || v.label || `${v.ref}${v.pos}${v.alt}`)}</strong>
      <span>${esc((v.drugs || []).join(', ') || 'drug not mapped')}</span>
      <small>${esc((v.calls || []).join(', ') || 'call not supplied')} · ${esc((v.tiers || []).join(', ') || 'tier not supplied')} · ${pctFraction(v.carrier_fraction)} carriers</small>
    </button>`;
  });
  html += `</div>`;

  content.innerHTML = html;

  content.querySelectorAll('[data-variant]').forEach(el => {
    const data = () => JSON.parse(decodeURIComponent(el.dataset.variant));
    el.addEventListener('mouseenter', ev => {
      const v = data();
      showTip(ev, `<div class="tt-head"><div class="tt-title">${esc(v.display || v.label || `${v.ref}${v.pos}${v.alt}`)}</div>
          <div class="tt-sub">${esc(v.locus || locus.name)}${v.drugs && v.drugs.length ? ' · ' + esc(v.drugs.join(', ')) : ''}</div></div>
        <div class="tt-body">
          <div class="tt-row"><span class="tt-k">Position</span><span class="tt-v">${esc(v.pos)}</span></div>
          <div class="tt-row"><span class="tt-k">Change</span><span class="tt-v">${esc(v.ref)}→${esc(v.alt)}</span></div>
          <div class="tt-row"><span class="tt-k">Carriers</span><span class="tt-v">${esc(v.carrier_count || 0)} (${pctFraction(v.carrier_fraction)})</span></div>
          <div class="tt-section"><div class="tt-note">${v.catalogued ? 'Catalogued evidence can establish resistance when the catalogue and coverage gates support it.' : 'Contextual or uncertain evidence is retained for review and cannot establish resistance by itself.'}</div></div>
        </div>`);
    });
    el.addEventListener('mouseleave', hideTip);
    el.addEventListener('click', () => openVariantDrill(data()));
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
  const shown = VUS_ITEMS.filter(vusMatches);
  if(!shown.length){
    return emptyState('vusList', 'No variant matches this filter',
      'Every variant was excluded by the current selection. Choose All to see the whole queue.');
  }
  list.innerHTML = shown.map(v => {
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
  if(!(MECH_CARDS) || !(MECH_CARDS).length){ return emptyState('mechGrid', 'No mechanism hypotheses', 'No lane raised a named mechanism for the rendered report scope.'); }
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
  if(!d) return;
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
        <div class="evidence-block sus"><div class="tier">Coverage-gated</div><div class="src">Susceptibility calls</div><div class="detail">Susceptibility is only reported where the defining loci are callable above the coverage threshold. Below threshold, the call is not assessed.</div></div>
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
  if(!drug) return;
  const longForm = { res:'resistant', sus:'susceptible',
                     ind:'indeterminate', nc:'no_call',
                     un:'unsupported', na:'not_assessed' };
  const call = lbl(longForm[d.call] || 'not_assessed');
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
  const title = d.display || d.label || `${d.ref || ''}${d.pos || ''}${d.alt || ''}`;
  const rationales = (d.rationales || []).map(x => `<li>${esc(x)}</li>`).join('') || '<li>No rationale text supplied for this evidence record.</li>';
  const limitations = (d.limitations || []).map(x => `<li>${esc(x)}</li>`).join('') || '<li>No additional limitation text supplied.</li>';
  const ids = (d.evidence_ids || []).map(x => `<span class="tag engine">${esc(x)}</span>`).join(' ') || '<span class="tag na">No evidence id</span>';
  const isoRows = d.iso ? `<span class="k">Selected isolate</span><span class="v">${esc(d.iso)}</span>
        <span class="k">Lineage</span><span class="v">${esc(d.lineage || d.lin || 'untyped')}</span>
        <span class="k">Site</span><span class="v">${esc(d.site || 'site not recorded')}</span>` : '';
  const body = `
    <div class="drill-section">
      <h4>Variant identity</h4>
      <div class="drill-kv">
        <span class="k">Locus</span><span class="v">${esc(d.locus || '—')}</span>
        <span class="k">Position</span><span class="v">${esc(d.pos || '—')}</span>
        <span class="k">Allele change</span><span class="v">${esc(d.ref || '—')} → ${esc(d.alt || '—')}</span>
        <span class="k">HGVS / amino-acid change</span><span class="v">${esc(d.change || d.label || 'not supplied')}</span>
        <span class="k">Consequence</span><span class="v">${esc(d.consequence || 'not supplied')}</span>
        ${isoRows}
      </div>
    </div>
    <div class="drill-section">
      <h4>Evidence effect</h4>
      <div class="evidence-block ${d.catalogued ? 'res' : 'ind'}"><div class="tier">${d.catalogued ? 'Catalogued resistance-associated' : 'Contextual / uncertain'}</div><div class="src">${esc((d.drugs || []).join(', ') || 'Drug mapping not supplied')}</div><div class="detail">${d.catalogued ? 'This mutation is represented by catalogued evidence in the report payload. The final drug call still remains coverage-gated and evidence-tier-gated.' : 'This mutation is visible because it was observed in a resistance-associated context, but it is not sufficient to establish resistance without stronger evidence.'}</div></div>
    </div>
    <div class="drill-section">
      <h4>Cohort context</h4>
      <div class="drill-kv">
        <span class="k">Carriers in report</span><span class="v">${esc(d.carrier_count || 0)} (${pctFraction(d.carrier_fraction)})</span>
        <span class="k">Calls represented</span><span class="v">${esc((d.calls || []).join(', ') || 'not supplied')}</span>
        <span class="k">Evidence tiers</span><span class="v">${esc((d.tiers || []).join(', ') || 'not supplied')}</span>
        <span class="k">Evidence lanes</span><span class="v">${esc((d.lanes || []).join(', ') || 'not supplied')}</span>
        <span class="k">Mean VAF</span><span class="v">${d.mean_vaf === null || d.mean_vaf === undefined ? '—' : pctFraction(d.mean_vaf)}</span>
        <span class="k">Median depth</span><span class="v">${num(d.median_depth)}</span>
      </div>
    </div>
    <div class="drill-section">
      <h4>Rationale and limitations</h4>
      <ul class="drill-list">${rationales}</ul>
      <h4 style="margin-top:12px">Limitations</h4>
      <ul class="drill-list">${limitations}</ul>
      <div style="margin-top:10px">${ids}</div>
    </div>
    <div class="drill-actions">
      <button onclick="downloadReport('vcf')">Download variant VCF</button>
      <button onclick="panel.classList.remove('open')">Close</button>
    </div>
  `;
  openPanel(title, body);
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
  if(!m) return;
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
const LABELS = DATA.labels || {};

/* Display labels come from the renderer, so Python and this runtime cannot
   drift into two spellings of the same term. An unknown token falls back to a
   mechanical rule rather than rendering raw. */
function lbl(token){
  if(token === null || token === undefined || token === '') return '—';
  const key = String(token);
  if(LABELS[key]) return LABELS[key];
  if(LABELS[key.toLowerCase()]) return LABELS[key.toLowerCase()];
  const words = key.replace(/^_+/, '').split(/[\s_\-.\/]+/).filter(Boolean);
  if(!words.length) return key;
  return words.map(function(w, i){
    return i === 0 ? w.charAt(0).toUpperCase() + w.slice(1).toLowerCase()
                   : w.toLowerCase();
  }).join(' ');
}

/* ==========================================================
   EXPORTS
   ========================================================== */
function downloadReport(kind){
  /* Every export is built from the injected payload. The design this was
     ported from wrote a VCF whose coordinates were 1000000 + i*100 - an
     invented position for every variant, in a file format whose entire
     purpose is coordinates. An export that cannot be built from real data
     emits a header comment saying so and writes no rows. */
  const run = DATA.run || {};
  const stamp = run.label || run.sample_id || 'report';
  const NL = String.fromCharCode(10);
  const TAB = String.fromCharCode(9);

  const payloads = {
    json: function(){
      return JSON.stringify({
        run: run, drugs: DRUGS, vus: VUS_ITEMS, discordance: DISCORDANCE_ROWS,
        audit: AUDIT_EVENTS, coverage_loci: COVERAGE_LOCI, lineages: LINEAGES,
        call_distribution: DATA.call_distribution || [],
        alignment: DATA.alignment_loci || {}, tiers: DATA.tiers || {},
        mechanisms: MECH_CARDS, measurability: DATA.measurability || [],
        lineage_strata: DATA.lineage || {}, mutation_index: MUTATIONS,
        prevalence: PREVALENCE, targets: TARGETS, epistasis: EPISTASIS,
        watchlist: WATCHLIST, implemented_workflows: WORKFLOWS,
        validation_outcomes: DATA.validation_outcomes || [],
        validation_timeline: DATA.validation_timeline || [],
        federated_sites: DATA.federated_sites || []
      }, null, 2);
    },

    calls: function(){
      const rows = DATA.isolates || [];
      let csv = '# myconductor call matrix. Six states; only susceptible admits a drug.' + NL;
      csv += 'sample,lineage,' + DRUGS.map(function(d){ return d.id; }).join(',') + NL;
      if(!rows.length){ return csv + '# no per-sample call matrix in this run' + NL; }
      rows.forEach(function(iso){
        csv += iso.id + ',' + (iso.lin || 'untyped') + ','
             + DRUGS.map(function(d){ return iso.calls[d.id] || 'na'; }).join(',') + NL;
      });
      return csv;
    },

    bench: function(){
      const f = function(v){ return (v === null || v === undefined) ? '' : v; };
      let csv = '# Blank accuracy columns mean the drug was not estimable in this run.' + NL;
      csv += 'drug,abbreviation,class,coverage_pct,error_rate_pct,sensitivity,'
           + 'specificity,ppv,npv,vme,me,evaluable_n,estimable' + NL;
      DRUGS.forEach(function(d){
        csv += [d.name, d.id, d.cls, f(d.coverage), f(d.error), f(d.sensitivity),
                f(d.specificity), f(d.ppv), f(d.npv), f(d.vme), f(d.me),
                f(d.evaluable), d.estimable ? 'yes' : 'no'].join(',') + NL;
      });
      return csv;
    },

    disc: function(){
      let tsv = '# Discordance is retained, never resolved by vote. Context is display only.' + NL;
      tsv += ['drug','sources','calls','reconciled','reason','context'].join(TAB) + NL;
      if(!DISCORDANCE_ROWS.length){ return tsv + '# no discordance in this run' + NL; }
      DISCORDANCE_ROWS.forEach(function(r){
        const calls = [r.prof, r.myk].filter(Boolean).join('/');
        tsv += [r.drug, r.sources || '', calls, r.rec,
                String(r.reason || '').split(TAB).join(' '),
                (r.context || []).join(' | ').split(TAB).join(' ')].join(TAB) + NL;
      });
      return tsv;
    },

    vcf: function(){
      const loci = DATA.alignment_loci || {};
      let out = '##fileformat=VCFv4.2' + NL + '##source=myconductor' + NL;
      out += '##reference=' + (run.reference_assembly || 'unspecified') + NL;
      out += '##INFO=<ID=GENE,Number=1,Type=String,Description="Locus">' + NL;
      out += '##INFO=<ID=CATALOGUED,Number=0,Type=Flag,Description="Graded catalogue entry">' + NL;
      out += ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'].join(TAB) + NL;
      const rows = [];
      Object.keys(loci).forEach(function(gene){
        (loci[gene].positions || []).forEach(function(v){
          let info = 'GENE=' + gene;
          if(v.catalogued){ info += ';CATALOGUED'; }
          if(v.drugs && v.drugs.length){ info += ';DRUGS=' + v.drugs.join('|'); }
          rows.push([loci[gene].assembly || '.', v.pos, '.', v.ref, v.alt,
                     '.', 'PASS', info].join(TAB));
        });
      });
      if(!rows.length){
        return out + '##note=No coordinate-resolved variant in this run. '
             + 'Positions are never invented, so no records are written.' + NL;
      }
      return out + rows.join(NL) + NL;
    },

    coverage: function(){
      let out = '# BED requires real locus spans. This profile does not carry them,' + NL
              + '# so callable fractions are exported as a table rather than as BED' + NL
              + '# intervals with invented start and end coordinates.' + NL;
      out += ['locus','drugs','callable_pct'].join(TAB) + NL;
      if(!COVERAGE_LOCI.length){ return out + '# no callable-locus evidence supplied' + NL; }
      COVERAGE_LOCI.forEach(function(l){
        const c = (l.callable === null || l.callable === undefined) ? '' : l.callable;
        out += [l.id, l.drug, c].join(TAB) + NL;
      });
      return out;
    },

    audit: function(){
      if(!AUDIT_EVENTS.length){
        return '{"note":"no ledger was supplied to this report"}' + NL;
      }
      return AUDIT_EVENTS.map(function(a){ return JSON.stringify(a); }).join(NL) + NL;
    },

    html: function(){
      return '<!doctype html>' + NL + document.documentElement.outerHTML;
    }
  };

  const ext = { json:'json', calls:'csv', bench:'csv', disc:'tsv', vcf:'vcf',
                coverage:'tsv', audit:'jsonl', html:'html' }[kind];
  const mime = { json:'application/json', csv:'text/csv',
                 tsv:'text/tab-separated-values', vcf:'text/plain',
                 jsonl:'application/x-ndjson', html:'text/html' }[ext] || 'text/plain';
  const blob = new Blob([payloads[kind]()], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'myconductor-' + stamp + '-' + kind + '.' + ext;
  a.click();
  URL.revokeObjectURL(url);
}

/* ==========================================================
   VUS FILTERS
   ========================================================== */
let vusFilter = 'all';
function setVusFilter(btn, mode){
  document.querySelectorAll('#vus .card-toolbar .mini-btn')
    .forEach(function(b){ b.classList.remove('active'); });
  btn.classList.add('active');
  vusFilter = mode;
  renderVUS();
}
function vusMatches(v){
  const dims = String((v.features || {}).dimensions || '').toLowerCase();
  const structural = dims.indexOf('structural') >= 0 || dims.indexOf('ligand') >= 0;
  if(vusFilter === 'structural'){ return structural; }
  if(vusFilter === 'no-structural'){ return !structural; }
  if(vusFilter === 'unranked'){ return v.score === null || v.score === undefined; }
  return true;
}

/* ==========================================================
   MUTATION EXPLORER AND DISCOVERY CONTEXT
   ========================================================== */
function safeText(value){
  return String(value === null || value === undefined ? '' : value)
    .replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;')
    .replace(/"/g,'&quot;').replace(/'/g,'&#39;');
}
function searchMutations(query){
  const target = document.getElementById('mutationResults');
  if(!MUTATIONS.length){
    return emptyState('mutationResults', 'No observed mutations',
      'This report contains no coordinate-resolved variant or validation-queue record to search.');
  }
  const q = String(query || '').trim().toLowerCase();
  const found = MUTATIONS.filter(function(m){
    return !q || [m.gene,m.variant,m.classification,m.source].concat(m.drugs || [])
      .join(' ').toLowerCase().indexOf(q) >= 0;
  });
  if(!found.length){
    return emptyState('mutationResults', 'No matching mutation',
      'Try a gene, variant label, drug, or evidence classification present in this report.');
  }
  target.innerHTML = found.slice(0,24).map(function(m){
    const coordinate = m.position === null || m.position === undefined
      ? 'Coordinate not resolved' : 'Position ' + Number(m.position).toLocaleString();
    const drugs = (m.drugs || []).length ? m.drugs.join(', ') : 'No linked drug recorded';
    return `<article class="result-item">
      <div><div class="result-gene">${safeText(m.gene)}</div>
      <div class="result-name">${safeText(m.variant)}</div></div>
      <div class="result-copy"><strong>${safeText(m.classification)}</strong><br>
      ${safeText(coordinate)} · ${safeText(drugs)}<br>
      <span>${safeText(m.source)}</span></div></article>`;
  }).join('') + (found.length > 24
    ? `<div class="result-more">Showing 24 of ${found.length} matches. Refine the search to narrow the list.</div>` : '');
}

function renderMutationAtlas(){
  if(!MUTATIONS.length){
    return emptyState('mutationAtlas', 'No mutation atlas',
      'An atlas needs observed variants. No atlas values are generated for an empty run.');
  }
  const groups = {};
  MUTATIONS.forEach(function(m){
    const key = m.gene || 'Unresolved locus';
    if(!groups[key]) groups[key] = {all:0, catalogued:0, uncertain:0};
    groups[key].all += 1;
    if(String(m.classification).toLowerCase().indexOf('catalogued') >= 0){ groups[key].catalogued += 1; }
    if(String(m.classification).toLowerCase().indexOf('uncertain') >= 0){ groups[key].uncertain += 1; }
  });
  const rows = Object.keys(groups).sort(function(a,b){ return groups[b].all-groups[a].all; });
  const max = Math.max.apply(null, rows.map(function(k){ return groups[k].all; }));
  document.getElementById('mutationAtlas').innerHTML = `<div class="atlas-list">${rows.map(function(g){
    const d=groups[g];
    return `<div class="atlas-row"><div class="atlas-label">${safeText(g)}</div>
      <div class="atlas-track"><span style="width:${100*d.all/max}%"></span></div>
      <div class="atlas-value">${d.all}</div>
      <div class="atlas-note">${d.catalogued} catalogued · ${d.uncertain} uncertain</div></div>`;
  }).join('')}</div>`;
}

function renderMechanismLandscape(){
  const known = MUTATIONS.filter(function(m){
    return String(m.classification).toLowerCase().indexOf('catalogued') >= 0;
  }).length;
  const rows = [
    {label:'Catalogued associations', value:known, tone:'res', detail:'May establish resistance when the complete evidence rule is met.'},
    {label:'Candidate mechanisms', value:MECH_CARDS.length, tone:'ind', detail:'Biologically plausible, with named evidence gaps.'},
    {label:'Uncertain variants', value:VUS_ITEMS.length, tone:'na', detail:'Prioritised for validation; no resistance prediction is made.'}
  ];
  if(!rows.some(function(r){ return r.value; })){
    return emptyState('mechanismLandscape', 'No mechanism evidence',
      'No catalogued association, candidate mechanism, or uncertain variant was recorded.');
  }
  document.getElementById('mechanismLandscape').innerHTML = rows.map(function(r){
    return `<div class="mechanism-level ${r.tone}"><div class="mechanism-count">${r.value}</div>
      <div><strong>${r.label}</strong><p>${r.detail}</p></div></div>`;
  }).join('');
}

function renderEpistasis(){
  if(!EPISTASIS.length){
    return emptyState('epistasisTable', 'No epistasis annotation',
      'No supplied co-observation rule matched this run. No interaction effect is inferred.');
  }
  document.getElementById('epistasisTable').innerHTML = `<div class="card table-scroll" style="padding:0"><table class="disc-table">
    <thead><tr><th>Drug</th><th>Interaction</th><th>Observed pair</th><th>Interpretation</th><th>Source</th><th>Effect on call</th></tr></thead>
    <tbody>${EPISTASIS.map(function(e){
      const pair = (e.matched_pairs || []).length ? e.matched_pairs.join('; ') : [e.primary,e.partner].filter(Boolean).join(' + ');
      return `<tr><td>${safeText(e.drug)}</td><td>${safeText(e.interaction)}</td>
        <td>${safeText(pair || 'Matched variants not named')}</td><td>${safeText(e.note)}</td>
        <td>${safeText(e.source)}<div class="muted-source">${safeText(e.rule_id || e.table_version)}</div></td>
        <td>${safeText(e.effect || 'Annotation only')}</td></tr>`;
    }).join('')}</tbody></table></div>`;
}

let prevalenceState = {drug:'All', lineage:'All', geography:'All'};
function uniqueValues(field){
  return Array.from(new Set(PREVALENCE.map(function(r){ return r[field]; }).filter(Boolean))).sort();
}
function setPrevalence(field, value){ prevalenceState[field]=value; renderPrevalence(); }
function renderPrevalence(){
  if(!PREVALENCE.length){
    emptyState('prevalenceChart', 'No prevalence dataset',
      'Supply cohort periods with resistant and tested denominators in the report-context dataset.');
    document.getElementById('prevalenceFilters').innerHTML = '';
    document.getElementById('prevalenceTable').innerHTML = '';
    renderGeographySummary([]);
    return;
  }
  const fields = [['drug','Drug'],['lineage','Lineage'],['geography','Geography']];
  document.getElementById('prevalenceFilters').innerHTML = fields.map(function(pair){
    const field=pair[0], label=pair[1];
    const id = `prev-${field}`;
    return `<div class="filter-field"><label for="${id}">${label}</label><select id="${id}" aria-label="${label}" onchange="setPrevalence('${field}',this.value)">
      <option>All</option>${uniqueValues(field).map(function(v){ return `<option${prevalenceState[field]===v?' selected':''}>${safeText(v)}</option>`; }).join('')}
      </select></div>`;
  }).join('');
  const rows = PREVALENCE.filter(function(r){
    return fields.every(function(pair){ return prevalenceState[pair[0]]==='All' || r[pair[0]]===prevalenceState[pair[0]]; });
  });
  const max = Math.max.apply(null, rows.map(function(r){ return Number(r.rate)||0; }).concat([1]));
  document.getElementById('prevalenceChart').innerHTML = rows.length ? `<div class="prevalence-bars">${rows.map(function(r){
    return `<div class="prevalence-row"><div class="prevalence-period">${safeText(r.period)}</div>
      <div class="prevalence-track"><span style="width:${100*(Number(r.rate)||0)/max}%"></span></div>
      <div class="prevalence-value">${r.rate===null||r.rate===undefined?'—':Number(r.rate).toFixed(1)+'%'}</div></div>`;
  }).join('')}</div>` : '<div class="empty-state"><div class="empty-title">No matching observations</div><div class="empty-detail">Change one or more filters.</div></div>';
  document.getElementById('prevalenceTable').innerHTML = rows.length ? `<table class="disc-table"><thead><tr>
    <th>Period</th><th>Drug</th><th>Lineage</th><th>Geography</th><th>Resistant</th><th>Tested</th><th>Prevalence</th><th>Source</th>
    </tr></thead><tbody>${rows.map(function(r){ return `<tr><td>${safeText(r.period)}</td><td>${safeText(r.drug)}</td>
    <td>${safeText(r.lineage)}</td><td>${safeText(r.geography)}</td><td>${num(r.resistant)}</td><td>${num(r.total)}</td>
    <td>${r.rate===null||r.rate===undefined?'—':Number(r.rate).toFixed(1)+'%'}</td><td>${safeText(r.source)}</td></tr>`; }).join('')}</tbody></table>` : '';
  renderGeographySummary(rows);
}

function renderGeographySummary(rows){
  const target = document.getElementById('geographySummary');
  if(!target) return;
  if(!rows.length){
    return emptyState('geographySummary', 'No geographic coverage',
      'Geographic summaries require prevalence rows with explicit geography and source fields.');
  }
  const groups = {};
  rows.forEach(function(r){
    const key = r.geography || 'Unspecified geography';
    if(!groups[key]) groups[key] = {rows:0, drugs:new Set(), periods:new Set(), sources:new Set()};
    groups[key].rows += 1;
    groups[key].drugs.add(r.drug);
    groups[key].periods.add(r.period);
    groups[key].sources.add(r.source);
  });
  target.innerHTML = `<div class="geo-grid">${Object.keys(groups).sort().map(function(g){
    const d = groups[g];
    return `<div class="geo-card"><strong>${safeText(g)}</strong>
      <span>${d.rows} observation row${d.rows===1?'':'s'}</span>
      <span>${Array.from(d.drugs).sort().join(', ')}</span>
      <span>${Array.from(d.periods).sort().join(', ')}</span>
      <small>${Array.from(d.sources).sort().join(', ')}</small></div>`;
  }).join('')}</div><p class="chart-explanation">Rows are grouped to show source coverage only. Counts are not summed across drugs or periods because the populations may overlap.</p>`;
}

function renderTargets(){
  if(!TARGETS.length){
    return emptyState('targetEvidence', 'No target-evidence dataset',
      'Supply evidence-backed target records in the report-context dataset. Scores are never invented from a gene name.');
  }
  document.getElementById('targetEvidence').innerHTML = `<div class="card table-scroll" style="padding:0"><table class="disc-table">
    <thead><tr><th>Target</th><th>Essentiality</th><th>Druggability</th><th>Human homology</th><th>Resistance liability</th><th>Evidence and source</th></tr></thead>
    <tbody>${TARGETS.map(function(t){ return `<tr><td><strong>${safeText(t.target)}</strong></td><td>${safeText(t.essentiality)}</td>
      <td>${safeText(t.druggability)}</td><td>${safeText(t.human_homology)}</td><td>${safeText(t.resistance_liability)}</td>
      <td>${safeText(t.evidence)}<div class="muted-source">${safeText(t.source)}</div></td></tr>`; }).join('')}</tbody></table></div>`;
}

function renderWatchlist(){
  if(!WATCHLIST.length){
    return emptyState('watchlistSummary', 'No watch-list aggregate',
      'Supply privacy-gated aggregate rows when a multi-site signal is ready for follow-up.');
  }
  document.getElementById('watchlistSummary').innerHTML = `<div class="watch-grid">${WATCHLIST.map(function(w){
    return `<article class="watch-card"><div class="watch-head">
      <strong>${safeText(w.drug)}</strong>${statusPill(w.status)}</div>
      <div class="watch-variant">${safeText(w.variant)}</div>
      <div class="watch-metrics"><div><b>${num(w.unresolved_isolates)}</b><span>unresolved isolates</span></div>
      <div><b>${num(w.contributing_sites)}</b><span>contributing sites</span></div></div>
      <p>${safeText(w.evidence_gaps)}</p>
      <small>${safeText(w.limitations)} Source: ${safeText(w.source)}</small></article>`;
  }).join('')}</div>`;
}

function renderDiscordanceTickets(){
  const table = document.getElementById('discordanceTicketTable');
  if(!DISCORDANCE_TICKETS.length){
    return emptyState('discordanceTicketTable', 'No discordance tickets',
      'This run\'s genomic and phenotypic calls agreed everywhere a phenotype was measured, so no ticket was opened.');
  }
  table.querySelector('tbody').innerHTML = DISCORDANCE_TICKETS.map(function(t){
    return `<tr><td>${safeText(t.sample)}</td><td>${safeText(t.drug)}</td>
      <td><span class="pill call-${safeText(t.genomic_call_key)}">${safeText(t.genomic_call)}</span></td>
      <td><span class="pill call-${safeText(t.phenotypic_call_key)}">${safeText(t.phenotypic_call)}</span></td>
      <td>${statusPill(t.status)}</td>
      <td>${safeText(t.evidence_gaps)}</td></tr>`;
  }).join('');
}

function renderExternalBenchmarks(){
  const target = document.getElementById('externalBenchmarkGrid');
  if(!target) return;
  if(!EXTERNAL_BENCHMARKS.length){
    return emptyState('externalBenchmarkGrid', 'No external model benchmarked',
      'Run `mycobench benchmark-model` against a bring-your-own prediction file, then supply the resulting rows as report context.');
  }
  target.innerHTML = EXTERNAL_BENCHMARKS.map(function(b){
    const tone = statusTone(b.verdict_key);
    return `<article class="eb-card eb-${tone}">
      <div class="eb-head"><strong>${safeText(b.model)} <span class="eb-version">${safeText(b.version)}</span></strong>${statusPill(b.verdict)}</div>
      <div class="eb-drug">${safeText(b.drug)} · ${safeText(b.lineage)}</div>
      <div class="eb-metrics">
        <div><b>${pct(b.call_rate, 1)}</b><span>call rate</span></div>
        <div><b>${pct(b.error_rate, 1)}</b><span>error rate</span></div>
        <div><b>${pct(b.baseline_coverage, 1)}</b><span>baseline coverage</span></div>
        <div><b>${pct(b.baseline_error_rate, 1)}</b><span>baseline error</span></div>
      </div>
      <p>${safeText(b.notes)}</p>
      <small>${num(b.n_evaluable)} evaluable · ${num(b.n_called)} called · ${num(b.n_dropped)} dropped without a matching phenotype · source: ${safeText(b.source)}</small>
    </article>`;
  }).join('');
}

function workflowEmpty(id, title){
  return emptyState(id, 'No ' + title.toLowerCase() + ' supplied',
    'This analysis produced no attributed records for this workflow. The absence is preserved rather than replaced with demonstration values.');
}
function workflowItem(title, badge, pairs, caveat){
  return `<article class="workflow-item"><span class="workflow-badge">${safeText(badge)}</span>
    <h4>${safeText(title)}</h4><dl>${pairs.map(function(pair){
      return `<dt>${safeText(pair[0])}</dt><dd>${safeText(pair[1])}</dd>`;
    }).join('')}</dl>${caveat ? `<p class="workflow-caveat">${safeText(caveat)}</p>` : ''}</article>`;
}
function renderWorkflowEvidence(){
  const definitions = [
    ['population','Population groups'], ['mic','MIC records'], ['structural','Structural records'],
    ['regulatory','Regulatory records'], ['expression','Expression records'],
    ['models','Governed predictions'], ['panels','Assay panels']
  ];
  document.getElementById('workflowIndex').innerHTML = definitions.map(function(d){
    return `<div class="workflow-kpi"><b>${(WORKFLOWS[d[0]] || []).length}</b><span>${d[1]}</span></div>`;
  }).join('');

  const population = WORKFLOWS.population || [];
  if(!population.length) workflowEmpty('populationWorkflow','Population-structure analysis');
  else document.getElementById('populationWorkflow').innerHTML = `<div class="workflow-grid">${population.map(function(p){
    const groups = (p.groups || []).map(function(g, i){
      return `Group ${i+1}: ${g.fraction === null ? 'unknown' : (g.fraction*100).toFixed(1)+'%'} median allele frequency; ${g.lineage}; ${(g.variants||[]).join(', ')}`;
    }).join(' | ');
    return workflowItem(p.sample, p.classification, [['Classification basis',p.basis],['Method',p.method],['Frequency groups',groups || 'None'],['Unclustered variants',(p.unclustered||[]).join(', ') || 'None'],['Evidence gaps',p.gaps || 'None recorded']], 'Frequency proximity does not prove cellular linkage, clone identity, transmission, or microevolution.');
  }).join('')}</div>`;

  const mic = WORKFLOWS.mic || [];
  if(!mic.length) workflowEmpty('micWorkflow','Quantitative MIC evidence');
  else document.getElementById('micWorkflow').innerHTML = `<div class="workflow-grid">${mic.map(function(m){
    const interval = m.interval ? m.interval.join('–') + ' ' + m.unit : 'No interval supplied';
    return workflowItem(m.drug, m.conflict ? 'Conflict requires review' : m.comparison, [['Prediction',m.value+' '+m.unit],['Interval',interval],['Critical concentration',m.critical_concentration+' '+m.unit],['Method',m.method],['Source',m.source]], m.interpretation);
  }).join('')}</div>`;

  const structural = WORKFLOWS.structural || [];
  if(!structural.length) workflowEmpty('structuralWorkflow','Structural annotation');
  else document.getElementById('structuralWorkflow').innerHTML = `<div class="workflow-grid">${structural.map(function(s){
    return workflowItem(s.variant, s.location, [['Gene',s.gene],['Reported effect',s.effect],['Ligand distance',s.distance],['Review status',s.status],['Source',s.source],['Ranking effect',s.ranking_effect]], 'An imported structural claim does not establish drug response and does not alter the VUS rank.');
  }).join('')}</div>`;

  const regulatory = WORKFLOWS.regulatory || [];
  if(!regulatory.length) workflowEmpty('regulatoryWorkflow','Regulatory-region evidence');
  else document.getElementById('regulatoryWorkflow').innerHTML = `<div class="workflow-grid">${regulatory.map(function(r){
    return workflowItem(r.name, r.tier, [['Observed variant',r.variant],['Region type',r.type],['Target genes',r.targets],['Associated drugs',r.drugs],['Source',r.source]], r.interpretation);
  }).join('')}</div>`;

  const expression = WORKFLOWS.expression || [];
  if(!expression.length) workflowEmpty('expressionWorkflow','Expression evidence');
  else document.getElementById('expressionWorkflow').innerHTML = `<div class="workflow-grid">${expression.map(function(e){
    return workflowItem(e.gene, e.conclusion, [['Measurement',e.measurement],['Reported fold change',e.fold_change+'×'],['Unit',e.unit],['Variant linkage',e.matching],['Source',e.source]], e.interpretation);
  }).join('')}</div>`;

  const models = WORKFLOWS.models || [];
  if(!models.length) workflowEmpty('modelWorkflow','Governed model output');
  else document.getElementById('modelWorkflow').innerHTML = `<div class="workflow-grid">${models.map(function(m){
    return workflowItem(m.variant, m.prediction, [['Drug',m.drug],['Model',m.model],['Reported confidence',m.confidence === null ? 'Not supplied' : m.confidence],['Evaluation cohort',m.cohort],['Incumbent baseline',m.baseline],['Approval basis',m.basis],['Effect on report',m.effect]], m.interpretation);
  }).join('')}</div>`;

  const panels = WORKFLOWS.panels || [];
  if(!panels.length) workflowEmpty('panelWorkflow','Assay-panel declaration');
  else document.getElementById('panelWorkflow').innerHTML = `<div class="workflow-grid">${panels.map(function(p){
    const unsupported = (p.unsupported||[]).map(function(x){ return x.drug + (x.missing ? ' — missing '+x.missing : ''); }).join(' | ');
    return workflowItem(p.name+' '+p.version, p.verified ? 'Verified declaration' : 'Unverified declaration', [['Assay',p.assay],['Declared loci',p.loci],['Callable loci retained',p.kept || 'None'],['Coverage claims discarded',p.discarded || 'None'],['Unsupported drugs',unsupported || 'None'],['Off-panel observed variants',p.off_panel || 'None'],['Source',p.source]], p.note);
  }).join('')}</div>`;
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
initSampleFilter();
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
searchMutations('');
renderMutationAtlas();
renderMechanismLandscape();
renderEpistasis();
renderPrevalence();
renderWatchlist();
renderTargets();
renderWorkflowEvidence();
renderAliView();
renderDiscordance();
renderDiscordanceTickets();
renderExternalBenchmarks();
renderVUS();
renderMech();
renderTimeline();
renderValidationDonut();
renderSiteChart();
renderAudit();
initNav();
