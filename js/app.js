/* =========================================================
   app.js — user interface for the Stalling Profile Explorer
   ========================================================= */
(function () {
  'use strict';

  const SD = window.SPEData;
  const PALETTE = ['#7b6cd9', '#e8836b', '#3fa88d', '#5b9bd5', '#d9a23b', '#c46fa0',
                   '#7fa65a', '#4c4f8f', '#e4a3b8', '#2f8a9e', '#a58a6b', '#8a8fa8'];
  const MAX_FACETS = 12;
  const PAGE = 25;
  const TABLE_META_PREF = ['motif.type', 'STR.motif.ssDNA', 'motif.length', 'pred.struct', 'periodicity'];
  const collator = new Intl.Collator('en', { numeric: true, sensitivity: 'base' });

  // ---------- helpers ----------
  const $ = (sel, el = document) => el.querySelector(sel);
  const $$ = (sel, el = document) => Array.from(el.querySelectorAll(sel));
  const fmtInt = n => Number(n).toLocaleString('en-GB', { maximumFractionDigits: 0 });
  const fmt = (x, d = 3) => (x == null || !Number.isFinite(x) ? '—' : x.toFixed(d));
  const esc = s => String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const tick = () => new Promise(r => setTimeout(r, 0));
  const debounce = (fn, ms = 160) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };
  const slug = s => String(s).replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '') || 'plot';
  const hexA = (hex, a) => { const n = parseInt(hex.slice(1), 16); return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`; };
  const csvCell = v => { const s = v == null ? '' : String(v); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; };
  function download(name, text) {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([text], { type: 'text/csv' }));
    a.download = name;
    document.body.appendChild(a); a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 0);
  }
  let toastTimer;
  function toast(msg, isError) {
    const el = $('#toast');
    el.textContent = msg;
    el.classList.toggle('is-error', !!isError);
    el.hidden = false;
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => { el.hidden = true; }, isError ? 9000 : 3500);
  }

  // ---------- state ----------
  const state = {
    D: null,
    active: [],
    filters: { search: '', minReads: 0, cat: {}, num: {} },
    settings: { norm: 'max', origin: 100, win: [-10, 10], sem: true, shade: true },
    tab: 'overview',
    filtered: [],
    dist: [],
    winFrac: {},
    overview: { facet: '', lines: 'sample' },
    heat: { a: '', b: '', sort: 'win', rows: 800 },
    table: { a: '', b: '', sortKey: 'win_a', sortDir: -1, page: 0, order: [] },
    detail: null,
    lastAgg: null,
  };

  const sampleColor = s => PALETTE[Math.max(0, state.D.samples.indexOf(s)) % PALETTE.length];
  function isValid(s, i) {
    const D = state.D;
    if (!s || !D.present[s] || D.present[s][i] !== 1) return false;
    const t = D.total[s][i];
    return t > 0 && t >= state.filters.minReads;
  }
  function denomFor(s, i, norm) {
    const D = state.D;
    return norm === 'max' ? D.maxv[s][i] : norm === 'sum' ? D.total[s][i] : 1;
  }
  const yTitle = norm => (norm === 'max' ? 'Normalised coverage (÷ max)' : norm === 'sum' ? 'Fraction of profile reads' : 'Read count');
  const winSorted = () => { const [a, b] = state.settings.win; return a <= b ? [a, b] : [b, a]; };
  const peakDist = (s, i) => state.dist[state.D.peak[s][i]];

  // ---------- derived quantities ----------
  function recomputeDerived() {
    const D = state.D;
    const o = state.settings.origin;
    state.dist = D.positions.map(p => o - p);
    state.colOrder = D.positions.map((_, k) => k).sort((a, b) => state.dist[a] - state.dist[b]);
    const [w0, w1] = winSorted();
    const inWin = [];
    state.dist.forEach((d, k) => { if (d >= w0 && d <= w1) inWin.push(k); });
    state.winFrac = {};
    for (const s of D.samples) {
      const P = D.prof[s], T = D.total[s];
      const out = new Float32Array(D.N);
      for (let i = 0; i < D.N; i++) {
        const t = T[i];
        if (!(t > 0)) { out[i] = NaN; continue; }
        let acc = 0;
        const off = i * D.L;
        for (const k of inWin) acc += P[off + k];
        out[i] = acc / t;
      }
      state.winFrac[s] = out;
    }
  }

  // ---------- filtering ----------
  function applyFilters() {
    const D = state.D;
    if (!D) return;
    const f = state.filters;
    const terms = f.search.trim().toLowerCase().split(/\s+/).filter(Boolean);
    const cats = Object.entries(f.cat);
    const nums = Object.entries(f.num);
    const out = [];
    outer: for (let i = 0; i < D.N; i++) {
      if (terms.length) { const h = D.hay[i]; for (const t of terms) if (!h.includes(t)) continue outer; }
      for (const [col, set] of cats) if (!set.has(D.mstr[col][i])) continue outer;
      for (const [col, [mn, mx]] of nums) {
        const v = Number(D.mstr[col][i]);
        if (!(v >= mn && v <= mx)) continue outer;
      }
      out.push(i);
    }
    state.filtered = out;
    state.table.page = 0;
    $('#filter-count').innerHTML = `<b>${fmtInt(out.length)}</b> of ${fmtInt(D.N)} sequences`;
    renderActive();
  }
  const applyFiltersSoon = debounce(applyFilters);

  // ---------- plot helpers ----------
  function baseLayout(extra = {}) {
    const axis = {
      zeroline: false, gridcolor: '#eef0f6', linecolor: '#dcdfea', showline: true, ticks: 'outside',
      tickcolor: '#dcdfea', ticklen: 4, tickfont: { size: 11, color: '#6a6f8c' }, automargin: true,
    };
    const L = {
      margin: { l: 58, r: 14, t: 12, b: 48 },
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: '#ffffff',
      font: { family: 'Manrope, system-ui, sans-serif', size: 12, color: '#22243b' },
      hoverlabel: { font: { family: 'Manrope, system-ui, sans-serif', size: 12, color: '#22243b' }, bordercolor: '#d9d3f5', bgcolor: '#ffffff' },
      hovermode: 'x unified', showlegend: false,
      xaxis: Object.assign({}, axis, { title: { text: 'Distance from STR start (nt)', font: { size: 12.5 } } }),
      yaxis: Object.assign({}, axis),
    };
    for (const k in extra) {
      if ((k === 'xaxis' || k === 'yaxis') && extra[k] && typeof extra[k] === 'object') L[k] = Object.assign({}, L[k], extra[k]);
      else L[k] = extra[k];
    }
    return L;
  }
  const plotConfig = name => ({
    responsive: true, displaylogo: false,
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'toggleSpikelines'],
    toImageButtonOptions: { format: 'svg', filename: slug(name) },
  });
  function windowShapes() {
    const shapes = [{ type: 'line', xref: 'x', yref: 'paper', x0: 0, x1: 0, y0: 0, y1: 1, line: { color: '#b4b7cc', width: 1, dash: 'dot' }, layer: 'below' }];
    if (state.settings.shade) {
      const [w0, w1] = winSorted();
      shapes.unshift({ type: 'rect', xref: 'x', yref: 'paper', x0: w0 - 0.5, x1: w1 + 0.5, y0: 0, y1: 1, fillcolor: 'rgba(232,131,107,0.12)', line: { width: 0 }, layer: 'below' });
    }
    return shapes;
  }
  function purgeIn(el) { $$('.js-plotly-plot', el).forEach(p => Plotly.purge(p)); if (el.classList.contains('js-plotly-plot')) Plotly.purge(el); }
  const emptyHTML = msg => `<div class="empty">${msg}</div>`;

  // ---------- grouping & aggregation ----------
  function groupsOf(idxs, col) {
    if (!col) return [{ name: 'All filtered sequences', idx: idxs }];
    const D = state.D;
    const vals = D.mstr[col];
    const m = new Map();
    for (const i of idxs) { const g = vals[i]; let a = m.get(g); if (!a) { a = []; m.set(g, a); } a.push(i); }
    const levels = D.metaColMap[col].levels || [...m.keys()].sort(collator.compare);
    return levels.filter(l => m.has(l)).map(l => ({ name: l, idx: m.get(l) }));
  }

  function aggregate(idxs, s) {
    const D = state.D, L = D.L, P = D.prof[s], norm = state.settings.norm;
    const sum = new Float64Array(L), sq = new Float64Array(L);
    let n = 0;
    for (const i of idxs) {
      if (!isValid(s, i)) continue;
      const d = denomFor(s, i, norm);
      if (!(d > 0)) continue;
      const off = i * L;
      for (let k = 0; k < L; k++) { const v = P[off + k] / d; sum[k] += v; sq[k] += v * v; }
      n++;
    }
    const mean = new Array(L), sem = new Array(L);
    for (let k = 0; k < L; k++) {
      if (!n) { mean[k] = null; sem[k] = null; continue; }
      const m = sum[k] / n;
      mean[k] = m;
      sem[k] = n > 1 ? Math.sqrt(Math.max(0, (sq[k] - n * m * m) / (n - 1)) / n) : 0;
    }
    return { mean, sem, n };
  }

  // =========================================================
  // OVERVIEW
  // =========================================================
  function renderOverview() {
    const D = state.D;
    const grid = $('#overview-grid');
    purgeIn(grid); purgeIn($('#overview-win'));
    const samples = state.active;
    const facet = state.overview.facet;
    const mode = facet ? state.overview.lines : 'sample';
    $('#o-lines').disabled = !facet;
    $('#overview-note').textContent = '';

    const bail = msg => {
      grid.className = 'facet-grid single'; grid.innerHTML = emptyHTML(msg);
      $('#overview-legend').innerHTML = ''; $('#overview-win').innerHTML = ''; $('#win-caption').textContent = '';
    };
    if (!samples.length) return bail('Select at least one sample in the sidebar.');
    if (!state.filtered.length) return bail('No sequences match the current filters.');

    let groups = groupsOf(state.filtered, facet);
    let note = '';
    if (groups.length > MAX_FACETS) {
      const totalGroups = groups.length;
      const lv = D.metaColMap[facet].levels;
      groups = groups.slice().sort((a, b) => b.idx.length - a.idx.length).slice(0, MAX_FACETS)
        .sort((a, b) => lv.indexOf(a.name) - lv.indexOf(b.name));
      note = `Showing the ${MAX_FACETS} largest of ${totalGroups} groups. Use the filters to pick others.`;
    }

    const agg = groups.map(g => samples.map(s => aggregate(g.idx, s)));
    state.lastAgg = { groups, samples, agg, facet };

    let panels, legend;
    if (mode === 'sample') {
      panels = groups.map((g, gi) => ({
        title: facet ? g.name : 'All filtered sequences',
        sub: `${fmtInt(g.idx.length)} sequences`,
        lines: samples.map((s, si) => ({ name: s, color: sampleColor(s), a: agg[gi][si] })),
      }));
      legend = samples.map(s => ({ name: s, color: sampleColor(s) }));
    } else {
      panels = samples.map((s, si) => ({
        title: s, sub: `${groups.length} groups`,
        lines: groups.map((g, gi) => ({ name: g.name, color: PALETTE[gi % PALETTE.length], a: agg[gi][si] })),
      }));
      legend = groups.map((g, gi) => ({ name: `${g.name} (${fmtInt(g.idx.length)})`, color: PALETTE[gi % PALETTE.length] }));
    }
    $('#overview-legend').innerHTML = legend.map(l =>
      `<span class="legend-item"><span class="legend-line" style="background:${l.color}"></span>${esc(l.name)}</span>`).join('');

    const norm = state.settings.norm;
    let ymax = 1.02;
    if (norm !== 'max') {
      ymax = 0;
      for (const p of panels) for (const ln of p.lines) for (let k = 0; k < D.L; k++) {
        const v = ln.a.mean[k];
        if (v != null) ymax = Math.max(ymax, v + (state.settings.sem ? ln.a.sem[k] : 0));
      }
      ymax = (ymax || 1) * 1.05;
    }

    grid.className = 'facet-grid' + (panels.length === 1 ? ' single' : '');
    grid.innerHTML = panels.map((p, k) => `
      <article class="card">
        <div class="card-head"><h3 class="card-title">${esc(p.title)}</h3><span class="card-meta">${p.sub}</span></div>
        <div class="card-body"><div class="plot" id="ov-plot-${k}"></div></div>
      </article>`).join('');

    const x = state.dist;
    const xRev = x.slice().reverse();
    const yt = yTitle(norm);
    panels.forEach((p, k) => {
      const traces = [];
      for (const ln of p.lines) {
        if (!ln.a.n) continue;
        if (state.settings.sem) {
          const up = ln.a.mean.map((m, j) => m + ln.a.sem[j]);
          const lo = ln.a.mean.map((m, j) => m - ln.a.sem[j]).reverse();
          traces.push({ x: x.concat(xRev), y: up.concat(lo), type: 'scatter', mode: 'lines', fill: 'toself',
            fillcolor: hexA(ln.color, 0.16), line: { width: 0 }, hoverinfo: 'skip', showlegend: false });
        }
        traces.push({ x, y: ln.a.mean, type: 'scatter', mode: 'lines', name: ln.name,
          line: { color: ln.color, width: 2 },
          hovertemplate: `${esc(ln.name)} (n=${fmtInt(ln.a.n)}): %{y:.3f}<extra></extra>` });
      }
      const layout = baseLayout({
        height: panels.length === 1 ? 420 : 300,
        yaxis: { title: { text: 'Mean ' + yt.charAt(0).toLowerCase() + yt.slice(1), font: { size: 12.5 } }, range: [0, ymax] },
        shapes: windowShapes(),
      });
      Plotly.newPlot(`ov-plot-${k}`, traces, layout, plotConfig(`profile_${p.title}`));
    });
    $('#overview-note').textContent = note;
    renderWinChart(groups, samples, facet);
  }

  function renderWinChart(groups, samples, facet) {
    const el = $('#overview-win');
    purgeIn(el);
    const [w0, w1] = winSorted();
    $('#win-caption').textContent = `For each sequence and sample: reads ending between ${w0} and ${w1} nt from the STR start, divided by all reads in that profile. Bars show the mean across sequences ± SEM. Uses raw counts, so it does not depend on the normalisation setting.`;
    const traces = samples.map(s => {
      const W = state.winFrac[s];
      const xs = [], ys = [], es = [], ns = [];
      for (const g of groups) {
        let sum = 0, sq = 0, c = 0;
        for (const i of g.idx) {
          if (!isValid(s, i)) continue;
          const v = W[i];
          if (!Number.isFinite(v)) continue;
          sum += v; sq += v * v; c++;
        }
        const m = c ? sum / c : null;
        const sd = c > 1 ? Math.sqrt(Math.max(0, (sq - c * m * m) / (c - 1))) : 0;
        xs.push(g.name); ys.push(m); es.push(c ? sd / Math.sqrt(c) : 0); ns.push(c);
      }
      return {
        type: 'bar', name: s, x: xs, y: ys, customdata: ns,
        marker: { color: sampleColor(s) },
        error_y: { type: 'data', array: es, visible: true, color: '#6a6f8c', thickness: 1.2, width: 3 },
        hovertemplate: `${esc(s)}<br>%{x}: %{y:.3f} (n=%{customdata})<extra></extra>`,
      };
    });
    const layout = baseLayout({
      height: 330, barmode: 'group', bargap: 0.25, bargroupgap: 0.08, hovermode: 'closest', showlegend: true,
      legend: { orientation: 'h', x: 0, y: 1.14, font: { family: 'Manrope, system-ui, sans-serif', size: 12 } },
      margin: { l: 58, r: 14, t: 34, b: 48 },
      xaxis: { title: { text: facet || '', font: { size: 12.5 } }, type: 'category', showgrid: false },
      yaxis: { title: { text: 'Fraction of reads in window', font: { size: 12.5 } }, rangemode: 'tozero' },
    });
    Plotly.newPlot(el, traces, layout, plotConfig('stall_window_fraction'));
  }

  function downloadOverviewCSV() {
    const A = state.lastAgg;
    if (!A) return;
    const rows = [['group', 'sample', 'dist', 'mean_norm', 'sem', 'n'].join(',')];
    A.groups.forEach((g, gi) => A.samples.forEach((s, si) => {
      const a = A.agg[gi][si];
      if (!a.n) return;
      state.dist.forEach((d, k) => rows.push([csvCell(g.name), csvCell(s), d, a.mean[k], a.sem[k], a.n].join(',')));
    }));
    download(`coverage_profile_${slug(A.facet || 'all')}_${state.settings.norm}.csv`, rows.join('\n'));
  }

  // =========================================================
  // HEATMAP
  // =========================================================
  function renderHeatmap() {
    const D = state.D;
    const el = $('#heat-plot');
    purgeIn(el);
    el.innerHTML = '';
    const a = state.heat.a;
    const b = state.heat.b && state.heat.b !== a ? state.heat.b : '';
    let norm = state.settings.norm;
    let normNote = '';
    if (b && norm === 'raw') { norm = 'sum'; normNote = ' Raw counts are not comparable between samples, so the difference uses fraction-of-reads normalisation.'; }

    const rows = state.filtered.filter(i => isValid(a, i) && (!b || isValid(b, i)));
    if (!rows.length) {
      el.innerHTML = emptyHTML('No sequences with reads in the chosen sample(s) match the current filters.');
      $('#heat-caption').textContent = '';
      return;
    }

    const W = state.winFrac;
    const byNum = (f, dir) => (x, y) => {
      const u = f(x), v = f(y);
      const nu = !Number.isFinite(u), nv = !Number.isFinite(v);
      if (nu || nv) return nu - nv;
      return dir * (u - v) || x - y;
    };
    const sk = state.heat.sort;
    let groups;
    const within = byNum(i => W[a][i], -1);
    if (sk.startsWith('meta:')) {
      groups = groupsOf(rows, sk.slice(5));
      groups.forEach(g => g.idx.sort(within));
    } else {
      let cmp;
      if (sk === 'peak') cmp = byNum(i => peakDist(a, i), -1);
      else if (sk === 'total') cmp = byNum(i => D.total[a][i], -1);
      else if (sk === 'delta' && b) cmp = byNum(i => W[b][i] - W[a][i], -1);
      else cmp = within;
      groups = [{ name: '', idx: rows.slice().sort(cmp) }];
    }

    // bin rows (within each group) when there are too many to draw one per sequence
    const maxRows = state.heat.rows;
    const total = rows.length;
    const binned = total > maxRows;
    const outRows = [];
    for (const g of groups) {
      const len = g.idx.length;
      let nb = binned ? Math.max(1, Math.round(maxRows * len / total)) : len;
      nb = Math.min(nb, len, Math.max(1, maxRows - outRows.length));
      for (let k = 0; k < nb; k++) {
        const s0 = Math.floor(k * len / nb), s1 = Math.floor((k + 1) * len / nb);
        outRows.push({ members: g.idx.slice(s0, s1), group: g.name });
      }
    }

    const cols = state.colOrder;
    const L = D.L;
    const val = (s, i, k) => {
      let v = D.prof[s][i * L + k] / denomFor(s, i, norm);
      if (norm === 'raw') v = Math.log10(1 + v);
      return v;
    };
    const z = [], text = [];
    let seen = 0;
    for (const r of outRows) {
      const acc = new Float64Array(L);
      for (const i of r.members) {
        for (let k = 0; k < L; k++) acc[k] += b ? val(b, i, k) - val(a, i, k) : val(a, i, k);
      }
      const m = r.members.length;
      z.push(cols.map(k => acc[k] / m));
      const label = m === 1
        ? D.ids[r.members[0]] + (r.group ? `, ${r.group}` : '')
        : `${fmtInt(m)} sequences (rows ${fmtInt(seen + 1)}–${fmtInt(seen + m)})${r.group ? ', ' + r.group : ''}`;
      seen += m;
      text.push(new Array(cols.length).fill(label));
    }

    // colour range: clip at the 98th percentile so a few extreme rows do not wash out the rest
    const flat = [];
    for (const row of z) for (const v of row) if (Number.isFinite(v)) flat.push(b ? Math.abs(v) : v);
    flat.sort((p, q) => p - q);
    const q98 = flat.length ? flat[Math.min(flat.length - 1, Math.floor(flat.length * 0.98))] : 1;
    let zmin, zmax, colorscale, cbTitle;
    if (b) {
      const lim = q98 || 1;
      zmin = -lim; zmax = lim;
      colorscale = [[0, '#2f8a73'], [0.25, '#a6dccb'], [0.5, '#ffffff'], [0.75, '#f6bcab'], [1, '#c9573e']];
      cbTitle = `B − A`;
    } else {
      zmin = 0; zmax = norm === 'max' ? 1 : (q98 || 1);
      colorscale = [[0, '#ffffff'], [0.12, '#eeeafb'], [0.35, '#b9aeec'], [0.65, '#7b6cd9'], [1, '#2e2675']];
      cbTitle = norm === 'max' ? '÷ max' : norm === 'sum' ? 'fraction' : 'log10(1+reads)';
    }

    const x = cols.map(k => state.dist[k]);
    const shapes = [];
    const tickvals = [], ticktext = [];
    if (groups.length > 1 || (groups[0] && groups[0].name)) {
      let start = 0;
      groups.forEach((g, gi) => {
        const nRows = outRows.filter(r => r.group === g.name).length;
        tickvals.push(start + (nRows - 1) / 2);
        ticktext.push(`${g.name} (${fmtInt(g.idx.length)})`);
        if (gi > 0) shapes.push({ type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: start - 0.5, y1: start - 0.5, line: { color: '#22243b', width: 1.2 } });
        start += nRows;
      });
    }
    if (state.settings.shade) {
      const [w0, w1] = winSorted();
      [w0 - 0.5, w1 + 0.5].forEach(xv => shapes.push({ type: 'line', xref: 'x', yref: 'paper', x0: xv, x1: xv, y0: 0, y1: 1, line: { color: '#e8836b', width: 1.5, dash: 'dash' } }));
    }

    const trace = {
      type: 'heatmap', x, y: outRows.map((_, k) => k), z, text,
      zmin, zmax, colorscale, zsmooth: false,
      hovertemplate: '%{text}<br>dist %{x} nt, value %{z:.3f}<extra></extra>',
      colorbar: { title: { text: cbTitle, side: 'right', font: { size: 11 } }, thickness: 12, len: 0.6, outlinewidth: 0, tickfont: { size: 10 } },
    };
    const layout = baseLayout({
      height: 660, hovermode: 'closest',
      margin: { l: 20, r: 10, t: 10, b: 50 },
      yaxis: {
        autorange: 'reversed', showgrid: false, showline: false, ticks: '', zeroline: false,
        tickmode: tickvals.length ? 'array' : 'auto', tickvals, ticktext,
        showticklabels: tickvals.length > 0, tickfont: { size: 11 },
      },
      xaxis: { showgrid: false },
      shapes,
    });
    Plotly.newPlot(el, [trace], layout, plotConfig(`heatmap_${a}${b ? '_vs_' + b : ''}`)).then(() => {
      el.on('plotly_click', ev => {
        const pt = ev && ev.points && ev.points[0];
        if (!pt) return;
        const r = outRows[Math.round(pt.y)];
        if (r) openDetail(r.members[Math.floor(r.members.length / 2)], true);
      });
    });

    const opt = $('#h-sort').selectedOptions[0];
    const sortLabel = opt ? opt.textContent.charAt(0).toLowerCase() + opt.textContent.slice(1) : '';
    $('#heat-caption').innerHTML =
      `<b>${fmtInt(total)}</b> sequences with reads in ${b ? `both <b>${esc(a)}</b> (A) and <b>${esc(b)}</b> (B), coloured by B − A` : `<b>${esc(a)}</b>`}` +
      (binned ? `, averaged into ${fmtInt(outRows.length)} rows (about ${(total / outRows.length).toFixed(1)} sequences per row)` : '') +
      `. Ordered by ${esc(sortLabel)}.${normNote}`;
  }

  // =========================================================
  // SEQUENCES TABLE
  // =========================================================
  function tableMetaCols() {
    const D = state.D;
    const usable = D.metaCols.filter(c => c.type !== 'seq').map(c => c.name);
    const pref = TABLE_META_PREF.filter(n => usable.includes(n));
    const rest = usable.filter(n => !pref.includes(n));
    return pref.concat(rest).slice(0, 3);
  }

  function tableColumns() {
    const D = state.D;
    const a = state.table.a;
    const b = state.table.b && state.table.b !== a ? state.table.b : '';
    const W = state.winFrac;
    const intFmt = v => (Number.isFinite(v) ? fmtInt(v) : '—');
    const cols = [{ key: 'id', label: 'Sequence ID', cls: 'id', get: i => D.ids[i] }];
    for (const c of tableMetaCols()) {
      cols.push({ key: 'meta:' + c, label: c, cls: 'meta', numericSort: D.metaColMap[c].numeric, get: i => D.mstr[c][i] });
    }
    cols.push({ key: 'tot_a', label: 'Reads A', num: true, get: i => (isValid(a, i) ? D.total[a][i] : NaN), fmt: intFmt });
    cols.push({ key: 'peak_a', label: 'Peak A (nt)', num: true, get: i => (isValid(a, i) ? peakDist(a, i) : NaN), fmt: v => (Number.isFinite(v) ? String(v) : '—') });
    cols.push({ key: 'win_a', label: 'Window A', num: true, get: i => (isValid(a, i) ? W[a][i] : NaN), fmt: v => fmt(v) });
    if (b) {
      cols.push({ key: 'tot_b', label: 'Reads B', num: true, get: i => (isValid(b, i) ? D.total[b][i] : NaN), fmt: intFmt });
      cols.push({ key: 'win_b', label: 'Window B', num: true, get: i => (isValid(b, i) ? W[b][i] : NaN), fmt: v => fmt(v) });
      cols.push({ key: 'delta', label: 'Δ B − A', num: true,
        get: i => (isValid(a, i) && isValid(b, i) ? W[b][i] - W[a][i] : NaN),
        fmt: v => (Number.isFinite(v) ? `<span class="${v > 0 ? 'delta-pos' : v < 0 ? 'delta-neg' : ''}">${v > 0 ? '+' : ''}${v.toFixed(3)}</span>` : '—') });
    }
    return cols;
  }

  function sortedOrder(cols) {
    let col = cols.find(c => c.key === state.table.sortKey);
    if (!col) { col = cols.find(c => c.key === 'win_a'); state.table.sortKey = 'win_a'; state.table.sortDir = -1; }
    const dir = state.table.sortDir;
    const rows = state.filtered.slice();
    if (col.num || col.numericSort) {
      const key = new Float64Array(state.D.N);
      for (const i of rows) { const v = Number(col.get(i)); key[i] = Number.isFinite(v) ? v : NaN; }
      rows.sort((x, y) => {
        const u = key[x], v = key[y];
        const nu = Number.isNaN(u), nv = Number.isNaN(v);
        if (nu || nv) return nu - nv;
        return dir * (u - v) || x - y;
      });
    } else {
      rows.sort((x, y) => dir * collator.compare(col.get(x), col.get(y)) || x - y);
    }
    return { rows, col };
  }

  function renderTable() {
    const cols = tableColumns();
    const { rows, col } = sortedOrder(cols);
    state.table.order = rows;
    const pages = Math.max(1, Math.ceil(rows.length / PAGE));
    state.table.page = Math.max(0, Math.min(state.table.page, pages - 1));
    const start = state.table.page * PAGE;
    const slice = rows.slice(start, start + PAGE);
    const dir = state.table.sortDir;

    $('#t-table thead').innerHTML = '<tr>' + cols.map(c =>
      `<th class="${c.num ? 'num' : ''}${c.key === col.key ? ' is-sorted' : ''}" data-key="${esc(c.key)}" title="Sort by ${esc(c.label)}">${esc(c.label)}${c.key === col.key ? (dir > 0 ? ' ↑' : ' ↓') : ''}</th>`).join('') + '</tr>';
    $('#t-table tbody').innerHTML = slice.length
      ? slice.map(i => `<tr data-i="${i}"${i === state.detail ? ' class="is-selected"' : ''}>` + cols.map(c => {
        const v = c.get(i);
        if (c.num) return `<td class="num">${c.fmt(v)}</td>`;
        const shown = v === 'NA' ? '<span class="na">NA</span>' : esc(v);
        return `<td class="${c.cls || ''}" title="${esc(v)}">${shown}</td>`;
      }).join('') + '</tr>').join('')
      : `<tr><td colspan="${cols.length}" class="na" style="padding:2rem;text-align:center">No sequences match the current filters.</td></tr>`;

    const b = state.table.b && state.table.b !== state.table.a ? state.table.b : '';
    $('#t-pager').innerHTML = `
      <span>${rows.length ? `${fmtInt(start + 1)}–${fmtInt(Math.min(start + PAGE, rows.length))} of ${fmtInt(rows.length)}` : '0 sequences'}
        <span class="pager-ab">A is <b>${esc(state.table.a)}</b>${b ? `, B is <b>${esc(b)}</b>` : ''}</span></span>
      <span class="pager-btns">
        <button class="btn-icon" data-page="first" ${state.table.page === 0 ? 'disabled' : ''} aria-label="First page">«</button>
        <button class="btn-icon" data-page="prev" ${state.table.page === 0 ? 'disabled' : ''} aria-label="Previous page">‹</button>
        <button class="btn-icon" data-page="next" ${state.table.page >= pages - 1 ? 'disabled' : ''} aria-label="Next page">›</button>
        <button class="btn-icon" data-page="last" ${state.table.page >= pages - 1 ? 'disabled' : ''} aria-label="Last page">»</button>
      </span>`;
  }

  function downloadTableCSV() {
    const D = state.D;
    const { rows } = sortedOrder(tableColumns());
    const metaNames = D.metaCols.filter(c => c.type !== 'seq').map(c => c.name);
    const header = ['seq.id', ...metaNames];
    for (const s of D.samples) header.push(`reads.${s}`, `peak.dist.${s}`, `window.frac.${s}`);
    const out = [header.map(csvCell).join(',')];
    for (const i of rows) {
      const r = [D.ids[i], ...metaNames.map(n => D.mstr[n][i])];
      for (const s of D.samples) {
        if (isValid(s, i)) r.push(D.total[s][i], peakDist(s, i), state.winFrac[s][i].toFixed(5));
        else r.push('', '', '');
      }
      out.push(r.map(csvCell).join(','));
    }
    const [w0, w1] = winSorted();
    download(`sequences_window_${w0}_${w1}.csv`, out.join('\n'));
  }

  // =========================================================
  // SEQUENCE DETAIL
  // =========================================================
  function openDetail(i, switchTab) {
    state.detail = i;
    try { history.replaceState(null, '', '#seq=' + encodeURIComponent(state.D.ids[i])); } catch (e) { /* not allowed on file:// in some browsers */ }
    if (switchTab && state.tab !== 'sequences') setTab('sequences');
    else { renderDetail(); highlightRow(); }
    $('#detail').scrollIntoView({ behavior: 'smooth', block: 'start' });
  }
  function closeDetail() {
    purgeIn($('#detail'));
    state.detail = null;
    $('#detail').hidden = true;
    try { history.replaceState(null, '', location.pathname + location.search); } catch (e) { /* ignore */ }
    highlightRow();
  }
  function highlightRow() {
    $$('#t-table tbody tr').forEach(tr => tr.classList.toggle('is-selected', Number(tr.dataset.i) === state.detail));
  }
  function stepDetail(delta) {
    const order = state.table.order;
    const p = order.indexOf(state.detail);
    if (p < 0) return;
    const next = order[p + delta];
    if (next === undefined) return;
    const page = Math.floor((p + delta) / PAGE);
    if (page !== state.table.page) { state.table.page = page; renderTable(); }
    openDetail(next, false);
  }

  function seqHTML(i) {
    const D = state.D;
    if (!D.seqCol) return '';
    const seq = String(D.mstr[D.seqCol][i] || '').toUpperCase();
    if (!seq || seq === 'NA') return '';
    const aligned = seq.length === D.L;
    const cls = new Array(seq.length).fill('');

    // highlight tandem copies of any motif given in the metadata
    const motifs = new Set();
    for (const c of D.metaCols) {
      if (c.type === 'seq' || !/motif/i.test(c.name)) continue;
      const v = String(D.mstr[c.name][i] || '').toUpperCase();
      if (/^[ACGTU]{1,12}$/.test(v)) motifs.add(v);
    }
    for (const m of motifs) {
      const re = new RegExp(`(?:${m}){${Math.max(2, Math.ceil(6 / m.length))},}`, 'g');
      let mt;
      while ((mt = re.exec(seq))) {
        for (let k = mt.index; k < mt.index + mt[0].length; k++) if (!cls[k].includes('str')) cls[k] += ' str';
        if (!mt[0].length) re.lastIndex++;
      }
    }
    let peakLabel = '';
    if (aligned) {
      const [w0, w1] = winSorted();
      state.dist.forEach((d, k) => { if (d >= w0 && d <= w1) cls[k] += ' win'; });
      const ps = state.active.find(s => isValid(s, i));
      if (ps) { cls[D.peak[ps][i]] += ' peak'; peakLabel = ps; }
    }

    const LINE = 50;
    let html = '';
    for (let s0 = 0; s0 < seq.length; s0 += LINE) {
      let bases = '';
      for (let k = s0; k < Math.min(seq.length, s0 + LINE); k++) {
        if (k > s0 && (k - s0) % 10 === 0) bases += ' ';
        const tip = aligned ? `pos ${D.positions[k]}, dist ${state.dist[k]}` : `nt ${k + 1}`;
        bases += `<span class="b nt-${seq[k]}${cls[k]}" title="${tip}">${seq[k]}</span>`;
      }
      html += `<div class="seq-line"><span class="seq-pos">${aligned ? D.positions[s0] : s0 + 1}</span><span>${bases}</span></div>`;
    }
    const key = [
      motifs.size ? '<span><i class="key-str"></i>STR repeat</span>' : '',
      aligned ? '<span><i class="key-win"></i>Stall window</span>' : '',
      peakLabel ? `<span><i class="key-peak"></i>Highest position in ${esc(peakLabel)}</span>` : '',
      !aligned ? `<span>The sequence is ${seq.length} nt but the profiles have ${D.L} positions, so bases are not aligned to the plot.</span>` : '',
    ].join('');
    return `<div class="seqbox">${html}</div><div class="seq-key">${key}</div>`;
  }

  function renderDetail() {
    const card = $('#detail');
    purgeIn(card);
    const i = state.detail;
    const D = state.D;
    if (i == null || !D) { card.hidden = true; return; }
    card.hidden = false;
    const order = state.table.order;
    const p = order.indexOf(i);
    const chips = D.metaCols.filter(c => c.type !== 'seq' && D.mstr[c.name][i] !== 'NA')
      .map(c => `<span class="chip"><b>${esc(c.name)}</b><span>${esc(D.mstr[c.name][i])}</span></span>`).join('');
    const metrics = state.active.map(s => {
      const ok = isValid(s, i);
      const has = D.present[s][i] && D.total[s][i] > 0;
      const body = ok
        ? `<b>${fmtInt(D.total[s][i])}</b> reads<br>peak at <b>${peakDist(s, i)}</b> nt<br>window fraction <b>${fmt(state.winFrac[s][i])}</b>`
        : has ? `${fmtInt(D.total[s][i])} reads — below the minimum` : 'No reads';
      return `<div class="metric"><div class="metric-name"><span class="swatch" style="background:${sampleColor(s)}"></span>${esc(s)}</div><div class="metric-vals">${body}</div></div>`;
    }).join('');

    card.innerHTML = `
      <div class="detail-head">
        <div>
          <div class="eyebrow">${p >= 0 ? `Sequence ${fmtInt(p + 1)} of ${fmtInt(order.length)}` : 'Sequence outside the current filters'}</div>
          <div class="detail-title">${esc(D.ids[i])}</div>
        </div>
        <div class="detail-actions">
          <button class="btn-icon" data-act="prev" ${p <= 0 ? 'disabled' : ''} aria-label="Previous sequence" title="Previous (←)">‹</button>
          <button class="btn-icon" data-act="next" ${p < 0 || p >= order.length - 1 ? 'disabled' : ''} aria-label="Next sequence" title="Next (→)">›</button>
          <button class="btn-icon" data-act="close" aria-label="Close" title="Close (Esc)">✕</button>
        </div>
      </div>
      <div class="chips">${chips}</div>
      <div class="card-body"><div class="plot" id="detail-plot"></div></div>
      <div class="metrics">${metrics}</div>
      ${seqHTML(i)}`;

    const norm = state.settings.norm;
    const traces = [];
    let nSample = 0;
    for (const s of state.active) {
      if (!(D.present[s][i] && D.total[s][i] > 0)) continue;
      nSample++;
      const d = denomFor(s, i, norm);
      const y = new Array(D.L);
      for (let k = 0; k < D.L; k++) y[k] = D.prof[s][i * D.L + k] / d;
      traces.push({ x: state.dist, y, type: 'scatter', mode: 'lines', name: s, line: { color: sampleColor(s), width: 1.8 },
        hovertemplate: `${esc(s)}: %{y:.3f}<extra></extra>` });
    }
    if (D.seqCol) {
      const seq = String(D.mstr[D.seqCol][i] || '').toUpperCase();
      if (seq.length === D.L) {
        traces.push({ x: state.dist, y: new Array(D.L).fill(0), type: 'scatter', mode: 'markers', marker: { opacity: 0, size: 1 },
          customdata: D.positions.map((pp, k) => [seq[k], pp]), hovertemplate: 'base %{customdata[0]}, pos %{customdata[1]}<extra></extra>', showlegend: false });
      }
    }
    if (!nSample) { $('#detail-plot').innerHTML = emptyHTML('No reads for this sequence in the selected samples.'); return; }
    Plotly.newPlot('detail-plot', traces, baseLayout({
      height: 320, showlegend: true,
      legend: { orientation: 'h', x: 0, y: 1.12, font: { family: 'Manrope, system-ui, sans-serif', size: 12 } },
      margin: { l: 58, r: 14, t: 30, b: 48 },
      yaxis: { title: { text: yTitle(norm), font: { size: 12.5 } }, rangemode: 'tozero' },
      shapes: windowShapes(),
    }), plotConfig(`profile_${D.ids[i]}`));
  }

  // =========================================================
  // TABS & RENDER DISPATCH
  // =========================================================
  function setTab(t) {
    state.tab = t;
    $$('.tab').forEach(b => { const on = b.dataset.tab === t; b.classList.toggle('is-active', on); b.setAttribute('aria-selected', String(on)); });
    $$('.view').forEach(v => v.classList.toggle('is-active', v.id === 'view-' + t));
    renderActive();
  }
  function renderActive() {
    if (!state.D || typeof Plotly === 'undefined') return;
    if (state.tab === 'overview') renderOverview();
    else if (state.tab === 'heatmap') renderHeatmap();
    else if (state.tab === 'sequences') { renderTable(); renderDetail(); }
  }
  const renderSoon = debounce(renderActive, 120);

  // =========================================================
  // DATASET SETUP
  // =========================================================
  function setDataset(D) {
    state.D = D;
    if (Number.isFinite(D.origin)) { state.settings.origin = D.origin; $('#s-origin').value = D.origin; }
    state.active = D.samples.slice();
    state.filters = { search: '', minReads: 0, cat: {}, num: {} };
    $('#f-search').value = ''; $('#f-minreads').value = 0;
    const cats = D.metaCols.filter(c => c.type === 'cat' && c.levels.length >= 2);
    const keep = state.overview.facet && cats.some(c => c.name === state.overview.facet);
    state.overview.facet = keep ? state.overview.facet : (cats.some(c => c.name === 'motif.type') ? 'motif.type' : '');
    state.heat.a = D.samples[0]; state.heat.b = '';
    state.table.a = D.samples[0]; state.table.b = D.samples[1] || '';
    state.table.sortKey = state.table.b ? 'delta' : 'win_a'; state.table.sortDir = -1; state.table.page = 0;
    state.detail = null;
    $('#detail').hidden = true;

    const label = D.isDemo ? 'Simulated demo dataset' : D.name;
    $('#dataset-label').innerHTML = `<span class="ds-name">${esc(label)}</span><span class="ds-chip">${fmtInt(D.N)} sequences</span><span class="ds-chip">${D.samples.length} sample${D.samples.length === 1 ? '' : 's'}</span><span class="ds-chip">${D.L} positions</span>`;
    $('#dataset-label').title = D.files.length ? D.files.map(f => `${f.name} (${f.kind}, ${fmtInt(f.rows)} rows)`).join('\n') : label;
    $('#demo-banner').hidden = !D.isDemo;
    document.title = D.isDemo ? 'Stalling Profile Explorer' : `${D.name} | Stalling Profile Explorer`;
    if (!D.isDemo && D.metaMissing) toast('Profiles loaded without a sequences file. Add sequences.tsv.gz to get metadata filters and sequences.');

    buildSampleList();
    buildMetaFilters();
    buildViewControls();
    recomputeDerived();
    applyFilters();

    const m = location.hash.match(/^#seq=(.+)$/);
    if (m) {
      const k = D.ids.indexOf(decodeURIComponent(m[1]));
      if (k >= 0) openDetail(k, true);
    }
  }

  function buildSampleList() {
    const D = state.D;
    $('#sample-list').innerHTML = D.samples.map(s => {
      let n = 0;
      for (let i = 0; i < D.N; i++) if (D.present[s][i] && D.total[s][i] > 0) n++;
      return `<label class="sample-row" title="${esc(s)}: ${fmtInt(n)} sequences with reads">
        <input type="checkbox" value="${esc(s)}" ${state.active.includes(s) ? 'checked' : ''}>
        <span class="swatch" style="background:${sampleColor(s)}"></span>
        <span class="sample-name">${esc(s)}</span><span class="sample-n">${fmtInt(n)}</span></label>`;
    }).join('');
  }

  function buildMetaFilters() {
    const D = state.D;
    const el = $('#meta-filters');
    const parts = [];
    for (const c of D.metaCols) {
      if (c.type === 'cat' && c.levels.length >= 2) {
        const counts = new Map();
        for (const v of D.mstr[c.name]) counts.set(v, (counts.get(v) || 0) + 1);
        parts.push(`<details class="ms" data-col="${esc(c.name)}">
          <summary><span class="ms-name">${esc(c.name)}</span><span class="ms-val">All</span></summary>
          <div class="ms-body">
            <div class="ms-actions"><button type="button" class="btn-link" data-act="all">Select all</button><button type="button" class="btn-link" data-act="none">Clear</button></div>
            ${c.levels.map(l => `<label class="ms-opt"><input type="checkbox" value="${esc(l)}" checked><span>${esc(l)}</span><span class="ms-n">${fmtInt(counts.get(l) || 0)}</span></label>`).join('')}
          </div></details>`);
      } else if (c.type === 'num') {
        parts.push(`<div class="field" data-num="${esc(c.name)}"><span class="num-label">${esc(c.name)}</span>
          <div class="pair"><input type="number" step="any" placeholder="${c.min}" data-b="min" aria-label="${esc(c.name)} minimum"><span>to</span><input type="number" step="any" placeholder="${c.max}" data-b="max" aria-label="${esc(c.name)} maximum"></div></div>`);
      }
    }
    el.innerHTML = parts.join('');

    $$('.ms', el).forEach(det => {
      const col = det.dataset.col;
      const update = () => {
        const boxes = $$('input', det);
        const sel = boxes.filter(b => b.checked).map(b => b.value);
        const out = $('.ms-val', det);
        if (sel.length === boxes.length) { delete state.filters.cat[col]; out.textContent = 'All'; }
        else {
          state.filters.cat[col] = new Set(sel);
          out.textContent = sel.length === 0 ? 'None' : sel.length === 1 ? sel[0] : `${sel.length} of ${boxes.length}`;
        }
        applyFiltersSoon();
      };
      det.addEventListener('change', update);
      det.addEventListener('click', e => {
        const act = e.target && e.target.dataset ? e.target.dataset.act : null;
        if (!act) return;
        e.preventDefault();
        $$('input', det).forEach(b => { b.checked = act === 'all'; });
        update();
      });
    });
    $$('[data-num]', el).forEach(box => {
      const col = box.dataset.num;
      box.addEventListener('input', () => {
        const mn = $('[data-b=min]', box).value, mx = $('[data-b=max]', box).value;
        if (mn === '' && mx === '') delete state.filters.num[col];
        else state.filters.num[col] = [mn === '' ? -Infinity : Number(mn), mx === '' ? Infinity : Number(mx)];
        applyFiltersSoon();
      });
    });
  }

  function buildViewControls() {
    const D = state.D;
    const cats = D.metaCols.filter(c => c.type === 'cat' && c.levels.length >= 2);
    $('#o-facet').innerHTML = '<option value="">Nothing (all together)</option>' +
      cats.map(c => `<option value="${esc(c.name)}">${esc(c.name)}</option>`).join('');
    $('#o-facet').value = state.overview.facet;
    $('#o-lines').value = state.overview.lines;
    const opts = D.samples.map(s => `<option value="${esc(s)}">${esc(s)}</option>`).join('');
    const none = '<option value="">— none —</option>';
    $('#h-a').innerHTML = opts; $('#h-a').value = state.heat.a;
    $('#h-b').innerHTML = none + opts; $('#h-b').value = state.heat.b;
    $('#t-a').innerHTML = opts; $('#t-a').value = state.table.a;
    $('#t-b').innerHTML = none + opts; $('#t-b').value = state.table.b;
    buildHeatSort();
  }
  function buildHeatSort() {
    const D = state.D;
    const cats = D.metaCols.filter(c => c.type === 'cat' && c.levels.length >= 2);
    const b = state.heat.b && state.heat.b !== state.heat.a;
    const options = [['win', 'Window fraction in A'], ['peak', 'Peak position in A'], ['total', 'Total reads in A']];
    if (b) options.push(['delta', 'Δ window fraction (B − A)']);
    cats.forEach(c => options.push(['meta:' + c.name, `Group: ${c.name}`]));
    if (!options.some(o => o[0] === state.heat.sort)) state.heat.sort = 'win';
    $('#h-sort').innerHTML = options.map(([v, l]) => `<option value="${esc(v)}">${esc(l)}</option>`).join('');
    $('#h-sort').value = state.heat.sort;
  }

  // =========================================================
  // LOADING
  // =========================================================
  function showLoading(msg, frac) {
    $('#loading').hidden = false;
    $('#loading-msg').textContent = msg;
    $('#loading-bar').style.width = frac == null ? '0%' : `${Math.round(frac * 100)}%`;
  }
  const hideLoading = () => { $('#loading').hidden = true; };

  async function loadFiles(fileList) {
    const files = Array.from(fileList || []);
    if (!files.length) return;
    try {
      const B = SD.createBuilder();
      for (const f of files) {
        showLoading(`Reading ${f.name}…`, 0);
        await tick();
        const text = await SD.bytesToText(await f.arrayBuffer());
        await SD.parseText(B, text, f.name, p => showLoading(`Parsing ${f.name}…`, p));
      }
      showLoading('Building dataset…', 1);
      await tick();
      const D = SD.finalize(B, { name: files.map(f => SD.fileStem(f.name)).join(' + '), origin: state.settings.origin });
      setDataset(D);
      toast(`Loaded ${fmtInt(D.N)} sequences across ${D.samples.length} sample${D.samples.length === 1 ? '' : 's'}.`);
    } catch (err) {
      console.error(err);
      toast(err.message || String(err), true);
    } finally {
      hideLoading();
    }
  }

  async function tryLocalManifest() {
    if (!/^https?:$/.test(location.protocol)) return false;
    let manifest;
    try {
      const r = await fetch('data/local/manifest.json', { cache: 'no-store' });
      if (!r.ok) return false;
      manifest = await r.json();
    } catch (e) { return false; }
    try {
      const B = SD.createBuilder();
      for (const fn of manifest.files || []) {
        showLoading(`Opening ${fn}…`, 0);
        const r = await fetch('data/local/' + fn, { cache: 'no-store' });
        if (!r.ok) throw new Error(`Could not read data/local/${fn} (${r.status}).`);
        const text = await SD.bytesToText(await r.arrayBuffer());
        await SD.parseText(B, text, fn, p => showLoading(`Parsing ${fn}…`, p));
      }
      showLoading('Building dataset…', 1);
      await tick();
      const origin = Number(manifest.origin);
      setDataset(SD.finalize(B, { name: manifest.name || 'Local data', origin: Number.isFinite(origin) ? origin : state.settings.origin }));
      return true;
    } catch (err) {
      console.error(err);
      toast(err.message || String(err), true);
      return false;
    } finally {
      hideLoading();
    }
  }

  // =========================================================
  // EVENTS
  // =========================================================
  function wire() {
    $$('.tab').forEach(b => b.addEventListener('click', () => setTab(b.dataset.tab)));

    const fileInput = $('#file-input');
    $('#btn-load').addEventListener('click', () => fileInput.click());
    $('#dropzone').addEventListener('click', () => fileInput.click());
    $('#dropzone').addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); fileInput.click(); } });
    fileInput.addEventListener('change', () => { loadFiles(fileInput.files); fileInput.value = ''; });

    let dragDepth = 0;
    const hasFiles = e => e.dataTransfer && Array.from(e.dataTransfer.types || []).includes('Files');
    window.addEventListener('dragenter', e => { if (!hasFiles(e)) return; e.preventDefault(); dragDepth++; $('#drop-overlay').hidden = false; });
    window.addEventListener('dragover', e => { if (hasFiles(e)) e.preventDefault(); });
    window.addEventListener('dragleave', e => { if (!hasFiles(e)) return; dragDepth = Math.max(0, dragDepth - 1); if (!dragDepth) $('#drop-overlay').hidden = true; });
    window.addEventListener('drop', e => {
      if (!hasFiles(e)) return;
      e.preventDefault(); dragDepth = 0; $('#drop-overlay').hidden = true;
      loadFiles(e.dataTransfer.files);
    });

    $('#sample-list').addEventListener('change', () => {
      state.active = $$('#sample-list input').filter(x => x.checked).map(x => x.value);
      renderActive();
    });
    $('#f-search').addEventListener('input', e => { state.filters.search = e.target.value; applyFiltersSoon(); });
    $('#f-minreads').addEventListener('input', e => { state.filters.minReads = Math.max(0, Number(e.target.value) || 0); renderSoon(); });
    $('#btn-reset').addEventListener('click', () => {
      state.filters = { search: '', minReads: 0, cat: {}, num: {} };
      $('#f-search').value = ''; $('#f-minreads').value = 0;
      buildMetaFilters();
      applyFilters();
    });

    $('#s-norm').addEventListener('change', e => { state.settings.norm = e.target.value; renderActive(); });
    $('#s-origin').addEventListener('input', e => {
      const v = Number(e.target.value);
      if (e.target.value === '' || !Number.isFinite(v)) return;
      state.settings.origin = v; recomputeDerived(); renderSoon();
    });
    const winInput = () => {
      const a = $('#s-win-from').value, b = $('#s-win-to').value;
      if (a === '' || b === '' || !Number.isFinite(Number(a)) || !Number.isFinite(Number(b))) return;
      state.settings.win = [Number(a), Number(b)]; recomputeDerived(); renderSoon();
    };
    $('#s-win-from').addEventListener('input', winInput);
    $('#s-win-to').addEventListener('input', winInput);
    $('#s-sem').addEventListener('change', e => { state.settings.sem = e.target.checked; renderActive(); });
    $('#s-shade').addEventListener('change', e => { state.settings.shade = e.target.checked; renderActive(); });

    $('#o-facet').addEventListener('change', e => { state.overview.facet = e.target.value; renderActive(); });
    $('#o-lines').addEventListener('change', e => { state.overview.lines = e.target.value; renderActive(); });
    $('#o-csv').addEventListener('click', downloadOverviewCSV);

    $('#h-a').addEventListener('change', e => { state.heat.a = e.target.value; buildHeatSort(); renderActive(); });
    $('#h-b').addEventListener('change', e => { state.heat.b = e.target.value; buildHeatSort(); renderActive(); });
    $('#h-sort').addEventListener('change', e => { state.heat.sort = e.target.value; renderActive(); });
    $('#h-rows').addEventListener('change', e => { state.heat.rows = Number(e.target.value); renderActive(); });

    $('#t-a').addEventListener('change', e => { state.table.a = e.target.value; state.table.page = 0; renderActive(); });
    $('#t-b').addEventListener('change', e => {
      state.table.b = e.target.value; state.table.page = 0;
      if (!state.table.b && ['delta', 'win_b', 'tot_b'].includes(state.table.sortKey)) state.table.sortKey = 'win_a';
      renderActive();
    });
    $('#t-csv').addEventListener('click', downloadTableCSV);
    $('#t-table thead').addEventListener('click', e => {
      const th = e.target.closest('th');
      if (!th) return;
      const key = th.dataset.key;
      if (state.table.sortKey === key) state.table.sortDir *= -1;
      else { state.table.sortKey = key; state.table.sortDir = (key === 'id' || key.startsWith('meta:')) ? 1 : -1; }
      state.table.page = 0;
      renderTable();
    });
    $('#t-table tbody').addEventListener('click', e => {
      const tr = e.target.closest('tr[data-i]');
      if (tr) openDetail(Number(tr.dataset.i), false);
    });
    $('#t-pager').addEventListener('click', e => {
      const btn = e.target.closest('[data-page]');
      if (!btn || btn.disabled) return;
      const pages = Math.max(1, Math.ceil(state.table.order.length / PAGE));
      const p = state.table.page;
      state.table.page = { first: 0, prev: p - 1, next: p + 1, last: pages - 1 }[btn.dataset.page];
      renderTable();
    });
    $('#detail').addEventListener('click', e => {
      const btn = e.target.closest('[data-act]');
      if (!btn || btn.disabled) return;
      if (btn.dataset.act === 'prev') stepDetail(-1);
      else if (btn.dataset.act === 'next') stepDetail(1);
      else if (btn.dataset.act === 'close') closeDetail();
    });

    document.addEventListener('keydown', e => {
      if (state.tab !== 'sequences' || state.detail == null) return;
      if (/^(INPUT|SELECT|TEXTAREA)$/.test(e.target.tagName)) return;
      if (e.key === 'ArrowRight') { e.preventDefault(); stepDetail(1); }
      else if (e.key === 'ArrowLeft') { e.preventDefault(); stepDetail(-1); }
      else if (e.key === 'Escape') closeDetail();
    });
  }

  async function init() {
    wire();
    if (typeof Plotly === 'undefined') {
      toast('The plotting library could not be loaded. Check your internet connection and reload.', true);
    }
    const loaded = await tryLocalManifest();
    if (!loaded) {
      showLoading('Generating demo data…', null);
      await tick();
      try { setDataset(window.SPEDemo.make()); } finally { hideLoading(); }
    }
  }

  init();
})();
