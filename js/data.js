/* =========================================================
   data.js — reading and assembling coverage-profile datasets.
   No DOM code here, so it can also be tested with Node.

   Accepted inputs (tab- or comma-separated, optionally .gz):
   - metadata:  seq.id + any other columns (one row per sequence)
   - wide:      seq.id, sample, p1 … pN   (one row per sequence x sample)
   - long:      seq.id, pos, depth, sample (or one *.depth column per sample)
   ========================================================= */
(function (root) {
  'use strict';

  const collator = new Intl.Collator('en', { numeric: true, sensitivity: 'base' });
  const tick = () => new Promise(r => setTimeout(r, 0));

  // ---------- bytes -> text (transparent gzip) ----------
  async function bytesToText(buf) {
    const u8 = buf instanceof Uint8Array ? buf : new Uint8Array(buf);
    if (u8.length > 2 && u8[0] === 0x1f && u8[1] === 0x8b) {
      if (typeof DecompressionStream === 'undefined') {
        throw new Error('This browser cannot read .gz files. Use a recent Chrome, Firefox or Safari, or decompress the file first.');
      }
      const stream = new Blob([u8]).stream().pipeThrough(new DecompressionStream('gzip'));
      return await new Response(stream).text();
    }
    return new TextDecoder().decode(u8);
  }

  // ---------- line splitting ----------
  function splitLine(line, sep) {
    if (line.indexOf('"') === -1) return line.split(sep);
    const out = [];
    let cur = '', q = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (q) {
        if (c === '"') {
          if (line[i + 1] === '"') { cur += '"'; i++; } else q = false;
        } else cur += c;
      } else if (c === '"') q = true;
      else if (c === sep) { out.push(cur); cur = ''; }
      else cur += c;
    }
    out.push(cur);
    return out;
  }

  function detectSep(header) {
    const tabs = (header.match(/\t/g) || []).length;
    const commas = (header.match(/,/g) || []).length;
    return tabs > 0 && tabs >= commas ? '\t' : ',';
  }

  function fileStem(name) {
    return String(name).split(/[\\/]/).pop()
      .replace(/\.gz$/i, '').replace(/\.(tsv|csv|txt)$/i, '');
  }

  const ID_NAMES = ['seq.id', 'seq_id', 'seqid', 'sequence_id', 'id'];

  function classify(headers) {
    const lower = headers.map(h => h.toLowerCase());
    let idCol = -1;
    for (const n of ID_NAMES) { idCol = lower.indexOf(n); if (idCol >= 0) break; }
    if (idCol < 0) return { kind: 'unknown' };
    const posCol = lower.indexOf('pos');
    const sampleCol = lower.indexOf('sample');
    const posLike = [];
    headers.forEach((h, j) => { if (/^p?\d+$/i.test(h)) posLike.push(j); });
    if (posCol >= 0) return { kind: 'long', idCol, posCol, sampleCol, lower };
    if (posLike.length >= 5) return { kind: 'wide', idCol, sampleCol, posLike, lower };
    return { kind: 'meta', idCol, lower };
  }

  // ---------- builder: accumulates several files ----------
  function createBuilder() {
    const B = { ids: new Map(), idList: [], meta: [], chunks: [], sampleOrder: [], sampleSet: new Set(), files: [] };
    B.idx = id => {
      let k = B.ids.get(id);
      if (k === undefined) { k = B.idList.length; B.ids.set(id, k); B.idList.push(id); B.meta.push(null); }
      return k;
    };
    B.addSample = s => { if (!B.sampleSet.has(s)) { B.sampleSet.add(s); B.sampleOrder.push(s); } };
    B.setMeta = (k, obj, overwrite) => {
      if (!B.meta[k]) B.meta[k] = {};
      const m = B.meta[k];
      for (const key in obj) if (overwrite || !(key in m)) m[key] = obj[key];
    };
    return B;
  }

  function cleanLine(s) {
    if (!s) return '';
    if (s.charCodeAt(s.length - 1) === 13) s = s.slice(0, -1);
    return s.trim() ? s : '';
  }

  async function parseText(B, text, fileName, onProgress) {
    if (text.charCodeAt(0) === 0xFEFF) text = text.slice(1);
    const lines = text.split('\n');
    let h = 0;
    while (h < lines.length && !lines[h].trim()) h++;
    if (h >= lines.length) throw new Error(`${fileName} is empty.`);
    const header = cleanLine(lines[h]);
    const sep = detectSep(header);
    const headers = splitLine(header, sep).map(s => s.trim());
    const cls = classify(headers);
    if (cls.kind === 'unknown') {
      throw new Error(`${fileName}: no "seq.id" column. Columns found: ${headers.slice(0, 6).join(', ')}${headers.length > 6 ? '…' : ''}`);
    }
    const stem = fileStem(fileName);
    const n = lines.length;
    const report = async i => { if (onProgress) onProgress(i / n); await tick(); };
    const get = (f, j) => (f[j] === undefined ? '' : f[j].trim());
    let rows = 0;

    if (cls.kind === 'meta') {
      const metaIdx = headers.map((_, j) => j).filter(j => j !== cls.idCol);
      for (let i = h + 1; i < n; i++) {
        const line = cleanLine(lines[i]);
        if (!line) continue;
        const f = splitLine(line, sep);
        const id = get(f, cls.idCol);
        if (!id) continue;
        const o = {};
        for (const j of metaIdx) o[headers[j]] = get(f, j);
        B.setMeta(B.idx(id), o, true);
        rows++;
        if ((i & 8191) === 0) await report(i);
      }
      B.files.push({ name: fileName, kind: 'meta', rows });
      return;
    }

    if (cls.kind === 'wide') {
      const posCols = cls.posLike;
      const positions = posCols.map(j => parseInt(headers[j].replace(/^p/i, ''), 10));
      const skip = new Set([cls.idCol, cls.sampleCol, ...posCols]);
      const metaIdx = headers.map((_, j) => j).filter(j => !skip.has(j));
      const chunk = { positions, entries: [] };
      for (let i = h + 1; i < n; i++) {
        const line = cleanLine(lines[i]);
        if (!line) continue;
        const f = splitLine(line, sep);
        const id = get(f, cls.idCol);
        if (!id) continue;
        const sample = cls.sampleCol >= 0 ? (get(f, cls.sampleCol) || stem) : stem;
        const k = B.idx(id);
        B.addSample(sample);
        const vals = new Float32Array(posCols.length);
        for (let j = 0; j < posCols.length; j++) {
          const v = +f[posCols[j]];
          vals[j] = Number.isFinite(v) ? v : 0;
        }
        chunk.entries.push({ sample, idx: k, vals });
        if (metaIdx.length && !B.meta[k]) {
          const o = {};
          for (const j of metaIdx) o[headers[j]] = get(f, j);
          B.setMeta(k, o, false);
        }
        rows++;
        if ((i & 4095) === 0) await report(i);
      }
      B.chunks.push(chunk);
      B.files.push({ name: fileName, kind: 'profiles (wide)', rows });
      return;
    }

    // ----- long format -----
    const lower = cls.lower;
    const skip = new Set([cls.idCol, cls.posCol]);
    if (cls.sampleCol >= 0) skip.add(cls.sampleCol);
    let depthCols;
    if (cls.sampleCol >= 0) {
      let d = lower.indexOf('depth');
      if (d < 0) d = lower.findIndex((x, j) => !skip.has(j) && /depth|count|reads/.test(x));
      if (d < 0) throw new Error(`${fileName}: found "pos" and "sample" but no "depth" column.`);
      depthCols = [{ j: d, name: null }];
    } else {
      depthCols = lower
        .map((x, j) => ({ x, j }))
        .filter(o => !skip.has(o.j) && /depth/.test(o.x))
        .map(o => ({ j: o.j, name: headers[o.j].replace(/[._-]?depth$/i, '') || stem }));
      if (!depthCols.length) throw new Error(`${fileName}: found "pos" but no depth column.`);
    }
    depthCols.forEach(dc => skip.add(dc.j));
    ['dist', 'norm', 'mean_norm', 'sem', 'n'].forEach(nm => { const j = lower.indexOf(nm); if (j >= 0) skip.add(j); });
    const metaIdx = headers.map((_, j) => j).filter(j => !skip.has(j));

    const recs = new Map();     // "sample\u0001idx" -> {sample, idx, arr}
    const posSeen = new Set();
    for (let i = h + 1; i < n; i++) {
      const line = cleanLine(lines[i]);
      if (!line) continue;
      const f = splitLine(line, sep);
      const id = get(f, cls.idCol);
      const pos = parseInt(f[cls.posCol], 10);
      if (!id || !Number.isFinite(pos) || pos < 0) continue;
      const k = B.idx(id);
      posSeen.add(pos);
      for (const dc of depthCols) {
        const sample = dc.name || get(f, cls.sampleCol) || stem;
        const key = sample + '\u0001' + k;
        let r = recs.get(key);
        if (!r) { r = { sample, idx: k, arr: new Float32Array(Math.max(160, pos + 1)) }; recs.set(key, r); B.addSample(sample); }
        if (pos >= r.arr.length) {
          const bigger = new Float32Array(Math.max(pos + 1, r.arr.length * 2));
          bigger.set(r.arr); r.arr = bigger;
        }
        const v = +f[dc.j];
        if (Number.isFinite(v)) r.arr[pos] += v;
      }
      if (metaIdx.length && !B.meta[k]) {
        const o = {};
        for (const j of metaIdx) o[headers[j]] = get(f, j);
        B.setMeta(k, o, false);
      }
      rows++;
      if ((i & 16383) === 0) await report(i);
    }
    const positions = [...posSeen].sort((a, b) => a - b);
    const chunk = { positions, entries: [] };
    for (const r of recs.values()) {
      const vals = new Float32Array(positions.length);
      for (let j = 0; j < positions.length; j++) vals[j] = r.arr[positions[j]] || 0;
      chunk.entries.push({ sample: r.sample, idx: r.idx, vals });
    }
    B.chunks.push(chunk);
    B.files.push({ name: fileName, kind: 'profiles (long)', rows });
  }

  // ---------- builder -> dataset ----------
  function finalize(B, opts = {}) {
    if (!B.chunks.length) {
      throw new Error('No coverage profiles found. Load the profiles file (seq.id, sample, p1 … p150) together with the sequences file.');
    }
    const posSet = new Set();
    for (const c of B.chunks) c.positions.forEach(p => posSet.add(p));
    const positions = [...posSet].sort((a, b) => a - b);
    const L = positions.length;
    const pIndex = new Map(positions.map((p, i) => [p, i]));

    const used = new Uint8Array(B.idList.length);
    for (const c of B.chunks) for (const e of c.entries) used[e.idx] = 1;
    const order = [];
    for (let k = 0; k < used.length; k++) if (used[k]) order.push(k);
    order.sort((x, y) => collator.compare(B.idList[x], B.idList[y]));
    const newIdx = new Int32Array(B.idList.length).fill(-1);
    order.forEach((old, ni) => { newIdx[old] = ni; });
    const N = order.length;

    const samples = B.sampleOrder.slice();
    const prof = {}, present = {};
    for (const s of samples) { prof[s] = new Float32Array(N * L); present[s] = new Uint8Array(N); }
    for (const c of B.chunks) {
      const map = c.positions.map(p => pIndex.get(p));
      for (const e of c.entries) {
        const ni = newIdx[e.idx];
        const P = prof[e.sample];
        const off = ni * L;
        for (let j = 0; j < map.length; j++) P[off + map[j]] += e.vals[j];
        present[e.sample][ni] = 1;
      }
    }
    return assemble({
      name: opts.name, isDemo: false, origin: opts.origin,
      positions, samples, prof, present,
      ids: order.map(o => B.idList[o]),
      meta: order.map(o => B.meta[o]),
      files: B.files.slice(),
      metaMissing: B.meta.every(m => m === null),
    });
  }

  // ---------- common post-processing ----------
  function assemble(o) {
    const { positions, samples, ids, prof } = o;
    const meta = o.meta || [];
    const N = ids.length, L = positions.length;
    const present = o.present || Object.fromEntries(samples.map(s => [s, new Uint8Array(N).fill(1)]));

    const total = {}, maxv = {}, peak = {};
    for (const s of samples) {
      const P = prof[s];
      const T = new Float64Array(N), M = new Float32Array(N), K = new Uint16Array(N);
      for (let i = 0; i < N; i++) {
        let t = 0, m = -Infinity, km = 0;
        const off = i * L;
        for (let k = 0; k < L; k++) {
          const v = P[off + k];
          t += v;
          if (v > m) { m = v; km = k; }
        }
        T[i] = t; M[i] = m > 0 ? m : 0; K[i] = km;
      }
      total[s] = T; maxv[s] = M; peak[s] = K;
    }

    // metadata columns, typed
    const colNames = [], seen = new Set();
    for (const m of meta) if (m) for (const k in m) if (!seen.has(k)) { seen.add(k); colNames.push(k); }
    const metaCols = [], mstr = {};
    let seqCol = null;
    for (const name of colNames) {
      const vals = new Array(N);
      const uniq = new Set();
      let numeric = true, nonEmpty = 0;
      for (let i = 0; i < N; i++) {
        let v = meta[i] ? meta[i][name] : undefined;
        v = v === undefined || v === null ? '' : String(v).trim();
        if (v === '' || v === 'NA' || v === 'NaN') v = 'NA';
        else { nonEmpty++; if (numeric && !Number.isFinite(Number(v))) numeric = false; }
        vals[i] = v;
        uniq.add(v);
      }
      if (!nonEmpty) continue;
      const lower = name.toLowerCase();
      let type;
      if (!seqCol && (lower === 'seq' || lower === 'sequence')) { type = 'seq'; seqCol = name; }
      else if (uniq.size <= 40) type = 'cat';
      else if (numeric) type = 'num';
      else type = 'text';
      const col = { name, type, numeric };
      if (type === 'cat') {
        col.levels = [...uniq].sort((a, b) => (a === 'NA' ? 1 : b === 'NA' ? -1 :
          numeric ? Number(a) - Number(b) : collator.compare(a, b)));
      }
      if (type === 'num') {
        let mn = Infinity, mx = -Infinity;
        for (const v of vals) if (v !== 'NA') { const x = Number(v); if (x < mn) mn = x; if (x > mx) mx = x; }
        col.min = mn; col.max = mx;
      }
      metaCols.push(col);
      mstr[name] = vals;
    }

    const searchCols = metaCols.filter(c => c.type !== 'num').map(c => mstr[c.name]);
    const hay = new Array(N);
    for (let i = 0; i < N; i++) {
      let s = ids[i];
      for (const v of searchCols) s += ' ' + v[i];
      hay[i] = s.toLowerCase();
    }

    return {
      name: o.name || 'Dataset', isDemo: !!o.isDemo, origin: o.origin,
      files: o.files || [], metaMissing: !!o.metaMissing,
      positions, L, N, samples, ids, metaCols,
      metaColMap: Object.fromEntries(metaCols.map(c => [c.name, c])),
      mstr, seqCol, prof, present, total, maxv, peak, hay,
    };
  }

  const api = { bytesToText, splitLine, detectSep, classify, createBuilder, parseText, finalize, assemble, fileStem };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.SPEData = api;
})(typeof window !== 'undefined' ? window : globalThis);
