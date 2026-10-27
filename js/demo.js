/* =========================================================
   demo.js — simulated dataset so the public site works
   without any real data. Deterministic (seeded), generated
   in the browser, so nothing unpublished is in the repo.
   ========================================================= */
(function (root) {
  'use strict';

  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function make() {
    const rnd = mulberry32(20261027);
    const L = 150, N = 900, START = 100;
    const BASES = 'ACGT';
    const COMP = { A: 'T', C: 'G', G: 'C', T: 'A' };
    const randSeq = n => { let s = ''; for (let i = 0; i < n; i++) s += BASES[(rnd() * 4) | 0]; return s; };
    const revcomp = s => s.split('').reverse().map(c => COMP[c]).join('');
    const normal = () => { let u = 0; while (!u) u = rnd(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rnd()); };
    const pick = list => { let r = rnd(), acc = 0; for (const [v, w] of list) { acc += w; if (r < acc) return v; } return list[list.length - 1][0]; };

    const samples = ['DMSO', 'PhenDC3_1uM', 'PhenDC3_5uM'];
    const dose = [0, 1, 2];
    const positions = Array.from({ length: L }, (_, k) => k + 1);
    const prof = {};
    samples.forEach(s => { prof[s] = new Float32Array(N * L); });
    const ids = [], meta = [];

    for (let i = 0; i < N; i++) {
      const type = pick([['STR', 0.55], ['G4', 0.14], ['I_motif', 0.1], ['Hairpin', 0.11], ['Control', 0.1]]);
      let insert, motif = 'NA', mlen = 'NA', struct = 'none', base, offset;
      if (type === 'STR') {
        const ml = 1 + ((rnd() * 6) | 0);
        do { motif = randSeq(ml); } while (ml > 1 && /^(.)\1*$/.test(motif));
        const len = rnd() < 0.5 ? 24 : 48;
        insert = motif.repeat(Math.ceil(len / ml)).slice(0, len);
        mlen = String(ml);
        const gc = (motif.match(/[GC]/g) || []).length / ml;
        const gRun = /GG|CC/.test(motif + motif) ? 0.12 : 0;
        base = 0.04 + 0.28 * gc * (ml <= 3 ? 1 : 0.6) + gRun + (len === 48 ? 0.06 : 0) + 0.08 * rnd();
        offset = -2;
      } else if (type === 'G4') {
        const loop = () => randSeq(1 + ((rnd() * 3) | 0)).replace(/G/g, 'T');
        insert = 'GGG' + loop() + 'GGG' + loop() + 'GGG' + loop() + 'GGG';
        struct = 'G4'; base = 0.22 + 0.2 * rnd(); offset = 3;
      } else if (type === 'I_motif') {
        const loop = () => randSeq(1 + ((rnd() * 3) | 0)).replace(/C/g, 'A');
        insert = 'CCC' + loop() + 'CCC' + loop() + 'CCC' + loop() + 'CCC';
        struct = 'i-motif'; base = 0.14 + 0.14 * rnd(); offset = 1;
      } else if (type === 'Hairpin') {
        const stem = randSeq(8 + ((rnd() * 5) | 0));
        insert = stem + 'GAAA' + revcomp(stem);
        struct = 'hairpin'; base = 0.18 + 0.22 * rnd(); offset = 0;
      } else {
        insert = randSeq(30); base = 0.02 + 0.04 * rnd(); offset = 0;
      }
      const seq = (randSeq(START - 1) + insert + randSeq(L)).slice(0, L);
      const id = `Seq_${String(i + 1).padStart(5, '0')}_${type}${type === 'STR' ? '_' + motif : ''}`;
      ids.push(id);
      meta.push({
        seq, 'motif.type': type, 'motif.length': mlen, 'STR.motif.ssDNA': motif,
        'insert.length': String(insert.length), 'pred.struct': struct,
      });

      const centre = START + offset + Math.round(normal() * 1.5);
      const depth = Math.exp(5.4 + 0.9 * normal());
      samples.forEach((s, si) => {
        const eff = type === 'G4' ? 1 + 0.85 * dose[si] : type === 'I_motif' ? 1 + 0.08 * dose[si] : 1;
        const p = Math.min(0.9, base * eff);
        const T = Math.max(1, Math.round(depth * (0.8 + 0.4 * rnd())));
        const P = prof[s];
        const off = i * L;
        for (let r = 0; r < T; r++) {
          const u = rnd();
          let pos;
          if (u < p) pos = centre + Math.round(normal() * 2.2);
          else if (u < p + 0.55 * (1 - p)) pos = L;
          else pos = 1 + ((rnd() * L) | 0);
          if (pos < 1) pos = 1; else if (pos > L) pos = L;
          P[off + pos - 1] += 1;
        }
      });
    }

    return root.SPEData.assemble({
      name: 'Simulated demo dataset', isDemo: true, origin: START,
      positions, samples, ids, meta, prof,
    });
  }

  root.SPEDemo = { make };
})(typeof window !== 'undefined' ? window : globalThis);
