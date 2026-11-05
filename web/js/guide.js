// Cómo se lee el mapa: un disco real con llamadas numeradas y tres figuras de uso.
import {
  D, t, onLang, fitCanvas, drawDiscBase, projectInto, magOf, MAG_RADIUS, CLASSES, CLASS_COLOR, CLASS_RES,
  DOMAINS, DOMAIN_COLOR, sprite, fmt, titleRes,
} from './core.js';

const $ = id => document.getElementById(id);

export function initGuide() {
  const fit = D.fits.catalytic;
  const ux = new Float32Array(D.n), uy = new Float32Array(D.n), hemi = new Uint8Array(D.n);
  projectInto(fit, 'equal', ux, uy, hemi);

  const inNorth = p => (D.posRows.get(p) || []).every(i => hemi[i] === 0);
  const centre = p => {
    const rows = D.posRows.get(p);
    return [d3.mean(rows, i => ux[i]), d3.mean(rows, i => uy[i])];
  };
  const spread = p => {
    const rows = D.posRows.get(p), [cx, cy] = centre(p);
    return d3.quantile(rows.map(i => Math.hypot(ux[i] - cx, uy[i] - cy)).sort(d3.ascending), .85);
  };
  // Ejemplos elegidos por los datos: nudo cerrado y nube abierta del disco norte
  const full = [...D.positions.values()].filter(p => p.n >= .9 * D.pdbs.length && inNorth(p.hmm_position));
  const near = D.positions.get(D.refPositions[0]);
  const [nx, ny] = centre(near.hmm_position);
  const apart = p => { const [x, y] = centre(p.hmm_position); return Math.hypot(x - nx, y - ny) > .4 && Math.hypot(x, y) > .4 && Math.hypot(x, y) < .86; };
  const tight = full.filter(p => p.distance > 12 && apart(p)).sort((a, b) => a.rmsf - b.rmsf)[0];
  const loose = full.filter(p => p.distance > 12 && p.rmsf <= 5 && apart(p)).sort((a, b) => b.rmsf - a.rmsf)[0];

  function dots(ctx, d, rows, { color, alpha = .5, radiusK = .7 } = {}) {
    ctx.globalCompositeOperation = 'lighter';
    for (const i of rows) {
      const m = magOf(fit.distance[i]);
      const sp = sprite(color ? color(i) : CLASS_COLOR[CLASSES[D.cls[i]]], Math.max(.5, MAG_RADIUS[m] * radiusK), d.dpr);
      ctx.globalAlpha = (m < 2 ? 1 : .7) * alpha;
      ctx.drawImage(sp, d.x + ux[i] * d.R * d.k + d.ox - sp.half, d.y + uy[i] * d.R * d.k + d.oy - sp.half, sp.css, sp.css);
    }
    ctx.globalAlpha = 1;
    ctx.globalCompositeOperation = 'source-over';
  }
  const northRows = [];
  for (let i = 0; i < D.n; i++) if (!hemi[i]) northRows.push(i);

  /* Disco grande con llamadas */
  const cv = $('guide-canvas');
  function drawGuide() {
    const { ctx, w, h, dpr } = fitCanvas(cv);
    ctx.clearRect(0, 0, w, h);
    if (w < 60) return;
    const R = (w / 2 - (w > 540 ? 28 : 21)) / 1.035;
    const d = { x: w / 2, y: h / 2, R, k: 1, ox: 0, oy: 0, dpr };
    drawDiscBase(ctx, d.x, d.y, R, { hemi: 'n' });
    dots(ctx, d, northRows, { alpha: .5, radiusK: Math.min(.8, w / 800) });

    const mark = (n, x, y, ring = 0) => {
      ctx.strokeStyle = '#ece4d2'; ctx.lineWidth = 1;
      if (ring) { ctx.setLineDash([1, 3]); ctx.lineCap = 'round'; ctx.beginPath(); ctx.arc(x, y, ring, 0, 7); ctx.stroke(); ctx.setLineDash([]); }
      const a = -Math.PI / 4, off = ring ? ring : 0;
      const bx = x + Math.cos(a) * (off + 16), by = y + Math.sin(a) * (off + 16);
      if (off || n !== 1) { ctx.beginPath(); ctx.moveTo(x + Math.cos(a) * off, y + Math.sin(a) * off); ctx.lineTo(bx, by); ctx.stroke(); }
      const px = n === 1 ? x : bx, py = n === 1 ? y : by;
      ctx.fillStyle = '#ece4d2';
      ctx.beginPath(); ctx.arc(px, py, 11, 0, 7); ctx.fill();
      ctx.fillStyle = '#080c16'; ctx.font = 'italic 600 15px "Bodoni Moda", serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillText(String(n), px, py + 1);
    };
    mark(1, d.x, d.y);
    const a2 = Math.PI * .86;   // lugar del limbo sin rótulo de grados
    mark(2, d.x + Math.cos(a2) * R, d.y + Math.sin(a2) * R);
    const [x3, y3] = centre(near.hmm_position);
    mark(3, d.x + x3 * R, d.y + y3 * R, 13);
    if (tight) { const [x, y] = centre(tight.hmm_position); mark(4, d.x + x * R, d.y + y * R, Math.max(9, spread(tight.hmm_position) * R + 6)); }
    if (loose) { const [x, y] = centre(loose.hmm_position); mark(5, d.x + x * R, d.y + y * R, Math.max(14, spread(loose.hmm_position) * R + 6)); }
  }

  /* Figuras de uso */
  function frameFor(rows, w, R) {
    const xs = d3.extent(rows, i => ux[i]), ys = d3.extent(rows, i => uy[i]);
    const k = Math.min(9, .62 * 2 / Math.max(xs[1] - xs[0], ys[1] - ys[0], .02));
    return { k, ox: -(xs[0] + xs[1]) / 2 * R * k, oy: -(ys[0] + ys[1]) / 2 * R * k };
  }
  function drawFilter() {
    const { ctx, w, h, dpr } = fitCanvas($('use-1'));
    ctx.clearRect(0, 0, w, h);
    if (w < 60) return;
    const R = w / 2 - 10, d = { x: w / 2, y: h / 2, R, k: 1, ox: 0, oy: 0, dpr };
    drawDiscBase(ctx, d.x, d.y, R, { hemi: 'n', labels: false });
    const rows = northRows.filter(i => CLASS_RES.arom.includes(D.residue[i]) && fit.distance[i] < 12);
    dots(ctx, d, rows, { alpha: .9, radiusK: 1 });
  }
  function drawClick() {
    const { ctx, w, h, dpr } = fitCanvas($('use-2'));
    ctx.clearRect(0, 0, w, h);
    if (w < 60) return;
    const R = w / 2 - 10, d = { x: w / 2, y: h / 2, R, k: 1, ox: 0, oy: 0, dpr };
    drawDiscBase(ctx, d.x, d.y, R, { hemi: 'n', labels: false });
    const Q = D.struct(D.fam.reference_pdb);
    const sel = Q.pos.indexOf(D.refPositions[0]);
    ctx.globalCompositeOperation = 'lighter';
    for (let i = 0; i < Q.n; i++) {
      if (Q.hemi[i]) continue;
      ctx.globalAlpha = .75; ctx.fillStyle = CLASS_COLOR[CLASSES[Q.cls[i]]];
      ctx.beginPath(); ctx.arc(d.x + Q.ux[i] * R, d.y + Q.uy[i] * R, MAG_RADIUS[magOf(Q.distance[i])], 0, 6.3); ctx.fill();
    }
    ctx.globalAlpha = 1; ctx.globalCompositeOperation = 'source-over';
    const P = i => [d.x + Q.ux[i] * R, d.y + Q.uy[i] * R];
    const [sx, sy] = P(sel);
    ctx.strokeStyle = '#d8b46a'; ctx.lineWidth = 1;
    for (let r = 0; r < Q.n; r++) {
      if (r === sel || Q.hemi[r]) continue;
      if (Math.hypot(Q.x[r] - Q.x[sel], Q.y[r] - Q.y[sel], Q.z[r] - Q.z[sel]) > 7) continue;
      const [x, y] = P(r);
      ctx.beginPath(); ctx.moveTo(sx, sy); ctx.lineTo(x, y); ctx.stroke();
      ctx.beginPath(); ctx.arc(x, y, 3.5, 0, 7); ctx.stroke();
    }
    ctx.strokeStyle = '#fff8e8';
    ctx.beginPath(); ctx.arc(sx, sy, 6, 0, 7);
    for (const [a, b] of [[0, -1], [0, 1], [-1, 0], [1, 0]]) { ctx.moveTo(sx + a * 9, sy + b * 9); ctx.lineTo(sx + a * 15, sy + b * 15); }
    ctx.stroke();
  }
  function drawCompare() {
    const { ctx, w, h, dpr } = fitCanvas($('use-3'));
    ctx.clearRect(0, 0, w, h);
    if (w < 60) return;
    const R = w / 2 - 10;
    const rows = (D.posRows.get(333) || []).filter(i => D.residue[i] === 'TRP' && !hemi[i]);
    if (!rows.length) return;
    const f = frameFor(rows, w, R);
    const d = { x: w / 2, y: h / 2, R, dpr, ...f };
    // Retícula local del zoom
    ctx.save();
    ctx.beginPath(); ctx.arc(w / 2, h / 2, R, 0, 7); ctx.clip();
    drawDiscBase(ctx, d.x + d.ox, d.y + d.oy, R * d.k, { hemi: 'n', labels: false, rim: false });
    ctx.restore();
    ctx.strokeStyle = 'rgba(216,180,106,.75)'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(w / 2, h / 2, R, 0, 7); ctx.stroke();
    dots(ctx, d, rows, { color: i => DOMAIN_COLOR[DOMAINS[D.dom[i]]], alpha: 1, radiusK: 1.9 });
  }

  /* Cierre: el disco sur, con todos sus residuos */
  function drawFinale() {
    const { ctx, w, h, dpr } = fitCanvas($('finale-canvas'));
    ctx.clearRect(0, 0, w, h);
    if (w < 60) return;
    const small = w < 861;
    const R = small ? w * .62 : Math.min(h * .78, w * .34);
    const d = { x: small ? w * .5 : w - R * .78, y: small ? h + R * .12 : h * .5, R, k: 1, ox: 0, oy: 0, dpr };
    drawDiscBase(ctx, d.x, d.y, R, { hemi: 's', labels: false });
    const south = [];
    for (let i = 0; i < D.n; i++) if (hemi[i]) south.push(i);
    dots(ctx, d, south, { alpha: .42, radiusK: .66 });
  }

  function steps() {
    $('how-steps').innerHTML = [1, 2, 3, 4, 5].map(i => `<li><b>${t(`h${i}_t`)}</b> ${t(`h${i}_p`)}</li>`).join('');
    const keys = ['fam', 'dist', 'res', 'pos', 'prot', 'dom', 'proj', 'fit', 'color', 'spec', 'zoom', 'table'];
    $('controls').innerHTML = keys.map(k => `<dt>${t('c2_' + k)}</dt><dd>${t('c2_' + k + '_p', { a: D.refPositions[0], b: D.refPositions[1], p: D.pdbs.length })}</dd>`).join('');
  }

  function copy() {
    steps();
    document.querySelectorAll('[data-t="fin_p"]').forEach(el => { el.textContent = t('fin_p', { p: D.pdbs.length }); });
    $('k3-p').textContent = t('k3_p', { res: titleRes(near.top_residue), pos: near.hmm_position, d: fmt(near.distance) });
    if (tight) $('k4-p').textContent = t('k4_p', { p: tight.hmm_position, res: titleRes(tight.top_residue), r: fmt(tight.rmsf, 2) });
    if (loose) $('k5-p').textContent = t('k5_p', { p: loose.hmm_position, r: fmt(loose.rmsf, 2) });
    $('k6-leg').innerHTML = CLASSES.map(c => `<span class="k" style="--c:${CLASS_COLOR[c]}"><i></i>${t('cls_' + c)}</span>`).join('');
  }
  const drawAll = () => { drawGuide(); drawFilter(); drawClick(); drawCompare(); drawFinale(); };
  copy();
  onLang(copy);
  drawAll();
  document.fonts?.ready.then(drawAll);
  let rt = 0;
  new ResizeObserver(() => { clearTimeout(rt); rt = setTimeout(drawAll, 80); }).observe($('read'));
}
