// Portada: un cielo que respira. Cada estrella es una posición de la familia y
// viaja sin pausas por los lugares que ocupa en cada proteína. Sin decoración: son datos.
import {
  D, t, lang, onLang, REDUCED, fitCanvas, drawDiscBase, radialFn, dirToDisc, magOf, MAG_RADIUS,
  CLASSES, CLASS_COLOR, sprite, engravedStar, spanningTree, clamp01, fmt, fmtInt, titleRes,
} from './core.js';

const STEP = 2.6;               // segundos por proteína, sin paradas
const KEY_LIMIT = 8;            // las constelaciones unen posiciones a menos de 8 Å
const LABELS = [165, 360, 333, 406, 121, 295];
const RIM_FADE = .06;           // una estrella se apaga al cruzar de un disco al otro

// Radio continuo según la distancia: sin saltos entre magnitudes
const RAD_D = [6, 9, 13, 18, 25, 40], RAD_R = [2.9, 2.2, 1.65, 1.25, .95, .72];
function radiusAt(d) {
  if (d <= RAD_D[0]) return RAD_R[0];
  for (let i = 1; i < RAD_D.length; i++) {
    if (d < RAD_D[i]) return RAD_R[i - 1] + (RAD_R[i] - RAD_R[i - 1]) * (d - RAD_D[i - 1]) / (RAD_D[i] - RAD_D[i - 1]);
  }
  return RAD_R[RAD_R.length - 1];
}
// Spline de Catmull-Rom: pasa por cada proteína y no se detiene en ella
const cr = (a, b, c, d, f) => .5 * (2 * b + (c - a) * f + (2 * a - 5 * b + 4 * c - d) * f * f + (3 * b - a - 3 * c + d) * f * f * f);

export function initHero() {
  const cv = document.getElementById('hero-canvas');
  const baseCv = document.getElementById('hero-base');
  const section = document.getElementById('top');
  const readout = document.getElementById('hero-readout');
  const nowEl = document.getElementById('hero-now');
  const pauseBtn = document.getElementById('hero-pause');
  const rad = radialFn('equal');
  const fit = D.fits.catalytic;
  const P = D.pdbs.length;
  const posList = [...D.posRows.keys()].sort((a, b) => a - b);
  const NP = posList.length;

  // Coordenadas por posición y proteína
  const raw = new Float32Array(NP * P * 3).fill(NaN);
  const rawC = new Uint8Array(NP * P);
  posList.forEach((pos, pi) => {
    for (const i of D.posRows.get(pos)) {
      const o = (pi * P + D.pdbI[i]) * 3;
      raw[o] = fit.x[i]; raw[o + 1] = fit.y[i]; raw[o + 2] = fit.z[i];
      rawC[pi * P + D.pdbI[i]] = D.cls[i];
    }
  });
  // Orden de paso: cada proteína va seguida de la más parecida, para que el cielo fluya
  const gap = (a, b) => {
    let s = 0, n = 0;
    for (let pi = 0; pi < NP; pi += 3) {
      const u = (pi * P + a) * 3, v = (pi * P + b) * 3;
      if (Number.isNaN(raw[u]) || Number.isNaN(raw[v])) continue;
      s += Math.hypot(raw[u] - raw[v], raw[u + 1] - raw[v + 1], raw[u + 2] - raw[v + 2]); n++;
    }
    return n > 20 ? s / n : 1e6;
  };
  const left = new Set(D.pdbs.map((_, i) => i));
  const start = D.pdbIdx.get(D.fam.reference_pdb) ?? 0;
  const order = [start]; left.delete(start);
  while (left.size) {
    const cur0 = order[order.length - 1];
    let best = -1, bd = Infinity;
    for (const c of left) { const g = gap(cur0, c); if (g < bd) { bd = g; best = c; } }
    order.push(best); left.delete(best);
  }
  const M = new Float32Array(NP * P * 3).fill(NaN);
  const C = new Uint8Array(NP * P);
  for (let pi = 0; pi < NP; pi++) for (let k = 0; k < P; k++) {
    const s = (pi * P + order[k]) * 3, d = (pi * P + k) * 3;
    M[d] = raw[s]; M[d + 1] = raw[s + 1]; M[d + 2] = raw[s + 2];
    C[pi * P + k] = rawC[pi * P + order[k]];
  }

  const keyIdx = posList.map((p, i) => (D.positions.get(p)?.distance < KEY_LIMIT && D.positions.get(p)?.n >= 0.9 * P) ? i : -1).filter(i => i >= 0);
  const centroid = pi => {
    let x = 0, y = 0, z = 0, n = 0;
    for (let k = 0; k < P; k++) { const o = (pi * P + k) * 3; if (!Number.isNaN(M[o])) { x += M[o]; y += M[o + 1]; z += M[o + 2]; n++; } }
    return [x / n, y / n, z / n];
  };
  const groups = { n: [], s: [] };
  keyIdx.forEach(pi => {
    const c = centroid(pi), q = dirToDisc(c[0], c[1], c[2], rad);
    groups[q.north ? 'n' : 's'].push({ pi, u: q.u, v: q.v });
  });
  const trees = { n: spanningTree(groups.n), s: spanningTree(groups.s) };
  const labelIdx = LABELS.map(p => posList.indexOf(p)).filter(i => i >= 0);

  let L, dpr, W, H, running = !REDUCED, visible = true, raf = 0;
  let clock = 0, last = 0, hover = null;
  const cur = new Array(NP);

  function layout() {
    const f = fitCanvas(cv); W = f.w; H = f.h; dpr = f.dpr;
    const small = W < 861;
    if (small) {
      // En pantalla estrecha el cielo va en su propio hueco, entre el texto y el pie: nunca debajo del texto.
      const box = section.getBoundingClientRect();
      const top = section.querySelector('.hero-copy').getBoundingClientRect().bottom - box.top;
      const foot = section.querySelector('.hero-foot').getBoundingClientRect().top - box.top;
      const R = Math.max(80, Math.min(W * .41, (foot - top) / 2 - 32));
      const main = { x: W * .5, y: (top + foot) / 2, R };
      L = { main, n: main, s: null };
    } else {
      // Se reserva sitio abajo para que su rótulo de 180° no se corte.
      const y = H * .52;
      const R = Math.min(H * .44, W * .31, (H - y - 28) / 1.035);
      const main = { x: W - R * 1.12 - Math.max(16, W * .02), y, R };
      const R2 = R * (H < 800 ? .24 : .285);
      const y2 = H - 58 - R2 * 1.14, gp = (R + R2) * 1.06, dy = y2 - main.y;
      L = { main, n: main, s: { x: main.x - Math.sqrt(Math.max(0, gp * gp - dy * dy)), y: y2, R: R2 } };
    }
    buildBase();
  }

  // Lienzo fijo aparte: limbo, retícula y polvo. Se pinta una vez, no en cada fotograma.
  function buildBase() {
    const { ctx: c } = fitCanvas(baseCv);
    c.clearRect(0, 0, W, H);
    drawDiscBase(c, L.n.x, L.n.y, L.n.R, { hemi: 'n' });
    if (L.s) drawDiscBase(c, L.s.x, L.s.y, L.s.R, { hemi: 's', labels: false });
    c.globalCompositeOperation = 'lighter';
    const { theta, phi, distance } = fit;
    for (let i = 0; i < D.n; i++) {
      const north = theta[i] < Math.PI / 2;
      const disc = north ? L.n : L.s; if (!disc) continue;
      const r = rad(north ? theta[i] : Math.PI - theta[i]);
      const m = magOf(distance[i]);
      const k = disc.R / L.main.R;
      const sp = sprite(CLASS_COLOR[CLASSES[D.cls[i]]], Math.max(.5, MAG_RADIUS[m] * (k < 1 ? .5 : .68)), dpr);
      c.globalAlpha = (m < 2 ? .6 : .36) * Math.min(1, disc.R / 170);
      c.drawImage(sp, disc.x + (north ? -1 : 1) * r * Math.sin(phi[i]) * disc.R - sp.half, disc.y - r * Math.cos(phi[i]) * disc.R - sp.half, sp.css, sp.css);
    }
    c.globalAlpha = 1;
    c.globalCompositeOperation = 'source-over';
  }

  function sampleSky(tt) {
    const base = Math.floor(tt), f = tt - base;
    const k0 = base % P, k1 = (k0 + 1) % P, km = (k0 + P - 1) % P, k2 = (k0 + 2) % P;
    for (let pi = 0; pi < NP; pi++) {
      const o = pi * P;
      const a = (o + k0) * 3, b = (o + k1) * 3;
      const ha = !Number.isNaN(M[a]), hb = !Number.isNaN(M[b]);
      let x, y, z, al = 1;
      if (ha && hb) {
        let m = (o + km) * 3, n = (o + k2) * 3;
        if (Number.isNaN(M[m])) m = a;
        if (Number.isNaN(M[n])) n = b;
        x = cr(M[m], M[a], M[b], M[n], f); y = cr(M[m + 1], M[a + 1], M[b + 1], M[n + 1], f); z = cr(M[m + 2], M[a + 2], M[b + 2], M[n + 2], f);
      } else if (ha) { x = M[a]; y = M[a + 1]; z = M[a + 2]; al = 1 - f; }
      else if (hb) { x = M[b]; y = M[b + 1]; z = M[b + 2]; al = f; }
      else { cur[pi] = null; continue; }
      const q = dirToDisc(x, y, z, rad);
      const disc = q.north ? L.n : L.s;
      if (!disc) { cur[pi] = null; continue; }
      al *= clamp01(Math.abs(y) / q.d / RIM_FADE);
      cur[pi] = {
        x: disc.x + q.u * disc.R, y: disc.y + q.v * disc.R, north: q.north, d: q.d, al, scale: disc.R / L.main.R,
        c0: C[o + (ha ? k0 : k1)], c1: C[o + (hb ? k1 : k0)], f, pi,
      };
    }
    return Math.round(tt) % P;
  }

  function resOf(pi, k) {
    const pos = posList[pi];
    for (const i of D.posRows.get(pos)) if (D.pdbI[i] === order[k]) return D.residue[i];
    return D.positions.get(pos)?.top_residue || '';
  }

  function dot(ctx, s, r, color, alpha) {
    ctx.globalAlpha = alpha; ctx.fillStyle = color;
    ctx.beginPath(); ctx.arc(s.x, s.y, r, 0, 6.3); ctx.fill();
  }
  function draw() {
    const ctx = cv.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    const k = sampleSky(clock);

    // Constelaciones sobre las estrellas reales
    ctx.lineWidth = 1;
    ctx.strokeStyle = '#d8b46a';
    for (const h of ['n', 's']) {
      for (const [a, b] of trees[h]) {
        const p = cur[groups[h][a].pi], q = cur[groups[h][b].pi];
        if (!p || !q || p.north !== q.north) continue;
        ctx.globalAlpha = .62 * Math.min(p.al, q.al);
        ctx.beginPath(); ctx.moveTo(p.x, p.y); ctx.lineTo(q.x, q.y); ctx.stroke();
      }
    }

    // Estrellas. El tipo de residuo cambia con un fundido; el radio cambia de forma continua.
    const big = Math.min(1.3, Math.max(.8, L.main.R / 340));
    for (let pi = 0; pi < NP; pi++) {
      const s = cur[pi]; if (!s) continue;
      const e = clamp01((10.5 - s.d) / 2.5);        // peso del dibujo grabado frente al punto
      const col0 = CLASS_COLOR[CLASSES[s.c0]], col1 = CLASS_COLOR[CLASSES[s.c1]];
      const same = s.c0 === s.c1;
      if (e < 1) {
        const r = radiusAt(s.d) * (s.scale < 1 ? .9 : 1.5);
        const a = s.al * (1 - e);
        if (same) dot(ctx, s, r, col0, a);
        else { dot(ctx, s, r, col0, a * (1 - s.f)); dot(ctx, s, r, col1, a * s.f); }
      }
      if (e > 0) {
        const r = Math.max(5, Math.min(11.5, 5 + (10 - s.d) * 1.2)) * (s.scale < 1 ? .62 : 1) * big;
        const a = s.al * e;
        if (same) engravedStar(ctx, s.x, s.y, r, col0, a);
        else { engravedStar(ctx, s.x, s.y, r, col0, a * (1 - s.f)); engravedStar(ctx, s.x, s.y, r, col1, a * s.f); }
        engravedStar(ctx, s.x, s.y, r * .34, '#fff8e8', a);
      }
    }
    ctx.globalAlpha = 1;

    // Rótulos de atlas
    ctx.fillStyle = '#ece4d2';
    ctx.textBaseline = 'middle';
    const fs = Math.max(13, Math.min(17, L.main.R / 21));
    for (const pi of labelIdx) {
      const s = cur[pi];
      if (!s || (s.scale < 1 && L.s.R < 120)) continue;
      ctx.font = `italic 500 ${s.scale < 1 ? 12 : fs}px "Bodoni Moda", serif`;
      const lab = `${titleRes(resOf(pi, k))} ${posList[pi]}`;
      const flip = s.x + 11 + ctx.measureText(lab).width > W - 8;
      ctx.textAlign = flip ? 'right' : 'left';
      ctx.globalAlpha = s.al;
      ctx.fillText(lab, s.x + (flip ? -11 : 11), s.y - 9);
    }
    ctx.globalAlpha = 1;

    if (hover) {
      ctx.strokeStyle = '#ece4d2'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(hover.x, hover.y, 8, 0, 7); ctx.stroke();
    }

    const code = D.pdbs[order[k]];
    if (code !== draw.code || lang !== draw.lang) {
      draw.code = code; draw.lang = lang;
      const org = D.organism(code);
      nowEl.innerHTML = `${t('hero_now')}: ${code}${org ? ' <i>' + org + '</i>' : ''}`;
    }
  }

  function frame(ts) {
    raf = 0;
    const dt = Math.min(.05, (ts - last) / 1000 || 0); last = ts;
    clock += dt / STEP;
    draw();
    if (running && visible) raf = requestAnimationFrame(frame);
  }
  function play() { if (!raf && running && visible) { last = performance.now(); raf = requestAnimationFrame(frame); } }

  function caption() {
    document.getElementById('hero-caption').textContent = t('hero_caption', { p: P, r: fmtInt(D.n) });
    pauseBtn.textContent = t(running ? 'hero_pause' : 'hero_play');
    pauseBtn.setAttribute('aria-pressed', String(!running));
    if (!running) draw();
  }

  pauseBtn.addEventListener('click', () => { running = !running; caption(); play(); });
  new IntersectionObserver(es => { visible = es[0].isIntersecting; play(); }, { threshold: .02 }).observe(section);
  document.addEventListener('visibilitychange', () => { visible = !document.hidden; play(); });

  // Lectura al pasar por una estrella
  section.addEventListener('pointermove', e => {
    if (e.pointerType === 'touch') return;
    const r = cv.getBoundingClientRect();
    const x = e.clientX - r.left, y = e.clientY - r.top;
    let best = null, bd = 18;
    for (const s of cur) { if (!s) continue; const d = Math.hypot(s.x - x, s.y - y); if (d < bd) { bd = d; best = s; } }
    hover = best;
    if (best) {
      const pos = posList[best.pi];
      const info = D.positions.get(pos);
      readout.hidden = false;
      readout.textContent = `${titleRes(info?.top_residue || '')} ${pos} · ${fmt(info?.distance)} Å · RMSF ${fmt(info?.rmsf, 2)} Å`;
      readout.style.left = Math.min(W - 230, best.x + 14) + 'px';
      readout.style.top = (best.y + 14) + 'px';
    } else readout.hidden = true;
    if (!running) draw();
  });
  section.addEventListener('pointerleave', () => { hover = null; readout.hidden = true; if (!running) draw(); });

  let rt = 0;
  const ro = new ResizeObserver(() => { clearTimeout(rt); rt = setTimeout(() => { layout(); draw(); }, 80); });
  [section, section.querySelector('.hero-copy'), section.querySelector('.hero-foot')].forEach(el => ro.observe(el));

  layout();
  caption();
  draw();
  document.fonts?.ready.then(() => { layout(); draw(); });
  onLang(() => { caption(); draw(); });
  play();
}
