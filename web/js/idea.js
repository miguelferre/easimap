// La idea: del enzima al mapa con el scroll. Un solo lienzo, cinco estados.
// s va de 0 a 4: proteína, bolsillo, esfera, discos, familia.
import {
  D, t, onLang, REDUCED, fitCanvas, drawDiscBase, radialFn, dirToDisc, magOf, CLASSES, CLASS_COLOR,
  MAG_RADIUS, sprite, ease, clamp01, fmtInt,
} from './core.js';

export function initIdea() {
  const cv = document.getElementById('idea-canvas');
  const stage = cv.parentElement;
  const stepsEl = document.getElementById('idea-steps');
  const steps = [...stepsEl.children];
  const tag = document.getElementById('idea-tag');
  const rad = radialFn('equal');
  const fit = D.fits.catalytic;
  const R0 = D.struct(D.fam.reference_pdb);   // la referencia: 3SCU, fuera del conjunto actual (D-012)
  const rows = d3.range(R0.n).sort((p, q) => (R0.auth[p] ?? 0) - (R0.auth[q] ?? 0));
  const N = rows.length;
  const px = rows.map(i => R0.x[i]), py = rows.map(i => R0.y[i]), pz = rows.map(i => R0.z[i]);
  const dist = rows.map(i => R0.distance[i]);
  const cls = rows.map(i => R0.cls[i]);
  const disc = rows.map(i => dirToDisc(R0.x[i], R0.y[i], R0.z[i], rad));
  const lig = (D.extras.ligand || []).filter(a => a.element !== 'H');
  const bonds = [];
  for (let i = 0; i < lig.length; i++) for (let j = i + 1; j < lig.length; j++) {
    if (Math.hypot(lig[i].x - lig[j].x, lig[i].y - lig[j].y, lig[i].z - lig[j].z) < 1.75) bonds.push([i, j]);
  }
  // Centro de la proteína: la cámara arranca aquí y viaja al azúcar
  const com = [px, py, pz].map(a => a.reduce((s, v) => s + v, 0) / N);

  let W, H, dpr, dust = null, s = 0, target = 0, raf = 0, yaw = -.6, spin = 0, visible = false, lastTs = 0;
  const Rs = 25;

  function size() {
    const f = fitCanvas(cv); W = f.w; H = f.h; dpr = f.dpr;
    buildDust();
  }
  function discs() {
    const small = W < 500;
    const R = Math.min(W / (small ? 4.5 : 4.75), H * .36);
    const gap = R * (small ? 1.1 : 1.2);
    return { n: { x: W / 2 - gap, y: H / 2 - 8, R }, s: { x: W / 2 + gap, y: H / 2 - 8, R } };
  }
  function buildDust() {
    const L = discs();
    dust = document.createElement('canvas');
    dust.width = cv.width; dust.height = cv.height;
    const c = dust.getContext('2d');
    c.setTransform(dpr, 0, 0, dpr, 0, 0);
    c.globalCompositeOperation = 'lighter';
    const { theta, phi, distance } = fit;
    for (let i = 0; i < D.n; i++) {
      const north = theta[i] < Math.PI / 2;
      const r = rad(north ? theta[i] : Math.PI - theta[i]);
      const d = north ? L.n : L.s;
      const m = magOf(distance[i]);
      const sp = sprite(CLASS_COLOR[CLASSES[D.cls[i]]], Math.max(.55, MAG_RADIUS[m] * .62), dpr);
      c.globalAlpha = (m < 2 ? .55 : .3) * Math.min(1, L.n.R / 190);
      c.drawImage(sp, d.x + (north ? -1 : 1) * r * Math.sin(phi[i]) * d.R - sp.half, d.y - r * Math.cos(phi[i]) * d.R - sp.half, sp.css, sp.css);
    }
  }

  // Cámara: rotación en yaw y tilt, perspectiva suave
  function cam(a, tilt, zoom, focus) {
    const ca = Math.cos(a), sa = Math.sin(a), ct = Math.cos(tilt), st = Math.sin(tilt);
    const base = Math.min(W, H) / 96 * zoom;
    return (x, y, z) => {
      x -= focus[0]; y -= focus[1]; z -= focus[2];
      const x1 = ca * x + sa * z, z1 = -sa * x + ca * z;
      const y1 = ct * y - st * z1, z2 = st * y + ct * z1;
      const p = 1 / (1 - z2 / 190);
      return [W / 2 + x1 * base * p, H / 2 - y1 * base * p, z2, p];
    };
  }

  function draw() {
    const ctx = cv.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    const s1 = ease(clamp01(s)), s2 = ease(clamp01(s - 1)), s3 = ease(clamp01(s - 2)), s4 = ease(clamp01(s - 3));
    const zoomFit = .84 * 48 / Rs;
    const zoom = (1 + s1 * .75) * (1 - s2) + zoomFit * s2;
    const focus = [com[0] * (1 - s1), com[1] * (1 - s1), com[2] * (1 - s1)];
    const P = cam(yaw + spin, .26, zoom, focus);
    const L = discs();

    // Posición de cada residuo en pantalla
    const scr = new Array(N);
    for (let i = 0; i < N; i++) {
      const d = dist[i] || 1;
      const k = 1 + (Rs / d - 1) * s2;
      const q = P(px[i] * k, py[i] * k, pz[i] * k);
      const dd = disc[i], dc = dd.north ? L.n : L.s;
      const tx = dc.x + dd.u * dc.R, ty = dc.y + dd.v * dc.R;
      scr[i] = [q[0] + (tx - q[0]) * s3, q[1] + (ty - q[1]) * s3, q[2], q[3]];
    }

    // Cadena
    const traceA = (1 - s2) * (1 - s1 * .35);
    if (traceA > .01) {
      ctx.globalAlpha = traceA * .5; ctx.strokeStyle = '#ece4d2'; ctx.lineWidth = 1;
      ctx.beginPath();
      for (let i = 0; i < N; i++) {
        if (i && (R0.auth[rows[i]] - R0.auth[rows[i - 1]]) > 3) ctx.moveTo(scr[i][0], scr[i][1]);
        else if (i) ctx.lineTo(scr[i][0], scr[i][1]); else ctx.moveTo(scr[i][0], scr[i][1]);
      }
      ctx.stroke(); ctx.globalAlpha = 1;
    }

    // Esfera
    const sphA = s2 * (1 - s3);
    const sphR = Rs * Math.min(W, H) / 96 * zoom;
    if (sphA > .01) {
      ctx.save(); ctx.globalAlpha = sphA;
      ctx.strokeStyle = 'rgba(216,180,106,.8)'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(W / 2, H / 2, sphR, 0, 7); ctx.stroke();
      ctx.strokeStyle = 'rgba(150,162,190,.3)';
      for (const lat of [-60, -30, 0, 30, 60]) {
        const pts = [];
        for (let a = 0; a <= 360; a += 6) {
          const rr = Rs * Math.cos(lat * Math.PI / 180), yy = Rs * Math.sin(lat * Math.PI / 180);
          pts.push(P(rr * Math.cos(a * Math.PI / 180), yy, rr * Math.sin(a * Math.PI / 180)));
        }
        ctx.beginPath(); pts.forEach((p, i) => { if (p[2] < -1) { ctx.moveTo(p[0], p[1]); } else { i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]); } });
        ctx.stroke();
      }
      // Rayos de visión hacia los residuos cercanos
      const o = P(0, 0, 0);
      ctx.strokeStyle = 'rgba(216,180,106,.42)'; ctx.beginPath();
      for (let i = 0; i < N; i++) if (dist[i] < 11) { ctx.moveTo(o[0], o[1]); ctx.lineTo(scr[i][0], scr[i][1]); }
      ctx.stroke();
      ctx.restore();
    }

    // Discos
    if (s3 > .01) {
      ctx.save(); ctx.globalAlpha = s3;
      drawDiscBase(ctx, L.n.x, L.n.y, L.n.R, { hemi: 'n' });
      drawDiscBase(ctx, L.s.x, L.s.y, L.s.R, { hemi: 's' });
      ctx.fillStyle = '#ece4d2'; ctx.textAlign = 'center';
      ctx.font = 'italic 500 16px "Bodoni Moda", serif';
      const small = L.n.R < 110;
      ctx.fillText(t('north'), L.n.x, L.n.y + L.n.R * 1.1 + (small ? 16 : 28));
      ctx.fillText(t('south'), L.s.x, L.s.y + L.s.R * 1.1 + (small ? 16 : 28));
      ctx.restore();
    }

    // Familia
    if (s4 > .01 && dust) {
      ctx.save(); ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.globalAlpha = s4; ctx.drawImage(dust, 0, 0); ctx.restore();
    }

    // Residuos de la proteína de referencia, del fondo al frente
    const ord = Array.from({ length: N }, (_, i) => i);
    if (s3 < .5) ord.sort((a, b) => scr[a][2] - scr[b][2]);
    const fade = 1 - s4 * .8;
    for (const i of ord) {
      const m = magOf(dist[i]);
      const shrink = (1 - s2) * 1 + s2 * 1;
      const base = s2 > .02 ? MAG_RADIUS[m] * (1 + s2 * .4) : 1.8 * Math.max(.6, scr[i][3]);
      const r = (s3 > .5 ? MAG_RADIUS[m] : base) * shrink;
      const sp = sprite(CLASS_COLOR[CLASSES[cls[i]]], Math.round(r * 4) / 4, dpr);
      const depth = s3 < .5 ? .45 + .55 * clamp01((scr[i][2] + 40) / 80) : 1;
      ctx.globalAlpha = fade * (depth + (1 - depth) * s3);
      ctx.drawImage(sp, scr[i][0] - sp.half, scr[i][1] - sp.half, sp.css, sp.css);
    }
    ctx.globalAlpha = 1;

    // Ligando real, en su sitio
    const ligA = clamp01(s * 1.6 - .25) * (1 - s3);
    if (ligA > .01 && lig.length) {
      const q = lig.map(a => P(a.x, a.y, a.z));
      ctx.globalAlpha = ligA; ctx.strokeStyle = '#d8b46a'; ctx.lineWidth = 2;
      ctx.beginPath();
      for (const [a, b] of bonds) { ctx.moveTo(q[a][0], q[a][1]); ctx.lineTo(q[b][0], q[b][1]); }
      ctx.stroke();
      // Punto de vista: carbono anomérico en el origen
      const o = P(0, 0, 0);
      ctx.fillStyle = '#fff8e8'; ctx.beginPath(); ctx.arc(o[0], o[1], 3.5, 0, 7); ctx.fill();
      ctx.strokeStyle = '#fff8e8'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(o[0], o[1], 8 + 3 * s2, 0, 7); ctx.stroke();
      ctx.globalAlpha = 1;
    }
  }

  function setTag() {
    const i = Math.max(0, Math.min(4, Math.round(s)));
    const txt = t('tag_' + (i + 1), { n: fmtInt(N), p: D.pdbs.length });
    if (tag.textContent !== txt) tag.textContent = txt;
    steps.forEach((el, k) => el.classList.toggle('is-on', k === i));
  }

  function progress() {
    // Cada paso tiene una meseta: el estado se queda quieto mientras se lee
    const first = steps[0].getBoundingClientRect();
    const last = steps[steps.length - 1].getBoundingClientRect();
    const mid = innerHeight * .55;
    const span = (last.top + last.height / 2) - (first.top + first.height / 2);
    const p = (mid - (first.top + first.height / 2)) / span * 4;
    const b = Math.floor(p), f = p - b;
    const plateau = clamp01((f - .3) / .4);
    return Math.max(0, Math.min(4, b + plateau));
  }

  function tick(ts) {
    raf = 0;
    target = progress();
    const dt = Math.min(.05, (ts - lastTs) / 1000 || 0); lastTs = ts;
    if (REDUCED) s = Math.round(target);
    else s += (target - s) * Math.min(1, dt * 7);
    if (!REDUCED && s < 2.5) spin += dt * .09;
    draw(); setTag();
    if (visible && (!REDUCED)) raf = requestAnimationFrame(tick);
  }
  function kick() { if (!raf) raf = requestAnimationFrame(tick); }

  new IntersectionObserver(es => { visible = es[0].isIntersecting; if (visible) { lastTs = performance.now(); kick(); } }, { rootMargin: '10% 0px' }).observe(stage.parentElement);
  addEventListener('scroll', () => { if (visible && REDUCED) kick(); }, { passive: true });
  let rt = 0;
  new ResizeObserver(() => { clearTimeout(rt); rt = setTimeout(() => { size(); kick(); }, 80); }).observe(stage);

  function copy() {
    document.querySelectorAll('[data-t="s5_p"]').forEach(el => { el.innerHTML = t('s5_p', { p: D.pdbs.length }); });
    setTag(); kick();
  }
  size();
  copy();
  onLang(copy);
  document.fonts?.ready.then(kick);
}
