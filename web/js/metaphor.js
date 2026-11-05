// Dos cielos: la misma figura desde la Tierra y desde el azúcar.
// A la izquierda, un esquema con estrellas reales de Orión y su entorno.
// A la derecha, datos reales de 3SCU: cada punto es un residuo del hemisferio norte.
import { D, t, onLang, fitCanvas, CLASSES, CLASS_COLOR, MAG_RADIUS, magOf, engravedStar, spanningTree, titleRes } from './core.js';

const $ = id => document.getElementById(id);

// Ascensión recta (horas), declinación (grados) y magnitud aparente, valores redondeados
const STARS = [
  ['Betelgeuse', 5.92, 7.4, .45], ['Rigel', 5.24, -8.2, .13], ['Bellatrix', 5.42, 6.35, 1.64], ['Mintaka', 5.53, -.3, 2.2],
  ['Alnilam', 5.6, -1.2, 1.7], ['Alnitak', 5.68, -1.9, 1.8], ['Saiph', 5.8, -9.7, 2.1],
  ['Sirius', 6.75, -16.7, -1.46], ['Aldebaran', 4.6, 16.5, .85], ['Procyon', 7.65, 5.2, .34],
];
const ORION = [[0, 2], [0, 5], [2, 3], [3, 4], [4, 5], [3, 1], [5, 6], [1, 6]];

function domeGeom(w, h) {
  const R = Math.min(w * .45, h * .78);
  return { cx: w / 2, base: h - Math.max(34, h * .13), R };
}
function drawDome(ctx, g) {
  const { cx, base, R } = g;
  ctx.save();
  ctx.lineWidth = 1;
  ctx.strokeStyle = 'rgba(150,162,190,.22)';
  for (const k of [.36, .72]) { ctx.beginPath(); ctx.ellipse(cx, base, R * k, R, 0, Math.PI, Math.PI * 2); ctx.stroke(); }
  for (const lat of [30, 60]) {
    const a = lat * Math.PI / 180, r = R * Math.cos(a);
    ctx.beginPath(); ctx.ellipse(cx, base - R * Math.sin(a), r, r * .14, 0, 0, Math.PI * 2); ctx.stroke();
  }
  ctx.strokeStyle = 'rgba(216,180,106,.8)';
  ctx.beginPath(); ctx.arc(cx, base, R, Math.PI, Math.PI * 2); ctx.stroke();
  ctx.beginPath(); ctx.ellipse(cx, base, R, R * .14, 0, 0, Math.PI * 2); ctx.stroke();
  ctx.restore();
}
function rays(ctx, x0, y0, pts) {
  ctx.save();
  ctx.strokeStyle = 'rgba(216,180,106,.38)'; ctx.lineWidth = 1;
  ctx.beginPath();
  for (const [x, y] of pts) { ctx.moveTo(x0, y0); ctx.lineTo(x, y); }
  ctx.stroke(); ctx.restore();
}
function label(ctx, txt, x, y, size, w) {
  ctx.font = `italic 500 ${size}px "Bodoni Moda", serif`;
  ctx.fillStyle = '#ece4d2'; ctx.textBaseline = 'middle';
  const flip = x + 9 + ctx.measureText(txt).width > w - 4;
  ctx.textAlign = flip ? 'right' : 'left';
  ctx.fillText(txt, x + (flip ? -9 : 9), y);
}

export function initMetaphor() {
  const R0 = D.struct(D.fam.reference_pdb);   // la referencia: 3SCU, fuera del conjunto actual (D-012)
  const rows = d3.range(R0.n).filter(i => R0.y[i] > 0);
  const lig = (D.extras.ligand || []).filter(a => a.element !== 'H');
  const bonds = [];
  for (let i = 0; i < lig.length; i++) for (let j = i + 1; j < lig.length; j++) {
    if (Math.hypot(lig[i].x - lig[j].x, lig[i].y - lig[j].y, lig[i].z - lig[j].z) < 1.75) bonds.push([i, j]);
  }

  function earth() {
    const { ctx, w, h } = fitCanvas($('sky-earth'));
    ctx.clearRect(0, 0, w, h);
    if (w < 60) return;
    const g = domeGeom(w, h);
    drawDome(ctx, g);
    const P = s => [g.cx - (s[1] - 5.9) * 15 / 44 * g.R * .8, g.base - g.R * .5 - (s[2] - 4) / 44 * g.R * .86];
    const q = STARS.map(P);
    rays(ctx, g.cx, g.base, [q[0], q[1], q[7]]);
    ctx.strokeStyle = 'rgba(216,180,106,.7)'; ctx.lineWidth = 1;
    ctx.beginPath();
    for (const [a, b] of ORION) { ctx.moveTo(q[a][0], q[a][1]); ctx.lineTo(q[b][0], q[b][1]); }
    ctx.stroke();
    STARS.forEach((s, i) => {
      const r = Math.max(4, 11 - s[3] * 2.6) * Math.min(1.2, Math.max(.75, g.R / 260));
      engravedStar(ctx, q[i][0], q[i][1], r, '#ece4d2');
    });
    const fs = Math.max(12, Math.min(16, g.R / 18));
    for (const i of [0, 1, 7, 8, 9]) label(ctx, t('star_' + STARS[i][0]), q[i][0], q[i][1] - 8, fs, w);
    // El suelo: una curva de la Tierra y una persona que mira
    ctx.strokeStyle = '#ece4d2'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(g.cx, g.base - 15, 3.2, 0, 7); ctx.moveTo(g.cx, g.base - 11.6); ctx.lineTo(g.cx, g.base); ctx.stroke();
    ctx.fillStyle = '#d8b46a';
    ctx.beginPath(); ctx.arc(g.cx, g.base, 2.4, 0, 7); ctx.fill();
  }

  function enzyme() {
    const { ctx, w, h } = fitCanvas($('sky-enzyme'));
    ctx.clearRect(0, 0, w, h);
    if (w < 60) return;
    const g = domeGeom(w, h);
    drawDome(ctx, g);
    const P = i => { const d = R0.distance[i]; return [g.cx + R0.x[i] / d * g.R, g.base - R0.y[i] / d * g.R]; };
    const key = rows.filter(i => R0.pos[i] > 0 && R0.distance[i] < 8.5);
    rays(ctx, g.cx, g.base, key.map(P));
    const pts = key.map(i => { const [u, v] = P(i); return { u, v }; });
    ctx.strokeStyle = 'rgba(216,180,106,.7)'; ctx.lineWidth = 1;
    ctx.beginPath();
    for (const [a, b] of spanningTree(pts)) { ctx.moveTo(pts[a].u, pts[a].v); ctx.lineTo(pts[b].u, pts[b].v); }
    ctx.stroke();
    const k = Math.min(1.25, Math.max(.8, g.R / 240));
    ctx.globalCompositeOperation = 'lighter';
    for (const i of rows) {
      const m = magOf(R0.distance[i]);
      if (m < 2) continue;
      const [x, y] = P(i);
      ctx.fillStyle = CLASS_COLOR[CLASSES[R0.cls[i]]];
      ctx.beginPath(); ctx.arc(x, y, MAG_RADIUS[m] * k, 0, 7); ctx.fill();
    }
    ctx.globalCompositeOperation = 'source-over';
    for (const i of rows) {
      const m = magOf(R0.distance[i]);
      if (m >= 2) continue;
      const [x, y] = P(i);
      engravedStar(ctx, x, y, (m ? 7.5 : 10.5) * k, CLASS_COLOR[CLASSES[R0.cls[i]]]);
    }
    const fs = Math.max(12, Math.min(16, g.R / 18));
    for (const p of [D.refPositions[0], 333, 121]) {
      const i = rows.find(r => R0.pos[r] === p);
      if (i === undefined) continue;
      const [x, y] = P(i);
      label(ctx, `${titleRes(R0.residue[i])} ${p}`, x, y - 8, fs, w);
    }
    // El suelo: la celopentaosa real, a escala propia, con el punto de vista en el origen
    if (lig.length) {
      const s = Math.max(1.5, g.R / 110);
      ctx.strokeStyle = '#d8b46a'; ctx.lineWidth = 1.4;
      ctx.beginPath();
      for (const [a, b] of bonds) { ctx.moveTo(g.cx + lig[a].x * s, g.base - lig[a].y * s); ctx.lineTo(g.cx + lig[b].x * s, g.base - lig[b].y * s); }
      ctx.stroke();
    }
    ctx.fillStyle = '#fff8e8'; ctx.strokeStyle = '#fff8e8'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(g.cx, g.base, 3, 0, 7); ctx.fill();
    ctx.beginPath(); ctx.arc(g.cx, g.base, 7, 0, 7); ctx.stroke();
  }

  function table() {
    const n = 7;
    const body = Array.from({ length: n }, (_, i) => `<tr><td>${t(`r${i + 1}a`)}</td><td>${t(`r${i + 1}b`, { p: D.pdbs.length })}</td></tr>`).join('');
    $('rosetta').innerHTML = `<thead><tr><th scope="col">${t('r_head_a')}</th><th scope="col">${t('r_head_b')}</th></tr></thead><tbody>${body}</tbody>`;
    $('sky-enzyme-cap').textContent = t('tw2_p', { n: rows.length });
  }

  const draw = () => { earth(); enzyme(); };
  table(); draw();
  onLang(() => { table(); draw(); });
  document.fonts?.ready.then(draw);
  let rt = 0;
  new ResizeObserver(() => { clearTimeout(rt); rt = setTimeout(draw, 80); }).observe($('metaphor'));
}
