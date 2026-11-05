// Núcleo común: idioma, datos, proyección y dibujo base del disco.
import { TEXTS } from './i18n.js';

export const REDUCED = matchMedia('(prefers-reduced-motion: reduce)').matches;

/* Idioma */
const langSubs = [];
function readLang() {
  const q = new URLSearchParams(location.search).get('lang');
  if (q === 'es' || q === 'en') {
    try { localStorage.setItem('easimap-lang', q); } catch (e) { /* sin persistencia */ }
    return q;
  }
  try { const s = localStorage.getItem('easimap-lang'); if (s === 'es' || s === 'en') return s; } catch (e) { /* almacenamiento bloqueado */ }
  return 'en';
}
export let lang = readLang();
export function t(key, vars) {
  let s = TEXTS[lang][key];
  if (s === undefined) s = TEXTS.en[key];
  if (s === undefined) return key;
  if (vars) for (const k in vars) s = s.split('{' + k + '}').join(vars[k]);
  return s;
}
export const onLang = f => { langSubs.push(f); };
export function setLang(l) {
  lang = l;
  try { localStorage.setItem('easimap-lang', l); } catch (e) { /* sin persistencia */ }
  applyStatic();
  langSubs.forEach(f => f());
}
export function applyStatic() {
  document.documentElement.lang = lang;
  document.title = t('doc_title');
  document.querySelectorAll('[data-t]').forEach(el => { el.innerHTML = t(el.dataset.t); });
  document.querySelectorAll('[data-t-aria]').forEach(el => el.setAttribute('aria-label', t(el.dataset.tAria)));
  document.querySelectorAll('[data-t-ph]').forEach(el => el.setAttribute('placeholder', t(el.dataset.tPh)));
  document.querySelectorAll('[data-lang]').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.lang === lang)));
}
export const fmt = (v, d = 1) => (v == null || Number.isNaN(v)) ? '–' : new Intl.NumberFormat(lang, { minimumFractionDigits: d, maximumFractionDigits: d }).format(v);
export const fmtInt = v => new Intl.NumberFormat(lang).format(v);

/* Clases de residuo y colores. Orden fijo; validado con validate_palette.js en oscuro. */
export const CLASSES = ['acid', 'base', 'polar', 'arom', 'gp', 'hyd'];
export const CLASS_COLOR = { acid: '#f2765a', base: '#6bb0f2', polar: '#5fd3b0', arom: '#f2c84b', gp: '#d886f0', hyd: '#aeb6c4' };
export const CLASS_RES = {
  acid: ['ASP', 'GLU'], base: ['ARG', 'HIS', 'LYS'], polar: ['ASN', 'CYS', 'GLN', 'SER', 'THR'],
  arom: ['PHE', 'TRP', 'TYR'], gp: ['GLY', 'PRO'], hyd: ['ALA', 'ILE', 'LEU', 'MET', 'VAL'],
};
export const ALIAS = { MSE: 'MET', CAS: 'CYS' };
const RES_CLASS = {};
CLASSES.forEach((c, i) => CLASS_RES[c].forEach(r => { RES_CLASS[r] = i; }));
export const classOf = r => RES_CLASS[ALIAS[r] || r] ?? 5;
export const DOMAINS = ['Eukaryota', 'Bacteria', 'Archaea', 'Unknown'];
export const DOMAIN_COLOR = { Eukaryota: '#f2c84b', Bacteria: '#6bb0f2', Archaea: '#f2765a', Unknown: '#7f8796' };
export const titleRes = r => r[0] + r.slice(1).toLowerCase();

/* Magnitud: la distancia se convierte en tamaño, como en un atlas estelar */
export const MAG_LIMITS = [6, 9, 13, 18, 25];
export const MAG_RADIUS = [2.9, 2.2, 1.65, 1.25, .95, .72];
export function magOf(d) { for (let i = 0; i < 5; i++) if (d < MAG_LIMITS[i]) return i; return 5; }

/* Proyección. El radio sale de d3; la orientación sigue el contrato:
   vista desde el centro, norte antihorario, sur horario, φ = 0 arriba. */
const HALF = Math.PI / 2;
export function radialFn(proj) {
  const raw = proj === 'eqd' ? d3.geoAzimuthalEquidistantRaw : d3.geoAzimuthalEqualAreaRaw;
  const rim = raw(HALF, 0)[0];
  return tt => raw(tt, 0)[0] / rim;
}
export function projectInto(fit, proj, ux, uy, hemi) {
  const rad = radialFn(proj);
  const { theta, phi } = fit;
  for (let i = 0; i < theta.length; i++) {
    const north = theta[i] < HALF;
    const r = rad(north ? theta[i] : Math.PI - theta[i]);
    ux[i] = (north ? -1 : 1) * r * Math.sin(phi[i]);
    uy[i] = -r * Math.cos(phi[i]);
    hemi[i] = north ? 0 : 1;
  }
}
export function dirToDisc(x, y, z, rad) {
  const d = Math.hypot(x, y, z) || 1;
  const th = Math.acos(Math.max(-1, Math.min(1, y / d)));
  const ph = Math.atan2(z, x);
  const north = th < HALF;
  const r = rad(north ? th : Math.PI - th);
  return { u: (north ? -1 : 1) * r * Math.sin(ph), v: -r * Math.cos(ph), north, d, th, ph };
}

/* Datos: solo los conjuntos actuales del pipeline (D-023). */
export const FAMILIES = ['GH1', 'GH5', 'GH16', 'GH18', 'GH20'];
export const PFAM = { GH1: 'PF00232', GH5: 'PF00150', GH16: 'PF00722', GH18: 'PF00704', GH20: 'PF00728' };
const famHash = (location.hash.match(/[#?&]fam=([a-z0-9]+)/i) || [])[1]?.toUpperCase();
export const FAM = FAMILIES.includes(famHash) ? famHash : 'GH1';
export const dataUrl = f => `data/${f}${window.EASIMAP_V ? '?v=' + window.EASIMAP_V : ''}`;
export const D = { ready: false };
export async function loadData(family = 'GH1') {
  const stem = family.toLowerCase();
  const [fam, extras] = await Promise.all([
    fetch(dataUrl(`${stem}.json`)).then(r => { if (!r.ok) throw new Error(stem + '.json ' + r.status); return r.json(); }),
    fetch(dataUrl(`${stem}.extras.json`)).then(r => (r.ok ? r.json() : {})).catch(() => ({})),
  ]);
  D.fam = fam; D.extras = extras; D.family = family;
  D.source = fam.source; D.generated = fam.generated;
  D.refPositions = fam.reference_positions;
  const c = fam.fits.catalytic;
  const n = c.pdb.length;
  D.n = n;
  const pdbs = []; const pdbIdx = new Map();
  D.pdbI = new Uint16Array(n); D.cls = new Uint8Array(n); D.dom = new Uint8Array(n); D.pos = new Int16Array(n);
  for (let i = 0; i < n; i++) {
    let k = pdbIdx.get(c.pdb[i]);
    if (k === undefined) { k = pdbs.length; pdbs.push(c.pdb[i]); pdbIdx.set(c.pdb[i], k); }
    D.pdbI[i] = k;
    D.cls[i] = classOf(c.residue[i]);
    D.dom[i] = Math.max(0, DOMAINS.indexOf(c.domain[i]));
    D.pos[i] = c.hmm_position[i] ?? 0;
  }
  D.pdbs = pdbs; D.pdbIdx = pdbIdx;
  D.residue = c.residue; D.auth = c.auth_seq_id; D.chain = c.chain; D.pdb = c.pdb;
  D.fits = {};
  for (const key of ['catalytic', 'core']) {
    const f = fam.fits[key];
    if (!f) continue;
    D.fits[key] = {
      theta: Float32Array.from(f.theta), phi: Float32Array.from(f.phi), distance: Float32Array.from(f.distance),
      x: Float32Array.from(f.x), y: Float32Array.from(f.y), z: Float32Array.from(f.z),
    };
  }
  D.positions = new Map(fam.positions.map(p => [p.hmm_position, p]));
  D.posRows = new Map();
  D.protRows = pdbs.map(() => []);
  for (let i = 0; i < n; i++) {
    D.protRows[D.pdbI[i]].push(i);
    const p = D.pos[i];
    if (p > 0) { let a = D.posRows.get(p); if (!a) D.posRows.set(p, a = []); a.push(i); }
  }
  D.proteinInfo = extras.proteins || {};
  D.organism = code => (D.proteinInfo[code] && D.proteinInfo[code].organism) || '';
  // Una proteína como columnas propias. Si no está en este conjunto, sale de extras (3SCU).
  const cache = new Map();
  D.struct = code => {
    if (cache.has(code)) return cache.get(code);
    const s = extras.reference_structures?.[code];
    let o;
    if (s) o = { ...s, pos: s.hmm_position.map(v => v ?? 0), auth: s.auth_seq_id };
    else {
      const f = D.fits.catalytic, rows = D.protRows[D.pdbIdx.get(code)];
      const col = a => rows.map(i => a[i]);
      o = { x: col(f.x), y: col(f.y), z: col(f.z), distance: col(f.distance), theta: col(f.theta), phi: col(f.phi), residue: col(D.residue), pos: col(D.pos), auth: col(D.auth) };
    }
    o.n = o.x.length;
    o.cls = o.residue.map(classOf);
    o.ux = new Float32Array(o.n); o.uy = new Float32Array(o.n); o.hemi = new Uint8Array(o.n);
    projectInto(o, 'equal', o.ux, o.uy, o.hemi);
    cache.set(code, o);
    return o;
  };
  D.ready = true;
  return D;
}

export async function loadFamilies() {
  D.families = await fetch(dataUrl('families.json')).then(r => (r.ok ? r.json() : [])).catch(() => []);
  return D.families;
}

/* Canvas con densidad de píxel */
export function fitCanvas(cv, maxDpr = 2) {
  const r = cv.getBoundingClientRect();
  const dpr = Math.min(window.devicePixelRatio || 1, maxDpr);
  const w = Math.max(1, Math.round(r.width * dpr)), h = Math.max(1, Math.round(r.height * dpr));
  if (cv.width !== w || cv.height !== h) { cv.width = w; cv.height = h; }
  const ctx = cv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, w: r.width, h: r.height, dpr };
}

/* Sprites: un punto por color y magnitud. Más rápido que arc y deja acumular luz. */
const spriteCache = new Map();
export function sprite(color, radius, dpr) {
  const key = color + radius + dpr;
  let s = spriteCache.get(key);
  if (s) return s;
  const size = Math.ceil((radius + 1.2) * 2 * dpr);
  s = document.createElement('canvas');
  s.width = s.height = size;
  const c = s.getContext('2d');
  c.fillStyle = color;
  c.beginPath(); c.arc(size / 2, size / 2, radius * dpr, 0, 7); c.fill();
  s.half = size / 2 / dpr; s.css = size / dpr;
  spriteCache.set(key, s);
  return s;
}

/* Disco de atlas: anillos de θ, radios de φ y limbo graduado */
export function drawDiscBase(ctx, cx, cy, R, { proj = 'equal', hemi = 'n', labels = true, rim = true, strength = 1 } = {}) {
  if (!(R > 0)) return;   // lienzo aún sin medidas, por ejemplo en una pestaña en segundo plano
  const rad = radialFn(proj);
  ctx.save();
  ctx.lineWidth = 1;
  ctx.strokeStyle = `rgba(150,162,190,${.20 * strength})`;
  for (const a of [15, 30, 45, 60, 75]) {
    const r = rad(a * Math.PI / 180) * R;
    ctx.setLineDash(a % 30 ? [1, 4] : []);
    ctx.beginPath(); ctx.arc(cx, cy, r, 0, 7); ctx.stroke();
  }
  ctx.setLineDash([]);
  const s = hemi === 'n' ? -1 : 1;
  for (let a = 0; a < 360; a += 30) {
    const rr = a % 90 ? rad(15 * Math.PI / 180) * R : 0;
    const sn = Math.sin(a * Math.PI / 180) * s, cs = -Math.cos(a * Math.PI / 180);
    ctx.beginPath(); ctx.moveTo(cx + sn * rr, cy + cs * rr); ctx.lineTo(cx + sn * R, cy + cs * R); ctx.stroke();
  }
  if (rim) {
    ctx.strokeStyle = `rgba(216,180,106,${.75 * strength})`;
    ctx.beginPath(); ctx.arc(cx, cy, R, 0, 7); ctx.stroke();
    const R2 = R * 1.035 + 2;
    ctx.strokeStyle = `rgba(216,180,106,${.42 * strength})`;
    ctx.beginPath(); ctx.arc(cx, cy, R2, 0, 7); ctx.stroke();
    ctx.beginPath();
    const step = R > 170 ? 2 : 5;
    for (let a = 0; a < 360; a += step) {
      const sn = Math.sin(a * Math.PI / 180), cs = -Math.cos(a * Math.PI / 180);
      const len = a % 30 === 0 ? 1 : (a % 10 === 0 ? .62 : .32);
      const r0 = R + (R2 - R) * (1 - len);
      ctx.moveTo(cx + sn * r0, cy + cs * r0); ctx.lineTo(cx + sn * R2, cy + cs * R2);
    }
    ctx.stroke();
    if (labels && R > 110) {
      ctx.fillStyle = `rgba(163,169,180,${.95 * strength})`;
      ctx.font = `400 ${R > 260 ? 11 : 9.5}px "Spline Sans Mono", monospace`;
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      const Rl = R2 + (R > 260 ? 15 : 11);
      for (let a = 0; a < 360; a += 30) {
        const sn = Math.sin(a * Math.PI / 180) * s, cs = -Math.cos(a * Math.PI / 180);
        ctx.fillText(a + '°', cx + sn * Rl, cy + cs * Rl);
      }
    }
  }
  ctx.restore();
}

/* Estrella grabada: cuatro puntas finas para las de primera magnitud */
export function engravedStar(ctx, x, y, r, color, alpha = 1) {
  ctx.globalAlpha = alpha;
  ctx.fillStyle = color;
  const k = r * .17;
  ctx.beginPath();
  ctx.moveTo(x, y - r); ctx.lineTo(x + k, y - k); ctx.lineTo(x + r, y); ctx.lineTo(x + k, y + k);
  ctx.lineTo(x, y + r); ctx.lineTo(x - k, y + k); ctx.lineTo(x - r, y); ctx.lineTo(x - k, y - k);
  ctx.closePath(); ctx.fill();
  ctx.globalAlpha = 1;
}

/* Árbol de expansión mínima sobre puntos: las constelaciones unen residuos reales */
export function spanningTree(pts) {
  const n = pts.length; if (n < 2) return [];
  const inT = new Array(n).fill(false); inT[0] = true;
  const edges = [];
  for (let e = 0; e < n - 1; e++) {
    let best = null, bd = Infinity;
    for (let i = 0; i < n; i++) if (inT[i]) for (let j = 0; j < n; j++) if (!inT[j]) {
      const d = Math.hypot(pts[i].u - pts[j].u, pts[i].v - pts[j].v);
      if (d < bd) { bd = d; best = [i, j]; }
    }
    inT[best[1]] = true; edges.push(best);
  }
  return edges;
}

export const ease = t => t < .5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
export const clamp01 = v => Math.max(0, Math.min(1, v));
