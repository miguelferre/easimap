// Explorador: dos discos en Canvas, filtros al momento, espectro de posiciones,
// selección con vecinos y enlace con Mol*.
import {
  D, t, lang, onLang, REDUCED, fitCanvas, drawDiscBase, radialFn, projectInto, magOf, MAG_LIMITS, MAG_RADIUS,
  CLASSES, CLASS_COLOR, CLASS_RES, ALIAS, DOMAINS, DOMAIN_COLOR, FAMILIES, PFAM, sprite, ease, fmt, fmtInt, titleRes,
} from './core.js';
import { molShow, molIsOn, molMessage } from './mol.js';

const DMAX = 45;              // el extremo del deslizador significa «sin límite»
const NEIGH = 7;              // Å
const DEFAULTS = { fit: 'catalytic', proj: 'equal', color: 'class', dmin: 0, dmax: 14, pos: '', pdb: '' };
const $ = id => document.getElementById(id);

export function initExplorer() {
  const S = { ...DEFAULTS, res: new Set(Object.values(CLASS_RES).flat()), doms: new Set(DOMAINS), sel: -1 };
  const allRes = [...S.res];
  const hasCore = !!D.fits.core;

  /* Buffers */
  const n = D.n;
  const ux = new Float32Array(n), uy = new Float32Array(n), hemi = new Uint8Array(n);
  const fromX = new Float32Array(n), fromY = new Float32Array(n), toX = new Float32Array(n), toY = new Float32Array(n);
  const mask = new Uint8Array(n);
  const mag = new Uint8Array(n);
  let vis = [[], []], visCount = 0, proteinsShown = 0;
  let trees = [null, null], tween = 0, raf = 0, samePos = new Set(), neighbours = [];

  function reproject(animate) {
    const f = D.fits[S.fit];
    fromX.set(ux); fromY.set(uy);
    const oldH = hemi.slice();
    projectInto(f, S.proj, toX, toY, hemi);
    for (let i = 0; i < n; i++) mag[i] = magOf(f.distance[i]);
    if (animate && !REDUCED) {
      for (let i = 0; i < n; i++) if (oldH[i] !== hemi[i]) { fromX[i] = toX[i]; fromY[i] = toY[i]; }
      const t0 = performance.now();
      tween = t0;
      const step = ts => {
        if (tween !== t0) return;
        const e = ease(Math.min(1, (ts - t0) / 440));
        for (let i = 0; i < n; i++) { ux[i] = fromX[i] + (toX[i] - fromX[i]) * e; uy[i] = fromY[i] + (toY[i] - fromY[i]) * e; }
        drawNow();
        if (e < 1) requestAnimationFrame(step); else { tween = 0; trees = [null, null]; }
      };
      requestAnimationFrame(step);
    } else { ux.set(toX); uy.set(toY); trees = [null, null]; }
  }

  /* Filtro */
  function parsePos(txt) {
    const out = [];
    for (const part of txt.split(/[,;\s]+/)) {
      if (!part) continue;
      const m = part.match(/^(\d+)(?:[-–](\d+))?$/);
      if (m) out.push([+m[1], +(m[2] ?? m[1])]);
    }
    return out.length ? out : null;
  }
  function apply({ keepSel = false, silent = false } = {}) {
    const f = D.fits[S.fit];
    const ranges = parsePos(S.pos);
    const pdbI = S.pdb ? D.pdbIdx.get(S.pdb) : -1;
    const hi = S.dmax >= DMAX ? Infinity : S.dmax;
    const vN = [], vS = [];
    const seen = new Set();
    for (let i = 0; i < n; i++) {
      const d = f.distance[i];
      let ok = d >= S.dmin && d <= hi && S.doms.has(DOMAINS[D.dom[i]]);
      if (ok) { const r = D.residue[i]; ok = S.res.has(ALIAS[r] || r); }
      if (ok && pdbI >= 0) ok = D.pdbI[i] === pdbI;
      if (ok && ranges) { const p = D.pos[i]; ok = false; for (const [a, b] of ranges) if (p >= a && p <= b) { ok = true; break; } }
      mask[i] = ok ? 1 : 0;
      if (ok) { (hemi[i] ? vS : vN).push(i); seen.add(D.pdbI[i]); }
    }
    // Lo cercano se pinta encima
    const byFar = (a, b) => f.distance[b] - f.distance[a];
    vN.sort(byFar); vS.sort(byFar);
    vis = [vN, vS]; visCount = vN.length + vS.length; proteinsShown = seen.size;
    trees = [null, null];
    if (!keepSel && S.sel >= 0 && !mask[S.sel]) select(-1);
    updateCount(); drawHist(); drawSpectrum();
    if ($('tbl-box').open) renderTable();
    if (!silent) writeHash();
    draw();
  }
  function updateCount() {
    $('count-n').textContent = fmtInt(visCount);
    $('count-of').textContent = t('of_total', { n: fmtInt(n), p: proteinsShown });
    const hi = S.dmax >= DMAX ? DMAX + '+' : fmt(S.dmax, S.dmax % 1 ? 1 : 0);
    $('f-dist-out').textContent = `${fmt(S.dmin, S.dmin % 1 ? 1 : 0)}–${hi} Å`;
  }

  /* Discos */
  const discs = [...document.querySelectorAll('.disc')].map((el, h) => {
    const cv = el.querySelector('canvas');
    return { el, cv, h, hemi: h ? 's' : 'n', tf: d3.zoomIdentity, w: 0, R: 0, cx: 0, cy: 0, dpr: 1 };
  });
  const discsEl = $('discs');
  discsEl.dataset.show = 'n';

  function sizeDiscs() {
    for (const d of discs) {
      if (!d.cv.offsetParent) continue;
      const f = fitCanvas(d.cv);
      d.w = f.w; d.dpr = f.dpr;
      d.R = (f.w / 2 - (f.w > 540 ? 28 : 21)) / 1.035;
      d.cx = d.cy = f.w / 2;
      d.zoom.extent([[0, 0], [f.w, f.w]]).translateExtent([[0, 0], [f.w, f.w]]);
    }
    trees = [null, null];
    drawHist(); drawSpectrum(); draw();
  }

  const colorOf = i => S.color === 'class' ? CLASS_COLOR[CLASSES[D.cls[i]]] : DOMAIN_COLOR[DOMAINS[D.dom[i]]];
  const screen = (d, i) => [d.tf.applyX(d.cx + ux[i] * d.R), d.tf.applyY(d.cy + uy[i] * d.R)];

  function refMarks(d, ctx, sizeK) {
    ctx.save();
    ctx.strokeStyle = '#ece4d2'; ctx.fillStyle = '#ece4d2'; ctx.lineWidth = 1;
    ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
    ctx.font = `italic 500 ${Math.round(12 + 2 * sizeK)}px "Bodoni Moda", serif`;
    for (const p of D.refPositions) {
      const rows = D.posRows.get(p) || [];
      let sx = 0, sy = 0, k = 0;
      for (const i of rows) if (hemi[i] === d.h) { sx += ux[i]; sy += uy[i]; k++; }
      if (k < rows.length / 2) continue;
      const x = d.tf.applyX(d.cx + sx / k * d.R), y = d.tf.applyY(d.cy + sy / k * d.R);
      const r = 8;
      ctx.beginPath(); ctx.moveTo(x, y - r); ctx.lineTo(x + r, y); ctx.lineTo(x, y + r); ctx.lineTo(x - r, y); ctx.closePath(); ctx.stroke();
      ctx.lineWidth = 3; ctx.strokeStyle = '#080c16'; ctx.strokeText(String(p), x + r + 4, y);
      ctx.lineWidth = 1; ctx.strokeStyle = '#ece4d2'; ctx.fillText(String(p), x + r + 4, y);
    }
    ctx.restore();
  }

  function drawDisc(d) {
    if (!d.w) return;
    const ctx = d.cv.getContext('2d');
    ctx.setTransform(d.dpr, 0, 0, d.dpr, 0, 0);
    ctx.clearRect(0, 0, d.w, d.w);
    const k = d.tf.k;
    drawDiscBase(ctx, d.tf.applyX(d.cx), d.tf.applyY(d.cy), d.R * k, { proj: S.proj, hemi: d.hemi });
    const list = vis[d.h];
    const sizeK = Math.max(.6, Math.min(1, d.w / 620));
    const zoomK = Math.pow(k, .55);
    const dense = visCount > 9000 ? .42 : (visCount > 2500 ? .55 : (visCount > 500 ? .8 : 1));
    const dimmed = S.sel >= 0;
    ctx.globalCompositeOperation = 'lighter';
    const X = d.tf.k, TX = d.tf.x, TY = d.tf.y;
    for (let j = 0; j < list.length; j++) {
      const i = list[j];
      const col = colorOf(i);
      const r = Math.round(MAG_RADIUS[mag[i]] * sizeK * zoomK * 4) / 4;
      const sp = sprite(col, Math.max(.5, r), d.dpr);
      const x = (d.cx + ux[i] * d.R) * X + TX, y = (d.cy + uy[i] * d.R) * X + TY;
      ctx.globalAlpha = dimmed ? (samePos.has(i) ? 1 : .26) : dense;
      ctx.drawImage(sp, x - sp.half, y - sp.half, sp.css, sp.css);
    }
    ctx.globalAlpha = 1;
    ctx.globalCompositeOperation = 'source-over';
    refMarks(d, ctx, sizeK);

    // Selección
    if (S.sel >= 0) {
      ctx.lineWidth = 1;
      if (hemi[S.sel] === d.h) {
        const [sx, sy] = screen(d, S.sel);
        ctx.strokeStyle = 'rgba(216,180,106,.8)';
        ctx.beginPath();
        for (const m of neighbours) if (hemi[m.i] === d.h) { const [x, y] = screen(d, m.i); ctx.moveTo(sx, sy); ctx.lineTo(x, y); }
        ctx.stroke();
      }
      ctx.strokeStyle = '#d8b46a';
      for (const m of neighbours) {
        if (hemi[m.i] !== d.h) continue;
        const [x, y] = screen(d, m.i);
        ctx.fillStyle = colorOf(m.i);
        ctx.beginPath(); ctx.arc(x, y, 2, 0, 7); ctx.fill();
        ctx.beginPath(); ctx.arc(x, y, 4.5, 0, 7); ctx.stroke();
      }
      if (hemi[S.sel] === d.h) {
        const [x, y] = screen(d, S.sel);
        ctx.strokeStyle = '#fff8e8';
        ctx.beginPath(); ctx.arc(x, y, 7, 0, 7);
        for (const [a, b] of [[0, -1], [0, 1], [-1, 0], [1, 0]]) { ctx.moveTo(x + a * 10, y + b * 10); ctx.lineTo(x + a * 18, y + b * 18); }
        ctx.stroke();
        ctx.fillStyle = '#fff8e8'; ctx.font = 'italic 500 15px "Bodoni Moda", serif'; ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
        const lab = `${titleRes(D.residue[S.sel])} ${D.pos[S.sel] || ''}`;
        const lx = x + 22 + ctx.measureText(lab).width > d.w ? x - 22 - ctx.measureText(lab).width : x + 22;
        ctx.fillText(lab, lx, y - 12);
      }
    }
    if (hoverRow >= 0 && hemi[hoverRow] === d.h) {
      const [x, y] = screen(d, hoverRow);
      ctx.strokeStyle = '#ece4d2'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(x, y, 6, 0, 7); ctx.stroke();
    }
  }
  function drawNow() { discs.forEach(drawDisc); }
  function draw() { if (!raf) raf = requestAnimationFrame(() => { raf = 0; drawNow(); }); }

  /* Zoom y gestos */
  for (const d of discs) {
    d.zoom = d3.zoom().scaleExtent([1, 20]).clickDistance(5)
      .filter(e => {
        if (e.type === 'wheel') return e.ctrlKey || e.metaKey;
        if (e.type === 'touchstart') return e.touches.length >= 2 || d.tf.k > 1.01;
        if (e.type === 'mousedown') return e.button === 0 && d.tf.k > 1.01;
        return !e.button;
      })
      .on('zoom', e => {
        d.tf = e.transform; trees = [null, null];
        d.cv.style.touchAction = d.tf.k > 1.01 ? 'none' : 'pan-y';
        d.cv.style.cursor = d.tf.k > 1.01 ? 'grab' : 'crosshair';
        drawDisc(d);
      });
    d3.select(d.cv).call(d.zoom).on('dblclick.zoom', null);
    d.cv.style.touchAction = 'pan-y';
  }
  const activeDiscs = () => discs.filter(d => d.cv.offsetParent);
  $('z-in').addEventListener('click', () => activeDiscs().forEach(d => d3.select(d.cv).transition().duration(REDUCED ? 0 : 200).call(d.zoom.scaleBy, 1.8)));
  $('z-out').addEventListener('click', () => activeDiscs().forEach(d => d3.select(d.cv).transition().duration(REDUCED ? 0 : 200).call(d.zoom.scaleBy, 1 / 1.8)));
  $('z-0').addEventListener('click', () => activeDiscs().forEach(d => d3.select(d.cv).transition().duration(REDUCED ? 0 : 200).call(d.zoom.transform, d3.zoomIdentity)));

  /* Punto bajo el cursor */
  let hoverRow = -1;
  function tree(d) {
    if (!trees[d.h]) trees[d.h] = d3.quadtree().x(i => d.cx + ux[i] * d.R).y(i => d.cy + uy[i] * d.R).addAll(vis[d.h]);
    return trees[d.h];
  }
  function pick(d, e, radius) {
    const r = d.cv.getBoundingClientRect();
    const [x, y] = d.tf.invert([e.clientX - r.left, e.clientY - r.top]);
    const i = tree(d).find(x, y, radius / d.tf.k);
    return i === undefined ? -1 : i;
  }
  const tip = $('tip');
  function showTip(d, i) {
    if (i < 0) { tip.hidden = true; return; }
    const f = D.fits[S.fit];
    const code = D.pdb[i], p = D.pos[i];
    tip.innerHTML = `<b>${titleRes(D.residue[i])} ${D.auth[i] ?? ''}</b> · ${code}<span>${p ? t('tip_pos', { p }) : t('tip_nopos')} · ${fmt(f.distance[i])} Å</span><span>${D.organism(code)}</span>`;
    tip.hidden = false;
    const [x, y] = screen(d, i);
    const host = discsEl.getBoundingClientRect(), r = d.cv.getBoundingClientRect();
    const tw = tip.offsetWidth, th = tip.offsetHeight;
    let left = r.left - host.left + x + 14, top = r.top - host.top + y + 14;
    if (left + tw > host.width) left = r.left - host.left + x - tw - 14;
    if (top + th > host.height) top = r.top - host.top + y - th - 14;
    tip.style.left = Math.max(0, left) + 'px'; tip.style.top = Math.max(0, top) + 'px';
  }
  for (const d of discs) {
    d.cv.addEventListener('pointermove', e => {
      if (e.pointerType === 'touch' || e.buttons) return;
      const i = pick(d, e, 12);
      if (i !== hoverRow) { hoverRow = i; showTip(d, i); draw(); }
    });
    d.cv.addEventListener('pointerleave', () => { hoverRow = -1; tip.hidden = true; draw(); });
    d.cv.addEventListener('click', e => {
      const touch = e.pointerType === 'touch' || matchMedia('(pointer: coarse)').matches;
      const i = pick(d, e, touch ? 24 : 12);
      select(i);
      if (touch) showTip(d, i);
    });
  }

  /* Selección */
  function select(i) {
    S.sel = i;
    samePos = new Set(); neighbours = [];
    const card = $('detail-card');
    if (i < 0) { card.hidden = true; drawSpectrum(); draw(); return; }
    const f = D.fits[S.fit];
    if (D.pos[i] > 0) (D.posRows.get(D.pos[i]) || []).forEach(r => samePos.add(r));
    samePos.add(i);
    for (const r of D.protRows[D.pdbI[i]]) {
      if (r === i) continue;
      const dd = Math.hypot(f.x[r] - f.x[i], f.y[r] - f.y[i], f.z[r] - f.z[i]);
      if (dd <= NEIGH) { neighbours.push({ i: r, d: dd }); samePos.add(r); }
    }
    neighbours.sort((a, b) => a.d - b.d);
    renderCard();
    card.hidden = false;
    if (molIsOn()) sendToMol();
    drawSpectrum(); draw();
  }
  function renderCard() {
    const i = S.sel; if (i < 0) return;
    const f = D.fits[S.fit];
    const code = D.pdb[i], p = D.pos[i], info = D.positions.get(p);
    const org = D.organism(code);
    $('d-eyebrow').textContent = `${code}${org ? ' · ' + org : ''}`;
    $('d-name').textContent = D.proteinInfo[code]?.name || '';
    $('d-title').textContent = `${titleRes(D.residue[i])} ${p || ''}`.trim();
    const deg = v => fmt(v * 180 / Math.PI) + '°';
    const rows = [
      [t('d_hmm'), p || t('d_none')],
      [t('d_pdbres'), `${D.residue[i]} ${D.auth[i] ?? '–'}`],
      [t('d_dist'), fmt(f.distance[i]) + ' Å'],
      [t('d_theta'), deg(f.theta[i])],
      [t('d_phi'), deg(f.phi[i])],
      [t('d_rmsf'), info ? fmt(info.rmsf, 2) + ' Å' : '–'],
      [t('d_n'), info ? info.n : '–'],
      [t('d_domain'), t('dom_' + DOMAINS[D.dom[i]])],
    ];
    $('d-facts').innerHTML = rows.map(([a, b]) => `<dt>${a}</dt><dd class="num">${b}</dd>`).join('');
    $('d-links').innerHTML = `${t('d_links')} <a href="https://www.rcsb.org/structure/${code}" rel="noopener">RCSB</a> · <a href="https://www.ebi.ac.uk/pdbe/entry/pdb/${code.toLowerCase()}" rel="noopener">PDBe</a> · <a href="https://www.ebi.ac.uk/interpro/entry/pfam/${PFAM[D.family]}/" rel="noopener">Pfam</a>`;
    const box = $('d-neigh');
    const shown = neighbours.slice(0, 10);
    box.innerHTML = '';
    if (!shown.length) box.innerHTML = `<span class="hint">${t('d_no_neigh')}</span>`;
    for (const m of shown) {
      const b = document.createElement('button');
      b.type = 'button'; b.className = 'chip'; b.setAttribute('aria-pressed', 'true');
      b.style.setProperty('--c', colorOf(m.i));
      b.innerHTML = `<i></i>${titleRes(D.residue[m.i])} ${D.pos[m.i] || D.auth[m.i]} · ${fmt(m.d)} Å`;
      b.addEventListener('click', () => select(m.i));
      box.appendChild(b);
    }
    if (!molIsOn()) { $('mol-msg').textContent = t('mol_size'); $('mol-load').hidden = false; }
  }
  function sendToMol() {
    const i = S.sel; if (i < 0) return;
    const code = D.pdb[i];
    molShow({ pdb: code, auth: D.auth[i], chain: D.chain[i] || D.proteinInfo[code]?.chain || 'A', residue: D.residue[i] });
  }
  $('mol-load').addEventListener('click', () => { $('mol-load').hidden = true; sendToMol(); });
  $('d-close').addEventListener('click', () => select(-1));

  /* Controles */
  function seg(host, items, get, set) {
    host.innerHTML = '';
    items.forEach(([v, key]) => {
      const b = document.createElement('button');
      b.type = 'button'; b.setAttribute('role', 'radio'); b.dataset.v = v; b.dataset.t = key;
      b.textContent = t(key);
      b.addEventListener('click', () => set(v));
      host.appendChild(b);
    });
    const sync = () => [...host.children].forEach(b => b.setAttribute('aria-checked', String(b.dataset.v === get())));
    sync();
    return sync;
  }
  const syncProj = seg($('sw-proj'), [['equal', 'proj_equal'], ['eqd', 'proj_eqd']], () => S.proj, v => { S.proj = v; syncProj(); reproject(true); apply({ keepSel: true }); });
  const syncFit = seg($('sw-fit'), [['catalytic', 'fit_catalytic'], ['core', 'fit_core']], () => S.fit, v => {
    if (v === 'core' && !hasCore) return;
    S.fit = v; syncFit(); reproject(true); apply({ keepSel: true }); if (S.sel >= 0) select(S.sel);
  });
  const syncColor = seg($('sw-color'), [['class', 'color_class'], ['domain', 'color_domain']], () => S.color, v => { S.color = v; syncColor(); renderLegend(); writeHash(); draw(); });
  const changeFam = v => {
    if (v === D.family) return;
    // Cambiar de familia limpia los filtros: posiciones y proteínas no valen en otra
    const q = new URLSearchParams();
    if (v !== 'GH1') q.set('fam', v.toLowerCase());
    if (S.proj !== DEFAULTS.proj) q.set('proj', S.proj);
    if (S.color !== DEFAULTS.color) q.set('color', S.color);
    const s = q.toString();
    history.replaceState(null, '', location.pathname + location.search + '#sky' + (s ? '?' + s : ''));
    location.reload();
  };
  seg($('sw-fam'), FAMILIES.map(f => [f, 'fam_label_' + f]), () => D.family, changeFam);
  [...$('sw-fam').children].forEach(b => { b.textContent = b.dataset.v; });
  $('sw-fit').hidden = !hasCore;
  const syncHemi = seg($('sw-hemi'), [['n', 'north'], ['s', 'south']], () => discsEl.dataset.show, v => { discsEl.dataset.show = v; syncHemi(); sizeDiscs(); });

  const dminI = $('f-dmin'), dmaxI = $('f-dmax');
  function onDist(e) {
    let a = +dminI.value, b = +dmaxI.value;
    if (a > b) { if (e.target === dminI) { b = a; dmaxI.value = b; } else { a = b; dminI.value = a; } }
    S.dmin = a; S.dmax = b; apply();
  }
  dminI.addEventListener('input', onDist); dmaxI.addEventListener('input', onDist);

  const resHost = $('f-res');
  function buildResidues() {
    resHost.innerHTML = '';
    CLASSES.forEach(c => {
      const g = document.createElement('div'); g.className = 'res-g'; g.style.setProperty('--c', CLASS_COLOR[c]);
      const head = document.createElement('button'); head.type = 'button'; head.innerHTML = `<i></i>${t('cls_' + c)}`;
      head.addEventListener('click', () => {
        const only = CLASS_RES[c].every(r => S.res.has(r)) && S.res.size === CLASS_RES[c].length;
        S.res = new Set(only ? allRes : CLASS_RES[c]); syncRes(); apply();
      });
      const box = document.createElement('div');
      CLASS_RES[c].forEach(r => {
        const b = document.createElement('button');
        b.type = 'button'; b.className = 'chip'; b.dataset.r = r; b.textContent = titleRes(r);
        b.addEventListener('click', () => { S.res.has(r) ? S.res.delete(r) : S.res.add(r); syncRes(); apply(); });
        box.appendChild(b);
      });
      g.append(head, box); resHost.appendChild(g);
    });
    syncRes();
  }
  const syncRes = () => resHost.querySelectorAll('.chip').forEach(b => b.setAttribute('aria-pressed', String(S.res.has(b.dataset.r))));
  $('f-res-all').addEventListener('click', () => { S.res = new Set(allRes); syncRes(); apply(); });

  const domHost = $('f-dom');
  function buildDomains() {
    domHost.innerHTML = '';
    const counts = DOMAINS.map(() => 0);
    D.pdbs.forEach((_, p) => { counts[D.dom[D.protRows[p][0]]]++; });
    DOMAINS.forEach((dm, k) => {
      const b = document.createElement('button');
      b.type = 'button'; b.className = 'chip'; b.dataset.d = dm; b.style.setProperty('--c', DOMAIN_COLOR[dm]);
      b.innerHTML = `<i></i>${t('dom_' + dm)} ${counts[k]}`;
      b.addEventListener('click', () => { S.doms.has(dm) ? S.doms.delete(dm) : S.doms.add(dm); syncDom(); apply(); });
      domHost.appendChild(b);
    });
    syncDom();
  }
  const syncDom = () => domHost.querySelectorAll('.chip').forEach(b => b.setAttribute('aria-pressed', String(S.doms.has(b.dataset.d))));

  const posI = $('f-pos');
  posI.addEventListener('input', () => { S.pos = posI.value; apply(); });

  const pdbI = $('f-pdb');
  function buildProteins() {
    const v = S.pdb;
    pdbI.innerHTML = `<option value="">${t('f_protein_all', { n: D.pdbs.length })}</option>` +
      [...D.pdbs].sort().map(c => { const o = D.organism(c); return `<option value="${c}">${c}${o ? ' · ' + o : ''}</option>`; }).join('');
    pdbI.value = v;
  }
  pdbI.addEventListener('change', () => { S.pdb = pdbI.value; apply(); });

  function syncAll() {
    dminI.value = S.dmin; dmaxI.value = S.dmax; posI.value = S.pos; pdbI.value = S.pdb;
    syncRes(); syncDom(); syncProj(); syncFit(); syncColor(); renderLegend();
  }
  function reset() {
    const prevFit = S.fit, prevProj = S.proj;
    Object.assign(S, DEFAULTS, { res: new Set(allRes), doms: new Set(DOMAINS) });
    select(-1);
    syncAll();
    if (prevFit !== S.fit || prevProj !== S.proj) reproject(true);
    brushClear();
    apply();
  }
  $('f-reset').addEventListener('click', reset);

  /* Leyenda */
  function renderLegend() {
    const host = $('legend');
    const items = S.color === 'class'
      ? CLASSES.map(c => [CLASS_COLOR[c], t('cls_' + c)])
      : DOMAINS.map(d => [DOMAIN_COLOR[d], t('dom_' + d)]);
    const dots = MAG_RADIUS.map((r, i) => `<circle cx="${6 + i * 13}" cy="7" r="${r}" fill="#ece4d2"/>`).join('');
    host.innerHTML = items.map(([c, l]) => `<span class="k" style="--c:${c}"><i></i>${l}</span>`).join('') +
      `<span class="k"><svg width="14" height="14" viewBox="-7 -7 14 14" aria-hidden="true"><path d="M0-6 6 0 0 6-6 0Z" fill="none" stroke="#ece4d2"/></svg>${t('ref_mark')}</span>` +
      `<span class="mag">${t('mag_near')} <svg width="76" height="14" aria-hidden="true">${dots}</svg> ${t('mag_far')}</span>`;
  }

  /* Histograma de distancias */
  const histCv = $('f-hist');
  const bins = new Array(DMAX).fill(0);
  { const dd = D.fits.catalytic.distance; for (let i = 0; i < n; i++) bins[Math.min(DMAX - 1, Math.floor(dd[i]))]++; }
  function drawHist() {
    if (!histCv.offsetParent) return;
    const { ctx, w, h } = fitCanvas(histCv);
    ctx.clearRect(0, 0, w, h);
    const mx = Math.max(...bins), bw = w / DMAX;
    bins.forEach((v, i) => {
      const on = i + 1 > S.dmin && i < (S.dmax >= DMAX ? Infinity : S.dmax);
      ctx.fillStyle = on ? '#d8b46a' : '#2b3347';
      const hh = Math.max(1, Math.sqrt(v / mx) * h);
      ctx.fillRect(i * bw, h - hh, Math.max(1, bw - 1), hh);
    });
  }

  /* Espectro */
  const specHost = $('spectrum');
  const specCv = specHost.querySelector('canvas'), specSvg = d3.select(specHost.querySelector('svg'));
  const allPos = [...D.positions.keys()].sort((a, b) => a - b);
  const pMin = allPos[0], pMax = allPos[allPos.length - 1];
  let sx = null, brush = null, brushing = false;
  if (!allPos.length) specHost.parentElement.hidden = true;
  function drawSpectrum() {
    if (!specCv.offsetParent || !allPos.length) return;
    const { ctx, w, h } = fitCanvas(specCv);
    ctx.clearRect(0, 0, w, h);
    sx = d3.scaleLinear([pMin - .5, pMax + .5], [0, w]);
    const ranges = parsePos(S.pos);
    const maxR = 8, base = h - 14;
    const bw = Math.max(1, w / (pMax - pMin + 1) - .4);
    for (const p of allPos) {
      const info = D.positions.get(p);
      const hh = Math.max(2, Math.min(1, info.rmsf / maxR) * (base - 2));
      const inR = !ranges || ranges.some(([a, b]) => p >= a && p <= b);
      const c = CLASS_COLOR[CLASSES[D.cls[(D.posRows.get(p) || [0])[0]] ?? 5]];
      const top = classColor(info.top_residue);
      ctx.fillStyle = top || c;
      ctx.globalAlpha = (inR ? 1 : .22) * (.4 + .6 * info.n / D.pdbs.length);
      ctx.fillRect(sx(p) - bw / 2, base - hh, bw, hh);
    }
    ctx.globalAlpha = 1;
    // Referencias y selección
    ctx.fillStyle = '#a3a9b4'; ctx.font = '10px "Spline Sans Mono", monospace'; ctx.textBaseline = 'alphabetic'; ctx.textAlign = 'center';
    for (const p of [100, 200, 300, 400]) if (p > pMin && p < pMax) ctx.fillText(p, sx(p), h - 2);
    ctx.strokeStyle = '#ece4d2'; ctx.lineWidth = 1;
    for (const p of D.refPositions) { ctx.beginPath(); ctx.moveTo(sx(p), 0); ctx.lineTo(sx(p), base); ctx.stroke(); }
    if (S.sel >= 0 && D.pos[S.sel] > 0) {
      ctx.strokeStyle = '#d8b46a'; ctx.beginPath(); ctx.moveTo(sx(D.pos[S.sel]), 0); ctx.lineTo(sx(D.pos[S.sel]), base); ctx.stroke();
    }
    if (!brush) {
      brush = d3.brushX().extent([[0, 0], [w, h]]).on('end', e => {
        if (!e.sourceEvent) return;
        if (!e.selection) {
          const [x] = d3.pointer(e.sourceEvent, specHost);
          const p = Math.round(sx.invert(x));
          S.pos = S.pos === String(p) ? '' : String(p);
        } else {
          const a = Math.round(sx.invert(e.selection[0])), b = Math.round(sx.invert(e.selection[1]));
          S.pos = a === b ? String(a) : `${a}-${b}`;
        }
        posI.value = S.pos; apply();
      });
      specSvg.call(brush);
      specSvg.selectAll('.overlay').attr('fill', 'transparent').style('cursor', 'crosshair');
    } else brush.extent([[0, 0], [w, h]]);
    if (!brushing) { brushing = true; specSvg.call(brush); brushing = false; }
  }
  function classColor(res) { const c = CLASSES.find(k => CLASS_RES[k].includes(ALIAS[res] || res)); return c ? CLASS_COLOR[c] : null; }
  function brushClear() { if (brush) specSvg.call(brush.move, null); }
  specHost.addEventListener('pointermove', e => {
    if (e.pointerType === 'touch' || !sx || e.buttons) return;
    const r = specHost.getBoundingClientRect();
    const p = Math.round(sx.invert(e.clientX - r.left));
    const info = D.positions.get(p);
    if (!info) { tip.hidden = true; return; }
    tip.innerHTML = `<b>${t('tip_pos', { p })}</b> · ${titleRes(info.top_residue)}<span>RMSF ${fmt(info.rmsf, 2)} Å · ${fmt(info.distance)} Å · n ${info.n}</span>`;
    tip.hidden = false;
    const host = discsEl.getBoundingClientRect();
    tip.style.left = Math.max(0, Math.min(host.width - tip.offsetWidth, e.clientX - host.left + 10)) + 'px';
    tip.style.top = (r.top - host.top - tip.offsetHeight - 4) + 'px';
  });
  specHost.addEventListener('pointerleave', () => { tip.hidden = true; });

  /* Tabla alternativa */
  function renderTable() {
    const f = D.fits[S.fit];
    const agg = new Map();
    for (const list of vis) for (const i of list) {
      const p = D.pos[i]; if (!p) continue;
      const a = agg.get(p) || { c: 0, s: 0, row: i };
      a.c++; a.s += f.distance[i]; if (f.distance[i] < f.distance[a.row]) a.row = i;
      agg.set(p, a);
    }
    const rows = [...agg].map(([p, a]) => ({ p, c: a.c, d: a.s / a.c, row: a.row, info: D.positions.get(p) })).sort((a, b) => a.d - b.d);
    $('tbl-note').textContent = t('tbl_note', { n: rows.length });
    $('tbl').innerHTML = `<thead><tr><th>${t('th_pos')}</th><th>${t('th_res')}</th><th class="num">${t('th_seen')}</th><th class="num">${t('th_dist')}</th><th class="num">${t('th_rmsf')}</th></tr></thead><tbody>` +
      rows.slice(0, 460).map(r => `<tr><td><button type="button" data-row="${r.row}">${r.p}</button></td><td>${r.info ? titleRes(r.info.top_residue) : '–'}</td><td class="num">${r.c}</td><td class="num">${fmt(r.d)}</td><td class="num">${r.info ? fmt(r.info.rmsf, 2) : '–'}</td></tr>`).join('') + '</tbody>';
  }
  $('tbl-box').addEventListener('toggle', () => { if ($('tbl-box').open) renderTable(); });
  $('tbl').addEventListener('click', e => { const b = e.target.closest('button[data-row]'); if (b) { select(+b.dataset.row); $('sky').scrollIntoView({ behavior: REDUCED ? 'auto' : 'smooth' }); } });

  /* CSV */
  $('csv').addEventListener('click', () => {
    const f = D.fits[S.fit];
    const out = ['pdb,auth_seq_id,hmm_position,residue,domain,theta,phi,distance,x,y,z'];
    for (const list of vis) for (const i of list) out.push([D.pdb[i], D.auth[i] ?? '', D.pos[i] || '', D.residue[i], DOMAINS[D.dom[i]], f.theta[i].toFixed(4), f.phi[i].toFixed(4), f.distance[i].toFixed(2), f.x[i].toFixed(2), f.y[i].toFixed(2), f.z[i].toFixed(2)].join(','));
    const url = URL.createObjectURL(new Blob([out.join('\n')], { type: 'text/csv' }));
    const a = document.createElement('a'); a.href = url; a.download = `easimap-${D.family.toLowerCase()}-${S.fit}.csv`; a.click();
    setTimeout(() => URL.revokeObjectURL(url), 2000);
  });

  /* Estado en la URL */
  let touched = false;
  function writeHash() {
    if (!touched) return;
    const q = new URLSearchParams();
    if (D.family !== 'GH1') q.set('fam', D.family.toLowerCase());
    if (S.dmin) q.set('dmin', S.dmin);
    if (S.dmax !== DEFAULTS.dmax) q.set('dmax', S.dmax);
    if (S.pos) q.set('pos', S.pos);
    if (S.pdb) q.set('pdb', S.pdb);
    if (S.fit !== DEFAULTS.fit) q.set('fit', S.fit);
    if (S.proj !== DEFAULTS.proj) q.set('proj', S.proj);
    if (S.color !== DEFAULTS.color) q.set('color', S.color);
    if (S.res.size !== allRes.length) q.set('res', [...S.res].join('.'));
    if (S.doms.size !== DOMAINS.length) q.set('dom', [...S.doms].join('.'));
    const s = q.toString();
    history.replaceState(null, '', location.pathname + location.search + '#sky' + (s ? '?' + s : ''));
  }
  function readHash() {
    const m = location.hash.match(/^#sky\?(.*)$/);
    if (!m) return false;
    const q = new URLSearchParams(m[1]);
    if (q.has('dmin')) S.dmin = +q.get('dmin') || 0;
    if (q.has('dmax')) S.dmax = +q.get('dmax') || DEFAULTS.dmax;
    if (q.has('pos')) S.pos = q.get('pos');
    if (q.has('pdb') && D.pdbIdx.has(q.get('pdb'))) S.pdb = q.get('pdb');
    if (q.get('fit') === 'core' && hasCore) S.fit = 'core';
    if (q.get('proj') === 'eqd') S.proj = 'eqd';
    if (q.get('color') === 'domain') S.color = 'domain';
    if (q.has('res')) S.res = new Set(q.get('res').split('.').filter(r => allRes.includes(r)));
    if (q.has('dom')) S.doms = new Set(q.get('dom').split('.').filter(r => DOMAINS.includes(r)));
    return true;
  }
  document.querySelector('.console').addEventListener('pointerdown', () => { touched = true; }, { capture: true });
  document.querySelector('.console').addEventListener('keydown', () => { touched = true; }, { capture: true });

  /* Encuadra lo visible cuando un filtro deja un grupo pequeño */
  function frameVisible() {
    for (const d of activeDiscs()) {
      const list = vis[d.h];
      const sel = d3.select(d.cv).transition().duration(REDUCED ? 0 : 500);
      if (list.length < 2) { sel.call(d.zoom.transform, d3.zoomIdentity); continue; }
      const xs = d3.extent(list, i => d.cx + toX[i] * d.R), ys = d3.extent(list, i => d.cy + toY[i] * d.R);
      const k = Math.max(1, Math.min(9, .62 * d.w / Math.max(xs[1] - xs[0], ys[1] - ys[0], 8)));
      const mx = (xs[0] + xs[1]) / 2, my = (ys[0] + ys[1]) / 2;
      sel.call(d.zoom.transform, d3.zoomIdentity.translate(d.w / 2, d.w / 2).scale(k).translate(-mx, -my));
    }
  }

  /* Puntos de entrada */
  function preset(p) {
    const prevFit = S.fit;
    const { frame, ...vals } = p;
    Object.assign(S, DEFAULTS, { res: new Set(allRes), doms: new Set(DOMAINS), proj: S.proj, fit: S.fit }, vals);
    select(-1);
    syncAll();
    if (prevFit !== S.fit) reproject(true);
    brushClear();
    touched = true;
    apply();
    if (p.frame) frameVisible(); else activeDiscs().forEach(d => d3.select(d.cv).call(d.zoom.transform, d3.zoomIdentity));
  }
  const tryOne = D.pdbIdx.has(D.fam.reference_pdb) ? D.fam.reference_pdb : (D.pdbIdx.has('1CBG') ? '1CBG' : [...D.pdbs].sort()[0]);
  const PRESETS = [
    ['try_1', { pos: D.refPositions.join(', '), dmax: DMAX }],
    ['try_2', { res: new Set(CLASS_RES.arom), dmax: 10 }],
    ['try_3', { pos: '333', res: new Set(['TRP']), dmax: DMAX, color: 'domain', frame: true }],
    ['try_4', { pdb: tryOne, dmax: DMAX }],
  ].filter(([k]) => k !== 'try_3' || D.family === 'GH1');
  function buildTry() {
    const host = $('try'); host.innerHTML = '';
    PRESETS.forEach(([k, p]) => {
      const b = document.createElement('button'); b.type = 'button'; b.textContent = t(k, { c: tryOne });
      b.addEventListener('click', () => preset(p));
      host.appendChild(b);
    });
  }

  /* Arranque */
  function rebuildText() {
    buildResidues(); buildDomains(); buildProteins(); buildTry(); renderLegend(); updateCount();
    [$('sw-proj'), $('sw-fit'), $('sw-color'), $('sw-hemi')].forEach(h => [...h.children].forEach(b => { b.textContent = t(b.dataset.t); }));
    $('sky-h').innerHTML = t('sky_title', { f: D.family });
    $('north-sub').textContent = t(D.family === 'GH5' ? 'north_sub_gh5' : 'north_sub');
    const euk = D.extras.stats?.by_domain?.Eukaryota ?? 0;
    const note = { GH5: ['note_gh5'], GH16: ['note_few'], GH18: ['note_num'], GH20: ['note_num', 'note_few'] }[D.family] || [];
    $('fam-note').hidden = !note.length;
    $('fam-note').textContent = note.map(k => t(k, { e: euk })).join(' ');
    if (S.sel >= 0) renderCard();
    if ($('tbl-box').open) renderTable();
    molMessage();
  }
  const hadHash = readHash();
  rebuildText();
  syncAll();
  reproject(false);
  apply({ silent: true });
  let rt = 0;
  new ResizeObserver(() => { clearTimeout(rt); rt = setTimeout(sizeDiscs, 60); }).observe(discsEl);
  sizeDiscs();
  document.fonts?.ready.then(draw);
  onLang(rebuildText);
  if (hadHash) {
    touched = true;
    setTimeout(() => { if (S.pos) frameVisible(); }, 60);
  }

  return { preset, state: S, select, vis: () => vis };
}
