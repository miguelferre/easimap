// Historias y método. Todas las cifras salen de los datos.
import { D, t, onLang, radialFn, fmt, fmtInt, DOMAINS, DOMAIN_COLOR, drawDiscBase } from './core.js';

const $ = id => document.getElementById(id);
const SHAPES = { Eukaryota: d3.symbolCircle, Bacteria: d3.symbolSquare, Archaea: d3.symbolTriangle, Unknown: d3.symbolDiamond };

export function initStories() {
  const st = D.extras.story333, off = D.extras.offsets;
  if (st) buildFig333(st);
  if (D.fam.horizon) buildHorizon();
  const copy = () => {
    if (st) copy333(st, off);
    if (off) buildOffset(off);
    if (D.fam.horizon) copyHorizon();
    buildMethod();
    if (st) legend333(st);
  };
  copy();
  onLang(copy);
  let rt = 0, lastW = $('fig-offset').clientWidth;
  new ResizeObserver(() => {
    if ($('fig-offset').clientWidth === lastW) return;
    lastW = $('fig-offset').clientWidth;
    clearTimeout(rt); rt = setTimeout(() => { if (off) buildOffset(off); }, 120);
  }).observe($('fig-offset'));
}

/* Posición 333 */
function copy333(st, off) {
  const A = st.groups.A.by_domain, B = st.groups.B.by_domain;
  const total = D.pdbs.length, nw = st.new_only;
  $('st1-p1').textContent = t('st1_p1', { n: st.n_trp, t: total });
  $('st1-p2').textContent = t('st1_p2', {
    sep: fmt(st.centroid_separation, 2), ae: A.Eukaryota || 0, bb: B.Bacteria || 0, ba: B.Archaea || 0, bu: B.Unknown || 0,
    k: st.fit_all.known, ta: fmt(st.groups.A.theta_mean_deg), tb: fmt(st.groups.B.theta_mean_deg),
  });
  $('st1-p3').textContent = t('st1_p3', {
    n: st.n_trp, nn: nw.n_trp, ig: nw.in_group, kn: nw.known, sep: fmt(nw.centroid_separation, 2), r: fmt(st.rmsf.current, 2),
  });
  $('st1-p4').textContent = st.chimera ? t('st1_p4', { c: st.chimera.pdb }) : '';
  const pTxt = st.p_perm < 1e-4 ? fmt(.0001, 4) : fmt(Math.max(st.p_perm, 1 / st.n_perm), 4);
  $('st1-limit').textContent = t('st1_limit', { p: pTxt, a: st.fit_all.in_group, k: st.fit_all.known, c: st.chimera?.pdb || '' });
  if (off) $('st1b-hyp').textContent = t('st1b_hyp', { k: off.pairs?.length ?? 0, h: off.hits.length, t: off.tested });
  $('st1b-p').textContent = t('st1b_p');
  const rows = [];
  for (const g of ['A', 'B']) for (const d of DOMAINS) { const c = st.groups[g].by_domain[d]; if (c) rows.push([t('grp_' + g.toLowerCase()), t('dom_' + d), c]); }
  let html = `<thead><tr><th>${t('th_group')}</th><th>${t('th_domain')}</th><th class="num">${t('th_count')}</th></tr></thead><tbody>` +
    rows.map(r => `<tr><td>${r[0]}</td><td>${r[1]}</td><td class="num">${r[2]}</td></tr>`).join('') + '</tbody>';
  if (off) html += `<thead><tr><th>${t('th_pos')}</th><th class="num">${t('th_same')}</th><th class="num">${t('th_back')}</th></tr></thead><tbody>` +
    off.table.map(r => `<tr><td>${r.p}</td><td class="num">${fmt(r.k0)}</td><td class="num">${fmt(r.km1)}</td></tr>`).join('') + '</tbody>';
  $('tbl-333').innerHTML = html;
  d3.select('#fig-333 svg').attr('aria-label', t('fig333_alt'));
}
function legend333(st) {
  const items = DOMAINS.filter(d => st.by_domain[d]).map(d => {
    const path = d3.symbol(SHAPES[d], 46)();
    return `<span class="k"><svg width="12" height="12" viewBox="-6 -6 12 12" aria-hidden="true"><path d="${path}" fill="${DOMAIN_COLOR[d]}"/></svg>${t('dom_' + d)} <b class="num">${st.by_domain[d]}</b></span>`;
  });
  $('leg-333').innerHTML = items.join('');
  d3.selectAll('#fig-333 text.g-a').text(t('grp_a'));
  d3.selectAll('#fig-333 text.g-b').text(t('grp_b'));
}
function buildFig333(st) {
  const rad = radialFn('equal');
  const W = 720, H = 470, pad = 54;
  const pts = st.points.map(p => {
    const r = rad(p.theta);
    return { ...p, u: -r * Math.sin(p.phi), v: -r * Math.cos(p.phi) };
  });
  const xs = d3.extent(pts, d => d.u), ys = d3.extent(pts, d => d.v);
  const k = Math.min((W - 2 * pad) / (xs[1] - xs[0]), (H - 2 * pad) / (ys[1] - ys[0]));
  const ox = W / 2 - k * (xs[0] + xs[1]) / 2, oy = H / 2 - k * (ys[0] + ys[1]) / 2;
  const X = u => ox + k * u, Y = v => oy + k * v;
  const svg = d3.select('#fig-333').append('svg').attr('viewBox', `0 0 ${W} ${H}`).attr('role', 'img').style('width', '100%').style('height', 'auto').style('display', 'block');

  // Retícula real del mapa en esta zona
  const grid = svg.append('g').attr('fill', 'none').attr('stroke', '#2b3347');
  const line = d3.line().x(d => d[0]).y(d => d[1]);
  for (let a = 20; a <= 90; a += 5) {
    const r = rad(a * Math.PI / 180), seg = [];
    for (let ph = 40; ph <= 140; ph += 1) seg.push([X(-r * Math.sin(ph * Math.PI / 180)), Y(-r * Math.cos(ph * Math.PI / 180))]);
    grid.append('path').attr('d', line(seg));
  }
  for (let ph = 40; ph <= 140; ph += 5) {
    const seg = [];
    for (let a = 20; a <= 90; a += 2) { const r = rad(a * Math.PI / 180); seg.push([X(-r * Math.sin(ph * Math.PI / 180)), Y(-r * Math.cos(ph * Math.PI / 180))]); }
    grid.append('path').attr('d', line(seg));
  }
  svg.append('clipPath').attr('id', 'c333').append('rect').attr('width', W).attr('height', H);
  grid.attr('clip-path', 'url(#c333)');

  // Envolvente de cada grupo: línea fina, como las figuras del atlas
  for (const g of ['A', 'B']) {
    const q = pts.filter(p => p.group === g).map(p => [X(p.u), Y(p.v)]);
    const hull = d3.polygonHull(q);
    if (!hull) continue;
    const cx = d3.mean(hull, d => d[0]), cy = d3.mean(hull, d => d[1]);
    const grown = hull.map(([x, y]) => { const dx = x - cx, dy = y - cy, l = Math.hypot(dx, dy) || 1; return [x + dx / l * 16, y + dy / l * 16]; });
    svg.append('path').attr('d', d3.line().curve(d3.curveCatmullRomClosed.alpha(.7))(grown))
      .attr('fill', 'none').attr('stroke', '#d8b46a').attr('stroke-width', 1).attr('stroke-dasharray', '1 4').attr('stroke-linecap', 'round');
    const top = grown.reduce((a, b) => (b[1] < a[1] ? b : a));
    const isA = g === 'A';
    svg.append('text').attr('class', isA ? 'g-a' : 'g-b').attr('x', top[0] + (isA ? 24 : -24)).attr('y', Math.max(22, top[1] + 12)).attr('text-anchor', isA ? 'start' : 'end')
      .attr('fill', '#ece4d2').style('font', `italic 500 ${innerWidth < 600 ? 38 : 22}px "Bodoni Moda", serif`);
  }

  const tip = d3.select('#fig-333').style('position', 'relative').append('div').attr('class', 'tip').attr('hidden', true);
  svg.append('g').selectAll('path').data(pts).join('path')
    .attr('d', d => d3.symbol(SHAPES[d.domain] || d3.symbolCircle, 74)())
    .attr('transform', d => `translate(${X(d.u)},${Y(d.v)})`)
    .attr('fill', d => DOMAIN_COLOR[d.domain]).attr('stroke', '#080c16').attr('stroke-width', 1.2)
    .on('pointerenter', (e, d) => {
      const r = e.currentTarget.ownerSVGElement.getBoundingClientRect(), s = r.width / W;
      tip.attr('hidden', null).html(`<b>${d.pdb}</b><span>${d.organism || ''}</span><span>${fmt(d.distance)} Å</span>`)
        .style('left', Math.min(r.width - 220, X(d.u) * s + 12) + 'px').style('top', (Y(d.v) * s + 12) + 'px');
    })
    .on('pointerleave', () => tip.attr('hidden', true));
  const ch = st.chimera && pts.find(p => p.pdb === st.chimera.pdb);
  if (ch) {
    svg.append('text').attr('x', X(ch.u) + 10).attr('y', Y(ch.v) - 8).attr('fill', '#ece4d2').style('font-size', innerWidth < 600 ? '22px' : '13px').text(ch.pdb);
  }
}

/* Desfase 333-337: dos puntos por posición sobre una escala común en Å */
function buildOffset(off) {
  const host = d3.select('#fig-offset'); host.selectAll('*').remove();
  const rows = off.table;
  const W = Math.max(300, Math.min(760, host.node().clientWidth || 640)), rowH = 40, top = 34, left = 44, right = 30, H = top + rows.length * rowH + 26;
  const x = d3.scaleLinear([0, 9], [left, W - right]);
  const svg = host.append('svg').attr('viewBox', `0 0 ${W} ${H}`).attr('role', 'img').attr('aria-label', t('off_axis')).style('width', '100%').style('height', 'auto').style('display', 'block');
  const axisY = top + rows.length * rowH + 4;
  svg.append('g').selectAll('line').data(d3.range(0, 10, 1)).join('line')
    .attr('x1', d => x(d)).attr('x2', d => x(d)).attr('y1', top - 8).attr('y2', axisY).attr('stroke', '#2b3347');
  svg.append('g').selectAll('text').data(d3.range(0, 10, 3)).join('text')
    .attr('x', d => x(d)).attr('y', axisY + 15).attr('text-anchor', 'middle').attr('fill', '#a3a9b4').style('font-size', '12px').text(d => d);
  host.append('figcaption').attr('class', 'hint').text(t('off_axis'));
  // Leyenda directa
  svg.append('circle').attr('cx', left + 5).attr('cy', 12).attr('r', 5.5).attr('fill', '#ece4d2');
  svg.append('text').attr('x', left + 16).attr('y', 16).attr('fill', '#ece4d2').style('font-size', '12px').text(t('off_same'));
  const lx = left + 20 + t('off_same').length * 7.6 + 18;
  svg.append('circle').attr('cx', lx).attr('cy', 12).attr('r', 5.5).attr('fill', 'none').attr('stroke', '#d8b46a').attr('stroke-width', 2);
  svg.append('text').attr('x', lx + 11).attr('y', 16).attr('fill', '#ece4d2').style('font-size', '12px').text(t('off_back'));
  rows.forEach((r, i) => {
    const y = top + i * rowH + rowH / 2;
    const better = r.km1 < r.k0;
    svg.append('text').attr('x', left - 12).attr('y', y + 4).attr('text-anchor', 'end').attr('fill', better ? '#ece4d2' : '#a3a9b4').style('font-size', '13px').text(r.p);
    svg.append('line').attr('x1', x(r.k0)).attr('x2', x(r.km1)).attr('y1', y).attr('y2', y).attr('stroke', better ? '#d8b46a' : '#465069');
    svg.append('circle').attr('cx', x(r.k0)).attr('cy', y).attr('r', 5.5).attr('fill', '#ece4d2');
    svg.append('circle').attr('cx', x(r.km1)).attr('cy', y).attr('r', 5.5).attr('fill', '#080c16').attr('stroke', '#d8b46a').attr('stroke-width', 2);
    const lo = Math.min(r.k0, r.km1), hi = Math.max(r.k0, r.km1);
    svg.append('text').attr('x', x(lo) - 10).attr('y', y + 4).attr('text-anchor', 'end').attr('fill', '#a3a9b4').style('font-size', '11px').text(lo > .6 ? fmt(lo) : '');
    svg.append('text').attr('x', x(hi) + 10).attr('y', y + 4).attr('fill', '#a3a9b4').style('font-size', '11px').text(fmt(hi));
  });
}

/* Horizonte: bandas con d3.contours, recortadas al disco */
const BANDS = [8, 12, 16, 22];
const HZ = ['#10141f', '#1d2130', '#3b3a3c', '#8a7248', '#f2d9a0'];
let hzStats = null;

function buildHorizon() {
  const { grid: G, north, south } = D.fam.horizon;
  const host = d3.select('#fig-horizon');
  hzStats = { cells: [0, 0, 0, 0, 0], perHemi: [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], best: null, total: 0 };
  [['n', north], ['s', south]].forEach(([h, vals], hi) => {
    // Fuera del disco la malla vale 0: se rellena con la celda válida más cercana para que el contorno no se rompa en el borde
    const v = Array.from(vals);
    const inside = i => vals[i] > 0;
    for (let j = 0; j < G; j++) for (let i = 0; i < G; i++) {
      const id = j * G + i;
      if (inside(id)) {
        const b = BANDS.findIndex(x => vals[id] < x);
        const k = b < 0 ? 4 : b;
        hzStats.cells[k]++; hzStats.perHemi[hi][k]++; hzStats.total++;
        if (!hzStats.best || vals[id] > hzStats.best.d) hzStats.best = { d: vals[id], h, i, j };
        continue;
      }
      let bd = 1e9, bv = 0;
      for (let q = 0; q < G * G; q++) if (inside(q)) { const dd = (q % G - i) ** 2 + (Math.floor(q / G) - j) ** 2; if (dd < bd) { bd = dd; bv = vals[q]; } }
      v[id] = bv;
    }
    const S = 300, R = 138, c = S / 2;
    const fig = host.append('figure');
    const svg = fig.append('svg').attr('viewBox', `0 0 ${S} ${S}`).attr('role', 'img').attr('data-h', h);
    const id = 'hz-' + h;
    svg.append('clipPath').attr('id', id).append('circle').attr('cx', c).attr('cy', c).attr('r', R);
    // Malla fina por interpolación bilineal y un suavizado leve: solo para dibujar las curvas.
    // La lectura al pasar usa siempre las celdas originales.
    const U = 5, F = G * U, fine = new Float64Array(F * F);
    for (let y = 0; y < F; y++) for (let x = 0; x < F; x++) {
      const gx = Math.min(G - 1.001, Math.max(0, (x + .5) / U - .5)), gy = Math.min(G - 1.001, Math.max(0, (y + .5) / U - .5));
      const i0 = Math.floor(gx), j0 = Math.floor(gy), fx = gx - i0, fy = gy - j0;
      fine[y * F + x] = v[j0 * G + i0] * (1 - fx) * (1 - fy) + v[j0 * G + i0 + 1] * fx * (1 - fy) + v[(j0 + 1) * G + i0] * (1 - fx) * fy + v[(j0 + 1) * G + i0 + 1] * fx * fy;
    }
    d3.blur2({ data: fine, width: F, height: F }, 3);
    const cont = d3.contours().size([F, F]).smooth(true).thresholds([0, ...BANDS])(fine);
    const sc = 2 * R / G, sf = 2 * R / F;
    const proj = d3.geoTransform({ point(x, y) { this.stream.point(c - R + x * sf, c - R + y * sf); } });
    const path = d3.geoPath(proj);
    const g = svg.append('g').attr('clip-path', `url(#${id})`);
    g.append('circle').attr('cx', c).attr('cy', c).attr('r', R).attr('fill', HZ[0]);
    // De cerca a lejos: cada umbral pinta lo que queda más lejos
    cont.forEach((m, i) => {
      if (i === 0) return;
      g.append('path').attr('d', path(m)).attr('fill', HZ[i]).attr('stroke', i > 2 ? '#d8b46a' : '#465069').attr('stroke-width', .5);
    });
    // Retícula del atlas encima: 30° y 60°
    const rr = radialFn('equal');
    for (const a of [30, 60]) svg.append('circle').attr('cx', c).attr('cy', c).attr('r', rr(a * Math.PI / 180) * R).attr('fill', 'none').attr('stroke', '#a3a9b4').attr('stroke-opacity', .28).attr('stroke-width', .6);
    for (let a = 0; a < 180; a += 30) { const sx = Math.sin(a * Math.PI / 180) * R, sy = Math.cos(a * Math.PI / 180) * R; svg.append('line').attr('x1', c - sx).attr('y1', c - sy).attr('x2', c + sx).attr('y2', c + sy).attr('stroke', '#a3a9b4').attr('stroke-opacity', .16).attr('stroke-width', .6); }
    svg.append('circle').attr('cx', c).attr('cy', c).attr('r', R).attr('fill', 'none').attr('stroke', '#d8b46a');
    svg.append('circle').attr('cx', c).attr('cy', c).attr('r', R + 6).attr('fill', 'none').attr('stroke', '#d8b46a').attr('stroke-opacity', .4).attr('stroke-width', .7);
    fig.append('figcaption');
    // Lectura al pasar
    const tip = host.selectAll('.tip').empty() ? host.append('div').attr('class', 'tip').attr('hidden', true) : host.select('.tip');
    svg.on('pointermove', e => {
      const r = svg.node().getBoundingClientRect();
      const px = (e.clientX - r.left) / r.width * S, py = (e.clientY - r.top) / r.height * S;
      if (Math.hypot(px - c, py - c) > R) { tip.attr('hidden', true); return; }
      const i = Math.min(G - 1, Math.max(0, Math.floor((px - (c - R)) / sc))), j = Math.min(G - 1, Math.max(0, Math.floor((py - (c - R)) / sc)));
      const val = vals[j * G + i];
      if (!(val > 0)) { tip.attr('hidden', true); return; }
      const hr = host.node().getBoundingClientRect();
      tip.attr('hidden', null).html(`<b>${fmt(val)} Å</b>`).style('left', (e.clientX - hr.left + 12) + 'px').style('top', (e.clientY - hr.top + 12) + 'px');
    }).on('pointerleave', () => tip.attr('hidden', true));
  });
  // Cruz en la vista más amplia
  const b = hzStats.best, S = 300, R = 138, c = S / 2, sc = 2 * R / G;
  const svg = d3.select(`#fig-horizon svg[data-h="${b.h}"]`);
  const x = c - R + (b.i + .5) * sc, y = c - R + (b.j + .5) * sc;
  svg.append('path').attr('d', `M${x - 8} ${y}H${x + 8}M${x} ${y - 8}V${y + 8}`).attr('stroke', '#080c16').attr('stroke-width', 1.6);
}
function copyHorizon() {
  const { grid: G } = D.fam.horizon;
  const s = hzStats;
  const pct = k => fmt(100 * k / s.total, 0) + ' %';
  $('st2-p2').textContent = t('st2_p2', { w: pct(s.cells[0] + s.cells[1]), o: pct(s.cells[4]), h: t(s.best.h === 'n' ? 'hz_north' : 'hz_south') });
  document.querySelectorAll('[data-t="st2_limit"]').forEach(el => { el.textContent = t('st2_limit', { g: G }); });
  const names = ['n', 's'];
  document.querySelectorAll('#fig-horizon figure').forEach((f, i) => {
    f.querySelector('figcaption').textContent = t(i ? 'south' : 'north');
    f.querySelector('svg').setAttribute('aria-label', t('hz_alt', { h: t(names[i] === 'n' ? 'hz_north' : 'hz_south') }));
  });
  const labels = [t('hz_less', { n: 8 }), '8–12 Å', '12–16 Å', '16–22 Å', t('hz_more', { n: 22 })];
  $('leg-horizon').innerHTML = `<span>${t('hz_legend')}</span>` + labels.map((l, i) => `<span class="k"><i style="--c:${HZ[i]};border-radius:0;outline:1px solid #465069"></i>${l}</span>`).join('') +
    `<span class="k">+ ${t('hz_widest')}</span>`;
  $('tbl-horizon').innerHTML = `<thead><tr><th>${t('th_band')}</th><th class="num">${t('th_share')}</th><th class="num">${t('th_north')}</th><th class="num">${t('th_south')}</th></tr></thead><tbody>` +
    labels.map((l, i) => `<tr><td>${l}</td><td class="num">${pct(s.cells[i])}</td><td class="num">${s.perHemi[0][i]}</td><td class="num">${s.perHemi[1][i]}</td></tr>`).join('') + '</tbody>';
}

/* Método */
function buildMethod() {
  const [a, b] = D.refPositions;
  const vars = { p: D.pdbs.length, r: fmtInt(D.n), a, b, c: D.extras.stats?.core_size ?? '' };
  $('m-data').textContent = t('m_data_pipe', vars);
  const R = D.extras.replica || {};
  $('m-replica').textContent = t('m_replica', { e: fmt(R.median_error, 5), p: fmtInt(R.pairs ?? 0), n: R.proteins, m: fmt(100 * (R.hmm_match?.same ?? 0) / (R.hmm_match?.total || 1)) });
  const fams = D.families || [];
  $('m-families').textContent = t('m_families', { n: fams.length });
  const keys = ['author', 'tutor', 'prior', 'related', 'struct', 'family', 'biology', 'afdb', 'viewer', 'maps', 'fonts'];
  $('m-list').innerHTML = keys.map(k => `<dt>${t('c_' + k)}</dt><dd>${t('c_' + k + '_v')}</dd>`).join('');
}
