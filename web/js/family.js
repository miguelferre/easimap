// Quién vive en este cielo: qué es GH1, qué hace y ocho cielos de proteínas reales.
import { D, t, onLang, fitCanvas, drawDiscBase, magOf, MAG_RADIUS, CLASSES, CLASS_COLOR, DOMAINS, fmt } from './core.js';

const $ = id => document.getElementById(id);
const SKIES = ['3SCU', '1CBG', '2E9L', '1WCG', '3WQ4', '3AHZ', '1GOW', '3W53'];
const linkOf = s => `visor.html#sky?pdb=${s.pdb}&amp;dmax=45`;

export function initFamily() {
  const fit = D.fits.catalytic;
  // Cada cielo es la propia entrada, otra de su grupo al 90 % (extras.skies)
  const skies = SKIES.map(c => ({ key: c, ...(D.extras.skies?.[c] || { pdb: c }) }));

  /* Esquema del corte, sobre los ejes del mapa: θ y distancia de 165 y 360 son medias reales */
  const mean = p => {
    const rows = D.posRows.get(p);
    return { th: d3.mean(rows, i => fit.theta[i]), d: d3.mean(rows, i => fit.distance[i]) };
  };
  function scheme() {
    const W = 560, H = 430, cx = 310, cy = 216, S = 38;
    const [pa, pn] = D.refPositions;
    const A = mean(pa), N = mean(pn);
    const at = q => [cx - q.d * S * Math.sin(q.th), cy - q.d * S * Math.cos(q.th)];
    const [ax, ay] = at(A), [nx, ny] = at(N);
    const ox = cx, oy = cy - 1.43 * S;
    const hex = (x, y, r, rot) => d3.range(6).map(i => { const a = rot + i * Math.PI / 3; return `${x + r * Math.cos(a)},${y + r * Math.sin(a)}`; }).join(' ');
    const lab = (x, y, a, b, anchor = 'start') => `<text x="${x}" y="${y}" text-anchor="${anchor}" class="s-t">${a}</text>` + (b ? `<text x="${x}" y="${y + 16}" text-anchor="${anchor}" class="s-m">${b}</text>` : '');
    $('fam-fig').innerHTML = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${t('sch_label')}">
      <line x1="24" x2="${W - 8}" y1="${cy}" y2="${cy}" class="s-eq"/>
      <text x="${W - 8}" y="${cy - 8}" text-anchor="end" class="s-m">${t('north')}</text>
      <text x="${W - 8}" y="${cy + 18}" text-anchor="end" class="s-m">${t('south')}</text>
      <polygon points="${hex(cx + 52, cy, 52, Math.PI)}" class="s-ring"/>
      <polygon points="${hex(ox + 62, oy - 66, 44, Math.PI / 2)}" class="s-ring"/>
      <line x1="${cx}" y1="${cy}" x2="${ox}" y2="${oy}" class="s-bond"/>
      <line x1="${ox}" y1="${oy}" x2="${ox + 62 - 38}" y2="${oy - 66 + 22}" class="s-bond"/>
      <path d="M${cx - 13} ${cy - 22}L${cx + 13} ${cy - 32}" class="s-cut"/>
      <circle cx="${ox}" cy="${oy}" r="6" class="s-o"/>
      <circle cx="${cx}" cy="${cy}" r="4" class="s-c"/><circle cx="${cx}" cy="${cy}" r="9" class="s-c2"/>
      <line x1="${ax}" y1="${ay}" x2="${ox - 9}" y2="${oy - 2}" class="s-arrow" marker-end="url(#fa)"/>
      <line x1="${nx}" y1="${ny}" x2="${cx - 7}" y2="${cy + 8}" class="s-arrow" marker-end="url(#fa)"/>
      <defs><marker id="fa" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0L8 4L0 8Z" fill="#d8b46a"/></marker></defs>
      <circle cx="${ax}" cy="${ay}" r="6" class="s-res"/><circle cx="${nx}" cy="${ny}" r="6" class="s-res"/>
      ${lab(ax - 10, ay - 26, `Glu ${pa}`, t('sch_acid'), 'middle')}
      ${lab(nx - 10, ny + 28, `Glu ${pn}`, t('sch_nuc'), 'middle')}
      ${lab(cx + 16, cy + 66, t("sch_c"), t("sch_c2"))}
      ${lab(ox + 15, oy + 1, t('sch_o'))}
      <text x="${cx + 52}" y="${cy + 5}" text-anchor="middle" class="s-m">−1</text>
      <text x="${ox + 62}" y="${oy - 61}" text-anchor="middle" class="s-m">+1</text>
    </svg>`;
    $('fam-cap').textContent = t('sch_cap', { a: fmt(A.th * 180 / Math.PI, 0), b: fmt(N.th * 180 / Math.PI, 0), d: fmt(A.d) });
  }

  /* Ocho cielos */
  function list() {
    $('fam-skies').innerHTML = skies.map(s => `<figure><a href="${linkOf(s)}" aria-label="${t('sk_' + s.key + '_t')}, ${s.pdb}"><canvas data-c="${s.pdb}" aria-hidden="true"></canvas></a><figcaption><h4>${t('sk_' + s.key + '_t')}</h4><p class="sk-code num">${s.pdb}${s.set === '2020' ? ' · ' + t('sk_tag_2020') : ''}</p><p>${t('sk_' + s.key + '_p')}</p></figcaption></figure>`).join('');
    draw();
  }
  function draw() {
    document.querySelectorAll('#fam-skies canvas').forEach(cv => {
      const { ctx, w, h } = fitCanvas(cv);
      ctx.clearRect(0, 0, w, h);
      if (w < 30) return;
      const R = w / 2 - 7, x0 = w / 2, y0 = h / 2;
      drawDiscBase(ctx, x0, y0, R, { hemi: 'n', labels: false });
      const Q = D.struct(cv.dataset.c);
      ctx.globalCompositeOperation = 'lighter';
      for (let i = 0; i < Q.n; i++) {
        if (Q.hemi[i]) continue;
        const m = magOf(Q.distance[i]);
        ctx.fillStyle = CLASS_COLOR[CLASSES[Q.cls[i]]];
        ctx.beginPath(); ctx.arc(x0 + Q.ux[i] * R, y0 + Q.uy[i] * R, MAG_RADIUS[m] * Math.min(1, w / 230), 0, 6.3); ctx.fill();
      }
      ctx.globalCompositeOperation = 'source-over';
    });
  }

  function copy() {
    const c = DOMAINS.map(() => 0);
    D.pdbs.forEach((_, p) => { c[D.dom[D.protRows[p][0]]]++; });
    $('fam-count').textContent = t('fam_count', { p: D.pdbs.length, e: c[0], b: c[1], a: c[2], u: c[3] });
    scheme(); list();
  }
  copy();
  onLang(copy);
  document.fonts?.ready.then(scheme);
  let rt = 0;
  new ResizeObserver(() => { clearTimeout(rt); rt = setTimeout(draw, 80); }).observe($('fam-skies'));
}
