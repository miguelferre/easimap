// Seis cielos: una entrada por familia, con su cielo real y lo que salió de verdad.
import { D, t, onLang, fitCanvas, drawDiscBase, CLASSES, CLASS_COLOR, MAG_RADIUS, fmt, fmtInt } from './core.js';

const $ = id => document.getElementById(id);
const link = f => (f.key === 'GH1' ? 'visor.html' : `visor.html#sky?fam=${f.file}`);

function decode(b64) {
  const s = atob(b64), n = s.length / 3, o = new Array(n);
  for (let i = 0; i < n; i++) o[i] = [s.charCodeAt(i * 3) / 127.5 - 1, s.charCodeAt(i * 3 + 1) / 127.5 - 1, s.charCodeAt(i * 3 + 2)];
  return o;
}

export function initSkies() {
  const fams = D.families || [];
  if (!fams.length) { $('six').hidden = true; return; }
  const pts = new Map(fams.map(f => [f.key, decode(f.sky)]));

  function draw() {
    document.querySelectorAll('#six-grid canvas').forEach(cv => {
      const { ctx, w, h } = fitCanvas(cv);
      ctx.clearRect(0, 0, w, h);
      if (w < 30) return;
      const R = w / 2 - 7, x0 = w / 2, y0 = h / 2;
      drawDiscBase(ctx, x0, y0, R, { hemi: 'n', labels: false });
      ctx.globalCompositeOperation = 'lighter';
      const k = Math.min(1, w / 240);
      for (const [u, v, c] of pts.get(cv.dataset.k)) {
        ctx.globalAlpha = .42;
        ctx.fillStyle = CLASS_COLOR[CLASSES[c >> 3]];
        ctx.beginPath(); ctx.arc(x0 + u * R, y0 + v * R, MAG_RADIUS[c & 7] * k * .8, 0, 6.3); ctx.fill();
      }
      ctx.globalAlpha = 1; ctx.globalCompositeOperation = 'source-over';
    });
  }

  const euk = f => f.by_domain.Eukaryota || 0;
  const pro = f => (f.by_domain.Bacteria || 0) + (f.by_domain.Archaea || 0);
  function result(f) {
    const c = f.candidates;
    if (f.key === 'GH1') return t('res_GH1');
    if (!c) return t('res_' + f.key);
    return t('res_' + f.key, { q: c.q05, n: c.family05.length, e: c.euk, p: c.family05[0] ?? '' });
  }
  function copy() {
    $('six-grid').innerHTML = fams.map(f => `<figure>
      <a href="${link(f)}" aria-label="${t('six_open', { f: f.key })}"><canvas data-k="${f.key}" aria-hidden="true"></canvas></a>
      <figcaption>
        <h3>${f.key}</h3>
        <p class="sk-code num">${t('six_stat', { n: f.structures, r: fmtInt(f.residues) })}</p>
        <p>${t('six_' + f.key)}</p>
        <p class="six-go"><a href="${link(f)}">${t('six_open', { f: f.key })}</a></p>
      </figcaption></figure>`).join('');
    $('six-table').innerHTML = `<thead><tr><th>${t('six_th_f')}</th><th class="num">${t('six_th_n')}</th><th class="num">${t('six_th_ep')}</th><th class="num">${t('six_th_rmsf')}</th><th>${t('six_th_res')}</th></tr></thead><tbody>` +
      fams.map(f => `<tr><td>${f.key}</td><td class="num" data-l="${t('six_th_n')}">${f.structures}</td><td class="num" data-l="${t('six_th_ep_s')}">${euk(f)} / ${pro(f)}</td><td class="num" data-l="RMSF">${fmt(f.median_rmsf, 2)} Å</td><td class="six-res">${result(f)}</td></tr>`).join('') + '</tbody>';
    const g1 = fams.find(f => f.key === 'GH1'), g5 = fams.find(f => f.key === 'GH5');
    $('six-frank').textContent = t('six_frank', { a: fmt(g5.median_rmsf, 2), b: fmt(g1.median_rmsf, 2) });
    draw();
  }
  copy();
  onLang(copy);
  let rt = 0;
  new ResizeObserver(() => { clearTimeout(rt); rt = setTimeout(draw, 80); }).observe($('six-grid'));
}
