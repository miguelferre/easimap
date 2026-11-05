// Arranque de las dos páginas: inicio y visor
import { applyStatic, loadData, loadFamilies, setLang, onLang, t, D, FAM } from './core.js';

const page = document.body.dataset.page;
applyStatic();
document.querySelectorAll('[data-lang]').forEach(b => b.addEventListener('click', () => setLang(b.dataset.lang)));

const bar = document.querySelector('.bar');
const onScroll = () => bar.classList.toggle('is-scrolled', scrollY > 24);
addEventListener('scroll', onScroll, { passive: true });
onScroll();

// La sección visible queda marcada en la navegación
const links = new Map([...document.querySelectorAll('.bar-nav a[href^="#"]')].map(a => [a.getAttribute('href').slice(1), a]));
const spy = new IntersectionObserver(es => es.forEach(e => {
  const a = links.get(e.target.id);
  if (a) { if (e.isIntersecting) a.setAttribute('aria-current', 'true'); else a.removeAttribute('aria-current'); }
}), { rootMargin: '-45% 0px -50% 0px' });
links.forEach((_, id) => { const s = document.getElementById(id); if (s) spy.observe(s); });

await new Promise(ok => { if (window.d3) ok(); else addEventListener('DOMContentLoaded', ok, { once: true }); });

const foot = () => { document.getElementById('foot-data').textContent = t('foot_data', { s: D.source, d: (D.generated || '').slice(0, 10) }); };

loadData(page === 'viewer' ? FAM : 'GH1').then(async () => {
  foot(); onLang(foot);
  if (page === 'viewer') {
    const title = () => { document.title = t('doc_title_viewer', { f: D.family }); };
    title(); onLang(title);
    const { initExplorer } = await import('./explorer.js');
    window.easimap = { D, explorer: initExplorer() };
  } else {
    await loadFamilies();
    const [{ initHero }, { initMetaphor }, { initFamily }, { initSkies }, { initIdea }, { initGuide }, { initStories }] = await Promise.all([
      import('./hero.js'), import('./metaphor.js'), import('./family.js'), import('./skies.js'), import('./idea.js'), import('./guide.js'), import('./stories.js'),
    ]);
    // Un fallo en una pieza no apaga las demás
    for (const init of [initHero, initMetaphor, initFamily, initSkies, initIdea, initGuide, initStories]) {
      try { init(); } catch (e) { console.error(init.name, e); }
    }
    window.easimap = { D };
    // El ancla llega antes que el alto real de las figuras
    if (location.hash) setTimeout(() => document.getElementById(location.hash.slice(1))?.scrollIntoView({ behavior: 'instant' }), 80);
  }
}).catch(e => {
  console.error(e);
  const c = document.getElementById('hero-caption') || document.getElementById('count-of');
  if (c) c.textContent = 'Data could not be loaded. No se pudieron cargar los datos.';
});
