// Mol* con carga diferida. Solo se descarga cuando alguien pide la estructura 3D.
import { t } from './core.js';

const VERSION = '5.13.1';
const BASE = `https://cdn.jsdelivr.net/npm/molstar@${VERSION}/build/viewer`;
let viewer = null, booting = null, loaded = null, last = null, msgKey = null, msgVars = null;
const box = () => document.getElementById('mol');
const say = (key, vars) => { msgKey = key; msgVars = vars; document.getElementById('mol-msg').textContent = key ? t(key, vars) : ''; };
export const molMessage = () => { if (msgKey) say(msgKey, msgVars); };
export const molIsOn = () => !!booting;

function loadAssets() {
  return new Promise((ok, no) => {
    const css = document.createElement('link');
    css.rel = 'stylesheet'; css.href = `${BASE}/theme/dark.css`; document.head.appendChild(css);
    const s = document.createElement('script');
    s.src = `${BASE}/molstar.js`; s.crossOrigin = 'anonymous'; s.onload = ok; s.onerror = () => no(new Error('molstar.js'));
    document.head.appendChild(s);
  });
}
function boot() {
  if (booting) return booting;
  box().classList.add('is-on');
  booting = loadAssets().then(() => window.molstar.Viewer.create(box(), {
    layoutIsExpanded: false, layoutShowControls: false, layoutShowSequence: false, layoutShowLog: false,
    layoutShowLeftPanel: false, viewportShowExpand: false, viewportShowSelectionMode: false,
    viewportShowAnimation: false, viewportShowTrajectoryControls: false, pdbProvider: 'rcsb',
    viewportBackgroundColor: '#0d1322', volumeStreamingDisabled: true,
  })).then(v => { viewer = v; return v; });
  return booting;
}

// Resalta el residuo y muestra su entorno en bolas y varillas
function focus(q) {
  const elements = { auth_asym_id: q.chain, auth_seq_id: q.auth };
  viewer.structureInteractivity({ elements, action: ['select'] });
  const plugin = viewer.plugin;
  const st = plugin.managers.structure.hierarchy.current.structures[0]?.cell?.obj?.data;
  if (st) {
    const loci = plugin.managers.structure.selection.getLoci(st);
    plugin.managers.structure.focus.setFromLoci(loci);
    plugin.managers.camera.focusLoci(loci, { minRadius: 11, extraRadius: 5, durationMs: 400 });
  }
}

export async function molShow(q) {
  last = q;
  if (q.auth == null) { say('mol_noid'); return; }
  say('mol_loading', { pdb: q.pdb });
  try {
    await boot();
    if (last !== q) return;
    if (loaded !== q.pdb) {
      await viewer.plugin.clear();
      await viewer.loadPdb(q.pdb.toLowerCase());
      loaded = q.pdb;
      if (last !== q) return;
    }
    focus(q);
    say('mol_ok', { res: q.residue, num: q.auth, pdb: q.pdb });
    window.__molState = { pdb: q.pdb, auth: q.auth, chain: q.chain };
  } catch (e) {
    console.warn('Mol*:', e);
    say('mol_fail');
  }
}
