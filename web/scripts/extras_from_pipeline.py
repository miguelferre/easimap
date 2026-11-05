#!/usr/bin/env python
"""Genera los extras de la web desde los ficheros del pipeline.

web/data/gh1.extras.json sale de gh1.json, el conjunto actual, con historia y desfases.
web/data/gh5, gh16, gh18 y gh20.extras.json llevan nombres, organismos, estadisticas y bandas.
web/data/families.json resume las cinco familias para el inicio, con un cielo pequeno de cada una.

Nombres y organismos salen de data/gh1.build.json y de RCSB (solo con codigos PDB).
Estadisticas, bandas de RMSF, historia de la 333 y desfases usan metodos.py. El ligando y
los puntos del 3SCU se miden en el PDB, porque el 3SCU da el marco del sustrato y no esta
en el conjunto actual (D-012). Nada sale de la tabla de 2020 (D-023).

Uso: python web/scripts/extras_from_pipeline.py
El resultado no depende de la fecha: la consulta a RCSB se guarda en una carpeta temporal.
"""
import base64
import collections
import json
import os
import sys
import tomllib
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import metodos as seed  # noqa: E402
from easimap import fetch as pipe_fetch, replica2020  # noqa: E402
from easimap.geometry import to_frame  # noqa: E402
from easimap.pfam_archive import archived_hmm  # noqa: E402
from easimap.structure import read_ligand_atoms  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "web" / "data"
OUT = Path(os.environ.get("EASIMAP_OUT", str(DATA))) / "gh1.extras.json"
GRAPHQL = "https://data.rcsb.org/graphql"
SKIES = ["3SCU", "1CBG", "2E9L", "1WCG", "3WQ4", "3AHZ", "1GOW", "3W53"]
SKY_MATES = {"3SCU": "2RGM", "1GOW": "7UZ1", "3AHZ": "3VII", "3WQ4": "3WQ5"}   # otra entrada de su grupo al 90 %
NOT_IN_TABLE = {"1DWA", "6JFP"}   # en la lista de 2020 de GH1.toml, sin filas en su tabla
CHIMERA = "4GXP"
HMM_MATCH = dict(same=27521, total=27820)   # numeracion de Pfam 32.0 frente a la tabla
QUERY = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id polymer_entities{
  rcsb_polymer_entity{pdbx_description}
  rcsb_polymer_entity_container_identifiers{auth_asym_ids}
  rcsb_polymer_entity_group_membership{aggregation_method group_id similarity_cutoff}
  rcsb_entity_source_organism{ncbi_scientific_name}}}}"""


def load(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def rcsb(codes):
    """Metadatos por entrada. Se guardan por conjunto de codigos para repetir la salida."""
    path = seed.CACHE / ("extras_rcsb_" + "_".join(sorted(codes))[:60] + f"_{len(codes)}.json")
    seed.CACHE.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        body = json.dumps({"query": QUERY, "variables": {"ids": codes}}).encode()
        seed.fetch(GRAPHQL, path, data=body, headers={"Content-Type": "application/json"})
    return {e["rcsb_id"]: e for e in load(path)["data"]["entries"]}


def entity_of(entry, chain):
    ents = entry["polymer_entities"]
    for e in ents:
        if chain in (e["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"] or []):
            return e
    return ents[0]


def group90(entity):
    for g in entity["rcsb_polymer_entity_group_membership"] or []:
        if g["aggregation_method"] == "sequence_identity" and g["similarity_cutoff"] == 90.0:
            return g["group_id"]
    return None


def organisms(entity):
    names = []
    for o in entity["rcsb_entity_source_organism"] or []:
        n = o.get("ncbi_scientific_name")
        if n and n not in names:
            names.append(n)
    return names


def rows_by_protein(f):
    out = collections.OrderedDict()
    for i, c in enumerate(f["pdb"]):
        out.setdefault(c, []).append(i)
    return out


def rmsf_table(f, idx):
    """RMSF por posicion alrededor del centroide, con 20 proteinas o mas, como el seed."""
    X = np.array([f["x"], f["y"], f["z"]]).T
    res = {}
    for p, rows in idx.items():
        if len(rows) >= 20:
            Y = X[rows]
            res[p] = float(np.sqrt(((Y - Y.mean(0)) ** 2).sum(1).mean()))
    return res


def reference_points(cfg):
    """Puntos y ligando de la estructura de referencia, medidos en el PDB con el marco del pipeline."""
    code = cfg["reference_pdb"]
    path = pipe_fetch.pdb_file(code)
    tables, _, _ = replica2020.run(cfg, path.parent, archived_hmm(cfg["pfam"], cfg["hmm_release"]), [code], local=True)
    q = tables[code]
    order = np.argsort(q.resid)   # la tabla sigue el orden de texto del pipeline: 10, 100, 101
    xyz = np.array(q.xyz)[order]
    pick = lambda a, n=None: [round(float(a[i]), n) if n else a[i] for i in order]  # noqa: E731
    points = dict(x=[round(float(v), 2) for v in xyz[:, 0]], y=[round(float(v), 2) for v in xyz[:, 1]], z=[round(float(v), 2) for v in xyz[:, 2]],
                  residue=pick(q.resname), hmm_position=[p or None for p in pick(q.hmm_position)], auth_seq_id=pick(q.resid),
                  distance=pick(q.distance, 2), theta=pick(q.theta, 4), phi=pick(q.phi, 4))
    origin, axes = replica2020.substrate_from_reference(path, cfg)
    ligand = []
    for seq in range(1, 6):   # la celopentaosa: BGC 1 a 5 de la cadena C
        for name, c in read_ligand_atoms(path, "BGC", "C", seq).items():
            x, y, z = to_frame(c, origin, axes)[0]
            ligand.append(dict(element=name[0], name=f"{name}.{seq}", x=round(float(x), 2), y=round(float(y), 2), z=round(float(z), 2)))
    return {code: points}, ligand


def main_gh1():
    G = load(DATA / "gh1.json")
    CFG = tomllib.loads((ROOT / "familias" / "GH1.toml").read_text(encoding="utf-8"))
    B = load(ROOT / "data" / "gh1.build.json")
    R = load(ROOT / "data" / "replica2020" / "resumen.json")
    GAP = load(ROOT / "data" / "gh1.gap334.json")
    cat, core = G["fits"]["catalytic"], G["fits"]["core"]
    entries = {e["pdb"]: e for e in B["entries"]}
    codes = list(rows_by_protein(cat))
    assert set(codes) == set(entries), "gh1.json y gh1.build.json no tienen las mismas estructuras"
    old = set(CFG["replica2020"]["proteins"]) - NOT_IN_TABLE   # solo para saber qué estructuras son nuevas

    # Nombres, organismos y grupos al 90 %
    wanted = sorted(set(codes) | set(SKIES))
    meta = {}
    for i in range(0, len(wanted), 40):
        meta.update(rcsb(wanted[i:i + 40]))
    idx_rows = rows_by_protein(cat)
    proteins, groups = {}, {}
    for c in codes:
        e = entries[c]
        ent = entity_of(meta[c], e["chain"])
        orgs = organisms(ent)
        proteins[c] = dict(domain=e["domain"], organism=" / ".join(orgs) or e["organism"], chain=e["chain"],
                           n_residues=len(idx_rows[c]), resolution=e["resolution"],
                           name=(ent["rcsb_polymer_entity"] or {}).get("pdbx_description") or "")
        groups.setdefault(group90(ent), []).append(c)

    # Los ocho cielos: la propia entrada u otra de su grupo al 90 %
    skies = {}
    for c in SKIES:
        pdb = c if c in entries else SKY_MATES[c]
        assert pdb in entries, pdb
        skies[c] = dict(pdb=pdb, same_entry=pdb == c)

    reference, ligand = reference_points(CFG)

    # Por posicion
    pos_idx = collections.defaultdict(list)
    for i, p in enumerate(cat["hmm_position"]):
        if p:
            pos_idx[int(p)].append(i)
    base, new = rmsf_table(cat, pos_idx), rmsf_table(core, pos_idx)
    X = np.array([cat["x"], cat["y"], cat["z"]]).T
    cdist = {p: float(np.linalg.norm(X[pos_idx[p]].mean(0))) for p in base}
    bands = []
    for lo, hi, lab in ((0, 8, "0-8"), (8, 12, "8-12"), (12, 20, "12-20"), (20, 30, "20-30"), (30, 1e9, "30+")):
        ps = [p for p in base if lo <= cdist[p] < hi]
        if ps:
            bands.append(dict(band=lab, n=len(ps), rmsf_catalytic=round(float(np.median([base[p] for p in ps])), 2),
                              rmsf_core=round(float(np.median([new[p] for p in ps])), 2)))
    seed.log("bandas", bands)

    # Historia de la 333, con los metodos del seed
    dom = {c: entries[c]["domain"] for c in codes}
    by_pos = collections.defaultdict(dict)
    for p, rows in pos_idx.items():
        for i in rows:
            by_pos[p][cat["pdb"][i]] = dict(xyz=X[i], res=cat["residue"][i], th=cat["theta"][i], ph=cat["phi"][i], d=cat["distance"][i])
    pts = [dict(pdb=m, th=float(r["th"]), ph=float(r["ph"]), d=float(r["d"]), xyz=r["xyz"], organism=proteins[m]["organism"], domain=dom[m])
           for m, r in sorted(by_pos[333].items()) if r["res"] == "TRP"]
    s333, lab = seed.story(pts, "conjunto actual")
    s333["points"] = [dict(pdb=p["pdb"], theta=round(p["th"], 4), phi=round(p["ph"], 4), distance=round(p["d"], 2), domain=p["domain"],
                           group="A" if g else "B", organism=p["organism"]) for p, g in zip(pts, lab)]

    def fits(points, labels):
        pairs = [(p, g) for p, g in zip(points, labels) if p["domain"] != "Unknown"]
        return dict(known=len(pairs), in_group=sum((p["domain"] == "Eukaryota") == bool(g) for p, g in pairs))

    s333["fit_all"] = fits(pts, lab)
    newp = [p for p in pts if p["pdb"] not in old]
    sn, ln = seed.story(newp, "solo nuevas")
    s333["new_only"] = dict(n_trp=sn["n_trp"], by_domain=sn["by_domain"], groups=sn["groups"], centroid_separation=sn["centroid_separation"],
                            agreement=sn["agreement"], n_known=sn["n_known"], p_perm=sn["p_perm"], **fits(newp, ln))
    th = {d: [p["th"] for p in pts if p["domain"] == d] for d in ("Eukaryota",)}
    th["Prokaryota"] = [p["th"] for p in pts if p["domain"] in ("Bacteria", "Archaea")]
    s333["theta_median_deg"] = {k: round(float(np.degrees(np.median(v))), 1) for k, v in th.items()}
    ch = [(p, g) for p, g in zip(pts, lab) if p["pdb"] == CHIMERA]
    s333["chimera"] = dict(pdb=CHIMERA, group="A" if ch[0][1] else "B", organism=proteins[CHIMERA]["organism"]) if ch else None
    pos = {q["hmm_position"]: q for q in G["positions"]}
    s333["rmsf"] = dict(current=pos[333]["rmsf"], n=pos[333]["n"])

    off = seed.offsets(by_pos, dom, "conjunto actual")
    off["pairs"] = [f'{q["eukaryote"]}/{q["prokaryote"]}' for q in GAP["pairs"]]

    # Replica y seleccion, tal como las dejo el pipeline
    counts = B["counts"]
    rv = R["variantes"]["pfam32.0_local"]
    gen = G["generated"]
    extras = dict(
        family="GH1", generated=gen, proteins=proteins,
        stats=dict(rows_out=len(cat["pdb"]), proteins=len(codes), positions=len(G["positions"]), core_size=len(B["core_positions"]),
                   by_domain=dict(collections.Counter(dom.values())), entities_annotated=next(v for k, v in counts.items() if "IPR001360" in k),
                   hmm=B["hmm"], generated=gen),
        rmsf_bands=dict(bands=bands, core_size=len(B["core_positions"])),
        story333=s333, offsets=off, skies=skies, reference_structures=reference, ligand=ligand,
        replica=dict(pairs=rv["pares"], proteins=rv["proteinas_reproducidas"], median_error=rv["error_mediana_A"],
                     not_reproduced=rv["en_tabla_2020_y_no_reproducidas"], glocal_median_error=R["variantes"]["actual_glocal"]["error_mediana_A"],
                     hmm_match=HMM_MATCH))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(extras, separators=(",", ":"), ensure_ascii=False, sort_keys=True, allow_nan=False), encoding="utf-8")
    seed.log(f"escrito {OUT} ({OUT.stat().st_size / 1e3:.0f} KB)")
    for k in ("n_trp", "by_domain", "centroid_separation", "agreement", "n_known", "p_perm", "fit_all", "theta_median_deg", "chimera", "rmsf"):
        seed.log(k, s333[k])
    seed.log("grupos", s333["groups"])
    seed.log("nuevas", {k: v for k, v in s333["new_only"].items()})
    seed.log("estadisticas", extras["stats"])
    seed.log("cielos", skies)
    seed.log("desfases", off["table"], off["hits"], off["tested"], off["n_euk"], off["n_pro"], off["pairs"])


def compact_replica(stem):
    path = ROOT / "data" / "replica2020" / f"resumen-{stem}.json"
    if not path.exists():
        return None
    r = load(path)
    v = r["variantes"]
    k = "pfam32.0_local" if v.get("pfam32.0_local", {}).get("pares") else "tabla_2020"
    b = v[k]
    return dict(variant=k, proteins_2020=r["proteinas_2020"], proteins=b["proteinas_reproducidas"], pairs=b["pares"],
                median_error=b["error_mediana_A"], not_reproduced=sorted(b.get("en_tabla_2020_y_no_reproducidas") or []))


def family_extras(stem):
    """Extras de una familia del pipeline: nombres, organismos, estadisticas y bandas de RMSF."""
    G = load(DATA / f"{stem}.json")
    B = load(ROOT / "data" / f"{stem}.build.json")
    cat, core = G["fits"]["catalytic"], G["fits"]["core"]
    entries = {e["pdb"]: e for e in B["entries"]}
    idx_rows = rows_by_protein(cat)
    codes = list(idx_rows)
    assert set(codes) == set(entries), stem
    meta = {}
    for i in range(0, len(codes), 40):
        meta.update(rcsb(sorted(codes)[i:i + 40]))
    proteins = {}
    for c in codes:
        e = entries[c]
        ent = entity_of(meta[c], e["chain"])
        proteins[c] = dict(domain=e["domain"], organism=" / ".join(organisms(ent)) or e["organism"], chain=e["chain"], n_residues=len(idx_rows[c]),
                           resolution=e["resolution"], name=(ent["rcsb_polymer_entity"] or {}).get("pdbx_description") or "")
    pos_idx = collections.defaultdict(list)
    for i, p in enumerate(cat["hmm_position"]):
        if p:
            pos_idx[int(p)].append(i)
    base, new = rmsf_table(cat, pos_idx), rmsf_table(core, pos_idx)
    X = np.array([cat["x"], cat["y"], cat["z"]]).T
    cdist = {p: float(np.linalg.norm(X[pos_idx[p]].mean(0))) for p in base}
    bands = []
    for lo, hi, lab in ((0, 8, "0-8"), (8, 12, "8-12"), (12, 20, "12-20"), (20, 30, "20-30"), (30, 1e9, "30+")):
        ps = [p for p in base if lo <= cdist[p] < hi]
        if ps:
            bands.append(dict(band=lab, n=len(ps), rmsf_catalytic=round(float(np.median([base[p] for p in ps])), 2),
                              rmsf_core=round(float(np.median([new[p] for p in ps])), 2)))
    core_size = len(B["core_positions"])
    extras = dict(family=G["family"], generated=G["generated"], proteins=proteins,
                  stats=dict(rows_out=len(cat["pdb"]), proteins=len(codes), positions=len(G["positions"]), core_size=core_size,
                             by_domain=dict(collections.Counter(e["domain"] for e in B["entries"])), hmm=B["hmm"], generated=G["generated"]),
                  rmsf_bands=dict(bands=bands, core_size=core_size), replica=compact_replica(stem))
    (OUT.parent / f"{stem}.extras.json").write_text(json.dumps(extras, separators=(",", ":"), ensure_ascii=False, sort_keys=True, allow_nan=False), encoding="utf-8")
    seed.log(f"{stem}: {len(codes)} estructuras, bandas {bands}")


def sky_bytes(cat):
    """Cielo norte de la familia en 3 bytes por residuo: u, v y clase con magnitud. Residuos a menos de 18 A."""
    cls = dict(ASP=0, GLU=0, ARG=1, HIS=1, LYS=1, ASN=2, CYS=2, GLN=2, SER=2, THR=2, PHE=3, TRP=3, TYR=3, GLY=4, PRO=4)
    alias = dict(MSE="MET", CAS="CYS")
    lim = (6, 9, 13, 18, 25)
    out = bytearray()
    for th, ph, d, r in zip(cat["theta"], cat["phi"], cat["distance"], cat["residue"]):
        if th >= np.pi / 2 or d >= 18:
            continue
        rad = np.sqrt(2) * np.sin(th / 2)
        u, v = -rad * np.sin(ph), -rad * np.cos(ph)
        m = next((i for i, x in enumerate(lim) if d < x), 5)
        out += bytes((int(round((u + 1) * 127.5)), int(round((v + 1) * 127.5)), cls.get(alias.get(r, r), 5) * 8 + m))
    return base64.b64encode(bytes(out)).decode()


def families():
    """Resumen de las cinco familias para la lamina del inicio. Todas las cifras salen de estos ficheros."""
    out = []
    for key, stem in (("GH1", "gh1"), ("GH5", "gh5"), ("GH16", "gh16"), ("GH18", "gh18"), ("GH20", "gh20")):
        G, E = load(DATA / f"{stem}.json"), load(DATA / f"{stem}.extras.json")
        cat = G["fits"]["catalytic"]
        cand = None
        cpath = ROOT / "data" / f"{stem}.candidatas.json"
        if cpath.exists():
            c = load(cpath)
            ps = c["posiciones"]
            cand = dict(tested=c["resumen"]["posiciones_probadas"], q05=sum(p["q_bh"] < .05 for p in ps), family05=sorted(p["hmm_position"] for p in ps if p["p_familia"] < .05),
                        euk=c["resumen"]["estructuras_a"], pro=c["resumen"]["estructuras_b"])
        out.append(dict(key=key, file=stem, source=G["source"], generated=G["generated"], reference_pdb=G["reference_pdb"], ligand=G["ligand"],
                        reference_positions=G["reference_positions"], structures=E["stats"]["proteins"], residues=len(cat["pdb"]), positions=len(G["positions"]),
                        by_domain=E["stats"]["by_domain"], core_size=E["stats"].get("core_size", 0),
                        median_rmsf=round(float(np.median([q["rmsf"] for q in G["positions"]])), 2), candidates=cand, sky=sky_bytes(cat)))
        seed.log(key, {k: v for k, v in out[-1].items() if k not in ("sky",)}, "cielo", len(out[-1]["sky"]))
    (OUT.parent / "families.json").write_text(json.dumps(out, separators=(",", ":"), ensure_ascii=False, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    main_gh1()
    for fam in ("gh5", "gh16", "gh18", "gh20"):
        family_extras(fam)
    families()
