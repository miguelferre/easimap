"""Métodos comunes de los extras: consulta con caché, grupos de la 333 y desfases."""
import collections
import os
import tempfile
import urllib.request
from pathlib import Path

import numpy as np

CACHE = Path(os.environ.get("EASIMAP_CACHE", str(Path(tempfile.gettempdir()) / "easimap-seed-cache")))
UA = {"User-Agent": "easimap-seed"}
RNG_SEED, N_PERM = 1, 20000


def log(*a):
    print(*a, flush=True)


def fetch(url, path, data=None, headers=None):
    if path.exists():
        return path.read_bytes()
    req = urllib.request.Request(url, data=data, headers=headers or UA)
    body = urllib.request.urlopen(req, timeout=90).read()
    path.write_bytes(body)
    return body


def two_means(X, th):
    c = X[[int(th.argmin()), int(th.argmax())]]
    l = None
    for _ in range(50):
        l = ((X[:, None] - c[None]) ** 2).sum(2).argmin(1)
        c = np.array([X[l == k].mean(0) for k in (0, 1)])
    return l == 0, c


def story(pts, label):
    """pts: lista de dicts pdb, th, ph, d, domain, xyz, organism. Como check5.py."""
    th = np.array([p["th"] for p in pts])
    X = np.array([p["xyz"] for p in pts])
    dom = np.array([p["domain"] for p in pts])
    lab, c = two_means(X, th)
    known = dom != "Unknown"
    euk = dom == "Eukaryota"
    obs = (lab[known] == euk[known]).mean()
    obs = float(max(obs, 1 - obs))
    rng = np.random.default_rng(RNG_SEED)
    perm = np.array([max(q, 1 - q) for q in ((lab[known] == rng.permutation(euk[known])).mean() for _ in range(N_PERM))])
    cnt = {g: collections.Counter(dom[lab == (g == "A")].tolist()) for g in ("A", "B")}
    res = dict(
        n_trp=len(pts), by_domain=dict(collections.Counter(dom.tolist())),
        groups={g: dict(n=int((lab == (g == "A")).sum()), by_domain=dict(cnt[g]),
                        theta_mean_deg=round(float(np.degrees(th[lab == (g == "A")].mean())), 1)) for g in ("A", "B")},
        centroid_separation=round(float(np.linalg.norm(c[0] - c[1])), 2),
        agreement=round(obs, 3), n_known=int(known.sum()),
        p_perm=round(float((perm >= obs).mean()), 5), n_perm=N_PERM, seed=RNG_SEED)
    log(f"  [{label}] TRP {len(pts)} {res['by_domain']} A={res['groups']['A']['n']} {cnt['A']} th={res['groups']['A']['theta_mean_deg']} "
        f"B={res['groups']['B']['n']} {cnt['B']} th={res['groups']['B']['theta_mean_deg']} sep={res['centroid_separation']} "
        f"acuerdo={res['agreement']} n={res['n_known']} p={res['p_perm']}")
    return res, lab


def offsets(by_pos, domof, label):
    euks = {m for m, d in domof.items() if d == "Eukaryota"}
    pro = {m for m, d in domof.items() if d in ("Bacteria", "Archaea")}

    def sep(p, k, mn):
        aa = [r["xyz"] for m, r in by_pos.get(p, {}).items() if m in euks]
        bb = [r["xyz"] for m, r in by_pos.get(p + k, {}).items() if m in pro]
        if len(aa) >= mn and len(bb) >= mn:
            return float(np.linalg.norm(np.median(aa, 0) - np.median(bb, 0)))
        return None

    table = []
    for p in range(333, 338):
        a, b = sep(p, 0, 5), sep(p, -1, 5)
        table.append(dict(p=p, k0=None if a is None else round(a, 1), km1=None if b is None else round(b, 1)))
    hits, tested = [], 0
    for p in sorted(by_pos):
        res = [(sep(p, k, 10), k) for k in (-1, 0, 1)]
        if all(r[0] is not None for r in res):
            tested += 1
            best = min(res)
            zero = [q for q in res if q[1] == 0][0][0]
            if best[1] != 0 and zero > 4 and best[0] < 3:
                hits.append([p, best[1], round(best[0], 1), round(zero, 1)])
    log(f"  [{label}] offsets {[(t['p'], t['k0'], t['km1']) for t in table]}")
    log(f"  [{label}] hits {len(hits)} de {tested}: {[h[0] for h in hits]}")
    return dict(table=table, hits=hits, tested=tested,
                n_euk=len(euks), n_pro=len(pro))
