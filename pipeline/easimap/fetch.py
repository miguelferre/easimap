"""Descargas con caché en disco: estructuras, HMM y mapeos de PDBe."""
from __future__ import annotations

from pathlib import Path
import gzip
import json
import time
import requests

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "cache"
RCSB_PDB = "https://files.rcsb.org/download/{code}.pdb"
RCSB_CIF = "https://files.rcsb.org/download/{code}.cif"
# Espejo del archivo wwPDB (mismo contenido). Se usa si files.rcsb.org devuelve error 5xx, como pasó el 2026-11-04.
WWPDB_PDB = "https://files.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb{lower}.ent.gz"
WWPDB_CIF = "https://files.wwpdb.org/pub/pdb/data/structures/all/mmCIF/{lower}.cif.gz"
INTERPRO_HMM = "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{acc}?annotation=hmm"
PDBE_MAP = "https://www.ebi.ac.uk/pdbe/api/mappings/{kind}/{code}"


def _get(url: str, retries: int = 6) -> requests.Response:
    for i in range(retries):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                r.raise_for_status()
        except requests.RequestException:
            if i == retries - 1:
                raise
        time.sleep(3 * (i + 1) ** 2)
    r.raise_for_status()
    return r


def _download(url: str, mirror: str) -> bytes:
    """Contenido de `url`; si falla, el del espejo (un .gz que se descomprime)."""
    try:
        return _get(url, retries=2).content
    except requests.RequestException:
        return gzip.decompress(_get(mirror).content)


def pdb_file(code: str, cache: Path = CACHE) -> Path:
    """Fichero PDB de RCSB (versión actual). Se descarga una vez."""
    dest = cache / "pdb" / f"{code.upper()}.pdb"
    if not dest.exists() or dest.stat().st_size == 0:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(_download(RCSB_PDB.format(code=code.upper()), WWPDB_PDB.format(lower=code.lower())))
    return dest


def cif_file(code: str, cache: Path = CACHE) -> Path:
    dest = cache / "cif" / f"{code.upper()}.cif"
    if not dest.exists() or dest.stat().st_size == 0:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(_download(RCSB_CIF.format(code=code.upper()), WWPDB_CIF.format(lower=code.lower())))
    return dest


def current_hmm(accession: str, cache: Path = CACHE) -> Path:
    """HMM vigente de Pfam, servido por InterPro."""
    dest = cache / "hmm" / f"{accession}.hmm"
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(gzip.decompress(_get(INTERPRO_HMM.format(acc=accession)).content).decode("latin-1"), encoding="latin-1")
    return dest


def pdbe_mapping(code: str, kind: str = "pfam", cache: Path = CACHE) -> dict:
    """Mapeo de PDBe (`pfam` o `interpro`) para una entrada: dominios con cadena y rango de residuos.

    Una entrada sin mapeo devuelve {}. Se guarda también, para no repetir la consulta.
    """
    dest = cache / "pdbe" / f"{code.lower()}.{kind}.json"
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        r = requests.get(PDBE_MAP.format(kind=kind, code=code.lower()), timeout=60)
        dest.write_text(r.text if r.status_code == 200 else "{}", encoding="utf-8")
    return json.loads(dest.read_text(encoding="utf-8"))
