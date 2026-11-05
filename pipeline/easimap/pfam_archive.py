"""Extrae un solo HMM de un lanzamiento antiguo de Pfam, sin guardar el fichero completo.

Pfam-A.hmm.gz pesa unos 270 MB. Se lee en flujo y se para al terminar el modelo pedido.
"""
from __future__ import annotations

from pathlib import Path
import time
import zlib
import requests

URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{release}/Pfam-A.hmm.gz"


def extract_hmms(wanted: dict[str, Path], release: str, chunk: int = 1 << 20, retries: int = 4) -> dict[str, Path]:
    """Descarga varios HMM (`{accession sin versión: destino}`) de Pfam `release` en una sola lectura del archivo.

    Se para en cuanto aparecen todos. Si la conexión falla, reintenta desde el principio.
    """
    for dest in wanted.values():
        dest.parent.mkdir(parents=True, exist_ok=True)
    last: Exception | None = None
    for attempt in range(retries):
        try:
            return _stream(dict(wanted), release, chunk)
        except (requests.RequestException, zlib.error) as exc:   # el FTP de EBI corta a veces el TLS
            last = exc
            time.sleep(5 * (attempt + 1))
    raise last  # type: ignore[misc]


def _stream(pending: dict[str, Path], release: str, chunk: int) -> dict[str, Path]:
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    buf = ""
    block: list[str] = []
    acc = ""
    done: dict[str, Path] = {}
    with requests.get(URL.format(release=release), stream=True, timeout=60) as r:
        r.raise_for_status()
        for raw in r.iter_content(chunk):
            buf += dec.decompress(raw).decode("latin-1")
            lines = buf.split("\n")
            buf = lines.pop()
            for line in lines:
                if line.startswith("HMMER3"):
                    block, acc = [], ""
                block.append(line)
                if line.startswith("ACC "):
                    acc = line.split()[1].split(".")[0]
                if line == "//":
                    if acc in pending:
                        pending[acc].write_text("\n".join(block) + "\n", encoding="latin-1")
                        done[acc] = pending.pop(acc)
                        if not pending:
                            return done
                    block = []
    raise LookupError(f"{sorted(pending)} no está en Pfam {release}")


def extract_hmm(accession: str, release: str, dest: Path, chunk: int = 1 << 20) -> Path:
    """Descarga el HMM `accession` (sin versión, p. ej. PF00232) de Pfam `release` (p. ej. 32.0) a `dest`."""
    return extract_hmms({accession: dest}, release, chunk)[accession]


def archived_hmm(accession: str, release: str, cache: Path | None = None) -> Path:
    """HMM de un lanzamiento antiguo de Pfam, en caché. La primera vez tarda unos cinco minutos (lee 270 MB en flujo)."""
    from . import fetch
    dest = (cache or fetch.CACHE) / "hmm" / f"{accession}_pfam{release}.hmm"
    if not dest.exists():
        extract_hmm(accession, release, dest)
    return dest
