"""Caché en disco cuya clave depende de todos los parámetros que cambian el resultado."""
from __future__ import annotations

import hashlib
import json
import pickle
from pathlib import Path
from typing import Callable

from . import fetch

CODE_VERSION = 3   # súbelo si cambia la lógica de lectura o de medida


def file_hash(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def key(**params) -> str:
    """Clave estable: hash del JSON ordenado de los parámetros y de la versión del código."""
    blob = json.dumps({"code_version": CODE_VERSION, **params}, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:20]


def memo(kind: str, params: dict, compute: Callable, root: Path = fetch.CACHE) -> object:
    """Devuelve el resultado guardado para estos parámetros o lo calcula y lo guarda."""
    dest = root / "build" / kind / f"{key(**params)}.pkl"
    if dest.exists():
        return pickle.loads(dest.read_bytes())
    value = compute()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(pickle.dumps(value))
    return value
