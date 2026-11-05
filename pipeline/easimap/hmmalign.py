"""Alineamiento de una secuencia a un HMM de Pfam y número de posición por residuo."""
from __future__ import annotations

from pathlib import Path
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence

_ALPHABET = Alphabet.amino()
_cache: dict[str, pyhmmer.plan7.HMM] = {}


def load_hmm(path: str | Path):
    key = str(path)
    if key not in _cache:
        with pyhmmer.plan7.HMMFile(key) as fh:
            _cache[key] = next(iter(fh))
    return _cache[key]


def assign_positions(hmm, seq: str, all_consensus_cols: bool = True) -> list[int]:
    """Posición HMM de cada residuo de `seq`; 0 si el residuo cae fuera de las columnas de consenso.

    all_consensus_cols=True da el número real de la columna del modelo (1..M).
    all_consensus_cols=False reproduce `hmmalign` de la línea de órdenes sin `--allcol`: las columnas
    de consenso sin ningún residuo desaparecen y la cuenta salta. Así numeraba 2020.
    """
    dseq = TextSequence(name=b"q", sequence=seq).digitize(_ALPHABET)
    msa = pyhmmer.hmmer.hmmalign(hmm, [dseq], all_consensus_cols=all_consensus_cols)
    aligned = msa.sequences[0].sequence
    ref = msa.reference
    if ref is None:
        raise RuntimeError("hmmalign no devolvió la línea RF")
    if isinstance(ref, bytes):
        ref = ref.decode()
    out: list[int] = []
    p = 0
    for ch, rf in zip(aligned, ref):
        if rf == "x":
            p += 1
        if ch not in "-.":
            out.append(p if rf == "x" else 0)
    if len(out) != len(seq):
        raise RuntimeError(f"alineamiento incoherente: {len(out)} residuos frente a {len(seq)}")
    return out


def assign_positions_local(hmm, seq: str) -> list[int]:
    """Posición HMM de cada residuo con alineamiento local de un solo dominio (el de mayor puntuación).

    Así numeraba 2020: con el HMM de Pfam 32.0 coincide en el 98,9 % de los residuos.
    Un residuo en columna de inserción, o fuera del dominio, recibe 0.
    """
    dseq = TextSequence(name=b"q", sequence=seq).digitize(_ALPHABET)
    best = None
    for hits in pyhmmer.hmmsearch(hmm, [dseq], E=1e6, Z=1, domE=1e6, incE=1e6, incdomE=1e6):
        for hit in hits:
            for dom in hit.domains:
                if best is None or dom.score > best.score:
                    best = dom
    out = [0] * len(seq)
    if best is None:
        return out
    al = best.alignment
    p, t = al.hmm_from, al.target_from - 1
    for hc, tc in zip(al.hmm_sequence, al.target_sequence):
        if hc == ".":          # inserción: residuo sin columna
            t += 1
        elif tc == "-":        # deleción: avanza el modelo
            p += 1
        else:                  # columna de consenso con residuo
            out[t] = p
            t += 1
            p += 1
    return out
