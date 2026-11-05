#!/usr/bin/env python
"""Pone un solo valor de version en index.html y visor.html.

El valor es un resumen del contenido de web/js, web/css y web/data. Los modulos, el estilo y
los datos se piden con ?v=<valor>, y asi un navegador no mezcla ficheros de dos subidas.
Hay que ejecutarlo antes de cada commit que cambie esos ficheros.

Uso: python web/scripts/stamp_version.py [--check]
"""
import hashlib
import re
import sys
from pathlib import Path

WEB = Path(__file__).resolve().parents[1]


def version():
    h = hashlib.sha1()
    for folder in ("js", "css", "data"):
        for f in sorted((WEB / folder).rglob("*")):
            if f.is_file():
                h.update(f.relative_to(WEB).as_posix().encode())
                h.update(f.read_bytes().replace(b"\r\n", b"\n"))
    return h.hexdigest()[:10]


def main():
    v = version()
    stale = []
    for name in ("index.html", "visor.html"):
        p = WEB / name
        s = p.read_text(encoding="utf-8")
        t = re.sub(r"window\.EASIMAP_V = '[0-9a-z]+'", f"window.EASIMAP_V = '{v}'", s)
        t = re.sub(r'\?v=[0-9a-z]+"', f'?v={v}"', t)
        if t != s:
            stale.append(name)
            if "--check" not in sys.argv:
                p.write_text(t, encoding="utf-8", newline="\n")
    print(v, "sin cambios" if not stale else ("desfasado: " if "--check" in sys.argv else "escrito: ") + ", ".join(stale))
    return 1 if stale and "--check" in sys.argv else 0


if __name__ == "__main__":
    sys.exit(main())
