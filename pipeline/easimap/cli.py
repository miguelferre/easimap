"""Línea de órdenes: `python -m easimap build GH1`."""
from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="easimap")
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="construye los datos de una familia con los métodos nuevos")
    b.add_argument("family")
    b.add_argument("--out", type=Path, help="salida; por defecto web/data/<familia>.json")
    b.add_argument("--limit", type=int, help="solo las primeras N entidades (pruebas)")
    r = sub.add_parser("replica", help="cifras de la réplica de 2020 (necesita la tabla de 2020)")
    r.add_argument("family")
    g = sub.add_parser("gap", help="prueba del desfase en las posiciones 334-336 de GH1")
    g.add_argument("family")
    args = ap.parse_args(argv)

    if args.cmd == "build":
        from .build import build
        out = args.out or ROOT / "web" / "data" / f"{args.family.lower()}.json"
        res = build(args.family, out=out, report_path=ROOT / "data" / f"{args.family.lower()}.build.json", limit=args.limit)
        for k, v in res["log"]["counts"].items():
            print(f"{k}: {v}")
        print("escrito", out)
    elif args.cmd == "gap":
        from .gap_check import run
        run(ROOT / "web" / "data" / f"{args.family.lower()}.json", ROOT / "data" / f"{args.family.lower()}.gap334.json")
        print("escrito data/%s.gap334.json" % args.family.lower())
    else:
        from .replica_report import report
        from .replica_report import default_out
        res = report(default_out(args.family), args.family)
        for k, v in res["variantes"].items():
            if "error" in v:
                print(k, "no se puede:", v["error"])
                continue
            print(k, v["proteinas_reproducidas"], "proteínas", v["pares"], "pares, mediana", f"{v['error_mediana_A']:.5f}", "Å")
    return 0
