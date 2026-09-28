"""Ejecuta el pipeline completo.

Uso:
    python run_pipeline.py config.yaml
    python run_pipeline.py config.yaml --noches noche1
    python run_pipeline.py --demo           # genera datos sintéticos y los procesa
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from instrastro import pipeline, synthetic  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", nargs="?", help="archivo YAML de configuración")
    parser.add_argument("--noches", nargs="*", help="procesar solo estas noches (por nombre)")
    parser.add_argument("--demo", metavar="CARPETA", nargs="?", const="demo",
                        help="genera un set sintético en CARPETA y lo procesa")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")

    if args.demo:
        _, config = synthetic.make_dataset(args.demo)
    elif args.config:
        config = args.config
    else:
        parser.error("indica un archivo de configuración o usa --demo")
    out = pipeline.run(config, noches=args.noches)
    print(f"Resultados en {out}")


if __name__ == "__main__":
    main()
