"""Command line: uv run python -m foldlab <step> <experiment folder>"""

from __future__ import annotations

import argparse
import importlib

from .config import Experiment

STEPS = {
    "prepare": "fetch the sequence and the experimental reference, write one FASTA per variant",
    "msa": "search for evolutionary relatives, build one alignment per variant",
    "predict": "run AlphaFold2 on every variant with identical settings",
    "analyze": "measure confidence, structural differences, noise floor, conservation",
    "report": "draw the figures, the 3D viewer and REPORT.md",
}


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="foldlab",
        description="Steps: " + "; ".join(f"{k} = {v}" for k, v in STEPS.items()),
    )
    parser.add_argument("step", choices=[*STEPS, "all"])
    parser.add_argument("experiment", help="experiment folder, e.g. experiments/001-ubiquitin-L67A")
    args = parser.parse_args()

    exp = Experiment.load(args.experiment)
    for step in STEPS if args.step == "all" else [args.step]:
        print(f"\n=== {step}: {STEPS[step]} ===")
        importlib.import_module(f"foldlab.{step}").run(exp)


if __name__ == "__main__":
    main()
