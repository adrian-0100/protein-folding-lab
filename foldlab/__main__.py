"""Command line: uv run python -m foldlab <step> <experiment folder>"""

from __future__ import annotations

import argparse
import importlib

from .config import Experiment

# step -> (what it does, the experiment.toml section that switches it on; None = always)
STEPS = {
    "select": ("choose the mutants by a written rule, from experimental structures only", "selection"),
    "prepare": ("fetch the sequence, the experimental structures and data, write one FASTA per variant", None),
    "msa": ("search for evolutionary relatives, build one alignment per variant", None),
    "score": ("score every measured mutant with sequence-based yardsticks", "scores"),
    "predict": ("run AlphaFold2 on every variant with identical settings", None),
    "analyze": ("measure confidence, structural differences, noise floor, agreement with experiment", None),
    "report": ("draw the figures, the 3D viewer and REPORT.md", None),
}
# `ready` stops before anything is computed: everything an experiment needs is fetched,
# chosen and checked, and nothing is predicted or scored.
GROUPS = {"all": list(STEPS), "ready": ["select", "prepare", "msa"]}


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="foldlab",
        description="Steps: " + "; ".join(f"{k} = {v[0]}" for k, v in STEPS.items())
        + ". 'ready' = select, prepare, msa. 'all' = every step that applies.",
    )
    parser.add_argument("step", choices=[*STEPS, *GROUPS])
    parser.add_argument("experiment", help="experiment folder, e.g. experiments/001-ubiquitin-L67A")
    args = parser.parse_args()

    for step in GROUPS.get(args.step, [args.step]):
        exp = Experiment.load(args.experiment)  # reloaded: a step may have added variants
        what, section = STEPS[step]
        if section and section not in exp.cfg:
            if args.step == step:
                raise SystemExit(f"'{step}' needs a [{section}] section in experiment.toml")
            continue
        print(f"\n=== {step}: {what} ===", flush=True)
        importlib.import_module(f"foldlab.{step}").run(exp)


if __name__ == "__main__":
    main()
