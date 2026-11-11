"""Step 3 - predict: run AlphaFold2 (through ColabFold) on every variant with identical settings.

Two ways to run, giving the same predictions:

  one at a time   (default) one ColabFold process per variant. Right for a laptop.
  in parallel     FOLDLAB_JOBS=N: one process per variant and network, N of them at once.
                  AlphaFold2 on a Mac runs on the CPU and one process only keeps a few cores
                  busy, so a machine with many cores is used by running many processes.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

from .config import COLABFOLD_ENV, WEIGHTS_DIR, Experiment


def _sha256(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _command(exp: Experiment, variant: str, out: Path, models: list[int]) -> list[str]:
    settings = exp.cfg["prediction"]
    command = [
        "colabfold_batch", str(exp.msa / f"{variant}.a3m"), str(out),
        "--model-type", settings["model_type"],
        "--num-recycle", str(settings["num_recycle"]),
        "--num-seeds", str(settings["num_seeds"]),
        "--random-seed", str(settings["random_seed"]),
        "--data", str(WEIGHTS_DIR),
        # Only these networks, in this order. Seeds and everything else are unchanged.
        "--num-models", str(len(models)), "--model-order", ",".join(map(str, models)),
    ]
    hint = exp.templates(variant)
    if hint == "search":
        # ColabFold asks its server for the best-matching structures in the PDB.
        command += ["--templates"]
    elif hint:
        # Hand the network exactly these structures. ColabFold writes index files next to them,
        # so every process gets a private copy: ten processes sharing one folder would collide.
        folder = out / "_templates"
        folder.mkdir(parents=True, exist_ok=True)
        for pdb in hint:
            shutil.copy(exp.dir / "templates" / f"{pdb}.pdb", folder / f"{pdb.lower()}.pdb")
        command += ["--templates", "--custom-template-path", str(folder)]
    return command


def _reuse(exp: Experiment, variant: str) -> None:
    """Copy predictions an earlier experiment already made with identical inputs and settings."""
    reuse = exp.cfg.get("reuse")
    if not reuse or variant not in reuse.get("variants", []):
        return
    source = Experiment.load(exp.dir / reuse["from"])
    for model in exp.models:
        there = source.predictions / variant / f"model_{model}"
        here = exp.predictions / variant / f"model_{model}"
        if _done(there, variant) and not _done(here, variant):
            shutil.copytree(there, here, dirs_exist_ok=True)
            print(f"  {variant} model {model}: copied from {source.cfg['experiment']['id']}")


def _done(folder: Path, variant: str) -> bool:
    return (folder / f"{variant}.done.txt").exists()


def variant_done(exp: Experiment, variant: str) -> bool:
    """True if every prediction of this variant is on disk, in either layout."""
    out = exp.predictions / variant
    return _done(out, variant) or all(_done(out / f"model_{m}", variant) for m in exp.models)


def _run_task(exp: Experiment, variant: str, model: int) -> tuple[str, float]:
    out = exp.predictions / variant / f"model_{model}"
    label = f"{variant} model {model}"
    if _done(out, variant):
        return label, 0.0
    out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with open(out / "stdout.log", "w") as log:
        subprocess.run(_command(exp, variant, out, [model]), check=True, env=COLABFOLD_ENV,
                       stdout=log, stderr=subprocess.STDOUT)
    seconds = time.time() - started
    print(f"  {label}: done in {seconds:.0f} s", flush=True)
    return label, seconds


def run(exp: Experiment) -> None:
    settings = exp.cfg["prediction"]
    jobs = int(os.environ.get("FOLDLAB_JOBS", "1"))
    exp.predictions.mkdir(parents=True, exist_ok=True)
    exp.results.mkdir(parents=True, exist_ok=True)
    per_variant = len(exp.models) * settings["num_seeds"]

    timings: dict[str, float] = {}
    started = time.time()
    for name in exp.variants:
        _reuse(exp, name)
    todo = [name for name in exp.variants if not variant_done(exp, name)]
    for name in exp.variants:
        if name not in todo:
            print(f"{name}: already predicted, skipping")

    if jobs <= 1:
        for name in todo:
            print(f"{name}: predicting {per_variant} structures "
                  f"({len(exp.models)} models x {settings['num_seeds']} seeds)...", flush=True)
            t0 = time.time()
            subprocess.run(_command(exp, name, exp.predictions / name, exp.models), check=True, env=COLABFOLD_ENV)
            timings[name] = round(time.time() - t0, 1)
            print(f"{name}: done in {timings[name]:.0f} s")
    elif todo:
        tasks = [
            (name, model) for name in todo for model in exp.models
            if not _done(exp.predictions / name / f"model_{model}", name)
        ]
        print(f"Predicting {len(tasks) * settings['num_seeds']} structures as {len(tasks)} processes "
              f"(variant x network), {min(jobs, len(tasks))} at a time...", flush=True)
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            for label, seconds in pool.map(lambda task: _run_task(exp, *task), tasks):
                timings[label] = round(seconds, 1)
    wall = round(time.time() - started, 1)
    if todo:
        print(f"all predictions done in {wall:.0f} s")

    # Provenance: enough to rerun this exactly, or to explain why a rerun differs.
    provenance_file = exp.results / "provenance.json"
    previous = json.loads(provenance_file.read_text()) if provenance_file.exists() else {}
    runs = previous.get("prediction_runs", [])
    if todo:
        runs.append({
            "recorded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "variants": todo,
            "parallel_processes": jobs,
            "wall_seconds": wall,
            "seconds_per_process": timings,
            "machine": os.environ.get("FOLDLAB_MACHINE") or platform.node(),
        })
    provenance = {
        "experiment": exp.cfg["experiment"]["id"],
        "software": {
            "colabfold": version("colabfold"),
            "alphafold-colabfold": version("alphafold-colabfold"),
            "jax": version("jax"),
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "weights": "alphafold_params_2021-07-14 (AlphaFold2-ptm)",
        "settings": settings,
        "inputs_sha256": {
            f"msa/{name}.a3m": _sha256(exp.msa / f"{name}.a3m") for name in exp.variants
        },
        # What each variant was handed as a template, and for a search what ColabFold picked.
        "templates": {
            name: {"asked": exp.templates(name),
                   "used": sorted({p.stem for p in (exp.predictions / name).rglob("templates_*/*.cif")}
                                  | {p.stem for p in (exp.predictions / name).rglob("_templates/*.pdb")})}
            for name in exp.variants if exp.templates(name)
        },
        "prediction_runs": runs,
    }
    # Keep what an earlier version of this file recorded, so no history is lost.
    for key in ("recorded_utc", "seconds_per_variant"):
        if key in previous:
            provenance[f"first_run_{key}"] = previous[key]
        elif f"first_run_{key}" in previous:
            provenance[f"first_run_{key}"] = previous[f"first_run_{key}"]
    provenance_file.write_text(json.dumps(provenance, indent=2) + "\n")
    print("wrote results/provenance.json")
