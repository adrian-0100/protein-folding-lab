"""Laboratory measurements to compare scores against: a deep mutational scan from ProteinGym.

A deep mutational scan measures, in one experiment, what each of a protein's possible single
mutations does (here: how well yeast grows when its only ubiquitin carries that mutation).
ProteinGym is a public collection of such scans, cleaned and labelled.
"""

from __future__ import annotations

import json
import random
import zipfile

import pandas as pd
import requests

from .config import DATA_DIR, MUTATION, Experiment

ARCHIVE = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip"


def truth_dir(exp: Experiment):
    return exp.dir / "truth"


def fetch(exp: Experiment) -> None:
    """Download the ProteinGym archive once per machine and copy this experiment's scans out of it."""
    wanted = [exp.cfg["truth"]["assay"], *filter(None, [exp.cfg["truth"].get("ceiling_assay")])]
    out = truth_dir(exp)
    out.mkdir(exist_ok=True)
    if all((out / f"{assay}.csv").exists() for assay in wanted):
        return
    archive = DATA_DIR / "proteingym" / "DMS_ProteinGym_substitutions.zip"
    if not archive.exists():
        archive.parent.mkdir(parents=True, exist_ok=True)
        print("  downloading the ProteinGym substitution scans (43 MB)...")
        with requests.get(ARCHIVE, stream=True, timeout=120) as response:
            response.raise_for_status()
            with open(archive, "wb") as fh:
                for block in response.iter_content(1 << 20):
                    fh.write(block)
    with zipfile.ZipFile(archive) as bundle:
        for assay in wanted:
            member = next(n for n in bundle.namelist() if n.endswith(f"/{assay}.csv"))
            (out / f"{assay}.csv").write_bytes(bundle.read(member))


def load(exp: Experiment, assay: str | None = None) -> pd.DataFrame:
    """One row per single mutant: mutant (e.g. L67A), position, measured score, and the harmful label.

    Higher score = fitter. `harmful` uses ProteinGym's own cut-off for this scan.
    """
    assay = assay or exp.cfg["truth"]["assay"]
    table = pd.read_csv(truth_dir(exp) / f"{assay}.csv")
    table = table[table.mutant.str.match(MUTATION.pattern)].copy()  # single mutants only
    table["position"] = table.mutant.str[1:-1].astype(int)
    table["harmful"] = table.DMS_score_bin == 0
    table = table[table.position.between(exp.first_residue, exp.last_residue)]
    return table.rename(columns={"DMS_score": "measured"})[["mutant", "position", "measured", "harmful"]]


def check_sequence(exp: Experiment, construct: str) -> None:
    """The scan must have been done on the protein we are about to score."""
    table = load(exp)
    wrong = sorted({m for m in table.mutant if construct[int(m[1:-1]) - exp.first_residue] != m[0]})
    if wrong:
        raise SystemExit(f"The scan's wild type differs from the construct at: {wrong[:8]}")
    print(f"Measured scan {exp.cfg['truth']['assay']}: {len(table)} single mutants at "
          f"{table.position.nunique()} positions; {int(table.harmful.sum())} harmful, "
          f"{int((~table.harmful).sum())} tolerated")
    print("  its wild-type residues match the construct")


def draw_subset(exp: Experiment) -> None:
    """Choose, at random but reproducibly, the mutants that get a structure prediction.

    Predicting a structure for every mutant in the scan would take days, so a sample is drawn:
    equal numbers of harmful and tolerated mutants, from a fixed random seed, before anything is
    predicted. The draw is written down and never repeated.
    """
    record = exp.selection_dir / "variants.json"
    if record.exists():
        chosen = json.loads(record.read_text())
        print(f"Subset: already drawn ({len(chosen['variants'])} variants in selection/variants.json)")
        return
    spec = exp.cfg["subset"]
    first, last = spec["positions"]
    table = load(exp)
    pool = table[table.position.between(first, last)]
    rng = random.Random(spec["random_seed"])
    picked = {}
    for label, harmful in (("harmful", True), ("tolerated", False)):
        names = sorted(pool[pool.harmful == harmful].mutant)  # sorted, so the draw does not depend on file order
        picked[label] = sorted(rng.sample(names, spec[label]), key=lambda m: (int(m[1:-1]), m))
    named = [m for m in spec.get("always_include", []) if m not in picked["harmful"] + picked["tolerated"]]
    missing = [m for m in named if m not in set(table.mutant)]
    if missing:
        raise SystemExit(f"always_include names mutants the scan did not measure: {missing}")
    exp.selection_dir.mkdir(exist_ok=True)
    record.write_text(json.dumps({
        "how": f"random.Random({spec['random_seed']}).sample over the sorted mutant names of each class, "
               f"positions {first}-{last}; plus the named cases",
        "random_sample": picked,
        "named_cases": named,
        # The baseline, then every chosen mutant under its own name.
        "variants": {"wt": [], **{m: [m] for m in picked["harmful"] + picked["tolerated"] + named}},
    }, indent=2) + "\n")
    print(f"Subset drawn: {spec['harmful']} harmful + {spec['tolerated']} tolerated at random "
          f"(seed {spec['random_seed']}), plus {named or 'no'} named cases -> selection/variants.json")
