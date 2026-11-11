"""Step 4 - analyze: turn predictions into numbers, and measure the noise before reading the signal.

Three kinds of comparison, all measured the same way:

  mutation    baseline vs mutant, same network and same seed. The only difference is the mutation.
  seed_noise  the same sequence predicted twice with different random seeds. Biologically
              meaningless variation: the smallest difference worth paying attention to.
  model_noise the same sequence predicted by two of AlphaFold2's five networks.

A mutation "shows up" in a metric only if it moves that metric more than seed noise does.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from itertools import combinations

import numpy as np
import pandas as pd

from . import structure as st
from .config import Experiment
from .msa import read_a3m
from .prepare import reference_sequence

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# metric -> (column, compare by absolute value?)
METRICS = {
    "Cα RMSD, folded core (Å)": ("rmsd_core", False),
    "All-atom RMSD around the site (Å)": ("rmsd_local_all_atom", False),
    "|ΔpLDDT| at the site": ("delta_plddt_site", True),
    "|ΔpLDDT| whole protein": ("delta_plddt_mean", True),
}


def conservation(exp: Experiment) -> pd.DataFrame:
    """Per position: how often relatives keep the same amino acid, and how varied the column is."""
    _, records = read_a3m((exp.msa / f"{exp.baseline}.a3m").read_text())
    query = records[0][1]
    # In a3m, lowercase letters are insertions relative to the query; dropping them leaves
    # every row the same length as the query. Identical rows are counted once.
    rows = {"".join(c for c in seq if not c.islower()) for _, seq in records}
    rows = [row for row in rows if len(row) == len(query)]
    table = []
    for i, residue in enumerate(query):
        column = Counter(row[i] for row in rows if row[i] in AMINO_ACIDS)
        total = sum(column.values())
        entropy = -sum((n / total) * math.log2(n / total) for n in column.values())
        table.append(
            {
                "position": exp.first_residue + i,
                "residue": residue,
                "msa_depth": total,
                "identity": column[residue] / total,
                "conservation": 1 - entropy / math.log2(20),
                # ties broken alphabetically, so the table is the same on every run
                "most_common_substitutions": " ".join(
                    f"{aa}:{n / total:.0%}"
                    for aa, n in sorted(column.items(), key=lambda kv: (-kv[1], kv[0]))[:4] if aa != residue
                ),
                # Full profile: how often each amino acid occupies this position among relatives.
                **{f"freq_{aa}": column[aa] / total for aa in AMINO_ACIDS},
            }
        )
    return pd.DataFrame(table)


def _stats(values) -> dict:
    values = np.asarray(list(values), dtype=float)
    return {
        "median": round(float(np.median(values)), 3),
        "min": round(float(values.min()), 3),
        "max": round(float(values.max()), 3),
        "n": int(values.size),
    }


def _loop_verdict(to_baseline: float, to_this: float, rule: dict) -> str:
    """Which crystal a prediction's moved residues resemble, by the rule fixed in the experiment README."""
    if to_this < rule["near"] and to_baseline > rule["far"]:
        return "follows this crystal"
    if to_baseline < rule["near"] and to_this > rule["far"]:
        return "stays with the baseline crystal"
    return "in between"


def run(exp: Experiment) -> None:
    exp.results.mkdir(parents=True, exist_ok=True)
    if "truth" in exp.cfg:
        # A different kind of experiment: scores against a laboratory scan, not structures against structures.
        from . import agreement
        return agreement.run(exp)
    offset = exp.first_residue - 1
    first, last = exp.cfg["analysis"]["core_residues"]
    core = list(range(first, last + 1))
    cutoff = exp.cfg["analysis"]["neighbor_cutoff"]
    # Optional second condition for calling something an effect (see the experiment README).
    minimum_effect = exp.cfg["analysis"].get("minimum_effect")
    moved_threshold = exp.cfg["analysis"].get("moved_threshold", 2.0)
    baseline = exp.baseline

    reference = st.load_atoms(exp.reference_pdb, exp.cfg["reference"]["chain"], exp.reference_offset)
    reference_core = [r for r in core if r in reference]
    predictions = {
        name: st.load_predictions(exp.predictions / name, name, offset) for name in exp.variants
    }
    if not all(predictions.values()):
        raise SystemExit("No predictions found - run the predict step first")

    # --- 1. Each prediction on its own -------------------------------------------------------
    rows = []
    for name, group in predictions.items():
        for p in group:
            lowest = min(p.plddt, key=p.plddt.get)
            _, _, to_reference = st.superpose(st.ca(p.atoms, reference_core), st.ca(reference, reference_core))
            rows.append(
                {
                    "variant": name, "rank": p.rank, "model": p.model, "seed": p.seed,
                    "plddt_mean": np.mean(list(p.plddt.values())),
                    "plddt_min": p.plddt[lowest], "plddt_min_residue": lowest,
                    "ptm": p.ptm,
                    "rmsd_core_to_reference": to_reference,
                }
            )
    per_prediction = pd.DataFrame(rows).sort_values(["variant", "rank"])
    per_prediction.round(3).to_csv(exp.results / "predictions.csv", index=False)

    # --- 2. Pairwise comparisons: mutation vs noise -------------------------------------------
    by_key = {name: {(p.model, p.seed): p for p in group} for name, group in predictions.items()}
    comparisons, deviations = [], []
    summary_mutants = {}
    # A variant that changes no residue (same sequence, different template) is examined at the
    # positions the experiment is about: every site mutated in any variant.
    every_site = sorted({m.position for muts in exp.variants.values() for m in muts})
    for mutant in exp.mutants:
        sites = [m.position for m in exp.variants[mutant]] or every_site
        local = st.neighbors(reference, sites, cutoff)

        pairs = [("mutation", by_key[baseline][key], by_key[mutant][key]) for key in sorted(by_key[baseline])]
        for name in (baseline, mutant):
            pairs += [
                ("seed_noise", a, b)
                for a, b in combinations(predictions[name], 2)
                if a.model == b.model and a.seed != b.seed
            ]
        pairs += [
            ("model_noise", a, b)
            for a, b in combinations(predictions[baseline], 2)
            if a.model != b.model and a.seed == b.seed
        ]
        for kind, a, b in pairs:
            result = st.compare(a, b, core, sites, local)
            deviation = result.pop("deviation")
            comparisons.append({"mutant": mutant, "kind": kind, "a": a.label, "b": b.label, **result})
            deviations += [
                {"mutant": mutant, "kind": kind, "position": r, "ca_deviation": d}
                for r, d in deviation.items()
            ]

        table = pd.DataFrame([c for c in comparisons if c["mutant"] == mutant])
        metrics = {}
        for label, (column, absolute) in METRICS.items():
            values = table[column].abs() if absolute else table[column]
            effect = values[table.kind == "mutation"]
            noise = values[table.kind == "seed_noise"]
            metrics[label] = {
                "mutation": _stats(effect),
                "seed_noise": _stats(noise),
                "model_noise": _stats(values[table.kind == "model_noise"]),
                # The decision rule, fixed before looking at the data (see the experiment README).
                "exceeds_seed_noise": bool(effect.median() > noise.max()),
                "exceeds_model_noise": bool(effect.median() > values[table.kind == "model_noise"].median()),
            }
            m = metrics[label]
            if not m["exceeds_seed_noise"]:
                m["verdict"] = "within noise"
            elif minimum_effect != "model_noise":
                m["verdict"] = "above noise"
            elif m["exceeds_model_noise"]:
                m["verdict"] = "above noise and above the difference between networks"
            else:
                m["verdict"] = "detectable but smaller than the difference between networks"
        summary_mutants[mutant] = {
            "sites": sites,
            "neighbors_within_cutoff": local,
            "signed_delta_plddt_site_median": round(float(table[table.kind == "mutation"].delta_plddt_site.median()), 3),
            "signed_delta_plddt_mean_median": round(float(table[table.kind == "mutation"].delta_plddt_mean.median()), 3),
            "metrics": metrics,
            "distance_to_functional_sites_A": {
                label: round(st.side_chain_distance(reference, sites, [r for r in residues if r not in sites]), 1)
                for label, residues in exp.cfg["analysis"].get("functional_sites", {}).items()
            },
        }
    pd.DataFrame(comparisons).round(3).to_csv(exp.results / "comparisons.csv", index=False)

    # --- 3. Per-residue table: evidence for choosing a site, and where things moved ------------
    sequence = reference_sequence(exp)
    accessibility = st.relative_accessibility(
        exp.reference_pdb, exp.cfg["reference"]["chain"], sequence, exp.reference_offset
    )
    per_residue = conservation(exp)
    per_residue["relative_accessibility"] = per_residue.position.map(accessibility)
    for name, group in predictions.items():
        values = pd.DataFrame([p.plddt for p in group])
        per_residue[f"plddt_{name}_mean"] = per_residue.position.map(values.mean())
        per_residue[f"plddt_{name}_min"] = per_residue.position.map(values.min())
        per_residue[f"plddt_{name}_max"] = per_residue.position.map(values.max())
    dev = pd.DataFrame(deviations)
    for mutant in exp.mutants:
        for kind in ("mutation", "seed_noise"):
            grouped = dev[(dev.mutant == mutant) & (dev.kind == kind)].groupby("position").ca_deviation
            per_residue[f"ca_dev_{mutant}_{kind}_median"] = per_residue.position.map(grouped.median())
            per_residue[f"ca_dev_{mutant}_{kind}_max"] = per_residue.position.map(grouped.max())

    # --- 3b. Where a crystal structure of the mutant exists: did the prediction move the way the
    #         experiment did? ------------------------------------------------------------------
    for mutant, (path, chain) in exp.variant_references.items():
        crystal = st.load_atoms(path, chain, exp.reference_offset)
        shared = [r for r in reference_core if r in crystal and "CA" in crystal[r]]
        rot, tran, crystal_rmsd = st.superpose(st.ca(crystal, shared), st.ca(reference, shared))
        overlaid = st.move(crystal, rot, tran)
        shift = {r: float(np.linalg.norm(overlaid[r]["CA"] - reference[r]["CA"])) for r in shared}
        per_residue[f"experimental_ca_dev_{mutant}"] = per_residue.position.map(shift)
        moved = sorted(r for r, d in shift.items() if d > moved_threshold)
        still = [r for r in shared if r not in moved]
        rows = per_residue[per_residue.position.isin(moved)]
        truth = {
            "pdb": path.stem,
            "crystal_vs_crystal_rmsd": round(crystal_rmsd, 3),
            "largest_shift": round(max(shift.values()), 2),
            "largest_shift_residue": max(shift, key=shift.get),
            "moved_threshold": moved_threshold,
            "moved_residues": moved,
            "predicted_largest_shift": round(float(per_residue[f"ca_dev_{mutant}_mutation_median"].max()), 2),
            "predicted_largest_shift_residue": int(
                per_residue.position[per_residue[f"ca_dev_{mutant}_mutation_median"].idxmax()]
            ),
            # Every prediction against each crystal, over all residues both crystals have.
            "rmsd_to_baseline_crystal": {
                name: _stats(st.rmsd_after_fit(p.atoms, reference, shared, shared) for p in predictions[name])
                for name in (baseline, mutant)
            },
            "rmsd_to_mutant_crystal": {
                name: _stats(st.rmsd_after_fit(p.atoms, crystal, shared, shared) for p in predictions[name])
                for name in (baseline, mutant)
            },
        }
        if moved:
            truth.update({
                "experimental_shift_of_moved_residues": _stats(shift[r] for r in moved),
                "predicted_shift_of_moved_residues": _stats(rows[f"ca_dev_{mutant}_mutation_median"]),
                "seed_noise_at_moved_residues": _stats(rows[f"ca_dev_{mutant}_seed_noise_max"]),
                # Overlay on the residues that stay put, then look only at the ones that move:
                # is the mutant prediction's version of them nearer the mutant crystal or the baseline one?
                "moved_residues_rmsd_to_baseline_crystal": _stats(
                    st.rmsd_after_fit(p.atoms, reference, still, moved) for p in predictions[mutant]),
                "moved_residues_rmsd_to_mutant_crystal": _stats(
                    st.rmsd_after_fit(p.atoms, crystal, still, moved) for p in predictions[mutant]),
            })
            # One line per prediction: which crystal do its moved residues resemble?
            verdict_rule = exp.cfg["analysis"].get("loop_verdict")
            each = []
            for p in sorted(predictions[mutant], key=lambda p: (p.model, p.seed)):
                to_baseline = st.rmsd_after_fit(p.atoms, reference, still, moved)
                to_this = st.rmsd_after_fit(p.atoms, crystal, still, moved)
                each.append({"model": p.model, "seed": p.seed, "to_baseline_crystal": round(to_baseline, 2),
                             "to_this_crystal": round(to_this, 2),
                             "verdict": _loop_verdict(to_baseline, to_this, verdict_rule) if verdict_rule else None})
            truth["per_prediction"] = each
            if verdict_rule:
                verdicts = {e["verdict"] for e in each}
                truth["verdict"] = verdicts.pop() if len(verdicts) == 1 else "split: " + ", ".join(sorted(verdicts))
        summary_mutants[mutant]["experimental_truth"] = truth
    per_residue.round(3).to_csv(exp.results / "per_residue.csv", index=False)

    ranked = per_residue.conservation.rank(ascending=False, method="min").astype(int)
    for mutant, info in summary_mutants.items():
        new_residue = {m.position: m.new for m in exp.variants[mutant]}
        info["site_evidence"] = [
            {
                "position": int(row.position),
                "residue": row.residue,
                "new_residue": new_residue[row.position],
                "new_residue_frequency_in_msa": round(row[f"freq_{new_residue[row.position]}"], 4),
                "identity_in_msa": round(row.identity, 3),
                "conservation": round(row.conservation, 3),
                "conservation_rank": f"{ranked[i]} of {len(per_residue)}",
                "relative_accessibility": round(row.relative_accessibility, 3),
                "most_common_substitutions": row.most_common_substitutions,
            }
            for i, row in per_residue.iterrows()
            if row.position in new_residue
        ]

    # --- 3c. How much do two crystal structures of the *same* protein differ? The experimental
    #         counterpart of seed noise: a mutant only "moves" if it moves more than this. --------
    experimental_noise = None
    if exp.replicates:
        chain = exp.cfg["reference"]["chain"]
        rmsds, largest = [], []
        for path in exp.replicates:
            other = st.load_atoms(path, chain, exp.reference_offset)
            shared = [r for r in reference_core if r in other and "CA" in other[r]]
            rot, tran, rmsd = st.superpose(st.ca(other, shared), st.ca(reference, shared))
            overlaid = st.move(other, rot, tran)
            rmsds.append(rmsd)
            largest.append(max(float(np.linalg.norm(overlaid[r]["CA"] - reference[r]["CA"])) for r in shared))
        experimental_noise = {"structures": [p.stem for p in exp.replicates],
                              "ca_rmsd": _stats(rmsds), "largest_shift": _stats(largest)}

    # --- 3d. A panel chosen by rule: one line per mutant, and the count the experiment is after. --
    panel = None
    if "selection" in exp.cfg and (exp.selection_dir / "panel.csv").exists():
        roles = pd.read_csv(exp.selection_dir / "panel.csv").set_index("mutation").role.to_dict()
        floor = experimental_noise["largest_shift"]["max"] if experimental_noise else moved_threshold
        lines = []
        for mutant, info in summary_mutants.items():
            truth = info.get("experimental_truth")
            if not truth or mutant not in roles:
                continue
            fold = info["metrics"]["Cα RMSD, folded core (Å)"]
            line = {"mutant": mutant, "role": roles[mutant], "pdb": truth["pdb"],
                    "experimental_largest_shift": truth["largest_shift"],
                    "predicted_largest_shift": truth["predicted_largest_shift"],
                    "fold_verdict": fold["verdict"], "captured": None, "false_alarm": None}
            if roles[mutant] == "mover" and truth["moved_residues"]:
                half = truth["predicted_shift_of_moved_residues"]["median"] >= 0.5 * truth["experimental_shift_of_moved_residues"]["median"]
                nearer = (truth["moved_residues_rmsd_to_mutant_crystal"]["median"]
                          < truth["moved_residues_rmsd_to_baseline_crystal"]["median"])
                line["captured"] = bool(half and nearer)
            elif roles[mutant] == "non-mover":
                line["false_alarm"] = bool(fold["exceeds_seed_noise"] and fold["exceeds_model_noise"]
                                           and truth["predicted_largest_shift"] > floor)
            lines.append(line)
        table = pd.DataFrame(lines)
        table.to_csv(exp.results / "panel.csv", index=False)
        movers, others = table[table.role == "mover"], table[table.role == "non-mover"]
        panel = {
            "noise_floor_largest_shift": floor,
            "movers": len(movers), "captured": int(movers.captured.sum()),
            "non_movers": len(others), "false_alarms": int(others.false_alarm.sum()),
            "rank_correlation_predicted_vs_experimental": round(float(
                table.experimental_largest_shift.corr(table.predicted_largest_shift, method="spearman")), 2),
            "mutants": lines,
        }

    # --- 4. Summary -----------------------------------------------------------------------------
    summary = {
        "experiment": exp.cfg["experiment"]["id"],
        "baseline": baseline,
        "msa_unique_sequences": int(per_residue.msa_depth.max()),
        "variants": {
            name: {
                "predictions": len(group),
                "plddt_mean": _stats(group.plddt_mean),
                "plddt_min": _stats(group.plddt_min),
                "ptm": _stats(group.ptm),
                "rmsd_core_to_reference": _stats(group.rmsd_core_to_reference),
            }
            for name in [baseline, *exp.mutants]
            for group in [per_prediction[per_prediction.variant == name]]
        },
        "mutants": summary_mutants,
    }
    if experimental_noise:
        summary["experimental_noise"] = experimental_noise
    if panel:
        summary["panel"] = panel
    (exp.results / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")

    # --- What just happened ----------------------------------------------------------------------
    base = summary["variants"][baseline]
    print(f"Positive control - does the predictor reproduce the experiment for '{baseline}'?")
    print(f"  Cα RMSD to {exp.cfg['reference']['pdb']} over residues {first}-{last}: "
          f"median {base['rmsd_core_to_reference']['median']} Å (max {base['rmsd_core_to_reference']['max']} Å)")
    for name, info in summary["variants"].items():
        print(f"  {name}: mean pLDDT median {info['plddt_mean']['median']}, pTM median {info['ptm']['median']}")
    for mutant, info in summary_mutants.items():
        print(f"\n{mutant} vs {baseline} (median of matched pairs) against seed noise (max):")
        for label, m in info["metrics"].items():
            print(f"  {label:36s} {m['mutation']['median']:7.3f}  vs noise max {m['seed_noise']['max']:7.3f}  -> {m['verdict']}")
        truth = info.get("experimental_truth")
        if truth:
            print(f"  crystal {truth['pdb']} vs baseline crystal: largest shift {truth['largest_shift']} Å at residue "
                  f"{truth['largest_shift_residue']}; predicted largest shift {truth['predicted_largest_shift']} Å "
                  f"at residue {truth['predicted_largest_shift_residue']}")
        if truth and truth.get("verdict"):
            print(f"  over the residues that move in {truth['pdb']}: {truth['verdict']}")
    if panel:
        print(f"\nPanel: {panel['captured']} of {panel['movers']} movers captured; "
              f"{panel['false_alarms']} of {panel['non_movers']} non-movers raised a false alarm")
    print("\nwrote results/predictions.csv, comparisons.csv, per_residue.csv, summary.json")
