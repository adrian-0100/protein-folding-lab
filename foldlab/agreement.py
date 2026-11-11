"""Analysis and report for variant-effect experiments: which score agrees with the laboratory?

Every score is turned into the same two numbers against the measured scan:

  rank correlation (Spearman)  rank the mutants by the score and by the measurement; how alike are the two
                               rankings? 1 = identical, 0 = unrelated.
  AUC                          pick one harmful and one tolerated mutant at random; how often does the score
                               put them in the right order? 0.5 = coin flip, 1 = always.

All scores are oriented so that higher = better tolerated.
"""

from __future__ import annotations

import json
import textwrap

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from . import structure as st
from . import truth
from .config import Experiment

SCORES = {
    "blosum62": "S0 chemistry only (BLOSUM62)",
    "alignment_log_odds": "S1 alignment log-odds",
    "language_model": "S2 language model (ESM-2)",
    "af2_plddt_site": "S3 AlphaFold2: change in pLDDT at the site",
}
AF2_EXTRA = {
    "af2_plddt_mean": "AlphaFold2: change in mean pLDDT",
    "af2_ca_shift": "AlphaFold2: minus the Cα shift of the fold",
}


def auc(score: pd.Series, harmful: pd.Series) -> float:
    """Chance that a tolerated mutant scores higher than a harmful one (ties count half)."""
    ranks = rankdata(score)
    tolerated = ~harmful.to_numpy()
    n_tol, n_harm = tolerated.sum(), (~tolerated).sum()
    return float((ranks[tolerated].sum() - n_tol * (n_tol + 1) / 2) / (n_tol * n_harm))


def structure_scores(exp: Experiment) -> pd.DataFrame:
    """S3 for every mutant that has predictions: mutant minus wild type, same network and seed."""
    offset = exp.first_residue - 1
    first, last = exp.cfg["analysis"]["core_residues"]
    core = list(range(first, last + 1))
    base = {(p.model, p.seed): p for p in st.load_predictions(exp.predictions / exp.baseline, exp.baseline, offset)}
    rows = []
    for name in exp.mutants:
        site = exp.variants[name][0].position
        pairs = [(base[(p.model, p.seed)], p) for p in st.load_predictions(exp.predictions / name, name, offset)
                 if (p.model, p.seed) in base]
        if not pairs:
            continue
        at_site = [m.plddt[site] - w.plddt[site] for w, m in pairs]
        overall = [np.mean(list(m.plddt.values())) - np.mean(list(w.plddt.values())) for w, m in pairs]
        shift = [st.superpose(st.ca(m.atoms, core), st.ca(w.atoms, core))[2] for w, m in pairs]
        rows.append({"mutant": name, "predictions": len(pairs),
                     "af2_plddt_site": np.mean(at_site), "af2_plddt_site_spread": np.std(at_site),
                     "af2_plddt_mean": np.mean(overall), "af2_ca_shift": -np.mean(shift)})
    return pd.DataFrame(rows)


def _agreement(table: pd.DataFrame, column: str) -> dict:
    rows = table.dropna(subset=[column])
    return {"n": len(rows), "spearman": round(float(spearmanr(rows[column], rows.measured).statistic), 3),
            "auc": round(auc(rows[column], rows.harmful), 3)}


def run(exp: Experiment) -> None:
    measured = truth.load(exp)
    scores = pd.read_csv(exp.dir / "scores" / "scores.csv")
    table = measured.merge(scores.drop(columns="position"), on="mutant", how="left")
    record = json.loads((exp.selection_dir / "variants.json").read_text())
    sample = set(record["random_sample"]["harmful"]) | set(record["random_sample"]["tolerated"])
    table["in_sample"] = table.mutant.isin(sample)
    af2 = structure_scores(exp)
    if len(af2):
        table = table.merge(af2, on="mutant", how="left")
    table.round(4).to_csv(exp.results / "scores_with_measurements.csv", index=False)

    available = [c for c in [*SCORES, *AF2_EXTRA] if c in table.columns]
    sequence_scores = [c for c in available if not c.startswith("af2_")]
    results = {
        "all_mutants": {c: _agreement(table, c) for c in sequence_scores},
        # The like-for-like comparison: every score on the mutants that also have a structure prediction.
        "random_sample": {c: _agreement(table[table.in_sample], c) for c in available},
    }

    # The ceiling: how well does a second, independent laboratory scan agree with the first?
    ceiling = None
    second = exp.cfg["truth"].get("ceiling_assay")
    if second:
        other = truth.load(exp, second).rename(columns={"measured": "second_scan"})[["mutant", "second_scan"]]
        both = table.merge(other, on="mutant")
        ceiling = {"assay": second, **_agreement(both, "second_scan")}

    named = table[table.mutant.isin(record.get("named_cases", []))]
    cases = []
    for row in named.itertuples():
        case = {"mutant": row.mutant, "measured": round(row.measured, 3), "harmful": bool(row.harmful),
                "measured_percentile": round(float((table.measured < row.measured).mean() * 100))}
        for c in available:
            value = getattr(row, c)
            pool = table[c].dropna()
            case[c] = None if pd.isna(value) else round(float(value), 3)
            # where this mutant falls among all mutants that have this score (0 = worst)
            case[f"{c}_percentile"] = None if pd.isna(value) or len(pool) < 20 else round(float((pool < value).mean() * 100))
        cases.append(case)

    # Decision rules, as written in the experiment README before the run.
    rules = {}
    everyone, drawn = results["all_mutants"], results["random_sample"]
    if "alignment_log_odds" in everyone:
        rules["H1 evolution carries the signal (S1 AUC >= 0.70 on all mutants)"] = everyone["alignment_log_odds"]["auc"] >= 0.70
    if {"alignment_log_odds", "language_model"} <= everyone.keys():
        rules["H2 the language model does at least as well (S2 rho >= S1 rho - 0.05)"] = (
            everyone["language_model"]["spearman"] >= everyone["alignment_log_odds"]["spearman"] - 0.05)
    if {"af2_plddt_site", "blosum62"} <= drawn.keys():
        s3, s0 = drawn["af2_plddt_site"]["auc"], drawn["blosum62"]["auc"]
        rules["H3 AlphaFold2's confidence does not separate them (S3 AUC < 0.65, or not above S0)"] = (
            True if s3 < 0.65 or s3 <= s0 else False if s3 >= 0.75 else "inconclusive at this sample size")

    summary = {"experiment": exp.cfg["experiment"]["id"], "assay": exp.cfg["truth"]["assay"],
               "mutants_measured": len(table), "harmful": int(table.harmful.sum()),
               "sample_size": int(table.in_sample.sum()), "agreement": results, "ceiling": ceiling,
               "named_cases": cases, "decision_rules": rules}
    (exp.results / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")

    print(f"{len(table)} measured mutants ({summary['harmful']} harmful); {summary['sample_size']} in the random sample")
    print(f"{'score':46s} {'all: rho':>9s} {'AUC':>6s}   {'sample: rho':>11s} {'AUC':>6s}")
    for c in available:
        a, b = everyone.get(c), drawn[c]
        left = f"{a['spearman']:9.2f} {a['auc']:6.2f}" if a else f"{'-':>9s} {'-':>6s}"
        print(f"{({**SCORES, **AF2_EXTRA})[c]:46s} {left}   {b['spearman']:11.2f} {b['auc']:6.2f}")
    if ceiling:
        print(f"{'ceiling: a second laboratory scan':46s} {ceiling['spearman']:9.2f} {ceiling['auc']:6.2f}")
    for rule, verdict in rules.items():
        print(f"  {rule}: {'holds' if verdict is True else 'fails' if verdict is False else verdict}")
    print("wrote results/scores_with_measurements.csv, summary.json")


def report(exp: Experiment) -> None:
    import matplotlib.pyplot as plt

    from .report import BASELINE, CONTEXT, INK, INK_2, MUTANT  # also applies the shared figure style

    results = exp.results
    summary = json.loads((results / "summary.json").read_text())
    table = pd.read_csv(results / "scores_with_measurements.csv")
    names = {**SCORES, **AF2_EXTRA}
    shown = [c for c in SCORES if c in table.columns]

    # Figure 1: each score against the measurement.
    fig, axes = plt.subplots(1, len(shown), figsize=(3.3 * len(shown), 3.4), sharey=True, squeeze=False)
    for ax, c in zip(axes[0], shown):
        rows = table.dropna(subset=[c])
        for harmful, color, label in ((False, BASELINE, "tolerated"), (True, MUTANT, "harmful")):
            part = rows[rows.harmful == harmful]
            ax.scatter(part[c], part.measured, s=9, color=color, alpha=0.55, linewidth=0, label=label)
        where = "all_mutants" if c in summary["agreement"]["all_mutants"] else "random_sample"
        stats = summary["agreement"][where][c]
        ax.set_title(f"{textwrap.fill(names[c], 30)}\nρ = {stats['spearman']:.2f}, AUC = {stats['auc']:.2f}, n = {stats['n']}",
                     fontsize=9)
        ax.set_xlabel("score (higher = better tolerated)")
    axes[0][0].set_ylabel("measured growth (higher = fitter)")
    axes[0][0].legend(loc="lower right", markerscale=2)
    fig.tight_layout()
    fig.savefig(results / "score_vs_measurement.png")
    plt.close(fig)

    # Figure 2: AUC of every score on the same sample, with the laboratory's own ceiling.
    drawn = summary["agreement"]["random_sample"]
    order = [c for c in SCORES if c in drawn]
    fig, ax = plt.subplots(figsize=(7, 0.55 * len(order) + 1.6))
    values = [drawn[c]["auc"] for c in order]
    ax.barh(range(len(order)), values, height=0.6, color=[MUTANT if c.startswith("af2_") else CONTEXT for c in order])
    for y, v in enumerate(values):
        ax.text(v + 0.01, y, f"{v:.2f}", va="center", fontsize=9, color=INK)
    ax.axvline(0.5, color=INK_2, linewidth=0.8, linestyle=(0, (2, 2)))
    ax.text(0.5, -0.75, "coin flip", ha="center", fontsize=8, color=INK_2)
    if summary["ceiling"]:
        ax.axvline(summary["ceiling"]["auc"], color=INK, linewidth=1)
        ax.text(summary["ceiling"]["auc"], -0.75, "a second lab scan", ha="center", fontsize=8, color=INK)
    ax.set_yticks(range(len(order)), [names[c] for c in order])
    ax.set_ylim(len(order) - 0.5, -1.1)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel(f"AUC on the {summary['sample_size']} randomly drawn mutants")
    ax.set_title("How often does each score order a harmful and a tolerated mutant correctly?")
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(results / "auc_by_score.png")
    plt.close(fig)

    lines = [f"# Results: {summary['experiment']}", "",
             "Generated by `foldlab report`. Interpretation belongs in `NOTEBOOK.md`.", "",
             f"**Question.** {exp.cfg['experiment']['question']}", "",
             f"Measured scan: `{summary['assay']}`, {summary['mutants_measured']} single mutants, "
             f"{summary['harmful']} of them harmful. Structure predictions were made for a random sample of "
             f"{summary['sample_size']}.", "",
             "## Agreement with the laboratory", "",
             "ρ: rank correlation with measured growth. AUC: how often the score orders a harmful and a "
             "tolerated mutant correctly (0.5 = coin flip).", "",
             "| Score | All mutants: ρ | AUC | n | Random sample: ρ | AUC | n |", "|---|---|---|---|---|---|---|"]
    everyone = summary["agreement"]["all_mutants"]
    for c, stats in drawn.items():
        a = everyone.get(c)
        left = f"{a['spearman']:.2f} | {a['auc']:.2f} | {a['n']}" if a else "- | - | -"
        lines.append(f"| {names[c]} | {left} | {stats['spearman']:.2f} | {stats['auc']:.2f} | {stats['n']} |")
    if summary["ceiling"]:
        c = summary["ceiling"]
        lines.append(f"| *Ceiling: a second laboratory scan (`{c['assay']}`)* | {c['spearman']:.2f} | {c['auc']:.2f} | {c['n']} | | | |")
    lines += ["", "![scores against the measurement](score_vs_measurement.png)", "",
              "![AUC by score](auc_by_score.png)", "", "## Decision rules", ""]
    lines += [f"- {rule}: **{'holds' if v is True else 'fails' if v is False else v}**"
              for rule, v in summary["decision_rules"].items()]
    if summary["named_cases"]:
        lines += ["", "## The named cases", "",
                  "Percentile: where the mutant falls among all mutants on that scale (0 = the worst).", "",
                  "| Mutant | Measured growth | Label | " + " | ".join(names[c] for c in drawn) + " |",
                  "|---|---|---|" + "---|" * len(drawn)]
        for case in summary["named_cases"]:
            cells = [("-" if case.get(c) is None else f"{case[c]:.2f}")
                     + ("" if case.get(f"{c}_percentile") is None else f" (pct {case[f'{c}_percentile']})") for c in drawn]
            lines.append(f"| {case['mutant']} | {case['measured']:.2f} (pct {case['measured_percentile']}) | "
                         f"{'harmful' if case['harmful'] else 'tolerated'} | " + " | ".join(cells) + " |")
    (results / "REPORT.md").write_text("\n".join(lines) + "\n")
    print("wrote results/score_vs_measurement.png, auc_by_score.png, REPORT.md")
