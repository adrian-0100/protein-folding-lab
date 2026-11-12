"""Step 5 - report: figures, an interactive 3D viewer, and REPORT.md with the numbers."""

from __future__ import annotations

import itertools
import json
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, PDBParser

from . import structure as st
from .analyze import METRICS
from .config import Experiment

# Colour follows the entity: the baseline is always blue, a mutant always orange,
# and anything that is context (noise, the experiment, other residues) is grey.
BASELINE, MUTANT, CONTEXT = "#2a78d6", "#eb6834", "#8a8983"
MUTANT_COLORS = [MUTANT, "#1baf7a", "#eda100"]   # 3D viewer only, where several mutants share one scene
SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "text.color": INK, "axes.labelcolor": INK_2, "xtick.color": INK_2, "ytick.color": INK_2,
        "axes.edgecolor": GRID, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
        "legend.frameon": False, "figure.dpi": 200,
    }
)


def _mark_sites(ax, sites: list[int], label: str) -> None:
    for site in sites:
        ax.axvline(site, color=INK_2, linewidth=0.8, linestyle=(0, (2, 2)), zorder=1)
    ax.annotate(label, (sites[0], 1.0), xycoords=("data", "axes fraction"), xytext=(0, 3),
                textcoords="offset points", ha="center", va="bottom", fontsize=9, color=INK)


def fig_confidence(exp, table, mutant, sites, out):
    fig, ax = plt.subplots(figsize=(8, 3.6))
    for name, color in ((exp.baseline, BASELINE), (mutant, MUTANT)):
        ax.fill_between(table.position, table[f"plddt_{name}_min"], table[f"plddt_{name}_max"],
                        color=color, alpha=0.18, linewidth=0)
        ax.plot(table.position, table[f"plddt_{name}_mean"], color=color, linewidth=2, label=name)
    for level, text in ((90, "very high"), (70, "confident")):
        ax.axhline(level, color=INK_2, linewidth=0.6)
        ax.text(table.position.max() + 0.6, level, text, va="center", fontsize=8, color=INK_2)
    _mark_sites(ax, sites, mutant)
    ax.set(xlabel="Residue position", ylabel="pLDDT (0-100)", ylim=(40, 100),
           xlim=(table.position.min(), table.position.max()))
    ax.set_title("Predicted confidence per residue", pad=18)
    ax.legend(loc="lower left", ncols=2, title="line = mean, band = range over all predictions",
              title_fontsize=8, alignment="left")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_deviation(table, mutant, sites, out):
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.fill_between(table.position, 0, table[f"ca_dev_{mutant}_seed_noise_max"], color=CONTEXT,
                    alpha=0.3, linewidth=0, label="seed noise (max)")
    ax.plot(table.position, table[f"ca_dev_{mutant}_mutation_median"], color=MUTANT, linewidth=2,
            label=f"{mutant} vs baseline (median of matched pairs)")
    _mark_sites(ax, sites, mutant)
    ax.set(xlabel="Residue position", ylabel="Cα displacement (Å)", ylim=(0, None),
           xlim=(table.position.min(), table.position.max()))
    ax.set_title("How far each residue moves after overlaying the folded core", pad=18)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_truth(table, mutant, sites, pdb, out):
    """The same quantity twice: how far each residue moves in the crystals, and in the predictions."""
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.fill_between(table.position, 0, table[f"ca_dev_{mutant}_seed_noise_max"], color=CONTEXT,
                    alpha=0.3, linewidth=0, label="prediction: seed noise (max)")
    ax.plot(table.position, table[f"experimental_ca_dev_{mutant}"], color=INK, linewidth=2,
            label=f"experiment: {pdb} vs wild-type crystal")
    ax.plot(table.position, table[f"ca_dev_{mutant}_mutation_median"], color=MUTANT, linewidth=2,
            label=f"prediction: {mutant} vs baseline (median of matched pairs)")
    _mark_sites(ax, sites, mutant)
    ax.set(xlabel="Residue position (UniProt numbering)", ylabel="Cα displacement (Å)", ylim=(0, None),
           xlim=(table.position.min(), table.position.max()))
    ax.set_title("Does the prediction move where the experiment moves?", pad=18)
    # Below the plot: the peaks can be anywhere, so no corner inside it is safe.
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.2), ncols=1)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_site_evidence(table, mutant, sites, out):
    fig, axes = plt.subplots(2, 1, figsize=(8, 4.6), sharex=True)
    colors = [MUTANT if p in sites else CONTEXT for p in table.position]
    panels = (
        ("conservation", "Conservation across relatives (1 = never changes)"),
        ("relative_accessibility", "Exposure to water in the crystal structure (0 = buried)"),
    )
    for ax, (column, title) in zip(axes, panels):
        # Chain ends can exceed the tabulated maximum exposure; cap the bar at 1.
        ax.bar(table.position, table[column].clip(upper=1), width=0.7, color=colors)
        ax.set_title(title)
        ax.set_ylim(0, 1.12)
        ax.set_yticks([0, 0.5, 1])
        ax.grid(axis="x", visible=False)
        for site in sites:
            value = table.loc[table.position == site, column].iloc[0]
            ax.annotate(f"{table.loc[table.position == site, 'residue'].iloc[0]}{site}: {value:.2f}",
                        (site, min(value, 1)), xytext=(0, 14), textcoords="offset points",
                        ha="center", va="bottom", fontsize=9, color=INK,
                        bbox={"facecolor": SURFACE, "edgecolor": "none", "pad": 1.5},
                        arrowprops={"arrowstyle": "-", "color": INK_2, "linewidth": 0.8})
    axes[1].set_xlabel("Residue position")
    axes[1].set_xlim(table.position.min() - 1, table.position.max() + 1)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_effect_vs_noise(comparisons, mutant, out):
    kinds = (("mutation", f"{mutant} vs baseline", MUTANT), ("seed_noise", "seed noise", CONTEXT),
             ("model_noise", "model noise", CONTEXT))
    fig, axes = plt.subplots(1, len(METRICS), figsize=(11, 3.0), sharey=True)
    rng = np.random.default_rng(0)
    for ax, (label, (column, absolute)) in zip(axes, METRICS.items()):
        for y, (kind, _, color) in enumerate(kinds):
            values = comparisons.loc[comparisons.kind == kind, column]
            values = values.abs() if absolute else values
            ax.scatter(values, y + rng.uniform(-0.16, 0.16, len(values)), s=26, color=color,
                       edgecolor=SURFACE, linewidth=0.8, zorder=3)
            ax.plot([values.median()] * 2, [y - 0.3, y + 0.3], color=INK, linewidth=1.5, zorder=4)
        ax.set_title(textwrap.fill(label, 22), fontsize=9.5)
        ax.set_xlim(0, None)
        ax.set_yticks(range(len(kinds)), [k[1] for k in kinds])
        ax.set_ylim(len(kinds) - 0.5, -0.5)
        ax.grid(axis="y", visible=False)
    fig.suptitle("Is the mutation's effect bigger than noise?  (dot = one pair of predictions, bar = median)",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def _aligned_pdb(source: Path, chain_id, offset, reference: st.Atoms, core: list[int], out: Path) -> str:
    """Save a copy of `source` overlaid on the reference, so all structures share one frame."""
    structure = PDBParser(QUIET=True).get_structure(source.stem, source)
    model = structure[0]
    chain = model[chain_id] if chain_id else next(iter(model))
    for other in [c for c in model if c.id != chain.id]:
        model.detach_child(other.id)
    for residue in [r for r in chain if r.id[0] != " "]:
        chain.detach_child(residue.id)
    atoms = st.load_atoms(source, chain_id, offset)
    shared = [r for r in core if r in reference and r in atoms]
    rot, tran, _ = st.superpose(st.ca(atoms, shared), st.ca(reference, shared))
    for atom in structure.get_atoms():
        atom.transform(rot, tran)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out))
    return out.read_text()


VIEWER = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"></script>
<style>
  :root { --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --line:#e6e5e0; }
  @media (prefers-color-scheme: dark) { :root { --surface:#1a1a19; --ink:#f0efec; --ink2:#c3c2b7; --line:#3a3a37; } }
  body { margin:0; background:var(--surface); color:var(--ink); font:15px/1.45 system-ui, sans-serif; }
  main { max-width:960px; margin:0 auto; padding:16px; }
  h1 { font-size:18px; margin:0 0 4px; } p { margin:0 0 12px; color:var(--ink2); }
  #viewer { position:relative; width:100%; height:min(68vh, 620px); border:1px solid var(--line); border-radius:8px; }
  .row { display:flex; flex-wrap:wrap; gap:8px 20px; margin:12px 0; align-items:center; }
  label { display:inline-flex; gap:6px; align-items:center; cursor:pointer; }
  .swatch { width:12px; height:12px; border-radius:3px; display:inline-block; }
  .legend span { margin-right:14px; white-space:nowrap; }
</style>
</head>
<body>
<main>
  <h1>__TITLE__</h1>
  <p>Drag to rotate, scroll to zoom. All structures are overlaid on the folded core of the experimental structure. The mutated residue is drawn as sticks.</p>
  <div class="row" id="toggles"></div>
  <div class="row">
    <label><input type="radio" name="mode" value="structure" checked> Colour by structure</label>
    <label><input type="radio" name="mode" value="plddt"> Colour predictions by confidence (pLDDT)</label>
  </div>
  <div id="viewer"></div>
  <div class="row legend" id="legend"></div>
</main>
<script>
const structures = __STRUCTURES__;
const bands = [[90, "#0053D6", "very high (>90)"], [70, "#65CBF3", "confident (70-90)"], [50, "#FFDB13", "low (50-70)"], [0, "#FF7D45", "very low (<50)"]];
const dark = matchMedia("(prefers-color-scheme: dark)").matches;
const viewer = $3Dmol.createViewer("viewer", { backgroundColor: dark ? "#1a1a19" : "#fcfcfb" });
structures.forEach(s => viewer.addModel(s.pdb, "pdb"));
const toggles = document.getElementById("toggles");
structures.forEach((s, i) => {
  toggles.insertAdjacentHTML("beforeend",
    `<label><input type="checkbox" data-i="${i}" checked><span class="swatch" style="background:${s.color}"></span>${s.label}</label>`);
});
const plddtColor = b => bands.find(([min]) => b >= min)[1];
function draw() {
  const mode = document.querySelector("input[name=mode]:checked").value;
  structures.forEach((s, i) => {
    const on = toggles.querySelector(`[data-i="${i}"]`).checked;
    const byConfidence = mode === "plddt" && s.prediction;
    const paint = byConfidence ? { colorfunc: a => plddtColor(a.b) } : { color: s.color };
    viewer.setStyle({ model: i }, on ? { cartoon: { ...paint, opacity: s.prediction ? 1 : 0.75 } } : {});
    if (on) viewer.addStyle({ model: i, resi: s.sites }, { stick: { color: s.color, radius: 0.22 } });
  });
  document.getElementById("legend").innerHTML = mode === "plddt"
    ? bands.map(([, c, t]) => `<span><span class="swatch" style="background:${c}"></span> ${t}</span>`).join("")
      + "<span>(the experimental structure keeps its grey)</span>"
    : "";
  viewer.render();
}
document.querySelectorAll("input").forEach(el => el.addEventListener("change", draw));
viewer.zoomTo();
draw();
</script>
</body>
</html>
"""


def _range(stats: dict, digits: int = 2) -> str:
    return f"{stats['median']:.{digits}f} ({stats['min']:.{digits}f}-{stats['max']:.{digits}f})"


def write_report(exp: Experiment, summary: dict, provenance: dict) -> None:
    cfg = exp.cfg
    pdb = cfg["reference"]["pdb"]
    first, last = cfg["analysis"]["core_residues"]
    lines = [
        f"# Results: {cfg['experiment']['id']}",
        "",
        "Generated by `foldlab report` from the files in `results/`. Do not edit by hand: "
        "interpretation belongs in `NOTEBOOK.md`.",
        "",
        f"**Question.** {cfg['experiment']['question']}",
        "",
        "## Does the predictor work on this protein? (positive control)",
        "",
        f"Each prediction is overlaid on the experimental structure {pdb} using the Cα atoms of residues "
        f"{first}-{last}. Values are median (min-max) over all predictions of that variant.",
        "",
        f"| Variant | Predictions | Mean pLDDT | Lowest pLDDT | pTM | Cα RMSD to {pdb} (Å) |",
        "|---|---|---|---|---|---|",
    ]
    for name, v in summary["variants"].items():
        lines.append(
            f"| {name} | {v['predictions']} | {_range(v['plddt_mean'], 1)} | {_range(v['plddt_min'], 1)} "
            f"| {_range(v['ptm'])} | {_range(v['rmsd_core_to_reference'])} |"
        )
    for mutant, info in summary["mutants"].items():
        lines += ["", f"## {mutant}: why this site?", "",
                  "| Residue | Same amino acid in relatives | Most common substitutions | Relatives that carry the new amino acid | Conservation (0-1) | Conservation rank | Exposure to water (0-1) |",
                  "|---|---|---|---|---|---|---|"]
        for site in info["site_evidence"]:
            lines.append(
                f"| {site['residue']}{site['position']} | {site['identity_in_msa']:.0%} | {site['most_common_substitutions']} "
                f"| {site['new_residue']}: {site['new_residue_frequency_in_msa']:.1%} | {site['conservation']:.2f} "
                f"| {site['conservation_rank']} | {site['relative_accessibility']:.2f} |"
            )
        lines += ["", f"Alignment depth: {summary['msa_unique_sequences']} unique sequences. "
                  f"Residues within {cfg['analysis']['neighbor_cutoff']} Å of the site in {pdb}: "
                  + ", ".join(str(r) for r in info["neighbors_within_cutoff"]) + ".", "",
                  "Closest side-chain approach to known functional residues:", ""]
        lines += [f"- {label}: {distance} Å" for label, distance in info["distance_to_functional_sites_A"].items()]
        lines += ["", "![site evidence](site_evidence_" + mutant + ".png)", "",
                  f"## {mutant} vs {summary['baseline']}: effect against noise", "",
                  "Decision rule (fixed in advance): the mutation shows up in a metric only if the median over "
                  "matched pairs is larger than the largest value seen between two seeds of the same sequence."
                  + (" To count as more than detectable it must also exceed the median difference between two "
                     "of AlphaFold2's networks." if cfg["analysis"].get("minimum_effect") == "model_noise" else ""), "",
                  "| Metric | Mutation, median (min-max) | Seed noise, median (min-max) | Model noise, median (min-max) | Verdict |",
                  "|---|---|---|---|---|"]
        for label, m in info["metrics"].items():
            verdict = m["verdict"] if m["verdict"] == "within noise" else f"**{m['verdict']}**"
            lines.append(f"| {label.replace('|', chr(92) + '|')} | {_range(m['mutation'])} | {_range(m['seed_noise'])} "
                         f"| {_range(m['model_noise'])} | {verdict} |")
        lines += ["", f"Signed change in pLDDT, {mutant} minus {summary['baseline']} (median of matched pairs): "
                  f"{info['signed_delta_plddt_site_median']:+.2f} at the site, "
                  f"{info['signed_delta_plddt_mean_median']:+.2f} over the whole protein.", "",
                  f"![effect vs noise](effect_vs_noise_{mutant}.png)", "",
                  f"![confidence](confidence_{mutant}.png)", "",
                  f"![displacement](displacement_{mutant}.png)"]
        truth = info.get("experimental_truth")
        if truth:
            base = summary["baseline"]
            lines += ["", f"## {mutant}: prediction against the crystal structure {truth['pdb']}", "",
                      f"What the experiment shows: overlaid on {pdb}, the crystal structure {truth['pdb']} differs by "
                      f"{truth['crystal_vs_crystal_rmsd']:.2f} Å (Cα RMSD), with the largest shift, "
                      f"{truth['largest_shift']:.1f} Å, at residue {truth['largest_shift_residue']}. "
                      + (f"Residues that move more than {truth['moved_threshold']:g} Å: "
                         + ", ".join(str(r) for r in truth["moved_residues"]) + "."
                         if truth["moved_residues"] else f"No residue moves more than {truth['moved_threshold']:g} Å."),
                      "",
                      f"What the prediction shows: the largest shift between matched {base} and {mutant} predictions "
                      f"is {truth['predicted_largest_shift']:.2f} Å, at residue {truth['predicted_largest_shift_residue']}.", ""]
            if truth["moved_residues"]:
                lines += ["| Over the residues that move in the experiment | Median (min-max), Å |", "|---|---|",
                          f"| Shift between the two crystal structures | {_range(truth['experimental_shift_of_moved_residues'])} |",
                          f"| Shift between {base} and {mutant} predictions | {_range(truth['predicted_shift_of_moved_residues'])} |",
                          f"| Largest seed noise there | {_range(truth['seed_noise_at_moved_residues'])} |",
                          f"| {mutant} predictions, distance from {pdb} ({base} crystal) | {_range(truth['moved_residues_rmsd_to_baseline_crystal'])} |",
                          f"| {mutant} predictions, distance from {truth['pdb']} | {_range(truth['moved_residues_rmsd_to_mutant_crystal'])} |",
                          ""]
            lines += [f"| Cα RMSD over the whole protein, median (min-max), Å | vs {pdb} ({base} crystal) | vs {truth['pdb']} |",
                      "|---|---|---|"]
            lines += [f"| {name} predictions | {_range(truth['rmsd_to_baseline_crystal'][name])} | {_range(truth['rmsd_to_mutant_crystal'][name])} |"
                      for name in (base, mutant)]
            if truth.get("per_prediction") and truth.get("verdict"):
                lines += ["", f"Each prediction, over the residues that move in {truth['pdb']} "
                              f"(distance in Å after overlaying the residues that do not move). "
                              f"Verdict for this variant: **{truth['verdict']}**.", "",
                          f"| Network | Seed | From the {base} crystal ({pdb}) | From {truth['pdb']} | Reading |",
                          "|---|---|---|---|---|"]
                lines += [f"| {e['model']} | {e['seed']} | {e['to_baseline_crystal']:.2f} | {e['to_this_crystal']:.2f} | {e['verdict']} |"
                          for e in truth["per_prediction"]]
            lines += ["", f"![experiment vs prediction](experiment_vs_prediction_{mutant}.png)"]

    # Experiments about templates: one table of sequence against what the network was handed.
    handed = {name: exp.templates(name) for name in summary["variants"] if exp.templates(name)}
    if handed:
        lines += ["", "## Sequence against template", "",
                  "Which crystal the moved residues of each variant's predictions resemble.", "",
                  "| Variant | Mutations | Template handed over | Scored against | Verdict |", "|---|---|---|---|---|"]
        for name in summary["variants"]:
            info = summary["mutants"].get(name, {})
            truth = info.get("experimental_truth") or {}
            mutations = ", ".join(str(m) for m in exp.variants[name]) or "none"
            hint = handed.get(name)
            hint = "none" if not hint else "database search" if hint == "search" else ", ".join(hint)
            lines.append(f"| {name} | {mutations} | {hint} | {truth.get('pdb', '-')} | {truth.get('verdict', '-')} |")

    # Experiments on a panel chosen by rule: the count, and the table behind it.
    noise = summary.get("experimental_noise")
    if noise:
        lines += ["", "## How much do crystals of the same protein differ?", "",
                  f"{len(noise['structures'])} other crystal structures of the unmutated protein "
                  f"({', '.join(noise['structures'])}) against {pdb}: Cα RMSD {_range(noise['ca_rmsd'])} Å, "
                  f"largest single Cα shift {_range(noise['largest_shift'])} Å. A mutant has moved only if it "
                  f"moves more than this."]
    panel = summary.get("panel")
    if panel:
        lines += ["", "## The panel", "",
                  f"**{panel['captured']} of {panel['movers']} movers captured. "
                  f"{panel['false_alarms']} of {panel['non_movers']} non-movers raised a false alarm.** "
                  f"Rank correlation between experimental and predicted largest shift, over the whole panel: "
                  f"{panel['rank_correlation_predicted_vs_experimental']:.2f}.", "",
                  "| Mutant | Role | Crystal | Largest shift in the crystal (Å) | Largest shift in the prediction (Å) | Fold-level verdict | Captured | False alarm |",
                  "|---|---|---|---|---|---|---|---|"]
        yes_no = lambda v: "-" if v is None else "yes" if v else "no"
        lines += [f"| {m['mutant']} | {m['role']} | {m['pdb']} | {m['experimental_largest_shift']:.2f} | "
                  f"{m['predicted_largest_shift']:.2f} | {m['fold_verdict']} | {yes_no(m['captured'])} | {yes_no(m['false_alarm'])} |"
                  for m in panel["mutants"]]
        lines += ["", "![panel](panel.png)"]
        fig_panel(panel, exp.results / "panel.png")
    software = provenance.get("software", {})
    lines += ["", "## Provenance", "",
              f"- ColabFold {software.get('colabfold')}, AlphaFold2 weights `{provenance.get('weights')}`, "
              f"JAX {software.get('jax')}, Python {software.get('python')}, {software.get('platform')}",
              "- Settings: " + ", ".join(f"{k}={v}" for k, v in provenance.get("settings", {}).items()),
              "- Every variant used the same alignment; only the query row differs (see `msa/`).",
              "- 3D view: open `results/viewer.html` in a browser.", ""]
    (exp.results / "REPORT.md").write_text("\n".join(lines))


def fig_panel(panel: dict, out) -> None:
    """One dot per mutant: how far the crystal moved against how far the prediction moved."""
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    table = pd.DataFrame(panel["mutants"])
    for role, color in (("mover", MUTANT), ("non-mover", BASELINE)):
        part = table[table.role == role]
        ax.scatter(part.experimental_largest_shift, part.predicted_largest_shift, s=46, color=color,
                   edgecolor=SURFACE, linewidth=1, label=role, zorder=3)
        for row in part.itertuples():
            ax.annotate(row.mutant, (row.experimental_largest_shift, row.predicted_largest_shift),
                        xytext=(5, 4), textcoords="offset points", fontsize=8, color=INK_2)
    top = max(table.experimental_largest_shift.max(), table.predicted_largest_shift.max()) * 1.1
    ax.plot([0, top], [0, top], color=INK_2, linewidth=0.8, linestyle=(0, (2, 2)), zorder=1)
    ax.text(top * 0.98, top * 0.93, "prediction = experiment", ha="right", fontsize=8, color=INK_2)
    ax.axvspan(0, panel["noise_floor_largest_shift"], color=CONTEXT, alpha=0.15, linewidth=0)
    ax.text(panel["noise_floor_largest_shift"] / 2, top * 0.5, "experimental\nnoise", ha="center", fontsize=8, color=INK_2)
    ax.set(xlim=(0, top), ylim=(0, top), xlabel="Largest Cα shift in the crystal (Å)",
           ylabel="Largest Cα shift in the prediction (Å)")
    ax.set_title("Does the prediction move when the crystal does?")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def run(exp: Experiment) -> None:
    if "truth" in exp.cfg:
        from . import agreement
        return agreement.report(exp)
    results = exp.results
    summary = json.loads((results / "summary.json").read_text())
    provenance_file = results / "provenance.json"
    provenance = json.loads(provenance_file.read_text()) if provenance_file.exists() else {}
    per_residue = pd.read_csv(results / "per_residue.csv")
    comparisons = pd.read_csv(results / "comparisons.csv")

    for mutant, info in summary["mutants"].items():
        sites = info["sites"]
        fig_site_evidence(per_residue, mutant, sites, results / f"site_evidence_{mutant}.png")
        fig_confidence(exp, per_residue, mutant, sites, results / f"confidence_{mutant}.png")
        fig_deviation(per_residue, mutant, sites, results / f"displacement_{mutant}.png")
        fig_effect_vs_noise(comparisons[comparisons.mutant == mutant], mutant,
                            results / f"effect_vs_noise_{mutant}.png")
        truth = info.get("experimental_truth")
        if truth:
            fig_truth(per_residue, mutant, sites, truth["pdb"], results / f"experiment_vs_prediction_{mutant}.png")
        print(f"wrote {5 if truth else 4} figures for {mutant}")

    # 3D viewer: the experiment, the top-ranked baseline prediction, and each mutant's matched prediction.
    aligned = results / "aligned"
    aligned.mkdir(exist_ok=True)
    offset = exp.first_residue - 1
    first, last = exp.cfg["analysis"]["core_residues"]
    core = list(range(first, last + 1))
    chain = exp.cfg["reference"]["chain"]
    reference = st.load_atoms(exp.reference_pdb, chain, exp.reference_offset)
    pdb = exp.cfg["reference"]["pdb"]
    top = min(st.load_predictions(exp.predictions / exp.baseline, exp.baseline, offset), key=lambda p: p.rank)
    structures = [
        {"label": f"{pdb} (experiment)", "color": CONTEXT, "prediction": False,
         "pdb": _aligned_pdb(exp.reference_pdb, chain, exp.reference_offset, reference, core, aligned / f"{pdb}.pdb")},
        {"label": f"{exp.baseline} (predicted, model {top.model} seed {top.seed})", "color": BASELINE, "prediction": True,
         "pdb": _aligned_pdb(top.path, None, offset, reference, core, aligned / f"{exp.baseline}.pdb")},
    ]
    sites = sorted({s for info in summary["mutants"].values() for s in info["sites"]})
    for mutant, color in zip(summary["mutants"], itertools.cycle(MUTANT_COLORS)):
        match = next(p for p in st.load_predictions(exp.predictions / mutant, mutant, offset)
                     if (p.model, p.seed) == (top.model, top.seed))
        structures.append(
            {"label": f"{mutant} (predicted, same model and seed)", "color": color, "prediction": True,
             "pdb": _aligned_pdb(match.path, None, offset, reference, core, aligned / f"{mutant}.pdb")})
    for mutant, (path, mutant_chain) in exp.variant_references.items():
        structures.append(
            {"label": f"{path.stem} ({mutant}, experiment)", "color": INK_2, "prediction": False,
             "pdb": _aligned_pdb(path, mutant_chain, exp.reference_offset, reference, core, aligned / f"{path.stem}.pdb")})
    # Each file keeps its own residue numbering: predictions start at 1, the experiment as deposited.
    for s in structures:
        shift = offset if s["prediction"] else exp.reference_offset
        s["sites"] = [site - shift for site in sites]
    html = (VIEWER.replace("__TITLE__", f"{exp.cfg['protein']['name']}: {exp.cfg['experiment']['id']}")
            .replace("__STRUCTURES__", json.dumps(structures)))
    (results / "viewer.html").write_text(html)
    print("wrote results/viewer.html")

    write_report(exp, summary, provenance)
    print("wrote results/REPORT.md")
