"""Step 0 (panel experiments) - select: choose mutants by a written rule, from experimental structures only.

Hand-picking which mutants to test is how a result gets biased. Here the PDB is surveyed for
every crystal structure of a single mutant of the protein, each is measured against the
reference crystal structure, and a rule fixed in experiment.toml decides which ones enter the
experiment. No prediction is involved at any point.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import requests

from . import structure as st
from .config import Experiment
from .prepare import _download, fetch_structure, UNIPROT

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
DETAILS = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id
 symmetry{space_group_name_H_M}
 rcsb_entry_info{resolution_combined}
 nonpolymer_entities{nonpolymer_comp{chem_comp{id}}}
 polymer_entities{entity_poly{pdbx_seq_one_letter_code pdbx_seq_one_letter_code_can}
  rcsb_polymer_entity_container_identifiers{auth_asym_ids}}}}"""


def survey_pdb(accession: str) -> list[dict]:
    """Every X-ray entry containing this UniProt protein, with what the rule needs to know about it."""
    query = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match", "value": accession}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}}]},
        "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 5000}},
    }
    found = requests.get(SEARCH, params={"json": json.dumps(query)}, timeout=60)
    found.raise_for_status()
    ids = [hit["identifier"] for hit in found.json()["result_set"]]
    entries = []
    for i in range(0, len(ids), 100):
        reply = requests.post(GRAPHQL, json={"query": DETAILS, "variables": {"ids": ids[i:i + 100]}}, timeout=120)
        reply.raise_for_status()
        entries += reply.json()["data"]["entries"]
    return entries


def run(exp: Experiment) -> None:
    rule = exp.cfg["selection"]
    out = exp.selection_dir
    record = out / "variants.json"
    if record.exists():
        print(f"Panel already selected ({len(json.loads(record.read_text())['variants']) - 1} mutants); "
              "delete selection/variants.json to select again")
        return
    out.mkdir(exist_ok=True)
    exp.inputs.mkdir(exist_ok=True)

    # The baseline sequence, and the reference crystal structure.
    accession = exp.cfg["protein"]["uniprot"]
    entry = exp.inputs / f"uniprot_{accession}.fasta"
    if not entry.exists():
        entry.write_text(_download(UNIPROT.format(accession=accession)))
    construct = exp.construct()
    reference_id, chain = exp.cfg["reference"]["pdb"], exp.cfg["reference"]["chain"]
    fetch_structure(exp.reference_pdb)
    reference = st.load_atoms(exp.reference_pdb, chain, exp.reference_offset)
    first, last = exp.cfg["analysis"]["core_residues"]
    core = [r for r in range(first, last + 1) if r in reference]

    # 1. Survey.
    entries = survey_pdb(accession)
    print(f"PDB: {len(entries)} X-ray entries contain {accession}")
    allowed = set(rule["allowed_bound_molecules"])
    rows, reasons = [], {}
    for e in entries:
        pdb = e["rcsb_id"]
        entity = e["polymer_entities"][0]
        sequence = entity["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "")
        n = min(len(sequence), len(construct))
        changed = [f"{a}{i + exp.first_residue}{b}" for i, (a, b) in enumerate(zip(construct[:n], sequence[:n])) if a != b]
        bound = sorted({x["nonpolymer_comp"]["chem_comp"]["id"] for x in (e["nonpolymer_entities"] or [])})
        resolution = (e["rcsb_entry_info"]["resolution_combined"] or [99.0])[0]
        if len(e["polymer_entities"]) != 1 or len(entity["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]) != 1:
            why = "more than one chain or protein"
        elif not len(construct) - 4 <= len(sequence) <= len(construct):
            why = "different length (a fusion or a fragment)"
        elif not e["symmetry"] or e["symmetry"]["space_group_name_H_M"] != rule["space_group"]:
            why = "different crystal form"
        elif len(changed) != 1:
            why = "not a single mutant" if changed else "unmutated"
        elif resolution > rule["max_resolution"]:
            why = f"resolution {resolution} Å"
        elif "(" in entity["entity_poly"]["pdbx_seq_one_letter_code"]:
            why = "chemically modified residue"
        elif set(bound) - allowed:
            why = "bound molecule: " + ", ".join(sorted(set(bound) - allowed))
        else:
            why = ""
        reasons[why or "qualifies"] = reasons.get(why or "qualifies", 0) + 1
        if len(changed) == 1:
            rows.append({"mutation": changed[0], "pdb": pdb, "resolution": resolution,
                         "bound_molecules": " ".join(bound), "excluded_because": why})
    for why, count in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"  {count:4d}  {why}")

    # 2. Measure every single mutant against the reference crystal (also the excluded ones: the
    #    table then shows what each exclusion cost).
    table = pd.DataFrame(rows)
    measured = []
    for row in table.itertuples():
        path = out / "pdb" / f"{row.pdb}.pdb"
        try:
            fetch_structure(path)
            atoms = st.load_atoms(path, chain, exp.reference_offset)
            shared = [r for r in core if r in atoms and "CA" in atoms[r]]
            rot, tran, rmsd = st.superpose(st.ca(atoms, shared), st.ca(reference, shared))
            moved = st.move(atoms, rot, tran)
            shift = {r: float(np.linalg.norm(moved[r]["CA"] - reference[r]["CA"])) for r in shared}
            worst = max(shift, key=shift.get)
            measured.append((round(rmsd, 3), round(shift[worst], 2), worst, sum(d > 1.0 for d in shift.values())))
        except Exception as problem:  # an entry that cannot be read is recorded, not fatal
            measured.append((None, None, None, None))
            table.loc[row.Index, "excluded_because"] = f"unreadable ({type(problem).__name__})"
    table[["ca_rmsd", "largest_shift", "at_residue", "residues_over_1A"]] = measured
    table["position"] = table.mutation.str[1:-1].astype(int)
    table = table.sort_values(["largest_shift", "resolution"], ascending=[False, True])
    table.to_csv(out / "survey.csv", index=False)

    # 3. Apply the rule. One structure per mutation: the best resolved of those that qualify.
    ok = table[table.excluded_because == ""].sort_values("resolution").drop_duplicates("mutation")
    movers = ok[ok.largest_shift > rule["mover_threshold"]].sort_values("largest_shift", ascending=False)
    movers = movers.head(rule["max_movers"])
    still = ok[ok.largest_shift < rule["non_mover_threshold"]]
    non_movers = []
    for mover in movers.itertuples():
        pool = still[~still.mutation.isin(non_movers)]
        # nearest in the sequence; ties go to the better resolved structure
        pick = pool.assign(distance=(pool.position - mover.position).abs()).sort_values(["distance", "resolution"]).iloc[0]
        non_movers.append(pick.mutation)
    chosen = pd.concat([movers.assign(role="mover"), ok[ok.mutation.isin(non_movers)].assign(role="non-mover")])
    chosen.to_csv(out / "panel.csv", index=False)

    record.write_text(json.dumps({
        "how": f"foldlab select: single mutants of the baseline in {rule['space_group']} from the PDB, "
               f"rule in experiment.toml [selection]; see survey.csv and panel.csv",
        "variants": {"wt": [], **{row.mutation: [row.mutation] for row in chosen.itertuples()}},
        "references": {row.mutation: {"pdb": row.pdb, "chain": chain} for row in chosen.itertuples()},
    }, indent=2) + "\n")

    print(f"\n{table.mutation.nunique()} different single mutations have a structure; "
          f"{len(ok)} qualify under the rule.")
    print(f"Movers (largest Cα shift > {rule['mover_threshold']} Å):")
    for row in movers.itertuples():
        print(f"  {row.mutation:7s} {row.pdb}  {row.largest_shift:5.2f} Å at residue {int(row.at_residue)}  ({row.resolution} Å)")
    print(f"Non-movers (< {rule['non_mover_threshold']} Å), the nearest one in the sequence to each mover:")
    for row in chosen[chosen.role == "non-mover"].itertuples():
        print(f"  {row.mutation:7s} {row.pdb}  {row.largest_shift:5.2f} Å  ({row.resolution} Å)")
    left_out = table[(table.excluded_because != "") & (table.largest_shift > rule["mover_threshold"])
                     & ~table.mutation.isin(ok.mutation)]
    left_out = left_out.sort_values("resolution").drop_duplicates("mutation")
    if len(left_out):
        print("Large movers the rule left out:")
        for row in left_out.itertuples():
            print(f"  {row.mutation:7s} {row.pdb}  {row.largest_shift:5.2f} Å  - {row.excluded_because}")
    print("wrote selection/survey.csv, panel.csv, variants.json")
