"""Step (variant-effect experiments) - score: how surprising is each mutation, by three yardsticks.

None of these predicts a structure. Each gives one number per mutation, higher = more acceptable:

  blosum62            chemistry only: how readily this kind of swap is tolerated in proteins in general
  alignment_log_odds  evolution of this protein: how often the new amino acid appears at this position
                      among its relatives, compared with the original
  language_model      a protein language model (ESM-2): hide the position, ask the model how likely
                      each amino acid is there given the rest of the sequence, compare new with original

The scores are computed for every mutant in the laboratory scan, without looking at the measurements.
"""

from __future__ import annotations

import math
from collections import Counter

import pandas as pd
from Bio.Align import substitution_matrices

from . import truth
from .config import Experiment
from .msa import read_a3m

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def blosum62(mutants: list[str]) -> list[float]:
    matrix = substitution_matrices.load("BLOSUM62")
    return [float(matrix[m[0]][m[-1]]) for m in mutants]


def alignment_log_odds(exp: Experiment, mutants: list[str]) -> list[float]:
    """log( (count of new + 1) / (count of original + 1) ) in that column of the alignment."""
    _, records = read_a3m((exp.msa / f"{exp.baseline}.a3m").read_text())
    query = records[0][1]
    rows = {"".join(c for c in seq if not c.islower()) for _, seq in records}  # drop insertions, deduplicate
    rows = [row for row in rows if len(row) == len(query)]
    columns = [Counter(row[i] for row in rows) for i in range(len(query))]
    scores = []
    for m in mutants:
        column = columns[int(m[1:-1]) - exp.first_residue]
        scores.append(math.log((column[m[-1]] + 1) / (column[m[0]] + 1)))
    return scores


def language_model(exp: Experiment, mutants: list[str], model_name: str) -> list[float]:
    """Masked-marginal score: log p(new) - log p(original) with that position hidden from the model."""
    import torch  # heavy and optional: only needed for this score
    from transformers import AutoTokenizer, EsmForMaskedLM

    sequence = exp.sequence(exp.baseline)
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmForMaskedLM.from_pretrained(model_name).eval()
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    model.to(device)
    positions = sorted({int(m[1:-1]) for m in mutants})
    encoded = tokenizer([sequence] * len(positions), return_tensors="pt")
    for row, position in enumerate(positions):
        # +1: the tokenizer puts a start token in front of the sequence
        encoded["input_ids"][row, position - exp.first_residue + 1] = tokenizer.mask_token_id
    log_p = {}
    with torch.no_grad():
        for start in range(0, len(positions), 16):
            batch = {k: v[start:start + 16].to(device) for k, v in encoded.items()}
            logits = model(**batch).logits.float().cpu()
            for offset, position in enumerate(positions[start:start + 16]):
                log_p[position] = torch.log_softmax(logits[offset, position - exp.first_residue + 1], dim=-1)
    index = {aa: tokenizer.convert_tokens_to_ids(aa) for aa in AMINO_ACIDS}
    print(f"  language model {model_name} on {device}: {len(positions)} positions scored")
    return [float(log_p[int(m[1:-1])][index[m[-1]]] - log_p[int(m[1:-1])][index[m[0]]]) for m in mutants]


def run(exp: Experiment) -> None:
    wanted = exp.cfg["scores"]
    mutants = sorted(truth.load(exp).mutant, key=lambda m: (int(m[1:-1]), m))
    table = pd.DataFrame({"mutant": mutants, "position": [int(m[1:-1]) for m in mutants]})
    if wanted.get("blosum62"):
        table["blosum62"] = blosum62(mutants)
    if wanted.get("alignment_log_odds"):
        table["alignment_log_odds"] = alignment_log_odds(exp, mutants)
    if wanted.get("language_model"):
        table["language_model"] = language_model(exp, mutants, wanted["language_model"])
    out = exp.dir / "scores"
    out.mkdir(exist_ok=True)
    table.round(4).to_csv(out / "scores.csv", index=False)
    print(f"Scored {len(table)} mutants with: {', '.join(c for c in table.columns if c not in ('mutant', 'position'))}")
    print("wrote scores/scores.csv (no measurements in this file: the comparison happens in analyze)")
