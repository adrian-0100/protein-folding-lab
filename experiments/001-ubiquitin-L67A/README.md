# Experiment 001: ubiquitin L67A

Written before the predictions were run. Results are in [results/REPORT.md](results/REPORT.md),
conclusions in [NOTEBOOK.md](NOTEBOOK.md).

## Question

Does replacing leucine 67 of human ubiquitin with alanine change the structure, or
the confidence, that AlphaFold2 predicts?

## Background

**Ubiquitin** is a 76-residue protein that cells attach to other proteins as a tag;
the best-known tag marks a protein for destruction by the proteasome. It is small,
folds on its own into one compact domain (a β-sheet wrapped around one α-helix),
has a high-quality experimental structure (PDB 1UBQ, 1.8 Å), and is among the most
conserved proteins known: yeast and human ubiquitin differ at 3 of 76 positions. That
makes it a good first protein: the predictor should get it right, and there is a
reference to check against.

**Leucine 67** sits on the last β-strand with its side chain pointing into the
hydrophobic core. **L67A** cuts that side chain down from four carbons to one. In
real proteins, this kind of large-to-small swap in the core leaves a cavity and
typically costs stability, while the folded shape barely moves (the classic
demonstration is the T4 lysozyme cavity series, Eriksson et al., *Science* 1992).

So the mutation was chosen to be one where biochemistry gives a clear expectation:
the protein should still fold into the same shape, but less stably. The analysis
checks the two premises (conserved, buried) from data instead of taking them on trust.

## Hypotheses

- **H1 (positive control).** AlphaFold2 reproduces the experimental structure of
  wild-type ubiquitin: Cα RMSD to 1UBQ below 1 Å over the folded core, mean pLDDT
  above 90.
- **H2 (structure).** L67A does not change the predicted fold: Cα RMSD between
  matched wild-type and mutant predictions stays within seed noise.
- **H3 (confidence).** If AlphaFold2 "notices" the damage to the core, pLDDT at and
  around residue 67 drops by more than seed noise.

H3 is the interesting one. Published tests on other proteins found that AlphaFold2's
output and pLDDT respond very little to single mutations, even destabilising ones
(Buel & Walters, *Nat. Struct. Mol. Biol.* 2022; Pak et al., *PLoS ONE* 2023). The
expectation is therefore that H3 fails. Either outcome is informative: a clear drop
would be a usable signal; no drop shows where this tool stops being useful.

## Method

All settings are in [experiment.toml](experiment.toml).

1. Sequence: residues 1-76 of UniProt P0CG48. Reference: PDB 1UBQ chain A.
2. One MSA search for the wild type. The mutant uses the same alignment with only
   the query row changed.
3. AlphaFold2-ptm through ColabFold, 5 models x 2 seeds = 10 predictions per
   variant, 3 recycles, no templates, no relaxation.
4. Comparisons, all measured identically:
   - **mutation**: wild type vs L67A, same model and same seed (10 matched pairs)
   - **seed noise**: same sequence, same model, different seed (10 pairs)
   - **model noise**: wild type, same seed, different model (20 pairs)

## Metrics

| Metric | What it captures |
|---|---|
| Cα RMSD to 1UBQ, residues 1-72 | Is the prediction right at all? |
| Cα RMSD between two predictions, residues 1-72 | Did the fold move? |
| All-atom RMSD of residues within 6 Å of residue 67 | Did the packing around the site rearrange? |
| ΔpLDDT at residue 67, and over the whole protein | Did the model's confidence change? |
| Conservation and solvent exposure of residue 67 | Was the site a sensible choice? |

Residues 73-76 are left out of RMSD: they form a flexible tail in solution.

## Decision rule

A metric shows an effect of the mutation only if its **median over the matched
pairs is larger than the largest value seen in seed noise**. Anything else is
reported as "within noise".

## What this experiment cannot tell us

- Whether L67A ubiquitin is less stable, misfolds, or works in a cell. That needs a
  stability predictor (a different kind of tool) and then a measurement.
- Anything general about AlphaFold2 and mutations. It is one mutation in one protein.
