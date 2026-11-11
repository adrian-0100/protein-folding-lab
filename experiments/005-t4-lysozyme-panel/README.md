# Experiment 005: from one case to a rate

**Status: prepared, not run (2026-11-10).** The PDB has been surveyed, the rule has chosen the panel, and structures and alignment are in place. Nothing has been predicted. This is the protocol. The selection rule, hypotheses
and decision rules are fixed here before any prediction is made.

## Question

Experiment 002 found one mutation whose crystal structure moves and whose
prediction does not. Is that typical? Across a panel of mutations with known
crystal structures, how often does AlphaFold2's prediction move when the
experiment moves, and stay put when it does not?

## Why this, why now

Two mutations are an anecdote. A dozen, chosen by a rule and not by hand, give a
rate, and a rate is what a report can stand on. This is also the experiment best
suited to a student: it has a selection step, a table, a plot, and a result that
is a count.

## Background

**The best-studied protein for this question.** T4 lysozyme is a 164-residue
enzyme from a virus that infects bacteria. From the 1980s onward, Brian Matthews's
laboratory made hundreds of single mutants of it, measured their stability, and
solved their crystal structures, to learn how proteins tolerate change. Despite the
name it is a different protein from human lysozyme (experiment 002), not a variant
of it.

**Surveyed for this proposal**, using only experimental structures, no predictions:

| | |
|---|---|
| Crystal structures containing T4 lysozyme in the PDB | 764 |
| ... that are the enzyme alone, one chain, in the standard crystal form (P 32 2 1) | 483 |
| ... of which single mutants of the standard laboratory background | 232 structures, 156 different mutations at 64 positions |
| Resolution of those | 1.0 to 2.5 Å, median 1.8 |

These figures are from the pilot survey made for the proposal. The survey the
experiment itself ran, with every structure and the reason each was kept or left
out, is [selection/survey.csv](selection/survey.csv).

**Same crystal form matters.** How molecules pack in a crystal can itself bend a
protein. Because all 232 structures share one crystal form, a difference between
a mutant and the reference is far more likely to come from the mutation than from
the packing. Experiment 002 could not say that as firmly.

**An experimental noise floor.** The reference protein itself has been
crystallised many times. Eight other structures of the *identical* protein differ
from the reference structure by a largest Cα shift of 0.30 to 1.11 Å (median 0.35).
So a mutant only "moves" if it moves clearly more than 1.1 Å. This is the
experimental counterpart of the seed noise used in 001 and 002, and it comes from
the laboratory.

**What the survey shows already.** Most mutations move nothing: the median largest
shift is about 0.5 Å, inside the noise floor. Only 8 of 156 exceed 1.5 Å, and 89 stay
under 0.5 Å. Proteins absorb most single changes without rearranging. The D67H
movement of experiment 002 is the exception, not the rule.

## Three traps, found while preparing this

1. **The laboratory protein is not the UniProt sequence.** Nearly every T4
   lysozyme structure differs from UniProt P00720 at positions 12 (Gly, not Arg)
   and 137 (Arg, not Ile), and most carry two more engineered changes, C54T and
   C97A, that remove the cysteines (the background called WT\*). The baseline here
   must be WT\*, taken from the reference structure 1L63, not from UniProt.
2. **The PDB's own mutation count is wrong for these entries.** It reports zero
   mutations for structures that have four. Mutations were found by comparing
   sequences directly.
3. **Some "mutants" carry more than a mutation.** Two of the eight largest movers
   have a chemical spin label attached to the mutated residue, and several
   structures of the cavity mutant L99A have a small molecule bound in the cavity,
   which moves the protein by itself (L99A ranges from 0.45 to 3.1 Å across its
   structures). These must be excluded by rule.

## Selection rule

Applied to the survey table, in this order, before any prediction:

1. Single mutation relative to WT\*; one chain; crystal form P 32 2 1.
2. Resolution 2.1 Å or better.
3. No modified residues, and no bound molecules other than ions and the reducing
   agents used in crystallisation (Cl, K, Na, phosphate, sulphate, BME, HED).
4. For a mutation with several qualifying structures, the best-resolved one.
5. **Movers:** every remaining mutation whose largest Cα shift exceeds 1.5 Å,
   up to six, largest first.
6. **Non-movers:** for each mover, the remaining mutation with largest shift under
   0.5 Å at the nearest residue position, without repeats. This pairs every mover
   with a mutation in the same part of the protein that did not move anything.

### What the rule chose

Run on 2026-11-10, before any prediction. Of 169 single mutations with a structure,
135 qualify. The panel is in [selection/panel.csv](selection/panel.csv):

| Mover | Crystal | Largest Cα shift | Paired non-mover | Crystal | Largest Cα shift |
|---|---|---|---|---|---|
| V111M | 1CV1 | 3.00 Å | V111A | 235L | 0.38 Å |
| V111I | 2L78 | 2.77 Å | T109D | 1L62 | 0.24 Å |
| L84A | 247L | 2.29 Å | A82S | 123L | 0.24 Å |
| V103M | 1CUQ | 1.57 Å | M102K | 1L54 | 0.35 Å |

Four movers, not six: that is all the rule found. Larger movers exist and were left
out for stated reasons: L118C and I9C carry a spin label, L99G has hexanediol bound,
A129W is at 2.5 Å resolution, and several others are in a different crystal form or
have more than one molecule in the crystal. Every one is listed with its reason in
the survey table.

One pairing is exactly the kind of contrast this experiment is for: at position 111,
methionine and isoleucine move the backbone by about 3 Å and alanine moves it by
0.4 Å.

## Hypotheses

- **H1 (positive control).** AlphaFold2 reproduces WT\*: Cα RMSD to 1L63 below 1 Å,
  mean pLDDT above 90.
- **H2.** For non-movers, the prediction does not move either.
- **H3.** For movers, the prediction does not follow the experiment: it stays at
  the WT\* conformation, as it did for D67H.

Expectation: all three hold.

## Decision rules

For each mutant, over the residues that shift more than 1 Å in its crystal,
measured as in experiment 002:

1. A mover is **captured** if the predicted shift of those residues is at least
   half the experimental shift (medians), and the mutant's predictions are closer
   to the mutant's own crystal than to the WT\* crystal there.
2. A non-mover is a **false alarm** if any fold-level metric passes the
   effect-against-noise rule of 002 (above seed noise and above the difference
   between networks) with a largest predicted shift over 1.1 Å.
3. H3 holds if at most one mover is captured. It fails if at least two thirds of
   the movers are captured (with four movers: three). Anything between is reported
   as mixed.
4. H2 holds if at most one non-mover is a false alarm.

Also reported, without a verdict attached: a plot of predicted against
experimental largest shift for all mutants in the panel, and the rank correlation
between them.

## Method

1. `select` (new step): rebuild the survey table from the PDB and apply the rule.
   The table and the chosen panel are saved in `selection/`.
2. `prepare`: WT\* sequence from 1L63; one FASTA per mutant; each mutant's crystal
   structure as its own reference; the eight other WT\* structures for the
   experimental noise floor.
3. `msa`: one search for WT\*.
4. `predict`: WT\* and the 8 mutants, 5 networks x 2 seeds, no templates.
   90 predictions of 164 residues.
5. `analyze`, `report`: as in 002, for every mutant, plus the panel summary.

## What is in place

Settings: [experiment.toml](experiment.toml).

| Piece | Where |
|---|---|
| `background` mutations applied before the variants | `foldlab/config.py`, `prepare.py` |
| The `select` step: survey the PDB, apply the rule, write the panel | `foldlab/select.py` |
| Experimental noise floor from replicate crystal structures | `foldlab/analyze.py` |
| Panel summary, table and figure | `foldlab/analyze.py`, `report.py` |

No new software. Tested so far: the panel summary on experiment 002 treated as a
panel of two. Not yet run on this panel.

To run it:

```bash
infra/remote-gpu-pc/remote run experiments/005-t4-lysozyme-panel
```

## Cost

90 predictions of 164 residues (nine sequences, five networks, two seeds).
Extrapolating from experiment 002 (130 residues, a little under a minute each at
ten at a time), roughly two hours on the remote machine. Not measured.

## Known before running

- **The movers are not independent.** Nearly all of them shift the same stretch,
  around residues 108-113, next to the cavity that the well-known L99A mutation
  opens. The panel will show how AlphaFold2 handles that one flexible region more
  than how it handles movement in general.
- **Few movers exist.** The rule found four. The count of captured movers will be
  out of four and should be read that way.
- **All of these structures were in AlphaFold2's training data**, as in 002.
- **Hexanediol and similar additives** sit on the border of rule 3. The rule
  excluded them, which cost the largest clean mover, L99G (4.7 Å).
- **Movement is not harm.** Stability was measured for most of these mutants and
  is a separate question, closer to experiment 003.

## What we will have at the end

A sentence of the form "AlphaFold2 followed the experiment in k of n mutations
that move the backbone, and raised a false alarm in m of n that do not", with the
table and figure behind it.
