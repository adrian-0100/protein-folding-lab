# Lab notebook: 004-lysozyme-templates

Numbers come from `results/`; what they mean is written here by a person.
Residue numbers are given as in the literature (mature enzyme), with the UniProt
number used by the pipeline in brackets where it helps: add 18. The tables in
`results/` use UniProt numbers.

---

## 2026-11-11 · First run

**Goal.** Find out why experiment 002 failed: hand AlphaFold2 the crystal structure
of D67H as a starting hint and see whether the two loops move, with the wild-type
sequence as control.

**What was done.**

1. Question, hypotheses and decision rules written in [README.md](README.md) on
   2026-11-10, before any prediction with a template existed.
2. `prepare`, `msa` (2026-11-10): same protein and references as 002; its alignment
   copied over (5,264 sequences); the two template structures, 1LZ1 and 1LYY,
   placed in `templates/`.
3. `predict` on the remote GPU PC, networks 1 and 2, two seeds. The 12 predictions
   without a template were copied from 002. The 28 new ones ran as 14 processes,
   ten at a time, 36 minutes in total. No process failed or stalled.
4. `analyze`, `report`.

**Were the templates really used?** Yes. This was the one piece of the pipeline
never run before, so it was checked before reading any result. Every process
records what it was handed: `found templates: ['1lyy_A']` for the two 1LYY
conditions, `['1lz1_A']` for the two 1LZ1 conditions, and twenty hits for each
search. And the predictions change with the template (below), which they could not
if it had been dropped.

**What the search found.** Twenty structures for each of the three sequences, the
same set in slightly different order (I56T's list has 1OUH where the others have
1YAN). By their titles all are engineered mutants of human lysozyme, most of them from
studies of stability and catalysis. None is 1LYY, 1LOZ (the I56T crystal) or 1LZ1. AlphaFold2
uses at most four templates; the output lists the twenty hits and does not say
which four it used.

**Deviation from plan.** None in the method. Two things after the run:

- The report called the crystal each condition is measured against "the wt_search
  crystal", "the wt_1LYY crystal" and so on. It is 1LYY every time. The wording was
  fixed in `foldlab/report.py` and `report` rerun on the laptop; `summary.json` is
  identical and only those label lines of `REPORT.md` changed.
- The report says "stays with the baseline crystal" where the protocol says
  "wild-type loops". Same rule, same thing.

**A detail worth knowing.** Predictions with a template are slower. Ten processes
side by side took 24 to 25 minutes each; without templates, in 002, the same jobs
took 13 to 16.

### Results

**Positive control: holds.** The wild-type predictions without a template are
0.41 Å (0.40-0.42) from 1LZ1, mean pLDDT 98.0, pTM 0.91.

**The table.** Thirteen residues sit more than 2 Å apart in the two crystals:
45-50, 68-73 and 103 (UniProt 63-68, 86-91, 121). Over those, after overlaying on
everything else, the two crystals are 6.15 Å apart. Each cell is the median of
four predictions (range in brackets), in Å.

| Template handed over | Sequence | From the wild-type crystal 1LZ1 | From the D67H crystal 1LYY | By the rule |
|---|---|---|---|---|
| None (from 002) | wild type | 0.95 (0.93-0.99) | 5.62 (5.61-5.66) | not scored in the report |
| None (from 002) | D67H | 0.82 (0.76-0.87) | 5.74 (5.72-5.77) | wild-type loops |
| Found by search | wild type | 0.60 (0.59-0.61) | 6.10 (6.09-6.10) | wild-type loops |
| Found by search | D67H | 0.52 (0.50-0.52) | 6.13 (6.13-6.14) | wild-type loops |
| 1LZ1 | wild type | 0.42 (0.38-0.46) | 5.97 (5.93-6.01) | wild-type loops |
| 1LZ1 | D67H | 0.39 (0.35-0.40) | 6.01 (5.97-6.05) | wild-type loops |
| **1LYY** | **wild type** | **1.40 (1.34-1.44)** | **5.26 (5.24-5.29)** | **wild-type loops** |
| **1LYY** | **D67H** | **1.28 (1.24-1.33)** | **5.32 (5.31-5.33)** | **wild-type loops** |

The first row was computed for this entry with the pipeline's own function, since
the report does not score the wild type without a template. In every cell the four
predictions agree, so every cell gets a verdict.

**H1, handed its own crystal the D67H sequence takes up the D67H loops: fails.**
The four predictions are 1.24 to 1.33 Å from the wild-type crystal and 5.31 to
5.33 Å from the crystal they were handed. The protocol expected H1 to hold.

**H2, handed the same crystal the wild-type sequence does not fully follow it:
holds**, though not in the way that was meant. It does not follow it at all, and
neither does D67H, so the control has nothing to be compared with.

**H3, with templates found by search D67H stays wild-type: holds.** 0.52 Å from
1LZ1, closer than without any template.

Of the four outcomes the protocol listed for the last row, this is the third:
wild-type loops for both sequences. Its reading there: "The network cannot or will
not leave the wild-type conformation."

**The templates are not ignored. They are outvoted.**

- Handed 1LZ1, the loops go from 0.95 to 0.42 Å from that crystal (wild-type
  sequence; 0.82 to 0.39 for D67H), and the whole protein from 0.41 to 0.20 Å.
- Handed 1LYY, the loops shift by about 0.5 Å from where they are without a
  template, in the direction of 1LYY: 0.36 Å closer to it for the wild-type
  sequence, 0.41 Å for D67H. That is 6 to 7% of the distance between the crystals.
- Once the template is fixed, the sequence makes little difference: the wild-type
  and D67H predictions with the same template have loops 0.15 to 0.24 Å apart.

**Confidence drops where template and alignment disagree, and almost only there.**

| Mean pLDDT | No template | 1LZ1 | 1LYY |
|---|---|---|---|
| Wild type, the 13 residues | 96.1 | 97.6 | 94.0 |
| D67H, the 13 residues | 95.8 | 97.0 | 91.6 |
| Wild type, everything else | 98.2 | 98.4 | 98.3 |
| D67H, everything else | 98.1 | 98.4 | 98.0 |

Handed 1LYY, the lowest value anywhere in the protein is at residue 49 (67), in the
first of the two loops: 86.0 for the wild-type sequence and 84.7 for D67H, against
93 without a template. At the mutated residue itself (D67H sequence) pLDDT is 96.7
without a template, 97.4 with 1LZ1 and 93.1 with 1LYY: it is not one of the 13,
and it is the only residue outside them that loses more than 2 points. For the
wild-type sequence nothing outside the 13 loses more than 0.7. The drop is larger
for the D67H sequence than for the wild type, in the loops and at the lowest
residue.

### Interpretation

For the question asked: the prediction follows neither the sequence nor the
structure it was handed. It follows the alignment. Handed the answer, AlphaFold2
moved the loops by half an ångström and lowered its confidence in them.

This rules out one explanation of 002. The network did not miss the D67H loops for
lack of a hint about where they go: it was given their exact coordinates and kept
the wild-type ones.

By the protocol's own table the reading is "cannot or will not". This run cannot
tell those two apart, and the protocol named the reason under "Known before
running": the alignment is the same in every condition, and it describes the
wild-type family. So in the last row the network had three sources: 5,264
relatives saying wild type, one structure saying D67H, and a sequence differing
from wild type at one position. The relatives won, for both sequences. Whether the
D67H conformation is out of the network's reach, or only out of reach while the
full alignment is present, is the experiment that follows from this one.

That is not a far-fetched possibility. Published work on getting AlphaFold2 to
produce alternative conformations of a protein does it by making the alignment
shallower, with or without templates (del Alamo et al., *eLife* 2022, on
transporters and receptors; cited from memory, to be checked against the paper
before it is relied on).

Two smaller things are new:

1. **pLDDT pointed at the right place for once.** In 002 the D67H prediction came
   with pLDDT 98 and no warning. Here, when the handed structure and the alignment
   disagree, confidence drops in the residues they disagree about, by up to 8.5
   points, and at the mutated residue; nowhere else by more than 2. An 85 is
   still "confident" on the usual scale.
2. **The sequence does register, faintly.** The drop is larger when the sequence
   is D67H, the one that agrees with the template. One changed letter is not
   enough to change the answer; it is enough to make the network less sure of it.

### Limitations

- Every condition used the full alignment. This is the limitation that decides
  how far the result can be read (see above).
- Two networks and two seeds per cell. All four agree in every cell, but it is
  four.
- The two 1LYY cells pass the "wild-type loops" rule narrowly: the rule asks for
  under 1.5 Å from 1LZ1 and the worst prediction is at 1.44. The label depends on
  a threshold set in advance; "6 to 7% of the way to 1LYY" does not.
- The wild-type row without a template is not in `results/`; the positive control
  for it is.
- For the search condition, which four of the twenty hits the network read is not
  recorded, and their loop conformations were not measured. That they resemble
  1LZ1 is inferred from the predictions moving towards it.
- Search results depend on the server's database on the day (2026-11-11).
- As in 002: a crystal structure is one snapshot, and nothing here bears on why
  the mutation causes disease.

### Next

1. **Repeat the 1LYY row with less alignment**: once with a few dozen sequences
   and once with the sequence alone, for both sequences. If the loops move, the
   conformation is within reach and the alignment is what holds the prediction in
   place. If they do not move even then, "cannot" becomes the better reading.
   Sixteen predictions; needs one new setting in the pipeline (alignment depth per
   variant). The wild-type control matters more here than it did in this run: with
   little or no alignment the network may simply copy whatever it is handed.
2. **Experiment 005** is prepared and independent of this one.
3. **The language model on I56T and D67H**, carried over from experiment 003.
