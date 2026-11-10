# Lab notebook: 002-lysozyme-amyloid

Numbers come from `results/`; what they mean is written here by a person.
Residue numbers are given as in the literature (mature enzyme), with the UniProt
number used by the pipeline in brackets where it helps: add 18.

---

## 2026-11-09 · First run

**Goal.** Test AlphaFold2 on two mutations whose real structures are known: one
that changes the crystal structure (D67H) and one that does not (I56T).

**What was done.**

1. Question, hypotheses and decision rules written in [README.md](README.md)
   before running. The decision rule was amended from experiment 001: an effect
   must exceed seed noise and the typical difference between two networks.
2. `prepare`: residues 19-148 of UniProt P61626; wild-type reference 1LZ1; mutant
   references 1LOZ (I56T) and 1LYY (D67H). The pipeline checked that each crystal
   carries its mutation and that the numbering lines up after the +18 shift.
3. `msa`: one search for the wild type: 5,264 aligned sequences, 5,066 unique.
4. `predict` on the remote GPU PC, 5 networks x 2 seeds per variant, 30 predictions,
   run as 15 processes, ten at a time.
5. `analyze`, `report`.

**What went wrong.** Four of the fifteen processes never started predicting. They
sat at 0% CPU for half an hour while the other eleven finished (in 13 to 16
minutes). Cause: ColabFold draws a plot of the alignment with matplotlib, which on
a Mac tries to open a real window. Over SSH with no screen that request can wait
forever, and with ten processes starting in the same second, four of them did. The
pipeline now tells matplotlib to draw to a file (`MPLBACKEND=Agg`). The four were
stopped and rerun alone (11 minutes); the eleven finished ones were kept. A
prediction does not depend on which other processes are running: experiment 001's
predictions came out bit-for-bit identical when rerun ten at a time.

**A detail worth knowing.** The mutants' alignments were built from the wild
type's by changing row 1 only. ColabFold then removes exact duplicates, and
reported one more usable sequence for each mutant (4,752 against 4,751). The most
likely reason: a database sequence identical to wild-type lysozyme is a duplicate
of the query in the wild-type run and a distinct relative in the mutant runs. The mutant runs
therefore see the wild-type sequence as one of 5,000 relatives. This is also what
would happen with a fresh search, and it is not a second variable in any practical
sense, but "only row 1 differs" is not literally the whole story.

### Results

**H1, positive control: holds.** Wild-type predictions sit 0.41 Å (0.39-0.42) from
the crystal structure 1LZ1, mean pLDDT 98.1, pTM 0.91.

**H2, I56T predicts no change of fold: holds.** Matched wild-type and I56T
predictions differ by 0.08 Å, detectable against seed noise (0.03) and smaller
than the difference between networks (0.12). The largest shift anywhere is 0.22 Å,
at the mutated residue. This agrees with the crystal, where I56T differs from wild
type by 0.28 Å. One metric did pass the rule: pLDDT at the site drops by 0.19, from
98.8 to 98.6. The rule calls it an effect because that number barely varies
between seeds or networks; it is too small to mean anything.

**H3, D67H predicts the loop movement: fails, by a wide margin.**

| Over the 13 residues that move more than 2 Å in the D67H crystal | |
|---|---|
| How far they move between the two crystals | 6.3 Å (2.1-9.5) |
| How far they move between wild-type and D67H predictions | 0.16 Å (0.06-0.57) |
| D67H predictions, distance from the **wild-type** crystal there | 0.79 Å |
| D67H predictions, distance from the **D67H** crystal there | 5.76 Å |

Over the whole protein, the D67H predictions are 0.38 Å from the wild-type crystal
and 1.83 Å from the D67H crystal. AlphaFold2 was given the D67H sequence and
returned the wild-type structure, with pLDDT 98. The figure
`results/experiment_vs_prediction_D67H.png` shows it at a glance: two experimental
peaks of 7.5 and 9.5 Å, and a predicted line that never leaves the floor.

Both decision conditions for H3 fail: the predicted shift is 2.5% of the
experimental one (the bar was half), and the predictions are closer to the wrong
crystal.

**But the network is not blind to D67H.** Around the mutated residue the predicted
atoms move by 0.61 Å (seed noise at most 0.08), and pLDDT there drops by 0.97, from
98.3 to 97.3. The largest predicted shifts are at residues 67 and 68 (0.72 and
0.57 Å), and both loops that move in the crystal show small bumps (0.12-0.29 Å):
clearly above seed noise in the second loop, only just in the first. So the response is in the right place and in the right direction
compared with I56T, which produced less on every local measure. It is about forty
times too small.

### Interpretation

For the question asked: no. AlphaFold2 does not predict the difference between
these two mutants that crystallography shows. For both it returns the wild-type
fold with full confidence.

Three things make this a stronger result than experiment 001:

1. There is experimental truth for the mutant, and the prediction can be scored
   against it: 5.8 Å off in the loops.
2. All three crystal structures were in AlphaFold2's training data (1984 and
   1997). Having seen the answer did not help.
3. The confidence score gave no warning. A reader shown only the D67H prediction
   and its pLDDT of 98 would have no reason to doubt it.

The faint local response is real and is the most interesting loose end: the
network's output does shift where the experiment shifts, only not nearly enough.

**The site evidence would have misled too.** By the two quick checks used to pick
a mutation, I56T looks damaging (fully buried, exposure 0.01; threonine in 0.1% of
relatives) and D67H looks harmless (on the surface, exposure 0.20; a variable
position where histidine appears in 5.7% of relatives). Both cause the same
disease, and it is the harmless-looking one that rearranges the structure.
Conservation and burial are a reasonable first filter and nothing more.

### Limitations

- Two mutations in one protein.
- No templates and no relaxation, as in 001. With templates switched on, the
  network could copy from a deposited structure; whether it would pick the mutant
  one is a separate question, not tested.
- A crystal structure is one snapshot, in one crystal form. The comparison takes
  the 1LYY loop positions as the truth for D67H.
- Crystal against crystal was overlaid using all residues, including the ones
  that move, which slightly understates how far the loops travel.
- Nothing here bears on why either mutation causes amyloid. The accepted
  explanation is reduced stability and partial unfolding, which a single predicted
  structure does not show.
- The first batch of predictions was interrupted (see above). The reruns used the
  same inputs, seeds and software.

### Next

1. **Ask a stability or variant-effect predictor** about I56T and D67H. Both
   reduce stability in reality; a tool built for that question should flag both.
2. **Try a predictor that uses the Apple GPU** (see `infra/remote-gpu-pc/README.md`). It
   is a different model, so the first thing to find out is whether it passes the
   same positive control, and then whether it handles D67H any differently.
3. **Rerun with templates enabled**, to see whether the output changes when the
   network is allowed to look at deposited structures.
4. **More mutations with known structures.** Two points do not make a trend; the
   PDB holds many lysozyme and T4 lysozyme mutants.
