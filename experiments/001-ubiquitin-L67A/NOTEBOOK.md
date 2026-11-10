# Lab notebook: 001-ubiquitin-L67A

One entry per working session. Numbers come from `results/`; what they mean is
written here by a person. Newest entry last.

---

## 2026-11-09 · Setup and first run

**Goal.** Get the pipeline working end to end on a protein whose answer is known,
and test whether AlphaFold2 responds to a single core mutation.

**What was done.**

1. Wrote the question, hypotheses and decision rule in [README.md](README.md)
   before running anything.
2. `prepare`: residues 1-76 of UniProt P0CG48; reference 1UBQ chain A. The
   reference sequence and numbering match the construct at all 76 positions.
3. `msa`: one MMseqs2 search (ColabFold server, UniRef + environmental) for the
   wild type: 21,577 aligned sequences, 20,566 unique. L67A uses the same file with
   only row 1 changed.
4. `predict`: AlphaFold2-ptm, 5 models x 2 seeds, 3 recycles, for each variant.
5. `analyze`, `report`.

**Where it ran.** The wild type was first predicted on the laptop (M3, 8 GB):
38 minutes for ten predictions. The laptop run of L67A was interrupted halfway.
Both variants were then run from scratch on the remote GPU PC (M1 Ultra, 128 GB)
through `infra/remote-gpu-pc/remote run`: 21 and 20 minutes. All results in `results/`
come from the remote run, so both variants share one machine. The laptop's
predictions are kept locally in `archive/001-laptop-run/` (not in the repository)
and were used for one extra check (below).

**Deviation from plan.** None in the method. The change of machine was not
planned; it is recorded here and in `results/provenance.json`.

### Results

**H1, positive control: holds.** The ten wild-type predictions overlay on the
crystal structure with a Cα RMSD of 0.38 Å (0.36-0.42) over residues 1-72, mean
pLDDT 95.3, pTM 0.83. Two crystal structures of the same protein typically differ
by about this much. The predictor works on this protein. The only low-confidence
region is the tail, residues 73-76 (pLDDT down to 48), which is flexible in
reality too.

**The site was a sensible choice.** In 1UBQ, leucine 67 has 0% of its surface
exposed to water: fully buried. Across relatives the position is leucine in 59%,
isoleucine in 22%, valine in 14%, phenylalanine in 2%: a large greasy side chain
in 97% of sequences. Alanine appears in 0.2%. Evolution almost never accepts this
substitution. (By plain conservation score it ranks 17th of 76, because the three
greasy residues count as different letters; "always hydrophobic" is the more
accurate description.) Its side chain comes within 3.8 Å of the Ile44 patch that
most of ubiquitin's partners bind.

**H2, the fold does not change: holds in substance, fails by the letter of the
rule.** Matched wild-type and L67A predictions differ by 0.085 Å (0.079-0.090 in
all ten pairs). Seed noise is 0.03 Å (0.02-0.03). So by the decision rule the
mutation does show up: the shift is about three times seed noise and perfectly
reproducible. It is also concentrated where it should be: the Cα of residue 67
moves 0.28 Å and its neighbours at 59-61 move 0.16-0.20 Å, against a seed noise
of 0.05 Å there.

But 0.085 Å is smaller than the difference between two of AlphaFold2's own five
networks (0.13 Å) and more than four times smaller than the prediction's distance
from the crystal structure (0.38 Å). No experiment could resolve it. The honest
statement is: *AlphaFold2 registers the mutation, and predicts the same fold.*

**H3, confidence drops at the site: fails.** pLDDT at residue 67 changes by -0.19
(median; 98.1 to 97.8 on a 0-100 scale), within seed noise (up to 0.32). Over the
whole protein: -0.10. One observation, not a conclusion: the change at the site
is negative in all ten matched pairs. That consistency suggests a real but
minuscule response; the pre-set rule compares magnitudes, and by that rule it is
noise.

**A ladder of differences**, all Cα RMSD over residues 1-72, which is the most
useful thing this run produced:

| Comparing | Typical difference |
|---|---|
| Same input, laptop vs remote | 0.009 Å |
| Same sequence, different random seed | 0.03 Å |
| **Wild type vs L67A, everything else equal** | **0.085 Å** |
| Same sequence, different AlphaFold2 network | 0.13 Å |
| Prediction vs crystal structure | 0.38 Å |

### Interpretation

AlphaFold2 predicts that ubiquitin L67A has the wild-type fold, with the same
confidence. Taken alone that is plausible: cavity-creating core mutations usually
do leave the folded structure nearly intact. What they change is stability, how
much of the time the protein spends folded, and nothing in this output measures
that. A pLDDT of 98 at a position where three of the four side-chain carbons
were just deleted shows that pLDDT is not a stability score.

So the answer to the question is: no meaningful change in predicted structure or
confidence, and this method could not have told us if the mutation were harmful.
That matches what Buel & Walters (2022) and Pak et al. (2023) reported on other
proteins.

### Limitations

- One mutation in one protein. Nothing here generalises.
- The decision rule tested whether an effect is *detectable*, not whether it is
  *large enough to matter*. It called 0.085 Å "above noise", which is true and
  unhelpful. The rule stays as written for this experiment; the next one adds a
  minimum effect size (see below).
- No relaxation step. Side chains around the new cavity are as the network left
  them. The all-atom RMSD near the site falls into two groups, a low one and one
  near 0.4 Å, in seed noise as well as in the mutation pairs. That looks like a
  neighbouring side chain flipping between two positions; not investigated.
- The alignment contains ubiquitin-like domains as well as ubiquitin, which makes
  every position look less conserved than it is within ubiquitin proper.
- Not verified here: that L67A is actually destabilising in ubiquitin. The
  expectation comes from the general behaviour of core mutations, not from a
  measurement I have read for this mutant.

### Next

1. **Amend the decision rule for future experiments**: an effect counts only if it
   exceeds seed noise *and* a stated minimum size (for RMSD, the model-to-model
   difference is a natural floor).
2. **Ask the question with a tool built for it.** Structure prediction answers
   "what shape"; a stability predictor or a protein language model answers "how
   much does this substitution hurt". Run one on L67A and on a surface control
   (for example K63A) and see whether it separates them.
3. **Choose a case with experimental truth for the mutant.** See
   `experiments/002-lysozyme-amyloid` and `docs/student-project.md`.
4. **The wet-lab experiment that would settle it**: express wild type and L67A,
   measure their melting temperatures. If L67A melts lower, the biochemical
   expectation was right and the prediction was blind to it.
