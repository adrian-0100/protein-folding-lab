# Lab notebook: 003-ubiquitin-variant-effect

Numbers come from `results/`; what they mean is written here by a person.
Mutants are named as in the scan: original amino acid, position in ubiquitin,
new amino acid. Everything is yeast ubiquitin unless it says otherwise.

---

## 2026-11-11 · First run

**Goal.** Put four scores next to a laboratory measurement of 1,195 ubiquitin
mutants and see which agrees with it: three scores built from sequences alone, and
AlphaFold2's change in confidence.

**What was done.**

1. Question, hypotheses and decision rules written in [README.md](README.md) on
   2026-11-10, before any score existed. The sample of 40 mutants for AlphaFold2
   was drawn the same day, by seed.
2. `prepare`, `msa` (done on 2026-11-10): residues 1-76 of UniProt P0CH08;
   reference 1UBQ with three declared differences; the Roscoe 2013 scan from
   ProteinGym; one search for the wild type, 21,413 aligned sequences.
3. `score`: BLOSUM62, alignment log-odds and ESM-2 for all 1,195 mutants. Under a
   minute, the language model on the remote machine's GPU.
4. `predict` on the remote GPU PC: wild type and 42 mutants, 5 networks, 1 seed,
   215 predictions as 215 processes, ten at a time. 2 h 32 min in total, about 7
   minutes per prediction while ten run side by side. No process failed or stalled.
5. `analyze`, `report`.

**Deviation from plan.** None in the method. After the run, the bar chart
`results/auc_by_score.png` turned out to have its title cut off; the drawing code
was fixed and `report` rerun on the laptop. `REPORT.md` and `summary.json` came out
identical, so no number changed.

**A check that is not in the report.** The wild-type predictions were overlaid on
the crystal structure for this entry, with the pipeline's own overlay function:
0.39 Å (0.38-0.43) over residues 1-72, mean pLDDT 95.2, pTM 0.83. That is the same
as human ubiquitin in experiment 001 (0.38 Å, 95.3, 0.83). The predictor works on
the yeast protein, and the three differences from the human crystal do not show.

### Results

AUC: the chance that a score puts a harmful and a tolerated mutant in the right
order. ρ: how alike the two rankings are.

| Score | All 1,195 mutants: ρ | AUC | The 40 sampled: ρ | AUC |
|---|---|---|---|---|
| S0 chemistry only (BLOSUM62) | 0.29 | 0.64 | 0.39 | 0.68 |
| S1 alignment log-odds | 0.42 | 0.72 | 0.38 | 0.69 |
| S2 language model (ESM-2) | 0.63 | 0.84 | 0.71 | 0.84 |
| S3 AlphaFold2, change in pLDDT at the site | | | 0.41 | 0.71 |
| Also reported: AlphaFold2, change in mean pLDDT | | | 0.57 | 0.78 |
| Also reported: AlphaFold2, shift of the fold | | | 0.55 | 0.72 |
| Ceiling: the second laboratory scan (1,138 shared mutants) | 0.84 | 0.97 | | |

The sample behaves like the whole: for S0 to S2 the AUC on the 40 is within 0.05
of the AUC on all mutants, well inside the ±0.09 expected from a sample this size.

**H1, evolution carries the signal: holds, narrowly.** The alignment score reaches
AUC 0.72 on all mutants; the bar was 0.70. It beats chemistry alone (0.64), not by
much.

**H2, the language model does at least as well: holds, by a wide margin.** ESM-2
reaches ρ 0.63 against 0.42 for the alignment score, and AUC 0.84. The rule only
asked it not to be worse by more than 0.05. It is also above the 0.4 to 0.5 that
such scores average over many proteins, as expected for a protein this well
represented among known sequences.

**H3, AlphaFold2's confidence does not separate them: inconclusive.** On the 40
sampled mutants S3 reaches AUC 0.71. The rule needed it below 0.65, or not above
S0 on the same sample (0.68), for H3 to hold, and 0.75 or more for H3 to fail. It
landed in between. The expectation written in the protocol was that H3 would hold.
It did not.

**The two other AlphaFold2 numbers do no worse.** Change in mean pLDDT reaches AUC
0.78 and ρ 0.57, the best of anything from AlphaFold2 and above the alignment
score on the same 40. It was listed in the protocol as "also reported", not as the
score the rule is applied to, and it stays that way: choosing it now, after seeing
that it did better, would be picking the winner after the race.

**What the AlphaFold2 signal looks like up close.**

- The fold never moves. The largest shift among the 42 mutants is 0.34 Å (L8P);
  the typical one is 0.05 Å for tolerated mutants and 0.11 Å for harmful ones. Two
  of AlphaFold2's own networks differ by 0.14 Å on the wild type.
- The confidence at the mutated residue drops by 0.6 points (median) for tolerated
  mutants and by 1.9 for harmful ones, on a scale of 100. The five networks
  disagree with each other by about 0.5.
- Every mutant loses confidence, the harmless ones too: mean pLDDT falls in all
  42, by 0.17 at least. "No effect" is not zero on this scale. A likely reason,
  not investigated: the alignment was built around the wild type, so any change
  makes the sequence fit it a little worse.
- The separation comes from the large drops. Ten of the 40 lose more than 3
  points; eight of those are harmful (the two that are not: D52F, K33C). Five of
  the eight introduce a glycine or a proline: I61G (-21.5, with the networks far
  apart), I36G, L71P, G35P, K48P.
- Eight of the 20 harmful mutants lose less than 1 point, no more than a tolerated
  one: F4H, K6W, L8P, D32R, R42A, L43A, I44K, V70M. Two kinds of miss are in that
  list. Leu8, Ile44 and Val70 form the surface patch that most of ubiquitin's
  partners bind: I44K is one of the most harmful mutants in the whole scan
  (growth -0.97, 4th percentile) and the language model scores it in its worst 1%,
  while AlphaFold2's confidence there changes by -0.84. And L43A removes most of
  a fully buried leucine (growth -0.88) for a change of -0.97: the same kind of
  mutation as L67A in experiment 001, with the same result.
- Burial does not explain which harmful mutants AlphaFold2 reacts to. Measured on
  the crystal, the ones with large drops are if anything more exposed (median
  0.41) than the ones with small drops (0.28).

**Where the language model is wrong.** Its worst misses among harmful mutants are
instructive: K48R (growth -0.93) looks mild on every sequence scale, because
arginine is the closest thing to lysine. But Lys48 is where the chains that send
proteins to be degraded are built, and arginine cannot be linked. G76C (growth
-0.85) is the last residue, the glycine that gets attached to other proteins; the
model prefers cysteine there to the natural glycine, while the plain alignment
score puts every change at position 76 at the very bottom of its ranking. In the
other direction, 15 of the 18 measured
changes at Phe45 are tolerated in the laboratory, and the model scores them as
severe.

**The named cases.**

| | Measured growth | Label | S1 | S2 | S3 |
|---|---|---|---|---|---|
| L67A | -0.35 | harmful | -5.6 (11th percentile) | -7.8 (34th) | -0.48 |
| K63A | -0.05 | tolerated | -1.9 (88th) | 0.0 (97th) | -0.83 |

L67A is harmful in a living cell, mildly: the cut-off is -0.2 and the scan goes
down to -1.5. Experiment 001 assumed this without a measurement; now there is one.
Both evolution-based scores put the two mutants in the right order. AlphaFold2
puts them in the wrong order, and neither value stands out from the 0.6 that a
typical tolerated mutant loses. That agrees with experiment 001, which found
-0.19 for L67A on the human protein and called it noise.

### Interpretation

For the question asked: the scores built from evolution agree with the laboratory,
and the language model agrees best. It ranked 1,195 mutants in under a minute and
orders a harmful and a tolerated mutant correctly 84 times in 100. That is the
practical answer to "what do I run to judge a mutation?".

It is still far from the ceiling. Two laboratories measuring the same mutants
agree at 0.97. The misses show why some of the gap may never close from sequence
alone: K48R is harmful for what Lys48 does, not for what it does to the fold.

AlphaFold2 is the part that did not go as expected, and experiments 001 and 002
need a correction because of it. They concluded that AlphaFold2's output barely
reacts to a mutation. For one mutant at a time that still stands: the typical drop
for a harmful mutant (1.9 points) is the same order as the drop every mutant
shows (0.6) and as the disagreement between networks (0.5), so no single value
can be read. Across 40 mutants, though, the drops are larger for harmful mutants
more often than not. The confidence score carries some information about harm.
What it picks up is mostly glycines and prolines, the substitutions that change
what the backbone can do. It does not pick up a cavity in the core (L43A, L67A)
or a change on the binding surface (I44K).

What this run cannot say is how much information. At 0.71 with ±0.09, AlphaFold2's
confidence could be as weak as chemistry alone or nearly as good as the language
model. It cost 2.5 hours of a 20-core machine for 42 mutants, against under a
minute for 1,195.

### Limitations

- One protein and one kind of measurement. Growth of yeast sums up everything
  ubiquitin has to do; it is not a measurement of stability.
- Forty mutants for S3. H3 is undecided, and the difference between S3 and any
  other score on the sample is within the sampling noise, the language model
  included if taken strictly (0.84 against 0.71, with ±0.09 on each).
- One seed per network. The spread between networks (0.5 points) is as large as
  the typical effect for a tolerated mutant.
- The sample leaves out positions 73-76, by design. That tail is where almost
  every substitution is harmful, and AlphaFold2 was not tested there.
- The result for mean pLDDT (0.78) was noticed after the fact. It is a lead, not
  a finding.
- S1 is the simplest possible use of an alignment: one position at a time, on an
  alignment that includes ubiquitin-like domains. More careful alignment-based
  methods exist and were not run, so "the language model beats the alignment"
  holds for this alignment score only.
- The language model was trained on natural sequences, ubiquitin's among them,
  many times over. It did not see this measurement, but it may do less well on a
  protein with few known relatives.
- "Harmful" is ProteinGym's cut-off at -0.2. Mutants just below it (D58Y, D32I)
  are harmful by a hair.

### Next

1. **Settle H3 on a fresh sample.** Draw 40 new mutants with a new seed and state
   beforehand which AlphaFold2 number is being tested. This is the only fair way
   to follow the mean-pLDDT lead: the first sample suggested it, so the first
   sample cannot confirm it. About 2.5 hours on the remote machine. Adding 40 to
   the existing sample and keeping S3 would probably leave the verdict where it
   is: at 0.71, a sample of 80 still falls between the two thresholds.
2. **Ask the language model about experiment 002's mutants**, I56T and D67H of
   human lysozyme. Both cause disease and AlphaFold2 saw neither. Seconds to run.
3. **Harmful to growth, or to folding?** ProteinGym also holds direct stability
   measurements (Tsuboyama et al. 2023). The same four scores against one of
   those would say whether any score, AlphaFold2's included, does better when
   the measurement is about the fold and not about growth.
4. **Experiments 004 and 005** are prepared and independent of this one.
