# Experiment 003: a tool built for the question

**Status: prepared, not run (2026-11-10).** The measurements are downloaded, the sample of 40 mutants is drawn, the alignment is built and the language model is installed on the remote machine. No score has been computed and no structure predicted. This is the protocol. Hypotheses and decision rules
are fixed here before any score is computed.

## Question

For ubiquitin there is a laboratory measurement of what almost every single
mutation does. Which computational score agrees with it: AlphaFold2's confidence,
or scores built from evolution alone?

## Why this, why now

Experiments 001 and 002 showed that AlphaFold2's output barely reacts to a single
mutation. That leaves the practical question open: if you want to know whether a
mutation is harmful, what *should* you run? Researchers use a different family of
tools for this, called variant-effect predictors. They do not predict a shape at
all. They estimate how surprising a substitution is, given what evolution has
tolerated in related proteins.

This experiment puts three such scores and AlphaFold2's confidence side by side
against measured data. It also closes a loose end from 001, which assumed that
L67A is harmful without evidence for it.

## Background

**The measurement.** Roscoe et al. (*J. Mol. Biol.* 2013) built yeast strains in
which the cell's only ubiquitin carried one chosen mutation, and measured how fast
each strain grew. This is a **deep mutational scan** (DMS): one experiment, more
than a thousand mutants. Slow growth means the mutation damages something
ubiquitin has to do. The data are distributed in ProteinGym, a public benchmark
collection, already cleaned and labelled.

Checked on download:

| | |
|---|---|
| Dataset | `RL40A_YEAST_Roscoe_2013`, in ProteinGym v1.3 substitutions (43 MB zip) |
| Mutants | 1,195 single mutants, positions 2-76 |
| Labels | 437 harmful, 758 tolerated (ProteinGym's cut-off, growth score below -0.2) |
| Includes | L67A and K63A, the mutation of experiment 001 and its proposed surface control |

**It is yeast ubiquitin.** Yeast and human ubiquitin differ at three positions:
19 (Ser in yeast, Pro in human), 24 (Asp, Glu) and 28 (Ser, Ala). Everything in
this experiment uses the yeast sequence (UniProt P0CH08, residues 1-76), because
that is the protein that was measured. The reference structure stays 1UBQ, the
human crystal structure, with those three differences declared.

**Growth is not stability.** A mutation can leave the protein perfectly folded and
still be harmful because it breaks a binding surface. The scan measures the sum of
everything. No score here is expected to match it perfectly, and that includes the
experiment itself: a second, independent scan of the same protein (Mavor et al.,
*eLife* 2016, run as a university class project) gives the ceiling on how well
anything can agree.

## The four scores

| | Score | What it knows | Cost |
|---|---|---|---|
| S0 | BLOSUM62 | Chemistry only: how often this kind of swap is tolerated in proteins in general. Knows nothing about ubiquitin | none |
| S1 | Alignment log-odds | How often the new amino acid appears at this position among ubiquitin's relatives, compared with the original. Computed from the same alignment the pipeline already builds | none |
| S2 | ESM-2 | A protein language model (650M parameters) trained to fill in hidden amino acids across millions of sequences. Score: how much less likely it finds the mutant residue than the original when that position is hidden | minutes on a laptop |
| S3 | AlphaFold2 | Change in pLDDT at the mutated residue, mutant minus wild type, averaged over the five networks. Also reported: change in mean pLDDT, and Cα shift | 2.5 h on the remote machine, for a subset |

S0 to S2 are computed for all 1,195 mutants. S3 needs a structure prediction per
mutant, so it is computed for a subset.

**The subset for S3.** Forty mutants: 20 drawn at random from the harmful ones and
20 from the tolerated ones, positions 2-72 only (the tail 73-76 is floppy and
low-confidence whatever the sequence), with a fixed random seed recorded in
`experiment.toml`. L67A and K63A are added as named cases and reported separately,
not counted in the statistics. The draw is made before any prediction is run.

## How agreement is measured

Two numbers, both descriptive, no significance tests.

- **Spearman correlation (ρ).** Rank the mutants by the score and by measured
  growth, and ask how similar the two rankings are. 1 is identical, 0 is unrelated.
- **AUC.** Pick one harmful and one tolerated mutant at random. AUC is the chance
  that the score puts them in the right order. 0.5 is a coin flip, 1.0 is perfect.

With 20 against 20, an AUC computed from the subset is itself noisy, by roughly
±0.09. This is why S0 to S2 are reported both on all 1,195 mutants and on the same
40, so that the comparison with S3 is like for like.

## Hypotheses

- **H1.** Evolution carries the signal: S1 separates harmful from tolerated
  mutations clearly better than chance.
- **H2.** The language model does at least as well as the plain alignment score.
- **H3.** AlphaFold2's confidence does not: on the subset, S3 is no better than
  the chemistry-only score S0.

Expectation: all three hold. Averaged over many proteins, the ProteinGym benchmark
reports ρ of roughly 0.4 to 0.5 for alignment-based and language-model scores. The
figure for this particular scan has deliberately not been looked up.

## Decision rules

1. H1 holds if S1 reaches AUC ≥ 0.70 on all mutants.
2. H2 holds if S2's ρ on all mutants is not lower than S1's by more than 0.05.
3. H3 holds if, on the subset, S3's AUC is below 0.65, or is not higher than S0's
   on the same subset. H3 fails if S3 reaches 0.75 or more.
   Between 0.65 and 0.75 and above S0 is reported as inconclusive at this sample
   size, with the recommendation to enlarge the subset.
4. The ceiling: ρ between the 2013 and 2016 scans on shared mutants is reported
   next to every score. No score is criticised for falling short of a number the
   experiment itself does not reach.

## Method

1. `prepare`: residues 1-76 of UniProt P0CH08; reference 1UBQ chain A, with the
   three declared differences; download and unpack the ProteinGym file.
2. `msa`: one search for the yeast wild type.
3. `score` (new step): S0, S1, S2 for every mutant in the scan.
4. `predict`: wild type and the 42 selected mutants, 5 networks, 1 seed, no
   templates, on the remote machine. 215 predictions.
5. `analyze`: S3 from the predictions; ρ and AUC for every score against the 2013
   scan, on all mutants and on the subset; the ceiling from the 2016 scan.
6. `report`: one scatter plot per score (score against measured growth), one bar
   chart of AUC with the ceiling marked, and a table of where L67A and K63A fall.

## What is in place

Settings: [experiment.toml](experiment.toml). The sample that was drawn, with the
seed and the method: [selection/variants.json](selection/variants.json).

| Piece | Where |
|---|---|
| Reference check that accepts declared differences | `foldlab/prepare.py` |
| Loader for a ProteinGym scan, and the random draw | `foldlab/truth.py` |
| The `score` step: BLOSUM62, alignment log-odds, ESM-2 | `foldlab/score.py` |
| Rank correlation and AUC against the scan, figures, report | `foldlab/agreement.py` |

The language model needs `torch` and `transformers`, kept optional so the base
pipeline stays light: `uv sync --extra plm`. They and the ESM-2 weights (2.4 GB)
are installed on the remote machine, where the model runs on the GPU.

Tested so far: the scoring and analysis code on made-up measurements, and the
language model on three lysozyme variants. Not yet run on this experiment's data.

To run it:

```bash
infra/remote-gpu-pc/remote run experiments/003-ubiquitin-variant-effect
```

## Cost

S0 to S2: minutes on the laptop. S3: 215 predictions of 76 residues, about 2.5
hours on the remote machine at ten at a time (from experiment 001's measured 42
seconds per prediction).

## Known before running

- **Growth is not stability**, as above. A score can be "wrong" about a mutation
  that is harmful for a reason no sequence model could know.
- **The language model may have seen ubiquitin.** It was trained on natural
  sequences, not on this measurement, so this is not leakage of the answer. But
  ubiquitin is one of the most-sequenced proteins there is, which favours S1 and S2.
- **One seed for S3.** The spread between the five networks is reported as the
  uncertainty of each S3 value. Experiment 001 measured seed noise on human
  ubiquitin at 0.3 pLDDT points or less.
- **A fair test of S3 needs the subset to be random.** That is why it is drawn by
  a seed, not picked by hand.

## What we will have at the end

A table that answers "what should I run to judge a mutation?", with AlphaFold2
included as one of the candidates, and a measured answer to whether L67A is
harmful in a living cell.
