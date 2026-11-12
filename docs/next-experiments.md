# The next three experiments

**003 and 004 are run (2026-11-11). 005 is prepared, not run.** Each experiment
has its full protocol in its own folder, written the way 001 and 002 were:
question, hypotheses and decision rules fixed before any result exists. As of
2026-11-10 everything they need was built, fetched and checked on the remote
machine. What is left is to run 005.

What 003 found, in one line: the language model agrees best with the laboratory
(AUC 0.84 on all 1,195 mutants), the plain alignment score passes its bar (0.72),
and AlphaFold2's confidence came out undecided on the sample of 40 (0.71). See its
[notebook](../experiments/003-ubiquitin-variant-effect/NOTEBOOK.md).

What 004 found, in one line: handed the D67H crystal, AlphaFold2 kept the
wild-type loops for both sequences (6 to 7% of the way to the crystal it was
given) and lowered its confidence in them; whether that is "cannot" or "will not
while the full alignment is present" is still open. See its
[notebook](../experiments/004-lysozyme-templates/NOTEBOOK.md).

## Where we are

| | Asked | Found |
|---|---|---|
| [001](../experiments/001-ubiquitin-L67A/NOTEBOOK.md) | Does AlphaFold2 react to a mutation in ubiquitin's core? | Barely: same fold, same confidence |
| [002](../experiments/002-lysozyme-amyloid/NOTEBOOK.md) | Does it predict the loop movement that the D67H crystal of lysozyme shows? | No: it returns the wild-type structure, 5.8 Å off in the loops, with full confidence |

Both say what AlphaFold2 does *not* do. They leave three questions, and each of
the next experiments answers one.

## The three experiments

| | Question | In one line | Predictions to run | Estimate on the remote machine |
|---|---|---|---|---|
| [003](../experiments/003-ubiquitin-variant-effect/README.md) | If not AlphaFold2, then what? | Score all 1,195 measured ubiquitin mutants with three evolution-based scores, and a random 40 of them with AlphaFold2's confidence, and see which agrees with the laboratory | 215, of 76 residues | 2.5 h |
| [004](../experiments/004-lysozyme-templates/README.md) | Why did 002 fail? | Hand the network the D67H crystal as a starting hint and see whether the prediction moves, with the wild-type sequence as control | 28, of 130 residues | 0.5 h |
| [005](../experiments/005-t4-lysozyme-panel/README.md) | Was 002 typical? | Eight T4 lysozyme mutants with crystal structures, chosen by rule: how often does the prediction move when the experiment does? | 90, of 164 residues | about 2 h (extrapolated) |

What each adds that the first two experiments did not have:

- **003** is the first with a *positive* result expected, and the first that uses a
  different kind of tool. It answers the question a student will actually have:
  what do I run to judge a mutation?
- **004** is the first that explains a result instead of only reporting one.
- **005** is the first with enough cases to state a rate, and the first with an
  experimental noise floor (the same protein crystallised nine times).

## Running them

Each is one command, from the laptop. It sends the experiment to the remote
machine, runs every remaining step there, and brings the results back. The one
that remains:

```bash
infra/remote-gpu-pc/remote run experiments/005-t4-lysozyme-panel
```

003 and 004 were run the same way on 2026-11-11. 003 took 2 h 32 min, as
estimated. 004 took 36 minutes for its 28 new predictions: predictions with a
template are slower than the half hour allowed for.

005 is about two hours, during which the remote machine's processor is fully
used: agree the run with whoever else uses it, and prefer the evening.

After each run: read `results/REPORT.md`, then write the experiment's `NOTEBOOK.md`.
The report gives numbers and verdicts by the pre-set rules; what they mean is
written by a person.

## What was prepared

| | 003 | 004 | 005 |
|---|---|---|---|
| Sequence and reference structure checked | yeast ubiquitin against human 1UBQ, three declared differences | as in 002 | laboratory background WT\* against 1L63 |
| Laboratory data | 1,195 measured mutants, plus a second independent scan | crystal of D67H (1LYY) | 8 mutant crystals, 8 replicate crystals of the unmutated protein |
| Chosen by rule, before any prediction | 20 harmful + 20 tolerated mutants at random (seed 0), plus L67A and K63A | - | 4 movers and 4 paired non-movers out of 135 qualifying mutations |
| Alignment | new search, 21,413 sequences | copied from 002 | new search, 4,326 sequences |
| Extra software | ESM-2 language model, on the remote's GPU | none | none |

New steps in the pipeline: `select` (choose a panel by rule), `score` (sequence-based
scores), and `ready`, which runs everything up to but not including the first
computation. That is what was run for all three:

```bash
uv run python -m foldlab ready experiments/<name>
```

### What has been tested, and what has not

Tested: every new analysis path, on experiment 002's real predictions or on made-up
measurements; that ColabFold accepts a crystal structure file as a template; that
the language model loads and scores on the remote machine.

Not tested, because testing it is running the experiment: the panel summary on
the real panel. If it fails on the first run, nothing is lost: finished
predictions are kept and the run can be repeated. Two more things were on this
list until they ran, and both worked on the first run: the agreement analysis on
the real measurements (003) and an AlphaFold2 prediction with a template (004).

## How the three fit the student project

| Student-guide option | Is |
|---|---|
| 1. Does a tool built for the question do better? | Experiment 003, on ubiquitin, where the measured answer exists |
| 2. Is it the same with more mutants? | Experiment 005 |
| 3. Does letting the network look at known structures change the answer? | Experiment 004 |
| 4. Bring it to physiology | Still open, and now with three more worked protocols to copy from |

005 is the natural one to hand to a student: the rule, the panel table and the
final count are all things a sixth-semester student can run and defend, and it is
prepared to the point where running it is one command.

## Considered and left out for now

- **A structure predictor that runs on the Mac's GPU** (Apple's SimpleFold, ports
  of OpenFold3). It would speed everything up, but it is a different model with
  its own installation risk, and the first thing it would need is its own positive
  controls. Better as its own experiment once these three are in: "does a
  different predictor see D67H?"
- **Physics-based stability predictors** (FoldX, Rosetta). The classic way to
  estimate how much a mutation destabilises a protein. Both need separate licences
  and installations, and 003 answers the same practical question with free tools.
- **TEM-1 β-lactamase and sickle-cell haemoglobin**, from the student guide.
  Still good projects; they add proteins, not new kinds of evidence.
- **A stability scan in place of the growth scan for 003.** ProteinGym also holds
  direct stability measurements for hundreds of small domains (Tsuboyama et al.
  2023). None is of a protein already in this repository, so ubiquitin was kept;
  it is the natural follow-up if 003 raises the question "harmful to growth, or to
  folding?".
