# The next three experiments

A proposal. Nothing here has been run. Each experiment has its full protocol in
its own folder, written the way 001 and 002 were: question, hypotheses and decision
rules fixed before any result exists.

## Where we are

| | Asked | Found |
|---|---|---|
| [001](../experiments/001-ubiquitin-L67A/NOTEBOOK.md) | Does AlphaFold2 react to a mutation in ubiquitin's core? | Barely: same fold, same confidence |
| [002](../experiments/002-lysozyme-amyloid/NOTEBOOK.md) | Does it predict the loop movement that the D67H crystal of lysozyme shows? | No: it returns the wild-type structure, 5.8 Å off in the loops, with full confidence |

Both say what AlphaFold2 does *not* do. They leave three questions, and each of
the proposals answers one.

## The three proposals

| | Question | In one line | New data | New software | Compute |
|---|---|---|---|---|---|
| [003](../experiments/003-ubiquitin-variant-effect/README.md) | If not AlphaFold2, then what? | Score all 1,195 measured ubiquitin mutants with three evolution-based scores and with AlphaFold2's confidence, and see which agrees with the laboratory | ProteinGym scan of ubiquitin, 43 MB | `torch`, `transformers`, ESM-2 weights (2.6 GB) | minutes on the laptop, plus 2.5 h remote |
| [004](../experiments/004-lysozyme-templates/README.md) | Why did 002 fail? | Hand the network the D67H crystal as a starting hint and see whether the prediction moves, with the wild-type sequence as control | none | none | 30 min remote |
| [005](../experiments/005-t4-lysozyme-panel/README.md) | Was 002 typical? | A panel of about twelve T4 lysozyme mutants with crystal structures, chosen by rule: how often does the prediction move when the experiment does? | about 250 PDB files, 40 MB | none | about 3 h remote |

What each adds that the first two experiments did not have:

- **003** is the first with a *positive* result expected, and the first that uses a
  different kind of tool. It answers the question a student will actually have:
  what do I run to judge a mutation?
- **004** is the first that explains a result instead of only reporting one.
- **005** is the first with enough cases to state a rate, and the first with an
  experimental noise floor (the same protein crystallised nine times).

## Recommended order

1. **004 first.** Half an hour of computing, no new downloads, and the smallest
   code change. It also exercises two pipeline changes (variants as tables, choice
   of networks) that the other two can build on.
2. **003 second.** The most new code and the only new dependencies, but it is the
   one with the most useful answer, and most of it runs on the laptop.
3. **005 third, or leave it for the student.** It is the most work to run and the
   easiest to hand over: the selection step, the panel table and the final count
   are all things a sixth-semester student can do and defend. It is option 2 in
   [the student guide](student-project.md).

They are independent. Any one can be run without the others.

## Everything needed

### Already checked

| | Result |
|---|---|
| Ubiquitin scan in ProteinGym | Present: 1,195 single mutants, positions 2-76, labelled harmful or tolerated, includes L67A and K63A. A second independent scan of the same protein is there too |
| It is yeast ubiquitin | Yes. Differs from human at positions 19, 24 and 28. UniProt P0CH08, residues 1-76 |
| `torch` and `transformers` next to the current environment | Resolve cleanly (dry run: 31 packages) |
| Template support in ColabFold | Works with `.pdb` files given by hand; only networks 1 and 2 read templates; the search programs run on this kind of machine |
| T4 lysozyme mutant structures | 156 different single mutations with a crystal structure in one crystal form; 8 extra structures of the unmutated protein for a noise floor |
| Human lysozyme variants in UniProt | Only three: I56T and D67H (disease) and T70N (listed with no disease). Too few for a variant-effect test, which is why 003 uses ubiquitin |

### To download

| What | Size | For |
|---|---|---|
| ProteinGym v1.3 substitutions archive | 43 MB | 003 |
| `torch`, `transformers` and their dependencies | about 1 GB | 003 |
| ESM-2 weights, `facebook/esm2_t33_650M_UR50D` | about 2.6 GB | 003 |
| T4 lysozyme structures from the PDB | 40 MB | 005 |
| Template search results | small, fetched by ColabFold | 004 |

All free, all from the projects' own servers.

### To build in the pipeline

Listed once here; each protocol says where its pieces go.

| Change | Needed by | Size |
|---|---|---|
| Variants written as tables, so a variant can carry settings besides its mutations | 004 | small |
| Choice of networks (`models = [1, 2]`) | 004 | small |
| Templates: search, or chosen structures, with a private template folder per process | 004 | small |
| A loops verdict per prediction and a sequence-by-template table | 004 | small |
| Reference check that accepts declared sequence differences | 003 | small |
| Loader for a ProteinGym scan | 003 | small |
| A `score` step: BLOSUM62, alignment log-odds, ESM-2 | 003 | medium |
| A second analysis path: rank correlation and AUC against a scan | 003 | medium |
| `background` mutations: the protein in the crystals is not always the UniProt sequence | 005 | small |
| A `select` step: survey the PDB, apply a written rule, write the panel | 005 | medium (a working draft exists) |
| Experimental noise floor from replicate crystal structures | 005 | small |
| Panel summary across many mutants | 005 | medium |

The pipeline is 1,360 lines today. All three together would add roughly half that
again.

### Computing

| | Predictions | Residues | Estimate on the remote machine |
|---|---|---|---|
| 003 | 215 | 76 | 2.5 h (from 001's measured rate) |
| 004 | 28 | 130 | 0.5 h (from 002's measured rate) |
| 005 | about 130 | 164 | about 3 h (extrapolated, not measured) |

About six hours in total, during which the remote machine's processor is fully
used. Each run should be agreed with whoever else uses it, and is best started in
the evening. The parts of 003 that do not need structure predictions run on the
laptop.

### Decisions that are yours

1. **Which to run, and in what order.** The recommendation is above.
2. **Whether 005 is ours or the student's.** If the student's, we stop at the
   protocol and the pilot table that already exist.
3. **Adding `torch` and 3.6 GB of downloads** to the laptop for 003.
4. **When the remote machine can be used at full load.**

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

## How the three fit the student project

| Student-guide option | Becomes |
|---|---|
| 1. Does a tool built for the question do better? | Experiment 003, on ubiquitin, where the measured answer exists |
| 2. Is it the same with more mutants? | Experiment 005 |
| 3. Does letting the network look at known structures change the answer? | Experiment 004 |
| 4. Bring it to physiology | Still open, and now with three more worked protocols to copy from |
