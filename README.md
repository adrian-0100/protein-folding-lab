# protein-folding-lab

A small, reproducible setup for running structure-prediction experiments the way a
research group does, written so that someone new to the subject can follow what
happens at each step and why.

Two experiments are already run. Read this page first, then them, in order:

1. [001-ubiquitin-L67A](experiments/001-ubiquitin-L67A/README.md): does AlphaFold2
   react to a mutation in the core of a small protein? (Barely.)
2. [002-lysozyme-amyloid](experiments/002-lysozyme-amyloid/README.md): two disease
   mutations whose real structures are known. Does it predict the one that moves?
   (No.)

Three more are written up and ready to run:
[docs/next-experiments.md](docs/next-experiments.md).

## The idea in five minutes

A protein is a chain of amino acids. There are 20 kinds, each written as one letter,
so a protein can be written as text: `MQIFVKTLTGK...`. That text is its **sequence**.

In a cell the chain does not stay stretched out. It folds, within milliseconds to
seconds, into one specific 3D shape, and that shape is what lets the protein do its
job: bind a hormone, cut a sugar, carry oxygen. The same sequence folds into the
same shape every time, so in principle the shape can be computed from the sequence.
Doing that computation is the "protein folding problem".

Measuring a shape in the laboratory (X-ray crystallography, cryo-electron
microscopy, NMR) takes months and sometimes fails. **AlphaFold2** is a neural
network that predicts the shape from the sequence in minutes. Its main trick is
evolution: it looks up thousands of relatives of your protein in other species and
lines them up (a **multiple sequence alignment**, MSA). Positions that never change
are important; pairs of positions that change *together* are probably touching in
3D. From those clues it outputs the coordinates of every atom plus a confidence
score for every residue (**pLDDT**, 0-100).

Three things the AI does not do, and that the whole workflow is built around:

1. It does not tell you whether the protein *works*. It gives a shape.
2. It does not measure anything. A prediction is a hypothesis that an experiment
   could prove wrong.
3. It is not equally good at everything. It is excellent at "what does this natural
   protein look like" and, as both experiments here show, close to blind to "what
   does this one mutation do".

Every term in bold, and a few more, is explained in [docs/concepts.md](docs/concepts.md).

## The workflow researchers use, and where it lives here

| # | What a researcher does | Why | Here |
|---|---|---|---|
| 1 | Write a question that can be answered with a number | "Does mutation X change metric Y" can be wrong; "understand protein X" cannot | the experiment's `README.md` and `experiment.toml` |
| 2 | Get the sequence and what is already known | Start from the curated record, not from a copy-pasted string | `prepare` (UniProt + the experimental structure from the PDB) |
| 3 | Find evolutionary relatives | The alignment is the predictor's main input | `msa` (MMseqs2 search through the ColabFold server) |
| 4 | Predict the structure | | `predict` (AlphaFold2 through ColabFold, run locally) |
| 5 | Compare variants, against controls | A difference means nothing until you know how big noise is | `analyze` |
| 6 | Look, summarise, write down | | `report`, then `NOTEBOOK.md` by hand |
| 7 | Cross-check with physics, then with a real experiment | The prediction only proposes | outside this repo; each experiment's README says what the next step would be |

## Run it

Requirements: [uv](https://docs.astral.sh/uv/) and about 4 GB of disk for the
AlphaFold2 weights. No GPU needed for small proteins.

```bash
uv sync
```

```bash
uv run python -m foldlab all experiments/001-ubiquitin-L67A
```

`all` runs the five steps in order. On a machine with many cores, put
`FOLDLAB_JOBS=10` in front to run ten predictions at a time; the results are
identical, only sooner.

`all` is not required: each step can be run alone (`prepare`, `msa`,
`predict`, `analyze`, `report`), and each skips work that is already on disk, so
rerunning `analyze` and `report` after changing the analysis takes seconds.

Then open, inside the experiment folder:

- `results/REPORT.md` - the numbers and figures
- `results/viewer.html` - the structures in 3D, in your browser
- `NOTEBOOK.md` - what was concluded and what was not

## What is in the box

```
foldlab/                  the pipeline (about 700 lines of Python, meant to be read)
  prepare.py   step 1     sequence from UniProt, reference structure from the PDB, one FASTA per variant
  msa.py       step 2     one database search, one alignment per variant
  predict.py   step 3     AlphaFold2 with identical settings for every variant; records versions and settings
  analyze.py   step 4     confidence, structural difference, noise floor, conservation, burial
  report.py    step 5     figures, 3D viewer, REPORT.md
  structure.py            geometry helpers: read a structure, overlay two, measure distances
experiments/
  001-ubiquitin-L67A/
    experiment.toml       the whole experiment as data: protein, mutation, settings
    README.md             question, hypothesis, plan and decision rule, written before the run
    NOTEBOOK.md           lab notebook: what was done, what came out, what it means
    inputs/ reference/ msa/ predictions/ results/     produced by the pipeline
  002-lysozyme-amyloid/   the same layout; adds crystal structures of the mutants themselves
docs/
  concepts.md             glossary for someone who has never done this
  student-project.md      how this maps onto a sixth-semester biomedicine project
infra/remote-gpu-pc/      run on the remote GPU PC: guest account setup, the `remote` command
weights/                  AlphaFold2 parameters (3.7 GB, downloaded once, not in git)
```

## The research habits this setup enforces

These are what separate "I ran AlphaFold" from an experiment.

- **A positive control.** Before trusting the predictor on a mutant, check it on
  something whose answer is known: the baseline prediction is overlaid on the
  experimental structure.
- **A noise floor.** Each sequence is predicted several times (five networks, two
  random seeds). Two predictions of the *same* sequence are never identical, and a
  mutation only counts if it moves a metric more than that.
- **One variable.** Every variant gets the same alignment, model, seeds and
  settings. The mutant's alignment is the baseline's with only the query row changed.
- **A decision rule fixed in advance.** Written in the experiment README before the
  run, so the conclusion cannot be bent to fit the result.
- **Evidence for the choice of mutation.** Conservation and burial are computed and
  reported, not asserted.
- **Provenance.** Software versions, weights, settings and input checksums are saved
  with the results.
- **A notebook.** Numbers are generated; conclusions are written by a person, next
  to their limitations.

## Starting your own experiment

1. Copy `experiments/001-ubiquitin-L67A/experiment.toml` into a new folder under `experiments/`.
2. Edit the protein (UniProt accession, residue range), the reference structure
   (PDB id, chain), and the mutation.
3. Write the experiment `README.md` (question, hypothesis, decision rule) *before* running.
4. `uv run python -m foldlab all experiments/<your folder>`.

`prepare` stops with a clear message if the mutation does not match the sequence or
if the reference structure is a different protein or numbered differently (set
`numbering_offset` under `[reference]` when the deposited structure numbers the
mature protein and UniProt numbers the precursor).

## What this machine can do

Predictions run on the CPU, and that is slow: ten predictions of ubiquitin (76
residues) took 38 minutes on this laptop (M3, 8 GB). Time and memory grow faster
than linearly with length, so anything much longer than 150 residues, and any
complex, belongs on a GPU (Google Colab's free ColabFold notebook runs the same
software).

There is also a guest account on a remote GPU PC (a Mac: M1 Ultra, 20 cores, 128 GB).
It is still CPU-only: AlphaFold2 has no Apple GPU backend. One prediction keeps about 6
of its 20 cores busy, so the way to use the machine is to run several predictions
side by side, which the pipeline does there by default. One command sends an
experiment over and brings the results back:

```bash
infra/remote-gpu-pc/remote run experiments/001-ubiquitin-L67A
```

Setting it up, sharing access with other people, and what the account can and
cannot touch are in [infra/remote-gpu-pc/README.md](infra/remote-gpu-pc/README.md).

The `msa` step sends the baseline sequence to the public ColabFold MMseqs2 server.
That is fine for published proteins; do not send unpublished sequences there.

## Credit

If results from this pipeline are used in a report, cite ColabFold (Mirdita et al.,
*Nature Methods* 2022), AlphaFold2 (Jumper et al., *Nature* 2021) and MMseqs2
(Steinegger & Söding, *Nature Biotechnology* 2017). ColabFold writes the exact
references into each prediction folder as `cite.bibtex`.
