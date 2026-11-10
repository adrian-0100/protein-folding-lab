# Experiment 004: what if AlphaFold2 is shown the answer?

**Status: proposed, not run.** This is the protocol. Hypotheses and decision rules
are fixed here before any prediction is made.

## Question

In experiment 002, AlphaFold2 predicted the wild-type structure for the D67H
mutant of lysozyme and missed a loop movement of up to 9.5 Å. If the network is
handed a known structure to start from, does the prediction follow the sequence or
the structure it was handed?

## Why this, why now

It is the cheapest experiment available and it explains the most. It reuses the
proteins, crystal structures and analysis of 002 and changes one setting. And it
separates two explanations of 002's result that predict different things:

- *The network cannot represent the D67H conformation at all.* Then even with the
  D67H crystal in hand, it will still return the wild-type loops.
- *The network can, but the sequence alone gives it no reason to.* Then with the
  D67H crystal in hand, it will move the loops.

## Background

**Templates.** AlphaFold2 can take, besides the sequence and its alignment, up to
four known structures of related proteins, called templates, as a starting hint.
Experiments 001 and 002 ran with templates off, the ColabFold default. With them
on, ColabFold searches a database of PDB structures and passes on the best matches.
It can also be given specific structures by hand.

**Only two of the five networks read templates.** AlphaFold2's networks 1 and 2
were trained with templates; 3, 4 and 5 ignore them (confirmed in the installed
code). So this experiment uses networks 1 and 2 only. That halves the predictions
and removes three networks that could not respond.

**The measuring stick already exists.** In the D67H crystal (1LYY), 13 residues sit
more than 2 Å from where they are in the wild-type crystal (1LZ1). Experiment 002
measured, over those residues, how far each prediction is from each crystal. The
D67H predictions were 0.79 Å from the wild-type crystal and 5.76 Å from their own.

## Design

Two sequences, four template conditions, networks 1 and 2, two seeds: four
predictions per cell.

| Template given | Wild-type sequence | D67H sequence | I56T sequence |
|---|---|---|---|
| None | from 002 | from 002 | from 002 |
| Whatever ColabFold's search finds | new | new | new |
| Wild-type crystal 1LZ1 only | new | new | |
| D67H crystal 1LYY only | new | new | |

The last row is the experiment. Its wild-type column is the control that makes it
interpretable: if the *wild-type* sequence also takes up the D67H loops when handed
1LYY, the network is copying the template and the sequence is not deciding
anything.

New predictions: 28, of 130 residues.

For the search condition, the templates ColabFold actually chose are recorded for
every run, since that is not under our control.

## Hypotheses

- **H1.** Handed the D67H crystal, the D67H sequence takes up the D67H loop
  conformation. Templates can steer the output.
- **H2.** Handed the same crystal, the wild-type sequence does *not* fully follow
  it, because its alignment still describes the wild-type protein.
- **H3.** With templates found by search, D67H stays wild-type, because nearly
  every lysozyme structure in the PDB has the wild-type loops.

Expectation: H1 and H3 hold. H2 is a real unknown. The four possible outcomes of
the last row:

| D67H sequence + 1LYY | Wild-type sequence + 1LYY | Reading |
|---|---|---|
| mutant loops | wild-type loops | The network uses sequence and template together. The strongest result |
| mutant loops | mutant loops | The network copies the template. The sequence is not deciding |
| wild-type loops | wild-type loops | The network cannot or will not leave the wild-type conformation |
| wild-type loops | mutant loops | Not expected. Would need explaining before anything else |

## Decision rules

Over the 13 residues that move in the crystal, after overlaying on the residues
that do not (the measure used in 002):

1. A prediction has **mutant loops** if it is under 1.5 Å from 1LYY and over 3 Å
   from 1LZ1.
2. It has **wild-type loops** if it is under 1.5 Å from 1LZ1 and over 3 Å from 1LYY.
3. Anything else is **in between**, and is reported with both distances.
4. A cell of the table gets a verdict only if all four of its predictions agree.

## Method

As experiment 002, same alignment files, same seeds. New in the pipeline: a
template setting per variant, and a choice of which networks to run.

## Planned `experiment.toml`

```toml
[experiment]
id = "004-lysozyme-templates"
question = "Handed a known structure to start from, does AlphaFold2's prediction follow the sequence or the structure?"

[protein]
name = "Lysozyme C (human)"
uniprot = "P61626"
residues = [19, 148]

[reference]
pdb = "1LZ1"
chain = "A"
numbering_offset = 18

[variants]
# name = mutations, and what the network is given as a template
wt            = { mutations = [] }
D67H          = { mutations = ["D85H"] }
I56T          = { mutations = ["I74T"] }
wt_search     = { mutations = [],       templates = "search" }
D67H_search   = { mutations = ["D85H"], templates = "search" }
I56T_search   = { mutations = ["I74T"], templates = "search" }
wt_1LZ1       = { mutations = [],       templates = ["1LZ1"] }
D67H_1LZ1     = { mutations = ["D85H"], templates = ["1LZ1"] }
wt_1LYY       = { mutations = [],       templates = ["1LYY"] }
D67H_1LYY     = { mutations = ["D85H"], templates = ["1LYY"] }

[reference.variants]
# every D67H condition is scored against the D67H crystal; wild-type ones too, for the control
D67H = { pdb = "1LYY", chain = "A" }
D67H_search = { pdb = "1LYY", chain = "A" }
D67H_1LZ1 = { pdb = "1LYY", chain = "A" }
D67H_1LYY = { pdb = "1LYY", chain = "A" }
wt_1LYY = { pdb = "1LYY", chain = "A" }

[prediction]
model_type = "alphafold2_ptm"
msa_mode = "mmseqs2_uniref_env"
models = [1, 2]             # the two networks that read templates
num_recycle = 3
num_seeds = 2
random_seed = 0
```

The three template-free variants are not recomputed: their predictions from
networks 1 and 2 are copied from experiment 002.

## What has to be built

| Piece | Where | Size |
|---|---|---|
| Variants written as tables, with an optional `templates` entry | `foldlab/config.py` | small |
| `models = [...]` in place of "the first N networks" | `foldlab/config.py`, `predict.py` | small |
| Passing templates to ColabFold: `--templates` for search, plus `--custom-template-path` for chosen structures | `foldlab/predict.py` | small |
| A private copy of the template folder for each process. ColabFold writes index files into it, and ten processes sharing one folder would collide | `foldlab/predict.py` | small |
| Recording which templates the search picked | `foldlab/predict.py`, `provenance.json` | small |
| A loops verdict per prediction, and the 2 by 2 table | `foldlab/analyze.py`, `report.py` | small |

No new software. Checked: ColabFold accepts `.pdb` files as custom templates, and
its template-search programs (`hhsearch`, `kalign`) are present and run on this
kind of machine.

## Cost

28 predictions of 130 residues. About half an hour on the remote machine; an
overnight run on the laptop would also do.

## Known before running

- **Handing over 1LYY is handing over the answer.** A "mutant loops" result in that
  cell does not show that AlphaFold2 can predict the D67H structure. It shows that
  the conformation is within its reach, which is the thing in doubt.
- **The alignment is unchanged**, so the network gets two sources that disagree
  when the sequence is wild type and the template is mutant, or the reverse. Which
  one wins is the measurement.
- **Two networks, two seeds.** Four predictions per cell is thin. Rule 4 is there
  so that a split result is reported as split.
- **Search results can change.** The template database behind ColabFold's server is
  updated from time to time, so that condition may not reproduce exactly later.

## What we will have at the end

A statement of *why* 002 came out as it did: either AlphaFold2 could not reach the
D67H structure, or it could and nothing in the sequence told it to.
