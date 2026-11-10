# Concepts

Short definitions, in the order you meet them when running an experiment.

## The molecule

**Amino acid / residue.** The building block of a protein. Twenty kinds, each with a
one-letter code (L = leucine, A = alanine). Once it is part of a chain it is called a
residue, and it is named by letter and position: L67 is the leucine at position 67.

**Backbone and side chain.** Every residue has the same backbone atoms, which link
into the chain, and a side chain that differs between the 20 kinds. Leucine's side
chain is a greasy four-carbon branch; alanine's is a single carbon.

**Cα (alpha carbon).** The backbone atom the side chain hangs from. One per residue,
so the trace of Cα atoms is the standard simplified picture of a protein's shape.

**Fold.** The 3D shape the chain settles into. Built from α-helices (spirals) and
β-strands (flat stretches that pair up into sheets).

**Hydrophobic core.** Greasy side chains avoid water, so they pack together in the
middle of the protein. This packing is the main force holding the fold together.
Residues there are "buried"; residues on the outside are "exposed".

**Wild type (WT).** The normal, unmutated sequence. The baseline of an experiment.

**Point mutation.** One residue swapped for another. Written `L67A`: leucine at
position 67 becomes alanine.

**Stability (ΔG) and ΔΔG.** A folded protein is constantly, briefly unfolding. ΔG
is how strongly it prefers the folded state; ΔΔG is how much a mutation changes
that. A destabilising mutation can leave the folded shape almost untouched while
making the protein spend more time unfolded, where the cell destroys it. Shape and
stability are different questions.

## The databases

**UniProt.** The curated catalogue of protein sequences. Every protein has an
accession (human polyubiquitin-C is `P0CG48`) with its sequence, function, known
variants and literature.

**PDB (Protein Data Bank).** The archive of experimentally measured structures.
Every entry has a four-character id (`1UBQ` is ubiquitin by X-ray crystallography).

**FASTA.** The plain-text format for sequences: a `>name` line, then the letters.

**Resolution (Å).** How fine the detail of an experimental structure is. An
ångström (Å) is 0.1 nanometres, about the size of an atom. Below 2 Å is good.

## The prediction

**Homologs.** Proteins in other species (or the same one) that descend from a
common ancestor of yours. Same job, similar sequence, same fold.

**MSA (multiple sequence alignment).** Your sequence and its homologs stacked so
that equivalent positions line up in columns. Stored here as `.a3m` files; row 1 is
your own sequence (the query).

**Conservation.** How little a column of the MSA varies. A position that is the
same in nearly every species is under strong selection: changing it usually breaks
something.

**Co-evolution.** When two columns change together across species, the two residues
are usually in contact: a change in one had to be compensated by the other. This is
the signal AlphaFold2 turns into 3D contacts.

**MMseqs2.** The fast search tool that finds the homologs. ColabFold sends the
query to a public MMseqs2 server instead of requiring terabytes of databases locally.

**AlphaFold2.** The neural network (DeepMind, 2021) that predicts atom coordinates
from a sequence and its MSA.

**ColabFold.** A faster, easier packaging of AlphaFold2: same network and weights,
MSA from the MMseqs2 server. What this project runs.

**Weights (parameters).** The numbers learned during training. AlphaFold2 ships
**five models**: five networks trained separately. Running all five gives five
slightly different answers, which is useful as a measure of uncertainty.

**Recycle.** The network feeds its own output back in and refines it. Three
recycles is the standard.

**Seed.** The starting value of the random number generator. It changes which
subset of a large MSA the network sees. Different seeds give slightly different
predictions of the same sequence; that spread is the **noise floor**.

**Rank.** ColabFold orders its predictions by confidence. `rank_001` is the one it
trusts most.

## Reading the output

**pLDDT.** Per-residue confidence, 0-100. Above 90: very high, side chains usually
right. 70-90: backbone right. 50-70: low. Below 50: often a region with no fixed
structure at all. It is the model's estimate of its own accuracy, not a
measurement of stability or flexibility. In the output PDB files it is stored in
the B-factor column.

**pTM.** A single 0-1 score for whether the overall fold is right. Above about 0.5
the fold is probably correct.

**PAE (predicted aligned error).** For every pair of residues, how uncertain their
relative position is. Matters for multi-domain proteins and complexes; not used in
experiment 001.

**Superposition.** Rotating and shifting one structure onto another so they overlap
as well as possible. Needed before any comparison because coordinates have no
absolute frame.

**RMSD.** Root-mean-square deviation: after superposition, the typical distance
between equivalent atoms, in Å. Two crystal structures of the same protein differ
by about 0.3-0.5 Å. Under 1 Å is "the same structure"; above 2-3 Å something moved.

**Relaxation.** A brief physics-based clean-up of the predicted coordinates to
remove atom clashes. Skipped here; it does not change the comparison.

## Experimental design

**Positive control.** A case where you already know the answer, run to prove the
method works before you trust it on the unknown.

**Matched pair.** The baseline and the mutant predicted by the same network with the
same seed. Everything is identical except the mutation, so any difference is caused
by the mutation.

**Noise floor.** How much a metric changes when nothing meaningful changes (here: a
different seed). A real effect has to be bigger than this.

**Decision rule.** The threshold for "yes, there is an effect", written before
looking at the data.
