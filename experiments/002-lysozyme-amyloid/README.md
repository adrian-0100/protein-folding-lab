# Experiment 002: the amyloid mutations of human lysozyme

Written before the predictions were run. Results are in [results/REPORT.md](results/REPORT.md),
conclusions in [NOTEBOOK.md](NOTEBOOK.md).

## Question

Two single mutations of human lysozyme cause the same disease. In the crystal, one
of them moves two loops by several ångströms and the other moves almost nothing.
Does AlphaFold2 predict that difference?

## Background

**Lysozyme** is a 130-residue enzyme in tears, saliva and macrophages that cuts the
cell wall of bacteria. Two inherited mutations, **I56T** and **D67H**, make it
deposit as amyloid fibrils in the kidneys, liver and spleen (hereditary systemic
amyloidosis; Pepys et al., *Nature* 1993).

Crystal structures exist for the wild type and for both mutants (Booth et al.,
*Nature* 1997). Measured with this repository's own code, overlaying each mutant
crystal on the wild-type crystal 1LZ1:

| Mutant crystal | Cα RMSD | Largest Cα shift | Residues shifting more than 2 Å |
|---|---|---|---|
| I56T (1LOZ) | 0.28 Å | 1.0 Å | none |
| D67H (1LYY) | 1.94 Å | 9.5 Å | 45-50, 68-73, 103 |

So there is an experimental answer for each mutant, and the two answers are
different. That is what experiment 001 lacked.

**Numbering.** The table and the literature count from the first residue of the
mature enzyme. UniProt counts from the precursor, which has an 18-residue signal
peptide in front. The pipeline and everything in `results/` use UniProt numbering:
add 18. I56T is I74T, D67H is D85H, and the loops are 63-68 and 86-91.

## Hypotheses

- **H1 (positive control).** AlphaFold2 reproduces wild-type lysozyme: Cα RMSD to
  1LZ1 below 1 Å, mean pLDDT above 90.
- **H2 (I56T).** The prediction shows no change of fold, as in the crystal.
- **H3 (D67H).** If AlphaFold2 captures what crystallography sees, the D67H
  predictions move the two loops, and over those residues they sit closer to the
  D67H crystal (1LYY) than to the wild-type crystal (1LZ1).

Expectation, from experiment 001 and from the literature it cited: H1 and H2 hold
and H3 fails, because the network returns close to the wild-type structure
whatever single residue is changed.

## Method

All settings are in [experiment.toml](experiment.toml). As in experiment 001:
residues 19-148 of UniProt P61626; one MSA search for the wild type, reused for
both mutants with only the query row changed; AlphaFold2-ptm through ColabFold,
5 networks x 2 seeds, 3 recycles, no templates, no relaxation. Thirty predictions
in total, run on the remote GPU PC.

New in this experiment: every mutant prediction is also compared with the
mutant's own crystal structure.

## Decision rules

Fixed here, before the run.

1. **Effect against noise** (amended after experiment 001). A metric shows an
   effect only if its median over matched pairs exceeds the largest seed noise,
   *and* exceeds the median difference between two of AlphaFold2's own networks.
   Above the first but not the second is reported as "detectable but small".
2. **H3 is supported** only if both hold:
   - the predicted shift of the residues that move in the D67H crystal is at
     least half of their experimental shift (medians), and
   - over those residues, the D67H predictions are closer to 1LYY than to 1LZ1.
3. **H2 is supported** if no fold-level metric for I56T passes rule 1.

## Known before running, and it matters

- **The network may have seen the answers.** AlphaFold2 was trained on PDB
  structures released up to 2018. 1LZ1, 1LOZ and 1LYY are from 1984 and 1997. If
  the D67H prediction does reproduce the loop movement, memorisation is a possible
  explanation. If it does not, that is the stronger result: even with the answer in
  its training data, a single changed residue was not enough to change the output.
- **A crystal is one snapshot.** The D67H loop positions could be influenced by how
  the molecules pack in that crystal. The original paper treats the movement as a
  property of the mutant; this experiment takes that as given.
- **Neither outcome explains the disease.** The accepted mechanism is that both
  mutants are less stable and partly unfold, and the amyloid forms from the
  unfolded part. A predicted folded structure, right or wrong, does not show that.
