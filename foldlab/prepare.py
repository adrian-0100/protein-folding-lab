"""Step 1 - prepare: get the sequence, the experimental reference, and write one FASTA per variant."""

from __future__ import annotations

import requests
from Bio.PDB import PDBParser, PPBuilder

from .config import Experiment

UNIPROT = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
RCSB = "https://files.rcsb.org/download/{pdb}.pdb"


def _download(url: str) -> str:
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    return response.text


def structure_sequence(path, chain_id: str, offset: int) -> dict[int, str]:
    """Residue number -> one-letter amino acid, read from an experimental structure."""
    chain = PDBParser(QUIET=True).get_structure("ref", path)[0][chain_id]
    residues: dict[int, str] = {}
    for peptide in PPBuilder().build_peptides(chain):
        for residue, letter in zip(peptide, peptide.get_sequence()):
            residues[residue.id[1] + offset] = letter
    return residues


def reference_sequence(exp: Experiment) -> dict[int, str]:
    return structure_sequence(exp.reference_pdb, exp.cfg["reference"]["chain"], exp.reference_offset)


def run(exp: Experiment) -> None:
    protein = exp.cfg["protein"]
    exp.inputs.mkdir(parents=True, exist_ok=True)
    exp.reference_dir.mkdir(parents=True, exist_ok=True)

    # 1. The full UniProt entry, kept untouched as the record of where the sequence came from.
    accession = protein["uniprot"]
    entry = exp.inputs / f"uniprot_{accession}.fasta"
    if not entry.exists():
        entry.write_text(_download(UNIPROT.format(accession=accession)))
    lines = entry.read_text().splitlines()
    full = "".join(lines[1:])
    print(f"UniProt {accession}: {lines[0][1:]}")
    print(f"  full entry is {len(full)} residues")

    # 2. The construct we actually fold: a slice of that entry.
    start, end = exp.first_residue, exp.last_residue
    construct = full[start - 1 : end]
    print(f"  using residues {start}-{end} ({len(construct)} aa): {construct}")

    # 3. One FASTA per variant. The baseline is the slice; mutants are the slice with changes.
    for name, mutations in exp.variants.items():
        sequence = exp.apply(construct, mutations)
        (exp.inputs / f"{name}.fasta").write_text(f">{name}\n{sequence}\n")
        label = ", ".join(str(m) for m in mutations) or "no mutations (baseline)"
        print(f"  wrote inputs/{name}.fasta - {label}")

    # 4. The experimental structure, and a check that it really is the same protein.
    pdb = exp.cfg["reference"]["pdb"]
    if not exp.reference_pdb.exists():
        exp.reference_pdb.write_text(_download(RCSB.format(pdb=pdb)))
    observed = reference_sequence(exp)
    mismatches = [
        f"{number}: construct {construct[number - start]} vs {pdb} {letter}"
        for number, letter in observed.items()
        if start <= number <= end and construct[number - start] != letter
    ]
    covered = sum(1 for number in observed if start <= number <= end)
    print(f"Reference {pdb}: {covered}/{len(construct)} residues observed in the experiment")
    if mismatches:
        raise SystemExit(
            f"Reference {pdb} does not match the construct at:\n  " + "\n  ".join(mismatches)
        )
    print("  sequence and numbering match the construct")

    # 5. Experimental structures of the mutants themselves, when someone has solved them.
    for name, (path, chain) in exp.variant_references.items():
        if not path.exists():
            path.write_text(_download(RCSB.format(pdb=path.stem)))
        expected = exp.apply(construct, exp.variants[name])
        observed = structure_sequence(path, chain, exp.reference_offset)
        wrong = [n for n, letter in observed.items() if start <= n <= end and expected[n - start] != letter]
        if wrong:
            raise SystemExit(f"{path.stem} is not the {name} variant: it differs at {wrong}")
        print(f"Reference for {name}: {path.stem}, {len(observed)} residues observed, carries the mutation")
