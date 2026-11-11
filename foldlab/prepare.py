"""Step 1 - prepare: get the sequence, the experimental structures, and write one FASTA per variant."""

from __future__ import annotations

import shutil

import requests
from Bio.PDB import PDBParser, PPBuilder

from . import truth
from .config import Experiment

UNIPROT = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
RCSB = "https://files.rcsb.org/download/{pdb}.pdb"


def _download(url: str) -> str:
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response.text


def fetch_structure(path) -> None:
    """Download a PDB entry to `path` (named <id>.pdb) unless it is already there."""
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(_download(RCSB.format(pdb=path.stem)))


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


def differences(exp: Experiment, observed: dict[int, str], sequence: str) -> list[str]:
    """Where a structure's residues differ from `sequence`, written like mutations: S19P."""
    start, end = exp.first_residue, exp.last_residue
    return [
        f"{sequence[number - start]}{number}{letter}"
        for number, letter in sorted(observed.items())
        if start <= number <= end and sequence[number - start] != letter
    ]


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
    print(f"UniProt {accession}: {lines[0][1:]}")
    print(f"  full entry is {len(''.join(lines[1:]))} residues")

    # 2. The construct we actually fold: a slice of that entry, plus any laboratory background.
    start, end = exp.first_residue, exp.last_residue
    construct = exp.construct()
    print(f"  using residues {start}-{end} ({len(construct)} aa): {construct}")
    if exp.background:
        print("  with the laboratory background " + ", ".join(str(m) for m in exp.background))

    # 3. Laboratory measurements, and the sample of mutants drawn from them (variant-effect experiments).
    if "truth" in exp.cfg:
        truth.fetch(exp)
        truth.check_sequence(exp, construct)
        if "subset" in exp.cfg:
            truth.draw_subset(exp)
            exp = Experiment.load(exp.dir)  # the draw adds variants

    # 4. One FASTA per variant. The baseline is the construct; mutants are the construct with changes.
    for name, mutations in exp.variants.items():
        sequence = exp.apply(construct, mutations)
        (exp.inputs / f"{name}.fasta").write_text(f">{name}\n{sequence}\n")
    listed = list(exp.variants.items())
    for name, mutations in listed[:12]:
        label = ", ".join(str(m) for m in mutations) or "no mutations"
        hint = exp.templates(name)
        print(f"  wrote inputs/{name}.fasta - {label}" + (f"; template: {hint}" if hint else ""))
    if len(listed) > 12:
        print(f"  ... and {len(listed) - 12} more")

    # 5. The experimental structure, and a check that it really is the same protein.
    pdb = exp.cfg["reference"]["pdb"]
    fetch_structure(exp.reference_pdb)
    observed = reference_sequence(exp)
    found = differences(exp, observed, construct)
    declared = [str(m) for m in exp.allowed_differences]
    covered = sum(1 for number in observed if start <= number <= end)
    print(f"Reference {pdb}: {covered}/{len(construct)} residues observed in the experiment")
    if sorted(found) != sorted(declared):
        raise SystemExit(
            f"Reference {pdb} differs from the construct at {found or 'no position'}, "
            f"but the experiment declares {declared or 'no differences'}"
        )
    print("  sequence and numbering match the construct"
          + (f", apart from the declared differences {', '.join(declared)}" if declared else ""))

    # 6. Other crystal structures of the same unmutated protein: the experimental noise floor.
    for path in exp.replicates:
        fetch_structure(path)
        chain = exp.cfg["reference"]["chain"]
        wrong = differences(exp, structure_sequence(path, chain, exp.reference_offset), construct)
        if wrong:
            raise SystemExit(f"Replicate {path.stem} is not the same protein as the construct: {wrong}")
    if exp.replicates:
        print(f"Replicates: {len(exp.replicates)} more structures of the identical protein "
              f"({', '.join(p.stem for p in exp.replicates)})")

    # 7. Experimental structures to score each variant against. Normally the crystal of that
    #    variant; a control may deliberately be scored against another variant's crystal.
    sequences = {name: exp.apply(construct, mutations) for name, mutations in exp.variants.items()}
    by_sequence: dict[str, str] = {}
    for name, sequence in sequences.items():
        by_sequence.setdefault(sequence, name)
    for name, (path, chain) in exp.variant_references.items():
        fetch_structure(path)
        observed = structure_sequence(path, chain, exp.reference_offset)
        if not differences(exp, observed, sequences[name]):
            continue
        owner = next((other for seq, other in by_sequence.items() if not differences(exp, observed, seq)), None)
        if owner is None:
            raise SystemExit(f"{path.stem} matches no variant of this experiment (asked for by {name})")
        print(f"  note: {name} is scored against {path.stem}, which is the structure of {owner}")
    if exp.variant_references:
        used = sorted({path.stem for path, _ in exp.variant_references.values()})
        print(f"Variant references: {len(exp.variant_references)} variants scored against {len(used)} "
              f"structure(s): {', '.join(used[:8])}{' ...' if len(used) > 8 else ''}")

    # 8. Structures to hand to the network as templates.
    wanted = sorted({pdb for name in exp.variants for pdb in (exp.templates(name) or []) if exp.templates(name) != "search"})
    for pdb in wanted:
        fetch_structure(exp.dir / "templates" / f"{pdb}.pdb")
    if wanted:
        print(f"Templates: {', '.join(wanted)} in templates/")

    # 9. Work an earlier experiment already did: same alignments, so the comparison stays like for like.
    reuse = exp.cfg.get("reuse")
    if reuse:
        source = Experiment.load(exp.dir / reuse["from"])
        exp.msa.mkdir(exist_ok=True)
        baseline_a3m = exp.msa / f"{exp.baseline}.a3m"
        if not baseline_a3m.exists():
            shutil.copy(source.msa / f"{source.baseline}.a3m", baseline_a3m)
            print(f"Alignment: copied from {source.cfg['experiment']['id']} (msa/{exp.baseline}.a3m)")
