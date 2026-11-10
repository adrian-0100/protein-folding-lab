"""Load an experiment folder: its experiment.toml plus the paths every step uses."""

from __future__ import annotations

import os
import re
import tomllib
from dataclasses import dataclass
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
WEIGHTS_DIR = PROJECT_ROOT / "weights"

# Environment for every ColabFold process we start. ColabFold draws a few plots with matplotlib,
# and on a Mac matplotlib's default is to open a real window. Over SSH, with no screen, that
# request can wait forever: four of fifteen processes hung that way in experiment 002.
# "Agg" draws straight to a file and never talks to the screen.
COLABFOLD_ENV = {**os.environ, "MPLBACKEND": "Agg"}

MUTATION = re.compile(r"^([A-Z])(\d+)([A-Z])$")


@dataclass(frozen=True)
class Mutation:
    original: str
    position: int  # numbered as in the UniProt entry
    new: str

    @classmethod
    def parse(cls, text: str) -> "Mutation":
        m = MUTATION.match(text)
        if not m:
            raise ValueError(f"Mutation {text!r} is not of the form L67A")
        return cls(m.group(1), int(m.group(2)), m.group(3))

    def __str__(self) -> str:
        return f"{self.original}{self.position}{self.new}"


@dataclass(frozen=True)
class Experiment:
    dir: Path
    cfg: dict

    @classmethod
    def load(cls, directory: str | Path) -> "Experiment":
        directory = Path(directory).resolve()
        with open(directory / "experiment.toml", "rb") as fh:
            return cls(directory, tomllib.load(fh))

    # --- folders -----------------------------------------------------------
    @property
    def inputs(self) -> Path:
        return self.dir / "inputs"

    @property
    def reference_dir(self) -> Path:
        return self.dir / "reference"

    @property
    def msa(self) -> Path:
        return self.dir / "msa"

    @property
    def predictions(self) -> Path:
        return self.dir / "predictions"

    @property
    def results(self) -> Path:
        return self.dir / "results"

    # --- protein -----------------------------------------------------------
    @property
    def first_residue(self) -> int:
        return self.cfg["protein"]["residues"][0]

    @property
    def last_residue(self) -> int:
        return self.cfg["protein"]["residues"][1]

    @property
    def reference_offset(self) -> int:
        """Added to the reference's residue numbers to get UniProt numbers (0 when they already agree).

        Needed when the deposited structure numbers the mature protein, e.g. after a signal
        peptide was cut off: mature residue 1 is then UniProt residue 1 + offset.
        """
        return self.cfg["reference"].get("numbering_offset", 0)

    @property
    def reference_pdb(self) -> Path:
        return self.reference_dir / f"{self.cfg['reference']['pdb']}.pdb"

    @property
    def variant_references(self) -> dict[str, tuple[Path, str]]:
        """Experimental structures of the mutants themselves, where they exist: name -> (file, chain)."""
        return {
            name: (self.reference_dir / f"{ref['pdb']}.pdb", ref["chain"])
            for name, ref in self.cfg["reference"].get("variants", {}).items()
        }

    # --- variants ----------------------------------------------------------
    @property
    def variants(self) -> dict[str, list[Mutation]]:
        return {
            name: [Mutation.parse(m) for m in muts]
            for name, muts in self.cfg["variants"].items()
        }

    @property
    def baseline(self) -> str:
        unmutated = [name for name, muts in self.variants.items() if not muts]
        if len(unmutated) != 1:
            raise ValueError("Exactly one variant must have no mutations (the baseline)")
        return unmutated[0]

    @property
    def mutants(self) -> list[str]:
        return [name for name, muts in self.variants.items() if muts]

    def sequence(self, variant: str) -> str:
        fasta = (self.inputs / f"{variant}.fasta").read_text().splitlines()
        return "".join(line.strip() for line in fasta if not line.startswith(">"))

    def apply(self, sequence: str, mutations: list[Mutation]) -> str:
        """Return `sequence` with the point mutations applied, checking each one."""
        residues = list(sequence)
        for mut in mutations:
            index = mut.position - self.first_residue
            if not 0 <= index < len(residues):
                raise ValueError(f"{mut}: position is outside the sequence")
            if residues[index] != mut.original:
                raise ValueError(
                    f"{mut}: position {mut.position} is {residues[index]}, not {mut.original}"
                )
            residues[index] = mut.new
        return "".join(residues)
