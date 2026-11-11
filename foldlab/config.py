"""Load an experiment folder: its experiment.toml plus the paths every step uses."""

from __future__ import annotations

import json
import os
import re
import tomllib
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
WEIGHTS_DIR = PROJECT_ROOT / "weights"
# Downloads shared between experiments (not in git).
DATA_DIR = PROJECT_ROOT / "data"

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

    @property
    def selection_dir(self) -> Path:
        """Where steps that choose variants by rule (select, the subset draw) record what they chose."""
        return self.dir / "selection"

    # --- protein -----------------------------------------------------------
    @property
    def first_residue(self) -> int:
        return self.cfg["protein"]["residues"][0]

    @property
    def last_residue(self) -> int:
        return self.cfg["protein"]["residues"][1]

    @property
    def background(self) -> list[Mutation]:
        """Changes that turn the UniProt sequence into the protein actually used in the laboratory.

        They are part of the baseline, not something under study: every variant carries them.
        """
        return [Mutation.parse(m) for m in self.cfg["protein"].get("background", [])]

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
    def allowed_differences(self) -> list[Mutation]:
        """Declared differences between the construct and the reference structure (construct -> reference).

        For when the best available structure is of a close relative, such as human ubiquitin
        standing in for yeast ubiquitin.
        """
        return [Mutation.parse(m) for m in self.cfg["reference"].get("allowed_differences", [])]

    @property
    def replicates(self) -> list[Path]:
        """Other experimental structures of the unmutated protein: an experimental noise floor."""
        return [self.reference_dir / f"{pdb}.pdb" for pdb in self.cfg["reference"].get("replicates", [])]

    # --- variants ----------------------------------------------------------
    @cached_property
    def _generated(self) -> dict:
        """Variants written by a step of the pipeline (a rule or a random draw), not by hand."""
        path = self.selection_dir / "variants.json"
        return json.loads(path.read_text()) if path.exists() else {}

    @cached_property
    def _variant_tables(self) -> dict[str, dict]:
        """Every variant as {mutations: [...], templates: None | "search" | [pdb ids]}."""
        tables = {}
        for name, spec in self.cfg.get("variants", {}).items():
            spec = {"mutations": spec} if isinstance(spec, list) else dict(spec)
            tables[name] = {"mutations": spec.get("mutations", []), "templates": spec.get("templates")}
        for name, mutations in self._generated.get("variants", {}).items():
            tables[name] = {"mutations": mutations, "templates": None}
        return tables

    @property
    def variants(self) -> dict[str, list[Mutation]]:
        return {name: [Mutation.parse(m) for m in t["mutations"]] for name, t in self._variant_tables.items()}

    def templates(self, variant: str) -> None | str | list[str]:
        """What the network is handed as a starting hint: nothing, a database search, or chosen structures."""
        return self._variant_tables[variant]["templates"]

    @property
    def baseline(self) -> str:
        plain = [n for n, t in self._variant_tables.items() if not t["mutations"] and not t["templates"]]
        if len(plain) != 1:
            raise ValueError("Exactly one variant must have no mutations and no templates (the baseline)")
        return plain[0]

    @property
    def mutants(self) -> list[str]:
        """Every variant that is compared with the baseline."""
        return [name for name in self._variant_tables if name != self.baseline]

    @property
    def variant_references(self) -> dict[str, tuple[Path, str]]:
        """Experimental structures to score a variant against, where they exist: name -> (file, chain)."""
        references = {**self.cfg["reference"].get("variants", {}), **self._generated.get("references", {})}
        return {
            name: (self.reference_dir / f"{ref['pdb']}.pdb", ref["chain"])
            for name, ref in references.items()
        }

    @property
    def models(self) -> list[int]:
        """Which of AlphaFold2's five networks to run."""
        settings = self.cfg["prediction"]
        return list(settings["models"]) if "models" in settings else list(range(1, settings["num_models"] + 1))

    # --- sequences ---------------------------------------------------------
    def construct(self) -> str:
        """The baseline sequence: the chosen slice of the UniProt entry, with the background applied."""
        entry = self.inputs / f"uniprot_{self.cfg['protein']['uniprot']}.fasta"
        full = "".join(entry.read_text().splitlines()[1:])
        return self.apply(full[self.first_residue - 1 : self.last_residue], self.background)

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
