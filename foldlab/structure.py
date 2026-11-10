"""Reading structures and comparing them. Everything here is plain geometry."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SVDSuperimposer import SVDSuperimposer

# residue number -> atom name -> xyz
Atoms = dict[int, dict[str, np.ndarray]]

PREDICTION_FILE = re.compile(r"_unrelaxed_rank_(\d+)_.*_model_(\d+)_seed_(\d+)\.pdb$")

# Largest possible solvent-exposed area per amino acid (Å², Tien et al. 2013).
MAX_ASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}


def _chain(path: Path, chain_id: str | None):
    model = PDBParser(QUIET=True).get_structure(path.stem, path)[0]
    chain = model[chain_id] if chain_id else next(iter(model))
    for residue in [r for r in chain if r.id[0] != " "]:  # drop waters and ligands
        chain.detach_child(residue.id)
    return chain


def load_atoms(path: Path, chain_id: str | None = None, offset: int = 0) -> Atoms:
    """Heavy atoms of one chain. `offset` shifts residue numbers (predictions start at 1)."""
    return {
        residue.id[1] + offset: {
            atom.get_name(): atom.coord.astype(float)
            for atom in residue
            if atom.element != "H"
        }
        for residue in _chain(path, chain_id)
    }


def relative_accessibility(path: Path, chain_id: str, sequence: dict[int, str], offset: int = 0) -> dict[int, float]:
    """Fraction of each residue's surface that water can touch: 0 = buried, 1 = fully exposed."""
    chain = _chain(path, chain_id)
    ShrakeRupley().compute(chain, level="R")
    return {
        residue.id[1] + offset: residue.sasa / MAX_ASA[sequence[residue.id[1] + offset]]
        for residue in chain
        if residue.id[1] + offset in sequence
    }


@dataclass
class Prediction:
    variant: str
    rank: int
    model: int
    seed: int
    path: Path
    atoms: Atoms
    plddt: dict[int, float]
    ptm: float

    @property
    def label(self) -> str:
        return f"{self.variant} model {self.model} seed {self.seed}"


def load_predictions(folder: Path, variant: str, offset: int) -> list[Prediction]:
    """Every prediction of one variant, whether ColabFold wrote them in one run or one run per network."""
    predictions = []
    for path in sorted(folder.rglob(f"{variant}_unrelaxed_rank_*.pdb")):
        match = PREDICTION_FILE.search(path.name)
        if not match:
            continue
        rank, model, seed = (int(g) for g in match.groups())
        scores = json.loads(Path(str(path).replace("_unrelaxed_", "_scores_")).with_suffix(".json").read_text())
        atoms = load_atoms(path, offset=offset)
        predictions.append(
            Prediction(
                variant, rank, model, seed, path, atoms,
                plddt={number: scores["plddt"][i] for i, number in enumerate(sorted(atoms))},
                ptm=scores["ptm"],
            )
        )
    if any(p.path.parent != folder for p in predictions):
        # A run split across processes only ranks within each process, so rank the whole set
        # here by mean confidence, best first: the rule ColabFold uses for single chains.
        predictions.sort(key=lambda p: -float(np.mean(list(p.plddt.values()))))
        for rank, prediction in enumerate(predictions, start=1):
            prediction.rank = rank
    return predictions


def ca(atoms: Atoms, residues: list[int]) -> np.ndarray:
    return np.array([atoms[r]["CA"] for r in residues])


def superpose(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Best rigid rotation + translation of `mobile` onto `target`; returns (rot, tran, RMSD in Å)."""
    sup = SVDSuperimposer()
    sup.set(target, mobile)
    sup.run()
    rot, tran = sup.get_rotran()
    return rot, tran, float(sup.get_rms())


def move(atoms: Atoms, rot: np.ndarray, tran: np.ndarray) -> Atoms:
    return {r: {name: xyz @ rot + tran for name, xyz in group.items()} for r, group in atoms.items()}


def rmsd_after_fit(mobile: Atoms, target: Atoms, fit: list[int], measure: list[int]) -> float:
    """Overlay `mobile` on `target` using the Cα of `fit`, then the Cα RMSD over `measure` (Å)."""
    rot, tran, _ = superpose(ca(mobile, fit), ca(target, fit))
    moved = ca(mobile, measure) @ rot + tran
    return float(np.sqrt(np.mean(np.sum((moved - ca(target, measure)) ** 2, axis=1))))


def neighbors(atoms: Atoms, sites: list[int], cutoff: float) -> list[int]:
    """Residues with any heavy atom within `cutoff` Å of any heavy atom of the site residues."""
    site_xyz = np.array([xyz for s in sites for xyz in atoms[s].values()])
    near = []
    for number, group in atoms.items():
        if number in sites:
            continue
        xyz = np.array(list(group.values()))
        if np.min(np.linalg.norm(xyz[:, None] - site_xyz[None], axis=-1)) <= cutoff:
            near.append(number)
    return sorted(near)


def side_chain_distance(atoms: Atoms, a: list[int], b: list[int]) -> float:
    """Closest approach between the side chains of two groups of residues, in Å.

    Backbone N, C and O are left out so that being next to each other in the chain does not
    count as being close; Cα is kept so that glycine, which has no side chain, still has a position.
    """
    pick = lambda group: np.array(
        [xyz for r in group for name, xyz in atoms[r].items() if name not in ("N", "C", "O")]
    )
    xa, xb = pick(a), pick(b)
    return float(np.min(np.linalg.norm(xa[:, None] - xb[None], axis=-1)))


def compare(a: Prediction, b: Prediction, core: list[int], sites: list[int], local: list[int]) -> dict:
    """Overlay b on a using the core Cα atoms, then measure how much moved and how confidence changed."""
    rot, tran, rmsd_core = superpose(ca(b.atoms, core), ca(a.atoms, core))
    moved = move(b.atoms, rot, tran)
    everyone = sorted(a.atoms)
    deviation = {r: float(np.linalg.norm(moved[r]["CA"] - a.atoms[r]["CA"])) for r in everyone}

    # Around the mutation, compare every atom the two structures share (side chains included).
    sq = [
        np.sum((moved[r][name] - xyz) ** 2)
        for r in local
        for name, xyz in a.atoms[r].items()
        if name in moved[r]
    ]
    mean = lambda p, residues: float(np.mean([p.plddt[r] for r in residues]))
    return {
        "rmsd_core": rmsd_core,
        "rmsd_local_all_atom": float(np.sqrt(np.mean(sq))),
        "max_ca_deviation": max(deviation.values()),
        "max_ca_deviation_residue": max(deviation, key=deviation.get),
        "delta_plddt_mean": mean(b, everyone) - mean(a, everyone),
        "delta_plddt_site": mean(b, sites) - mean(a, sites),
        "delta_plddt_local": mean(b, local) - mean(a, local),
        "deviation": deviation,
    }
