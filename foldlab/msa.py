"""Step 2 - msa: search for evolutionary relatives once, then reuse that alignment for every variant.

The alignment (MSA) is AlphaFold2's main source of information. If each variant
got its own search, two things would differ between runs: the mutation and the
alignment. So the baseline is searched once and each mutant gets a copy of that
alignment where only the first row (the query itself) carries the mutation.
"""

from __future__ import annotations

import shutil
import subprocess

from .config import COLABFOLD_ENV, Experiment


def read_a3m(text: str) -> tuple[list[str], list[tuple[str, str]]]:
    """Split an .a3m file into its comment lines and its (header, sequence) records."""
    comments: list[str] = []
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith("#"):
            comments.append(line)
        elif line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line, []
        elif line.strip():
            chunks.append(line.strip().replace("\x00", ""))
    if header is not None:
        records.append((header, "".join(chunks)))
    return comments, records


def write_a3m(comments: list[str], records: list[tuple[str, str]]) -> str:
    lines = list(comments)
    for header, sequence in records:
        lines += [header, sequence]
    return "\n".join(lines) + "\n"


def run(exp: Experiment) -> None:
    exp.msa.mkdir(parents=True, exist_ok=True)
    baseline = exp.baseline
    baseline_a3m = exp.msa / f"{baseline}.a3m"

    if not baseline_a3m.exists():
        search_dir = exp.msa / "_search"
        print(f"Searching sequence databases for relatives of '{baseline}' (ColabFold MMseqs2 server)...")
        subprocess.run(
            [
                "colabfold_batch",
                "--msa-only",
                "--msa-mode",
                exp.cfg["prediction"]["msa_mode"],
                str(exp.inputs / f"{baseline}.fasta"),
                str(search_dir),
            ],
            check=True,
            env=COLABFOLD_ENV,
        )
        shutil.copy(search_dir / f"{baseline}.a3m", baseline_a3m)

    comments, records = read_a3m(baseline_a3m.read_text())
    query_header, query = records[0]
    expected = exp.sequence(baseline)
    if query != expected:
        raise SystemExit("The first row of the alignment is not the baseline sequence")
    print(f"msa/{baseline}.a3m: {len(records)} aligned sequences (row 1 is the query)")

    for name in exp.mutants:
        mutated = exp.apply(query, exp.variants[name])
        (exp.msa / f"{name}.a3m").write_text(
            write_a3m(comments, [(query_header, mutated)] + records[1:])
        )
        print(f"msa/{name}.a3m: same alignment, query row carries {name}")
