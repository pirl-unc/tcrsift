"""Regenerate packaged ``canonical_constants*.fasta`` from pyensembl.

Run as::

    python -m tcrsift._regen_canonical_constants [--species human|mouse] [--release 112] [--out PATH]

Fetches each canonical TCR-constant transcript's ``protein_sequence`` from
pyensembl, strips the leading ``X`` splice-boundary placeholder, and writes a
bare-mature FASTA.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .species import normalize_species, species_to_pyensembl

logger = logging.getLogger(__name__)

# Canonical transcripts and their UniProt cross-references.
CANONICAL_TRANSCRIPTS_BY_SPECIES: dict[str, dict[str, dict[str, str]]] = {
    "human": {
        "TRAC": {"ensembl": "ENST00000611116", "uniprot": "P01848"},
        "TRBC1": {"ensembl": "ENST00000633705", "uniprot": "P01850"},
        "TRBC2": {"ensembl": "ENST00000466254", "uniprot": "A0A5B9"},
    },
    "mouse": {
        "TRAC": {"ensembl": "ENSMUST00000198398", "uniprot": "P01849"},
        "TRBC1": {"ensembl": "ENSMUST00000192856", "uniprot": "P01852"},
        "TRBC2": {"ensembl": "ENSMUST00000103299", "uniprot": "A0A075B5J4"},
    },
}

# Backwards-compatible alias for human transcripts.
CANONICAL_TRANSCRIPTS: dict[str, dict[str, str]] = CANONICAL_TRANSCRIPTS_BY_SPECIES["human"]


def _resolve_protein_sequence(genome, tx) -> str | None:
    """Return ``tx.protein_sequence``, falling back to version-stripped
    ``protein_sequences.fasta_dictionary`` lookup for Ensembl >=112 where FASTA
    headers carry version suffixes (e.g. ``ENSMUSP00000143014.2``).
    """
    seq = tx.protein_sequence
    if seq:
        return seq
    pid = getattr(tx, "protein_id", None)
    if not pid:
        return None
    fdict = genome.protein_sequences.fasta_dictionary
    if pid in fdict:
        return fdict[pid]
    prefix = f"{pid}."
    for k, v in fdict.items():
        if k.startswith(prefix):
            return v
    return None


def fetch_canonical_constants(
    release: int = 110,
    *,
    species: str = "human",
) -> dict[str, str]:
    """Fetch canonical TCR constant-region AA sequences from pyensembl.

    Strips the leading ``X`` (splice-boundary partial codon) that pyensembl
    emits and returns the bare mature C-region for all three genes.
    """
    try:
        import pyensembl
    except ImportError as e:
        raise RuntimeError(
            "pyensembl is required to regenerate canonical_constants.fasta. "
            "Install with: pip install pyensembl"
        ) from e

    sp = normalize_species(species)
    transcripts = CANONICAL_TRANSCRIPTS_BY_SPECIES[sp]
    genome = pyensembl.EnsemblRelease(release, species=species_to_pyensembl(sp))
    genome.download()
    genome.index()

    out: dict[str, str] = {}
    for gene, ids in transcripts.items():
        tx = genome.transcript_by_id(ids["ensembl"])
        seq = _resolve_protein_sequence(genome, tx)
        if not seq:
            raise RuntimeError(
                f"No protein_sequence for {gene} ({ids['ensembl']}) in "
                f"{sp} Ensembl release {release}"
            )
        if not seq.startswith("X"):
            raise RuntimeError(
                f"Expected leading 'X' splice placeholder on {gene} "
                f"({ids['ensembl']}); got start {seq[:10]!r}. Check whether "
                f"Ensembl's convention changed in release {release}."
            )
        seq = seq[1:]
        # Human Ensembl ENST00000466254 encodes TRBC2*03 (K at mature pos 9,
        # 0-indexed 8), whereas ~99% of humans carry TRBC2*01 (E at pos 9,
        # NCBI AAA60662). Patch to TRBC2*01 so regenerating canonical_constants.fasta
        # never reverts human TRBC2*01 to *03 mislabeled as allele:01.
        if sp == "human" and gene == "TRBC2" and len(seq) > 8 and seq[8] == "K":
            seq = seq[:8] + "E" + seq[9:]
        out[gene] = seq
    return out


def write_fasta(
    constants: dict[str, str],
    path: Path,
    release: int,
    *,
    species: str = "human",
) -> None:
    """Write ``constants`` to ``path`` with provenance comments."""
    sp = normalize_species(species)
    transcripts = CANONICAL_TRANSCRIPTS_BY_SPECIES[sp]
    uniprot_summary = ", ".join(
        f"{info['uniprot']} ({gene})" for gene, info in transcripts.items()
    )
    lines = [
        f"# tcrsift canonical {sp} TCR constant-region AA sequences",
        f"# Source: Ensembl release {release} protein_sequence for each canonical transcript",
        f"# Verified against UniProt {uniprot_summary}",
        "#",
        "# All sequences are BARE MATURE protein (leading splice-boundary 'X' stripped).",
        "#",
        "# To regenerate this file:",
        f"#   python -m tcrsift._regen_canonical_constants --species {sp} --release {release}",
    ]
    for gene, seq in constants.items():
        ids = transcripts[gene]
        header = (
            f">{gene}|UniProt:{ids['uniprot']}|Ensembl:{ids['ensembl']}|"
            f"release:{release}|allele:01"
        )
        lines.append(header)
        for i in range(0, len(seq), 60):
            lines.append(seq[i : i + 60])
    path.write_text("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--species",
        choices=("human", "mouse"),
        default="human",
        help="Target species (default: human)",
    )
    parser.add_argument("--release", type=int, default=None, help="Ensembl release")
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output FASTA path",
    )
    args = parser.parse_args(argv)
    sp = normalize_species(args.species)
    release = args.release if args.release is not None else (110 if sp == "human" else 112)
    default_name = (
        "canonical_constants.fasta"
        if sp == "human"
        else "canonical_constants_mouse.fasta"
    )
    out_path = args.out or (Path(__file__).parent / "refseqs" / default_name)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    constants = fetch_canonical_constants(release=release, species=sp)
    write_fasta(constants, out_path, release=release, species=sp)
    logger.info("Wrote %s (%d genes)", out_path, len(constants))
    return 0


if __name__ == "__main__":
    sys.exit(main())
