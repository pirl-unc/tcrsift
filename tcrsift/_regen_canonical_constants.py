"""Regenerate ``tcrsift/data/canonical_constants.fasta`` from pyensembl.

Run as::

    python -m tcrsift._regen_canonical_constants [--release 110] [--out PATH]

Fetches each canonical TCR-constant transcript's ``protein_sequence``
from pyensembl, strips the leading ``X`` splice-boundary placeholder,
and writes a FASTA with the prepend-E convention for β chains.

The output is the source of truth for :data:`tcrsift.assemble.HUMAN_TRAC_AA`,
:data:`tcrsift.assemble.HUMAN_TRBC1_AA`, and :data:`tcrsift.assemble.HUMAN_TRBC2_AA`
at runtime. See #100 for the motivation; previous hand-typed strings drifted
from UniProt at seven positions across the three genes.

CI runs the same fetch and asserts the packaged FASTA hasn't moved
(see ``tests/test_assemble_constants.py::test_canonical_matches_pyensembl``).
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

logger = logging.getLogger(__name__)

# Canonical transcripts and their UniProt cross-references. Stable across
# Ensembl releases at least through 114 (last verified 2026-05-28).
CANONICAL_TRANSCRIPTS: dict[str, dict[str, str]] = {
    "TRAC":  {"ensembl": "ENST00000611116", "uniprot": "P01848"},
    "TRBC1": {"ensembl": "ENST00000633705", "uniprot": "P01850"},
    "TRBC2": {"ensembl": "ENST00000466254", "uniprot": "A0A5B9"},
}


def fetch_canonical_constants(release: int = 110) -> dict[str, str]:
    """Fetch canonical TCR constant-region AA sequences from pyensembl.

    Strips the leading ``X`` (splice-boundary partial codon) that pyensembl
    emits and returns the bare mature C-region for all three genes.
    The J→C junction residue (universally E for β chains, J-dependent for
    α chains) is no longer baked into the canonical sequence; it's read
    per-clone from the CellRanger contig at assembly time (see #105).
    Pre-2.0 versions prepended a synthetic E to β chains.

    Parameters
    ----------
    release : int
        Ensembl release to query. Defaults to 110 — matches the rest of
        the tcrsift CLI's pyensembl-release default.

    Returns
    -------
    dict[str, str]
        ``{"TRAC": ..., "TRBC1": ..., "TRBC2": ...}``.
    """
    try:
        import pyensembl
    except ImportError as exc:
        raise RuntimeError(
            "pyensembl is required for regenerating canonical constants. "
            "Reinstall tcrsift with its core dependencies."
        ) from exc
    ens = pyensembl.EnsemblRelease(release)
    out: dict[str, str] = {}
    for gene, ids in CANONICAL_TRANSCRIPTS.items():
        tx = ens.transcript_by_id(ids["ensembl"])
        prot = tx.protein_sequence
        if prot is None:
            raise RuntimeError(
                f"pyensembl returned no protein_sequence for {gene} ({ids['ensembl']})"
            )
        if prot.startswith("X"):
            prot = prot[1:]
        out[gene] = prot
    return out


def write_fasta(constants: dict[str, str], path: Path, release: int) -> None:
    """Write the canonical-constants FASTA with provenance headers."""
    lines = [
        "# tcrsift canonical TCR constant-region AA sequences",
        f"# Source: Ensembl release {release} protein_sequence for each canonical transcript",
        "# Verified against UniProt P01848 (TRAC), P01850 (TRBC1), A0A5B9 (TRBC2)",
        "#",
        "# pyensembl's protein_sequence emits a leading 'X' for the J-C splice-boundary",
        "# partial codon. We strip that X and instead retain tcrsift's prepend-E convention",
        "# for β-chain C-regions (the splice junction residue is E for the canonical",
        "# TRBJ-TRBC1/2 pairs). α-chain stores the bare mature sequence (no prepend).",
        "#",
        "# To regenerate this file:",
        "#   python -m tcrsift._regen_canonical_constants",
    ]
    for gene, seq in constants.items():
        ids = CANONICAL_TRANSCRIPTS[gene]
        header = (
            f">{gene}|UniProt:{ids['uniprot']}|Ensembl:{ids['ensembl']}|"
            f"release:{release}|prepended-E:{str(gene.startswith('TRBC')).lower()}"
        )
        lines.append(header)
        # 60-char-wrapped FASTA body
        for i in range(0, len(seq), 60):
            lines.append(seq[i : i + 60])
    path.write_text("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=int, default=110, help="Ensembl release")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "refseqs" / "canonical_constants.fasta",
        help="Output FASTA path",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    constants = fetch_canonical_constants(release=args.release)
    write_fasta(constants, args.out, release=args.release)
    logger.info("Wrote %s (%d genes)", args.out, len(constants))
    return 0


if __name__ == "__main__":
    sys.exit(main())
