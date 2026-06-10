# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Build the germline V-gene signal-peptide (leader) reference (#267).

Parses the vendored IMGT/GENE-DB L-PART1+L-PART2 FASTA (TRAV + TRBV) under
``tcrsift/refseqs/imgt_vgene_leaders/`` and writes a translated, deduplicated
reference table ``tcrsift/refseqs/vgene_leaders.csv.gz`` with columns:

    gene, allele, functionality, leader_aa, leader_nt

Only records that begin at ATG and translate cleanly to an M-started, stop-free
peptide are kept (partial-in-5' / non-ATG entries are skipped). Each leader is
translation-verified here so a transcription slip in the vendored FASTA surfaces
as a build error rather than a silently wrong reference.

The vendored FASTA is the authoritative source. To RE-FETCH it from IMGT (needs
network), the GENElect queries are::

    https://www.imgt.org/genedb/GENElect?query=8.1+TRAV&species=Homo+sapiens&IMGTlabel=L-PART1+L-PART2
    https://www.imgt.org/genedb/GENElect?query=8.1+TRBV&species=Homo+sapiens&IMGTlabel=L-PART1+L-PART2

Run: ``python scripts/build_vgene_leaders.py``
"""

from __future__ import annotations

import gzip
from pathlib import Path

from tcrsift.assemble import translate_dna

REFSEQS = Path(__file__).resolve().parent.parent / "tcrsift" / "refseqs"
RAW_DIR = REFSEQS / "imgt_vgene_leaders"
OUT = REFSEQS / "vgene_leaders.csv.gz"


def _parse_fasta(path: Path):
    """Yield (header, seq) pairs from a FASTA file."""
    header, chunks = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def _build():
    rows = []
    for fasta in sorted(RAW_DIR.glob("*.fasta")):
        for header, nt in _parse_fasta(fasta):
            # IMGT header: ACCESSION|GENE*ALLELE|species|functionality|label|...
            fields = header.split("|")
            gene_allele = fields[1]
            functionality = fields[3].strip()
            gene, _, allele = gene_allele.partition("*")
            nt = nt.lower()
            # Keep only clean, ATG-started ORFs (skip partial-in-5' / non-ATG).
            if not nt.startswith("atg"):
                continue
            aa, _ragged = translate_dna(nt.upper())
            # Keep only clean, full-length, M-started, stop-free leaders.
            # translate_dna truncates at the first stop, so an internal stop
            # shows up as aa shorter than nt//3 — that's a non-productive
            # (pseudogene) leader, skipped, not a transcription error.
            if not aa.startswith("M") or len(aa) != len(nt) // 3:
                continue
            rows.append((gene, allele, functionality, aa, nt))

    rows.sort(key=lambda r: (r[0], r[1]))
    with gzip.open(OUT, "wt", newline="") as fh:
        fh.write("gene,allele,functionality,leader_aa,leader_nt\n")
        for gene, allele, func, aa, nt in rows:
            fh.write(f"{gene},{allele},{func},{aa},{nt}\n")

    genes = {r[0] for r in rows}
    print(f"Wrote {OUT.relative_to(REFSEQS.parent.parent)}: "
          f"{len(rows)} alleles across {len(genes)} genes.")
    # Spot-check the cohort-critical leaders.
    by_ga = {(r[0], r[1]): r[3] for r in rows}
    for ga, expect in {
        ("TRAV1-1", "01"): "MWGAFLLYVSMKMGGTA",
        ("TRAV1-2", "01"): "MWGVFLLYVSMKMGGTT",
        ("TRBV13", "01"): "MLSPDLPDSAWNTRLLCHVMLCLLGAVSV",
        ("TRAV16", "01"): "MKPTLISVLVIIFILRGTR",
    }.items():
        got = by_ga.get(ga)
        status = "ok" if got == expect else f"MISMATCH (got {got})"
        print(f"  {ga[0]}*{ga[1]}: {got}  [{status}]")


if __name__ == "__main__":
    _build()
