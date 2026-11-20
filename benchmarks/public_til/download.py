#!/usr/bin/env python3
"""Download the public data used by the TIL prioritization benchmark.

    python benchmarks/public_til/download.py --data-dir /path/to/cache

About 1.4 GB. Nothing here needs controlled access:

* Caushi et al., Nature 2021 (NSCLC, anti-PD-1): GEO GSE176021 per-sample
  Cell Ranger matrices + ``filtered_contig_annotations.csv``; antigen labels
  from the paper's Supplementary Table 8 (Springer, open).
* Meng et al., Sci Transl Med 2023 (PDAC): GEO GSE254250 Cell Ranger h5
  matrices plus the CD8 metadata table carrying per-cell clonotype ids and
  tumor-reactive / non-reactive calls.
"""

from __future__ import annotations

import argparse
import gzip
import io
import tarfile
import urllib.request
from pathlib import Path

GEO = "https://ftp.ncbi.nlm.nih.gov/geo"

# Tumor (and metastasis) samples of the five patients whose MANA-specific
# clonotypes appear in tumor single-cell data. (GSM, sample name).
CAUSHI_SAMPLES = [
    ("GSM5352889", "MD01-004_tumor_1"),
    *[(f"GSM53528{n}", f"MD043-011_mettumor_{i}") for n, i in ((93, 1), (95, 2), (96, 3), (97, 4), (98, 5))],
    *[(f"GSM5352{897 + i}", f"MD043-011_tumor_{i}") for i in (2, 3, 4, 5)],
    *[(f"GSM53529{26 + i}", f"NY016-014_tumor_{i}") for i in (1, 2, 3)],
    *[(f"GSM53529{42 + i}", f"NY016-025_tumor_{i}") for i in (1, 2, 3, 4)],
    *[(f"GSM53529{58 + i}", f"MD01-005_tumor_{i}") for i in range(2, 10)],
]
# GEO misspells one archive's suffix.
VDJ_SUFFIX = {"NY016-025_tumor_3": "vdt"}
CAUSHI_LABELS = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03752-4/"
                 "MediaObjects/41586_2021_3752_MOESM8_ESM.xlsx")

MENG = f"{GEO}/series/GSE254nnn/GSE254250/suppl"
MENG_FILES = [f"GSE254250_{s}_GEX_filtered_feature_bc_matrix.h5"
              for s in ("249", "262", "282", "309", "418", "432", "Pool")]
MENG_FILES.append("GSE254250_PDAC_9sample_CD8_RNA_metadata_matrix.csv.gz")


def fetch(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"  {url}")
    request = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
    with urllib.request.urlopen(request) as response, open(tmp, "wb") as out:
        while chunk := response.read(1 << 20):
            out.write(chunk)
    tmp.rename(dest)
    return dest


def extract_flat(archive: Path, dest: Path, wanted: tuple[str, ...]) -> None:
    """Copy only the named regular files out of an archive, ignoring its paths."""
    dest.mkdir(parents=True, exist_ok=True)
    data = archive.read_bytes()
    while data[:2] == b"\x1f\x8b":  # some GEO VDJ archives are gzipped twice
        data = gzip.decompress(data)
    with tarfile.open(fileobj=io.BytesIO(data)) as tar:
        for member in tar:
            name = Path(member.name).name
            if member.isfile() and name in wanted and not (dest / name).exists():
                (dest / name).write_bytes(tar.extractfile(member).read())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True, type=Path)
    root = parser.parse_args().data_dir

    print("Caushi 2021 (GSE176021)")
    fetch(CAUSHI_LABELS, root / "caushi" / "supplementary_tables.xlsx")
    for gsm, name in CAUSHI_SAMPLES:
        base = f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_{name}"
        raw = root / "caushi" / "raw"
        gex = fetch(f"{base}.tar.gz", raw / f"{name}.tar.gz")
        vdj = fetch(f"{base}.{VDJ_SUFFIX.get(name, 'vdj')}.tar.gz", raw / f"{name}.vdj.tar.gz")
        sample = root / "caushi" / "samples" / name
        extract_flat(gex, sample / "gex" / "filtered_feature_bc_matrix",
                     ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"))
        extract_flat(vdj, sample / "vdj", ("filtered_contig_annotations.csv",))

    print("Meng 2023 (GSE254250)")
    for name in MENG_FILES:
        fetch(f"{MENG}/{name}", root / "meng" / name)
    print(f"done: {root}")


if __name__ == "__main__":
    main()
