#!/usr/bin/env python3
"""Download the public data used by the TIL prioritization benchmark.

    python benchmarks/public_til/download.py --data-dir /path/to/cache

About 2.5 GB. Nothing here needs controlled access:

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
import re
import tarfile
import urllib.request
from pathlib import Path

GEO = "https://ftp.ncbi.nlm.nih.gov/geo"

# Every sample (tumor, metastasis, lymph node, adjacent normal lung) of the
# eight patients with antigen-labelled clonotypes in Supplementary Table 8.
CAUSHI_SERIES = f"{GEO}/series/GSE176nnn/GSE176021/suppl"
CAUSHI_PATIENTS = ("MD01-004", "MD01-005", "MD043-003", "MD043-008",
                   "MD043-011", "NY016-007", "NY016-014", "NY016-025")
ARCHIVE = re.compile(r"(GSM\d+)_(([A-Z0-9]+-\d+)_[A-Za-z]+_\d+)\.(vd[a-z]\.)?tar\.gz")
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
    listing = fetch(f"{CAUSHI_SERIES}/filelist.txt", root / "caushi" / "filelist.txt").read_text()
    archives = {}  # sample name -> {"gex": url, "vdj": url}; GEO misspells one ".vdj" suffix
    for gsm, name, patient, vdj in ARCHIVE.findall(listing):
        if patient in CAUSHI_PATIENTS:
            url = f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_{name}.{vdj}tar.gz"
            archives.setdefault(name, {})["vdj" if vdj else "gex"] = url
    raw = root / "caushi" / "raw"
    for name, urls in sorted(archives.items()):
        if set(urls) != {"gex", "vdj"}:
            continue  # prioritization needs both
        sample = root / "caushi" / "samples" / name
        extract_flat(fetch(urls["gex"], raw / f"{name}.tar.gz"), sample / "gex" / "filtered_feature_bc_matrix",
                     ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"))
        extract_flat(fetch(urls["vdj"], raw / f"{name}.vdj.tar.gz"), sample / "vdj",
                     ("filtered_contig_annotations.csv",))

    print("Meng 2023 (GSE254250)")
    for name in MENG_FILES:
        fetch(f"{MENG}/{name}", root / "meng" / name)
    print(f"done: {root}")


if __name__ == "__main__":
    main()
