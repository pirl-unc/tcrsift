#!/usr/bin/env python3
"""Download the public data used by the TIL prioritization benchmark.

    python benchmarks/public_til/download.py --data-dir /path/to/cache

About 3 GB. Nothing here needs controlled access:

* Caushi et al., Nature 2021 (NSCLC, anti-PD-1): GEO GSE176021 per-sample
  Cell Ranger matrices + ``filtered_contig_annotations.csv``; antigen labels
  from the paper's Supplementary Table 8 (Springer, open).
* Eberhardt et al., Nature 2021 (HPV+ head and neck cancer): GEO GSE180268
  matrices and contig tables for HPV-tetramer-sorted and bulk PD-1+ CD8 TIL.
* Meng et al., Sci Transl Med 2023 (PDAC): GEO GSE254250 Cell Ranger h5
  matrices plus the CD8 metadata table carrying per-cell clonotype ids and
  tumor-reactive / non-reactive calls.
* Chiffelle et al., Immunity 2024 (melanoma, TIL therapy): baseline-tumor
  count tables and contig tables from GEO GSE222448 (deposited with Barras et
  al. 2024); tumor-reactive / non-reactive calls from the paper's Table S3.
* Pulliam et al., Cell Rep Med 2024 (Merkel cell carcinoma): GEO GSE227708
  runs Z1249 and Z1308, the two with labelled clones in a tumor sample whose
  hashtags match the deposited sample sheet; labels from the paper's Table S4.
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import io
import re
import socket
import tarfile
import urllib.error
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

# Eberhardt 2021: the three patients with a bulk PD-1+ CD8 TIL sample next to
# their HPV-tetramer-sorted samples.
EBERHARDT_SERIES = f"{GEO}/series/GSE180nnn/GSE180268/suppl"
EBERHARDT_PATIENTS = ("HPV7", "HPV34", "HPV51")
EBERHARDT_FILE = re.compile(
    r"(GSM\d+)_((HPV\d+)_\w+?)_(barcodes\.tsv|features\.tsv|matrix\.mtx|all_contig_annotations\.csv)\.gz")

MENG = f"{GEO}/series/GSE254nnn/GSE254250/suppl"
MENG_FILES = [f"GSE254250_{s}_GEX_filtered_feature_bc_matrix.h5"
              for s in ("249", "262", "282", "309", "418", "432", "Pool")]
MENG_FILES.append("GSE254250_PDAC_9sample_CD8_RNA_metadata_matrix.csv.gz")

# Chiffelle 2024: 5' CD45+ baseline-tumor (T0) libraries of the 11 patients
# with tested clonotypes. GEO holds dense gene-by-cell count tables.
CHIFFELLE_SERIES = f"{GEO}/series/GSE222nnn/GSE222448/suppl"
CHIFFELLE_PATIENTS = (1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13)
CHIFFELLE_COUNTS = re.compile(r"(GSM\d+)_patient(\d+)_(T0\d?)_CD45-counts\.tsv\.gz")
CHIFFELLE_CONTIGS = re.compile(r"(GSM\d+)_all_contig_annotations_Patient(\d+)_(T0\d?)\.csv\.gz")
CHIFFELLE_LABELS = "https://ars.els-cdn.com/content/image/1-s2.0-S1074761324004138-mmc3.xlsx"

PULLIAM_SERIES = f"{GEO}/series/GSE227nnn/GSE227708/suppl"
PULLIAM_RUNS = ("Z1249", "Z1308")
PULLIAM_FILE = re.compile(
    r"(GSM\d+)_(\w+?)_(barcodes\.tsv|features\.tsv|matrix\.mtx|all_contig_annotations\.csv)\.gz")
PULLIAM_LABELS = "https://ars.els-cdn.com/content/image/1-s2.0-S2666379124000351-mmc5.xlsx"


@contextlib.contextmanager
def _ipv4_only():
    """Force AF_INET resolution when a dual-stack mirror (e.g. NCBI GEO IPv6) returns HTTP 403."""
    orig = socket.getaddrinfo

    def _af_inet(host, port, family=0, type=0, proto=0, flags=0):
        return orig(host, port, socket.AF_INET, type, proto, flags)

    socket.getaddrinfo = _af_inet
    try:
        yield
    finally:
        socket.getaddrinfo = orig


def fetch(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"  {url}")
    request = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response, open(tmp, "wb") as out:
            while chunk := response.read(1 << 20):
                out.write(chunk)
    except urllib.error.HTTPError as exc:
        tmp.unlink(missing_ok=True)
        if exc.code != 403:
            raise
        with _ipv4_only(), urllib.request.urlopen(request, timeout=60) as response, open(tmp, "wb") as out:
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


def filtered_contigs(archive: Path, target: Path) -> None:
    """Write Cell Ranger's filtered contig table (called cells, confident
    contigs) from a gzipped ``all_contig_annotations.csv``."""
    if target.exists():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    lines = gzip.decompress(archive.read_bytes()).decode().splitlines()
    header = lines[0].split(",")
    cell, confident = header.index("is_cell"), header.index("high_confidence")
    kept = [line for line in lines[1:]
            if line.split(",")[cell].lower() == "true" and line.split(",")[confident].lower() == "true"]
    target.write_text("\n".join([lines[0], *kept]) + "\n")


def counts_table_to_matrix(table: Path, dest: Path) -> None:
    """Turn a dense gene-by-cell count table (first line: cell names, then one
    line per gene) into a Cell Ranger style ``filtered_feature_bc_matrix``.
    Cell names are ``BARCODE-library``; barcodes are written as ``BARCODE-1``."""
    if (dest / "matrix.mtx.gz").exists():
        return
    dest.mkdir(parents=True, exist_ok=True)
    genes, entries = [], []
    with gzip.open(table, "rt") as handle:
        cells = [name.split("-")[0] + "-1" for name in handle.readline().rstrip("\n").split("\t") if name]
        for row, line in enumerate(handle, start=1):
            gene, *values = line.rstrip("\n").split("\t")
            genes.append(gene)
            entries += [f"{row} {col} {value}" for col, value in enumerate(values, start=1) if value != "0"]
    with gzip.open(dest / "barcodes.tsv.gz", "wt") as out:
        out.write("\n".join(cells) + "\n")
    with gzip.open(dest / "features.tsv.gz", "wt") as out:
        out.write("".join(f"{gene}\t{gene}\tGene Expression\n" for gene in genes))
    tmp = dest / "matrix.mtx.gz.part"
    with gzip.open(tmp, "wt", compresslevel=3) as out:
        out.write("%%MatrixMarket matrix coordinate integer general\n")
        out.write(f"{len(genes)} {len(cells)} {len(entries)}\n")
        out.write("\n".join(entries) + "\n")
    tmp.rename(dest / "matrix.mtx.gz")


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

    print("Eberhardt 2021 (GSE180268)")
    listing = fetch(f"{EBERHARDT_SERIES}/filelist.txt", root / "eberhardt" / "filelist.txt").read_text()
    for gsm, name, patient, kind in EBERHARDT_FILE.findall(listing):
        if patient not in EBERHARDT_PATIENTS:
            continue
        filename = f"{gsm}_{name}_{kind}.gz"
        archive = fetch(f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{filename}", root / "eberhardt" / "raw" / filename)
        sample = root / "eberhardt" / "samples" / name
        if kind.startswith("all_contig"):
            # GEO has the unfiltered contig table; apply Cell Ranger's own
            # definition of the filtered one (called cells, confident contigs).
            filtered_contigs(archive, sample / "vdj" / "filtered_contig_annotations.csv")
        else:
            target = sample / "gex" / "filtered_feature_bc_matrix" / f"{kind}.gz"
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(archive.read_bytes())

    print("Meng 2023 (GSE254250)")
    for name in MENG_FILES:
        fetch(f"{MENG}/{name}", root / "meng" / name)

    print("Chiffelle 2024 (GSE222448)")
    fetch(CHIFFELLE_LABELS, root / "chiffelle" / "table_s3.xlsx")
    listing = fetch(f"{CHIFFELLE_SERIES}/filelist.txt", root / "chiffelle" / "filelist.txt").read_text()
    for pattern, kind in ((CHIFFELLE_COUNTS, "gex"), (CHIFFELLE_CONTIGS, "vdj")):
        for match in pattern.finditer(listing):
            gsm, patient, biopsy = match.groups()
            if int(patient) not in CHIFFELLE_PATIENTS:
                continue
            archive = fetch(f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{match.group(0)}",
                            root / "chiffelle" / "raw" / match.group(0))
            sample = root / "chiffelle" / "samples" / f"P{patient}_tumor_{biopsy}"
            if kind == "gex":
                counts_table_to_matrix(archive, sample / "gex" / "filtered_feature_bc_matrix")
            else:
                filtered_contigs(archive, sample / "vdj" / "filtered_contig_annotations.csv")

    print("Pulliam 2024 (GSE227708)")
    fetch(PULLIAM_LABELS, root / "pulliam" / "table_s4.xlsx")
    fetch(f"{PULLIAM_SERIES}/GSE227708_meta_data.csv.gz", root / "pulliam" / "meta_data.csv.gz")
    listing = fetch(f"{PULLIAM_SERIES}/filelist.txt", root / "pulliam" / "filelist.txt").read_text()
    for gsm, run, kind in PULLIAM_FILE.findall(listing):
        if run in PULLIAM_RUNS:
            fetch(f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_{run}_{kind}.gz",
                  root / "pulliam" / "runs" / run / f"{kind}.gz")
    print(f"done: {root}")


if __name__ == "__main__":
    main()
