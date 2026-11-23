#!/usr/bin/env python3
"""Download the open-access human datasets for the paired blood-tumor, ACT, and CD4+ TIL benchmark.

    python benchmarks/public_blood_cd4/download.py --data-dir /path/to/blood_cd4_bench

About 450 MB standalone (~250 MB when reusing --til-cache-dir). No controlled access required:

1. Pulliam et al., Cell Reports Medicine 2024 (PMID:38340723; GEO GSE227708):
   Hashed peripheral blood and matched Merkel cell carcinoma (MCPyV-driven) tumor
   scRNA-seq + scTCR-seq runs (Z1210, Z1249, Z1308, Z1493, Z1493B) with pMHC
   multimer-verified MCPyV antitumor vs CMV/EBV/Influenza control viral TCRs in
   Table S4. Evaluates `--context blood-tumor` (`NeoTCR_PBL`) and blood-to-tumor
   clonal overlap.
2. Chiffelle et al., Immunity 2024 / Barras et al., Science Immunology 2024
   (PMID:39321802; GEO GSE222448):
   Paired baseline pre-ACT (`T0`) and Day-30 post-ACT (`T30`) melanoma tumor
   biopsies across Patients 1, 5, 7, 8, 10, and 13, linked to autologous
   tumor-reactive (`TR`) vs non-reactive (`NTR`) TCRs and ACT infusion-product
   expansion (`Clonally expanded in ACTP`) in Supplementary Table S3.
3. Veatch et al., Cancer Cell 2022 (PMID:35413271; GEO GSE198265):
   Human melanoma scRNA/AbSeq + TCR profiling across Patients X197, X198, X205,
   and X422 with experimentally verified autologous neoantigen-specific (`CD4_neo`,
   e.g. BRAF V600E, ERAP1, SNRPB, HTATSF1, ITGB5) and melanoma self-antigen
   (`CD4_self`) CD4+ TIL clonotypes.
4. Oh et al., Cell 2020 (PMID:32497499; GEO GSE149652):
   Whole-transcriptome 10x scRNA-seq + paired α/β TCR sequencing of human bladder
   cancer CD4+ TILs and matched adjacent non-malignant bladder CD4+ T cells across
   7 patients (Anti-PD-L1 A-D, Chemo, Untreated A-B; 13 libraries), with annotated
   cytotoxic (`CD4GZMB`, `CD4GZMK`), `CD4CXCL13`, and `FOXP3+/IL2RAhi` Treg states.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import re
import shutil
import socket
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

GEO = "https://ftp.ncbi.nlm.nih.gov/geo"

PULLIAM_SERIES = f"{GEO}/series/GSE227nnn/GSE227708/suppl"
PULLIAM_RUNS = ("Z1210", "Z1249", "Z1308", "Z1493", "Z1493B")
PULLIAM_FILE = re.compile(
    r"(GSM\d+)_(\w+?)_(barcodes\.tsv|features\.tsv|matrix\.mtx|all_contig_annotations\.csv)\.gz"
)
PULLIAM_LABELS = "https://ars.els-cdn.com/content/image/1-s2.0-S2666379124000351-mmc5.xlsx"

CHIFFELLE_SERIES = f"{GEO}/series/GSE222nnn/GSE222448/suppl"
CHIFFELLE_ACT_PATIENTS = (1, 5, 7, 8, 10, 13)
CHIFFELLE_COUNTS = re.compile(r"(GSM\d+)_patient(\d+)_(T(?:0\d?|30))_CD45-counts\.tsv\.gz")
CHIFFELLE_CONTIGS = re.compile(r"(GSM\d+)_all_contig_annotations_Patient(\d+)_(T(?:0\d?|30))\.csv\.gz")
CHIFFELLE_LABELS = "https://ars.els-cdn.com/content/image/1-s2.0-S1074761324004138-mmc3.xlsx"

VEATCH_FILES = [
    ("X197", "gex", "GSM5942454", "GSM5942454_X197_TCR_DBEC_MolsPerCell.csv.gz"),
    ("X197", "vdj", "GSM5942469", "GSM5942469_X197_TCR_VDJ_annotated-update.csv.gz"),
    ("X198", "gex", "GSM5942455", "GSM5942455_Combined_X198_TCR_DBEC_MolsPerCell.csv.gz"),
    ("X198", "vdj", "GSM5942470", "GSM5942470_X198_annotated_VDJ_update.csv.gz"),
    ("X198", "tag", "GSM5942470", "GSM5942470_X198_TCR_Sample_Tag_Calls.csv.gz"),
    ("X205", "gex", "GSM5942456", "GSM5942456_X205_DBEC_MolsPerCell.csv.gz"),
    ("X205", "vdj", "GSM5942456", "GSM5942456_X205_annotated_VDJ_update.csv.gz"),
    ("X422", "gex", "GSM5942457", "GSM5942457_Combined_X422_TCR_DBEC_MolsPerCell.csv.gz"),
    ("X422", "vdj", "GSM5942471", "GSM5942471_X422_TCR_VDJ_annotated-update.csv.gz"),
    ("X422", "tag", "GSM5942457", "GSM5942457_X422_TCR_Sample_Tag_Calls.csv.gz"),
]

OH_SAMPLES = [
    ("AntiPDL1_A_tumor", "AntiPDL1_A", "tumor", "GSM4506643", "Anti-PD-L1_A_tumor_CD4"),
    ("AntiPDL1_B_tumor", "AntiPDL1_B", "tumor", "GSM4506645", "Anti-PD-L1_B_tumor_CD4"),
    ("AntiPDL1_B_normal", "AntiPDL1_B", "normal", "GSM4506646", "Anti-PD-L1_B_normal_CD4"),
    ("AntiPDL1_C_tumor", "AntiPDL1_C", "tumor", "GSM4506648", "Anti-PD-L1_C_tumor_CD4"),
    ("AntiPDL1_C_normal", "AntiPDL1_C", "normal", "GSM4506649", "Anti-PD-L1_C_normal_CD4"),
    ("AntiPDL1_D_tumor", "AntiPDL1_D", "tumor", "GSM4506651", "Anti-PD-L1_D_tumor_CD4"),
    ("AntiPDL1_D_normal", "AntiPDL1_D", "normal", "GSM4506652", "Anti-PD-L1_D_normal_CD4"),
    ("Chemo_tumor", "Chemo", "tumor", "GSM4506654", "Chemo_tumor_CD4"),
    ("Chemo_normal", "Chemo", "normal", "GSM4506655", "Chemo_normal_CD4"),
    ("Untreated_A_tumor", "Untreated_A", "tumor", "GSM4506658", "Untreated_A_tumor_CD4"),
    ("Untreated_A_normal", "Untreated_A", "normal", "GSM4506659", "Untreated_A_normal_CD4"),
    ("Untreated_B_tumor", "Untreated_B", "tumor", "GSM4506661", "Untreated_B_tumor_CD4"),
    ("Untreated_B_normal", "Untreated_B", "normal", "GSM4506662", "Untreated_B_normal_CD4"),
]


@contextlib.contextmanager
def _ipv4_only():
    orig = socket.getaddrinfo

    def _af_inet(host, port, family=0, type=0, proto=0, flags=0):
        return orig(host, port, socket.AF_INET, type, proto, flags)

    socket.getaddrinfo = _af_inet
    try:
        yield
    finally:
        socket.getaddrinfo = orig


def _verify_gz(path: Path) -> None:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass


def fetch(url: str, dest: Path, fallback_existing: Path | None = None) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    if fallback_existing is not None and fallback_existing.exists() and fallback_existing.stat().st_size > 0:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(fallback_existing, dest)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"  {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "tcrsift-benchmark/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp, open(tmp, "wb") as out:
            while chunk := resp.read(1 << 20):
                out.write(chunk)
        _verify_gz(tmp)
    except (urllib.error.URLError, TimeoutError, OSError, EOFError) as exc:
        tmp.unlink(missing_ok=True)
        if isinstance(exc, urllib.error.HTTPError) and exc.code not in (403, 502, 503, 504):
            raise
        with _ipv4_only(), urllib.request.urlopen(req, timeout=60) as resp, open(tmp, "wb") as out:
            while chunk := resp.read(1 << 20):
                out.write(chunk)
        _verify_gz(tmp)
    tmp.rename(dest)
    return dest


def filtered_contigs(archive: Path, target: Path) -> None:
    if target.exists() and target.stat().st_size > 0:
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    lines = gzip.decompress(archive.read_bytes()).decode("utf-8", errors="replace").splitlines()
    header = lines[0].split(",")
    cell, confident = header.index("is_cell"), header.index("high_confidence")
    kept = [
        line
        for line in lines[1:]
        if line.split(",")[cell].lower() == "true" and line.split(",")[confident].lower() == "true"
    ]
    target.write_text("\n".join([lines[0], *kept]) + "\n")


def counts_table_to_matrix(table: Path, dest: Path, sep: str = "\t") -> None:
    """Convert dense gene-by-cell count table to 10x Matrix Market format."""
    if (dest / "matrix.mtx.gz").exists():
        return
    dest.mkdir(parents=True, exist_ok=True)
    genes, entries = [], []
    with gzip.open(table, "rt", encoding="utf-8", errors="replace") as handle:
        raw_cells = [c.strip() for c in handle.readline().rstrip("\n").split(sep)]
        # First header field in CSV count tables may be 'index' or empty.
        if sep == "," and raw_cells and raw_cells[0].lower() in ("", "index", "gene"):
            raw_cells = raw_cells[1:]
        else:
            raw_cells = [c for c in raw_cells if c]
        cells = [name.split("-")[0] + "-1" for name in raw_cells]
        for row, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            gene, *values = line.rstrip("\n").split(sep)
            genes.append(gene)
            for col, value in enumerate(values, start=1):
                if value not in ("0", "0.0", ""):
                    iv = int(round(float(value)))
                    if iv > 0:
                        entries.append(f"{row} {col} {iv}")
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


def oh_paired_tcr_to_contigs(tcr_gz: Path, target: Path) -> None:
    """Convert Oh et al. 2020 paired TCR matrix CSV into Cell Ranger filtered_contig_annotations.csv."""
    if target.exists() and target.stat().st_size > 0:
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(tcr_gz, compression="gzip")
    fieldnames = [
        "barcode", "is_cell", "contig_id", "high_confidence", "length", "chain",
        "v_gene", "d_gene", "j_gene", "c_gene", "full_length", "productive",
        "cdr3", "cdr3_nt", "reads", "umis", "raw_clonotype_id",
    ]
    rows = []
    for idx, r in df.iterrows():
        bc_raw = str(r.get("cell.barcode") or "").strip()
        if not bc_raw or bc_raw.lower() == "nan":
            continue
        bc = bc_raw.split("-")[0] + "-1"
        a_cdr3 = str(r.get("alphaCDR3") or "").strip()
        b_cdr3 = str(r.get("betaCDR3") or "").strip()
        cid = f"clonotype_{idx}"
        if a_cdr3 and a_cdr3.lower() != "nan":
            rows.append({
                "barcode": bc,
                "is_cell": "True",
                "contig_id": f"{bc}_contig_1",
                "high_confidence": "True",
                "length": 500,
                "chain": "TRA",
                "v_gene": str(r.get("alphaVGene") or "None"),
                "d_gene": "None",
                "j_gene": str(r.get("alphaJGene") or "None"),
                "c_gene": "TRAC",
                "full_length": "True",
                "productive": "True",
                "cdr3": a_cdr3,
                "cdr3_nt": str(r.get("alphaClonalSeq") or ""),
                "reads": int(r.get("cloneCount") or 100),
                "umis": 10,
                "raw_clonotype_id": cid,
            })
        if b_cdr3 and b_cdr3.lower() != "nan":
            b_j = str(r.get("betaJGene") or "None")
            b_c = "TRBC2" if b_j.upper().startswith("TRBJ2") else "TRBC1"
            rows.append({
                "barcode": bc,
                "is_cell": "True",
                "contig_id": f"{bc}_contig_2",
                "high_confidence": "True",
                "length": 500,
                "chain": "TRB",
                "v_gene": str(r.get("betaVGene") or "None"),
                "d_gene": str(r.get("betaDGene") or "None"),
                "j_gene": b_j,
                "c_gene": b_c,
                "full_length": "True",
                "productive": "True",
                "cdr3": b_cdr3,
                "cdr3_nt": str(r.get("betaClonalSeq") or ""),
                "reads": int(r.get("cloneCount") or 100),
                "umis": 10,
                "raw_clonotype_id": cid,
            })
    with open(target, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument(
        "--til-cache-dir",
        type=Path,
        default=None,
        help="Optional existing public_til cache directory to reuse shared Pulliam/Chiffelle files",
    )
    args = parser.parse_args()
    root = args.data_dir
    til_cache = args.til_cache_dir

    print("1/4 Pulliam 2024 paired Blood + Tumor Merkel cell carcinoma (GSE227708)")
    fetch(
        PULLIAM_LABELS,
        root / "pulliam" / "table_s4.xlsx",
        (til_cache / "pulliam" / "table_s4.xlsx") if til_cache else None,
    )
    fetch(
        f"{PULLIAM_SERIES}/GSE227708_meta_data.csv.gz",
        root / "pulliam" / "meta_data.csv.gz",
        (til_cache / "pulliam" / "meta_data.csv.gz") if til_cache else None,
    )
    listing = fetch(
        f"{PULLIAM_SERIES}/filelist.txt",
        root / "pulliam" / "filelist.txt",
        (til_cache / "pulliam" / "filelist.txt") if til_cache else None,
    ).read_text()
    for gsm, run, kind in PULLIAM_FILE.findall(listing):
        if run in PULLIAM_RUNS:
            fb = (til_cache / "pulliam" / "runs" / run / f"{kind}.gz") if til_cache else None
            fetch(
                f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_{run}_{kind}.gz",
                root / "pulliam" / "runs" / run / f"{kind}.gz",
                fb,
            )

    print("2/4 Chiffelle/Barras 2024 longitudinal pre-ACT T0 & post-ACT T30 melanoma (GSE222448)")
    fetch(
        CHIFFELLE_LABELS,
        root / "chiffelle" / "table_s3.xlsx",
        (til_cache / "chiffelle" / "table_s3.xlsx") if til_cache else None,
    )
    c_listing = fetch(
        f"{CHIFFELLE_SERIES}/filelist.txt",
        root / "chiffelle" / "filelist.txt",
        (til_cache / "chiffelle" / "filelist.txt") if til_cache else None,
    ).read_text()
    for pattern, kind in ((CHIFFELLE_COUNTS, "gex"), (CHIFFELLE_CONTIGS, "vdj")):
        for match in pattern.finditer(c_listing):
            gsm, patient, biopsy = match.groups()
            if int(patient) not in CHIFFELLE_ACT_PATIENTS:
                continue
            sample_dir = root / "chiffelle" / "samples" / f"P{patient}_tumor_{biopsy}"
            existing_sdir = (
                (til_cache / "chiffelle" / "samples" / f"P{patient}_tumor_{biopsy}")
                if til_cache
                else None
            )
            if (
                existing_sdir is not None
                and (existing_sdir / "gex" / "filtered_feature_bc_matrix" / "matrix.mtx.gz").exists()
                and (existing_sdir / "vdj" / "filtered_contig_annotations.csv").exists()
                and not sample_dir.exists()
            ):
                shutil.copytree(existing_sdir, sample_dir)
                continue
            fb_raw = (til_cache / "chiffelle" / "raw" / match.group(0)) if til_cache else None
            archive = fetch(
                f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{match.group(0)}",
                root / "chiffelle" / "raw" / match.group(0),
                fb_raw,
            )
            if kind == "gex":
                counts_table_to_matrix(archive, sample_dir / "gex" / "filtered_feature_bc_matrix", sep="\t")
            else:
                filtered_contigs(archive, sample_dir / "vdj" / "filtered_contig_annotations.csv")

    print("3/4 Veatch 2022 human melanoma neoantigen-specific CD4+ TILs (GSE198265)")
    for patient, _kind, gsm, fname in VEATCH_FILES:
        fetch(
            f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{fname}",
            root / "veatch" / patient / fname,
        )

    print("4/4 Oh 2020 human bladder cancer paired α/β CD4+ TILs + normal (GSE149652)")
    for sample_id, _patient, _tissue, gsm, prefix in OH_SAMPLES:
        raw_dir = root / "oh" / "raw"
        sdir = root / "oh" / "samples" / sample_id
        tcr_gz = fetch(
            f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_{prefix}_TCR_paired_matrice.csv.gz",
            raw_dir / f"{gsm}_{prefix}_TCR_paired_matrice.csv.gz",
        )
        ci_gz = fetch(
            f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_{prefix}_droplet_cellinfo_matrice.csv.gz",
            raw_dir / f"{gsm}_{prefix}_droplet_cellinfo_matrice.csv.gz",
        )
        cnt_gz = fetch(
            f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_{prefix}_droplet_count_matrice.csv.gz",
            raw_dir / f"{gsm}_{prefix}_droplet_count_matrice.csv.gz",
        )
        sdir.mkdir(parents=True, exist_ok=True)
        if not (sdir / "cellinfo.csv.gz").exists():
            shutil.copyfile(ci_gz, sdir / "cellinfo.csv.gz")
        counts_table_to_matrix(cnt_gz, sdir / "gex" / "filtered_feature_bc_matrix", sep=",")
        oh_paired_tcr_to_contigs(tcr_gz, sdir / "vdj" / "filtered_contig_annotations.csv")

    print(f"done: {root}")


if __name__ == "__main__":
    main()
