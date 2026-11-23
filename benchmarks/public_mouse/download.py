#!/usr/bin/env python3
"""Download the open-access public murine datasets used by the mouse TCRsift benchmark.

    python benchmarks/public_mouse/download.py --data-dir /path/to/mouse_bench

About 800 MB total across three open-access NCBI GEO series (no controlled access):

1. Multi-tissue murine MC38 colon adenocarcinoma (GEO GSE168944):
   Matched Tumor, Peripheral Blood, Draining Lymph Node (dLN), and Spleen 10x
   scRNA-seq + paired scTCR-seq across biological replicate mice (Mouse01 &
   Mouse02; 8 paired GEX + VDJ libraries: GSM5172668-GSM5172675 and
   GSM5172684-GSM5172691).
2. Kasmani et al., Life Science Alliance 2022 (PMID:35667687; GEO GSE200608):
   NOD mouse pancreatic islet and spleen CD8+ T cells FACS-sorted with H-2Kd
   IGRP206-214 / NRP-V7 (KYNKANVFL) MHC-I tetramers vs islet tetramer-negative
   CD8+ T cells across two biological rounds (6 paired GEX + VDJ libraries:
   GSM6038943-GSM6038954). Also represented in VDJdb MusMusculus (780 entries),
   enabling direct circular vs study-held-out database evaluation.
3. Khatun et al., J Exp Med 2021 (PMID:33201171; GEO GSE158896):
   C57BL/6 mouse splenic virus-specific CD4+ T cells FACS-sorted with I-Ab
   LCMV GP66-77 (DIYKGVYQFKSV) MHC-II tetramers at Day 10 post-acute LCMV
   Armstrong infection across five mice (M1-M5; 5 paired GEX + VDJ libraries:
   GSM4815363-GSM4815372).
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import shutil
import socket
import urllib.error
import urllib.request
from pathlib import Path

GEO = "https://ftp.ncbi.nlm.nih.gov/geo"

# 1. GSE168944 — Multi-compartment MC38 murine tumor model (Mouse01 & Mouse02)
MC38_SAMPLES = [
    # (sample_id, mouse_id, tissue, gex_gsm, gex_tag, vdj_gsm, vdj_tag)
    ("Mouse01_Tumor", "Mouse01", "tumor", "GSM5172668", "DJC01", "GSM5172684", "DJT01"),
    ("Mouse02_Tumor", "Mouse02", "tumor", "GSM5172669", "DJC02", "GSM5172685", "DJT02"),
    ("Mouse01_Blood", "Mouse01", "blood", "GSM5172670", "DJC04", "GSM5172686", "DJT04"),
    ("Mouse02_Blood", "Mouse02", "blood", "GSM5172671", "DJC05", "GSM5172687", "DJT05"),
    ("Mouse01_Spleen", "Mouse01", "spleen", "GSM5172672", "DJC09", "GSM5172688", "DJT09"),
    ("Mouse02_Spleen", "Mouse02", "spleen", "GSM5172673", "DJC10", "GSM5172689", "DJT10"),
    ("Mouse01_dLN", "Mouse01", "dLN", "GSM5172674", "DJC12", "GSM5172690", "DJT12"),
    ("Mouse02_dLN", "Mouse02", "dLN", "GSM5172675", "DJC13", "GSM5172691", "DJT13"),
]

# 2. GSE200608 — Kasmani et al. 2022 NOD mouse islet & spleen IGRP tetramer+/- CD8 T cells
KASMANI_SAMPLES = [
    # (sample_id, round_id, tissue, sort_group, gex_gsm, gex_prefix, vdj_gsm, vdj_prefix)
    ("R1_Islet_TetNeg", "Round1", "islet", "tet_neg", "GSM6038943", "Round_1_GEX_Isl_tet_N", "GSM6038949", "Round_1_TCR_Isl_tet_N"),
    ("R1_Islet_TetPos", "Round1", "islet", "tet_pos", "GSM6038944", "Round_1_GEX_Isl_tet_P", "GSM6038950", "Round_1_TCR_Isl_tet_P"),
    ("R1_Spleen_TetPos", "Round1", "spleen", "tet_pos", "GSM6038945", "Round_1_GEX_Spl_tet_P", "GSM6038951", "Round_1_TCR_Spl_tet_P"),
    ("R2_Islet_TetNeg", "Round2", "islet", "tet_neg", "GSM6038946", "Round_2_GEX_Isl_tet_N", "GSM6038952", "Round_2_TCR_Isl_tet_N"),
    ("R2_Islet_TetPos", "Round2", "islet", "tet_pos", "GSM6038947", "Round_2_GEX_Isl_tet_P", "GSM6038953", "Round_2_TCR_Isl_tet_P"),
    ("R2_Spleen_TetPos", "Round2", "spleen", "tet_pos", "GSM6038948", "Round_2_GEX_Spl_tet_P", "GSM6038954", "Round_2_TCR_Spl_tet_P"),
]

# 3. GSE158896 — Khatun et al. 2021 LCMV Armstrong GP66:I-Ab tetramer+ murine CD4 T cells
KHATUN_SAMPLES = [
    ("M1_GP66_CD4", "M1", "GSM4815363", "M1_GEX", "GSM4815368", "M1_TCR"),
    ("M2_GP66_CD4", "M2", "GSM4815364", "M2_GEX", "GSM4815369", "M2_TCR"),
    ("M3_GP66_CD4", "M3", "GSM4815365", "M3_GEX", "GSM4815370", "M3_TCR"),
    ("M4_GP66_CD4", "M4", "GSM4815366", "M4_GEX", "GSM4815371", "M4_TCR"),
    ("M5_GP66_CD4", "M5", "GSM4815367", "M5_GEX", "GSM4815372", "M5_TCR"),
]


@contextlib.contextmanager
def _ipv4_only():
    """Force AF_INET resolution when a dual-stack GEO mirror returns HTTP 403 or hangs."""
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


def fetch(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
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


def _gsm_url(gsm: str, filename: str) -> str:
    return f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{filename}"


def _genes_to_features_gz(genes_gz: Path, features_gz: Path) -> None:
    """Convert Cell Ranger v2 2-column genes.tsv.gz to 3-column features.tsv.gz."""
    if features_gz.exists() and features_gz.stat().st_size > 0:
        return
    features_gz.parent.mkdir(parents=True, exist_ok=True)
    lines = gzip.decompress(genes_gz.read_bytes()).decode("utf-8", errors="replace").splitlines()
    out_lines = []
    for ln in lines:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) >= 3:
            out_lines.append(ln)
        elif len(parts) == 2:
            out_lines.append(f"{parts[0]}\t{parts[1]}\tGene Expression")
        else:
            out_lines.append(f"{parts[0]}\t{parts[0]}\tGene Expression")
    with gzip.open(features_gz, "wt", compresslevel=3) as fh:
        fh.write("\n".join(out_lines) + "\n")


def _gunzip_to(src_gz: Path, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src_gz, "rb") as src, open(dest, "wb") as dst:
        shutil.copyfileobj(src, dst)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True, type=Path)
    root = parser.parse_args().data_dir

    print("1/3 MC38 multi-tissue murine tumor cohort (GSE168944)")
    for sample_id, _mouse, _tissue, gex_gsm, gex_tag, vdj_gsm, vdj_tag in MC38_SAMPLES:
        raw_dir = root / "mc38" / "raw"
        sdir = root / "mc38" / "samples" / sample_id
        mtx_dir = sdir / "gex" / "filtered_feature_bc_matrix"
        mtx_dir.mkdir(parents=True, exist_ok=True)

        bc_raw = fetch(
            _gsm_url(gex_gsm, f"{gex_gsm}_{gex_tag}_barcodes.tsv.gz"),
            raw_dir / f"{gex_gsm}_{gex_tag}_barcodes.tsv.gz",
        )
        gn_raw = fetch(
            _gsm_url(gex_gsm, f"{gex_gsm}_{gex_tag}_genes.tsv.gz"),
            raw_dir / f"{gex_gsm}_{gex_tag}_genes.tsv.gz",
        )
        mx_raw = fetch(_gsm_url(gex_gsm, f"{gex_gsm}_{gex_tag}_matrix.mtx.gz"), raw_dir / f"{gex_gsm}_{gex_tag}_matrix.mtx.gz")
        vdj_raw = fetch(_gsm_url(vdj_gsm, f"{vdj_gsm}_{vdj_tag}_filtered_contig_annotations.csv.gz"), raw_dir / f"{vdj_gsm}_{vdj_tag}_filtered_contig_annotations.csv.gz")

        if not (mtx_dir / "barcodes.tsv.gz").exists():
            shutil.copyfile(bc_raw, mtx_dir / "barcodes.tsv.gz")
        if not (mtx_dir / "matrix.mtx.gz").exists():
            shutil.copyfile(mx_raw, mtx_dir / "matrix.mtx.gz")
        _genes_to_features_gz(gn_raw, mtx_dir / "features.tsv.gz")
        _gunzip_to(vdj_raw, sdir / "vdj" / "filtered_contig_annotations.csv")

    print("2/3 Kasmani 2022 NOD mouse islet & spleen IGRP tetramer CD8 cohort (GSE200608)")
    for sample_id, _round, _tissue, _sort, gex_gsm, gex_pref, vdj_gsm, vdj_pref in KASMANI_SAMPLES:
        raw_dir = root / "kasmani" / "raw"
        sdir = root / "kasmani" / "samples" / sample_id
        mtx_dir = sdir / "gex" / "filtered_feature_bc_matrix"
        mtx_dir.mkdir(parents=True, exist_ok=True)

        for kind in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"):
            fname = f"{gex_gsm}_{gex_pref}_{kind}"
            src = fetch(_gsm_url(gex_gsm, fname), raw_dir / fname)
            if not (mtx_dir / kind).exists():
                shutil.copyfile(src, mtx_dir / kind)
        vdj_fname = f"{vdj_gsm}_{vdj_pref}_filtered_contig_annotations.csv.gz"
        vdj_raw = fetch(_gsm_url(vdj_gsm, vdj_fname), raw_dir / vdj_fname)
        _gunzip_to(vdj_raw, sdir / "vdj" / "filtered_contig_annotations.csv")

    print("3/3 Khatun 2021 LCMV GP66:I-Ab tetramer+ murine CD4 cohort (GSE158896)")
    for sample_id, _mouse, gex_gsm, gex_pref, vdj_gsm, vdj_pref in KHATUN_SAMPLES:
        raw_dir = root / "khatun" / "raw"
        sdir = root / "khatun" / "samples" / sample_id
        mtx_dir = sdir / "gex" / "filtered_feature_bc_matrix"
        mtx_dir.mkdir(parents=True, exist_ok=True)

        bc_raw = fetch(_gsm_url(gex_gsm, f"{gex_gsm}_{gex_pref}_barcodes.tsv.gz"), raw_dir / f"{gex_gsm}_{gex_pref}_barcodes.tsv.gz")
        gn_raw = fetch(_gsm_url(gex_gsm, f"{gex_gsm}_{gex_pref}_genes.tsv.gz"), raw_dir / f"{gex_gsm}_{gex_pref}_genes.tsv.gz")
        mx_raw = fetch(_gsm_url(gex_gsm, f"{gex_gsm}_{gex_pref}_matrix.mtx.gz"), raw_dir / f"{gex_gsm}_{gex_pref}_matrix.mtx.gz")
        vdj_raw = fetch(_gsm_url(vdj_gsm, f"{vdj_gsm}_{vdj_pref}_filtered_contig_annotations.csv.gz"), raw_dir / f"{vdj_gsm}_{vdj_pref}_filtered_contig_annotations.csv.gz")

        if not (mtx_dir / "barcodes.tsv.gz").exists():
            shutil.copyfile(bc_raw, mtx_dir / "barcodes.tsv.gz")
        if not (mtx_dir / "matrix.mtx.gz").exists():
            shutil.copyfile(mx_raw, mtx_dir / "matrix.mtx.gz")
        _genes_to_features_gz(gn_raw, mtx_dir / "features.tsv.gz")
        _gunzip_to(vdj_raw, sdir / "vdj" / "filtered_contig_annotations.csv")

    print(f"done: {root}")


if __name__ == "__main__":
    main()
