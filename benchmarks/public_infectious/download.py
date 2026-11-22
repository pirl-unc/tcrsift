#!/usr/bin/env python3
"""Download open-access public infectious-disease scRNA+scTCR datasets (~750 MB).

    python benchmarks/public_infectious/download.py --data-dir /path/to/cache

Datasets downloaded (no controlled access required):

1. Fischer et al., Nature Communications 2021 (COVID-19 PBMC reverse phenotyping,
   GEO GSE171037): 4 paired unstimulated and SARS-CoV-2 Spike peptide-stimulated
   libraries (GSM5216975-GSM5216978) with GEX .h5, filtered contig .csv.gz, and
   all_contig_annotations.json.gz.
2. 10x Genomics Multiplexed 50-Dextramer Peripheral Blood CD8+ T cells (Donors 1
   and 2): filtered_feature_bc_matrix.h5 (GEX + 50 pMHC dextramers + CITE-seq),
   all_contig_annotations.csv, and binarized_matrix.csv.
3. Minervina et al., Nature Immunology 2022 (SARS-CoV-2 infection, vaccination,
   and breakthrough CD8+ T cells): Zenodo 6232103 aggregated 10x GEX+VDJ archive
   plus curated antigen-specific cell table (pogorely/COVID_vax_CD8).
"""

from __future__ import annotations

import argparse
import urllib.request
import zipfile
from pathlib import Path

GEO = "https://ftp.ncbi.nlm.nih.gov/geo"
FISCHER_SUPPL = f"{GEO}/series/GSE171nnn/GSE171037/suppl"
FISCHER_SAMPLES = {
    "GT6_unstim": ("GSM5216975", "Sample_20063a009_01"),
    "GT6_covid_stim": ("GSM5216976", "Sample_20063a010_01"),
    "GT_unstim": ("GSM5216977", "Sample_20063a011_01"),
    "GT_covid_stim": ("GSM5216978", "Sample_20063a012_01"),
}

TENX_BASE = "https://cf.10xgenomics.com/samples/cell-vdj/3.0.2"
TENX_DONORS = {
    "donor1": "vdj_v1_hs_aggregated_donor1",
    "donor2": "vdj_v1_hs_aggregated_donor2",
}

MINERVINA_ZIP = "https://zenodo.org/api/records/6232103/files/libs_aggregate_rev.zip/content"
MINERVINA_TSV = "https://raw.githubusercontent.com/pogorely/COVID_vax_CD8/main/cd8_only_dextr_rev_clean.tsv"


def fetch(url: str, dest: Path, fallback_url: str | None = None) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for candidate in (url, fallback_url):
        if not candidate:
            continue
        try:
            print(f"  {candidate}")
            req = urllib.request.Request(candidate, headers={"User-Agent": "tcrsift-benchmark/1.0"})
            with urllib.request.urlopen(req, timeout=60) as resp, open(tmp, "wb") as out:
                while chunk := resp.read(1 << 20):
                    out.write(chunk)
            tmp.rename(dest)
            return dest
        except Exception as exc:
            tmp.unlink(missing_ok=True)
            if candidate == fallback_url or fallback_url is None:
                raise RuntimeError(f"Failed to download {candidate}: {exc}") from exc
    return dest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", required=True, type=Path)
    root = parser.parse_args().data_dir

    print("1. Fischer et al. 2021 (GSE171037 — COVID-19 PBMC unstim vs Spike stim)")
    for label, (gsm, sample_tag) in FISCHER_SAMPLES.items():
        sdir = root / "fischer_2021_covid" / label
        for suffix in (
            "filtered_feature_bc_matrix.h5",
            "filtered_contig_annotations.csv.gz",
            "all_contig_annotations.json.gz",
        ):
            fname = f"{gsm}_{sample_tag}_{suffix}"
            sample_url = f"{GEO}/samples/{gsm[:7]}nnn/{gsm}/suppl/{fname}"
            series_url = f"{FISCHER_SUPPL}/{fname}"
            fetch(sample_url, sdir / fname, fallback_url=series_url)

    print("2. 10x Genomics Multiplexed 50-Dextramer Donors 1 & 2")
    for donor, prefix in TENX_DONORS.items():
        ddir = root / "tenx_dextramer" / donor
        fetch(
            f"{TENX_BASE}/{prefix}/{prefix}_sc5p_v2_hs_B_T_all_feature_bc_matrix_filtered.h5",
            ddir / f"{prefix}_filtered_feature_bc_matrix.h5",
            fallback_url=f"{TENX_BASE}/{prefix}/{prefix}_filtered_feature_bc_matrix.h5",
        )
        for suffix in ("all_contig_annotations.csv", "binarized_matrix.csv"):
            fname = f"{prefix}_{suffix}"
            fetch(f"{TENX_BASE}/{prefix}/{fname}", ddir / fname)

    print("3. Minervina et al. 2022 (Zenodo 6232103 — SARS-CoV-2 CD8 dextramers)")
    mdir = root / "minervina_2022_covid_flu"
    fetch(MINERVINA_TSV, mdir / "cd8_only_dextr_rev_clean.tsv")
    zip_path = fetch(MINERVINA_ZIP, mdir / "libs_aggregate_rev.zip")
    target_h5 = mdir / "extracted" / "libs_aggregate_rev" / "filtered_feature_bc_matrix.h5"
    if not target_h5.exists():
        print("  Extracting libs_aggregate_rev.zip ...")
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(mdir / "extracted")

    print(f"Done: {root}")


if __name__ == "__main__":
    main()
