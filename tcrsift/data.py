# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import annotations

from os import PathLike
from os.path import join

import pandas as pd
import scanpy as sc


def annotate_combined_df(
    df,
    tcell_min_read_ratio_cd4_vs_cd8: int = 3,
    min_cd3_reads: int = 0,
    min_percent_mt: float = 2,
    max_percent_mt: float = 8,
    min_genes: int = 250,
    max_genes: int = 15_000,
    min_reads: int = 500,
    max_reads: int = 100_000,
):
    """
    Load a TSV file with multiple samples of combined gene expression and VDJ sequence data.

    multisample_tsv_path
        Expected columns:
        - "barcode" or "Barcode": nucleotide barcode for each cell
        - "sample" or "Sample": antigen specific expansion pool name
        - "cdr3_aa1" or "CDR3_alpha": amino acid sequence of TCR alpha chain
        - "cdr3_aa2" or "CDR3_beta": amino acid sequence of TCR beta chain
        - "nCount_RNA": number of reads per cell
        - "nFeature_RNA": number of genes detected per cell

        --- Gene expression columns ---
        - "ENSG00000167286.10" (or CD3D): read count of CD3 delta chain
        - "ENSG00000198851.10" (or CD3E): read count of CD3 epsilon chain
        - "ENSG00000160654.11" (or "CD3G"): read count of CD3 gamma chain
        - "ENSG00000010610.10" (or "CD4"):  read count of CD4
        - "ENSG00000153563.16" (or "CD8A"): read count of CD8 alpha chain
        - "ENSG00000172116.23" (or "CD8B"): read count of CD8 beta chain

    tcell_annotation_ratio_cutoff
        Ratio of CD4 to CD8 RNA reads to consider a cell confidently CD4+ or CD8+.


    """
    df = df.rename(
        columns={
            "ENSG00000167286.10": "CD3D",
            "ENSG00000198851.10": "CD3E",
            "ENSG00000160654.11": "CD3G",
            "ENSG00000010610.10": "CD4",
            "ENSG00000153563.16": "CD8A",
            "ENSG00000172116.23": "CD8B",
            "barcode": "Barcode",
            "Row.names": "Barcode",
            "sample": "Sample",
            "cdr3_aa1": "CDR3_alpha",
            "cdr3_aa2": "CDR3_beta",
        }
    )

    if "Peptide_Number" not in df.columns:
        df["Peptide_Number"] = [s.split("_")[-1] for s in df["Sample"]]
    df["filter:has_alpha"] = ~df["CDR3_alpha"].isnull()
    df["filter:has_beta"] = ~df["CDR3_beta"].isnull()

    df["CDR3_alpha"] = df["CDR3_alpha"].fillna("")
    df["CDR3_beta"] = df["CDR3_beta"].fillna("")
    df["CTaa"] = df["CDR3_alpha"] + "_" + df["CDR3_beta"]
    df["CDR3a/b"] = df["CDR3_alpha"] + "/" + df["CDR3_beta"]
    df["CTaa_pairs"] = [
        [aa + "_" + bb for aa in a.split(";") for bb in b.split(";")]
        for a, b in zip(df["CDR3_alpha"], df["CDR3_beta"])
    ]
    df["CDR3a/b_pairs"] = [
        [aa + "_" + bb for aa in a.split(";") for bb in b.split(";")]
        for a, b in zip(df["CDR3_alpha"], df["CDR3_beta"])
    ]
    df["Cell_ID"] = df.Peptide_Number + "-" + df.Barcode
    df["CD3"] = df["CD3D"] + df["CD3E"] + df["CD3G"]
    df["CD8"] = df["CD8A"] + df["CD8B"]
    df["Both_CD4_and_CD8"] = (df.CD4 > 1) & (df.CD8 > 1)

    df["Seq"] = df.CDR3_alpha.fillna("") + "_" + df.CDR3_beta.fillna("")
    df["CDR3_alpha_missing"] = (
        df.CDR3_alpha.isnull()
        | (df.CDR3_alpha.str.len() == 0)
        | (df.CDR3_alpha == "NA")
    )
    df["CDR3_beta_missing"] = (
        df.CDR3_beta.isnull() | (df.CDR3_beta.str.len() == 0) | (df.CDR3_beta == "NA")
    )
    df["TCR_complete"] = ~df.CDR3_alpha_missing & ~df.CDR3_beta_missing

    tcell_type = [
        (
            "Confident CD8+"
            if ((1 + cd8) / (1 + cd4)) > tcell_min_read_ratio_cd4_vs_cd8
            else (
                "Confident CD4+"
                if ((1 + cd4) / (1 + cd8)) > tcell_min_read_ratio_cd4_vs_cd8
                else (
                    "Likely CD8+"
                    if (cd8 > 0 and cd4 == 0)
                    else "Likely CD4+" if (cd4 > 0 and cd8 == 0) else "Unknown"
                )
            )
        )
        for (cd8, cd4) in zip(df.CD8, df.CD4)
    ]
    df["Tcell_type"] = pd.Categorical(tcell_type, categories=sorted(set(tcell_type)))
    df["confident_and_complete"] = (
        df["Tcell_type"].str.startswith("Confident")
    ) & df.TCR_complete
    df["filter:cd3_reads"] = df.CD3 >= min_cd3_reads
    df["filter:percent.mt"] = (df["percent.mt"] > min_percent_mt) & (
        df["percent.mt"] < max_percent_mt
    )
    df["filter:num_reads"] = (df["nCount_RNA"] > min_reads) & (
        df["nCount_RNA"] < max_reads
    )

    df["filter:num_genes"] = (df["nFeature_RNA"] > min_genes) & (
        df["nFeature_RNA"] < max_genes
    )

    filter_cols = [col for col in df.columns if col.startswith("filter:")]
    df["filter:all"] = df[filter_cols].all(axis=1)
    df["filtered_confident_and_complete"] = (
        df["confident_and_complete"] & df["filter:all"]
    )
    return df


def load_combined_tsv(multisample_tsv_path: str | PathLike, **kwargs):
    """
    Load a TSV file with multiple samples of combined gene expression and VDJ sequence data.

    multisample_tsv_path
        Expected columns:
        - "barcode" or "Barcode": nucleotide barcode for each cell
        - "sample" or "Sample": antigen specific expansion pool name
        - "cdr3_aa1" or "CDR3_alpha": amino acid sequence of TCR alpha chain
        - "cdr3_aa2" or "CDR3_beta": amino acid sequence of TCR beta chain
        - "nCount_RNA": number of reads per cell
        - "nFeature_RNA": number of genes detected per cell

        --- Gene expression columns ---
        - "ENSG00000167286.10" (or CD3D): read count of CD3 delta chain
        - "ENSG00000198851.10" (or CD3E): read count of CD3 epsilon chain
        - "ENSG00000160654.11" (or "CD3G"): read count of CD3 gamma chain
        - "ENSG00000010610.10" (or "CD4"):  read count of CD4
        - "ENSG00000153563.16" (or "CD8A"): read count of CD8 alpha chain
        - "ENSG00000172116.23" (or "CD8B"): read count of CD8 beta chain

    tcell_annotation_ratio_cutoff
        Ratio of CD4 to CD8 RNA reads to consider a cell confidently CD4+ or CD8+.


    """
    return annotate_combined_df(pd.read_csv(multisample_tsv_path, sep="\t"), **kwargs)


def load_sample_cellranger_count_and_vdj_outputs(
    gene_expression__dir: str | PathLike,
    vdj_dir: str | PathLike,
    vdj_clonotypes_filename="clonotypes.csv",
    vdj_annotations_filename="all_contig_annotations.csv",
):
    gene_expression_matrix_dir = join(
        gene_expression__dir, "filtered_feature_bc_matrix"
    )
    vdj_clonotypes_csv_path = join(vdj_dir, vdj_clonotypes_filename)
    vdj_annotations_csv_path = join(vdj_dir, vdj_annotations_filename)
    gene_expression_data = sc.read_10x_mtx(
        gene_expression_matrix_dir, var_names="gene_ids"
    )
    df_clonotypes = pd.read_csv(vdj_clonotypes_csv_path)
    df_annotations = pd.read_csv(vdj_annotations_csv_path)
    print("Loaded data")
    return (gene_expression_data, df_clonotypes, df_annotations)


def run_from_sample_sheet_for_cellranger_count_and_vdj_outputs(
    sample_sheet_csv_path: str | PathLike,
):
    """
    Load CellRanger 5'GEX+3'VDJ data for each sample based on paths in a sample sheet CSV file.
    Assumes you have run 'cellranger vdj' and 'cellranger count' on the samples.

    Parameters
    ----------
    sample_sheet_csv_path
        CSV file with columns 'sample', 'gene_expression_dir', 'vdj_dir'.
    """
    print("Running from sample sheet")
    sample_sheet = pd.read_csv(sample_sheet_csv_path)
    for _, row in sample_sheet.iterrows():
        print(row)
        (gene_expression_data, df_clonotypes, df_annotations) = (
            load_sample_cellranger_count_and_vdj_outputs(
                row["gene_expression_matrix_dir"],
                row["vdj_clonotypes_csv_path"],
                row["vdj_annotations_csv_path"],
            )
        )
