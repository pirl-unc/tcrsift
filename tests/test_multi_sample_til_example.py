"""Focused tests for the auditable multi-sample TIL example."""

import anndata as ad
import numpy as np
import pandas as pd

from examples import multi_sample_til
from tcrsift.clonotype import aggregate_clonotypes, build_clone_sample_long


def test_mart1_match_covers_native_and_altered_epitopes():
    df = pd.DataFrame(
        {
            "db_epitope": [
                "EAAGIGILTV",
                "AAGIGILTV",
                "ELAGIGILTV",
                "NLVPMVATV",
            ]
        }
    )
    assert multi_sample_til._known_mart1_mask(df).tolist() == [
        True,
        True,
        True,
        False,
    ]


def test_trav12_2_flag_is_allele_insensitive_and_exact():
    df = pd.DataFrame(
        {
            "alpha_v_gene": [
                "TRAV12-2*01",
                "trav12-2",
                "TRAV12-1*01",
                None,
            ]
        }
    )
    assert multi_sample_til._trav12_2_mask(df).tolist() == [
        True,
        True,
        False,
        False,
    ]


def test_same_clone_is_not_collapsed_across_patients():
    obs = pd.DataFrame(
        {
            "sample": ["sample_1", "sample_2"],
            "patient_id": ["patient_1", "patient_2"],
            "CDR3_alpha": ["CAVR", "CAVR"],
            "CDR3_beta": ["CASSF", "CASSF"],
            "TRA_1_umis": [3, 3],
            "TRB_1_umis": [3, 3],
        },
        index=["cell_1", "cell_2"],
    )
    adata = ad.AnnData(X=np.zeros((2, 0)), obs=obs)
    clonotypes, clone_sample = multi_sample_til._aggregate_within_patients(
        adata,
        aggregate_clonotypes,
        build_clone_sample_long,
    )

    assert len(clonotypes) == 2
    assert clonotypes["donor"].tolist() == ["patient_1", "patient_2"]
    assert clone_sample[["donor", "CDR3ab"]].drop_duplicates().shape[0] == 2


def test_parser_defaults_keep_heuristic_filters_auditable(tmp_path):
    args = multi_sample_til.create_parser().parse_args(["samples.yaml", "-o", str(tmp_path)])
    assert args.exclude_known_viral
    assert args.exclude_known_mart1
    assert not args.exclude_trav12_2
    assert args.exclude_public_quantile is None
