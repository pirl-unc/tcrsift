"""Tests for the packaged sample-sheet TIL prioritization workflow."""

import runpy
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tcrsift import til_prioritize
from tcrsift.cli import create_parser, main
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
    assert til_prioritize._known_mart1_mask(df).tolist() == [
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
    assert til_prioritize._trav12_2_mask(df).tolist() == [
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
    clonotypes, clone_sample = til_prioritize._aggregate_within_patients(
        adata,
        aggregate_clonotypes,
        build_clone_sample_long,
    )

    assert len(clonotypes) == 2
    assert clonotypes["donor"].tolist() == ["patient_1", "patient_2"]
    assert clone_sample[["donor", "CDR3ab"]].drop_duplicates().shape[0] == 2


def test_parser_defaults_keep_heuristic_filters_auditable(tmp_path):
    args = create_parser().parse_args(["til-prioritize", "samples.yaml", "-o", str(tmp_path)])
    assert args.exclude_known_viral
    assert args.exclude_known_mart1
    assert not args.exclude_trav12_2
    assert args.exclude_public_quantile is None


def test_parser_accepts_example_options(tmp_path):
    args = create_parser().parse_args([
        "til-prioritize", "samples.csv", "-o", str(tmp_path),
        "--min-cells", "5", "--min-frequency", "0.01",
        "--signature-quantile", "0.95", "--min-signature-support", "2",
        "--vdjdb", "vdjdb.tsv", "--iedb", "iedb.tsv", "--cedar", "cedar.tsv",
        "--database-match", "strict_ab", "--no-exclude-known-viral",
        "--no-exclude-known-mart1", "--exclude-trav12-2",
        "--exclude-public-quantile", "0.9", "--verbose",
    ])
    assert args.sample_sheet == Path("samples.csv")
    assert args.output_dir == tmp_path
    assert (args.min_cells, args.min_frequency) == (5, 0.01)
    assert (args.signature_quantile, args.min_signature_support) == (0.95, 2)
    assert (args.vdjdb, args.iedb, args.cedar) == tuple(map(Path, ["vdjdb.tsv", "iedb.tsv", "cedar.tsv"]))
    assert args.database_match == "strict_ab"
    assert not args.exclude_known_viral and not args.exclude_known_mart1
    assert args.exclude_trav12_2 and args.exclude_public_quantile == 0.9
    assert args.verbose


def _til_cells():
    """Two patients sharing an expanded TCR; full gene universe for real scoring."""
    from tcrsift.signature_methods import NEOANTIGEN_SIGNATURES, SIGNATURES

    registry = {**SIGNATURES, **NEOANTIGEN_SIGNATURES}
    genes = sorted({g for name in til_prioritize.SIGNATURE_NAMES for g in registry[name].all_genes})
    genes += [f"background_{i}" for i in range(500)]
    rng = np.random.default_rng(4)
    X = rng.poisson(np.linspace(1, 12, len(genes)), size=(12, len(genes))).astype(float) + 1
    expanded = np.array([True] * 4 + [False] * 2)  # repeated in both patients
    for gene in ["PRF1", "GZMB", "CXCL13", "ENTPD1"]:
        X[:, genes.index(gene)] = np.tile(np.where(expanded, 100, 1), 2)
    obs = pd.DataFrame({
        "sample": ["sample_1"] * 6 + ["sample_2"] * 6,
        "patient_id": ["patient_1"] * 6 + ["patient_2"] * 6,
        "CDR3_alpha": (["CAVSDGGSQGNLIF"] * 4 + ["CAVNAGGSQGNLIF"] * 2) * 2,
        "CDR3_beta": (["CASSLGQAYEQYF"] * 4 + ["CASSLAGAYEQYF"] * 2) * 2,
        "TRA_1_umis": 3,
        "TRB_1_umis": 3,
        "CD3": 30,
        "CD4": [0] * 6 + [20] * 6,
        "CD8": [20] * 6 + [0] * 6,
    }, index=[f"cell_{i}" for i in range(12)])
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))


@pytest.mark.parametrize("entrypoint", ["cli", "example"])
def test_workflow_writes_scored_and_selected_clones(tmp_path, monkeypatch, capsys, entrypoint):
    import tcrsift

    # Loading is covered by loader tests. Exercise the actual scoring,
    # clone aggregation, publicness, annotation, filtering, and CSV writes.
    cells = _til_cells()

    def load_samples(path):
        assert path == Path("samples.yaml")
        return cells.copy()

    monkeypatch.setattr(tcrsift, "load_samples", load_samples)
    argv = ["samples.yaml", "-o", str(tmp_path), "--min-cells", "3"]
    if entrypoint == "cli":
        main(["til-prioritize", *argv])
    else:
        script = Path(__file__).resolve().parents[1] / "examples" / "multi_sample_til.py"
        monkeypatch.setattr(sys, "argv", [str(script), *argv])
        runpy.run_path(str(script), run_name="__main__")

    candidates = pd.read_csv(tmp_path / "candidate_clones.csv")
    audit = pd.read_csv(tmp_path / "all_scored_clones.csv")
    scores = pd.read_csv(tmp_path / "clone_sample_scores.csv")
    assert len(candidates) == 2
    assert candidates["donor"].tolist() == ["patient_1", "patient_2"]
    assert candidates["CDR3ab"].tolist() == ["CAVSDGGSQGNLIF_CASSLGQAYEQYF"] * 2
    assert candidates["cell_count"].tolist() == [4, 4]
    assert candidates["selected_for_review"].all()
    assert len(audit) == len(scores) == 4
    assert audit["selected_for_review"].tolist() == [True, True, False, False]
    for name in til_prioritize.SIGNATURE_NAMES:
        assert np.isfinite(scores[f"signature_{name}"]).all()
    assert "Wrote 2 candidates" in capsys.readouterr().out


@pytest.mark.parametrize("option,value", [
    ("--min-cells", "0"), ("--min-frequency", "1.1"),
    ("--signature-quantile", "0"), ("--min-signature-support", "7"),
    ("--exclude-public-quantile", "0"),
])
def test_invalid_threshold_fails_before_loading(tmp_path, monkeypatch, caplog, option, value):
    import tcrsift

    def unexpected_load(*args, **kwargs):
        pytest.fail("Invalid thresholds must fail before loading data")

    monkeypatch.setattr(tcrsift, "load_samples", unexpected_load)
    with pytest.raises(SystemExit) as error:
        main(["til-prioritize", "samples.yaml", "-o", str(tmp_path), option, value])
    assert error.value.code == 1
    assert option in caplog.text


def test_single_sample_has_clear_error(tmp_path, monkeypatch, caplog):
    import tcrsift

    cells = _til_cells()[:6].copy()
    monkeypatch.setattr(tcrsift, "load_samples", lambda path: cells)
    with pytest.raises(SystemExit) as error:
        main(["til-prioritize", "samples.yaml", "-o", str(tmp_path)])
    assert error.value.code == 1
    assert "at least two named samples" in caplog.text
