"""Tests for the packaged sample-sheet TIL prioritization workflow."""

import json
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
    assert not args.exclude_known_viral
    assert args.exclude_known_mart1
    assert not args.exclude_trav12_2
    assert args.exclude_public_quantile is None
    assert args.tcell_type == "cd8"
    assert args.context == "solid-tumor"
    assert args.max_clones == 200
    assert args.signature_quantile == 0.9
    assert args.min_signature_score == 0.0


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
    genes = sorted({g for signature in registry.values() for g in signature.all_genes})
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
        "TRA_1_reads": 30,
        "TRB_1_reads": 30,
        "CD3": 30,
        "CD4": [0] * 6 + [20] * 6,
        "CD8": [20] * 6 + [0] * 6,
    }, index=[f"cell_{i}" for i in range(12)])
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))


def _mock_samples(tmp_path, monkeypatch, cells):
    from tcrsift import loader

    sheet = tmp_path / "samples.yaml"
    sheet.write_text("samples:\n" + "".join(
        f"  - sample: {sample}\n    vdj_dir: vdj\n    gex_dir: gex\n"
        for sample in cells.obs["sample"].unique()
    ))
    def load_samples(sample_sheet, **kwargs):
        assert [s.sample for s in sample_sheet] == list(cells.obs["sample"].unique())
        assert kwargs["min_mito_pct"] == 0
        return cells.copy()
    monkeypatch.setattr(loader, "load_samples", load_samples)
    return str(sheet)


@pytest.mark.parametrize("entrypoint", ["cli", "example"])
def test_workflow_writes_scored_and_selected_clones(tmp_path, monkeypatch, capsys, entrypoint):
    # Loading is covered by loader tests. Exercise the actual scoring,
    # clone aggregation, publicness, annotation, filtering, and CSV writes.
    cells = _til_cells()

    sheet = _mock_samples(tmp_path, monkeypatch, cells)
    argv = [sheet, "-o", str(tmp_path), "--min-cells", "3", "--tcell-type", "both"]
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
    assert candidates.selection_rank.tolist() == [1, 2]
    assert len(audit) == len(scores) == 4
    assert audit["selected_for_review"].sum() == 2
    config = json.loads((tmp_path / "prioritization.json").read_text())
    for name in config["resolved_signatures"]:
        from tcrsift.prioritize import SIGNATURE_LINEAGE
        valid = scores.lineage.eq(SIGNATURE_LINEAGE[name]) if name in SIGNATURE_LINEAGE else np.ones(len(scores), bool)
        assert np.isfinite(scores.loc[valid, f"signature_{name}"]).all()
        assert scores.loc[~valid, f"signature_{name}"].isna().all()
    pass_cols = [f"signature_{name}_passes_cutoff" for name in config["resolved_signatures"]]
    support = scores.groupby(["donor", "CDR3ab"])[pass_cols].max().sum(axis=1)
    assert audit.set_index(["donor", "CDR3ab"]).signature_support_count.to_dict() == support.to_dict()
    assert "Wrote 2 candidates" in capsys.readouterr().out


@pytest.mark.parametrize("option,value", [
    ("--min-cells", "0"), ("--min-frequency", "1.1"),
    ("--signature-quantile", "1.1"), ("--min-signature-support", "999"),
    ("--exclude-public-quantile", "0"),
    ("--max-clones", "-1"), ("--min-expression", "GZMB=nan"),
    ("--min-vdj-reads", "-1"), ("--min-alpha-cdr3-length", "-1"),
    ("--min-signature-score", "nan"), ("--min-signature-score", "inf"),
])
def test_invalid_threshold_fails_before_loading(tmp_path, monkeypatch, caplog, option, value):
    from tcrsift import loader

    def unexpected_load(*args, **kwargs):
        pytest.fail("Invalid thresholds must fail before loading data")

    monkeypatch.setattr(loader, "load_samples", unexpected_load)
    with pytest.raises(SystemExit) as error:
        main(["til-prioritize", "samples.yaml", "-o", str(tmp_path), option, value])
    assert error.value.code == 1
    assert option in caplog.text


def test_single_sample_and_default_cd8_are_supported(tmp_path, monkeypatch):
    cells = _til_cells()[:6].copy()
    sheet = _mock_samples(tmp_path, monkeypatch, cells)
    main(["prioritize", sheet, "-o", str(tmp_path), "--max-clones", "1"])
    candidates = pd.read_csv(tmp_path / "candidate_clones.csv")
    assert len(candidates) == 1
    assert candidates.Tcell_type_consensus.str.contains("CD8").all()
    assert candidates.selection_rank.tolist() == [1]


@pytest.mark.parametrize("exclude_viral,expected", [(False, 1), (True, 0)])
def test_mart1_excluded_by_default_and_viral_opt_in(tmp_path, monkeypatch, exclude_viral, expected):
    from tcrsift import annotate

    cells = _til_cells()[:6].copy()
    sheet = _mock_samples(tmp_path, monkeypatch, cells)

    def annotate_clonotypes(frame, **kwargs):
        result = frame.copy()
        mart1 = result.CDR3_alpha.eq("CAVSDGGSQGNLIF")
        result["db_epitope"] = np.where(mart1, "EAAGIGILTV", "NLVPMVATV")
        result["is_viral"] = ~mart1
        return result

    monkeypatch.setattr(annotate, "annotate_clonotypes", annotate_clonotypes)
    options = ["--exclude-known-viral"] if exclude_viral else []
    main(["prioritize", sheet, "-o", str(tmp_path), "--signatures", "Cytolytic",
          "--signature-quantile", "0", "--min-signature-score", "-100", *options])
    candidates = pd.read_csv(tmp_path / "candidate_clones.csv")
    audit = pd.read_csv(tmp_path / "all_scored_clones.csv")
    assert len(candidates) == expected
    mart1 = audit[audit.known_mart1_match]
    assert len(mart1) == 1 and not mart1.selected_for_review.any()
    assert mart1.excluded_reason.tolist() == ["known_MART1"]
    if not exclude_viral:
        assert candidates.known_viral_match.all()


def test_short_cdr3_is_audited_and_excluded(tmp_path, monkeypatch):
    cells = _til_cells()[:6].copy()
    cells.obs.loc[cells.obs.CDR3_alpha.eq("CAVSDGGSQGNLIF"), "CDR3_alpha"] = "CAVF"
    sheet = _mock_samples(tmp_path, monkeypatch, cells)
    main(["prioritize", sheet, "-o", str(tmp_path), "--signatures", "Cytolytic",
          "--min-alpha-cdr3-length", "8", "--signature-quantile", "0", "--min-signature-score", "-100"])
    audit = pd.read_csv(tmp_path / "all_scored_clones.csv")
    short = audit[audit.alpha_cdr3_length == 4]
    assert len(short) == 1
    assert short.excluded_reason.tolist() == ["short_alpha_cdr3"]
    assert not short.selected_for_review.any()
    assert len(pd.read_csv(tmp_path / "candidate_clones.csv")) == 1


def test_cli_budget_is_total_even_with_multiple_patients(tmp_path, monkeypatch):
    sheet = _mock_samples(tmp_path, monkeypatch, _til_cells())
    main(["prioritize", sheet, "-o", str(tmp_path), "--tcell-type", "both", "--max-clones", "1"])
    candidates = pd.read_csv(tmp_path / "candidate_clones.csv")
    audit = pd.read_csv(tmp_path / "all_scored_clones.csv")
    assert len(candidates) == 1
    assert candidates.selection_rank.tolist() == [1]
    assert audit.eligible_for_review.sum() > len(candidates)
