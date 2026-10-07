"""Behavioral tests for context presets, shared cell gates, and balanced selection."""

import numpy as np
import pandas as pd
import pytest

from tcrsift.cli import create_parser
from tcrsift.clonotype import aggregate_clonotypes, build_clone_sample_long
from tcrsift.phenotype import phenotype_cells
from tcrsift.prioritize import (
    CONTEXT_SIGNATURES,
    filter_cells,
    gene_filter_mask,
    resolve_signatures,
    select_round_robin,
    signature_pass_mask,
)
from tests.test_til_prioritize import _til_cells


@pytest.mark.parametrize("context", CONTEXT_SIGNATURES)
def test_presets_respect_lineage(context):
    cd8 = resolve_signatures(context, "cd8")
    cd4 = resolve_signatures(context, "cd4")
    assert "NeoTCR4" not in cd8
    assert not {"NeoTCR8", "NeoTCR_PBL", "MANAscore"}.intersection(cd4)
    assert {"AcuteActivation", "Cytolytic", "Differentiated"}.issubset(cd8)


def test_contexts_and_manual_overrides():
    blood = resolve_signatures("blood", "cd8")
    heme = resolve_signatures("heme", "cd8")
    assert blood == heme
    assert "CirculatingMemory" in blood
    assert not {"TumorReactive", "NeoTCR8", "NeoTCR_PBL", "MANAscore"}.intersection(blood)
    assert "NeoTCR_PBL" in resolve_signatures("blood-tumor", "cd8")
    mpe = resolve_signatures("mpe", "both")
    assert "TumorReactive" not in mpe and {"NeoTCR8", "NeoTCR4"}.issubset(mpe)
    assert resolve_signatures("solid-tumor", "cd8", ["cytolytic", "AIM", "Cytolytic"], ["aim"]) == ["Cytolytic"]
    with pytest.raises(ValueError, match="incompatible"):
        resolve_signatures("generic", "cd8", ["NeoTCR4"])
    with pytest.raises(ValueError, match="Unknown signature"):
        resolve_signatures("blood", "cd8", ["typo"])
    with pytest.raises(ValueError, match="No signatures"):
        resolve_signatures("blood", "cd8", ["AIM"], ["AIM"])


def _ranking_tables():
    clones = pd.DataFrame({
        "donor": ["p1"] * 4 + ["p2"] * 4,
        "CDR3ab": list("ABCD") * 2,
        "eligible_for_review": [True, True, True, False] * 2,
    })
    rows = []
    for donor in ("p1", "p2"):
        for sample in ("s1", "s2"):
            for clone, x, y in zip("ABCD", [9, 8, 1, 99], [1, 8, 9, 99]):
                rows.append({"donor": donor, "sample": sample, "lineage": "cd8", "CDR3ab": clone,
                             "signature_X": x if sample == "s1" else y, "signature_Y": y,
                             "signature_X_percentile": 0.9, "signature_Y_percentile": 0.9,
                             "cells": 2, "frequency": 0.25})
    return clones, pd.DataFrame(rows)


def test_round_robin_covers_samples_deduplicates_and_uses_one_total_budget():
    clones, scores = _ranking_tables()
    out = select_round_robin(clones, scores, ["X", "Y"], 4)
    selected = out[out.selected_for_review]
    assert selected.groupby("donor").CDR3ab.agg(list).to_dict() == {"p1": ["A", "C"], "p2": ["A", "C"]}
    assert selected.selected_sample.tolist() == ["s1", "s2"] * 2
    assert selected.selection_rank.tolist() == [1, 2, 3, 4]
    # Input row order never breaks ties or changes selection.
    shuffled = select_round_robin(clones, scores.sample(frac=1, random_state=7), ["X", "Y"], 4)
    pd.testing.assert_frame_equal(out, shuffled)
    unlimited = select_round_robin(clones, scores, ["X", "Y"], 0)
    assert unlimited[unlimited.selected_for_review].groupby("donor").CDR3ab.agg(list).tolist() == [list("ACB")] * 2
    capped = select_round_robin(clones, scores, ["X", "Y"], 3)
    assert capped.selected_for_review.sum() == 3
    assert capped[capped.selected_for_review].donor.tolist() == ["p1", "p1", "p2"]


def test_round_robin_takes_top_of_distinct_signatures_and_handles_empty_lists():
    clones, scores = _ranking_tables()
    scores = scores[scores["sample"] == "s1"].copy()
    selected = select_round_robin(clones, scores, ["X", "Y"], 4)
    assert selected[selected.selected_for_review].selected_signature.tolist() == ["X", "X", "Y", "Y"]
    scores["signature_X"] = np.nan
    selected = select_round_robin(clones, scores, ["X", "Y"], 2)
    assert selected[selected.selected_for_review].CDR3ab.tolist() == ["C", "C"]
    selected = select_round_robin(clones, scores, ["X", "Y"], 10, quantile=0.95)
    assert not selected.selected_for_review.any()


def test_default_budget_is_200_across_patients():
    clones = pd.DataFrame({"donor": np.repeat(["p1", "p2"], 1000),
                           "CDR3ab": [str(i) for i in range(1000)] * 2,
                           "eligible_for_review": True})
    scores = clones.assign(sample="s", lineage="cd8", cells=2, frequency=1/1000,
                           signature_X=list(range(1, 1001)) * 2)
    scores["signature_X_percentile"] = scores.groupby("donor").signature_X.rank(pct=True)
    selected = select_round_robin(clones, scores, ["X"])
    assert selected.selected_for_review.sum() == 200
    assert selected[selected.selected_for_review].groupby("donor").size().to_dict() == {"p1": 100, "p2": 100}
    assert selected[selected.selected_for_review].selection_rank.tolist() == list(range(1, 201))


def test_exhausted_signature_yields_to_remaining_lists_without_padding():
    clones = pd.DataFrame({"donor": "p1", "CDR3ab": list("ABCDE"), "eligible_for_review": True})
    scores = clones.assign(sample="s", lineage="cd8", cells=2, frequency=0.2,
                           signature_X=[5, 0, -1, -2, -3], signature_Y=[1, 5, 4, 3, -1])
    for name in ("X", "Y"):
        scores[f"signature_{name}_percentile"] = scores[f"signature_{name}"].rank(pct=True)
    selected = select_round_robin(clones, scores, ["X", "Y"], 200, quantile=0)
    chosen = selected[selected.selected_for_review]
    assert chosen.CDR3ab.tolist() == list("ABCD")
    assert chosen.selected_signature.tolist() == ["X", "Y", "Y", "Y"]
    assert selected.loc[selected.CDR3ab.eq("E"), "selected_for_review"].tolist() == [False]


def test_cutoff_rejects_flat_nonfinite_negative_and_singleton_lists():
    clones, scores = _ranking_tables()
    scores["signature_X"] = 5.0
    assert not signature_pass_mask(scores, "X", quantile=0).any()
    scores["signature_X"] = [-1, -2, np.nan, -3] * 4
    assert not signature_pass_mask(scores, "X", quantile=0).any()
    assert signature_pass_mask(scores, "X", quantile=0, min_score=-2).sum() == 4
    singleton = scores.iloc[:1].copy()
    singleton["signature_X"] = 8.0
    assert not signature_pass_mask(singleton, "X", quantile=0).any()
    scores["signature_X"] = np.inf
    assert not signature_pass_mask(scores, "X", quantile=0).any()


def test_score_and_percentile_must_qualify_in_the_same_stratum():
    scores = pd.DataFrame({"donor": "p", "sample": ["s1", "s1", "s2", "s2"], "lineage": "cd8",
                           "CDR3ab": ["A", "B", "A", "B"],
                           "signature_X": [-1, -2, 1, 2], "signature_X_percentile": [1, .5, .5, 1]})
    assert signature_pass_mask(scores, "X").tolist() == [False, False, False, True]


def _args(*extra):
    return create_parser().parse_args(["prioritize", "samples.yaml", "-o", "out", *extra])


def test_cell_filters_and_low_umi_override_share_frequency_denominator():
    cells = phenotype_cells(_til_cells())
    cells.obs.loc["cell_0", "TRA_1_reads"] = 1
    cells.obs.loc["cell_1", "TRA_1_umis"] = 1
    cells.obs.loc["cell_2", "CD3"] = 0
    cells, counts, _ = filter_cells(cells, _args("--min-vdj-reads", "10", "--min-vdj-umis", "1"))
    assert cells.obs_names.tolist() == ["cell_1", "cell_3", "cell_4", "cell_5"]
    assert counts[-1]["after"] == 4
    clones = aggregate_clonotypes(cells, min_umi=1).set_index("CDR3ab")
    long = build_clone_sample_long(cells, min_umi=1).set_index("CDR3ab")
    assert clones.cell_count.to_dict() == long.cells.to_dict()
    assert long.frequency.sum() == 1
    assert clones.max_frequency.to_dict() == long.frequency.to_dict()


def test_expression_thresholds_use_log_cp10k_and_both_excludes_unknown():
    cells = _til_cells()
    cells.obs.loc["cell_0", ["CD4", "CD8"]] = 0
    cells = phenotype_cells(cells)
    kept, _, _ = filter_cells(cells, _args("--tcell-type", "both", "--min-expression", "GZMB=3"))
    assert kept.obs_names.tolist() == ["cell_1", "cell_2", "cell_3", "cell_6", "cell_7", "cell_8", "cell_9"]
    assert set(kept.obs.lineage) == {"cd4", "cd8"}
    low, _, _ = filter_cells(cells, _args("--max-expression", "GZMB=1"))
    assert low.obs_names.tolist() == ["cell_4", "cell_5"]
    with pytest.raises(KeyError, match="NOT_A_GENE"):
        filter_cells(cells, _args("--min-expression", "NOT_A_GENE=1"))


def test_requested_read_gate_fails_on_missing_read_data():
    cells = phenotype_cells(_til_cells())
    del cells.obs["TRA_1_reads"]
    with pytest.raises(ValueError, match="TRA_1_reads"):
        filter_cells(cells, _args("--min-vdj-reads", "1"))


def test_exact_gene_gates_handle_alleles_and_missing_calls():
    clones = pd.DataFrame({"alpha_v_gene": ["TRAV12-2*01", "TRAV12-1", None],
                           "beta_v_gene": ["TRBV1", "TRBV1", "TRBV1"]})
    assert gene_filter_mask(clones, [], ["trav12-2"], "v").tolist() == [False, True, False]
    assert gene_filter_mask(clones, ["TRAV12-2"], [], "v").tolist() == [True, False, False]
