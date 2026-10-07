"""Behavioral calibration checks against known nulls and planted signals."""

import numpy as np
import pandas as pd
import pytest

from tcrsift.prioritize import select_round_robin, signature_pass_mask
from tcrsift.signature_background import calibrate_signature_background


def _tables(values, size=10):
    values = np.asarray(values)
    cells = pd.DataFrame({
        "donor": "p", "sample": "s", "lineage": "cd8",
        "CDR3ab": np.repeat(np.arange(len(values) // size).astype(str), size),
        "signature_X": values, "signature_Y": 7 + 3 * values,
    })
    keys = ["donor", "sample", "lineage", "CDR3ab"]
    groups = cells.groupby(keys, observed=True)
    scores = groups[["signature_X", "signature_Y"]].mean().join(groups.size().rename("cells")).reset_index()
    for name in ("X", "Y"):
        scores[f"signature_{name}_percentile"] = scores[f"signature_{name}"].rank(pct=True)
    return scores, cells


def test_exact_small_population_does_not_force_a_positive():
    scores, cells = _tables([0, 0, 1, 1], size=2)
    result, audit = calibrate_signature_background(scores, cells, ["X"])
    assert audit.method.tolist() == ["exact"]
    assert audit.draws.tolist() == [6]
    assert result.signature_X_noise_floor.tolist() == [1, 1]
    assert result.signature_X_background_tail_probability.tolist() == [1, 1 / 6]
    assert not signature_pass_mask(result, "X", quantile=0).any()


def test_noise_only_population_has_small_tail_instead_of_fixed_top_decile():
    rng = np.random.default_rng(81)
    scores, cells = _tables(rng.normal(size=10000))
    result, audit = calibrate_signature_background(scores, cells, ["X", "Y"])
    # N(0,1) means over 10 cells have a 99th percentile near 2.326/sqrt(10).
    assert abs(audit.iloc[0].noise_floor - 2.326 / np.sqrt(10)) < .08
    passing = signature_pass_mask(result, "X")
    assert 1 <= passing.sum() <= 25  # ~1% per comparison, explicitly not FDR control
    # Native scale and offset cannot affect qualification.
    np.testing.assert_allclose(result.signature_Y_noise_floor, 7 + 3 * result.signature_X_noise_floor)
    assert passing.equals(signature_pass_mask(result, "Y"))


def test_detects_planted_clone_expression_with_matching_cell_count():
    rng = np.random.default_rng(19)
    values = rng.normal(size=10000)
    values[:100] += 4
    scores, cells = _tables(values)
    result, _ = calibrate_signature_background(scores, cells, ["X"])
    passing = signature_pass_mask(result, "X")
    assert passing[result.CDR3ab.astype(int) < 10].all()
    assert passing.sum() < 40


def test_floor_accounts_for_clone_size():
    scores, cells = _tables(np.random.default_rng(7).normal(size=1000), size=1)
    cells.loc[:99, "CDR3ab"] = "large"
    groups = cells.groupby(["donor", "sample", "lineage", "CDR3ab"], observed=True)
    scores = groups[["signature_X"]].mean().join(groups.size().rename("cells")).reset_index()
    result, audit = calibrate_signature_background(scores, cells, ["X"])
    floors = audit.set_index("clone_cells").noise_floor
    assert floors[100] < floors[1] / 5
    assert audit.set_index("clone_cells").loc[1, "method"] == "exact"
    assert result.loc[result.CDR3ab.eq("large"), "signature_X_noise_floor"].iloc[0] == floors[100]


def test_strata_are_isolated_and_input_order_and_other_signatures_do_not_change_draws():
    scores, cells = _tables(np.random.default_rng(41).normal(size=200))
    one, one_audit = calibrate_signature_background(scores, cells, ["X"])
    other_cells = cells.assign(donor="q", sample="other", lineage="cd4", signature_X=99)
    other_scores = scores.assign(donor="q", sample="other", lineage="cd4", signature_X=99)
    both, audit = calibrate_signature_background(
        pd.concat([scores, other_scores], ignore_index=True).sample(frac=1, random_state=2),
        pd.concat([cells, other_cells], ignore_index=True).sample(frac=1, random_state=3), ["Y", "X"],
    )
    for suffix in ("noise_floor", "background_tail_probability"):
        col = f"signature_X_{suffix}"
        pd.testing.assert_series_equal(one[col], both.loc[both.donor.eq("p"), col].sort_index())
    assert audit.loc[(audit.donor == "p") & (audit.signature == "X"), "noise_floor"].iloc[0] == one_audit.noise_floor.iloc[0]
    assert not signature_pass_mask(both[both.donor.eq("q")], "X", quantile=0).any()


@pytest.mark.parametrize("value", [0.0, .1, np.nan, np.inf])
def test_flat_and_unavailable_signatures_never_pass(value):
    scores, cells = _tables(np.repeat(value, 30), size=3)
    result, audit = calibrate_signature_background(scores, cells, ["X"])
    assert not signature_pass_mask(result, "X", quantile=0).any()
    if np.isfinite(value):
        assert result.signature_X_background_tail_probability.eq(1).all()
    else:
        assert result.signature_X_noise_floor.isna().all()
        assert audit.status.eq("unavailable").all()


def test_cell_population_must_match_clone_counts():
    scores, cells = _tables(np.arange(20))
    with pytest.raises(ValueError, match="match the cells"):
        calibrate_signature_background(scores, cells.iloc[1:], ["X"])


def test_background_floor_can_be_negative_and_manual_floor_is_additional():
    scores, cells = _tables(np.r_[np.repeat(-1, 10), np.repeat(-10, 990)])
    result, _ = calibrate_signature_background(scores, cells, ["X"])
    assert signature_pass_mask(result, "X").sum() == 1
    assert not signature_pass_mask(result, "X", min_score=0).any()
    assert not signature_pass_mask(result, "X", cutoff="legacy").any()


def test_ten_of_twenty_five_high_cells_qualify_despite_low_clone_mean():
    values = np.random.default_rng(52).normal(size=2500)
    values[:25] = 2  # a clone with a uniformly raised mean
    values[25:50] = [10] * 10 + [-10] * 15  # strong subset diluted in the mean
    scores, cells = _tables(values, size=25)
    result, _ = calibrate_signature_background(scores, cells, ["X"])
    mixed = result.CDR3ab.eq("1")
    assert result.loc[mixed, "signature_X_high_cells"].tolist() == [10]
    assert result.loc[mixed, "signature_X_high_fraction"].tolist() == [.4]
    assert not signature_pass_mask(result, "X")[mixed].any()
    assert signature_pass_mask(result, "X", evidence="subset")[mixed].all()
    assert signature_pass_mask(result, "X", evidence="either")[mixed].all()
    assert not signature_pass_mask(result, "X", evidence="subset", min_high_fraction=.5)[mixed].any()
    assert not signature_pass_mask(result, "X", evidence="subset", min_high_cells=11)[mixed].any()
    assert not signature_pass_mask(result, "X", evidence="subset", min_score=0)[mixed].any()
    clones = result[["donor", "CDR3ab"]].assign(eligible_for_review=True)
    selected = select_round_robin(clones, result.assign(frequency=.01), ["X"], 2, evidence="either")
    chosen = selected[selected.selected_for_review]
    assert chosen.CDR3ab.tolist() == ["0", "1"]
    assert chosen.selected_evidence.tolist() == ["mean", "subset"]
    assert chosen.iloc[1].selected_cells == 25
    assert chosen.iloc[1].selected_high_cells == 10
    assert chosen.iloc[1].selected_high_fraction == .4
    pd.testing.assert_frame_equal(selected, select_round_robin(
        clones, result.assign(frequency=.01).sample(frac=1, random_state=4), ["X"], 2, evidence="either",
    ))


def test_subset_tail_matches_exact_combinatorial_probability():
    from math import comb

    scores, cells = _tables(np.arange(20), size=5)
    result, audit = calibrate_signature_background(scores, cells, ["X"], high_quantile=.5)
    row = result.loc[result.CDR3ab.eq("3")].iloc[0]
    assert row.signature_X_high_cells == 5
    assert audit.background_high_cells.tolist() == [9]
    expected = comb(9, 5) / comb(20, 5)
    assert row.signature_X_subset_tail_probability == pytest.approx(expected)
    assert signature_pass_mask(result, "X", evidence="subset").sum() == 1


def test_subset_requires_several_high_cells_and_no_majority_is_required():
    values = np.random.default_rng(92).normal(size=1000)
    values[:10] = [100] + [-100] * 9
    scores, cells = _tables(values)
    result, _ = calibrate_signature_background(scores, cells, ["X"])
    assert result.loc[result.CDR3ab.eq("0"), "signature_X_high_cells"].tolist() == [1]
    assert not signature_pass_mask(result, "X", evidence="subset")[result.CDR3ab.eq("0")].any()


@pytest.mark.parametrize("value", [0., .1, np.nan, np.inf])
def test_flat_unavailable_subset_cannot_qualify(value):
    scores, cells = _tables(np.repeat(value, 100), size=10)
    result, _ = calibrate_signature_background(scores, cells, ["X"])
    assert not signature_pass_mask(result, "X", evidence="either", quantile=0).any()


def test_subset_evidence_cannot_silently_fall_back_to_legacy():
    scores, _ = _tables(np.arange(20), size=5)
    with pytest.raises(ValueError, match="Subset evidence requires"):
        signature_pass_mask(scores, "X", cutoff="legacy", evidence="either")


def test_subset_support_cannot_combine_count_and_enrichment_from_different_samples():
    scores = pd.DataFrame({
        "donor": "p", "sample": ["a", "b"], "lineage": "cd8", "CDR3ab": "A",
        "signature_X": [1, 2], "signature_X_high_cells": [10, 1],
        "signature_X_high_fraction": [.4, .1], "signature_X_subset_percentile": 1.,
        "signature_X_subset_tail_probability": [.2, .0001],
    })
    assert not signature_pass_mask(scores, "X", evidence="subset").any()
