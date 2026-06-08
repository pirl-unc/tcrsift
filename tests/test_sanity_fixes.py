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

"""Regression tests for the sanity-check bug-fix batch.

Each test pins a concrete bug fixed in this batch; most fail on the prior code.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from tcrsift.assemble import _pick_best_allele
from tcrsift.candidate import select_candidates
from tcrsift.filter import assign_tiers_threshold, filter_clonotypes_threshold
from tcrsift.report import _format_selection_item
from tcrsift.selection import select_freq_prism_per_condition
from tcrsift.signature_methods import call_positive
from tcrsift.unify import compute_condition_statistics, merge_experiments
from tcrsift.validation import (
    TCRsiftValidationError,
    canonicalize_clonotype_columns,
    validate_clonotype_df,
)


class TestConditionSubstringConflation:
    def test_prefix_substring_conditions_do_not_conflate(self):
        # pool1 must not also pull in condition_pool10's fractions (#sanity).
        df = pd.DataFrame({
            "CDR3_pair": ["a", "b"],
            "condition_pool1.frac": [0.5, 0.0],
            "condition_pool10.frac": [0.1, 0.2],
        })
        out = compute_condition_statistics(
            df, ["pool1", "pool10"], verbose=False, show_progress=False,
        )
        assert out["combined.pool1.frac.sum"].tolist() == [0.5, 0.0]
        assert out["combined.pool10.frac.sum"].tolist() == [0.1, 0.2]


class TestMergeExperimentsDuplicateName:
    def test_duplicate_name_raises(self):
        df = pd.DataFrame(
            {"CDR3_pair": ["a/b"], "CDR3_alpha": ["a"], "CDR3_beta": ["b"]}
        )
        with pytest.raises(TCRsiftValidationError, match="Duplicate experiment name"):
            merge_experiments([(df, "X"), (df, "X")], verbose=False)


class TestSelectCandidatesNaN:
    def test_nan_scored_clones_not_selected_when_pool_smaller_than_top_n(self):
        # Pool of 3 tier3 clones, only 1 with a real score, top_n=3.
        df = pd.DataFrame({
            "tier": ["tier3", "tier3", "tier3"],
            "signature_x": [0.9, np.nan, np.nan],
        })
        out = select_candidates(df, signatures=("x",), top_n=3)
        assert out["tier3plus_top_by_x_signature"].tolist() == [True, False, False]
        assert out["is_selected"].tolist() == [True, False, False]


class TestCallPositiveNaN:
    @pytest.mark.parametrize("method,kw", [
        ("quantile", {"quantile": 0.5}),
        ("otsu", {}),
        ("gap", {"search_top": 0.8}),
    ])
    def test_single_nan_does_not_empty_positive_set(self, method, kw):
        # A clear high cluster (5, 6) above a baseline (0s), with one NaN cell.
        scores = pd.Series([0.0, 0.0, 0.0, np.nan, 5.0, 6.0])
        pos = call_positive(scores, method=method, **kw)
        assert pos.sum() >= 1  # not silently emptied by the NaN
        assert not bool(pos.iloc[3])  # the NaN cell is never positive

    def test_all_nan_returns_all_false(self):
        scores = pd.Series([np.nan, np.nan])
        pos = call_positive(scores, method="quantile", quantile=0.5)
        assert pos.sum() == 0


class TestPickBestAlleleByFraction:
    def test_picks_higher_fraction_not_higher_count(self):
        # Allele A (len 12) matches 12/12 = 1.0; allele B (len 20) matches
        # 13/20 = 0.65 — B has the higher absolute count (13 > 12) but the
        # lower fraction. The picker must choose A.
        pool = {"A": "Y" * 12, "B": "Y" * 20}
        contig = "Y" * 13 + "X" * 7  # len 20
        allele, score, _ = _pick_best_allele(contig, "TRBC1", pool)
        assert allele == "A"
        assert score == pytest.approx(1.0)


class TestFilterMaxConditionsSentinel:
    def _df(self, n_conditions):
        n = len(n_conditions)
        return pd.DataFrame({
            "CDR3ab": [f"CARA{i}_CASB{i}" for i in range(n)],
            "CDR3_alpha": [f"CARA{i}" for i in range(n)],
            "CDR3_beta": [f"CASB{i}" for i in range(n)],
            "cell_count": [5] * n,
            "max_frequency": [0.1] * n,
            "n_conditions": n_conditions,
        })

    def test_default_none_disables_filter(self):
        out = filter_clonotypes_threshold(self._df([1, 9]), min_cells=2)
        assert len(out) == 2

    def test_explicit_value_applies(self):
        out = filter_clonotypes_threshold(self._df([1, 9]), min_cells=2, max_conditions=2)
        assert len(out) == 1

    def test_explicit_999_is_not_a_magic_disable(self):
        # A clone in 1000 conditions must be dropped by max_conditions=999
        # (old code treated 999 as "off" and kept it).
        out = filter_clonotypes_threshold(
            self._df([5, 1000]), min_cells=2, max_conditions=999,
        )
        assert len(out) == 1


class TestAssignTiersStrictnessOrder:
    def test_strictest_qualifying_tier_wins_for_nonlexical_names(self):
        # "loose" sorts before "strict" lexically but is the more permissive
        # tier; the clone qualifies for both and must get the stricter one.
        tier_defs = {
            "loose": {"min_cells": 2, "min_frequency": 0.0},
            "strict": {"min_cells": 10, "min_frequency": 0.01},
        }
        df = pd.DataFrame({"cell_count": [10], "max_frequency": [0.02]})
        out = assign_tiers_threshold(df, tier_definitions=tier_defs)
        assert out["tier"].iloc[0] == "strict"


def _prism_feat():
    rows = []
    for cond in ("A", "B"):
        for i in range(6):
            rows.append({
                "CDR3ab": f"{cond}_C{i}",
                "condition": cond,
                "frequency": 0.02 - i * 0.001,
                "ppost_alpha": 10 ** (-i),
                "ppost_beta": 10 ** (-i),
                "antigen_response_score": float(i),
                "naive_score": float(i),
                "coverage": float(10 - i),
            })
    return pd.DataFrame(rows)


class TestPrismRescueAndTieBreak:
    def test_rescue_marks_thin_condition_picks(self):
        feat = _prism_feat()
        # Gate so condition A keeps few clones, forcing rescue.
        out = select_freq_prism_per_condition(
            feat, condition_col="condition", gate=0.018,
            top_freq=4, top_prism=0, rescue_target=4, rescue_rank_col="coverage",
        )
        assert "rescued" in out.columns
        assert out["rescued"].any()  # at least one sub-threshold clone rescued

    def test_no_rescue_when_target_zero(self):
        feat = _prism_feat()
        out = select_freq_prism_per_condition(
            feat, condition_col="condition", gate=0.0, top_freq=3, top_prism=0,
        )
        assert not out["rescued"].any()

    def test_cdr3ab_only_tie_break_is_deterministic(self):
        feat = _prism_feat()
        base = select_freq_prism_per_condition(
            feat, condition_col="condition", gate=0.0, top_freq=3, top_prism=0,
            tie_break=["cdr3ab"],
        )
        shuffled = feat.sample(frac=1.0, random_state=1).reset_index(drop=True)
        out = select_freq_prism_per_condition(
            shuffled, condition_col="condition", gate=0.0, top_freq=3, top_prism=0,
            tie_break=["cdr3ab"],
        )
        pd.testing.assert_frame_equal(
            out.reset_index(drop=True), base.reset_index(drop=True)
        )


class TestCanonicalizeClonotypeColumns:
    @pytest.mark.parametrize("alias,canonical", [
        ("CDR3a", "CDR3_alpha"),
        ("CDR3b", "CDR3_beta"),
        ("cdr3_alpha", "CDR3_alpha"),
        ("TRB_cdr3", "CDR3_beta"),
        ("cdr3.beta", "CDR3_beta"),
        ("cdr3ab", "CDR3ab"),
    ])
    def test_aliases_renamed(self, alias, canonical):
        df = pd.DataFrame({alias: ["X"]})
        out = canonicalize_clonotype_columns(df)
        assert canonical in out.columns
        assert alias not in out.columns

    def test_does_not_overwrite_existing_canonical(self):
        df = pd.DataFrame({"CDR3_beta": ["canon"], "CDR3b": ["synonym"]})
        out = canonicalize_clonotype_columns(df)
        # The real canonical column is preserved; the synonym is left as-is.
        assert out["CDR3_beta"].tolist() == ["canon"]
        assert "CDR3b" in out.columns

    def test_validate_accepts_synonym_columns(self):
        # CDR3b alone is enough to pass validation once canonicalized.
        df = pd.DataFrame({"CDR3b": ["CASSL"], "cell_count": [3]})
        out = validate_clonotype_df(df, for_filtering=True)
        assert "CDR3_beta" in out.columns

    def test_no_synonyms_means_no_rename(self):
        # build_paired_id=False isolates the rename behavior: canonical-only
        # input has nothing to rename and is returned unchanged.
        df = pd.DataFrame({"CDR3_alpha": ["a"], "CDR3_beta": ["b"]})
        out = canonicalize_clonotype_columns(df, build_paired_id=False)
        assert list(out.columns) == ["CDR3_alpha", "CDR3_beta"]


class TestAddObsColumnsNoFragmentation:
    def test_adds_many_columns_without_performance_warning(self):
        import anndata as ad
        from pandas.errors import PerformanceWarning

        from tcrsift.loader import _add_obs_columns

        adata = ad.AnnData(np.zeros((4, 2), dtype="float32"))
        new_cols = {f"c{i}": np.arange(4) for i in range(150)}
        with warnings.catch_warnings():
            warnings.simplefilter("error", PerformanceWarning)
            _add_obs_columns(adata, new_cols)
        assert "c149" in adata.obs.columns
        assert adata.obs.shape[1] >= 150

    def test_overwrites_existing_and_broadcasts_scalars(self):
        import anndata as ad

        from tcrsift.loader import _add_obs_columns

        adata = ad.AnnData(np.zeros((3, 2), dtype="float32"))
        adata.obs["x"] = [1, 2, 3]
        _add_obs_columns(adata, {"x": [9, 9, 9], "label": "donorA"})
        assert adata.obs["x"].tolist() == [9, 9, 9]
        assert adata.obs["label"].tolist() == ["donorA"] * 3


class TestCanonicalBuildsPairedId:
    def test_builds_cdr3ab_from_alpha_beta_when_absent(self):
        df = pd.DataFrame({"CDR3_alpha": ["CAA", "CAB"], "CDR3_beta": ["CASA", "CASB"]})
        out = canonicalize_clonotype_columns(df)
        assert out["CDR3ab"].tolist() == ["CAA_CASA", "CAB_CASB"]

    def test_builds_from_synonyms(self):
        df = pd.DataFrame({"CDR3a": ["CAA"], "CDR3b": ["CASB"]})
        out = canonicalize_clonotype_columns(df)
        assert out["CDR3ab"].tolist() == ["CAA_CASB"]

    def test_does_not_clobber_existing_cdr3ab(self):
        df = pd.DataFrame({
            "CDR3_alpha": ["CAA"], "CDR3_beta": ["CASB"], "CDR3ab": ["existing"],
        })
        out = canonicalize_clonotype_columns(df)
        assert out["CDR3ab"].tolist() == ["existing"]

    def test_beta_only_not_paired(self):
        df = pd.DataFrame({"CDR3_beta": ["CASB"]})
        out = canonicalize_clonotype_columns(df)
        assert "CDR3ab" not in out.columns

    def test_opt_out(self):
        df = pd.DataFrame({"CDR3_alpha": ["a"], "CDR3_beta": ["b"]})
        out = canonicalize_clonotype_columns(df, build_paired_id=False)
        assert "CDR3ab" not in out.columns

    def test_does_not_mutate_caller(self):
        df = pd.DataFrame({"CDR3_alpha": ["a"], "CDR3_beta": ["b"]})
        canonicalize_clonotype_columns(df)
        assert "CDR3ab" not in df.columns


class TestFormatSelectionItem:
    @pytest.mark.parametrize("raw,expected", [
        ("AIM+=freq#6(0.90%)", "AIM+ freq#6 (0.90%)"),
        ("IFN+CTY-=freq#12(0.25%)", "IFN+CTY- freq#12 (0.25%)"),
        ("  CTY-=prism#3(0.10%)  ", "CTY- prism#3 (0.10%)"),
        ("already formatted", "already formatted"),
        ("", ""),
    ])
    def test_format(self, raw, expected):
        assert _format_selection_item(raw) == expected
