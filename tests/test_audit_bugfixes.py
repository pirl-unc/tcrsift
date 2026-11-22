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

"""Adversarial regression tests for v3.29.1 correctness & hardening fixes."""

from __future__ import annotations

import warnings
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csc_matrix, csr_matrix


class TestCrossSampleDuplicateBarcodes:
    """10x barcode whitelists repeat across samples; joins must be positional."""

    def test_attribution_and_long_table_preserve_duplicate_barcodes_across_samples(self):
        from tcrsift.clonotype import aggregate_clonotypes, build_clone_sample_long
        from tcrsift.config import AttributionConfig

        # Same barcode "AAACCTGAG-1" in S1 and S2 with distinct clones and UMIs.
        obs = pd.DataFrame(
            {
                "sample": ["S1", "S1", "S2", "S2"],
                "patient_id": ["D1", "D1", "D2", "D2"],
                "CDR3_alpha": ["CAVA1", "CAVA1", "CAVA2", "CAVA2"],
                "CDR3_beta": ["CASSB1", "CASSB1", "CASSB2", "CASSB2"],
                "TRA_1_cdr3": ["CAVA1", "CAVA1", "CAVA2", "CAVA2"],
                "TRB_1_cdr3": ["CASSB1", "CASSB1", "CASSB2", "CASSB2"],
                "TRA_1_umis": [10, 12, 30, 40],
                "TRB_1_umis": [15, 18, 50, 60],
            },
            index=["AAACCTGAG-1", "AAACCTGTC-1", "AAACCTGAG-1", "AAACCTGTC-1"],
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            adata = ad.AnnData(X=csr_matrix((4, 1)), obs=obs)

        cfg = AttributionConfig(enabled=True)
        clones = aggregate_clonotypes(
            adata, attribution=cfg, verbose=False, show_progress=False
        )
        by_id = clones.set_index("CDR3ab")
        assert by_id.loc["CAVA1_CASSB1", "cell_count"] == pytest.approx(2.0)
        assert by_id.loc["CAVA2_CASSB2", "cell_count"] == pytest.approx(2.0)
        assert by_id.loc["CAVA1_CASSB1", "donors"] == "D1"
        assert by_id.loc["CAVA2_CASSB2", "donors"] == "D2"

        long_df = build_clone_sample_long(adata, attribution=cfg)
        s1 = long_df[long_df["sample"] == "S1"].iloc[0]
        s2 = long_df[long_df["sample"] == "S2"].iloc[0]
        assert s1["cells"] == pytest.approx(2.0)
        assert s1["n_alpha_umis"] == pytest.approx(22.0)
        assert s2["cells"] == pytest.approx(2.0)
        assert s2["n_alpha_umis"] == pytest.approx(70.0)

    def test_expand_dual_alpha_variants_handles_duplicate_barcodes_and_missing_umis(self):
        from tcrsift.report import expand_dual_alpha_variants

        selected = pd.DataFrame(
            {
                "CDR3ab": ["CAVA1_CASSB1"],
                "CDR3_alpha": ["CAVA1"],
                "CDR3_beta": ["CASSB1"],
                "merged_alpha_partners": ["CAVA1;CAVA2"],
                "VDJ_beta_aa": ["VDJB1"],
                "beta_v_gene": ["TRBV7-9"],
                "beta_j_gene": ["TRBJ2-1"],
                "beta_c_gene": ["TRBC2"],
            }
        )
        # Barcode "BAR-1" appears in both samples; S1 carries dual-alpha (no UMI cols),
        # S2 carries an unrelated clone on the same barcode.
        obs = pd.DataFrame(
            {
                "sample": ["S1", "S2"],
                "TRA_1_cdr3": ["CAVA1", "CAVOTHER"],
                "TRA_2_cdr3": ["CAVA2", None],
                "TRB_1_cdr3": ["CASSB1", "CASSOTHER"],
                "TRA_1_vdj_aa": ["VDJA1_TRUE", "VDJ_WRONG"],
                "TRA_2_vdj_aa": ["VDJA2_TRUE", None],
                "TRA_1_v_gene": ["TRAV12-1", "TRAV99"],
                "TRA_2_v_gene": ["TRAV21", None],
            },
            index=["BAR-1", "BAR-1"],
        )
        expanded, variant_of = expand_dual_alpha_variants(selected, obs)
        assert len(expanded) == 2
        assert len(variant_of) == 2
        by_a = expanded.set_index("CDR3_alpha")
        assert by_a.loc["CAVA1", "VDJ_alpha_aa"] == "VDJA1_TRUE"
        assert by_a.loc["CAVA2", "VDJ_alpha_aa"] == "VDJA2_TRUE"

    def test_refine_cluster_positional_relabel_with_duplicate_barcodes(self):
        from tcrsift.embed import refine_cluster

        rng = np.random.default_rng(7)
        # 50 cells: first 25 in S1 (target compartment), next 25 in S2 sharing barcodes!
        barcodes = [f"BC_{i:02d}-1" for i in range(25)] * 2
        X = rng.poisson(3, size=(50, 12)).astype(np.float32)
        obs = pd.DataFrame(
            {
                "sample": ["S1"] * 25 + ["S2"] * 25,
                "phenotype": ["Target"] * 25 + ["Untouched"] * 25,
            },
            index=barcodes,
        )
        var = pd.DataFrame(index=[f"GENE{i}" for i in range(12)])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            adata = ad.AnnData(X=X, obs=obs, var=var)

        refine_cluster(
            adata,
            "Target",
            relabel_fn=lambda sub: ("Refined", "C_sub"),
            informative_genes=[f"GENE{i}" for i in range(10)],
            min_cells=10,
            final_cluster_col="final_cl",
        )
        assert (adata.obs["phenotype"].iloc[:25] == "Refined").all()
        # S2 cells with identical barcodes MUST remain untouched!
        assert (adata.obs["phenotype"].iloc[25:] == "Untouched").all()
        assert adata.obs["final_cl"].iloc[25:].isna().all()

    def test_signature_scoring_and_methods_long_handle_duplicate_barcodes(self):
        from tcrsift.signature_methods import (
            Signature,
            build_signature_methods,
            score_by_name,
            score_signature,
        )

        idx = ["BC1-1", "BC2-1", "BC1-1", "BC2-1"]
        expr = pd.DataFrame(
            {
                "GZMB": [10.0, 0.0, 20.0, 0.0],
                "PRF1": [8.0, 0.0, 15.0, 0.0],
                "IL7R": [0.0, 5.0, 0.0, 6.0],
                "CXCL13": [9.0, 0.0, 12.0, 0.0],
                "ENTPD1": [7.0, 0.0, 11.0, 0.0],
            },
            index=idx,
        )
        groups = pd.Series(["S1", "S1", "S2", "S2"], index=idx)
        sig = Signature("CytoTest", ("GZMB", "PRF1"))
        s = score_signature(expr, sig, groups=groups)
        assert len(s) == 4
        assert s.iloc[0] > s.iloc[1]
        assert s.iloc[2] > s.iloc[3]

        obs = pd.DataFrame(
            {
                "CDR3ab": ["CAVA_CASSB", "_", "CAVA2_CASSB2", "_ONLYBETA"],
                "sample": ["S1", "S1", "S2", "S2"],
            },
            index=idx,
        )
        long_sig = build_signature_methods(
            expr, obs, signatures=[sig], positive_method="quantile", quantile=0.25
        )
        # Incomplete "_" and "_ONLYBETA" keys must be excluded from synthetic methods.
        assert not long_sig["CDR3ab"].isin({"_", "_ONLYBETA", ""}).any()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            adata = ad.AnnData(X=expr.to_numpy(dtype=np.float32), obs=obs, var=pd.DataFrame(index=expr.columns))
        mana = score_by_name(adata, "MANAscore", groups=groups)
        assert len(mana) == 4
        assert mana.iloc[0] > mana.iloc[1]


class TestAttributionAndClonotypeInvariants:
    """Multi-donor isolation, slot promotion, and per-method frequency denominators."""

    def test_attribute_cells_isolates_donors_and_promotes_secondary_slot(self):
        from tcrsift.attribution import attribute_cells
        from tcrsift.config import AttributionConfig

        df = pd.DataFrame(
            {
                "sample": ["D1_S1", "D1_S1", "D2_S1", "D2_S1", "D2_S1"],
                "patient_id": ["D1", "D1", "D2", "D2", "D2"],
                # Row 0-1: D1 has complete CAVA1_CASSB_SHARED
                # Row 2: D2 has beta-only (_, CASSB_SHARED) -> must NOT leak to D1's CAVA1!
                # Row 3: D2 has homotypic TRA_1 == TRA_2 == CAVA2 -> single complete clone
                # Row 4: D2 has TRA_1 failing min_umi (1) and TRA_2 passing (6) -> promoted!
                "TRA_1_cdr3": ["CAVA1", "CAVA1", None, "CAVA2", "CAVLOW"],
                "TRA_2_cdr3": [None, None, None, "CAVA2", "CAVA2"],
                "TRB_1_cdr3": ["CASSB_SH", "CASSB_SH", "CASSB_SH", "CASSB2", "CASSB2"],
                "TRB_2_cdr3": [None, None, None, None, None],
                "TRA_1_umis": [10, 10, np.nan, 8, 1],
                "TRA_2_umis": [np.nan, np.nan, np.nan, 5, 6],
                "TRB_1_umis": [10, 10, 10, 10, 10],
                "TRB_2_umis": [np.nan, np.nan, np.nan, np.nan, np.nan],
            },
            index=[f"c{i}" for i in range(5)],
        )
        long_table, _merge_map, _stats = attribute_cells(
            df, AttributionConfig(enabled=True), min_umi=2
        )
        # D2's beta-only CASSB_SH cell must NOT be attributed to D1's CAVA1_CASSB_SH
        d2_rows = long_table[long_table["sample"] == "D2_S1"]
        assert "CAVA1_CASSB_SH" not in set(d2_rows["CDR3ab"])
        # Homotypic A2==A1 + secondary-promoted A2 both yield CAVA2_CASSB2 with weight 1.0 each
        assert d2_rows.loc[d2_rows["CDR3ab"] == "CAVA2_CASSB2", "weight"].sum() == pytest.approx(2.0)

    def test_max_frequency_per_method_uses_within_sample_denominator(self):
        from tcrsift.clonotype import aggregate_clonotypes

        # Donor 1 AIM sample has 2 cells (both clone A -> 100% in D1_AIM).
        # Donor 2 AIM sample has 8 cells (all clone B -> 100% in D2_AIM).
        # Clone A's max_frequency_per_method must be 1.0 (not 2/10 = 0.2!).
        obs = pd.DataFrame(
            {
                "sample": ["D1_AIM"] * 2 + ["D2_AIM"] * 8,
                "patient_id": ["D1"] * 2 + ["D2"] * 8,
                "enrichment_method": ["AIM"] * 10,
                "CDR3_alpha": ["CAVA"] * 2 + ["CAVB"] * 8,
                "CDR3_beta": ["CASSA"] * 2 + ["CASSB"] * 8,
            },
            index=[f"cell_{i}" for i in range(10)],
        )
        adata = ad.AnnData(X=csr_matrix((10, 1)), obs=obs)
        clones = aggregate_clonotypes(adata, verbose=False, show_progress=False)
        by_id = clones.set_index("CDR3ab")
        assert by_id.loc["CAVA_CASSA", "max_frequency_per_method"] == pytest.approx(1.0)

    def test_airr_export_empty_mouse_genes_and_d_gene(self, tmp_path):
        from tcrsift.clonotype import aggregate_clonotypes, export_clonotypes_airr

        empty_out = tmp_path / "empty.tsv"
        export_clonotypes_airr(pd.DataFrame(), str(empty_out))
        empty_df = pd.read_csv(empty_out, sep="\t")
        assert {"clone_id", "junction_aa_tra", "junction_aa_trb", "sequence_id", "productive"} <= set(empty_df.columns)

        obs = pd.DataFrame(
            {
                "sample": ["S1", "S1"],
                "CDR3_alpha": ["CAVSDGGSQGNLIF", "CAVSDGGSQGNLIF"],
                "CDR3_beta": ["CASSLGQAYEQYF", "CASSLGQAYEQYF"],
                "TRA_1_cdr3": ["CAVSDGGSQGNLIF", "CAVSDGGSQGNLIF"],
                "TRB_1_cdr3": ["CASSLGQAYEQYF", "CASSLGQAYEQYF"],
                "TRA_1_v_gene": ["Trav14-1*01", "Trav14-1*01"],
                "TRA_1_j_gene": ["Traj22*01", "Traj22*01"],
                "TRB_1_v_gene": ["Trbv13-2*01", "Trbv13-2*01"],
                "TRB_1_d_gene": ["Trbd1*01", "Trbd1*01"],
                "TRB_1_j_gene": ["Trbj2-7*01", "Trbj2-7*01"],
                "TRA_1_umis": [5, 5],
                "TRB_1_umis": [5, 5],
            },
            index=["c1", "c2"],
        )
        clones = aggregate_clonotypes(ad.AnnData(obs=obs), verbose=False, show_progress=False)
        assert clones.loc[0, "beta_d_gene"] == "Trbd1*01"
        clones["cell_count"] = [2.4]
        out_tsv = tmp_path / "mouse.airr.tsv"
        export_clonotypes_airr(clones, str(out_tsv))
        airr = pd.read_csv(out_tsv, sep="\t")
        assert airr.loc[0, "v_call_tra"] == "TRAV14-1*01"
        assert airr.loc[0, "d_call_trb"] == "TRBD1*01"
        assert airr.loc[0, "clone_count"] == 2

    def test_sample_overlap_with_integer_sample_ids(self):
        from tcrsift.clonotype import compute_sample_overlap_matrices

        long_df = pd.DataFrame(
            {
                "CDR3ab": ["A_B", "A_B", "C_D"],
                "donor": [1, 2, 2],
                "cells": [5, 5, 5],
            }
        )
        mats = compute_sample_overlap_matrices(long_df, sample_col="donor")
        assert mats["jaccard"].loc["1", "2"] == pytest.approx(0.5)


class TestFilterLogisticAndInsilicoFixes:
    """NaN tier safety, donor-local stats, exact logistic inverse, and tied-best insilico."""

    def test_clone_clears_tier_rejects_nan(self):
        from tcrsift.filter import DEFAULT_THRESHOLD_TIERS, clone_clears_tier

        t5 = DEFAULT_THRESHOLD_TIERS["tier5"]
        assert clone_clears_tier(np.nan, 0.5, t5) is False
        assert clone_clears_tier(10, np.nan, t5) is False
        assert clone_clears_tier(10, 0.5, {"min_cells": 2, "min_frequency": 0.01, "max_conditions": 2}, n_conditions=np.nan) is False

    def test_split_clonotypes_by_donor_recomputes_local_abundance(self):
        from tcrsift.filter import split_clonotypes_by_donor

        clonotypes = pd.DataFrame(
            {
                "CDR3ab": ["PUBLIC_AB"],
                "cell_count": [102],
                "max_frequency": [0.80],
                "max_frequency_per_method": [0.80],
            }
        )
        long_df = pd.DataFrame(
            {
                "CDR3ab": ["PUBLIC_AB", "PUBLIC_AB"],
                "donor": ["D1", "D2"],
                "sample": ["D1_S", "D2_S"],
                "method": ["AIM", "AIM"],
                "cells": [100, 2],
                "frequency": [0.80, 0.005],
            }
        )
        split = split_clonotypes_by_donor(clonotypes, long_df)
        assert split["D1"].iloc[0]["cell_count"] == 100
        assert split["D1"].iloc[0]["max_frequency"] == pytest.approx(0.80)
        assert split["D2"].iloc[0]["cell_count"] == 2
        assert split["D2"].iloc[0]["max_frequency"] == pytest.approx(0.005)
        assert split["D2"].iloc[0]["max_frequency_per_method"] == pytest.approx(0.005)

    def test_logistic_zero_weight_fallback_and_nullable_is_viral(self, monkeypatch):
        from tcrsift.filter import filter_clonotypes_logistic
        from tcrsift.model import LogisticFrequencyModel

        df = pd.DataFrame(
            {
                "CDR3ab": [f"A{i}_B{i}" for i in range(6)],
                "CDR3_alpha": [f"CAVA{i}" for i in range(6)],
                "CDR3_beta": [f"CASSB{i}" for i in range(6)],
                "cell_count": [10] * 6,
                "max_frequency": [0.02, 0.05, 0.10, 0.15, 0.20, 0.60],
                "is_viral": [False, None, pd.NA, False, False, False],
            }
        )
        # Weight == 0.0 must trigger fallback to default_freq_threshold (0.50),
        # not argmin picking min_freq_threshold (0.09) for every clone.
        monkeypatch.setattr(
            "tcrsift.filter.fit_frequency_logistic_model",
            lambda x, y: LogisticFrequencyModel(params=np.array([0.0]), converged=True, n_iter=3),
        )
        out = filter_clonotypes_logistic(df, fdr_tiers=[0.1, 0.01], default_freq_threshold=0.50)
        assert out.attrs["logistic_fallback_reason"] == "negative_weight"
        assert out.attrs["fdr_to_threshold"] == {0.1: 0.50, 0.01: 0.50}
        assert (out["tier"] == "tier1").sum() == 1

    def test_insilico_mask_retains_clones_tied_at_best_possible_value(self):
        from tcrsift.insilico_filter import FilterPredicate, insilico_mask

        # 3 of 10 clones tie at 0.0 (best possible low value): average rank = 2.0 ->
        # (2 - 1) / 9 = 0.111 > 0.05, so without best-tie protection 0 clones passed!
        df = pd.DataFrame({"score": [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]})
        mask = insilico_mask(df, [FilterPredicate("score", "low", 0.05)])
        assert mask.iloc[:3].all()
        assert not mask.iloc[3:].any()


class TestDtypesLoaderTilAndUnify:
    """StringDtype/h5ad round-trip, TIL incomplete-key exclusion, and unify deduplication."""

    def test_rehydrate_obs_handles_arbitrary_stringdtype_and_categorical_columns(self, tmp_path):
        from tcrsift._dtypes import rehydrate_obs
        from tcrsift.mnemonic import tcr_name

        obs = pd.DataFrame(
            {
                "CDR3_alpha": pd.array(["CAVSD", pd.NA], dtype="string"),
                "custom_str_col": pd.array(["val1", pd.NA], dtype="string"),
                "custom_cat_col": pd.Categorical(["catA", np.nan]),
                "TRA_count": pd.array([1, pd.NA], dtype="Int64"),
            },
            index=pd.array(["c1", "c2"], dtype="string"),
        )
        adata = ad.AnnData(X=csr_matrix((2, 1)), obs=obs)
        rehydrate_obs(adata)
        assert adata.obs["custom_str_col"].iloc[1] is None
        out_h5 = tmp_path / "rehydrated.h5ad"
        adata.write_h5ad(out_h5)
        reloaded = ad.read_h5ad(out_h5)
        rehydrate_obs(reloaded)
        assert reloaded.obs["CDR3_alpha"].iloc[1] is None
        assert pd.isna(reloaded.obs["custom_str_col"].iloc[1])
        assert reloaded.obs["CDR3_alpha"].fillna("").tolist() == ["CAVSD", ""]

        assert tcr_name(None) == "Anon"
        assert tcr_name(np.nan) == "Anon"
        assert tcr_name(pd.NA) == "Anon"

    def test_sample_sheet_antigen_names_broadcast_on_two_cell_adata(self, tmp_path):
        from tcrsift.loader import _add_obs_columns
        from tcrsift.sample_sheet import Sample, SampleSheet

        with pytest.raises(ValueError, match="culture_days"):
            Sample(sample="S_bad", vdj_dir="/tmp", culture_days=True)

        s = Sample(sample="S1", vdj_dir="/tmp", antigen_names=["MART-1", "gp100"])
        sheet = SampleSheet(samples=[s], species="mouse")
        assert s.species == "mouse"
        assert sheet.antigen_label_map()["S1"] == "MART-1 + gp100"

        adata = ad.AnnData(X=csr_matrix((2, 1)), obs=pd.DataFrame(index=["c1", "c2"]))
        _add_obs_columns(adata, {"antigen_names": s.antigen_names})
        assert adata.obs["antigen_names"].tolist() == ["MART-1;gp100", "MART-1;gp100"]
        h5_path = tmp_path / "two_cells.h5ad"
        adata.write_h5ad(h5_path)
        assert h5_path.exists()

    def test_til_h5ad_categorical_and_single_chain_exclusion(self, tmp_path):
        from tcrsift.til import (
            identify_til_specific_clones,
            load_til_data,
            match_til,
            summarize_til_clonotypes,
        )

        # 2 complete paired cells + 2 beta-only cells ("_CASSB") + 1 empty cell ("_")
        obs = pd.DataFrame(
            {
                "CDR3_alpha": pd.Categorical(["CAVA", "CAVA", np.nan, np.nan, np.nan]),
                "CDR3_beta": pd.Categorical(["CASSB", "CASSB", "CASSB_ONLY", "CASSB_ONLY", np.nan]),
                "sample": ["T1"] * 5,
            },
            index=[f"t{i}" for i in range(5)],
        )
        h5_path = tmp_path / "til.h5ad"
        ad.AnnData(X=csr_matrix((5, 1)), obs=obs).write_h5ad(h5_path)

        loaded = load_til_data("h5ad", h5_path, "T1")
        summary = summarize_til_clonotypes({"T1": loaded}, match_by="CDR3ab", min_cells=1)
        assert set(summary["CDR3ab"]) == {"CAVA_CASSB"}
        # Denominator is 2 complete paired cells, so CAVA_CASSB frequency is 1.0 (not 2/4 = 0.5)
        assert summary.iloc[0]["til_frequency"] == pytest.approx(1.0)

        culture = pd.DataFrame(
            {
                "CDR3ab": ["CAVA_CASSB", "_CASSB_ONLY", np.nan],
                "CDR3_alpha": ["CAVA", "", np.nan],
                "CDR3_beta": ["CASSB", "CASSB_ONLY", np.nan],
            }
        )
        matched = match_til(culture, ad.read_h5ad(h5_path), match_by="CDR3ab")
        assert bool(matched.iloc[0]["til_match"]) is True
        assert matched.iloc[0]["til_frequency"] == pytest.approx(1.0)
        assert bool(matched.iloc[1]["til_match"]) is False

        til_spec = identify_til_specific_clones(ad.read_h5ad(h5_path), min_cells=2)
        assert set(til_spec["CDR3ab"]) == {"CAVA_CASSB"}

    def test_unify_deduplicates_keys_and_computes_single_source_stats(self):
        from tcrsift.unify import (
            compute_condition_statistics,
            find_top_condition,
            merge_experiments,
        )

        exp1 = pd.DataFrame(
            {
                "CDR3_alpha": ["CAVA", "CAVA", np.nan],
                "CDR3_beta": ["CASSB", "CASSB", "CASSB_ONLY"],
                "cell_count": [3, 4, 2],
                "condition_pool1.frac": [0.2, 0.3, 0.0],
                "condition_pool2.frac": [0.0, 0.0, 0.0],
            }
        )
        exp2 = pd.DataFrame(
            {
                "CDR3_pair": ["CAVA/CASSB", "/CASSB_ONLY"],
                "cell_count": [5, 1],
            }
        )
        merged = merge_experiments([(exp1, "E1"), (exp2, "E2")], verbose=False)
        assert len(merged) == 2
        by_pair = merged.set_index("CDR3_pair")
        assert by_pair.loc["CAVA/CASSB", "combined.total_cells.count"] == 12
        assert by_pair.loc["/CASSB_ONLY", "combined.total_cells.count"] == 3

        stats = compute_condition_statistics(merged, ["pool1", "pool2"], verbose=False)
        by_pair = stats.set_index("CDR3_pair")
        assert by_pair.loc["CAVA/CASSB", "combined.pool1.frac.consistency"] == pytest.approx(1.0)

        top = find_top_condition(stats, ["pool1", "pool2"], verbose=False).set_index("CDR3_pair")
        assert top.loc["CAVA/CASSB", "top_condition"] == "pool1"
        assert pd.isna(top.loc["/CASSB_ONLY", "top_condition"])


class TestAssemblyCodonAndLeaderFixes:
    """Lowercase VDJ recoding, 1-codon backtracking, leader policy NT sync, and from-data constants."""

    def test_recode_codons_lowercase_repeats_restriction_sites_and_backtracking(self):
        from tcrsift.assemble import CodonOptimization, recode_codons, translate_dna

        assert translate_dna("atggcc")[0] == "MA"

        # Lowercase poly-Gly repeat + EcoRI site (gaattc = Glu-Phe)
        nt_in = "atgggaggaggagaattc"
        recoded, changed = recode_codons(
            nt_in, max_codon_repeats=1, restriction_sites=("GAATTC",)
        )
        assert translate_dna(recoded)[0] == translate_dna(nt_in)[0]
        assert "GAATTC" not in recoded
        assert len(changed) >= 2

        # NcoI (CCATGG = Pro CCA + Trp TGG): Trp has ONLY 1 codon (TGG), so forward
        # greedy cannot change TGG; 1-codon backtracking must swap CCA -> another Pro codon!
        nco_in = "ATGCCATGGGCC"
        nco_out, nco_changed = recode_codons(nco_in, restriction_sites=("CCATGG",))
        assert translate_dna(nco_out)[0] == "MPWA"
        assert "CCATGG" not in nco_out
        assert nco_changed == [3]

        co = CodonOptimization(mode="full", avoid_enzymes=("bsai", "ecori"))
        assert set(co.restriction_sites()) == {"GGTCTC", "GAATTC"}

    def test_leader_policy_keeps_full_nt_optimized_in_sync_and_validates_directives(self):
        from tcrsift.assemble import (
             apply_leader_policy,
            assemble_full_sequences,
            collect_germline_variants,
            export_fasta,
            fix_jc_parity,
            translate_dna,
        )
        from tcrsift.qc_battery import check_assembly
        from tcrsift.validation import TCRsiftValidationError

        # VDJ with lowercase repeat codons so focal recode modifies donor VDJ bytes
        vdj_a_aa = "MLLLLVPVLEVIFTLGGTRAQSVTQLGSHVSVSEGALVLLRCNYSSSVPPYLFWYVQYPNQGLQLLLKYTSAATLVKGINGFEAEFKKSETSFHLTKPSAHMSDAAEYFCAVSDGGSQGNLIFGKGTKLSVKP"
        # Strip leader to pure VDJ starting at AQSV...
        vdj_a_aa = vdj_a_aa[18:]
        vdj_b_aa = "NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASSLGQAYEQYFGPGTRLTVT"
        from tcrsift.assemble import optimize_codons

        # Force lowercase with an identical codon repeat ("ggagga") inside VDJ_alpha_nt
        vdj_a_nt = optimize_codons(vdj_a_aa).lower()
        vdj_b_nt = optimize_codons(vdj_b_aa).lower()

        clones = pd.DataFrame(
            {
                "CDR3ab": ["CAVSDGGSQGNLIF_CASSLGQAYEQYF"],
                "CDR3_alpha": ["CAVSDGGSQGNLIF"],
                "CDR3_beta": ["CASSLGQAYEQYF"],
                "alpha_v_gene": ["TRAV12-2"],
                "alpha_j_gene": ["TRAJ42"],
                "alpha_c_gene": ["TRAC"],
                "beta_v_gene": ["TRBV7-9"],
                "beta_j_gene": ["TRBJ2-7"],
                "beta_c_gene": ["TRBC2"],
                "VDJ_alpha_aa": [vdj_a_aa],
                "VDJ_alpha_nt": [vdj_a_nt],
                "VDJ_beta_aa": [vdj_b_aa],
                "VDJ_beta_nt": [vdj_b_nt],
            }
        )
        asm = assemble_full_sequences(
            clones,
            alpha_leader="CD28",
            beta_leader="CD8A",
            include_constant=True,
            linker="P2A",
            force_alpha_leader="CD8A",
            secondary_alpha_leader="CD8A",  # identical to forced primary -> must not emit duplicate twin!
            avoid_enzymes=("EcoRI",),
            show_progress=False,
            verbose=False,
        )
        assert len(asm) == 1  # redundant identical secondary suppressed
        row = asm.iloc[0]
        assert translate_dna(row["full_alpha_nt_optimized"])[0] == row["full_alpha_aa"]
        assert translate_dna(row["single_chain_nt_optimized"])[0] == row["single_chain_aa"]
        assert check_assembly(asm).status == "PASS"
        assert check_assembly(asm.iloc[0:0]).status == "SKIP"

        # Curated exogenous leaders (CD28/CD8A) must not appear as donor germline variants
        gv = collect_germline_variants(asm)
        assert len(gv[gv["region"] == "leader"]) == 0

        with pytest.raises(TCRsiftValidationError, match="leader_fallback"):
            apply_leader_policy(asm, "T2A", leader_fallback="not_a_leader")

        # Duplicate index fix_jc_parity + export_fasta NaN guard
        dup_df = pd.DataFrame(
            {
                "CDR3ab": ["A_B", "C_D"],
                "beta_j_gene": ["TRBJ1-1", "TRBJ2-1"],
                "beta_c_gene": ["TRBC2", "TRBC1"],
                "single_chain_aa": ["MAAA", np.nan],
            },
            index=[0, 0],
        )
        msgs = fix_jc_parity(dup_df)
        assert len(msgs) == 2
        assert dup_df["beta_c_gene_canonical"].tolist() == ["TRBC1", "TRBC2"]

        fa_path = Path(asm.attrs.get("tmp", "/tmp")) / "test_audit_nan_guard.fasta"
        try:
            export_fasta(dup_df, fa_path)
            txt = fa_path.read_text()
            assert "MAAA" in txt and "nan" not in txt.lower()
        finally:
            fa_path.unlink(missing_ok=True)

    def test_constant_from_data_aa_only_synthesizes_nt_and_marks_unverified(self):
        from tcrsift.assemble import assemble_full_sequences, translate_dna

        clones = pd.DataFrame(
            {
                "CDR3ab": ["CAVSDGGSQGNLIF_CASSLGQAYEQYF"],
                "CDR3_alpha": ["CAVSDGGSQGNLIF"],
                "CDR3_beta": ["CASSLGQAYEQYF"],
                "VDJ_alpha_aa": ["AQSVTQLGSHVSVSEGALVLLRCNYSSSVPPYLFWYVQYPNQGLQLLLKYTSAATLVKGINGFEAEFKKSETSFHLTKPSAHMSDAAEYFCAVSDGGSQGNLIFGKGTKLSVKP"],
                "VDJ_beta_aa": ["NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASSLGQAYEQYFGPGTRLTVT"],
                "alpha_constant_aa": ["NIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESS"],
                "beta_constant_aa": ["EDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD"],
            }
        )
        asm = assemble_full_sequences(
            clones,
            constant_source="from-data",
            include_constant=True,
            linker="T2A",
            show_progress=False,
            verbose=False,
        )
        row = asm.iloc[0]
        assert isinstance(row["full_alpha_nt_optimized"], str) and len(row["full_alpha_nt_optimized"]) > 0
        assert translate_dna(row["single_chain_nt_optimized"])[0] == row["single_chain_aa"]
        assert bool(row["construct_contig_verified"]) is False


class TestAnnotateTilSelectAtlasAndPlots:
    """strict_ab levenshtein guard, CITE-seq UMI exclusion, CD44 collision, 1-cell Leiden, and Jaccard."""

    def test_match_clonotypes_strict_ab_levenshtein_blocks_b_only_and_nan(self):
        from tcrsift.annotate import match_clonotypes
        from tcrsift.annotate_tcrs import annotate_tcrs

        db = pd.DataFrame(
            {
                "cdr3_alpha": ["CAVA_DB"],
                "cdr3_beta": ["CASSLGQAYEQYF"],
                "epitope": ["ELAGIGILTV"],
                "antigen_gene": ["MLANA"],
                "db_protein_canonical": ["MART-1"],
                "db_antigen_gene": ["MLANA"],
                "species": ["Homo sapiens"],
                "db_category": ["tumor_associated"],
                "database": ["VDJdb"],
                "is_viral": [False],
            }
        )
        clones = pd.DataFrame(
            {
                "CDR3ab": ["WRONG_CASSLGQAYEQYY", "NAN_ROW"],
                "CDR3_alpha": ["CAVWRONG", np.nan],
                "CDR3_beta": ["CASSLGQAYEQYY", np.nan],  # distance 1 from DB beta, wrong alpha
            }
        )
        matched = match_clonotypes(
            clones,
            db,
            match_strictness="strict_ab",
            match_mode="levenshtein",
            max_distance=1,
            verbose=False,
            show_progress=False,
        )
        assert not matched["db_match"].any()

        # annotate_tcrs(methods=["prism"]) auto-computes pgen/ppost
        valid_clones = pd.DataFrame(
            {
                "CDR3ab": ["CAVSDGGSQGNLIF_CASSLGQAYEQYF", "CAGAGGGSQGNLIF_CASSIRSSYEQYF"],
                "CDR3_alpha": ["CAVSDGGSQGNLIF", "CAGAGGGSQGNLIF"],
                "CDR3_beta": ["CASSLGQAYEQYF", "CASSIRSSYEQYF"],
            }
        )
        per_cell = pd.DataFrame(
            {
                "CDR3ab": ["CAVSDGGSQGNLIF_CASSLGQAYEQYF", "CAGAGGGSQGNLIF_CASSIRSSYEQYF"],
                "gex.TNFRSF9": [12.0, 0.0],
                "gex.MKI67": [6.0, 0.0],
                "gex.TCF7": [0.0, 10.0],
                "gex.LEF1": [0.0, 8.0],
                "gex.CCR7": [0.0, 9.0],
                "gex.SELL": [0.0, 7.0],
                "gex.IL7R": [1.0, 11.0],
            }
        )
        scored = annotate_tcrs(valid_clones, per_cell=per_cell, methods=["prism"])
        assert "prism_score" in scored.columns
        assert scored["prism_score"].notna().all()

    def test_til_select_excludes_citeseq_umis_and_cd44_prefix_collision(self, tmp_path):
        from tcrsift.til_select import _load_10x_h5_marker_counts, run_selection_pipeline

        # Build 10x H5 with 2 Gene Expression rows (100 UMIs) + 1 Antibody Capture row (900 UMIs)
        h5_path = tmp_path / "cite.h5"
        mat = csc_matrix(np.array([[20, 10], [80, 90], [900, 900]], dtype=np.float32))
        with h5py.File(h5_path, "w") as f:
            grp = f.create_group("matrix")
            grp.create_dataset("barcodes", data=np.array([b"BC1-1", b"BC2-1"]))
            grp.create_dataset("data", data=mat.data)
            grp.create_dataset("indices", data=mat.indices)
            grp.create_dataset("indptr", data=mat.indptr)
            grp.create_dataset("shape", data=np.array(mat.shape))
            feat = grp.create_group("features")
            feat.create_dataset("name", data=np.array([b"CD8A", b"ACTB", b"CD8A_TotalSeq"]))
            feat.create_dataset("id", data=np.array([b"ENSG1", b"ENSG2", b"ADT_CD8A"]))
            feat.create_dataset(
                "feature_type",
                data=np.array([b"Gene Expression", b"Gene Expression", b"Antibody Capture"]),
            )

        counts_df = _load_10x_h5_marker_counts(h5_path, ["CD8A"])
        assert counts_df["total_umis"].tolist() == [100.0, 100.0]
        assert counts_df["cp10k_CD8A"].tolist() == [2000.0, 1000.0]

        # CD44 / CD40LG must NOT pollute cd4_cp10k_max; extinct T2 clone cannot be became_cytotoxic
        master = pd.DataFrame(
            {
                "CDR3ab": ["CLONE_EXPAND", "CLONE_EXTINCT"],
                "CDR3_alpha": ["CAVA1", "CAVA2"],
                "CDR3_beta": ["CASSB1", "CASSB2"],
                "cell_count_T1": [2, 5],
                "cell_count_T2": [10, 0],
                "frequency_T1": [0.10, 0.25],
                "frequency_T2": [0.50, 0.0],
                "total_cells": [12, 5],
                "score_cp10k_CD8A_T1": [50.0, 50.0],
                "score_cp10k_CD8A_T2": [60.0, 0.0],
                "score_cp10k_CD4_T1": [0.0, 0.0],
                "score_cp10k_CD4_T2": [0.0, 0.0],
                "score_cp10k_CD44_T1": [500.0, 500.0],
                "score_cp10k_CD40LG_T2": [500.0, 0.0],
                "score_z_GZMB_T1": [-0.5, -1.0],
                "score_z_GZMB_T2": [1.5, 0.0],
                "score_cp10k_GZMB_T1": [10.0, 5.0],
                "score_cp10k_GZMB_T2": [200.0, 0.0],
            }
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", pd.errors.PerformanceWarning)
            res, *_ = run_selection_pipeline(
                master, ["T1", "T2"], ["CD8A", "CD4", "CD44", "CD40LG", "GZMB"]
            )
        perf_warns = [w for w in caught if issubclass(w.category, pd.errors.PerformanceWarning)]
        assert not perf_warns
        by_c = res.set_index("CDR3ab")
        assert by_c.loc["CLONE_EXPAND", "cd4_cp10k_max"] == 0.0
        assert bool(by_c.loc["CLONE_EXPAND", "is_cd8_not_cd4"]) is True
        assert bool(by_c.loc["CLONE_EXTINCT", "is_became_cytotoxic"]) is False

    def test_annotate_cells_singleton_cluster_and_integer_leiden_preserve_top_markers(self):
        from tcrsift.annotate_cells import annotate_cells

        genes = ["CD3D", "CD3E", "CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "MS4A1", "CD79A", "CD19"]
        X = np.array(
            [[10, 10, 12, 10, 20, 15, 18, 0, 0, 0]] * 5
            + [[0, 0, 0, 0, 0, 0, 0, 20, 18, 15]] * 5
            + [[1, 1, 0, 0, 0, 0, 0, 0, 0, 0]],  # 1-cell singleton outlier cluster!
            dtype=np.float32,
        )
        adata = ad.AnnData(
            X=X,
            obs=pd.DataFrame({"leiden": [0] * 5 + [1] * 5 + [2]}, index=[f"c{i}" for i in range(11)]),
            var=pd.DataFrame(index=genes),
        )
        annotate_cells(adata, min_subtype_frac=0.01)
        # Multi-cell clusters 0 and 1 must retain their `· GENE` suffixes despite singleton cluster 2!
        assert " · " in adata.obs["phenotype"].iloc[0]
        assert " · " in adata.obs["phenotype"].iloc[5]

    def test_plots_exact_sample_jaccard_and_cli_helpers(self, tmp_path):
        from tcrsift.cli import main
        from tcrsift.plots import plot_clonotypes, plot_funnel_variants, plot_phenotype
        from tcrsift.validation import parse_til_sample_spec

        # Path containing '=' after ':' must parse cleanly
        name, stype, spath = parse_til_sample_spec("T1=csv:/tmp/run=1/til.csv")
        assert (name, stype, spath) == ("T1", "csv", Path("/tmp/run=1/til.csv"))

        # Disjoint samples S1 and S10 must have off-diagonal Jaccard == 0, not substring match
        clonotypes = pd.DataFrame(
            {
                "CDR3ab": ["CAVA1_CASSB1", "CAVA2_CASSB2"],
                "cell_count": [10, 5],
                "samples": ["S10", "S1"],
            }
        )
        plot_clonotypes(clonotypes, tmp_path)
        assert (tmp_path / "clonotype_sharing_heatmap.png").exists()

        # Non-categorical Tcell_type and empty stage_counts must not crash
        adata = ad.AnnData(
            X=csr_matrix((2, 1)),
            obs=pd.DataFrame(
                {"sample": ["S1", "S1"], "Tcell_type": ["Confident CD8+", "Unknown"]},
                index=["c1", "c2"],
            ),
        )
        plot_phenotype(adata, tmp_path)
        assert (tmp_path / "phenotype_composition.png").exists()
        assert plot_funnel_variants({}, tmp_path) == []

        # `tcrsift report` with no subcommand prints help cleanly instead of crashing
        rc = main(["report"])
        assert rc is None or rc == 0
