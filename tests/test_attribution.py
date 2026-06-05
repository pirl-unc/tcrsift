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

"""Tests for partial-information clonotype attribution (#176)."""

from __future__ import annotations

import pandas as pd
import pytest

from tcrsift.attribution import attribute_cells
from tcrsift.config import AttributionConfig


def _cells(rows):
    """Build a cell DataFrame; fills the chain columns attribute_cells reads.

    Each row dict may set a1/a2/b1/b2 (CDR3 strings, "" = absent) and sample;
    UMIs default to 5 (passing) for any present chain.
    """
    cols = ["CDR3_alpha", "TRA_2_cdr3", "CDR3_beta", "TRB_2_cdr3",
            "TRA_1_umis", "TRA_2_umis", "TRB_1_umis", "TRB_2_umis", "sample"]
    data = {c: [] for c in cols}
    idx = []
    for i, r in enumerate(rows):
        idx.append(r.get("bc", f"c{i}"))
        a1, a2 = r.get("a1", ""), r.get("a2", "")
        b1, b2 = r.get("b1", ""), r.get("b2", "")
        data["CDR3_alpha"].append(a1)
        data["TRA_2_cdr3"].append(a2)
        data["CDR3_beta"].append(b1)
        data["TRB_2_cdr3"].append(b2)
        data["TRA_1_umis"].append(r.get("au1", 5 if a1 else 0))
        data["TRA_2_umis"].append(r.get("au2", 5 if a2 else 0))
        data["TRB_1_umis"].append(r.get("bu1", 5 if b1 else 0))
        data["TRB_2_umis"].append(r.get("bu2", 5 if b2 else 0))
        data["sample"].append(r.get("sample", "S1"))
    return pd.DataFrame(data, index=idx)


def _on(**kw):
    kw.setdefault("enabled", True)
    return AttributionConfig(**kw)


class TestCompleteCells:
    def test_complete_single_weight_one(self):
        df = _cells([{"a1": "A", "b1": "B"}, {"a1": "A", "b1": "B"}])
        long, merge = attribute_cells(df, _on())
        assert merge == {}
        assert set(long["CDR3ab"]) == {"A_B"}
        assert long["weight"].sum() == pytest.approx(2.0)
        # Each cell contributes exactly 1.0.
        per_cell = long.groupby("cell_barcode")["weight"].sum()
        assert per_cell.to_dict() == pytest.approx({"c0": 1.0, "c1": 1.0})


class TestBetaShare:
    def test_alpha_dropout_splits_proportional_3_to_1(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},   # A1_B prior = 3
            {"a1": "A2", "b1": "B"},   # A2_B prior = 1
            {"bc": "drop", "a1": "", "b1": "B"},  # alpha-dropout, beta B
        ])
        long, _ = attribute_cells(df, _on())
        drop = long[long["cell_barcode"] == "drop"].set_index("CDR3ab")["weight"]
        assert drop["A1_B"] == pytest.approx(0.75)
        assert drop["A2_B"] == pytest.approx(0.25)
        assert drop.sum() == pytest.approx(1.0)

    def test_alpha_dropout_no_match_dropped(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"bc": "drop", "a1": "", "b1": "OTHER"},  # beta with no complete clone
        ])
        long, _ = attribute_cells(df, _on())
        assert "drop" not in set(long["cell_barcode"])

    def test_beta_sharing_disabled_drops(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"bc": "drop", "a1": "", "b1": "B"},
        ])
        long, _ = attribute_cells(df, _on(beta_sharing=False))
        assert "drop" not in set(long["cell_barcode"])

    def test_uniform_weighting_ignores_priors(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A2", "b1": "B"},
            {"bc": "drop", "a1": "", "b1": "B"},
        ])
        long, _ = attribute_cells(df, _on(share_weighting="uniform"))
        drop = long[long["cell_barcode"] == "drop"].set_index("CDR3ab")["weight"]
        assert drop["A1_B"] == pytest.approx(0.5)
        assert drop["A2_B"] == pytest.approx(0.5)


class TestDualAlpha:
    def test_recurrent_dual_alpha_merges(self):
        # Genuine allelic inclusion: the dual is the beta's DOMINANT config
        # (3 of 4 paired cells carry both alphas), so it merges (#198 gate passes).
        df = _cells([
            {"a1": "A1", "b1": "B"},   # one allele-dropout single (A1_B prior 1 -> canonical)
            {"bc": "d1", "a1": "A1", "a2": "A2", "b1": "B"},
            {"bc": "d2", "a1": "A2", "a2": "A1", "b1": "B"},
            {"bc": "d3", "a1": "A1", "a2": "A2", "b1": "B"},
        ])
        long, merge = attribute_cells(df, _on(dual_alpha_min_cells=2))
        assert merge.get("A2_B") == "A1_B"
        assert merge.get("A1_B") == "A1_B"
        # Dual cells attributed wholly to the merged clone.
        for bc in ("d1", "d2", "d3"):
            rows = long[long["cell_barcode"] == bc]
            assert list(rows["CDR3ab"]) == ["A1_B"]
            assert rows["weight"].iloc[0] == pytest.approx(1.0)
            assert rows["kind"].iloc[0] == "dual_alpha_merge"
        # A2_B no longer appears as its own clone (merged away).
        assert "A2_B" not in set(long["CDR3ab"])


class TestDualAlphaDominanceGate:
    """#198: only merge when the dual is the beta's dominant configuration, not
    a public/promiscuous beta's chance recurrence."""

    def _promiscuous(self):
        # Beta B pairs with 5 distinct alphas; the (A1,A2,B) dual recurs twice
        # but is a small minority of B's cells -> must NOT merge.
        rows = [{"a1": f"X{k}", "b1": "B"} for k in range(5)]
        rows += [{"a1": f"X{k}", "b1": "B"} for k in range(5)]   # more singles
        rows += [
            {"bc": "d1", "a1": "X0", "a2": "X1", "b1": "B"},
            {"bc": "d2", "a1": "X1", "a2": "X0", "b1": "B"},
        ]
        return _cells(rows)

    def test_promiscuous_beta_not_merged(self):
        long, merge = attribute_cells(self._promiscuous(), _on(dual_alpha_min_cells=2))
        # No merge: the dual is 2 of 12 B-cells (~17% < 0.5 majority).
        assert merge == {}

    def test_min_fraction_zero_restores_recurrence_only(self):
        # Opting out of the dominance gate reverts to bare recurrence (pre-#198).
        long, merge = attribute_cells(
            self._promiscuous(), _on(dual_alpha_min_cells=2, dual_alpha_min_fraction=0)
        )
        assert merge != {}

    def test_max_beta_partners_gate(self):
        # Even a dominant dual is rejected when the beta is too promiscuous.
        df = _cells([
            {"bc": "d1", "a1": "A1", "a2": "A2", "b1": "B"},
            {"bc": "d2", "a1": "A2", "a2": "A1", "b1": "B"},
            {"a1": "A3", "b1": "B"},  # a 3rd partner -> 3 distinct alphas on B
        ])
        # Dominant (2/3 cells) but 3 partners > max 2 -> no merge.
        _long, merge = attribute_cells(
            df, _on(dual_alpha_min_cells=2, dual_alpha_max_beta_partners=2)
        )
        assert merge == {}

    def test_singleton_dual_alpha_splits(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},   # A1_B prior 3
            {"a1": "A2", "b1": "B"},   # A2_B prior 1
            {"bc": "d1", "a1": "A1", "a2": "A2", "b1": "B"},  # singleton -> split
        ])
        long, merge = attribute_cells(df, _on(dual_alpha_min_cells=2))
        assert merge == {}
        rows = long[long["cell_barcode"] == "d1"].set_index("CDR3ab")["weight"]
        assert rows["A1_B"] == pytest.approx(0.75)
        assert rows["A2_B"] == pytest.approx(0.25)

    def test_merge_disabled_splits(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A2", "b1": "B"},
            {"bc": "d1", "a1": "A1", "a2": "A2", "b1": "B"},
            {"bc": "d2", "a1": "A1", "a2": "A2", "b1": "B"},
        ])
        long, merge = attribute_cells(df, _on(dual_alpha_merge=False))
        assert merge == {}
        assert (long[long["cell_barcode"] == "d1"]["kind"] == "dual_alpha_split").all()

    def test_dual_alpha_split_disabled_collapses_to_primary(self):
        # #181: dual_alpha_split=False collapses a non-merged dual-alpha to its
        # primary pair (full weight), instead of splitting it as a doublet.
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A2", "b1": "B"},
            {"bc": "d1", "a1": "A1", "a2": "A2", "b1": "B"},  # singleton dual-alpha
        ])
        long, _ = attribute_cells(
            df, _on(dual_alpha_min_cells=2, dual_alpha_split=False)
        )
        rows = long[long["cell_barcode"] == "d1"]
        assert len(rows) == 1
        assert rows["CDR3ab"].iloc[0] == "A1_B"
        assert rows["weight"].iloc[0] == pytest.approx(1.0)
        assert rows["kind"].iloc[0] == "primary"

    def test_split_flags_are_orthogonal(self):
        # #181: dual_alpha_split=False with dual_beta_split=True -> dual-alpha
        # collapses to primary while dual-beta still splits.
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A2", "b1": "B"},
            {"a1": "A", "b1": "B1"},
            {"a1": "A", "b1": "B1"},
            {"a1": "A", "b1": "B2"},
            {"bc": "da", "a1": "A1", "a2": "A2", "b1": "B"},   # dual-alpha singleton
            {"bc": "db", "a1": "A", "b1": "B1", "b2": "B2"},   # dual-beta
        ])
        long, _ = attribute_cells(
            df, _on(dual_alpha_min_cells=2, dual_alpha_split=False, dual_beta_split=True)
        )
        da = long[long["cell_barcode"] == "da"]
        assert len(da) == 1 and da["kind"].iloc[0] == "primary"
        db = long[long["cell_barcode"] == "db"]
        assert set(db["CDR3ab"]) == {"A_B1", "A_B2"}


class TestDualBeta:
    def test_no_phantom_zero_weight_clone(self):
        # The B2 partner has no complete prior, so a proportional split sends it
        # 0 weight — it must not surface as a phantom clone.
        df = _cells([
            {"a1": "A", "b1": "B1"},
            {"bc": "d1", "a1": "A", "b1": "B1", "b2": "B2"},  # B2 prior = 0
        ])
        long, _ = attribute_cells(df, _on())
        assert "A_B2" not in set(long["CDR3ab"])
        assert (long["weight"] > 0).all()
        # The cell's full weight still lands on the real clone.
        assert long[long["cell_barcode"] == "d1"]["weight"].sum() == pytest.approx(1.0)

    def test_dual_beta_splits(self):
        df = _cells([
            {"a1": "A", "b1": "B1"},
            {"a1": "A", "b1": "B1"},
            {"a1": "A", "b1": "B1"},   # A_B1 prior 3
            {"a1": "A", "b1": "B2"},   # A_B2 prior 1
            {"bc": "d1", "a1": "A", "b1": "B1", "b2": "B2"},
        ])
        long, _ = attribute_cells(df, _on())
        rows = long[long["cell_barcode"] == "d1"].set_index("CDR3ab")["weight"]
        assert rows["A_B1"] == pytest.approx(0.75)
        assert rows["A_B2"] == pytest.approx(0.25)


class TestWeightConservation:
    def test_every_cell_sums_to_one_or_zero(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A2", "b1": "B"},
            {"bc": "drop_ok", "a1": "", "b1": "B"},       # -> 1.0
            {"bc": "drop_no", "a1": "", "b1": "ZZZ"},     # -> 0.0 (no match)
            {"bc": "dual", "a1": "A1", "a2": "A2", "b1": "B"},  # -> 1.0
            {"bc": "nochain", "a1": "", "b1": ""},        # -> 0.0
            {"bc": "alpha_only", "a1": "A1", "b1": ""},   # -> 0.0 (no beta)
        ])
        long, _ = attribute_cells(df, _on())
        sums = long.groupby("cell_barcode")["weight"].sum()
        for bc in ("c0", "c1", "drop_ok", "dual"):
            assert sums.get(bc, 0.0) == pytest.approx(1.0), bc
        for bc in ("drop_no", "nochain", "alpha_only"):
            assert sums.get(bc, 0.0) == pytest.approx(0.0), bc

    def test_umi_gate_excludes_low_umi_chain(self):
        # Beta present but below UMI gate -> treated as no beta -> dropped.
        df = _cells([
            {"a1": "A", "b1": "B"},
            {"bc": "low", "a1": "A", "b1": "B", "bu1": 1},
        ])
        long, _ = attribute_cells(df, _on(), min_umi=2)
        assert "low" not in set(long["cell_barcode"])


def _adata(rows):
    """Wrap _cells rows in a minimal AnnData for aggregate_clonotypes."""
    import anndata as ad
    from scipy.sparse import csr_matrix

    obs = _cells(rows)
    X = csr_matrix((len(obs), 1))
    return ad.AnnData(X=X, obs=obs)


class TestAggregateIntegration:
    def test_off_path_byte_identical_and_integer(self):
        from tcrsift.clonotype import aggregate_clonotypes

        rows = [
            {"a1": "A", "b1": "B"},
            {"a1": "A", "b1": "B"},
            {"a1": "C", "b1": "D"},
        ]
        base = aggregate_clonotypes(_adata(rows), show_progress=False, verbose=False)
        off = aggregate_clonotypes(
            _adata(rows), attribution=_on(enabled=False), show_progress=False, verbose=False
        )
        # Disabled attribution == no attribution, byte-identical.
        pd.testing.assert_frame_equal(base, off)
        # OFF keeps integer cell counts.
        assert base["cell_count"].dtype.kind in "iu"
        assert base.set_index("CDR3ab")["cell_count"].to_dict() == {"A_B": 2, "C_D": 1}

    def test_on_path_weighted_counts_and_beta_share(self):
        from tcrsift.clonotype import aggregate_clonotypes

        rows = [
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},   # A1_B prior 3
            {"a1": "A2", "b1": "B"},   # A2_B prior 1
            {"a1": "", "b1": "B"},     # alpha-dropout -> 0.75/0.25
        ]
        clo = aggregate_clonotypes(
            _adata(rows), attribution=_on(), show_progress=False, verbose=False
        )
        counts = clo.set_index("CDR3ab")["cell_count"].to_dict()
        # A1_B: 3 + 0.75 = 3.75 ; A2_B: 1 + 0.25 = 1.25
        assert counts["A1_B"] == pytest.approx(3.75)
        assert counts["A2_B"] == pytest.approx(1.25)
        # Total weight conserved (5 cells, all attributable).
        assert clo["cell_count"].sum() == pytest.approx(5.0)

    def test_per_sample_frequency_sum_is_one(self):
        from tcrsift.clonotype import build_clone_sample_long

        rows = [
            {"a1": "A1", "b1": "B", "sample": "S1"},
            {"a1": "A1", "b1": "B", "sample": "S1"},
            {"a1": "A2", "b1": "B", "sample": "S1"},
            {"a1": "", "b1": "B", "sample": "S1"},   # alpha-dropout, shared
            {"a1": "A1", "a2": "A2", "b1": "B", "sample": "S1"},  # dual-alpha
        ]
        long = build_clone_sample_long(_adata(rows), attribution=_on())
        freq_sum = long[long["sample"] == "S1"]["frequency"].sum()
        assert freq_sum == pytest.approx(1.0)

    def test_handle_doublets_remove_superseded(self):
        from tcrsift.clonotype import aggregate_clonotypes

        rows = [
            {"a1": "A1", "b1": "B"},
            {"a1": "A2", "b1": "B"},
            # A multi-chain doublet cell that handle_doublets='remove' would drop.
            {"bc": "dbl", "a1": "A1", "a2": "A2", "b1": "B"},
        ]
        adata = _adata(rows)
        adata.obs["multi_chain"] = [False, False, True]
        clo = aggregate_clonotypes(
            adata, handle_doublets="remove", attribution=_on(dual_alpha_merge=False),
            show_progress=False, verbose=False,
        )
        # The doublet cell's weight is redistributed, not dropped: total = 3.0.
        assert clo["cell_count"].sum() == pytest.approx(3.0)


class TestRepresentativeVdjConsistency:
    """#184: the representative VDJ for a chain must come from a cell+slot whose
    CDR3 equals the clone's canonical CDR3 — primary slot normally, secondary
    slot for a merged dual-alpha whose canonical alpha is a cell's TRA_2."""

    def test_picks_primary_slot_matching_cdr3(self):
        from tcrsift.clonotype import _pick_chain_representative

        df = pd.DataFrame({
            "TRA_1_cdr3": ["CAAA", "CAAA", "CDDD"],
            "TRA_2_cdr3": ["", "CDDD", ""],
            "TRA_1_umis": [5, 5, 99],  # the CDDD cell has the most UMIs
            "TRA_2_umis": [0, 3, 0],
            "TRA_1_vdj_aa": ["xxCAAAxx", "xxCAAAxx", "xxCDDDxx"],
            "TRA_2_vdj_aa": ["", "xxCDDDxx", ""],
        })
        prefix, rep = _pick_chain_representative(df, "TRA", "CAAA")
        # Despite the CDDD cell having more UMIs, the rep carries CAAA.
        assert prefix == "TRA_1"
        assert "CAAA" in rep["TRA_1_vdj_aa"] and "CDDD" not in rep["TRA_1_vdj_aa"]

    def test_falls_to_secondary_slot_for_merged_alpha(self):
        from tcrsift.clonotype import _pick_chain_representative

        # Canonical alpha CEEE only appears as the cell's TRA_2 (allelic
        # inclusion); the rep + VDJ must come from the secondary slot.
        df = pd.DataFrame({
            "TRA_1_cdr3": ["CAAA"],
            "TRA_2_cdr3": ["CEEE"],
            "TRA_1_umis": [5],
            "TRA_2_umis": [3],
            "TRA_1_vdj_aa": ["xxCAAAxx"],
            "TRA_2_vdj_aa": ["xxCEEExx"],
        })
        prefix, rep = _pick_chain_representative(df, "TRA", "CEEE")
        assert prefix == "TRA_2"
        assert "CEEE" in rep["TRA_2_vdj_aa"]

    def test_fallback_when_no_cdr3_columns(self):
        from tcrsift.clonotype import _pick_chain_representative

        df = pd.DataFrame({"TRA_1_umis": [5, 9]})  # no _cdr3 columns
        prefix, rep = _pick_chain_representative(df, "TRA", "CAAA")
        assert prefix == "TRA_1" and rep is not None

    def test_merged_clone_vdj_contains_canonical_cdr3(self):
        import anndata as ad
        from scipy.sparse import csr_matrix

        from tcrsift.clonotype import aggregate_clonotypes

        # a1=CAAA (prior 3, canonical) and a2=CDDD (prior 1) share beta CBBB;
        # two recurrent dual-alpha cells merge them. The CDDD cell has the most
        # alpha UMIs, so pre-#184 the merged clone's VDJ_alpha would be CDDD's.
        def cell(a1, a2, au1, au2=0):
            return {
                "CDR3_alpha": a1, "TRA_2_cdr3": a2, "CDR3_beta": "CBBB", "TRB_2_cdr3": "",
                "TRA_1_cdr3": a1, "TRA_1_umis": au1, "TRA_2_umis": au2,
                "TRB_1_cdr3": "CBBB", "TRB_1_umis": 5,
                "TRA_1_vdj_aa": f"GG{a1}WW" if a1 else "",
                "TRA_2_vdj_aa": f"GG{a2}WW" if a2 else "",
                "TRB_1_vdj_aa": "GGCBBBWW",
                "sample": "S1",
            }
        rows = [
            cell("CAAA", "", 5),
            cell("CAAA", "", 5),
            cell("CAAA", "", 5),
            cell("CDDD", "", 99),               # most UMIs, wrong alpha
            cell("CAAA", "CDDD", 5, au2=4),     # dual-alpha
            cell("CAAA", "CDDD", 5, au2=4),     # dual-alpha (recurs -> merge)
        ]
        obs = pd.DataFrame(rows, index=[f"c{i}" for i in range(len(rows))])
        adata = ad.AnnData(X=csr_matrix((len(obs), 1)), obs=obs)
        clo = aggregate_clonotypes(
            adata, attribution=_on(dual_alpha_min_cells=2),
            show_progress=False, verbose=False,
        )
        merged = clo[clo["CDR3ab"] == "CAAA_CBBB"].iloc[0]
        assert "CAAA" in merged["VDJ_alpha_aa"]
        assert "CDDD" not in merged["VDJ_alpha_aa"]


class TestDisabledAndGuards:
    def test_cdr3b_only_returns_empty(self):
        df = _cells([{"a1": "A", "b1": "B"}])
        long, merge = attribute_cells(df, _on(), group_by="CDR3b_only")
        assert long.empty and merge == {}

    def test_empty_input(self):
        df = _cells([])
        long, merge = attribute_cells(df, _on())
        assert long.empty and merge == {}
