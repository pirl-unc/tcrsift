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

"""Tests for #118 (allele-NaN reason column + uniform qc_warnings),
#119 (cohort-level novel-allele detection), and #120 (per-clone
observed AA + divergence positions + polymorphic-position
tolerance).
"""

from __future__ import annotations

import pandas as pd

from tcrsift.assemble import (
    ALLELE_REASON_AUTO_DETECTED,
    ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION,
    ALLELE_REASON_INVALID_OVERRIDE,
    ALLELE_REASON_NO_CONTIG,
    ALLELE_REASON_OVERRIDDEN,
    ALLELE_REASON_SPARSE_CONTIG,
    ALLELE_REASON_VALUES,
    HUMAN_TRAC_AA,
    HUMAN_TRBC1_AA,
    _divergence_positions,
    _polymorphic_positions,
    allele_audit_report,
    assemble_full_sequences,
    back_translate,
    detect_novel_alleles,
)

# ---------------------------------------------------------------------------
# Pure-function helpers
# ---------------------------------------------------------------------------


class TestPolymorphicPositions:
    def test_empty_pool(self):
        assert _polymorphic_positions({}) == set()

    def test_single_allele(self):
        assert _polymorphic_positions({"01": "ABCDE"}) == set()

    def test_two_alleles_differ_at_pos_1(self):
        pool = {"01": "ABCDE", "02": "AXCDE"}
        assert _polymorphic_positions(pool) == {1}

    def test_truncates_at_shortest_allele(self):
        pool = {"01": "ABCDE", "02": "ABCD"}
        # Length 4; pos 0..3 all agree → no polymorphic positions.
        assert _polymorphic_positions(pool) == set()


class TestDivergencePositions:
    def test_no_divergence(self):
        assert _divergence_positions("ABCDE", "ABCDE") == []

    def test_single_divergence(self):
        # Pos 1-indexed: pos 3 in "ABCDE" vs "ABXDE" differs.
        assert _divergence_positions("ABCDE", "ABXDE") == [(3, "C", "X")]

    def test_caps_at_max_positions(self):
        a = "A" * 20
        b = "B" * 20
        # Default cap is 15.
        divs = _divergence_positions(a, b)
        assert len(divs) == 15
        assert all(d[0] <= 15 for d in divs)


# ---------------------------------------------------------------------------
# End-to-end via assemble_full_sequences
# ---------------------------------------------------------------------------


def _vdj_fixture():
    vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
    vdj_beta = "CASS" + "G" * 60 + "VETA"
    return vdj_alpha, vdj_beta


def _base_clone(vdj_alpha, vdj_beta, **extra):
    row = {
        "CDR3ab": "c1",
        "CDR3_alpha": vdj_alpha,
        "CDR3_beta": vdj_beta,
        "VDJ_alpha_aa": vdj_alpha,
        "VDJ_beta_aa": vdj_beta,
        "VDJ_alpha_nt": back_translate(vdj_alpha),
        "VDJ_beta_nt": back_translate(vdj_beta),
        "alpha_c_gene": "TRAC",
        "beta_c_gene": "TRBC1",
        "beta_j_gene": "TRBJ1-1",
        "samples": "S1",
    }
    row.update(extra)
    return row


class TestAlleleReasonNoContig:
    """Without contigs the picker has nothing to score; reason
    should be ``no_contig`` (#118)."""

    def test_no_contigs_dir(self):
        a, b = _vdj_fixture()
        df = pd.DataFrame([_base_clone(a, b)])
        out = assemble_full_sequences(
            df, alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        for chain in ("alpha", "beta"):
            assert out[f"{chain}_allele_called"].iloc[0] is None
            assert out[f"{chain}_allele_called_reason"].iloc[0] == ALLELE_REASON_NO_CONTIG


class TestAlleleReasonExplicitOverride:
    def test_valid_override(self):
        a, b = _vdj_fixture()
        df = pd.DataFrame([_base_clone(a, b, beta_c_gene="TRBC2", beta_j_gene="TRBJ2-1")])
        out = assemble_full_sequences(
            df, alpha_leader=None, beta_leader=None,
            trbc2_allele="03",
            verbose=False, show_progress=False,
        )
        assert out["beta_allele_called"].iloc[0] == "TRBC2*03"
        assert out["beta_allele_called_reason"].iloc[0] == ALLELE_REASON_OVERRIDDEN

    def test_invalid_override_emits_qc_warning(self):
        a, b = _vdj_fixture()
        df = pd.DataFrame([_base_clone(a, b, beta_c_gene="TRBC2", beta_j_gene="TRBJ2-1")])
        out = assemble_full_sequences(
            df, alpha_leader=None, beta_leader=None,
            trbc2_allele="99-not-a-real-allele",
            verbose=False, show_progress=False,
        )
        assert out["beta_allele_called"].iloc[0] is None
        assert out["beta_allele_called_reason"].iloc[0] == ALLELE_REASON_INVALID_OVERRIDE
        # qc_warnings populated.
        warnings = out["qc_warnings"].iloc[0]
        assert any("invalid allele override" in w for w in warnings)


class TestAlleleReasonSparseContig:
    """When the contig covers fewer than ``_DEFAULT_MIN_PICKER_POSITIONS``
    codons, the reason is ``sparse_contig`` and qc_warnings is populated
    (the gap #118 reports)."""

    def test_sparse_contig_populates_warning(self, tmp_path):
        a, b = _vdj_fixture()
        # Contig with only 4 codons past VDJ — below the 10-codon floor.
        leader_aa = "M" + "A" * 19
        sparse_beta_c_nt = back_translate(HUMAN_TRBC1_AA[:4])  # 12 nt
        beta_contig = back_translate(leader_aa) + back_translate(b) + sparse_beta_c_nt

        df = pd.DataFrame([_base_clone(
            a, b,
            beta_contig_ids="contig_b1",
        )])
        contig_dir = tmp_path / "S1"
        contig_dir.mkdir(parents=True)
        (contig_dir / "filtered_contig.fasta").write_text(
            f">contig_b1\n{beta_contig}\n"
        )

        out = assemble_full_sequences(
            df, contigs_dir=str(tmp_path),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        assert out["beta_allele_called"].iloc[0] is None
        assert out["beta_allele_called_reason"].iloc[0] == ALLELE_REASON_SPARSE_CONTIG
        warnings = out["qc_warnings"].iloc[0]
        assert any(
            "contig covers only" in w and "codons" in w
            for w in warnings
        ), f"sparse warning missing, got: {warnings}"


class TestAlleleReasonAutoDetected:
    def test_clean_contig_match(self, tmp_path):
        a, b = _vdj_fixture()
        leader_aa = "M" + "A" * 19
        # Real CellRanger contigs put the J→C junction codon (E for
        # β, GAG) at the very start of the bytes past VDJ; the
        # canonical AA starts AFTER that codon. Mirror that here so
        # the picker (which translates contig past junction) sees
        # HUMAN_TRBC1_AA[:15] in frame against the same allele.
        beta_c_nt = "GAG" + back_translate(HUMAN_TRBC1_AA[:15])
        beta_contig = back_translate(leader_aa) + back_translate(b) + beta_c_nt

        df = pd.DataFrame([_base_clone(a, b, beta_contig_ids="contig_b1")])
        contig_dir = tmp_path / "S1"
        contig_dir.mkdir(parents=True)
        (contig_dir / "filtered_contig.fasta").write_text(
            f">contig_b1\n{beta_contig}\n"
        )
        out = assemble_full_sequences(
            df, contigs_dir=str(tmp_path),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        assert out["beta_allele_called_reason"].iloc[0] == ALLELE_REASON_AUTO_DETECTED
        assert out["beta_allele_called"].iloc[0] is not None


class TestAlleleReasonValuesEnum:
    def test_all_values_are_strings(self):
        for v in ALLELE_REASON_VALUES:
            assert isinstance(v, str)
            assert v.strip() == v


# ---------------------------------------------------------------------------
# #120 — observed AA + divergence positions
# ---------------------------------------------------------------------------


class TestObservedConstantAaStart:
    def test_no_contig_means_none(self):
        a, b = _vdj_fixture()
        df = pd.DataFrame([_base_clone(a, b)])
        out = assemble_full_sequences(
            df, alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        assert out["beta_observed_constant_aa_start"].iloc[0] is None

    def test_captures_contig_translation(self, tmp_path):
        a, b = _vdj_fixture()
        leader_aa = "M" + "A" * 19
        # Junction codon (E for β) + canonical mature AA — mirrors
        # real CellRanger contig layout.
        beta_c_nt = "GAG" + back_translate(HUMAN_TRBC1_AA[:15])
        beta_contig = back_translate(leader_aa) + back_translate(b) + beta_c_nt
        df = pd.DataFrame([_base_clone(a, b, beta_contig_ids="contig_b1")])
        contig_dir = tmp_path / "S1"
        contig_dir.mkdir(parents=True)
        (contig_dir / "filtered_contig.fasta").write_text(
            f">contig_b1\n{beta_contig}\n"
        )
        out = assemble_full_sequences(
            df, contigs_dir=str(tmp_path),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        obs = out["beta_observed_constant_aa_start"].iloc[0]
        # observed_constant_aa_start is the translation of contig past
        # the junction codon, capped at 15 residues — exactly
        # HUMAN_TRBC1_AA[:15] here.
        assert obs == HUMAN_TRBC1_AA[:15]


# ---------------------------------------------------------------------------
# #119 — cohort novel-allele detection
# ---------------------------------------------------------------------------


class TestDetectNovelAlleles:
    def _make_cohort(self, n_clones: int, n_v_genes: int,
                     n_samples: int, divergent_fraction: float):
        """Build a synthetic full_sequences-shaped DataFrame with
        ``n_clones`` chains, ``divergent_fraction`` of them carrying
        a polymorphism at TRAC mature position 3 (N→K)."""
        rows = []
        n_divergent = int(round(n_clones * divergent_fraction))
        for i in range(n_clones):
            is_divergent = i < n_divergent
            v_gene = f"TRAV{1 + (i % n_v_genes)}"
            sample = f"S{1 + (i % n_samples)}"
            # Build observed AA — replace canonical TRAC[2] with K
            # when divergent.
            canonical = HUMAN_TRAC_AA[:15]
            if is_divergent:
                # Note: observed_constant_aa_start is the contig past
                # the JUNCTION CODON. The junction codon for α can be
                # N (varies by J), but for this test we'll use the
                # canonical without junction; divergence at mature pos
                # 3 maps to observed position 3 (1-indexed: 3rd char).
                observed = canonical[:2] + "K" + canonical[3:]
            else:
                observed = canonical
            rows.append({
                "alpha_c_gene": "TRAC",
                "alpha_c_gene_canonical": "TRAC",
                "alpha_allele_called": None if is_divergent else "TRAC*01",
                "alpha_allele_called_reason": (
                    ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION
                    if is_divergent else ALLELE_REASON_AUTO_DETECTED
                ),
                "alpha_observed_constant_aa_start": observed,
                "alpha_v_gene": v_gene,
                "samples": sample,
            })
        return pd.DataFrame(rows)

    def test_high_freq_spread_classified_as_novel(self):
        # 50% of clones (heterozygous-like), 5 V-genes, 3 samples → novel.
        df = self._make_cohort(
            n_clones=200, n_v_genes=5, n_samples=3,
            divergent_fraction=0.50,
        )
        result = detect_novel_alleles(df, min_pct=0.05, min_v_spread=3, min_samples=2)
        assert not result.empty
        # Find the row matching the N→K substitution.
        nk_rows = result[
            (result["expected_aa"] == "N") & (result["observed_aa"] == "K")
        ]
        assert len(nk_rows) == 1
        assert nk_rows.iloc[0]["verdict"] == "novel_allele_candidate"

    def test_low_freq_classified_as_artifact(self):
        df = self._make_cohort(
            n_clones=200, n_v_genes=1, n_samples=1,
            divergent_fraction=0.02,  # only 2% — below threshold
        )
        result = detect_novel_alleles(df, min_pct=0.05, min_v_spread=3, min_samples=2)
        if result.empty:
            # No detected divergence; that's fine — 2% × 1 V-gene
            # may register as a single artifact row.
            return
        assert all(result["verdict"] == "likely_artifact")

    def test_v_concentrated_high_freq_is_artifact(self):
        # 50% divergence but all in ONE V-gene → artifact.
        df = self._make_cohort(
            n_clones=200, n_v_genes=1, n_samples=3,
            divergent_fraction=0.50,
        )
        result = detect_novel_alleles(df, min_pct=0.05, min_v_spread=3, min_samples=2)
        nk_rows = result[
            (result["expected_aa"] == "N") & (result["observed_aa"] == "K")
        ]
        assert len(nk_rows) == 1
        assert nk_rows.iloc[0]["verdict"] == "likely_artifact"

    def test_empty_df_returns_empty(self):
        df = pd.DataFrame({"alpha_allele_called_reason": []})
        result = detect_novel_alleles(df)
        assert result.empty


class TestAlleleAuditReport:
    def test_report_renders_basic_breakdown(self):
        a, b = _vdj_fixture()
        df = pd.DataFrame([
            _base_clone(a, b),  # no_contig
            _base_clone(a, b),  # no_contig
        ])
        out = assemble_full_sequences(
            df, alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        report = allele_audit_report(out)
        assert "[allele audit]" in report
        assert "alpha chain" in report
        assert ALLELE_REASON_NO_CONTIG in report
