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
    ALLELE_REASON_DIVERGENT_CONTIG,
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


class TestAlleleReasonDivergentContig:
    """E2E for ``divergent_contig`` reason — contig translates but
    disagrees at every comparable position vs all packaged alleles."""

    def test_all_zero_score_yields_divergent_reason(self, tmp_path):
        a, b = _vdj_fixture()
        leader_aa = "M" + "A" * 19
        # Contig long enough to clear the sparse floor (≥10 codons)
        # but with a totally different AA — translates to W's, which
        # don't match anywhere in HUMAN_TRBC1_AA's prefix.
        # Junction codon (E for β) then 12 W codons.
        beta_c_nt = "GAG" + ("TGG" * 12)  # E + WWWWWWWWWWWW
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
        assert out["beta_allele_called"].iloc[0] is None
        assert out["beta_allele_called_reason"].iloc[0] == ALLELE_REASON_DIVERGENT_CONTIG
        warnings = out["qc_warnings"].iloc[0]
        assert any("disagrees with every packaged allele" in w for w in warnings)


class TestAlleleReasonDivergentAtPolymorphicPosition:
    """E2E for the #120 ask-4 behavior: when the contig disagrees with
    the best-fit allele at an allele-distinguishing position (TRBC2
    mature pos 9 between *01 and *03), the picker punts to NaN
    instead of silently flipping to the sibling on overall score."""

    @staticmethod
    def _trbc2_contig_with_residue_at_pos9(residue: str) -> str:
        from tcrsift.assemble import HUMAN_CONSTANT_ALLELES

        # Build a TRBC2 contig that matches *01 everywhere EXCEPT
        # at mature pos 9, where we inject ``residue``.
        trbc2_01 = HUMAN_CONSTANT_ALLELES["TRBC2"]["01"]
        canonical_prefix = trbc2_01[:8] + residue + trbc2_01[9:15]
        # Junction codon (E for β) + 15 mature codons.
        return "GAG" + back_translate(canonical_prefix)

    def test_disagrees_with_both_alleles_at_pos9_punts(self, tmp_path):
        a, b = _vdj_fixture()
        leader_aa = "M" + "A" * 19
        # Inject Q at pos 9 — disagrees with both *01 (E) and *03 (K).
        beta_c_nt = self._trbc2_contig_with_residue_at_pos9("Q")
        beta_contig = back_translate(leader_aa) + back_translate(b) + beta_c_nt

        df = pd.DataFrame([_base_clone(
            a, b,
            beta_c_gene="TRBC2",
            beta_j_gene="TRBJ2-1",
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
        assert (
            out["beta_allele_called_reason"].iloc[0]
            == ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION
        )
        # qc_warnings populated and names the position.
        warnings = out["qc_warnings"].iloc[0]
        assert any(
            "allele-distinguishing position" in w and "9" in w
            for w in warnings
        ), f"polymorphic warning missing or unspecific, got: {warnings}"
        # The observed AA carries the donor's actual residue at pos 9
        # (1-indexed pos 9 = 0-indexed idx 8).
        observed = out["beta_observed_constant_aa_start"].iloc[0]
        assert observed[8] == "Q"

    def test_allele_alternatives_preserved_when_polymorphic_punts(self, tmp_path):
        # The picker keeps ``allele_alternatives`` populated even when
        # it refuses to commit — that lets cohort audits see which
        # alleles were close.
        a, b = _vdj_fixture()
        leader_aa = "M" + "A" * 19
        beta_c_nt = self._trbc2_contig_with_residue_at_pos9("Q")
        beta_contig = back_translate(leader_aa) + back_translate(b) + beta_c_nt
        df = pd.DataFrame([_base_clone(
            a, b,
            beta_c_gene="TRBC2",
            beta_j_gene="TRBJ2-1",
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
        alts = out["beta_allele_alternatives"].iloc[0]
        assert isinstance(alts, str)
        # Both packaged TRBC2 alleles should appear in the audit.
        assert "TRBC2*01" in alts
        assert "TRBC2*03" in alts

    def test_constant_bytes_use_default_when_punted(self, tmp_path):
        """The #120 docstring note: punt-to-NaN punts the CALL, not
        the bytes. The user's ``_constant_nt_optimized`` still uses
        the default allele's canonical residue at the polymorphic
        position. Lock that in so a future refactor that accidentally
        drops the bytes too gets caught."""
        a, b = _vdj_fixture()
        leader_aa = "M" + "A" * 19
        beta_c_nt = self._trbc2_contig_with_residue_at_pos9("Q")
        beta_contig = back_translate(leader_aa) + back_translate(b) + beta_c_nt
        df = pd.DataFrame([_base_clone(
            a, b,
            beta_c_gene="TRBC2",
            beta_j_gene="TRBJ2-1",
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
        # constant_aa was built from the DEFAULT allele (HUMAN_TRBC2_AA)
        # which has E at mature pos 9 (idx 8). With the J-junction E
        # prepended, that maps to constant_aa[9] (idx 9 = junction +
        # idx 8 of mature).
        from tcrsift.assemble import HUMAN_TRBC2_AA
        constant_aa = out["beta_constant_aa"].iloc[0]
        # Junction E + mature TRBC2 → constant_aa[9] == HUMAN_TRBC2_AA[8].
        assert constant_aa[9] == HUMAN_TRBC2_AA[8]


class TestAlleleReasonFrameError:
    """When the contig's first codon past the junction doesn't
    translate to the picked allele's expected first residue, the
    picker punts with ``frame_error``."""

    def test_misaligned_first_codon_yields_frame_error(self, tmp_path):
        # Fixture rationale: contig past junction is (Q,
        # canonical[1:14]). With min_picker_positions=10 and
        # 13/14 AA matches, _pick_best_allele returns TRBC1*01 with
        # score 0.929 — definitely commits to a best_allele.
        # _verify_first_codon_frame then translates "CAG" → Q and
        # compares to HUMAN_TRBC1_AA[0] = E → returns False →
        # ALLELE_REASON_FRAME_ERROR fires deterministically.
        a, b = _vdj_fixture()
        leader_aa = "M" + "A" * 19
        from tcrsift.assemble import HUMAN_TRBC1_AA
        # Junction codon E + (Q at first position) + canonical[1:14].
        beta_c_nt = "GAG" + "CAG" + back_translate(HUMAN_TRBC1_AA[1:14])
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
        assert out["beta_allele_called"].iloc[0] is None
        assert out["beta_allele_called_reason"].iloc[0] == "frame_error"
        warnings = out["qc_warnings"].iloc[0]
        assert any(
            "frame error" in w.lower() or "frame error upstream" in w
            for w in warnings
        ), f"frame_error warning missing, got: {warnings}"


class TestAutoDetectedDivergencesFeedNovelAlleles:
    """#120 ask 3 plus #119: when the picker calls an allele but the
    contig still disagrees at non-distinguishing positions, those
    stored divergences should feed ``detect_novel_alleles``. This is
    the heterozygous TRAC N→K case from the issue — TRAC has one
    packaged allele, so the picker commits to TRAC*01 even when the
    contig has K instead of N; the cohort aggregator should still
    surface the polymorphism."""

    def test_stored_divergences_aggregate(self):
        # Build a 50-clone cohort with 25% K-at-position-3 across
        # 5 V-genes and 3 samples — should hit
        # novel_allele_candidate.
        rows = []
        for i in range(50):
            is_variant = i < 13  # 26%
            v_gene = f"TRAV{1 + (i % 5)}"
            sample = f"S{1 + (i % 3)}"
            if is_variant:
                divs_str = "3:N->K"
            else:
                divs_str = None
            rows.append({
                "alpha_c_gene": "TRAC",
                "alpha_c_gene_canonical": "TRAC",
                "alpha_allele_called": "TRAC*01",
                "alpha_allele_called_reason": ALLELE_REASON_AUTO_DETECTED,
                "alpha_allele_divergence_positions": divs_str,
                "alpha_observed_constant_aa_start": (
                    HUMAN_TRAC_AA[:2] + "K" + HUMAN_TRAC_AA[3:15]
                    if is_variant else HUMAN_TRAC_AA[:15]
                ),
                "alpha_v_gene": v_gene,
                "samples": sample,
            })
        df = pd.DataFrame(rows)
        result = detect_novel_alleles(
            df, min_pct=0.05, min_v_spread=3, min_samples=2,
        )
        nk_rows = result[
            (result["expected_aa"] == "N") & (result["observed_aa"] == "K")
        ]
        assert len(nk_rows) == 1, (
            f"expected 1 N→K row, got: {result.to_dict('records')}"
        )
        assert nk_rows.iloc[0]["verdict"] == "novel_allele_candidate"
        assert nk_rows.iloc[0]["n_clones"] == 13


class TestStoredAndRecomputedCoexist:
    """Cohorts can carry both stored divergences (from auto_detected
    clones where the picker called an allele but the contig still
    disagreed at non-distinguishing positions) AND recomputed
    divergences (from no-call clones where the picker punted). At
    the same ``(gene, position, expected_aa, observed_aa)``, both
    sources must aggregate into a single row."""

    def test_same_substitution_from_both_sources_aggregates(self):
        rows = []
        # 10 auto_detected clones with stored N→K divergence at TRAC pos 3.
        for i in range(10):
            rows.append({
                "alpha_c_gene": "TRAC",
                "alpha_c_gene_canonical": "TRAC",
                "alpha_allele_called": "TRAC*01",
                "alpha_allele_called_reason": ALLELE_REASON_AUTO_DETECTED,
                "alpha_allele_divergence_positions": "3:N->K",
                "alpha_observed_constant_aa_start": (
                    HUMAN_TRAC_AA[:2] + "K" + HUMAN_TRAC_AA[3:15]
                ),
                "alpha_v_gene": f"TRAV{1 + (i % 5)}",
                "samples": f"S{1 + (i % 2)}",
            })
        # 5 no-call clones with same observed substitution; recompute
        # against the default canonical (TRAC) yields the same N→K row.
        for i in range(5):
            rows.append({
                "alpha_c_gene": "TRAC",
                "alpha_c_gene_canonical": "TRAC",
                "alpha_allele_called": None,
                "alpha_allele_called_reason": ALLELE_REASON_DIVERGENT_CONTIG,
                "alpha_allele_divergence_positions": None,
                "alpha_observed_constant_aa_start": (
                    HUMAN_TRAC_AA[:2] + "K" + HUMAN_TRAC_AA[3:15]
                ),
                "alpha_v_gene": f"TRAV{6 + (i % 5)}",
                "samples": f"S{3 + (i % 2)}",
            })
        df = pd.DataFrame(rows)
        result = detect_novel_alleles(
            df, min_pct=0.05, min_v_spread=3, min_samples=2,
            min_cohort_size=0,
        )
        nk_rows = result[
            (result["expected_aa"] == "N") & (result["observed_aa"] == "K")
        ]
        # Single row aggregating both sources — n_clones = 10 + 5 = 15.
        assert len(nk_rows) == 1, (
            f"stored + recomputed should aggregate into one row; got "
            f"{nk_rows.to_dict('records')}"
        )
        assert nk_rows.iloc[0]["n_clones"] == 15
        # V-gene + sample union across both sources.
        assert nk_rows.iloc[0]["n_v_genes"] == 10  # TRAV1..TRAV10
        assert nk_rows.iloc[0]["n_samples"] == 4   # S1..S4


class TestMinCohortSizeGate:
    """``min_cohort_size`` skips the aggregator on small cohorts so
    ``allele_audit_report`` doesn't pay an unnecessary tax in the
    interactive assemble path."""

    def test_below_threshold_returns_empty(self):
        # 5 clones — well below the default min_cohort_size of 20.
        rows = [{
            "alpha_c_gene": "TRAC",
            "alpha_c_gene_canonical": "TRAC",
            "alpha_allele_called": None,
            "alpha_allele_called_reason": ALLELE_REASON_DIVERGENT_CONTIG,
            "alpha_allele_divergence_positions": None,
            "alpha_observed_constant_aa_start": (
                HUMAN_TRAC_AA[:2] + "K" + HUMAN_TRAC_AA[3:15]
            ),
            "alpha_v_gene": f"TRAV{i + 1}",
            "samples": f"S{i + 1}",
        } for i in range(5)]
        df = pd.DataFrame(rows)
        result = detect_novel_alleles(df, min_cohort_size=20)
        assert result.empty

    def test_zero_threshold_evaluates_anyway(self):
        # 5 clones but min_cohort_size=0 → runs.
        rows = [{
            "alpha_c_gene": "TRAC",
            "alpha_c_gene_canonical": "TRAC",
            "alpha_allele_called": None,
            "alpha_allele_called_reason": ALLELE_REASON_DIVERGENT_CONTIG,
            "alpha_allele_divergence_positions": None,
            "alpha_observed_constant_aa_start": (
                HUMAN_TRAC_AA[:2] + "K" + HUMAN_TRAC_AA[3:15]
            ),
            "alpha_v_gene": f"TRAV{i + 1}",
            "samples": f"S{i + 1}",
        } for i in range(5)]
        df = pd.DataFrame(rows)
        result = detect_novel_alleles(df, min_cohort_size=0)
        assert not result.empty


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
