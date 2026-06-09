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

"""Tests for #113: multi-allele canonical with per-clone auto-detect.

The new architecture (Option D from the issue thread):
  - Multi-allele FASTA ships per-allele mature C-region sequences
    (TRAC*01, TRBC1*01, TRBC2*01, TRBC2*03).
  - ``_pick_best_allele`` scores each packaged allele against the
    contig translation past VDJ and picks the best.
  - ``assemble_full_sequences`` takes per-gene ``{gene}_allele=...``
    kwargs: ``"auto"`` (default) runs the picker; an explicit allele
    label like ``"01"`` or ``"03"`` forces that allele.
  - Audit columns recorded per clone:
      * ``{chain}_allele_called`` — e.g. ``"TRBC2*01"`` or None.
      * ``{chain}_allele_score`` — float ∈ [0, 1] or None.
      * ``{chain}_allele_alternatives`` — ``"allele:score;..."`` for
        the other compatible alleles, sorted by score desc.
"""

from __future__ import annotations

import pandas as pd

from tcrsift.assemble import (
    HUMAN_CONSTANT_ALLELES,
    HUMAN_PREFERRED_CODONS,
    HUMAN_TRAC_AA,
    HUMAN_TRBC1_AA,
    HUMAN_TRBC2_AA,
    _pick_best_allele,
    _score_allele_against_contig,
    assemble_full_sequences,
)

# ---------------------------------------------------------------------------
# Packaging shape


class TestMultiAlleleFasta:
    def test_all_three_genes_present(self):
        assert set(HUMAN_CONSTANT_ALLELES) == {"TRAC", "TRBC1", "TRBC2"}

    def test_trac_has_only_01(self):
        """No human TRAC*02-protein form has been documented; ship *01 only."""
        assert set(HUMAN_CONSTANT_ALLELES["TRAC"]) == {"01"}

    def test_trbc1_has_only_01(self):
        """TRBC1*02 protein-level references are non-human (gorilla / ferret /
        orangutan per IMGT); ship *01 only until a human *02 reference
        is verified."""
        assert set(HUMAN_CONSTANT_ALLELES["TRBC1"]) == {"01"}

    def test_trbc2_has_both_protein_forms(self):
        """Two distinct human TRBC2 protein forms — *01 (E at pos 9,
        encoded by IMGT *01/*02/*04) and *03 (K at pos 9)."""
        assert set(HUMAN_CONSTANT_ALLELES["TRBC2"]) == {"01", "03"}

    def test_default_human_trbc2_aa_is_01_form(self):
        """As of 2.3, ``HUMAN_TRBC2_AA`` is the *01-protein form (the
        major-allele in humans). Pre-2.3 it was *03."""
        assert HUMAN_TRBC2_AA == HUMAN_CONSTANT_ALLELES["TRBC2"]["01"]
        assert HUMAN_TRBC2_AA[8] == "E"

    def test_trbc2_alleles_differ_only_at_pos_9(self):
        """The whole point of carrying both protein forms is that they
        differ at exactly one mature position."""
        a01 = HUMAN_CONSTANT_ALLELES["TRBC2"]["01"]
        a03 = HUMAN_CONSTANT_ALLELES["TRBC2"]["03"]
        assert len(a01) == len(a03)
        diffs = [i for i in range(len(a01)) if a01[i] != a03[i]]
        assert diffs == [8]
        assert a01[8] == "E"
        assert a03[8] == "K"


# ---------------------------------------------------------------------------
# Scoring + picker unit tests


class TestScoreAlleleAgainstContig:
    def test_full_agreement(self):
        n_agree, n_compared = _score_allele_against_contig("EDLNK", "EDLNK")
        assert (n_agree, n_compared) == (5, 5)

    def test_one_mismatch(self):
        n_agree, n_compared = _score_allele_against_contig("EDLXK", "EDLNK")
        assert (n_agree, n_compared) == (4, 5)

    def test_contig_shorter_than_allele(self):
        n_agree, n_compared = _score_allele_against_contig(
            "EDL", "EDLNKVFPP",
        )
        assert (n_agree, n_compared) == (3, 3)

    def test_allele_shorter_than_contig(self):
        n_agree, n_compared = _score_allele_against_contig(
            "EDLNKVFPPEXTRA", "EDLNK",
        )
        assert (n_agree, n_compared) == (5, 5)

    def test_empty_inputs(self):
        assert _score_allele_against_contig("", "EDLNK") == (0, 0)
        assert _score_allele_against_contig("EDLNK", "") == (0, 0)


class TestPickBestAllele:
    def test_picks_trbc2_01_when_contig_has_e_at_pos_9(self):
        """Donor contig codes E at pos 9 (the most common case) → *01 wins."""
        contig_aa = "DLKNVFPPE" + "VAVFEPSEAE"  # E at pos 9 (idx 8)
        best, score, all_scores = _pick_best_allele(
            contig_aa, "TRBC2",
        )
        assert best == "01"
        assert score == 1.0
        # Both alleles scored; *03 has a single mismatch at pos 9.
        assert all_scores["01"] == 1.0
        assert all_scores["03"] < 1.0

    def test_picks_trbc2_03_when_contig_has_k_at_pos_9(self):
        contig_aa = "DLKNVFPPK" + "VAVFEPSEAE"  # K at pos 9
        best, score, all_scores = _pick_best_allele(
            contig_aa, "TRBC2",
        )
        assert best == "03"
        assert score == 1.0
        assert all_scores["03"] == 1.0
        assert all_scores["01"] < 1.0

    def test_picks_trbc1_01_when_only_one_allele(self):
        """TRBC1 has only one packaged allele; picker returns it."""
        contig_aa = HUMAN_TRBC1_AA[:15]
        best, score, _ = _pick_best_allele(contig_aa, "TRBC1")
        assert best == "01"
        assert score == 1.0

    def test_returns_no_call_when_coverage_too_thin(self):
        """When the contig translation is shorter than ``min_compared``,
        the picker returns ``None`` so the caller can fall back to the
        user's default without confusing the no-decision case with a
        confident first-allele pick."""
        best, score, all_scores = _pick_best_allele(
            "DL", "TRBC2", min_compared=5,
        )
        assert best is None
        assert score == 0.0
        assert all_scores == {}

    def test_5_codon_contig_ties_under_default_min(self):
        """The default ``min_compared`` of 10 is past the deepest
        TRBC2 distinguishing position (pos 9, idx 8). With only 5
        codons of contig coverage the picker's window doesn't include
        the distinguishing residue, so it correctly no-decides instead
        of silently breaking the tie by FASTA order."""
        contig_aa = "DLKNV"  # 5 residues; agrees with both *01 and *03
        best, score, all_scores = _pick_best_allele(contig_aa, "TRBC2")
        assert best is None
        assert score == 0.0
        assert all_scores == {}

    def test_10_codon_contig_decides_with_default_min(self):
        """At exactly 10 codons (the default floor) the picker can
        distinguish *01 (E at idx 8) from *03 (K at idx 8)."""
        # *01-form (E at pos 9 / idx 8) for the first 10 residues:
        contig_aa = HUMAN_CONSTANT_ALLELES["TRBC2"]["01"][:10]
        best, score, _ = _pick_best_allele(contig_aa, "TRBC2")
        assert best == "01"
        assert score == 1.0


# ---------------------------------------------------------------------------
# End-to-end through assemble_full_sequences


def _build_clone(beta_contig_seq: str, *, alpha_leader_aa: str = "M" + "A" * 19):
    """Minimal fixture: one β-only clone with a synthetic VDJ + the
    given contig past VDJ."""
    vdj_beta = "CASS" + "G" * 60 + "VETA"
    vdj_beta_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_beta)
    leader_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in alpha_leader_aa)
    contig = leader_nt + vdj_beta_nt + beta_contig_seq
    vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
    df = pd.DataFrame([{
        "CDR3ab": "c1",
        "CDR3_alpha": vdj_alpha, "CDR3_beta": vdj_beta,
        "VDJ_alpha_aa": vdj_alpha, "VDJ_beta_aa": vdj_beta,
        "VDJ_alpha_nt": "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_alpha),
        "VDJ_beta_nt": vdj_beta_nt + "A",  # +1 overshoot for #91 trim
        "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC2",
        "beta_j_gene": "TRBJ2-1",  # TRBJ2 → TRBC2 per J-family parity (#90)
        "samples": "S1",
        "beta_contig_ids": "contig_b1",
    }])
    return df, contig


def _make_contigs_dir(tmp_path, contig_seq):
    sample_dir = tmp_path / "S1"
    sample_dir.mkdir(parents=True, exist_ok=True)
    (sample_dir / "filtered_contig.fasta").write_text(
        f">contig_b1\n{contig_seq}\n"
    )
    return tmp_path


class TestE2eAlleleAutoDetect:
    """End-to-end: the auto-picker chooses the right allele based on
    the donor's contig, and the audit columns are populated."""

    def _contig_for_trbc2(self, allele: str, junction_codon: str = "GAG"):
        """Build the post-VDJ contig bytes that encode the given allele's
        first ~15 mature residues. Includes the universal-E junction
        codon at the front (``junction_codon='GAG'``)."""
        canonical = HUMAN_CONSTANT_ALLELES["TRBC2"][allele][:15]
        canonical_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in canonical)
        return junction_codon + canonical_nt

    def test_auto_picks_trbc2_01_for_e_at_pos_9_donor(self, tmp_path):
        """B1-2/B1-3-style donor: 100% TRBC2*01 (E at pos 9)."""
        df, contig = _build_clone(self._contig_for_trbc2("01"))
        contigs_dir = _make_contigs_dir(tmp_path, contig)
        out = assemble_full_sequences(
            df, contigs_dir=str(contigs_dir),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        assert out["beta_allele_called"].iloc[0] == "TRBC2*01"
        assert out["beta_allele_score"].iloc[0] == 1.0
        # The assembled constant carries E at mature pos 9 — the
        # junction-prepend means stored idx 0 is the junction E, so
        # the C-exon body starts at idx 1; pos 9 of mature TRBC2 is at
        # idx 9 of the assembled constant.
        assert out["beta_constant_aa"].iloc[0][9] == "E"

    def test_auto_picks_trbc2_03_for_k_at_pos_9_donor(self, tmp_path):
        """Minority-allele donor: TRBC2*03 (K at pos 9)."""
        df, contig = _build_clone(self._contig_for_trbc2("03"))
        contigs_dir = _make_contigs_dir(tmp_path, contig)
        out = assemble_full_sequences(
            df, contigs_dir=str(contigs_dir),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        assert out["beta_allele_called"].iloc[0] == "TRBC2*03"
        assert out["beta_allele_score"].iloc[0] == 1.0
        # The assembled constant carries K at pos 9 — the contig's
        # actual sequence, no longer stomped to E (the #113 bug).
        assert out["beta_constant_aa"].iloc[0][9] == "K"

    def test_alternatives_column_records_other_allele_scores(self, tmp_path):
        """The audit trail surfaces the runner-up alleles so users can
        see how close the call was."""
        df, contig = _build_clone(self._contig_for_trbc2("01"))
        contigs_dir = _make_contigs_dir(tmp_path, contig)
        out = assemble_full_sequences(
            df, contigs_dir=str(contigs_dir),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        alternatives = out["beta_allele_alternatives"].iloc[0]
        # Winner first, runner-up second.
        assert "TRBC2*01:" in alternatives
        assert "TRBC2*03:" in alternatives
        # Winner has the higher score.
        parts = dict(p.split(":") for p in alternatives.split(";"))
        assert float(parts["TRBC2*01"]) > float(parts["TRBC2*03"])


class TestE2eExplicitAlleleOverride:
    """User-supplied allele kwargs take precedence over auto-detect."""

    def test_forced_trbc2_03_even_when_contig_says_01(self, tmp_path):
        """Edge case: user knows their donor carries the minority allele
        but their contig was sparse or noisy. ``trbc2_allele='03'``
        forces the picker's hand."""
        df, contig = _build_clone(
            # Contig says *01 (E at pos 9), but user forces *03.
            "GAG" + "".join(
                HUMAN_PREFERRED_CODONS[r]
                for r in HUMAN_CONSTANT_ALLELES["TRBC2"]["01"][:15]
            ),
        )
        contigs_dir = _make_contigs_dir(tmp_path, contig)
        out = assemble_full_sequences(
            df, contigs_dir=str(contigs_dir),
            trbc2_allele="03",
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        assert out["beta_allele_called"].iloc[0] == "TRBC2*03"
        assert out["beta_constant_aa"].iloc[0][9] == "K"

    def test_forced_trbc2_01_when_contig_says_03(self, tmp_path):
        df, contig = _build_clone(
            "GAG" + "".join(
                HUMAN_PREFERRED_CODONS[r]
                for r in HUMAN_CONSTANT_ALLELES["TRBC2"]["03"][:15]
            ),
        )
        contigs_dir = _make_contigs_dir(tmp_path, contig)
        out = assemble_full_sequences(
            df, contigs_dir=str(contigs_dir),
            trbc2_allele="01",
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        assert out["beta_allele_called"].iloc[0] == "TRBC2*01"
        assert out["beta_constant_aa"].iloc[0][9] == "E"


class TestE2eNoContigFallback:
    """When no contigs are loaded, the picker can't run; the assembly
    uses the default canonical (*01 for TRBC2 since 2.3)."""

    def test_no_contigs_uses_trbc2_01_default(self):
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha, "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha, "VDJ_beta_aa": vdj_beta,
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC2",
            "beta_j_gene": "TRBJ2-1",  # TRBJ2 → TRBC2 per J-family parity (#90)
            "samples": "S1",
        }])
        out = assemble_full_sequences(
            df, alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        # No allele call (picker couldn't run); default TRBC2*01 used.
        assert out["beta_allele_called"].iloc[0] is None
        # The constant_aa is the default TRBC2 (= *01) = ``HUMAN_TRBC2_AA``,
        # now prefixed with the invariant β junction E from the canonical
        # fallback (no contig available) so the chain isn't 1 aa short (#235).
        assert out["beta_constant_aa"].iloc[0] == "E" + HUMAN_TRBC2_AA
        assert out["beta_junction_residue"].iloc[0] == "E"


class TestE2eAlphaTracPicker:
    """The picker runs for α-chain TRAC too (single packaged allele).
    Audit columns should reflect the call."""

    def test_alpha_records_trac_01(self, tmp_path):
        # Mirror the β fixture but populate the α contig with the
        # universal "AAT" junction (N) + back-translated TRAC.
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_alpha_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_alpha)
        leader_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in "M" + "A" * 19)
        trac_start_nt = "".join(
            HUMAN_PREFERRED_CODONS[r] for r in HUMAN_TRAC_AA[:15]
        )
        alpha_contig = leader_nt + vdj_alpha_nt + "AAT" + trac_start_nt
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        vdj_beta_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_beta)
        beta_contig = (
            leader_nt + vdj_beta_nt + "GAG"
            + "".join(HUMAN_PREFERRED_CODONS[r] for r in HUMAN_TRBC2_AA[:15])
        )
        sample_dir = tmp_path / "S1"
        sample_dir.mkdir(parents=True)
        (sample_dir / "filtered_contig.fasta").write_text(
            f">contig_a1\n{alpha_contig}\n>contig_b1\n{beta_contig}\n"
        )
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha, "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha, "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": vdj_alpha_nt + "A",
            "VDJ_beta_nt": vdj_beta_nt + "A",
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC2",
            "beta_j_gene": "TRBJ2-1",
            "samples": "S1",
            "alpha_contig_ids": "contig_a1",
            "beta_contig_ids": "contig_b1",
        }])
        out = assemble_full_sequences(
            df, contigs_dir=str(tmp_path),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        assert out["alpha_allele_called"].iloc[0] == "TRAC*01"
        assert out["alpha_allele_score"].iloc[0] == 1.0
        # Only one allele packaged for TRAC, so alternatives is just *01.
        assert "TRAC*01" in out["alpha_allele_alternatives"].iloc[0]


class TestInvalidAlleleOverride:
    """User-supplied invalid allele labels (typos, non-existent
    alleles) must surface as QC warnings — not silently fall through
    to the default."""

    def _basic_fixture(self, tmp_path):
        leader_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in "M" + "A" * 19)
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        vdj_beta_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_beta)
        contig = (
            leader_nt + vdj_beta_nt + "GAG"
            + "".join(HUMAN_PREFERRED_CODONS[r] for r in HUMAN_TRBC2_AA[:15])
        )
        sample_dir = tmp_path / "S1"
        sample_dir.mkdir(parents=True)
        (sample_dir / "filtered_contig.fasta").write_text(
            f">contig_b1\n{contig}\n"
        )
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha, "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha, "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_alpha),
            "VDJ_beta_nt": vdj_beta_nt + "A",
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC2",
            "beta_j_gene": "TRBJ2-1",
            "samples": "S1",
            "beta_contig_ids": "contig_b1",
        }])
        return df

    def test_invalid_label_emits_qc_warning(self, tmp_path):
        df = self._basic_fixture(tmp_path)
        out = assemble_full_sequences(
            df, contigs_dir=str(tmp_path),
            trbc2_allele="99",  # not in the packaged pool
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        # Falls back to default; no allele call recorded.
        assert out["beta_allele_called"].iloc[0] is None
        # QC warning surfaces the invalid label.
        qc = out["qc_warnings"].iloc[0] or []
        assert any(
            "invalid allele override" in m and "99" in m for m in qc
        ), f"expected invalid-allele warning, got {qc}"

    def test_short_unzero_padded_label_treated_as_invalid(self, tmp_path):
        """Pin: ``trbc2_allele='3'`` is NOT the same as ``'03'`` —
        catches the common typo."""
        df = self._basic_fixture(tmp_path)
        out = assemble_full_sequences(
            df, contigs_dir=str(tmp_path),
            trbc2_allele="3",  # typo for "03"
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        qc = out["qc_warnings"].iloc[0] or []
        assert any("invalid allele override" in m for m in qc)


class TestConstantSourceFromDataIgnoresAlleleKwargs:
    """``constant_source='from-data'`` reads the constant from the
    input frame and bypasses the canonical FASTA entirely. Allele
    overrides on that combination are silently meaningless; we emit
    a top-level warning so the mismatch surfaces."""

    def test_warning_when_combined(self, caplog):
        import logging

        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha, "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha, "VDJ_beta_aa": vdj_beta,
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC2",
            "beta_j_gene": "TRBJ2-1",
            "samples": "S1",
            "alpha_constant_aa": HUMAN_TRAC_AA,
            "beta_constant_aa": HUMAN_TRBC2_AA,
        }])
        with caplog.at_level(logging.WARNING, logger="tcrsift.assemble"):
            assemble_full_sequences(
                df,
                constant_source="from-data",
                trbc2_allele="03",  # ignored, should warn
                alpha_leader=None, beta_leader=None,
                verbose=False, show_progress=False,
            )
        assert any(
            "from-data" in r.message and "ignores" in r.message
            for r in caplog.records
        )


class TestConstantSourceDecoupledFromAlleleCall:
    """The ``{chain}_constant_source`` audit string and the new
    ``{chain}_allele_called`` are decoupled — the source string just
    names the gene plus the contig/canonical blend, while the allele
    info lives only in the new audit columns. Downstream code parsing
    ``_constant_source`` to identify "which canonical" must read
    ``_allele_called`` instead."""

    def test_source_string_does_not_contain_allele(self, tmp_path):
        # Build the same fixture as the auto-detect test for *01.
        leader_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in "M" + "A" * 19)
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        vdj_beta_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_beta)
        contig = (
            leader_nt + vdj_beta_nt + "GAG"
            + "".join(HUMAN_PREFERRED_CODONS[r] for r in HUMAN_TRBC2_AA[:15])
        )
        sample_dir = tmp_path / "S1"
        sample_dir.mkdir(parents=True)
        (sample_dir / "filtered_contig.fasta").write_text(
            f">contig_b1\n{contig}\n"
        )
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha, "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha, "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_alpha),
            "VDJ_beta_nt": vdj_beta_nt + "A",
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC2",
            "beta_j_gene": "TRBJ2-1",
            "samples": "S1",
            "beta_contig_ids": "contig_b1",
        }])
        out = assemble_full_sequences(
            df, contigs_dir=str(tmp_path),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        source = out["beta_constant_source"].iloc[0]
        assert "TRBC2" in source  # the gene
        # No allele label in source — allele info is in its own column.
        assert "*01" not in source and "*03" not in source
        # The allele info is in its own column.
        assert out["beta_allele_called"].iloc[0] == "TRBC2*01"


class TestNTFrameVerification:
    """The picker's NT-level first-codon frame check (#115 review)
    catches frame-off-by-N bugs upstream where AA scoring would
    falsely succeed."""

    def test_first_codon_translates_correctly(self):
        from tcrsift.assemble import _verify_first_codon_frame

        trbc2_01 = HUMAN_CONSTANT_ALLELES["TRBC2"]["01"]
        # First codon of TRBC2 mature is D (codon GAT or GAC).
        # We use the codon-optimized choice from HUMAN_PREFERRED_CODONS.
        first_codon_nt = HUMAN_PREFERRED_CODONS[trbc2_01[0]]
        assert _verify_first_codon_frame(first_codon_nt, trbc2_01) is True

    def test_frame_shift_caught(self):
        """A frame-off-by-1 contig has a different first codon — the
        check rejects it even though downstream AA scoring on the
        shifted bytes might be high."""
        from tcrsift.assemble import _verify_first_codon_frame

        # Allele expects D at residue 0; contig's 'first codon' is
        # actually ATG (M) — frame error.
        assert _verify_first_codon_frame("ATG" + "X" * 30, "DLKN") is False

    def test_too_short_contig_fails_check(self):
        from tcrsift.assemble import _verify_first_codon_frame

        # < 3 nt cannot form a codon.
        assert _verify_first_codon_frame("GA", "DLKN") is False
        assert _verify_first_codon_frame("", "DLKN") is False

    def test_empty_allele_fails_check(self):
        from tcrsift.assemble import _verify_first_codon_frame

        assert _verify_first_codon_frame("GAT", "") is False
