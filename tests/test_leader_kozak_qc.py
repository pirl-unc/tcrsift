# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Signal-peptide QC + Kozak start-selection helpers (#263).

The Kozak cases are the real TRAV16 over-capture: three in-frame ATGs whose
contexts the user scored by hand — M0/M13 are weak 5'-UTR AUGs, M35 is the
adequate (-3 purine) true start. The longest-ORF rule grabs M0 (54-aa leader);
Kozak ranking must instead prefer M35 (~19-aa textbook signal peptide).
"""

from __future__ import annotations

from tcrsift.assemble import (
    DEFAULT_LEADERS,
    KOZAK_ADEQUATE,
    SP_LENGTH_MAX,
    SP_LENGTH_MIN,
    _has_h_region,
    _kozak_correct_leader,
    _kozak_score,
    leader_qc,
)


class TestLeaderQC:
    def test_curated_defaults_all_ok(self):
        # The hand-verified leaders must all pass (M-start, in-window, h-region).
        for name, spec in DEFAULT_LEADERS.items():
            assert leader_qc(spec["aa"]) == "ok", name

    def test_missing_and_no_met(self):
        assert leader_qc(None) == "missing"
        assert leader_qc("") == "missing"
        assert leader_qc("ALPVTALLLLLLLLA") == "no_met"

    def test_internal_stop(self):
        assert leader_qc("MALL*LLLLAARP") == "internal_stop"

    def test_length_window(self):
        assert leader_qc("M" + "L" * 5) == "too_short"        # 6 aa
        assert leader_qc("M" + "L" * 60) == "too_long"        # 61 aa
        # The 54-aa over-capture from B1-2 is flagged too_long.
        assert leader_qc("M" + "A" * 53) == "too_long"

    def test_no_h_region(self):
        # M-start, in-window length, but no hydrophobic core → flagged.
        assert leader_qc("MKKEKDEKRNKQDE") == "no_h_region"

    def test_window_constants_sane(self):
        assert SP_LENGTH_MIN == 12 and SP_LENGTH_MAX == 25


class TestHRegion:
    def test_hydrophobic_core_detected(self):
        assert _has_h_region("MALPVTALLLPLALLLHAARP")  # CD8A leader

    def test_polar_stretch_rejected(self):
        assert not _has_h_region("MKKEKDEKRNKQDE")


class TestKozakScoring:
    """Real TRAV16 contexts (user-scored by hand): 0 weak · >=2 adequate."""

    # bases -6..-1 then ATG then +4; ATG's A sits at index 6.
    # M0 : CTATAC ATG C  (-3=T, +4=C) weak     → 0
    # M13: GTGTCT ATG T  (-3=T, +4=T) weak     → 0
    # M35: TCAGAC ATG A  (-3=G, +4=A) adequate → 2 (-3 purine, +4 not G)
    M0 = "CTATAC" + "ATG" + "C"
    M13 = "GTGTCT" + "ATG" + "T"
    M35 = "TCAGAC" + "ATG" + "A"

    def _score(self, ctx):
        return _kozak_score(ctx, 6)

    def test_trav16_scores_match_hand_calibration(self):
        assert self._score(self.M0) == 0      # fails both key positions
        assert self._score(self.M13) == 0     # fails both
        assert self._score(self.M35) == 2     # -3=G purine (adequate), +4=A

    def test_only_m35_clears_the_adequate_bar(self):
        assert self._score(self.M35) >= KOZAK_ADEQUATE
        assert self._score(self.M0) < KOZAK_ADEQUATE
        assert self._score(self.M13) < KOZAK_ADEQUATE

    def test_strong_context_scores_higher(self):
        # -3=A (preferred purine) AND +4=G → strong (>=4).
        strong = _kozak_score("AAA" + "A" + "CC" + "ATG" + "G", 6)
        assert strong >= 4

    def test_a_preferred_over_g_at_minus3(self):
        a3 = _kozak_score("AAA" + "A" + "CC" + "ATG" + "C", 6)
        g3 = _kozak_score("AAA" + "G" + "CC" + "ATG" + "C", 6)
        assert a3 > g3 >= KOZAK_ADEQUATE

    def test_incomplete_5p_context_returns_negative(self):
        # ATG hard against the 5' end → -1, never chosen over a real start.
        assert _kozak_score("ATGGCC", 0) == -1


class TestLeakyScan:
    """Leaky-scan picks the first ADEQUATE in-frame AUG, not the first one.

    Reconstructs the TRAV16 collapse: three in-frame Mets at aa offsets 0/13/35,
    where M0/M13 are weak (pyrimidine -3, non-G +4) and only M35 has a -3 purine
    — so the 54-aa over-capture must collapse to the 19-aa real leader.
    """

    # Filler codon with a PYRIMIDINE first base (Leu, CTG): an ATG immediately
    # after it has -3='C' (pyrimidine) and, with the next CTG, +4='C' → weak.
    L = "CTG"
    # One Thr codon (ACT, purine first base) placed right before M35 gives it
    # the -3='A' purine that clears the adequate bar.
    THR = "ACT"
    VDJ_NT = "GAATTCGAATTCGAA"  # translates to EFEFE
    VDJ_AA = "EFEFE"

    def _build_contig(self):
        pre = "TTTTTT"                       # M0 upstream: -3='T' (weak)
        nt = (
            pre
            + "ATG" + self.L * 12            # M0  + 12 L           (k=0)
            + "ATG" + self.L * 20            # M13 + 20 L           (k=13)
            + self.THR                        #   Thr right before M35 (k=34)
            + "ATG" + self.L * 18            # M35 + 18 L           (k=35)
            + self.VDJ_NT
        )
        return nt

    def test_leaky_scan_collapses_to_m35(self):
        from tcrsift.assemble import find_longest_orf

        nt = self._build_contig()
        translated, offset, _ = find_longest_orf(nt)
        orig_leader = translated.split(self.VDJ_AA)[0]
        assert len(orig_leader) == 54  # the over-capture (M0 start)

        leader_aa, leader_nt, score, kind = _kozak_correct_leader(
            nt, offset, orig_leader
        )
        # First adequate-context AUG is M35 → 19-aa leader, not 54.
        assert kind == "kozak"
        assert len(leader_aa) == 19
        assert leader_aa == orig_leader[35:]
        assert score >= KOZAK_ADEQUATE
        # NT trimmed in lockstep with the AA.
        assert len(leader_nt) == 19 * 3

    def test_no_adequate_aug_flags_weak_kozak(self):
        # A single weak-context ATG (pyrimidine -3, non-G +4): nothing clears the
        # bar → weak_kozak, and we keep the longest-ORF leader (visible, flagged).
        from tcrsift.assemble import find_longest_orf

        nt = "TTTTTT" + "ATG" + self.L * 18 + self.VDJ_NT
        translated, offset, _ = find_longest_orf(nt)
        orig_leader = translated.split(self.VDJ_AA)[0]
        leader_aa, leader_nt, score, kind = _kozak_correct_leader(
            nt, offset, orig_leader
        )
        assert kind == "weak_kozak"
        assert leader_aa == orig_leader  # fell back to the longest-ORF start


class TestFromContigIntegration:
    """End-to-end: leaky-scan QC columns + curated fallback through assembly."""

    L = "CTG"
    THR = "ACT"
    VDJ_NT = "GAATTCGAATTCGAA"
    VDJ_AA = "EFEFE"

    # Over-capture contig: M0/M13 weak, M35 adequate → leaky-scan picks M35.
    OVERCAPTURE = (
        "TTTTTT"
        + "ATG" + "CTG" * 12
        + "ATG" + "CTG" * 20
        + "ACT"
        + "ATG" + "CTG" * 18
        + "GAATTCGAATTCGAA"
    )
    # Weak-only contig: one weak ATG, nothing clears the bar → weak_kozak.
    WEAK_ONLY = "TTTTTT" + "ATG" + "CTG" * 18 + "GAATTCGAATTCGAA"

    def _row(self, contig_id):
        import pandas as pd

        return pd.DataFrame([{
            "CDR3ab": "EFEFE_CASS",
            "CDR3_alpha": "EFEFE", "CDR3_beta": "CASS",
            "VDJ_alpha_aa": self.VDJ_AA, "VDJ_alpha_nt": self.VDJ_NT,
            "VDJ_beta_aa": "CASS" + "G" * 40 + "VETA",
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1",
            "alpha_j_gene": "TRAJ45", "beta_j_gene": "TRBJ1-1",
            "samples": "S1", "alpha_contig_ids": contig_id,
        }])

    def _assemble(self, tmp_path, contig_id, contig_seq, **kw):
        from tcrsift.assemble import assemble_full_sequences

        d = tmp_path / "S1" / "vdj_t"
        d.mkdir(parents=True)
        (d / "filtered_contig.fasta").write_text(f">{contig_id}\n{contig_seq}\n")
        out = assemble_full_sequences(
            self._row(contig_id), contigs_dir=str(tmp_path),
            sample_name_from="grandparent", alpha_leader="from_contig",
            beta_leader=None, include_constant=False,
            verbose=False, show_progress=False, **kw,
        )
        return out.iloc[0]

    def test_overcapture_collapses_to_19aa_ok(self, tmp_path):
        r = self._assemble(tmp_path, "c1", self.OVERCAPTURE)
        assert r["alpha_leader_aa"] == "M" + "L" * 18      # 19-aa real leader
        assert r["alpha_leader_qc"] == "ok"
        assert r["alpha_leader_source"] == "contig_kozak"
        assert r["alpha_leader_kozak_score"] >= 2

    def test_weak_only_flagged_weak_kozak(self, tmp_path):
        r = self._assemble(tmp_path, "c2", self.WEAK_ONLY)
        assert r["alpha_leader_qc"] == "weak_kozak"
        assert r["alpha_leader_source"] == "contig_weak_fallback"

    def test_leader_fallback_substitutes_curated_sp(self, tmp_path):
        from tcrsift.assemble import DEFAULT_LEADERS

        r = self._assemble(tmp_path, "c3", self.WEAK_ONLY, leader_fallback="CD8A")
        assert r["alpha_leader_aa"] == DEFAULT_LEADERS["CD8A"]["aa"]
        assert r["alpha_leader_qc"] == "curated_fallback"
        assert r["alpha_leader_source"] == "curated_fallback:CD8A"
