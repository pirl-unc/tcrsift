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
    SP_LENGTH_HARD,
    SP_LENGTH_MAX,
    SP_LENGTH_MIN,
    _classify_leader,
    _has_h_region,
    _kozak_correct_leader,
    _kozak_score,
    leader_qc,
)

# Minimal aa→codon map for building synthetic leader contigs in tests.
_CODON = {
    "M": "ATG", "L": "CTG", "S": "AGC", "P": "CCC", "D": "GAC", "A": "GCC",
    "W": "TGG", "N": "AAC", "T": "ACC", "R": "CGC", "C": "TGC", "V": "GTG",
    "G": "GGC", "E": "GAA", "F": "TTC", "K": "AAG", "Y": "TAC", "H": "CAC",
    "Q": "CAG", "I": "ATC",
}


def _nt(aa: str) -> str:
    return "".join(_CODON[a] for a in aa)


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

    def test_overcapture_reselects_to_m35(self):
        from tcrsift.assemble import find_longest_orf

        nt = self._build_contig()
        translated, offset, _ = find_longest_orf(nt)
        orig_leader = translated.split(self.VDJ_AA)[0]
        assert len(orig_leader) == 54  # the over-capture (M0 start)

        leader_aa, leader_nt, score, source = _kozak_correct_leader(
            nt, offset, orig_leader
        )
        # Over-capture (>25) → first adequate in-range AUG is M35 → 19-aa leader.
        assert source == "contig_kozak_reselected"
        assert len(leader_aa) == 19
        assert leader_aa == orig_leader[35:]
        assert score >= KOZAK_ADEQUATE
        assert len(leader_nt) == 19 * 3  # NT trimmed in lockstep

    def test_in_range_leader_with_weak_start_is_kept_unchanged(self):
        # REGRESSION GUARD: a 19-aa (in-range) leader whose true start is weak
        # Kozak but which HAS an adequate downstream Met must NOT be shortened —
        # only over-captures (>25) are reselected. Keeps a correct leader intact.
        from tcrsift.assemble import find_longest_orf

        nt = (
            "TTTTTT"
            + "ATG" + self.L * 5      # weak M0 (-3='T', +4='C')        (k=0)
            + self.THR                 #   Thr → next Met gets -3 purine  (k=6)
            + "ATG" + self.L * 11     # adequate Met at k=7 (12-aa tail)
            + self.VDJ_NT
        )
        translated, offset, _ = find_longest_orf(nt)
        orig_leader = translated.split(self.VDJ_AA)[0]
        assert len(orig_leader) == 19 <= 25  # in range
        leader_aa, leader_nt, score, source = _kozak_correct_leader(
            nt, offset, orig_leader
        )
        assert source == "contig"          # kept, NOT reselected
        assert leader_aa == orig_leader    # unchanged 19-aa leader
        assert score < KOZAK_ADEQUATE      # its start really is weak (surfaced as QC)

    def test_overcapture_without_in_range_adequate_start_is_kept(self):
        # 40-aa over-capture, single weak ATG, no nested SP at all → kept
        # (the caller flags it too_long); we never fabricate a worse leader.
        from tcrsift.assemble import find_longest_orf

        nt = "TTTTTT" + "ATG" + self.L * 39 + self.VDJ_NT
        translated, offset, _ = find_longest_orf(nt)
        orig_leader = translated.split(self.VDJ_AA)[0]
        assert len(orig_leader) == 40
        leader_aa, leader_nt, score, source = _kozak_correct_leader(
            nt, offset, orig_leader
        )
        assert source == "contig"
        assert leader_aa == orig_leader  # kept the over-capture, unfixable

    # Real TRAV1-2 SP nested at the tail of a 34-aa over-capture; its start is
    # weak-Kozak (why leaky-scan skips it) but it has a clear h-region.
    TRAV1_2_SP = "MWGVFLLYVSMKMGGTT"
    TRAV1_2_SP_NT = (
        "ATGTGGGGCGTGTTCCTGCTGTACGTGAGCATGAAGATGGGCGGCACCACC"
    )

    def test_structural_trim_recovers_weak_kozak_nested_sp(self):
        # 34-aa over-capture = M + 16 filler + 17-aa TRAV1-2 SP. No adequate
        # Kozak start in range → structural (h-region) trim recovers the 17-aa SP.
        from tcrsift.assemble import find_longest_orf

        nt = (
            "TTTTTT" + "ATG" + self.L * 16   # weak M0 + 16 L (no nested Met)
            + self.TRAV1_2_SP_NT             # 17-aa SP, weak-Kozak start at k=17
            + self.VDJ_NT
        )
        translated, offset, _ = find_longest_orf(nt)
        orig_leader = translated.split(self.VDJ_AA)[0]
        assert len(orig_leader) == 34
        leader_aa, leader_nt, score, source = _kozak_correct_leader(
            nt, offset, orig_leader
        )
        assert source == "contig_hregion_trimmed"
        assert leader_aa == self.TRAV1_2_SP   # MWGVFLLYVSMKMGGTT
        assert len(leader_aa) == 17
        assert score < KOZAK_ADEQUATE         # its start really is weak-Kozak


class TestLengthBanding:
    """_classify_leader length bands: ok / weak_kozak_start / long_leader / too_long."""

    OK20 = "M" + "L" * 19  # 20 aa, hydrophobic h-region

    def test_in_range_adequate_start_is_ok(self):
        assert _classify_leader(self.OK20, 3, "contig") == "ok"

    def test_in_range_weak_start_flagged(self):
        assert _classify_leader(self.OK20, 0, "contig") == "weak_kozak_start"

    def test_26_to_30_is_accepted_long_leader(self):
        # TRBV13's real 29-aa leader → long_leader, NOT an error.
        assert _classify_leader("M" + "L" * 28, 0, "contig") == "long_leader"
        assert _classify_leader("M" + "L" * 25, 0, "contig") == "long_leader"  # 26 aa

    def test_over_hard_max_is_too_long(self):
        assert _classify_leader("M" + "L" * 30, 0, "contig") == "too_long"  # 31 aa
        assert SP_LENGTH_HARD == 30

    def test_hregion_trim_and_reselect_sources(self):
        assert _classify_leader("M" + "L" * 16, 0, "contig_hregion_trimmed") == "hregion_trimmed"
        # reselected leader judged purely on shape (in-range → ok)
        assert _classify_leader("M" + "L" * 16, 0, "contig_kozak_reselected") == "ok"


class TestFromContigIntegration:
    """End-to-end: leaky-scan QC columns + curated fallback through assembly."""

    L = "CTG"
    VDJ_NT = "GAATTCGAATTCGAA"
    VDJ_AA = "EFEFE"

    # Over-capture: M0/M13 weak, M35 adequate, 54 aa → reselect to M35 (19 aa).
    OVERCAPTURE = (
        "TTTTTT"
        + "ATG" + "CTG" * 12
        + "ATG" + "CTG" * 20
        + "ACT"
        + "ATG" + "CTG" * 18
        + "GAATTCGAATTCGAA"
    )
    # In-range (19 aa) leader with a weak start → kept + flagged weak_kozak_start.
    INRANGE_WEAK = "TTTTTT" + "ATG" + "CTG" * 18 + "GAATTCGAATTCGAA"
    # Unfixable 40-aa over-capture (no in-range adequate start) → too_long.
    OVERCAPTURE_UNFIXABLE = "TTTTTT" + "ATG" + "CTG" * 39 + "GAATTCGAATTCGAA"

    def _row(self, contig_id, contig_ids=None):
        import pandas as pd

        return pd.DataFrame([{
            "CDR3ab": "EFEFE_CASS",
            "CDR3_alpha": "EFEFE", "CDR3_beta": "CASS",
            "VDJ_alpha_aa": self.VDJ_AA, "VDJ_alpha_nt": self.VDJ_NT,
            "VDJ_beta_aa": "CASS" + "G" * 40 + "VETA",
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1",
            "alpha_j_gene": "TRAJ45", "beta_j_gene": "TRBJ1-1",
            "samples": "S1", "alpha_contig_ids": contig_ids or contig_id,
        }])

    def _assemble(self, tmp_path, contigs, contig_ids=None, **kw):
        # contigs: list of (id, seq); contig_ids: the row's alpha_contig_ids.
        from tcrsift.assemble import assemble_full_sequences

        d = tmp_path / "S1" / "vdj_t"
        d.mkdir(parents=True)
        (d / "filtered_contig.fasta").write_text(
            "".join(f">{cid}\n{seq}\n" for cid, seq in contigs)
        )
        row = self._row(contigs[0][0], contig_ids=contig_ids or contigs[0][0])
        out = assemble_full_sequences(
            row, contigs_dir=str(tmp_path),
            sample_name_from="grandparent", alpha_leader="from_contig",
            beta_leader=None, include_constant=False,
            verbose=False, show_progress=False, **kw,
        )
        return out.iloc[0]

    def test_overcapture_collapses_to_19aa_ok(self, tmp_path):
        r = self._assemble(tmp_path, [("c1", self.OVERCAPTURE)])
        assert r["alpha_leader_aa"] == "M" + "L" * 18      # 19-aa real leader
        assert r["alpha_leader_qc"] == "ok"
        assert r["alpha_leader_source"] == "contig_kozak_reselected"
        assert r["alpha_leader_kozak_score"] >= 2
        assert r["alpha_leader_support"] == "1/1"

    def test_in_range_weak_start_kept_and_flagged(self, tmp_path):
        r = self._assemble(tmp_path, [("c2", self.INRANGE_WEAK)])
        assert r["alpha_leader_aa"] == "M" + "L" * 18      # kept, not shortened
        assert r["alpha_leader_qc"] == "weak_kozak_start"
        assert r["alpha_leader_source"] == "contig"

    def test_unfixable_overcapture_flagged_too_long(self, tmp_path):
        r = self._assemble(tmp_path, [("c3", self.OVERCAPTURE_UNFIXABLE)])
        assert r["alpha_leader_qc"] == "too_long"
        assert r["alpha_leader_source"] == "contig"

    def test_leader_fallback_substitutes_only_implausible(self, tmp_path):
        from tcrsift.assemble import DEFAULT_LEADERS

        # too_long is a bad SP → rejected → switched to the configured fallback.
        r = self._assemble(
            tmp_path, [("c4", self.OVERCAPTURE_UNFIXABLE)], leader_fallback="CD8A"
        )
        assert r["alpha_leader_aa"] == DEFAULT_LEADERS["CD8A"]["aa"]
        assert r["alpha_leader_qc"] == "substituted"
        assert r["alpha_leader_source"] == "curated:CD8A"

    def test_leader_fallback_does_not_touch_weak_kozak_start(self, tmp_path):
        # weak_kozak_start is plausible (in-range, well-shaped) → NOT substituted
        # even when a fallback is configured. The kept leader survives.
        r = self._assemble(
            tmp_path, [("c5", self.INRANGE_WEAK)], leader_fallback="CD8A"
        )
        assert r["alpha_leader_aa"] == "M" + "L" * 18
        assert r["alpha_leader_qc"] == "weak_kozak_start"

    def test_multi_contig_consensus_support(self, tmp_path):
        # Two contigs both yielding the same reselected leader → support 2/2.
        r = self._assemble(
            tmp_path,
            [("c6a", self.OVERCAPTURE), ("c6b", self.OVERCAPTURE)],
            contig_ids="c6a;c6b",
        )
        assert r["alpha_leader_aa"] == "M" + "L" * 18
        assert r["alpha_leader_support"] == "2/2"

    # TRAV1-2: 34-aa over-capture (M + 16 filler + real 17-aa SP), weak-Kozak
    # nested start → structural h-region trim recovers MWGVFLLYVSMKMGGTT.
    TRAV1_2 = (
        "TTTTTT" + "ATG" + "CTG" * 16
        + "ATGTGGGGCGTGTTCCTGCTGTACGTGAGCATGAAGATGGGCGGCACCACC"
        + "GAATTCGAATTCGAA"
    )
    # TRBV13: real 29-aa leader (IMGT TRBV13*01, donor-variant), no nested SP.
    TRBV13 = "GCAGCA" + _nt("MLSPDLPDSAWNTRLLCRVMLCLLGAGSV") + "GAATTCGAATTCGAA"

    def test_trav1_2_overcapture_structural_trim(self, tmp_path):
        r = self._assemble(tmp_path, [("c7", self.TRAV1_2)])
        assert r["alpha_leader_aa"] == "MWGVFLLYVSMKMGGTT"   # real 17-aa SP
        assert r["alpha_leader_qc"] == "hregion_trimmed"
        assert r["alpha_leader_source"] == "contig_hregion_trimmed"

    def test_trbv13_long_leader_accepted_not_trimmed(self, tmp_path):
        r = self._assemble(tmp_path, [("c8", self.TRBV13)])
        # Genuinely-long 29-aa native leader → kept, accepted (not an error).
        assert r["alpha_leader_aa"] == "MLSPDLPDSAWNTRLLCRVMLCLLGAGSV"
        assert r["alpha_leader_qc"] == "long_leader"
        assert r["alpha_leader_source"] == "contig"

    def test_fallback_leaves_long_leader_alone(self, tmp_path):
        # long_leader is accepted (not implausible) → curated fallback won't touch
        # it even when configured.
        r = self._assemble(tmp_path, [("c9", self.TRBV13)], leader_fallback="CD8A")
        assert r["alpha_leader_aa"] == "MLSPDLPDSAWNTRLLCRVMLCLLGAGSV"
        assert r["alpha_leader_qc"] == "long_leader"
