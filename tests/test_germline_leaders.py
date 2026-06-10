# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Germline-anchored signal-peptide selection (#267).

The germline V-leader reference is an alignment anchor: it picks the start/length
of the leader, but the SHIPPED leader keeps the donor's native contig residues
(preserving allelic variants). The anchor is confidence-gated — gene absent /
weird allele / too divergent → returns None and the caller uses the Kozak +
h-region heuristic (#263).
"""

from __future__ import annotations

import pandas as pd

from tcrsift.leaders import (
    _diff_vs_germline,
    characterize_divergence,
    germline_anchor_leader,
    germline_leader_nt,
    germline_vgene_leaders,
)


class TestReference:
    def test_cohort_genes_present(self):
        for gene in ("TRAV1-1", "TRAV1-2", "TRBV13", "TRAV16"):
            assert germline_vgene_leaders(gene), gene

    def test_trav1_2_canonical_leader(self):
        alleles = dict((a, aa) for a, _f, aa in germline_vgene_leaders("TRAV1-2"))
        assert alleles["01"] == "MWGVFLLYVSMKMGGTT"

    def test_trbv13_both_alleles(self):
        alleles = dict((a, aa) for a, _f, aa in germline_vgene_leaders("TRBV13"))
        assert alleles["01"] == "MLSPDLPDSAWNTRLLCHVMLCLLGAVSV"
        assert alleles["03"] == "MLSPDLPDSAWNTRLLCRVMLCLLGAGSV"  # donor's allele

    def test_unknown_gene_empty(self):
        assert germline_vgene_leaders("TRBVZZ") == []
        assert germline_vgene_leaders(None) == []

    def test_allele_suffix_stripped(self):
        # A V call carrying an allele still resolves to the gene's leaders.
        assert germline_vgene_leaders("TRAV1-2*01") == germline_vgene_leaders("TRAV1-2")

    def test_dv_gene_resolves_both_forms(self):
        # IMGT keys dual α/δ genes with a slash; a bare-prefix call must also
        # resolve (review finding — else the germline layer silently skips them).
        slashed = germline_vgene_leaders("TRAV14/DV4")
        assert slashed
        assert germline_vgene_leaders("TRAV14") == slashed
        assert germline_vgene_leaders("TRAV14/DV4*02") == slashed

    def test_multi_length_alleles_present(self):
        # TRBV8-2 has both 18-aa (*01) and 19-aa (*02/*03) leaders.
        lens = {len(aa) for _a, _f, aa in germline_vgene_leaders("TRBV8-2")}
        assert {18, 19} <= lens


class TestDiff:
    def test_substitutions_formatted(self):
        # 1-based residue numbering (M = 1): germline H at 18 → donor R,
        # germline V at 27 → donor G. (The donor's informal 0-based labels were
        # H17/V26 — same two positions.)
        donor = "MLSPDLPDSAWNTRLLCRVMLCLLGAGSV"
        germ = "MLSPDLPDSAWNTRLLCHVMLCLLGAVSV"
        assert _diff_vs_germline(donor, germ) == "H18R;V27G"

    def test_identical_empty(self):
        assert _diff_vs_germline("MWGVFLLYVSMKMGGTT", "MWGVFLLYVSMKMGGTT") == ""


class TestAnchor:
    def test_overcapture_recovers_suffix(self):
        # 34-aa over-capture (17 prepended + real 17-aa SP) → germline anchors
        # the tail; donor residues kept.
        over = "M" + "A" * 16 + "MWGVFLLYVSMKMGGTT"
        out = germline_anchor_leader(over, "TRAV1-2")
        assert out is not None
        leader_aa, allele, g_aa, identity, diff = out
        assert leader_aa == "MWGVFLLYVSMKMGGTT"
        assert identity == 1.0 and diff == ""

    def test_donor_allele_matched_exactly(self):
        # TRBV13 donor (R17/G26) matches *03 exactly, not *01-with-variants.
        donor = "MLSPDLPDSAWNTRLLCRVMLCLLGAGSV"
        out = germline_anchor_leader(donor, "TRBV13")
        assert out is not None
        leader_aa, allele, g_aa, identity, diff = out
        assert leader_aa == donor  # native residues preserved
        assert allele == "03" and identity == 1.0 and diff == ""

    def test_variant_vs_closest_allele_reported(self):
        # A donor 1 aa off from every allele still anchors (>=0.80) and the diff
        # is surfaced. Take TRAV1-2*01 and mutate position 5 F→Y.
        donor = "MWGVYLLYVSMKMGGTT"
        out = germline_anchor_leader(donor, "TRAV1-2")
        assert out is not None
        leader_aa, allele, g_aa, identity, diff = out
        assert leader_aa == donor
        assert diff == "F5Y" and identity > 0.9

    def test_unknown_gene_returns_none(self):
        assert germline_anchor_leader("M" + "A" * 16 + "MWGVFLLYVSMKMGGTT", "TRBVZZ") is None

    def test_too_divergent_returns_none(self):
        # All-K leader shares almost nothing with any TRBV13 allele.
        assert germline_anchor_leader("M" + "K" * 28, "TRBV13") is None

    def test_contig_shorter_than_germline_returns_none(self):
        # A 10-aa leader can't suffix-anchor a 29-aa germline.
        assert germline_anchor_leader("M" + "L" * 9, "TRBV13") is None

    def test_germline_insertion_falls_back_not_truncated(self):
        # A donor germline 1-aa INSERTION lengthens the SP past every reference
        # allele (TRAV1-2's alleles are all 17 aa). The length-matched suffix of
        # the 18-aa donor drops the start Met, so the M-start guard rejects every
        # allele -> returns None, and the caller falls back to the Kozak/h-region
        # heuristic which keeps the donor's longer SP. Guards against silently
        # truncating a real insertion by mis-anchoring it to a shorter allele.
        donor = "MAWGVFLLYVSMKMGGTT"  # TRAV1-2*01 (MWGVFLLYVSMKMGGTT) + 'A' after the start Met
        assert len(donor) == 18
        assert germline_anchor_leader(donor, "TRAV1-2") is None
        # Residual edge (documented, not yet caught): a 1-aa insertion whose
        # residue-2 is itself Met makes the dropped-Met suffix still start with M
        # and match the shorter allele -> would mis-anchor. Requires a duplicated
        # start Met AND high suffix identity, which is biologically rare.

    def test_in_range_leader_anchored_unchanged(self):
        # A leader already at the germline length (no over-capture) anchors and
        # is returned unchanged — the anchor isn't only for over-captures.
        out = germline_anchor_leader("MWGVFLLYVSMKMGGTT", "TRAV1-2")
        assert out is not None
        leader_aa, allele, _g, identity, diff = out
        assert leader_aa == "MWGVFLLYVSMKMGGTT" and identity == 1.0 and diff == ""

    def test_multi_length_allele_picks_matching_length(self):
        # Donor matching the 18-aa TRBV8-2*01 anchors to that allele/length, not
        # a 19-aa allele — the best-identity suffix wins.
        a01 = dict((a, aa) for a, _f, aa in germline_vgene_leaders("TRBV8-2"))["01"]
        assert len(a01) == 18
        out = germline_anchor_leader("MM" + a01, "TRBV8-2")  # 2-aa over-capture
        assert out is not None
        leader_aa, allele, _g, identity, _d = out
        assert leader_aa == a01 and allele == "01" and identity == 1.0


class TestEndToEndAssembly:
    """from_contig assembly emits germline columns + falls back when off-reference."""

    L = "CTG"
    VDJ_NT = "GAATTCGAATTCGAA"
    VDJ_AA = "EFEFE"
    # 34-aa over-capture whose tail is the real TRAV1-2 SP MWGVFLLYVSMKMGGTT.
    TRAV1_2_CONTIG = (
        "TTTTTT" + "ATG" + "CTG" * 16
        + "ATGTGGGGCGTGTTCCTGCTGTACGTGAGCATGAAGATGGGCGGCACCACC"
        + "GAATTCGAATTCGAA"
    )

    def _assemble(self, tmp_path, contig_seq, v_gene, **kw):
        from tcrsift.assemble import assemble_full_sequences

        d = tmp_path / "S1" / "vdj_t"
        d.mkdir(parents=True)
        (d / "filtered_contig.fasta").write_text(f">c1\n{contig_seq}\n")
        row = pd.DataFrame([{
            "CDR3ab": "EFEFE_CASS",
            "CDR3_alpha": "EFEFE", "CDR3_beta": "CASS",
            "VDJ_alpha_aa": self.VDJ_AA, "VDJ_alpha_nt": self.VDJ_NT,
            "VDJ_beta_aa": "CASS" + "G" * 40 + "VETA",
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1",
            "alpha_j_gene": "TRAJ45", "beta_j_gene": "TRBJ1-1",
            "alpha_v_gene": v_gene,
            "samples": "S1", "alpha_contig_ids": "c1",
        }])
        out = assemble_full_sequences(
            row, contigs_dir=str(tmp_path), sample_name_from="grandparent",
            alpha_leader="from_contig", beta_leader=None, include_constant=False,
            verbose=False, show_progress=False, **kw,
        )
        return out.iloc[0]

    def test_germline_anchor_trims_overcapture(self, tmp_path):
        r = self._assemble(tmp_path, self.TRAV1_2_CONTIG, "TRAV1-2")
        assert r["alpha_leader_aa"] == "MWGVFLLYVSMKMGGTT"   # native 17-aa SP
        assert r["alpha_leader_qc"] == "ok"
        assert r["alpha_leader_source"] == "contig_germline_anchored"
        assert r["alpha_germline_allele"] == "TRAV1-2*01"
        assert r["alpha_germline_identity"] == 1.0
        assert r["alpha_leader_vs_germline"] == "identical"

    def test_off_reference_gene_falls_back_to_heuristic(self, tmp_path):
        # V gene not in the reference → no germline anchor; the heuristic runs
        # (over-capture with a weak-Kozak h-region tail → hregion_trimmed).
        r = self._assemble(tmp_path, self.TRAV1_2_CONTIG, "TRBVZZ-9")
        assert r["alpha_leader_source"] != "contig_germline_anchored"
        assert r["alpha_leader_aa"] == "MWGVFLLYVSMKMGGTT"  # heuristic recovers it too
        assert "alpha_germline_allele" not in r or pd.isna(r.get("alpha_germline_allele"))

    # In-range leader at exactly the germline length (no over-capture).
    TRAV1_2_INRANGE = (
        "GCAGCA"
        + "ATGTGGGGCGTGTTCCTGCTGTACGTGAGCATGAAGATGGGCGGCACCACC"
        + "GAATTCGAATTCGAA"
    )

    def test_in_range_leader_germline_anchored(self, tmp_path):
        r = self._assemble(tmp_path, self.TRAV1_2_INRANGE, "TRAV1-2")
        assert r["alpha_leader_aa"] == "MWGVFLLYVSMKMGGTT"
        assert r["alpha_leader_source"] == "contig_germline_anchored"
        assert r["alpha_leader_qc"] == "ok"

    def test_wrong_gene_divergent_falls_back_but_flags_divergence(self, tmp_path):
        # Right contig, WRONG V gene (TRBV13): germline ANCHOR fails (29-aa suffix
        # doesn't start with M) → heuristic picks the start. But the universal
        # germline COMPARISON still runs and flags the divergence — the shipped
        # 17-aa leader vs TRBV13's 29-aa germline → length_mismatch, low identity.
        r = self._assemble(tmp_path, self.TRAV1_2_CONTIG, "TRBV13")
        assert r["alpha_leader_source"] != "contig_germline_anchored"
        assert str(r["alpha_germline_allele"]).startswith("TRBV13*")
        assert "length_mismatch" in r["alpha_leader_vs_germline"]
        assert r["alpha_germline_identity"] < 0.8


class TestCentralRecorder:
    """Every SP path flows through _record_leader + shared QC (review ask)."""

    def test_curated_default_leader_is_qcd(self, tmp_path):
        # A non-from_contig default leader (CD28) now carries source + QC, routed
        # through the same recorder as the contig/germline paths.
        from tcrsift.assemble import assemble_full_sequences

        row = pd.DataFrame([{
            "CDR3ab": "EFEFE_CASS",
            "CDR3_alpha": "EFEFE", "CDR3_beta": "CASS",
            "VDJ_alpha_aa": "EFEFE", "VDJ_alpha_nt": "GAATTCGAATTCGAA",
            "VDJ_beta_aa": "CASS" + "G" * 40 + "VETA",
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1",
            "alpha_j_gene": "TRAJ45", "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
        }])
        out = assemble_full_sequences(
            row, alpha_leader="CD28", beta_leader=None, include_constant=False,
            verbose=False, show_progress=False,
        ).iloc[0]
        from tcrsift.assemble import DEFAULT_LEADERS

        assert out["alpha_leader_aa"] == DEFAULT_LEADERS["CD28"]["aa"]
        assert out["alpha_leader_source"] == "curated_default"
        assert out["alpha_leader_qc"] == "ok"

    def test_curated_default_compared_to_native_germline(self, tmp_path):
        # The universal comparison runs for a curated default too: with a known V
        # gene, the (synthetic) CD28 leader is compared to the native germline so
        # its divergence is visible — exactly the "look it up + check difference"
        # signal, not only for germline-anchored leaders.
        from tcrsift.assemble import assemble_full_sequences

        row = pd.DataFrame([{
            "CDR3ab": "EFEFE_CASS",
            "CDR3_alpha": "EFEFE", "CDR3_beta": "CASS",
            "VDJ_alpha_aa": "EFEFE", "VDJ_alpha_nt": "GAATTCGAATTCGAA",
            "VDJ_beta_aa": "CASS" + "G" * 40 + "VETA",
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1",
            "alpha_j_gene": "TRAJ45", "beta_j_gene": "TRBJ1-1",
            "alpha_v_gene": "TRAV1-2",
            "samples": "S1",
        }])
        out = assemble_full_sequences(
            row, alpha_leader="CD28", beta_leader=None, include_constant=False,
            verbose=False, show_progress=False,
        ).iloc[0]
        # Native germline surfaced; CD28 (18 aa) differs from the 17-aa native SP
        # — flagged as a length divergence (insertion relative to germline), low id.
        assert out["alpha_germline_allele"] == "TRAV1-2*01"
        assert out["alpha_germline_leader_aa"] == "MWGVFLLYVSMKMGGTT"
        assert out["alpha_leader_vs_germline"] != "identical"
        assert "insertion" in out["alpha_leader_vs_germline"]
        assert out["alpha_germline_identity"] < 0.5

    def test_germline_columns_always_present(self, tmp_path):
        # Off-reference gene → germline columns exist but are NaN (uniform schema).
        from tcrsift.assemble import assemble_full_sequences

        row = pd.DataFrame([{
            "CDR3ab": "EFEFE_CASS",
            "CDR3_alpha": "EFEFE", "CDR3_beta": "CASS",
            "VDJ_alpha_aa": "EFEFE", "VDJ_alpha_nt": "GAATTCGAATTCGAA",
            "VDJ_beta_aa": "CASS" + "G" * 40 + "VETA",
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1",
            "alpha_j_gene": "TRAJ45", "beta_j_gene": "TRBJ1-1",
            "alpha_v_gene": "TRBVZZ-9",  # not in reference
            "samples": "S1",
        }])
        out = assemble_full_sequences(
            row, alpha_leader="CD28", beta_leader=None, include_constant=False,
            verbose=False, show_progress=False,
        ).iloc[0]
        for col in (
            "alpha_germline_allele", "alpha_germline_leader_aa",
            "alpha_germline_identity", "alpha_leader_vs_germline",
        ):
            assert col in out.index  # column present
            assert pd.isna(out[col])  # but NaN (no reference match)


class TestDivergenceCharacterization:
    def test_internal_deletion(self):
        germ = "MLLLLLLLGPGISLLLPGSLAGSGL"  # TRBV20-1*01, 25 aa
        donor = "MLLLLLLLGPGSGL"            # shared N+C, internal gap
        d = characterize_divergence(donor, germ)
        assert d.startswith("internal_deletion:")
        assert "ISLLLPGSLA" in d  # the missing chunk is reported

    def test_5p_truncation(self):
        germ = "MABCDEFGHIJKLMNP"
        donor = germ[-8:]  # clean C-terminal suffix
        assert characterize_divergence(donor, germ).startswith("5p_truncation")

    def test_insertion(self):
        assert characterize_divergence("M" + "A" * 20, "M" + "A" * 17).startswith("insertion")

    def test_same_length_substitution(self):
        # equal length → position diff, not a length descriptor
        assert characterize_divergence("MWGVYLLYVSMKMGGTT", "MWGVFLLYVSMKMGGTT") == "F5Y"


class TestLeaderPolicy:
    """apply_leader_policy: keep SP-sound+consistent contig leaders; switch
    bad/inconsistent ones to germline; force / secondary-construct options (#270)."""

    def _germline(self):
        g_aa = dict((a, aa) for a, _f, aa in germline_vgene_leaders("TRBV20-1"))["01"]
        return g_aa, germline_leader_nt("TRBV20-1", "01")

    def _df(self, beta_leaders, *, sp_ok=True):
        # beta: TRBV20-1, divergent contig leaders; alpha: matches germline.
        from tcrsift.assemble import signal_peptide_features
        g_aa, _g_nt = self._germline()
        rows = []
        for i, lead in enumerate(beta_leaders):
            lead_nt = "ATG" + "CTC" * (len(lead) - 1)
            body_aa, body_nt = "CASSF" + "G" * 30, "TGC" * 35
            rows.append({
                "CDR3ab": f"CAA_{i}_CASSF",
                "alpha_v_gene": "TRAV1-2", "beta_v_gene": "TRBV20-1",
                "alpha_leader_aa": "MWGVFLLYVSMKMGGTT", "alpha_leader_source": "contig",
                "alpha_germline_leader_aa": "MWGVFLLYVSMKMGGTT",
                "alpha_germline_allele": "TRAV1-2*01", "alpha_germline_identity": 1.0,
                "alpha_leader_qc": "ok", "alpha_sp_features_ok": True,
                "full_alpha_aa": "MWGVFLLYVSMKMGGTT" + "CAVF" + "G" * 20,
                "beta_leader_aa": lead, "beta_leader_nt": lead_nt,
                "beta_leader_source": "contig", "beta_leader_qc": "ok",
                "beta_sp_features_ok": sp_ok if sp_ok is not True
                else signal_peptide_features(lead)["features_ok"],
                "beta_germline_leader_aa": g_aa, "beta_germline_allele": "TRBV20-1*01",
                "beta_germline_identity": 0.5,  # divergent from germline
                "full_beta_aa": lead + body_aa, "full_beta_nt": lead_nt + body_nt,
            })
        return pd.DataFrame(rows)

    def test_default_keeps_consistent_sound_divergent(self):
        # SP-sound + consistent-divergent (likely germline indel) → KEPT by default.
        from tcrsift.assemble import apply_leader_policy
        df = self._df(["MLLLLLLLGPGSGL", "MLLLLLLLGPGSGL"])
        out = apply_leader_policy(df, "T2A")
        assert len(out) == 2  # no twins, no switch
        assert (out["beta_leader_aa"] == "MLLLLLLLGPGSGL").all()

    def test_switches_inconsistent_divergent_to_germline(self):
        # divergent but NOT consistent across the gene → switch to germline.
        from tcrsift.assemble import apply_leader_policy
        g_aa, _ = self._germline()
        df = self._df(["MLLLLLLLGPGSGL", "MLLLLLLLGPAAAAGSGL"])  # 2 different
        out = apply_leader_policy(df, "T2A")
        assert (out["beta_leader_aa"] == g_aa).all()  # both switched to germline
        assert (out["beta_leader_source"] == "germline_reference_leader").all()
        assert out["full_beta_aa"].iloc[0].startswith(g_aa)

    def test_switches_bad_sp_to_germline(self):
        # bad SP (features_ok False), even if consistent → switch to germline.
        from tcrsift.assemble import apply_leader_policy
        g_aa, _ = self._germline()
        df = self._df(["MKKKKEEDDKKEE", "MKKKKEEDDKKEE"], sp_ok=False)  # no h-core
        out = apply_leader_policy(df, "T2A")
        assert (out["beta_leader_aa"] == g_aa).all()

    def test_secondary_germline_emits_twin(self):
        # --leader-secondary germline: keep primary, ALSO emit a germline twin.
        from tcrsift.assemble import apply_leader_policy
        g_aa, g_nt = self._germline()
        df = self._df(["MLLLLLLLGPGSGL", "MLLLLLLLGPGSGL"])
        out = apply_leader_policy(df, "T2A", secondary_beta="germline")
        assert len(out) == 4
        prim = out[out["leader_variant"] == "primary"]
        sec = out[out["leader_variant"] == "secondary:germline"]
        assert len(prim) == 2 and len(sec) == 2
        assert (prim["beta_leader_aa"] == "MLLLLLLLGPGSGL").all()  # primary kept
        for _, s in sec.iterrows():
            assert s["beta_leader_aa"] == g_aa
            assert s["full_beta_aa"].startswith(g_aa)
            assert g_aa in s["single_chain_aa"]

    def test_force_beta_leader(self):
        from tcrsift.assemble import apply_leader_policy
        g_aa, _ = self._germline()
        df = self._df(["MLLLLLLLGPGSGL", "MLLLLLLLGPGSGL"])
        out = apply_leader_policy(df, "T2A", force_beta="germline")
        assert (out["beta_leader_aa"] == g_aa).all()
        assert out["beta_leader_source"].iloc[0].startswith("forced_")
        # alpha untouched by a beta-only force
        assert (out["alpha_leader_aa"] == "MWGVFLLYVSMKMGGTT").all()

    def test_force_curated_all_chains(self):
        from tcrsift.assemble import DEFAULT_LEADERS, apply_leader_policy
        df = self._df(["MLLLLLLLGPGSGL"])
        out = apply_leader_policy(df, "T2A", force_alpha="CD8A", force_beta="CD8A")
        assert (out["beta_leader_aa"] == DEFAULT_LEADERS["CD8A"]["aa"]).all()
        assert (out["alpha_leader_aa"] == DEFAULT_LEADERS["CD8A"]["aa"]).all()
