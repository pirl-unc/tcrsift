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
    germline_anchor_leader,
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
