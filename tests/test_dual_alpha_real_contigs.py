# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Focal tests on REAL dual-α data (B1-3, cell CGCCGTCGTATGGAGA-1).

This clone is the one that blocked the sarah-tcrseq B1-3 regen: a dual-α clone
whose non-dominant α variant fell to no_contig → blanket-N, even though BOTH α
are fully sequenced through the TRAC C-region. Sequences below are extracted
verbatim from B1-3/per_sample_outs/AIMpos-3/vdj_t/filtered_contig.fasta and the
phenotyped.h5ad obs (TRA_1 = CAVNDSGGGADGLTF/TRAJ45, TRA_2 = CAGQLGKTQGGSEKLVF/
TRAJ57). Both real contigs cover TRAC, so the correct per-α wiring must read the
J→C junction (Y) from each α's OWN contig.
"""

from __future__ import annotations

import pandas as pd

from tcrsift.assemble import assemble_full_sequences

# --- real extracted data (cell CGCCGTCGTATGGAGA-1, sample AIMpos-3) -----------
A2_CDR3, A2_J = "CAVNDSGGGADGLTF", "TRAJ45"
A2_VDJ_AA = "AQSVTQLGSHVSVSEGALVLLRCNYSSSVPPYLFWYVQYPNQGLQLLLKYTSAATLVKGINGFEAEFKKSETSFHLTKPSAHMSDAAEYFCAVNDSGGGADGLTFGKGTHLIIQP"
A2_VDJ_NT = "GCCCAGTCGGTGACCCAGCTTGGCAGCCACGTCTCTGTCTCTGAAGGAGCCCTGGTTCTGCTGAGGTGCAACTACTCATCGTCTGTTCCACCATATCTCTTCTGGTATGTGCAATACCCCAACCAAGGACTCCAGCTTCTCCTGAAGTACACATCAGCGGCCACCCTGGTTAAAGGCATCAACGGTTTTGAGGCTGAATTTAAGAAGAGTGAAACCTCCTTCCACCTGACGAAACCCTCAGCCCATATGAGCGACGCGGCTGAGTACTTCTGTGCTGTGAATGATTCAGGAGGAGGTGCTGACGGACTCACCTTTGGCAAAGGGACTCATCTAATCATCCAGCCCT"
A2_CONTIG_ID = "CGCCGTCGTATGGAGA-1_contig_1"
A2_CONTIG = "TTTTGAAACCCTTCAAAGGCAGAGACTTGTCCAGCCTAACCTGCCTGCTGCTCCTAGCTCCTGAGGCTCAGGGCCCTTGGCTTCTGTCCGCTCTGCTCAGGGCCCTCCAGCGTGGCCACTGCTCAGCCATGCTCCTGCTGCTCGTCCCAGTGCTCGAGGTGATTTTTACCCTGGGAGGAACCAGAGCCCAGTCGGTGACCCAGCTTGGCAGCCACGTCTCTGTCTCTGAAGGAGCCCTGGTTCTGCTGAGGTGCAACTACTCATCGTCTGTTCCACCATATCTCTTCTGGTATGTGCAATACCCCAACCAAGGACTCCAGCTTCTCCTGAAGTACACATCAGCGGCCACCCTGGTTAAAGGCATCAACGGTTTTGAGGCTGAATTTAAGAAGAGTGAAACCTCCTTCCACCTGACGAAACCCTCAGCCCATATGAGCGACGCGGCTGAGTACTTCTGTGCTGTGAATGATTCAGGAGGAGGTGCTGACGGACTCACCTTTGGCAAAGGGACTCATCTAATCATCCAGCCCTATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACT"

A1_CDR3, A1_J = "CAGQLGKTQGGSEKLVF", "TRAJ57"
A1_VDJ_AA = "GQQLNQSPQSMFIQEGEDVSMNCTSSSIFNTWLWYKQDPGEGPVLLIALYKAGELTSNGRLTAQFGITRKDSFLNISASIPSDVGIYFCAGQLGKTQGGSEKLVFGKGTKLTVNP"
A1_VDJ_NT = "GGTCAACAGCTGAATCAGAGTCCTCAATCTATGTTTATCCAGGAAGGAGAAGATGTCTCCATGAACTGCACTTCTTCAAGCATATTTAACACCTGGCTATGGTACAAGCAGGACCCTGGGGAAGGTCCTGTCCTCTTGATAGCCTTATATAAGGCTGGTGAATTGACCTCAAATGGAAGACTGACTGCTCAGTTTGGTATAACCAGAAAGGACAGCTTCCTGAATATCTCAGCATCCATACCTAGTGATGTAGGCATCTACTTCTGTGCTGGGCAGCTGGGGAAAACTCAGGGCGGATCTGAAAAGCTGGTCTTTGGAAAGGGAACGAAACTGACAGTAAACCCAT"
A1_CONTIG_ID = "CGCCGTCGTATGGAGA-1_contig_3"
A1_CONTIG = "GATTCAGGAAATAATTCTTTGCTGATAAGGATGCTCCTTGAACATTTATTAATAATCTTGTGGATGCAGCTGACATGGGTCAGTGGTCAACAGCTGAATCAGAGTCCTCAATCTATGTTTATCCAGGAAGGAGAAGATGTCTCCATGAACTGCACTTCTTCAAGCATATTTAACACCTGGCTATGGTACAAGCAGGACCCTGGGGAAGGTCCTGTCCTCTTGATAGCCTTATATAAGGCTGGTGAATTGACCTCAAATGGAAGACTGACTGCTCAGTTTGGTATAACCAGAAAGGACAGCTTCCTGAATATCTCAGCATCCATACCTAGTGATGTAGGCATCTACTTCTGTGCTGGGCAGCTGGGGAAAACTCAGGGCGGATCTGAAAAGCTGGTCTTTGGAAAGGGAACGAAACTGACAGTAAACCCATATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACT"

BETA = "CASSISASGGTEQFF"  # shared β of the real clone (assembled canonically here)


def _write_contigs(tmp_path, *contigs):
    """Write (contig_id, seq) pairs as a CellRanger-style fasta under sample S1."""
    d = tmp_path / "S1" / "vdj_t"
    d.mkdir(parents=True)
    (d / "filtered_contig.fasta").write_text(
        "".join(f">{cid}\n{seq}\n" for cid, seq in contigs)
    )
    return tmp_path


def _alpha_row(cdr3, vdj_aa, vdj_nt, j_gene, contig_id):
    """A single-α assembleable row (β is a stub; we assert only on α)."""
    return {
        "CDR3ab": f"{cdr3}_{BETA}",
        "CDR3_alpha": cdr3, "CDR3_beta": BETA,
        "VDJ_alpha_aa": vdj_aa, "VDJ_alpha_nt": vdj_nt,
        "VDJ_beta_aa": "CASS" + "G" * 40 + "VETA",
        "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1",
        "alpha_j_gene": j_gene, "beta_j_gene": "TRBJ1-1",
        "samples": "S1", "alpha_contig_ids": contig_id,
    }


def _assemble(tmp_path, row, contigs):
    # Lay contigs out CellRanger-style (<root>/S1/vdj_t/filtered_contig.fasta)
    # and load via cellranger_dir so the contig is keyed by sample "S1"
    # (matching the row's samples), mirroring the real run invocation.
    cdir = _write_contigs(tmp_path, *contigs)
    out = assemble_full_sequences(
        pd.DataFrame([row]), cellranger_dir=str(cdir), sample_name_from="parent",
        alpha_leader=None, beta_leader=None, verbose=False, show_progress=False,
    )
    return out.iloc[0]


class TestRealDualAlphaContigs:
    def test_alpha2_reads_junction_from_own_contig(self, tmp_path):
        r = _assemble(
            tmp_path,
            _alpha_row(A2_CDR3, A2_VDJ_AA, A2_VDJ_NT, A2_J, A2_CONTIG_ID),
            [(A2_CONTIG_ID, A2_CONTIG)],
        )
        assert r["alpha_junction_residue_source"] == "contig"
        assert r["alpha_junction_residue"] == "Y"

    def test_alpha1_reads_junction_from_own_contig(self, tmp_path):
        r = _assemble(
            tmp_path,
            _alpha_row(A1_CDR3, A1_VDJ_AA, A1_VDJ_NT, A1_J, A1_CONTIG_ID),
            [(A1_CONTIG_ID, A1_CONTIG)],
        )
        assert r["alpha_junction_residue_source"] == "contig"
        assert r["alpha_junction_residue"] == "Y"

    def test_the_bug_wrong_alpha_contig_falls_back(self, tmp_path):
        # The exact failure: the α1 variant pointed at α2's contig (the merged
        # clone's dominant-α list). α1's VDJ isn't in α2's contig → no_contig →
        # canonical fallback (blanket N). This is what the wiring fix prevents.
        r = _assemble(
            tmp_path,
            _alpha_row(A1_CDR3, A1_VDJ_AA, A1_VDJ_NT, A1_J, A2_CONTIG_ID),
            [(A2_CONTIG_ID, A2_CONTIG)],  # only the WRONG (α2) contig present
        )
        assert r["alpha_junction_residue_source"] == "canonical_fallback"
        assert r["alpha_junction_residue"] == "N"

    def test_case_insensitive_lowercase_vdj_nt(self, tmp_path):
        # Loader writes lowercase vdj_nt; contig is uppercase (#235 secondary).
        row = _alpha_row(A2_CDR3, A2_VDJ_AA, A2_VDJ_NT.lower(), A2_J, A2_CONTIG_ID)
        r = _assemble(tmp_path, row, [(A2_CONTIG_ID, A2_CONTIG)])
        assert r["alpha_junction_residue_source"] == "contig"
        assert r["alpha_junction_residue"] == "Y"

    def test_both_real_contigs_cover_trac(self):
        # Sanity on the raw data: each α's contig contains its VDJ AND continues
        # into the TRAC C-region (IQNPD start) — i.e. coverage was never missing.
        for vdj_nt, contig in ((A2_VDJ_NT, A2_CONTIG), (A1_VDJ_NT, A1_CONTIG)):
            i = contig.upper().find(vdj_nt.upper())
            assert i >= 0, "VDJ not found in its own contig"
            # canonical TRAC mature start IQNPD = ...ATCCAGAACCCTGAC...
            assert "ATCCAGAACCCTGAC" in contig.upper()[i:]
