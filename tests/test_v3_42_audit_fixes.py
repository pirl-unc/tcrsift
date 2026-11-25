"""Regression tests for v3.42.0 cross-tool audit fixes.

Covers:
1. Codon-aligned germline V-gene leaders (fixes 67-nt TRAV41*01 frameshift).
2. Rejection of unspliced-intron / non-canonical J->C junction codons.
3. Cell Ranger 3 TRAJ37 over-extended VDJ_alpha_aa ending in D + bare TRAC seam.
4. 1-2 nt non-codon contig stubs past VDJ not emitting false load-bearing errors.
5. Bundled murine k-mer Pgen background models (`kmer_pgen_mouse_{alpha,beta}.npz`).
"""

from __future__ import annotations

import pandas as pd

from tcrsift import seqprob
from tcrsift.assemble import (
    HUMAN_TRAC_AA,
    HUMAN_TRBC1_AA,
    assemble_full_sequences,
    back_translate,
    translate_dna,
    validate_sequences,
)
from tcrsift.assemble.leaders import _load_human_reference, germline_leader_nt


def test_all_packaged_human_leaders_are_codon_aligned():
    """Every leader in vgene_leaders.csv.gz (including TRAV41*01) must be 3*len(aa)."""
    by_gene, nt_idx = _load_human_reference()
    assert by_gene, "Expected bundled human leader reference"
    for gene, recs in by_gene.items():
        for allele, _func, leader_aa in recs:
            leader_nt = nt_idx.get((gene, allele))
            assert leader_nt is not None, f"Missing NT for {gene}*{allele}"
            assert len(leader_nt) == 3 * len(leader_aa), (
                f"{gene}*{allele}: len(leader_nt)={len(leader_nt)} != 3*{len(leader_aa)}"
            )
            tr, ragged = translate_dna(leader_nt.upper())
            assert not ragged
            assert tr == leader_aa, f"{gene}*{allele}: {tr!r} != {leader_aa!r}"

    trav41_nt = germline_leader_nt("TRAV41", "01", species="human")
    assert trav41_nt is not None
    assert len(trav41_nt) == 66


def test_trav41_germline_leader_no_frameshift(tmp_path):
    """Assembling a TRAV41 clone with germline leader fallback must stay in-frame."""
    vdj_a = "KNEVEQSPQNLTAQEGEFITINCSYSVGISALHWLQQHPGGGIVSLFMLSSGKKKHGRLIATINIQEKHSSLHITASHPRDSAVYICAVRLQRGGGADGLTFGKGTHLIIQP"
    vdj_b = "DGGITQSPKYLFRKEGQNVTLSCEQNLNHDAMYWYRQDPGQGLRLIYYSQIVNDFQKGDIAEGYSVSREKKESFPLTVTSAQKNPTAFYLCASSIRSSYEQYFGPGTRLTVT"
    df = pd.DataFrame(
        [
            {
                "CDR3ab": "CAVRLQRGGGADGLTF_CASSIRSSYEQYF",
                "CDR3_alpha": "CAVRLQRGGGADGLTF",
                "CDR3_beta": "CASSIRSSYEQYF",
                "alpha_v_gene": "TRAV41",
                "alpha_j_gene": "TRAJ45",
                "alpha_c_gene": "TRAC",
                "beta_v_gene": "TRBV19",
                "beta_j_gene": "TRBJ2-7",
                "beta_c_gene": "TRBC2",
                "VDJ_alpha_aa": vdj_a,
                "VDJ_alpha_nt": back_translate(vdj_a).lower(),
                "VDJ_beta_aa": vdj_b,
                "VDJ_beta_nt": back_translate(vdj_b).lower(),
                "samples": "S1",
                "alpha_contig_ids": "missing_a",
                "beta_contig_ids": "missing_b",
            }
        ]
    )
    (tmp_path / "S1").mkdir()
    (tmp_path / "S1" / "filtered_contig.fasta").write_text(">dummy\nATGAAA\n")
    out = assemble_full_sequences(
        df,
        contigs_dir=str(tmp_path),
        alpha_leader="from_contig",
        beta_leader="from_contig",
        leader_fallback="germline",
        species="human",
        linker="T2A",
        verbose=False,
        show_progress=False,
    )
    assert out["alpha_leader_source"].iloc[0] == "germline_reference_leader"
    assert len(out["alpha_leader_nt"].iloc[0]) == 66
    lb = [m for m in validate_sequences(out, species="human") if m.severity == "load_bearing"]
    assert lb == [], f"Unexpected load-bearing errors: {lb}"


def test_unspliced_intron_and_traj37_overextended_seam_and_1nt_stub(tmp_path):
    """Unspliced intron (G), TRAJ37 over-extended D+bare TRAC, and 1-nt stubs all validate cleanly."""
    vdj_b = "DGGITQSPKYLFRKEGQNVTLSCEQNLNHDAMYWYRQDPGQGLRLIYYSQIVNDFQKGDIAEGYSVSREKKESFPLTVTSAQKNPTAFYLCASSIRSSYEQYFGPGTRLTVT"
    vdj_b_nt = back_translate(vdj_b)

    # Row 0: TRAJ44 followed by unspliced intron GGTAGGTAAC... (would have read 'G' -> GIQNPDPAV)
    vdj_a_intron = "QKEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMFIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVNTGTASKLTFGTGTRLQVTL"
    vdj_a_intron_nt = back_translate(vdj_a_intron)

    # Row 1: TRAJ37 where VDJ_alpha_aa already ends in 'D' and contig continues directly with bare TRAC (IQNPDP...)
    vdj_a_j37 = "QKEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMFIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVSGGSNYKLTFGQGTTLQVKPD"
    vdj_a_j37_nt = back_translate(vdj_a_j37)

    # Row 2: 1-nt stub ('A') past VDJ on alpha and 2-nt stub ('GA') on beta
    vdj_a_stub = "QKEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMFIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVIQGAQKLVFGQGTRLTINP"
    vdj_a_stub_nt = back_translate(vdj_a_stub)

    sdir = tmp_path / "S1"
    sdir.mkdir()
    (sdir / "filtered_contig.fasta").write_text(
        f">a_intron\n{vdj_a_intron_nt}GGTAGGTAACAGAAACCG\n"
        f">a_j37\n{vdj_a_j37_nt}{back_translate(HUMAN_TRAC_AA[:10])}\n"
        f">a_stub\n{vdj_a_stub_nt}A\n"
        f">b_ok\n{vdj_b_nt}{back_translate('E' + HUMAN_TRBC1_AA[:10])}\n"
        f">b_stub\n{vdj_b_nt}GA\n"
    )

    df = pd.DataFrame(
        [
            {
                "CDR3ab": "CAVNTGTASKLTF_CASSIRSSYEQYF",
                "CDR3_alpha": "CAVNTGTASKLTF",
                "CDR3_beta": "CASSIRSSYEQYF",
                "alpha_v_gene": "TRAV12-2",
                "alpha_j_gene": "TRAJ44",
                "alpha_c_gene": "TRAC",
                "beta_v_gene": "TRBV19",
                "beta_j_gene": "TRBJ1-1",
                "beta_c_gene": "TRBC1",
                "VDJ_alpha_aa": vdj_a_intron,
                "VDJ_alpha_nt": vdj_a_intron_nt.lower(),
                "VDJ_beta_aa": vdj_b,
                "VDJ_beta_nt": vdj_b_nt.lower(),
                "samples": "S1",
                "alpha_contig_ids": "a_intron",
                "beta_contig_ids": "b_ok",
            },
            {
                "CDR3ab": "CAVSGGSNYKLTF_CASSIRSSYEQYF",
                "CDR3_alpha": "CAVSGGSNYKLTF",
                "CDR3_beta": "CASSIRSSYEQYF",
                "alpha_v_gene": "TRAV12-2",
                "alpha_j_gene": "TRAJ37",
                "alpha_c_gene": "TRAC",
                "beta_v_gene": "TRBV19",
                "beta_j_gene": "TRBJ1-1",
                "beta_c_gene": "TRBC1",
                "VDJ_alpha_aa": vdj_a_j37,
                "VDJ_alpha_nt": vdj_a_j37_nt.lower(),
                "VDJ_beta_aa": vdj_b,
                "VDJ_beta_nt": vdj_b_nt.lower(),
                "samples": "S1",
                "alpha_contig_ids": "a_j37",
                "beta_contig_ids": "b_ok",
            },
            {
                "CDR3ab": "CAVIQGAQKLVF_CASSIRSSYEQYF",
                "CDR3_alpha": "CAVIQGAQKLVF",
                "CDR3_beta": "CASSIRSSYEQYF",
                "alpha_v_gene": "TRAV12-2",
                "alpha_j_gene": "TRAJ54",
                "alpha_c_gene": "TRAC",
                "beta_v_gene": "TRBV19",
                "beta_j_gene": "TRBJ1-1",
                "beta_c_gene": "TRBC1",
                "VDJ_alpha_aa": vdj_a_stub,
                "VDJ_alpha_nt": vdj_a_stub_nt.lower(),
                "VDJ_beta_aa": vdj_b,
                "VDJ_beta_nt": vdj_b_nt.lower(),
                "samples": "S1",
                "alpha_contig_ids": "a_stub",
                "beta_contig_ids": "b_stub",
            },
        ]
    )

    out = assemble_full_sequences(
        df,
        contigs_dir=str(tmp_path),
        species="human",
        linker="T2A",
        verbose=False,
        show_progress=False,
    )
    # Row 0: intron rejected -> germline TRAJ44 'D' + IQNPDPAV
    assert out["alpha_constant_aa"].iloc[0].startswith("DIQNPDPAV")
    assert out["alpha_junction_residue_source"].iloc[0] == "j_inferred"

    # Row 1: TRAJ37 D shifted cleanly to constant seam (single D + single I, no IIQNPDPAV or DDIQNPDPAV)
    assert out["vdj_alpha_aa"].iloc[1].endswith("LQVKP")
    assert out["alpha_constant_aa"].iloc[1].startswith("DIQNPDPAV")
    assert "LQVKPDIIQNP" not in out["full_alpha_aa"].iloc[1]
    assert "LQVKPDIQNPDPAV" in out["full_alpha_aa"].iloc[1]
    assert out["alpha_junction_residue_source"].iloc[1] == "contig"

    # Row 2: 1-2 nt stubs stored as None rather than flagged as out-of-frame
    assert out["alpha_constant_nt_contig"].iloc[2] is None
    assert out["beta_constant_nt_contig"].iloc[2] is None

    lb = [m for m in validate_sequences(out, species="human") if m.severity == "load_bearing"]
    assert lb == [], f"Unexpected load-bearing errors: {lb}"


def test_packaged_mouse_kmer_pgen_models_load_and_score_mouse_motifs():
    """Bundled mouse k-mer Pgen tables load with strict_species=True and score murine V/J motifs."""
    ma = seqprob.load_background_model("alpha", "kmer", "pgen", species="mouse", strict_species=True)
    mb = seqprob.load_background_model("beta", "kmer", "pgen", species="mouse", strict_species=True)
    hb = seqprob.load_background_model("beta", "kmer", "pgen", species="human", strict_species=True)
    assert ma.n_train == 500_000
    assert mb.n_train == 500_000

    # Mouse-specific Trbv3 (CGAR) + Trbj1-5 (NNQAPLF) / Trbj2-4 (QNTLYF) motifs
    mouse_cdr3s = ["CGARNTYNNQAPLF", "CASSPPGYNNQAPLF", "CGARDWGGQNTLYF"]
    m_scores = mb.log_prob(mouse_cdr3s)
    h_scores = hb.log_prob(mouse_cdr3s)
    assert all(m > h + 10.0 for m, h in zip(m_scores, h_scores)), (
        f"Expected mouse Pgen model to score murine motifs much higher than human: {m_scores} vs {h_scores}"
    )
