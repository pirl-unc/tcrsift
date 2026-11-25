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

"""Regression tests for v3.37.4 assembly and pipeline correctness fixes."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from tcrsift.assemble import (
    apply_leader_policy,
    assemble_full_sequences,
    build_assembly_qc_report,
    enforce_contig_fidelity,
    translate_dna,
    validate_sequences,
)
from tcrsift.assemble.vregion import germline_compare_vregion
from tcrsift.pipeline import _common_contig_root


def _sample_clone_df(*, with_nt: bool = True) -> pd.DataFrame:
    row = {
        "CDR3ab": "CAVRDSNYQLIW_CASSLGGNTEVFF",
        "CDR3_alpha": "CAVRDSNYQLIW",
        "CDR3_beta": "CASSLGGNTEVFF",
        "VDJ_alpha_aa": "AQKVTQAQTEISVVEKEDVTLDCVYETRDTTYYLFWYKQPPSGELVFLIRRNSFDEQNEISGRYSWNFQKSTSSFNFTITASQVVDSAVYFCAVRDSNYQLIWGAGTKLIIKP",
        "VDJ_beta_aa": "DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDDSGMPEDRFSAKMPNASFSTLKIQPSEPRDSAVYFCASSLGGNTEVFFGKGTRLTVV",
        "alpha_v_gene": "TRAV12-2",
        "alpha_j_gene": "TRAJ33",
        "alpha_c_gene": "TRAC",
        "beta_v_gene": "TRBV19",
        "beta_j_gene": "TRBJ1-1",
        "beta_c_gene": "TRBC1",
    }
    if with_nt:
        from tcrsift.assemble import optimize_codons

        row["VDJ_alpha_nt"] = optimize_codons(row["VDJ_alpha_aa"])
        row["VDJ_beta_nt"] = optimize_codons(row["VDJ_beta_aa"])
    return pd.DataFrame([row])


def test_aa_only_vdj_does_not_drop_vdj_from_full_nt_or_single_chain_nt():
    df = _sample_clone_df(with_nt=False)
    out = assemble_full_sequences(
        df,
        include_constant=True,
        alpha_leader="CD28",
        beta_leader="CD8A",
        linker="T2A",
        show_progress=False,
    )
    full_a_aa = out.loc[0, "full_alpha_aa"]
    full_a_nt = out.loc[0, "full_alpha_nt"]
    sc_aa = out.loc[0, "single_chain_aa"]
    sc_nt = out.loc[0, "single_chain_nt"]

    tr_a, _ = translate_dna(full_a_nt)
    tr_sc, _ = translate_dna(sc_nt)
    assert tr_a.rstrip("*") == full_a_aa.rstrip("*")
    assert tr_sc.rstrip("*") == sc_aa.rstrip("*")
    assert "CAVRDSNYQLIW" in tr_a
    assert "CASSLGGNTEVFF" in tr_sc


def test_p2a_and_f2a_linker_recoded_to_avoid_bsmbi():
    df = _sample_clone_df(with_nt=True)
    for linker_name in ("P2A", "F2A"):
        out = assemble_full_sequences(
            df,
            include_constant=True,
            alpha_leader="CD28",
            beta_leader="CD8A",
            linker=linker_name,
            avoid_enzymes=("BsmBI",),
            codon_optimization="focal",
            show_progress=False,
        )
        sc_opt = str(out.loc[0, "single_chain_nt_optimized"]).upper()
        assert "CGTCTC" not in sc_opt
        assert "GAGACG" not in sc_opt
        tr_opt, _ = translate_dna(sc_opt)
        assert tr_opt.rstrip("*") == out.loc[0, "single_chain_aa"].rstrip("*")


def test_include_constant_false_passes_validate_and_qc_report():
    df = _sample_clone_df(with_nt=True)
    out = assemble_full_sequences(
        df,
        include_constant=False,
        alpha_leader="CD28",
        beta_leader="CD8A",
        linker=None,
        show_progress=False,
    )
    lb = [w for w in validate_sequences(out) if w.severity == "load_bearing"]
    assert lb == []
    report = build_assembly_qc_report(out)
    assert report.passed is True


def test_from_data_constant_source_passes_fidelity_and_skips_unpaired_chain():
    df = _sample_clone_df(with_nt=True)
    df.loc[0, "VDJ_alpha_aa"] = None
    df.loc[0, "VDJ_alpha_nt"] = None
    df.loc[0, "alpha_constant_aa"] = "IQNPDPAVYQL"
    df.loc[0, "beta_constant_aa"] = "DLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFFPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSVSYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDF"
    out = assemble_full_sequences(
        df,
        include_constant=True,
        constant_source="from-data",
        alpha_leader=None,
        beta_leader="CD8A",
        linker=None,
        show_progress=False,
    )
    assert pd.isna(out.loc[0, "alpha_constant_source"]) if "alpha_constant_source" in out.columns else True
    assert pd.isna(out.loc[0, "alpha_c_gene_canonical"]) if "alpha_c_gene_canonical" in out.columns else True
    note = enforce_contig_fidelity(out, allow_canonical_fallback=True, context="test")
    assert "from-data" in note
    assert "germline-inferred" not in note


def test_germline_compare_vregion_flags_truncated_germline_length_mismatch(monkeypatch):
    from tcrsift.assemble import vregion

    # Donor framework is 95 aa; best germline allele is 92 aa matching all 92 prefix residues.
    donor_fr = "A" * 95
    monkeypatch.setattr(
        vregion,
        "germline_vregion",
        lambda *a, **kw: [("01", "F", "A" * 92)],
    )
    res = germline_compare_vregion(donor_fr + "CASSF", "CASSF", "TRBV19")
    assert res is not None
    _allele, _g_overlap, identity, diff = res
    assert identity < 1.0
    assert diff != "identical"
    assert "len:95vs92" in diff


def test_apply_leader_policy_handles_duplicate_index_and_all_nan_float_cols():
    df = _sample_clone_df(with_nt=True)
    out = assemble_full_sequences(
        df,
        include_constant=True,
        alpha_leader="CD28",
        beta_leader="CD8A",
        linker="T2A",
        show_progress=False,
    )
    dup = pd.concat([out, out], ignore_index=False)  # duplicate index [0, 0]
    dup["alpha_leader_source"] = "missing_from_contig"
    dup["alpha_leader_aa"] = float("nan")  # all-NaN float64 column (CSV reload)
    dup["alpha_leader_nt"] = float("nan")
    res = apply_leader_policy(dup, "T2A", leader_fallback="CD28", species="human")
    assert len(res) == 2
    assert list(res.index) == [0, 0]
    assert all(isinstance(v, str) and v.startswith("M") for v in res["alpha_leader_aa"])


def test_mouse_germline_leader_respects_custom_data_dir(tmp_path):
    # Write a minimal mouse IMGT leader cache under tmp_path/imgt_trav_leader_mouse/TRAV_LEADER_mouse.html
    imgt_dir = tmp_path / "imgt_trav_leader_mouse"
    imgt_dir.mkdir(parents=True)
    fasta = (
        "<pre>\n"
        ">X00001|TRAV1*01|Mus musculus|F|L-PART1+L-PART2|1..54|54 nt|\n"
        "atgctgcgcctgctgctggccctgaacctgttccccagcatccaggtgaccggc\n"
        "</pre>\n"
    )
    (imgt_dir / "TRAV_LEADER_mouse.html").write_text(fasta)

    df = _sample_clone_df(with_nt=True)
    df.loc[0, "alpha_v_gene"] = "TRAV1"
    df.loc[0, "alpha_c_gene"] = "Trac"
    df.loc[0, "beta_c_gene"] = "Trbc1"
    out = assemble_full_sequences(
        df,
        species="mouse",
        include_constant=True,
        alpha_leader="mCD28",
        beta_leader="mCD8B",
        force_alpha_leader="germline",
        linker="T2A",
        data_dir=tmp_path,
        show_progress=False,
    )
    assert out.loc[0, "alpha_leader_source"] == "forced_germline_reference_leader"
    assert out.loc[0, "alpha_germline_allele"] == "TRAV1*01"


def test_common_contig_root_resolves_relative_vdj_dirs_against_base_dir(tmp_path):
    sheet_dir = tmp_path / "study"
    s1 = sheet_dir / "samples" / "S1" / "vdj"
    s2 = sheet_dir / "samples" / "S2" / "vdj"
    s1.mkdir(parents=True)
    s2.mkdir(parents=True)
    root = _common_contig_root(["samples/S1/vdj", "samples/S2/vdj"], base_dir=sheet_dir)
    assert root is not None
    assert Path(root) == (sheet_dir / "samples").resolve()
