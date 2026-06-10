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

"""Tests for the native independent QC battery (#279)."""

from __future__ import annotations

import pandas as pd

from tcrsift.qc_battery import (
    EXPECTED_T2A,
    check_assembly,
    check_cdr3_ref,
    check_contigs,
    check_source,
    check_synth,
    contig_cdr3_map,
    qc_summary_markdown,
    run_qc_battery,
)

# Independent codon table (same as the module) to back-translate test AA into NT
# so the round-trip check passes for a hand-built construct.
_B = "TCAG"
_COD = [a + b + c for a in _B for b in _B for c in _B]
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_AA2COD = {}
for cod, aa in zip(_COD, _AA):
    _AA2COD.setdefault(aa, cod)


def _bt(aa: str) -> str:
    """Back-translate AA → NT (one codon per residue)."""
    return "".join(_AA2COD[a] for a in aa)


def _good_row():
    """A single, internally-consistent β-T2A-α construct that passes A–E."""
    a_lead, a_vdj, a_const = "MALPV", "CAVRDF", "IQNP"
    b_lead, b_vdj, b_const = "MGTRL", "CASSLF", "EDLN"
    full_a = a_lead + a_vdj + a_const
    full_b = b_lead + b_vdj + b_const
    sc = full_b + EXPECTED_T2A + full_a
    return {
        "CDR3ab": "CASSLF_CAVRDF",
        "selected_clone": "clone1",
        "cell_count": 10.0,
        # leaders / VDJ / constants — AA + matching NT
        "alpha_leader_aa": a_lead, "alpha_leader_nt": _bt(a_lead),
        "beta_leader_aa": b_lead, "beta_leader_nt": _bt(b_lead),
        "VDJ_alpha_aa": a_vdj, "VDJ_alpha_nt": _bt(a_vdj),
        "VDJ_beta_aa": b_vdj, "VDJ_beta_nt": _bt(b_vdj),
        "alpha_constant_aa": a_const, "beta_constant_aa": b_const,
        "full_alpha_aa": full_a, "full_alpha_nt": _bt(full_a),
        "full_beta_aa": full_b, "full_beta_nt": _bt(full_b),
        "single_chain_aa": sc, "single_chain_nt": _bt(sc),
        # CDR3 + gene names
        "CDR3_alpha": a_vdj, "CDR3_beta": b_vdj,
        "alpha_v_gene": "TRAV1-2", "alpha_j_gene": "TRAJ33", "alpha_c_gene": "TRAC",
        "beta_v_gene": "TRBV20-1", "beta_j_gene": "TRBJ2-1", "beta_c_gene": "TRBC2",
        "alpha_junction_residue": "C", "beta_junction_residue": "C",
        # synthesis flags
        "construct_contig_verified": True,
        "synth_duplicate_construct": False,
        "synth_alpha_beta_swap": False,
        "alpha_contig_ids": "a1", "beta_contig_ids": "b1",
    }


def _good_df():
    return pd.DataFrame([_good_row()])


def test_all_checks_pass_on_consistent_construct():
    d = _good_df()
    cmap = {"a1": "CAVRDF", "b1": "CASSLF"}
    ct = d[["CDR3ab", "alpha_v_gene", "alpha_j_gene", "beta_v_gene",
            "beta_j_gene", "CDR3_alpha", "CDR3_beta", "cell_count"]]
    results = run_qc_battery(d, clonotypes=ct, contig_map=cmap)
    assert [r.status for r in results] == ["PASS"] * 5


def test_assembly_fails_on_recomposition_break():
    d = _good_df()
    d.loc[0, "full_alpha_aa"] = "MWRONGSEQ"  # no longer leader+vdj+const
    assert check_assembly(d).status == "FAIL"


def test_assembly_fails_on_wrong_linker():
    d = _good_df()
    row = d.loc[0]
    sc = row["full_beta_aa"] + "GGSGGSGG" + row["full_alpha_aa"]  # valid AA, wrong linker
    d.loc[0, "single_chain_aa"] = sc
    d.loc[0, "single_chain_nt"] = _bt(sc)
    assert check_assembly(d).status == "FAIL"


def test_cdr3_traj35_c_ending_allowed():
    d = _good_df()
    # CDR3 ending in C with a TRAJ35 call is non-canonical-but-correct.
    d.loc[0, "CDR3_alpha"] = "CAVRDC"
    d.loc[0, "VDJ_alpha_aa"] = "MALPVCAVRDCIQNP"
    d.loc[0, "alpha_j_gene"] = "TRAJ35"
    res = check_cdr3_ref(d)
    assert res.status == "PASS"
    assert "non-canonical" in res.detail


def test_cdr3_c_ending_without_traj35_fails():
    d = _good_df()
    d.loc[0, "CDR3_alpha"] = "CAVRDC"
    d.loc[0, "VDJ_alpha_aa"] = "MALPVCAVRDCIQNP"
    d.loc[0, "alpha_j_gene"] = "TRAJ33"  # canonical J, C-ending is wrong
    assert check_cdr3_ref(d).status == "FAIL"


def test_cdr3_invalid_gene_name_fails():
    d = _good_df()
    d.loc[0, "beta_v_gene"] = "IGHV1"  # not a TRBV
    assert check_cdr3_ref(d).status == "FAIL"


def test_source_skip_without_clonotypes():
    assert check_source(_good_df(), None).status == "SKIP"


def test_source_fails_on_gene_mismatch():
    d = _good_df()
    ct = d[["CDR3ab", "beta_v_gene", "cell_count"]].copy()
    ct.loc[0, "beta_v_gene"] = "TRBV9"  # disagrees with the selected row
    assert check_source(d, ct).status == "FAIL"


def test_source_dual_alpha_secondary_not_keyed_is_ok():
    d = _good_df()
    ct = d[["CDR3ab", "beta_v_gene", "cell_count"]].iloc[0:0]  # empty index
    res = check_source(d, ct)
    assert res.status == "PASS"
    assert "dual-alpha secondaries" in res.detail


def test_synth_fails_on_duplicate_single_chain():
    d = pd.concat([_good_df(), _good_df()], ignore_index=True)  # identical SC
    assert check_synth(d).status == "FAIL"


def test_contigs_skip_without_map():
    assert check_contigs(_good_df(), None).status == "SKIP"
    assert check_contigs(_good_df(), {}).status == "SKIP"


def test_contigs_fails_on_orphan_id():
    d = _good_df()
    assert check_contigs(d, {"b1": "CASSLF"}).status == "FAIL"  # a1 missing


def test_contigs_fails_on_cdr3_not_in_any_contig():
    d = _good_df()
    cmap = {"a1": "TOTALLYDIFFERENT", "b1": "CASSLF"}
    assert check_contigs(d, cmap).status == "FAIL"


def test_summary_markdown_all_pass():
    d = _good_df()
    cmap = {"a1": "CAVRDF", "b1": "CASSLF"}
    ct = d[["CDR3ab", "alpha_v_gene", "alpha_j_gene", "beta_v_gene",
            "beta_j_gene", "CDR3_alpha", "CDR3_beta", "cell_count"]]
    results = run_qc_battery(d, clonotypes=ct, contig_map=cmap)
    md, ok = qc_summary_markdown({"B1-2": results}, deliverable="deliv")
    assert ok
    assert "ALL PASS" in md
    assert "B1-2" in md
    assert "TRAJ35" in md  # edge-case note present


def test_summary_markdown_skip_does_not_fail():
    d = _good_df()
    results = run_qc_battery(d)  # no clonotypes, no contigs -> C & E SKIP
    md, ok = qc_summary_markdown({"B1-2": results})
    assert ok  # SKIP must not flip overall to FAIL
    assert "SKIP" in md


def test_contig_cdr3_map_reads_annotations(tmp_path):
    d = tmp_path / "per_sample_outs" / "S1" / "vdj_t"
    d.mkdir(parents=True)
    pd.DataFrame({"contig_id": ["c1", "c2"], "cdr3": ["CASSF", "CAVRF"],
                  "extra": [1, 2]}).to_csv(
        d / "filtered_contig_annotations.csv", index=False)
    m = contig_cdr3_map(cellranger_dir=str(tmp_path / "per_sample_outs"))
    assert m == {"c1": "CASSF", "c2": "CAVRF"}


def test_contig_cdr3_map_empty_without_dir():
    assert contig_cdr3_map() == {}
