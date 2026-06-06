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

"""Tests for the native selected-clones report (#188)."""

from __future__ import annotations

import matplotlib
import pandas as pd

matplotlib.use("Agg")

from tcrsift.assemble import back_translate  # noqa: E402
from tcrsift.report import build_selected_report, expand_dual_alpha_variants  # noqa: E402


def _obs_with_dual_alpha(a1, a2, beta):
    """A cell that co-expresses both alphas + the beta, plus a decoy."""
    return pd.DataFrame([
        {
            "TRA_1_cdr3": a1, "TRA_2_cdr3": a2, "TRB_1_cdr3": beta,
            "TRA_1_umis": 9, "TRA_2_umis": 7, "TRB_1_umis": 9,
            "TRA_1_vdj_aa": a1, "TRA_1_vdj_nt": back_translate(a1),
            "TRA_1_v_gene": "TRAV1", "TRA_1_j_gene": "TRAJ1", "TRA_1_c_gene": "TRAC",
            "TRA_2_vdj_aa": a2, "TRA_2_vdj_nt": back_translate(a2),
            "TRA_2_v_gene": "TRAV2", "TRA_2_j_gene": "TRAJ2", "TRA_2_c_gene": "TRAC",
        },
        {  # decoy: different beta, shouldn't be picked
            "TRA_1_cdr3": a1, "TRA_2_cdr3": a2, "TRB_1_cdr3": "CZZZF",
            "TRA_1_umis": 99, "TRA_2_umis": 99, "TRB_1_umis": 99,
            "TRA_1_vdj_aa": a1, "TRA_1_vdj_nt": back_translate(a1),
            "TRA_2_vdj_aa": a2, "TRA_2_vdj_nt": back_translate(a2),
        },
    ])


class TestExpandDualAlphaVariants:
    def test_emits_both_alpha_variants(self):
        a1, a2, beta = "CAAAF", "CBBBF", "CASSBETAF"
        selected = pd.DataFrame([{
            "CDR3ab": f"{a1}_{beta}", "CDR3_alpha": a1, "CDR3_beta": beta,
            "merged_alpha_partners": f"{a1};{a2}",
        }])
        obs = _obs_with_dual_alpha(a1, a2, beta)
        out, variant_of = expand_dual_alpha_variants(selected, obs)
        assert set(out["CDR3ab"]) == {f"{a1}_{beta}", f"{a2}_{beta}"}
        # Each variant's alpha VDJ matches its own CDR3.
        for _, r in out.iterrows():
            assert r["CDR3_alpha"] in str(r["VDJ_alpha_aa"])
        assert variant_of == {f"{a1}_{beta}": f"{a1}_{beta}", f"{a2}_{beta}": f"{a1}_{beta}"}
        assert set(out["selected_clone"]) == {f"{a1}_{beta}"}

    def test_non_merged_passthrough(self):
        selected = pd.DataFrame([{
            "CDR3ab": "CAAAF_CBETAF", "CDR3_alpha": "CAAAF", "CDR3_beta": "CBETAF",
            "merged_alpha_partners": "",
        }])
        out, variant_of = expand_dual_alpha_variants(selected, pd.DataFrame())
        assert list(out["CDR3ab"]) == ["CAAAF_CBETAF"]
        assert variant_of == {}

    def test_noop_without_obs(self):
        selected = pd.DataFrame([{
            "CDR3ab": "X_Y", "CDR3_alpha": "X", "CDR3_beta": "Y",
            "merged_alpha_partners": "A;B",
        }])
        out, variant_of = expand_dual_alpha_variants(selected, None)
        assert list(out["CDR3ab"]) == ["X_Y"] and variant_of == {}

    def test_no_coexpressing_cell_passthrough(self):
        # merged clone but no cell co-expresses both alphas with this beta.
        selected = pd.DataFrame([{
            "CDR3ab": "CAAAF_CBETAF", "CDR3_alpha": "CAAAF", "CDR3_beta": "CBETAF",
            "merged_alpha_partners": "CAAAF;CBBBF",
        }])
        obs = pd.DataFrame([{
            "TRA_1_cdr3": "CAAAF", "TRA_2_cdr3": "", "TRB_1_cdr3": "CBETAF",
            "TRA_1_umis": 5, "TRA_2_umis": 0, "TRB_1_umis": 5,
        }])
        out, variant_of = expand_dual_alpha_variants(selected, obs)
        assert list(out["CDR3ab"]) == ["CAAAF_CBETAF"] and variant_of == {}


def _assembleable_clone(cdr3ab, a, b, **extra):
    row = {
        "CDR3ab": cdr3ab,
        "CDR3_alpha": a, "CDR3_beta": b,
        "VDJ_alpha_aa": a, "VDJ_beta_aa": b,
        "VDJ_alpha_nt": back_translate(a), "VDJ_beta_nt": back_translate(b),
        "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1", "beta_j_gene": "TRBJ1-1",
        "samples": "S1",
    }
    row.update(extra)
    return row


class TestBuildSelectedReport:
    def test_writes_outputs_and_provenance(self, tmp_path):
        a = "CASS" + "A" * 40 + "VLF"
        b = "CASS" + "G" * 40 + "VEF"
        clonotypes = pd.DataFrame([_assembleable_clone("c1", a, b)])
        selected = pd.DataFrame([{"CDR3ab": "c1", "selection_rule": "shared", "global_rank": 1}])

        out = build_selected_report(
            selected, clonotypes, tmp_path,
            provenance_cols=["selection_rule", "global_rank"],
        )
        assert (tmp_path / "selected_clones.csv").exists()
        assert (tmp_path / "selected_clones_sequences.pdf").exists()
        assert (tmp_path / "selected_clones_qc.txt").exists()
        # Provenance carried through to the assembled output.
        written = pd.read_csv(tmp_path / "selected_clones.csv")
        assert "selection_rule" in written.columns
        assert "single_chain_aa" in out.columns

    def test_dual_alpha_emits_two_constructs(self, tmp_path):
        a1 = "CASS" + "A" * 40 + "VLF"
        a2 = "CASS" + "C" * 40 + "VLF"
        b = "CASS" + "G" * 40 + "VEF"
        clonotypes = pd.DataFrame([
            _assembleable_clone(f"{a1}_{b}", a1, b, merged_alpha_partners=f"{a1};{a2}"),
        ])
        selected = pd.DataFrame([{"CDR3ab": f"{a1}_{b}", "selection_rule": "private"}])
        obs = pd.DataFrame([{
            "TRA_1_cdr3": a1, "TRA_2_cdr3": a2, "TRB_1_cdr3": b,
            "TRA_1_umis": 9, "TRA_2_umis": 7, "TRB_1_umis": 9,
            "TRA_1_vdj_aa": a1, "TRA_1_vdj_nt": back_translate(a1),
            "TRA_1_v_gene": "TRAV1", "TRA_1_j_gene": "TRAJ1", "TRA_1_c_gene": "TRAC",
            "TRA_2_vdj_aa": a2, "TRA_2_vdj_nt": back_translate(a2),
            "TRA_2_v_gene": "TRAV2", "TRA_2_j_gene": "TRAJ2", "TRA_2_c_gene": "TRAC",
        }])
        out = build_selected_report(
            selected, clonotypes, tmp_path, obs=obs, provenance_cols=["selection_rule"],
        )
        # Two constructs (one per alpha variant) of the same selected clone.
        assert len(out) == 2
        assert set(out["CDR3ab"]) == {f"{a1}_{b}", f"{a2}_{b}"}


class TestCoverAndCombine:
    def _assembleable(self, cdr3ab, a, b):
        return _assembleable_clone(cdr3ab, a, b)

    def test_cover_page_prepended(self, tmp_path):
        from pypdf import PdfReader

        a = "CASS" + "A" * 40 + "VLF"
        b = "CASS" + "G" * 40 + "VEF"
        clonotypes = pd.DataFrame([_assembleable_clone("c1", a, b)])
        selected = pd.DataFrame([{"CDR3ab": "c1", "selection_rule": "shared"}])
        out = build_selected_report(
            selected, clonotypes, tmp_path, provenance_cols=["selection_rule"], cover=True,
        )
        assert len(out) == 1
        pdf = PdfReader(str(tmp_path / "selected_clones_sequences.pdf"))
        # cover page + 1 construct page.
        assert len(pdf.pages) == 2
        assert "T2A" in pdf.pages[0].extract_text()

    def test_no_cover(self, tmp_path):
        from pypdf import PdfReader

        a = "CASS" + "A" * 40 + "VLF"
        b = "CASS" + "G" * 40 + "VEF"
        clonotypes = pd.DataFrame([_assembleable_clone("c1", a, b)])
        selected = pd.DataFrame([{"CDR3ab": "c1"}])
        build_selected_report(selected, clonotypes, tmp_path, cover=False)
        pdf = PdfReader(str(tmp_path / "selected_clones_sequences.pdf"))
        assert len(pdf.pages) == 1  # just the construct page

    def test_combine_selected_pdfs(self, tmp_path):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from pypdf import PdfReader

        from tcrsift.report import combine_selected_pdfs

        # Two tiny per-donor "selected" PDFs (2 pages each).
        donor_pdfs = []
        for name in ("B1-2", "B1-3"):
            d = tmp_path / name
            d.mkdir()
            p = d / "selected_clones_sequences.pdf"
            from matplotlib.backends.backend_pdf import PdfPages
            with PdfPages(p) as pp:
                for _ in range(2):
                    fig = plt.figure()
                    pp.savefig(fig)
                    plt.close(fig)
            donor_pdfs.append(p)
        out = tmp_path / "cohort.pdf"
        combine_selected_pdfs(donor_pdfs, out, labels=["B1-2", "B1-3"])
        pages = PdfReader(str(out)).pages
        # cover + legend + (title + 2 pages) x 2 = 2 + 6 = 8.
        assert len(pages) == 8
