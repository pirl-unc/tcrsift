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
            "TRA_1_contig_id": "cellA_contig_1",
            "TRA_2_vdj_aa": a2, "TRA_2_vdj_nt": back_translate(a2),
            "TRA_2_v_gene": "TRAV2", "TRA_2_j_gene": "TRAJ2", "TRA_2_c_gene": "TRAC",
            "TRA_2_contig_id": "cellA_contig_2",
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


class TestDualAlphaContigWiring:
    """Each dual-α variant must point at its OWN α's contig (the bug behind a
    dual-α clone falling to no_contig/blanket-N even though the second α was
    fully sequenced). The merged clone's alpha_contig_ids only carries the
    dominant α's contigs; the expansion must refresh it per variant from the
    rep cell's matching slot. Hammer the edge cases."""

    def _clone(self, a1, a2, beta, **extra):
        # The merged clone carries ONLY the dominant α's contig list — the value
        # each variant must NOT keep.
        d = {
            "CDR3ab": f"{a1}_{beta}", "CDR3_alpha": a1, "CDR3_beta": beta,
            "merged_alpha_partners": f"{a1};{a2}",
            "alpha_contig_ids": "DOMINANT_only_c1;DOMINANT_only_c2",
        }
        d.update(extra)
        return pd.DataFrame([d])

    def test_each_variant_gets_its_own_slot_contig(self):
        a1, a2, beta = "CAAAF", "CBBBF", "CASSBETAF"
        obs = _obs_with_dual_alpha(a1, a2, beta)  # a1->TRA_1(c1), a2->TRA_2(c2)
        out, _ = expand_dual_alpha_variants(self._clone(a1, a2, beta), obs)
        cid = dict(zip(out["CDR3ab"], out["alpha_contig_ids"]))
        assert cid[f"{a1}_{beta}"] == "cellA_contig_1"
        assert cid[f"{a2}_{beta}"] == "cellA_contig_2"
        # Neither variant kept the merged clone's dominant-only list.
        assert "DOMINANT_only" not in "".join(map(str, cid.values()))

    def test_slot_order_independence(self):
        # Rep cell has the alphas in the OPPOSITE slots; mapping must follow the
        # cdr3→slot map, not row position.
        a1, a2, beta = "CAAAF", "CBBBF", "CASSBETAF"
        obs = pd.DataFrame([{
            "TRA_1_cdr3": a2, "TRA_2_cdr3": a1, "TRB_1_cdr3": beta,
            "TRA_1_umis": 9, "TRA_2_umis": 9, "TRB_1_umis": 9,
            "TRA_1_vdj_aa": a2, "TRA_1_vdj_nt": back_translate(a2),
            "TRA_1_j_gene": "TRAJ2", "TRA_1_contig_id": "slot1_contig",
            "TRA_2_vdj_aa": a1, "TRA_2_vdj_nt": back_translate(a1),
            "TRA_2_j_gene": "TRAJ1", "TRA_2_contig_id": "slot2_contig",
        }])
        out, _ = expand_dual_alpha_variants(self._clone(a1, a2, beta), obs)
        cid = dict(zip(out["CDR3ab"], out["alpha_contig_ids"]))
        jg = dict(zip(out["CDR3ab"], out["alpha_j_gene"]))
        # a1 lives in TRA_2 here → must get slot2's contig + J gene.
        assert cid[f"{a1}_{beta}"] == "slot2_contig" and jg[f"{a1}_{beta}"] == "TRAJ1"
        assert cid[f"{a2}_{beta}"] == "slot1_contig" and jg[f"{a2}_{beta}"] == "TRAJ2"

    def test_obs_without_contig_id_columns_clears_to_none(self):
        # Older h5ad lacking TRA_*_contig_id: no per-α contig is recoverable, so
        # the variant's alpha_contig_ids is cleared to None — NOT left as the
        # merged clone's dominant-α list (which would feed the WRONG α into
        # extraction). None → extraction falls back cleanly.
        a1, a2, beta = "CAAAF", "CBBBF", "CASSBETAF"
        obs = pd.DataFrame([{
            "TRA_1_cdr3": a1, "TRA_2_cdr3": a2, "TRB_1_cdr3": beta,
            "TRA_1_umis": 9, "TRA_2_umis": 7, "TRB_1_umis": 9,
            "TRA_1_vdj_aa": a1, "TRA_1_vdj_nt": back_translate(a1),
            "TRA_2_vdj_aa": a2, "TRA_2_vdj_nt": back_translate(a2),
        }])
        out, _ = expand_dual_alpha_variants(self._clone(a1, a2, beta), obs)
        assert out["alpha_contig_ids"].isna().all()

    def test_one_slot_missing_contig_id_is_per_variant(self):
        # Only TRA_1 has a contig id; the TRA_2 variant must NOT inherit TRA_1's
        # (or the dominant-α clone list) — it clears to None.
        a1, a2, beta = "CAAAF", "CBBBF", "CASSBETAF"
        obs = pd.DataFrame([{
            "TRA_1_cdr3": a1, "TRA_2_cdr3": a2, "TRB_1_cdr3": beta,
            "TRA_1_umis": 9, "TRA_2_umis": 7, "TRB_1_umis": 9,
            "TRA_1_vdj_aa": a1, "TRA_1_vdj_nt": back_translate(a1),
            "TRA_1_contig_id": "only_a1_contig",
            "TRA_2_vdj_aa": a2, "TRA_2_vdj_nt": back_translate(a2),
            "TRA_2_contig_id": float("nan"),  # a2's contig unknown
        }])
        out, _ = expand_dual_alpha_variants(self._clone(a1, a2, beta), obs)
        cid = dict(zip(out["CDR3ab"], out["alpha_contig_ids"]))
        assert cid[f"{a1}_{beta}"] == "only_a1_contig"
        # a2 has no contig → None (never a1's contig, never the dominant list).
        assert pd.isna(cid[f"{a2}_{beta}"])

    def test_multi_cell_contigs_unioned_per_alpha(self):
        # Each α's contig list spans ALL co-expressing cells (not just the rep
        # cell), de-duplicated — preserving the multi-contig consensus / max
        # C-region coverage the merged clone's list had.
        a1, a2, beta = "CAAAF", "CBBBF", "CASSBETAF"
        common = dict(
            TRB_1_cdr3=beta, TRA_1_cdr3=a1, TRA_2_cdr3=a2,
            TRA_1_vdj_aa=a1, TRA_1_vdj_nt=back_translate(a1),
            TRA_2_vdj_aa=a2, TRA_2_vdj_nt=back_translate(a2),
        )
        obs = pd.DataFrame([
            {**common, "TRA_1_umis": 9, "TRA_2_umis": 9, "TRB_1_umis": 9,
             "TRA_1_contig_id": "cellA_a1", "TRA_2_contig_id": "cellA_a2"},
            {**common, "TRA_1_umis": 3, "TRA_2_umis": 3, "TRB_1_umis": 3,
             "TRA_1_contig_id": "cellB_a1", "TRA_2_contig_id": "cellA_a2"},  # dup a2
        ])
        out, _ = expand_dual_alpha_variants(self._clone(a1, a2, beta), obs)
        cid = dict(zip(out["CDR3ab"], out["alpha_contig_ids"]))
        # a1 seen in two cells → both contigs; a2 duplicated → de-duped to one.
        assert set(cid[f"{a1}_{beta}"].split(";")) == {"cellA_a1", "cellB_a1"}
        assert cid[f"{a2}_{beta}"] == "cellA_a2"

    def test_three_alphas_not_expanded(self):
        # >2 α in merged_alpha_partners: expansion only handles the 2-α case;
        # the clone passes through unexpanded (no bogus contig assignment).
        a1, a2, a3, beta = "CAAAF", "CBBBF", "CCCCF", "CASSBETAF"
        clone = pd.DataFrame([{
            "CDR3ab": f"{a1}_{beta}", "CDR3_alpha": a1, "CDR3_beta": beta,
            "merged_alpha_partners": f"{a1};{a2};{a3}",
            "alpha_contig_ids": "kept",
        }])
        obs = _obs_with_dual_alpha(a1, a2, beta)
        out, variant_of = expand_dual_alpha_variants(clone, obs)
        assert list(out["CDR3ab"]) == [f"{a1}_{beta}"] and variant_of == {}
        assert out["alpha_contig_ids"].iloc[0] == "kept"


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
            allow_canonical_fallback=True,  # synthetic fixture: no contigs
        )
        assert (tmp_path / "selected_clones.csv").exists()
        assert (tmp_path / "selected_clones_sequences.pdf").exists()
        assert (tmp_path / "selected_clones_qc.txt").exists()
        # Provenance carried through to the assembled output.
        written = pd.read_csv(tmp_path / "selected_clones.csv")
        assert "selection_rule" in written.columns
        assert "single_chain_aa" in out.columns

    def test_emits_qc_summary(self, tmp_path):
        # #279: the independent QC battery writes qc-summary.md on every
        # selected-report build. With no contig dir, check E SKIPs; A/B/C/D run
        # on the assembled frame. This synthetic fixture has no contigs, so the
        # construct is a canonical fallback and check D correctly flags it as
        # not contig-verified — exactly the integrity signal the battery is for.
        a = "CASS" + "A" * 40 + "VLF"
        b = "CASS" + "G" * 40 + "VEF"
        clonotypes = pd.DataFrame([_assembleable_clone("c1", a, b)])
        selected = pd.DataFrame([{"CDR3ab": "c1", "selection_rule": "shared"}])
        out_dir = tmp_path / "B1-2"  # donor-named dir → cohort header
        build_selected_report(
            selected, clonotypes, out_dir,
            provenance_cols=["selection_rule"],
            allow_canonical_fallback=True,
        )
        summary = out_dir / "qc-summary.md"
        assert summary.exists()
        text = summary.read_text()
        assert "QC summary" in text
        assert "B1-2" in text
        assert "A. Assembly integrity" in text
        # Assembly integrity (A) must hold even for a canonical-fallback build.
        assert "A. Assembly integrity | PASS" in text
        # No contigs → E SKIPs (not a failure).
        assert "E. Raw cellranger contigs | SKIP" in text
        # Edge-case allow-list note is always present.
        assert "TRAJ35" in text

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
            allow_canonical_fallback=True,  # synthetic fixture: no contigs
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
            allow_canonical_fallback=True,  # synthetic fixture: no contigs
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
        build_selected_report(selected, clonotypes, tmp_path, cover=False, allow_canonical_fallback=True)
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
        # (these donor PDFs are blank figures, not legend covers, so nothing is
        # stripped — matches the no-cover path.)
        assert len(pages) == 8

    def test_combine_strips_redundant_donor_legend(self, tmp_path):
        # #combined-pdf: a donor PDF built WITH its own legend cover must not have
        # that cover repeated once the cohort legend is shown up front.
        from pypdf import PdfReader

        from tcrsift.report import combine_selected_pdfs

        a = "CASS" + "A" * 40 + "VLF"
        b = "CASS" + "G" * 40 + "VEF"
        donor_pdfs = []
        for name in ("B1-2", "B1-3"):
            d = tmp_path / name
            clonotypes = pd.DataFrame([_assembleable_clone("c1", a, b)])
            selected = pd.DataFrame([{"CDR3ab": "c1", "selection_rule": "shared"}])
            build_selected_report(
                selected, clonotypes, d, provenance_cols=["selection_rule"],
                cover=True, allow_canonical_fallback=True,
            )
            donor_pdfs.append(d / "selected_clones_sequences.pdf")
        # each donor PDF = legend cover + 1 construct page (2 pages).
        assert all(len(PdfReader(str(p)).pages) == 2 for p in donor_pdfs)
        out = tmp_path / "cohort.pdf"
        combine_selected_pdfs(donor_pdfs, out, labels=["B1-2", "B1-3"])
        pages = PdfReader(str(out)).pages
        # cover + legend + per donor [title + 1 construct] x2 = 2 + 4 = 6
        # (each donor's embedded legend cover stripped, not 8).
        assert len(pages) == 6
        # only ONE legend in the whole file (the cohort one).
        n_legend = sum(
            1 for p in pages if "ribosomal" in (p.extract_text() or "")
        )
        assert n_legend == 1

    def test_combine_keeps_legend_when_no_cohort_legend(self, tmp_path):
        # With include_legend=False the donor covers are the only legends and must
        # be kept (don't strip what isn't duplicated).
        from pypdf import PdfReader

        from tcrsift.report import combine_selected_pdfs

        a = "CASS" + "A" * 40 + "VLF"
        b = "CASS" + "G" * 40 + "VEF"
        d = tmp_path / "B1-2"
        clonotypes = pd.DataFrame([_assembleable_clone("c1", a, b)])
        selected = pd.DataFrame([{"CDR3ab": "c1", "selection_rule": "shared"}])
        build_selected_report(
            selected, clonotypes, d, provenance_cols=["selection_rule"],
            cover=True, allow_canonical_fallback=True,
        )
        p = d / "selected_clones_sequences.pdf"
        out = tmp_path / "cohort.pdf"
        combine_selected_pdfs([p], out, labels=["B1-2"], include_legend=False)
        pages = PdfReader(str(out)).pages
        # cohort title + donor title + [legend cover + construct] = 4 (cover kept)
        assert len(pages) == 4
        assert any("ribosomal" in (pg.extract_text() or "") for pg in pages)


def test_selection_breakdown_cols_match_helper():
    # The PDF skip-set must stay in sync with the helper's emitted columns:
    # every per-condition breakdown column except the formatted `selection_detail`
    # is CSV-only (#selection-cols). Guards against drift when a column is added.
    from tcrsift.report import _SELECTION_BREAKDOWN_COLS
    from tcrsift.selection import SELECTION_SUMMARY_COLS

    assert _SELECTION_BREAKDOWN_COLS == set(SELECTION_SUMMARY_COLS) - {"selection_detail"}


def test_attach_public_db_annotation():
    # Curated public-DB columns are surfaced on the selected clones; a dual-α
    # variant inherits its parent's annotation; absent annotations → empty cols.
    from tcrsift.report import _PUBLIC_DB_COLS, _attach_public_db_annotation
    assembled = pd.DataFrame({"CDR3ab": ["c1", "c2", "c1_variant"], "x": [1, 2, 3]})
    anno = pd.DataFrame({
        "CDR3ab": ["c1", "c2"],
        "is_viral": [True, False],
        "db_category": ["viral", "tumor"],
        "db_epitope": ["GLCTLVAML", "ELAGIGILTV"],
        "db_species": ["EBV", "human"],
        "db_database": ["IEDB", "CEDAR"],
    })
    variant_of = {"c1_variant": "c1"}  # variant inherits c1's annotation
    out = _attach_public_db_annotation(assembled.copy(), anno, variant_of)
    assert out.loc[out.CDR3ab == "c1", "db_epitope"].iloc[0] == "GLCTLVAML"
    assert out.loc[out.CDR3ab == "c1_variant", "db_category"].iloc[0] == "viral"  # inherited
    assert out.loc[out.CDR3ab == "c2", "is_viral"].iloc[0] == False  # noqa: E712
    # schema stability: all curated cols present even with no annotations
    out2 = _attach_public_db_annotation(assembled.copy(), None, {})
    assert all(c in out2.columns for c in _PUBLIC_DB_COLS)


def test_emit_frequency_by_condition(tmp_path):
    from tcrsift.report import _emit_frequency_by_condition
    df = pd.DataFrame({
        "CDR3ab": ["c1", "c2"],
        "frequency_per_condition": ["AIM⁺=0.90%;CTY⁻=0.34%", "AIM⁺=0.50%"],
    })
    pivot = _emit_frequency_by_condition(df, tmp_path, name="B1-2")
    assert pivot is not None
    assert set(pivot.columns) == {"AIM⁺", "CTY⁻"}
    assert pivot.loc["c1", "AIM⁺"] == 0.90 and pivot.loc["c1", "CTY⁻"] == 0.34
    assert pd.isna(pivot.loc["c2", "CTY⁻"])  # c2 not in CTY⁻
    assert (tmp_path / "selected_frequency_by_condition.csv").exists()
    assert (tmp_path / "selected_frequency_heatmap.png").exists()


def test_emit_frequency_by_condition_skips_without_column(tmp_path):
    from tcrsift.report import _emit_frequency_by_condition
    df = pd.DataFrame({"CDR3ab": ["c1"]})  # no frequency_per_condition
    assert _emit_frequency_by_condition(df, tmp_path) is None


def test_frequency_heatmap_carries_annotation_strip(tmp_path):
    # #selected-anno: the freq heatmap's public-DB annotation strip is driven by
    # db_category on the assembled frame.
    from tcrsift.plots import plot_selected_frequency_heatmap
    pivot = pd.DataFrame(
        {"AIM⁺": [0.9, 0.5, 0.1], "CTY⁻": [0.3, None, 0.2]},
        index=pd.Index(["c1", "c2", "c3"], name="CDR3ab"),
    )
    anno = pd.Series({"c1": "viral", "c2": None, "c3": "tumor_self"})
    out = plot_selected_frequency_heatmap(
        pivot, tmp_path / "hm.png", title="B1-2", annotations=anno,
    )
    assert out is not None and out.exists() and out.stat().st_size > 0
    # also renders fine with no annotations
    assert plot_selected_frequency_heatmap(pivot, tmp_path / "hm2.png").exists()
