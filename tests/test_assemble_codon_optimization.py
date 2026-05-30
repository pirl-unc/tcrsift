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

"""Tests for #116 — motif-aware codon optimization + dual stops +
NT-triad (assembly / optimized / blend).

Covers:

* :func:`optimize_codons` motif avoidance — homopolymer runs and
  common Type-II restriction sites do not appear in the output.
* :func:`optimize_codons` AA round-trip — translates back to the
  input.
* :func:`optimize_codons` determinism — same input → same output.
* :func:`stop_codons_nt` validation + concatenation.
* ``stop_codons`` plumbing through ``assemble_full_sequences`` —
  default dual stops appear at the tail of every codon-optimized
  constant, single-stop preserves pre-2.4 behavior, empty tuple
  omits stops entirely, invalid stops raise ``ValueError``.
* NT-triad columns (``_constant_nt_assembly`` /
  ``_constant_nt_optimized`` / ``_constant_nt``) all appear on
  output and obey their documented round-trip invariants.
* :func:`_add_single_chain` strips ALL trailing stops, not just one
  (single-chain construct stays in frame with the 2A linker under
  dual-stop default).
* The picker stop-codon bug from #116: user-overridden alleles
  and auto-picked alleles both end with the configured stops
  (pre-#116 they silently dropped them).
"""

from __future__ import annotations

import pandas as pd
import pytest

from tcrsift.assemble import (
    DEFAULT_FORBIDDEN_MOTIFS,
    HUMAN_CODON_ALTERNATIVES,
    HUMAN_PREFERRED_CODONS,
    HUMAN_TRBC1_AA,
    assemble_full_sequences,
    back_translate,
    get_constant_region_sequences,
    optimize_codons,
    stop_codons_nt,
    translate_dna,
)


class TestOptimizeCodonsBasics:
    """Core invariants of the in-house codon optimizer."""

    def test_translates_back_to_input_aa(self):
        aa = "MASTKLDEFGHIKLMNPQRSTVWY"
        nt = optimize_codons(aa)
        translated, ragged = translate_dna(nt)
        assert not ragged
        assert translated == aa

    def test_translates_back_with_stop(self):
        aa = "MAS*"
        nt = optimize_codons(aa)
        # Stop should map to TAA / TGA / TAG.
        assert nt[-3:] in {"TAA", "TGA", "TAG"}

    def test_deterministic(self):
        aa = "M" + ("ASKLDE" * 30)  # long enough that motif logic kicks in
        assert optimize_codons(aa) == optimize_codons(aa)

    def test_no_motifs_in_canonical_constants(self):
        # The real motivation: the actual canonical TCR constants —
        # which are mixed-AA sequences — come out clean.
        from tcrsift.assemble import HUMAN_CONSTANT_REGIONS_AA
        for gene, aa in HUMAN_CONSTANT_REGIONS_AA.items():
            nt = optimize_codons(aa)
            for motif in DEFAULT_FORBIDDEN_MOTIFS:
                assert motif not in nt, (
                    f"forbidden motif {motif!r} appeared in optimized "
                    f"{gene}: {aa!r}"
                )

    def test_avoids_homopolymer_at_modest_run(self):
        # Greedy avoidance works for realistic AA contexts (≤4 same
        # AA in a row, which is what TCR constants actually have).
        # Pathological runs of 8+ same AA can defeat any greedy
        # strategy without lookahead — by design.
        for aa in ("MKKKK", "MAAAA", "MGGGG", "MPPPP"):
            nt = optimize_codons(aa)
            for motif in ("AAAAA", "GGGGG", "CCCCC", "TTTTT"):
                assert motif not in nt, (
                    f"avoidable homopolymer {motif!r} appeared in "
                    f"optimized {aa!r}: {nt!r}"
                )

    def test_unknown_residue_emits_nnn(self):
        nt = optimize_codons("MX")
        # M is fixed (ATG), X has no entry → falls back to "NNN".
        assert nt == "ATGNNN"

    def test_prefix_motif_already_present_is_not_blamed_on_new_codon(self):
        # Prefix already contains GGGGG. Picking ANY G-starting codon
        # for the next G would "extend" the run, but the motif is
        # pre-existing — we can't rewrite the prefix. Optimizer falls
        # back to CAI-best (GGC) rather than rejecting all options.
        nt = optimize_codons("G", prefix_nt="GGGGG")
        # CAI-best for G — picker accepts since GGGGG was already in
        # the prefix.
        assert nt == "GGC"

    def test_avoids_introducing_motif_at_codon_boundary(self):
        # Prefix ends in "CCCC" but no CCCCC yet. The next codon must
        # not introduce one: codons for I (ATC/ATT/ATA) all start
        # with A, so any is fine. Test with K instead: alternatives
        # AAG / AAA — A-starting. Pick a sequence where the boundary
        # actually matters: prefix "CCCC", AA "P" (alternatives start
        # with C). Should NOT pick CCC (would form CCCCCCC) — should
        # pick CCT/CCA/CCG.
        nt = optimize_codons("P", prefix_nt="CCCC")
        assert nt != "CCC", (
            "optimizer picked CCC after CCCC prefix → CCCCCCC, "
            "should have avoided motif"
        )
        assert nt in {"CCT", "CCA", "CCG"}

    def test_forbidden_motifs_can_be_disabled(self):
        # Pass empty tuple → no motif checks → equivalent to back_translate.
        aa = "KKKKKK"
        nt = optimize_codons(aa, forbidden_motifs=())
        assert nt == "AAG" * 6  # CAI-best for every K


class TestExpressionQualityMotifs:
    """The default motif set extends past restriction sites and
    homopolymer runs to include polyA + splice donor consensus
    (#116 expression-quality additions)."""

    @pytest.mark.parametrize("motif,kind", [
        ("AATAAA", "polyA canonical"),
        ("ATTAAA", "polyA cryptic"),
        ("GTAAGT", "5' splice donor canonical"),
        ("GTGAGT", "5' splice donor variant"),
    ])
    def test_canonical_constants_free_of_motif(self, motif, kind):
        # The default optimization of every canonical constant must
        # avoid these expression-killing motifs.
        seqs = get_constant_region_sequences()
        for gene, nt in seqs.items():
            assert motif not in nt, (
                f"{kind} motif {motif!r} appeared in optimized {gene}"
            )

    def test_polya_motif_actively_avoided(self):
        # N-K is the natural source of AATAAA: N→AAT, K→AAA → "AATAAA".
        # Optimizer should pick N→AAC and/or K→AAG to break the motif.
        nt = optimize_codons("NK")
        assert "AATAAA" not in nt
        # Round-trip check (the avoidance shouldn't break translation).
        translated, _ = translate_dna(nt)
        assert translated == "NK"

    def test_splice_donor_avoided_at_VS_junction(self):
        # V-S can spell GTAAGT when V picks GTA and S picks AGT.
        # CAI-best is V→GTG, S→AGC → GTGAGC, which is also a 5' splice
        # donor variant. With the new GTGAGT in forbidden_motifs the
        # optimizer should reject GTGAGT and pick alternatives.
        nt = optimize_codons("VS")
        assert "GTAAGT" not in nt
        assert "GTGAGT" not in nt
        translated, _ = translate_dna(nt)
        assert translated == "VS"


class TestGcWindowing:
    """GC-content balancing (#116). Sliding window of 60 nt should be
    biased toward 40–65% GC where alternatives permit."""

    def test_gc_target_disabled_recovers_cai_pick(self):
        # With GC balancing disabled, K should always pick CAI-best AAG.
        nt = optimize_codons("KKKK", gc_target_range=(0.0, 1.0))
        # CAI-best K is AAG. (Without GC bias, all 4 K's get AAG.)
        assert nt == "AAGAAGAAGAAG"

    def test_default_gc_keeps_canonical_in_range(self):
        # The codon-optimized canonical constants should sit inside
        # the GC target range when averaged.
        seqs = get_constant_region_sequences()
        for gene, nt in seqs.items():
            # Drop trailing stops before measuring.
            body = nt
            while body[-3:] in {"TAA", "TAG", "TGA"}:
                body = body[:-3]
            gc = body.count("G") + body.count("C")
            frac = gc / len(body)
            assert 0.35 <= frac <= 0.70, (
                f"{gene} optimized GC out of broad range: {frac:.2%}"
            )

    def test_gc_invalid_range_raises(self):
        with pytest.raises(ValueError, match="gc_target_range"):
            optimize_codons("M", gc_target_range=(0.7, 0.3))


class TestOptimizeCodonsAvoidsRestrictionSites:
    """Restriction-site avoidance is the practical motivation for
    optimization — if these sites appear in synthesized constructs
    they break common cloning workflows."""

    @pytest.mark.parametrize("site,length", [
        ("GAATTC", 6),   # EcoRI
        ("GGATCC", 6),   # BamHI
        ("GCGGCCGC", 8), # NotI
        ("CTCGAG", 6),   # XhoI
    ])
    def test_no_restriction_site_in_canonical_constants(self, site, length):
        # The full canonical constants, fully optimized, should
        # not contain these common cloning sites.
        seqs = get_constant_region_sequences()
        for gene, nt in seqs.items():
            assert site not in nt, (
                f"restriction site {site} appeared in optimized {gene} NT"
            )


class TestStopCodonsNt:
    """Validation + concatenation helper for stop codons."""

    def test_default_dual_stop(self):
        assert stop_codons_nt(("TAA", "TGA")) == "TAATGA"

    def test_single_stop(self):
        assert stop_codons_nt(("TAA",)) == "TAA"

    def test_empty_returns_empty_string(self):
        assert stop_codons_nt(()) == ""

    def test_triple_stop(self):
        assert stop_codons_nt(("TAA", "TAG", "TGA")) == "TAATAGTGA"

    def test_rejects_non_stop_codon(self):
        with pytest.raises(ValueError, match="must be one of"):
            stop_codons_nt(("TAA", "GCC"))

    def test_rejects_lowercase(self):
        # Strict — we don't silently upper-case; downstream pipelines
        # are case-sensitive and a lower-case codon usually means the
        # user mixed up DNA / RNA conventions.
        with pytest.raises(ValueError, match="must be one of"):
            stop_codons_nt(("taa",))


class TestGetConstantRegionSequencesStops:
    """``get_constant_region_sequences`` should respect ``stop_codons``."""

    def test_default_has_dual_stop_tail(self):
        seqs = get_constant_region_sequences()
        for gene, nt in seqs.items():
            assert nt.endswith("TAATGA"), (
                f"{gene} should end with default dual stop (TAATGA); "
                f"got tail {nt[-12:]!r}"
            )

    def test_single_stop_arg(self):
        seqs = get_constant_region_sequences(stop_codons=("TAG",))
        for gene, nt in seqs.items():
            assert nt.endswith("TAG")
            # Should NOT end with two stops.
            second_to_last = nt[-6:-3]
            assert second_to_last not in {"TAA", "TGA", "TAG"}, (
                f"{gene} should have only one trailing stop; "
                f"got tail {nt[-12:]!r}"
            )

    def test_no_stops(self):
        seqs = get_constant_region_sequences(stop_codons=())
        for gene, nt in seqs.items():
            tail = nt[-3:]
            assert tail not in {"TAA", "TAG", "TGA"}, (
                f"{gene} should not have a trailing stop with stop_codons=(); "
                f"got tail {tail!r}"
            )

    def test_translates_back_to_canonical(self):
        # Each NT (minus stops) must translate to the canonical AA.
        from tcrsift.assemble import HUMAN_CONSTANT_REGIONS_AA
        seqs = get_constant_region_sequences()
        for gene, nt in seqs.items():
            body = nt
            while body[-3:] in {"TAA", "TAG", "TGA"}:
                body = body[:-3]
            translated, ragged = translate_dna(body)
            assert not ragged, f"{gene} NT not multiple-of-3 after stop strip"
            assert translated == HUMAN_CONSTANT_REGIONS_AA[gene]


class TestStopCodonsCli:
    """``assemble_full_sequences`` propagates and validates
    ``stop_codons``."""

    @staticmethod
    def _minimal_clonotype():
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        return pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha,
            "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha,
            "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": back_translate(vdj_alpha),
            "VDJ_beta_nt": back_translate(vdj_beta),
            "alpha_c_gene": "TRAC",
            "beta_c_gene": "TRBC1",
            "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
        }])

    def test_default_dual_stop_on_constant_nt(self):
        out = assemble_full_sequences(
            self._minimal_clonotype(),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        # Both _constant_nt (blend) and _constant_nt_optimized
        # should end in the dual stop.
        for chain in ("alpha", "beta"):
            for col_suffix in ("_constant_nt", "_constant_nt_optimized"):
                nt = out[f"{chain}{col_suffix}"].iloc[0]
                assert nt.endswith("TAATGA"), (
                    f"{chain}{col_suffix} should end with TAATGA; "
                    f"got tail {nt[-12:]!r}"
                )

    def test_single_stop_propagates(self):
        out = assemble_full_sequences(
            self._minimal_clonotype(),
            alpha_leader=None, beta_leader=None,
            stop_codons=("TAA",),
            verbose=False, show_progress=False,
        )
        for chain in ("alpha", "beta"):
            nt = out[f"{chain}_constant_nt_optimized"].iloc[0]
            assert nt.endswith("TAA")
            assert nt[-6:-3] not in {"TAA", "TAG", "TGA"}

    def test_no_stops(self):
        out = assemble_full_sequences(
            self._minimal_clonotype(),
            alpha_leader=None, beta_leader=None,
            stop_codons=(),
            verbose=False, show_progress=False,
        )
        for chain in ("alpha", "beta"):
            nt = out[f"{chain}_constant_nt_optimized"].iloc[0]
            assert nt[-3:] not in {"TAA", "TAG", "TGA"}

    def test_invalid_stop_raises(self):
        with pytest.raises(ValueError, match="must be one of"):
            assemble_full_sequences(
                self._minimal_clonotype(),
                stop_codons=("TAA", "ATG"),
                alpha_leader=None, beta_leader=None,
                verbose=False, show_progress=False,
            )


class TestNtTriadColumns:
    """The new ``_assembly`` / ``_optimized`` / ``_constant_nt`` triad."""

    @staticmethod
    def _df_no_contig():
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        return pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha,
            "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha,
            "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": back_translate(vdj_alpha),
            "VDJ_beta_nt": back_translate(vdj_beta),
            "alpha_c_gene": "TRAC",
            "beta_c_gene": "TRBC1",
            "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
        }])

    def test_columns_present_no_contig(self):
        out = assemble_full_sequences(
            self._df_no_contig(),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        for chain in ("alpha", "beta"):
            assert f"{chain}_constant_nt_assembly" in out.columns
            assert f"{chain}_constant_nt_optimized" in out.columns
            assert f"{chain}_constant_nt" in out.columns

    def test_assembly_none_without_contig(self):
        # Without contigs there is no donor NT to put in _assembly.
        out = assemble_full_sequences(
            self._df_no_contig(),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        for chain in ("alpha", "beta"):
            assert out[f"{chain}_constant_nt_assembly"].iloc[0] is None
            # full_{chain}_nt_assembly likewise None.
            assert out[f"full_{chain}_nt_assembly"].iloc[0] is None

    def test_optimized_present_without_contig(self):
        out = assemble_full_sequences(
            self._df_no_contig(),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        for chain in ("alpha", "beta"):
            opt = out[f"{chain}_constant_nt_optimized"].iloc[0]
            assert isinstance(opt, str)
            assert opt.endswith("TAATGA")

    def test_assembly_and_optimized_with_contig(self, tmp_path):
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        leader_aa = "M" + "A" * 19
        vdj_beta_nt = back_translate(vdj_beta)
        # Contig: leader + VDJ + first 8 codons of canonical TRBC1.
        # First codon serves as the junction E (HUMAN_TRBC1_AA[0]).
        c_region_start_nt = back_translate(HUMAN_TRBC1_AA[:8])
        beta_contig_seq = (
            back_translate(leader_aa) + vdj_beta_nt + c_region_start_nt
        )

        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha,
            "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha,
            "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": back_translate(vdj_alpha),
            "VDJ_beta_nt": vdj_beta_nt,
            "alpha_c_gene": "TRAC",
            "beta_c_gene": "TRBC1",
            "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
            "beta_contig_ids": "contig_b1",
        }])
        contig_dir = tmp_path / "S1"
        contig_dir.mkdir(parents=True)
        (contig_dir / "filtered_contig.fasta").write_text(
            f">contig_b1\n{beta_contig_seq}\n"
        )

        out = assemble_full_sequences(
            df, contigs_dir=str(tmp_path),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )

        # _assembly is the raw contig bytes past VDJ — exactly
        # c_region_start_nt for this fixture.
        assembly = out["beta_constant_nt_assembly"].iloc[0]
        assert assembly == c_region_start_nt
        # _optimized is the codon-optimized canonical + dual stops.
        opt = out["beta_constant_nt_optimized"].iloc[0]
        assert opt.endswith("TAATGA")
        # The two should differ (assembly has no stops; optimized has them).
        assert opt != assembly

    def test_optimized_translates_to_constant_aa(self):
        out = assemble_full_sequences(
            self._df_no_contig(),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        for chain in ("alpha", "beta"):
            nt = out[f"{chain}_constant_nt_optimized"].iloc[0]
            aa = out[f"{chain}_constant_aa"].iloc[0]
            body = nt
            while body[-3:] in {"TAA", "TAG", "TGA"}:
                body = body[:-3]
            translated, _ = translate_dna(body)
            assert translated == aa


class TestSingleChainDualStopStripping:
    """The 2A construct requires the upstream β CDS to have NO trailing
    stop. With dual stops the strip logic must remove ALL of them."""

    def test_single_chain_no_residual_stop_before_linker(self):
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha,
            "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha,
            "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": back_translate(vdj_alpha),
            "VDJ_beta_nt": back_translate(vdj_beta),
            "alpha_c_gene": "TRAC",
            "beta_c_gene": "TRBC1",
            "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
        }])
        out = assemble_full_sequences(
            df, alpha_leader=None, beta_leader=None,
            linker="T2A", verbose=False, show_progress=False,
        )
        sc_nt = out["single_chain_nt"].iloc[0]
        # Linker should appear once, with NO upstream stop right before it.
        from tcrsift.assemble import LINKERS
        linker_nt = LINKERS["T2A"]["dna"]
        # Find where the linker starts in the single chain.
        idx = sc_nt.find(linker_nt)
        assert idx > 0
        # The 3 nt just before the linker must NOT be a stop codon.
        prev_codon = sc_nt[idx - 3 : idx]
        assert prev_codon not in {"TAA", "TAG", "TGA"}, (
            f"single_chain_nt has a residual stop {prev_codon!r} just before "
            f"the T2A linker — strip_stop_codon_dna failed to remove all "
            f"trailing stops from the dual-stop β CDS."
        )
        # And the single-chain NT in frame translates to the AA.
        sc_aa = out["single_chain_aa"].iloc[0]
        translated, _ = translate_dna(sc_nt)
        assert translated == sc_aa, (
            "single_chain_nt does not translate to single_chain_aa — "
            "frame broke at the β→linker junction"
        )


class TestSingleChainTriad:
    """The single-chain (β-2A-α) cassette is emitted in three flavors
    matching the constant-region triad (#116): blend / optimized /
    assembly. Synthesis users want ``_optimized``; QC users want
    ``_assembly``; back-compat callers keep ``_constant_nt`` and
    ``single_chain_nt`` unchanged."""

    @staticmethod
    def _df_no_contig():
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        return pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha,
            "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha,
            "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": back_translate(vdj_alpha),
            "VDJ_beta_nt": back_translate(vdj_beta),
            "alpha_c_gene": "TRAC",
            "beta_c_gene": "TRBC1",
            "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
        }])

    def test_triad_columns_present(self):
        out = assemble_full_sequences(
            self._df_no_contig(),
            alpha_leader=None, beta_leader=None,
            linker="T2A",
            verbose=False, show_progress=False,
        )
        for col in (
            "single_chain_nt",
            "single_chain_nt_optimized",
            "single_chain_nt_assembly",
        ):
            assert col in out.columns

    def test_optimized_translates_to_single_chain_aa(self):
        out = assemble_full_sequences(
            self._df_no_contig(),
            alpha_leader=None, beta_leader=None,
            linker="T2A",
            verbose=False, show_progress=False,
        )
        nt = out["single_chain_nt_optimized"].iloc[0]
        aa = out["single_chain_aa"].iloc[0]
        translated, _ = translate_dna(nt)
        assert translated == aa

    def test_assembly_is_none_when_no_contig(self):
        # Without contigs, full_*_nt_assembly is None for both chains,
        # so the single-chain assembly must be None too.
        out = assemble_full_sequences(
            self._df_no_contig(),
            alpha_leader=None, beta_leader=None,
            linker="T2A",
            verbose=False, show_progress=False,
        )
        assert out["single_chain_nt_assembly"].iloc[0] is None

    def test_optimized_ends_with_dual_stop(self):
        # The β CDS in the cassette has all stops stripped (the
        # construct is one ORF across the linker), but the ALPHA
        # constant's stops survive at the 3' end.
        out = assemble_full_sequences(
            self._df_no_contig(),
            alpha_leader=None, beta_leader=None,
            linker="T2A",
            verbose=False, show_progress=False,
        )
        nt = out["single_chain_nt_optimized"].iloc[0]
        assert nt.endswith("TAATGA"), (
            f"single_chain_nt_optimized should end in dual stop; tail "
            f"is {nt[-12:]!r}"
        )

    def test_optimized_has_no_internal_stop_before_linker(self):
        # The 2A linker AA must appear once. The NT just before the
        # linker NT must NOT be a stop codon — those should all be
        # stripped from β.
        out = assemble_full_sequences(
            self._df_no_contig(),
            alpha_leader=None, beta_leader=None,
            linker="T2A",
            verbose=False, show_progress=False,
        )
        from tcrsift.assemble import LINKERS
        linker_nt = LINKERS["T2A"]["dna"]
        for col in ("single_chain_nt", "single_chain_nt_optimized"):
            nt = out[col].iloc[0]
            idx = nt.find(linker_nt)
            assert idx > 0
            assert nt[idx - 3 : idx] not in {"TAA", "TAG", "TGA"}, (
                f"{col} has residual stop before T2A — "
                f"strip_stop_codon_dna missed a trailing stop"
            )

    def test_optimized_differs_from_blend_at_constant_only(self):
        # The blend uses canonical-optimized C-region bytes when no
        # contig is provided, so the blend equals the optimized
        # variant in this no-contig fixture. Verify the equivalence
        # — it locks in the back-compat invariant.
        out = assemble_full_sequences(
            self._df_no_contig(),
            alpha_leader=None, beta_leader=None,
            linker="T2A",
            verbose=False, show_progress=False,
        )
        assert (
            out["single_chain_nt"].iloc[0]
            == out["single_chain_nt_optimized"].iloc[0]
        ), (
            "Without contigs the blend and the optimized cassette "
            "must be byte-identical — both use the optimized canonical "
            "for the C region."
        )

    def test_assembly_emitted_when_both_chains_have_contig(self, tmp_path):
        # When BOTH α and β have contig bytes past the J→C junction,
        # the assembly cassette must be a string, end with the
        # donor's actual bytes from the α tail, and be in-frame +
        # premature-stop-free across the whole construct.
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        leader_aa = "M" + "A" * 19
        from tcrsift.assemble import (
            HUMAN_TRAC_AA,
            HUMAN_TRBC1_AA,
        )

        beta_c_start_nt = back_translate(HUMAN_TRBC1_AA[:8])
        alpha_c_start_nt = back_translate(HUMAN_TRAC_AA[:8])
        beta_contig_seq = (
            back_translate(leader_aa)
            + back_translate(vdj_beta)
            + beta_c_start_nt
        )
        alpha_contig_seq = (
            back_translate(leader_aa)
            + back_translate(vdj_alpha)
            + alpha_c_start_nt
        )

        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha,
            "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha,
            "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": back_translate(vdj_alpha),
            "VDJ_beta_nt": back_translate(vdj_beta),
            "alpha_c_gene": "TRAC",
            "beta_c_gene": "TRBC1",
            "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
            "alpha_contig_ids": "contig_a1",
            "beta_contig_ids": "contig_b1",
        }])
        contig_dir = tmp_path / "S1"
        contig_dir.mkdir(parents=True)
        (contig_dir / "filtered_contig.fasta").write_text(
            f">contig_a1\n{alpha_contig_seq}\n>contig_b1\n{beta_contig_seq}\n"
        )

        out = assemble_full_sequences(
            df, contigs_dir=str(tmp_path),
            alpha_leader=None, beta_leader=None,
            linker="T2A",
            verbose=False, show_progress=False,
        )
        sc_asm = out["single_chain_nt_assembly"].iloc[0]
        assert isinstance(sc_asm, str), (
            "single_chain_nt_assembly should be a string when both "
            f"chains have contig coverage, got {type(sc_asm).__name__}"
        )
        # In-frame.
        assert len(sc_asm) % 3 == 0, (
            f"single_chain_nt_assembly length {len(sc_asm)} not divisible "
            "by 3 — frame broken"
        )
        # No mid-chain stops in the coding portion.
        translated, ragged = translate_dna(sc_asm)
        assert not ragged
        assert "*" not in translated[:-1], (
            "single_chain_nt_assembly has a premature stop in the "
            "coding region"
        )


class TestPickerStopBugRegression:
    """#116: the picker overrides used to call ``back_translate`` without
    appending stops, so user-overridden alleles (and any allele picked
    by auto-detect) silently dropped the stops while default-allele
    rows kept them. Verify dual stops are present in both branches."""

    def test_explicit_allele_override_keeps_stops(self):
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha,
            "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha,
            "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": back_translate(vdj_alpha),
            "VDJ_beta_nt": back_translate(vdj_beta),
            "alpha_c_gene": "TRAC",
            "beta_c_gene": "TRBC2",
            "beta_j_gene": "TRBJ2-1",
            "samples": "S1",
        }])
        # Force a non-default allele to exercise the override branch.
        out = assemble_full_sequences(
            df, alpha_leader=None, beta_leader=None,
            trbc2_allele="03",
            verbose=False, show_progress=False,
        )
        nt = out["beta_constant_nt_optimized"].iloc[0]
        assert nt.endswith("TAATGA"), (
            "Explicit allele override branch dropped stops (#116 regression)"
        )
        # And the blended _constant_nt also has stops.
        nt_blend = out["beta_constant_nt"].iloc[0]
        assert nt_blend.endswith("TAATGA"), (
            "Blended constant_nt dropped stops after explicit allele override"
        )


class TestBackTranslateUnchanged:
    """Back-compat: external callers using ``back_translate`` directly
    must continue to get the naive single-best codon output."""

    def test_back_translate_is_deterministic_naive(self):
        aa = "MASTKLDEFGHIKLMNPQRSTVWY"
        nt = back_translate(aa)
        # Every codon is the most-common one for its AA.
        assert nt == "".join(HUMAN_PREFERRED_CODONS[r] for r in aa)

    def test_back_translate_unknown_residue_is_nnn(self):
        assert back_translate("X") == "NNN"

    def test_codon_alternatives_first_is_preferred(self):
        # Sanity: HUMAN_CODON_ALTERNATIVES[aa][0] should equal
        # HUMAN_PREFERRED_CODONS[aa]. Otherwise back_translate +
        # optimize_codons could diverge for boring inputs and confuse
        # callers comparing the two.
        for aa, preferred in HUMAN_PREFERRED_CODONS.items():
            alternatives = HUMAN_CODON_ALTERNATIVES.get(aa, [])
            assert alternatives, f"no alternatives for {aa}"
            assert alternatives[0] == preferred, (
                f"{aa}: preferred {preferred} != first alternative "
                f"{alternatives[0]}"
            )
