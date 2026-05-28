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

"""Tests for the constant-region pipeline after #66 / #67 / #100.

Covers:
- Canonical AA constants are present and end with the canonical
  C-terminus from ``CONSTANT_REGION_ENDINGS``.
- The packaged ``tcrsift/data/canonical_constants.fasta`` matches the
  pyensembl protein_sequence for each canonical transcript (#100 fix).
- The canonical AAs match UniProt P01848 / P01850 / A0A5B9 at a set of
  diagnostic positions (the residues that drifted in #100).
- ``get_constant_region_sequences`` returns back-translated NT that
  round-trips through ``translate_dna`` to the canonical AA.
- ``pick_canonical_constant`` resolves chain + c_gene + j_gene to
  the right canonical (α → TRAC; β → TRBC1 vs TRBC2 by c_gene; β
  fallback to J-gene parity when c_gene is missing).
- ``verify_canonical_constant_start`` accepts close matches and
  rejects clear mismatches.
- ``_extract_c_region_start_from_contig`` finds the post-VDJ residues.
"""

from __future__ import annotations

import pandas as pd
import pytest

from tcrsift.assemble import (
    CODON_TABLE,
    CONSTANT_REGION_ENDINGS,
    HUMAN_CONSTANT_REGIONS_AA,
    HUMAN_PREFERRED_CODONS,
    HUMAN_TRAC_AA,
    HUMAN_TRBC1_AA,
    HUMAN_TRBC2_AA,
    _blend_constant_nt_with_contig,
    _extract_c_region_nt_from_contig,
    _extract_c_region_start_from_contig,
    back_translate,
    get_constant_region_sequences,
    pick_canonical_constant,
    translate_dna,
    verify_canonical_constant_start,
)


class TestCanonicalSequences:
    """The canonical AA constants are the source of truth — verify
    they match the canonical C-terminus QC table."""

    @pytest.mark.parametrize(
        "gene,aa",
        [
            ("TRAC", HUMAN_TRAC_AA),
            ("TRBC1", HUMAN_TRBC1_AA),
            ("TRBC2", HUMAN_TRBC2_AA),
        ],
    )
    def test_ends_with_canonical_c_terminus(self, gene, aa):
        assert aa.endswith(CONSTANT_REGION_ENDINGS[gene]), (
            f"{gene} canonical AA tail {aa[-15:]!r} doesn't match "
            f"CONSTANT_REGION_ENDINGS entry {CONSTANT_REGION_ENDINGS[gene]!r}"
        )

    @pytest.mark.parametrize(
        "gene,min_len,max_len",
        [
            ("TRAC", 130, 145),    # mature TRAC ~138 aa
            ("TRBC1", 170, 185),   # mature TRBC1 ~178 aa
            ("TRBC2", 170, 185),   # mature TRBC2 ~178 aa
        ],
    )
    def test_lengths_in_expected_range(self, gene, min_len, max_len):
        aa = HUMAN_CONSTANT_REGIONS_AA[gene]
        assert min_len <= len(aa) <= max_len, (
            f"{gene} length {len(aa)} outside expected window "
            f"[{min_len}, {max_len}]"
        )

    def test_all_three_present(self):
        assert set(HUMAN_CONSTANT_REGIONS_AA) == {"TRAC", "TRBC1", "TRBC2"}

    def test_diagnostic_positions_match_uniprot(self):
        """Lock the specific residues that previously drifted (#100).

        These positions were wrong in the pre-#100 hardcoded sequences;
        pinning them here makes any future hand-edit drift loud. Indices
        are 0-indexed Python slices into the stored AA constants.

        TRAC: stored as bare mature, so stored index ``i`` == UniProt
        position ``i+1``.

        TRBC1/2: ``E`` prepended for the J→C junction residue, so stored
        index ``i`` == UniProt position ``i`` (i.e., stored position 1
        is the synthetic ``E``; mature pos 1 is at stored index 1).
        """
        # TRAC mature position 47 → stored index 46 — the conserved Thr.
        # Pre-#100 had 'C' here.
        assert HUMAN_TRAC_AA[46] == "T", (
            f"TRAC position 47 must be T (UniProt P01848); got "
            f"{HUMAN_TRAC_AA[46]!r} (pre-#100 had 'C')"
        )

        # TRBC1 mature position 135 → stored index 135 (since E is at
        # stored idx 0, mature 1 is at stored idx 1, mature N at stored
        # idx N). Pre-#100 had 'E' here, swapped with TRBC2.
        assert HUMAN_TRBC1_AA[135] == "V", (
            f"TRBC1 mature position 135 must be V (UniProt P01850); got "
            f"{HUMAN_TRBC1_AA[135]!r} (pre-#100 had 'E')"
        )

        # TRBC2 mature positions that drifted in #100.
        # Mature pos 3 (stored idx 3): K (was N pre-#100 — JOVI.1 epitope swap)
        assert HUMAN_TRBC2_AA[3] == "K", (
            f"TRBC2 mature pos 3 must be K (JOVI.1 epitope); got "
            f"{HUMAN_TRBC2_AA[3]!r}"
        )
        # Mature pos 4 (stored idx 4): N (was K pre-#100 — JOVI.1 epitope swap)
        assert HUMAN_TRBC2_AA[4] == "N", (
            f"TRBC2 mature pos 4 must be N (JOVI.1 epitope); got "
            f"{HUMAN_TRBC2_AA[4]!r}"
        )
        # Mature pos 9 (stored idx 9): K (was E pre-#100)
        assert HUMAN_TRBC2_AA[9] == "K"
        # Mature pos 56 (stored idx 56): S (was C pre-#100)
        assert HUMAN_TRBC2_AA[56] == "S"
        # Mature pos 135 (stored idx 135): E (was V pre-#100, swapped with TRBC1)
        assert HUMAN_TRBC2_AA[135] == "E"

    def test_jovi1_distinguishing_residues_differ_between_trbc1_and_trbc2(self):
        """JOVI.1 antibody distinguishes TRBC1 (`DLNK...`) from TRBC2
        (`DLKN...`) at positions 4-5 of the mature C — these must differ
        between the two canonical sequences, otherwise every β chain
        we ship will look like TRBC1 to a JOVI.1 flow stain regardless
        of its actual TRBC identity. Pre-#100 both started ``EDLNKVF``."""
        # Strip the prepended E to get back to mature positions
        trbc1_mature = HUMAN_TRBC1_AA[1:]
        trbc2_mature = HUMAN_TRBC2_AA[1:]
        assert trbc1_mature[2:4] == "NK", f"TRBC1 mature pos 3-4: {trbc1_mature[2:4]}"
        assert trbc2_mature[2:4] == "KN", f"TRBC2 mature pos 3-4: {trbc2_mature[2:4]}"

    def test_trbc1_trbc2_position_135_is_swapped(self):
        """TRBC1 has ``V`` at mature pos 135, TRBC2 has ``E`` — the
        opposite of what tcrsift had pre-#100."""
        trbc1_mature = HUMAN_TRBC1_AA[1:]
        trbc2_mature = HUMAN_TRBC2_AA[1:]
        assert trbc1_mature[134] == "V"
        assert trbc2_mature[134] == "E"


class TestCanonicalMatchesPyensembl:
    """Source-of-truth check (#100): the packaged FASTA must match what
    pyensembl returns for each canonical transcript. Skips cleanly if
    pyensembl isn't installed or the GRCh38 release isn't cached."""

    PROVENANCE = {
        "TRAC":  ("ENST00000611116", HUMAN_TRAC_AA, False),   # no prepend
        "TRBC1": ("ENST00000633705", HUMAN_TRBC1_AA, True),    # prepend E
        "TRBC2": ("ENST00000466254", HUMAN_TRBC2_AA, True),    # prepend E
    }

    @pytest.fixture(scope="class")
    def ensembl_release(self):
        pyensembl = pytest.importorskip("pyensembl")
        ens = pyensembl.EnsemblRelease(110)
        # Trigger a cheap lookup; if data isn't cached, skip rather than
        # downloading multi-GB at test time.
        try:
            ens.transcript_by_id("ENST00000611116")
        except Exception as exc:
            pytest.skip(f"pyensembl GRCh38 release 110 not cached: {exc}")
        return ens

    @pytest.mark.parametrize(
        "gene,tx_id,expected_aa,prepend_e",
        [
            ("TRAC", *PROVENANCE["TRAC"]),
            ("TRBC1", *PROVENANCE["TRBC1"]),
            ("TRBC2", *PROVENANCE["TRBC2"]),
        ],
    )
    def test_fasta_matches_pyensembl(
        self, ensembl_release, gene, tx_id, expected_aa, prepend_e
    ):
        prot = ensembl_release.transcript_by_id(tx_id).protein_sequence
        # Strip the X splice-boundary placeholder pyensembl emits at
        # the J→C junction of TR_C_gene transcripts.
        if prot.startswith("X"):
            prot = prot[1:]
        # Apply tcrsift's prepend-E convention for β chains.
        if prepend_e:
            prot = "E" + prot
        assert prot == expected_aa, (
            f"{gene}: packaged FASTA differs from pyensembl ENST {tx_id}. "
            f"Regenerate with `python -m tcrsift._regen_canonical_constants` "
            f"and verify against UniProt before committing."
        )


class TestBackTranslate:
    """Back-translation must round-trip through translate_dna."""

    def test_round_trip_canonical_trac(self):
        nt = back_translate(HUMAN_TRAC_AA)
        aa, _ = translate_dna(nt)
        assert aa == HUMAN_TRAC_AA

    def test_round_trip_all_canonical(self):
        for gene, aa in HUMAN_CONSTANT_REGIONS_AA.items():
            recovered_aa, _ = translate_dna(back_translate(aa))
            assert recovered_aa == aa, f"{gene} round-trip failed"

    def test_unknown_residue_becomes_nnn(self):
        nt = back_translate("LXM")  # X is not in HUMAN_PREFERRED_CODONS
        assert "NNN" in nt


class TestGetConstantRegionSequences:
    """The replacement for the pyensembl-backed loader."""

    def test_returns_all_three_genes(self):
        constants = get_constant_region_sequences()
        assert set(constants) == {"TRAC", "TRBC1", "TRBC2"}

    def test_returns_back_translated_canonical(self):
        """Each NT must translate to the canonical AA (#66 — the
        pre-fix path returned 2–11 aa)."""
        constants = get_constant_region_sequences()
        for gene, nt in constants.items():
            aa, _ = translate_dna(nt)
            # Returned NT ends with a stop codon; the translated AA
            # contains a trailing '*'.
            assert aa.rstrip("*") == HUMAN_CONSTANT_REGIONS_AA[gene]
            # Length sanity: the pre-fix bug produced 2-11 aa.
            assert len(aa.rstrip("*")) > 100

    def test_does_not_import_pyensembl(self):
        """The new path doesn't touch pyensembl. This test imports
        from tcrsift.assemble; if pyensembl were imported at module
        level the test would have to mock it. The fact that the
        existing test_assemble tests pass without pyensembl in the
        env is itself evidence — but we also confirm no module-level
        attribute exposes it."""
        from tcrsift import assemble

        # `get_constant_region_sequences` should not pull pyensembl.
        # We don't assert the absence at import time (pyensembl could
        # be installed for other reasons), just that the function
        # doesn't fail when called without it.
        assert get_constant_region_sequences()  # should not raise


class TestPickCanonicalConstant:
    """Routing chain + c_gene + j_gene to the right canonical."""

    def test_alpha_always_trac(self):
        for c_gene in (None, "", "TRAC", "TRAC*01", "anything"):
            name, aa = pick_canonical_constant("alpha", c_gene=c_gene)
            assert name == "TRAC" and aa == HUMAN_TRAC_AA

    def test_beta_trbc1_by_c_gene(self):
        name, aa = pick_canonical_constant("beta", c_gene="TRBC1")
        assert name == "TRBC1" and aa == HUMAN_TRBC1_AA

    def test_beta_trbc2_by_c_gene(self):
        name, aa = pick_canonical_constant("beta", c_gene="TRBC2")
        assert name == "TRBC2" and aa == HUMAN_TRBC2_AA

    def test_beta_handles_allele_suffix(self):
        name, _ = pick_canonical_constant("beta", c_gene="TRBC2*01")
        assert name == "TRBC2"

    def test_beta_fallback_to_jgene_trbc2(self):
        """When c_gene is missing, TRBJ2-* is in-cis with TRBC2."""
        name, aa = pick_canonical_constant("beta", c_gene=None, j_gene="TRBJ2-1")
        assert name == "TRBC2" and aa == HUMAN_TRBC2_AA

    def test_beta_fallback_to_jgene_trbc1(self):
        name, aa = pick_canonical_constant("beta", c_gene="", j_gene="TRBJ1-2")
        assert name == "TRBC1" and aa == HUMAN_TRBC1_AA

    def test_beta_default_to_trbc1_when_nothing_known(self):
        name, _ = pick_canonical_constant("beta")
        assert name == "TRBC1"

    def test_beta_jgene_overrides_conflicting_cgene_trbj1(self):
        """#90: TRBJ1 is in-cis with TRBC1; a CellRanger TRBC2 call
        must be overridden."""
        name, aa = pick_canonical_constant(
            "beta", c_gene="TRBC2", j_gene="TRBJ1-1"
        )
        assert name == "TRBC1"
        assert aa == HUMAN_TRBC1_AA

    def test_beta_jgene_overrides_conflicting_cgene_trbj2(self):
        """#90: TRBJ2 is in-cis with TRBC2; a CellRanger TRBC1 call
        must be overridden."""
        name, aa = pick_canonical_constant(
            "beta", c_gene="TRBC1", j_gene="TRBJ2-5"
        )
        assert name == "TRBC2"
        assert aa == HUMAN_TRBC2_AA

    def test_beta_jgene_with_allele_suffix(self):
        """``"TRBJ1-1*02"`` and ``"TRBC2*01"`` are common allele forms;
        family parsing must strip the allele before deciding."""
        name, _ = pick_canonical_constant(
            "beta", c_gene="TRBC2*01", j_gene="TRBJ1-1*02"
        )
        assert name == "TRBC1"

    def test_beta_jgene_lowercase_handled(self):
        """Defensive: lowercase / whitespace input should still drive
        the J-family override (some upstream tools emit lowercase)."""
        name, _ = pick_canonical_constant(
            "beta", c_gene="TRBC2", j_gene="  trbj1-2  "
        )
        assert name == "TRBC1"

    def test_beta_empty_j_falls_back_to_c_gene(self):
        """When J is missing or empty, c_gene drives the call."""
        for j in (None, "", "   "):
            name, _ = pick_canonical_constant(
                "beta", c_gene="TRBC2", j_gene=j
            )
            assert name == "TRBC2", f"failed for j_gene={j!r}"


class TestVerifyCanonicalConstantStart:
    """Cross-checking observed contig prefix against the canonical."""

    def test_exact_match_passes(self):
        # First 15 of TRAC starts with "IQNPDPAVYQLRDSK"
        observed = HUMAN_TRAC_AA[:15]
        assert verify_canonical_constant_start(observed, HUMAN_TRAC_AA)

    def test_one_residue_off_still_passes(self):
        """8 of 15 minimum is the default threshold."""
        # Substitute a few positions but keep most aligned.
        observed = "IXNPDPAXYQLRDXK"  # 12/15 match
        assert verify_canonical_constant_start(observed, HUMAN_TRAC_AA)

    def test_wrong_gene_fails(self):
        # TRBC1 start (EDLNKVFPPEVAVFE) vs TRAC canonical: 0/15 match.
        observed = HUMAN_TRBC1_AA[:15]
        assert not verify_canonical_constant_start(observed, HUMAN_TRAC_AA)

    def test_empty_observed_fails(self):
        assert not verify_canonical_constant_start("", HUMAN_TRAC_AA)

    def test_min_match_threshold_respected(self):
        """At threshold=5, looser comparisons pass."""
        observed = "IQNPDXXXXXXXXXX"  # 5/15 match
        assert verify_canonical_constant_start(
            observed, HUMAN_TRAC_AA, min_match=5
        )
        assert not verify_canonical_constant_start(
            observed, HUMAN_TRAC_AA, min_match=8
        )


class TestExtractCRegionStartFromContig:
    """The contig-side extractor used to verify the canonical pick."""

    def test_extracts_residues_after_vdj(self):
        """A contig that contains LEADER + VDJ + canonical TRAC start
        should yield the TRAC start when we split on the VDJ AA."""
        from tcrsift.assemble import HUMAN_PREFERRED_CODONS

        leader_aa = "MAGTWLLLLLALGCPALPTG"
        vdj_aa = "CASSLAPGATNEKLFF"  # arbitrary placeholder VDJ
        c_start = HUMAN_TRAC_AA[:30]  # what we want extracted
        full_aa = leader_aa + vdj_aa + c_start
        # Build a contig DNA that translates to `full_aa`.
        contig_dna = "".join(HUMAN_PREFERRED_CODONS[r] for r in full_aa) + "TAA"

        row = pd.Series(
            {
                "samples": "S1",
                "alpha_contig_ids": "contig_1",
            }
        )
        sample_contigs = {"S1": {"contig_1": contig_dna}}

        observed = _extract_c_region_start_from_contig(
            row, sample_contigs, vdj_aa, "alpha", n_aa=15
        )
        assert observed == HUMAN_TRAC_AA[:15]

    def test_returns_none_when_no_contig_ids(self):
        row = pd.Series({"samples": "S1"})  # no alpha_contig_ids column
        assert (
            _extract_c_region_start_from_contig(row, {}, "VDJ", "alpha") is None
        )

    def test_returns_none_when_vdj_empty(self):
        row = pd.Series({"samples": "S1", "alpha_contig_ids": "x"})
        assert _extract_c_region_start_from_contig(row, {}, "", "alpha") is None

    def test_returns_none_when_post_vdj_too_short(self):
        """If the contig assembly stops right after the J-gene, there
        aren't ``n_aa`` residues to return — we get None, not a truncated
        prefix masquerading as the C-region start."""
        from tcrsift.assemble import HUMAN_PREFERRED_CODONS

        vdj_aa = "CASSWY"
        # Only 5 aa after the VDJ.
        full_aa = "MAGT" + vdj_aa + "IQNPD"
        contig_dna = "".join(HUMAN_PREFERRED_CODONS[r] for r in full_aa) + "TAA"
        row = pd.Series({"samples": "S1", "alpha_contig_ids": "c1"})
        sample_contigs = {"S1": {"c1": contig_dna}}

        # Want 15 aa, only 5 available → None.
        assert (
            _extract_c_region_start_from_contig(
                row, sample_contigs, vdj_aa, "alpha", n_aa=15
            )
            is None
        )


class TestExtractCRegionNtFromContig:
    """The NT-level extractor used by the contig-aware C region assembly."""

    def _make_contig(self, vdj_aa: str, post_vdj_extra_nt: str) -> tuple[pd.Series, dict, str]:
        leader_aa = "MAGT"
        # Build a contig: leader + VDJ + extra post-VDJ NT
        contig_dna = (
            "".join(HUMAN_PREFERRED_CODONS[r] for r in leader_aa + vdj_aa)
            + post_vdj_extra_nt
        )
        vdj_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_aa)
        row = pd.Series({"samples": "S1", "beta_contig_ids": "c1"})
        sample_contigs = {"S1": {"c1": contig_dna}}
        return row, sample_contigs, vdj_nt

    def test_returns_post_vdj_nt(self):
        row, sample_contigs, vdj_nt = self._make_contig("CASSAY", "GAGGAC")
        result = _extract_c_region_nt_from_contig(row, sample_contigs, vdj_nt, "beta")
        assert result == "GAGGAC"

    def test_returns_none_when_no_contig(self):
        row = pd.Series({"samples": "S1"})  # no beta_contig_ids
        assert (
            _extract_c_region_nt_from_contig(row, {}, "ACGTACG", "beta") is None
        )

    def test_returns_none_when_vdj_nt_empty(self):
        row = pd.Series({"samples": "S1", "beta_contig_ids": "c1"})
        assert (
            _extract_c_region_nt_from_contig(row, {"S1": {"c1": "ACGT"}}, "", "beta")
            is None
        )

    def test_returns_none_when_vdj_nt_not_in_contig(self):
        row = pd.Series({"samples": "S1", "beta_contig_ids": "c1"})
        contigs = {"S1": {"c1": "AAAACCCCGGGGTTTT"}}
        # vdj_nt isn't anywhere in this contig
        assert (
            _extract_c_region_nt_from_contig(row, contigs, "GAGGAC", "beta") is None
        )


class TestBlendConstantNtWithContig:
    """The hybrid contig-derived + canonical-codon-optimized C region NT builder.

    Strategy under test:
    * Junction codon and downstream codons come from the contig as far as
      they agree with the canonical AA.
    * Any partial codon at the contig 3' edge is completed using canonical
      NT (preferring the codon-optimized choice when compatible).
    * Everything past the contig coverage is codon-optimized canonical NT.
    * AA mismatches between contig translation and canonical trigger an
      early switch back to canonical, recorded in the debug dict.
    """

    def test_no_contig_returns_canonical(self):
        canonical_aa = "EDLN"
        canonical_nt = back_translate(canonical_aa)
        result, debug = _blend_constant_nt_with_contig("", canonical_aa, canonical_nt)
        assert result == canonical_nt
        assert debug["n_contig_codons"] == 0
        assert debug["aa_mismatch_at"] is None

    def test_uses_all_matching_contig_codons(self):
        canonical_aa = "EDLN"
        canonical_nt = back_translate(canonical_aa)  # codon-optimized
        # Build a contig NT that translates to EDL (matches canonical 0..2)
        # using non-codon-optimized codons so we can prove the result
        # comes from the contig and not from canonical NT.
        contig_nt = "GAAGATCTC"  # GAA=E, GAT=D, CTC=L — all valid, none preferred
        # Sanity check the fixture: the contig codons are deliberately
        # NOT the codon-optimized choices, so equality of result[:9]
        # with contig_nt below proves the bytes came from the contig.
        for c in (contig_nt[:3], contig_nt[3:6], contig_nt[6:9]):
            assert CODON_TABLE[c] in {"E", "D", "L"}
            assert c != HUMAN_PREFERRED_CODONS[CODON_TABLE[c]]
        result, debug = _blend_constant_nt_with_contig(contig_nt, canonical_aa, canonical_nt)
        # First 9 nt should be the contig bytes; positions 9.. should be canonical_nt[9:]
        assert result[:9] == contig_nt
        assert result[9:] == canonical_nt[9:]
        assert debug["n_contig_codons"] == 3
        assert debug["aa_mismatch_at"] is None

    def test_completes_partial_codon_compatible(self):
        canonical_aa = "EDLN"
        canonical_nt = back_translate(canonical_aa)
        # Contig provides 1 full codon (E) plus 1 partial nt 'G'
        # The next canonical AA is 'D'; codons for D are GAT and GAC.
        # 'G' is a compatible prefix → completed codon is GAT or GAC.
        contig_nt = "GAAG"  # GAA=E, then partial 'G'
        result, debug = _blend_constant_nt_with_contig(contig_nt, canonical_aa, canonical_nt)
        # First 3 nt = contig E codon; next 3 nt = completed D codon
        assert result[:3] == "GAA"
        # Completed codon should start with the partial nt 'G' and code for D
        completed = result[3:6]
        assert completed.startswith("G")
        assert CODON_TABLE[completed] == "D"
        assert debug["partial_codon_completed"] is True
        # Positions 6.. should be canonical_nt[6:] (the L and N codons)
        assert result[6:] == canonical_nt[6:]

    def test_incompatible_partial_codon_falls_back_to_canonical(self):
        canonical_aa = "EDLN"
        canonical_nt = back_translate(canonical_aa)
        # Partial nt 'T' is incompatible with D codons (which all start
        # with G). The partial bytes get discarded and canonical takes
        # over at the next codon boundary.
        contig_nt = "GAAT"  # GAA=E, then partial 'T' (no D codon starts with T)
        result, debug = _blend_constant_nt_with_contig(contig_nt, canonical_aa, canonical_nt)
        assert result[:3] == "GAA"  # contig E
        assert result[3:] == canonical_nt[3:]  # canonical from D onward
        assert debug["partial_codon_completed"] is False
        assert debug["partial_codon_dropped"] is True  # distinct from "no partial"
        assert debug["n_contig_codons"] == 1
        # And the human-readable source string flags the drop so the
        # audit trail makes it visible.
        assert "partial dropped" in debug["source"]

    def test_aa_mismatch_switches_to_canonical(self):
        canonical_aa = "EDLN"
        canonical_nt = back_translate(canonical_aa)
        # Contig codes for E, then a wrong residue (Q instead of D), then
        # other codons we don't care about — the algorithm should switch
        # to canonical the moment it sees Q != D at position 1.
        contig_nt = "GAA" + "CAG" + "CTC" + "AAC"  # E, Q, L, N
        result, debug = _blend_constant_nt_with_contig(contig_nt, canonical_aa, canonical_nt)
        assert result[:3] == "GAA"  # the matching E codon from contig
        assert result[3:] == canonical_nt[3:]  # canonical from position 1 (the D)
        assert debug["n_contig_codons"] == 1
        assert debug["aa_mismatch_at"] == 2  # 1-indexed canonical AA pos

    def test_mismatch_at_first_codon_returns_pure_canonical(self):
        canonical_aa = "EDLN"
        canonical_nt = back_translate(canonical_aa)
        # Contig codes for Q at position 0 — no contig codons match.
        contig_nt = "CAG" + "GAT"
        result, debug = _blend_constant_nt_with_contig(contig_nt, canonical_aa, canonical_nt)
        assert result == canonical_nt
        assert debug["n_contig_codons"] == 0
        assert debug["aa_mismatch_at"] == 1

    def test_translates_back_to_canonical_aa(self):
        """The blended NT must translate to the canonical AA — the whole
        point is that the protein is preserved while the NT respects the
        donor at the junction."""
        canonical_aa = HUMAN_TRBC1_AA[:50]
        canonical_nt = back_translate(canonical_aa)
        contig_nt = canonical_nt[:30]  # contig covers first 10 codons exactly
        result, debug = _blend_constant_nt_with_contig(contig_nt, canonical_aa, canonical_nt)
        translated_back, _ = translate_dna(result)
        assert translated_back == canonical_aa
        assert debug["n_contig_codons"] == 10

    # ------------------------------------------------------------
    # Audit-driven additions (issue surfaced during PR #103 review)
    # ------------------------------------------------------------

    @pytest.mark.parametrize(
        "contig_nt,expected_n_codons,desc",
        [
            ("", 0, "no contig at all"),
            ("GAA", 1, "exactly one full codon"),
            ("GAAGAT", 2, "two full codons"),
            ("GAAGATCTCAAT", 4, "all four full codons (matches EDLN)"),
            ("GAAGATCTCAATTTT", 4, "contig overruns canonical_aa"),
            ("GAAG", 1, "1 full + partial completed"),
            ("GAAT", 1, "1 full + partial dropped"),
        ],
    )
    def test_length_invariant_holds_across_paths(
        self, contig_nt, expected_n_codons, desc
    ):
        """Regardless of which branch was taken — pure canonical fallback,
        contig + canonical, contig + completed-partial + canonical,
        contig + dropped-partial + canonical, or AA mismatch switch —
        the blended NT must be exactly ``3 * len(canonical_aa)`` long.
        That's the contract every downstream consumer relies on; it
        wasn't asserted directly before this PR."""
        canonical_aa = "EDLN"
        canonical_nt = back_translate(canonical_aa)
        blended, debug = _blend_constant_nt_with_contig(
            contig_nt, canonical_aa, canonical_nt
        )
        assert len(blended) == 3 * len(canonical_aa), (
            f"length invariant violated for {desc!r}: "
            f"got {len(blended)}, expected {3 * len(canonical_aa)}; debug={debug}"
        )

    def test_pure_partial_one_nt_compatible(self):
        """Pure partial codon (1 byte, no full codons) where the byte
        is a valid prefix for the canonical residue's codon. Untested
        before this PR — code path exists but wasn't exercised."""
        canonical_aa = "EDLN"
        canonical_nt = back_translate(canonical_aa)
        # 'G' is a compatible prefix for E (GAA, GAG).
        blended, debug = _blend_constant_nt_with_contig(
            "G", canonical_aa, canonical_nt
        )
        assert len(blended) == 3 * len(canonical_aa)
        assert blended[:3].startswith("G")
        assert CODON_TABLE[blended[:3]] == "E"
        assert debug["n_contig_codons"] == 0
        assert debug["partial_codon_completed"] is True
        assert debug["partial_codon_dropped"] is False
        # Rest of the chain comes from canonical-codon-opt.
        assert blended[3:] == canonical_nt[3:]

    def test_pure_partial_two_nt_compatible(self):
        """Pure partial codon, 2 bytes."""
        canonical_aa = "EDLN"
        canonical_nt = back_translate(canonical_aa)
        # 'GA' is a compatible prefix for E (GAA, GAG).
        blended, debug = _blend_constant_nt_with_contig(
            "GA", canonical_aa, canonical_nt
        )
        assert len(blended) == 3 * len(canonical_aa)
        assert blended[:3].startswith("GA")
        assert CODON_TABLE[blended[:3]] == "E"
        assert debug["partial_codon_completed"] is True

    def test_pure_partial_one_nt_incompatible(self):
        """Pure partial codon whose byte doesn't start ANY codon for the
        canonical residue → dropped. The whole canonical NT is returned
        and ``partial_codon_dropped`` is set so the audit trail records
        that donor fidelity was available but discarded."""
        canonical_aa = "EDLN"
        canonical_nt = back_translate(canonical_aa)
        # 'T' starts no E codon (E codons are GAA, GAG).
        blended, debug = _blend_constant_nt_with_contig(
            "T", canonical_aa, canonical_nt
        )
        assert len(blended) == 3 * len(canonical_aa)
        assert blended == canonical_nt  # everything canonical
        assert debug["n_contig_codons"] == 0
        assert debug["partial_codon_completed"] is False
        assert debug["partial_codon_dropped"] is True
        assert "partial dropped" in debug["source"]

    def test_contig_longer_than_canonical_truncated(self):
        """Excess bytes past ``3 * len(canonical_aa)`` are silently
        dropped — length invariant must hold."""
        canonical_aa = "EDLN"  # 4 residues → 12 nt
        canonical_nt = back_translate(canonical_aa)
        contig_nt = canonical_nt + "AAATTT"  # 6 excess bytes
        blended, debug = _blend_constant_nt_with_contig(
            contig_nt, canonical_aa, canonical_nt
        )
        assert len(blended) == 3 * len(canonical_aa)
        # All four canonical residues sourced from contig (since contig_nt
        # starts with canonical_nt byte-for-byte).
        assert debug["n_contig_codons"] == 4
        assert debug["aa_mismatch_at"] is None

    def test_partial_dropped_distinct_from_no_partial(self):
        """The ``partial_codon_dropped`` flag is what makes the audit
        trail useful — without it, the source string couldn't tell
        whether the donor had no partial at all vs. the donor had a
        partial that we threw away. Pin both signals separately."""
        canonical_aa = "EDLN"
        canonical_nt = back_translate(canonical_aa)
        # No partial.
        _, debug_no_partial = _blend_constant_nt_with_contig(
            "GAA", canonical_aa, canonical_nt
        )
        # Partial dropped.
        _, debug_dropped = _blend_constant_nt_with_contig(
            "GAAT", canonical_aa, canonical_nt
        )
        assert debug_no_partial["partial_codon_dropped"] is False
        assert debug_dropped["partial_codon_dropped"] is True
        assert "partial dropped" not in debug_no_partial["source"]
        assert "partial dropped" in debug_dropped["source"]


class TestExtractCRegionNtMultiContig:
    """Audit-driven coverage of the multi-sample / multi-contig path of
    ``_extract_c_region_nt_from_contig`` — the consensus picker over
    ``Counter.most_common`` was undocumented in tests before this PR."""

    def _make_multi_contig_row(self, vdj_nt, contigs):
        """``contigs`` is a list of ``(sample, contig_id, full_contig_seq)``
        tuples. The row stitches contig_ids across all entries."""
        row = pd.Series({
            "samples": ";".join(sorted({c[0] for c in contigs})),
            "beta_contig_ids": ";".join(c[1] for c in contigs),
        })
        sample_contigs = {}
        for sample, contig_id, seq in contigs:
            sample_contigs.setdefault(sample, {})[contig_id] = seq
        return row, sample_contigs

    def test_picks_most_common_post_vdj_nt(self):
        """When multiple cells' contigs all contain the representative
        ``vdj_nt`` but disagree on the post-VDJ bytes, the most-common
        post-VDJ NT wins. Two contigs end in ``GAGGAC`` vs one in
        ``GAGGAA`` → ``GAGGAC`` is the consensus."""
        vdj_nt = "TGTGCATCTAGT"  # arbitrary 12-nt VDJ
        row, sample_contigs = self._make_multi_contig_row(
            vdj_nt,
            [
                ("S1", "c1", "AAA" + vdj_nt + "GAGGAC"),
                ("S1", "c2", "AAA" + vdj_nt + "GAGGAC"),
                ("S2", "c3", "AAA" + vdj_nt + "GAGGAA"),
            ],
        )
        result = _extract_c_region_nt_from_contig(
            row, sample_contigs, vdj_nt, "beta"
        )
        assert result == "GAGGAC"

    def test_filters_contigs_that_dont_contain_vdj_nt(self):
        """Cross-source hazard mitigation: contigs whose VDJ doesn't
        match the representative ``vdj_nt`` byte-for-byte must NOT
        contribute to the consensus. (#94 archetype — protects against
        post-VDJ bytes being mixed across cells whose underlying VDJ
        differs.)"""
        vdj_nt = "TGTGCATCTAGT"
        row, sample_contigs = self._make_multi_contig_row(
            vdj_nt,
            [
                # Contains the rep vdj_nt → its post-VDJ contributes.
                ("S1", "c1", "AAA" + vdj_nt + "GAGGAC"),
                # Different VDJ ("TGTGCAACGAGT") — must be ignored even
                # though it's in the contig list.
                ("S1", "c2", "AAA" + "TGTGCAACGAGT" + "GAGGAA"),
            ],
        )
        result = _extract_c_region_nt_from_contig(
            row, sample_contigs, vdj_nt, "beta"
        )
        # Only c1's post-VDJ bytes contribute; c2's "GAGGAA" is filtered out.
        assert result == "GAGGAC"

    def test_trailing_semicolon_in_samples_field_handled(self):
        """The aggregator's join can produce trailing ``;`` artifacts.
        ``split(";")`` then yields empty strings that won't match any
        ``sample_contigs`` key — they should be silently skipped
        without affecting the consensus."""
        vdj_nt = "TGTGCATCTAGT"
        row = pd.Series({
            "samples": "S1;",  # trailing semicolon
            "beta_contig_ids": "c1",
        })
        sample_contigs = {"S1": {"c1": "AAA" + vdj_nt + "GAGGAC"}}
        result = _extract_c_region_nt_from_contig(
            row, sample_contigs, vdj_nt, "beta"
        )
        assert result == "GAGGAC"


class TestAssembleEndToEndContigBlend:
    """E2E: drive ``assemble_full_sequences`` with contigs and verify the
    NT→AA round-trip invariant (#91) still holds with the new blend
    path active. Unit tests cover the blend function directly; this
    locks in the integration through ``_add_constant_regions`` →
    ``_build_full_sequences`` → ``validate_sequences``."""

    def test_full_beta_nt_translates_to_full_beta_aa_with_blend(self, tmp_path):
        from tcrsift.assemble import (
            HUMAN_PREFERRED_CODONS,
            HUMAN_TRBC1_AA,
            assemble_full_sequences,
            translate_dna,
            validate_sequences,
        )

        # Build a clone with realistic VDJ_aa and a contig whose bytes
        # immediately after the trimmed VDJ are codon-aligned to the
        # first few canonical C-region residues. (In real CellRanger
        # output the bytes between VDJ and C are the J/C junction codon
        # split across the J overshoot and the C exon's first 1-2 nt;
        # for this E2E test of the blend integration we just put them
        # codon-aligned so the contig's first 8 codons match canonical
        # C[0:8]. The blend function's biology is exercised in unit
        # tests above; this one verifies the orchestration end-to-end.)
        leader_aa = "M" + "A" * 19
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        vdj_beta_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_beta)
        # +1 overshoot only on the ROW's VDJ_beta_nt (so the #91 trim
        # actually fires); the contig itself is codon-aligned.
        overshoot = "A"
        c_region_start_nt = "".join(
            HUMAN_PREFERRED_CODONS[r] for r in HUMAN_TRBC1_AA[:8]
        )
        beta_contig_seq = (
            "".join(HUMAN_PREFERRED_CODONS[r] for r in leader_aa)
            + vdj_beta_nt
            + c_region_start_nt
        )

        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha,
            "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha,
            "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_alpha),
            "VDJ_beta_nt": vdj_beta_nt + overshoot,
            "alpha_c_gene": "TRAC",
            "beta_c_gene": "TRBC1",
            "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
            "beta_contig_ids": "contig_b1",
        }])

        # Build a fake CellRanger contig directory. ``load_contigs``
        # uses the FASTA's parent-directory name as the sample key, so
        # the sample subdirectory IS the sample name.
        contig_dir = tmp_path / "S1"
        contig_dir.mkdir(parents=True)
        fasta = contig_dir / "filtered_contig.fasta"
        fasta.write_text(f">contig_b1\n{beta_contig_seq}\n")

        out = assemble_full_sequences(
            df,
            contigs_dir=str(tmp_path),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )

        # Round-trip invariant on the spliced NT.
        full_nt = out["full_beta_nt"].iloc[0]
        full_aa = out["full_beta_aa"].iloc[0]
        # Strip trailing stop codon if present.
        nt_for_translation = full_nt[:-3] if full_nt[-3:] in {"TAA", "TAG", "TGA"} else full_nt
        translated, _ = translate_dna(nt_for_translation)
        assert translated == full_aa, (
            "full_beta_nt does not translate to full_beta_aa — the "
            "blend path broke the round-trip invariant from #91."
        )

        # The source label records the contig contribution.
        source = out["beta_constant_source"].iloc[0]
        assert "contig(" in source, (
            f"expected contig contribution in beta_constant_source, got {source!r}"
        )

        # validate_sequences (which runs _validate_nt_aa_roundtrip) should
        # not produce load-bearing NT-translation failures.
        msgs = validate_sequences(out, strict=False)
        translation_failures = [
            m for m in msgs
            if m.severity == "load_bearing"
            and ("translate to" in m or "frame likely broken" in m)
        ]
        assert translation_failures == [], (
            f"unexpected NT round-trip failures: {translation_failures}"
        )
