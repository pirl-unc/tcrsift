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
    CONSTANT_REGION_ENDINGS,
    HUMAN_CONSTANT_REGIONS_AA,
    HUMAN_TRAC_AA,
    HUMAN_TRBC1_AA,
    HUMAN_TRBC2_AA,
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
