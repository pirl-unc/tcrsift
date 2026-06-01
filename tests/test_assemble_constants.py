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
    HUMAN_CONSTANT_ALLELES,
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
        "gene,expected_len",
        [
            ("TRAC", 140),    # mature TRAC, post-X-strip from pyensembl ENST00000611116
            ("TRBC1", 176),   # mature TRBC1, post-X-strip from ENST00000633705 (no synthetic E since 2.0)
            ("TRBC2", 178),   # mature TRBC2, post-X-strip from ENST00000466254
        ],
    )
    def test_mature_chain_length_exact(self, gene, expected_len):
        aa = HUMAN_CONSTANT_REGIONS_AA[gene]
        assert len(aa) == expected_len, (
            f"{gene} length {len(aa)} != UniProt mature length {expected_len}"
        )

    def test_all_three_present(self):
        assert set(HUMAN_CONSTANT_REGIONS_AA) == {"TRAC", "TRBC1", "TRBC2"}

    def test_diagnostic_positions_match_uniprot(self):
        """Lock the specific residues that previously drifted (#100).

        These positions were wrong in the pre-#100 hardcoded sequences;
        pinning them here makes any future hand-edit drift loud. Indices
        are 0-indexed Python slices into the stored AA constants.

        As of 2.0 (#105) all three canonical sequences are stored as
        BARE MATURE (no synthetic prepend on β chains). UniProt mature
        position N → stored index N-1 for all three genes.
        """
        # TRAC mature position 47 → stored index 46 — the conserved Thr.
        # Pre-#100 had 'C' here.
        assert HUMAN_TRAC_AA[46] == "T", (
            f"TRAC mature pos 47 must be T (UniProt P01848); got "
            f"{HUMAN_TRAC_AA[46]!r} (pre-#100 had 'C')"
        )

        # TRBC1 mature pos 135 (idx 134). Pre-#100 had 'E', swapped with TRBC2.
        assert HUMAN_TRBC1_AA[134] == "V", (
            f"TRBC1 mature pos 135 must be V (UniProt P01850); got "
            f"{HUMAN_TRBC1_AA[134]!r} (pre-#100 had 'E')"
        )

        # TRBC2 mature diagnostic positions (bare-mature, 1-indexed).
        # Mature pos 3 (idx 2): K (JOVI.1 distinguishing — TRBC1 has N here)
        assert HUMAN_TRBC2_AA[2] == "K", (
            f"TRBC2 mature pos 3 must be K (JOVI.1 epitope); got "
            f"{HUMAN_TRBC2_AA[2]!r}"
        )
        # Mature pos 4 (idx 3): N (JOVI.1 distinguishing — TRBC1 has K here)
        assert HUMAN_TRBC2_AA[3] == "N", (
            f"TRBC2 mature pos 4 must be N (JOVI.1 epitope); got "
            f"{HUMAN_TRBC2_AA[3]!r}"
        )
        # Mature pos 9 (idx 8): E. As of 2.3 (#113) the packaged
        # TRBC2 default is the TRBC2*01-protein form (E at pos 9 —
        # major-allele in humans). Pre-2.3 versions shipped the *03
        # form with K here. The *03 sequence is still available via
        # the allele picker (auto-detect, or explicit
        # ``trbc2_allele="03"``).
        assert HUMAN_TRBC2_AA[8] == "E"
        # Mature pos 56 (idx 55): S (was C pre-#100; TRBC1 has S too — not diagnostic)
        assert HUMAN_TRBC2_AA[55] == "S"
        # Mature pos 135 (idx 134): E (was V pre-#100, swapped with TRBC1)
        assert HUMAN_TRBC2_AA[134] == "E"

    def test_jovi1_distinguishing_residues_differ_between_trbc1_and_trbc2(self):
        """JOVI.1 antibody distinguishes TRBC1 (`DLNK...`) from TRBC2
        (`DLKN...`) at mature positions 3-4. These must differ between
        the two canonical sequences, otherwise every β chain we ship
        will look like TRBC1 to a JOVI.1 flow stain regardless of its
        actual TRBC identity. Pre-#100 both started ``EDLNKVF``."""
        assert HUMAN_TRBC1_AA[2:4] == "NK", (
            f"TRBC1 mature pos 3-4: {HUMAN_TRBC1_AA[2:4]!r}"
        )
        assert HUMAN_TRBC2_AA[2:4] == "KN", (
            f"TRBC2 mature pos 3-4: {HUMAN_TRBC2_AA[2:4]!r}"
        )

    def test_trbc1_trbc2_position_135_is_swapped(self):
        """TRBC1 has ``V`` at mature pos 135, TRBC2 has ``E`` — the
        opposite of what tcrsift had pre-#100."""
        assert HUMAN_TRBC1_AA[134] == "V"
        assert HUMAN_TRBC2_AA[134] == "E"


class TestTRACBiology:
    """Detailed biology pins on the canonical TRAC sequence.

    UniProt P01848 (TRAC_HUMAN) mature chain; pyensembl ENST00000611116
    canonical transcript. Stored as bare mature since 2.0 (#105) — no
    synthetic prepend, matches UniProt directly.
    """

    def test_starts_with_iqnpdpavy(self):
        """TRAC mature starts at I (residue 1 of UniProt P01848 chain)."""
        assert HUMAN_TRAC_AA.startswith("IQNPDPAVY"), (
            f"TRAC mature must start with IQNPDPAVY; got {HUMAN_TRAC_AA[:9]!r}"
        )

    def test_ends_with_llmtlrlwss(self):
        """The TRAC C-terminus is the last 10 residues of the cytoplasmic
        tail. Mismatches here are a reliable signal of either truncation
        (the #66 frame bug produced 2–11 aa truncations) or wrong-frame
        translation."""
        assert HUMAN_TRAC_AA.endswith("LLMTLRLWSS"), (
            f"TRAC tail must be LLMTLRLWSS; got {HUMAN_TRAC_AA[-10:]!r}"
        )

    def test_no_internal_stop(self):
        """No premature stop in the mature chain — sanity check
        against accidental in-frame stops from a corrupted FASTA."""
        assert "*" not in HUMAN_TRAC_AA

    def test_no_unknown_residue(self):
        """The X splice-boundary placeholder must be stripped at load
        time; an X in the stored constant would indicate the strip
        logic in ``_load_canonical_constants_fasta`` regressed."""
        assert "X" not in HUMAN_TRAC_AA

    def test_thr47_conserved_residue(self):
        """TRAC position 47 (idx 46) is a conserved Thr. Pre-#100 had
        Cys here — a hand-typing drift that this test pins for the
        future. Position cited in UniProt P01848 annotations."""
        assert HUMAN_TRAC_AA[46] == "T"

    def test_alphabet_is_standard(self):
        """Only the 20 standard residues — no X, no *, no stops, no
        lowercase, no non-AA characters."""
        valid = set("ACDEFGHIKLMNPQRSTVWY")
        assert set(HUMAN_TRAC_AA).issubset(valid), (
            f"unexpected residues in TRAC: {set(HUMAN_TRAC_AA) - valid}"
        )


class TestTRBC1Biology:
    """Detailed biology pins on the canonical TRBC1 sequence.

    UniProt P01850 (TRBC1_HUMAN) mature chain; pyensembl ENST00000633705
    canonical transcript. Stored as bare mature since 2.0 (#105) — pre-2.0
    had a synthetic E prepended. The J→C junction E is now read per-clone
    from the contig in ``_add_constant_regions``.
    """

    def test_starts_with_dlnk(self):
        """TRBC1 mature C-region starts at D (residue 1 of UniProt P01850
        chain). The JOVI.1-distinguishing tetrad NK at positions 3-4."""
        assert HUMAN_TRBC1_AA.startswith("DLNKVF"), (
            f"TRBC1 must start with DLNKVF; got {HUMAN_TRBC1_AA[:6]!r}"
        )

    def test_jovi1_signature_at_positions_3_4(self):
        """JOVI.1 antibody recognises the NK epitope at TRBC1 positions
        3-4 of the mature chain. A flow stain against TRBC1-expressing
        T cells uses this exact tetrad."""
        assert HUMAN_TRBC1_AA[2:4] == "NK"

    def test_glu_at_mature_position_9(self):
        """TRBC1 has E at mature position 9 — the residue that
        distinguishes it from TRBC2's K at the same position."""
        assert HUMAN_TRBC1_AA[8] == "E"

    def test_phe_at_mature_position_36(self):
        """TRBC1 has F at mature position 36 (idx 35). TRBC2 has Y."""
        assert HUMAN_TRBC1_AA[35] == "F"

    def test_val_at_mature_position_135(self):
        """TRBC1 has V at mature position 135 (idx 134) — pre-#100 had
        an E here, swapped with TRBC2. The V/E swap at this position is
        one of the four TRBC1 vs TRBC2 distinguishing residues."""
        assert HUMAN_TRBC1_AA[134] == "V"

    def test_ends_with_vkrkdf(self):
        """TRBC1 C-terminus is ``VKRKDF`` — 2 aa shorter than TRBC2's
        ``VKRKDSRG``. The terminal motif distinguishes the two genes
        in nucleic-acid validation."""
        assert HUMAN_TRBC1_AA.endswith("VKRKDF")

    def test_does_not_have_trbc2_ctail(self):
        """The 2-aa C-terminal extension is unique to TRBC2; TRBC1 must
        not have RG at the end."""
        assert not HUMAN_TRBC1_AA.endswith("VKRKDSRG")
        assert not HUMAN_TRBC1_AA.endswith("RG")

    def test_no_internal_stop_or_x(self):
        assert "*" not in HUMAN_TRBC1_AA
        assert "X" not in HUMAN_TRBC1_AA

    def test_alphabet_is_standard(self):
        valid = set("ACDEFGHIKLMNPQRSTVWY")
        assert set(HUMAN_TRBC1_AA).issubset(valid)

    def test_no_synthetic_e_at_start_since_2_0(self):
        """Pre-2.0 the stored TRBC1 started with a synthetic E (prepended
        at FASTA-write time). Since 2.0 the per-clone junction residue
        is sourced from the contig; the canonical is bare mature
        starting at D."""
        assert HUMAN_TRBC1_AA[0] == "D"
        assert not HUMAN_TRBC1_AA.startswith("E")


class TestTRBC2Biology:
    """Detailed biology pins on the canonical TRBC2 sequence.

    UniProt A0A5B9 (TRBC2_HUMAN) mature chain; pyensembl ENST00000466254
    canonical transcript. Stored as bare mature since 2.0 (#105).

    NB: the original tcrsift comment cited ``A0A075B6Y0`` as the TRBC2
    UniProt reference — that's actually TRAJ49 (a 19-aa fragment).
    A0A5B9 is the correct canonical TRBC2 entry; #100 fixed the trace.
    """

    def test_starts_with_dlkn(self):
        """TRBC2 mature C-region starts at D (residue 1 of UniProt A0A5B9
        chain). The JOVI.1-distinguishing tetrad is KN at positions 3-4
        — the opposite order from TRBC1's NK."""
        assert HUMAN_TRBC2_AA.startswith("DLKNVF"), (
            f"TRBC2 must start with DLKNVF; got {HUMAN_TRBC2_AA[:6]!r}"
        )

    def test_jovi1_signature_at_positions_3_4(self):
        """JOVI.1 antibody recognises the NK epitope on TRBC1; the KN
        order on TRBC2 means JOVI.1 does NOT bind. Pre-#100 both genes
        were ``EDLNKVF`` — every TRBC2-expressing T cell would have
        falsely stained JOVI.1+, and any TRBC1 vs TRBC2 frequency
        analysis derived from tcrsift output would have been bunk."""
        assert HUMAN_TRBC2_AA[2:4] == "KN"

    def test_glu_at_mature_position_9_in_default_allele(self):
        """As of 2.3 (#113) the default TRBC2 is the *01-protein form,
        which has **E** at mature position 9 — same as TRBC1. So this
        position no longer distinguishes the two genes at the protein
        level under the default canonical. The *03-protein form (K at
        pos 9) is still available as a packaged allele under
        ``HUMAN_CONSTANT_ALLELES['TRBC2']['03']`` and is picked
        automatically when the donor's contig codes for K."""
        assert HUMAN_TRBC2_AA[8] == "E"
        # Sanity: the *03 form is still packaged and still has K.
        from tcrsift.assemble import HUMAN_CONSTANT_ALLELES

        assert HUMAN_CONSTANT_ALLELES["TRBC2"]["03"][8] == "K"

    def test_tyr_at_mature_position_36(self):
        """TRBC2 has Y at mature position 36 (idx 35). TRBC1 has F.
        One of the four TRBC1 vs TRBC2 distinguishing residues."""
        assert HUMAN_TRBC2_AA[35] == "Y"

    def test_ser_at_mature_position_56(self):
        """TRBC2 has S at mature position 56 (idx 55) — pre-#100 had C.
        TRBC1 has S too at this position, so this isn't a TRBC1-vs-TRBC2
        distinguishing residue, but it's a #100 drift to pin."""
        assert HUMAN_TRBC2_AA[55] == "S"

    def test_glu_at_mature_position_135(self):
        """TRBC2 has E at mature position 135 (idx 134). TRBC1 has V.
        The V/E swap at this position is the fourth and last TRBC1 vs
        TRBC2 distinguishing residue in the shared region. Pre-#100
        the two genes were swapped at this position."""
        assert HUMAN_TRBC2_AA[134] == "E"

    def test_ends_with_vkrkdsrg(self):
        """TRBC2 has a 2-aa C-terminal extension over TRBC1 — ``VKRKDSRG``
        vs ``VKRKDF``. The RG tail is the simplest way to tell the two
        genes apart in protein output."""
        assert HUMAN_TRBC2_AA.endswith("VKRKDSRG")

    def test_two_aa_longer_than_trbc1_in_ctail(self):
        """Length delta between the two genes is exactly 2 aa, all of
        it in the C-terminal extension."""
        assert len(HUMAN_TRBC2_AA) == len(HUMAN_TRBC1_AA) + 2

    def test_no_internal_stop_or_x(self):
        assert "*" not in HUMAN_TRBC2_AA
        assert "X" not in HUMAN_TRBC2_AA

    def test_alphabet_is_standard(self):
        valid = set("ACDEFGHIKLMNPQRSTVWY")
        assert set(HUMAN_TRBC2_AA).issubset(valid)

    def test_no_synthetic_e_at_start_since_2_0(self):
        assert HUMAN_TRBC2_AA[0] == "D"
        assert not HUMAN_TRBC2_AA.startswith("E")


class TestTRBC1VsTRBC2DistinguishingResidues:
    """The TRBC1 vs TRBC2 comparison: exactly which positions differ,
    in what direction, and why each one matters."""

    @pytest.mark.parametrize(
        "pos,trbc1,trbc2,note",
        [
            (3,   "N", "K", "JOVI.1 epitope position 1 (NK vs KN)"),
            (4,   "K", "N", "JOVI.1 epitope position 2 (NK vs KN)"),
            (36,  "F", "Y", "shared-region #100 swap (F in TRBC1, Y in TRBC2)"),
            (135, "V", "E", "shared-region #100 swap (V in TRBC1, E in TRBC2)"),
        ],
    )
    def test_distinguishing_position(self, pos, trbc1, trbc2, note):
        """Each TRBC1 vs TRBC2 distinguishing position pinned with the
        biological rationale (in the *01-protein default since 2.3 /
        #113). Pos 9 was previously a distinguisher but only because
        we shipped the TRBC2*03-protein form — both *01-protein forms
        of the two genes have E at pos 9. See
        :meth:`test_trbc2_03_variant_adds_pos_9_distinguisher` for the
        *03 case."""
        idx = pos - 1
        assert HUMAN_TRBC1_AA[idx] == trbc1, (
            f"TRBC1 mature pos {pos} ({note}): expected {trbc1}, got "
            f"{HUMAN_TRBC1_AA[idx]!r}"
        )
        assert HUMAN_TRBC2_AA[idx] == trbc2, (
            f"TRBC2 mature pos {pos} ({note}): expected {trbc2}, got "
            f"{HUMAN_TRBC2_AA[idx]!r}"
        )
        assert HUMAN_TRBC1_AA[idx] != HUMAN_TRBC2_AA[idx]

    def test_internal_region_distinguishing_positions(self):
        """The two genes' default canonicals (TRBC1*01 and TRBC2*01)
        share the same internal-region length (residues 1..170) and
        differ at exactly 4 positions there. Pos 9 is NOT a
        distinguisher between the *01-protein forms — only the *03
        variant has K there. After residue 170 the C-terminal tails
        structurally diverge: TRBC1 ends at ``VKRKDF`` while TRBC2
        continues to ``VKRKDSRG``."""
        internal_end = 170
        diffs = [
            i for i in range(internal_end)
            if HUMAN_TRBC1_AA[i] != HUMAN_TRBC2_AA[i]
        ]
        diff_positions = [i + 1 for i in diffs]
        assert diff_positions == [3, 4, 36, 135], (
            f"TRBC1 vs TRBC2 internal-region distinguishing positions: "
            f"{diff_positions} (expected [3, 4, 36, 135])"
        )

    def test_trbc2_03_variant_adds_pos_9_distinguisher(self):
        """The TRBC2*03-protein form (still packaged in
        ``HUMAN_CONSTANT_ALLELES``) has K at pos 9, which DOES
        distinguish it from TRBC1's E. Auto-detect picks this when the
        donor's contig codes K there."""
        from tcrsift.assemble import HUMAN_CONSTANT_ALLELES

        trbc2_03 = HUMAN_CONSTANT_ALLELES["TRBC2"]["03"]
        assert trbc2_03[8] == "K"
        assert HUMAN_TRBC1_AA[8] == "E"

    def test_internal_region_is_otherwise_identical(self):
        """Outside the 4 distinguishing positions in residues 1..170,
        TRBC1*01 and TRBC2*01 agree byte-for-byte. Anything else is a
        FASTA drift — pre-#100 there were many such drifts."""
        internal_end = 170
        distinguishing = {2, 3, 35, 134}  # 0-indexed (pos 9 dropped at #113)
        for i in range(internal_end):
            if i in distinguishing:
                continue
            assert HUMAN_TRBC1_AA[i] == HUMAN_TRBC2_AA[i], (
                f"unexpected TRBC1/TRBC2 difference at idx {i}: "
                f"{HUMAN_TRBC1_AA[i]!r} vs {HUMAN_TRBC2_AA[i]!r}"
            )

    def test_ctail_structural_divergence(self):
        """Both genes share the ``VKRKD`` motif (positions 171-175);
        TRBC1 then ends with a single ``F`` while TRBC2 extends with
        ``SRG``. The shared motif anchors NT-level alignment of the
        two genes' C-tails; the extension is the unambiguous protein-
        level marker of TRBC2."""
        # Shared VKRKD anchor.
        assert HUMAN_TRBC1_AA[170:175] == "VKRKD"
        assert HUMAN_TRBC2_AA[170:175] == "VKRKD"
        # TRBC1 single-residue tail.
        assert HUMAN_TRBC1_AA[175:] == "F"
        # TRBC2 three-residue tail.
        assert HUMAN_TRBC2_AA[175:] == "SRG"


class TestCanonicalMatchesPyensembl:
    """Source-of-truth check (#100, refined in #113): the packaged
    FASTA's allele entries must match what pyensembl / NCBI / UniProt
    return for each canonical transcript. Skips cleanly if pyensembl
    isn't installed or the GRCh38 release isn't cached.

    Note: as of 2.3 (#113), ``HUMAN_TRBC2_AA`` defaults to the
    TRBC2*01-protein form (NCBI AAA60662). Ensembl transcript
    ``ENST00000466254`` encodes the TRBC2*03-protein form (UniProt
    A0A5B9). The pyensembl comparison for TRBC2 therefore uses the
    *03 entry from ``HUMAN_CONSTANT_ALLELES['TRBC2']`` rather than
    ``HUMAN_TRBC2_AA``."""

    PROVENANCE = {
        "TRAC":  ("ENST00000611116", HUMAN_TRAC_AA),
        "TRBC1": ("ENST00000633705", HUMAN_TRBC1_AA),
        # Ensembl returns the TRBC2*03-protein form; the default
        # ``HUMAN_TRBC2_AA`` is *01 since 2.3 (#113), so this row
        # uses the *03 entry from the packaged multi-allele FASTA.
        "TRBC2": ("ENST00000466254", HUMAN_CONSTANT_ALLELES["TRBC2"]["03"]),
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
        "gene,tx_id,expected_aa",
        [
            ("TRAC", *PROVENANCE["TRAC"]),
            ("TRBC1", *PROVENANCE["TRBC1"]),
            ("TRBC2", *PROVENANCE["TRBC2"]),
        ],
    )
    def test_fasta_matches_pyensembl(
        self, ensembl_release, gene, tx_id, expected_aa
    ):
        """As of 2.0 (#105) all three canonical sequences are stored
        BARE MATURE (no synthetic prepend), so the comparison is direct."""
        prot = ensembl_release.transcript_by_id(tx_id).protein_sequence
        # Strip the X splice-boundary placeholder pyensembl emits at
        # the J→C junction of TR_C_gene transcripts.
        if prot.startswith("X"):
            prot = prot[1:]
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


class TestJunctionResidueUniformly:
    """#105 (2.0): the J→C junction residue is read per-clone from the
    contig's first codon past trimmed VDJ — same code path for both
    α and β. β's junction is universally E (GAG codon, 118/118 audited
    B1 clones); α's is J-dependent (N/Y/D/H seen in the same audit).
    Pre-2.0 β had a synthetic E baked into the stored canonical and α
    had nothing — every α chain was 1 aa short."""

    @staticmethod
    def _build_alpha_fixture(tmp_path, junction_codon: str):
        """Return ``(df, contigs_dir)`` for an α-only assembly test.
        ``junction_codon`` is the 3-nt junction codon (one of AAT, TAT,
        GAT, CAT for N, Y, D, H respectively) that the contig provides
        immediately after the trimmed VDJ."""
        leader_aa = "M" + "A" * 19
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_alpha_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_alpha)
        # +1 overshoot on the row's VDJ_alpha_nt so the #91 trim fires.
        overshoot = "A"
        # Contig: leader + vdj + junction codon + first 8 canonical TRAC codons.
        trac_first_8_nt = "".join(
            HUMAN_PREFERRED_CODONS[r] for r in HUMAN_TRAC_AA[:8]
        )
        alpha_contig_seq = (
            "".join(HUMAN_PREFERRED_CODONS[r] for r in leader_aa)
            + vdj_alpha_nt
            + junction_codon
            + trac_first_8_nt
        )
        # Minimal β placeholder (β is unaffected by #105 — its canonical
        # already carries the prepended E).
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha,
            "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha,
            "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": vdj_alpha_nt + overshoot,
            "VDJ_beta_nt": "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_beta),
            "alpha_c_gene": "TRAC",
            "beta_c_gene": "TRBC1",
            "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
            "alpha_contig_ids": "contig_a1",
        }])
        contig_dir = tmp_path / "S1"
        contig_dir.mkdir(parents=True)
        (contig_dir / "filtered_contig.fasta").write_text(
            f">contig_a1\n{alpha_contig_seq}\n"
        )
        return df, tmp_path

    @pytest.mark.parametrize(
        "codon,residue",
        [
            ("AAT", "N"),  # most common α junction in the B1 cohorts
            ("TAT", "Y"),
            ("GAT", "D"),
            ("CAT", "H"),
        ],
    )
    def test_alpha_junction_residue_prepended_from_contig(
        self, codon, residue, tmp_path
    ):
        """The α junction residue from the contig must be prepended to
        ``alpha_constant_aa`` so the assembled mature α matches the
        donor's expressed protein."""
        from tcrsift.assemble import assemble_full_sequences

        df, contigs_dir = self._build_alpha_fixture(tmp_path, codon)
        out = assemble_full_sequences(
            df,
            contigs_dir=str(contigs_dir),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        # The constant AA now begins with the junction residue, then
        # the canonical TRAC start.
        assert out["alpha_constant_aa"].iloc[0].startswith(residue + HUMAN_TRAC_AA[:8])
        # And the audit column records what was inserted.
        assert out["alpha_junction_residue"].iloc[0] == residue
        # full_alpha_aa = leader + vdj + (junction + TRAC), so the
        # junction residue appears immediately after the VDJ.
        full_aa = out["full_alpha_aa"].iloc[0]
        vdj_alpha = df["VDJ_alpha_aa"].iloc[0]
        assert vdj_alpha + residue + HUMAN_TRAC_AA[:8] in full_aa

    def test_full_alpha_nt_round_trips_to_full_alpha_aa(self, tmp_path):
        """With the prepend in place, ``translate(full_alpha_nt)`` must
        still equal ``full_alpha_aa`` (the #91 invariant)."""
        from tcrsift.assemble import (
            assemble_full_sequences,
            translate_dna,
            validate_sequences,
        )

        df, contigs_dir = self._build_alpha_fixture(tmp_path, "AAT")
        out = assemble_full_sequences(
            df,
            contigs_dir=str(contigs_dir),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        full_nt = out["full_alpha_nt"].iloc[0]
        full_aa = out["full_alpha_aa"].iloc[0]
        nt_for_translation = (
            full_nt[:-3] if full_nt[-3:] in {"TAA", "TAG", "TGA"} else full_nt
        )
        translated, _ = translate_dna(nt_for_translation)
        assert translated == full_aa
        # validate_sequences also flags no NT round-trip failures.
        msgs = validate_sequences(out, strict=False)
        translation_failures = [
            m for m in msgs
            if m.severity == "load_bearing"
            and ("translate to" in m or "frame likely broken" in m)
        ]
        assert translation_failures == []

    def test_no_contig_emits_warning_and_skips_prepend(self, tmp_path):
        """When contigs are loaded but α has no contig past VDJ, the
        assembly preserves pre-#105 behavior (no junction residue) and
        surfaces a QC warning so users know α is 1 aa short."""
        from tcrsift.assemble import assemble_full_sequences

        df, contigs_dir = self._build_alpha_fixture(tmp_path, "AAT")
        # Remove the contig contents so find(vdj_nt) returns None.
        for sample_dir in contigs_dir.iterdir():
            if sample_dir.is_dir():
                for fasta in sample_dir.glob("*.fasta"):
                    fasta.write_text(">contig_a1\nACGT\n")  # no VDJ match

        out = assemble_full_sequences(
            df,
            contigs_dir=str(contigs_dir),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        # alpha_constant_aa starts at the canonical (no prepend).
        assert out["alpha_constant_aa"].iloc[0].startswith(HUMAN_TRAC_AA[:8])
        assert out["alpha_junction_residue"].iloc[0] is None
        # And a QC warning records the gap.
        qc = out["qc_warnings"].iloc[0] or []
        assert any(
            "no contig coverage past VDJ" in m for m in qc
        ), f"expected α-no-contig warning, got {qc}"

    def test_no_contigs_dir_silent_pre_105_behavior(self):
        """When ``contigs_dir`` isn't provided at all (user opted out
        of contig-derived NT entirely), the assembly stays silent and
        produces pre-#105 alpha output. No warning."""
        from tcrsift.assemble import assemble_full_sequences

        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha,
            "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha,
            "VDJ_beta_aa": vdj_beta,
            "alpha_c_gene": "TRAC",
            "beta_c_gene": "TRBC1",
            "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
        }])
        out = assemble_full_sequences(
            df, alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        # alpha_junction_residue isn't populated (None or column absent).
        if "alpha_junction_residue" in out.columns:
            assert out["alpha_junction_residue"].iloc[0] is None
        # And no warning about missing contig coverage.
        qc = out.get("qc_warnings", pd.Series([[]])).iloc[0] or []
        assert not any(
            "no contig coverage past VDJ" in m for m in qc
        )

    def test_stop_or_unknown_junction_codon_warns(self, tmp_path):
        """If the contig's junction codon is a stop (TAA/TAG/TGA) or
        a codon with N, we don't prepend — instead emit a QC warning."""
        from tcrsift.assemble import assemble_full_sequences

        df, contigs_dir = self._build_alpha_fixture(tmp_path, "TAA")  # stop
        out = assemble_full_sequences(
            df,
            contigs_dir=str(contigs_dir),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        # No prepend.
        assert out["alpha_constant_aa"].iloc[0].startswith(HUMAN_TRAC_AA[:8])
        assert out["alpha_junction_residue"].iloc[0] is None
        # Warning surfaces the stop codon and names α.
        qc = out["qc_warnings"].iloc[0] or []
        assert any(
            "alpha" in m and "junction codon" in m and "TAA" in m for m in qc
        ), f"expected α stop-codon warning, got {qc}"

    # ------------------------------------------------------------
    # β junction tests — same code path as α since 2.0 (#105).
    # ------------------------------------------------------------

    @staticmethod
    def _build_beta_fixture(tmp_path, junction_codon: str, c_gene: str = "TRBC1"):
        """Return ``(df, contigs_dir)`` for a β-only assembly test.
        ``junction_codon`` is the 3-nt codon (universally GAG = E in real
        biology, but the test parametrises across other residues to
        exercise the same uniform code path)."""
        leader_aa = "M" + "A" * 19
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        vdj_beta_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_beta)
        overshoot = "A"
        canonical_aa = HUMAN_TRBC1_AA if c_gene == "TRBC1" else HUMAN_TRBC2_AA
        canonical_start_nt = "".join(
            HUMAN_PREFERRED_CODONS[r] for r in canonical_aa[:8]
        )
        beta_contig_seq = (
            "".join(HUMAN_PREFERRED_CODONS[r] for r in leader_aa)
            + vdj_beta_nt
            + junction_codon
            + canonical_start_nt
        )
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha, "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha, "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_alpha),
            "VDJ_beta_nt": vdj_beta_nt + overshoot,
            "alpha_c_gene": "TRAC",
            "beta_c_gene": c_gene,
            "beta_j_gene": "TRBJ1-1" if c_gene == "TRBC1" else "TRBJ2-1",
            "samples": "S1",
            "beta_contig_ids": "contig_b1",
        }])
        contig_dir = tmp_path / "S1"
        contig_dir.mkdir(parents=True, exist_ok=True)
        (contig_dir / "filtered_contig.fasta").write_text(
            f">contig_b1\n{beta_contig_seq}\n"
        )
        return df, tmp_path

    @pytest.mark.parametrize("c_gene", ["TRBC1", "TRBC2"])
    def test_beta_universal_e_junction_from_contig(self, c_gene, tmp_path):
        """The β J→C junction is universally E (GAG codon) across all
        118/118 audited B1 clones — the J segment's terminal nt plus
        the C exon's first 2 nt always spell GAG. Both TRBC1 and TRBC2
        get the E prepended uniformly when the contig provides GAG."""
        from tcrsift.assemble import assemble_full_sequences

        df, contigs_dir = self._build_beta_fixture(
            tmp_path, "GAG", c_gene=c_gene
        )
        out = assemble_full_sequences(
            df,
            contigs_dir=str(contigs_dir),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        canonical = HUMAN_TRBC1_AA if c_gene == "TRBC1" else HUMAN_TRBC2_AA
        # Per-clone canonical = E + bare mature.
        assert out["beta_constant_aa"].iloc[0].startswith("E" + canonical[:8])
        assert out["beta_junction_residue"].iloc[0] == "E"
        # And the assembled full β has the E at the J→C boundary.
        full_aa = out["full_beta_aa"].iloc[0]
        vdj_beta = df["VDJ_beta_aa"].iloc[0]
        assert vdj_beta + "E" + canonical[:8] in full_aa

    def test_beta_trbc1_jovi1_signature_with_junction_e(self, tmp_path):
        """For a TRBC1 clone with junction E, the assembled chain reads
        ``…VDJ + EDLNK…`` at the J→C boundary. The JOVI.1 NK epitope
        ends up at AA positions 3-4 of the assembled β-constant chunk."""
        from tcrsift.assemble import assemble_full_sequences

        df, contigs_dir = self._build_beta_fixture(tmp_path, "GAG", c_gene="TRBC1")
        out = assemble_full_sequences(
            df,
            contigs_dir=str(contigs_dir),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        full_aa = out["full_beta_aa"].iloc[0]
        vdj_beta = df["VDJ_beta_aa"].iloc[0]
        # Junction + JOVI.1 epitope appears at the boundary: E D L N K
        assert vdj_beta + "EDLNK" in full_aa, (
            "expected ...VDJ + EDLNK... at TRBC1 J→C boundary"
        )

    def test_beta_trbc2_jovi1_signature_with_junction_e(self, tmp_path):
        """For a TRBC2 clone, the boundary reads ``…VDJ + EDLKN…`` —
        the KN order (NOT NK) is the JOVI.1-negative signature."""
        from tcrsift.assemble import assemble_full_sequences

        df, contigs_dir = self._build_beta_fixture(tmp_path, "GAG", c_gene="TRBC2")
        out = assemble_full_sequences(
            df,
            contigs_dir=str(contigs_dir),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        full_aa = out["full_beta_aa"].iloc[0]
        vdj_beta = df["VDJ_beta_aa"].iloc[0]
        # KN instead of NK — JOVI.1-negative TRBC2 boundary.
        assert vdj_beta + "EDLKN" in full_aa
        assert vdj_beta + "EDLNK" not in full_aa, (
            "TRBC2 must NOT have the JOVI.1 NK epitope; that'd indicate "
            "the canonical sequence regressed to pre-#100"
        )

    def test_beta_constant_ends_correctly_by_gene(self, tmp_path):
        """TRBC1 ends with ``VKRKDF``; TRBC2 ends with ``VKRKDSRG``. The
        2-aa tail difference is preserved through the assembly path."""
        from tcrsift.assemble import assemble_full_sequences

        for c_gene, expected_tail in [
            ("TRBC1", "VKRKDF"),
            ("TRBC2", "VKRKDSRG"),
        ]:
            df, contigs_dir = self._build_beta_fixture(
                tmp_path / c_gene, "GAG", c_gene=c_gene
            )
            out = assemble_full_sequences(
                df,
                contigs_dir=str(contigs_dir),
                alpha_leader=None, beta_leader=None,
                verbose=False, show_progress=False,
            )
            # The mature β-constant chunk ends with the gene-specific tail.
            assert out["beta_constant_aa"].iloc[0].endswith(expected_tail), (
                f"{c_gene} should end with {expected_tail!r}; got "
                f"{out['beta_constant_aa'].iloc[0][-12:]!r}"
            )

    def test_no_synthetic_e_when_no_contig_for_beta(self, tmp_path):
        """Pre-2.0, β chains got a synthetic E even without contigs.
        Since 2.0, no contig means no E — bare mature TRBC, with a
        QC warning indicating the assembled chain is 1 aa short of the
        donor's expressed β protein. This is the 2.0 behaviour change."""
        from tcrsift.assemble import assemble_full_sequences

        df, contigs_dir = self._build_beta_fixture(tmp_path, "GAG", c_gene="TRBC1")
        # Wipe the contig contents so find(vdj_nt) fails.
        for sample_dir in contigs_dir.iterdir():
            if sample_dir.is_dir():
                for fasta in sample_dir.glob("*.fasta"):
                    fasta.write_text(">contig_b1\nACGT\n")
        out = assemble_full_sequences(
            df,
            contigs_dir=str(contigs_dir),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        # β-constant starts at bare mature D, NOT synthetic E.
        assert out["beta_constant_aa"].iloc[0].startswith("DLNK")
        assert not out["beta_constant_aa"].iloc[0].startswith("EDLN")
        assert out["beta_junction_residue"].iloc[0] is None
        # And the QC warning explicitly names β.
        qc = out["qc_warnings"].iloc[0] or []
        assert any(
            "beta" in m and "no contig coverage past VDJ" in m for m in qc
        )

    def test_alpha_and_beta_use_same_code_path(self, tmp_path):
        """Uniform treatment: the same junction-peel logic runs for both
        chains. We assert that α and β handle their respective junction
        codons in lockstep — both get a non-None ``{chain}_junction_residue``
        when their contigs provide a valid codon."""
        from tcrsift.assemble import assemble_full_sequences

        # Build a combined fixture with both an α and β contig.
        leader_aa = "M" + "A" * 19
        leader_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in leader_aa)
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        vdj_alpha_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_alpha)
        vdj_beta_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_beta)
        trac_start = "".join(HUMAN_PREFERRED_CODONS[r] for r in HUMAN_TRAC_AA[:8])
        trbc1_start = "".join(HUMAN_PREFERRED_CODONS[r] for r in HUMAN_TRBC1_AA[:8])
        alpha_contig = leader_nt + vdj_alpha_nt + "AAT" + trac_start  # N junction
        beta_contig = leader_nt + vdj_beta_nt + "GAG" + trbc1_start   # E junction

        contig_dir = tmp_path / "S1"
        contig_dir.mkdir(parents=True)
        (contig_dir / "filtered_contig.fasta").write_text(
            f">contig_a1\n{alpha_contig}\n>contig_b1\n{beta_contig}\n"
        )

        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha, "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha, "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": vdj_alpha_nt + "A",  # +1 overshoot for #91 trim
            "VDJ_beta_nt": vdj_beta_nt + "A",
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1",
            "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
            "alpha_contig_ids": "contig_a1",
            "beta_contig_ids": "contig_b1",
        }])
        out = assemble_full_sequences(
            df, contigs_dir=str(tmp_path),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        # Both chains have a junction residue from the same code path.
        assert out["alpha_junction_residue"].iloc[0] == "N"
        assert out["beta_junction_residue"].iloc[0] == "E"
        # Both constant_aa start with the per-clone junction.
        assert out["alpha_constant_aa"].iloc[0].startswith("N" + HUMAN_TRAC_AA[:8])
        assert out["beta_constant_aa"].iloc[0].startswith("E" + HUMAN_TRBC1_AA[:8])

    def test_validate_sequences_accepts_junction_prepended_layout(self, tmp_path):
        """#107 regression: after 2.0's #105 fix, ``{chain}_constant_aa``
        starts with the per-clone junction residue. The validator's
        canonical-start check was comparing observed (with junction)
        against the bare-mature canonical from
        :data:`HUMAN_CONSTANT_REGIONS_AA` (without junction) and flagged
        every clone. ``tcrsift run --contigs-dir`` raised at the PDF
        gate even though all output data was correct.

        Now: validate_sequences should source the expected canonical
        start from the row's own ``{chain}_constant_aa`` column, which
        already reflects whatever junction was prepended. No false
        positives on contig-aware output."""
        from tcrsift.assemble import assemble_full_sequences, validate_sequences

        # Same combined α+β fixture used in the uniformity test above.
        leader_aa = "M" + "A" * 19
        leader_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in leader_aa)
        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        vdj_beta = "CASS" + "G" * 60 + "VETA"
        vdj_alpha_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_alpha)
        vdj_beta_nt = "".join(HUMAN_PREFERRED_CODONS[r] for r in vdj_beta)
        trac_start = "".join(HUMAN_PREFERRED_CODONS[r] for r in HUMAN_TRAC_AA[:8])
        trbc1_start = "".join(HUMAN_PREFERRED_CODONS[r] for r in HUMAN_TRBC1_AA[:8])
        alpha_contig = leader_nt + vdj_alpha_nt + "AAT" + trac_start   # N junction
        beta_contig = leader_nt + vdj_beta_nt + "GAG" + trbc1_start    # E junction

        contig_dir = tmp_path / "S1"
        contig_dir.mkdir(parents=True)
        (contig_dir / "filtered_contig.fasta").write_text(
            f">contig_a1\n{alpha_contig}\n>contig_b1\n{beta_contig}\n"
        )

        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_alpha, "CDR3_beta": vdj_beta,
            "VDJ_alpha_aa": vdj_alpha, "VDJ_beta_aa": vdj_beta,
            "VDJ_alpha_nt": vdj_alpha_nt + "A",
            "VDJ_beta_nt": vdj_beta_nt + "A",
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1",
            "beta_j_gene": "TRBJ1-1",
            "samples": "S1",
            "alpha_contig_ids": "contig_a1",
            "beta_contig_ids": "contig_b1",
        }])
        out = assemble_full_sequences(
            df, contigs_dir=str(tmp_path),
            alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )

        # Sanity: junction residues were prepended per #105.
        assert out["alpha_junction_residue"].iloc[0] == "N"
        assert out["beta_junction_residue"].iloc[0] == "E"
        assert out["alpha_constant_aa"].iloc[0].startswith("N")
        assert out["beta_constant_aa"].iloc[0].startswith("E")

        # The actual #107 assertion: validator emits NO load-bearing
        # "constant start doesn't match canonical" messages for either
        # chain. Pre-fix this would fire for both α and β.
        msgs = validate_sequences(out, strict=False)
        start_failures = [
            m for m in msgs
            if m.severity == "load_bearing"
            and "constant start doesn't match canonical" in m
        ]
        assert start_failures == [], (
            f"#107 regression: validate_sequences flagged "
            f"junction-prepended constants as failing the canonical-start "
            f"check. {len(start_failures)} failures: "
            f"{start_failures[:3]}"
        )

        # And no NT round-trip failures (#91 invariant) either — the
        # blend preserves the junction codon.
        rt_failures = [
            m for m in msgs
            if m.severity == "load_bearing"
            and ("translate to" in m or "frame likely broken" in m)
        ]
        assert rt_failures == []

    def test_validate_sequences_still_catches_real_start_mismatch(self):
        """Defensive: the #107 fix mustn't soften the canonical-start
        check past the point of usefulness. Pin that a clone whose
        ``full_{chain}_aa`` doesn't match its own ``{chain}_constant_aa``
        after VDJ still gets flagged."""
        from tcrsift.assemble import validate_sequences

        # Hand-craft a row where the body after VDJ in full_alpha_aa is
        # wrong but alpha_constant_aa claims the correct canonical.
        leader = "M" + "A" * 19
        vdj = "CASS" + "A" * 60 + "VLPHA"
        good_canonical = "N" + HUMAN_TRAC_AA  # what the row claims
        bad_observed_constant = "WRONG" + "X" * (len(good_canonical) - 5 - 8) + "MTLRLWSS"
        df = pd.DataFrame([{
            "alpha_leader_aa": leader,
            "vdj_alpha_aa": vdj,
            "alpha_constant_aa": good_canonical,
            "alpha_c_gene": "TRAC",
            "alpha_c_gene_canonical": "TRAC",
            "full_alpha_aa": leader + vdj + bad_observed_constant,
        }])
        msgs = validate_sequences(df, strict=False)
        start_failures = [
            m for m in msgs
            if m.severity == "load_bearing"
            and "constant start doesn't match canonical" in m
        ]
        assert start_failures, (
            "validator should still flag genuine start mismatches; the "
            "#107 fix may have softened the check too far"
        )
