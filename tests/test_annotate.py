"""
Tests for TCR annotation using public databases.
"""

import numpy as np
import pandas as pd
import pytest

from tcrsift.annotate import (
    _flag_viral,
    annotate_clonotypes,
    classify_category,
    get_annotation_summary,
    load_cedar,
    load_databases,
    load_iedb,
    load_vdjdb,
    match_clonotypes,
)
from tcrsift.validation import TCRsiftValidationError


class TestFlagViral:
    """Tests for viral flagging."""

    def test_flag_cmv_viral(self):
        """CMV should be flagged as viral."""
        df = pd.DataFrame({"species": ["CMV", "Cytomegalovirus", "Human CMV"]})
        result = _flag_viral(df)
        assert all(result)

    def test_flag_ebv_viral(self):
        """EBV should be flagged as viral."""
        df = pd.DataFrame({"species": ["EBV", "Epstein-Barr virus"]})
        result = _flag_viral(df)
        assert all(result)

    def test_flag_hiv_viral(self):
        """HIV should be flagged as viral."""
        df = pd.DataFrame({"species": ["HIV", "Human Immunodeficiency Virus"]})
        result = _flag_viral(df)
        assert all(result)

    def test_non_viral_not_flagged(self):
        """Non-viral species should not be flagged."""
        df = pd.DataFrame({"species": ["Homo sapiens", "Mouse", "Self-antigen"]})
        result = _flag_viral(df)
        assert not any(result)

    def test_missing_species_not_flagged(self):
        """Missing species column should not be flagged."""
        df = pd.DataFrame({"epitope": ["PEPTIDE1", "PEPTIDE2"]})
        result = _flag_viral(df)
        assert not any(result)

    def test_null_species_not_flagged(self):
        """Null species should not be flagged."""
        df = pd.DataFrame({"species": [None, np.nan, ""]})
        result = _flag_viral(df)
        assert not any(result)


class TestClassifyCategory:
    """Tests for the species/antigen → category classifier."""

    def test_viral(self):
        cats = classify_category(
            pd.Series(["Human cytomegalovirus", "Epstein-Barr virus"]),
            pd.Series(["pp65", "EBNA-3"]),
        )
        assert (cats == "viral").all()

    def test_bacterial(self):
        cats = classify_category(
            pd.Series(["Mycobacterium tuberculosis", "Listeria monocytogenes"]),
            pd.Series(["ESAT-6", "LLO"]),
        )
        assert (cats == "bacterial").all()

    def test_self_homo_sapiens(self):
        cats = classify_category(
            pd.Series(["Homo sapiens (human)", "Homo sapiens"]),
            pd.Series(["beta-2-microglobulin", "insulin"]),
        )
        assert (cats == "self").all()

    def test_tumor_self_overrides_species(self):
        """MART-1 / NY-ESO-1 / MAGE peptides are Homo sapiens but should
        not bucket as plain self — they're tumor-associated self antigens."""
        cats = classify_category(
            pd.Series(["Homo sapiens", "Homo sapiens", "Homo sapiens"]),
            pd.Series(["MLANA (MART-1)", "NY-ESO-1", "MAGE-A3"]),
        )
        assert (cats == "tumor_self").all()

    def test_unknown_species(self):
        cats = classify_category(
            pd.Series(["", None]),
            pd.Series(["", None]),
        )
        assert (cats == "unknown").all()

    def test_other_known_but_unclassified(self):
        cats = classify_category(
            pd.Series(["Plasmodium falciparum"]),
            pd.Series(["CSP"]),
        )
        assert (cats == "other").all()

    def test_tumor_self_patterns_use_word_boundaries(self):
        """Short tumor-self tokens (her2, tert, psa, wt1) are matched at
        word boundaries — substrings inside unrelated names must not
        trigger the tumor_self override."""
        # Antigens that *contain* the short tokens as substrings but are
        # not actually tumor-self antigens. They should fall through to
        # "self" (Homo sapiens) rather than getting overridden.
        cats = classify_category(
            pd.Series(["Homo sapiens"] * 4),
            pd.Series([
                "interterritorial protein",  # contains "tert"
                "altered self",              # contains "tert" inside "altered"
                "lipase A",                  # nothing tumor-related
                "phosphatase",               # contains "psa" as substring
            ]),
        )
        assert (cats == "self").all()

    def test_tumor_self_patterns_match_standalone(self):
        """The same short tokens *do* fire when they appear as standalone
        antigen names (the actual VDJdb/IEDB convention)."""
        cats = classify_category(
            pd.Series(["Homo sapiens"] * 4),
            pd.Series(["HER2", "TERT", "PSA", "WT1"]),
        )
        assert (cats == "tumor_self").all()


class TestLoadVdjdb:
    """Tests for loading VDJdb.

    Test data is sampled from the real ``vdjdb_full.txt`` and
    ``vdjdb.txt`` files in the release zip so we exercise the actual
    on-disk shapes and sequence conventions (β-CDR3s start ``CAS``,
    α-CDR3s start ``CAV``/``CAA``, etc.) rather than synthetic stand-ins.
    """

    @pytest.fixture
    def mock_vdjdb_file(self, temp_dir):
        """Create a mock VDJdb file in the paired-chain format.

        Uses two real VDJdb rows (InfluenzaB NS1 / HLA-B*07:02
        epitope HPNGYKSLSTL) so the schema and CDR3 motifs are
        authentic.
        """
        vdjdb_path = temp_dir / "vdjdb_full.txt"
        df = pd.DataFrame(
            {
                "cdr3.alpha": ["CAYRSAGSGGSNYKLTF", "CAYRSATGYALNF"],
                "cdr3.beta": ["CASSLMTNQPQHF", "CASSLYGSVQNEQFF"],
                "antigen.epitope": ["HPNGYKSLSTL", "HPNGYKSLSTL"],
                "antigen.gene": ["NS1", "NS1"],
                "antigen.species": ["InfluenzaB", "InfluenzaB"],
                "mhc.a": ["HLA-B*07:02", "HLA-B*07:02"],
                "species": ["HomoSapiens", "HomoSapiens"],  # donor — to be dropped
            }
        )
        df.to_csv(vdjdb_path, sep="\t", index=False)
        return vdjdb_path

    @pytest.fixture
    def mock_vdjdb_dir(self, temp_dir, mock_vdjdb_file):
        """Create a mock VDJdb directory."""
        return temp_dir

    def test_load_from_file(self, mock_vdjdb_file):
        """Load VDJdb from file."""
        result = load_vdjdb(mock_vdjdb_file)

        assert len(result) == 2
        assert "cdr3_beta" in result.columns
        assert "cdr3_alpha" in result.columns
        assert "epitope" in result.columns
        assert "is_viral" in result.columns
        assert result["database"].iloc[0] == "VDJdb"

    def test_load_from_directory(self, mock_vdjdb_dir):
        """Load VDJdb from directory."""
        result = load_vdjdb(mock_vdjdb_dir)
        assert len(result) == 2

    def test_viral_flagging(self, mock_vdjdb_file):
        """VDJdb entries should have viral flag."""
        result = load_vdjdb(mock_vdjdb_file)
        # Both fixture rows are Influenza B (viral).
        assert all(result["is_viral"])

    def test_load_picks_vdjdb_full_when_present(self, temp_dir):
        """Directory-mode load picks the canonical ``vdjdb_full.txt``
        and ignores any unrelated files (sidecar metadata,
        partial/processed views) in the same directory (#45 bug 1).

        Uses a real VDJdb row (InfluenzaB NS1 / HLA-B*07:02).
        """
        full = pd.DataFrame(
            {
                "cdr3.alpha": ["CAYRSAGSGGSNYKLTF"],
                "cdr3.beta": ["CASSLMTNQPQHF"],
                "antigen.epitope": ["HPNGYKSLSTL"],
                "antigen.gene": ["NS1"],
                "antigen.species": ["InfluenzaB"],
                "species": ["HomoSapiens"],  # donor — should be dropped
                "mhc.a": ["HLA-B*07:02"],
            }
        )
        full.to_csv(temp_dir / "vdjdb_full.txt", sep="\t", index=False)
        # A metadata sidecar that alphabetizes before vdjdb_full.txt.
        pd.DataFrame(
            {"name": ["complex.id"], "type": ["txt"], "comment": ["..."]}
        ).to_csv(temp_dir / "vdjdb.meta.txt", sep="\t", index=False)

        result = load_vdjdb(temp_dir)
        assert len(result) == 1
        assert result.iloc[0]["cdr3_alpha"] == "CAYRSAGSGGSNYKLTF"
        assert result.iloc[0]["cdr3_beta"] == "CASSLMTNQPQHF"

    def test_load_falls_back_to_vdjdb_txt_when_full_missing(self, temp_dir):
        """Priority-2 path: when ``vdjdb_full.txt`` isn't present but
        ``vdjdb.txt`` is, the long-format file must be picked from a
        directory load (and then gene-filtered by the long normalizer).
        Uses a real TRB row from VDJdb (InfluenzaA M1 / GILGFVFTL).
        """
        pd.DataFrame(
            {
                "complex.id": [0],
                "gene": ["TRB"],
                "cdr3": ["CASSIVGGNEQFF"],
                "antigen.epitope": ["GILGFVFTL"],
                "antigen.species": ["InfluenzaA"],
            }
        ).to_csv(temp_dir / "vdjdb.txt", sep="\t", index=False)
        # Sidecars that mustn't be picked.
        pd.DataFrame({"name": ["foo"]}).to_csv(
            temp_dir / "vdjdb.meta.txt", sep="\t", index=False
        )

        result = load_vdjdb(temp_dir)
        assert len(result) == 1
        assert result.iloc[0]["cdr3_beta"] == "CASSIVGGNEQFF"

    def test_load_falls_back_to_vdjdb_slim_when_others_missing(self, temp_dir):
        """Priority-3 path: when only ``vdjdb.slim.txt`` is present,
        the slim file should be picked. Uses a real TRB row.
        """
        pd.DataFrame(
            {
                "complex.id": [0],
                "gene": ["TRB"],
                "cdr3": ["CASSMRSTGELFF"],
                "antigen.epitope": ["GILGFVFTL"],
                "antigen.species": ["InfluenzaA"],
            }
        ).to_csv(temp_dir / "vdjdb.slim.txt", sep="\t", index=False)
        # Sidecars present; must not be picked over the canonical slim.
        pd.DataFrame({"name": ["foo"]}).to_csv(
            temp_dir / "vdjdb.meta.txt", sep="\t", index=False
        )

        result = load_vdjdb(temp_dir)
        assert len(result) == 1
        assert result.iloc[0]["cdr3_beta"] == "CASSMRSTGELFF"

    def test_load_maps_cdr3_dot_beta_to_cdr3_beta(self, temp_dir):
        """``vdjdb_full.txt`` uses ``cdr3.beta`` (paired-chain row format),
        not ``cdr3`` (slim/long format). Without ``cdr3.beta`` in the
        mapping, ``match_clonotypes`` crashes with KeyError when the
        user explicitly points at the full file (#45 bug 1)."""
        path = temp_dir / "vdjdb_full.txt"
        pd.DataFrame(
            {
                "cdr3.alpha": ["CAYRSAGSGGSNYKLTF"],
                "cdr3.beta": ["CASSLMTNQPQHF"],
                "antigen.epitope": ["HPNGYKSLSTL"],
                "antigen.species": ["InfluenzaB"],
            }
        ).to_csv(path, sep="\t", index=False)

        result = load_vdjdb(path)
        assert "cdr3_beta" in result.columns
        assert result["cdr3_beta"].iloc[0] == "CASSLMTNQPQHF"

    def test_donor_species_does_not_collide_with_antigen_species(self, temp_dir):
        """VDJdb's full export carries both a donor ``species`` column
        (T-cell origin, usually ``HomoSapiens``) and ``antigen.species``
        (the epitope source organism, e.g. ``InfluenzaB``). The loader
        uses the latter for matching/category, so the former must not
        leak through as a same-name column collision (#45 bug 1).
        """
        path = temp_dir / "vdjdb_full.txt"
        pd.DataFrame(
            {
                "cdr3.alpha": ["CAYRSAGSGGSNYKLTF"],
                "cdr3.beta": ["CASSLMTNQPQHF"],
                "antigen.epitope": ["HPNGYKSLSTL"],
                "antigen.species": ["InfluenzaB"],
                "species": ["HomoSapiens"],
            }
        ).to_csv(path, sep="\t", index=False)

        result = load_vdjdb(path)
        # Only one column named ``species``, and it carries the antigen
        # source organism — not the donor species.
        assert (result.columns == "species").sum() == 1
        assert result["species"].iloc[0] == "InfluenzaB"

    def test_load_raises_when_no_canonical_file_present(self, temp_dir):
        """A directory without any of the canonical names (only
        sidecar/metadata files) must raise with a hint telling the
        user how to recover (#45 bug 1)."""
        pd.DataFrame({"name": ["foo"]}).to_csv(
            temp_dir / "vdjdb.meta.txt", sep="\t", index=False
        )
        pd.DataFrame({"name": ["bar"]}).to_csv(
            temp_dir / "vdjdb.slim.meta.txt", sep="\t", index=False
        )

        with pytest.raises(
            TCRsiftValidationError, match="No canonical VDJdb data file"
        ):
            load_vdjdb(temp_dir)

    def test_long_format_filters_to_trb_to_avoid_alpha_leakage(self, temp_dir):
        """The slim/long VDJdb file uses one row per chain — ``cdr3``
        holds an α-CDR3 when ``gene=TRA`` and a β-CDR3 when
        ``gene=TRB``. A naive ``cdr3 → cdr3_beta`` mapping puts α
        sequences into ``cdr3_beta`` and causes wrong-chain matches
        downstream. The loader must filter to TRB rows so
        ``cdr3_beta`` only carries actual β-CDR3s.

        Uses real VDJdb sequences: ``CAASGGYQKVTF`` is a real α-CDR3
        (HIV-1, gene=TRA); ``CASSIVGGNEQFF`` is a real β-CDR3 (InfluenzaA
        M1 / GILGFVFTL).
        """
        path = temp_dir / "vdjdb.txt"
        pd.DataFrame(
            {
                "complex.id": [0, 0],
                "gene": ["TRA", "TRB"],
                "cdr3": [
                    "CAASGGYQKVTF",   # α — must NOT land in cdr3_beta
                    "CASSIVGGNEQFF",  # β — should land in cdr3_beta
                ],
                "antigen.epitope": ["KAFSPEVIPMF", "GILGFVFTL"],
                "antigen.species": ["HIV-1", "InfluenzaA"],
            }
        ).to_csv(path, sep="\t", index=False)

        result = load_vdjdb(path)

        # Only the TRB row survives.
        assert len(result) == 1
        assert result.iloc[0]["cdr3_beta"] == "CASSIVGGNEQFF"
        # Crucially: the α-CDR3 sequence did NOT leak into cdr3_beta.
        assert "CAASGGYQKVTF" not in set(result["cdr3_beta"].dropna())
        # cdr3_alpha column present but empty — schema stability.
        assert "cdr3_alpha" in result.columns
        assert result["cdr3_alpha"].isna().all()

    def test_unrecognized_schema_raises_clear_error(self, temp_dir):
        """A file with neither paired nor long-format VDJdb columns
        should raise a validation error naming the columns it saw,
        rather than silently producing an empty/garbage result."""
        path = temp_dir / "vdjdb_weird.txt"
        pd.DataFrame(
            {"random_col": ["foo"], "another_col": ["bar"]}
        ).to_csv(path, sep="\t", index=False)

        with pytest.raises(TCRsiftValidationError, match="Unrecognized VDJdb schema"):
            load_vdjdb(path)


class TestLoadIedb:
    """Tests for loading IEDB."""

    @pytest.fixture
    def mock_iedb_file(self, temp_dir):
        """Create a mock IEDB file."""
        iedb_path = temp_dir / "iedb.tsv"
        df = pd.DataFrame(
            {
                "Chain 2 CDR3 Curated": ["CASSLGQAYEQYF"],
                "Chain 1 CDR3 Curated": ["CAVSDGGSQGNLIF"],
                "Epitope - Name": ["pp65"],
                "Epitope - Source Organism Name": ["CMV"],
            }
        )
        df.to_csv(iedb_path, sep="\t", index=False)
        return iedb_path

    def test_load_iedb(self, mock_iedb_file):
        """Load IEDB file."""
        result = load_iedb(mock_iedb_file)

        assert len(result) == 1
        assert "cdr3_beta" in result.columns
        assert result["database"].iloc[0] == "IEDB"

    def test_load_iedb_v3_hierarchical_header(self, temp_dir):
        """Parse the current v3 export shape: comma-separated, two-row header.

        The header starts ``Receptor,Receptor,Receptor,...`` and fields
        live one level down (``CDR3 Curated`` under ``Chain 1`` etc.).
        Synthesised here as a small CSV so the test doesn't depend on
        a 99MB cache file.
        """
        path = temp_dir / "tcr_full_v3.csv"
        content = (
            # Row 1 — section names (repeated per field)
            "Receptor,Receptor,Epitope,Epitope,Epitope,Assay,"
            "Chain 1,Chain 1,Chain 2,Chain 2\n"
            # Row 2 — field names
            "IEDB Receptor ID,Type,Name,Source Molecule,Source Organism,"
            "MHC Allele Names,Type,CDR3 Curated,Type,CDR3 Curated\n"
            # Data — one αβ row that should be picked up
            "1,alphabeta,NLVPMVATV,pp65,Human cytomegalovirus,HLA-A*02:01,"
            "alpha,CAVSDGGSQGNLIF,beta,CASSLGQAYEQYF\n"
            # A gammadelta row that should be dropped (not αβ-matchable)
            "2,gammadelta,SOMEPEP,foo,Homo sapiens,HLA-A*02:01,"
            "gamma,CAVGAMMA,delta,CASSDELTA\n"
        )
        path.write_text(content)

        result = load_iedb(path)

        assert len(result) == 1
        row = result.iloc[0]
        assert row["cdr3_alpha"] == "CAVSDGGSQGNLIF"
        assert row["cdr3_beta"] == "CASSLGQAYEQYF"
        assert row["epitope"] == "NLVPMVATV"
        assert row["antigen_gene"] == "pp65"
        assert row["species"] == "Human cytomegalovirus"
        assert row["mhc_allele"] == "HLA-A*02:01"
        assert bool(row["is_viral"]) is True
        assert row["database"] == "IEDB"

    def test_load_iedb_v3_falls_back_to_calculated_cdr3(self, temp_dir):
        """If a chain has no ``CDR3 Curated`` value but ``CDR3 Calculated``
        is present, the calculated value is used. Mirrors the priority
        VDJdb applies to its curated-over-inferred fields."""
        path = temp_dir / "tcr_full_v3.csv"
        content = (
            "Receptor,Receptor,Epitope,Chain 1,Chain 1,Chain 1,"
            "Chain 2,Chain 2,Chain 2\n"
            "IEDB Receptor ID,Type,Name,Type,CDR3 Curated,CDR3 Calculated,"
            "Type,CDR3 Curated,CDR3 Calculated\n"
            "1,alphabeta,EPI1,alpha,,CAVCALCALPHA,beta,CASSCURBETA,CASSCALCBETA\n"
        )
        path.write_text(content)

        result = load_iedb(path)

        assert len(result) == 1
        assert result.iloc[0]["cdr3_alpha"] == "CAVCALCALPHA"
        # Curated wins over Calculated when both are present.
        assert result.iloc[0]["cdr3_beta"] == "CASSCURBETA"


class TestIedbEpitopeOverride:
    """Tests for the optional epitope-table override (#54).

    Synthetic mini-files mirror the real on-disk shapes:
    - Receptor file: hierarchical CSV with ``Receptor``/``Epitope``/
      ``Assay``/``Chain 1``/``Chain 2`` sections.
    - Epitope file: hierarchical CSV with ``Epitope ID`` (single-level)
      + ``Epitope`` / ``Related Object`` sections.

    The receptor uses ``https://www.iedb.org/epitope/N`` IRIs; the
    epitope file uses ``http://www.iedb.org/epitope/N``. The override
    must normalize the URL scheme before joining.
    """

    @pytest.fixture
    def receptor_csv(self, temp_dir):
        path = temp_dir / "tcr_full_v3.csv"
        content = (
            # Row 1 — section names
            "Receptor,Receptor,Epitope,Epitope,Epitope,Epitope,Assay,"
            "Chain 1,Chain 1,Chain 2,Chain 2\n"
            # Row 2 — field names
            "IEDB Receptor ID,Type,IEDB IRI,Name,Source Molecule,Source Organism,"
            "MHC Allele Names,Type,CDR3 Curated,Type,CDR3 Curated\n"
            # Row A: long receptor name to be replaced by shorter epitope name
            "1,alphabeta,https://www.iedb.org/epitope/42,LLFGYPVYV,"
            "transcriptional activator Tax,Human T-cell leukemia virus type I,"
            "HLA-A*02:01,alpha,CAVALPHA1,beta,CASSBETA1\n"
            # Row B: receptor SM is empty; epitope table fills it in
            "2,alphabeta,https://www.iedb.org/epitope/43,EEYLKAWTF,,,"
            "HLA-A*02:01,alpha,CAVALPHA2,beta,CASSBETA2\n"
            # Row C: epitope table has no record for this IRI; receptor value retained
            "3,alphabeta,https://www.iedb.org/epitope/999,NLVPMVATV,"
            "65 kDa phosphoprotein,Human cytomegalovirus,"
            "HLA-A*02:01,alpha,CAVALPHA3,beta,CASSBETA3\n"
        )
        path.write_text(content)
        return path

    @pytest.fixture
    def epitope_csv(self, temp_dir):
        path = temp_dir / "epitope_full_v3.csv"
        # Real epitope file has many columns; we synthesize the minimum
        # ``load_iedb_epitope_lookup`` needs (top-level "Epitope ID" with
        # "IEDB IRI" plus "Epitope" with "Source Molecule" / "Source Organism").
        content = (
            # Row 1 — section names. Notice the http:// (vs receptor's https://)
            "Epitope ID,Epitope,Epitope,Epitope\n"
            # Row 2 — field names
            "IEDB IRI,Name,Source Molecule,Source Organism\n"
            # Match for receptor row A: shorter canonical name.
            "http://www.iedb.org/epitope/42,LLFGYPVYV,Protein Tax-1,"
            "Human T-cell leukemia virus type I\n"
            # Match for receptor row B: fills the blank.
            "http://www.iedb.org/epitope/43,EEYLKAWTF,EBNA-3,Epstein-Barr virus\n"
        )
        path.write_text(content)
        return path

    def test_override_shortens_canonical_names(self, receptor_csv, epitope_csv):
        """Receptor row whose Source Molecule is long ("transcriptional
        activator Tax") gets overridden with the epitope-table's
        shorter form ("Protein Tax-1") — the core win for #54."""
        result = load_iedb(receptor_csv, epitope_path=epitope_csv)
        rows = {r["cdr3_beta"]: r for _, r in result.iterrows()}
        assert rows["CASSBETA1"]["antigen_gene"] == "Protein Tax-1"

    def test_override_fills_blank_when_epitope_has_data(self, receptor_csv, epitope_csv):
        """When the receptor row's Source Molecule is empty but the
        epitope table has a value for that IRI, the blank is filled."""
        result = load_iedb(receptor_csv, epitope_path=epitope_csv)
        rows = {r["cdr3_beta"]: r for _, r in result.iterrows()}
        assert rows["CASSBETA2"]["antigen_gene"] == "EBNA-3"
        assert rows["CASSBETA2"]["species"] == "Epstein-Barr virus"

    def test_override_keeps_receptor_value_when_epitope_missing(
        self, receptor_csv, epitope_csv
    ):
        """An IRI absent from the epitope table keeps its receptor value
        — the override is non-destructive."""
        result = load_iedb(receptor_csv, epitope_path=epitope_csv)
        rows = {r["cdr3_beta"]: r for _, r in result.iterrows()}
        assert rows["CASSBETA3"]["antigen_gene"] == "65 kDa phosphoprotein"
        assert rows["CASSBETA3"]["species"] == "Human cytomegalovirus"

    def test_without_epitope_path_returns_receptor_values_unchanged(
        self, receptor_csv
    ):
        """When ``epitope_path`` isn't passed, the receptor values pass
        through unchanged — load_iedb's default behavior is unaffected."""
        result = load_iedb(receptor_csv)  # no epitope_path
        rows = {r["cdr3_beta"]: r for _, r in result.iterrows()}
        assert rows["CASSBETA1"]["antigen_gene"] == "transcriptional activator Tax"

    def test_https_to_http_iri_normalization(self, receptor_csv, epitope_csv):
        """Receptor IRIs are ``https://``, epitope IRIs are ``http://``.
        The join must normalize before matching — confirms the
        ``_normalize_iedb_iri`` step is wired in correctly."""
        result = load_iedb(receptor_csv, epitope_path=epitope_csv)
        # If normalization weren't applied, no row would override.
        assert any(result["antigen_gene"] == "Protein Tax-1")


class TestLoadCedar:
    """Tests for loading CEDAR."""

    @pytest.fixture
    def mock_cedar_file(self, temp_dir):
        """Create a mock CEDAR file."""
        cedar_path = temp_dir / "cedar.tsv"
        df = pd.DataFrame(
            {
                "cdr3_b_aa": ["CASSLGQAYEQYF"],
                "cdr3_a_aa": ["CAVSDGGSQGNLIF"],
                "epitope_sequence": ["NLVPMVATV"],
                "organism": ["CMV"],
            }
        )
        df.to_csv(cedar_path, sep="\t", index=False)
        return cedar_path

    def test_load_cedar(self, mock_cedar_file):
        """Load CEDAR file."""
        result = load_cedar(mock_cedar_file)

        assert len(result) == 1
        assert "cdr3_beta" in result.columns
        assert result["database"].iloc[0] == "CEDAR"


class TestLoadDatabases:
    """Tests for combined database loading."""

    def test_load_combined(self, temp_dir):
        """Load multiple databases."""
        # Create mock files
        vdjdb_path = temp_dir / "vdjdb.tsv"
        iedb_path = temp_dir / "iedb.tsv"

        pd.DataFrame(
            {
                "cdr3.alpha": ["CAVTEST1"],
                "cdr3.beta": ["CASSTEST1"],
                "antigen.species": ["CMV"],
            }
        ).to_csv(vdjdb_path, sep="\t", index=False)

        pd.DataFrame(
            {
                "Chain 2 CDR3 Curated": ["CASSTEST2"],
                "Chain 1 CDR3 Curated": ["CAVTEST2"],
                "Epitope - Source Organism Name": ["EBV"],
            }
        ).to_csv(iedb_path, sep="\t", index=False)

        result = load_databases(vdjdb_path=vdjdb_path, iedb_path=iedb_path)

        assert len(result) == 2
        assert set(result["database"].unique()) == {"VDJdb", "IEDB"}

    def test_load_raises_without_any_db(self):
        """Should raise if no database provided."""
        with pytest.raises(ValueError, match="At least one database"):
            load_databases()


class TestMatchClonotypes:
    """Tests for clonotype matching."""

    def test_match_by_cdr3ab(self, sample_clonotypes_df, sample_database_df):
        """Match by both alpha and beta chains."""
        result = match_clonotypes(
            sample_clonotypes_df,
            sample_database_df,
            match_by="CDR3ab",
        )

        assert "db_match" in result.columns
        assert "db_epitope" in result.columns
        assert "is_viral" in result.columns

        # First clone should match (same CDR3ab)
        assert result["db_match"].iloc[0]

    def test_match_by_cdr3b_only(self, sample_clonotypes_df, sample_database_df):
        """Match by beta chain only."""
        result = match_clonotypes(
            sample_clonotypes_df,
            sample_database_df,
            match_by="CDR3b_only",
        )

        # More matches expected with beta-only matching
        assert result["db_match"].sum() >= 1

    def test_no_match_annotation(self, sample_clonotypes_df):
        """Clonotypes without matches should have empty annotations."""
        # Database with non-matching sequences
        non_matching_db = pd.DataFrame(
            {
                "cdr3_beta": ["CASSNOMATCH"],
                "cdr3_alpha": ["CAVNOMATCH"],
                "epitope": ["NOTFOUND"],
                "species": ["Unknown"],
                "database": ["TestDB"],
                "is_viral": [False],
            }
        )

        result = match_clonotypes(sample_clonotypes_df, non_matching_db, match_by="CDR3ab")

        assert result["db_match"].sum() == 0
        assert (~result["is_viral"]).all()

    def test_skips_with_warning_when_database_missing_cdr3_columns(
        self, sample_clonotypes_df, caplog
    ):
        """When a database loader silently produces a DataFrame without
        `cdr3_alpha`/`cdr3_beta` (e.g. format drift in IEDB's v3 CSV
        export — the column-mapping in load_iedb expects flat headers
        but v3 ships hierarchical ones), match_clonotypes should log a
        clear warning and return the clonotypes unannotated rather than
        crashing with KeyError (#46)."""
        import logging

        # Simulate a load_iedb output where the column-rename map missed
        # every header — the result is a DataFrame with raw IEDB headers
        # and no `cdr3_alpha`/`cdr3_beta`.
        malformed_db = pd.DataFrame(
            {
                "Chain 1": ["foo"],
                "Chain 2": ["bar"],
                "database": ["IEDB"],
                "is_viral": [False],
            }
        )

        with caplog.at_level(logging.WARNING, logger="tcrsift.annotate"):
            result = match_clonotypes(
                sample_clonotypes_df, malformed_db, match_by="CDR3ab"
            )

        # Annotation columns are present and empty (initialized then
        # returned without matching).
        assert "db_match" in result.columns
        assert not result["db_match"].any()
        # Warning names the missing columns and points at the cause.
        warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert any("missing required columns" in m for m in warnings)
        assert any("cdr3_alpha" in m for m in warnings)

    def test_structured_columns_present_and_populated(self, sample_clonotypes_df):
        """match_clonotypes emits the structured-annotation columns added
        in #48 and populates them with the matched values, not just placeholders."""
        sample = sample_clonotypes_df.iloc[0]
        db = pd.DataFrame(
            {
                "cdr3_beta": [sample["CDR3_beta"]],
                "cdr3_alpha": [sample["CDR3_alpha"]],
                "epitope": ["NLVPMVATV"],
                "antigen_gene": ["pp65"],
                "species": ["Human cytomegalovirus"],
                "mhc_allele": ["HLA-A*02:01"],
                "database": ["VDJdb"],
                "is_viral": [True],
            }
        )

        result = match_clonotypes(sample_clonotypes_df, db, match_by="CDR3ab")

        for col in ("db_protein", "db_mhc", "db_category", "db_match_strength"):
            assert col in result.columns

        matched = result[result["db_match"]]
        assert len(matched) >= 1
        first = matched.iloc[0]
        assert first["db_protein"] == "pp65"
        assert first["db_mhc"] == "HLA-A*02:01"
        assert first["db_epitope"] == "NLVPMVATV"
        assert first["db_category"] == "viral"
        assert first["db_match_strength"] == "ab"

    def test_match_strictness_strict_ab_no_fallback(self, sample_clonotypes_df):
        """strict_ab should *not* fall back to β-only matching. A clone
        whose β matches but α doesn't should come back unmatched."""
        # DB has the β CDR3 but a different α CDR3 → only β matches.
        db = pd.DataFrame(
            {
                "cdr3_beta": ["CASSLGQAYEQYF"],
                "cdr3_alpha": ["CAVDIFFERENT"],
                "epitope": ["NLV"],
                "antigen_gene": ["pp65"],
                "species": ["Human cytomegalovirus"],
                "mhc_allele": ["HLA-A*02:01"],
                "database": ["VDJdb"],
                "is_viral": [True],
            }
        )

        strict = match_clonotypes(
            sample_clonotypes_df, db, match_strictness="strict_ab"
        )
        # Legacy ab_with_partial would mark this clone as matched (partial);
        # strict_ab must not.
        assert strict["db_match"].sum() == 0
        assert not strict["db_match_partial"].any()

        partial = match_clonotypes(
            sample_clonotypes_df, db, match_strictness="ab_with_partial"
        )
        # The legacy mode picks up the β-only fallback.
        assert partial["db_match"].sum() >= 1
        assert partial["db_match_partial"].sum() >= 1
        # Partial matches are labeled b_only in db_match_strength.
        matched = partial[partial["db_match"]]
        assert (matched["db_match_strength"] == "b_only").all()

    def test_match_strictness_b_only_matches_strength(self, sample_clonotypes_df):
        """b_only strictness should mark all matches as db_match_strength=b_only,
        without setting db_match_partial (that flag is reserved for the αβ
        fallback path, not an explicit β-only run)."""
        db = pd.DataFrame(
            {
                "cdr3_beta": ["CASSLGQAYEQYF"],
                "cdr3_alpha": [None],
                "epitope": ["NLV"],
                "antigen_gene": ["pp65"],
                "species": ["Human cytomegalovirus"],
                "mhc_allele": ["HLA-A*02:01"],
                "database": ["VDJdb"],
                "is_viral": [True],
            }
        )

        result = match_clonotypes(
            sample_clonotypes_df, db, match_strictness="b_only"
        )
        matched = result[result["db_match"]]
        assert len(matched) >= 1
        assert (matched["db_match_strength"] == "b_only").all()
        # Not a fallback — explicit β-only mode.
        assert not matched["db_match_partial"].any()

    def test_match_strictness_invalid_raises(self, sample_clonotypes_df, sample_database_df):
        """An unknown strictness should produce a clear validation error."""
        with pytest.raises(TCRsiftValidationError, match="match_strictness"):
            match_clonotypes(
                sample_clonotypes_df,
                sample_database_df,
                match_strictness="bogus",
            )

    def test_match_strictness_overrides_match_by(self, sample_clonotypes_df):
        """When both knobs are set, ``match_strictness`` wins over the
        legacy ``match_by``. Locking this in so the precedence branch in
        match_clonotypes can't quietly regress."""
        # DB has only a β CDR3 match — α doesn't match.
        db = pd.DataFrame(
            {
                "cdr3_beta": ["CASSLGQAYEQYF"],
                "cdr3_alpha": ["CAVDIFFERENT"],
                "epitope": ["NLV"],
                "antigen_gene": ["pp65"],
                "species": ["Human cytomegalovirus"],
                "mhc_allele": ["HLA-A*02:01"],
                "database": ["VDJdb"],
                "is_viral": [True],
            }
        )
        # match_by="CDR3ab" (legacy) would fall back to β-only; pair it
        # with match_strictness="strict_ab" and the explicit param must
        # win → no match.
        result = match_clonotypes(
            sample_clonotypes_df,
            db,
            match_by="CDR3ab",
            match_strictness="strict_ab",
        )
        assert result["db_match"].sum() == 0

    def test_db_category_populated_from_species_and_protein(self, sample_clonotypes_df):
        """db_category should derive from species + antigen, with tumor-self
        antigens (MART-1) overriding the species-derived bucket."""
        # Two DB entries: one viral (CMV/pp65), one tumor-self (Homo sapiens,
        # but MART-1 should override "self" → "tumor_self").
        sample_a = sample_clonotypes_df.iloc[0]
        sample_b = sample_clonotypes_df.iloc[1] if len(sample_clonotypes_df) > 1 else sample_a
        db = pd.DataFrame(
            {
                "cdr3_beta": [sample_a["CDR3_beta"], sample_b["CDR3_beta"]],
                "cdr3_alpha": [sample_a["CDR3_alpha"], sample_b["CDR3_alpha"]],
                "epitope": ["NLV", "ELAGIGILTV"],
                "antigen_gene": ["pp65", "MLANA"],
                "species": ["Human cytomegalovirus", "Homo sapiens"],
                "mhc_allele": ["HLA-A*02:01", "HLA-A*02:01"],
                "database": ["VDJdb", "IEDB"],
                "is_viral": [True, False],
            }
        )

        result = match_clonotypes(sample_clonotypes_df, db, match_by="CDR3ab")
        matched = result[result["db_match"]]
        assert len(matched) >= 1
        # Categories should be assigned per the curated tables.
        cats = set(matched["db_category"].dropna())
        assert cats.issubset({"viral", "tumor_self"})

    def test_viral_flag_propagation(self, sample_clonotypes_df, sample_database_df):
        """Viral flag should propagate from database."""
        result = match_clonotypes(
            sample_clonotypes_df,
            sample_database_df,
            match_by="CDR3ab",
        )

        # Matched clones should inherit viral status
        matched = result[result["db_match"]]
        if len(matched) > 0:
            # All matches in our test db are viral
            assert all(matched["is_viral"])


class TestAnnotateClonotypes:
    """Tests for main annotation function."""

    @pytest.fixture
    def mock_db_paths(self, temp_dir):
        """Create mock database files using the real VDJdb paired schema."""
        vdjdb_path = temp_dir / "vdjdb_full.txt"
        pd.DataFrame(
            {
                "cdr3.alpha": ["CAVSDGGSQGNLIF"],
                "cdr3.beta": ["CASSLGQAYEQYF"],
                "antigen.epitope": ["NLV"],
                "antigen.species": ["CMV"],
            }
        ).to_csv(vdjdb_path, sep="\t", index=False)

        return {"vdjdb_path": vdjdb_path}

    def test_annotate_clonotypes(self, sample_clonotypes_df, mock_db_paths):
        """Full annotation pipeline."""
        result = annotate_clonotypes(
            sample_clonotypes_df,
            vdjdb_path=mock_db_paths["vdjdb_path"],
        )

        assert "db_match" in result.columns
        assert len(result) == len(sample_clonotypes_df)

    def test_exclude_viral(self, sample_clonotypes_df, mock_db_paths):
        """Exclude viral clones option."""
        result = annotate_clonotypes(
            sample_clonotypes_df,
            vdjdb_path=mock_db_paths["vdjdb_path"],
            exclude_viral=True,
        )

        # Viral clones should be excluded
        assert not any(result.get("is_viral", pd.Series([False])))

    def test_no_database_paths_returns_default_annotation_columns(self, sample_clonotypes_df):
        """Annotation should be optional when no database paths are provided.

        The short-circuit must produce the *full* annotation column set
        (legacy + structured columns from #48) so downstream code can
        rely on column presence regardless of whether annotation ran.
        """
        result = annotate_clonotypes(sample_clonotypes_df)

        assert len(result) == len(sample_clonotypes_df)
        expected_columns = {
            "db_match",
            "db_match_partial",
            "db_match_strength",
            "db_epitope",
            "db_protein",
            "db_species",
            "db_mhc",
            "db_category",
            "db_database",
            "is_viral",
        }
        assert expected_columns.issubset(result.columns)
        assert (~result["db_match"]).all()
        assert (~result["is_viral"]).all()


class TestGetAnnotationSummary:
    """Tests for annotation summary."""

    def test_summary_basic(self, sample_clonotypes_df, sample_database_df):
        """Summary should have basic stats."""
        annotated = match_clonotypes(sample_clonotypes_df, sample_database_df)
        summary = get_annotation_summary(annotated)

        assert "total" in summary
        assert "matched" in summary
        assert "viral" in summary

    def test_summary_database_breakdown(self, sample_clonotypes_df, sample_database_df):
        """Summary should have database breakdown."""
        annotated = match_clonotypes(sample_clonotypes_df, sample_database_df)
        summary = get_annotation_summary(annotated)

        if "database_breakdown" in summary:
            for db in ["VDJdb", "IEDB", "CEDAR"]:
                assert db in summary["database_breakdown"]
