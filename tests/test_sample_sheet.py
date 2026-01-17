"""
Tests for sample sheet parsing and validation.
"""

import pytest
import pandas as pd
from pathlib import Path

from tcrsift.sample_sheet import (
    Sample,
    SampleSheet,
    load_sample_sheet,
    validate_sample_sheet,
    ANTIGEN_TYPE_TCELL_EXPECTATIONS,
    VALID_ANTIGEN_TYPES,
    VALID_SOURCES,
)


class TestSample:
    """Tests for the Sample dataclass."""

    def test_sample_with_gex_dir(self):
        """Sample with only gex_dir should be valid."""
        sample = Sample(sample="S1", gex_dir="/path/to/gex")
        assert sample.sample == "S1"
        assert sample.gex_dir == "/path/to/gex"
        assert sample.vdj_dir is None

    def test_sample_with_vdj_dir(self):
        """Sample with only vdj_dir should be valid."""
        sample = Sample(sample="S1", vdj_dir="/path/to/vdj")
        assert sample.sample == "S1"
        assert sample.vdj_dir == "/path/to/vdj"
        assert sample.gex_dir is None

    def test_sample_with_both_dirs(self):
        """Sample with both dirs should be valid."""
        sample = Sample(sample="S1", gex_dir="/path/to/gex", vdj_dir="/path/to/vdj")
        assert sample.gex_dir == "/path/to/gex"
        assert sample.vdj_dir == "/path/to/vdj"

    def test_sample_requires_data_dir(self):
        """Sample without gex_dir, vdj_dir, or sct_path should raise error."""
        with pytest.raises(ValueError, match="must have at least gex_dir, vdj_dir, or sct_path"):
            Sample(sample="S1")

    def test_invalid_antigen_type_raises(self):
        """Invalid antigen type should raise error."""
        with pytest.raises(ValueError, match="Invalid antigen_type"):
            Sample(sample="S1", vdj_dir="/path", antigen_type="invalid_type")

    def test_valid_antigen_types(self):
        """All valid antigen types should work."""
        for antigen_type in VALID_ANTIGEN_TYPES:
            sample = Sample(sample="S1", vdj_dir="/path", antigen_type=antigen_type)
            assert sample.antigen_type == antigen_type

    def test_invalid_source_raises(self):
        """Invalid source should raise error."""
        with pytest.raises(ValueError, match="Invalid source"):
            Sample(sample="S1", vdj_dir="/path", source="invalid_source")

    def test_valid_sources(self):
        """All valid sources should work."""
        for source in VALID_SOURCES:
            sample = Sample(sample="S1", vdj_dir="/path", source=source)
            assert sample.source == source

    def test_invalid_tcell_type_expected_raises(self):
        """Invalid tcell_type_expected should raise error."""
        with pytest.raises(ValueError, match="Invalid tcell_type_expected"):
            Sample(sample="S1", vdj_dir="/path", tcell_type_expected="invalid")

    def test_invalid_pre_sorted_raises(self):
        """Invalid pre_sorted should raise error."""
        with pytest.raises(ValueError, match="Invalid pre_sorted"):
            Sample(sample="S1", vdj_dir="/path", pre_sorted="invalid")

    def test_invalid_mhc_blocking_raises(self):
        """Invalid mhc_blocking should raise error."""
        with pytest.raises(ValueError, match="Invalid mhc_blocking"):
            Sample(sample="S1", vdj_dir="/path", mhc_blocking="invalid")


class TestSampleGetExpectedTcellType:
    """Tests for get_expected_tcell_type method."""

    def test_direct_specification_takes_priority(self):
        """Direct tcell_type_expected should override all."""
        sample = Sample(
            sample="S1",
            vdj_dir="/path",
            antigen_type="short_peptide",  # would give CD8
            tcell_type_expected="CD4",
        )
        assert sample.get_expected_tcell_type() == "CD4"

    def test_pre_sorted_overrides_antigen_type(self):
        """Pre-sorting should override antigen type."""
        sample = Sample(
            sample="S1",
            vdj_dir="/path",
            antigen_type="long_peptide",  # would give mixed
            pre_sorted="CD8",
        )
        assert sample.get_expected_tcell_type() == "CD8"

    def test_mhc_blocking_inference(self):
        """MHC blocking should infer opposite T cell type."""
        # MHC-I blocking means CD8 responses blocked -> expect CD4
        sample1 = Sample(sample="S1", vdj_dir="/path", mhc_blocking="MHC-I")
        assert sample1.get_expected_tcell_type() == "CD4"

        # MHC-II blocking means CD4 responses blocked -> expect CD8
        sample2 = Sample(sample="S2", vdj_dir="/path", mhc_blocking="MHC-II")
        assert sample2.get_expected_tcell_type() == "CD8"

    def test_antigen_type_expectations(self):
        """Antigen type should give expected T cell type."""
        for antigen_type, expected in ANTIGEN_TYPE_TCELL_EXPECTATIONS.items():
            sample = Sample(sample="S1", vdj_dir="/path", antigen_type=antigen_type)
            assert sample.get_expected_tcell_type() == expected

    def test_no_expectation(self):
        """Sample without hints should return None."""
        sample = Sample(sample="S1", vdj_dir="/path")
        assert sample.get_expected_tcell_type() is None


class TestSampleTypeChecks:
    """Tests for is_tetramer_or_sct and is_til methods."""

    def test_is_tetramer_by_source(self):
        """Source=tetramer should be detected."""
        sample = Sample(sample="S1", vdj_dir="/path", source="tetramer")
        assert sample.is_tetramer_or_sct() is True
        assert sample.is_til() is False

    def test_is_sct_by_source(self):
        """Source=sct should be detected."""
        sample = Sample(sample="S1", vdj_dir="/path", source="sct")
        assert sample.is_tetramer_or_sct() is True

    def test_is_tetramer_by_antigen_type(self):
        """Tetramer and SCT antigen types should be detected."""
        sample1 = Sample(sample="S1", vdj_dir="/path", antigen_type="tetramer_mhc1")
        assert sample1.is_tetramer_or_sct() is True

        sample2 = Sample(sample="S2", vdj_dir="/path", antigen_type="tetramer_mhc2")
        assert sample2.is_tetramer_or_sct() is True

        sample3 = Sample(sample="S3", vdj_dir="/path", antigen_type="sct")
        assert sample3.is_tetramer_or_sct() is True

    def test_is_til(self):
        """TIL source should be detected."""
        sample = Sample(sample="S1", vdj_dir="/path", source="til")
        assert sample.is_til() is True
        assert sample.is_tetramer_or_sct() is False

    def test_culture_sample(self):
        """Culture sample should not be tetramer or TIL."""
        sample = Sample(sample="S1", vdj_dir="/path", source="culture")
        assert sample.is_tetramer_or_sct() is False
        assert sample.is_til() is False


class TestSampleSheet:
    """Tests for SampleSheet class."""

    def test_len(self):
        """Test __len__."""
        ss = SampleSheet(samples=[
            Sample(sample="S1", vdj_dir="/path1"),
            Sample(sample="S2", vdj_dir="/path2"),
        ])
        assert len(ss) == 2

    def test_iter(self):
        """Test __iter__."""
        samples = [
            Sample(sample="S1", vdj_dir="/path1"),
            Sample(sample="S2", vdj_dir="/path2"),
        ]
        ss = SampleSheet(samples=samples)
        assert list(ss) == samples

    def test_getitem(self):
        """Test __getitem__."""
        samples = [
            Sample(sample="S1", vdj_dir="/path1"),
            Sample(sample="S2", vdj_dir="/path2"),
        ]
        ss = SampleSheet(samples=samples)
        assert ss[0].sample == "S1"
        assert ss[1].sample == "S2"

    def test_get_sample(self):
        """Test get_sample by name."""
        ss = SampleSheet(samples=[
            Sample(sample="S1", vdj_dir="/path1"),
            Sample(sample="S2", vdj_dir="/path2"),
        ])
        assert ss.get_sample("S1").sample == "S1"
        assert ss.get_sample("S2").sample == "S2"
        assert ss.get_sample("S3") is None

    def test_get_culture_samples(self):
        """Test filtering culture samples."""
        ss = SampleSheet(samples=[
            Sample(sample="C1", vdj_dir="/path1", source="culture"),
            Sample(sample="T1", vdj_dir="/path2", source="til"),
            Sample(sample="C2", vdj_dir="/path3", source="culture"),
        ])
        culture = ss.get_culture_samples()
        assert len(culture) == 2
        assert all(s.source == "culture" for s in culture)

    def test_get_til_samples(self):
        """Test filtering TIL samples."""
        ss = SampleSheet(samples=[
            Sample(sample="C1", vdj_dir="/path1", source="culture"),
            Sample(sample="T1", vdj_dir="/path2", source="til"),
        ])
        til = ss.get_til_samples()
        assert len(til) == 1
        assert til[0].sample == "T1"

    def test_get_tetramer_samples(self):
        """Test filtering tetramer samples."""
        ss = SampleSheet(samples=[
            Sample(sample="C1", vdj_dir="/path1", source="culture"),
            Sample(sample="Tet1", vdj_dir="/path2", source="tetramer"),
            Sample(sample="SCT1", vdj_dir="/path3", source="sct"),
        ])
        tetramer = ss.get_tetramer_samples()
        assert len(tetramer) == 2

    def test_to_dataframe(self):
        """Test conversion to DataFrame."""
        ss = SampleSheet(samples=[
            Sample(sample="S1", vdj_dir="/path1", antigen_type="short_peptide"),
            Sample(sample="S2", vdj_dir="/path2", source="til"),
        ])
        df = ss.to_dataframe()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert "sample" in df.columns
        assert "expected_tcell_type" in df.columns
        assert df.loc[0, "expected_tcell_type"] == "CD8"


class TestLoadSampleSheet:
    """Tests for loading sample sheets."""

    def test_load_csv(self, sample_csv_sample_sheet):
        """Test loading CSV sample sheet."""
        ss = load_sample_sheet(sample_csv_sample_sheet)
        assert len(ss) == 2
        assert ss[0].sample == "S1"
        assert ss[0].antigen_type == "short_peptide"

    def test_load_yaml(self, sample_yaml_sample_sheet):
        """Test loading YAML sample sheet."""
        ss = load_sample_sheet(sample_yaml_sample_sheet)
        assert len(ss) == 2
        assert ss[0].sample == "S1"

    def test_nonexistent_file_raises(self):
        """Test loading nonexistent file raises error."""
        with pytest.raises(FileNotFoundError):
            load_sample_sheet("/nonexistent/path.csv")

    def test_unsupported_format_raises(self, temp_dir):
        """Test unsupported format raises error."""
        bad_file = temp_dir / "samples.json"
        bad_file.write_text("{}")
        with pytest.raises(ValueError, match="Unsupported sample sheet format"):
            load_sample_sheet(bad_file)


class TestValidateSampleSheet:
    """Tests for sample sheet validation."""

    def test_duplicate_sample_names(self):
        """Test detection of duplicate sample names."""
        ss = SampleSheet(samples=[
            Sample(sample="S1", vdj_dir="/path1"),
            Sample(sample="S1", vdj_dir="/path2"),  # duplicate
        ])
        warnings = validate_sample_sheet(ss)
        assert any("Duplicate sample names" in w for w in warnings)

    def test_conflicting_long_peptide_cd8(self):
        """Test warning for long peptide expecting CD8."""
        ss = SampleSheet(samples=[
            Sample(
                sample="S1",
                vdj_dir="/path",
                antigen_type="long_peptide",
                tcell_type_expected="CD8",
            ),
        ])
        warnings = validate_sample_sheet(ss)
        assert any("Long peptides typically favor CD4" in w for w in warnings)

    def test_conflicting_short_peptide_cd4(self):
        """Test warning for short peptide expecting CD4."""
        ss = SampleSheet(samples=[
            Sample(
                sample="S1",
                vdj_dir="/path",
                antigen_type="short_peptide",
                tcell_type_expected="CD4",
            ),
        ])
        warnings = validate_sample_sheet(ss)
        assert any("Short peptides typically bind MHC-I" in w for w in warnings)

    def test_nonexistent_paths_warned(self, temp_dir):
        """Test warning for nonexistent paths."""
        ss = SampleSheet(samples=[
            Sample(sample="S1", gex_dir="/nonexistent/gex", vdj_dir="/nonexistent/vdj"),
        ])
        warnings = validate_sample_sheet(ss)
        assert any("gex_dir does not exist" in w for w in warnings)
        assert any("vdj_dir does not exist" in w for w in warnings)

    def test_valid_sample_sheet_no_warnings(self, temp_dir):
        """Valid sample sheet should have no warnings about conflicts."""
        gex_dir = temp_dir / "gex"
        vdj_dir = temp_dir / "vdj"
        gex_dir.mkdir()
        vdj_dir.mkdir()

        ss = SampleSheet(samples=[
            Sample(
                sample="S1",
                gex_dir=str(gex_dir),
                vdj_dir=str(vdj_dir),
                antigen_type="short_peptide",
            ),
        ])
        warnings = validate_sample_sheet(ss)
        # Should only have no warnings (paths exist, consistent metadata)
        assert len(warnings) == 0
