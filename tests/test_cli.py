"""
Tests for CLI commands.
"""

import pandas as pd
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
import argparse

from tcrsift.cli import (
    cmd_annotate_gex,
    create_parser,
)
from tcrsift.sample_sheet import Sample, SampleSheet, load_sample_sheet


class TestAnnotateGexParser:
    """Tests for annotate-gex CLI parser configuration."""

    def test_parser_has_annotate_gex_command(self):
        """Parser should include annotate-gex command."""
        parser = create_parser()
        # Parse valid annotate-gex command to verify it exists
        args = parser.parse_args([
            "annotate-gex",
            "-i", "input.csv",
            "-o", "output.csv",
            "--gex-file", "matrix.h5",
        ])
        assert args.command == "annotate-gex"
        assert hasattr(args, 'func')

    def test_annotate_gex_required_args(self):
        """Test that required arguments are enforced."""
        parser = create_parser()

        # Should fail without required args
        with pytest.raises(SystemExit):
            parser.parse_args(["annotate-gex"])

    def test_annotate_gex_all_options(self):
        """Test all annotate-gex options parse correctly."""
        parser = create_parser()
        args = parser.parse_args([
            "annotate-gex",
            "-i", "input.csv",
            "-o", "output.csv",
            "--gex-file", "matrix.h5",
            "--barcode-col", "cell_barcode",
            "--genes", "CD3D,CD4,CD8A",
            "--prefix", "expr",
            "--no-qc",
            "--aggregate",
            "--group-col", "clone_id",
            "--cd4-cd8-counts",
            "--verbose",
        ])

        assert args.input == "input.csv"
        assert args.output == "output.csv"
        assert args.gex_file == "matrix.h5"
        assert args.barcode_col == "cell_barcode"
        assert args.genes == "CD3D,CD4,CD8A"
        assert args.prefix == "expr"
        assert args.no_qc is True
        assert args.aggregate is True
        assert args.group_col == "clone_id"
        assert args.cd4_cd8_counts is True
        assert args.verbose is True


class TestAnnotateGexCommand:
    """Tests for annotate-gex command function."""

    @pytest.fixture
    def cells_csv(self, temp_dir):
        """Create a cells CSV file for testing."""
        df = pd.DataFrame({
            "barcode": ["bc1", "bc2", "bc3", "bc4", "bc5"],
            "CDR3_pair": ["A/B", "A/B", "A/B", "C/D", "C/D"],
            "sample": ["S1", "S1", "S1", "S2", "S2"],
        })
        path = temp_dir / "cells.csv"
        df.to_csv(path, index=False)
        return path

    @pytest.fixture
    def cells_with_gex_csv(self, temp_dir):
        """Create a cells CSV file with GEX columns for testing."""
        df = pd.DataFrame({
            "barcode": ["bc1", "bc2", "bc3", "bc4", "bc5"],
            "CDR3_pair": ["A/B", "A/B", "A/B", "C/D", "C/D"],
            "gex.CD4": [10.0, 0.0, 5.0, 15.0, 0.0],
            "gex.CD8": [0.0, 20.0, 0.0, 0.0, 10.0],
        })
        path = temp_dir / "cells_with_gex.csv"
        df.to_csv(path, index=False)
        return path

    def test_annotate_gex_missing_barcode_col_warning(self, temp_dir, capsys):
        """Test warning when barcode column is missing."""
        # Create CSV without barcode column
        df = pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D"],
            "cell_count": [3, 2],
        })
        input_path = temp_dir / "no_barcode.csv"
        df.to_csv(input_path, index=False)
        output_path = temp_dir / "output.csv"

        args = argparse.Namespace(
            input=str(input_path),
            output=str(output_path),
            gex_file="fake.h5",
            barcode_col="barcode",
            genes=None,
            prefix="gex",
            no_qc=False,
            aggregate=False,
            group_col="CDR3_pair",
            cd4_cd8_counts=False,
            verbose=False,
        )

        cmd_annotate_gex(args)

        captured = capsys.readouterr()
        assert "Warning: Barcode column 'barcode' not found" in captured.out

    def test_annotate_gex_aggregate_only(self, cells_with_gex_csv, temp_dir):
        """Test aggregation without GEX augmentation."""
        output_path = temp_dir / "aggregated.csv"

        args = argparse.Namespace(
            input=str(cells_with_gex_csv),
            output=str(output_path),
            gex_file="fake.h5",
            barcode_col="nonexistent",  # Will skip augmentation
            genes=None,
            prefix="gex",
            no_qc=False,
            aggregate=True,
            group_col="CDR3_pair",
            cd4_cd8_counts=False,
            verbose=False,
        )

        cmd_annotate_gex(args)

        result = pd.read_csv(output_path)
        assert len(result) == 2  # 2 clonotypes
        assert "CDR3_pair" in result.columns
        assert "total_cells.count" in result.columns
        assert result[result["CDR3_pair"] == "A/B"]["total_cells.count"].iloc[0] == 3

    def test_annotate_gex_cd4_cd8_counts_without_augmentation(self, temp_dir, capsys):
        """Test CD4/CD8 counts warning when no augmentation."""
        # Create CSV without GEX columns
        df = pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D"],
            "cell_count": [3, 2],
        })
        input_path = temp_dir / "no_gex.csv"
        df.to_csv(input_path, index=False)
        output_path = temp_dir / "output.csv"

        args = argparse.Namespace(
            input=str(input_path),
            output=str(output_path),
            gex_file="fake.h5",
            barcode_col="barcode",  # Missing
            genes=None,
            prefix="gex",
            no_qc=False,
            aggregate=False,
            group_col="CDR3_pair",
            cd4_cd8_counts=True,  # Requested but can't do it
            verbose=False,
        )

        cmd_annotate_gex(args)

        captured = capsys.readouterr()
        assert "Cannot compute CD4/CD8 counts without GEX augmentation" in captured.out

    def test_annotate_gex_custom_genes_parsing(self, capsys, temp_dir):
        """Test custom gene list parsing."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B"],
            "cell_count": [1],
        })
        input_path = temp_dir / "input.csv"
        df.to_csv(input_path, index=False)
        output_path = temp_dir / "output.csv"

        args = argparse.Namespace(
            input=str(input_path),
            output=str(output_path),
            gex_file="fake.h5",
            barcode_col="barcode",
            genes="GENE1, GENE2, GENE3",  # With spaces
            prefix="gex",
            no_qc=False,
            aggregate=False,
            group_col="CDR3_pair",
            cd4_cd8_counts=False,
            verbose=False,
        )

        cmd_annotate_gex(args)

        captured = capsys.readouterr()
        assert "Using custom gene list: ['GENE1', 'GENE2', 'GENE3']" in captured.out


class TestAutoDetectTilSamples:
    """Tests for auto-detection of TIL samples in run command."""

    @pytest.fixture
    def sample_sheet_with_til(self, temp_dir):
        """Create a sample sheet with TIL samples."""
        yaml_content = """
samples:
  - sample: "Culture_Pool1"
    vdj_dir: "/data/culture/vdj"
    source: "culture"
  - sample: "Culture_Pool2"
    vdj_dir: "/data/culture2/vdj"
    source: "culture"
  - sample: "Patient1_TIL"
    vdj_dir: "/data/til/vdj"
    source: "til"
  - sample: "Patient1_TIL2"
    vdj_dir: "/data/til2/vdj"
    source: "til"
"""
        path = temp_dir / "samples_with_til.yaml"
        path.write_text(yaml_content)
        return path

    @pytest.fixture
    def sample_sheet_without_til(self, temp_dir):
        """Create a sample sheet without TIL samples."""
        yaml_content = """
samples:
  - sample: "Culture_Pool1"
    vdj_dir: "/data/culture/vdj"
    source: "culture"
  - sample: "Tetramer_Sample"
    vdj_dir: "/data/tetramer/vdj"
    source: "tetramer"
"""
        path = temp_dir / "samples_no_til.yaml"
        path.write_text(yaml_content)
        return path

    def test_get_til_samples_detects_til(self, sample_sheet_with_til):
        """Test that TIL samples are detected from sample sheet."""
        sample_sheet = load_sample_sheet(sample_sheet_with_til)
        til_samples = sample_sheet.get_til_samples()

        assert len(til_samples) == 2
        til_names = [s.sample for s in til_samples]
        assert "Patient1_TIL" in til_names
        assert "Patient1_TIL2" in til_names

    def test_get_til_samples_empty_when_no_til(self, sample_sheet_without_til):
        """Test that empty list returned when no TIL samples."""
        sample_sheet = load_sample_sheet(sample_sheet_without_til)
        til_samples = sample_sheet.get_til_samples()

        assert len(til_samples) == 0

    def test_sample_is_til_method(self):
        """Test Sample.is_til() method."""
        til_sample = Sample(sample="TIL", vdj_dir="/path", source="til")
        culture_sample = Sample(sample="Culture", vdj_dir="/path", source="culture")
        no_source_sample = Sample(sample="Unknown", vdj_dir="/path")

        assert til_sample.is_til() is True
        assert culture_sample.is_til() is False
        assert no_source_sample.is_til() is False

    def test_sample_sheet_get_culture_samples(self, sample_sheet_with_til):
        """Test filtering culture samples."""
        sample_sheet = load_sample_sheet(sample_sheet_with_til)
        culture_samples = sample_sheet.get_culture_samples()

        assert len(culture_samples) == 2
        assert all(s.source == "culture" for s in culture_samples)

    def test_auto_detect_til_logic(self, sample_sheet_with_til):
        """Test the auto-detection logic used in run command."""
        # This mimics the logic in cmd_run
        sample_sheet = load_sample_sheet(sample_sheet_with_til)
        til_samples = sample_sheet.get_til_samples()

        if til_samples:
            til_sample_names = [s.sample for s in til_samples]
        else:
            til_sample_names = []

        assert til_sample_names == ["Patient1_TIL", "Patient1_TIL2"]


class TestSampleSheetSourceTypes:
    """Additional tests for sample sheet source type detection."""

    def test_tetramer_source(self):
        """Test tetramer source detection."""
        sample = Sample(sample="Tet", vdj_dir="/path", source="tetramer")
        assert sample.is_tetramer_or_sct() is True
        assert sample.is_til() is False

    def test_sct_source(self):
        """Test SCT source detection."""
        sample = Sample(sample="SCT", vdj_dir="/path", source="sct")
        assert sample.is_tetramer_or_sct() is True
        assert sample.is_til() is False

    def test_til_antigen_type_does_not_affect_til_check(self):
        """Test that is_til only checks source, not antigen_type."""
        # TIL identification should be by source, not antigen type
        sample = Sample(sample="S1", vdj_dir="/path", antigen_type="short_peptide", source="til")
        assert sample.is_til() is True

    def test_get_tetramer_samples(self):
        """Test filtering tetramer and SCT samples."""
        ss = SampleSheet(samples=[
            Sample(sample="C1", vdj_dir="/path1", source="culture"),
            Sample(sample="T1", vdj_dir="/path2", source="til"),
            Sample(sample="Tet1", vdj_dir="/path3", source="tetramer"),
            Sample(sample="SCT1", vdj_dir="/path4", source="sct"),
        ])
        tetramer = ss.get_tetramer_samples()

        assert len(tetramer) == 2
        assert all(s.is_tetramer_or_sct() for s in tetramer)
