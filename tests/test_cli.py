"""
Tests for CLI commands.
"""

import argparse
from unittest.mock import patch

import pandas as pd
import pytest

from tcrsift.cli import (
    cmd_annotate_gex,
    cmd_run,
    create_parser,
)
from tcrsift.sample_sheet import Sample, SampleSheet, load_sample_sheet
from tcrsift.validation import (
    TCRsiftValidationError,
    validate_annotate_gex_args,
    validate_assemble_args,
    validate_cli_conditional_requirement,
    validate_cli_mutually_exclusive,
    validate_run_args,
)


class TestAnnotateGexParser:
    """Tests for annotate-gex CLI parser configuration."""

    def test_parser_has_annotate_gex_command(self):
        """Parser should include annotate-gex command."""
        parser = create_parser()
        # Parse valid annotate-gex command to verify it exists
        args = parser.parse_args(
            [
                "annotate-gex",
                "-i",
                "input.csv",
                "-o",
                "output.csv",
                "--gex-file",
                "matrix.h5",
            ]
        )
        assert args.command == "annotate-gex"
        assert hasattr(args, "func")

    def test_annotate_gex_required_args(self):
        """Test that required arguments are enforced."""
        parser = create_parser()

        # Should fail without required args
        with pytest.raises(SystemExit):
            parser.parse_args(["annotate-gex"])

    def test_annotate_gex_all_options(self):
        """Test all annotate-gex options parse correctly."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "annotate-gex",
                "-i",
                "input.csv",
                "-o",
                "output.csv",
                "--gex-file",
                "matrix.h5",
                "--barcode-col",
                "cell_barcode",
                "--genes",
                "CD3D,CD4,CD8A",
                "--prefix",
                "expr",
                "--no-qc",
                "--aggregate",
                "--group-col",
                "CDR3ab",
                "--cd4-cd8-counts",
                "--verbose",
            ]
        )

        assert args.input == "input.csv"
        assert args.output == "output.csv"
        assert args.gex_file == "matrix.h5"
        assert args.barcode_col == "cell_barcode"
        assert args.genes == "CD3D,CD4,CD8A"
        assert args.prefix == "expr"
        assert args.no_qc is True
        assert args.aggregate is True
        assert args.group_col == "CDR3ab"
        assert args.cd4_cd8_counts is True
        assert args.verbose is True


class TestAnnotateParser:
    """Tests for annotate CLI parser configuration."""

    def test_annotate_match_by_option(self):
        """Annotate command should accept --match-by."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "annotate",
                "-i",
                "input.csv",
                "-o",
                "output.csv",
                "--match-by",
                "CDR3b_only",
            ]
        )
        assert args.match_by == "CDR3b_only"


class TestDataDownloadParser:
    """Tests for `tcrsift data download` accepting multiple databases (#44).

    The original argument was a single-value `--db`, so repeating it
    (`--db vdjdb --db iedb`) silently downloaded only the last one.
    Both space-separated (`--db vdjdb iedb`) and repeated-flag forms
    must produce the same list.
    """

    def test_data_download_no_db_flag_yields_none(self):
        parser = create_parser()
        args = parser.parse_args(["data", "download"])
        assert args.db is None  # handler interprets as "all"

    def test_data_download_space_separated_dbs(self):
        parser = create_parser()
        args = parser.parse_args(["data", "download", "--db", "vdjdb", "iedb"])
        assert args.db == ["vdjdb", "iedb"]

    def test_data_download_repeated_flag_dbs(self):
        """Reproduces the original bug surface: `--db vdjdb --db iedb`
        used to drop the first value. With `action=extend`, both end
        up in the list (#44)."""
        parser = create_parser()
        args = parser.parse_args(
            ["data", "download", "--db", "vdjdb", "--db", "iedb"]
        )
        assert args.db == ["vdjdb", "iedb"]

    def test_data_download_mixed_forms(self):
        parser = create_parser()
        args = parser.parse_args(
            ["data", "download", "--db", "vdjdb", "iedb", "--db", "cedar"]
        )
        assert args.db == ["vdjdb", "iedb", "cedar"]

    def test_data_clear_repeated_flag_dbs(self):
        """Same fix applied to `data clear --db` (#44 bonus)."""
        parser = create_parser()
        args = parser.parse_args(
            ["data", "clear", "--db", "vdjdb", "--db", "iedb"]
        )
        assert args.db == ["vdjdb", "iedb"]

    def test_data_download_unknown_db_rejected(self):
        """`choices=` constraint still applies — unknown names fail
        argparse validation rather than silently doing nothing."""
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["data", "download", "--db", "bogus_db"])


class TestDataDownloadCommand:
    """Tests for the cmd_data_download handler behavior (#44)."""

    def test_download_iterates_each_db(self, monkeypatch, tmp_path):
        """Multiple `--db` values must trigger one download call each,
        in the order given. Previously only the last value reached the
        download path."""
        from tcrsift.cli import cmd_data_download

        called: list[str] = []

        def fake_download(name, data_dir=None, force=False):
            called.append(name)
            return tmp_path / name / "fake"

        monkeypatch.setattr("tcrsift.datacache.download_database", fake_download)

        args = argparse.Namespace(
            db=["vdjdb", "iedb"], cache_dir=str(tmp_path), force=False
        )
        cmd_data_download(args)

        assert called == ["vdjdb", "iedb"]

    def test_download_deduplicates_repeats(self, monkeypatch, tmp_path):
        """A name repeated on the command line is downloaded only once.
        Order of first mention is preserved."""
        from tcrsift.cli import cmd_data_download

        called: list[str] = []

        def fake_download(name, data_dir=None, force=False):
            called.append(name)
            return tmp_path / name / "fake"

        monkeypatch.setattr("tcrsift.datacache.download_database", fake_download)

        args = argparse.Namespace(
            db=["vdjdb", "iedb", "vdjdb"], cache_dir=str(tmp_path), force=False
        )
        cmd_data_download(args)

        assert called == ["vdjdb", "iedb"]

    def test_download_no_db_downloads_all_with_url(self, monkeypatch, tmp_path):
        """Default behavior (no --db) still downloads every database
        spec that has a non-None download_url."""
        from tcrsift.cli import cmd_data_download
        from tcrsift.datacache import DATABASES

        called: list[str] = []

        def fake_download(name, data_dir=None, force=False):
            called.append(name)
            return tmp_path / name / "fake"

        monkeypatch.setattr("tcrsift.datacache.download_database", fake_download)

        args = argparse.Namespace(db=None, cache_dir=str(tmp_path), force=False)
        cmd_data_download(args)

        expected = [
            name for name, spec in DATABASES.items() if spec.download_url is not None
        ]
        assert called == expected


class TestAnnotateGexCommand:
    """Tests for annotate-gex command function."""

    @pytest.fixture
    def cells_csv(self, temp_dir):
        """Create a cells CSV file for testing."""
        df = pd.DataFrame(
            {
                "barcode": ["bc1", "bc2", "bc3", "bc4", "bc5"],
                "CDR3_pair": ["A/B", "A/B", "A/B", "C/D", "C/D"],
                "sample": ["S1", "S1", "S1", "S2", "S2"],
            }
        )
        path = temp_dir / "cells.csv"
        df.to_csv(path, index=False)
        return path

    @pytest.fixture
    def cells_with_gex_csv(self, temp_dir):
        """Create a cells CSV file with GEX columns for testing."""
        df = pd.DataFrame(
            {
                "barcode": ["bc1", "bc2", "bc3", "bc4", "bc5"],
                "CDR3_pair": ["A/B", "A/B", "A/B", "C/D", "C/D"],
                "gex.CD4": [10.0, 0.0, 5.0, 15.0, 0.0],
                "gex.CD8": [0.0, 20.0, 0.0, 0.0, 10.0],
            }
        )
        path = temp_dir / "cells_with_gex.csv"
        df.to_csv(path, index=False)
        return path

    @pytest.fixture
    def fake_gex_file(self, temp_dir):
        """Create a fake GEX file that passes validation."""
        path = temp_dir / "fake_matrix.h5"
        path.write_bytes(b"fake h5 content")
        return path

    def test_annotate_gex_missing_barcode_col_warning(self, temp_dir, fake_gex_file, capsys):
        """Test warning when barcode column is missing."""
        # Create CSV without barcode column
        df = pd.DataFrame(
            {
                "CDR3_pair": ["A/B", "C/D"],
                "cell_count": [3, 2],
            }
        )
        input_path = temp_dir / "no_barcode.csv"
        df.to_csv(input_path, index=False)
        output_path = temp_dir / "output.csv"

        args = argparse.Namespace(
            input=str(input_path),
            output=str(output_path),
            gex_file=str(fake_gex_file),
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

    def test_annotate_gex_aggregate_only(self, cells_with_gex_csv, temp_dir, fake_gex_file):
        """Test aggregation without GEX augmentation."""
        output_path = temp_dir / "aggregated.csv"

        args = argparse.Namespace(
            input=str(cells_with_gex_csv),
            output=str(output_path),
            gex_file=str(fake_gex_file),
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

    def test_annotate_gex_cd4_cd8_counts_without_augmentation(
        self, temp_dir, fake_gex_file, capsys
    ):
        """Test CD4/CD8 counts warning when no augmentation."""
        # Create CSV without GEX columns
        df = pd.DataFrame(
            {
                "CDR3_pair": ["A/B", "C/D"],
                "cell_count": [3, 2],
            }
        )
        input_path = temp_dir / "no_gex.csv"
        df.to_csv(input_path, index=False)
        output_path = temp_dir / "output.csv"

        args = argparse.Namespace(
            input=str(input_path),
            output=str(output_path),
            gex_file=str(fake_gex_file),
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

    def test_annotate_gex_custom_genes_parsing(self, capsys, temp_dir, fake_gex_file):
        """Test custom gene list parsing."""
        df = pd.DataFrame(
            {
                "CDR3_pair": ["A/B"],
                "cell_count": [1],
            }
        )
        input_path = temp_dir / "input.csv"
        df.to_csv(input_path, index=False)
        output_path = temp_dir / "output.csv"

        args = argparse.Namespace(
            input=str(input_path),
            output=str(output_path),
            gex_file=str(fake_gex_file),
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


class TestRunCommandTilSamples:
    """Tests for run command behavior with explicit TIL samples."""

    def test_run_with_explicit_til_samples_does_not_crash(
        self, tmp_path, monkeypatch
    ):
        """cmd_run should work when til_samples are explicitly set."""
        import anndata as ad
        import pandas as pd

        # Sample sheet with culture + TIL
        sample_sheet_path = tmp_path / "samples.yaml"
        sample_sheet_path.write_text(
            """
samples:
  - sample: "Culture1"
    vdj_dir: "/data/culture/vdj"
    source: "culture"
  - sample: "Patient1_TIL"
    vdj_dir: "/data/til/vdj"
    source: "til"
"""
        )

        # Minimal AnnData returned by patched load_samples
        obs = pd.DataFrame(
            {
                "sample": ["Culture1"],
                "source": ["culture"],
            },
            index=["cell1"],
        )
        adata = ad.AnnData(obs=obs)

        def fake_load_samples(*_args, **_kwargs):
            return adata

        def fake_phenotype_cells(adata_in, *_args, **_kwargs):
            return adata_in

        def fake_filter_by_tcell_type(adata_in, *_args, **_kwargs):
            return adata_in

        def fake_aggregate_clonotypes(*_args, **_kwargs):
            return pd.DataFrame(
                {
                    "CDR3ab": ["A_B"],
                    "CDR3_alpha": ["A"],
                    "CDR3_beta": ["B"],
                    "cell_count": [1],
                }
            )

        def fake_filter_clonotypes(df, *_args, **_kwargs):
            return df.assign(tier="tier1")

        def fake_split_by_tier(df, *_args, **_kwargs):
            return {"tier1": df}

        def fake_load_til_samples(*_args, **_kwargs):
            return {}

        # Patch heavy functions
        monkeypatch.setattr("tcrsift.loader.load_samples", fake_load_samples)
        monkeypatch.setattr("tcrsift.phenotype.phenotype_cells", fake_phenotype_cells)
        monkeypatch.setattr("tcrsift.phenotype.filter_by_tcell_type", fake_filter_by_tcell_type)
        monkeypatch.setattr("tcrsift.clonotype.aggregate_clonotypes", fake_aggregate_clonotypes)
        monkeypatch.setattr("tcrsift.filter.filter_clonotypes", fake_filter_clonotypes)
        monkeypatch.setattr("tcrsift.filter.split_by_tier", fake_split_by_tier)
        monkeypatch.setattr("tcrsift.til.load_til_samples", fake_load_til_samples)

        output_dir = tmp_path / "out"

        args = argparse.Namespace(
            sample_sheet=str(sample_sheet_path),
            output_dir=str(output_dir),
            config=None,
            # Disable optional outputs
            generate_plots=False,
            generate_report=False,
            # Ensure assembly is skipped
            no_leaders=True,
            include_constant=False,
            single_chain=False,
            # Explicit TIL samples
            til_samples="Patient1_TIL",
            # Required defaults for other args
            min_genes=None,
            max_genes=None,
            min_counts=None,
            max_counts=None,
            min_mito_pct=None,
            max_mito_pct=None,
            cd4_cd8_ratio=None,
            min_cd3_reads=None,
            group_by=None,
            handle_doublets=None,
            min_umi=None,
            tcell_type=None,
            method=None,
            min_cells=None,
            min_frequency=None,
            require_complete=None,
            fdr_tiers=None,
            vdjdb_path=None,
            iedb_path=None,
            cedar_path=None,
            match_by=None,
            exclude_viral=None,
            flag_only=None,
            til_match_by=None,
            min_til_cells=None,
            alpha_leader=None,
            beta_leader=None,
            leaders_from_contigs=False,
            contigs_dir=None,
            linker=None,
            constant_source=None,
            skip_plots=None,
            verbose=False,
        )

        # Should not raise (previously UnboundLocalError)
        cmd_run(args)

    def test_run_excludes_til_cells_from_clonotyping(self, tmp_path, monkeypatch):
        """Culture clonotyping should exclude TIL cells when present in sample sheet."""
        import anndata as ad
        import pandas as pd

        sample_sheet_path = tmp_path / "samples.yaml"
        sample_sheet_path.write_text(
            """
samples:
  - sample: "Culture1"
    vdj_dir: "/data/culture/vdj"
    source: "culture"
  - sample: "Patient1_TIL"
    vdj_dir: "/data/til/vdj"
    source: "til"
"""
        )

        obs = pd.DataFrame(
            {
                "sample": ["Culture1", "Patient1_TIL"],
                "source": ["culture", "til"],
            },
            index=["cell1", "cell2"],
        )
        adata = ad.AnnData(obs=obs)

        def fake_load_samples(*_args, **_kwargs):
            return adata

        def fake_phenotype_cells(adata_in, *_args, **_kwargs):
            return adata_in

        def fake_filter_by_tcell_type(adata_in, *_args, **_kwargs):
            return adata_in

        def fake_aggregate_clonotypes(adata_in, *_args, **_kwargs):
            # Ensure TIL cells were filtered out before clonotyping
            assert "til" not in adata_in.obs["source"].unique()
            return pd.DataFrame(
                {
                    "CDR3ab": ["A_B"],
                    "CDR3_alpha": ["A"],
                    "CDR3_beta": ["B"],
                    "cell_count": [1],
                }
            )

        def fake_filter_clonotypes(df, *_args, **_kwargs):
            return df.assign(tier="tier1")

        def fake_split_by_tier(df, *_args, **_kwargs):
            return {"tier1": df}

        # Patch heavy functions
        monkeypatch.setattr("tcrsift.loader.load_samples", fake_load_samples)
        monkeypatch.setattr("tcrsift.phenotype.phenotype_cells", fake_phenotype_cells)
        monkeypatch.setattr("tcrsift.phenotype.filter_by_tcell_type", fake_filter_by_tcell_type)
        monkeypatch.setattr("tcrsift.clonotype.aggregate_clonotypes", fake_aggregate_clonotypes)
        monkeypatch.setattr("tcrsift.filter.filter_clonotypes", fake_filter_clonotypes)
        monkeypatch.setattr("tcrsift.filter.split_by_tier", fake_split_by_tier)
        monkeypatch.setattr("tcrsift.til.load_til_samples", lambda *_a, **_k: {})

        output_dir = tmp_path / "out"

        args = argparse.Namespace(
            sample_sheet=str(sample_sheet_path),
            output_dir=str(output_dir),
            config=None,
            generate_plots=False,
            generate_report=False,
            no_leaders=True,
            include_constant=False,
            single_chain=False,
            til_samples="Patient1_TIL",
            min_genes=None,
            max_genes=None,
            min_counts=None,
            max_counts=None,
            min_mito_pct=None,
            max_mito_pct=None,
            cd4_cd8_ratio=None,
            min_cd3_reads=None,
            group_by=None,
            handle_doublets=None,
            min_umi=None,
            tcell_type=None,
            method=None,
            min_cells=None,
            min_frequency=None,
            require_complete=None,
            fdr_tiers=None,
            vdjdb_path=None,
            iedb_path=None,
            cedar_path=None,
            match_by=None,
            exclude_viral=None,
            flag_only=None,
            til_match_by=None,
            min_til_cells=None,
            alpha_leader=None,
            beta_leader=None,
            leaders_from_contigs=False,
            contigs_dir=None,
            linker=None,
            constant_source=None,
            skip_plots=None,
            verbose=False,
        )

        cmd_run(args)

    def test_run_with_direct_til_sample_specs(self, tmp_path, monkeypatch):
        """cmd_run should support repeatable --til-sample specs without TIL in sample sheet."""
        import anndata as ad
        import pandas as pd

        sample_sheet_path = tmp_path / "samples.yaml"
        sample_sheet_path.write_text(
            """
samples:
  - sample: "Culture1"
    vdj_dir: "/data/culture/vdj"
    source: "culture"
"""
        )

        # Create a valid direct TIL CSV path for early validation
        til_csv = tmp_path / "til.csv"
        til_csv.write_text("CDR3_alpha,CDR3_beta\nA,B\n")

        obs = pd.DataFrame(
            {
                "sample": ["Culture1"],
                "source": ["culture"],
            },
            index=["cell1"],
        )
        adata = ad.AnnData(obs=obs)

        def fake_load_samples(*_args, **_kwargs):
            return adata

        def fake_phenotype_cells(adata_in, *_args, **_kwargs):
            return adata_in

        def fake_filter_by_tcell_type(adata_in, *_args, **_kwargs):
            return adata_in

        def fake_aggregate_clonotypes(*_args, **_kwargs):
            return pd.DataFrame(
                {
                    "CDR3ab": ["A_B"],
                    "CDR3_alpha": ["A"],
                    "CDR3_beta": ["B"],
                    "cell_count": [1],
                }
            )

        def fake_filter_clonotypes(df, *_args, **_kwargs):
            return df.assign(tier="tier1")

        def fake_split_by_tier(df, *_args, **_kwargs):
            return {"tier1": df}

        def fake_load_til_specs(specs):
            assert specs == [f"T1=csv:{til_csv}"]
            return {"T1": pd.DataFrame({"CDR3_alpha": ["A"], "CDR3_beta": ["B"], "sample": ["T1"]})}

        def fake_match_til(df, *_args, **_kwargs):
            return df.assign(
                til_match=False,
                til_samples="",
                til_cell_count=0,
                til_frequency=0.0,
            )

        monkeypatch.setattr("tcrsift.loader.load_samples", fake_load_samples)
        monkeypatch.setattr("tcrsift.phenotype.phenotype_cells", fake_phenotype_cells)
        monkeypatch.setattr("tcrsift.phenotype.filter_by_tcell_type", fake_filter_by_tcell_type)
        monkeypatch.setattr("tcrsift.clonotype.aggregate_clonotypes", fake_aggregate_clonotypes)
        monkeypatch.setattr("tcrsift.filter.filter_clonotypes", fake_filter_clonotypes)
        monkeypatch.setattr("tcrsift.filter.split_by_tier", fake_split_by_tier)
        monkeypatch.setattr("tcrsift.til.load_til_specs", fake_load_til_specs)
        monkeypatch.setattr("tcrsift.til.match_til", fake_match_til)

        output_dir = tmp_path / "out"

        args = argparse.Namespace(
            sample_sheet=str(sample_sheet_path),
            output_dir=str(output_dir),
            config=None,
            generate_plots=False,
            generate_report=False,
            no_leaders=True,
            include_constant=False,
            single_chain=False,
            til_samples=None,
            til_sample=[f"T1=csv:{til_csv}"],
            min_genes=None,
            max_genes=None,
            min_counts=None,
            max_counts=None,
            min_mito_pct=None,
            max_mito_pct=None,
            cd4_cd8_ratio=None,
            min_cd3_reads=None,
            group_by=None,
            handle_doublets=None,
            min_umi=None,
            tcell_type=None,
            method=None,
            min_cells=None,
            min_frequency=None,
            require_complete=None,
            fdr_tiers=None,
            vdjdb_path=None,
            iedb_path=None,
            cedar_path=None,
            match_by=None,
            exclude_viral=None,
            flag_only=None,
            til_match_by=None,
            min_til_cells=None,
            alpha_leader=None,
            beta_leader=None,
            leaders_from_contigs=False,
            contigs_dir=None,
            linker=None,
            constant_source=None,
            skip_plots=None,
            verbose=False,
        )

        cmd_run(args)


class TestRunEmitsPerSampleTier:
    """`tcrsift run` enriches the emitted clone_sample_long.csv with a
    per-sample abundance tier via tcrsift.selection (#125 increment)."""

    def test_clone_sample_long_has_per_sample_tier_column(self, tmp_path, monkeypatch):
        import anndata as ad
        import numpy as np

        sample_sheet_path = tmp_path / "samples.yaml"
        sample_sheet_path.write_text(
            """
samples:
  - sample: "S1"
    vdj_dir: "/data/s1/vdj"
    source: "culture"
  - sample: "S2"
    vdj_dir: "/data/s2/vdj"
    source: "culture"
"""
        )

        # Two samples, two clones with real CDR3 columns so the live
        # build_clone_sample_long produces per-(clone, sample) rows.
        rows = (
            [("S1", "CAVA", "CASS_A")] * 12
            + [("S1", "CAVB", "CASS_B")] * 3
            + [("S2", "CAVA", "CASS_A")] * 4
        )
        n = len(rows)
        obs = pd.DataFrame(
            {
                "sample": [r[0] for r in rows],
                "source": ["culture"] * n,
                "CDR3_alpha": [r[1] for r in rows],
                "CDR3_beta": [r[2] for r in rows],
            },
            index=[f"cell{i}" for i in range(n)],
        )
        adata = ad.AnnData(X=np.zeros((n, 1), dtype=np.float32), obs=obs)

        monkeypatch.setattr("tcrsift.loader.load_samples", lambda *a, **k: adata)
        monkeypatch.setattr(
            "tcrsift.phenotype.phenotype_cells", lambda a, *args, **k: a
        )
        monkeypatch.setattr(
            "tcrsift.phenotype.filter_by_tcell_type", lambda a, *args, **k: a
        )
        monkeypatch.setattr(
            "tcrsift.clonotype.aggregate_clonotypes",
            lambda *a, **k: pd.DataFrame(
                {
                    "CDR3ab": ["CAVA_CASS_A", "CAVB_CASS_B"],
                    "CDR3_alpha": ["CAVA", "CAVB"],
                    "CDR3_beta": ["CASS_A", "CASS_B"],
                    "cell_count": [16, 3],
                }
            ),
        )
        monkeypatch.setattr(
            "tcrsift.filter.filter_clonotypes", lambda df, *a, **k: df.assign(tier="tier1")
        )
        monkeypatch.setattr(
            "tcrsift.filter.split_by_tier", lambda df, *a, **k: {"tier1": df}
        )
        monkeypatch.setattr("tcrsift.til.load_til_samples", lambda *a, **k: {})

        output_dir = tmp_path / "out"
        args = argparse.Namespace(
            sample_sheet=str(sample_sheet_path),
            output_dir=str(output_dir),
            config=None,
            generate_plots=False,
            generate_report=False,
            no_leaders=True,
            include_constant=False,
            single_chain=False,
            til_samples=None,
            min_genes=None, max_genes=None, min_counts=None, max_counts=None,
            min_mito_pct=None, max_mito_pct=None, cd4_cd8_ratio=None,
            min_cd3_reads=None, group_by=None, handle_doublets=None, min_umi=None,
            tcell_type=None, method=None, min_cells=None, min_frequency=None,
            require_complete=None, fdr_tiers=None, vdjdb_path=None, iedb_path=None,
            cedar_path=None, match_by=None, exclude_viral=None, flag_only=None,
            til_match_by=None, min_til_cells=None, alpha_leader=None,
            beta_leader=None, leaders_from_contigs=False, contigs_dir=None,
            linker=None, constant_source=None, skip_plots=None, verbose=False,
        )

        cmd_run(args)

        long_csv = output_dir / "data" / "clone_sample_long.csv"
        assert long_csv.exists(), "run should emit clone_sample_long.csv for >=2 samples"
        long_df = pd.read_csv(long_csv)
        assert "per_sample_tier" in long_df.columns
        # CAVA in S1: 12 cells / 15 in-sample = 0.8 freq -> tier1.
        a_s1 = long_df[
            (long_df["CDR3ab"] == "CAVA_CASS_A") & (long_df["sample"] == "S1")
        ]
        assert a_s1["per_sample_tier"].iloc[0] == "tier1"


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
        ss = SampleSheet(
            samples=[
                Sample(sample="C1", vdj_dir="/path1", source="culture"),
                Sample(sample="T1", vdj_dir="/path2", source="til"),
                Sample(sample="Tet1", vdj_dir="/path3", source="tetramer"),
                Sample(sample="SCT1", vdj_dir="/path4", source="sct"),
            ]
        )
        tetramer = ss.get_tetramer_samples()

        assert len(tetramer) == 2
        assert all(s.is_tetramer_or_sct() for s in tetramer)


# =============================================================================
# CLI Validation Tests
# =============================================================================


class TestCliConditionalRequirement:
    """Tests for validate_cli_conditional_requirement function."""

    def test_no_error_when_condition_not_met(self):
        """Should not raise when condition is not triggered."""
        args = argparse.Namespace(
            leaders_from_contigs=False,
            contigs_dir=None,
        )
        # Should not raise
        validate_cli_conditional_requirement(
            args,
            required_arg="contigs_dir",
            condition_args=["leaders_from_contigs"],
            condition_description="when using --leaders-from-contigs",
        )

    def test_no_error_when_requirement_satisfied(self):
        """Should not raise when condition is met and requirement satisfied."""
        args = argparse.Namespace(
            leaders_from_contigs=True,
            contigs_dir="/path/to/contigs",
        )
        # Should not raise
        validate_cli_conditional_requirement(
            args,
            required_arg="contigs_dir",
            condition_args=["leaders_from_contigs"],
            condition_description="when using --leaders-from-contigs",
        )

    def test_error_when_requirement_missing(self):
        """Should raise when condition met but requirement missing."""
        args = argparse.Namespace(
            leaders_from_contigs=True,
            contigs_dir=None,
        )
        with pytest.raises(TCRsiftValidationError) as exc_info:
            validate_cli_conditional_requirement(
                args,
                required_arg="contigs_dir",
                condition_args=["leaders_from_contigs"],
                condition_description="when using --leaders-from-contigs",
            )
        assert "--contigs-dir is required" in str(exc_info.value)
        assert "--leaders-from-contigs" in str(exc_info.value)

    def test_condition_values_matching(self):
        """Should check specific values when condition_values provided."""
        args = argparse.Namespace(
            alpha_leader="from_contig",
            beta_leader="CD8A",
            contigs_dir=None,
        )
        with pytest.raises(TCRsiftValidationError) as exc_info:
            validate_cli_conditional_requirement(
                args,
                required_arg="contigs_dir",
                condition_args=["alpha_leader", "beta_leader"],
                condition_values=["from_contig"],
                condition_description="when using leader=from_contig",
            )
        assert "--contigs-dir is required" in str(exc_info.value)
        assert "--alpha-leader=from_contig" in str(exc_info.value)

    def test_condition_values_not_matching(self):
        """Should not raise when condition_values don't match."""
        args = argparse.Namespace(
            alpha_leader="CD28",
            beta_leader="CD8A",
            contigs_dir=None,
        )
        # Should not raise - neither leader is "from_contig"
        validate_cli_conditional_requirement(
            args,
            required_arg="contigs_dir",
            condition_args=["alpha_leader", "beta_leader"],
            condition_values=["from_contig"],
            condition_description="when using leader=from_contig",
        )

    def test_multiple_conditions_any_triggers(self):
        """Should trigger when any condition is met."""
        args = argparse.Namespace(
            alpha_leader="CD28",
            beta_leader="from_contig",  # This triggers
            contigs_dir=None,
        )
        with pytest.raises(TCRsiftValidationError) as exc_info:
            validate_cli_conditional_requirement(
                args,
                required_arg="contigs_dir",
                condition_args=["alpha_leader", "beta_leader"],
                condition_values=["from_contig"],
                condition_description="when using leader=from_contig",
            )
        assert "--beta-leader=from_contig" in str(exc_info.value)

    def test_empty_string_treated_as_missing(self):
        """Empty string should be treated as missing."""
        args = argparse.Namespace(
            leaders_from_contigs=True,
            contigs_dir="",
        )
        with pytest.raises(TCRsiftValidationError):
            validate_cli_conditional_requirement(
                args,
                required_arg="contigs_dir",
                condition_args=["leaders_from_contigs"],
                condition_description="when using --leaders-from-contigs",
            )


class TestCliMutuallyExclusive:
    """Tests for validate_cli_mutually_exclusive function."""

    def test_no_error_when_none_set(self):
        """Should not raise when no args in group are set."""
        args = argparse.Namespace(
            no_leaders=False,
            leaders_from_contigs=False,
        )
        validate_cli_mutually_exclusive(
            args,
            arg_names=["no_leaders", "leaders_from_contigs"],
            group_description="leader options",
        )

    def test_no_error_when_one_set(self):
        """Should not raise when only one arg is set."""
        args = argparse.Namespace(
            no_leaders=True,
            leaders_from_contigs=False,
        )
        validate_cli_mutually_exclusive(
            args,
            arg_names=["no_leaders", "leaders_from_contigs"],
            group_description="leader options",
        )

    def test_error_when_multiple_set(self):
        """Should raise when multiple mutually exclusive args set."""
        args = argparse.Namespace(
            no_leaders=True,
            leaders_from_contigs=True,
        )
        with pytest.raises(TCRsiftValidationError) as exc_info:
            validate_cli_mutually_exclusive(
                args,
                arg_names=["no_leaders", "leaders_from_contigs"],
                group_description="leader options",
            )
        assert "Cannot use" in str(exc_info.value)
        assert "--no-leaders" in str(exc_info.value)
        assert "--leaders-from-contigs" in str(exc_info.value)


class TestValidateAssembleArgs:
    """Tests for validate_assemble_args function."""

    def test_valid_args_default_leaders(self):
        """Should pass with default leader options."""
        args = argparse.Namespace(
            input="clonotypes.csv",
            output="assembled.csv",
            alpha_leader="CD28",
            beta_leader="CD8A",
            leaders_from_contigs=False,
            contigs_dir=None,
        )
        # Should not raise
        validate_assemble_args(args)

    def test_valid_args_with_contigs(self):
        """Should pass when from_contig used with contigs_dir."""
        args = argparse.Namespace(
            input="clonotypes.csv",
            output="assembled.csv",
            alpha_leader="from_contig",
            beta_leader="from_contig",
            leaders_from_contigs=False,
            contigs_dir="/path/to/contigs",
        )
        # Mock directory existence
        with patch("tcrsift.validation.validate_directory_exists"):
            validate_assemble_args(args)

    def test_error_leaders_from_contigs_without_dir(self):
        """Should fail when --leaders-from-contigs used without --contigs-dir."""
        args = argparse.Namespace(
            input="clonotypes.csv",
            output="assembled.csv",
            alpha_leader="CD28",
            beta_leader="CD8A",
            leaders_from_contigs=True,
            contigs_dir=None,
        )
        with pytest.raises(TCRsiftValidationError) as exc_info:
            validate_assemble_args(args)
        assert "--contigs-dir is required" in str(exc_info.value)
        assert "--leaders-from-contigs" in str(exc_info.value)

    def test_error_alpha_from_contig_without_dir(self):
        """Should fail when alpha_leader=from_contig without contigs_dir."""
        args = argparse.Namespace(
            input="clonotypes.csv",
            output="assembled.csv",
            alpha_leader="from_contig",
            beta_leader="CD8A",
            leaders_from_contigs=False,
            contigs_dir=None,
        )
        with pytest.raises(TCRsiftValidationError) as exc_info:
            validate_assemble_args(args)
        assert "--contigs-dir is required" in str(exc_info.value)

    def test_error_beta_from_contig_without_dir(self):
        """Should fail when beta_leader=from_contig without contigs_dir."""
        args = argparse.Namespace(
            input="clonotypes.csv",
            output="assembled.csv",
            alpha_leader="CD28",
            beta_leader="from_contig",
            leaders_from_contigs=False,
            contigs_dir=None,
        )
        with pytest.raises(TCRsiftValidationError) as exc_info:
            validate_assemble_args(args)
        assert "--contigs-dir is required" in str(exc_info.value)


class TestValidateRunArgs:
    """Tests for validate_run_args function."""

    @pytest.fixture
    def sample_sheet_file(self, tmp_path):
        """Create a sample sheet file for testing."""
        ss = tmp_path / "samples.yaml"
        ss.write_text("samples:\n  - sample: S1\n    vdj_dir: /path\n")
        return ss

    def test_valid_args_default(self, sample_sheet_file):
        """Should pass with default options."""
        args = argparse.Namespace(
            sample_sheet=str(sample_sheet_file),
            output_dir="/output",
            alpha_leader=None,
            beta_leader=None,
            leaders_from_contigs=False,
            contigs_dir=None,
        )
        validate_run_args(args)

    def test_error_leaders_from_contigs_without_dir(self, sample_sheet_file):
        """Should fail when --leaders-from-contigs used without --contigs-dir."""
        args = argparse.Namespace(
            sample_sheet=str(sample_sheet_file),
            output_dir="/output",
            alpha_leader=None,
            beta_leader=None,
            leaders_from_contigs=True,
            contigs_dir=None,
        )
        with pytest.raises(TCRsiftValidationError) as exc_info:
            validate_run_args(args)
        assert "--contigs-dir is required" in str(exc_info.value)

    def test_error_missing_sample_sheet(self, tmp_path):
        """Should fail when sample sheet doesn't exist."""
        args = argparse.Namespace(
            sample_sheet=str(tmp_path / "nonexistent.yaml"),
            output_dir="/output",
            alpha_leader=None,
            beta_leader=None,
            leaders_from_contigs=False,
            contigs_dir=None,
        )
        with pytest.raises(TCRsiftValidationError) as exc_info:
            validate_run_args(args)
        assert "does not exist" in str(exc_info.value)

    def test_error_conflicting_til_sources(self, sample_sheet_file):
        """Should fail when both --til-samples and --til-sample are provided."""
        args = argparse.Namespace(
            sample_sheet=str(sample_sheet_file),
            output_dir="/output",
            alpha_leader=None,
            beta_leader=None,
            leaders_from_contigs=False,
            contigs_dir=None,
            til_samples="TIL1",
            til_sample=["TIL1=csv:/tmp/til.csv"],
        )
        with pytest.raises(TCRsiftValidationError) as exc_info:
            validate_run_args(args)
        assert "Conflicting TIL options" in str(exc_info.value)


class TestValidateAnnotateGexArgs:
    """Tests for validate_annotate_gex_args function."""

    @pytest.fixture
    def gex_file(self, tmp_path):
        """Create a fake GEX file for testing."""
        gex = tmp_path / "matrix.h5"
        gex.write_bytes(b"fake h5 content")
        return gex

    def test_valid_args(self, gex_file):
        """Should pass with valid gex_file."""
        args = argparse.Namespace(
            gex_file=str(gex_file),
        )
        validate_annotate_gex_args(args)

    def test_error_missing_gex_file(self, tmp_path):
        """Should fail when gex_file doesn't exist."""
        args = argparse.Namespace(
            gex_file=str(tmp_path / "nonexistent.h5"),
        )
        with pytest.raises(TCRsiftValidationError) as exc_info:
            validate_annotate_gex_args(args)
        assert "does not exist" in str(exc_info.value)


class TestAssembleParserHelpStrings:
    """Tests for assemble command help string organization."""

    def test_assemble_has_required_group(self):
        """Assemble parser should have a 'required arguments' group."""
        parser = create_parser()

        # Check the parser has expected structure
        subparsers = parser._subparsers._group_actions[0].choices
        asm_parser = subparsers["assemble"]

        group_titles = [g.title for g in asm_parser._action_groups]
        assert "required arguments" in group_titles
        assert "leader peptide options" in group_titles
        assert "sequence options" in group_titles

    def test_assemble_description_shows_conditional_requirement(self):
        """Assemble parser description should mention conditional requirements."""
        parser = create_parser()
        subparsers = parser._subparsers._group_actions[0].choices
        asm_parser = subparsers["assemble"]

        assert "CONDITIONALLY REQUIRED" in asm_parser.description
        assert "--contigs-dir" in asm_parser.description

    def test_match_til_has_required_group(self):
        """Match-til parser should have a 'required arguments' group."""
        parser = create_parser()
        subparsers = parser._subparsers._group_actions[0].choices
        til_parser = subparsers["match-til"]

        group_titles = [g.title for g in til_parser._action_groups]
        assert "required arguments" in group_titles

    def test_annotate_gex_has_required_group(self):
        """Annotate-gex parser should have a 'required arguments' group."""
        parser = create_parser()
        subparsers = parser._subparsers._group_actions[0].choices
        gex_parser = subparsers["annotate-gex"]

        group_titles = [g.title for g in gex_parser._action_groups]
        assert "required arguments" in group_titles


class TestCliEarlyValidation:
    """Integration tests for early CLI validation."""

    def test_assemble_fails_early_with_missing_contigs_dir(self, tmp_path, capsys):
        """Assemble should fail before loading data when contigs_dir missing."""
        # Create a valid input file
        input_csv = tmp_path / "clonotypes.csv"
        pd.DataFrame(
            {
                "CDR3_alpha": ["CAVX"],
                "CDR3_beta": ["CASSX"],
            }
        ).to_csv(input_csv, index=False)

        parser = create_parser()
        args = parser.parse_args(
            [
                "assemble",
                "-i",
                str(input_csv),
                "-o",
                str(tmp_path / "out.csv"),
                "--leaders-from-contigs",  # This requires --contigs-dir
            ]
        )

        # Should fail at validation, not at data loading
        with pytest.raises(TCRsiftValidationError) as exc_info:
            args.func(args)

        assert "--contigs-dir is required" in str(exc_info.value)

    def test_assemble_validates_contigs_dir_exists(self, tmp_path):
        """Assemble should validate that contigs_dir exists."""
        input_csv = tmp_path / "clonotypes.csv"
        pd.DataFrame(
            {
                "CDR3_alpha": ["CAVX"],
                "CDR3_beta": ["CASSX"],
            }
        ).to_csv(input_csv, index=False)

        parser = create_parser()
        args = parser.parse_args(
            [
                "assemble",
                "-i",
                str(input_csv),
                "-o",
                str(tmp_path / "out.csv"),
                "--leaders-from-contigs",
                "--contigs-dir",
                str(tmp_path / "nonexistent_dir"),
            ]
        )

        with pytest.raises(TCRsiftValidationError) as exc_info:
            args.func(args)

        assert "does not exist" in str(exc_info.value)


class TestMatchTilParser:
    """Tests for match-til CLI parser configuration."""

    def test_parser_has_match_til_command(self):
        """Parser should include match-til command."""
        parser = create_parser()
        # Parse valid match-til command to verify it exists
        args = parser.parse_args(
            [
                "match-til",
                "-i",
                "clonotypes.csv",
                "-o",
                "matched.csv",
                "--til-csv",
                "til.csv",
            ]
        )
        assert args.command == "match-til"
        assert hasattr(args, "func")

    def test_match_til_required_args(self):
        """Test that required arguments are enforced."""
        parser = create_parser()

        # Should fail without -i and -o
        with pytest.raises(SystemExit):
            parser.parse_args(["match-til"])

    def test_match_til_all_til_source_options(self):
        """Test that all TIL source options parse correctly."""
        parser = create_parser()

        # Test --til-h5ad
        args1 = parser.parse_args(
            [
                "match-til",
                "-i",
                "clonotypes.csv",
                "-o",
                "matched.csv",
                "--til-h5ad",
                "til.h5ad",
            ]
        )
        assert args1.til_h5ad == "til.h5ad"
        assert args1.til_csv is None

        # Test --til-csv
        args2 = parser.parse_args(
            [
                "match-til",
                "-i",
                "clonotypes.csv",
                "-o",
                "matched.csv",
                "--til-csv",
                "til.csv",
            ]
        )
        assert args2.til_csv == "til.csv"
        assert args2.til_h5ad is None

        # Test --til-vdj-dir
        args3 = parser.parse_args(
            [
                "match-til",
                "-i",
                "clonotypes.csv",
                "-o",
                "matched.csv",
                "--til-vdj-dir",
                "/path/to/vdj",
            ]
        )
        assert args3.til_vdj_dir == "/path/to/vdj"

        # Test --sample-sheet
        args4 = parser.parse_args(
            [
                "match-til",
                "-i",
                "clonotypes.csv",
                "-o",
                "matched.csv",
                "-s",
                "samples.yaml",
            ]
        )
        assert args4.sample_sheet == "samples.yaml"

        # Test repeatable --til-sample specs
        args5 = parser.parse_args(
            [
                "match-til",
                "-i",
                "clonotypes.csv",
                "-o",
                "matched.csv",
                "--til-sample",
                "T1=csv:/path/to/til_t1.csv",
                "--til-sample",
                "T2=vdj:/path/to/til_t2_vdj_outs",
            ]
        )
        assert args5.til_sample == [
            "T1=csv:/path/to/til_t1.csv",
            "T2=vdj:/path/to/til_t2_vdj_outs",
        ]

    def test_match_til_matching_options(self):
        """Test matching options parse correctly."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "match-til",
                "-i",
                "clonotypes.csv",
                "-o",
                "matched.csv",
                "--til-csv",
                "til.csv",
                "--match-by",
                "CDR3b_only",
                "--min-til-cells",
                "5",
            ]
        )
        assert args.match_by == "CDR3b_only"
        assert args.min_til_cells == 5


class TestTilClonotypeParser:
    """Tests for til-clonotype CLI parser configuration."""

    def test_parser_has_til_clonotype_command(self):
        """Parser should include til-clonotype command."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "til-clonotype",
                "-o",
                "til_clonotypes.csv",
                "--til-csv",
                "til.csv",
            ]
        )
        assert args.command == "til-clonotype"
        assert hasattr(args, "func")

    def test_til_clonotype_accepts_repeatable_til_sample(self):
        """til-clonotype should parse repeatable --til-sample options."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "til-clonotype",
                "-o",
                "til_clonotypes.csv",
                "--til-sample",
                "T1=csv:/path/to/til_t1.csv",
                "--til-sample",
                "T2=h5ad:/path/to/til_t2.h5ad",
            ]
        )
        assert args.til_sample == [
            "T1=csv:/path/to/til_t1.csv",
            "T2=h5ad:/path/to/til_t2.h5ad",
        ]


class TestTilSelectParser:
    """Tests for til-select CLI parser configuration."""

    def test_parser_has_til_select_command(self):
        """Parser should include til-select command."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "til-select",
                "--data-dir",
                "data",
            ]
        )
        assert args.command == "til-select"
        assert hasattr(args, "func")

    def test_til_select_accepts_samples_alias(self):
        """til-select should accept v2-style --samples mappings."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "til-select",
                "--samples",
                "T1=consensus_annotations.T1.csv,clonotypes.T1.csv",
                "T2=consensus_annotations.T2.csv,clonotypes.T2.csv",
            ]
        )
        assert args.samples == [
            "T1=consensus_annotations.T1.csv,clonotypes.T1.csv",
            "T2=consensus_annotations.T2.csv,clonotypes.T2.csv",
        ]

    def test_til_select_accepts_inputs_alias(self):
        """til-select should treat --inputs as alias for --samples."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "til-select",
                "--inputs",
                "T1=consensus_annotations.T1.csv,clonotypes.T1.csv",
            ]
        )
        assert args.samples == [
            "T1=consensus_annotations.T1.csv,clonotypes.T1.csv",
        ]


class TestValidateMatchTilArgs:
    """Tests for match-til argument validation."""

    def test_no_til_source_raises(self, tmp_path):
        """Should raise when no TIL source is provided."""
        from tcrsift.validation import validate_match_til_args

        # Create an args-like object with no TIL sources
        args = argparse.Namespace(
            sample_sheet=None,
            til_h5ad=None,
            til_csv=None,
            til_vdj_dir=None,
            til_sample=None,
        )

        with pytest.raises(TCRsiftValidationError, match="No TIL data source specified"):
            validate_match_til_args(args)

    def test_multiple_til_sources_raises(self, tmp_path):
        """Should raise when multiple TIL sources are provided."""
        from tcrsift.validation import validate_match_til_args

        # Create temp files so they "exist" for validation
        til_csv = tmp_path / "til.csv"
        til_csv.write_text("CDR3_alpha,CDR3_beta\nCAV,CAS\n")

        til_h5ad = tmp_path / "til.h5ad"
        til_h5ad.write_text("dummy")  # Not a real h5ad but just for path existence

        args = argparse.Namespace(
            sample_sheet=None,
            til_h5ad=str(til_h5ad),
            til_csv=str(til_csv),
            til_vdj_dir=None,
            til_sample=None,
        )

        with pytest.raises(TCRsiftValidationError, match="Multiple TIL data sources"):
            validate_match_til_args(args)

    def test_single_til_csv_valid(self, tmp_path):
        """Single til_csv source should be valid."""
        from tcrsift.validation import validate_match_til_args

        til_csv = tmp_path / "til.csv"
        til_csv.write_text("CDR3_alpha,CDR3_beta\nCAV,CAS\n")

        args = argparse.Namespace(
            sample_sheet=None,
            til_h5ad=None,
            til_csv=str(til_csv),
            til_vdj_dir=None,
            til_sample=None,
        )

        # Should not raise
        validate_match_til_args(args)

    def test_nonexistent_til_csv_raises(self, tmp_path):
        """Should raise when til_csv file doesn't exist."""
        from tcrsift.validation import validate_match_til_args

        args = argparse.Namespace(
            sample_sheet=None,
            til_h5ad=None,
            til_csv=str(tmp_path / "nonexistent.csv"),
            til_vdj_dir=None,
            til_sample=None,
        )

        with pytest.raises(TCRsiftValidationError, match="does not exist"):
            validate_match_til_args(args)

    def test_nonexistent_til_vdj_dir_raises(self, tmp_path):
        """Should raise when til_vdj_dir doesn't exist."""
        from tcrsift.validation import validate_match_til_args

        args = argparse.Namespace(
            sample_sheet=None,
            til_h5ad=None,
            til_csv=None,
            til_vdj_dir=str(tmp_path / "nonexistent_dir"),
            til_sample=None,
        )

        with pytest.raises(TCRsiftValidationError, match="does not exist"):
            validate_match_til_args(args)

    def test_til_sample_specs_valid(self, tmp_path):
        """Repeatable --til-sample specs should validate."""
        from tcrsift.validation import validate_match_til_args

        til_csv = tmp_path / "til.csv"
        til_csv.write_text("CDR3_alpha,CDR3_beta\nCAV,CAS\n")
        vdj_dir = tmp_path / "vdj_outs"
        vdj_dir.mkdir()

        args = argparse.Namespace(
            sample_sheet=None,
            til_h5ad=None,
            til_csv=None,
            til_vdj_dir=None,
            til_sample=[f"T1=csv:{til_csv}", f"T2=vdj:{vdj_dir}"],
        )

        validate_match_til_args(args)

    def test_til_sample_specs_invalid_format_raises(self, tmp_path):
        """Invalid --til-sample format should raise a clear error."""
        from tcrsift.validation import validate_match_til_args

        args = argparse.Namespace(
            sample_sheet=None,
            til_h5ad=None,
            til_csv=None,
            til_vdj_dir=None,
            til_sample=["not_a_valid_spec"],
        )

        with pytest.raises(TCRsiftValidationError, match="Invalid --til-sample spec"):
            validate_match_til_args(args)

    def test_til_sample_specs_duplicate_names_raise(self, tmp_path):
        """Duplicate sample names in --til-sample should raise."""
        from tcrsift.validation import validate_match_til_args

        til1 = tmp_path / "til1.csv"
        til2 = tmp_path / "til2.csv"
        til1.write_text("CDR3_alpha,CDR3_beta\nCAV,CAS\n")
        til2.write_text("CDR3_alpha,CDR3_beta\nCAV,CAS\n")

        args = argparse.Namespace(
            sample_sheet=None,
            til_h5ad=None,
            til_csv=None,
            til_vdj_dir=None,
            til_sample=[f"T1=csv:{til1}", f"T1=csv:{til2}"],
        )

        with pytest.raises(TCRsiftValidationError, match="Duplicate TIL sample name"):
            validate_match_til_args(args)

    def test_cmd_match_til_rejects_multiple_source_modes(self, tmp_path):
        """match-til command should fail early when multiple source modes are provided."""
        til_csv = tmp_path / "til.csv"
        til_csv.write_text("CDR3_alpha,CDR3_beta\nCAV,CAS\n")

        parser = create_parser()
        args = parser.parse_args(
            [
                "match-til",
                "-i",
                "clonotypes.csv",
                "-o",
                "matched.csv",
                "--til-csv",
                str(til_csv),
                "--til-sample",
                f"T1=csv:{til_csv}",
            ]
        )

        with pytest.raises(TCRsiftValidationError, match="Multiple TIL data sources"):
            args.func(args)
