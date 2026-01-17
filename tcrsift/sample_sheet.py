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
"""
Sample sheet parsing for TCRsift.

Supports both CSV and YAML formats for specifying samples and their metadata.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# Antigen types and their expected T cell responses
ANTIGEN_TYPE_TCELL_EXPECTATIONS = {
    "short_peptide": "CD8",      # 8-11aa, MHC-I direct binding
    "long_peptide": "mixed",     # 15-25+aa, requires processing, favors CD4
    "peptide_pool": "mixed",     # pool of peptides, mixed responses
    "minigene": "mixed",         # requires processing, favors CD4
    "minigene_library": "mixed", # library of minigenes, mixed responses
    "whole_protein": "mixed",    # requires processing, favors CD4
    "mrna": "mixed",             # mRNA encoding antigen(s), requires processing
    "tetramer_mhc1": "CD8",      # MHC-I tetramer selection (single antigen)
    "tetramer_mhc2": "CD4",      # MHC-II tetramer selection (single antigen)
    "sct": "CD8",                # single-chain trimer (pMHC-I alpha-B2M-peptide fusion)
}

VALID_ANTIGEN_TYPES = set(ANTIGEN_TYPE_TCELL_EXPECTATIONS.keys())
VALID_SOURCES = {"culture", "tetramer", "sct", "til"}
VALID_TCELL_TYPES = {"CD4", "CD8", "mixed", None}
VALID_MHC_BLOCKING = {"MHC-I", "MHC-II", None}
VALID_PRE_SORTED = {"CD4", "CD8", None}


@dataclass
class Sample:
    """Represents a single sample with its metadata."""
    sample: str
    gex_dir: str | None = None
    vdj_dir: str | None = None
    antigen_type: str | None = None
    antigen_description: str | None = None
    culture_days: int | None = None
    tcell_type_expected: str | None = None
    pre_sorted: str | None = None
    mhc_blocking: str | None = None
    source: str = "culture"
    # Antigen fields - what was given to APCs (protein, long peptide, minigene, etc.)
    antigen_name: str | None = None  # name of source antigen (e.g., "PRAME", "CMV pp65")
    antigen_sequence: str | None = None  # sequence of source antigen (may be long)
    # Epitope fields - minimal peptide that binds MHC
    epitope_sequence: str | None = None  # minimal peptide AA sequence (e.g., "SLLQHLIGL")
    mhc_allele: str | None = None  # MHC restriction (e.g., "HLA-A*02:01")
    # Pool/library fields - lists for multiple antigens
    antigen_names: list | None = None  # list of source antigen names
    antigen_sequences: list | None = None  # list of source antigen sequences
    epitope_sequences: list | None = None  # list of minimal epitope sequences (if known)
    # Other metadata
    tissue: str | None = None
    patient_id: str | None = None
    # Experiment grouping for multi-source unification
    experiment: str | None = None  # experiment name for unification (e.g., "TIL", "Culture", "SCT")
    # SCT-specific fields
    sct_path: str | None = None  # path to SCT Excel file
    sct_sheet: str = "Cell"  # sheet name in SCT Excel file
    # Standalone GEX file (for augmentation without full CellRanger output)
    gex_path: str | None = None  # path to 10x filtered_feature_bc_matrix.h5 file

    def __post_init__(self):
        # Validate at least one data source
        has_cellranger = self.gex_dir or self.vdj_dir
        has_sct = self.sct_path is not None
        if not has_cellranger and not has_sct:
            raise ValueError(
                f"Sample '{self.sample}' must have at least gex_dir, vdj_dir, or sct_path"
            )

        # Validate antigen type
        if self.antigen_type and self.antigen_type not in VALID_ANTIGEN_TYPES:
            raise ValueError(
                f"Invalid antigen_type '{self.antigen_type}' for sample '{self.sample}'. "
                f"Valid types: {VALID_ANTIGEN_TYPES}"
            )

        # Validate source
        if self.source not in VALID_SOURCES:
            raise ValueError(
                f"Invalid source '{self.source}' for sample '{self.sample}'. "
                f"Valid sources: {VALID_SOURCES}"
            )

        # Validate tcell_type_expected
        if self.tcell_type_expected and self.tcell_type_expected not in VALID_TCELL_TYPES:
            raise ValueError(
                f"Invalid tcell_type_expected '{self.tcell_type_expected}' for sample '{self.sample}'. "
                f"Valid types: {VALID_TCELL_TYPES}"
            )

        # Validate pre_sorted
        if self.pre_sorted and self.pre_sorted not in VALID_PRE_SORTED:
            raise ValueError(
                f"Invalid pre_sorted '{self.pre_sorted}' for sample '{self.sample}'. "
                f"Valid values: {VALID_PRE_SORTED}"
            )

        # Validate mhc_blocking
        if self.mhc_blocking and self.mhc_blocking not in VALID_MHC_BLOCKING:
            raise ValueError(
                f"Invalid mhc_blocking '{self.mhc_blocking}' for sample '{self.sample}'. "
                f"Valid values: {VALID_MHC_BLOCKING}"
            )

    def get_expected_tcell_type(self) -> str | None:
        """
        Determine expected T cell type based on antigen type, blocking, and sorting.

        Returns the most specific expectation available.
        """
        # Direct specification takes priority
        if self.tcell_type_expected:
            return self.tcell_type_expected

        # Pre-sorting is definitive
        if self.pre_sorted:
            return self.pre_sorted

        # MHC blocking tells us what's NOT expected
        if self.mhc_blocking == "MHC-I":
            return "CD4"  # CD8 responses blocked
        if self.mhc_blocking == "MHC-II":
            return "CD8"  # CD4 responses blocked

        # Antigen type gives us expectations
        if self.antigen_type:
            return ANTIGEN_TYPE_TCELL_EXPECTATIONS.get(self.antigen_type)

        return None

    def is_tetramer_or_sct(self) -> bool:
        """Check if this sample is from tetramer or SCT selection."""
        return self.source in {"tetramer", "sct"} or (
            self.antigen_type is not None and (
                self.antigen_type.startswith("tetramer_") or self.antigen_type == "sct"
            )
        )

    def is_til(self) -> bool:
        """Check if this sample is TIL data."""
        return self.source == "til"

    def is_sct_data(self) -> bool:
        """Check if this sample is SCT platform data."""
        return self.source == "sct" or self.sct_path is not None


@dataclass
class SampleSheet:
    """Collection of samples with their metadata."""
    samples: list[Sample] = field(default_factory=list)

    def __len__(self):
        return len(self.samples)

    def __iter__(self):
        return iter(self.samples)

    def __getitem__(self, idx):
        return self.samples[idx]

    def get_sample(self, name: str) -> Sample | None:
        """Get a sample by name."""
        for s in self.samples:
            if s.sample == name:
                return s
        return None

    def get_culture_samples(self) -> list[Sample]:
        """Get all culture samples (not TIL, tetramer, or SCT)."""
        return [s for s in self.samples if s.source == "culture"]

    def get_til_samples(self) -> list[Sample]:
        """Get all TIL samples."""
        return [s for s in self.samples if s.is_til()]

    def get_tetramer_samples(self) -> list[Sample]:
        """Get all tetramer/SCT samples."""
        return [s for s in self.samples if s.is_tetramer_or_sct()]

    def get_sct_samples(self) -> list[Sample]:
        """Get all SCT platform samples."""
        return [s for s in self.samples if s.is_sct_data()]

    def get_samples_by_experiment(self) -> dict[str, list[Sample]]:
        """
        Group samples by experiment name.

        Returns a dictionary mapping experiment names to lists of samples.
        Samples without an experiment field are grouped under 'default'.
        """
        groups: dict[str, list[Sample]] = {}
        for s in self.samples:
            exp_name = s.experiment or "default"
            if exp_name not in groups:
                groups[exp_name] = []
            groups[exp_name].append(s)
        return groups

    def get_experiment_names(self) -> list[str]:
        """Get unique experiment names, preserving order."""
        seen = set()
        names = []
        for s in self.samples:
            exp_name = s.experiment or "default"
            if exp_name not in seen:
                seen.add(exp_name)
                names.append(exp_name)
        return names

    def to_dataframe(self) -> pd.DataFrame:
        """Convert sample sheet to a pandas DataFrame."""
        records = []
        for s in self.samples:
            record = {
                "sample": s.sample,
                "gex_dir": s.gex_dir,
                "vdj_dir": s.vdj_dir,
                "gex_path": s.gex_path,
                "sct_path": s.sct_path,
                "sct_sheet": s.sct_sheet,
                "experiment": s.experiment,
                "antigen_type": s.antigen_type,
                "antigen_description": s.antigen_description,
                "antigen_name": s.antigen_name,
                "epitope_sequence": s.epitope_sequence,
                "mhc_allele": s.mhc_allele,
                "culture_days": s.culture_days,
                "tcell_type_expected": s.tcell_type_expected,
                "pre_sorted": s.pre_sorted,
                "mhc_blocking": s.mhc_blocking,
                "source": s.source,
                "expected_tcell_type": s.get_expected_tcell_type(),
            }
            records.append(record)
        return pd.DataFrame(records)


def load_sample_sheet(path: str | Path) -> SampleSheet:
    """
    Load a sample sheet from CSV or YAML file.

    Parameters
    ----------
    path : str or Path
        Path to the sample sheet file (.csv, .tsv, .yaml, or .yml)

    Returns
    -------
    SampleSheet
        Parsed sample sheet with all samples
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Sample sheet not found: {path}")

    suffix = path.suffix.lower()

    if suffix in {".yaml", ".yml"}:
        return _load_yaml_sample_sheet(path)
    elif suffix in {".csv", ".tsv"}:
        return _load_csv_sample_sheet(path, sep="," if suffix == ".csv" else "\t")
    else:
        raise ValueError(f"Unsupported sample sheet format: {suffix}. Use .csv, .tsv, .yaml, or .yml")


def _load_yaml_sample_sheet(path: Path) -> SampleSheet:
    """Load sample sheet from YAML file."""
    try:
        import yaml
    except ImportError:
        raise ImportError("PyYAML is required to load YAML sample sheets. Install with: pip install pyyaml")

    with open(path) as f:
        data = yaml.safe_load(f)

    if not isinstance(data, dict) or "samples" not in data:
        raise ValueError("YAML sample sheet must have a 'samples' key with a list of samples")

    samples = []
    for sample_data in data["samples"]:
        # Handle nested peptide_pool if present
        if "peptide_pool" in sample_data and isinstance(sample_data["peptide_pool"], list):
            # Keep as list of dicts
            pass

        # Convert None strings to actual None
        for key, value in sample_data.items():
            if value == "null" or value == "":
                sample_data[key] = None

        samples.append(Sample(**sample_data))

    return SampleSheet(samples=samples)


def _load_csv_sample_sheet(path: Path, sep: str = ",") -> SampleSheet:
    """Load sample sheet from CSV/TSV file."""
    df = pd.read_csv(path, sep=sep)

    # Required column
    if "sample" not in df.columns:
        raise ValueError("Sample sheet must have a 'sample' column")

    # At least one of gex_dir or vdj_dir must be present
    if "gex_dir" not in df.columns and "vdj_dir" not in df.columns:
        raise ValueError("Sample sheet must have at least 'gex_dir' or 'vdj_dir' column")

    samples = []
    for _, row in df.iterrows():
        sample_data = {}
        for col in df.columns:
            value = row[col]
            # Convert NaN to None
            if pd.isna(value):
                value = None
            # Convert numeric strings for culture_days
            elif col == "culture_days" and value is not None:
                value = int(value)
            sample_data[col] = value

        samples.append(Sample(**sample_data))

    return SampleSheet(samples=samples)


def validate_sample_sheet(sample_sheet: SampleSheet) -> list[str]:
    """
    Validate a sample sheet and return any warnings.

    Parameters
    ----------
    sample_sheet : SampleSheet
        The sample sheet to validate

    Returns
    -------
    list[str]
        List of warning messages (empty if no warnings)
    """
    warnings = []

    # Check for duplicate sample names
    names = [s.sample for s in sample_sheet.samples]
    if len(names) != len(set(names)):
        duplicates = [n for n in names if names.count(n) > 1]
        warnings.append(f"Duplicate sample names found: {set(duplicates)}")

    # Check for samples with conflicting expectations
    for sample in sample_sheet.samples:
        expected = sample.get_expected_tcell_type()

        # Warn if long peptide but expecting CD8
        if sample.antigen_type == "long_peptide" and expected == "CD8":
            warnings.append(
                f"Sample '{sample.sample}': Long peptides typically favor CD4+ responses, "
                "but CD8+ is expected. Consider if this is correct."
            )

        # Warn if short peptide but expecting CD4
        if sample.antigen_type == "short_peptide" and expected == "CD4":
            warnings.append(
                f"Sample '{sample.sample}': Short peptides typically bind MHC-I for CD8+ responses, "
                "but CD4+ is expected. Consider if this is correct."
            )

        # Check for path existence
        if sample.gex_dir and not Path(sample.gex_dir).exists():
            warnings.append(f"Sample '{sample.sample}': gex_dir does not exist: {sample.gex_dir}")

        if sample.vdj_dir and not Path(sample.vdj_dir).exists():
            warnings.append(f"Sample '{sample.sample}': vdj_dir does not exist: {sample.vdj_dir}")

        if sample.sct_path and not Path(sample.sct_path).exists():
            warnings.append(f"Sample '{sample.sample}': sct_path does not exist: {sample.sct_path}")

        if sample.gex_path and not Path(sample.gex_path).exists():
            warnings.append(f"Sample '{sample.sample}': gex_path does not exist: {sample.gex_path}")

    return warnings
