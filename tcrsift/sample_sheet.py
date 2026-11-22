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
    "short_peptide": "CD8",  # 8-11aa, MHC-I direct binding
    "long_peptide": "mixed",  # 15-25+aa, requires processing, favors CD4
    "peptide_pool": "mixed",  # pool of peptides, mixed responses
    "minigene": "mixed",  # requires processing, favors CD4
    "minigene_library": "mixed",  # library of minigenes, mixed responses
    "whole_protein": "mixed",  # requires processing, favors CD4
    "mrna": "mixed",  # mRNA encoding antigen(s), requires processing
    "tetramer_mhc1": "CD8",  # MHC-I tetramer selection (single antigen)
    "tetramer_mhc2": "CD4",  # MHC-II tetramer selection (single antigen)
    "sct": "CD8",  # single-chain trimer (pMHC-I alpha-B2M-peptide fusion)
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
    # Free-form enrichment / selection method label. Distinct from `source`
    # (which is a fixed pipeline-source enum); used to stratify clones across
    # selection strategies within a donor (e.g. AIMpos, IFNpos, tetpos).
    enrichment_method: str | None = None
    # Free-form timepoint label for longitudinal designs (e.g. "D7", "D14",
    # "baseline", "post-IL2"). Distinct from `culture_days` (a single int);
    # used to group samples that are different timepoints of the same study
    # arm. (#9)
    timepoint: str | None = None
    # Free-form APC label for multi-APC stim designs (e.g. "mDC", "B-LCL",
    # "autologous_PBMC", "K562-A2-CD80"). Distinct from `tissue`
    # (anatomical) and `antigen_type` (what was presented). (#9)
    apc_type: str | None = None
    # Experiment grouping for multi-source unification
    experiment: str | None = None  # experiment name for unification (e.g., "TIL", "Culture", "SCT")
    # SCT-specific fields
    sct_path: str | None = None  # path to SCT Excel file
    sct_sheet: str = "Cell"  # sheet name in SCT Excel file
    # Standalone GEX file (for augmentation without full CellRanger output)
    gex_path: str | None = None  # path to 10x filtered_feature_bc_matrix.h5 file
    # TIL-specific input fields (alternative to vdj_dir for TIL samples)
    til_csv: str | None = None  # path to CSV with CDR3_alpha, CDR3_beta columns
    til_h5ad: str | None = None  # path to pre-processed h5ad file
    # Organism species ("human" or "mouse"; None inherits sheet/CLI auto-detection)
    species: str | None = None

    def __post_init__(self):
        if self.species is not None:
            from .species import normalize_species

            self.species = normalize_species(self.species, allow_auto=True)
        # Coerce culture_days consistently for both CSV and YAML inputs (which
        # may deliver str/float). Round rather than silently truncate, and fail
        # clearly on non-numeric instead of raising a bare ValueError later.
        if isinstance(self.culture_days, bool):
            raise ValueError(
                f"culture_days must be numeric for sample '{self.sample}', "
                f"got {self.culture_days!r}"
            )
        if self.culture_days is not None and not isinstance(self.culture_days, int):
            try:
                self.culture_days = int(round(float(self.culture_days)))
            except (TypeError, ValueError):
                raise ValueError(
                    f"culture_days must be numeric for sample '{self.sample}', "
                    f"got {self.culture_days!r}"
                )

        # Validate at least one data source
        has_cellranger = self.gex_dir or self.vdj_dir
        has_sct = self.sct_path is not None
        has_til_data = self.til_csv is not None or self.til_h5ad is not None
        if not has_cellranger and not has_sct and not has_til_data:
            raise ValueError(
                f"Sample '{self.sample}' must have at least gex_dir, vdj_dir, sct_path, til_csv, or til_h5ad"
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
            self.antigen_type is not None
            and (self.antigen_type.startswith("tetramer_") or self.antigen_type == "sct")
        )

    def is_til(self) -> bool:
        """Check if this sample is TIL data."""
        return self.source == "til"

    def is_sct_data(self) -> bool:
        """Check if this sample is SCT platform data."""
        return self.source == "sct" or self.sct_path is not None

    def get_til_data_source(self) -> tuple[str, str] | None:
        """
        Get the TIL data source type and path.

        Returns
        -------
        tuple[str, str] | None
            Tuple of (source_type, path) where source_type is one of:
            - "h5ad": Pre-processed AnnData file
            - "csv": CSV file with CDR3_alpha, CDR3_beta columns
            - "vdj_dir": CellRanger VDJ directory
            Returns None if this is not a TIL sample or has no data source.
        """
        if not self.is_til():
            return None

        if self.til_h5ad:
            return ("h5ad", self.til_h5ad)
        elif self.til_csv:
            return ("csv", self.til_csv)
        elif self.vdj_dir:
            return ("vdj_dir", self.vdj_dir)
        return None


@dataclass
class SampleSheet:
    """Collection of samples with their metadata."""

    samples: list[Sample] = field(default_factory=list)
    # Cohort-level flag (#26): set true when donors share antigen + MHC +
    # cohort and a unified FDR ranking across donors is biologically valid.
    # Read from the top-level YAML key `donors_share_antigen: true`.
    # Drives the auto-resolution of `--fdr-scope` to `global` instead of
    # `per-donor` for multi-donor cohorts.
    donors_share_antigen: bool = False
    species: str | None = None

    def __post_init__(self):
        from .species import normalize_species
        from .validation import TCRsiftValidationError

        if self.species is not None:
            self.species = normalize_species(self.species, allow_auto=True)
        explicit_sample_species = {
            s.species
            for s in self.samples
            if getattr(s, "species", None) not in (None, "auto")
        }
        if len(explicit_sample_species) > 1:
            raise TCRsiftValidationError(
                f"SampleSheet mixes multiple explicit species: {sorted(explicit_sample_species)}",
                hint="Run human and mouse cohorts in separate sample sheets.",
            )
        if len(explicit_sample_species) == 1:
            only_sp = next(iter(explicit_sample_species))
            if self.species in (None, "auto"):
                self.species = only_sp
            elif self.species != only_sp:
                raise TCRsiftValidationError(
                    f"SampleSheet species={self.species!r} conflicts with sample species={only_sp!r}",
                    hint="Ensure top-level and per-sample species agree.",
                )
        if self.species not in (None, "auto"):
            for s in self.samples:
                if s.species is None:
                    s.species = self.species

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

    def antigen_label_map(self, *, key: str = "sample") -> dict[str, str]:
        """Map ``{condition token → antigen label}`` from each sample's
        ``antigen_name`` (#317), for routing plot/table condition labels
        through :func:`tcrsift.format.pretty_antigen` so a pool/condition is
        never shown by a bare token.

        ``key`` is the :class:`Sample` attribute used as the condition token
        (default the ``sample`` name; e.g. ``"enrichment_method"`` to key by
        sort). The label is opaque free text (peptide pool, minigene, protein,
        …). Samples with no ``antigen_name`` are omitted; a multi-antigen
        ``antigen_names`` list is joined with " + " when the scalar
        ``antigen_name`` is absent.
        """
        out: dict[str, str] = {}
        for s in self.samples:
            token = getattr(s, key, None)
            if token is None:
                continue
            label = s.antigen_name
            if not label and getattr(s, "antigen_names", None):
                raw_names = s.antigen_names
                if isinstance(raw_names, str):
                    parts = [p.strip() for p in raw_names.split(";") if p.strip()]
                    label = " + ".join(parts)
                else:
                    label = " + ".join(str(a) for a in raw_names if a)
            if label:
                out[str(token)] = str(label)
        return out

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
                "patient_id": s.patient_id,
                "enrichment_method": s.enrichment_method,
                "timepoint": s.timepoint,
                "apc_type": s.apc_type,
                "tissue": s.tissue,
                "species": s.species,
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
        raise ValueError(
            f"Unsupported sample sheet format: {suffix}. Use .csv, .tsv, .yaml, or .yml"
        )


def _build_sample(data: dict) -> Sample:
    """Construct a Sample from a row/record, rejecting unknown columns.

    A typo'd or stray sample-sheet column would otherwise reach ``Sample(**data)``
    as an unexpected keyword and raise an opaque ``TypeError``; this surfaces a
    clear, actionable error instead. Shared by the YAML and CSV loaders.
    """
    from dataclasses import fields as _fields

    valid = {f.name for f in _fields(Sample)}
    unknown = sorted(k for k in data if k not in valid)
    if unknown:
        raise ValueError(
            f"Unknown sample-sheet column(s): {unknown}. "
            f"Valid columns: {sorted(valid)}"
        )
    # A blank cell means "not given": let the field's default apply instead of
    # passing None over it (a blank ``source`` was rejected as invalid).
    data = {key: value for key, value in data.items() if value is not None}
    for key in ("sample", "patient_id"):  # YAML reads a bare 12 as an integer
        if key in data:
            data[key] = str(data[key])
    return Sample(**data)


def _reject_duplicate_samples(samples: list) -> None:
    names = pd.Series([s.sample for s in samples])
    repeated = sorted(names[names.duplicated()].unique())
    if repeated:
        raise ValueError(f"Sample sheet lists sample name(s) more than once: {repeated}")


def _load_yaml_sample_sheet(path: Path) -> SampleSheet:
    """Load sample sheet from YAML file."""
    try:
        import yaml
    except ImportError:
        raise ImportError(
            "PyYAML is required to load YAML sample sheets. Install with: pip install pyyaml"
        )

    with open(path) as f:
        data = yaml.safe_load(f)

    if not isinstance(data, dict) or "samples" not in data:
        raise ValueError("YAML sample sheet must have a 'samples' key with a list of samples")

    raw_samples = data["samples"]
    if raw_samples is None:
        raw_samples = []
    elif not isinstance(raw_samples, list):
        raise ValueError("YAML sample sheet 'samples' key must map to a list of samples")

    samples = []
    for raw_item in raw_samples:
        if not isinstance(raw_item, dict):
            raise ValueError(f"Each sample entry in YAML must be a mapping, got {type(raw_item).__name__}")
        sample_data = dict(raw_item)
        # Translate legacy peptide_pool list if present before _build_sample checks columns
        if "peptide_pool" in sample_data:
            pool = sample_data.pop("peptide_pool")
            if isinstance(pool, list) and pool:
                names = []
                epitopes = []
                for entry in pool:
                    if isinstance(entry, dict):
                        if entry.get("name"):
                            names.append(str(entry["name"]))
                        if entry.get("sequence") or entry.get("epitope"):
                            epitopes.append(str(entry.get("sequence") or entry.get("epitope")))
                    elif entry:
                        names.append(str(entry))
                if names and not sample_data.get("antigen_names"):
                    sample_data["antigen_names"] = names
                if epitopes and not sample_data.get("epitope_sequences"):
                    sample_data["epitope_sequences"] = epitopes

        # Convert None strings to actual None
        for key, value in list(sample_data.items()):
            if value == "null" or value == "":
                sample_data[key] = None

        samples.append(_build_sample(sample_data))

    _reject_duplicate_samples(samples)
    # Sheet-level flags (#26). Accept any truthy YAML value — bare `true`
    # parses to bool, but quoted forms like "true" / "yes" parse to
    # strings; users shouldn't have their flag silently ignored over
    # quoting style.
    raw_dsa = data.get("donors_share_antigen")
    if isinstance(raw_dsa, str):
        donors_share_antigen = raw_dsa.strip().lower() in (
            "true", "yes", "y", "1", "on",
        )
    else:
        donors_share_antigen = bool(raw_dsa)
    raw_species = data.get("species")
    sheet_species = None
    if raw_species is not None and str(raw_species).strip():
        from .species import normalize_species

        sheet_species = normalize_species(str(raw_species), allow_auto=True)
    return SampleSheet(
        samples=samples,
        donors_share_antigen=donors_share_antigen,
        species=sheet_species,
    )


def _load_csv_sample_sheet(path: Path, sep: str = ",") -> SampleSheet:
    """Load sample sheet from CSV/TSV file."""
    # Force ID-bearing columns to string so pandas type inference doesn't read
    # e.g. sample "0012" as int 12 (leading zeros dropped) or "12" as 12.0,
    # which would collide distinct samples and diverge from the YAML loader.
    # (pandas ignores dtype keys for columns not present in the file.)
    df = pd.read_csv(path, sep=sep, dtype={"sample": str, "patient_id": str})

    # Required column
    if "sample" not in df.columns:
        raise ValueError("Sample sheet must have a 'sample' column")

    # At least one supported data-source column must be present
    source_cols = {"gex_dir", "vdj_dir", "gex_path", "sct_path", "til_csv", "til_h5ad"}
    if not (source_cols & set(df.columns)):
        raise ValueError(
            "Sample sheet must have at least 'gex_dir', 'vdj_dir', 'gex_path', "
            "'sct_path', 'til_csv', or 'til_h5ad' column"
        )

    samples = []
    for _, row in df.iterrows():
        sample_data = {}
        for col in df.columns:
            value = row[col]
            # Convert NaN to None; culture_days coercion is centralized in
            # Sample.__post_init__ so CSV and YAML behave identically.
            if pd.isna(value):
                value = None
            sample_data[col] = value

        samples.append(_build_sample(sample_data))

    _reject_duplicate_samples(samples)
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

        if sample.til_csv and not Path(sample.til_csv).exists():
            warnings.append(f"Sample '{sample.sample}': til_csv does not exist: {sample.til_csv}")

        if sample.til_h5ad and not Path(sample.til_h5ad).exists():
            warnings.append(f"Sample '{sample.sample}': til_h5ad does not exist: {sample.til_h5ad}")

    return warnings
