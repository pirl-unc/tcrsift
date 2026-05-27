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
Input validation for TCRsift.

Provides clear, actionable error messages when inputs don't meet requirements.
"""

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional, Union

import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

logger = logging.getLogger(__name__)


class TCRsiftValidationError(ValueError):
    """Custom exception for validation errors with clear messages."""

    def __init__(self, message: str, hint: Optional[str] = None):
        self.hint = hint
        full_message = message
        if hint:
            full_message += f"\n\nHint: {hint}"
        super().__init__(full_message)


def validate_file_exists(
    path: Union[str, Path],
    file_description: str = "file",
) -> Path:
    """
    Validate that a file exists and is readable.

    Parameters
    ----------
    path : str or Path
        Path to validate
    file_description : str
        Description for error messages (e.g., "sample sheet", "VDJ directory")

    Returns
    -------
    Path
        Validated Path object

    Raises
    ------
    TCRsiftValidationError
        If file doesn't exist or isn't readable
    """
    path = Path(path)

    if not path.exists():
        raise TCRsiftValidationError(
            f"The {file_description} does not exist: {path}",
            hint=f"Check that the path is correct. Current working directory: {Path.cwd()}",
        )

    if path.is_dir():
        raise TCRsiftValidationError(
            f"Expected a file for {file_description}, but got a directory: {path}",
            hint="Provide the path to a file, not a directory.",
        )

    return path


def validate_directory_exists(
    path: Union[str, Path],
    dir_description: str = "directory",
    required_files: Optional[list[str]] = None,
) -> Path:
    """
    Validate that a directory exists and optionally contains required files.

    Parameters
    ----------
    path : str or Path
        Path to validate
    dir_description : str
        Description for error messages
    required_files : list of str, optional
        Files that must exist in the directory

    Returns
    -------
    Path
        Validated Path object
    """
    path = Path(path)

    if not path.exists():
        raise TCRsiftValidationError(
            f"The {dir_description} does not exist: {path}",
            hint=f"Check that the path is correct. Current working directory: {Path.cwd()}",
        )

    if not path.is_dir():
        raise TCRsiftValidationError(
            f"Expected a directory for {dir_description}, but got a file: {path}",
            hint="Provide the path to a directory, not a file.",
        )

    if required_files:
        missing = [f for f in required_files if not (path / f).exists()]
        if missing:
            available = [f.name for f in path.iterdir()][:10]
            raise TCRsiftValidationError(
                f"The {dir_description} is missing required files: {missing}",
                hint=f"Available files in directory: {available}. "
                "Make sure this is the correct CellRanger output directory.",
            )

    return path


def validate_dataframe(
    df: Any,
    name: str = "DataFrame",
    required_columns: Optional[list[str]] = None,
    min_rows: int = 0,
) -> pd.DataFrame:
    """
    Validate that input is a DataFrame with required columns.

    Parameters
    ----------
    df : Any
        Object to validate
    name : str
        Name for error messages
    required_columns : list of str, optional
        Columns that must be present
    min_rows : int
        Minimum number of rows required

    Returns
    -------
    pd.DataFrame
        Validated DataFrame
    """
    if not isinstance(df, pd.DataFrame):
        raise TCRsiftValidationError(
            f"Expected {name} to be a pandas DataFrame, but got {type(df).__name__}",
            hint="Make sure you're passing a DataFrame, not a file path or other object.",
        )

    if len(df) < min_rows:
        raise TCRsiftValidationError(
            f"{name} has {len(df)} rows, but at least {min_rows} are required",
            hint="Check that your input data is not empty and contains valid entries.",
        )

    if required_columns:
        missing = [col for col in required_columns if col not in df.columns]
        if missing:
            available = list(df.columns)[:20]
            raise TCRsiftValidationError(
                f"{name} is missing required columns: {missing}",
                hint=f"Available columns: {available}. "
                "Check that your data has the expected format.",
            )

    return df


def validate_anndata(
    adata: Any,
    name: str = "AnnData",
    required_obs_columns: Optional[list[str]] = None,
    min_cells: int = 0,
) -> "ad.AnnData":
    """
    Validate AnnData object.

    Parameters
    ----------
    adata : Any
        Object to validate
    name : str
        Name for error messages
    required_obs_columns : list of str, optional
        Columns required in adata.obs
    min_cells : int
        Minimum number of cells required

    Returns
    -------
    ad.AnnData
        Validated AnnData
    """
    import anndata as ad

    if not isinstance(adata, ad.AnnData):
        raise TCRsiftValidationError(
            f"Expected {name} to be an AnnData object, but got {type(adata).__name__}",
            hint="Load your data with scanpy or anndata first, or use tcrsift.load_samples().",
        )

    if adata.n_obs < min_cells:
        raise TCRsiftValidationError(
            f"{name} has {adata.n_obs} cells, but at least {min_cells} are required",
            hint="Check that your input data loaded correctly and contains cells.",
        )

    if required_obs_columns:
        missing = [col for col in required_obs_columns if col not in adata.obs.columns]
        if missing:
            available = list(adata.obs.columns)[:20]
            raise TCRsiftValidationError(
                f"{name} is missing required columns in .obs: {missing}",
                hint=f"Available columns in .obs: {available}. "
                "Did you run the previous pipeline steps?",
            )

    return adata


def validate_cellranger_vdj_dir(path: Union[str, Path]) -> Path:
    """
    Validate a CellRanger VDJ output directory.

    Parameters
    ----------
    path : str or Path
        Path to CellRanger VDJ output directory

    Returns
    -------
    Path
        Validated path
    """
    path = Path(path)

    # Check directory exists
    path = validate_directory_exists(path, "CellRanger VDJ directory")

    # Check for key files
    contig_files = [
        "filtered_contig_annotations.csv",
        "all_contig_annotations.csv",
    ]

    found_contig = any((path / f).exists() for f in contig_files)
    if not found_contig:
        available = [f.name for f in path.iterdir()][:15]
        raise TCRsiftValidationError(
            f"No contig annotations file found in CellRanger VDJ directory: {path}",
            hint=f"Expected one of: {contig_files}. "
            f"Available files: {available}. "
            "This should be the 'outs' directory from cellranger vdj.",
        )

    return path


def validate_cellranger_gex_dir(path: Union[str, Path]) -> Path:
    """
    Validate a CellRanger GEX output directory.

    Parameters
    ----------
    path : str or Path
        Path to CellRanger GEX output directory

    Returns
    -------
    Path
        Validated path
    """
    path = Path(path)

    # Check directory exists
    path = validate_directory_exists(path, "CellRanger GEX directory")

    # Check for key files/directories
    matrix_locations = [
        "filtered_feature_bc_matrix",
        "filtered_feature_bc_matrix.h5",
        "raw_feature_bc_matrix",
    ]

    found_matrix = any((path / f).exists() for f in matrix_locations)
    if not found_matrix:
        available = [f.name for f in path.iterdir()][:15]
        raise TCRsiftValidationError(
            f"No gene expression matrix found in CellRanger GEX directory: {path}",
            hint=f"Expected one of: {matrix_locations}. "
            f"Available files: {available}. "
            "This should be the 'outs' directory from cellranger count.",
        )

    return path


def validate_cdr3_sequence(
    seq: str,
    chain: str = "unknown",
    strict: bool = False,
) -> bool:
    """
    Validate a CDR3 amino acid sequence.

    Parameters
    ----------
    seq : str
        CDR3 sequence to validate
    chain : str
        Chain type for error messages ("alpha" or "beta")
    strict : bool
        If True, raise on invalid; if False, return bool

    Returns
    -------
    bool
        True if valid

    Notes
    -----
    This function checks:
    - Characters are valid amino acid letters (ACDEFGHIKLMNPQRSTVWY) or stop codon (*)
    - Warns (debug level) if sequence doesn't start with C
    - Warns (debug level) if sequence length is outside 5-30 aa range
    """
    if pd.isna(seq) or seq == "":
        return True  # Missing is OK

    # Valid amino acid characters
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY*")

    invalid_chars = set(seq.upper()) - valid_aa
    if invalid_chars:
        msg = f"CDR3 {chain} sequence contains invalid characters: {invalid_chars} in '{seq}'"
        if strict:
            raise TCRsiftValidationError(
                msg,
                hint="CDR3 sequences should contain only standard amino acid letters.",
            )
        logger.warning(msg)
        return False

    # Check typical CDR3 patterns
    seq_upper = seq.upper()
    if chain == "alpha" and not seq_upper.startswith("C"):
        logger.debug(f"CDR3 alpha doesn't start with C: {seq}")
    elif chain == "beta" and not seq_upper.startswith("C"):
        logger.debug(f"CDR3 beta doesn't start with C: {seq}")

    # Check length bounds (typical CDR3 is 5-30 aa)
    if len(seq) < 5:
        logger.debug(f"CDR3 {chain} unusually short ({len(seq)} aa): {seq}")
    elif len(seq) > 30:
        logger.debug(f"CDR3 {chain} unusually long ({len(seq)} aa): {seq}")

    return True


def validate_cdr3_dataframe(
    df: pd.DataFrame,
    alpha_col: str = "CDR3_alpha",
    beta_col: str = "CDR3_beta",
    strict: bool = False,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Validate CDR3 sequences in a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing CDR3 sequences
    alpha_col : str
        Column name for alpha chain CDR3
    beta_col : str
        Column name for beta chain CDR3
    strict : bool
        If True, raise on first invalid sequence

    Returns
    -------
    tuple[pd.DataFrame, list[str]]
        - DataFrame with 'cdr3_valid' column added
        - List of warning messages for invalid sequences

    Examples
    --------
    >>> df, warnings = validate_cdr3_dataframe(clonotypes)
    >>> if warnings:
    ...     print(f"Found {len(warnings)} invalid CDR3 sequences")
    >>> valid_df = df[df['cdr3_valid']]
    """
    warnings = []
    valid_mask = pd.Series(True, index=df.index)

    for col, chain in [(alpha_col, "alpha"), (beta_col, "beta")]:
        if col not in df.columns:
            continue

        for idx, seq in df[col].items():
            if not validate_cdr3_sequence(seq, chain=chain, strict=strict):
                valid_mask[idx] = False
                warnings.append(f"Row {idx}: Invalid CDR3 {chain} sequence: {seq}")

    df = df.copy()
    df["cdr3_valid"] = valid_mask

    return df, warnings


# =============================================================================
# Safe Calculation Utilities
# =============================================================================


def safe_divide(
    numerator: float,
    denominator: float,
    default: float = 0.0,
) -> float:
    """
    Safely divide two numbers, returning a default value if denominator is zero.

    Parameters
    ----------
    numerator : float
        The numerator
    denominator : float
        The denominator
    default : float
        Value to return if denominator is zero (default: 0.0)

    Returns
    -------
    float
        Result of division or default value

    Examples
    --------
    >>> safe_divide(10, 2)
    5.0
    >>> safe_divide(10, 0)
    0.0
    >>> safe_divide(10, 0, default=float('nan'))
    nan
    """
    if denominator == 0:
        return default
    return numerator / denominator


def safe_percentage(
    part: float,
    total: float,
    default: float = 0.0,
) -> float:
    """
    Safely calculate a percentage, returning default if total is zero.

    Parameters
    ----------
    part : float
        The part (numerator)
    total : float
        The total (denominator)
    default : float
        Value to return if total is zero (default: 0.0)

    Returns
    -------
    float
        Percentage (0-100) or default value

    Examples
    --------
    >>> safe_percentage(25, 100)
    25.0
    >>> safe_percentage(10, 0)
    0.0
    """
    if total == 0:
        return default
    return (part / total) * 100


def safe_mode(series: pd.Series, default: Any = None) -> Any:
    """
    Safely get the mode of a pandas Series, returning default if empty.

    This is useful when aggregating data where some groups may be empty
    or have no valid values.

    Parameters
    ----------
    series : pd.Series
        Series to get mode from
    default : Any
        Value to return if series is empty or has no mode (default: None)

    Returns
    -------
    Any
        The mode value or default

    Examples
    --------
    >>> safe_mode(pd.Series(['A', 'A', 'B']))
    'A'
    >>> safe_mode(pd.Series([]))
    None
    >>> safe_mode(pd.Series([np.nan, np.nan]), default='Unknown')
    'Unknown'
    """
    if series is None or len(series) == 0:
        return default

    # Drop NA values before computing mode
    clean_series = series.dropna()
    if len(clean_series) == 0:
        return default

    mode_result = clean_series.mode()
    if len(mode_result) == 0:
        return default

    return mode_result.iloc[0]


def pick_representative_cell(
    clone_df: pd.DataFrame,
    *,
    umi_cols: tuple[str, ...] = ("TRA_1_umis", "TRB_1_umis"),
) -> pd.Series | None:
    """Pick one row from a per-cell DataFrame to act as the canonical
    source for coupled per-contig columns (#94).

    Ranks rows by summed UMI evidence across ``umi_cols`` (default:
    α+β UMI sums) and returns the top-ranked row as a :class:`pd.Series`.
    Ties are broken by the row's position in ``clone_df`` (first wins).
    When no ``umi_cols`` are present in the frame, falls back to the
    first row — that's still single-cell-sourced, just ranked by a
    weaker (positional) signal.

    Use this whenever aggregating multiple correlated columns from a
    group of cells (CDR3 AA, CDR3 NT, V/J/C gene calls, contig IDs,
    leader, constant region — anything that should describe a single
    biological molecule). Never call :func:`safe_mode` independently on
    coupled columns: it picks per-column modes with lex tie-breaking,
    so AA and NT can end up sourced from different cells and produce
    an internally-inconsistent (AA, NT) pair.

    Parameters
    ----------
    clone_df
        Per-cell rows for a single clonotype (or any grouping where
        you want to copy multiple correlated columns from one source).
    umi_cols
        Columns to sum for the UMI-evidence ranking. Defaults to the
        CellRanger α+β UMI columns. Pass a single-element tuple to
        rank by one chain only, or supply a different ranking proxy.

    Returns
    -------
    pd.Series | None
        The representative row, or ``None`` when ``clone_df`` is empty.
    """
    if clone_df is None or len(clone_df) == 0:
        return None
    available = [c for c in umi_cols if c in clone_df.columns]
    if available:
        score = clone_df[available].fillna(0).sum(axis=1)
        rep_idx = score.idxmax()
    else:
        rep_idx = clone_df.index[0]
    return clone_df.loc[rep_idx]


def validate_sample_sheet_entry(
    entry: dict,
    index: int,
) -> list[str]:
    """
    Validate a single sample sheet entry.

    Parameters
    ----------
    entry : dict
        Sample entry to validate
    index : int
        Entry index for error messages

    Returns
    -------
    list of str
        List of warning messages (empty if valid)
    """
    warnings = []

    # Required fields
    if "sample_name" not in entry and "name" not in entry:
        raise TCRsiftValidationError(
            f"Sample sheet entry {index + 1} is missing 'sample_name' or 'name' field",
            hint="Each sample must have a name. Check your sample sheet format.",
        )

    # Check paths exist
    for path_field in ["vdj_path", "gex_path"]:
        if path_field in entry:
            path = Path(entry[path_field])
            if not path.exists():
                warnings.append(f"Sample {index + 1}: {path_field} does not exist: {path}")

    # Check for unusual values
    sample_name = entry.get("sample_name") or entry.get("name", "")
    if not sample_name.strip():
        warnings.append(f"Sample {index + 1} has empty name")

    return warnings


def validate_clonotype_df(
    df: pd.DataFrame,
    for_filtering: bool = False,
    for_annotation: bool = False,
    for_assembly: bool = False,
) -> pd.DataFrame:
    """
    Validate a clonotype DataFrame for specific operations.

    Parameters
    ----------
    df : pd.DataFrame
        Clonotype DataFrame to validate
    for_filtering : bool
        Validate for filtering operations
    for_annotation : bool
        Validate for annotation operations
    for_assembly : bool
        Validate for assembly operations

    Returns
    -------
    pd.DataFrame
        Validated DataFrame
    """
    # Basic validation
    df = validate_dataframe(df, "Clonotype DataFrame", min_rows=1)

    # Check for clone identifier
    if "CDR3ab" not in df.columns and "CDR3_beta" not in df.columns:
        raise TCRsiftValidationError(
            "Clonotype DataFrame must have 'CDR3ab' or 'CDR3_beta' column",
            hint="Make sure you're using output from tcrsift.aggregate_clonotypes().",
        )

    if for_filtering:
        if "cell_count" not in df.columns:
            raise TCRsiftValidationError(
                "Filtering requires 'cell_count' column in clonotype DataFrame",
                hint="This column is created by aggregate_clonotypes(). "
                "Make sure you're using the correct input file.",
            )

    if for_annotation:
        cdr3_cols = ["CDR3_alpha", "CDR3_beta"]
        if not any(col in df.columns for col in cdr3_cols):
            raise TCRsiftValidationError(
                "Annotation requires CDR3 sequence columns (CDR3_alpha, CDR3_beta)",
                hint="Make sure your clonotype file has CDR3 sequence columns.",
            )

    if for_assembly:
        required = ["CDR3_alpha", "CDR3_beta"]
        missing = [col for col in required if col not in df.columns]
        if missing:
            raise TCRsiftValidationError(
                f"Assembly requires columns: {required}. Missing: {missing}",
                hint="Full-length assembly needs both alpha and beta CDR3 sequences.",
            )

    return df


def validate_numeric_param(
    value: Any,
    name: str,
    min_value: Optional[float] = None,
    max_value: Optional[float] = None,
) -> float:
    """
    Validate a numeric parameter.

    Parameters
    ----------
    value : Any
        Value to validate
    name : str
        Parameter name for error messages
    min_value : float, optional
        Minimum allowed value
    max_value : float, optional
        Maximum allowed value

    Returns
    -------
    float
        Validated value
    """
    try:
        value = float(value)
    except (TypeError, ValueError):
        raise TCRsiftValidationError(
            f"Parameter '{name}' must be a number, got: {value} ({type(value).__name__})",
            hint="Check that you're passing a valid number.",
        )

    if min_value is not None and value < min_value:
        raise TCRsiftValidationError(
            f"Parameter '{name}' must be >= {min_value}, got: {value}",
            hint=f"Adjust the value to be at least {min_value}.",
        )

    if max_value is not None and value > max_value:
        raise TCRsiftValidationError(
            f"Parameter '{name}' must be <= {max_value}, got: {value}",
            hint=f"Adjust the value to be at most {max_value}.",
        )

    return value


def log_validation_summary(
    n_valid: int,
    n_total: int,
    item_name: str = "items",
    warnings: Optional[list[str]] = None,
):
    """
    Log a summary of validation results.

    Parameters
    ----------
    n_valid : int
        Number of valid items
    n_total : int
        Total number of items
    item_name : str
        Name of items for message
    warnings : list of str, optional
        Warning messages to log
    """
    pct = n_valid / n_total * 100 if n_total > 0 else 0
    logger.info(f"Validated {n_valid}/{n_total} {item_name} ({pct:.1f}% passed)")

    if warnings:
        for w in warnings[:10]:  # Limit to first 10 warnings
            logger.warning(w)
        if len(warnings) > 10:
            logger.warning(f"... and {len(warnings) - 10} more warnings")


# =============================================================================
# CLI Argument Validation
# =============================================================================


def parse_til_sample_spec(spec: str) -> tuple[str, str, Path]:
    """
    Parse a --til-sample specification.

    Expected formats
    ----------------
    - NAME=TYPE:PATH
    - TYPE:PATH  (sample name inferred from PATH stem)

    TYPE aliases
    ------------
    - csv
    - h5ad
    - vdj, vdj_dir, vdj-dir  -> vdj_dir
    """
    if not isinstance(spec, str) or not spec.strip():
        raise TCRsiftValidationError(
            "Empty --til-sample specification",
            hint="Use format NAME=TYPE:PATH, e.g. T1=csv:/path/to/til.csv",
        )

    spec = spec.strip()
    sample_name = None
    source_spec = spec

    if "=" in spec:
        sample_name, source_spec = spec.split("=", 1)
        sample_name = sample_name.strip()
        if not sample_name:
            raise TCRsiftValidationError(
                f"Invalid --til-sample spec: '{spec}'",
                hint="Sample name before '=' cannot be empty.",
            )

    if ":" not in source_spec:
        raise TCRsiftValidationError(
            f"Invalid --til-sample spec: '{spec}'",
            hint="Expected TYPE:PATH after optional NAME=. Example: T1=csv:/path/to/til.csv",
        )

    source_type_raw, path_raw = source_spec.split(":", 1)
    source_type_raw = source_type_raw.strip().lower()
    path_raw = path_raw.strip()

    source_aliases = {
        "csv": "csv",
        "h5ad": "h5ad",
        "vdj": "vdj_dir",
        "vdj_dir": "vdj_dir",
        "vdj-dir": "vdj_dir",
    }
    source_type = source_aliases.get(source_type_raw)
    if source_type is None:
        raise TCRsiftValidationError(
            f"Invalid TIL source type in --til-sample: '{source_type_raw}'",
            hint="Use one of: csv, h5ad, vdj (or vdj_dir).",
        )

    if not path_raw:
        raise TCRsiftValidationError(
            f"Missing path in --til-sample spec: '{spec}'",
            hint="Expected non-empty path after TYPE:.",
        )

    path = Path(path_raw).expanduser()

    if sample_name is None:
        sample_name = path.stem if path.suffix else path.name
        if not sample_name:
            sample_name = "TIL"

    return sample_name, source_type, path


def validate_til_sample_specs(specs: list[str]) -> list[tuple[str, str, Path]]:
    """
    Validate and parse --til-sample specifications.

    Returns parsed tuples: (sample_name, source_type, path).
    """
    if not specs:
        raise TCRsiftValidationError(
            "No --til-sample values provided",
            hint="Provide at least one spec, e.g. --til-sample T1=csv:/path/to/til.csv",
        )

    parsed_specs: list[tuple[str, str, Path]] = []
    seen_names: set[str] = set()

    for spec in specs:
        sample_name, source_type, path = parse_til_sample_spec(spec)

        if sample_name in seen_names:
            raise TCRsiftValidationError(
                f"Duplicate TIL sample name in --til-sample specs: '{sample_name}'",
                hint="Use unique sample names, e.g. T1=..., T2=...",
            )
        seen_names.add(sample_name)

        if source_type in {"csv", "h5ad"}:
            validate_file_exists(path, f"TIL {source_type} file for sample '{sample_name}'")
        elif source_type == "vdj_dir":
            validate_directory_exists(path, f"TIL VDJ directory for sample '{sample_name}'")

        parsed_specs.append((sample_name, source_type, path))

    return parsed_specs


def validate_cli_conditional_requirement(
    args,
    required_arg: str,
    condition_args: list[str],
    condition_values: Optional[list] = None,
    condition_description: str = "",
) -> None:
    """
    Validate that a conditionally required CLI argument is provided.

    Use this for arguments that are only required when certain other
    arguments are set to specific values.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments
    required_arg : str
        Name of the argument that is conditionally required
    condition_args : list of str
        Names of arguments that trigger the requirement
    condition_values : list, optional
        If provided, the condition is only met when any condition_arg
        equals any of these values. If None, any truthy value triggers it.
    condition_description : str
        Human-readable description of when the argument is required

    Raises
    ------
    TCRsiftValidationError
        If the condition is met but the required argument is missing

    Examples
    --------
    # --contigs-dir required when --leaders-from-contigs is set
    validate_cli_conditional_requirement(
        args,
        required_arg="contigs_dir",
        condition_args=["leaders_from_contigs"],
        condition_description="when using --leaders-from-contigs",
    )

    # --contigs-dir required when --alpha-leader or --beta-leader is "from_contig"
    validate_cli_conditional_requirement(
        args,
        required_arg="contigs_dir",
        condition_args=["alpha_leader", "beta_leader"],
        condition_values=["from_contig"],
        condition_description="when using leader='from_contig'",
    )
    """
    # Check if any condition is met
    condition_met = False
    triggering_args = []

    for cond_arg in condition_args:
        value = getattr(args, cond_arg, None)
        if condition_values is not None:
            # Check if value matches any of the specified values
            if value in condition_values:
                condition_met = True
                triggering_args.append(f"--{cond_arg.replace('_', '-')}={value}")
        else:
            # Any truthy value triggers the condition
            if value:
                condition_met = True
                triggering_args.append(f"--{cond_arg.replace('_', '-')}")

    if not condition_met:
        return  # Condition not met, no requirement to check

    # Check if required argument is provided
    required_value = getattr(args, required_arg, None)
    if required_value is None or required_value == "":
        arg_flag = f"--{required_arg.replace('_', '-')}"
        raise TCRsiftValidationError(
            f"{arg_flag} is required {condition_description}",
            hint=f"You specified {', '.join(triggering_args)}, which requires {arg_flag}. "
            f"Either provide {arg_flag} or use a different option.",
        )


def validate_cli_mutually_exclusive(
    args,
    arg_names: list[str],
    group_description: str = "",
) -> None:
    """
    Validate that at most one of a set of mutually exclusive arguments is set.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments
    arg_names : list of str
        Names of mutually exclusive arguments
    group_description : str
        Human-readable description of the group

    Raises
    ------
    TCRsiftValidationError
        If more than one argument in the group is set
    """
    set_args = []
    for arg_name in arg_names:
        value = getattr(args, arg_name, None)
        if value is not None and value is not False:
            set_args.append(f"--{arg_name.replace('_', '-')}")

    if len(set_args) > 1:
        raise TCRsiftValidationError(
            f"Cannot use {' and '.join(set_args)} together",
            hint=f"These options are mutually exclusive{': ' + group_description if group_description else ''}. "
            "Choose only one.",
        )


def validate_assemble_args(args) -> None:
    """
    Validate arguments for the assemble command.

    Checks conditional requirements specific to sequence assembly.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments for assemble command

    Raises
    ------
    TCRsiftValidationError
        If required arguments are missing
    """
    # Check if --leaders-from-contigs requires --contigs-dir
    validate_cli_conditional_requirement(
        args,
        required_arg="contigs_dir",
        condition_args=["leaders_from_contigs"],
        condition_description="when using --leaders-from-contigs",
    )

    # Check if --alpha-leader=from_contig or --beta-leader=from_contig requires --contigs-dir
    validate_cli_conditional_requirement(
        args,
        required_arg="contigs_dir",
        condition_args=["alpha_leader", "beta_leader"],
        condition_values=["from_contig"],
        condition_description="when using --alpha-leader=from_contig or --beta-leader=from_contig",
    )

    # Validate contigs-dir exists if provided
    if getattr(args, "contigs_dir", None):
        validate_directory_exists(args.contigs_dir, "contigs directory")


def validate_annotate_gex_args(args) -> None:
    """
    Validate arguments for the annotate-gex command.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments for annotate-gex command

    Raises
    ------
    TCRsiftValidationError
        If required arguments are missing or invalid
    """
    # Validate GEX file exists
    validate_file_exists(args.gex_file, "gene expression file (--gex-file)")


def validate_run_args(args) -> None:
    """
    Validate arguments for the run command (unified pipeline).

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments for run command

    Raises
    ------
    TCRsiftValidationError
        If required arguments are missing
    """
    # Check if --leaders-from-contigs requires --contigs-dir
    validate_cli_conditional_requirement(
        args,
        required_arg="contigs_dir",
        condition_args=["leaders_from_contigs"],
        condition_description="when using --leaders-from-contigs",
    )

    # Check if --alpha-leader=from_contig or --beta-leader=from_contig requires --contigs-dir
    validate_cli_conditional_requirement(
        args,
        required_arg="contigs_dir",
        condition_args=["alpha_leader", "beta_leader"],
        condition_values=["from_contig"],
        condition_description="when using --alpha-leader=from_contig or --beta-leader=from_contig",
    )

    # Validate sample sheet exists
    validate_file_exists(args.sample_sheet, "sample sheet (--sample-sheet)")

    # Prevent ambiguous TIL source selection in run command
    if getattr(args, "til_samples", None) and getattr(args, "til_sample", None):
        raise TCRsiftValidationError(
            "Conflicting TIL options: --til-samples and --til-sample",
            hint="Use either --til-samples (sample names from sample sheet) "
            "or repeat --til-sample (direct source specs), not both.",
        )

    # Validate contigs-dir exists if provided
    if getattr(args, "contigs_dir", None):
        validate_directory_exists(args.contigs_dir, "contigs directory (--contigs-dir)")

    # Validate direct TIL sample specs if provided
    if getattr(args, "til_sample", None):
        validate_til_sample_specs(args.til_sample)


def validate_match_til_args(args) -> None:
    """
    Validate arguments for the match-til command.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments for match-til command

    Raises
    ------
    TCRsiftValidationError
        If required arguments are missing or invalid
    """
    # Check that at least one TIL data source is provided
    til_sources = ["sample_sheet", "til_h5ad", "til_csv", "til_vdj_dir", "til_sample"]
    provided_sources = []
    for src in til_sources:
        value = getattr(args, src, None)
        if src == "til_sample":
            if value:
                provided_sources.append(src)
        elif value:
            provided_sources.append(src)

    if not provided_sources:
        raise TCRsiftValidationError(
            "No TIL data source specified",
            hint="Provide one of: --sample-sheet, --til-h5ad, --til-csv, --til-vdj-dir, or --til-sample",
        )

    # Validate that only one source is provided (to avoid confusion)
    if len(provided_sources) > 1:
        flags = [f"--{src.replace('_', '-')}" for src in provided_sources]
        raise TCRsiftValidationError(
            f"Multiple TIL data sources specified: {', '.join(flags)}",
            hint="Provide only one TIL data source. For multiple TIL samples, use --sample-sheet or repeat --til-sample.",
        )

    # Validate the provided source exists
    source = provided_sources[0]
    source_path = getattr(args, source)

    if source == "sample_sheet":
        validate_file_exists(source_path, "TIL sample sheet (--sample-sheet)")
    elif source == "til_h5ad":
        validate_file_exists(source_path, "TIL h5ad file (--til-h5ad)")
    elif source == "til_csv":
        validate_file_exists(source_path, "TIL CSV file (--til-csv)")
    elif source == "til_vdj_dir":
        validate_directory_exists(source_path, "TIL VDJ directory (--til-vdj-dir)")
    elif source == "til_sample":
        validate_til_sample_specs(source_path)
