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

    return True


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
                warnings.append(
                    f"Sample {index + 1}: {path_field} does not exist: {path}"
                )

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
    if "clone_id" not in df.columns and "CDR3_beta" not in df.columns:
        raise TCRsiftValidationError(
            "Clonotype DataFrame must have 'clone_id' or 'CDR3_beta' column",
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
        cdr3_cols = ["CDR3_alpha", "CDR3_beta", "CDR3_beta"]
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
