# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
"""Clone-size-matched background tails for expression-state prioritization.

This competitive reference asks whether a clone exceeds random groups of cells
from its own sample/lineage. It is not a technical-noise model or a classifier
of antigen specificity. See docs/user-guide/til-signatures.md.
"""

from __future__ import annotations

import hashlib
import itertools
import json

import numpy as np
import pandas as pd

STRATUM = ["donor", "sample", "lineage"]


def _small_population(n, size, limit):
    """Whether all subsets fit the draw budget, without huge binomial ints."""
    count = 1
    for k in range(1, min(size, n - size) + 1):
        count = count * (n - k + 1) // k
        if count > limit:
            return False
    return True


def _background_means(matrix, size, draws, seed):
    """Draw cell subsets without replacement; enumerate small populations."""
    n = len(matrix)
    if size == 1:
        return matrix.copy(), "exact"
    if size == n:
        return matrix.mean(axis=0, keepdims=True), "exact"
    exact = _small_population(n, size, draws)
    rng = np.random.default_rng(seed)
    subsets = (itertools.combinations(range(n), size) if exact else
               (rng.choice(n, size, replace=False) for _ in range(draws)))
    # One subset at a time bounds memory even for a very large clone. Reuse
    # draws across signatures, preserving all observed within-cell covariance.
    means = [matrix[np.asarray(indices)].mean(axis=0) for indices in subsets]
    return np.asarray(means), "exact" if exact else "monte_carlo"


def validate_background_parameters(quantile, draws, seed):
    """Reject unusable calibration settings before loading expression data."""
    if not 0.5 < quantile < 1:
        raise ValueError("--signature-background-quantile must be between 0.5 and 1 (exclusive)")
    if draws < 1000 or draws * (1 - quantile) < 10 - 1e-9:
        raise ValueError("--signature-background-draws must be >= 1000 and leave >= 10 draws in the upper tail")
    if seed < 0:
        raise ValueError("--signature-background-seed must be >= 0")


def calibrate_signature_background(scores, per_cell, signatures, *, quantile=0.99, draws=2000, seed=0):
    """Return clone scores with floors/tail probabilities and a calibration table.

    Cells must be the same population used to form the clone means and counts.
    Each distinct clone size receives its own reference within a stratum. Tail
    probabilities include ties and use a +1 correction for Monte Carlo draws;
    they are unadjusted, competitive probabilities, not antigen-reactivity FDRs.
    """
    validate_background_parameters(quantile, draws, seed)
    result = scores.copy()
    columns = [f"signature_{name}" for name in signatures]
    for col in columns:
        result[f"{col}_noise_floor"] = np.nan
        result[f"{col}_background_tail_probability"] = np.nan
    cell_groups = per_cell.groupby(STRATUM, observed=True)
    records = []
    for key, group in result.groupby(STRATUM, observed=True, sort=True):
        cells = cell_groups.get_group(key).sort_index()
        counts = cells.groupby("CDR3ab", observed=True).size()
        if not group.CDR3ab.map(counts).eq(group.cells).all():
            raise ValueError("Background cells must match the cells used for clone scores and counts")
        matrix = cells[columns].to_numpy(dtype=float)
        finite = np.isfinite(matrix).all(axis=0)
        for size, same_size in group.groupby("cells", observed=True, sort=True):
            size = int(size)
            token = json.dumps([seed, *map(str, key), size]).encode()
            local_seed = int.from_bytes(hashlib.sha256(token).digest()[:8], "little")
            means, method = _background_means(matrix, size, draws, local_seed)
            for j, (name, col) in enumerate(zip(signatures, columns)):
                null = means[:, j]
                floor = float(np.quantile(null, quantile, method="higher")) if finite[j] else np.nan
                result.loc[same_size.index, f"{col}_noise_floor"] = floor
                if finite[j]:
                    # Include numerical ties so an all-equal program cannot
                    # acquire a tiny tail probability from reduction roundoff.
                    ordered = np.sort(null)
                    observed = same_size[col].to_numpy(dtype=float)
                    tolerance = 1e-12 + 1e-9 * np.abs(observed)
                    exceed = len(null) - np.searchsorted(ordered, observed - tolerance, side="left")
                    correction = int(method == "monte_carlo")
                    tail = (exceed + correction) / (len(null) + correction)
                    result.loc[same_size.index, f"{col}_background_tail_probability"] = tail
                records.append(dict(zip(STRATUM, map(str, key))) | {
                    "signature": name, "clone_cells": size, "background_cells": len(cells),
                    "quantile": quantile, "noise_floor": floor, "draws": len(null),
                    "method": method, "status": "ok" if finite[j] else "unavailable",
                })
    return result, pd.DataFrame(records)
