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

"""One-parameter ridge-regularized logistic frequency model for FDR tiers."""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class LogisticFrequencyModel:
    """Small logistic model compatible with the needs of threshold selection."""

    params: np.ndarray
    converged: bool
    n_iter: int

    def predict(self, x: np.ndarray | list[float] | pd.Series) -> np.ndarray:
        values = np.asarray(x, dtype=float)
        logits = np.clip(values * float(self.params[0]), -500.0, 500.0)
        return 1.0 / (1.0 + np.exp(-logits))


def fit_frequency_logistic_model(
    x: np.ndarray | pd.Series,
    y: np.ndarray | pd.Series,
    max_iter: int = 100,
    tol: float = 1e-8,
    ridge: float = 1e-6,
) -> LogisticFrequencyModel:
    """
    Fit a one-parameter logistic model with a small ridge penalty.

    The original implementation used statsmodels for `Logit(y, x)` with no
    intercept. This local implementation matches that functional form without
    requiring an additional dependency stack.
    """
    x_values = np.asarray(x, dtype=float).reshape(-1)
    y_values = np.asarray(y, dtype=float).reshape(-1)

    if x_values.size == 0:
        raise ValueError("Cannot fit logistic model on empty inputs")
    if x_values.shape != y_values.shape:
        raise ValueError("x and y must have the same shape")
    if not np.isfinite(x_values).all() or not np.isfinite(y_values).all():
        raise ValueError("x and y must contain only finite values")
    if np.unique(y_values).size < 2:
        raise ValueError("Need both positive and negative examples to fit logistic model")

    weight = 0.0
    converged = False
    n_iter = 0

    for n_iter in range(1, max_iter + 1):
        logits = np.clip(x_values * weight, -500.0, 500.0)
        probs = 1.0 / (1.0 + np.exp(-logits))

        gradient = float(np.sum(x_values * (y_values - probs)) - ridge * weight)
        hessian = float(-np.sum((x_values**2) * probs * (1.0 - probs)) - ridge)

        if abs(gradient) < tol:
            converged = True
            break

        if abs(hessian) < tol:
            logger.warning(
                "Logistic fit aborting at iter %d: near-singular Hessian (%.3g); "
                "weight=%.6g may be unconverged.",
                n_iter,
                hessian,
                weight,
            )
            break

        step = gradient / hessian
        weight -= step

        if abs(step) < tol:
            converged = True
            break

    return LogisticFrequencyModel(
        params=np.asarray([weight], dtype=float),
        converged=converged,
        n_iter=n_iter,
    )


def analytical_logistic_threshold(
    weight: float,
    fdr: float,
    min_value: float = 0.0,
    max_frequency: float = 1.0,
) -> float:
    """Exact frequency threshold x* solving sigmoid(weight * x*) = 1 - fdr.

    Returns ``float('inf')`` when ``weight <= 0``, ``fdr <= 0``, or the exact
    threshold exceeds ``max_frequency`` (so unreachable FDR tiers never clamp
    onto the largest observed clone).
    """
    if weight <= 0 or not np.isfinite(weight) or fdr <= 0:
        return float("inf")
    if fdr >= 1.0:
        return float(min_value)
    x_star = float(np.log((1.0 - fdr) / fdr) / weight)
    if x_star > max_frequency:
        return float("inf")
    return max(float(min_value), x_star)
