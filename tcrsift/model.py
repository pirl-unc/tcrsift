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

from __future__ import annotations

import logging
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

DEFAULT_FDRS = (0.15, 0.1, 0.01, 0.001, 0.0001)
DEFAULT_THRESHOLD_COLORS = (
    "grey",
    "pink",
    "orange",
    "red",
    "purple",
    "blue",
    "green",
    "brown",
    "cyan",
    "magenta",
)


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


def count_at_threshold(df: pd.DataFrame, freq_threshold: float) -> int:
    return int(
        (
            (df["specificity_description"] == "Single-culture")
            & (df["max_freq"] >= freq_threshold)
        ).sum()
    )


def calc_threshold(
    df: pd.DataFrame,
    x_plot: np.ndarray,
    y_plot: np.ndarray,
    fdr: float = 0.1,
    min_value: float = 0.0,
) -> tuple[float, int]:
    y_target = 1.0 - fdr
    threshold_idx = int(np.argmin(np.abs(y_plot - y_target)))
    freq_threshold = max(float(min_value), float(x_plot[threshold_idx]))
    n = count_at_threshold(df, freq_threshold)
    return freq_threshold, n


def calc_thresholds_and_counts(
    df: pd.DataFrame,
    fdrs: tuple[float, ...] | list[float] = DEFAULT_FDRS,
    min_freq_threshold: float = 0.09,
    default_freq_threshold: float = 0.5,
    only_avoid_viral: bool = True,
) -> tuple[dict[float, float], dict[float, int], LogisticFrequencyModel]:
    fdr_values = tuple(fdrs)
    fdr_to_threshold: dict[float, float] = {}
    threshold_to_count: dict[float, int] = {}

    target_above_min_freq = (df["max_freq"] > min_freq_threshold).values

    if only_avoid_viral:
        target = target_above_min_freq & (df["specificity_description"] != "Viral").values
    else:
        target = target_above_min_freq & (df["specificity_description"] == "Single-culture").values

    try:
        model = fit_frequency_logistic_model(df["max_freq"].values, target.astype(float))
    except ValueError:
        model = LogisticFrequencyModel(
            params=np.asarray([0.0], dtype=float),
            converged=False,
            n_iter=0,
        )
        weight = float(model.params[0])
    else:
        weight = float(model.params[0])

    if weight <= 0:
        for fdr in fdr_values:
            fdr_to_threshold[fdr] = float(default_freq_threshold)
        threshold_to_count[float(default_freq_threshold)] = count_at_threshold(
            df, float(default_freq_threshold)
        )
        return fdr_to_threshold, threshold_to_count, model

    x_plot = np.linspace(df["max_freq"].min(), df["max_freq"].max(), 10000)
    y_plot = model.predict(x_plot)
    for fdr in fdr_values:
        threshold, n = calc_threshold(
            df,
            x_plot,
            y_plot,
            fdr,
            min_value=min_freq_threshold,
        )
        fdr_to_threshold[fdr] = threshold
        threshold_to_count[threshold] = n

    return fdr_to_threshold, threshold_to_count, model


def annotate_plot_with_thresholds_and_counts(
    df: pd.DataFrame,
    ax: plt.Axes,
    model: LogisticFrequencyModel,
    fdr_to_threshold: dict[float, float],
    threshold_to_count: dict[float, int],
    preferred_fdr: float = 0.15,
    colors: tuple[str, ...] | list[str] = DEFAULT_THRESHOLD_COLORS,
) -> None:
    x_plot = np.linspace(df["max_freq"].min(), df["max_freq"].max(), 10000)
    y_plot = model.predict(x_plot)
    ax.plot(x_plot, y_plot, color="green", linewidth=1, linestyle=":")

    fdrs = reversed(sorted(fdr_to_threshold))
    for i, (fdr, color) in enumerate(zip(fdrs, colors)):
        threshold = fdr_to_threshold[fdr]
        n = threshold_to_count[threshold]
        ax.axvline(
            threshold,
            color=color,
            linestyle="--",
            alpha=0.6 if fdr == preferred_fdr else 0.4,
        )
        ax.text(
            x=threshold + 0.03,
            y=0.3 + i * 0.12,
            s=f"{100 * fdr:g}% FDR\nthreshold\n= {threshold:0.2}%\n(n={n})",
            fontweight="bold" if fdr == preferred_fdr else "medium",
            alpha=0.9 if fdr == preferred_fdr else 0.8,
        )
