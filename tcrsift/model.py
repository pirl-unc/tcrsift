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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm


def count_at_threshold(df : pd.DataFrame, freq_threshold : float) -> int:
    return  (
        (df["specificity_description"] == "Single-culture") & (df["max_freq"] >= freq_threshold)
    ).sum()

def calc_threshold(df : pd.DataFrame, x_plot : np.ndarray, y_plot : np.ndarray, fdr : float=0.1, min_value : float=0.0):
    y_target = 1.0 - fdr

    threshold_idx = np.argmin(np.abs(y_plot - y_target))
    freq_threshold = max(min_value, x_plot[threshold_idx])
    n = count_at_threshold(df, freq_threshold)
    return freq_threshold, n

def calc_thresholds_and_counts(
        df : pd.DataFrame,
        fdrs : list[float] = [0.15, 0.1, 0.01, 0.001,  0.0001],
        min_freq_threshold : float = 0.09,
        default_freq_threshold : float = 0.5,
        only_avoid_viral : bool =True) -> tuple[dict[float, float], dict[float, int], sm.Logit]:
    fdr_to_threshold = {}
    threshold_to_count = {}

     # Fit logistic regression
    target_above_min_freq = (df["max_freq"] > min_freq_threshold).values

    if only_avoid_viral:
        target = target_above_min_freq & (df["specificity_description"] != "Viral").values
    else:
        target = target_above_min_freq & (df["specificity_description"] == "Single-culture").values

    model = sm.Logit(target, df["max_freq"].values).fit()
    weight = model.params[0]
    if weight < 0:
        print(f"-- Data too messy to choose thresholds adaptively, defaulting to {default_freq_threshold:0.2f}")
        for fdr in fdrs:
            fdr_to_threshold[fdr] = default_freq_threshold
            threshold_to_count[default_freq_threshold] = count_at_threshold(df, default_freq_threshold)
    else:
        x_plot = np.linspace(df["max_freq"].min(), df["max_freq"].max(), 10000)
        y_plot = model.predict(x_plot)
        for fdr in fdrs:
            t, n = calc_threshold(df, x_plot, y_plot, fdr, min_value=min_freq_threshold)
            fdr_to_threshold[fdr] = t
            threshold_to_count[t] = n
    return fdr_to_threshold, threshold_to_count, model

def annotate_plot_with_thresholds_and_counts(
        df : pd.DataFrame,
        ax : plt.Axes,
        model : sm.Logit,
        fdr_to_threshold : dict[float, float],
        threshold_to_count : dict[float, int],
        preferred_fdr : float = 0.15,
        colors=['grey', 'pink', 'orange', 'red', 'purple', 'blue', 'green', 'brown',  'cyan', 'magenta']):
    for i, (fdr, color) in enumerate(zip(fdr_to_threshold.keys(), colors)):
        x_plot = np.linspace(df["max_freq"].min(), df["max_freq"].max(), 10000)
        y_plot = model.predict(x_plot)
        ax.plot(x_plot, y_plot, color="green", linewidth=1, linestyle=':')
        fdrs = reversed(sorted(fdr_to_threshold.keys()))
        for i, (fdr, color) in enumerate(zip(fdrs, colors)):
            t = fdr_to_threshold[fdr]
            n = threshold_to_count[t]
            ax.axvline(t, color=color, linestyle="--", alpha=0.6 if fdr == preferred_fdr else 0.4)
            ax.text(
                x=t+0.03,
                y=0.3 + i * 0.12,
                s=f"{100*fdr:g}% FDR\nthreshold\n= {t:0.2}%\n(n={n})",
                fontweight='bold' if fdr==preferred_fdr else 'medium',
                alpha=0.9 if fdr==preferred_fdr else 0.8)
