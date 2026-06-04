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

"""Per-condition signature-consistency QC (#161).

Does each enrichment sort behave like its label? The single most diagnostic
cohort signal: each sort should be enriched for the GEX signature its gate
implies, consistently across donors. A sort that *isn't* — or a donor whose
whole sort×signature pattern is an outlier — is a failed/mislabelled sort or
a weird donor, flagged at a glance instead of via a forensic dig.

- **Within-donor:** each sort's expected signature(s) are enriched
  (z-scored across the donor's cells) vs that donor's baseline.
- **Cross-donor:** a donor whose per-(sort, signature) enrichment matrix
  correlates poorly with the leave-one-out cohort consensus is flagged.

The sort→signature map is **config-driven** (the biology isn't hardcoded).
NB: ``CTYneg`` is a **Cell Trace Yellow dilution** gate — it selects
*divided / expanded* cells, so its signature is Proliferation, **not** low
cytolytic. Compound sorts (``AIMpos_CTYneg``) require all components.
"""

from __future__ import annotations

import logging
import re

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Default sort → [(signature, direction)] map. Direction is "high" (enriched)
# or "low" (depleted). Sorts not listed (e.g. tetpos — a tetramer/specificity
# gate with no clean GEX correlate) are skipped for the within-donor check.
DEFAULT_SORT_SIGNATURE_MAP: dict[str, list[tuple[str, str]]] = {
    "AIMpos": [("AcuteActivation", "high")],        # AIM markers (4-1BB/OX40)
    "CTYneg": [("Proliferation", "high")],          # CTY dilution = divided/expanded
    "IFNpos": [("AntigenExperienced", "high")],     # IFN-γ secretion = effector
}

# strip a trailing donor suffix like "-2" / "_B1-2" is left to the caller via
# donor_col; here we only strip a trailing "-<digits>".
_DONOR_SUFFIX = re.compile(r"-\d+$")


def parse_sort_label(method: str) -> list[str]:
    """Split a sort label into its base component sorts.

    ``"AIMpos_CTYneg-2"`` → ``["AIMpos", "CTYneg"]`` (donor suffix stripped,
    compound split on ``_``).
    """
    if not isinstance(method, str) or not method:
        return []
    base = _DONOR_SUFFIX.sub("", method.strip())
    return [tok for tok in base.split("_") if tok]


def _base_method(method: str) -> str:
    """Strip a trailing donor suffix so a sort aligns across donors.

    ``"AIMpos_CTYneg-2"`` → ``"AIMpos_CTYneg"`` (keeps the compound label).
    """
    return _DONOR_SUFFIX.sub("", str(method).strip())


def sort_signature_consistency(
    df: pd.DataFrame,
    *,
    method_col: str = "enrichment_method",
    donor_col: str = "donor",
    signature_cols: list[str] | None = None,
    signature_prefix: str = "signature_",
    sort_signature_map: dict[str, list[tuple[str, str]]] | None = None,
    zscore_within_donor: bool = True,
    enrich_threshold: float = 0.0,
    cross_donor_min_corr: float = 0.5,
) -> pd.DataFrame:
    """Per-(donor, sort) signature-consistency QC table + warnings (#161).

    ``df`` is a per-cell (or per-clone) frame with ``method_col``,
    ``donor_col`` and per-signature score columns (``signature_<name>`` by
    default; signature names match :data:`DEFAULT_SORT_SIGNATURE_MAP` values
    without the prefix). Scores are z-scored within each donor (so enrichment
    is relative to that donor's baseline) and averaged per sort.

    Returns one row per (donor, sort) with: ``n_cells``, the mean z of each
    expected signature, ``within_donor_consistent`` (all expected signatures
    pass their direction at ``enrich_threshold``), ``donor_pattern_corr``
    (correlation of the donor's sort×signature matrix to the leave-one-out
    cohort consensus; NaN with < 3 donors) + ``donor_outlier``, and a
    ``warning`` string. Sorts with no mapped signature are reported (n_cells)
    but not consistency-flagged.
    """
    sort_map = sort_signature_map or DEFAULT_SORT_SIGNATURE_MAP
    if method_col not in df.columns:
        raise ValueError(f"sort_signature_consistency: missing {method_col!r}")
    if signature_cols is None:
        signature_cols = [c for c in df.columns if c.startswith(signature_prefix)]
    if not signature_cols:
        raise ValueError(
            "sort_signature_consistency: no signature score columns "
            f"(prefix {signature_prefix!r}) found")
    sig_name = {c: c[len(signature_prefix):] if c.startswith(signature_prefix)
                else c for c in signature_cols}
    have_donor = donor_col in df.columns
    work = df.copy()
    if not have_donor:
        work[donor_col] = "(all)"

    # z-score signatures within each donor (transform avoids the grouping-
    # column apply deprecation), then mean per (donor, sort).
    if zscore_within_donor:
        for c in signature_cols:
            work[c] = work.groupby(donor_col, observed=True)[c].transform(
                lambda v: (v - v.mean()) / v.std(ddof=0)
                if v.std(ddof=0) and v.std(ddof=0) > 0 else 0.0)
    agg = (
        work.groupby([donor_col, method_col], observed=True)[signature_cols]
        .mean().reset_index()
    )
    agg["n_cells"] = (
        work.groupby([donor_col, method_col], observed=True).size()
        .reset_index(name="n_cells")["n_cells"]
    )
    # Base sort label (donor suffix stripped) so sorts align across donors.
    agg["_base"] = agg[method_col].map(_base_method)

    # cross-donor: per-donor flattened (sort×signature) vector vs LOO consensus.
    pivot = agg.set_index([donor_col, "_base"])[signature_cols]
    pivot = pivot[~pivot.index.duplicated(keep="first")]
    donors = list(agg[donor_col].unique())
    corr = {}
    if len(donors) >= 3:
        for d in donors:
            mine = pivot.xs(d, level=donor_col)
            others = pivot.drop(index=d, level=donor_col)
            consensus = others.groupby(level="_base", observed=True).mean()
            common = mine.index.intersection(consensus.index)
            if len(common) >= 2:
                a = mine.loc[common].to_numpy().ravel()
                b = consensus.loc[common].to_numpy().ravel()
                mask = np.isfinite(a) & np.isfinite(b)
                corr[d] = (float(np.corrcoef(a[mask], b[mask])[0, 1])
                           if mask.sum() >= 2 and np.std(a[mask]) > 0
                           and np.std(b[mask]) > 0 else np.nan)
            else:
                corr[d] = np.nan
    else:
        corr = {d: np.nan for d in donors}

    rows = []
    for _, r in agg.iterrows():
        donor, method = r[donor_col], r[method_col]
        expected = []
        for token in parse_sort_label(str(method)):
            expected.extend(sort_map.get(token, []))
        row = {donor_col: donor, method_col: method, "n_cells": int(r["n_cells"])}
        problems = []
        for sig, direction in expected:
            col = next((c for c in signature_cols if sig_name[c] == sig), None)
            if col is None:
                continue
            z = float(r[col])
            row[f"z_{sig}"] = z
            ok = z >= enrich_threshold if direction == "high" else z <= -enrich_threshold
            if not ok:
                problems.append(
                    f"{sig} not {direction} (z={z:+.2f})")
        row["expected_signatures"] = ", ".join(f"{s}:{d}" for s, d in expected)
        row["within_donor_consistent"] = (not problems) if expected else True
        dc = corr.get(donor, np.nan)
        row["donor_pattern_corr"] = dc
        row["donor_outlier"] = bool(dc == dc and dc < cross_donor_min_corr)
        warn = []
        if expected and problems:
            warn.append("sort not enriched as labelled: " + "; ".join(problems))
        if row["donor_outlier"]:
            warn.append(f"donor sort-pattern outlier (consensus corr {dc:.2f})")
        row["warning"] = "; ".join(warn)
        rows.append(row)
    return pd.DataFrame(rows)


def sort_signature_consistency_from_adata(
    adata,
    *,
    method_col: str = "enrichment_method",
    donor_col: str = "donor",
    sort_signature_map: dict[str, list[tuple[str, str]]] | None = None,
    **kwargs,
) -> pd.DataFrame:
    """Per-cell signature scores from an AnnData → :func:`sort_signature_consistency`.

    Scores each signature named in the (default) sort→signature map per cell
    (via the shared expression resolver, ignoring genes absent from the
    matrix), assembles a per-cell frame with ``method_col`` + ``donor_col``,
    and runs the consistency QC. Returns an empty frame when the method
    column is absent or no mapped signature's genes are measurable.
    """
    from .signature_methods import SIGNATURES, score_signature

    sort_map = sort_signature_map or DEFAULT_SORT_SIGNATURE_MAP
    obs = adata.obs
    if method_col not in obs.columns:
        return pd.DataFrame()
    donor = donor_col if donor_col in obs.columns else (
        "patient_id" if "patient_id" in obs.columns else None)
    needed = {sig for pairs in sort_map.values() for sig, _ in pairs}
    out = pd.DataFrame(index=obs.index)
    out[method_col] = obs[method_col].to_numpy()
    if donor is not None:
        out["donor"] = obs[donor].to_numpy()
    scored = 0
    for name in needed:
        sig = SIGNATURES.get(name)
        if sig is None:
            continue
        from .signature_methods import expression_frame_from_adata
        expr = expression_frame_from_adata(adata, sig.all_genes, on_missing="ignore")
        if expr.shape[1] == 0:
            continue
        out[f"signature_{name}"] = score_signature(
            expr, sig, on_missing="ignore").to_numpy()
        scored += 1
    if scored == 0:
        return pd.DataFrame()
    return sort_signature_consistency(
        out, method_col=method_col, donor_col="donor",
        sort_signature_map=sort_map, **kwargs)
