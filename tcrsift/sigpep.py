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

"""Signal-peptide feature analysis — von Heijne tripartite model (#270).

Scores a leader against the canonical signal-peptide features so a contig SP
that diverges from germline (e.g. an internal-deletion that may be a rare
germline indel) can be kept on the evidence that it still has *everything
expected of a signal peptide* — rather than substituted just because its length
differs from the reference. The defining feature is the hydrophobic h-region;
the n-region charge, c-region cleavage context, and helix-breaker are recorded
as corroborating detail.

Tripartite model:
  - **n-region** (1-5 aa): begins with the initiator Met; net positive
    (Lys/Arg) — the "positive-inside" rule.
  - **h-region / hydrophobic core** (~7-15 aa): a strongly hydrophobic
    (L/A/V/I/F/M) stretch — SRP-recognized, the single most critical element.
    Must clear a hydrophobicity threshold over a >=7-residue window, no charged
    residue inside.
  - **c-region** (3-7 aa): more polar, carries the cleavage site; small/uncharged
    at -3 and -1 (von Heijne (-3,-1) rule), Pro disfavoured at -3..-1; a
    helix-breaker (Pro/Gly) near the h/c boundary positions cleavage.
  - Overall ~15-30 aa, start Met.

``signal_peptide_features`` returns the measured features plus a ``features_ok``
gate (the defining subset: Met start, sane length, a strong >=7-aa h-core with
no internal charge) and a ``score`` (fraction of all criteria met, including the
softer cleavage / charge / helix-breaker preferences that real TCR leaders don't
always satisfy — e.g. -1=L in TRBV20-1).
"""

from __future__ import annotations

# Kyte-Doolittle hydropathy (higher = more hydrophobic).
_KD = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
_CHARGED = frozenset("DEKR")
_SMALL = frozenset("AGSCT")  # accepted at cleavage -3 / -1 (von Heijne)
_POSITIVE = frozenset("KR")
_NEGATIVE = frozenset("DE")

# h-core acceptance: a >=7-residue window with mean Kyte-Doolittle above this.
H_CORE_MIN_LEN = 7
H_CORE_MIN_KD = 1.6
SP_LEN_MIN = 12
SP_LEN_MAX = 30


def _mean_kd(seq: str) -> float:
    return sum(_KD.get(a, 0.0) for a in seq) / len(seq) if seq else 0.0


def _find_h_core(leader: str) -> tuple[int, int, float] | None:
    """Most hydrophobic contiguous window (length 7-15) → (start, end, mean_KD),
    or None when the leader is too short to contain a 7-aa core."""
    n = len(leader)
    if n < H_CORE_MIN_LEN:
        return None
    best = None
    for i in range(n - H_CORE_MIN_LEN + 1):
        for j in range(i + H_CORE_MIN_LEN, min(i + 16, n) + 1):
            m = _mean_kd(leader[i:j])
            if best is None or m > best[2]:
                best = (i, j, m)
    return best


def signal_peptide_features(leader: str | None) -> dict:
    """Measure von Heijne SP features for a leader. See module docstring.

    Always returns the full key set (``features_ok`` False / fields None for a
    missing or too-short leader) so callers can emit a uniform schema.
    """
    keys_ok = (
        "met_start", "length_ok", "h_core_ok", "h_core_no_charge",
        "n_positive", "cleavage_ok", "helix_breaker",
    )
    if not isinstance(leader, str) or not leader:
        return {
            "length": 0, "h_core": "", "h_core_len": 0, "h_core_kd": 0.0,
            "n_region": "", "n_charge": 0, "c_region": "",
            "cleavage_m3": "", "cleavage_m1": "",
            **{k: False for k in keys_ok},
            "features_ok": False, "score": 0.0,
        }

    n = len(leader)
    met_start = leader.startswith("M")
    length_ok = SP_LEN_MIN <= n <= SP_LEN_MAX

    core = _find_h_core(leader)
    if core is None:
        h_start, h_end, h_kd, h_core = n, n, 0.0, ""
    else:
        h_start, h_end, h_kd = core
        h_core = leader[h_start:h_end]
    h_core_ok = bool(core) and (h_end - h_start) >= H_CORE_MIN_LEN and h_kd >= H_CORE_MIN_KD
    h_core_no_charge = not (set(h_core) & _CHARGED)

    n_region = leader[:h_start]
    n_charge = sum(a in _POSITIVE for a in n_region) - sum(a in _NEGATIVE for a in n_region)
    n_positive = n_charge >= 0  # positive-inside (Met-only n-region scores 0 = ok)

    c_region = leader[h_end:]
    m1 = leader[-1] if n >= 1 else ""
    m3 = leader[-3] if n >= 3 else ""
    cleavage_ok = (m1 in _SMALL) and (m3 in _SMALL) and ("P" not in leader[-3:])
    # Helix-breaker (Pro/Gly) at/after the h/c boundary positions the cleavage.
    helix_breaker = bool(set(leader[max(h_end - 1, 0):]) & set("PG"))

    gate = met_start and length_ok and h_core_ok and h_core_no_charge
    crit = (met_start, length_ok, h_core_ok, h_core_no_charge,
            n_positive, cleavage_ok, helix_breaker)
    return {
        "length": n,
        "h_core": h_core, "h_core_len": len(h_core), "h_core_kd": round(h_kd, 2),
        "n_region": n_region, "n_charge": n_charge,
        "c_region": c_region, "cleavage_m3": m3, "cleavage_m1": m1,
        "met_start": met_start, "length_ok": length_ok,
        "h_core_ok": h_core_ok, "h_core_no_charge": h_core_no_charge,
        "n_positive": n_positive, "cleavage_ok": cleavage_ok,
        "helix_breaker": helix_breaker,
        "features_ok": gate,
        "score": round(sum(crit) / len(crit), 3),
    }


def sp_features_summary(f: dict) -> str:
    """Compact, auditable one-line summary of the SP features."""
    if not f.get("length"):
        return "no_leader"
    return (
        f"M{'' if f['met_start'] else '!'} len={f['length']} "
        f"n[{f['n_region']} {f['n_charge']:+d}] "
        f"h[{f['h_core']} {f['h_core_len']}aa KD{f['h_core_kd']:+.1f}] "
        f"c[{f['c_region']}] cleav[-3={f['cleavage_m3']} -1={f['cleavage_m1']}"
        f"{':ok' if f['cleavage_ok'] else ''}]"
        f"{' helixbrk' if f['helix_breaker'] else ''}"
    )
