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

"""Lightweight, zero-dependency TCRdist similarity metric and search index.

Implements the Dash et al. (*Nature* 2017) / Mayer-Blackwell et al. (*eLife* 2021)
TCR distance metric natively in pure Python + NumPy:

- **V-gene CDR1 / CDR2 / CDR2.5 ("pMHC") loop distance** (weight 1): fixed-width
  IMGT-aligned loop comparison using ``tcrsift/refseqs/vgene_cdr_loops.csv.gz``
  (covering human and mouse ``TRAV`` and ``TRBV`` genes).
- **BLOSUM62-weighted trimmed CDR3 distance** (weight 3, trimming conserved 3
  N-terminal and 2 C-terminal anchor residues): O(L) single-gap centered
  alignment with substitution cost ``min(4, max(0, 4 - BLOSUM62[a, b]))`` and
  gap penalty ``4`` (``12`` after the 3x CDR3 weight).
- **Exact admissible-bound search index** (:class:`TCRDistIndex`): prunes DB
  V-gene buckets by ``v_dist_a + v_dist_b <= radius`` and CDR3 length bounds
  before vectorized NumPy BLOSUM62 scoring — returning exact matches at radius
  ``<= max_distance`` without ``parasail`` or ``tcrdist3``.
"""

from __future__ import annotations

import functools
import gzip
from collections import defaultdict
from importlib.resources import files as _pkg_files

import numpy as np
import pandas as pd

from .genes import canonicalize_gene
from .species import normalize_species

# Standard 20 amino acids + gap/unknown sentinel (index 20).
_AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
_AA_TO_IDX: dict[str, int] = {aa: i for i, aa in enumerate(_AA_ORDER)}
_GAP_IDX = 20
_UNK_IDX = 21

# Dash et al. (2017) BLOSUM62 substitution distance: min(4, max(0, 4 - BLOSUM62[a, b])).
# Rows/columns follow _AA_ORDER = "ARNDCQEGHILKMFPSTWYV", plus gap (20) and unknown (21).
_BLOSUM62_DIAG_ROWS = (
    # A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    (4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0),    # A
    (-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3),    # R
    (-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3),        # N
    (-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3),   # D
    (0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1),# C
    (-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2),       # Q
    (-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2),      # E
    (0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3),  # G
    (-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3),    # H
    (-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3),   # I
    (-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1),   # L
    (-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2),    # K
    (-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1),    # M
    (-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1),    # F
    (-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2),# P
    (1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2),       # S
    (0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0),   # T
    (-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3),# W
    (-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1),  # Y
    (0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4),    # V
)

DEFAULT_GAP_PENALTY = 4
DEFAULT_CDR3_WEIGHT = 3
DEFAULT_NTRIM = 3
DEFAULT_CTRIM = 2
DEFAULT_PAIRED_RADIUS = 48
DEFAULT_BETA_RADIUS = 12


def _build_cost_matrix(gap_penalty: int = DEFAULT_GAP_PENALTY) -> np.ndarray:
    mat = np.full((22, 22), gap_penalty, dtype=np.int16)
    for i in range(20):
        for j in range(20):
            mat[i, j] = 0 if i == j else min(4, max(0, 4 - _BLOSUM62_DIAG_ROWS[i][j]))
    # Gap-to-gap is 0; gap-to-AA is gap_penalty; unknown-to-unknown is 0.
    mat[_GAP_IDX, _GAP_IDX] = 0
    mat[_UNK_IDX, _UNK_IDX] = 0
    return mat


_COST_MATRIX = _build_cost_matrix(DEFAULT_GAP_PENALTY)


def _encode_seq(seq: str) -> np.ndarray:
    arr = np.empty(len(seq), dtype=np.int8)
    for i, ch in enumerate(seq):
        if ch in (".", "-"):
            arr[i] = _GAP_IDX
        else:
            arr[i] = _AA_TO_IDX.get(ch, _UNK_IDX)
    return arr


def trim_cdr3(seq: str, ntrim: int = DEFAULT_NTRIM, ctrim: int = DEFAULT_CTRIM) -> str:
    """Trim conserved N- and C-terminal CDR3 anchor residues (default 3 and 2)."""
    if not isinstance(seq, str):
        return ""
    s = seq.strip().upper()
    if len(s) > ntrim + ctrim:
        return s[ntrim : len(s) - ctrim]
    return s


def cdr3_distance(
    seq1: str,
    seq2: str,
    *,
    weight: int = DEFAULT_CDR3_WEIGHT,
    gap_penalty: int = DEFAULT_GAP_PENALTY,
    ntrim: int = DEFAULT_NTRIM,
    ctrim: int = DEFAULT_CTRIM,
) -> int:
    """Compute Dash/TCRdist BLOSUM62 single-gap CDR3 distance between two CDR3 strings."""
    s1 = trim_cdr3(seq1, ntrim=ntrim, ctrim=ctrim)
    s2 = trim_cdr3(seq2, ntrim=ntrim, ctrim=ctrim)
    if s1 == s2:
        return 0
    a1 = _encode_seq(s1)
    a2 = _encode_seq(s2)
    return int(_cdr3_dist_encoded(a1, a2, gap_penalty=gap_penalty) * weight)


def _cdr3_dist_encoded(a1: np.ndarray, a2: np.ndarray, gap_penalty: int = DEFAULT_GAP_PENALTY) -> int:
    """Unweighted single-gap centered BLOSUM62 alignment distance between encoded trimmed CDR3s."""
    n1, n2 = len(a1), len(a2)
    if n1 == n2:
        if n1 == 0:
            return 0
        return int(_COST_MATRIX[a1, a2].sum())
    if n1 > n2:
        a1, a2 = a2, a1
        n1, n2 = n2, n1
    g = n2 - n1
    if n1 == 0:
        return g * gap_penalty
    # Prefix costs aligned to a2[:n1]; suffix costs aligned to a2[g:]
    pref_step = _COST_MATRIX[a1, a2[:n1]]
    suf_step = _COST_MATRIX[a1, a2[g:]]
    pref = np.empty(n1 + 1, dtype=np.int32)
    pref[0] = 0
    np.cumsum(pref_step, out=pref[1:])
    suf = np.empty(n1 + 1, dtype=np.int32)
    suf[n1] = 0
    suf[:n1] = np.cumsum(suf_step[::-1])[::-1]
    return int((pref + suf).min()) + g * gap_penalty


def _cdr3_dist_batch_same_len(q_arr: np.ndarray, ref_mat: np.ndarray) -> np.ndarray:
    """Vectorized unweighted BLOSUM62 distance from one query (L,) to ref_mat (M, L)."""
    if ref_mat.shape[1] == 0:
        return np.zeros(ref_mat.shape[0], dtype=np.int32)
    return _COST_MATRIX[ref_mat, q_arr].sum(axis=1, dtype=np.int32)


def _cdr3_dist_batch_diff_len(
    q_arr: np.ndarray, ref_mat: np.ndarray, gap_penalty: int = DEFAULT_GAP_PENALTY
) -> np.ndarray:
    """Vectorized unweighted single-gap BLOSUM62 distance from q_arr (Lq,) to ref_mat (M, Lr)."""
    lq = len(q_arr)
    lr = ref_mat.shape[1]
    if lq == lr:
        return _cdr3_dist_batch_same_len(q_arr, ref_mat)
    if lq < lr:
        g = lr - lq
        if lq == 0:
            return np.full(ref_mat.shape[0], g * gap_penalty, dtype=np.int32)
        pref_step = _COST_MATRIX[ref_mat[:, :lq], q_arr]
        suf_step = _COST_MATRIX[ref_mat[:, g:], q_arr]
        n = lq
    else:
        g = lq - lr
        if lr == 0:
            return np.full(ref_mat.shape[0], g * gap_penalty, dtype=np.int32)
        pref_step = _COST_MATRIX[ref_mat, q_arr[:lr]]
        suf_step = _COST_MATRIX[ref_mat, q_arr[g:]]
        n = lr
    m = ref_mat.shape[0]
    pref = np.zeros((m, n + 1), dtype=np.int32)
    np.cumsum(pref_step, axis=1, out=pref[:, 1:])
    suf = np.zeros((m, n + 1), dtype=np.int32)
    suf[:, :n] = np.cumsum(suf_step[:, ::-1], axis=1)[:, ::-1]
    return (pref + suf).min(axis=1) + g * gap_penalty


def _loop_dist(loop1: str, loop2: str, gap_penalty: int = DEFAULT_GAP_PENALTY) -> int:
    """Fixed-width IMGT loop distance (CDR1, CDR2, or CDR2.5)."""
    if loop1 == loop2:
        return 0
    n = min(len(loop1), len(loop2))
    a1 = _encode_seq(loop1[:n])
    a2 = _encode_seq(loop2[:n])
    d = int(_COST_MATRIX[a1, a2].sum())
    if len(loop1) != len(loop2):
        d += abs(len(loop1) - len(loop2)) * gap_penalty
    return d


@functools.lru_cache(maxsize=4)
def _load_v_loop_table(species: str = "human") -> dict[tuple[str, str], tuple[str, str, str]]:
    """Load ``{(chain, canonical_v_gene): (cdr1, cdr2, cdr25)}`` for ``species``."""
    sp = normalize_species(species)
    ref = _pkg_files("tcrsift.refseqs").joinpath("vgene_cdr_loops.csv.gz")
    with ref.open("rb") as raw, gzip.open(raw, "rt") as fh:
        df = pd.read_csv(fh, keep_default_na=False)
    sub = df[df["species"] == sp]
    out: dict[tuple[str, str], tuple[str, str, str]] = {}
    for row in sub.itertuples(index=False):
        out[(str(row.chain), str(row.gene))] = (str(row.cdr1), str(row.cdr2), str(row.cdr25))
    return out


@functools.lru_cache(maxsize=8)
def _pairwise_v_distances(species: str, chain: str) -> dict[str, dict[str, int]]:
    """Precompute exact pairwise V-gene CDR1+CDR2+CDR2.5 distance matrix for ``(species, chain)``."""
    loops = _load_v_loop_table(species)
    genes = sorted(g for (ch, g) in loops if ch == chain)
    dist: dict[str, dict[str, int]] = {g: {} for g in genes}
    for i, g1 in enumerate(genes):
        l1 = loops[(chain, g1)]
        dist[g1][g1] = 0
        for g2 in genes[i + 1 :]:
            l2 = loops[(chain, g2)]
            d = _loop_dist(l1[0], l2[0]) + _loop_dist(l1[1], l2[1]) + _loop_dist(l1[2], l2[2])
            dist[g1][g2] = d
            dist[g2][g1] = d
    return dist


def v_gene_distance(
    v1: str | None, v2: str | None, chain: str, *, species: str = "human"
) -> int | None:
    """IMGT CDR1+CDR2+CDR2.5 BLOSUM62 distance between two V genes, or None if unmapped."""
    g1 = canonicalize_gene(v1)
    g2 = canonicalize_gene(v2)
    if not g1 or not g2:
        return None
    if g1 == g2:
        return 0
    pw = _pairwise_v_distances(normalize_species(species), chain.lower())
    row = pw.get(g1)
    if row is None:
        return None
    return row.get(g2)


def tcr_distance(
    *,
    cdr3_a1: str | None = None,
    v_a1: str | None = None,
    cdr3_b1: str | None = None,
    v_b1: str | None = None,
    cdr3_a2: str | None = None,
    v_a2: str | None = None,
    cdr3_b2: str | None = None,
    v_b2: str | None = None,
    species: str = "human",
) -> int:
    """Compute TCRdist between two TCRs (paired αβ, β-only, or α-only).

    Missing V-gene arguments contribute 0 to the V-loop term (graceful CDR3-only
    degradation); mismatched/unrecognized V genes fall back to exact canonical
    name equality (0 if equal, 30 if unequal and unmapped).
    """
    total = 0
    for chain, c1, v1, c2, v2 in (
        ("alpha", cdr3_a1, v_a1, cdr3_a2, v_a2),
        ("beta", cdr3_b1, v_b1, cdr3_b2, v_b2),
    ):
        if c1 and c2:
            total += cdr3_distance(c1, c2)
        if v1 and v2:
            vd = v_gene_distance(v1, v2, chain, species=species)
            if vd is not None:
                total += vd
            elif canonicalize_gene(v1) != canonicalize_gene(v2):
                total += 30
    return total


class TCRDistIndex:
    """Exact radius-search index over a reference TCR database using admissible V+length bounds.

    Parameters
    ----------
    db_positions : np.ndarray
        Positional row indices into the reference DataFrame.
    cdr3_beta : Sequence[str]
        Reference CDR3β amino-acid strings.
    v_beta : Sequence[str | None], optional
        Reference TRBV gene calls.
    cdr3_alpha : Sequence[str | None], optional
        Reference CDR3α amino-acid strings (required for ``paired=True``).
    v_alpha : Sequence[str | None], optional
        Reference TRAV gene calls.
    species : str, default ``"human"``
        Organism species (``"human"`` or ``"mouse"``).
    paired : bool, default ``True``
        When True, indexes paired αβ entries (joint distance ``d_alpha + d_beta``);
        when False, indexes β-only entries (``d_beta``).
    """

    def __init__(
        self,
        db_positions: np.ndarray,
        cdr3_beta: list[str],
        *,
        v_beta: list[str | None] | None = None,
        cdr3_alpha: list[str | None] | None = None,
        v_alpha: list[str | None] | None = None,
        species: str = "human",
        paired: bool = True,
    ) -> None:
        self.species = normalize_species(species)
        self.paired = paired
        self._v_pw_a = _pairwise_v_distances(self.species, "alpha") if paired else {}
        self._v_pw_b = _pairwise_v_distances(self.species, "beta")
        self.has_v_alpha = paired and v_alpha is not None and any(
            canonicalize_gene(v) in self._v_pw_a for v in v_alpha
        )
        self.has_v_beta = v_beta is not None and any(
            canonicalize_gene(v) in self._v_pw_b for v in v_beta
        )

        # Group reference rows by (va_key, vb_key) -> (len_a, len_b) -> vectorized arrays.
        # DB rows lacking a mapped V gene on a chain are kept in the "" bucket so
        # mixed databases (e.g. VDJdb + V-less IEDB/CEDAR rows) never drop exact
        # CDR3 matches, while fuzzy matches still require mapped V-gene loops.
        raw_buckets: dict[
            tuple[str, str], dict[tuple[int, int], dict[tuple[str, str], list[int]]]
        ] = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        exact_trimmed: dict[tuple[str, str], list[int]] = defaultdict(list)

        for i, pos in enumerate(db_positions):
            b_raw = cdr3_beta[i]
            if not isinstance(b_raw, str) or not b_raw:
                continue
            tb = trim_cdr3(b_raw)
            vb_c = canonicalize_gene(v_beta[i]) if (v_beta is not None and self.has_v_beta) else ""
            if not vb_c or vb_c not in self._v_pw_b:
                vb_c = ""

            if paired:
                a_raw = cdr3_alpha[i] if cdr3_alpha is not None else None
                if not isinstance(a_raw, str) or not a_raw:
                    continue
                ta = trim_cdr3(a_raw)
                va_c = (
                    canonicalize_gene(v_alpha[i])
                    if (v_alpha is not None and self.has_v_alpha)
                    else ""
                )
                if not va_c or va_c not in self._v_pw_a:
                    va_c = ""
                raw_buckets[(va_c, vb_c)][(len(ta), len(tb))][(ta, tb)].append(int(pos))
                exact_trimmed[(ta, tb)].append(int(pos))
            else:
                raw_buckets[("", vb_c)][(0, len(tb))][("", tb)].append(int(pos))
                exact_trimmed[("", tb)].append(int(pos))

        # Freeze into compact NumPy matrices for fast vectorized query scoring.
        self._buckets: dict[
            tuple[str, str], list[tuple[int, int, np.ndarray | None, np.ndarray, list[np.ndarray]]]
        ] = {}
        self._v_keys_by_a: dict[str, list[str]] = defaultdict(list)
        for vkey, len_map in raw_buckets.items():
            va_k, vb_k = vkey
            self._v_keys_by_a[va_k].append(vb_k)
            blocks = []
            for (la, lb), seq_map in len_map.items():
                keys = list(seq_map.keys())
                m = len(keys)
                mat_b = np.empty((m, lb), dtype=np.int8)
                mat_a = np.empty((m, la), dtype=np.int8) if paired else None
                pos_groups: list[np.ndarray] = []
                for idx_k, (ta, tb) in enumerate(keys):
                    if paired and mat_a is not None:
                        mat_a[idx_k] = _encode_seq(ta)
                    mat_b[idx_k] = _encode_seq(tb)
                    pos_groups.append(np.asarray(seq_map[(ta, tb)], dtype=np.int64))
                blocks.append((la, lb, mat_a, mat_b, pos_groups))
            self._buckets[vkey] = blocks

        self._va_keys = list(self._v_keys_by_a.keys())
        self._exact_trimmed_pos: dict[tuple[str, str], np.ndarray] = {
            k: np.sort(np.asarray(v, dtype=np.int64)) for k, v in exact_trimmed.items()
        }

    def _candidate_v_buckets(
        self,
        v_alpha: str | None,
        v_beta: str | None,
        radius: int,
    ) -> tuple[list[tuple[int, str, str, bool, bool]], bool]:
        """Resolve admissible ``(v_dist, va_k, vb_k, exact_a_only, exact_b_only)`` buckets.

        Returns ``(buckets, needs_exact_fallback)`` where ``needs_exact_fallback`` is
        True if a V-annotated chain has a missing/unmapped query V gene (in which
        case exact trimmed-CDR3 lookup preserves exact hits without fuzzy false
        positives).
        """
        if self.has_v_alpha:
            qva = canonicalize_gene(v_alpha)
            va_row = self._v_pw_a.get(qva) if qva else None
            if va_row is None:
                return [], True
            va_cands: list[tuple[str, int, bool]] = []
            for va_k in self._va_keys:
                if va_k == "":
                    va_cands.append(("", 0, True))
                else:
                    d_va = va_row.get(va_k)
                    if d_va is not None and d_va <= radius:
                        va_cands.append((va_k, d_va, False))
        else:
            va_cands = [("", 0, False)] if "" in self._v_keys_by_a else []

        if self.has_v_beta:
            qvb = canonicalize_gene(v_beta)
            vb_row = self._v_pw_b.get(qvb) if qvb else None
            if vb_row is None:
                return [], True
        else:
            vb_row = None

        bucket_cands: list[tuple[int, str, str, bool, bool]] = []
        for va_k, dva, exact_a in va_cands:
            for vb_k in self._v_keys_by_a[va_k]:
                if vb_row is not None:
                    if vb_k == "":
                        dvb, exact_b = 0, True
                    else:
                        d_vb = vb_row.get(vb_k)
                        if d_vb is None:
                            continue
                        dvb, exact_b = d_vb, False
                else:
                    dvb, exact_b = 0, False
                vd = dva + dvb
                if vd <= radius:
                    bucket_cands.append((vd, va_k, vb_k, exact_a, exact_b))
        bucket_cands.sort(key=lambda t: t[0])
        return bucket_cands, False

    def query_nearest(
        self,
        cdr3_beta: str,
        *,
        v_beta: str | None = None,
        cdr3_alpha: str | None = None,
        v_alpha: str | None = None,
        max_distance: int | None = None,
    ) -> tuple[int, np.ndarray] | None:
        """Find all nearest-distance DB row positions within ``max_distance``.

        Returns ``(best_distance, np.ndarray_of_db_positions)`` or ``None`` if no
        reference entry lies within ``max_distance``.
        """
        if not cdr3_beta or (self.paired and not cdr3_alpha):
            return None
        radius = (
            (DEFAULT_PAIRED_RADIUS if self.paired else DEFAULT_BETA_RADIUS)
            if max_distance is None
            else int(max_distance)
        )
        tb = trim_cdr3(cdr3_beta)
        qb = _encode_seq(tb)
        lqb = len(qb)

        if self.paired:
            ta = trim_cdr3(cdr3_alpha or "")
            qa = _encode_seq(ta)
            lqa = len(qa)
        else:
            ta = ""
            qa = None
            lqa = 0

        bucket_cands, needs_exact_fallback = self._candidate_v_buckets(v_alpha, v_beta, radius)
        if needs_exact_fallback:
            exact_pos = self._exact_trimmed_pos.get((ta, tb))
            return (0, exact_pos) if exact_pos is not None else None

        best_d = radius + 1
        best_pos: list[np.ndarray] = []
        gap_unit = DEFAULT_GAP_PENALTY * DEFAULT_CDR3_WEIGHT  # 12

        for vd, va_k, vb_k, exact_a, exact_b in bucket_cands:
            if vd > best_d:
                break
            for la, lb, mat_a, mat_b, pos_groups in self._buckets[(va_k, vb_k)]:
                if (exact_b and lb != lqb) or (self.paired and exact_a and la != lqa):
                    continue
                len_lb = (abs(lqb - lb) + (abs(lqa - la) if self.paired else 0)) * gap_unit
                if vd + len_lb > best_d:
                    continue
                rem_cdr3_unweighted = (best_d - vd) // DEFAULT_CDR3_WEIGHT
                db_raw = _cdr3_dist_batch_diff_len(qb, mat_b)
                if not self.paired:
                    if exact_b:
                        keep = np.flatnonzero(db_raw == 0)
                        if len(keep) == 0:
                            continue
                        tot = np.full(len(keep), vd, dtype=np.int32)
                        hits = keep
                        min_block = vd
                    else:
                        tot = vd + db_raw * DEFAULT_CDR3_WEIGHT
                        min_block = int(tot.min())
                        if min_block > best_d:
                            continue
                        hits = np.flatnonzero(tot == min_block)
                else:
                    # Prune rows where beta alone already exceeds remaining budget
                    b_limit = 0 if exact_b else rem_cdr3_unweighted
                    keep = np.flatnonzero(db_raw <= b_limit)
                    if len(keep) == 0:
                        continue
                    assert qa is not None and mat_a is not None
                    da_raw = _cdr3_dist_batch_diff_len(qa, mat_a[keep])
                    if exact_a:
                        ok_a = np.flatnonzero(da_raw == 0)
                        if len(ok_a) == 0:
                            continue
                        keep = keep[ok_a]
                        da_raw = da_raw[ok_a]
                    tot = vd + (da_raw + db_raw[keep]) * DEFAULT_CDR3_WEIGHT
                    min_block = int(tot.min())
                    if min_block > best_d:
                        continue
                    hits = keep[np.flatnonzero(tot == min_block)]

                if min_block < best_d:
                    best_d = min_block
                    best_pos = [pos_groups[h] for h in hits]
                else:
                    best_pos.extend(pos_groups[h] for h in hits)

        if not best_pos or best_d > radius:
            return None
        return best_d, np.sort(np.concatenate(best_pos))
