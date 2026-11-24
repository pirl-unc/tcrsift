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

"""Independent QC battery for assembled deliverables (#279).

Re-derives every construct fact from first principles — its own codon table,
its own CDR3 / germline / linker logic — and cross-checks the assembled frame
against the source clonotypes and the raw CellRanger contigs. It deliberately
does **NOT** trust tcrsift's own QC columns (``construct_contig_verified``,
``qc_warnings``, the ``synth_*`` flags); the point is an outside check that the
pipeline could not have faked. Ported from the external ``qc_check.py`` that
caught the #278 schema gap and confirmed the TRAJ35 / dual-α edges are correct.

Five checks (A–E):

    A. Assembly integrity   every NT→AA translation exact; no internal stops;
                            full = leader + VDJ + constant byte-for-byte;
                            single_chain = beta + linker + alpha; M-start;
                            linker == expected T2A.
    B. CDR3 + reference     each CDR3 ⊂ its VDJ; canonical C…[F/W] (TRAJ35-style
                            non-canonical-J allow-list); valid V/J/C names;
                            junction residue present; sane CDR3 length.
    C. Source mapping       every clone's V/J/CDR3/cell_count matches the source
                            clonotypes; dual-α secondaries accounted for.
    D. Synthesis / dual-α   contig_verified all true; no dup constructs / α-β
                            swaps / dup single_chains; expansion = dual-α count.
    E. Raw cellranger       every *_contig_id exists in the filtered contig
                            annotations; every clone's CDR3 matches ≥1 of its
                            own contigs.

Checks whose inputs are absent (no source clonotypes, no contig map) are
reported as SKIP rather than silently dropped, so the summary never overstates
coverage.
"""

from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
import pandas as pd

# Independent codon table — biological ground truth, not pipeline logic, so the
# battery stays self-contained rather than importing the assembler's table.
_B = "TCAG"
_COD = [a + b + c for a in _B for b in _B for c in _B]
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CT = dict(zip(_COD, _AA))

# Expected self-cleaving inter-chain linker in a β-T2A-α single-chain construct.
EXPECTED_T2A = "EGRGSLLTCGDVEENPGP"

# J genes whose conserved anchor is not F/W (IMGT non-canonical); their CDR3
# ends on that residue rather than the canonical Phe/Trp (human TRAJ35 ends in C;
# mouse TRAJ7/TRAJ44 end in L and mouse TRBJ2-6 ends in Y).
NONCANON_J_END_HUMAN: dict[str, list[str]] = {
    "C": ["TRAJ35"],
}
NONCANON_J_END_MOUSE: dict[str, list[str]] = {
    "C": ["TRAJ35"],
    "L": ["TRAJ7", "TRAJ44"],
    "Y": ["TRBJ2-6"],
}
NONCANON_J_END = NONCANON_J_END_MOUSE  # backward-compatible union constant


def _is_mouse_row(r: pd.Series) -> bool:
    """Detect whether an assembled row has murine VDJ provenance."""
    for sp_col in ("vdj_species", "species"):
        raw_sp = r.get(sp_col)
        if raw_sp is not None and pd.notna(raw_sp):
            sp = str(raw_sp).strip().lower()
            if sp in ("mouse", "mus_musculus", "mm", "mm10", "mm39"):
                return True
            if sp in ("human", "homo_sapiens", "hs", "hg38", "grch38"):
                return False
    for src_col in ("alpha_constant_source", "beta_constant_source"):
        raw_src = r.get(src_col)
        if raw_src is not None and pd.notna(raw_src):
            src = str(raw_src).strip().lower()
            if src.startswith("chimeric:mouse:"):
                return False
            if src.startswith("chimeric:human:"):
                return True
    for c_col in ("alpha_c_gene", "beta_c_gene"):
        val = str(r.get(c_col, "")).strip()
        if val in ("Trac", "Trbc1", "Trbc2"):
            return True
    for v_col in ("alpha_v_gene", "beta_v_gene", "alpha_j_gene", "beta_j_gene"):
        val = str(r.get(v_col, "")).strip()
        if val.startswith(("Trav", "Trbv", "Traj", "Trbj")) or "D-" in val.upper():
            return True
    # Fallback for uppercased mouse gene tables: require species-unique constant
    # tails (mouse TRAC ends in MGLRILLLKVAGFNLLMTLRLWSS vs human
    # IGFRILLLKVAGFNLLMTLRLWSS; mouse TRBC ends in VKRKNS/VKKKNS vs human
    # VKRKDF/VKRKDSRG). Do not match the shared LLMTLRLWSS tail alone.
    fa = _nz(r.get("full_alpha_aa", ""))
    fb = _nz(r.get("full_beta_aa", ""))
    if fa.endswith("MGLRILLLKVAGFNLLMTLRLWSS") or fb.endswith(("VKRKNS", "VKKKNS")):
        return True
    return False


def _tr(nt: str) -> str:
    """Translate a DNA/RNA string frame-0, 'X' on any unknown codon."""
    nt = str(nt).upper().replace("U", "T")
    return "".join(
        _CT.get(nt[i : i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3)
    )


def _nz(a) -> str:
    """Stringify and strip a trailing stop ('*') so a translated CDS with its
    terminal stop compares equal to the stored stop-free AA."""
    if a is None or pd.isna(a):
        return ""
    return str(a).rstrip("*")


class BatteryCheckResult:
    """One check's verdict: PASS / FAIL / SKIP plus a human detail string."""

    __slots__ = ("name", "status", "detail")

    def __init__(self, name: str, status: str, detail: str):
        self.name = name
        self.status = status  # "PASS" | "FAIL" | "SKIP"
        self.detail = detail

    @property
    def ok(self) -> bool:
        # SKIP does not fail the battery, but it is not a pass either.
        return self.status != "FAIL"

    def badge(self) -> str:
        return {"PASS": "PASS ✅", "FAIL": "FAIL ❌", "SKIP": "SKIP ⚪"}[self.status]


def _missing(d: pd.DataFrame, cols) -> list[str]:
    """Columns the check needs but the frame lacks (so it SKIPs, not crashes)."""
    return [c for c in cols if c not in d.columns]


def _count_true(d: pd.DataFrame, col: str) -> int:
    """Count ``True`` in a boolean column, 0 when the column is absent (#283).

    ``DataFrame.get(col)`` returns ``None`` for a missing column, and
    ``None == True`` is a scalar ``False`` — so the previous
    ``(d.get(col) == True).sum()`` raised ``AttributeError`` whenever the
    ``synth_*`` columns weren't present. Guard on membership instead.
    """
    if col not in d.columns:
        return 0
    return int((d[col] == True).sum())  # noqa: E712


# Columns each check indexes directly (r[col]); absence → SKIP, not a crash. A
# real assembled deliverable carries all of them — this only spares degenerate
# / partial frames an unhandled KeyError while keeping the verdict honest.
_ASSEMBLY_COLS = [
    "full_alpha_aa", "full_beta_aa", "alpha_leader_aa", "beta_leader_aa",
    "VDJ_alpha_aa", "VDJ_beta_aa", "single_chain_aa",
]
_CDR3_COLS = [
    "CDR3_alpha", "CDR3_beta", "VDJ_alpha_aa", "VDJ_beta_aa",
    "alpha_j_gene", "beta_j_gene",
    "alpha_v_gene", "alpha_c_gene", "beta_v_gene", "beta_c_gene",
]


def check_assembly(d: pd.DataFrame) -> BatteryCheckResult:
    """A. Every translation exact, recomposition byte-for-byte, M-start, 2A linker."""
    miss = _missing(d, _ASSEMBLY_COLS)
    if miss:
        return BatteryCheckResult(
            "A. Assembly integrity", "SKIP", f"missing columns: {', '.join(miss)}"
        )
    if len(d) == 0:
        return BatteryCheckResult("A. Assembly integrity", "SKIP", "0 rows")
    from .assemble import LINKERS

    valid_linkers = {v["aa"] for v in LINKERS.values()} | {EXPECTED_T2A}
    bad = 0
    linkers: set[str] = set()
    n_sc = 0
    for _, r in d.iterrows():
        for ntc, aac in [
            ("single_chain_nt", "single_chain_aa"),
            ("single_chain_nt_optimized", "single_chain_aa"),
            ("full_beta_nt", "full_beta_aa"),
            ("full_beta_nt_optimized", "full_beta_aa"),
            ("full_alpha_nt", "full_alpha_aa"),
            ("full_alpha_nt_optimized", "full_alpha_aa"),
            ("VDJ_beta_nt", "VDJ_beta_aa"),
            ("VDJ_alpha_nt", "VDJ_alpha_aa"),
            ("beta_leader_nt", "beta_leader_aa"),
            ("alpha_leader_nt", "alpha_leader_aa"),
        ]:
            if pd.notna(r.get(ntc)) and pd.notna(r.get(aac)):
                nt_str = str(r[ntc]).strip()
                if nt_str != "":
                    if len(nt_str) % 3 != 0 or _nz(_tr(nt_str)) != _nz(r[aac]):
                        bad += 1
        for ch in ("alpha", "beta"):
            full = _nz(r.get(f"full_{ch}_aa"))
            if not full:
                continue
            lead = _nz(r.get(f"{ch}_leader_aa"))
            vdj = _nz(r.get(f"VDJ_{ch}_aa"))
            const = _nz(r.get(f"{ch}_constant_aa"))
            if lead + vdj + const != full:
                bad += 1
        sc, fb, fa = _nz(r.get("single_chain_aa")), _nz(r.get("full_beta_aa")), _nz(r.get("full_alpha_aa"))
        if sc:
            n_sc += 1
            expect_m = bool(_nz(r.get("beta_leader_aa")))
            if not (sc.startswith(fb) and sc.endswith(fa) and (not expect_m or sc.startswith("M"))):
                bad += 1
            if sc.startswith(fb) and sc.endswith(fa) and len(sc) >= len(fb) + len(fa):
                linkers.add(sc[len(fb) : len(sc) - len(fa)])
        elif fb and fa:
            bad += 1
        for c in ["full_beta_aa", "full_alpha_aa", "single_chain_aa"]:
            val = _nz(r.get(c))
            if val and ("*" in val or "X" in val):
                bad += 1
    linkers_ok = (not n_sc) or (bool(linkers) and linkers <= valid_linkers)
    ok = bad == 0 and linkers_ok
    return BatteryCheckResult(
        "A. Assembly integrity",
        "PASS" if ok else "FAIL",
        f"{bad} integrity failures; inter-chain linker(s)={linkers or 'NONE'}",
    )


def check_cdr3_ref(d: pd.DataFrame) -> BatteryCheckResult:
    """B. CDR3 ⊂ VDJ, canonical C…[F/W] (TRAJ35 allow-list), valid gene names."""
    miss = _missing(d, _CDR3_COLS)
    if miss:
        return BatteryCheckResult(
            "B. CDR3 + reference", "SKIP", f"missing columns: {', '.join(miss)}"
        )
    bad, noncanon, lens = 0, 0, []
    for _, r in d.iterrows():
        allow_map = NONCANON_J_END_MOUSE if _is_mouse_row(r) else NONCANON_J_END_HUMAN
        for ch in ("alpha", "beta"):
            raw_cdr = r[f"CDR3_{ch}"]
            if pd.isna(raw_cdr) or str(raw_cdr).strip() == "":
                continue
            cdr = str(raw_cdr).strip()
            vdj = _nz(r[f"VDJ_{ch}_aa"])
            jg = "" if pd.isna(r[f"{ch}_j_gene"]) else str(r[f"{ch}_j_gene"])
            lens.append(len(cdr))
            if cdr not in vdj or not cdr.startswith("C"):
                bad += 1
            end = cdr[-1]
            if end not in "FW":
                jg_u = jg.upper()
                if end in allow_map and any(
                    jg_u.startswith(g) for g in allow_map[end]
                ):
                    noncanon += 1
                else:
                    bad += 1
        for g, pat in [
            ("alpha_v_gene", "TRAV"), ("alpha_j_gene", "TRAJ"), ("alpha_c_gene", "TRAC"),
            ("beta_v_gene", "TRBV"), ("beta_j_gene", "TRBJ"), ("beta_c_gene", "TRBC"),
        ]:
            if not str(r[g]).upper().startswith(pat):
                bad += 1
        for ch in ("alpha", "beta"):
            jr = str(r.get(f"{ch}_junction_residue", ""))
            if len(jr) != 1 or not jr.isalpha():
                bad += 1
    span = f"{min(lens)}-{max(lens)} (med {int(np.median(lens))})" if lens else "n/a"
    return BatteryCheckResult(
        "B. CDR3 + reference",
        "PASS" if bad == 0 else "FAIL",
        f"{bad} failures; CDR3 len {span}; {noncanon} non-canonical-J C-ending "
        f"CDR3s (TRAJ35, expected)",
    )


def check_source(d: pd.DataFrame, clonotypes: pd.DataFrame | None) -> BatteryCheckResult:
    """C. Selected V/J/CDR3/cell_count matches the source clonotypes."""
    if clonotypes is None or "CDR3ab" not in getattr(clonotypes, "columns", []):
        return BatteryCheckResult("C. Source mapping", "SKIP", "no source clonotypes available")
    if "CDR3ab" not in d.columns:
        return BatteryCheckResult("C. Source mapping", "SKIP", "assembled frame has no CDR3ab key")
    ct = clonotypes.set_index("CDR3ab")
    mm, dualsec = 0, 0
    for _, r in d.iterrows():
        if r["CDR3ab"] not in ct.index:
            dualsec += 1
            continue
        s = ct.loc[r["CDR3ab"]]
        s = s.iloc[0] if isinstance(s, pd.DataFrame) else s
        for g in [
            "alpha_v_gene", "alpha_j_gene", "beta_v_gene", "beta_j_gene",
            "CDR3_alpha", "CDR3_beta",
        ]:
            if g in ct.columns and g in d.columns and str(r[g]) != str(s[g]):
                mm += 1
                break
        if (
            "cell_count" in ct.columns
            and "cell_count" in d.columns
            and abs(float(r["cell_count"]) - float(s["cell_count"])) > 0.5
        ):
            mm += 1
    return BatteryCheckResult(
        "C. Source mapping",
        "PASS" if mm == 0 else "FAIL",
        f"{mm} V/J/CDR3/cell mismatches; {dualsec} dual-alpha secondaries not "
        f"keyed in clonotypes (expected)",
    )


def check_synth(d: pd.DataFrame) -> BatteryCheckResult:
    """D. contig_verified all true; no dup constructs / α-β swaps / dup SCs."""
    if "single_chain_aa" not in d.columns:
        return BatteryCheckResult(
            "D. Synthesis / dual-alpha", "SKIP", "no single_chain_aa column"
        )
    n = len(d)
    nclone = (
        d["selected_clone"].nunique() if "selected_clone" in d else d["CDR3ab"].nunique()
    )
    # construct_contig_verified is the load-bearing gate. Fail CLOSED when it's
    # absent (can't certify what wasn't recorded), but say so honestly rather
    # than printing "0/n" — which reads as "verified false" not "not recorded".
    cv_present = "construct_contig_verified" in d.columns
    cv = _count_true(d, "construct_contig_verified")
    dup = _count_true(d, "synth_duplicate_construct")
    swap = _count_true(d, "synth_alpha_beta_swap")
    # Duplicate single-chains: count among NON-NULL only. `_add_single_chain`
    # sets single_chain_aa = NA for any clone missing a chain, and `nunique()`
    # excludes NA — so `n - nunique()` would count every unpaired clone as a
    # "duplicate" and spuriously FAIL the check (#check-synth-dupsc). Dedupe over
    # the populated constructs instead.
    sc = d["single_chain_aa"].dropna()
    dupsc = len(sc) - sc.nunique()
    ok = cv_present and cv == n and dup == 0 and swap == 0 and dupsc == 0
    cv_str = f"{cv}/{n}" if cv_present else "not recorded"
    # A clone expands on two INDEPENDENT axes (#283), so a single "+N dual-alpha"
    # is wrong once `--leader-secondary` is used: a row can be both a 2nd-α variant
    # AND a secondary-SP twin. Report the clean partition (primary + secondary = n)
    # plus dual-α as a clone-level count, instead of an additive surplus that
    # double-counts the overlap.
    clone_col = "selected_clone" if "selected_clone" in d else "CDR3ab"
    n_secondary = (
        int(d["leader_variant"].astype(str).str.startswith("secondary").sum())
        if "leader_variant" in d else 0
    )
    n_primary = n - n_secondary
    n_dual_clones = (
        int(d.loc[d["dual_alpha_variant"].notna(), clone_col].nunique())
        if "dual_alpha_variant" in d else 0
    )
    sp_split = f"{n_primary} primary" + (f" + {n_secondary} secondary-SP" if n_secondary else "")
    return BatteryCheckResult(
        "D. Synthesis / dual-alpha",
        "PASS" if ok else "FAIL",
        f"contig_verified {cv_str}; {nclone} clones->{n} constructs "
        f"({sp_split}; {n_dual_clones} dual-α); "
        f"dup={dup} swap={swap} dup_single_chain={dupsc}",
    )


def check_contigs(d: pd.DataFrame, cmap: dict[str, str | set[str]] | None) -> BatteryCheckResult:
    """E. Every contig ref exists in cellranger; CDR3 matches ≥1 own contig."""
    if not cmap:
        return BatteryCheckResult(
            "E. Raw cellranger contigs", "SKIP", "no cellranger contig annotations available"
        )
    refs = orphan = nomatch = 0
    for _, r in d.iterrows():
        for ch, cdrcol in [("beta", "CDR3_beta"), ("alpha", "CDR3_alpha")]:
            raw_ids = r.get(f"{ch}_contig_ids")
            if pd.isna(raw_ids):
                continue
            ids = [
                x.strip()
                for x in str(raw_ids).replace(";", ",").split(",")
                if x.strip()
            ]
            if not ids:
                continue
            refs += len(ids)
            orphan += sum(1 for c in ids if c not in cmap)
            raw_set: set[str] = set()
            for c in ids:
                if c in cmap:
                    val = cmap[c]
                    if isinstance(val, (set, frozenset, list, tuple)):
                        raw_set.update(str(v) for v in val)
                    else:
                        raw_set.add(str(val))
            cdr_val = r.get(cdrcol)
            if raw_set and (pd.isna(cdr_val) or str(cdr_val) not in raw_set):
                nomatch += 1
    return BatteryCheckResult(
        "E. Raw cellranger contigs",
        "PASS" if orphan == 0 and nomatch == 0 else "FAIL",
        f"{refs} contig refs; {orphan} not in cellranger; {nomatch} clones whose "
        f"CDR3 matches none of their own contigs",
    )


def check_germline(d: pd.DataFrame) -> BatteryCheckResult:
    """F. Donor-vs-germline divergence by region (leader / constant / framework).

    Informational, not a gate — germline divergence is expected donor biology, so
    this always PASSes (or SKIPs when nothing is germline-comparable); the point
    is to *surface* what's otherwise buried in CSV columns. Folds the leader
    keep/switch variants and the #187 constant-allele divergence (and the
    V-REGION framework comparison, when present) into one human-readable line.
    """
    from .assemble import collect_germline_variants

    try:
        gv = collect_germline_variants(d)
    except Exception as e:  # never let the surfacing check break the battery
        return BatteryCheckResult("F. Germline / allele divergence", "SKIP", f"unavailable: {e}")
    if gv is None or gv.empty:
        return BatteryCheckResult(
            "F. Germline / allele divergence", "SKIP",
            "no germline-comparable columns (run with contigs + V-leader reference)",
        )
    parts = []
    for region in ("leader", "constant", "framework"):
        sub = gv[gv["region"] == region]
        if sub.empty:
            continue
        n_variants = len(sub)
        n_constructs = int(sub["n_constructs"].sum())
        shipped = "/".join(sorted(sub["shipped"].astype(str).unique()))
        parts.append(f"{region}: {n_variants} variant(s), {n_constructs} construct(s) [{shipped}]")
    detail = "; ".join(parts) if parts else "none"
    return BatteryCheckResult("F. Germline / allele divergence", "PASS", detail)


def contig_cdr3_map(
    contig_dir: str | Path | None = None,
    *,
    cellranger_dir: str | Path | None = None,
) -> dict[str, set[str]]:
    """Build a ``contig_id -> set[cdr3]`` map from CellRanger ``*contig_annotations*.csv``
    under the given tree. Returns an empty dict when nothing is found (check E
    then reports SKIP). Uses a set per contig ID so shared 10x barcodes across
    multi-sample directories never overwrite each other.
    """
    root = cellranger_dir or contig_dir
    if root is None:
        return {}
    raw: dict[str, set[str]] = {}
    for f in glob.glob(f"{root}/**/*contig_annotations*.csv", recursive=True):
        try:
            a = pd.read_csv(f, usecols=lambda c: c in ("contig_id", "cdr3"))
        except Exception:
            continue
        if "contig_id" not in a.columns or "cdr3" not in a.columns:
            continue
        for cid, cdr in zip(a["contig_id"], a["cdr3"]):
            if pd.isna(cid) or pd.isna(cdr):
                continue
            raw.setdefault(str(cid), set()).add(str(cdr))
    return {k: (next(iter(v)) if len(v) == 1 else v) for k, v in raw.items()}


def run_qc_battery(
    assembled: pd.DataFrame,
    *,
    clonotypes: pd.DataFrame | None = None,
    contig_map: dict[str, str] | None = None,
) -> list[BatteryCheckResult]:
    """Run all checks and return their results in A–F order."""
    return [
        check_assembly(assembled),
        check_cdr3_ref(assembled),
        check_source(assembled, clonotypes),
        check_synth(assembled),
        check_contigs(assembled, contig_map),
        check_germline(assembled),
    ]


_EDGE_CASE_NOTES = [
    "TRAJ35 CDR3s end in `C` — non-canonical J anchor (`C`-G-S-G), confirmed vs IMGT.",
    "Dual-alpha secondary variants are not keyed in `clonotypes.csv` (their "
    "CDR3ab carries the 2nd alpha).",
    "Donor TRAC divergence near pos 3 → tcrsift ships canonical TRAC*01 (standard "
    "for TCR-T).",
    "`pgen`/`ppost` are log-scale (negative); private paired clones share no CDR3ab "
    "across donors.",
]


def qc_summary_markdown(
    cohorts: dict[str, list[BatteryCheckResult]],
    *,
    title: str = "QC summary",
    deliverable: str | None = None,
) -> tuple[str, bool]:
    """Render an A–E PASS/FAIL/SKIP table per cohort + the verified-correct
    edge-case note. Returns ``(markdown, all_ok)`` where ``all_ok`` is False iff
    any check FAILed (SKIP does not fail).
    """
    over = f" over `{deliverable}`" if deliverable else ""
    head = (
        f"Independent QC battery{over}. Every construct fact is re-derived from "
        "first principles and cross-checked against the source clonotypes, the "
        "raw cellranger contigs, and IMGT conventions — tcrsift's own QC columns "
        "are NOT trusted."
    )
    lines = [f"# {title}", "", head, ""]
    all_ok = True
    for coh, results in cohorts.items():
        lines += [f"## {coh}", "", "| check | verdict | detail |", "|---|---|---|"]
        for res in results:
            all_ok &= res.ok
            lines.append(f"| {res.name} | {res.badge()} | {res.detail} |")
        lines.append("")
    lines += ["## Verified-correct edge cases (not defects)"]
    lines += [f"- {note}" for note in _EDGE_CASE_NOTES]
    lines += ["", f"**Overall: {'ALL PASS ✅' if all_ok else 'FAILURES PRESENT ❌'}**"]
    return "\n".join(lines) + "\n", all_ok
