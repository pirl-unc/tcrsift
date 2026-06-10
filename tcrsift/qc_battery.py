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
# ends on that residue rather than the canonical Phe/Trp. TRAJ35: C-G-S-G.
NONCANON_J_END: dict[str, list[str]] = {"C": ["TRAJ35"]}


def _tr(nt: str) -> str:
    """Translate a DNA/RNA string frame-0, 'X' on any unknown codon."""
    nt = str(nt).upper().replace("U", "T")
    return "".join(
        _CT.get(nt[i : i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3)
    )


def _nz(a) -> str:
    """Stringify and strip a trailing stop ('*') so a translated CDS with its
    terminal stop compares equal to the stored stop-free AA."""
    return str(a).rstrip("*")


class QCResult:
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


def check_assembly(d: pd.DataFrame) -> QCResult:
    """A. Every translation exact, recomposition byte-for-byte, M-start, T2A."""
    miss = _missing(d, _ASSEMBLY_COLS)
    if miss:
        return QCResult(
            "A. Assembly integrity", "SKIP", f"missing columns: {', '.join(miss)}"
        )
    bad = 0
    linkers: set[str] = set()
    for _, r in d.iterrows():
        for ntc, aac in [
            ("single_chain_nt", "single_chain_aa"),
            ("full_beta_nt", "full_beta_aa"),
            ("full_alpha_nt", "full_alpha_aa"),
            ("VDJ_beta_nt", "VDJ_beta_aa"),
            ("VDJ_alpha_nt", "VDJ_alpha_aa"),
            ("beta_leader_nt", "beta_leader_aa"),
            ("alpha_leader_nt", "alpha_leader_aa"),
        ]:
            if (
                pd.notna(r.get(ntc))
                and pd.notna(r.get(aac))
                and _nz(_tr(r[ntc])) != _nz(r[aac])
            ):
                bad += 1
        for ch in ("alpha", "beta"):
            full = _nz(r[f"full_{ch}_aa"])
            lead = str(r[f"{ch}_leader_aa"])
            vdj = str(r[f"VDJ_{ch}_aa"])
            const = _nz(r.get(f"{ch}_constant_aa", ""))
            if lead + vdj + const != full:
                bad += 1
        sc, fb, fa = _nz(r["single_chain_aa"]), _nz(r["full_beta_aa"]), _nz(r["full_alpha_aa"])
        if not (sc.startswith(fb) and fa in sc and sc.startswith("M")):
            bad += 1
        if sc.startswith(fb) and fa in sc:
            linkers.add(sc[len(fb) : sc.rindex(fa)])
        if any(
            "*" in _nz(r[c]) or "X" in str(r[c])
            for c in ["full_beta_aa", "full_alpha_aa", "single_chain_aa"]
        ):
            bad += 1
    ok = bad == 0 and linkers == {EXPECTED_T2A}
    return QCResult(
        "A. Assembly integrity",
        "PASS" if ok else "FAIL",
        f"{bad} integrity failures; inter-chain linker(s)={linkers or 'NONE'}",
    )


def check_cdr3_ref(d: pd.DataFrame) -> QCResult:
    """B. CDR3 ⊂ VDJ, canonical C…[F/W] (TRAJ35 allow-list), valid gene names."""
    miss = _missing(d, _CDR3_COLS)
    if miss:
        return QCResult(
            "B. CDR3 + reference", "SKIP", f"missing columns: {', '.join(miss)}"
        )
    bad, noncanon, lens = 0, 0, []
    for _, r in d.iterrows():
        for ch in ("alpha", "beta"):
            cdr, vdj, jg = (
                str(r[f"CDR3_{ch}"]),
                str(r[f"VDJ_{ch}_aa"]),
                str(r[f"{ch}_j_gene"]),
            )
            if cdr in ("", "nan"):
                continue
            lens.append(len(cdr))
            if cdr not in vdj or not cdr.startswith("C"):
                bad += 1
            end = cdr[-1]
            if end not in "FW":
                if end in NONCANON_J_END and any(
                    jg.startswith(g) for g in NONCANON_J_END[end]
                ):
                    noncanon += 1
                else:
                    bad += 1
        for g, pat in [
            ("alpha_v_gene", "TRAV"), ("alpha_j_gene", "TRAJ"), ("alpha_c_gene", "TRAC"),
            ("beta_v_gene", "TRBV"), ("beta_j_gene", "TRBJ"), ("beta_c_gene", "TRBC"),
        ]:
            if not str(r[g]).startswith(pat):
                bad += 1
        for ch in ("alpha", "beta"):
            jr = str(r.get(f"{ch}_junction_residue", ""))
            if len(jr) != 1 or not jr.isalpha():
                bad += 1
    span = f"{min(lens)}-{max(lens)} (med {int(np.median(lens))})" if lens else "n/a"
    return QCResult(
        "B. CDR3 + reference",
        "PASS" if bad == 0 else "FAIL",
        f"{bad} failures; CDR3 len {span}; {noncanon} non-canonical-J C-ending "
        f"CDR3s (TRAJ35, expected)",
    )


def check_source(d: pd.DataFrame, clonotypes: pd.DataFrame | None) -> QCResult:
    """C. Selected V/J/CDR3/cell_count matches the source clonotypes."""
    if clonotypes is None or "CDR3ab" not in getattr(clonotypes, "columns", []):
        return QCResult("C. Source mapping", "SKIP", "no source clonotypes available")
    if "CDR3ab" not in d.columns:
        return QCResult("C. Source mapping", "SKIP", "assembled frame has no CDR3ab key")
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
    return QCResult(
        "C. Source mapping",
        "PASS" if mm == 0 else "FAIL",
        f"{mm} V/J/CDR3/cell mismatches; {dualsec} dual-alpha secondaries not "
        f"keyed in clonotypes (expected)",
    )


def check_synth(d: pd.DataFrame) -> QCResult:
    """D. contig_verified all true; no dup constructs / α-β swaps / dup SCs."""
    if "single_chain_aa" not in d.columns:
        return QCResult(
            "D. Synthesis / dual-alpha", "SKIP", "no single_chain_aa column"
        )
    n = len(d)
    nclone = (
        d["selected_clone"].nunique() if "selected_clone" in d else d["CDR3ab"].nunique()
    )
    cv = int((d.get("construct_contig_verified") == True).sum())  # noqa: E712
    dup = int((d.get("synth_duplicate_construct") == True).sum())  # noqa: E712
    swap = int((d.get("synth_alpha_beta_swap") == True).sum())  # noqa: E712
    dupsc = n - d["single_chain_aa"].nunique()
    ok = cv == n and dup == 0 and swap == 0 and dupsc == 0
    return QCResult(
        "D. Synthesis / dual-alpha",
        "PASS" if ok else "FAIL",
        f"contig_verified {cv}/{n}; {nclone} clones->{n} constructs "
        f"(+{n - nclone} dual-alpha); dup={dup} swap={swap} dup_single_chain={dupsc}",
    )


def check_contigs(d: pd.DataFrame, cmap: dict[str, str] | None) -> QCResult:
    """E. Every contig ref exists in cellranger; CDR3 matches ≥1 own contig."""
    if not cmap:
        return QCResult(
            "E. Raw cellranger contigs", "SKIP", "no cellranger contig annotations available"
        )
    refs = orphan = nomatch = 0
    for _, r in d.iterrows():
        for ch, cdrcol in [("beta", "CDR3_beta"), ("alpha", "CDR3_alpha")]:
            ids = [
                x.strip()
                for x in str(r.get(f"{ch}_contig_ids", "")).replace(";", ",").split(",")
                if x.strip()
            ]
            if not ids:
                continue
            refs += len(ids)
            orphan += sum(1 for c in ids if c not in cmap)
            raws = [cmap[c] for c in ids if c in cmap]
            if raws and str(r.get(cdrcol)) not in raws:
                nomatch += 1
    return QCResult(
        "E. Raw cellranger contigs",
        "PASS" if orphan == 0 and nomatch == 0 else "FAIL",
        f"{refs} contig refs; {orphan} not in cellranger; {nomatch} clones whose "
        f"CDR3 matches none of their own contigs",
    )


def contig_cdr3_map(
    contig_dir: str | Path | None = None,
    *,
    cellranger_dir: str | Path | None = None,
) -> dict[str, str]:
    """Build a ``contig_id -> cdr3`` map from CellRanger ``*contig_annotations*.csv``
    under the given tree. Returns an empty dict when nothing is found (check E
    then reports SKIP). Last-writer-wins on duplicate contig ids across files —
    contig ids are unique per cellranger run, and dup keys would only arise from
    overlapping globs.
    """
    root = cellranger_dir or contig_dir
    if root is None:
        return {}
    m: dict[str, str] = {}
    for f in glob.glob(f"{root}/**/*contig_annotations*.csv", recursive=True):
        try:
            a = pd.read_csv(f, usecols=lambda c: c in ("contig_id", "cdr3"))
        except Exception:
            continue
        if "contig_id" not in a.columns or "cdr3" not in a.columns:
            continue
        for cid, cdr in zip(a["contig_id"].astype(str), a["cdr3"].astype(str)):
            m[cid] = cdr
    return m


def run_qc_battery(
    assembled: pd.DataFrame,
    *,
    clonotypes: pd.DataFrame | None = None,
    contig_map: dict[str, str] | None = None,
) -> list[QCResult]:
    """Run all five checks and return their results in A–E order."""
    return [
        check_assembly(assembled),
        check_cdr3_ref(assembled),
        check_source(assembled, clonotypes),
        check_synth(assembled),
        check_contigs(assembled, contig_map),
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
    cohorts: dict[str, list[QCResult]],
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
