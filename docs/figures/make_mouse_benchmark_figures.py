#!/usr/bin/env python3
"""Draw the murine (`--species mouse`) benchmark figures from ``benchmarks/public_mouse/results``.

    python docs/figures/make_mouse_benchmark_figures.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from make_figures import AQUA, BLUE, INK, INK2, NEUTRAL, ORANGE, OUT, SURFACE

RESULTS = Path(__file__).resolve().parents[2] / "benchmarks" / "public_mouse" / "results"
SIG_LABELS = {
    "clone_size": "Clone size (baseline)",
    "TumorReactive": "TumorReactive (MGI)",
    "TumorReactive_HumanCXCL13": "TumorReactive (human CXCL13)",
    "Differentiated": "Differentiated",
    "AntigenExperienced": "AntigenExperienced",
    "Cytolytic": "Cytolytic",
    "CirculatingMemory": "CirculatingMemory",
    "AcuteActivation": "AcuteActivation",
    "Proliferation": "Proliferation",
    "AIM": "AIM",
}


def save_tight(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        OUT / f"{name}.svg",
        bbox_inches="tight",
        pad_inches=0.08,
        metadata={"Date": None},
    )
    plt.close(fig)


def percent_axis(ax, axis: str = "x", decimals: int = 0) -> None:
    getattr(ax, f"{axis}axis").set_major_formatter(
        matplotlib.ticker.PercentFormatter(1.0, decimals=decimals)
    )


def dataset_designs() -> None:
    columns = ["Cohort", "Sequenced", "Biological ground truth", "Benchmark role"]
    rows = [
        (
            "MC38 tumor\n(GSE168944)",
            [
                "2 C57BL/6 mice across\n4 matched organs:\nTumor, Blood, dLN, Spleen",
                "8 paired 10x GEX +\nVDJ libraries\n(CD8 + CD4 T cells)",
                "4-compartment clonal\ntracking (tumor-enriched\nvs lymphoid-shared)",
                "MGI TumorReactive,\nblood-tumor tracking,\nand mouse assembly",
            ],
        ),
        (
            "Kasmani 2022\nNOD diabetes\n(GSE200608)",
            [
                "NOD mice across 2 rounds:\npancreatic islets and\nmatched spleen",
                "6 paired 10x GEX +\nVDJ libraries on\nsorted CD8 T cells",
                "H-2Kd IGRP (KYNKANVFL)\ntetramer+ in islet vs\nspleen, plus islet Tet-",
                "Target-organ vs splenic\nsister clones & VDJdb\ncircular vs held-out",
            ],
        ),
        (
            "Khatun 2021\nLCMV GP66 CD4\n(GSE158896)",
            [
                "5 C57BL/6 mice (M1-M5)\nat Day 10 post-acute\nLCMV Armstrong",
                "5 paired 10x GEX +\nVDJ libraries on\nsplenic CD4 T cells",
                "I-Ab GP66-77\n(DIYKGVYQFKSV) MHC-II\ntetramer+ CD4 T cells",
                "Mouse CD4 phenotyping,\npublic GP66 sharing,\nand CD4 assembly",
            ],
        ),
    ]
    fig, ax = plt.subplots(figsize=(11, 4.9), constrained_layout=True)
    width, height, x0, dx, dy = 2.25, 1.05, 1.85, 2.55, 1.4
    for c, name in enumerate(columns):
        ax.text(
            x0 + c * dx + width / 2,
            0.42,
            name,
            ha="center",
            va="center",
            fontsize=9.5,
            fontweight="bold",
            color=INK,
        )
    for r, (study, cells) in enumerate(rows):
        y = -r * dy - 0.25
        ax.text(
            1.6,
            y - height / 2,
            study,
            ha="right",
            va="center",
            fontsize=9.3,
            fontweight="bold",
            color=INK,
        )
        for c, words in enumerate(cells):
            x = x0 + c * dx
            last = c == len(cells) - 1
            ax.add_patch(
                matplotlib.patches.FancyBboxPatch(
                    (x, y - height),
                    width,
                    height,
                    boxstyle="round,pad=0,rounding_size=0.08",
                    facecolor=BLUE if last else SURFACE,
                    edgecolor=BLUE if last else NEUTRAL,
                    linewidth=1.3,
                )
            )
            ax.text(
                x + width / 2,
                y - height / 2,
                words,
                ha="center",
                va="center",
                fontsize=8.5,
                color=SURFACE if last else INK,
                linespacing=1.25,
            )
            if not last:
                ax.annotate(
                    "",
                    xy=(x + dx - 0.04, y - height / 2),
                    xytext=(x + width + 0.04, y - height / 2),
                    arrowprops={"arrowstyle": "->", "color": INK2, "lw": 1.1},
                )
    ax.set_xlim(0, x0 + 3 * dx + width + 0.1)
    ax.set_ylim(-(len(rows) - 1) * dy - height - 0.4, 0.7)
    ax.axis("off")
    save_tight(fig, "mouse_benchmark_datasets")


def mouse_auroc_forest() -> None:
    df = pd.read_csv(RESULTS / "auroc_mouse.csv")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.7), constrained_layout=True)

    panels = [
        (
            ax1,
            df[df["dataset"].eq("mc38") & df["contrast"].eq("tumor_til_vs_dln_spleen_nontumor")],
            "MC38: tumor TILs vs matched dLN & spleen",
        ),
        (
            ax2,
            df[df["dataset"].eq("kasmani") & df["contrast"].eq("islet_tet_pos_vs_spleen_tet_pos")],
            "Kasmani IGRP Tet+: pancreatic islet vs spleen",
        ),
    ]
    for ax, sub, title in panels:
        sub = sub.sort_values("auroc").reset_index(drop=True)
        y_trans = ax.get_yaxis_transform()
        for i, row in sub.iterrows():
            is_base = row["signature"] == "clone_size"
            is_hcx = row["signature"] == "TumorReactive_HumanCXCL13"
            color = NEUTRAL if is_base else ORANGE if is_hcx else BLUE
            ax.plot([row["auroc_low"], row["auroc_high"]], [i, i], color=color, lw=2, solid_capstyle="round")
            ax.scatter(
                [row["auroc"]],
                [i],
                s=52,
                color=SURFACE if is_hcx else color,
                edgecolor=color if is_hcx else SURFACE,
                linewidth=1.6 if is_hcx else 1.3,
                zorder=3,
            )
            ax.text(1.02, i, f"{row['auroc']:.2f}", transform=y_trans, va="center", fontsize=8.5, color=INK)
        ax.axvline(0.5, color=INK2, lw=1, ls=(0, (4, 3)))
        ax.text(0.49, len(sub) - 0.4, "chance", ha="right", va="center", fontsize=8.2, color=INK2)
        r0 = sub.iloc[0]
        ax.set_yticks(range(len(sub)), [SIG_LABELS.get(s, s) for s in sub["signature"]])
        ax.set_xlim(0.04, 1.02)
        ax.set_ylim(-0.6, len(sub) - 0.3)
        ax.set_xlabel("AUROC (95% bootstrap interval)")
        ax.set_title(f"{title}\n({int(r0['n_positive'])} vs {int(r0['n_negative'])} clones)", fontsize=9.8)
        ax.grid(axis="x", lw=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)

    save_tight(fig, "mouse_auroc_mc38_kasmani")


def mouse_vdjdb_recovery() -> None:
    df = pd.read_csv(RESULTS / "vdjdb_mouse_recovery.csv")
    rules = ["strict_ab", "ab_with_partial", "b_only", "b_only_lev1"]
    rule_labels = ["strict_ab", "ab_with_partial", "b_only", "levenshtein\nb_only"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.4), constrained_layout=True)

    # Panel A: Kasmani IGRP circular vs held-out recall
    kas_circ = df[df["cohort"].eq("kasmani_igrp_nod") & df["database"].eq("vdjdb_full_circular")].set_index("strictness")
    kas_held = df[df["cohort"].eq("kasmani_igrp_nod") & df["database"].eq("vdjdb_heldout_no_kasmani")].set_index("strictness")
    xs = np.arange(len(rules))
    w = 0.34
    c_vals = [float(kas_circ.loc[r, "igrp_recall"]) for r in rules]
    h_vals = [float(kas_held.loc[r, "igrp_recall"]) for r in rules]

    ax1.bar(
        xs - w / 2,
        c_vals,
        width=w,
        facecolor=SURFACE,
        edgecolor=ORANGE,
        linewidth=1.8,
        linestyle="--",
        label="Full VDJdb (circular: PMID:35667687 in DB)",
    )
    ax1.bar(
        xs + w / 2,
        h_vals,
        width=w,
        facecolor=BLUE,
        edgecolor=BLUE,
        linewidth=1.2,
        label="Study held out (independent)",
    )
    for x, cv, hv in zip(xs, c_vals, h_vals):
        ax1.text(x - w / 2, cv + 0.015, f"{cv:.0%}", ha="center", va="bottom", fontsize=8.0, color=ORANGE)
        ax1.text(x + w / 2, hv + 0.015, f"{hv:.0%}", ha="center", va="bottom", fontsize=8.0, color=INK)

    percent_axis(ax1, "y", decimals=0)
    ax1.set_xticks(xs, rule_labels)
    ax1.set_ylim(0, max(max(c_vals), 0.25) * 1.22)
    ax1.set_ylabel("IGRP (KYNKANVFL) recall")
    n_q = int(kas_held.iloc[0]["n_query_clones"])
    ax1.set_title(f"Kasmani NOD IGRP CD8 clones (n={n_q})")
    ax1.grid(axis="y", lw=0.6)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=8.1, loc="upper left")

    # Panel B: Khatun GP66 CD4 LCMV recall vs precision (independent)
    kha = df[df["cohort"].eq("khatun_gp66_cd4")].set_index("strictness")
    r_vals = [float(kha.loc[r, "igrp_recall"]) for r in rules]
    p_vals = [
        float(kha.loc[r, "igrp_precision_among_matched"])
        if pd.notna(kha.loc[r, "igrp_precision_among_matched"])
        else 0.0
        for r in rules
    ]
    ax2.bar(xs - w / 2, r_vals, width=w, color=BLUE, label="LCMV recall")
    ax2.bar(xs + w / 2, p_vals, width=w, color=AQUA, label="LCMV precision among matched")
    for x, rv, pv, r in zip(xs, r_vals, p_vals, rules):
        ax2.text(x - w / 2, rv + 0.018, f"{rv:.1%}", ha="center", va="bottom", fontsize=7.9, color=INK)
        if int(kha.loc[r, "n_matched"]) > 0:
            ax2.text(x + w / 2, pv + 0.018, f"{pv:.0%}", ha="center", va="bottom", fontsize=7.9, color=INK)
        else:
            ax2.text(x + w / 2, 0.02, "0 hits", ha="center", va="bottom", fontsize=7.8, color=INK2)

    percent_axis(ax2, "y", decimals=0)
    ax2.set_xticks(xs, rule_labels)
    ax2.set_ylim(0, 1.15)
    ax2.set_ylabel("rate across unique GP66 CD4 clones")
    n_k = int(kha.iloc[0]["n_query_clones"])
    ax2.set_title(f"Khatun LCMV GP66 CD4 clones (n={n_k}, independent)")
    ax2.grid(axis="y", lw=0.6)
    ax2.set_axisbelow(True)
    ax2.legend(frameon=False, fontsize=8.1, loc="upper left")

    save_tight(fig, "mouse_vdjdb_recovery")


def mouse_assembly_and_shortlist() -> None:
    sh = pd.read_csv(RESULTS / "shortlist_mc38.csv").set_index("dataset")
    asm = pd.read_csv(RESULTS / "mouse_assembly.csv")
    summary = json.loads((RESULTS / "summary.json").read_text())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2), constrained_layout=True)

    # Panel A: Shortlist capture rates in MC38 Tumor & Blood
    cats = ["MC38 Tumor\n(--context solid-tumor)", "MC38 Blood\n(--context blood)"]
    keys = ["mc38_tumor", "mc38_blood"]
    pos_rates = [
        sh.loc[k, "tumor_enriched_shortlisted"] / max(sh.loc[k, "tumor_enriched_total"], 1)
        for k in keys
    ]
    neg_rates = [
        sh.loc[k, "lymphoid_shared_shortlisted"] / max(sh.loc[k, "lymphoid_shared_total"], 1)
        for k in keys
    ]
    xs = np.arange(2)
    w = 0.32
    ax1.bar(xs - w / 2, pos_rates, width=w, color=BLUE, label="Tumor-enriched / tumor-infiltrating")
    ax1.bar(xs + w / 2, neg_rates, width=w, color=NEUTRAL, label="Lymphoid-shared / blood-only")
    for x, k, pr, nr in zip(xs, keys, pos_rates, neg_rates):
        r = sh.loc[k]
        ax1.text(
            x - w / 2,
            pr + 0.015,
            f"{int(r['tumor_enriched_shortlisted'])}/{int(r['tumor_enriched_total'])}",
            ha="center",
            va="bottom",
            fontsize=8.2,
            color=INK,
        )
        ax1.text(
            x + w / 2,
            nr + 0.015,
            f"{int(r['lymphoid_shared_shortlisted'])}/{int(r['lymphoid_shared_total'])}",
            ha="center",
            va="bottom",
            fontsize=8.2,
            color=INK2,
        )
    percent_axis(ax1, "y", decimals=0)
    ax1.set_xticks(xs, cats)
    ax1.set_ylim(0, max(max(pos_rates), max(neg_rates), 0.25) * 1.25)
    ax1.set_ylabel("fraction shortlisted by default prioritize")
    ax1.set_title("MC38 multi-organ shortlist selectivity")
    ax1.grid(axis="y", lw=0.6)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=8.2, loc="upper right")

    # Panel B: Murine Trac + Trbc1/Trbc2 construct assembly breakdown
    labels_b = [
        f"MC38 CD8\n(n={int(asm.loc[0, 'n_candidates'])})",
        f"Khatun GP66 CD4\n(n={int(asm.loc[1, 'n_candidates'])})",
    ]
    trbc1_f = (asm["trbc1_constructs"] / asm["n_candidates"]).to_numpy()
    trbc2_f = (asm["trbc2_constructs"] / asm["n_candidates"]).to_numpy()
    xb = np.arange(len(asm))
    ax2.bar(xb, trbc1_f, width=0.45, color=BLUE, label="Trac + Trbc1 single-chain construct")
    ax2.bar(xb, trbc2_f, width=0.45, bottom=trbc1_f, color=AQUA, label="Trac + Trbc2 single-chain construct")
    for x, row in enumerate(asm.itertuples()):
        tot = row.complete_rate
        ax2.text(
            x,
            tot + 0.02,
            f"{row.n_complete_single_chain}/{row.n_candidates} (100%)\n{row.median_single_chain_aa_len} aa median",
            ha="center",
            va="bottom",
            fontsize=8.1,
            color=INK,
        )
    percent_axis(ax2, "y", decimals=0)
    ax2.set_xticks(xb, labels_b)
    ax2.set_ylim(0, 1.22)
    ax2.set_ylabel("assembled synthesis-ready constructs")
    n_err = summary["mc38"]["assembly_validation_errors"] + summary["khatun"]["assembled_cd4_validation_errors"]
    ax2.set_title(f"Murine germline V+J + C assembly ({n_err} validation errors)")
    ax2.grid(axis="y", lw=0.6)
    ax2.set_axisbelow(True)
    ax2.legend(frameon=False, fontsize=8.2, loc="lower right")

    save_tight(fig, "mouse_assembly_and_compartments")


if __name__ == "__main__":
    dataset_designs()
    mouse_auroc_forest()
    mouse_vdjdb_recovery()
    mouse_assembly_and_shortlist()
    print("wrote mouse benchmark figures")
