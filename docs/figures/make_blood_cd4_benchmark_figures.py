#!/usr/bin/env python3
"""Draw the paired blood-tumor, post-ACT, and CD4+ TIL benchmark figures.

    python docs/figures/make_blood_cd4_benchmark_figures.py
"""

from __future__ import annotations

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

RESULTS = Path(__file__).resolve().parents[2] / "benchmarks" / "public_blood_cd4" / "results"
SIG_LABELS = {
    "clone_size": "Clone size (baseline)",
    "NeoTCR_PBL": "NeoTCR_PBL",
    "NeoTCR4": "NeoTCR4",
    "NeoTCR8": "NeoTCR8",
    "MANAscore": "MANAscore",
    "TumorReactive": "TumorReactive",
    "CirculatingMemory": "CirculatingMemory",
    "Differentiated": "Differentiated",
    "AntigenExperienced": "AntigenExperienced",
    "Cytolytic": "Cytolytic",
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
    columns = ["Cohort", "Sequenced", "Experimental ground truth", "Benchmark role"]
    rows = [
        (
            "Pulliam 2024\nMCC blood + tumor\n(GSE227708)",
            [
                "4 Merkel cell carcinoma\npatients (10 blood & 2\ntumor libraries >=20 cells)",
                "Hashed 10x GEX + VDJ\nacross pre- and on-\nanti-PD-(L)1 draws",
                "pMHC multimers for\nMCPyV tumor virus vs\nCMV, EBV, and Influenza",
                "--context blood-tumor\n(NeoTCR_PBL) & tumor-\ninfiltrating clone recall",
            ],
        ),
        (
            "Chiffelle 2024\nPre/Post-ACT TIL\n(GSE222448)",
            [
                "6 melanoma TIL-ACT\npatients (P1, P5, P7,\nP8, P10, P13)",
                "Paired baseline T0 (8)\nand Day-30 post-ACT\nT30 (6) tumor biopsies",
                "Cloned TCRs tested vs\nautologous melanoma\nlines (TR vs NTR)",
                "Post-infusion T30 vs T0\nsignature stability &\nACTP persistence",
            ],
        ),
        (
            "Veatch 2022\nMelanoma CD4\n(GSE198265)",
            [
                "4 melanoma patients\n(X197, X198, X205,\nX422) fresh TILs",
                "BD Rhapsody 397-gene\ntargeted mRNA + AbSeq\n+ dominant TCR Vβ",
                "Cloned CD4 TCRs screened\nagainst autologous\nneoantigens & self",
                "Verified CD4_neo vs\nuntested/self & FOXP3+\nTreg audit",
            ],
        ),
        (
            "Oh 2020\nBladder CD4\n(GSE149652)",
            [
                "7 bladder cancer\npatients (7 tumors +\n6 matched normals)",
                "Whole-transcriptome\n10x GEX + paired α/β\nTCR on sorted CD4+",
                "Curated CD4CXCL13,\ncytotoxic GZMB/GZMK,\nFOXP3+ Treg & normal",
                "Whole-transcriptome\nNeoTCR4 (40/40), pre_sorted\n& --exclude-treg",
            ],
        ),
    ]
    fig, ax = plt.subplots(figsize=(11.2, 6.2), constrained_layout=True)
    width, height, x0, dx, dy = 2.25, 1.05, 1.85, 2.55, 1.38
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
            fontsize=9.2,
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
                fontsize=8.4,
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
    save_tight(fig, "blood_cd4_benchmark_datasets")


def _paired_forest(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    order: list[str],
    title: str,
    label1: str,
    label2: str,
    out_name: str,
    color2: str = AQUA,
    xlim: tuple[float, float] = (0.04, 1.02),
) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
    off = 0.15
    y_t = ax.get_yaxis_transform()
    for i, sig in enumerate(order):
        r1, r2 = df1.loc[sig], df2.loc[sig]
        c1 = NEUTRAL if sig == "clone_size" else BLUE
        c2 = INK2 if sig == "clone_size" else color2
        ax.plot([r1.auroc_low, r1.auroc_high], [i + off, i + off], color=c1, lw=2, solid_capstyle="round")
        ax.scatter([r1.auroc], [i + off], s=48, color=c1, edgecolor=SURFACE, linewidth=1.2, zorder=3)
        ax.plot([r2.auroc_low, r2.auroc_high], [i - off, i - off], color=c2, lw=2, solid_capstyle="round")
        ax.scatter([r2.auroc], [i - off], s=48, marker="s", color=c2, edgecolor=SURFACE, linewidth=1.2, zorder=3)
        ax.text(1.02, i, f"{r1.auroc:.2f} / {r2.auroc:.2f}", transform=y_t, va="center", fontsize=8.6, color=INK)

    ax.axvline(0.5, color=INK2, lw=1, ls=(0, (4, 3)))
    ax.set_yticks(range(len(order)), [SIG_LABELS.get(s, s) for s in order])
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.6, len(order) - 0.3)
    ax.set_xlabel("AUROC (95% bootstrap interval)")
    ax.set_title(title, fontsize=10.2)
    ax.grid(axis="x", lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    handles = [
        matplotlib.lines.Line2D([0], [0], color=BLUE, marker="o", lw=2, ms=5.8, label=label1),
        matplotlib.lines.Line2D([0], [0], color=color2, marker="s", lw=2, ms=5.8, label=label2),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8.4, loc="lower right")
    save_tight(fig, out_name)


def pulliam_and_act_aurocs() -> None:
    df = pd.read_csv(RESULTS / "auroc_blood_cd4.csv")

    p_vir = df[
        df["dataset"].eq("pulliam_blood")
        & df["subset"].eq("all_clones_ge1")
        & df["contrast"].eq("MCPyV_vs_CMV_EBV_Flu")
    ].set_index("signature")
    p_tum = df[
        df["dataset"].eq("pulliam_blood")
        & df["subset"].eq("matched_tumor_donors_ge1")
        & df["contrast"].eq("tumor_infiltrating_ab_vs_blood_only")
    ].set_index("signature")
    sigs_p = [s for s in p_vir.index if s in p_tum.index]
    order_p = sorted(sigs_p, key=lambda s: float(p_vir.loc[s, "auroc"]))
    _paired_forest(
        p_vir,
        p_tum,
        order_p,
        "Pulliam MCC peripheral blood (--context blood-tumor)",
        f"MCPyV vs CMV/EBV/Flu ({int(p_vir.iloc[0]['n_positive'])} vs {int(p_vir.iloc[0]['n_negative'])})",
        f"In matched tumor vs blood-only ({int(p_tum.iloc[0]['n_positive'])} vs {int(p_tum.iloc[0]['n_negative'])})",
        "blood_cd4_auroc_pulliam",
        color2=AQUA,
        xlim=(0.04, 1.02),
    )

    c_t0 = df[
        df["dataset"].eq("chiffelle_act_6patients")
        & df["subset"].eq("T0_pre_ACT")
        & df["contrast"].eq("TR_vs_NTR")
    ].set_index("signature")
    c_t30 = df[
        df["dataset"].eq("chiffelle_act_6patients")
        & df["subset"].eq("T30_post_ACT")
        & df["contrast"].eq("TR_vs_NTR")
    ].set_index("signature")
    sigs_c = [s for s in c_t0.index if s in c_t30.index]
    order_c = sorted(sigs_c, key=lambda s: 0.5 * (float(c_t0.loc[s, "auroc"]) + float(c_t30.loc[s, "auroc"])))
    _paired_forest(
        c_t0,
        c_t30,
        order_c,
        "Chiffelle melanoma: Pre-ACT T0 vs Day-30 Post-ACT T30",
        f"Pre-ACT T0 ({int(c_t0.iloc[0]['n_positive'])} TR vs {int(c_t0.iloc[0]['n_negative'])} NTR)",
        f"Post-ACT T30 ({int(c_t30.iloc[0]['n_positive'])} TR vs {int(c_t30.iloc[0]['n_negative'])} NTR)",
        "blood_cd4_auroc_chiffelle_act",
        color2=AQUA,
        xlim=(0.04, 1.02),
    )


def cd4_til_aurocs() -> None:
    df = pd.read_csv(RESULTS / "auroc_blood_cd4.csv")

    v_all = df[
        df["dataset"].eq("veatch_melanoma_cd4") & df["contrast"].eq("CD4_neo_vs_untested_cd4")
    ].set_index("signature")
    v_nt = df[
        df["dataset"].eq("veatch_melanoma_cd4") & df["contrast"].eq("CD4_neo_vs_non_treg_untested")
    ].set_index("signature")
    sigs_v = [s for s in v_all.index if s in v_nt.index]
    order_v = sorted(sigs_v, key=lambda s: float(v_nt.loc[s, "auroc"]))
    _paired_forest(
        v_all,
        v_nt,
        order_v,
        "Veatch melanoma CD4+ TIL: verified CD4_neo clones",
        f"vs all untested CD4 ({int(v_all.iloc[0]['n_positive'])} vs {int(v_all.iloc[0]['n_negative'])})",
        f"vs non-Treg CD4 ({int(v_nt.iloc[0]['n_positive'])} vs {int(v_nt.iloc[0]['n_negative'])})",
        "blood_cd4_auroc_veatch",
        color2=AQUA,
        xlim=(0.2, 1.02),
    )

    o_cm = df[
        df["dataset"].eq("oh_bladder_cd4") & df["contrast"].eq("CD4CXCL13_GZMB_GZMK_vs_CD4CM_TH17")
    ].set_index("signature")
    o_tr = df[
        df["dataset"].eq("oh_bladder_cd4") & df["contrast"].eq("CD4CXCL13_GZMB_GZMK_vs_CD4IL2RAHI_Treg")
    ].set_index("signature")
    sigs_o = [s for s in o_cm.index if s in o_tr.index]
    order_o = sorted(sigs_o, key=lambda s: float(o_cm.loc[s, "auroc"]))
    _paired_forest(
        o_cm,
        o_tr,
        order_o,
        "Oh bladder CD4+ TIL: CD4CXCL13 + cytotoxic effectors",
        f"vs CD4CM/TH17 bystanders ({int(o_cm.iloc[0]['n_positive'])} vs {int(o_cm.iloc[0]['n_negative'])})",
        f"vs FOXP3+ CD4IL2RAHI Treg ({int(o_tr.iloc[0]['n_positive'])} vs {int(o_tr.iloc[0]['n_negative'])})",
        "blood_cd4_auroc_oh",
        color2=ORANGE,
        xlim=(0.05, 1.02),
    )


def treg_and_presorted_ablation() -> None:
    oh = pd.read_csv(RESULTS / "shortlist_oh_cd4.csv").set_index("run")
    ve = pd.read_csv(RESULTS / "shortlist_veatch_cd4.csv").set_index("run")

    # Figure: Oh 2020 shortlist counts across pre_sorted & --exclude-treg settings
    fig1, ax1 = plt.subplots(figsize=(7.0, 4.4), constrained_layout=True)
    runs_oh = [
        ("no_presorted_cd4_mrna_only", "No pre_sorted\n(CD4 mRNA > 0)"),
        ("presorted_default", "pre_sorted=CD4\n(default)"),
        ("presorted_exclude_treg_0.50", "pre_sorted=CD4\n--exclude-treg 0.50"),
        ("presorted_exclude_treg_0.25", "pre_sorted=CD4\n--exclude-treg 0.25"),
    ]
    xs = np.arange(len(runs_oh))
    w = 0.34
    eff_c = [int(oh.loc[k, "eff_cxcl13_cyto_shortlisted"]) for k, _ in runs_oh]
    treg_c = [int(oh.loc[k, "treg_il2rahi_shortlisted"]) for k, _ in runs_oh]
    ax1.bar(xs - w / 2, eff_c, width=w, color=BLUE, label="CD4CXCL13 + cytotoxic GZMB/GZMK")
    ax1.bar(xs + w / 2, treg_c, width=w, color=ORANGE, label="CD4IL2RAHI (FOXP3+ Treg)")
    for x, ev, tv in zip(xs, eff_c, treg_c):
        ax1.text(x - w / 2, ev + 0.6, str(ev), ha="center", va="bottom", fontsize=8.5, color=INK)
        ax1.text(x + w / 2, tv + 0.6, str(tv), ha="center", va="bottom", fontsize=8.5, color=ORANGE)
    ax1.set_xticks(xs, [lbl for _, lbl in runs_oh], fontsize=8.5)
    ax1.set_ylim(0, max(max(eff_c), max(treg_c), 10) * 1.25)
    ax1.set_ylabel("shortlisted CD4 clones across 7 bladder tumors")
    ax1.set_title("Oh 2020: pre_sorted=CD4 rescue and --exclude-treg", fontsize=10.2)
    ax1.grid(axis="y", lw=0.6)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=8.4, loc="upper left")
    save_tight(fig1, "blood_cd4_oh_presorted_treg")

    # Figure: Veatch 2022 CD4_neo vs FOXP3+ Treg shortlist counts
    fig2, ax2 = plt.subplots(figsize=(6.2, 4.3), constrained_layout=True)
    runs_ve = [
        ("default_keep_treg", "Default\n(--no-exclude-treg)"),
        ("exclude_treg_0.50", "With --exclude-treg\n(FOXP3 >= 50%)"),
    ]
    xv = np.arange(len(runs_ve))
    neo_v = [int(ve.loc[k, "cd4_neo_shortlisted"]) for k, _ in runs_ve]
    treg_v = [int(ve.loc[k, "treg_shortlisted"]) for k, _ in runs_ve]
    ax2.bar(
        xv - w / 2,
        neo_v,
        width=w,
        color=BLUE,
        label=f"Verified CD4_neo (of {int(ve.iloc[0]['cd4_neo_ge2'])} >=2-cell)",
    )
    ax2.bar(xv + w / 2, treg_v, width=w, color=ORANGE, label="FOXP3+ Treg-like clones")
    for x, nv, tv in zip(xv, neo_v, treg_v):
        ax2.text(x - w / 2, nv + 1.0, str(nv), ha="center", va="bottom", fontsize=8.6, color=INK)
        ax2.text(x + w / 2, tv + 1.0, str(tv), ha="center", va="bottom", fontsize=8.6, color=ORANGE)
    ax2.set_xticks(xv, [lbl for _, lbl in runs_ve])
    ax2.set_ylim(0, max(max(neo_v), max(treg_v), 15) * 1.25)
    ax2.set_ylabel("shortlisted qualifying clones")
    ax2.set_title("Veatch 2022: --exclude-treg removes FOXP3+ Treg picks", fontsize=10.2)
    ax2.grid(axis="y", lw=0.6)
    ax2.set_axisbelow(True)
    ax2.legend(frameon=False, fontsize=8.4, loc="upper right")
    save_tight(fig2, "blood_cd4_veatch_treg")


if __name__ == "__main__":
    dataset_designs()
    pulliam_and_act_aurocs()
    cd4_til_aurocs()
    treg_and_presorted_ablation()
    print("wrote blood & CD4 benchmark figures")
