#!/usr/bin/env python3
"""Draw the infectious-disease and blood benchmark figures from ``benchmarks/public_infectious/results``.

    python scripts/figures/make_infectious_benchmark_figures.py

Every data panel plots empirical measurements written by
``benchmarks/public_infectious/run.py`` on public single-cell VDJ + GEX cohorts.
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

RESULTS = Path(__file__).resolve().parents[2] / "benchmarks" / "public_infectious" / "results"
SIG_LABELS = {
    "clone_size": "Clone size (baseline)",
    "AntigenExperienced": "AntigenExperienced",
    "CirculatingMemory": "CirculatingMemory",
    "AcuteActivation": "AcuteActivation",
    "Differentiated": "Differentiated",
    "Cytolytic": "Cytolytic",
    "Proliferation": "Proliferation",
    "NeoTCR_PBL": "NeoTCR_PBL",
    "TumorReactive": "TumorReactive",
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
    """Schematic of the three infectious-disease cohorts and benchmark labels."""
    columns = ["Cohort", "Sequenced", "How specificity was established", "Benchmark role"]
    rows = [
        (
            "10x Genomics\ndextramer blood",
            [
                "2 healthy donors\n(Donors 1 and 2,\nresting blood)",
                "Single-cell RNA + TCR\n+ CITE-seq + 50 pMHC\ndextramers on CD8 T cells",
                "44 viral and self pMHC\ndextramers calibrated\nagainst 6 negative controls",
                "Viral memory vs MART-1\nnaive, multimer QC,\nand VDJdb matching",
            ],
        ),
        (
            "Minervina 2022\nSARS-CoV-2",
            [
                "Longitudinal blood:\npre-exposure, infection,\nvaccination, breakthrough",
                "Single-cell RNA + TCR\nof dextramer-sorted\nCD8 T cells",
                "18 HLA-matched\nSARS-CoV-2 Spike and\nnon-Spike dextramers",
                "State trajectories across\n5 exposure stages and\nVDJdb matching",
            ],
        ),
        (
            "Fischer 2021\nCOVID-19 PBMC",
            [
                "2 hospitalized severe\nCOVID-19 donors\n(GT_2 and GT_6)",
                "Paired unstimulated and\n16 h Spike-peptide\nstimulated PBMC RNA + TCR",
                "Within-clone activation\nshift between paired\nunstim and stim wells",
                "Persistent clones:\nabundance fold-change\nvs cytokine/AIM shift",
            ],
        ),
    ]
    fig, ax = plt.subplots(figsize=(11, 4.9), constrained_layout=True)
    width, height, x0, dx, dy = 2.25, 1.05, 1.75, 2.6, 1.4
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
            1.5,
            y - height / 2,
            study,
            ha="right",
            va="center",
            fontsize=9.5,
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
                fontsize=8.6,
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
    save_tight(fig, "benchmark_infectious_datasets")


def paired_donor_forest(
    df: pd.DataFrame,
    comparison: str,
    title: str,
    out_name: str,
    xlim: tuple[float, float] = (0.0, 1.0),
) -> None:
    """Two-donor AUROC forest plot (10x Donor 1 and Donor 2)."""
    sub = df[df["comparison"].eq(comparison)].copy()
    d1 = sub[sub["dataset"].eq("tenx_donor1")].set_index("signature")
    d2 = sub[sub["dataset"].eq("tenx_donor2")].set_index("signature")
    sigs = [s for s in d1.index if s in d2.index]
    order = sorted(sigs, key=lambda s: 0.5 * (d1.loc[s, "auroc"] + d2.loc[s, "auroc"]))

    fig, ax = plt.subplots(figsize=(7.6, 4.8), constrained_layout=True)
    offset = 0.14
    y_trans = ax.get_yaxis_transform()
    for i, sig in enumerate(order):
        r1, r2 = d1.loc[sig], d2.loc[sig]
        c1 = NEUTRAL if sig == "clone_size" else BLUE
        c2 = INK2 if sig == "clone_size" else AQUA
        ax.plot([r1.auroc_low, r1.auroc_high], [i + offset, i + offset], color=c1, lw=2, solid_capstyle="round")
        ax.scatter([r1.auroc], [i + offset], s=48, color=c1, edgecolor=SURFACE, linewidth=1.3, zorder=3)
        ax.plot([r2.auroc_low, r2.auroc_high], [i - offset, i - offset], color=c2, lw=2, solid_capstyle="round")
        ax.scatter([r2.auroc], [i - offset], s=48, marker="s", color=c2, edgecolor=SURFACE, linewidth=1.3, zorder=3)
        ax.text(
            1.02,
            i,
            f"{r1.auroc:.2f} / {r2.auroc:.2f}",
            transform=y_trans,
            va="center",
            fontsize=8.6,
            color=INK,
        )

    ax.axvline(0.5, color=INK2, lw=1, ls=(0, (4, 3)))
    ax.text(0.49, len(order) - 0.45, "chance", ha="right", va="center", fontsize=8.5, color=INK2)
    ax.set_yticks(range(len(order)), [SIG_LABELS.get(s, s) for s in order])
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.6, len(order) - 0.3)
    ax.set_xlabel("AUROC (95% bootstrap interval)")
    ax.set_title(title)
    ax.grid(axis="x", lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)

    n1_pos, n1_neg = int(d1.iloc[0]["n_positive"]), int(d1.iloc[0]["n_negative"])
    n2_pos, n2_neg = int(d2.iloc[0]["n_positive"]), int(d2.iloc[0]["n_negative"])
    handles = [
        matplotlib.lines.Line2D([0], [0], color=BLUE, marker="o", lw=2, ms=6, label=f"Donor 1 ({n1_pos} vs {n1_neg})"),
        matplotlib.lines.Line2D([0], [0], color=AQUA, marker="s", lw=2, ms=6, label=f"Donor 2 ({n2_pos} vs {n2_neg})"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8.6, loc="lower right")
    save_tight(fig, out_name)


def blood_signature_aurocs() -> None:
    df = pd.read_csv(RESULTS / "blood_signature_auroc.csv")
    paired_donor_forest(
        df,
        "viral_memory_vs_naive_mart1",
        "Resting blood: viral memory (Flu/EBV/CMV) vs naive MART-1 clones",
        "benchmark_infectious_auroc_memory_vs_naive",
        xlim=(0.0, 1.0),
    )
    paired_donor_forest(
        df,
        "expanded_viral_vs_unlabeled_expanded",
        "Expanded blood clones (>=2 cells): viral dextramer+ vs unlabeled",
        "benchmark_infectious_auroc_expanded_blood",
        xlim=(0.2, 1.0),
    )
    paired_donor_forest(
        df,
        "influenza_m1_vs_ebv_expanded",
        "Memory subtype (>=2 cells): Influenza M1 (>0.5) vs EBV (<0.5)",
        "benchmark_infectious_auroc_flu_vs_ebv",
        xlim=(0.1, 1.0),
    )


def exposure_trajectories() -> None:
    exp = pd.read_csv(RESULTS / "exposure_state_profiles.csv").set_index("donor_category")
    auroc_df = pd.read_csv(RESULTS / "blood_signature_auroc.csv")
    recall = auroc_df[
        auroc_df["dataset"].eq("minervina_2022")
        & auroc_df["comparison"].eq("acute_post_vax_dose2_vs_resting_convalescent")
    ].sort_values("auroc")

    stages = [
        ("naive", "Pre-exposure"),
        ("post_inf", "Post-infection\n(convalescent)"),
        ("inf_after_1", "Infection +\n1 vaccine dose"),
        ("inf_after_2", "Acute recall\npost-dose 2"),
        ("break", "Breakthrough\ninfection"),
    ]
    present = [(k, lbl) for k, lbl in stages if k in exp.index]
    xs = np.arange(len(present))

    fig1, ax1 = plt.subplots(figsize=(7.6, 4.3), constrained_layout=True)
    series = [
        ("mean_Proliferation", "Proliferation", NEUTRAL, "o", "--"),
        ("mean_AcuteActivation", "AcuteActivation", ORANGE, "o", "-"),
        ("mean_Differentiated", "Differentiated", BLUE, "s", "-"),
        ("mean_TumorReactive", "TumorReactive", INK2, "P", "-."),
        ("mean_AIM", "AIM", INK, "D", "--"),
        ("mean_CirculatingMemory", "CirculatingMemory", AQUA, "^", "-"),
    ]
    for col, label, color, marker, ls in series:
        if col not in exp.columns:
            continue
        ys = [float(exp.loc[k, col]) for k, _ in present]
        ax1.plot(xs, ys, color=color, marker=marker, ls=ls, lw=2, ms=5.5, label=label)

    ax1.axhline(0.0, color=INK2, lw=0.9, ls=(0, (4, 3)))
    ax1.set_xticks(xs, [lbl for _, lbl in present])
    ax1.set_ylabel("mean signature score across SARS-CoV-2 CD8 T cells")
    ax1.set_title("SARS-CoV-2 CD8 T-cell state trajectories (Minervina 2022)")
    ax1.grid(axis="y", lw=0.6)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=8.4, ncol=2, loc="upper left")
    save_tight(fig1, "benchmark_infectious_exposure_states")

    fig2, ax2 = plt.subplots(figsize=(6.8, 4.2), constrained_layout=True)
    y_trans = ax2.get_yaxis_transform()
    for i, row in enumerate(recall.itertuples()):
        base = row.signature == "clone_size"
        color = NEUTRAL if base else BLUE
        ax2.plot([row.auroc_low, row.auroc_high], [i, i], color=color, lw=2, solid_capstyle="round")
        ax2.scatter([row.auroc], [i], s=52, color=color, edgecolor=SURFACE, linewidth=1.4, zorder=3)
        ax2.text(1.02, i, f"{row.auroc:.2f}", transform=y_trans, va="center", fontsize=8.8, color=INK)
    ax2.axvline(0.5, color=INK2, lw=1, ls=(0, (4, 3)))
    ax2.text(0.49, len(recall) - 0.45, "chance", ha="right", va="center", fontsize=8.5, color=INK2)
    r0 = recall.iloc[0]
    ax2.set_yticks(range(len(recall)), [SIG_LABELS.get(s, s) for s in recall["signature"]])
    ax2.set_xlim(0.2, 1.0)
    ax2.set_ylim(-0.6, len(recall) - 0.3)
    ax2.set_xlabel("AUROC (95% bootstrap interval)")
    ax2.set_title(
        f"Acute post-dose-2 recall vs convalescent clones ({r0.n_positive} vs {r0.n_negative})"
    )
    ax2.grid(axis="x", lw=0.6)
    ax2.set_axisbelow(True)
    ax2.tick_params(axis="y", length=0)
    ax2.spines["left"].set_visible(False)
    save_tight(fig2, "benchmark_infectious_auroc_recall")


def vdjdb_precision_recall() -> None:
    df = pd.read_csv(RESULTS / "vdjdb_matching.csv")
    df["rule"] = df["match_mode"] + ":" + df["match_strictness"]
    rule_order = [
        ("exact:strict_ab", "strict_ab"),
        ("exact:ab_with_partial", "ab_with_partial"),
        ("exact:b_only", "b_only"),
        ("levenshtein:b_only", "levenshtein b_only"),
    ]

    fig1, ax1 = plt.subplots(figsize=(6.8, 4.5), constrained_layout=True)

    # 10x Donor 1 + Donor 2 pooled exact-epitope precision vs recall
    tenx = df[df["dataset"].isin(["tenx_donor1", "tenx_donor2"])].copy()
    pooled = (
        tenx.groupby(["db_scope", "circular", "rule"], as_index=False)[
            ["n_query_clones", "n_matched", "n_exact_epitope"]
        ]
        .sum()
    )
    pooled["precision"] = pooled["n_exact_epitope"] / pooled["n_matched"].where(pooled["n_matched"] > 0, np.nan)
    pooled["recall"] = pooled["n_exact_epitope"] / pooled["n_query_clones"]

    panel_a_offsets = {
        ("full_vdjdb", "strict_ab"): (-8, -11, "right"),
        ("full_vdjdb", "levenshtein b_only"): (-7, 7, "right"),
        ("leave_study_out", "strict_ab"): (6, 6, "left"),
        ("leave_study_out", "ab_with_partial"): (6, 7, "left"),
        ("leave_study_out", "b_only"): (6, -12, "left"),
        ("leave_study_out", "levenshtein b_only"): (-7, -13, "right"),
    }
    for scope, label, color, ls, hollow in (
        ("full_vdjdb", "Full VDJdb (circular: 10x in DB)", ORANGE, (0, (3, 2)), True),
        ("leave_study_out", "10x study held out (independent)", BLUE, "-", False),
    ):
        sub = pooled[pooled["db_scope"].eq(scope)].set_index("rule")
        pts = [(sub.loc[r, "recall"], sub.loc[r, "precision"], short) for r, short in rule_order if r in sub.index]
        ax1.plot([p[0] for p in pts], [p[1] for p in pts], color=color, ls=ls, lw=1.8, zorder=2)
        for rec, prec, short in pts:
            ax1.scatter(
                [rec],
                [prec],
                s=64,
                color=SURFACE if hollow else color,
                edgecolor=color if hollow else SURFACE,
                linewidth=1.8 if hollow else 1.4,
                zorder=3,
            )
            off = panel_a_offsets.get((scope, short))
            if off is not None:
                ax1.annotate(
                    short,
                    (rec, prec),
                    textcoords="offset points",
                    xytext=(off[0], off[1]),
                    ha=off[2],
                    fontsize=8.4,
                    color=INK,
                )
        ax1.plot([], [], color=color, ls=ls, marker="o", mfc=SURFACE if hollow else color, mec=color, lw=1.8, label=label)

    percent_axis(ax1, "x", decimals=0)
    percent_axis(ax1, "y", decimals=0)
    ax1.set_xlim(0.0, 1.03)
    ax1.set_ylim(0.55, 1.05)
    ax1.set_xlabel("exact-epitope recall")
    ax1.set_ylabel("exact-epitope precision")
    ax1.set_title("10x Donors 1 & 2 (Flu, EBV, CMV, MART-1): VDJdb matching")
    ax1.grid(lw=0.6)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=8.5, loc="lower left")
    save_tight(fig1, "benchmark_infectious_vdjdb_10x")

    # Minervina 2022 SARS-CoV-2 virus precision vs recall
    fig2, ax2 = plt.subplots(figsize=(6.8, 4.5), constrained_layout=True)
    minv = df[df["dataset"].eq("minervina_2022_sarscov2")].copy()
    n_sars_query = int(minv["n_query_clones"].iloc[0])
    panel_b_offsets = {
        ("full_vdjdb", "strict_ab"): (-6, 7, "right"),
        ("full_vdjdb", "ab_with_partial"): (6, 7, "left"),
        ("full_vdjdb", "b_only"): (6, -11, "left"),
        ("full_vdjdb", "levenshtein b_only"): (6, 6, "left"),
        ("leave_study_out", "strict_ab"): (6, 6, "left"),
        ("leave_study_out", "ab_with_partial"): (6, 7, "left"),
        ("leave_study_out", "b_only"): (6, -12, "left"),
        ("leave_study_out", "levenshtein b_only"): (6, 6, "left"),
    }
    for scope, label, color, ls, hollow in (
        ("full_vdjdb", "Full VDJdb (circular: study in DB)", ORANGE, (0, (3, 2)), True),
        ("leave_study_out", "PMID:35383307 held out (independent)", BLUE, "-", False),
    ):
        sub = minv[minv["db_scope"].eq(scope)].set_index("rule")
        pts = [
            (sub.loc[r, "virus_recall"], sub.loc[r, "virus_precision"], short)
            for r, short in rule_order
            if r in sub.index
        ]
        ax2.plot([p[0] for p in pts], [p[1] for p in pts], color=color, ls=ls, lw=1.8, zorder=2)
        for rec, prec, short in pts:
            ax2.scatter(
                [rec],
                [prec],
                s=64,
                color=SURFACE if hollow else color,
                edgecolor=color if hollow else SURFACE,
                linewidth=1.8 if hollow else 1.4,
                zorder=3,
            )
            off = panel_b_offsets.get((scope, short), (6, 6, "left"))
            ax2.annotate(
                short,
                (rec, prec),
                textcoords="offset points",
                xytext=(off[0], off[1]),
                ha=off[2],
                fontsize=8.4,
                color=INK,
            )
        ax2.plot([], [], color=color, ls=ls, marker="o", mfc=SURFACE if hollow else color, mec=color, lw=1.8, label=label)

    percent_axis(ax2, "x", decimals=0)
    percent_axis(ax2, "y", decimals=0)
    ax2.set_xlim(-0.01, 0.36)
    ax2.set_ylim(0.12, 1.06)
    ax2.set_xlabel("SARS-CoV-2 recall")
    ax2.set_ylabel("SARS-CoV-2 species precision")
    ax2.set_title(f"Minervina 2022 ({n_sars_query:,} SARS-CoV-2 clones): VDJdb matching")
    ax2.grid(lw=0.6)
    ax2.set_axisbelow(True)
    ax2.legend(frameon=False, fontsize=8.5, loc="lower right")
    save_tight(fig2, "benchmark_infectious_vdjdb_minervina")


def multimer_concordance_qc() -> None:
    df = pd.read_csv(RESULTS / "multimer_calling.csv")
    y_col = "mean_memory_clr_diff" if "mean_memory_clr_diff" in df.columns else None
    if y_col is None:
        df["mean_memory_clr_diff"] = df["mean_cd45ro_log1p"] - df["mean_cd45ra_log1p"]
        y_col = "mean_memory_clr_diff"

    fig, ax = plt.subplots(figsize=(7.8, 4.6), constrained_layout=True)
    custom_labels = {
        ("tenx_donor1", "KLGGALQAK"): ("KLGGALQAK (D1, n={n})", (10, 10), "left"),
        ("tenx_donor2", "KLGGALQAK"): ("KLGGALQAK (D2, n={n})", (10, -13), "left"),
        ("tenx_donor1", "GILGFVFTL"): ("GILGFVFTL (D1, n={n})", (-10, 10), "right"),
        ("tenx_donor2", "GILGFVFTL"): ("GILGFVFTL (D2, n={n})", (-12, 14), "right"),
        ("tenx_donor2", "RAKFKQLL"): ("RAKFKQLL (D2, n={n})", (-12, -16), "right"),
        ("tenx_donor1", "IVTDFSVIK"): ("IVTDFSVIK (D1, n={n})", (-12, 14), "right"),
        ("tenx_donor2", "GLCTLVAML"): ("GLCTLVAML (D2, n={n})", (-10, -13), "right"),
        ("tenx_donor2", "FLYALALLL"): ("FLYALALLL (D2, n={n})", (-12, -14), "right"),
    }

    for row in df.itertuples():
        is_sticky = bool(row.is_known_sticky)
        is_hiv_bg = "HIV" in row.multimer
        color = ORANGE if is_sticky else NEUTRAL if is_hiv_bg else BLUE
        marker = "s" if row.dataset == "tenx_donor2" else "o"
        size = np.clip(28 + 2.0 * row.expanded_clones_with_any_call, 42, 210)
        y_val = getattr(row, y_col)
        ax.scatter(
            [row.median_clone_concordance],
            [y_val],
            s=size,
            marker=marker,
            color=color,
            edgecolor=SURFACE,
            linewidth=1.3,
            alpha=0.92,
            zorder=3,
        )
        pep = row.multimer.split("_")[0]
        spec = custom_labels.get((row.dataset, pep))
        if spec is not None:
            tmpl, xy_off, ha = spec
            ax.annotate(
                tmpl.format(n=row.expanded_clones_with_any_call),
                (row.median_clone_concordance, y_val),
                textcoords="offset points",
                xytext=xy_off,
                ha=ha,
                fontsize=8.0,
                color=INK,
                arrowprops={"arrowstyle": "-", "color": INK2, "lw": 0.5, "shrinkA": 0, "shrinkB": 3},
            )

    ax.axvline(0.5, color=INK2, lw=1.1, ls=(0, (4, 3)))
    ax.axhline(0.0, color=NEUTRAL, lw=0.9, ls=":")
    percent_axis(ax, "x", decimals=0)
    ax.set_xlim(0.10, 1.07)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin - 0.18, ymax + 0.32)
    ax.text(
        0.515,
        ax.get_ylim()[1] - 0.08,
        "50% concordance\nQC floor",
        ha="left",
        va="top",
        fontsize=8.2,
        color=INK2,
    )
    ax.set_xlabel("median clonal multimer concordance across expanded clones")
    ax.set_ylabel("CITE-seq memory shift: CLR(CD45RO) - CLR(CD45RA)")
    ax.set_title("Clonal concordance separates viral dextramers from sticky KLGGALQAK")
    ax.grid(lw=0.6)
    ax.set_axisbelow(True)

    handles = [
        matplotlib.lines.Line2D([0], [0], color=BLUE, marker="o", lw=0, ms=6.5, label="Viral multimer (Donor 1 o / Donor 2 s)"),
        matplotlib.lines.Line2D([0], [0], color=ORANGE, marker="o", lw=0, ms=6.5, label="Sticky KLGGALQAK_IE-1_CMV"),
        matplotlib.lines.Line2D([0], [0], color=NEUTRAL, marker="s", lw=0, ms=6.5, label="HIV Gag background (HIV-neg donor)"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8.2, loc="upper left")
    save_tight(fig, "benchmark_infectious_multimer_qc")


def fischer_stimulation() -> None:
    df = pd.read_csv(RESULTS / "fischer_stimulation.csv")
    fig, ax = plt.subplots(figsize=(7.4, 4.5), constrained_layout=True)

    reactive = df["delta_effector_cytokines"] > 0.5
    bystander = ~reactive
    ax.scatter(
        df.loc[bystander, "log2_fold_change"],
        df.loc[bystander, "delta_effector_cytokines"],
        s=25 + 1.4 * (df.loc[bystander, "stim_cells"] + df.loc[bystander, "control_cells"]).clip(upper=80),
        color=NEUTRAL,
        edgecolor=SURFACE,
        linewidth=0.9,
        alpha=0.85,
        label=f"Unchanged / bystander clones ({int(bystander.sum())})",
        zorder=2,
    )
    ax.scatter(
        df.loc[reactive, "log2_fold_change"],
        df.loc[reactive, "delta_effector_cytokines"],
        s=35 + 1.8 * (df.loc[reactive, "stim_cells"] + df.loc[reactive, "control_cells"]).clip(upper=80),
        color=BLUE,
        edgecolor=SURFACE,
        linewidth=1.2,
        alpha=0.95,
        label=f"Cytokine-induced clones, Δ > 0.5 ({int(reactive.sum())})",
        zorder=3,
    )

    top4 = df[reactive].nlargest(4, "delta_effector_cytokines").sort_values(
        "delta_effector_cytokines", ascending=False
    )
    offsets = [(7, 6), (7, -11), (7, 5), (7, -10)]
    for row, xy_off in zip(top4.itertuples(), offsets):
        ax.annotate(
            f"{row.CDR3_beta} ({row.donor}, {row.stim_cells} vs {row.control_cells})",
            (row.log2_fold_change, row.delta_effector_cytokines),
            textcoords="offset points",
            xytext=xy_off,
            fontsize=8.0,
            color=INK,
        )

    ax.axhline(0.5, color=BLUE, lw=1.0, ls=(0, (4, 3)))
    ax.axvline(0.0, color=INK2, lw=0.9, ls=":")
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax + 0.35)
    ax.set_xlabel("clone abundance log2 fold-change (16 h Spike stim vs unstim)")
    ax.set_ylabel("within-clone Δ effector-cytokine transcripts (IFNG, TNF, IL2, CCL4)")
    ax.set_title("16 h Spike restimulation induces cytokine mRNAs without clonal division")
    ax.grid(lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8.4, loc="upper right")
    save_tight(fig, "benchmark_infectious_fischer_stim")


if __name__ == "__main__":
    dataset_designs()
    blood_signature_aurocs()
    exposure_trajectories()
    vdjdb_precision_recall()
    multimer_concordance_qc()
    fischer_stimulation()
    print("wrote infectious benchmark figures")
