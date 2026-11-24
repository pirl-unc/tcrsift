#!/usr/bin/env python3
"""Draw the head-to-head comparison SVG figures from ``benchmarks/head_to_head/results``.

    python docs/figures/make_overview_benchmark_figures.py

Produces six single-panel SVGs in ``docs/assets/figures/``:
1. ``h2h_prioritize_auroc.svg`` — clone-level AUROC across 5 human ex vivo cohorts.
2. ``h2h_prioritize_cd4_treg.svg`` — CD4 shortlist Treg fraction on Veatch 2022 and Oh 2020.
3. ``h2h_assembly_stitchr.svg`` — construct assembly completion vs ``stitchr`` / ``thimble 1.3.1``.
4. ``h2h_pgen_olga.svg`` — cross-donor public sharing AUROC and speed vs ``OLGA 1.3.0``.
5. ``h2h_multimer_counts.svg`` — sticky ``KLGGALQAK`` / HIV Gag false positives vs canonical viral clones.
6. ``h2h_multimer_precision.svg`` — cognate epitope precision on 10x-held-out ``strict_ab`` VDJdb.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from make_figures import AQUA, BLUE, INK, INK2, NEUTRAL, ORANGE, OUT, SURFACE

H2H_RESULTS = Path(__file__).resolve().parents[2] / "benchmarks" / "head_to_head" / "results"


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


def h2h_prioritize_figures() -> None:
    """Single-panel AUROC and CD4 Treg shortlist figures across 5 ex vivo cohorts."""
    df = pd.read_csv(H2H_RESULTS / "h2h_prioritize_baselines.csv")

    cohort_order = [
        ("meng_pdac_cd8", "Meng 2023\nPDAC CD8 TIL"),
        ("chiffelle_melanoma_cd8", "Chiffelle 2024\nMelanoma CD8 TIL"),
        ("pulliam_mcc_blood_cd8", "Pulliam 2024\nMCC Blood CD8"),
        ("veatch_melanoma_cd4", "Veatch 2022\nMelanoma CD4 TIL"),
        ("oh_bladder_cd4", "Oh 2020\nBladder CD4 TIL"),
    ]
    method_specs = [
        ("clone_size", "Clone size", NEUTRAL),
        ("raw_CXCL13", "Single-gene CXCL13", ORANGE),
        ("scanpy_score_genes_NeoTCR", "Scanpy score_genes (NeoTCR)", AQUA),
        ("scanpy_composite_topk", "Scanpy composite top-K", INK2),
        ("tcrsift", "TCRsift default", BLUE),
    ]

    # Figure 1: Grouped horizontal bars of AUROC across the 5 cohorts
    fig1, ax1 = plt.subplots(figsize=(7.6, 5.2), constrained_layout=True)
    y_centers = np.arange(len(cohort_order))[::-1]
    bar_h = 0.145
    offsets = np.linspace(2 * bar_h, -2 * bar_h, len(method_specs))

    for m_idx, (m_key, m_label, color) in enumerate(method_specs):
        vals = []
        ys = []
        for c_idx, (c_key, _) in enumerate(cohort_order):
            sub = df[df["cohort"].eq(c_key)]
            if m_key == "tcrsift":
                row = sub[sub["method"].str.startswith("tcrsift")].iloc[0]
            else:
                row = sub[sub["method"].eq(m_key)].iloc[0]
            vals.append(float(row["auroc"]))
            ys.append(y_centers[c_idx] + offsets[m_idx])
        ax1.barh(
            ys,
            vals,
            height=bar_h * 0.88,
            color=color,
            edgecolor=SURFACE,
            linewidth=0.6,
            label=m_label,
        )
        if m_key == "tcrsift":
            for yv, v in zip(ys, vals):
                ax1.text(
                    v + 0.012,
                    yv,
                    f"{v:.2f}",
                    va="center",
                    ha="left",
                    fontsize=8.3,
                    fontweight="bold",
                    color=BLUE,
                )

    ax1.axvline(0.5, color=INK2, lw=1.0, ls=(0, (4, 3)))
    ax1.set_yticks(y_centers, [lbl for _, lbl in cohort_order])
    ax1.set_xlim(0.0, 0.98)
    ax1.set_xlabel("Clone-level AUROC")
    ax1.set_title("Ex vivo ranking AUROC across 5 human cohorts", fontsize=10.2)
    ax1.grid(axis="x", lw=0.6)
    ax1.set_axisbelow(True)
    ax1.spines["left"].set_visible(False)
    ax1.tick_params(axis="y", length=0)
    ax1.legend(frameon=False, fontsize=8.4, loc="lower right")
    save_tight(fig1, "h2h_prioritize_auroc")

    # Figure 2: CD4 Treg shortlist contamination on Veatch 2022 & Oh 2020
    fig2, ax2 = plt.subplots(figsize=(6.8, 4.4), constrained_layout=True)
    cd4_cohorts = [
        ("veatch_melanoma_cd4", "Veatch 2022\nMelanoma CD4 (FOXP3+ Treg)"),
        ("oh_bladder_cd4", "Oh 2020\nBladder CD4 (CD4IL2RAHI Treg)"),
    ]
    xs = np.arange(len(cd4_cohorts))
    bw = 0.15
    b_offsets = np.linspace(-2 * bw, 2 * bw, len(method_specs))
    ymax = 0.0
    for m_idx, (m_key, m_label, color) in enumerate(method_specs):
        treg_fracs = []
        for c_key, _ in cd4_cohorts:
            sub = df[df["cohort"].eq(c_key)]
            if m_key == "tcrsift":
                row = sub[sub["method"].str.startswith("tcrsift")].iloc[0]
            else:
                row = sub[sub["method"].eq(m_key)].iloc[0]
            tf = float(row["treg_shortlist_fraction"])
            treg_fracs.append(tf)
            ymax = max(ymax, tf)
        ax2.bar(
            xs + b_offsets[m_idx],
            treg_fracs,
            width=bw * 0.9,
            color=color,
            edgecolor=SURFACE,
            linewidth=0.6,
            label=m_label,
        )
        for xv, tf in zip(xs + b_offsets[m_idx], treg_fracs):
            ax2.text(
                xv,
                tf + 0.012,
                f"{tf:.0%}",
                ha="center",
                va="bottom",
                fontsize=8.2,
                fontweight="bold" if m_key == "tcrsift" else "normal",
                color=BLUE if m_key == "tcrsift" else INK,
            )

    percent_axis(ax2, "y", decimals=0)
    ax2.set_xticks(xs, [lbl for _, lbl in cd4_cohorts])
    ax2.set_ylim(0, max(ymax * 1.28, 0.35))
    ax2.set_ylabel("Treg fraction in top-K shortlist (lower is better)")
    ax2.set_title("CD4 shortlist Treg fraction (--exclude-treg)", fontsize=10.2)
    ax2.grid(axis="y", lw=0.6)
    ax2.set_axisbelow(True)
    ax2.legend(frameon=False, fontsize=8.2, loc="upper right")
    save_tight(fig2, "h2h_prioritize_cd4_treg")


def h2h_assembly_and_pgen_figures() -> None:
    """Single-panel construct assembly vs stitchr 1.3.1 and Pgen/Ppost vs OLGA 1.3.0."""
    asm = pd.read_csv(H2H_RESULTS / "h2h_assembly_vs_stitchr.csv")
    pgen = pd.read_csv(H2H_RESULTS / "h2h_pgen_vs_olga.csv")

    # Figure 3: Construct assembly completion rate & variant/error annotations
    fig1, ax1 = plt.subplots(figsize=(7.8, 4.9), constrained_layout=True)
    asm_bars = [
        ("minervina_2022_human_cr6_vdj", "tcrsift_donor_contig", "Minervina CR6\nTCRsift donor contig", BLUE),
        ("minervina_2022_human_cr6_vdj", "stitchr_thimble_v1.3.1", "Minervina CR6\nstitchr 1.3.1 (*01)", NEUTRAL),
        ("fischer_2021_5prime_truncated_human", "tcrsift_donor_contig_only", "Fischer 5'-trunc\nTCRsift strict contig", INK2),
        ("fischer_2021_5prime_truncated_human", "tcrsift_germline_vdj_fallback", "Fischer 5'-trunc\nTCRsift --germline-vdj-fallback", BLUE),
        ("fischer_2021_5prime_truncated_human", "stitchr_thimble_v1.3.1", "Fischer 5'-trunc\nstitchr 1.3.1", NEUTRAL),
        ("mouse_cellranger_101_mc38_khatun", "tcrsift_mouse_germline_fallback", "Mouse CR\nTCRsift --species mouse", BLUE),
        ("mouse_cellranger_101_mc38_khatun", "stitchr_thimble_v1.3.1_default", "Mouse CR\nstitchr default", ORANGE),
        ("mouse_cellranger_101_mc38_khatun", "stitchr_thimble_v1.3.1_skip_c_checks", "Mouse CR\nstitchr -sc (skip C check)", NEUTRAL),
    ]
    ys = np.arange(len(asm_bars))[::-1]
    rates = []
    notes = []
    colors = []
    labels = []
    for cohort, tool, lbl, col in asm_bars:
        r = asm[asm["cohort"].eq(cohort) & asm["tool"].eq(tool)].iloc[0]
        rate = float(r["complete_rate"])
        rates.append(rate)
        colors.append(col)
        n_tot = int(r["n_clones"])
        head, tail = lbl.split("\n", 1)
        labels.append(f"{head} (n={n_tot})\n{tail}")
        poly = int(r["donor_vdj_polymorphic_clones_preserved"])
        misc = int(r["mismatched_trbc_shipped"])
        if tool == "tcrsift_donor_contig" and poly > 0:
            notes.append(f"{rate:.1%} ({poly} donor V-variants kept)")
        elif misc > 0:
            notes.append(f"{rate:.1%} ({misc} Trbj/Trbc mismatches)")
        else:
            notes.append(f"{rate:.1%} ({int(r['n_complete'])}/{int(r['n_clones'])})")

    ax1.barh(ys, rates, height=0.64, color=colors, edgecolor=SURFACE, linewidth=0.8)
    for yv, rv, txt, col in zip(ys, rates, notes, colors):
        ax1.text(
            min(rv + 0.015, 1.02),
            yv,
            txt,
            va="center",
            ha="left",
            fontsize=8.2,
            fontweight="bold" if col == BLUE else "normal",
            color=ORANGE if "mismatches" in txt else INK,
        )
    percent_axis(ax1, "x", decimals=0)
    ax1.set_xticks([0.0, 0.25, 0.50, 0.75, 1.0])
    ax1.set_yticks(ys, labels, fontsize=8.4)
    ax1.set_xlim(0, 1.48)
    ax1.set_xlabel("Complete single-chain constructs built")
    ax1.set_title("Construct assembly: TCRsift vs stitchr (thimble 1.3.1)", fontsize=10.2)
    ax1.grid(axis="x", lw=0.6)
    ax1.set_axisbelow(True)
    ax1.spines["left"].set_visible(False)
    ax1.tick_params(axis="y", length=0)
    save_tight(fig1, "h2h_assembly_stitchr")

    # Figure 4: Cross-donor public sharing AUROC & throughput vs OLGA 1.3.0
    fig2, ax2 = plt.subplots(figsize=(6.8, 4.5), constrained_layout=True)
    pgen_cohorts = [
        ("human_blood_donors_1500", "Human blood\n(n=1,500 clones)"),
        ("mouse_khatun_gp66_1781", "Mouse GP66 CD4\n(n=1,781 clones)"),
    ]
    models = [
        ("tcrsift_kmer_pgen_beta", "TCRsift k-mer Pgen", AQUA),
        ("tcrsift_kmer_ppost_beta", "TCRsift k-mer Ppost", BLUE),
        ("olga_exact_vdj_pgen_beta", "OLGA 1.3.0 exact DP", NEUTRAL),
    ]
    xs = np.arange(len(pgen_cohorts))
    bw = 0.24
    offs = [-bw, 0.0, bw]
    for m_idx, (m_key, m_lbl, col) in enumerate(models):
        aucs = []
        spds = []
        r_vals = []
        for c_key, _ in pgen_cohorts:
            r = pgen[pgen["cohort"].eq(c_key) & pgen["model"].eq(m_key)].iloc[0]
            aucs.append(float(r["public_sharing_auroc"]))
            spds.append(float(r["clones_per_sec"]))
            r_vals.append(r["pearson_r_vs_olga_beta"])
        ax2.bar(
            xs + offs[m_idx],
            aucs,
            width=bw * 0.9,
            color=col,
            edgecolor=SURFACE,
            linewidth=0.8,
            label=m_lbl,
        )
        for xv, av, sp, rv in zip(xs + offs[m_idx], aucs, spds, r_vals):
            spd_txt = f"{sp / 1000:.0f}k/s" if sp >= 1000 else f"{sp:.0f}/s"
            extra = f"\n(r={float(rv):.2f})" if m_key == "tcrsift_kmer_pgen_beta" and pd.notna(rv) else ""
            ax2.text(
                xv,
                av + 0.012,
                f"{av:.2f}\n{spd_txt}{extra}",
                ha="center",
                va="bottom",
                fontsize=8.0,
                color=INK,
                linespacing=1.15,
            )

    ax2.axhline(0.5, color=INK2, lw=1.0, ls=(0, (4, 3)))
    ax2.set_xticks(xs, [lbl for _, lbl in pgen_cohorts])
    ax2.set_ylim(0.4, 1.06)
    ax2.set_ylabel("Cross-donor public sharing AUROC")
    ax2.set_title("Publicness AUROC and throughput: TCRsift k-mer vs OLGA 1.3.0", fontsize=10.2)
    ax2.grid(axis="y", lw=0.6)
    ax2.set_axisbelow(True)
    ax2.legend(frameon=False, fontsize=8.4, loc="lower right")
    save_tight(fig2, "h2h_pgen_olga")


def h2h_multimer_figures() -> None:
    """Single-panel pMHC multimer clone count and VDJdb precision figures on 10x Donors 1 & 2."""
    df = pd.read_csv(H2H_RESULTS / "h2h_multimer_callers.csv").set_index("caller")
    caller_order = [
        ("10x_per_cell_threshold", "10x per-cell\nthreshold", NEUTRAL),
        ("clr_z3_per_cell", "CLR z >= 3\nper-cell", INK2),
        ("gmm_bimodal_per_cell", "2-comp GMM\nper-cell", AQUA),
        ("tcrsift_clonal_concordance", "TCRsift clonal\nconcordance", BLUE),
    ]

    xs = np.arange(len(caller_order))
    sticky = [int(df.loc[k, "sticky_KLGGALQAK_fp_clones"]) for k, _, _ in caller_order]
    hiv = [int(df.loc[k, "hiv_gag_fp_clones"]) for k, _, _ in caller_order]
    canon = [int(df.loc[k, "canonical_flu_ebv_cmv_clones"]) for k, _, _ in caller_order]
    precs = [float(df.loc[k, "vdjdb_cognate_precision"]) for k, _, _ in caller_order]

    # Figure 5: Sticky KLGGALQAK / HIV Gag FPs vs canonical viral clones
    fig1, ax1 = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    w = 0.36
    ax1.bar(
        xs - w / 2,
        sticky,
        width=w,
        color=ORANGE,
        edgecolor=SURFACE,
        label="Sticky CMV KLGGALQAK + HIV Gag FPs",
    )
    ax1.bar(xs - w / 2, hiv, bottom=sticky, width=w, color="#b84a14", edgecolor=SURFACE)
    ax1.bar(
        xs + w / 2,
        canon,
        width=w,
        color=BLUE,
        edgecolor=SURFACE,
        label="Canonical Flu / EBV / CMV clones",
    )
    for x, s, h, c in zip(xs, sticky, hiv, canon):
        tot_fp = s + h
        ax1.text(
            x - w / 2,
            tot_fp + 18,
            f"{tot_fp:,}",
            ha="center",
            va="bottom",
            fontsize=8.4,
            fontweight="bold",
            color=ORANGE,
        )
        ax1.text(x + w / 2, c + 18, f"{c:,}", ha="center", va="bottom", fontsize=8.4, color=BLUE)

    ax1.set_xticks(xs, [lbl for _, lbl, _ in caller_order])
    ax1.set_ylim(0, max(max(s + h for s, h in zip(sticky, hiv)), max(canon)) * 1.18)
    ax1.set_ylabel("Called expanded clones (>=2 cells, Donors 1 + 2)")
    ax1.set_title("Sticky KLGGALQAK / HIV Gag false positives vs canonical viral clones", fontsize=10.0)
    ax1.grid(axis="y", lw=0.6)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=8.4, loc="upper right")
    save_tight(fig1, "h2h_multimer_counts")

    # Figure 6: Independent cognate precision on 10x-held-out strict_ab VDJdb
    fig2, ax2 = plt.subplots(figsize=(6.8, 4.3), constrained_layout=True)
    bars = ax2.bar(
        xs,
        precs,
        width=0.52,
        color=[col for _, _, col in caller_order],
        edgecolor=SURFACE,
        linewidth=0.9,
    )
    for rect, (k, _, col), pv in zip(bars, caller_order, precs):
        n_ok = int(df.loc[k, "vdjdb_cognate_epitope_concordant"])
        n_m = int(df.loc[k, "vdjdb_strict_matched_clones"])
        exact_p = (n_ok / n_m) if n_m > 0 else pv
        ax2.text(
            rect.get_x() + rect.get_width() / 2,
            pv + 0.015,
            f"{exact_p:.1%}\n({n_ok}/{n_m})",
            ha="center",
            va="bottom",
            fontsize=8.4,
            fontweight="bold" if col == BLUE else "normal",
            color=BLUE if col == BLUE else INK,
        )
    percent_axis(ax2, "y", decimals=0)
    ax2.set_xticks(xs, [lbl for _, lbl, _ in caller_order])
    ax2.set_ylim(0, 1.16)
    ax2.set_ylabel("Cognate epitope precision vs held-out VDJdb")
    ax2.set_title("Cognate precision on 10x-held-out strict_ab VDJdb", fontsize=10.2)
    ax2.grid(axis="y", lw=0.6)
    ax2.set_axisbelow(True)
    save_tight(fig2, "h2h_multimer_precision")


def main() -> None:
    h2h_prioritize_figures()
    h2h_assembly_and_pgen_figures()
    h2h_multimer_figures()
    print(f"Wrote 6 head-to-head SVGs to {OUT}")


if __name__ == "__main__":
    main()
