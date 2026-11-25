# Ground-Truth & Comparator Provenance (`benchmarks/head_to_head`)

This directory contains reproducible head-to-head comparisons of `tcrsift` against eight open-source bioinformatics tools and baselines across all six algorithmic stages of single-cell TCR construct generation, antigen-specificity annotation, and transcriptional prioritization:

1. **Ex vivo TIL & Cancer-Blood Clonotype Prioritization (`h2h_prioritize_baselines.csv`)**:
   - **Comparators**:
     - `clone_size`: Within-donor clonal abundance rank (`size_percentile`, breaking ties by `cell_count`).
     - `raw_CXCL13`: Per-clone mean `log1p(CP10K)` *CXCL13* single-gene expression.
     - `scanpy_score_genes_NeoTCR`: Standard `scanpy.tl.score_genes(..., ctrl_size=50, n_bins=25)` (Satija et al. 2015 binned-control module score, subtracting `genes_down` when signed) on published `NeoTCR8` (CD8 TIL), `NeoTCR4` (CD4 TIL), or `NeoTCR_PBL` (blood-tumor) gene sets, taking top-$K$ expanded (`>=2`-cell) clones at the exact shortlist size $K$ produced by `tcrsift`.
     - `pyucell_NeoTCR`: **`pyucell` `v0.8.0`** (rank-based Mann-Whitney U UCell scoring, Andreatta & Carmona 2021) on the exact same signed `NeoTCR8` / `NeoTCR4` / `NeoTCR_PBL` gene sets (`max_rank=min(1500, n_vars)`, `w_neg=1.0`), taking top-$K$ expanded clones.
     - `scanpy_composite_topk`: Equal-weighted mean of `scanpy.tl.score_genes` across all context preset signatures (`>=2`-cell clones, top-$K$ at matched budget $K$).
     - `tcrsift_solid_tumor` / `tcrsift_blood_tumor` / `tcrsift_solid_tumor_cd4_excl_treg`: Explicit context presets (`--context solid-tumor`, `--context blood-tumor`, or `--context solid-tumor --tcell-type cd4 --exclude-treg`), running within-sample gene-wise z-scoring (`score_by_name`, including `NeoTCR8`/`NeoTCR4`/`NeoTCR_PBL`/`MANAscore` built into `tcrsift`) + empirical cell-resampling null calibration (`calibrate_signature_background`) + priority-tiered round-robin (`select_round_robin`).
   - **Cohorts**: Meng et al. 2023 PDAC CD8+ TIL (`GSE254250`, `PMID:37992152`), Chiffelle et al. 2024 Melanoma baseline `T0` CD8+ TIL (`GSE222448`, `PMID:39321802`), Pulliam et al. 2024 MCC peripheral blood (`GSE227708`, `PMID:38340723`), Veatch et al. 2022 Melanoma CD4+ TIL (`GSE198265`, `PMID:35413271`), and Oh et al. 2020 Bladder CD4+ TIL (`GSE149652`, `PMID:32497499`).
   - **Why `predicTCR` and `TRTpred` are documented rather than executed**: `predicTCR` (Tan et al., *Nat Biotechnol* 2024) model weights are proprietary behind `predictcr.com`; `TRTpred` (Pétremand et al., *Nat Biotechnol* 2024; `Harari-Lab/TRTpred`) is an R/Seurat+Bioconductor package under a custom license trained on melanoma cohorts overlapping public melanoma test sets. All baselines in `run.py` are open-source Python packages, untrained on the evaluation cohorts, and reproducible offline.

2. **Full-Length TCR Construct Assembly (`h2h_assembly_vs_stitchr.csv`)**:
   - **Comparator**: **`stitchr` / `thimble` `v1.3.1`** (Heather et al., *Nucleic Acids Res* 2022, `DOI:10.1093/nar/gkac190`; `pip install stitchr`, initialized with `stitchrdl -s HUMAN` and `stitchrdl -s MOUSE`) using `Linker=T2A`, `Link_order=BA`.
   - **Cohorts**:
     - `minervina_2022_human_cr6_vdj`: All `46` shortlisted paired α/β human clones from `minervina_blood_default` (`Zenodo 6232103`, Cell Ranger 6 `fwr1_nt`–`fwr4_nt` donor VDJ sequences).
     - `fischer_2021_truncated_human`: All `160` paired α/β clones from Fischer et al. 2021 (`GSE171037`, `PMID:34301934`) where `98.7%–99.7%` of contigs are 5'-truncated (testing default auto-rescue `tcrsift assemble` vs fail-closed `--no-germline-vdj-fallback` vs `thimble -s HUMAN`).
     - `mouse_cellranger_101`: All `101` shortlisted murine clones (`41` MC38 CD8+ from `GSE168944` + `60` Khatun GP66 CD4+ from `GSE158896`, `PMID:33201171`) fed with raw Cell Ranger V/J/C annotations.

3. **Synthesis Codon Optimization & Hazard Clearance (`h2h_codon_opt_vs_dnachisel.csv`)**:
   - **Comparator**: **`DNAChisel` `v3.2.16`** (Zulkower & Rosser, *Bioinformatics* 2020; `DnaOptimizationProblem` with `EnforceTranslation`, `AvoidPattern` for `EcoRI`/`BamHI`/`XhoI`/`BsaI`/`BsmBI` and 9-nt homopolymers, `EnforceGCContent(0.30, 0.70, 50)`, and `CodonOptimize(species="h_sapiens")`) vs `tcrsift assemble` `codon_optimization="focal"` (default) and `"full"` on the `46` Minervina 2022 human single-chain constructs.

4. **Held-Out Public TCR Database Matching (`h2h_db_matching_vs_scirpy.csv`)**:
   - **Comparator**: **`scirpy` `v0.26.1`** (Sturm et al., *Bioinformatics* 2020; `ir.pp.ir_dist` + `ir.tl.ir_query` with `metric="identity"` paired/`VDJ` and `metric="levenshtein", cutoff=1` paired) vs `tcrsift annotate` (`strict_ab` exact, `ab_with_partial` exact, and `strict_ab` paired Levenshtein-1 `ab_near`) on study-held-out VDJdb (`vdjdb-2026-11-06`) across Minervina 2022 SARS-CoV-2 (`7,317` paired clones) and 10x Healthy Donor 1 canonical viral multimer clones (`246` expanded clones).

5. **TCR Generation Probability & Post-Selection Publicness (`h2h_pgen_vs_olga.csv`)**:
   - **Comparators**: **`OLGA` `v1.3.0`** (Sethna et al., *Bioinformatics* 2019; exact V(D)J recombination `Pgen`) and **`soNNia` `v0.4.0`** (Isacchini et al., *PNAS* 2021; `Sonia(ppost_model="human_T_beta"|"mouse_T_beta")` pretrained post-selection `Ppost = Q * Pgen`).
   - **Cohorts**: `1,500` unique paired α/β human clones across 10x Healthy Donors 1 & 2 and Minervina 2022, and all `1,781` unique `CDR3_beta` clones across 5 inbred C57BL/6 mice (`M1`–`M5`, Khatun 2021 `GSE158896`).

6. **Multiplexed pMHC Dextramer Calling (`h2h_multimer_callers.csv`)**:
   - **Comparators** (on expanded `>=2`-cell paired α/β CD8+ clones in 10x Healthy Donors 1 & 2 across `44` pMHC dextramers + `6` negative controls):
     - `10x_per_cell_threshold`: `UMI >= 10` & `>= 5x max(NC, 1)` single-binder in `>=1` cell.
     - `clr_z3_per_cell`: Centered log-ratio across the `44` pMHCs + `z >= 3.0` & `UMI >= 5` in `>=1` cell.
     - `muon_dsb_per_cell`: **`muon` `v0.2.1`** (`muon.prot.pp.dsb`, Mulè et al., *Nat Commun* 2022) using TCR-negative cells in the same library as background + the `6` `*_NC` multimers as `isotype_controls`, calling single-dominant binders with `DSB >= 8.0` & `UMI >= 10`.
     - `gmm_bimodal_per_cell`: Per-dextramer 2-component Gaussian Mixture Model (`log1p(UMI)` posterior `>= 0.95` + `UMI >= 10` + `> 5x max(NC, 1)`).
     - `tcrsift_clonal_concordance`: `call_multimer_specificity` + `summarize_clone_multimer_specificity` (`UMI >= 10`, `>= 5x max(NC, 1)`, `min_specific_cells=2`, `min_concordance=0.5`).
