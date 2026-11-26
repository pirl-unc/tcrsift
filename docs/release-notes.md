# Release notes

Release changelog for `tcrsift`, including behavior changes in `tcrsift
prioritize`, `tcrsift run`, and construct assembly that can alter results.
Record the package version alongside every analysis; the
[prioritization guide](user-guide/til-signatures.md) documents the current
release. Full release history is on
[GitHub](https://github.com/pirl-unc/tcrsift/releases) and
[PyPI](https://pypi.org/project/tcrsift/#history).

## 3.43.1

Refined native TCRdist matching, repeat-free full-ORF codon optimization, code
simplification, and API/methods documentation:

- **Mixed V-gene database support & zip-safe resource loading (`tcrsift.tcrdist`)**:
  `TCRDistIndex` now retains V-gene-less database rows (common in combined
  VDJdb + IEDB + CEDAR tables) and queries with a missing/unmapped V gene via
  exact trimmed-CDR3 matching rather than dropping exact hits, while continuing
  to require mapped V-gene loops for all fuzzy (`d > 0`) matches. Switched
  packaged `vgene_cdr_loops.csv.gz` loading to zip-safe `ref.open("rb")`.
- **Distance-aware `db_match_score` for `match_mode="tcrdist"`**:
  `_compute_match_score` now decays smoothly with TCRdist (`max(0.4, 1 - 0.5 *
  distance / scale)`, `scale = 48` paired / `24` β-only) instead of applying a
  flat `0.7x` penalty across all `*_near` hits, while keeping `exact` and
  `levenshtein` scores unchanged.
- **Repeat-free `codon_optimization="full"` & shared prefix optimizer**:
  Added direct-repeat avoidance (`max_direct_repeat=20`, `avoid_kmers`,
  `max_codon_repeats`) to `optimize_codons` and unified the duplicated
  `focal`/`full` prefix optimization logic in `construct.py` and
  `leader_policy.py` into `optimize_or_recode_prefix`. Optimizing β before α
  and passing β k-mers into α eliminates cross-chain ≥20 bp FR2/FR4 repeats
  (`46/46 = 100%` repeat-free, max `17 bp` vs `19 bp` in `DNAChisel`).
- **Documentation updates**: Added `tcrsift.tcrdist` API docs to
  `docs/api/annotate.md`, codon-optimization docs to `docs/api/assemble.md`,
  and refreshed `docs/how-it-works.md`, `docs/methods-section.md`, and
  `docs/benchmarks-overview.md`.

## 3.43.0

Added native, zero-dependency **TCRdist** (`--match-mode tcrdist`) similarity
scoring and public TCR database matching:

- **Native `tcrsift.tcrdist` module (pure Python + NumPy, human & mouse)**:
  Implements the Dash (*Nature* 2017) / `tcrdist3` (*eLife* 2021) distance
  metric — IMGT V-gene CDR1/CDR2/CDR2.5 loop distance (weight 1) plus
  BLOSUM62-weighted trimmed CDR3 single-gap alignment (weight 3, `ntrim=3`,
  `ctrim=2`, `gap_penalty=4`) — backed by a compact 4.3 KB bundled V-gene loop
  table (`tcrsift/refseqs/vgene_cdr_loops.csv.gz`, 448 canonicalized human and
  mouse `TRAV`/`TRBV` genes) and an exact admissible-bound branch-and-bound
  index (`TCRDistIndex`). Exports `cdr3_distance`, `v_gene_distance`,
  `tcr_distance`, and `TCRDistIndex` lazily from `tcrsift`.
- **Database matching via `--match-mode tcrdist`**: Wired into
  `match_clonotypes`, `annotate_clonotypes`, `AnnotateConfig`, `tcrsift annotate`,
  and `tcrsift run`. Supports `strict_ab`, `ab_with_partial`, and `b_only`
  strictness, allows `--max-distance 0` (exact V-gene loops + trimmed CDR3),
  auto-defaults `--max-distance` to `48` (paired) or `12` (`b_only`), and
  gracefully degrades to BLOSUM62 trimmed-CDR3 distance when V-gene columns are
  absent.
- **Head-to-head parity & 15.8x speedup vs `tcrdist3 0.3`**: On study-held-out
  VDJdb (`vdjdb-2026-11-06`), native `tcrsift --match-mode tcrdist` (`r <= 48`)
  recovers **567/607 (93.4% precision, 13.5% recall)** on Minervina 2022
  verified SARS-CoV-2 clones ($n=4{,}212$) in **12.9 s** end-to-end (vs
  557/596 = 93.5% in **204.1 s** for `tcrdist3 0.3`, **15.8x faster**, with no
  `parasail`/`tcrdist3` dependency) and **143/145 (98.6% precision, 58.1%
  recall)** on 10x Donor 1 canonical viral clones ($n=246$) in **8.0 s**
  (**115/115 = 100.0%** at `r <= 24`).

## 3.42.0

Fixed four assembly and murine `Pgen` bugs discovered during a cross-tool audit
and added `tcrdist3 0.3` to the head-to-head benchmark suite:

- **Codon-aligned germline V-gene leaders (`TRAV41*01` frameshift fix)**:
  Trimmed trailing 1–2 nt IMGT exon-boundary overhangs in
  `scripts/build_vgene_leaders.py`, `tcrsift/refseqs/vgene_leaders.csv.gz`, and
  `tcrsift/assemble/leaders.py` (14 alleles including functional `TRAV41*01`,
  previously 67 nt = `3*22 + 1`), preventing downstream reading-frame shifts
  when germline leader fallback is used.
- **Validated α/β J→C junction seams & Cell Ranger 3 `TRAJ37` boundary rescue**:
  `_add_constant_regions` now validates contig junction codons against
  `_VALID_JUNCTION_RESIDUES` (`{N, Y, D, H}` for α, `{E}` for β), rejecting
  unspliced intronic C-tail reads in favor of germline `_traj_junction_residue`,
  and cleanly shifts over-extended Cell Ranger 3 `TRAJ37` junction residues
  (`...LQVKPD` + bare `ATCCAGAAC...`) so `alpha_constant_aa` starts with
  `DIQNPDPAV` rather than `IIQNPDPAV`. Non-codon 1–2 nt stubs past VDJ are no
  longer written to `*_constant_nt_contig` or flagged as load-bearing frame
  errors. On Fischer 2021 (160 5'-truncated clones), load-bearing validation
  errors dropped from **29 to 0**.
- **Bundled murine k-mer `Pgen` models (`kmer_pgen_mouse_{alpha,beta}.npz`)**:
  Generated offline from 500,000 `OLGA` `mouse_T_{alpha,beta}` draws (~30 KB
  each, Apache-2.0 NumPy tables only). On Khatun 2021 murine GP66 CD4 TCRs
  ($n=1{,}781$), correlation with `OLGA mouse_T_beta` jumped from $r = 0.198$
  to **$r = 0.910$** ($\rho = 0.904$, exceeding pretrained `soNNia mouse_T_beta`
  at $r = 0.877$) while running **266x faster than `OLGA`** and **846x faster
  than `soNNia`**.
- **Head-to-head benchmark audit & `tcrdist3 0.3` comparator**: Added
  `tcrdist3 0.3` (via isolated `tcrdist_worker.py` + `Dockerfile.tcrdist`),
  aligned Stage 2a Minervina evaluation to the exact 4,212 dextramer-verified
  SARS-CoV-2 clones with uniform modal-species aggregation (`tcrdist3` paired
  `r <= 48` reaches 557/596 = 93.5% precision, 13.2% recall vs 210/226 = 92.9%
  for `tcrsift` paired Lev-1 at 14.9s), calibrated `muon` DSB to the standard
  `>= 3.0` cutoff (101/103 = 98.1% precision, 369 sticky `KLGGALQAK` FPs) with
  `n_distinct_epitopes_called`, and added `n_repeat_lt20` / `max_repeat_bp` to
  the `DNAChisel` codon-optimization table.

## 3.41.0

Added paired αβ Levenshtein-1 database matching (`ab_near`) and expanded the
head-to-head benchmark suite with five additional open-source Python tools:

- **Paired αβ Levenshtein-1 matching (`ab_near`)**: `match_clonotypes(...,
  match_strictness="strict_ab", match_mode="levenshtein")` now performs O(L)
  deletion-canonical paired αβ neighbor search (matching paired database
  entries within edit distance ≤ 1 on each chain, `db_match_strength="ab_near"`),
  rather than silently acting as an exact-only no-op under `strict_ab`. On
  study-held-out VDJdb, `strict_ab` Levenshtein-1 matches the exact same 292
  Minervina clones and 110 10x Donor 1 clones (99.1% cognate precision, 2x
  recall over exact) as `scirpy 0.26.1` while running **31.7x faster** (10.2s vs
  323.9s on 7,317 clones).
- **Five new external open-source comparators in `benchmarks/head_to_head`**:
  Added **`DNAChisel 3.2.16`** (Stage 1b synthesis codon optimization),
  **`scirpy 0.26.1`** (Stage 2a held-out VDJdb matching), **`muon 0.2.1` DSB**
  (Stage 2b pMHC multimer background denoising), **`soNNia 0.4.0`** (Stage 2c
  pretrained `Sonia` post-selection `Ppost`), and **`pyucell 0.8.0`** (Stage 3
  rank-based UCell signature scoring) alongside `stitchr`, `OLGA`, and `Scanpy`
  with two new SVG figures (`h2h_db_matching_scirpy.svg` and
  `h2h_codon_opt_dnachisel.svg`).

## 3.40.1

Regenerated all five benchmark suites on `3.40.1` and restructured all five
benchmark documentation pages:

- **High-level results & dataset motivation on every benchmark page**: Added
  scannable **High-level / At-a-glance results** scorecard tables and explicit
  **Why these datasets were included** motivation tables to
  `docs/benchmarks-overview.md`, `docs/benchmarks.md`,
  `docs/benchmarks-infectious.md`, `docs/benchmarks-blood-cd4.md`, and
  `docs/benchmarks-mouse.md`.
- **Full 5-suite benchmark & figure regeneration**: Re-ran `public_til`,
  `public_infectious`, `public_blood_cd4`, `public_mouse`, and `head_to_head`
  (plus all 6 SVG figure generators) on `3.40.1` so every benchmark CSV,
  `summary.json`, and documentation figure is synchronized on the current
  release.

## 3.40.0

Enabled automatic truncation detection and IMGT germline VDJ fallback by default
across construct assembly:

- **`germline_vdj_fallback=True` by default**: `assemble_full_sequences()`,
  `AssembleConfig`, `PrioritizeConfig`, `tcrsift assemble`, `tcrsift run`,
  `tcrsift prioritize`, and `tcrsift report selected` now automatically detect
  per-chain truncation (`partial_VDJ_*`) or bare `(v_gene, CDR3, j_gene)` rows
  and rescue them against cached IMGT V-REGION and J-REGION references.
- **Intact donor contigs untouched**: Full-length donor VDJ sequences
  (`VDJ_{chain}_aa` present) remain a complete no-op (`{chain}_vdj_source =
  "contig_full"`), preserving non-`*01` donor V-region alleles verbatim and
  never injecting back-translated codons into AA-only intact inputs. Every
  germline-imputed chain continues to stamp `{chain}_vdj_source`,
  `{chain}_vdj_germline_allele`, `qc_warnings`, and
  `construct_contig_verified=False`.
- **Fail-closed strict opt-out**: Pass `--no-germline-vdj-fallback` (or
  `germline_vdj_fallback=False`) to disable germline VDJ imputation and require
  5'-complete donor contigs. When truncated or missing chains cannot be rescued
  because IMGT references are uncached, assembly logs a warning pointing to
  `tcrsift data download --group germline` (or `germline_mouse`) /
  `--auto-download-germline`.
- **Head-to-head assembly benchmark updated**: Regenerated
  `benchmarks/head_to_head` (`h2h_assembly_vs_stitchr.csv` and
  `h2h_assembly_stitchr.svg`) showing default `tcrsift assemble` achieving 100%
  completion across Minervina CR6 (46/46, 9 donor V-variants kept), Fischer
  5'-truncated (160/160, 32 partial 5' contigs rescued), and Mouse Cell Ranger
  (101/101, 52 `Trbc` miscalls fixed) alongside `--no-germline-vdj-fallback`
  strict mode.

## 3.39.2

Cleaned up benchmark and prioritization documentation:

- Removed the out-of-place CLI cheatsheet wall from `docs/benchmarks-overview.md`
  so the benchmark page leads straight into Stage 1 (Construct Assembly) →
  Stage 2 (Identifying Expected TCRs) → Stage 3 (Expression Signatures).
- Rewrote `docs/user-guide/til-signatures.md` examples and context-preset table
  in plain, concise prose focused on the 5 primary contexts (`solid-tumor`,
  `blood-tumor`, `peptide-stim`, `blood`, `generic`).

## 3.39.1

Empirically tightened context presets, simplified mode guide, and regenerated
head-to-head benchmarks:

- **Removed unvalidated `post-act` context**: Empirical inspection of Chiffelle
  Day-30 post-ACT (`T30`) biopsies confirmed that `AntigenExperienced` (0.33),
  `Differentiated` (0.45), and `Cytolytic` (0.48) are all below chance after
  rapid TIL expansion and perform worse than `TumorReactive` (0.61). Removed
  `post-act` so every primary context in `PRIMARY_CONTEXTS` (`solid-tumor`,
  `blood-tumor`, `blood`, `peptide-stim`, `generic`) is empirically supported.
- **Simplified mode-selection guide**: Replaced verbose multi-column taxonomy in
  `docs/benchmarks-overview.md` and `docs/user-guide/til-signatures.md` with a
  concise 4-step decision rule focused strictly on choices that change what
  TCRsift computes.
- **Updated head-to-head benchmark labels & outputs**: Renamed `tcrsift_default`
  in `benchmarks/head_to_head/run.py`, `PROVENANCE.md`, and
  `h2h_prioritize_baselines.csv` to explicit preset names
  (`tcrsift_solid_tumor`, `tcrsift_blood_tumor`,
  `tcrsift_solid_tumor_cd4_excl_treg`) and regenerated `benchmarks/head_to_head`
  results and SVGs on `3.39.1`.

## 3.39.0

Explicit biological-context & assay ontology, new `post-act` and `peptide-stim`
prioritization contexts, and restructured benchmark overview:

- **Context & assay ontology (`tcrsift.PRIMARY_CONTEXTS`, `tcrsift.CONTEXT_ONTOLOGY`)**:
  Organized `tcrsift prioritize --context` into 6 orthogonal primary biological
  regimes — `solid-tumor` (fresh TIL), `post-act` (new: post-ACT / infused TIL
  biopsies where chronic *CXCL13*/*CD39* residency is attenuated; Tier 1
  `AntigenExperienced`, `Cytolytic`, `Differentiated`), `blood-tumor` (cancer
  blood), `blood` (infection/vaccine/healthy/heme circulation), `peptide-stim`
  (new: 6–24 h unpaired restimulation; Tier 1 `AIM`, `AcuteActivation`,
  `Proliferation`, `Cytolytic`; paired-control restimulation remains
  `tcrsift run --method stimulation` and multi-day culture remains
  `tcrsift run --method threshold|logistic`), and `generic` (6 core state panels
  only) — while keeping `heme` (alias for `blood`) and `mpe` (non-epithelial
  effusion variant) for backward compatibility.
- **CD4 TIL Treg warning**: `tcrsift prioritize` now logs an explicit warning
  recommending `--exclude-treg` when `FOXP3+` Treg-like clones are detected
  under `--tcell-type cd4`/`both` in `solid-tumor` or `mpe` contexts with
  `--exclude-treg` disabled.
- **Restructured benchmark overview & clearer figure labels**: Reordered
  `docs/benchmarks-overview.md` from objective ground truth outward —
  **Stage 1: Construct assembly correctness → Stage 2: Identifying expected
  antigen-specific TCRs (VDJdb `strict_ab`, pMHC multimers, `Pgen`/`Ppost` vs
  `OLGA`, 16 h restimulation) → Stage 3: Context-dependent expression signatures
  & selection heuristics** — and updated the prioritization table and SVG
  figures (`h2h_prioritize_auroc.svg`, `h2h_prioritize_cd4_treg.svg`) to label
  exact `--context` presets rather than `"tcrsift_default"` and clarify that
  `NeoTCR8`/`NeoTCR4`/`NeoTCR_PBL`/`MANAscore` are built-in `tcrsift`
  signatures.

## 3.38.1

Regenerated all five empirical benchmark suites and all documentation figures on
the post-`3.38.0` codebase:

- **Full benchmark verification**: Re-ran `benchmarks/public_til/`,
  `benchmarks/public_infectious/`, `benchmarks/public_blood_cd4/`,
  `benchmarks/public_mouse/`, and `benchmarks/head_to_head/` end-to-end against
  VDJdb `2026-11-06`, IMGT human/mouse references, `stitchr/thimble v1.3.1`,
  `OLGA v1.3.0`, and `Scanpy 1.12.3`. Every scientific metric (AUROCs, recall,
  precision, clone/shortlist counts, assembly completeness, donor polymorphism
  preservation, Pgen/Ppost correlation, and multimer concordance) matched
  `3.37.1` identically.
- **Figures & runner fix**: Regenerated all SVG figures via all six scripts in
  `docs/figures/` (including tightened `bbox_inches="tight"` bounds on the five
  methods/schematic SVGs from `make_figures.py`), and fixed
  `benchmarks/public_til/run.py` to count entries from the relational VDJdb
  layout (`chains.tsv.gz` + `records.tsv.gz`) via `load_vdjdb`.

## 3.38.0

Dead-code removal, internal helper deduplication, and `til-select` cleanup:

- **Dead-code removal**: Removed unreferenced internal functions and unused
  constants across `annotate.py` (`_pick_vdjdb_file`), `pgen_models.py`
  (`_observed_seqs`), `til.py` (`_clean_cdr3_col`), `assemble/constants.py`
  (`_count_trbc_overrides`), `assemble/vregion.py`
  (`jregion_reference_is_available`), `validation.py`
  (`validate_cellranger_vdj_dir`, `validate_cellranger_gex_dir`), `genes.py`
  (`CD3_GENES`, `CD8_GENES`, `local_pyensembl_releases`), `mnemonic.py`
  (unused vowel/cluster tables and unreachable type check), `pgen.py` (dead
  `_MOUSE_TRAV_ALIASES` lookup inside `_j_usage_log`), and `unify.py` (unused
  tracking lists).
- **Helper deduplication**: Unified reverse-complement in `assemble/qc.py`
  onto `assemble.codons._revcomp`; unified IMGT functionality ranking in
  `assemble/germline.py` onto `assemble.vregion._func_rank_key`; replaced
  duplicate boolean/mode wrappers in `til_select.py` with
  `_schema.truthy_flag_series` and `validation.most_common`; simplified
  `til_prioritize._ensure_clone_key` onto `_schema.paired_cdr3ab_key`; and
  removed redundant `ImportError` guards for guaranteed core dependencies
  (`reportlab`, `pyyaml`, `matplotlib`).
- **`til-select` output & CLI cleanup**: Stopped emitting 14 empty placeholder
  PNGs into `fig_dir` and `all_plots.pdf`, and hid legacy no-op compatibility
  options (`--trend-increase-ratio-min`, `--trend-decrease-ratio-max`,
  `--pyensembl-release`) from `tcrsift til-select --help` while preserving CLI
  backward compatibility.

## 3.37.5

Core I/O, QC, GEX, multimer, report, `til-select`, and plotting bugfixes:

- **SCT & GEX**: Fixed `sct.aggregate_sct` inverted boolean aggregation
  (`any=False, all=True`) by normalizing with `truthy_flag_series`; gated CP10K
  normalization in `gex.compute_signature_scores_per_cell` on the full matrix
  rather than the queried signature subset; resolved `GENE_SYNONYMS` pairs via
  `sym_to_var` in `gex.augment_with_gex`; guarded empty DataFrames in
  `gex.compute_cd4_cd8_counts`; and guarded 0-cell sparse `ZeroDivisionError`
  and `NaN` column selection in `genes.resolve_adata_var_map`.
- **QC, attribution & h5ad serialization**: Excluded `NaN` secondary-chain
  slots (`TRA_2_*`, `TRB_2_*`) from falsely failing `qc.validate_clonotypes`,
  fixed `Categorical` `.fillna("")` crashes in `qc.py` and `attribution.py`,
  added fallback scoring for small/targeted panels in `cross_lineage_doublets`,
  and stringified mixed-type object columns in `_dtypes.rehydrate_obs` so
  `write_h5ad` never fails on mixed-scalar `obs` columns.
- **Multimers, reports & datacache**: Fixed `_DEFAULT_NC_RE` to recognize
  `NC_1_TotalSeqC`, `dCODE_NC_1_PE`, and `NC_ctrl` while still excluding
  NCBI accessions; implemented sample-prefix stripping (`S1_ACGT-1` → `ACGT-1`)
  in `call_multimers_for_adata`; wrapped `_image_reader` headers in `pdf_safe`,
  validated `extra_images` existence in `bundle_figure_pdf`, and emitted
  `selected_clone` on the early-return path of `expand_dual_alpha_variants`;
  converted `zipfile.BadZipFile` into `DownloadError` in `datacache._extract_zip`.
- **CLI, `til-select` & plots**: Excluded overlapping columns when merging
  `--clone-scores` in `tcrsift select --prism`; allowed `tcrsift annotate-gex`
  to fall back to `CDR3ab` or skip signature aggregation when `CDR3_pair` is
  absent; stopped `--exhaustion-genes` default from masking `--enrichment-genes`;
  fixed `dominant_specificity_from_master` double-negative lineage assignment
  and deduplication, `select_top_monotonic_increase_clones(..., require_positive_all_timepoints=True)`,
  `plot_annotated_heatmap` annotation tags, duplicate-`CDR3ab` crashes in
  `plots/ranks.py` and `plots/annotations.py`, empty/missing-`cell_count`
  handling in `plot_filter`, and `'<NA>'` leakage in `plot_provenance`.

## 3.37.4

Assembly and pipeline correctness fixes:

- **AA-only VDJ nucleotide assembly**: Back-translate missing `vdj_*_nt` /
  `*_leader_nt` segments in `_build_full_sequences` when AA sequences are
  present so `full_*_nt` and `single_chain_nt` never splice `leader_nt` +
  `constant_nt` directly while omitting the V(D)J region.
- **Mouse germline leader cache forwarding**: Propagated `data_dir` /
  `--cache-dir` through `_assemble_clone`, `_extract_leader_from_contigs_single`,
  `_record_leader`, and `apply_leader_policy` so custom-cached murine IMGT
  leader references resolve cleanly.
- **2A linker enzyme avoidance**: Recoded `P2A`/`F2A` linker DNA against
  `avoid_enzymes` (including β–linker and linker–α seams) in
  `single_chain_nt_optimized`, eliminating the reverse-strand `BsmBI` site.
- **Assembly QC & fidelity gates**: Made `validate_sequences` and
  `build_assembly_qc_report` recognize VDJ-only (`include_constant=False`)
  assemblies without failing 200-aa / TRAC/TRBC constant-tail checks; reported
  `from-data` constants accurately in `enforce_contig_fidelity` and skipped
  unpaired chains in `_add_constant_from_row`; flagged length discrepancies in
  `germline_compare_vregion` instead of labeling truncated germlines
  `"identical"`; and fixed non-unique index and all-NaN `float64` assignment
  errors in `apply_leader_policy`.
- **Pipeline (`run` & `prioritize --assemble`)**: Applied
  `--min-til-cells-per-donor` after `match_til` populates
  `max_til_cells_per_donor` in Step 6, resolved relative sample-sheet `vdj_dir`
  paths against the sample sheet's directory, and guarded stimulation scoring in
  `run_pipeline` against double-`normalize_total`/`log1p` on pre-logged `h5ad`
  inputs.

## 3.37.3

Dependency hardening, Pandas 3 / `harmonypy 2.0` compatibility, and optional
extra guards:

- **Packaging (`pyproject.toml`)**: Declared direct core dependencies `scipy`
  and `h5py`, added `<3.0.0` upper bound on `pandas`, removed unused `leidenalg`
  from `[atlas]` and `[full]` (Leiden uses `flavor="igraph"`), and updated the
  `openpyxl` rationale comment and installation docs.
- **Pandas 3 / PyArrow `StringDtype` compatibility**: Rewrote
  `annotate._match_viral_species_lower` with RE2-compatible word boundaries
  (`\bflu\b`, `\binfluenza\b`) so PyArrow-backed string columns no longer raise
  `ArrowInvalid`, and eliminated `combine_first` empty-series `FutureWarning`s
  in `_load_iedb_receptor_export`.
- **Atlas embedding & `harmonypy 2.0.0`**: Guarded `nclust >= 2` and
  `verbose=False` in `embed_cells` so `harmonypy 2.0.0` no longer crashes when
  `n_obs < 60`, added actionable `ImportError` hints pointing to
  `pip install 'tcrsift[atlas]'` for `harmonypy` and `igraph` (including in
  `tcrsift cells embed`), narrowed Pearson-residual warning suppression to
  `RuntimeWarning`, and used `np.nanmean` in `gd_cd8_relabel` so `pd.NA` in
  `ab_col` no longer mislabels γδ T sub-clusters as CD8.
- **Pgen / plotting optional guards**: Guarded optional `tcrpeg` import in
  `pgen_models._bundled_observed_beta`, filtered `NaN` CDR3 rows before reading
  V/J genes in `pgen_models.train_model`, and corrected the `matplotlib-venn` /
  `tcrsift[plots]` hint in `plots/selection.py`.

## 3.37.2

Bug-fix release addressing 25 correctness, safety-control, CLI/config, plotting,
and edge-case issues across QC, assembly, prioritization, annotation, filtering,
and documentation:

- **QC battery & species inference**: Replaced `qc_battery._is_mouse_row` with
  public `infer_row_vdj_species(row) -> "human" | "mouse"` using species-unique
  TRAC/TRBC tails (`MGLRILLLKVAGFNLLMTLRLWSS` / `VKRKNS` / `VKKKNS`) plus
  explicit species and chimeric constant-source metadata so human constructs no
  longer fall into mouse J-anchor rules. Updated `check_assembly` to tolerate
  CellRanger v7+'s standard +1 nt J→C splice-donor overhang on raw `VDJ_*_nt`
  columns while keeping assembled `vdj_*_nt`, `full_*_nt`, and
  `single_chain_nt` strictly in-frame.
- **Assembly & codon optimization**: Extended `optimize_codons` with bounded
  multi-codon lookahead so 3-codon restriction sites (such as `NcoI` `CCATGG`
  across `T-M-E` in mouse `TRAC`) are eliminated, and propagated
  `forbidden_motifs` into explicit-allele override and chimeric constant
  branches in `_add_constant_regions`.
- **Prioritization & pipeline safety controls**: Added
  `allow_canonical_fallback` to `PrioritizeConfig` so
  `--no-allow-canonical-fallback` is preserved; passed `clone_scores=til_matched`
  to `select_from_clone_sample_long` in `run_pipeline` so `percentile_rank`
  selection rules rank by actual scores instead of `inf`; fixed
  `select_round_robin` non-unique index corruption; normalized `tcell_type` case
  in `resolve_signatures`.
- **Cell & TCR annotation / signatures**: Resolved `gene_symbols` in
  `annotate_cells.top_markers` on Ensembl-indexed AnnData so clonal TCR/Ig
  segments are excluded and cluster suffixes emit readable gene symbols; ignored
  blank `antigen_gene` / `db_protein_canonical` strings in
  `annotate._annotate_match`; aligned `score_by_name` `zscore`/`weighted_z`
  implicit-`log1p` gating with `geneset_enrichment` to avoid double-logging
  log-normalized matrices lacking `uns["log1p"]`.
- **CLI, config, and validation**: Added `--sample-sheet` support to
  `tcrsift assemble` and `tcrsift report selected` (and accepted file paths in
  `load_contigs(..., sample_name_from="sheet")`), enabled sample-sheet `vdj_dir`
  auto-discovery for `--sample-name-from sheet` in `run_pipeline`, cleared the
  mutually exclusive `contigs_dir` / `cellranger_dir` counterpart on CLI
  override, and allowed sample-sheet `vdj_dir`s in `validate_run_args`.
- **Edge-case hardening & plots**: Fixed empty/β-only DB crash in
  `til_select.match_clonotypes`, duplicate-index crash in
  `filter_clonotypes_stimulation`, raw `KeyError` on β-only tables in
  `filter_clonotypes_threshold`, zero-division on all-doublet input in
  `aggregate_clonotypes`, empty-frame `KeyError` in `calculate_clone_frequencies`,
  `NaN` clone-ID `IndexError` in `select_freq_prism_per_condition`, missing
  `"sample"` crash in `identify_til_specific_clones`, 0-cell crashes in
  `phenotype.py`, `publicness_col` coloring in `plot_pgen_distribution`,
  duplicate-`CDR3ab` crashes in `plots/ranks.py`, missing-`CDR3ab` guard in
  `plot_vgene_usage_by_method`, and stale `reports`/`pdfkit` extra docs in
  `README.md` and `docs/getting-started/installation.md`.

## 3.37.1

Re-ran all five empirical benchmark suites (`public_til`, `public_infectious`,
`public_blood_cd4`, `public_mouse`, and `head_to_head`) end-to-end on the
post-3.35/3.37 codebase and refreshed all benchmark figures and documentation:

- All TIL (`public_til`), blood/post-ACT/CD4 (`public_blood_cd4`), murine
  (`public_mouse`), and head-to-head prioritization / multimer-caller CSVs
  reproduced identically to prior releases.
- Updated Fischer 2021 5'-partial contig rescue counts (`33 → 32` rescued
  partial contigs, `25 → 23` J→C verified, `160/160` complete constructs)
  reflecting the `3.35.0` FR4-tail guard, and updated 10x multimer summary
  counts (`multimer_calling.csv`: `1,019 → 39` sticky `KLGGALQAK` clones, `0`
  confident HIV Gag false positives) reflecting deterministic tie-breaking.
- Fixed nested `vdjdb.txt` discovery in `benchmarks/public_til/run.py` so
  VDJdb provenance stamps reliably into `summary.json`.

## 3.37.0

Dead-code removal, assembly subpackage consolidation, and CLI/config/CI cleanup:

- **Dead & duplicated code removed:**
  - Removed extinct pre-v1 `tcrsift/data.py` (`annotate_combined_df`) and unused
    legacy plotting/threshold helpers (and eager `matplotlib` import) from
    `tcrsift/model.py`, retaining `LogisticFrequencyModel`,
    `fit_frequency_logistic_model`, and `analytical_logistic_threshold`.
  - Deduplicated `tcrsift/til_select.py` IEDB/CEDAR loaders by delegating
    directly to `tcrsift.annotate.load_iedb` and `tcrsift.annotate.load_cedar`.
  - Renamed `tcrsift.qc_battery.QCResult` to `BatteryCheckResult` to eliminate
    the class-name collision with `tcrsift.qc.QCResult`.
- **Assembly subpackage consolidation:**
  - Moved `leaders.py`, `sigpep.py`, and `vregion.py` into `tcrsift/assemble/`
    (preserving all public re-exports on `tcrsift.assemble` and maintaining a
    strict acyclic import DAG).
  - Moved offline reference-regeneration helper `_regen_canonical_constants.py`
    out of the installed wheel into `scripts/regen_canonical_constants.py`.
- **CLI, config & CI polish:**
  - `tcrsift report` and `tcrsift pgen` without a subcommand now exit non-zero
    (`1`) instead of silently succeeding with `0`.
  - `TCRsiftConfig._from_dict` treats empty YAML sections (`load:`) as defaults
    and raises `TCRsiftValidationError` (instead of raw `TypeError`) for unknown
    nested keys; `cmd_generate_config` streams YAML directly to stdout without a
    temp file; `tcrsift prioritize --assemble` now exposes
    `--allow-canonical-fallback` / `--no-allow-canonical-fallback`; `cells qc`
    displays actual default thresholds in `--help`; `data download` respects
    `--verbose` instead of forcing debug logging.
  - Added a `build` (`python -m build` + `twine check`) CI job to
    `.github/workflows/tests.yml`, removed redundant `ruff` runs across test
    matrix legs, and dropped `torch`/`tcrpeg` from the `lint` job; removed
    obsolete `pylintrc` and `requirements-dev.txt` and aligned `CONTRIBUTING.md`
    and `RELEASING.md`.

## 3.36.1

Targeted bug fixes, documentation rendering fixes, and benchmark figure polish
from the holistic post-3.36.0 repository audit:

- **Top-level API & CLI fixes:**
  - Fixed `tcrsift.match_clonotypes` to export `tcrsift.annotate.match_clonotypes`
    rather than the stripped-down legacy copy in `tcrsift.til_select`.
  - Wired `apply_cd3_gate` into standalone `tcrsift phenotype` so it matches
    `tcrsift run` and its CLI documentation.
  - Normalized string `"none"` to `None` in `AssembleConfig.__post_init__` so
    YAML `alpha_leader: none` / `beta_leader: none` no longer crashes at Step 7.
  - Registered `--germline-vdj-fallback`, `--auto-download-germline`, and
    `--cache-dir` on `tcrsift report selected`; fixed `validate_run_args` and
    `validate_assemble_args` when `--no-assemble` / `--no-leaders` or YAML
    `from_contig` settings are used.
  - Fixed `reconstruct_germline_vdj_chain(..., species="mouse")` to inherit
    `codon_species` from `species` when unspecified; moved `tcrsift.mnemonic`'s
    module docstring ahead of `from __future__ import annotations`; added
    `python -m tcrsift` (`tcrsift/__main__.py`) and `oncoref` to run provenance.
- **Scripts, packaging & CI:**
  - Fixed `deploy.sh --dry-run <version>` so it no longer mutates `version.py`,
    commits, or pushes during dry runs; guarded `set -e` in `develop.sh` so
    sourcing it cannot kill interactive shells; aligned `lint.sh` targets with CI.
  - Scoped `.github/workflows/docs.yml` concurrency per-PR so PR doc builds
    cannot cancel live `main` GitHub Pages deploys.
  - Cleaned `pyproject.toml` extras (dropped unused `pdfkit` and duplicate
    `reportlab`, added `matplotlib-venn` to `plots`/`full`, aligned docs URL).
- **Docs & benchmark figures:**
  - Enabled `pymdownx.arithmatex` + MathJax 3 in `mkdocs.yml` so equations in
    `methods-section.md` and `benchmarks-overview.md` render properly.
  - Updated `docs/methods-section.md`, `docs/user-guide/atlas.md`,
    `docs/api/insilico_filter.md`, `docs/api/signatures.md`, and
    `docs/benchmarks-overview.md` to match current `3.36` code behavior.
  - Regenerated benchmark SVGs to eliminate legend-over-bar collisions, label
    overlaps on `benchmark_infectious_multimer_qc.svg`, and right-edge clipping.

## 3.36.0

Internal package organization cleanup (no CLI flag, scoring, or selection
behavior changes):

- Converted `tcrsift.cli`, `tcrsift.assemble`, and `tcrsift.plots` from flat
  modules plus underscore-prefixed helper files (`_cli_commands.py`,
  `_assemble_*.py`, `_plots_*.py`) into cohesive domain subpackages:
  - `tcrsift.cli` (`pipeline`, `load`, `cells`, `filter`, `til`, `annotate`,
    `assemble`, `report`, `data`, `common`, `__main__`), pairing parser
    registration and command handlers per workflow domain.
  - `tcrsift.assemble` (`construct`, `contigs`, `leader_policy`, `constants`,
    `codons`, `germline`, `qc`, `alleles`).
  - `tcrsift.plots` (`state`, `repertoire`, `annotations`, `ranks`, `atlas`,
    `selection`, `reports`).
- Preserved all existing public re-exports on `tcrsift.cli`, `tcrsift.assemble`,
  and `tcrsift.plots`, plus `python -m tcrsift.cli`.

## 3.35.0

Whole-library bug-fix release addressing correctness, statistical edge cases,
assembly QC, and I/O safety across all pipeline stages, plus a standalone
manuscript Methods section ([methods-section.md](methods-section.md)):

- **Loading, QC, phenotyping, and clonotyping.** Guarded mitochondrial QC
  (`min_mito_pct` / `max_mito_pct`) when no mitochondrial genes are present in
  the expression matrix; made marker extraction and GEX signature scoring safe
  under duplicate `var_names` and already-log-transformed matrices; fixed
  `species="auto"` gene-ID resolution; removed mouse V-genes (`TRAV12-2`,
  `TRBV28`) from the human-only species detector; enforced numeric
  `min_cd3_reads` in `apply_cd3_gate` and synthesized `CD3` from
  `CD3D`/`CD3E`/`CD3G`; added an expression-ratio guard (`max_imbalance_ratio`)
  to `cd4_cd8_doublet_mask`; fixed multi-sample contig ID collision tracking in
  `qc_battery`; and handled string booleans, missing chains, and per-donor
  method recovery cleanly.
- **Scoring, filtering, selection, TIL matching, and attribution.** Replaced
  grid-clamped logistic FDR thresholds with analytical thresholds
  (`analytical_logistic_threshold`) and non-misleading fallback tiers; fixed
  `subset_percentile` to rank `-hypergeom.logsf` with `method="max"` (avoiding
  float64 tail underflow and midpoint tie bias); tightened `RECEPTOR_GENE_RE`
  and added receptor-gene anti-leakage checks across signature scorers; fixed
  `fdr_scope="auto"` multi-donor cohort detection in `run_pipeline`; moved
  `config.selection` outside the assembly guard so `--no-assemble` still runs
  selection; propagated `patient_id` in `load_til_samples`; folded dual-alpha
  merge counts into attribution priors; and guarded constant score columns
  against passing 100% of clones in `insilico_filter`.
- **Database annotation, probability models, and cache safety.** Fixed IEDB /
  CEDAR UniProt tumor-self classification (MART-1, HER2, PRAME, gp100, TRP-1,
  Tyrosinase, NY-BR-1) and word-bounded influenza/bacterial matching so
  *Haemophilus influenzae* is not misclassified as viral influenza; added Zip
  Slip validation, thread-safe IPv4 fallback, and atomic `force=True` staging
  in `datacache`; fixed log-space insertion probabilities and `pd.NA` handling
  in `pgen`; set `train_model` default k-mer order to `2` matching `seqprob`;
  and fixed `unify` (`Likely_CD8` symmetry, dotted condition names, tie
  handling) and `cohort` empty-beta overlap.
- **Construct assembly and reporting.** Fixed premature-stop detection in
  `_nt_contig_in_frame_no_premature_stop`, bare-mature 1-aa junction offset
  validation, and uppercase mouse `TRAJ` detection in `_assemble_qc`; added
  `NcoI`/`BsaI`/`BsmBI` to forbidden synthesis motifs and guarded recoding
  against introducing restriction sites; skipped leader/constant assembly on
  missing chains of unpaired clones; fixed `TRBC01`/`TRBC02` allele parsing
  and human `TRBC2*01` regeneration; made `tcrsift.annotate_cells` and
  `tcrsift.annotate_tcrs` callable even after submodule import; supported
  CSV-round-tripped stringified `sample_frequencies` in plotting; fixed
  `pretty_sample` on non-marker hyphenated names (`Patient1-TIL`, `MART-1`);
  and preserved multi-condition provenance in `build_selected_report`.

## 3.34.1

Documentation cleanup across the overview, TIL, infectious disease, paired
blood-tumor/CD4, and murine benchmark pages: simplified prose and tables,
removed the `FM1`–`FM8` scorecard layout from the overview page, and split
side-by-side two-panel benchmark plots into individual captioned figures. No
code, scoring, or selection changes.

## 3.34.0

Added `tcrsift prioritize --assemble`, tiered `tcrsift --help` output,
`from_contig` germline leader fallback on 5'-truncated contigs, and the
cross-cohort benchmark overview with external tool comparisons:

- **Construct assembly in `prioritize` / `til-prioritize` (`--assemble`)
  and `from_contig` leader fallback.** Added `--assemble`,
  `--germline-vdj-fallback`, and `--single-chain-linker` (`T2A` default, `P2A`,
  `E2A`, `F2A`, or `none`) to `tcrsift prioritize`, `tcrsift til-prioritize`,
  and `PrioritizeConfig`. When `--assemble` is passed, shortlisted candidates
  are assembled into full-length α, β, and single-chain constructs
  (`candidate_constructs.csv` and `candidate_constructs.fasta`, plus `assembly`
  and `contig_verified_constructs` metrics in `prioritization.json`),
  automatically detecting sample `*contig*.fasta(.gz)` and
  `*contig_annotations*.json(.gz)` files. Also fixed `apply_leader_policy` and
  `_resolve_named_leader` in `tcrsift/assemble.py` so `from_contig` chains whose
  contigs are 5'-truncated or absent receive their V-gene's `M...` germline
  (or species default) signal peptide under `leader_fallback` rather than
  emitting leaderless chains. Shortlist selection logic and default cutoffs are
  unchanged.
- **Tiered `tcrsift --help` CLI organization.** Grouped all 25 top-level
  subcommands into "Primary workflows" (`prioritize`, `run`, `assemble`,
  `data`), "Exploration & annotation" (`cells`, `annotate`, `report`), and
  "Modular pipeline steps & legacy commands" while keeping all existing
  subcommands and arguments unchanged.
- **Cross-cohort overview and head-to-head external benchmarks.** Added
  [Benchmark overview and head-to-head comparisons](benchmarks-overview.md)
  and `benchmarks/head_to_head/`, comparing `tcrsift` against
  `scanpy.tl.score_genes` (`NeoTCR8`, `NeoTCR4`, `NeoTCR_PBL`) and single-gene
  `CXCL13` across 5 human TIL and cancer-blood cohorts, `stitchr` / `thimble`
  1.3.1 across human Cell Ranger 6 donor contigs (9/46 = 19.6% donor
  V-framework variants preserved), 5'-truncated human deposits (160/160 = 100%
  complete `M...` constructs vs 135/160 = 84.4% in `stitchr`), and murine Cell
  Ranger shortlists (101/101 = 100% with 52 `Trbc1`/`Trbc2` miscalls
  corrected vs 78/101 = 77.2% in `stitchr`), `OLGA` 1.3.0 exact V(D)J
  dynamic programming ($r = 0.954$, 451x faster on human blood), and per-cell,
  CLR `z >= 3`, and bimodal GMM pMHC multimer callers (96.3% sticky `KLGGALQAK`
  false-positive reduction and 100.0% [109/109] cognate epitope precision on
  10x-held-out `strict_ab` VDJdb).

## 3.33.0

Empirical validation and real-data fixes across human `--context blood-tumor`
(`NeoTCR_PBL`), longitudinal post-ACT (`T30`) TILs, human `--tcell-type cd4`
(`NeoTCR4`, `pre_sorted`, and `--exclude-treg`), and murine `--species mouse`
workflows:

- **FACS `pre_sorted` coreceptor-dropout rescue (`Sample.pre_sorted`).**
  `load_sample` now propagates `pre_sorted` (`"CD4"` or `"CD8"`, case-insensitive)
  into `adata.obs["pre_sorted"]`, and `phenotype_cells` assigns `CD3`-passing
  cells with zero `CD4` and zero `CD8` UMIs to `"Likely CD4+"` / `"Likely CD8+"`
  on FACS pre-sorted libraries while still rejecting cells expressing the
  opposite coreceptor (`8.7x` cell retention, `708` → `6,144` cells by rescuing
  `88.5%` *CD4* mRNA dropout on human CD4+ 10x libraries in Oh et al. 2020;
  `+26.2%` cell retention, `9,046` → `11,413` cells on murine `GP66:I-Ab` CD4+
  sorts in Khatun et al. 2021). Unsorted samples are unchanged.
- **`FOXP3+` regulatory T-cell auditing and `--exclude-treg`.** `prioritize` and
  `til-prioritize` now compute per-clone raw-count `foxp3_pos_fraction` and
  `is_treg_like` (human `FOXP3` / mouse `Foxp3`, default threshold
  `--treg-foxp3-fraction 0.5`), append `treg_FOXP3` to `risk_flags`, and exclude
  `is_treg_like` clones from candidate selection when `--exclude-treg` is passed
  (removing `61/61` `FOXP3+` Treg picks on Veatch et al. 2022 melanoma CD4+ TILs
  while retaining `100%` of shortlisted neoantigen/self-antigen CD4+ clones).
- **Cell Ranger `hg19` (GENCODE v19) ↔ GRCh38 gene symbol synonyms.** Extended
  `GENE_SYNONYMS` in `tcrsift/genes.py` with 18 official HGNC/GENCODE symbol
  aliases (`ADGRG1`/`GPR56`, `LINC01943`/`RP11-291B21.2`, `PLPP1`/`PPAP2A`,
  `SLF1`/`ANKRD32`, `JAML`/`AMICA1`, `FAAP20`/`C1ORF86`, `FYB1`/`FYB`,
  `BEX3`/`NGFRAP1`, `CYTOR`/`LINC00152`, `ATP5PD`/`ATP5H`, `VSIR`/`C10ORF54`,
  `MCUB`/`CCDC109B`, `MTERF4`/`MTERFD2`, `SARAF`/`TMEM66`,
  `MIR4435-2HG`/`MIR4435-1HG`, `LINC01480`/`AP001055.6`, `RP11-47I22.3`, and
  `CTB-61M7.2`), resolving `40/40` (`100%`) `NeoTCR4` genes on legacy `hg19`
  matrices instead of failing strict gene resolution.
- **Murine Ensembl `-dv` V-gene normalization and static IMGT V-QUEST germline
  downloads.** Normalized Cell Ranger mouse Ensembl shared α/δ V-gene names
  (`TRAV13-4-dv7` → `TRAV13-4/DV7`, `TRAV16d-dv11` → `TRAV16D/DV11`) in
  `normalize_vgene` and `canonicalize_gene` (raising murine germline V+J
  construct assembly to `101/101` = `100%`), switched all eight human and mouse
  IMGT `V-REGION` and `J-REGION` download specs in `datacache.py` from IMGT's
  `GENElect` CGI to static IMGT V-QUEST reference FASTAs while adding session
  cookie + automatic IPv4 fallback to `_default_urlopen` (fixing
  `imgt_trav_leader_mouse` / `imgt_trbv_leader_mouse` redirects as well), and
  pre-filtered `host_species` in `match_clonotypes` before full-table species
  canonicalization.
- **New public benchmark suites.** Added reproducible evaluations and
  provenance documentation for
  [paired blood-tumor, post-ACT, and CD4+ TILs](benchmarks-blood-cd4.md)
  (`benchmarks/public_blood_cd4/`) and
  [murine `--species mouse` models](benchmarks-mouse.md)
  (`benchmarks/public_mouse/`).

## 3.32.0

Priority-tiered round-robin scheduling, single-pass `strict_ab` database matching
in `prioritize`, IMGT `V-REGION` + `J-REGION` germline fallback for 5'-truncated
contigs, first-class `--method stimulation` and `--call-multimers` workflows,
lightweight default installation (`tcrsift[tcrpeg]`), and regenerated `3.32.0`
benchmarks:

- **Priority-tiered round-robin scheduling (`--round-robin-mode tiered`, default).**
  Context presets (`solid-tumor`, `mpe`, `blood-tumor`) now place dedicated
  disease-reactivity signatures (`TumorReactive`, `MANAscore`, `NeoTCR8`,
  `NeoTCR4`, `NeoTCR_PBL`) first in `CONTEXT_SIGNATURES` and schedule all
  qualifying clones from those Tier-1 lists before backfilling any remaining
  `--max-clones` budget from general effector, activation, or proliferation
  panels (Tier 2). Pass `--round-robin-mode flat` to rotate all preset signatures
  with equal weight from round 1.
- **`--database-match` defaults to `strict_ab` in `prioritize` / `til-prioritize`.**
  Because every clone scored by `prioritize` carries paired α/β chains, database
  annotations and viral/MART-1 exclusions now default to exact paired α+β
  matching (`strict_ab`) in a single pass rather than `ab_with_partial`, avoiding
  false-positive β-only collisions with public viral entries on private
  repertoires. Pass `--database-match ab_with_partial` or `b_only` to opt into
  β-only database matching.
- **IMGT `V-REGION` + `J-REGION` germline fallback (`--germline-vdj-fallback`).**
  Added cached human and mouse IMGT `TRAJ`/`TRBJ` `J-REGION` references
  (`tcrsift data download --group germline` / `--group germline_mouse`) and
  `reconstruct_germline_vdj` for 5'-truncated contigs (`partial_VDJ_*`) and bare
  `(v_gene, CDR3, j_gene)` tables (`construct_contig_verified = False`, stamped in
  `{chain}_vdj_source`, `{chain}_vdj_germline_allele`, and `qc_warnings`). Also
  fixed `_build_full_sequences` so chains missing VDJ sequences never emit
  `leader + constant` chimeras.
- **First-class `--method stimulation` and `--call-multimers`.** Added
  `filter_clonotypes_stimulation` (`tcrsift filter --method stimulation` and
  `tcrsift run --method stimulation`, combining per-donor expansion `FDR`/`log2FC`
  with paired within-clone `delta_signature_AIM` / `delta_signature_AcuteActivation`
  shifts over `source: control`) and `call_multimers_for_adata` (`--call-multimers`
  in `prioritize` and `run`).
- **Lightweight core install (`tcrsift[tcrpeg]`) and modularized internals.**
  Moved `torch` and `tcrpeg` to `pip install "tcrsift[tcrpeg]"` while keeping
  bundled NumPy k-mer `Pgen`/`Ppost` models in core, adopted PEP 639 SPDX
  licensing (`setuptools>=77`), and split `assemble.py`, `plots.py`, and `cli.py`
  into focused internal submodules with full backwards-compatible re-exports.

## 3.31.0

VDJdb 2026+ relational-archive compatibility, programmatic pipeline/prioritization
APIs, and regenerated benchmark evaluations with explicit ground-truth provenance:

- **VDJdb 2026-11+ relational release (`chains.tsv.gz` + `records.tsv.gz` and
  `vdjdb.parquet`) support.** Starting with the `2026-11-06` GitHub release,
  `antigenomics/vdjdb-db` ships normalized `chains.tsv.gz`, `records.tsv.gz`, and
  `vdjdb.parquet` inside a versioned `vdjdb-YYYY-MM-DD/` subdirectory rather than
  root-level `vdjdb_full.txt` / `vdjdb.txt` files. Both `tcrsift data download
  --db vdjdb` (`cached_path`) and `load_vdjdb` now resolve relational TSV/gzip
  bundles, `.parquet` tables, and versioned release subdirectories alongside
  legacy flat TSV files, reconstructing paired α/β receptors (`200,180` paired
  out of `280,894` `TRB` entries in `vdjdb-2026-11-06`) via `record_id`.
- **Programmatic `run_pipeline` and `PrioritizeConfig` APIs.** Extracted the
  end-to-end `tcrsift run` workflow into `tcrsift.pipeline.run_pipeline`
  (returning `PipelineResult`), added the typed `PrioritizeConfig` dataclass for
  `run_til_prioritize`, unified CDR3/clone-key/boolean schema normalization in
  `tcrsift._schema` (fixing single-chain `_beta` / `alpha_` key leakage in
  `read_donor_clones`), and auto-derived `TCRsiftConfig` flat-key mappings from
  dataclass definitions.
- **Regenerated TIL and infectious-disease benchmarks on 3.31.0.** Re-ran both
  [`benchmarks/public_til/`](benchmarks.md) and
  [`benchmarks/public_infectious/`](benchmarks-infectious.md) end-to-end with
  machine-readable `ground_truth_provenance` metadata and default `--context
  blood` shortlist/cutoff evaluations.

## 3.30.0

Database-cache, VDJ-loading, and public-deposit compatibility fixes uncovered by
the [infectious-disease and blood benchmark](benchmarks-infectious.md):

- **Cached VDJdb resolves paired α/β entries (`vdjdb_full.txt`).** The VDJdb
  release archive contains both `vdjdb_full.txt` (paired α/β rows) and
  `vdjdb.txt` (one chain per row). `cached_path("vdjdb")` previously returned
  `vdjdb.txt` directly, losing α-chain pairing under `prioritize`'s default
  `--database-match ab_with_partial` (degrading paired VDJdb entries into
  β-only matches during default MART-1 exclusion and `--exclude-known-viral`)
  and giving zero matches under `--database-match strict_ab` when VDJdb was read
  from the `tcrsift data` cache. `cached_path("vdjdb")` now returns the cache
  directory so `vdjdb_full.txt` is preferred, and standalone
  `vdjdb.txt`/`vdjdb.slim.txt` files reconstruct paired α/β entries from
  `complex.id > 0` rows.
- **5'-partial Cell Ranger contigs (`productive == "None"`) with complete V/J
  genes and in-frame CDR3 are retained.** In 5'-truncated VDJ libraries where
  reads start downstream of the V-leader `ATG` (`full_length == False`), Cell
  Ranger marks `productive` as `"None"` rather than `"True"`, while true
  frameshift or stop-codon contigs are marked `"False"`. Retaining `is_cell`,
  `high_confidence` contigs with `productive == "None"`, valid V and J gene
  calls, and a stop-free canonical `C...[FW]` CDR3 recovers 5'-partial cohorts
  (such as GSE171037) while changing ≤0.01% of contigs in standard full-length
  libraries.
- **Compressed, sample-prefixed, and pre-Cell Ranger 6 inputs.** VDJ and GEX
  loaders accept `.csv.gz` contig tables, unambiguous sample-prefixed GEO/10x
  contig and `.h5` filenames, `Sample.gex_path`, and sibling
  `*contig_annotations*.json(.gz)` files to populate full-length `VDJ_aa` and
  `VDJ_nt` when `fwr1..fwr4` CSV columns are absent.

## 3.28.0

Bug fixes affecting receptor assignment and shortlist selection. Re-run
analyses generated with earlier versions.

- **Non-productive contigs are dropped on load.** Cell Ranger's filtered contig
  table retains non-productive contigs, and the loader previously ranked a
  cell's contigs by UMI count alone. A high-UMI non-productive contig could
  therefore displace a cell's productive alpha or beta chain or leave the cell
  without a CDR3 (7% of primary chains in the public lung dataset were
  non-productive). Applies to every command that loads Cell Ranger VDJ output,
  including `run`.
- **Barcodes are matched exactly unless no exact matches exist.** Stripping
  `-1`/`-2` suffixes per barcode could previously pair GEX cell `X-2` with the
  TCR of an unrelated VDJ cell `X-1`. Core-barcode fallback now runs only when
  a sample has zero exact barcode matches, skips ambiguous core barcodes, and
  warns when nothing matches.
- **`.h5` expression matrices are indexed by Ensembl ID**, matching matrix
  directories. Duplicate gene symbols previously caused `.h5`-only samples to
  fail on load, and mixing `.h5` and matrix-directory samples duplicated the
  gene axis.
- **Mid-round budget sharing and `--max-clones 0` fixes.** A list whose sole
  remaining clone had already been selected by another list could still consume
  a turn, violating prefix consistency across `--max-clones` budgets. In
  addition, `--max-clones 0` was treated as a cap equal to the eligible clone
  count, altering `selection_rank` and `selected_signature` relative to a
  non-binding cap.
- **V and J gene filters with repeated barcodes across samples.** When a
  clone's highest-UMI cell shared its barcode string with another cell in a
  different sample, representative gene lookup returned multiple rows and
  caused `--exclude-trav12-2`, `--exclude-v-gene`, and `--include-v-gene` to
  misclassify the clone.
- **Sample-sheet parsing.** Blank `source` cells fall back to the default
  instead of raising an error, numeric YAML sample names are coerced to strings,
  and duplicate sample names are rejected on load.

## 3.27.0

Two shortlist-changing updates motivated by the
[public TIL benchmark](benchmarks.md#4-what-do-the-cell-filters-drop):

- **Default `--min-cd3` lowered from 10 to 3**, matching `run` and `phenotype`.
  Across public TIL datasets, a floor of 10 discarded 65% to 92% of CD8 cells
  without improving ranking among surviving clones. Because more cells and
  clones enter scoring, signature percentiles, background floors, and
  shortlists change. Pass `--min-cd3 10` to reproduce earlier filtering.
- **Fair mid-round budget allocation.** When `--max-clones` expires partway
  through a round, remaining turns are assigned to the signature (and within
  it, the sample) with the fewest picks so far in that round. Previously, partial
  rounds were truncated in preset order, which could starve signatures listed
  near the end of a preset when qualifying lists outnumbered the budget. Full
  rounds retain the existing order.

## 3.26.0

Loader and scoring fixes uncovered while benchmarking on public datasets
([benchmarks](benchmarks.md)). Existing runs that completed without error
produce unchanged scores and shortlists.

- Multi-sample runs reading `filtered_feature_bc_matrix/` directories preserve
  the gene symbols supplied with the matrices instead of failing with missing
  signature genes under `--context solid-tumor`.
- Signature lookup resolves common HGNC symbol renames via a built-in alias
  table (such as `PLA2G16`/`PLAAT3` and `CARS`/`CARS1`), allowing `NeoTCR8`,
  `NeoTCR4`, and `NeoTCR_PBL` to score across older and newer Cell Ranger
  references. Genuinely absent genes still raise an error.
- When only `all_contig_annotations.csv` is present, the VDJ loader filters to
  called cells and high-confidence contigs with valid chain calls.
- `prioritize` logs a warning whenever any single non-lineage cell filter
  removes more than half of the remaining cells.

## 3.25.0

Background cutoff mode defaults to `--signature-evidence either`, allowing
clones to qualify via elevated mean expression or hypergeometric enrichment of
high-scoring cells, each contributing its own ranked list within the shared
200-clone budget. Output tables add high-cell counts, fractions, exact tail
probabilities, and `selected_evidence`. Pass `--signature-evidence mean` for
3.24.0 behavior; `legacy` cutoff mode continues to use mean scores only.

## 3.24.0

Replaced the fixed `score > 0` floor with a size-matched empirical resampling
floor (`--signature-cutoff background`), keeping the top-decile rank gate and
200-clone cap. Pass `--signature-cutoff legacy` for the previous rule. NeoTCR
module scoring now always reads normalized `.X` even when `.raw` is present,
and bare-DataFrame fallback scoring applies `log1p` by default (pass
`log1p=False` if inputs are already log-transformed).

## 3.23.1

Changed `--max-clones` from 100 per patient to **200 across the entire run**
with global `selection_rank`. Required percentile ≥ 0.90 and score > 0,
skipping numerically flat lists and yielding turns from exhausted lists without
padding below cutoffs.

## 3.23.0 (changes from 3.22)

Unified `prioritize` and `til-prioritize` around context presets, CD8 selection
by default, MART-1 exclusion on, viral exclusion off, no minimum mitochondrial
floor, and single-sample support.
