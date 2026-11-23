# Prioritization release notes

Behavior changes in `tcrsift prioritize` and `til-prioritize` that can alter a
shortlist. Record the package version alongside every analysis; the
[prioritization guide](user-guide/til-signatures.md) documents the current
release. Full release history is on
[GitHub](https://github.com/pirl-unc/tcrsift/releases) and
[PyPI](https://pypi.org/project/tcrsift/#history).

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
