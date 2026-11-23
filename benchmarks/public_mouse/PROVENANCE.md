# Ground-Truth Data Provenance (`benchmarks/public_mouse`)

This document records the primary open-access data deposits, sample-to-accession mappings, experimental assays establishing biological ground truth, preprocessing and clonal aggregation rules in `run.py`, reference database versions, and circularity controls for all three murine (`--species mouse`) cohorts evaluated in `docs/benchmarks-mouse.md`.

---

## Summary Table

| Cohort ID | Study | Biological Context | Expression & VDJ Deposit | Ground-Truth Label Source | Assay Defining Positives / Contrasts | Circularity & Controls |
| --- | --- | --- | --- | --- | --- | --- |
| `mc38` | Multi-tissue murine MC38 colon adenocarcinoma ([GSE168944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168944)) | Syngeneic C57BL/6 MC38 subcutaneous tumor, peripheral blood, draining lymph node (`dLN`), and spleen (`Mouse01` & `Mouse02`) | GEO `GSE168944`: `GSM5172668`–`GSM5172675` (GEX `barcodes.tsv.gz`, `genes.tsv.gz`, `matrix.mtx.gz`) paired with `GSM5172684`–`GSM5172691` (VDJ `filtered_contig_annotations.csv.gz`) | Exact paired `CDR3ab` clonal tracking across 4 anatomical compartments within the same mouse (`Tumor`, `Blood`, `dLN`, `Spleen`) | Multi-organ paired scRNA-seq + scTCR-seq comparing tumor-exclusive expanded clones (`>=5` tumor cells and `0` lymphoid cells) against lymphoid-shared non-enriched clones (`>=2` tumor cells, `>=1` lymphoid cell, and `<=2x` tumor vs lymphoid frequency), plus blood-to-tumor clonal overlap | Independent of all TCRsift expression signatures and absent from VDJdb; also serves as a C57BL/6 out-of-strain specificity control against NOD `IGRP` VDJdb lookups |
| `kasmani` | Kasmani, Ciecko, Brown, Petrova, Gorski, Chen & Cui, *Life Science Alliance* 2022 ([DOI:10.26508/lsa.202201503](https://doi.org/10.26508/lsa.202201503); PMID:35667687) | Autoimmune diabetes NOD/ShiLtJ mouse pancreatic islets and spleen (`Round1` & `Round2`) | GEO [GSE200608](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200608): `GSM6038943`–`GSM6038948` (GEX) paired with `GSM6038949`–`GSM6038954` (VDJ `filtered_contig_annotations.csv.gz`) | Physical FACS gating with `H-2Kd IGRP206-214 / NRP-V7` (`KYNKANVFL`) MHC-I tetramers across pancreatic islets (`Islet_TetPos` vs `Islet_TetNeg`) and spleen (`Spleen_TetPos`) | (A) Same cognate `IGRP` (`KYNKANVFL`) specificity across target organ (`Islet_TetPos`) vs circulating lymphoid reservoir (`Spleen_TetPos`), including identical sister `CDR3ab` clones; (B) within-islet `TetPos` vs `TetNeg` CD8+ T cells | **Ingested into VDJdb (`PMID:35667687`, 780 `MusMusculus` records).** Evaluated both on unmodified VDJdb (**circular**) and with all `PMID:35667687` records held out (**independent**). Independent of all TCRsift expression signatures |
| `khatun` | Khatun, Kasmani, Zander et al., *J Exp Med* 2021 ([DOI:10.1084/jem.20200650](https://doi.org/10.1084/jem.20200650); PMID:33201171) | Day-10 acute LCMV Armstrong infection splenic CD4+ T cells across 5 C57BL/6 mice (`M1`–`M5`) | GEO [GSE158896](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158896): `GSM4815363`–`GSM4815367` (GEX) paired with `GSM4815368`–`GSM4815372` (VDJ `filtered_contig_annotations.csv.gz`) | Physical FACS sorting of `CD4+ CD44hi` `I-Ab GP66-77` (`DIYKGVYQFKSV`) MHC class II tetramer-binding T cells across 5 genetically identical mice | Virus-specific murine CD4+ T-cell phenotyping (`--species mouse --tcell-type cd4`), Th1 vs Tfh/AIM clonal state selection, cross-mouse public receptor sharing, VDJdb `LCMV` lookup, and murine CD4 full-length construct assembly | Independent of all TCRsift signatures and absent from VDJdb (`0` records with `PMID:33201171`, verified in `summary.json`) |

---

## Cohort-by-Cohort Construction & Format Notes

### 1. Multi-Compartment MC38 Colon Adenocarcinoma (`mc38` — `GSE168944`)
- **Library pairing (`Mouse01` and `Mouse02`)**:
  - `Mouse01_Tumor`: `GSM5172668` (`DJC01` GEX) + `GSM5172684` (`DJT01` VDJ)
  - `Mouse02_Tumor`: `GSM5172669` (`DJC02` GEX) + `GSM5172685` (`DJT02` VDJ)
  - `Mouse01_Blood`: `GSM5172670` (`DJC04` GEX) + `GSM5172686` (`DJT04` VDJ)
  - `Mouse02_Blood`: `GSM5172671` (`DJC05` GEX) + `GSM5172687` (`DJT05` VDJ)
  - `Mouse01_Spleen`: `GSM5172672` (`DJC09` GEX) + `GSM5172688` (`DJT09` VDJ)
  - `Mouse02_Spleen`: `GSM5172673` (`DJC10` GEX) + `GSM5172689` (`DJT10` VDJ)
  - `Mouse01_dLN`: `GSM5172674` (`DJC12` GEX) + `GSM5172690` (`DJT12` VDJ)
  - `Mouse02_dLN`: `GSM5172675` (`DJC13` GEX) + `GSM5172691` (`DJT13` VDJ)
- **Cell Ranger v2 `genes.tsv.gz` normalization**: `download.py` converts 2-column `genes.tsv.gz` (`ENSMUSG... \t Symbol`) into standard 3-column `features.tsv.gz` (`ENSMUSG... \t Symbol \t Gene Expression`) without altering counts or barcodes.
- **Species-adapted MGI `TumorReactive` ablation**: Human CD8+ tumor-reactive TILs transcribe high levels of *CXCL13*, whereas murine CD8+ TILs rarely express *Cxcl13* and instead upregulate *Ccl3*, *Ccl4*, *Entpd1*, *Tox*, and *Itgae* (`MOUSE_SIGNATURES["TumorReactive"]`). `run.py` measures single-cell *Cxcl13* vs *Ccl3*/*Ccl4* detection rates on MC38 CD8+ TILs and compares MGI `TumorReactive` directly against naive human `CXCL13`-containing `TumorReactive`.

### 2. Kasmani et al. 2022 NOD Islet & Spleen `IGRP` Tetramer CD8+ Cohort (`kasmani` — `GSE200608`)
- **Library pairing (`Round1` and `Round2`)**:
  - `R1_Islet_TetNeg`: `GSM6038943` + `GSM6038949`
  - `R1_Islet_TetPos`: `GSM6038944` + `GSM6038950`
  - `R1_Spleen_TetPos`: `GSM6038945` + `GSM6038951`
  - `R2_Islet_TetNeg`: `GSM6038946` + `GSM6038952`
  - `R2_Islet_TetPos`: `GSM6038947` + `GSM6038953`
  - `R2_Spleen_TetPos`: `GSM6038948` + `GSM6038954`
- **Between-compartment baseline normalization**: Because `tcrsift prioritize` ranks clones *within* each `(donor, sample, lineage)` stratum (giving every sample an internal mean percentile of `0.50`), comparing clones sorted from two separate physical compartments (`Islet_TetPos` vs `Spleen_TetPos`, or `Islet_TetPos` vs `Islet_TetNeg`) requires scoring cells against a shared per-round baseline (`groups=adata_cd8.obs["patient_id"]`) before ranking clones within each round. Identical sister `CDR3ab` clonotypes captured in both `Spleen_TetPos` and `Islet_TetPos` of the same round provide an internal receptor-matched control for target-organ antigen engagement.
- **Circularity control in VDJdb**: `PMID:35667687` contributes `780` `MusMusculus` `TRB` rows in `vdjdb-2026-11-06`. `run.py` reports both unmodified VDJdb (**circular** self-retrieval) and study-held-out VDJdb (all `35667687` rows removed prior to indexing), alongside out-of-strain C57BL/6 MC38 clones to measure false-positive `IGRP` collision rates.

### 3. Khatun et al. 2021 LCMV `GP66:I-Ab` Murine CD4+ Cohort (`khatun` — `GSE158896`)
- **Library pairing (`M1`–`M5`)**: `GSM4815363`–`GSM4815367` (`M1_GEX`–`M5_GEX`) paired with `GSM4815368`–`GSM4815372` (`M1_TCR`–`M5_TCR`).
- **5'-partial VDJ contigs & germline construct assembly**: All three GEO deposits (`GSE168944`, `GSE200608`, `GSE158896`) distribute pre-Cell Ranger 6 18-column `filtered_contig_annotations.csv` files (`v_gene`, `j_gene`, `c_gene`, `cdr3`, `cdr3_nt`, without `fwr1_nt`..`fwr4_nt`). Full-length murine `Trac` and `Trbc1`/`Trbc2` single-chain constructs are assembled via `assemble_full_sequences(..., species="mouse", germline_vdj_fallback=True)` against cached IMGT `Mus_musculus` `V-REGION` and `J-REGION` references (`tcrsift data download --group germline_mouse`) and validated with `validate_sequences(..., species="mouse")`.

### 4. Reference Databases
- **VDJdb**: `vdjdb-2026-11-06` (`14,008` `MusMusculus` `TRB` rows; `9,647` paired α/β rows), filtered via `host_species_filter="same"` for `species="mouse"`.
- **IMGT Mouse Germline References**: Static IMGT V-QUEST reference directory FASTAs (`https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/{TRAV,TRBV,TRAJ,TRBJ}.fasta`), cached via `tcrsift data download --group germline_mouse`.
