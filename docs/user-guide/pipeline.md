# Pipeline Overview

TCRsift processes single-cell TCR data through a series of steps, each building on the previous.

## Pipeline Steps

```
┌─────────────────┐
│  Sample Sheet   │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Load Data     │  ← CellRanger VDJ + GEX
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Phenotype     │  ← CD4/CD8 classification
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Clonotype     │  ← Aggregate by CDR3
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│    Filter       │  ← Tiered selection
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Annotate      │  ← VDJdb, IEDB, CEDAR
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Match TIL     │  ← Optional TIL matching
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Assemble      │  ← Full-length sequences
└─────────────────┘
```

## 1. Load Data

Parses CellRanger VDJ and GEX outputs:

- **VDJ**: `filtered_contig_annotations.csv` for TCR sequences
- **GEX**: `filtered_feature_bc_matrix/` for gene expression

```bash
tcrsift load --sample-sheet samples.yaml -o loaded.h5ad
```

The output is an AnnData object with:

- Cell barcodes as observations
- Gene expression in `X`
- VDJ annotations in `obs` columns

## 2. Phenotype Cells

Classifies cells as CD4+ or CD8+ based on gene expression:

```bash
tcrsift phenotype -i loaded.h5ad -o phenotyped.h5ad --cd4-cd8-ratio 3.0
```

Classification categories:

| Category | Criteria |
|----------|----------|
| Confident CD8+ | CD8/CD4 ratio > threshold |
| Confident CD4+ | CD4/CD8 ratio > threshold |
| Likely CD8+ | CD8 > 0 and CD4 = 0 |
| Likely CD4+ | CD4 > 0 and CD8 = 0 |
| Unknown | Similar expression levels |

## 3. Aggregate Clonotypes

Groups cells by CDR3 sequences:

```bash
tcrsift clonotype -i phenotyped.h5ad -o clonotypes.csv --group-by CDR3ab
```

Grouping options:

- `CDR3ab`: Match by both alpha and beta chains (stricter)
- `CDR3b_only`: Match by beta chain only (more permissive)

## 4. Filter Clonotypes

Applies tiered filtering to prioritize antigen-specific clones:

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method threshold --tcell-type cd8
```

See [Filtering Strategies](filtering.md) for detailed options.

## 5. Annotate Clonotypes

Matches against public TCR databases:

```bash
tcrsift annotate -i filtered/tier1.csv -o annotated.csv \
    --vdjdb /path/to/vdjdb \
    --iedb /path/to/iedb
```

Annotations include:

- Known epitope specificity
- Viral vs tumor antigens
- Database source

## 6. Match TIL (Optional)

For tumor studies, match culture clonotypes against TIL:

```bash
tcrsift match-til -i annotated.csv --til-csv til_clonotypes.csv -o matched.csv

# Multiple TIL samples without a sample sheet:
tcrsift match-til -i annotated.csv -o matched.csv \
  --til-sample T1=csv:/path/to/til_t1.csv \
  --til-sample T2=h5ad:/path/to/til_t2.h5ad
```

This identifies clones that:

- Were expanded in culture AND present in tumor
- Are TIL-specific (not in culture)

TIL samples are excluded from culture aggregation/filtering and only used for matching.

For TIL-only analysis (no culture input), use:

```bash
tcrsift til-clonotype -o til_clonotypes.csv \
  --til-sample T1=csv:/path/to/til_t1.csv \
  --til-sample T2=h5ad:/path/to/til_t2.h5ad
```

For TIL-only 10x VDJ+GEX timepoint prioritization (CD8 + enrichment + immunogenic masks), use:

```bash
tcrsift til-select \
  --data-dir /path/to/til_timepoint_data \
  --rank-by marker_score_z_mean
```

`til-select` runs in v2-compatible CSV mode by default, so with the same data/options
it reproduces legacy `harmonize_abtcr_timepoints.py` CSV outputs. Figures and PDFs may
still differ byte-for-byte.

Expected per-timepoint files in `--data-dir`:
- `consensus_annotations.<TP>.csv`
- `clonotypes.<TP>.csv`
- `filtered_contig_annotations.<TP>.csv`
- `sample_filtered_feature_bc_matrix.<TP>.h5`

## 7. Assemble Full Sequences

Build full-length TCR sequences:

```bash
tcrsift assemble -i annotated.csv -o sequences.csv \
    --include-constant \
    --linker T2A \
    --fasta sequences.fasta
```

Output includes:

- Leader peptide (from contigs)
- Variable region (VDJ)
- Constant region (canonical TRAC / TRBC1 / TRBC2 — see below)
- Single-chain construct (beta-T2A-alpha)

### How the constant region NT is built

The C region's amino-acid sequence comes from a packaged canonical FASTA
(`tcrsift/refseqs/canonical_constants.fasta`, sourced from pyensembl GRCh38
release 110 and cross-checked against UniProt P01848 / P01850 / A0A5B9 —
see issue #100 for the provenance story).

The C region's **nucleotide** sequence is built as a hybrid of donor-real
contig bytes and codon-optimized canonical:

1. **J→C junction codon and surrounding C-region codons come from the
   CellRanger contig where possible.** The contig retains the donor's
   actual NT at the boundary (which the #91 fix trimmed off `vdj_*_nt`
   at the clonotype-aggregation step). When contigs are available, the
   assembler locates `vdj_{chain}_nt` in the contig and copies the bytes
   immediately past it for as many codons as agree with the canonical
   AA. This keeps the assembled NT faithful to what's actually in the
   donor's cells at the J→C boundary.
2. **Incomplete codons at the contig 3' edge are completed using the
   canonical reference.** If the contig provides 1-2 nt of an otherwise-
   incomplete codon, the assembler picks a codon that:
   - starts with the contig's partial bytes (preserving donor fidelity), AND
   - codes for the canonical AA at that position.

   When the codon-optimized canonical codon matches both constraints,
   it's used; otherwise any compatible codon is chosen. If no codon can
   satisfy both (donor partial bytes are incompatible with the canonical
   residue), the partial bytes are discarded and the canonical codon
   takes that position.
3. **Everything past the contig's coverage uses codon-optimized canonical
   NT.** The deep C region is invariant across donors and benefits from
   codon optimization for downstream synthesis.

The breakdown is recorded in `{chain}_constant_source`, e.g.:

```
canonical:TRBC1 (contig-verified) [contig(8 codons, partial completed) + canonical-codon-opt]
```

If the contig translation disagrees with the canonical AA at some position,
the assembler switches to canonical-codon-optimized at that point and adds
a QC warning naming the AA position of the divergence. The protein is
always canonical; only the NT can carry donor-specific bytes at the boundary.

To opt out and force codon-optimized canonical for the whole C region —
e.g., when reproducing pre-1.3 assemblies — omit `--contigs-dir` so no
contigs are loaded. `_blend_constant_nt_with_contig` then sees no contig
NT and returns the canonical-codon-optimized sequence unchanged.

## Running the Complete Pipeline

Use `tcrsift run` to execute all steps:

```bash
tcrsift run \
    --sample-sheet samples.yaml \
    --output-dir results/ \
    --vdjdb /path/to/vdjdb \
    --tcell-type cd8 \
    --method threshold \
    # report generation is enabled by default; use --no-report to disable
```

This creates a complete output directory with all intermediate files and a summary report.
