# Future Plans and Known Limitations

This document tracks known limitations and planned future enhancements for TCRsift.

## Species Limitations

**Current Status:** TCRsift currently only supports **human** TCR data.

### Human-Specific Components

The following components are hardcoded for human:

1. **Constant Region Sequences** (`assemble.py`)
   - Uses `pyensembl` with human genome (GRCh38, release 93)
   - Fetches TRAC, TRBC1, TRBC2 sequences from Ensembl
   - Constant region QC checks use human sequence endings

2. **Default Leader Sequences** (`assemble.py`)
   - CD8A, CD28, IgK, TRAC, TRBC leaders are human sequences
   - DNA codons are human codon-optimized

3. **Gene Names** (`loader.py`, `clonotype.py`)
   - Assumes human TCR gene nomenclature (TRAV, TRBV, etc.)
   - CD marker genes (CD3D, CD3E, CD4, CD8A, CD8B) are human

4. **Public Database Matching** (`annotate.py`)
   - VDJdb, IEDB, CEDAR matching assumes human TCR format
   - CDR3 matching uses human conventions

### Multi-Species Support Roadmap

To support other species (mouse, non-human primates, etc.):

#### Phase 1: Configuration-Based Species Selection

```python
# Proposed API
from tcrsift import TCRsiftConfig

config = TCRsiftConfig(species="mouse")  # or "human" (default)
adata = load_samples(sample_sheet, config=config)
```

#### Phase 2: Species-Specific Reference Data

1. **Create species modules:**
   ```
   tcrsift/
     species/
       __init__.py
       human.py      # Current hardcoded data
       mouse.py      # Mouse-specific sequences
       base.py       # Abstract base class
   ```

2. **Each species module provides:**
   - Constant region sequences (TRAC, TRBC1, TRBC2 equivalents)
   - Default leader peptide sequences
   - Gene name patterns (TRAV vs Trav)
   - CD marker gene names
   - Ensembl release and species name

3. **Example species module structure:**
   ```python
   # tcrsift/species/mouse.py

   ENSEMBL_RELEASE = 93
   ENSEMBL_SPECIES = "mus_musculus"

   CONSTANT_REGIONS = {
       "TRAC": "...",
       "TRBC1": "...",
       "TRBC2": "...",
   }

   DEFAULT_LEADERS = {
       "CD8A": {"aa": "...", "dna": "..."},
       # ...
   }

   CD_MARKERS = {
       "CD3": ["Cd3d", "Cd3e", "Cd3g"],
       "CD4": ["Cd4"],
       "CD8": ["Cd8a", "Cd8b1"],
   }

   TCR_GENE_PATTERN = r"Tr[ab][vdjc]\d+"
   ```

#### Phase 3: Database Updates

1. Update VDJdb/IEDB/CEDAR loaders to filter by species
2. Add species column to annotation outputs
3. Support cross-species TCR comparison (humanized mouse TCRs, etc.)

### Contributing Species Support

If you need support for a specific species:

1. Open an issue describing your use case
2. Provide reference sequences if available (constant regions, leaders)
3. Consider contributing a species module (see `tcrsift/species/` structure above)

## Other Planned Features

### Near-Term

- [ ] Batch effect correction for multi-sample analyses
- [ ] Integration with Seurat/Scanpy pipelines
- [ ] HTML report generation improvements

### Medium-Term

- [ ] TCR clustering by sequence similarity
- [ ] Structural modeling integration (AlphaFold, TCRmodel)
- [ ] Neoantigen binding prediction integration

### Long-Term

- [ ] Web interface for non-programmers
- [ ] Cloud deployment support
- [ ] Real-time analysis of sequencing data

## Feedback

Please open issues at https://github.com/pirl-unc/tcrsift/issues for:
- Species support requests
- Feature requests
- Bug reports
