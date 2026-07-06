# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Canonical T-cell gene-set signatures.

Single source of truth for the gene sets used across TIL selection
(``til_select.py``) and the per-sample signature scatter
(``plots.py``). Pulled out of ``til_select.py`` so they can be reused
for selecting other kinds of T cells — antigen-response screens,
exhaustion-state phenotyping, healthy-donor panels — without inheriting
the TIL-specific scoring code.

All gene symbols are HGNC (human). When working with Ensembl IDs,
translate via the 10x ``features.tsv.gz`` mapping.

Signatures grouped by intent:

- :data:`EFFECTOR_GENES_HGNC` — cytotoxic-effector panel
  (IFNG, GZMB, PRF1, GNLY, NKG7). Used in cytotoxicity readouts.
  ``ACTIVATION_GENES_HGNC`` is a deprecated alias (#142): this set is
  effector differentiation, *not* immediate-early activation (CD69,
  NR4A1-3, EGR1/2, FOS, JUN) — the old name invited that conflation.
- :data:`NAIVE_STEM_GENES_HGNC` — naïve / stem-memory program
  (TCF7, LEF1, CCR7, SELL, IL7R, CD27, CD28). The "down" pole of the
  effector−naïve differentiation contrast (#141), the best axis for
  separating antigen-expanded clones from naïve-like bystanders.
- :data:`ANTIGEN_RESPONSE_GENES_HGNC` — focal 2-gene recent-antigen
  marker: TNFRSF9 (4-1BB / CD137, the AIM-assay marker per Wölfl
  2007, Frentsch 2005, Bacher 2013) + MKI67 (Ki-67, proliferation).
  Kept as the back-compat subset of :data:`AIM_GENES_HGNC`. Per the
  expansion-correlation analysis (#303), both members positively track
  *culture* clonal expansion (TNFRSF9 ρ≈0.28, MKI67 ρ≈0.46) despite
  TNFRSF9's sparsity — sparse-but-specific, not a reason to demote it.
- :data:`AIM_GENES_HGNC` — activation-induced-marker panel for the
  in-vitro restimulation / antigen-response axis: TNFRSF9 (4-1BB),
  TNFRSF4 (OX40), IL2RA (CD25) + MKI67 (proliferation). These track
  *culture* expansion but sit near zero in TIL (#303), so this is the
  culture-compartment activation readout, not an in-vivo state.
- :data:`CYTOLYTIC_GENES_HGNC` — minimal canonical cytotoxic
  effector pair (PRF1 perforin, GZMB granzyme B). Used as the
  effector readout in Caushi 2021, Krishna 2021, Hanada 2022. Tracks
  *TIL* expansion (ρ≈0.29–0.34) but is flat/negative in culture (#303).
- :data:`EXHAUSTION_GENES_HGNC` — canonical co-inhibitory /
  chronic-antigen exhausted-T-cell set (PDCD1, TOX, LAG3, HAVCR2,
  TIGIT). CTLA4 was dropped (#303): it is an activation/Treg marker,
  not specific to the exhausted state. Every member positively tracks
  *TIL* expansion in the H37 tumor — the expanded intratumoral clones
  are the chronically-stimulated ones. Residency/tumor-reactive
  extension: ENTPD1 (CD39), CXCL13, ITGAE (CD103).
- :data:`TUMOR_REACTIVE_GENES_HGNC` — TIL-resident tumor-reactive
  phenotype: CXCL13 (Workel 2019, Cohen 2022, Veatch 2022) + ENTPD1
  / CD39 (Duhen 2018, Simoni 2018, Thommen 2018).
- :data:`EXPANSION_CORE_GENES_HGNC` — cross-compartment expansion
  core (MKI67, TNFRSF9, EGR2, IFNG, CXCL13, HAVCR2): the markers whose
  clone-level expression tracks expansion in *both* the peptide-culture
  and TIL compartments (#303), so the one panel usable to score
  "expanding & antigen-engaged" regardless of compartment.
- :data:`MARKER_PANEL_HGNC` — the default per-clone GEX scoring panel.
  Widened (#303) to the union of the functional panels above (plus
  CD4/CD8 lineage) so every panel is computable from a single
  extraction pass. Not a single-intent signature (so it is
  intentionally absent from :data:`T_CELL_SIGNATURES`); it is the union
  display panel shown per clone. ``til_select`` re-exports it as
  ``MARKER_GENES_DEFAULT`` for back-compat.

Compartment note (#303): the panels mean different things by
compartment. AIM markers track *in-vitro culture* expansion; cytolytic
and exhaustion track *in-vivo TIL* state/expansion; the expansion core
tracks both. Detectability (% of cells expressing) is *not* the right
yardstick — CD69 is the most detectable activation marker yet
anti-correlates with expansion, while sparse TNFRSF9 works.

:data:`T_CELL_SIGNATURES` is a snake-case name → tuple dict for
convenient iteration.
"""

from __future__ import annotations

EFFECTOR_GENES_HGNC: tuple[str, ...] = ("IFNG", "GZMB", "PRF1", "GNLY", "NKG7")
# Deprecated alias (#142): mis-named — this is cytotoxic effector, not
# immediate-early activation. Kept for back-compat; prefer EFFECTOR_GENES_HGNC.
ACTIVATION_GENES_HGNC: tuple[str, ...] = EFFECTOR_GENES_HGNC
NAIVE_STEM_GENES_HGNC: tuple[str, ...] = (
    "TCF7", "LEF1", "CCR7", "SELL", "IL7R", "CD27", "CD28",
)
ANTIGEN_RESPONSE_GENES_HGNC: tuple[str, ...] = ("TNFRSF9", "MKI67")
# Activation-induced markers (#303): the culture-compartment antigen-response
# axis. ANTIGEN_RESPONSE_GENES_HGNC above is its back-compat 2-gene subset.
AIM_GENES_HGNC: tuple[str, ...] = ("TNFRSF9", "TNFRSF4", "IL2RA", "MKI67")
CYTOLYTIC_GENES_HGNC: tuple[str, ...] = ("PRF1", "GZMB")
# Co-inhibitory / chronic-antigen exhaustion set (#303). CTLA4 dropped:
# it is an activation/Treg marker, not specific to the exhausted state.
EXHAUSTION_GENES_HGNC: tuple[str, ...] = (
    "PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT",
)
TUMOR_REACTIVE_GENES_HGNC: tuple[str, ...] = ("CXCL13", "ENTPD1")
# Cross-compartment expansion core (#303): tracks clonal expansion in BOTH
# the peptide-culture and TIL compartments.
EXPANSION_CORE_GENES_HGNC: tuple[str, ...] = (
    "MKI67", "TNFRSF9", "EGR2", "IFNG", "CXCL13", "HAVCR2",
)

# Default per-clone GEX scoring panel. Widened (#303) to the union of the
# functional panels (plus CD4/CD8 lineage and the immediate-early NR4A/GITR
# context markers) so every panel is computable from one extraction pass. A
# union display panel, not a single intent signature — kept out of
# T_CELL_SIGNATURES below on purpose.
MARKER_PANEL_HGNC: tuple[str, ...] = (
    # lineage
    "CD4", "CD8A", "CD8B",
    # activation / AIM (+ immediate-early context: GITR, NR4A, EGR2)
    "MKI67", "TNFRSF9", "TNFRSF4", "IL2RA", "TNFRSF18", "NR4A1", "NR4A2",
    "EGR2", "CTLA4",
    # cytolytic / effector
    "PRF1", "GZMB", "GZMK", "NKG7", "GNLY", "IFNG",
    # exhaustion / chronic antigen (+ residency / tumor-reactive)
    "PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "ENTPD1", "CXCL13", "ITGAE",
)

T_CELL_SIGNATURES: dict[str, tuple[str, ...]] = {
    "effector":          EFFECTOR_GENES_HGNC,
    "activation":        ACTIVATION_GENES_HGNC,  # deprecated alias of effector
    "naive_stem":        NAIVE_STEM_GENES_HGNC,
    "antigen_response":  ANTIGEN_RESPONSE_GENES_HGNC,
    "aim":               AIM_GENES_HGNC,
    "cytolytic":         CYTOLYTIC_GENES_HGNC,
    "exhaustion":        EXHAUSTION_GENES_HGNC,
    "tumor_reactive":    TUMOR_REACTIVE_GENES_HGNC,
    "expansion_core":    EXPANSION_CORE_GENES_HGNC,
}


# --------------------------------------------------------------------------- #
# Neoantigen-reactivity signatures (#309)
# --------------------------------------------------------------------------- #
# These published "neoantigen-reactivity" signatures are NOT interchangeable
# weighted gene sums — each has a distinct structure, so the registry in
# ``signature_methods`` records each one's genes, per-gene sign, input units,
# and scoring method rather than pretending they share a format.
#
# MANAscore (Zeng/Smith, Nat Commun 2025; PMID 39900903) is the only small
# fixed-gene model, and even it has NO published per-gene coefficients — the
# weights live in a trained RF+linear ensemble. The reproducible part is the
# gene *directions* (+CXCL13, +ENTPD1, -IL7R) and *input units*
# (log-normalized). We carry that as a transparent signed-z proxy:
#   score = (z(CXCL13) + z(ENTPD1) - z(IL7R)) / sqrt(3)   on log1p CP10K.
MANASCORE_UP_HGNC: tuple[str, ...] = ("CXCL13", "ENTPD1")
MANASCORE_DOWN_HGNC: tuple[str, ...] = ("IL7R",)
# Signed per-gene weights for the ``weighted_z`` proxy (all unit magnitude;
# only the sign is published).
MANASCORE_WEIGHTS_HGNC: dict[str, float] = {
    "CXCL13": +1.0,
    "ENTPD1": +1.0,
    "IL7R": -1.0,
}

# NeoTCR8 (Lowery/Rosenberg, Science 2022; PMID 35113651) — 243-gene CD8
# neoantigen-reactive set (Table S10). An UNWEIGHTED gene set scored by rank
# enrichment (scGSEA / score_genes) — there are no published per-gene weights.
NEOTCR8_GENES_HGNC: tuple[str, ...] = (
    "ATP10D", "GZMB", "ENTPD1", "KIR2DL4", "LAYN", "HTRA1", "CD70",
    "CXCR6", "HMOX1", "ADGRG1", "LRRN3", "ACP5", "CTSW", "GALNT2",
    "LINC01480", "CARS", "LAG3", "TOX", "PTPRCAP", "ASB2", "ITGB7",
    "PTMS", "CD8A", "GPR68", "NSMCE1", "ABI3", "SLC1A4", "PLEKHF1",
    "CD8B", "LINC01871", "CCL4", "NKG7", "CLIC3", "NDFIP2", "PLPP1",
    "PCED1B", "CXCL13", "PDCD1", "PRF1", "HLA-DMA", "GPR25", "CD9",
    "TIGIT", "HLA-DRB5", "SYTL3", "SLF1", "NEK1", "CASP1", "SMC4",
    "TSEN54", "PLSCR1", "GNPTAB", "HLA-DPB1", "PLEKHA1", "ARHGAP9",
    "ALOX5AP", "SH3BP1", "NCF4", "NELL2", "GATA3", "PPM1M", "TNFRSF1A",
    "AC022706.1", "MCM5", "HLA-DRB1", "TNFSF10", "TRIM21", "HDLBP",
    "ERN1", "CALHM2", "SASH3", "ACTA2", "MAST4", "CAPG", "MPST",
    "IGFLR1", "GZMA", "CD27", "ITGAE", "SLA2", "RHOC", "COMMD8",
    "MYO1G", "SP140", "PHPT1", "CD2BP2", "PLEKHO1", "STAM", "MRPL16",
    "IL2RB", "ID2", "TESPA1", "GOLGA8B", "MIS18BP1", "VAMP5", "DAPK2",
    "HLA-DPA1", "TSG101", "IL4R", "CCND2", "CTSC", "TRAF3IP3", "NLRC3",
    "ORAI3", "GNLY", "MIR155HG", "CARD16", "CD82", "ECH1", "JAML",
    "EEF1G", "ETFB", "DAXX", "RBM4", "HCST", "RAB27A", "YPEL2",
    "CHST12", "ARPC1B", "PDIA4", "PDIA6", "AC243960.1", "TBC1D10C",
    "PTPN6", "PYCARD", "BST2", "BTN3A2", "MTG1", "MLEC", "DUSP4",
    "GSDMD", "SLAMF1", "IFI6", "PCID2", "GIMAP1", "ITGA1", "CSNK2B",
    "CDK2AP2", "MYO1F", "AC004687.1", "PTTG1", "APOBEC3C", "TSPAN14",
    "MOB3A", "STXBP2", "LCP2", "PLA2G16", "LINC00649", "CST7", "TADA3",
    "SIT1", "APOBEC3G", "SUSD3", "CD3G", "CCL5", "CDC25B", "TNFRSF1B",
    "HMGN3", "THEMIS", "ASF1A", "CTNNB1", "FIBP", "CCDC85B", "POLR3GL",
    "GIMAP6", "ARL6IP1", "CALCOCO2", "CCPG1", "KLRB1", "ACAA2", "ISG15",
    "EIF4A1", "CAT", "MANF", "XAB2", "GRINA", "GLO1", "LSM2", "SLFN5",
    "FKBP1A", "AKNA", "TAP1", "LMO4", "APEH", "C12orf75", "TMEM14A",
    "DNPH1", "C17orf49", "NUDT5", "MGAT1", "CCDC69", "EIF4EBP1", "PDHB",
    "ARL3", "UCP2", "IFI35", "HSBP1", "LYST", "MRFAP1L1", "ITGAL",
    "AIP", "RASAL3", "CAPN1", "ITGB1", "RBPJ", "LBH", "DYNLL1", "NME2",
    "MT1F", "SYNGR2", "ABTB1", "ZGPAT", "CD63", "ILK", "SKA2",
    "TMEM204", "ACO2", "HOPX", "CRIP1", "OXNAD1", "CCS", "GRAP2",
    "GSTO1", "HADHB", "IL16", "PIN4", "CUEDC2", "CALM3", "SAMSN1",
    "HM13", "SNAP23", "LPCAT4", "FAAP20", "EFHD2", "PRDX3", "CCM2",
    "C22orf39", "SDHA", "ARRDC1", "MAP4K1", "NDUFA13", "IL27RA",
    "C14orf119"
)

# NeoTCR4 (Lowery/Rosenberg, Science 2022; PMID 35113651) — 40-gene CD4
# neoantigen-reactive set (Table S10). An UNWEIGHTED gene set scored by rank
# enrichment — no published per-gene weights.
NEOTCR4_GENES_HGNC: tuple[str, ...] = (
    "CXCL13", "HMOX1", "ETV7", "ADGRG1", "PDCD1", "ENTPD1", "CCDC50",
    "TOX", "CD4", "TIGIT", "TNFRSF18", "NMB", "MYL6B", "AHI1", "MAF",
    "IFNG", "LAG3", "CXCR6", "IGFLR1", "DUSP4", "ACP5", "LINC01943",
    "LIMS1", "BATF", "PCED1B", "ITGAL", "YPEL2", "MAL", "PPT1", "ELMO1",
    "MIS18BP1", "TMEM173", "ADI1", "SLA", "GALM", "LBH", "SECISBP2L",
    "CTSB", "C17orf49", "CORO1B"
)

# NeoTCR_PBL (Yossef/Rosenberg, Cancer Cell 2023; PMID 38039963) — 151-gene
# signature of circulating (peripheral-blood) neoantigen-reactive CD8 T cells:
# the genes upregulated in cluster C9 at avg_log2FC >= 0.5 (Table S2D). An
# UNWEIGHTED gene set scored by rank enrichment — no published per-gene weights.
# Verbatim from the supplement (descending log2FC order); a couple of entries are
# Ensembl clone IDs / lncRNAs (AC119396.1, MIR4435-2HG, CYTOR, NEAT1) that only
# score where present in the matrix — kept to preserve fidelity to the published
# set. Distinct from Lowery's NeoTCR8/4 (a separate blood-derived signature).
NEOTCRPBL_GENES_HGNC: tuple[str, ...] = (
    "PASK", "MT2A", "MT1E", "S100A4", "SELL", "LIME1", "DENND10",
    "ALOX5AP", "COTL1", "ITM2A", "MT1X", "EIF3A", "AQP3", "SMC4",
    "UBXN11", "YWHAB", "S100A11", "S1PR4", "MYO1G", "LIMS1", "CHN1",
    "RBPJ", "TMSB10", "CRIP1", "LEF1", "ITGB1", "TMEM123", "CD52",
    "GYPC", "ANXA2", "SELPLG", "NPDC1", "LGALS3", "PAG1", "ANXA5",
    "UCP2", "CORO1B", "CD82", "CD55", "NEAT1", "ITGB7", "VIM",
    "CDC25B", "TLE5", "KLF2", "S100A10", "S100A6", "ACTG1", "TPM4",
    "ITGAE", "MIR4435-2HG", "CNN2", "FYB1", "HLA-DQB1", "FXYD5", "SESN3",
    "HLA-DRB1", "BIN1", "ISG20", "TIGIT", "RHOG", "MED15", "RNASET2",
    "BEX3", "CAPZB", "CORO1A", "LIMD2", "MT1F", "PFN1", "ACTB",
    "SLC1A4", "NAP1L4", "LSP1", "SUSD3", "ICAM2", "TNFRSF25", "CYTOR",
    "GPSM3", "NR3C1", "TTN", "SPOCK2", "GSTK1", "CD5", "PLSCR3",
    "LTB", "OCIAD2", "LEPROTL1", "PRDX1", "NME2", "NDUFA12", "RAB37",
    "SH3KBP1", "NIBAN1", "SPN", "FKBP5", "NDUFB9", "SHMT2", "ZYX",
    "SMAP2", "CSGALNACT1", "TRADD", "ATP5PD", "HLA-DRB5", "ELOVL5", "ARHGDIB",
    "CTLA4", "R3HDM4", "GPR171", "RCSD1", "MYL6B", "TSPO", "IFITM1",
    "CAMK4", "CD27", "TKT", "PPP2R5C", "VSIR", "CD7", "PAICS",
    "S1PR1", "VOPP1", "RASA3", "MPRIP", "MCUB", "EMB", "GPI",
    "RASGRP2", "MTERF4", "DUSP16", "SIT1", "CARHSP1", "ZMIZ1", "GMFG",
    "RAC2", "PLEC", "DGKA", "SLC25A5", "PFKL", "OPTN", "ETS1",
    "MSC", "ARID5B", "HLA-DPA1", "MFNG", "P2RY8", "SARAF", "MAP3K1",
    "TRAPPC5", "AC119396.1", "SLC16A3", "CD3D"
)


# --------------------------------------------------------------------------- #
# Cell-type / T-state registries for the per-cell annotator (#312)
# --------------------------------------------------------------------------- #
# The generalized cell-typing gene sets, shared by the per-cell annotator
# (:mod:`tcrsift.annotate_cells`). ``T_CELL_SIGNATURES`` above stays the
# FUNCTIONAL-program registry (effector/exhaustion/…); these are cell-LINEAGE
# and T-STATE registries. Tumor / cancer-testis-antigen typing is deliberately
# OUT of scope (belongs in oncoref, #310) and is not included — tissue
# restriction is done with ``annotate_clusters(allowed_types=...)``.
#
# TCR/Ig CONSTANT genes are lineage evidence (kept); clonal V/J segments are
# excluded elsewhere so structure/labels aren't driven by clonotype.
CELL_TYPE_SIGNATURES: dict[str, list[str]] = {
    # --- lymphoid ---
    "T cell": ["CD3D", "CD3E", "CD3G", "CD247", "CD2", "CD5", "CD7",
               "TRAC", "TRBC1", "TRBC2"],
    "NK cell": ["NKG7", "KLRD1", "KLRF1", "NCR1", "NCR3", "GNLY", "NCAM1",
                "KLRC1", "KLRK1", "CD160", "CD244", "FCGR3A", "EOMES",
                "KIR2DL3", "KIR3DL1", "TYROBP", "PRF1"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "Plasma cell": ["MZB1", "XBP1", "SDC1", "PRDM1", "DERL3", "TNFRSF17", "JCHAIN"],
    # --- myeloid ---
    "Macrophage": ["CD68", "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "APOE",
                   "APOC1", "TREM2", "MARCO", "SELENOP"],
    "Monocyte": ["FCN1", "VCAN", "S100A8", "S100A9", "S100A12", "CD14", "LYZ"],
    "Dendritic cell": ["CLEC9A", "XCR1", "CD1C", "FCER1A", "LILRA4", "CLEC10A", "LAMP3"],
    # KIT (CD117) is NOT mast-specific (HSPC/ILC2 express it): require the mast
    # granule genes — tryptases, carboxypeptidase CPA3, FcERI-beta, HDC.
    "Mast cell": ["TPSAB1", "TPSB2", "CPA3", "MS4A2", "HDC"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2", "G0S2", "FUT4", "CEACAM3", "IFITM2"],
    # --- stroma / vasculature ---
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA",
                   "PDGFRB", "COL6A1", "COL6A2", "COL6A3", "FBLN1", "MFAP5",
                   "SFRP2", "GSN"],
    "Smooth muscle": ["ACTA2", "MYH11", "TAGLN", "CNN1", "MYL9", "DES",
                      "PLN", "LMOD1", "MYLK"],
    "Pericyte": ["RGS5", "NOTCH3", "PDGFRB", "MCAM", "KCNJ8", "HIGD1B",
                 "CSPG4", "ABCC9"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5", "RAMP2", "EGFL7",
                    "PLVAP", "ERG", "FLI1"],
    "Platelet/megakaryocyte": ["PPBP", "PF4", "ITGA2B", "GP9", "GP1BB",
                               "TUBB1", "NRGN", "GNG11"],
    "Erythroid": ["GYPA", "ALAS2", "CA1", "SLC4A1", "KLF1", "GATA1"],
}

# Cell types that can appear in a Ficoll-isolated PBMC / peptide-stimulated
# culture — pass to ``annotate_clusters(allowed_types=...)`` so a cultured
# moDC/macrophage cluster can't win a stromal/granulocyte signature it shares
# an activation program with. None = consider every type (whole-tissue TIL).
PBMC_CULTURE_TYPES: frozenset = frozenset({
    "T cell", "NK cell", "B cell", "Plasma cell",
    "Macrophage", "Monocyte", "Dendritic cell",
})

# T-cell sub-states, one layer below the lineage call. A cluster's state is the
# argmax of these, subject to the gates in DEFAULT_GATES (a shared master TF is
# not enough — the defining effector cytokine must clear a floor).
T_STATE_SIGNATURES: dict[str, list[str]] = {
    "naive/Tcm": ["CCR7", "SELL", "TCF7", "LEF1", "IL7R"],
    "effector/cytotoxic": ["GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1",
                           "GNLY", "NKG7", "FGFBP2", "KLRG1"],
    "Trm": ["ITGAE", "ZNF683", "CXCR6", "ITGA1", "CD69"],
    "Tex/CXCL13+": ["CXCL13", "PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT",
                    "ENTPD1", "LAYN"],
    "Tfh": ["CXCR5", "BCL6", "ICOS", "TOX2", "CD200"],
    "Th1": ["TBX21", "IFNG", "CXCR3", "IL12RB2"],
    "Th2": ["GATA3", "IL13", "IL4", "IL5", "IL17RB"],
    "Th17": ["RORC", "IL17A", "IL17F", "CCR6", "IL23R", "IL22"],
    "proliferating": ["MKI67", "TOP2A", "STMN1", "TYMS"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "CCR8"],
    "IFN-stimulated": ["ISG15", "IFIT1", "IFIT3", "MX1", "OAS1"],
}

# B-cell sub-states (applied within B-cell clusters).
B_STATE_SIGNATURES: dict[str, list[str]] = {
    "naive B": ["TCL1A", "IGHD", "FCER2", "IL4R"],
    "memory B": ["CD27", "TNFRSF13B", "CD80", "AIM2"],
    "germinal-center B": ["BCL6", "RGS13", "AICDA", "MEF2B", "LRMP"],
}

# CD8/CD4 for T-lineage; CD3 = pan-T lineage genes that are TRULY T-specific
# (absent in NK) to disambiguate cytotoxic CD8 from NK; NKspec = NK receptors
# cytotoxic CD8 LACKS (the positive NK signal).
LINEAGE_GENES: dict[str, list[str]] = {
    "CD8": ["CD8A", "CD8B"],
    "CD4": ["CD4"],
    "CD3": ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "CD5"],
    "NKspec": ["GNLY", "KLRF1", "NCR1", "NCAM1", "FCGR3A", "KLRC1", "CD160"],
}

# Rare-but-important subtypes that keep their distinguishing-gene suffix even
# when a small fraction of their lineage (exempt from the subtype-merge rule).
PROTECTED_SUBTYPES: frozenset = frozenset({"pDC", "cDC1"})
