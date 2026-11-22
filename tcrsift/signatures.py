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

"""Canonical T-cell and single-cell lineage gene-set signatures (human & mouse).

Single source of truth for the gene sets used across clone prioritization,
TIL selection, per-sample scatter plots, and single-cell cluster annotation.

Gene symbols are stored uppercase to match :func:`tcrsift.genes.adata_symbol_array`
uniformly across human HGNC (``GZMB``, ``CD8B``) and mouse MGI (``Gzmb`` → ``GZMB``,
``Cd8b1`` → ``CD8B1``). Species-specific registries handle genes that lack 1:1
orthologs across species (for example, human ``GNLY``, ``GZMH``, ``FGFBP2``,
``KLRF1``, ``NCR3``, ``FCN1``, ``CD1C``, ``LILRA4``, ``FCGR3B`` have no mouse
ortholog, and mouse NK cells express ``Ncr1``/``Klrb1c``/``Itga2``/``Klra*``
rather than ``NCAM1``/``KIR*``).
"""

from __future__ import annotations

from .species import normalize_species

EFFECTOR_GENES_HGNC: tuple[str, ...] = ("IFNG", "GZMB", "PRF1", "GNLY", "NKG7")
# Deprecated alias (#142): cytotoxic effector, not immediate-early activation.
ACTIVATION_GENES_HGNC: tuple[str, ...] = EFFECTOR_GENES_HGNC
NAIVE_STEM_GENES_HGNC: tuple[str, ...] = (
    "TCF7", "LEF1", "CCR7", "SELL", "IL7R", "CD27", "CD28",
)
ANTIGEN_RESPONSE_GENES_HGNC: tuple[str, ...] = ("TNFRSF9", "MKI67")
AIM_GENES_HGNC: tuple[str, ...] = ("TNFRSF9", "TNFRSF4", "IL2RA", "MKI67")
CYTOLYTIC_GENES_HGNC: tuple[str, ...] = ("PRF1", "GZMB")
EXHAUSTION_GENES_HGNC: tuple[str, ...] = (
    "PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT",
)
TUMOR_REACTIVE_GENES_HGNC: tuple[str, ...] = ("CXCL13", "ENTPD1")
EXPANSION_CORE_GENES_HGNC: tuple[str, ...] = (
    "MKI67", "TNFRSF9", "EGR2", "IFNG", "CXCL13", "HAVCR2",
)

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
    "activation":        ACTIVATION_GENES_HGNC,
    "naive_stem":        NAIVE_STEM_GENES_HGNC,
    "antigen_response":  ANTIGEN_RESPONSE_GENES_HGNC,
    "aim":               AIM_GENES_HGNC,
    "cytolytic":         CYTOLYTIC_GENES_HGNC,
    "exhaustion":        EXHAUSTION_GENES_HGNC,
    "tumor_reactive":    TUMOR_REACTIVE_GENES_HGNC,
    "expansion_core":    EXPANSION_CORE_GENES_HGNC,
}

# Mouse MGI equivalents (uppercase for adata_symbol_array matching).
# Granulysin (GNLY) has no mouse ortholog; mouse tumor-reactive CD8 TILs in
# syngeneic models rarely express Cxcl13 and instead co-express Entpd1 (CD39),
# Pdcd1 (PD-1), Havcr2 (TIM-3), and Itgae (CD103).
EFFECTOR_GENES_MGI: tuple[str, ...] = ("IFNG", "GZMB", "PRF1", "NKG7")
TUMOR_REACTIVE_GENES_MGI: tuple[str, ...] = ("ENTPD1", "PDCD1", "HAVCR2", "ITGAE")
EXPANSION_CORE_GENES_MGI: tuple[str, ...] = (
    "MKI67", "TNFRSF9", "EGR2", "IFNG", "ENTPD1", "HAVCR2",
)
MARKER_PANEL_MGI: tuple[str, ...] = (
    "CD4", "CD8A", "CD8B1",
    "MKI67", "TNFRSF9", "TNFRSF4", "IL2RA", "TNFRSF18", "NR4A1", "NR4A2",
    "EGR2", "CTLA4",
    "PRF1", "GZMB", "GZMA", "GZMK", "GZMC", "NKG7", "IFNG",
    "PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "ENTPD1", "ITGAE", "CCL3", "CCL4",
)

MOUSE_T_CELL_SIGNATURES: dict[str, tuple[str, ...]] = {
    "effector":          EFFECTOR_GENES_MGI,
    "activation":        EFFECTOR_GENES_MGI,
    "naive_stem":        NAIVE_STEM_GENES_HGNC,
    "antigen_response":  ANTIGEN_RESPONSE_GENES_HGNC,
    "aim":               AIM_GENES_HGNC,
    "cytolytic":         CYTOLYTIC_GENES_HGNC,
    "exhaustion":        EXHAUSTION_GENES_HGNC,
    "tumor_reactive":    TUMOR_REACTIVE_GENES_MGI,
    "expansion_core":    EXPANSION_CORE_GENES_MGI,
}


# --------------------------------------------------------------------------- #
# Neoantigen-reactivity signatures (#309)
# --------------------------------------------------------------------------- #
MANASCORE_UP_HGNC: tuple[str, ...] = ("CXCL13", "ENTPD1")
MANASCORE_DOWN_HGNC: tuple[str, ...] = ("IL7R",)
MANASCORE_WEIGHTS_HGNC: dict[str, float] = {
    "CXCL13": +1.0,
    "ENTPD1": +1.0,
    "IL7R": -1.0,
}

# NeoTCR8 (Lowery/Rosenberg, Science 2022; PMID 35113651) — 243-gene human CD8 set.
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
    "C14orf119",
)

# NeoTCR4 (Lowery/Rosenberg, Science 2022; PMID 35113651) — 40-gene human CD4 set.
NEOTCR4_GENES_HGNC: tuple[str, ...] = (
    "CXCL13", "HMOX1", "ETV7", "ADGRG1", "PDCD1", "ENTPD1", "CCDC50",
    "TOX", "CD4", "TIGIT", "TNFRSF18", "NMB", "MYL6B", "AHI1", "MAF",
    "IFNG", "LAG3", "CXCR6", "IGFLR1", "DUSP4", "ACP5", "LINC01943",
    "LIMS1", "BATF", "PCED1B", "ITGAL", "YPEL2", "MAL", "PPT1", "ELMO1",
    "MIS18BP1", "TMEM173", "ADI1", "SLA", "GALM", "LBH", "SECISBP2L",
    "CTSB", "C17orf49", "CORO1B",
)

# NeoTCR_PBL (Yossef/Rosenberg, Cancer Cell 2023; PMID 38039963) — 151-gene human CD8 PBL set.
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
    "TRAPPC5", "AC119396.1", "SLC16A3", "CD3D",
)


# --------------------------------------------------------------------------- #
# Cell-type / T-state / B-state registries for per-cell annotation (#312)
# --------------------------------------------------------------------------- #

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

PBMC_CELL_TYPE_SIGNATURES: dict[str, list[str]] = CELL_TYPE_SIGNATURES

_SOLID_TUMOR_PARENCHYMAL_SIGNATURES: dict[str, list[str]] = {
    "Mesothelial": ["MSLN", "UPK3B", "WT1", "CALB2", "PDPN", "LRRN4", "KRT19"],
    "Osteoclast": ["ACP5", "CTSK", "MMP9", "OSCAR", "ATP6V0D2", "DCSTAMP", "CALCR"],
    "Skeletal muscle": ["MYH1", "MYH2", "ACTN2", "TNNT3", "TNNC2", "CKM", "TTN",
                        "MYLPF"],
    "Adipocyte": ["ADIPOQ", "LEP", "PLIN1", "FABP4", "CFD", "GPD1", "LIPE"],
    "Schwann/nerve": ["PLP1", "MPZ", "S100B", "SOX10", "NRXN1", "CDH19", "GPM6B"],
}

SOLID_TUMOR_CELL_TYPE_SIGNATURES: dict[str, list[str]] = {
    **CELL_TYPE_SIGNATURES,
    **_SOLID_TUMOR_PARENCHYMAL_SIGNATURES,
}

# Mouse (Mus musculus) cell-type registry (MGI symbols uppercased).
# Replaces human-only receptors (GNLY, KLRF1, NCR3, KIR*, FCN1, S100A12, CD1C,
# LILRA4, FCGR3B) and avoids NCAM1 on NK cells (neural/stromal in mouse).
MOUSE_CELL_TYPE_SIGNATURES: dict[str, list[str]] = {
    **CELL_TYPE_SIGNATURES,
    "NK cell": ["NCR1", "KLRB1C", "KLRB1B", "ITGA2", "KLRA1", "KLRA4", "KLRA7",
                "KLRA8", "KLRD1", "KLRK1", "KLRC1", "NKG7", "PRF1", "EOMES",
                "TYROBP", "FCGR3", "CD244A", "CD160"],
    "Monocyte": ["LY6C2", "CCR2", "LYZ2", "VCAN", "S100A8", "S100A9", "CD14",
                 "CHIL3", "PLAC8", "FCGR4", "NR4A1", "MS4A6C"],
    "Dendritic cell": ["ITGAX", "ZBTB46", "FLT3", "CLEC9A", "XCR1", "SIRPA",
                       "CD209A", "MGL2", "CLEC10A", "SIGLECH", "BST2", "CCR9",
                       "CCR7", "FSCN1", "LAMP3"],
    "Mast cell": ["TPSAB1", "TPSB2", "MCPT4", "CMA1", "CPA3", "MS4A2", "HDC"],
    "Neutrophil": ["LY6G", "S100A8", "S100A9", "CSF3R", "CXCR2", "RETNLG",
                   "G0S2", "MMP9", "IL1B", "ACOD1", "IFITM2"],
    "Erythroid": ["GYPA", "ALAS2", "CAR1", "SLC4A1", "KLF1", "GATA1"],
}

MOUSE_SOLID_TUMOR_CELL_TYPE_SIGNATURES: dict[str, list[str]] = {
    **MOUSE_CELL_TYPE_SIGNATURES,
    **_SOLID_TUMOR_PARENCHYMAL_SIGNATURES,
}

PBMC_CULTURE_TYPES: frozenset = frozenset({
    "T cell", "NK cell", "B cell", "Plasma cell",
    "Macrophage", "Monocyte", "Dendritic cell",
})

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

MOUSE_T_STATE_SIGNATURES: dict[str, list[str]] = {
    "naive/Tcm": ["CCR7", "SELL", "TCF7", "LEF1", "IL7R"],
    "effector/cytotoxic": ["GZMA", "GZMB", "GZMK", "GZMC", "PRF1", "NKG7",
                           "KLRG1", "FASL"],
    "Trm": ["ITGAE", "ZFP683", "CXCR6", "ITGA1", "CD69"],
    "Tex": ["PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "ENTPD1", "LAYN",
            "CCL3", "CCL4"],
    "Tfh": ["CXCR5", "BCL6", "ICOS", "TOX2", "CD200"],
    "Th1": ["TBX21", "IFNG", "CXCR3", "IL12RB2"],
    "Th2": ["GATA3", "IL13", "IL4", "IL5", "IL17RB"],
    "Th17": ["RORC", "IL17A", "IL17F", "CCR6", "IL23R", "IL22"],
    "proliferating": ["MKI67", "TOP2A", "STMN1", "TYMS"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "CCR8"],
    "IFN-stimulated": ["ISG15", "IFIT1", "IFIT3", "MX1", "OAS1A", "OAS1G"],
}

B_STATE_SIGNATURES: dict[str, list[str]] = {
    "naive B": ["TCL1A", "IGHD", "FCER2", "IL4R"],
    "memory B": ["CD27", "TNFRSF13B", "CD80", "AIM2"],
    "germinal-center B": ["BCL6", "RGS13", "AICDA", "MEF2B", "LRMP"],
}

MOUSE_B_STATE_SIGNATURES: dict[str, list[str]] = {
    "naive B": ["TCL1", "IGHD", "FCER2A", "IL4RA", "SELL"],
    "memory B": ["CD80", "PDCD1LG2", "NT5E", "CCR6", "ZBTB32", "AIM2", "TNFRSF13B"],
    "germinal-center B": ["BCL6", "RGS13", "AICDA", "MEF2B", "IRAG2"],
}

LINEAGE_GENES: dict[str, list[str]] = {
    "CD8": ["CD8A", "CD8B"],
    "CD4": ["CD4"],
    "CD3": ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "CD5"],
    "NKspec": ["GNLY", "KLRF1", "NCR1", "NCAM1", "FCGR3A", "KLRC1", "CD160"],
}

MOUSE_LINEAGE_GENES: dict[str, list[str]] = {
    "CD8": ["CD8A", "CD8B1"],
    "CD4": ["CD4"],
    "CD3": ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "CD5"],
    "NKspec": ["NCR1", "KLRB1C", "ITGA2", "KLRA4", "KLRA7", "KLRA8", "FCGR3"],
}

PROTECTED_SUBTYPES: frozenset = frozenset({"pDC", "cDC1"})


def cell_type_signatures_for_species(
    species: str = "human",
    *,
    solid_tumor: bool = False,
) -> dict[str, list[str]]:
    """Return cell-type signature registry for ``species``."""
    if normalize_species(species) == "mouse":
        return MOUSE_SOLID_TUMOR_CELL_TYPE_SIGNATURES if solid_tumor else MOUSE_CELL_TYPE_SIGNATURES
    return SOLID_TUMOR_CELL_TYPE_SIGNATURES if solid_tumor else CELL_TYPE_SIGNATURES


def t_state_signatures_for_species(species: str = "human") -> dict[str, list[str]]:
    """Return T-cell sub-state signature registry for ``species``."""
    return MOUSE_T_STATE_SIGNATURES if normalize_species(species) == "mouse" else T_STATE_SIGNATURES


def b_state_signatures_for_species(species: str = "human") -> dict[str, list[str]]:
    """Return B-cell sub-state signature registry for ``species``."""
    return MOUSE_B_STATE_SIGNATURES if normalize_species(species) == "mouse" else B_STATE_SIGNATURES


def lineage_genes_for_species(species: str = "human") -> dict[str, list[str]]:
    """Return CD8/CD4/CD3/NKspec lineage marker blocks for ``species``."""
    return MOUSE_LINEAGE_GENES if normalize_species(species) == "mouse" else LINEAGE_GENES


def marker_panel_for_species(species: str = "human") -> tuple[str, ...]:
    """Return default per-clone marker panel for ``species``."""
    return MARKER_PANEL_MGI if normalize_species(species) == "mouse" else MARKER_PANEL_HGNC
