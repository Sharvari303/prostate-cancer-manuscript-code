"""
config.py — Single source of truth for the PCa survival analysis pipeline.

All gene lists, cohort IDs, endpoints, thresholds, colors, and paths are
defined here. No other module should hardcode any of these values.

Biological axes (derived from SHAP analysis on PhysiCell ABM surrogate):
  Molecular axis   — PTEN, MDM2, MDM4, AR, CDKN1A, BAX, ATM
                     Alteration = CNA (GISTIC) + somatic mutation + SV
                     7 cohorts combined, OS endpoint

  Androgen axis    — SLCO2B1, SLCO1B3, AKR1C3
                     mRNA z-scores, per cohort independently

  Adhesion/motility axis — CDH1, EPCAM (epithelial, low is bad)
                           VIM, CDH2, MMP9, CXCL8, SNAI1, RHOA, ITGB1, FN1 (mesenchymal, high is bad)
                           mRNA z-scores, per cohort independently

Cohort pools:
  Molecular (7 cohorts, OS): TCGA_PRAD, SU2C, MCTP, MSK_2025, PIK3R1_MSK, MSK_2024, MCSPC_MSK
    MSKCC excluded (no OS); IDH_MUTANT excluded (biologically atypical)
  Expression (4 cohorts, per-cohort): TCGA_PRAD, SU2C, MSKCC, MCTP
    Endpoint: DFS for TCGA_PRAD and MSKCC; OS for SU2C and MCTP
    z-scores are never pooled across cohorts (different normalization baselines)
"""

from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────────────────────────

ROOT_DIR    = Path(__file__).resolve().parent.parent
DATA_DIR    = ROOT_DIR / "data"
CACHE_DIR   = DATA_DIR / "cache"
OUTPUT_DIR  = ROOT_DIR / "outputs"
FIGURES_DIR = OUTPUT_DIR / "figures"
TABLES_DIR  = OUTPUT_DIR / "tables"
LOG_DIR     = ROOT_DIR / "logs"

# Ensure all directories exist
for d in [CACHE_DIR, TABLES_DIR, LOG_DIR,
          FIGURES_DIR / "task1_molecular_individual",
          FIGURES_DIR / "task2_molecular_or",
          FIGURES_DIR / "task3_androgen_mrna_individual",
          FIGURES_DIR / "task3b_androgen_cna",
          FIGURES_DIR / "task4_androgen_or",
          FIGURES_DIR / "task6_adhesion_mrna_individual",
          FIGURES_DIR / "task6b_adhesion_cna",
          FIGURES_DIR / "task7_adhesion_or",
          FIGURES_DIR / "task_pfs_tcga",
          # Crowding / mechanobiology axes (module8) — mRNA expression, TCGA-PRAD only
          FIGURES_DIR / "task8_crowding_axis1",
          FIGURES_DIR / "task8_crowding_axis2",
          FIGURES_DIR / "task8_crowding_axis3",
          FIGURES_DIR / "task8_crowding_axis4",
          FIGURES_DIR / "task8_crowding_axis5",
          FIGURES_DIR / "task8_crowding_axis6",
          FIGURES_DIR / "task8_crowding_composite",
          FIGURES_DIR / "task8_crowding_forest"]:
    d.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# API
# ─────────────────────────────────────────────────────────────────────────────

CBIOPORTAL_BASE_URL = "https://www.cbioportal.org/api"
API_TIMEOUT_SECONDS = 30
API_RETRY_ATTEMPTS  = 3
API_RETRY_DELAY     = 2.0   # seconds between retries
API_POLITE_DELAY    = 0.5   # seconds between consecutive calls

# Per-study mRNA z-score profile overrides.
# Needed where the auto-selected profile (first match to "Zscore") is wrong:
#   MSKCC: "mrna_median_Zscores" is actually log2 raw expression (Affymetrix);
#          "mrna_zbynorm" is the true z-score profile (vs normal prostate tissue).
#   SU2C:  "mrna_seq_fpkm_capture_Zscores" is targeted capture (201 patients);
#          "mrna_seq_fpkm_polya_Zscores" is polyA whole-transcriptome (262 patients,
#          80 with OS data vs 65 from capture) — standard for CRPC cohorts.
MRNA_PROFILE_OVERRIDES = {
    "prad_mskcc":      "prad_mskcc_mrna_zbynorm",
    "prad_su2c_2019":  "prad_su2c_2019_mrna_seq_fpkm_polya_Zscores",
}


# ─────────────────────────────────────────────────────────────────────────────
# DATASETS
# ─────────────────────────────────────────────────────────────────────────────

DATASETS = {
    # ===== mRNA COHORTS — have both expression z-scores and survival =====
    # Analyzed per-cohort independently for expression axes (never pooled).
    # OS events are low in TCGA/MSKCC primary disease → use DFS as primary for those.
    "TCGA_PRAD": {
        "study_id":          "prad_tcga_pan_can_atlas_2018",
        "label":             "TCGA-PRAD",
        "stage":             "Primary",
        "n_expected":        494,
        "has_mrna":          True,
        # OS: 10/494 events (2%) — too few for OS KM; DFS: 30/494 (6%) — use as primary
        "primary_endpoint":  "DFS",
        "in_molecular_pool": True,
    },
    "SU2C": {
        "study_id":          "prad_su2c_2019",
        "label":             "SU2C/PCF CRPC",
        "stage":             "Metastatic CRPC",
        "n_expected":        444,
        "has_mrna":          True,
        # OS: 76/128 events (59%) — good OS data
        "primary_endpoint":  "OS",
        "in_molecular_pool": True,
    },
    "MSKCC": {
        "study_id":          "prad_mskcc",
        "label":             "MSKCC",
        "stage":             "Primary + Metastatic",
        "n_expected":        198,
        "has_mrna":          True,
        # OS: 0 events — DFS only; excluded from molecular cohort pool
        "primary_endpoint":  "DFS",
        "in_molecular_pool": False,   # no OS data — excluded from combined molecular KM
    },
    "MCTP": {
        "study_id":          "prad_mich",
        "label":             "MCTP Michigan",
        "stage":             "Metastatic",
        "n_expected":        119,
        "has_mrna":          True,    # 31 patients have both OS and mRNA z-scores
        # OS: 48/48 events (100%); mRNA subset n=31
        "primary_endpoint":  "OS",
        "in_molecular_pool": True,
    },
    # ===== CNA + MUTATION ONLY COHORTS — molecular axis pool =====
    "MSK_2025": {
        "study_id":          "prad_msk_2025",
        "label":             "MSK 2025",
        "stage":             "Mixed",
        "n_expected":        120,
        "has_mrna":          False,
        "primary_endpoint":  "OS",
        "in_molecular_pool": True,
    },
    "PIK3R1_MSK": {
        "study_id":          "prad_pik3r1_msk_2021",
        "label":             "PIK3R1 MSK 2021",
        "stage":             "Mixed",
        "n_expected":        1417,
        "has_mrna":          False,
        "primary_endpoint":  "OS",
        "in_molecular_pool": True,
    },
    "MSK_2024": {
        "study_id":          "prostate_msk_2024",
        "label":             "MSK Clin Cancer Res 2024",
        "stage":             "Mixed",
        "n_expected":        2260,
        "has_mrna":          False,
        "primary_endpoint":  "OS",
        "in_molecular_pool": True,
    },
    "MCSPC_MSK": {
        "study_id":          "prad_mcspc_mskcc_2020",
        "label":             "Metastatic CSPC MSK",
        "stage":             "Metastatic CSPC",
        "n_expected":        424,
        "has_mrna":          False,
        "primary_endpoint":  "OS",
        "in_molecular_pool": True,
    },
    # ===== EXCLUDED FROM MOLECULAR POOL — kept for reference/separate analyses =====
    "IDH_MUTANT": {
        "study_id":          "prad_idhmut_msk_2025",
        "label":             "IDH Mutant MSK 2024",
        "stage":             "IDH-mutant",
        "n_expected":        99,
        "has_mrna":          False,
        # Excluded: IDH-mutant PCa is a biologically distinct rare subtype;
        # pooling with general PCa would confound PTEN/MDM2/MDM4/AR alteration rates
        "primary_endpoint":  "OS",
        "in_molecular_pool": False,
    },
}

# Combined label - includes only cohorts with mRNA for expression analysis
COMBINED_LABEL = "TCGA-PRAD + SU2C + MSKCC + MCTP (mRNA cohorts)"
COMBINED_MOLECULAR_LABEL = "7-Cohort Molecular Pool (OS endpoint)"

# ─────────────────────────────────────────────────────────────────────────────
# FOCUSED COHORT SETS
# ─────────────────────────────────────────────────────────────────────────────

# 7-cohort molecular OS pool: all have CNA+mut data and OS endpoint
# Excludes: MSKCC (no OS), IDH_MUTANT (biologically atypical rare subtype)
MOLECULAR_COHORT_KEYS = [k for k, v in DATASETS.items() if v["in_molecular_pool"]]

# mRNA cohorts: analyzed independently per cohort, never pooled across studies.
# z-scores are normalized within each cohort — cross-study pooling is invalid.
# TCGA/MSKCC use DFS as primary (insufficient OS events); SU2C/MCTP use OS.
EXPRESSION_COHORT_KEYS = ["TCGA_PRAD", "SU2C", "MSKCC", "MCTP"]


# ─────────────────────────────────────────────────────────────────────────────
# SURVIVAL ENDPOINTS
# Per-dataset because TCGA and SU2C use different primary endpoints
# ─────────────────────────────────────────────────────────────────────────────

ENDPOINTS = {
    # TCGA-PRAD: primary disease, only 10/494 OS events (2%) — use DFS as primary.
    # OS is kept as secondary for completeness but should not be used for mRNA KM.
    "TCGA_PRAD": {
        "primary": {
            "time":   "DFS_MONTHS",
            "event":  "DFS_EVENT",
            "label":  "Disease-Free Survival",
            "short":  "DFS",
        },
        "secondary": {
            "time":   "OS_MONTHS",
            "event":  "OS_EVENT",
            "label":  "Overall Survival",
            "short":  "OS",
        },
    },
    # SU2C: metastatic CRPC, 76/128 OS events — OS is the appropriate endpoint
    "SU2C": {
        "primary": {
            "time":   "OS_MONTHS",
            "event":  "OS_EVENT",
            "label":  "Overall Survival",
            "short":  "OS",
        },
        "secondary": {
            "time":   "PFS_MONTHS",
            "event":  "PFS_EVENT",
            "label":  "Progression-Free Survival",
            "short":  "PFS",
        },
    },
    # MSKCC: no OS data — DFS only. Excluded from molecular OS pool.
    # Used only for mRNA expression DFS analysis.
    "MSKCC": {
        "primary": {
            "time":   "DFS_MONTHS",
            "event":  "DFS_EVENT",
            "label":  "Disease-Free Survival",
            "short":  "DFS",
        },
    },
    # MCTP: 48/48 OS events (100%); 31 patients have mRNA + OS
    "MCTP": {
        "primary": {
            "time":   "OS_MONTHS",
            "event":  "OS_EVENT",
            "label":  "Overall Survival",
            "short":  "OS",
        },
    },
    # Default endpoint for any cohort not explicitly listed (OS)
    "DEFAULT": {
        "primary": {
            "time":   "OS_MONTHS",
            "event":  "OS_EVENT",
            "label":  "Overall Survival",
            "short":  "OS",
        },
    },
    # MCSPC_MSK uses non-standard OS column names from cBioPortal
    "MCSPC_MSK": {
        "primary": {
            "time":   "OS_SMP_MONTHS",
            "event":  "OS_EVENT",
            "label":  "Overall Survival",
            "short":  "OS",
            "status_col": "SURVIVAL_STATUS",
        },
    },
    # Endpoint for 7-cohort combined molecular analysis (OS throughout)
    "MOLECULAR_COHORT_COMBINED": {
        "primary": {
            "time":   "OS_MONTHS",
            "event":  "OS_EVENT",
            "label":  "Overall Survival",
            "short":  "OS",
        },
    },
    "COMBINED_MOLECULAR": {
        "primary": {
            "time":   "OS_MONTHS",
            "event":  "OS_EVENT",
            "label":  "Overall Survival",
            "short":  "OS",
        },
    },
}

# Raw string values from cBioPortal that map to event = 1
EVENT_POSITIVE_STRINGS = {
    # MCSPC_MSK uses SURVIVAL_STATUS → renamed to OS_STATUS; its values are
    # plain "Dead" / "Alive" (no "0:"/"1:" prefix), so both must be listed here
    "OS_STATUS":      ["1:DECEASED", "DECEASED", "Dead", "DEAD"],
    "DFS_STATUS":     ["1:Recurred/Progressed", "1:RECURRENCE",
                       "Recurred", "Progressed"],
    "PFS_STATUS":     ["1:PROGRESSION",   "Progressed"],
}

EVENT_NEGATIVE_STRINGS = {
    # Explicit censored values — anything listed here maps to 0 instead of NaN.
    # Critical for MCSPC_MSK whose SURVIVAL_STATUS uses "Alive" (no "0:" prefix).
    "OS_STATUS":  ["0:LIVING", "LIVING", "Alive", "ALIVE"],
    "DFS_STATUS": ["0:DiseaseFree", "DiseaseFree"],
    "PFS_STATUS": ["0:CENSORED",   "CENSORED"],
}


# ─────────────────────────────────────────────────────────────────────────────
# ENTREZ GENE IDs
# Required by cBioPortal API — maps Hugo symbol → Entrez ID
# ─────────────────────────────────────────────────────────────────────────────

ENTREZ_IDS = {
    # Molecular axis
    "PTEN":    5728,
    "MDM2":    4193,
    "MDM4":    4194,
    "AR":       367,
    # AR activity axis
    "KLK3":    3816,
    "TMPRSS2": 7113,
    "NKX3.1":  4824,
    # Androgen uptake axis (mRNA + CNA analysis)
    "SLCO2B1": 11309,
    "SLCO1B3": 28234,
    "AKR1C3":  8644,
    # Adhesion-Motility axis (mRNA + CNA analysis)
    "CDH1":     999,   # E-cadherin; epithelial retention, loss is bad
    "EPCAM":   4072,   # Epithelial cell adhesion molecule; loss is bad
    "VIM":     7431,   # Vimentin; mesenchymal marker, gain is bad
    "CDH2":    1000,   # N-cadherin; mesenchymal switch, gain is bad
    "MMP9":    4318,   # Matrix metalloproteinase-9; invasion, gain is bad
    "CXCL8":   3576,   # IL-8; pro-invasive chemokine, gain is bad
    # ACTB removed: housekeeping gene, constitutively expressed,
    # z-score near zero in all samples — not an informative motility marker
    "SNAI1":   6615,   # Snail; canonical EMT transcription factor, represses CDH1, gain is bad
    "RHOA":     387,   # Rho GTPase A; contractility/invasion, gain is bad
    "ITGB1":   3688,   # Integrin beta-1; ECM attachment, gain is bad
    "FN1":     2335,   # Fibronectin-1; mesenchymal ECM, gain is bad
    # DNA damage response axis
    "CDKN1A":  1026,
    "BAX":      581,
    "ATM":      472,
    # Apoptosis / survival signaling axis
    "BCL2":     596,
    "RAF1":    5894,
    "CASP9":    842,
    "AKT1":     207,
    # ─── CROWDING / MECHANOBIOLOGY AXES (mRNA expression; see CROWDING_AXES) ───
    # Axis 1 — Hippo/YAP-TAZ
    "YAP1":   10413,
    "WWTR1":  25937,   # TAZ
    "LATS1":   9113,
    "LATS2":  26524,
    # Axis 2 — mechanosensitive ion channels
    "PIEZO1":  9780,
    "TRPV4":  59341,
    # Axis 3 — phosphoinositide kinases
    "PIP4K2B": 8396,
    "PIP5K1C":23396,
    "PIK3CA":  5290,
    # Axis 4 — cell cycle (CDKN1A 1026 already defined above — reused, not re-added)
    "CDKN1B":  1027,   # p27
    "CCND1":    595,   # Cyclin D1
    "MKI67":   4288,   # Ki-67
    # Axis 5 — hypoxia
    "HIF1A":   3091,
    "EPAS1":   2034,   # HIF2A
    "VEGFA":   7422,
    # Axis 6 — PIP2 / cytoskeletal mechanosensing (novel)
    "ANXA1":    301,
    "KPNA4":   3840,
    "RAE1":    8480,
    "PLS1":    5357,
    "ZYX":     7791,
    "ARF1":     375,
    "CDC42":    998,
    "EZR":     7430,
    "LAMP2":   3920,
    "HMOX1":   3162,
    "FLNC":    2318,
    "VAMP3":   9341,   # in PIP2 composite formula (not in brief per-gene table)
}


# ─────────────────────────────────────────────────────────────────────────────
# GENE AXIS DEFINITIONS
# Each gene includes its data type, biological mechanism, and direction of harm
# ─────────────────────────────────────────────────────────────────────────────

MOLECULAR_GENES = {
    # DNA-level structural events — permanent, hardwired genomic lesions
    # Explicit separation: CNA (copy number) vs Mutation (point variants)
    # Different biological origins, detection properties, and frequencies
    # Mechanism links to cellular MHS model SHAP top features
    "PTEN": {
        "alteration_type": "suppressor",
        "mechanism":      "PI3K_AKT_dysregulation",
        "direction":      "loss_is_bad",
        "shap_rank":      "Top 3 all models",
        "note":           "Primary SA variable in cellular model",
        "confidence":     "very_high",
        # ─── DATA SOURCES: Explicit hierarchy ───
        "data_sources": {
            "CNA": {
                "primary":     True,      # CNA is main driver for PTEN loss
                "description": "Copy number alterations (GISTIC scores)",
                "thresholds":  {
                    "strict":  -2,  # homozygous deletion only
                    "relaxed": -1,  # includes heterozygous deletion
                },
            },
            "mutation": {
                "primary":     False,     # secondary
                "description": "Somatic point mutations",
                "types":       ["Nonsense_Mutation", "Frame_Shift_Del",
                               "Frame_Shift_Ins", "Splice_Site"],
            },
            "sv": {
                "primary":     False,
                "description": "Structural variants (fusions, translocations)",
            },
            "mrna": {
                "primary":     False,
                "description": "mRNA expression z-score (±2 threshold)",
            },
        },
        "alteration_logic": "CNA_OR_mutation_OR_sv_OR_mrna",
    },
    "MDM2": {
        "alteration_type": "oncogene",
        "mechanism":      "p53_suppression",
        "direction":      "gain_is_bad",
        "shap_rank":      "Top 5 RF and SVM models",
        "note":           "Co-alteration with PTEN drives double-hit p53/AKT",
        "confidence":     "high",
        # ─── DATA SOURCES: Explicit hierarchy ───
        "data_sources": {
            "CNA": {
                "primary":     True,      # CNA is main driver for MDM2 gain
                "description": "Copy number alterations (GISTIC scores)",
                "thresholds":  {
                    "strict":  +2,
                    "relaxed": +1,
                },
            },
            "mutation": {
                "primary":     False,     # mutations rare in MDM2
                "description": "Somatic point mutations",
                "types":       ["Missense_Mutation"],
            },
            "sv": {
                "primary":     False,
                "description": "Structural variants",
            },
        },
        "alteration_logic": "CNA_OR_mutation_OR_sv",
    },
    "MDM4": {
        "alteration_type": "oncogene",
        "mechanism":      "p53_transcription_inhibition",
        "direction":      "gain_is_bad",
        "shap_rank":      "#1 in NN model",
        "note":           "Top SHAP feature — borderline KM expected "
                          "in primary disease, stronger in CRPC",
        "confidence":     "high",
        # ─── DATA SOURCES: Explicit hierarchy ───
        "data_sources": {
            "CNA": {
                "primary":     True,      # CNA is main driver
                "description": "Copy number alterations (GISTIC scores)",
                "thresholds":  {
                    "strict":  +2,
                    "relaxed": +1,
                },
            },
            "mutation": {
                "primary":     False,
                "description": "Somatic point mutations",
                "types":       ["Missense_Mutation"],
            },
            "sv": {
                "primary":     False,
                "description": "Structural variants",
            },
        },
        "alteration_logic": "CNA_OR_mutation_OR_sv",
    },
    "AR": {
        "alteration_type": "oncogene",
        "mechanism":      "androgen_independent_activation",
        "direction":      "gain_is_bad",
        "shap_rank":      "Top 3 all models",
        "note":           "Rare in primary TCGA (~1% CNA), ~42-52% in CRPC. "
                          "CNA amplification is the dominant and unambiguous mechanism. "
                          "Mutations excluded: protein-change data not stored in master; "
                          "mutations add <5% patients beyond CNA alone.",
        "confidence":     "high",
        # ─── DATA SOURCES ───
        "data_sources": {
            "CNA": {
                "primary":     True,      # Amplification is the main, high-frequency AR driver
                "description": "Copy number amplification (GISTIC scores)",
                "thresholds":  {
                    "strict":  +2,   # homozygous/high-level amplification
                    "relaxed": +1,   # includes low-level gain
                },
            },
        },
        "alteration_logic": "CNA_only",
    },
    "BCL2": {
        "alteration_type": "oncogene",
        "mechanism":      "apoptosis_suppression",
        "direction":      "gain_is_bad",
        "shap_rank":      "Apoptosis axis",
        "note":           "Anti-apoptotic; overexpression blocks intrinsic apoptosis",
        "confidence":     "high",
        "data_sources": {
            "CNA": {
                "primary":     True,
                "description": "Copy number alterations (GISTIC scores)",
                "thresholds":  {"strict": +2, "relaxed": +1},
            },
            "mutation": {
                "primary":     False,
                "description": "Somatic point mutations",
                "types":       ["Missense_Mutation"],
            },
        },
        "alteration_logic": "CNA_OR_mutation",
    },
    "RAF1": {
        "alteration_type": "oncogene",
        "mechanism":      "MAPK_pathway_activation",
        "direction":      "gain_is_bad",
        "shap_rank":      "MAPK/ERK axis",
        "note":           "RAF1 amplification activates MAPK/ERK proliferation signaling",
        "confidence":     "high",
        "data_sources": {
            "CNA": {
                "primary":     True,
                "description": "Copy number alterations (GISTIC scores)",
                "thresholds":  {"strict": +2, "relaxed": +1},
            },
            "mutation": {
                "primary":     False,
                "description": "Somatic point mutations",
                "types":       ["Missense_Mutation"],
            },
        },
        "alteration_logic": "CNA_OR_mutation",
    },
    "CASP9": {
        "alteration_type": "suppressor",
        "mechanism":      "intrinsic_apoptosis_initiation",
        "direction":      "loss_is_bad",
        "shap_rank":      "Apoptosis axis",
        "note":           "Initiator caspase; loss blocks mitochondrial apoptosis pathway",
        "confidence":     "high",
        "data_sources": {
            "CNA": {
                "primary":     True,
                "description": "Copy number alterations (GISTIC scores)",
                "thresholds":  {"strict": -2, "relaxed": -1},
            },
            "mutation": {
                "primary":     False,
                "description": "Somatic point mutations",
                "types":       ["Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins"],
            },
        },
        "alteration_logic": "CNA_OR_mutation",
    },
    "AKT1": {
        "alteration_type": "oncogene",
        "mechanism":      "PI3K_AKT_survival_signaling",
        "direction":      "gain_is_bad",
        "shap_rank":      "PI3K/AKT axis",
        "note":           "Downstream of PTEN; AKT1 amplification/mutation drives survival and proliferation",
        "confidence":     "high",
        "data_sources": {
            "CNA": {
                "primary":     True,
                "description": "Copy number alterations (GISTIC scores)",
                "thresholds":  {"strict": +2, "relaxed": +1},
            },
            "mutation": {
                "primary":     False,
                "description": "Somatic point mutations",
                "types":       ["Missense_Mutation"],
            },
        },
        "alteration_logic": "CNA_OR_mutation",
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# AR ACTIVITY AXIS (NEW)
# Measures functional androgen receptor signaling (independent of AR alterations)
# Data: mRNA z-scores of AR target genes
# Maps to: AR activation / proliferation rate SA in PhysiCell ABM
# ─────────────────────────────────────────────────────────────────────────────
# CRITICAL: Many CRPC tumors have WT AR + very high AR activity
# AR alterations (CNA/mutation) ≠ AR activity (transcriptional output)
# Literature: Cato et al., Takeda et al. show AR activity > AR status in CRPC

AR_ACTIVITY_GENES = {
    "KLK3": {  # Prostate-specific antigen (PSA)
        "data_type":  "mRNA_expression",
        "mechanism":  "ar_target_gene_canonical",
        "direction":  "high_is_bad",
        "note":       "Canonical AR target gene. mRNA levels proxy for serum PSA. "
                      "High = strong AR-driven proliferation and secretion",
        "confidence": "very_high",
        "evidence":   ["Classical AR target in prostate cancer",
                       "Measured clinically as serum PSA",
                       "Strong association with AR activity"],
    },
    "TMPRSS2": {
        "data_type":  "mRNA_expression",
        "mechanism":  "ar_target_gene_fusion_hotspot",
        "direction":  "high_is_bad",
        "note":       "TMPRSS2-ERG fusions are oncogenic; mRNA high indicates "
                      "either AR-driven transcription or fusion activity. "
                      "Proxy for AR-regulated genes",
        "confidence": "very_high",
        "evidence":   ["TMPRSS2-ERG fusion common in prostate cancer",
                       "High mRNA = strong AR activity or fusion-driven",
                       "Drives ETS transcriptional program"],
    },
    "NKX3.1": {
        "data_type":  "mRNA_expression",
        "mechanism":  "ar_target_gene_paradoxical",
        "direction":  "high_is_bad",
        "note":       "Paradoxical: high in CRPC despite androgen deprivation. "
                      "Indicates sustained AR activity independent of circulating androgens. "
                      "Marker of intracrine/paracrine AR signaling",
        "confidence": "high",
        "evidence":   ["AR target in normal prostate",
                       "Retained/upregulated in CRPC despite ADT",
                       "Indicates AR axis not shut down in resistant disease"],
    },
}

ANDROGEN_GENES = {
    # Androgen uptake and synthesis — transcriptional adaptation to ADT
    # Distinct from AR ACTIVITY (see AR_ACTIVITY_GENES above)
    # Data: mRNA z-scores (RNA-seq)
    # Maps to: testosterone uptake rate SA in PhysiCell ABM
    # Key distinction:
    #   AR_ACTIVITY_GENES → Measures AR signaling output (KLK3, TMPRSS2, NKX3.1)
    #   ANDROGEN_GENES → Measures AR input/supply (SLCO, AKR)
    "SLCO2B1": {
        "data_type":  "mRNA_expression",
        "mechanism":  "androgen_precursor_transport",
        "direction":  "high_is_bad",
        "note":       "Solute carrier organic anion transporter family member 2B1. "
                      "Upregulated under ADT; proxy for high cellular uptake rate in ABM",
        "confidence": "high",
    },
    "SLCO1B3": {
        "data_type":  "mRNA_expression",
        "mechanism":  "testosterone_transport",
        "direction":  "high_is_bad",
        "note":       "Testosterone transporter; CRPC enriched. "
                      "Enables import of circulating androgens in castrate setting",
        "confidence": "high",
    },
    "AKR1C3": {
        "data_type":  "mRNA_expression",
        "mechanism":  "intratumoral_androgen_synthesis",
        "direction":  "high_is_bad",
        "note":       "Aldo-keto reductase converts 4-androstenedione → testosterone. "
                      "Intracrine synthesis circumvents ADT. Cited in manuscript",
        "confidence": "high",
    },
}

ADHESION_MOTILITY_GENES = {
    # 10-gene adhesion/motility axis — mRNA z-score analysis per cohort
    # Epithelial retention (LOW is bad — loss of epithelial identity promotes invasion)
    "CDH1": {
        "data_type": "mRNA_expression",
        "direction": "low_is_bad",
        "sub_group": "epithelial",
        "note":      "E-cadherin; loss drives PCa invasion",
    },
    "EPCAM": {
        "data_type": "mRNA_expression",
        "direction": "low_is_bad",
        "sub_group": "epithelial",
        "note":      "Epithelial cell adhesion molecule; loss promotes EMT",
    },
    # Mesenchymal/invasion (HIGH is bad — promotes invasion and metastasis)
    "VIM": {
        "data_type": "mRNA_expression",
        "direction": "high_is_bad",
        "sub_group": "mesenchymal",
        "note":      "Vimentin; mesenchymal intermediate filament marker",
    },
    "CDH2": {
        "data_type": "mRNA_expression",
        "direction": "high_is_bad",
        "sub_group": "mesenchymal",
        "note":      "N-cadherin; cadherin switch drives metastasis",
    },
    "MMP9": {
        "data_type": "mRNA_expression",
        "direction": "high_is_bad",
        "sub_group": "mesenchymal",
        "note":      "Matrix metalloproteinase-9; ECM degradation and invasion",
    },
    "CXCL8": {
        "data_type": "mRNA_expression",
        "direction": "high_is_bad",
        "sub_group": "mesenchymal",
        "note":      "IL-8; pro-invasive chemokine; promotes angiogenesis and invasion",
    },
    # ACTB removed — housekeeping gene; constitutively expressed; z-score near 0
    # in all samples; not a selective motility marker.
    "SNAI1": {
        "data_type": "mRNA_expression",
        "direction": "high_is_bad",
        "sub_group": "mesenchymal",
        "note":      "Snail transcription factor; canonical EMT driver; directly represses CDH1",
    },
    "RHOA": {
        "data_type": "mRNA_expression",
        "direction": "high_is_bad",
        "sub_group": "mesenchymal",
        "note":      "Rho GTPase A; promotes actomyosin contractility and invasion",
    },
    "ITGB1": {
        "data_type": "mRNA_expression",
        "direction": "high_is_bad",
        "sub_group": "mesenchymal",
        "note":      "Integrin beta-1; ECM attachment and PCa bone metastasis",
    },
    "FN1": {
        "data_type": "mRNA_expression",
        "direction": "high_is_bad",
        "sub_group": "mesenchymal",
        "note":      "Fibronectin-1; mesenchymal ECM component; promotes invasion",
    },
}

# Sub-group lists for adhesion/motility axis (mRNA)
# ACTB removed (housekeeping gene); SNAI1 added (canonical EMT driver)
ADHESION_EPITHELIAL_GENES   = ["CDH1", "EPCAM"]
ADHESION_MESENCHYMAL_GENES  = ["VIM", "CDH2", "MMP9", "CXCL8", "SNAI1", "RHOA", "ITGB1", "FN1"]

# Individual genes (CNA/mutation analysis — molecular axis extensions)
INDIVIDUAL_GENES = {
    "CDKN1A": {
        "alteration_type": "suppressor",
        "direction":       "loss_is_bad",
        "description":     "p21/WAF1; p53-induced cell cycle inhibitor",
    },
    "BAX": {
        "alteration_type": "suppressor",
        "direction":       "loss_is_bad",
        "description":     "BCL2 associated X; p53-induced apoptosis regulator",
    },
    "ATM": {
        "alteration_type": "suppressor",
        "direction":       "loss_is_bad",
        "description":     "Ataxia telangiectasia mutated kinase; DNA damage sensor",
    },
}

# Androgen genes analyzed at DNA level (CNA + mut + SV) on 7-cohort molecular pool — Task 3b
ANDROGEN_CNA_GENES = {
    "SLCO2B1": {"alteration_type": "oncogene", "direction": "gain_is_bad"},
    "SLCO1B3": {"alteration_type": "oncogene", "direction": "gain_is_bad"},
    "AKR1C3":  {"alteration_type": "oncogene", "direction": "gain_is_bad"},
}

# Adhesion/motility genes analyzed at DNA level (CNA + mut + SV) on 7-cohort molecular pool — Task 6b
# ACTB removed (housekeeping, CNA/mut not informative); SNAI1 added
ADHESION_CNA_GENES = {
    "CDH1":  {"alteration_type": "suppressor", "direction": "loss_is_bad"},
    "EPCAM": {"alteration_type": "suppressor", "direction": "loss_is_bad"},
    "VIM":   {"alteration_type": "oncogene",   "direction": "gain_is_bad"},
    "CDH2":  {"alteration_type": "oncogene",   "direction": "gain_is_bad"},
    "MMP9":  {"alteration_type": "oncogene",   "direction": "gain_is_bad"},
    "CXCL8": {"alteration_type": "oncogene",   "direction": "gain_is_bad"},
    "SNAI1": {"alteration_type": "oncogene",   "direction": "gain_is_bad"},
    "RHOA":  {"alteration_type": "oncogene",   "direction": "gain_is_bad"},
    "ITGB1": {"alteration_type": "oncogene",   "direction": "gain_is_bad"},
    "FN1":   {"alteration_type": "oncogene",   "direction": "gain_is_bad"},
}

# All gene lists for convenience
ALL_MOLECULAR_GENES         = list(MOLECULAR_GENES.keys())
ALL_AR_ACTIVITY_GENES       = list(AR_ACTIVITY_GENES.keys())
ALL_ANDROGEN_GENES          = list(ANDROGEN_GENES.keys())
ALL_ADHESION_MOTILITY_GENES = list(ADHESION_MOTILITY_GENES.keys())
ALL_ADHESION_GENES          = ALL_ADHESION_MOTILITY_GENES
ALL_INDIVIDUAL_GENES        = list(INDIVIDUAL_GENES.keys())
ALL_ANDROGEN_CNA_GENES      = list(ANDROGEN_CNA_GENES.keys())
ALL_ADHESION_CNA_GENES      = list(ADHESION_CNA_GENES.keys())
ALL_EXPRESSION_GENES        = ALL_ANDROGEN_GENES + ALL_ADHESION_MOTILITY_GENES
# All genes needing CNA/mut downloads (molecular + individual + androgen CNA + adhesion CNA)
ALL_CNA_GENES = (list(MOLECULAR_GENES.keys()) + ALL_INDIVIDUAL_GENES +
                 ALL_ANDROGEN_CNA_GENES + ALL_ADHESION_CNA_GENES)


# ─────────────────────────────────────────────────────────────────────────────
# CROWDING / MECHANOBIOLOGY AXES (Axes 1–6)
# Clinical anchor for the Fig 5 crowding result (PTEN-deleted clones proliferate
# under physical crowding). Source: PI brief PCa_Crowding_6Axes_KM_Brief.pdf,
# built from Nukpezah lab thesis Ch.3 (2025, "Nukpezah et al., in preparation").
#
# DESIGN (decided 2026-06-14):
#   - ALL six axes are mRNA EXPRESSION (RNA-seq z-scores), NOT DNA/genetic calls.
#     They route through module5_flag_expression, never module4.
#   - The ONLY DNA variable is the PTEN stratifier: deep deletion (CNA == -2),
#     column PTEN_DEEPDEL (NOT PTEN_ALT_STRICT, which folds in mutation/SV).
#   - Cohort: TCGA-PRAD only. Endpoints: DFS primary + PFS secondary, NO OS
#     (TCGA-PRAD has ~10/494 OS events = 2%, effectively event-free).
#   - Stratified evidence = Cox gene_z * PTEN_DEEPDEL interaction term (all ~490
#     pts); PTEN-stratified KM curves are descriptive only (flag arms < 20).
#   - Dual-data genes CDKN1A/AR: crowding uses mRNA ONLY; their existing CNA
#     analyses are untouched and coexist as distinct columns.
#   - Splits: median default; quartile only for ANXA1 and the PIP2 composite.
#   - FDR: Benjamini-Hochberg across Axis 6 genes, q < 0.10 (exploratory).
#
# direction:  high_is_bad | low_is_bad
# pten_stratify: whether to fit the gene * PTEN_DEEPDEL interaction
# ─────────────────────────────────────────────────────────────────────────────

CROWDING_AXES = {
    # ===== Axis 1 — Hippo/YAP-TAZ (contact inhibition master regulator) =====
    "YAP1":    {"axis": 1, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "published",
                "note": "Nuclear YAP overgrowth; AKT inhibits LATS under PTEN loss"},
    "WWTR1":   {"axis": 1, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "published",
                "note": "TAZ; YAP paralog, activated by PTEN loss"},
    "LATS1":   {"axis": 1, "data_type": "mRNA_expression", "direction": "low_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "published",
                "note": "Low LATS1 = YAP constitutively active = proliferates through crowding"},
    "LATS2":   {"axis": 1, "data_type": "mRNA_expression", "direction": "low_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "published",
                "note": "As LATS1; interaction with PTEN loss"},

    # ===== Axis 2 — Mechanosensitive ion channels =====
    "PIEZO1":  {"axis": 2, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "seminovel",
                "note": "Compressive mechanosensor; novel in PCa (known in breast)"},
    "TRPV4":   {"axis": 2, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "seminovel",
                "note": "Activates PI3K/AKT->p27 bypass; amplified in PTEN loss"},

    # ===== Axis 3 — Phosphoinositide kinases =====
    "PIP4K2B": {"axis": 3, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "seminovel",
                "note": "Mechanoresponsive PIP2 kinase; YAP nuclear retention"},
    "PIP5K1C": {"axis": 3, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "seminovel",
                "note": "Generates PIP2 substrate pool; max PIP2->PIP3 in PTEN loss"},
    "PIK3CA":  {"axis": 3, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "context",
                "note": "Known oncogene; PI3K-AKT context for crowding resistance"},

    # ===== Axis 4 — Cell cycle arrest genes =====
    "CDKN1B":  {"axis": 4, "data_type": "mRNA_expression", "direction": "low_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "published",
                "note": "p27; crowding arrest bypassed in PTEN-deleted"},
    "CDKN1A":  {"axis": 4, "data_type": "mRNA_expression", "direction": "low_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "published",
                "note": "p21 (mRNA expression — distinct from molecular CNA CDKN1A)"},
    "CCND1":   {"axis": 4, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "published",
                "note": "Cyclin D1; maintained proliferation despite crowding"},
    "MKI67":   {"axis": 4, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": False, "split": "median",   "priority": "published",
                "note": "Ki-67; proliferation phenotype confirmation, no PTEN strat"},

    # ===== Axis 5 — Hypoxia =====
    "HIF1A":   {"axis": 5, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "published",
                "note": "Crowded-tumor hypoxia; PTEN loss prevents HIF1A arrest"},
    "EPAS1":   {"axis": 5, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "published",
                "note": "HIF2A; relevant to CRPC"},
    "VEGFA":   {"axis": 5, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": False, "split": "median",   "priority": "published",
                "note": "Angiogenic output of crowding/hypoxia; no PTEN strat"},

    # ===== Axis 6 — PIP2 / cytoskeletal mechanosensing (NOVEL — headline) =====
    "ANXA1":   {"axis": 6, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "quartile", "priority": "highest",
                "novel": True, "note": "Dominant PIP2 mechanotype predictor; NO PCa data"},
    "KPNA4":   {"axis": 6, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "highest",
                "novel": True, "note": "p=1.7e-7 HR=3.09 across cancers; imports crowding TFs"},
    "RAE1":    {"axis": 6, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "highest",
                "novel": True, "note": "p=2.6e-7 HR=3.25 across cancers; mRNA export"},
    "PLS1":    {"axis": 6, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "high",
                "novel": True, "note": "PIP2-dependent actin bundler"},
    "ZYX":     {"axis": 6, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "high",
                "novel": True, "note": "Focal adhesion force sensor"},
    "ARF1":    {"axis": 6, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "high",
                "novel": True, "note": "PIP2-Golgi coupling; mechanosensitive trafficking"},
    "CDC42":   {"axis": 6, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "high",
                "novel": True, "note": "PIP2-activated actin nucleation; maintained under crowding"},
    "EZR":     {"axis": 6, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": False, "split": "median",   "priority": "medium",
                "novel": True, "note": "PIP2-activated membrane-cytoskeleton linker"},
    "LAMP2":   {"axis": 6, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "medium",
                "novel": True, "note": "Lysosomal survival under crowding+hypoxia (22RV1 Fig 3.12)"},
    "HMOX1":   {"axis": 6, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": True,  "split": "median",   "priority": "medium",
                "novel": True, "note": "Oxidative stress defense; crowding->hypoxia->HMOX1 in CRPC"},
    "FLNC":    {"axis": 6, "data_type": "mRNA_expression", "direction": "high_is_bad",
                "pten_stratify": False, "split": "median",   "priority": "medium",
                "novel": True, "note": "Actin cross-linker / force sensor"},
}

# Cohort and endpoint scope for crowding axes
CROWDING_COHORT_KEY = "TCGA_PRAD"          # single cohort (n≈494)
CROWDING_ENDPOINTS  = ["DFS", "PFS"]       # primary, secondary — NO OS (event-free)
CROWDING_FDR_METHOD = "fdr_bh"             # Benjamini-Hochberg
CROWDING_FDR_Q      = 0.10                 # exploratory discovery threshold

# Convenience lists
ALL_CROWDING_GENES = list(CROWDING_AXES.keys())
CROWDING_AXIS6_GENES = [g for g, v in CROWDING_AXES.items() if v["axis"] == 6]
CROWDING_GENES_PTEN_STRATIFIED = [g for g, v in CROWDING_AXES.items() if v["pten_stratify"]]


# ─────────────────────────────────────────────────────────────────────────────
# CNA THRESHOLDS
# ─────────────────────────────────────────────────────────────────────────────

CNA_THRESHOLDS = {
    "strict":  {"suppressor": -2, "oncogene": +2},
    "relaxed": {"suppressor": -1, "oncogene": +1},
}

# Default threshold used in main analysis
DEFAULT_CNA_THRESHOLD = "strict"


# ─────────────────────────────────────────────────────────────────────────────
# EXPRESSION SPLIT STRATEGIES
# How to split continuous mRNA z-scores into High / Low groups
# ─────────────────────────────────────────────────────────────────────────────

EXPRESSION_SPLITS = {
    "median": {
        "label":       "Median split",
        "description": "High = above median, Low = below median",
        "primary":     True,   # used in main figures
    },
    "quartile": {
        "label":       "Quartile split",
        "description": "High = top 25%, Low = bottom 25%",
        "primary":     False,  # used in sensitivity check figures
        "upper_q":     0.75,
        "lower_q":     0.25,
    },
}

DEFAULT_EXPRESSION_SPLIT = "median"


# ─────────────────────────────────────────────────────────────────────────────
# COMPOSITE AXIS SCORE DEFINITIONS
# How to combine individual gene z-scores into a single axis score
# Direction: positive score = more aggressive
# ─────────────────────────────────────────────────────────────────────────────

COMPOSITE_SCORES = {
    "ar_activity_score": {
        "label":       "AR Signaling Activity Score",
        "description": "Mean z-score of KLK3 + TMPRSS2 + NKX3.1. "
                       "Measures AR transcriptional activity independent of AR alterations. "
                       "High = strong AR-driven proliferation. "
                       "Proxy for proliferation rate SA in ABM.",
        "genes":       ["KLK3", "TMPRSS2", "NKX3.1"],
        "weights":     [1, 1, 1],       # equal weighting
        "weighting_justification": (
            "Equal weights maintain interpretability and avoid overfitting. "
            "Each gene represents distinct aspects of AR activity: "
            "KLK3 (canonical target), TMPRSS2 (fusion hotspot/AR target), "
            "NKX3.1 (paradoxical elevation in CRPC). Unbiased by cohort-specific variance."
        ),
        "direction":   "high_is_bad",
        "note":        "CRITICAL: AR activity ≠ AR alterations. Many CRPC tumors have "
                       "WT AR + very high activity. Expected to outperform AR CNA/mutation "
                       "in SU2C cohort.",
    },
    "androgen_uptake_score": {
        "label":       "Androgen Uptake Score",
        "description": "Mean z-score of SLCO2B1 + SLCO1B3 + AKR1C3. "
                       "Higher = more aggressive androgen sequestration (uptake + synthesis). "
                       "Proxy for high cellular androgen uptake rate in ABM.",
        "genes":       ["SLCO2B1", "SLCO1B3", "AKR1C3"],
        "weights":     [1, 1, 1],       # equal weighting
        "weighting_justification": (
            "Equal weights for uptake genes: SLCO2B1 (adrenal precursor), "
            "SLCO1B3 (testosterone), AKR1C3 (intracrine synthesis). "
            "Each represents distinct metabolic step."
        ),
        "direction":   "high_is_bad",
        "note":        "Androgen supply axis (how much androgen entering/made). "
                       "Distinct from AR_ACTIVITY (how much AR is active).",
    },
    "adhesion_motility_score": {
        "label":       "Adhesion & Motility Score",
        "description": "mesenchymal_mean - epithelial_mean. Higher = more invasive phenotype.",
        # epithelial_retention = mean(CDH1, EPCAM) — subtracted
        "genes_epithelial":  ["CDH1", "EPCAM"],
        # mesenchymal/invasion = mean(VIM, CDH2, MMP9, CXCL8, SNAI1, RHOA, ITGB1, FN1) — added
        "genes_mesenchymal": ["VIM", "CDH2", "MMP9", "CXCL8", "SNAI1", "RHOA", "ITGB1", "FN1"],
        "direction":   "high_is_bad",
    },
    "combined_tme_score": {
        "label":       "Combined TME Aggressiveness Score",
        "description": "AR Activity + Androgen Uptake + Adhesion-Motility Score. "
                       "Captures compounded proliferative + metabolic + mechanical "
                       "TME pressure modeled in ABM.",
        "components":  ["ar_activity_score", "androgen_uptake_score", "adhesion_motility_score"],
        "direction":   "high_is_bad",
        "note":        "Composite of all three functional axes. Higher = more aggressive phenotype.",
    },
    "pip2_trafficking_score": {
        "label":       "PIP2 Trafficking Score",
        "description": "Mean z-score of ANXA1 + ARF1 + CDC42 + EZR + VAMP3. "
                       "Nukpezah thesis PIP2 interactome signature (top 5). "
                       "Dichotomized top vs bottom quartile. Higher = crowding-resistant "
                       "PIP2-mediated mechanical adaptation. Axis 6 headline result.",
        "genes":       ["ANXA1", "ARF1", "CDC42", "EZR", "VAMP3"],
        "weights":     [1, 1, 1, 1, 1],   # equal weighting
        "weighting_justification": (
            "Equal weights for the top-5 PIP2 interactome proteins from the thesis "
            "mutual-information feature selection. VAMP3 retained per composite formula "
            "though absent from the brief's per-gene KM table."
        ),
        "split":       "quartile",        # top vs bottom quartile per brief
        "direction":   "high_is_bad",
        "pten_stratify": True,
        "note":        "MOST NOVEL: no published study links PIP2 interactome expression "
                       "to PCa survival. Tested with gene * PTEN_DEEPDEL interaction.",
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# STATISTICAL SETTINGS
# ─────────────────────────────────────────────────────────────────────────────

STATS = {
    "min_group_size":          20,    # minimum patients per KM arm
    "min_cox_events":          10,    # minimum events for stable Cox HR estimate
    "significance_threshold":  0.05,
    "pairwise_correction":     "bonferroni",
    "confidence_interval":     0.95,
}

# Age column names used by each cohort in cBioPortal raw data → standardized to AGE
AGE_COLUMN_ALIASES = {
    "AGE_AT_DIAGNOSIS":        "AGE",   # SU2C
    "AGE_CURRENT":             "AGE",   # PIK3R1_MSK
    "CURRENT_AGE":             "AGE",   # MSK_2024
    "CURRENT_AGE_DEID":        "AGE",   # IDH_MUTANT
    "AGE_AT_START_OF_PARPI":   "AGE",   # MSK_2025
    "AGE_AT_SAMPLE_COLLECTION":"AGE",   # MCSPC_MSK
    "DXAGE":                   "AGE",   # IDH_MUTANT alternate
}


# ─────────────────────────────────────────────────────────────────────────────
# PLOT AESTHETICS
# Consistent across all figures
# ─────────────────────────────────────────────────────────────────────────────

COLORS = {
    # Binary groups (altered/high vs unaltered/low)
    "altered":     "#E63946",
    "unaltered":   "#457B9D",
    "high":        "#E63946",
    "low":         "#2A9D8F",

    # Composite risk score groups
    "score_0":     "#2196F3",
    "score_1":     "#FF9800",
    "score_2plus": "#F44336",

    # Per-axis accent colors
    "molecular":   "#6D2B7E",
    "androgen":    "#E76F51",
    "adhesion":    "#2A9D8F",
    "combined":    "#9B2226",

    # Co-alteration subgroups
    "no_alt":          "#2196F3",
    "pten_only":       "#4CAF50",
    "mdm_only":        "#FF9800",
    "pten_mdm":        "#F44336",
    "ar_other":        "#9C27B0",

    # Dataset markers
    "tcga":  "#1565C0",
    "su2c":  "#B71C1C",
}

PLOT_STYLE = {
    "figure_dpi":       300,
    "font_family":      "DejaVu Sans",   # Arial on Windows/Mac locally
    "title_fontsize":   18,     # Figure titles (13 → 18, +30%)
    "label_fontsize":   17,     # X/Y axis labels (13 → 17, +30%)
    "tick_fontsize":    16,     # Axis tick numbers (12 → 16, +30%)
    "legend_fontsize":  14,     # Legend text (11 → 14, +30%)
    "linewidth":        2.5,    # KM curve thickness (increased from 2.0)
    "ci_alpha":         0.15,    # confidence interval shading transparency
    "formats":          ["pdf", "svg"],   # output formats per figure
}

