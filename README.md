# PCa Survival Analysis — KM Pipeline

Kaplan-Meier and Cox Proportional Hazards survival analysis across prostate cancer cohorts sourced from cBioPortal. Analyzes genetic/intrinsic and extrinsic heterogeneity drivers: (1) molecular alterations, (2a) androgen uptake/synthesis, (2b) adhesion/motility regulators, and (3) crowding / mechanobiology axes.

---

## Cohorts and Genes by Axis

### Axis 1: Molecular Alterations (Intrinsic Heterogeneity)
**Cohorts:** TCGA-PRAD (primary PCa), SU2C/PCF CRPC (metastatic CRPC), MCTP Michigan (primary PCa), MSK 2025, PIK3R1 MSK 2021, MSK Clin Cancer Res 2024, Metastatic CSPC MSK (metastatic CSPC)

**Genes:** PTEN, MDM2, MDM4, AR, CDKN1A, BAX, ATM

### Axis 2a: Androgen Uptake & Synthesis (Extrinsic Heterogeneity)
**Cohorts:** TCGA-PRAD (primary PCa), MSKCC (primary PCa), SU2C/PCF CRPC (metastatic CRPC)

**Genes:** SLCO2B1, SLCO1B3, AKR1C3

### Axis 2b: Adhesion / Motility/ (Extrinsic Heterogeneity)
**Cohorts:** TCGA-PRAD (primary PCa), MSKCC (primary PCa), SU2C/PCF CRPC (metastatic CRPC)

**Genes:** CDH1, EPCAM, VIM, CDH2, MMP9, CXCL8, SNAI1, RHOA, ITGB1, FN1

### Axis 3: Crowding / Mechanobiology (clinical anchor for Fig 5 crowding result)
**Cohort:** TCGA-PRAD (primary PCa)

**Genes (6 sub-axes):** Hippo/YAP (YAP1, WWTR1, LATS1, LATS2); mechanosensitive ion channels (PIEZO1, TRPV4); phosphoinositide kinases (PIP4K2B, PIP5K1C, PIK3CA); cell cycle (CDKN1A, CDKN1B, CCND1, MKI67); hypoxia (HIF1A, EPAS1, VEGFA); PIP2 / cytoskeletal (ANXA1, KPNA4, RAE1, PLS1, ZYX, ARF1, CDC42, EZR, LAMP2, HMOX1, FLNC) — plus a PIP2 composite score (mean z of ANXA1, ARF1, CDC42, EZR, VAMP3)

---

## Kaplan-Meier Curve Generation

### Molecular Alterations (Axis 1)
Genetic differences derived from **CNA (copy number alteration), SV (structural variation), and mutation** data pooled across 7 cohorts. Patients stratified as altered vs. wildtype based on presence of homozygous deletions (suppressors: CNA ≤ −2), amplifications (oncogenes: CNA ≥ +2), or deleterious mutations. Primary endpoint: overall survival (OS).

### mRNA Expression (Axes 2a, 2b)
Expression differences stratified **within each cohort independently** using z-score normalization. Patients stratified as high (z > 1.0) vs. rest using the ZSCORE threshold. Primary endpoint: DFS for TCGA-PRAD and MSKCC; OS for SU2C/PCF CRPC.

### Crowding / Mechanobiology (Axis 3)
mRNA expression in TCGA-PRAD only (DFS primary, PFS secondary; OS omitted as event-poor). Per-gene KM by median split (quartile for ANXA1 and the PIP2 composite). PTEN dependence tested via a Cox `gene_z × PTEN_DEEPDEL` interaction term (deep deletion, CNA == −2). Benjamini-Hochberg FDR within test families (q < 0.10). Run by `module8_crowding.py`.

---

## Cox Regression

Univariate and multivariate Cox proportional hazards regression performed for all genes and composite alterations. Multivariate models include cohort indicators (for pooled molecular models) and age where available (MSKCC excluded due to missing age data). Hazard ratios (HR), 95% confidence intervals (CI), and p-values reported for each predictor.

---

## Environment Setup

**Python Version:** Python 3.8 or higher

Create and activate a conda environment:
```bash
conda create -n km_prostate python=3.9
conda activate km_prostate
cd km_analysis
pip install -r requirements.txt
```

Or, recreate from an existing environment file (if available):
```bash
conda env create -f environment.yml
conda activate km_prostate
```

---

## How to Run

Full pipeline (re-downloads all data):
```bash
cd km_analysis/scripts
python run_all.py --refresh
```

Resume from a specific module (uses cached data):
```bash
python run_all.py --start module6_comprehensive_km
```

Run one module only:
```bash
python run_all.py --only module7
```

---

## Key Thresholds

| Parameter | Value |
|-----------|-------|
| CNA threshold (strict) | Suppressor ≤ −2 (homozygous del), Oncogene ≥ +2 (amplification) |
| mRNA z-score cutoff | z > 1.0 (ZSCORE split) |
| Min patients per KM arm | 20 |
| Min events for Cox | 10 |
| Significance threshold | p < 0.05 |
| Crowding axes FDR (Axis 6) | Benjamini-Hochberg, q < 0.10 |
| PTEN stratifier (crowding) | Deep deletion, CNA == −2 (`PTEN_DEEPDEL`) |

---

## Directory Structure

```
km_analysis/
├── scripts/
│   ├── config.py                    # Gene lists, cohort IDs, thresholds, paths
│   ├── module2_download.py          # Fetch clinical, CNA, mutation, SV, mRNA from cBioPortal
│   ├── module3_clean.py             # Data cleaning and standardization
│   ├── module4_flag_molecular.py    # Molecular alteration flags (CNA/mut/SV)
│   ├── module5_flag_expression.py   # Expression stratification (z-score, median, quartile)
│   ├── module6_comprehensive_km.py  # Kaplan-Meier curves
│   ├── module7_cox_regression.py    # Cox regression (univariate + multivariate)
│   ├── module8_crowding.py          # Crowding axes 1-6 (mRNA KM + PTEN Cox interaction + FDR)
│   ├── run_all.py                   # Master pipeline runner
│   └── utils/
│       ├── km_engine.py             # KaplanMeierCurve class, log-rank test, plotting
│       ├── api_client.py            # cBioPortal API wrappers
│       ├── axes.py                  # Composite score helpers
│       └── logger.py                # Logging setup
├── data/
│   └── cache/                       # Intermediate CSVs
├── outputs/
│   ├── figures/                     # KM plots (PDF/SVG)
│   └── tables/
│       ├── km_statistics.csv        # KM statistics (N, events, p-value)
│       ├── median_os_table.csv      # Median OS/DFS per arm
│       ├── cox_results.csv          # Cox regression results
│       ├── crowding_km_statistics.csv    # Crowding axes KM stats + q-values
│       └── crowding_cox_interaction.csv  # Crowding gene × PTEN_DEEPDEL interaction
└── logs/                            # Pipeline logs
```
*Generated with assistance from Claude AI coding agent.*