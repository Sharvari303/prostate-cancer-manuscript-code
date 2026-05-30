"""
utils/axes.py
─────────────────────────────────────────────────────────────────────────────
Shared axis computation functions for KM and Cox regression analysis.

Gene panels (sourced from config.py — do not hardcode here):
  Epithelial (low_is_bad):  CDH1, EPCAM
  Mesenchymal (high_is_bad): VIM, CDH2, MMP9, CXCL8, SNAI1, RHOA, ITGB1, FN1
  Androgen uptake:           SLCO2B1, SLCO1B3, AKR1C3

Composite score formula:
  adhesion_motility_score = mean(mesenchymal genes) - mean(epithelial genes)
  High score = aggressive invasive/motile phenotype.

Previously included genes now removed:
  ACTB — housekeeping gene; z-score uniformly near 0; not informative
  KRT8, KRT18, MYH9, ACTA2, TUBA1B, TUBB3 — never downloaded; silently returned 0;
    removed to eliminate silent silent fallback behaviour

Missing data policy (STRICT):
  - If a patient has NO data for any gene in an axis → NaN (excluded from analysis)
  - If a patient has data for SOME genes in an OR analysis → flag from available genes
  - NaN patients are dropped at the plot/model level via .dropna()
  - Patients are NEVER silently assigned to wild-type / low group
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logger import get_logger

log = get_logger("axes")


# ─────────────────────────────────────────────────────────────────────────────
# MOLECULAR AXIS
# ─────────────────────────────────────────────────────────────────────────────

def get_molecular_genes(df):
    """
    Extract individual molecular alterations (NaN-safe).

    Each gene flag is:
      - 1  if the patient has a confirmed alteration
      - 0  if the patient has confirmed no alteration (data present, not altered)
      - NaN if no data source is available for this patient for this gene

    Requires: PTEN_CNA, PTEN_MUT, MDM2_CNA, MDM2_MUT, MDM4_CNA, AR_CNA, AR_MUT
    """
    genes = {}

    # PTEN loss: CNA <= -2 OR truncating mutation
    if 'PTEN_CNA' in df.columns or 'PTEN_MUT' in df.columns:
        cna_avail = 'PTEN_CNA' in df.columns
        mut_avail = 'PTEN_MUT' in df.columns
        cna_flag = (df['PTEN_CNA'] <= -2) if cna_avail else pd.Series(False, index=df.index)
        mut_flag = (df['PTEN_MUT'] == 1) if mut_avail else pd.Series(False, index=df.index)
        altered = cna_flag | mut_flag
        # Restore NaN where no data existed at all
        any_data = (df['PTEN_CNA'].notna() if cna_avail else pd.Series(False, index=df.index)) | \
                   (df['PTEN_MUT'].notna() if mut_avail else pd.Series(False, index=df.index))
        altered = altered.where(any_data, other=np.nan)
        genes['PTEN_loss'] = altered

    # MDM2 gain: CNA >= +2 OR missense mutation
    if 'MDM2_CNA' in df.columns or 'MDM2_MUT' in df.columns:
        cna_avail = 'MDM2_CNA' in df.columns
        mut_avail = 'MDM2_MUT' in df.columns
        cna_flag = (df['MDM2_CNA'] >= +2) if cna_avail else pd.Series(False, index=df.index)
        mut_flag = (df['MDM2_MUT'] == 1) if mut_avail else pd.Series(False, index=df.index)
        altered = cna_flag | mut_flag
        any_data = (df['MDM2_CNA'].notna() if cna_avail else pd.Series(False, index=df.index)) | \
                   (df['MDM2_MUT'].notna() if mut_avail else pd.Series(False, index=df.index))
        altered = altered.where(any_data, other=np.nan)
        genes['MDM2_gain'] = altered

    # MDM4 gain: CNA >= +2
    if 'MDM4_CNA' in df.columns:
        altered = (df['MDM4_CNA'] >= +2)
        any_data = df['MDM4_CNA'].notna()
        altered = altered.where(any_data, other=np.nan)
        genes['MDM4_gain'] = altered

    # AR gain: CNA >= +2 OR gain-of-function mutation (mutation-primary)
    if 'AR_CNA' in df.columns or 'AR_MUT' in df.columns:
        cna_avail = 'AR_CNA' in df.columns
        mut_avail = 'AR_MUT' in df.columns
        cna_flag = (df['AR_CNA'] >= +2) if cna_avail else pd.Series(False, index=df.index)
        mut_flag = (df['AR_MUT'] == 1) if mut_avail else pd.Series(False, index=df.index)
        altered = cna_flag | mut_flag
        any_data = (df['AR_CNA'].notna() if cna_avail else pd.Series(False, index=df.index)) | \
                   (df['AR_MUT'].notna() if mut_avail else pd.Series(False, index=df.index))
        altered = altered.where(any_data, other=np.nan)
        genes['AR_gain'] = altered

    return genes


def get_individual_genes(df):
    """
    Extract individual supplementary gene alterations (CDKN1A, BAX, ATM) — NaN-safe.

    Loss-of-function: CNA <= -2 OR truncating mutation.
    Returns dict of gene_name → Series (1/0/NaN).
    """
    genes = {}
    gene_configs = {
        'CDKN1A': ('CDKN1A_CNA', 'CDKN1A_MUT', 'loss'),
        'BAX':    ('BAX_CNA',    'BAX_MUT',    'loss'),
        'ATM':    ('ATM_CNA',    'ATM_MUT',    'loss'),
    }
    for gene_name, (cna_col, mut_col, direction) in gene_configs.items():
        cna_avail = cna_col in df.columns
        mut_avail = mut_col in df.columns
        if not cna_avail and not mut_avail:
            continue
        cna_flag = (df[cna_col] <= -2) if cna_avail else pd.Series(False, index=df.index)
        mut_flag = (df[mut_col] == 1) if mut_avail else pd.Series(False, index=df.index)
        altered = cna_flag | mut_flag
        any_data = (df[cna_col].notna() if cna_avail else pd.Series(False, index=df.index)) | \
                   (df[mut_col].notna() if mut_avail else pd.Series(False, index=df.index))
        altered = altered.where(any_data, other=np.nan)
        genes[gene_name] = altered
    return genes


def get_molecular_burden(df):
    """Count alterations: 0, 1, 2, or 3+ genes altered (PTEN/MDM2/MDM4/AR).

    Returns NaN for patients where all genes are NaN (no data).
    """
    mol_genes = get_molecular_genes(df)
    if not mol_genes:
        return pd.Series(np.nan, index=df.index)
    # Sum with skipna=False so that if all-NaN → NaN; partial data sums available flags
    flags = pd.DataFrame(mol_genes)
    # Convert True/False to 1/0 for summation, keep NaN as NaN
    flags = flags.apply(pd.to_numeric, errors='coerce')
    # Patient-level sum: count how many genes are altered (True=1)
    # skipna=True here: if patient has some NaN genes, sum the ones we know about
    # A patient with ALL NaN gets NaN
    has_any_data = flags.notna().any(axis=1)
    burden = flags.sum(axis=1, skipna=True)
    burden = burden.where(has_any_data, other=np.nan)
    return burden


def get_molecular_binary(df):
    """Binary: any alteration vs none (PTEN/MDM2/MDM4/AR) — NaN-safe.

    Returns:
      1   = at least one gene altered
      0   = all genes with data are unaltered
      NaN = no molecular data available for this patient
    """
    mol_genes = get_molecular_genes(df)
    if not mol_genes:
        return pd.Series(np.nan, index=df.index)
    flags = pd.DataFrame(mol_genes).apply(pd.to_numeric, errors='coerce')
    has_any_data = flags.notna().any(axis=1)
    any_altered = flags.any(axis=1, skipna=True)  # True if any gene is 1
    result = any_altered.astype(float)
    result = result.where(has_any_data, other=np.nan)
    return result


def get_all_molecular_or_flag(df):
    """
    OR-combined: any of PTEN/MDM2/MDM4/AR/CDKN1A/BAX/ATM altered — NaN-safe.

    Includes both molecular axis genes and individual genes.
    Patient excluded only if they have NO data for ANY of the 7 genes.
    """
    mol_genes = get_molecular_genes(df)
    indiv_genes = get_individual_genes(df)
    all_flags = {**mol_genes, **indiv_genes}
    if not all_flags:
        return pd.Series(np.nan, index=df.index)
    flags = pd.DataFrame(all_flags).apply(pd.to_numeric, errors='coerce')
    has_any_data = flags.notna().any(axis=1)
    any_altered = flags.any(axis=1, skipna=True)
    result = any_altered.astype(float)
    result = result.where(has_any_data, other=np.nan)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# ANDROGEN AXIS
# ─────────────────────────────────────────────────────────────────────────────

def get_androgen_score(df):
    """Androgen uptake score: mean(SLCO2B1, SLCO1B3, AKR1C3).

    Returns NaN if no androgen z-score columns available for a patient.
    """
    cols = ['SLCO2B1_ZSCORE', 'SLCO1B3_ZSCORE', 'AKR1C3_ZSCORE']
    available = [c for c in cols if c in df.columns]
    if not available:
        return pd.Series(np.nan, index=df.index)
    score = df[available].mean(axis=1)  # NaN if all values NaN for that row
    return score


def get_androgen_binary(df):
    """High vs Low androgen uptake (composite score > 1.0, ~top 16% expression).

    Returns NaN if patient has no androgen z-score data.
    """
    score = get_androgen_score(df)
    result = (score > 1.0).astype(float)
    result = result.where(score.notna(), other=np.nan)
    return result


def get_androgen_or_flag(df):
    """
    OR-combined androgen flag: any of SLCO2B1/SLCO1B3/AKR1C3 has z > 1.0 — NaN-safe.

    Returns:
      1   = at least one gene upregulated (z > 1.0)
      0   = all available genes at z <= 1.0
      NaN = no androgen gene data available for this patient
    """
    cols = ['SLCO2B1_ZSCORE', 'SLCO1B3_ZSCORE', 'AKR1C3_ZSCORE']
    available = [c for c in cols if c in df.columns]
    if not available:
        return pd.Series(np.nan, index=df.index)
    flags = pd.DataFrame({c: (df[c] > 1.0) for c in available})
    has_any_data = df[available].notna().any(axis=1)
    any_high = flags.any(axis=1, skipna=True)
    result = any_high.astype(float)
    result = result.where(has_any_data, other=np.nan)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# ADHESION & MOTILITY AXIS
# ─────────────────────────────────────────────────────────────────────────────
# Composite score:
#   epithelial_retention = mean(CDH1, EPCAM)                          — SUBTRACTED
#   mesenchymal          = mean(VIM, CDH2, MMP9, CXCL8, SNAI1,        — ADDED
#                               RHOA, ITGB1, FN1)
#   score = mesenchymal_mean - epithelial_mean
#
# High score = aggressive invasive/motile phenotype.
#
# Gene panel matches ADHESION_EPITHELIAL_GENES / ADHESION_MESENCHYMAL_GENES in config.py.
# Removed: ACTB (housekeeping), KRT8/KRT18/MYH9/ACTA2/TUBA1B/TUBB3 (never downloaded).
# ─────────────────────────────────────────────────────────────────────────────

_EPITHELIAL_COLS = ['CDH1_ZSCORE', 'EPCAM_ZSCORE']
_MESENCHYMAL_COLS = [
    'VIM_ZSCORE', 'CDH2_ZSCORE', 'MMP9_ZSCORE', 'CXCL8_ZSCORE',
    'SNAI1_ZSCORE', 'RHOA_ZSCORE', 'ITGB1_ZSCORE', 'FN1_ZSCORE',
]

# Direction of harm for each gene (for OR-combined flag)
# "high_is_bad": z > 1.0 → bad direction
# "low_is_bad":  z < -1.0 → bad direction (epithelial gene loss)
_ADHESION_GENE_DIRECTIONS = {
    'CDH1_ZSCORE':   'low_is_bad',
    'EPCAM_ZSCORE':  'low_is_bad',
    'VIM_ZSCORE':    'high_is_bad',
    'CDH2_ZSCORE':   'high_is_bad',
    'MMP9_ZSCORE':   'high_is_bad',
    'CXCL8_ZSCORE':  'high_is_bad',
    'SNAI1_ZSCORE':  'high_is_bad',
    'RHOA_ZSCORE':   'high_is_bad',
    'ITGB1_ZSCORE':  'high_is_bad',
    'FN1_ZSCORE':    'high_is_bad',
}


def get_adhesion_motility_score(df):
    """Adhesion & Motility composite score.

    Formula: mesenchymal_mean - epithelial_mean
    where each sub-score is the mean of available genes in that group.
    Returns NaN if no mesenchymal columns are available (epithelial defaults to 0).
    """
    epith_avail = [c for c in _EPITHELIAL_COLS if c in df.columns]
    mesen_avail = [c for c in _MESENCHYMAL_COLS if c in df.columns]

    if not mesen_avail:
        log.warning("No mesenchymal columns found for adhesion_motility_score")
        return pd.Series(np.nan, index=df.index)

    epith_score = df[epith_avail].mean(axis=1) if epith_avail else pd.Series(0.0, index=df.index)
    mesen_score = df[mesen_avail].mean(axis=1)

    return mesen_score - epith_score


def get_adhesion_motility_binary(df):
    """High vs Low adhesion-motility (composite score > 1.0, ~top 16% invasive phenotype).

    Returns NaN if patient has no adhesion/motility data.
    """
    score = get_adhesion_motility_score(df)
    result = (score > 1.0).astype(float)
    result = result.where(score.notna(), other=np.nan)
    return result


def get_adhesion_motility_or_flag(df):
    """
    OR-combined adhesion/motility flag — NaN-safe.

    A patient is flagged (1) if ANY gene shows expression in its harmful direction:
      - Epithelial genes (CDH1/KRT8/KRT18/EPCAM): z < -1.0  (loss = bad)
      - Mesenchymal/motility genes: z > 1.0  (gain = bad)

    Returns:
      1   = at least one gene dysregulated in harmful direction
      0   = no gene dysregulated (all available genes in normal range)
      NaN = no adhesion/motility gene data available for this patient
    """
    all_cols = list(_ADHESION_GENE_DIRECTIONS.keys())
    available = [c for c in all_cols if c in df.columns]
    if not available:
        return pd.Series(np.nan, index=df.index)

    flags = pd.DataFrame(index=df.index)
    for col in available:
        direction = _ADHESION_GENE_DIRECTIONS[col]
        if direction == 'low_is_bad':
            flags[col] = (df[col] < -1.0)
        else:
            flags[col] = (df[col] > 1.0)

    has_any_data = df[available].notna().any(axis=1)
    any_dysregulated = flags.any(axis=1, skipna=True)
    result = any_dysregulated.astype(float)
    result = result.where(has_any_data, other=np.nan)
    return result

