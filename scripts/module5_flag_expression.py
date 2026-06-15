"""
module5_flag_expression.py — Assigns High/Low groups from mRNA z-scores.

Runs for 4 cohorts only (TCGA_PRAD, SU2C, MSKCC, MCTP). Each cohort is
processed independently — z-scores are never pooled across cohorts because
each study uses a different normalization reference population.

Three split strategies per gene, computed within each cohort:
  MEDIAN    — above/below cohort median
  QUARTILE  — top 25% vs bottom 25% (middle 50% excluded)
  ZSCORE    — z > 1.0 vs rest

Patients with missing z-scores always get NaN and are excluded from analysis.

Also computes OR-combined flags (any gene dysregulated) and composite
axis scores (ADHESION_MOTILITY_SCORE, ANDROGEN_SCORE).

Output: data/cache/{KEY}_flagged_expression.csv  (one per expression cohort)
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    ANDROGEN_GENES, ADHESION_MOTILITY_GENES,
    ADHESION_EPITHELIAL_GENES, ADHESION_MESENCHYMAL_GENES,
    AR_ACTIVITY_GENES,
    CROWDING_AXES, COMPOSITE_SCORES,
    EXPRESSION_COHORT_KEYS,
    CACHE_DIR, DATASETS,
)
from utils.logger import get_logger

log = get_logger("module5_flag_expression")

# Crowding axes are mRNA expression — they get the same {GENE}_ZSCORE → {GENE}_GROUP_*
# treatment as the other expression genes. CDKN1A appears here (mRNA) and also in the
# molecular axis (CNA); the two live in distinct columns and never collide.
ALL_EXPR_GENES = {**AR_ACTIVITY_GENES, **ANDROGEN_GENES, **ADHESION_MOTILITY_GENES,
                  **CROWDING_AXES}

SPLIT_NAMES = ["MEDIAN", "QUARTILE", "ZSCORE"]


# ─────────────────────────────────────────────────────────────────────────────
# NaN-SAFE SPLITTING
# Returns a Series with values "High", "Low", or NaN.
# Patients with missing z-score always get NaN — never assigned to a group.
# For quartile split: middle 50% patients also get NaN (excluded from KM).
# ─────────────────────────────────────────────────────────────────────────────

def _split_median(series):
    med = series.median()
    result = pd.Series(np.nan, index=series.index, dtype=object)
    valid  = series.notna()
    result[valid & (series > med)]  = "High"
    result[valid & (series <= med)] = "Low"
    return result


def _split_quartile(series):
    valid = series.notna()
    q75   = series[valid].quantile(0.75)
    q25   = series[valid].quantile(0.25)
    result = pd.Series(np.nan, index=series.index, dtype=object)
    result[valid & (series >= q75)] = "High"
    result[valid & (series <= q25)] = "Low"
    # Middle 50% patients: remain NaN (excluded from KM comparison)
    return result


def _split_zscore(series):
    result = pd.Series(np.nan, index=series.index, dtype=object)
    valid  = series.notna()
    result[valid & (series > 1.0)]  = "High"
    result[valid & (series <= 1.0)] = "Low"
    return result


_SPLITTERS = {
    "MEDIAN":   _split_median,
    "QUARTILE": _split_quartile,
    "ZSCORE":   _split_zscore,
}


def add_expression_splits(df):
    """
    Adds {GENE}_GROUP_{SPLIT} columns for all expression genes.
    Values: "High", "Low", or NaN.
    """
    log.info("Adding NaN-safe expression group columns...")
    for gene in ALL_EXPR_GENES:
        zscore_col = f"{gene}_ZSCORE"
        if zscore_col not in df.columns:
            log.warning(f"  {zscore_col} not found — skipping {gene}")
            continue
        series = pd.to_numeric(df[zscore_col], errors="coerce")
        n_valid = series.notna().sum()
        log.info(f"  {gene}: {n_valid} valid z-scores")
        for split_name, splitter in _SPLITTERS.items():
            col = f"{gene}_GROUP_{split_name}"
            df[col] = splitter(series)
            n_high = (df[col] == "High").sum()
            n_low  = (df[col] == "Low").sum()
            n_nan  = df[col].isna().sum()
            log.info(f"    {split_name}: High={n_high}, Low={n_low}, excluded={n_nan}")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# OR-COMBINED FLAGS (NaN-safe)
# ─────────────────────────────────────────────────────────────────────────────

def _or_combined(df, genes, split_name, bad_value):
    """
    OR-combined flag across a set of genes for one split strategy.
    bad_value: "High" for all-high-is-bad genes, "Low" for all-low-is-bad genes.

    Returns Series: 1 = any gene dysregulated, 0 = all genes have data & none dysregulated, NaN = excluded.
    Wildtype (0) requires ALL genes to have non-NaN GROUP values AND none meets the criterion.
    """
    cols = [f"{g}_GROUP_{split_name}" for g in genes if f"{g}_GROUP_{split_name}" in df.columns]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    flags     = df[cols]
    any_bad   = (flags == bad_value).any(axis=1)
    all_data  = flags.notna().all(axis=1)
    result    = np.where(any_bad, 1, np.where(all_data, 0, np.nan))
    return pd.Series(result, index=df.index)


def add_or_combined_expression(df):
    """
    Adds OR-combined flags for androgen and adhesion/motility axes per split.
    Androgen OR: any of SLCO2B1/SLCO1B3/AKR1C3 HIGH (all high_is_bad)
    Adhesion OR: any epithelial gene LOW OR any mesenchymal gene HIGH
    """
    androgen_genes = list(ANDROGEN_GENES.keys())
    epith_genes    = ADHESION_EPITHELIAL_GENES
    mesen_genes    = ADHESION_MESENCHYMAL_GENES

    for split_name in SPLIT_NAMES:
        # Androgen OR — any high = dysregulated
        df[f"ANDROGEN_OR_{split_name}"] = _or_combined(
            df, androgen_genes, split_name, "High")

        # Adhesion OR — directional: epith LOW or mesen HIGH
        epith_cols = [f"{g}_GROUP_{split_name}" for g in epith_genes
                      if f"{g}_GROUP_{split_name}" in df.columns]
        mesen_cols = [f"{g}_GROUP_{split_name}" for g in mesen_genes
                      if f"{g}_GROUP_{split_name}" in df.columns]
        all_cols   = epith_cols + mesen_cols

        if all_cols:
            flags     = df[all_cols]
            any_epith_low  = (df[epith_cols] == "Low").any(axis=1)  if epith_cols else pd.Series(False, index=df.index)
            any_mesen_high = (df[mesen_cols] == "High").any(axis=1) if mesen_cols else pd.Series(False, index=df.index)
            any_bad   = any_epith_low | any_mesen_high
            all_data  = flags.notna().all(axis=1)
            result    = np.where(any_bad, 1, np.where(all_data, 0, np.nan))
            df[f"ADHESION_OR_{split_name}"] = pd.Series(result, index=df.index)

        n_and = (df[f"ANDROGEN_OR_{split_name}"] == 1).sum()
        n_adh = (df.get(f"ADHESION_OR_{split_name}", pd.Series()) == 1).sum()
        log.info(f"  OR {split_name}: androgen dysreg={n_and}, adhesion dysreg={n_adh}")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# COMPOSITE SCORE
# ─────────────────────────────────────────────────────────────────────────────

def add_composite_scores(df):
    """
    Adhesion-Motility composite score: mesenchymal_mean - epithelial_mean.
    """
    epith_cols = [f"{g}_ZSCORE" for g in ADHESION_EPITHELIAL_GENES
                  if f"{g}_ZSCORE" in df.columns]
    mesen_cols = [f"{g}_ZSCORE" for g in ADHESION_MESENCHYMAL_GENES
                  if f"{g}_ZSCORE" in df.columns]

    epith_mean = df[epith_cols].mean(axis=1) if epith_cols else pd.Series(np.nan, index=df.index)
    mesen_mean = df[mesen_cols].mean(axis=1) if mesen_cols else pd.Series(np.nan, index=df.index)
    df["ADHESION_MOTILITY_SCORE"] = mesen_mean - epith_mean

    log.info(f"  ADHESION_MOTILITY_SCORE: "
             f"valid={df['ADHESION_MOTILITY_SCORE'].notna().sum()}, "
             f"mean={df['ADHESION_MOTILITY_SCORE'].mean():.3f}")

    androgen_cols = [f"{g}_ZSCORE" for g in ANDROGEN_GENES
                     if f"{g}_ZSCORE" in df.columns]
    df["ANDROGEN_SCORE"] = (df[androgen_cols].mean(axis=1)
                            if androgen_cols else pd.Series(np.nan, index=df.index))
    log.info(f"  ANDROGEN_SCORE: "
             f"valid={df['ANDROGEN_SCORE'].notna().sum()}, "
             f"mean={df['ANDROGEN_SCORE'].mean():.3f}")

    # PIP2 trafficking composite (Axis 6 headline): mean z of ANXA1+ARF1+CDC42+EZR+VAMP3.
    # Dichotomized top vs bottom quartile in module8 (middle 50% → NaN). Equal weights.
    pip2_cfg  = COMPOSITE_SCORES["pip2_trafficking_score"]
    pip2_cols = [f"{g}_ZSCORE" for g in pip2_cfg["genes"] if f"{g}_ZSCORE" in df.columns]
    if pip2_cols:
        df["PIP2_TRAFFICKING_SCORE"] = df[pip2_cols].mean(axis=1)
        # Quartile group column (top=High, bottom=Low, middle 50%=NaN)
        df["PIP2_TRAFFICKING_GROUP_QUARTILE"] = _split_quartile(df["PIP2_TRAFFICKING_SCORE"])
        log.info(f"  PIP2_TRAFFICKING_SCORE: valid={df['PIP2_TRAFFICKING_SCORE'].notna().sum()} "
                 f"(from {len(pip2_cols)}/{len(pip2_cfg['genes'])} genes), "
                 f"mean={df['PIP2_TRAFFICKING_SCORE'].mean():.3f}")
    else:
        df["PIP2_TRAFFICKING_SCORE"] = np.nan
        df["PIP2_TRAFFICKING_GROUP_QUARTILE"] = np.nan
        log.warning("  PIP2_TRAFFICKING_SCORE: no component z-scores found — set NaN")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def flag_all(molecular_dfs=None, force_refresh=False):
    """Run expression flagging for mRNA-eligible cohorts only.

    Expression splits (median, quartile, z-score) are computed WITHIN each cohort
    using that cohort's own distribution. The COMBINED dataset is excluded because
    pooling z-scores across studies with different normalization references is invalid.

    Eligible cohorts: EXPRESSION_COHORT_KEYS = TCGA_PRAD, SU2C, MSKCC, MCTP
    """
    results = {}

    # Explicitly log which cohorts are being skipped and why
    all_keys = list(DATASETS.keys())
    skipped  = [k for k in all_keys if k not in EXPRESSION_COHORT_KEYS]
    log.info(f"Expression flagging — eligible cohorts: {EXPRESSION_COHORT_KEYS}")
    log.info(f"Skipped (no mRNA or not per-cohort eligible): {skipped}")
    log.info("NOTE: COMBINED dataset intentionally excluded — cross-cohort z-score "
             "pooling is invalid (different normalization baselines).")

    for key in EXPRESSION_COHORT_KEYS:
        out_path = CACHE_DIR / f"{key}_flagged_expression.csv"

        if out_path.exists() and not force_refresh:
            log.info(f"Loading expression flags from cache: {out_path.name}")
            results[key] = pd.read_csv(out_path)
            continue

        if molecular_dfs and key in molecular_dfs:
            df = molecular_dfs[key].copy()
        else:
            mol_path = CACHE_DIR / f"{key}_flagged_molecular.csv"
            if not mol_path.exists():
                raise FileNotFoundError(
                    f"Missing molecular flagged file: {mol_path}. "
                    f"Run module4_flag_molecular.py first.")
            df = pd.read_csv(mol_path)

        label = DATASETS[key]["label"]
        n_valid_mrna = df[[c for c in df.columns if '_ZSCORE' in c]].notna().any(axis=1).sum()
        log.info(f"\n{'='*55}\nExpression flagging: {label}  "
                 f"(n={len(df)}, with mRNA={n_valid_mrna})\n{'='*55}")

        # PTEN_DEEPDEL (crowding stratifier) flows in from the molecular flagged frame.
        # Warn loudly if absent — module4 must be re-run after adding it.
        if "PTEN_DEEPDEL" not in df.columns:
            log.warning("  PTEN_DEEPDEL missing — re-run module4_flag_molecular --refresh "
                        "for crowding (module8) PTEN-stratified analysis.")

        df = add_expression_splits(df)
        df = add_or_combined_expression(df)
        df = add_composite_scores(df)

        df.to_csv(out_path, index=False)
        log.info(f"Saved: {out_path}")
        results[key] = df

    return results


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--refresh", action="store_true")
    args = p.parse_args()
    flag_all(force_refresh=args.refresh)
    log.info("Expression flagging complete.")
