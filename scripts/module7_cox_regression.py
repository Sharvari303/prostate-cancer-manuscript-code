"""
module7_cox_regression.py — Univariate + Multivariate Cox regression for all tasks.

Min events per model: 10.

Molecular pool analyses:
  7 cohorts combined, OS endpoint
  Covariates: COHORT_CODE (integer) + AGE (where available)

Per-cohort mRNA analyses:
  4 mRNA cohorts, cohort-specific endpoint (OS or DFS per config.ENDPOINTS)
  Covariates: AGE (where available); no cohort term (single cohort)
  Runs both continuous z-score Cox and binary-per-split Cox.

OR-combined flags are binary (0/1/NaN).

Output: outputs/tables/cox_results.csv
"""
import sys
import warnings
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    DATASETS, MOLECULAR_COHORT_KEYS, EXPRESSION_COHORT_KEYS,
    MOLECULAR_GENES, INDIVIDUAL_GENES,
    ANDROGEN_GENES, ANDROGEN_CNA_GENES,
    ADHESION_MOTILITY_GENES, ADHESION_CNA_GENES,
    ENDPOINTS,
    CACHE_DIR, TABLES_DIR, STATS,
    DEFAULT_CNA_THRESHOLD,
)
from utils.logger import get_logger

log = get_logger("module7_cox_regression")

MIN_EVENTS = STATS["min_cox_events"]   # 10

# Default endpoint for molecular pool (all 7 cohorts use OS)
TIME_COL   = "OS_MONTHS"
EVENT_COL  = "OS_EVENT"


def _get_endpoint(cohort_key):
    """Return (time_col, event_col, endpoint_short) for a cohort."""
    ep = ENDPOINTS.get(cohort_key, ENDPOINTS["DEFAULT"])["primary"]
    return ep["time"], ep["event"], ep["short"]

SPLIT_NAMES = ["MEDIAN", "QUARTILE", "ZSCORE"]


# ─────────────────────────────────────────────────────────────────────────────
# UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

def _to_binary(series):
    """Convert True/False/NaN (string or bool from CSV) to 1.0/0.0/NaN float."""
    def _map(v):
        if pd.isna(v):
            return np.nan
        if v is True or str(v).strip().lower() == "true":
            return 1.0
        if v is False or str(v).strip().lower() == "false":
            return 0.0
        try:
            return float(v)
        except (ValueError, TypeError):
            return np.nan
    return series.map(_map)


def _clean_endpoint(df, time_col=TIME_COL, event_col=EVENT_COL):
    """Coerce and filter survival columns for any endpoint (OS or DFS)."""
    df = df.copy()
    df[time_col]  = pd.to_numeric(df.get(time_col,  pd.Series(dtype=float)), errors="coerce")
    df[event_col] = pd.to_numeric(df.get(event_col, pd.Series(dtype=float)), errors="coerce")
    df = df.dropna(subset=[time_col, event_col])
    df = df[df[time_col] > 0].copy()
    return df


def _clean_os(df):
    """Legacy alias — uses OS endpoint."""
    return _clean_endpoint(df, TIME_COL, EVENT_COL)


def _group_to_cox_binary(series, direction):
    """
    Convert GROUP string column ("High"/"Low"/NaN) to Cox-ready binary.
    Bad group is always encoded as 1 so HR > 1 means worse prognosis.
    high_is_bad: High=1, Low=0
    low_is_bad:  Low=1,  High=0
    """
    result = pd.Series(np.nan, index=series.index, dtype=float)
    if direction == "high_is_bad":
        result[series == "High"] = 1.0
        result[series == "Low"]  = 0.0
    else:
        result[series == "Low"]  = 1.0
        result[series == "High"] = 0.0
    return result


# ─────────────────────────────────────────────────────────────────────────────
# COX ENGINE
# ─────────────────────────────────────────────────────────────────────────────

def _run_cox(df, feature_col, covariates=None, time_col=None, event_col=None):
    """
    Fit univariate or multivariate Cox model.
    Returns dict with HR, CI_lower, CI_upper, pval, N, events — or None.
    Requires at least MIN_EVENTS events in the clean subset.
    """
    try:
        from lifelines import CoxPHFitter
    except ImportError:
        log.error("lifelines not installed — cannot run Cox regression")
        return None

    tc = time_col  or TIME_COL
    ec = event_col or EVENT_COL

    cols = [tc, ec, feature_col] + (covariates or [])
    df_c = df[cols].copy()
    df_c = df_c.dropna()
    df_c = df_c[df_c[tc] > 0]
    df_c[ec] = pd.to_numeric(df_c[ec], errors="coerce")
    df_c = df_c.dropna(subset=[ec])

    n_events = int((df_c[ec] == 1).sum())
    if n_events < MIN_EVENTS or len(df_c) < 20:
        return None

    # Drop constant columns (causes singular matrix)
    for col in [feature_col] + (covariates or []):
        if df_c[col].nunique() < 2:
            return None

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cph = CoxPHFitter()
            cph.fit(df_c, duration_col=tc, event_col=ec,
                    show_progress=False)

        hr      = float(np.exp(cph.params_[feature_col]))
        ci_lo   = float(np.exp(cph.confidence_intervals_.loc[feature_col, "95% lower-bound"]))
        ci_hi   = float(np.exp(cph.confidence_intervals_.loc[feature_col, "95% upper-bound"]))
        pval    = float(cph.summary.loc[feature_col, "p"])
        return {
            "HR": hr, "CI_lower": ci_lo, "CI_upper": ci_hi,
            "pval": pval, "N": len(df_c), "events": n_events,
        }
    except Exception as e:
        log.debug(f"    Cox failed: {e}")
        return None


def _record(rows, task, gene, cohort, split, predictor_type,
            label_bad, covariates, uni, multi, endpoint="OS"):
    """Append one result row — handles None uni/multi gracefully."""
    def _v(d, k): return d[k] if d else np.nan
    covars_str = ",".join(covariates) if covariates else "none"
    rows.append({
        "task":          task,
        "gene":          gene,
        "cohort":        cohort,
        "split":         split,
        "endpoint":      endpoint,
        "predictor":     predictor_type,
        "label_bad":     label_bad,
        "covariates":    covars_str,
        "N":             _v(uni, "N"),
        "events":        _v(uni, "events"),
        "HR_uni":        _v(uni, "HR"),
        "CI_uni_lo":     _v(uni, "CI_lower"),
        "CI_uni_hi":     _v(uni, "CI_upper"),
        "p_uni":         _v(uni, "pval"),
        "HR_multi":      _v(multi, "HR"),
        "CI_multi_lo":   _v(multi, "CI_lower"),
        "CI_multi_hi":   _v(multi, "CI_upper"),
        "p_multi":       _v(multi, "pval"),
    })
    if uni:
        log.info(f"    HR={uni['HR']:.2f} [{uni['CI_lower']:.2f}–{uni['CI_upper']:.2f}] "
                 f"p={uni['pval']:.3f}  n={uni['N']} ev={uni['events']}")


def _cox_pair(df, feat_col, label_bad, task, gene, cohort, split,
              predictor_type, covariates, rows,
              time_col=None, event_col=None, endpoint="OS"):
    """Run univariate + multivariate Cox for one predictor; record result."""
    uni   = _run_cox(df, feat_col, covariates=None,
                     time_col=time_col, event_col=event_col)
    multi = (_run_cox(df, feat_col, covariates=covariates,
                      time_col=time_col, event_col=event_col)
             if covariates else None)
    _record(rows, task, gene, cohort, split, predictor_type,
            label_bad, covariates, uni, multi, endpoint=endpoint)


# ─────────────────────────────────────────────────────────────────────────────
# DATA LOADING
# ─────────────────────────────────────────────────────────────────────────────

def _load_molecular_pool():
    """Stack flagged molecular CSVs for the 7-cohort molecular OS pool.

    Uses MOLECULAR_COHORT_KEYS (excludes MSKCC and IDH_MUTANT).
    All cohorts use OS endpoint.
    """
    dfs = []
    for i, key in enumerate(MOLECULAR_COHORT_KEYS):
        path = CACHE_DIR / f"{key}_flagged_molecular.csv"
        if not path.exists():
            log.warning(f"  Missing: {path.name}")
            continue
        df = pd.read_csv(path, low_memory=False)
        df["COHORT_KEY"]  = key
        df["COHORT_CODE"] = i
        dfs.append(df)
    if not dfs:
        raise FileNotFoundError("No flagged molecular files found for 7-cohort set")
    combined = pd.concat(dfs, ignore_index=True)
    combined = _clean_os(combined)
    combined["AGE"] = pd.to_numeric(combined.get("AGE", pd.Series(dtype=float)), errors="coerce")
    log.info(f"  7-cohort: {len(combined)} patients, AGE available={combined['AGE'].notna().sum()}")
    return combined


def _molecular_pool_covariates(df):
    """Build covariate list for 7-cohort combined model."""
    covars = ["COHORT_CODE"]
    if df["AGE"].notna().sum() >= 50:
        covars.append("AGE")
    return covars


def _load_mrna_cohort(key):
    """Load flagged expression CSV for one mRNA cohort using its primary endpoint.

    Returns (df, time_col, event_col, ep_label) or (None, ...) if missing.
    TCGA_PRAD and MSKCC use DFS; SU2C and MCTP use OS.
    """
    path = CACHE_DIR / f"{key}_flagged_expression.csv"
    if not path.exists():
        return None, None, None, None
    time_col, event_col, ep_label = _get_endpoint(key)
    df = pd.read_csv(path, low_memory=False)
    df = _clean_endpoint(df, time_col, event_col)
    df["AGE"] = pd.to_numeric(df.get("AGE", pd.Series(dtype=float)), errors="coerce")
    log.info(f"  {key}: {len(df)} patients with valid {ep_label} "
             f"({time_col}/{event_col})")
    return df, time_col, event_col, ep_label


def _mrna_cohort_covariates(df):
    """Build covariate list for single-cohort mRNA model (AGE only if available)."""
    if df["AGE"].notna().sum() >= 20:
        return ["AGE"]
    return []


# ─────────────────────────────────────────────────────────────────────────────
# TASK 1 — Individual molecular genes, 7 cohorts combined
# ─────────────────────────────────────────────────────────────────────────────

def cox_task1(mol_pool, covars, rows):
    log.info("\n--- Cox Task 1: Individual molecular genes ---")
    all_genes = list(MOLECULAR_GENES.keys()) + list(INDIVIDUAL_GENES.keys())
    for gene in all_genes:
        col = f"{gene}_ALT_{DEFAULT_CNA_THRESHOLD.upper()}"
        if col not in mol_pool.columns:
            continue
        mol_pool[f"_feat_{gene}"] = _to_binary(mol_pool[col])
        log.info(f"  {gene}")
        _cox_pair(mol_pool, f"_feat_{gene}", f"{gene} Altered",
                  "task1", gene, "7_cohorts", "", "binary_alteration",
                  covars, rows)


# ─────────────────────────────────────────────────────────────────────────────
# TASK 2 — OR-combined (PTEN+AR+MDM4), 7 cohorts combined
# ─────────────────────────────────────────────────────────────────────────────

def cox_task2(mol_pool, covars, rows):
    log.info("\n--- Cox Task 2: OR-combined PTEN+AR+MDM4 ---")
    col = "MOLECULAR_OR"
    if col not in mol_pool.columns:
        log.warning("  MOLECULAR_OR not found")
        return
    mol_pool["_feat_mol_or"] = pd.to_numeric(mol_pool[col], errors="coerce")
    _cox_pair(mol_pool, "_feat_mol_or", "Any Altered (PTEN/AR/MDM4)",
              "task2", "PTEN_AR_MDM4", "7_cohorts", "", "binary_or",
              covars, rows)


# ─────────────────────────────────────────────────────────────────────────────
# TASK 3 — Individual androgen mRNA genes, per cohort
# Continuous z-score + binary per split
# ─────────────────────────────────────────────────────────────────────────────

def cox_task3(rows):
    log.info("\n--- Cox Task 3: Individual androgen mRNA per cohort ---")
    for cohort_key in EXPRESSION_COHORT_KEYS:
        df, time_col, event_col, ep_label = _load_mrna_cohort(cohort_key)
        if df is None:
            continue
        covars = _mrna_cohort_covariates(df)
        cohort_label = DATASETS[cohort_key]["label"]
        log.info(f"  Cohort: {cohort_label}")

        for gene in ANDROGEN_GENES:
            # Continuous z-score
            zcol = f"{gene}_ZSCORE"
            if zcol in df.columns:
                df["_feat_z"] = pd.to_numeric(df[zcol], errors="coerce")
                _cox_pair(df, "_feat_z", f"High {gene} (continuous z)",
                          "task3", gene, cohort_key, "ZSCORE_CONTINUOUS",
                          "continuous_zscore", covars, rows,
                          time_col=time_col, event_col=event_col,
                          endpoint=ep_label)

            # Binary per split
            for split in SPLIT_NAMES:
                gcol = f"{gene}_GROUP_{split}"
                if gcol not in df.columns:
                    continue
                df["_feat_b"] = _group_to_cox_binary(df[gcol], "high_is_bad")
                log.info(f"    {gene} {split}")
                _cox_pair(df, "_feat_b", f"High {gene}",
                          "task3", gene, cohort_key, split,
                          "binary_group", covars, rows,
                          time_col=time_col, event_col=event_col,
                          endpoint=ep_label)


# ─────────────────────────────────────────────────────────────────────────────
# TASK PFS TCGA — All androgen + adhesion mRNA genes vs PFS in TCGA-PRAD
# ─────────────────────────────────────────────────────────────────────────────

def cox_task_pfs_tcga(rows):
    log.info("\n--- Cox Task PFS TCGA: All androgen + adhesion mRNA vs PFS — TCGA-PRAD ---")
    path = CACHE_DIR / "TCGA_PRAD_flagged_expression.csv"
    if not path.exists():
        log.warning("  TCGA_PRAD_flagged_expression.csv missing — skipping")
        return
    df = pd.read_csv(path, low_memory=False)
    df = _clean_endpoint(df, "PFS_MONTHS", "PFS_EVENT")
    df["AGE"] = pd.to_numeric(df.get("AGE", pd.Series(dtype=float)), errors="coerce")
    log.info(f"  TCGA_PRAD PFS: {len(df)} patients")
    covars = _mrna_cohort_covariates(df)

    # Androgen genes
    for gene in ANDROGEN_GENES:
        zcol = f"{gene}_ZSCORE"
        if zcol in df.columns:
            df["_feat_z"] = pd.to_numeric(df[zcol], errors="coerce")
            _cox_pair(df, "_feat_z", f"High {gene} (continuous z)",
                      "task_pfs_tcga", gene, "TCGA_PRAD", "ZSCORE_CONTINUOUS",
                      "continuous_zscore", covars, rows,
                      time_col="PFS_MONTHS", event_col="PFS_EVENT", endpoint="PFS")
        for split in SPLIT_NAMES:
            gcol = f"{gene}_GROUP_{split}"
            if gcol not in df.columns:
                continue
            df["_feat_b"] = _group_to_cox_binary(df[gcol], "high_is_bad")
            log.info(f"    {gene} {split} PFS")
            _cox_pair(df, "_feat_b", f"High {gene}",
                      "task_pfs_tcga", gene, "TCGA_PRAD", split,
                      "binary_group", covars, rows,
                      time_col="PFS_MONTHS", event_col="PFS_EVENT", endpoint="PFS")

    # Adhesion/motility genes
    for gene, cfg in ADHESION_MOTILITY_GENES.items():
        direction = cfg["direction"]
        label_bad = f"Low {gene} (epithelial loss)" if direction == "low_is_bad" else f"High {gene}"
        zcol = f"{gene}_ZSCORE"
        if zcol in df.columns:
            df["_feat_z"] = pd.to_numeric(df[zcol], errors="coerce")
            _cox_pair(df, "_feat_z", f"{label_bad} (continuous z)",
                      "task_pfs_tcga", gene, "TCGA_PRAD", "ZSCORE_CONTINUOUS",
                      "continuous_zscore", covars, rows,
                      time_col="PFS_MONTHS", event_col="PFS_EVENT", endpoint="PFS")
        for split in SPLIT_NAMES:
            gcol = f"{gene}_GROUP_{split}"
            if gcol not in df.columns:
                continue
            df["_feat_b"] = _group_to_cox_binary(df[gcol], direction)
            log.info(f"    {gene} {split} PFS")
            _cox_pair(df, "_feat_b", label_bad,
                      "task_pfs_tcga", gene, "TCGA_PRAD", split,
                      "binary_group", covars, rows,
                      time_col="PFS_MONTHS", event_col="PFS_EVENT", endpoint="PFS")


# ─────────────────────────────────────────────────────────────────────────────
# TASK 3b — Individual androgen CNA genes, 7 cohorts combined
# ─────────────────────────────────────────────────────────────────────────────

def cox_task3b(mol_pool, covars, rows):
    log.info("\n--- Cox Task 3b: Individual androgen CNA ---")
    for gene in ANDROGEN_CNA_GENES:
        col = f"{gene}_CNA_ALT"
        if col not in mol_pool.columns:
            continue
        mol_pool[f"_feat_{gene}"] = _to_binary(mol_pool[col])
        log.info(f"  {gene}")
        _cox_pair(mol_pool, f"_feat_{gene}", f"{gene} Altered (Gain/Mut)",
                  "task3b", gene, "7_cohorts", "", "binary_alteration",
                  covars, rows)


# ─────────────────────────────────────────────────────────────────────────────
# TASK 4 — OR-combined androgen mRNA, per cohort × 3 splits
# ─────────────────────────────────────────────────────────────────────────────

def cox_task4(rows):
    log.info("\n--- Cox Task 4: OR-combined androgen mRNA per cohort ---")
    for cohort_key in EXPRESSION_COHORT_KEYS:
        df, time_col, event_col, ep_label = _load_mrna_cohort(cohort_key)
        if df is None:
            continue
        covars = _mrna_cohort_covariates(df)
        cohort_label = DATASETS[cohort_key]["label"]
        log.info(f"  Cohort: {cohort_label}")

        for split in SPLIT_NAMES:
            col = f"ANDROGEN_OR_{split}"
            if col not in df.columns:
                continue
            df["_feat_or"] = pd.to_numeric(df[col], errors="coerce")
            log.info(f"    {split}")
            _cox_pair(df, "_feat_or", "Any Androgen Gene Upregulated",
                      "task4", "ANDROGEN_OR", cohort_key, split,
                      "binary_or", covars, rows,
                      time_col=time_col, event_col=event_col,
                      endpoint=ep_label)


# ─────────────────────────────────────────────────────────────────────────────
# TASK 6 — Individual adhesion/motility mRNA genes, per cohort
# Continuous z-score + binary per split (direction-aware coding)
# ─────────────────────────────────────────────────────────────────────────────

def cox_task6(rows):
    log.info("\n--- Cox Task 6: Individual adhesion/motility mRNA per cohort ---")
    for cohort_key in EXPRESSION_COHORT_KEYS:
        df, time_col, event_col, ep_label = _load_mrna_cohort(cohort_key)
        if df is None:
            continue
        covars = _mrna_cohort_covariates(df)
        cohort_label = DATASETS[cohort_key]["label"]
        log.info(f"  Cohort: {cohort_label}")

        for gene, cfg in ADHESION_MOTILITY_GENES.items():
            direction = cfg["direction"]

            # Continuous z-score (raw z, HR interpretation depends on direction)
            zcol = f"{gene}_ZSCORE"
            if zcol in df.columns:
                df["_feat_z"] = pd.to_numeric(df[zcol], errors="coerce")
                label_note = "high is bad" if direction == "high_is_bad" else "low is bad"
                _cox_pair(df, "_feat_z", f"{gene} z-score ({label_note})",
                          "task6", gene, cohort_key, "ZSCORE_CONTINUOUS",
                          "continuous_zscore", covars, rows,
                          time_col=time_col, event_col=event_col,
                          endpoint=ep_label)

            # Binary per split — bad group coded as 1
            for split in SPLIT_NAMES:
                gcol = f"{gene}_GROUP_{split}"
                if gcol not in df.columns:
                    continue
                df["_feat_b"] = _group_to_cox_binary(df[gcol], direction)
                label_bad = (f"Low {gene} (epithelial loss)"
                             if direction == "low_is_bad"
                             else f"High {gene}")
                log.info(f"    {gene} {split}")
                _cox_pair(df, "_feat_b", label_bad,
                          "task6", gene, cohort_key, split,
                          "binary_group", covars, rows,
                          time_col=time_col, event_col=event_col,
                          endpoint=ep_label)


# ─────────────────────────────────────────────────────────────────────────────
# TASK 6b — Individual adhesion/motility CNA genes, 7 cohorts combined
# ─────────────────────────────────────────────────────────────────────────────

def cox_task6b(mol_pool, covars, rows):
    log.info("\n--- Cox Task 6b: Individual adhesion/motility CNA ---")
    for gene, cfg in ADHESION_CNA_GENES.items():
        col = f"{gene}_CNA_ALT"
        if col not in mol_pool.columns:
            continue
        mol_pool[f"_feat_{gene}"] = _to_binary(mol_pool[col])
        alt_note = "loss" if cfg["direction"] == "loss_is_bad" else "gain"
        log.info(f"  {gene}")
        _cox_pair(mol_pool, f"_feat_{gene}", f"{gene} Altered ({alt_note})",
                  "task6b", gene, "7_cohorts", "", "binary_alteration",
                  covars, rows)


# ─────────────────────────────────────────────────────────────────────────────
# TASK 7 — OR-combined adhesion/motility mRNA, per cohort × 3 splits
# ─────────────────────────────────────────────────────────────────────────────

def cox_task7(rows):
    log.info("\n--- Cox Task 7: OR-combined adhesion/motility mRNA per cohort ---")
    for cohort_key in EXPRESSION_COHORT_KEYS:
        df, time_col, event_col, ep_label = _load_mrna_cohort(cohort_key)
        if df is None:
            continue
        covars = _mrna_cohort_covariates(df)
        cohort_label = DATASETS[cohort_key]["label"]
        log.info(f"  Cohort: {cohort_label}")

        for split in SPLIT_NAMES:
            col = f"ADHESION_OR_{split}"
            if col not in df.columns:
                continue
            df["_feat_or"] = pd.to_numeric(df[col], errors="coerce")
            log.info(f"    {split}")
            _cox_pair(df, "_feat_or", "Any Adhesion/Motility Gene Dysregulated",
                      "task7", "ADHESION_OR", cohort_key, split,
                      "binary_or", covars, rows,
                      time_col=time_col, event_col=event_col,
                      endpoint=ep_label)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    log.info("\n" + "="*70)
    log.info("MODULE 7: COX PROPORTIONAL HAZARDS REGRESSION")
    log.info("="*70)

    rows = []

    # Load 7-cohort combined once (reused for tasks 1, 2, 3b, 6b)
    log.info("\nLoading 7-cohort combined data...")
    try:
        mol_pool    = _load_molecular_pool()
        covars = _molecular_pool_covariates(mol_pool)
        log.info(f"  7-cohort covariates: {covars}")

        cox_task1(mol_pool,  covars, rows)
        cox_task2(mol_pool,  covars, rows)
        cox_task3b(mol_pool, covars, rows)
        cox_task6b(mol_pool, covars, rows)
    except FileNotFoundError as e:
        log.error(f"  7-cohort load failed: {e}")

    # Per-cohort mRNA tasks (load inside each function)
    cox_task3(rows)
    cox_task_pfs_tcga(rows)
    cox_task4(rows)
    cox_task6(rows)
    cox_task7(rows)

    # Save
    out_path = TABLES_DIR / "cox_results.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)
    log.info(f"\nSaved: {out_path}  ({len(rows)} rows)")

    log.info("\n" + "="*70)
    log.info("MODULE 7 COMPLETE")
    log.info("="*70)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--refresh", action="store_true")
    args = p.parse_args()
    main()
    log.info("Cox regression complete.")
