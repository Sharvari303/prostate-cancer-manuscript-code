"""
module3_clean.py — Standardizes and validates all cohort data.

Renames cohort-specific column names to standard labels (e.g.
"Overall Survival (Months)" → OS_MONTHS), converts string status
values to binary events (1=event, 0=censored), and removes rows
that have no valid survival endpoint at all.

Also produces a combined CSV stacking all cohorts for reference.

Output: data/cache/{KEY}_clean.csv  (one per cohort)
        data/cache/combined_clean.csv
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATASETS, EVENT_POSITIVE_STRINGS, EVENT_NEGATIVE_STRINGS, CACHE_DIR, AGE_COLUMN_ALIASES
from utils.logger import get_logger

log = get_logger("module3_clean")

COLUMN_ALIASES = {
    # Survival endpoints
    "OS_MONTHS": "OS_MONTHS",
    "Overall Survival (Months)": "OS_MONTHS",
    # MCSPC_MSK uses OS_SMP_MONTHS for overall survival time
    "OS_SMP_MONTHS": "OS_MONTHS",
    "OS_STATUS": "OS_STATUS",
    "Overall Survival Status": "OS_STATUS",
    # MCSPC_MSK uses SURVIVAL_STATUS instead of OS_STATUS
    "SURVIVAL_STATUS": "OS_STATUS",
    "DFS_MONTHS": "DFS_MONTHS",
    "Disease Free (Survival) (Months)": "DFS_MONTHS",
    "Disease Free Months": "DFS_MONTHS",
    "DFS_STATUS": "DFS_STATUS",
    "Disease Free Status": "DFS_STATUS",
    "PFS_MONTHS": "PFS_MONTHS",
    "Progression Free Survival (Months)": "PFS_MONTHS",
    "PFS_STATUS": "PFS_STATUS",
    "Progression Free Status": "PFS_STATUS",
    # Age columns — standardize all variants to AGE for consistent Cox covariate use
    **AGE_COLUMN_ALIASES,
}

ENDPOINT_PAIRS = [
    ("OS_MONTHS",  "OS_EVENT",  "Overall Survival"),
    ("DFS_MONTHS", "DFS_EVENT", "Disease-Free Survival"),
    ("PFS_MONTHS", "PFS_EVENT", "Progression-Free Survival"),
]


def _standardize_columns(df):
    renames = {}
    # Initialize with existing column names so we never rename two columns to
    # the same target (e.g. CURRENT_AGE_DEID + DXAGE both → AGE).
    seen_targets = set(df.columns)
    for c in df.columns:
        if c in COLUMN_ALIASES:
            target = COLUMN_ALIASES[c]
            if target == c:
                continue
            if target in seen_targets:
                continue
            renames[c] = target
            seen_targets.add(target)
    return df.rename(columns=renames)


def _status_to_event(series, status_col):
    positives = EVENT_POSITIVE_STRINGS.get(status_col, [])
    negatives = EVENT_NEGATIVE_STRINGS.get(status_col, [])
    def _map(val):
        if pd.isna(val):
            return np.nan
        s = str(val).strip()
        if s in positives: return 1
        if s in negatives: return 0
        if s.startswith("1:"): return 1
        if s.startswith("0:"): return 0
        try: return int(float(s))
        except ValueError: return np.nan
    return series.map(_map)


def _add_event_columns(df):
    status_map = {"OS_STATUS": "OS_EVENT",
                  "DFS_STATUS": "DFS_EVENT",
                  "PFS_STATUS": "PFS_EVENT"}
    for scol, ecol in status_map.items():
        if scol in df.columns:
            df[ecol] = _status_to_event(df[scol], scol)
            n_ev = (df[ecol] == 1).sum()
            n_ce = (df[ecol] == 0).sum()
            log.info(f"  {ecol}: events={n_ev}, censored={n_ce}, "
                     f"missing={df[ecol].isna().sum()}")
    return df


def _deduplicate_suffixed_columns(df):
    """
    Handles duplicate columns from merges with _x and _y suffixes.
    For genes present in both AR activity and androgen/adhesion fetches,
    keep the _x version (AR activity priority) and drop _y.
    Returns cleaned dataframe with single _ZSCORE columns.
    """
    cols_to_drop = []

    for col in df.columns:
        if col.endswith('_x') and col[:-2] + '_y' in df.columns:
            base_col = col[:-2]  # Remove _x suffix
            if base_col.endswith('_ZSCORE'):
                x_vals = df[col].copy()
                y_vals = df[col[:-2] + '_y'].copy()
                merged = x_vals.fillna(y_vals)
                df[base_col] = merged
                cols_to_drop.append(col)
                cols_to_drop.append(col[:-2] + '_y')

    df = df.drop(columns=cols_to_drop, errors='ignore')
    return df


def _coerce_numeric(df):
    for col in df.columns:
        if any(p in col for p in ["_MONTHS", "_CNA", "_ZSCORE"]):
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _remove_invalid_rows(df):
    n_start = len(df)
    # Zero out negative times
    for col in [c for c in df.columns if "_MONTHS" in c]:
        df.loc[df[col] < 0, col] = np.nan
    # Keep rows with at least one valid endpoint
    has_valid = pd.Series(False, index=df.index)
    for tc, ec, _ in ENDPOINT_PAIRS:
        if tc in df.columns and ec in df.columns:
            has_valid |= (df[tc].notna() & df[ec].notna() & (df[tc] > 0))
    df = df[has_valid].copy()
    log.info(f"  Rows: {n_start} → {len(df)} (removed {n_start-len(df)})")
    return df


def clean_dataset(master_df, dataset_key):
    label = DATASETS[dataset_key]["label"]
    log.info(f"\n{'='*55}\nCleaning: {label}\n{'='*55}")
    df = master_df.copy()
    df = _standardize_columns(df)
    df = _deduplicate_suffixed_columns(df)
    df = _coerce_numeric(df)
    df = _add_event_columns(df)
    df = _remove_invalid_rows(df)
    if "DATASET" not in df.columns:
        df["DATASET"] = label
    # Endpoint summary
    for tc, ec, lbl in ENDPOINT_PAIRS:
        if tc in df.columns and ec in df.columns:
            v = df[tc].notna() & df[ec].notna()
            n_ev = (df.loc[v, ec] == 1).sum()
            med  = df.loc[v, tc].median()
            log.info(f"  {lbl:<30} n={v.sum()}, "
                     f"events={n_ev}, median={med:.1f}mo")
    return df


def combine_datasets(clean_dfs):
    # Guard: deduplicate columns in each DataFrame before concat to prevent
    # "Reindexing only valid with uniquely valued Index objects" crash.
    deduped = {}
    for key, df in clean_dfs.items():
        if df.columns.duplicated().any():
            dup = df.columns[df.columns.duplicated()].tolist()
            log.warning(f"  {key}: dropping duplicate columns: {dup}")
            df = df.loc[:, ~df.columns.duplicated(keep='first')]
        deduped[key] = df
    combined = pd.concat(list(deduped.values()), axis=0, ignore_index=True)
    log.info(f"\nCombined: {len(combined)} total patients")
    for key, df in clean_dfs.items():
        lbl = DATASETS[key]["label"]
        log.info(f"  {lbl}: {(combined['DATASET']==lbl).sum()}")
    return combined


def run_cleaning(master_dfs=None, force_refresh=False):
    combined_cache = CACHE_DIR / "combined_clean.csv"
    if combined_cache.exists() and not force_refresh:
        log.info("Loading clean datasets from cache...")
        results = {}
        for key in DATASETS:
            p = CACHE_DIR / f"{key}_clean.csv"
            results[key] = pd.read_csv(p)
        results["COMBINED"] = pd.read_csv(combined_cache)
        return results

    if master_dfs is None:
        master_dfs = {}
        for key in DATASETS:
            p = CACHE_DIR / f"{key}_master.csv"
            if not p.exists():
                raise FileNotFoundError(
                    f"Missing: {p}. Run module2_download.py first.")
            master_dfs[key] = pd.read_csv(p)

    clean_dfs = {}
    for key, df in master_dfs.items():
        clean = clean_dataset(df, key)
        clean_dfs[key] = clean
        out = CACHE_DIR / f"{key}_clean.csv"
        clean.to_csv(out, index=False)
        log.info(f"Saved: {out}")

    combined = combine_datasets(clean_dfs)
    combined.to_csv(combined_cache, index=False)
    log.info(f"Saved: {combined_cache}")

    clean_dfs["COMBINED"] = combined
    return clean_dfs


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--refresh", action="store_true")
    args = p.parse_args()
    results = run_cleaning(force_refresh=args.refresh)
    log.info("Done.")
