"""
module6_comprehensive_km.py — Kaplan-Meier survival analyses for all tasks.

Reads flagged CSVs produced by module4 (molecular flags) and module5 (expression flags).

Endpoint routing (defined in config.ENDPOINTS):
  Molecular pool (7 cohorts combined): OS throughout
  mRNA cohorts (per-cohort):
    TCGA_PRAD → DFS   SU2C → OS   MSKCC → DFS   MCTP → OS

Tasks:
  1:  Individual molecular genes (PTEN/MDM2/MDM4/AR/CDKN1A/BAX/ATM)
      7 cohorts combined, binary altered vs wildtype (OS)
  2:  OR-combined (PTEN+AR+MDM4), 7 cohorts combined (OS)
  3:  Individual androgen mRNA genes (SLCO2B1/SLCO1B3/AKR1C3)
      per mRNA cohort × 3 splits (MEDIAN/QUARTILE/ZSCORE)
  3pfs: AKR1C3 mRNA vs PFS — TCGA-PRAD × 3 splits (PFS endpoint)
  3b: Individual androgen CNA genes, 7 cohorts combined (OS)
  4:  OR-combined androgen mRNA, per mRNA cohort × 3 splits
  6:  Individual adhesion/motility mRNA genes, per mRNA cohort × 3 splits
      direction-aware: epithelial LOW is bad, mesenchymal HIGH is bad
  6b: Individual adhesion/motility CNA genes, 7 cohorts combined (OS)
  7:  OR-combined adhesion/motility mRNA, per mRNA cohort × 3 splits

Missing data rule: NaN → excluded (never assigned to wildtype).
Min group size: 20 patients per KM arm.

Output:
  outputs/tables/km_statistics.csv   — one row per figure
  outputs/tables/median_os_table.csv — one row per KM arm per figure
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    DATASETS, MOLECULAR_COHORT_KEYS, EXPRESSION_COHORT_KEYS,
    MOLECULAR_GENES, INDIVIDUAL_GENES,
    ANDROGEN_GENES, ANDROGEN_CNA_GENES,
    ADHESION_MOTILITY_GENES, ADHESION_CNA_GENES,
    ENDPOINTS,
    CACHE_DIR, TABLES_DIR,
    COLORS, PLOT_STYLE, STATS,
    DEFAULT_CNA_THRESHOLD,
)
from utils.km_engine import (
    apply_plot_style, save_figure,
    logrank_test, format_pvalue,
    plot_two_group_km,
)
from utils.logger import get_logger

log = get_logger("module6_comprehensive_km")

MIN_GROUP  = STATS["min_group_size"]   # 20 patients per arm

# Default for molecular pool (OS throughout for all 7 cohorts)
TIME_COL   = "OS_MONTHS"
EVENT_COL  = "OS_EVENT"

SPLIT_NAMES = ["MEDIAN", "QUARTILE", "ZSCORE"]


def _get_endpoint(cohort_key):
    """Return (time_col, event_col, endpoint_label) for a given cohort.

    Uses config.ENDPOINTS primary entry per cohort. Falls back to OS.
    This ensures TCGA_PRAD and MSKCC use DFS while SU2C/MCTP use OS.
    """
    ep = ENDPOINTS.get(cohort_key, ENDPOINTS["DEFAULT"])["primary"]
    return ep["time"], ep["event"], ep["short"]
SPLIT_LABELS = {
    "MEDIAN":   "Median split",
    "QUARTILE": "Quartile split (top 25% vs bottom 25%)",
    "ZSCORE":   "z > 1.0 vs rest",
}


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
    """Coerce survival columns to numeric and drop rows missing time or event.

    Accepts explicit time_col/event_col so mRNA cohorts using DFS can pass
    DFS_MONTHS/DFS_EVENT rather than always using OS.
    """
    df = df.copy()
    df[time_col]  = pd.to_numeric(df.get(time_col,  pd.Series(dtype=float)), errors="coerce")
    df[event_col] = pd.to_numeric(df.get(event_col, pd.Series(dtype=float)), errors="coerce")
    df = df.dropna(subset=[time_col, event_col])
    df = df[df[time_col] > 0].copy()
    return df


def _clean_os(df):
    """Legacy alias — uses OS endpoint (for molecular pool tasks)."""
    return _clean_endpoint(df, TIME_COL, EVENT_COL)


# ─────────────────────────────────────────────────────────────────────────────
# CORE KM RUNNER
# ─────────────────────────────────────────────────────────────────────────────

def _run_binary_km(df_a, df_b, label_a, label_b,
                   title, figure_id, task, gene, cohort, split,
                   subdir, stats_rows, median_rows,
                   color_a=None, color_b=None,
                   time_col=None, event_col=None, endpoint_label=None):
    """
    One binary KM comparison: df_a = 'bad/altered' arm, df_b = 'good/wildtype' arm.
    Skips if either arm has < MIN_GROUP patients with valid endpoint data.
    Saves figure; appends one row to stats_rows and two rows to median_rows.

    time_col/event_col: override default OS columns for cohorts using DFS.
    endpoint_label: short string for axis labels (e.g. "OS" or "DFS").
    """
    tc = time_col  or TIME_COL
    ec = event_col or EVENT_COL
    ep = endpoint_label or "OS"

    df_a_v = _clean_endpoint(df_a, tc, ec)
    df_b_v = _clean_endpoint(df_b, tc, ec)
    na, nb = len(df_a_v), len(df_b_v)

    if na < MIN_GROUP or nb < MIN_GROUP:
        log.warning(f"  {figure_id}: skip — n_a={na}, n_b={nb} (min={MIN_GROUP})")
        return

    pval = logrank_test(
        df_a_v[tc].values, df_a_v[ec].values,
        df_b_v[tc].values, df_b_v[ec].values,
    )

    apply_plot_style()
    fig, ax = plt.subplots(figsize=(9, 6))

    km_a, km_b = plot_two_group_km(
        ax, df_a_v, df_b_v, tc, ec,
        label_a, label_b,
        color_high=color_a or COLORS["high"],
        color_low=color_b  or COLORS["low"],
        title=title,
        pval=pval,
        ylabel=f"{ep} Probability",
    )

    save_figure(fig, figure_id, subdir)   # closes figure

    if km_a is None or km_b is None:
        return

    pval_str = format_pvalue(pval)
    log.info(f"  {figure_id}: n_a={na}(ev={km_a.n_events}), "
             f"n_b={nb}(ev={km_b.n_events}), {pval_str}")

    stats_rows.append({
        "figure_id": figure_id, "task": task, "gene": gene,
        "cohort": cohort, "split": split, "endpoint": ep,
        "N_group1": na, "N_group2": nb,
        "events_group1": km_a.n_events, "events_group2": km_b.n_events,
        "pval": pval,
    })
    for lbl, km, n in [(label_a, km_a, na), (label_b, km_b, nb)]:
        median_rows.append({
            "figure_id": figure_id, "task": task, "gene": gene,
            "cohort": cohort, "split": split, "endpoint": ep,
            "group": lbl, "N": n, "events": km.n_events,
            "median_months": km.median_survival,
            "CI_lower": km.median_ci_lower, "CI_upper": km.median_ci_upper,
        })


# ─────────────────────────────────────────────────────────────────────────────
# DATA LOADING
# ─────────────────────────────────────────────────────────────────────────────

def _load_molecular_pool():
    """Load and stack flagged molecular CSVs for the 7-cohort molecular OS pool.

    All 7 cohorts use OS endpoint. MSKCC and IDH_MUTANT are excluded
    (no OS and biologically atypical, respectively — see config.DATASETS).
    """
    dfs = []
    for i, key in enumerate(MOLECULAR_COHORT_KEYS):
        path = CACHE_DIR / f"{key}_flagged_molecular.csv"
        if not path.exists():
            log.warning(f"  Missing: {path.name} — skipping {key}")
            continue
        df = pd.read_csv(path, low_memory=False)
        df["COHORT_KEY"]  = key
        df["COHORT_CODE"] = i
        dfs.append(df)
    if not dfs:
        raise FileNotFoundError("No flagged molecular files found for 7-cohort set")
    combined = pd.concat(dfs, ignore_index=True)
    combined = _clean_os(combined)
    log.info(f"  7-cohort combined: {len(combined)} patients with valid OS")
    return combined


def _load_mrna_cohorts():
    """Load flagged expression CSVs for the 4 mRNA cohorts.

    Each cohort is loaded with its own primary endpoint (from config.ENDPOINTS):
      TCGA_PRAD → DFS_MONTHS / DFS_EVENT  (OS has only 10 events)
      MSKCC     → DFS_MONTHS / DFS_EVENT  (no OS data)
      SU2C      → OS_MONTHS  / OS_EVENT
      MCTP      → OS_MONTHS  / OS_EVENT

    Returns dict: cohort_key → (df, time_col, event_col, endpoint_label)
    """
    dfs = {}
    for key in EXPRESSION_COHORT_KEYS:
        path = CACHE_DIR / f"{key}_flagged_expression.csv"
        if not path.exists():
            log.warning(f"  Missing: {path.name} — skipping {key}")
            continue
        df = pd.read_csv(path, low_memory=False)
        time_col, event_col, ep_label = _get_endpoint(key)
        df_clean = _clean_endpoint(df, time_col, event_col)
        log.info(f"  {key}: {len(df_clean)} patients with valid {ep_label} "
                 f"({time_col}/{event_col})")
        dfs[key] = (df_clean, time_col, event_col, ep_label)
    return dfs


# ─────────────────────────────────────────────────────────────────────────────
# TASK 1 — Individual molecular genes, 7 cohorts combined
# 7 genes: PTEN, MDM2, MDM4, AR, CDKN1A, BAX, ATM
# ─────────────────────────────────────────────────────────────────────────────

def run_task1(mol_pool, stats_rows, median_rows):
    log.info("\n" + "="*60)
    log.info("TASK 1: Individual molecular genes — 7 cohorts combined OS")
    log.info("="*60)

    all_genes = list(MOLECULAR_GENES.keys()) + list(INDIVIDUAL_GENES.keys())

    for gene in all_genes:
        col = f"{gene}_ALT_{DEFAULT_CNA_THRESHOLD.upper()}"
        if col not in mol_pool.columns:
            log.warning(f"  {col} not found — skipping {gene}")
            continue

        flags = _to_binary(mol_pool[col])
        df_alt = mol_pool[flags == 1.0].copy()
        df_wt  = mol_pool[flags == 0.0].copy()

        title     = f"{gene} — Altered vs Wildtype | Combined Molecular Cohorts (OS)"
        figure_id = f"{gene.lower()}_7cohorts_os"

        _run_binary_km(
            df_alt, df_wt,
            f"{gene} Altered", f"{gene} Wildtype",
            title=title, figure_id=figure_id,
            task="task1", gene=gene, cohort="7_cohorts", split="",
            subdir="task1_molecular_individual",
            stats_rows=stats_rows, median_rows=median_rows,
            color_a=COLORS["altered"], color_b=COLORS["unaltered"],
        )


# ─────────────────────────────────────────────────────────────────────────────
# TASK 2 — OR-combined (PTEN + AR + MDM4), 7 cohorts combined
# 1 figure total
# ─────────────────────────────────────────────────────────────────────────────

def run_task2(mol_pool, stats_rows, median_rows):
    log.info("\n" + "="*60)
    log.info("TASK 2: OR-combined (PTEN+AR+MDM4) — 7 cohorts combined OS")
    log.info("="*60)

    col = "MOLECULAR_OR"
    if col not in mol_pool.columns:
        log.warning(f"  {col} not found — skipping task2")
        return

    flags  = pd.to_numeric(mol_pool[col], errors="coerce")
    df_alt = mol_pool[flags == 1].copy()
    df_wt  = mol_pool[flags == 0].copy()

    _run_binary_km(
        df_alt, df_wt,
        "Any Altered (PTEN/AR/MDM4)", "All Wildtype",
        title="PTEN OR AR OR MDM4 — Any Altered vs All Wildtype | Combined Molecular Cohorts (OS)",
        figure_id="pten_ar_mdm4_or_7cohorts_os",
        task="task2", gene="PTEN_AR_MDM4", cohort="7_cohorts", split="",
        subdir="task2_molecular_or",
        stats_rows=stats_rows, median_rows=median_rows,
        color_a=COLORS["altered"], color_b=COLORS["unaltered"],
    )


# ─────────────────────────────────────────────────────────────────────────────
# TASK 3 — Individual androgen mRNA genes, per cohort × 3 splits
# 3 cohorts × 3 genes × 3 splits = 27 figures
# ─────────────────────────────────────────────────────────────────────────────

def run_task3(mrna_dfs, stats_rows, median_rows):
    log.info("\n" + "="*60)
    log.info("TASK 3: Individual androgen mRNA — per cohort × 3 splits")
    log.info("="*60)

    for cohort_key, (df, time_col, event_col, ep_label) in mrna_dfs.items():
        cohort_label = DATASETS[cohort_key]["label"]
        log.info(f"\n  Cohort: {cohort_label}")

        for gene in ANDROGEN_GENES:
            for split in SPLIT_NAMES:
                col = f"{gene}_GROUP_{split}"
                if col not in df.columns:
                    log.warning(f"    {col} not found")
                    continue

                df_high = df[df[col] == "High"].copy()
                df_low  = df[df[col] == "Low"].copy()

                figure_id = (f"{gene.lower()}_{cohort_key.lower()}"
                             f"_{split.lower()}_{ep_label.lower()}")
                title     = (f"{gene} — Androgen Uptake | {SPLIT_LABELS[split]}\n"
                             f"{cohort_label} ({ep_label})")

                _run_binary_km(
                    df_high, df_low,
                    f"High {gene}", f"Low {gene}",
                    title=title, figure_id=figure_id,
                    task="task3", gene=gene,
                    cohort=cohort_key, split=split,
                    subdir="task3_androgen_mrna_individual",
                    stats_rows=stats_rows, median_rows=median_rows,
                    color_a=COLORS["high"], color_b=COLORS["low"],
                    time_col=time_col, event_col=event_col, endpoint_label=ep_label,
                )


# ─────────────────────────────────────────────────────────────────────────────
# TASK 3b — Individual androgen CNA genes, 7 cohorts combined
# 3 genes: SLCO2B1, SLCO1B3, AKR1C3 — gain is bad
# ─────────────────────────────────────────────────────────────────────────────

def run_task3b(mol_pool, stats_rows, median_rows):
    log.info("\n" + "="*60)
    log.info("TASK 3b: Individual androgen CNA — 7 cohorts combined OS")
    log.info("="*60)

    for gene in ANDROGEN_CNA_GENES:
        col = f"{gene}_CNA_ALT"
        if col not in mol_pool.columns:
            log.warning(f"  {col} not found — skipping {gene}")
            continue

        flags  = _to_binary(mol_pool[col])
        df_alt = mol_pool[flags == 1.0].copy()
        df_wt  = mol_pool[flags == 0.0].copy()

        title     = f"{gene} CNA/Mutation — Altered vs Wildtype | Combined Molecular Cohorts (OS)"
        figure_id = f"{gene.lower()}_cna_7cohorts_os"

        _run_binary_km(
            df_alt, df_wt,
            f"{gene} Altered (Gain/Mut)", f"{gene} Wildtype",
            title=title, figure_id=figure_id,
            task="task3b", gene=gene, cohort="7_cohorts", split="",
            subdir="task3b_androgen_cna",
            stats_rows=stats_rows, median_rows=median_rows,
            color_a=COLORS["altered"], color_b=COLORS["unaltered"],
        )


# ─────────────────────────────────────────────────────────────────────────────
# TASK 4 — OR-combined androgen mRNA, per cohort × 3 splits
# 3 cohorts × 3 splits = 9 figures
# ─────────────────────────────────────────────────────────────────────────────

def run_task4(mrna_dfs, stats_rows, median_rows):
    log.info("\n" + "="*60)
    log.info("TASK 4: OR-combined androgen mRNA — per cohort × 3 splits")
    log.info("="*60)

    for cohort_key, (df, time_col, event_col, ep_label) in mrna_dfs.items():
        cohort_label = DATASETS[cohort_key]["label"]
        log.info(f"\n  Cohort: {cohort_label}")

        for split in SPLIT_NAMES:
            col = f"ANDROGEN_OR_{split}"
            if col not in df.columns:
                log.warning(f"    {col} not found")
                continue

            flags   = pd.to_numeric(df[col], errors="coerce")
            df_up   = df[flags == 1].copy()
            df_none = df[flags == 0].copy()

            figure_id = (f"androgen_or_{cohort_key.lower()}"
                         f"_{split.lower()}_{ep_label.lower()}")
            title     = (f"Any Androgen Gene Upregulated (SLCO2B1/SLCO1B3/AKR1C3)\n"
                         f"{SPLIT_LABELS[split]} | {cohort_label} ({ep_label})")

            _run_binary_km(
                df_up, df_none,
                "Any Upregulated", "None Upregulated",
                title=title, figure_id=figure_id,
                task="task4", gene="ANDROGEN_OR",
                cohort=cohort_key, split=split,
                subdir="task4_androgen_or",
                stats_rows=stats_rows, median_rows=median_rows,
                color_a=COLORS["high"], color_b=COLORS["low"],
                time_col=time_col, event_col=event_col, endpoint_label=ep_label,
            )


# ─────────────────────────────────────────────────────────────────────────────
# TASK 6 — Individual adhesion/motility mRNA genes, per cohort × 3 splits
# 3 cohorts × 10 genes × 3 splits = 90 figures
# Direction-aware: epithelial LOW is bad, mesenchymal HIGH is bad
# ─────────────────────────────────────────────────────────────────────────────

def run_task6(mrna_dfs, stats_rows, median_rows):
    log.info("\n" + "="*60)
    log.info("TASK 6: Individual adhesion/motility mRNA — per cohort × 3 splits")
    log.info("="*60)

    for cohort_key, (df, time_col, event_col, ep_label) in mrna_dfs.items():
        cohort_label = DATASETS[cohort_key]["label"]
        log.info(f"\n  Cohort: {cohort_label}")

        for gene, cfg in ADHESION_MOTILITY_GENES.items():
            direction = cfg["direction"]

            for split in SPLIT_NAMES:
                col = f"{gene}_GROUP_{split}"
                if col not in df.columns:
                    log.warning(f"    {col} not found")
                    continue

                if direction == "low_is_bad":
                    df_bad  = df[df[col] == "Low"].copy()
                    df_good = df[df[col] == "High"].copy()
                    label_bad  = f"Low {gene} (epithelial loss)"
                    label_good = f"High {gene}"
                    color_bad  = COLORS["altered"]
                    color_good = COLORS["unaltered"]
                    note = "epithelial loss"
                else:  # high_is_bad
                    df_bad  = df[df[col] == "High"].copy()
                    df_good = df[df[col] == "Low"].copy()
                    label_bad  = f"High {gene}"
                    label_good = f"Low {gene}"
                    color_bad  = COLORS["high"]
                    color_good = COLORS["low"]
                    note = "mesenchymal gain"

                figure_id = (f"{gene.lower()}_{cohort_key.lower()}"
                             f"_{split.lower()}_{ep_label.lower()}")
                title     = (f"{gene} ({note}) | {SPLIT_LABELS[split]}\n"
                             f"{cohort_label} ({ep_label})")

                _run_binary_km(
                    df_bad, df_good,
                    label_bad, label_good,
                    title=title, figure_id=figure_id,
                    task="task6", gene=gene,
                    cohort=cohort_key, split=split,
                    subdir="task6_adhesion_mrna_individual",
                    stats_rows=stats_rows, median_rows=median_rows,
                    color_a=color_bad, color_b=color_good,
                    time_col=time_col, event_col=event_col, endpoint_label=ep_label,
                )


# ─────────────────────────────────────────────────────────────────────────────
# TASK 6b — Individual adhesion/motility CNA genes, 7 cohorts combined
# 10 genes — CDH1/EPCAM: loss is bad; VIM/CDH2/MMP9/CXCL8/ACTB/RHOA/ITGB1/FN1: gain is bad
# ─────────────────────────────────────────────────────────────────────────────

def run_task6b(mol_pool, stats_rows, median_rows):
    log.info("\n" + "="*60)
    log.info("TASK 6b: Individual adhesion/motility CNA — 7 cohorts combined OS")
    log.info("="*60)

    for gene, cfg in ADHESION_CNA_GENES.items():
        col = f"{gene}_CNA_ALT"
        if col not in mol_pool.columns:
            log.warning(f"  {col} not found — skipping {gene}")
            continue

        flags  = _to_binary(mol_pool[col])
        df_alt = mol_pool[flags == 1.0].copy()
        df_wt  = mol_pool[flags == 0.0].copy()

        alt_note  = "loss" if cfg["direction"] == "loss_is_bad" else "gain"
        figure_id = f"{gene.lower()}_cna_7cohorts_os"
        title     = f"{gene} CNA/Mutation ({alt_note}) — Altered vs Wildtype | Combined Molecular Cohorts (OS)"

        _run_binary_km(
            df_alt, df_wt,
            f"{gene} Altered ({alt_note})", f"{gene} Wildtype",
            title=title, figure_id=figure_id,
            task="task6b", gene=gene, cohort="7_cohorts", split="",
            subdir="task6b_adhesion_cna",
            stats_rows=stats_rows, median_rows=median_rows,
            color_a=COLORS["altered"], color_b=COLORS["unaltered"],
        )


# ─────────────────────────────────────────────────────────────────────────────
# TASK 7 — OR-combined adhesion/motility mRNA, per cohort × 3 splits
# 3 cohorts × 3 splits = 9 figures
# ─────────────────────────────────────────────────────────────────────────────

def run_task7(mrna_dfs, stats_rows, median_rows):
    log.info("\n" + "="*60)
    log.info("TASK 7: OR-combined adhesion/motility mRNA — per cohort × 3 splits")
    log.info("="*60)

    for cohort_key, (df, time_col, event_col, ep_label) in mrna_dfs.items():
        cohort_label = DATASETS[cohort_key]["label"]
        log.info(f"\n  Cohort: {cohort_label}")

        for split in SPLIT_NAMES:
            col = f"ADHESION_OR_{split}"
            if col not in df.columns:
                log.warning(f"    {col} not found")
                continue

            flags     = pd.to_numeric(df[col], errors="coerce")
            df_dysreg = df[flags == 1].copy()
            df_normal = df[flags == 0].copy()

            figure_id = (f"adhesion_or_{cohort_key.lower()}"
                         f"_{split.lower()}_{ep_label.lower()}")
            title     = (f"Adhesion/Motility OR-combined — Any Gene Dysregulated\n"
                         f"{SPLIT_LABELS[split]} | {cohort_label} ({ep_label})")

            _run_binary_km(
                df_dysreg, df_normal,
                "Any Gene Dysregulated", "None Dysregulated",
                title=title, figure_id=figure_id,
                task="task7", gene="ADHESION_OR",
                cohort=cohort_key, split=split,
                subdir="task7_adhesion_or",
                stats_rows=stats_rows, median_rows=median_rows,
                color_a=COLORS["high"], color_b=COLORS["low"],
                time_col=time_col, event_col=event_col, endpoint_label=ep_label,
            )


# ─────────────────────────────────────────────────────────────────────────────
# TASK PFS TCGA — All androgen + adhesion mRNA genes vs PFS in TCGA-PRAD
# Mirrors task3_pfs but covers all expression genes (13 total × 3 splits).
# Loads from cache directly to avoid DFS endpoint filter truncating PFS cohort.
# ─────────────────────────────────────────────────────────────────────────────

def run_task_pfs_tcga(stats_rows, median_rows):
    log.info("\n" + "="*60)
    log.info("TASK PFS TCGA: All androgen + adhesion mRNA genes vs PFS — TCGA-PRAD × 3 splits")
    log.info("="*60)

    cohort_label = DATASETS["TCGA_PRAD"]["label"]
    path = CACHE_DIR / "TCGA_PRAD_flagged_expression.csv"
    if not path.exists():
        log.warning("  TCGA_PRAD_flagged_expression.csv missing — skipping task_pfs_tcga")
        return

    df_raw = pd.read_csv(path, low_memory=False)
    df_pfs = _clean_endpoint(df_raw, "PFS_MONTHS", "PFS_EVENT")
    log.info(f"  TCGA_PRAD PFS: {len(df_pfs)} patients with valid PFS")

    # Androgen genes
    for gene in ANDROGEN_GENES:
        for split in SPLIT_NAMES:
            col = f"{gene}_GROUP_{split}"
            if col not in df_pfs.columns:
                log.warning(f"    {col} not found — skipping")
                continue

            df_high = df_pfs[df_pfs[col] == "High"].copy()
            df_low  = df_pfs[df_pfs[col] == "Low"].copy()

            figure_id = f"{gene.lower()}_tcga_prad_{split.lower()}_pfs"
            title     = (f"{gene} — Androgen Uptake | {SPLIT_LABELS[split]}\n"
                         f"{cohort_label} (PFS)")

            _run_binary_km(
                df_high, df_low,
                f"High {gene}", f"Low {gene}",
                title=title, figure_id=figure_id,
                task="task_pfs_tcga", gene=gene,
                cohort="TCGA_PRAD", split=split,
                subdir="task_pfs_tcga",
                stats_rows=stats_rows, median_rows=median_rows,
                color_a=COLORS["high"], color_b=COLORS["low"],
                time_col="PFS_MONTHS", event_col="PFS_EVENT", endpoint_label="PFS",
            )

    # Adhesion/motility genes
    for gene, cfg in ADHESION_MOTILITY_GENES.items():
        direction = cfg["direction"]
        for split in SPLIT_NAMES:
            col = f"{gene}_GROUP_{split}"
            if col not in df_pfs.columns:
                log.warning(f"    {col} not found — skipping")
                continue

            if direction == "low_is_bad":
                df_bad  = df_pfs[df_pfs[col] == "Low"].copy()
                df_good = df_pfs[df_pfs[col] == "High"].copy()
                label_bad  = f"Low {gene} (epithelial loss)"
                label_good = f"High {gene}"
                color_bad  = COLORS["altered"]
                color_good = COLORS["unaltered"]
                note = "epithelial loss"
            else:
                df_bad  = df_pfs[df_pfs[col] == "High"].copy()
                df_good = df_pfs[df_pfs[col] == "Low"].copy()
                label_bad  = f"High {gene}"
                label_good = f"Low {gene}"
                color_bad  = COLORS["high"]
                color_good = COLORS["low"]
                note = "mesenchymal gain"

            figure_id = f"{gene.lower()}_tcga_prad_{split.lower()}_pfs"
            title     = (f"{gene} ({note}) | {SPLIT_LABELS[split]}\n"
                         f"{cohort_label} (PFS)")

            _run_binary_km(
                df_bad, df_good,
                label_bad, label_good,
                title=title, figure_id=figure_id,
                task="task_pfs_tcga", gene=gene,
                cohort="TCGA_PRAD", split=split,
                subdir="task_pfs_tcga",
                stats_rows=stats_rows, median_rows=median_rows,
                color_a=color_bad, color_b=color_good,
                time_col="PFS_MONTHS", event_col="PFS_EVENT", endpoint_label="PFS",
            )


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    log.info("\n" + "="*70)
    log.info("MODULE 6: KM ANALYSIS — 8 TASKS")
    log.info("="*70)

    stats_rows  = []
    median_rows = []

    log.info("\nLoading 7-cohort combined data (molecular flags)...")
    mol_pool = _load_molecular_pool()

    log.info("\nLoading mRNA cohort data (expression flags)...")
    mrna_dfs = _load_mrna_cohorts()

    run_task1(mol_pool,  stats_rows, median_rows)
    run_task2(mol_pool,  stats_rows, median_rows)
    run_task3(mrna_dfs,      stats_rows, median_rows)
    run_task3b(mol_pool,     stats_rows, median_rows)
    run_task_pfs_tcga(       stats_rows, median_rows)
    run_task4(mrna_dfs,   stats_rows, median_rows)
    run_task6(mrna_dfs,   stats_rows, median_rows)
    run_task6b(mol_pool, stats_rows, median_rows)
    run_task7(mrna_dfs,   stats_rows, median_rows)

    stats_path  = TABLES_DIR / "km_statistics.csv"
    median_path = TABLES_DIR / "median_os_table.csv"
    pd.DataFrame(stats_rows).to_csv(stats_path,  index=False)
    pd.DataFrame(median_rows).to_csv(median_path, index=False)
    log.info(f"\nSaved: {stats_path}  ({len(stats_rows)} rows)")
    log.info(f"Saved: {median_path} ({len(median_rows)} rows)")

    log.info("\n" + "="*70)
    log.info("MODULE 6 COMPLETE")
    log.info("="*70)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--refresh", action="store_true")
    args = p.parse_args()
    main()
    log.info("KM analysis complete.")
