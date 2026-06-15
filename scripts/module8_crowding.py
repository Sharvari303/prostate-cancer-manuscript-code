"""
module8_crowding.py — Kaplan-Meier + Cox analyses for the crowding /
mechanobiology axes (Axes 1–6) anchoring the Fig 5 crowding result.

Scope (decided 2026-06-14, see config.CROWDING_AXES header):
  - mRNA EXPRESSION only (RNA-seq z-scores). PTEN is the ONLY DNA variable,
    used purely as the deep-deletion stratifier (PTEN_DEEPDEL, CNA == -2).
  - Cohort: TCGA-PRAD only. Endpoints: DFS primary + PFS secondary. NO OS
    (TCGA-PRAD is effectively OS-event-free: ~10/494 = 2%).
  - Stratified evidence = Cox  gene_z * PTEN_DEEPDEL  interaction term (uses all
    patients). PTEN-stratified KM curves are produced for figures but are
    DESCRIPTIVE; arms below STATS.min_group_size are flagged, not trusted.
  - Splits: median default; quartile for ANXA1 and the PIP2 composite.
  - FDR: Benjamini-Hochberg across Axis 6 gene p-values; q < 0.10 (exploratory).

Reads:  data/cache/TCGA_PRAD_flagged_expression.csv  (from module5; must contain
        {GENE}_ZSCORE, {GENE}_GROUP_*, PIP2_TRAFFICKING_*, PTEN_DEEPDEL)
Writes: outputs/tables/crowding_km_statistics.csv
        outputs/tables/crowding_cox_interaction.csv
        outputs/figures/task8_crowding_axis{1..6}/*.pdf|svg
        outputs/figures/task8_crowding_composite/*.pdf|svg
        outputs/figures/task8_crowding_forest/*.pdf|svg
"""
import sys
import warnings
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    CROWDING_AXES, COMPOSITE_SCORES,
    CROWDING_COHORT_KEY, CROWDING_ENDPOINTS,
    CROWDING_FDR_Q, CROWDING_AXIS6_GENES,
    DATASETS, ENDPOINTS,
    CACHE_DIR, TABLES_DIR, STATS, COLORS, PLOT_STYLE,
)
from utils.km_engine import (
    apply_plot_style, save_figure,
    plot_two_group_km, fit_km, validate_and_clean,
    logrank_test, format_pvalue,
)
from utils.logger import get_logger

log = get_logger("module8_crowding")

MIN_GROUP  = STATS["min_group_size"]    # 20
MIN_EVENTS = STATS["min_cox_events"]    # 10

# Endpoint column names (TCGA-PRAD has both)
ENDPOINT_COLS = {
    "DFS": ("DFS_MONTHS", "DFS_EVENT", "Disease-Free Survival"),
    "PFS": ("PFS_MONTHS", "PFS_EVENT", "Progression-Free Survival"),
    "OS":  ("OS_MONTHS",  "OS_EVENT",  "Overall Survival"),
}


# ─────────────────────────────────────────────────────────────────────────────
# DATA
# ─────────────────────────────────────────────────────────────────────────────

def _load_crowding_frame():
    """Load the TCGA-PRAD expression frame with crowding columns + PTEN_DEEPDEL."""
    path = CACHE_DIR / f"{CROWDING_COHORT_KEY}_flagged_expression.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing. Run module2→module5 (with --refresh) first so the "
            f"crowding mRNA z-scores and PTEN_DEEPDEL are present.")
    df = pd.read_csv(path, low_memory=False)
    if "PTEN_DEEPDEL" not in df.columns:
        log.warning("PTEN_DEEPDEL absent — re-run module4 --refresh then module5 "
                    "--refresh. PTEN-stratified Cox interaction will be skipped.")
    else:
        df["PTEN_DEEPDEL"] = pd.to_numeric(df["PTEN_DEEPDEL"], errors="coerce")
    return df


def _clean(df, time_col, event_col):
    out = df.copy()
    out[time_col]  = pd.to_numeric(out.get(time_col,  pd.Series(dtype=float)), errors="coerce")
    out[event_col] = pd.to_numeric(out.get(event_col, pd.Series(dtype=float)), errors="coerce")
    out = out.dropna(subset=[time_col, event_col])
    out = out[out[time_col] > 0]
    return out


# ─────────────────────────────────────────────────────────────────────────────
# COX (univariate + interaction)
# ─────────────────────────────────────────────────────────────────────────────

def _cox_univariate_z(df, zcol, time_col, event_col):
    """Univariate Cox on continuous z-score. Returns dict or None.

    HR is per +1 z unit. For low_is_bad genes a protective (HR<1) direction is
    expected on the raw z; interpretation is handled by the caller/label.
    """
    try:
        from lifelines import CoxPHFitter
    except ImportError:
        return None
    d = df[[time_col, event_col, zcol]].copy()
    d[zcol] = pd.to_numeric(d[zcol], errors="coerce")
    d = d.dropna()
    d = d[d[time_col] > 0]
    n_ev = int((d[event_col] == 1).sum())
    if n_ev < MIN_EVENTS or len(d) < MIN_GROUP or d[zcol].nunique() < 2:
        return None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cph = CoxPHFitter()
            cph.fit(d, duration_col=time_col, event_col=event_col)
        return {
            "HR":       float(np.exp(cph.params_[zcol])),
            "CI_lower": float(np.exp(cph.confidence_intervals_.loc[zcol, "95% lower-bound"])),
            "CI_upper": float(np.exp(cph.confidence_intervals_.loc[zcol, "95% upper-bound"])),
            "pval":     float(cph.summary.loc[zcol, "p"]),
            "N":        len(d),
            "events":   n_ev,
        }
    except Exception as e:
        log.debug(f"    univariate Cox failed for {zcol}: {e}")
        return None


def _cox_interaction(df, zcol, time_col, event_col, strat_col="PTEN_DEEPDEL"):
    """Cox model: z + PTEN_DEEPDEL + z:PTEN_DEEPDEL.

    Returns dict with main/PTEN/interaction HR+p, or None. The interaction term
    p-value is the primary stratified evidence ('does PTEN loss amplify risk').
    """
    try:
        from lifelines import CoxPHFitter
    except ImportError:
        return None
    if strat_col not in df.columns:
        return None
    d = df[[time_col, event_col, zcol, strat_col]].copy()
    d[zcol]      = pd.to_numeric(d[zcol], errors="coerce")
    d[strat_col] = pd.to_numeric(d[strat_col], errors="coerce")
    d = d.dropna()
    d = d[d[time_col] > 0]
    # need events and variation in both terms
    n_ev = int((d[event_col] == 1).sum())
    if (n_ev < MIN_EVENTS or len(d) < MIN_GROUP
            or d[zcol].nunique() < 2 or d[strat_col].nunique() < 2):
        return None
    d["_interaction"] = d[zcol] * d[strat_col]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cph = CoxPHFitter()
            cph.fit(d, duration_col=time_col, event_col=event_col)
        def _hr(term):  return float(np.exp(cph.params_[term]))
        def _p(term):   return float(cph.summary.loc[term, "p"])
        return {
            "HR_main":        _hr(zcol),
            "p_main":         _p(zcol),
            "HR_pten":        _hr(strat_col),
            "p_pten":         _p(strat_col),
            "HR_interaction": _hr("_interaction"),
            "p_interaction":  _p("_interaction"),
            "N":              len(d),
            "events":         n_ev,
            "n_pten_loss":    int((d[strat_col] == 1).sum()),
        }
    except Exception as e:
        log.debug(f"    interaction Cox failed for {zcol}: {e}")
        return None


# ─────────────────────────────────────────────────────────────────────────────
# PER-GENE KM
# ─────────────────────────────────────────────────────────────────────────────

def _group_col_for(gene, split):
    return f"{gene}_GROUP_{split.upper()}"


def _km_for_gene(df, gene, cfg, endpoint, km_rows, fig=True):
    """KM (high vs low) for one gene on one endpoint. Direction-aware labels.

    Records a row in km_rows with N/events per arm, log-rank p, and the
    univariate continuous-z Cox HR. Plots the figure when fig=True.
    """
    time_col, event_col, ep_label = ENDPOINT_COLS[endpoint]
    split = cfg["split"]                      # 'median' or 'quartile'
    direction = cfg["direction"]
    gcol = _group_col_for(gene, split)
    zcol = f"{gene}_ZSCORE"

    if gcol not in df.columns:
        log.warning(f"  {gene}: {gcol} missing — skipping")
        return
    d = _clean(df, time_col, event_col)
    if len(d) < MIN_GROUP:
        log.warning(f"  {gene} [{endpoint}]: only {len(d)} with endpoint — skipping")
        return

    # "High"/"Low" arms. Bad arm depends on direction (for color + labeling only).
    df_high = d[d[gcol] == "High"]
    df_low  = d[d[gcol] == "Low"]
    n_high, n_low = len(df_high), len(df_low)
    ev_high = int((df_high[event_col] == 1).sum())
    ev_low  = int((df_low[event_col]  == 1).sum())

    pval = None
    ch, okh = validate_and_clean(df_high, time_col, event_col, "high")
    cl, okl = validate_and_clean(df_low,  time_col, event_col, "low")
    if okh and okl:
        pval = logrank_test(ch[time_col], ch[event_col], cl[time_col], cl[event_col])

    cox = _cox_univariate_z(d, zcol, time_col, event_col) if zcol in d.columns else None

    bad = "High" if direction == "high_is_bad" else "Low"
    km_rows.append({
        "axis":        cfg["axis"],
        "gene":        gene,
        "endpoint":    endpoint,
        "split":       split,
        "direction":   direction,
        "bad_group":   bad,
        "n_high":      n_high, "events_high": ev_high,
        "n_low":       n_low,  "events_low":  ev_low,
        "logrank_p":   pval,
        "cox_z_HR":    cox["HR"]       if cox else np.nan,
        "cox_z_CI_lo": cox["CI_lower"] if cox else np.nan,
        "cox_z_CI_hi": cox["CI_upper"] if cox else np.nan,
        "cox_z_p":     cox["pval"]     if cox else np.nan,
        "cox_z_N":     cox["N"]        if cox else np.nan,
        "cox_z_events":cox["events"]   if cox else np.nan,
        "low_arm_below_min": (n_high < MIN_GROUP) or (n_low < MIN_GROUP),
        "novel":       cfg.get("novel", False),
        "priority":    cfg.get("priority", ""),
    })

    if fig and (okh or okl):
        apply_plot_style()
        f, ax = plt.subplots(figsize=(6, 5))
        split_note = "median" if split == "median" else "top vs bottom quartile"
        title = f"{gene} ({split_note}) — TCGA-PRAD {endpoint}"
        # color the worse-prognosis arm red
        c_high = COLORS["altered"] if direction == "high_is_bad" else COLORS["low"]
        c_low  = COLORS["low"]     if direction == "high_is_bad" else COLORS["altered"]
        plot_two_group_km(ax, df_high, df_low, time_col, event_col,
                          label_high=f"{gene} High", label_low=f"{gene} Low",
                          color_high=c_high, color_low=c_low,
                          title=title, pval=pval, ylabel=f"{ep_label} Probability")
        save_figure(f, f"{gene.lower()}_{endpoint.lower()}", f"task8_crowding_axis{cfg['axis']}")


def _km_pten_stratified(df, gene, cfg, endpoint):
    """Descriptive PTEN-stratified KM figure (4 arms collapsed to 2 panels).

    NOT the primary evidence — that is the Cox interaction. Arms below
    MIN_GROUP are still drawn but the figure is for visualization only.
    """
    if "PTEN_DEEPDEL" not in df.columns:
        return
    time_col, event_col, ep_label = ENDPOINT_COLS[endpoint]
    split = cfg["split"]
    gcol  = _group_col_for(gene, split)
    if gcol not in df.columns:
        return
    d = _clean(df, time_col, event_col)
    d = d[d["PTEN_DEEPDEL"].notna()]
    if len(d) < MIN_GROUP:
        return

    apply_plot_style()
    f, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, pten_val, pten_lab in [(axes[0], 1.0, "PTEN deep-del"),
                                   (axes[1], 0.0, "PTEN intact")]:
        sub = d[d["PTEN_DEEPDEL"] == pten_val]
        df_high = sub[sub[gcol] == "High"]
        df_low  = sub[sub[gcol] == "Low"]
        pval = None
        ch, okh = validate_and_clean(df_high, time_col, event_col, "h")
        cl, okl = validate_and_clean(df_low,  time_col, event_col, "l")
        if okh and okl:
            pval = logrank_test(ch[time_col], ch[event_col], cl[time_col], cl[event_col])
        c_high = COLORS["altered"] if cfg["direction"] == "high_is_bad" else COLORS["low"]
        c_low  = COLORS["low"]     if cfg["direction"] == "high_is_bad" else COLORS["altered"]
        plot_two_group_km(ax, df_high, df_low, time_col, event_col,
                          label_high=f"{gene} High", label_low=f"{gene} Low",
                          color_high=c_high, color_low=c_low,
                          title=f"{pten_lab} (n={len(sub)})", pval=pval,
                          ylabel=f"{ep_label} Probability")
    f.suptitle(f"{gene} stratified by PTEN deep-deletion — TCGA-PRAD {endpoint} "
               f"(descriptive; see Cox interaction)", fontsize=13)
    save_figure(f, f"{gene.lower()}_{endpoint.lower()}_pten_strat",
                f"task8_crowding_axis{cfg['axis']}")


# ─────────────────────────────────────────────────────────────────────────────
# PIP2 COMPOSITE (headline)
# ─────────────────────────────────────────────────────────────────────────────

def _km_composite(df, endpoint, km_rows, inter_rows):
    cfg = COMPOSITE_SCORES["pip2_trafficking_score"]
    score_col = "PIP2_TRAFFICKING_SCORE"
    group_col = "PIP2_TRAFFICKING_GROUP_QUARTILE"
    if score_col not in df.columns or group_col not in df.columns:
        log.warning("  PIP2 composite columns missing — skipping")
        return
    time_col, event_col, ep_label = ENDPOINT_COLS[endpoint]
    d = _clean(df, time_col, event_col)
    df_high = d[d[group_col] == "High"]   # top quartile
    df_low  = d[d[group_col] == "Low"]    # bottom quartile

    pval = None
    ch, okh = validate_and_clean(df_high, time_col, event_col, "h")
    cl, okl = validate_and_clean(df_low,  time_col, event_col, "l")
    if okh and okl:
        pval = logrank_test(ch[time_col], ch[event_col], cl[time_col], cl[event_col])

    cox   = _cox_univariate_z(d, score_col, time_col, event_col)
    inter = _cox_interaction(d, score_col, time_col, event_col)

    km_rows.append({
        "axis": 6, "gene": "PIP2_COMPOSITE", "endpoint": endpoint,
        "split": "quartile", "direction": "high_is_bad", "bad_group": "High",
        "n_high": len(df_high), "events_high": int((df_high[event_col] == 1).sum()),
        "n_low":  len(df_low),  "events_low":  int((df_low[event_col]  == 1).sum()),
        "logrank_p": pval,
        "cox_z_HR":  cox["HR"]       if cox else np.nan,
        "cox_z_CI_lo": cox["CI_lower"] if cox else np.nan,
        "cox_z_CI_hi": cox["CI_upper"] if cox else np.nan,
        "cox_z_p":   cox["pval"]     if cox else np.nan,
        "cox_z_N":   cox["N"]        if cox else np.nan,
        "cox_z_events": cox["events"] if cox else np.nan,
        "low_arm_below_min": (len(df_high) < MIN_GROUP) or (len(df_low) < MIN_GROUP),
        "novel": True, "priority": "highest",
    })
    if inter:
        inter_rows.append({"gene": "PIP2_COMPOSITE", "axis": 6, "endpoint": endpoint, **inter})

    if okh or okl:
        apply_plot_style()
        f, ax = plt.subplots(figsize=(6.5, 5.2))
        plot_two_group_km(ax, df_high, df_low, time_col, event_col,
                          label_high="PIP2 score High (top quartile)",
                          label_low="PIP2 score Low (bottom quartile)",
                          color_high=COLORS["altered"], color_low=COLORS["low"],
                          title=f"PIP2 Trafficking Score — TCGA-PRAD {endpoint}\n(Axis 6 headline)",
                          pval=pval, ylabel=f"{ep_label} Probability")
        save_figure(f, f"pip2_composite_{endpoint.lower()}", "task8_crowding_composite")


# ─────────────────────────────────────────────────────────────────────────────
# FOREST PLOT
# ─────────────────────────────────────────────────────────────────────────────

def _forest_plot(km_df, endpoint):
    """HR (continuous-z univariate Cox) ± 95% CI for all genes, one endpoint."""
    sub = km_df[(km_df["endpoint"] == endpoint) & km_df["cox_z_HR"].notna()].copy()
    if sub.empty:
        return
    sub = sub.sort_values(["axis", "gene"])
    apply_plot_style()
    h = max(4, 0.32 * len(sub))
    f, ax = plt.subplots(figsize=(7, h))
    y = np.arange(len(sub))[::-1]
    for yi, (_, r) in zip(y, sub.iterrows()):
        sig = (r["cox_z_p"] is not None) and (r["cox_z_p"] < STATS["significance_threshold"])
        color = "#C62828" if sig else "#555555"
        ax.plot([r["cox_z_CI_lo"], r["cox_z_CI_hi"]], [yi, yi], color=color, lw=1.5)
        ax.plot(r["cox_z_HR"], yi, "o", color=color, ms=5)
    ax.axvline(1.0, color="grey", ls="--", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"A{int(r['axis'])} {r['gene']}" for _, r in sub.iterrows()],
                       fontsize=8)
    ax.set_xlabel("HR per +1 z (univariate Cox)")
    ax.set_xscale("log")
    ax.set_title(f"Crowding axes — {endpoint} hazard ratios (TCGA-PRAD)", fontsize=12)
    save_figure(f, f"forest_{endpoint.lower()}", "task8_crowding_forest")


# ─────────────────────────────────────────────────────────────────────────────
# FDR CORRECTION (Benjamini-Hochberg, per test family)
#
# The crowding analysis is an exploratory screen, so every p-value is corrected
# for multiple testing within a coherent test family (the q-value denominator =
# number of tests in its family). Pooling unrelated families would over-penalise.
#   Family 1  KM, median split     : per-gene high-vs-low log-rank, all genes × {DFS,PFS}
#   Family 2  KM, quartile split    : ANXA1 + PIP2 composite quartile KM × {DFS,PFS}
#   Family 3  gene × PTEN interaction: Cox interaction term, all genes × {DFS,PFS}
# DFS and PFS are pooled within a family (conservative).
# ─────────────────────────────────────────────────────────────────────────────

def _bh_within(df, pcol, qcol, family_mask):
    """BH-correct pcol within rows where family_mask is True; write q to qcol."""
    try:
        from statsmodels.stats.multitest import multipletests
    except ImportError:
        log.warning("statsmodels not installed — q-values left NaN "
                    "(pip install statsmodels)")
        return df
    sub = df[family_mask & df[pcol].notna()]
    if len(sub) >= 1:
        q = multipletests(sub[pcol].values, method="fdr_bh")[1]
        df.loc[sub.index, qcol] = q
    return df


def _apply_fdr(km_df, inter_df):
    """Add q_logrank_BH / q_interaction_BH + fdr_significant to both tables.

    Benjamini-Hochberg FDR is applied within each test family (KM-median,
    KM-quartile, gene x PTEN interaction) so the q-value denominator is the
    number of tests answering the same question; DFS and PFS are pooled within
    a family. q < CROWDING_FDR_Q (0.10) flags an exploratory discovery.
    """
    km_df["q_logrank_BH"] = pd.NA
    km_df = _bh_within(km_df, "logrank_p", "q_logrank_BH", km_df["split"] == "median")    # Family 1
    km_df = _bh_within(km_df, "logrank_p", "q_logrank_BH", km_df["split"] == "quartile")  # Family 2
    km_df["fdr_significant"] = pd.to_numeric(km_df["q_logrank_BH"], errors="coerce") < CROWDING_FDR_Q

    if not inter_df.empty:
        inter_df["q_interaction_BH"] = pd.NA
        inter_df = _bh_within(inter_df, "p_interaction", "q_interaction_BH",
                              inter_df["p_interaction"].notna())                          # Family 3
        inter_df["fdr_significant"] = pd.to_numeric(
            inter_df["q_interaction_BH"], errors="coerce") < CROWDING_FDR_Q
    return km_df, inter_df


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    log.info("\n" + "=" * 70)
    log.info("MODULE 8: CROWDING / MECHANOBIOLOGY AXES (1–6) — TCGA-PRAD")
    log.info("=" * 70)

    df = _load_crowding_frame()
    log.info(f"Loaded {CROWDING_COHORT_KEY}: {len(df)} patients")
    if "PTEN_DEEPDEL" in df.columns:
        log.info(f"  PTEN deep-del: {int((df['PTEN_DEEPDEL'] == 1).sum())}, "
                 f"intact: {int((df['PTEN_DEEPDEL'] == 0).sum())}")

    km_rows, inter_rows = [], []

    for endpoint in CROWDING_ENDPOINTS:        # DFS, PFS
        log.info(f"\n--- Endpoint: {endpoint} ---")
        for gene, cfg in CROWDING_AXES.items():
            log.info(f"  {gene} (axis {cfg['axis']})")
            _km_for_gene(df, gene, cfg, endpoint, km_rows, fig=True)
            # Cox interaction (primary stratified evidence)
            if cfg["pten_stratify"] and f"{gene}_ZSCORE" in df.columns:
                d = _clean(df, *ENDPOINT_COLS[endpoint][:2])
                inter = _cox_interaction(d, f"{gene}_ZSCORE", *ENDPOINT_COLS[endpoint][:2])
                if inter:
                    inter_rows.append({"gene": gene, "axis": cfg["axis"],
                                       "endpoint": endpoint, **inter})
                # descriptive stratified figure
                _km_pten_stratified(df, gene, cfg, endpoint)
        # composite
        _km_composite(df, endpoint, km_rows, inter_rows)

    km_df = pd.DataFrame(km_rows)
    inter_df = pd.DataFrame(inter_rows)

    # BH-FDR within test families (KM-median, KM-quartile, interaction)
    km_df, inter_df = _apply_fdr(km_df, inter_df)

    # Save tables
    km_path = TABLES_DIR / "crowding_km_statistics.csv"
    km_df.to_csv(km_path, index=False)
    log.info(f"\nSaved: {km_path}  ({len(km_df)} rows)")

    inter_path = TABLES_DIR / "crowding_cox_interaction.csv"
    inter_df.to_csv(inter_path, index=False)
    log.info(f"Saved: {inter_path}  ({len(inter_df)} rows)")

    # Forest plots
    for endpoint in CROWDING_ENDPOINTS:
        _forest_plot(km_df, endpoint)

    # Headline summary — FDR survivors (q < threshold) across both families
    log.info("\n--- FDR survivors (q < %.2f) ---" % CROWDING_FDR_Q)
    km_hits = km_df[km_df.get("fdr_significant", False) == True]
    log.info("KM (per-gene high vs low):")
    if km_hits.empty:
        log.info("  none")
    else:
        for _, r in km_hits.iterrows():
            log.info(f"  {r['gene']:<14} {r['endpoint']} {r['split']:<8} "
                     f"logrank p={r['logrank_p']:.4g} q={float(r['q_logrank_BH']):.4g}")
    log.info("Interaction (gene × PTEN):")
    ix_hits = inter_df[inter_df.get("fdr_significant", False) == True] if not inter_df.empty else inter_df
    if inter_df.empty or ix_hits.empty:
        log.info("  none")
    else:
        for _, r in ix_hits.iterrows():
            log.info(f"  {r['gene']:<14} {r['endpoint']} HR={r['HR_interaction']:.2f} "
                     f"p={r['p_interaction']:.4g} q={float(r['q_interaction_BH']):.4g}")

    log.info("\n" + "=" * 70)
    log.info("MODULE 8 COMPLETE")
    log.info("=" * 70)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--refresh", action="store_true")
    args = p.parse_args()
    main()
    log.info("Crowding analysis complete.")
