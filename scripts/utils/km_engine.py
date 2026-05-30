"""
utils/km_engine.py
─────────────────────────────────────────────────────────────────────────────
Core Kaplan-Meier engine — backed by lifelines for all statistical computation.

Provides:
  - KaplanMeierCurve   : thin wrapper around lifelines.KaplanMeierFitter
  - logrank_test        : two-group log-rank test (lifelines)
  - multivariate_logrank: k-group log-rank test (lifelines)
  - plot_two_group_km   : plots two KM curves on matplotlib axes
  - plot_multigroup_km  : plots multiple KM curves
  - save_figure         : saves figure in all configured formats
─────────────────────────────────────────────────────────────────────────────
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test as _ll_logrank
from lifelines.statistics import multivariate_logrank_test as _ll_multivariate_logrank

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import STATS, PLOT_STYLE, COLORS, FIGURES_DIR
from utils.logger import get_logger

log = get_logger("km_engine")


# ─────────────────────────────────────────────────────────────────────────────
# MATPLOTLIB STYLE
# ─────────────────────────────────────────────────────────────────────────────

def apply_plot_style():
    plt.rcParams.update({
        "font.family":       PLOT_STYLE.get("font_family", "DejaVu Sans"),
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.titlesize":    PLOT_STYLE["title_fontsize"],
        "axes.labelsize":    PLOT_STYLE["label_fontsize"],
        "xtick.labelsize":   PLOT_STYLE["tick_fontsize"],
        "ytick.labelsize":   PLOT_STYLE["tick_fontsize"],
        "legend.fontsize":   PLOT_STYLE["legend_fontsize"],
        "figure.dpi":        150,
    })


# ─────────────────────────────────────────────────────────────────────────────
# KAPLAN-MEIER CURVE CLASS
# Thin wrapper around lifelines.KaplanMeierFitter that exposes the same
# interface used throughout module6 (timeline, survival, ci_lower, ci_upper,
# median_survival, n_total, n_events, label, plot()).
# ─────────────────────────────────────────────────────────────────────────────

class KaplanMeierCurve:
    """
    Wraps lifelines.KaplanMeierFitter.
    CI: Greenwood formula with log(-log) transformation (lifelines default).
    Median CI: Brookmeyer-Crowley (lifelines default).
    """

    def __init__(self, times, events, label=""):
        self.label    = label
        self.n_total  = len(times)
        self.n_events = int(np.sum(events))

        kmf = KaplanMeierFitter(alpha=0.05)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            kmf.fit(
                np.asarray(times,  dtype=float),
                np.asarray(events, dtype=float),
                label=label,
            )
        self._kmf = kmf

        # Expose arrays in the same shape as before for plot()
        self.timeline = np.concatenate([[0.0], kmf.survival_function_.index.values])
        self.survival = np.concatenate([[1.0], kmf.survival_function_.values.flatten()])
        self.ci_lower = np.concatenate([[1.0], kmf.confidence_interval_.iloc[:, 0].values])
        self.ci_upper = np.concatenate([[1.0], kmf.confidence_interval_.iloc[:, 1].values])

        self.median_survival  = float(kmf.median_survival_time_)
        # Brookmeyer-Crowley CI for median via lifelines percentile method
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                self.median_ci_lower = float(kmf.percentile(0.5,
                    interpolation="midpoint"))
            except Exception:
                self.median_ci_lower = np.nan
        # lifelines does not expose median CI bounds directly;
        # approximate via the CI curve crossing 0.5 (same as Brookmeyer-Crowley)
        ci_lo_arr = self.ci_lower
        ci_hi_arr = self.ci_upper
        idx_lo = np.searchsorted(-ci_hi_arr, -0.5)
        idx_hi = np.searchsorted(-ci_lo_arr, -0.5)
        self.median_ci_lower = self.timeline[idx_lo] if idx_lo < len(self.timeline) else np.nan
        self.median_ci_upper = self.timeline[idx_hi] if idx_hi < len(self.timeline) else np.nan

    def plot(self, ax, color, linewidth=2.0, ci_alpha=0.15, step_where="post"):
        """Draws step function + CI shading on axes."""
        ax.step(self.timeline, self.survival,
                where=step_where, color=color,
                linewidth=linewidth, label=self.label)
        ax.fill_between(self.timeline, self.ci_lower, self.ci_upper,
                        step=step_where, color=color, alpha=ci_alpha)


# ─────────────────────────────────────────────────────────────────────────────
# LOG-RANK TESTS
# ─────────────────────────────────────────────────────────────────────────────

def logrank_test(times_a, events_a, times_b, events_b):
    """
    Two-sample log-rank test via lifelines.
    Returns p-value (float) or None if test cannot be performed.
    """
    t_a = np.asarray(times_a,  dtype=float)
    e_a = np.asarray(events_a, dtype=float)
    t_b = np.asarray(times_b,  dtype=float)
    e_b = np.asarray(events_b, dtype=float)

    if len(t_a) < 2 or len(t_b) < 2:
        return None
    if e_a.sum() + e_b.sum() == 0:
        return None

    try:
        result = _ll_logrank(t_a, t_b, e_a, e_b)
        return float(result.p_value)
    except Exception as exc:
        log.debug(f"logrank_test failed: {exc}")
        return None


def multivariate_logrank_test(df, group_col, time_col, event_col):
    """
    K-sample log-rank test via lifelines.
    Returns overall p-value or None.
    """
    groups = df[group_col].dropna().unique()
    if len(groups) < 2:
        return None

    valid = df[[group_col, time_col, event_col]].dropna()
    valid = valid[valid[time_col] > 0]

    if len(valid) < 10:
        return None

    try:
        result = _ll_multivariate_logrank(
            valid[time_col], valid[group_col], valid[event_col]
        )
        return float(result.p_value)
    except Exception as exc:
        log.debug(f"multivariate_logrank_test failed: {exc}")
        return None


# ─────────────────────────────────────────────────────────────────────────────
# VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

def check_group_size(df, label="group"):
    n     = len(df) if hasattr(df, '__len__') else 0
    min_n = STATS["min_group_size"]
    if n < min_n:
        log.warning(f"  {label}: n={n} is below minimum ({min_n}). "
                    f"KM curve may be unreliable.")
    return n >= min_n


def validate_and_clean(df, time_col, event_col, label=""):
    """
    Validates and returns a clean copy of the DataFrame suitable for KM.
    Returns (clean_df, is_valid).
    """
    if time_col not in df.columns or event_col not in df.columns:
        log.warning(f"  {label}: missing columns {time_col}/{event_col}")
        return pd.DataFrame(), False

    clean = df[[time_col, event_col]].copy()
    clean[time_col]  = pd.to_numeric(clean[time_col],  errors="coerce")
    clean[event_col] = pd.to_numeric(clean[event_col], errors="coerce")
    clean = clean.dropna()
    clean = clean[clean[time_col] > 0]

    if len(clean) < 2:
        log.warning(f"  {label}: only {len(clean)} valid rows")
        return clean, False

    n_events = (clean[event_col] == 1).sum()
    if n_events == 0:
        log.warning(f"  {label}: 0 events — cannot fit KM")
        return clean, False

    return clean, True


def fit_km(df, time_col, event_col, label=""):
    """
    Fits KaplanMeierCurve. Returns curve object or None.
    """
    clean, ok = validate_and_clean(df, time_col, event_col, label)
    if not ok:
        return None
    return KaplanMeierCurve(
        clean[time_col].values,
        clean[event_col].values,
        label=label,
    )


# ─────────────────────────────────────────────────────────────────────────────
# P-VALUE FORMATTING
# ─────────────────────────────────────────────────────────────────────────────

def format_pvalue(p):
    if p is None:
        return "p = N/A"
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}"


def annotate_pvalue(ax, pval, x=0.05, y=0.08):
    txt   = format_pvalue(pval)
    sig   = pval is not None and pval < STATS["significance_threshold"]
    color = "#C62828" if sig else "#424242"
    ax.text(x, y, txt,
            transform=ax.transAxes,
            fontsize=14, fontstyle="italic", color=color,
            bbox=dict(boxstyle="round,pad=0.3",
                      facecolor="white", alpha=0.85,
                      edgecolor=color, linewidth=0.8))


# ─────────────────────────────────────────────────────────────────────────────
# AT-RISK TABLE
# ─────────────────────────────────────────────────────────────────────────────

def add_at_risk_table(ax, curves_and_dfs, time_col, colors, time_points=None):
    """
    Adds a simple at-risk count row below the KM plot.
    curves_and_dfs : list of (KaplanMeierCurve, DataFrame, time_col, event_col, color)
    """
    if time_points is None:
        all_times = []
        for _, df, tc, _ in curves_and_dfs:
            if tc in df.columns:
                all_times.extend(df[tc].dropna().values)
        if not all_times:
            return
        max_t       = np.percentile(all_times, 95)
        time_points = np.linspace(0, max_t, 6).astype(int)

    y_base = -0.22
    ax.text(-0.02, y_base, "At risk:",
            transform=ax.transAxes, fontsize=7, ha="right", va="center")

    for row_idx, (curve, df, tc, ec, color) in enumerate(curves_and_dfs):
        y   = y_base - row_idx * 0.07
        lbl = curve.label.split("(")[0].strip() if curve else ""
        ax.text(-0.02, y, lbl, transform=ax.transAxes,
                fontsize=6, ha="right", va="center", color=color)
        if tc in df.columns:
            t_arr = pd.to_numeric(df[tc], errors="coerce").values
            for t in time_points:
                n = int(np.sum(t_arr >= t))
                ax.text(t / ax.get_xlim()[1] if ax.get_xlim()[1] > 0 else 0.5,
                        y, str(n), transform=ax.transAxes,
                        fontsize=6, ha="center", va="center", color=color)


# ─────────────────────────────────────────────────────────────────────────────
# PLOTTING FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def _style_axes(ax, title, xlabel, ylabel):
    ax.set_title(title, fontsize=PLOT_STYLE["title_fontsize"],
                 fontweight="bold", pad=8)
    ax.set_xlabel(xlabel, fontsize=PLOT_STYLE["label_fontsize"])
    ax.set_ylabel(ylabel, fontsize=PLOT_STYLE["label_fontsize"])
    ax.set_ylim(-0.05, 1.05)
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{x:.0%}"))
    ax.legend(fontsize=PLOT_STYLE["legend_fontsize"],
              loc="upper right", framealpha=0.8)


def plot_two_group_km(ax, df_high, df_low,
                      time_col, event_col,
                      label_high, label_low,
                      color_high=None, color_low=None,
                      title="", pval=None,
                      show_at_risk=False,
                      ylabel="Survival Probability"):
    """
    Plots two KM curves on a given axes object.
    Returns (curve_high, curve_low).
    """
    color_high = color_high or COLORS["high"]
    color_low  = color_low  or COLORS["low"]
    lw         = PLOT_STYLE["linewidth"]
    ci_alpha   = PLOT_STYLE["ci_alpha"]

    n_high = len(df_high.dropna(subset=[time_col, event_col]))
    n_low  = len(df_low.dropna(subset=[time_col, event_col]))

    km_high = fit_km(df_high, time_col, event_col, label=f"{label_high} (n={n_high})")
    km_low  = fit_km(df_low,  time_col, event_col, label=f"{label_low} (n={n_low})")

    if km_high is None and km_low is None:
        ax.text(0.5, 0.5, "Insufficient data",
                ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_title(title)
        return None, None

    if km_high:
        km_high.plot(ax, color=color_high, linewidth=lw, ci_alpha=ci_alpha)
    if km_low:
        km_low.plot(ax,  color=color_low,  linewidth=lw, ci_alpha=ci_alpha)

    if pval is None:
        clean_h, ok_h = validate_and_clean(df_high, time_col, event_col, "high")
        clean_l, ok_l = validate_and_clean(df_low,  time_col, event_col, "low")
        if ok_h and ok_l:
            pval = logrank_test(
                clean_h[time_col], clean_h[event_col],
                clean_l[time_col], clean_l[event_col])

    annotate_pvalue(ax, pval)
    _style_axes(ax, title, "Time (Months)", ylabel)

    for km, color in [(km_high, color_high), (km_low, color_low)]:
        if km and not np.isnan(km.median_survival):
            ax.axvline(km.median_survival, color=color,
                       linestyle="--", linewidth=0.8, alpha=0.6)

    return km_high, km_low


def plot_multigroup_km(ax, group_dfs, time_col, event_col,
                       colors, title="", pval=None,
                       show_at_risk=False,
                       ylabel="Survival Probability"):
    """
    Plots multiple KM curves on one axes.
    group_dfs : list of (label, DataFrame) tuples
    """
    lw       = PLOT_STYLE["linewidth"]
    ci_alpha = PLOT_STYLE["ci_alpha"]
    curves   = []

    for (label, df), color in zip(group_dfs, colors):
        n  = len(df.dropna(subset=[time_col, event_col]))
        km = fit_km(df, time_col, event_col, label=f"{label} (n={n})")
        if km:
            km.plot(ax, color=color, linewidth=lw, ci_alpha=ci_alpha)
            curves.append(km)

    annotate_pvalue(ax, pval)
    _style_axes(ax, title, "Time (Months)", ylabel)
    return curves


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE SAVING
# ─────────────────────────────────────────────────────────────────────────────

def save_figure(fig, filename_stem, subdir):
    out_dir = FIGURES_DIR / subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    for fmt in PLOT_STYLE["formats"]:
        out_path = out_dir / f"{filename_stem}.{fmt}"
        try:
            fig.savefig(out_path, dpi=PLOT_STYLE["figure_dpi"],
                        bbox_inches="tight", format=fmt)
            log.info(f"  Saved: {out_path.name}")
        except Exception as e:
            log.error(f"  Failed to save {out_path.name}: {e}")

    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# ALIASES
# ─────────────────────────────────────────────────────────────────────────────

def run_logrank(df_a, df_b, time_col, event_col):
    clean_a, ok_a = validate_and_clean(df_a, time_col, event_col, "group_a")
    clean_b, ok_b = validate_and_clean(df_b, time_col, event_col, "group_b")
    if not ok_a or not ok_b:
        return None
    return logrank_test(
        clean_a[time_col].values, clean_a[event_col].values,
        clean_b[time_col].values, clean_b[event_col].values,
    )


def run_multivariate_logrank(df, group_col, time_col, event_col):
    return multivariate_logrank_test(df, group_col, time_col, event_col)
