import os
import math
from datetime import datetime
from functools import lru_cache

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, f_oneway, levene, bartlett, mannwhitneyu
from scipy.optimize import curve_fit

#####################
# --- Run-directory layout  ---
# Each simulation run writes its output to a directory whose name is
#   f'{RUN_DIR_PREFIX}{Run_ID + RUN_DIR_ID_OFFSET}'
# located under RUN_BASE. These are the values for this project's runs; anyone
# running the script from a fresh checkout can override them via environment
# variables instead of editing the source:
#   export PCA_RUN_BASE=/path/to/my/runs
#   export PCA_RUN_DIR_PREFIX=my_run_        # e.g. dirs my_run_0, my_run_1, ...
#   export PCA_RUN_DIR_ID_OFFSET=0           # 0 if Run_ID already matches dir index
#####################


# ----------- CONFIGURATION -----------
MINUTES_PER_MONTH = 43200  # 30 days * 24h * 60min


RUN_BASE          = os.environ.get('PCA_RUN_BASE', '/ocean/projects/mcb200052p/sskemkar') ##update this to your run base path
RUN_DIR_PREFIX    = os.environ.get('PCA_RUN_DIR_PREFIX', 'PCa_ABM_dp_') ##update this to your run directory prefix
RUN_DIR_ID_OFFSET = int(os.environ.get('PCA_RUN_DIR_ID_OFFSET', '-1')) ##update this to your run directory numbering offset 

# Shared plot styling (used by every plotting function)
LABEL_SIZE = 20
TITLE_SIZE = 22
TICK_SIZE  = 18
# Muted, high-contrast palette suitable for manuscripts
COLORS = ['#1b2a49', '#c0392b', '#2e7d32', '#6a0dad', '#b8860b', '#1a6b8a',
          '#7f3f00', '#004d40', '#880e4f', '#37474f']

# Update metric labels to match your analysis_over_time.csv
METRIC_LABELS = {
    'Total_Cells': 'Total Cells',
    'Alive_PTEN_Deleted': 'PTEN R Cells',
    'Alive_PTEN_Normal': 'PTEN S Cells',
    'PTEN_Ratio': 'PTEN R/PTEN S Ratio',
    'Testosterone_Depletion_Radius': 'Testosterone Depletion Radius',
    'C_avg_all': 'Avg Clustering Index (All Cells)',
    'C_avg_PTEN_normal': 'Avg Clustering Index (PTEN Normal)',
    'C_avg_PTEN_deleted': 'Avg Clustering Index (PTEN Deleted)',
    'Front_Radius_All': 'Front Radius (All Cells)',
    'Front_Radius_PTEN_Normal': 'Front Radius (PTEN Normal)',
    'Front_Radius_PTEN_Deleted': 'Front Radius (PTEN Deleted)',
    'Front_Speed_All': 'Front Speed (All Cells)',
    'Front_Speed_PTEN_Normal': 'Front Speed (PTEN Normal)',
    'Front_Speed_PTEN_Deleted': 'Front Speed (PTEN Deleted)'
}


def load_masterlist(masterlist_path):
    """Load masterlist CSV file"""
    return pd.read_csv(masterlist_path)


def run_id_to_dir(run_id, prefix=RUN_DIR_PREFIX, id_offset=RUN_DIR_ID_OFFSET):
    """Map a masterlist Run_ID to its simulation directory name.

    Each run's outputs live in a directory named f'{prefix}{Run_ID + id_offset}'.
    The v2 masterlist is 1-based while the directories are 0-based, so the
    default id_offset is -1 (PCa_ABM_dp_0 == Run_ID 1).

    When someone else runs this from GitHub, their batch job may name the run
    directories differently. Override the naming without editing this function
    via the PCA_RUN_DIR_PREFIX / PCA_RUN_DIR_ID_OFFSET environment variables
    (or by passing prefix/id_offset explicitly). The mapping convention lives
    here alone.
    """
    return f'{prefix}{int(run_id) + id_offset}'


@lru_cache(maxsize=None)
def _read_run_csv(run_base_dir, run_id, analysis_rel_path):
    """Read one run's analysis CSV, cached so repeated aggregations reuse it."""
    path = os.path.join(run_base_dir, run_id_to_dir(run_id), analysis_rel_path)
    if not os.path.exists(path):
        print(f"Missing: {path}")
        return None
    return pd.read_csv(path)


def filter_masterlist(masterlist_df, constant_conditions):
    """Filter masterlist by constant conditions"""
    df = masterlist_df.copy()
    for key, value in constant_conditions.items():
        if key in df.columns:
            try:
                # Try numeric comparison first (handles int/float column vs int/float value)
                num_value = float(value)
                df = df[pd.to_numeric(df[key], errors='coerce') == num_value]
            except (ValueError, TypeError):
                # Fall back to string comparison for non-numeric values
                df = df[df[key].astype(str) == str(value)]
        else:
            print(f"Warning: Column '{key}' not found in masterlist. Available columns: {df.columns.tolist()}")
    return df


def get_scanning_groups(masterlist_df, scanning_variable):
    """Group runs by scanning variable(s) and return dict: {composition_label: [run_ids]}"""
    # Ensure scanning_variable is a list
    if isinstance(scanning_variable, str):
        scanning_variable = [scanning_variable]

    # Check if scanning variables exist in the dataframe
    for var in scanning_variable:
        if var not in masterlist_df.columns:
            print(f"Warning: Scanning variable '{var}' not found in filtered masterlist.")
            print(f"Available columns: {masterlist_df.columns.tolist()}")
            return {}

    # Group by scanning variable(s). Passing a list always yields tuple keys.
    groups = {}
    for name, group in masterlist_df.groupby(scanning_variable):
        label = ','.join(f"{var}={val}" for var, val in zip(scanning_variable, name))
        groups[label] = group['Run_ID'].astype(str).tolist()

    return groups


def load_run_data(run_base_dir, run_id, analysis_rel_path, max_time_min=None):
    """Load analysis data for a single run (returns a defensive copy)"""
    df = _read_run_csv(run_base_dir, run_id, analysis_rel_path)
    if df is None:
        return None
    if max_time_min is not None:
        df = df[df['time_min'] <= max_time_min]
    return df.reset_index(drop=True)


def aggregate_runs(run_base_dir, run_ids, analysis_rel_path, metric, final_time_months):
    """Aggregate a metric across multiple runs. Returns (times, mean, std, data)."""
    max_time_min = final_time_months * MINUTES_PER_MONTH
    dfs = []
    for run_id in run_ids:
        df = load_run_data(run_base_dir, run_id, analysis_rel_path, max_time_min)
        if df is not None:
            # Compute derived metrics if needed
            if metric == 'PTEN_Ratio':
                df['PTEN_Ratio'] = df['Alive_PTEN_Deleted'] / (df['Alive_PTEN_Normal'] + 1e-8)
            dfs.append(df)

    if not dfs:
        return None

    # Just use time points from first run (assuming all identical)
    times = dfs[0]['time_min'].values
    # Stack data directly without interpolation
    data = np.array([df[metric].values for df in dfs])
    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)
    return times, mean, std, data


def final_timepoint_values(groups, run_base_dir, analysis_rel_path, metric, final_time_months):
    """For each group, return the per-run metric values at the final timepoint.

    Returns dict: {label: np.ndarray of values}. Shared by the plotting and
    statistical routines so runs are aggregated once per (group, metric).
    """
    data_dict = {}
    for label, run_ids in groups.items():
        result = aggregate_runs(run_base_dir, run_ids, analysis_rel_path, metric, final_time_months)
        if result is None:
            continue
        times, _mean, _std, data = result
        idx = np.argmin(np.abs(times / MINUTES_PER_MONTH - final_time_months))
        data_dict[label] = data[:, idx]
    return data_dict


def save_legend(handles, labels, output_dir, timestamp, plot_type):
    """Save a shared legend as its own figure (drawn once per plot-type).

    The scan-group labels are identical across every metric/cohort of a given
    plot-type, so a single legend file documents them all -- keeping it off the
    individual plots so they render at full size.
    """
    if not handles:
        return
    fig_leg = plt.figure(figsize=(4, 0.5 + 0.4 * len(labels)))
    fig_leg.legend(handles, labels, loc='center', fontsize=TICK_SIZE, frameon=False)
    fig_leg.gca().axis('off')
    fname = os.path.join(output_dir, f"legend_{plot_type}_{timestamp}.png")
    fig_leg.savefig(fname, dpi=300, bbox_inches='tight')
    plt.close(fig_leg)
    print(f"Saved legend: {fname}")


def plot_temporal(groups, run_base_dir, analysis_rel_path, metric, constant_conditions, final_time_months, output_dir, timestamp, ylims=None):
    """Plot temporal profile (line plot with error bands).

    Returns (handles, labels) for the shared legend, which is saved separately
    (see save_legend) rather than drawn on the crowded plot itself.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    if not groups:
        print(f"No groups to plot for {metric}")
        plt.close(fig)
        return [], []

    for i, (label, run_ids) in enumerate(groups.items()):
        result = aggregate_runs(run_base_dir, run_ids, analysis_rel_path, metric, final_time_months)
        if result is None:
            continue
        times, mean, std, data = result
        color = COLORS[i % len(COLORS)]
        ax.plot(times / MINUTES_PER_MONTH, mean, label=label, color=color, linewidth=2.5)
        ax.fill_between(times / MINUTES_PER_MONTH, mean - std, mean + std, alpha=0.15, color=color)

    ax.set_xlabel('Time (months)', fontsize=LABEL_SIZE)
    ax.set_ylabel(METRIC_LABELS.get(metric, metric), fontsize=LABEL_SIZE)
    ax.set_title(f'{METRIC_LABELS.get(metric, metric)} vs Time', fontsize=TITLE_SIZE)
    ax.tick_params(axis='both', labelsize=TICK_SIZE, width=1.5, length=6)
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if ylims and metric in ylims:
        ax.set_ylim(ylims[metric])

    # Legend is saved separately (save_legend) so it doesn't shrink the plot.
    handles, labels = ax.get_legend_handles_labels()

    plt.tight_layout()
    fname = f"{metric}__{constant_conditions}_temporal_{timestamp}.png"
    plt.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches='tight')
    plt.close(fig)
    return handles, labels


def plot_violin(groups, run_base_dir, analysis_rel_path, metric, constant_conditions, final_time_months, output_dir, timestamp, ylims=None):
    """Plot distribution at final time point (violin plot)"""
    fig, ax = plt.subplots(figsize=(10, 6))

    if not groups:
        print(f"No groups to plot for {metric}")
        return

    data_dict = final_timepoint_values(groups, run_base_dir, analysis_rel_path, metric, final_time_months)
    violin_labels = [label for label, values in data_dict.items() if len(values) > 0]
    violin_data = [data_dict[label] for label in violin_labels]

    if not violin_data:
        print(f"No data to plot violin for {metric}")
        return

    sns.violinplot(data=violin_data, ax=ax)
    ax.set_xticks(range(len(violin_labels)))
    ax.set_xticklabels(violin_labels, rotation=45, ha='right', fontsize=TICK_SIZE)
    ax.set_xlabel('Tissue Composition', fontsize=LABEL_SIZE)
    ax.set_ylabel(METRIC_LABELS.get(metric, metric), fontsize=LABEL_SIZE)
    ax.set_title(f'{METRIC_LABELS.get(metric, metric)} at {final_time_months} months', fontsize=TITLE_SIZE)
    ax.tick_params(axis='both', labelsize=TICK_SIZE, width=1.5, length=6)
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if ylims and metric in ylims:
        ax.set_ylim(ylims[metric])

    plt.tight_layout()
    fname = f"{metric}_{constant_conditions}_violin_{timestamp}.png"
    plt.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches='tight')
    plt.show()


def plot_metrics(
    masterlist_df,
    run_base_dir,
    analysis_rel_path,
    constant_conditions,
    scanning_variable,
    metrics,
    plot_type='temporal',
    final_time_months=15,
    output_dir='datanalaysis_output',
    ylims=None
):
    """Main function to generate plots based on specified conditions"""
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    print(f"Filtering by constant conditions: {constant_conditions}")
    filtered_df = filter_masterlist(masterlist_df, constant_conditions)

    if filtered_df.empty:
        print("No runs match the specified conditions!")
        return

    print(f"Found {len(filtered_df)} matching runs")
    print(f"Grouping by scanning variables: {scanning_variable}")
    groups = get_scanning_groups(filtered_df, scanning_variable)

    if not groups:
        print("No valid groups found after scanning!")
        return

    print(f"Found {len(groups)} groups: {list(groups.keys())}")

    # Accept a single plot type or a list, so violin + temporal run together.
    plot_types = [plot_type] if isinstance(plot_type, str) else list(plot_type)
    for pt in plot_types:
        if pt not in ('temporal', 'violin'):
            raise ValueError("plot_type must be 'temporal' or 'violin' (or a list of these)")

    all_stats_rows = []
    temporal_legend = ([], [])  # captured once, saved as a shared legend file

    for metric in metrics:
        for pt in plot_types:
            print(f"Plotting {metric} ({pt})")
            if pt == 'temporal':
                handles, labels = plot_temporal(groups, run_base_dir, analysis_rel_path, metric, constant_conditions, final_time_months, output_dir, timestamp, ylims=ylims)
                if handles:
                    temporal_legend = (handles, labels)
            else:  # violin
                plot_violin(groups, run_base_dir, analysis_rel_path, metric, constant_conditions, final_time_months, output_dir, timestamp, ylims=ylims)

        # --- Statistical analysis (metric-level, independent of plot type) ---
        print(f"Statistical analysis for {metric} at {final_time_months} months:")
        rows = statistical_analysis_groups(groups, run_base_dir, analysis_rel_path, metric, final_time_months)
        all_stats_rows.extend(rows)
        if 'violin' in plot_types:
            print(f"Adjacent group comparisons for {metric}:")
            adj_rows = statistical_analysis_adjacent(groups, run_base_dir, analysis_rel_path, metric, final_time_months)
            all_stats_rows.extend(adj_rows)

    # Save the shared temporal legend once (labels are the same for every metric).
    if 'temporal' in plot_types:
        save_legend(*temporal_legend, output_dir, timestamp, 'temporal')

    stats_df = pd.DataFrame(all_stats_rows)
    stats_fname = os.path.join(output_dir, f"stats_{timestamp}.csv")
    stats_df.to_csv(stats_fname, index=False)
    print(f"Saved stats: {stats_fname}")


def _stat_row(metric, **overrides):
    """Build one stats-CSV row with empty defaults, overriding the fields given."""
    row = {'metric': metric, 'group': '', 'n': '', 'mean': '', 'std': '', 'median': '',
           'test': '', 'stat': np.nan, 'p_value': np.nan, 'p_bonferroni': '', 'comparison': ''}
    row.update(overrides)
    return row


def statistical_analysis_groups(groups, run_base_dir, analysis_rel_path, metric, final_time_months):
    """
    For each group, collect metric values at the final timepoint and perform statistical tests.
    Returns a list of result dicts for saving to CSV.
    """
    data_dict = final_timepoint_values(groups, run_base_dir, analysis_rel_path, metric, final_time_months)

    rows = []
    if len(data_dict) < 2:
        print(f"Not enough groups for statistical comparison for {metric}.")
        return rows

    values = list(data_dict.values())
    labels = list(data_dict.keys())

    # Per-group descriptive stats
    for label, vals in data_dict.items():
        rows.append(_stat_row(metric, group=label, n=len(vals),
                              mean=np.mean(vals), std=np.std(vals), median=np.median(vals)))

    if len(values) == 2:
        stat, p = ttest_ind(values[0], values[1])
        comparison = f"{labels[0]} vs {labels[1]}"
        print(f"T-test for {metric}: stat={stat:.3f}, p={p:.3e} ({comparison})")
        rows.append(_stat_row(metric, test='ttest_ind', stat=stat, p_value=p, comparison=comparison))
    else:
        stat, p = f_oneway(*values)
        print(f"ANOVA for {metric}: stat={stat:.3f}, p={p:.3e} (groups: {labels})")
        rows.append(_stat_row(metric, test='one_way_ANOVA', stat=stat, p_value=p, comparison=str(labels)))

    stat_lev, p_lev = levene(*values)
    stat_bart, p_bart = bartlett(*values)
    print(f"Levene's test for {metric}: stat={stat_lev:.3f}, p={p_lev:.3e}")
    print(f"Bartlett's test for {metric}: stat={stat_bart:.3f}, p={p_bart:.3e}")
    rows.append(_stat_row(metric, test='levene', stat=stat_lev, p_value=p_lev, comparison=str(labels)))
    rows.append(_stat_row(metric, test='bartlett', stat=stat_bart, p_value=p_bart, comparison=str(labels)))

    return rows


def _sort_key(label):
    """Extract trailing numeric value from a group label for sorting.
    e.g. 'PTEN_null=15' -> 15.0, 'Cell_cell_adhesion_multiplier=0.1' -> 0.1
    Falls back to alphabetical sort if no numeric value found.
    """
    try:
        return float(label.split('=')[-1])
    except ValueError:
        return label


def _bonferroni(p, n):
    """Bonferroni-corrected p-value, clamped to 1.0."""
    return min(p * n, 1.0)


def statistical_analysis_adjacent(groups, run_base_dir, analysis_rel_path, metric, final_time_months):
    """
    For each pair of adjacent groups (sorted numerically by scanning variable value),
    run Mann-Whitney U (means) and Levene's test (variances) with Bonferroni correction.
    Works for any scanning variable (PTEN_null, adhesion multiplier, uptake rate, etc.).
    Returns a list of result dicts for saving to CSV.
    """
    data_dict = final_timepoint_values(groups, run_base_dir, analysis_rel_path, metric, final_time_months)

    # Sort labels numerically by the scanning variable value
    labels = sorted(data_dict.keys(), key=_sort_key)
    rows = []

    if len(labels) < 2:
        print(f"Not enough groups for adjacent comparison for {metric}.")
        return rows

    n_pairs = len(labels) - 1  # number of adjacent pairs for Bonferroni

    for i in range(n_pairs):
        label_a, label_b = labels[i], labels[i + 1]
        vals_a, vals_b = data_dict[label_a], data_dict[label_b]
        comparison = f"{label_a} vs {label_b}"

        # Mann-Whitney U (mean comparison)
        stat_mw, p_mw = mannwhitneyu(vals_a, vals_b, alternative='two-sided')
        p_mw_corrected = _bonferroni(p_mw, n_pairs)
        print(f"  Mann-Whitney U [{metric}] {comparison}: U={stat_mw:.3f}, p={p_mw:.3e}, p_bonf={p_mw_corrected:.3e}")
        rows.append(_stat_row(metric, test='mannwhitneyu', stat=stat_mw, p_value=p_mw,
                              p_bonferroni=p_mw_corrected, comparison=comparison))

        # Levene's test (variance comparison)
        stat_lev, p_lev = levene(vals_a, vals_b)
        p_lev_corrected = _bonferroni(p_lev, n_pairs)
        print(f"  Levene      [{metric}] {comparison}: stat={stat_lev:.3f}, p={p_lev:.3e}, p_bonf={p_lev_corrected:.3e}")
        rows.append(_stat_row(metric, test='levene_adjacent', stat=stat_lev, p_value=p_lev,
                              p_bonferroni=p_lev_corrected, comparison=comparison))

    return rows


def _sigmoid(x, L, x0, k, b):
    """4-parameter sigmoid: L / (1 + exp(-k*(x - x0))) + b"""
    return L / (1.0 + np.exp(-k * (x - x0))) + b


def plot_clustering_vs_adhesion(
    masterlist_df,
    run_base_dir,
    analysis_rel_path,
    cohorts,
    metric='C_avg_PTEN_deleted',
    androgen_condition='High',
    scenario='Scenario2',
    uptake_rate_multipliers=None,
    final_time_months=15,
    output_dir='datanalaysis_output'
):
    """
    Plot a clustering index metric vs log10 adhesion strength for multiple cohorts.
    Produces one figure per uptake_rate_multiplier value, with a sigmoid fit overlaid
    on each cohort's data points.

    Parameters
    ----------
    uptake_rate_multipliers : list of float, optional
        Values to iterate over. If None, one figure is made with no uptake filter.
    """
    os.makedirs(output_dir, exist_ok=True)

    colors = {'BR': '#1f77b4', 'CTRL': '#ff7f0e', 'TR': '#2ca02c'}
    markers = {'BR': 'o', 'CTRL': 's', 'TR': '^'}
    metric_label = METRIC_LABELS.get(metric, metric)

    uptake_list = uptake_rate_multipliers if uptake_rate_multipliers is not None else [None]

    for uptake in uptake_list:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        fig, ax = plt.subplots(figsize=(12, 7))
        uptake_label = f'uptake={uptake}' if uptake is not None else 'all_uptakes'
        print(f"\n--- {metric} | {uptake_label} ---")

        sigmoid_rows = []

        for cohort in cohorts:
            print(f"  Cohort: {cohort}")

            filter_conds = {
                'Cohort': cohort,
                'Androgen_condition': androgen_condition,
                'Scenario': scenario,
            }
            if uptake is not None:
                filter_conds['Uptake_rate_multiplier'] = uptake

            filtered_df = filter_masterlist(masterlist_df, filter_conds)
            if filtered_df.empty:
                print(f"    No runs found")
                continue

            # Convert adhesion column to numeric once, then iterate its levels.
            filtered_df = filtered_df.assign(
                Cell_cell_adhesion_multiplier=pd.to_numeric(
                    filtered_df['Cell_cell_adhesion_multiplier'], errors='coerce'))
            adhesion_values = sorted(
                filtered_df['Cell_cell_adhesion_multiplier'].dropna().unique())

            clustering_means, clustering_stds, log_adhesion_values = [], [], []

            for adhesion in adhesion_values:
                subset_df = filtered_df[filtered_df['Cell_cell_adhesion_multiplier'] == adhesion]
                if subset_df.empty:
                    continue

                run_ids = subset_df['Run_ID'].astype(str).tolist()
                result = aggregate_runs(run_base_dir, run_ids, analysis_rel_path,
                                        metric, final_time_months)
                if result is None:
                    continue

                times, mean, std, data = result
                final_value = mean[-1]
                final_std = std[-1]
                clustering_means.append(final_value)
                clustering_stds.append(final_std)
                log_adhesion_values.append(math.log10(adhesion))
                print(f"    adhesion={adhesion:.2g} (log10={math.log10(adhesion):.2f}): "
                      f"{metric} = {final_value:.3f} ± {final_std:.3f}")

            if not clustering_means:
                continue

            x = np.array(log_adhesion_values)
            y = np.array(clustering_means)
            color = colors.get(cohort, None)

            # Data points with error bars
            ax.errorbar(x, y, yerr=clustering_stds,
                        marker=markers.get(cohort, 'o'), linestyle='none',
                        linewidth=2.5, markersize=12, label=cohort,
                        color=color, capsize=7, elinewidth=2, alpha=0.9)

            # Sigmoid fit
            try:
                y_range = y.max() - y.min()
                x_mid = x.mean()
                p0 = [y_range, x_mid, 1.0, y.min()]
                popt, pcov = curve_fit(_sigmoid, x, y, p0=p0, maxfev=10000)
                perr = np.sqrt(np.diag(pcov))
                x_fit = np.linspace(x.min(), x.max(), 300)
                y_fit = _sigmoid(x_fit, *popt)
                ax.plot(x_fit, y_fit, linestyle='--', linewidth=2.5, color=color, alpha=0.7)
                L, x0, k, b = popt
                print(f"    Sigmoid fit: L={L:.3f}±{perr[0]:.3f}, x0={x0:.3f}±{perr[1]:.3f}, "
                      f"k={k:.3f}±{perr[2]:.3f}, b={b:.3f}±{perr[3]:.3f}")
                sigmoid_rows.append({
                    'cohort': cohort, 'metric': metric,
                    'uptake_rate_multiplier': uptake,
                    'L': L, 'L_err': perr[0],
                    'x0': x0, 'x0_err': perr[1],
                    'k': k, 'k_err': perr[2],
                    'b': b, 'b_err': perr[3],
                })
            except Exception as e:
                print(f"    Sigmoid fit failed for {cohort}: {e}")

        # Formatting
        ax.set_xlabel('log10(Cell-cell Adhesion Multiplier)', fontsize=20, fontweight='bold')
        ax.set_ylabel(metric_label, fontsize=20, fontweight='bold')
        ax.set_title(f'{metric_label} vs Adhesion Strength\n'
                     f'({androgen_condition} AR, {scenario}, {uptake_label})',
                     fontsize=22, fontweight='bold')
        ax.tick_params(axis='both', labelsize=18, width=1.5, length=6)
        ax.legend(fontsize=16, loc='best', framealpha=0.9)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()

        fname = f"{metric}_vs_adhesion_{uptake_label}_{timestamp}.png"
        plt.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches='tight')
        print(f"  Saved: {os.path.join(output_dir, fname)}")
        plt.show()

        # Save sigmoid parameters to CSV
        if sigmoid_rows:
            sig_df = pd.DataFrame(sigmoid_rows)
            sig_fname = os.path.join(output_dir, f"{metric}_sigmoid_params_{uptake_label}_{timestamp}.csv")
            sig_df.to_csv(sig_fname, index=False)
            print(f"  Saved sigmoid params: {sig_fname}")


# ----------- CONFIGURATION FOR RESULT 1 -----------
if __name__ == '__main__':


    _SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    MASTERLIST  = os.environ.get(
        'PCA_MASTERLIST',
        os.path.join(_SCRIPT_DIR, 'ABMruns_masterlist_prostatecancer_v2.csv')) ##update this to your masterlist path
    ANALYSIS    = 'output/analysis_over_time.csv'
    METRICS     = ['Total_Cells', 'Alive_PTEN_Deleted', 'Alive_PTEN_Normal', 'PTEN_Ratio', 'C_avg_all', 'C_avg_PTEN_deleted']
    FINAL_TIME  = 15
    # One or both of 'violin' / 'temporal'. A list produces both in one run.
    PLOT_TYPE   = ['violin', 'temporal']
    OUTPUT_ROOT = 'results/result_3' ##update this to your output root path

    # --- LOOP_OVER: variable used to split into separate plots (one plot file per value) ---
    # Each entry: (column_name, value) — one plot_metrics call per entry
    # Example: loop over cohorts
    LOOP_OVER = [
        ('Cohort', 'BR'),
        ('Cohort', 'CTRL'),
        ('Cohort', 'TR'),
    ]
    # Example: loop over androgen conditions instead
    # LOOP_OVER = [
    #     ('Androgen_condition', 'High'),
    #     ('Androgen_condition', 'Low'),
    # ]

    # --- CONSTANT_CONDITIONS: fixed for all plots (do NOT include the LOOP_OVER key here) ---
    CONSTANT_CONDITIONS = {
        'Scenario'           : 'Scenario2',
        'Androgen_condition' : 'High',
        'Uptake_rate_multiplier': 1
    }
    # Example: if looping over cohorts and scanning adhesion:
    # CONSTANT_CONDITIONS = {
    #     'Scenario': 'Scenario2',
    #     'Androgen_condition': 'High',
    #     'Uptake_rate_multiplier': '1',
    # }

    # --- SCAN_VAR: variable whose values become separate lines within each plot ---
    SCAN_VAR = ['Cell_cell_adhesion_multiplier']
    # Example for adhesion SA:  SCAN_VAR = ['Cell_cell_adhesion_multiplier']
    # Example for uptake SA:    SCAN_VAR = ['Uptake_rate_multiplier']

    # Load the masterlist once and reuse it everywhere below.
    print(f"Loading masterlist from: {MASTERLIST}")
    masterlist_df = load_masterlist(MASTERLIST)

    # --- compute global y-limits across ALL loop values and ALL scan groups ---
    # Temporal needs the full time-series mean +/- std range; violin only needs
    # the final timepoint. The full-series range contains the final-timepoint
    # values, so when temporal is active we use it for both plot types.
    print("Computing global y-limits...")
    _plot_types = [PLOT_TYPE] if isinstance(PLOT_TYPE, str) else PLOT_TYPE
    use_full_series = 'temporal' in _plot_types
    ylims = {}
    for metric in METRICS:
        all_vals = []
        for loop_col, loop_val in LOOP_OVER:
            conditions = {**CONSTANT_CONDITIONS, loop_col: loop_val}
            filtered_df = filter_masterlist(masterlist_df, conditions)
            groups = get_scanning_groups(filtered_df, SCAN_VAR)
            for label, run_ids in groups.items():
                result = aggregate_runs(RUN_BASE, run_ids, ANALYSIS, metric, FINAL_TIME)
                if result is not None:
                    times, mean, std, data = result
                    if use_full_series:
                        # full time series mean +/- std (covers temporal & violin)
                        all_vals.append(mean + std)
                        all_vals.append(mean - std)
                    else:
                        # violin only: values at the final timepoint
                        idx = np.argmin(np.abs(times / MINUTES_PER_MONTH - FINAL_TIME))
                        all_vals.append(data[:, idx])
        if all_vals:
            all_vals_arr = np.concatenate(all_vals)
            ymin = np.min(all_vals_arr)
            ymax = np.max(all_vals_arr)
            pad  = 0.05 * (ymax - ymin)
            ylims[metric] = (ymin - pad, ymax + pad)

    # --- loop and plot ---
    for loop_col, loop_val in LOOP_OVER:
        print(f"\n=== {loop_col}: {loop_val} ===")
        conditions = {**CONSTANT_CONDITIONS, loop_col: loop_val}
        plot_metrics(
            masterlist_df       = masterlist_df,
            run_base_dir        = RUN_BASE,
            analysis_rel_path   = ANALYSIS,
            constant_conditions = conditions,
            scanning_variable   = SCAN_VAR,
            metrics             = METRICS,
            plot_type           = PLOT_TYPE,
            final_time_months   = FINAL_TIME,
            output_dir          = f'{OUTPUT_ROOT}/{loop_val}',
            ylims               = ylims
        )

    # --- clustering index vs adhesion strength plots ---
    # One figure per (metric x uptake_rate_multiplier); cohorts as separate lines with sigmoid fits
    UPTAKE_RATES = sorted(
        pd.to_numeric(masterlist_df['Uptake_rate_multiplier'], errors='coerce').dropna().unique()
    )
    for clust_metric in ['C_avg_all', 'C_avg_PTEN_deleted']:
        plot_clustering_vs_adhesion(
            masterlist_df            = masterlist_df,
            run_base_dir             = RUN_BASE,
            analysis_rel_path        = ANALYSIS,
            cohorts                  = [loop_val for _, loop_val in LOOP_OVER],
            metric                   = clust_metric,
            androgen_condition       = CONSTANT_CONDITIONS['Androgen_condition'],
            scenario                 = CONSTANT_CONDITIONS['Scenario'],
            uptake_rate_multipliers  = UPTAKE_RATES,
            final_time_months        = FINAL_TIME,
            output_dir               = f'{OUTPUT_ROOT}/clustering_vs_adhesion',
        )
