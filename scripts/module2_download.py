"""
module2_download.py — Downloads all data from cBioPortal for each cohort.

For each cohort, fetches 13 data types and left-merges everything onto
clinical data (clinical is the anchor — patients without survival data
are excluded). Results saved as one master CSV per cohort.

Data fetched per cohort:
  Clinical          survival endpoints (OS, DFS)
  CNA               GISTIC2 discrete scores for all gene axes
  Mutations         somatic mutations (binary flag + mutation type)
  mRNA z-scores     molecular genes, AR activity genes, androgen + adhesion genes
  Structural vars   gene fusions/rearrangements for molecular genes

Mutation wildtype rule: a patient is called wildtype (0) only if confirmed
present in the mutation profile. Patients not in the profile stay NaN.

Output: data/cache/{KEY}_master.csv  (one per cohort)
"""

import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import (
    DATASETS, ENTREZ_IDS,
    ALL_MOLECULAR_GENES, ALL_AR_ACTIVITY_GENES, ALL_EXPRESSION_GENES,
    ALL_INDIVIDUAL_GENES, ALL_ANDROGEN_CNA_GENES, ALL_ADHESION_CNA_GENES,
    ALL_CROWDING_GENES,
    CACHE_DIR,
)
from utils.api_client import (
    get_study_info, fetch_clinical_data,
    fetch_cna_data, fetch_mutation_data, fetch_mrna_data,
    fetch_sv_data,
)
from utils.logger import get_logger

log = get_logger("module2_download")


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _entrez_for(gene_list):
    return [ENTREZ_IDS[g] for g in gene_list if g in ENTREZ_IDS]


def _merge_left_all(base_df, *other_dfs):
    """Left-merges all dfs onto base using patientId."""
    result = base_df.copy()
    result["patientId"] = result["patientId"].astype(str)
    for df in other_dfs:
        if df is not None and len(df) > 0 and "patientId" in df.columns:
            df = df.copy()
            df["patientId"] = df["patientId"].astype(str)
            result = result.merge(df, on="patientId", how="left")
    return result


def _pivot_mutations_to_wide(mut_df, genes):
    """
    Long mutation records → wide binary flags per gene.
    Returns DataFrame with columns: patientId, {GENE}_MUT,
    {GENE}_MUT_TYPE.
    """
    if mut_df.empty:
        log.warning("Empty mutation dataframe — all flags will be 0")
        return pd.DataFrame(columns=["patientId"])

    rows = {}
    for gene in genes:
        gene_muts = mut_df[mut_df["gene"] == gene]
        for patient in gene_muts["patientId"].unique():
            if patient not in rows:
                rows[patient] = {"patientId": patient}
            rows[patient][f"{gene}_MUT"] = 1
            pt_types = gene_muts[
                gene_muts["patientId"] == patient
            ]["mutationType"].values
            rows[patient][f"{gene}_MUT_TYPE"] = (
                pt_types[0] if len(pt_types) else "")

    if not rows:
        return pd.DataFrame(columns=["patientId"])

    wide = pd.DataFrame(list(rows.values()))
    for gene in genes:
        col = f"{gene}_MUT"
        if col not in wide.columns:
            wide[col] = 0
        else:
            wide[col] = wide[col].fillna(0).astype(int)
    return wide


def _log_column_summary(df):
    """Prints completeness for all key columns."""
    key_cols = (
        [c for c in df.columns if "MONTHS" in c or "STATUS" in c] +
        [c for c in df.columns if "_CNA" in c] +
        [c for c in df.columns if "_MUT" in c and "TYPE" not in c] +
        [c for c in df.columns if "_ZSCORE" in c]
    )
    log.info("  Column completeness:")
    for col in key_cols:
        n     = df[col].notna().sum()
        pct   = 100 * n / len(df)
        log.info(f"    {col:<35} {n:>4}/{len(df)} ({pct:.0f}%)")


# ─────────────────────────────────────────────────────────────────────────────
# SINGLE DATASET PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def download_dataset(dataset_key, force_refresh=False):
    """
    Downloads and assembles all data for one dataset.

    Parameters
    ----------
    dataset_key   : key in config.DATASETS e.g. "TCGA_PRAD"
    force_refresh : bypass cache and re-fetch from API

    Returns
    -------
    master DataFrame
    """
    cfg      = DATASETS[dataset_key]
    study_id = cfg["study_id"]
    label    = cfg["label"]

    log.info(f"\n{'='*60}")
    log.info(f"Dataset: {label}  ({study_id})")
    log.info(f"{'='*60}")

    # Check assembled master cache first
    master_cache = CACHE_DIR / f"{dataset_key}_master.csv"
    if master_cache.exists() and not force_refresh:
        log.info(f"Loading master from cache: {master_cache.name}")
        df = pd.read_csv(master_cache)
        log.info(f"  {len(df)} patients, {len(df.columns)} columns")
        return df

    # Verify study exists
    info = get_study_info(study_id)
    log.info(f"Study confirmed: {info['name']} "
             f"({info['allSampleCount']} samples)")

    # Step 1: Clinical
    log.info("Step 1/13 — Clinical data")
    clinical = fetch_clinical_data(study_id, force_refresh)
    clinical["DATASET"] = label

    # Step 2: CNA (optional — some cohorts lack CNA data)
    log.info("Step 2/13 — CNA data (molecular axis)")
    try:
        cna = fetch_cna_data(
            study_id,
            entrez_ids=_entrez_for(ALL_MOLECULAR_GENES),
            gene_names=ALL_MOLECULAR_GENES,
            force_refresh=force_refresh,
            cache_suffix="cna_molecular",
        )
    except ValueError as e:
        log.warning(f"  CNA data not available: {e}")
        cna = pd.DataFrame(columns=["patientId"])

    # Step 3: Mutations (optional — some cohorts may have limited mutations)
    log.info("Step 3/13 — Mutation data (molecular axis)")
    try:
        mut_raw, mol_profiled_patients = fetch_mutation_data(
            study_id,
            entrez_ids=_entrez_for(ALL_MOLECULAR_GENES),
            gene_names=ALL_MOLECULAR_GENES,
            force_refresh=force_refresh,
            cache_suffix="mutations_molecular",
        )
        mut_wide = _pivot_mutations_to_wide(mut_raw, ALL_MOLECULAR_GENES)
        mol_mut_available = True
    except (ValueError, KeyError) as e:
        log.warning(f"  Mutation data not available or incomplete: {e}")
        mut_wide = pd.DataFrame(columns=["patientId"])
        mol_mut_available = False
        mol_profiled_patients = set()

    # Step 4: CNA data for individual genes (CDKN1A, BAX, ATM)
    log.info("Step 4/13 — CNA data (individual genes: CDKN1A, BAX, ATM)")
    try:
        cna_indiv = fetch_cna_data(
            study_id,
            entrez_ids=_entrez_for(ALL_INDIVIDUAL_GENES),
            gene_names=ALL_INDIVIDUAL_GENES,
            force_refresh=force_refresh,
            cache_suffix="cna_individual",
        )
    except ValueError as e:
        log.warning(f"  Individual genes CNA data not available: {e}")
        cna_indiv = pd.DataFrame(columns=["patientId"])

    # Step 5: Mutations for individual genes (CDKN1A, BAX, ATM)
    log.info("Step 5/13 — Mutation data (individual genes: CDKN1A, BAX, ATM)")
    try:
        mut_indiv_raw, indiv_profiled_patients = fetch_mutation_data(
            study_id,
            entrez_ids=_entrez_for(ALL_INDIVIDUAL_GENES),
            gene_names=ALL_INDIVIDUAL_GENES,
            force_refresh=force_refresh,
            cache_suffix="mutations_individual",
        )
        mut_indiv_wide = _pivot_mutations_to_wide(mut_indiv_raw, ALL_INDIVIDUAL_GENES)
        indiv_mut_available = True
    except (ValueError, KeyError) as e:
        log.warning(f"  Individual genes mutation data not available: {e}")
        mut_indiv_wide = pd.DataFrame(columns=["patientId"])
        indiv_mut_available = False
        indiv_profiled_patients = set()

    # Step 6: mRNA z-scores for molecular axis genes
    log.info("Step 6/13 — mRNA z-scores (molecular axis genes)")
    try:
        mrna_mol = fetch_mrna_data(
            study_id,
            entrez_ids=_entrez_for(ALL_MOLECULAR_GENES),
            gene_names=ALL_MOLECULAR_GENES,
            force_refresh=force_refresh,
            cache_suffix="mrna_molecular",
        )
    except ValueError as e:
        log.warning(f"  mRNA data not available: {e}")
        mrna_mol = pd.DataFrame(columns=["patientId"])

    # Step 7: mRNA z-scores for AR activity genes
    log.info("Step 7/13 — mRNA z-scores (AR activity: KLK3, TMPRSS2, NKX3.1)")
    try:
        mrna_ar_activity = fetch_mrna_data(
            study_id,
            entrez_ids=_entrez_for(ALL_AR_ACTIVITY_GENES),
            gene_names=ALL_AR_ACTIVITY_GENES,
            force_refresh=force_refresh,
            cache_suffix="mrna_ar_activity",
        )
    except ValueError as e:
        log.warning(f"  AR activity mRNA not available: {e}")
        mrna_ar_activity = pd.DataFrame(columns=["patientId"])

    # Step 8: mRNA expression (androgen uptake + adhesion-motility axis)
    log.info("Step 8/13 — mRNA expression (androgen uptake + adhesion-motility: 14 genes)")
    try:
        mrna = fetch_mrna_data(
            study_id,
            entrez_ids=_entrez_for(ALL_EXPRESSION_GENES),
            gene_names=ALL_EXPRESSION_GENES,
            force_refresh=force_refresh,
        )
    except ValueError as e:
        log.warning(f"  Expression mRNA not available: {e}")
        mrna = pd.DataFrame(columns=["patientId"])

    # Step 9: Structural variants
    log.info("Step 9/13 — Structural variants (molecular axis genes)")
    try:
        sv = fetch_sv_data(
            study_id,
            entrez_ids=_entrez_for(ALL_MOLECULAR_GENES),
            gene_names=ALL_MOLECULAR_GENES,
            force_refresh=force_refresh,
        )
    except (ValueError, KeyError) as e:
        log.warning(f"  Structural variants not available: {e}")
        sv = pd.DataFrame(columns=["patientId"])

    # Step 10: CNA data for androgen CNA genes (SLCO2B1, SLCO1B3, AKR1C3)
    log.info("Step 10/13 — CNA data (androgen CNA genes: SLCO2B1, SLCO1B3, AKR1C3)")
    try:
        cna_androgen = fetch_cna_data(
            study_id,
            entrez_ids=_entrez_for(ALL_ANDROGEN_CNA_GENES),
            gene_names=ALL_ANDROGEN_CNA_GENES,
            force_refresh=force_refresh,
            cache_suffix="cna_androgen",
        )
    except ValueError as e:
        log.warning(f"  Androgen CNA data not available: {e}")
        cna_androgen = pd.DataFrame(columns=["patientId"])

    # Step 11: Mutations for androgen CNA genes
    log.info("Step 11/13 — Mutation data (androgen CNA genes)")
    try:
        mut_androgen_raw, androgen_profiled_patients = fetch_mutation_data(
            study_id,
            entrez_ids=_entrez_for(ALL_ANDROGEN_CNA_GENES),
            gene_names=ALL_ANDROGEN_CNA_GENES,
            force_refresh=force_refresh,
            cache_suffix="mutations_androgen",
        )
        mut_androgen_wide = _pivot_mutations_to_wide(mut_androgen_raw, ALL_ANDROGEN_CNA_GENES)
        androgen_mut_available = True
    except (ValueError, KeyError) as e:
        log.warning(f"  Androgen mutation data not available: {e}")
        mut_androgen_wide = pd.DataFrame(columns=["patientId"])
        androgen_mut_available = False
        androgen_profiled_patients = set()

    # Step 12: CNA data for adhesion/motility CNA genes
    log.info("Step 12/13 — CNA data (adhesion/motility CNA genes)")
    try:
        cna_adhesion = fetch_cna_data(
            study_id,
            entrez_ids=_entrez_for(ALL_ADHESION_CNA_GENES),
            gene_names=ALL_ADHESION_CNA_GENES,
            force_refresh=force_refresh,
            cache_suffix="cna_adhesion",
        )
    except ValueError as e:
        log.warning(f"  Adhesion CNA data not available: {e}")
        cna_adhesion = pd.DataFrame(columns=["patientId"])

    # Step 13: Mutations for adhesion/motility CNA genes
    log.info("Step 13/13 — Mutation data (adhesion/motility CNA genes)")
    try:
        mut_adhesion_raw, adhesion_profiled_patients = fetch_mutation_data(
            study_id,
            entrez_ids=_entrez_for(ALL_ADHESION_CNA_GENES),
            gene_names=ALL_ADHESION_CNA_GENES,
            force_refresh=force_refresh,
            cache_suffix="mutations_adhesion",
        )
        mut_adhesion_wide = _pivot_mutations_to_wide(mut_adhesion_raw, ALL_ADHESION_CNA_GENES)
        adhesion_mut_available = True
    except (ValueError, KeyError) as e:
        log.warning(f"  Adhesion mutation data not available: {e}")
        mut_adhesion_wide = pd.DataFrame(columns=["patientId"])
        adhesion_mut_available = False
        adhesion_profiled_patients = set()

    # Step 14: mRNA z-scores for crowding / mechanobiology axes (Axes 1–6)
    # Expression only; PTEN stratifier comes from CNA (already fetched in Step 2).
    # Analyzed for TCGA-PRAD only (module8), but fetched wherever mRNA exists.
    log.info(f"Step 14/14 — mRNA z-scores (crowding axes: {len(ALL_CROWDING_GENES)} genes)")
    try:
        mrna_crowding = fetch_mrna_data(
            study_id,
            entrez_ids=_entrez_for(ALL_CROWDING_GENES),
            gene_names=ALL_CROWDING_GENES,
            force_refresh=force_refresh,
            cache_suffix="mrna_crowding",
        )
    except ValueError as e:
        log.warning(f"  Crowding mRNA not available: {e}")
        mrna_crowding = pd.DataFrame(columns=["patientId"])

    # Assemble master
    log.info("Assembling master dataframe...")
    master = _merge_left_all(
        clinical, cna, mut_wide, cna_indiv, mut_indiv_wide,
        mrna_mol, mrna_ar_activity, mrna, sv,
        cna_androgen, mut_androgen_wide,
        cna_adhesion, mut_adhesion_wide,
        mrna_crowding,
    )

    # Mutation flag columns: 0-fill ONLY for patients confirmed to be in the
    # mutation profile (i.e. sequenced). Absent from mutation records AND in the
    # profile → true wildtype (0). Not in the profile at all → NaN (unknown).
    def _fill_mut_flags(master, genes, profiled_patients):
        if not profiled_patients:
            # Profile unavailable for this cohort — all mutation flags stay NaN
            for gene in genes:
                col = f"{gene}_MUT"
                if col not in master.columns:
                    master[col] = pd.NA
            return master
        is_profiled = master["patientId"].isin(profiled_patients)
        n_profiled  = is_profiled.sum()
        n_not       = (~is_profiled).sum()
        log.info(f"  Mutation fill: {n_profiled} profiled (→0 if absent), "
                 f"{n_not} not profiled (→NaN)")
        for gene in genes:
            col = f"{gene}_MUT"
            if col not in master.columns:
                master[col] = pd.NA
            master.loc[is_profiled, col] = (
                master.loc[is_profiled, col].fillna(0).astype(int)
            )
            # Patients not in the profile retain NaN — never assigned wildtype
        return master

    master = _fill_mut_flags(master, ALL_MOLECULAR_GENES,   mol_profiled_patients)
    master = _fill_mut_flags(master, ALL_INDIVIDUAL_GENES,  indiv_profiled_patients)
    master = _fill_mut_flags(master, ALL_ANDROGEN_CNA_GENES, androgen_profiled_patients)
    master = _fill_mut_flags(master, ALL_ADHESION_CNA_GENES, adhesion_profiled_patients)

    log.info(f"Master assembled: "
             f"{len(master)} patients × {len(master.columns)} columns")
    _log_column_summary(master)

    master.to_csv(master_cache, index=False)
    log.info(f"Saved: {master_cache}")
    return master


# ─────────────────────────────────────────────────────────────────────────────
# ALL DATASETS
# ─────────────────────────────────────────────────────────────────────────────

def download_all(force_refresh=False):
    """
    Downloads all datasets.

    Returns
    -------
    dict of {dataset_key: master_DataFrame}
    """
    results = {}
    for key in DATASETS:
        results[key] = download_dataset(key, force_refresh=force_refresh)

    log.info("\n" + "="*60)
    log.info("ALL DOWNLOADS COMPLETE")
    log.info("="*60)
    for key, df in results.items():
        log.info(f"  {DATASETS[key]['label']:<25} "
                 f"{len(df):>4} patients, "
                 f"{len(df.columns):>3} columns")
    return results


# ─────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Download cBioPortal data for KM analysis"
    )
    parser.add_argument(
        "--refresh", action="store_true",
        help="Force re-fetch from API, bypassing local cache"
    )
    parser.add_argument(
        "--dataset", default="all",
        choices=list(DATASETS.keys()) + ["all"],
        help="Which dataset to download (default: all)"
    )
    args = parser.parse_args()

    if args.dataset == "all":
        download_all(force_refresh=args.refresh)
    else:
        download_dataset(args.dataset, force_refresh=args.refresh)
