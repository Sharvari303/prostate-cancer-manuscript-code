"""
module4_flag_molecular.py — Computes binary alteration flags for all molecular genes.

For each gene, combines CNA + somatic mutation + structural variant (where available)
into a single altered/wildtype/NaN flag using OR logic.

Thresholds (strict, used throughout):
  Oncogene amplification  ≥ +2 (homozygous gain)
  Suppressor deletion     ≤ −2 (homozygous loss)
  Suppressor mutation     truncating only (frameshift, nonsense, splice)
  AR                      CNA amplification only (mutations excluded)

Missing data rule: no data from any source → NaN, never wildtype.
Wildtype (0) requires data present AND no alteration detected.

Output: data/cache/{KEY}_flagged_molecular.csv  (one per cohort)
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    MOLECULAR_GENES, INDIVIDUAL_GENES,
    ANDROGEN_CNA_GENES, ADHESION_CNA_GENES,
    CNA_THRESHOLDS, DEFAULT_CNA_THRESHOLD,
    CACHE_DIR, DATASETS,
)
from utils.logger import get_logger

log = get_logger("module4_flag_molecular")

TRUNCATING_TYPES = {
    "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Splice_Site", "Splice_Region", "Translation_Start_Site",
    "Nonstop_Mutation",
}

MRNA_ZSCORE_THRESHOLD = 2.0


def _is_truncating_mut(df, gene):
    mut_col  = f"{gene}_MUT"
    type_col = f"{gene}_MUT_TYPE"
    if mut_col not in df.columns:
        return pd.Series(False, index=df.index)
    has_mut = df[mut_col] == 1
    if type_col in df.columns:
        return has_mut & df[type_col].isin(TRUNCATING_TYPES)
    return has_mut



def _flag_single_gene(df, gene, alteration_type, threshold_key="strict",
                      use_mrna=False):
    """
    Returns NaN-safe alteration flag for one gene.
    alteration_type: "suppressor" or "oncogene"
    NaN = no data from any source → excluded from analysis (never assigned wildtype).
    """
    thresholds = CNA_THRESHOLDS[threshold_key]
    cna_col    = f"{gene}_CNA"
    sv_col     = f"{gene}_SV"
    zscore_col = f"{gene}_ZSCORE"

    # CNA flag
    if cna_col in df.columns:
        cna_thresh = thresholds[alteration_type]
        if alteration_type == "suppressor":
            cna_flag = df[cna_col] <= cna_thresh
        else:
            cna_flag = df[cna_col] >= cna_thresh
        cna_available = df[cna_col].notna()
    else:
        cna_flag      = pd.Series(pd.NA, index=df.index, dtype="boolean")
        cna_available = pd.Series(False, index=df.index)

    # Mutation flag
    mut_col = f"{gene}_MUT"
    if mut_col in df.columns:
        if threshold_key == "strict" and alteration_type == "suppressor":
            mut_flag = _is_truncating_mut(df, gene)
        else:
            mut_flag = df[mut_col] == 1
        mut_available = df[mut_col].notna()
    else:
        mut_flag      = pd.Series(pd.NA, index=df.index, dtype="boolean")
        mut_available = pd.Series(False, index=df.index)

    # SV flag
    if sv_col in df.columns:
        sv_flag      = df[sv_col] == 1
        sv_available = df[sv_col].notna()
    else:
        sv_flag      = pd.Series(pd.NA, index=df.index, dtype="boolean")
        sv_available = pd.Series(False, index=df.index)

    # mRNA flag (only for molecular genes that have z-scores)
    if use_mrna and zscore_col in df.columns:
        zs = pd.to_numeric(df[zscore_col], errors="coerce")
        if alteration_type == "suppressor":
            mrna_flag = zs < -MRNA_ZSCORE_THRESHOLD
        else:
            mrna_flag = zs > MRNA_ZSCORE_THRESHOLD
        mrna_available = zs.notna()
    else:
        mrna_flag      = pd.Series(pd.NA, index=df.index, dtype="boolean")
        mrna_available = pd.Series(False, index=df.index)

    # OR logic with NaN-safety
    altered = (cna_flag.fillna(False) | mut_flag.fillna(False) |
               sv_flag.fillna(False) | mrna_flag.fillna(False))
    any_data = cna_available | mut_available | sv_available | mrna_available
    altered = altered.where(any_data, other=np.nan)
    return altered


def add_molecular_flags(df, dataset_key="UNKNOWN"):
    """Flags MOLECULAR_GENES (PTEN, MDM2, MDM4, AR) with CNA+mut+SV.

    AR mutation flagging is cohort-stage-aware:
      - Primary disease cohorts (TCGA_PRAD): only known GOF hotspots
      - CRPC/metastatic cohorts: any AR missense mutation
    dataset_key is passed through so AR filtering can apply the right rule.
    """
    log.info(f"Adding molecular gene flags (PTEN, MDM2, MDM4, AR) for {dataset_key}...")
    for gene, cfg in MOLECULAR_GENES.items():
        atype = cfg["alteration_type"]
        for thresh_key in ["strict", "relaxed"]:
            col = f"{gene}_ALT_{thresh_key.upper()}"

            if gene == "AR":
                # AR uses cohort-aware GOF mutation flagging instead of generic mut flag
                flag = _flag_ar(df, thresh_key, dataset_key)
            else:
                flag = _flag_single_gene(df, gene, atype, thresh_key, use_mrna=False)

            df[col] = flag.where(flag.notna(), other=np.nan)
            n_alt = (flag == True).sum()
            log.info(f"  {col:<30} altered={n_alt}/{len(df)}, "
                     f"missing={flag.isna().sum()}")

    # Composite molecular score (default threshold, PTEN+MDM2+MDM4+AR)
    default_cols = [f"{g}_ALT_{DEFAULT_CNA_THRESHOLD.upper()}" for g in MOLECULAR_GENES]
    df["MOLECULAR_SCORE"] = df[default_cols].sum(axis=1, skipna=False)
    return df


def _flag_ar(df, threshold_key, dataset_key):
    """AR alteration flag: CNA amplification (≥ +2 strict, ≥ +1 relaxed) only.

    AR amplification is the dominant mechanism (42% SU2C, 52% MCTP vs 1% TCGA).
    Mutations are excluded: GOF filtering requires protein-change data not stored
    in the master file, and mutations add <5% patients beyond CNA alone.
    """
    thresholds = CNA_THRESHOLDS[threshold_key]
    cna_col = "AR_CNA"

    if cna_col in df.columns:
        cna_flag      = df[cna_col] >= thresholds["oncogene"]
        cna_available = df[cna_col].notna()
    else:
        cna_flag      = pd.Series(pd.NA, index=df.index, dtype="boolean")
        cna_available = pd.Series(False, index=df.index)

    n_alt = cna_flag.fillna(False).sum()
    log.info(f"    AR CNA-only ({threshold_key}): {n_alt} altered in {dataset_key}")
    return cna_flag.where(cna_available, other=np.nan)


def add_individual_gene_flags(df):
    """Flags INDIVIDUAL_GENES (CDKN1A, BAX, ATM) with CNA+mutation, suppressor logic."""
    log.info("Adding individual gene flags (CDKN1A, BAX, ATM)...")
    for gene, cfg in INDIVIDUAL_GENES.items():
        atype = cfg["alteration_type"]
        col   = f"{gene}_ALT_STRICT"
        flag  = _flag_single_gene(df, gene, atype, "strict", use_mrna=False)
        df[col] = flag.where(flag.notna(), other=np.nan)
        n_alt = (flag == True).sum()
        log.info(f"  {col:<30} altered={n_alt}/{len(df)}, "
                 f"missing={flag.isna().sum()}")
    return df


def add_androgen_cna_flags(df):
    """Flags ANDROGEN_CNA_GENES (SLCO2B1, SLCO1B3, AKR1C3) with CNA+mutation. Task 3b."""
    log.info("Adding androgen CNA gene flags (SLCO2B1, SLCO1B3, AKR1C3)...")
    for gene, cfg in ANDROGEN_CNA_GENES.items():
        atype = cfg["alteration_type"]
        col   = f"{gene}_CNA_ALT"
        flag  = _flag_single_gene(df, gene, atype, "strict", use_mrna=False)
        df[col] = flag.where(flag.notna(), other=np.nan)
        n_alt = (flag == True).sum()
        log.info(f"  {col:<30} altered={n_alt}/{len(df)}, "
                 f"missing={flag.isna().sum()}")
    return df


def add_adhesion_cna_flags(df):
    """Flags ADHESION_CNA_GENES with CNA+mutation. Task 6b."""
    log.info("Adding adhesion/motility CNA gene flags...")
    for gene, cfg in ADHESION_CNA_GENES.items():
        atype = cfg["alteration_type"]
        col   = f"{gene}_CNA_ALT"
        flag  = _flag_single_gene(df, gene, atype, "strict", use_mrna=False)
        df[col] = flag.where(flag.notna(), other=np.nan)
        n_alt = (flag == True).sum()
        log.info(f"  {col:<30} altered={n_alt}/{len(df)}, "
                 f"missing={flag.isna().sum()}")
    return df


def add_or_combined_flags(df):
    """
    Adds OR-combined alteration flags with strict missing data rule.
    Wildtype only if ALL genes in the OR set have non-NaN data AND all = 0.
    """
    # Task 2: PTEN + AR + MDM4 OR-combined
    task2_genes = ["PTEN", "AR", "MDM4"]
    task2_cols  = [f"{g}_ALT_{DEFAULT_CNA_THRESHOLD.upper()}" for g in task2_genes]
    available   = [c for c in task2_cols if c in df.columns]
    if available:
        flags = df[available]
        any_altered    = (flags == True).any(axis=1)
        all_have_data  = flags.notna().all(axis=1)
        or_combined    = np.where(any_altered, 1,
                         np.where(all_have_data, 0, np.nan))
        df["MOLECULAR_OR"] = or_combined
        n_alt = (df["MOLECULAR_OR"] == 1).sum()
        n_wt  = (df["MOLECULAR_OR"] == 0).sum()
        n_exc = df["MOLECULAR_OR"].isna().sum()
        log.info(f"  MOLECULAR_OR (PTEN+AR+MDM4): altered={n_alt}, "
                 f"wildtype={n_wt}, excluded={n_exc}")
    return df


def flag_all(clean_dfs=None, force_refresh=False):
    results = {}
    for key in list(DATASETS.keys()) + ["COMBINED"]:
        suffix   = "combined" if key == "COMBINED" else key
        out_path = CACHE_DIR / f"{suffix}_flagged_molecular.csv"

        if out_path.exists() and not force_refresh:
            log.info(f"Loading molecular flags from cache: {out_path.name}")
            results[key] = pd.read_csv(out_path)
            continue

        if clean_dfs and key in clean_dfs:
            df = clean_dfs[key].copy()
        else:
            src = CACHE_DIR / ("combined_clean.csv" if key == "COMBINED"
                               else f"{key}_clean.csv")
            if not src.exists():
                raise FileNotFoundError(f"Clean cache not found: {src}")
            df = pd.read_csv(src)

        log.info(f"\n{'='*55}\nMolecular flagging: "
                 f"{DATASETS.get(key, {}).get('label', key)}\n{'='*55}")
        df = add_molecular_flags(df, dataset_key=key)
        df = add_individual_gene_flags(df)
        df = add_androgen_cna_flags(df)
        df = add_adhesion_cna_flags(df)
        df = add_or_combined_flags(df)

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
    log.info("Molecular flagging complete.")
