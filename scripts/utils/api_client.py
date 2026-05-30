"""
utils/api_client.py
─────────────────────────────────────────────────────────────────────────────
cBioPortal REST API client.
Handles:
  - HTTP requests with retry logic and timeout
  - Local CSV caching (avoids re-fetching on every run)
  - Profile discovery (handles naming differences across datasets)
  - Polite rate limiting between calls
─────────────────────────────────────────────────────────────────────────────
"""

import time
import requests
import pandas as pd
from pathlib import Path
from typing import Any, Dict, List, Optional

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import (
    CBIOPORTAL_BASE_URL,
    API_TIMEOUT_SECONDS,
    API_RETRY_ATTEMPTS,
    API_RETRY_DELAY,
    API_POLITE_DELAY,
    CACHE_DIR,
    MRNA_PROFILE_OVERRIDES,
)
from utils.logger import get_logger

log = get_logger("api_client")


# ─────────────────────────────────────────────────────────────────────────────
# CORE HTTP
# ─────────────────────────────────────────────────────────────────────────────

def _get(endpoint: str,
         params: Optional[Dict] = None) -> Any:
    """
    GET request with retry logic.
    Returns parsed JSON response.
    """
    url = f"{CBIOPORTAL_BASE_URL}/{endpoint.lstrip('/')}"

    for attempt in range(1, API_RETRY_ATTEMPTS + 1):
        try:
            log.debug(f"GET {url} params={params} (attempt {attempt})")
            r = requests.get(url,
                             params=params,
                             timeout=API_TIMEOUT_SECONDS)
            r.raise_for_status()
            time.sleep(API_POLITE_DELAY)
            return r.json()

        except requests.exceptions.HTTPError as e:
            log.warning(f"HTTP error on attempt {attempt}: {e}")
        except requests.exceptions.ConnectionError as e:
            log.warning(f"Connection error on attempt {attempt}: {e}")
        except requests.exceptions.Timeout:
            log.warning(f"Timeout on attempt {attempt}")

        if attempt < API_RETRY_ATTEMPTS:
            time.sleep(API_RETRY_DELAY * attempt)

    raise RuntimeError(f"API GET failed after {API_RETRY_ATTEMPTS} "
                       f"attempts: {url}")


def _post(endpoint: str, body: Dict, params: Optional[Dict] = None) -> Any:
    """
    POST request with retry logic.
    Returns parsed JSON response.
    """
    url = f"{CBIOPORTAL_BASE_URL}/{endpoint.lstrip('/')}"

    for attempt in range(1, API_RETRY_ATTEMPTS + 1):
        try:
            log.debug(f"POST {url} (attempt {attempt})")
            r = requests.post(url,
                              json=body,
                              params=params,
                              timeout=API_TIMEOUT_SECONDS)
            r.raise_for_status()
            time.sleep(API_POLITE_DELAY)
            return r.json()

        except requests.exceptions.HTTPError as e:
            log.warning(f"HTTP error on attempt {attempt}: {e}")
        except requests.exceptions.ConnectionError as e:
            log.warning(f"Connection error on attempt {attempt}: {e}")
        except requests.exceptions.Timeout:
            log.warning(f"Timeout on attempt {attempt}")

        if attempt < API_RETRY_ATTEMPTS:
            time.sleep(API_RETRY_DELAY * attempt)

    raise RuntimeError(f"API POST failed after {API_RETRY_ATTEMPTS} "
                       f"attempts: {url}")


# ─────────────────────────────────────────────────────────────────────────────
# CACHING
# ─────────────────────────────────────────────────────────────────────────────

def _cache_path(study_id: str, data_type: str) -> Path:
    return CACHE_DIR / f"{study_id}_{data_type}.csv"


def _load_cache(cache_file: Path) -> Optional[pd.DataFrame]:
    if cache_file.exists():
        log.info(f"  Cache hit: {cache_file.name}")
        return pd.read_csv(cache_file)
    return None


def _save_cache(df: pd.DataFrame, cache_file: Path) -> None:
    df.to_csv(cache_file, index=False)
    log.info(f"  Cached to: {cache_file.name}")


# ─────────────────────────────────────────────────────────────────────────────
# PROFILE DISCOVERY
# Handles naming inconsistencies across cBioPortal datasets
# ─────────────────────────────────────────────────────────────────────────────

def get_molecular_profiles(study_id: str) -> Dict[str, str]:
    """
    Returns dict of {profile_id: alteration_type} for a study.
    """
    data = _get(f"studies/{study_id}/molecular-profiles")
    return {
        p["molecularProfileId"]: p["molecularAlterationType"]
        for p in data
    }


def find_profile(study_id: str,
                 target_type: str,
                 preferred_suffix: Optional[str] = None) -> str:
    """
    Finds the molecular profile ID for a given alteration type.

    Parameters
    ----------
    study_id        : cBioPortal study ID
    target_type     : e.g. "COPY_NUMBER_ALTERATION", "MUTATION_EXTENDED",
                           "MRNA_EXPRESSION"
    preferred_suffix: optional suffix hint to resolve ambiguity
                      e.g. "Zscores", "gistic"

    Returns
    -------
    profile_id string
    """
    profiles = get_molecular_profiles(study_id)
    log.debug(f"Available profiles for {study_id}: "
              f"{list(profiles.keys())}")

    # Filter by alteration type
    matches = {pid: ptype
               for pid, ptype in profiles.items()
               if target_type in ptype}

    if not matches:
        raise ValueError(
            f"No profile of type '{target_type}' found "
            f"for study '{study_id}'.\n"
            f"Available: {profiles}"
        )

    if len(matches) == 1:
        pid = list(matches.keys())[0]
        log.debug(f"Profile resolved: {pid}")
        return pid

    # Multiple matches — use suffix hint
    if preferred_suffix:
        for pid in matches:
            if preferred_suffix.lower() in pid.lower():
                log.debug(f"Profile resolved via suffix hint: {pid}")
                return pid

    # Fallback: return first match and warn
    pid = list(matches.keys())[0]
    log.warning(f"Multiple profiles matched for {target_type} "
                f"in {study_id}. Using: {pid}. "
                f"All matches: {list(matches.keys())}")
    return pid


# ─────────────────────────────────────────────────────────────────────────────
# STUDY INFO
# ─────────────────────────────────────────────────────────────────────────────

def get_study_info(study_id: str) -> Dict:
    """Returns metadata dict for a study."""
    return _get(f"studies/{study_id}")


def get_sample_list_id(study_id: str) -> str:
    """Returns the 'all samples' list ID for a study."""
    return f"{study_id}_all"


def get_profiled_patients(study_id: str, profile_id: str) -> set:
    """
    Returns the set of patientIds confirmed to be sequenced in this study.

    Uses the study's _sequenced sample list via /studies/{id}/samples, which
    returns patientId alongside sampleId. This is the correct endpoint —
    /molecular-profiles/{id}/samples does not exist in the cBioPortal API.
    """
    sequenced_list_id = f"{study_id}_sequenced"
    data = _get(f"studies/{study_id}/samples",
                params={"sampleListId": sequenced_list_id,
                        "projection":   "SUMMARY",
                        "pageSize":     10000})
    patient_ids = {item["patientId"] for item in data if "patientId" in item}
    log.debug(f"  Profiled patients ({sequenced_list_id}): {len(patient_ids)}")
    return patient_ids


# ─────────────────────────────────────────────────────────────────────────────
# CLINICAL DATA
# ─────────────────────────────────────────────────────────────────────────────

def fetch_clinical_data(study_id: str,
                        force_refresh: bool = False) -> pd.DataFrame:
    """
    Fetches all patient-level clinical attributes.
    Returns wide-format DataFrame: one row per patient,
    one column per clinical attribute.
    """
    cache_file = _cache_path(study_id, "clinical")
    if not force_refresh:
        cached = _load_cache(cache_file)
        if cached is not None:
            return cached

    log.info(f"Fetching clinical data for {study_id}...")
    data = _get(
        f"studies/{study_id}/clinical-data",
        params={"clinicalDataType": "PATIENT",
                "projection":       "SUMMARY"}
    )

    if not data:
        raise ValueError(f"No clinical data returned for {study_id}")

    df = pd.DataFrame(data)
    df = (df
          .pivot(index="patientId",
                 columns="clinicalAttributeId",
                 values="value")
          .reset_index())

    log.info(f"  Clinical data: {len(df)} patients, "
             f"{len(df.columns)} attributes")
    _save_cache(df, cache_file)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# CNA DATA
# ─────────────────────────────────────────────────────────────────────────────

def fetch_cna_data(study_id: str,
                   entrez_ids: List[int],
                   gene_names: List[str],
                   force_refresh: bool = False,
                   cache_suffix: str = "cna") -> pd.DataFrame:
    """
    Fetches discrete CNA (GISTIC) scores.
    Returns wide DataFrame: one row per patient,
    one column per gene named {GENE}_CNA.
    """
    cache_file = _cache_path(study_id, cache_suffix)
    if not force_refresh:
        cached = _load_cache(cache_file)
        if cached is not None:
            return cached

    log.info(f"Fetching CNA data for {study_id}...")
    profile_id = find_profile(study_id,
                               "COPY_NUMBER_ALTERATION",
                               preferred_suffix="gistic")

    body = {
        "sampleListId":  get_sample_list_id(study_id),
        "entrezGeneIds": entrez_ids,
    }

    data = _post(
        f"molecular-profiles/{profile_id}/discrete-copy-number/fetch",
        body,
        params={"discreteCopyNumberEventType": "ALL"},
    )

    if not data:
        log.warning(f"No CNA data returned for {study_id}")
        return pd.DataFrame(columns=["patientId"])

    entrez_to_hugo = dict(zip(entrez_ids, gene_names))

    records = []
    for item in data:
        patient_id = item.get("patientId") or item.get("uniquePatientKey")
        hugo       = (item["gene"]["hugoGeneSymbol"]
                      if "gene" in item
                      else entrez_to_hugo.get(item.get("entrezGeneId"), "UNKNOWN"))
        alteration = item.get("alteration", item.get("value"))
        records.append({"patientId": patient_id,
                        "gene":      hugo,
                        "CNA":       alteration})

    df = pd.DataFrame(records)
    df = (df
          .pivot_table(index="patientId",
                       columns="gene",
                       values="CNA",
                       aggfunc="first")
          .reset_index())

    # Rename columns: GENE → GENE_CNA
    df.columns = (["patientId"] +
                  [f"{c}_CNA" for c in df.columns if c != "patientId"])

    log.info(f"  CNA data: {len(df)} patients")
    _save_cache(df, cache_file)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# MUTATION DATA
# ─────────────────────────────────────────────────────────────────────────────

def fetch_mutation_data(study_id: str,
                        entrez_ids: List[int],
                        gene_names: List[str],
                        force_refresh: bool = False,
                        cache_suffix: str = "mutations"):
    """
    Fetches somatic mutation data.

    Returns
    -------
    df : long-format DataFrame, one row per mutation event.
         Columns: patientId, gene, mutationType, variantClassification,
                  proteinChange, mutationStatus
    profiled_patients : set of patientIds confirmed sequenced in this profile.
         Use this to distinguish true wildtype (sequenced, no mutation) from
         missing data (not sequenced). Only patients in this set should have
         their mutation flag filled to 0; all others remain NaN.
    """
    cache_file        = _cache_path(study_id, cache_suffix)
    profiled_cache    = _cache_path(study_id, cache_suffix + "_profiled")

    if not force_refresh:
        cached = _load_cache(cache_file)
        profiled_cached = _load_cache(profiled_cache)
        if cached is not None and profiled_cached is not None:
            profiled_patients = set(profiled_cached["patientId"].tolist())
            return cached, profiled_patients

    log.info(f"Fetching mutation data for {study_id}...")
    profile_id = find_profile(study_id,
                               "MUTATION_EXTENDED",
                               preferred_suffix="mutations")

    # Fetch the set of patients confirmed profiled in this mutation panel
    profiled_patients = get_profiled_patients(study_id, profile_id)
    log.info(f"  Patients in mutation profile: {len(profiled_patients)}")

    records = []
    for gene, eid in zip(gene_names, entrez_ids):
        data = _get(
            f"molecular-profiles/{profile_id}/mutations",
            params={
                "sampleListId":   get_sample_list_id(study_id),
                "entrezGeneId":   eid,
                "projection":     "DETAILED",
            }
        )
        for item in data:
            patient_id = (item.get("patientId") or
                          item.get("uniquePatientKey", ""))
            records.append({
                "patientId":             patient_id,
                "gene":                  gene,
                "mutationType":          item.get("mutationType", ""),
                "variantClassification": item.get("variantClassification", ""),
                "proteinChange":         item.get("proteinChange", ""),
                "mutationStatus":        item.get("mutationStatus", ""),
            })

    df = pd.DataFrame(records) if records else pd.DataFrame(
        columns=["patientId", "gene", "mutationType",
                 "variantClassification", "proteinChange",
                 "mutationStatus"])

    log.info(f"  Mutation data: {len(df)} mutation events across "
             f"{df['patientId'].nunique() if len(df) else 0} patients")

    _save_cache(df, cache_file)
    _save_cache(pd.DataFrame({"patientId": list(profiled_patients)}),
                profiled_cache)
    return df, profiled_patients


# ─────────────────────────────────────────────────────────────────────────────
# STRUCTURAL VARIANT DATA
# ─────────────────────────────────────────────────────────────────────────────

def fetch_sv_data(study_id: str,
                  entrez_ids: List[int],
                  gene_names: List[str],
                  force_refresh: bool = False) -> pd.DataFrame:
    """
    Fetches structural variant data for a list of genes.
    Returns wide DataFrame: one row per patient,
    one column per gene named {GENE}_SV (1 if any SV, 0 otherwise).
    """
    cache_file = _cache_path(study_id, "sv")
    if not force_refresh:
        cached = _load_cache(cache_file)
        if cached is not None:
            return cached

    log.info(f"Fetching structural variant data for {study_id}...")
    try:
        profile_id = find_profile(study_id, "STRUCTURAL_VARIANT")
    except ValueError as e:
        log.warning(f"  No SV profile for {study_id}: {e}")
        return pd.DataFrame(columns=["patientId"])

    data = _post(
        "structural-variant/fetch",
        {"molecularProfileIds": [profile_id],
         "entrezGeneIds":       entrez_ids}
    )

    if not data:
        log.warning(f"No SV data returned for {study_id}")
        df = pd.DataFrame(columns=["patientId"])
        _save_cache(df, cache_file)
        return df

    records = []
    for item in data:
        patient_id = item.get("patientId")
        # A gene is involved if it appears at either breakpoint
        site1 = item.get("site1EntrezGeneId")
        site2 = item.get("site2EntrezGeneId")
        entrez_to_hugo = dict(zip(entrez_ids, gene_names))
        for eid in [site1, site2]:
            if eid in entrez_to_hugo:
                records.append({"patientId": patient_id,
                                 "gene":      entrez_to_hugo[eid]})

    if not records:
        log.warning(f"No SV records matched target genes for {study_id}")
        df = pd.DataFrame(columns=["patientId"])
        _save_cache(df, cache_file)
        return df

    long_df = pd.DataFrame(records).drop_duplicates()
    long_df["SV"] = 1
    wide = (long_df
            .pivot_table(index="patientId", columns="gene",
                         values="SV", aggfunc="max")
            .reset_index())
    wide.columns = (["patientId"] +
                    [f"{c}_SV" for c in wide.columns if c != "patientId"])

    log.info(f"  SV data: {len(wide)} patients with SVs")
    _save_cache(wide, cache_file)
    return wide


# ─────────────────────────────────────────────────────────────────────────────
# mRNA EXPRESSION DATA
# ─────────────────────────────────────────────────────────────────────────────

def fetch_mrna_data(study_id: str,
                    entrez_ids: List[int],
                    gene_names: List[str],
                    force_refresh: bool = False,
                    cache_suffix: str = "mrna") -> pd.DataFrame:
    """
    Fetches mRNA expression z-scores.
    Returns wide DataFrame: one row per patient,
    one column per gene named {GENE}_ZSCORE.
    """
    cache_file = _cache_path(study_id, cache_suffix)
    if not force_refresh:
        cached = _load_cache(cache_file)
        if cached is not None:
            return cached

    log.info(f"Fetching mRNA expression data for {study_id}...")
    if study_id in MRNA_PROFILE_OVERRIDES:
        profile_id = MRNA_PROFILE_OVERRIDES[study_id]
        log.info(f"  Using explicit profile override: {profile_id}")
    else:
        profile_id = find_profile(study_id,
                                   "MRNA_EXPRESSION",
                                   preferred_suffix="Zscores")

    body = {
        "sampleListId":  get_sample_list_id(study_id),
        "entrezGeneIds": entrez_ids,
    }

    data = _post(
        f"molecular-profiles/{profile_id}/molecular-data/fetch",
        body
    )

    if not data:
        log.warning(f"No mRNA data returned for {study_id}")
        return pd.DataFrame(columns=["patientId"])

    entrez_to_hugo = dict(zip(entrez_ids, gene_names))

    records = []
    for item in data:
        patient_id = item.get("patientId") or item.get("uniquePatientKey")
        hugo       = (item["gene"]["hugoGeneSymbol"]
                      if "gene" in item
                      else entrez_to_hugo.get(item.get("entrezGeneId"), "UNKNOWN"))
        value      = item.get("value")
        records.append({"patientId": patient_id,
                        "gene":      hugo,
                        "zscore":    value})

    df = pd.DataFrame(records)
    df = (df
          .pivot_table(index="patientId",
                       columns="gene",
                       values="zscore",
                       aggfunc="first")
          .reset_index())

    # Rename columns: GENE → GENE_ZSCORE
    df.columns = (["patientId"] +
                  [f"{c}_ZSCORE"
                   for c in df.columns if c != "patientId"])

    log.info(f"  mRNA data: {len(df)} patients")
    _save_cache(df, cache_file)
    return df


def fetch_raw_mrna_data(study_id: str,
                        entrez_ids: List[int],
                        gene_names: List[str],
                        force_refresh: bool = False,
                        cache_suffix: str = "mrna_raw") -> pd.DataFrame:
    """
    Fetches RAW mRNA expression data (not z-scores).
    Returns wide DataFrame: one row per patient,
    one column per gene named {GENE}_RAW.

    Used for harmonization across cohorts.
    """
    cache_file = _cache_path(study_id, cache_suffix)
    if not force_refresh:
        cached = _load_cache(cache_file)
        if cached is not None:
            return cached

    log.info(f"Fetching RAW mRNA expression data for {study_id}...")

    # Find raw mRNA profile (not z-scores)
    try:
        profile_id = find_profile(study_id,
                                   "MRNA_EXPRESSION",
                                   preferred_suffix=None)  # Get any, but not Zscore
        # If we got a Zscore profile, try to find non-Zscore version
        profiles = get_molecular_profiles(study_id)
        mrna_profiles = {pid: ptype for pid, ptype in profiles.items()
                        if "MRNA_EXPRESSION" in ptype}

        # Prefer profiles WITHOUT Zscore
        non_zscore = {pid for pid in mrna_profiles if "Zscore" not in pid.lower()}
        if non_zscore:
            profile_id = list(non_zscore)[0]
            log.debug(f"Using raw mRNA profile: {profile_id}")
        else:
            log.warning(f"No raw mRNA profile found for {study_id}, using z-score profile")

    except ValueError:
        log.warning(f"Could not find raw mRNA profile for {study_id}")
        return pd.DataFrame(columns=["patientId"])

    body = {
        "sampleListId":  get_sample_list_id(study_id),
        "entrezGeneIds": entrez_ids,
    }

    data = _post(
        f"molecular-profiles/{profile_id}/molecular-data/fetch",
        body
    )

    if not data:
        log.warning(f"No raw mRNA data returned for {study_id}")
        return pd.DataFrame(columns=["patientId"])

    entrez_to_hugo = dict(zip(entrez_ids, gene_names))

    records = []
    for item in data:
        patient_id = item.get("patientId") or item.get("uniquePatientKey")
        hugo       = (item["gene"]["hugoGeneSymbol"]
                      if "gene" in item
                      else entrez_to_hugo.get(item.get("entrezGeneId"), "UNKNOWN"))
        value      = item.get("value")
        records.append({"patientId": patient_id,
                        "gene":      hugo,
                        "raw":       value})

    df = pd.DataFrame(records)
    df = (df
          .pivot_table(index="patientId",
                       columns="gene",
                       values="raw",
                       aggfunc="first")
          .reset_index())

    # Rename columns: GENE → GENE_RAW
    df.columns = (["patientId"] +
                  [f"{c}_RAW"
                   for c in df.columns if c != "patientId"])

    log.info(f"  Raw mRNA data: {len(df)} patients")
    _save_cache(df, cache_file)
    return df
