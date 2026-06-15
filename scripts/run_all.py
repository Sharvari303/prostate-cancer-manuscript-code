"""
run_all.py — Master pipeline runner.

Runs all 6 modules in order. Each module reads cached outputs from the
previous step, so you can resume mid-pipeline with --start or run a
single module with --only.

Usage:
  python run_all.py                    # full pipeline
  python run_all.py --refresh          # re-download everything from cBioPortal
  python run_all.py --start module6    # resume from a specific module
  python run_all.py --only module7     # run one module only
"""
import sys
import time
import argparse
import traceback
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import CACHE_DIR
from utils.logger import get_logger

log = get_logger("run_all")

MODULE_ORDER = [
    "module2_download",
    "module3_clean",
    "module4_flag_molecular",
    "module5_flag_expression",
    "module6_comprehensive_km",
    "module7_cox_regression",
    "module8_crowding",
]


def import_module(name):
    """Dynamically imports a module by name."""
    import importlib
    return importlib.import_module(name)


def run_pipeline(start_from=None, only=None, force_refresh=False):
    log.info("="*60)
    log.info(f"KM SURVIVAL ANALYSIS PIPELINE")
    log.info(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log.info("="*60)

    modules_to_run = MODULE_ORDER
    if start_from:
        idx = next((i for i, m in enumerate(MODULE_ORDER)
                    if start_from in m), None)
        if idx is None:
            log.error(f"Module '{start_from}' not found.")
            return
        modules_to_run = MODULE_ORDER[idx:]
    if only:
        modules_to_run = [m for m in MODULE_ORDER if only in m]
        if not modules_to_run:
            log.error(f"Module '{only}' not found.")
            return

    # State passed between modules to avoid redundant disk reads
    state = {
        "master_dfs":    None,
        "clean_dfs":     None,
        "molecular_dfs": None,
        "flagged_dfs":   None,
        "summaries":     {},
    }

    results = {}
    for mod_name in modules_to_run:
        log.info(f"\n{'─'*60}")
        log.info(f"Running: {mod_name}")
        log.info(f"{'─'*60}")
        t0 = time.time()

        try:
            mod = import_module(mod_name)

            # ── Module 2: Download ──────────────────────────────────
            if mod_name == "module2_download":
                state["master_dfs"] = mod.download_all(
                    force_refresh=force_refresh)

            # ── Module 3: Clean ─────────────────────────────────────
            elif "module3" in mod_name:
                state["clean_dfs"] = mod.run_cleaning(
                    master_dfs=state["master_dfs"],
                    force_refresh=force_refresh)

            # ── Module 4: Molecular flags ───────────────────────────
            elif "module4" in mod_name:
                state["molecular_dfs"] = mod.flag_all(
                    clean_dfs=state["clean_dfs"],
                    force_refresh=force_refresh)

            # ── Module 5: Expression flags ──────────────────────────
            elif "module5" in mod_name:
                state["flagged_dfs"] = mod.flag_all(
                    molecular_dfs=state["molecular_dfs"],
                    force_refresh=force_refresh)

            # ── Module 6: Comprehensive KM plots (Tier 1/2/3) ──────────────────
            elif "module6_comprehensive_km" in mod_name:
                mod.main()
                log.info("  57 KM figures + statistics table generated")

            # ── Module 7: Cox regression (univariate + multivariate) ──────────
            elif "module7_cox_regression" in mod_name:
                mod.main()
                log.info("  Cox regression results table generated")

            # ── Module 8: Crowding / mechanobiology axes (KM + Cox interaction) ─
            elif "module8_crowding" in mod_name:
                mod.main()
                log.info("  Crowding axes KM + Cox interaction tables generated")

            elapsed = time.time() - t0
            log.info(f"✓ {mod_name} completed in {elapsed:.1f}s")
            results[mod_name] = "SUCCESS"

        except Exception as e:
            elapsed = time.time() - t0
            log.error(f"✗ {mod_name} FAILED after {elapsed:.1f}s")
            log.error(f"  Error: {e}")
            log.debug(traceback.format_exc())
            results[mod_name] = f"FAILED: {e}"
            # Continue to next module rather than aborting
            continue

    # Final summary
    log.info(f"\n{'='*60}")
    log.info("PIPELINE COMPLETE — Summary")
    log.info("="*60)
    for mod, status in results.items():
        icon = "✓" if status == "SUCCESS" else "✗"
        log.info(f"  {icon} {mod:<45} {status}")
    log.info(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="KM Survival Analysis Pipeline"
    )
    parser.add_argument(
        "--refresh", action="store_true",
        help="Force re-fetch all data from cBioPortal API")
    parser.add_argument(
        "--start", default=None, metavar="MODULE",
        help="Resume pipeline from this module "
             "(e.g. module6_km_molecular_individual)")
    parser.add_argument(
        "--only", default=None, metavar="MODULE",
        help="Run only one specific module")
    args = parser.parse_args()

    run_pipeline(
        start_from=args.start,
        only=args.only,
        force_refresh=args.refresh,
    )
