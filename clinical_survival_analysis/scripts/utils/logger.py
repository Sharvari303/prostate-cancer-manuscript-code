"""
utils/logger.py
─────────────────────────────────────────────────────────────────────────────
Shared logger for all pipeline modules.
Writes to both console and a timestamped log file.
─────────────────────────────────────────────────────────────────────────────
"""

import logging
import sys
from datetime import datetime
from pathlib import Path


def get_logger(name: str, log_dir: Path = None) -> logging.Logger:
    """
    Returns a logger that writes to console + log file.

    Parameters
    ----------
    name     : module name, e.g. "module2_download"
    log_dir  : directory for log file (uses config.LOG_DIR if None)
    """
    if log_dir is None:
        # Import here to avoid circular import
        from config import LOG_DIR
        log_dir = LOG_DIR

    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # Avoid adding duplicate handlers on reimport
    if logger.handlers:
        return logger

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler — INFO and above
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    # File handler — DEBUG and above
    timestamp = datetime.now().strftime("%Y%m%d")
    log_file  = log_dir / f"{timestamp}_pipeline.log"
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger
