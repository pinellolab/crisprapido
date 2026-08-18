#!/usr/bin/env python3
import runpy
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
runpy.run_path(
    str(SCRIPT_DIR.parent / "chr2_100_guides" / "collect_peak_rss.py"),
    run_name="__main__",
)

