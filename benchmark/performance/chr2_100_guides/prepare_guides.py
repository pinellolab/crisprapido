#!/usr/bin/env python3
import os
import runpy
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.environ["BENCH_PACKAGE_DIR"] = str(SCRIPT_DIR)
runpy.run_path(str(SCRIPT_DIR.parent / "chr22_500_guides" / "prepare_guides.py"), run_name="__main__")
