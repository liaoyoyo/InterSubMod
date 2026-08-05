#!/usr/bin/env python3
"""Deprecated CLI redirect to the topic-owned matched-normal analyzer.

The maintained implementation is under the 20260715 cooccurrence-validation
research topic. This wrapper exists only to preserve the accidental top-level
path without creating a second implementation.
"""

from __future__ import annotations

import runpy
from pathlib import Path


TARGET = (
    Path(__file__).resolve().parents[1]
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
    / "scripts"
    / "analyze_matched_normal_candidate_controls.py"
)


if __name__ == "__main__":
    runpy.run_path(str(TARGET), run_name="__main__")
