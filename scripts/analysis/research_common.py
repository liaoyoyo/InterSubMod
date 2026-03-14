#!/usr/bin/env python3
"""Shared helpers for research-round analysis scripts."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_json(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Dict) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)


def parse_variant_counts(path: Path) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    if not path.exists():
        return counts
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            try:
                counts[key.strip()] = int(value.strip())
            except ValueError:
                continue
    return counts


def as_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def to_float(value, default: float = math.nan) -> float:
    if value is None:
        return default
    try:
        text = str(value).strip()
        if not text:
            return default
        return float(text)
    except (TypeError, ValueError):
        return default


def to_int(value, default: int = 0) -> int:
    if value is None:
        return default
    try:
        text = str(value).strip()
        if not text:
            return default
        return int(float(text))
    except (TypeError, ValueError):
        return default


def safe_div(num: float, den: float) -> float:
    if den == 0:
        return 0.0
    return num / den


def compute_metrics(tp: int, fp: int, truth_total: int) -> Dict[str, float]:
    fn = max(truth_total - tp, 0)
    precision = safe_div(tp, tp + fp)
    recall = safe_div(tp, truth_total)
    f1 = safe_div(2 * precision * recall, precision + recall) if (precision + recall) else 0.0
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "calls_total": tp + fp,
        "truth_total": truth_total,
    }


def infer_platform(sample: str) -> str:
    sample_upper = (sample or "").upper()
    if "DORADO" in sample_upper:
        return "ONT_Dorado"
    if sample_upper == "HCC1395":
        return "ONT_5kHz"
    return "ONT"


def infer_truth_set(sample: str) -> str:
    sample_upper = (sample or "").upper()
    if sample_upper.startswith("HCC1395"):
        return "SEQC2_v1.2.1_HC_SNV"
    if sample_upper == "COLO829":
        return "NYGC"
    return "orthogonal-tools"


def infer_caller_name(vcf_path: str) -> str:
    text = (vcf_path or "").lower()
    if "deepsomatic" in text:
        return "DeepSomatic"
    if "clairs-to" in text or "clairs_to" in text:
        return "ClairS-TO"
    if "clairs" in text:
        return "ClairS"
    return "CallerInput"


def load_tsv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def write_tsv_rows(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict]) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def markdown_table(columns: Sequence[str], rows: Sequence[Dict]) -> str:
    if not rows:
        return "_No rows_"
    lines = []
    lines.append("| " + " | ".join(columns) + " |")
    lines.append("| " + " | ".join(["---"] * len(columns)) + " |")
    for row in rows:
        values = [str(row.get(col, "")) for col in columns]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)
