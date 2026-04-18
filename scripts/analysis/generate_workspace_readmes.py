#!/usr/bin/env python3
"""Generate README.md for observation workspaces that lack one.

Extracts metadata from:
  1. round_context.json (if exists): title, script, observation_id, description
  2. Directory name: date, topic keywords
  3. File listing: figure count, data file count, key files

Usage:
    python3 scripts/analysis/generate_workspace_readmes.py [--dry-run]
"""

from __future__ import annotations

import argparse
import json
import os
import re
from datetime import datetime
from pathlib import Path
from typing import Optional

WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces"
)

# Workspace status classification based on naming suffixes
SUFFIX_STATUS = {
    "_before_hp_fix": ("superseded", "已被 HP integer tag fix 後版本取代"),
    "_smoke": ("smoketest", "測試用工作區"),
    "_smoketest_": ("smoketest", "測試用工作區"),
}

# Known workspace metadata for workspaces without round_context.json
KNOWN_METADATA: dict[str, dict] = {
    "20260317_cross_sample_methylation_observation_workspace": {
        "title": "Cross-sample methylation observation",
        "status": "historical",
        "note": "早期跨樣本甲基化觀察，已被 O-series 取代",
    },
    "20260317_quadrant_analysis_workspace": {
        "title": "Quadrant analysis",
        "status": "historical",
        "note": "早期象限分析，已被 O01 取代",
    },
    "20260318_cross_sample_quadrant_methylation_workspace": {
        "title": "Cross-sample quadrant methylation analysis",
        "status": "historical",
        "note": "早期跨樣本象限甲基化，已被 O-series 取代",
    },
    "20260318_per_sample_tp_fp_analysis": {
        "title": "Per-sample TP/FP analysis",
        "status": "historical",
        "note": "早期每樣本 TP/FP 分析",
    },
    "20260319_to_fp_normal_solved_analysis": {
        "title": "TO FP normal-solved analysis",
        "status": "historical",
        "note": "TO FP 以 normal BAM 解決的分析",
    },
    "20260324_methylation_cluster_allele_cooccurrence_HCC1395": {
        "title": "Methylation-cluster-allele co-occurrence (HCC1395)",
        "status": "historical",
        "note": "HCC1395 甲基化-分群-等位基因共現分析",
    },
}


def parse_date_from_name(name: str) -> str:
    """Extract YYYYMMDD from workspace name and format as YYYY-MM-DD."""
    m = re.match(r"(\d{8})_", name)
    if m:
        d = m.group(1)
        return f"{d[:4]}-{d[4:6]}-{d[6:8]}"
    return "unknown"


def classify_status(name: str) -> tuple[str, str]:
    """Determine workspace status from name suffixes."""
    for suffix, (status, note) in SUFFIX_STATUS.items():
        if suffix in name:
            return status, note
    if name in KNOWN_METADATA:
        return KNOWN_METADATA[name].get("status", "active"), ""
    return "active", ""


def find_generation_script(name: str) -> Optional[str]:
    """Try to find the Python script that generated this workspace."""
    scripts_dir = Path(__file__).resolve().parent
    # Try common prefixes
    candidates = []

    # O-series
    m = re.search(r"_O(\d+[a-z]?)(?:v\d)?_", name)
    if m:
        obs_id = m.group(1)
        candidates.extend(scripts_dir.glob(f"build_observation_O{obs_id}*.py"))

    # LOH round
    m = re.search(r"_loh_round(\d+)_", name)
    if m:
        candidates.extend(scripts_dir.glob(f"build_loh_round{m.group(1)}*.py"))

    # Germline FP
    if "germline_fp" in name:
        candidates.extend(scripts_dir.glob("build_germline_fp_*.py"))

    # Generic topic match
    topic_parts = name.split("_")[1:]  # skip date
    for part in topic_parts:
        if len(part) > 4 and part not in ("before", "after", "post", "analysis", "workspace"):
            candidates.extend(scripts_dir.glob(f"build_*{part}*.py"))

    if candidates:
        # Deduplicate and return first match
        seen = set()
        unique = []
        for c in candidates:
            if c.name not in seen:
                seen.add(c.name)
                unique.append(c)
        return f"scripts/analysis/{unique[0].name}"
    return None


def list_workspace_contents(ws_path: Path) -> dict:
    """Summarize workspace contents."""
    figures = sorted(ws_path.glob("figures/*.png")) if (ws_path / "figures").is_dir() else []
    png_root = sorted(ws_path.glob("*.png"))
    all_figs = figures + png_root

    data_files = []
    for ext in ("*.tsv", "*.csv", "*.tsv.gz"):
        data_files.extend(ws_path.glob(ext))
        if (ws_path / "data").is_dir():
            data_files.extend((ws_path / "data").glob(ext))
    data_files = sorted(set(data_files))

    md_files = sorted(
        f for f in ws_path.glob("*.md")
        if f.name not in ("README.md",)
    )

    has_context = (ws_path / "round_context.json").exists()
    has_report = any(
        (ws_path / name).exists()
        for name in ("observation_report.md", "round_summary.md", "觀察.md")
    )

    return {
        "figure_count": len(all_figs),
        "data_file_count": len(data_files),
        "md_file_count": len(md_files),
        "has_context": has_context,
        "has_report": has_report,
        "key_files": [f.name for f in data_files[:10]],
        "figure_names": [f.name for f in all_figs[:8]],
    }


def load_round_context(ws_path: Path) -> dict:
    """Load round_context.json if it exists."""
    ctx_path = ws_path / "round_context.json"
    if ctx_path.exists():
        try:
            with open(ctx_path) as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError):
            pass
    return {}


def generate_readme(ws_path: Path) -> str:
    """Generate README.md content for a workspace."""
    name = ws_path.name
    date_str = parse_date_from_name(name)
    status, status_note = classify_status(name)
    ctx = load_round_context(ws_path)
    contents = list_workspace_contents(ws_path)
    known = KNOWN_METADATA.get(name, {})

    # Title
    title = ctx.get("title") or known.get("title") or name
    obs_id = ctx.get("observation_id", "")

    # Script
    script = ctx.get("script") or find_generation_script(name)
    if script and "/" in script:
        script = script.split("/InterSubMod/")[-1] if "/InterSubMod/" in script else script

    # Description
    description = ctx.get("description") or known.get("note") or ""

    # Build README
    lines = [f"# {title}"]
    lines.append("")
    lines.append(f"- **建立日期**: {date_str}")
    lines.append(f"- **狀態**: {status}")
    if obs_id:
        lines.append(f"- **觀察 ID**: {obs_id}")
    if script:
        lines.append(f"- **生成腳本**: `{script}`")
    lines.append(f"- **輸入資料**: Master dataset (`all_region_rows.tsv.gz`)")

    # Status warning
    if status == "superseded":
        lines.append("")
        lines.append(f"> **警告**: {status_note or '此工作區已被取代，不應引用其結果。'}")
        # Try to find the replacement
        base = name.replace("_before_hp_fix", "")
        post = base + "_post_to_hp_fix"
        if (ws_path.parent / post).exists():
            lines.append(f"> **替代版本**: `{post}`")
        elif (ws_path.parent / base).exists() and base != name:
            lines.append(f"> **替代版本**: `{base}`")
    elif status == "smoketest":
        lines.append("")
        lines.append("> **注意**: 此為測試工作區，數據不完整，不應引用。")
    elif status == "historical":
        lines.append("")
        lines.append(f"> **注意**: {known.get('note', '此為早期探索工作區，已被後續系統性觀察取代。')}")

    if description:
        lines.append("")
        lines.append(f"## 說明")
        lines.append("")
        lines.append(description)

    # Contents summary
    lines.append("")
    lines.append("## 產出摘要")
    lines.append("")
    lines.append(f"| 項目 | 數量 |")
    lines.append(f"|------|------|")
    lines.append(f"| 圖表 (PNG) | {contents['figure_count']} |")
    lines.append(f"| 數據表 (TSV/CSV) | {contents['data_file_count']} |")
    lines.append(f"| 報告 (MD) | {contents['md_file_count']} |")

    if contents["key_files"]:
        lines.append("")
        lines.append("## 主要檔案")
        lines.append("")
        for f in contents["key_files"]:
            lines.append(f"- `{f}`")

    # Metadata from round_context
    if ctx.get("samples"):
        lines.append("")
        lines.append(f"## 涵蓋樣本")
        lines.append("")
        lines.append(", ".join(ctx["samples"]))

    if ctx.get("modes"):
        lines.append(f"")
        lines.append(f"**模式**: {', '.join(ctx['modes'])}")

    lines.append("")
    lines.append("---")
    lines.append(f"*自動生成於 {datetime.now().strftime('%Y-%m-%d')} by `generate_workspace_readmes.py`*")
    lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Generate workspace READMEs")
    parser.add_argument("--dry-run", action="store_true", help="Print without writing")
    parser.add_argument("--force", action="store_true", help="Overwrite existing READMEs")
    args = parser.parse_args()

    if not WORKSPACE_ROOT.is_dir():
        print(f"ERROR: Workspace root not found: {WORKSPACE_ROOT}")
        return 1

    generated = 0
    skipped = 0
    for ws_dir in sorted(WORKSPACE_ROOT.iterdir()):
        if not ws_dir.is_dir() or ws_dir.name.startswith("."):
            continue
        readme_path = ws_dir / "README.md"
        if readme_path.exists() and not args.force:
            skipped += 1
            continue

        content = generate_readme(ws_dir)
        if args.dry_run:
            print(f"\n{'='*60}")
            print(f"[DRY RUN] Would write: {readme_path}")
            print(f"{'='*60}")
            print(content[:500])
            print("...")
        else:
            readme_path.write_text(content, encoding="utf-8")
            print(f"  Generated: {ws_dir.name}/README.md")
        generated += 1

    print(f"\nSummary: {generated} generated, {skipped} skipped (already have README)")
    return 0


if __name__ == "__main__":
    exit(main())
