#!/usr/bin/env python3
"""從既有已驗證 run 抽取 canonical 數字，凍結成回歸比對用的 baseline。

原則：所有數字由本腳本從真實檔案讀出，**絕不手打**。
缺任何來源檔即 exit 2，不產出半套 baseline。

用法:
    python3 pipeline/lineage/freeze_baseline.py                    # 產生 baseline
    python3 pipeline/lineage/freeze_baseline.py --check            # 只驗證來源可讀
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BIG7 = Path("/big7_disk/liaoyoyo2001/big7_disk_output")

SCHEMA_NAME = "intersubmod.lineage_pipeline_baseline"
SCHEMA_VERSION = "1.0.0"

SOURCES = {
    "strict_linkage": REPO
    / "research/20260723_production_exact_ps_strict_read_linkage"
    / "20260723_exactPS嚴格ReadLinkage全資料集報告_01/data/all7_report_data.json",
    "topology": BIG7
    / "synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples"
    / "all7_strict_guard1000_v1/summary/all7_summary.json",
}

DATASETS = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]

# 從 strict_linkage 抽的 per-sample 欄位
LINKAGE_FIELDS = [
    "candidate_loci_S",
    "canonical_molecule_rows",
    "all_components",
    "active_node_memberships",
    "active_unique_loci",
    "HP1_W",
    "HP2_W",
    "W_k_gt12",
    "W_memberships",
    "W_span_gt_50kb",
]

# 從 topology 抽的 per-sample 欄位
TOPOLOGY_FIELDS = [
    "groups_total",
    "mutation_bearing_units",
    "mutation_family_complete_units",
    "mutation_family_abstain_units",
    "objective_certified_units",
    "objective_abstain_units",
    "no_active_alt_units",
    "ranked_units",
    "best_tree_unique_units",
    "best_tree_tied_units",
    "recurrence_required_units",
]

TOPOLOGY_TOTALS = [
    "groups_total",
    "mutation_bearing_units",
    "mutation_family_complete_units",
    "mutation_family_abstain_units",
    "objective_certified_units",
    "objective_abstain_units",
    "no_active_alt_units",
    "ranked_units",
    "best_tree_unique_units",
    "best_tree_tied_units",
    "recurrence_required_units",
    "zero_denominator_units",
    "topology_runtime_seconds",
]


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_sources() -> dict[str, dict]:
    missing = [name for name, p in SOURCES.items() if not p.is_file()]
    if missing:
        for name in missing:
            print(f"MISSING SOURCE: {name} -> {SOURCES[name]}", file=sys.stderr)
        sys.exit(2)
    return {name: json.loads(p.read_text(encoding="utf-8")) for name, p in SOURCES.items()}


def build_baseline(raw: dict[str, dict]) -> dict:
    linkage = {e["dataset"]: e for e in raw["strict_linkage"].get("datasets", [])}
    topo = {s["sample"]: s for s in raw["topology"].get("samples", [])}

    absent = [d for d in DATASETS if d not in linkage or d not in topo]
    if absent:
        print(f"MISSING DATASETS in sources: {absent}", file=sys.stderr)
        sys.exit(2)

    samples = {}
    for name in DATASETS:
        lk, tp = linkage[name], topo[name]
        entry = {f: lk.get(f) for f in LINKAGE_FIELDS}
        entry.update({f: tp.get(f) for f in TOPOLOGY_FIELDS})
        entry["W_total"] = (lk.get("HP1_W") or 0) + (lk.get("HP2_W") or 0)
        entry["active_k_distribution"] = tp.get("active_k_distribution", {})
        el = tp.get("solver_elapsed_microseconds", {})
        entry["solver_elapsed_us_sum"] = el.get("sum")
        entry["solver_elapsed_us_max"] = el.get("max")
        entry["search_nodes_sum"] = (tp.get("search_nodes") or {}).get("sum")
        samples[name] = entry

    totals_src = raw["topology"].get("totals", {})
    return {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "provenance": {
            name: {
                "path": str(p),
                "sha256": sha256_of(p),
                "size_bytes": p.stat().st_size,
            }
            for name, p in SOURCES.items()
        },
        "source_all_pass": {
            "strict_linkage": raw["strict_linkage"].get("all_pass"),
            "topology": raw["topology"].get("all_pass"),
        },
        "datasets": DATASETS,
        "samples": samples,
        "cohort_totals": {f: totals_src.get(f) for f in TOPOLOGY_TOTALS},
        "cohort_W_total": sum(samples[d]["W_total"] for d in DATASETS),
        "guards": (topo["HCC1395"].get("configured_guards") or {}),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="只驗證來源可讀，不寫檔")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "baseline" / "cohort_baseline.json",
    )
    args = ap.parse_args()

    raw = load_sources()
    if args.check:
        for name, p in SOURCES.items():
            print(f"OK  {name}: {p} ({p.stat().st_size} bytes)")
        return 0

    baseline = build_baseline(raw)

    if not all(baseline["source_all_pass"].values()):
        print(f"REFUSE: source all_pass is not true: {baseline['source_all_pass']}", file=sys.stderr)
        return 3

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(baseline, ensure_ascii=False, indent=1, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"baseline written: {args.output}")
    print(f"  datasets      : {len(baseline['datasets'])}")
    print(f"  cohort W      : {baseline['cohort_W_total']}")
    print(f"  cohort groups : {baseline['cohort_totals']['groups_total']}")
    print(f"  runtime (s)   : {baseline['cohort_totals']['topology_runtime_seconds']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
