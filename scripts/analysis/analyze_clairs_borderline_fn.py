#!/usr/bin/env python3
"""Analyze ClairS borderline FN variants (Pool B) and evaluate rescue rules.

Pool B = ClairS non-PASS variants (excluding NonSomatic) that overlap the
truth set. These are variants ClairS generated but filtered out, yet they
are true somatic mutations (FN from caller perspective).

Outputs:
  pool_b_fn_stats.tsv             - Overall Pool B FN count and baseline metrics
  pool_b_fn_distribution.tsv      - GQ/AF/FILTER distribution by group
  pool_b_caller_rescue_rules.tsv  - Caller-only rules F1 evaluation
  pool_b_combined_rescue_rules.tsv - Caller+methyl combined rules (requires --significance-csv)
  pool_b_joined_features.tsv      - Joined feature table (caller + methyl features)
  pool_b_summary.md               - Human-readable summary
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

from research_common import compute_metrics, infer_platform, to_float, to_int, write_tsv_rows


RuleFunc = Callable[[Dict[str, object]], bool]

CALLER_RULE_FIELDS = [
    "rule_id",
    "description",
    "tp_rescued",
    "fp_reintroduced",
    "rescue_precision",
    "new_tp",
    "new_fp",
    "precision",
    "recall",
    "f1",
    "delta_f1",
    "baseline_f1",
]

COMBINED_RULE_FIELDS = CALLER_RULE_FIELDS  # same schema

DISTRIBUTION_FIELDS = [
    "group",
    "count",
    "median_gq",
    "median_af",
    "median_dp",
    "filter_lowqual_only_frac",
    "filter_varcluster_frac",
    "filter_strandbias_frac",
    "filter_nonsomatic_frac",
]

STATS_FIELDS = [
    "key",
    "value",
]

JOINED_FIELDS = [
    "region_key",
    "group",
    "gq",
    "af",
    "dp",
    "filter",
    "filter_tags",
    "VerificationClass",
    "AlleleDelta",
    "CramersV",
    "PairwiseMedianDist",
    "HPMergedDelta",
    "HPMergedSig",
    "NumReads",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pool-b-fn-vcf", required=True,
                   help="Pool B FN VCF: non-PASS variants that ARE in truth set")
    p.add_argument("--pool-b-non-fn-vcf", required=True,
                   help="Pool B non-FN VCF: non-PASS variants NOT in truth set (true FP/germline)")
    p.add_argument("--significance-csv",
                   help="InterSubMod significance_summary.csv for Pool B FN (intersubmod_tp/)")
    p.add_argument("--significance-csv-non-fn",
                   help="InterSubMod significance_summary.csv for Pool B non-FN (intersubmod_fp/); "
                        "enables accurate combined rule FP estimation")
    p.add_argument("--baseline-tp", type=int, required=True,
                   help="ClairS-TO PASS baseline TP count")
    p.add_argument("--baseline-fp", type=int, required=True,
                   help="ClairS-TO PASS baseline FP count")
    p.add_argument("--truth-total", type=int, required=True,
                   help="Truth total SNV count (for recall / F1)")
    p.add_argument("--sample", default="HCC1395", help="Sample name")
    p.add_argument("--output-dir", required=True, help="Output directory")
    return p.parse_args()


def _open_vcf(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_vcf(path: Path) -> List[Dict[str, object]]:
    """Parse VCF; return list of variant dicts with caller-level fields."""
    rows: List[Dict[str, object]] = []
    with _open_vcf(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _vid, ref, alt, qual_str, filt, info = parts[:8]
            rest = parts[8:]

            info_map: Dict[str, str] = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_map[k] = v

            fmt_map: Dict[str, str] = {}
            if len(rest) >= 2:
                fmt_map = dict(zip(rest[0].split(":"), rest[1].split(":")))

            # AF: prefer FORMAT AF/VAF, fallback to INFO
            af = math.nan
            for key in ("AF", "VAF"):
                if key in fmt_map:
                    af = to_float(fmt_map[key].split(",")[0])
                    break
                if key in info_map:
                    af = to_float(info_map[key].split(",")[0])
                    break

            dp = to_int(fmt_map.get("DP", "0"))
            gq = to_int(fmt_map.get("GQ", "0"))

            # filter tags as frozenset
            filter_tags = frozenset(t for t in filt.split(";") if t and t not in {".", "PASS"})

            rows.append({
                "region_key": f"{chrom}:{pos}:{ref}:{alt}",
                "gq": gq,
                "af": af,
                "dp": dp,
                "filter": filt,
                "filter_tags": filter_tags,
                "qual": to_float(qual_str),
            })
    return rows


def load_significance(path: Path) -> Dict[str, Dict[str, object]]:
    """Load significance_summary.csv; key = Chr:Pos:Ref:Alt."""
    result: Dict[str, Dict[str, object]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = f"{row['Chr']}:{row['Pos']}:{row['Ref']}:{row['Alt']}"
            result[key] = row
    return result


def _median(values: Sequence[float]) -> str:
    clean = sorted(v for v in values if not math.isnan(v))
    if not clean:
        return ""
    n = len(clean)
    mid = n // 2
    if n % 2:
        return f"{clean[mid]:.4f}"
    return f"{((clean[mid - 1] + clean[mid]) / 2.0):.4f}"


def build_distribution_row(group: str, variants: List[Dict]) -> Dict:
    if not variants:
        return {f: "" for f in DISTRIBUTION_FIELDS} | {"group": group, "count": 0}
    gqs = [float(v["gq"]) for v in variants]
    afs = [float(v["af"]) for v in variants if not math.isnan(float(v["af"]))]
    dps = [float(v["dp"]) for v in variants]
    n = len(variants)
    lowqual_only = sum(1 for v in variants if v["filter_tags"] == frozenset({"LowQual"}))
    varcluster = sum(1 for v in variants if "VariantCluster" in v["filter_tags"])
    strandbias = sum(1 for v in variants if "StrandBias" in v["filter_tags"])
    nonsomatic = sum(1 for v in variants if "NonSomatic" in v["filter_tags"])
    return {
        "group": group,
        "count": n,
        "median_gq": _median(gqs),
        "median_af": _median(afs),
        "median_dp": _median(dps),
        "filter_lowqual_only_frac": f"{lowqual_only / n:.4f}",
        "filter_varcluster_frac": f"{varcluster / n:.4f}",
        "filter_strandbias_frac": f"{strandbias / n:.4f}",
        "filter_nonsomatic_frac": f"{nonsomatic / n:.4f}",
    }


def _build_caller_rules() -> Dict[str, tuple[str, RuleFunc]]:
    """Returns dict of rule_id -> (description, callable)."""
    return {
        "gq_ge_15": ("GQ >= 15", lambda v: v["gq"] >= 15),
        "gq_ge_20": ("GQ >= 20", lambda v: v["gq"] >= 20),
        "gq_ge_25": ("GQ >= 25", lambda v: v["gq"] >= 25),
        "af_ge_0_10": ("AF >= 0.10", lambda v: not math.isnan(float(v["af"])) and float(v["af"]) >= 0.10),
        "af_ge_0_15": ("AF >= 0.15", lambda v: not math.isnan(float(v["af"])) and float(v["af"]) >= 0.15),
        "af_ge_0_20": ("AF >= 0.20", lambda v: not math.isnan(float(v["af"])) and float(v["af"]) >= 0.20),
        "gq15_and_af015": (
            "GQ >= 15 AND AF >= 0.15",
            lambda v: v["gq"] >= 15 and not math.isnan(float(v["af"])) and float(v["af"]) >= 0.15,
        ),
        "lowqual_only": (
            "FILTER = LowQual only (no other tags)",
            lambda v: v["filter_tags"] == frozenset({"LowQual"}),
        ),
        "lowqual_only_and_gq15": (
            "FILTER = LowQual only AND GQ >= 15",
            lambda v: v["filter_tags"] == frozenset({"LowQual"}) and v["gq"] >= 15,
        ),
        "no_varcluster": (
            "VariantCluster NOT in FILTER",
            lambda v: "VariantCluster" not in v["filter_tags"],
        ),
        "no_varcluster_and_gq15": (
            "VariantCluster NOT in FILTER AND GQ >= 15",
            lambda v: "VariantCluster" not in v["filter_tags"] and v["gq"] >= 15,
        ),
        "no_strandbias_and_gq15": (
            "StrandBias NOT in FILTER AND GQ >= 15",
            lambda v: "StrandBias" not in v["filter_tags"] and v["gq"] >= 15,
        ),
    }


def _build_combined_rules() -> Dict[str, tuple[str, RuleFunc]]:
    """Combined caller + methylation rescue rules.

    Variants without methylation data get VerificationClass='NA'.
    """
    METHYL_CLASSES_STRONG = {"Strong", "Subclone"}

    def _af(v):
        return float(v["af"]) if not math.isnan(float(v["af"])) else 0.0

    def _pd(v):
        return float(v.get("PairwiseMedianDist") or 0.0)

    def _ad(v):
        return float(v.get("AlleleDelta") or 0.0)

    def _cls(v):
        return v.get("VerificationClass", "NA")

    def _hp_sig(v):
        return str(v.get("HPMergedSig", "")).strip().lower() in {"true", "1", "yes"}

    return {
        # Caller-only baselines
        "gq_ge_20": ("GQ >= 20 (caller only)", lambda v: v["gq"] >= 20),
        "lowqual_only_gq15": (
            "FILTER=LowQual only AND GQ >= 15 (caller only)",
            lambda v: v["filter_tags"] == frozenset({"LowQual"}) and v["gq"] >= 15,
        ),
        # Methyl-only
        "class_strong_or_subclone": (
            "VerificationClass in {Strong, Subclone}",
            lambda v: _cls(v) in METHYL_CLASSES_STRONG,
        ),
        "pairwise_ge_020": (
            "PairwiseMedianDist >= 0.20",
            lambda v: _pd(v) >= 0.20,
        ),
        # Combined
        "gq15_and_not_noise": (
            "GQ >= 15 AND VerificationClass != Noise",
            lambda v: v["gq"] >= 15 and _cls(v) not in {"Noise", "NA"},
        ),
        "gq15_and_strong_or_subclone": (
            "GQ >= 15 AND VerificationClass in {Strong, Subclone}",
            lambda v: v["gq"] >= 15 and _cls(v) in METHYL_CLASSES_STRONG,
        ),
        "gq15_and_allele_delta_low": (
            "GQ >= 15 AND AlleleDelta < 0.15 (low AD = not artifact-like)",
            lambda v: v["gq"] >= 15 and _ad(v) < 0.15,
        ),
        "lowqual_and_pairwise_ge_015": (
            "FILTER=LowQual only AND PairwiseMedianDist >= 0.15",
            lambda v: v["filter_tags"] == frozenset({"LowQual"}) and _pd(v) >= 0.15,
        ),
        "any_gq_and_hp_sig": (
            "GQ >= 10 AND HPMergedSig = True",
            lambda v: v["gq"] >= 10 and _hp_sig(v),
        ),
        "gq15_and_pairwise_ge_015": (
            "GQ >= 15 AND PairwiseMedianDist >= 0.15",
            lambda v: v["gq"] >= 15 and _pd(v) >= 0.15,
        ),
    }


def evaluate_rules(
    rules: Dict[str, tuple[str, RuleFunc]],
    fn_variants: List[Dict],
    non_fn_variants: List[Dict],
    baseline_tp: int,
    baseline_fp: int,
    truth_total: int,
    baseline_f1: float,
) -> List[Dict]:
    rows = []
    for rule_id, (desc, rule_fn) in rules.items():
        tp_rescued = sum(1 for v in fn_variants if rule_fn(v))
        fp_reintroduced = sum(1 for v in non_fn_variants if rule_fn(v))
        new_tp = baseline_tp + tp_rescued
        new_fp = baseline_fp + fp_reintroduced
        metrics = compute_metrics(new_tp, new_fp, truth_total)
        delta = metrics["f1"] - baseline_f1
        rescue_prec = tp_rescued / (tp_rescued + fp_reintroduced) if (tp_rescued + fp_reintroduced) else float("nan")
        rows.append({
            "rule_id": rule_id,
            "description": desc,
            "tp_rescued": tp_rescued,
            "fp_reintroduced": fp_reintroduced,
            "rescue_precision": f"{rescue_prec:.4f}" if not math.isnan(rescue_prec) else "",
            "new_tp": new_tp,
            "new_fp": new_fp,
            "precision": f"{metrics['precision']:.6f}",
            "recall": f"{metrics['recall']:.6f}",
            "f1": f"{metrics['f1']:.6f}",
            "delta_f1": f"{delta:.6f}",
            "baseline_f1": f"{baseline_f1:.6f}",
        })
    return sorted(rows, key=lambda r: float(r["delta_f1"]), reverse=True)


def join_methylation(
    variants: List[Dict],
    sig_map: Dict[str, Dict],
    group: str,
) -> List[Dict]:
    result = []
    for v in variants:
        key = v["region_key"]
        sig = sig_map.get(key, {})
        result.append({
            "region_key": key,
            "group": group,
            "gq": v["gq"],
            "af": f"{v['af']:.4f}" if not math.isnan(float(v["af"])) else "",
            "dp": v["dp"],
            "filter": v["filter"],
            "filter_tags": ";".join(sorted(v["filter_tags"])),
            "VerificationClass": sig.get("VerificationClass", ""),
            "AlleleDelta": sig.get("AlleleDelta", ""),
            "CramersV": sig.get("CramersV", ""),
            "PairwiseMedianDist": sig.get("PairwiseMedianDist", ""),
            "HPMergedDelta": sig.get("HPMergedDelta", ""),
            "HPMergedSig": sig.get("HPMergedSig", ""),
            "NumReads": sig.get("NumReads", ""),
        })
    return result


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    pool_b_fn = parse_vcf(Path(args.pool_b_fn_vcf))
    pool_b_non_fn = parse_vcf(Path(args.pool_b_non_fn_vcf))

    baseline_metrics = compute_metrics(args.baseline_tp, args.baseline_fp, args.truth_total)
    baseline_f1 = baseline_metrics["f1"]

    # Load methylation data if provided
    sig_map: Dict[str, Dict] = {}       # for Pool B FN (TP side)
    sig_map_non_fn: Dict[str, Dict] = {}  # for Pool B non-FN (FP side)
    if args.significance_csv:
        sig_path = Path(args.significance_csv)
        if sig_path.exists():
            sig_map = load_significance(sig_path)
            print(f"[analyze_clairs_borderline_fn] Loaded {len(sig_map)} regions from significance CSV (FN)")
        else:
            print(f"[analyze_clairs_borderline_fn] WARNING: significance CSV not found: {sig_path}")
    if args.significance_csv_non_fn:
        sig_path_nfn = Path(args.significance_csv_non_fn)
        if sig_path_nfn.exists():
            sig_map_non_fn = load_significance(sig_path_nfn)
            print(f"[analyze_clairs_borderline_fn] Loaded {len(sig_map_non_fn)} regions from significance CSV (non-FN)")
        else:
            print(f"[analyze_clairs_borderline_fn] WARNING: non-FN significance CSV not found: {sig_path_nfn}")

    fn_keys = {v["region_key"] for v in pool_b_fn}
    fn_with_methyl = [v for v in pool_b_fn if v["region_key"] in sig_map]
    fn_without_methyl = [v for v in pool_b_fn if v["region_key"] not in sig_map]
    methyl_coverage = len(fn_with_methyl) / len(pool_b_fn) if pool_b_fn else 0.0
    non_fn_with_methyl = [v for v in pool_b_non_fn if v["region_key"] in sig_map_non_fn]
    non_fn_methyl_coverage = len(non_fn_with_methyl) / len(pool_b_non_fn) if pool_b_non_fn else 0.0

    # Stats file
    stats_rows = [
        {"key": "sample", "value": args.sample},
        {"key": "platform", "value": infer_platform(args.sample)},
        {"key": "pool_b_fn_count", "value": len(pool_b_fn)},
        {"key": "pool_b_non_fn_count", "value": len(pool_b_non_fn)},
        {"key": "pool_b_total", "value": len(pool_b_fn) + len(pool_b_non_fn)},
        {"key": "fn_with_methyl_coverage", "value": f"{methyl_coverage:.4f}"},
        {"key": "fn_with_methyl_count", "value": len(fn_with_methyl)},
        {"key": "non_fn_with_methyl_coverage", "value": f"{non_fn_methyl_coverage:.4f}"},
        {"key": "non_fn_with_methyl_count", "value": len(non_fn_with_methyl)},
        {"key": "baseline_tp", "value": args.baseline_tp},
        {"key": "baseline_fp", "value": args.baseline_fp},
        {"key": "truth_total", "value": args.truth_total},
        {"key": "baseline_precision", "value": f"{baseline_metrics['precision']:.6f}"},
        {"key": "baseline_recall", "value": f"{baseline_metrics['recall']:.6f}"},
        {"key": "baseline_f1", "value": f"{baseline_f1:.6f}"},
        {"key": "go_nogo_threshold_fn", "value": "30"},
        {"key": "go_nogo_result", "value": "GO" if len(pool_b_fn) >= 30 else "NO-GO (Pool B FN < 30)"},
        {"key": "methyl_coverage_threshold", "value": "0.40"},
        {"key": "methyl_coverage_sufficient", "value": "YES" if methyl_coverage >= 0.40 else "NO"},
    ]
    write_tsv_rows(output_dir / "pool_b_fn_stats.tsv", STATS_FIELDS, stats_rows)

    # Distribution
    dist_rows = [
        build_distribution_row("pool_b_fn", pool_b_fn),
        build_distribution_row("pool_b_non_fn", pool_b_non_fn),
        build_distribution_row("pool_b_fn_with_methyl", fn_with_methyl),
        build_distribution_row("pool_b_fn_without_methyl", fn_without_methyl),
    ]
    write_tsv_rows(output_dir / "pool_b_fn_distribution.tsv", DISTRIBUTION_FIELDS, dist_rows)

    # Caller-only rescue rules
    caller_rules = _build_caller_rules()
    caller_rule_rows = evaluate_rules(
        caller_rules, pool_b_fn, pool_b_non_fn,
        args.baseline_tp, args.baseline_fp, args.truth_total, baseline_f1,
    )
    write_tsv_rows(output_dir / "pool_b_caller_rescue_rules.tsv", CALLER_RULE_FIELDS, caller_rule_rows)

    # Joined features table (all pool_b_fn + pool_b_non_fn with methyl)
    non_fn_sig = sig_map_non_fn if sig_map_non_fn else sig_map
    joined_rows = (
        join_methylation(pool_b_fn, sig_map, "pool_b_fn")
        + join_methylation(pool_b_non_fn, non_fn_sig, "pool_b_non_fn")
    )
    write_tsv_rows(output_dir / "pool_b_joined_features.tsv", JOINED_FIELDS, joined_rows)

    # Combined rules (only if methylation available for at least FN side)
    combined_rule_rows = []
    if sig_map:
        def _enrich(variants: List[Dict], lookup: Dict[str, Dict]) -> List[Dict]:
            enriched = []
            for v in variants:
                sig = lookup.get(v["region_key"], {})
                enriched.append({
                    **v,
                    "VerificationClass": sig.get("VerificationClass", "NA"),
                    "AlleleDelta": to_float(sig.get("AlleleDelta"), 0.0),
                    "CramersV": to_float(sig.get("CramersV"), 0.0),
                    "PairwiseMedianDist": to_float(sig.get("PairwiseMedianDist"), 0.0),
                    "HPMergedDelta": to_float(sig.get("HPMergedDelta"), 0.0),
                    "HPMergedSig": sig.get("HPMergedSig", ""),
                })
            return enriched

        # Use dedicated sig_map for each group; fall back to sig_map if non-fn not provided
        non_fn_lookup = sig_map_non_fn if sig_map_non_fn else sig_map
        combined_rules = _build_combined_rules()
        combined_rule_rows = evaluate_rules(
            combined_rules,
            _enrich(pool_b_fn, sig_map),
            _enrich(pool_b_non_fn, non_fn_lookup),
            args.baseline_tp, args.baseline_fp, args.truth_total, baseline_f1,
        )
        if not sig_map_non_fn:
            print("[analyze_clairs_borderline_fn] WARNING: no --significance-csv-non-fn provided; "
                  "combined rule FP estimates may be underestimated")
        write_tsv_rows(output_dir / "pool_b_combined_rescue_rules.tsv", COMBINED_RULE_FIELDS, combined_rule_rows)

    # Markdown summary
    best_caller = caller_rule_rows[0] if caller_rule_rows else {}
    best_combined = combined_rule_rows[0] if combined_rule_rows else {}

    md = [
        "# Pool B FN Analysis Summary",
        "",
        f"- Sample: `{args.sample}` / Platform: `{infer_platform(args.sample)}`",
        f"- Truth total: `{args.truth_total}`",
        f"- Baseline (ClairS-TO PASS): TP=`{args.baseline_tp}`, FP=`{args.baseline_fp}`, F1=`{baseline_f1:.6f}`",
        "",
        "## Pool B FN 數量（Go/No-Go）",
        "",
        f"- Pool B FN: **{len(pool_b_fn)}** 個",
        f"- Pool B non-FN (過濾的 FP/germline): {len(pool_b_non_fn)} 個",
        f"- 甲基化覆蓋率: {methyl_coverage:.1%} ({len(fn_with_methyl)}/{len(pool_b_fn)})",
        f"- **Go/No-Go**: {'GO → 繼續實驗 2-4' if len(pool_b_fn) >= 30 else 'NO-GO → Pool B FN < 30，考慮 Pool A-light'}",
        "",
        "## 最佳 Caller-Only 救援規則",
        "",
    ]
    if best_caller:
        md += [
            f"- Rule: `{best_caller['rule_id']}`",
            f"- Description: {best_caller['description']}",
            f"- TP rescued: {best_caller['tp_rescued']} / FP reintroduced: {best_caller['fp_reintroduced']}",
            f"- Rescue precision: {best_caller['rescue_precision']}",
            f"- F1: {best_caller['f1']} (delta: {best_caller['delta_f1']})",
        ]
    if best_combined:
        md += [
            "",
            "## 最佳 Combined (Caller + Methyl) 救援規則",
            "",
            f"- Rule: `{best_combined['rule_id']}`",
            f"- Description: {best_combined['description']}",
            f"- TP rescued: {best_combined['tp_rescued']} / FP reintroduced: {best_combined['fp_reintroduced']}",
            f"- Rescue precision: {best_combined['rescue_precision']}",
            f"- F1: {best_combined['f1']} (delta: {best_combined['delta_f1']})",
        ]
    md += [
        "",
        "## 輸出檔案",
        "",
        "- `pool_b_fn_stats.tsv` - 整體統計與 Go/No-Go 判斷",
        "- `pool_b_fn_distribution.tsv` - GQ/AF/FILTER 分布比較",
        "- `pool_b_caller_rescue_rules.tsv` - Caller-only 規則 F1 評估",
        "- `pool_b_combined_rescue_rules.tsv` - Caller+Methyl 聯合規則 F1 評估",
        "- `pool_b_joined_features.tsv` - 完整特徵表（含甲基欄位）",
    ]
    (output_dir / "pool_b_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print(f"[analyze_clairs_borderline_fn] Pool B FN count: {len(pool_b_fn)}")
    print(f"[analyze_clairs_borderline_fn] Pool B non-FN count: {len(pool_b_non_fn)}")
    print(f"[analyze_clairs_borderline_fn] Methyl coverage: {methyl_coverage:.1%}")
    print(f"[analyze_clairs_borderline_fn] Go/No-Go: {'GO' if len(pool_b_fn) >= 30 else 'NO-GO'}")
    print(f"[analyze_clairs_borderline_fn] Results written to: {output_dir}")


if __name__ == "__main__":
    main()
