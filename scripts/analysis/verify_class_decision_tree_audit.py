#!/usr/bin/env python3
"""
verify_class_decision_tree_audit.py
VerificationClass decision tree quantitative audit across 7 samples.

Covers:
  Step 2: VerificationClass decision tree flow (TP/FP counts per branch)
  Step 4: Significant field vs VerificationClass comparison
  Step 5: QualityScore penalty audit (LOH=-25, coverage anomalies)
  Step 6: LOH vs TP/FP distribution analysis
  Step 7: New feature discriminability (LOH_HP_Signal, HP_Coverage_Symmetry,
          PairwiseMedianDist stratified by LOH)
  Step 8: PERMANOVA penalty in HeuristicScore (x0.7 if ClusterPermanovaP>0.05)

Output:
  - Console: all tables
  - docs/methodology/20260323_VerificationClass決策樹跨樣本量化_01.md
"""

import argparse
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    CURRENT_CLASSES_V2,
    SchemaContractError,
    read_evidence,
    select_current_view,
)

# ── Data Paths ──────────────────────────────────────────────────────────────

BASE = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical"

SAMPLE_PATHS = {
    "HCC1395":        "20260314_HCC1395_paired_full_full_complete_matrix",
    "HCC1395_DORADO": "20260315_HCC1395_DORADO_paired_full_full_complete_matrix",
    "COLO829":        "20260315_COLO829_paired_full_full_complete_matrix",
    "H1437":          "20260315_H1437_paired_full_full_complete_matrix",
    "H2009":          "20260315_H2009_paired_full_full_complete_matrix",
    "HCC1937":        "20260315_HCC1937_paired_full_full_complete_matrix",
    "HCC1954":        "20260315_HCC1954_paired_full_full_complete_matrix",
}

# Also include HCC1395 pileup n=99 as a reference (large FP set)
EXTRA_PATHS = {
    "HCC1395_pileup_n99": {
        "tp": "/tmp/n99_full/tp/significance_summary.csv",
        "fp": "/tmp/n99_full/fp/significance_summary.csv",
    }
}

OUTPUT_DOC = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/20260323_VerificationClass決策樹跨樣本量化_01.md"

CURRENT_NOISE_CLASSES = (
    "Noise_Uniform",
    "Noise_Chaotic",
    "Noise_Uncorrelated",
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def safe_float(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def safe_bool(v):
    """Parse 'True'/'False'/'1'/'0' strings to bool."""
    return str(v).strip().lower() in ("true", "1", "yes")


def apply_verification_contract(frame):
    """Validate C2 + typed evidence; this audit forbids unknown current classes."""
    current = select_current_view(frame)
    if current.unknown_counts:
        raise SchemaContractError(
            f"decision-tree audit: unknown VerificationClass values: {current.unknown_counts}"
        )
    evidence = read_evidence(frame)
    selected = frame.copy()
    selected["_VerificationClass_C2"] = current.values
    selected["_LabelFirstSupport_E"] = evidence["LabelFirstSupport"]
    selected["_ClusterFirstSupport_E"] = evidence["ClusterFirstSupport"]
    selected["_EvidencePath_E"] = evidence["EvidencePath"]
    selected["_VerificationSourceField"] = current.field
    selected["_VerificationSchemaStatus"] = current.schema_status
    return selected, {
        "verification": current.metadata(),
        "evidence_fields": [
            "LabelFirstSupport",
            "ClusterFirstSupport",
            "WithinHPSupport",
            "DispersionWarning",
            "EvidencePath",
            "EvidenceDerivation",
        ],
    }


def load_csv(path):
    if not path or not os.path.exists(path):
        return []
    frame = pd.read_csv(path, low_memory=False)
    try:
        selected, _ = apply_verification_contract(frame)
    except SchemaContractError as exc:
        raise SchemaContractError(f"{path}: {exc}") from exc
    return [enrich_row(row) for row in selected.to_dict("records")]


def auroc(scores_tp, scores_fp):
    """Wilcoxon-Mann-Whitney AUROC (TP > FP direction)."""
    if not scores_tp or not scores_fp:
        return float("nan")
    n_tp = len(scores_tp)
    n_fp = len(scores_fp)
    combined = [(s, 1) for s in scores_tp] + [(s, 0) for s in scores_fp]
    combined.sort(key=lambda x: x[0])
    rank_sum = 0
    for rank, (_, label) in enumerate(combined, start=1):
        if label == 1:
            rank_sum += rank
    u = rank_sum - n_tp * (n_tp + 1) / 2
    return u / (n_tp * n_fp)


def precision(tp_count, fp_count):
    total = tp_count + fp_count
    return tp_count / total if total > 0 else float("nan")


def fmt_pct(n, total):
    if total == 0:
        return "0.0%"
    return f"{n/total:.1%}"


# ── Per-row Feature Computation ───────────────────────────────────────────────

def enrich_row(row):
    """Add derived features to a row dict."""
    hp1n = safe_float(row.get("HP1FamilyN", 0))
    hp2n = safe_float(row.get("HP2FamilyN", 0))
    max_hp = max(hp1n, hp2n)
    min_hp = min(hp1n, hp2n)
    row["_HP_Coverage_Symmetry"] = min_hp / (max_hp + 1)  # 0=asymmetric, 1=symmetric

    ad = abs(safe_float(row.get("AlleleDelta", 0)))
    hp_ratio = safe_float(row.get("HP_Ratio", 0.5))
    hp_imbalance = abs(hp_ratio - 0.5) * 2  # 0=balanced, 1=fully imbalanced
    loh = safe_bool(row.get("Potential_LOH", False))
    row["_HP_Imbalance"] = hp_imbalance
    row["_LOH_HP_Signal"] = hp_imbalance * ad if loh else ad

    row["_passed_gate"] = safe_bool(row.get("PassedGating", False))
    if "_VerificationClass_C2" not in row:
        raise SchemaContractError("decision-tree audit: row lacks validated C2 class")
    vc = row["_VerificationClass_C2"]
    row["_label_sig"] = bool(row["_LabelFirstSupport_E"])
    row["_cluster_sig"] = bool(row["_ClusterFirstSupport_E"])
    row["_vc"] = vc
    row["_significant"] = safe_bool(row.get("Significant", False))
    row["_loh"] = loh
    row["_quality_score"] = safe_float(row.get("Quality_Score", 0))
    row["_quality_tier"] = row.get("Quality_Tier", "")
    row["_perm_f"] = safe_float(row.get("ClusterPermanovaF", 0))
    row["_perm_p"] = safe_float(row.get("ClusterPermanovaP", 1))
    row["_perm_valid"] = safe_bool(row.get("ClusterPermanovaValid", False))
    row["_heuristic"] = safe_float(row.get("HeuristicScore", 0))
    row["_global_p"] = safe_float(row.get("GlobalP", 1))
    row["_cramers_v"] = safe_float(row.get("CramersV", 0))
    row["_allele_delta"] = ad
    row["_pairwise_mean"] = safe_float(row.get("PairwiseMeanDist", 0))
    row["_pairwise_median"] = safe_float(row.get("PairwiseMedianDist", 0))
    row["_suggest_filter"] = safe_bool(row.get("SuggestFilter", False))
    return row


# ── Decision Tree Flow Analysis ───────────────────────────────────────────────

def decision_tree_flow(tp_rows, fp_rows):
    """Compute TP/FP counts at each branch of the decision tree."""
    result = {}
    n_tp = len(tp_rows)
    n_fp = len(fp_rows)
    result["total"] = {"tp": n_tp, "fp": n_fp}

    # Gate split
    tp_pg = [r for r in tp_rows if r["_passed_gate"]]
    fp_pg = [r for r in fp_rows if r["_passed_gate"]]
    tp_ng = [r for r in tp_rows if not r["_passed_gate"]]
    fp_ng = [r for r in fp_rows if not r["_passed_gate"]]
    result["passed_gate_true"]  = {"tp": len(tp_pg), "fp": len(fp_pg)}
    result["passed_gate_false"] = {"tp": len(tp_ng), "fp": len(fp_ng)}

    # Canonical VerificationClass v2 branches
    result["verification_classes"] = {}
    for vc in CURRENT_CLASSES_V2:
        tp_vc = [r for r in tp_rows if r["_vc"] == vc]
        fp_vc = [r for r in fp_rows if r["_vc"] == vc]
        result["verification_classes"][vc] = {"tp": len(tp_vc), "fp": len(fp_vc)}

    # Evidence support comes from typed E fields, never from class names.
    for key, predicate in (
        ("cluster_first_support_true", lambda r: r["_cluster_sig"]),
        ("label_first_support_true", lambda r: r["_label_sig"]),
        ("bidirectional_support_true", lambda r: r["_cluster_sig"] and r["_label_sig"]),
    ):
        result[key] = {
            "tp": sum(1 for r in tp_rows if predicate(r)),
            "fp": sum(1 for r in fp_rows if predicate(r)),
        }

    # Significant field
    tp_sig = [r for r in tp_rows if r["_significant"]]
    fp_sig = [r for r in fp_rows if r["_significant"]]
    result["significant_true"]  = {"tp": len(tp_sig), "fp": len(fp_sig)}
    result["significant_false"] = {"tp": n_tp - len(tp_sig), "fp": n_fp - len(fp_sig)}

    # C2 non-noise classes (explicit exact enum set; no fallback/default-to-Noise).
    tp_nn = [r for r in tp_rows if r["_vc"] not in CURRENT_NOISE_CLASSES]
    fp_nn = [r for r in fp_rows if r["_vc"] not in CURRENT_NOISE_CLASSES]
    result["vc_not_noise"] = {"tp": len(tp_nn), "fp": len(fp_nn)}

    return result


def vc_auroc_features(tp_rows, fp_rows):
    """AUROC of key features for TP/FP discrimination."""
    metrics = {}
    features = {
        "GlobalP_inv":        (lambda r: -r["_global_p"],   True),
        "CramersV":           (lambda r: r["_cramers_v"],    True),
        "HeuristicScore":     (lambda r: r["_heuristic"],    True),
        "AlleleDelta":        (lambda r: r["_allele_delta"], False),  # FP is higher
        "PairwiseMeanDist":   (lambda r: r["_pairwise_mean"], False),
        "LOH_HP_Signal":      (lambda r: r["_LOH_HP_Signal"], False),  # FP is higher
        "HP_Coverage_Symmetry": (lambda r: r["_HP_Coverage_Symmetry"], True),
        "ClusterPermanovaF":  (lambda r: r["_perm_f"] if r["_perm_valid"] else 0.0, True),
    }
    for name, (fn, tp_higher) in features.items():
        tp_s = [fn(r) for r in tp_rows]
        fp_s = [fn(r) for r in fp_rows]
        if tp_higher:
            auc = auroc(tp_s, fp_s)
        else:
            # FP is higher → flip: FP-discriminating = 1-AUROC when TP=1 direction
            auc = 1.0 - auroc(tp_s, fp_s)
        metrics[name] = auc
    return metrics


# ── QualityScore Audit ─────────────────────────────────────────────────────────

def quality_audit(tp_rows, fp_rows):
    """QualityScore distribution by LOH and VerificationClass."""
    def bucket(rows, condition, label):
        subset = [r for r in rows if condition(r)]
        if not subset:
            return {"n": 0, "mean": float("nan"), "p25": float("nan"),
                    "median": float("nan"), "p75": float("nan")}
        scores = sorted(r["_quality_score"] for r in subset)
        n = len(scores)
        return {
            "n": n,
            "mean": sum(scores) / n,
            "p25": scores[n // 4],
            "median": scores[n // 2],
            "p75": scores[3 * n // 4],
        }

    result = {}
    result["tp_loh_true"]  = bucket(tp_rows, lambda r: r["_loh"], "TP LOH=True")
    result["tp_loh_false"] = bucket(tp_rows, lambda r: not r["_loh"], "TP LOH=False")
    result["fp_loh_true"]  = bucket(fp_rows, lambda r: r["_loh"], "FP LOH=True")
    result["fp_loh_false"] = bucket(fp_rows, lambda r: not r["_loh"], "FP LOH=False")

    # TP in canonical C2 noise classes by quality tier
    tp_noise = [r for r in tp_rows if r["_vc"] in CURRENT_NOISE_CLASSES]
    for tier in ("High", "Medium", "Low", "VeryLow", ""):
        key = f"tp_noise_tier_{tier or 'unknown'}"
        subset = [r for r in tp_noise if r["_quality_tier"] == tier]
        result[key] = len(subset)

    # LOH rate by canonical C2 VerificationClass.
    result["loh_rate_by_verification_class"] = {}
    for vc in CURRENT_CLASSES_V2:
        tp_vc = [r for r in tp_rows if r["_vc"] == vc]
        fp_vc = [r for r in fp_rows if r["_vc"] == vc]
        tp_loh = sum(1 for r in tp_vc if r["_loh"])
        fp_loh = sum(1 for r in fp_vc if r["_loh"])
        result["loh_rate_by_verification_class"][vc] = {
            "tp": tp_loh / len(tp_vc) if tp_vc else float("nan"),
            "fp": fp_loh / len(fp_vc) if fp_vc else float("nan"),
        }

    return result


# ── PERMANOVA Penalty Audit ────────────────────────────────────────────────────

def permanova_penalty_audit(tp_rows, fp_rows):
    """How many TP (passed_gate) have ClusterPermanovaP > 0.05 → get x0.7 penalty."""
    tp_pg = [r for r in tp_rows if r["_passed_gate"]]
    fp_pg = [r for r in fp_rows if r["_passed_gate"]]

    tp_perm_valid = [r for r in tp_pg if r["_perm_valid"]]
    fp_perm_valid = [r for r in fp_pg if r["_perm_valid"]]

    tp_penalized = [r for r in tp_perm_valid if r["_perm_p"] > 0.05]
    fp_penalized = [r for r in fp_perm_valid if r["_perm_p"] > 0.05]

    # Compare heuristic score distribution for penalized vs not
    def hs_stats(rows):
        if not rows:
            return {"n": 0, "mean": float("nan"), "median": float("nan")}
        scores = sorted(r["_heuristic"] for r in rows)
        n = len(scores)
        return {"n": n, "mean": sum(scores) / n, "median": scores[n // 2]}

    return {
        "tp_passed_gate": len(tp_pg),
        "fp_passed_gate": len(fp_pg),
        "tp_perm_valid": len(tp_perm_valid),
        "fp_perm_valid": len(fp_perm_valid),
        "tp_penalized_n": len(tp_penalized),
        "fp_penalized_n": len(fp_penalized),
        "tp_penalized_rate": len(tp_penalized) / len(tp_perm_valid) if tp_perm_valid else float("nan"),
        "fp_penalized_rate": len(fp_penalized) / len(fp_perm_valid) if fp_perm_valid else float("nan"),
        "tp_penalized_hs": hs_stats(tp_penalized),
        "tp_not_penalized_hs": hs_stats([r for r in tp_perm_valid if r["_perm_p"] <= 0.05]),
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def load_sample(sample, tp_path, fp_path):
    tp_raw = load_csv(tp_path)
    fp_raw = load_csv(fp_path)
    if not tp_raw:
        print(f"  [WARN] {sample}: no TP data at {tp_path}")
        return None, None
    tp_rows = tp_raw
    fp_rows = fp_raw
    print(f"  {sample}: TP={len(tp_rows)}, FP={len(fp_rows)}")
    return tp_rows, fp_rows


def collect_samples(include_extra=True):
    samples = {}
    print("\n=== Loading sample data ===")
    for sample, run_id in SAMPLE_PATHS.items():
        tp_path = f"{BASE}/{sample}/paired_full/{run_id}/intersubmod_tp/significance_summary.csv"
        fp_path = f"{BASE}/{sample}/paired_full/{run_id}/intersubmod_fp/significance_summary.csv"
        tp, fp = load_sample(sample, tp_path, fp_path)
        if tp is not None:
            samples[sample] = {"tp": tp, "fp": fp}

    if include_extra:
        for name, paths in EXTRA_PATHS.items():
            tp, fp = load_sample(name, paths["tp"], paths["fp"])
            if tp is not None:
                samples[name] = {"tp": tp, "fp": fp}

    return samples


def print_decision_tree(sample, flow, n_tp_total, n_fp_total):
    print(f"\n--- {sample} Decision Tree ---")
    total = flow["total"]
    print(f"  ALL: TP={total['tp']} FP={total['fp']}")
    pg = flow["passed_gate_true"]
    ng = flow["passed_gate_false"]
    print(f"  PassedGating=True:  TP={pg['tp']} ({fmt_pct(pg['tp'],total['tp'])}) "
          f"FP={pg['fp']} ({fmt_pct(pg['fp'],total['fp'])}) "
          f"Prec={precision(pg['tp'],pg['fp']):.1%}")
    print(f"  PassedGating=False: TP={ng['tp']} ({fmt_pct(ng['tp'],total['tp'])}) "
          f"FP={ng['fp']} ({fmt_pct(ng['fp'],total['fp'])})")

    print(f"  VerificationClass v2 breakdown:")
    for vc in CURRENT_CLASSES_V2:
        b = flow["verification_classes"][vc]
        prec = precision(b["tp"], b["fp"])
        print(f"    {vc:24s}: TP={b['tp']:5d} ({fmt_pct(b['tp'],total['tp'])}) "
              f"FP={b['fp']:5d} ({fmt_pct(b['fp'],total['fp'])}) "
              f"Prec={prec:.1%}")

    print("  Typed evidence support:")
    for key in ("cluster_first_support_true", "label_first_support_true", "bidirectional_support_true"):
        branch = flow[key]
        print(
            f"    {key:30s}: TP={branch['tp']:5d} FP={branch['fp']:5d} "
            f"Prec={precision(branch['tp'], branch['fp']):.1%}"
        )

    sig = flow["significant_true"]
    nn  = flow["vc_not_noise"]
    print(f"  Significant=True:   TP={sig['tp']:5d} ({fmt_pct(sig['tp'],total['tp'])}) "
          f"FP={sig['fp']:5d} ({fmt_pct(sig['fp'],total['fp'])}) "
          f"Prec={precision(sig['tp'],sig['fp']):.1%}")
    print(f"  C2 non-noise:       TP={nn['tp']:5d} ({fmt_pct(nn['tp'],total['tp'])}) "
          f"FP={nn['fp']:5d} ({fmt_pct(nn['fp'],total['fp'])}) "
          f"Prec={precision(nn['tp'],nn['fp']):.1%}")


def print_auroc_table(sample, aucs):
    print(f"\n--- {sample} AUROC (TP-discriminating direction) ---")
    for feat, auc in sorted(aucs.items(), key=lambda x: -x[1]):
        bar = "█" * int(auc * 20)
        print(f"  {feat:30s}: {auc:.3f} {bar}")


def print_quality_audit(sample, qa):
    print(f"\n--- {sample} QualityScore Audit ---")
    for key in ("tp_loh_true", "tp_loh_false", "fp_loh_true", "fp_loh_false"):
        d = qa[key]
        if d["n"] == 0:
            print(f"  {key}: n=0")
        else:
            print(f"  {key}: n={d['n']} mean={d['mean']:.1f} median={d['median']:.0f} "
                  f"p25={d['p25']:.0f} p75={d['p75']:.0f}")

    print(f"  TP in C2 noise classes by Quality_Tier:")
    for tier in ("High", "Medium", "Low", "VeryLow", "unknown"):
        key = f"tp_noise_tier_{tier}"
        print(f"    {tier}: {qa.get(key, 0)}")

    print(f"  LOH rate per C2 class:")
    for vc in CURRENT_CLASSES_V2:
        rates = qa["loh_rate_by_verification_class"][vc]
        tr, fr = rates["tp"], rates["fp"]
        print(
            f"    {vc:24s}: TP_LOH%={tr:.1%}  FP_LOH%={fr:.1%}"
            if not (math.isnan(tr) or math.isnan(fr))
            else f"    {vc}: N/A"
        )


def print_permanova_audit(sample, pa):
    print(f"\n--- {sample} PERMANOVA Penalty (x0.7) Audit ---")
    print(f"  TP passed_gate={pa['tp_passed_gate']}, perm_valid={pa['tp_perm_valid']}, "
          f"penalized={pa['tp_penalized_n']} ({pa['tp_penalized_rate']:.1%})")
    print(f"  FP passed_gate={pa['fp_passed_gate']}, perm_valid={pa['fp_perm_valid']}, "
          f"penalized={pa['fp_penalized_n']} ({pa['fp_penalized_rate']:.1%})")
    phs = pa["tp_penalized_hs"]
    nhs = pa["tp_not_penalized_hs"]
    if phs["n"] > 0 and nhs["n"] > 0:
        print(f"  TP penalized HeuristicScore: n={phs['n']} mean={phs['mean']:.2f} median={phs['median']:.2f}")
        print(f"  TP NOT penalized HeuristicScore: n={nhs['n']} mean={nhs['mean']:.2f} median={nhs['median']:.2f}")


def build_markdown(all_results):
    """Generate the methodology .md document."""
    lines = []
    lines += [
        "<!--",
        "建立時間: 2026-03-23",
        "目標: VerificationClass 決策樹跨樣本量化分析（Steps 2/4/5/6/7/8）",
        "處理範圍: 7 個樣本 + HCC1395 pileup n=99",
        "關聯檔案:",
        "  - src/core/SignificanceAnalyzer.cpp (304-337)",
        "  - src/core/GlobalTest.cpp (137-142)",
        "  - scripts/pipeline/steps/03_filter_analysis.py",
        "  - scripts/analysis/verify_class_decision_tree_audit.py",
        "-->",
        "",
        "# VerificationClass 決策樹跨樣本量化分析",
        "",
        "## 1. 定義",
        "",
        "### VerificationClass v2 與 evidence 契約",
        "```",
        "current_class = VerificationClass (VerificationSchemaVersion == 2)",
        "label_support = typed LabelFirstSupport",
        "cluster_support = typed ClusterFirstSupport",
        "non_noise = current_class not in {Noise_Uniform, Noise_Chaotic, Noise_Uncorrelated}",
        "# evidence is never reconstructed from a class name",
        "```",
        "",
        "### 欄位來源（significance_summary.csv）",
        "| 欄位 | 說明 |",
        "|------|------|",
        "| PassedGating | passed_gate 結果 |",
        "| CramersV | Cramér's V 效應量 |",
        "| GlobalP | 全局 p-value（AlleleDelta chi-sq） |",
        "| VerificationSchemaVersion | 必須逐列為 2 |",
        "| VerificationClass | canonical 11-class C2 taxonomy |",
        "| LabelFirstSupport | typed label-first evidence |",
        "| ClusterFirstSupport | typed cluster-first evidence |",
        "| EvidencePath | canonical evidence path enum |",
        "| Significant | Python 層判斷（候選去除） |",
        "| Potential_LOH | LOH 旗標 |",
        "| Quality_Score | QualityScore (0-100) |",
        "| ClusterPermanovaF | PERMANOVA F 統計量 |",
        "| ClusterPermanovaP | PERMANOVA p-value |",
        "",
    ]

    # Decision tree Mermaid diagrams
    lines += [
        "## 2. 決策樹流量（跨樣本）",
        "",
    ]

    for sample, data in all_results.items():
        flow = data["flow"]
        total = flow["total"]
        if total["fp"] == 0:
            continue  # Skip samples with no FP (can't show precision)

        pg  = flow["passed_gate_true"]
        ng  = flow["passed_gate_false"]
        cfs = flow["cluster_first_support_true"]
        lfs = flow["label_first_support_true"]
        bidir = flow["bidirectional_support_true"]
        sig = flow["significant_true"]
        nn  = flow["vc_not_noise"]

        lines += [
            f"### {sample}",
            f"",
            f"```mermaid",
            f"flowchart TD",
            f'  ALL["全量<br>TP={total["tp"]:,} FP={total["fp"]:,}"]',
            f'  ALL --> CFS["ClusterFirstSupport=True<br>TP={cfs["tp"]:,} FP={cfs["fp"]:,}<br>Prec={precision(cfs["tp"],cfs["fp"]):.0%}"]',
            f'  ALL --> LFS["LabelFirstSupport=True<br>TP={lfs["tp"]:,} FP={lfs["fp"]:,}<br>Prec={precision(lfs["tp"],lfs["fp"]):.0%}"]',
            f'  ALL --> BI["Both evidence arms=True<br>TP={bidir["tp"]:,} FP={bidir["fp"]:,}<br>Prec={precision(bidir["tp"],bidir["fp"]):.0%}"]',
            f"```",
            f"",
            f"| 判斷標準 | TP | FP | TP% | FP% | Precision |",
            f"|---------|----|----|-----|-----|-----------|",
            f"| 全量 | {total['tp']:,} | {total['fp']:,} | 100% | 100% | {precision(total['tp'],total['fp']):.1%} |",
            f"| PassedGating=True | {pg['tp']:,} | {pg['fp']:,} | {fmt_pct(pg['tp'],total['tp'])} | {fmt_pct(pg['fp'],total['fp'])} | {precision(pg['tp'],pg['fp']):.1%} |",
            f"| ClusterFirstSupport=True | {cfs['tp']:,} | {cfs['fp']:,} | {fmt_pct(cfs['tp'],total['tp'])} | {fmt_pct(cfs['fp'],total['fp'])} | {precision(cfs['tp'],cfs['fp']):.1%} |",
            f"| LabelFirstSupport=True | {lfs['tp']:,} | {lfs['fp']:,} | {fmt_pct(lfs['tp'],total['tp'])} | {fmt_pct(lfs['fp'],total['fp'])} | {precision(lfs['tp'],lfs['fp']):.1%} |",
            f"| Both evidence arms=True | {bidir['tp']:,} | {bidir['fp']:,} | {fmt_pct(bidir['tp'],total['tp'])} | {fmt_pct(bidir['fp'],total['fp'])} | {precision(bidir['tp'],bidir['fp']):.1%} |",
            f"| Significant=True | {sig['tp']:,} | {sig['fp']:,} | {fmt_pct(sig['tp'],total['tp'])} | {fmt_pct(sig['fp'],total['fp'])} | {precision(sig['tp'],sig['fp']):.1%} |",
            f"| C2 non-noise | {nn['tp']:,} | {nn['fp']:,} | {fmt_pct(nn['tp'],total['tp'])} | {fmt_pct(nn['fp'],total['fp'])} | {precision(nn['tp'],nn['fp']):.1%} |",
        ]
        for vc in CURRENT_CLASSES_V2:
            branch = flow["verification_classes"][vc]
            lines.append(
                f"| C2 `{vc}` | {branch['tp']:,} | {branch['fp']:,} | "
                f"{fmt_pct(branch['tp'],total['tp'])} | {fmt_pct(branch['fp'],total['fp'])} | "
                f"{precision(branch['tp'],branch['fp']):.1%} |"
            )
        lines.append("")

    # AUROC summary table
    lines += [
        "## 3. 特徵辨別力 AUROC（跨樣本）",
        "",
        "AUROC 為 TP-discriminating 方向（>0.5 表示 TP 此特徵值更高，<0.5 表示 FP 更高）",
        "",
    ]

    feature_names = ["GlobalP_inv", "CramersV", "HeuristicScore", "AlleleDelta",
                     "PairwiseMeanDist", "LOH_HP_Signal", "HP_Coverage_Symmetry", "ClusterPermanovaF"]
    header = "| 特徵 | " + " | ".join(s for s in all_results.keys() if all_results[s]["flow"]["total"]["fp"] > 5) + " |"
    sep = "|------" + "|------" * sum(1 for s in all_results if all_results[s]["flow"]["total"]["fp"] > 5) + "|"
    lines += [header, sep]
    for feat in feature_names:
        vals = []
        for sample, data in all_results.items():
            if data["flow"]["total"]["fp"] <= 5:
                continue
            auc = data["auroc"].get(feat, float("nan"))
            if math.isnan(auc):
                vals.append("N/A")
            else:
                vals.append(f"{auc:.3f}")
        lines.append(f"| {feat} | " + " | ".join(vals) + " |")
    lines.append("")

    # Step 4: Significant vs canonical C2 non-noise
    lines += [
        "## 4. Significant 欄位評估",
        "",
        "比較 Significant=True 與 C2 non-noise 作為分類準則的差異：",
        "",
        "| 樣本 | Sig=True TP% | Sig=True FP% | Sig=True Prec | C2 non-noise TP% | C2 non-noise FP% | C2 non-noise Prec |",
        "|------|-------------|-------------|--------------|--------------|--------------|---------------|",
    ]
    for sample, data in all_results.items():
        flow = data["flow"]
        total = flow["total"]
        if total["fp"] == 0:
            continue
        sig = flow["significant_true"]
        nn  = flow["vc_not_noise"]
        lines.append(
            f"| {sample} | {fmt_pct(sig['tp'],total['tp'])} | {fmt_pct(sig['fp'],total['fp'])} | "
            f"{precision(sig['tp'],sig['fp']):.1%} | "
            f"{fmt_pct(nn['tp'],total['tp'])} | {fmt_pct(nn['fp'],total['fp'])} | "
            f"{precision(nn['tp'],nn['fp']):.1%} |"
        )
    lines.append("")

    # Step 5: QualityScore audit
    lines += [
        "## 5. QualityScore 審查",
        "",
        "### LOH 懲罰影響（-25分）",
        "",
        "| 樣本 | TP LOH=T mean | TP LOH=F mean | FP LOH=T mean | FP LOH=F mean |",
        "|------|--------------|--------------|--------------|--------------|",
    ]
    for sample, data in all_results.items():
        qa = data["quality"]
        t_lt = qa["tp_loh_true"]
        t_lf = qa["tp_loh_false"]
        f_lt = qa["fp_loh_true"]
        f_lf = qa["fp_loh_false"]
        lines.append(
            f"| {sample} | {t_lt['mean']:.0f}(n={t_lt['n']}) | "
            f"{t_lf['mean']:.0f}(n={t_lf['n']}) | "
            f"{f_lt['mean']:.0f}(n={f_lt['n']}) | "
            f"{f_lf['mean']:.0f}(n={f_lf['n']}) |"
        )
    lines.append("")

    # Step 6: LOH rate by selected C2 classes
    lines += [
        "### LOH 比率（selected C2 classes）",
        "",
        "| 樣本 | Strong_Bidirectional TP/FP LOH% | ClusterFirstOnly TP/FP LOH% | C2 noise TP LOH% |",
        "|------|------------------|---------------------|--------------|",
    ]
    for sample, data in all_results.items():
        qa = data["quality"]
        total = data["flow"]["total"]
        if total["fp"] == 0:
            continue
        rates = qa["loh_rate_by_verification_class"]
        strong = rates["Strong_Bidirectional"]
        cluster = rates["ClusterFirstOnly"]
        tp_rows = data["flow"]["total"]["tp"]
        noise_tp_count = sum(
            data["flow"]["verification_classes"][vc]["tp"] for vc in CURRENT_NOISE_CLASSES
        )
        lines.append(
            f"| {sample} | "
            f"TP={strong['tp']:.0%} / FP={strong['fp']:.0%} | "
            f"TP={cluster['tp']:.0%} / FP={cluster['fp']:.0%} | "
            f"{fmt_pct(noise_tp_count, tp_rows)} |"
        )
    lines.append("")

    # Step 8: PERMANOVA penalty
    lines += [
        "## 8. PERMANOVA 懲罰（×0.7 if ClusterPermanovaP>0.05）",
        "",
        "| 樣本 | TP passed_gate | TP perm_valid | TP penalized% | FP penalized% | TP penalized mean_HS | TP normal mean_HS |",
        "|------|--------------|--------------|--------------|--------------|---------------------|-----------------|",
    ]
    for sample, data in all_results.items():
        pa = data["permanova"]
        total = data["flow"]["total"]
        if total["fp"] == 0:
            continue
        p_hs = pa["tp_penalized_hs"]
        n_hs = pa["tp_not_penalized_hs"]
        lines.append(
            f"| {sample} | {pa['tp_passed_gate']} | {pa['tp_perm_valid']} | "
            f"{pa['tp_penalized_rate']:.0%} | {pa['fp_penalized_rate']:.0%} | "
            f"{p_hs['mean']:.2f} | {n_hs['mean']:.2f} |"
        )

    lines += [
        "",
        "## 結論",
        "",
        "### Step 2: VerificationClass 決策樹",
        "- **Pending**: 需結合 3+ 樣本數據確認 passed_gate 雙重閾值合理性",
        "",
        "### Step 4: Significant 欄位",
        "- **Pending**: 對比 Significant=True vs C2 non-noise 的 Precision 差異確認是否需要移除",
        "",
        "### Step 5: QualityScore LOH 懲罰",
        "- **Pending**: 確認 LOH TP 的 QualityScore 是否顯著低於非 LOH TP",
        "",
        "### Step 6: LOH 在各 VC class 的比例",
        "- **Pending**: 確認 FP 是否有更高的 LOH 比率（若有，LOH 可作為 FP 指標）",
        "",
        "### Step 8: PERMANOVA 懲罰",
        "- **Pending**: 確認懲罰比率與 HeuristicScore 影響程度",
        "",
    ]

    return "\n".join(lines)


def main():
    global BASE, OUTPUT_DOC
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", default=BASE, help="Canonical sample root")
    parser.add_argument("--output-doc", default=OUTPUT_DOC)
    parser.add_argument(
        "--skip-extra",
        action="store_true",
        help="Explicitly exclude the optional /tmp n=99 reference inputs.",
    )
    args = parser.parse_args()
    BASE = args.base
    OUTPUT_DOC = args.output_doc
    try:
        samples = collect_samples(include_extra=not args.skip_extra)
    except SchemaContractError as exc:
        raise SystemExit(f"[decision-tree][schema-contract] {exc}") from exc
    if not samples:
        print("ERROR: No sample data loaded")
        sys.exit(1)

    all_results = {}

    print("\n" + "=" * 70)
    print("STEP 2: VerificationClass Decision Tree Quantitative Analysis")
    print("=" * 70)

    for sample, data in samples.items():
        tp_rows = data["tp"]
        fp_rows = data["fp"]
        flow = decision_tree_flow(tp_rows, fp_rows)
        aucs = vc_auroc_features(tp_rows, fp_rows)
        qa   = quality_audit(tp_rows, fp_rows)
        pa   = permanova_penalty_audit(tp_rows, fp_rows)

        all_results[sample] = {
            "flow": flow,
            "auroc": aucs,
            "quality": qa,
            "permanova": pa,
        }

        print_decision_tree(sample, flow, len(tp_rows), len(fp_rows))
        if flow["total"]["fp"] > 5:
            print_auroc_table(sample, aucs)
        print_quality_audit(sample, qa)
        print_permanova_audit(sample, pa)

    # Cross-sample summary
    print("\n" + "=" * 70)
    print("CROSS-SAMPLE SUMMARY: VerificationClass Precision")
    print("=" * 70)
    print(
        f"{'Sample':25s} {'StrongBi Prec':13s} {'ClusterOnly':12s} "
        f"{'CFS Prec':10s} {'Noise TP%':10s} {'SigTrue Prec':12s} {'C2!=N Prec':10s}"
    )
    for sample, data in all_results.items():
        flow = data["flow"]
        total = flow["total"]
        if total["fp"] == 0:
            print(f"  {sample}: FP=0, skipping precision")
            continue
        s = flow["verification_classes"]["Strong_Bidirectional"]
        sc = flow["verification_classes"]["ClusterFirstOnly"]
        cfs = flow["cluster_first_support_true"]
        sig = flow["significant_true"]
        nn  = flow["vc_not_noise"]
        noise_tp = sum(flow["verification_classes"][vc]["tp"] for vc in CURRENT_NOISE_CLASSES)
        noise_tp_rate = noise_tp / total["tp"] if total["tp"] > 0 else 0
        print(f"  {sample:25s} "
              f"{precision(s['tp'],s['fp']):11.1%}  "
              f"{precision(sc['tp'],sc['fp']):10.1%}  "
              f"{precision(cfs['tp'],cfs['fp']):8.1%}  "
              f"{noise_tp_rate:8.1%}  "
              f"{precision(sig['tp'],sig['fp']):10.1%}  "
              f"{precision(nn['tp'],nn['fp']):8.1%}")

    # Write markdown
    os.makedirs(os.path.dirname(OUTPUT_DOC), exist_ok=True)
    md_content = build_markdown(all_results)
    with open(OUTPUT_DOC, "w") as f:
        f.write(md_content)
    print(f"\n[OK] Methodology doc written to: {OUTPUT_DOC}")


if __name__ == "__main__":
    main()
