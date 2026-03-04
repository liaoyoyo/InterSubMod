#!/usr/bin/env python3
"""
evaluate_dorado_purity_filter_components.py

目的：
1) 在 HCC1395_DORADO purity 結果上，分解評估：
   - QUAL 條件：QUAL < q
   - 甲基條件：AD > a AND CV < c AND VAF < v
   - OR 組合：QUAL 或 甲基條件
2) 驗證目前固定門檻與舊門檻的實際影響（TP/FP/FN/P/R/F1）
3) 進行門檻搜尋，找出全域/分 purity 最佳方案

輸入：
- /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395_DORADO/purity_t*_20260213_dorado_purity_full/

輸出：
- components_eval.tsv
- grid_eval_top20.tsv
- best_global.json
"""

from __future__ import annotations

import csv
import gzip
import itertools
import json
import math
from pathlib import Path
from typing import Callable, Dict, List, Tuple


TRUTH_TOTAL = 39447
BASE = Path(
    "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395_DORADO"
)
RUN_TAG = "20260213_dorado_purity_full"

PURITY_DIRS = sorted(BASE.glob(f"purity_t*_{RUN_TAG}"))
OUT_DIR = BASE / "purity_runs" / RUN_TAG / "component_analysis"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_vcf_features(vcf_path: Path) -> Dict[Tuple[str, int], Dict[str, float]]:
    """讀取 QUAL + VAF (chr,pos)"""
    features: Dict[Tuple[str, int], Dict[str, float]] = {}
    opener = gzip.open if vcf_path.suffix == ".gz" else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            chrom = p[0]
            pos = int(p[1])
            qual = 0.0 if p[5] == "." else float(p[5])

            fmt = p[8].split(":")
            sample = p[9].split(":")
            d = dict(zip(fmt, sample))
            vaf = 0.0
            if "VAF" in d:
                try:
                    vaf = float(d["VAF"])
                except ValueError:
                    vaf = 0.0
            elif "AF" in d:
                try:
                    vaf = float(d["AF"])
                except ValueError:
                    vaf = 0.0
            elif "AD" in d:
                try:
                    ad = d["AD"].split(",")
                    if len(ad) >= 2:
                        r = int(ad[0])
                        a = int(ad[1])
                        t = r + a
                        vaf = a / t if t > 0 else 0.0
                except Exception:
                    vaf = 0.0
            features[(chrom, pos)] = {"qual": qual, "vaf": vaf}
    return features


def read_summary(csv_path: Path) -> List[dict]:
    rows = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
    return rows


def metric(tp: int, fp: int, truth_total: int = TRUTH_TOTAL) -> dict:
    fn = truth_total - tp
    p = tp / (tp + fp) if tp + fp > 0 else 0.0
    r = tp / truth_total if truth_total > 0 else 0.0
    f1 = 2 * p * r / (p + r) if p + r > 0 else 0.0
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": p,
        "recall": r,
        "f1": f1,
    }


def row_values(row: dict, feats: Dict[Tuple[str, int], Dict[str, float]]) -> dict:
    chrom = row["Chr"]
    pos = int(row["Pos"])
    ad = abs(float(row.get("AlleleDelta", 0.0)))
    cv = float(row.get("CramersV", 0.0))
    q = feats.get((chrom, pos), {}).get("qual")
    vaf = feats.get((chrom, pos), {}).get("vaf")
    return {"ad": ad, "cv": cv, "qual": q, "vaf": vaf}


def eval_filter(
    tp_rows: List[dict],
    fp_rows: List[dict],
    tp_feat: Dict[Tuple[str, int], Dict[str, float]],
    fp_feat: Dict[Tuple[str, int], Dict[str, float]],
    baseline_tp: int,
    baseline_fp: int,
    rule: Callable[[dict], bool],
) -> dict:
    tp_removed = 0
    fp_removed = 0
    for r in tp_rows:
        if rule(row_values(r, tp_feat)):
            tp_removed += 1
    for r in fp_rows:
        if rule(row_values(r, fp_feat)):
            fp_removed += 1
    m = metric(baseline_tp - tp_removed, baseline_fp - fp_removed)
    m["tp_removed"] = tp_removed
    m["fp_removed"] = fp_removed
    return m


def get_purity_name(path: Path) -> str:
    # purity_t19_n29_20260213... -> t19_n29
    name = path.name
    return name.split("_", 2)[1] + "_" + name.split("_", 3)[2]


def main() -> None:
    # 收集資料
    samples = []
    for d in PURITY_DIRS:
        mpath = d / "metrics.json"
        tp_csv = d / "intersubmod_tp" / "significance_summary.csv"
        fp_csv = d / "intersubmod_fp" / "significance_summary.csv"
        tp_vcf = d / "longphase_s" / "filtered_snv_tp.vcf.gz"
        fp_vcf = d / "longphase_s" / "filtered_snv_fp.vcf.gz"
        if not all(p.exists() for p in [mpath, tp_csv, fp_csv, tp_vcf, fp_vcf]):
            continue
        metrics = json.load(open(mpath))
        samples.append(
            {
                "name": get_purity_name(d),
                "dir": d,
                "baseline": metrics["baseline"],
                "tp_rows": read_summary(tp_csv),
                "fp_rows": read_summary(fp_csv),
                "tp_feat": parse_vcf_features(tp_vcf),
                "fp_feat": parse_vcf_features(fp_vcf),
            }
        )

    # 1) 組件分解：你指定的條件
    # A: QUAL < 0.75
    # B: AD > 0.25 AND CV < 0.05 AND VAF < 0.24
    # C: A OR B
    def rule_qual(x):
        return x["qual"] is not None and x["qual"] < 0.75

    def rule_acv(x):
        return (x["ad"] > 0.25) and (x["cv"] < 0.05) and (x["vaf"] is not None and x["vaf"] < 0.24)

    def rule_or(x):
        return rule_qual(x) or rule_acv(x)

    # D: 目前程式預設（QUAL停用）
    # AD > 0.15 AND CV < 0.03 AND VAF < 0.15
    def rule_current(x):
        return (x["ad"] > 0.15) and (x["cv"] < 0.03) and (x["vaf"] is not None and x["vaf"] < 0.15)

    component_rows: List[dict] = []
    for s in samples:
        base = s["baseline"]
        base_f1 = float(base["f1"])
        for tag, rule in [
            ("QUAL_only_q0.75", rule_qual),
            ("ACV_only_a0.25_c0.05_v0.24", rule_acv),
            ("QUAL_OR_ACV", rule_or),
            ("CURRENT_a0.15_c0.03_v0.15", rule_current),
        ]:
            res = eval_filter(
                s["tp_rows"],
                s["fp_rows"],
                s["tp_feat"],
                s["fp_feat"],
                int(base["tp"]),
                int(base["fp"]),
                rule,
            )
            component_rows.append(
                {
                    "purity": s["name"],
                    "rule": tag,
                    "baseline_f1": round(base_f1, 4),
                    "final_f1": round(res["f1"], 4),
                    "f1_delta": round(res["f1"] - base_f1, 4),
                    "tp_removed": res["tp_removed"],
                    "fp_removed": res["fp_removed"],
                    "final_tp": res["tp"],
                    "final_fp": res["fp"],
                    "final_fn": res["fn"],
                    "final_precision": round(res["precision"], 4),
                    "final_recall": round(res["recall"], 4),
                }
            )

    comp_tsv = OUT_DIR / "components_eval.tsv"
    with open(comp_tsv, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "purity",
                "rule",
                "baseline_f1",
                "final_f1",
                "f1_delta",
                "tp_removed",
                "fp_removed",
                "final_tp",
                "final_fp",
                "final_fn",
                "final_precision",
                "final_recall",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(component_rows)

    # 2) 門檻搜尋（OR 組合）：QUAL<q OR (AD>a & CV<c & VAF<v)
    q_vals = [0.60, 0.65, 0.70, 0.75, 0.80]
    a_vals = [0.15, 0.20, 0.25, 0.30, 0.35]
    c_vals = [0.02, 0.03, 0.05, 0.07]
    v_vals = [0.10, 0.15, 0.20, 0.24, 0.30]

    grid_rows = []
    for q, a, c, v in itertools.product(q_vals, a_vals, c_vals, v_vals):
        deltas = []
        per = {}
        for s in samples:
            base = s["baseline"]
            base_f1 = float(base["f1"])

            def rule(x, q=q, a=a, c=c, v=v):
                cond_q = (x["qual"] is not None and x["qual"] < q)
                cond_m = (x["ad"] > a) and (x["cv"] < c) and (x["vaf"] is not None and x["vaf"] < v)
                return cond_q or cond_m

            res = eval_filter(
                s["tp_rows"],
                s["fp_rows"],
                s["tp_feat"],
                s["fp_feat"],
                int(base["tp"]),
                int(base["fp"]),
                rule,
            )
            delta = res["f1"] - base_f1
            deltas.append(delta)
            per[s["name"]] = {
                "delta": round(delta, 4),
                "f1": round(res["f1"], 4),
                "tp_removed": res["tp_removed"],
                "fp_removed": res["fp_removed"],
            }

        avg = sum(deltas) / len(deltas)
        min_delta = min(deltas)
        grid_rows.append(
            {
                "q": q,
                "a": a,
                "c": c,
                "v": v,
                "avg_delta": avg,
                "min_delta": min_delta,
                "all_non_negative": int(min_delta >= 0.0),
                "detail": json.dumps(per, ensure_ascii=False),
            }
        )

    # 排序：優先 all_non_negative, 再 avg_delta
    grid_rows.sort(key=lambda x: (x["all_non_negative"], x["avg_delta"]), reverse=True)

    top20_tsv = OUT_DIR / "grid_eval_top20.tsv"
    with open(top20_tsv, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["q", "a", "c", "v", "avg_delta", "min_delta", "all_non_negative", "detail"],
            delimiter="\t",
        )
        w.writeheader()
        for r in grid_rows[:20]:
            rr = r.copy()
            rr["avg_delta"] = round(rr["avg_delta"], 6)
            rr["min_delta"] = round(rr["min_delta"], 6)
            w.writerow(rr)

    best = grid_rows[0]
    with open(OUT_DIR / "best_global.json", "w") as f:
        json.dump(
            {
                "best": {
                    "q": best["q"],
                    "a": best["a"],
                    "c": best["c"],
                    "v": best["v"],
                    "avg_delta": round(best["avg_delta"], 6),
                    "min_delta": round(best["min_delta"], 6),
                    "all_non_negative": bool(best["all_non_negative"]),
                    "detail": json.loads(best["detail"]),
                },
                "note": "目標函式：先最大化 all_non_negative，再最大化 avg_delta",
            },
            f,
            indent=2,
            ensure_ascii=False,
        )

    # 3) ACV-only 更廣門檻搜尋（不含 QUAL）
    a_vals2 = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
    c_vals2 = [0.005, 0.01, 0.02, 0.03, 0.05, 0.07]
    v_vals2 = [0.05, 0.10, 0.15, 0.20, 0.24, 0.30]
    acv_rows = []
    per_purity_best = {}
    for a, c, v in itertools.product(a_vals2, c_vals2, v_vals2):
        deltas = []
        per = {}
        for s in samples:
            base = s["baseline"]
            base_f1 = float(base["f1"])

            def rule(x, a=a, c=c, v=v):
                return (x["ad"] > a) and (x["cv"] < c) and (x["vaf"] is not None and x["vaf"] < v)

            res = eval_filter(
                s["tp_rows"],
                s["fp_rows"],
                s["tp_feat"],
                s["fp_feat"],
                int(base["tp"]),
                int(base["fp"]),
                rule,
            )
            delta = res["f1"] - base_f1
            deltas.append(delta)
            per[s["name"]] = {
                "delta": round(delta, 4),
                "f1": round(res["f1"], 4),
                "tp_removed": res["tp_removed"],
                "fp_removed": res["fp_removed"],
            }
            # 每個 purity 的最佳參數（最大 delta）
            prev = per_purity_best.get(s["name"])
            cand = {"a": a, "c": c, "v": v, "delta": delta, "f1": res["f1"], "tp_removed": res["tp_removed"], "fp_removed": res["fp_removed"]}
            if prev is None or cand["delta"] > prev["delta"]:
                per_purity_best[s["name"]] = cand

        avg = sum(deltas) / len(deltas)
        min_delta = min(deltas)
        acv_rows.append(
            {
                "a": a,
                "c": c,
                "v": v,
                "avg_delta": avg,
                "min_delta": min_delta,
                "all_non_negative": int(min_delta >= 0.0),
                "detail": json.dumps(per, ensure_ascii=False),
            }
        )

    acv_rows.sort(key=lambda x: (x["all_non_negative"], x["avg_delta"]), reverse=True)
    acv_top20_tsv = OUT_DIR / "acv_only_grid_top20.tsv"
    with open(acv_top20_tsv, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["a", "c", "v", "avg_delta", "min_delta", "all_non_negative", "detail"],
            delimiter="\t",
        )
        w.writeheader()
        for r in acv_rows[:20]:
            rr = r.copy()
            rr["avg_delta"] = round(rr["avg_delta"], 6)
            rr["min_delta"] = round(rr["min_delta"], 6)
            w.writerow(rr)

    with open(OUT_DIR / "acv_only_best_global.json", "w") as f:
        b = acv_rows[0]
        json.dump(
            {
                "best": {
                    "a": b["a"],
                    "c": b["c"],
                    "v": b["v"],
                    "avg_delta": round(b["avg_delta"], 6),
                    "min_delta": round(b["min_delta"], 6),
                    "all_non_negative": bool(b["all_non_negative"]),
                    "detail": json.loads(b["detail"]),
                },
                "per_purity_best": {
                    k: {
                        "a": round(vv["a"], 4),
                        "c": round(vv["c"], 4),
                        "v": round(vv["v"], 4),
                        "delta": round(vv["delta"], 6),
                        "f1": round(vv["f1"], 6),
                        "tp_removed": vv["tp_removed"],
                        "fp_removed": vv["fp_removed"],
                    }
                    for k, vv in sorted(per_purity_best.items())
                },
                "note": "ACV-only: AD>a AND CV<c AND VAF<v",
            },
            f,
            indent=2,
            ensure_ascii=False,
        )

    print(f"[OK] components: {comp_tsv}")
    print(f"[OK] grid top20: {top20_tsv}")
    print(f"[OK] best: {OUT_DIR / 'best_global.json'}")
    print(f"[OK] acv grid top20: {acv_top20_tsv}")
    print(f"[OK] acv best: {OUT_DIR / 'acv_only_best_global.json'}")


if __name__ == "__main__":
    main()
