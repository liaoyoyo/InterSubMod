#!/usr/bin/env python3
"""lineage 管線回歸驗證器 — 改動後跑這支，自動判斷結果是否合理。

三層檢查：
  L1 守恆    : 結果自身必須滿足的算術恆等式（違反 = 一定有 bug）
  L2 baseline: 與凍結基準逐欄比對（超出容差 = 需要人工解釋）
  L3 已知問題: 追蹤 KNOWN_ISSUES 登記的異常是否仍存在 / 是否惡化

用法:
    # 驗證 baseline 自身（自檢，確認驗證器邏輯正確）
    python3 pipeline/lineage/verify.py --self-test

    # 驗證一份新結果
    python3 pipeline/lineage/verify.py --candidate <new_summary.json> [--tolerance 0.0]

離開碼: 0=PASS  1=L2/L3 有差異需確認  2=L1 守恆違反(必有 bug)  3=輸入錯誤
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_BASELINE = HERE / "baseline" / "cohort_baseline.json"

# ── L1 守恆規則：(名稱, 檢查函式) — 對每個 sample 執行 ──────────────────
def _ak_sum(s: dict) -> int:
    return sum(int(v) for v in (s.get("active_k_distribution") or {}).values())


CONSERVATION = [
    (
        "C1 active_k 總和 == groups_total",
        lambda s: _ak_sum(s) == s["groups_total"],
        lambda s: f"{_ak_sum(s)} vs {s['groups_total']}",
    ),
    (
        "C2 active_k[0] == no_active_alt_units",
        lambda s: int((s.get("active_k_distribution") or {}).get("0", 0)) == s["no_active_alt_units"],
        lambda s: f"{(s.get('active_k_distribution') or {}).get('0')} vs {s['no_active_alt_units']}",
    ),
    (
        "C3 groups_total - active_k[0] == mutation_bearing_units",
        lambda s: s["groups_total"] - int((s.get("active_k_distribution") or {}).get("0", 0))
        == s["mutation_bearing_units"],
        lambda s: f"{s['groups_total'] - int((s.get('active_k_distribution') or {}).get('0', 0))} vs {s['mutation_bearing_units']}",
    ),
    (
        "C4 objective_certified + objective_abstain == groups_total",
        lambda s: s["objective_certified_units"] + s["objective_abstain_units"] == s["groups_total"],
        lambda s: f"{s['objective_certified_units']} + {s['objective_abstain_units']} vs {s['groups_total']}",
    ),
    (
        "C5 mutation_family_complete + abstain == mutation_bearing_units",
        lambda s: s["mutation_family_complete_units"] + s["mutation_family_abstain_units"]
        == s["mutation_bearing_units"],
        lambda s: f"{s['mutation_family_complete_units']} + {s['mutation_family_abstain_units']} vs {s['mutation_bearing_units']}",
    ),
    (
        "C6 best_tree_unique + best_tree_tied == ranked_units",
        lambda s: s["best_tree_unique_units"] + s["best_tree_tied_units"] == s["ranked_units"],
        lambda s: f"{s['best_tree_unique_units']} + {s['best_tree_tied_units']} vs {s['ranked_units']}",
    ),
    (
        "C7 HP1_W + HP2_W == W_total",
        lambda s: (s.get("HP1_W") or 0) + (s.get("HP2_W") or 0) == s.get("W_total"),
        lambda s: f"{s.get('HP1_W')} + {s.get('HP2_W')} vs {s.get('W_total')}",
    ),
    (
        "C8 active_k 最大值 <= 12 (k>12 必已切分)",
        lambda s: all(int(k) <= 12 for k in (s.get("active_k_distribution") or {})),
        lambda s: f"max_k={max((int(k) for k in (s.get('active_k_distribution') or {})), default=0)}",
    ),
]

# ── L3 已知問題（與 docs/KNOWN_ISSUES.md 同步）────────────────────────
KNOWN_ISSUES = [
    {
        "id": "D1",
        "desc": "COLO829 的 groups_total 少於 W_total（其他樣本皆增加）",
        "sample": "COLO829",
        "probe": lambda s: s["groups_total"] - s["W_total"],
        "baseline_value": -14,
        "severity": "high",
        "note": "表示有 W 未進入 topology 階段，來源未查明",
    },
    {
        "id": "D5",
        "desc": "tie 比例偏高，影響 read 標籤的唯一性",
        "sample": None,  # cohort 層級
        "probe": None,
        "baseline_value": None,
        "severity": "high",
        "note": "cohort best_tree_tied / ranked",
    },
]


def load(path: Path) -> dict:
    if not path.is_file():
        print(f"NOT FOUND: {path}", file=sys.stderr)
        sys.exit(3)
    return json.loads(path.read_text(encoding="utf-8"))


def check_conservation(samples: dict) -> list[str]:
    failures = []
    for name, s in samples.items():
        for label, fn, detail in CONSERVATION:
            try:
                ok = fn(s)
            except (KeyError, TypeError) as exc:
                failures.append(f"  {name:16s} {label}  -> 欄位缺失 {exc}")
                continue
            if not ok:
                failures.append(f"  {name:16s} {label}  -> {detail(s)}")
    return failures


def compare_to_baseline(base: dict, cand: dict, tol: float) -> list[str]:
    diffs = []
    bs, cs = base["samples"], cand.get("samples", {})
    for name in base["datasets"]:
        if name not in cs:
            diffs.append(f"  {name:16s} 候選結果缺少此樣本")
            continue
        b, c = bs[name], cs[name]
        for field, bv in b.items():
            if isinstance(bv, dict) or bv is None:
                continue
            cv = c.get(field)
            if cv is None:
                diffs.append(f"  {name:16s} {field:36s} baseline={bv} candidate=缺失")
                continue
            if isinstance(bv, (int, float)) and isinstance(cv, (int, float)):
                if bv == cv:
                    continue
                if bv != 0 and abs(cv - bv) / abs(bv) <= tol:
                    continue
                pct = "n/a" if bv == 0 else f"{(cv - bv) / bv * 100:+.2f}%"
                diffs.append(f"  {name:16s} {field:36s} {bv} -> {cv}  ({pct})")
            elif bv != cv:
                diffs.append(f"  {name:16s} {field:36s} {bv} -> {cv}")
    return diffs


def check_known_issues(base: dict, cand: dict) -> list[str]:
    out = []
    for issue in KNOWN_ISSUES:
        if issue["probe"] is None or issue["sample"] is None:
            continue
        s = cand.get("samples", {}).get(issue["sample"])
        if s is None:
            out.append(f"  [{issue['id']}] 樣本 {issue['sample']} 不在候選結果中")
            continue
        try:
            now = issue["probe"](s)
        except (KeyError, TypeError):
            out.append(f"  [{issue['id']}] 無法評估（欄位缺失）")
            continue
        was = issue["baseline_value"]
        if now == was:
            out.append(f"  [{issue['id']}] 仍存在，未變化 ({now})  — {issue['desc']}")
        elif was is not None and abs(now) < abs(was):
            out.append(f"  [{issue['id']}] 改善 {was} -> {now}  — {issue['desc']}")
        else:
            out.append(f"  [{issue['id']}] 惡化/變動 {was} -> {now}  — {issue['desc']}")
    return out


def report_cohort(base: dict) -> None:
    t = base["cohort_totals"]
    ranked = t.get("ranked_units") or 0
    tied = t.get("best_tree_tied_units") or 0
    print("cohort 概況（baseline）:")
    print(f"  W_total        : {base['cohort_W_total']}")
    print(f"  groups_total   : {t.get('groups_total')}")
    print(f"  ranked_units   : {ranked}")
    print(f"  tied / ranked  : {tied}/{ranked} = {tied / ranked * 100:.1f}%" if ranked else "  tied/ranked: n/a")
    print(f"  abstain        : {t.get('objective_abstain_units')}")
    print(f"  runtime (s)    : {t.get('topology_runtime_seconds')}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    ap.add_argument("--candidate", type=Path, help="新結果（與 baseline 同 schema）")
    ap.add_argument("--self-test", action="store_true", help="以 baseline 自身當候選，驗證器自檢")
    ap.add_argument("--tolerance", type=float, default=0.0, help="相對容差，預設 0（完全相等）")
    args = ap.parse_args()

    base = load(args.baseline)
    if args.self_test:
        cand = base
        print(f"=== SELF-TEST（baseline 自身）===\n")
    elif args.candidate:
        cand = load(args.candidate)
        print(f"=== 驗證 {args.candidate} ===\n")
    else:
        ap.error("需要 --candidate 或 --self-test")

    print(f"baseline: {args.baseline}")
    for name, prov in base.get("provenance", {}).items():
        print(f"  source {name}: sha256={prov['sha256'][:16]}… ({prov['size_bytes']} bytes)")
    print()
    report_cohort(base)
    print()

    exit_code = 0

    print("── L1 守恆檢查 " + "─" * 50)
    fails = check_conservation(cand.get("samples", {}))
    if fails:
        print(f"❌ {len(fails)} 項守恆違反：")
        for f in fails:
            print(f)
        exit_code = 2
    else:
        n = len(cand.get("samples", {})) * len(CONSERVATION)
        print(f"✅ 全過（{len(cand.get('samples', {}))} 樣本 × {len(CONSERVATION)} 規則 = {n} 項）")

    print("\n── L2 baseline 比對 " + "─" * 45)
    diffs = compare_to_baseline(base, cand, args.tolerance)
    if diffs:
        print(f"⚠ {len(diffs)} 項差異（容差 {args.tolerance}）：")
        for d in diffs[:40]:
            print(d)
        if len(diffs) > 40:
            print(f"  … 另有 {len(diffs) - 40} 項")
        exit_code = max(exit_code, 1)
    else:
        print("✅ 與 baseline 完全一致")

    print("\n── L3 已知問題追蹤 " + "─" * 46)
    notes = check_known_issues(base, cand)
    if notes:
        for n in notes:
            print(n)
    else:
        print("（無可自動評估的已知問題）")

    print("\n" + "=" * 64)
    print({0: "PASS", 1: "DIFF — 需人工確認", 2: "FAIL — 守恆違反，必有 bug"}.get(exit_code, "?"))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
