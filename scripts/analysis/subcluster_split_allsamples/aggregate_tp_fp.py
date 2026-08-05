#!/usr/bin/env python3
"""TP vs FP 切割總帳對照 — 讀 result_<S>.json (TP) + result_<S>_FP.json (FP),
產出每樣本 TP/FP 的 A=cansplit% 與 C=clean% 對照 + json。零捏造,只讀真值。
用法: aggregate_tp_fp.py <RESULT_DIR>
輸出: <RESULT_DIR>/tp_fp_comparison.json + TP_FP_comparison.md
"""
import sys, json, os

RD = sys.argv[1] if len(sys.argv) > 1 else "."
ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H2009", "H1437", "HCC1937", "HCC1954"]

def load(p):
    return json.load(open(p)) if os.path.exists(p) else None

rows = []
for s in ORDER:
    tp = load(f"{RD}/result_{s}.json")
    fp = load(f"{RD}/result_{s}_FP.json")
    if tp is None:
        continue
    row = {"sample": s, "tp": tp, "fp": fp}
    rows.append(row)

def pct(r, k):
    N = r.get("N", 0)
    return round(100 * r.get(k, 0) / N, 1) if N else None

out = []
lines = []
lines.append("# TP vs FP 切割總帳對照（canonical paired_full longphase-s, tumor-only）")
lines.append("")
lines.append("> A=cansplit(best_k≥2) / C=clean PERMANOVA(任一軸,扣dispersion)。FP region 數小者(N<100)比例為低N噪音,僅供參考。")
lines.append("")
lines.append("| 樣本 | TP N | TP A% | TP C% | TP Jac | FP N | FP A% | FP C% | FP Jac | A%差(TP−FP) | C%差(TP−FP) |")
lines.append("|------|------|-------|-------|--------|------|-------|-------|--------|-------------|-------------|")
for r in rows:
    tp, fp = r["tp"], r["fp"]
    tpA, tpC = pct(tp, "A_cansplit"), pct(tp, "C_clean_permanova")
    if fp:
        fpA, fpC = pct(fp, "A_cansplit"), pct(fp, "C_clean_permanova")
        dA = round(tpA - fpA, 1) if (tpA is not None and fpA is not None) else None
        dC = round(tpC - fpC, 1) if (tpC is not None and fpC is not None) else None
        lowN = " ⚠低N" if fp.get("N", 0) < 100 else ""
        lines.append(f"| {r['sample']} | {tp['N']} | {tpA} | {tpC} | {tp.get('jaccard')} "
                     f"| {fp['N']}{lowN} | {fpA} | {fpC} | {fp.get('jaccard')} | {dA} | {dC} |")
    else:
        lines.append(f"| {r['sample']} | {tp['N']} | {tpA} | {tpC} | {tp.get('jaccard')} | MISSING | — | — | — | — | — |")
    out.append({"sample": r["sample"],
                "tp": {k: tp.get(k) for k in ("N", "A_cansplit", "C_clean_permanova", "AnC", "jaccard", "A_pct", "C_pct", "cant", "insuff", "clear")},
                "fp": ({k: fp.get(k) for k in ("N", "A_cansplit", "C_clean_permanova", "AnC", "jaccard", "A_pct", "C_pct", "cant", "insuff", "clear")} if fp else None)})

# pooled FP (N加權) — 因個別 FP 低N,給合併視角
fp_tot_N = sum(r["fp"]["N"] for r in rows if r["fp"])
fp_tot_A = sum(r["fp"]["A_cansplit"] for r in rows if r["fp"])
fp_tot_C = sum(r["fp"]["C_clean_permanova"] for r in rows if r["fp"])
tp_tot_N = sum(r["tp"]["N"] for r in rows)
tp_tot_A = sum(r["tp"]["A_cansplit"] for r in rows)
tp_tot_C = sum(r["tp"]["C_clean_permanova"] for r in rows)
lines.append("")
lines.append("## 合併（pooled，所有樣本相加）")
lines.append("")
lines.append(f"- **TP pooled**: N={tp_tot_N}, A=cansplit {tp_tot_A} ({round(100*tp_tot_A/tp_tot_N,1)}%), C=clean {tp_tot_C} ({round(100*tp_tot_C/tp_tot_N,1)}%)")
if fp_tot_N:
    lines.append(f"- **FP pooled**: N={fp_tot_N}, A=cansplit {fp_tot_A} ({round(100*fp_tot_A/fp_tot_N,1)}%), C=clean {fp_tot_C} ({round(100*fp_tot_C/fp_tot_N,1)}%)")
out_summary = {
    "tp_pooled": {"N": tp_tot_N, "A": tp_tot_A, "A_pct": round(100*tp_tot_A/tp_tot_N,1), "C": tp_tot_C, "C_pct": round(100*tp_tot_C/tp_tot_N,1)},
    "fp_pooled": {"N": fp_tot_N, "A": fp_tot_A, "A_pct": round(100*fp_tot_A/fp_tot_N,1) if fp_tot_N else None, "C": fp_tot_C, "C_pct": round(100*fp_tot_C/fp_tot_N,1) if fp_tot_N else None},
    "per_sample": out,
}
json.dump(out_summary, open(f"{RD}/tp_fp_comparison.json", "w"), indent=1)
open(f"{RD}/TP_FP_comparison.md", "w").write("\n".join(lines))
print("\n".join(lines))
print(f"\n[tp_fp -> {RD}/tp_fp_comparison.json]")
