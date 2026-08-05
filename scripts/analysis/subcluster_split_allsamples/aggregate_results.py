#!/usr/bin/env python3
"""彙整全樣本 result_<S>.json → 跨樣本切割總帳對照表 (markdown + json)。
零捏造: 只讀 result_*.json 真值。缺檔則該樣本標 MISSING。
用法: aggregate_results.py <RESULT_DIR>
輸出: <RESULT_DIR>/SUMMARY_all_samples.md + summary_all_samples.json
"""
import sys, json, glob, os

RESULT_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H2009", "H1437", "HCC1937", "HCC1954"]

results = {}
for f in glob.glob(f"{RESULT_DIR}/result_*.json"):
    r = json.load(open(f))
    results[r["sample"]] = r

samples = [s for s in SAMPLE_ORDER if s in results] + [s for s in results if s not in SAMPLE_ORDER]

def g(s, k):
    return results[s].get(k, "—")

lines = []
lines.append("# 全樣本 subcluster 切割總帳對照（canonical paired_full longphase-s tagged BAM, tumor-only）")
lines.append("")
lines.append(f"> 樣本數: {len(samples)} | binary: feat/summary-nreadsvalid @ 2d5692f | ISM: w=5000, BERNOULLI, C_min=3, SKIP")
lines.append("> A=cansplit(best_k≥2) / C=clean PERMANOVA(任一軸,扣dispersion) — 邏輯逐字同 permanova_clean_and_4group.py Q1")
lines.append("")

# ===== 主表: A/C/Jaccard Venn =====
lines.append("## 主表 — A/C clean-location Venn（每樣本）")
lines.append("")
lines.append("| 樣本 | N | A=cansplit | C=clean PERMANOVA | A∩C | A−C | C−A | Jaccard |")
lines.append("|------|---|-----------|-------------------|-----|-----|-----|---------|")
for s in samples:
    lines.append(
        f"| {s} | {g(s,'N')} "
        f"| {g(s,'A_cansplit')} ({g(s,'A_pct')}%) "
        f"| {g(s,'C_clean_permanova')} ({g(s,'C_pct')}%) "
        f"| {g(s,'AnC')} ({g(s,'AnC_pct')}%) "
        f"| {g(s,'A_minus_C')} ({g(s,'A_minus_C_pct')}%) "
        f"| {g(s,'C_minus_A')} ({g(s,'C_minus_A_pct')}%) "
        f"| {g(s,'jaccard')} |"
    )
lines.append("")

# ===== 切割品質總帳 =====
lines.append("## 切割品質總帳（cant vs cansplit + 分帶；分帶為依定義重建，邊界±容差）")
lines.append("")
lines.append("| 樣本 | N | insuff(n<6) | cant(切不出) | cansplit | clear(V≥.7+sig) | mid(.5-.7) | weak(.3-.5) | vweak(<.3) | degen |")
lines.append("|------|---|-------------|--------------|----------|-----------------|------------|-------------|------------|-------|")
def pct(s, k):
    N = results[s].get("N", 0)
    v = results[s].get(k, 0)
    return f"{v} ({round(100*v/N,1)}%)" if N else str(v)
for s in samples:
    lines.append(
        f"| {s} | {g(s,'N')} | {pct(s,'insuff')} | {pct(s,'cant')} | {pct(s,'cansplit')} "
        f"| {pct(s,'clear')} | {g(s,'mid')} | {g(s,'weak')} | {g(s,'vweak')} | {g(s,'degen')} |"
    )
lines.append("")

# ===== runtime =====
lines.append("## 執行資訊")
lines.append("")
lines.append("| 樣本 | driver elapsed(s) | matched significance |")
lines.append("|------|-------------------|----------------------|")
for s in samples:
    lines.append(f"| {s} | {g(s,'elapsed_s')} | {g(s,'matched_significance')} |")
lines.append("")

out_md = f"{RESULT_DIR}/SUMMARY_all_samples.md"
open(out_md, "w").write("\n".join(lines))
json.dump({s: results[s] for s in samples}, open(f"{RESULT_DIR}/summary_all_samples.json", "w"), indent=1)
print(f"[aggregate] {len(samples)} samples → {out_md}")
print("\n".join(lines))
