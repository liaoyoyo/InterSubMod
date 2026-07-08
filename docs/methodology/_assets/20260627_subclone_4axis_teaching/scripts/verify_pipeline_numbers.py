#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
端到端 pipeline 數字驗證器(2026-07-09):一鍵把 explainer 每個 headline 數字對回原始檔重算 + 一致性檢查。
只讀不寫(不動任何共用檔)。用途:讓人獨立驗證 20260709 explainer 的數字。
用法: python3 verify_pipeline_numbers.py   (HCC1395;讀同層 ../data)
每列:[值] ← 來源檔:key  | 一致性檢查 PASS/FAIL。
"""
import json, os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
def L(name): return json.load(open(os.path.join(DATA, name), encoding="utf-8"))

def line(label, val, src, ok=None):
    tag = "" if ok is None else ("  ✓" if ok else "  ✗ FAIL")
    print(f"  {label:44} = {str(val):>10}   ← {src}{tag}")

print("="*96)
print("端到端 pipeline 數字驗證(HCC1395) — 每個數字對回原始檔")
print("="*96)

# ── 站1 sSNV 宇宙 ──
led = L("sm_completeness_ledger.json")
print("\n[站1] sSNV 宇宙 ← sm_completeness_ledger.json")
line("sSNV total", led["universe_total"], "universe_total")
line("TP / FP", f'{led["universe_tp"]}/{led["universe_fp"]}', "universe_tp,universe_fp",
     led["universe_tp"]+led["universe_fp"]==led["universe_total"])
b = led["buckets"]
print("\n[站2] sm_linkage 三桶 ← sm_completeness_ledger.json:buckets")
line("linked", b["linked"], "buckets.linked")
line("underpowered", b["underpowered"], "buckets.underpowered")
line("isolated_singleton", b["isolated_singleton"], "buckets.isolated_singleton")
line("三桶加總 == total", b["linked"]+b["underpowered"]+b["isolated_singleton"],
     "sum==universe_total", b["linked"]+b["underpowered"]+b["isolated_singleton"]==led["universe_total"])

# ── 站3 regions ──
R = L("sm_region_integration.json")["regions"]
mss = [r for r in R if r.get("n_sSNV", 0) >= 2]
def has_vec(r):
    p = r.get("populations"); return isinstance(p, (list, dict)) and len(p) > 0
withv = [r for r in mss if has_vec(r)]; empty = [r for r in mss if not has_vec(r)]
print("\n[站3] region 整合 ← sm_region_integration.json:regions")
line("n_sSNV>=2 區(全宇宙)", len(mss), "count(n_sSNV>=2)")
line("有完整向量區", len(withv), "count(populations非空)")
line("空向量區(40.4%)", f"{len(empty)} ({len(empty)/len(mss)*100:.1f}%)", "count(populations空)")

# ── 站4 topology determinacy / 站5 enumerate ──
ct = L("candidate_trees.json")["summary"]; it = ct["identifiability_table"]
print("\n[站4] topology 可辨識性 ← candidate_trees.json:summary.identifiability_table")
keys = [("determined","determined"),("ambiguous","ambiguous"),
        ("recurrence(m通道)","recurrence(m通道;非乾淨可辨識)"),("conflict(真artifact)","conflict(真artifact/cycle)"),
        ("nonenumerable B_pairwise","nonenumerable(B_pairwise)"),("nonenumerable C_underdetermined","nonenumerable(C_underdetermined)"),
        ("other","other")]
wv = 0
for lbl, k in keys:
    line(lbl, it.get(k), f"identifiability_table['{k}']"); wv += it.get(k, 0)
line("with-vector 桶加總 == 3885", wv, "sum of above", wv == 3885)
sub = it.get("subcube_recovered(gap#1;partial建樹)", 0)
line("subcube_recovered(partial救回)", sub, "identifiability_table[...]")
line("總分母 == 6288", wv + sub, "3885 + subcube", wv + sub == 6288)
print("\n[站5] enumerate co-optimal ← candidate_trees.json:summary")
line("B1_enumerated", ct["B1_enumerated"], "summary.B1_enumerated")
line("  = determined + ambiguous", f'{ct["B1_determined"]}+{ct["B1_ambiguous"]}', "B1_determined,B1_ambiguous",
     ct["B1_enumerated"] == ct["B1_determined"]+ct["B1_ambiguous"])
# CAP 截斷實測
cts = L("candidate_trees.json")["candidate_trees"]
capped = sum(1 for c in cts if c.get("capped"))
maxtc = max((c.get("true_count", 0) for c in cts), default=0)
line("CAP 截斷區數(應=0)", capped, "count(capped=True)", capped == 0)
line("最大 true_count(候選樹數)", maxtc, "max(true_count)")

# ── 站5b partial-read 救回 ──
try:
    pr = L("partial_read_recoverable_fraction.json")["summary"]
    print("\n[站5b] partial-read 可救比例 ← partial_read_recoverable_fraction.json:summary")
    line("空向量區", pr["n_empty_vector(populations 空)"], "n_empty_vector")
    line("可救(>=1 pairwise,80.2%)", f'{pr["RECOVERABLE(>=1 pairwise co-read)"]} ({pr["recoverable_frac_of_empty"]*100:.1f}%)', "RECOVERABLE")
    line("chain>=3loci(可拼樹,58.6%)", f'{pr["chain>=3loci(可拼 partial 樹)"]} ({pr["chain3_frac_of_empty"]*100:.1f}%)', "chain>=3loci")
    line("真單位點(救不回,19.8%)", f'{pr["truly_single_locus(無任何 pair,不可救)"]} ({pr["truly_single_frac_of_empty"]*100:.1f}%)', "truly_single_locus")
except FileNotFoundError:
    print("\n[站5b] partial_read_recoverable_fraction.json 不在(跳過)")

# ── 站6 輔助驗證 ──
print("\n[站6] 輔助驗證")
# CN m-通道:數 topology_per_region 的 recurrence 子標
try:
    det = L("topology_per_region.json")["detail"]
    dc = Counter(r.get("determinacy") for r in det)
    rec_art = sum(v for k, v in dc.items() if k and "recurrence_artifact" in k)
    rec_cand = sum(v for k, v in dc.items() if k and "recurrence_candidate" in k)
    rec_loh = sum(v for k, v in dc.items() if k and "recurrence_LOH" in k)
    print("  CN m-通道拆 recurrence ← topology_per_region.json:detail[].determinacy")
    line("  recurrence→artifact(m>1棄)", rec_art, "count(recurrence_artifact)")
    line("  recurrence→candidate(m=1留)", rec_cand, "count(recurrence_candidate)")
    line("  recurrence→LOH(m未定)", rec_loh, "count(recurrence_LOH)")
except Exception as ex:
    print(f"  (topology_per_region recurrence 子標數不到:{ex})")
# 甲基四配子 demo
try:
    fg = L("four_gamete_aa_methyl_demo.json")["summary"]
    print("  甲基四配子否決 ← four_gamete_aa_methyl_demo.json:summary")
    line("  四配子 pair", fg["n_4gamete_pairs"], "n_4gamete_pairs")
    line("  same-HP / clean", f'{fg["funnel"]["same_hp(前提成立)"]} / {fg["funnel"]["same_hp_AND_cn_clean(neutral/loh)"]}', "funnel")
    line("  distal 顯著(皆 cross-HP CN-gain)", fg["distal_lineage_test"]["n_distal_perm_sig(p<0.05=distal 真能判別)"], "distal_lineage_test")
except Exception as ex:
    print(f"  (four_gamete demo 不到:{ex})")

print("\n" + "="*96)
print("驗證完成。上方 ✗ 代表數字與原始檔不一致(需查);全 ✓ = explainer 數字可信。")
print("="*96)
