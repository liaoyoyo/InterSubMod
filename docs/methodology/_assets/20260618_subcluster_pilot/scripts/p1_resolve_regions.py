#!/usr/bin/env python3
"""[P1 候選 region dir 解析] cis_candidates.json 的每個 (chrom,pos,set) → 確定構造 region dir 路徑
(window = anchor±5000) → 檢查 methylation.csv + reads.tsv 存在 → land cis_candidates_resolved.json。
任何 cis-test 設計都需此映射；驗證矩陣可達 + 抓 missing。純讀。"""
import json, os, glob

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
KEEP = f"{WT}/output/_phylo_wg_full"
cands = json.load(open(f"{A}/cis_candidates.json"))


def region_dir(ch, pos, s):
    # records_v6 pos = window start = anchor-5000; SNV anchor = pos+5000; window = pos..pos+10000
    anchor = pos + 5000
    sl = s.lower()
    return f"{KEEP}/{ch}_{s}/filtered_snv_{sl}_{ch}/{ch}/{ch}_{anchor}/{ch}_{pos}_{pos + 10000}"


def resolve(c):
    ch, pos, s = c["chrom"], int(c["pos"]), c["set"]
    rd = region_dir(ch, pos, s)
    if os.path.exists(f"{rd}/methylation/methylation.csv"):
        return rd
    # fallback: glob by window-start = pos (robust to anchor/edge variation)
    g = glob.glob(f"{KEEP}/{ch}_{s}/filtered_snv_{s.lower()}_{ch}/{ch}/{ch}_*/{ch}_{pos}_*/methylation/methylation.csv")
    return os.path.dirname(os.path.dirname(g[0])) if g else None


found = 0
missing = []
for c in cands:
    rd = resolve(c)
    c["region_dir"] = rd
    ok = rd is not None and os.path.exists(f"{rd}/reads/reads.tsv")
    c["resolved"] = ok
    if ok:
        found += 1
    else:
        missing.append(f"{c['chrom']}:{c['pos']}:{c['set']}")
print(f"resolved {found}/{len(cands)}")
print(f"missing ({len(missing)}) sample: {missing[:15]}")
# 解析率分 TP/FP
tp = [c for c in cands if c["set"] == "TP"]
fp = [c for c in cands if c["set"] == "FP"]
print(f"  TP resolved {sum(c['resolved'] for c in tp)}/{len(tp)} | FP resolved {sum(c['resolved'] for c in fp)}/{len(fp)}")
json.dump(cands, open(f"{A}/cis_candidates_resolved.json", "w"), ensure_ascii=False, indent=1)
print("[-> cis_candidates_resolved.json]")
