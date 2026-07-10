#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""capped 密集區「聚集 + 衝突檢查」重新框架（2026-07-10 使用者提案）。

背景：capped(676) 是 partial-heavy 高 k 密集區。完整枚舉所有最小 Steiner 樹在組合上爆炸，
  且 V4/V5 oracle 本身用冪集、對這些區也不可行（無法 brute-force 驗證完整）。
使用者洞見：這些區不需完整列出所有可能樹——**只需確認覆蓋的共現是否有衝突（four-gamete）**。
  無衝突 ⟺ perfect phylogeny 存在（Gusfield 定理）= 拓撲乾淨巢狀、僅 order/completion 組合多。

本腳本對每個 capped 區：① 蒐集全部觀測 pattern（full+partial）② four-gamete 衝突檢查（O(k²·patterns)，可驗證）
  → 分類 compatible（無衝突·乾淨巢狀存在）vs conflict（有衝突·需 homoplasy）。
只讀不寫源碼。輸出 capped_cluster_conflict.json。用法：python3 capped_cluster_conflict_check.py
"""
import json, os, glob
from collections import Counter
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.normpath(os.path.join(HERE, "..", "data"))
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"


def four_gamete(patterns, k):
    """對每對位點：有覆蓋(非X)的 read 中若 RR/RA/AR/AA 四型全現 → 衝突(incompatible)。回 (衝突對數, 可測對數)。"""
    conf = 0; testable = 0
    for i, j in combinations(range(k), 2):
        seen = set()
        for pat in patterns:
            a, b = pat[i], pat[j]
            if a in "RA" and b in "RA":
                seen.add(a + b)
        if seen:
            testable += 1
        if len(seen) == 4:
            conf += 1
    return conf, testable


def main():
    lay = json.load(open(f"{PILOT}/layered_reconstruction_HCC1395.json"))
    capped = [u for u in lay["detail"] if u.get("is_lineage") and "capped" in u["L1_class"]]
    mreg = {}
    for f in sorted(glob.glob(f"{PILOT}/mlhp_part_*.json")):
        for g in json.load(open(f)).get("groups", []):
            mreg[f"{g['chrom']}:{g['start']}-{g['end']}"] = g

    total = 0; compatible = 0; conflict = 0
    distinct = Counter(); conflict_regions = []
    for u in capped:
        g = mreg.get(u["region"])
        if not g:
            continue
        fam = u["family"]
        pops = (g.get("populations_by_hp", {}) or {}).get(fam, {}) or {}
        subs = (g.get("subread_groups_by_hp", {}) or {}).get(fam, {}) or {}
        k = len(g.get("positions", []))
        if k < 2:
            continue
        total += 1
        pats = list(set(pops) | set(subs))
        distinct[min(len(pats), 8)] += 1
        c, _ = four_gamete(pats, k)
        if c == 0:
            compatible += 1
        else:
            conflict += 1
            if len(conflict_regions) < 20:
                conflict_regions.append({"region": u["region"], "fam": fam, "k": k, "n_conflict_pairs": c})

    out = {
        "n_capped_analyzed": total,
        "compatible(無衝突·perfect-phylogeny 存在)": compatible,
        "compatible_pct": round(compatible / total * 100, 1) if total else 0,
        "conflict(有衝突·需 homoplasy)": conflict,
        "conflict_pct": round(conflict / total * 100, 1) if total else 0,
        "distinct_pattern_dist": {str(k if k < 8 else "8+"): v for k, v in sorted(distinct.items())},
        "conflict_region_examples": conflict_regions,
        "method": "four-gamete 相容性（Gusfield：無衝突 ⟺ perfect phylogeny 存在）；O(k²·patterns) 可完整驗證，非枚舉",
        "interpretation": "compatible = 拓撲乾淨巢狀、樹存在，只是 order/completion 組合多（無必要完整枚舉）；conflict = 需 recurrence（多為 CN artifact，同 recurrence 類送 CN 裁決）",
    }
    json.dump({"summary": out}, open(f"{D}/capped_cluster_conflict.json", "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
