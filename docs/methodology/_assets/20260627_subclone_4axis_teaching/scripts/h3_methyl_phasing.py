#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Q6: H3(unphased) reads 用甲基距離試定相 HP1 vs HP2 vs 真獨立。
理論:H3 = 沒連到 germline 而未定相;若只是定相失敗→甲基應接近 HP1 或 HP2(同染色體 germline 甲基);
      若 H3 甲基同時遠離 HP1+HP2 = 真第三群(罕見:HP1+HP2+HP3 同時)。
前提:該區 HP1 vs HP2 甲基要可分(germline-ASM 存在)才能定相。
每 read 平均 β;比 mean(H3) 距 mean(HP1)/mean(HP2)。輸出 h3_methyl_phasing.json。compute batch(§13.0)。
"""
import json, os
from collections import defaultdict
import numpy as np
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
MAPQ = 20; MINR = 4; READ_CAP = 600
det = json.load(open(f"{DATA}/topology_per_region.json"))["detail"]
# 含 H3 的區(haplotypes 帶 '3' 或 node_hp 有 H3)
targets = [r for r in det if "3" in r.get("haplotypes", "") or "H3" in str(r.get("node_hp", {}).values())]

def read_meth(a):
    try: mb = a.modified_bases
    except Exception: return None
    qr = {q: rr for q, rr in a.get_aligned_pairs(matches_only=True)}
    m = {}
    if mb:
        for k, lst in mb.items():
            if k[2] != "m": continue
            for qp, mlq in lst:
                rr = qr.get(qp)
                if rr is not None: m[rr] = mlq / 255.0
    return m

def hp_norm(a):
    if not a.has_tag("HP"): return 0
    try: return int(str(a.get_tag("HP")).split("-")[0])
    except Exception: return 0

tb = pysam.AlignmentFile(TBAM, "rb")
out = []
for r in targets:
    chrom, s = r["chrom"], r["start"]; e = int(r["region"].split("-")[-1])
    bins = {1: [], 2: [], 3: []}  # HP -> per-read mean β
    n = 0
    for a in tb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        h = hp_norm(a)
        if h not in (1, 2, 3): continue
        m = read_meth(a)
        if not m or len(m) < 3: continue
        bins[h].append(float(np.mean(list(m.values())))); n += 1
        if n >= READ_CAP: break
    n1, n2, n3 = len(bins[1]), len(bins[2]), len(bins[3])
    rec = {"region": r["region"], "cn": r["cn"], "haplotypes": r.get("haplotypes", "?"),
           "n_HP1": n1, "n_HP2": n2, "n_H3": n3}
    if n3 < MINR:
        rec["verdict"] = "H3 reads 帶甲基不足(<4)→無法測"; out.append(rec); continue
    m3 = float(np.mean(bins[3])); rec["beta_H3"] = round(m3, 3)
    if n1 >= MINR and n2 >= MINR:
        m1, m2 = float(np.mean(bins[1])), float(np.mean(bins[2]))
        sep = abs(m1 - m2)
        rec["beta_HP1"], rec["beta_HP2"], rec["HP1_HP2_sep"] = round(m1, 3), round(m2, 3), round(sep, 3)
        d1, d2 = abs(m3 - m1), abs(m3 - m2)
        if sep < 0.1:
            rec["verdict"] = "HP1≈HP2 甲基不可分(無 germline-ASM)→甲基無法定相 H3"
        elif min(d1, d2) > sep:  # H3 離兩者都比 HP1-HP2 間距遠
            rec["verdict"] = f"H3 遠離 HP1+HP2(d1={round(d1,3)},d2={round(d2,3)})→可能真第三群"
            rec["assign"] = "distinct"
        else:
            rec["assign"] = "HP1" if d1 < d2 else "HP2"
            rec["verdict"] = f"H3 甲基近 {rec['assign']}(d1={round(d1,3)},d2={round(d2,3)})→很可能只是定相失敗、屬該HP"
    elif n1 >= MINR or n2 >= MINR:
        hp = 1 if n1 >= MINR else 2
        mh = float(np.mean(bins[hp])); rec[f"beta_HP{hp}"] = round(mh, 3)
        rec["verdict"] = f"只有 HP{hp} 可比(另一 HP reads 不足);H3 距 HP{hp}={round(abs(m3-mh),3)}"
    else:
        rec["verdict"] = "HP1/HP2 reads 皆不足→無基線可比"
    out.append(rec)

from collections import Counter
testable = [x for x in out if x.get("assign")]
summary = {"n_h3_regions": len(targets), "n_with_h3_meth": sum(1 for x in out if x["n_H3"] >= MINR),
           "n_testable(HP1+HP2+H3)": len(testable),
           "assign_dist": dict(Counter(x.get("assign", x["verdict"][:12]) for x in out)),
           "verdict_summary": dict(Counter(
               ("近HP(定相失敗)" if "很可能只是定相失敗" in x["verdict"]
                else "真第三群" if "真第三群" in x["verdict"]
                else "HP不可分" if "不可分" in x["verdict"]
                else "資料不足") for x in out)),
           "note": "近HP=H3 其實是 HP1/HP2 只是沒定相(支持用戶理論);真第三群=H3 甲基遠離兩HP(罕見);HP不可分=該區無 germline-ASM 甲基無法定相。"}
json.dump({"summary": summary, "regions": out}, open(f"{DATA}/h3_methyl_phasing.json", "w"), ensure_ascii=False, indent=1)
print("H3 PHASING DONE", json.dumps(summary, ensure_ascii=False))
