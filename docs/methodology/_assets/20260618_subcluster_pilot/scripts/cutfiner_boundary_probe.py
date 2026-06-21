#!/usr/bin/env python3
"""檢視「要不要切更細」邊界的合理性：join 兩份既有結果——
  cluster_redesign.json   : 每 k 的對齊式 confidence(CONFIRMED/STRUCTURE/...)  ← 現行「切更細」判準
  kprofile_stability.json : 每 k 的 chance-corrected stability-excess(real vs within-1-group null) ← 「是否真實」
對每個 k 攤開：『真實(excess>門檻)?』×『可歸因(對齊可靠)?』2×2，找出兩判準分歧處 = 邊界脆弱點。
純讀 json，不重算。"""
import json
from collections import Counter
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
EXC=0.10  # stability-excess 門檻（真實結構，扣 null 後）

red={r["key"]:r for r in json.load(open(f"{A}/cluster_redesign.json"))["loci"]}
stab={(r["chrom"]+"_"+r["pos"]):r for r in json.load(open(f"{A}/kprofile_stability.json"))["loci"]}

cells=Counter()  # (real?, attributed?) -> count of (locus,k)
rows=[]
for key,R in red.items():
    S=stab.get(key)
    if not S: continue
    exc_by_k={p["k"]:p["stab_excess"] for p in S["per_k"]}
    conf_by_k={res["k_core"]:res["confidence"] for res in R["new_resolutions"]}
    for k in sorted(set(exc_by_k)|set(conf_by_k)):
        exc=exc_by_k.get(k); conf=conf_by_k.get(k,"not_listed")
        real = (exc is not None and exc>=EXC)
        attributed = (conf=="CONFIRMED")
        cells[(real,attributed)] += 1
        # 只記分歧 / 關鍵情形
        if (real and not attributed) or (attributed and exc is not None and exc<EXC):
            rows.append((key,k,exc,conf,real,attributed))

print("=== 2×2  (locus,k) 計數：『真實 excess≥0.10』 × 『可歸因 CONFIRMED』===")
print(f"  真實 & 可歸因      (該切、且知道為何): {cells[(True,True)]}")
print(f"  真實 & 不可歸因    (該切、但novel/單樣本無力歸因): {cells[(True,False)]}")
print(f"  不真實 & 可歸因    (對齊說CONFIRMED但扣null後不真→邊界假陽性?): {cells[(False,True)]}")
print(f"  不真實 & 不可歸因  (不該切): {cells[(False,False)]}")
print()
print("=== 分歧樣本（real但unattributed，或 CONFIRMED但excess<0.10）前 18 ===")
for key,k,exc,conf,real,attr in sorted(rows,key=lambda x:-(x[2] or -9))[:18]:
    flag = "REAL_unattributed" if (real and not attr) else "CONFIRMED_but_low_excess"
    print(f"  {key:16s} k={k} excess={exc} conf={conf:10s} → {flag}")
json.dump({"cells":{f"real={r},attr={a}":c for (r,a),c in cells.items()},"divergences":len(rows)},
          open(f"{A}/cutfiner_boundary_probe.json","w"),indent=1)
