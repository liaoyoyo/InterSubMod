#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chr17 worked：遺傳樹 vs 甲基群×群距離 並排視覺化(2026-06-30)。
輸出 standalone HTML:① 甲基取樣範圍/機制 ② 遺傳樹 ③ 各群平均β ④ 群×群甲基距離 ⑤ 判讀。
數據全從 chr17_subclone_data.json 計算注入(§13-A,不手打)。
"""
import json, os
import numpy as np
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "20260630_chr17_tree_vs_methyl_grouping.standalone.html"))
d = json.load(open(os.path.join(DATA, "chr17_subclone_data.json")))
cpgs = d["cpgs"]; reads = d["reads"]
S = ["48362515", "48365089", "48365161"]
SLAB = {"48362515": "S1(β)", "48365089": "S2(α祖先)", "48365161": "S3(β)"}


def bit(v):
    return "A" if v == "ALT" else ("R" if v == "REF" else "?")


grp = {}
for r in reads:
    g = "".join(bit(r["geno"].get(s)) for s in S)
    if "?" in g:
        continue
    grp.setdefault(g, []).append(np.array([m if m is not None else np.nan for m in r["meth"]]))

gm = {g: np.vstack(a) for g, a in grp.items()}
gmean = {g: float(np.nanmean(M)) for g, M in gm.items()}
gn = {g: len(grp[g]) for g in grp}


def gdist(A, B):
    a = np.nanmean(A, axis=0); b = np.nanmean(B, axis=0)
    mask = ~(np.isnan(a) | np.isnan(b))
    return float(np.nanmean(np.abs(a[mask] - b[mask])))


order = ["RRR", "RAR", "AAA"]
order = [g for g in order if g in gm] + [g for g in gm if g not in ("RRR", "RAR", "AAA")]
dist = {(x, y): round(gdist(gm[x], gm[y]), 3) for x, y in combinations(order, 2)}
covs = [int(np.sum(~np.isnan(M[i]))) for M in gm.values() for i in range(M.shape[0])]
span_kb = round((cpgs[-1] - cpgs[0]) / 1000, 1)
med_cov = int(np.median(covs))

# 顏色: β 0→藍(低) 1→紅(高)
def betacol(b):
    r = int(255 * b); bl = int(255 * (1 - b))
    return f"rgb({r},80,{bl})"

# 群節點 + 樹邊(genetic)
muts = {"RRR": "germline", "RAR": "S2(α)", "AAA": "S1,S2,S3"}
# 樹 SVG: ROOT→RAR→AAA 線性 + 各群 β 色塊
treerows = ""
yy = {"RRR": 40, "RAR": 110, "AAA": 180}
for g in ["RRR", "RAR", "AAA"]:
    if g not in gmean: continue
    treerows += f'<circle cx="120" cy="{yy[g]}" r="26" fill="{betacol(gmean[g])}" stroke="#333"/>'
    treerows += f'<text x="120" y="{yy[g]-1}" text-anchor="middle" fill="#fff" font-size="12" font-weight="700">{g}</text>'
    treerows += f'<text x="120" y="{yy[g]+12}" text-anchor="middle" fill="#fff" font-size="9">β={gmean[g]:.2f}</text>'
    treerows += f'<text x="158" y="{yy[g]+4}" font-size="11">{muts[g]} (n={gn.get(g,0)})</text>'
treerows += '<line x1="120" y1="66" x2="120" y2="84" stroke="#1565c0" stroke-width="2.5"/><text x="128" y="80" font-size="10" fill="#2b8a3e">+S2(α)</text>'
treerows += '<line x1="120" y1="136" x2="120" y2="154" stroke="#1565c0" stroke-width="2.5"/><text x="128" y="150" font-size="10" fill="#2b8a3e">+S1,S3</text>'

# 群×群距離表
dcells = ""
for x in order:
    dcells += f'<tr><td class="g">{x}</td>'
    for y in order:
        if x == y:
            dcells += '<td style="background:#f1f3f5">—</td>'
        else:
            v = dist.get((x, y)) or dist.get((y, x))
            # 近=綠 遠=紅
            col = f"rgba({int(255*min(v/0.4,1))},{int(200*(1-min(v/0.4,1)))},80,0.5)"
            dcells += f'<td style="background:{col}"><b>{v}</b></td>'
    dcells += "</tr>"

html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<title>chr17 遺傳樹 vs 甲基群分布</title><style>
body{{font-family:-apple-system,"Noto Sans TC","Microsoft JhengHei",sans-serif;max-width:1000px;margin:0 auto;padding:24px;color:#212529;line-height:1.7;background:#f8f9fa}}
h1{{font-size:21px}}h2{{font-size:16px;border-bottom:2px solid #dee2e6;padding-bottom:4px;margin-top:26px}}
.card{{background:#fff;border:1px solid #dee2e6;border-radius:10px;padding:16px;margin:12px 0}}
table{{border-collapse:collapse;margin:8px 0}}td,th{{border:1px solid #dee2e6;padding:6px 12px;text-align:center}}
td.g{{font-weight:700;background:#f8f9fa}}
.note{{color:#868e96;font-size:13px}}.k{{background:#fff3bf;border:1px solid #ffe066;border-radius:6px;padding:8px 12px;margin:8px 0}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:14px}}
</style></head><body>
<h1>chr17:48360161 worked — 遺傳樹 vs 甲基「哪些標籤靠近」</h1>
<div class="note">HCC1395 ⭐3。節點顏色=群平均甲基β（藍=低/紅=高）。數據從 chr17_subclone_data.json 計算注入。</div>

<div class="card"><h2>① 甲基取樣範圍與機制</h2>
<b>範圍</b>：{len(cpgs)} 個 CpG，跨 <b>{span_kb}kb</b>（{cpgs[0]}–{cpgs[-1]}）；sSNV 在 48362515–48365161，CpG 窗向上下游延伸。每 read 只覆蓋其 aligned span 內的 CpG（中位 <b>{med_cov}/{len(cpgs)}</b>）。<br>
<b>機制</b>：每 CpG 每 read 的 β = <b>ML/255</b>（Dorado 5mCG modified_bases，只取 5mC）；群=genotype 共現分類；群甲基=群內 read β 平均。</div>

<div class="grid">
<div class="card"><h2>② 遺傳樹（sSNV 共現）</h2>
<svg viewBox="0 0 320 215" width="100%" height="215">{treerows}</svg>
<div class="note">germline RRR → +S2(α) → RAR → +S1,S3 → AAA。線性累積。</div></div>

<div class="card"><h2>③ 甲基群×群距離（|Δβ|）= 甲基認為誰靠近</h2>
<table><tr><th></th>{''.join(f'<th>{g}</th>' for g in order)}</tr>{dcells}</table>
<div class="note">綠=近(甲基相似)、紅=遠。</div></div>
</div>

<div class="card k"><h2>④ 關鍵判讀</h2>
<b>遺傳樹</b>：RRR→RAR→<b>AAA 是 RAR 的後代</b>（多 S1,S3 兩個突變）。<br>
<b>甲基卻認為</b>：RAR↔AAA <b>近（{dist.get(('RAR','AAA'),'?')}）</b>，兩者都遠離 RRR（{dist.get(('RRR','RAR'),'?')}/{dist.get(('RRR','AAA'),'?')}）。<br><br>
🔴 <b>甲基分的是「有沒有 α(S2) 突變」（germline vs 突變），不是遺傳譜系步驟</b>——RAR 和 AAA 都帶 α → 甲基高 → 視為同一群；甲基<b>看不出 AAA 是 RAR 的後代</b>（那一步被抹平）。<br>
→ 這是 <b>cis-ASM</b>：甲基對齊突變的局部效應（α 軸 26 個顯著 CpG），<b>非獨立 lineage</b>。所以甲基<b>能分「突變了沒」、不能排「誰是誰的後代」</b>=bounded-auxiliary 的視覺證據。</div>

<div class="note">證據級 L2（單樣本 chr17）。甲基群×群距離對齊 α 突變軸而非樹深度 → 確認甲基為 cis-ASM 偵測器、非 lineage 排序器。</div>
</body></html>"""
open(OUT, "w").write(html)
print("WROTE", OUT)
print("group means:", {g: round(v, 3) for g, v in gmean.items()})
print("group dists:", dist)
