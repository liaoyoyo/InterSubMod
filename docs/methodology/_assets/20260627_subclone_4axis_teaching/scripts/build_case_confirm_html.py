#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HCC1395 清楚案例綜合確認 HTML(2026-06-30):chr17 旗艦案例 = 演化樹 + 甲基熱圖(read×CpG) +
距離熱圖(read×read,UPGMA序) + 分類標籤條(genotype/lineage/HP/甲基cluster) + 漏斗 + residual spotlight。
數據全從 chr17_tree_data.json + chr17_subclone_data.json + methyl_auxiliary_annotation.json 注入(§13-A)。
一次觀察確認「甲基對齊 genotype-cis 非 lineage」。
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "20260630_HCC1395_case_confirmation.standalone.html"))
t = json.load(open(os.path.join(DATA, "chr17_tree_data.json")))
sd = json.load(open(os.path.join(DATA, "chr17_subclone_data.json")))
aux = json.load(open(os.path.join(DATA, "methyl_auxiliary_annotation.json")))["summary"]["b_hidden_substructure"]
reads = sd["reads"]; cpgs = sd["cpgs"]
S = ["48362515", "48365089", "48365161"]


def bit(v):
    return "A" if v == "ALT" else ("R" if v == "REF" else "?")


for r in reads:
    r["gt"] = "".join(bit(r["geno"].get(s)) for s in S)
    r["lin"] = r.get("lineage", "?").split("_")[0]


def betacol(b):
    if b is None or (isinstance(b, float) and np.isnan(b)):
        return "#e9ecef"
    b = max(0.0, min(1.0, b))
    return f"rgb({int(60+195*b)},{int(80*(1-abs(b-0.5)*2))},{int(255-195*b)})"


# ---- 甲基熱圖: reads 排序(lineage→genotype),cols=CpG ----
ordr = sorted(range(len(reads)), key=lambda i: (reads[i]["lin"], reads[i]["gt"], i))
CW = 6.5; RH = 7
linc = {"L0": "#adb5bd", "L1": "#4dabf7", "L2": "#e8590c"}
gtc = {"RRR": "#adb5bd", "RAR": "#4dabf7", "AAA": "#e8590c"}
hm = ""
for row, i in enumerate(ordr):
    r = reads[i]; y = row * RH
    for c in range(len(cpgs)):
        hm += f'<rect x="{c*CW}" y="{y}" width="{CW}" height="{RH}" fill="{betacol(r["meth"][c])}"/>'
# 左標籤條: lineage / genotype / HP / methyl-cluster
labx = len(cpgs) * CW + 4
clu = t.get("clusters_k2", [])
lab = ""
for row, i in enumerate(ordr):
    r = reads[i]; y = row * RH
    lab += f'<rect x="{labx}" y="{y}" width="10" height="{RH}" fill="{linc.get(r["lin"],"#ccc")}"/>'
    lab += f'<rect x="{labx+11}" y="{y}" width="10" height="{RH}" fill="{gtc.get(r["gt"],"#ccc")}"/>'
    cv = clu[i] if i < len(clu) else 0
    lab += f'<rect x="{labx+22}" y="{y}" width="10" height="{RH}" fill="{"#7048e8" if cv==1 else "#82c91e"}"/>'
hm_h = len(ordr) * RH

# ---- 距離熱圖: UPGMA 序 ----
dm = np.array(t["distance_matrix"]); oo = t["read_order_upgma"]
dmax = float(np.nanmax(dm)) or 1.0
DC = 6.5
dh = ""
for a, i in enumerate(oo):
    for b, j in enumerate(oo):
        v = dm[i][j] / dmax
        g = int(255 * (1 - v))
        dh += f'<rect x="{b*DC}" y="{a*DC}" width="{DC}" height="{DC}" fill="rgb({g},{g},{255 if a==b else g})"/>'
dh_sz = len(oo) * DC

# ---- 群平均β + 群×群距離 ----
from itertools import combinations
gm = {}
for r in reads:
    if "?" in r["gt"]:
        continue
    gm.setdefault(r["gt"], []).append(np.array([m if m is not None else np.nan for m in r["meth"]]))
gmean = {g: round(float(np.nanmean(np.vstack(a))), 3) for g, a in gm.items()}
gn = {g: len(v) for g, v in gm.items()}


def gd(x, y):
    a = np.nanmean(np.vstack(gm[x]), axis=0); b = np.nanmean(np.vstack(gm[y]), axis=0)
    m = ~(np.isnan(a) | np.isnan(b)); return round(float(np.nanmean(np.abs(a[m] - b[m]))), 3)


go = [g for g in ["RRR", "RAR", "AAA"] if g in gm]
gdist = {(x, y): gd(x, y) for x, y in combinations(go, 2)}

# ---- axis bar ----
axc = t["axis_sig_cpg_count"]
axbar = ""
mx = max(axc.values()) or 1
for k, v in sorted(axc.items(), key=lambda x: -x[1]):
    col = "#e03131" if "α" in k or "@48365089" in k else ("#868e96" if "HP" in k else "#1971c2")
    axbar += f'<div style="display:flex;align-items:center;gap:6px;margin:2px 0"><div style="width:160px;font-size:11px">{k}</div><div style="background:{col};height:14px;width:{int(180*v/mx)}px"></div><b style="font-size:11px">{v}</b></div>'

cx = t["cross_clu2_x_lineage"]
# residual spotlight
res = [h for h in json.load(open(os.path.join(DATA, "methyl_auxiliary_annotation.json"))).get("hidden_candidates", []) if h["flag"] == "residual-candidate"]
res_same = [h for h in res if h.get("genotype") == "RR"]  # germline 也雙峰
res_html = ""
for h in res[:6]:
    rr = " 🔴germline也雙峰" if h.get("genotype") == "RR" else ""
    res_html += f'<tr><td class="mono">{h["region"]}</td><td>{h.get("genotype")}{rr}</td><td>{h.get("cn")}</td><td>{h.get("means")}</td><td>{h.get("sizes")}</td></tr>'

gmeanrow = " ｜ ".join(f'{g} β={gmean[g]}(n={gn[g]})' for g in go)
gdrow = " ｜ ".join(f'd({x},{y})={v}' for (x, y), v in gdist.items())

html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<title>HCC1395 清楚案例綜合確認</title><style>
body{{font-family:-apple-system,"Noto Sans TC","Microsoft JhengHei",sans-serif;max-width:1180px;margin:0 auto;padding:24px;color:#212529;line-height:1.6;background:#f8f9fa}}
h1{{font-size:22px}}h2{{font-size:16px;border-bottom:2px solid #dee2e6;padding-bottom:4px;margin-top:24px}}
.card{{background:#fff;border:1px solid #dee2e6;border-radius:10px;padding:16px;margin:12px 0}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:14px}}
table{{border-collapse:collapse;font-size:12px}}td,th{{border:1px solid #dee2e6;padding:4px 8px}}.mono{{font-family:monospace}}
.note{{color:#868e96;font-size:12px}}.k{{background:#fff3bf;border:1px solid #ffe066;border-radius:6px;padding:10px 14px;margin:10px 0}}
.lg{{display:inline-block;width:11px;height:11px;border-radius:2px;vertical-align:middle;margin-right:3px}}
</style></head><body>
<h1>HCC1395 清楚案例綜合確認 — chr17:48360161 旗艦</h1>
<div class="note">⭐3 單樣本。一次觀察「遺傳樹 vs 甲基」。色：甲基β 藍=低/紅=高。數據從 chr17_tree_data + chr17_subclone_data + methyl_auxiliary 注入。</div>

<div class="grid">
<div class="card"><h2>① 癌症演化樹（sSNV 共現）</h2>
<svg viewBox="0 0 300 200" width="100%" height="190">
<circle cx="60" cy="35" r="22" fill="{betacol(gmean.get('RRR'))}" stroke="#333"/><text x="60" y="39" text-anchor="middle" fill="#fff" font-size="11" font-weight="700">RRR</text>
<text x="92" y="38" font-size="11">germline 根 β={gmean.get('RRR')}</text>
<line x1="60" y1="57" x2="60" y2="83" stroke="#1565c0" stroke-width="2"/><text x="66" y="74" font-size="9" fill="#2b8a3e">+S2(α)</text>
<circle cx="60" cy="105" r="22" fill="{betacol(gmean.get('RAR'))}" stroke="#333"/><text x="60" y="109" text-anchor="middle" fill="#fff" font-size="11" font-weight="700">RAR</text>
<text x="92" y="108" font-size="11">+α祖先 β={gmean.get('RAR')}</text>
<line x1="60" y1="127" x2="60" y2="153" stroke="#1565c0" stroke-width="2"/><text x="66" y="144" font-size="9" fill="#2b8a3e">+S1,S3</text>
<circle cx="60" cy="175" r="22" fill="{betacol(gmean.get('AAA'))}" stroke="#333"/><text x="60" y="179" text-anchor="middle" fill="#fff" font-size="11" font-weight="700">AAA</text>
<text x="92" y="178" font-size="11">+β後代 β={gmean.get('AAA')}</text>
</svg>
<div class="note">線性累積。節點色=群平均甲基。</div></div>

<div class="card"><h2>② 群平均甲基 + 群×群距離</h2>
<p style="font-size:13px">{gmeanrow}</p>
<p style="font-size:13px"><b>{gdrow}</b></p>
<div class="k" style="font-size:12px">🔑 RAR↔AAA 近、都遠離 RRR → 甲基分「有無α突變」非譜系步驟（AAA 是 RAR 後代但甲基視為近）。</div>
</div>
</div>

<div class="grid">
<div class="card"><h2>③ 甲基熱圖（read × {len(cpgs)} CpG）</h2>
<div class="note"><span class="lg" style="background:#adb5bd"></span>L0/RRR <span class="lg" style="background:#4dabf7"></span>L1/RAR <span class="lg" style="background:#e8590c"></span>L2/AAA ｜ 標籤條:lineage|genotype|甲基cluster</div>
<svg viewBox="0 0 {labx+40} {hm_h}" width="100%" height="{min(hm_h,360)}">{hm}{lab}</svg>
<div class="note">每列=1 read(依 lineage→genotype 排)。RRR(上)偏藍/低、RAR+AAA(下)偏紅/高 → 甲基隨 α 突變跳;RAR vs AAA 色相近=分不出後代。</div></div>

<div class="card"><h2>④ read×read 甲基距離熱圖（UPGMA 序）</h2>
<svg viewBox="0 0 {dh_sz} {dh_sz}" width="100%" height="{min(dh_sz,340)}">{dh}</svg>
<div class="note">亮=近、暗=遠。甲基分群（k=2）對應「突變 vs germline」而非細分 lineage。</div></div>
</div>

<div class="grid">
<div class="card"><h2>⑤ 分類標籤:甲基cluster × 遺傳lineage 交叉</h2>
<table><tr><th>甲基群＼lineage</th><th>L0(根)</th><th>L1(α)</th><th>L2(α+β)</th><th>other</th></tr>
<tr><td class="mono">cluster 1</td><td>{cx.get('1',{}).get('L0',0)}</td><td>{cx.get('1',{}).get('L1',0)}</td><td>{cx.get('1',{}).get('L2',0)}</td><td>{cx.get('1',{}).get('other',0)}</td></tr>
<tr><td class="mono">cluster 2</td><td>{cx.get('2',{}).get('L0',0)}</td><td>{cx.get('2',{}).get('L1',0)}</td><td>{cx.get('2',{}).get('L2',0)}</td><td>{cx.get('2',{}).get('other',0)}</td></tr></table>
<div class="k" style="font-size:12px">🔴 甲基 cluster 1 同含 L1=19 + L2=19 → <b>甲基分群 ≠ 遺傳 lineage</b>（無法分 RAR/AAA）。cluster 2 ≈ germline 根。</div></div>

<div class="card"><h2>⑥ 甲基對齊哪個軸（顯著 CpG 數）</h2>
{axbar}
<div class="note">α(genotype)軸最強(26)、lineage 軸弱(9)、<b>HP1-vs-HP1-1(subclone)軸=0</b> → 甲基=cis-ASM 存在性偵測器,非 lineage 排序器。</div></div>
</div>

<div class="card"><h2>⑦ 全基因組漏斗（3,416 sSNV 群）</h2>
<table><tr><th>類別</th><th>群數</th><th>%</th><th>判定</th></tr>
<tr><td>unimodal(無結構)</td><td>{aux['flag_dist'].get('unimodal',0)}</td><td>{round(100*aux['flag_dist'].get('unimodal',0)/aux['n_groups_tested'],1)}%</td><td>✅ 負篩可信</td></tr>
<tr><td>bimodal→對齊 CN-gain</td><td>{aux['flag_dist'].get('cn-flagged',0)}</td><td>—</td><td>confound(multiplicity)</td></tr>
<tr><td>bimodal→對齊 HP</td><td>{aux['flag_dist'].get('hp-explained',0)}</td><td>—</td><td>confound(germline HP)</td></tr>
<tr><td>bimodal→residual(未解釋)</td><td>{aux['flag_dist'].get('residual-candidate',0)}</td><td>{round(100*aux['flag_dist'].get('residual-candidate',0)/aux['n_groups_tested'],2)}%</td><td>⚠ L3,但多為 intrinsic(見⑧)</td></tr>
<tr><td><b>確認 subclone-specific 甲基</b></td><td><b>0</b></td><td><b>0%</b></td><td>❌ 無一通過</td></tr></table></div>

<div class="card"><h2>⑧ residual spotlight — 為何 residual ≠ subclone</h2>
<table><tr><th>region</th><th>genotype</th><th>cn</th><th>雙峰means</th><th>sizes</th></tr>{res_html}</table>
<div class="k" style="font-size:12px">🔴 chr2:207617973 連 <b>germline 群「RR」（無突變）也雙峰</b> → 甲基分裂非突變驅動 = <b>intrinsic 表觀異質性，非 subclone</b>（無遺傳佐證=double-dip）。residual 中 germline-also-bimodal 數={len(res_same)}。</div></div>

<div class="card k"><h2>🎯 一句話確認</h2>
甲基 = <b>bounded-auxiliary（chr17 視覺鐵證 + 全基因組漏斗 + 跨3癌種一致）</b>：<br>
✅ 能說「無隱藏次結構」(95.8% unimodal) + 錨定「germline-根 vs 突變-衍生」(α 軸 cis-ASM)<br>
❌ 不能排序衍生(RAR/AAA 甲基相近、HP軸 0 sig CpG) / 不能確認 subclone(0 通過、residual 是 intrinsic)<br>
→ 確認需 single-cell。</div>
<div class="note">證據級 L2(單樣本 chr17 視覺 + 全基因組漏斗)。</div>
</body></html>"""
open(OUT, "w").write(html)
print("WROTE", OUT, f"({os.path.getsize(OUT)//1024}KB)")
print("group means:", gmean, "| dists:", gdist)
