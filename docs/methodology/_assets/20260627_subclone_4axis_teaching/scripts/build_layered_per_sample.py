#!/usr/bin/env python3
"""
build_layered_per_sample.py — 分層工作站「可攜版」driver(2026-07-07)
複用既有 topology_workstation/ 可攜慣例(每樣本一檔 + index.html),但內容=新分層模型
(L0 HP家族 → L1 sSNV枚舉樹 → L2 CN → L3甲基;欄位=multi-HP% / region-determinacy / lineage-determinacy)。
不覆蓋舊 topology_workstation/(6 樣本尚無新資料);輸出到 layered_workstation/。

輸入:每樣本 layered_region_view_{S}.json(SM_ML_GLOB 上游 mlhp + layered driver 產)。
  HCC1395=20260618_subcluster_pilot/;6 樣本=multisample_subclone/{S}/(待 5b 跑)。
用法:python3 build_layered_per_sample.py
"""
import os, json, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "build_layered_workstation.py")
PY = "python3"
OUTDIR = os.path.normpath(os.path.join(HERE, "..", "..", "layered_workstation"))
os.makedirs(OUTDIR, exist_ok=True)

# 樣本 → region-view json 路徑(有才建;無=待 5b)
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
SAMPLES = [
    ("HCC1395", "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/layered_region_view_HCC1395.json", "SEQC2 truth·主樣本(5khz)"),
    ("COLO829", f"{MSROOT}/COLO829/layered_region_view_COLO829.json", "melanoma·外部真值"),
    ("H1437",   f"{MSROOT}/H1437/layered_region_view_H1437.json", "LUAD"),
    ("H2009",   f"{MSROOT}/H2009/layered_region_view_H2009.json", "LUAD"),
    ("HCC1395_DORADO", f"{MSROOT}/HCC1395_DORADO/layered_region_view_HCC1395_DORADO.json", "同細胞株·DORADO basecaller"),
    ("HCC1937", f"{MSROOT}/HCC1937/layered_region_view_HCC1937.json", "breast"),
    ("HCC1954", f"{MSROOT}/HCC1954/layered_region_view_HCC1954.json", "breast"),
]

def human(b):
    for u in ("B", "KB", "MB"):
        if b < 1024: return f"{b:.0f} {u}"
        b /= 1024
    return f"{b:.1f} GB"

rows = []
for name, rv, note in SAMPLES:
    if not os.path.exists(rv):
        rows.append({"name": name, "note": note, "pending": True}); continue
    out = os.path.join(OUTDIR, name + ".html")
    env = dict(os.environ, SM_RV=rv, SM_OUT=out)
    r = subprocess.run([PY, BUILD], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(f"[{name}] BUILD FAIL:\n{r.stdout.decode()[:500]}"); rows.append({"name": name, "note": note, "pending": True, "err": True}); continue
    C = json.load(open(rv, encoding="utf-8"))["census"]
    nreg = C["n_regions"]; rd = C["region_determinacy"]; mult = C["hp_multiplicity"]; L1 = C["L1"]
    nlin = L1["n_lineage_units"]
    rows.append({"name": name, "note": note, "pending": False,
                 "ssnv": C.get("U1_sSNV_somatic_total", 0), "nreg": nreg,
                 "multiHP": mult.get("2", 0), "det": rd.get("all_determined", 0),
                 "amb": rd.get("has_ambiguous", 0), "cap": rd.get("has_capped", 0),
                 "rec": rd.get("has_recurrence", 0), "nlin": nlin,
                 "avgtree": (C.get("U5_trees", {}).get("sum_ntrees_noncapped", 0) / nlin) if nlin else 0,
                 "v17": L1.get("all_V1V7_pass", False),
                 "size": os.path.getsize(out)})
    print(f"[{name}] OK {human(os.path.getsize(out))}")

# ===== index.html(新分層欄位) =====
def pctcell(v, tot, col):
    p = 100 * v / tot if tot else 0
    return f'<td class="n"><span style="color:{col}">{p:.1f}%</span><div class="bar" style="width:{max(2,p*0.6):.0f}px;background:{col}"></div><div class="note">{v:,}</div></td>'

trh = ""
for r in rows:
    if r["pending"]:
        trh += f'<tr style="opacity:.5"><td><b>{r["name"]}</b><div class="note">{r["note"]}</div></td><td colspan="8" class="note">⏳ 待 5b：跑 mlhp multilocus + layered driver（尚無 layered_region_view）</td><td class="note">—</td></tr>'
        continue
    trh += ('<tr>'
        f'<td><a href="{r["name"]}.html"><b>{r["name"]}</b></a><div class="note">{r["note"]}</div></td>'
        f'<td class="n">{r["ssnv"]:,}</td><td class="n">{r["nreg"]:,}</td>'
        + pctcell(r["multiHP"], r["nreg"], "#db2777")
        + pctcell(r["det"], r["nreg"], "#16a34a")
        + pctcell(r["amb"], r["nreg"], "#d97706")
        + pctcell(r["cap"], r["nreg"], "#9333ea")
        + f'<td class="n">{r["rec"]}</td>'
        + f'<td class="n">{r["avgtree"]:.2f}</td>'
        + f'<td class="note">{"✅" if r["v17"] else "⚠"} {human(r["size"])}</td>'
        + f'<td><a class="btn" href="{r["name"]}.html">開啟 ▶</a></td></tr>')

ndone = sum(1 for r in rows if not r["pending"])
INDEX = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>分層樹重建工作站 — 主頁(per-HP-家族 · 枚舉全最小樹)</title><style>
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Segoe UI","Noto Sans TC","Microsoft JhengHei",sans-serif;color:#212529;background:#f8f9fa;line-height:1.5}}
.wrap{{max-width:1120px;margin:0 auto;padding:22px}}h1{{font-size:22px;margin:0 0 4px}}h2{{font-size:15px;color:#1971c2;margin:20px 0 8px}}
.sub{{color:#868e96;font-size:13px;margin:0 0 14px}}table{{border-collapse:collapse;width:100%;font-size:12.5px;background:#fff;border-radius:8px;overflow:hidden}}
th,td{{padding:6px 9px;text-align:left;border-bottom:1px solid #f1f3f5}}th{{background:#f1f3f5;font-size:11px}}td.n,th.n{{text-align:right;font-variant-numeric:tabular-nums}}
.note{{color:#868e96;font-size:10.5px}}a{{color:#1971c2;text-decoration:none}}a:hover{{text-decoration:underline}}
.btn{{background:#1971c2;color:#fff;padding:4px 11px;border-radius:6px;font-size:11.5px;white-space:nowrap}}.btn:hover{{background:#1864ab}}
.card{{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:11px 15px;margin:10px 0;font-size:12px}}
.ribbon{{background:#fffbf5;border:1px solid #ffe0a3;border-radius:6px;padding:8px 12px;margin:10px 0;font-size:12px;color:#a37200}}
.bar{{display:inline-block;height:7px;border-radius:2px;margin-left:5px;vertical-align:middle}}.mono{{font-family:ui-monospace,monospace}}
</style></head><body><div class="wrap">
<h1>🧬 分層樹重建工作站 — 主頁（per-HP-家族 · 枚舉全最小樹）</h1>
<p class="sub">L0 HP家族 → L1 sSNV枚舉 → L2 CN → L3甲基 · 樹 per germline 家族分開建 · 每檔獨立可攜 · by build_layered_per_sample.py</p>
<div class="ribbon">🏷️ <b>主分母=region(7100)</b> · region 確定⟺全 germline lineage 唯一樹 · 甲基 L3 事後不 rank 樹 · 資料模型:<span class="mono">InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md</span> · {ndone}/7 樣本已建（其餘待 5b）</div>

<h2>📊 跨樣本分層比較表</h2>
<table><tr><th>樣本</th><th class="n">somatic<br>sSNV</th><th class="n">總區</th><th class="n">multi-HP%<br>(雙親代各一樹)</th><th class="n">region<br>all-determined%</th><th class="n">has-<br>ambiguous%</th><th class="n">has-<br>capped%</th><th class="n">recur</th><th class="n">avg樹/<br>lineage</th><th>V1-7·檔</th><th></th></tr>
{trh}
</table>

<div class="card" style="background:#f4f6ff;border-color:#c3d0f5"><b>🔭 讀法（新分層模型）</b><div style="line-height:1.85;margin-top:4px">
• <b>multi-HP%</b>=兩 germline 家族(1&2)都帶 mutation → 各建一樹的區比例（舊 pooled 混淆 allelic/clonal 的比例）<br>
• <b>region all-determined%</b>=該區<b>所有</b> germline lineage 都唯一樹（多-HP 需雙確定,故嚴於 lineage-unit）<br>
• <b>has-ambiguous%</b>=≥1 家族有多棵等機率最小樹（枚舉全集,「定不出來即答案」）<br>
• <b>has-capped%</b>=太密(>4 隱藏祖先)無法窮舉 → NP-hard-dense 誠實 flag<br>
• <b>avg樹/lineage</b>=每 lineage-unit 平均枚舉樹數（determined=1；ambiguous≥2）<br>
🔴 <b>分母鐵則</b>：region% 與 lineage-unit% 分母不同不可比；跨樣本只比比例不比絕對數（絕對受深度混淆）
</div></div>

<div class="card"><b>📖 欄位/單位定義</b><div class="note" style="line-height:1.85;margin-top:3px">
• <b>somatic sSNV</b>=census somatic==True（重建骨幹；HCC1395=23,810，非總 census 35,332）<br>
• <b>總區</b>=multilocus 分析群(≤8 sSNV 窗)；linkage 全 span 版另計(同批區,位置99.9%重疊)<br>
• 完整單位階層(U1 sSNV→U6 node)+ 每比例分子÷分母 → 見資料模型 spec §1-§2<br>
• 每樣本檔內：dashboard 4 表 + region 瀏覽器(逐 lineage 枚舉樹 SVG + L0-L3 判斷軌跡)+ V1-V7 驗證
</div></div>
<p class="note" style="margin-top:16px">重生：<span class="mono">python3 build_layered_per_sample.py</span> · 單樣本 HTML by build_layered_workstation.py(SM_RV/SM_OUT)</p>
</div></body></html>"""
open(os.path.join(OUTDIR, "index.html"), "w", encoding="utf-8").write(INDEX)
print(f"\nOK index.html + {ndone} 樣本 → {OUTDIR}")
