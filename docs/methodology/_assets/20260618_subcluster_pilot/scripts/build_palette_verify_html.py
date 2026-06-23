#!/usr/bin/env python3
"""[色盤修正驗證 HTML] 修正前(舊·錯色,已部署 FULL dashboard figs) vs 修正後(canonical, fixed phylo_cpp_render)。
少量(~10)真實圖 → 依規約3 用 inline base64(單檔可攜)。legend 用 canonical 色(示範已修)。供肉眼確認後再全部重跑。"""
import base64, json, os

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
pairs = json.load(open("/tmp/pcr_verify/pairs.json"))


def b64(p):
    return base64.b64encode(open(p, "rb").read()).decode() if os.path.exists(p) else ""


# canonical legend swatches (示範: 用對的顏色)
LEG = [("HP1", "#60a5fa"), ("HP1-1", "#1e3a8a"), ("HP2", "#c084fc"), ("HP2-1", "#6b21a8"),
       ("ALT", "#dc2626"), ("REF", "#fbbf24"), ("群1", "#db2777"), ("群2", "#0d9488")]
leg_html = "".join(f'<span class="sw"><i style="background:{c}"></i>{n}</span>' for n, c in LEG)

rows = []
for p in pairs:
    rows.append(f"""<div class="card">
  <div class="hd"><b>{p['chrom']}:{p['snv']}</b> <span class="tag">{p['set']}</span> <span class="tag">coarse={p['coarse']}群</span></div>
  <div class="pair">
    <figure><figcaption class="bad">修正前 · 舊(錯色: HP1=紅撞ALT)</figcaption><img src="data:image/png;base64,{b64(p['old'])}"></figure>
    <figure><figcaption class="good">修正後 · canonical(HP1=藍/HP2=紫, 群=teal-pink)</figcaption><img src="data:image/png;base64,{b64(p['new'])}"></figure>
  </div></div>""")

HTML = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>色盤修正驗證 — before/after</title><style>
:root{{--bg:#0f1115;--card:#1a1d24;--bd:#2a2f3a;--tx:#e6e8ec;--mut:#9aa3b2;--acc:#d97757;--bad:#f87171;--good:#34d399}}
body{{margin:0;background:var(--bg);color:var(--tx);font:14px/1.55 system-ui,'Noto Sans CJK TC',sans-serif}}
.wrap{{max-width:1500px;margin:0 auto;padding:18px}}
h1{{font-size:20px;margin:0 0 4px}}.sub{{color:var(--mut);font-size:13px;margin-bottom:10px}}
.banner{{background:#13261c;border:1px solid var(--good);border-radius:8px;padding:10px 14px;margin:10px 0;font-size:13px}}
.legend{{display:flex;flex-wrap:wrap;gap:12px;background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px 14px;margin:10px 0}}
.sw{{display:inline-flex;align-items:center;gap:5px;font-size:12.5px}}.sw i{{width:14px;height:14px;border-radius:3px;display:inline-block;border:1px solid #0006}}
.card{{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px;margin:14px 0}}
.hd{{margin-bottom:6px}}.tag{{font-size:11px;color:var(--mut);border:1px solid var(--bd);border-radius:9px;padding:2px 7px;margin-left:4px}}
.pair{{display:grid;grid-template-columns:1fr 1fr;gap:10px}}
figure{{margin:0}}figcaption{{font-size:12px;margin-bottom:3px;font-weight:600}}
.bad{{color:var(--bad)}}.good{{color:var(--good)}}
img{{width:100%;border-radius:5px;background:#fff;display:block}}
.ft{{color:var(--mut);font-size:12px;border-top:1px solid var(--bd);margin-top:18px;padding-top:10px}}
</style></head><body><div class="wrap">
<h1>色盤修正驗證 — 修正前(舊·錯色) vs 修正後(canonical)</h1>
<div class="sub">用 fixed <code>phylo_cpp_render.py</code> 重畫 {len(pairs)} 個代表位點(coarse≥2)，與已部署 FULL dashboard 的舊圖並列。確認顏色正確後再全部重跑(~53min)。</div>
<div class="banner">✅ <b>應看到</b>: 修正後 <b>HP 欄 = 藍/紫族</b>(不再紅)、<b>ALT 欄 = 紅/琥珀</b>、<b>C++ coarse 群欄 + 樹枝 = teal/pink</b>(不撞語意軸)。修正前 HP1 是紅色(撞 ALT 紅)。</div>
<div class="legend"><b style="font-size:12px;align-self:center">canonical 圖例:</b>{leg_html}</div>
{''.join(rows)}
<div class="ft">來源: 舊圖 = obs_ws/cpp_wg/figs/(已部署, 待重畫) · 新圖 = fixed phylo_cpp_render.py(import ism_heatmap_std 的 H.HP_COL/H.ALT_COL/H.CLUSTER_COL) · 單樣本 HCC1395 · 圖 inline base64(規約3少量真實圖)。確認 ✓ → 全基因組重跑重畫 10707 圖 + regen FULL dashboard。</div>
</div></body></html>"""

out = f"{A}/20260623_palette_fix_verification_preview.standalone.html"
open(out, "w").write(HTML)
print(f"wrote {os.path.basename(out)}  ({len(HTML)} bytes, {len(pairs)} pairs)")
