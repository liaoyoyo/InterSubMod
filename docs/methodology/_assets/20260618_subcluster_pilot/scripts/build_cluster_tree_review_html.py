#!/usr/bin/env python3
"""[演化樹配色 強對比版 review HTML] 清楚定義 + 色票 + 5 深樹 mockup(舊扁平 vs 新)。
色票由 cluster_tree_palette 即時算。inline base64。"""
import base64, json, os, sys
sys.path.insert(0, os.path.dirname(__file__)); import cluster_tree_palette as CT

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
mock = json.load(open("/tmp/cluster_tree_mockup/mockup.json"))


def b64(p):
    return base64.b64encode(open(p, "rb").read()).decode() if os.path.exists(p) else ""


def swrow(items):
    return "".join(f'<span class="sw"><i style="background:{c}"></i>{n}</span>' for n, c in items)


CC = CT.cluster_tree_color
left = [(p, CC(p)) for p in ["1-1", "1-1-1", "1-1-2", "1-1-1-1", "1-1-1-2"]]       # root1 左支 magenta 家族
right = [(p, CC(p)) for p in ["1-2", "1-2-1", "1-2-2"]]                            # root1 右支 cyan 家族
root2 = [(p, CC(p)) for p in ["2", "2-1"]]
exc = [("other", CC("other")), ("outlier", CC("outlier"))]
HP = [("HP1", "#60a5fa"), ("HP1-1", "#1e3a8a"), ("HP2", "#c084fc"), ("HP2-1", "#6b21a8"), ("HP3", "#14b8a6")]
RES = [("ALT", "#dc2626"), ("REF", "#fbbf24"), ("tumor", "#f97316"), ("normal", "#22c55e")]

cards = "".join(f"""<div class="card"><div class="hd"><b>{m['title']}</b> <span class="tag">標籤 {', '.join(m['terms'])}</span></div>
<img src="data:image/png;base64,{b64(m['fig'])}"></div>""" for m in mock)

HTML = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>演化樹配色(強對比) — 定義 + review</title><style>
:root{{--bg:#0f1115;--card:#1a1d24;--bd:#2a2f3a;--tx:#e6e8ec;--mut:#9aa3b2;--acc:#d97757;--good:#34d399}}
body{{margin:0;background:var(--bg);color:var(--tx);font:14px/1.6 system-ui,'Noto Sans CJK TC',sans-serif}}
.wrap{{max-width:1500px;margin:0 auto;padding:18px}}
h1{{font-size:20px;margin:0 0 4px}}h2{{font-size:15px;border-left:4px solid var(--acc);padding-left:8px;margin:20px 0 8px}}
.sub{{color:var(--mut);font-size:13px}}
.def{{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:12px 16px;margin:8px 0}}
.def ol{{margin:6px 0 6px 18px;padding:0}}.def li{{margin:4px 0}}.def code{{background:#0f1115;padding:1px 5px;border-radius:4px;font-size:12px}}
.legend{{display:flex;flex-wrap:wrap;gap:11px;background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px 14px;margin:6px 0;align-items:center}}
.sw{{display:inline-flex;align-items:center;gap:5px;font-size:12px}}.sw i{{width:16px;height:16px;border-radius:3px;display:inline-block;border:1px solid #0006}}
.ok{{background:#13261c;border:1px solid var(--good);border-radius:8px;padding:9px 14px;margin:8px 0;font-size:13px}}.ok b{{color:var(--good)}}
.warn{{background:#2a2410;border:1px solid #d97706;border-radius:8px;padding:9px 14px;margin:8px 0;font-size:13px}}.warn b{{color:#fbbf24}}
.card{{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px;margin:12px 0}}
.hd{{margin-bottom:6px}}.tag{{font-size:11px;color:var(--mut);border:1px solid var(--bd);border-radius:9px;padding:2px 7px;margin-left:4px}}
img{{width:100%;border-radius:5px;background:#fff;display:block}}
.ft{{color:var(--mut);font-size:12px;border-top:1px solid var(--bd);margin-top:18px;padding-top:10px}}
b.k{{color:var(--acc)}}
</style></head><body><div class="wrap">
<h1>演化樹 (phylo cluster) tree-aware 配色 — 強對比版</h1>
<div class="sub">取代扁平 CLUSTER_COL。HP 不動。本版把<b>第一次分裂改成跳到不同色相錨點(magenta vs cyan)</b>,讓兩大支一眼分得開。</div>

<h2>① 清楚定義(SoT, <code>cluster_tree_palette.py</code>)</h2>
<div class="def">
<b>標籤語意</b>: <code>coarse_label</code>=遞迴二元樹路徑。<code>segs=label.split('-')</code>;<code>segs[0]</code>=根群,之後每位∈{{1,2}}=該層左/右子;深度=segs 長度。例 <code>1-2-1</code>=根1→右子→左子(深3)。
<br><b>配色三維(全 deterministic,同 label 跨 region 顏色固定)</b>:
<ol>
<li><b class="k">色相 HUE</b> = <b>第一次分裂(segs[1])直接跳到不同色相錨點</b>(明顯對比,非同色深淺):根1 左子=<b>magenta</b> / 右子=<b>cyan</b>。更深的分裂(segs[2:])=在錨點鄰域內<b>微幅 hue 漂移</b>(左−/右+,逐層減半)→ 表親可分、仍歸該大支。</li>
<li><b class="k">亮度 L</b> = 隨<b>深度遞減</b>(深=暗);冗餘且 luminance 故<b>色盲安全</b>;兄弟 -2 微亮。</li>
<li><b class="k">飽和 S</b> = 固定 {CT.SAT}(鮮明)。</li>
<li><b class="k">例外</b>(other/outlier/edge) = 中性<b>灰 #9ca3af</b>。</li>
</ol>
<b>錨點色相</b>: 根1 {{左 magenta {CT.ANCHOR[('1','1')]}, 右 cyan {CT.ANCHOR[('1','2')]}}}、根2 {{rose / violet-magenta}}。
<br><b>護欄</b>: 錨點只用 cluster-reserved {{magenta, rose, cyan}},脫離 HP 藍/紫/teal + ALT 紅/REF 琥珀/T-N 橘綠。
</div>

<div class="warn">⚠ <b>對比天花板(誠實)</b>: 扣掉 HP(藍/紫/teal)+保留軸(紅/琥珀/橘/綠)後,真正「明顯不同」的色相只剩 <b>~3 區(magenta / rose / cyan)</b>。兩大支用 magenta vs cyan = 最大對比;3 支以上靠 rose + 亮度階。<b>要再多色就得放寬護欄</b>(讓 cluster 借用 HP 的藍/紫,接受與 HP 欄較近)。唯一近鄰: cyan ↔ HP3-teal(HP3 罕見、不同欄)。</div>

<h2>② 色票(即時算)</h2>
<div class="sub">根1 <b>左支 1-1-*</b>(magenta 家族, 深=暗):</div><div class="legend">{swrow(left)}</div>
<div class="sub">根1 <b>右支 1-2-*</b>(cyan 家族):</div><div class="legend">{swrow(right)}</div>
<div class="sub">根2(rose) + 例外(灰):</div><div class="legend">{swrow(root2 + exc)}</div>
<div class="sub" style="margin-top:8px">撞色檢查 — HP(不動): <span class="legend" style="display:inline-flex;padding:5px 9px">{swrow(HP)}</span></div>
<div class="sub">保留軸: <span class="legend" style="display:inline-flex;padding:5px 9px">{swrow(RES)}</span></div>

<div class="ok">✅ 第一次分裂 1-1 vs 1-2 = <b>magenta vs cyan</b>(一眼分得開);子層用亮度+微漂顯示深度與兄弟;根2=rose。與 HP(藍/紫)區隔。</div>

<h2>③ 5 個深樹 mockup — 左 舊扁平 vs 右 新強對比(HP 欄不動=撞色檢查)</h2>
<div class="sub">右側: 上半 magenta(1-1 支) / 下半 cyan(1-2 支), 一眼兩大群; 各支內深度=暗淺。</div>
{cards}

<div class="ft">來源: cluster_tree_palette.py(定義+色) / render_cluster_tree_mockup.py(圖) / workflow w9bt81b3l(裁決) · 單樣本 HCC1395 · inline base64。
<br>🔸 還不夠明顯 → 我可 (a)放寬護欄借更多色相(接受與 HP 較近) (b)再加 rose/第4色 (c)調飽和/亮度。確認 OK → 改 SoT + 全部重畫。</div>
</div></body></html>"""

out = f"{A}/20260623_cluster_tree_palette_review.standalone.html"
open(out, "w").write(HTML)
print(f"wrote {os.path.basename(out)} ({len(HTML)} bytes, {len(mock)} mockups)")
