#!/usr/bin/env python3
"""One-off: byte-copy the 3 verified SVGs from method_explainer into v3 decks as appendix.

Guarantees SVG fidelity (no hand-transcription of coordinates/numbers).
- presenter_v3.html: 2 new appendix slides (A=FIG1 svgA+svgB, B=FIG2), denominators 12->14
- report_v3.html: 1 appendix section id=appx (svgA+svgB+FIG2)
All numbers inside SVGs are pre-verified true values (grep'd from JSON sources).
"""
import re
import sys

D = "/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/in_progress/20260603_methyl_haplotype_story_deck_v2"
ME = f"{D}/method_explainer_dbeta_vs_ism_v1.html"
PRES = f"{D}/presenter_v3.html"
REP = f"{D}/report_v3.html"

me = open(ME, encoding="utf-8").read()
svgs = re.findall(r"<svg\b.*?</svg>", me, re.S)
assert len(svgs) == 3, f"expected 3 svg, got {len(svgs)}"
svgA, svgB, svgFig2 = svgs

def sized(svg):
    # add responsive sizing without touching internal coords/numbers
    return svg.replace("<svg ", '<svg style="width:100%;height:auto;max-width:780px;display:block;margin:8px auto" ', 1)

svgA, svgB, svgFig2 = sized(svgA), sized(svgB), sized(svgFig2)

# ---------- presenter: 2 appendix slides ----------
pres = open(PRES, encoding="utf-8").read()
appx_a = (
    '\n<!-- 附錄 A · 方法圖示 FIG1 -->\n<section class="slide">\n'
    '  <span class="tag">附錄 A · 方法圖示</span><span class="snum">13 / 14</span>\n'
    '  <h1 class="title">Δβ 怎麼算：兩層平均 → 擅長方向型、漏掉交叉型</h1>\n'
    '  <div class="sub">同步自 method_explainer；聚合數字皆已驗證真值，read 甲基標記為示意放大</div>\n'
    '  <div class="body" style="display:block">\n'
    f'    {svgA}\n    {svgB}\n'
    '    <div class="figcap">情境 A（方向型・BRCA2 真實）多數 CpG 一致地差 → 平均後仍明顯 → Δβ=−0.122 抓到；'
    '情境 B（交叉型・合成示意）前後反向 → 兩層平均抵消（Δβ≈0、Wilcoxon p=0.438 漏），'
    '但 read pattern 不同 → 盲分群 ARI=0.544 抓到。「篩選平均」是優點也是盲點 → 需 max|Δ|/ARI/PERMANOVA 雙鍵（任務 #19）。</div>\n'
    '  </div>\n</section>\n'
)
appx_b = (
    '\n<!-- 附錄 B · 方法圖示 FIG2 -->\n<section class="slide">\n'
    '  <span class="tag">附錄 B · 方法圖示</span><span class="snum">14 / 14</span>\n'
    "  <h1 class=\"title\">舊 ISM（Cramér's V）vs 新 Δβ：同樣有 tag，為何舊方法抓不到</h1>\n"
    '  <div class="sub">不是缺資料（tag/normal 模組都在 binary 內），是方法學（指標設計）問題 → 任務 #19</div>\n'
    '  <div class="body" style="display:block">\n'
    f'    {svgFig2}\n'
    '    <div class="figcap">舊＝無監督分群的類別關聯（Cramér&#39;s V，混淆 germline-allelic+copy、tumor-only、無方向）；'
    '新＝方向性、somatic-controlled、normal-anchored、copy-partitioned 的 Δβ → 把 somatic-cis 訊號隔離出來並誠實標出它小。</div>\n'
    '  </div>\n</section>\n'
)
# insert after last </section> (inside .slides container, before container close)
idx = pres.rfind("</section>")
assert idx != -1
cut = idx + len("</section>")
pres = pres[:cut] + appx_a + appx_b + pres[cut:]
# bump denominators 12 -> 14 (12 snum spans + 1 nav initial text)
n_before = pres.count(" / 12</span>")
pres = pres.replace(" / 12</span>", " / 14</span>")
open(PRES, "w", encoding="utf-8").write(pres)

# ---------- report: 1 appendix section ----------
rep = open(REP, encoding="utf-8").read()
appx = (
    '\n<!-- 附錄 A · 方法圖示 -->\n'
    '<h2 id="appx"><span class="n">附A</span>附錄 A · 方法圖示（Δβ 算法 + 舊 ISM vs 新）</h2>\n'
    "<p>同步自 method_explainer。聚合數字（β/Δβ/d_cis/d_within/silhouette/ARI/Cramér's V）皆為已驗證真值，read 甲基標記為示意放大。</p>\n"
    '<h3>圖 A1 · Δβ = 兩層平均 → 擅長方向型、漏掉交叉型</h3>\n'
    f"{svgA}\n{svgB}\n"
    '<p style="font-size:13px;color:#555">情境 A（方向型・BRCA2 真實）平均後仍明顯 → Δβ=−0.122 抓到；'
    '情境 B（交叉型・合成示意）前後反向 → 兩層平均抵消（Δβ≈0、Wilcoxon p=0.438 漏），但盲分群 ARI=0.544 抓到。</p>\n'
    "<h3>圖 A2 · 舊 ISM（Cramér's V）vs 新 Δβ：為何同樣有 tag，舊方法抓不到</h3>\n"
    f"{svgFig2}\n"
    '<p style="font-size:13px;color:#555">不是缺資料（tag/normal 模組都在 binary 內），是方法學（指標設計）問題 → 任務 #19。</p>\n'
)
marker = '<div class="foot" id="prov">'
assert marker in rep, "report foot marker not found"
rep = rep.replace(marker, appx + marker, 1)
open(REP, "w", encoding="utf-8").write(rep)

print(f"OK svgA/B/Fig2 lengths: {len(svgA)}/{len(svgB)}/{len(svgFig2)}")
print(f"presenter denominators bumped: {n_before} occurrences ' / 12</span>' -> ' / 14</span>'")
print("presenter appendix slides inserted: 2 (附錄 A FIG1, 附錄 B FIG2)")
print("report appendix section inserted: 1 (id=appx)")
