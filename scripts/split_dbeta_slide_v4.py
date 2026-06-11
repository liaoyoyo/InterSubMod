#!/usr/bin/env python3
"""Fix presenter_v4 only display issue: §2 Δβ slide overflows (svgA+svgB stacked = 940px > 648).
Split into 2 slides (情境A directional / 情境B blind-spot), renumber 15->16. Lossless SVG reuse."""
import re

F = "/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/in_progress/20260603_methyl_haplotype_story_deck_v2/presenter_v4.html"
html = open(F, encoding="utf-8").read()
secs = re.findall(r"<section\b.*?</section>", html, re.S)
assert len(secs) == 15, f"expected 15, got {len(secs)}"

idx = next(i for i, s in enumerate(secs) if "Δβ 怎麼算" in s)
fig1 = secs[idx]
svgs = re.findall(r"<svg\b.*?</svg>", fig1, re.S)
assert len(svgs) == 2, f"Δβ slide should have 2 svg, got {len(svgs)}"
svgA, svgB = svgs

slide_a = f"""<section class="slide">
  <span class="tag">§2 · Δβ 方法</span><span class="snum">NN / 16</span>
  <h1 class="title">Δβ 怎麼算：兩層平均 → 方向型（BRCA2）抓得到</h1>
  <div class="sub">同步自 method_explainer；聚合數字皆已驗證真值，read 甲基標記為示意放大</div>
  <div class="body" style="display:block">
    {svgA}
    <div class="figcap">情境 A（方向型・BRCA2 真實）：多數 CpG 一致地差 → 兩層平均後仍明顯 → Δβ=−0.122 抓到；normal-anchored cis-test 分 cis（−0.142）vs drift（−0.022）。</div>
  </div>
</section>"""

slide_b = f"""<section class="slide">
  <span class="tag">§2 · Δβ 盲點</span><span class="snum">NN / 16</span>
  <h1 class="title">Δβ 盲點：交叉型被平均抵消，盲分群 ARI 才抓得到</h1>
  <div class="sub">「篩選平均」是優點也是盲點 → 需 max|Δ| / ARI / PERMANOVA 雙鍵（任務 #19）</div>
  <div class="body" style="display:block">
    {svgB}
    <div class="figcap">情境 B（交叉型・合成示意，非 BRCA2）：前後反向 → 兩層平均抵消（Δβ≈0、Wilcoxon p=0.438 漏），但 read pattern 不同 → 盲分群 ARI=0.544 抓到。</div>
  </div>
</section>"""

new_secs = secs[:idx] + [slide_a, slide_b] + secs[idx + 1:]
assert len(new_secs) == 16

# renumber snum 01-16 / 16
out = []
for i, sec in enumerate(new_secs, start=1):
    out.append(re.sub(r'<span class="snum">[^<]*</span>', f'<span class="snum">{i:02d} / 16</span>', sec, count=1))

first = html.index("<section")
last = html.rindex("</section>") + len("</section>")
head, tail = html[:first], html[last:]
head = head.replace("共 15 張（封面＋§1–§5）", "共 16 張（封面＋§1–§5）")
open(F, "w", encoding="utf-8").write(head + "\n".join(out) + "\n" + tail)
print(f"split done: 15->16 slides; Δβ split at idx {idx} -> {idx+1}(A)/{idx+2}(B)")
print(f"sections now={len(re.findall(chr(60)+'section', head + chr(10).join(out) + tail))}")
