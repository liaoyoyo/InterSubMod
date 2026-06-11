#!/usr/bin/env python3
"""Build presenter_v4.html: re-architect v3's verified sections into 5 focused sections.

User decisions (2026-06-04):
 - NEW v4, keep v3 untouched
 - phasing RESEARCH AXIS completely removed (drop C5; strip C7 phasing line/row; drop prov phasing row)
   (KEEP technical mentions: C1 label-flip 'phasing run', C2b external-tool '甲基→phasing' niche)
 - §4 framed as Productive-Failure methodology value

Reuses v3 section blocks losslessly (incl. IGV base64 + 3 SVGs). Authors 2 NEW slides.
All numbers are pre-verified true values (this session). No number production here.
"""
import re

D = "/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/in_progress/20260603_methyl_haplotype_story_deck_v2"
SRC = f"{D}/presenter_v3.html"
OUT = f"{D}/presenter_v4.html"

html = open(SRC, encoding="utf-8").read()
secs = re.findall(r"<section\b.*?</section>", html, re.S)
assert len(secs) == 14, f"expected 14 sections, got {len(secs)}"

def find(marker):
    hits = [s for s in secs if marker in s]
    assert len(hits) == 1, f"marker {marker!r} matched {len(hits)}"
    return hits[0]

cover = find('class="slide cover"')
c1    = find("真實 cis-candidate（BRCA2）")
c2    = find("四步可重現流程")
c2b   = find("ISM 填的是「用既有 tag")
c3    = find("其他位點稀有且需嚴格 gate")
c4    = find("為什麼難：訊號混亂")
c4b   = find("外部文獻為「困難」命名")
c6    = find("高特異、低召回、可調靈敏度")
c7    = find("其他甲基方向：機會存在")
c8    = find("論文定位：ISM 甲基工具")
prov  = find("可驗證來源、證據分級")
fig1  = find("Δβ 怎麼算：兩層平均")
fig2  = find("為何舊方法抓不到")
# c5 (phasing panorama) intentionally NOT selected -> dropped

# ---- strip phasing RESEARCH AXIS from reused sections ----
c7 = re.sub(r'\s*<li class="live">⚠ 但<b>純 phasing.*?</li>', "", c7)
c7 = re.sub(r"\s*<tr><td>純 phasing 延長 somatic read</td>.*?</tr>", "", c7, flags=re.S)
c7 = c7.replace("機會 vs 困難：甲基版受限訊號；純 phasing 版才是強路",
                "機會 vs 困難：甲基方向多受限訊號（AUC≈0.505）；列出畫清哪條值得做、哪條不夠")
c7 = c7.replace("但我要明確切開：如果是用純 phasing 或 germline-het（不靠甲基）去延長 somatic read，那是完全不同的事，難度價值都高，是我們的 LIVE 主軸——別把這兩條混在一起。", "")
assert "LIVE 主軸" not in c7 and "純 phasing" not in c7, "C7 phasing axis not fully stripped"

prov = re.sub(r"\s*<tr><td>phasing 主軸</td>.*?</tr>", "", prov, flags=re.S)
prov = prov.replace("v3 主軸＝ISM 甲基工具；phasing 整合於 C5。",
                    "v4 主軸＝ISM 甲基-haplotype ASM（純甲基焦點，phasing 軸不在本 deck）。")
assert "phasing 主軸" not in prov, "prov phasing row not removed"

# ---- cover lede -> 5-section agenda (drop 全景/phasing) ----
cover = cover.replace(
    "從一個真實個案（BRCA2）出發，說明方法、推廣、困難、全景與工具價值——稀有但真實",
    "5 段：①個案存在 ②如何找到＋方法差異 ③為何困難 ④快速驗證→哪些 possible ⑤數據驗證——稀有但真實、誠實負結果")

# ---- NEW slide: §3 why-not-detected-before ----
new_whynot = """<section class="slide">
  <span class="tag">TAG</span><span class="snum">NN / 15</span>
  <h1 class="title">為何之前偵測不到：舊方法看得到訊號，但不是 somatic 特異</h1>
  <div class="sub">分群＋tag 關聯確實存在（Cramér's V=0.80），卻被三源混淆</div>
  <div class="body">
    <div class="col-pts">
      <ul class="pts">
        <li class="info">舊 ISM 分群 ↔ tag 關聯<b>確實存在</b>（Cramér's V=0.80, p=0.0, passed gating）</li>
        <li class="dead">但<b>非 somatic 特異</b>：混 germline-allelic(HP1↔HP2)＋copy＋n_reads</li>
        <li>缺 normal-anchor＋缺方向性 → BRCA2 淹沒在<b>數千高-V 位點</b>中不突出</li>
      </ul>
      <div class="caveat">⚠ epipolymorphism 當「<b>困難真實</b>」可引（Landan 2012 / Derrien 2021）；當「<b>判別特徵</b>」是 n_reads artifact（AUC 0.845→0.530 校正後），<b>兩者勿混</b>。</div>
    </div>
    <div class="col-fig"><div class="figwrap">
      <table class="tbl">
        <tr><th>混淆源</th><th>為何看不到 somatic</th><th class="c">新方法如何隔離</th></tr>
        <tr><td>germline-allelic</td><td>HP1↔HP2 常見 ASM 也高 V</td><td class="c">HP1 vs HP1-1 同單倍型</td></tr>
        <tr><td>copy-number</td><td>LOH/CN 造成假甲基差</td><td class="c">copy-partition（d_within）</td></tr>
        <tr><td>n_reads（覆蓋）</td><td>O11 異質性 0.845→0.530</td><td class="c">coverage-rarefied</td></tr>
        <tr><td>無 normal 基線</td><td>分不出 cis vs drift</td><td class="c">normal-anchored cis-test</td></tr>
        <tr><td>類別關聯非方向</td><td>無法排序/定方向</td><td class="c">方向性 Δβ（可排序）</td></tr>
      </table>
      <div class="figcap">「偵測不到 somatic ASM」＝舊指標看得到訊號但混淆；新方法逐源 de-confound（視覺見 §2 方法差異）</div>
    </div></div>
  </div>
  <details class="note"><summary>講稿</summary><div class="script">這頁回答「為什麼之前偵測不到」。重點是：舊 ISM 的分群跟 tag 關聯其實看得到——Cramér's V 高達 0.80、通過 gating。問題不是沒訊號，是這個訊號不是 somatic 特異：它同時混了 germline-allelic（HP1 跟 HP2 之間常見的 ASM 也會高 V）、copy-number（LOH/CN 造成的假甲基差）、還有 n_reads 覆蓋（這正是 O11 epipolymorphism 的陷阱，異質性 AUC 看起來 0.845，做 n_reads 校正後掉到 0.530，是 artifact）。再加上舊輸出沒有 normal 基線、也沒有方向性，所以 BRCA2 這種真實位點就淹沒在數千個高-V 位點裡。新方法是逐一把這些混淆源拆掉：HP1 vs HP1-1 控住 germline、copy-partition 拆 copy、coverage-rarefied 殺 n_reads、normal-anchored 分 cis/drift、方向性 Δβ 可排序——這就是為何同樣的 tag 資料，新方法抓得到而舊方法漏掉。一句話：不是缺資料，是方法學。</div></details>
</section>"""

# ---- NEW slide: §4 quick-verify verdict map (Productive-Failure core) ----
new_verdict = """<section class="slide">
  <span class="tag">TAG</span><span class="snum">NN / 15</span>
  <h1 class="title">快速驗證每個甲基突破想法 → 我們先替你關掉的坑</h1>
  <div class="sub">Productive Failure：6 個問題快速驗證 → 1 稀有 anchor＋5 關閉；possible 只剩 characterization/tooling</div>
  <div class="body">
    <div class="col-pts">
      <ul class="pts">
        <li class="info">把每個「想用甲基突破」的問題都<b>快速驗證</b>，不留模糊待辦</li>
        <li class="live">方法學貢獻＝<b>把死路標清楚</b>，reviewer 會問的先自己關掉</li>
        <li class="hold">possible 只剩 <b>characterization＋tooling</b>（非 filter/subclone/two-hit）</li>
      </ul>
      <div class="caveat">誠實口徑：6 問題中只有 Q1–Q3（定位）是 PARTIAL（稀有真實）；Q3b/Q4/Q5/Q6 全 NEGATIVE/DEAD；filter（Q6）已 concluded，<b>不列待辦</b>。</div>
    </div>
    <div class="col-fig"><div class="figwrap">
      <table class="tbl">
        <tr><th>想用甲基突破的問題</th><th class="c">快速驗證</th><th class="c">依據</th></tr>
        <tr><td>Q1 哪些 somatic 位點有甲基差異</td><td class="c">🟡 PARTIAL</td><td class="c">BRCA2＋regime-A 15/23；0.34%</td></tr>
        <tr><td>Q2 甲基差異與哪種 tag 關聯</td><td class="c">🟡 PARTIAL</td><td class="c">HP-axis modest；多為 copy</td></tr>
        <tr><td>Q3 定位明顯差異位點</td><td class="c">🟡 定位</td><td class="c">regime-A 可輸出（稀有）</td></tr>
        <tr><td>Q3b 輔助 phase / tag / HP3 / unphase</td><td class="c">🔴 NEG</td><td class="c">甲基判別 AUC=0.505</td></tr>
        <tr><td>Q4 分析 subclone（多分群/演化樹）</td><td class="c">🔴 NEG</td><td class="c">不過 germline-het baseline</td></tr>
        <tr><td>Q5 二次打擊＋發生順序</td><td class="c">🔴 NEG</td><td class="c">region-level 無 two-hit 訊號</td></tr>
        <tr><td>Q6 甲基提升 F1（降 FP）</td><td class="c">⚫ DEAD</td><td class="c">LOSO 100% circular（held-out −0.00012）</td></tr>
      </table>
      <div class="figcap">6 問題＝1 稀有 anchor（Q1–3 PARTIAL）＋5 關閉 → 把每條死路都驗證並標清＝productive-failure 方法學貢獻</div>
    </div></div>
  </div>
  <details class="note"><summary>講稿</summary><div class="script">這是這份報告的核心誠實段。我們把每一個「想用甲基去突破」的問題都做了快速驗證，而不是留一堆模糊待辦。結果：Q1 哪些 somatic 位點有甲基差異——PARTIAL，BRCA2 加 regime-A 過嚴格 gate 只剩 15/23，全基因組 strong-ASM 才 0.34%，稀有；Q2 跟哪種 tag 關聯——PARTIAL，HP-axis 有但 modest，copy-partition 顯示多數其實是 copy；Q3 定位明顯位點——可以輸出，但稀有；Q3b 想用甲基輔助 phase、tag、HP3、unphase——NEGATIVE，甲基判別力全域 AUC 只有 0.505，接近隨機；Q4 subclone 多分群、演化樹——NEGATIVE，過不了 germline-het baseline；Q5 二次打擊跟發生順序——NEGATIVE，region-level 沒有訊號；Q6 用甲基提升 F1、降 FP——DEAD，LOSO 證實 100% circularity，已 concluded，不再開。所以六個問題其實是一個稀有但真實的 anchor 加上五條關閉的路。這正是 productive-failure 的價值：把每條死路都驗證清楚、把 reviewer 會質疑的點先自己關掉，剩下真正 possible 的只有 characterization 跟 tooling。</div></details>
</section>"""

# ---- v4 order + §N tags ----
order = [cover, c1, c2, fig1, fig2, c2b, c6, c4, new_whynot, c4b, new_verdict, c3, c7, c8, prov]
tags  = [None, "§1 · 個案存在", "§2 · 如何找到", "§2 · Δβ 方法", "§2 · 方法差異",
         "§2 · 工具地景", "§2 · 工具定位", "§3 · 為何困難", "§3 · 為何之前偵測不到",
         "§3 · 文獻誠實", "§4 · 快速驗證 verdict", "§4 · 推廣誠實", "§4 · 哪些 possible",
         "§4 · 前路 roadmap", "§5 · 數據與驗證"]
assert len(order) == len(tags) == 15

out_secs = []
for i, (sec, tag) in enumerate(zip(order, tags), start=1):
    s = re.sub(r'<span class="snum">[^<]*</span>', f'<span class="snum">{i:02d} / 15</span>', sec, count=1)
    if tag:
        s = re.sub(r'<span class="tag">[^<]*</span>', f'<span class="tag">{tag}</span>', s, count=1)
    out_secs.append(s)

# ---- reassemble ----
first = html.index("<section")
last = html.rindex("</section>") + len("</section>")
head_part, tail = html[:first], html[last:]
head_part = head_part.replace(
    "共 12 張（封面＋C1–C8＋附錄）· v3 主軸＝ISM 甲基工具",
    "共 15 張（封面＋§1–§5）· v4 主軸＝ISM 甲基-haplotype ASM（純甲基）")
head_part = head_part.replace(
    "v3 主軸＝ISM 甲基工具(C1 BRCA2 hero→C6 工具定義)，phasing 整合於 C5 全景。",
    "v4 主軸＝ISM 甲基-haplotype ASM 5 段(個案/方法/困難/verdict/數據)，phasing 軸不放。")

new_html = head_part + "\n".join(out_secs) + "\n" + tail
open(OUT, "w", encoding="utf-8").write(new_html)

# ---- self-check ----
n_sec = len(re.findall(r"<section\b", new_html))
n_slidecls = len(re.findall(r'class="slide', new_html))
axis = sum(new_html.count(x) for x in ["phasing 主軸", "phasing：n=7", "LIVE 主軸", "phasing 整合於 C5", "純 phasing"])
print(f"OUT={OUT}")
print(f"sections={n_sec}  slide-class={n_slidecls}  (expect 15/15)")
print(f"phasing-AXIS leftover mentions (expect 0)={axis}")
print(f"d_within present={new_html.count('d_within')}  verdict-map present={'快速驗證每個甲基突破' in new_html}")
print(f"§ tags={len(re.findall(chr(167)+'[1-5] ·', new_html))}")
