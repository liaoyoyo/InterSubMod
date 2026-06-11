#!/usr/bin/env python3
"""Build report_v4.html: re-architect report_v3's verified prose into 5 sections.

Mirrors presenter_v4. User decisions: NEW v4 (keep v3) / phasing axis removed / §4 productive-failure.
Reuses report_v3 prose+SVG losslessly. Authors 2 NEW sections. No number production.
"""
import re

D = "/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/in_progress/20260603_methyl_haplotype_story_deck_v2"
SRC = f"{D}/report_v3.html"
OUT = f"{D}/report_v4.html"

html = open(SRC, encoding="utf-8").read()

FOOT = '<div class="foot" id="prov">'
h2_first = re.search(r'<h2 id=', html).start()
foot_at = html.index(FOOT)
head = html[:h2_first]
middle = html[h2_first:foot_at]
foot_tail = html[foot_at:]

sections = [s for s in re.split(r'(?=<h2 id=)', middle) if s.strip()]
def find(idstr):
    hits = [s for s in sections if f'<h2 id="{idstr}"' in s]
    assert len(hits) == 1, f'id {idstr} matched {len(hits)}'
    return hits[0]

def relabel(sec, new_h2):
    return re.sub(r'<h2\b.*?</h2>', new_h2, sec, count=1, flags=re.S)

# ---- strip phasing RESEARCH AXIS ----
def strip_c7(s):
    s = re.sub(r'\s*<tr><td><b>純 phasing / germline-het（非甲基）延長 somatic read</b></td>.*?</tr>', '', s, flags=re.S)
    s = s.replace('；但用<b>純 phasing/germline-het（非甲基）</b>延長 read 是難度價值完全不同的 LIVE 主軸', '')
    return s
def strip_c8(s):
    return s.replace(' + 純 phasing read extension', '')

s0   = find('s0')
c1   = relabel(find('c1'),  '<h2 id="s1"><span class="n">§1</span>個案存在：突變 × haplotype × 甲基差異 × 生物機制 共存</h2>')
c2   = relabel(find('c2'),  '<h2 id="s2"><span class="n">§2</span>如何找到並驗證這個位點（四步流程 + 工具地景）</h2>')
appx = relabel(find('appx'),"<h2 id=\"appx\"><span class=\"n\">§2</span>方法圖示：Δβ 算法 vs 舊 ISM（Cramér's V）</h2>")
c6   = relabel(find('c6'),  '<h2 id="s2c"><span class="n">§2</span>ISM 作為甲基分析工具的價值（能做/不能做）</h2>')
c4   = relabel(find('c4'),  '<h2 id="s3"><span class="n">§3</span>困難與文獻：為什麼甲基研究這麼難</h2>')
c3   = relabel(find('c3'),  '<h2 id="s4b"><span class="n">§4</span>推廣：其他位點稀有、其他樣本 private</h2>')
c7   = relabel(strip_c7(find('c7')), '<h2 id="s4c"><span class="n">§4</span>其他甲基方向（productive-failure 盤點）</h2>')
c8   = relabel(strip_c8(find('c8')), '<h2 id="s4d"><span class="n">§4</span>前路 + verdict</h2>')
refs = relabel(find('refs'),'<h2 id="s5"><span class="n">§5</span>參考文獻（已驗證）</h2>')
# c5 (phasing panorama) dropped
assert 'LIVE 主軸' not in c7 and '純 phasing' not in c7, 'c7 axis not stripped'
assert '純 phasing read extension' not in c8, 'c8 axis not stripped'

new_whynot = """<h2 id="s3b"><span class="n">§3</span>為何之前偵測不到：舊方法看得到訊號，卻不是 somatic 特異</h2>
<p>「之前為什麼偵測不到 somatic ASM」的精確答案是<b>方法學，不是缺資料</b>。舊 ISM 用 region 內無監督甲基分群 + <b>Cramér's V</b> 量「分群是否關聯單倍型」——這個關聯<b>確實存在</b>（V=0.80, p=0.0, 通過 gating），但它<b>不是 somatic 特異</b>：同時混了 germline-allelic（HP1↔HP2 之間常見的 ASM 也高 V）、copy-number（LOH/CN 造成假甲基差）、以及 n_reads 覆蓋（O11 epipolymorphism 的陷阱：異質性 AUC 看似 0.845，做 n_reads 校正後掉到 0.530，是 artifact）。再加上舊輸出沒有 normal 基線、沒有方向性，BRCA2 這種真實位點就淹沒在數千個高-V 位點裡。</p>
<table>
  <tr><th>混淆源</th><th>為何看不到 somatic</th><th>新方法如何逐源隔離</th></tr>
  <tr><td>germline-allelic</td><td>HP1↔HP2 常見 ASM 也高 V</td><td>HP1 vs HP1-1（同單倍型只差 somatic allele）</td></tr>
  <tr><td>copy-number</td><td>LOH/CN 造成假甲基差</td><td>copy-partition（純 allele d_within=−0.023）</td></tr>
  <tr><td>n_reads（覆蓋）</td><td>O11 異質性 AUC 0.845→0.530</td><td>coverage-rarefied</td></tr>
  <tr><td>無 normal 基線</td><td>分不出 cis vs drift</td><td>normal-anchored cis-test</td></tr>
  <tr><td>類別關聯（非方向）</td><td>無法排序 / 定方向</td><td>方向性 Δβ（可全基因組排序）</td></tr>
</table>
<div class="caveat">⚠ epipolymorphism 當「<b>困難真實</b>」可引（Landan 2012 / Derrien 2021）；當「<b>判別特徵</b>」是 n_reads artifact（O11，校正後 AUC 0.530），兩者勿混。一句話：不是缺資料（tag/normal 模組都在 binary 內），是<b>方法學（指標設計）</b> → 任務 #19。</div>
"""

new_verdict = """<h2 id="s4"><span class="n">§4</span>快速驗證每個甲基突破想法 → 哪些是 possible 甲基研究</h2>
<p>本報告把每一個「想用甲基去突破」的問題都做了<b>快速驗證</b>，而不是留模糊待辦。這正是 productive-failure 的價值：<b>把每條死路都驗證清楚、把 reviewer 會質疑的點先自己關掉</b>。結果是六個問題＝一個稀有但真實的 anchor（Q1–Q3 定位）＋五條關閉的路。</p>
<table>
  <tr><th>想用甲基突破的問題</th><th class="c">快速驗證</th><th>依據</th></tr>
  <tr><td>Q1 哪些 somatic 位點有甲基差異</td><td class="c"><span class="rel-p">🟡 PARTIAL</span></td><td>BRCA2＋regime-A 15/23；strong-ASM 0.34%（稀有）</td></tr>
  <tr><td>Q2 甲基差異與哪種 tag 關聯</td><td class="c"><span class="rel-p">🟡 PARTIAL</span></td><td>HP-axis modest；copy-partition 顯示多為 copy</td></tr>
  <tr><td>Q3 定位明顯差異位點</td><td class="c"><span class="rel-p">🟡 定位</span></td><td>regime-A 可輸出（稀有）</td></tr>
  <tr><td>Q3b 用甲基輔助 phase / tag / HP3 / unphase</td><td class="c"><span class="rel-c">🔴 NEGATIVE</span></td><td>甲基判別力全域 AUC=0.505（≈隨機）</td></tr>
  <tr><td>Q4 分析 subclone（多分群 / 演化樹）</td><td class="c"><span class="rel-c">🔴 NEGATIVE</span></td><td>不過 germline-het baseline</td></tr>
  <tr><td>Q5 二次打擊 + 發生順序</td><td class="c"><span class="rel-c">🔴 NEGATIVE</span></td><td>region-level 無 two-hit 訊號</td></tr>
  <tr><td>Q6 甲基提升 F1（降 FP）</td><td class="c"><span class="rel-c">⚫ DEAD</span></td><td>LOSO 100% circular（held-out −0.00012），已 concluded</td></tr>
</table>
<div class="caveat">誠實口徑：只有 Q1–Q3（定位）是 PARTIAL（稀有真實）；Q3b/Q4/Q5/Q6 全 NEGATIVE/DEAD。<b>possible 只剩 characterization ＋ tooling</b>（非 filter / subclone / two-hit）；filter（Q6）已 concluded，不列待辦。</div>
"""

order = [s0, c1, c2, appx, c6, c4, new_whynot, new_verdict, c3, c7, c8, refs]
new_middle = "\n".join(order)

# ---- rebuild TOC ----
new_toc = """<nav class="toc">
  <h4>目錄 · v4</h4>
  <a href="#s0">§0 一句話主張</a>
  <a href="#s1">§1 個案存在 (BRCA2)</a>
  <a href="#s2">§2 如何找到 + 方法差異</a>
  <a href="#s3">§3 為何困難</a>
  <a href="#s4">§4 快速驗證 → 哪些 possible</a>
  <a href="#s5">§5 參考文獻</a>
  <a href="#prov">Provenance</a>
</nav>"""
head = re.sub(r'<nav class="toc">.*?</nav>', lambda m: new_toc, head, count=1, flags=re.S)

# ---- hero + head-comment updates (drop phasing axis) ----
head = head.replace(
    "8 章報告 — 從一個真實個案出發，說明方法、推廣、困難、全景、工具價值與其他方向",
    "5 段報告 — 從一個真實個案出發：①個案存在 ②如何找到＋方法差異 ③為何困難 ④快速驗證→哪些 possible ⑤數據驗證")
head = head.replace("characterization · v3", "characterization · v4")
head = head.replace("provenance_note: v3 主軸＝ISM 甲基工具；phasing 整合於 C5。",
                    "provenance_note: v4 主軸＝ISM 甲基-haplotype ASM 5 段；phasing 軸不放。")

# ---- prov footer: strip phasing axis ----
foot_tail = foot_tail.replace(
    "v3 主軸＝ISM 甲基-read-群聚工具（C1 BRCA2 hero → C6 工具定義）；phasing 整合於 C5 全景（強正能力）。",
    "v4 主軸＝ISM 甲基-haplotype ASM（純甲基焦點，phasing 軸不在本 deck）；5 段＝個案/方法/困難/快速驗證/數據。")

new_html = head + new_middle + "\n" + foot_tail
open(OUT, "w", encoding="utf-8").write(new_html)

# ---- self-check ----
axis = sum(new_html.count(x) for x in ["LIVE 主軸", "phasing 主軸", "phasing 整合於 C5", "純 phasing read extension", "C5 全景"])
print(f"OUT={OUT}")
print(f"h2 count={len(re.findall(r'<h2 ', new_html))} (s0,§1,§2x3,§3x2,§4x4,§5 = 11)")
print(f"phasing-AXIS leftover (expect 0)={axis}")
print(f"verdict-map={'快速驗證每個甲基突破' in new_html}  why-not={'為何之前偵測不到' in new_html}")
appx_kept = 'id="appx"' in new_html
print(f"#appx anchor kept={appx_kept}  svg count={new_html.count('<svg')}")
print(f"d_within={new_html.count('d_within')}")
