#!/usr/bin/env python3
"""Illustrated case gallery: archive heatmap + this round's evidence plot, per locus.

Reads figure_cases.json (per-read methylation means for 9 representative loci,
extracted from the 2026-01 w5000 archive) and embeds the matching archive
cluster heatmaps as base64. No external requests, no emoji.
"""
import base64, io, json, math, os, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "20260726_methyl_case_gallery.data.json")
OUT = os.path.join(HERE, "20260726_methyl_case_gallery.standalone.html")
with open(DATA) as fh:
    CASES = json.load(fh)

C_ALT, C_REF, C_BAD, C_OK, C_MUT = "#0072B2", "#5b6672", "#B2182B", "#2E7D32", "#7b848f"

TAGS = {
    "M2": ("M2 僅 ALT 內多群", C_BAD,
           "ALT 群分裂成兩群且效應量大，同位點 REF 群沒有 —— 這是最像 subclone 的形狀。"),
    "M3": ("M3 區域性多群", "#8c6d1f",
           "ALT 與 REF 都分裂 —— 整個區域的分子彼此不同，與帶不帶突變無關。"),
    "M1": ("M1 ALT/REF 位移", "#4a6fa5",
           "兩群沒有分裂，但整體平均有位移。"),
    "M0": ("M0 無標記", C_OK,
           "既沒分裂也沒位移 —— 佔全體 87.3%，是常態。"),
}


def opt_split(v, minc=3):
    v = sorted(v); n = len(v)
    if n < 2 * minc: return None
    best = None
    for k in range(minc, n - minc + 1):
        lo, hi = v[:k], v[k:]
        w = (len(lo) * st.pstdev(lo) ** 2 + len(hi) * st.pstdev(hi) ** 2) / n
        if best is None or w < best[0]: best = (w, k)
    k = best[1]
    return (v[:k], v[k:], (v[k] + v[k - 1]) / 2.0)


def dotplot(case):
    """Two rows of per-read mean methylation: ALT above, REF below."""
    W, H = 700, 226
    L, R, T = 74, 22, 34
    pw = W - L - R
    sx = lambda v: L + v * pw
    p = ['<svg viewBox="0 0 {0} {1}" width="100%" role="img" aria-label="per-read methylation">'.format(W, H)]
    for i in range(6):
        v = i / 5.0
        p.append('<line x1="{0:.1f}" y1="{1}" x2="{0:.1f}" y2="{2}" stroke="currentColor" stroke-opacity=".10"/>'.format(sx(v), T - 6, H - 34))
        p.append('<text x="{0:.1f}" y="{1}" text-anchor="middle" font-size="10.5" fill="{2}">{3:.1f}</text>'.format(sx(v), H - 18, C_MUT, v))
    p.append('<text x="{0:.0f}" y="{1}" text-anchor="middle" font-size="11" fill="{2}">每條 read 在共同 CpG 集上的平均甲基化</text>'.format(L + pw / 2, H - 4, C_MUT))
    Q95 = 3.857
    rows = [("ALT", case["alt"], C_ALT, float(case["row"]["sep_alt"])),
            ("REF", case["ref"], C_REF, float(case["row"]["sep_ref"]))]
    for ri, (lab, vals, col, sepv) in enumerate(rows):
        y0 = T + ri * 84
        p.append('<text x="{0}" y="{1}" text-anchor="end" font-size="12.5" font-weight="700" fill="{2}">{3}</text>'.format(L - 12, y0 + 30, col, lab))
        p.append('<text x="{0}" y="{1}" text-anchor="end" font-size="10" fill="{2}">n={3}</text>'.format(L - 12, y0 + 44, C_MUT, len(vals)))
        p.append('<rect x="{0}" y="{1}" width="{2}" height="62" rx="4" fill="currentColor" fill-opacity=".03"/>'.format(L, y0, pw))
        sp = opt_split(vals)
        if sp:
            lo, hi, thr = sp
            gap = st.mean(hi) - st.mean(lo)
            if gap >= 0.20 and sepv >= Q95:
                p.append('<line x1="{0:.1f}" y1="{1}" x2="{0:.1f}" y2="{2}" stroke="{3}" stroke-width="2" stroke-dasharray="5 3"/>'.format(sx(thr), y0 + 2, y0 + 60, C_BAD))
                p.append('<text x="{0:.1f}" y="{1}" text-anchor="middle" font-size="10.5" font-weight="700" fill="{2}">Δβ={3:.2f}</text>'.format(sx(thr), y0 - 4, C_BAD, gap))
            for grp in (lo, hi):
                m = st.mean(grp)
                p.append('<line x1="{0:.1f}" y1="{1}" x2="{0:.1f}" y2="{2}" stroke="{3}" stroke-width="2.4" stroke-opacity=".55"/>'.format(sx(m), y0 + 6, y0 + 56, col))
        # deterministic vertical jitter, densest values spread most
        for i, v in enumerate(vals):
            yy = y0 + 12 + (i * 37 % 39)
            p.append('<circle cx="{0:.1f}" cy="{1}" r="4" fill="{2}" fill-opacity=".72" stroke="{2}" stroke-opacity=".95" stroke-width=".7"/>'.format(sx(v), yy, col))
    p.append("</svg>")
    return "".join(p)


def img64(path, width=880):
    from PIL import Image
    im = Image.open(path).convert("RGB")
    if im.width > width:
        im = im.resize((width, int(im.height * width / im.width)), Image.LANCZOS)
    buf = io.BytesIO(); im.save(buf, "PNG", optimize=True)
    return base64.b64encode(buf.getvalue()).decode(), im.size


def card(c):
    r = c["row"]
    tl, tc, tdesc = TAGS[c["tag"]]
    b64, size = img64(c["heatmap"])
    sp = opt_split(c["alt"])
    gap = (st.mean(sp[1]) - st.mean(sp[0])) if sp else 0
    spr = opt_split(c["ref"])
    gapr = (st.mean(spr[1]) - st.mean(spr[0])) if spr else 0
    read = READINGS[(c["chrom"], c["pos"])]
    return """
<section class="case">
<div class="chead">
  <h3>{chrom}:{pos} <span class="hp">germline HP{hp}</span></h3>
  <span class="chip" style="border-color:{tc};color:{tc}">{tl}</span>
  <span class="chip">{ls}</span>
  <span class="chip">CN {cn}</span>
  <span class="chip">局部密度 {dens}</span>
</div>
<div class="cgrid">
  <div>
    <div class="ctitle">本輪證據：每條 read 的平均甲基化</div>
    {dot}
    <table class="mini"><tbody>
      <tr><td>ALT 群分離度 S</td><td class="num">{sa}</td><td>門檻 null-95th = 3.86</td></tr>
      <tr><td>ALT 兩子群 Δβ</td><td class="num">{gap:+.3f}</td><td>{gapnote}</td></tr>
      <tr><td>REF 群分離度 S（對照）</td><td class="num">{sr}</td><td>子群 Δβ {gapr:+.3f}</td></tr>
      <tr><td>ALT−REF 平均位移</td><td class="num">{dar}</td><td>p={par} · 對 null 超額 {exc}</td></tr>
      <tr><td>共同 CpG / 全窗</td><td class="num">{nc} / {nt}</td><td>管線 optimal_k={k} · gating {gt}</td></tr>
    </tbody></table>
  </div>
  <div>
    <div class="ctitle">原始熱圖（archive · UPGMA）</div>
    <img src="data:image/png;base64,{b64}" alt="cluster heatmap {chrom}:{pos}" width="{w}" height="{h}">
    <div class="cap">Allele 欄深灰=ALT、淺灰=REF；HP 欄藍=HP1、橘=HP2-1、黃=HP2；白色=該 read 未覆蓋該 CpG。</div>
  </div>
</div>
<div class="reading"><b>判讀：</b>{read}</div>
</section>""".format(
        chrom=c["chrom"], pos=c["pos"], hp=c["hp"], tc=tc, tl=tl,
        ls=r["lineage_state"], cn=r["cn_state"], dens=r["local_density_50kb"],
        dot=dotplot(c), sa=r["sep_alt"], sr=r["sep_ref"], gap=gap, gapr=gapr,
        gapnote="達 M2 效應量門檻 0.20" if abs(gap) >= 0.2 else "未達 0.20 門檻",
        dar=r["d_altref"], par=r["p_altref"], exc=r["excess_vs_null"],
        nc=c["n_common"], nt=c["n_cpg"], k=r["pipeline_optimal_k"],
        gt="通過" if r["pipeline_passed_gating"] == "1" else "未過",
        b64=b64, w=size[0], h=size[1], read=read)


READINGS = {
    ("chr7", 36599703): "ALT 的 26 條 read 明顯分成低甲基與高甲基兩坨（Δβ=0.59），REF 的 16 條則集中。"
    "這是全基因組效應量最大的 M2 之一。但看右圖：主導 UPGMA 樹狀分群的是 <b>HP 欄的顏色分塊</b>"
    "（上半橘/黃=HP2 家族、下半藍=HP1），不是 Allele 欄 —— 熱圖上「看起來分兩群」主要來自 germline 單倍型，"
    "本輪的 within-HP 設計已把它扣掉，剩下的才是左圖。",
    ("chr3", 159027373): "ALT 27 條分裂（Δβ=0.48）而 REF 41 條沒有，形狀乾淨。"
    "但共同 CpG 只有 8 個（全窗 37 個）—— 平均項數少時，兩坨的分離會被放大，這是本例最大的不確定來源。",
    ("chr7", 80055524): "ALT 25 條分裂、REF 76 條集中。REF 深度是 ALT 的 3 倍，"
    "深度不對等本身就會讓深度低的一側看起來更分散，判讀時要扣掉這一層。",
    ("chr7", 106280304): "<b>這是關鍵對照。</b>ALT 與 REF <b>都</b>分裂成兩群 —— 整個區域的分子彼此不同，"
    "與帶不帶突變無關。這種形狀在全基因組佔 1.3%，而它證明「分裂」本身不足以推論譜系。",
    ("chr17", 69884800): "ALT 55 條、REF 41 條都分裂，且分裂位置相近。"
    "若只看 ALT 這一列會誤判成 subclone；把 REF 畫在同一張圖上，結論立刻反轉。",
    ("chr20", 29104598): "沒有分裂，但 ALT 整體比 REF 位移。位移的統計顯著性來自 read 數多，"
    "效應量本身不大 —— M1 全基因組佔 7.5%，且其骨幹富集在密度分層後消失。",
    ("chr4", 133868344): "位移方向清楚但兩側各只有十幾條 read。這種深度下單一位點無法定論，"
    "只能作為候選清單的一列。",
    ("chr7", 159316686): "ALT 138 條、REF 59 條，兩列都是單一連續分布 —— 這才是常態（87.3%）。"
    "深度這麼高仍看不到分裂，說明前面幾例的分裂不是深度不足造成的假象。",
    ("chr7", 153978134): "另一個高深度的 M0。把它與上面的 M2 並排看，"
    "可以校準「什麼程度的離散才算真的分裂」。",
}

cards = "".join(card(c) for c in CASES)
n_m2 = sum(1 for c in CASES if c["tag"] == "M2")

HTML = """<title>甲基案例圖解 — sSNV-lineage × 甲基標記層</title>
<style>
:root{{--bg:#fbfbfc;--fg:#1b1f24;--mut:#5b6672;--card:#fff;--line:#e3e6ea;--bad:#B2182B;--ok:#2E7D32;--a:#0072B2}}
@media (prefers-color-scheme:dark){{:root{{--bg:#14171a;--fg:#e8eaed;--mut:#9aa4af;--card:#1c2024;
--line:#2b3137;--bad:#f08a80;--ok:#7cc47f;--a:#5aa9dd}}}}
:root[data-theme=dark]{{--bg:#14171a;--fg:#e8eaed;--mut:#9aa4af;--card:#1c2024;--line:#2b3137;--bad:#f08a80;--ok:#7cc47f;--a:#5aa9dd}}
:root[data-theme=light]{{--bg:#fbfbfc;--fg:#1b1f24;--mut:#5b6672;--card:#fff;--line:#e3e6ea;--bad:#B2182B;--ok:#2E7D32;--a:#0072B2}}
body{{background:var(--bg);color:var(--fg);font:15px/1.7 -apple-system,"Noto Sans CJK TC","Microsoft JhengHei",sans-serif;margin:0;padding:32px 20px 80px}}
.wrap{{max-width:1080px;margin:0 auto}}
h1{{font-size:26px;margin:0 0 6px}}
h2{{font-size:19px;margin:44px 0 10px;padding-top:16px;border-top:1px solid var(--line)}}
h3{{font-size:16.5px;margin:0}}
.sub{{color:var(--mut);font-size:13px;margin:0 0 24px}}
.guide{{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:18px 22px;margin:20px 0}}
.guide h4{{margin:0 0 8px;font-size:14.5px}}
.guide ul{{margin:6px 0 14px;padding-left:20px}}
.guide li{{font-size:13.5px;margin-bottom:4px}}
.case{{background:var(--card);border:1px solid var(--line);border-radius:11px;padding:18px 20px;margin:20px 0}}
.chead{{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:14px}}
.hp{{font-size:12.5px;font-weight:500;color:var(--mut)}}
.chip{{font-size:11.5px;font-weight:600;padding:2px 9px;border-radius:11px;border:1px solid var(--line);color:var(--mut)}}
.cgrid{{display:grid;grid-template-columns:1fr 1fr;gap:20px}}
@media (max-width:900px){{.cgrid{{grid-template-columns:1fr}}}}
.ctitle{{font-size:12px;font-weight:700;letter-spacing:.04em;color:var(--mut);margin-bottom:6px}}
.cgrid img{{max-width:100%;height:auto;border:1px solid var(--line);border-radius:6px;background:#fff}}
.cap{{font-size:11.5px;color:var(--mut);margin-top:6px;line-height:1.55}}
table.mini{{width:100%;border-collapse:collapse;font-size:12.5px;margin-top:10px}}
table.mini td{{padding:5px 8px;border-bottom:1px solid var(--line)}}
table.mini td:last-child{{color:var(--mut);font-size:11.5px}}
.num{{text-align:right;font-variant-numeric:tabular-nums;font-weight:600}}
.reading{{margin-top:14px;font-size:13.5px;line-height:1.68;border-left:3px solid var(--a);padding-left:13px;color:var(--fg)}}
.warn{{background:var(--card);border:1px solid var(--line);border-left:5px solid var(--bad);border-radius:10px;padding:16px 20px;margin:20px 0}}
.warn .vt{{font-size:11px;font-weight:800;letter-spacing:.1em;color:var(--bad);margin-bottom:6px}}
code{{background:rgba(127,127,127,.13);padding:1px 5px;border-radius:3px;font-size:12.5px}}
</style>
<div class="wrap">
<h1>甲基案例圖解 —— 每個位點兩張圖並排</h1>
<p class="sub">HCC1395 · 2026-01 w5000 archive · 9 個代表位點涵蓋 4 種標記狀態 ·
所有數值取自 <code>ssnv_lineage_x_methyl_annotation.tsv</code></p>

<div class="warn"><div class="vt">先看這個，否則會誤讀</div>
<p style="margin:0">右側原始熱圖裡最醒目的分塊，多數來自 <b>HP 欄（germline 單倍型）</b>而不是 <b>Allele 欄</b>。
germline 甲基單倍型差異是真的，但它<b>不是 subclone</b>。左側是本輪的圖：已固定在同一個 germline HP 家族內，
只比 ALT vs REF —— 兩張圖看到的東西不一樣，這正是重點。</p></div>

<div class="guide">
<h4>怎麼讀左圖（本輪證據）</h4>
<ul>
<li>每個點 = 一條 tumor read，橫軸 = 它在<b>共同 CpG 集</b>上的平均甲基化（0 到 1）。
用共同 CpG 是為了避免「不同 read 覆蓋不同 CpG」偽造出分群。</li>
<li>上排藍點 = 帶突變的 ALT read；下排灰點 = 同一個 germline HP 內的 REF read，<b>它是內建對照</b>。</li>
<li>紅色虛線 = 最佳 2-切分位置，<b>只在該列同時滿足分離度 ≥ null-95th（3.86）與兩子群 Δβ ≥ 0.20 時才畫</b>
（＝標記層的旗標條件）。實心短豎線 = 各子群平均，無論有沒有畫虛線都會標。</li>
<li><b>判斷準則：只有「上排分裂、下排不分裂」才值得注意。</b>兩排都分裂 = 區域屬性。</li>
</ul>
<h4>怎麼讀右圖（archive 原始熱圖）</h4>
<ul>
<li>矩陣：行 = read（依 UPGMA 樹狀圖排序）、列 = CpG（依基因組位置）；
藍 = 低甲基、紅 = 高甲基、<b>白 = 該 read 沒覆蓋到該 CpG</b>。</li>
<li>四條註記欄由左到右：<b>HP</b>（藍=1、橘=2-1、黃=2、灰=0）、<b>Strand</b>、<b>Source</b>（紅=Tumor）、
<b>Allele</b>（深灰=ALT、淺灰=REF）。</li>
<li>大量白色 = read 覆蓋範圍差異很大，這是原始熱圖最容易誤導的地方。</li>
</ul>
<h4>四種標記狀態</h4>
<ul>
{taglist}
</ul>
</div>

<h2>案例（{n} 個）</h2>
{cards}

<h2>把 9 個案例合起來看</h2>
<p>M2（{n_m2} 例在此）在全基因組佔 <b>3.8%</b>，而它在 sSNV 骨幹位點與孤立位點上出現的頻率相同
（3.7% vs 3.6%，OR=1.03, p=0.854）。M3 這種「ALT 與 REF 都分裂」的形狀證明分裂本身不專屬於突變等位；
把 ALT 群的分離度與同位點 REF 群相比，全基因組是 <b>23.9% vs 23.6%</b>（配對 p=0.29）。</p>
<p>所以這些圖的正確用途是：<b>挑案例、看資料長相、寫 limitation</b>；
而不是拿其中任何一張宣稱該位點屬於某個 subclone。完整統計與標記層欄位定義見
<code>20260726_ssnv_lineage_x_methyl.standalone.html</code>。</p>
</div>
""".format(
    cards=cards, n=len(CASES), n_m2=n_m2,
    taglist="".join('<li><b style="color:{0}">{1}</b> —— {2}</li>'.format(c, t, d)
                    for t, c, d in TAGS.values()))

with open(OUT, "w") as fh:
    fh.write(HTML)
print("wrote", OUT, os.path.getsize(OUT), "bytes")
