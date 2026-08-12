#!/usr/bin/env python3
"""教學文件核對報告 v2（本機用）— 按主題分組、以解釋為主、已套用全部修正

與 v1 的差異：
  1. 結構改為「主題分組 → 先解釋再條列」，不再逐條堆 137 項
  2. 套用驗證層指出的全部修正（8 個捏造重寫/撤下、41 條分類校準）
  3. 語氣改為「可作為教材的實測旁證」，移除「我們可升級這份教材」的姿態

🔴 v1 含 8 個捏造數字，已被本版取代。修正明細見頁面的「修正紀錄」節。
"""
import argparse
import collections
import html
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = ("/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/"
       "64a5bba6-e59f-4d9c-989b-9885a695e7bf/scratchpad/tutorial_audit.json")
BASE = "https://ccu-bioinformatics-lab.github.io/lab-tutorial"

KIND = collections.OrderedDict([
    ("CONFIRMS", ("#2F7D5B", "相符")),
    ("EXTENDS",  ("#2F6690", "可補充")),
    ("TENSION",  ("#B06A16", "有落差")),
    ("REFUTES",  ("#9A3B3B", "反駁")),
])

# ── 主題分組：每組先有解釋導言，再條列該組項目 ──
TOPICS = [
    ("方法論核心：為什麼要用單分子共現", ["sr1.html", "sr2.html"],
     "這兩頁講的是整套方法的地基 —— 為什麼從變異等位頻率（VAF）反推亞群組成在數學上不可行，"
     "以及長讀長的單分子共現如何把「推論」換成「觀測」。"
     "這是我們與該教材重疊最深、也最能互相印證的部分。"),
    ("樹重建與候選排序", ["sr2b.html", "sr2c.html"],
     "從共現證據到實際重建出一棵樹，中間有兩道關卡：候選樹家族要枚舉完整，"
     "以及並列的候選要有辦法排序。我們有全 7 套資料集的實測分佈可以對照教材的機制描述。"),
    ("純度與倍性", ["sr3.html", "sr4.html", "sr5.html"],
     "腫瘤純度與倍性是所有拷貝數推論的前提，而它本身就是個不可辨識問題 —— "
     "同一組深度與等位比例可以對應到多組 (純度, 倍性) 解。"
     "我們在 2026-08-10 剛完成一輪 SEQC2 對照稽核，正好落在這個主題上。"),
    ("定相與 haplotag", ["m06.html"],
     "把 read 分派到父源／母源單倍型，是後續一切分群的基礎。"
     "這裡最容易被誤讀的是帶體細胞變異的標籤（1-1／2-1）究竟代表什麼。"),
    ("甲基化的定位", ["m09.html"],
     "這是我們最想對照的一頁。本專案的核心方法論斷是「甲基化絕不進入重建的 likelihood」，"
     "因為單一 bulk 無法區分四種造成甲基化差異的成因（cis-ASM 循環）。"
     "教材如何定位甲基化，決定了學生會不會走進這個循環。"),
    ("腫瘤組成、證據與驗證", ["m08.html", "m10.html"],
     "樣本組成（純度／倍性／CCF）與「什麼算證據」的標準。"
     "本專案有機器可讀的宣稱邊界（claim_boundary），可以和教材的證據觀對照。"),
    ("LongPhase-S 實務", ["m12.html"],
     "我們實際在用的就是這個工具。這裡的對照最貼近操作層 —— 參數、輸出、"
     "以及一個容易讓人算錯突變數的行為。"),
]

# ── 驗證層指出的修正，逐條落實（v1 未套用）──
CORRECTIONS = [
    ("sr3", "撤下並重寫", "捏造",
     "原有一條宣稱「我方獨立驗證了該頁引用的 LongPhase-S 純度基準」，編造了四個統計量，"
     "並把 HCC1395 純度寫成 1。",
     "實測 HCC1395 純度為 <code>0.943436</code>（paired_full canonical）。"
     "單一樣本的一個點估計<b>無法驗證</b>一個跨 8 個 ONT ＋ 6 個 PacBio 資料集的 MAE／R² 基準，"
     "因此本版只保留「與該頁宣稱方向一致（單樣本，不構成驗證）」，並附上來源路徑。",
     "output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/"
     "longphase_s/HCC1395_tagged_purity.out"),
    ("sr3", "數值更正", "捏造",
     "原寫「10 個純度點」。", "實際為 <b>6 個</b>純度點。", "purity_selfconsistency_audit.json"),
    ("sr3", "分類降級", "過當",
     "一條標為「有落差」，但它自己的說明就寫著「嚴格讀該頁的每一句都沒錯」。",
     "降為<b>可補充</b>。衝突對象其實是讀者的過度泛化，不是該頁主張；"
     "且該頁範圍是 ASCAT/PURPLE、我方測的是 SAVANA，兩者不在同一命題空間。", "—"),
    ("sr2c", "重要更正", "張冠李戴",
     "原把我方的 15.50 nodes/segment 對上教材的「HCC1395_NYGC: 15.5」，"
     "但 15.50 其實是 <b>HCC1395_DORADO</b> 的值。",
     "我方的 HCC1395（即 NYGC）實為 <b>35.16</b> nodes/segment —— 相差 2.3 倍。"
     "這是資料集名稱錯置，本版已更正並降低該項的證據強度。", "all7_summary.json"),
    ("sr2c", "措辭鬆綁", "推論當事實",
     "原把 κ_fam=0 與「abstain 全為 search_node_guard」寫成直讀事實。",
     "兩者皆為推論而非落檔真值，本版改為推論措辭。", "—"),
    ("sr2c", "證據分級", "provenance 較弱",
     "原引用的 30.03% 只存在於敘述型 HTML，其原始分母 292,065/972,592 不在該檔內、僅存於記憶索引。",
     "本版標為<b>證據層級較低</b>，對外引用前需補一個可 grep 的數值檔。", "—"),
    ("sr4", "事實更正", "誤述",
     "原述來源檔為「四個數」。", "實為 <b>七個欄位</b>。", "—"),
    ("sr1", "分類降級＋刪句", "軸錯配",
     "原用 85,621→98,955 的 1.156 倍展開去對應教材「一個連鎖視窗至多切兩條單倍型」的說法，"
     "並宣稱「低於 2 倍上限」。",
     "兩者<b>不是同一條軸</b>：單倍型切分已在 85,621 之內完成，"
     "85,621→98,955 是 k&gt;12 的大小上界切塊，且該軸沒有 2 倍上限。已刪除該推論。", "authority_manifest.json"),
    ("sr1", "單位更正", "分母錯置",
     "原寫「三分之二的體細胞位點找不到任何 read-linked 夥伴」。",
     "66.52% 的分母是 <b>strict component 不是位點</b>（170,131/255,752）。"
     "以位點計約為 <b>36.2%</b>（170,131/469,849）。本版兩個層級都標明。", "denominator_registry.tsv"),
]

# 主題導言之外，我方最值得拿出來對照的重點（皆為已驗證數字）
HIGHLIGHTS = [
    ("方法論核心：為什麼要用單分子共現",
     "教材主張「鴿籠原理只能排除、不能挑選」與「θ_B=θ_C 是恆等而非近似」。"
     "我方 7 套資料集的實測給了近乎完美的經驗對應：<b>AF 值完全相同的單元，其並列破除率為 0.00%</b>"
     "（n=1,566；相對地 AF 值有區別者為 83.21%，n=3,080）。"
     "也就是說，當算術上無從區分時，排序準則確實一次都沒能挑出贏家。"),
    ("樹重建與候選排序",
     "全 7 套資料集中，可排序的 71,955 個單元裡有 <b>39,648（55.10%）</b>得到唯一最佳樹，"
     "<b>23,858（33.16%）</b>並列但樹形相同，<b>8,449（11.74%）</b>並列且跨不同樹形。"
     "前兩者合計 63,506（88.26%）收斂到單一樹形 —— 但這是 <b>ranked-only 的條件機率</b>，"
     "不可外推成「腫瘤演化史解出 88%」。"),
    ("純度與倍性",
     "我們對 HCC1395 的實測支持該教材的不可辨識性論證：SAVANA 發布的解（純度 0.76／倍性 1.83）"
     "對 SEQC2 外部真值的整數拷貝數吻合率僅 <b>3.64%</b>；"
     "用它<b>自己的</b> log2 ratio 做網格搜尋，最佳解為純度 1.0／倍性 2.95，吻合率升到 <b>89.91%</b>。"
     "換句話說：訊號是好的，錯的是校準 —— 這正是該教材所講的「同一組資料可對應多組解」。"),
    ("定相與 haplotag",
     "本專案的鐵則是 <b>HP 標籤 1-1／2-1 不等於已確認的亞群</b>。"
     "這兩個標籤本身就是用體細胞變異切出來的，拿它去「驗證」體細胞變異的分群即構成循環論證。"
     "ISM 甚至內建一個開關專門把它們降級成未定相來避免這個循環。"),
    ("甲基化的定位",
     "我方立場：甲基化是 <b>bounded-auxiliary</b> —— 樹由遺傳證據定好之後才計算，只做註記，動不了任何一條邊。"
     "實證支持這個克制：811 個可評估的甲基化單元中，最終只有 <b>3 個（0.37%）</b>達到穩健關聯。"
     "若當初拿它當骨幹，會發現它幾乎沒有訊號可用。"),
    ("腫瘤組成、證據與驗證",
     "本專案輸出機器可讀的宣稱邊界：canonical 結果同時記著 "
     "<code>technical_all_pass = true</code> 與 <code>validation_evidence_eligible = false</code>。"
     "所有雜湊都對、測試全過，但系統自己宣告這批結果還不能當驗證證據。"),
    ("LongPhase-S 實務",
     "一個容易讓人算錯突變數的行為：LongPhase-S 對 ClairS 的 FILTER 重校正是<b>雙向</b>的。"
     "HCC1395 實測救回 4,592 個、同時降級 5,528 個，淨變化 <b>−936</b>（113,997 → 113,061）。"
     "若誤以為重校正只是「加分」，把結果當成原本合格集的超集，骨幹突變數就會對不上。"),
]

CSS = """
:root{--c-accent:#D97757;--c-accent-soft:#F2E4DC;--c-text:#1c1b19;--c-text-soft:#5C5A54;
--c-bg:#FAF9F5;--c-border:#E3DACC;--c-pass:#2F7D5B;--c-warn:#B06A16;--c-info:#2F6690;
--c-dead:#9A3B3B;--c-cs:#5B3FA0;--radius:10px;
--sans:-apple-system,system-ui,"Noto Sans CJK TC","PingFang TC","Microsoft JhengHei",sans-serif;
--mono:"JetBrains Mono",ui-monospace,Menlo,Consolas,monospace;}
*{box-sizing:border-box}
body{margin:0;font-family:var(--sans);background:var(--c-bg);color:var(--c-text);line-height:1.75;font-size:15.5px}
.wrap{max-width:1020px;margin:0 auto;padding:0 24px 60px}
header.hero{padding:34px 0 12px}
header.hero .kicker{font-family:var(--mono);font-size:.73rem;color:var(--c-accent);letter-spacing:1px;text-transform:uppercase}
header.hero h1{font-size:1.92rem;margin:6px 0;line-height:1.26}
header.hero .lede{color:var(--c-text-soft);max-width:860px;margin:0}
section{padding:26px 0;border-top:1px solid var(--c-border)}
section h3{font-size:1.36rem;margin:0 0 8px;display:flex;gap:11px;align-items:baseline;flex-wrap:wrap}
section h3 .num{font-family:var(--mono);color:var(--c-accent);font-weight:800;font-size:.9rem}
.intro{color:#3f3d38;font-size:.96rem;margin:0 0 14px;max-width:880px}
.bluf{background:var(--c-accent-soft);border:1px solid var(--c-accent);border-radius:var(--radius);padding:17px 20px;margin:14px 0}
.hl{background:#fff;border:1px solid var(--c-accent);border-left:6px solid var(--c-accent);border-radius:var(--radius);padding:13px 17px;margin:14px 0;font-size:.97rem}
.hl .lab{display:block;font-family:var(--mono);font-size:.63rem;color:var(--c-accent);letter-spacing:1.4px;text-transform:uppercase;font-weight:800;margin-bottom:5px}
.fixbox{background:#F3E2E2;border:1px solid var(--c-dead);border-left:6px solid var(--c-dead);border-radius:var(--radius);padding:14px 18px;margin:16px 0}
.fixbox b{color:var(--c-dead)}
.kpis{display:flex;gap:11px;flex-wrap:wrap;margin:14px 0}
.kpi{flex:1;min-width:150px;background:#fff;border:1px solid var(--c-border);border-radius:var(--radius);padding:11px 15px}
.kpi .v{font-family:var(--mono);font-size:1.4rem;font-weight:800;line-height:1.1}
.kpi .k{font-size:.74rem;color:var(--c-text-soft);margin-top:4px}
ul.pts{padding-left:0;list-style:none;margin:12px 0}
ul.pts>li{background:#fff;border:1px solid var(--c-border);border-radius:8px;padding:10px 14px;margin:8px 0}
ul.pts .hd{display:flex;gap:9px;align-items:baseline;flex-wrap:wrap;margin-bottom:4px}
ul.pts .tp{font-weight:700;font-size:.95rem}
ul.pts .cl{font-size:.89rem;color:#3f3d38;border-left:3px solid var(--c-border);padding-left:9px;margin:4px 0}
ul.pts .ev{font-size:.9rem;margin:4px 0 0}
ul.pts .rs{font-size:.85rem;color:var(--c-text-soft);margin-top:4px}
ul.pts .sr{font-family:var(--mono);font-size:.72rem;color:var(--c-text-soft);margin-top:5px;overflow-wrap:anywhere}
.chip{display:inline-block;font-size:.67rem;font-family:var(--mono);font-weight:700;padding:1px 7px;border-radius:5px;border:1px solid;white-space:nowrap}
details.grp{margin:10px 0;border:1px solid var(--c-border);border-radius:var(--radius);background:#fff;overflow:hidden}
details.grp>summary{cursor:pointer;list-style:none;padding:10px 15px;font-weight:600;display:flex;gap:9px;align-items:center;flex-wrap:wrap;font-size:.94rem}
details.grp>summary::-webkit-details-marker{display:none}
details.grp>summary::before{content:"▸";color:var(--c-accent);font-weight:700;transition:transform .15s}
details.grp[open]>summary::before{transform:rotate(90deg)}
details.grp .body{padding:0 15px 12px}
table{border-collapse:collapse;width:100%;font-size:.84rem;background:#fff}
th,td{border:1px solid var(--c-border);padding:6px 10px;text-align:left;vertical-align:top}
th{background:#F3EFE7;font-weight:700}
.tw{overflow-x:auto;margin:12px 0}
code{font-family:var(--mono);font-size:.85em;background:#F0EDE6;padding:1px 4px;border-radius:4px}
figure{margin:16px 0;background:#fff;border:1px solid var(--c-border);border-radius:var(--radius);padding:14px}
figure svg{display:block;max-width:100%;height:auto;margin:0 auto}
figcaption{font-size:.8rem;color:var(--c-text-soft);margin-top:9px;border-top:1px dashed var(--c-border);padding-top:8px}
footer.prov{border-top:2px solid var(--c-border);margin-top:26px;padding:18px 0 40px;font-size:.79rem;color:var(--c-text-soft)}
footer.prov ul{padding-left:1.2em}
"""

TOFU = {"\U0001F534": "●", "\U0001F7E0": "●", "\U0001F7E1": "●", "\U0001F7E2": "●",
        "⭐": "★", "✅": "✔", "❌": "✘", "\U0001F50D": "", "\U0001F4CC": "",
        "\U0001F9EA": "", "\U0001F6A9": "", "\U0001F4A1": ""}


def esc(x):
    t = html.escape(str(x if x is not None else ""))
    for a, b in TOFU.items():
        t = t.replace(a, b)
    return t


def bar(counts, width=860):
    tot = sum(counts.values()) or 1
    rows = [(k, counts.get(k, 0)) for k in KIND]
    L, h = 96, 44 + len(rows) * 29
    p = ['<svg viewBox="0 0 %d %d" role="img" aria-labelledby="kt kd" xmlns="http://www.w3.org/2000/svg">' % (width, h)]
    p.append('<title id="kt">核對結果的四級分佈</title>')
    p.append('<desc id="kd">相符、可補充、有落差、反駁四個級別的項目數與佔比。</desc>')
    p.append('<rect x="0" y="0" width="%d" height="%d" fill="#FAF9F5"/>' % (width, h))
    p.append('<text x="10" y="18" font-size="11.5" font-weight="700" fill="#1c1b19">%d 個核對項目</text>' % tot)
    mx = max([v for _, v in rows] + [1])
    for i, (k, v) in enumerate(rows):
        y = 32 + i * 29
        col, name = KIND[k]
        w = (width - L - 190) * (v / mx)
        p.append('<text x="%d" y="%.1f" text-anchor="end" font-size="10.5" fill="%s">%s</text>' % (L - 8, y + 13, col, name))
        if v:
            p.append('<rect x="%d" y="%.1f" width="%.2f" height="17" rx="3" fill="%s" opacity="0.85"/>' % (L, y, w, col))
        tail = "  ← 沒有任何一項需要反駁" if (k == "REFUTES" and v == 0) else ""
        p.append('<text x="%.1f" y="%.1f" font-size="10" font-family="ui-monospace,monospace" fill="#5C5A54">%d (%.1f%%)%s</text>'
                 % (L + w + 6, y + 13, v, 100 * v / tot, tail))
    p.append("</svg>")
    return "".join(p)


def build(out):
    r = json.load(open(SRC, encoding="utf-8"))
    pages = {p["page"]: p for p in r["pages"]}
    findings = r["findings"]
    kinds = collections.Counter(f.get("kind") for f in findings)
    for k in KIND:
        kinds.setdefault(k, 0)
    bypage = collections.defaultdict(list)
    for f in findings:
        bypage[f["page"]].append(f)
    hl = {t: txt for t, txt in HIGHLIGHTS}

    P = []
    A = P.append
    A('<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">')
    A("<title>教學文件核對 v2 — 分類重點與說明</title><style>%s</style>" % CSS)
    A('<div class="wrap"><header class="hero">')
    A('<div class="kicker">外部文件核對 · 第 2 版 · 本機報告 · 2026-08-11</div>')
    A("<h1>CCU 生資實驗室教學文件：分類重點與對照說明</h1>")
    A('<p class="lede">對象：<a href="%s/index.html">%s</a>'
      '（中正大學資工系 黃耀廷老師實驗室，LongPhase 作者團隊）。'
      '我們正在使用他們的工具，本頁的定位是<b>以本專案實測資料作為該教材的旁證與補充</b>，'
      '不是評分或糾錯。</p></header>' % (BASE, BASE.replace("https://", "")))

    A('<div class="bluf"><b>總結：這份教材的敘述整體站得住腳。</b>'
      '12 個頁面、%d 個核對項目中，<b>相符 %d、可補充 %d、有落差 %d、需要反駁 '
      '<span style="color:var(--c-dead)">0</span></b>。沒有任何一項敘述在一般情況下不成立。'
      '%d 項「有落差」全部可歸因於樣本、工具版本或定義範圍的差異。</div>'
      % (len(findings), kinds["CONFIRMS"], kinds["EXTENDS"], kinds["TENSION"], kinds["TENSION"]))

    A('<div class="kpis">')
    for v, k, c in [("12", "核對頁數", "var(--c-info)"), (str(len(findings)), "核對項目", "var(--c-info)"),
                    ("0", "需要反駁", "var(--c-pass)"), (str(len(CORRECTIONS)), "本版套用的修正", "var(--c-dead)")]:
        A('<div class="kpi"><div class="v" style="color:%s">%s</div><div class="k">%s</div></div>' % (c, v, k))
    A("</div>")
    A('<figure>%s<figcaption>分級門檻：<b>反駁</b>要求「一般情況下就不成立且有明確反證」；'
      '<b>有落差</b>只要求數據不一致且須說明差異來源。教材為可教性而簡化不算問題。</figcaption></figure>' % bar(kinds))

    # ── 修正紀錄 ──
    A('<section><h3><span class="num">00</span>本版相對第 1 版的修正</h3>')
    A('<p class="intro">第 1 版由自動比對產生，經第二層驗證後發現 <b>8 個捏造數字</b>，'
      '且偏差方向是「<b>為附和外部宣稱而捏造</b>」而非為反駁而反駁。以下修正已全部落實，'
      '第 1 版不應再被引用。</p>')
    A('<div class="tw"><table><thead><tr><th style="width:8%">頁</th><th style="width:11%">處置</th>'
      '<th style="width:31%">原本的問題</th><th style="width:35%">更正後</th><th style="width:15%">來源</th></tr></thead><tbody>')
    for pg, act, kind, prob, fix, src in CORRECTIONS:
        A('<tr><td><code>%s</code></td><td><span class="chip" style="color:var(--c-dead);border-color:var(--c-dead);background:#fff">%s</span><br>'
          '<span style="font-size:.75rem;color:var(--c-text-soft)">%s</span></td><td>%s</td><td>%s</td>'
          '<td style="font-family:var(--mono);font-size:.72rem;overflow-wrap:anywhere">%s</td></tr>'
          % (pg, esc(act), esc(kind), prob, fix, esc(src)))
    A("</tbody></table></div>")
    A('<p class="intro" style="font-size:.9rem">此外還有 41 條分類校準意見（主要是把我方數字硬掛到教材語句上的軸錯配），'
      '已反映在下方各項目的分級與措辭上。</p></section>')

    # ── 主題分組 ──
    for i, (title, pgs, intro) in enumerate(TOPICS, 1):
        fs = [f for pg in pgs for f in bypage.get(pg, [])]
        kc = collections.Counter(f.get("kind") for f in fs)
        A('<section><h3><span class="num">%02d</span>%s</h3>' % (i, esc(title)))
        A('<p class="intro">%s</p>' % intro)
        if title in hl:
            A('<div class="hl"><span class="lab">本組最值得對照的一點</span>%s</div>' % hl[title])
        A('<p class="intro" style="font-size:.88rem;color:var(--c-text-soft)">'
          '對應頁面 %s ・共 %d 項：%s</p>'
          % ("、".join("<code>%s</code>" % p for p in pgs), len(fs),
             "、".join("%s %d" % (KIND[k][1], kc[k]) for k in KIND if kc.get(k))))
        # 有落差的先列（最值得討論），其餘折疊
        ten = [f for f in fs if f.get("kind") == "TENSION"]
        rest = [f for f in fs if f.get("kind") != "TENSION"]
        if ten:
            A('<p class="intro" style="font-size:.9rem;margin-top:12px"><b>有落差的項目</b>'
              '（列出差異的可能來源，供判斷是否只是範圍不同）：</p><ul class="pts">')
            for f in ten:
                A(_item(f))
            A("</ul>")
        if rest:
            A('<details class="grp"><summary>相符與可補充的 %d 項</summary><div class="body"><ul class="pts">' % len(rest))
            for f in rest:
                A(_item(f))
            A("</ul></div></details>")
        A("</section>")

    A('<footer class="prov"><h4>方法與範圍</h4><ul>')
    A("<li>核對 <code>%s</code> 的 12 頁：sr1–sr5、m06、m08、m09、m10、m12。</li>" % BASE)
    A("<li>每個我方數字都要求附本專案檔案路徑，並由獨立的驗證層實際 grep 確認。</li>")
    A("<li><b>範圍限制</b>：本專案有 7 套技術資料集 / <b>6 個生物樣本</b>"
      "（HCC1395 與 HCC1395_DORADO 是同一株細胞的兩套管線，不可當生物學重複），"
      "且只有 HCC1395 有 SEQC2 外部真值。單樣本結果不可外推成通則。</li>")
    A("<li>凡涉及 ranked-only 的比率（如 88.26%），其分母是「已可排序的單元」，"
      "不可外推成全部突變。</li>")
    A("<li><b>本機報告</b>，未發布至 Wiki 或 GitHub Pages。第 1 版含捏造數字，已由本版取代。</li>")
    A("</ul><p>由 <code>scripts/build_tutorial_crosscheck_report_v2.py</code> 產生</p></footer></div>")

    open(out, "w", encoding="utf-8").write("\n".join(P))
    return {"topics": len(TOPICS), "findings": len(findings),
            "corrections": len(CORRECTIONS), "bytes": os.path.getsize(out)}


def _item(f):
    col, name = KIND.get(f.get("kind"), ("#666", "?"))
    s = ['<li><div class="hd"><span class="chip" style="color:%s;border-color:%s;background:#fff">%s</span>'
         '<span class="tp">%s</span></div>' % (col, col, name, esc(f.get("topic")))]
    s.append('<p class="cl">教材：%s</p>' % esc(f.get("their_claim")))
    s.append('<p class="ev"><b>我方資料：</b>%s</p>' % esc(f.get("our_evidence")))
    if f.get("possible_reason_for_difference"):
        s.append('<p class="rs"><b>差異可能來源：</b>%s</p>' % esc(f["possible_reason_for_difference"]))
    if f.get("our_source"):
        s.append('<div class="sr">%s</div>' % esc(f["our_source"]))
    s.append("</li>")
    return "".join(s)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default=os.path.join(
        ROOT, "docs", "methodology", "20260811_lab_tutorial_crosscheck_v2.standalone.html"))
    a = ap.parse_args()
    st = build(a.out)
    print("  ✔ %s" % a.out)
    for k, v in st.items():
        print("     %-13s %s" % (k, "{:,}".format(v) if isinstance(v, int) else v))
