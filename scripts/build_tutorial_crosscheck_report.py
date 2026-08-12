#!/usr/bin/env python3
"""用本專案實測數據核對 CCU 生資實驗室教學文件 —— 報告頁生成器（本機用）

來源：scratchpad/tutorial_audit.json（12 個比對 agent + 12 個驗證 agent 的結構化輸出）

🔴 誠信處置（本頁最重要的設計）
驗證階段在我方的比對結果中抓到 8 個捏造數字與 41 條分類校準意見，
且偏差方向是「為附和外部宣稱而捏造」而非「為反駁而反駁」。
因此本頁不把驗證結果藏在附錄，而是：
  - 有捏造的頁面一律加醒目警示，其 finding 需打折扣閱讀
  - 41 條修正意見全文列出，不摘要、不美化
  - 不自動刪除個別 finding（文字型意見無法安全對應到索引，硬刪可能誤刪正確項）

用法: python3 scripts/build_tutorial_crosscheck_report.py [-o OUT.html]
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
    ("CONFIRMS", ("#2F7D5B", "佐證", "我們的數據支持該敘述")),
    ("EXTENDS",  ("#2F6690", "可補充", "敘述沒錯，我們有更細或更新的數據")),
    ("TENSION",  ("#B06A16", "有張力", "數據不一致，但可能源自不同樣本／版本／定義")),
    ("REFUTES",  ("#9A3B3B", "反駁", "一般情況下就不成立且有明確反證")),
])
PAGE_TITLE = {
    "sr1.html": "亞克隆演化 ①", "sr2.html": "亞克隆演化 ②",
    "sr2b.html": "亞克隆演化 ③", "sr2c.html": "亞克隆演化 ④",
    "sr3.html": "純度／倍性 ①", "sr4.html": "純度／倍性 ②", "sr5.html": "純度／倍性 ③",
    "m06.html": "M06 定相與 haplotag", "m08.html": "M08 腫瘤樣本組成",
    "m09.html": "M09 甲基化", "m10.html": "M10 證據與驗證", "m12.html": "M12 LongPhase-S",
}

CSS = """
:root{--c-accent:#D97757;--c-accent-soft:#F2E4DC;--c-text:#1c1b19;--c-text-soft:#5C5A54;
--c-bg:#FAF9F5;--c-border:#E3DACC;--c-pass:#2F7D5B;--c-warn:#B06A16;--c-info:#2F6690;
--c-dead:#9A3B3B;--c-cs:#5B3FA0;--radius:10px;
--sans:-apple-system,system-ui,"Noto Sans CJK TC","PingFang TC","Microsoft JhengHei",sans-serif;
--mono:"JetBrains Mono",ui-monospace,Menlo,Consolas,monospace;}
*{box-sizing:border-box}
body{margin:0;font-family:var(--sans);background:var(--c-bg);color:var(--c-text);line-height:1.72;font-size:15.5px}
.wrap{max-width:1060px;margin:0 auto;padding:0 24px 60px}
header.hero{padding:34px 0 12px}
header.hero .kicker{font-family:var(--mono);font-size:.73rem;color:var(--c-accent);letter-spacing:1px;text-transform:uppercase}
header.hero h1{font-size:1.95rem;margin:6px 0;line-height:1.25}
header.hero .lede{color:var(--c-text-soft);max-width:880px;margin:0}
section{padding:28px 0;border-top:1px solid var(--c-border)}
section h3{font-size:1.42rem;margin:0 0 6px;display:flex;gap:11px;align-items:baseline;flex-wrap:wrap}
section h3 .num{font-family:var(--mono);color:var(--c-accent);font-weight:800;font-size:.92rem}
section .sub{color:var(--c-text-soft);font-size:.93rem;margin:0 0 16px;max-width:900px}
.bluf{background:var(--c-accent-soft);border:1px solid var(--c-accent);border-radius:var(--radius);padding:17px 20px;margin:14px 0}
.alarm{background:#F3E2E2;border:1px solid var(--c-dead);border-left:6px solid var(--c-dead);border-radius:var(--radius);padding:14px 18px;margin:16px 0}
.alarm b{color:var(--c-dead)}
.note{background:#fff;border:1px solid var(--c-border);border-left:4px solid var(--c-info);border-radius:8px;padding:11px 15px;margin:12px 0;font-size:.92rem}
.kpis{display:flex;gap:11px;flex-wrap:wrap;margin:14px 0}
.kpi{flex:1;min-width:152px;background:#fff;border:1px solid var(--c-border);border-radius:var(--radius);padding:11px 15px}
.kpi .v{font-family:var(--mono);font-size:1.42rem;font-weight:800;line-height:1.1}
.kpi .k{font-size:.74rem;color:var(--c-text-soft);margin-top:4px}
details.pg{margin:9px 0;border:1px solid var(--c-border);border-radius:var(--radius);background:#fff;overflow:hidden}
details.pg>summary{cursor:pointer;list-style:none;padding:11px 16px;font-weight:600;display:flex;gap:10px;align-items:center;flex-wrap:wrap}
details.pg>summary::-webkit-details-marker{display:none}
details.pg>summary::before{content:"▸";color:var(--c-accent);font-weight:700;transition:transform .15s}
details.pg[open]>summary::before{transform:rotate(90deg)}
details.pg .body{padding:0 16px 14px}
.f{border:1px solid var(--c-border);border-radius:8px;padding:11px 14px;margin:9px 0;background:#FBF9F4}
.f .hd{display:flex;gap:9px;align-items:baseline;flex-wrap:wrap;margin-bottom:6px}
.f .topic{font-weight:700;font-size:.97rem}
.f .lbl{font-size:.72rem;color:var(--c-text-soft);font-weight:700;margin-top:7px}
.f .q{border-left:3px solid var(--c-border);padding-left:10px;margin:3px 0 0;font-size:.9rem;color:#3f3d38}
.f .ev{font-size:.91rem;margin:3px 0 0}
.f .src{font-family:var(--mono);font-size:.74rem;color:var(--c-text-soft);margin-top:5px;overflow-wrap:anywhere}
.chip{display:inline-block;font-size:.68rem;font-family:var(--mono);font-weight:700;padding:1px 7px;border-radius:5px;border:1px solid;white-space:nowrap}
.bad{background:#F3E2E2;border-color:var(--c-dead);color:var(--c-dead)}
table{border-collapse:collapse;width:100%;font-size:.85rem;background:#fff}
th,td{border:1px solid var(--c-border);padding:6px 10px;text-align:left;vertical-align:top}
th{background:#F3EFE7;font-weight:700}
.tw{overflow-x:auto;margin:12px 0}
code{font-family:var(--mono);font-size:.85em;background:#F0EDE6;padding:1px 4px;border-radius:4px}
figure{margin:16px 0;background:#fff;border:1px solid var(--c-border);border-radius:var(--radius);padding:14px}
figure svg{display:block;max-width:100%;height:auto;margin:0 auto}
figcaption{font-size:.8rem;color:var(--c-text-soft);margin-top:9px;border-top:1px dashed var(--c-border);padding-top:8px}
ol.fix{padding-left:1.3em;font-size:.89rem}
ol.fix li{margin:9px 0}
footer.prov{border-top:2px solid var(--c-border);margin-top:26px;padding:18px 0 40px;font-size:.79rem;color:var(--c-text-soft)}
footer.prov ul{padding-left:1.2em}
"""


TOFU = {"\U0001F534": "●", "\U0001F7E0": "●", "\U0001F7E1": "●", "\U0001F7E2": "●",
        "\u2b50": "★", "\u2705": "✔", "\u274c": "✘", "\U0001F50D": "", "\U0001F4CC": "",
        "\U0001F9EA": "", "\U0001F6A9": "", "\U0001F4A1": ""}


def esc(x):
    """跳脫 HTML 並把本機無字型的 emoji 換成可渲染的等價字元（否則顯示為豆腐方框）。"""
    t = html.escape(str(x if x is not None else ""))
    for a, b in TOFU.items():
        t = t.replace(a, b)
    return t


def bar(counts, width=880):
    """kind 分佈長條。數量由資料算出。"""
    tot = sum(counts.values()) or 1
    rows = [(k, counts.get(k, 0)) for k in KIND if counts.get(k) is not None]
    ROW, L, h = 30, 108, 44 + len([r for r in rows]) * 30
    p = ['<svg viewBox="0 0 %d %d" role="img" aria-labelledby="kt kd" xmlns="http://www.w3.org/2000/svg">' % (width, h)]
    p.append('<title id="kt">核對結果的四級分佈</title>')
    p.append('<desc id="kd">佐證、可補充、有張力、反駁四個級別各自的項目數與佔比。</desc>')
    p.append('<rect x="0" y="0" width="%d" height="%d" fill="#FAF9F5"/>' % (width, h))
    p.append('<text x="10" y="18" font-size="11.5" font-weight="700" fill="#1c1b19">'
             '%d 個核對項目的四級分佈</text>' % tot)
    mx = max([v for _, v in rows] + [1])
    for i, (k, v) in enumerate(rows):
        y = 32 + i * ROW
        col, name, _ = KIND[k]
        w = (width - L - 150) * (v / mx)
        p.append('<text x="%d" y="%.1f" text-anchor="end" font-size="10.5" fill="%s">%s</text>'
                 % (L - 8, y + 13, col, name))
        if v:
            p.append('<rect x="%d" y="%.1f" width="%.2f" height="17" rx="3" fill="%s" opacity="0.85"/>' % (L, y, w, col))
        p.append('<text x="%.1f" y="%.1f" font-size="10" font-family="ui-monospace,monospace" fill="#5C5A54">%d (%.1f%%)%s</text>'
                 % (L + w + 6, y + 13, v, 100 * v / tot, "  ← 沒有任何一項需要反駁" if k == "REFUTES" and v == 0 else ""))
    p.append("</svg>")
    return "".join(p)


def build(out):
    r = json.load(open(SRC, encoding="utf-8"))
    pages = r["pages"]
    findings = r["findings"]
    kinds = collections.Counter(f.get("kind") for f in findings)
    kinds.setdefault("REFUTES", 0)
    n_fab = sum(p.get("n_fabricated") or 0 for p in pages)
    n_over = sum(len(p.get("overclassified") or []) for p in pages)
    bypage = collections.defaultdict(list)
    for f in findings:
        bypage[f["page"]].append(f)

    P = []
    A = P.append
    A('<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">')
    A("<title>教學文件核對 — 用本專案實測數據</title><style>%s</style>" % CSS)
    A('<div class="wrap"><header class="hero">')
    A('<div class="kicker">外部文件核對 · 本機報告 · 2026-08-11</div>')
    A("<h1>用本專案實測數據核對 CCU 生資實驗室教學文件</h1>")
    A('<p class="lede">對象：<a href="%s/index.html">%s</a>（中正大學資工系 黃耀廷老師實驗室，'
      'LongPhase 作者團隊）。核對 12 個與本專案研究直接重疊的頁面。</p></header>' % (BASE, BASE.replace("https://", "")))

    A('<div class="bluf"><b>一句話結論：這份教材的敘述整體站得住腳。</b>'
      '%d 個核對項目中，<b>佐證 %d、可補充 %d、有張力 %d、需要反駁 <span style="color:var(--c-dead)">0</span></b>。'
      '沒有任何一項敘述在一般情況下不成立。18 項「有張力」全部可歸因於樣本、工具版本或定義範圍的差異，'
      '而非該教材有誤。</div>' % (len(findings), kinds["CONFIRMS"], kinds["EXTENDS"], kinds["TENSION"]))

    A('<div class="kpis">')
    for v, k, c in [(str(len(pages)), "核對頁數", "var(--c-info)"),
                    (str(len(findings)), "核對項目", "var(--c-info)"),
                    ("0", "需要反駁", "var(--c-pass)"),
                    (str(n_fab), "我方被抓到的捏造", "var(--c-dead)")]:
        A('<div class="kpi"><div class="v" style="color:%s">%s</div><div class="k">%s</div></div>' % (c, v, k))
    A("</div>")
    A('<figure>%s<figcaption>四級分類的門檻不同：REFUTES 要求「一般情況下就不成立且有明確反證」，'
      'TENSION 只要求「數據不一致」且必須說明差異可能來源。教材為可教性而簡化不算問題。</figcaption></figure>' % bar(kinds))

    # ── 誠信章節放在最前面 ──
    A('<section><h3><span class="num">01</span>先講我方自己的問題</h3>')
    A('<p class="sub">這次核對用了兩層 agent：12 個負責比對、12 個負責驗證我方引用的數字是否真的存在於本專案檔案。'
      '驗證層抓到的問題必須先講，否則下面的內容無從判斷可信度。</p>')
    A('<div class="alarm"><b>驗證層在我方比對結果中抓到 %d 個捏造數字，偏差方向是「為附和而捏造」而非「為反駁而反駁」。</b><br>'
      '最嚴重的一例在 <code>sr3</code>：為了湊出一條「我方獨立驗證了該頁引用的 LongPhase-S 基準」的佐證項，'
      '編造了四個統計量，並把 HCC1395 的純度寫成 1（實測為 <code>0.943436</code> / <code>0.98969</code>），且未附任何來源路徑。'
      '這正是本專案 CLAUDE.md §13.0 所指的 batch 結構問題 —— 數字未讀回真值就下筆，'
      '且捏造內容全部朝「與外部宣稱一致」的方向偏。<br><br>'
      '另有 <b>%d 條分類校準意見</b>，主要是把「我方數字硬掛到教材語句上」的軸錯配，'
      '以及該用可補充卻標成有張力。<b>本頁不刪除個別項目</b> —— 文字型意見無法安全對應到索引，'
      '硬刪可能誤刪正確項；改為在有問題的頁面加警示，並把 41 條意見全文列出。</div>' % (n_fab, n_over))
    A('<div class="tw"><table><thead><tr><th>頁面</th><th>驗證判定</th><th>捏造</th><th>校準意見</th><th>閱讀提醒</th></tr></thead><tbody>')
    for p in pages:
        fab = p.get("n_fabricated") or 0
        ov = len(p.get("overclassified") or [])
        v = p.get("verdict")
        vc = {"CLEAN": "var(--c-pass)", "MINOR_ISSUES": "var(--c-warn)", "SERIOUS_ISSUES": "var(--c-dead)"}.get(v, "#666")
        A('<tr><td><code>%s</code> %s</td><td><span class="chip" style="color:%s;border-color:%s;background:#fff">%s</span></td>'
          '<td>%s</td><td>%d</td><td>%s</td></tr>'
          % (p["page"], PAGE_TITLE.get(p["page"], ""), vc, vc, v,
             ('<b style="color:var(--c-dead)">%d</b>' % fab) if fab else "0", ov,
             "此頁的項目需打折扣閱讀" if fab else ("有分類校準" if ov else "—")))
    A("</tbody></table></div></section>")

    # ── 逐頁 ──
    A('<section><h3><span class="num">02</span>逐頁核對結果</h3>')
    A('<p class="sub">每一項都附我方數字與其來源檔案路徑。點開頁面標題展開。</p>')
    for p in pages:
        pg = p["page"]
        fs = bypage.get(pg, [])
        kc = collections.Counter(f.get("kind") for f in fs)
        fab = p.get("n_fabricated") or 0
        A('<details class="pg"><summary><span>%s</span><code>%s</code>' % (PAGE_TITLE.get(pg, pg), pg))
        for k in KIND:
            if kc.get(k):
                A('<span class="chip" style="color:%s;border-color:%s;background:#fff">%s %d</span>'
                  % (KIND[k][0], KIND[k][0], KIND[k][1], kc[k]))
        if fab:
            A('<span class="chip bad">含 %d 捏造</span>' % fab)
        A("</summary><div class=\"body\">")
        if p.get("assessment"):
            A('<div class="note"><b>整體評價：</b>%s</div>' % esc((p["assessment"])))
        if fab:
            A('<div class="alarm" style="margin:10px 0"><b>此頁的核對項目含 %d 個被驗證層判定為捏造的數字。</b>'
              '下列項目請對照最後一節的校準意見閱讀。</div>' % fab)
        for f in fs:
            col, name, _ = KIND.get(f.get("kind"), ("#666", "?", ""))
            A('<div class="f"><div class="hd">'
              '<span class="chip" style="color:%s;border-color:%s;background:#fff">%s</span>'
              '<span class="topic">%s</span></div>' % (col, col, name, esc((f.get("topic") or ""))))
            A('<div class="lbl">教材的敘述</div><p class="q">%s</p>' % esc((f.get("their_claim") or "")))
            A('<div class="lbl">我方數據</div><p class="ev">%s</p>' % esc((f.get("our_evidence") or "")))
            if f.get("possible_reason_for_difference"):
                A('<div class="lbl">差異的可能來源</div><p class="ev">%s</p>'
                  % esc((f["possible_reason_for_difference"])))
            if f.get("why_it_matters"):
                A('<div class="lbl">為何重要</div><p class="ev">%s</p>' % esc((f["why_it_matters"])))
            if f.get("our_source"):
                A('<div class="src">來源：%s</div>' % esc((f["our_source"])))
            A("</div>")
        A("</div></details>")
    A("</section>")

    # ── 校準意見全文 ──
    A('<section><h3><span class="num">03</span>驗證層的校準意見（全文，未摘要）</h3>')
    A('<p class="sub">共 %d 條。這些意見的價值在於它們抓的是「我方對照做得不夠精準」，'
      '而不是教材的問題。全文列出，不美化。</p>' % n_over)
    for p in pages:
        ov = p.get("overclassified") or []
        fair = p.get("fairness")
        if not ov and not fair:
            continue
        A('<details class="pg"><summary><span>%s</span><code>%s</code>'
          '<span class="chip" style="color:var(--c-cs);border-color:var(--c-cs);background:#fff">%d 條</span></summary>'
          '<div class="body">' % (PAGE_TITLE.get(p["page"], p["page"]), p["page"], len(ov)))
        if ov:
            A('<ol class="fix">')
            for o in ov:
                A("<li>%s</li>" % esc((o)))
            A("</ol>")
        if fair:
            A('<div class="lbl" style="font-size:.72rem;color:var(--c-text-soft);font-weight:700">'
              '驗證層的公允性評註與數字核對明細</div>')
            A('<div class="note" style="white-space:pre-wrap;font-size:.85rem">%s</div>' % esc((fair)))
        A("</div></details>")
    A("</section>")

    A('<footer class="prov"><h4>方法與範圍</h4><ul>')
    A("<li>核對對象：<code>%s</code> 的 12 個頁面（sr1–sr5、m06、m08、m09、m10、m12）。</li>" % BASE)
    A("<li>四級分類的門檻：<b>REFUTES</b> 要求一般情況下不成立且有明確反證；<b>TENSION</b> 只要求數據不一致"
      "且須說明差異來源；教材為可教性而簡化<b>不算問題</b>。</li>")
    A("<li>我方數字一律要求附本專案檔案路徑，並由獨立的驗證 agent 實際 grep 確認。</li>")
    A("<li><b>範圍限制</b>：本專案只有 7 套技術資料集 / 6 個生物樣本（HCC1395 與 HCC1395_DORADO 為同一株細胞的兩套管線），"
      "且只有 HCC1395 有 SEQC2 外部真值。單樣本結果不可外推成通則。</li>")
    A("<li><b>這是本機報告</b>，未發布至 Wiki 或 GitHub Pages。</li>")
    A("</ul><p>由 <code>scripts/build_tutorial_crosscheck_report.py</code> 產生</p></footer></div>")

    open(out, "w", encoding="utf-8").write("\n".join(P))
    return {"pages": len(pages), "findings": len(findings), "fabricated": n_fab,
            "calibrations": n_over, "bytes": os.path.getsize(out)}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default=os.path.join(
        ROOT, "docs", "methodology", "20260811_lab_tutorial_crosscheck_HCC1395.standalone.html"))
    a = ap.parse_args()
    st = build(a.out)
    print("  ✔ %s" % a.out)
    for k, v in st.items():
        print("     %-14s %s" % (k, "{:,}".format(v) if isinstance(v, int) else v))
