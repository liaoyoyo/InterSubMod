#!/usr/bin/env python3
"""Subclone 重建軟體全景報告生成器（本機用）

來源：scratchpad/tools_survey_clean.json（6 類調研 agent + 6 類引用查證 agent）

🔴 引用紀律
調研層產出的每一筆文獻／repo，都由獨立的查證層用 PubMed / Crossref 反查。
查出 13 個問題，其中 4 個是 PMID/DOI 錯配（識別碼指向完全不相干的論文）——
這是文獻調研最危險也最難察覺的一型。本頁對所有非 EXISTS 的項目一律標紅，
並在最前面獨立成節，**不可直接引用未通過查證者**。

用法: python3 scripts/build_subclone_tool_landscape.py [-o OUT.html]
"""
import argparse
import collections
import html
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = ("/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/"
       "64a5bba6-e59f-4d9c-989b-9885a695e7bf/scratchpad/tools_survey_clean.json")

FIT = collections.OrderedDict([
    ("GOOD", ("#2F7D5B", "適用")),
    ("PARTIAL", ("#B06A16", "有條件")),
    ("POOR", ("#9A3B3B", "不適用")),
    ("NOT_APPLICABLE", ("#6B6B66", "不適用（前提不符）")),
])
CAT_NAME = {
    "vaf_cluster": "VAF-based 亞克隆分群", "phylo": "演化樹重建",
    "cn_purity": "拷貝數與純度前處理", "longread": "長讀長專用與單分子方法",
    "sc_multi": "單細胞與多區域", "bench": "評測與比較文獻",
}
TOFU = {"\U0001F534": "●", "\U0001F7E0": "●", "\U0001F7E1": "●", "\U0001F7E2": "●",
        "⭐": "★", "✅": "✔", "❌": "✘", "\U0001F50D": "", "\U0001F4CC": "", "\U0001F9EA": ""}

CSS = """
:root{--c-accent:#D97757;--c-accent-soft:#F2E4DC;--c-text:#1c1b19;--c-text-soft:#5C5A54;
--c-bg:#FAF9F5;--c-border:#E3DACC;--c-pass:#2F7D5B;--c-warn:#B06A16;--c-info:#2F6690;
--c-dead:#9A3B3B;--c-cs:#5B3FA0;--radius:10px;
--sans:-apple-system,system-ui,"Noto Sans CJK TC","PingFang TC","Microsoft JhengHei",sans-serif;
--mono:"JetBrains Mono",ui-monospace,Menlo,Consolas,monospace;}
*{box-sizing:border-box}
body{margin:0;font-family:var(--sans);background:var(--c-bg);color:var(--c-text);line-height:1.74;font-size:15.5px}
.wrap{max-width:1060px;margin:0 auto;padding:0 24px 60px}
header.hero{padding:34px 0 12px}
header.hero .kicker{font-family:var(--mono);font-size:.73rem;color:var(--c-accent);letter-spacing:1px;text-transform:uppercase}
header.hero h1{font-size:1.92rem;margin:6px 0;line-height:1.26}
header.hero .lede{color:var(--c-text-soft);max-width:870px;margin:0}
section{padding:26px 0;border-top:1px solid var(--c-border)}
section h3{font-size:1.36rem;margin:0 0 8px;display:flex;gap:11px;align-items:baseline;flex-wrap:wrap}
section h3 .num{font-family:var(--mono);color:var(--c-accent);font-weight:800;font-size:.9rem}
.intro{color:#3f3d38;font-size:.95rem;margin:0 0 14px;max-width:880px}
.bluf{background:var(--c-accent-soft);border:1px solid var(--c-accent);border-radius:var(--radius);padding:17px 20px;margin:14px 0}
.alarm{background:#F3E2E2;border:1px solid var(--c-dead);border-left:6px solid var(--c-dead);border-radius:var(--radius);padding:14px 18px;margin:16px 0}
.alarm b{color:var(--c-dead)}
.pick{background:#fff;border:1px solid var(--c-pass);border-left:6px solid var(--c-pass);border-radius:var(--radius);padding:13px 17px;margin:14px 0}
.pick .lab{display:block;font-family:var(--mono);font-size:.63rem;color:var(--c-pass);letter-spacing:1.4px;text-transform:uppercase;font-weight:800;margin-bottom:5px}
.kpis{display:flex;gap:11px;flex-wrap:wrap;margin:14px 0}
.kpi{flex:1;min-width:148px;background:#fff;border:1px solid var(--c-border);border-radius:var(--radius);padding:11px 15px}
.kpi .v{font-family:var(--mono);font-size:1.4rem;font-weight:800;line-height:1.1}
.kpi .k{font-size:.73rem;color:var(--c-text-soft);margin-top:4px}
table{border-collapse:collapse;width:100%;font-size:.83rem;background:#fff}
th,td{border:1px solid var(--c-border);padding:6px 9px;text-align:left;vertical-align:top}
th{background:#F3EFE7;font-weight:700}
.tw{overflow-x:auto;margin:12px 0}
.chip{display:inline-block;font-size:.66rem;font-family:var(--mono);font-weight:700;padding:1px 6px;border-radius:5px;border:1px solid;white-space:nowrap}
details.grp{margin:9px 0;border:1px solid var(--c-border);border-radius:var(--radius);background:#fff;overflow:hidden}
details.grp>summary{cursor:pointer;list-style:none;padding:10px 15px;font-weight:600;display:flex;gap:9px;align-items:center;flex-wrap:wrap;font-size:.94rem}
details.grp>summary::-webkit-details-marker{display:none}
details.grp>summary::before{content:"▸";color:var(--c-accent);font-weight:700;transition:transform .15s}
details.grp[open]>summary::before{transform:rotate(90deg)}
details.grp .body{padding:0 15px 12px}
.tool{border:1px solid var(--c-border);border-radius:8px;padding:10px 13px;margin:8px 0;background:#FBF9F4}
.tool .hd{display:flex;gap:8px;align-items:baseline;flex-wrap:wrap;margin-bottom:4px}
.tool .nm{font-weight:700;font-size:.98rem}
.tool .lb{font-size:.71rem;color:var(--c-text-soft);font-weight:700;margin-top:6px}
.tool p{margin:2px 0;font-size:.89rem}
.tool .cite{font-family:var(--mono);font-size:.72rem;color:var(--c-text-soft);margin-top:5px;overflow-wrap:anywhere}
code{font-family:var(--mono);font-size:.85em;background:#F0EDE6;padding:1px 4px;border-radius:4px}
figure{margin:16px 0;background:#fff;border:1px solid var(--c-border);border-radius:var(--radius);padding:14px}
figure svg{display:block;max-width:100%;height:auto;margin:0 auto}
figcaption{font-size:.79rem;color:var(--c-text-soft);margin-top:9px;border-top:1px dashed var(--c-border);padding-top:8px}
footer.prov{border-top:2px solid var(--c-border);margin-top:26px;padding:18px 0 40px;font-size:.78rem;color:var(--c-text-soft)}
footer.prov ul{padding-left:1.2em}
"""


def esc(x):
    t = html.escape(str(x if x is not None else ""))
    for a, b in TOFU.items():
        t = t.replace(a, b)
    return t


def fitbar(counts, width=840):
    tot = sum(counts.values()) or 1
    L, h = 128, 42 + len(FIT) * 29
    p = ['<svg viewBox="0 0 %d %d" role="img" aria-labelledby="ft fd" xmlns="http://www.w3.org/2000/svg">' % (width, h)]
    p.append('<title id="ft">工具對本專案資料的適配分佈</title>')
    p.append('<desc id="fd">81 個工具依適用、有條件、不適用三級對本專案 ONT 單一 bulk 資料的適配分佈。</desc>')
    p.append('<rect x="0" y="0" width="%d" height="%d" fill="#FAF9F5"/>' % (width, h))
    p.append('<text x="10" y="17" font-size="11.5" font-weight="700" fill="#1c1b19">'
             '%d 個工具對「ONT 長讀長 + 單一 bulk + CN-altered 95%%」的適配</text>' % tot)
    mx = max(list(counts.values()) + [1])
    for i, (k, (col, name)) in enumerate(FIT.items()):
        v = counts.get(k, 0)
        y = 30 + i * 29
        w = (width - L - 150) * (v / mx)
        p.append('<text x="%d" y="%.1f" text-anchor="end" font-size="10.5" fill="%s">%s</text>' % (L - 8, y + 13, col, name))
        if v:
            p.append('<rect x="%d" y="%.1f" width="%.2f" height="17" rx="3" fill="%s" opacity="0.85"/>' % (L, y, w, col))
        p.append('<text x="%.1f" y="%.1f" font-size="10" font-family="ui-monospace,monospace" fill="#5C5A54">%d (%.1f%%)</text>'
                 % (L + w + 6, y + 13, v, 100 * v / tot))
    p.append("</svg>")
    return "".join(p)


def build(out):
    r = json.load(open(SRC, encoding="utf-8"))
    cats = r["categories"]
    tools = [t for c in cats for t in (c.get("tools") or [])]
    bench = [(c["key"], b) for c in cats for b in (c.get("benchmarks") or [])]
    fitc = collections.Counter(t.get("fit_for_our_data") for t in tools)
    bad = [(c["key"], x) for c in cats for x in (c.get("cite_results") or []) if x.get("status") != "EXISTS"]
    n_checked = sum(c.get("cite_checked") or 0 for c in cats)

    P, A = [], None
    A = P.append
    A('<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">')
    A("<title>Subclone 重建軟體全景與適配分析</title><style>%s</style>" % CSS)
    A('<div class="wrap"><header class="hero">')
    A('<div class="kicker">工具全景調研 · 本機報告 · 2026-08-12</div>')
    A("<h1>Subclone 重建軟體：現況、限制、評測，以及哪個適合我們的資料</h1>")
    A('<p class="lede">六大類、%d 個工具、%d 篇評測文獻。'
      '每個工具都對照本專案的資料特性（<b>ONT 長讀長 · 單一 bulk · 6 個生物樣本 · CN-altered 95.44%%</b>）'
      '做適配判斷。</p></header>' % (len(tools), len(bench)))

    A('<div class="bluf"><b>總結：沒有任何現成工具「正好適合」我們的資料，這不是壞消息。</b>'
      '主流工具幾乎全部建立在三個我們不滿足的前提上 —— 短讀長、多樣本、以及可靠的 CN／purity 前處理。'
      '%d 個工具中被判定為<b>直接適用</b>的僅 %d 個；'
      '大多數落在「有條件可用（當正交參照或上游前處理）」。'
      '這正是本專案存在的理由：<b>ONT 單一 bulk ＋ 單分子共現 ＋ haplotag</b> 這個組合目前是空白。</div>'
      % (len(tools), fitc.get("GOOD", 0)))

    A('<div class="kpis">')
    for v, k, c in [(str(len(tools)), "調研工具數", "var(--c-info)"),
                    (str(len(bench)), "評測文獻", "var(--c-info)"),
                    (str(n_checked), "引用查證筆數", "var(--c-cs)"),
                    (str(len(bad)), "引用有問題", "var(--c-dead)")]:
        A('<div class="kpi"><div class="v" style="color:%s">%s</div><div class="k">%s</div></div>' % (c, v, k))
    A("</div>")
    A('<figure>%s<figcaption>適配判斷的依據是「工具的前提假設 vs 我們的資料特性」，'
      '不是工具本身的好壞。很多被判「不適用」的工具在它自己的場景裡是領域標準。</figcaption></figure>'
      % fitbar(fitc))

    # ── 引用查證（誠信優先，放最前）──
    A('<section><h3><span class="num">00</span>先講引用可靠度</h3>')
    A('<p class="intro">文獻調研最大的風險是引用幻覺。本次每一筆文獻與 repo 都由獨立的查證層'
      '用 PubMed / Crossref 反查，共驗 <b>%d 筆</b>，抓到 <b>%d 個問題</b>。</p>' % (n_checked, len(bad)))
    A('<div class="alarm"><b>其中最危險的一型是「識別碼錯配」—— 論文確實存在，但給的 PMID/DOI 指向完全不相干的研究。</b>'
      '例如 CITUP 的 DOI 實際指向一篇自然語言處理論文、Orchard 的 PMID 指向果蠅病毒研究、'
      'DeCiFer 的 PMID 指向多巴胺神經科學論文。這種錯誤在校稿時極難察覺，'
      '但審稿人一查就破。<b>下表所列項目在補正前不可直接引用。</b></div>')
    A('<div class="tw"><table><thead><tr><th style="width:10%">類別</th><th style="width:22%">項目</th>'
      '<th style="width:13%">狀態</th><th style="width:55%">問題與更正</th></tr></thead><tbody>')
    for ck, x in bad:
        st = x.get("status")
        col = "var(--c-dead)" if st == "NOT_FOUND" else "var(--c-warn)"
        note = x.get("note") or x.get("corrected_citation") or ""
        A('<tr><td><code>%s</code></td><td>%s</td>'
          '<td><span class="chip" style="color:%s;border-color:%s;background:#fff">%s</span></td>'
          '<td>%s</td></tr>' % (ck, esc(x.get("name")), col, col, esc(st), esc(note)[:600]))
    A("</tbody></table></div>")
    A('<p class="intro" style="font-size:.88rem">其餘 %d 筆通過查證。'
      '整體而言長讀長那一類的引用品質最好（21 筆全通過）。</p></section>'
      % (n_checked - len(bad)))

    # ── 逐類 ──
    for i, c in enumerate(cats, 1):
        ts = c.get("tools") or []
        bs = c.get("benchmarks") or []
        fc = collections.Counter(t.get("fit_for_our_data") for t in ts)
        A('<section><h3><span class="num">%02d</span>%s</h3>' % (i, esc(CAT_NAME.get(c["key"], c["key"]))))
        if c.get("summary"):
            A('<p class="intro">%s</p>' % esc(c["summary"]))
        A('<p class="intro" style="font-size:.87rem;color:var(--c-text-soft)">%d 個工具 ・ %s</p>'
          % (len(ts), " ・ ".join("%s %d" % (FIT[k][1], fc[k]) for k in FIT if fc.get(k))))
        if c.get("verdict"):
            A('<div class="pick"><span class="lab">對我們資料的裁決</span>%s</div>' % esc(c["verdict"]))
        if ts:
            A('<details class="grp"><summary>%d 個工具的逐項說明</summary><div class="body">' % len(ts))
            for t in ts:
                col, nm = FIT.get(t.get("fit_for_our_data"), ("#666", "?"))
                A('<div class="tool"><div class="hd">'
                  '<span class="chip" style="color:%s;border-color:%s;background:#fff">%s</span>'
                  '<span class="nm">%s</span>%s</div>'
                  % (col, col, nm, esc(t.get("name")),
                     ('<span style="font-size:.75rem;color:var(--c-text-soft)">%s</span>' % esc(t.get("year"))) if t.get("year") else ""))
                A("<p>%s</p>" % esc(t.get("what_it_does")))
                for lab, key in [("需要什麼輸入", "input_requirements"), ("優點", "strengths"),
                                 ("限制", "limitations"), ("對我們資料", "fit_reason"), ("維護狀態", "still_maintained")]:
                    if t.get(key):
                        A('<div class="lb">%s</div><p>%s</p>' % (lab, esc(t[key])))
                if t.get("citation") or t.get("repo"):
                    A('<div class="cite">%s %s</div>' % (esc(t.get("citation")), esc(t.get("repo"))))
                A("</div>")
            A("</div></details>")
        if bs:
            A('<details class="grp"><summary>%d 篇評測文獻</summary><div class="body">' % len(bs))
            for b in bs:
                A('<div class="tool"><div class="hd"><span class="nm">%s</span></div>' % esc(b.get("name")))
                for lab, key in [("比較了什麼", "what_it_compared"), ("關鍵發現", "key_finding"),
                                 ("該評測本身的限制", "caveat")]:
                    if b.get(key):
                        A('<div class="lb">%s</div><p>%s</p>' % (lab, esc(b[key])))
                if b.get("citation"):
                    A('<div class="cite">%s</div>' % esc(b["citation"]))
                A("</div>")
            A("</div></details>")
        if c.get("conflicts"):
            A('<details class="grp"><summary style="color:var(--c-dead)">與本專案既有結論的衝突與需修正處</summary>'
              '<div class="body"><p style="font-size:.89rem;white-space:pre-wrap">%s</p></div></details>'
              % esc(c["conflicts"]))
        A("</section>")

    A('<footer class="prov"><h4>方法與範圍</h4><ul>')
    A("<li>六類分別由獨立 agent 調研，再由第二層 agent 用 PubMed / Crossref 反查每一筆文獻與 repo 的存在性。</li>")
    A("<li>適配判斷的基準：本專案為 <b>ONT 長讀長</b>、<b>單一 bulk</b>、7 套技術資料集 / 6 個生物樣本、"
      "HCC1395 unit 層 <b>CN-altered 95.44%</b>（CN 中性且無 LOH 僅 4.56%）、"
      "有 germline haplotag 與 5mC、僅 HCC1395 有 SEQC2 外部真值。</li>")
    A("<li><b>本頁不做工具好壞的排名</b> —— 判斷的是「工具前提 vs 我們資料」的匹配度。"
      "被判不適用者多半在其原生場景裡是領域標準。</li>")
    A("<li>標為 NOT_FOUND / WRONG_DETAILS 的引用<b>在補正前不可使用</b>。</li>")
    A("<li><b>本機報告</b>，未發布至 Wiki 或 GitHub Pages。</li>")
    A("</ul><p>由 <code>scripts/build_subclone_tool_landscape.py</code> 產生</p></footer></div>")

    open(out, "w", encoding="utf-8").write("\n".join(P))
    return {"tools": len(tools), "benchmarks": len(bench), "cite_checked": n_checked,
            "cite_bad": len(bad), "fit_good": fitc.get("GOOD", 0), "bytes": os.path.getsize(out)}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default=os.path.join(
        ROOT, "docs", "method_comparison", "20260812_subclone_tool_landscape_and_fit.standalone.html"))
    a = ap.parse_args()
    st = build(a.out)
    print("  ✔ %s" % a.out)
    for k, v in st.items():
        print("     %-13s %s" % (k, "{:,}".format(v) if isinstance(v, int) else v))
