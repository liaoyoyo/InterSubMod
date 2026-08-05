#!/usr/bin/env python3
"""Render the within-HP near/distal stratification report to standalone HTML.

Every displayed number is injected from the sibling .data.json; a missing key
raises immediately rather than rendering a placeholder (CLAUDE.md 13-A).
No emoji (this machine has no emoji font) and no external requests.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "20260726_within_hp_near_distal.data.json")
OUT = os.path.join(HERE, "20260726_within_hp_near_distal.standalone.html")

with open(DATA) as fh:
    D = json.load(fh)


def req(path):
    """Fetch a dotted key path; refuse (raise) when absent."""
    cur = D
    for part in path.split("."):
        if isinstance(cur, list):
            cur = cur[int(part)]
        elif part in cur:
            cur = cur[part]
        else:
            raise KeyError("data.json missing required key: " + path)
    return cur


C_TP, C_TPN = "#0072B2", "#9ecae1"
C_FP, C_FPN = "#D55E00", "#F0B27A"
C_NEG, C_INK, C_MUT = "#B2182B", "#1b1f24", "#5b6672"


def mid(label):
    lo, hi = label.split("-")
    return (int(lo) + int(hi)) / 2.0


# --------------------------------------------------------------- chart A / B
def line_chart(series, ymax, title, ylab, height=330, excess=False):
    W, H = 780, height
    L, R, T, B = 66, 182, 26, 52
    pw, ph = W - L - R, H - T - B

    def sx(v):
        return L + v / 5200.0 * pw

    def sy(v):
        return T + ph - (v / ymax) * ph

    p = ['<svg viewBox="0 0 {0} {1}" width="100%" role="img" aria-label="{2}">'.format(W, H, title)]
    p.append('<rect x="{0}" y="{1}" width="{2}" height="{3}" fill="none" stroke="currentColor" stroke-opacity=".16"/>'.format(L, T, pw, ph))
    steps = 5
    for i in range(steps + 1):
        v = ymax * i / steps
        y = sy(v)
        p.append('<line x1="{0}" y1="{1:.1f}" x2="{2}" y2="{1:.1f}" stroke="currentColor" stroke-opacity=".09"/>'.format(L, y, L + pw))
        p.append('<text x="{0}" y="{1:.1f}" text-anchor="end" font-size="11" fill="{2}">{3:.3f}</text>'.format(L - 8, y + 3.5, C_MUT, v))
    for lab in [s["bin"] for s in series[0]["pts"]]:
        x = sx(mid(lab))
        p.append('<line x1="{0:.1f}" y1="{1}" x2="{0:.1f}" y2="{2}" stroke="currentColor" stroke-opacity=".07"/>'.format(x, T, T + ph))
        p.append('<text x="{0:.1f}" y="{1}" text-anchor="middle" font-size="10.5" fill="{2}">{3}</text>'.format(x, T + ph + 17, C_MUT, lab.replace("-", "–")))
    p.append('<text x="{0:.0f}" y="{1}" text-anchor="middle" font-size="11.5" fill="{2}">CpG 與 sSNV 的距離 (bp)</text>'.format(L + pw / 2, T + ph + 40, C_MUT))
    p.append('<text transform="translate(15,{0:.0f}) rotate(-90)" text-anchor="middle" font-size="11.5" fill="{1}">{2}</text>'.format(T + ph / 2, C_MUT, ylab))

    # shaded excess bands (observed minus null), drawn first
    if not excess:
        for s in series:
            if not s.get("band"):
                continue
            up = s["pts"]
            dn = s["band"]
            pts = ["{0:.1f},{1:.1f}".format(sx(mid(a["bin"])), sy(a["v"])) for a in up]
            pts += ["{0:.1f},{1:.1f}".format(sx(mid(b["bin"])), sy(b["v"])) for b in reversed(dn)]
            p.append('<polygon points="{0}" fill="{1}" fill-opacity=".13"/>'.format(" ".join(pts), s["color"]))

    for s in series:
        pts = " ".join("{0:.1f},{1:.1f}".format(sx(mid(a["bin"])), sy(a["v"])) for a in s["pts"])
        dash = ' stroke-dasharray="6 4"' if s.get("dash") else ""
        p.append('<polyline points="{0}" fill="none" stroke="{1}" stroke-width="{2}"{3} stroke-linejoin="round"/>'.format(pts, s["color"], s.get("w", 2.6), dash))
        for a in s["pts"]:
            p.append('<circle cx="{0:.1f}" cy="{1:.1f}" r="3.6" fill="{2}"/>'.format(sx(mid(a["bin"])), sy(a["v"]), s["color"]))

    ly = T + 6
    p.append('<rect x="{0}" y="{1}" width="{2}" height="{3}" rx="6" fill="currentColor" fill-opacity=".04" stroke="currentColor" stroke-opacity=".12"/>'.format(L + pw + 14, ly, R - 26, 24 * len(series) + 14))
    for i, s in enumerate(series):
        yy = ly + 20 + i * 24
        dash = ' stroke-dasharray="6 4"' if s.get("dash") else ""
        p.append('<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="{3}" stroke-width="2.8"{4}/>'.format(L + pw + 26, yy, L + pw + 52, s["color"], dash))
        p.append('<text x="{0}" y="{1}" font-size="11.5" fill="currentColor">{2}</text>'.format(L + pw + 58, yy + 4, s["name"]))
    p.append("</svg>")
    return "".join(p)


tp_c, fp_c = req("step3.tp"), req("step3.fp")
chart_a = line_chart(
    [
        {"name": "TP 觀察", "color": C_TP, "pts": [{"bin": c["bin"], "v": c["obs"]} for c in tp_c],
         "band": [{"bin": c["bin"], "v": c["null"]} for c in tp_c]},
        {"name": "TP 置換 null", "color": C_TPN, "dash": True, "w": 2.2,
         "pts": [{"bin": c["bin"], "v": c["null"]} for c in tp_c]},
        {"name": "FP 觀察（對照）", "color": C_FP, "pts": [{"bin": c["bin"], "v": c["obs"]} for c in fp_c],
         "band": [{"bin": c["bin"], "v": c["null"]} for c in fp_c]},
        {"name": "FP 置換 null", "color": C_FPN, "dash": True, "w": 2.2,
         "pts": [{"bin": c["bin"], "v": c["null"]} for c in fp_c]},
    ],
    0.08, "distance decay observed vs null", "|Δβ| 中位（每個 CpG）", 360)

chart_b = line_chart(
    [
        {"name": "TP 超額", "color": C_TP, "pts": [{"bin": c["bin"], "v": c["exc"]} for c in tp_c]},
        {"name": "FP 超額（對照）", "color": C_FP, "pts": [{"bin": c["bin"], "v": c["exc"]} for c in fp_c]},
    ],
    0.02, "excess over null", "超額 = 觀察 − null", 300, excess=True)


# ------------------------------------------------------------------ chart C
def matched_chart():
    s2 = req("step2")
    groups = [
        ("未配對 CpG 數", s2["raw_near"], s2["raw_distal"], s2["raw_p"],
         "near 只有 {0:.0f} 個 CpG、distal 有 {1:.0f} 個".format(s2["cpg_near"], s2["cpg_distal"])),
        ("配對 CpG 數後", s2["mat_near"], s2["mat_distal"], s2["mat_p"], "兩側取相同 CpG 數"),
    ]
    W, H = 780, 300
    L, T, ph, gw, bw = 66, 26, 190, 300, 74
    p = ['<svg viewBox="0 0 {0} {1}" width="100%" role="img" aria-label="matched near distal">'.format(W, H)]
    ymax = 0.04
    for i in range(5):
        v = ymax * i / 4
        y = T + ph - v / ymax * ph
        p.append('<line x1="{0}" y1="{1:.1f}" x2="{2}" y2="{1:.1f}" stroke="currentColor" stroke-opacity=".09"/>'.format(L, y, L + 660))
        p.append('<text x="{0}" y="{1:.1f}" text-anchor="end" font-size="11" fill="{2}">{3:.3f}</text>'.format(L - 8, y + 3.5, C_MUT, v))
    for gi, (gname, nv, dv, pv, note) in enumerate(groups):
        gx = L + 60 + gi * gw
        for bi, (lab, val, col) in enumerate((("near ≤1 kb", nv, C_TP), ("distal ≥2 kb", dv, C_MUT))):
            x = gx + bi * (bw + 18)
            h = val / ymax * ph
            p.append('<rect x="{0:.0f}" y="{1:.1f}" width="{2}" height="{3:.1f}" rx="3" fill="{4}" fill-opacity="{5}"/>'.format(x, T + ph - h, bw, h, col, .85 if bi == 0 else .45))
            p.append('<text x="{0:.0f}" y="{1:.1f}" text-anchor="middle" font-size="12" font-weight="600" fill="currentColor">{2:.4f}</text>'.format(x + bw / 2, T + ph - h - 7, val))
            p.append('<text x="{0:.0f}" y="{1}" text-anchor="middle" font-size="11" fill="{2}">{3}</text>'.format(x + bw / 2, T + ph + 17, C_MUT, lab))
        ratio = nv / dv
        verdict = "顯著" if pv < 0.05 else "不顯著"
        col = C_NEG if pv < 0.05 else "#2E7D32"
        p.append('<text x="{0:.0f}" y="{1}" text-anchor="middle" font-size="12.5" font-weight="700" fill="currentColor">{2}</text>'.format(gx + bw + 9, T + ph + 40, gname))
        p.append('<text x="{0:.0f}" y="{1}" text-anchor="middle" font-size="11.5" fill="{2}">比值 {3:.2f} · p={4:.2g} · {5}</text>'.format(gx + bw + 9, T + ph + 58, col, ratio, pv, verdict))
        p.append('<text x="{0:.0f}" y="{1}" text-anchor="middle" font-size="10.5" fill="{2}">{3}</text>'.format(gx + bw + 9, T + ph + 74, C_MUT, note))
    p.append('<text transform="translate(15,{0:.0f}) rotate(-90)" text-anchor="middle" font-size="11.5" fill="{1}">|Δβ| 中位（每個 位點×HP 單元）</text>'.format(T + ph / 2, C_MUT))
    p.append("</svg>")
    return "".join(p)


chart_c = matched_chart()

# ------------------------------------------------------------------- gates
GATES = [
    ("關卡 1", "germline ASM", "由設計排除",
     "只比較<b>同一個 germline HP 家族內部</b>的 ALT vs REF read。germline 單倍型專一甲基化在兩組中完全相同，相減即抵消。",
     "可測單元 {0:,} 個（位點×HP），其中 {1:,} 個 p&lt;0.05（{2:.1f}%）".format(req("step1.n_unit"), req("step1.n_sig"), req("step1.pct_sig")), "ok"),
    ("關卡 2", "CpG 取樣噪音", "已排除",
     "近端 CpG 數天生較少（中位 {0:.0f} 個 vs {1:.0f} 個），平均後噪音較大，會偽裝成「近端效應較強」。取相同 CpG 數重算即可檢驗。".format(req("step2.cpg_near"), req("step2.cpg_distal")),
     "比值由 {0:.2f}（p={1:.1e}）降到 {2:.2f}（p={3:.2f}，不顯著）".format(req("step2.raw_near") / req("step2.raw_distal"), req("step2.raw_p"), req("step2.mat_near") / req("step2.mat_distal"), req("step2.mat_p")), "ok"),
    ("關卡 3", "somatic-cis 足跡", "已排除",
     "點源的 cis 足跡必須<b>隨距離衰減</b>。以置換 null 校正後，超額值在 0.5 kb 到 5 kb 之間完全持平。",
     "超額 {0:+.4f} → {1:+.4f}（跨 5 個距離帶，無衰減）".format(tp_c[0]["exc"], tp_c[-1]["exc"]), "ok"),
    ("關卡 4", "somatic lineage 殘餘", "未通過",
     "殘存的均勻超額若源於腫瘤譜系，體細胞 TP 位點應<b>高於</b>非體細胞 FP 對照。實測相反。",
     "TP 超額 {0:+.4f}（比值 {1:.2f}）< FP 超額 {2:+.4f}（比值 {3:.2f}）".format(req("step4.tp.exc"), req("step4.tp.ratio"), req("step4.fp.exc"), req("step4.fp.ratio")), "bad"),
]

gate_html = []
for tag, name, verdict, why, ev, kind in GATES:
    cls = "gate ok" if kind == "ok" else "gate bad"
    gate_html.append(
        '<div class="{0}"><div class="gtop"><span class="gtag">{1}</span>'
        '<span class="gname">{2}</span><span class="gv">{3}</span></div>'
        '<p class="gwhy">{4}</p><p class="gev">{5}</p></div>'.format(cls, tag, name, verdict, why, ev))
gate_html = "".join(gate_html)

rows = "".join(
    '<tr><td class="mono">{0}</td><td>HP{1}</td><td class="num">{2:+.3f}</td>'
    '<td class="num">{3:.1e}</td><td class="num">{4} / {5}</td></tr>'.format(
        t["loc"], t["hp"], t["d"], t["p"], t["n_alt"], t["n_ref"])
    for t in req("step1.top"))

HTML = """<title>within-HP near vs distal 甲基分層 — HCC1395</title>
<style>
:root{{--bg:#fbfbfc;--fg:#1b1f24;--mut:#5b6672;--card:#fff;--line:#e3e6ea;
--ok:#2E7D32;--bad:#B2182B;--tp:#0072B2;--fp:#D55E00;--accent:#0072B2}}
@media (prefers-color-scheme:dark){{:root{{--bg:#14171a;--fg:#e8eaed;--mut:#9aa4af;
--card:#1c2024;--line:#2b3137;--ok:#7cc47f;--bad:#f08a80;--tp:#5aa9dd;--fp:#f0985c;--accent:#5aa9dd}}}}
:root[data-theme=dark]{{--bg:#14171a;--fg:#e8eaed;--mut:#9aa4af;--card:#1c2024;--line:#2b3137;
--ok:#7cc47f;--bad:#f08a80;--tp:#5aa9dd;--fp:#f0985c;--accent:#5aa9dd}}
:root[data-theme=light]{{--bg:#fbfbfc;--fg:#1b1f24;--mut:#5b6672;--card:#fff;--line:#e3e6ea;
--ok:#2E7D32;--bad:#B2182B;--tp:#0072B2;--fp:#D55E00;--accent:#0072B2}}
body{{background:var(--bg);color:var(--fg);font:15px/1.7 -apple-system,"Noto Sans CJK TC","Microsoft JhengHei",sans-serif;
margin:0;padding:32px 20px 80px}}
.wrap{{max-width:940px;margin:0 auto}}
h1{{font-size:26px;line-height:1.35;margin:0 0 6px}}
h2{{font-size:19px;margin:44px 0 6px;padding-top:16px;border-top:1px solid var(--line)}}
h2 .lv{{font-size:11px;font-weight:700;letter-spacing:.06em;color:var(--mut);
border:1px solid var(--line);border-radius:4px;padding:2px 6px;margin-right:8px;vertical-align:2px}}
.sub{{color:var(--mut);font-size:13px;margin:0 0 26px}}
.verdict{{background:var(--card);border:1px solid var(--line);border-left:5px solid var(--bad);
border-radius:10px;padding:18px 22px;margin:22px 0}}
.verdict .vt{{font-size:11px;font-weight:800;letter-spacing:.1em;color:var(--bad);margin-bottom:6px}}
.verdict p{{margin:0;font-size:16.5px;line-height:1.65}}
.kpis{{display:grid;grid-template-columns:repeat(auto-fit,minmax(196px,1fr));gap:12px;margin:20px 0}}
.kpi{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:14px 16px}}
.kpi .k{{font-size:11.5px;color:var(--mut);margin-bottom:5px}}
.kpi .v{{font-size:23px;font-weight:700;font-variant-numeric:tabular-nums;letter-spacing:-.01em}}
.kpi .n{{font-size:11.5px;color:var(--mut);margin-top:4px}}
.gate{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:14px 17px;margin:10px 0}}
.gate.ok{{border-left:5px solid var(--ok)}} .gate.bad{{border-left:5px solid var(--bad)}}
.gtop{{display:flex;align-items:center;gap:10px;flex-wrap:wrap}}
.gtag{{font-size:11px;font-weight:700;letter-spacing:.06em;color:var(--mut)}}
.gname{{font-size:15.5px;font-weight:700;flex:1}}
.gv{{font-size:12px;font-weight:700;padding:2px 9px;border-radius:11px;border:1px solid currentColor}}
.gate.ok .gv{{color:var(--ok)}} .gate.bad .gv{{color:var(--bad)}}
.gwhy{{margin:8px 0 6px;font-size:14px;color:var(--fg)}}
.gev{{margin:0;font-size:13px;color:var(--mut);font-variant-numeric:tabular-nums}}
figure{{margin:18px 0;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px 14px;overflow-x:auto}}
figcaption{{font-size:13.5px;color:var(--mut);margin-top:10px;padding:0 4px}}
figcaption b{{color:var(--fg)}}
table{{width:100%;border-collapse:collapse;font-size:13.5px;margin:12px 0}}
th,td{{padding:7px 10px;border-bottom:1px solid var(--line);text-align:left}}
th{{font-size:11.5px;color:var(--mut);font-weight:600;letter-spacing:.03em}}
.num{{text-align:right;font-variant-numeric:tabular-nums}}
.mono{{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px}}
details{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:12px 16px;margin:10px 0}}
summary{{cursor:pointer;font-weight:650;font-size:14.5px}}
details p,details li{{font-size:13.5px;color:var(--fg)}}
.note{{font-size:13px;color:var(--mut);border-left:3px solid var(--line);padding-left:12px;margin:14px 0}}
code{{background:rgba(127,127,127,.13);padding:1px 5px;border-radius:3px;font-size:12.5px}}
</style>
<div class="wrap">
<h1>within-HP 條件下的 near vs distal 甲基分層</h1>
<p class="sub">HCC1395 · {archive_short} · 產生日期 {date} · 全部數字由
<code>{datafile}</code> 注入</p>

<div class="verdict"><div class="vt">結論 — 陰性</div>
<p>在同一個 germline HP 家族內比較 ALT vs REF read，甲基差異<b>不隨距離衰減</b>（否定點源 cis 足跡），
但殘存的均勻超額在<b>非體細胞的 FP 對照位點上反而更大</b>。
因此 within-HP 的 ALT/REF 甲基差異<b>沒有可歸因於 somatic lineage 的成分</b>。</p></div>

<div class="kpis">
<div class="kpi"><div class="k">超額值跨距離變化</div><div class="v">{exc_lo:+.4f} → {exc_hi:+.4f}</div>
<div class="n">0.5 kb 至 5 kb，完全持平</div></div>
<div class="kpi"><div class="k">CpG 數配對後 near/distal</div><div class="v">{mat_ratio:.2f}</div>
<div class="n">配對前 {raw_ratio:.2f}（p={raw_p:.0e}）</div></div>
<div class="kpi"><div class="k">TP 對 null 的比值</div><div class="v">{tp_ratio:.2f}</div>
<div class="n">體細胞真陽性位點</div></div>
<div class="kpi"><div class="k">FP 對照對 null 的比值</div><div class="v">{fp_ratio:.2f}</div>
<div class="n">非體細胞 — 反而更高</div></div>
</div>

<h2><span class="lv">L1</span>四道關卡：三道通過，第四道不通過</h2>
<p class="sub" style="margin-bottom:14px">每道關卡各排除一個可能解釋。前三道成功把混淆源移除，
第四道用來檢驗「剩下的是不是 lineage」—— 而它沒有通過。</p>
{gates}

<h2><span class="lv">L2</span>證據一：甲基差異不隨距離衰減，否定點源 cis 足跡</h2>
<figure>{chart_a}
<figcaption>實線為觀察值，虛線為<b>置換 null</b>（把同一批 read 的 ALT/REF 標籤隨機重排後重算）。
陰影 = 超額。原始觀察值隨距離<i>微升</i>，但 null 同步上升（窗口邊緣覆蓋較低、噪音較大）——
兩者相減後的超額才是訊號，而它是平的。<b>點源 cis 足跡必然衰減，此處沒有衰減。</b>
同時注意橘色 FP 對照的陰影帶比藍色 TP <b>更寬</b>。</figcaption></figure>

<figure>{chart_b}
<figcaption>只畫超額值，把上圖的結論獨立出來。TP 五個距離帶為
{exc_list}，變異小於 0.001；FP 對照在<b>每一個距離帶都高於 TP</b>。</figcaption></figure>

<h2><span class="lv">L2</span>證據二：near &gt; distal 是 CpG 取樣噪音，不是生物訊號</h2>
<figure>{chart_c}
<figcaption>初步分析顯示 near 的 |Δβ| 明顯大於 distal（比值 {raw_ratio:.2f}，p={raw_p:.1e}），看似支持 cis 足跡。
但每個單元的 near 只有中位 {cpg_near:.0f} 個 CpG、distal 有 {cpg_distal:.0f} 個 —— 平均項數少一半以上，
離散度自然較大。<b>兩側取相同 CpG 數重算後差異消失</b>（比值 {mat_ratio:.2f}，p={mat_p:.2f}）。
這是本輪最容易誤判的一步。</figcaption></figure>

<h2><span class="lv">L2</span>證據三：體細胞位點沒有超過它自己的陰性對照</h2>
<table><thead><tr><th>位點集合</th><th class="num">位點池</th><th class="num">可測單元</th>
<th class="num">CpG 數</th><th class="num">觀察</th><th class="num">null</th>
<th class="num">超額</th><th class="num">比值</th></tr></thead><tbody>
<tr><td><b style="color:var(--tp)">TP</b> 體細胞真陽性</td><td class="num">{tp_pool:,}</td><td class="num">{tp_units:,}</td>
<td class="num">{tp_ncpg:,}</td><td class="num">{tp_obs:.4f}</td><td class="num">{tp_null:.4f}</td>
<td class="num">{tp_exc:+.4f}</td><td class="num">{tp_ratio:.2f}</td></tr>
<tr><td><b style="color:var(--fp)">FP</b> 非體細胞對照</td><td class="num">{fp_pool:,}</td><td class="num">{fp_units:,}</td>
<td class="num">{fp_ncpg:,}</td><td class="num">{fp_obs:.4f}</td><td class="num">{fp_null:.4f}</td>
<td class="num">{fp_exc:+.4f}</td><td class="num">{fp_ratio:.2f}</td></tr>
</tbody></table>
<p class="note">若均勻超額源自腫瘤譜系，TP 應高於 FP。實際 <b>TP − FP = {delta_exc:+.4f}</b>（比值差 {delta_ratio:+.2f}），
方向與 lineage 假說相反。FP 的絕對噪音較高（null {fp_null:.4f} vs {tp_null:.4f}），
所以應以<b>比值</b>比較 —— TP 在兩個指標上都沒有勝出。</p>

<h2><span class="lv">L2</span>個別位點仍有大效應，但無法歸因</h2>
<p class="sub" style="margin-bottom:8px">下列為 within-HP 檢定中 |Δ| 最大且顯著的位點。
它們是真實的資料特徵，可作<b>描述性註記</b>；但依上述三項證據，
不能宣稱這些差異來自 subclone 譜系。</p>
<table><thead><tr><th>位點</th><th>germline HP</th><th class="num">Δβ (ALT−REF)</th>
<th class="num">p</th><th class="num">ALT / REF read 數</th></tr></thead><tbody>{rows}</tbody></table>

<h2><span class="lv">L3</span>方法、定義與限制</h2>

<details><summary>名詞與座標軸定義</summary>
<ul>
<li><b>germline HP 家族</b>：HP tag 破折號前的第一碼。<code>1</code> 與 <code>1-1</code> 同族、
<code>2</code> 與 <code>2-1</code> 同族。同族內比較 ALT vs REF 時，germline 單倍型專一甲基化（ASM）在兩組完全相同，
因此自動抵消 —— 這正是專案既有的「HP1-1 vs HP1-2 = clonal 而非 ASM」定義。</li>
<li><b>Δβ</b>：同一個 CpG 上，ALT read 的平均甲基化值減去 REF read 的平均值。範圍 −1 到 +1，
本文一律取絕對值中位數（分布右偏，平均值受極端值主導）。</li>
<li><b>置換 null</b>：保持 read 集合與兩組人數不變，隨機重排 ALT/REF 標籤後重算 |Δβ|。
它回答的是「ALT/REF 這個切法，是否比同一批 read 的任意切法差更多」。</li>
<li><b>超額</b> = 觀察中位 − null 中位。<b>比值</b> = 觀察中位 ÷ null 中位（尺度無關，用於跨集合比較）。</li>
<li><b>near / distal</b>：CpG 距 sSNV ≤1,000 bp 為 near、≥2,000 bp 為 distal（中間留白避免邊界效應）。</li>
</ul></details>

<details><summary>統計設計與門檻</summary>
<ul>
<li>窗口 {window}；每組至少 {minr} 條 tumor read；每個 CpG 兩組各至少 {mincov} 條有效讀值。</li>
<li>檢定：{stat}。距離曲線的比較單位是 CpG，near/distal 的比較單位是 位點×HP 單元。</li>
<li>每個單元做 {nperm} 次置換；在 10<sup>5</sup> 量級的 CpG 上 null 中位已穩定。這是<b>集合層</b>的 null，
不是逐位點的顯著性校正。</li>
<li>取樣：位點以固定亂數種子（20260726）隨機打亂後掃描至額度為止，非全量。
TP 掃描 {tp_scanned:,} 個位點、FP {fp_scanned:,} 個。</li>
</ul></details>

<details><summary>限制（必讀）</summary>
<ul>
<li><b>FP 不是純 germline-het null。</b>FP 集合富含 artifact，絕對噪音較高。它是<b>方向正確</b>的對照
（若 lineage 有貢獻，TP 應超過 FP），但不能取代正式的 germline-het 對照。</li>
<li><b>單一樣本、單一資料線</b>（HCC1395、2026-01 w5000 archive）。跨樣本未驗證。</li>
<li><b>覆蓋深度偏誤</b>：要求同一 HP 內各有 ≥{minr} 條 ALT 與 REF read，
{tp_yield:.1f}% 的 TP 位點才可測 —— 結果偏向高覆蓋區。</li>
<li><b>CN 未在本輪控制</b>（依當輪指示暫置）。CN gain 會同時影響 ALT read 比例與甲基分布，
惟本輪的關鍵比較是 TP vs FP，兩者都受 CN 影響，方向性結論不因此翻轉。</li>
<li>本分析回答的是「ALT/REF 這個切法是否帶額外甲基訊號」，<b>不是</b>「這些區域有沒有甲基差異」。
後者為真（見既有 ASM 結果），但那屬 germline haplotype 軸。</li>
</ul></details>

<details><summary>資料來源與可重現性</summary>
<p>archive：<code>{archive}</code></p>
<p>來源 JSON（已隨報告封存於 <code>sources/</code>）：{srcs}</p>
<p>重建：<code>python3 assemble_data.py &amp;&amp; python3 build_page.py</code></p>
<p>本頁所有數字皆由 <code>{datafile}</code> 注入；缺少必要 key 時 builder 直接拋錯，不渲染占位值。</p>
</details>

<h2><span class="lv">L3</span>這對既有結論的影響</h2>
<p>與 <code>project_methylation_use_exhausted_bounded_auxiliary</code>（甲基用途已窮盡 = bounded-auxiliary）
及 <code>project_methyl_allele_axis_backbone_coenrichment_refuted</code>（allele 軸 × 骨幹共同富集已否決）
<b>一致</b>，並補上先前缺的一環：那兩份記錄把 allele 軸訊號歸因於 somatic-cis 足跡，
本輪直接檢定該假設，發現<b>連 cis 足跡也不成立</b> —— 距離曲線是平的，
而殘餘的均勻成分在非體細胞對照上更強。結論方向不變（甲基不作為 subclone 譜系證據），
但歸因需修正為「一般性的 allele-split 效應」，而非「突變局部 cis 足跡」。</p>
<p class="note">重新開啟條件：具備正式 germline-het null、CN-neutral 且足量的樣本，
且在<b>固定局部突變密度</b>下 TP 的比值仍顯著高於對照。</p>
</div>
""".format(
    archive_short=os.path.basename(req("meta.archive")),
    archive=req("meta.archive"),
    date=req("meta.date"),
    datafile=os.path.basename(DATA),
    srcs="、".join("<code>{0}</code>".format(s) for s in req("meta.sources")),
    gates=gate_html, chart_a=chart_a, chart_b=chart_b, chart_c=chart_c, rows=rows,
    exc_lo=tp_c[0]["exc"], exc_hi=tp_c[-1]["exc"],
    exc_list="、".join("{0:+.4f}".format(c["exc"]) for c in tp_c),
    mat_ratio=req("step2.mat_near") / req("step2.mat_distal"),
    raw_ratio=req("step2.raw_near") / req("step2.raw_distal"),
    raw_p=req("step2.raw_p"), mat_p=req("step2.mat_p"),
    cpg_near=req("step2.cpg_near"), cpg_distal=req("step2.cpg_distal"),
    tp_ratio=req("step4.tp.ratio"), fp_ratio=req("step4.fp.ratio"),
    tp_pool=req("step4.tp.pool"), fp_pool=req("step4.fp.pool"),
    tp_units=req("step4.tp.units"), fp_units=req("step4.fp.units"),
    tp_ncpg=req("step4.tp.n_cpg"), fp_ncpg=req("step4.fp.n_cpg"),
    tp_obs=req("step4.tp.obs"), fp_obs=req("step4.fp.obs"),
    tp_null=req("step4.tp.null"), fp_null=req("step4.fp.null"),
    tp_exc=req("step4.tp.exc"), fp_exc=req("step4.fp.exc"),
    delta_exc=req("step4.delta_excess"), delta_ratio=req("step4.delta_ratio"),
    tp_scanned=req("step4.tp.scanned"), fp_scanned=req("step4.fp.scanned"),
    tp_yield=req("step4.tp_yield"),
    window=req("design.window"), minr=req("design.min_reads_per_group"),
    mincov=req("design.min_cov_per_cpg"), nperm=req("design.n_perm"),
    stat=req("design.stat"),
)

with open(OUT, "w") as fh:
    fh.write(HTML)
print("wrote", OUT, os.path.getsize(OUT), "bytes")
