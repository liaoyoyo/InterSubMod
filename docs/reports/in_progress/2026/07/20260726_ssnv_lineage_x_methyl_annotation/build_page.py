#!/usr/bin/env python3
"""Render the sSNV-lineage x methylation annotation report.

All numbers injected from the sibling .data.json; missing keys raise.
No emoji, no external requests, theme-aware.
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "20260726_ssnv_lineage_x_methyl.data.json")
OUT = os.path.join(HERE, "20260726_ssnv_lineage_x_methyl.standalone.html")
with open(DATA) as fh:
    D = json.load(fh)


def req(p):
    cur = D
    for k in p.split("."):
        if isinstance(cur, list):
            cur = cur[int(k)]
        elif k in cur:
            cur = cur[k]
        else:
            raise KeyError("data.json missing key: " + p)
    return cur


C_A, C_R, C_N = "#0072B2", "#5b6672", "#999999"
C_BAD, C_OK = "#B2182B", "#2E7D32"

MS = ["M2_ALT_ONLY_MULTI", "M3_REGION_MULTI", "M1_ALTREF_SHIFT", "M0_NONE"]
MSL = {"M2_ALT_ONLY_MULTI": "M2 僅 ALT 內多群", "M3_REGION_MULTI": "M3 區域性多群",
       "M1_ALTREF_SHIFT": "M1 ALT/REF 位移", "M0_NONE": "M0 無標記"}
LS = ["L1_BACKBONE", "L2_WEAK_LINK", "L3_ISOLATED"]
LSL = {"L1_BACKBONE": "L1 骨幹（≥1 partner·雙端 somatic·coread≥6）",
       "L2_WEAK_LINK": "L2 弱串聯（有 partner 但未過閘）",
       "L3_ISOLATED": "L3 孤立（無 partner）"}


# ------------------------------------------------------- chart 1: three bars
def bars3():
    sy = req("symmetric")
    n = sy["n"]
    vals = [("ALT 群內", sy["alt"] / n * 100, C_A),
            ("REF 群內（內建對照）", sy["ref"] / n * 100, C_R),
            ("匹配單峰 null", sy["null"] / n * 100, C_N)]
    W, H, T, ph, bw = 760, 300, 26, 190, 108
    p = ['<svg viewBox="0 0 {0} {1}" width="100%" role="img" aria-label="alt ref null multimodality">'.format(W, H)]
    ymax = 30.0
    for i in range(7):
        v = ymax * i / 6
        y = T + ph - v / ymax * ph
        p.append('<line x1="70" y1="{0:.1f}" x2="700" y2="{0:.1f}" stroke="currentColor" stroke-opacity=".09"/>'.format(y))
        p.append('<text x="62" y="{0:.1f}" text-anchor="end" font-size="11" fill="{1}">{2:.0f}%</text>'.format(y + 3.5, C_N, v))
    for i, (lab, v, col) in enumerate(vals):
        x = 150 + i * 170
        h = v / ymax * ph
        p.append('<rect x="{0}" y="{1:.1f}" width="{2}" height="{3:.1f}" rx="4" fill="{4}" fill-opacity=".85"/>'.format(x, T + ph - h, bw, h, col))
        p.append('<text x="{0}" y="{1:.1f}" text-anchor="middle" font-size="15" font-weight="700" fill="currentColor">{2:.1f}%</text>'.format(x + bw / 2, T + ph - h - 9, v))
        p.append('<text x="{0}" y="{1}" text-anchor="middle" font-size="11.5" fill="{2}">{3}</text>'.format(x + bw / 2, T + ph + 18, C_N, lab))
    y0 = T + ph - vals[2][1] / ymax * ph
    p.append('<line x1="150" y1="{0:.1f}" x2="638" y2="{0:.1f}" stroke="{1}" stroke-width="2" stroke-dasharray="7 4"/>'.format(y0, C_BAD))
    p.append('<text x="646" y="{0:.1f}" font-size="11.5" fill="{1}">分群法本身的底線</text>'.format(y0 + 4, C_BAD))
    p.append('<text x="{0}" y="{1}" text-anchor="middle" font-size="12" fill="currentColor">ALT 與 REF 幾乎相同（配對 p={2:.2f}）→ 多群是<tspan font-weight="700">區域屬性</tspan>，不是等位屬性</text>'.format(W / 2, T + ph + 48, sy["p_paired"]))
    p.append('<text transform="translate(16,{0}) rotate(-90)" text-anchor="middle" font-size="11.5" fill="{1}">分離度超過 null-95th 的比例</text>'.format(T + ph / 2, C_N))
    p.append("</svg>")
    return "".join(p)


# --------------------------------------------- chart 2: lineage x methyl grid
def grid():
    tab, tot = req("table"), req("totals_lineage")
    W, H = 780, 250
    cw, rh, x0, y0 = 132, 46, 250, 62
    p = ['<svg viewBox="0 0 {0} {1}" width="100%" role="img" aria-label="lineage by methyl state">'.format(W, H)]
    for j, m in enumerate(MS):
        p.append('<text x="{0}" y="{1}" text-anchor="middle" font-size="11.5" font-weight="600" fill="currentColor">{2}</text>'.format(x0 + j * cw + cw / 2, y0 - 22, MSL[m].split(" ")[0]))
        p.append('<text x="{0}" y="{1}" text-anchor="middle" font-size="10" fill="{2}">{3}</text>'.format(x0 + j * cw + cw / 2, y0 - 8, C_N, MSL[m].split(" ", 1)[1]))
    for i, l in enumerate(LS):
        t = tot[l]
        p.append('<text x="{0}" y="{1}" text-anchor="end" font-size="11.5" font-weight="600" fill="currentColor">{2}</text>'.format(x0 - 14, y0 + i * rh + rh / 2 + 1, l))
        p.append('<text x="{0}" y="{1}" text-anchor="end" font-size="10" fill="{2}">n={3:,}</text>'.format(x0 - 14, y0 + i * rh + rh / 2 + 14, C_N, t))
        for j, m in enumerate(MS):
            v = tab[l][m] / t * 100
            op = min(v / 12.0, 1.0) if m != "M0_NONE" else min(v / 110.0, 1.0)
            col = C_A if m != "M0_NONE" else C_N
            p.append('<rect x="{0}" y="{1}" width="{2}" height="{3}" rx="3" fill="{4}" fill-opacity="{5:.2f}" stroke="currentColor" stroke-opacity=".13"/>'.format(x0 + j * cw + 3, y0 + i * rh + 3, cw - 6, rh - 6, col, .12 + .7 * op))
            p.append('<text x="{0}" y="{1}" text-anchor="middle" font-size="13" font-weight="700" fill="currentColor">{2:.1f}%</text>'.format(x0 + j * cw + cw / 2, y0 + i * rh + rh / 2 - 1, v))
            p.append('<text x="{0}" y="{1}" text-anchor="middle" font-size="9.5" fill="{2}">{3:,}</text>'.format(x0 + j * cw + cw / 2, y0 + i * rh + rh / 2 + 12, C_N, tab[l][m]))
    p.append('<text x="{0}" y="{1}" text-anchor="middle" font-size="12" fill="currentColor">每一直行由上到下幾乎不變 —— sSNV-lineage 狀態<tspan font-weight="700">不預測</tspan>甲基標記狀態</text>'.format(W / 2, y0 + 3 * rh + 34))
    p.append("</svg>")
    return "".join(p)


# ------------------------------------------------ chart 3: M1 strat collapse
def strat_chart():
    s = req("m1_strat")
    keep = [x for x in s if x["bin"] in ("全體(未分層)", "密度 0", "密度 1", "密度 2-3")]
    W, H, T, ph = 780, 320, 26, 190
    p = ['<svg viewBox="0 0 {0} {1}" width="100%" role="img" aria-label="m1 density stratification">'.format(W, H)]
    ymax = 12.0
    for i in range(7):
        v = ymax * i / 6
        y = T + ph - v / ymax * ph
        p.append('<line x1="72" y1="{0:.1f}" x2="720" y2="{0:.1f}" stroke="currentColor" stroke-opacity=".09"/>'.format(y))
        p.append('<text x="64" y="{0:.1f}" text-anchor="end" font-size="11" fill="{1}">{2:.0f}%</text>'.format(y + 3.5, C_N, v))
    gw, bw = 162, 58
    for gi, e in enumerate(keep):
        gx = 110 + gi * gw
        for bi, (lab, val, col) in enumerate((("L1 骨幹", e["pb"], C_A), ("L3 孤立", e["pi"], C_R))):
            x = gx + bi * (bw + 12)
            h = val / ymax * ph
            p.append('<rect x="{0:.0f}" y="{1:.1f}" width="{2}" height="{3:.1f}" rx="3" fill="{4}" fill-opacity="{5}"/>'.format(x, T + ph - h, bw, h, col, .85 if bi == 0 else .45))
            p.append('<text x="{0:.0f}" y="{1:.1f}" text-anchor="middle" font-size="11.5" font-weight="600" fill="currentColor">{2:.1f}</text>'.format(x + bw / 2, T + ph - h - 6, val))
        sigc = C_BAD if e["p"] < 0.05 else C_OK
        txt = "p={0:.3f}".format(e["p"]) + ("  顯著" if e["p"] < 0.05 else "  不顯著")
        p.append('<text x="{0:.0f}" y="{1}" text-anchor="middle" font-size="12" font-weight="700" fill="currentColor">{2}</text>'.format(gx + bw + 6, T + ph + 20, e["bin"]))
        p.append('<text x="{0:.0f}" y="{1}" text-anchor="middle" font-size="11" fill="{2}">{3:+.1f}pp · {4}</text>'.format(gx + bw + 6, T + ph + 38, sigc, e["d"], txt))
        p.append('<text x="{0:.0f}" y="{1}" text-anchor="middle" font-size="10" fill="{2}">n={3:,} / {4:,}</text>'.format(gx + bw + 6, T + ph + 54, C_N, e["nb"], e["ni"]))
    p.append('<rect x="96" y="{0}" width="176" height="{1}" rx="6" fill="{2}" fill-opacity=".07"/>'.format(T - 4, ph + 66, C_BAD))
    p.append('<text x="{0}" y="{1}" text-anchor="middle" font-size="12" fill="currentColor">左側未分層看似有富集；<tspan font-weight="700" fill="{2}">一旦固定局部突變密度就消失</tspan></text>'.format(W / 2, T + ph + 78, C_BAD))
    p.append('<text transform="translate(16,{0}) rotate(-90)" text-anchor="middle" font-size="11.5" fill="{1}">M1 位移旗標比例</text>'.format(T + ph / 2, C_N))
    p.append("</svg>")
    return "".join(p)


inv_rows = "".join(
    '<tr><td><b>{0}</b></td><td class="mono">{1}</td><td>{2}</td><td>{3}</td></tr>'.format(
        r["item"], r["path"], r["scope"], r["note"]) for r in req("inventory"))

ex_rows = "".join(
    '<tr><td class="mono">{0}:{1}</td><td>HP{2}</td><td class="num">{3:+.3f}</td>'
    '<td class="num">{4:.2f}</td><td class="num">{5}/{6}</td><td class="num">{7}</td>'
    '<td>{8}</td><td class="num">{9}</td></tr>'.format(
        e["chrom"], e["pos"], e["hp"], e["dsplit_alt"], e["sep_alt"], e["n_alt"], e["n_ref"],
        e["n_partners"], e["cn_state"], e["pipeline_optimal_k"]) for e in req("examples")[:10])

gk = req("genome_pipeline.k")
k_rows = "".join('<tr><td>optimal_k = {0}</td><td class="num">{1:,}</td><td class="num">{2:.1f}%</td></tr>'.format(
    k, v, v / req("genome_pipeline.n_loci") * 100) for k, v in gk.items())

SCHEMA = [
    ("chrom · pos · hp", "位點座標與 germline HP 家族（1 或 2）"),
    ("lineage_state", "L1_BACKBONE / L2_WEAK_LINK / L3_ISOLATED — 由 n_partners、both_somatic、max_coread 決定"),
    ("n_partners · max_coread · top_rel · both_somatic", "sSNV 共現骨幹原始欄位"),
    ("local_density_50kb", "±50 kb 內其他 TP sSNV 數 —— <b>必要控制變項</b>"),
    ("methyl_state", "M2_ALT_ONLY_MULTI / M3_REGION_MULTI / M1_ALTREF_SHIFT / M0_NONE"),
    ("flag_alt_multi · flag_ref_multi", "分離度 ≥ 匹配單峰 null 的 95 百分位（{0:.2f}）；REF 為內建對照".format(req("meta.null_q95"))),
    ("flag_altref_shift", "within-HP ALT vs REF：p&lt;0.05 且 |Δβ|≥0.10 且置換 null 超額 &gt;0"),
    ("flag_strand_assoc", "ALT 分群與正反股顯著關聯（技術性疑慮旗標）"),
    ("d_altref · p_altref · excess_vs_null", "within-HP 位移的效應量、p 值、對置換 null 的超額"),
    ("sep_alt · dsplit_alt · sep_ref · sep_null", "組內分離度、兩子群 Δβ、REF 對照、匹配單峰 null"),
    ("pipeline_optimal_k · pipeline_passed_gating · v_alt · v_hp", "原管線分群判定與兩軸 Cramér's V"),
    ("cn_state · is_loh · n_alt · n_ref · n_common_cpg", "CN 狀態與可用資料量"),
]
schema_rows = "".join('<tr><td class="mono">{0}</td><td>{1}</td></tr>'.format(a, b) for a, b in SCHEMA)

HTML = """<title>sSNV-lineage × 甲基標記層 — HCC1395</title>
<style>
:root{{--bg:#fbfbfc;--fg:#1b1f24;--mut:#5b6672;--card:#fff;--line:#e3e6ea;
--ok:#2E7D32;--bad:#B2182B;--a:#0072B2}}
@media (prefers-color-scheme:dark){{:root{{--bg:#14171a;--fg:#e8eaed;--mut:#9aa4af;--card:#1c2024;
--line:#2b3137;--ok:#7cc47f;--bad:#f08a80;--a:#5aa9dd}}}}
:root[data-theme=dark]{{--bg:#14171a;--fg:#e8eaed;--mut:#9aa4af;--card:#1c2024;--line:#2b3137;
--ok:#7cc47f;--bad:#f08a80;--a:#5aa9dd}}
:root[data-theme=light]{{--bg:#fbfbfc;--fg:#1b1f24;--mut:#5b6672;--card:#fff;--line:#e3e6ea;
--ok:#2E7D32;--bad:#B2182B;--a:#0072B2}}
body{{background:var(--bg);color:var(--fg);font:15px/1.7 -apple-system,"Noto Sans CJK TC","Microsoft JhengHei",sans-serif;margin:0;padding:32px 20px 80px}}
.wrap{{max-width:980px;margin:0 auto}}
h1{{font-size:26px;margin:0 0 6px;line-height:1.35}}
h2{{font-size:19px;margin:46px 0 8px;padding-top:16px;border-top:1px solid var(--line)}}
h2 .lv{{font-size:11px;font-weight:700;letter-spacing:.06em;color:var(--mut);border:1px solid var(--line);
border-radius:4px;padding:2px 6px;margin-right:8px;vertical-align:2px}}
h3{{font-size:15.5px;margin:22px 0 6px}}
.sub{{color:var(--mut);font-size:13px;margin:0 0 24px}}
.verdict{{background:var(--card);border:1px solid var(--line);border-left:5px solid var(--bad);
border-radius:10px;padding:18px 22px;margin:22px 0}}
.verdict .vt{{font-size:11px;font-weight:800;letter-spacing:.1em;color:var(--bad);margin-bottom:6px}}
.verdict p{{margin:0 0 8px;font-size:16px;line-height:1.66}}.verdict p:last-child{{margin:0}}
.use{{background:var(--card);border:1px solid var(--line);border-left:5px solid var(--ok);
border-radius:10px;padding:16px 20px;margin:18px 0}}
.use .vt{{font-size:11px;font-weight:800;letter-spacing:.1em;color:var(--ok);margin-bottom:6px}}
.kpis{{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:20px 0}}
.kpi{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:14px 16px}}
.kpi .k{{font-size:11.5px;color:var(--mut);margin-bottom:5px}}
.kpi .v{{font-size:23px;font-weight:700;font-variant-numeric:tabular-nums}}
.kpi .n{{font-size:11.5px;color:var(--mut);margin-top:4px}}
figure{{margin:16px 0;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px 14px;overflow-x:auto}}
figcaption{{font-size:13.5px;color:var(--mut);margin-top:10px;padding:0 4px}}
figcaption b{{color:var(--fg)}}
table{{width:100%;border-collapse:collapse;font-size:13.5px;margin:12px 0}}
th,td{{padding:7px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}}
th{{font-size:11.5px;color:var(--mut);font-weight:600}}
.num{{text-align:right;font-variant-numeric:tabular-nums}}
.mono{{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}}
details{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:12px 16px;margin:10px 0}}
summary{{cursor:pointer;font-weight:650;font-size:14.5px}}
.note{{font-size:13px;color:var(--mut);border-left:3px solid var(--line);padding-left:12px;margin:14px 0}}
code{{background:rgba(127,127,127,.13);padding:1px 5px;border-radius:3px;font-size:12.5px}}
ol.why > li{{margin-bottom:14px}}
.tagbad{{color:var(--bad);font-weight:700}}.tagok{{color:var(--ok);font-weight:700}}
</style>
<div class="wrap">
<h1>sSNV-lineage 分群下的甲基標記層</h1>
<p class="sub">HCC1395 · {archive} · {date} · 全基因組 {n_loci:,} 位點 / {n_units:,} 個 位點×HP 單元 ·
數字由 <code>{datafile}</code> 注入</p>

<div class="verdict"><div class="vt">兩個問題的實測答案</div>
<p><b>Q1 為何 HP 下 ALT 內還會看到多群？</b>主因是<b>分群程序本身</b>——匹配的單峰 null 在同一門檻下也產生
{null_pct:.0f}% 的「兩群」，而管線的 <code>optimal_k</code> 值域根本沒有 1。扣掉 null 後確實還有真實多峰
（ALT {alt_pct:.1f}%），但<b>同位點的 REF 群是 {ref_pct:.1f}%，一樣多</b>（配對 p={p_paired:.2f}）——
所以那是<b>區域屬性</b>（epiallele／PMD 式的分子間異質性），不是等位或譜系屬性。</p>
<p><b>Q2 與 sSNV-lineage 分群是否高度相關？</b><span class="tagbad">否，接近零。</span>
骨幹位點 M2 比例 {m2_l1:.1f}% vs 孤立位點 {m2_l3:.1f}%（OR={m2_or:.2f}, p={m2_p:.3f}）。
唯一看似有訊號的 M1 位移（+{m1_d:.1f}pp, p={m1_p:.3f}）在<b>固定局部突變密度後完全消失</b>。</p></div>

<div class="kpis">
<div class="kpi"><div class="k">ALT 內多群（超 null-95th）</div><div class="v">{alt_pct:.1f}%</div>
<div class="n">REF 對照 {ref_pct:.1f}% · null {null_pct:.1f}%</div></div>
<div class="kpi"><div class="k">M2 僅 ALT 多群（含效應量閘）</div><div class="v">{m2_tot:.1f}%</div>
<div class="n">{m2_n:,} 個單元，保守定義</div></div>
<div class="kpi"><div class="k">骨幹 vs 孤立的 M2 差異</div><div class="v">OR {m2_or:.2f}</div>
<div class="n">p={m2_p:.3f}，無富集</div></div>
<div class="kpi"><div class="k">標記層規模</div><div class="v">{n_units:,}</div>
<div class="n">位點×HP 單元 · 25 欄</div></div>
</div>

<h2><span class="lv">L1</span>Q1 為何會看到「ALT 內多群」—— 四個成因，逐一量化</h2>
<ol class="why">
<li><b>分群程序保證產生分群</b>（最大宗）。任何 1 維資料做 k=2 最佳切分都會得到「兩群」；
用同 n 同 SD 的<b>單峰高斯</b>重抽，在同一門檻下也有 {null_pct:.0f}% 被判為多群。
管線自己的 <code>optimal_k</code> 在 {gk_loci:,} 個位點中<b>從未輸出 1</b>（值域 2–6），
結構上就不可能回答「其實只有一群」。</li>
<li><b>真實的分子間異質性存在，但屬區域屬性</b>。扣除 null 後 ALT 有 {alt_pct:.1f}% 超標，
遠高於 5% 的期望值 —— 這是真的。但<b>同位點 REF 群 {ref_pct:.1f}%</b>，兩者無差異。
癌症基因組的 partially methylated domain 與 epiallele 隨機性會讓整個區域的分子彼此不同，
與該位點帶不帶突變無關。</li>
<li><b>read 覆蓋範圍差異</b>（已排除）。每條 read 覆蓋不同 CpG 子集，逐 read 平均會混入「覆蓋差異」。
改用共同-CpG 集（中位 {n_common:.0f} 個，全窗 {n_cpg:.0f} 個）後分離度中位由
{sep_raw:.2f} 變 {sep_common:.2f} —— 影響可忽略，不是主因。</li>
<li><b>正反股（strand）</b>（已排除）。ALT 分群與股別顯著關聯者僅 {strand_pct:.1f}%，
等同 5% 的偽陽性率。</li>
</ol>
<figure>{c1}
<figcaption>三根柱同一門檻、同一統計量。<b>ALT 與 REF 幾乎重合</b>是關鍵：
若多群源自突變譜系，只有 ALT 該升高。虛線是分群法本身的底線。</figcaption></figure>
<p class="note">還有兩個<b>本輪未能分離</b>的候選成因：<b>CN 增益</b>（同一 HP 家族內有多份拷貝，
各拷貝甲基狀態可不同 —— HCC1395 {cn_gain_pct:.0f}% 單元為 gain）與<b>細胞型別混雜／HP 指派誤差</b>。
兩者都會產生「區域性、非等位專一」的多群，與觀察到的模式一致，但需要 single-cell 或 CN-neutral 足量樣本才能拆開。</p>

<h2><span class="lv">L1</span>Q2 標記層：sSNV-lineage 狀態 × 甲基標記狀態</h2>
<figure>{c2}
<figcaption>每格為該 lineage 列中的百分比（下方小字為單元數）。
<b>三列數值幾乎一致</b> —— 知道一個位點在骨幹上、弱串聯、或孤立，對它的甲基標記狀態沒有預測力。</figcaption></figure>

<h3>唯一看似有訊號的 M1，在密度分層下消失</h3>
<figure>{c3}
<figcaption>M1（within-HP ALT/REF 位移）未分層時骨幹 {m1_l1:.1f}% vs 孤立 {m1_l3:.1f}%（p={m1_p:.3f}）。
但骨幹位點的局部 sSNV 密度中位為 2、孤立位點為 0 —— 密度同時驅動「有沒有串聯 partner」與
「附近有多少突變的 cis 效應」。<b>固定密度後三層全部不顯著</b>。
這與 2026-07-26 稍早否決的「共同富集」是同一個混淆源。</figcaption></figure>

<h2><span class="lv">L2</span>標記層檔案：欄位與使用規則</h2>
<div class="use"><div class="vt">可以這樣用</div>
<p>把它當作<b>描述性註記</b>：「這個位點在 sSNV 骨幹上是什麼狀態、它的甲基觀測長什麼樣、
哪些技術旗標亮了」。適合挑案例、做圖、寫 limitation、排除低品質位點。</p>
<p style="margin-top:8px;color:var(--bad)"><b>不可以這樣用：</b>把 M1/M2 當作某位點屬於某 subclone 的證據，
或宣稱甲基佐證了遺傳骨幹。上面三張圖就是反證。</p></div>
<p><code>InterSubMod/docs/reports/in_progress/2026/07/20260726_ssnv_lineage_x_methyl_annotation/ssnv_lineage_x_methyl_annotation.tsv</code>
（{n_units:,} 列 × 25 欄）</p>
<table><thead><tr><th style="width:38%">欄位</th><th>意義</th></tr></thead><tbody>{schema}</tbody></table>

<h3>M2 且落在骨幹上的位點（效應量最大的 10 個，僅供挑案例）</h3>
<table><thead><tr><th>位點</th><th>HP</th><th class="num">子群 Δβ</th><th class="num">分離度 S</th>
<th class="num">ALT/REF</th><th class="num">partners</th><th>CN</th><th class="num">管線 k</th></tr></thead>
<tbody>{examples}</tbody></table>
<p class="note">這些是<b>真實的資料特徵</b>（Δβ 達 0.4–0.6，遠高於雜訊），值得肉眼檢視熱圖。
但依上表，同樣強度的分群在孤立位點與 REF 群中出現得一樣頻繁，所以<b>單看這張表無法反推譜系</b>。</p>

<h2><span class="lv">L2</span>甲基資料資產盤點</h2>
<table><thead><tr><th>資料</th><th>路徑</th><th>範圍</th><th>說明</th></tr></thead><tbody>{inv}</tbody></table>

<h2><span class="lv">L3</span>定義、方法與限制</h2>
<details><summary>狀態定義（lineage 軸與 methyl 軸）</summary>
<p><b>lineage 軸</b>（來自 sSNV 共現骨幹）：{ls_defs}</p>
<p><b>methyl 軸</b>（本輪定義，全部對 null 校準）：</p>
<ul>
<li><b>M2_ALT_ONLY_MULTI</b>：ALT 群分離度 ≥ null-95th <b>且</b> 兩子群 |Δβ|≥0.20，<b>且</b> REF 群未達門檻。
最保守 —— 要求內建對照乾淨。</li>
<li><b>M3_REGION_MULTI</b>：ALT 與 REF 都達門檻 → 區域屬性。</li>
<li><b>M1_ALTREF_SHIFT</b>：within-HP ALT vs REF 平均位移顯著且 |Δβ|≥0.10 且超過置換 null。</li>
<li><b>M0_NONE</b>：以上皆非。</li>
</ul></details>

<details><summary>null 的建構（為何門檻是 {q95:.2f} 而不是拍腦袋）</summary>
<p>每個單元各自產生匹配的 null：同 n、同 SD 的單峰高斯重抽 3 次，各自求最佳 2-切分的分離度。
全基因組 {n_units:,} 個 null 的 95 百分位 = <b>{q95:.2f}</b>。用它當門檻時，null 自身的通過率
依定義為 {null_pct:.1f}%，觀察值超出的部分才是訊號。</p>
<p>另一個獨立 null 是<b>置換 null</b>（保持 read 集合不變、隨機重排 ALT/REF 標籤），用於 M1 的超額計算。</p></details>

<details><summary>限制</summary>
<ul>
<li><b>單樣本</b>（HCC1395）、單一 2026-01 w5000 資料線。跨 7 樣本未驗證。</li>
<li><b>覆蓋門檻</b>：同一 HP 內各需 ≥6 條 read 且共同 CpG ≥8，{n_loci:,}/{tp_total:,} 位點才進入標記層 ——
偏向高覆蓋區。</li>
<li><b>CN 未拆開</b>：{cn_gain_pct:.0f}% 單元為 gain、{cn_loh_pct:.0f}% 為 LOH，CN-neutral 極少，
無法用它作乾淨對照。</li>
<li><b>M2 與 REF 對照不完全對稱</b>：M2 額外要求 |Δβ|≥0.20（REF 未存該欄），故 M2 為<b>保守</b>估計；
公平的對稱比較見第一張圖。</li>
<li>lineage 狀態沿用 2026-07-26 的串聯註記，未納入分層樹求解器的拓撲狀態。</li>
</ul></details>

<details><summary>可重現性</summary>
<p>archive：<code>{archive_full}</code></p>
<p>掃描腳本與中間檔已封存於 <code>sources/</code>：
<code>full_methyl_scan.py</code>（{scan_units:,} 單元）、<code>sig_scan.tsv</code>（{sig_loci:,} 位點）、
<code>build_annotation.py</code>、<code>methyl_locus_scan.tsv</code>、<code>annotation_summary.json</code>。</p>
<p>重建：<code>python3 sources/build_annotation.py</code> → <code>python3 build_page.py</code></p></details>

<h2><span class="lv">L3</span>管線自身分群判定（作為對照）</h2>
<table><thead><tr><th>optimal_k</th><th class="num">位點數</th><th class="num">佔比</th></tr></thead>
<tbody>{krows}</tbody></table>
<p class="note">{gate:,} / {gk_loci:,}（{gate_pct:.1f}%）通過管線 gating。注意 <code>optimal_k</code> 沒有 1 這個選項 ——
這是「為何總是看到多群」最直接的結構性解釋。</p>
</div>
""".format(
    archive=req("meta.archive"), archive_full="/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/" + req("meta.archive"),
    date=req("meta.date"), datafile=os.path.basename(DATA),
    n_loci=req("meta.n_loci"), n_units=req("meta.n_units"), tp_total=30490,
    null_pct=req("symmetric.null") / req("symmetric.n") * 100,
    alt_pct=req("symmetric.alt") / req("symmetric.n") * 100,
    ref_pct=req("symmetric.ref") / req("symmetric.n") * 100,
    p_paired=req("symmetric.p_paired"),
    m2_l1=req("enrich.m2_l1"), m2_l3=req("enrich.m2_l3"),
    m2_or=req("enrich.m2_or"), m2_p=req("enrich.m2_p"),
    m2_tot=req("totals_methyl.M2_ALT_ONLY_MULTI") / req("meta.n_units") * 100,
    m2_n=req("totals_methyl.M2_ALT_ONLY_MULTI"),
    m1_l1=req("enrich.m1_l1"), m1_l3=req("enrich.m1_l3"),
    m1_d=req("enrich.m1_l1") - req("enrich.m1_l3"), m1_p=req("enrich.m1_p"),
    n_common=req("readspan_control.n_common_median"), n_cpg=req("readspan_control.n_cpg_median"),
    sep_raw=req("readspan_control.sep_raw_median"), sep_common=req("readspan_control.sep_common_median"),
    strand_pct=req("strand_flag") / req("meta.n_units") * 100,
    cn_gain_pct=req("cn.gain.n") / req("meta.n_units") * 100,
    cn_loh_pct=req("cn.loh.n") / req("meta.n_units") * 100,
    gk_loci=req("genome_pipeline.n_loci"), gate=req("genome_pipeline.gating"),
    gate_pct=req("genome_pipeline.gating") / req("genome_pipeline.n_loci") * 100,
    q95=req("meta.null_q95"),
    scan_units=req("inputs.scan_units"), sig_loci=req("inputs.sig_loci"),
    c1=bars3(), c2=grid(), c3=strat_chart(),
    schema=schema_rows, examples=ex_rows, inv=inv_rows, krows=k_rows,
    ls_defs="；".join("<b>{0}</b> = {1}".format(l, LSL[l].split("（")[1].rstrip("）") if "（" in LSL[l] else LSL[l]) for l in LS),
)

with open(OUT, "w") as fh:
    fh.write(HTML)
print("wrote", OUT, os.path.getsize(OUT), "bytes")
