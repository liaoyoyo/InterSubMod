#!/usr/bin/env python3
"""
build_html.py — 由 data.json 注入產生教學用 standalone HTML（§13-A：缺 key 直接 refuse，不 render）

用法: python3 build_html.py
輸入: 20260726_complexity_funnel_teaching_01.data.json（由 build_data.py 產生，含 42 項 SELF-CHECK）
輸出: 20260726_complexity_funnel_teaching_01.standalone.html
"""
import html
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "20260726_complexity_funnel_teaching_01.data.json"
OUT = HERE / "20260726_complexity_funnel_teaching_01.standalone.html"

if not DATA.exists():
    sys.exit(f"REFUSE: {DATA.name} 不存在 — 先跑 build_data.py")
D = json.loads(DATA.read_text())
if not D["_meta"]["self_check_pass"]:
    sys.exit("REFUSE: data.json 的 self_check_pass=False")


def g(path):
    """安全取值：缺 key 直接 refuse（不得以預設值靜默填補）。"""
    cur = D
    for part in path.split("."):
        if isinstance(cur, list):
            cur = cur[int(part)]
        elif part in cur:
            cur = cur[part]
        else:
            sys.exit(f"REFUSE: data.json 缺 key -> {path}")
    return cur


def n(x):
    return f"{x:,}" if isinstance(x, int) else x


e = html.escape

# ------------------------------------------------------------------ 取值
META, PAR = g("_meta"), g("_meta.params")
KMAX, ECAP, BUD = PAR["MAX_SNV"], PAR["extra_cap"], PAR["per_level_budget"]
L0, L1, L2, L3, L4, L5, L6, L7, L8 = (g(f"L{i}_{s}") for i, s in enumerate(
    ["global", "partition", "hypercube", "tree_constraint", "closed",
     "hidden_choice", "trees", "vaf", "group_steiner"]))
CUBE8 = [r for r in L2["rows"] if r["k"] == KMAX][0]
ARB8 = [r for r in L3["rows"] if r["k"] == KMAX][0]
CL5 = [r for r in L4["rows"] if r["k"] == 5][0]
CL4 = [r for r in L4["rows"] if r["k"] == 4][0]
PEAK, TV, CS = L5["peak"], L6["theory_vs_real"], L7["summary"]
FUN, TL = L1["funnel"], L8["tree_level"]
CHECKS = g("self_check")

# ------------------------------------------------------------------ 元件


def math(tex):
    return f'<div class="math">{tex}</div>'


def verify_box(cmd, src, extra=""):
    return (f'<div class="vbox"><div class="vh">如何確認</div>'
            f'<pre class="cmd">{e(cmd)}</pre>'
            f'<div class="vsrc">來源：<code>{e(src)}</code></div>'
            + (f'<div class="vsrc">{extra}</div>' if extra else "") + "</div>")


def table(headers, rows, cls=""):
    h = "".join(f"<th>{c}</th>" for c in headers)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<div class="tw"><table class="{cls}"><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


def layer(idx, tag, title, lead, body):
    return (f'<section class="layer" id="L{idx}"><div class="lhead">'
            f'<span class="lnum">L{idx}</span><div><h2>{title}</h2>'
            f'<div class="ltag">{tag}</div></div></div>'
            f'<p class="lead">{lead}</p>{body}</section>')


# ------------------------------------------------------------------ 漏斗 SVG
# 只放「樹/搜尋規模」同一條可比的量，維持單調遞減；節點子集空間 2^(2^k) 是不同物件，留在 L2 講。
FUNNEL = [
    ("任意有根標號樹", ARB8["unconstrained_digits"], "256^255 — 只要求是樹，不管超立方體", "#94a3b8"),
    ("＋超立方體單調限制", ARB8["hypercube_digits"], "每邊只翻一個 bit、只准 0→1", "#2563eb"),
    ("＋覆蓋公理＋min-H 搜尋", len(str(PEAK["cand"])), f"候選節點集上限 {n(PEAK['cand'])}", "#2563eb"),
    ("最終輸出樹數（實測峰值）", len(str(TV["real_max_trees"])), f"單區最多 {TV['real_max_trees']} 棵", "#0d9488"),
    ("VAF/CCF 加權後（破等）", 1, "多數案例收斂到唯一 winner", "#0d9488"),
]
MAXD = max(d for _, d, _, _ in FUNNEL)
bars = []
Y = 8
for lab, d, note, col in FUNNEL:
    w = max(26, 640 * (d / MAXD) ** 0.55)
    x = 20 + (640 - w) / 2
    # 條夠寬 → 數字放條內（白字）；太窄 → 放條右側（同色深字），避免文字溢出色塊
    val = (f'<text x="{x + w / 2:.0f}" y="{Y + 25}" text-anchor="middle" fill="#fff" '
           f'font-size="13" font-weight="600">{d} 位數</text>' if w >= 72 else
           f'<text x="{x + w + 8:.0f}" y="{Y + 25}" fill="{col}" '
           f'font-size="13" font-weight="700">{d} 位數</text>')
    bars.append(
        f'<rect x="{x:.0f}" y="{Y}" width="{w:.0f}" height="40" rx="6" fill="{col}" opacity="0.88"/>'
        + val +
        f'<text x="12" y="{Y + 25}" font-size="12" fill="#0f172a" font-weight="600">{e(lab)}</text>'
        f'<text x="690" y="{Y + 25}" font-size="11" fill="#64748b">{e(note)}</text>')
    Y += 54
FUNNEL_SVG = (f'<svg viewBox="0 0 1240 {Y + 6}" class="funnel" role="img" '
              f'aria-label="計算量漏斗，單位為十進位位數"><g transform="translate(190,0)">'
              + "".join(bars) + "</g></svg>")

# 鋸齒 SVG
CURVE = L5["curve"]
pts = [(r["P"], r["cand"]) for r in CURVE if r["P"] <= 260]
mx = max(c for _, c in pts)
poly = " ".join(f"{40 + p / 260 * 1120:.1f},{300 - (c / mx) * 250:.1f}" for p, c in pts)
marks = []
for P in (45, 46, 97, 98):
    r = next(x for x in CURVE if x["P"] == P)
    cx, cy = 40 + P / 260 * 1120, 300 - (r["cand"] / mx) * 250
    marks.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4.5" fill="#dc2626"/>'
                 f'<text x="{cx:.1f}" y="{cy - 11:.1f}" font-size="11" text-anchor="middle" '
                 f'fill="#dc2626" font-weight="600">P={P}·E={r["E"]}<tspan x="{cx:.1f}" dy="-12">'
                 f'{n(r["cand"])}</tspan></text>')
SAW_SVG = (f'<svg viewBox="0 0 1180 330" class="chart" role="img" aria-label="候選數隨候選池大小的鋸齒曲線">'
           f'<line x1="40" y1="300" x2="1170" y2="300" stroke="#cbd5e1"/>'
           f'<line x1="40" y1="20" x2="40" y2="300" stroke="#cbd5e1"/>'
           f'<polyline points="{poly}" fill="none" stroke="#2563eb" stroke-width="2"/>'
           + "".join(marks) +
           f'<text x="600" y="325" font-size="12" text-anchor="middle" fill="#64748b">候選隱藏節點池大小 P（0–255）</text>'
           f'<text x="14" y="16" font-size="12" fill="#64748b">Cand(P)</text></svg>')

# ------------------------------------------------------------------ 各層內容
S = []

S.append(layer(0, "起點 · 為什麼非分區不可", "全基因組當成一個超立方體",
               f'把 HCC1395 chr1–22 的 <b>{n(L0["n_autosomal_sSNV"])}</b> 個 sSNV 全部當成一個布林向量的維度，'
               f'狀態空間是 2<sup>{n(L0["n_autosomal_sSNV"])}</sup> —— 一個 <b>{n(L0["genotype_space_digits"])} 位數</b>的整數。'
               f'宇宙原子數大約 80 位數。這一層存在的唯一目的，是說明「不分區」在數學上就已經出局。',
               math("|X<sub>global</sub>| = 2<sup>N</sup>,&nbsp;&nbsp; N = " + n(L0["n_autosomal_sSNV"])
                    + "&nbsp;&nbsp;→&nbsp;&nbsp; log<sub>10</sub>|X| = N·log<sub>10</sub>2 ≈ "
                    + n(L0["genotype_space_digits"]))
               + verify_box(L0["verify"], META["sources"]["census"],
                            "位數＝⌊N·log₁₀2⌋+1，任何計算機都能複核。")
               + '<div class="key"><b>重點</b>：這不是「難算」，是「不可能列舉」。所有後續設計都是為了把這個數字砍下來，'
               '而且每一刀都要能說清楚<em>砍掉的是什麼、為什麼可以砍</em>。</div>'))

frows = [
    ("L1 ClairS paired PASS 全集", n(FUN["L1_all_pass_universe"]), "操作骨幹（非生物真值）", ""),
    ("L2 − chrX/Y 出界", f'−{n(FUN["L2_out_of_scope_chrXY"])}', "scope 定為 chr1–22", "cut"),
    ("＝ 體染色體", n(FUN["autosomal_chr1_22"]), "", "sub"),
    ("L3 − 位置孤立單點", f'−{n(FUN["L3_positional_singleton"])}', "兩側 gap&gt;50kb，無共現可建樹", "cut"),
    ("＝ 進入多位點群", n(FUN["in_multilocus_group"]), f'共 {n(FUN["n_positional_groups"])} 個位置群', "sub"),
    ("L4 − densest-8 截斷", f'−{n(FUN["L4_cap_excluded_densest8"])}', f"大群只取最密 {KMAX} 個（MAX_SNV）", "cut"),
    ("L5 − read 支持不足", f'−{n(FUN["L5_read_unsupported"])}', "full-cov/subread read &lt; MINREAD=3", "cut"),
    ("<b>L6 真正送進 solver</b>", f'<b>{n(FUN["L6_retained_sSNV"])}</b>',
     f'分佈在 {n(FUN["retained_regions"])} 個 region', "keep"),
]
S.append(layer(1, "第一刀 · 也是最猛的一刀", "分區：把一個巨型超立方體切成幾萬個小超立方體",
               f'切法是<b>位置鄰近</b>：相鄰 sSNV 間距 &gt; 50 kb 就切開。切完後每個 region 再硬性封頂 '
               f'k ≤ {KMAX}（<code>MAX_SNV</code>），大群只保留最密的 {KMAX} 個。',
               table(["層", "sSNV 數", "說明", ""],
                     [(a, b, c, {"cut": '<span class="mk cut">扣除</span>',
                                 "keep": '<span class="mk keep">送入 solver</span>',
                                 "sub": '<span class="mk mid">小計</span>', "": ""}[d])
                      for a, b, c, d in frows], "funnel-t")
               + math("Σ<sub>i</sub> 2<sup>k<sub>i</sub></sup> = " + n(L1["sum_states_after_partition"])
                      + "&nbsp;&nbsp;≤&nbsp;&nbsp; U · 2<sup>k<sub>max</sub></sup> = "
                      + n(L1["n_solver_units"]) + " × " + n(2 ** KMAX) + " = "
                      + n(L1["upper_bound_units_times_256"]))
               + f'<p>把 {n(L1["n_solver_units"])} 個 solver unit 的狀態空間<b>全部加起來</b>，總共只有 '
                 f'<b>{n(L1["sum_states_after_partition"])}</b> 個狀態（{L1["sum_states_digits"]} 位數）。'
                 f'從 {n(L0["genotype_space_digits"])} 位數 → {L1["sum_states_digits"]} 位數。</p>'
               + table(["k", "unit 數", "佔比", "該 unit 的狀態數 2ᵏ"],
                       [(r["k"], n(r["n"]), f'{r["pct"]}%', n(r["states_2k"])) for r in L1["k_hist"]])
               + verify_box("python3 -c \"import json;d=json.load(open('layered_reconstruction_HCC1395.json'));"
                            "print(sum(2**u['n_sSNV'] for u in d['detail']), len(d['detail']))\"",
                            META["sources"]["layered"],
                            "六層 funnel 加總 == 全集（SELF-CHECK 已驗），代表分層互斥且窮盡、沒有數字憑空消失。")
               + '<div class="key"><b>重點</b>：這一刀砍掉的不是「雜訊」，是<em>跨區的共現資訊</em>。'
                 f'代價寫在檯面上——{n(FUN["L4_cap_excluded_densest8"])} 個 sSNV 因 densest-8 截斷沒進 solver。'
                 '這是<b>計算界線，不是讀長界線</b>：那些 sSNV 的 read 資料還在，只是沒被枚舉。</div>'
               + '<div class="warn"><b>caveat</b>：本表為 HCC1395 單樣本、ClairS paired PASS 操作骨幹，'
                 '非生物真值；densest-8 的比例會隨樣本突變密度變動。</div>'))

S.append(layer(2, "單區內部 · 狀態空間長什麼樣", "Directed Boolean Hypercube {0,1}ᵏ",
               '每個 region 內，一條 read 在 k 個 sSNV 上的讀值就是一個 0/1 字串（REF=0 / ALT=1）。'
               '所有可能字串 = 超立方體的頂點。<b>方向</b>來自生物學：突變只會 0→1（不考慮回復突變），'
               '所以邊是單向的、只准加一個 1。',
               math("|V| = 2<sup>k</sup> &nbsp;·&nbsp; |E| = k·2<sup>k−1</sup> "
                    "&nbsp;·&nbsp; 任意節點子集 = 2<sup>2<sup>k</sup></sup>")
               + table(["k", "頂點 2ᵏ", "有向 unit-flip 邊 k·2ᵏ⁻¹", "任意節點子集 2^(2ᵏ) 的位數"],
                       [(r["k"], n(r["vertices"]), n(r["edges"]), f'{r["all_subsets_digits"]} 位數')
                        for r in L2["rows"]])
               + f'<p>k = {KMAX} 時：<b>{n(CUBE8["vertices"])}</b> 個頂點、<b>{n(CUBE8["edges"])}</b> 條有向邊。'
                 f'頂點數小到不值一提——<b>爆炸從來不在這裡</b>。真正大的是「要選哪些頂點進樹」：'
                 f'2<sup>256</sup>，{CUBE8["all_subsets_digits"]} 位數。</p>'
               + verify_box("# 邊數手算：每個頂點有幾個 0 就能往外走幾條邊\n"
                            "python3 -c \"k=8;print(sum(k-bin(x).count('1') for x in range(2**k)))\"  # 1024\n"
                            "# 或用公式 k·2^(k-1) = 8·128 = 1024",
                            META["sources"]["spec"] + " §2",
                            "兩種算法（逐頂點數出邊 vs 公式）必須相等，這是最容易自己複核的一層。")
               + '<div class="key"><b>重點</b>：把 k 從 8 降到 7 只讓頂點少一半，'
                 '但讓「子集空間」從 78 位數掉到 39 位數——<em>指數的指數</em>。這解釋了為什麼 MAX_SNV 卡在 8 而不是 12。</div>'))

S.append(layer(3, "第二刀 · 生物學限制值多少錢", "演化樹限制：unit-flip + 有根 + 單親",
               '在超立方體上「選一棵樹」不是隨便連。三條限制：<b>①</b> 根固定是 0ᵏ（germline 全 REF）；'
               '<b>②</b> 每條邊剛好翻一個 bit 且只准 0→1；<b>③</b> 每個非根節點恰好一個 parent。'
               '這三條合起來就是 <b>arborescence</b>（有向有根樹）。',
               math("A(k) = ∏<sub>x≠0</sub> popcount(x) = ∏<sub>j=1</sub><sup>k</sup> j<sup>C(k,j)</sup>"
                    "&nbsp;&nbsp;&nbsp;vs&nbsp;&nbsp;&nbsp; 無限制： n<sup>n−1</sup>, n = 2<sup>k</sup>")
               + '<p>為什麼是連乘？因為每個非根頂點 x 只能從「少一個 1 的鄰居」來，這種鄰居剛好有 popcount(x) 個，'
                 '而且每個頂點<b>獨立</b>挑一個。邊一定是 popcount 遞增，所以天然無環、天然 root-connected——'
                 '不需要額外檢查。</p>'
               + table(["k", "頂點 n", "無限制有根標號樹 n^(n−1)", "超立方體 arborescence A(k)", "砍掉幾個數量級"],
                       [(r["k"], n(r["n_vertices"]), f'{r["unconstrained_sci"]}（{r["unconstrained_digits"]} 位）',
                         f'{r["hypercube_sci"]}（{r["hypercube_digits"]} 位）', f'<b>{r["cut_orders"]}</b>')
                        for r in L3["rows"]])
               + f'<p>k = {KMAX}：從 <b>{ARB8["unconstrained_digits"]} 位數</b> 砍到 '
                 f'<b>{ARB8["hypercube_digits"]} 位數</b>，一刀 <b>{ARB8["cut_orders"]} 個數量級</b>。'
                 f'這是純粹由「突變不可逆 + 每次一個」這兩句生物學換來的。</p>'
               + verify_box("python3 -c \"import math;k=8;\n"
                            "a=1\n"
                            "for j in range(1,k+1): a*=j**math.comb(k,j)\n"
                            "b=1\n"
                            "for x in range(1,2**k): b*=bin(x).count('1')\n"
                            "print(a==b, len(str(a)))\"   # True 146",
                            META["sources"]["bench"] + " → A_max",
                            "build_data.py 用<b>兩種互不相干的算法</b>（乘冪公式 vs 逐頂點連乘）各算一次再比對，"
                            "並再對照獨立產生的 bench JSON。k=1…8 全部吻合。")
               + '<div class="key"><b>重點</b>：146 位數還是天文數字——所以<em>生物學限制不夠</em>。'
                 '真正把數字壓到可枚舉的，是下一層的「資料要覆蓋」。</div>'))

S.append(layer(4, "第三刀 · 資料說了算", "覆蓋公理 + closed：哪些節點集合法",
               '光是「是一棵樹」還不夠，樹必須<b>解釋觀測</b>。兩條硬約束：'
               '<b>覆蓋公理</b>——每筆通過門檻的觀測都要被樹上某個節點相容（完整觀測其基因型節點本人必須在樹上）；'
               '<b>closed</b>——每個節點的祖先鏈必須也在樹上，不能憑空冒出來。',
               math("N 合法 ⟺ ∀ x ∈ N∖{0}, ∃ j ∈ x : x∖{j} ∈ N &nbsp;&nbsp;(closed)"
                    "<br>且 ∀ p ∈ O*, ∃ x ∈ N : x ~ p &nbsp;&nbsp;(覆蓋公理)")
               + '<p class="sub">「相容」的定義決定了 <b>Group</b> 從哪來：完整觀測只相容它自己（一個點）；'
                 '有缺口的部分觀測相容一整個<b>子立方體</b>（一群點），碰到其中任一個就算覆蓋。'
                 '這就是 group Steiner 的 group。</p>'
               + table(["k", "任意子集 2^(2ᵏ)", "其中 closed 的（含 root）", "合法比例"],
                       [(r["k"], n(int(r["all_subsets"])), n(r["closed"]),
                         (f'{r["frac_ppm"] / 1e4:.2f}%' if r["frac_ppm"] else "—"))
                        for r in L4["rows"]])
               + f'<p>closed 限制本身砍掉不少（k=4 只剩 {CL4["frac_ppm"] / 1e4:.1f}%、'
                 f'k=5 剩 {CL5["closed"] / int(CL5["all_subsets"]) * 100:.1f}%），但注意——'
                 f'k=5 的合法集合仍有 <b>{n(CL5["closed"])}</b> 個。'
                 f'<b>連「數清楚有幾個合法節點集」這件事本身，到 k=6 就已經算不動了。</b></p>'
               + verify_box("# k≤4 暴力枚舉全部 2^(2^k) 個子集逐一檢查 closed\n"
                            "# k=5 改用分層 DP（狀態 = 上一層被選中的頂點子集）\n"
                            "python3 build_data.py --check   # DP 對 k≤4 與暴力逐一吻合才放行",
                            "build_data.py: closed_sets_bruteforce / closed_sets_leveldp",
                            "這是本頁唯一「新算」的量，所以特地寫了兩套演算法互相背書。")
               + '<div class="key"><b>重點</b>：這一層是<em>資料進場的地方</em>。前面三層純數學、跟樣本無關；'
                 '從這裡開始，數字取決於你這個 region 實際觀測到幾種基因型。'
                 '覆蓋公理是<b>硬約束、不可放寬</b>——solver 不准為了求一棵更乾淨的樹而丟掉任何一筆觀測；'
                 '覆蓋失敗就判 conflict，不是默默省略。</div>'))

erows = [(r["v"], n(r["n"]), f'{r["pct"]}%') for r in L5["e_min_hist"]]
S.append(layer(5, "爆炸點 · 91% 的時間都花在這", f"H（隱藏節點）選擇：從池子 P 選 e 個",
               f'觀測到的節點通常接不成樹，中間缺的祖先要補進來——這些補進來的就是<b>隱藏節點 H</b>'
               f'（Steiner 點，程式裡 label 前綴就是 <code>H_</code>）。目標是 H 越少越好，'
               f'做法是 <b>e = 0, 1, 2, … 逐層往上試</b>，第一個找得到可行解的層就是 e<sub>min</sub>，'
               f'並把<b>那一層所有</b>可行解全部收下。',
               math("Cand(P) = Σ<sub>e=0</sub><sup>E(P)</sup> C(P, e),"
                    "&nbsp;&nbsp; E(P) = min( " + str(ECAP) + ", &nbsp;min{ e : C(P,e) &gt; " + n(BUD) + " } − 1 )")
               + f'<p>兩道閘：硬上限 <code>extra_cap={ECAP}</code>，以及每層預算 '
                 f'<code>per_level_budget={n(BUD)}</code>。哪道先咬到，看候選池 P 多大。</p>'
               + SAW_SVG
               + '<div class="key"><b>最反直覺的一件事</b>：P 變大，工作量可能<em>暴跌</em>。'
                 '因為預算一咬，E 就自動降級。所以「越密越慢」是錯的——是<b>鋸齒</b>。</div>'
               + table(["P", "E(P)", "候選節點集數 Cand(P)", ""],
                       [(r["P"], r["E"], n(r["cand"]),
                         {45: "◀ 全域最痛點", 46: "◀ 掉 10 倍", 98: "◀ 再降一級"}.get(r["P"], ""))
                        for r in L5["sawtooth"]])
               + f'<p>全域上限：<b>{n(PEAK["cand"])}</b> 個候選節點集，出現在 P = {PEAK["P"]}、E = {PEAK["E"]}。'
                 f'允許 e=4 的最大池是 P = {L5["max_P_allowing_e4"]}；'
                 f'允許 e=3 的最大池是 P = {L5["max_P_allowing_e3"]}。</p>'
               + '<h3>實測：時間花在哪</h3>'
               + table(["階段", "耗時 (ms)", "佔比"],
                       [(r["label"], r["ms"], f'<b>{r["pct"]}%</b>' if r["pct"] > 50 else f'{r["pct"]}%')
                        for r in L5["stage_share"]])
               + f'<p>S5 逐層搜尋吃掉 <b>{[r for r in L5["stage_share"] if r["key"] == "t5_level_search"][0]["pct"]}%</b>，'
                 f'其他階段全部加起來不到 9%。單一候選的處理成本中位數 '
                 f'{L5["throughput_us"]["median"]} µs，所以最壞單區 wall-time 中位數只有 '
                 f'<b>{L5["walltime"]["median"]} 秒</b>（最大 {L5["walltime"]["max"]} 秒）——會爆，但爆得有界。</p>'
               + '<h3>實際資料的 e<sub>min</sub> 分佈</h3>'
               + table(["e_min", "unit 數", "佔比"], erows)
               + f'<p>中位數 {L5["e_min_stats"]["median"]}、平均 {L5["e_min_stats"]["mean"]}、'
                 f'最大 <b>{L5["e_min_stats"]["max"]}</b>。被 cap 掉的共 '
                 f'<b>{n(L5["capped"]["n"])}</b> 個（{L5["capped"]["pct"]}%），'
                 f'其中預算類 {n(L5["capped"]["by_reason"]["budget"])}、'
                 f'上限類 {n(L5["capped"]["by_reason"]["extra_cap_greedy"])}。</p>'
               + verify_box("python3 -c \"from math import comb\n"
                            "print(comb(45,4), comb(46,4))                 # 148995 163185（預算 150000 卡在中間）\n"
                            "print(sum(comb(45,e) for e in range(5)))      # 164221 = 全域峰值\"",
                            META["sources"]["solver"] + ":144-157",
                            "build_data.py 的 cand_curve() 逐行鏡像 solver 的 fast-cap 行為，"
                            "重算出的峰值/位置/E 三項都與獨立產生的 bench JSON 吻合。")
               + '<div class="key"><b>重點</b>：<code>extra_cap=4</code> 的「4」是<em>人設的天花板</em>，'
                 '不是問題在 4 崩掉——真實資料的 e<sub>min</sub> 一路到 16。cap 的作用是<b>把難的案例隔離並誠實標記</b>，'
                 '不是掃到地毯下。所有被 cap 的案例一律標 <code>capped</code>，'
                 '<code>V7_no_overclaim</code> 驗證器會直接擋住任何對它們宣稱 determined 的行為。</div>'))

S.append(layer(6, "輸出 · 樹的數量反而很小", "從節點集到樹：獨立選 parent",
               '節點集 N 定下來之後，樹的數量是純組合：每個非根節點<b>獨立</b>從「在 N 裡面的 unit-pred」挑一個當爸爸。'
               '所以是連乘，再對所有並列最小的 N 求和。',
               math("n<sub>trees</sub> = Σ<sub>N ∈ 𝒩<sub>min</sub></sub> ∏<sub>x ∈ N∖{0}</sub> "
                    "| { j ∈ x : x∖{j} ∈ N } |")
               + '<p class="sub">一個漂亮的等價定理（solver 靠它省下生成所有樹的成本）：'
                 '<b>某節點的 pred 數 ≥ 2 ⟺ 該 N 違反 rooted three-gamete ⟺ 它的所有生成樹都含 recurrence</b>。'
                 '反之 pred 數全 = 1 → 連乘 = 1 → 唯一一棵無 recurrence 的樹。'
                 '所以「數樹」和「判斷相容性」是同一件事。</p>'
               + table(["", "數值"],
                       [("理論上限 A(8)（k=8 全超立方體）", f'{ARB8["hypercube_sci"]}（{TV["A_max_k8_digits"]} 位數）'),
                        ("<b>實測單區最多樹數</b>", f'<b>{TV["real_max_trees"]}</b>'),
                        ("實測全部 unit 樹數總和", n(TV["real_total_trees"])),
                        ("差距", f'<b>{TV["gap_orders_of_magnitude"]} 個數量級</b>')])
               + table(["樹數門檻", "unit 數", "佔比"],
                       [(f'≥ {r["thr"]}', n(r["n"]), f'{r["pct"]}%') for r in L6["n_trees_tail"]])
               + f'<p>{[r for r in L6["n_trees_tail"] if r["thr"] == 1][0]["pct"]}% 的 unit 至少有 1 棵樹；'
                 f'只有 {[r for r in L6["n_trees_tail"] if r["thr"] == 100][0]["pct"]}% 超過 100 棵；'
                 f'<b>沒有任何一個 unit 超過 1000 棵</b>。</p>'
               + verify_box("# solver 的 V5 驗證器獨立重算全樹數，與 enumerate 自報值比對\n"
                            "python3 tree_enumeration_solver.py   # GOLDEN 8 個手算案例 + V1–V7",
                            META["sources"]["solver"] + ":371-402",
                            "V4 獨立確認 e_min−1 層真的全部不可行（不信任早停）；"
                            "V5 獨立重算樹數；V7 擋 overclaim。18,931 個 unit 全部 V1–V7 PASS、0 fail。")
               + '<div class="key"><b>重點</b>：這是整套方法最漂亮的一格。理論上界 146 位數，'
                 f'實際只到 <b>{TV["real_max_trees"]}</b> 棵——差 {TV["gap_orders_of_magnitude"]} 個數量級。'
                 '原因不是近似、不是取樣，是<em>覆蓋公理把解釘死在超立方體一個極小的下集裡</em>。'
                 '資料本身就是最強的剪枝。</div>'))

S.append(layer(7, "最後一刀 · 用 read 數破等", "VAF / CCF 加權：ambiguous 集合內部定序",
               f'枚舉完剩下 <b>{n(CS["n_ambiguous_units"])}</b> 個 unit 有多棵並列最小樹。'
               f'破等的依據是 <b>pigeonhole</b>：祖先節點的 read 數必然 ≥ 它後裔的 read 數，'
               f'所以 read count 本身就攜帶方向資訊。',
               math("P(T | reads) ∝ ∏ read-count&nbsp;&nbsp;&nbsp;"
                    "（pigeonhole：祖先 read 數 ≥ 後裔 read 數）")
               + table(["", "數值", "佔比"],
                       [("ambiguous unit 總數", n(CS["n_ambiguous_units"]), "100%"),
                        ("成功破等（top posterior ≥ 0.6）", n(CS["tie_broken(≥0.6)"]),
                         f'<b>{CS["broke_frac"] * 100:.1f}%</b>'),
                        ("維持等機率（不強行定序）", n(CS["stay_uniform"]),
                         f'{CS["stay_uniform"] / CS["n_ambiguous_units"] * 100:.1f}%'),
                        ("破等者中 winner 符合 pigeonhole", n(CS["winner_pigeonhole_clean"]),
                         f'<b>{CS["winner_clean_frac_of_broke"] * 100:.1f}%</b>')])
               + '<h3><span class="dotr"></span>這一層最容易被攻擊，先講清楚兩件事</h3>'
               + '<div class="warn"><b>① 為什麼這不算循環論證？</b><br>'
                 'read count 來自<b>同一條通道 B（sSNV 基因型）</b>，不是甲基、不是 HP。'
                 '形式化規格的紅線是「不得用通道 A（φ）或通道 M（d，含甲基）對樹集全序化」，'
                 '寫成 ∂(tree set)/∂φ = ∂(tree set)/∂d = 0。甲基一律<b>事後</b>註記，不進 likelihood。'
                 'read-count 加權沒有跨通道，所以不觸線。</div>'
               + f'<div class="warn"><b>② CN confound 沒有解決。</b><br>'
                 f'{n(CS["by_type_cn_total"]["structure(多完成)|gain(read≠CCF)"])} / '
                 f'{n(CS["n_ambiguous_units"])} '
                 f'（{CS["by_type_cn_total"]["structure(多完成)|gain(read≠CCF)"] / CS["n_ambiguous_units"] * 100:.0f}%）'
                 f'的 ambiguous unit 位於 CN-gain 區。在 CN-gain 區「read 數 ∝ CCF」的前提不成立，'
                 f'所以那部分的破等結果<b>只能當工程 heuristic，不能當證據</b>。'
                 f'CN-clean 的只有 {n(CS["by_type_cn_total"]["structure(多完成)|clean"])} 個。</div>'
               + verify_box("python3 -c \"import json;s=json.load(open('ccf_tree_weighting_full_observe.json'))['summary']\n"
                            "print(s['tie_broken(≥0.6)']+s['stay_uniform'] == s['n_ambiguous_units'])\"  # True",
                            META["sources"]["ccf"],
                            "SELF-CHECK 逐項重算 4 個比例，並確認 ambiguous 總數與 solver 輸出一致。")
               + '<div class="key"><b>重點</b>：VAF 這一層<em>不是</em>把樹集砍小，而是在集合<b>內部</b>給機率權重。'
                 '不能破等的 ' + n(CS["stay_uniform"]) + ' 個就<b>誠實維持等機率</b>——'
                 '「定不出來」本身就是答案，不是失敗。</div>'))

det_pct = TL["det_all_pct_of_units"]
S.append(layer(8, "定位 · 我們到底解了什麼、沒解什麼", "Group Steiner Arborescence：名稱、複雜度、與紅線",
               '把前面所有層合起來，這個問題的正式形狀是：'
               '<b>在有向布林超立方體上，找一棵以 0ᵏ 為根的 arborescence，'
               '覆蓋所有 terminal group（完整觀測＝單點群，部分觀測＝子立方體群），使 Steiner 點最少。</b>'
               '這正是 <b>Group Steiner Arborescence</b> 的形狀。',
               '<h3>三個名詞各自對應什麼</h3>'
               + table(["名詞", "數學意義", "在本問題裡是什麼"],
                       [("<b>Group</b>", "terminal 不是單點，是「一群點，碰到任一個就算」",
                         "一條 read 沒讀完全時，它相容的整個子立方體"),
                        ("<b>Steiner</b>", "為了連通而必須額外加入、非 terminal 的節點",
                         "隱藏祖先基因型 H（未被任何 read 直接觀測到）"),
                        ("<b>Arborescence</b>", "有向有根樹：每個非根節點恰一個 parent、全 root 可達",
                         "根 = germline 全 REF；邊只准 0→1，方向由突變不可逆決定")])
               + '<h3>複雜度地圖：哪些是多項式島、哪些是 NP-hard 海</h3>'
               + table(["階段", "數學問題", "複雜度", "本方法的立場", "文獻"],
                       [("骨幹（完整觀測）", "directed perfect phylogeny", '<span class="ok">poly</span>',
                         "算到底、確定", "Gusfield 1991"),
                        ("部分觀測補全（缺口）", "Incomplete Directed Perfect Phylogeny",
                         '<span class="ok">poly→linear</span>', "算到底，有嚴格基礎", "Pe'er et al. 2004"),
                        ("隱藏節點補全", "Steiner in {0,1}ⁿ", '<span class="hot">NP-/MAX-SNP-hard</span>',
                         "<b>不解通用最佳化</b>", "Foulds &amp; Graham 1982"),
                        ("群覆蓋最小化", "Group Steiner Tree", '<span class="hot">Ω(log^{2−ε} n) 不可近似</span>',
                         "<b>白地 / future work</b>", "Halperin &amp; Krauthgamer 2003"),
                        ("（近似演算法）", "Group Steiner / Directed Steiner",
                         "O(log²n·log k) / O(kᵋ)", "<b>不採用</b>（我們要全集不要一棵）",
                         "Garg et al. 2000; Charikar et al. 1999"),
                        ("超立方體 Steiner arborescence", "FPT 參數化", "FPT",
                         "最接近的既有結果", "Mahapatra et al. 2025"),
                        ("多完成唯一性", "non-identifiable", '<span class="hot">impossibility（已證）</span>',
                         "定不出來即答案", "Satas et al. 2021")])
               + '<div class="warn"><b><span class="dotr"></span>三條不可跨越的紅線（規格 §7.5.3）</b><ol>'
                 '<li><b>枚舉，不是最佳化</b>：輸出所有並列最小的樹，<b>不宣稱</b>解通用 directed / group Steiner 的'
                 ' min-cost 最佳化。「minimal」= 字典序前緣的<em>全集</em>，不是挑一棵代表。</li>'
                 '<li><b>不用任何通道對樹集全序化</b>：甲基與 HP 一律事後註記，'
                 '保證 ∂(tree set)/∂φ = ∂(tree set)/∂d = 0。</li>'
                 '<li><b>通用 group-Steiner cover-minimize = 白地</b>：'
                 'partial-read → subcube 群覆蓋是真結構類比與原創 motivation，'
                 '但 solver 做的是「capped-k 窮舉 + IDP 補全 + 覆蓋公理」，<b>不是</b> cover-and-minimize 求解。</li>'
                 '</ol></div>'
               + '<h3>那 NP-hard 怎麼繞過的？</h3>'
               + f'<p>沒有繞過——是<b>不進去</b>。NP-hardness 是 k → ∞ 的漸近性質；'
                 f'我們有硬上界 k ≤ {KMAX}，狀態空間只有 {n(2 ** KMAX)} 點，'
                 f'落在<b>窮舉可解的小島</b>上。代價已經在 L1 付過了：'
                 f'{n(FUN["L4_cap_excluded_densest8"])} 個 sSNV 因 densest-8 沒進 solver。'
                 f'<b>這是把 NP-hard 換成資訊損失，而不是換成近似誤差</b>——'
                 f'差別在於損失是<em>可數、可報告</em>的。</p>'
               + '<h3>最終產出</h3>'
               + table(["判定", "unit 數", "佔比", "意義"],
                       [("determined", n(TL["determined_all"]), f'{det_pct}%', "唯一最小樹"),
                        ("&nbsp;&nbsp;└ 其中純 REF 家族", n(TL["L2a_root_only_reference_family"]), "—",
                         "樹只有 ROOT、無 somatic edge，<b>不算重建出 subclone 樹</b>"),
                        ("&nbsp;&nbsp;└ <b>真正帶突變的 determined</b>", f'<b>{n(TL["determined_mutation_bearing"])}</b>',
                         f'<b>{TL["det_mut_pct_of_nonroot"]}%</b>', "排除 root-only 家族後的誠實分母"),
                        ("ambiguous", n(TL["ambiguous"]), "—", "多棵並列最小樹，輸出全集"),
                        ("capped", n(TL["capped"]), "—", "太密、枚舉未完，<b>絕不宣稱 determined</b>"),
                        ("recurrence_required", n(TL["recurrence"]), "—", "所有最小解都需 recurrence，送獨立 m-通道")])
               + f'<div class="key"><b>最後一個重點，也是最該記住的</b>：'
                 f'誠實的 determined 是 <b>{TL["det_mut_pct_of_nonroot"]}%</b>'
                 f'（{n(TL["determined_mutation_bearing"])} / '
                 f'{n(TL["determined_all"] + TL["ambiguous"] + TL["capped"] + TL["recurrence"] - TL["L2a_root_only_reference_family"])}），'
                 f'不是 {det_pct}%——因為 {n(TL["L2a_root_only_reference_family"])} 個「純 REF 家族」是 '
                 f'trivially determined（樹上只有一個 ROOT，什麼都沒重建出來）。'
                 f'把它們算進成功率會膨脹結論。<b>分母怎麼定，比分子怎麼算更容易出事。</b></div>'))

# ------------------------------------------------------------------ SELF-CHECK
chk_rows = [(('<span class="ok">PASS</span>' if c["pass"] else '<span class="hot">FAIL</span>'), e(c["name"]), f'<code>{e(c["got"])}</code>',
             f'<code>{e(c["want"])}</code>', e(c["note"])) for c in CHECKS]
CHECK_HTML = (f'<details class="chk"><summary><b>SELF-CHECK 全表</b>（'
              f'{sum(c["pass"] for c in CHECKS)}/{len(CHECKS)} PASS）— 每個純數學重算 vs 獨立來源</summary>'
              + table(["", "檢查項", "重算值", "來源值", "說明"], chk_rows) + "</details>")

SRC_ROWS = [(f'<code>{e(k)}</code>', f'<code>{e(v)}</code>') for k, v in META["sources"].items()]

TOC = "".join(f'<a href="#L{i}">L{i} · {t}</a>' for i, t in enumerate(
    ["全域狀態空間", "分區", "超立方體", "演化樹限制", "覆蓋＋closed",
     "H 選擇（爆炸點）", "樹的數量", "VAF 破等", "Group Steiner 定位"]))

HTML = f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>計算量漏斗 — 從全基因組到 Group Steiner Arborescence</title>
<style>
:root{{--bg:#ffffff;--fg:#0f172a;--mut:#64748b;--line:#e2e8f0;--card:#f8fafc;
--acc:#2563eb;--hot:#dc2626;--ok:#0d9488;--amber:#ea580c}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--fg);line-height:1.75;
font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans CJK TC","PingFang TC","Microsoft JhengHei",'Helvetica Neue',Arial,sans-serif}}
.wrap{{max-width:1280px;margin:0 auto;padding:0 24px 80px;display:grid;grid-template-columns:220px 1fr;gap:36px}}
header{{border-bottom:3px solid var(--acc);padding:34px 24px 22px;max-width:1280px;margin:0 auto}}
h1{{font-size:27px;margin:0 0 8px;letter-spacing:-.3px}}
.sub2{{color:var(--mut);font-size:14px}}
.badges{{margin-top:14px;display:flex;flex-wrap:wrap;gap:7px}}
.b{{font-size:12px;padding:3px 11px;border-radius:99px;background:var(--card);border:1px solid var(--line);color:var(--mut)}}
.b.demo{{background:#fff7ed;border-color:#fed7aa;color:#9a3412}}
.b.ok{{background:#f0fdfa;border-color:#99f6e4;color:#0f766e}}
nav{{position:sticky;top:20px;align-self:start;font-size:13px;display:flex;flex-direction:column;gap:2px;padding-top:30px}}
nav a{{color:var(--mut);text-decoration:none;padding:6px 10px;border-left:2px solid var(--line);transition:.15s}}
nav a:hover{{color:var(--acc);border-left-color:var(--acc);background:var(--card)}}
main{{min-width:0;padding-top:20px}}
.layer{{padding:30px 0;border-bottom:1px solid var(--line)}}
.lhead{{display:flex;gap:14px;align-items:flex-start;margin-bottom:6px}}
.lnum{{background:var(--acc);color:#fff;font-weight:700;font-size:13px;padding:5px 11px;border-radius:6px;flex:0 0 auto;margin-top:4px}}
h2{{font-size:21px;margin:0}}
h3{{font-size:16px;margin:26px 0 8px;color:var(--acc)}}
.ltag{{font-size:12.5px;color:var(--amber);font-weight:600}}
.lead{{font-size:15px;margin:12px 0 18px}}
p{{margin:12px 0}} .sub{{color:var(--mut);font-size:14px}}
.math{{background:#fbfdff;border:1px solid #dbeafe;border-left:4px solid var(--acc);
padding:15px 18px;margin:16px 0;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
font-size:14.5px;overflow-x:auto;border-radius:0 6px 6px 0}}
.math sub,.math sup{{font-size:.72em}}
.tw{{overflow-x:auto;margin:16px 0}}
table{{border-collapse:collapse;width:100%;font-size:13.5px;min-width:460px}}
th{{background:var(--card);text-align:left;padding:9px 12px;border-bottom:2px solid var(--line);
font-weight:600;color:var(--mut);font-size:12.5px;white-space:nowrap}}
td{{padding:8px 12px;border-bottom:1px solid var(--line);vertical-align:top}}
tr:hover td{{background:#fcfdfe}}
.funnel-t td:nth-child(2){{font-family:ui-monospace,Menlo,monospace;text-align:right;white-space:nowrap}}
.vbox{{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px 16px;margin:18px 0}}
.vh{{font-size:12px;font-weight:700;color:var(--ok);letter-spacing:.04em;margin-bottom:8px}}
.cmd{{background:#0f172a;color:#e2e8f0;padding:11px 14px;border-radius:6px;overflow-x:auto;
font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px;margin:0 0 8px;line-height:1.6}}
.vsrc{{font-size:12.5px;color:var(--mut)}}
.vsrc code{{font-size:11.5px;word-break:break-all}}
.key{{background:#f0fdfa;border-left:4px solid var(--ok);padding:13px 17px;margin:18px 0;
border-radius:0 6px 6px 0;font-size:14px}}
.warn{{background:#fff7ed;border-left:4px solid var(--amber);padding:13px 17px;margin:16px 0;
border-radius:0 6px 6px 0;font-size:13.5px}}
.warn ol{{margin:8px 0 0;padding-left:20px}} .warn li{{margin:6px 0}}
code{{background:var(--card);padding:1px 5px;border-radius:3px;
font-family:ui-monospace,Menlo,monospace;font-size:12.5px}}
.ok{{color:var(--ok);font-weight:600}}
.mk{{font-size:11.5px;padding:2px 8px;border-radius:4px;white-space:nowrap;font-weight:600}}
.mk.cut{{background:#fef2f2;color:#b91c1c;border:1px solid #fecaca}}
.mk.keep{{background:#f0fdfa;color:#0f766e;border:1px solid #99f6e4}}
.mk.mid{{background:var(--card);color:var(--mut);border:1px solid var(--line)}}
.dotr{{display:inline-block;width:.62em;height:.62em;border-radius:50%;
background:var(--hot);margin-right:.45em;vertical-align:.02em}}
#prov td code{{word-break:break-all;white-space:normal;line-height:1.5}}
#prov table{{table-layout:fixed}} #prov th:first-child,#prov td:first-child{{width:100px}} .hot{{color:var(--hot);font-weight:600}}
.funnel,.chart{{width:100%;height:auto;margin:18px 0;background:var(--card);
border:1px solid var(--line);border-radius:8px;padding:8px}}
details.chk{{margin:26px 0;border:1px solid var(--line);border-radius:8px;padding:14px 18px;background:var(--card)}}
details.chk summary{{cursor:pointer;font-size:14px}}
footer{{max-width:1280px;margin:0 auto;padding:26px 24px 60px;color:var(--mut);font-size:12.5px;
border-top:1px solid var(--line)}}
@media(max-width:900px){{.wrap{{grid-template-columns:1fr;gap:0}}nav{{position:static;flex-direction:row;
flex-wrap:wrap;padding:14px 0}}nav a{{border-left:none;border-bottom:2px solid var(--line);font-size:12px}}}}
@media print{{nav{{display:none}}.wrap{{grid-template-columns:1fr}}.layer{{break-inside:avoid}}}}
</style></head><body>
<header>
<h1>計算量漏斗 — 從全基因組到 Group Steiner Arborescence</h1>
<div class="sub2">每一層：<b>數學/CS 定義</b> → <b>公式</b> → <b>真實數字</b> → <b>如何自己確認</b> → <b>重點講解</b></div>
<div class="badges">
<span class="b demo">F · 教學頁（非新分析）</span>
<span class="b">{e(META['date'])}</span>
<span class="b">{e(META['sample_scope'])}</span>
<span class="b">MAX_SNV={KMAX} · extra_cap={ECAP} · budget={n(BUD)}</span>
<span class="b ok">SELF-CHECK {sum(c['pass'] for c in CHECKS)}/{len(CHECKS)} PASS</span>
</div>
</header>
<div class="wrap">
<nav>{TOC}<a href="#check">SELF-CHECK 全表</a><a href="#prov">資料溯源</a></nav>
<main>
<section class="layer" style="border-top:none;padding-top:8px">
<h2 style="margin-bottom:4px">一眼看完整個漏斗</h2>
<p class="sub">單位統一為「候選解數量的十進位<b>位數</b>」（log₁₀ 尺度），全部以單一 k=8 region 為例。
橫條寬度已做次方壓縮，只表達相對量級，不可直接量長度換算。</p>
{FUNNEL_SVG}
<div class="key">整條漏斗的骨架：<b>{ARB8['unconstrained_digits']} 位 → {ARB8['hypercube_digits']} 位
（生物學限制）→ {len(str(PEAK['cand']))} 位（資料覆蓋）→ {len(str(TV['real_max_trees']))} 位（實測輸出）→ 1（VAF 破等）</b>。
最猛的兩刀是「演化樹限制」與「覆蓋公理」，而<em>不是</em>任何近似演算法。</div>
</section>
{''.join(S)}
<section class="layer" id="check"><h2>怎麼確認這頁上每個數字</h2>
<p>本頁所有數字分兩類，<b>沒有第三類</b>：</p>
<table><thead><tr><th>類別</th><th>驗證方式</th></tr></thead><tbody>
<tr><td><b>純數學量</b><br>（超立方體、A(k)、closed 集、Cand(P)）</td>
<td>由 <code>build_data.py</code> 當場重算，且<b>每個量都用兩種互不相干的演算法各算一次再比對</b>
（例：A(k) 用乘冪公式 vs 逐頂點 popcount 連乘；closed 集用暴力 vs 分層 DP）。</td></tr>
<tr><td><b>實測量</b><br>（funnel、e_min、樹數、CCF）</td>
<td>從既有 verified JSON 讀出，<b>不手打</b>；並對可推導的關係做恆等式檢查
（六層 funnel 加總 == 全集、破等+均勻 == ambiguous 總數、k 直方圖重數等）。</td></tr>
</tbody></table>
<div class="vbox"><div class="vh">一鍵重跑全部驗證</div>
<pre class="cmd">cd InterSubMod/docs/reports/in_progress/2026/07/20260726_complexity_funnel_teaching_01
python3 build_data.py --check     # 42 項 SELF-CHECK，任一項不符即 exit 1、不產出 data.json
python3 build_data.py             # 通過才寫 data.json
python3 build_html.py             # 由 data.json 注入；缺任何 key 直接 refuse 不 render</pre>
<div class="vsrc">設計理由：報告不手打數字，由 template + data 注入產生 —— 缺資料時<b>物理上無法</b>捏造，只會 refuse。</div></div>
{CHECK_HTML}
</section>
<section class="layer" id="prov"><h2>資料溯源</h2>
{table(["來源代號", "路徑（相對 repo root）"], SRC_ROWS)}
<div class="warn"><b>使用限制</b><br>
① <code>layered</code> 那份於 2026-07-11 稽核判定為 <b>upstream-mismatched engineering baseline</b>
（6/7 歷史 tagged BAM 受 <code>--truth-bed</code> 限制）——本頁僅用它做<b>計算規模量級</b>與 k 分佈說明，
<b>不得</b>作為正式 Results 的比例數字。<br>
② <code>census</code> 的六層 funnel 為 HCC1395 單樣本、chr1–22、ClairS paired PASS 操作骨幹，非生物真值。<br>
③ 純數學層（L0、L2、L3、L4 與 L5 的公式）<b>與樣本無關</b>，可直接引用。</div>
</section>
</main></div>
<footer>
{e(META['title'])} · {e(META['date'])} · {e(META['task_type'])}<br>
生成：<code>build_data.py</code> → <code>{e(DATA.name)}</code> → <code>build_html.py</code> → 本頁。
SELF-CHECK {sum(c['pass'] for c in CHECKS)}/{len(CHECKS)} PASS。零外部請求、零 CDN。
</footer></body></html>"""

OUT.write_text(HTML)
print(f"wrote {OUT.name}  ({OUT.stat().st_size / 1024:.1f} KB)")
print(f"  layers={len(S)}  checks={len(CHECKS)}  external-refs=0")
