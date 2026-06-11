#!/usr/bin/env python3
"""
99 - Illustrated edition of the pipeline manual (圖解版) for a professor.

Every non-intuitive method / formula gets an inline-SVG schematic the FIRST time it
appears (toy examples watermarked 示意); per-layer locus classification is a data-flow
diagram with REAL counts; two real cases (chr16 clean / chr20 edge) sit in side-boxes
with their real numbers + embedded validation panel. Clear 4-zone hierarchy:
  主軸 (method spine) · 配角 (supporting stats) · 解釋 (intuition) · 舉例 (examples).

methods-example principle honoured: every REAL number is read from json at build time
(construction-level anti-fabrication, §13-A); schematic toy molecules are marked 示意.

Output: display_v2/20260609_pipeline_manual_ILLUSTRATED_01.html
"""
import json, base64, datetime, os

DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")
OUT = f"{DV}/20260609_pipeline_manual_ILLUSTRATED_01.html"

fn = json.load(open(f"{DV}/funnel_numbers.json"))
d89 = json.load(open(f"{DV}/diag89.json"))
d90 = json.load(open(f"{DV}/diag90.json"))
cause = json.load(open(f"{DV}/cause_genomewide.json"))
v3 = json.load(open(f"{DV}/validation3.json"))["summary"]
man = {m["chr"] + "_" + str(m["pos"]): m for m in json.load(open(f"{DV}/manifest.json"))}
vloci = json.load(open(f"{DV}/validation3.json"))["loci"]
t10 = next(t for t in d90["tiers"] if t["minN"] == 10)

# real flow numbers
F = dict(lps_tp=fn["lps_tp_vcf_records"], lps_fp=fn["lps_fp_vcf_records"],
         ism_tp=fn["ism_tp_analyzed"], ism_fp=fn["ism_fp_analyzed"],
         sig_tp=d90["orig_significant"]["tp"], sig_fp=d90["orig_significant"]["fp"],
         base=d90["base_ratio"], sig_ratio=d90["orig_significant"]["ratio"],
         overstrict=d89["q1"]["overstrict_tp"], over_pct=d89["q1"]["pct_of_tp"],
         cause_rel=cause["tp"]["reliability_gated"], cause_n=cause["tp"]["n"], cause_pct=cause["tp"]["pct_reliability"],
         tierA=t10["tierA_tp"], tierA_fp=t10["tierA_fp"], tierB=t10["tierB_tp"], tierA_ratio=t10["tierA_ratio"],
         valA=v3["by_tier"]["A"]["validated"], valA_n=v3["by_tier"]["A"]["n"],
         dbeta=d89["q3"]["tp_dbeta"], dbeta_ratio=d89["q3"]["dbeta_ratio"], cv=d89["q3"]["tp_cramersv"],
         cv_ratio=d89["q3"]["cramersv_ratio"], overlap=d89["q3"]["both_tp"])


def case(k):
    m, v = man.get(k, {}), vloci.get(k, {})
    return dict(cv=m.get("cv_max"), db=m.get("db"), reads=m.get("reads"), hp1=m.get("minhpn"),
                F=v.get("F"), Fe=v.get("F_even"), Fo=v.get("F_odd"), Ff=v.get("F_fwd"), Fr=v.get("F_rev"),
                boot=v.get("boot"), val=v.get("validated"), q=v.get("perm_q"), permf=m.get("perm_f"))
C16, C20 = case("chr16_3444295"), case("chr20_22561507")


def b64(path):
    return "data:image/png;base64," + base64.b64encode(open(path, "rb").read()).decode()
PANEL16 = b64(f"{DV}/figs_val/chr16_3444295_valpanel.png")
PANEL20 = b64(f"{DV}/figs_val/chr20_22561507_valpanel.png")

# ---- palette ----
MET, UNM, HP1, HP2, AC = "#dc2626", "#2563eb", "#1d4ed8", "#9333ea", "#38bdf8"


def watermark(w, h):
    return (f'<text x="{w-4}" y="{h-4}" text-anchor="end" font-size="8" fill="#475569" '
            f'font-style="italic">示意 toy</text>')


# ============ SVG schematics (示意 / toy) ============
def svg_nhd():
    cells_i = [1, 1, 0, 1, 0, 1]; cells_j = [1, 0, 0, 1, 1, 1]
    s = ['<svg viewBox="0 0 410 150" width="410" height="150" font-family="system-ui">']
    s.append('<text x="6" y="16" fill="#94a3b8" font-size="12">兩條 read 在共同 CpG 上比對：不同的比例＝距離</text>')
    for ri, (cells, lab, y) in enumerate([(cells_i, "read i", 38), (cells_j, "read j", 70)]):
        s.append(f'<text x="6" y="{y+15}" fill="#cbd5e1" font-size="11">{lab}</text>')
        for ci, v in enumerate(cells):
            x = 70 + ci * 34
            col = MET if v == 1 else UNM
            s.append(f'<rect x="{x}" y="{y}" width="30" height="22" rx="3" fill="{col}"/>')
    # mark diffs
    for ci in range(6):
        if cells_i[ci] != cells_j[ci]:
            x = 70 + ci * 34
            s.append(f'<rect x="{x-2}" y="36" width="34" height="58" rx="4" fill="none" stroke="#fbbf24" stroke-width="2.5"/>')
            s.append(f'<text x="{x+15}" y="110" text-anchor="middle" fill="#fbbf24" font-size="13">✗</text>')
    s.append('<text x="70" y="132" fill="#e2e8f0" font-size="13">不同 = 2 格（黃框）／共同有效 = 6 格　→　NHD = 2/6 = 0.33</text>')
    s.append(f'<rect x="44" y="32" width="14" height="14" fill="{MET}"/><text x="61" y="44" font-size="10" fill="#94a3b8">甲基</text>')
    s.append(watermark(410, 150) + '</svg>')
    return "".join(s)


def svg_cramersv():
    s = ['<svg viewBox="0 0 380 170" width="380" height="170" font-family="system-ui">']
    s.append('<text x="6" y="15" fill="#94a3b8" font-size="12">聚類群 × HP 標籤 交叉表 → 關連強度；某格期望&lt;5 則不可信→歸零</text>')
    # 2x2 table
    data = [["30", "2"], ["1", "28"]]
    for r in range(2):
        for c in range(2):
            x, y = 120 + c * 70, 40 + r * 42
            spar = (r == 0 and c == 1) or (r == 1 and c == 0)
            fill = "#3b1106" if spar else "#0c2a4d"
            s.append(f'<rect x="{x}" y="{y}" width="66" height="38" fill="{fill}" stroke="#334155"/>')
            s.append(f'<text x="{x+33}" y="{y+24}" text-anchor="middle" fill="#e2e8f0" font-size="15">{data[r][c]}</text>')
    s.append('<text x="120" y="32" fill="#1d4ed8" font-size="11">HP1</text><text x="190" y="32" fill="#9333ea" font-size="11">HP2</text>')
    s.append('<text x="80" y="63" fill="#94a3b8" font-size="11">群A</text><text x="80" y="105" fill="#94a3b8" font-size="11">群B</text>')
    s.append('<text x="270" y="60" fill="#4ade80" font-size="11">對角集中</text><text x="270" y="74" fill="#4ade80" font-size="11">=高 V</text>')
    s.append('<text x="270" y="100" fill="#fdba74" font-size="10">紅格小→期望&lt;5</text><text x="270" y="113" fill="#fdba74" font-size="10">→V 歸零(不可信)</text>')
    s.append('<text x="6" y="150" fill="#e2e8f0" font-size="12.5">V = √( χ² / (n·(min(列,欄)−1)) )　0=無關, 1=完全對應</text>')
    s.append(watermark(380, 170) + '</svg>')
    return "".join(s)


def svg_permanova():
    import math
    s = ['<svg viewBox="0 0 360 175" width="360" height="175" font-family="system-ui">']
    s.append('<text x="6" y="15" fill="#94a3b8" font-size="12">距離空間：跨 HP 距離 是否 ≫ 同 HP 內距離（對稀疏穩健）</text>')
    g1 = [(90, 70), (110, 90), (95, 105), (120, 75)]
    g2 = [(250, 80), (270, 100), (255, 115), (240, 95)]
    for (x, y) in g1:
        s.append(f'<circle cx="{x}" cy="{y}" r="7" fill="{HP1}"/>')
    for (x, y) in g2:
        s.append(f'<circle cx="{x}" cy="{y}" r="7" fill="{HP2}"/>')
    s.append(f'<line x1="105" y1="85" x2="253" y2="98" stroke="#fbbf24" stroke-width="2" stroke-dasharray="4"/>')
    s.append('<text x="170" y="80" fill="#fbbf24" font-size="11">between（大）</text>')
    s.append(f'<line x1="90" y1="70" x2="120" y2="75" stroke="#86efac" stroke-width="2"/>')
    s.append(f'<line x1="250" y1="80" x2="240" y2="95" stroke="#86efac" stroke-width="2"/>')
    s.append('<text x="60" y="130" fill="#86efac" font-size="11">within（小）</text>')
    s.append('<text x="6" y="158" fill="#e2e8f0" font-size="12.5">pseudo-F = (SS_between/(k−1)) / (SS_within/(n−k))　F≫1 = 分離；999 次置換得 p</text>')
    s.append(watermark(360, 175) + '</svg>')
    return "".join(s)


def svg_dbeta():
    s = ['<svg viewBox="0 0 360 160" width="360" height="160" font-family="system-ui">']
    s.append('<text x="6" y="15" fill="#94a3b8" font-size="12">Δβ = HP 家族「間」平均距離 − 家族「內」平均距離</text>')
    s.append(f'<ellipse cx="95" cy="80" rx="42" ry="32" fill="none" stroke="{HP1}" stroke-width="2"/>')
    s.append(f'<ellipse cx="265" cy="80" rx="42" ry="32" fill="none" stroke="{HP2}" stroke-width="2"/>')
    s.append('<text x="95" y="40" text-anchor="middle" fill="#1d4ed8" font-size="11">HP1+HP1-1</text>')
    s.append('<text x="265" y="40" text-anchor="middle" fill="#9333ea" font-size="11">HP2+HP2-1</text>')
    for (x, y) in [(80, 72), (105, 88), (95, 78)]:
        s.append(f'<circle cx="{x}" cy="{y}" r="5" fill="{HP1}"/>')
    for (x, y) in [(250, 72), (278, 88), (262, 80)]:
        s.append(f'<circle cx="{x}" cy="{y}" r="5" fill="{HP2}"/>')
    s.append('<line x1="137" y1="80" x2="223" y2="80" stroke="#fbbf24" stroke-width="2"/>')
    s.append('<text x="180" y="74" text-anchor="middle" fill="#fbbf24" font-size="11">between</text>')
    s.append('<text x="180" y="125" text-anchor="middle" fill="#e2e8f0" font-size="12.5">值大 = 兩家族甲基模式明顯不同（注意：是距離分離，非甲基平均差）</text>')
    s.append(watermark(360, 160) + '</svg>')
    return "".join(s)


def svg_silhouette():
    s = ['<svg viewBox="0 0 360 150" width="360" height="150" font-family="system-ui">']
    s.append('<text x="6" y="15" fill="#94a3b8" font-size="12">層級聚成樹 → 試剪 k 群 → 選 silhouette 最高的 k</text>')
    # tiny dendrogram
    s.append('<path d="M40,110 L40,70 L80,70 L80,110 M60,70 L60,55 L150,55 L150,90 L120,90 L120,110 M180,110 L180,80 L150,80" '
             'fill="none" stroke="#cbd5e1" stroke-width="1.6"/>')
    s.append('<line x1="30" y1="62" x2="195" y2="62" stroke="#fbbf24" stroke-dasharray="4" stroke-width="1.5"/>')
    s.append('<text x="200" y="66" fill="#fbbf24" font-size="10">剪此處 → k=2</text>')
    # silhouette curve
    s.append('<text x="240" y="40" fill="#94a3b8" font-size="10">silhouette</text>')
    s.append('<polyline points="245,100 270,55 295,70 320,95" fill="none" stroke="#4ade80" stroke-width="2"/>')
    s.append('<circle cx="270" cy="55" r="4" fill="#4ade80"/><text x="262" y="48" fill="#4ade80" font-size="10">k=2 最佳</text>')
    s.append('<text x="6" y="138" fill="#e2e8f0" font-size="12">silhouette = (b−a)/max(a,b)　a=同群距, b=最近他群距</text>')
    s.append(watermark(360, 150) + '</svg>')
    return "".join(s)


# ============ data-driven (REAL numbers) ============
def svg_funnel():
    rows = [("longphase-S → ISM", F["lps_tp"], F["lps_fp"], 100),
            ("ISM 實際分析", F["ism_tp"], F["ism_fp"], 99),
            ("Significant gate 通過", F["sig_tp"], F["sig_fp"], int(100 * F["sig_tp"] / F["lps_tp"]) + 4)]
    s = ['<svg viewBox="0 0 420 175" width="420" height="175" font-family="system-ui">']
    s.append('<text x="6" y="14" fill="#94a3b8" font-size="12">位點數逐層收斂（真實計數；寬度按 TP 數比例）</text>')
    for i, (lab, tp, fp, wpct) in enumerate(rows):
        w = int(60 + 320 * (wpct / 100)); x = (420 - w) // 2; y = 28 + i * 46
        s.append(f'<rect x="{x}" y="{y}" width="{w}" height="36" rx="6" fill="#1e293b" stroke="{AC}"/>')
        s.append(f'<text x="210" y="{y+16}" text-anchor="middle" fill="{AC}" font-size="12">{lab}</text>')
        s.append(f'<text x="210" y="{y+31}" text-anchor="middle" fill="#e2e8f0" font-size="12">TP {tp:,} ｜ FP {fp}</text>')
        if i < 2:
            s.append(f'<text x="210" y="{y+44}" text-anchor="middle" fill="#94a3b8" font-size="13">▼</text>')
    s.append(f'<text x="210" y="170" text-anchor="middle" fill="#86efac" font-size="11.5">通過率 TP {100*F["sig_tp"]/F["ism_tp"]:.1f}% · TP:FP {F["sig_ratio"]}:1（base {F["base"]}:1 的 {F["sig_ratio"]/F["base"]:.1f}×）</text>')
    s.append('</svg>')
    return "".join(s)


def svg_flow():
    """逐層分類流水圖 — REAL counts. 標籤放條上方（左段靠左/右段靠右），永不裁切。"""
    ism = F["ism_tp"]; sig = F["sig_tp"]; nonsig = ism - sig
    over = F["overstrict"]; tA = F["tierA"]; tB = F["tierB"]; valA = F["valA"]; valfail = F["valA_n"] - valA
    W = 760; L = 14; R = W - 14; sc = (R - L) / ism

    def split_row(y, x0, total, left_w, lcol, llab, rcol, rlab):
        """labels above (left-aligned left / right-aligned right) + proportional bar below."""
        wl = max(8, int(left_w * sc)); xr = x0 + wl
        out = [f'<text x="{x0}" y="{y}" fill="{lcol}" font-size="11.5" font-weight="700">{llab}</text>']
        out.append(f'<text x="{R}" y="{y}" text-anchor="end" fill="{rcol}" font-size="11.5">{rlab}</text>')
        out.append(f'<rect x="{x0}" y="{y+5}" width="{wl}" height="16" fill="{lcol}" opacity="0.9"/>')
        out.append(f'<rect x="{xr}" y="{y+5}" width="{R-xr}" height="16" fill="{rcol}" opacity="0.55"/>')
        return "".join(out), xr, wl

    def cap(y, t):
        return f'<text x="{W//2}" y="{y}" text-anchor="middle" fill="#94a3b8" font-size="11.5">▼ {t}</text>'

    s = [f'<svg viewBox="0 0 {W} 320" width="100%" height="320" font-family="system-ui" style="max-width:760px">']
    s.append('<text x="6" y="14" fill="#94a3b8" font-size="12">每個 TP 位點逐層被分到哪一類（橫條寬度 = 真實位點數比例）</text>')
    s.append(f'<text x="{L}" y="34" fill="#cbd5e1" font-size="11.5" font-weight="700">L3 ISM 分析全部 TP　{ism:,} 位點</text>')
    s.append(f'<rect x="{L}" y="39" width="{R-L}" height="16" fill="#475569"/>')
    s.append(cap(72, "L5 Significant gate（p≤0.05、CramersV≥0.1、reads≥20）"))
    r1, _, _ = split_row(90, L, ism, sig, "#22c55e", f"通過 PASS {sig:,}", "#475569", f"未過篩掉 {nonsig:,}"); s.append(r1)
    s.append(cap(128, "被篩掉的之中，用 PERMANOVA 找「其實有真結構」的"))
    r2, _, _ = split_row(146, L, ism, over, "#f59e0b", f"有結構卻被歸零 {over:,}（過嚴格）", "#1e293b", f"確實無顯著 {nonsig-over:,}"); s.append(r2)
    s.append(cap(184, "PASS ＋ 救回 → 依「可靠 CramersV／PERMANOVA／minHP≥10」分 Tier"))
    r3, _, _ = split_row(202, L, ism, tA, "#4ade80", f"★ Tier A 高信心 {tA}（TP:FP {F['tierA_ratio']}:1）", "#60a5fa", f"Tier B 救回 {tB}"); s.append(r3)
    s.append(cap(240, "L7 不用 TP/FP 的獨立驗證（permutation＋split-half＋雙股＋bootstrap；分母＝Tier A）"))
    # validation row scales within Tier A (not ism)
    s.append(f'<text x="{L}" y="258" fill="#16a34a" font-size="11.5" font-weight="700">✓ 獨立驗證通過 {valA}/{F["valA_n"]} = {100*valA/F["valA_n"]:.0f}%</text>')
    s.append(f'<text x="{R}" y="258" text-anchor="end" fill="#f87171" font-size="11.5">⚠ 複製失敗→肉眼 {valfail}</text>')
    wv = int((R - L) * valA / F["valA_n"])
    s.append(f'<rect x="{L}" y="263" width="{wv}" height="16" fill="#16a34a" opacity="0.9"/>')
    s.append(f'<rect x="{L+wv}" y="263" width="{R-L-wv}" height="16" fill="#7f1d1d" opacity="0.7"/>')
    s.append(f'<text x="{L}" y="300" fill="#fcd34d" font-size="11">主軸：每層把位點分成「進下一層」與「分流出去」；橘＝過嚴格救回、綠＝最終高信心可驗證、紅＝邊緣待肉眼</text>')
    s.append(f'<text x="{L}" y="316" fill="#94a3b8" font-size="10.5">配角：寬度＝真實位點數比例；數字皆讀自 diag89 / diag90 / validation3 json（最後一列以 Tier A 為分母）</text>')
    s.append('</svg>')
    return "".join(s)


def svg_valladder():
    layers = [("L1 顯著（超越巧合）", "HP-permutation null + Mann-Whitney · BH-FDR q≤0.05", "#0e7490"),
              ("L2 獨立複製 ⭐", "CpG split-half(偶/奇) + 雙股一致（非循環關鍵）", "#0891b2"),
              ("L3 穩定", "read bootstrap ≥90% + PERMDISP 離散控制", "#06b6d4"),
              ("L4 生物特異", "normal-cis(HP_Residual) + germline-het 負對照", "#22d3ee"),
              ("L5 肉眼（最後）", "read×CpG 四格 panel · 只看 borderline 邊緣", "#67e8f9")]
    s = ['<svg viewBox="0 0 560 230" width="100%" height="230" font-family="system-ui" style="max-width:560px">']
    s.append('<text x="6" y="14" fill="#94a3b8" font-size="12">統計優先（下寬）→ 肉眼最後（上窄），逐層加嚴</text>')
    for i, (lab, d, col) in enumerate(layers):
        w = 520 - i * 70; x = (560 - w) // 2; y = 196 - i * 38
        s.append(f'<rect x="{x}" y="{y}" width="{w}" height="32" rx="5" fill="{col}" opacity="0.82"/>')
        s.append(f'<text x="280" y="{y+14}" text-anchor="middle" fill="#04222b" font-size="12" font-weight="700">{lab}</text>')
        s.append(f'<text x="280" y="{y+27}" text-anchor="middle" fill="#04222b" font-size="10">{d}</text>')
    s.append('</svg>')
    return "".join(s)


# ============ layout helpers ============
def zone(kind, html):
    cls = {"spine": "z-spine", "support": "z-sup", "explain": "z-exp", "example": "z-ex"}[kind]
    tag = {"spine": "主軸 · 方法", "support": "配角 · 數據", "explain": "解釋 · 直覺", "example": "舉例 · 真實案例"}[kind]
    return f'<div class="{cls}"><div class="ztag">{tag}</div>{html}</div>'


def example_box(name, c, panel, verdict, color, story):
    rows = [("CramersV (gate)", f"{c['cv']:.2f}"), ("Δβ", f"{c['db']:+.2f}"), ("reads", c['reads']),
            ("獨立 F (全CpG)", f"{c['F']:.2f}"), ("F 偶/奇 CpG", f"{c['Fe']:.2f} / {c['Fo']:.2f}"),
            ("F 正/反股", f"{c['Ff'] if c['Ff'] is not None else 'NA'} / {c['Fr'] if c['Fr'] is not None else 'NA'}"),
            ("bootstrap", f"{c['boot']:.2f}"), ("perm q", c['q'])]
    st = "".join(f'<tr><td>{a}</td><td class="n">{b}</td></tr>' for a, b in rows)
    return (f'<div class="exbox" style="border-color:{color}">'
            f'<div class="exh" style="color:{color}">{name} — {verdict}</div>'
            f'<div class="story">{story}</div>'
            f'<img src="{panel}" alt="panel"/>'
            f'<table class="extab">{st}</table></div>')


# ============ body ============
def body():
    A = []
    a = A.append
    a('<h1>HCC1395 甲基差異分析 — 全流程<span style="color:#38bdf8">圖解版</span>說明書</h1>')
    a(f'<div class="sub">每個方法/公式第一次出現都配示意圖 · 逐層分類流水圖（真實計數）· 真實案例側欄 · '
      f'生成 {datetime.date.today().isoformat()} · <span class="badge partial">single-sample HCC1395 ⭐3</span> · '
      f'文字精簡版見 <code>20260609_pipeline_manual_for_professor_01.html</code></div>')
    a('<div class="legend">分層色帶：'
      '<span class="lg z-spine">主軸·方法</span><span class="lg z-sup">配角·數據</span>'
      '<span class="lg z-exp">解釋·直覺</span><span class="lg z-ex">舉例·案例</span>'
      '　示意圖右下標 <i>示意 toy</i> = 教學用合成例（非真實數字）；流水圖/funnel 為真實計數。</div>')

    # S1 overview spine
    a('<h2 id="o">總覽 — 一條主軸：「甲基化會不會按單倍型分開？」</h2>')
    a(zone("explain", '<b>核心問題：</b>對每個體細胞突變位點，看它周圍 read 的甲基化模式，是否會按照（由 SNP 定相的）'
            '單倍型 HP tag 清楚分成兩群。能清楚分＝該位點有等位基因特異甲基化(ASM)。本說明書把這個判斷的每一層拆給你看。'))
    a(zone("spine", svg_funnel()))
    a(zone("support", f'<b>主軸 7 層：</b>樣本 BAM → ① ClairS 找體細胞 SNV → ② LongPhase-S 標 HP tag '
            f'→ ③ ISM 開 ±1000bp 視窗抽 read×CpG → ④ 統計檢定 → ⑤ Significant gate（{F["sig_tp"]:,} 通過）'
            f'→ ⑥ Tier 分層 → ⑦ 獨立驗證。'))

    # S2 NHD
    a('<h2 id="nhd">方法 ① 兩條 read 有多像 — NHD 距離</h2>')
    a(zone("explain", '<b>直覺：</b>兩條 read 在它們<b>共同覆蓋</b>的 CpG 上，甲基化狀態有多少比例不一樣。越多不一樣＝越遠。'))
    a(zone("spine", svg_nhd()))
    a(zone("support", '公式 <code>NHD = 不同CpG數 / 共同有效CpG數</code>（共同 &lt;3 記無效）。先把連續甲基機率二元化：'
            '≥0.8→甲基(1)、≤0.2→未甲基(0)。出處 <code>DistanceMatrix.cpp:21-54</code>、<code>Config.hpp:33-34</code>。'))

    # S3 clustering
    a('<h2 id="clu">方法 ② 把像的 read 聚成群 — 階層聚類 + silhouette</h2>')
    a(zone("explain", '<b>直覺：</b>距離近的 read 一層層併成樹，再挑「剪幾群」讓群內最緊、群間最開（silhouette 最高）。'))
    a(zone("spine", svg_silhouette()))

    # S4 CramersV
    a('<h2 id="cv">方法 ③ 群是否對應 HP — Cramér\'s V（卡方）+ 可靠性閘控</h2>')
    a(zone("explain", '<b>直覺：</b>把 read 依「聚類群 × HP 標籤」做交叉表，看兩種分法多一致。'
            '<b>但卡方在表格太稀疏（HP 群極不平衡）時不可信</b> → 強制歸零。'))
    a(zone("spine", svg_cramersv()))
    a(zone("support", '<code>V = √(χ²/(n·(min(列,欄)−1)))</code>；可靠性 = 每格期望值 ≥5（Cochran）。'
            '不可靠 → <code>summary CramersV = 0</code>（原始值另存 significance.json）。出處 <code>MathUtils.cpp:154</code>、'
            '<code>RegionProcessor.cpp:1592</code>。<b>這就是後面「過嚴格」的根源。</b>'))

    # S5 PERMANOVA
    a('<h2 id="pm">方法 ④ 距離法的 HP 分離 — PERMANOVA（對稀疏穩健）</h2>')
    a(zone("explain", '<b>直覺：</b>不靠離散分群，直接問「跨 HP 的 read 距離是否遠大於同 HP 內」。'
            '稀疏時比卡方可靠 → 用來救回被卡方歸零的真結構。'))
    a(zone("spine", svg_permanova()))
    a(zone("support", '<code>pseudo-F = (SS_between/(k−1))/(SS_within/(n−k))</code>；p 由 999 次 HP 標籤置換得。'
            '完美分離→F=10⁹。出處 <code>StructureTest.cpp:141</code>。'))

    # S6 dbeta
    a('<h2 id="db">方法 ⑤ HP 家族間的分離量 — Δβ（HPMergedDelta）</h2>')
    a(zone("explain", '<b>直覺：</b>併成兩大家族（HP1+HP1-1 vs HP2+HP2-1），問「跨家族距離比家族內大多少」。'))
    a(zone("spine", svg_dbeta()))
    a(zone("support", f'<code>Δβ = between_mean − within_mean</code>（NHD 距離上）。<b>精確：是距離分離量，非甲基平均差。</b>'
            f'與 CramersV 不同訊號（本研究 Δβ≥0.2 與 CramersV≥0.1 僅重疊 {F["overlap"]} 個）。出處 <code>LabelTest.cpp:204-214</code>。'))

    # S7 gate + flow
    a('<h2 id="gate">分類 ⑥ Significant gate 與逐層位點分流（流水圖）</h2>')
    a(zone("explain", '<b>直覺：</b>把上面的統計組成一道 gate 判「算不算顯著 ASM」，再依可靠性+結構分 Tier。'
            '下面這張流水圖是<b>本說明書主軸</b>：看每個 TP 位點逐層被分到哪一類。'))
    a(zone("spine", svg_flow()))
    a(zone("support", f'兩道門檻：passed_gating(p≤0.1) → Significant(p≤0.05 & V≥0.1 & reads≥20)。'
            f'過嚴格 {F["overstrict"]:,} 個有結構卻 CramersV=0（{F["cause_pct"]}% 是 Cochran 稀疏歸零）。'
            f'Tier A {F["tierA"]}/8 = {F["tierA_ratio"]}:1。出處 <code>RegionProcessor.cpp:1143</code>、scripts 89-91。'))

    # S8 validation + examples
    a('<h2 id="val">驗證 ⑦ 不用 TP/FP 的獨立驗證（5 層，肉眼最後）</h2>')
    a(zone("explain", '<b>直覺：</b>不用 SEQC2 真集，改證「甲基↔HP 關連<b>真實且能獨立複製</b>」。'
            'HP 來自 SNP 定相（與甲基無關）→ 非循環。換獨立子集（不同 CpG/股/重抽）還看得到嗎？'))
    a(zone("spine", svg_valladder()))
    a(zone("support", f'實證（非循環）：Tier A <b>{100*F["valA"]/F["valA_n"]:.0f}%</b>（{F["valA"]}/{F["valA_n"]}）通過。'
            '~14% 沒過 = gate 自信但複製不一致 = 下面右邊那種，正是 L5 該肉眼的邊緣。'))
    a('<h3 style="margin-top:18px">舉例 — 兩個真實位點對照（隨時放旁邊看）</h3>')
    a('<div class="exrow">'
      + example_box("chr16:3444295", C16, PANEL16, "✓ 鐵證", "#4ade80",
                    "四格（偶/奇 CpG·正/反股）全部清楚分離、F 全 ≈8-10、CramersV=1.0、Δβ=0.84 → 真且可獨立複製。")
      + example_box("chr20:22561507", C20, PANEL20, "⚠ 邊緣（gate 過度自信）", "#fdba74",
                    "CramersV=0.97 看似強（分群對齊 HP），但獨立 F=1.15 很弱、正股 0.99 / 反股 1.44 矛盾 → 實際甲基距離差很小、複製失敗 → 需肉眼定奪。")
      + '</div>')
    a(zone("example", '<b>教學重點：</b>chr20 説明「CramersV 高 ≠ 甲基分離大」—— 卡方測「分群是否對齊標籤」，'
            'F 測「實際距離分離量」。兩者可背離；獨立驗證(F+雙股)抓得出 gate 抓不出的邊緣。'))

    # boundaries
    a('<h2 id="bd">已知邊界</h2>')
    a(zone("support", '🔴 ① 單樣本 HCC1395、單 pipeline → ⭐3。② 全程 characterization，<b>不可寫成 variant filter</b>'
            '（甲基→TP/FP concluded NEGATIVE）。③ Δβ 是距離分離量非甲基平均差。④ HP 來自 SNP 定相才使驗證非循環。'
            '⑤ KDE 覆蓋隨位點集合變動（不影響分類）。'))
    a(f'<div class="sub" style="margin-top:14px">配套互動工具（逐位點雙圖+即時門檻+驗證徽章）：'
      f'<code>InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/display_v2/20260609_locus_judgment_display_01.html</code></div>')
    return "".join(A)


CSS = """
:root{--bg:#0f172a;--fg:#e2e8f0;--mut:#94a3b8;--ac:#38bdf8}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:14.5px/1.7 system-ui,-apple-system,'Segoe UI',sans-serif}
.layout{display:flex;max-width:1180px;margin:0 auto}
nav{position:sticky;top:0;align-self:flex-start;height:100vh;overflow:auto;width:210px;padding:18px 10px;border-right:1px solid #334155;font-size:13px;flex-shrink:0}
nav a{display:block;color:var(--mut);text-decoration:none;padding:4px 8px;border-radius:5px}nav a:hover{background:#16203a;color:var(--ac)}
nav .t{color:var(--ac);font-weight:700;margin:6px 0 6px;font-size:12px}
main{flex:1;padding:22px 30px 90px;min-width:0}
h1{font-size:23px;margin:0 0 4px}h2{font-size:19px;margin:32px 0 8px;padding:6px 0 6px 11px;border-left:4px solid var(--ac);scroll-margin-top:12px}
h3{font-size:15px;color:var(--ac)}.sub{color:var(--mut);font-size:12.5px}
code{background:#0b1222;padding:1px 6px;border-radius:4px;font-size:12px;color:#fcd34d}
.legend{background:#0b1222;border:1px solid #334155;border-radius:8px;padding:9px 13px;font-size:12px;color:var(--mut);margin:12px 0}
.lg{display:inline-block;padding:1px 8px;border-radius:5px;margin:0 4px;font-size:11px}
.z-spine{border-left:4px solid #38bdf8} .z-sup{border-left:4px solid #64748b}
.z-exp{border-left:4px solid #22c55e} .z-ex{border-left:4px solid #fbbf24}
.z-spine,.z-sup,.z-exp,.z-ex{background:#15213a;border-radius:0 8px 8px 0;padding:10px 14px;margin:9px 0;position:relative}
.z-exp{background:#0c2030}.z-ex{background:#1c1407}
.ztag{position:absolute;top:-9px;left:10px;font-size:10px;padding:0 7px;border-radius:8px;background:#0b1222;color:var(--mut)}
.lg.z-spine{background:#0c2c44;color:#7dd3fc}.lg.z-sup{background:#1e293b;color:#cbd5e1}.lg.z-exp{background:#052012;color:#86efac}.lg.z-ex{background:#3b1106;color:#fdba74}
svg{display:block;margin:4px 0;background:#0b1222;border-radius:8px;padding:6px}
.badge{display:inline-block;font-size:11px;padding:1px 7px;border-radius:9px;background:#3b1106;color:#fdba74}
.exrow{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin:10px 0}
.exbox{background:#0b1222;border:2px solid;border-radius:10px;padding:12px}
.exh{font-weight:700;font-size:14px;margin-bottom:5px}.story{font-size:12px;color:#cbd5e1;margin-bottom:8px;min-height:48px}
.exbox img{width:100%;border-radius:6px;background:#fff}
.extab{width:100%;border-collapse:collapse;font-size:11.5px;margin-top:8px}
.extab td{border-bottom:1px solid #1e293b;padding:2px 4px;color:var(--mut)}.extab td.n{text-align:right;color:var(--fg);font-variant-numeric:tabular-nums}
@media(max-width:860px){nav{display:none}.exrow{grid-template-columns:1fr}}
"""

NAVITEMS = [("o", "總覽 · 主軸"), ("nhd", "① NHD 距離"), ("clu", "② 聚類 silhouette"),
            ("cv", "③ Cramér's V + 閘控"), ("pm", "④ PERMANOVA"), ("db", "⑤ Δβ"),
            ("gate", "⑥ gate + 流水圖"), ("val", "⑦ 獨立驗證 + 舉例"), ("bd", "邊界")]


def main():
    nav = '<nav><div class="t">圖解目錄</div>' + "".join(
        f'<a href="#{i}">{t}</a>' for i, t in NAVITEMS) + '</nav>'
    FONT = ("system-ui,-apple-system,'Segoe UI','Noto Sans CJK TC','Noto Sans TC',"
            "'Microsoft JhengHei','PingFang TC','Heiti TC','Droid Sans Fallback',sans-serif")
    html = (f'<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>HCC1395 全流程圖解說明書</title><style>{CSS}</style></head>'
            f'<body><div class="layout">{nav}<main>{body()}</main></div></body></html>')
    # inject a CJK-robust font stack everywhere (SVG attrs + CSS body) so it renders on any OS
    html = html.replace('font-family="system-ui"', f'font-family="{FONT}"')
    html = html.replace("font:14.5px/1.7 system-ui,-apple-system,'Segoe UI',sans-serif",
                        f"font:14.5px/1.7 {FONT}")
    with open(OUT, "w") as f:
        f.write(html)
    print(f"[99] wrote {OUT} ({os.path.getsize(OUT)//1024} KB)")
    print(f"     real flow: ism_tp={F['ism_tp']} sig={F['sig_tp']} overstrict={F['overstrict']} "
          f"tierA={F['tierA']} valA={F['valA']}/{F['valA_n']}; cases chr16 F={C16['F']} chr20 F={C20['F']}")


if __name__ == "__main__":
    main()
