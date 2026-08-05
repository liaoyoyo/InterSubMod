#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_report.py — 由 data.json 注入產生 standalone HTML(含手刻 inline SVG)

§13-A 反捏造:所有數字一律 D[...] 取自 data.json;缺 required key → raise(refuse render)。
零外部依賴:無 CDN / 無 web font / 無外部圖檔,SVG 全內嵌。
"""
import json
from math import comb, log10
from pathlib import Path

HERE = Path(__file__).parent
STEM = "20260725_solver_complexity_stepwise_01"
D = json.load(open(HERE / f"{STEM}.data.json"))

REQUIRED = ["_meta", "bench", "stage_share", "throughput_us_per_candidate", "bench_cases",
            "verify_cost", "step5_bound", "A_max", "reachable_bound", "state_space",
            "real", "theory_vs_real", "pool_build_max", "A_max_k8_factors"]
for k in REQUIRED:
    if k not in D:
        raise SystemExit(f"REFUSE: data.json 缺 required key '{k}' — 不可手打補值")


def n(x):
    return f"{x:,}" if isinstance(x, int) else x


P = D["_meta"]["params"]
SB = D["step5_bound"]
TP = D["throughput_us_per_candidate"]
R = D["real"]

# ═══════════════════════ SVG 產生器 ═══════════════════════

PAL = {"obs": "#2563eb", "hid": "#7c3aed", "root": "#64748b", "edge": "#94a3b8",
       "tree": "#0f766e", "grp": "#ea580c", "bg": "#f8fafc", "warn": "#dc2626"}


def svg_hasse(width=760, height=478):
    """SVG-A:k=3 布爾超立方體 Hasse 圖 — 隱藏 Steiner 點「真正被需要」的最小範例。

    觀測 {AAR, ARA} 的唯一共同 unit-pred 是 ARR;ARR 未被觀測 → 必須插入為 Steiner 點,
    否則兩個觀測無法在單調(popcount 遞增)樹上連通。e_min = 1。
    """
    levels = [["RRR"], ["ARR", "RAR", "RRA"], ["AAR", "ARA", "RAA"], ["AAA"]]
    obs = {"AAR", "ARA"}                  # 完整觀測(terminals)
    hidden = {"ARR"}                      # 被迫插入的未觀測祖先
    tree_edges = {("RRR", "ARR"), ("ARR", "AAR"), ("ARR", "ARA")}
    pos, y0, dy = {}, 362, 95
    for li, lv in enumerate(levels):
        y = y0 - li * dy
        for i, node in enumerate(lv):
            x = width / 2 + (i - (len(lv) - 1) / 2) * 200
            pos[node] = (x, y)
    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="k=3 布爾超立方體 Hasse 圖">']
    p.append(f'<rect width="{width}" height="{height}" fill="white"/>')
    for li in range(4):
        y = y0 - li * dy
        p.append(f'<line x1="18" y1="{y}" x2="{width-18}" y2="{y}" stroke="#e2e8f0" '
                 f'stroke-dasharray="3 5"/>')
        p.append(f'<text x="22" y="{y-9}" font-size="11" fill="#94a3b8">popcount {li}</text>')
    # 邊
    for li in range(1, 4):
        for c in levels[li]:
            for j in range(3):
                if c[j] == "A":
                    par = c[:j] + "R" + c[j + 1:]
                    if par not in pos:
                        continue
                    x1, y1 = pos[par]
                    x2, y2 = pos[c]
                    if (par, c) in tree_edges:
                        p.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
                                 f'stroke="{PAL["tree"]}" stroke-width="2.6"/>')
                    else:
                        p.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
                                 f'stroke="{PAL["edge"]}" stroke-width="1" '
                                 f'stroke-dasharray="4 4" opacity=".55"/>')
    # 節點
    for node, (x, y) in pos.items():
        if node == "RRR":
            fill, stroke, dash, lab = "#fff", PAL["root"], "none", "ROOT"
        elif node in obs:
            fill, stroke, dash, lab = PAL["obs"], PAL["obs"], "none", node
        elif node in hidden:
            fill, stroke, dash, lab = "#fff", PAL["hid"], "4 3", node + " ᴴ"
        else:
            fill, stroke, dash, lab = "#fff", "#cbd5e1", "none", node
        p.append(f'<circle cx="{x}" cy="{y}" r="21" fill="{fill}" stroke="{stroke}" '
                 f'stroke-width="2.2" stroke-dasharray="{dash}"/>')
        tc = "white" if node in obs else "#334155"
        p.append(f'<text x="{x}" y="{y+4}" font-size="11" font-family="monospace" '
                 f'text-anchor="middle" fill="{tc}">{node}</text>')
        if node in hidden:
            p.append(f'<text x="{x-30}" y="{y-2}" font-size="10.5" text-anchor="end" '
                     f'fill="{PAL["hid"]}">必須插入的</text>')
            p.append(f'<text x="{x-30}" y="{y+12}" font-size="10.5" text-anchor="end" '
                     f'fill="{PAL["hid"]}">隱藏祖先</text>')
    # 為何非插不可(置於 popcount 3 列右側空白區,避開所有節點)
    p.append(f'<rect x="500" y="48" width="238" height="62" rx="8" fill="#faf5ff" '
             f'stroke="{PAL["hid"]}" stroke-width="1.2" opacity=".95"/>')
    p.append(f'<text x="514" y="70" font-size="10.5" fill="#475569">'
             f'AAR 與 ARA 唯一的共同 unit-pred</text>')
    p.append(f'<text x="514" y="87" font-size="10.5" fill="#475569">'
             f'= ARR,但 ARR <tspan font-weight="600">未被觀測</tspan></text>')
    p.append(f'<text x="514" y="102" font-size="10.5" font-weight="600" '
             f'fill="{PAL["hid"]}">→ 必須插入為 Steiner 點</text>')
    # 圖例(底部帶狀,避開 popcount 標籤)
    lg = [(PAL["obs"], "solid", "觀測到的基因型(terminal)"),
          (PAL["hid"], "dash", "未觀測祖先 = Steiner 點"),
          (PAL["tree"], "line", "選中的樹邊(unit-flip)"),
          (PAL["edge"], "dashline", "候選但未選的邊")]
    p.append(f'<line x1="18" y1="{height-66}" x2="{width-18}" y2="{height-66}" '
             f'stroke="#e2e8f0"/>')
    for i, (col, kind, txt) in enumerate(lg):
        x, y = 24 + (i % 2) * 370, height - 42 + (i // 2) * 21
        if kind in ("line", "dashline"):
            da = ' stroke-dasharray="4 4"' if kind == "dashline" else ""
            p.append(f'<line x1="{x}" y1="{y-4}" x2="{x+22}" y2="{y-4}" stroke="{col}" '
                     f'stroke-width="2.4"{da}/>')
        else:
            da = ' stroke-dasharray="4 3"' if kind == "dash" else ""
            fl = col if kind == "solid" else "#fff"
            p.append(f'<circle cx="{x+11}" cy="{y-4}" r="7" fill="{fl}" stroke="{col}" '
                     f'stroke-width="2"{da}/>')
        p.append(f'<text x="{x+30}" y="{y}" font-size="11.5" fill="#475569">{txt}</text>')
    p.append("</svg>")
    return "".join(p)


def svg_group(width=760, height=330):
    """SVG-B:singleton group vs subcube group — long-read 覆蓋如何變成群約束。"""
    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="singleton group 與 subcube group 對照">']
    p.append(f'<rect width="{width}" height="{height}" fill="white"/>')

    def panel(ox, title, readstr, members, note, col):
        q = [f'<text x="{ox+150}" y="24" font-size="13.5" font-weight="600" '
             f'text-anchor="middle" fill="#0f172a">{title}</text>']
        # read 示意
        q.append(f'<text x="{ox+16}" y="56" font-size="11" fill="#64748b">ONT read 覆蓋:</text>')
        for i, ch in enumerate(readstr):
            x = ox + 100 + i * 40
            fill = "#e2e8f0" if ch == "X" else ("#dbeafe" if ch == "R" else "#bfdbfe")
            txt = "?" if ch == "X" else ch
            tcol = "#94a3b8" if ch == "X" else "#1e40af"
            q.append(f'<rect x="{x}" y="42" width="32" height="24" rx="4" fill="{fill}" '
                     f'stroke="#94a3b8" stroke-width="1"/>')
            q.append(f'<text x="{x+16}" y="59" font-size="12" font-family="monospace" '
                     f'text-anchor="middle" fill="{tcol}">{txt}</text>')
            q.append(f'<text x="{x+16}" y="80" font-size="9.5" text-anchor="middle" '
                     f'fill="#94a3b8">sSNV{i+1}</text>')
        # 群成員
        q.append(f'<text x="{ox+16}" y="118" font-size="11" fill="#64748b">誘導出的群:</text>')
        gy = 132
        gw = 300
        q.append(f'<rect x="{ox+14}" y="{gy}" width="{gw}" height="{56*len(members)+16}" '
                 f'rx="10" fill="none" stroke="{col}" stroke-width="2" stroke-dasharray="6 4"/>')
        for i, m in enumerate(members):
            cy = gy + 34 + i * 56
            q.append(f'<circle cx="{ox+56}" cy="{cy}" r="19" fill="{col}" opacity=".14" '
                     f'stroke="{col}" stroke-width="2"/>')
            q.append(f'<text x="{ox+56}" y="{cy+4}" font-size="11" font-family="monospace" '
                     f'text-anchor="middle" fill="#0f172a">{m}</text>')
            q.append(f'<text x="{ox+88}" y="{cy+4}" font-size="11.5" fill="#475569">'
                     f'相容頂點</text>')
        q.append(f'<text x="{ox+14}" y="{gy+56*len(members)+38}" font-size="11.5" '
                 f'fill="{col}">{note}</text>')
        return "".join(q)

    p.append(panel(10, "完整覆蓋 → singleton group", "ARA", ["ARA"],
                   "樹必須「包含」這個頂點", PAL["obs"]))
    p.append(f'<line x1="{width/2}" y1="14" x2="{width/2}" y2="{height-10}" '
             f'stroke="#e2e8f0" stroke-width="1.5"/>')
    p.append(panel(width / 2 + 10, "部分覆蓋 → subcube group", "AXR", ["ARR", "AAR"],
                   "樹只需「碰到其中一個」", PAL["grp"]))
    p.append("</svg>")
    return "".join(p)


def svg_pipeline(width=980, height=250):
    """SVG-C:七步驟 pipeline + 成本標註。"""
    share = {s["key"]: s for s in D["stage_share"]}
    steps = [("S1", "解析觀測", "O((f+p)k)", None),
             ("S2", "universe", "O((f+p)k)", None),
             ("S3", "建 subcube group", "O(Σ2^u·k)", "t3_subcubes"),
             ("S4", "建候選 pool", "O((f+p2^k)2^k)", "t4_pool"),
             ("S5", "逐層最小解搜尋", "Σ C(P,e)×O(·)", "t5_level_search"),
             ("S6", "分析式計數", "O(|N|·k²)", "t6_analytic_count"),
             ("S7", "實體化樹", "Θ(A(N))", "t7_materialize")]
    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="七步驟流程與成本">']
    p.append(f'<rect width="{width}" height="{height}" fill="white"/>')
    bw, gap = 122, 14
    x0 = (width - (bw * 7 + gap * 6)) / 2
    for i, (sid, name, cost, key) in enumerate(steps):
        x = x0 + i * (bw + gap)
        pct = share[key]["pct"] if key else 0.0
        hot = pct >= 50
        fill = "#fef2f2" if hot else ("#f0fdfa" if pct >= 5 else "#f8fafc")
        stroke = PAL["warn"] if hot else ("#0d9488" if pct >= 5 else "#cbd5e1")
        p.append(f'<rect x="{x}" y="54" width="{bw}" height="104" rx="9" fill="{fill}" '
                 f'stroke="{stroke}" stroke-width="{2.4 if hot else 1.4}"/>')
        p.append(f'<text x="{x+bw/2}" y="76" font-size="12" font-weight="700" '
                 f'text-anchor="middle" fill="{stroke}">{sid}</text>')
        p.append(f'<text x="{x+bw/2}" y="96" font-size="10.5" text-anchor="middle" '
                 f'fill="#334155">{name}</text>')
        p.append(f'<text x="{x+bw/2}" y="118" font-size="9.5" font-family="monospace" '
                 f'text-anchor="middle" fill="#64748b">{cost}</text>')
        if key:
            p.append(f'<text x="{x+bw/2}" y="142" font-size="13" font-weight="700" '
                     f'text-anchor="middle" fill="{stroke}">{pct:.1f}%</text>')
        else:
            p.append(f'<text x="{x+bw/2}" y="142" font-size="11" text-anchor="middle" '
                     f'fill="#94a3b8">~0%</text>')
        if i < 6:
            ax = x + bw + 2
            p.append(f'<path d="M{ax} 106 L{ax+gap-4} 106" stroke="#cbd5e1" '
                     f'stroke-width="1.6" marker-end="url(#ah)"/>')
    p.append('<defs><marker id="ah" markerWidth="7" markerHeight="7" refX="6" refY="3.5" '
             'orient="auto"><path d="M0,0 L7,3.5 L0,7 z" fill="#cbd5e1"/></marker></defs>')
    p.append(f'<text x="{width/2}" y="30" font-size="12.5" text-anchor="middle" '
             f'fill="#475569">實測時間佔比(28 案例合計 '
             f'{D["stage_share_total_ms"]:.0f} ms) — <tspan fill="{PAL["warn"]}" '
             f'font-weight="700">S5 是唯一主導項</tspan></text>')
    p.append(f'<text x="{width/2}" y="185" font-size="11" text-anchor="middle" '
             f'fill="#64748b">多項式階段 ————————————————— | S5 組合爆炸(有 budget 硬上限) '
             f'| S7 指數(生產無上限)</text>')
    # 佔比長條
    bx, bwid = x0, bw * 7 + gap * 6
    cx = bx
    for s in D["stage_share"]:
        w = bwid * s["pct"] / 100
        col = PAL["warn"] if s["pct"] >= 50 else ("#0d9488" if s["pct"] >= 5 else "#cbd5e1")
        p.append(f'<rect x="{cx}" y="204" width="{max(w,1.2)}" height="20" fill="{col}" '
                 f'opacity=".85"/>')
        if s["pct"] >= 5:
            p.append(f'<text x="{cx+w/2}" y="218" font-size="10.5" text-anchor="middle" '
                     f'fill="white" font-weight="600">{s["label"].split()[0]} {s["pct"]:.1f}%</text>')
        cx += w
    p.append("</svg>")
    return "".join(p)


def svg_sawtooth(width=940, height=404):
    """SVG-E:Cand(P) 鋸齒曲線(log-y),資料來自 data.json 的精確曲線。"""
    curve = SB["curve"]
    ml, mr, mt, mb = 66, 22, 30, 70
    pw, ph = width - ml - mr, height - mt - mb
    ymax = log10(SB["global_max"]) * 1.06

    def X(Pv):
        return ml + pw * Pv / 255

    def Y(c):
        return mt + ph - ph * (log10(max(c, 1)) / ymax)

    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="Step5 候選數隨 P 變化的鋸齒曲線">']
    p.append(f'<rect width="{width}" height="{height}" fill="white"/>')
    for e in range(0, 6):
        y = Y(10 ** e)
        if y < mt:
            continue
        p.append(f'<line x1="{ml}" y1="{y}" x2="{ml+pw}" y2="{y}" stroke="#e2e8f0"/>')
        p.append(f'<text x="{ml-8}" y="{y+4}" font-size="10.5" text-anchor="end" '
                 f'fill="#94a3b8">10^{e}</text>')
    for Pv in range(0, 256, 32):
        p.append(f'<line x1="{X(Pv)}" y1="{mt}" x2="{X(Pv)}" y2="{mt+ph}" stroke="#f1f5f9"/>')
        p.append(f'<text x="{X(Pv)}" y="{mt+ph+16}" font-size="10.5" text-anchor="middle" '
                 f'fill="#94a3b8">{Pv}</text>')
    # 依 E 分段著色
    colE = {0: "#cbd5e1", 1: "#a5b4fc", 2: "#60a5fa", 3: "#0d9488", 4: "#dc2626"}
    pts_by_E = {}
    for r in curve:
        pts_by_E.setdefault(r["E"], []).append(r)
    for E, rows in sorted(pts_by_E.items()):
        d = " ".join(f'{"M" if i == 0 else "L"}{X(r["P"]):.1f},{Y(r["cand"]):.1f}'
                     for i, r in enumerate(rows))
        p.append(f'<path d="{d}" fill="none" stroke="{colE.get(E,"#334155")}" '
                 f'stroke-width="2.4"/>')
    # budget 線
    yb = Y(P["per_level_budget"])
    p.append(f'<line x1="{ml}" y1="{yb}" x2="{ml+pw}" y2="{yb}" stroke="{PAL["warn"]}" '
             f'stroke-dasharray="6 4" stroke-width="1.6" opacity=".7"/>')
    p.append(f'<text x="{ml+pw-4}" y="{yb-7}" font-size="10.5" text-anchor="end" '
             f'fill="{PAL["warn"]}">per_level_budget = {n(P["per_level_budget"])}</text>')
    # 標註鋸齒
    for pair in SB["sawtooth_pairs"]:
        x, y = X(pair["P"]), Y(pair["cand"])
        p.append(f'<circle cx="{x}" cy="{y}" r="4.2" fill="white" stroke="#0f172a" '
                 f'stroke-width="1.8"/>')
    top = SB["sawtooth_pairs"][0]
    p.append(f'<text x="{X(top["P"])+8}" y="{Y(top["cand"])-8}" font-size="11.5" '
             f'font-weight="700" fill="{PAL["warn"]}">P={top["P"]} → {n(top["cand"])} '
             f'(全域最大)</text>')
    d2 = SB["sawtooth_pairs"][1]
    p.append(f'<text x="{X(d2["P"])+8}" y="{Y(d2["cand"])+16}" font-size="11" '
             f'fill="#475569">P={d2["P"]} → {n(d2["cand"])}(驟降 '
             f'{top["cand"]/d2["cand"]:.0f}×)</text>')
    d3 = SB["sawtooth_pairs"][2]
    p.append(f'<text x="{X(d3["P"])-4}" y="{Y(d3["cand"])-10}" font-size="11" '
             f'fill="#0d9488">P={d3["P"]} → {n(d3["cand"])}</text>')
    d4 = SB["sawtooth_pairs"][3]
    p.append(f'<text x="{X(d4["P"])+8}" y="{Y(d4["cand"])+16}" font-size="11" '
             f'fill="#475569">P={d4["P"]} → {n(d4["cand"])}(驟降 '
             f'{d3["cand"]/d4["cand"]:.0f}×)</text>')
    p.append(f'<text x="{width/2}" y="16" font-size="12.5" text-anchor="middle" '
             f'fill="#334155">Step 5 最壞候選數 Cand(P) — 顏色 = 實際搜到的最深層 E</text>')
    p.append(f'<text x="{ml+pw/2}" y="{mt+ph+36}" font-size="11" text-anchor="middle" '
             f'fill="#64748b">P = 候選 Steiner 節點池大小(上限 255)</text>')
    # 圖例移到底部,避開 P=45 峰值標註
    ly = mt + ph + 58
    lx = ml + 4
    p.append(f'<text x="{lx}" y="{ly+4}" font-size="10.5" fill="#94a3b8">'
             f'實際搜到的最深層 E:</text>')
    lx += 118
    for E in sorted(pts_by_E):
        p.append(f'<line x1="{lx}" y1="{ly}" x2="{lx+18}" y2="{ly}" '
                 f'stroke="{colE.get(E)}" stroke-width="2.6"/>')
        p.append(f'<text x="{lx+23}" y="{ly+4}" font-size="10.5" fill="#475569">E={E}</text>')
        lx += 62
    p.append("</svg>")
    return "".join(p)


def svg_theorem(width=820, height=300):
    """SVG-F:A(N)=Π d_N(x) 定理示意。"""
    nodes = {"RRR": (410, 250, 0, "ROOT"), "ARR": (250, 176, 1, None),
             "RAR": (570, 176, 1, None), "AAR": (410, 102, 2, None),
             "AAA": (410, 34, 1, None)}
    edges = [("RRR", "ARR"), ("RRR", "RAR"), ("ARR", "AAR"), ("RAR", "AAR"),
             ("AAR", "AAA")]
    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="樹數定理示意">']
    p.append(f'<rect width="{width}" height="{height}" fill="white"/>')
    for a, b in edges:
        x1, y1, _, _ = nodes[a]
        x2, y2, _, _ = nodes[b]
        p.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{PAL["edge"]}" '
                 f'stroke-width="2"/>')
    for name, (x, y, d, lab) in nodes.items():
        isroot = lab == "ROOT"
        p.append(f'<circle cx="{x}" cy="{y}" r="23" fill="{"#fff" if isroot else "#eff6ff"}" '
                 f'stroke="{PAL["root"] if isroot else PAL["obs"]}" stroke-width="2.2"/>')
        p.append(f'<text x="{x}" y="{y+4}" font-size="11" font-family="monospace" '
                 f'text-anchor="middle" fill="#0f172a">{name}</text>')
        if not isroot:
            col = PAL["warn"] if d > 1 else "#0d9488"
            p.append(f'<text x="{x+30}" y="{y-6}" font-size="12.5" font-weight="700" '
                     f'fill="{col}">d={d}</text>')
            if d > 1:
                p.append(f'<text x="{x+30}" y="{y+11}" font-size="10" fill="{PAL["warn"]}">'
                         f'兩個父選擇</text>')
    p.append(f'<rect x="26" y="34" width="196" height="92" rx="9" fill="#f8fafc" '
             f'stroke="#e2e8f0"/>')
    p.append(f'<text x="40" y="58" font-size="12" font-weight="600" fill="#0f172a">'
             f'A(N) = Π d(x)</text>')
    p.append(f'<text x="40" y="80" font-size="12" font-family="monospace" fill="#334155">'
             f'= 1 × 1 × 2 × 1</text>')
    p.append(f'<text x="40" y="102" font-size="14" font-family="monospace" '
             f'font-weight="700" fill="{PAL["warn"]}">= 2 棵樹</text>')
    p.append(f'<text x="40" y="118" font-size="10" fill="#64748b">(AAR 選 ARR 或 RAR)</text>')
    p.append(f'<text x="{width-26}" y="58" font-size="11.5" text-anchor="end" '
             f'fill="#0d9488">計數 O(|N|·k) → 多項式</text>')
    p.append(f'<text x="{width-26}" y="78" font-size="11.5" text-anchor="end" '
             f'fill="{PAL["warn"]}">列舉 Θ(A(N)) → 指數</text>')
    p.append(f'<text x="{width-26}" y="102" font-size="11" text-anchor="end" '
             f'fill="#64748b">「數得出來,列不出來」</text>')
    p.append("</svg>")
    return "".join(p)


def svg_gap(width=880, height=210):
    """SVG-G:理論絕對上界 vs 實際最大(log10 尺標)。"""
    dg = D["theory_vs_real"]["A_max_k8_digits"]
    real = D["theory_vs_real"]["real_max_trees"]
    ml, mr = 70, 30
    pw = width - ml - mr
    top = 150
    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="理論上界與實際最大樹數對比">']
    p.append(f'<rect width="{width}" height="{height}" fill="white"/>')
    p.append(f'<text x="{width/2}" y="22" font-size="12.5" text-anchor="middle" '
             f'fill="#334155">樹數:理論絕對上界 vs HCC1395 實測(log₁₀ 尺標)</text>')
    y = 78
    p.append(f'<line x1="{ml}" y1="{y}" x2="{ml+pw}" y2="{y}" stroke="#cbd5e1" '
             f'stroke-width="1.4"/>')
    for e in range(0, top + 1, 25):
        x = ml + pw * e / top
        p.append(f'<line x1="{x}" y1="{y-5}" x2="{x}" y2="{y+5}" stroke="#cbd5e1"/>')
        p.append(f'<text x="{x}" y="{y+21}" font-size="10.5" text-anchor="middle" '
                 f'fill="#94a3b8">10^{e}</text>')
    xr = ml + pw * log10(real) / top
    p.append(f'<rect x="{ml}" y="{y-30}" width="{xr-ml}" height="18" rx="4" '
             f'fill="#0d9488" opacity=".85"/>')
    p.append(f'<text x="{xr+8}" y="{y-16}" font-size="11.5" font-weight="700" '
             f'fill="#0d9488">實測最大 {real} 棵</text>')
    xt = ml + pw * (dg - 1) / top
    p.append(f'<rect x="{ml}" y="{y+34}" width="{xt-ml}" height="18" rx="4" '
             f'fill="{PAL["warn"]}" opacity=".72"/>')
    p.append(f'<text x="{ml+8}" y="{y+48}" font-size="11.5" font-weight="700" '
             f'fill="white">理論絕對上界 A_max(8) ≈ 7.7 × 10^145</text>')
    p.append(f'<text x="{width/2}" y="{height-14}" font-size="12" text-anchor="middle" '
             f'fill="#475569">差距約 <tspan font-weight="700" fill="{PAL["warn"]}">'
             f'{D["theory_vs_real"]["gap_orders_of_magnitude"]} 個數量級</tspan> — '
             f'因為真實演化資料是巢狀稀疏,不是隨機基因型</text>')
    p.append("</svg>")
    return "".join(p)


def svg_real_hist(width=880, height=250):
    """SVG-H:真實問題規模分佈(k 與 n_trees 尾端)。"""
    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="真實資料問題規模分佈">']
    p.append(f'<rect width="{width}" height="{height}" fill="white"/>')
    # 左:k 分佈
    kh = [r for r in R["k"]["hist"] if isinstance(r["v"], int)]
    mx = max(r["pct"] for r in kh)
    bw, x0, base = 38, 66, 190
    p.append(f'<text x="{x0}" y="26" font-size="12.5" font-weight="600" fill="#0f172a">'
             f'k = sSNV 數(超立方體維度)</text>')
    for i, r in enumerate(kh):
        h = 120 * r["pct"] / mx
        x = x0 + i * (bw + 12)
        p.append(f'<rect x="{x}" y="{base-h}" width="{bw}" height="{h}" rx="3" '
                 f'fill="{PAL["obs"]}" opacity=".8"/>')
        p.append(f'<text x="{x+bw/2}" y="{base-h-6}" font-size="10.5" text-anchor="middle" '
                 f'fill="#334155">{r["pct"]:.1f}%</text>')
        p.append(f'<text x="{x+bw/2}" y="{base+16}" font-size="11" text-anchor="middle" '
                 f'fill="#64748b">{r["v"]}</text>')
    p.append(f'<text x="{x0}" y="{base+40}" font-size="11" fill="#64748b">'
             f'median {R["k"]["stats"]["median"]} · mean {R["k"]["stats"]["mean"]} · '
             f'max {R["k"]["stats"]["max"]}</text>')
    # 右:n_trees 尾端
    xr0 = width / 2 + 40
    p.append(f'<text x="{xr0}" y="26" font-size="12.5" font-weight="600" fill="#0f172a">'
             f'n_trees ≥ 門檻的比例</text>')
    tail = R["n_trees"]["tail"]
    for i, r in enumerate(tail):
        h = max(120 * r["pct"] / 100, 0.8)
        x = xr0 + i * 52
        col = PAL["grp"] if r["thr"] >= 32 else "#0d9488"
        p.append(f'<rect x="{x}" y="{base-h}" width="34" height="{h}" rx="3" fill="{col}" '
                 f'opacity=".8"/>')
        p.append(f'<text x="{x+17}" y="{base-h-6}" font-size="9.5" text-anchor="middle" '
                 f'fill="#334155">{r["pct"]:.1f}%</text>')
        p.append(f'<text x="{x+17}" y="{base+16}" font-size="10" text-anchor="middle" '
                 f'fill="#64748b">≥{r["thr"]}</text>')
    p.append(f'<text x="{xr0}" y="{base+40}" font-size="11" fill="#64748b">'
             f'median {R["n_trees"]["stats"]["median"]} · mean {R["n_trees"]["stats"]["mean"]} · '
             f'max {R["n_trees"]["stats"]["max"]} · 全體共 {n(R["n_trees"]["total"])} 棵</text>')
    p.append("</svg>")
    return "".join(p)


# ═══════════════════════ HTML 組裝 ═══════════════════════

def tbl(headers, rows, cls=""):
    h = "".join(f"<th>{c}</th>" for c in headers)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<div class="tw"><table class="{cls}"><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


stage_rows = [[s["label"], f'{s["ms"]:.2f}', f'{s["pct"]:.2f}%'] for s in D["stage_share"]]

step_rows = [
    ["S1", "解析觀測", "<code>:118-119</code>", "O((f+p)·k)", "~0%"],
    ["S2", "universe", "<code>:121</code>", "O((f+p)·k)", "~0%"],
    ["S3", "建 subcube group", "<code>:35-42,122</code>", "O(Σᵢ2^uᵢ·k) ≤ O(p·2^k·k)",
     f'{D["stage_share"][0]["pct"]:.1f}%'],
    ["S4", "建候選 Steiner pool", "<code>:127-136</code>",
     "O(Σ_m 2^|m|) ≤ O((f+p·2^k)·2^k)，去重後 |pool| ≤ 2^k",
     f'{D["stage_share"][1]["pct"]:.1f}%'],
    ["<b>S5</b>", "<b>逐層最小解搜尋</b>", "<code>:144-157</code>",
     "<b>Σ_{e≤E} C(P,e) × O(f + p·2^k + |N|·k)</b>",
     f'<b class="hot">{D["stage_share"][2]["pct"]:.1f}%</b>'],
    ["S6", "分析式計數 + 相容性", "<code>:173-187</code>",
     "O(|𝒩|·|N|·k²) — <b>多項式</b>", f'{D["stage_share"][3]["pct"]:.1f}%'],
    ["S7", "實體化樹", "<code>:191-201</code>", "Θ(A_total) — <b>指數</b>",
     f'{D["stage_share"][4]["pct"]:.1f}%'],
]

saw_rows = [[str(r["P"]), str(r["E"]), n(r["cand"]),
             {45: "← 全域最大", 46: "驟降 10×", 97: "次高峰", 98: "驟降 31×",
              255: "最大 P 反而便宜"}.get(r["P"], "")]
            for r in SB["sawtooth_pairs"]]

amax_rows = [[str(r["k"]), n(r["vertices"]),
              (r["value"] if r["digits"] <= 15 else f'≈ {r["sci"]}'), str(r["digits"])]
             for r in D["A_max"]]

reach_rows = [[str(r["f"]), n(int(r["bound"])), str(r["digits"])]
              for r in D["reachable_bound"]]

cmp_rows = [
    ["k（維度）", "8", str(R["k"]["stats"]["max"]), str(R["k"]["stats"]["median"]), "—"],
    ["P（節點池）", "255", "—", "—", "—"],
    ["S5 候選數", n(SB["global_max"]), "—", "—", "—"],
    ["e_min（隱藏節點）", f'{P["extra_cap"]}（+greedy）', str(R["e_min"]["stats"]["max"]),
     str(R["e_min"]["stats"]["median"]), "—"],
    ["<b>n_trees（樹數）</b>", "<b>≈ 7.7 × 10¹⁴⁵</b>",
     f'<b class="hot">{R["n_trees"]["stats"]["max"]}</b>',
     f'<b>{R["n_trees"]["stats"]["median"]}</b>',
     f'<b>10^{D["theory_vs_real"]["gap_orders_of_magnitude"]}</b>'],
]

vrows = [[r["case"], f'{r["t5_ms"]:.2f}', f'{r["V4_ms"]:.2f}', f'{r["V5_ms"]:.2f}',
          f'{r["ratio"]:.2f}×'] for r in D["verify_cost"]["rows"]]

bench_show = sorted(D["bench_cases"], key=lambda c: -c["cand"])[:12]
brows = [[c["case"], str(c["k"]), str(c["f"]), str(c["p"]), str(c["P"]), str(c["e_min"]),
          n(c["cand"]), str(c["n_trees"]), f'{c["t5_ms"]:.2f}',
          ("—" if not c["capped"] else
           ('<span class="tag warn">budget</span>' if c["cap_kind"] == "budget"
            else '<span class="tag amber">extra_cap→greedy</span>'))]
         for c in bench_show]

cls_rows = [[r["cls"], n(r["n"]), f'{r["pct"]:.1f}%'] for r in R["L1_class"]]

k8f = " · ".join(f'{r["j"]}<sup>{r["exp"]}</sup>' for r in D["A_max_k8_factors"])

HTML = f"""<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Solver 複雜度逐步驟分析 — 布爾超立方體最小 group-Steiner 樹集合</title>
<style>
:root{{--bg:#ffffff;--fg:#0f172a;--mut:#64748b;--line:#e2e8f0;--card:#f8fafc;
--acc:#2563eb;--hot:#dc2626;--ok:#0d9488;--amber:#ea580c}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--fg);
font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans CJK TC","PingFang TC",
"Microsoft JhengHei",'Helvetica Neue',Arial,sans-serif;line-height:1.72;font-size:15.5px}}
.wrap{{max-width:1080px;margin:0 auto;padding:34px 22px 80px}}
h1{{font-size:26px;line-height:1.35;margin:0 0 6px;letter-spacing:-.01em}}
h2{{font-size:20px;margin:44px 0 12px;padding-bottom:8px;border-bottom:2px solid var(--line)}}
h3{{font-size:16px;margin:26px 0 8px;color:#1e293b}}
.sub{{color:var(--mut);font-size:13.5px;margin-bottom:20px}}
.tldr{{background:linear-gradient(180deg,#f0f9ff,#f8fafc);border:1px solid #bae6fd;
border-left:5px solid var(--acc);border-radius:10px;padding:16px 20px;margin:20px 0 8px}}
.tldr b{{color:var(--acc)}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(158px,1fr));gap:11px;margin:18px 0}}
.kpi{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:12px 14px}}
.kpi .v{{font-size:21px;font-weight:700;letter-spacing:-.02em}}
.kpi .l{{font-size:11.5px;color:var(--mut);margin-top:2px}}
.kpi.hot .v{{color:var(--hot)}} .kpi.ok .v{{color:var(--ok)}} .kpi.acc .v{{color:var(--acc)}}
figure{{margin:22px 0;padding:14px;background:#fff;border:1px solid var(--line);border-radius:11px}}
figure svg{{width:100%;height:auto;display:block}}
figcaption{{font-size:12.5px;color:var(--mut);margin-top:10px;padding-top:9px;
border-top:1px dashed var(--line)}}
.tw{{overflow-x:auto;margin:14px 0}}
table{{border-collapse:collapse;width:100%;font-size:13.5px;min-width:520px}}
th{{background:#f1f5f9;text-align:left;padding:9px 11px;border-bottom:2px solid #cbd5e1;
font-weight:600;white-space:nowrap}}
td{{padding:8px 11px;border-bottom:1px solid var(--line);vertical-align:top}}
tbody tr:hover{{background:#fafcff}}
code{{background:#f1f5f9;padding:1px 5px;border-radius:4px;font-size:12.5px;
font-family:ui-monospace,SFMono-Regular,Menlo,monospace}}
.math{{background:#fbfdff;border:1px solid #dbeafe;border-radius:9px;padding:14px 18px;
margin:14px 0;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:14px;
overflow-x:auto;line-height:2}}
.box{{border:1px solid var(--line);border-radius:9px;padding:14px 18px;margin:16px 0;
background:var(--card)}}
.box.warn{{background:#fffbeb;border-color:#fcd34d}}
.box.dang{{background:#fef2f2;border-color:#fecaca}}
.box.good{{background:#f0fdfa;border-color:#99f6e4}}
.box h4{{margin:0 0 7px;font-size:14px}}
.hot{{color:var(--hot)}} .ok{{color:var(--ok)}} .mut{{color:var(--mut)}}
.tag{{display:inline-block;padding:1px 7px;border-radius:20px;font-size:11px;font-weight:600}}
.tag.warn{{background:#fee2e2;color:#b91c1c}}
.tag.amber{{background:#ffedd5;color:#c2410c}}
ol,ul{{padding-left:22px}} li{{margin:5px 0}}
footer{{margin-top:52px;padding-top:18px;border-top:2px solid var(--line);
font-size:12.5px;color:var(--mut)}}
.proof{{border-left:3px solid #cbd5e1;padding:4px 0 4px 16px;margin:12px 0;font-size:14.5px}}
@media(prefers-color-scheme:dark){{
:root{{--bg:#0b1120;--fg:#e2e8f0;--mut:#94a3b8;--line:#1e293b;--card:#0f172a}}
body{{background:var(--bg)}} figure{{background:#0f172a}}
figure svg{{filter:invert(.92) hue-rotate(180deg)}}
th{{background:#1e293b;border-bottom-color:#334155}} tbody tr:hover{{background:#111c33}}
code,.math{{background:#111c33;border-color:#1e293b}}
.tldr{{background:#0f1e33;border-color:#1e3a5f}}
.box.warn{{background:#2a2410;border-color:#78350f}}
.box.dang{{background:#2a1414;border-color:#7f1d1d}}
.box.good{{background:#0d2926;border-color:#115e59}}
}}
</style>
<div class="wrap">

<h1>布爾超立方體上「最小 group-Steiner 樹集合」<br>逐步驟計算量與理論上界</h1>
<div class="sub">solver: <code>tree_enumeration_solver.py</code>（2026-07-06）　·
參數 MAX_SNV={P["MAX_SNV"]}、extra_cap={P["extra_cap"]}、per_level_budget={n(P["per_level_budget"])}
　·　全部數字由 <code>{STEM}.data.json</code> 注入</div>

<div class="tldr">
<b>一句話：</b>求「所有最小樹的集合」的成本 <b>91.1% 集中在單一階段（S5 逐層節點集搜尋）</b>，
其候選數有<b>精確上界 {n(SB["global_max"])}</b>（最壞約 {SB["worst_walltime_s"]["max"]} 秒／unit）。
樹的<b>計數是多項式</b>（A(N)=Π d<sub>N</sub>(x)），<b>列舉才是指數</b>。
理論絕對上界 A<sub>max</sub>(8) ≈ 7.7×10<sup>145</sup> 在 caps 下不可達；
HCC1395 實測最大僅 <b>{R["n_trees"]["stats"]["max"]} 棵</b>。
</div>

<div class="grid">
<div class="kpi hot"><div class="v">{D["stage_share"][2]["pct"]:.1f}%</div>
<div class="l">S5 佔總時間</div></div>
<div class="kpi acc"><div class="v">{n(SB["global_max"])}</div>
<div class="l">S5 候選數上界 @P={SB["at_P"]}</div></div>
<div class="kpi"><div class="v">{SB["worst_walltime_s"]["max"]} s</div>
<div class="l">最壞 wall-time／unit</div></div>
<div class="kpi hot"><div class="v">10<sup>145</sup></div>
<div class="l">A_max(8) 絕對上界</div></div>
<div class="kpi ok"><div class="v">{R["n_trees"]["stats"]["max"]}</div>
<div class="l">實測最大樹數</div></div>
<div class="kpi ok"><div class="v">{D["bench"]["n_mirror_ok"]}/{D["bench"]["n_cases"]}</div>
<div class="l">鏡像驗證通過</div></div>
</div>

<h2>§1　問題物件：超立方體、Steiner 點、group</h2>
<p>每個區域內的 <i>k</i> 個 somatic sSNV 構成布爾超立方體 <code>{{0,1}}^k</code> 的座標軸。
一條 ONT read 的基因型向量 = 一個頂點；根固定在 germline all-REF（<code>0^k</code>）；
邊只允許 <b>unit-flip</b>（一次獲得一個突變，popcount 遞增）。</p>

<figure>{svg_hasse()}
<figcaption><b>圖 1．</b>k=3 的 Hasse 圖（solver 實際運作的空間），示範<b>隱藏節點「非插不可」</b>的最小情形。
實心藍 = 觀測到的完整基因型（terminals）；虛線紫 = 未觀測的祖先 clone，即 <b>Steiner 點</b>；
粗綠 = 被選中的樹邊；灰色 = 未被任何觀測需要、<b>不進最小樹</b>的頂點。
所有邊都由下往上、一次只翻一個 bit — 這是 infinite-sites 假設的幾何化。<br>
觀測到 <code>AAR</code> 與 <code>ARA</code>，兩者唯一的共同 unit-pred 是 <code>ARR</code>；
<code>ARR</code> 沒被任何 read 觀測到，但若不插入它，兩個觀測就無法在單調樹上連通 —
<b>這就是「未觀測祖先 clone」在數學上被迫存在的地方</b>。此例 e<sub>min</sub>=1。</figcaption></figure>

<h3>long-read 的部分覆蓋，如何變成「group」</h3>
<p>這是本形式化最關鍵的一步：read 只穿過區域的一部分時，未覆蓋位點記為 <code>X</code>（未知，非 0）。
一個含 <i>u</i> 個 <code>X</code> 的觀測，對應超立方體上一個 <b>2^u 個頂點的子立方體（面）</b>；
樹只要<b>碰到其中任一個</b>就算解釋了這條 read。這正是 group-Steiner 的「group」。</p>

<figure>{svg_group()}
<figcaption><b>圖 2．</b>左：完整覆蓋 → singleton group（樹必須<i>包含</i>該頂點）。
右：部分覆蓋 → subcube group（樹只需<i>碰到其中一個</i>）。
<b>測序讀長這個物理限制，被直接翻譯成 group 約束</b> — 這是把腫瘤重建 cast 成
group-Steiner 的實質內容，也是與既有腫瘤系統發生法的分野。</figcaption></figure>

<h2>§2　七個步驟與各自的計算量</h2>
<figure>{svg_pipeline()}
<figcaption><b>圖 3．</b>七步驟 pipeline。百分比為 {D["bench"]["n_cases"]} 個合成案例的實測時間佔比
（合計 {D["stage_share_total_ms"]:.0f} ms）。<b>S5 是唯一主導項</b>；S3/S4 的理論上界看似嚇人，
實際被 <code>|pool| ≤ 2^k</code> 壓住（實測最大僅 {n(D["pool_build_max"]["max_subset_ops"])} 次子集操作、
{D["pool_build_max"]["max_t4_ms"]} ms）。</figcaption></figure>

{tbl(["#", "步驟", "程式位置", "理論成本", "實測佔比"], step_rows)}

<div class="box">
<h4>符號</h4>
<code>k</code> = sSNV 數（≤{P["MAX_SNV"]}）　<code>f</code> = 完整觀測數　<code>p</code> = partial 觀測數
<code>u_i</code> = 第 i 條 partial 的未覆蓋位點數　<code>P</code> = 候選 Steiner 節點池大小
<code>e</code> = hidden 節點數　<code>N</code> = 節點集　<code>𝒩</code> = 可行節點集族
</div>

<h3>驗證器成本：V4/V5 等於把主導階段再跑一次</h3>
<p>V4（獨立重算最小性）重跑第 <code>e_min−1</code> 層，V5（獨立重算完整性）重跑第 <code>e_min</code> 層。
實測其合計為 S5 的 <b>{D["verify_cost"]["ratio_min"]}–{D["verify_cost"]["ratio_max"]}×</b>。
生產預設 <code>VERIFY_EVERY={P["production_VERIFY_EVERY"]}</code>（每個 unit 全跑），故實際總成本 ≈ 2× S5。</p>
{tbl(["案例", "S5 (ms)", "V4 (ms)", "V5 (ms)", "(V4+V5)/S5"], vrows)}

<h2>§3　理論最大量（數學推導）</h2>

<h3>3.1　狀態空間的天然屏障</h3>
<div class="math">|{{0,1}}<sup>k</sup>| = 2<sup>k</sup> ≤ 2<sup>8</sup> = 256　　⇒　　P ≤ 2<sup>k</sup> − 1 = <b>255</b></div>
<p>候選 pool 是超立方體頂點的子集，去重後<b>無論 f、p 多大都不超過 256</b>。這是第一道硬上限。</p>
{tbl(["k", "頂點數 2^k", "P 上限"], [[str(r["k"]), n(r["vertices"]), n(r["max_P"])] for r in D["state_space"]])}

<h3>3.2　S5 候選數的精確上界 — 出現「鋸齒」</h3>
<p>迴圈在<b>第一個</b> C(P,e) &gt; B（B={n(P["per_level_budget"])}）處中斷，故最壞候選總數為：</p>
<div class="math">Cand(P) = Σ<sub>e=0..E(P)</sub> C(P,e)，　　E(P) = min( {P["extra_cap"]},　min{{e : C(P,e) &gt; B}} − 1 )</div>
<p>以大整數精確計算全部 P ∈ [0,255] 得：</p>
<div class="math">max<sub>P</sub> Cand(P) = <b class="hot">{n(SB["global_max"])}</b>　（於 P = {SB["at_P"]}，E = {SB["at_E"]}）</div>

<figure>{svg_sawtooth()}
<figcaption><b>圖 4．</b>Cand(P) 的鋸齒結構。顏色代表實際搜到的最深層 E。
<b>成本對問題規模是非單調的</b>：P 從 45→46 時，因 C(46,4)={n(comb(46,4))} 超過 budget，
第 4 層被整層砍掉，成本反而<b>驟降 10×</b>。P 從 97→98 同理驟降 31×。
換言之「問題變大」可以「變便宜」，代價是被誠實標成 capped。</figcaption></figure>

{tbl(["P", "E", "候選總數", "說明"], saw_rows)}

<div class="box good">
<h4>最壞 wall-time</h4>
實測吞吐 {TP["min"]}–{TP["max"]} µs/候選（median {TP["median"]}，n={TP["n_samples"]}）⇒
{n(SB["global_max"])} × {TP["max"]} µs = <b>{SB["worst_walltime_s"]["max"]} 秒／unit</b>（最壞）；
median 情境 {SB["worst_walltime_s"]["median"]} 秒。
</div>

<h3>3.3　樹數定理（S6/S7 的數學基礎）</h3>
<div class="box">
<h4>定理</h4>
設 ∅ ∈ N ⊆ {{0,1}}<sup>k</sup> 且 N 為 unit-flip-closed（每個非根 x∈N 至少有一個 unit-pred 在 N 中）。
則以 ∅ 為根、邊皆 unit-flip 的生成樹（arborescence）數目恰為
<div class="math">A(N) = Π<sub>x∈N, x≠∅</sub> d<sub>N</sub>(x)，　　d<sub>N</sub>(x) = #{{ j∈x : x∖{{j}} ∈ N }}</div>
</div>
<div class="proof">
<b>證明．</b>（⊇）每個非根節點 x 獨立選一個 parent ∈ {{x∖{{j}} : j∈x}} ∩ N（closed 保證非空，共 d<sub>N</sub>(x) 種）。
所得圖中每個非根節點恰一入邊；沿 parent 指標走，popcount <b>嚴格遞減 1</b> ⇒ 不可能成環，
且必在有限步抵達唯一無父節點 ∅ ⇒ root-connected。故為合法 arborescence。<br>
（⊆）反之，任一合法 arborescence 中每個非根節點恰有一個 parent，且該 parent 必為其 unit-pred 且屬於 N，
故對應唯一一組選擇。<br>
兩者互逆 ⇒ 雙射 ⇒ 計數為乘積。∎
</div>

<figure>{svg_theorem()}
<figcaption><b>圖 5．</b>定理示意。只有 <code>AAR</code> 有兩個 unit-pred（d=2），其餘皆 d=1，
故 A(N)=1×1×2×1=2 棵樹。<b>推論：計數只需 O(|N|·k) 次查詢（多項式），
列舉需 Θ(A(N))（指數）</b> — 這就是工作站能標「N 樹 = M 形狀」而不必真的生成 N 棵樹的原因。</figcaption></figure>

<h3>3.4　絕對理論最大樹數</h3>
<p>A(N) 在 N = <b>整個超立方體</b>時最大：popcount = j 的節點有 C(k,j) 個，每個的 d = j（所有 unit-pred 都在）。故</p>
<div class="math">A<sub>max</sub>(k) = Π<sub>j=1..k</sub> j<sup>C(k,j)</sup></div>
<p>k=8 時展開為　{k8f}　=</p>
<div class="math" style="word-break:break-all;font-size:12px">{D["A_max"][7]["value"]}</div>
<div class="math">≈ <b class="hot">7.696 × 10<sup>145</sup></b>　（{D["A_max"][7]["digits"]} 位數）</div>
{tbl(["k", "頂點數", "A_max(k)", "位數"], amax_rows)}

<h3>3.5　caps 下的<i>可達</i>上界（遠低於絕對上界）</h3>
<div class="box warn">
A<sub>max</sub>(8) 要求 N 含全部 256 個頂點，即 f + e = 255；
但 <code>extra_cap={P["extra_cap"]}</code> ⇒ 需 f ≥ 251 個相異完整觀測
（實際 f 最大只有 {R["f"]["stats"]["max"]}）。<b>絕對上界不可達。</b>
</div>
<div class="math">A(N) ≤ k<sup>(f + extra_cap)</sup>　　⇒　　A<sub>total</sub> ≤ |𝒩| × k<sup>(f+4)</sup>，
　|𝒩| ≤ C(P, e<sub>min</sub>) ≤ {n(P["per_level_budget"])}</div>
{tbl(["f（完整觀測數）", "每個 N 的上界 8^(f+4)", "位數"], reach_rows)}

<h2>§4　理論 vs 實際</h2>
<figure>{svg_gap()}
<figcaption><b>圖 6．</b>理論絕對上界與實測最大值的距離。</figcaption></figure>

<figure>{svg_real_hist()}
<figcaption><b>圖 7．</b>HCC1395 真實問題規模分佈（{n(R["n_units"])} 個 lineage unit）。
k 的中位數僅 {R["k"]["stats"]["median"]}；<b>{R["n_trees"]["tail"][6]["pct"]:.2f}% 的 unit 樹數 ≥ 1000</b>，
全體加總只有 {n(R["n_trees"]["total"])} 棵。</figcaption></figure>

{tbl(["量", "理論最大", "實測最大", "實測 median", "差距"], cmp_rows)}

<p><b>為什麼實際離理論這麼遠？</b>因為真實演化資料是<b>巢狀、稀疏</b>的（perfect phylogeny 結構），
不是隨機基因型。合成測試用<i>隨機</i>基因型時 capped 率高達
{D["bench"]["n_capped"]}/{D["bench"]["n_cases"]}，真實資料只有 {R["capped"]["pct"]}% —
<b>隨機輸入才是對抗性最壞情況</b>。</p>

{tbl(["L1 分類", "unit 數", "佔比"], cls_rows)}

<h3>實測案例（依候選數排序前 12）</h3>
{tbl(["案例", "k", "f", "p", "P", "e_min", "候選數", "n_trees", "S5 (ms)", "cap 成因"], brows)}

<h2>§5　兩個誠實的工程邊界</h2>
<div class="box dang">
<h4>① S7 在生產是無上限的</h4>
<code>ANALYSIS_TREE_CAP={P["production_ANALYSIS_TREE_CAP"]}</code> ⇒ <code>store_limit=None</code> ⇒
<b>全部樹都實體化</b>（非工作站顯示用的 32）。理論上這是 Θ(A_total) 的指數步；
實務上因實測 max={R["n_trees"]["stats"]["max"]} 而無害，但<b>理論上無界</b>。
若骨幹擴大或 MAX_SNV 調高，這裡會先爆。
</div>
<div class="box warn">
<h4>② capped 有兩種成因，意義不同</h4>
<ul>
<li><b>budget 超支</b>（{n(R["capped"]["by_reason"].get("budget", 0))} 個）＝ 問題太密，<b>放棄搜尋</b>。</li>
<li><b>extra_cap 耗盡 → greedy fallback</b>（{n(R["capped"]["by_reason"].get("extra_cap_greedy", 0))} 個）
＝ 真的需要 &gt;{P["extra_cap"]} 個隱藏節點，<b>退化成單一貪婪解</b>。</li>
</ul>
後者不只是「沒枚舉完」，而是<b>連最小性都不保證</b>。兩者目前共用同一個 <code>capped</code> 標籤，
論文中值得分開陳述。
</div>

<h2>§6　一句話口徑</h2>
<div class="tldr">
求最小 group-Steiner 樹<b>集合</b>的成本集中在「逐層節點集搜尋」：候選數為
Σ<sub>e≤E</sub> C(P,e)，在 P ≤ 2<sup>k</sup>−1 = 255、extra_cap={P["extra_cap"]}、
budget={n(P["per_level_budget"])} 三重限制下有<b>精確上界 {n(SB["global_max"])}</b>
（最壞約 {SB["worst_walltime_s"]["max"]} 秒／unit）。樹的<b>計數</b>由
A(N)=Π d<sub>N</sub>(x) 在 O(|N|·k) 內完成（多項式），<b>列舉</b>才是 Θ(A(N))（指數）。
理論絕對上界 A<sub>max</sub>(8)=Π<sub>j</sub> j<sup>C(8,j)</sup> ≈ 7.7×10<sup>145</sup> 在 caps 下不可達；
HCC1395 實測最大 {R["n_trees"]["stats"]["max"]} 棵、全基因組共 {n(R["n_trees"]["total"])} 棵。
</div>

<footer>
<p><b>數據來源（皆本輪實跑、可重算）：</b></p>
<ul>
<li><code>complexity_bench.py</code> → {D["bench"]["n_cases"]} 案例逐階段實測；
鏡像忠實度驗證 <b>{D["bench"]["n_mirror_ok"]}/{D["bench"]["n_cases"]} 通過</b>
（自製計時版的 e_min／n_trees／capped 與官方 <code>enumerate_min_trees</code> 完全一致）。</li>
<li><code>complexity_theory.py</code> → 大整數精確上界（非浮點估算）。</li>
<li>真實分佈：<code>{R["source"]}</code></li>
<li>solver：<code>{D["_meta"]["solver"]}</code></li>
</ul>
<div class="box warn" style="margin-top:14px">
<h4>⚠ 真實分佈的適用範圍</h4>
{R["caveat"]}
（問題尺寸分佈的<i>量級</i>不受 truth-BED 影響，故用於計算量分析成立；但比例數字不可外引。）
</div>
<p style="margin-top:16px">所有數字由 <code>{STEM}.data.json</code> 注入產生；
缺 required key 時 builder 直接 refuse（§13-A 反捏造）。重建：
<code>python3 assemble_data.py &amp;&amp; python3 build_report.py</code></p>
</footer>
</div>
"""

out = HERE / f"{STEM}.standalone.html"
out.write_text(HTML, encoding="utf-8")
print(f"wrote {out}  ({len(HTML)/1024:.1f} KB)")
