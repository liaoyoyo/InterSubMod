#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成題目聚焦母版的 3 張待補 SVG（schematic 流程/概念示意，無分析真值）。
   風格對齊 Fig 1（methods-example）：bg #FAF9F5，藍=germline/橘=somatic，5mC red/blue ramp。
   所有示意分子打 synthetic 浮水印；含 HP1-1 的圖加 label-flip 註（lint C5）。"""
import html

C = dict(bg="#FAF9F5", ink="#141413", soft="#6B6B66", line="#E4E2DA",
         germ="#2F5597", germL="#8FAADC", som="#C55A11", somL="#ED9B6A",
         meth="#C0392B", unmeth="#2D6CB5", grey="#C9C6BD", warn="#9a6212",
         teal="#1B8A86", gold="#B8860B")

def esc(s): return html.escape(str(s), quote=True)
def T(x, y, s, sz=10, col="ink", a="start", w=None, it=False):
    extra = (f' font-weight="{w}"' if w else '') + (' font-style="italic"' if it else '')
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{sz}" fill="{C.get(col,col)}" text-anchor="{a}"{extra}>{esc(s)}</text>'
def R(x, y, w, h, fill, stroke=None, rx=0, sw=1.2, dash=None):
    st = f' stroke="{C.get(stroke,stroke)}" stroke-width="{sw}"' if stroke else ''
    da = f' stroke-dasharray="{dash}"' if dash else ''
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" fill="{C.get(fill,fill)}"{st}{da}/>'
def L(x1, y1, x2, y2, col="ink", sw=1.2, dash=None):
    da = f' stroke-dasharray="{dash}"' if dash else ''
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{C.get(col,col)}" stroke-width="{sw}"{da}/>'
def arrow(x, y, col="soft"):  # right arrow head at (x,y)
    return f'<path d="M{x-9:.1f},{y-6:.1f} L{x:.1f},{y:.1f} L{x-9:.1f},{y+6:.1f}" fill="{C.get(col,col)}"/>'

def wrap(title, sub, body, W, H):
    head = (T(12, 26, title, 16, "ink", w=800) + T(12, 46, sub, 11.5, "soft"))
    foot = T(W-12, H-10, "示意 schematic — 分子為 synthetic，非真實資料", 9, "warn", a="end")
    return (f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" role="img" '
            f'aria-label="{esc(title)}"><title>{esc(title)}</title><desc>{esc(sub)}</desc>'
            f'<rect width="{W}" height="{H}" fill="{C["bg"]}"/>\n{head}\n{body}\n{foot}\n</svg>')

# 5mC ramp cell：H=實心紅, L=空心藍框, miss=灰
def mcell(x, y, s, v):
    if v == "H": return R(x, y, s, s, "meth", rx=2)
    if v == "L": return R(x, y, s, s, "#fff", stroke="unmeth", rx=2, sw=1.6)
    return R(x, y, s, s, "grey", rx=2)

# ---------------- Fig 2：單位點一條龍 ----------------
def fig2():
    W, H = 1180, 430
    s = []
    # latent groups: A=(R1,R3,R5) germline-like, B=(R2,R4,R6) somatic-like
    rows = ["R1","R2","R3","R4","R5","R6"]
    grp = {"R1":"A","R3":"A","R5":"A","R2":"B","R4":"B","R6":"B"}
    M = {  # 5 CpG pattern：A 群偏甲基, B 群偏去甲基
        "R1":["H","H","L",".","L"], "R3":["H","H",".","L","L"], "R5":["H",".","H","L","."],
        "R2":["L","L","H","H","."], "R4":[".","L","H","H","L"], "R6":["L",".","H",".","H"],
    }
    # Panel ①：read×CpG
    ax, ay, cs = 60, 92, 26
    s.append(T(ax, ay-14, "① read×CpG 甲基矩陣", 11.5, "ink", w=700))
    for j in range(5): s.append(T(ax+34+j*cs+cs/2, ay-2, f"C{j+1}", 9, "soft", a="middle"))
    for i, r in enumerate(rows):
        ry = ay + i*cs
        col = "germ" if grp[r]=="A" else "som"
        s.append(T(ax+28, ry+cs/2+3, r, 9.5, col, a="end", w=700))
        for j in range(5):
            s.append(mcell(ax+34+j*cs, ry, cs-3, M[r][j]))
    s.append(arrow(ax+34+5*cs+34, ay+3*cs))
    # Panel ②：距離矩陣
    bx, by, bs = 470, 92, 26
    s.append(T(bx, by-14, "② read–read 距離矩陣", 11.5, "ink", w=700))
    # schematic distances: within-group small, cross-group large
    def dist(a,b):
        if a==b: return None
        return 0.12 if grp[a]==grp[b] else 0.82
    for j,r in enumerate(rows): s.append(T(bx+30+j*bs+bs/2, by-2, r, 8.5, "soft", a="middle"))
    for i,a in enumerate(rows):
        s.append(T(bx+26, by+i*bs+bs/2+3, a, 9, "soft", a="end"))
        for j,b in enumerate(rows):
            d = dist(a,b)
            if d is None:
                s.append(R(bx+30+j*bs, by+i*bs, bs-3, bs-3, "#efeee8", rx=2))
            else:
                op = 0.18 + d*0.72
                s.append(f'<rect x="{bx+30+j*bs:.1f}" y="{by+i*bs:.1f}" width="{bs-3:.1f}" height="{bs-3:.1f}" rx="2" fill="{C["teal"]}" opacity="{op:.2f}"/>')
    s.append(T(bx+30, by+6*bs+16, "深=遠 / 淺=近（組內近、跨組遠）", 9, "soft"))
    s.append(arrow(bx+30+6*bs+30, by+3*bs))
    # Panel ③：UPGMA 分群 + tree
    tx, ty = 880, 92
    s.append(T(tx, ty-14, "③ UPGMA 分群 → 2 群", 11.5, "ink", w=700))
    # two clusters
    gA = ["R1","R3","R5"]; gB = ["R2","R4","R6"]
    yy = ty
    def leaf(label, y, col):
        return T(tx+150, y+4, label, 10, col, a="start", w=700)
    # group A box
    s.append(R(tx, ty, 250, 78, "#eef3f9", stroke="germ", rx=8, sw=1.4))
    s.append(T(tx+10, ty+18, "群 1 ≈ HP1（germline）", 10.5, "germ", w=700))
    for k,r in enumerate(gA): s.append(T(tx+18, ty+36+k*15, r, 9.5, "germ"))
    s.append(R(tx, ty+96, 250, 78, "#fdf1e8", stroke="som", rx=8, sw=1.4))
    s.append(T(tx+10, ty+114, "群 2 ≈ HP1-1（somatic）", 10.5, "som", w=700))
    for k,r in enumerate(gB): s.append(T(tx+18, ty+132+k*15, r, 9.5, "som"))
    # mini dendrogram on the right
    dx = tx+170
    s.append(L(dx, ty+39, dx+30, ty+39, "germ", 1.6))
    s.append(L(dx, ty+135, dx+30, ty+135, "som", 1.6))
    s.append(L(dx+30, ty+39, dx+30, ty+135, "ink", 1.4))
    s.append(L(dx+30, ty+87, dx+55, ty+87, "ink", 1.6))
    s.append(T(dx+58, ty+90, "分歧", 9, "soft"))
    # legend + label-flip + caption
    ly = 372
    s.append(R(60, ly, 14, 14, "meth", rx=2)); s.append(T(80, ly+11, "5mC 甲基", 10, "ink"))
    s.append(R(180, ly, 14, 14, "#fff", stroke="unmeth", rx=2, sw=1.6)); s.append(T(200, ly+11, "去甲基", 10, "ink"))
    s.append(R(280, ly, 14, 14, "grey", rx=2)); s.append(T(300, ly+11, "未覆蓋", 10, "ink"))
    s.append(T(60, ly+30, "▣ label-flip：HP1-1 ≡ HP1 germline 鏈上帶 somatic ALT 的子單倍型（非第三 haplotype）；甲基矩陣值為示意。", 10, "warn", w=700))
    return wrap("圖 2 · 單一位點上 ISM 怎麼跑（一條龍）",
                "read×CpG 甲基矩陣 → read–read 距離矩陣 → UPGMA 分群（示意一個位點）", "\n".join(s), W, H)

# ---------------- Fig 5：軸 C vs 軸 A ----------------
def fig5():
    W, H = 1180, 430
    s = []
    # 中間：問題不同
    s.append(T(W/2, 78, "同一位點，兩種問法 — 互補，不是「更好」", 13.5, "ink", a="middle", w=800))
    # 左：軸 A 率差
    ax = 60
    s.append(R(ax, 100, 480, 250, "#fff", stroke="line", rx=10))
    s.append(T(ax+18, 128, "軸 A · per-position 率差", 13, "teal", w=800))
    s.append(T(ax+18, 148, "modkit dmr · DSS · methylKit", 10.5, "soft"))
    s.append(T(ax+18, 174, "把 reads 塌成「每個位點的甲基 %」，比兩組率差", 10.5, "ink"))
    # two rate bars HP1 vs HP1-1
    by = 196; bw = 300
    for k,(lab,rate,col) in enumerate([("HP1", 0.62, "germ"), ("HP1-1", 0.50, "som")]):
        yy = by + k*46
        s.append(T(ax+18, yy+16, lab, 10.5, col, w=700))
        s.append(R(ax+90, yy, bw, 22, "#eee", rx=4))
        s.append(R(ax+90, yy, bw*rate, 22, col, rx=4))
        s.append(T(ax+90+bw*rate+8, yy+16, f"{int(rate*100)}%", 10, col, w=700))
    s.append(T(ax+18, by+104, "輸出：Δ 率差（一個數）", 10.5, "ink", w=700))
    s.append(T(ax+18, by+124, "→ 看不到 read 內部結構", 10, "warn"))
    # 右：軸 C 結構
    cx = 640
    s.append(R(cx, 100, 480, 250, "#fff", stroke="line", rx=10))
    s.append(T(cx+18, 128, "軸 C · read-to-read 結構（ISM）", 13, "gold", w=800))
    s.append(T(cx+18, 148, "ISM · cvlr · ASMS · qFDRP", 10.5, "soft"))
    s.append(T(cx+18, 174, "保留 read×CpG，問「reads 是否分成可分離亞群」", 10.5, "ink"))
    # mini matrix → 2 clusters
    mx, my, cs = cx+30, 196, 16
    pat = [["H","L","H","L"],["H","L","L","L"],["L","H","L","H"],["L","H","H","H"]]
    rc = ["germ","germ","som","som"]
    for i in range(4):
        s.append(T(mx-6, my+i*cs+cs-3, f"R{i+1}", 8, rc[i], a="end"))
        for j in range(4):
            s.append(mcell(mx+j*cs, my+i*cs, cs-3, pat[i][j]))
    s.append(arrow(mx+4*cs+24, my+2*cs))
    # clusters
    clx = mx+4*cs+34
    s.append(R(clx, my-4, 150, 30, "#eef3f9", stroke="germ", rx=6)); s.append(T(clx+8, my+15, "群1 ≈ HP1", 9.5, "germ", w=700))
    s.append(R(clx, my+34, 150, 30, "#fdf1e8", stroke="som", rx=6)); s.append(T(clx+8, my+53, "群2 ≈ HP1-1", 9.5, "som", w=700))
    s.append(T(cx+18, by+104, "輸出：PERMANOVA 結構顯著性 + 距離矩陣", 10.5, "ink", w=700))
    s.append(T(cx+18, by+124, "→ 保留誰跟誰像（結構）", 10, "ink"))
    # 中央連接 + 誠實註
    s.append(T(W/2, 372, "✅ modkit 已交叉驗證率差層 r≈0.98（同方向同量級）　🔴 discordant-case（率差≈0 但結構顯著）未跑＝硬傷　·　鐵則：組合無人佔 ≠ 更好", 10, "ink", a="middle"))
    s.append(T(60, 404, "▣ label-flip：HP1-1 ≡ HP1 上帶 somatic ALT 的子單倍型（非第三 haplotype）。", 9.5, "warn", w=700))
    return wrap("圖 5 · ISM 站哪條軸：read 結構（軸 C） vs 率差（軸 A）",
                "與 modkit/DSS 互補；ISM 保留 read-level 結構並做顯著性檢定", "\n".join(s), W, H)

# ---------------- Fig 6：顯著判定流程（示意）----------------
def fig6():
    W, H = 1180, 420
    s = []
    s.append(T(12, 70, "雙向驗證：① 切群 → 看群是否對齊標籤　② 依標籤 → 看群內外是否顯著", 12.5, "ink", w=700))
    # 左：切群 k=2/3/6
    ax = 60
    s.append(R(ax, 92, 360, 250, "#fff", stroke="line", rx=10))
    s.append(T(ax+16, 116, "① 依組距離切群（k 掃描）", 12, "teal", w=800))
    for k,(kk,yy) in enumerate([("k=2",150),("k=3",206),("k=6",262)]):
        s.append(T(ax+24, yy, kk, 11, "ink", w=700))
        # mini bars representing clusters
        n = int(kk.split("=")[1])
        cw = min(40, 300/n)
        cols = ["germ","som","teal","gold","germL","somL"]
        for c in range(n):
            s.append(R(ax+80+c*(cw+4), yy-12, cw, 16, cols[c%6], rx=3))
    s.append(T(ax+24, 300, "→ 每種 k 都檢查群與標籤的關聯", 10, "soft"))
    s.append(arrow(ax+360+30, 217))
    # 右：4 標籤軸關聯（示意值）
    bx = 510
    s.append(R(bx, 92, 610, 250, "#fff", stroke="line", rx=10))
    s.append(T(bx+16, 116, "② 依標籤驗證群內外顯著（Fisher + Cramér's V）", 12, "gold", w=800))
    axes = [("Allele", 1.00, "germ", "分群跟著等位走"),
            ("HP", 0.98, "germ", "分群跟著單倍型走"),
            ("Strand", 0.00, "grey", "分群不跟股別走（好）"),
            ("Source", 1.00, "som", "分群跟著 tumor/normal")]
    by = 146; bw = 300
    for k,(lab,v,col,note) in enumerate(axes):
        yy = by + k*44
        s.append(T(bx+20, yy+15, lab, 11, "ink", w=700))
        s.append(R(bx+110, yy, bw, 20, "#eee", rx=4))
        s.append(R(bx+110, yy, max(bw*v,3), 20, col, rx=4))
        tag = "高顯著" if v>=0.5 else "不顯著"
        s.append(T(bx+110+bw+10, yy+15, f"{int(v*100)}% {tag}", 10, col if v>=0.5 else "soft", w=700))
        s.append(T(bx+110+bw+96, yy+15, note, 9.5, "soft"))
    # watermark big
    s.append(T(W/2, 384, "⚠ 示意 schematic — 100/98/0/100 是說明流程的舉例值，非真實位點結果", 12, "warn", a="middle", w=800))
    return wrap("圖 6 · 顯著判定流程（雙向驗證，示意）",
                "切群↔標籤雙向檢查；Strand 不顯著＝甲基分群不跟定序股別走", "\n".join(s), W, H)

if __name__ == "__main__":
    import os
    out = os.path.dirname(os.path.abspath(__file__))
    for name, fn in [("fig2_single_locus_pipeline", fig2),
                     ("fig5_axisC_vs_axisA", fig5),
                     ("fig6_significance_schematic", fig6)]:
        p = os.path.join(out, name + ".svg")
        open(p, "w", encoding="utf-8").write(fn())
        print("wrote", p)
