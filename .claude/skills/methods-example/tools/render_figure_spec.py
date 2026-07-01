#!/usr/bin/env python3
"""render_figure_spec.py — methods-example 物件組合渲染器（§13-A 構造防捏造）.

讀 figure_spec.json + 同目錄 data.json → 對每個 primitive 從 data_ref 注入真值 →
組合成單一 SVG（+ HTML 預覽 + provenance 稽核表）.

核心鐵則（與專案 scripts/fill_report.py 同模式）:
  * 凡 primitive 參數 *_ref 指向 data 的 `verified.*` 節點 → 該節點必須有非 null `value`,
    否則 **REFUSE**（exit 3）—— 數字物理上只能來自 data_ref, 無法手打捏造.
  * `verified.*` 節點若 value 故意為 null（如 max_abs_delta TODO）→ 渲染成「方法說明」而非假數字.
  * `schematic.* / synthetic` 區為示意（synthetic flag）, 打 watermark, 與真值物理隔離.

用法:
  python3 render_figure_spec.py <figure_spec.json> [-o out.svg] [--html] [--audit]
退碼: 0 OK / 2 spec 或 data 讀取錯 / 3 REFUSE（缺 verified 真值）.
"""
import argparse
import html
import json
import os
import sys


class Refuse(Exception):
    pass


# ---------- data ref resolution ----------
def resolve(data, dotted):
    """走訪 dotted path（支援 a.b.c）。缺 key → Refuse。"""
    node = data
    for key in dotted.split("."):
        if not isinstance(node, dict) or key not in node:
            raise Refuse(f"missing data path: {dotted!r} (broke at {key!r})")
        node = node[key]
    return node


def req_val(data, dotted):
    """verified 真值: 取 node.value, None → Refuse。回 (value, src)。"""
    node = resolve(data, dotted)
    if not (isinstance(node, dict) and "value" in node):
        raise Refuse(f"{dotted!r} is not a verified value-node (need .value + .src)")
    if node["value"] is None:
        raise Refuse(f"{dotted!r} value is null — 先 LOCK-AND-GATHER 落檔真值再渲染")
    return node["value"], node.get("src", "?")


def opt_val(data, dotted):
    """副鍵/可空真值: value=null → 回 (None, src 說明)，不 Refuse。"""
    node = resolve(data, dotted)
    if isinstance(node, dict) and "value" in node:
        return node["value"], node.get("src", "?")
    return None, "?"


PROV = []  # provenance 稽核: (metric, value, src)


def track(metric, value, src):
    PROV.append((metric, value, src))
    return value


def fmt(v):
    if isinstance(v, float):
        s = f"{v:.3f}".rstrip("0").rstrip(".")
        return s if s not in ("-0", "") else "0"
    return str(v)


# ---------- svg helpers ----------
def esc(s):
    return html.escape(str(s), quote=True)


def circle(cx, cy, r, filled, col_meth, col_unmeth):
    if filled:
        return f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{col_meth}"/>'
    return f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="#fff" stroke="{col_unmeth}" stroke-width="2"/>'


# ---------- primitive renderers (P1-P5; 回 (svg, height)) ----------
def p_read_cpg_matrix(data, prim, sh, y):
    cs = sh["color_scale"]
    groups = resolve(data, prim["params"]["groups_ref"])
    cpg_labels = resolve(data, prim["params"]["cpg_labels_ref"])
    snv_after = resolve(data, prim["params"]["snv_after_ref"])
    note = resolve(data, prim["params"]["snv_note_ref"]) if "snv_note_ref" in prim["params"] else \
        resolve(data, prim["params"].get("synthetic_note_ref", "schematic.note"))
    intu = resolve(data, prim["intuition_ref"]) if prim.get("intuition_ref") else None

    x0, col0, dx = 132, 150, 42
    ncol = len(cpg_labels)
    cols = [col0 + i * dx for i in range(ncol)]
    snv_x = (cols[snv_after - 1] + cols[snv_after]) // 2
    s = []
    yy = y
    if intu:
        s.append(f'<text x="10" y="{yy}" font-size="12.5" font-style="normal" fill="{cs["ink"]}">'
                 f'<tspan font-weight="700">直覺：</tspan>{esc(intu)}</text>')
        yy += 22
    # axis + CpG ticks
    s.append(f'<line x1="{col0}" y1="{yy+4}" x2="{cols[-1]}" y2="{yy+4}" stroke="#bbb"/>')
    for cx, lab in zip(cols, cpg_labels):
        s.append(f'<line x1="{cx}" y1="{yy}" x2="{cx}" y2="{yy+8}" stroke="#999"/>'
                 f'<text x="{cx}" y="{yy-2}" font-size="8" fill="#999" text-anchor="middle">{esc(lab)}</text>')
    # SNV marker
    s.append(f'<polygon points="{snv_x},{yy-12} {snv_x-6},{yy-22} {snv_x+6},{yy-22}" fill="{cs["meth"]}"/>'
             f'<text x="{snv_x}" y="{yy-26}" font-size="9" fill="{cs["meth"]}" text-anchor="middle" font-weight="700">G&gt;A</text>')
    yy += 18
    snv_top = yy - 4
    for g in groups:
        gcol = cs.get(g["color_key"], cs["ink"])
        s.append(f'<text x="{x0-6}" y="{yy+6}" font-size="10.5" font-weight="700" fill="{gcol}" text-anchor="end">{esc(g["label"])}</text>')
        s.append(f'<text x="{x0-6}" y="{yy+18}" font-size="8.5" fill="{cs["soft"]}" text-anchor="end">{esc(g["sublabel"])}</text>')
        for read in g["reads"]:
            cells = read.split()
            s.append(f'<line x1="{col0}" y1="{yy}" x2="{cols[-1]}" y2="{yy}" stroke="#cfcabd" stroke-width="6" stroke-linecap="round"/>')
            for cx, st in zip(cols, cells):
                s.append(circle(cx, yy, 5, st == "m", cs["meth"], cs["unmeth"]))
            yy += 16
        yy += 6
    # SNV vertical guide
    s.insert(0, f'<line x1="{snv_x}" y1="{snv_top}" x2="{snv_x}" y2="{yy-22}" stroke="{cs["meth"]}" stroke-dasharray="2 3" opacity=".45"/>')
    s.append(f'<text x="{col0}" y="{yy}" font-size="8.5" fill="{cs["synthetic"]}">（{esc(note)}）</text>')
    yy += 14
    return "\n".join(s), yy - y


def p_beta_bar(data, prim, sh, y):
    cs = sh["color_scale"]
    s = [f'<text x="10" y="{y}" font-size="12" font-weight="700" fill="{cs["ink"]}">{esc(prim.get("title",""))}</text>']
    yy = y + 14
    maxw = 300
    for b in prim["params"]["bars"]:
        v, src = req_val(data, b["value_ref"])
        track(b["label"], v, src)
        col = cs.get(b["color_key"], cs["ink"])
        w = int(maxw * max(0.0, min(1.0, float(v))))
        s.append(f'<text x="124" y="{yy+11}" font-size="10" fill="{col}" text-anchor="end">{esc(b["label"])}</text>')
        s.append(f'<rect x="132" y="{yy+2}" width="{w}" height="13" rx="3" fill="{col}" opacity=".55"/>'
                 f'<text x="{132+w+6}" y="{yy+12}" font-size="11" font-weight="800" fill="{col}">{fmt(v)}</text>')
        yy += 22
    return "\n".join(s), yy - y + 4


def p_dbeta_step(data, prim, sh, y):
    cs = sh["color_scale"]
    mean, msrc = req_val(data, prim["params"]["mean_ref"]); track("Δβ_mean", mean, msrc)
    n, nsrc = req_val(data, prim["params"]["n_ref"]); track("n_shared_cpg", n, nsrc)
    p, psrc = req_val(data, prim["params"]["p_ref"]); track("p_somatic", p, psrc)
    maxabs, _ = opt_val(data, prim["params"].get("maxabs_ref", "verified._none"))
    s = [f'<text x="10" y="{y}" font-size="12" font-weight="700" fill="{cs["ink"]}">{esc(prim.get("title",""))}</text>']
    yy = y + 18
    s.append(f'<rect x="10" y="{yy}" width="{sh.get("width",760)-20}" height="44" rx="6" fill="#F3F6F3" stroke="#cfe0d2"/>')
    s.append(f'<text x="22" y="{yy+19}" font-size="13" font-weight="800" fill="{cs["meth"]}">主鍵 Δβ = 平均(Δ) = {fmt(mean)}'
             f'<tspan font-size="10.5" font-weight="400" fill="#3a3a35">　paired Wilcoxon over ~{fmt(n)} CpG, p{"&lt;0.001" if float(p)==0.0 else "="+fmt(p)}</tspan></text>')
    if maxabs is None:
        s.append(f'<text x="22" y="{yy+37}" font-size="10.5" fill="#1f6b3f"><tspan font-weight="700">副鍵 max|Δ|</tspan>'
                 f'：方法同時保留 per-CpG 最大絕對差作副鍵，抓「交叉/位置型」（mean≈0 但 pattern 大）盲點（此位點 max|Δ| 未在 verified JSON，故只列方法，不顯示未驗證數字）。</text>')
    else:
        s.append(f'<text x="22" y="{yy+37}" font-size="10.5" fill="#1f6b3f"><tspan font-weight="700">副鍵 max|Δ|</tspan> = {fmt(maxabs)}（抓交叉/位置型盲點）</text>')
    yy += 54
    return "\n".join(s), yy - y


def p_normal_cis_triplet(data, prim, sh, y):
    cs = sh["color_scale"]
    nb, s0 = req_val(data, prim["params"]["normal_ref"]); track("normalHP1 β", nb, s0)
    t1, s1 = req_val(data, prim["params"]["tumorHP1_ref"]); track("tumorHP1 β", t1, s1)
    t11, s2 = req_val(data, prim["params"]["tumorHP11_ref"]); track("tumorHP1-1 β", t11, s2)
    dcis, s3 = req_val(data, prim["params"]["dcis_ref"]); track("d_cis", dcis, s3)
    pcis, _ = req_val(data, prim["params"]["pcis_ref"])
    ddr, s4 = req_val(data, prim["params"]["ddrift_ref"]); track("d_drift", ddr, s4)
    pdr, _ = req_val(data, prim["params"]["pdrift_ref"])
    dw, s5 = req_val(data, prim["params"]["dwithin_ref"]); track("d_within", dw, s5)
    pw, _ = req_val(data, prim["params"]["pwithin_ref"])
    s = [f'<text x="10" y="{y}" font-size="12" font-weight="700" fill="{cs["ink"]}">{esc(prim.get("title",""))}</text>']
    yy = y + 16
    maxw = 250
    rows = [("normal HP1", nb, cs["normal"]), ("tumor HP1", t1, cs["hp1"]), ("tumor HP1-1", t11, cs["hp11"])]
    for lab, v, col in rows:
        w = int(maxw * max(0.0, min(1.0, float(v))))
        s.append(f'<text x="124" y="{yy+11}" font-size="10" fill="{col}" text-anchor="end">{esc(lab)}</text>')
        s.append(f'<rect x="132" y="{yy+2}" width="{w}" height="13" rx="3" fill="{col}" opacity=".6"/>'
                 f'<text x="{132+w+6}" y="{yy+12}" font-size="10.5" font-weight="800" fill="{col}">{fmt(v)}</text>')
        yy += 20
    ax = 432
    s.append(f'<text x="{ax}" y="{y+30}" font-size="10.5" fill="#3a3a35"><tspan font-weight="700">d_drift</tspan> = {fmt(ddr)} (p={fmt(pdr)}, NS) ← 小</text>')
    s.append(f'<text x="{ax}" y="{y+50}" font-size="10.5" fill="{cs["meth"]}"><tspan font-weight="700">d_cis</tspan> = {fmt(dcis)} (p{"&lt;0.001" if float(pcis)==0.0 else "="+fmt(pcis)}) ← 大</text>')
    s.append(f'<text x="{ax}" y="{y+70}" font-size="10.5" font-weight="700" fill="{cs["synthetic"]}">⚠ 純 allele cis d_within = {fmt(dw)} (p={fmt(pw)}) 小但真</text>')
    return "\n".join(s), yy - y + 8


# ---------- P-READTRACK-LEGEND / hap_split_track (U1; deck 慣例：藍 germline / 橘 somatic) ----------
def p_readtrack_legend(data, prim, sh, y):
    cs = sh["color_scale"]
    items = resolve(data, prim["params"]["items_ref"])
    s = [f'<text x="10" y="{y+11}" font-size="10" font-weight="700" fill="{cs["soft"]}">圖例</text>']
    x = 56
    for it in items:
        col = cs.get(it["color_key"], cs["ink"])
        s.append(f'<circle cx="{x}" cy="{y+7}" r="5.5" fill="{col}"/>'
                 f'<text x="{x+11}" y="{y+11}" font-size="10.5" fill="{cs["ink"]}">{esc(it["label"])}</text>')
        x += 168
    return "\n".join(s), 24


def p_hap_split_track(data, prim, sh, y):
    cs = sh["color_scale"]
    detail = prim["params"].get("detail", "tick")  # tick(預設,隱藏 base 降載) / base(變異格才顯字母,強調哪種突變) / sequence(全序列,僅特例)
    germ = resolve(data, prim["params"]["germline_ref"])
    som = resolve(data, prim["params"]["somatic_ref"])
    cols = resolve(data, prim["params"]["cols_ref"])
    vcol = {"1": cs.get("hp1", "#2F5597"), "2": cs.get("hp2", "#8FAADC"),
            "s": cs.get("hp11", "#C55A11"), "e": cs.get("seqerr", "#7F7F7F")}
    x0, col0, dx = 150, 172, 40
    colx = [col0 + i * dx for i in range(len(cols))]
    s, st = [], {"yy": y}

    def draw_group(grp, layer_label):
        s.append(f'<text x="{x0-12}" y="{st["yy"]+2}" font-size="10" font-weight="700" fill="{cs["soft"]}" text-anchor="end">{esc(layer_label)}</text>')
        for g in grp:
            gcol = cs.get(g["color_key"], cs["ink"])
            s.append(f'<text x="{x0-12}" y="{st["yy"]+18}" font-size="10.5" font-weight="800" fill="{gcol}" text-anchor="end">{esc(g["label"])}</text>')
            if g.get("sublabel"):
                s.append(f'<text x="{x0-12}" y="{st["yy"]+29}" font-size="8" fill="{cs["soft"]}" text-anchor="end">{esc(g["sublabel"])}</text>')
            for read in g["reads"]:
                marks = (read.get("marks", "") or "").split()
                bases = (read.get("bases", "") or "").split()
                yy = st["yy"]
                x1, x2, tip = col0 - 14, colx[-1] + 14, colx[-1] + 24
                s.append(f'<polygon points="{x1},{yy+4} {x2},{yy+4} {tip},{yy+12} {x2},{yy+20} {x1},{yy+20}" fill="{gcol}" opacity=".14"/>')
                for i, cx in enumerate(colx):
                    m = marks[i] if i < len(marks) else "."
                    b = bases[i] if i < len(bases) else ""
                    if m != ".":  # 變異一律顯（位置 + 類別色）= L0 必要
                        r = 5.5 if detail == "tick" else 8
                        s.append(f'<circle cx="{cx}" cy="{yy+12}" r="{r}" fill="{vcol.get(m, cs["ink"])}"/>')
                        if detail in ("base", "sequence") and b:  # base 字母 = L2，僅強調「哪種突變」時
                            s.append(f'<text x="{cx}" y="{yy+16}" font-size="10" font-weight="800" fill="#fff" text-anchor="middle">{esc(b)}</text>')
                    elif detail == "sequence" and b:  # 背景 base = L3 噪音，預設隱藏；僅 full-sequence 顯（淡）
                        s.append(f'<text x="{cx}" y="{yy+16}" font-size="9" fill="#c4c2ba" text-anchor="middle">{esc(b)}</text>')
                st["yy"] += 24
            st["yy"] += 8

    draw_group(germ, "Germline")
    if som:  # 純 phasing 教學件(t1)可給空 somatic → 跳過分隔線+somatic 層
        s.append(f'<line x1="{x0-40}" y1="{st["yy"]}" x2="{colx[-1]+28}" y2="{st["yy"]}" stroke="{cs.get("line","#ccc")}" stroke-dasharray="4 3"/>')
        s.append(f'<text x="{x0-40}" y="{st["yy"]-3}" font-size="8.5" fill="{cs["soft"]}">↓ haplotag → 依 somatic 變異細分子標</text>')
        st["yy"] += 16
        draw_group(som, "Somatic")
    return "\n".join(s), st["yy"] - y + 4


# ---------- U2 igv_pileup（真實結果風格；read=lane rect，SNV=base 色直 tick）----------
def p_igv_pileup(data, prim, sh, y):
    cs = sh["color_scale"]
    lanes = resolve(data, prim["params"]["lanes_ref"])
    ncol = prim["params"].get("n_cols", 24)
    W = sh.get("width", 760)
    x0 = 150
    w = W - x0 - 20
    base_col = {"A": "#3a8a3a", "C": "#2D6CB5", "G": "#9a6a2a", "T": "#C0392B"}
    s, yy = [], y
    s.append(f'<text x="{x0-10}" y="{yy+9}" font-size="9" fill="{cs["soft"]}" text-anchor="end">coverage</text>')
    for i in range(ncol):
        h = 5 + (i * 5) % 13
        cx = x0 + i * (w / ncol)
        s.append(f'<rect x="{cx:.1f}" y="{yy+13-h}" width="{w/ncol-1:.1f}" height="{h}" fill="#cfcabd"/>')
    yy += 20
    locus_top, locus_x = yy, x0 + w * 0.5
    for lane in lanes:
        lcol = cs.get(lane["color_key"], "#bbb")
        nreads = lane.get("n_reads", 4)
        s.append(f'<text x="{x0-10}" y="{yy+10}" font-size="9.5" font-weight="700" fill="{lcol}" text-anchor="end">{esc(lane["label"])}</text>')
        for r in range(nreads):
            s.append(f'<rect x="{x0}" y="{yy+r*6}" width="{w:.1f}" height="4.5" rx="2" fill="{lcol}" opacity=".5"/>')
        for vc in lane.get("variants", []):
            vx = x0 + (vc["col"] / ncol) * w
            s.append(f'<rect x="{vx:.1f}" y="{yy}" width="2.6" height="{nreads*6-2}" fill="{base_col.get(vc.get("base","A"),"#444")}"/>')
        yy += nreads * 6 + 8
    s.insert(0, f'<line x1="{locus_x:.1f}" y1="{locus_top}" x2="{locus_x:.1f}" y2="{yy-8}" stroke="{cs.get("hp11","#C55A11")}" stroke-dasharray="2 3" opacity=".55"/>')
    return "\n".join(s), yy - y


# ---------- U3 cpg_beta_matrix（甲基化軸：cell = red>0.5 / blue<0.5 ramp）+ W 視窗 ----------
def p_cpg_beta_matrix(data, prim, sh, y):
    cs = sh["color_scale"]
    groups = resolve(data, prim["params"]["groups_ref"])
    ncpg = prim["params"].get("n_cpg", 8)
    win = prim["params"].get("window")
    col0, dx, x0 = 172, 30, 150
    colx = [col0 + i * dx for i in range(ncpg)]
    s, yy = [], y
    if win:
        wx1, wx2 = colx[win[0]] - 13, colx[win[1]] + 13
        s.append(f'<text x="{(wx1+wx2)/2:.0f}" y="{yy+8}" font-size="9" font-weight="700" fill="{cs["synthetic"]}" text-anchor="middle">W 視窗</text>')
        s.append(f'<path d="M{wx1},{yy+20} L{wx1},{yy+13} L{wx2},{yy+13} L{wx2},{yy+20}" fill="none" stroke="{cs["synthetic"]}" stroke-width="1.4"/>')
        yy += 22
    for cx in colx:
        s.append(f'<text x="{cx}" y="{yy+8}" font-size="8" fill="{cs["soft"]}" text-anchor="middle">CpG</text>')
    yy += 12
    for g in groups:
        gcol = cs.get(g.get("color_key", "hp1"), cs["ink"])
        s.append(f'<text x="{x0-10}" y="{yy+11}" font-size="10" font-weight="700" fill="{gcol}" text-anchor="end">{esc(g["label"])}</text>')
        for row in g["rows"]:
            for cx, v in zip(colx, [float(x) for x in row.split()]):
                col, op = (cs["meth"], (v - 0.5) * 1.5 + 0.3) if v >= 0.5 else (cs["unmeth"], (0.5 - v) * 1.5 + 0.3)
                s.append(f'<rect x="{cx-12}" y="{yy}" width="24" height="13" rx="2" fill="{col}" opacity="{min(op,1):.2f}"/>')
            yy += 16
        yy += 6
    return "\n".join(s), yy - y


# ---------- U5 grouped_bar（跨樣本 method 對照；zoomed 軸必標 caveat）----------
def p_grouped_bar(data, prim, sh, y):
    cs = sh["color_scale"]
    samples = resolve(data, prim["params"]["samples_ref"])
    methods = resolve(data, prim["params"]["methods_ref"])
    zoomed = prim["params"].get("zoomed", False)
    ymin, ymax = prim["params"].get("ymin", 0.0), prim["params"].get("ymax", 1.0)
    x0, W = 150, sh.get("width", 760)
    plot_w, plot_h = W - x0 - 30, 150
    base_y = y + 14 + plot_h
    s = []
    for gv in range(5):
        gy = base_y - gv / 4 * plot_h
        s.append(f'<line x1="{x0}" y1="{gy:.1f}" x2="{x0+plot_w}" y2="{gy:.1f}" stroke="#eee" stroke-dasharray="2 2"/>')
        s.append(f'<text x="{x0-6}" y="{gy+3:.1f}" font-size="8.5" fill="{cs["soft"]}" text-anchor="end">{ymin+gv/4*(ymax-ymin):.2f}</text>')
    group_w = plot_w / len(samples)
    bw = group_w * 0.66 / len(methods)
    for si, samp in enumerate(samples):
        gx = x0 + si * group_w + group_w * 0.17
        for mi, m in enumerate(methods):
            v = samp["values"].get(m["name"])
            if v is None:
                continue
            bh = (v - ymin) / (ymax - ymin) * plot_h
            bx = gx + mi * bw
            mcol = cs.get(m["color_key"], "#888")
            s.append(f'<rect x="{bx:.1f}" y="{base_y-bh:.1f}" width="{bw-2:.1f}" height="{bh:.1f}" fill="{mcol}"/>')
            s.append(f'<text x="{bx+bw/2:.1f}" y="{base_y-bh-2:.1f}" font-size="8" fill="{mcol}" text-anchor="middle">{v:.3f}</text>')
        s.append(f'<text x="{gx+group_w*0.33:.1f}" y="{base_y+13}" font-size="9" fill="{cs["ink"]}" text-anchor="middle">{esc(samp["name"])}</text>')
    yy = base_y + 20
    if zoomed:
        s.append(f'<text x="{x0}" y="{yy+9}" font-size="9" font-weight="700" fill="{cs["synthetic"]}">⚠ y 軸非零起點（{ymin}–{ymax}）放大差異，注意感知偏誤</text>')
        yy += 14
    return "\n".join(s), yy - y


# ---------- U4 facet_grid（小倍數 metric×sample；邊緣才標軸 + 共用 legend）----------
def p_facet_grid(data, prim, sh, y):
    cs = sh["color_scale"]
    rows = resolve(data, prim["params"]["rows_ref"])      # [{metric, cells:[{baseline,method}]}]
    cols = resolve(data, prim["params"]["col_labels_ref"])  # [sample]
    W, x0 = sh.get("width", 760), 150
    grid_w = W - x0 - 20
    cw, ph = grid_w / len(cols), 52
    s, yy = [], y
    for ci, c in enumerate(cols):
        s.append(f'<text x="{x0+ci*cw+cw/2:.0f}" y="{yy+8}" font-size="9" font-weight="700" fill="{cs["ink"]}" text-anchor="middle">{esc(c)}</text>')
    yy += 14
    for row in rows:
        s.append(f'<text x="{x0-8}" y="{yy+ph/2:.0f}" font-size="9.5" font-weight="700" fill="{cs["soft"]}" text-anchor="end">{esc(row["metric"])}</text>')
        for ci, cell in enumerate(row["cells"]):
            px, pw, by, ph_in = x0 + ci * cw + 6, cw - 12, yy + ph - 8, ph - 16
            s.append(f'<rect x="{px:.1f}" y="{yy+2}" width="{pw:.1f}" height="{ph-4}" fill="none" stroke="#eee"/>')
            for key, col, off in [("baseline", cs["unmeth"], 4), ("method", cs["hp11"], pw / 2 + 2)]:
                bh = float(cell.get(key, 0)) * ph_in
                s.append(f'<rect x="{px+off:.1f}" y="{by-bh:.1f}" width="{pw/2-6:.1f}" height="{bh:.1f}" fill="{col}"/>')
        yy += ph + 4
    s.append(f'<circle cx="{x0+10}" cy="{yy+6}" r="5" fill="{cs["unmeth"]}"/><text x="{x0+20}" y="{yy+10}" font-size="9.5" fill="{cs["ink"]}">baseline</text>')
    s.append(f'<circle cx="{x0+120}" cy="{yy+6}" r="5" fill="{cs["hp11"]}"/><text x="{x0+130}" y="{yy+10}" font-size="9.5" fill="{cs["ink"]}">method</text>')
    return "\n".join(s), yy + 18 - y


# ---------- U6 stacked_bar（組成隨 param sweep；schematic idealized）----------
def p_stacked_bar(data, prim, sh, y):
    cs = sh["color_scale"]
    xs = resolve(data, prim["params"]["x_values_ref"])
    series = resolve(data, prim["params"]["series_ref"])
    W, x0 = sh.get("width", 760), 150
    plot_w, plot_h = W - x0 - 120, 150
    base_y = y + 16 + plot_h
    n, bw = len(xs), (W - x0 - 120) / len(xs) * 0.6
    s = [f'<text x="{x0-6}" y="{y+22}" font-size="8.5" fill="{cs["soft"]}" text-anchor="end">100%</text>',
         f'<text x="{x0-6}" y="{base_y+3}" font-size="8.5" fill="{cs["soft"]}" text-anchor="end">0%</text>']
    for i, xv in enumerate(xs):
        bx, acc = x0 + i * (plot_w / n) + (plot_w / n - bw) / 2, 0.0
        for ser in series:
            h = float(ser["values"][i]) * plot_h
            s.append(f'<rect x="{bx:.1f}" y="{base_y-acc-h:.1f}" width="{bw:.1f}" height="{h:.1f}" fill="{cs.get(ser["color_key"],"#888")}"/>')
            acc += h
        s.append(f'<text x="{bx+bw/2:.1f}" y="{base_y+12}" font-size="8.5" fill="{cs["ink"]}" text-anchor="middle">{esc(str(xv))}</text>')
    s.append(f'<text x="{x0+plot_w/2:.0f}" y="{base_y+26}" font-size="9" fill="{cs["soft"]}" text-anchor="middle">{esc(prim["params"].get("x_label",""))}</text>')
    ly = y + 16
    for ser in series:
        s.append(f'<rect x="{x0+plot_w+16}" y="{ly}" width="10" height="10" fill="{cs.get(ser["color_key"],"#888")}"/>'
                 f'<text x="{x0+plot_w+30}" y="{ly+9}" font-size="9" fill="{cs["ink"]}">{esc(ser["name"])}</text>')
        ly += 16
    return "\n".join(s), base_y + 30 - y


# ---------- U7 loh_ideogram（單染色體 band + centromere + interval 軌）----------
def p_loh_ideogram(data, prim, sh, y):
    cs = sh["color_scale"]
    chrom = prim["params"].get("chrom", "chr13")
    length = float(prim["params"].get("length_mb", 115))
    cen = float(prim["params"].get("centromere_mb", 17))
    tracks = resolve(data, prim["params"]["tracks_ref"])
    W, x0 = sh.get("width", 760), 150
    w = W - x0 - 20

    def mb2x(mb):
        return x0 + float(mb) / length * w

    cx = mb2x(cen)
    s, yy = [], y
    s.append(f'<text x="{x0-10}" y="{yy+16}" font-size="10" font-weight="700" fill="{cs["ink"]}" text-anchor="end">{esc(chrom)}</text>')
    s.append(f'<rect x="{x0}" y="{yy+5}" width="{cx-x0:.1f}" height="15" rx="6" fill="#d8d5cc"/>')
    s.append(f'<rect x="{cx:.1f}" y="{yy+5}" width="{x0+w-cx:.1f}" height="15" rx="6" fill="#d8d5cc"/>')
    s.append(f'<circle cx="{cx:.1f}" cy="{yy+12}" r="4" fill="#a8a59c"/>')
    yy += 28
    for tr in tracks:
        s.append(f'<text x="{x0-10}" y="{yy+10}" font-size="9" fill="{cs["soft"]}" text-anchor="end">{esc(tr["label"])}</text>')
        s.append(f'<line x1="{x0}" y1="{yy+6}" x2="{x0+w}" y2="{yy+6}" stroke="#eee"/>')
        for iv in tr.get("intervals", []):
            ix1, ix2 = mb2x(iv["start"]), mb2x(iv["end"])
            s.append(f'<rect x="{ix1:.1f}" y="{yy}" width="{max(ix2-ix1,2):.1f}" height="12" rx="2" fill="{cs.get(iv.get("color_key","loh"),"#8E44AD")}" opacity=".8"/>')
        yy += 18
    return "\n".join(s), yy - y


# ---------- metro_map（地鐵圖：手放 grid 座標 + 線站序；粗彩線 + 轉乘站大圈 + 平行 offset）----------
def p_metro_map(data, prim, sh, y):
    cs = sh["color_scale"]
    stations = resolve(data, prim["params"]["stations_ref"])
    lines = resolve(data, prim["params"]["lines_ref"])
    dx = prim["params"].get("col_dx", 70)
    dy = prim["params"].get("row_dy", 56)
    x0, y0 = 80, y + 36
    off, nl = 5.0, len(lines)

    def sxy(sid):
        st = stations[sid]
        return x0 + st["col"] * dx, y0 + st["row"] * dy

    member = {}
    for li, ln in enumerate(lines):
        for sid in ln["stations"]:
            member.setdefault(sid, []).append(li)
    s = []
    # 線 legend（頂列）
    lx = x0
    for ln in lines:
        col = cs.get(ln["color_key"], "#888")
        dash = ' stroke-dasharray="6 4"' if ln.get("future") else ''
        s.append(f'<line x1="{lx}" y1="{y+12}" x2="{lx+22}" y2="{y+12}" stroke="{col}" stroke-width="5"{dash}/>')
        s.append(f'<text x="{lx+27}" y="{y+15}" font-size="9" fill="{cs["ink"]}">{esc(ln["name"])}</text>')
        lx += 32 + len(ln["name"]) * 9 + 12
    # 線（每線垂直 offset → 共段平行）
    for li, ln in enumerate(lines):
        col = cs.get(ln["color_key"], "#888")
        oy = (li - (nl - 1) / 2.0) * off
        pts = " ".join(f"{sxy(sid)[0]:.1f},{sxy(sid)[1]+oy:.1f}" for sid in ln["stations"])
        dash = ' stroke-dasharray="7 5"' if ln.get("future") else ''
        s.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="6" stroke-linejoin="round" stroke-linecap="round" opacity="{0.5 if ln.get("future") else 0.92}"{dash}/>')
    # 站 + label
    maxy = y0
    for sid, st in stations.items():
        x, yy_ = sxy(sid)
        maxy = max(maxy, yy_)
        lis = member.get(sid, [])
        if len(lis) >= 2:  # 轉乘站 = 大白圈黑框
            s.append(f'<circle cx="{x:.1f}" cy="{yy_:.1f}" r="9" fill="#fff" stroke="{cs["ink"]}" stroke-width="2.5"/>')
            sy = yy_
        else:
            col = cs.get(lines[lis[0]]["color_key"], "#888") if lis else "#888"
            sy = yy_ + ((lis[0] - (nl - 1) / 2.0) * off if lis else 0)
            s.append(f'<circle cx="{x:.1f}" cy="{sy:.1f}" r="5.5" fill="#fff" stroke="{col}" stroke-width="2.5"/>')
        lab = st.get("label", "")
        if st.get("label_pos") == "right":
            s.append(f'<text x="{x+12:.1f}" y="{sy+3:.1f}" font-size="8.5" font-weight="700" fill="{cs["ink"]}">{esc(lab)}</text>')
        else:  # backbone 站 label 旋轉 -38°（EPI2ME 慣例，避重疊）
            s.append(f'<text transform="rotate(-38 {x:.1f} {yy_-13:.1f})" x="{x:.1f}" y="{yy_-13:.1f}" font-size="8.5" font-weight="700" fill="{cs["ink"]}" text-anchor="start">{esc(lab)}</text>')
    return "\n".join(s), maxy - y + 40


# ---------- pipeline_diagram（EPI2ME drawio 風 banded swimlane + scope 註解區塊）----------
def p_pipeline_diagram(data, prim, sh, y):
    cs = sh["color_scale"]
    bands = resolve(data, prim["params"]["bands_ref"])
    W = sh.get("width", 1100)
    x0, lblw = 12, 184
    box_x0 = x0 + lblw + 8
    avail = W - box_x0 - 16
    box_h, gap = 34, 26
    s, yy = [], y
    for bi, band in enumerate(bands):
        bcol = cs.get(band.get("color_key", "normal"), "#888")
        lblcol = bcol if band.get("color_key") not in (None, "normal") else cs["ink"]
        stages = band["stages"]
        n = len(stages)
        bw = max(58, min(150, (avail - gap * (n - 1)) / max(n, 1)))
        band_h = box_h + 30
        # scope block（範圍註解意義區塊）= 半透明背景 + 左側 label + note
        s.append(f'<rect x="{x0}" y="{yy}" width="{W-2*x0}" height="{band_h}" rx="9" fill="{bcol}" opacity=".07"/>')
        s.append(f'<text x="{x0+10}" y="{yy+20}" font-size="11.5" font-weight="800" fill="{lblcol}">{esc(band["label"])}</text>')
        if band.get("note"):
            s.append(f'<text x="{x0+10}" y="{yy+37}" font-size="8.5" fill="{cs["soft"]}">{esc(band["note"])}</text>')
        cy = yy + band_h / 2
        fsz = 9.5 if bw >= 92 else (8.5 if bw >= 72 else 8)
        bx, prev_right = box_x0, None
        for st in stages:
            future = st.get("status") == "future"
            dash = ' stroke-dasharray="4 3"' if future else ''
            inkc = cs["soft"] if future else cs["ink"]
            s.append(f'<rect x="{bx:.1f}" y="{cy-box_h/2:.1f}" width="{bw:.1f}" height="{box_h}" rx="6" fill="#fff" stroke="{bcol}" stroke-width="1.4"{dash}/>')
            sub = st.get("sub")
            ty = (cy - 1) if sub else (cy + 3)
            s.append(f'<text x="{bx+bw/2:.1f}" y="{ty:.1f}" font-size="{fsz}" font-weight="700" fill="{inkc}" text-anchor="middle">{esc(st["label"])}</text>')
            if sub:
                s.append(f'<text x="{bx+bw/2:.1f}" y="{cy+11:.1f}" font-size="8" fill="{cs["soft"]}" text-anchor="middle">{esc(sub)}</text>')
            if prev_right is not None:
                s.append(f'<line x1="{prev_right:.1f}" y1="{cy:.1f}" x2="{bx-4:.1f}" y2="{cy:.1f}" stroke="{bcol}" stroke-width="2.5"/>')
                s.append(f'<path d="M{bx-4:.1f},{cy-3.5:.1f} L{bx:.1f},{cy:.1f} L{bx-4:.1f},{cy+3.5:.1f}" fill="{bcol}"/>')
            prev_right = bx + bw
            bx += bw + gap
        yy += band_h + 22
        if bi < len(bands) - 1:  # 向下 connector 到下一 band
            dx = box_x0 + 28
            s.append(f'<line x1="{dx}" y1="{yy-22}" x2="{dx}" y2="{yy-5}" stroke="{cs["soft"]}" stroke-width="2"/>')
            s.append(f'<path d="M{dx-3.5},{yy-9} L{dx},{yy-5} L{dx+3.5},{yy-9}" fill="{cs["soft"]}"/>')
    return "\n".join(s), yy - y + 2


# ---------- B1 gene_model_track（exon/intron/UTR + strand）----------
def p_gene_model_track(data, prim, sh, y):
    cs = sh["color_scale"]
    feats = resolve(data, prim["params"]["features_ref"])
    length = float(prim["params"].get("length", 100))
    strand = prim["params"].get("strand", "+")
    W, x0 = sh.get("width", 760), 150
    w = W - x0 - 30

    def sx(p):
        return x0 + float(p) / length * w

    cy = y + 16
    s = [f'<text x="{x0-10}" y="{cy+4}" font-size="10" font-weight="700" fill="{cs["ink"]}" text-anchor="end">{esc(prim["params"].get("label","gene"))}</text>',
         f'<line x1="{x0}" y1="{cy}" x2="{x0+w}" y2="{cy}" stroke="{cs.get("ink","#444")}" stroke-width="1.2"/>']
    for i in range(1, 6):
        ax, d = x0 + i * w / 6, (4 if strand == "+" else -4)
        s.append(f'<path d="M{ax-d:.1f},{cy-3} L{ax+d:.1f},{cy} L{ax-d:.1f},{cy+3}" fill="none" stroke="#999" stroke-width="1"/>')
    for f in feats:
        x1, x2 = sx(f["start"]), sx(f["end"])
        if f["type"] == "exon":
            s.append(f'<rect x="{x1:.1f}" y="{cy-8}" width="{max(x2-x1,2):.1f}" height="16" rx="1.5" fill="{cs.get("hp1","#2F5597")}"/>')
        elif f["type"] == "utr":
            s.append(f'<rect x="{x1:.1f}" y="{cy-4}" width="{max(x2-x1,2):.1f}" height="8" rx="1" fill="{cs.get("hp2","#8FAADC")}"/>')
    s.append(f'<text x="{x0+w}" y="{cy+24}" font-size="8.5" fill="{cs["soft"]}" text-anchor="end">exon=實心框 · UTR=細框 · intron=線 · 箭頭=strand ({strand})</text>')
    return "\n".join(s), 42


# ---------- B2 coverage_track（depth area + ploidy 線）----------
def p_coverage_track(data, prim, sh, y):
    cs = sh["color_scale"]
    depths = [float(v) for v in resolve(data, prim["params"]["depths_ref"])]
    ploidy = prim["params"].get("ploidy_line")
    W, x0 = sh.get("width", 760), 150
    w, h = W - x0 - 30, 70
    base_y, mx, n = y + 14 + h, max(depths) or 1, len(depths)
    bw = w / n
    s = [f'<text x="{x0-10}" y="{y+22}" font-size="10" font-weight="700" fill="{cs["ink"]}" text-anchor="end">{esc(prim["params"].get("label","coverage"))}</text>',
         f'<text x="{x0-10}" y="{base_y+3}" font-size="8" fill="{cs["soft"]}" text-anchor="end">0</text>',
         f'<text x="{x0-10}" y="{y+18}" font-size="8" fill="{cs["soft"]}" text-anchor="end">{int(mx)}×</text>']
    for i, d in enumerate(depths):
        bh = d / mx * h
        s.append(f'<rect x="{x0+i*bw:.1f}" y="{base_y-bh:.1f}" width="{bw:.1f}" height="{bh:.1f}" fill="{cs.get("normal","#999")}" opacity=".55"/>')
    if ploidy is not None:
        py = base_y - float(ploidy) / mx * h
        s.append(f'<line x1="{x0}" y1="{py:.1f}" x2="{x0+w}" y2="{py:.1f}" stroke="{cs.get("hp11","#C55A11")}" stroke-dasharray="4 3"/>')
        s.append(f'<text x="{x0+w}" y="{py-2:.1f}" font-size="8" fill="{cs.get("hp11","#C55A11")}" text-anchor="end">expected ploidy</text>')
    return "\n".join(s), h + 18


# ---------- B3 methyl_lollipop_track（單軌 5mC，filled/open 形狀冗餘）----------
def p_methyl_lollipop_track(data, prim, sh, y):
    cs = sh["color_scale"]
    sites = resolve(data, prim["params"]["sites_ref"])
    length = float(prim["params"].get("length", 100))
    W, x0 = sh.get("width", 760), 150
    w, base_y = W - x0 - 30, y + 40

    def sx(p):
        return x0 + float(p) / length * w

    s = [f'<text x="{x0-10}" y="{base_y+4}" font-size="10" font-weight="700" fill="{cs["ink"]}" text-anchor="end">{esc(prim["params"].get("label","5mC"))}</text>',
         f'<line x1="{x0}" y1="{base_y}" x2="{x0+w}" y2="{base_y}" stroke="#cfcabd" stroke-width="2"/>']
    for st in sites:
        px, top = sx(st["pos"]), base_y - 22
        s.append(f'<line x1="{px:.1f}" y1="{base_y}" x2="{px:.1f}" y2="{top}" stroke="#bbb" stroke-width="1"/>')
        if st.get("meth"):
            s.append(f'<circle cx="{px:.1f}" cy="{top}" r="5.5" fill="{cs["meth"]}"/>')
        else:
            s.append(f'<circle cx="{px:.1f}" cy="{top}" r="5.5" fill="#fff" stroke="{cs["unmeth"]}" stroke-width="2"/>')
    s.append(f'<text x="{x0+w}" y="{base_y+18}" font-size="8.5" fill="{cs["soft"]}" text-anchor="end">實心=甲基化 / 空心=未甲基化（形狀冗餘，不單靠顏色）</text>')
    return "\n".join(s), 62


# ---------- B7 variant_class_legendcard（germline/somatic 概念圖卡）----------
def p_variant_class_legendcard(data, prim, sh, y):
    cs = sh["color_scale"]
    items = resolve(data, prim["params"]["items_ref"])
    W, x0 = sh.get("width", 760), 150
    s = [f'<rect x="{x0}" y="{y}" width="{W-x0-30}" height="{len(items)*30+14}" rx="8" fill="#fff" stroke="{cs.get("line","#E4E2DA")}"/>']
    yy = y + 24
    for it in items:
        col = cs.get(it["color_key"], cs["ink"])
        s.append(f'<circle cx="{x0+22}" cy="{yy-4}" r="7" fill="{col}"/>')
        s.append(f'<text x="{x0+40}" y="{yy}" font-size="11.5" font-weight="700" fill="{col}">{esc(it["label"])}</text>')
        s.append(f'<text x="{x0+150}" y="{yy}" font-size="10.5" fill="{cs["ink"]}">{esc(it.get("desc",""))}</text>')
        yy += 30
    return "\n".join(s), yy - y + 6


RENDERERS = {
    "metro_map": p_metro_map,
    "pipeline_diagram": p_pipeline_diagram,
    "hap_split_track": p_hap_split_track,
    "readtrack_legend": p_readtrack_legend,
    "gene_model_track": p_gene_model_track,
    "coverage_track": p_coverage_track,
    "methyl_lollipop_track": p_methyl_lollipop_track,
    "variant_class_legendcard": p_variant_class_legendcard,
    "igv_pileup": p_igv_pileup,
    "cpg_beta_matrix": p_cpg_beta_matrix,
    "grouped_bar": p_grouped_bar,
    "facet_grid": p_facet_grid,
    "stacked_bar": p_stacked_bar,
    "loh_ideogram": p_loh_ideogram,
    "read_cpg_matrix": p_read_cpg_matrix,
    "beta_bar": p_beta_bar,
    "dbeta_step": p_dbeta_step,
    "normal_cis_triplet": p_normal_cis_triplet,
}


def render(spec, data):
    sh = spec["shared"]
    cs = sh["color_scale"]
    W = sh.get("width", 760)
    frags, y = [], 96
    for prim in spec["primitives"]:
        fn = RENDERERS.get(prim["type"])
        if fn is None:
            raise Refuse(f"primitive type {prim['type']!r} 尚未實作（C1 pilot 覆蓋 P1-P5；P6/P7 = C2/C3 待擴）")
        frag, h = fn(data, prim, sh, y + 18)
        frags.append(frag)
        y += h + 34
        frags.append(f'<line x1="10" y1="{y-18}" x2="{W-10}" y2="{y-18}" stroke="#eee"/>')

    # header
    title = esc(spec.get("title", ""))
    sub = esc(spec.get("subtitle", ""))
    head = [
        f'<text x="10" y="26" font-size="16" font-weight="800" fill="{cs["ink"]}">{title}</text>',
        f'<text x="10" y="46" font-size="11.5" fill="{cs["soft"]}">{sub}</text>',
    ]
    if spec.get("annotations", {}).get("methyl_legend", True):  # 甲基化圖預設顯示 5mC legend；haplotype 圖關閉用 readtrack_legend
        head.append(
            f'<g transform="translate(10,58)">'
            f'<circle cx="8" cy="8" r="5" fill="{cs["meth"]}"/><text x="18" y="12" font-size="10.5" fill="{cs["ink"]}">5mC high（甲基）</text>'
            f'<circle cx="150" cy="8" r="5" fill="#fff" stroke="{cs["unmeth"]}" stroke-width="2"/><text x="160" y="12" font-size="10.5" fill="{cs["ink"]}">5mC low（去甲基）</text>'
            f'<text x="300" y="12" font-size="10.5" fill="{cs["ink"]}">β＝該群在某 CpG 的甲基比例</text></g>')
    # annotations footer
    ann = spec.get("annotations", {})
    foot = []
    fy = y + 4
    if ann.get("label_flip_ref"):
        lf = esc(resolve(data, ann["label_flip_ref"]))
        foot.append(f'<text x="10" y="{fy}" font-size="10.5" font-weight="700" fill="#7a4">▣ label-flip：{lf}</text>')
        fy += 18
    if ann.get("single_sample_caveat_ref"):
        sc = esc(resolve(data, ann["single_sample_caveat_ref"]))
        foot.append(f'<text x="10" y="{fy}" font-size="10.5" fill="{cs["synthetic"]}">⚠ 誠實標註：{sc}</text>')
        fy += 18
    if ann.get("process_schematic"):
        foot.append(f'<text x="10" y="{fy}" font-size="9.5" fill="{cs["soft"]}">※ 流程 schematic：模組名與常數取自 ISM 源碼（附 file:line，可 grep），本圖無分析數值。</text>')
    else:
        foot.append(f'<text x="10" y="{fy}" font-size="9.5" fill="{cs["soft"]}">※ 聚合數字皆 verified 真值（{len(PROV)} 項，附 src 可 grep）；read 甲基點為示意（synthetic）。</text>')
    fy += 16
    if ann.get("synthetic_watermark"):
        foot.append(f'<text x="{W-12}" y="{fy}" font-size="8.5" fill="#ccc" text-anchor="end">schematic reads = synthetic / 示意</text>')
        fy += 8

    H = fy + 12
    body = "\n".join(head + frags + foot)
    return (f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" role="img" '
            f'aria-label="{title}"><title>{title}</title><desc>{sub}</desc>'
            # CJK font-family on all <text> → 中文 renders correctly in BROWSER (真交付物；瀏覽器逐字 cascade).
            # ⚠ cairosvg raster 預覽的豆腐是系統缺字型問題非本碼可修：cairosvg 2.8.2 只取 stack 第一個名字、無
            # per-glyph fallback，且本機 fontconfig 把未裝的 Noto 回退到 DejaVu(無CJK)；Droid Sans Fallback 有
            # CJK 但缺 latin → 混排無單字型解。乾淨 raster 需 `apt install fonts-noto-cjk`（有 latin+CJK），或用瀏覽器開 .html。
            f'<style>text{{font-family:"Noto Sans CJK TC","Source Han Sans TW","Droid Sans Fallback","WenQuanYi Zen Hei",system-ui,sans-serif}}</style>'
            f'<rect width="{W}" height="{H}" fill="#FAF9F5"/>\n{body}\n</svg>')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("-o", "--out")
    ap.add_argument("--html", action="store_true", help="同時輸出 .html 預覽 wrapper")
    ap.add_argument("--audit", action="store_true", help="印 provenance 稽核表")
    a = ap.parse_args()
    try:
        spec = json.load(open(a.spec, encoding="utf-8"))
        ddir = os.path.dirname(os.path.abspath(a.spec))
        data = json.load(open(os.path.join(ddir, spec["data_file"]), encoding="utf-8"))
    except Exception as e:
        print(f"[ERR] 讀 spec/data 失敗: {e}", file=sys.stderr); sys.exit(2)
    try:
        svg = render(spec, data)
    except Refuse as e:
        print(f"[REFUSE] {e}\n  → 補齊 verified 真值（落檔 + Read 回）再渲染；不可手打捏造數字（§13-A）。", file=sys.stderr)
        sys.exit(3)
    out = a.out or os.path.splitext(a.spec)[0].replace("figure_spec", "method_explainer") + ".svg"
    open(out, "w", encoding="utf-8").write(svg)
    print(f"[OK] SVG → {out}  ({len(PROV)} verified 真值注入)")
    if a.html:
        hp = out[:-4] + ".html"
        open(hp, "w", encoding="utf-8").write(
            f'<!DOCTYPE html><meta charset="utf-8"><title>{esc(spec.get("title",""))}</title>'
            f'<body style="background:#FAF9F5;font-family:system-ui;margin:0;padding:20px">{svg}</body>')
        print(f"[OK] HTML → {hp}")
    if a.audit:
        print("\n=== provenance 稽核（metric → value → src）===")
        for m, v, s in PROV:
            print(f"  {m:18s} = {fmt(v):>8s}   ⇐ {s}")


if __name__ == "__main__":
    main()
