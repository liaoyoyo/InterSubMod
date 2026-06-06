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


RENDERERS = {
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
        f'<g transform="translate(10,58)">'
        f'<circle cx="8" cy="8" r="5" fill="{cs["meth"]}"/><text x="18" y="12" font-size="10.5" fill="{cs["ink"]}">5mC high（甲基）</text>'
        f'<circle cx="150" cy="8" r="5" fill="#fff" stroke="{cs["unmeth"]}" stroke-width="2"/><text x="160" y="12" font-size="10.5" fill="{cs["ink"]}">5mC low（去甲基）</text>'
        f'<text x="300" y="12" font-size="10.5" fill="{cs["ink"]}">β＝該群在某 CpG 的甲基比例</text></g>',
    ]
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
    foot.append(f'<text x="10" y="{fy}" font-size="9.5" fill="{cs["soft"]}">※ 聚合數字皆 verified 真值（{len(PROV)} 項，附 src 可 grep）；read 甲基點為示意（synthetic）。</text>')
    fy += 16
    if ann.get("synthetic_watermark"):
        foot.append(f'<text x="{W-12}" y="{fy}" font-size="8.5" fill="#ccc" text-anchor="end">schematic reads = synthetic / 示意</text>')
        fy += 8

    H = fy + 12
    body = "\n".join(head + frags + foot)
    return (f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" role="img" '
            f'aria-label="{title}"><title>{title}</title><desc>{sub}</desc>'
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
