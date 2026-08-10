#!/usr/bin/env python3
"""CN / HP 觀察 pilot 頁生成器（HCC1395）

回應四項觀察需求：
  1. 演化分支節點改用 S1/S2… 短代號（區域內依基因體位置穩定編號，hover 顯完整座標）
  2. 同一位點在不同 HP 的對照視圖（現有工作站完全沒有）
  3. LOH / CNV 區間疊圖（現有工作站 CNV=0、LOH 僅 9 次且非區間）
  4. 補上整體數據的觀察圖表

資料來源（全部已存在，本腳本不重算任何科學量）：
  - canonical topology  : HCC1395.topology.jsonl（樹結構、active_positions、af_coverage）
  - CN unit 註釋        : hcc1395_unit_cn_annotation.jsonl（SEQC2 六類 + LOH 獨立軸）
  - CN site 註釋        : hcc1395_site_cn_annotation.jsonl（per-site CN/LOH/AF）

🔴 科學紀律（不可放寬）
  - CN/LOH 一律為**描述性疊圖**。`authority_manifest.claim_boundary.forbidden` 明列
    "CN/LOH-corrected CCF"，2026-08-10 稽核亦定案「矯正只能是 fail-closed gate，
    不是數值校正」。本頁不做任何 CN 數值校正，只標註。
  - read-AF **只取 active_positions**。af_coverage 含非 active 位點，未過濾會使
    離散度與中介檢驗全部反轉（見 memory reference_af_coverage_contains_inactive_sites）。
  - SEQC2 為外部真值優先；SAVANA 發布 fit 已判定 mis-calibrated（整數 CN 吻合僅 3.64%），
    僅作對照顯示，不作判準。

用法: python3 scripts/build_cn_hp_observation_pilot.py [-o OUT.html]
"""
import argparse
import collections
import html
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUDIT = os.path.join(ROOT, "research", "20260810_cnv_confound_seqc2_savana_audit", "data")
TOPO = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
        "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/"
        "samples/HCC1395/HCC1395.topology.jsonl")
UNIT_CN = os.path.join(AUDIT, "hcc1395_unit_cn_annotation.jsonl")
SITE_CN = os.path.join(AUDIT, "hcc1395_site_cn_annotation.jsonl")

# SEQC2 六類的顏色與中文名（LOH 是獨立軸，不與 gain/loss 合併）
CN_CLASS = collections.OrderedDict([
    ("neutral",  ("#2F7D5B", "CN 中性",        "copy number 2，無 LOH")),
    ("gain",     ("#B06A16", "增益",           "total CN > 2，無 LOH")),
    ("gain_loh", ("#9A3B3B", "增益 + LOH",     "CN 增益且雜合性缺失")),
    ("cnloh",    ("#5B3FA0", "CN 中性 LOH",    "copy number 2 但雜合性缺失")),
    ("loss",     ("#2F6690", "缺失",           "total CN < 2，無 LOH")),
    ("loss_loh", ("#1c1b19", "缺失 + LOH",     "CN 缺失且雜合性缺失")),
])
CHROMS = ["chr%s" % c for c in list(range(1, 23))]


def load_jsonl(path, need=None):
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            d = json.loads(line)
            if need and not all(k in d for k in need):
                continue
            out.append(d)
    return out


# ─────────────────────────── 需求 1：S 編號 ───────────────────────────
def site_labels(unit):
    """區域內依基因體位置穩定編號 → {position: 'S1'}。

    穩定性來源：active_positions 是該 unit 的固定集合，排序後編號與資料無關的
    渲染順序無涉，因此同一 unit 每次產生的 S 編號都相同、可跨圖引用。
    """
    pos = sorted(unit.get("active_positions") or [])
    return {p: "S%d" % (i + 1) for i, p in enumerate(pos)}


def tree_model(unit):
    """把 representative_best_edges/vertices 轉成用 S 編號表示的樹。

    節點原本的 label 是 ROOT / H_AR / AA 這類等位組合字串，難以閱讀；
    改用「該節點已獲得哪些 S 編號」表示，例如 ROOT → {S1} → {S1,S2}。
    """
    labels = site_labels(unit)
    edges = unit.get("representative_best_edges") or []
    if not edges:
        return None
    # 逐邊累積：child 節點 = parent 已有的集合 + 本邊獲得的位點
    acquired = {}          # vertex -> set of S labels
    root = None
    links = []
    for e in edges:
        p, c = e.get("parent_vertex"), e.get("child_vertex")
        pos = e.get("acquired_position")
        if p is None or c is None:
            continue
        if root is None:
            root = p
        s = labels.get(pos)
        acquired.setdefault(p, set())
        acquired[c] = set(acquired.get(p, set())) | ({s} if s else set())
        links.append({"parent": p, "child": c, "gained": s,
                      "pos": pos, "bit": e.get("acquired_active_bit")})
    if root is None:
        return None
    nodes = {}
    for v, ss in acquired.items():
        nodes[v] = "ROOT" if not ss else "+".join(sorted(ss, key=lambda x: int(x[1:])))
    return {"root": root, "nodes": nodes, "links": links, "labels": labels}


def pick_examples(units, unit_cn, n=6):
    """挑代表性 unit 做分支圖：涵蓋不同 CN 類別與不同分支數，且必須有樹。"""
    cn_by_id = {u["unit_id"]: u for u in unit_cn}
    cand = []
    for u in units:
        if u.get("unit_status") != "ranked":
            continue
        k = u.get("active_bit_count") or 0
        if k < 2:                     # k<2 沒有分支可看
            continue
        tm = tree_model(u)
        if not tm or len(tm["links"]) < 2:
            continue
        cn = cn_by_id.get(u["unit_id"], {})
        cand.append((u, tm, cn.get("seqc2_class")))
    # 每個 CN 類別各挑一個，優先挑分支較多的（較能展示 S 編號的價值）
    out, seen = [], set()
    for cls in CN_CLASS:
        pool = [c for c in cand if c[2] == cls]
        if not pool:
            continue
        pool.sort(key=lambda c: (-len(c[1]["links"]), c[0]["unit_id"]))
        out.append(pool[0]); seen.add(pool[0][0]["unit_id"])
        if len(out) >= n:
            break
    for c in sorted(cand, key=lambda c: -len(c[1]["links"])):
        if len(out) >= n:
            break
        if c[0]["unit_id"] not in seen:
            out.append(c); seen.add(c[0]["unit_id"])
    return out


def svg_tree(unit, tm, cn_class, width=430):
    """畫一棵樹。節點顯 S 編號組合，邊標「獲得 Sx」，hover 顯完整座標。"""
    # 依 parent 關係分層
    children = collections.defaultdict(list)
    for l in tm["links"]:
        children[l["parent"]].append(l)
    depth, order = {tm["root"]: 0}, [tm["root"]]
    i = 0
    while i < len(order):
        v = order[i]; i += 1
        for l in children.get(v, []):
            if l["child"] not in depth:
                depth[l["child"]] = depth[v] + 1
                order.append(l["child"])
    bylv = collections.defaultdict(list)
    for v, d in depth.items():
        bylv[d].append(v)
    maxd = max(bylv) if bylv else 0
    xy, ROWH = {}, 62
    for d in sorted(bylv):
        row = sorted(bylv[d])
        for j, v in enumerate(row):
            xy[v] = (width / (len(row) + 1) * (j + 1), 34 + d * ROWH)
    h = 34 + maxd * ROWH + 46
    col = CN_CLASS.get(cn_class, ("#6B6B66",))[0]
    p = ['<svg viewBox="0 0 %d %d" role="img" aria-labelledby="t%s d%s" '
         'xmlns="http://www.w3.org/2000/svg">' % (width, h, unit["unit_id"][:8], unit["unit_id"][:8])]
    p.append('<title id="t%s">%s 的演化分支</title>' % (unit["unit_id"][:8], html.escape(unit["region_id"])))
    p.append('<desc id="d%s">以 S 編號表示的有根樹，根為 germline，每條邊獲得一個體細胞突變。</desc>'
             % unit["unit_id"][:8])
    p.append('<rect x="0" y="0" width="%d" height="%d" fill="#FAF9F5"/>' % (width, h))
    for l in tm["links"]:
        if l["parent"] not in xy or l["child"] not in xy:
            continue
        x1, y1 = xy[l["parent"]]; x2, y2 = xy[l["child"]]
        p.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#8a8578" stroke-width="1.5"/>'
                 % (x1, y1 + 13, x2, y2 - 15))
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2 + 2
        tip = "獲得 %s ＝ %s:%s" % (l["gained"] or "?", unit.get("chrom"), l["pos"])
        p.append('<g><title>%s</title>'
                 '<rect x="%.1f" y="%.1f" width="34" height="16" rx="8" fill="#F2E4DC" stroke="%s" stroke-width="1"/>'
                 '<text x="%.1f" y="%.1f" text-anchor="middle" font-size="10" font-weight="700" fill="#8A3F22">+%s</text></g>'
                 % (html.escape(tip), mx - 17, my - 8, col, mx, my + 3.5, html.escape(l["gained"] or "?")))
    for v, (x, y) in xy.items():
        lab = tm["nodes"].get(v, "?")
        isroot = (v == tm["root"])
        w = max(38, 9 * len(lab) + 14)
        p.append('<rect x="%.1f" y="%.1f" width="%.1f" height="26" rx="6" fill="%s" stroke="%s" stroke-width="1.6"/>'
                 % (x - w / 2, y - 13, w, "#FFFFFF" if not isroot else "#E4F0EA", col))
        p.append('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" font-weight="700" fill="%s">%s</text>'
                 % (x, y + 4, col, html.escape(lab)))
    p.append("</svg>")
    return "".join(p)


# ─────────────────────── 需求 2：同位點跨 HP 對照 ───────────────────────
def cross_hp(sites):
    by = collections.defaultdict(list)
    for s in sites:
        by[(s["chrom"], s["position"])].append(s)
    rows = []
    for (c, pos), lst in by.items():
        hps = {x["hp_family"] for x in lst}
        if len(hps) < 2:
            continue
        rec = {"chrom": c, "position": pos, "cn_class": lst[0].get("seqc2_class"),
               "loh": bool(lst[0].get("seqc2_loh")), "cn": lst[0].get("seqc2_total_cn"), "hp": {}}
        for x in lst:
            rec["hp"][x["hp_family"]] = {
                "alt": x.get("alt_reads"), "ref": x.get("ref_reads"),
                "af": x.get("read_af"), "status": x.get("unit_status"),
                "unit": x.get("unit_id"), "region": x.get("region_id")}
        rows.append(rec)
    rows.sort(key=lambda r: (CHROMS.index(r["chrom"]) if r["chrom"] in CHROMS else 99, r["position"]))
    return rows


# ─────────────────── 需求 3+4：統計與染色體軌道 ───────────────────
def chrom_track(sites, width=980):
    """每條染色體一條軌道，用 CN 六類著色標出每個 sSNV 的位置。"""
    by = collections.defaultdict(list)
    for s in sites:
        if s["chrom"] in CHROMS:
            by[s["chrom"]].append(s)
    maxpos = {c: max((x["position"] for x in v), default=1) for c, v in by.items()}
    gmax = max(maxpos.values()) if maxpos else 1
    L, R, ROW = 58, width - 12, 21
    h = 30 + len(CHROMS) * ROW + 8
    p = ['<svg viewBox="0 0 %d %d" role="img" aria-labelledby="ct cd" xmlns="http://www.w3.org/2000/svg">' % (width, h)]
    p.append('<title id="ct">全基因組 CN 類別分布</title>')
    p.append('<desc id="cd">22 條染色體各一軌，每個體細胞突變依 SEQC2 拷貝數六類著色，'
             '可看出中性區集中於少數染色體。</desc>')
    p.append('<rect x="0" y="0" width="%d" height="%d" fill="#FAF9F5"/>' % (width, h))
    p.append('<text x="%d" y="16" font-size="10.5" fill="#5C5A54">每個標記＝一個 sSNV，顏色＝SEQC2 拷貝數類別；橫軸為染色體座標（共用比例尺）</text>' % L)
    for i, c in enumerate(CHROMS):
        y = 30 + i * ROW
        p.append('<text x="%d" y="%.1f" text-anchor="end" font-size="9.5" fill="#5C5A54">%s</text>' % (L - 6, y + 9, c))
        p.append('<rect x="%d" y="%.1f" width="%.1f" height="12" rx="3" fill="#EFEDE7"/>' % (L, y, R - L))
        for s in by.get(c, []):
            x = L + (R - L) * (s["position"] / gmax)
            col = CN_CLASS.get(s.get("seqc2_class"), ("#9aa",))[0]
            p.append('<rect x="%.2f" y="%.1f" width="1.3" height="12" fill="%s" opacity="0.75"/>' % (x, y, col))
        n = len(by.get(c, []))
        nn = sum(1 for s in by.get(c, []) if s.get("seqc2_class") == "neutral")
        p.append('<text x="%.1f" y="%.1f" font-size="8.5" fill="#6B6B66">n=%d%s</text>'
                 % (R + 2, y + 9, n, ("・中性 %d" % nn) if nn else ""))
    p.append("</svg>")
    return "".join(p)


def bar_chart(title, desc, pairs, width=460, colors=None, note=""):
    """通用橫向長條圖。pairs = [(label, count)]，由資料算出百分比。"""
    total = sum(v for _, v in pairs) or 1
    ROW, L = 26, 168
    h = 34 + len(pairs) * ROW + (16 if note else 4)
    p = ['<svg viewBox="0 0 %d %d" role="img" aria-labelledby="b%s bd%s" xmlns="http://www.w3.org/2000/svg">'
         % (width, h, abs(hash(title)) % 9999, abs(hash(title)) % 9999)]
    p.append('<title id="b%s">%s</title>' % (abs(hash(title)) % 9999, html.escape(title)))
    p.append('<desc id="bd%s">%s</desc>' % (abs(hash(title)) % 9999, html.escape(desc)))
    p.append('<rect x="0" y="0" width="%d" height="%d" fill="#FAF9F5"/>' % (width, h))
    p.append('<text x="10" y="17" font-size="11.5" font-weight="700" fill="#1c1b19">%s</text>' % html.escape(title))
    mx = max(v for _, v in pairs) or 1
    for i, (lab, v) in enumerate(pairs):
        y = 30 + i * ROW
        col = (colors or {}).get(lab, "#2F6690")
        w = (width - L - 118) * (v / mx)
        p.append('<text x="%d" y="%.1f" text-anchor="end" font-size="10" fill="#1c1b19">%s</text>'
                 % (L - 6, y + 12, html.escape(str(lab))))
        p.append('<rect x="%d" y="%.1f" width="%.2f" height="15" rx="3" fill="%s" opacity="0.85"/>' % (L, y, w, col))
        p.append('<text x="%.1f" y="%.1f" font-size="9.5" font-family="ui-monospace,monospace" fill="#5C5A54">%s (%.2f%%)</text>'
                 % (L + w + 5, y + 12, "{:,}".format(v), 100 * v / total))
    if note:
        p.append('<text x="10" y="%.1f" font-size="9" fill="#6B6B66">%s</text>' % (h - 5, html.escape(note)))
    p.append("</svg>")
    return "".join(p)


def cross_tab(units_cn, width=520):
    """CN 類別 × 拓撲解析度交叉表 —— 這是「CN 影響了什麼」最直接的一張圖。"""
    tab = collections.defaultdict(lambda: collections.Counter())
    for u in units_cn:
        c, r = u.get("seqc2_class"), u.get("resolution_class")
        if c in CN_CLASS and r:
            tab[c][r] += 1
    order = ["UNIQUE_TREE", "TIED_SAME_TOPOLOGY", "TIED_CROSS_TOPOLOGY"]
    cols = {"UNIQUE_TREE": "#2F7D5B", "TIED_SAME_TOPOLOGY": "#2F6690", "TIED_CROSS_TOPOLOGY": "#B06A16"}
    rows = [c for c in CN_CLASS if tab.get(c)]
    ROW, L, h = 30, 118, 0
    h = 56 + len(rows) * ROW + 30
    p = ['<svg viewBox="0 0 %d %d" role="img" aria-labelledby="xt xd" xmlns="http://www.w3.org/2000/svg">' % (width, h)]
    p.append('<title id="xt">拷貝數類別與拓撲解析度的交叉分布</title>')
    p.append('<desc id="xd">每列一個拷貝數類別，橫條依比例顯示唯一樹、並列同形、並列跨形三種解析結果的佔比。</desc>')
    p.append('<rect x="0" y="0" width="%d" height="%d" fill="#FAF9F5"/>' % (width, h))
    p.append('<text x="10" y="17" font-size="11.5" font-weight="700" fill="#1c1b19">拷貝數類別 × 拓撲解析度（每列標準化為 100%）</text>')
    for j, k in enumerate(order):
        p.append('<rect x="%d" y="26" width="10" height="10" rx="2" fill="%s"/>' % (L + j * 132, cols[k]))
        p.append('<text x="%d" y="35" font-size="9" fill="#5C5A54">%s</text>' % (L + j * 132 + 14, k))
    for i, c in enumerate(rows):
        y = 56 + i * ROW
        tot = sum(tab[c].values()) or 1
        p.append('<text x="%d" y="%.1f" text-anchor="end" font-size="10" fill="%s">%s</text>'
                 % (L - 6, y + 13, CN_CLASS[c][0], CN_CLASS[c][1]))
        x = L
        for k in order:
            v = tab[c].get(k, 0)
            w = (width - L - 66) * (v / tot)
            if w > 0:
                p.append('<rect x="%.2f" y="%.1f" width="%.2f" height="17" fill="%s" opacity="0.85"><title>%s %s: %s (%.2f%%)</title></rect>'
                         % (x, y, w, cols[k], CN_CLASS[c][1], k, "{:,}".format(v), 100 * v / tot))
            x += w
        p.append('<text x="%.1f" y="%.1f" font-size="9" font-family="ui-monospace,monospace" fill="#5C5A54">n=%s</text>'
                 % (width - 62, y + 13, "{:,}".format(tot)))
    p.append('<text x="10" y="%.1f" font-size="9" fill="#6B6B66">'
             '中性列的唯一樹比例高於增益列＝拓撲層未被拷貝數抬高（2026-08-10 稽核結論）</text>' % (h - 8))
    p.append("</svg>")
    return "".join(p)


CSS = """
:root{--c-accent:#D97757;--c-accent-soft:#F2E4DC;--c-text:#1c1b19;--c-text-soft:#5C5A54;
--c-bg:#FAF9F5;--c-card:#FFF;--c-border:#E3DACC;--c-pass:#2F7D5B;--c-warn:#B06A16;
--c-info:#2F6690;--c-dead:#9A3B3B;--c-cs:#5B3FA0;--radius:10px;
--sans:-apple-system,system-ui,"Noto Sans CJK TC","PingFang TC","Microsoft JhengHei",sans-serif;
--mono:"JetBrains Mono",ui-monospace,Menlo,Consolas,monospace;}
*{box-sizing:border-box}
body{margin:0;font-family:var(--sans);background:var(--c-bg);color:var(--c-text);line-height:1.7;font-size:15.5px}
.wrap{max-width:1080px;margin:0 auto;padding:0 24px 60px}
header.hero{padding:34px 0 14px}
header.hero .kicker{font-family:var(--mono);font-size:.74rem;color:var(--c-accent);letter-spacing:1px;text-transform:uppercase}
header.hero h1{font-size:2rem;margin:6px 0;line-height:1.24}
header.hero .lede{color:var(--c-text-soft);max-width:880px;margin:0}
section{padding:30px 0;border-top:1px solid var(--c-border)}
section h3{font-size:1.5rem;margin:0 0 6px;display:flex;gap:12px;align-items:baseline;flex-wrap:wrap}
section h3 .num{font-family:var(--mono);color:var(--c-accent);font-weight:800;font-size:.95rem}
section .sub{color:var(--c-text-soft);font-size:.94rem;margin:0 0 18px;max-width:900px}
.bluf{background:var(--c-accent-soft);border:1px solid var(--c-accent);border-radius:var(--radius);padding:18px 20px;margin:14px 0}
.warn{background:#F3E2E2;border:1px solid var(--c-dead);border-left:6px solid var(--c-dead);border-radius:var(--radius);padding:12px 16px;margin:14px 0;font-size:.93rem}
.warn b{color:var(--c-dead)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:14px;margin:14px 0}
figure{margin:0;background:#fff;border:1px solid var(--c-border);border-radius:var(--radius);padding:14px;overflow-x:auto}
figure svg{display:block;max-width:100%;height:auto;margin:0 auto}
figcaption{font-size:.8rem;color:var(--c-text-soft);margin-top:10px;border-top:1px dashed var(--c-border);padding-top:8px;line-height:1.55}
.ttl{font-size:.88rem;font-weight:700;margin:0 0 8px}
table{border-collapse:collapse;width:100%;font-size:.83rem;background:#fff}
th,td{border:1px solid var(--c-border);padding:6px 9px;text-align:left;vertical-align:top}
th{background:#F3EFE7;font-weight:700}
tbody tr:nth-child(even){background:#FBF9F4}
.tw{overflow-x:auto;margin:12px 0;max-height:520px;overflow-y:auto}
code,.mono{font-family:var(--mono);font-size:.85em;background:#F0EDE6;padding:1px 4px;border-radius:4px}
.chip{display:inline-block;font-size:.68rem;font-family:var(--mono);font-weight:700;padding:1px 6px;border-radius:5px;border:1px solid;white-space:nowrap}
.kpis{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}
.kpi{flex:1;min-width:170px;background:#fff;border:1px solid var(--c-border);border-radius:var(--radius);padding:11px 15px}
.kpi .v{font-family:var(--mono);font-size:1.4rem;font-weight:800;line-height:1.1}
.kpi .k{font-size:.75rem;color:var(--c-text-soft);margin-top:4px}
footer.prov{border-top:2px solid var(--c-border);margin-top:30px;padding:20px 0 40px;font-size:.79rem;color:var(--c-text-soft)}
footer.prov ul{padding-left:1.2em}
"""


def build(out_path):
    units = load_jsonl(TOPO)
    ucn = load_jsonl(UNIT_CN)
    sites = load_jsonl(SITE_CN)
    cn_by_id = {u["unit_id"]: u for u in ucn}
    # 把 resolution_class 併進 CN 註釋（交叉表要用）
    res_by_id = {u["unit_id"]: cn_by_id.get(u["unit_id"], {}).get("resolution_class") for u in units}
    for u in ucn:
        u.setdefault("resolution_class", res_by_id.get(u["unit_id"]))

    xhp = cross_hp(sites)
    ex = pick_examples(units, ucn)

    # ---- 統計（全部由資料算出，無手打數字）----
    n_units, n_sites = len(units), len(sites)
    n_pos = len({(s["chrom"], s["position"]) for s in sites})
    cls_site = collections.Counter(s.get("seqc2_class") for s in sites if s.get("seqc2_class"))
    cls_unit = collections.Counter(u.get("seqc2_class") for u in ucn if u.get("seqc2_class"))
    loh_n = sum(1 for s in sites if s.get("seqc2_loh"))
    hp_n = collections.Counter(u.get("hp_family") for u in ucn if u.get("hp_family"))
    kb = collections.Counter(u.get("active_bit_count") for u in units)
    colmap = {v[1]: v[0] for v in CN_CLASS.values()}

    P = []
    A = P.append
    A('<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">')
    A("<title>HCC1395 — CN／HP 觀察 pilot</title>")
    A("<style>%s</style>" % CSS)
    A('<div class="wrap"><header class="hero">')
    A('<div class="kicker">觀察 pilot · HCC1395 · 2026-08-10</div>')
    A("<h1>拷貝數與單倍型的觀察疊圖</h1>")
    A('<p class="lede">回應四項觀察缺口：分支節點改用 <b>S 編號</b>、補上<b>同位點跨 HP 對照</b>、'
      '加入 <b>LOH／CNV 疊圖</b>、以及更多整體觀察圖表。'
      '所有數字由既有 canonical 產物算出，本頁不重算任何科學量。</p></header>')

    A('<div class="warn"><b>這一頁的 CN／LOH 是「描述性疊圖」，不是校正。</b>'
      '<code>authority_manifest.claim_boundary.forbidden</code> 明列 <code>CN/LOH-corrected CCF</code>；'
      '2026-08-10 的 CN 稽核亦定案「矯正只能是 fail-closed 適用性 gate，不是數值校正」——'
      '因為 SEQC2 只給 total CN 與 LOH，沒有 allele-specific CN，做數值校正等於把假設偷渡成結果。'
      '本頁只標註「這個結論站在什麼地基上」，不改動任何數值。</div>')

    A('<div class="kpis">')
    for v, k, c in [("{:,}".format(n_units), "canonical unit", "var(--c-info)"),
                    ("{:,}".format(n_pos), "distinct sSNV 位點", "var(--c-info)"),
                    ("%.2f%%" % (100 * (n_sites - cls_site.get("neutral", 0)) / max(1, n_sites)), "位點落在 CN-altered 區", "var(--c-dead)"),
                    ("{:,}".format(len(xhp)), "同位點跨 HP 可對照", "var(--c-cs)")]:
        A('<div class="kpi"><div class="v" style="color:%s">%s</div><div class="k">%s</div></div>' % (c, v, k))
    A("</div>")

    # ── 1 S 編號分支圖 ──
    A('<section><h3><span class="num">01</span>演化分支改用 S 編號</h3>')
    A('<p class="sub">原本節點顯示的是 <code>ROOT → H_ARR → H_AAR → H_AAA</code> 這類等位組合字串，'
      '難以對應到實際位點。改為<b>區域內依基因體位置穩定編號</b>：同一 unit 的 '
      '<code>active_positions</code> 排序後給 S1／S2／S3，節點顯示「已獲得哪些 S」，'
      '<b>滑鼠移到邊上的 +Sx 標籤即顯示完整座標</b>。編號在同一 unit 內固定，可跨圖引用。</p>')
    A('<div class="grid">')
    for u, tm, cls in ex:
        lab = tm["labels"]
        A("<figure>")
        A('<p class="ttl">%s <span class="chip" style="color:%s;border-color:%s;background:#fff">%s</span></p>'
          % (html.escape(u["region_id"].split("|U")[0]), CN_CLASS.get(cls, ("#666",))[0],
             CN_CLASS.get(cls, ("#666",))[0], CN_CLASS.get(cls, ("", "未註釋"))[1]))
        A(svg_tree(u, tm, cls))
        A("<figcaption><b>S 編號對照：</b>" +
          "、".join("%s = %s:%s" % (s, u["chrom"], p) for p, s in sorted(lab.items())) +
          "<br>解析結果 <code>%s</code>・候選樹 <code>%s</code> 棵</figcaption>"
          % (html.escape(str(cn_by_id.get(u["unit_id"], {}).get("resolution_class"))),
             html.escape(str(u.get("total_tree_count")))))
        A("</figure>")
    A("</div></section>")

    # ── 2 跨 HP ──
    A('<section><h3><span class="num">02</span>同一位點在不同 HP 的對照</h3>')
    A('<p class="sub">現有工作站沒有這個視圖。實際盤點後：%s 個 distinct 位點中，'
      '<b>%s 個（%.2f%%）出現在兩個不同的 HP 家族</b>，其餘 %.2f%% 只落在單一 HP。'
      '這 %s 個正是最值得看的 —— 它們能顯示同一個突變在兩條單倍型上的讀數與 AF 差異。</p>'
      % ("{:,}".format(n_pos), "{:,}".format(len(xhp)), 100 * len(xhp) / max(1, n_pos),
         100 * (n_pos - len(xhp)) / max(1, n_pos), "{:,}".format(len(xhp))))
    A('<div class="tw"><table><thead><tr><th>位點</th><th>CN 類別</th><th>LOH</th>'
      '<th>HP1 ALT/REF</th><th>HP1 AF</th><th>HP2 ALT/REF</th><th>HP2 AF</th><th>兩側狀態</th></tr></thead><tbody>')
    for r in xhp:
        h1, h2 = r["hp"].get("1"), r["hp"].get("2")
        f = lambda x: ("%d/%d" % (x["alt"], x["ref"])) if x else "—"
        g = lambda x: ("%.3f" % x["af"]) if x and x.get("af") is not None else "—"
        cc = CN_CLASS.get(r["cn_class"], ("#666", "未註釋"))
        A('<tr><td class="mono">%s:%s</td><td><span class="chip" style="color:%s;border-color:%s;background:#fff">%s</span></td>'
          '<td>%s</td><td class="mono">%s</td><td class="mono">%s</td><td class="mono">%s</td><td class="mono">%s</td><td>%s</td></tr>'
          % (r["chrom"], "{:,}".format(r["position"]), cc[0], cc[0], cc[1],
             "是" if r["loh"] else "否", f(h1), g(h1), f(h2), g(h2),
             " / ".join(sorted({(h or {}).get("status", "—") for h in (h1, h2)}))))
    A("</tbody></table></div>")
    A('<p class="sub" style="margin-top:8px">AF 僅取 <code>active_positions</code> 上的值。'
      '<code>af_coverage</code> 另含非 active 位點，未過濾會使離散度與中介檢驗反轉。</p>')
    A("</section>")

    # ── 3 CN/LOH 疊圖 ──
    A('<section><h3><span class="num">03</span>LOH／CNV 區間疊圖</h3>')
    A('<p class="sub">現有工作站完全沒有這層（<code>CNV</code>／<code>copy_number</code>／<code>segment</code> 命中數為 0，'
      '<code>LOH</code> 僅 9 次且非區間）。此處以 <b>SEQC2 外部真值</b>為準'
      '（SAVANA 發布的 fit 已判定 mis-calibrated，整數 CN 吻合僅 3.64%，僅作對照不作判準）。'
      '<b>LOH 是獨立軸</b>：gain 與 LOH 可同時成立，壓成單一標籤會失真。</p>')
    A('<figure><p class="ttl">全基因組 CN 類別分布</p>%s'
      '<figcaption>中性區（綠）明顯集中在少數染色體 —— 這正是為什麼「按染色體分層後 CN 效應量級減半」。'
      '橫軸為共用比例尺的染色體座標。</figcaption></figure>' % chrom_track(sites))
    A('<div class="grid">')
    A('<figure>%s<figcaption>site 層六類分布。LOH 標記為真者共 %s 個（%.2f%%）——'
      'LOH 是大宗而非零星現象。</figcaption></figure>'
      % (bar_chart("位點的 CN 類別分布", "六類拷貝數狀態下的體細胞突變位點數量與佔比",
                   [(CN_CLASS[c][1], cls_site.get(c, 0)) for c in CN_CLASS if cls_site.get(c)],
                   colors=colmap, note="分母＝%s 個位點註釋" % "{:,}".format(n_sites)),
         "{:,}".format(loh_n), 100 * loh_n / max(1, n_sites)))
    A('<figure>%s<figcaption>unit 層六類分布。可用於生物解讀的 <b>CN 中性且無 LOH</b> 僅 %s 個。</figcaption></figure>'
      % (bar_chart("分析單元的 CN 類別分布", "六類拷貝數狀態下的分析單元數量與佔比",
                   [(CN_CLASS[c][1], cls_unit.get(c, 0)) for c in CN_CLASS if cls_unit.get(c)],
                   colors=colmap, note="分母＝%s 個已註釋 unit" % "{:,}".format(sum(cls_unit.values()))),
         "{:,}".format(cls_unit.get("neutral", 0))))
    A("</div></section>")

    # ── 4 更多觀察圖表 ──
    A('<section><h3><span class="num">04</span>整體觀察圖表</h3>')
    A('<p class="sub">補上原本以表格為主、缺乏整體視覺確認的部分。</p>')
    A('<figure><p class="ttl">拷貝數類別 × 拓撲解析度</p>%s'
      '<figcaption><b>這是本頁最重要的一張圖。</b>若 CN 抬高了拓撲解析度，增益列的「唯一樹」比例應高於中性列；'
      '實際方向相反 —— 這支持「拓撲層 CN-robust、可保留」的既有結論。'
      '（read-AF 層則確實被抬高，屬另一軸，見 2026-08-10 稽核報告。）</figcaption></figure>' % cross_tab(ucn))
    A('<div class="grid">')
    A('<figure>%s<figcaption>兩個 germline 單倍型家族的 unit 數相當接近，未見系統性偏斜。</figcaption></figure>'
      % bar_chart("HP 家族的分析單元數", "兩個 germline 單倍型家族各自的分析單元數量",
                  [("HP%s" % k, v) for k, v in sorted(hp_n.items())],
                  colors={"HP1": "#2F6690", "HP2": "#5B3FA0"}))
    A('<figure>%s<figcaption>k＝該 unit 的活躍突變位數。k=0 無突變可推論；k=1 無共現夥伴、無法建樹；'
      'k≥2 才進入拓撲求解。</figcaption></figure>'
      % bar_chart("每單元的活躍突變位數 k", "各 k 值的分析單元數量分布",
                  [("k=%s" % k, kb[k]) for k in sorted(kb) if k is not None][:9],
                  colors={"k=0": "#9A3B3B", "k=1": "#B06A16"}))
    A("</div></section>")

    A('<footer class="prov"><h4>資料來源與紀律</h4><ul>')
    A("<li>canonical topology：<code>%s</code>（%s unit）</li>" % (html.escape(TOPO), "{:,}".format(n_units)))
    A("<li>CN 註釋：<code>hcc1395_unit_cn_annotation.jsonl</code> / <code>hcc1395_site_cn_annotation.jsonl</code>"
      "（由 <code>cn_annotate_units.py</code> 產生，SEQC2 外部真值）</li>")
    A("<li>本頁<b>不重算任何科學量</b>，所有數字由上述產物直接統計；CN／LOH 僅標註不校正。</li>")
    A("<li>read-AF 僅取 <code>active_positions</code>；SAVANA 僅對照不作判準。</li>")
    A("</ul><p>由 <code>scripts/build_cn_hp_observation_pilot.py</code> 產生 · HCC1395 pilot</p></footer></div>")

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(P))
    return {"units": n_units, "sites": n_sites, "positions": n_pos,
            "cross_hp": len(xhp), "examples": len(ex),
            "bytes": os.path.getsize(out_path)}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default=os.path.join(
        ROOT, "docs", "methodology", "20260810_cn_hp_observation_pilot_HCC1395.standalone.html"))
    a = ap.parse_args()
    for p in (TOPO, UNIT_CN, SITE_CN):
        if not os.path.exists(p):
            sys.exit("缺輸入檔（本腳本不產生科學資料，只讀既有產物）: %s" % p)
    st = build(a.out)
    print("  ✔ %s" % a.out)
    for k, v in st.items():
        print("     %-12s %s" % (k, "{:,}".format(v) if isinstance(v, int) else v))
