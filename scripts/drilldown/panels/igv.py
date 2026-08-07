#!/usr/bin/env python3
"""IGV 式 read 對齊圖：每條 read 一條水平長條，依 lineage 分組，
標出各 sSNV 的 REF/ALT 與 CpG 甲基。

概念取自手刻的 docs/methodology/20260625_chr17_igv_subclone_explainer_01.standalone.html
（1080×528 SVG、2,655 個 rect）。這裡把它自動化。

## 四條資料鏈（workflow 已逐 region 實測）

    read 座標    ISM <region>/reads/reads.tsv 的 start / end（全基因組座標）
    各 sSNV 等位  同檔 alt_support（該窗的位點）；跨窗合併用 read_name
    lineage 類別  同檔 lineage_path（來自 lineage_tagged BAM 的 lv 標籤）
    CpG 甲基     <region>/methylation/methylation.csv（欄名即 CpG 座標）

## 三個必須守住的界線

1. **lineage 類別在 k≥7 完全不可得** —— 那 153 個 region 全是
   family_incomplete / ABSTAIN_RESOURCE_LIMIT，solver 根本沒建樹，
   不是標籤漏掉。那種 region 退回用 pattern 基因型分組，並明確標示。
2. **ISM 只在每個 sSNV ±1kb 開窗**，寬 region 的甲基是斷續的
   （實測 chr4 span 10,140 只有 53.0% 有 CpG 資料）。無資料段留白，
   不用色階填滿假裝有測。
3. **read 沒有 lineage 標籤 ≠ 不屬於任何 lineage** —— 它可能只是沒跨到
   足夠位點。歸入「未標記」組並標明數量。
"""
from __future__ import annotations

import csv
import os
from collections import OrderedDict, defaultdict
from pathlib import Path

from pngenc import hex_to_rgb          # noqa: F401  （色碼轉換共用）

import sys
_STD = Path(__file__).resolve().parent.parent.parent / "ism_heatmap_std.py"
sys.path.insert(0, str(_STD.parent))
import ism_heatmap_std as H            # noqa: E402

# 甲基 β → 色（沿用視覺規約 SoT，不另立色盤）
def beta_color(b):
    return H.rdbu_hex(b, 12)


ALT_COL = {"ALT": "#dc2626", "REF": "#fbbf24"}
READ_BG = {"tumor": "#3b4a5a", "normal": "#4a5a3b"}
NA_GREY = "#9ca3af"


def _read_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def _meth_of(locus_dir: Path):
    """回傳 {read_id: {cpg_pos: beta}} 與 CpG 座標清單。"""
    p = locus_dir / "methylation" / "methylation.csv"
    if not p.is_file():
        return {}, []
    with open(p, newline="") as fh:
        rd = csv.reader(fh)
        head = next(rd)
        cpgs = []
        for c in head[1:]:
            try:
                cpgs.append(int(c))
            except ValueError:
                cpgs.append(None)
        out = {}
        for r in rd:
            if not r:
                continue
            vals = {}
            for i, x in enumerate(r[1:]):
                x = x.strip()
                if x in ("", "NA", "nan") or i >= len(cpgs) or cpgs[i] is None:
                    continue
                try:
                    vals[cpgs[i]] = float(x)
                except ValueError:
                    pass
            out[r[0]] = vals
    return out, [c for c in cpgs if c is not None]


def collect(ism_root, chrom, positions, locus_dir_fn):
    """跨該 region 的所有 sSNV 窗，把 read 資訊合併起來。

    合併鍵是 read_name —— read_id 是**窗內索引**，跨窗不通用。
    """
    reads = OrderedDict()          # read_name -> {...}
    cpg_all = {}                   # cpg_pos -> {read_name: beta}
    covered = []                   # 有 ISM 窗的 sSNV
    windows = []

    for pos in positions:
        ld = locus_dir_fn(ism_root, chrom, pos)
        if ld is None:
            continue
        covered.append(pos)
        rows = _read_rows(ld / "reads" / "reads.tsv")
        if not rows:
            continue
        # 窗的實際範圍（供圖上標示哪些區段有甲基資料）
        try:
            lo, hi = ld.name.split("_")[1:3]
            windows.append((int(lo), int(hi)))
        except (ValueError, IndexError):
            pass

        meth, _cpgs = _meth_of(ld)
        for r in rows:
            nm = (r.get("read_name") or "").strip()
            if not nm:
                continue
            e = reads.get(nm)
            if e is None:
                e = reads[nm] = {
                    "start": None, "end": None, "hp": (r.get("hp") or "").strip(),
                    "strand": (r.get("strand") or "").strip(),
                    "tn": H.tn_of(r.get("is_tumor") or ""),
                    "lin": (r.get("lineage_path") or "").strip(),
                    "lst": (r.get("lineage_status") or "").strip(),
                    "alleles": {},
                }
            try:
                s, t = int(r["start"]), int(r["end"])
                e["start"] = s if e["start"] is None else min(e["start"], s)
                e["end"] = t if e["end"] is None else max(e["end"], t)
            except (KeyError, TypeError, ValueError):
                pass
            a = (r.get("alt_support") or "").strip()
            if a in ("ALT", "REF"):
                e["alleles"][pos] = a
            if not e["lin"]:
                lp = (r.get("lineage_path") or "").strip()
                if lp and lp != ".":
                    e["lin"] = lp
                    e["lst"] = (r.get("lineage_status") or "").strip()

            rid = (r.get("read_id") or "").strip()
            for cp, b in (meth.get(rid) or {}).items():
                cpg_all.setdefault(cp, {})[nm] = b

    return reads, cpg_all, covered, windows


def group_reads(reads, positions):
    """依 lineage 標籤分組；沒有標籤的退回用等位 pattern 分組。

    回傳 [(組名, 說明, [read_name…])]，並標明該組是「lineage 標籤」
    還是「等位 pattern」—— 兩者證據等級完全不同。
    """
    by_lin = defaultdict(list)
    unlabeled = []
    for nm, e in reads.items():
        lp = e["lin"]
        if lp and lp not in (".", ""):
            # 多 block 時是逗號分隔，取第一個
            by_lin[lp.split(",")[0]].append(nm)
        else:
            unlabeled.append(nm)

    groups = []
    for lp in sorted(by_lin, key=lambda x: (x.count("-"), x)):
        note = "lineage 標籤"
        if lp.endswith("+"):
            note = "lineage 標籤（子樹斷言：此節點或其後代）"
        groups.append((lp, note, by_lin[lp], "lineage"))

    if unlabeled:
        # 沒有 lineage 標籤者依等位 pattern 再分，讓圖上仍看得出結構
        by_pat = defaultdict(list)
        for nm in unlabeled:
            al = reads[nm]["alleles"]
            pat = "".join("A" if al.get(p) == "ALT" else
                          "R" if al.get(p) == "REF" else "X" for p in positions)
            by_pat[pat].append(nm)
        for pat in sorted(by_pat, key=lambda p: (-len(by_pat[p]), p)):
            groups.append((pat, "等位 pattern（無 lineage 標籤）", by_pat[pat], "pattern"))
    return groups


def render_svg(reads, cpg_all, positions, covered, windows, groups,
               chrom, focus_pos=None, width=1000, row_h=6, gap=3, zoom=True):
    """手刻 inline SVG。read 數通常數百條以內，SVG 體積可接受
    （不像甲基矩陣那樣 read×CpG 會爆到 1.5 MB）。

    zoom=True 時 x 軸只涵蓋 sSNV 跨度 ±20% —— ONT read 常長達 40 kb，
    而 sSNV 跨度可能只有幾百 bp，用全 read 跨度畫會讓所有位點擠成一條線
    （實測 525 bp / 90 kb）。超出視窗的 read 段畫成截斷箭頭，不假裝沒有。
    """
    if not any(g[2] for g in groups):
        return "<div class='capability-off'>此 region 沒有可畫的 read。</div>"

    rlo = min((e["start"] for e in reads.values() if e["start"] is not None),
              default=positions[0])
    rhi = max((e["end"] for e in reads.values() if e["end"] is not None),
              default=positions[-1])
    if zoom:
        pad = max((positions[-1] - positions[0]) * 0.25, 300)
        lo, hi = positions[0] - pad, positions[-1] + pad
    else:
        lo, hi = rlo, rhi
    span = max(hi - lo, 1)
    padL, padR, padT = 112, 16, 50
    plot = width - padL - padR

    def X(p):
        return padL + plot * ((min(max(p, lo), hi) - lo) / span)

    LBL = 15                      # 每組標籤佔的高度
    height = padT + sum(LBL + len(g[2]) * row_h + gap for g in groups) + 46
    S = [f"<svg viewBox='0 0 {width} {height}' width='100%' "
         f"style='background:#11161f' font-family='ui-monospace,monospace' font-size='9'>"]

    for (a, b) in windows:
        if b < lo or a > hi:
            continue
        S.append(f"<rect x='{X(a):.1f}' y='{padT - 8}' "
                 f"width='{max(X(b) - X(a), 1):.1f}' "
                 f"height='{height - padT - 34}' fill='#1a2230'/>")

    S.append(f"<line x1='{padL}' y1='{padT - 14}' x2='{width - padR}' y2='{padT - 14}' "
             f"stroke='#4a5568'/>")
    for t in range(5):
        p = lo + span * t / 4
        S.append(f"<line x1='{X(p):.1f}' y1='{padT - 18}' x2='{X(p):.1f}' y2='{padT - 14}' "
                 f"stroke='#4a5568'/>")
        S.append(f"<text x='{X(p):.1f}' y='{padT - 22}' text-anchor='middle' "
                 f"fill='#98a09a'>{p / 1e6:.4f}M</text>")

    # sSNV 垂直線 + S 標籤（錯開避免重疊）
    for i, p in enumerate(positions):
        col = "#facc15" if p == focus_pos else "#8b95a5"
        S.append(f"<line x1='{X(p):.1f}' y1='{padT - 14}' x2='{X(p):.1f}' y2='{height - 32}' "
                 f"stroke='{col}' stroke-width='{1.8 if p == focus_pos else .9}' "
                 f"stroke-dasharray='2 2' opacity='.7'/>")
        yy = height - 20 + (i % 2) * 11
        S.append(f"<text x='{X(p):.1f}' y='{yy}' text-anchor='middle' fill='{col}' "
                 f"font-weight='{700 if p == focus_pos else 400}'>S{i + 1}</text>")

    y = padT
    for gname, gnote, members, gkind in groups:
        gcol = "#db2777" if gkind == "lineage" else "#6b7280"
        S.append(f"<line x1='0' y1='{y - 2}' x2='{width}' y2='{y - 2}' "
                 f"stroke='#2a3441' stroke-width='.6'/>")
        S.append(f"<text x='4' y='{y + 9}' fill='{gcol}' font-weight='700'>"
                 f"{gname[:15]}</text>")
        S.append(f"<text x='4' y='{y + 19}' fill='#667069' font-size='7.5'>"
                 f"{len(members)} 條 · {'lineage' if gkind == 'lineage' else 'pattern'}</text>")
        y += LBL
        for nm in members:
            e = reads[nm]
            if e["start"] is None:
                y += row_h
                continue
            x0, x1 = X(e["start"]), X(e["end"])
            bg = READ_BG["tumor"] if e["tn"] == "1" else READ_BG["normal"]
            S.append(f"<rect x='{x0:.1f}' y='{y}' width='{max(x1 - x0, 1):.1f}' "
                     f"height='{row_h - 1}' fill='{bg}' rx='1'><title>{nm[:8]}… "
                     f"{'tumor' if e['tn'] == '1' else 'normal'} · HP{e['hp']} · "
                     f"{e['end'] - e['start']:,} bp · {chrom}:{e['start']:,}-{e['end']:,}"
                     f"</title></rect>")
            # 被視窗截斷的兩端畫箭頭，表示 read 還有延續
            if e["start"] < lo:
                S.append(f"<polygon points='{padL},{y + row_h / 2} {padL + 5},{y} "
                         f"{padL + 5},{y + row_h - 1}' fill='{bg}'/>")
            if e["end"] > hi:
                xe = width - padR
                S.append(f"<polygon points='{xe},{y + row_h / 2} {xe - 5},{y} "
                         f"{xe - 5},{y + row_h - 1}' fill='{bg}'/>")
            for cp, per in cpg_all.items():
                b = per.get(nm)
                if b is None or not (lo <= cp <= hi) or not (e["start"] <= cp <= e["end"]):
                    continue
                S.append(f"<rect x='{X(cp) - 1:.1f}' y='{y}' width='2' "
                         f"height='{row_h - 1}' fill='{beta_color(b)}'/>")
            for p, a in e["alleles"].items():
                if not (lo <= p <= hi):
                    continue
                S.append(f"<rect x='{X(p) - 2:.1f}' y='{y - .5}' width='4' "
                         f"height='{row_h}' fill='{ALT_COL[a]}'/>")
            y += row_h
        y += gap

    S.append("</svg>")
    return "".join(S)


def legend_html(reads, cpg_all, positions, covered, groups, zoom=True):
    """圖例與判讀界線。這張圖很容易被過度解讀，界線必須在圖旁邊而非別處。"""
    n_lin = sum(len(g[2]) for g in groups if g[3] == "lineage")
    n_pat = sum(len(g[2]) for g in groups if g[3] == "pattern")
    n_t = sum(1 for e in reads.values() if e["tn"] == "1")
    miss = [p for p in positions if p not in covered]
    return (
        "<div class='legend' style='margin:.4rem 0'>"
        "<span class='bar-key'><i class='swatch' style='background:#3b4a5a'></i>tumor read</span>"
        "<span class='bar-key'><i class='swatch' style='background:#4a5a3b'></i>normal read</span>"
        "<span class='bar-key'><i class='swatch' style='background:#dc2626'></i>該位點 ALT</span>"
        "<span class='bar-key'><i class='swatch' style='background:#fbbf24'></i>該位點 REF</span>"
        "<span class='bar-key'><i class='swatch' style='background:#2166ac'></i>CpG β→0</span>"
        "<span class='bar-key'><i class='swatch' style='background:#b2182b'></i>CpG β→1</span>"
        "<span class='bar-key'><i class='swatch' style='background:#1a2230'></i>有 ISM 甲基窗</span>"
        "</div>"
        "<div class='denom' style='line-height:1.75'>"
        f"共 <b>{len(reads):,}</b> 條 read（tumor {n_t:,} / normal {len(reads) - n_t:,}）、"
        f"<b>{len(cpg_all):,}</b> 個 CpG、<b>{len(positions)}</b> 個 sSNV。<br>"
        f"分組：<b>{n_lin:,}</b> 條有 lineage 標籤（粉紅組名），"
        f"<b>{n_pat:,}</b> 條沒有、退回用等位 pattern 分組（灰色組名）。"
        "<b>沒有 lineage 標籤不等於不屬於任何 lineage</b> —— "
        "多半只是那條 read 沒跨到足夠位點。<br>"
        + (f"⚠ <b>{len(miss)} 個 sSNV 沒有 ISM 窗</b>（"
           + "、".join(f"{p:,}" for p in miss[:4]) + "…）：那幾條垂直線上不會有甲基色塊，"
           "是沒測不是沒甲基。<br>" if miss else "")
        + ("深色底的橫向區段是有 ISM 甲基資料的窗（每個 sSNV ±1 kb）；"
           "窗外留白代表<b>沒測</b>，不是沒甲基。<br>")
        + ("x 軸已縮放到 sSNV 跨度 ±25%；超出視窗的 read 段畫成箭頭。<br>" if zoom else "")
        + "<b>組名帶 <code>+</code> 者是子樹斷言</b>（此節點或其後代），"
        "與不帶 <code>+</code> 的<b>不可相加</b>。"
        "</div>")


def build_one(args):
    """給 ProcessPoolExecutor 用的單 region 產圖。回傳 (region_id, dict|None)。"""
    ism_root, chrom, region_id, positions = args
    import sys as _s
    from pathlib import Path as _P
    _here = _P(__file__).resolve().parent
    for _p in (str(_here), str(_here.parent), str(_here.parent / "sources")):
        if _p not in _s.path:
            _s.path.insert(0, _p)
    import ism as _ism
    try:
        reads, cpg, cov, wins = collect(ism_root, chrom, positions, _ism.locus_dir)
    except Exception:                                  # noqa: BLE001
        return region_id, None
    if len(reads) < 4 or not cov or len(reads) > 600:
        return region_id, None
    g = group_reads(reads, positions)
    return region_id, {
        "svg": render_svg(reads, cpg, positions, cov, wins, g, chrom),
        "legend": legend_html(reads, cpg, positions, cov, g),
        "n": len(reads), "cpg": len(cpg), "cov": len(cov),
    }


def build_many(ism_root, regions, workers=0, log=print):
    """對合格 region 平行產圖。回傳 {region_id: dict}。

    合格條件：2–8 個 sSNV、跨度 ≤50 kb —— 更大的 region read 數會讓 SVG 過胖，
    且 k≥7 的 region solver 已 abstain（沒有 lineage 分組可畫）。
    """
    import concurrent.futures as cf
    import os as _os

    tasks = []
    for r in regions:
        ap = r.get("ap") or []
        if 2 <= len(ap) <= 8 and (ap[-1] - ap[0]) <= 50000:
            tasks.append((str(ism_root), r["c"], r["id"], list(ap)))
    if not tasks:
        return {}
    if not workers:
        workers = max(1, min(24, (_os.cpu_count() or 4) - 2))
    log(f"  IGV 圖：{len(tasks):,} 個合格 region（2–8 sSNV、跨度 ≤50 kb），平行度 {workers}")

    out, done = {}, 0
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        for rid, v in ex.map(build_one, tasks, chunksize=16):
            done += 1
            if v:
                out[rid] = v
            if done % 2000 == 0:
                log(f"    已處理 {done:,} / {len(tasks):,}（產出 {len(out):,}）")
    log(f"  IGV 圖產出 {len(out):,} / {len(tasks):,} 個 region")
    return out
