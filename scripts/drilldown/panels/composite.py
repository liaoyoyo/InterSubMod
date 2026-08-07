#!/usr/bin/env python3
"""從 ISM 的逐位點產物合成雙面板圖：左 read×CpG 甲基、右 read×read 距離。

read 的排序用 clustering/leaf_order.txt —— 那是甲基距離分群的葉序。
這很重要：依 HP 排序（既有預渲染 PNG 的做法）看不出「甲基自己怎麼分」，
依葉序排才能讓分群結構浮現，再對照左側色帶看它跟哪個標籤軸對齊。

側欄六軌，順序沿用 ism_heatmap_std 的慣例（cluster 最靠熱圖）：
    HP │ ALT │ T/N │ Strand │ lineage │ cluster
lineage 軌需要 sha256(read_name) → qname_sha256 的 join；join 不到的塗灰
並在圖說標明 n/N，不假裝有資料。
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

from pngenc import Canvas, hex_to_rgb

# 視覺規約 SoT。找不到就 fail，不自己另立色盤 —— 那會造成第四套配色。
_ISM_STD = Path(__file__).resolve().parent.parent.parent / "ism_heatmap_std.py"
sys.path.insert(0, str(_ISM_STD.parent))
import ism_heatmap_std as H            # noqa: E402

NA = hex_to_rgb("#d1d5db")
BG = "#11161f"
GAP = 4


def _rdbu(v):
    return hex_to_rgb(H.rdbu_hex(v, 12))


def _magma(t):
    return hex_to_rgb(H.magma_hex(t, 16))


def read_matrix(path):
    """回傳 (read_ids, cpg_positions, rows) —— rows[i][j] 是 β 或 None。"""
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        cpgs = [c for c in header[1:]]
        ids, rows = [], []
        for r in rd:
            if not r:
                continue
            ids.append(r[0])
            vals = []
            for x in r[1:]:
                x = x.strip()
                if x in ("", "NA", "nan", "NaN"):
                    vals.append(None)
                else:
                    try:
                        vals.append(float(x))
                    except ValueError:
                        vals.append(None)
            rows.append(vals)
    return ids, cpgs, rows


def read_square(path):
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        ids = [c for c in header[1:]]
        m = {}
        for r in rd:
            if not r:
                continue
            m[r[0]] = [None if x.strip() in ("", "NA", "nan") else float(x) for x in r[1:]]
    return ids, m


def read_reads_tsv(path):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rid = (r.get("read_id") or "").strip()
            if not rid:
                continue
            out[rid] = {
                "hp": (r.get("hp") or "").strip(),
                "alt_support": (r.get("alt_support") or "").strip(),
                "is_tumor": (r.get("is_tumor") or "").strip(),
                "strand": (r.get("strand") or "").strip(),
                "name": (r.get("read_name") or "").strip(),
            }
    return out


def leaf_order(path):
    try:
        return [ln.strip() for ln in open(path) if ln.strip()]
    except OSError:
        return []


def flat_clusters(linkage_path, leaves, k: int):
    """從 linkage matrix 做平面切割。◆ 儀表板自算 —— 不是 C++ 產物。

    ⚠ 本專案的 linkage_matrix.csv 有**顯式** new_cluster_id 欄
    （cluster_i, cluster_j, distance, new_cluster_id, size），
    不是 scipy 的隱式 n+i 編號。曾按 scipy 慣例解讀而算出 k=41（實際 optimal_k=3）。

    leaves = 依 leaf_order 的 read 清單；回傳 {read: cluster_index}。
    """
    if k < 2 or len(leaves) < 2:
        return None
    merges = []
    try:
        with open(linkage_path, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                merges.append((int(float(r["cluster_i"])), int(float(r["cluster_j"])),
                               int(float(r["new_cluster_id"]))))
    except (OSError, ValueError, KeyError, TypeError):
        return None
    if not merges:
        return None

    # 葉節點 0..n-1 對應 leaf_order 的順序
    n = len(leaves)
    parent = {}

    def find(a):
        while parent.get(a, a) != a:
            parent[a] = parent.get(parent[a], parent[a])
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
        return ra

    take = max(0, len(merges) - (k - 1))
    for i, j, nid in merges[:take]:
        parent[nid] = union(i, j)

    roots, out = {}, {}
    for idx, rid in enumerate(leaves):
        r = find(idx)
        if r not in roots:
            roots[r] = len(roots)
        out[rid] = roots[r]
    return out


def sidebars(ids, reads, lineage_of, clusters):
    """回傳 [(label, [rgb per read])]。順序：HP ALT T/N Strand lineage cluster"""
    def col(fn):
        return [hex_to_rgb(fn(i)) for i in ids]

    bars = [
        ("HP", col(lambda i: H.HP_COL.get(reads.get(i, {}).get("hp", ""), H.NA_GREY))),
        ("ALT", col(lambda i: H.ALT_COL.get(reads.get(i, {}).get("alt_support", ""), H.NA_GREY))),
        ("T/N", col(lambda i: H.TN_COL.get(H.tn_of(reads.get(i, {}).get("is_tumor", "")), H.NA_GREY))),
        ("Strand", col(lambda i: H.STRAND_COL.get(reads.get(i, {}).get("strand", ""), H.NA_GREY))),
    ]
    if lineage_of:
        bars.append(("lineage", col(lambda i: H.cluster_tree_color(lineage_of.get(i, "")))))
    if clusters:
        cmap = {}
        for c in sorted(set(clusters)):
            if c < 0:
                continue
            cmap[c] = H.CLUSTER_COL[len(cmap) % len(H.CLUSTER_COL)]
        bars.append(("cluster",
                     [NA if c < 0 else hex_to_rgb(cmap[c]) for c in clusters]))
    return bars


def build(locus_dir: Path, out_png: Path, lineage_map=None, cell_h: int = 1):
    """合成一張雙面板 PNG。回傳 dict（含尺寸、read/CpG 數、lineage join 率），
    資料不足回 None。"""
    meth = locus_dir / "methylation" / "methylation.csv"
    if not meth.is_file():
        return None
    ids, cpgs, rows = read_matrix(meth)
    if len(ids) < 3 or not cpgs:
        return None

    reads = read_reads_tsv(locus_dir / "reads" / "reads.tsv")

    # 距離矩陣（metric 目錄名不固定，取第一個）
    dist_ids, dist = [], {}
    dd = locus_dir / "distance"
    if dd.is_dir():
        for metric in sorted(p for p in dd.iterdir() if p.is_dir()):
            mp = metric / "matrix.csv"
            if mp.is_file():
                dist_ids, dist = read_square(mp)
                break

    # 排序用甲基分群的葉序 —— 依 HP 排看不出甲基自己怎麼分。
    # ⚠ leaf_order.txt 存的是 read NAME（UUID），methylation.csv 用的是 read_id 整數，
    #   兩者交集為 0，必須經 reads.tsv 的 read_id↔read_name 轉換（曾直接比對而全數落空）。
    name_to_id = {}
    for rid, info in reads.items():
        if info.get("name"):
            name_to_id[info["name"]] = rid
    idset = set(ids)
    lo_names = leaf_order(locus_dir / "clustering" / "leaf_order.txt")
    order = [name_to_id[nm] for nm in lo_names
             if nm in name_to_id and name_to_id[nm] in idset]
    clustered = set(order)
    # 沒進分群的 read 排在最後 —— 它們不是被丟掉，是分群階段沒收進去
    tail = [i for i in ids if i not in clustered]
    n_clustered = len(order)
    order = order + tail
    idx_of = {rid: i for i, rid in enumerate(ids)}

    k = None
    sig = locus_dir / "clustering" / "significance.json"
    if sig.is_file():
        try:
            k = int(json.loads(sig.read_text()).get("optimal_k") or 0)
        except (OSError, ValueError, TypeError):
            k = None
    cmap_by_read = flat_clusters(locus_dir / "clustering" / "linkage_matrix.csv",
                                 order[:n_clustered], k or 0)
    clusters = None
    if cmap_by_read:
        # 未進分群者給 -1，畫成灰 —— 不併進任何一群
        clusters = [cmap_by_read.get(rid, -1) for rid in order]

    lin_of, lin_hit = {}, 0
    if lineage_map:
        for rid in order:
            nm = reads.get(rid, {}).get("name", "")
            if nm:
                lab = lineage_map.get(hashlib.sha256(nm.encode()).hexdigest())
                if lab:
                    lin_of[rid] = lab
                    lin_hit += 1

    bars = sidebars(order, reads, lin_of if lin_of else None, clusters)
    nbar = len(bars)
    SB, nrow = 3, len(order)
    ncol, ndist = len(cpgs), len(dist_ids)
    ch = max(1, cell_h)

    w = nbar * SB + ncol + (GAP + nbar * SB + ndist if ndist else 0)
    h = nrow * ch
    cv = Canvas(w, h, BG)

    for r, rid in enumerate(order):
        y = r * ch
        for b, (_lab, colors) in enumerate(bars):
            cv.rect(b * SB, y, SB - 1, ch, colors[r])
        src = rows[idx_of[rid]]
        x0 = nbar * SB
        for c in range(ncol):
            v = src[c] if c < len(src) else None
            cv.rect(x0 + c, y, 1, ch, NA if v is None else _rdbu(v))

        if ndist:
            x1 = nbar * SB + ncol + GAP
            for b, (_lab, colors) in enumerate(bars):
                cv.rect(x1 + b * SB, y, SB - 1, ch, colors[r])
            drow = dist.get(rid)
            x2 = x1 + nbar * SB
            if drow:
                dpos = {d: i for i, d in enumerate(dist_ids)}
                mx = max((v for v in drow if v is not None), default=1.0) or 1.0
                for c, cid in enumerate(order):
                    j = dpos.get(cid)
                    v = drow[j] if (j is not None and j < len(drow)) else None
                    cv.rect(x2 + c, y, 1, ch, NA if v is None else _magma(min(v / mx, 1.0)))

    out_png.parent.mkdir(parents=True, exist_ok=True)
    nbytes = cv.save(out_png)
    return {
        "bytes": nbytes, "w": w, "h": h,
        "reads": nrow, "cpgs": ncol, "hasDist": bool(ndist),
        "bars": [b[0] for b in bars],
        "barW": SB, "gap": GAP,
        "lineageHit": lin_hit, "lineageTotal": nrow,
        "clusterK": (len({c for c in clusters if c >= 0}) if clusters else None),
        "optimalK": k,
        "clustered": n_clustered, "notClustered": nrow - n_clustered,
    }
