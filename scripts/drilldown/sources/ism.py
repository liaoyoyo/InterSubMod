#!/usr/bin/env python3
"""ISM 甲基層（擴充能力）。缺了就把甲基面板降級，不影響拓撲骨幹。

一個 ISM run 的形狀：
    <root>/<chrom>/<chrom>/<chrom>/<chrom>_<pos>/<chrom>_<start>_<end>/
        metadata.txt
        methylation/methylation.csv      read × CpG 的 β 矩陣
        distance/<METRIC>/matrix.csv     read × read 距離
        clustering/leaf_order.txt         甲基距離分群的葉序
        clustering/significance.json      含 optimal_k
        reads/reads.tsv                   每條 read 的 hp / alt_support / strand …
    <root>/<chrom>/significance_summary.csv   逐位點的各軸 p 值

⚠ 兩個必須標明的事：
  1. 窗大小決定 CpG 數，直接影響分群強度。不同 run 的窗不同就不可互比。
  2. 「未檢定」（*Valid=false 或群數 <2）是第三態，不等於「不顯著」。
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

from capability import (Capability, Registry, finalize, make, probe,
                        probe_count, probe_exists, probe_linkage)

# 各軸在 significance_summary.csv 的欄名。舊 binary 沒有 Lineage* 欄，
# 探測時會發現並把 lineage 軸標成不可用而非靜默當成「不顯著」。
AXES = [
    {"id": "hp", "title": "HP（merged）", "p": "HPMergedP", "eff": "HPMergedDelta",
     "valid": None, "n": None},
    {"id": "hpfine", "title": "HP（fine）", "p": "HPFineP", "eff": "HPFineF",
     "valid": None, "n": "HPFineNGroups"},
    {"id": "allele", "title": "ALT / REF", "p": "AlleleP", "eff": "AlleleDelta",
     "valid": None, "n": None},
    {"id": "cluster", "title": "甲基自身分群", "p": "ClusterPermanovaP", "eff": "ClusterPermanovaF",
     "valid": "ClusterPermanovaValid", "n": None},
    {"id": "lineage", "title": "lineage", "p": "LineagePermanovaP", "eff": "LineagePermanovaF",
     "valid": "LineagePermanovaValid", "n": "LineageNGroups"},
]

SUMMARY_KEEP = ["Chr", "Pos", "NumReads", "NumCpGs", "ClusterDispersionWarn"]


def _f(v):
    try:
        x = float(v)
        return x if x == x else None      # NaN → None
    except (TypeError, ValueError):
        return None


def _i(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _truthy(v):
    return str(v).strip().lower() in ("1", "true", "yes")


def load(reg: Registry, root, l1, chroms) -> Capability:
    """root = ISM run 根目錄（底下每條染色體一個子目錄）。"""
    cap = make(reg, "ism_dirs", "ISM 甲基層（逐位點產物）",
               enables=["methyl_panel", "filter:axis", "cooccur:ism_k", "layer:ismRing"])
    root = Path(root) if root else None
    if not root or not root.is_dir():
        cap.state = "ABSENT"
        probe(cap, "S0", False, f"ISM run 目錄不存在：{root}")
        return cap
    probe(cap, "S0", True, f"{root}")
    cap.paths.append(_dirref(root))

    # S1：逐染色體讀 significance_summary.csv，順帶偵測有哪些軸
    axis_present, rows_by_pos, n_rows = {}, {}, 0
    window_bp = None
    have_chrom = []
    for c in chroms:
        s = root / c / "significance_summary.csv"
        if not s.is_file():
            continue
        have_chrom.append(c)
        with open(s, newline="") as fh:
            rd = csv.DictReader(fh)
            cols = set(rd.fieldnames or [])
            for a in AXES:
                if a["p"] in cols:
                    axis_present[a["id"]] = True
            for r in rd:
                n_rows += 1
                pos = _i(r.get("Pos"))
                if pos is None:
                    continue
                rec = {k: r.get(k) for k in SUMMARY_KEEP if k in r}
                for a in AXES:
                    if a["id"] not in axis_present:
                        continue
                    rec[a["id"] + "_p"] = _f(r.get(a["p"]))
                    rec[a["id"] + "_eff"] = _f(r.get(a["eff"]))
                    rec[a["id"] + "_n"] = _i(r.get(a["n"])) if a["n"] else None
                    rec[a["id"] + "_valid"] = _truthy(r.get(a["valid"])) if a["valid"] else None
                rows_by_pos[(c, pos)] = rec

    if not have_chrom:
        cap.state = "MALFORMED"
        probe(cap, "S1", False, f"{root} 底下找不到任何 <chrom>/significance_summary.csv")
        return cap

    missing_axes = [a["title"] for a in AXES if a["id"] not in axis_present]
    probe(cap, "S1", True,
          f"{len(have_chrom)} 條染色體有 summary；可用軸 {len(axis_present)}/{len(AXES)}"
          + (f"；缺 {', '.join(missing_axes)}" if missing_axes else ""))
    if missing_axes:
        cap.state = "PARTIAL"
        cap.reason = ("這個 run 的 binary 沒有算 " + "、".join(missing_axes) +
                      " 軸。那些軸顯示為「未檢定」，<b>不等於不顯著</b>。")

    # 窗大小：直接影響 CpG 數與分群強度，不同 run 不可互比
    rp = root / have_chrom[0] / "run_params.json"
    if rp.is_file():
        try:
            window_bp = json.loads(rp.read_text(encoding="utf-8"))["region"]["window_size_bp"]
        except (OSError, ValueError, KeyError):
            window_bp = None
    cap.counts["window_size_bp"] = window_bp
    cap.counts["window_span_bp"] = (window_bp * 2 + 1) if window_bp else None
    probe(cap, "S2", True,
          f"summary 共 {n_rows:,} 列；窗 ±{window_bp:,} bp"
          f"（span {window_bp * 2 + 1:,} bp）" if window_bp else f"summary 共 {n_rows:,} 列")
    probe_count(cap, n_rows, None, "significance_summary 列")

    # S3：與硬核心的連結 —— 有多少 sSNV 對得上 ISM 位點
    chrom_names = l1["chroms"]
    hit, hit_dir = set(), set()
    dir_cache = {}
    for i in range(l1["n"]):
        c = chrom_names[l1["chrom"][i]]
        p = _abs_pos(l1, i)
        if (c, p) in rows_by_pos:
            hit.add(i)
        d = dir_cache.get(c)
        if d is None:
            d = _locus_dir_set(root, c)
            dir_cache[c] = d
        if p in d:
            hit_dir.add(i)

    probe_linkage(cap, len(hit), l1["n"], "sSNV 在 significance_summary 中有對應列")
    cap.counts["locus_dir_hit"] = len(hit_dir)
    probe(cap, "S3", True,
          f"有逐位點產物目錄的 sSNV：{len(hit_dir):,} / {l1['n']:,}"
          f"（{len(hit_dir) / max(l1['n'], 1) * 100:.1f}%）")

    cap.payload = {
        "rows": rows_by_pos,
        "hit": hit,
        "hit_dir": hit_dir,
        "axes": [a for a in AXES if a["id"] in axis_present],
        "missing_axes": missing_axes,
        "window_bp": window_bp,
        "root": str(root),
    }
    return finalize(cap)


def _abs_pos(l1, i):
    """L1 存的是染色體內 delta，這裡還原絕對座標。"""
    if not hasattr(_abs_pos, "_cache") or _abs_pos._cache_id != id(l1):
        acc, prev, out = 0, -1, []
        for j in range(l1["n"]):
            if l1["chrom"][j] != prev:
                acc = l1["dpos"][j]
                prev = l1["chrom"][j]
            else:
                acc += l1["dpos"][j]
            out.append(acc)
        _abs_pos._cache = out
        _abs_pos._cache_id = id(l1)
    return _abs_pos._cache[i]


def _locus_dir_set(root: Path, chrom: str) -> set:
    """ISM 的巢狀是 <root>/<chrom>/<chrom>/<chrom>/<chrom>_<pos>/。
    層數會因 CLI 的 -o 而異，所以往下探到出現 <chrom>_<pos> 為止。"""
    cur = root / chrom
    for _ in range(4):
        if not cur.is_dir():
            return set()
        names = os.listdir(cur)
        pos = set()
        for n in names:
            if n.startswith(chrom + "_"):
                tail = n[len(chrom) + 1:]
                if tail.isdigit():
                    pos.add(int(tail))
        if pos:
            return pos
        nxt = cur / chrom
        if nxt.is_dir():
            cur = nxt
        else:
            return set()
    return set()


def _dirref(p: Path):
    from capability import FileRef
    return FileRef(path=str(p), exists=True, size_bytes=0, sha256="", mtime="")


def locus_dir(root: Path, chrom: str, pos: int):
    """回傳某位點的產物目錄（含視窗子目錄），找不到回 None。"""
    cur = Path(root) / chrom
    for _ in range(4):
        cand = cur / f"{chrom}_{pos}"
        if cand.is_dir():
            subs = [d for d in cand.iterdir() if d.is_dir()]
            return subs[0] if subs else cand
        nxt = cur / chrom
        if not nxt.is_dir():
            return None
        cur = nxt
    return None
