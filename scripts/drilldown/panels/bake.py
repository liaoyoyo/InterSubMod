#!/usr/bin/env python3
"""把 ISM 的逐位點產物批次合成 PNG，並產出 HTML 端要用的 manifest。

輸出落在 <out>/panels/<chrom>/<chrom>_<pos>.png，HTML 用相對路徑連結。
預設 --figs-mode copy 的精神：輸出夾自足、可整包搬走。

manifest 記錄每張圖的尺寸與各區塊的像素邊界，讓前端能疊一層 SVG 標軸與圖例
（PNG 本身沒有文字 —— 1 px = 1 cell 的圖放大後文字會糊）。
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

import composite


def load_lineage_map(assign_path, paths_path) -> dict:
    """qname_sha256 → lineage_path。給側欄的 lineage 軌用。

    assignments 給 read→region，paths 給 region→vertex 的階層標籤；
    一條 read 可能落多個 block，取第一個（前端會標明）。
    """
    if not assign_path or not Path(assign_path).is_file():
        return {}
    region_label = {}
    if paths_path and Path(paths_path).is_file():
        with gzip.open(paths_path, "rt", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rid = row.get("region_id")
                lp = (row.get("lineage_path") or "").strip()
                if rid and lp and rid not in region_label:
                    region_label[rid] = lp
    out = {}
    with gzip.open(assign_path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            q = (row.get("qname_sha256") or "").strip()
            rid = (row.get("region_id") or "").strip()
            if q and q not in out:
                lab = region_label.get(rid)
                if lab:
                    out[q] = lab
    return out


def bake(ism_root, out_dir: Path, loci, lineage_map=None, cell_h: int = 2,
         limit: int = 0, log=print) -> dict:
    """loci: [(chrom, pos), ...]。回傳 manifest。"""
    import ism as src_ism           # 借用它的目錄探測，避免兩套路徑推導邏輯

    ism_root = Path(ism_root)
    panels = out_dir / "panels"
    manifest, made, skipped, total_bytes = {}, 0, 0, 0

    for n, (chrom, pos) in enumerate(loci):
        if limit and made >= limit:
            break
        ld = src_ism.locus_dir(ism_root, chrom, pos)
        if ld is None:
            skipped += 1
            continue
        png = panels / chrom / f"{chrom}_{pos}.png"
        try:
            info = composite.build(ld, png, lineage_map=lineage_map, cell_h=cell_h)
        except Exception as exc:                     # noqa: BLE001
            log(f"  [warn] {chrom}:{pos} 合成失敗 {type(exc).__name__}: {exc}")
            skipped += 1
            continue
        if not info:
            skipped += 1
            continue
        info["file"] = f"panels/{chrom}/{chrom}_{pos}.png"
        manifest.setdefault(chrom, {})[str(pos)] = info
        made += 1
        total_bytes += info["bytes"]
        if made % 2000 == 0:
            log(f"  已產 {made:,} 張（{total_bytes / 2**20:.0f} MB）")

    return {"panels": manifest, "made": made, "skipped": skipped,
            "bytes": total_bytes, "cellH": cell_h}


def write_manifest(out_dir: Path, mani: dict, chroms) -> dict:
    """逐染色體分片，與 L2/L4 同一套載入機制。"""
    data = out_dir / "data"
    data.mkdir(parents=True, exist_ok=True)
    out = {}
    for c in chroms:
        rows = mani["panels"].get(c, {})
        p = data / f"L5.{c}.js"
        body = ("window.__DD=window.__DD||{};window.__DD.L5=window.__DD.L5||{};"
                f"window.__DD.L5[{json.dumps(c)}]="
                f"{json.dumps(rows, ensure_ascii=False, separators=(',', ':'))};")
        p.write_text(body, encoding="utf-8")
        out[c] = {"file": f"data/L5.{c}.js", "loci": len(rows), "bytes": p.stat().st_size}
    return out


_ = hashlib
