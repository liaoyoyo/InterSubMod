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
import os
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


def _one(args):
    """單一位點的合成。給 ProcessPoolExecutor 用，所以必須是模組層函式。"""
    ism_root, out_dir, chrom, pos, lineage_map, cell_h = args
    import ism as src_ism
    ld = src_ism.locus_dir(Path(ism_root), chrom, pos)
    if ld is None:
        return chrom, pos, None, "no-dir"
    png = Path(out_dir) / "panels" / chrom / f"{chrom}_{pos}.png"
    tpng = Path(out_dir) / "panels" / chrom / f"{chrom}_{pos}.T.png"
    try:
        info = composite.build(ld, png, lineage_map=lineage_map, cell_h=cell_h)
    except Exception as exc:                          # noqa: BLE001
        return chrom, pos, None, f"{type(exc).__name__}: {exc}"
    if not info:
        return chrom, pos, None, "insufficient"
    info["file"] = f"panels/{chrom}/{chrom}_{pos}.png"
    # tumor-only 版本。失敗（例如 T read 不足 3 條）不影響主版本，只是沒有切換。
    try:
        tinfo = composite.build(ld, tpng, lineage_map=lineage_map, cell_h=cell_h,
                                tumor_only=True)
        if tinfo:
            tinfo["file"] = f"panels/{chrom}/{chrom}_{pos}.T.png"
            info["tumorOnly"] = tinfo
    except Exception:                                  # noqa: BLE001
        pass
    return chrom, pos, info, None


def bake(ism_root, out_dir: Path, loci, lineage_map=None, cell_h: int = 2,
         limit: int = 0, workers: int = 0, log=print) -> dict:
    """loci: [(chrom, pos), ...]。回傳 manifest。

    逐位點獨立，所以用 process pool 平行 —— 單執行緒實測約 0.9 s/位點，
    16,304 個要 4 小時；平行後降到十幾分鐘。
    """
    import concurrent.futures as cf

    if limit:
        loci = loci[:limit]
    if not workers:
        workers = max(1, min(24, (os.cpu_count() or 4) - 2))

    # 各 chrom 目錄先建好，避免子行程同時 mkdir 競爭
    for c in {c for c, _ in loci}:
        (Path(out_dir) / "panels" / c).mkdir(parents=True, exist_ok=True)

    manifest, made, skipped, total_bytes = {}, 0, 0, 0
    reasons = {}
    tasks = [(str(ism_root), str(out_dir), c, p, lineage_map, cell_h) for c, p in loci]

    log(f"  平行度 {workers}，共 {len(tasks):,} 個位點")
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        for chrom, pos, info, why in ex.map(_one, tasks, chunksize=32):
            if info is None:
                skipped += 1
                reasons[why] = reasons.get(why, 0) + 1
                continue
            manifest.setdefault(chrom, {})[str(pos)] = info
            made += 1
            total_bytes += info["bytes"]
            if made % 2000 == 0:
                log(f"  已產 {made:,} 張（{total_bytes / 2**20:.0f} MB）")

    if reasons:
        log("  略過原因：" + "、".join(f"{k}×{v:,}" for k, v in sorted(reasons.items())))
    return {"panels": manifest, "made": made, "skipped": skipped,
            "bytes": total_bytes, "cellH": cell_h, "skipReasons": reasons}


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
