#!/usr/bin/env python3
"""把探測到的資料組成前端要的形狀：shell 內嵌的 boot payload + 逐染色體分片。

分片用 `<script src>` 而非 fetch —— Chrome 對 file:// 頁面的 fetch 一律以 origin
null 阻擋，但 subresource script 可以正常載入。所以每個分片是一行賦值：
    window.__DD.L2.chr7 = {...};
"""
from __future__ import annotations

import json
from pathlib import Path

# 篩選維度的定義。每個維度的 keys 決定 UI 上有哪些勾選項，
# 初始一律全勾（預設全顯示），使用者取消勾選才會收斂。
K_COLORS = ["#c9c8bd", "#176b58", "#2f8f76", "#5aa892", "#8bbfae", "#b66e20",
            "#c98a45", "#a94336", "#6b5592"]


def json_for_html(obj) -> str:
    """`</` 會提前關閉 script 標籤，必須跳脫。"""
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def build_dims(topo_cap, extra=None) -> list:
    """從硬核心可得的欄位組出基本維度。擴充層（ISM / lineage / BED）在 extra 加。"""
    regions = topo_cap.payload["regions"]
    l1 = topo_cap.payload["l1"]

    ks = sorted({min(k, 8) for k in l1["k"]})
    dims = [
        {
            "id": "k", "title": "位點數 k", "src": "topology", "srcLabel": "topology",
            "keys": [{"v": str(k), "label": (f"k={k}" if k < 8 else "k≥8"),
                      "color": K_COLORS[min(k, len(K_COLORS) - 1)]} for k in ks],
        },
        {
            "id": "unique", "title": "拓撲是否唯一解", "src": "topology", "srcLabel": "topology",
            "keys": [
                {"v": "all", "label": "所屬 region 全為唯一解", "color": "#176b58"},
                {"v": "some", "label": "部分 region 唯一解", "color": "#b66e20"},
                {"v": "none", "label": "皆非唯一解", "color": "#a94336"},
            ],
        },
        {
            "id": "ranked", "title": "unit 是否 ranked", "src": "topology", "srcLabel": "topology",
            "keys": [{"v": "yes", "label": "有 ranked region", "color": "#176b58"},
                     {"v": "no", "label": "無", "color": "#c9c8bd"}],
        },
        {
            "id": "multiRegion", "title": "是否跨多 region", "src": "topology", "srcLabel": "topology",
            "keys": [{"v": "yes", "label": "跨多 region", "color": "#6b5592"},
                     {"v": "no", "label": "單一 region", "color": "#c9c8bd"}],
        },
        {
            "id": "chrom", "title": "染色體", "src": "topology", "srcLabel": "topology",
            "keys": [{"v": c, "label": c} for c in l1["chroms"]],
        },
    ]
    dims.extend(extra or [])
    _ = regions
    return dims


LAYERS = [
    {"id": "density", "title": "1 Mb 密度帶", "on": True,
     "hint": "全顯示時擁擠處看不出多寡，底下補一條分箱灰階"},
    {"id": "ismRing", "title": "有 ISM 甲基資料的位點加圈", "on": True,
     "hint": "需要 ISM 能力；缺時此開關停用"},
    {"id": "treeCandidates", "title": "樹：候選聯集（藍）", "on": True,
     "hint": "需要 candidates 能力"},
    {"id": "treeMethyl", "title": "樹：甲基投影（橘）", "on": True,
     "hint": "事後投影，不參與拓撲推論"},
    {"id": "treeLineageLabel", "title": "樹：階層編號 HP2-1-1", "on": True,
     "hint": "需要 lineage 能力"},
]


def build_boot(sample: str, reg, dims: list, layers: list = None) -> dict:
    topo = reg.get("topology")
    return {
        "sample": sample,
        "l1": topo.payload["l1"],
        "capabilities": reg.matrix(),
        "dims": dims,
        "layers": layers if layers is not None else LAYERS,
    }


def write_shards(out_dir: Path, topo_cap, chroms: list, ism_cap=None) -> dict:
    """L2（region + 代表樹）與 L4（ISM 逐位點統計）逐染色體分片，同一個檔。

    分片而非全內嵌，是因為 build_exact_ps_layered_workstation 的 H2009 單頁
    188 MB 已證明全內嵌撐不住 genome-scale。"""
    data_dir = out_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    regions = topo_cap.payload["regions"]

    by_chrom = {}
    for r in regions:
        by_chrom.setdefault(r["c"], []).append(r)

    ism_rows = (ism_cap.payload or {}).get("rows", {}) if (ism_cap and ism_cap.usable) else {}
    ism_by_chrom = {}
    for (c, pos), rec in ism_rows.items():
        ism_by_chrom.setdefault(c, {})[str(pos)] = rec

    manifest = {}
    for c in chroms:
        rows = by_chrom.get(c, [])
        l4 = ism_by_chrom.get(c, {})
        path = data_dir / f"L2.{c}.js"
        body = "window.__DD=window.__DD||{};"
        body += "window.__DD.L2=window.__DD.L2||{};window.__DD.L4=window.__DD.L4||{};"
        body += f"window.__DD.L2[{json.dumps(c)}]={json_for_html(rows)};"
        body += f"window.__DD.L4[{json.dumps(c)}]={json_for_html(l4)};"
        path.write_text(body, encoding="utf-8")
        manifest[c] = {"file": f"data/L2.{c}.js", "regions": len(rows),
                       "ism_loci": len(l4), "bytes": path.stat().st_size}
    return manifest
