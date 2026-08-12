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


# 甲基軸的顏色語意。刻意不含「甲基自身分群」——那個軸是循環論證（見 C10），
# 把它畫上去會讓幾乎每個位點都亮起來，等於沒有資訊。
AXIS_BITS = [("hp", 1, "HP", "#285f8f"), ("allele", 2, "ALT/REF", "#a94336"),
             ("lineage", 4, "lineage", "#6b5592")]
AXIS_NONE = 8          # 有測但都不顯著
AXIS_UNTESTED = 16     # 有 ISM summary 列，但三個全域軸都無法檢定

# 全域軸與 detail-only 軸必須明確分開。舊 standalone 頁會同時列出
# allele / HP 軸，但沒有一個可多樣本重現的定義表；改由 payload 作 SoT。
AXIS_DEFINITIONS = [
    {
        "id": "hp", "label": "HP（merged）", "field": "HPMergedP",
        "scope": "全域軸道 + 位點細節",
        "definition": "比較 merged haplotype 分組的甲基距離；raw p ≤ 0.05 只作探索標記。",
    },
    {
        "id": "allele", "label": "ALT / REF", "field": "AlleleP",
        "scope": "全域軸道 + 位點細節",
        "definition": "比較 ALT-support 與 REF-support read 的甲基距離；是關聯，不是因果。",
    },
    {
        "id": "lineage", "label": "lineage", "field": "LineagePermanovaP",
        "scope": "全域軸道 + 位點細節",
        "definition": "比較 read lineage 分組；lineage vertex 不等於 clone / CCF。",
    },
    {
        "id": "hpfine", "label": "HP（fine）", "field": "HPFineP",
        "scope": "只在位點細節",
        "definition": "較細 HP 分組，與 merged HP 相關；不重複編入全域單一類別軸。",
    },
    {
        "id": "cluster", "label": "甲基自身分群", "field": "ClusterPermanovaP",
        "scope": "只在位點細節，不納入結論",
        "definition": "分群標籤與檢定使用同一距離矩陣，屬循環論證軸。",
    },
]


def fill_axis_codes(l1: dict, ism_cap) -> None:
    """逐 sSNV 算一個 byte：沿哪些軸顯著。放進 L1 讓全基因組圖能直接畫，
    不必等分片載入。"""
    n = l1["n"]
    codes = [0] * n
    if not (ism_cap and ism_cap.usable and ism_cap.payload):
        l1["axis"] = codes
        return
    rows = ism_cap.payload["rows"]
    chroms = l1["chroms"]
    acc, prev = 0, -1
    for i in range(n):
        if l1["chrom"][i] != prev:
            acc = l1["dpos"][i]
            prev = l1["chrom"][i]
        else:
            acc += l1["dpos"][i]
        rec = rows.get((chroms[l1["chrom"][i]], acc))
        if rec is None:
            continue
        code, tested = 0, False
        for aid, bit, _lab, _col in AXIS_BITS:
            p = rec.get(aid + "_p")
            v = rec.get(aid + "_valid")
            k = rec.get(aid + "_n")
            if v is False or (k is not None and k < 2) or p is None:
                continue
            tested = True
            if p <= 0.05:
                code |= bit
        if tested and code == 0:
            code = AXIS_NONE
        elif not tested:
            # 不能與「沒有 summary 列」共用 0，否則 missingness 會被誤畫成
            # ISM 完全缺席。
            code = AXIS_UNTESTED
        codes[i] = code
    l1["axis"] = codes


def observation_scope(l1: dict, ism_cap) -> dict:
    """全基因組觀察的分母與 missingness ledger。

    這是從舊 standalone 頁「總體 → 符合條件 → 頁面展示」漏斗汲取的可泛化
    部分；不帶入該頁尚未驗證的科學結論或 HCC 常數。
    """
    n = int(l1.get("n", 0) or 0)
    defs = []
    available_axes = set()
    if ism_cap and ism_cap.usable and ism_cap.payload:
        available_axes = {a.get("id") for a in (ism_cap.payload.get("axes") or [])}
    for row in AXIS_DEFINITIONS:
        item = dict(row)
        item["available_in_run"] = row["id"] in available_axes
        defs.append(item)

    base = {
        "claim_ceiling": "OBSERVATION_ONLY_NO_TRUTH_SET",
        "universe": {"num": n, "den": n, "unit": "sSNV", "available": True,
                     "label": "topology 觀察母體"},
        "axis_definitions": defs,
        "alpha": 0.05,
        "multiplicity": "raw p ≤ 0.05；未在此全域觀察層做 cohort-wide 多重檢定校正",
        "legacy_crosswalk": ("本頁不把舊 standalone 的 A_ALLELE_clean / HP 類別映射到現行軸；"
                             "兩者母體、檢定與分類 contract 不同，只沿用其漸進展開與"
                             "明示選案的資訊架構。"),
        "missingness": ("ISM 可得性受 coverage / CpG gate 影響，不可假設為隨機缺失；"
                        "必須搭配自檢頁的「ISM 可得性 × k」前哨。"),
    }
    if not (ism_cap and ism_cap.usable and ism_cap.payload):
        reason = (ism_cap.reason if ism_cap else "未提供 ISM 能力") or "ISM 能力不可用"
        base.update({
            "available": False, "reason": reason,
            "summary_join": {"available": False, "num": None, "den": n},
            "axis_testable": {"available": False, "num": None, "den": n},
            "axis_significant_raw": {"available": False, "num": None, "den": None},
            "detail_available": {"available": False, "num": None, "den": n},
        })
        return base

    codes = list(l1.get("axis") or [0] * n)
    summary_n = sum(1 for c in codes if c != 0)
    testable_n = sum(1 for c in codes if c not in (0, AXIS_UNTESTED))
    sig_mask = sum(bit for _aid, bit, _lab, _col in AXIS_BITS)
    sig_n = sum(1 for c in codes if c & sig_mask)
    detail_n = len(ism_cap.payload.get("hit_dir") or set())
    base.update({
        "available": True,
        "summary_join": {"available": True, "num": summary_n, "den": n,
                         "label": "有 significance summary 列"},
        "axis_testable": {"available": True, "num": testable_n, "den": summary_n,
                          "label": "至少一個預指定非循環軸可檢定"},
        "axis_significant_raw": {"available": True, "num": sig_n, "den": testable_n,
                                 "label": "至少一個全域軸 raw p ≤ 0.05"},
        "detail_available": {"available": True, "num": detail_n, "den": n,
                             "label": "有逐位點細節目錄"},
    })
    return base


def axis_dim(ism_cap) -> list:
    """甲基軸做成篩選維度 —— 使用者要能「只看沿 lineage 分群的位點」。"""
    if not (ism_cap and ism_cap.usable):
        return [{
            "id": "axis", "title": "甲基沿哪個軸顯著", "src": "ism", "srcLabel": "ISM",
            "disabled": True, "disabledReason": "缺 ISM 能力",
            "keys": [{"v": "none", "label": "（不可用）"}],
        }]
    keys = [{"v": "0", "label": "無 ISM summary 列", "color": "#e7e5da"},
            {"v": str(AXIS_UNTESTED), "label": "有 summary，但全域軸皆未檢定",
             "color": "#d7a85b"},
            {"v": str(AXIS_NONE), "label": "至少一軸有檢定，但都不顯著",
             "color": "#c9c8bd"}]
    keys += [{"v": str(bit), "label": lab + " 顯著", "color": col}
             for _aid, bit, lab, col in AXIS_BITS]
    keys.append({"v": "multi", "label": "多軸同時顯著", "color": "#17231e"})
    return [{
        "id": "axis", "title": "甲基沿哪個軸顯著", "src": "ism", "srcLabel": "ISM",
        "keys": keys,
        "note": ("全域只編碼 HP merged、ALT/REF、lineage；HP fine 留在位點細節。"
                 "不含「甲基自身分群」——該軸是循環論證，見自檢 C10。"
                 "此處使用 raw p ≤ 0.05，僅供探索式下鑽。"),
    }]


LAYERS = [
    {"id": "density", "title": "1 Mb 密度帶", "on": False,
     "hint": "擁擠處看不出多寡時可開；與甲基軸軌道同時開會偏擠，故預設關"},
    {"id": "axisTrack", "title": "甲基軸軌道（顯著軸以顏色標示）", "on": True,
     "hint": "每條染色體下方一條細軌，顏色 = 該位點的甲基沿哪個軸顯著"},
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


def write_shards(out_dir: Path, topo_cap, chroms: list, ism_cap=None, mlhp_cap=None,
                 igv_fn=None) -> dict:
    """L2（region + 代表樹）與 L4（ISM 逐位點統計）逐染色體分片，同一個檔。

    分片而非全內嵌，是因為 build_exact_ps_layered_workstation 的 H2009 單頁
    188 MB 已證明全內嵌撐不住 genome-scale。"""
    data_dir = out_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    regions = topo_cap.payload["regions"]

    by_chrom = {}
    for r in regions:
        by_chrom.setdefault(r["c"], []).append(r)

    mlhp = (mlhp_cap.payload or {}).get("by_region", {}) if (mlhp_cap and mlhp_cap.usable) else {}
    if mlhp:
        for r in regions:
            m = mlhp.get(r["id"])
            if m:
                r["pat"] = m          # read pattern 分布，read state matrix 用

    ism_rows = (ism_cap.payload or {}).get("rows", {}) if (ism_cap and ism_cap.usable) else {}
    ism_by_chrom = {}
    for (c, pos), rec in ism_rows.items():
        ism_by_chrom.setdefault(c, {})[str(pos)] = rec

    manifest = {}
    for c in chroms:
        rows = by_chrom.get(c, [])
        l4 = ism_by_chrom.get(c, {})
        # IGV 圖**不進分片** —— 每張 SVG 約 73 KB，5,410 張全內嵌會讓
        # 分片爆到 1.4 GB（chr7 單片 138 MB），正是要避免的 H2009 失敗模式。
        # 改成逐 region 獨立 .js，點到該 region 才載入；分片只留一個布林旗標。
        if igv_fn:
            for r in rows:
                v = igv_fn(c, r)
                if v:
                    r["igvOk"] = 1
        path = data_dir / f"L2.{c}.js"
        body = "window.__DD=window.__DD||{};"
        body += "window.__DD.L2=window.__DD.L2||{};window.__DD.L4=window.__DD.L4||{};"
        body += f"window.__DD.L2[{json.dumps(c)}]={json_for_html(rows)};"
        body += f"window.__DD.L4[{json.dumps(c)}]={json_for_html(l4)};"
        path.write_text(body, encoding="utf-8")
        manifest[c] = {"file": f"data/L2.{c}.js", "regions": len(rows),
                       "ism_loci": len(l4), "bytes": path.stat().st_size}
    return manifest
