#!/usr/bin/env python3
"""共出現：維度之間的交集。

只放真的有結構的組合，不硬湊。每個矩陣都帶 Pearson 殘差（觀察−期望）/√期望，
因為單看 n 或列% 看不出「這個組合是不是超出隨機預期」。

n < 30 的格子不顯示百分比與殘差（NCHS 標準）—— 小樣本的比例會騙人。
"""
from __future__ import annotations

import math
from collections import Counter

MIN_N_FOR_PCT = 30


def _matrix(pairs, xs, ys):
    """pairs: [(x, y), ...]。回傳含計數、期望、Pearson 殘差的矩陣。"""
    cnt = Counter(pairs)
    total = sum(cnt.values())
    rowsum = Counter(); colsum = Counter()
    for (x, y), v in cnt.items():
        colsum[x] += v
        rowsum[y] += v

    cells = []
    for y in ys:
        row = []
        for x in xs:
            n = cnt.get((x, y), 0)
            exp = (colsum[x] * rowsum[y] / total) if total else 0
            resid = (n - exp) / math.sqrt(exp) if exp > 0 else None
            row.append({
                "n": n,
                "rowPct": (n / rowsum[y] * 100) if rowsum[y] >= MIN_N_FOR_PCT else None,
                "colPct": (n / colsum[x] * 100) if colsum[x] >= MIN_N_FOR_PCT else None,
                "resid": round(resid, 2) if resid is not None else None,
                "small": n < MIN_N_FOR_PCT,
            })
        cells.append(row)
    return {"xs": xs, "ys": ys, "cells": cells, "total": total,
            "colsum": [colsum[x] for x in xs], "rowsum": [rowsum[y] for y in ys]}


def build(reg) -> list:
    topo = reg.get("topology")
    regions = topo.payload["regions"]
    out = []

    # ── 1. k × 唯一解率 ────────────────────────────────────────────
    # 實測有明顯崖壁：k=1 96.5% → k=5 46.6% → k=6 驟降 4.5% → k≥7 全 0%
    by_k = {}
    for r in regions:
        k = r["k"]
        d = by_k.setdefault(k, [0, 0])
        d[0] += 1
        if r["tu"]:
            d[1] += 1
    rows = []
    for k in sorted(by_k):
        n, u = by_k[k]
        rows.append({"k": k, "n": n, "unique": u,
                     "rate": (u / n * 100) if n >= MIN_N_FOR_PCT else None,
                     "small": n < MIN_N_FOR_PCT})
    out.append({
        "id": "k_unique", "title": "位點數 k × 唯一解率",
        "kind": "bars", "rows": rows,
        "denom": f"分母 = {len(regions):,} 個 region",
        "src": "topology", "srcLabel": "topology",
        "note": "k=0 是「無 active ALT」，不是樹解不出來 —— 那些 region 根本沒有要解的問題。"
                "唯一解率隨 k 下降是預期的（候選空間變大），要看的是<b>下降有沒有崖壁</b>。",
        "caveat": "n < 30 的 k 不顯示比率（NCHS 標準），只給原始計數。",
    })

    # ── 2. unit_status × 是否有代表樹 ──────────────────────────────
    xs = sorted({r["us"] for r in regions})
    ys = ["有樹", "無樹"]
    m = _matrix([(r["us"], "有樹" if r["v"] else "無樹") for r in regions], xs, ys)
    out.append({
        "id": "status_tree", "title": "unit_status × 是否有代表樹",
        "kind": "matrix", "matrix": m,
        "denom": f"分母 = {len(regions):,} 個 region",
        "src": "topology", "srcLabel": "topology",
        "note": "這張圖應該是<b>完全塊狀</b>的 —— 每個 status 要嘛全有樹、要嘛全無樹。"
                "若出現混合格，代表 topology 的 status 與實際產物不一致，是資料問題。",
        "caveat": "",
    })

    # ── 3. k × 同分並列棵數 ────────────────────────────────────────
    def tie_bucket(r):
        try:
            t = int(r["tc"])
        except (TypeError, ValueError):
            return "未知"
        if t <= 1:
            return "1（唯一）"
        if t <= 3:
            return "2–3"
        if t <= 10:
            return "4–10"
        return "≥11"
    kx = [str(min(r["k"], 8)) for r in regions]
    ks = sorted({k for k in kx}, key=lambda s: int(s))
    ty = ["1（唯一）", "2–3", "4–10", "≥11", "未知"]
    m2 = _matrix([(str(min(r["k"], 8)), tie_bucket(r)) for r in regions], ks, ty)
    out.append({
        "id": "k_tie", "title": "位點數 k × 同分並列棵數",
        "kind": "matrix", "matrix": m2,
        "denom": f"分母 = {len(regions):,} 個 region",
        "src": "topology", "srcLabel": "topology",
        "note": "k 越大、同分並列越多是結構限制（候選空間隨 k 指數成長），"
                "所以<b>對角趨勢本身不是發現</b>。有訊息的是同一個 k 內的離散程度。",
        "caveat": "",
    })

    return out


def build_selection_sentinel(reg) -> dict:
    """ISM 可得性 × k —— 選擇偏誤前哨。

    若 ISM 可得性隨 k 變化，則任何「k × 甲基」的結論都被選擇偏誤污染。
    這張圖必須放在自檢層置頂，是所有甲基結論的前提檢查。
    """
    ism = reg.get("ism_dirs")
    topo = reg.get("topology")
    if not (ism and ism.usable and ism.payload and "hit" in (ism.payload or {})):
        return {
            "id": "ism_k_sentinel", "available": False,
            "title": "ISM 可得性 × k（選擇偏誤前哨）",
            "reason": "缺 ISM 能力。在接上之前，任何「拓撲 × 甲基」的結論都無從檢查是否被選擇偏誤污染。",
        }
    l1 = topo.payload["l1"]
    hit = ism.payload["hit"]          # set of L1 index
    xs = sorted({str(min(k, 8)) for k in l1["k"]}, key=lambda s: int(s))
    pairs = [(str(min(l1["k"][i], 8)), "有 ISM" if i in hit else "無 ISM")
             for i in range(l1["n"])]
    m = _matrix(pairs, xs, ["有 ISM", "無 ISM"])
    return {
        "id": "ism_k_sentinel", "available": True,
        "title": "ISM 可得性 × k（選擇偏誤前哨）",
        "kind": "matrix", "matrix": m,
        "denom": f"分母 = {l1['n']:,} 個 sSNV",
        "src": "derived", "srcLabel": "儀表板自算",
        "note": "ISM 的覆蓋不是隨機的 —— 缺的位點是被 ISM 自己的 coverage / CpG gate 濾掉的，"
                "與 CpG 密度、覆蓋深度相關。<b>若各 k 的可得性明顯不同，"
                "則所有「k × 甲基」的結論都被選擇偏誤污染</b>，必須在結論裡聲明母體。",
        "caveat": "此矩陣的母體是全部 sSNV；下游甲基面板的母體則只有「通過 ISM gate」那批。",
    }
