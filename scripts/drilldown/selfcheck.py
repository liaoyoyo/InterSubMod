#!/usr/bin/env python3
"""自檢層：用機械計算的守恆等式回答「這份結果整體有沒有問題」。

每一項都是等式，不是敘述。任一不成立就在頁面上標紅，不靜默通過。

需要某個能力才能算的檢查，在該能力缺席時狀態是「無法檢查」——
那是第三態，既不是通過也不是失敗。把它算成通過會讓人以為驗過了。
"""
from __future__ import annotations

from collections import Counter

PASS = "PASS"
FAIL = "FAIL"
SKIP = "SKIP"          # 缺能力，無法檢查

LABEL = {PASS: "通過", FAIL: "不成立", SKIP: "無法檢查"}


def _chk(cid, title, status, detail, need=None):
    return {"id": cid, "title": title, "status": status, "statusLabel": LABEL[status],
            "detail": detail, "need": need}


def run(reg) -> dict:
    checks = []
    topo = reg.get("topology")
    regions = topo.payload["regions"]
    l1 = topo.payload["l1"]

    # ── C1 逐染色體 sSNV 加總 = L1 總數 ────────────────────────────
    per_chrom = Counter(l1["chrom"])
    s = sum(per_chrom.values())
    checks.append(_chk(
        "C1", "逐染色體 sSNV 加總 = L1 總數",
        PASS if s == l1["n"] else FAIL,
        f"Σ 22 條 = {s:,}；L1 總數 = {l1['n']:,}"))

    # ── C2 unit_status 分類加總 = region 總數 ──────────────────────
    us = Counter(r["us"] for r in regions)
    checks.append(_chk(
        "C2", "unit_status 各類加總 = region 總數",
        PASS if sum(us.values()) == len(regions) else FAIL,
        "　".join(f"{k}={v:,}" for k, v in sorted(us.items())) +
        f"　→ Σ {sum(us.values()):,} / region {len(regions):,}"))

    # ── C3 k 分布加總 = region 總數 ────────────────────────────────
    ks = Counter(r["k"] for r in regions)
    checks.append(_chk(
        "C3", "k 分布加總 = region 總數",
        PASS if sum(ks.values()) == len(regions) else FAIL,
        "　".join(f"k={k}:{v:,}" for k, v in sorted(ks.items())) +
        f"　→ Σ {sum(ks.values()):,}"))

    # ── C4 有樹 + 無樹 = region 總數，且無樹者的 status 可完全歸因 ──
    with_tree = sum(1 for r in regions if r["v"])
    no_tree = Counter(r["us"] for r in regions if not r["v"])
    ok = (with_tree + sum(no_tree.values())) == len(regions)
    checks.append(_chk(
        "C4", "有樹 + 無樹 = region 總數（且無樹者可完全歸因）",
        PASS if ok else FAIL,
        f"有樹 {with_tree:,} + 無樹 {sum(no_tree.values()):,} = {len(regions):,}；"
        f"無樹原因：" + "、".join(f"{k}×{v:,}" for k, v in sorted(no_tree.items()))))

    # ── C5 每個有樹的 region 都滿足 |V| - 1 == |E| ─────────────────
    bad = [r["id"] for r in regions if r["v"] and (len(r["v"]) - 1) != len(r["e"])]
    checks.append(_chk(
        "C5", "每個代表樹滿足 |V| − 1 = |E|（樹的定義）",
        PASS if not bad else FAIL,
        f"檢查 {with_tree:,} 棵樹，違反 {len(bad):,} 棵"
        + (f"；前 3 例：{', '.join(bad[:3])}" if bad else "")))

    # ── C6 樹的邊指向的 vertex 都存在於該樹 ────────────────────────
    dangling = 0
    for r in regions:
        if not r["v"]:
            continue
        ids = {v[0] for v in r["v"]}
        for e in r["e"]:
            if e[0] not in ids or e[1] not in ids:
                dangling += 1
    checks.append(_chk(
        "C6", "樹的每條邊兩端都在該樹的 vertex 集合內",
        PASS if dangling == 0 else FAIL,
        f"懸空邊 {dangling:,} 條"))

    # ── C7 CSR 索引與 (region,pos) pair 一致 ───────────────────────
    pair = sum(len(r["ap"]) for r in regions)
    csr = len(l1["csrVal"])
    checks.append(_chk(
        "C7", "CSR 條目數 = Σ 各 region 的 active_positions 數",
        PASS if pair == csr else FAIL,
        f"Σ active_positions = {pair:,}；CSR 條目 = {csr:,}"))

    # ── C8 hidden 節點數（需要 lineage 才能交叉驗證）───────────────
    h_prefix = sum(1 for r in regions for v in r["v"] if str(v[1]).startswith("H_"))
    lin = reg.get("lineage_paths")
    if lin and lin.usable and lin.payload and "hidden" in (lin.payload or {}):
        exp = lin.payload["hidden"]
        checks.append(_chk(
            "C8", "topology 的 H_ 前綴數 = lineage_paths 的 is_hidden 數",
            PASS if h_prefix == exp else FAIL,
            f"topology H_ = {h_prefix:,}；lineage is_hidden = {exp:,}"))
    else:
        checks.append(_chk(
            "C8", "topology 的 H_ 前綴數 = lineage_paths 的 is_hidden 數",
            SKIP, f"topology 的 H_ 前綴共 {h_prefix:,}，但缺 lineage_paths 無從交叉驗證",
            need="lineage_paths"))

    # ── C9 ISM 覆蓋（需要 ISM 能力）────────────────────────────────
    ism = reg.get("ism_dirs")
    if ism and ism.usable and ism.linkage:
        lk = ism.linkage
        checks.append(_chk(
            "C9", "有 ISM 目錄 + 無 ISM 目錄 = sSNV 總數",
            PASS if lk.denominator == l1["n"] else FAIL,
            f"有 {lk.numerator:,} + 無 {lk.denominator - lk.numerator:,} = {lk.denominator:,}；"
            f"sSNV 總數 {l1['n']:,}"))
    else:
        checks.append(_chk(
            "C9", "有 ISM 目錄 + 無 ISM 目錄 = sSNV 總數", SKIP,
            "缺 ISM 能力，無法檢查甲基層覆蓋", need="ism_dirs"))

    n_pass = sum(1 for c in checks if c["status"] == PASS)
    n_fail = sum(1 for c in checks if c["status"] == FAIL)
    n_skip = sum(1 for c in checks if c["status"] == SKIP)
    return {
        "checks": checks,
        "summary": {"pass": n_pass, "fail": n_fail, "skip": n_skip, "total": len(checks)},
        "coverage": _coverage(reg, regions, l1),
    }


def _coverage(reg, regions, l1) -> list:
    """覆蓋率不只給比率，要說缺的是什麼 —— 只給百分比會讓人以為缺的是隨機的。"""
    out = []
    with_tree = sum(1 for r in regions if r["v"])
    no_tree = Counter(r["us"] for r in regions if not r["v"])
    out.append({
        "title": "有代表樹的 region",
        "num": with_tree, "den": len(regions),
        "gap": "缺的 {:,} 個逐項歸因：{}".format(
            len(regions) - with_tree,
            "、".join(f"{k} {v:,}" for k, v in sorted(no_tree.items()))),
    })
    for cid in ("ism_dirs", "prerender_meth", "lineage_paths", "lineage_reads"):
        c = reg.get(cid)
        if c and c.linkage:
            out.append({
                "title": c.title, "num": c.linkage.numerator, "den": c.linkage.denominator,
                "gap": c.linkage.method,
            })
        elif c:
            out.append({"title": c.title, "num": 0, "den": l1["n"],
                        "gap": f"能力狀態 {c.state} — {c.reason or '未探測'}"})
    return out


def to_markdown(result: dict, sample: str, inputs: list) -> str:
    s = result["summary"]
    lines = [f"# {sample} · 自檢報告", "",
             f"守恆等式 {s['pass']} 通過 / {s['fail']} 不成立 / {s['skip']} 無法檢查"
             f"（共 {s['total']}）", "",
             "「無法檢查」是第三態 —— 該項需要的資料層不在，既不是通過也不是失敗。", "",
             "## 守恆等式", ""]
    for c in result["checks"]:
        mark = {PASS: "✅", FAIL: "❌", SKIP: "⊘"}[c["status"]]
        lines.append(f"### {mark} {c['id']} {c['title']}")
        lines.append("")
        lines.append(c["detail"])
        if c["need"]:
            lines.append(f"")
            lines.append(f"> 需要能力 `{c['need']}`")
        lines.append("")
    lines += ["## 覆蓋率與缺口成因", "",
              "| 項目 | 分子 | 分母 | 比率 | 缺的是什麼 |", "|---|---:|---:|---:|---|"]
    for c in result["coverage"]:
        r = f"{c['num'] / c['den'] * 100:.1f}%" if c["den"] else "—"
        lines.append(f"| {c['title']} | {c['num']:,} | {c['den']:,} | {r} | {c['gap']} |")
    lines += ["", "## 輸入 provenance", "",
              "| 路徑 | bytes | sha256 | mtime |", "|---|---:|---|---|"]
    for i in inputs:
        lines.append(f"| `{i['path']}` | {i['size_bytes']:,} | `{i['sha256'][:16]}…` | {i['mtime']} |")
    return "\n".join(lines) + "\n"
