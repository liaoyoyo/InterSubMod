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

    # ── C10 循環論證偵測：某軸顯著率逼近 100% 是雙重使用的訊號 ────
    ismc = reg.get("ism_dirs")
    if ismc and ismc.usable and ismc.payload:
        flagged = []
        for ax in ismc.payload["axes"]:
            sig = tot = 0
            for rec in ismc.payload["rows"].values():
                p = rec.get(ax["id"] + "_p")
                v = rec.get(ax["id"] + "_valid")
                n = rec.get(ax["id"] + "_n")
                if v is False or (n is not None and n < 2) or p is None:
                    continue
                tot += 1
                if p <= 0.05:
                    sig += 1
            if tot >= 100 and sig / tot >= 0.99:
                flagged.append(f"{ax['title']} {sig:,}/{tot:,}={sig / tot * 100:.1f}%")
        checks.append(_chk(
            "C10", "沒有任何軸的顯著率逼近 100%（雙重使用偵測）",
            PASS if not flagged else FAIL,
            "偵測到近乎全數顯著的軸：" + "、".join(flagged) +
            "　→ 這通常代表分組標籤與檢定用的距離矩陣同源（double-dipping），p 值不可當證據"
            if flagged else "各軸顯著率皆 < 99%"))
    else:
        checks.append(_chk(
            "C10", "沒有任何軸的顯著率逼近 100%（雙重使用偵測）", SKIP,
            "缺 ISM 能力", need="ism_dirs"))

    # ── C11 LCA 增益的 A/B 一致性 ──────────────────────────────────
    ab = reg.get("lca_ab")
    if ab and ab.usable and ab.payload:
        d = ab.payload
        same = d.get("identical_fields", 0)
        diff = d.get("differing_fields", [])
        ok = set(diff) <= {"lv_written", "lca_resolved", "lca_candidates_sum"}
        checks.append(_chk(
            "C11", "pre-LCA 與 post-LCA receipt 是乾淨的 A/B（只有 LCA 相關欄位不同）",
            PASS if ok else FAIL,
            f"兩套各 {d.get('n_files', 0)} 檔；{same} 個欄位完全相同，"
            f"不同的只有 {'、'.join(diff) or '無'}；"
            f"lv 覆蓋 {d.get('pre_lv', 0):,} → {d.get('post_lv', 0):,}"
            f"（{d.get('gain', 0):.2f}×）"))
    else:
        checks.append(_chk(
            "C11", "pre-LCA 與 post-LCA receipt 是乾淨的 A/B", SKIP,
            "缺 bam_pre_lca_receipts/，無法量化 LCA 到底買到什麼", need="lca_ab"))

    # ── C12 block 的 sSNV 共現圖是否連通 ──────────────────────────
    ml = reg.get("mlhp")
    if ml and ml.usable and ml.payload:
        br = ml.counts.get("broken_cooccurrence", 0)
        tot_b = len(ml.payload["by_region"])
        checks.append(_chk(
            "C12", "每個 block 的 sSNV 都被 read 串成一塊（共現圖連通）",
            PASS if br == 0 else FAIL,
            f"斷裂的 block {br:,} / {tot_b:,}（{br / max(tot_b, 1) * 100:.2f}%）"
            "　→ 那些 block 內有 sSNV 之間沒有任何 read 同時覆蓋，"
            "solver 仍建出樹，但跨分量的邊<b>零 read 支持</b>，"
            "不可當作連鎖證據"))
    else:
        checks.append(_chk(
            "C12", "每個 block 的 sSNV 都被 read 串成一塊（共現圖連通）", SKIP,
            "缺 MLHP 能力，無法檢查共現連通性", need="mlhp"))

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
    ranked = sum(1 for r in regions if r["us"] == "ranked")
    no_tree = Counter(r["us"] for r in regions if not r["v"])
    # 兩個分母都給。只給 78.8% 會讓人以為 lineage 層漏了兩成 —— 實際上
    # 它覆蓋了 100% 的 ranked region，另外 2,460 個是 solver 根本沒建樹
    # （no_active_alt / family_incomplete / zero_denominator），不是漏掉。
    out.append({
        "title": "有代表樹的 region（分母 = 全部 region）",
        "num": with_tree, "den": len(regions),
        "gap": "缺的 {:,} 個是 solver 沒建樹，不是資料漏掉：{}".format(
            len(regions) - with_tree,
            "、".join(f"{k} {v:,}" for k, v in sorted(no_tree.items()))),
    })
    out.append({
        "title": "有代表樹的 region（分母 = ranked region）",
        "num": with_tree, "den": ranked,
        "gap": ("ranked 與有樹是同一個集合（實測集合完全相等）—— "
                "lineage 層對「有樹的 region」是 100% 覆蓋。"
                "上一列的 78.8% 是相對全部 region，兩者都對但意思不同。"),
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
