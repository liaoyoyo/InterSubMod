#!/usr/bin/env python3
"""strict endpoint edges（擴充能力）：block 內 sSNV 兩兩之間的真實 read 連鎖證據。

這是 component 形成的**權威資料源**（region_linkage_rule =
"strict fixed-R/A endpoint-pair support connected component"）。

🔴 為什麼一定要用這個而不是 MLHP 的 pattern 表：
   subread_groups_by_hp 是「exact-PS projected molecule support after
   O/D/S/L/X→X marginalization」—— 跨越兩個位點的 read 若在其他位點是
   O/D/S/L，投影後整條會出現 X，pattern 表因此看不到那條連鎖。
   曾用 pattern 表算連通性，得到「16.16% 的 block 斷裂」的錯誤結論；
   用本檔的邊重算是 **0 / 11,590**，全部連通，設計本來就是對的。

每條邊帶：support_total、RR/RA/AR/AA 拆解、gap_bp、thresholds_passed。
thresholds_passed 讓我們能回答「門檻改成 T 會多／少多少邊」而不必重跑管線。
"""
from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from pathlib import Path

from capability import Capability, Registry, finalize, make, probe, probe_linkage

REQUIRED = {"chrom", "hp_family", "phase_set", "pos_i1", "pos_j1",
            "support_total", "passes_primary_threshold"}


def load(reg: Registry, root, regions, mlhp_cap=None) -> Capability:
    cap = make(reg, "strict_edges", "strict endpoint edges（read 連鎖證據）",
               enables=["edge_support", "selfcheck:C12", "threshold_whatif"])
    root = Path(root) if root else None
    if not root or not root.is_dir():
        cap.state = "ABSENT"
        probe(cap, "S0", False, f"chromosomes 目錄不存在：{root}")
        return cap

    files = sorted(root.glob("chr*/strict_regions/*.endpoint_edges.tsv.gz"))
    if not files:
        cap.state = "ABSENT"
        probe(cap, "S0", False, f"{root} 底下找不到 */strict_regions/*.endpoint_edges.tsv.gz")
        return cap
    probe(cap, "S0", True, f"{len(files)} 個染色體的 endpoint_edges")

    by_ps = defaultdict(dict)      # (chrom, hp, ps) -> {(a,b): edge}
    n_edge = n_pass = 0
    thresholds = defaultdict(int)
    with gzip.open(files[0], "rt", newline="") as fh:
        head = set((csv.DictReader(fh, delimiter="\t").fieldnames or []))
    missing = sorted(REQUIRED - head)
    if missing:
        cap.state = "MALFORMED"
        probe(cap, "S1", False, f"缺必要欄位：{', '.join(missing)}")
        return cap
    probe(cap, "S1", True, f"{len(head)} 欄，必要欄位齊備")

    for f in files:
        with gzip.open(f, "rt", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                n_edge += 1
                for t in (r.get("thresholds_passed") or "").split(","):
                    if t.strip().isdigit():
                        thresholds[int(t)] += 1
                if r.get("passes_primary_threshold") != "true":
                    continue
                n_pass += 1
                a, b = int(r["pos_i1"]), int(r["pos_j1"])
                if a > b:
                    a, b = b, a
                by_ps[(r["chrom"], str(r["hp_family"]), str(r["phase_set"]))][(a, b)] = {
                    "n": int(r.get("support_total") or 0),
                    "rr": int(r.get("support_RR") or 0),
                    "ra": int(r.get("support_RA") or 0),
                    "ar": int(r.get("support_AR") or 0),
                    "aa": int(r.get("support_AA") or 0),
                    "gap": int(r.get("gap_bp") or 0),
                }
    probe(cap, "S2", True,
          f"邊 {n_edge:,}；通過主門檻 {n_pass:,}（{n_pass / max(n_edge, 1) * 100:.1f}%）")
    probe(cap, "S2", True,
          "門檻 what-if：" + "　".join(
              f"T≥{t} {thresholds[t]:,}（{thresholds[t] / max(thresholds.get(3, 1), 1):.2f}×）"
              for t in sorted(thresholds)))

    # S3：逐 block 檢查連通性 —— 這才是正確的檢查方式
    # ⚠ 用**全部** positions（含非 active），不是 topology 的 active_positions。
    #   component 是用全部位點的邊形成的；若橋接兩個 active 位點的是一個
    #   非 active 位點，只看 active 子集會誤判成斷裂（實測差 394 個 block）。
    mlhp_pos = {}
    if mlhp_cap and mlhp_cap.usable and mlhp_cap.payload:
        for rid, m in mlhp_cap.payload["by_region"].items():
            if m.get("pos"):
                mlhp_pos[rid] = m["pos"]

    broken, checked, no_edge = 0, 0, 0
    broken_active = 0
    broken_ids = []
    for r in regions:
        pos = mlhp_pos.get(r["id"]) or r.get("ap") or []
        if len(pos) < 2:
            continue
        checked += 1
        E = by_ps.get((r["c"], str(r["hp"]), str(r["ps"])))
        if not E:
            no_edge += 1
            continue
        S = set(pos)
        adj = {p: set() for p in pos}
        for (a, b) in E:
            if a in S and b in S:
                adj[a].add(b)
                adj[b].add(a)
        seen, nc = set(), 0
        for p in pos:
            if p in seen:
                continue
            nc += 1
            st = [p]
            while st:
                x = st.pop()
                if x in seen:
                    continue
                seen.add(x)
                st.extend(adj[x] - seen)
        if nc > 1:
            broken += 1
            if len(broken_ids) < 5:
                broken_ids.append(r["id"])

        # 另外記錄：只看 active 子集時會不會斷。那不是資料問題，
        # 是「橋接者是非 active 位點」，但畫樹時只有 active 進去，值得標明。
        act = [p for p in (r.get("ap") or []) if p in S]
        if len(act) >= 2:
            sa = set(act)
            seen2, nc2 = set(), 0
            for p in act:
                if p in seen2:
                    continue
                nc2 += 1
                st = [p]
                while st:
                    x = st.pop()
                    if x in seen2:
                        continue
                    seen2.add(x)
                    st.extend((adj[x] & sa) - seen2)
            if nc2 > 1:
                broken_active += 1

    probe_linkage(cap, checked - broken, checked,
                  "k≥2 的 block 其**全部** sSNV 由通過門檻的邊連通")
    probe(cap, "S3", True,
          f"只看 active 子集時斷裂的 block：{broken_active:,} / {checked:,}"
          f"（{broken_active / max(checked, 1) * 100:.2f}%）—— "
          f"那不是資料問題，是橋接者為非 active 位點；畫樹時只有 active 進去")

    cap.counts.update({"edges": n_edge, "edges_passed": n_pass,
                       "blocks_checked": checked, "blocks_broken": broken,
                       "blocks_broken_active_only": broken_active,
                       "blocks_no_edge": no_edge,
                       "thresholds": dict(thresholds)})
    cap.payload = {"by_ps": by_ps, "thresholds": dict(thresholds),
                   "broken": broken, "checked": checked, "broken_ids": broken_ids,
                   "broken_active": broken_active,
                   "n_edge": n_edge, "n_pass": n_pass}
    return finalize(cap)
