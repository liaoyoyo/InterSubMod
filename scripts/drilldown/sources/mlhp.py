#!/usr/bin/env python3
"""MLHP（擴充能力）：提供 read state matrix 需要的 pattern 分布。

topology.jsonl 只有樹與 AF，沒有「哪些 read 觀察到哪個 R/A/X pattern」。
那在 MLHP 的 groups[] 裡：
    populations_by_hp      {"2": {"AA": 4}}            完整覆蓋的 pattern
    subread_groups_by_hp   {"2": {"AX":3,"XA":32,...}} 部分覆蓋的 pattern
    n_full_cov_reads / reads_by_hp / projected_molecule_rows

沒有這層，read state matrix 就畫不出來 —— 面板會顯示「不可用 + 原因」而不是空白。
"""
from __future__ import annotations

import json
from pathlib import Path

from capability import (MALFORMED, Capability, Registry, finalize, make, probe,
                        probe_exists, probe_linkage)


def load(reg: Registry, path, region_ids: set, expected_sample=None) -> Capability:
    cap = make(reg, "mlhp", "MLHP（read pattern 分布）",
               enables=["read_state_matrix", "pattern_stats"])
    if not path or not Path(path).is_file():
        cap.state = "ABSENT"
        probe(cap, "S0", False, f"MLHP 不存在：{path}")
        return cap
    if not probe_exists(cap, path):
        return cap

    try:
        doc = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        cap.state = "MALFORMED"
        probe(cap, "S1", False, f"無法解析：{type(exc).__name__}: {exc}")
        return cap

    if expected_sample:
        observed = str(doc.get("sample") or doc.get("sample_id") or "").strip()
        if observed != str(expected_sample):
            cap.state = MALFORMED
            probe(cap, "S1", False,
                  f"MLHP sample 不符：--sample={expected_sample}，"
                  f"文件 sample={observed or '(缺)'}；為避免跨樣本混用已禁用此層")
            return cap
        probe(cap, "S1", True, f"sample 身分一致：{observed}")

    groups = doc.get("groups")
    if not isinstance(groups, list) or not groups:
        cap.state = "MALFORMED"
        probe(cap, "S1", False, "缺 groups[] 或為空")
        return cap
    need = {"region_id", "populations_by_hp", "subread_groups_by_hp", "positions"}
    missing = sorted(need - set(groups[0]))
    if missing:
        cap.state = "MALFORMED"
        probe(cap, "S1", False, f"groups[0] 缺欄位：{', '.join(missing)}")
        return cap
    probe(cap, "S1", True, f"groups {len(groups):,} 筆，必要欄位齊備")

    def _components(pats, k):
        """把 pattern 轉成 sSNV 共現圖，回傳 (分量數, 各分量大小)。

        🔴 這裡算出的「斷裂」**不等於**「沒有 read 串連」——曾經這樣解讀，錯了。

        groups 是 exact-PS 投影後的 pattern 表：read 在任一位點是 O/D/S/L 時，
        投影會讓整條變成 X。因此一條真的橫跨兩個位點的 read，只要在別的位點
        落在那些狀態，就不會出現在 pattern 表的共現裡。用這張表判連通性會
        系統性低估，當初據此得出「16.16% 的 block 斷裂、跨分量的樹邊零 read
        支持」的結論，實測是錯的。

        判連通性的權威來源是 endpoint_edges（見 sources/strict_edges.py）：
        用通過門檻的邊重算，11,590 / 11,590 全部連通，設計本來就是對的。

        本函式保留下來，是因為「投影後的 pattern 表看不到多少連鎖」本身
        是有意義的統計量（它決定 read state matrix 能顯示什麼），
        但它**只描述投影表**，不可外推成 read 支持度的結論。
        """
        if k < 2:
            return 1, [k]
        adj = {i: set() for i in range(k)}
        for pat in pats:
            idx = [i for i, ch in enumerate(pat) if ch != "X"]
            for a in range(len(idx)):
                for b in range(a + 1, len(idx)):
                    adj[idx[a]].add(idx[b])
                    adj[idx[b]].add(idx[a])
        seen, comps = set(), []
        for i in range(k):
            if i in seen:
                continue
            stack, cur = [i], set()
            while stack:
                x = stack.pop()
                if x in seen:
                    continue
                seen.add(x)
                cur.add(x)
                stack.extend(adj[x] - seen)
            comps.append(sorted(cur))
        return len(comps), comps

    by_region = {}
    n_broken = 0
    for g in groups:
        rid = g.get("region_id")
        if not rid:
            continue
        _pats = {}
        for _hp, _m in (g.get("populations_by_hp") or {}).items():
            _pats.update(_m)
        for _hp, _m in (g.get("subread_groups_by_hp") or {}).items():
            _pats.update(_m)
        _k = g.get("n_sSNV") or len(g.get("positions") or [])
        _nc, _comps = _components(_pats, _k)
        if _nc > 1:
            n_broken += 1

        by_region[rid] = {
            # 共現連通性：nc>1 表示這個 block 的 sSNV 沒有被 read 串成一塊
            "nComp": _nc,
            "comps": _comps if _nc > 1 else None,
            # 完整位點清單（含非 active）。pattern 字串的字元數對應這個，
            # 不是 topology 的 active_positions —— 兩者長度可能不同。
            "pos": g.get("positions") or [],
            "nSnv": g.get("n_sSNV"),
            "cov": g.get("col_coverage_by_hp") or {},
            "full": g.get("populations_by_hp") or {},
            "part": g.get("subread_groups_by_hp") or {},
            "nFull": g.get("n_full_cov_reads"),
            "readsByHp": g.get("reads_by_hp") or {},
            "projRows": g.get("projected_molecule_rows"),
            "treeSupported": g.get("tree_supported_molecule_block_incidences"),
            "patternCount": g.get("tree_supported_pattern_count"),
            "minRead": g.get("min_read"),
            "cn": g.get("cn"),
            "cnAvail": g.get("cn_data_available"),
            "vafEligible": g.get("vaf_eligible"),
        }

    hit = len(region_ids & set(by_region)) if region_ids else len(by_region)
    if not probe_linkage(cap, hit, len(region_ids) if region_ids else len(by_region),
                         "region_id 對得上 topology", min_rate=0.50,
                         fail_state=MALFORMED):
        return cap
    cap.counts["groups"] = len(by_region)
    cap.counts["broken_cooccurrence"] = n_broken
    probe(cap, "S2", True,
          f"投影後 pattern 表看不到共現的 block：{n_broken:,} / {len(by_region):,}"
          f"（{n_broken / max(len(by_region), 1) * 100:.2f}%）—— "
          f"這是 O/D/S/L→X 投影造成的視野限制，"
          f"**不是**缺乏 read 支持（連通性以 strict endpoint edges 為準，該處 100% 連通）")
    cap.payload = {"by_region": by_region}
    return finalize(cap)
