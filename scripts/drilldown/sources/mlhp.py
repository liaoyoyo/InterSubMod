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

from capability import Capability, Registry, finalize, make, probe, probe_exists, probe_linkage


def load(reg: Registry, path, region_ids: set) -> Capability:
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

        一個 block 理應是「被 read 串連起來的 sSNV 集合」。若共現圖斷裂，
        代表那些位點之間沒有任何 read 同時覆蓋 —— solver 仍會建出一棵樹，
        但跨分量的邊完全沒有 read 支持。實測 16.16% 的 block 有這個情形。
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
    probe_linkage(cap, hit, len(region_ids) if region_ids else len(by_region),
                  "region_id 對得上 topology")
    cap.counts["groups"] = len(by_region)
    cap.counts["broken_cooccurrence"] = n_broken
    probe(cap, "S2", True,
          f"共現圖斷裂的 block：{n_broken:,} / {len(by_region):,}"
          f"（{n_broken / max(len(by_region), 1) * 100:.2f}%）—— "
          f"那些 block 的 sSNV 之間沒有 read 串連，跨分量的樹邊零 read 支持")
    cap.payload = {"by_region": by_region}
    return finalize(cap)
