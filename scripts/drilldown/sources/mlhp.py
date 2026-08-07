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

    by_region = {}
    for g in groups:
        rid = g.get("region_id")
        if not rid:
            continue
        by_region[rid] = {
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
    cap.payload = {"by_region": by_region}
    return finalize(cap)
