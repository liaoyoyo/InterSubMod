#!/usr/bin/env python3
"""硬核心：topology.jsonl。缺這個就沒有骨幹可畫，直接 refuse。

產出兩層：
    L1  全基因組 sSNV 索引 —— 每個 (chrom, pos) 一筆，columnar 送給前端
    L2  region 層 —— 代表樹的 vertices/edges、逐位點 AF、各種 status

一個 sSNV 可能落在多個 region（實測 19,849 個 distinct pos 對 20,060 個
(region,pos) pair），所以 L1→L2 用 CSR 表示多對多，前端點到多 region 的位點時
列出全部讓使用者自己選，不替他挑。
"""
from __future__ import annotations

import json
from pathlib import Path

from capability import (CORE, MALFORMED, Capability, Registry, file_ref, finalize,
                        first_jsonl, make, open_maybe_gzip, probe, probe_count,
                        probe_exists, probe_schema)

SCHEMA = "intersubmod.exact_ps_cpp_topology_af.unit"

# 每一筆都必須有的欄位。探測時檢查首筆。
REQUIRED_FIELDS = {
    "region_id", "chrom", "active_positions", "hp_family", "phase_set",
    "unit_id", "block_id", "active_bit_count", "unit_status",
}

# 只有「跑出樹」的 region 才有這些欄位 —— unit_status 為 zero_denominator /
# no_active_alt / family_incomplete 的 region 沒有 representative tree。
# 因此它們不能列進 REQUIRED_FIELDS（曾誤列，導致首筆恰為無樹 region 的樣本
# 被判成 MALFORMED），改在掃完後統計有樹比例。
TREE_FIELDS = {"representative_best_vertices", "representative_best_edges"}

# 22 條體染色體。chrX/Y 不在拓撲 scope 內（census 定 chr1-22），
# 若 jsonl 出現其他 contig 會被記進 skipped 而不是靜默丟掉。
CHROMS = [f"chr{i}" for i in range(1, 23)]
CHROM_INDEX = {c: i for i, c in enumerate(CHROMS)}

CHROM_LEN = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468,
}


def _int(v, default=0):
    try:
        return int(v)
    except (TypeError, ValueError):
        return default


def _compact_region(rec: dict, idx: int) -> dict:
    """把 topology 的 44 欄壓成前端要的形狀。欄名縮短是因為這會乘上 11,590 筆，
    長欄名在 JSON 裡是純粹的浪費。"""
    verts = rec.get("representative_best_vertices") or []
    edges = rec.get("representative_best_edges") or []
    af = rec.get("af_coverage") or []
    return {
        "i": idx,
        "id": rec.get("region_id", ""),
        "c": rec.get("chrom", ""),
        "ps": str(rec.get("phase_set", "")),
        "hp": str(rec.get("hp_family", "")),
        "u": rec.get("unit_id", ""),
        "b": rec.get("block_id", ""),
        "k": _int(rec.get("active_bit_count")),
        "k0": _int(rec.get("original_bit_count")),
        "ap": list(rec.get("active_positions") or []),
        # vertex：[編號, label]。label 以 H_ 起頭者是 solver 補入的隱藏節點
        "v": [[_int(v.get("vertex")), str(v.get("label", ""))] for v in verts],
        # edge：[父, 子, 獲得位置, edge score]
        "e": [[_int(e.get("parent_vertex")), _int(e.get("child_vertex")),
               _int(e.get("acquired_position")), str(e.get("edge_score_fraction", ""))]
              for e in edges],
        "af": [[_int(a.get("position")), _int(a.get("alt_reads")), _int(a.get("ref_reads"))]
               for a in af],
        "tu": bool(rec.get("best_tree_unique")),
        "tc": str(rec.get("best_tree_tie_count", "")),
        "tt": str(rec.get("total_tree_count", "")),
        "vs": _int(rec.get("best_vertex_set_count")),
        "ms": _int(rec.get("minimum_vertex_set_count")),
        "us": rec.get("unit_status", ""),
        "fs": rec.get("family_status", ""),
        "os": rec.get("objective_status", ""),
        "rs": rec.get("read_af_status", ""),
        "rr": bool(rec.get("recurrence_required")),
        "rc": _int(rec.get("recurrence_compatible_vertex_set_count")),
        "ri": _int(rec.get("recurrence_incompatible_vertex_set_count")),
        "sn": _int(rec.get("search_nodes")),
        "reason": rec.get("reason", ""),
    }


def load(reg: Registry, topology_path, receipt_path=None, expected_sample=None) -> Capability:
    cap = make(reg, "topology", "拓撲骨幹（topology.jsonl）", tier=CORE,
               enables=["genome", "region", "tree", "filter:k", "filter:hp",
                        "filter:unit_status", "cooccur:k"])

    if not probe_exists(cap, topology_path):
        return cap

    head = first_jsonl(topology_path)
    if not probe_schema(cap, head, SCHEMA):
        return cap
    if expected_sample:
        observed = str(head.get("sample") or "").strip()
        if observed != str(expected_sample):
            cap.state = MALFORMED
            probe(cap, "S1", False,
                  f"topology sample 不符：--sample={expected_sample}，首筆 sample={observed or '(缺)'}")
            return cap
        probe(cap, "S1", True, f"sample 身分一致：{observed}")
    missing = sorted(REQUIRED_FIELDS - set(head))
    if missing:
        cap.state = "MALFORMED"
        probe(cap, "S1", False, f"首筆缺必要欄位：{', '.join(missing)}")
        return cap

    regions, skipped = [], {}
    with_tree = 0
    # pos_map: (chrom, pos) -> [region 索引…]，之後轉成 CSR
    pos_map = {}
    with open_maybe_gzip(topology_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if expected_sample and str(rec.get("sample") or "").strip() != str(expected_sample):
                cap.state = MALFORMED
                probe(cap, "S1", False,
                      f"topology 內混入異樣本：region={rec.get('region_id', '?')}，"
                      f"sample={rec.get('sample') or '(缺)'}，期望 {expected_sample}")
                return cap
            chrom = rec.get("chrom", "")
            if chrom not in CHROM_INDEX:
                skipped[chrom] = skipped.get(chrom, 0) + 1
                continue
            idx = len(regions)
            row = _compact_region(rec, idx)
            if row["v"]:
                with_tree += 1
            regions.append(row)
            for p in row["ap"]:
                pos_map.setdefault((chrom, p), []).append(idx)

    if skipped:
        probe(cap, "S1", True,
              "略過非 chr1-22 的 contig：" + ", ".join(f"{k}×{v}" for k, v in sorted(skipped.items())))

    expected = None
    if receipt_path:
        rp = Path(receipt_path)
        if rp.is_file():
            cap.paths.append(file_ref(rp))
            try:
                rec = json.loads(rp.read_text(encoding="utf-8"))
                for key in ("units", "unit_count", "records", "rows"):
                    if isinstance(rec.get(key), int):
                        expected = rec[key]
                        break
            except (OSError, ValueError):
                expected = None
    probe_count(cap, len(regions), expected, "region")

    cap.counts["region_with_tree"] = with_tree
    cap.counts["region_without_tree"] = len(regions) - with_tree
    no_tree_status = {}
    for r in regions:
        if not r["v"]:
            no_tree_status[r["us"]] = no_tree_status.get(r["us"], 0) + 1
    cap.counts["no_tree_by_status"] = no_tree_status
    probe(cap, "S2", True,
          f"有代表樹的 region {with_tree:,} / {len(regions):,}"
          f"（無樹者：{', '.join(f'{k}×{v:,}' for k, v in sorted(no_tree_status.items())) or '無'}）")

    # L1：sSNV 索引。排序讓前端可以直接二分搜尋，不用自己再排一次。
    keys = sorted(pos_map.keys(), key=lambda t: (CHROM_INDEX[t[0]], t[1]))
    l1 = _build_l1(keys, pos_map, regions)

    cap.counts["sSNV"] = len(keys)
    cap.counts["region_pos_pair"] = sum(len(v) for v in pos_map.values())
    cap.counts["multi_region_sSNV"] = sum(1 for v in pos_map.values() if len(v) > 1)
    probe(cap, "S3", True,
          f"sSNV {len(keys):,}；(region,pos) pair {cap.counts['region_pos_pair']:,}；"
          f"跨多 region 的 sSNV {cap.counts['multi_region_sSNV']:,}")

    cap.payload = {"regions": regions, "l1": l1, "pos_index": {k: i for i, k in enumerate(keys)}}
    return finalize(cap)


def _build_l1(keys, pos_map, regions) -> dict:
    """columnar 版 sSNV 索引。位置用染色體內 delta 存，因為排序後相鄰位置很近，
    delta 的數值範圍小很多，JSON 體積差幾倍。"""
    chrom_col, dpos, k_col, flags = [], [], [], []
    csr_off, csr_val = [0], []
    prev_chrom, prev_pos = None, 0

    for (chrom, pos) in keys:
        ci = CHROM_INDEX[chrom]
        chrom_col.append(ci)
        if chrom != prev_chrom:
            dpos.append(pos)
            prev_chrom, prev_pos = chrom, pos
        else:
            dpos.append(pos - prev_pos)
            prev_pos = pos

        idxs = pos_map[(chrom, pos)]
        csr_val.extend(idxs)
        csr_off.append(len(csr_val))

        # 這個 sSNV 的 k 取所屬 region 的最大值（同一位點在不同 region 可能 k 不同）
        k_col.append(max(regions[i]["k"] for i in idxs))

        f = 0
        if any(regions[i]["tu"] for i in idxs):
            f |= 1 << 0          # 至少一個 region 是唯一解
        if all(regions[i]["tu"] for i in idxs):
            f |= 1 << 1          # 全部 region 都是唯一解
        if any(regions[i]["us"] == "ranked" for i in idxs):
            f |= 1 << 2
        if len(idxs) > 1:
            f |= 1 << 3          # 跨多 region
        if any(regions[i]["rr"] for i in idxs):
            f |= 1 << 4          # 需要 recurrence
        flags.append(f)

    return {
        "n": len(keys),
        "chrom": chrom_col,
        "dpos": dpos,
        "k": k_col,
        "flags": flags,
        "csrOff": csr_off,
        "csrVal": csr_val,
        "chroms": CHROMS,
        "chromLen": [CHROM_LEN[c] for c in CHROMS],
        "flagBits": ["anyUnique", "allUnique", "anyRanked", "multiRegion", "recurrence"],
        # 甲基軸碼：由 ISM 能力事後填入（見 build_drilldown_dashboard.py）。
        # 0 = 無 ISM 資料；bit0 HP 顯著 / bit1 ALT 顯著 / bit2 lineage 顯著 /
        # bit3 有測但都不顯著。刻意不含「甲基自身分群」—— 那個軸是循環論證。
        "axis": [],
    }
