"""unit_lineage_paths.tsv.gz —— LongLineage 解出的每個樹頂點一列。

為什麼需要這一層：自檢 C8 要交叉驗證「topology 裡 `H_` 前綴的節點數」
與「lineage_paths 裡 is_hidden=true 的節點數」是否一致。兩者是同一件事
（solver 為使樹合法而補入、reads 未直接觀察到的中間祖先）由兩支不同程式
各自寫出，對得上才代表兩邊看的是同一棵樹。

⚠ 這支模組是補接線用的。selfcheck.py 從一開始就 `reg.get("lineage_paths")`，
但沒有任何 source 註冊過這個能力 —— C8 因此永遠停在「無法檢查」，
看起來像資料不全，其實是自己的接線缺口。2026-08-09 補上。
"""

from __future__ import annotations

import gzip
from pathlib import Path

from capability import EXTENSION, Capability, Registry, finalize, make, probe, probe_linkage

REQUIRED = {"region_id", "vertex", "vertex_label", "is_hidden", "depth", "lineage_path"}


def load(reg: Registry, paths_dir, region_ids=None) -> Capability:
    cap = make(reg, "lineage_paths", "lineage 頂點表（unit_lineage_paths）",
               tier=EXTENSION, enables=["selfcheck:C8", "panel:tree_depth"])
    d = Path(paths_dir) if paths_dir else None
    if not d or not d.is_dir():
        cap.state = "ABSENT"
        probe(cap, "S0", False, f"paths 目錄不存在：{d}")
        return cap
    files = sorted(d.glob("*.unit_lineage_paths.tsv.gz"))
    if not files:
        cap.state = "ABSENT"
        probe(cap, "S0", False, f"{d} 下沒有 *.unit_lineage_paths.tsv.gz")
        return cap
    total_bytes = sum(f.stat().st_size for f in files)
    probe(cap, "S0", True, f"{len(files)} 檔、{total_bytes:,} bytes")

    hidden = 0
    rows = 0
    regions: set[str] = set()
    depth_max = 0
    try:
        for f in files:
            with gzip.open(f, "rt", encoding="utf-8") as fh:
                header = fh.readline().rstrip("\n").split("\t")
                missing = REQUIRED - set(header)
                if missing:
                    cap.state = "MALFORMED"
                    probe(cap, "S1", False, f"{f.name} 缺欄位：{', '.join(sorted(missing))}")
                    return cap
                idx = {c: i for i, c in enumerate(header)}
                i_rid, i_hid, i_dep = idx["region_id"], idx["is_hidden"], idx["depth"]
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) <= i_dep:
                        continue
                    rows += 1
                    regions.add(p[i_rid])
                    # 欄位是字串 "true"/"false" —— 直接當布林用會全部成真
                    if p[i_hid].strip().lower() == "true":
                        hidden += 1
                    try:
                        depth_max = max(depth_max, int(p[i_dep]))
                    except ValueError:
                        pass
    except (OSError, ValueError) as exc:
        cap.state = "MALFORMED"
        probe(cap, "S1", False, f"讀取失敗：{exc}")
        return cap

    probe(cap, "S1", True, f"必要欄位齊備（{len(REQUIRED)} 欄）")
    probe(cap, "S2", True,
          f"頂點 {rows:,} 個、region {len(regions):,} 個、"
          f"補入節點 {hidden:,}（{hidden / max(rows, 1) * 100:.1f}%）、最深 depth {depth_max}")

    if region_ids:
        probe_linkage(cap, len(regions & set(region_ids)), len(region_ids),
                      "region_id 對得上 topology")

    cap.counts["vertices"] = rows
    cap.counts["hidden"] = hidden
    cap.payload = {"hidden": hidden, "vertices": rows,
                   "regions": len(regions), "depth_max": depth_max}
    return finalize(cap)
