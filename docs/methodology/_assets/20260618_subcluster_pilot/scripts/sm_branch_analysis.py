"""[分支型態 + 相鄰區域一致性] 讀 sm_region_integration.json。
1. 分支拓樸 taxonomy: 每個有結構區域 → 標準型態(linear/bifurcation/sibling-pair/full-tree...) + modal。
2. 相鄰區域 shape 一致性: 同 chrom 連續區域是否同 shape；observed vs shuffle-null（空間自相關）。
3. VAF-level 一致性: 全基因組區域主 VAF 是否有一致峰（= subclone CCF 跨區重現）。
輸出 sm_branch_analysis.json。np.random seed 固定可重算。
"""
import json
import numpy as np
from collections import Counter, defaultdict

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
np.random.seed(7)
reg = json.load(open(f"{A}/sm_region_integration.json"))["regions"]


def topo_label(r):
    """標準分支型態。"""
    nn, nest, sib, dep = r["n_nodes"], r["n_nested"], r["n_sibling"], r["max_depth"]
    if r["has_cycle"]:
        return "inconsistent"
    if nest == 0 and sib == 0:
        return "co_linked_single_lineage" if r.get("n_colinked_merges", 0) >= 1 else "no_structure"
    if nest >= 1 and sib >= 1:
        return "full_tree(branch+depth)"
    if sib >= 1 and nest == 0:
        return "sibling_only(2branch)" if sib == 1 else "sibling_multi(>=3branch)"
    # nested only
    if dep >= 2:
        return "linear_chain(depth>=2)"
    # depth 1 nested: 1 ancestor -> k descendants
    # out-degree of root
    outdeg = Counter()
    for u, v in r["nested_edges"]:
        outdeg[u] += 1
    mo = max(outdeg.values()) if outdeg else 1
    return "single_nested(A->B)" if mo == 1 else f"fan_out(A->{mo})"


# --- 1. taxonomy ---
struct = [r for r in reg if r["tree_shape"] not in ("no_confirmed_structure", "inconsistent")]
labels = Counter(topo_label(r) for r in struct)
ns = len(struct)
# branch-bearing only (有實際分支: sibling>=1 或 fan_out 或 full_tree)
branch = [r for r in struct if r["n_sibling"] >= 1 or any(c >= 2 for c in [Counter(u for u, v in r["nested_edges"])[k] for k in set(u for u, v in r["nested_edges"])] or [0])]

# --- 2. 相鄰 shape 一致性 ---
bychrom = defaultdict(list)
for r in reg:
    bychrom[r["chrom"]].append(r)
adj_pairs = 0
adj_same = 0
adj_dist = []
shapes_seq = []
for c, rs in bychrom.items():
    rs = sorted(rs, key=lambda x: x["start"])
    for i in range(len(rs) - 1):
        adj_pairs += 1
        if rs[i]["tree_shape"] == rs[i + 1]["tree_shape"]:
            adj_same += 1
        adj_dist.append(rs[i + 1]["start"] - rs[i]["end"])
obs_rate = adj_same / adj_pairs
# null: shuffle shapes within chrom
null_rates = []
for _ in range(200):
    s = 0
    p = 0
    for c, rs in bychrom.items():
        sh = [r["tree_shape"] for r in rs]
        np.random.shuffle(sh)
        for i in range(len(sh) - 1):
            p += 1
            if sh[i] == sh[i + 1]:
                s += 1
    null_rates.append(s / p)
null_mean = float(np.mean(null_rates))
null_sd = float(np.std(null_rates))
z = (obs_rate - null_mean) / null_sd if null_sd > 0 else 0

# --- 3. VAF-level 一致性（全基因組區域主 VAF 峰）---
# 每區域取各 node VAF；CN-clean(loh/neutral) 才用
clean_vafs = []
for r in reg:
    if r["cn"] in ("loh", "neutral") and r.get("vaf"):
        clean_vafs.extend(r["vaf"].values())
hist, edges = np.histogram(clean_vafs, bins=20, range=(0, 1))
peaks = sorted([(round(edges[i], 2), int(hist[i])) for i in range(20)], key=lambda x: -x[1])[:5]

out = {
    "n_structured_regions": ns,
    "branch_topology_taxonomy": dict(labels.most_common()),
    "branch_topology_proportion": {k: round(v / ns, 3) for k, v in labels.most_common()},
    "modal_pattern": labels.most_common(1)[0],
    "adjacent_consistency": {
        "n_adjacent_pairs": adj_pairs,
        "observed_same_shape_rate": round(obs_rate, 3),
        "null_mean": round(null_mean, 3), "null_sd": round(null_sd, 4), "z_score": round(z, 1),
        "median_adjacent_gap_bp": int(np.median(adj_dist)),
        "interpretation": "obs >> null → 相鄰區域結構正相關(空間相干); obs≈null → 各區域獨立(無空間結構)",
    },
    "vaf_consistency_clean": {
        "n_clean_vaf_values": len(clean_vafs),
        "top5_vaf_peaks(bin_start,count)": peaks,
        "note": "若集中在少數峰 = subclone CCF 跨區一致; 若均勻 = 無一致 subclone 層級",
    },
}
json.dump(out, open(f"{A}/sm_branch_analysis.json", "w"), ensure_ascii=False, indent=1)
print("=== 1. 分支型態 taxonomy (structured n=%d) ===" % ns)
for k, v in labels.most_common():
    print(f"  {k:30s} {v:>5} ({100*v/ns:5.1f}%)")
print(f"\n=== 2. 相鄰區域 shape 一致性 ===")
print(f"  adjacent pairs={adj_pairs}  observed same-shape={obs_rate:.3f}  null={null_mean:.3f}±{null_sd:.4f}  z={z:.1f}")
print(f"  median adjacent gap={int(np.median(adj_dist)):,}bp")
print(f"\n=== 3. CN-clean VAF 峰 (subclone CCF 一致性) ===")
print(f"  n_clean_vaf={len(clean_vafs)}  top peaks(VAF:count): {peaks}")
print(f"-> sm_branch_analysis.json")
