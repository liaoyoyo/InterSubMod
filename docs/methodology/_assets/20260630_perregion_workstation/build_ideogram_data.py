#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HG38 ideogram 資料（census-based，精準）。
isolated/underpowered/linked 逐位點來自 census_categories.json（sm_linkage_genomewide census：
isolated=n_partners_le50k==0；underpowered=partner>0 & realized==0；linked=realized>0）。
tree 區 markers + shape 來自 topology_per_region.json。iso/under binned 2Mb 密度。
§13: 從 census + detail 算，計數須對齊 tier_r(8320/5458/21554)。"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CAT = "/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/e33ba939-69f7-4671-a8fd-ecdf6e8ebe55/scratchpad/census_categories.json"
TOPO = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json"
OUT = os.path.join(HERE, "data", "ideogram_data.json")
BIN = 2_000_000

HG38 = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
        "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
        "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
        "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468}

cat = json.load(open(CAT, encoding="utf-8"))   # {isolated:[...], underpowered:[...], linked:[...]}
topo = json.load(open(TOPO, encoding="utf-8"))
SHAPE = {"full_tree": "F", "linear_nested": "S", "sibling_only": "S", "co_linked_lineage": "S",
         "inconsistent": "I", "no_confirmed_structure": "N"}

trees_by_chrom = {}
for r in topo["detail"]:
    trees_by_chrom.setdefault(r["chrom"], []).append([r["start"], SHAPE.get(r.get("tree_shape"), "N")])
for c in trees_by_chrom:
    trees_by_chrom[c].sort()


def bin_positions(keys, chrom, nbin):
    b = [0] * nbin
    for k in keys:
        c, p = k.rsplit(":", 1)
        if c == chrom:
            idx = int(p) // BIN
            if 0 <= idx < nbin:
                b[idx] += 1
    return b


iso_keys = cat.get("isolated", [])
und_keys = cat.get("underpowered", [])
per_chrom = {}
tot = {"isolated": len(iso_keys), "underpowered": len(und_keys), "linked": len(cat.get("linked", []))}
for chrom in HG38:
    nbin = HG38[chrom] // BIN + 1
    iso_b = bin_positions(iso_keys, chrom, nbin)
    und_b = bin_positions(und_keys, chrom, nbin)
    per_chrom[chrom] = {"len": HG38[chrom], "bin": BIN,
                        "trees": trees_by_chrom.get(chrom, []),
                        "iso_bins": iso_b, "und_bins": und_b,
                        "n_tree": len(trees_by_chrom.get(chrom, [])),
                        "n_iso": sum(iso_b), "n_und": sum(und_b)}

json.dump({"source": "census_categories + topology_per_region (GRCh38)", "bin": BIN,
           "totals": tot, "per_chrom": per_chrom}, open(OUT, "w"), ensure_ascii=False)
print(f"DONE -> {OUT}  isolated={tot['isolated']} underpowered={tot['underpowered']} "
      f"linked={tot['linked']} (對照 tier_r 8320/5458/21554)")
