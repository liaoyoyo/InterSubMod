#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Per-sample read-driven 新版數據 deriver（零 BAM 重讀）。
從現成 sm_linkage_genomewide.json 的 census key 衍生 isolated/underpowered/linked
（規則同 HCC1395 build_ideogram_data.py：isolated=n_partners_le50k==0；
 underpowered=partner>0 & realized==0；linked=realized>0），
+ topology_per_region.json 樹 markers → 產 census_categories.json + ideogram_data.json。
用法: build_sample_ideogram.py <sample_dir>   (寫回該 dir)
§13: 純衍生，全部數字來自磁碟現成檔，無捏造。"""
import json
import os
import sys

SDIR = sys.argv[1]
LINK = os.path.join(SDIR, "sm_linkage_genomewide.json")
TOPO = os.path.join(SDIR, "topology_per_region.json")
OUT_CAT = os.path.join(SDIR, "census_categories.json")
OUT_IDE = os.path.join(SDIR, "ideogram_data.json")
BIN = 2_000_000

HG38 = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
        "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
        "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
        "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468}
SHAPE = {"full_tree": "F", "linear_nested": "S", "sibling_only": "S", "co_linked_lineage": "S",
         "inconsistent": "I", "no_confirmed_structure": "N"}

link = json.load(open(LINK, encoding="utf-8"))
census = link["census"]

iso, und, lnk = [], [], []
for k, v in census.items():
    npar = v.get("n_partners_le50k", 0) or 0
    nrl = v.get("n_realized_links", 0) or 0
    if npar == 0:
        iso.append(k)
    elif nrl == 0:
        und.append(k)
    else:
        lnk.append(k)
json.dump({"isolated": iso, "underpowered": und, "linked": lnk},
          open(OUT_CAT, "w"), ensure_ascii=False)

topo = json.load(open(TOPO, encoding="utf-8"))
trees_by_chrom = {}
for r in topo["detail"]:
    trees_by_chrom.setdefault(r["chrom"], []).append([r["start"], SHAPE.get(r.get("tree_shape"), "N")])
for c in trees_by_chrom:
    trees_by_chrom[c].sort()

# --- artifact 分類（>50kb 跨距 / >8 sSNV 窗外塌縮，從 span 衍生，零 BAM）---
# 依據：截斷區(>8 sSNV)的 tree 由窗內 multilocus(cap 8)建；span>100kb 時無單條 ONT read 能端到端
# 跨越該區(read N50 10-30kb、ultralong 罕>100kb)→ 「樹」無法被 read 端到端實現＝窗口聚合 artifact，
# 結構不可信。span<=100kb 但 >8 sSNV＝dense_capped：read 可能真跨→rd_perregion 才有定樹增益(候選)。
SPAN_READ_BOUND = 100_000
artifacts = []
for r in topo["detail"]:
    if not r.get("truncated"):
        continue
    cls = "window_aggregation" if r["span"] > SPAN_READ_BOUND else "dense_capped"
    artifacts.append({"region": r["region"], "chrom": r["chrom"], "start": r["start"],
                      "span": r["span"], "n_sSNV": r["n_sSNV"], "cn": r.get("cn"),
                      "klass": cls, "has_cycle": bool(r.get("has_cycle")),
                      "cycle_cause": r.get("cycle_cause"), "tree_shape": r.get("tree_shape")})
art_summary = {"n_truncated": len(artifacts),
               "n_window_aggregation": sum(1 for a in artifacts if a["klass"] == "window_aggregation"),
               "n_dense_capped": sum(1 for a in artifacts if a["klass"] == "dense_capped"),
               "n_cycle": sum(1 for r in topo["detail"] if r.get("has_cycle")),
               "span_read_bound": SPAN_READ_BOUND}


def bin_positions(keys, chrom, nbin):
    b = [0] * nbin
    for kk in keys:
        c, p = kk.rsplit(":", 1)
        if c == chrom:
            idx = int(p) // BIN
            if 0 <= idx < nbin:
                b[idx] += 1
    return b


per_chrom = {}
tot = {"isolated": len(iso), "underpowered": len(und), "linked": len(lnk)}
for chrom in HG38:
    nbin = HG38[chrom] // BIN + 1
    iso_b = bin_positions(iso, chrom, nbin)
    und_b = bin_positions(und, chrom, nbin)
    per_chrom[chrom] = {"len": HG38[chrom], "bin": BIN,
                        "trees": trees_by_chrom.get(chrom, []),
                        "iso_bins": iso_b, "und_bins": und_b,
                        "n_tree": len(trees_by_chrom.get(chrom, [])),
                        "n_iso": sum(iso_b), "n_und": sum(und_b)}

json.dump({"source": "sm_linkage census + topology_per_region (GRCh38, derived no-BAM)",
           "bin": BIN, "totals": tot, "per_chrom": per_chrom,
           "artifacts_summary": art_summary, "artifacts": artifacts},
          open(OUT_IDE, "w"), ensure_ascii=False)
print(f"{os.path.basename(SDIR)}: census_total={len(census)} "
      f"isolated={tot['isolated']} underpowered={tot['underpowered']} linked={tot['linked']} "
      f"sum={sum(tot.values())} | tree_regions={len(topo['detail'])} | "
      f"artifact: trunc={art_summary['n_truncated']} "
      f"window_agg={art_summary['n_window_aggregation']} dense_capped={art_summary['n_dense_capped']} "
      f"cycle={art_summary['n_cycle']}")
