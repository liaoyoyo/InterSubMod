#!/usr/bin/env python3
"""per-region read-driven (HCC1395) — 把 read-driven 串接歸到 pipeline 區域。
逐區: 多-ALT read 數 / distinct combos / max chain → 與 topology n_clusters 交叉確認 + 建逐區 HTML。
argv[1]=chroms(逗號) argv[2]=out。§13: 純分析→落 JSON。"""
import gzip
import json
import sys
import os
from collections import defaultdict
from bisect import bisect_right
import pysam

TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
VD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
TOPO = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json"
MAPQ_MIN = 20
CHROMS = sys.argv[1].split(",") if len(sys.argv) > 1 else [f"chr{c}" for c in range(1, 23)]
OUT = sys.argv[2] if len(sys.argv) > 2 else "/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/e33ba939-69f7-4671-a8fd-ecdf6e8ebe55/scratchpad/rd_perregion.json"

topo = json.load(open(TOPO, encoding="utf-8"))
regions_by_chrom = defaultdict(list)
for r in topo["detail"]:
    regions_by_chrom[r["chrom"]].append((r["start"], r["start"] + r["span"], r["region"]))
for c in regions_by_chrom:
    regions_by_chrom[c].sort()


def load_union(chrom):
    snvs = {}
    for src in ("tp", "fp"):
        f = f"{VD}/filtered_snv_{src}_{chrom}.vcf.gz"
        try:
            with gzip.open(f, "rt") as fh:
                for ln in fh:
                    if ln.startswith("#"):
                        continue
                    p = ln.split("\t")
                    pos = int(p[1])
                    if pos not in snvs:
                        snvs[pos] = (p[3].upper(), p[4].strip().upper())
        except FileNotFoundError:
            pass
    return snvs


def region_of(pos, starts, regs):
    i = bisect_right(starts, pos) - 1
    if i >= 0 and regs[i][0] <= pos <= regs[i][1]:
        return regs[i][2]
    return None


bam = pysam.AlignmentFile(TBAM, "rb")
per_region = defaultdict(lambda: {"multi_alt": 0, "max_chain": 0, "combos": set()})
for chrom in CHROMS:
    snvs = load_union(chrom)
    if not snvs:
        continue
    regs = regions_by_chrom.get(chrom, [])
    starts = [r[0] for r in regs]
    positions = sorted(snvs)
    pos2reg = {p: region_of(p, starts, regs) for p in positions}
    read_vec = defaultdict(dict)
    for pos in positions:
        ref, alt, = snvs[pos][0], snvs[pos][1]
        for col in bam.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=0, stepper="samtools"):
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                aln = pr.alignment
                if aln.mapping_quality < MAPQ_MIN:
                    continue
                base = aln.query_sequence[pr.query_position].upper()
                read_vec[aln.query_name][pos] = "ALT" if base == alt else ("REF" if base == ref else "O")
    for rn, vec in read_vec.items():
        alts = [p for p, v in vec.items() if v == "ALT"]
        if len(alts) < 2:
            continue
        byreg = defaultdict(list)
        for p in alts:
            rg = pos2reg.get(p)
            if rg is not None:
                byreg[rg].append(p)
        for rg, plist in byreg.items():
            if len(plist) >= 2:
                d = per_region[rg]
                d["multi_alt"] += 1
                d["max_chain"] = max(d["max_chain"], len(plist))
                d["combos"].add(frozenset(plist))
    del read_vec
    print(f"[{chrom}] regions={len(regs)} regions_with_chains={sum(1 for rg in per_region if rg.startswith(chrom+':'))}", flush=True)
bam.close()
out = {rg: {"rd_multi_alt": d["multi_alt"], "rd_max_chain": d["max_chain"], "rd_combos": len(d["combos"])}
       for rg, d in per_region.items()}
json.dump(out, open(OUT, "w"), ensure_ascii=False)
print(f"DONE regions_with_chains={len(out)} -> {OUT}")
