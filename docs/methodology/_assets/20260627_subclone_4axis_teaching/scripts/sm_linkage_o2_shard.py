#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""sm_linkage O2 分片 worker(2026-07-01):
處理一個位置分片 [lo,hi](邊界落在 gap>50kb 的群間隙 → 不切斷任何 linkage group → byte-identical)。
O2: per_read_calls 改 per-group 單次 region pileup(已驗 397 group + chr22 端到端 7/7 byte-identical、chr8 超寬群 2.98-6.19x)。
byte-identical 保證: main() 完全不動,只(a)monkeypatch per_read_calls=region 版(輸出等價) (b)load_union 過濾到 [lo,hi]。
argv: chrom lo hi out_path。env: SM_TBAM/SM_NBAM/SM_VD/SM_VCF_MODE。記 wall+peak RSS 到 out.bench。"""
import os, sys, glob, time, json, resource
from collections import defaultdict
SCD = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCD)
import sm_linkage_genomewide as M
import pysam

def per_read_calls_o2(bam, chrom, group):
    """O2: group span 單次 region pileup,只在 sSNV 位點 tabulate。與逐位點版 byte-identical。"""
    calls = defaultdict(dict); hp = {}; ps = {}
    refalt = {pos: (ref, alt) for pos, ref, alt, _ in group}
    posset = set(refalt); lo = min(posset); hi = max(posset)
    for col in bam.pileup(chrom, lo - 1, hi, truncate=True, min_base_quality=0, stepper="samtools"):
        pos = col.reference_pos + 1
        if pos not in posset: continue
        ref, alt = refalt[pos]
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None: continue
            aln = pr.alignment
            if aln.mapping_quality < M.MAPQ_MIN: continue
            rn = aln.query_name
            base = aln.query_sequence[pr.query_position].upper()
            calls[rn][pos] = "REF" if base == ref else ("ALT" if base == alt else "OTHER")
            if aln.has_tag("HP"): hp[rn] = aln.get_tag("HP")
            if aln.has_tag("PS"): ps[rn] = aln.get_tag("PS")
    return calls, hp, ps

def main():
    chrom, lo, hi, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    M.per_read_calls = per_read_calls_o2                       # O2 patch
    _orig = M.load_union
    M.load_union = lambda c: [s for s in _orig(c) if lo <= s[0] <= hi]   # 過濾到分片
    t0 = time.time()
    M.main([chrom], out)
    wall = time.time() - t0
    peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0   # Linux ru_maxrss = KB
    d = json.load(open(out))
    json.dump({"chrom": chrom, "lo": lo, "hi": hi, "wall_sec": round(wall, 1),
               "peak_rss_mb": round(peak_mb, 1), "n_pairs": d["n_pairs_recorded"],
               "n_groups_multi": d["n_groups_multi"], "n_singletons": d["n_singletons"]},
              open(out + ".bench", "w"))
    print(f"SHARD {chrom}:{lo}-{hi} wall={wall:.1f}s peak={peak_mb:.0f}MB pairs={d['n_pairs_recorded']}", flush=True)

if __name__ == "__main__":
    main()
