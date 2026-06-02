#!/usr/bin/env python3
"""
Per-read 甲基萃取 + HP tag → Read×CpG 矩陣（唯讀，不改 BAM）。

用途：給定 region，抽出所有 overlapping read 的 per-CpG 5mC 機率 + HP tag，
產出 (1) read×CpG 矩陣 TSV (2) per-read metadata TSV (HP/MAPQ/strand/n_CpG)。
下游由 heatmap.py 畫 clustermap（HP 側欄）。

設計原則（對齊 design_02 §6+R4/R8）：
- 甲基值用 ML qual (0-255) → 機率 ml/255；缺值 = NaN（不是 0，避免 0/未測混淆）。
- HP tag: HP:Z:{1,2,1-1,2-1,3}；無 tag → "unphase"。
- 只取 primary (-F 0x900 等效：skip secondary+supplementary)。
- CpG 參考座標：用 read 對 reference 的對位 (get_aligned_pairs) 把 query mod pos → ref pos。
- 限 5mC（mod code 'm' / 27551）。
"""
import sys, argparse, csv
import pysam
import numpy as np

def parse_region(s):
    chrom, rest = s.split(":")
    a, b = rest.replace(",", "").split("-")
    return chrom, int(a), int(b)

def get_hp(read):
    try:
        v = read.get_tag("HP")
        return str(v)
    except KeyError:
        return "unphase"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--region", required=True, help="chr:start-end (1-based, comma ok)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-cpg", type=int, default=3, help="read 至少 N 個 CpG 才納入")
    ap.add_argument("--max-reads", type=int, default=4000, help="安全上限")
    args = ap.parse_args()

    chrom, start, end = parse_region(args.region)
    bam = pysam.AlignmentFile(args.bam, "rb")

    # 收集所有 read 的 (ref_pos -> prob)，再對齊成共同 CpG 欄
    seen_ids = set()
    read_rows = []      # list of dict: read_id, hp, mapq, strand, methyl={ref_pos:prob}
    all_positions = set()
    n_seen = 0
    n_no_mod = 0

    for read in bam.fetch(chrom, start, end):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        n_seen += 1
        if n_seen > args.max_reads:
            break
        # pysam modified_bases: {(base, strand, mod): [(qpos, qual), ...]}
        mods = read.modified_bases or {}
        # 找 5mC key（base 'C', mod 'm' 或 27551）
        mc_calls = None
        for key, calls in mods.items():
            base, strand_flag, mod = key[0], key[1], key[2]
            if base in ("C", b"C") and mod in ("m", 27551):
                mc_calls = calls
                break
        if not mc_calls:
            n_no_mod += 1
            continue
        # query pos -> ref pos 對位表
        q2r = {}
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            q2r[qpos] = rpos
        methyl = {}
        for qpos, qual in mc_calls:
            rpos = q2r.get(qpos)
            if rpos is None:
                continue
            if rpos < start or rpos > end:
                continue
            methyl[rpos] = qual / 255.0
        if len(methyl) < args.min_cpg:
            continue
        if read.query_name in seen_ids:
            continue
        seen_ids.add(read.query_name)
        read_rows.append({
            "read_id": read.query_name,
            "hp": get_hp(read),
            "mapq": read.mapping_quality,
            "strand": "-" if read.is_reverse else "+",
            "methyl": methyl,
        })
        all_positions.update(methyl.keys())

    bam.close()

    positions = sorted(all_positions)
    if not read_rows:
        sys.stderr.write(f"[WARN] no reads with >={args.min_cpg} CpG in {args.region}\n")
        # 仍寫空檔以利下游判斷
    # 寫矩陣
    mat_path = args.out_prefix + "_matrix.tsv"
    meta_path = args.out_prefix + "_meta.tsv"
    with open(mat_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["read_id"] + [f"{chrom}:{p}" for p in positions])
        for r in read_rows:
            row = [r["read_id"]]
            for p in positions:
                v = r["methyl"].get(p)
                row.append("" if v is None else f"{v:.4f}")
            w.writerow(row)
    with open(meta_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["read_id", "hp", "mapq", "strand", "n_cpg"])
        for r in read_rows:
            w.writerow([r["read_id"], r["hp"], r["mapq"], r["strand"], len(r["methyl"])])

    # 摘要
    from collections import Counter
    hp_counts = Counter(r["hp"] for r in read_rows)
    print(f"region={args.region}")
    print(f"  reads_seen(primary)={min(n_seen,args.max_reads)}  no_mod_tag={n_no_mod}  reads_kept(>={args.min_cpg}CpG)={len(read_rows)}")
    print(f"  CpG positions={len(positions)}")
    print(f"  HP distribution: {dict(hp_counts)}")
    print(f"  -> {mat_path}")
    print(f"  -> {meta_path}")

if __name__ == "__main__":
    main()
