#!/usr/bin/env python3
"""
T2 嚴格驗證 + 邊緣分類：H3 read 的 germline 覆蓋狀況 + 甲基可指派性。

longphase-S H3 = 帶 somatic 變異、somatic 一致(≥0.6)、但 germline 不一致/缺(norHPsim<0.6)。
用戶點名兩類特殊處理：(a) 沒經過任何 germline (b) 經過不一致 germline。

每 H3 read 用 germline VCF 數覆蓋的 het 數 + 一致性 → 分類：
  no_germline (0 het) / inconsistent (≥2 het 但 H1/H2 票分歧 <0.8) / weak (1 het 或弱一致)
然後用 1-1 vs 2-1 的甲基質心（同窗）把 H3 read 指派 H1-1/H2-1 → 信心 margin（各類比較）。
⚠ H3 無 germline ground truth → 報「可指派性/信心分布」而非 accuracy；正控=同窗 1-1/2-1 held-out 自我判別。
唯讀。每 chr 落 rigor_t2_{chrom}.json。
"""
import sys, json, argparse, bisect
import numpy as np
import pysam

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GVCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
SOMATIC_VCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/somatic_pass.vcf.gz"

def get_hp(read):
    try: return str(read.get_tag("HP"))
    except KeyError: return "unphase"

def load_hets(chrom):
    vcf = pysam.VariantFile(GVCF); hl = {}
    for rec in vcf.fetch(chrom):
        s0 = rec.samples[0]; gt = s0.get("GT")
        if gt is None or len(gt) != 2 or gt[0] == gt[1] or None in gt: continue
        if len(rec.alleles) < 2 or len(rec.ref) != 1 or len(rec.alts[0]) != 1: continue
        hl[rec.pos - 1] = (rec.ref.upper(), rec.alts[0].upper(), gt)
    return hl

def read_meth(read, ws, we, min_cpg=5):
    mods = read.modified_bases or {}; mc = None
    for k, calls in mods.items():
        if k[0] in ("C", b"C") and k[2] in ("m", 27551): mc = calls; break
    if not mc: return None
    q2r = {q: r for q, r in read.get_aligned_pairs(matches_only=True)}
    meth = {}
    for qpos, qual in mc:
        rpos = q2r.get(qpos)
        if rpos is None or rpos < ws or rpos > we: continue
        meth[rpos] = qual / 255.0
    return meth if len(meth) >= min_cpg else None

def germline_class(read, hl, hp_sorted):
    """數 H3 read 覆蓋的 germline het + 一致性 → 分類。"""
    qseq = read.query_sequence
    if qseq is None: return "no_germline", 0, 0.0
    r2q = {rp: qp for qp, rp in read.get_aligned_pairs(matches_only=True)}
    votes = {1: 0, 2: 0}
    lo = bisect.bisect_left(hp_sorted, read.reference_start); hi = bisect.bisect_right(hp_sorted, read.reference_end)
    for rp in hp_sorted[lo:hi]:
        qp = r2q.get(rp)
        if qp is None: continue
        base = qseq[qp].upper(); r, a, g = hl[rp]; al = [r, a]; h1, h2 = al[g[0]], al[g[1]]
        if base == h1 and base != h2: votes[1] += 1
        elif base == h2 and base != h1: votes[2] += 1
    tot = votes[1] + votes[2]
    if tot == 0: return "no_germline", 0, 0.0
    sim = max(votes[1], votes[2]) / tot
    if tot >= 2 and sim < 0.8: return "inconsistent_germline", tot, round(sim, 3)
    return "weak_or_consistent", tot, round(sim, 3)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True); ap.add_argument("--max-sites", type=int, default=300)
    ap.add_argument("--win", type=int, default=2000); ap.add_argument("--min-anchor", type=int, default=8)
    ap.add_argument("--out-dir", default="."); args = ap.parse_args()
    hl = load_hets(args.chrom); hp_sorted = sorted(hl)
    vcf = pysam.VariantFile(SOMATIC_VCF)
    sites = [rec.pos - 1 for rec in vcf.fetch(args.chrom) if rec.alts and len(rec.ref) == 1 and len(rec.alts[0]) == 1]
    vcf.close()
    step = max(1, len(sites) // (args.max_sites * 3)) if sites else 1
    bam = pysam.AlignmentFile(BAM, "rb")
    h3_class = {"no_germline": 0, "inconsistent_germline": 0, "weak_or_consistent": 0}
    h3_assignable = {"no_germline": [], "inconsistent_germline": [], "weak_or_consistent": []}
    n_h3 = n_sites_used = 0
    for c in sites[::step]:
        if n_sites_used >= args.max_sites: break
        ws, we = c - args.win, c + args.win
        g11, g21, h3reads = [], [], []
        seen = set()
        for read in bam.fetch(args.chrom, max(0, ws), we):
            if read.is_secondary or read.is_supplementary or read.is_unmapped: continue
            if read.query_name in seen: continue
            hp = get_hp(read)
            if hp not in ("1-1", "2-1", "3"): continue
            m = read_meth(read, ws, we)
            if m is None: continue
            seen.add(read.query_name)
            if hp == "1-1": g11.append(m)
            elif hp == "2-1": g21.append(m)
            else: h3reads.append((read, m))
        if not h3reads: continue
        n_sites_used += 1
        # 1-1/2-1 質心（指派 H3 用）
        have_centroids = len(g11) >= args.min_anchor and len(g21) >= args.min_anchor
        if have_centroids:
            allp = set()
            for m in g11 + g21: allp.update(m.keys())
            pos = [p for p in allp if sum(1 for m in g11 + g21 if p in m) >= 0.3 * (len(g11) + len(g21))]
            def cen(grp):
                return {p: np.mean([m[p] for m in grp if p in m]) for p in pos if any(p in m for m in grp)}
            c11, c21 = cen(g11), cen(g21)
        for read, m in h3reads:
            n_h3 += 1
            cls, ntot, sim = germline_class(read, hl, hp_sorted)
            h3_class[cls] += 1
            if have_centroids:
                sh = [p for p in m if p in c11 and p in c21]
                if len(sh) >= 5:
                    d11 = np.sqrt(np.mean([(m[p] - c11[p]) ** 2 for p in sh]))
                    d21 = np.sqrt(np.mean([(m[p] - c21[p]) ** 2 for p in sh]))
                    margin = abs(d11 - d21) / (d11 + d21 + 1e-9)
                    h3_assignable[cls].append(round(float(margin), 4))
    bam.close()
    def summ(lst):
        return {"n_with_centroid": len(lst), "margin_median": round(float(np.median(lst)), 4) if lst else None,
                "frac_high_conf(margin>0.1)": round(float(np.mean([1 if x > 0.1 else 0 for x in lst])), 4) if lst else None}
    out = {"chrom": args.chrom, "n_somatic_sites": len(sites), "n_sites_with_H3": n_sites_used, "n_H3_reads": n_h3,
           "H3_germline_class": h3_class,
           "H3_assignability_by_class": {k: summ(v) for k, v in h3_assignable.items()}}
    outp = f"{args.out_dir}/rigor_t2_{args.chrom}.json"
    json.dump(out, open(outp, "w"), indent=1, ensure_ascii=False)
    print(f"[rigor-t2] {args.chrom}: H3 reads={n_h3} | class={h3_class} | "
          f"no_germline assignable margin_med={out['H3_assignability_by_class']['no_germline']['margin_median']} -> {outp}")

if __name__ == "__main__":
    main()
