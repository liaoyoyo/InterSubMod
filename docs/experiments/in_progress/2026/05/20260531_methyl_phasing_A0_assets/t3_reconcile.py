#!/usr/bin/env python3
"""
調和 local-allele(POSITIVE) vs rigor_t3 HP-tag(NEGATIVE)，並查 confound。
對 somatic 位點的 HP1 ALT(亞群) vs REF(母本) 覆蓋 read：
 (1) matched Δβ(ALT,REF) vs 同窗噪音(split REF) → 用 rigor_t3 的 metric 測；若也正 → 差異是「標籤」(HP-tag 被污染) 非 metric。
 (2) farCpG 版 Δβ（排除突變 ±100bp）。
 (3) confound：ALT/REF 兩群 per-read n_cpg + read 長度中位數（差太多→覆蓋 confound）。
唯讀。落 t3_reconcile_{chrom}.json。
"""
import sys, json, argparse
import numpy as np
import pysam
from scipy.stats import wilcoxon

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
SVCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/somatic_pass.vcf.gz"

def get_hp(r):
    try: return str(r.get_tag("HP"))
    except KeyError: return "u"

def dbeta(g0, g1, excl=None, pos0=None):
    pos = set()
    for m in g0 + g1: pos.update(m.keys())
    if excl is not None: pos = {p for p in pos if abs(p - pos0) > excl}
    pos = [p for p in pos if sum(1 for m in g0 + g1 if p in m) >= 0.3 * (len(g0) + len(g1))]
    if len(pos) < 5: return None
    def cen(g):
        return np.array([np.mean([m[p] for m in g if p in m]) if any(p in m for m in g) else np.nan for p in pos])
    c0, c1 = cen(g0), cen(g1); mask = ~np.isnan(c0) & ~np.isnan(c1)
    return float(np.mean(np.abs(c0[mask] - c1[mask]))) if mask.sum() >= 5 else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True); ap.add_argument("--max-sites", type=int, default=200)
    ap.add_argument("--win", type=int, default=4000); ap.add_argument("--min", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260611); ap.add_argument("--out-dir", default="."); a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    vcf = pysam.VariantFile(SVCF)
    sites = [(r.pos - 1, r.ref.upper(), r.alts[0].upper()) for r in vcf.fetch(a.chrom)
             if r.alts and len(r.ref) == 1 and len(r.alts[0]) == 1]; vcf.close()
    step = max(1, len(sites) // (a.max_sites * 3)) if sites else 1
    bam = pysam.AlignmentFile(BAM, "rb")
    d_real, d_noise, d_far, ncpg_alt, ncpg_ref, len_alt, len_ref = [], [], [], [], [], [], []
    for pos0, ref, alt in sites[::step]:
        if len(d_real) >= a.max_sites: break
        ws, we = pos0 - a.win, pos0 + a.win
        ALT, REF, lA, lR = [], [], [], []
        seen = set()
        for read in bam.fetch(a.chrom, max(0, ws), we):
            if read.is_secondary or read.is_supplementary or read.is_unmapped: continue
            if read.query_name in seen or get_hp(read) != "1": continue
            qseq = read.query_sequence
            if qseq is None: continue
            pairs = read.get_aligned_pairs(matches_only=True)
            r2q = {rp: qp for qp, rp in pairs}
            if pos0 not in r2q: continue
            base = qseq[r2q[pos0]].upper()
            if base not in (ref, alt): continue
            mods = read.modified_bases or {}; mc = None
            for k, calls in mods.items():
                if k[0] in ("C", b"C") and k[2] in ("m", 27551): mc = calls; break
            if not mc: continue
            q2r = {q: r for q, r in pairs}; meth = {}
            for qp, ql in mc:
                rp = q2r.get(qp)
                if rp is None or rp < ws or rp > we: continue
                meth[rp] = ql / 255.0
            if len(meth) < 5: continue
            seen.add(read.query_name)
            rl = read.reference_length or 0
            if base == alt: ALT.append(meth); lA.append(rl)
            else: REF.append(meth); lR.append(rl)
        if len(ALT) < a.min or len(REF) < a.min: continue
        dr = dbeta(ALT, REF);
        if dr is None: continue
        idx = rng.permutation(len(REF)); h = len(REF) // 2
        dn = dbeta([REF[i] for i in idx[:h]], [REF[i] for i in idx[h:]])
        df = dbeta(ALT, REF, excl=100, pos0=pos0)
        d_real.append(dr)
        if dn is not None: d_noise.append(dn)
        if df is not None: d_far.append(df)
        ncpg_alt.append(np.median([len(m) for m in ALT])); ncpg_ref.append(np.median([len(m) for m in REF]))
        len_alt.append(np.median(lA)); len_ref.append(np.median(lR))
    bam.close()
    # matched Δβ(real)−noise 同窗配對（real 與 noise 同一組 sites 子集；取共同 index）
    n = min(len(d_real), len(d_noise))
    paired = [d_real[i] - d_noise[i] for i in range(n)]
    out = {"chrom": a.chrom, "n_windows": len(d_real),
           "dbeta_real_median": round(float(np.median(d_real)), 4) if d_real else None,
           "dbeta_noise_median": round(float(np.median(d_noise)), 4) if d_noise else None,
           "dbeta_farCpG_median": round(float(np.median(d_far)), 4) if d_far else None,
           "matched_real_minus_noise_median": round(float(np.median(paired)), 4) if paired else None,
           "matched_wilcoxon_p": float(wilcoxon(paired).pvalue) if len(paired) >= 10 else None,
           "ncpg_alt_median": round(float(np.median(ncpg_alt)), 1) if ncpg_alt else None,
           "ncpg_ref_median": round(float(np.median(ncpg_ref)), 1) if ncpg_ref else None,
           "readlen_alt_median": int(np.median(len_alt)) if len_alt else None,
           "readlen_ref_median": int(np.median(len_ref)) if len_ref else None}
    json.dump(out, open(f"{a.out_dir}/t3_reconcile_{a.chrom}.json", "w"), indent=1, ensure_ascii=False)
    print(f"[recon] {a.chrom} n={out['n_windows']}: Δβ real={out['dbeta_real_median']} noise={out['dbeta_noise_median']} far={out['dbeta_farCpG_median']} "
          f"matched(real-noise)={out['matched_real_minus_noise_median']} p={out['matched_wilcoxon_p']} | "
          f"nCpG ALT={out['ncpg_alt_median']}/REF={out['ncpg_ref_median']} len ALT={out['readlen_alt_median']}/REF={out['readlen_ref_median']}")

if __name__ == "__main__":
    main()
