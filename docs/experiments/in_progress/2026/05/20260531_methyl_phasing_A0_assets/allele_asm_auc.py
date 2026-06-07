#!/usr/bin/env python3
"""
單 SNP REF-vs-ALT allele-specific methylation AUC — tumor vs normal 對等對照。

目的（回答用戶：HP 群內外甲基差異「不是 copy、是真 haplotype ASM」如何驗證）：
  在每個 germline het SNP，把覆蓋該位點的 read 依「帶 REF / 帶 ALT」分兩群（= 兩個 allele），
  算兩群的 per-read 甲基 anchor AUC。這是 allele-specific methylation 的最基本定義：
    - 不需跨 SNP 相位、不靠 longphase HP tag → tumor / normal 完全對等、互相獨立。
    - normal 二倍體 copy-clean、絕大多數位點非 ASM → normal 的 AUC 分布 = 「無 ASM 無 copy」的誠實 null。
    - copy 只影響深度不影響單 read 甲基；若 normal random-het AUC 掉回 ~0.5-0.6 → 方法乾淨、tumor 升高是真 ASM；
      若 normal 也 ~0.97 → 方法樂觀 artifact（與 copy 無關）。

一份腳本同時交付三件：
  (1) 決定性 normal 對照（--bam 指 normal BL）
  (2) CN=2-only 分層（每區記 SEQC2 status，aggregate 時 filter neutral）
  (3) depth-matched downsample（每區把兩群降到等量 read 再算 AUC → 排除 P-06 深度 confound）

唯讀 BAM。tumor/normal 用同一 VCF+chrom+max-sites+seed → 抽到完全相同的 SNP 位置（可配對比較）。
數字全落 JSON；報告須 Read JSON 後才引用（feedback_no_fabricated_numbers_in_reports）。
"""
import sys, json, argparse, bisect
import numpy as np
import pysam
from sklearn.metrics import roc_auc_score

VCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
TUMOR_BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
NORMAL_BAM = "/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam"
CNVDIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
GAIN_LOSS_LOH = CNVDIR + "/ngs_benchmark_cnvs_gain_loss_loh.bed"

# ---------- SEQC2 BED ----------
def load_bed_intervals(path, val_col=3):
    d = {}
    try:
        with open(path) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) <= val_col or not p[0].startswith("chr"):
                    continue
                try:
                    s, e = int(p[1]), int(p[2])
                except ValueError:
                    continue
                d.setdefault(p[0], []).append((s, e, p[val_col].strip()))
    except FileNotFoundError:
        pass
    for c in d:
        d[c].sort()
    return d

def query_interval(intervals, pos):
    if not intervals:
        return None
    starts = [x[0] for x in intervals]
    i = bisect.bisect_right(starts, pos) - 1
    for j in range(max(0, i - 3), min(len(intervals), i + 2)):
        s, e, v = intervals[j]
        if s <= pos <= e:
            return v
    return None

# ---------- 甲基 anchor AUC（沿用 shuffle_control 已驗證邏輯） ----------
def loo_centroid_auc(D, lab):
    lab = np.asarray(lab); n = len(lab)
    m0 = (lab == 0); m1 = (lab == 1); n0 = m0.sum(); n1 = m1.sum()
    if n0 < 2 or n1 < 2:
        return None
    sum0 = D[:, m0].sum(axis=1); sum1 = D[:, m1].sum(axis=1)
    denom0 = np.where(m0, n0 - 1, n0).astype(float); denom1 = np.where(m1, n1 - 1, n1).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        scores = sum0 / denom0 - sum1 / denom1
    v = np.isfinite(scores)
    if v.sum() < 4 or len(set(lab[v].tolist())) < 2:
        return None
    try:
        a = roc_auc_score(lab[v], scores[v])
    except Exception:
        return None
    return max(a, 1 - a)

def shuffle_p95(D, lab, K, rng):
    sh = []
    for _ in range(K):
        sl = rng.permutation(lab)
        if len(set(sl)) < 2:
            continue
        a = loo_centroid_auc(D, sl)
        if a is not None:
            sh.append(a)
    return float(np.percentile(sh, 95)) if sh else None

def build_D(M):
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            mask = ~np.isnan(M[i]) & ~np.isnan(M[j])
            d = np.sqrt(np.mean((M[i, mask] - M[j, mask]) ** 2)) if mask.sum() >= 3 else np.nan
            D[i, j] = D[j, i] = d
    med = np.nanmedian(D[D > 0]) if np.any(D > 0) else 1.0
    D[np.isnan(D)] = med
    return D

# ---------- 單 SNP allele 標籤 + 甲基矩陣 ----------
def region_allele_matrix(bam, chrom, snp_pos0, ref, alt, win, min_anchor, min_cpg):
    """以 het SNP 為錨：read 依該位鹼基分 REF(0)/ALT(1) 兩群，收 ±win CpG 甲基。回 (M, lab) 或 None。"""
    start = max(0, snp_pos0 - win); end = snp_pos0 + win
    rows, seen, allpos = [], set(), set()
    for read in bam.fetch(chrom, start, end):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if read.query_name in seen:
            continue
        seq = read.query_sequence
        if seq is None:
            continue
        q2r = {}
        snp_q = None
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            q2r[qpos] = rpos
            if rpos == snp_pos0:
                snp_q = qpos
        if snp_q is None:
            continue
        base = seq[snp_q].upper()
        if base == ref:
            allele = 0
        elif base == alt:
            allele = 1
        else:
            continue
        mods = read.modified_bases or {}
        mc = None
        for k, calls in mods.items():
            if k[0] in ("C", b"C") and k[2] in ("m", 27551):
                mc = calls; break
        if not mc:
            continue
        meth = {}
        for qpos, qual in mc:
            rpos = q2r.get(qpos)
            if rpos is None or rpos < start or rpos > end:
                continue
            meth[rpos] = qual / 255.0
        if len(meth) < min_cpg:
            continue
        seen.add(read.query_name)
        rows.append((allele, meth)); allpos.update(meth.keys())
    n0 = sum(1 for a, _ in rows if a == 0); n1 = sum(1 for a, _ in rows if a == 1)
    if n0 < min_anchor or n1 < min_anchor:
        return None
    positions = sorted(allpos)
    cov = {p: sum(1 for _, m in rows if p in m) for p in positions}
    positions = [p for p in positions if cov[p] >= 0.3 * len(rows)]
    if len(positions) < min_cpg:
        return None
    M = np.full((len(rows), len(positions)), np.nan)
    lab = np.zeros(len(rows), dtype=int)
    for i, (a, m) in enumerate(rows):
        lab[i] = a
        for j, p in enumerate(positions):
            if p in m:
                M[i, j] = m[p]
    keep = np.array([np.sum(~np.isnan(M[i])) >= min_cpg for i in range(M.shape[0])])
    M, lab = M[keep], lab[keep]
    if (lab == 0).sum() < min_anchor or (lab == 1).sum() < min_anchor:
        return None
    return M, lab

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, choices=["tumor", "normal"])
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--max-sites", type=int, default=150)
    ap.add_argument("--win", type=int, default=2000)
    ap.add_argument("--min-anchor", type=int, default=10)
    ap.add_argument("--min-cpg", type=int, default=5)
    ap.add_argument("--shuffleK", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20260604)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    bam_path = TUMOR_BAM if args.sample == "tumor" else NORMAL_BAM
    rng = np.random.default_rng(args.seed)
    gll = load_bed_intervals(GAIN_LOSS_LOH).get(args.chrom, [])

    # 收集 het SNV（雙等位、單鹼基）位置 — 與 sample 無關，tumor/normal 用同一組
    vcf = pysam.VariantFile(VCF)
    hets = []
    for rec in vcf.fetch(args.chrom):
        if rec.alts is None or len(rec.alts) != 1:
            continue
        if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue
        gt = rec.samples[0].get("GT")
        if not gt or None in gt or len(set(gt)) != 2:
            continue
        hets.append((rec.pos, rec.ref.upper(), rec.alts[0].upper()))
    vcf.close()
    # 均勻抽樣（確定性）：每隔 step 取一個，直到通過篩選達 max-sites
    step = max(1, len(hets) // (args.max_sites * 6)) if hets else 1
    cand = hets[::step]

    bam = pysam.AlignmentFile(bam_path, "rb")
    regions = []
    for pos, ref, alt in cand:
        if len(regions) >= args.max_sites:
            break
        res = region_allele_matrix(bam, args.chrom, pos - 1, ref, alt, args.win, args.min_anchor, args.min_cpg)
        if res is None:
            continue
        M, lab = res
        D = build_D(M)
        auc = loo_centroid_auc(D, lab)
        if auc is None:
            continue
        sp95 = shuffle_p95(D, lab, args.shuffleK, rng)
        # depth-matched：兩群降到等量 read
        n0 = int((lab == 0).sum()); n1 = int((lab == 1).sum())
        k = min(n0, n1)
        auc_dm = None
        if k >= args.min_anchor:
            idx0 = rng.choice(np.where(lab == 0)[0], k, replace=False)
            idx1 = rng.choice(np.where(lab == 1)[0], k, replace=False)
            idx = np.sort(np.concatenate([idx0, idx1]))
            Dm = D[np.ix_(idx, idx)]; lm = lab[idx]
            auc_dm = loo_centroid_auc(Dm, lm)
        status = query_interval(gll, pos) or "neutral"
        regions.append({
            "pos": pos, "ref": ref, "alt": alt,
            "auc": round(float(auc), 4),
            "auc_depthmatched": round(float(auc_dm), 4) if auc_dm is not None else None,
            "shuffle_p95": round(float(sp95), 4) if sp95 is not None else None,
            "sig_over_shuffle": bool(sp95 is not None and auc > sp95),
            "n_ref": n0, "n_alt": n1, "n_cpg": int(M.shape[1]),
            "status": status,
        })
    bam.close()

    aucs = [r["auc"] for r in regions]
    out = {
        "sample": args.sample, "chrom": args.chrom, "bam": bam_path,
        "n_het_candidates": len(hets), "n_regions": len(regions),
        "auc_median": round(float(np.median(aucs)), 4) if aucs else None,
        "auc_mean": round(float(np.mean(aucs)), 4) if aucs else None,
        "frac_sig_over_shuffle": round(float(np.mean([1 if r["sig_over_shuffle"] else 0 for r in regions])), 4) if regions else None,
        "regions": regions,
    }
    outp = f"{args.out_dir}/allele_asm_{args.sample}_{args.chrom}.json"
    with open(outp, "w") as f:
        json.dump(out, f, indent=1)
    print(f"[allele_asm] {args.sample} {args.chrom}: {len(regions)} regions, auc_median={out['auc_median']}, frac_sig={out['frac_sig_over_shuffle']} -> {outp}")

if __name__ == "__main__":
    main()
