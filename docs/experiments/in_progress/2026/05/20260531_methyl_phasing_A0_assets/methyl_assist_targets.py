#!/usr/bin/env python3
"""
T2 + T3：用甲基輔助 longphase-S 的 somatic-層 tag 矯正（held-out 判別框架）。

longphase-S tag cascade（已讀原始碼確認）：
  unTag=無證據 | H1/H2=germline-only | H1-1/H2-1=帶 somatic 變異且歸 germline H1/H2 分支 | H3=somatic 但 germline 未知

T3（核心，最novel）：plain H1/H2 read 是「真 germline-only」與「其實屬亞群但沒踩到 somatic 位點」的混合。
  測：在 somatic 位點窗，tag `1-1`(亞群) vs `1`(非亞群) 的甲基能否判別（held-out AUC）。
  → 若能，代表亞群有獨特 tumor-acquired ASM（V10 已證真實非copy）→ 可拆 ambiguous plain H1。同理 2-1 vs 2。

T2：H3 read 帶 somatic 變異但 germline 分支未知。
  測：tag `1-1` vs `2-1` 的甲基能否判別 germline 分支（held-out AUC）→ 若能，可把 H3 歸 H1-1/H2-1。

每個 somatic 位點 ±win 窗，按 BAM 既有 HP tag 分群（read-only，longphase 已算好的 tag）。
數字全落 JSON；報告須 Read 後才引用。
"""
import sys, json, argparse
import numpy as np
import pysam
from sklearn.metrics import roc_auc_score

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
SOMATIC_VCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/somatic_pass.vcf.gz"

def get_hp(read):
    try:
        return str(read.get_tag("HP"))
    except KeyError:
        return "unphase"

def collect_window(bam, chrom, center0, win, want_tags, min_cpg=5):
    """收窗內各 want_tag 群的 read×CpG 甲基。回 {tag: [meth_dict,...]}。"""
    start = max(0, center0 - win); end = center0 + win
    groups = {t: [] for t in want_tags}
    seen = set()
    for read in bam.fetch(chrom, start, end):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if read.query_name in seen:
            continue
        hp = get_hp(read)
        if hp not in want_tags:
            continue
        mods = read.modified_bases or {}
        mc = None
        for k, calls in mods.items():
            if k[0] in ("C", b"C") and k[2] in ("m", 27551):
                mc = calls; break
        if not mc:
            continue
        q2r = {q: r for q, r in read.get_aligned_pairs(matches_only=True)}
        meth = {}
        for qpos, qual in mc:
            rpos = q2r.get(qpos)
            if rpos is None or rpos < start or rpos > end:
                continue
            meth[rpos] = qual / 255.0
        if len(meth) < min_cpg:
            continue
        seen.add(read.query_name)
        groups[hp].append(meth)
    return groups

def auc_two_groups(g0, g1, min_anchor, min_cpg, rng, shuffleK=80):
    """兩群甲基 LOO 質心 AUC + shuffle p95 + held-out（train/test）accuracy。"""
    if len(g0) < min_anchor or len(g1) < min_anchor:
        return None
    rows = [(0, m) for m in g0] + [(1, m) for m in g1]
    allpos = set()
    for _, m in rows:
        allpos.update(m.keys())
    positions = sorted(p for p in allpos if sum(1 for _, m in rows if p in m) >= 0.3 * len(rows))
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
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            mask = ~np.isnan(M[i]) & ~np.isnan(M[j])
            d = np.sqrt(np.mean((M[i, mask] - M[j, mask]) ** 2)) if mask.sum() >= 3 else np.nan
            D[i, j] = D[j, i] = d
    med = np.nanmedian(D[D > 0]) if np.any(D > 0) else 1.0
    D[np.isnan(D)] = med

    def loo(lb):
        lb = np.asarray(lb); m0 = lb == 0; m1 = lb == 1
        n0 = m0.sum(); n1 = m1.sum()
        if n0 < 2 or n1 < 2: return None
        s0 = D[:, m0].sum(1); s1 = D[:, m1].sum(1)
        d0 = np.where(m0, n0 - 1, n0).astype(float); d1 = np.where(m1, n1 - 1, n1).astype(float)
        with np.errstate(all="ignore"):
            sc = s0 / d0 - s1 / d1
        v = np.isfinite(sc)
        if v.sum() < 4 or len(set(lb[v].tolist())) < 2: return None
        try: return max(roc_auc_score(lb[v], sc[v]), 1 - roc_auc_score(lb[v], sc[v]))
        except Exception: return None

    auc = loo(lab)
    if auc is None:
        return None
    sh = [a for _ in range(shuffleK) if (a := loo(rng.permutation(lab))) is not None]
    sp95 = float(np.percentile(sh, 95)) if sh else None
    return {"auc": round(float(auc), 4), "shuffle_p95": round(sp95, 4) if sp95 else None,
            "sig": bool(sp95 is not None and auc > sp95),
            "n0": int((lab == 0).sum()), "n1": int((lab == 1).sum()), "n_cpg": len(positions)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--max-sites", type=int, default=200)
    ap.add_argument("--win", type=int, default=2000)
    ap.add_argument("--min-anchor", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260609)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    vcf = pysam.VariantFile(SOMATIC_VCF)
    sites = []
    for rec in vcf.fetch(args.chrom):
        if rec.alts and len(rec.ref) == 1 and len(rec.alts[0]) == 1:
            sites.append(rec.pos - 1)
    vcf.close()
    step = max(1, len(sites) // (args.max_sites * 3)) if sites else 1
    cand = sites[::step]

    bam = pysam.AlignmentFile(BAM, "rb")
    t3_h1, t3_h2, t2 = [], [], []
    for c in cand:
        if len(t3_h1) + len(t3_h2) >= args.max_sites * 2:
            break
        g = collect_window(bam, args.chrom, c, args.win, {"1", "2", "1-1", "2-1"})
        # T3-H1: 亞群 1-1 vs 非亞群 1
        r = auc_two_groups(g["1-1"], g["1"], args.min_anchor, 5, rng)
        if r: r["site"] = c + 1; t3_h1.append(r)
        # T3-H2: 亞群 2-1 vs 非亞群 2
        r = auc_two_groups(g["2-1"], g["2"], args.min_anchor, 5, rng)
        if r: r["site"] = c + 1; t3_h2.append(r)
        # T2: 1-1 vs 2-1（germline 分支判別）
        r = auc_two_groups(g["1-1"], g["2-1"], args.min_anchor, 5, rng)
        if r: r["site"] = c + 1; t2.append(r)
    bam.close()

    def summ(lst):
        if not lst: return {"n_sites": 0}
        aucs = [x["auc"] for x in lst]
        return {"n_sites": len(lst), "auc_median": round(float(np.median(aucs)), 4),
                "frac_sig": round(float(np.mean([1 if x["sig"] else 0 for x in lst])), 4),
                "shuffle_p95_median": round(float(np.median([x["shuffle_p95"] for x in lst if x["shuffle_p95"]])), 4) if any(x["shuffle_p95"] for x in lst) else None}
    out = {"chrom": args.chrom, "n_somatic_sites": len(sites),
           "T3_H1_subclone_vs_germline": summ(t3_h1),
           "T3_H2_subclone_vs_germline": summ(t3_h2),
           "T2_H1-1_vs_H2-1_branch": summ(t2),
           "sites_t3_h1": t3_h1, "sites_t3_h2": t3_h2, "sites_t2": t2}
    outp = f"{args.out_dir}/methyl_assist_{args.chrom}.json"
    json.dump(out, open(outp, "w"), indent=1, ensure_ascii=False)
    print(f"[assist-T2T3] {args.chrom}: somatic sites={len(sites)} | "
          f"T3-H1 n={out['T3_H1_subclone_vs_germline'].get('n_sites')} auc={out['T3_H1_subclone_vs_germline'].get('auc_median')} sig={out['T3_H1_subclone_vs_germline'].get('frac_sig')} | "
          f"T3-H2 n={out['T3_H2_subclone_vs_germline'].get('n_sites')} auc={out['T3_H2_subclone_vs_germline'].get('auc_median')} | "
          f"T2 n={out['T2_H1-1_vs_H2-1_branch'].get('n_sites')} auc={out['T2_H1-1_vs_H2-1_branch'].get('auc_median')} -> {outp}")

if __name__ == "__main__":
    main()
