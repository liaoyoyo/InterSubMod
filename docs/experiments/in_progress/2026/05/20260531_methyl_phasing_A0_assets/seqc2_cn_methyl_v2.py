#!/usr/bin/env python3
"""
v2：在 seqc2_cn_methyl.py 基礎上加三件事（用戶 2026-06-01 a/b/c）：
  (a) CpG-SNP 排除控制：剔除「會被 het SNP 改變的 CpG 欄」，重算 anchor AUC
      → 比 raw AUC，看 95% 訊號扣掉 CpG-SNP artifact 後還剩多少
  (b) 速度：per-chr 內部不再每區開 samtools subprocess（已改 pysam）；外層 38 平行
  (c) LOH 個案細節：每個 LOH region 輸出 raw/cpgsnp-excluded AUC + delta + n，
      供事後篩「哪些 LOH 反而分得開」+ 異常標記

CpG-SNP 排除邏輯（保守）：
  - 載入該區所有 germline het SNP（含本 phased VCF）。
  - 一個 CpG（ref C 在 +strand / G 在 -strand）若其 C 位、或其下游 G 位（CpG 的 G）落在 het SNP 上
    → 該 SNP 可能破壞/創造 CpG 或讓 allele 序列本身不同 → 標為 cpgsnp-tainted，重算時剔除該欄。
  - 同時剔除「位置正好是 het SNP」的 CpG（最直接的 pseudo-ASM 來源）。

輸出每區同時含 anchor_auc(raw) 與 anchor_auc_cpgsnp_excl，可比較。
唯讀。每 chr 落 seqc2_cn_methyl_v2_{chrom}.json。
"""
import sys, json, argparse, bisect
import numpy as np
import pysam
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GVCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
REF = "/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
CNVDIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
AD = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets"

def load_bed(path, val_col=3, as_int=False):
    d = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= val_col or not p[0].startswith("chr"):
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            v = p[val_col].strip()
            if as_int:
                try: v = int(v)
                except ValueError: continue
            d.setdefault(p[0], []).append((s, e, v))
    for c in d: d[c].sort()
    return d

def query_iv(iv, pos):
    if not iv: return None
    starts = [x[0] for x in iv]
    i = bisect.bisect_right(starts, pos) - 1
    for j in range(max(0, i-3), min(len(iv), i+2)):
        s, e, v = iv[j]
        if s <= pos <= e: return v
    return None

def get_hp(read):
    try: return str(read.get_tag("HP"))
    except KeyError: return "unphase"

def loo_centroid_auc(D, lab):
    lab = np.asarray(lab); n = len(lab)
    m0 = (lab == 0); m1 = (lab == 1)
    n0 = m0.sum(); n1 = m1.sum()
    if n0 < 2 or n1 < 2: return None
    sum0 = D[:, m0].sum(axis=1); sum1 = D[:, m1].sum(axis=1)
    denom0 = np.where(m0, n0-1, n0).astype(float)
    denom1 = np.where(m1, n1-1, n1).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        scores = sum0/denom0 - sum1/denom1
    v = np.isfinite(scores)
    if v.sum() < 4 or len(set(lab[v].tolist())) < 2: return None
    try: a = roc_auc_score(lab[v], scores[v])
    except Exception: return None
    return max(a, 1-a)

def pairwise_D(M, min_overlap=3):
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            mask = ~np.isnan(M[i]) & ~np.isnan(M[j])
            d = np.sqrt(np.mean((M[i, mask]-M[j, mask])**2)) if mask.sum() >= min_overlap else np.nan
            D[i, j] = D[j, i] = d
    med = np.nanmedian(D[D > 0]) if np.any(D > 0) else 1.0
    D[np.isnan(D)] = med
    return D

def shuffle_p95(D, lab, K, rng):
    sh = []
    for _ in range(K):
        sl = rng.permutation(lab)
        if len(set(sl)) < 2: continue
        a = loo_centroid_auc(D, sl)
        if a is not None: sh.append(a)
    return float(np.percentile(sh, 95)) if sh else None

def auc_from_matrix(M, lab, K, rng):
    """給 read×CpG matrix + lab → (auc, p95)。"""
    keep = np.array([np.sum(~np.isnan(M[i])) >= 3 for i in range(M.shape[0])])
    M2, lab2 = M[keep], lab[keep]
    if (lab2 == 0).sum() < 5 or (lab2 == 1).sum() < 5 or M2.shape[1] < 3:
        return None, None, int((lab2 == 0).sum()), int((lab2 == 1).sum())
    D = pairwise_D(M2)
    auc = loo_centroid_auc(D, lab2)
    p95 = shuffle_p95(D, lab2, K, rng) if auc is not None else None
    return auc, p95, int((lab2 == 0).sum()), int((lab2 == 1).sum())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--max-per-status", type=int, default=25)
    ap.add_argument("--K", type=int, default=30)
    ap.add_argument("--win", type=int, default=2000)
    args = ap.parse_args()
    chrom = args.chrom
    rng = np.random.RandomState(20260601)

    type_iv = load_bed(CNVDIR+"/ngs_benchmark_cnvs_gain_loss_loh.bed", 3, False).get(chrom, [])
    gain_iv = load_bed(CNVDIR+"/ngs_benchmark_cnv_gain_cn.bed", 3, True).get(chrom, [])
    loss_iv = load_bed(CNVDIR+"/ngs_benchmark_cnv_loss_cn.bed", 3, True).get(chrom, [])

    bam = pysam.AlignmentFile(BAM, "rb")
    vcf = pysam.VariantFile(GVCF)
    ref = pysam.FastaFile(REF)

    # 收 het SNP positions（全 chr，供 region 內 CpG-SNP 判定 + 取樣）
    het_pos = []
    het_set = set()
    for rec in vcf.fetch(chrom):
        gt = None
        for s in rec.samples.values(): gt = s.get("GT"); break
        if gt is None or len(gt) != 2 or gt[0] == gt[1] or None in gt: continue
        het_pos.append(rec.pos)
        het_set.add(rec.pos)  # 1-based
    het_sorted = sorted(het_pos)
    order = list(het_pos); rng.shuffle(order)

    def het_in_window(s, e):
        lo = bisect.bisect_left(het_sorted, s)
        hi = bisect.bisect_right(het_sorted, e)
        return set(het_sorted[lo:hi])

    results = []
    sc = {"gain": 0, "loss": 0, "loh": 0, "neutral": 0}
    for pos in order:
        if all(v >= args.max_per_status for v in sc.values()):
            break
        typ = query_iv(type_iv, pos)
        status = typ if typ in ("gain", "loss", "loh") else "neutral"
        if sc[status] >= args.max_per_status:
            continue
        start, end = pos - args.win, pos + args.win
        # 該區 het SNP（0-based ref 座標集，用於 CpG-SNP 判定）
        win_hets = het_in_window(start, end)
        het0 = set(p - 1 for p in win_hets)  # 0-based

        # 收 read × CpG
        rows, seen, allpos = [], set(), set()
        for read in bam.fetch(chrom, start, end):
            if read.is_secondary or read.is_supplementary or read.is_unmapped: continue
            if read.query_name in seen: continue
            hp = get_hp(read)
            if hp not in ("1", "2"): continue
            mods = read.modified_bases or {}
            mc = None
            for k, calls in mods.items():
                if k[0] in ("C", b"C") and k[2] in ("m", 27551):
                    mc = calls; break
            if not mc: continue
            q2r = {q: r for q, r in read.get_aligned_pairs(matches_only=True)}
            meth = {}
            for qpos, qual in mc:
                rpos = q2r.get(qpos)
                if rpos is None or rpos < start or rpos > end: continue
                meth[rpos] = qual / 255.0
            if len(meth) < 5: continue
            seen.add(read.query_name)
            rows.append((hp, meth)); allpos.update(meth.keys())

        n1 = sum(1 for h, _ in rows if h == "1")
        n2 = sum(1 for h, _ in rows if h == "2")
        if n1 < 10 or n2 < 10: continue
        positions = sorted(allpos)
        cov = {p: sum(1 for _, m in rows if p in m) for p in positions}
        positions = [p for p in positions if cov[p] >= 0.3 * len(rows)]
        if len(positions) < 5: continue

        # 判定每個 CpG 欄是否 cpgsnp-tainted
        # CpG: ref pos 0-based. C 在 +strand → ref[pos]=C 且 ref[pos+1]=G。
        # tainted 若: pos 或 pos+1 落在 het SNP；或 pos-1(若該 CpG 由 -strand G 表示)落 het。
        tainted = set()
        for p in positions:
            # p 為 0-based ref（pysam aligned_pairs 給 0-based）
            if p in het0 or (p + 1) in het0 or (p - 1) in het0:
                tainted.add(p)
                continue
            # 額外：檢查 ref 序列該位是否真 CpG，且鄰位 het 破壞之（保守已由上 cover）
        clean_pos = [p for p in positions if p not in tainted]

        def build_M(pos_list):
            M = np.full((len(rows), len(pos_list)), np.nan)
            lab = np.zeros(len(rows), dtype=int)
            rowmean = np.full(len(rows), np.nan)
            for i, (h, m) in enumerate(rows):
                lab[i] = 0 if h == "1" else 1
                vals = []
                for j, p in enumerate(pos_list):
                    if p in m:
                        M[i, j] = m[p]; vals.append(m[p])
                if vals: rowmean[i] = np.mean(vals)
            return M, lab, rowmean

        M_raw, lab, rowmean = build_M(positions)
        auc_raw, p95_raw, kn1, kn2 = auc_from_matrix(M_raw, lab, args.K, np.random.RandomState(int(pos) % 99999))
        if auc_raw is None: continue

        # CpG-SNP excluded
        n_tainted = len(tainted)
        if clean_pos and len(clean_pos) >= 5:
            M_clean, lab_c, _ = build_M(clean_pos)
            auc_clean, p95_clean, _, _ = auc_from_matrix(M_clean, lab_c, args.K, np.random.RandomState(int(pos) % 99999 + 7))
        else:
            auc_clean, p95_clean = None, None

        # GMM 雙峰
        x = rowmean[~np.isnan(rowmean)].reshape(-1, 1)
        bic = None
        if len(x) >= 12:
            try:
                bic = float(GaussianMixture(1, random_state=0).fit(x).bic(x) -
                            GaussianMixture(2, random_state=0).fit(x).bic(x))
            except Exception: bic = None

        cn = query_iv(gain_iv, pos) if status == "gain" else (query_iv(loss_iv, pos) if status == "loss" else None)
        try:
            covarr = bam.count_coverage(chrom, start, end, quality_threshold=0)
            depth = float(np.array(covarr, dtype=np.int64).sum(axis=0).mean())
        except Exception:
            depth = None

        sc[status] += 1
        results.append({
            "chrom": chrom, "pos": pos, "seqc2_status": status, "seqc2_cn": cn,
            "n_hp1": n1, "n_hp2": n2, "n_cpg_raw": len(positions),
            "n_cpg_tainted": n_tainted, "n_cpg_clean": len(clean_pos),
            "anchor_auc_raw": round(auc_raw, 4),
            "shuffle_p95_raw": round(p95_raw, 4) if p95_raw is not None else None,
            "anchor_auc_cpgsnp_excl": round(auc_clean, 4) if auc_clean is not None else None,
            "shuffle_p95_excl": round(p95_clean, 4) if p95_clean is not None else None,
            "auc_drop_after_excl": round(auc_raw - auc_clean, 4) if auc_clean is not None else None,
            "mean_depth": round(depth, 1) if depth is not None else None,
            "gmm_bic_diff": round(bic, 2) if bic is not None else None,
        })
    bam.close()

    out = {"chrom": chrom, "n_regions": len(results), "status_count": sc,
           "params": {"win": args.win, "K": args.K, "max_per_status": args.max_per_status},
           "regions": results}
    outpath = f"{AD}/seqc2_cn_methyl_v2_{chrom}.json"
    json.dump(out, open(outpath, "w"), ensure_ascii=False, indent=2)
    print(f"[v2] {chrom}: {len(results)} regions status={sc} -> {outpath}")

if __name__ == "__main__":
    main()
