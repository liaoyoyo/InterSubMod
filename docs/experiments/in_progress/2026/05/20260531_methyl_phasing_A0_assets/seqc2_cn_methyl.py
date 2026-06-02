#!/usr/bin/env python3
"""
SEQC2 CN/LOH × 甲基分群一致性分析（Q1+Q2+Q3 共用資料基礎）。

回答用戶三問：
  Q1 甲基分群 AUC 是否與 SEQC2 CN-gain/LOH 狀態高度相關？（若是→95% 訊號疑為 CN confound）
  Q2 高甲基分離是否集中在 CN-gain/低 anchor 區（倍體+alignment artifact）？
  Q3 甲基(coverage/分群)能否預測幾倍體？

方法（每個 germline het SNP ±2kb window）：
  - 甲基 anchor AUC（HP1 vs HP2 leave-one-out 質心，對稱化）+ K-shuffle p95（沿用已驗證法，可比 null）
  - SEQC2 標註：type ∈ {gain,loss,loh,neutral} + CN 整數值（gain_cn/loss_cn BED）
  - local mean coverage（samtools depth；驗證 SEQC2 CN 在本 BAM 成立 = Q3 的 CN proxy）
  - GMM bimodality（1 vs 2 component BIC 差）+ diptest（per-read mean methylation）= Q2 multimodal

唯讀 BAM。每 chr 落一個 JSON：seqc2_cn_methyl_{chrom}.json。
數字全由本腳本算出落檔；報告須 Read JSON 後才引用（feedback_no_fabricated_numbers_in_reports）。
"""
import sys, json, argparse, bisect, subprocess
import numpy as np
import pysam
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture
try:
    import diptest as _diptest
    HAVE_DIP = True
except ImportError:
    HAVE_DIP = False

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GERMLINE_VCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
CNVDIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
GAIN_LOSS_LOH = CNVDIR + "/ngs_benchmark_cnvs_gain_loss_loh.bed"
GAIN_CN = CNVDIR + "/ngs_benchmark_cnv_gain_cn.bed"
LOSS_CN = CNVDIR + "/ngs_benchmark_cnv_loss_cn.bed"

# ---------- SEQC2 BED 載入（per-chr sorted intervals） ----------
def load_bed_intervals(path, val_col=3, val_is_int=False):
    """回傳 {chrom: [(start,end,val),...] sorted by start}。"""
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
            if val_is_int:
                try:
                    v = int(v)
                except ValueError:
                    continue
            d.setdefault(p[0], []).append((s, e, v))
    for c in d:
        d[c].sort()
    return d

def query_interval(intervals, pos):
    """intervals=[(s,e,v)] sorted；回傳覆蓋 pos 的 v（第一個），無則 None。"""
    if not intervals:
        return None
    starts = [x[0] for x in intervals]
    i = bisect.bisect_right(starts, pos) - 1
    # 往前看幾個（區間可能重疊）
    for j in range(max(0, i - 3), min(len(intervals), i + 2)):
        s, e, v = intervals[j]
        if s <= pos <= e:
            return v
    return None

# ---------- 甲基 anchor AUC（沿用 shuffle_control 邏輯） ----------
def get_hp(read):
    try:
        return str(read.get_tag("HP"))
    except KeyError:
        return "unphase"

def region_D_lab(bam, chrom, start, end, min_anchor=10, min_cpg=5):
    rows, seen, allpos = [], set(), set()
    for read in bam.fetch(chrom, start, end):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if read.query_name in seen:
            continue
        hp = get_hp(read)
        if hp not in ("1", "2"):
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
        rows.append((hp, meth))
        allpos.update(meth.keys())
    n1 = sum(1 for h, _ in rows if h == "1")
    n2 = sum(1 for h, _ in rows if h == "2")
    if n1 < min_anchor or n2 < min_anchor:
        return None
    positions = sorted(allpos)
    cov = {p: sum(1 for _, m in rows if p in m) for p in positions}
    positions = [p for p in positions if cov[p] >= 0.3 * len(rows)]
    if len(positions) < min_cpg:
        return None
    M = np.full((len(rows), len(positions)), np.nan)
    lab = np.zeros(len(rows), dtype=int)
    rowmean = np.full(len(rows), np.nan)
    for i, (h, m) in enumerate(rows):
        lab[i] = 0 if h == "1" else 1
        vals = []
        for j, p in enumerate(positions):
            if p in m:
                M[i, j] = m[p]; vals.append(m[p])
        if vals:
            rowmean[i] = np.mean(vals)
    keep = np.array([np.sum(~np.isnan(M[i])) >= min_cpg for i in range(M.shape[0])])
    M, lab, rowmean = M[keep], lab[keep], rowmean[keep]
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
    return D, lab, int((lab == 0).sum()), int((lab == 1).sum()), len(positions), rowmean

def loo_centroid_auc(D, lab):
    """向量化 leave-one-out 質心距離 AUC（同邏輯，O(n²) 一次矩陣運算，shuffle K 次才不會慢）。"""
    lab = np.asarray(lab)
    n = len(lab)
    m0 = (lab == 0); m1 = (lab == 1)
    n0 = m0.sum(); n1 = m1.sum()
    if n0 < 2 or n1 < 2:
        return None
    # 對 read i：到 group0 其他成員的均距 = (sum_{k in g0} D[i,k] - (i in g0 ? D[i,i]=0 : 0)) / (g0 大小 - 自身)
    sum0 = D[:, m0].sum(axis=1)  # 含自身(D[i,i]=0 不影響和)
    sum1 = D[:, m1].sum(axis=1)
    # leave-one-out：若 i 屬 g0，分母 n0-1；否則 n0
    denom0 = np.where(m0, n0 - 1, n0).astype(float)
    denom1 = np.where(m1, n1 - 1, n1).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        d0 = sum0 / denom0
        d1 = sum1 / denom1
    scores = d0 - d1  # 越大越像 group1
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

def bimodality(rowmean):
    """per-read mean methylation 的多模態：GMM BIC(2)−BIC(1) + diptest p。"""
    x = rowmean[~np.isnan(rowmean)].reshape(-1, 1)
    if len(x) < 12:
        return None, None
    try:
        g1 = GaussianMixture(1, random_state=0).fit(x)
        g2 = GaussianMixture(2, random_state=0).fit(x)
        bic_diff = float(g1.bic(x) - g2.bic(x))  # >0 表 2-comp 較優 = 多模態
    except Exception:
        bic_diff = None
    dp = None
    if HAVE_DIP:
        try:
            _, dp = _diptest.diptest(x.ravel())
            dp = float(dp)
        except Exception:
            dp = None
    return bic_diff, dp

def region_mean_depth(bam, chrom, start, end):
    """平均覆蓋（CN proxy）— 用已開啟的 pysam handle in-process count_coverage（比 subprocess samtools depth 快數十倍，避免大染色體卡死）。"""
    try:
        cov = bam.count_coverage(chrom, start, end, quality_threshold=0)
        # cov = 4 arrays (A,C,G,T) per base；相加 = per-base depth
        per_base = np.array(cov, dtype=np.int64).sum(axis=0)
        return float(per_base.mean()) if len(per_base) else 0.0
    except Exception:
        return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out-dir", default="/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets")
    ap.add_argument("--max-per-status", type=int, default=40, help="每類 SEQC2 status 最多取樣 region 數")
    ap.add_argument("--K", type=int, default=50)
    ap.add_argument("--win", type=int, default=2000)
    args = ap.parse_args()
    chrom = args.chrom
    rng = np.random.RandomState(20260601)

    type_iv = load_bed_intervals(GAIN_LOSS_LOH, 3, False).get(chrom, [])
    gain_cn_iv = load_bed_intervals(GAIN_CN, 3, True).get(chrom, [])
    loss_cn_iv = load_bed_intervals(LOSS_CN, 3, True).get(chrom, [])

    bam = pysam.AlignmentFile(BAM, "rb")
    vcf = pysam.VariantFile(GERMLINE_VCF)

    # 收集 het SNP
    hets = []
    for rec in vcf.fetch(chrom):
        gt = None
        for s in rec.samples.values():
            gt = s.get("GT"); break
        if gt is None or len(gt) != 2 or gt[0] == gt[1] or None in gt:
            continue
        hets.append(rec.pos)
    rng.shuffle(hets)

    results = []
    status_count = {"gain": 0, "loss": 0, "loh": 0, "neutral": 0}
    for pos in hets:
        if all(v >= args.max_per_status for v in status_count.values()):
            break
        typ = query_interval(type_iv, pos)
        status = typ if typ in ("gain", "loss", "loh") else "neutral"
        if status_count[status] >= args.max_per_status:
            continue
        rd = region_D_lab(bam, chrom, pos - args.win, pos + args.win)
        if rd is None:
            continue
        D, lab, n1, n2, ncpg, rowmean = rd
        auc = loo_centroid_auc(D, lab)
        if auc is None:
            continue
        p95 = shuffle_p95(D, lab, args.K, rng)
        bic_diff, dip_p = bimodality(rowmean)
        cn = None
        if status == "gain":
            cn = query_interval(gain_cn_iv, pos)
        elif status == "loss":
            cn = query_interval(loss_cn_iv, pos)
        depth = region_mean_depth(bam, chrom, pos - args.win, pos + args.win)
        status_count[status] += 1
        results.append({
            "chrom": chrom, "pos": pos, "seqc2_status": status, "seqc2_cn": cn,
            "anchor_auc": round(auc, 4),
            "shuffle_p95": round(p95, 4) if p95 is not None else None,
            "auc_minus_p95": round(auc - p95, 4) if p95 is not None else None,
            "n_hp1": n1, "n_hp2": n2, "n_cpg": ncpg,
            "mean_depth": round(depth, 1) if depth is not None else None,
            "gmm_bic_diff": round(bic_diff, 2) if bic_diff is not None else None,
            "dip_p": round(dip_p, 4) if dip_p is not None else None,
        })
    bam.close()

    out = {
        "chrom": chrom, "n_regions": len(results),
        "status_count": status_count,
        "params": {"win": args.win, "K": args.K, "max_per_status": args.max_per_status},
        "regions": results,
    }
    outpath = f"{args.out_dir}/seqc2_cn_methyl_{chrom}.json"
    json.dump(out, open(outpath, "w"), ensure_ascii=False, indent=2)
    print(f"[seqc2] {chrom}: {len(results)} regions  status={status_count}  -> {outpath}")

if __name__ == "__main__":
    main()
