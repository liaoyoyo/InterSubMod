#!/usr/bin/env python3
"""
germline-het null：BRCA2 anchor AUC=0.866 是真訊號還是 ASM 放大假象？

設計：在隨機 germline het 區（非已知 ASM / 非 imprinting）用「與 BRCA2 完全相同」
的方法算 HP1-vs-HP2 甲基 anchor AUC 分布。

判讀：
  - null AUC 分布中位數 ~0.8+ → 甲基普遍跟 HP 相關（救 unphase 好消息）
  - null AUC 分布中位數 ~0.5-0.6 → BRCA2 0.866 是特異強 ASM，不可外推

方法（與 heatmap_and_separation.py 一致）：
  - 每個 region：抽 primary read 的 per-CpG 5mC + HP tag
  - 只用 HP1/HP2 anchor read（germline-phased，有真 HP）
  - leave-one-out 質心：每 read 對 HP1 群均距 vs HP2 群均距 → score → ROC AUC
  - pairwise-complete 距離（只比共同覆蓋 CpG）

NULL region 選擇：
  - 從 germline_phased_merged.vcf.gz 取 het SNP
  - 隨機抽 N 個，各取 ±2kb 窗
  - 排除：BRCA2 區、已知 imprinting DMR、chrX/chrY（單套）
  - 要求每窗 HP1>=10 且 HP2>=10 anchor read（否則跳過）

唯讀；全部結果寫 germline_het_null_results.json。
"""
import sys, json, random
import numpy as np
import pysam
from scipy.spatial.distance import pdist
from sklearn.metrics import roc_auc_score, silhouette_score

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GERMLINE_VCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/germline_het_null_results.json"

# 排除區（hg38）：BRCA2 + 已知 imprinting DMR
EXCLUDE = [
    ("chr13", 32_310_000, 32_320_000),   # BRCA2
    ("chr11", 1_995_000, 2_006_000),     # H19/IGF2
    ("chr15", 24_950_000, 24_960_000),   # SNRPN
    ("chr20", 58_885_000, 58_915_000),   # GNAS
]
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]

def in_exclude(chrom, pos):
    for c, a, b in EXCLUDE:
        if chrom == c and a <= pos <= b:
            return True
    return False

def get_hp(read):
    try:
        return str(read.get_tag("HP"))
    except KeyError:
        return "unphase"

def region_anchor_auc(bam, chrom, start, end, min_anchor=10, min_cpg=5):
    """回傳 (auc, silhouette, n_hp1, n_hp2, n_cpg) 或 None。"""
    rows = []  # (hp, {refpos:prob})
    seen = set()
    allpos = set()
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
    # 過濾稀疏欄（>=30% 覆蓋）
    cov = {p: sum(1 for _, m in rows if p in m) for p in positions}
    positions = [p for p in positions if cov[p] >= 0.3 * len(rows)]
    if len(positions) < min_cpg:
        return None
    M = np.full((len(rows), len(positions)), np.nan)
    lab = np.zeros(len(rows), dtype=int)
    for i, (h, m) in enumerate(rows):
        lab[i] = 0 if h == "1" else 1
        for j, p in enumerate(positions):
            if p in m:
                M[i, j] = m[p]
    # 過濾 <min_cpg 覆蓋 read
    keep = np.array([np.sum(~np.isnan(M[i])) >= min_cpg for i in range(M.shape[0])])
    M, lab = M[keep], lab[keep]
    if (lab == 0).sum() < min_anchor or (lab == 1).sum() < min_anchor:
        return None
    # pairwise-complete dist
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            mask = ~np.isnan(M[i]) & ~np.isnan(M[j])
            d = np.sqrt(np.mean((M[i, mask] - M[j, mask]) ** 2)) if mask.sum() >= 3 else np.nan
            D[i, j] = D[j, i] = d
    med = np.nanmedian(D[D > 0]) if np.any(D > 0) else 1.0
    D[np.isnan(D)] = med
    # leave-one-out anchor AUC
    scores = []
    for i in range(n):
        o1 = [k for k in range(n) if k != i and lab[k] == 0]
        o2 = [k for k in range(n) if k != i and lab[k] == 1]
        if not o1 or not o2:
            scores.append(np.nan); continue
        scores.append(np.mean([D[i, k] for k in o1]) - np.mean([D[i, k] for k in o2]))
    scores = np.array(scores)
    v = ~np.isnan(scores)
    if v.sum() < 2 * min_anchor or len(set(lab[v])) < 2:
        return None
    try:
        auc = roc_auc_score(lab[v], scores[v]); auc = max(auc, 1 - auc)
    except Exception:
        return None
    try:
        sil = float(silhouette_score(D, lab, metric="precomputed"))
    except Exception:
        sil = None
    return dict(auc=float(auc), silhouette=sil, n_hp1=int((lab == 0).sum()),
               n_hp2=int((lab == 1).sum()), n_cpg=len(positions),
               region=f"{chrom}:{start}-{end}")

def main():
    n_target = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    seed = 20260601
    rng = random.Random(seed)
    bam = pysam.AlignmentFile(BAM, "rb")
    vcf = pysam.VariantFile(GERMLINE_VCF)

    # 收集 autosome het SNP（GT 0|1 或 1|0），分散抽樣
    print("[null] 掃 germline het SNP（每染色體抽樣）...", file=sys.stderr)
    candidates = []
    per_chr_cap = 2000
    chr_counts = {c: 0 for c in AUTOSOMES}
    for rec in vcf.fetch():
        c = rec.chrom
        if c not in AUTOSOMES:
            continue
        if chr_counts[c] >= per_chr_cap:
            continue
        # het 判定
        gt = None
        for s in rec.samples.values():
            gt = s.get("GT"); break
        if gt is None or len(gt) != 2 or gt[0] == gt[1] or None in gt:
            continue
        if in_exclude(c, rec.pos):
            continue
        candidates.append((c, rec.pos))
        chr_counts[c] += 1
    print(f"[null] 候選 het SNP={len(candidates)}", file=sys.stderr)
    rng.shuffle(candidates)

    results = []
    tried = 0
    for c, pos in candidates:
        if len(results) >= n_target:
            break
        tried += 1
        if tried > n_target * 30:
            break
        r = region_anchor_auc(bam, c, pos - 2000, pos + 2000)
        if r is not None:
            results.append(r)
            print(f"[null] {len(results)}/{n_target} {r['region']} AUC={r['auc']:.3f} "
                  f"HP1={r['n_hp1']} HP2={r['n_hp2']} cpg={r['n_cpg']}", file=sys.stderr)
    bam.close()

    aucs = [r["auc"] for r in results]
    summary = {
        "n_null_regions": len(results),
        "n_tried": tried,
        "seed": seed,
        "auc_distribution": {
            "min": float(np.min(aucs)) if aucs else None,
            "p25": float(np.percentile(aucs, 25)) if aucs else None,
            "median": float(np.median(aucs)) if aucs else None,
            "mean": float(np.mean(aucs)) if aucs else None,
            "p75": float(np.percentile(aucs, 75)) if aucs else None,
            "max": float(np.max(aucs)) if aucs else None,
        },
        "frac_above_0.58": float(np.mean([a > 0.58 for a in aucs])) if aucs else None,
        "frac_above_0.70": float(np.mean([a > 0.70 for a in aucs])) if aucs else None,
        "frac_above_0.80": float(np.mean([a > 0.80 for a in aucs])) if aucs else None,
        "brca2_reference_auc": 0.866,
        "brca2_percentile_in_null": (
            float(np.mean([a < 0.866 for a in aucs])) if aucs else None),
        "regions": results,
    }
    json.dump(summary, open(OUT, "w"), ensure_ascii=False, indent=2)
    print(json.dumps({k: v for k, v in summary.items() if k != "regions"},
                     ensure_ascii=False, indent=2))
    print(f"[null] DONE -> {OUT}")

if __name__ == "__main__":
    main()
