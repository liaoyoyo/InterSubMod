#!/usr/bin/env python3
"""
Shuffle-label 控制：判定 germline-het null AUC median 0.974 是真 label-dependent 訊號，
還是 leave-one-out 質心 AUC 方法本身的 leak / 對稱化偏差。

設計（與 germline_het_null.py 同一批 region，seed 20260601，可重現）：
  對每個 null 區的同一距離矩陣 D：
    - real_auc      = loo_centroid_auc(D, 真 HP 標籤)   [對稱化]
    - shuffle_auc   = K 次 loo_centroid_auc(D, 打亂標籤) [對稱化] 的 mean / p95
  判讀：
    - real >> shuffle  → 0.974 是真 label-dependent（甲基隨 HP 變）→ 對 rescue 是好消息
    - real ≈ shuffle   → 方法 leak（AUC 與真標籤無關仍高）→ 前面所有 AUC 打折

額外診斷：同時報 raw（未對稱化）real AUC，量化對稱化 max(auc,1-auc) 的上偏程度。

⚠ 注意：shuffle 控制能抓「與真標籤無關的方法 leak」，但**抓不到 CpG-SNP pseudo-ASM**
（那是與真標籤相關的 artifact，shuffle 後也會掉 → 看似真訊號）。CpG-SNP 排除是另一支控制（下一步）。

唯讀。結果寫 shuffle_control_results.json。
"""
import sys, json, random
import numpy as np
import pysam
from sklearn.metrics import roc_auc_score

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GERMLINE_VCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/shuffle_control_results.json"

EXCLUDE = [
    ("chr13", 32_310_000, 32_320_000),
    ("chr11", 1_995_000, 2_006_000),
    ("chr15", 24_950_000, 24_960_000),
    ("chr20", 58_885_000, 58_915_000),
]
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]

def in_exclude(chrom, pos):
    return any(chrom == c and a <= pos <= b for c, a, b in EXCLUDE)

def get_hp(read):
    try:
        return str(read.get_tag("HP"))
    except KeyError:
        return "unphase"

def region_D_lab(bam, chrom, start, end, min_anchor=10, min_cpg=5):
    """回傳 (D, lab, n_hp1, n_hp2, n_cpg) 或 None。與 germline_het_null.py 同邏輯。"""
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
    for i, (h, m) in enumerate(rows):
        lab[i] = 0 if h == "1" else 1
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
    return D, lab, int((lab == 0).sum()), int((lab == 1).sum()), len(positions)

def loo_centroid_auc(D, lab, symmetrize=True):
    n = len(lab)
    scores = np.full(n, np.nan)
    for i in range(n):
        o1 = [k for k in range(n) if k != i and lab[k] == 0]
        o2 = [k for k in range(n) if k != i and lab[k] == 1]
        if not o1 or not o2:
            continue
        scores[i] = np.mean([D[i, k] for k in o1]) - np.mean([D[i, k] for k in o2])
    v = ~np.isnan(scores)
    if v.sum() < 4 or len(set(lab[v])) < 2:
        return None
    try:
        auc = roc_auc_score(lab[v], scores[v])
    except Exception:
        return None
    return max(auc, 1 - auc) if symmetrize else auc

def main():
    n_target = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    seed = 20260601
    rng = random.Random(seed)
    nprng = np.random.RandomState(seed)
    bam = pysam.AlignmentFile(BAM, "rb")
    vcf = pysam.VariantFile(GERMLINE_VCF)

    print("[shuffle] 掃 het SNP（同 null 選法）...", file=sys.stderr)
    candidates, chr_counts = [], {c: 0 for c in AUTOSOMES}
    for rec in vcf.fetch():
        c = rec.chrom
        if c not in AUTOSOMES or chr_counts[c] >= 2000:
            continue
        gt = None
        for s in rec.samples.values():
            gt = s.get("GT"); break
        if gt is None or len(gt) != 2 or gt[0] == gt[1] or None in gt:
            continue
        if in_exclude(c, rec.pos):
            continue
        candidates.append((c, rec.pos)); chr_counts[c] += 1
    rng.shuffle(candidates)
    print(f"[shuffle] 候選={len(candidates)}", file=sys.stderr)

    results, tried = [], 0
    for c, pos in candidates:
        if len(results) >= n_target:
            break
        tried += 1
        if tried > n_target * 30:
            break
        rd = region_D_lab(bam, c, pos - 2000, pos + 2000)
        if rd is None:
            continue
        D, lab, n1, n2, ncpg = rd
        real = loo_centroid_auc(D, lab, symmetrize=True)
        real_raw = loo_centroid_auc(D, lab, symmetrize=False)
        if real is None:
            continue
        sh = []
        for _ in range(K):
            sl = nprng.permutation(lab)
            if len(set(sl)) < 2:
                continue
            a = loo_centroid_auc(D, sl, symmetrize=True)
            if a is not None:
                sh.append(a)
        if not sh:
            continue
        results.append({
            "region": f"{c}:{pos-2000}-{pos+2000}",
            "n_hp1": n1, "n_hp2": n2, "n_cpg": ncpg,
            "real_auc_sym": float(real),
            "real_auc_raw": float(real_raw) if real_raw is not None else None,
            "shuffle_auc_sym_mean": float(np.mean(sh)),
            "shuffle_auc_sym_p95": float(np.percentile(sh, 95)),
            "real_minus_shuffle_p95": float(real - np.percentile(sh, 95)),
        })
        print(f"[shuffle] {len(results)}/{n_target} {c}:{pos} "
              f"real={real:.3f} raw={real_raw:.3f} shuf_mean={np.mean(sh):.3f} "
              f"shuf_p95={np.percentile(sh,95):.3f} Δ={real-np.percentile(sh,95):+.3f}",
              file=sys.stderr)
    bam.close()

    reals = [r["real_auc_sym"] for r in results]
    reals_raw = [r["real_auc_raw"] for r in results if r["real_auc_raw"] is not None]
    shuf_means = [r["shuffle_auc_sym_mean"] for r in results]
    deltas = [r["real_minus_shuffle_p95"] for r in results]
    summary = {
        "n_regions": len(results), "n_tried": tried, "K_shuffles": K, "seed": seed,
        "real_auc_sym": {"median": float(np.median(reals)), "mean": float(np.mean(reals))} if reals else None,
        "real_auc_raw": {"median": float(np.median(reals_raw)), "mean": float(np.mean(reals_raw)),
                          "min": float(np.min(reals_raw)), "max": float(np.max(reals_raw))} if reals_raw else None,
        "shuffle_auc_sym": {"median": float(np.median(shuf_means)), "mean": float(np.mean(shuf_means)),
                             "max": float(np.max(shuf_means))} if shuf_means else None,
        "real_minus_shuffle_p95": {"median": float(np.median(deltas)), "mean": float(np.mean(deltas)),
                                    "min": float(np.min(deltas)), "frac_positive": float(np.mean([d > 0 for d in deltas]))} if deltas else None,
        "regions": results,
    }
    json.dump(summary, open(OUT, "w"), ensure_ascii=False, indent=2)
    print(json.dumps({k: v for k, v in summary.items() if k != "regions"}, ensure_ascii=False, indent=2))
    print(f"[shuffle] DONE -> {OUT}")

if __name__ == "__main__":
    main()
