#!/usr/bin/env python3
"""
外推驗證（補最後一哩）：用 PS phase-block 當獨立 ground truth，直接驗
「甲基能否救『假裝成 unphase』的 read」。

前序問題：所有 anchor AUC 都用「有 germline 證據的 read」自證，unphase read
本身無 ground truth，「甲基能分 anchor」→「能救 unphase」是外推。

本腳本的解法（held-out simulation，模擬真 unphase 救援）：
  1. 在一個大 PS block 內，用「該 block 的 germline het SNP」給每條 read 定真 HP（ground truth）。
  2. 把 read 隨機分 train/test 兩半。
  3. test read：**完全遮住其 germline 證據**（當作 unphase），只給甲基。
  4. 用 train read 的甲基（含真 HP 標籤）建 HP1/HP2 甲基 profile（質心）。
  5. test read 用甲基距離預測 HP，對照 step1 的真 HP → 算救援正確率 (accuracy) + AUC。
  6. 對照：shuffle train 標籤 → null accuracy（應掉到 ~0.5）。

這比 anchor leave-one-out 更接近真實救援：predictor 從『其他 read』學，
被預測 read 自己『沒有 germline 證據』，正是 unphase read 的處境。

唯讀。每 chr 落 extrapolation_{chrom}.json。
"""
import sys, json, argparse, bisect
import numpy as np
import pysam
from sklearn.metrics import roc_auc_score

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GVCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
AD = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets"

def load_phased_hets(chrom):
    """回傳 [(pos0, ref, alt, hap)] hap: 0|1→alt在H2, 1|0→alt在H1。只取單一最大 PS block。"""
    vcf = pysam.VariantFile(GVCF)
    by_ps = {}
    for rec in vcf.fetch(chrom):
        s0 = rec.samples[0]
        gt = s0.get("GT"); ps = s0.get("PS")
        if gt is None or len(gt) != 2 or gt[0] == gt[1] or None in gt or ps is None:
            continue
        if len(rec.alleles) < 2: continue
        ref, alt = rec.ref, rec.alts[0]
        if len(ref) != 1 or len(alt) != 1: continue  # 只 SNP
        # phased: gt=(0,1) 表 H1=ref,H2=alt; (1,0) 表 H1=alt,H2=ref
        by_ps.setdefault(ps, []).append((rec.pos - 1, ref, alt, gt))
    if not by_ps:
        return None, None
    best_ps = max(by_ps, key=lambda k: len(by_ps[k]))
    return best_ps, by_ps[best_ps]

def read_truehp_and_meth(bam, chrom, hets, win_start, win_end, rng, min_het_votes=3, min_cpg=5):
    """對每條 read：用 hets 定真 HP（多數票），收 per-CpG 甲基。
    回傳 list of dict: read_id, true_hp(1/2), het_positions_used(set), meth{rpos:prob}, hp_votes."""
    het_lookup = {p: (ref, alt, gt) for p, ref, alt, gt in hets}
    het_pos_sorted = sorted(het_lookup.keys())
    out = []
    seen = set()
    for read in bam.fetch(chrom, win_start, win_end):
        if read.is_secondary or read.is_supplementary or read.is_unmapped: continue
        if read.query_name in seen: continue
        # read 覆蓋的 het 位點 → 看 read 在該位是 ref 還 alt → 投 H1/H2
        try:
            pairs = read.get_aligned_pairs(matches_only=True, with_seq=False)
        except Exception:
            continue
        qseq = read.query_sequence
        if qseq is None: continue
        r2q = {}
        for qpos, rpos in pairs:
            r2q[rpos] = qpos
        votes = {1: 0, 2: 0}
        used = set()
        lo = bisect.bisect_left(het_pos_sorted, read.reference_start)
        hi = bisect.bisect_right(het_pos_sorted, read.reference_end)
        for rp in het_pos_sorted[lo:hi]:
            qp = r2q.get(rp)
            if qp is None: continue
            base = qseq[qp].upper()
            ref, alt, gt = het_lookup[rp]
            # gt (a,b): H1 allele = alleles[a], H2 allele = alleles[b]; alleles=[ref,alt]
            alleles = [ref, alt]
            h1_allele = alleles[gt[0]]; h2_allele = alleles[gt[1]]
            if base == h1_allele and base != h2_allele:
                votes[1] += 1; used.add(rp)
            elif base == h2_allele and base != h1_allele:
                votes[2] += 1; used.add(rp)
        total = votes[1] + votes[2]
        if total < min_het_votes: continue
        # 真 HP = 多數，且要夠純（>=0.8 一致）
        maxhp = 1 if votes[1] >= votes[2] else 2
        if max(votes[1], votes[2]) / total < 0.8: continue
        # 甲基
        mods = read.modified_bases or {}
        mc = None
        for k, calls in mods.items():
            if k[0] in ("C", b"C") and k[2] in ("m", 27551):
                mc = calls; break
        if not mc: continue
        q2r = {q: r for q, r in pairs}
        meth = {}
        for qpos, qual in mc:
            rpos = q2r.get(qpos)
            if rpos is None: continue
            meth[rpos] = qual / 255.0
        if len(meth) < min_cpg: continue
        seen.add(read.query_name)
        out.append({"read_id": read.query_name, "true_hp": maxhp,
                    "votes": dict(votes), "meth": meth})
    return out

def predict_acc(reads, K_shuffle, rng):
    """train/test split → 甲基質心預測 test 真 HP → accuracy + AUC + shuffle null。"""
    n = len(reads)
    if n < 20: return None
    idx = rng.permutation(n)
    half = n // 2
    train = [reads[i] for i in idx[:half]]
    test = [reads[i] for i in idx[half:]]
    if sum(r["true_hp"] == 1 for r in train) < 4 or sum(r["true_hp"] == 2 for r in train) < 4:
        return None
    if sum(r["true_hp"] == 1 for r in test) < 4 or sum(r["true_hp"] == 2 for r in test) < 4:
        return None
    # 共同 CpG 欄
    allpos = set()
    for r in train + test: allpos.update(r["meth"].keys())
    positions = sorted(allpos)
    cov = {p: sum(1 for r in train+test if p in r["meth"]) for p in positions}
    positions = [p for p in positions if cov[p] >= 0.3 * (len(train)+len(test))]
    if len(positions) < 5: return None
    def vec(r):
        return np.array([r["meth"].get(p, np.nan) for p in positions])
    # 質心（train，忽略 NaN）
    def centroid(group):
        mat = np.array([vec(r) for r in group])
        return np.nanmean(mat, axis=0)
    def acc_for_labels(train_lab):
        c1 = centroid([train[i] for i in range(len(train)) if train_lab[i] == 1])
        c2 = centroid([train[i] for i in range(len(train)) if train_lab[i] == 2])
        correct = 0; total = 0; scores = []; ys = []
        for r in test:
            v = vec(r)
            mask = ~np.isnan(v) & ~np.isnan(c1) & ~np.isnan(c2)
            if mask.sum() < 3: continue
            d1 = np.sqrt(np.mean((v[mask]-c1[mask])**2))
            d2 = np.sqrt(np.mean((v[mask]-c2[mask])**2))
            pred = 1 if d1 < d2 else 2
            correct += (pred == r["true_hp"]); total += 1
            scores.append(d2 - d1)  # 越大越像 H1
            ys.append(1 if r["true_hp"] == 1 else 0)
        if total < 8 or len(set(ys)) < 2: return None, None, total
        acc = correct / total
        try: auc = roc_auc_score(ys, scores); auc = max(auc, 1-auc)
        except Exception: auc = None
        return acc, auc, total
    real_lab = [r["true_hp"] for r in train]
    acc, auc, ntest = acc_for_labels(real_lab)
    if acc is None: return None
    # shuffle null
    null_accs = []
    for _ in range(K_shuffle):
        sl = list(rng.permutation(real_lab))
        if len(set(sl)) < 2: continue
        a, _, _ = acc_for_labels(sl)
        if a is not None: null_accs.append(a)
    return {"accuracy": round(acc, 4), "auc": round(auc, 4) if auc else None,
            "n_test": ntest, "n_train": len(train),
            "null_acc_mean": round(float(np.mean(null_accs)), 4) if null_accs else None,
            "null_acc_p95": round(float(np.percentile(null_accs, 95)), 4) if null_accs else None,
            "acc_minus_null_p95": round(acc - np.percentile(null_accs, 95), 4) if null_accs else None}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--n-blocks", type=int, default=15, help="取幾個大 PS block")
    ap.add_argument("--block-span", type=int, default=2000, help="每 block 取中心 ±span 的 read（須夠小讓 read 共享 CpG）")
    ap.add_argument("--windows-per-block", type=int, default=3, help="每 block 取幾個局部窗測試")
    ap.add_argument("--K", type=int, default=30)
    args = ap.parse_args()
    chrom = args.chrom
    rng = np.random.RandomState(20260602)

    # 取多個大 PS block（不只最大一個）
    vcf = pysam.VariantFile(GVCF)
    by_ps = {}
    for rec in vcf.fetch(chrom):
        s0 = rec.samples[0]; gt = s0.get("GT"); ps = s0.get("PS")
        if gt is None or len(gt) != 2 or gt[0] == gt[1] or None in gt or ps is None: continue
        if len(rec.ref) != 1 or len(rec.alts[0]) != 1: continue
        by_ps.setdefault(ps, []).append((rec.pos-1, rec.ref, rec.alts[0], gt))
    big_blocks = sorted(by_ps.items(), key=lambda kv: -len(kv[1]))[:args.n_blocks]

    bam = pysam.AlignmentFile(BAM, "rb")
    results = []
    for ps, hets in big_blocks:
        positions = sorted(h[0] for h in hets)
        # 在 block 內取多個局部窗（±block_span 小窗，read 才共享 CpG）
        lo, hi = positions[0], positions[-1]
        if hi - lo < 2 * args.block_span:
            centers = [int(np.median(positions))]
        else:
            centers = [int(c) for c in np.linspace(lo + args.block_span, hi - args.block_span,
                                                   args.windows_per_block)]
        for center in centers:
            ws, we = center - args.block_span, center + args.block_span
            reads = read_truehp_and_meth(bam, chrom, hets, ws, we, rng)
            n1 = sum(r["true_hp"] == 1 for r in reads); n2 = sum(r["true_hp"] == 2 for r in reads)
            res = predict_acc(reads, args.K, rng) if len(reads) >= 20 else None
            if res:
                res.update({"ps_block": ps, "center": center, "n_reads": len(reads),
                            "n_hp1": n1, "n_hp2": n2, "n_block_hets": len(hets)})
                results.append(res)
    bam.close()

    accs = [r["accuracy"] for r in results]
    deltas = [r["acc_minus_null_p95"] for r in results if r["acc_minus_null_p95"] is not None]
    out = {"chrom": chrom, "n_blocks_tested": len(results),
           "accuracy_median": round(float(np.median(accs)), 4) if accs else None,
           "accuracy_mean": round(float(np.mean(accs)), 4) if accs else None,
           "null_acc_median": round(float(np.median([r["null_acc_mean"] for r in results if r["null_acc_mean"]])), 4) if results else None,
           "frac_acc_above_null_p95": round(float(np.mean([d > 0 for d in deltas])), 3) if deltas else None,
           "blocks": results}
    outpath = f"{AD}/extrapolation_{chrom}.json"
    json.dump(out, open(outpath, "w"), ensure_ascii=False, indent=2)
    print(f"[extrap] {chrom}: {len(results)} blocks, acc_median={out['accuracy_median']}, "
          f"null_median={out['null_acc_median']}, frac>null_p95={out['frac_acc_above_null_p95']} -> {outpath}")

if __name__ == "__main__":
    main()
