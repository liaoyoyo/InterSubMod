#!/usr/bin/env python3
"""
聚焦驗證：在「已驗證明顯分群」的位點，甲基資訊夠不夠輔助 tag？（無循環設計）

回答用戶：「用已驗證有明顯分群的位點的 read，確認是否有足夠資訊輔助 tag 結果」。

無循環關鍵（與 V6 不同）：
  - 分離度（separability）在 **train read** 上量（LOO-AUC）→ 獨立判斷「此位點是否明顯分群」。
  - assist 正確率在 **disjoint 的 test read** 上量（test read 遮 germline 證據當 unphase，只用甲基預測）。
  - 兩者用不同 read 子集 → 「用分離位點測 assist」不循環（分離 verify 與 assist measure 不共享 read）。

每 PS-block window：
  1. 用 block het SNP 給每條 read 定真 HP（多數票，≥80% 純）= ground truth。
  2. train/test 隨機各半。
  3. train_sep = train read 的 LOO 質心 AUC（此位點甲基分群明顯度，獨立於 test）。
  4. 從 train 建 HP1/HP2 甲基質心 → 預測 test read（遮 germline）。
  5. 記 test assist accuracy + 每條 test read 信心 margin → 高信心 yield + 高信心 acc。
  6. 彙總按 train_sep 分層：verified-separable (≥0.85) vs 不可分 → 比 assist。
唯讀。每 chr 落 assist_tag_{chrom}.json。
"""
import sys, json, argparse, bisect
import numpy as np
import pysam
from sklearn.metrics import roc_auc_score

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GVCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"

def load_ps_blocks(chrom, min_hets=20):
    """回傳 {ps: [(pos0,ref,alt,gt)]}（phased het SNV，按 PS 分組，僅留 ≥min_hets 的 block）。"""
    vcf = pysam.VariantFile(GVCF)
    by_ps = {}
    for rec in vcf.fetch(chrom):
        s0 = rec.samples[0]
        gt = s0.get("GT"); ps = s0.get("PS")
        if gt is None or len(gt) != 2 or gt[0] == gt[1] or None in gt or ps is None:
            continue
        if len(rec.alleles) < 2 or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue
        by_ps.setdefault(ps, []).append((rec.pos - 1, rec.ref.upper(), rec.alts[0].upper(), gt))
    vcf.close()
    return {k: v for k, v in by_ps.items() if len(v) >= min_hets}

def read_truehp_and_meth(bam, chrom, hets, win_start, win_end, min_het_votes=2, min_cpg=5):
    het_lookup = {p: (ref, alt, gt) for p, ref, alt, gt in hets}
    het_pos_sorted = sorted(het_lookup.keys())
    out, seen = [], set()
    for read in bam.fetch(chrom, max(0, win_start), win_end):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if read.query_name in seen:
            continue
        try:
            pairs = read.get_aligned_pairs(matches_only=True)
        except Exception:
            continue
        qseq = read.query_sequence
        if qseq is None:
            continue
        r2q = {rpos: qpos for qpos, rpos in pairs}
        votes = {1: 0, 2: 0}
        lo = bisect.bisect_left(het_pos_sorted, read.reference_start)
        hi = bisect.bisect_right(het_pos_sorted, read.reference_end)
        for rp in het_pos_sorted[lo:hi]:
            qp = r2q.get(rp)
            if qp is None:
                continue
            base = qseq[qp].upper()
            ref, alt, gt = het_lookup[rp]
            alleles = [ref, alt]
            h1, h2 = alleles[gt[0]], alleles[gt[1]]
            if base == h1 and base != h2:
                votes[1] += 1
            elif base == h2 and base != h1:
                votes[2] += 1
        total = votes[1] + votes[2]
        if total < min_het_votes:
            continue
        maxhp = 1 if votes[1] >= votes[2] else 2
        if max(votes[1], votes[2]) / total < 0.8:
            continue
        mods = read.modified_bases or {}
        mc = None
        for k, calls in mods.items():
            if k[0] in ("C", b"C") and k[2] in ("m", 27551):
                mc = calls; break
        if not mc:
            continue
        q2r = {q: r for q, r in pairs}
        meth = {}
        for qpos, qual in mc:
            rpos = q2r.get(qpos)
            if rpos is None or rpos < win_start or rpos > win_end:
                continue
            meth[rpos] = qual / 255.0
        if len(meth) < min_cpg:
            continue
        seen.add(read.query_name)
        out.append({"true_hp": maxhp, "meth": meth})
    return out

def centroid(reads, positions):
    M = np.full((len(reads), len(positions)), np.nan)
    for i, r in enumerate(reads):
        for j, p in enumerate(positions):
            if p in r["meth"]:
                M[i, j] = r["meth"][p]
    return np.nanmean(M, axis=0) if len(reads) else None

def rmse(a_meth, c, positions):
    vals, cv = [], []
    for j, p in enumerate(positions):
        if p in a_meth and not np.isnan(c[j]):
            vals.append(a_meth[p]); cv.append(c[j])
    if len(vals) < 3:
        return None
    return np.sqrt(np.mean((np.array(vals) - np.array(cv)) ** 2))

def loo_train_sep(train, positions):
    """train read 的 LOO 質心 AUC（獨立分離度）。"""
    n = len(train)
    if n < 8:
        return None
    scores, labs = [], []
    for i in range(n):
        rest = [train[j] for j in range(n) if j != i]
        c1 = centroid([r for r in rest if r["true_hp"] == 1], positions)
        c2 = centroid([r for r in rest if r["true_hp"] == 2], positions)
        if c1 is None or c2 is None:
            continue
        d1 = rmse(train[i]["meth"], c1, positions); d2 = rmse(train[i]["meth"], c2, positions)
        if d1 is None or d2 is None:
            continue
        scores.append(d2 - d1); labs.append(0 if train[i]["true_hp"] == 1 else 1)
    if len(set(labs)) < 2 or len(labs) < 6:
        return None
    try:
        a = roc_auc_score(labs, scores); return max(a, 1 - a)
    except Exception:
        return None

def window_eval(reads, rng, margin_thr=0.10):
    n = len(reads)
    if n < 20:
        return None
    idx = rng.permutation(n); half = n // 2
    train = [reads[i] for i in idx[:half]]; test = [reads[i] for i in idx[half:]]
    if sum(r["true_hp"] == 1 for r in train) < 4 or sum(r["true_hp"] == 2 for r in train) < 4:
        return None
    if sum(r["true_hp"] == 1 for r in test) < 4 or sum(r["true_hp"] == 2 for r in test) < 4:
        return None
    allpos = set()
    for r in train + test:
        allpos.update(r["meth"].keys())
    positions = sorted(p for p in allpos if sum(1 for r in train + test if p in r["meth"]) >= 0.3 * (len(train) + len(test)))
    if len(positions) < 5:
        return None
    train_sep = loo_train_sep(train, positions)
    c1 = centroid([r for r in train if r["true_hp"] == 1], positions)
    c2 = centroid([r for r in train if r["true_hp"] == 2], positions)
    if c1 is None or c2 is None:
        return None
    correct = tot = conf_n = conf_correct = 0
    for r in test:
        d1 = rmse(r["meth"], c1, positions); d2 = rmse(r["meth"], c2, positions)
        if d1 is None or d2 is None:
            continue
        pred = 1 if d1 < d2 else 2
        margin = abs(d1 - d2) / (d1 + d2 + 1e-9)
        ok = (pred == r["true_hp"])
        tot += 1; correct += int(ok)
        if margin >= margin_thr:
            conf_n += 1; conf_correct += int(ok)
    if tot < 4:
        return None
    return {
        "train_sep": round(float(train_sep), 4) if train_sep is not None else None,
        "test_acc": round(correct / tot, 4),
        "n_train": len(train), "n_test": tot,
        "conf_yield": round(conf_n / tot, 4),
        "conf_acc": round(conf_correct / conf_n, 4) if conf_n else None,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--max-blocks", type=int, default=60)
    ap.add_argument("--win", type=int, default=2000)
    ap.add_argument("--blocks-per-ps", type=int, default=6)
    ap.add_argument("--seed", type=int, default=20260609)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    blocks = load_ps_blocks(args.chrom)
    bam = pysam.AlignmentFile(BAM, "rb")
    results = []
    # 取最大的幾個 PS block，每個 block 內取數個 het 中心窗
    for ps in sorted(blocks, key=lambda k: -len(blocks[k])):
        if len(results) >= args.max_blocks:
            break
        hets = blocks[ps]
        centers = [h[0] for h in hets]
        step = max(1, len(centers) // args.blocks_per_ps)
        for c in centers[::step]:
            if len(results) >= args.max_blocks:
                break
            reads = read_truehp_and_meth(bam, args.chrom, hets, c - args.win, c + args.win)
            ev = window_eval(reads, rng)
            if ev:
                ev["ps"] = ps; ev["center"] = c + 1
                results.append(ev)
    bam.close()
    sep = [r for r in results if r["train_sep"] is not None and r["train_sep"] >= 0.85]
    nsep = [r for r in results if r["train_sep"] is not None and r["train_sep"] < 0.85]
    def med(rs, k):
        xs = [r[k] for r in rs if r.get(k) is not None]
        return round(float(np.median(xs)), 4) if xs else None
    out = {
        "chrom": args.chrom, "n_windows": len(results),
        "verified_separable_train_sep>=0.85": {
            "n": len(sep), "test_acc_median": med(sep, "test_acc"),
            "conf_yield_median": med(sep, "conf_yield"), "conf_acc_median": med(sep, "conf_acc"),
        },
        "not_separable_train_sep<0.85": {
            "n": len(nsep), "test_acc_median": med(nsep, "test_acc"),
            "conf_yield_median": med(nsep, "conf_yield"), "conf_acc_median": med(nsep, "conf_acc"),
        },
        "overall_test_acc_median": med(results, "test_acc"),
        "windows": results,
    }
    outp = f"{args.out_dir}/assist_tag_{args.chrom}.json"
    json.dump(out, open(outp, "w"), indent=1, ensure_ascii=False)
    s = out["verified_separable_train_sep>=0.85"]; ns = out["not_separable_train_sep<0.85"]
    print(f"[assist] {args.chrom}: {len(results)} win | SEPARABLE(n={s['n']}) test_acc={s['test_acc_median']} conf_yield={s['conf_yield_median']} conf_acc={s['conf_acc_median']} | NOT-SEP(n={ns['n']}) test_acc={ns['test_acc_median']} -> {outp}")

if __name__ == "__main__":
    main()
