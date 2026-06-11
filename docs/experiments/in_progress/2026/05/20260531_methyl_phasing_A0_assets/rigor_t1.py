#!/usr/bin/env python3
"""
T1 嚴格驗證 + 反例特徵化：unphase→H1/H2 救援在哪裡失敗？
held-out（train 量分離度、test 量 assist，無循環），每 window 記富特徵：
  train_sep / test_acc / n_cpg / n_train / n_test / mean_cov / SEQC2 status(LOH/gain/neutral)
→ 特徵化「失敗 window（test_acc<0.7）」與「不可分 window（train_sep<0.85）」的成因。
唯讀。每 chr 落 rigor_t1_{chrom}.json。
"""
import sys, json, argparse, bisect
import numpy as np
import pysam
from sklearn.metrics import roc_auc_score

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GVCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
CNVDIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
GLL = CNVDIR + "/ngs_benchmark_cnvs_gain_loss_loh.bed"

def load_bed(path):
    d = {}
    try:
        for line in open(path):
            p = line.rstrip("\n").split("\t")
            if len(p) <= 3 or not p[0].startswith("chr"): continue
            try: s, e = int(p[1]), int(p[2])
            except ValueError: continue
            d.setdefault(p[0], []).append((s, e, p[3].strip()))
    except FileNotFoundError: pass
    for c in d: d[c].sort()
    return d

def qstatus(iv, pos):
    if not iv: return "neutral"
    st = [x[0] for x in iv]; i = bisect.bisect_right(st, pos) - 1
    for j in range(max(0, i - 3), min(len(iv), i + 2)):
        s, e, v = iv[j]
        if s <= pos <= e: return v
    return "neutral"

def load_ps_blocks(chrom, min_hets=20):
    vcf = pysam.VariantFile(GVCF); by_ps = {}
    for rec in vcf.fetch(chrom):
        s0 = rec.samples[0]; gt = s0.get("GT"); ps = s0.get("PS")
        if gt is None or len(gt) != 2 or gt[0] == gt[1] or None in gt or ps is None: continue
        if len(rec.alleles) < 2 or len(rec.ref) != 1 or len(rec.alts[0]) != 1: continue
        by_ps.setdefault(ps, []).append((rec.pos - 1, rec.ref.upper(), rec.alts[0].upper(), gt))
    vcf.close()
    return {k: v for k, v in by_ps.items() if len(v) >= min_hets}

def reads_truehp(bam, chrom, hets, ws, we, min_votes=2, min_cpg=5):
    hl = {p: (r, a, g) for p, r, a, g in hets}; hp_sorted = sorted(hl)
    out, seen = [], set()
    for read in bam.fetch(chrom, max(0, ws), we):
        if read.is_secondary or read.is_supplementary or read.is_unmapped: continue
        if read.query_name in seen: continue
        qseq = read.query_sequence
        if qseq is None: continue
        pairs = read.get_aligned_pairs(matches_only=True)
        r2q = {rp: qp for qp, rp in pairs}
        votes = {1: 0, 2: 0}
        lo = bisect.bisect_left(hp_sorted, read.reference_start); hi = bisect.bisect_right(hp_sorted, read.reference_end)
        for rp in hp_sorted[lo:hi]:
            qp = r2q.get(rp)
            if qp is None: continue
            base = qseq[qp].upper(); r, a, g = hl[rp]; al = [r, a]
            h1, h2 = al[g[0]], al[g[1]]
            if base == h1 and base != h2: votes[1] += 1
            elif base == h2 and base != h1: votes[2] += 1
        tot = votes[1] + votes[2]
        if tot < min_votes: continue
        mh = 1 if votes[1] >= votes[2] else 2
        if max(votes[1], votes[2]) / tot < 0.8: continue
        mods = read.modified_bases or {}; mc = None
        for k, calls in mods.items():
            if k[0] in ("C", b"C") and k[2] in ("m", 27551): mc = calls; break
        if not mc: continue
        q2r = {q: r for q, r in pairs}; meth = {}
        for qpos, qual in mc:
            rpos = q2r.get(qpos)
            if rpos is None or rpos < ws or rpos > we: continue
            meth[rpos] = qual / 255.0
        if len(meth) < min_cpg: continue
        seen.add(read.query_name); out.append({"true_hp": mh, "meth": meth})
    return out

def cen(reads, pos):
    M = np.full((len(reads), len(pos)), np.nan)
    for i, r in enumerate(reads):
        for j, p in enumerate(pos):
            if p in r["meth"]: M[i, j] = r["meth"][p]
    return np.nanmean(M, axis=0) if len(reads) else None

def rmse(m, c, pos):
    v, cv = [], []
    for j, p in enumerate(pos):
        if p in m and not np.isnan(c[j]): v.append(m[p]); cv.append(c[j])
    return np.sqrt(np.mean((np.array(v) - np.array(cv)) ** 2)) if len(v) >= 3 else None

def loo(train, pos):
    n = len(train)
    if n < 8: return None
    sc, lb = [], []
    for i in range(n):
        rest = [train[j] for j in range(n) if j != i]
        c1 = cen([r for r in rest if r["true_hp"] == 1], pos); c2 = cen([r for r in rest if r["true_hp"] == 2], pos)
        if c1 is None or c2 is None: continue
        d1 = rmse(train[i]["meth"], c1, pos); d2 = rmse(train[i]["meth"], c2, pos)
        if d1 is None or d2 is None: continue
        sc.append(d2 - d1); lb.append(0 if train[i]["true_hp"] == 1 else 1)
    if len(set(lb)) < 2 or len(lb) < 6: return None
    try: a = roc_auc_score(lb, sc); return max(a, 1 - a)
    except Exception: return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True); ap.add_argument("--max-blocks", type=int, default=60)
    ap.add_argument("--win", type=int, default=2000); ap.add_argument("--seed", type=int, default=20260609)
    ap.add_argument("--out-dir", default="."); args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    gll = load_bed(GLL).get(args.chrom, [])
    blocks = load_ps_blocks(args.chrom)
    bam = pysam.AlignmentFile(BAM, "rb"); rows = []
    for ps in sorted(blocks, key=lambda k: -len(blocks[k])):
        if len(rows) >= args.max_blocks: break
        hets = blocks[ps]; centers = [h[0] for h in hets]; step = max(1, len(centers) // 6)
        for c in centers[::step]:
            if len(rows) >= args.max_blocks: break
            reads = reads_truehp(bam, args.chrom, hets, c - args.win, c + args.win)
            n = len(reads)
            if n < 20: continue
            idx = rng.permutation(n); h = n // 2; train = [reads[i] for i in idx[:h]]; test = [reads[i] for i in idx[h:]]
            if min(sum(r["true_hp"] == 1 for r in train), sum(r["true_hp"] == 2 for r in train),
                   sum(r["true_hp"] == 1 for r in test), sum(r["true_hp"] == 2 for r in test)) < 4: continue
            allp = set()
            for r in train + test: allp.update(r["meth"].keys())
            pos = sorted(p for p in allp if sum(1 for r in train + test if p in r["meth"]) >= 0.3 * (len(train) + len(test)))
            if len(pos) < 5: continue
            ts = loo(train, pos)
            c1 = cen([r for r in train if r["true_hp"] == 1], pos); c2 = cen([r for r in train if r["true_hp"] == 2], pos)
            if c1 is None or c2 is None: continue
            corr = tot = 0
            for r in test:
                d1 = rmse(r["meth"], c1, pos); d2 = rmse(r["meth"], c2, pos)
                if d1 is None or d2 is None: continue
                tot += 1; corr += int((1 if d1 < d2 else 2) == r["true_hp"])
            if tot < 4: continue
            rows.append({"acc": round(corr / tot, 4), "train_sep": round(float(ts), 4) if ts is not None else None,
                         "n_cpg": len(pos), "n_train": len(train), "n_test": tot, "n_reads": n,
                         "status": qstatus(gll, c + 1)})
    bam.close()
    def grp(pred):
        sub = [r for r in rows if pred(r)]
        if not sub: return {"n": 0}
        return {"n": len(sub), "acc_median": round(float(np.median([r["acc"] for r in sub])), 4),
                "n_cpg_median": int(np.median([r["n_cpg"] for r in sub])),
                "n_reads_median": int(np.median([r["n_reads"] for r in sub])),
                "train_sep_median": round(float(np.median([r["train_sep"] for r in sub if r["train_sep"] is not None])), 4) if any(r["train_sep"] is not None for r in sub) else None}
    # 失敗 vs 成功；可分 vs 不可分；按 status
    out = {"chrom": args.chrom, "n_windows": len(rows),
           "FAIL_acc<0.7": grp(lambda r: r["acc"] < 0.7),
           "OK_acc>=0.9": grp(lambda r: r["acc"] >= 0.9),
           "separable_sep>=0.85": grp(lambda r: r["train_sep"] is not None and r["train_sep"] >= 0.85),
           "notsep_sep<0.85": grp(lambda r: r["train_sep"] is not None and r["train_sep"] < 0.85),
           "by_status": {s: grp(lambda r, s=s: r["status"] == s) for s in ["neutral", "gain", "loss", "loh"]},
           "rows": rows}
    outp = f"{args.out_dir}/rigor_t1_{args.chrom}.json"
    json.dump(out, open(outp, "w"), indent=1, ensure_ascii=False)
    f = out["FAIL_acc<0.7"]; o = out["OK_acc>=0.9"]
    print(f"[rigor-t1] {args.chrom}: win={len(rows)} | FAIL(n={f['n']}) n_cpg_med={f.get('n_cpg_median')} n_reads_med={f.get('n_reads_median')} sep_med={f.get('train_sep_median')} "
          f"| OK(n={o['n']}) n_cpg_med={o.get('n_cpg_median')} sep_med={o.get('train_sep_median')} -> {outp}")

if __name__ == "__main__":
    main()
