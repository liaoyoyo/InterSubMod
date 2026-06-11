#!/usr/bin/env python3
"""
T3 局部 allele 版（回應用戶）：用「該 somatic 位點本身的 ALT/REF」當亞群/母本標籤，
比全域 HP tag(1-1 vs 1) 更乾淨。測兩件事：

GATE（必要條件）：蓋到位點的 HP1 reads，ALT(亞群 A1) vs REF(母本 B1) 甲基能不能分？
  held-out 分類 AUC + raw Δβ + 同窗噪音地板(隨機對半 B1)。若分不開 → 下面的「偏向」=雜訊。
USER Q（偏向）：靠近位點但【沒蓋到】該位點的 plain HP1 read(Q1)，其甲基距離是否更偏向 A1(亞群質心)？
  報 frac(closer to A1) vs 50%(null)。只有 GATE 通過才有意義。
HP2/A2/B2/Q2 同樣統計。唯讀。每 chr 落 t3_local_allele_{chrom}.json。
"""
import sys, json, argparse
import numpy as np
import pysam
from sklearn.metrics import roc_auc_score

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
SOMATIC_VCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/somatic_pass.vcf.gz"

def get_hp(read):
    try: return str(read.get_tag("HP"))
    except KeyError: return "unphase"

def read_info(read, pos0, ws, we, ref, alt, min_cpg=5):
    """回 (covers_site, allele('ALT'/'REF'/None), meth dict) 或 None。"""
    qseq = read.query_sequence
    if qseq is None: return None
    pairs = read.get_aligned_pairs(matches_only=True)
    r2q = {}; allele = None
    for qpos, rpos in pairs:
        r2q[rpos] = qpos
    covers = read.reference_start <= pos0 <= read.reference_end
    if pos0 in r2q:
        base = qseq[r2q[pos0]].upper()
        allele = "ALT" if base == alt else ("REF" if base == ref else None)
    mods = read.modified_bases or {}; mc = None
    for k, calls in mods.items():
        if k[0] in ("C", b"C") and k[2] in ("m", 27551): mc = calls; break
    if not mc: return None
    q2r = {q: r for q, r in pairs}; meth = {}
    for qpos, qual in mc:
        rpos = q2r.get(qpos)
        if rpos is None or rpos < ws or rpos > we: continue
        meth[rpos] = qual / 255.0
    if len(meth) < min_cpg: return None
    return (covers, allele, meth)

def centroid(group, pos):
    M = np.full((len(group), len(pos)), np.nan)
    for i, m in enumerate(group):
        for j, p in enumerate(pos):
            if p in m: M[i, j] = m[p]
    return np.nanmean(M, axis=0) if len(group) else None

def rmse(m, c, pos):
    v, cv = [], []
    for j, p in enumerate(pos):
        if p in m and not np.isnan(c[j]): v.append(m[p]); cv.append(c[j])
    return np.sqrt(np.mean((np.array(v) - np.array(cv)) ** 2)) if len(v) >= 3 else None

def heldout_auc(A, B, rng, min_anchor):
    """A(亞群) vs B(母本) held-out 分類 AUC（train 建質心、test 評）。"""
    if len(A) < min_anchor or len(B) < min_anchor: return None
    def split(g):
        idx = rng.permutation(len(g)); h = len(g) // 2; return [g[i] for i in idx[:h]], [g[i] for i in idx[h:]]
    trA, teA = split(A); trB, teB = split(B)
    if min(len(trA), len(teA), len(trB), len(teB)) < 3: return None
    allp = set()
    for m in trA + trB: allp.update(m.keys())
    pos = [p for p in allp if sum(1 for m in trA + trB if p in m) >= 0.3 * (len(trA) + len(trB))]
    if len(pos) < 5: return None
    cA, cB = centroid(trA, pos), centroid(trB, pos)
    if cA is None or cB is None: return None
    sc, lb = [], []
    for m in teA:
        dA, dB = rmse(m, cA, pos), rmse(m, cB, pos)
        if dA is not None and dB is not None: sc.append(dB - dA); lb.append(1)
    for m in teB:
        dA, dB = rmse(m, cA, pos), rmse(m, cB, pos)
        if dA is not None and dB is not None: sc.append(dB - dA); lb.append(0)
    if len(set(lb)) < 2 or len(lb) < 6: return None
    try: return roc_auc_score(lb, sc)
    except Exception: return None

def lean(Q, A, B, min_anchor):
    """ambiguous Q reads：frac 更接近 A(亞群質心) than B(母本質心)。"""
    if len(A) < min_anchor or len(B) < min_anchor or not Q: return None, 0
    allp = set()
    for m in A + B: allp.update(m.keys())
    pos = [p for p in allp if sum(1 for m in A + B if p in m) >= 0.3 * (len(A) + len(B))]
    if len(pos) < 5: return None, 0
    cA, cB = centroid(A, pos), centroid(B, pos)
    closerA = tot = 0
    for m in Q:
        dA, dB = rmse(m, cA, pos), rmse(m, cB, pos)
        if dA is None or dB is None: continue
        tot += 1; closerA += int(dA < dB)
    return (closerA / tot if tot else None), tot

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True); ap.add_argument("--max-sites", type=int, default=250)
    ap.add_argument("--win", type=int, default=4000); ap.add_argument("--min-anchor", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260611); ap.add_argument("--excl-dist", type=int, default=100); ap.add_argument("--out-dir", default="."); args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    vcf = pysam.VariantFile(SOMATIC_VCF)
    sites = [(rec.pos - 1, rec.ref.upper(), rec.alts[0].upper()) for rec in vcf.fetch(args.chrom)
             if rec.alts and len(rec.ref) == 1 and len(rec.alts[0]) == 1]
    vcf.close()
    step = max(1, len(sites) // (args.max_sites * 3)) if sites else 1
    bam = pysam.AlignmentFile(BAM, "rb"); rows = []
    for pos0, ref, alt in sites[::step]:
        if len(rows) >= args.max_sites: break
        ws, we = pos0 - args.win, pos0 + args.win
        A1, B1, A2, B2, Q1, Q2 = [], [], [], [], [], []
        seen = set()
        for read in bam.fetch(args.chrom, max(0, ws), we):
            if read.is_secondary or read.is_supplementary or read.is_unmapped: continue
            if read.query_name in seen: continue
            hp = get_hp(read)
            if hp not in ("1", "2"): continue
            info = read_info(read, pos0, ws, we, ref, alt)
            if info is None: continue
            seen.add(read.query_name)
            covers, allele, meth = info
            if covers and allele == "ALT": (A1 if hp == "1" else A2).append(meth)
            elif covers and allele == "REF": (B1 if hp == "1" else B2).append(meth)
            elif not covers: (Q1 if hp == "1" else Q2).append(meth)
        rec = {"site": pos0 + 1, "nA1": len(A1), "nB1": len(B1), "nA2": len(A2), "nB2": len(B2), "nQ1": len(Q1), "nQ2": len(Q2)}
        a = heldout_auc(A1, B1, rng, args.min_anchor)
        if a is not None: rec["gate_auc_h1"] = round(float(a), 4)
        a = heldout_auc(A2, B2, rng, args.min_anchor)
        if a is not None: rec["gate_auc_h2"] = round(float(a), 4)
        # CONTROL：排除 somatic 位點 ±D bp 的 CpG（測 GATE 是否為 somatic-CpG artifact）
        D = args.excl_dist
        far = lambda g: [{p: v for p, v in m.items() if abs(p - pos0) > D} for m in g]
        A1f, B1f, A2f, B2f = far(A1), far(B1), far(A2), far(B2)
        a = heldout_auc([m for m in A1f if len(m) >= 5], [m for m in B1f if len(m) >= 5], rng, args.min_anchor)
        if a is not None: rec["gate_auc_h1_farCpG"] = round(float(a), 4)
        a = heldout_auc([m for m in A2f if len(m) >= 5], [m for m in B2f if len(m) >= 5], rng, args.min_anchor)
        if a is not None: rec["gate_auc_h2_farCpG"] = round(float(a), 4)
        # 噪音地板：同群(母本 B)隨機對半當假 A/B → held-out AUC（純抽樣噪音）
        if len(B1) >= 2 * args.min_anchor:
            idx = rng.permutation(len(B1)); h = len(B1) // 2
            nf = heldout_auc([B1[i] for i in idx[:h]], [B1[i] for i in idx[h:]], rng, args.min_anchor)
            if nf is not None: rec["gate_noise_h1"] = round(float(nf), 4)
        if len(B2) >= 2 * args.min_anchor:
            idx = rng.permutation(len(B2)); h = len(B2) // 2
            nf = heldout_auc([B2[i] for i in idx[:h]], [B2[i] for i in idx[h:]], rng, args.min_anchor)
            if nf is not None: rec["gate_noise_h2"] = round(float(nf), 4)
        l, n = lean(Q1, A1, B1, args.min_anchor)
        if l is not None: rec["lean_h1"] = round(l, 4); rec["lean_h1_n"] = n
        l, n = lean(Q2, A2, B2, args.min_anchor)
        if l is not None: rec["lean_h2"] = round(l, 4); rec["lean_h2_n"] = n
        if any(k in rec for k in ("gate_auc_h1", "gate_auc_h2", "lean_h1", "lean_h2")):
            rows.append(rec)
    bam.close()
    def med(key):
        xs = [r[key] for r in rows if key in r]
        return round(float(np.median(xs)), 4) if xs else None
    # lean 加權平均（read 數加權）
    def wlean(lk, nk):
        num = sum(r[lk] * r[nk] for r in rows if lk in r); den = sum(r[nk] for r in rows if lk in r)
        return round(num / den, 4) if den else None, den
    wl1, n1 = wlean("lean_h1", "lean_h1_n"); wl2, n2 = wlean("lean_h2", "lean_h2_n")
    out = {"chrom": args.chrom, "n_somatic_sites": len(sites), "n_windows": len(rows),
           "GATE_subclone_vs_ancestral_heldout_auc": {"h1_median": med("gate_auc_h1"), "h2_median": med("gate_auc_h2"),
               "h1_NOISE_median": med("gate_noise_h1"), "h2_NOISE_median": med("gate_noise_h2"),
               "h1_farCpG_median": med("gate_auc_h1_farCpG"), "h2_farCpG_median": med("gate_auc_h2_farCpG"),
               "h1_n": sum(1 for r in rows if "gate_auc_h1" in r), "h2_n": sum(1 for r in rows if "gate_auc_h2" in r)},
           "USERQ_ambiguous_lean_to_subclone": {"h1_weighted_frac_closer_to_subclone": wl1, "h1_n_reads": n1,
               "h2_weighted_frac_closer_to_subclone": wl2, "h2_n_reads": n2, "null_expectation": 0.5},
           "rows": rows}
    outp = f"{args.out_dir}/t3_local_allele_{args.chrom}.json"
    json.dump(out, open(outp, "w"), indent=1, ensure_ascii=False)
    g = out["GATE_subclone_vs_ancestral_heldout_auc"]; u = out["USERQ_ambiguous_lean_to_subclone"]
    print(f"[t3-local] {args.chrom}: win={len(rows)} | GATE H1={g['h1_median']}(noise {g['h1_NOISE_median']}, farCpG {g['h1_farCpG_median']}) H2={g['h2_median']}(noise {g['h2_NOISE_median']}, farCpG {g['h2_farCpG_median']}) "
          f"| 偏向亞群 H1={u['h1_weighted_frac_closer_to_subclone']} H2={u['h2_weighted_frac_closer_to_subclone']} vs 0.5 -> {outp}")

if __name__ == "__main__":
    main()
