#!/usr/bin/env python3
"""
T3 嚴格反駁/確認：subclone(1-1) 甲基是否與 ancestral(1) 可區分？用戶前提的必要條件。

用戶 T3 前提：亞群 read（有/沒踩 somatic 位點）甲基一致 → 可用有踩到的校正沒踩到的。
此前提的【必要條件】= 亞群(1-1) 甲基 ≠ 母本(1) 甲基（否則無從區分）。

決定性測試（不只 AUC，小群+對稱化會膨脹 null）：
  - raw Δβ(1-1 vs 1) = 兩群質心在共享 CpG 的平均甲基差（實際 effect size）
  - 對照 raw Δβ(1 vs 2) = 已知真實 haplotype 差（V10 正控）── 在【同窗】比較
  - directional AUC（不對稱化）+ 各群 n + group-1 甲基雙峰 BIC（隱藏亞群該雙峰）
  判讀：若 Δβ(1-1,1) ≈ 0 ≪ Δβ(1,2) → 必要條件不成立（亞群無獨特甲基）→ T3 NEGATIVE 確認。
唯讀。每 chr 落 rigor_t3_{chrom}.json。
"""
import sys, json, argparse
import numpy as np
import pysam
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
SOMATIC_VCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/somatic_pass.vcf.gz"

def get_hp(read):
    try: return str(read.get_tag("HP"))
    except KeyError: return "unphase"

def collect(bam, chrom, center0, win, tags, min_cpg=5):
    start = max(0, center0 - win); end = center0 + win
    groups = {t: [] for t in tags}; seen = set()
    for read in bam.fetch(chrom, start, end):
        if read.is_secondary or read.is_supplementary or read.is_unmapped: continue
        if read.query_name in seen: continue
        hp = get_hp(read)
        if hp not in tags: continue
        mods = read.modified_bases or {}; mc = None
        for k, calls in mods.items():
            if k[0] in ("C", b"C") and k[2] in ("m", 27551): mc = calls; break
        if not mc: continue
        q2r = {q: r for q, r in read.get_aligned_pairs(matches_only=True)}
        meth = {}
        for qpos, qual in mc:
            rpos = q2r.get(qpos)
            if rpos is None or rpos < start or rpos > end: continue
            meth[rpos] = qual / 255.0
        if len(meth) < min_cpg: continue
        seen.add(read.query_name); groups[hp].append(meth)
    return groups

def _dbeta(g0, g1):
    """兩群在共享 CpG 的平均 |Δβ|（per-CpG 質心差）+ pattern-based LOO 質心 AUC。"""
    pos = set()
    for m in g0 + g1: pos.update(m.keys())
    pos = [p for p in pos if sum(1 for m in g0 + g1 if p in m) >= 0.3 * (len(g0) + len(g1))]
    if len(pos) < 5: return None
    def cen(g):
        v = []
        for p in pos:
            xs = [m[p] for m in g if p in m]
            v.append(np.mean(xs) if xs else np.nan)
        return np.array(v)
    c0, c1 = cen(g0), cen(g1)
    mask = ~np.isnan(c0) & ~np.isnan(c1)
    if mask.sum() < 5: return None
    dbeta = float(np.mean(np.abs(c0[mask] - c1[mask])))  # 平均 |Δβ| per CpG（pattern）
    # pattern-based separability：每 read 到兩質心 RMSE 距離，scores=d0-d1
    rows = [(0, m) for m in g0] + [(1, m) for m in g1]
    sc, lb = [], []
    for a, m in rows:
        v0 = [(m[p] - c0[i]) ** 2 for i, p in enumerate(pos) if p in m and not np.isnan(c0[i])]
        v1 = [(m[p] - c1[i]) ** 2 for i, p in enumerate(pos) if p in m and not np.isnan(c1[i])]
        if len(v0) < 3 or len(v1) < 3: continue
        sc.append(np.sqrt(np.mean(v0)) - np.sqrt(np.mean(v1))); lb.append(a)
    auc = None
    if len(set(lb)) == 2 and len(lb) >= 6:
        try: a = roc_auc_score(lb, sc); auc = max(a, 1 - a)
        except Exception: auc = None
    return {"dbeta": round(dbeta, 4), "pat_auc": round(float(auc), 4) if auc is not None else None, "n_cpg": int(mask.sum())}

def centroid_diff(g0, g1):
    return _dbeta(g0, g1)

def noise_floor(g, rng):
    """同群隨機對半 → Δβ 噪音地板（純抽樣噪音，非生物訊號）。"""
    if len(g) < 16: return None
    idx = rng.permutation(len(g)); h = len(g) // 2
    return _dbeta([g[i] for i in idx[:h]], [g[i] for i in idx[h:]])

def bimod(g):
    rm = np.array([np.mean(list(m.values())) for m in g]).reshape(-1, 1)
    if len(rm) < 12: return None
    try:
        b = float(GaussianMixture(1, random_state=0).fit(rm).bic(rm) - GaussianMixture(2, random_state=0).fit(rm).bic(rm))
        return round(b, 2)
    except Exception:
        return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--max-sites", type=int, default=250)
    ap.add_argument("--win", type=int, default=2000)
    ap.add_argument("--min-anchor", type=int, default=8)
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--seed", type=int, default=20260609)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    vcf = pysam.VariantFile(SOMATIC_VCF)
    sites = [rec.pos - 1 for rec in vcf.fetch(args.chrom) if rec.alts and len(rec.ref) == 1 and len(rec.alts[0]) == 1]
    vcf.close()
    step = max(1, len(sites) // (args.max_sites * 3)) if sites else 1
    bam = pysam.AlignmentFile(BAM, "rb")
    rows = []
    for c in sites[::step]:
        if len(rows) >= args.max_sites: break
        g = collect(bam, args.chrom, c, args.win, {"1", "2", "1-1", "2-1"})
        rec = {"site": c + 1, "n_1": len(g["1"]), "n_2": len(g["2"]), "n_11": len(g["1-1"]), "n_21": len(g["2-1"])}
        # T3-H1: 亞群 vs 母本（同 haplotype）
        if len(g["1-1"]) >= args.min_anchor and len(g["1"]) >= args.min_anchor:
            r = centroid_diff(g["1-1"], g["1"])
            if r: rec["t3_h1"] = r
        # T3-H2
        if len(g["2-1"]) >= args.min_anchor and len(g["2"]) >= args.min_anchor:
            r = centroid_diff(g["2-1"], g["2"])
            if r: rec["t3_h2"] = r
        # 正控：H1 vs H2（不同 haplotype，已知真實）— 同窗
        if len(g["1"]) >= args.min_anchor and len(g["2"]) >= args.min_anchor:
            r = centroid_diff(g["1"], g["2"])
            if r: rec["ctrl_h1_vs_h2"] = r
        # 噪音地板：同群隨機對半（純抽樣噪音）
        nf = noise_floor(g["1"], rng)
        if nf: rec["noise_floor_1"] = nf
        nf2 = noise_floor(g["2"], rng)
        if nf2: rec["noise_floor_2"] = nf2
        # group-1 雙峰（隱藏亞群？）
        if len(g["1"]) >= 12:
            rec["bimod_group1"] = bimod(g["1"])
        if "t3_h1" in rec or "t3_h2" in rec or "ctrl_h1_vs_h2" in rec:
            rows.append(rec)
    bam.close()

    def med(key, sub):
        xs = [r[key][sub] for r in rows if key in r and r[key].get(sub) is not None]
        return round(float(np.median(xs)), 4) if xs else None
    out = {"chrom": args.chrom, "n_somatic_sites": len(sites), "n_windows": len(rows),
           "T3_H1_subclone_vs_ancestral": {"dbeta_median": med("t3_h1", "dbeta"), "pat_auc_median": med("t3_h1", "pat_auc"),
                                            "n_windows": sum(1 for r in rows if "t3_h1" in r)},
           "T3_H2_subclone_vs_ancestral": {"dbeta_median": med("t3_h2", "dbeta"), "pat_auc_median": med("t3_h2", "pat_auc"),
                                            "n_windows": sum(1 for r in rows if "t3_h2" in r)},
           "CTRL_H1_vs_H2_positive": {"dbeta_median": med("ctrl_h1_vs_h2", "dbeta"), "pat_auc_median": med("ctrl_h1_vs_h2", "pat_auc"),
                                       "n_windows": sum(1 for r in rows if "ctrl_h1_vs_h2" in r)},
           "NOISE_FLOOR_random_split": {"dbeta_median_g1": med("noise_floor_1", "dbeta"), "dbeta_median_g2": med("noise_floor_2", "dbeta"),
                                        "pat_auc_median_g1": med("noise_floor_1", "pat_auc")},
           "bimod_group1_frac_strong(BIC>10)": round(float(np.mean([1 if (r.get("bimod_group1") or 0) > 10 else 0 for r in rows if "bimod_group1" in r])), 4) if any("bimod_group1" in r for r in rows) else None,
           "n11_median": int(np.median([r["n_11"] for r in rows])) if rows else 0,
           "rows": rows}
    outp = f"{args.out_dir}/rigor_t3_{args.chrom}.json"
    json.dump(out, open(outp, "w"), indent=1, ensure_ascii=False)
    t3 = out["T3_H1_subclone_vs_ancestral"]; ct = out["CTRL_H1_vs_H2_positive"]; nf = out["NOISE_FLOOR_random_split"]
    print(f"[rigor-t3] {args.chrom}: win={len(rows)} | T3-H1 Δβ={t3['dbeta_median']} patAUC={t3['pat_auc_median']} (n={t3['n_windows']}) "
          f"| 正控 H1vsH2 Δβ={ct['dbeta_median']} patAUC={ct['pat_auc_median']} | 噪音地板 Δβ={nf['dbeta_median_g1']} patAUC={nf['pat_auc_median_g1']} | n11={out['n11_median']} -> {outp}")

if __name__ == "__main__":
    main()
