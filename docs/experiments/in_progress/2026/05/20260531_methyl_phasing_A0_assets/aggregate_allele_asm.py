#!/usr/bin/env python3
"""
彙總 allele-ASM AUC：tumor vs normal 配對 → 裁決「甲基 allele 差異是 copy / 方法樂觀 / 真 haplotype ASM」。

核心比較：
  1. 整體 tumor vs normal AUC median（normal=copy-clean 二倍體；若 ≈ → 高 AUC 非 copy、非 tumor-特異）
  2. 配對 delta = tumor_auc − normal_auc（同位點）；median ≈0 → 無 tumor-特異/copy-特異增益
  3. by tumor-CN-status：neutral(CN=2 兩樣本皆二倍體=最乾淨對照) / gain / loh
  4. depth-matched AUC 分布（排除 P-06 read-count confound）
  5. normal 作生物學 null：normal random-het AUC 分布（copy-clean、多數非 ASM）
"""
import json, glob, os
import numpy as np

AD = os.path.dirname(os.path.abspath(__file__))

def load(sample):
    out = {}
    for fp in glob.glob(f"{AD}/allele_asm_{sample}_chr*.json"):
        d = json.load(open(fp))
        for r in d["regions"]:
            out[(d["chrom"], r["pos"])] = r
    return out

T = load("tumor"); N = load("normal")
common = sorted(set(T) & set(N))

def med(xs):
    xs = [x for x in xs if x is not None]
    return round(float(np.median(xs)), 4) if xs else None

def summ(regs):
    a = [r["auc"] for r in regs]
    adm = [r["auc_depthmatched"] for r in regs if r["auc_depthmatched"] is not None]
    sp = [r["shuffle_p95"] for r in regs if r["shuffle_p95"] is not None]
    sig = [1 if r["sig_over_shuffle"] else 0 for r in regs]
    return {
        "n": len(regs), "auc_median": med(a), "auc_mean": round(float(np.mean(a)), 4) if a else None,
        "auc_depthmatched_median": med(adm),
        "shuffle_p95_median": med(sp),
        "frac_sig_over_shuffle": round(float(np.mean(sig)), 4) if sig else None,
    }

res = {
    "n_tumor_regions": len(T), "n_normal_regions": len(N), "n_paired_loci": len(common),
    "overall": {"tumor": summ(list(T.values())), "normal": summ(list(N.values()))},
    "by_tumor_status": {},
    "paired": {}, "by_chrom": {},
}

# by tumor CN status (normal at SAME loci as comparison)
for st in ["neutral", "gain", "loss", "loh"]:
    tlist = [r for r in T.values() if r["status"] == st]
    if not tlist:
        continue
    # normal 在 tumor 為此 status 的同位點
    nlist = [N[k] for k in common if T[k]["status"] == st]
    res["by_tumor_status"][st] = {
        "tumor": summ(tlist),
        "normal_same_loci": summ(nlist) if nlist else None,
    }

# 配對 delta（同位點 tumor − normal）
deltas = [T[k]["auc"] - N[k]["auc"] for k in common]
deltas_neu = [T[k]["auc"] - N[k]["auc"] for k in common if T[k]["status"] == "neutral"]
deltas_loh = [T[k]["auc"] - N[k]["auc"] for k in common if T[k]["status"] == "loh"]
res["paired"] = {
    "n": len(common),
    "delta_tumor_minus_normal_median": med(deltas),
    "delta_mean": round(float(np.mean(deltas)), 4) if deltas else None,
    "frac_tumor_gt_normal": round(float(np.mean([1 if d > 0 else 0 for d in deltas])), 4) if deltas else None,
    "neutral_delta_median": med(deltas_neu), "neutral_n": len(deltas_neu),
    "loh_delta_median": med(deltas_loh), "loh_n": len(deltas_loh),
}

for c in sorted({k[0] for k in T}):
    tl = [r for k, r in T.items() if k[0] == c]
    nl = [r for k, r in N.items() if k[0] == c]
    res["by_chrom"][c] = {"tumor_auc_median": summ(tl)["auc_median"], "tumor_n": len(tl),
                          "normal_auc_median": summ(nl)["auc_median"] if nl else None, "normal_n": len(nl)}

outp = AD + "/allele_asm_aggregate.json"
with open(outp, "w") as f:
    json.dump(res, f, indent=1, ensure_ascii=False)

print(f"paired loci: {len(common)}")
print(f"OVERALL  tumor auc_med={res['overall']['tumor']['auc_median']} (dm={res['overall']['tumor']['auc_depthmatched_median']}, sig={res['overall']['tumor']['frac_sig_over_shuffle']})")
print(f"         normal auc_med={res['overall']['normal']['auc_median']} (dm={res['overall']['normal']['auc_depthmatched_median']}, sig={res['overall']['normal']['frac_sig_over_shuffle']})")
print(f"PAIRED   delta(tumor-normal) median={res['paired']['delta_tumor_minus_normal_median']}  frac tumor>normal={res['paired']['frac_tumor_gt_normal']}")
print(f"         neutral(CN=2) delta_med={res['paired']['neutral_delta_median']} (n={res['paired']['neutral_n']}); loh delta_med={res['paired']['loh_delta_median']} (n={res['paired']['loh_n']})")
for st, v in res["by_tumor_status"].items():
    nn = v["normal_same_loci"]["auc_median"] if v["normal_same_loci"] else None
    print(f"  status {st:7s}: tumor auc_med={v['tumor']['auc_median']} (n={v['tumor']['n']})  normal@same={nn}")
print(f"-> {outp}")
