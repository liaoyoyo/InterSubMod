#!/usr/bin/env python3
"""
imprinting 正控：已知 imprinted DMR 的 per-read 甲基雙峰性（不需 HP anchor）。

邏輯：imprinted DMR 一條親代 allele 甲基化、一條未甲基化（parent-of-origin）→ per-read 平均甲基
應呈 ~50:50 雙峰。這是「甲基測量本身可在已知 ASM 位點乾淨分兩群」的正控，
獨立於 copy number、獨立於 germline het anchor（解 V2 imprinting INCONCLUSIVE 的 anchor 死結）。

讀既有 *_matrix.tsv（read×CpG），算 per-read mean methylation → GMM 1 vs 2 comp BIC + 雙峰分裂比。
⚠ caveat：此處 matrix 來自 TUMOR BAM；tumor 可能有 loss-of-imprinting(LOI)，故 NEGATIVE 不必然=方法失敗。
   最乾淨正控應在 normal（見 --normal-extract 由另一步補）。
"""
import sys, json, glob, os
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

def analyze(matrix_path, min_cpg=5):
    m = pd.read_csv(matrix_path, sep="\t", index_col=0)
    # per-read 平均甲基（只算非 NaN），要求每 read ≥ min_cpg 個 CpG
    n_per_read = m.notna().sum(axis=1)
    keep = n_per_read >= min_cpg
    rm = m[keep].mean(axis=1, skipna=True).dropna().values
    n = len(rm)
    if n < 12:
        return {"n_reads_used": int(n), "verdict": "TOO_FEW_READS"}
    x = rm.reshape(-1, 1)
    g1 = GaussianMixture(1, random_state=0).fit(x)
    g2 = GaussianMixture(2, random_state=0).fit(x)
    bic_diff = float(g1.bic(x) - g2.bic(x))  # >10 = 強雙峰
    means = sorted(g2.means_.ravel().tolist())
    labels = g2.predict(x)
    # 兩群比例 + 兩群中心甲基差
    frac_hi = float(np.mean(labels == np.argmax(g2.means_.ravel())))
    delta_centers = float(abs(means[1] - means[0]))
    # 簡單分裂：低甲基 read 比例（mean<0.3）vs 高甲基（>0.7）
    frac_low = float(np.mean(rm < 0.3)); frac_high = float(np.mean(rm > 0.7)); frac_mid = float(np.mean((rm >= 0.3) & (rm <= 0.7)))
    verdict = "BIMODAL" if (bic_diff > 10 and delta_centers > 0.25) else ("WEAK_BIMODAL" if bic_diff > 2 else "UNIMODAL")
    return {
        "n_reads_used": int(n),
        "bic_diff_2vs1": round(bic_diff, 2),
        "gmm_center_low": round(means[0], 3), "gmm_center_high": round(means[1], 3),
        "delta_centers": round(delta_centers, 3),
        "frac_in_high_cluster": round(frac_hi, 3),
        "frac_low_meth(<0.3)": round(frac_low, 3), "frac_mid": round(frac_mid, 3), "frac_high_meth(>0.7)": round(frac_high, 3),
        "rm_mean": round(float(np.mean(rm)), 3), "rm_std": round(float(np.std(rm)), 3),
        "verdict": verdict,
    }

def main():
    AD = os.path.dirname(os.path.abspath(__file__))
    targets = {}
    for fp in sorted(glob.glob(AD + "/*_matrix.tsv")):
        base = os.path.basename(fp).replace("_matrix.tsv", "")
        if any(k in base for k in ("imprint", "h19", "gnas", "snrpn")):
            targets[base] = fp
    out = {"note": "imprinted DMR per-read methylation bimodality (positive control, no anchor needed)",
           "source": "TUMOR BAM matrices (caveat: tumor LOI possible)", "loci": {}}
    for name, fp in targets.items():
        out["loci"][name] = analyze(fp)
    outp = AD + "/imprint_bimodality.json"
    with open(outp, "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    for name, r in out["loci"].items():
        print(f"  {name}: {r.get('verdict')} n={r.get('n_reads_used')} bic={r.get('bic_diff_2vs1')} "
              f"centers={r.get('gmm_center_low')}/{r.get('gmm_center_high')} Δ={r.get('delta_centers')}")
    print(f"-> {outp}")

if __name__ == "__main__":
    main()
