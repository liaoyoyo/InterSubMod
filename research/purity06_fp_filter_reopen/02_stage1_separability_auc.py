#!/usr/bin/env python3
"""
Stage 1 — Separability: can ISM features separate TP vs FP at 60% purity?
Reopen of the 4x-NEGATIVE "ISM/methyl -> delete FP -> raise F1" direction,
new precondition C3 = 60% purity (HCC1395 DORADO t30_n20 subsample, MM/ML tagged).

Read-only on existing ISM outputs (no BAM, no re-run). Zero-cost characterization.

Pre-registered (strict reopen, scientific-rigor 8.3.1):
  H1: axis-B / ISM read-level feature AUC(TP vs FP) > 0.58 AND survives within-AF-bin confound.
  Ceiling context (research/.../ceiling_math.json): even a PERFECT FP oracle only yields
  dF1 = +0.0031 at this stage; to clear Cohen-small +0.01 under TP-loss<=2% needs single-feature
  AUC ~0.82-0.91 (here even IMPOSSIBLE at 60% longphase). So Stage1 is characterization only:
  it cannot rescue Stage2, but quantifies HOW separable TP/FP actually are.

Output: research/purity06_fp_filter_reopen/stage1_auc.json + stdout table.
"""
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

BASE = ("/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/"
        "s-pure-pileup/HCC1395_DORADO/purity_t30_n20_20260213_dorado_purity_full")
TP_CSV = f"{BASE}/intersubmod_tp/significance_summary.csv"
FP_CSV = f"{BASE}/intersubmod_fp/significance_summary.csv"
OUT = "research/purity06_fp_filter_reopen/stage1_auc.json"

SEED = 42
N_BOOT = 1000
AUC_SIG = 0.58          # screening bar (project convention)
AUC_NEEDED_F1 = 0.82    # from ceiling_math: min single-feature AUC to clear +0.01 (best case)

ID_COLS = {"RegionID", "Chr", "Pos", "Ref", "Alt", "label"}
# non-ordinal / categorical text columns to skip for AUC
CAT_COLS = {"DominantLabel", "VerificationClass", "LOH_Subtype", "Coverage_Category",
            "Quality_Tier", "UnassignedDir", "Significant", "SuggestFilter",
            "PassedGating", "HPMergedSig", "HPFineSig", "AlleleSig", "Potential_LOH"}

def boot_auc_ci(y, x, n=N_BOOT, seed=SEED):
    m = ~np.isnan(x)
    y2, x2 = y[m], x[m]
    nfp = int(y2.sum())
    if len(np.unique(y2)) < 2 or nfp < 5:
        return None
    auc = roc_auc_score(y2, x2)
    rng = np.random.default_rng(seed)
    idx = np.arange(len(y2))
    aucs = []
    for _ in range(n):
        s = rng.choice(idx, len(idx), replace=True)
        if len(np.unique(y2[s])) < 2:
            continue
        aucs.append(roc_auc_score(y2[s], x2[s]))
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return float(auc), float(lo), float(hi), int(m.sum()), nfp

def main():
    tp = pd.read_csv(TP_CSV)
    fp = pd.read_csv(FP_CSV)
    tp["label"] = 0   # TP = negative class (protect)
    fp["label"] = 1   # FP = positive class (want to detect/delete)
    df = pd.concat([tp, fp], ignore_index=True)
    y = df["label"].values.astype(float)

    feats = [c for c in df.columns if c not in ID_COLS and c not in CAT_COLS]
    rows = []
    for c in feats:
        x = pd.to_numeric(df[c], errors="coerce").values.astype(float)
        if np.all(np.isnan(x)):
            continue
        r = boot_auc_ci(y, x)
        if r is None:
            continue
        auc, lo, hi, n_used, n_fp = r
        oriented = max(auc, 1 - auc)
        rows.append({
            "feature": c,
            "AUC_raw": round(auc, 4),
            "AUC_oriented": round(oriented, 4),
            "CI95": [round(lo, 3), round(hi, 3)],
            "direction": "FP>TP" if auc >= 0.5 else "FP<TP",
            "n_used": n_used, "n_FP_used": n_fp,
        })
    rows.sort(key=lambda d: -d["AUC_oriented"])

    n_sig = sum(1 for r in rows if r["AUC_oriented"] >= AUC_SIG)
    best = rows[0] if rows else None
    summary = {
        "purity": "60% (DORADO t30_n20)", "stage": "LongPhase (ISM input set)",
        "n_TP": int((y == 0).sum()), "n_FP": int((y == 1).sum()),
        "n_features_scored": len(rows),
        "n_features_AUC_ge_0.58": n_sig,
        "best_feature": best,
        "AUC_needed_for_F1_plus0.01": AUC_NEEDED_F1,
        "verdict_stage1": ("separable" if best and best["AUC_oriented"] >= AUC_SIG else "not-separable"),
        "note": ("Stage1 measures separability only. Per ceiling_math.json the F1 outcome is "
                 "fixed (perfect oracle +0.0031, realistic negative) regardless of this AUC, "
                 "because FP density=0.89% (base-rate) and FN=10901 dominate (recall-bound)."),
    }
    out = {"summary": summary, "features": rows}
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print(f"=== Stage 1 separability @ 60% purity (TP={summary['n_TP']} FP={summary['n_FP']}) ===")
    print(f"{'feature':28} {'AUC_or':>7} {'CI95':>16} {'dir':>7}")
    for r in rows[:18]:
        print(f"{r['feature']:28} {r['AUC_oriented']:>7.3f} "
              f"[{r['CI95'][0]:.2f},{r['CI95'][1]:.2f}]".rjust(17) + f" {r['direction']:>7}")
    print(f"\nfeatures scored={len(rows)}  AUC>=0.58: {n_sig}  "
          f"(need AUC>={AUC_NEEDED_F1} just to MAYBE clear +0.01 -> {'NONE reach it' if (not best or best['AUC_oriented']<AUC_NEEDED_F1) else 'some reach'})")
    print(f"[written] {OUT}")

if __name__ == "__main__":
    main()
