#!/usr/bin/env python3
"""
HANDOFF 3b/3c/3d — tumor HP-AUC < normal 配對 delta + confound 分層 (CNV/LOH/coverage/HP-balance) + 顯著性關係。
回答：tumor<normal 在「乾淨」diploid 層是否仍成立 → H_somatic vs H_cnv。
read-only（讀 significance_summary.csv），數字落 JSON 供 HTML 報告引用。
依據：docs/.../method_audit/HANDOFF_tumor_vs_normal_hpauc.md §3b/3c/3d/4
"""
import json, sys
import pandas as pd
import numpy as np

CSV = sys.argv[1] if len(sys.argv) > 1 else "/tmp/wg_hpauc/significance_summary.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "method_audit/confound_stratify.json"

cols = ["Chr","Pos","NumReads","NReadsValid","HP_AUC_Normal","HP_AUC_Tumor","CramersV",
        "Potential_LOH","Coverage_Category","NTumorReads","NNormalReads","Tumor_HP1","Tumor_HP2","Significant"]
df = pd.read_csv(CSV, usecols=cols, low_memory=False)

# normalize types
for c in ["HP_AUC_Normal","HP_AUC_Tumor","CramersV"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")
for c in ["Tumor_HP1","Tumor_HP2","NTumorReads","NNormalReads","NumReads"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
df["LOH"] = df["Potential_LOH"].astype(str).str.lower() == "true"
df["Sig"] = df["Significant"].astype(str).str.lower() == "true"

def stats(d):
    if len(d) == 0:
        return {"n": 0}
    delta = (d["HP_AUC_Normal"] - d["HP_AUC_Tumor"])
    return {
        "n": int(len(d)),
        "median_delta": round(float(delta.median()), 4),
        "mean_delta": round(float(delta.mean()), 4),
        "frac_delta_gt0": round(float((delta > 0).mean()), 3),
        "frac_delta_gt0.2": round(float((delta > 0.2).mean()), 3),
        "median_tumor_auc": round(float(d["HP_AUC_Tumor"].median()), 4),
        "median_normal_auc": round(float(d["HP_AUC_Normal"].median()), 4),
    }

# ---- both-HP valid subset (HANDOFF 3b) ----
both = df[(df["Tumor_HP1"] > 0) & (df["Tumor_HP2"] > 0) &
         (df["HP_AUC_Normal"] >= 0) & (df["HP_AUC_Tumor"] >= 0)].copy()

out = {
    "input": CSV,
    "n_regions_total": int(len(df)),
    "3b_paired_delta_both_hp": stats(both),
}

# ---- 3c-CNV: 分 Coverage_Category ----
out["3c_by_coverage_category"] = {
    cat: stats(both[both["Coverage_Category"] == cat])
    for cat in ["Normal","Elevated","CNV_Gain","CNV_Loss","High_Copy","Low"]
}

# ---- 3c-LOH ----
out["3c_by_loh"] = {
    "LOH_true": stats(both[both["LOH"]]),
    "LOH_false": stats(both[~both["LOH"]]),
}

# ---- 🔴 CLEAN layer: diploid Normal + non-LOH ----
clean = both[(both["Coverage_Category"] == "Normal") & (~both["LOH"])]
out["3c_clean_diploid_nonloh"] = stats(clean)

# coverage control inside clean layer (避免 tumor read 少→AUC 估計差)
out["3c_clean_coverage_control"] = {
    "NTumorReads_ge30": stats(clean[clean["NTumorReads"] >= 30]),
    "NTumorReads_ge50": stats(clean[clean["NTumorReads"] >= 50]),
    "NTumorReads_ge80": stats(clean[clean["NTumorReads"] >= 80]),
}

# HP balance control inside clean layer (極不平衡→AUC 不穩)
mn = np.minimum(clean["Tumor_HP1"], clean["Tumor_HP2"])
mx = np.maximum(clean["Tumor_HP1"], clean["Tumor_HP2"]).replace(0, np.nan)
bal = (mn / mx)
out["3c_clean_hp_balance_ge0.33"] = stats(clean[bal >= 0.33])
out["3c_clean_hp_balance_ge0.5"] = stats(clean[bal >= 0.5])

# ---- 3d: tumor AUC vs 顯著性 / CramersV (within both-HP) ----
out["3d_significance_relation"] = {
    "median_tumor_auc_when_sig": round(float(both[both["Sig"]]["HP_AUC_Tumor"].median()), 4),
    "median_tumor_auc_when_notsig": round(float(both[~both["Sig"]]["HP_AUC_Tumor"].median()), 4),
    "median_cramersv_tumor_auc_lt0.5": round(float(both[both["HP_AUC_Tumor"] < 0.5]["CramersV"].median()), 4),
    "median_cramersv_tumor_auc_ge0.7": round(float(both[both["HP_AUC_Tumor"] >= 0.7]["CramersV"].median()), 4),
}

# ---- 4: 案例選取 ----
# 乾淨 somatic 候選: normal≥0.8 & tumor≤0.5 & 乾淨 diploid + 充足覆蓋, 每染色體取最大 delta 1 個
cand = clean[(clean["HP_AUC_Normal"] >= 0.8) & (clean["HP_AUC_Tumor"] <= 0.5) &
             (clean["NTumorReads"] >= 50) & (bal >= 0.33)].copy()
cand["delta"] = cand["HP_AUC_Normal"] - cand["HP_AUC_Tumor"]
cand = cand.sort_values("delta", ascending=False)
clean_cases = []
seen_chr = set()
for _, r in cand.iterrows():
    if r["Chr"] in seen_chr:
        continue
    seen_chr.add(r["Chr"])
    clean_cases.append({"chr": r["Chr"], "pos": int(r["Pos"]),
                        "normal_auc": round(float(r["HP_AUC_Normal"]),3), "tumor_auc": round(float(r["HP_AUC_Tumor"]),3),
                        "delta": round(float(r["delta"]),3), "NTumorReads": int(r["NTumorReads"]),
                        "Tumor_HP1": int(r["Tumor_HP1"]), "Tumor_HP2": int(r["Tumor_HP2"]),
                        "cov": r["Coverage_Category"], "loh": bool(r["LOH"]), "sig": bool(r["Sig"])})
    if len(clean_cases) >= 6:
        break
out["clean_somatic_candidate_count"] = int(len(cand))
out["clean_somatic_cases"] = clean_cases

# 對照案例: CNV_Gain 或 LOH 的 tumor<<normal (confound 版)
contrast = both[((both["Coverage_Category"] == "CNV_Gain") | (both["LOH"])) &
                (both["HP_AUC_Normal"] >= 0.8) & (both["HP_AUC_Tumor"] <= 0.5) &
                (both["NTumorReads"] >= 50)].copy()
contrast["delta"] = contrast["HP_AUC_Normal"] - contrast["HP_AUC_Tumor"]
contrast = contrast.sort_values("delta", ascending=False)
out["contrast_confound_cases"] = [
    {"chr": r["Chr"], "pos": int(r["Pos"]), "normal_auc": round(float(r["HP_AUC_Normal"]),3),
     "tumor_auc": round(float(r["HP_AUC_Tumor"]),3), "delta": round(float(r["delta"]),3),
     "cov": r["Coverage_Category"], "loh": bool(r["LOH"]), "NTumorReads": int(r["NTumorReads"])}
    for _, r in contrast.head(3).iterrows()
]

json.dump(out, open(OUT, "w"), indent=2, ensure_ascii=False)
print(f"[confound_stratify] wrote {OUT}")
print(f"  both-HP n={out['3b_paired_delta_both_hp']['n']} median_delta={out['3b_paired_delta_both_hp']['median_delta']}")
print(f"  CLEAN diploid+nonLOH: n={out['3c_clean_diploid_nonloh']['n']} median_delta={out['3c_clean_diploid_nonloh']['median_delta']} frac>0={out['3c_clean_diploid_nonloh']['frac_delta_gt0']}")
print(f"  clean+cov>=50: median_delta={out['3c_clean_coverage_control']['NTumorReads_ge50']['median_delta']}")
print(f"  clean+balance>=0.5: median_delta={out['3c_clean_hp_balance_ge0.5']['median_delta']}")
print(f"  clean somatic candidates={out['clean_somatic_candidate_count']} cases={len(clean_cases)}")
