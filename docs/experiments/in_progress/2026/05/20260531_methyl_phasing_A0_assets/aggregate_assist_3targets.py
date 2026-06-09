#!/usr/bin/env python3
"""彙總 T1(assist_tag) + T2/T3(methyl_assist) 三 target → 一份判別表 + 3-panel 圖。數字全讀 JSON。"""
import json, glob, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK TC", "Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
AD = os.path.dirname(os.path.abspath(__file__))

# ---- T1: assist_tag (unphase->H1/H2 at separable loci) ----
sep_acc, nsep_acc, sep_yield, sep_cacc, sep_n, nsep_n = [], [], [], [], 0, 0
for f in glob.glob(AD + "/assist_tag_chr*.json"):
    d = json.load(open(f)); s = d["verified_separable_train_sep>=0.85"]; ns = d["not_separable_train_sep<0.85"]
    if s["test_acc_median"] is not None: sep_acc.append(s["test_acc_median"]); sep_n += s["n"]
    if s["conf_yield_median"] is not None: sep_yield.append(s["conf_yield_median"])
    if s["conf_acc_median"] is not None: sep_cacc.append(s["conf_acc_median"])
    if ns["test_acc_median"] is not None: nsep_acc.append(ns["test_acc_median"]); nsep_n += ns["n"]
T1 = {"separable_n": sep_n, "separable_test_acc": round(float(np.median(sep_acc)), 4) if sep_acc else None,
      "separable_conf_yield": round(float(np.median(sep_yield)), 4) if sep_yield else None,
      "separable_conf_acc": round(float(np.median(sep_cacc)), 4) if sep_cacc else None,
      "notsep_n": nsep_n, "notsep_test_acc": round(float(np.median(nsep_acc)), 4) if nsep_acc else None}

# ---- T2/T3: methyl_assist ----
def pool(key):
    aucs, sigs, sp, n = [], [], [], 0
    for f in glob.glob(AD + "/methyl_assist_chr*.json"):
        d = json.load(open(f)); blk = d[key]
        if blk.get("n_sites", 0) == 0: continue
        # 用 per-site list 算 pooled median
        listkey = {"T3_H1_subclone_vs_germline": "sites_t3_h1", "T3_H2_subclone_vs_germline": "sites_t3_h2", "T2_H1-1_vs_H2-1_branch": "sites_t2"}[key]
        for s in d.get(listkey, []):
            aucs.append(s["auc"]); sigs.append(1 if s["sig"] else 0)
            if s["shuffle_p95"] is not None: sp.append(s["shuffle_p95"])
        n += blk["n_sites"]
    return {"n_sites": n, "auc_median": round(float(np.median(aucs)), 4) if aucs else None,
            "frac_sig": round(float(np.mean(sigs)), 4) if sigs else None,
            "shuffle_p95_median": round(float(np.median(sp)), 4) if sp else None}
T2 = pool("T2_H1-1_vs_H2-1_branch")
T3H1 = pool("T3_H1_subclone_vs_germline"); T3H2 = pool("T3_H2_subclone_vs_germline")

res = {"T1_unphase_to_H1H2": T1, "T2_H3_to_branch": T2, "T3_H1_subclone": T3H1, "T3_H2_subclone": T3H2}
json.dump(res, open(AD + "/assist_3targets_aggregate.json", "w"), indent=1, ensure_ascii=False)

# ---- 圖 ----
fig, ax = plt.subplots(figsize=(9.5, 5))
labels = ["T1\nunphase→H1/H2\n(可分位點)", "T2\nH3→H1-1/H2-1\n(germline 分支)", "T3-H1\n1-1 vs 1\n(亞群 vs 母本)", "T3-H2\n2-1 vs 2\n(亞群 vs 母本)"]
vals = [T1["separable_test_acc"], T2["auc_median"], T3H1["auc_median"], T3H2["auc_median"]]
nulls = [0.5, T2["shuffle_p95_median"], T3H1["shuffle_p95_median"], T3H2["shuffle_p95_median"]]
colors = ["#15803D", "#15803D", "#C2410C", "#C2410C"]
x = np.arange(len(labels))
ax.bar(x - 0.2, vals, 0.4, color=colors, alpha=0.85, label="real (AUC/acc)")
ax.bar(x + 0.2, nulls, 0.4, color="#9CA3AF", alpha=0.7, label="null (shuffle p95 / chance)")
for i, (v, nl) in enumerate(zip(vals, nulls)):
    if v is not None: ax.text(i - 0.2, v + 0.01, f"{v:.2f}", ha="center", fontsize=9)
ax.axhline(0.5, ls=":", c="gray", lw=1)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("判別力 (AUC / held-out accuracy)"); ax.set_ylim(0, 1.05)
ax.set_title("甲基輔助 longphase-S tag 矯正：T1/T2 可行（分不同 haplotype），T3 弱（同 haplotype 內亞群）\n機制：甲基訊號是 germline-haplotype 層級，非 subclone 層級")
ax.legend(fontsize=9, loc="lower right")
fig.tight_layout(); fig.savefig(AD + "/fig_assist_3targets.png", dpi=130, bbox_inches="tight"); plt.close()

print("=== 3-target 彙總 ===")
print(f"T1 unphase→H1/H2 (可分位點 n={T1['separable_n']}): held-out acc={T1['separable_test_acc']} conf_yield={T1['separable_conf_yield']} conf_acc={T1['separable_conf_acc']} | 不可分 acc={T1['notsep_test_acc']}")
print(f"T2 H3→branch (n={T2['n_sites']}): auc={T2['auc_median']} (null {T2['shuffle_p95_median']}) frac_sig={T2['frac_sig']}")
print(f"T3-H1 亞群vs母本 (n={T3H1['n_sites']}): auc={T3H1['auc_median']} (null {T3H1['shuffle_p95_median']}) frac_sig={T3H1['frac_sig']}")
print(f"T3-H2 亞群vs母本 (n={T3H2['n_sites']}): auc={T3H2['auc_median']} (null {T3H2['shuffle_p95_median']}) frac_sig={T3H2['frac_sig']}")
print("saved fig_assist_3targets.png + assist_3targets_aggregate.json")
