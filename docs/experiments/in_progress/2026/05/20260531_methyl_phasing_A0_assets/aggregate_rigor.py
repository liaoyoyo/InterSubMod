#!/usr/bin/env python3
"""彙總 rigor t1/t2/t3 跨染色體 → 一份反例/邊緣/嚴格驗證表 + 圖。數字全讀 JSON（池化 per-window rows）。"""
import json, glob, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK TC", "Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
AD = os.path.dirname(os.path.abspath(__file__))

def med(xs):
    xs = [x for x in xs if x is not None]
    return round(float(np.median(xs)), 4) if xs else None

# ---- T3: pool per-window dbeta（subclone / ctrl / noise） ----
t3_sub, t3_ctrl, t3_noise, t3_sub_auc, t3_ctrl_auc, t3_noise_auc = [], [], [], [], [], []
for f in glob.glob(AD + "/rigor_t3_chr*.json"):
    d = json.load(open(f))
    for r in d["rows"]:
        if "t3_h1" in r: t3_sub.append(r["t3_h1"]["dbeta"]); t3_sub_auc.append(r["t3_h1"].get("pat_auc"))
        if "t3_h2" in r: t3_sub.append(r["t3_h2"]["dbeta"]); t3_sub_auc.append(r["t3_h2"].get("pat_auc"))
        if "ctrl_h1_vs_h2" in r: t3_ctrl.append(r["ctrl_h1_vs_h2"]["dbeta"]); t3_ctrl_auc.append(r["ctrl_h1_vs_h2"].get("pat_auc"))
        if "noise_floor_1" in r: t3_noise.append(r["noise_floor_1"]["dbeta"]); t3_noise_auc.append(r["noise_floor_1"].get("pat_auc"))
        if "noise_floor_2" in r: t3_noise.append(r["noise_floor_2"]["dbeta"]); t3_noise_auc.append(r["noise_floor_2"].get("pat_auc"))
T3 = {"n_windows": len(t3_sub),
      "subclone_dbeta": med(t3_sub), "ctrl_haplotype_dbeta": med(t3_ctrl), "noise_floor_dbeta": med(t3_noise),
      "subclone_excess_over_noise": round((med(t3_sub) or 0) - (med(t3_noise) or 0), 4),
      "ctrl_excess_over_noise": round((med(t3_ctrl) or 0) - (med(t3_noise) or 0), 4),
      "subclone_pat_auc": med(t3_sub_auc), "ctrl_pat_auc": med(t3_ctrl_auc), "noise_pat_auc": med(t3_noise_auc)}

# ---- T1: pool failure vs ok ----
fail, ok, sep, nsep = [], [], [], []
by_status = {}
for f in glob.glob(AD + "/rigor_t1_chr*.json"):
    d = json.load(open(f))
    for r in d["rows"]:
        if r["acc"] < 0.7: fail.append(r)
        elif r["acc"] >= 0.9: ok.append(r)
        if r["train_sep"] is not None:
            (sep if r["train_sep"] >= 0.85 else nsep).append(r)
        by_status.setdefault(r["status"], []).append(r)
def grp(rs):
    if not rs: return {"n": 0}
    return {"n": len(rs), "acc_median": med([r["acc"] for r in rs]), "n_cpg_median": int(np.median([r["n_cpg"] for r in rs])),
            "n_reads_median": int(np.median([r["n_reads"] for r in rs])),
            "train_sep_median": med([r["train_sep"] for r in rs if r["train_sep"] is not None])}
T1 = {"FAIL_acc<0.7": grp(fail), "OK_acc>=0.9": grp(ok),
      "separable_sep>=0.85": grp(sep), "notsep_sep<0.85": grp(nsep),
      "by_status": {s: grp(rs) for s, rs in by_status.items()}}

# ---- T2: pool H3 class + assignability ----
cls = {"no_germline": 0, "inconsistent_germline": 0, "weak_or_consistent": 0}
assignmargin = {"no_germline": [], "inconsistent_germline": [], "weak_or_consistent": []}
n_h3 = 0
for f in glob.glob(AD + "/rigor_t2_chr*.json"):
    d = json.load(open(f))
    n_h3 += d["n_H3_reads"]
    for k in cls: cls[k] += d["H3_germline_class"].get(k, 0)
T2 = {"n_H3_reads": n_h3, "H3_germline_class": cls,
      "H3_class_pct": {k: round(100 * v / max(1, sum(cls.values())), 1) for k, v in cls.items()},
      "H3_assignability_by_class": {}}
# assignability：從各 chr JSON 的彙總欄取（已是 median/frac）
for f in glob.glob(AD + "/rigor_t2_chr*.json"):
    d = json.load(open(f))
    for k, v in d["H3_assignability_by_class"].items():
        T2["H3_assignability_by_class"].setdefault(k, {"margin": [], "highconf": [], "n": 0})
        if v.get("margin_median") is not None: T2["H3_assignability_by_class"][k]["margin"].append(v["margin_median"])
        if v.get("frac_high_conf(margin>0.1)") is not None: T2["H3_assignability_by_class"][k]["highconf"].append(v["frac_high_conf(margin>0.1)"])
        T2["H3_assignability_by_class"][k]["n"] += v.get("n_with_centroid", 0)
for k in list(T2["H3_assignability_by_class"]):
    a = T2["H3_assignability_by_class"][k]
    T2["H3_assignability_by_class"][k] = {"n_with_centroid": a["n"], "margin_median": med(a["margin"]), "frac_high_conf": med(a["highconf"])}

res = {"T1_failure_analysis": T1, "T2_H3_edge": T2, "T3_subclone_rigor": T3}
json.dump(res, open(AD + "/rigor_aggregate.json", "w"), indent=1, ensure_ascii=False)

# ---- 圖：T3 Δβ 三柱對照 ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))
bars = [T3["noise_floor_dbeta"], T3["subclone_dbeta"], T3["ctrl_haplotype_dbeta"]]
labs = ["噪音地板\n(同群隨機對半)", "T3 亞群vs母本\n(同 haplotype)", "正控 H1vsH2\n(不同 haplotype, V10)"]
cols = ["#9CA3AF", "#C2410C", "#15803D"]
ax1.bar(range(3), bars, color=cols, alpha=0.85)
for i, v in enumerate(bars):
    ax1.text(i, v + 0.003, f"{v:.3f}", ha="center", fontsize=10)
ax1.set_xticks(range(3)); ax1.set_xticklabels(labs, fontsize=9); ax1.set_ylabel("raw |Δβ| per-CpG (median)")
ax1.set_title(f"T3 嚴格反駁：亞群 Δβ 僅略超噪音(+{T3['subclone_excess_over_noise']})\n遠不及 haplotype 訊號(+{T3['ctrl_excess_over_noise']}) → 太弱不可用")
# T1 failure
g = T1["FAIL_acc<0.7"]; o = T1["OK_acc>=0.9"]
ax2.bar([0, 1], [g.get("train_sep_median", 0) or 0, o.get("train_sep_median", 0) or 0], color=["#C2410C", "#15803D"], alpha=0.85)
ax2.text(0, (g.get("train_sep_median") or 0) + 0.01, f"sep={g.get('train_sep_median')}\nn_cpg={g.get('n_cpg_median')}\nn={g['n']}", ha="center", fontsize=9)
ax2.text(1, (o.get("train_sep_median") or 0) + 0.01, f"sep={o.get('train_sep_median')}\nn_cpg={o.get('n_cpg_median')}\nn={o['n']}", ha="center", fontsize=9)
ax2.set_xticks([0, 1]); ax2.set_xticklabels(["T1 FAIL (acc<0.7)", "T1 OK (acc>=0.9)"], fontsize=9)
ax2.set_ylabel("train_sep median"); ax2.set_ylim(0, 1.15)
ax2.set_title("T1 失敗成因 = 低 train_sep(本質分不開)\n非低覆蓋/CpG（失敗組 n_cpg 不低）")
fig.tight_layout(); fig.savefig(AD + "/fig_rigor.png", dpi=130, bbox_inches="tight"); plt.close()

print("=== rigor 彙總 ===")
print(f"T3: 噪音 Δβ={T3['noise_floor_dbeta']} | 亞群 Δβ={T3['subclone_dbeta']}(+{T3['subclone_excess_over_noise']}) | 正控 Δβ={T3['ctrl_haplotype_dbeta']}(+{T3['ctrl_excess_over_noise']}) (n_win={T3['n_windows']})")
print(f"T1: FAIL n={T1['FAIL_acc<0.7']['n']} sep={T1['FAIL_acc<0.7'].get('train_sep_median')} n_cpg={T1['FAIL_acc<0.7'].get('n_cpg_median')} | OK n={T1['OK_acc>=0.9']['n']} sep={T1['OK_acc>=0.9'].get('train_sep_median')}")
print(f"   by_status: " + " ".join(f"{s}:acc={v.get('acc_median')}(n={v['n']})" for s, v in T1['by_status'].items()))
print(f"T2: H3={T2['n_H3_reads']} class%={T2['H3_class_pct']}")
for k, v in T2["H3_assignability_by_class"].items(): print(f"   {k}: n={v['n_with_centroid']} margin={v['margin_median']} highconf={v['frac_high_conf']}")
print("saved fig_rigor.png + rigor_aggregate.json")
