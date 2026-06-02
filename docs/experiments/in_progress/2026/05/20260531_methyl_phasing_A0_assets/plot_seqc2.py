#!/usr/bin/env python3
"""畫 SEQC2 × 甲基 4 張圖（數據全來自 seqc2_cn_methyl_chr*.json，不產生新統計）。"""
import json, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK TC", "Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

AD = os.path.dirname(os.path.abspath(__file__))
rows = []
for fp in glob.glob(AD + "/seqc2_cn_methyl_chr*.json"):
    rows.extend(json.load(open(fp))["regions"])

C = {"gain": "#C2410C", "loss": "#1E3A8A", "loh": "#A16207", "neutral": "#6B6B66"}

# ── 圖1: anchor AUC by status (boxplot) ──
fig, ax = plt.subplots(figsize=(7, 4.5))
order = ["neutral", "gain", "loss", "loh"]
data = [[r["anchor_auc"] for r in rows if r["seqc2_status"] == s] for s in order]
bp = ax.boxplot(data, labels=[f"{s}\n(n={len(d)})" for s, d in zip(order, data)],
                patch_artist=True, widths=0.6)
for patch, s in zip(bp["boxes"], order):
    patch.set_facecolor(C[s]); patch.set_alpha(0.6)
ax.axhline(0.58, ls="--", c="gray", lw=1, label="AUC=0.58 門檻")
ax.set_ylabel("甲基 anchor AUC (HP1 vs HP2)")
ax.set_title("圖1: 甲基分離 AUC × SEQC2 CN 狀態\nLOH 區顯著較低 (p=0.0019)；gain≈neutral (p=0.70)")
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(AD + "/fig_seqc2_auc_by_status.png", dpi=130); plt.close()

# ── 圖2: CN 整數 vs AUC (scatter, 看是否相關) ──
fig, ax = plt.subplots(figsize=(7, 4.5))
cn_rows = [r for r in rows if isinstance(r["seqc2_cn"], (int, float)) and float(r["seqc2_cn"]).is_integer()]
xs = [int(r["seqc2_cn"]) + np.random.RandomState(int(r["pos"]) % 1000).uniform(-0.15, 0.15) for r in cn_rows]
ys = [r["anchor_auc"] for r in cn_rows]
ax.scatter(xs, ys, c="#D97757", alpha=0.5, s=22)
ax.set_xlabel("SEQC2 copy number (整數)")
ax.set_ylabel("甲基 anchor AUC")
ax.set_title("圖2: 倍體 vs 甲基分離 — 無相關 (Spearman rho=0.035, p=0.71)\n→ 甲基分得開『不是』因為 CN 高 (否證 CN-confound 假設)")
fig.tight_layout(); fig.savefig(AD + "/fig_seqc2_cn_vs_auc.png", dpi=130); plt.close()

# ── 圖3: CN vs depth (CN proxy 有效性, Q3) ──
fig, ax = plt.subplots(figsize=(7, 4.5))
dc = [(int(r["seqc2_cn"]), r["mean_depth"]) for r in cn_rows if r["mean_depth"] is not None]
cns = sorted(set(c for c, _ in dc))
bydc = [[d for c, d in dc if c == cn] for cn in cns]
bp = ax.boxplot(bydc, labels=[f"CN{cn}\n(n={len(b)})" for cn, b in zip(cns, bydc)], patch_artist=True, widths=0.6)
for patch in bp["boxes"]:
    patch.set_facecolor("#15803D"); patch.set_alpha(0.5)
ax.set_xlabel("SEQC2 copy number"); ax.set_ylabel("本 BAM 平均 coverage")
ax.set_title("圖3: copy number vs 實測 coverage — 強單調 (Spearman rho=0.923, p<1e-50)\n→ coverage 是有效 CN proxy；甲基不需要也能由 depth 讀倍體")
fig.tight_layout(); fig.savefig(AD + "/fig_seqc2_cn_vs_depth.png", dpi=130); plt.close()

# ── 圖4: bimodality (GMM BIC) by status (Q2) ──
fig, ax = plt.subplots(figsize=(7, 4.5))
data2 = [[r["gmm_bic_diff"] for r in rows if r["seqc2_status"] == s and r["gmm_bic_diff"] is not None] for s in order]
bp = ax.boxplot(data2, labels=[f"{s}\n(n={len(d)})" for s, d in zip(order, data2)], patch_artist=True, widths=0.6)
for patch, s in zip(bp["boxes"], order):
    patch.set_facecolor(C[s]); patch.set_alpha(0.6)
ax.axhline(10, ls="--", c="gray", lw=1, label="BIC diff=10 (多模態門檻)")
ax.set_ylabel("GMM BIC(1) − BIC(2)  (>0 = 甲基雙峰)")
ax.set_title("圖4: 甲基多模態 × CN 狀態 (Q2)\ngain/loss 較常雙峰 (~50%) vs neutral (21%)")
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(AD + "/fig_seqc2_bimodality.png", dpi=130); plt.close()

print("4 figs saved:")
for f in ["fig_seqc2_auc_by_status.png", "fig_seqc2_cn_vs_auc.png", "fig_seqc2_cn_vs_depth.png", "fig_seqc2_bimodality.png"]:
    p = AD + "/" + f
    print(f"  {f}: {'OK' if os.path.exists(p) else 'MISSING'} ({os.path.getsize(p)//1024 if os.path.exists(p) else 0}KB)")
