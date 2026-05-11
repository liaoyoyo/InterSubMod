#!/usr/bin/env python3
"""Step 8: Cross-sample PCA on Coverage_Category distribution (stale vs KDE-fixed).

Context: Original O1-O10 Fig 10 claim "COLO829 isolated in PCA → depth is the
first driver" may be artifact of stale 75× compressing CovM to 0.39 and flooding
CNV_Loss category. Under KDE-fixed baseline, COLO829 CNV_Loss dropped from ~60%
(est) to 3.93%, so its Category vector should no longer be an outlier.

Method:
  - Feature vector per sample = [CNV_Loss%, Low%, Normal%, Elevated%, CNV_Gain%, High_Copy%]
  - PCA on stale (paired mode) → plot PC1/PC2 with sample labels
  - PCA on fixed (paired_pileup + paired_full averaged) → same
  - Compare COLO829 distance to nearest-neighbor (isolation metric)
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIXED = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/kde_rerun_B_14combos/all_region_rows_kde_B_tp.tsv.gz"
STALE = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
         "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz")
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/kde_fix_validation/outputs/step8_pca")
FIG = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260420_kde_fix_acceptance")

SAMPLES = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H2009", "H1437", "COLO829"]
CATS = ["CNV_Loss", "Low", "Normal", "Elevated", "CNV_Gain", "High_Copy"]


def category_vector(df, sample, mode_filter):
    sub = df[(df["sample"] == sample) & df["mode"].isin(mode_filter)]
    dist = sub["Coverage_Category"].value_counts(normalize=True)
    return np.array([100.0 * dist.get(c, 0.0) for c in CATS])


def pca_2d(X):
    Xc = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    pc = Xc @ Vt.T
    var_expl = (S ** 2) / (S ** 2).sum()
    return pc[:, :2], var_expl[:2]


def nn_distance(pc, idx):
    d = np.linalg.norm(pc - pc[idx], axis=1)
    d[idx] = np.inf
    return float(d.min())


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    print("Loading masters ...")
    cols = ["sample", "mode", "Coverage_Category"]
    stale_chunks = []
    for chk in pd.read_csv(STALE, sep="\t", usecols=cols + ["truth_label"],
                           chunksize=200_000, low_memory=False):
        m = (chk["truth_label"] == "TP") & (chk["sample"].isin(SAMPLES)) & (chk["mode"] == "paired")
        if m.any():
            stale_chunks.append(chk[m][cols].copy())
    st = pd.concat(stale_chunks, ignore_index=True)
    fx = pd.read_csv(FIXED, sep="\t", usecols=cols, low_memory=False)
    print(f"Stale TP paired rows: {len(st):,}  Fixed rows: {len(fx):,}")

    X_stale = np.vstack([category_vector(st, s, ["paired"]) for s in SAMPLES])
    X_fixed = np.vstack([category_vector(fx, s, ["paired_pileup", "paired_full"]) for s in SAMPLES])

    pc_s, var_s = pca_2d(X_stale)
    pc_f, var_f = pca_2d(X_fixed)

    colo_idx = SAMPLES.index("COLO829")
    d_colo_stale = nn_distance(pc_s, colo_idx)
    d_colo_fixed = nn_distance(pc_f, colo_idx)

    other_dists_stale = [nn_distance(pc_s, i) for i in range(len(SAMPLES)) if i != colo_idx]
    other_dists_fixed = [nn_distance(pc_f, i) for i in range(len(SAMPLES)) if i != colo_idx]

    isolation_ratio_stale = d_colo_stale / np.median(other_dists_stale)
    isolation_ratio_fixed = d_colo_fixed / np.median(other_dists_fixed)

    # Feature vector table
    rows = []
    for i, s in enumerate(SAMPLES):
        for k, c in enumerate(CATS):
            rows.append({
                "sample": s, "category": c,
                "stale_pct": X_stale[i, k], "fixed_pct": X_fixed[i, k],
                "delta_pp": X_fixed[i, k] - X_stale[i, k],
            })
    pd.DataFrame(rows).to_csv(OUT / "category_vectors.tsv", sep="\t", index=False, float_format="%.3f")

    pca_rows = []
    for i, s in enumerate(SAMPLES):
        pca_rows.append({
            "sample": s,
            "stale_PC1": pc_s[i, 0], "stale_PC2": pc_s[i, 1],
            "fixed_PC1": pc_f[i, 0], "fixed_PC2": pc_f[i, 1],
        })
    pd.DataFrame(pca_rows).to_csv(OUT / "pca_coords.tsv", sep="\t", index=False, float_format="%.3f")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    colors = plt.cm.tab10.colors
    for ax, pc, var, title in [
        (axes[0], pc_s, var_s, f"STALE baseline 75× (COLO829 NN-dist={d_colo_stale:.1f}, ratio={isolation_ratio_stale:.2f}×)"),
        (axes[1], pc_f, var_f, f"KDE-FIXED per-sample (COLO829 NN-dist={d_colo_fixed:.1f}, ratio={isolation_ratio_fixed:.2f}×)"),
    ]:
        for i, s in enumerate(SAMPLES):
            marker = "X" if s == "COLO829" else "o"
            size = 200 if s == "COLO829" else 120
            ax.scatter(pc[i, 0], pc[i, 1], c=[colors[i]], s=size, marker=marker,
                       edgecolors="black", linewidths=0.8, zorder=3)
            ax.annotate(s, (pc[i, 0], pc[i, 1]), fontsize=9,
                        xytext=(5, 5), textcoords="offset points")
        ax.set_xlabel(f"PC1 ({100*var[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({100*var[1]:.1f}%)")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
        ax.axhline(0, color="grey", linewidth=0.5)
        ax.axvline(0, color="grey", linewidth=0.5)
    fig.suptitle("Cross-sample Coverage_Category PCA: stale vs KDE-fixed baseline", y=1.00, fontsize=11)
    fig.tight_layout()
    out_fig = FIG / "fig7_cross_sample_pca.png"
    fig.savefig(out_fig, dpi=130, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_fig}")

    # Summary
    print("\n" + "=" * 80)
    print("COLO829 isolation in Coverage_Category PCA")
    print("=" * 80)
    print(f"  Stale  NN-dist = {d_colo_stale:.2f}  (median of others: {np.median(other_dists_stale):.2f})")
    print(f"  Fixed  NN-dist = {d_colo_fixed:.2f}  (median of others: {np.median(other_dists_fixed):.2f})")
    print(f"  Isolation ratio stale→fixed: {isolation_ratio_stale:.2f}× → {isolation_ratio_fixed:.2f}×")
    if isolation_ratio_fixed < 1.5 * isolation_ratio_stale and isolation_ratio_fixed < 3.0:
        verdict = "COLO829 no longer extreme outlier after KDE fix"
    elif isolation_ratio_fixed < isolation_ratio_stale * 0.8:
        verdict = "COLO829 isolation significantly reduced (artifact was dominant driver)"
    else:
        verdict = "COLO829 isolation persists (not purely baseline artifact)"
    print(f"  Verdict: {verdict}")

    md = OUT / "step8_summary.md"
    with open(md, "w") as f:
        f.write("# Step 8: 跨樣本 Coverage_Category PCA — baseline fix 下游影響\n\n")
        f.write("## 背景\n\n")
        f.write("原 O1-O10 Fig 10 結論：「COLO829 在 PCA 中孤立 → Depth 是第一驅動因素」\n")
        f.write("疑慮：stale 75× 讓 COLO829 CovM 壓縮至 0.39，造成 CNV_Loss 泛濫，可能是 PCA 孤立的 artifact。\n\n")
        f.write("## 特徵向量（per-sample, 6-D Coverage_Category %）\n\n")
        f.write("### COLO829 vector\n\n")
        f.write("| Category | Stale % | Fixed % | Δpp |\n")
        f.write("|----------|--------:|--------:|----:|\n")
        ci = SAMPLES.index("COLO829")
        for k, c in enumerate(CATS):
            f.write(f"| {c} | {X_stale[ci, k]:.2f} | {X_fixed[ci, k]:.2f} | {X_fixed[ci, k]-X_stale[ci, k]:+.2f} |\n")
        f.write("\n### PCA 孤立度對比\n\n")
        f.write("| 指標 | Stale 75× | KDE-fixed | 變化 |\n")
        f.write("|------|----------:|----------:|------|\n")
        f.write(f"| COLO829 NN-dist | {d_colo_stale:.2f} | {d_colo_fixed:.2f} | "
                f"{100*(d_colo_fixed-d_colo_stale)/d_colo_stale:+.1f}% |\n")
        f.write(f"| 其他樣本 NN-dist 中位 | {np.median(other_dists_stale):.2f} | "
                f"{np.median(other_dists_fixed):.2f} | — |\n")
        f.write(f"| COLO829 isolation ratio | {isolation_ratio_stale:.2f}× | "
                f"{isolation_ratio_fixed:.2f}× | — |\n")
        f.write(f"| PC1/PC2 var explained (stale) | {100*var_s[0]:.1f}% + {100*var_s[1]:.1f}% | — | — |\n")
        f.write(f"| PC1/PC2 var explained (fixed) | — | {100*var_f[0]:.1f}% + {100*var_f[1]:.1f}% | — |\n")
        f.write(f"\n**Verdict**: {verdict}\n\n")
        f.write("## 結論\n\n")
        if isolation_ratio_fixed < isolation_ratio_stale * 0.8:
            f.write("- 「Depth 是第一驅動」屬於 stale baseline artifact\n")
            f.write("- KDE fix 後 COLO829 Category 分佈已向其他樣本靠攏\n")
            f.write("- 原 O1-O10 Fig 10 結論需**降級**：從「depth driver」改為「stale baseline artifact」\n")
        elif isolation_ratio_fixed < 3.0:
            f.write("- COLO829 孤立度大幅降低；原結論部分成立但需修正\n")
            f.write("- Depth 仍是一個因素，但 baseline artifact 貢獻顯著\n")
        else:
            f.write("- COLO829 仍相對孤立；baseline artifact 並非唯一原因\n")
            f.write("- 原「depth driver」結論仍成立（需額外說明 baseline 影響）\n")

    print(f"Wrote {md}")


if __name__ == "__main__":
    main()
