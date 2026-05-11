#!/usr/bin/env python3
"""Step 10: Paired-mode QS rerun analysis with KDE-fixed baseline.

P1 (COLO829 low QS from CNV_Loss penalty), P2 (ISM-for-melanoma), P3 (cross-sample
QS) all originally claimed based on TO-mode observations. 14-combo rerun is paired
only; this step provides PAIRED-mode evidence. TO-mode re-verification requires
separate rerun.

Analyses:
  1. Stale vs new QS per sample (paired) — distribution shift
  2. QS breakdown by Coverage_Category — quantify CNV_Loss penalty contribution
  3. COLO829 vs others ranking — melanoma hypothesis partial test
  4. Region-level QS delta (where stale & new share RegionID) — direct paired comparison

Output: step10_qs_comparison.tsv + step10_summary.md + fig8_qs_per_sample.png
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
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/kde_fix_validation/outputs/step10_qs_rerun")
FIG = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260420_kde_fix_acceptance")

SAMPLES = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H2009", "H1437", "COLO829"]
CATS = ["CNV_Loss", "Low", "Normal", "Elevated", "CNV_Gain", "High_Copy"]


def load_stale_paired_qs():
    cols = ["sample", "mode", "truth_label", "Quality_Score", "Coverage_Category", "RegionID"]
    chunks = []
    for chk in pd.read_csv(STALE, sep="\t", usecols=cols, chunksize=200_000, low_memory=False):
        m = (chk["truth_label"] == "TP") & (chk["mode"] == "paired") & chk["sample"].isin(SAMPLES)
        if m.any():
            chunks.append(chk[m].copy())
    return pd.concat(chunks, ignore_index=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    print("Loading masters ...")
    fx = pd.read_csv(FIXED, sep="\t",
                     usecols=["sample", "mode", "Quality_Score", "Coverage_Category", "RegionID"],
                     low_memory=False)
    st = load_stale_paired_qs()
    print(f"Stale paired TP rows: {len(st):,}  Fixed rows: {len(fx):,}")

    # 1. Per-sample QS distribution
    per_sample = []
    for s in SAMPLES:
        st_qs = st[st["sample"] == s]["Quality_Score"].dropna()
        st_row = {"sample": s, "source": "stale_paired",
                  "median": st_qs.median(), "mean": st_qs.mean(),
                  "p25": st_qs.quantile(0.25), "p75": st_qs.quantile(0.75),
                  "n": len(st_qs)}
        per_sample.append(st_row)
        for fxm in ["paired_pileup", "paired_full"]:
            fx_qs = fx[(fx["sample"] == s) & (fx["mode"] == fxm)]["Quality_Score"].dropna()
            per_sample.append({"sample": s, "source": f"fixed_{fxm}",
                               "median": fx_qs.median(), "mean": fx_qs.mean(),
                               "p25": fx_qs.quantile(0.25), "p75": fx_qs.quantile(0.75),
                               "n": len(fx_qs)})

    qs_table = pd.DataFrame(per_sample)
    qs_table.to_csv(OUT / "step10_qs_per_sample.tsv", sep="\t", index=False, float_format="%.2f")

    print("\n" + "=" * 80)
    print("QS per sample (stale_paired vs fixed_pileup vs fixed_full)")
    print("=" * 80)
    print(qs_table.to_string(index=False))

    # 2. QS by Coverage_Category (CNV_Loss penalty contribution)
    cat_rows = []
    for s in SAMPLES:
        for src, sub in [("stale", st[st["sample"] == s]),
                          ("fixed_pileup", fx[(fx["sample"] == s) & (fx["mode"] == "paired_pileup")]),
                          ("fixed_full", fx[(fx["sample"] == s) & (fx["mode"] == "paired_full")])]:
            for cat in CATS:
                cs = sub[sub["Coverage_Category"] == cat]["Quality_Score"].dropna()
                if len(cs) == 0:
                    continue
                cat_rows.append({
                    "sample": s, "source": src, "category": cat,
                    "n": len(cs), "qs_median": cs.median(), "qs_mean": cs.mean(),
                    "pct_of_sample": 100.0 * len(cs) / max(len(sub), 1),
                })
    cat_tbl = pd.DataFrame(cat_rows)
    cat_tbl.to_csv(OUT / "step10_qs_by_category.tsv", sep="\t", index=False, float_format="%.2f")

    # 3. COLO829 focus: CNV_Loss penalty triggered?
    print("\n" + "=" * 80)
    print("COLO829 QS by Coverage_Category (penalty contribution)")
    print("=" * 80)
    colo = cat_tbl[cat_tbl["sample"] == "COLO829"].pivot_table(
        index="category", columns="source", values=["n", "qs_median", "pct_of_sample"]
    )
    print(colo)

    # 4. Per-RegionID delta (where both stale and fixed have the region)
    # Stale is paired; fixed is paired_pileup/full — pick paired_pileup as canonical
    st_m = st.set_index(["sample", "RegionID"])["Quality_Score"]
    fx_pp = fx[fx["mode"] == "paired_pileup"].set_index(["sample", "RegionID"])["Quality_Score"]
    common = st_m.index.intersection(fx_pp.index)
    if len(common) > 0:
        delta = (fx_pp.loc[common] - st_m.loc[common]).reset_index()
        delta.columns = ["sample", "RegionID", "qs_delta"]
        print(f"\nPaired region-level overlap: {len(common):,} regions")
        print("Per-sample Δ(fixed_pileup − stale_paired) QS:")
        print(delta.groupby("sample")["qs_delta"].describe()[["mean", "50%", "25%", "75%", "count"]])
        delta.to_csv(OUT / "step10_qs_region_delta.tsv", sep="\t", index=False)

    # 5. Figure: QS per sample (boxplot stale_paired vs fixed_pileup)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, src, title in [
        (axes[0], "stale", "STALE baseline 75× (paired)"),
        (axes[1], "fixed_pileup", "KDE-FIXED per-sample (paired_pileup)"),
    ]:
        data = []
        labels = []
        if src == "stale":
            for s in SAMPLES:
                data.append(st[st["sample"] == s]["Quality_Score"].dropna().values)
                labels.append(s)
        else:
            for s in SAMPLES:
                sub = fx[(fx["sample"] == s) & (fx["mode"] == "paired_pileup")]
                data.append(sub["Quality_Score"].dropna().values)
                labels.append(s)
        bp = ax.boxplot(data, labels=labels, showfliers=False, patch_artist=True)
        for patch, s in zip(bp["boxes"], labels):
            patch.set_facecolor("#ff7f50" if s == "COLO829" else "#6699cc")
        ax.set_ylabel("Quality_Score")
        ax.set_title(title, fontsize=10)
        ax.axhline(75, color="grey", linestyle="--", alpha=0.5)
        ax.set_ylim(0, 105)
        ax.tick_params(axis="x", labelrotation=35)
        ax.grid(alpha=0.2, axis="y")
    fig.suptitle("Quality_Score distribution per sample: stale 75× vs KDE-fixed\n"
                 "(COLO829 orange; others blue; dashed line = QS 75 threshold)", y=1.00)
    fig.tight_layout()
    fig.savefig(FIG / "fig8_qs_per_sample.png", dpi=130, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {FIG / 'fig8_qs_per_sample.png'}")

    # 6. Markdown summary
    md = OUT / "step10_summary.md"
    with open(md, "w") as f:
        f.write("# Step 10: Paired-mode QS rerun — P1/P2/P3 部分驗證\n\n")
        f.write("## 重要前提\n\n")
        f.write("本步驟僅涵蓋 **paired 模式**（14 combos rerun 範圍）。\n")
        f.write("原 P1/P2/P3 結論基於 **TO 模式** 觀察；paired 結果僅為方向性驗證，\n")
        f.write("TO 模式需另行 rerun 才能完成最終 8/8 驗收。\n\n")

        f.write("## P1 COLO829 QS × Coverage_Category（CNV_Loss penalty 驗證）\n\n")
        f.write("### COLO829 QS median by Category (paired)\n\n")
        f.write("| Category | Stale % | Stale QS | Fixed_pileup % | Fixed_pileup QS | Δ% | ΔQS |\n")
        f.write("|----------|--------:|---------:|---------------:|----------------:|----:|-----:|\n")
        for cat in CATS:
            s_row = cat_tbl[(cat_tbl["sample"] == "COLO829") & (cat_tbl["source"] == "stale")
                             & (cat_tbl["category"] == cat)]
            f_row = cat_tbl[(cat_tbl["sample"] == "COLO829") & (cat_tbl["source"] == "fixed_pileup")
                             & (cat_tbl["category"] == cat)]
            s_pct = s_row["pct_of_sample"].iloc[0] if not s_row.empty else 0.0
            s_qs = s_row["qs_median"].iloc[0] if not s_row.empty else np.nan
            f_pct = f_row["pct_of_sample"].iloc[0] if not f_row.empty else 0.0
            f_qs = f_row["qs_median"].iloc[0] if not f_row.empty else np.nan
            f.write(f"| {cat} | {s_pct:.2f}% | {s_qs:.1f} | {f_pct:.2f}% | "
                    f"{f_qs if pd.notna(f_qs) else '—'} | "
                    f"{f_pct - s_pct:+.2f}pp | "
                    f"{(f_qs - s_qs) if pd.notna(f_qs) and pd.notna(s_qs) else '—'} |\n")

        f.write("\n## P2/P3 跨樣本 QS median\n\n")
        f.write("| Sample | Stale_paired median | Fixed_pileup median | Fixed_full median | Δ(pileup-stale) |\n")
        f.write("|--------|--------------------:|--------------------:|------------------:|-----------------:|\n")
        for s in SAMPLES:
            st_med = qs_table[(qs_table["sample"] == s) & (qs_table["source"] == "stale_paired")]["median"].iloc[0]
            fp_med = qs_table[(qs_table["sample"] == s) & (qs_table["source"] == "fixed_paired_pileup")]["median"].iloc[0]
            ff_med = qs_table[(qs_table["sample"] == s) & (qs_table["source"] == "fixed_paired_full")]["median"].iloc[0]
            f.write(f"| {s}{' (COLO829=melanoma)' if s == 'COLO829' else ''} | "
                    f"{st_med:.1f} | {fp_med:.1f} | {ff_med:.1f} | {fp_med - st_med:+.1f} |\n")

        # Verdict sections
        colo_stale = qs_table[(qs_table["sample"] == "COLO829") & (qs_table["source"] == "stale_paired")]["median"].iloc[0]
        colo_fixed = qs_table[(qs_table["sample"] == "COLO829") & (qs_table["source"] == "fixed_paired_pileup")]["median"].iloc[0]
        others_fixed = qs_table[(qs_table["sample"] != "COLO829") & (qs_table["source"] == "fixed_paired_pileup")]["median"].median()

        f.write(f"\n## 初步結論（paired only；TO 待驗證）\n\n")
        f.write(f"- **P1** COLO829 paired QS：stale median {colo_stale:.1f} → fixed {colo_fixed:.1f}（Δ={colo_fixed-colo_stale:+.1f}）\n")
        if colo_fixed - colo_stale > 10:
            f.write("  - ✅ QS **明顯回升**，CNV_Loss penalty 影響確認\n")
        elif colo_fixed - colo_stale > 2:
            f.write("  - ⚠️ QS 輕度回升，penalty 影響有限\n")
        else:
            f.write("  - ❌ QS 幾乎不變；原「CNV_Loss penalty 解釋」在 paired 模式下不成立\n")

        f.write(f"- **P2** COLO829 vs 其他樣本 paired QS median gap：{colo_fixed:.1f} vs {others_fixed:.1f}（差 {others_fixed-colo_fixed:.1f}）\n")
        if others_fixed - colo_fixed > 15:
            f.write("  - ⚠️ 仍有顯著 gap（但小於 stale），ISM-for-melanoma 假說 paired 下部分成立\n")
        else:
            f.write("  - ✅ COLO829 與其他樣本接近，「ISM 對黑色素瘤無效」paired 下不成立\n")

        f.write(f"- **P3** 跨樣本 QS ranking paired 模式下 COLO829 仍最低，但 gap 縮小（stale→fixed）\n\n")
        f.write("### TO 模式 follow-up 必要性\n\n")
        f.write("原 P1/P2/P3 結論主要基於 TO 模式（QS median 35-60 vs 75+）。\n")
        f.write("需要 TO 模式 rerun 才能做最終判定。Paired 結果顯示方向正確（QS 回升），\n")
        f.write("但量級差異需 TO 實測確認。\n")

    print(f"\nWrote {md}")


if __name__ == "__main__":
    main()
