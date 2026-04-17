"""
B.1-2 質疑驗證：HPFineNGroups 飽和效應檢查

核心問題：N≥4+NR≥80 TP rate 89.1% 是 NGroups 訊號還是 NR 驅動的飽和？

驗證步驟：
  S1: 讀取 748K 資料並拆分 mode × LOH
  S2: NR bin × NGroups 雙層交叉 TP rate 表
  S3: NR≥80 內 NGroups=3 vs 4 TP rate 差異（核心）
  S4: NR 匹配 Fisher exact 檢驗
  S5: Paired mode 並行複驗
  S6: 7 樣本分別 NR≥80 內 NGroups=3 vs 4

執行：python scripts/01_saturation_check.py
輸出：data/*.tsv, figures/*.png
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[3]
DATA_TSV = ROOT / "output" / "synthesis" / "observation_workspaces" / \
           "20260327_loh_round1_cross_sample_audit" / "all_region_rows.tsv.gz"
OUT = Path(__file__).resolve().parents[1]
DATA_OUT = OUT / "data"
FIG_OUT = OUT / "figures"

NR_BINS = [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100), (100, 150), (150, 500)]
NR_LABELS = [f"{a}-{b}" for a, b in NR_BINS]
NGROUPS_VALUES = [1, 2, 3, 4]
MIN_N = 30
ALPHA = 0.05


def nr_bin(nr: int) -> str | None:
    for (a, b), label in zip(NR_BINS, NR_LABELS):
        if a <= nr < b:
            return label
    return None


def load_data() -> pd.DataFrame:
    cols = ["sample", "mode", "truth_label", "NumReads", "HPFineNGroups",
            "to_loh_bed_hit", "caller_af", "HP_Ratio", "Potential_LOH",
            "Coverage_Multiple"]
    print(f"[S1] Loading {DATA_TSV} ...", flush=True)
    df = pd.read_csv(DATA_TSV, sep="\t", usecols=cols,
                     dtype={"sample": "category", "mode": "category",
                            "truth_label": "category",
                            "to_loh_bed_hit": "category"})
    print(f"[S1] Loaded rows={len(df):,}", flush=True)
    # Normalize
    df["is_tp"] = (df["truth_label"] == "TP").astype(int)
    df["is_loh"] = df["to_loh_bed_hit"].astype(str).str.lower().isin(["true", "1"])
    df["nr_bin"] = df["NumReads"].apply(nr_bin)
    df = df.dropna(subset=["nr_bin"])
    df["HPFineNGroups"] = df["HPFineNGroups"].astype(int)
    return df


def s2_crosstab(df: pd.DataFrame, mode: str, loh: bool) -> pd.DataFrame:
    sub = df[(df["mode"] == mode) & (df["is_loh"] == loh)]
    rows = []
    for nr_label in NR_LABELS:
        for ng in NGROUPS_VALUES:
            bin_df = sub[(sub["nr_bin"] == nr_label) & (sub["HPFineNGroups"] == ng)]
            n = len(bin_df)
            tp = bin_df["is_tp"].sum()
            rate = tp / n if n > 0 else np.nan
            rows.append({
                "mode": mode,
                "LOH": loh,
                "nr_bin": nr_label,
                "HPFineNGroups": ng,
                "n": n,
                "tp": int(tp),
                "fp": n - int(tp),
                "tp_rate": rate,
            })
    return pd.DataFrame(rows)


def s3_nr80_saturation(df: pd.DataFrame, mode: str, loh: bool) -> dict:
    """NR≥80 內 NGroups=3 vs 4 的 Fisher exact 檢驗。"""
    sub = df[(df["mode"] == mode) & (df["is_loh"] == loh) & (df["NumReads"] >= 80)]
    n3 = sub[sub["HPFineNGroups"] == 3]
    n4 = sub[sub["HPFineNGroups"] == 4]
    tp3, fp3 = int(n3["is_tp"].sum()), len(n3) - int(n3["is_tp"].sum())
    tp4, fp4 = int(n4["is_tp"].sum()), len(n4) - int(n4["is_tp"].sum())
    rate3 = tp3 / len(n3) if len(n3) else np.nan
    rate4 = tp4 / len(n4) if len(n4) else np.nan
    odds, p = fisher_exact([[tp4, fp4], [tp3, fp3]], alternative="greater")
    return {
        "mode": mode,
        "LOH": loh,
        "filter": "NR>=80",
        "n3": len(n3), "tp3": tp3, "rate3": rate3,
        "n4": len(n4), "tp4": tp4, "rate4": rate4,
        "delta_pp": (rate4 - rate3) * 100 if not np.isnan(rate4 - rate3) else np.nan,
        "odds_ratio": odds,
        "fisher_p_one_sided": p,
    }


def s4_matched_nr(df: pd.DataFrame, mode: str, loh: bool) -> pd.DataFrame:
    """每個 NR bin 內 NGroups=3 vs 4 的 Fisher exact。"""
    rows = []
    sub = df[(df["mode"] == mode) & (df["is_loh"] == loh)]
    for nr_label in NR_LABELS:
        bin_df = sub[sub["nr_bin"] == nr_label]
        n3 = bin_df[bin_df["HPFineNGroups"] == 3]
        n4 = bin_df[bin_df["HPFineNGroups"] == 4]
        if len(n3) < MIN_N or len(n4) < MIN_N:
            continue
        tp3, fp3 = int(n3["is_tp"].sum()), len(n3) - int(n3["is_tp"].sum())
        tp4, fp4 = int(n4["is_tp"].sum()), len(n4) - int(n4["is_tp"].sum())
        rate3 = tp3 / len(n3)
        rate4 = tp4 / len(n4)
        odds, p = fisher_exact([[tp4, fp4], [tp3, fp3]], alternative="greater")
        rows.append({
            "mode": mode, "LOH": loh, "nr_bin": nr_label,
            "n3": len(n3), "rate3": rate3,
            "n4": len(n4), "rate4": rate4,
            "delta_pp": (rate4 - rate3) * 100,
            "odds_ratio": odds,
            "fisher_p_one_sided": p,
            "signif": "*" if p < ALPHA else "",
        })
    return pd.DataFrame(rows)


def s6_per_sample(df: pd.DataFrame, mode: str, loh: bool) -> pd.DataFrame:
    """每個樣本 NR≥80 內 NGroups=3 vs 4。"""
    rows = []
    samples = df["sample"].cat.categories if hasattr(df["sample"], "cat") else df["sample"].unique()
    for sample in samples:
        sub = df[(df["sample"] == sample) & (df["mode"] == mode) &
                 (df["is_loh"] == loh) & (df["NumReads"] >= 80)]
        n3 = sub[sub["HPFineNGroups"] == 3]
        n4 = sub[sub["HPFineNGroups"] == 4]
        if len(n3) < 10 or len(n4) < 10:
            rows.append({
                "sample": sample, "mode": mode, "LOH": loh,
                "n3": len(n3), "n4": len(n4),
                "rate3": np.nan, "rate4": np.nan,
                "delta_pp": np.nan, "fisher_p": np.nan,
                "direction": "INSUFFICIENT",
            })
            continue
        tp3, fp3 = int(n3["is_tp"].sum()), len(n3) - int(n3["is_tp"].sum())
        tp4, fp4 = int(n4["is_tp"].sum()), len(n4) - int(n4["is_tp"].sum())
        rate3 = tp3 / len(n3); rate4 = tp4 / len(n4)
        _, p = fisher_exact([[tp4, fp4], [tp3, fp3]], alternative="greater")
        rows.append({
            "sample": sample, "mode": mode, "LOH": loh,
            "n3": len(n3), "rate3": rate3,
            "n4": len(n4), "rate4": rate4,
            "delta_pp": (rate4 - rate3) * 100,
            "fisher_p": p,
            "direction": "POS" if rate4 > rate3 else "NEG" if rate4 < rate3 else "FLAT",
        })
    return pd.DataFrame(rows)


def plot_heatmap(df_crosstab: pd.DataFrame, out_path: Path) -> None:
    """NR × NGroups TP rate heatmap (TO NonLOH + Paired NonLOH + TO LOH 三聯)."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    configs = [("to", False, "TO NonLOH"),
               ("to", True, "TO LOH"),
               ("paired", False, "Paired NonLOH")]
    for ax, (mode, loh, title) in zip(axes, configs):
        sub = df_crosstab[(df_crosstab["mode"] == mode) & (df_crosstab["LOH"] == loh)]
        pivot_rate = sub.pivot(index="nr_bin", columns="HPFineNGroups", values="tp_rate")
        pivot_n = sub.pivot(index="nr_bin", columns="HPFineNGroups", values="n")
        pivot_rate = pivot_rate.reindex(NR_LABELS)
        pivot_n = pivot_n.reindex(NR_LABELS)
        im = ax.imshow(pivot_rate.values, cmap="RdYlGn", vmin=0.3, vmax=1.0, aspect="auto")
        for i in range(pivot_rate.shape[0]):
            for j in range(pivot_rate.shape[1]):
                rate = pivot_rate.values[i, j]
                n = pivot_n.values[i, j]
                if pd.notna(rate) and n >= MIN_N:
                    ax.text(j, i, f"{rate:.2f}\nn={int(n)}",
                            ha="center", va="center", fontsize=8,
                            color="black" if rate > 0.6 else "white")
        ax.set_xticks(range(len(NGROUPS_VALUES)))
        ax.set_xticklabels([f"N={n}" for n in NGROUPS_VALUES])
        ax.set_yticks(range(len(NR_LABELS)))
        ax.set_yticklabels(NR_LABELS)
        ax.set_xlabel("HPFineNGroups")
        ax.set_ylabel("NumReads bin")
        ax.set_title(title)
        plt.colorbar(im, ax=ax, label="TP rate")
    plt.suptitle("S2: TP rate heatmap (NR × NGroups)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def plot_per_sample_forest(per_sample_df: pd.DataFrame, out_path: Path) -> None:
    sub = per_sample_df[per_sample_df["direction"] != "INSUFFICIENT"]
    if sub.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, mode in zip(axes, ["to", "paired"]):
        mode_sub = sub[sub["mode"] == mode]
        y = np.arange(len(mode_sub))
        colors = ["green" if r > 0 else "red" for r in mode_sub["delta_pp"]]
        ax.barh(y, mode_sub["delta_pp"], color=colors, alpha=0.7)
        ax.axvline(0, color="black", lw=0.8)
        ax.axvline(5, color="orange", ls="--", lw=0.8, label="+5pp threshold")
        ax.axvline(10, color="red", ls="--", lw=0.8, label="+10pp threshold")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{s} (n3={n3},n4={n4})"
                            for s, n3, n4 in zip(mode_sub["sample"], mode_sub["n3"], mode_sub["n4"])])
        ax.set_xlabel("Δ TP rate (N4-N3) in pp")
        ax.set_title(f"S6: {mode} NonLOH NR≥80, per-sample N4-N3")
        ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    DATA_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    df = load_data()

    # S1 分佈摘要
    print("\n[S1] Distribution summary:")
    summary = df.groupby(["mode", "is_loh", "truth_label"], observed=True).size()
    print(summary)

    # S2
    print("\n[S2] NR × NGroups crosstab ...", flush=True)
    crosstab_rows = []
    for mode in ["to", "paired"]:
        for loh in [False, True]:
            crosstab_rows.append(s2_crosstab(df, mode, loh))
    crosstab = pd.concat(crosstab_rows, ignore_index=True)
    crosstab.to_csv(DATA_OUT / "nr_ngroups_crosstab.tsv", sep="\t", index=False)

    # S3
    print("[S3] NR≥80 saturation test ...", flush=True)
    s3_rows = []
    for mode in ["to", "paired"]:
        for loh in [False, True]:
            s3_rows.append(s3_nr80_saturation(df, mode, loh))
    s3_df = pd.DataFrame(s3_rows)
    s3_df.to_csv(DATA_OUT / "s3_nr80_saturation.tsv", sep="\t", index=False)
    print(s3_df.to_string(index=False))

    # S4 matched
    print("\n[S4] NR-matched Fisher exact ...", flush=True)
    s4_rows = []
    for mode in ["to", "paired"]:
        for loh in [False, True]:
            df_matched = s4_matched_nr(df, mode, loh)
            s4_rows.append(df_matched)
    s4 = pd.concat(s4_rows, ignore_index=True)
    s4.to_csv(DATA_OUT / "nr_matched_test.tsv", sep="\t", index=False)
    print(s4.to_string(index=False))

    # S6
    print("\n[S6] Per-sample NR≥80 ...", flush=True)
    s6_rows = []
    for mode in ["to", "paired"]:
        for loh in [False, True]:
            s6_rows.append(s6_per_sample(df, mode, loh))
    s6 = pd.concat(s6_rows, ignore_index=True)
    s6.to_csv(DATA_OUT / "per_sample_nr80_ngroups.tsv", sep="\t", index=False)
    print(s6.to_string(index=False))

    # Figures
    print("\n[Fig] heatmap ...", flush=True)
    plot_heatmap(crosstab, FIG_OUT / "01_nr_ngroups_tp_rate_heatmap.png")
    plot_per_sample_forest(s6, FIG_OUT / "02_nr_matched_delta_forest.png")

    # Final verdict
    print("\n=== FINAL VERDICT (TO NonLOH NR≥80) ===")
    to_nonloh = s3_df[(s3_df["mode"] == "to") & (s3_df["LOH"] == False)].iloc[0]
    delta = to_nonloh["delta_pp"]
    p_val = to_nonloh["fisher_p_one_sided"]
    print(f"rate3={to_nonloh['rate3']:.4f} (n={to_nonloh['n3']})")
    print(f"rate4={to_nonloh['rate4']:.4f} (n={to_nonloh['n4']})")
    print(f"delta_pp={delta:.2f} pp, fisher_p={p_val:.3e}")
    if delta < 5:
        print("VERDICT: SATURATION CONFIRMED — 結論翻轉為 NEGATIVE")
    elif delta < 10:
        print("VERDICT: WEAK MARKER — 結論降級為 ⭐2")
    else:
        print("VERDICT: TRUE SUBCLONE MARKER — 結論保留 ⭐3")

    return 0


if __name__ == "__main__":
    sys.exit(main())
