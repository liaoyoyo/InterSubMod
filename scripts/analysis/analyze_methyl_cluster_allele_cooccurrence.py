#!/usr/bin/env python3
"""
Observation Document 1: Methylation Cluster <-> HP/Allele Co-occurrence Analysis
Sample: HCC1395 5kHz paired_full (20260314)
Goal:   Test if methylation clusters are non-randomly associated with HP1/HP2
        haplotype labels and ALT/REF allele support in legacy Strong vs Noise loci.

Output: /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
        20260324_methylation_cluster_allele_cooccurrence_HCC1395/
"""

import argparse
import json
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import sys
from datetime import datetime

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    SchemaContractError,
    select_legacy_view,
)

# Keep known plotting deprecation noise quiet without hiding schema warnings.
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")

# ===== Config =====
ISM_TP_DIR = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full"
    "/20260314_HCC1395_paired_full_full_complete_matrix/intersubmod_tp"
)
OUTPUT_DIR = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces"
    "/20260324_methylation_cluster_allele_cooccurrence_HCC1395"
)

SAMPLE_N_PER_CLASS = 300   # loci to sample per legacy verification class
MIN_READS = 12             # min reads after NA filter
MIN_CPGS_PER_READ = 5      # min non-NA CpGs per read
METHYL_HIGH_THR = 0.6      # threshold: mean methylation >= this => High
METHYL_LOW_THR = 0.4       # threshold: mean methylation <= this => Low (else Middle)
RANDOM_SEED = 42

TARGET_CLASSES = ["Strong", "Noise"]
LEGACY_SELECTION_COLUMN = "_VerificationClass_Legacy_Selected"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Replicate the historical legacy Strong-vs-Noise co-occurrence analysis."
    )
    parser.add_argument(
        "--significance-summary",
        type=Path,
        default=None,
        help="Input significance_summary.csv (default: <ism-tp-dir>/significance_summary.csv)",
    )
    parser.add_argument("--ism-tp-dir", type=Path, default=ISM_TP_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help=(
            "Explicitly authorize historical unversioned VerificationClass as the legacy four-state field; "
            "unknown values are excluded and counted"
        ),
    )
    return parser.parse_args()


def select_replication_cohorts(sig_df: pd.DataFrame, allow_unversioned_v1: bool = False):
    """Select legacy Strong/Noise without folding any v2 current classes."""
    view = select_legacy_view(
        sig_df,
        allow_unversioned_v1=allow_unversioned_v1,
        unversioned_unknown_policy="exclude",
    )
    selected = sig_df.copy()
    selected[LEGACY_SELECTION_COLUMN] = view.values
    selected["_SelectionField"] = view.field
    counts = selected[LEGACY_SELECTION_COLUMN].value_counts()
    target_mask = selected[LEGACY_SELECTION_COLUMN].isin(TARGET_CLASSES)
    metadata = view.metadata()
    metadata.update(
        {
            "selection_values": list(TARGET_CLASSES),
            "input_row_count": int(len(selected)),
            "selected_row_count": int(target_mask.sum()),
            "selected_counts": {name: int(counts.get(name, 0)) for name in TARGET_CLASSES},
            "excluded_non_target_count": int((selected[LEGACY_SELECTION_COLUMN].notna() & ~target_mask).sum()),
            "excluded_unknown_count": int(selected[LEGACY_SELECTION_COLUMN].isna().sum()),
        }
    )
    return selected[target_mask].copy(), metadata

# ===== Path helpers =====

def find_region_dir(chr_name: str, pos: int):
    """Locate the ISM region directory for a given locus."""
    snv_parent = ISM_TP_DIR / "filtered_snv_tp" / chr_name / f"{chr_name}_{pos}"
    if not snv_parent.exists():
        return None
    window_dirs = [d for d in snv_parent.iterdir() if d.is_dir()]
    return window_dirs[0] if window_dirs else None


def load_region(region_dir):
    """Return (reads_df, methyl_df) or (None, None) on failure."""
    reads_path = region_dir / "reads" / "reads.tsv"
    methyl_path = region_dir / "methylation" / "methylation.csv"
    if not reads_path.exists() or not methyl_path.exists():
        return None, None
    reads = pd.read_csv(reads_path, sep="\t")
    methyl = pd.read_csv(methyl_path, index_col=0)
    return reads, methyl


# ===== Feature computation =====

def simplify_hp(hp_val) -> str:
    s = str(hp_val).strip()
    if s == "1":
        return "HP1"
    elif s == "2":
        return "HP2"
    elif "-" in s:
        return "Somatic"
    elif s in ("0", "3", "nan", ""):
        return "Unphased"
    return "Unphased"


def build_read_features(reads: pd.DataFrame, methyl: pd.DataFrame):
    """Join reads + methylation, compute per-read mean methylation."""
    if methyl.empty or reads.empty:
        return None

    methyl_mean = methyl.mean(axis=1, skipna=True)
    methyl_cov = methyl.notna().sum(axis=1)

    df = reads.copy()
    df["mean_methyl"] = methyl_mean.values
    df["n_cpg_covered"] = methyl_cov.values

    # Filter: enough CpG coverage and valid methylation
    df = df[(df["n_cpg_covered"] >= MIN_CPGS_PER_READ) & df["mean_methyl"].notna()]
    if len(df) < MIN_READS:
        return None

    # Methylation class: High / Low / Middle
    df["methyl_cls"] = "Middle"
    df.loc[df["mean_methyl"] >= METHYL_HIGH_THR, "methyl_cls"] = "High"
    df.loc[df["mean_methyl"] <= METHYL_LOW_THR, "methyl_cls"] = "Low"

    # Binary methylation for Fisher test
    df["methyl_bin"] = (df["mean_methyl"] >= 0.5).map({True: "High", False: "Low"})

    # Simplified HP group
    df["hp_group"] = df["hp"].apply(simplify_hp)

    return df


# ===== Association statistics =====

def fisher_association(df: pd.DataFrame, group_col: str, methyl_col: str = "methyl_bin"):
    """Compute Fisher's exact test + Cramér's V for group vs methylation association."""
    valid_groups = {"HP1", "HP2"} if "HP" in group_col else {"ALT", "REF"}
    sub = df[df[group_col].isin(valid_groups)].copy()
    if len(sub) < 10:
        return None

    ct = pd.crosstab(sub[group_col], sub[methyl_col])
    if ct.shape != (2, 2):
        return None

    try:
        odds, p = stats.fisher_exact(ct.values)
        chi2 = stats.chi2_contingency(ct.values)[0]
        n = ct.values.sum()
        v = np.sqrt(chi2 / (n * (min(ct.shape) - 1))) if n > 0 else 0.0
        # Fraction of High-methyl reads per group
        frac = sub.groupby(group_col)["methyl_bin"].apply(lambda x: (x == "High").mean())
        return dict(
            n_reads=len(sub),
            p_value=p,
            odds_ratio=odds,
            cramers_v=v,
            contingency=ct,
            frac_high=frac.to_dict(),
        )
    except Exception:
        return None


# ===== Per-locus analysis =====

def analyze_locus(row: pd.Series):
    """Return a result dict for one locus row."""
    region_dir = find_region_dir(row["Chr"], int(row["Pos"]))
    if region_dir is None:
        return None

    reads, methyl = load_region(region_dir)
    if reads is None:
        return None

    df = build_read_features(reads, methyl)
    if df is None:
        return None

    hp_res = fisher_association(df, "hp_group")
    al_res = fisher_association(df, "alt_support")

    n_hp1 = (df["hp_group"] == "HP1").sum()
    n_hp2 = (df["hp_group"] == "HP2").sum()
    n_alt = (df["alt_support"] == "ALT").sum()
    n_ref = (df["alt_support"] == "REF").sum()

    rec = {
        "region_id": row["RegionID"],
        "chr": row["Chr"],
        "pos": row["Pos"],
        "selection_field": row["_SelectionField"],
        "selection_value": row[LEGACY_SELECTION_COLUMN],
        "dominant_label": row.get("DominantLabel", ""),
        "hp_merged_delta": row.get("HPMergedDelta", np.nan),
        "allele_delta": row.get("AlleleDelta", np.nan),
        "n_reads_total": len(reads),
        "n_reads_analyzed": len(df),
        "n_hp1": n_hp1,
        "n_hp2": n_hp2,
        "n_alt": n_alt,
        "n_ref": n_ref,
        "pct_high_methyl": (df["methyl_bin"] == "High").mean(),
        "mean_methyl_hp1": df.loc[df["hp_group"] == "HP1", "mean_methyl"].mean(),
        "mean_methyl_hp2": df.loc[df["hp_group"] == "HP2", "mean_methyl"].mean(),
        "mean_methyl_alt": df.loc[df["alt_support"] == "ALT", "mean_methyl"].mean(),
        "mean_methyl_ref": df.loc[df["alt_support"] == "REF", "mean_methyl"].mean(),
    }

    if hp_res:
        rec.update({
            "hp_n_phased": hp_res["n_reads"],
            "hp_p": hp_res["p_value"],
            "hp_odds": hp_res["odds_ratio"],
            "hp_cramers_v": hp_res["cramers_v"],
            "hp1_frac_high": hp_res["frac_high"].get("HP1", np.nan),
            "hp2_frac_high": hp_res["frac_high"].get("HP2", np.nan),
        })
    else:
        rec.update({"hp_n_phased": 0, "hp_p": np.nan, "hp_odds": np.nan,
                    "hp_cramers_v": np.nan, "hp1_frac_high": np.nan, "hp2_frac_high": np.nan})

    if al_res:
        rec.update({
            "allele_n": al_res["n_reads"],
            "allele_p": al_res["p_value"],
            "allele_odds": al_res["odds_ratio"],
            "allele_cramers_v": al_res["cramers_v"],
            "alt_frac_high": al_res["frac_high"].get("ALT", np.nan),
            "ref_frac_high": al_res["frac_high"].get("REF", np.nan),
        })
    else:
        rec.update({"allele_n": 0, "allele_p": np.nan, "allele_odds": np.nan,
                    "allele_cramers_v": np.nan, "alt_frac_high": np.nan, "ref_frac_high": np.nan})

    return rec


# ===== Main analysis =====

def run_analysis(sig_df: pd.DataFrame) -> pd.DataFrame:
    all_results = []
    for cls in TARGET_CLASSES:
        subset = sig_df[sig_df[LEGACY_SELECTION_COLUMN] == cls]
        if len(subset) > SAMPLE_N_PER_CLASS:
            subset = subset.sample(SAMPLE_N_PER_CLASS, random_state=RANDOM_SEED)
        print(f"\n[{cls}] Analyzing {len(subset)} loci...")
        ok, skip = 0, 0
        for i, (_, row) in enumerate(subset.iterrows()):
            if i % 50 == 0:
                print(f"  {i}/{len(subset)}...", end="\r", flush=True)
            res = analyze_locus(row)
            if res:
                all_results.append(res)
                ok += 1
            else:
                skip += 1
        print(f"  [{cls}] OK={ok}  Skipped={skip}     ")
    return pd.DataFrame(all_results)


# ===== Summary statistics =====

def summary_stats(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cls in TARGET_CLASSES:
        sub = results[results["selection_value"] == cls]
        if sub.empty:
            continue

        def pct_sig(col, alpha=0.05):
            valid = sub[col].dropna()
            return (valid < alpha).mean() * 100 if len(valid) > 0 else np.nan

        def median_q(col):
            return sub[col].median()

        rows.append({
            "class": cls,
            "n_loci": len(sub),
            # HP association
            "hp_sig_pct": pct_sig("hp_p"),
            "hp_cramers_v_median": median_q("hp_cramers_v"),
            "hp_methyl_delta_median": (sub["hp1_frac_high"] - sub["hp2_frac_high"]).abs().median(),
            # Allele association
            "allele_sig_pct": pct_sig("allele_p"),
            "allele_cramers_v_median": median_q("allele_cramers_v"),
            "allele_methyl_delta_median": (sub["alt_frac_high"] - sub["ref_frac_high"]).abs().median(),
            # Overall methylation
            "mean_methyl_median": median_q("pct_high_methyl"),
        })
    return pd.DataFrame(rows)


# ===== Plots =====

def make_plots(results: pd.DataFrame):
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Methylation Cluster ↔ HP/Allele Association\nHCC1395 5kHz (Legacy Strong vs Noise)",
                 fontsize=13, fontweight="bold")

    palette = {"Strong": "#e74c3c", "Noise": "#3498db"}
    classes = ["Strong", "Noise"]

    # 1. HP Cramér's V distribution
    ax = axes[0, 0]
    for cls in classes:
        sub = results[(results["selection_value"] == cls) & results["hp_cramers_v"].notna()]
        ax.hist(sub["hp_cramers_v"], bins=30, alpha=0.6, label=cls,
                color=palette[cls], density=True)
    ax.set_xlabel("HP Cramér's V")
    ax.set_ylabel("Density")
    ax.set_title("HP1/HP2 ↔ Methylation Association Strength")
    ax.legend()

    # 2. Allele Cramér's V distribution
    ax = axes[0, 1]
    for cls in classes:
        sub = results[(results["selection_value"] == cls) & results["allele_cramers_v"].notna()]
        ax.hist(sub["allele_cramers_v"], bins=30, alpha=0.6, label=cls,
                color=palette[cls], density=True)
    ax.set_xlabel("Allele Cramér's V")
    ax.set_title("ALT/REF ↔ Methylation Association Strength")
    ax.legend()

    # 3. HP Fisher -log10(p) violin
    ax = axes[0, 2]
    plot_data = []
    for cls in classes:
        sub = results[(results["selection_value"] == cls) & results["hp_p"].notna()]
        log_p = -np.log10(sub["hp_p"].clip(lower=1e-15))
        for v in log_p:
            plot_data.append({"class": cls, "-log10(p)": v})
    pdata = pd.DataFrame(plot_data)
    if not pdata.empty:
        sns.violinplot(data=pdata, x="class", y="-log10(p)", palette=palette,
                       order=classes, ax=ax, cut=0)
    ax.axhline(-np.log10(0.05), color="gray", linestyle="--", label="p=0.05")
    ax.set_title("HP Fisher Test -log10(p)")
    ax.legend(fontsize=8)

    # 4. HP methylation delta (HP1 vs HP2 fraction high)
    ax = axes[1, 0]
    for cls in classes:
        sub = results[results["selection_value"] == cls].copy()
        delta = (sub["hp1_frac_high"] - sub["hp2_frac_high"]).abs().dropna()
        ax.hist(delta, bins=30, alpha=0.6, label=cls, color=palette[cls], density=True)
    ax.set_xlabel("|HP1_HighFrac - HP2_HighFrac|")
    ax.set_title("HP Methylation Fraction Delta")
    ax.legend()

    # 5. Allele methylation delta (ALT vs REF fraction high)
    ax = axes[1, 1]
    for cls in classes:
        sub = results[results["selection_value"] == cls].copy()
        delta = (sub["alt_frac_high"] - sub["ref_frac_high"]).abs().dropna()
        ax.hist(delta, bins=30, alpha=0.6, label=cls, color=palette[cls], density=True)
    ax.set_xlabel("|ALT_HighFrac - REF_HighFrac|")
    ax.set_title("Allele Methylation Fraction Delta")
    ax.legend()

    # 6. Scatter: HP Cramér's V vs Allele Cramér's V
    ax = axes[1, 2]
    for cls in classes:
        sub = results[results["selection_value"] == cls].dropna(subset=["hp_cramers_v", "allele_cramers_v"])
        ax.scatter(sub["hp_cramers_v"], sub["allele_cramers_v"],
                   alpha=0.3, s=15, label=cls, color=palette[cls])
    ax.set_xlabel("HP Cramér's V")
    ax.set_ylabel("Allele Cramér's V")
    ax.set_title("HP vs Allele Association (per locus)")
    ax.legend()

    plt.tight_layout()
    out_path = OUTPUT_DIR / "methylation_cluster_association_summary.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Plot saved: {out_path}")


# ===== Observation document =====

def write_observation_doc(results: pd.DataFrame, summary: pd.DataFrame, metadata):
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    n_strong = len(results[results["selection_value"] == "Strong"])
    n_noise = len(results[results["selection_value"] == "Noise"])

    strong = summary[summary["class"] == "Strong"].iloc[0] if n_strong > 0 else {}
    noise = summary[summary["class"] == "Noise"].iloc[0] if n_noise > 0 else {}

    doc = f"""<!--
建立時間: {now}
目標: 觀察 legacy Strong/Noise loci 中甲基 cluster 與 HP/Allele 標籤的共存關係
處理範圍: HCC1395 5kHz paired_full (20260314)
關聯檔案:
  - scripts/analysis/analyze_methyl_cluster_allele_cooccurrence.py
  - {ISM_TP_DIR}
-->

# 觀察文件 1：甲基 Cluster ↔ HP/Allele 共存分析

**樣本**：HCC1395 5kHz paired_full
**分析時間**：{now}
**歷史 cohort 欄位**：`{metadata['selection_field']}`（schema status=`{metadata['schema_status']}`）
**歷史 cohort 值**：`Strong`、`Noise`
**未知值計數**：`{metadata['unknown_counts']}`（不納入 cohort）
**分析 loci**：legacy Strong={n_strong}，legacy Noise={n_noise}（各最多 {SAMPLE_N_PER_CLASS} 個，random seed={RANDOM_SEED}）

---

## 方法

1. 讀取 `significance_summary.csv`，按明確的 legacy four-state 欄位分組抽樣；不折疊 current v2 classes
2. 每個 locus 合併 `reads.tsv`（hp, alt_support）與 `methylation.csv`（per-CpG 機率）
3. 每個 read 計算 mean methylation（忽略 NA，需 ≥{MIN_CPGS_PER_READ} 個 CpG）
4. Binary 分類：mean ≥ 0.5 → High，< 0.5 → Low
5. HP1/HP2 × High/Low → Fisher's exact test + Cramér's V
6. ALT/REF × High/Low → Fisher's exact test + Cramér's V

---

## 結果總覽

| 指標 | Legacy Strong ({n_strong} loci) | Legacy Noise ({n_noise} loci) |
|------|------|------|
| **HP 顯著比例**（p<0.05） | {strong.get('hp_sig_pct', float('nan')):.1f}% | {noise.get('hp_sig_pct', float('nan')):.1f}% |
| **HP Cramér's V 中位數** | {strong.get('hp_cramers_v_median', float('nan')):.3f} | {noise.get('hp_cramers_v_median', float('nan')):.3f} |
| **HP 甲基 delta 中位數** | {strong.get('hp_methyl_delta_median', float('nan')):.3f} | {noise.get('hp_methyl_delta_median', float('nan')):.3f} |
| **Allele 顯著比例**（p<0.05） | {strong.get('allele_sig_pct', float('nan')):.1f}% | {noise.get('allele_sig_pct', float('nan')):.1f}% |
| **Allele Cramér's V 中位數** | {strong.get('allele_cramers_v_median', float('nan')):.3f} | {noise.get('allele_cramers_v_median', float('nan')):.3f} |
| **Allele 甲基 delta 中位數** | {strong.get('allele_methyl_delta_median', float('nan')):.3f} | {noise.get('allele_methyl_delta_median', float('nan')):.3f} |

---

## 詮釋框架

### HP 關聯（Haplotype-specific methylation）
- **高 Cramér's V**（>0.3）+ **p < 0.05** → HP1/HP2 有不同甲基模式
- 可能原因：(a) germline ASM（HP 決定的正常甲基差異），(b) tumor-specific clonal methylation change
- Legacy Strong > legacy Noise 代表歷史 four-state cohort 捕捉到較強的甲基-allele 共存；不得解讀為 current v2 class 比較

### Allele 關聯（Somatic allele-specific methylation）
- **高 Cramér's V** + **p < 0.05** → ALT/REF 等位基因有不同甲基模式
- 更直接反映體細胞 SNV 與甲基化狀態的共存
- 如果 Allele delta > HP delta，代表甲基變化比 haplotype 分離更緊密地與 somatic mutation 共存

---

## 後續行動建議

1. **如果 Strong Cramér's V 顯著高於 Noise**：
   - 進入方向 E（ML）：以 HP Cramér's V、Allele Cramér's V 作為訓練特徵
   - 現有 significance_summary.csv 已有 AlleleDelta、HPMergedDelta，可直接用

2. **如果 Allele 關聯強於 HP 關聯**：
   - 支持「somatic allele 與甲基共存」的核心論點
   - 可直接寫成第一個可主張的研究發現

3. **如果 Strong 與 Noise 差異不顯著**：
   - 代表 ISM 現有 Strong/Noise 分類不是甲基-allele 共存的有效代理
   - 需要重新定義「有意義的甲基共存」的判斷標準

---

## 圖表

- `methylation_cluster_association_summary.png`：6 panel 綜合圖
- `results_per_locus.csv`：每個 locus 的完整統計
- `summary_stats.csv`：legacy Strong vs Noise 匯總
- `analysis_metadata.json`：selection field/value、schema status、warning 與 unknown counts

"""
    doc_path = OUTPUT_DIR / "20260324_methylation_cluster_allele_cooccurrence_obs_01.md"
    doc_path.write_text(doc, encoding="utf-8")
    print(f"Observation doc: {doc_path}")


# ===== Entry point =====

def main():
    global ISM_TP_DIR, OUTPUT_DIR
    args = parse_args()
    ISM_TP_DIR = args.ism_tp_dir
    OUTPUT_DIR = args.output_dir
    significance_summary = args.significance_summary or (ISM_TP_DIR / "significance_summary.csv")

    print("=" * 60)
    print("Methylation Cluster <-> HP/Allele Co-occurrence Analysis")
    print("HCC1395 5kHz paired_full")
    print("=" * 60)

    sig_df = pd.read_csv(significance_summary)
    print(f"\nLoaded significance_summary: {len(sig_df)} loci")
    try:
        sig_df, metadata = select_replication_cohorts(
            sig_df,
            allow_unversioned_v1=args.allow_unversioned_v1,
        )
    except SchemaContractError as exc:
        print(f"ERROR: legacy verification schema contract failed: {exc}", file=sys.stderr)
        sys.exit(2)

    print(
        "Selection: "
        f"field={metadata['selection_field']} values={metadata['selection_values']} "
        f"schema_status={metadata['schema_status']} counts={metadata['selected_counts']}"
    )
    for message in metadata["warnings"]:
        print(f"WARNING: {message}", file=sys.stderr)
    if metadata["unknown_counts"]:
        print(
            f"WARNING: excluded unknown legacy values: {metadata['unknown_counts']}",
            file=sys.stderr,
        )

    # Do not create partial analysis artifacts until the cohort schema is valid.
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    results = run_analysis(sig_df)
    print(f"\nTotal analyzed loci: {len(results)}")

    # Save per-locus results
    results.to_csv(OUTPUT_DIR / "results_per_locus.csv", index=False)
    print(f"Per-locus results saved.")

    summary = summary_stats(results)
    summary.to_csv(OUTPUT_DIR / "summary_stats.csv", index=False)

    print("\n===== Summary =====")
    print(summary.to_string(index=False))

    make_plots(results)
    write_observation_doc(results, summary, metadata)

    metadata.update(
        {
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "significance_summary": str(significance_summary.resolve()),
            "ism_tp_dir": str(ISM_TP_DIR.resolve()),
            "output_dir": str(OUTPUT_DIR.resolve()),
            "sample_n_per_class": SAMPLE_N_PER_CLASS,
            "random_seed": RANDOM_SEED,
            "analyzed_counts": {
                name: int((results["selection_value"] == name).sum()) for name in TARGET_CLASSES
            },
        }
    )
    metadata_path = OUTPUT_DIR / "analysis_metadata.json"
    metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Analysis metadata: {metadata_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
