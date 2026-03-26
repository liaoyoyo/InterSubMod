#!/usr/bin/env python3
"""Build round-level markdown summary and plots from generated research outputs."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from research_common import markdown_table, read_json


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", required=True, help="Round sample bundle directory")
    return parser.parse_args()


def plot_benchmark(df: pd.DataFrame, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.barplot(data=df, x="method", y="f1", ax=ax)
    ax.set_ylim(0, min(1.0, max(df["f1"].astype(float).max() + 0.05, 0.2)))
    ax.set_title("Benchmark F1 Comparison")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_agreement_counts(df: pd.DataFrame, out_png: Path) -> None:
    order = df["agreement_type"].value_counts().index.tolist()
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.countplot(data=df, x="agreement_type", order=order, ax=ax)
    ax.set_title("Agreement Types")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_class_shift(df: pd.DataFrame, out_png: Path) -> None:
    shift_df = df[df["class_shift"] != "unchanged"].copy()
    fig, ax = plt.subplots(figsize=(8, 4))
    if shift_df.empty:
        ax.text(0.5, 0.5, "No class shifts", ha="center", va="center")
        ax.axis("off")
    else:
        order = shift_df["class_shift"].value_counts().index.tolist()
        sns.countplot(data=shift_df, x="class_shift", order=order, ax=ax)
        ax.tick_params(axis="x", rotation=30)
        ax.set_title("Class Shifts")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def write_failure_diagnosis(sample_dir: Path, context: Dict, metrics: Dict, agreement_df: pd.DataFrame, haplotag_df: pd.DataFrame) -> None:
    diagnosis_path = sample_dir / "failure_diagnosis.md"
    if diagnosis_path.exists():
        return

    filtered = metrics.get("filtered", {})
    baseline = metrics.get("baseline", {})
    lines: List[str] = []
    lines.append(f"# Failure Diagnosis: {context.get('sample', sample_dir.name)}")
    lines.append("")
    lines.append("## 現象")
    lines.append("")
    lines.append(f"- Baseline F1: `{baseline.get('f1', '')}`")
    lines.append(f"- Filtered F1: `{filtered.get('f1', '')}`")
    lines.append(f"- F1 delta: `{metrics.get('improvement', {}).get('f1_delta', '')}`")
    if not haplotag_df.empty:
        lines.append(f"- Sample HP assign rate: `{haplotag_df.iloc[0].get('hp_assign_rate', '')}`")
    lines.append("")
    lines.append("## 可能原因")
    lines.append("")
    if float(metrics.get("improvement", {}).get("f1_delta", 0.0)) <= 0:
        lines.append("- 目前規則沒有帶來正向 F1 增益，需檢查 filter 是否過嚴或特徵本身不足。")
    if not haplotag_df.empty and float(haplotag_df.iloc[0].get("hp_assign_rate", 0.0)) < 0.05:
        lines.append("- HP assign rate 偏低，可能使 label-first HP 驗證不具可解釋性。")
    conflict_count = int((agreement_df["agreement_type"] == "conflict").sum()) if not agreement_df.empty else 0
    if conflict_count > 0:
        lines.append(f"- 存在 `{conflict_count}` 個 label/cluster 衝突 region，需抽樣做 diagnostics。")
    if conflict_count == 0 and float(metrics.get("improvement", {}).get("f1_delta", 0.0)) > 0:
        lines.append("- 目前沒有明顯失敗訊號，可轉向擴大樣本或做敏感度矩陣。")
    lines.append("")
    lines.append("## 下一個可驗證修改點")
    lines.append("")
    lines.append("- 調整 read filter，確認是否因過嚴而降低 label-first 可用性。")
    lines.append("- 比較 BERNOULLI / NHD / L1 距離法是否改變 upgrade/conflict 分布。")
    lines.append("- 對 conflict 或 label_upgrade region 匯出 diagnostics 圖與表。")
    diagnosis_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_context(sample_dir: Path) -> Dict:
    for name in ("round_context.json", "run_context.json"):
        path = sample_dir / name
        if path.exists():
            return read_json(path)
    return {}


def main() -> None:
    args = parse_args()
    sample_dir = Path(args.sample_dir).resolve()
    plots_dir = sample_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    context = load_context(sample_dir)
    metrics = read_json(sample_dir / "metrics.json")

    benchmark_path = sample_dir / "benchmark_comparison.tsv"
    design_path = sample_dir / "method_design_validation.tsv"
    agreement_path = sample_dir / "label_cluster_agreement.tsv"
    haplotag_path = sample_dir / "haplotag_qc.tsv"

    benchmark_df = pd.read_csv(benchmark_path, sep="\t", low_memory=False) if benchmark_path.exists() else pd.DataFrame()
    design_df = pd.read_csv(design_path, sep="\t", low_memory=False) if design_path.exists() else pd.DataFrame()
    agreement_df = pd.read_csv(agreement_path, sep="\t", low_memory=False) if agreement_path.exists() else pd.DataFrame()
    haplotag_df = pd.read_csv(haplotag_path, sep="\t", low_memory=False) if haplotag_path.exists() else pd.DataFrame()

    plot_links: List[str] = []
    if not benchmark_df.empty:
        plot_benchmark(benchmark_df, plots_dir / "benchmark_f1.png")
        plot_links.append("- [benchmark_f1.png](plots/benchmark_f1.png)")
    if not agreement_df.empty:
        plot_agreement_counts(agreement_df, plots_dir / "agreement_types.png")
        plot_class_shift(agreement_df, plots_dir / "class_shifts.png")
        plot_links.append("- [agreement_types.png](plots/agreement_types.png)")
        plot_links.append("- [class_shifts.png](plots/class_shifts.png)")

    write_failure_diagnosis(sample_dir, context, metrics, agreement_df, haplotag_df)

    filtered = metrics.get("filtered", {})
    baseline = metrics.get("baseline", {})
    improvement = metrics.get("improvement", {})
    upgrade_count = int((agreement_df["agreement_type"] == "label_upgrade").sum()) if not agreement_df.empty else 0
    conflict_count = int((agreement_df["agreement_type"] == "conflict").sum()) if not agreement_df.empty else 0

    lines: List[str] = []
    lines.append(f"# Round Summary: {context.get('sample', sample_dir.name)}")
    lines.append("")
    lines.append("## 基本資訊")
    lines.append("")
    lines.append(f"- 時間: `{context.get('created_at', '')}`")
    lines.append(f"- 樣本: `{context.get('sample', '')}`")
    lines.append(f"- 樣本集合: `{context.get('sample_set', '')}`")
    lines.append(f"- 平台: `{context.get('platform', '')}`")
    lines.append(f"- 測試重點: `{context.get('test_focus', '')}`")
    lines.append(f"- 改動內容: `{context.get('changes', '')}`")
    lines.append(f"- 來源 run: `{context.get('source_run_dir', '')}`")
    lines.append("")
    lines.append("## 主要結果")
    lines.append("")
    lines.append(f"- Baseline F1: `{baseline.get('f1', '')}`")
    lines.append(f"- Filtered F1: `{filtered.get('f1', '')}`")
    lines.append(f"- F1 delta: `{improvement.get('f1_delta', '')}`")
    if not haplotag_df.empty:
        lines.append(f"- HP assign rate: `{haplotag_df.iloc[0].get('hp_assign_rate', '')}`")
    lines.append(f"- Label upgrades: `{upgrade_count}`")
    lines.append(f"- Label/cluster conflicts: `{conflict_count}`")
    lines.append("")
    lines.append("## Benchmark Table")
    lines.append("")
    if not benchmark_df.empty:
        lines.append(
            markdown_table(
                ["method", "tp", "fp", "fn", "precision", "recall", "f1", "delta_f1_vs_baseline"],
                benchmark_df.to_dict(orient="records"),
            )
        )
    else:
        lines.append("_No benchmark comparison found_")
    lines.append("")
    lines.append("## 輸出檔案")
    lines.append("")
    lines.append("- [benchmark_comparison.tsv](benchmark_comparison.tsv)")
    lines.append("- [benchmark_comparison.md](benchmark_comparison.md)")
    lines.append("- [method_design_validation.tsv](method_design_validation.tsv)")
    lines.append("- [label_cluster_agreement.tsv](label_cluster_agreement.tsv)")
    lines.append("- [failure_diagnosis.md](failure_diagnosis.md)")
    if plot_links:
        lines.append("")
        lines.append("## 圖片")
        lines.append("")
        lines.extend(plot_links)
    if not design_df.empty:
        lines.append("")
        lines.append("## 方法學觀察")
        lines.append("")
        lines.append(f"- 設計驗證列數: `{len(design_df.index)}`")
        lines.append(f"- agreement 列數: `{len(agreement_df.index)}`")
    lines.append("")
    lines.append("## 初步結論")
    lines.append("")
    lines.append(f"- `{context.get('conclusion', '待補結論')}`")
    lines.append(f"- 下一步: `{context.get('next_step', '挑選 conflict / upgrade region 做 diagnostics')}`")

    summary_path = sample_dir / "round_summary.md"
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[build_round_dashboard] Wrote {summary_path}")


if __name__ == "__main__":
    main()
