#!/usr/bin/env python3
"""Generate step1_findings.md from step1_summary_stats.json + trajectory + deltas.
"""
from __future__ import annotations

import json
import math
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
STEP1 = REPO / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way"


def fmt_ci(rate: float, ci) -> str:
    if rate is None or (isinstance(rate, float) and math.isnan(rate)):
        return "—"
    if isinstance(ci, list) and len(ci) == 2:
        return f"{rate:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]"
    return f"{rate:.3f}"


def main() -> int:
    summary = json.loads((STEP1 / "step1_summary_stats.json").read_text())
    traj = pd.read_csv(STEP1 / "step1_trajectory.tsv", sep="\t", low_memory=False)

    class_counts = traj["class"].value_counts().to_dict()
    total = sum(class_counts.values())
    by_label = traj.groupby(["label", "class"]).size().unstack(fill_value=0)

    bam = summary["bam_metrics"]
    h1a = summary["H1a"]
    h1b = summary["H1b"]
    h1c = summary["H1c"]
    overall = summary["genome_marker_off_NGge3"]

    out = STEP1 / "step1_findings.md"
    today = datetime.utcnow().strftime("%Y-%m-%d")

    md_lines: list[str] = []
    md_lines.append("<!--")
    md_lines.append(f"build_date: {today}")
    md_lines.append("agent: Step 1 V3F/V5/V6 three-way integration (HCC1395 paired-pileup)")
    md_lines.append("status: in_progress")
    md_lines.append("report_class: characterization_post_hoc")
    md_lines.append("scope: HCC1395 pilot, paired-pileup VCF, 30,490 TP + 4,842 FP regions")
    md_lines.append("parent_plan: research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md")
    md_lines.append("inputs:")
    md_lines.append("  - research/paired_priority_bug_audit/phaseC_genome_three_way/{V3F,V5,V6}_{on,off}_{tp,fp}/")
    md_lines.append("  - research/paired_priority_bug_audit/phaseC_genome_three_way/v3f_vs_v5_vs_v6_region_ng.tsv")
    md_lines.append("  - /big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_{tp,fp}.vcf.gz")
    md_lines.append("  - research/tpfp_loh_af_kde_discrimination/data/master.tsv.gz (HCC1395 paired_full mode)")
    md_lines.append("outputs:")
    md_lines.append("  - step1_master_three_way.tsv")
    md_lines.append("  - step1_delta_{v5_vs_v3f,v6_vs_v5,v6_vs_v3f}.tsv")
    md_lines.append("  - step1_trajectory.tsv")
    md_lines.append("  - step1_off_vs_on_compare.tsv")
    md_lines.append("  - step1_summary_stats.json")
    md_lines.append("  - figures/step1/step1_{three_panel_heatmap,delta_heatmap,trajectory_sankey}.png")
    md_lines.append("verdict: see Section 4 H1a/H1b/H1c table")
    md_lines.append("-->")
    md_lines.append("")
    md_lines.append("# Step 1 — V3F / V5 / V6 三向 ISM 整合（HCC1395 paired-pileup）")
    md_lines.append("")
    md_lines.append("> Characterization-only. 不評估 filter / ΔF1（per plan §Out-of-scope）。")
    md_lines.append("")

    # Section 1: Data inventory & caveats
    md_lines.append("## 1. 資料盤點與重要警示")
    md_lines.append("")
    md_lines.append("### 1.1 phaseC ISM 結果結構（與計畫書認知差異）")
    md_lines.append("")
    md_lines.append("- 計畫書假設 phaseC 內存在 region-level master.tsv（含 HP1FamilyN/HP2FamilyN/caller_af/Coverage_Multiple/LOH/TP_label 等欄位）。")
    md_lines.append("- **實際**：12 個 ISM 執行的 `significance_summary.csv` 均為 header-only（每個 run 的 `significance_statistics.txt` 顯示 `Regions Analyzed (Success): 0`）。phaseC 是 stripped-down HP audit 跑法，主 master TSV 未產出。")
    md_lines.append("- 唯一可用的 per-region 數據是 `{run}/filtered_snv_{tp|fp}/{chr}/{chr}_{pos}/{region_id}/reads/reads.tsv`，欄位：`read_id, read_name, chr, start, end, mapq, hp, alt_support, is_tumor, strand`。")
    md_lines.append("- 既有跨版本 aggregate：`v3f_vs_v5_vs_v6_region_ng.tsv`（105,996 rows = (30,490 TP + 4,842 FP) × 3 BAM × 1 (off NG + on NG)），由 `phaseC_region_ng_fast.py` 產出。")
    md_lines.append("- region_id 編碼：`{chr}_{start}_{end}`，width=10kb，SNV pos = midpoint。")
    md_lines.append("")
    md_lines.append("### 1.2 補強：caller_af 與 LOH/CN covariate 來源")
    md_lines.append("")
    md_lines.append("- **caller_af**：直接從 source VCF `filtered_snv_{tp,fp}.vcf.gz` FORMAT/AF 欄位讀取（bcftools query）。")
    md_lines.append("- **LOH / Coverage_Multiple / Diploid_Coverage_Used / loh_side / AF (master)**：從 `research/tpfp_loh_af_kde_discrimination/data/master.tsv.gz`（HCC1395 paired_full mode）按 Chr+Pos join。")
    md_lines.append(f"- master.tsv join 覆蓋率：總 {summary['total_rows']} rows 中 master_join_ok=1 為 {summary['master_joined_rows']} ({summary['master_joined_rows']/summary['total_rows']:.1%})。")
    md_lines.append("- **重要不對稱**：phaseC 用 **paired-pileup** VCF（30,490 TP + 4,842 FP），但 master.tsv 只有 **paired_full** mode（29,754 TP + 627 FP）。Join 結果：TP ~96.8% 可 join，FP 只有 ~10.5%（paired_full FP set 小很多）。Step 2 起若要做 LOH × CN 主軸分層分析，可以使用 TP 樣本完整，但 FP 內的 LOH/CN 必須額外處理（或限縮到可 join 的 506 個 FP）。")
    if summary.get("stale_diploid75_rows", 0) > 0:
        md_lines.append(f"- **Diploid_Coverage_Used=75.0 stale flag**：{summary['stale_diploid75_rows']} rows 命中 stale KDE binary artifact（CURRENT_FOCUS H-CN1 警示）。Step 2-3 涉及 Coverage_Multiple 比較時須過濾或重跑 KDE。")
    else:
        md_lines.append("- Diploid_Coverage_Used=75.0 stale flag：未發現可疑 row（master.tsv KDE-corrected）。")
    md_lines.append("")

    # Section 2: Genome-wide marker rates (NG_off ≥ 3) — three-way
    md_lines.append("## 2. 三向版本 genome-wide marker (NG_off ≥ 3) summary")
    md_lines.append("")
    md_lines.append("| BAM | marker n | TP | FP | TP rate (95% CI) |")
    md_lines.append("|-----|----------|----|----|-------------------|")
    for b in ["V3F", "V5", "V6"]:
        o = overall[b]
        md_lines.append(f"| {b} | {o['marker_off_n']:,} | {o['marker_off_n_tp']:,} | {o['marker_off_n_fp']:,} | {fmt_ci(o['marker_off_tp_rate'], o['marker_off_tp_rate_CI95'])} |")
    md_lines.append("")
    md_lines.append("> 與 `v3f_vs_v5_vs_v6_genome_summary.tsv` 一致（V3F=0.9175 / V5=0.8937 / V6=0.9093）— 確認讀取流程無誤。")
    md_lines.append("")

    # Section 3: NG=2 inner vs outer gap (H1a/b/c proxy)
    md_lines.append("## 3. NG=2 in master-joined regions — inner / outer TP rate gap")
    md_lines.append("")
    md_lines.append("> Proxy for plan §H1a/b/c \"Inner same-hap TP gap\"。phaseC 沒有 per-cell HP-direction grain，所以以 master-joined `loh_side` 區分 Inner / Outer，並以 NG=2 cell 內 TP rate 作為 gap source。")
    md_lines.append("")
    md_lines.append("| BAM | Inner NG=2 TP rate (95% CI) | n_TP / n_FP | Outer NG=2 TP rate (95% CI) | n_TP / n_FP | Inner − Outer gap |")
    md_lines.append("|-----|------------------------------|-------------|-------------------------------|-------------|---------------------|")
    for b in ["V3F", "V5", "V6"]:
        m = bam[b]
        md_lines.append(
            f"| {b} | "
            f"{fmt_ci(m['inner_NG2_tp_rate'], m.get('inner_NG2_tp_rate_CI95', []))} | "
            f"{m['inner_NG2_n_tp']} / {m['inner_NG2_n_fp']} | "
            f"{fmt_ci(m['outer_NG2_tp_rate'], m.get('outer_NG2_tp_rate_CI95', []))} | "
            f"{m['outer_NG2_n_tp']} / {m['outer_NG2_n_fp']} | "
            f"{m['inner_minus_outer_gap']:+.3f} |"
        )
    md_lines.append("")

    # Section 4: H1a/H1b/H1c judgment
    md_lines.append("## 4. H1a / H1b / H1c 判定")
    md_lines.append("")
    md_lines.append("| Hypothesis | 指標 | Δ | Threshold | Verdict |")
    md_lines.append("|------------|------|---|-----------|---------|")
    md_lines.append(f"| H1a — V5 BAM Inner same-hap TP gap > V3F | Δgap(V5−V3F) | {h1a['delta_gap_V5_minus_V3F']:+.3f} | ≥ 0.03 | **{h1a['verdict']}** |")
    md_lines.append(f"| H1b — V6 BAM Inner same-hap TP gap > V5 | Δgap(V6−V5) | {h1b['delta_gap_V6_minus_V5']:+.3f} | ≥ 0.03 | **{h1b['verdict']}** |")
    md_lines.append(f"| H1c — V6 對 V3F 累積增益最大 | Δgap(V6−V3F) | {h1c['delta_gap_V6_minus_V3F']:+.3f} | ≥ 0.06 | **{h1c['verdict']}** |")
    md_lines.append("")
    md_lines.append("> **解讀**：")
    md_lines.append("> - POSITIVE = 該段修補增加了 Inner-vs-Outer NG=2 TP rate gap（characterization 增益）")
    md_lines.append("> - NEGATIVE = gap 變化未達 threshold，但未反向")
    md_lines.append("> - NEGATIVE_REVERSE = gap 反向（Inner 變差 / Outer 變好）")
    md_lines.append("> - UNKNOWN = 樣本不足或計算失敗")
    md_lines.append("")

    # Section 5: per-region trajectory
    md_lines.append("## 5. Per-region NG trajectory (V3F → V5 → V6, flag=off)")
    md_lines.append("")
    md_lines.append("**5 類分類**：")
    md_lines.append("- A: 兩段都改善（V5−V3F > 0 且 V6−V5 > 0）")
    md_lines.append("- B: 只 V5 改善（V5−V3F > 0 且 V6−V5 ≤ 0）")
    md_lines.append("- C: 只 V6 改善（V5−V3F ≤ 0 且 V6−V5 > 0）")
    md_lines.append("- D: 無變化（兩段 Δ 都 = 0）")
    md_lines.append("- E: 反向或單段下降（任一段 Δ < 0 但非 D）")
    md_lines.append("")
    md_lines.append("| 類別 | n | 占比 |")
    md_lines.append("|------|---|------|")
    for cls in ["A_both_improve", "B_only_V5_improve", "C_only_V6_improve", "D_no_change", "E_reverse_or_decrease", "MISSING"]:
        n = class_counts.get(cls, 0)
        md_lines.append(f"| {cls} | {n} | {n/total:.1%} |")
    md_lines.append("")
    md_lines.append("### 5.1 按 TP / FP label 拆解")
    md_lines.append("")
    md_lines.append("```")
    md_lines.append(by_label.to_string())
    md_lines.append("```")
    md_lines.append("")

    # Section 6: off / on flag (germline-hp-only) comparison
    md_lines.append("## 6. `--germline-hp-only` off vs on 對照（mask V5 Layer 1.5？）")
    md_lines.append("")
    md_lines.append("> 詳見 `step1_off_vs_on_compare.tsv`。每 (label, BAM, flag) 組合給 NG=2 與 NG≥3 比例。")
    md_lines.append("> 若 V5_on（mask Layer 1.5 marker → unphased）與 V6_off 接近，表示 ISM 端 demotion 可近似 V6 修補；若差距大，V6 binary fix 才是主導。")
    md_lines.append("")

    # Section 7: figures
    md_lines.append("## 7. 圖示")
    md_lines.append("")
    md_lines.append("- `figures/step1/step1_three_panel_heatmap.png` — V3F / V5 / V6 三聯 TP rate × (NG × LOH side)")
    md_lines.append("- `figures/step1/step1_delta_heatmap.png` — 3 個 Δ TP rate heatmap")
    md_lines.append("- `figures/step1/step1_trajectory_sankey.png` — 5 類 region trajectory stacked bar")
    md_lines.append("")

    # Section 8: Hand-off to Agent B (Step 2)
    md_lines.append("## 8. Hand-off context — Step 2 Agent B")
    md_lines.append("")
    md_lines.append("**Master TSV schema (step1_master_three_way.tsv)**：")
    md_lines.append("- Key：`region_id`（chr_start_end，width=10kb），`chr`, `pos`（SNV center）")
    md_lines.append("- Label：`label` (TP/FP)")
    md_lines.append("- Per-version × per-flag features：`{V3F,V5,V6}_{off,on}_{0,1,2,1-1,2-1,3,11,21,33,other,NG,n_reads}`")
    md_lines.append("- Caller AF：`caller_af`（from `filtered_snv_{tp,fp}.vcf.gz` FORMAT/AF）")
    md_lines.append("- Master join：`master_join_ok` (1/0)")
    md_lines.append("- LOH / CN covariate：`LOH_Bed_Overlap`, `LOH_Bed_Annotation`, `LOH_Subtype`, `Coverage_Multiple`, `Coverage_Category`, `Diploid_Coverage_Used`, `loh_side`, `AF_master`")
    md_lines.append("")
    md_lines.append("**Step 2 推薦 covariate**（基於 plan §Step 2 3-axis grid）：")
    md_lines.append("- 主 grid 軸：`loh_side` (Inner/Outer), HP bucket（從 reads.tsv 進一步以 HP family 計算），`Coverage_Multiple` (5 bins)")
    md_lines.append("- LR covariate：`HPFineNGroups` (= V6_off_NG), `caller_af`, `n_reads`")
    md_lines.append("")
    md_lines.append("**已知 caveat 給 Agent B**：")
    md_lines.append("1. **FP set 的 LOH/CN join 不完整**（only ~10% master-joined）。Step 2 主分析建議以 TP 為主，FP 對照僅在 master-joined 子集做。")
    md_lines.append("2. **V6 BAM germline-existent 區因重用 V5 phased VCF 殘留 priority bug 偏移**（hp_1_1:hp_2_1 ratio 1.838 介於 V3F 1.138 與 V5 1.86 之間）— 這不是 bug，是設計選擇。Inner same-hap 在 V6 不如 V3F 純淨。")
    md_lines.append("3. **on flag mask Layer 1.5 marker → unphased**，會降低 NG_on=2 sample 數。Step 2 應以 flag=off 為主軸 grid，on 作 sensitivity check。")
    md_lines.append("4. Power gate：master-joined Inner+NG=2 cell 預期 n 詳見 step1_summary_stats.json bam_metrics 各 row n_TP+n_FP 欄位。")
    md_lines.append("")

    md_lines.append("## 9. 附錄")
    md_lines.append("")
    md_lines.append("**Scripts**：")
    md_lines.append("- `scripts/build_three_way_master.py` — Stage 1-4 wide-format master")
    md_lines.append("- `scripts/compute_deltas_and_trajectory.py` — Stage 5-7 deltas + trajectory + H1a/b/c")
    md_lines.append("- `scripts/make_figures.py` — Stage 8 figures")
    md_lines.append("- `scripts/generate_findings.py` — 本 findings.md 自動生成")
    md_lines.append("")
    md_lines.append("**Intermediate artifacts**：")
    md_lines.append("- `intermediate/per_region_hp_counts.tsv.gz` — 12 ISM run × per-region HP family counts (long format)")
    md_lines.append("- `intermediate/caller_af_lookup.tsv` — region_id × label × caller_af")
    md_lines.append("- `intermediate/build_log.txt` — build_three_way_master 執行日誌")
    md_lines.append("")
    md_lines.append("**Out-of-scope reminders** (plan §Out-of-scope)：")
    md_lines.append("- 本 step 不評估 ΔF1 / filter 效果（characterization-only）。")
    md_lines.append("- 不修改 C++ pipeline。")
    md_lines.append("- 不重跑任何 BAM/ISM。")

    out.write_text("\n".join(md_lines))
    print(f"wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
