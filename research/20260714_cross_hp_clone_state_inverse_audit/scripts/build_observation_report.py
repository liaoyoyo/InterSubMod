#!/usr/bin/env python3
"""Build an auditable Markdown record and a standalone professor-facing HTML report."""

from __future__ import annotations

import argparse
import html
import json
from datetime import datetime
from pathlib import Path
from typing import Any


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def fmt(value: int | float) -> str:
    return f"{value:,.3f}" if isinstance(value, float) else f"{value:,}"


def pct(numerator: int, denominator: int, digits: int = 2) -> str:
    return "—" if denominator == 0 else f"{100 * numerator / denominator:.{digits}f}%"


def esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


def find_record(records: list[dict[str, Any]], dataset: str, region: str) -> dict[str, Any]:
    for record in records:
        if record["dataset"] == dataset and record["region"] == region:
            return record
    raise KeyError(f"Missing candidate {dataset} {region}")


def markdown_sample_table(summaries: list[dict[str, Any]]) -> str:
    lines = [
        "| Dataset | Biological ID | Regions | 雙 primary HP | Topo enumerated / complete | 單樹 stored / complete | Observed strict | CN contract |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for summary in summaries:
        counts = summary["counts"]
        cn = summary["copy_number_contract"]
        cn_text = f"{cn['source']}（coarse）" if cn["availability"] == "measured" else "unavailable"
        lines.append(
            f"| {summary['dataset']} | {summary['biological_id']} | "
            f"{counts['regions_total']:,} | {counts['regions_two_primary_hp']:,} | "
            f"{counts['regions_direct_sister_shape_invariant']:,} / "
            f"{counts['regions_direct_sister_shape_invariant_pair_complete']:,} | "
            f"{counts['regions_direct_sister_tree_unique']:,} / "
            f"{counts['regions_direct_sister_tree_unique_pair_complete']:,} | "
            f"{counts['regions_observed_hp1_direct_only_hp2_sister_only']:,} | {cn_text} |"
        )
    return "\n".join(lines)


def make_markdown(
    audit: dict[str, Any], comparison: dict[str, Any], ascn: dict[str, Any]
) -> str:
    counts = audit["aggregate_counts"]
    sets = comparison["exact_region_set_comparisons"]
    stable = ascn["category_stable_pass_counts"]
    cn_mother = ascn["dual_hp_mother_set_screen"]
    created = datetime.now().astimezone().isoformat(timespec="seconds")
    return f"""<!--
建立時間: {created}
目標: 驗證現有 7-dataset 資料是否能觀察跨 HP direct/sister、同位點 collision，並條件式反推 clone mixture
處理範圍: chr1-22；LongPhase-S PASS canonical 7 datasets / 6 biological samples；ClairS PASS sensitivity
關聯檔案: InterSubMod/research/20260714_cross_hp_clone_state_inverse_audit/results/cross_hp_candidate_audit.json；InterSubMod/research/20260714_cross_hp_clone_state_inverse_audit/results/candidate_regions.tsv
-->

# 跨 HP clone-state 反推：現有資料可行性、候選觀察與驗證紀錄

用 SCQA + Claim–Evidence–Limit：先回答能不能，再拆出資料證據、反推條件與目前不能越過的界線。

> **TL;DR — 可以做候選結構 census 與定點 raw-read 驗證；目前不能從 bulk ONT 的兩個 HP marginals 直接確認四個 biological clones。** 目前枚舉結果有 {counts['regions_direct_sister_shape_invariant']:,} 個 direct/sister shape-invariant，其中 {counts['regions_direct_sister_shape_invariant_pair_complete']:,} 個兩側 candidate set 完整；{counts['regions_direct_sister_tree_unique']:,} 個目前各側只儲存一棵樹，其中 {counts['regions_direct_sister_tree_unique_pair_complete']:,} 個才是 analysis-complete tree-unique。固定兩-sSNV五成分 catalog 完整反推仍為 0（影響：高，信心：高）。

## 1. 直答判定

| 問題 | 判定 | 數據與理由 |
|---|---|---|
| 能否觀察每區域、每 HP 的 read states 與候選樹？ | **可以，但 canonical 有 threshold** | 7/7 canonical 完成；可重跑 MINREAD≥3 的 thresholded R/A/X states、counts 與 trees；低支持 raw states 需由 BAM bounded re-extraction 取得 |
| 能否找出一側 direct、另一側 sister？ | **可以** | 目前枚舉 shape-invariant {counts['regions_direct_sister_shape_invariant']:,}（analysis-complete {counts['regions_direct_sister_shape_invariant_pair_complete']:,}）；目前單樹 {counts['regions_direct_sister_tree_unique']:,}（analysis-complete {counts['regions_direct_sister_tree_unique_pair_complete']:,}） |
| 能否觀察同一 sSNV 在 HP1/HP2 都有 ALT？ | **可以初篩** | 至少一位點 {counts['regions_two_primary_hp_any_cross_hp_same_site_alt']:,}；全部位點 {counts['regions_two_primary_hp_all_sites_cross_hp_same_site_alt']:,} |
| 固定五成分模型在數學上可否反推？ | **有條件可以** | 固定 catalog 時 design rank=5；需 constrained likelihood |
| 現有真實區域是否達完整反推門檻？ | **沒有** | k=2 + unique direct/sister + all-site collision + full catalog + nonnegative = 0 |
| 能否宣稱 biological clone 或 joint clone tree？ | **不可以** | bulk HP marginals 無 cell pairing；沒有正交 clone truth |

## 2. 資料範圍與單位

- 7 datasets = 6 biological samples；HCC1395 與 HCC1395_DORADO 是同一細胞系的兩個技術資料版本。
- Canonical：20260713_layered_reconstruction_v3_raw_all_lps_pass_v5，LongPhase-S recalibrated PASS，chr1–22。
- 638,259 raw/all records → 582,820 LongPhase-S PASS → 469,849 chr1–22 biallelic sSNVs → 194,149 retained sSNVs → 51,815 regions。
- 跨 HP 比較母集合：22,779 個同時有 mutation-bearing HP1 與 HP2 primary lineage 的 regions。
- 玩具例 G1 是 germline anchor，不進 somatic state vector；G1-S1-S2 在 canonical 的 k 是 2，不是 3。

## 3. 平行證據鏈

### Solver topology

- 22,779 雙 primary HP regions。
- 18,588 兩側 candidate sets complete，占雙 HP {pct(counts['regions_pair_complete'], counts['regions_two_primary_hp'])}。
- 目前枚舉／儲存的結果有 {counts['regions_direct_sister_shape_invariant']:,} 個維持 direct-only / sister-only；其中 {counts['regions_direct_sister_shape_invariant_pair_complete']:,} 個兩側 candidate set 完整，才可作 exhaustive shape-invariant 解讀。
- 目前各側只儲存一棵樹為 {counts['regions_direct_sister_tree_unique']:,} 個；其中 {counts['regions_direct_sister_tree_unique_pair_complete']:,} 個同時 candidate-set complete，才可稱 analysis-complete tree-unique。

### Cross-HP same-site ALT

- 2,287 個至少有一個 sSNV 在 HP1/HP2 都有 ALT≥3。
- 86 個區域的所有 sSNV 都有雙 HP ALT≥3；其中 k=2 為 60。
- 與 direct+sister 取交集：{counts['regions_direct_sister_shape_invariant_pair_complete']:,} 個 complete topology 中 {counts['regions_direct_sister_shape_invariant_pair_complete_any_collision']:,} 個有 any collision；{counts['regions_direct_sister_tree_unique_pair_complete']:,} 個 complete unique 中 {counts['regions_direct_sister_tree_unique_pair_complete_any_collision']:,} 個；all-site collision 為 0。

54 與 86 是兩個平行集合，**不可寫成 86 → 54 的線性漏斗**。

## 4. 兩種 direct/sister 定義

| 定義 | 單位 | 數量 | 能說什麼 |
|---|---|---:|---|
| Solver shape-invariant（目前枚舉） | 已生成／儲存候選樹 | {counts['regions_direct_sister_shape_invariant']:,} | discovery metric；含 capped，不保證 exhaustive |
| Solver shape-invariant（analysis-complete） | 兩側完整 candidate sets | {counts['regions_direct_sister_shape_invariant_pair_complete']:,} | 所有可行候選均維持粗拓撲類別 |
| Solver 單樹（目前儲存） | 每 HP output | {counts['regions_direct_sister_tree_unique']:,} | 含 capped；n_trees=1 不必然等於真正唯一 |
| Solver tree-unique（analysis-complete） | 兩側完整 candidate sets | {counts['regions_direct_sister_tree_unique_pair_complete']:,} | 兩側各只有一棵可行樹 |
| Observed full-state geometry broad | thresholded full reads | 14 | HP1 有 subset chain；HP2 有 incomparable pair |
| Observed strict HP1 direct-only + HP2 sister-only | thresholded full reads | 7 | observed 幾何符合方向；不代表 joint cell pairing |

7 個 observed strict regions 全為 single non-missing PS 且 pair-complete，但 canonical coarse CN-neutral proxy 為 0。此處 single-PS 只表示非缺值 PS 的 distinct count=1，不表示每條 alignment 都帶 PS，也不證明跨區域 phase continuity；例如 H2009 HP1 為 64/75、HP2 為 77/77，HP1 缺 PS 的 11 條全為 XXA，移除後不改變 major-PS sensitivity 判讀。

## 5. 各 dataset

{markdown_sample_table(audit['sample_summaries'])}

7/7 region-view ↔ MLHP group conservation 通過：missing join=0、duplicate key=0。

## 6. 固定五成分模型

~~~
              direct HP      sister HP
Background       RR             RR
C1                I             RR
C2               AA             RR
C3               RR              U
C4                I              V

pi2 = p_direct(AA)
pi3 = p_sister(U)
pi4 = p_sister(V)
pi0 = p_direct(RR) - p_sister(U)
pi1 = p_direct(I)  - p_sister(V)
~~~

合成 contract test 已精確回復 pi=(0.10, 0.20, 0.25, 0.15, 0.30)。其中 π0 是 RR|RR background component；沒有外部 purity 與 normal/tumor model 時，不能直接命名為 normal contamination。若允許完整 3×3 HP1×HP2 pairing，兩組 marginals 仍有 4 個不可識別自由度，所以 catalog 必須被明示且與替代模型比較。

## 7. 為何目前是 0，而不是沒有 subclone

1. {counts['regions_direct_sister_shape_invariant_pair_complete']:,} 個 analysis-complete solver direct+sister 候選的 k 從 3 起，沒有 k=2。
2. 60 個 k=2 all-site collision regions 多為 direct+direct（38）或 sister+sister（7）。
3. 現成 populations_by_hp 與 subread_groups_by_hp 已套 MINREAD=3，低頻 states 被摘要刪除。
4. Full reads 也存在於 subread_groups；兩者相加會重複計數。
5. HP1/HP2 是 bulk marginals，沒有 cell barcode。

所以 0 的正確意思是：**thresholded canonical summary 中，固定且很窄的五成分 catalog 沒有 inference-ready region**；不是證明生物上不存在此類 clone。

## 8. CN、PS 與 bounded raw-read 驗證

- External SAVANA 宣告 gate 在 22,779 個 dataset-regions 中篩出 {cn_mother['magnitude_pass']:,}（{pct(cn_mother['magnitude_pass'], cn_mother['total'])}）1+1 magnitude candidates；加入目前 solution-stability heuristic 後為 {cn_mother['stable_pass']:,}（{pct(cn_mother['stable_pass'], cn_mother['total'])}）。H2009 分別為 {cn_mother['by_dataset']['H2009']['magnitude_pass']:,} / {cn_mother['by_dataset']['H2009']['stable_pass']:,}。這是 dataset-region engineering screen，HCC1395/DORADO 共用 CN 來源，不是 6 個獨立 biological samples 的 prevalence；也不能 orient HP1/HP2 或排除 subclonal CNA。
- 與本分析交集後，{counts['regions_direct_sister_shape_invariant_pair_complete']:,} 個 analysis-complete topology-invariant 中只有 {stable['topology_invariant_analysis_complete']} 個通過 CN magnitude + solution-stability screen；{counts['regions_direct_sister_tree_unique_pair_complete']:,} 個 analysis-complete tree-unique 為 {stable['tree_unique_analysis_complete']}；7 個 observed strict 為 {stable['observed_hp1_direct_only_hp2_sister_only']}。
- 唯一交集是 H2009 chr1:120007237-120040749：total=2.0329、minor=0.8958、major=1.1371。ranked output 目前只報一個 purity/ploidy solution，因此通過本次 heuristic；這不是整個 solution-space uniqueness 或正交 purity validation。該區仍只有 topology invariant，HP1 exact tree 有兩解。
- HCC1937 chr20:58489593-58518912 raw re-extraction：306 alignment exposures、299 QNAME、sidecar exact join 306/306、missing/conflict=0；baseline 加四個 denominator/PS sensitivity views 都維持 solver HP1 direct-only + HP2 sister-only，但 CN unavailable。
- H2009 候選 raw re-extraction：183 個至少在一個候選 sSNV 可判讀的 exposures（182 primary + 1 supplementary）、182 QNAME、sidecar exact join 183/183；baseline 加四個 sensitivity views 都維持 solver HP1 兩棵 exact trees／一種 direct-only shape + HP2 一棵 sister-only tree。低頻 `<3` states 全為 partial，沒有隱藏 full state 推翻判讀；observed thresholded full-state geometry 兩側則都是 single-only，direct/sister 來自 partial constraints 與 hidden-node solver。
- 上述 306 / 183 是 sSNV-state re-extraction exposures；下節的 H2009 251 是 methyl BAM 中 primary、MAPQ≥20、interval-overlap alignment identities，不要求 sSNV 可判讀，兩者不可互換分母。

## 9. Backbone sensitivity

- Analysis-complete topology invariant：canonical {counts['regions_direct_sister_shape_invariant_pair_complete']}、ClairS PASS {sets['direct_sister_shape_invariant_analysis_complete']['right_n']}；exact overlap {sets['direct_sister_shape_invariant_analysis_complete']['exact_intersection_n']}，Jaccard {sets['direct_sister_shape_invariant_analysis_complete']['jaccard']:.3f}。Discovery metric 另為 54 / 47。
- Analysis-complete tree-unique：{counts['regions_direct_sister_tree_unique_pair_complete']} vs {sets['direct_sister_tree_unique_analysis_complete']['right_n']}；overlap {sets['direct_sister_tree_unique_analysis_complete']['exact_intersection_n']}，Jaccard {sets['direct_sister_tree_unique_analysis_complete']['jaccard']:.3f}。Stored metric 另為 17 / 14。
- k=2 all-site collision：60 vs 69；overlap 48，Jaccard {sets['two_site_all_sites_cross_hp_collision']['jaccard']:.3f}。
- 固定 catalog inference-ready：兩個 backbone 都是 0。

## 10. Methyl 與 VAF

- 最新 canonical 保存 max_vaf 與 HP-specific REF/ALT counts，但沒有用 VAF 對候選樹做 ranking。
- VAF 是 joint read-state model 的 marginal；不應把同一批 reads 產生的 VAF 再當獨立分數重複加權。VAF/read fraction 是 molecule fraction，不是 cancer cell fraction（CCF）；轉成 CCF 仍需 purity、allele-specific CN、mutant multiplicity 與 normal contamination 校正。
- 54/54 topology-invariant candidates 的 bounded region 都能 exact join methyl BAM/HP sidecar；49/54 的 HP1/HP2 MM+ML rate ≥95%，29/54 兩個 HP 都有 full-span methyl reads。
- H2009 指定區域 primary MAPQ≥20 interval-overlap 為 251 alignment identities，exact join 251/251；HP1 94/94、HP2 96/96 有 MM+ML；full-span 10+5=15/15 具有區域內 modification entries，以 `alignment.modified_bases_forward` 解碼共 3,155 entries（不是 `modified_bases` 的 3,173）。
- 這只證明 methyl read-level data 可對齊；layered L3 仍是 not_evaluated，尚未做 K=1/K=2、matched-normal cis-ASM 或 clone clustering。
- KDE 對本 read/CN/methyl mixture 問題是 N/A，因為這裡不使用 Coverage_Multiple。
- Methyl 只能在 genotype/HP 固定後作 residual corroboration，不可單獨確認 clone 或 tree edge。

## 11. 下一輪最小驗證

1. 輸出未截斷 region_read_state：alignment id、QNAME、HP、PS、R/A/X、callable mask、品質與來源 digest。
2. 輸出 region_state_count：alignment-identity 與 primary-QNAME 兩種 denominator，full/partial 分開。
3. 以 H2009 做 primary 1+1 pilot；HCC1395/DORADO 作技術再現；H1437 作 CN-solution sensitivity。
4. 比較固定與替代 HP1×HP2 catalogs，使用 constrained MLE/EM、bootstrap與 abstain rule。
5. Methyl 只在 genetic-state anchored reads 內做 K=1 vs K=2 stability。
6. biological confirmation 需 single-cell、cell barcode、clone isolate、多區域或 longitudinal evidence。

## 12. 可直接用於報告的句子

> 在最新 LongPhase-S PASS canonical 的 7 個 datasets（6 個 biological samples）中，我們可以利用 ONT read-level sSNV 共現與 HP 分層，篩出跨 homolog 的 regional mutation-state candidates。51,815 個 regions 中有 22,779 個同時具有 HP1/HP2 mutation-bearing lineages；目前枚舉結果有 54 個呈現一側 direct、另一側 sister，其中 45 個 candidate-set complete；17 個目前各側只儲存一棵樹，其中 13 個才是 analysis-complete tree-unique。然而，沒有任何區域同時符合固定兩-sSNV五成分模型的完整反推條件，因此目前結果支持 regional mutation-state architecture candidates 的定點驗證，不支持已確認四個 cellular clones 或唯一 joint clone tree。

## 13. 證據與執行

輸入：

- /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5
- /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v6

輸出：

- InterSubMod/research/20260714_cross_hp_clone_state_inverse_audit/results/cross_hp_candidate_audit.json
- InterSubMod/research/20260714_cross_hp_clone_state_inverse_audit/results/candidate_regions.tsv
- InterSubMod/research/20260714_cross_hp_clone_state_inverse_audit/results/backbone_sensitivity_comparison.json
- InterSubMod/research/20260714_cross_hp_clone_state_inverse_audit/results/ascn_candidate_overlap.json
- InterSubMod/research/20260714_cross_hp_clone_state_inverse_audit/results/bounded_read_methyl_validation_receipt.json
- InterSubMod/research/20260714_cross_hp_clone_state_inverse_audit/results/independent_review_receipt.json

驗證：3/3 synthetic contract tests PASS；7/7 region conservation PASS；兩個 bounded candidate raw-read sensitivities不翻轉；canonical 與 sensitivity run 均 SUCCEEDED；3/3 independent regression reviewers PASS。
"""


def html_sample_rows(summaries: list[dict[str, Any]]) -> str:
    rows = []
    for summary in summaries:
        counts = summary["counts"]
        cn = summary["copy_number_contract"]
        cn_text = (
            f"{esc(cn['source'])}<small>coarse；非正式 1+1</small>"
            if cn["availability"] == "measured"
            else "<span class='muted'>unavailable</span>"
        )
        rows.append(
            f"""<tr>
<th scope="row">{esc(summary['dataset'])}<small>{esc(summary['biological_id'])}</small></th>
<td>{fmt(counts['regions_total'])}</td>
<td>{fmt(counts['regions_two_primary_hp'])}<small>{pct(counts['regions_two_primary_hp'], counts['regions_total'])} / total</small></td>
<td>{fmt(counts['regions_direct_sister_shape_invariant'])}<small>{fmt(counts['regions_direct_sister_shape_invariant_pair_complete'])} complete</small></td>
<td>{fmt(counts['regions_direct_sister_tree_unique'])}<small>{fmt(counts['regions_direct_sister_tree_unique_pair_complete'])} complete</small></td>
<td>{fmt(counts['regions_observed_hp1_direct_only_hp2_sister_only'])}</td>
<td>{cn_text}</td>
</tr>"""
        )
    return "".join(rows)


def sample_bar_svg(summaries: list[dict[str, Any]]) -> str:
    width, left, top, row_h, plot_w = 940, 132, 38, 48, 720
    maximum = max(row["counts"]["regions_total"] for row in summaries)
    height = top + row_h * len(summaries) + 38
    items = [
        f'<svg class="sample-bars" viewBox="0 0 {width} {height}" role="img" aria-label="各 dataset 全部與雙 primary HP regions">',
        '<rect x="108" y="8" width="14" height="10" rx="2" class="bar-total"/>',
        '<text x="132" y="18" class="svg-legend">全部 regions</text>',
        '<rect x="276" y="8" width="14" height="10" rx="2" class="bar-dual"/>',
        '<text x="300" y="18" class="svg-legend">雙 primary HP</text>',
    ]
    for index, summary in enumerate(summaries):
        y = top + index * row_h
        counts = summary["counts"]
        total_w = plot_w * counts["regions_total"] / maximum
        dual_w = plot_w * counts["regions_two_primary_hp"] / maximum
        items.extend(
            [
                f'<text x="0" y="{y + 17}" class="svg-label">{esc(summary["dataset"])}</text>',
                f'<rect x="{left}" y="{y}" width="{total_w:.2f}" height="22" rx="3" class="bar-total"/>',
                f'<rect x="{left}" y="{y + 7}" width="{dual_w:.2f}" height="8" rx="2" class="bar-dual"/>',
                f'<text x="{left + total_w + 8:.2f}" y="{y + 16}" class="svg-value">{fmt(counts["regions_total"])}</text>',
            ]
        )
    items.append("</svg>")
    return "".join(items)


def flow_svg(counts: dict[str, int]) -> str:
    total = counts["regions_total"]
    dual = counts["regions_two_primary_hp"]
    complete = counts["regions_pair_complete"]
    shape = counts["regions_direct_sister_shape_invariant_pair_complete"]
    unique = counts["regions_direct_sister_tree_unique_pair_complete"]
    any_alt = counts["regions_two_primary_hp_any_cross_hp_same_site_alt"]
    all_alt = counts["regions_two_primary_hp_all_sites_cross_hp_same_site_alt"]
    k2 = counts["regions_two_primary_hp_two_site_all_sites_cross_hp_same_site_alt"]
    intersect = counts["regions_direct_sister_shape_invariant_pair_complete_any_collision"]
    unique_intersect = counts["regions_direct_sister_tree_unique_pair_complete_any_collision"]
    return f"""<svg class="flow-svg" viewBox="0 0 1100 610" role="img" aria-labelledby="flow-title flow-desc">
<title id="flow-title">跨 HP 候選平行篩選</title><desc id="flow-desc">拓撲與同位點支持為兩條平行證據鏈，最後取交集。</desc>
<defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="7" refY="3" orient="auto"><path d="M0,0 L0,6 L8,3 z" class="arrow-head"/></marker></defs>
<text x="50" y="42" class="svg-kicker">共同母集合</text>
<rect x="50" y="62" width="250" height="105" rx="10" class="flow-root"/><text x="75" y="94" class="flow-label">全部 regions</text><text x="75" y="135" class="flow-number">{fmt(total)}</text>
<path d="M300 114 H380" class="flow-line" marker-end="url(#arrow)"/>
<rect x="390" y="62" width="270" height="105" rx="10" class="flow-root"/><text x="415" y="94" class="flow-label">雙 primary HP</text><text x="415" y="135" class="flow-number">{fmt(dual)}</text><text x="415" y="155" class="flow-note">{pct(dual,total)} / 全部</text>
<rect x="780" y="55" width="270" height="115" rx="10" class="flow-intersection"/><text x="805" y="83" class="flow-label">A ∩ B（非線性漏斗）</text><text x="805" y="112" class="flow-note">complete Topo + any：<tspan class="strong">{intersect}</tspan></text><text x="805" y="137" class="flow-note">complete unique + any：<tspan class="strong">{unique_intersect}</tspan></text><text x="805" y="160" class="flow-note">complete Topo + all：<tspan class="strong">0</tspan></text>
<text x="50" y="230" class="svg-kicker">A · SOLVER TOPOLOGY</text>
<rect x="50" y="250" width="250" height="92" rx="10" class="flow-a"/><text x="75" y="282" class="flow-label">兩側 candidate complete</text><text x="75" y="320" class="flow-number small">{fmt(complete)}</text><text x="185" y="320" class="flow-note">{pct(complete,dual)} / 雙 HP</text>
<path d="M300 296 H370" class="flow-line a" marker-end="url(#arrow)"/>
<rect x="380" y="250" width="250" height="92" rx="10" class="flow-a"/><text x="405" y="282" class="flow-label">complete direct+sister</text><text x="405" y="320" class="flow-number small">{shape}</text>
<path d="M630 296 H700" class="flow-line a" marker-end="url(#arrow)"/>
<rect x="710" y="250" width="250" height="92" rx="10" class="flow-a"/><text x="735" y="282" class="flow-label">complete tree-unique</text><text x="735" y="320" class="flow-number small">{unique}</text><text x="815" y="320" class="flow-note">{pct(unique,shape)} / Topo</text>
<text x="50" y="400" class="svg-kicker">B · SAME-SITE ALT</text>
<rect x="50" y="420" width="250" height="92" rx="10" class="flow-b"/><text x="75" y="452" class="flow-label">至少一位點雙 HP ALT</text><text x="75" y="490" class="flow-number small">{fmt(any_alt)}</text><text x="170" y="490" class="flow-note">{pct(any_alt,dual)} / 雙 HP</text>
<path d="M300 466 H370" class="flow-line b" marker-end="url(#arrow)"/>
<rect x="380" y="420" width="250" height="92" rx="10" class="flow-b"/><text x="405" y="452" class="flow-label">全部位點雙 HP ALT</text><text x="405" y="490" class="flow-number small">{all_alt}</text><text x="485" y="490" class="flow-note">{pct(all_alt,any_alt)} / any</text>
<path d="M630 466 H700" class="flow-line b" marker-end="url(#arrow)"/>
<rect x="710" y="420" width="250" height="92" rx="10" class="flow-b"/><text x="735" y="452" class="flow-label">其中 k=2</text><text x="735" y="490" class="flow-number small">{k2}</text><text x="815" y="490" class="flow-note">{pct(k2,all_alt)} / all-site</text>
</svg>"""


def candidate_card(title: str, record: dict[str, Any], verdict: str, reason: str) -> str:
    return f"""<article class="candidate-card">
<div class="candidate-head">{esc(record['dataset'])}<code>{esc(record['region'])}</code></div>
<h4>{esc(title)}</h4><dl>
<div><dt>k</dt><dd>{record['n_sSNV']}</dd></div><div><dt>CN</dt><dd>{esc(record['cn_label'])}</dd></div>
<div><dt>PS</dt><dd>{record['n_unique_phase_sets']}</dd></div><div><dt>collision</dt><dd>{len(record['collision_positions'])}/{record['n_sSNV']}</dd></div>
<div><dt>HP1</dt><dd>{esc(', '.join(record['hp1_shape_classes']))}</dd></div><div><dt>HP2</dt><dd>{esc(', '.join(record['hp2_shape_classes']))}</dd></div>
<div><dt>Observed HP1</dt><dd>{esc(record['observed_full_geometry_hp1']['geometry_class'])}</dd></div><div><dt>Observed HP2</dt><dd>{esc(record['observed_full_geometry_hp2']['geometry_class'])}</dd></div>
</dl><p class="verdict-line"><strong>{esc(verdict)}</strong> — {esc(reason)}</p>
<details><summary>thresholded full-read states</summary><pre>HP1 {esc(json.dumps(record['hp1_full_counts_retained'], ensure_ascii=False, sort_keys=True))}
HP2 {esc(json.dumps(record['hp2_full_counts_retained'], ensure_ascii=False, sort_keys=True))}</pre></details></article>"""


def make_html(
    audit: dict[str, Any], comparison: dict[str, Any], ascn: dict[str, Any], css: str
) -> str:
    counts = audit["aggregate_counts"]
    summaries = audit["sample_summaries"]
    records = audit["candidate_records"]
    sets = comparison["exact_region_set_comparisons"]
    stable = ascn["category_stable_pass_counts"]
    cn_mother = ascn["dual_hp_mother_set_screen"]
    hcc1937 = find_record(records, "HCC1937", "chr20:58489593-58518912")
    h2009 = find_record(records, "H2009", "chr1:120007237-120040749")
    hcc1395 = find_record(records, "HCC1395", "chr17:21303695-21320366")
    created = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
    return f"""<!doctype html>
<html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="description" content="跨 HP clone-state 反推可行性、7-dataset 候選觀察與驗證紀錄">
<title>跨 HP clone-state 反推｜可行性與觀察紀錄</title><style>{css}</style></head>
<body><a class="skip" href="#main">跳到主要內容</a><div id="progress" aria-hidden="true"></div>
<header class="hero"><div class="hero-grid"><div><div class="stamp">InterSubMod · Comprehensive validation · 2026-07-14</div>
<h1>跨 HP clone-state<br>能反推到哪裡？</h1>
<p>以最新 7-dataset canonical，分開檢查 read-state、solver 拓撲、跨 HP 同位點與固定 clone catalog；把能觀察與已確認切開。</p></div>
<aside class="hero-aside"><strong>0</strong><small>符合完整五成分反推門檻的真實 regions</small>
<div class="hero-note">不是沒有 subclone。<br>是目前沒有一個區域同時通過所有可識別性、CN、PS 與 read-state gates。</div></aside></div></header>
<div class="layout"><nav aria-label="報告目錄"><div class="nav-title">Observation ledger</div><ol>
<li><a href="#answer">01 直答結論</a></li><li><a href="#scope">02 範圍與單位</a></li>
<li><a href="#funnel">03 平行證據鏈</a></li><li><a href="#definitions">04 定義</a></li>
<li><a href="#samples">05 樣本數據</a></li><li><a href="#inverse">06 反推模型</a></li>
<li><a href="#cases">07 候選個案</a></li><li><a href="#sensitivity">08 敏感度</a></li>
<li><a href="#limits">09 限制與下一步</a></li><li><a href="#claim">10 報告句子</a></li>
<li><a href="#evidence">11 證據鏈</a></li></ol></nav><main id="main">

<section id="answer"><div class="section-head"><div class="section-no">01 / ANSWER</div><h2>可以篩候選；<br>還不能確認四個 clones。</h2></div>
<p class="lead">現有資料已足以做 region × HP mutation-state census、找出 direct/sister 拓撲與同位點跨 HP ALT；但 bulk reads 沒有 cell pairing，未截斷 state table與逐區域 methyl 驗證也尚未成為 canonical，因此 biological clone tree 仍是 PROBE。</p>
<div class="verdict-grid"><article class="verdict go"><span class="tag">GO · VERIFIED</span><strong>候選篩選</strong><p>7/7 canonical 完成；MINREAD≥3 的 thresholded R/A/X states、counts、trees 可重跑；低支持 raw states 由 bounded BAM re-extraction 補查。</p></article>
<article class="verdict probe"><span class="tag">PROBE · CONDITIONAL</span><strong>比例反推</strong><p>固定五成分 catalog 時可識別；真實資料需 constrained likelihood 與替代模型比較。</p></article>
<article class="verdict no"><span class="tag">NO-GO · CURRENT CLAIM</span><strong>clone 確認</strong><p>不能把兩個 HP marginals 宣稱為同一細胞的 homolog pairing 或唯一 joint clone tree。</p></article></div>
<div class="kpis"><div class="kpi"><strong>{fmt(counts['regions_total'])}</strong><span>全部 regions</span><em>chr1–22 · 7 datasets</em></div>
<div class="kpi"><strong>{fmt(counts['regions_two_primary_hp'])}</strong><span>雙 primary HP</span><em>{pct(counts['regions_two_primary_hp'],counts['regions_total'])} / total</em></div>
<div class="kpi"><strong>{counts['regions_direct_sister_shape_invariant']}</strong><span>目前 shape-invariant</span><em>{counts['regions_direct_sister_shape_invariant_pair_complete']} analysis-complete</em></div>
<div class="kpi"><strong>{counts['regions_direct_sister_tree_unique']}</strong><span>目前各側單樹</span><em>{counts['regions_direct_sister_tree_unique_pair_complete']} analysis-complete</em></div></div>
<div class="callout"><strong>最容易誤解：</strong> 54 是 solver 拓撲集合；86 是 all-site cross-HP ALT 集合。兩者為平行證據，不是 86 再篩成 54。</div></section>

<section id="scope"><div class="section-head"><div class="section-no">02 / SCOPE</div><h2>先固定母集合，<br>每個數字才有意義。</h2></div>
<div class="two-col"><article class="definition"><h3>資料範圍</h3><p><b>Canonical</b>：LongPhase-S recalibrated FILTER=PASS，chr1–22，run v5，SUCCEEDED。</p>
<p>638,259 raw/all → 582,820 LPS PASS → 469,849 autosomal biallelic sSNVs → 194,149 retained → 51,815 regions。</p></article>
<article class="definition"><h3>7 不等於 7 個病人</h3><p><b>7 datasets / 6 biological samples</b>。HCC1395 與 HCC1395_DORADO 是同一細胞系兩套資料，不可在 biological denominator 重複計數。</p></article></div>
<details><summary>為何玩具例是 k=2，而不是 k=3？</summary><div><p>G1 是 germline SNP，用來決定 HP family；樹的 state vector 只放 somatic sSNVs。因此 G1-S1-S2 的正式 query 是 k=2。</p></div></details>
<details><summary>Region、state tree、clone 的差別</summary><div><p>Region 是 sSNV connected component；HP lineage 是依 germline HP family 分層的 evidence；state 是 R/A/X；tree 是符合 constraints 的 mutation-state topology；clone 是 cell population。前四者可算，第五者不能直接觀察。</p></div></details></section>

<section id="funnel"><div class="section-head"><div class="section-no">03 / EVIDENCE</div><h2>兩條證據鏈，<br>最後才取交集。</h2></div>
<div class="flow-wrap">{flow_svg(counts)}</div><p class="lead">Topology 描述候選演化幾何；same-site collision 描述同一 genomic sSNV 在兩個 HP 都有 ALT evidence。只有交集才接近固定跨 HP clone catalog。</p></section>

<section id="definitions"><div class="section-head"><div class="section-no">04 / DEFINITIONS</div><h2>同一句 direct+sister，<br>其實有兩個 estimator。</h2></div>
<div class="two-col"><article class="definition"><h3>A · Solver topology</h3><p>目前枚舉結果中，一側全為 direct-only、另一側全為 sister-only，稱 shape-invariant。只有兩側 candidate set complete 時才可作 exhaustive 解讀。</p><strong>{counts['regions_direct_sister_shape_invariant']} enumerated / {counts['regions_direct_sister_shape_invariant_pair_complete']} complete · {counts['regions_direct_sister_tree_unique']} stored-single / {counts['regions_direct_sister_tree_unique_pair_complete']} complete-unique</strong></article>
<article class="definition"><h3>B · Observed full-state geometry</h3><p>只看 thresholded full reads。mutation sets 有 proper subset 稱 direct；互不包含稱 sister。</p><strong>{counts['regions_observed_hp1_direct_hp2_sister_broad']} broad · {counts['regions_observed_hp1_direct_only_hp2_sister_only']} strict</strong></article></div>
<div class="callout">兩個 estimator 不能混稱 clone 數。Solver 納入 partial subcube constraints；observed geometry 只看完整且 count≥3 的 states。</div>
<details><summary>Cross-HP same-site collision</summary><div><p>同一 sSNV position 在 HP1 與 HP2 的 ALT count 都 ≥3。它是候選旗標，仍可能來自 LOH/gain、HP switch、mapping artifact 或 caller error。</p></div></details>
<details><summary>single-PS 的精確意思</summary><div><p>只表示非缺值 PS 的 distinct count=1，不表示每條 alignment 都帶 PS，也不證明跨區域 phase continuity。H2009 HP1 為 64/75、HP2 為 77/77；HP1 缺 PS 的 11 條全是 XXA，移除後 major-PS sensitivity 不翻轉。</p></div></details>
<details><summary>為何 populations 與 subreads 不能相加？</summary><div><p>每條 read 先進 subread_groups，完整覆蓋時又進 populations；full key 會重複呈現。反推時 full counts 只能取 obs_populations。</p></div></details></section>

<section id="samples"><div class="section-head"><div class="section-no">05 / DATASETS</div><h2>每個 dataset 都有數，<br>但 CN 能力不相同。</h2></div>
{sample_bar_svg(summaries)}
<div class="table-wrap"><table><caption>7 datasets 同口徑 census；Topo 與單樹欄為目前枚舉值，small text 為 analysis-complete</caption><thead><tr><th>Dataset / biological ID</th><th>Regions</th><th>雙 HP</th><th>Topo enumerated</th><th>單樹 stored</th><th>Observed strict</th><th>Canonical CN</th></tr></thead><tbody>{html_sample_rows(summaries)}</tbody></table></div>
<p class="callout"><strong>Canonical neutral 不等於 allele-specific 1+1。</strong> External SAVANA 可做 magnitude screen，但 major/minor 沒有 orient 到 HP1/HP2，也不能排除 subclonal CNA。</p></section>

<section id="inverse"><div class="section-head"><div class="section-no">06 / INVERSE</div><h2>固定 catalog 時可解；<br>放開 catalog 就不唯一。</h2></div>
<div class="formula"><pre>direct HP      sister HP
Background RR  RR
C1      I      RR
C2      AA     RR
C3      RR     U
C4      I      V

π₂ = pᴅ(AA)
π₃ = pₛ(U)
π₄ = pₛ(V)
π₀ = pᴅ(RR) − pₛ(U)
π₁ = pᴅ(I)  − pₛ(V)</pre>
<div class="assumption"><h3>成立前提</h3><p>固定五種 components、CN=1+1、無 HP error、無 state-specific bias、使用未截斷 per-molecule counts。</p><p>合成測試 design rank 5/5；π=(.10,.20,.25,.15,.30) 精確回復。π₀ 是 RR|RR background；無外部 purity prior 時不能直接叫 normal fraction。</p></div></div>
<div class="matrix"><div class="yes"><strong>數學模型</strong><p>固定 catalog 可條件式識別。</p></div><div class="conditional"><strong>真實估計</strong><p>需 MLE/EM、bootstrap與 alternative catalog comparison。</p></div><div class="nope"><strong>Cell pairing</strong><p>HP marginals 無法產生同一細胞的 homolog 配對。</p></div></div>
<h3>為何真實資料是 0？</h3><p>{counts['regions_direct_sister_shape_invariant_pair_complete']} 個 analysis-complete direct+sister invariant regions 的 k 從 3 起；60 個 k=2 all-site collision regions 多是 direct+direct 或 sister+sister。兩條證據鏈尚未在固定玩具 catalog 交會。</p></section>

<section id="cases"><div class="section-head"><div class="section-no">07 / CASES</div><h2>三個最接近的區域，<br>各自卡在不同 gate。</h2></div>
<div class="candidate-grid">
{candidate_card('結構最穩健',hcc1937,'PROBE','Raw re-extraction與四種 sensitivity不翻轉；但 CN unavailable。')}
{candidate_card('唯一 asCN 交集',h2009,'PROBE','Solver HP1 direct-only同形兩樹、HP2 sister-only一樹；observed full-state兩側皆single-only。')}
{candidate_card('Solver 單樹 + collision',hcc1395,'NO-GO for fraction','Solver HP2為sister-only，但displayed full-state geometry為single-only；且CN=LOH。')}
</div>
<p class="callout"><strong>分母不可混用：</strong> H2009 的 183 是至少一個候選 sSNV 可判讀的 alignment exposures（182 primary + 1 supplementary）；251 是 methyl BAM 的 primary、MAPQ≥20、interval-overlap alignment identities，不要求 sSNV 可判讀。</p>
<h3>External allele-specific CN screen</h3>
<p>22,779 個 dataset-regions 中，{cn_mother['magnitude_pass']:,}（{pct(cn_mother['magnitude_pass'],cn_mother['total'])}）通過宣告的 1+1 magnitude gate；加入目前 heuristic stability 後為 {cn_mother['stable_pass']:,}（{pct(cn_mother['stable_pass'],cn_mother['total'])}）。H2009 為 {cn_mother['by_dataset']['H2009']['magnitude_pass']:,} / {cn_mother['by_dataset']['H2009']['stable_pass']:,}。HCC1395/DORADO 共用 CN，這不是獨立 biological prevalence。與候選交集後：{counts['regions_direct_sister_shape_invariant_pair_complete']} complete topology 僅 {stable['topology_invariant_analysis_complete']} 個；{counts['regions_direct_sister_tree_unique_pair_complete']} complete unique 為 {stable['tree_unique_analysis_complete']}；7 observed strict 為 {stable['observed_hp1_direct_only_hp2_sister_only']}。</p>
<details><summary>H2009 唯一交集的 CN 數據</summary><div><p>chr1_seg20：108,270,001–121,500,000；total 2.0329、minor 0.8958、major 1.1371；13.23 Mb、1,323 bins、6,901 het SNPs。ranked output 只報一列（0.95 / 2.22），因此通過本次 heuristic；不是全 solution-space uniqueness 或 orthogonal purity validation。這仍是 magnitude screen，不提供 HP orientation，也不能排除 subclonal CNA。</p></div></details></section>

<section id="sensitivity"><div class="section-head"><div class="section-no">08 / SENSITIVITY</div><h2>結論 0 穩健，<br>個案 membership 仍會變。</h2></div>
<div class="sensitivity"><article><h3>Complete topology</h3><div class="pair"><span>Canonical</span><strong>{counts['regions_direct_sister_shape_invariant_pair_complete']}</strong></div><div class="pair"><span>ClairS PASS</span><strong>{sets['direct_sister_shape_invariant_analysis_complete']['right_n']}</strong></div><p>overlap {sets['direct_sister_shape_invariant_analysis_complete']['exact_intersection_n']} · Jaccard {sets['direct_sister_shape_invariant_analysis_complete']['jaccard']:.3f}<br>discovery metric 54 / 47</p></article>
<article><h3>Complete tree-unique</h3><div class="pair"><span>Canonical</span><strong>{counts['regions_direct_sister_tree_unique_pair_complete']}</strong></div><div class="pair"><span>ClairS PASS</span><strong>{sets['direct_sister_tree_unique_analysis_complete']['right_n']}</strong></div><p>overlap {sets['direct_sister_tree_unique_analysis_complete']['exact_intersection_n']} · Jaccard {sets['direct_sister_tree_unique_analysis_complete']['jaccard']:.3f}<br>stored metric 17 / 14</p></article>
<article><h3>k=2 all-site</h3><div class="pair"><span>Canonical</span><strong>60</strong></div><div class="pair"><span>ClairS PASS</span><strong>69</strong></div><p>overlap 48 · Jaccard {sets['two_site_all_sites_cross_hp_collision']['jaccard']:.3f}</p></article></div>
<p class="callout"><strong>解讀：</strong> 固定 catalog inference-ready 在兩個 backbones 都是 0；analysis-complete exact Jaccard 約 0.53–0.56，單一 region 仍需 raw-read 複核。</p></section>

<section id="limits"><div class="section-head"><div class="section-no">09 / LIMITS</div><h2>缺的不是更多圖，<br>是四個資料契約。</h2></div>
<ol class="steps"><li><div><strong>未截斷 read-state table</strong><p>保留 count&lt;3、alignment identity、QNAME、HP、PS、state、callable mask與品質。</p></div></li>
<li><div><strong>Molecule denominator sensitivity</strong><p>同時報 alignment-identity 與 primary-QNAME，full/partial 分開。</p></div></li>
<li><div><strong>CN 與 HP error</strong><p>整合 1+1 magnitude、purity、LOH、switch/confusion；否則只報 molecule fractions。</p></div></li>
<li><div><strong>Model comparison + orthogonal evidence</strong><p>固定與替代 catalogs 做 likelihood/bootstrapping；methyl只作 anchored auxiliary；最終需 single-cell或多區域 evidence。</p></div></li></ol>
<details><summary>VAF 放在哪裡？</summary><div><p>VAF 是 joint read-state model 的 marginal。候選模型應同時預測 per-site VAF並校準；不要把同一批 reads導出的 VAF再當獨立分數重複加權。VAF/read fraction 是 molecule fraction，不是 cancer cell fraction；轉成 CCF 需 purity、allele-specific CN、mutant multiplicity 與 normal contamination 校正。Canonical目前沒有 per-tree VAF ranking。</p></div></details>
<details><summary>Methyl 目前能做到哪裡？</summary><div><p>54/54 topology candidates 均可 bounded exact join；49/54 兩個 HP 的 MM/ML rate ≥95%，29/54 兩側都有 full-span methyl reads。H2009 指定區域 251/251 exact join，HP1 94/94、HP2 96/96 有 MM+ML；full-span 15/15 有區域內 modification entries，以 <code>alignment.modified_bases_forward</code> 解碼 3,155 entries（HP1 2,086 + HP2 1,069）。這只證明資料可對齊；layered L3 仍是 not_evaluated。穩定 methyl K=2 也只能叫 candidate epigenetic subpopulation，不能單獨確認 clone。</p></div></details></section>

<section id="claim"><div class="section-head"><div class="section-no">10 / CLAIM</div><h2>可直接用於教授報告的句子。</h2></div>
<p class="claim">在最新 LongPhase-S PASS canonical 的 7 個 datasets（6 個 biological samples）中，我們可以利用 ONT read-level sSNV 共現與 HP 分層，篩出跨 homolog 的 regional mutation-state candidates。51,815 個 regions 中有 22,779 個同時具有 HP1/HP2 mutation-bearing lineages；目前枚舉結果有 54 個呈現一側 direct、另一側 sister，其中 45 個 candidate-set complete；17 個目前各側只儲存一棵樹，其中 13 個才是 analysis-complete tree-unique。然而，沒有任何區域同時符合固定兩-sSNV五成分模型的完整反推條件，因此目前結果支持 regional mutation-state architecture candidates 的定點驗證，不支持已確認四個 cellular clones 或唯一 joint clone tree。</p></section>

<section id="evidence"><div class="section-head"><div class="section-no">11 / PROVENANCE</div><h2>每個數字有來源，<br>每個限制有出口。</h2></div>
<div class="evidence-list"><p><strong>Primary input</strong><br>/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5</p>
<p><strong>Sensitivity input</strong><br>/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v6</p>
<p><strong>Machine-readable outputs</strong><br>../results/cross_hp_candidate_audit.json<br>../results/candidate_regions.tsv<br>../results/backbone_sensitivity_comparison.json<br>../results/ascn_candidate_overlap.json<br>../results/bounded_read_methyl_validation_receipt.json<br>../results/independent_review_receipt.json</p>
<p><strong>Verification</strong><br>3/3 synthetic tests PASS · 7/7 region conservation PASS · two bounded raw-read candidate sensitivities不翻轉 · canonical/sensitivity SUCCEEDED · 3/3 independent regression reviewers PASS</p></div>
<p>產生時間：{esc(created)}。這是 observation record，不改寫 canonical pipeline。</p></section>
</main></div><footer><div><strong>Claim ceiling</strong><p>Candidate cross-HP regional mutation-state pattern；not confirmed cell-level HP pairing, clone identity, or biological lineage.</p></div></footer>
<script>const bar=document.getElementById("progress");const update=()=>{{const h=document.documentElement;const d=h.scrollHeight-h.clientHeight;bar.style.width=(d?100*h.scrollTop/d:0)+"%"}};addEventListener("scroll",update,{{passive:true}});update();</script>
</body></html>"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--ascn", required=True, type=Path)
    parser.add_argument("--style", required=True, type=Path)
    parser.add_argument("--out-md", required=True, type=Path)
    parser.add_argument("--out-html", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    audit = load_json(args.audit)
    comparison = load_json(args.comparison)
    ascn = load_json(args.ascn)
    css = args.style.read_text(encoding="utf-8")
    markdown = make_markdown(audit, comparison, ascn)
    standalone = make_html(audit, comparison, ascn, css)
    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_html.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.write_text(markdown, encoding="utf-8")
    args.out_html.write_text(standalone, encoding="utf-8")
    print(f"WROTE {args.out_md.resolve()} ({args.out_md.stat().st_size:,} bytes)")
    print(f"WROTE {args.out_html.resolve()} ({args.out_html.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
