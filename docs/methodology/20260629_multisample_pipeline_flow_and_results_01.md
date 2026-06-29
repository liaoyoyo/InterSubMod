---
title: 7 ONT 樣本 clone/subclone 全流程 + 基因註釋 + 多分頁工作站 — 流程與結果
date: 2026-06-29
status: in_progress
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/multisample_summary.json, docs/methodology/_assets/20260629_multisample_topology_workstation.standalone.html
evidence_tier: L2 (7 樣本機械重建 + 跨樣本一致;單 pipeline 封頂 ⭐3)
---

# 7 ONT 樣本 clone/subclone 全流程與結果

## TL;DR（1 行）
7 個 ONT 癌症細胞株（ClairS full-model + LongPhase-S tag）全跑到 clone/subclone 樹結構，每區補 COSMIC/GENCODE/DGIdb 基因註釋，整合成單一多分頁互動工作站 HTML。

## 摘要（3 行）
- **為何**：把 HCC1395 單樣本 subclone 重建擴到全部 7 樣本，驗證 sSNV 共現骨幹方法的跨樣本/跨癌種一致性，並為每個區域補基因意義（名稱/啟動子/oncogene-TSG/用藥）。
- **怎麼做**：複用既有 8-step pipeline（per-sample 化 + 48 核平行 orchestrator），下游接 COSMIC CGC + Cell Lines Project + GENCODE + DGIdb 註釋，HTML 用 data-swap 分頁複用 20260628 工作站。
- **結果**：7 樣本全產出（with-vector 區 1,636–4,768；A_determined 比例 50–66%），cycle（incompatible）跨樣本 0–265（與 CN-gain multiplicity 相關），癌基因區 60–236、用藥區 264–983。

---

## §1 Pipeline 流程（8 step + 註釋 + HTML）

```
輸入(per sample):
  tagged tumor BAM (canonical/.../longphase_s/{S}_tagged.bam, ClairS full + LongPhase-S somatic tag)
  normal BAM (/big8_disk/data/{S}/...)
  ClairS VCF (filtered_snv tp/fp) + germline phased VCF
        │
        ▼
[1] sm_linkage_genomewide.py  (per-chrom ×22, 平行)
        sSNV × HP × read 單分子共現掃描 → part_chr{1..22}.json
        │
        ▼
[2] sm_merge_parts.py         合併 part_chr*.json → sm_linkage_genomewide.json (pairs + census)
        │
        ▼
[3] sm_pair_lists.py          2×2 (RR/RA/AR/AA) 配對統計 → sm_pair_lists_counts.json
        │
        ▼
[4] sm_multilocus_combinations.py  (per-chrom ×22, 平行)
        多位點基因型向量族群 → ml_part_chr{1..22}.json
        │
        ▼
[5] sm_region_integration.py  pairs + ml_part_* + CN → 區域整合 → sm_region_integration.json
        │
        ▼
[6] topology_analysis.py      每區拓樸判定 → topology_per_region.json
        determinacy: A_determined / A_ambiguous / B_pairwise / C_underdetermined
        topology_type: single / linear / branched / germline_only / no_genotype_vectors
        has_cycle (incompatible = 上游 pairwise graph 有環, 多為 CN-gain multiplicity)
        │
        ▼
[7] candidate_scoring.py      每區評分 + 佇列 3 欄(why_conflict/parsimony/methyl_applicability) → candidate_scoring.json
        │
        ▼
[8] region_gene_annotation.py 每區 → 基因名/啟動子/癌症角色/用藥 → region_gene_annotation.json
        GENCODE v46(基因 body + TSS±2kb 啟動子) + DGIdb(用藥) + COSMIC CGC v104(oncogene/TSG/tier)
        │
        ▼
build_topology_workstation.py  7 樣本 → 單一多分頁 HTML(data-swap, 複用 20260628 工作站)
```

### 平行化（48 核）
`run_multisample_parallel.py`：4 階段（P1 linkage 全域池 → P2 merge+pair_lists per-sample → P3 multilocus 全域池 → P4 region+topology+candidate+gene per-sample），per-chrom（22-way）×4 樣本 = 88 jobs、40 並發。**修正了原本沿用 5-way split 只用 5–7 核的低度利用**（load 7→46，~6-8× 加速）。

---

## §2 7 樣本結果彙總（數字源：`multisample_summary.json`）

| 樣本 | 總區 | 有向量區 | A_det | branched | linear | cycle | 癌基因區 | 用藥區 | 佇列 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| HCC1395 | 6,772 | 3,885 | 1,812 | 1,113 | 754 | 12 | 159 | 623 | 2,118 |
| COLO829 | 7,063 | 3,786 | 1,231 | 1,076 | 158 | 1 | 132 | 565 | 2,588 |
| H1437 | 8,317 | 4,768 | 2,161 | 1,431 | 775 | 65 | 160 | 799 | 2,708 |
| H2009 | 9,605 | 4,243 | 2,049 | 1,403 | 909 | 265 | 236 | 983 | 2,471 |
| HCC1395_DORADO | 3,807 | 2,379 | 1,206 | 707 | 516 | 4 | 92 | 343 | 1,316 |
| HCC1937 | 1,942 | 1,636 | 1,017 | 613 | 425 | 0 | 60 | 264 | 634 |
| HCC1954 | 3,095 | 1,979 | 1,148 | 920 | 234 | 1 | 63 | 265 | 851 |

- **canonical 分母**：determinacy 比例算在「有向量區」上，非「總區」（無向量區 = 無基因型向量、不可建樹）。
- **A_determined 比例**：47–62%（HCC1937 最高 62%、H1437/COLO829 較低）— 單分子向量唯一決定拓樸的區。
- **cycle（incompatible）**：H2009 265 最高、HCC1937 0 最低 — incompatible 多為 CN-gain multiplicity artifact（非真衝突），與樣本 CN 負荷相關。

---

## §3 問題與修復紀錄

### 🔴 P-1 跨-run merge 污染（已修復 + 量化）
- **症狀**：H2009 merge 時 `sm_linkage_genomewide.json` = 325MB（其餘樣本 ~200MB）。
- **根因**：先前的舊序列 run（`run_all_remaining_samples.sh`, 5-way split）卡在 H2009 part_2（311G BAM, 跑 4.5hr 未完）。改用新平行 orchestrator（per-chrom）後，**舊 5-way 進程未被殺、繼續在背景把 `part_1/3/4/5.json`（5-way 命名）寫回 H2009 目錄**。merge 的 glob `sm_linkage_gw_part_*.json` 同時吃到 5-way（part_N）+ per-chrom（part_chrN）兩種命名 → chr1-22 部分被重複計入。
- **修復**：殺殘留 5-way 進程 → 移除 5-way 污染檔（保留乾淨 part_chr*）→ H2009 從 merge 起乾淨重做。
- **驗證**：乾淨 merge = **203.3MB**（vs 污染 325MB，移除 122MB / 38% 重複）；重做後 4,243 區全唯一（0 重複）、ml_part 22/22、無 error。
- **教訓**：並行/序列 run 切換時，舊 run 的背景 compute 進程必須先確認終止；merge glob 應綁定唯一命名 pattern 避免跨-run 混入。

### P-2 5-way 核心低度利用（已修，見 §1 平行化）

---

## §4 COSMIC 資料整合

- **Cancer Gene Census v104**（768 基因）→ 每區 oncogene/TSG/tier/相關癌種標籤（已接入 step 8）。
- **Cell Lines Project v104**：6 樣本全有外部真值（命名對映 `COLO829→COLO-829` / `H1437→NCI-H1437` / `H2009→NCI-H2009`），per-gene CN（`TOTAL_CN`/`MINOR_ALLELE`）+ 突變。可補目前 cn=unknown + 交叉驗證區域突變（待整合為註釋層）。
- 其他下載：StructuralVariants/Breakpoints（SV 軸非循環錨）、MutantCensus、Fusion、Resistance、Classification、Genes、CompleteDifferentialMethylation（甲基軸）。
- 存放：repo 外 `/big7_disk/liaoyoyo2001/gene_annotation/`（GENCODE/DGIdb/CGC）+ `cosmic_v104/`（其餘）。

---

## §5 誠實邊界

- **證據層級 ⭐3（single-pipeline 封頂）**：7 樣本同一 pipeline、同一參數，跨樣本一致是**內部一致性**非獨立平台確認；升 ⭐4 需 single-cell/multi-region 正交。
- **甲基 = bounded-auxiliary**：骨幹是 sSNV 共現；甲基僅負向篩選 + L3 旗標，不參與本流程的分群/排序（見 `project_methylation_use_exhausted_bounded_auxiliary`）。
- **CN=unknown**：本輪未跑 SAVANA；incompatible/multiplicity 判讀力下降（誠實標 cn=unknown）。COSMIC CLP CN 為 gene-level（SNP-array），可作交叉檢核，pipeline 級 segment CN 仍待 SAVANA。
- **truth set 各樣本不同**（NYGC/SEQC2/GIAB）→ TP/FP 標記僅觀察、不進前處理。
- **reconstruction = regional partition 非單一腫瘤系統發生樹**；每區獨立判定拓樸，不跨區串成全基因組譜系。

## 輸出檔
- 多分頁工作站：`InterSubMod/docs/methodology/_assets/20260629_multisample_topology_workstation.standalone.html`（32MB, 7 分頁）
- 彙總：`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/multisample_summary.json`
- per-sample：`/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone/{S}/{topology_per_region,candidate_scoring,region_gene_annotation}.json`
