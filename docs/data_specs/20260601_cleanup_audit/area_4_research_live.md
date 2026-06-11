---
title: Cleanup Audit Area [4/6] — research/ LIVE projects
status: in_progress
owner: cleanup-audit-agent
created: 2026-06-01
scope: research/ active-project subset (read-only audit; NO delete/move)
sizing_note: per-project du; paired_priority_bug_audit total = timeout(large) due to genome-wide per-region ISM trees
---

# Area [4/6] — research/ 現役專案 清理稽核

唯讀稽核。**未刪除/搬移任何檔案**（刪/移為 Hard Gate，需用戶確認後由主 agent 執行）。

## 最高原則套用
- 結案 NEGATIVE ≠ 可刪。SAFE_DELETE 僅給重跑可重生的中間檔 / 明確重複副本 / 純 transient。
- 撐起已發佈結論的原始證據 → 至少 ARCHIVE。不確定 → NEEDS_USER_DECISION（保守）。
- `autoresearch / data_registry / _template` = **KEEP，零清理建議**（活躍狀態機 + SoT；CLAUDE.md `research/` 規範列為「不可 git rm 清單」）。

## 母專案 LIVE 判定（讀 manifest/README/00_*.md 交叉確認）

| 專案 | bytes | 母專案判定 | 依據 |
|------|------:|-----------|------|
| tpfp_loh_af_kde_discrimination | 33,493,123 | **LIVE** KEEP | 00_INDEX status=in_progress；论文主軸 LOH×AF×KDE phasing |
| loh_subclone_af | 4,831,420 | **LIVE** KEEP | manifest status=completed/POSITIVE；memory paired 層保留 |
| loh_subclone_af_paired | 4,698,490 | **LIVE** KEEP | manifest 同上；paired 對照 |
| loh_investigation | 46,293,563 | **LIVE** KEEP | LOH 判定機制研究區塊（论文主軸）|
| germline_asm_analysis | 2,675,182 | **LIVE** KEEP | ZAR1L/BRCA2 ASM ⭐3 characterization（全 figures）|
| tsg_promoter_asm_reviewer | 68,300,374 | **LIVE** KEEP | ASM ⭐3 reviewer survey（memory: validated v2 TSV）|
| B2_loh_skeptical_check | 40,400 | **LIVE** KEEP | ASM ⭐3 skeptical control（tiny）|
| selfphasing_v6_production | 42,512 | **LIVE** KEEP | V6 production manifest（指向外部 binary）|
| v6_bam_tpfp_hp_loh_cn | 203,106,032 | **LIVE** KEEP | V6 BAM TP/FP × HP×LOH×CN characterization |
| paired_priority_bug_audit | timeout(large) | **LIVE** KEEP | V6 priority-bug audit（撐 V6 validated 報告）|
| hku_collaboration | 19,847,201 | **LIVE** KEEP | HKU handoff 交付紀錄（5/24 已交付，保留）|
| thread_d_paper | 30,858 | **LIVE** KEEP | TP-enriched phasing tool paper draft（论文主軸 G6）|
| v5_provenance_followup | 2,260,614,371 | **LIVE** KEEP（results）| CLAUDE.md 列為「不可 git rm」validated project 範例 |
| autoresearch | 239,731 | **KEEP（狀態機）** | evidence_ledger.jsonl + hypothesis_queue.json SoT |
| data_registry | 11,994 | **KEEP（狀態機）** | 跨 project 共用資料索引 |
| _template | 19,750 | **KEEP（狀態機）** | init-research 範本 |

**所有母專案 LIVE 確認 → 整專案 KEEP。** 以下表格僅針對 **專案內部可回收的 intermediate/重複/superseded** 標 ARCHIVE / NEEDS_USER_DECISION（即使母專案 LIVE）。

---

## Artifact 級判定表（專案內可回收項）

| path | 大小 | bytes | purpose | trust_tier | conclusion_status | cleanup_verdict | reclaimable_bytes | referenced_by | rationale |
|------|-----:|------:|---------|-----------|-------------------|-----------------|------------------:|---------------|-----------|
| research/v5_provenance_followup/T1_2_read_level_audit/vote_dump_baseline_genome.tsv.gz | 744M | 779,000,000* | baseline 全基因組 read-level vote dump（priority-bug audit 原始輸入）| PRE-FIX | CONCLUDED_POSITIVE | **NEEDS_USER_DECISION→ARCHIVE** | 2,219,541,005 (三檔合計) | T1_2_priority_bug_mechanism_report.md §10（重現法）；validated 報告引用「結論」非「dump」 | audit 已定案（4-path PASS, V3F/V5 修正率 100%）；gitignored 純中間檔，可由 3 binaries 重生；2026-05-10 cleanup proposal 已建議封存待用戶確認 finalize |
| research/v5_provenance_followup/T1_2_read_level_audit/vote_dump_v5_genome.tsv.gz | 687M | (含上合計) | V5 HEAD 全基因組 vote dump | CURRENT | CONCLUDED_POSITIVE | **NEEDS_USER_DECISION→ARCHIVE** | (合計內) | 同上 | 同上 |
| research/v5_provenance_followup/T1_2_read_level_audit/vote_dump_v3f_genome.tsv.gz | 687M | (含上合計) | V3F 全基因組 vote dump | SUPERSEDED(by V5) | CONCLUDED_POSITIVE | **NEEDS_USER_DECISION→ARCHIVE** | (合計內) | 同上 | 同上 |
| research/v5_provenance_followup/T1_2_read_level_audit/vote_dump_*_chr19.tsv.gz (3 檔) | 40M | 41,010,669 | chr19 子集 vote dump（全基因組版的子集）| PRE-FIX/CURRENT | CONCLUDED_POSITIVE | NEEDS_USER_DECISION | 41,010,669 | 同上 | genome 版已覆蓋；chr19 子集為早期 pilot 殘留，可隨 genome dump 一同封存 |
| research/paired_priority_bug_audit/phaseC_genome_three_way/{baseline,V3F,V5,V6}_{on,off}_{tp,fp}/ | timeout(large) | -1 | baseline/V3F/V5/V6 × on/off × tp/fp 全基因組 per-region ISM 輸出（chr1 單一 config 即 2624 region 子目錄）| 混 PRE-FIX+CURRENT | CONCLUDED_POSITIVE | **NEEDS_USER_DECISION→ARCHIVE** | -1 (估 tens of GB；本區最大回收標的) | docs/.../20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md 等引用「summary」非「per-region」| 重跑可重生的 ISM pipeline 輸出；結論已落 validated .md（08_phaseD findings / V6 characterization）；16 config × 3 phaseC 變體 × 數千 region 大量複製 |
| research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance/ | timeout(large) | -1 | 同上加顯著性版本（部分重疊）| 同上 | CONCLUDED_POSITIVE | NEEDS_USER_DECISION | -1 | 同上 | 與無 significance 版高度重疊；可能其一 superseded |
| research/paired_priority_bug_audit/phaseC_v6_4sample_with_significance/ | timeout(large) | -1 | V6 4 樣本顯著性 per-region 輸出 | CURRENT | CONCLUDED_POSITIVE | NEEDS_USER_DECISION | -1 | V6 characterization 報告 | 重生可得；保守待用戶定 |
| research/paired_priority_bug_audit/{phaseD_v6_5sample, v6c_phaseB_runs, v6c_phaseB_v5bam_runs, v6c_phaseB_v6bam_runs}/ | timeout(large) | -1 | phaseB/D 各版 per-region 跑批中間輸出 | 混 | CONCLUDED_POSITIVE | NEEDS_USER_DECISION | -1 | 04/05/08 findings .md | 中間跑批；結論已成文 |
| research/v6_bam_tpfp_hp_loh_cn/step4_cross_sample_extension/intermediate/file_lists/*_on/off_tp/fp.txt (16 檔) | 83M | 86,891,610 | per-sample read/region 路徑清單（建 master 用的 file enumeration）| CURRENT | LIVE | **ARCHIVE** | 86,891,610 | step4_findings.md | 純路徑列舉中間檔，可由 pipeline 重生；on/off 成對重複；母專案 LIVE 但此項為 disposable intermediate |
| research/v6_bam_tpfp_hp_loh_cn/step4_cross_sample_extension/per_sample_master/*_v6_master.tsv (4 檔) | 40M | ~41,000,000 | H2009/H1437/HCC1954/HCC1937 v6 master TSV（cross-sample 分析核心輸入）| CURRENT | LIVE | **KEEP** | 0 | step4_findings.md（核心數據）| 跨樣本分析直接消費；可重生但為當前分析骨幹，保留 |
| research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_master_augmented.tsv | 26M | ~27,000,000 | methyl filter pilot augmented master | CURRENT | CONCLUDED（methyl filter DEAD ⭐2 L4）| NEEDS_USER_DECISION | ~27,000,000 | step5_findings.md | methyl-as-FP-filter 方向已 DEAD，但 Phase2 LOSO NEGATIVE 為 G6 論文 methods 證據 → 至少 ARCHIVE，**勿 SAFE_DELETE**；母專案 LIVE 故保守 |
| research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/intermediate/per_region_hp_counts.tsv.gz | 3.2M | ~3,300,000 | per-region HP counts 中間檔 | CURRENT | LIVE | KEEP | 0 | step1 deltas | 小；step1 分析直接用 |
| research/tsg_promoter_asm_reviewer/data/hg38.ncbiRefSeq.gtf.gz | 40M | 40,000,000* | NCBI RefSeq hg38 基因註解（公開下載）| NA(reference) | NA | **ARCHIVE** | 41,902,342 (含 cpg) | scripts/17_genome_asm_analysis.py | 公開 reference annotation，可重新下載（非分析輸出）；佔 tsg 專案 ~60%；可移至共用 reference 或外部 |
| research/tsg_promoter_asm_reviewer/data/cpgIslandExt_hg38.txt.gz | 702K | 702,736 | UCSC CpG island annotation（公開下載）| NA(reference) | NA | ARCHIVE | (含上合計) | scripts | 公開可重下；同上 |
| research/tsg_promoter_asm_reviewer/msa.log | 470K | 480,969 | MSA build 過程 log（transient）| NA | TRANSIENT | **SAFE_DELETE** | 480,969 | none | 純建構 log，重跑可重生，無證據價值 |
| research/tsg_promoter_asm_reviewer/genome_survey/ (v1) | 3.0M | 2,992,736 | 第一版全基因組 ASM survey（疑被 genome_survey_v2 取代）| SUPERSEDED? | CONCLUDED | NEEDS_USER_DECISION | 2,992,736 | scripts/17_genome_asm_analysis.py（仍引用）| memory 註記「validated v2 TSV」「buggy −0.054 是 MSA Level1 artifact」→ v1 疑為 buggy 前版；但 script 仍引用，需用戶確認 v1 是否 superseded |
| research/loh_investigation/figures/{loh_round2,loh_round3,loh_round4,phase1a_round3}/ | 2.6M | 2,557,755 | LOH 圖表 round2-4 + phase1a（多輪迭代版本）| 部分 SUPERSEDED | LIVE | NEEDS_USER_DECISION | 2,557,755 | none（scoped grep 未見直接引用）| 多輪 figure 迭代；早期 round 疑被後續取代，但母專案 LIVE 论文主軸 → 保守；建議用戶確認哪輪為權威 |
| research/tpfp_loh_af_kde_discrimination/08_archive_to_crosssample.md | 8.0K | 8,000 | 08 報告（INDEX 標 superseded）| SUPERSEDED | CONCLUDED | KEEP | 0 | 00_INDEX（標 superseded by 07）| 8KB 微檔；保留為紀錄不值得刪 |

\* genome dump 個別 bytes 為估值；三檔精確合計 = 2,219,541,005 bytes（du -bc 實測）。tsg gtf+cpg 合計 = 42,605,078 bytes（實測）。

---

## 反常事項（anomalies）

1. **v5 vote_dump .gitignore 註解 + 2026-05-10 cleanup proposal 嚴重低估大小**：兩處皆寫 ~26–40 MB，但 genome dumps 實測 **2.1 GB**（chr19 時代估值未更新）。封存決策的尺寸假設過期。
2. **paired_priority_bug_audit 總量無法 du 完**（>13 min 仍跑，因 phaseC 單一 config 的 filtered_snv_tp/chr1 即 2624 region 子目錄）。實際為本區**最大回收質量**（保守估 tens of GB 的 regeneratable ISM per-region 輸出），但 du 超時故 bytes=-1，未納入精確總計。
3. **phaseC 三變體高度重疊**：`phaseC_genome_three_way` vs `..._with_significance` 為同基底加/不加顯著性，疑其一 superseded — 需用戶確認以避免雙倍佔用。
4. **tsg genome_survey v1 vs v2 supersession 無定論**：v2 命名 + memory「validated v2 TSV / buggy 前版」暗示 v1 superseded，但 script 17 仍引用 v1 → 衝突，標 NEEDS_USER_DECISION。
5. **step5 methyl filter DEAD 但屬論文 methods 證據**：母專案 v6 LIVE 內含已 DEAD 的 methyl-as-FP-filter pilot；其 NEGATIVE 結論是 G6 論文 methods 段落證據 → 即使 DEAD 也只能 ARCHIVE，禁 SAFE_DELETE。

## 總計（不含 paired_priority_bug_audit timeout 項）

- total_scanned（15 個可精確 du 專案）= **2,644,245,001 bytes**（~2.64 GB）。paired_priority_bug_audit 另計 timeout(large)，估 tens of GB。
- KEEP（母專案 + 核心數據）≈ 上述大部分。
- ARCHIVE（明確 regeneratable/reference/重複）= v6 file_lists 86,891,610 + tsg gtf+cpg 42,605,078 = **129,496,688 bytes**。
- SAFE_DELETE（純 transient）= tsg msa.log **480,969 bytes**。
- NEEDS_USER_DECISION（保守，含 v5 dumps 2.26G + step5 27M + tsg survey v1 3.0M + loh figures 2.6M + chr19 dumps 41M + paired phaseC timeout）= **2,293,604,165 bytes 已知 + paired timeout(large 未計)**。
