<!--
建立時間: 2026-01-12 00:00
更新時間: 2026-03-25 11:20
狀態: validated
資料來源:
  - docs/standards/20260228_文件命名與狀態管理規範_01.md
  - docs/standards/20260228_output軟連結與版本控管規範_01.md
  - scripts/analysis/check_ai_agent_readiness.sh
-->

# 當前目標

## 1. 目前狀態

1. docs 重整已完成核心落地：
   - 命名：`YYYYMMDD_主題_流水號.md`
   - 報告分層：`reports/validated|finalized`
   - 實驗分層：`experiments/in_progress|validated|finalized`
2. `output/` 入口已固定為軟連結：
   - `output -> /big7_disk/liaoyoyo2001/big7_disk_output`
3. Knowledge MCP 已接入：
   - `.mcp.json` 指向 `/big8_disk/liaoyoyo2001/knowledge/scripts/mcp/knowledge_server.py`

## 2. AI Agent 主要入口

1. docs 導航：`docs/README.md`
2. 研究歷史索引：`docs/experiments/INDEX.md`（已試驗方向、成功/失敗結論、建議後續）
3. Agent 手冊：`docs/references/manual/20260301_AI_Agent_快速操作手冊_01.md`
4. 健康檢查：`scripts/analysis/check_ai_agent_readiness.sh`
5. 文件規範：`docs/standards/README.md`
6. 研究文件主 Agent：[/big8_disk/liaoyoyo2001/InterSubMod/.claude/agents/intersubmod-weekly-research-agent.md](/big8_disk/liaoyoyo2001/InterSubMod/.claude/agents/intersubmod-weekly-research-agent.md)
7. 研究文件使用手冊：[/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究報告Agent與Skills使用手冊_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究報告Agent與Skills使用手冊_01.md)

## 3. 當前進行中

### 主線切換摘要（2026-03-24 確認）

1. 方法學審查已完成，`label-first / cluster-first` 的合理性、TO rescue 邊界、`paired_persistent_final_fp` 機制都已有 closeout。
2. 當前最穩定的高層流程分工已定調為：
   - 第一層：`caller-first`
   - 第二層：`methylation-support`
   - 第三層：`annotation / QC / artifact triage`
3. 下一階段突破方向優先序已正式收斂為：
   - `Phase 1`：方向 E，`ML read classification`
   - `Phase 2`：方向 A+D，`normal methylation reference + CN/Purity-aware`
   - `Phase 3`：方向 B，`gene-level / mechanism-level evidence integration`
   - `Phase 4`：方向 C，`CpG 功能分層與智慧選點`
4. 樣本角色同步調整為：
   - `HCC1395 5kHz paired/TO`：discovery 與 stress-test
   - `HCC1395_DORADO paired/TO`：同 biological sample 的 cross-platform validation
   - 其他 pure paired 樣本：穩定性與可攜性檢查
5. 近期不再作為主線重複驗證的項目：
   - 單純重做 `paired pure` 合理性驗證
   - `H006 window_bp=1000`
   - `AlleleDelta standalone filter`
   - 以 gate 微調直接解決 germline ASM FP

### Phase 1 目前啟動點（2026-03-25）

1. `Phase 1` 已正式進入研究啟動，不再停留在 roadmap 標題層。
2. 目前最優先任務不是直接訓練模型，而是先建立 `per-read training/export layer`。
3. 已確認既有資產已足夠支撐第一版資料層：
   - `rescue_joined_features.tsv`
   - `reads.tsv`
   - `methylation.csv`
   - `distance/.../matrix.csv`
   - 既有 `TO support` 與 `paired old output` diagnostics
4. 目前最主要缺口：
   - 統一的 exporter
   - 清楚的 label schema
   - `5kHz / DORADO` harmonization 與 split 規則
5. 本輪啟動文件：
   - [20260325_Phase1_ML_read_classification研究啟動與資料缺口盤點_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260325_Phase1_ML_read_classification研究啟動與資料缺口盤點_01.md)
   - [20260325_Phase1_label_schema與harmonization規格_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260325_Phase1_label_schema與harmonization規格_01.md)
6. 本輪已落地的 Phase 1 腳本與 pilot：
   - 腳本：[export_phase1_read_training_table.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_phase1_read_training_table.py)
   - pilot 輸出：[/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_read_training_exporter_pilot](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_read_training_exporter_pilot)
   - 目前最小驗證結果：`HCC1395 5kHz TO`、`2 regions`、`77 read rows`、`0 missing regions`
7. 本輪新增的 Phase 1 baseline manifest：
   - 腳本：[build_phase1_training_manifest.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_phase1_training_manifest.py)
   - 完整 manifest：[/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_v1)
   - 四象限 smoke：[/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_smoke4](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_smoke4)
   - 目前結果：
     - 四個 baseline dataset 共 `141,014` 個 regions 已納入 manifest
     - smoke resolve `8/8` 成功、`missing=0`
     - smoke read export 共 `1,701` 筆 `read × region` rows
8. 本輪新增的 manifest-driven shard export：
   - 腳本：[export_phase1_manifest_shard.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_phase1_manifest_shard.py)
   - smoke4：[/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_smoke4](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_smoke4)
   - sample80：[/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample80](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample80)
   - 目前結果：
     - smoke4：`8/8` resolve、`1,701` read rows
     - sample80：`80/80` resolve、`14,165` read rows
9. 本輪已收斂的任務定義：
   - `Phase 1A` 已正式定義為 `within-tumor alt-support read classification`
   - exporter schema 已內建：
     - `dataset_role`
     - `harmonization_group`
     - `phase1a_task`
     - `phase1a_region_label`
     - `phase1a_read_label`
     - `phase1b_ready`
     - `phase1b_blocker`
   - `Phase 1B` 目前不可直接開始：
     - sample80 的 `14,165` 筆 reads 中 `is_tumor` 全為 `1`
     - 這代表現在的 read export 仍是 tumor-only universe
     - code-level 檢查已定位到：
       - [RegionProcessor.cpp](/big7_disk/liaoyoyo2001/InterSubMod/src/core/RegionProcessor.cpp) 會初始化 `normal_reader`
       - 但 `process_single_region(...)` 目前只接收 `tumor_reader`
       - `normal_reader` 在目前程式內沒有後續使用點
10. `Phase 1A split manifest` 已完成：
   - 腳本：[build_phase1a_split_manifest.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_phase1a_split_manifest.py)
   - 輸出：[/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_split_manifest_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_split_manifest_v1)
   - 目前結果：
     - discovery regions：`70,463`
     - external validation regions：`70,551`
11. `Phase 1A` 的第一版 benchmark 已完成：
   - baseline table：
     [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_head_to_head_baseline_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_head_to_head_baseline_v1)
   - round 1 文件：
     [20260325_Phase1A_read_classifier_benchmark_round1_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260325_Phase1A_read_classifier_benchmark_round1_01.md)
   - sample80 benchmark：
     [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample80_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample80_v1)
   - sample200 benchmark：
     [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample200_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample200_v1)
   - 關鍵結論：
     - `context-rich` 模型已可穩定完成 `Phase 1A`
     - `pure methyl baseline` 明顯不足
     - `paired-pure` 顯著較容易，`to-pure` 是主要困難子任務
     - `to-pure` 初版 failure diagnosis 已完成：
       [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_to_pure_failure_diagnosis_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_to_pure_failure_diagnosis_v1)
12. `Phase 1A` 的更嚴格 incremental test 與 multi-bio validation 已完成：
   - round 2 文件：
     [20260325_Phase1A_incremental_test與multi_bio_validation_round2_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260325_Phase1A_incremental_test與multi_bio_validation_round2_01.md)
   - mode-mixed sample400：
     [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample400](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample400)
     [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample400_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample400_v1)
     [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_sample400_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_sample400_v1)
   - paired-only multi-bio sample637：
     [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_paired_multibio_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_paired_multibio_v1)
     [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_paired_multibio_sample637](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_paired_multibio_sample637)
     [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_paired_multibio_sample637_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_paired_multibio_sample637_v1)
     [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_paired_multibio_sample637_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_paired_multibio_sample637_v1)
   - 目前最重要結果：
     - `methyl + context` 不是全域有效
     - `sample400` mode-mixed 測試中：
       - discovery holdout `delta F1=-0.0088`
       - external validation `delta F1=-0.0206`
       - 負增益主要由 `to-pure` 驅動
     - `paired-only multi-bio external validation` 中：
       - `context-only F1=0.8611`
       - `methyl+context F1=0.8722`
       - `delta F1=+0.0112`
       - bootstrap `95% CI=[+0.0044, +0.0188]`
     - sample-level 目前已明確支持增益的 paired samples：
       - `HCC1954`
       - `COLO829`
       - `HCC1395_DORADO paired`
     - `H1437 / HCC1937 / H2009` 仍存在異質性，不能直接宣稱所有樣本都提升
     - bucket-level diagnostics 顯示：
       - `sample400 to-pure` 的退步主要集中在 `TP REF` 誤判與 `FP Subclone / Strong`
       - `H2009` 的負向表現主要集中在 `FP + Weak`
       - `HCC1954 / COLO829` 的正增益來自明確的 `Subclone / Weak / Strong` buckets
13. 目前已確認的下一批擴充任務：
   - 將 `paired-pure multi-bio` 升級成 `Phase 1A` 的主要驗證軸
   - `to-pure` 與 `paired-pure` 分開建模，不再共用同一個 `methyl 是否有效` 結論
   - 針對 `H2009 / HCC1937 / H1437` 做 dataset-specific diagnostics
   - 補 `5kHz / DORADO / PAO / Google` 的 group-aware normalization 策略
   - 規劃 normal-read export layer，作為 `Phase 1B` 前置
   - 視需要回收 `DORADO candidate-specific` 路徑 provenance，但這已不是 blocking issue

> `2026-03-24` 新確認（方法學審查 + 研究收尾 + 全域分析報告）：
> - **研究方法與相關文獻突破方向全域分析報告完成**
>   - 整合 12 篇外部核心研究 + ISM 全實驗歷程，確認突破方向優先序：E（ML read classification）→ A+D（Normal ref + CN/Purity-aware）→ B（Gene-level integration）→ C（CpG 功能分層）
>   - 確認 ISM 獨特定位：「以 somatic SNV 為錨點 + read-level methylation clustering + 雙向顯著性驗證」尚無完整對標
>   - 詳見：[20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md)
> - **方法學審查全部完成**（16份觀察文件，Steps 1-8 + 補充分析）
>   - 詳見：[20260324_方法學審查全域結論報告_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/20260324_方法學審查全域結論報告_01.md)
> - **H006 window_bp=1000 TO 調整分析完成（理論 Rejected）**
>   - window=1000 → NumCpGs 降至 ~15（5x 減少）；ISM FN 覆蓋率從 7% 降至 ~1.4%
>   - PassedGating TP/FP 比值在低 CpG 條件下 = 0.98（無改善）；Tagged BAM 不可用
>   - 詳見：[20260324_H006_窗口大小TO調整分析_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/20260324_H006_窗口大小TO調整分析_01.md)
> - **paired_persistent_final_fp 深度追蹤完成（研究方向關閉）**
>   - 45個跨平台共享 FP：suggest_filter=False（所有 rule sweep 零觸發）
>   - 機制：(A) Strong/Weak FP(43%)=germline ASM，ISM 正確偵測但無法區分胚系/體細胞；(B) Noise FP(49%)=caller FP，無特徵可過濾
>   - 改善需要 normal sample 甲基化參考，ISM 架構下不可改善
>   - 詳見：[20260324_paired_persistent_final_fp_深度追蹤_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260324_paired_persistent_final_fp_深度追蹤_01.md)
> - **binary_methyl_high/low 死碼修正完成**（RegionProcessor.cpp/hpp，F1 delta=0，107/107 tests pass）

> `2026-03-23` 新確認：
> - 已完成 residual FP deep-dive workspace：[20260323_to_residual_fp_deep_dive](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260323_to_residual_fp_deep_dive)
> - `raw_absent` 可以再細分，但主切法不是 methylation hard filter，而是 `shared exact hotspot` vs `platform-specific tail`
> - `HCC1395` 與 `HCC1395_DORADO` 的 `raw_absent` 各自有 `11430 / 11424` 個，其中 `11220` 個為 cross-platform exact shared hotspot，兩邊 shared 比例都超過 `98.16%`
> - 針對整個 `raw_absent`、`shared raw_absent` 與 `persistent` 的同層級 rule sweep 皆無正向規則；只有 `HCC1395 raw_absent_platform_specific` 出現極小正訊號：`ad_alt<=3`，`11 FP / 7 TP`，`delta F1=+0.000032`
> - `paired_persistent_final_fp` 已完成深度追蹤（見 2026-03-24 條目），確認為 irreducible FP，研究方向關閉
> - 詳細報告：[20260323_TO_residual_FP_deep_dive_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260323_TO_residual_FP_deep_dive_01.md)

> `2026-03-22` 新確認：
> - 已完成 TO FP provenance 與 paired 對照 observation workspace：[20260322_to_fp_provenance_analysis](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260322_to_fp_provenance_analysis)
> - `HCC1395` 與 `HCC1395_DORADO` 的 TO raw FP universe 中，`PON / NonSomatic` caller veto 佔比都超過 `99.48%`，且本次檢查顯示 `caller_pon_filtered` 100% 帶顯式 `PoN_1..PoN_4` tag
> - `LongPhase-TO` 在這兩個 big7 TO pilot 上對 truth-scoped raw FP 的額外移除為 `0`；目前 TO final rule 只額外移除 `HCC1395 8 FP / 4 TP`、`DORADO 6 FP / 15 TP`
> - TO final residual FP 中，`paired` 可解比例仍超過 `99%`，主因幾乎全是 `paired_raw_absent`，不是 `normal_alt_support`；真正 hard blind spots 收斂為 `HCC1395 87` 與 `DORADO 77` 個 `paired_persistent_final_fp`
> - 詳細報告：[20260322_TO_FP來源分解與paired對照分析_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260322_TO_FP來源分解與paired對照分析_01.md)
> - 摘要：[20260322_TO_FP來源分解摘要_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260322_TO_FP來源分解摘要_01.md)


1. `2026-03-07 ~ 2026-03-24` 的純樣本與 TO 基線已完成，可作為後續新方法的固定對照：
   - 第一主 paired benchmark：`HCC1395 5kHz`
   - 第一 cross-platform benchmark：`HCC1395 ONT_Dorado + HCC1395BL`
   - TO discovery / validation：`HCC1395 5kHz TO` 與 `HCC1395_DORADO TO`
2. 這輪基線研究已明確回答的核心問題：
   - 甲基訊號目前無法取代 caller，較適合作為 `support / annotation / triage`
   - `label-first` 與 `cluster-first` 可保留為解釋框架，但不是下一階段唯一主研究目標
   - 真正需要突破的是 `read-level modeling`、`normal baseline`、`CN/purity-aware correction` 與更高層證據整合
3. `2026-03-07` 新版完整欄位 paired pure 正式 rerun 已完成：
   - `HCC1395 5kHz`：`LongPhase-S F1=0.8522` → `InterSubMod F1=0.8532`（`+0.0010`）
   - `HCC1395_DORADO`：`LongPhase-S F1=0.8592` → `InterSubMod F1=0.8590`（`-0.0002`）
   - 跨 7 個 pure paired 樣本比較後，`低 VAF + 偏高 AlleleDelta + 不穩定全域結構` 只在 `HCC1395 5kHz` 顯示明顯正增益
   - `Noise->Strong` 也幾乎只在 `HCC1395` 系列明顯出現，尤其 `5kHz` 最值得警戒
   - 新增 refined label 分析後，`Strong + 低 VAF + 高 AlleleDelta` 在 `HCC1395 5kHz` 可移除 `68` 個 FP、僅誤傷 `1` 個 TP（`F1 +0.000806`），但 `DORADO` 與其他樣本不支持將其升級為全域規則
   - `Noise->Strong` 本身不適合直接刪除，因為在 `5kHz` 會移除 `1659` 個 TP、僅移除 `51` 個 FP
   - 反向救援分析顯示：在 `HCC1395 5kHz`，`LongPhase-S` 丟掉的 caller TP 中，有 `45` 筆可用單純 `GQ>=20` 救回且不重新引入精確對齊集合中的 FP（`F1 +0.000739`）
   - 目前的 InterSubMod `VerificationClass / PairwiseDist` 對 `LongPhase-S` 的 TP rescue 沒有超過 caller `GQ` 的表現，甲基訊號現階段更適合用於 FP triage 與 Strong 細分類
4. 研究流程與紀錄制度化：
   - 研究藍圖：`docs/research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md`
   - 執行計畫：`docs/plans/2026/03/20260307_純樣本甲基研究執行計畫_01.md`
   - 研究手冊：`docs/references/manual/20260307_研究推進與實驗觀察手冊_01.md`
   - 最新進度報告：`docs/experiments/in_progress/2026/03/20260307_純樣本round執行與進度報告_01.md`
   - 候選規則與 shift 跨樣本分析：`docs/experiments/in_progress/2026/03/20260307_低VAF高AlleleDelta與shift跨樣本分析_01.md`
   - Strong 細分類與 haplotagged samtools 驗證：`docs/experiments/in_progress/2026/03/20260307_Strong細分類與samtools驗證分析_01.md`
   - LongPhase 救回 TP 與甲基救援分析：`docs/experiments/in_progress/2026/03/20260308_LongPhase救回TP與甲基救援分析_01.md`
   - LongPhase-TO 空間需求與中間產物確認：`docs/experiments/in_progress/2026/03/20260307_LongPhaseTO_空間需求與中間產物確認_01.md`
   - HCC1395 5kHz TO pilot 啟動與執行紀錄：`docs/experiments/in_progress/2026/03/20260307_HCC1395_5kHz_TO_pilot啟動與執行紀錄_01.md`
   - HCC1395 5kHz TO borderline rescue 特徵證據鏈整理：`docs/reports/validated/2026/03/20260308_HCC1395_5kHz_TO_borderline_rescue特徵證據鏈整理_01.md`
5. `2026-03-07` 新確認：
   - 真正要做 `LongPhase-TO + InterSubMod` 驗證時，仍需 `haplotagged BAM`；目前 repo pipeline `02_intersubmod.sh` 直接吃 `--tagged-bam`
   - `benchmark_split_snv_vcf.sh` 已補上 plain `.vcf` 自動 `bgzip + tabix`；`run_longphase_to_intersubmod_pilot.sh` 已補上 `--tag-prefix`
   - `20260307_hcc1395_to_pilot_1` 已完成：
     - `ClairS-TO`
     - `LongPhase-TO phase`
     - `LongPhase-TO haplotag`
     - `InterSubMod`
     - benchmark / design validation / round summary
   - `HCC1395 5kHz tumor-only` 本輪結果：
     - `ClairS-TO F1=0.7127`
     - `LongPhase-TO F1=0.7127`
     - `InterSubMod F1=0.7130`
     - `InterSubMod` 移除 `36` 個 FP、誤傷 `3` 個 TP，`F1 delta=+0.0003`
   - 在這個 `5kHz TO` pilot 上，`LongPhase-TO phase` 本身沒有改變 benchmark call set；目前正增益來自後續甲基過濾層
   - TO 場景下 `label-first` 仍非空：
     - `label_upgrade=2756`
     - `conflict=3227`
     - `consistent_strong=1471`
   - `5kHz TO tagged BAM` 實測大小為 `260G`，`.bai` 約 `113M`
   - 截至本輪主流程完成，磁碟大約剩：
     - `/big8_disk` `228G`
     - `/bip8_disk` `207G`
   - 因此目前**不需要立刻手動清理**，但在進入下一個 `200G+` 級 tagged BAM 任務前，應先將 `5kHz TO tagged BAM` 搬到 `Archive_pending_delete/`
   - `HCC1395_DORADO` 已補齊 TO baseline benchmark（`F1=0.7226`），仍是下一個最合理的第二階段驗證樣本
   - `2026-03-08` TO 特徵回查已確認：
     - `label_upgrade=2756` 中 `TP=2526`、`FP=230`，`FP rate=8.35%`
     - `conflict=3227` 中 `TP=2465`、`FP=762`，`FP rate=23.61%`
     - `Weak->Strong` 與 `Noise->Strong` 在 TO 下都以 `TP` 為主，不適合直接當 artifact 刪除規則
     - TO 下真正仍有效的 `FP` 特徵是 `low VAF + high AlleleDelta`，而且這批幾乎不落在 `label_upgrade/conflict`
     - 該特徵在 TO 下主要落在 `cluster_only + Strong->Weak + 低 GQ/QUAL`
     - `low VAF + high AlleleDelta + low CramersV` 在 `5kHz TO` 可移除 `36 FP`、僅誤傷 `2 TP`，`F1 +0.000261`
   - `2026-03-08` TO diagnostics 與 scheme 調整新確認：
     - 本輪 `5kHz TO` 是 `ClairS-TO pileup` 路線，不是 full model；目前結論先只對 pileup 成立
     - top `Strong->Weak + lowVAF/highAD` 可疑 FP 中有 `4/5` 聚集在 `chr11` 的 `1402 bp` 區段
     - `HP skew` 本身不足以區分 FP/TP；TP 對照也可呈現極端單一 HP
     - 比較有辨識力的是 `高 AlleleDelta + 低 PairwiseMedianDist + 很低 alt fraction + 低 QUAL/GQ`
     - 現行 core 的 `Strong` 會混入「只有 allele sig、沒有 HP sig」的位點；這批不適合整包刪除（`5104 TP / 2042 FP`）
     - 新增 `cluster_plus_weak_label` 後，原先的 TO 可疑子集更精確地落在：
       - `cluster_plus_weak_label + Strong->Weak + lowVAF/highAD`
     - 分類語意調整目前主要提升的是可解釋性，尚未超過既有 `low VAF + high AlleleDelta + low CramersV` 的 F1 增益
   - `2026-03-08` ClairS 邊緣 TP rescue 與甲基輔助評估第一輪已完成（**只用 final VCF，不用 pileup candidate/tensor**）：
     - `HCC1395 5kHz paired`：
       - `caller_lost_tp=1052`（`candidate_eligible=920`）
       - `caller_removed_fp=12974`（`candidate_eligible=5266`）
       - 依規劃排序的最佳 safe rule：`candidate_gq_ge_15`，`106 TP / 75 FP`，`F1 +0.000825`
       - 若只看最佳 `delta F1`：`gq>=20`，`59 TP / 8 FP`，`F1 +0.000871`
       - 甲基 summary 目前只覆蓋 `111/1052` lost TP 與 `802/12974` removed FP；本輪 `methylation support/veto` 尚未超過 caller-only
     - `HCC1395_DORADO paired`：
       - `caller_lost_tp=1489`（`candidate_eligible=1122`）
       - `caller_removed_fp=2533`（`candidate_eligible=1658`）
       - 依規劃排序的最佳 safe rule：`candidate_gq_ge_15`，`97 TP / 88 FP`，`F1 +0.000502`
       - 若只看最佳 `delta F1`：`gq>=20`，`53 TP / 7 FP`，`F1 +0.000782`
       - 尚缺 `DORADO paired` 的 candidate-specific InterSubMod summary，因此本輪不能把甲基 rescue 解讀成有效 negative result
     - `HCC1395 5kHz TO`：
       - `caller_lost_tp=2108`（`candidate_eligible=773`）
       - `caller_removed_fp=2720008`（`candidate_eligible=298`，`pon_hit_fraction=0.999019`）
       - 最佳 caller-only safe rule：`qual>=10 or gq>=20`，`491 TP / 118 FP`，`F1 +0.006824`
       - 第二名：`gq>=15`，`365 TP / 71 FP`，`F1 +0.005233`
       - candidate-specific InterSubMod 已補跑完成：
         - `candidate_lost_tp export=773`、`candidate_removed_fp export=298`，missing 均為 `0`
         - 可分析 region：`lost_tp=675/773 (87.32%)`、`removed_fp=213/298 (71.48%)`
       - 最佳甲基 support：
         - `PairwiseMedianDist >= 0.20`：`300 TP / 68 FP`，`F1 +0.004219`
       - 最佳 label-aware support：
         - `agreement_positive`：`148 TP / 25 FP`，`F1 +0.002163`
         - 比單純 `Strong/Subclone`（`149 TP / 30 FP`）更乾淨
       - `artifact veto` 對 TO rescue 幫助有限：
         - `lowVAF/highAlleleDelta/lowCramersV` 與 `gq>=15` 組合沒有額外效果
         - `combined_artifact_veto` 只少 `6 FP`，但同時少 `57 TP`
     - **結論**：TO 下甲基資料可以幫 borderline TP rescue，但目前仍是次要 support 訊號，尚未超過最佳 caller-only 規則
     - `2026-03-09` 已補齊 `HCC1395_DORADO paired` candidate-specific InterSubMod：
       - scratch 改寫到 `/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/`
       - 重建 `HCC1395_DORADO_tagged.bam` 成功，大小 `223G`；`bam.bai` 已完成
       - `caller_lost_tp` 可分析 `93/1122 (8.29%)`
       - `caller_removed_fp` 可分析 `326/1658 (19.66%)`
	     - `gq>=10 + PairwiseMedianDist>=0.20`：`10 TP / 36 FP`，`F1 -0.000280`
	     - `gq>=10 + agreement_positive`：`36 TP / 72 FP`，`F1 -0.000298`
	     - `gq>=10 + Strong/Subclone`：`25 TP / 38 FP`，`F1 -0.000059`
	     - `caller-only` 仍最佳：`gq>=15` 為 `97 TP / 88 FP`，`F1 +0.000502`；`gq>=20` 為最佳 `delta F1`（`+0.000782`）
	     - **結論**：`5kHz TO` 已成立的甲基 rescue support 規則，目前**不能升級成跨樣本規則**；`DORADO paired` 新證據反而支持 caller-first 穩定、methylation-support 不穩定
	     - `2026-03-09` feature-space 深入分析已完成：
	       - 新腳本：[analyze_methylation_rescue_feature_space.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_methylation_rescue_feature_space.py)
	       - 主報告：[20260309_5kHz_TO與DORADO_paired_甲基rescue特徵空間分析_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260309_5kHz_TO與DORADO_paired_甲基rescue特徵空間分析_01.md)
	       - `HCC1395 5kHz TO`：
	         - `Quality_Score>=60`、`PairwiseMedianDist>=0.15/0.20`、`agreement_positive` 都能從 baseline 救回 TP
	         - 但全部 **沒有超過 `gq>=10` caller gate**
	         - `allele_assign_rate>=0.99` 與 `CramersV<=0.05` 看似數值高，但屬技術 proxy，不應升級為正式甲基 rescue 特徵
	       - `HCC1395_DORADO paired`：
	         - support 方向與 `5kHz TO` 不同；正向訊號偏向 `高 hp_assign_rate + 低 Pairwise + 較高 Quality_Score`
	         - 單特徵最佳增量為 `gq>=15 + hp_assign_rate>=0.99`：`50 TP / 15 FP`，相對 gate `F1 +0.000132`
	         - 雙特徵最佳增量為 `gq>=15 + Pairwise<=0.20 + hp_assign>=0.99`：`46 TP / 12 FP`，相對 gate `F1 +0.000103`
	         - 但整體仍 **沒有超過 `gq>=20` caller-only**
	     - `2026-03-11` GQ 與甲基 rescue 系統性驗證已完成：
	       - 新腳本：[analyze_gq_methylation_rescue_matrix.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_gq_methylation_rescue_matrix.py)
	       - 主報告：[20260311_GQ與甲基rescue系統性驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_GQ與甲基rescue系統性驗證_01.md)
	       - 主輸出目錄：[/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix/](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix/)
	       - `GQ` 仍是最穩定的 caller-first rescue 訊號，但不能直接解讀成 `ClairS-TO` 應全域放寬；在 `HCC1395 5kHz TO` 中，`caller_removed_fp=2,720,008`，真正 `candidate_eligible` 只有 `298`
	       - `HCC1395 5kHz TO`：
	         - `gq>=10`：`499 TP / 119 FP`，`F1 +0.006943`
	         - `gq>=15`：`365 TP / 71 FP`，`F1 +0.005233`
	         - 甲基單獨有效但不超過 caller：`Quality_Score>=60` 為 `556 TP / 162 FP`，`F1 +0.007466`；`PairwiseMedianDist>=0.20` 為 `408 TP / 132 FP`，`F1 +0.005374`
	         - `gq>=10 + 甲基` 目前沒有超過 `gq>=10 only`；`gq>=10 + Pairwise>=0.20` 為 `300 TP / 68 FP`，相對 gate `delta F1=-0.002724`
	       - `HCC1395 5kHz paired`：
	         - `gq>=10` 為負增益（`183 TP / 616 FP`, `F1 -0.004459`）
	         - paired 最穩定仍是 `gq>=15` 或 `gq>=20`
	         - 甲基單獨與 `gq>=10 + 甲基` 目前都不形成正向 rescue
	       - `HCC1395_DORADO paired`：
	         - caller-first 仍以 `gq>=15 / gq>=20` 為主
	         - 甲基單獨最佳僅 `hp_assign_rate>=0.99`，`68 TP / 87 FP`，`F1 +0.000041`
	         - `gq>=10 + quality>=60 + hp_assign>=0.99` 雖比 `gq>=10` gate 更好（`delta F1 vs gate=+0.001961`），但整體仍低於 `gq>=20`
	       - `PairwiseMedianDist` 呈現明顯區間性與樣本依賴：
	         - `5kHz TO` 在 `[0.15,0.25)` 為正向 support，但 `>=0.25` 不再更好
	         - `DORADO paired` 方向相反，較低 pairwise 反而更合理
	       - `low VAF + high AlleleDelta (+ low CramersV)` 仍較適合後段 artifact triage，不適合提前升級成 TP rescue 主規則
	     - `2026-03-11` 已補齊 `HCC1395_DORADO TO` candidate-specific InterSubMod：
	       - 主報告：[20260311_HCC1395_DORADO_TO_candidate_specific甲基rescue驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_HCC1395_DORADO_TO_candidate_specific甲基rescue驗證_01.md)
	       - round 目錄：[/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue)
	       - 這輪採用 `candidate-window subset BAM + subset haplotag`，避免生成 `200G+` full tagged BAM；整個 round 體積約 `3.1G`
	       - candidate-specific coverage 已達 `100%`：
	         - `caller_lost_tp=247/247`
	         - `caller_removed_fp=94/94`
	       - `DORADO TO` 的 caller-only rescue 仍成立：
	         - `gq>=10`：`40 TP / 11 FP`，`F1 +0.000540`
	         - `gq>=15`：`33 TP / 9 FP`，`F1 +0.000446`
	       - 甲基單特徵不是負效果：
	         - `Quality_Score>=60`：`208 TP / 77 FP`，`F1 +0.002620`
	         - `PairwiseMedianDist<=0.20`：`173 TP / 60 FP`，`F1 +0.002217`
	         - `hp_assign_rate>=0.99`：`122 TP / 41 FP`，`F1 +0.001577`
	       - 但在固定 caller gate 後，甲基目前仍只作 support：
	         - `gq>=10 + Quality>=60`：相對 gate `delta F1=-0.000064`
	         - `gq>=10 + Pairwise<=0.20`：`-0.000142`
	         - `gq>=10 + agreement_positive`：`-0.000430`
	       - 這代表 `TO` 下跨 `5kHz` 與 `DORADO` 目前都支持同一個高層結論：
	         - 第一層：`caller-first`
	         - 第二層：`methylation-support`
	         - 獨立旁路：`artifact triage`
	       - `PairwiseMedianDist` 的方向仍具樣本依賴：
	         - `5kHz TO` 偏高 pairwise 為正向 support
	         - `DORADO TO` 偏低 pairwise（`<=0.20`）較合理
	         - 因此不可升級成所有 TO 樣本的全域規則
	     - `2026-03-11` 已補做 `TO support` 特徵的 read-level diagnostics：
	       - 主報告：[20260311_TO_support特徵_readlevel_diagnostics_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md)
	       - 執行腳本：[run_to_support_feature_diagnostics.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_to_support_feature_diagnostics.py)
	       - 主輸出目錄：[/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics)
	       - `5kHz TO` 與 `DORADO TO` 都已補上代表性 `TP/FP` 的 matrix heatmap 與 `samtools` snapshot
	       - `Quality_Score` 在兩個 TO 資料集中都不是負訊號，但 `FP` 也可同樣很高；目前最合理定位是 `soft support / ranking / annotation`
	       - `PairwiseMedianDist` 的方向依然顯示明顯 dataset-dependence：
	         - `5kHz TO` 偏高 pairwise 的代表性 TP 常帶有較高結構異質性
	         - `DORADO TO` 偏低 pairwise 的代表性 TP/FP 都存在，顯示其更適合作 dataset-aware annotation
	       - `hp_assign_rate` 的代表性 `TP/FP` 都常達 `1.0`，而且 `FP` 常伴隨極高 alt fraction 與極端 haplotype skew；目前應降階為 `phase/QC annotation`
	       - 更新後的最合理流程分工：
	         - 第一層：`caller-first`
	         - 第二層：`Quality_Score` 等甲基 support
	         - annotation/QC：`PairwiseMedianDist`、`hp_assign_rate`、`NA HP fraction`、`haplotype skew`
	     - `2026-03-11` 已補做四象限整合矩陣與可比性盤點：
	       - 主報告：[20260311_HCC1395_四象限甲基rescue整合報告_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md)
	       - 整合腳本：[build_hcc1395_cross_quadrant_matrix.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_hcc1395_cross_quadrant_matrix.py)
	       - 主輸出目錄：[/bip8_disk/liaoyoyo2001/InterSubMod_out/output/research_rounds/20260311_hcc1395_cross_quadrant_matrix](/bip8_disk/liaoyoyo2001/InterSubMod_out/output/research_rounds/20260311_hcc1395_cross_quadrant_matrix)
	       - 目前已把 `HCC1395 5kHz / HCC1395_DORADO × paired / TO` 的 benchmark、candidate-specific coverage、特徵矩陣、交集/正交性與 TO snapshot scope 差異放進同一個可比框架
	       - 已確認：
	         - paired raw caller 有 `pileup_filter.vcf` 與 `full_alignment_filter.vcf`，但目前只完成 availability 與 phase 2 規劃，尚未做四象限 direct benchmark
	         - TO raw caller 目前只有 pileup-only 路線；不能假設已有 full model 可直接比較
	         - `5kHz TO full tagged BAM` 與 `DORADO TO candidate-window subset tagged BAM` 只能做 read-level 方向比較，不能做絕對值硬比較
	       - 目前最合理高層分工仍維持：
	         - 第一層：`caller-first`
	         - 第二層：`methylation-support`
	         - 獨立旁路：`artifact triage`
	         - 補充註記：`annotation / QC`
	     - 本輪正式報告：
	       - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_ClairS邊緣TP_rescue與甲基輔助評估_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_ClairS邊緣TP_rescue與甲基輔助評估_01.md)
	       - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_HCC1395_5kHz_TO_candidate_specific甲基rescue驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_HCC1395_5kHz_TO_candidate_specific甲基rescue驗證_01.md)
	       - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260309_HCC1395_DORADO_paired_candidate_specific甲基rescue驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260309_HCC1395_DORADO_paired_candidate_specific甲基rescue驗證_01.md)
	       - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260309_5kHz_TO與DORADO_paired_甲基rescue特徵空間分析_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260309_5kHz_TO與DORADO_paired_甲基rescue特徵空間分析_01.md)
	       - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_GQ與甲基rescue系統性驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_GQ與甲基rescue系統性驗證_01.md)
	       - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_HCC1395_DORADO_TO_candidate_specific甲基rescue驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_HCC1395_DORADO_TO_candidate_specific甲基rescue驗證_01.md)
	       - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md)
	       - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_phase2_paired_raw與feature_interval_orthogonality分析_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_phase2_paired_raw與feature_interval_orthogonality分析_01.md)
	     - `2026-03-11` phase 2 已完成：
	       - 新腳本：[run_phase2_paired_model_feature_analysis.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase2_paired_model_feature_analysis.py)
	       - 主輸出目錄：[/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis)
	       - `paired raw pileup/full` 直接 benchmark 已確認：
	         - `5kHz paired` 中 raw `pileup/full` 都比 merged output 有更高 recall，但 FP 暴增，F1 分別掉到 `0.7896 / 0.7422`
	         - `DORADO paired` 中 raw `full` 優於 raw `pileup`，但仍低於 merged output（`0.8545 < 0.8565`）
	       - finer feature interval 顯示：
	         - `5kHz TO` 的 `PairwiseMedianDist` 正向 support 帶集中在 `0.18~0.20`
	         - `DORADO TO` 則集中在 `0.12~0.15`
	         - `Quality_Score` 在兩個 TO 資料集都偏向 `55~60` 區間較乾淨
	       - orthogonality 分析顯示：
	         - paired 目前沒有清楚的甲基正交補強，仍以 `GQ` 為主
	         - TO 中 `GQ / Pairwise / hp_assign / agreement / Strong` 間有小幅 union 增益，但仍屬 support，不足以取代 caller-first
	     - `2026-03-11` 已將 phase 2 finer interval 回接到 analysis-layer annotation：
	       - 新腳本：[build_phase2_annotation_layer.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_phase2_annotation_layer.py)
	       - 主輸出目錄：[/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer)
	       - 已補出四個 dataset 的 annotated candidates、annotation summary、policy evaluation 與 overlap summary
	       - `Quality_Score` 現階段最適合升級為 `support annotation`
	       - `PairwiseMedianDist` 現階段最適合升級為 `dataset-aware support annotation`
	       - `hp_assign_rate` 應明確降階為 `phase/QC annotation`
	       - annotation tier 實測後，目前沒有任何 policy 穩定超過最佳 caller-first gate：
	         - `HCC1395 5kHz TO`：最佳仍是 `gq>=10`
	         - `HCC1395 5kHz paired`：最佳仍是 `gq>=15`
	         - `HCC1395_DORADO TO`：最佳仍是 `gq>=10`
	         - `HCC1395_DORADO paired`：最佳反而是更嚴格的 `gq>=20`
	       - 高層結論正式收斂為：
	         - 第一層：`caller-first`
	         - 第二層：`methylation-support`
	         - 第三層：`annotation / QC / artifact triage`
       - 本輪正式報告：
         - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_phase2_finer_interval回接annotation_layer驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_phase2_finer_interval回接annotation_layer驗證_01.md)
         - [/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md)
   - `2026-03-08` ClairS 邊緣 FN 探勘（Pool B）與甲基救援分析已完成（**實驗 25**）：
     - `Pool B FN`（ClairS non-PASS ∩ truth set）：**840 個**
     - 策略 B 補跑後甲基覆蓋率：**99.9%（839/840）**
     - **最佳 Caller-only**：`no_varcluster`（F1 delta=+0.003391）；最乾淨：`no_varcluster_and_gq15`（precision=71.9%，F1 delta=+0.001452）
     - **最佳 Combined（Caller+Methyl）**：`gq15_and_allele_delta_low`（F1 delta=+0.002702），僅比 Caller-only 好 +0.000386
     - **甲基唯一規則全數為負**：Pool B non-FN 含大量 germline 變異，Strong rate 比 FN 更高（44.7% vs 30.8%）
     - **AlleleDelta 具 somatic 特異性**：FN 中位 AD=0.010，non-FN 中位 AD=0.065，但 precision 仍只 ~32%
     - **結論**：甲基訊號不適合單獨用於 Pool B rescue；現行標準無需調整，問題在 Pool B non-FN 本身混入 germline
     - **不推薦升為正式規則**；`no_varcluster_and_gq15` 建議在 `HCC1395_DORADO TO` 做交叉驗證
     - 正式報告：`docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md`
   - `2026-03-12` 未完項 closing plan 已建立並完成三條 closeout：
     - closing plan：[/big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260312_未完項closing_plan_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260312_未完項closing_plan_01.md)
     - round 0 盤點：[/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260312_TO雙模型與snapshot_scope盤點_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260312_TO雙模型與snapshot_scope盤點_01.md)
     - validated closeout：[/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260312_TO雙模型可得性closeout_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260312_TO雙模型可得性closeout_01.md)
     - same-scope control：[/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260312_TO_snapshot_scope_same_scope_control_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260312_TO_snapshot_scope_same_scope_control_01.md)
     - Pool B closeout：[/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260312_PoolB_FN_integration_closeout_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260312_PoolB_FN_integration_closeout_01.md)
     - 本輪 3 條 closing workstream 已全部完成：
       - `TO pileup vs full model`：已正式 closeout，現有 source tree 中尚未確認 full model 產物；目前 TO caller 證據應視為 `pileup-route only`
       - `5kHz TO` vs `DORADO TO` snapshot：`5kHz TO` same-scope control 已完成，證明 subset 技術本身不造成指標偏移；目前不可做跨 dataset 絕對值硬比較的主因是 dataset/platform 差異，不是 subset 技術
       - `Pool B FN`：已完成 integration closeout；再次支持 `caller-first`，甲基只提供弱補強，且 `AlleleDelta` 不可跨母體全域化
     - `2026-03-12` same-scope control 已定案：
       - `HCC1395 5kHz TO` 對相同 18 個代表性 region 補做 `candidate-window subset snapshot`
       - `full tagged BAM` vs `subset tagged BAM` 比較結果為 `18/18 fully identical`
       - `read_count`、`target_alt_fraction`、`na_hp_fraction`、`collapsed_hp_balance_delta` 等全部 `max_abs_delta = 0.0`
       - 結論：`subset snapshot` 可安全用於同 dataset 內的 read-level ranking 與 diagnostics；目前 `5kHz TO` 與 `DORADO TO` 仍不可做跨 dataset 絕對值硬比，主因是 dataset/platform 差異，而不是 subset 技術本身
     - `2026-03-12` Pool B integration closeout 已定案：
       - `Pool B FN = 840`
       - 原始 caller-only 最佳規則為 `no_varcluster`（`428 TP / 390 FP`, `F1 +0.003391`）
       - 最乾淨的 caller-side 規則為 `no_varcluster_and_gq15`（`115 TP / 45 FP`, `precision=71.9%`）
       - 最佳 combined 規則為 `gq15_and_allele_delta_low`（`431 TP / 473 FP`, `F1 +0.002702`），僅比相近 caller gate 多 `+0.000386`，仍低於原始 caller-only 最佳
       - `Pool B` 中 `AlleleDelta` 與 TO kept-set / artifact triage 的方向不同，不能全域化
6. docs 維護工作仍持續進行，但現已降為研究主線的次要維護事項：
   - 降低舊語彙殘留（例如舊輸出路徑與舊文件分層名稱）
   - 收斂腳本路徑設定（優先使用可配置 `OUTPUT_ROOT` 或 `output/` 入口）
   - 清理不必要空目錄，降低 Agent 探索噪音
7. `2026-03-10` 研究文件 Agent 與 Skills 第一版已落地：
   - 主 Agent：[/big8_disk/liaoyoyo2001/InterSubMod/.claude/agents/intersubmod-weekly-research-agent.md](/big8_disk/liaoyoyo2001/InterSubMod/.claude/agents/intersubmod-weekly-research-agent.md)
   - Skills：
     - [/home/liaoyoyo2001/.codex/skills/intersubmod-context-synthesizer/SKILL.md](/home/liaoyoyo2001/.codex/skills/intersubmod-context-synthesizer/SKILL.md)
     - [/home/liaoyoyo2001/.codex/skills/intersubmod-weekly-report-writer/SKILL.md](/home/liaoyoyo2001/.codex/skills/intersubmod-weekly-report-writer/SKILL.md)
     - [/home/liaoyoyo2001/.codex/skills/intersubmod-report-prompt-refiner/SKILL.md](/home/liaoyoyo2001/.codex/skills/intersubmod-report-prompt-refiner/SKILL.md)
   - 手冊與規格：
     - [/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_InterSubMod研究文件Agent規格_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_InterSubMod研究文件Agent規格_01.md)
     - [/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究脈絡整理Skill規格_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究脈絡整理Skill規格_01.md)
     - [/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_週報生成Skill規格_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_週報生成Skill規格_01.md)
     - [/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_指令修正與偏好收斂Skill規格_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_指令修正與偏好收斂Skill規格_01.md)
     - [/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_個人化研究寫作偏好設定規格_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_個人化研究寫作偏好設定規格_01.md)
     - [/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究報告Agent與Skills使用手冊_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究報告Agent與Skills使用手冊_01.md)
   - 目前狀態：可直接用來產生週報、主題式證據鏈與 prompt refinement；下一輪值得補上真實 eval prompts 與 changelog
8. `2026-03-10` 本週主線週報與週報指令手冊已完成：
   - 本週主線週報：[/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260310_研究主線週報_20260305_20260310_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260310_研究主線週報_20260305_20260310_01.md)
   - 週報指令手冊：[/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究週報撰寫指令與skill草案_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究週報撰寫指令與skill草案_01.md)
   - 可直接作為之後研究週報、主題式證據鏈與 prompt refinement 的第一版標準入口

## 4. 阻塞與風險

1. `archive/deep/` 為 immutable 快照，保留歷史失效連結（不回寫）。
2. 活躍腳本仍有部分硬編碼歷史路徑，可能影響跨專案重用。
3. `purity-aware` / subsample 結果容易被 tumor-normal 組織甲基差異混淆，暫不作主證據。
4. `label-first` 與 `cluster-first` 雖已有輸出，但方法學正當性與適用邊界仍需驗證。
5. 若只報 F1 而不附 `truth_total/TP/FP/FN`，容易高估方法改善幅度或誤判結論。
6. `pileup` 與 `full model` 的 TO 輸入若混用，會直接改變 caller VCF 與 downstream 結論；後續報告必須明確標示來源。
7. `TO` 的 borderline rescue 若直接拿 baseline kept-set 的 InterSubMod summary 來做 join，會得到錯誤結論；必須使用 candidate-specific `lost_tp / removed_fp` run。

## 5. 每次研究啟動必查

1. `docs/CURRENT_FOCUS.md`
2. `docs/experiments/INDEX.md`
3. `docs/README.md`
4. `docs/research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md`
5. `docs/references/manual/20260307_研究推進與實驗觀察手冊_01.md`
6. `/big8_disk/liaoyoyo2001/InterSubMod/.claude/agents/intersubmod-weekly-research-agent.md`
7. `/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究報告Agent與Skills使用手冊_01.md`
