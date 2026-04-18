<!--
建立時間: 2026-01-12 00:00
更新時間: 2026-04-17 22:30
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
   - 資訊分層：Active / Recent / Archive（詳見 AGENTS.md）
2. `output/` 入口已固定為軟連結：
   - `output -> /big7_disk/liaoyoyo2001/big7_disk_output`
3. Knowledge MCP 已接入：
   - `.mcp.json` 指向 `/big8_disk/liaoyoyo2001/knowledge/scripts/mcp/knowledge_server.py`

## 2. AI Agent 主要入口

1. docs 導航：`docs/README.md`
2. 研究歷史索引：`docs/experiments/INDEX.md`
3. 研究全景：`docs/reports/research_landscape/00_INDEX.md`
4. Agent 手冊：`docs/references/manual/20260301_AI_Agent_快速操作手冊_01.md`
5. 健康檢查：`scripts/analysis/check_ai_agent_readiness.sh`
6. 文件規範：`docs/standards/README.md`

## 3. 當前進行中

### Phase 2 研究方向（優先序）

1. **Phase 2 方向 A+D**：Normal Methylation Reference + CN/Purity-aware correction
   - **Phase A-D 程式碼已完成（2026-04-13）**：Normal BAM integration, Sample ASM, Somatic HP ASM, LOH BED annotation, Cross-region subclone analysis
   - HCC1395 全基因體驗證通過：Sample ASM 97.3% sig, Normal Baseline 100% valid, LOH concordance 94.1%, 4-group subclone
   - [驗證報告](experiments/validated/2026/04/20260413_Phase_BCD_Dual_BAM_Validation_01.md)
   - 待執行：7 samples 全量驗證、haplotag 重跑後的 Phase 2A 正式分析

2. **Phase 2 方向 B**：Gene-level / mechanism-level evidence integration
   - 目前狀態：待 Phase 2A 全量分析完成後啟動

3. **Phase 2 方向 C**：CpG 功能分層與智慧選點
   - 目前狀態：待規劃

### Self-Phasing 修正後重跑

- PON-only phasing 已驗證：LOH.bed 不變、somatic bias 消除、N50 +99.7%、phased rate +23.6pp
- **P0-3 LOH.bed 生成機制已確認（2026-04-11）**：LOH.bed 使用 VCF AF/VAF（`PhasingGraph.cpp:1817`，VAF >= 0.8 → HOM），ISM hp_ratio 使用 BAM HP tags（`ReadParser.cpp:123`）。兩套系統使用不同數據源，解釋 Jaccard=1.0 與 62% ISM LOH 消失的矛盾。
- 待執行：haplotag + ISM 全量重跑（7 samples × paired + TO）
- 重跑後可啟動：Phase 2A normal methylation reference baseline

### Phase 1A 已鎖定

- 最優方案：`methyl+context`，paired-pure multi-bio external validation delta F1=+0.0112, CI=[+0.0044,+0.0188]
- TO 模式下甲基化增益為負（delta F1=-0.0206）
- LOH feature 對 read-level 分類無用（LOH 是 region-level）
- 已正式關閉 TO single-sample post-hoc germline FP 識別

### LOH 雙定義與特徵探索全面關閉（2026-04-06 確定性結論）

- SEQC2 外部驗證 Jaccard=0.928：LOH.bed 與 FDA 金標準幾乎完全吻合
- LOH 區域 10/10 filter 策略全失敗（FP rate=0.239 < Non-LOH 0.338，LOH 是 TP-enriched）
- Non-LOH max AUC<0.58，多特徵 Voting AUC=0.577
- cnLOH 表面 0.587 是 Simpson's Paradox（per-sample mean=0.50）
- QS mode-aware 已實作：TO 模式停用 LOH penalty 與 verify bonus
- **結論：ISM 價值在 read-level epigenetic characterization，非 variant filter**

### Coverage_Multiple GC 校正與甲基化-CN 驗證 — 全 NO-GO（2026-04-11 確定性結論）

- **TO Pipeline 數據**：ClairS-TO → LongPhase-TO → ISM（TP=28,383, FP=11,830, F1=0.7127）
- **GC-Content 校正**：delta-r = -0.0002（Go/No-Go 門檻 ≥0.03）；ONT 5kHz GC bias 極小（98.7% regions 變化<5%）
- **甲基化-CN 相關**：所有 HP-free 特徵 residualized |r| < 0.07（甲基化對 CN 完全無感）
- **HPFineNGroups-CN**：raw r=0.495 → residualized r=0.160（68% 是 NumReads confound）
- **KDE auto-estimation 確認正確**：CN 分類準確度 6.2%→43.8%（已實作，`--expected-coverage` CLI 可用）
- **Coverage_Multiple 作為 CN proxy 已足夠**：TO r=0.827，Paired r=0.831
- 腳本：`scripts/analysis/gc_correction_to_validation.py`, `scripts/analysis/methylation_cn_validation.py`
- 圖表：Figures 28-33（`research/seqc2_cnv_stratification/figures/`）
- **結論：GC 校正不需實作；甲基化無法驗證 CN；Coverage_Multiple 現有精度已足夠**

### SEQC2 CNV 分層觀察 — CNV zone-aware filter 關閉（2026-04-10 確定性結論）

- **SEQC2 正交驗證 CNV truth set**（6 callers × 21 replicates × 3 technologies）用於 HCC1395 分層
- **Coverage_Multiple 驗證**：與 SEQC2 真實 CN Pearson r=0.831，ISM read-depth 代理可信
- **HCC1395 單樣本假象**：Gain+LOH zone AUC=0.782（AlleleDelta），但跨 7 樣本驗證 mean AUC ≤ 0.641
- **Simpson's Paradox 排除**：CNV 非 Quality_Score pooling 問題根源（分層後 QS diff=-0.042，pooling 反而有利）
- **Zone 排除策略全不可行**：所有排除策略 trade-off 低於 break-even（如排除 CN_Loss 移除 45% FP 但損失 11% TP）
- **Gain+LOH 根因**：CN=3+LOH 環境造成最高 FP rate 12.9%（allele imbalance 誤導 caller），FP rate 隨 CN 增加反而下降
- **15 張圖表 + 5 TSV 數據檔**：報告見 `docs/experiments/in_progress/2026/04/20260409_SEQC2_CNV分層觀察_01.md`
- **結論：CNV 不是特徵空間耗盡的根因；zone-aware filter 無可行策略，正式關閉**

### R1-R5 特徵設計研究（2026-04-07 確認）

- **R1**: CramersV 93% 為零 = 2×2 框架缺陷；HPMergedDelta 多群時反向；HPFineNGroups 已克服（AUC +0.125）
- **R2**: Excess groups 概念有效（跨子集統一 +0.059）但子集內無新信號，不需修改 C++
- **R3**: 結構清楚子集 AUC 反而下降 → **確認是 identifiability 問題而非特徵設計問題**
- **R4**: HPFineNGroups N≥4 + NR≥80 → TP rate 89.1%；低 AF (0.1-0.2) 信號最強（+50pp）
- **R5**: PairwiseMeanDist 與 HPFineN 正交（Spearman=0.07），微弱獨立增量
- **HPFineNGroups 確認為 somatic heterogeneity 標記** — 7/7 一致，residualized AUC=0.617，不能用於 filter 但有明確生物學價值

### Option C 雙路測試 NEGATIVE（2026-04-07 確定性結論）

- **架構確認**：cluster_labels 已是 HP-free（甲基化 hierarchical clustering），ClusterPermanova 被 passed_gate 擋住（TO 22% 有效）
- ClusterPermanovaF AUC=0.512（n=90,572）— 純甲基化 cluster 品質完全隨機
- HP-free 5 features combo AUC=0.564 vs HP-dependent combo AUC=0.598 vs 全部 AUC=0.607
- HP-free 特徵僅增加 +0.009 AUC — 無實用價值
- **結論：所有區分力來自 HP tags，純甲基化 clustering 無法突破 identifiability problem**
- C++ 代碼修改取消，確認正確策略為 PON-only phasing 上游修正

### O9 FN 特徵觀察 NO-GO（2026-04-08 確定性結論）

- 7 samples × 2 modes (Paired+TO)，122,790 FN regions 完整 ISM 執行
- **HP-free 甲基化特徵全 AUC < 0.53** — 純甲基化空間 FN 與 TP 不可區分
- 最強信號 LabelAllelePermanovaF=0.664（Paired）/ 0.636（TO）是 AF 代理非甲基化
- **TO Quality_Score AUC=0.338（嚴重反轉）** — FN 的 QS 比 TP 高，Cohen's d=-0.671
- FN VerificationClass: 55% Noise, 23% Weak, 18% Strong — 與 TP 分布相似
- **NO-GO：ISM 無法 rescue FN，甲基化空間 FN≡TP**

### TO-pure 獨立建模 NEGATIVE（2026-04-08 確定性結論）

- LOSO 7-fold：HP-free AUC=0.53、All-ISM=0.60-0.64、Caller-only=0.63、ISM+Caller=0.66
- caller_af (AUC=0.654) 單獨超越全部 ISM 特徵組合
- ISM 在 Caller 之上僅增加 +0.003（LR）~ +0.030（RF）
- HP-free ISM 完全隨機（AUC~0.53），HP-dependent +0.07~0.10 但可能循環
- **TO 模式 ISM 近乎無效，caller_af 是唯一有效判別器**

### 獨立分析完成（2026-04-11）

- **PON 跨樣本移除率驗證 ✓**：7 樣本 raw PON rate 95.16-98.81%（mean 97.77%），refined FP-level 全 > 98%。結論穩定度 3→4。H2009 最低（95.16%）因高突變負荷非 PON 失效。
- **H2009 負向根因 ✓**：Paired FP rate 僅 0.06%（86/132,994），改進天花板 +0.0004 F1。76.7% FP 在 LOH、89.5% Noise class。根因=caller 已近乎完美，ISM 只能誤傷 TP。甲基化特徵區分力反而優於平均。

### LOH / CN / AF 結論速查

完整三維度統合見 [07_LOH_CN_AF_研究總整理](reports/research_landscape/07_LOH_CN_AF_研究總整理.md)。

| 維度 | Filter 方向 | Characterization 方向 |
|------|-----------|---------------------|
| **LOH** | ❌ 全面關閉（10/10 策略失敗） | ✅ Subclone AF×Methylation 雙模式確認 |
| **CN** | ❌ Zone 排除全不可行 | ✅ Coverage_Multiple r=0.831 代理可信 |
| **AF** | ⚠️ 唯一有效信號但來自 caller | ✅ AF gradient 預測 NGroups |

### Zone-Aware Confidence Framework（2026-04-17 完整驗證）

- **構想文件**：[Zone-Aware Confidence Framework](concepts/2026/04/20260417_Zone_Aware_Confidence_Framework_01.md)
- **5 Zone 定義**：Z1 LOH Subclonal Active, Z2 High Somatic Hetero, Z3 Complete LOH, Z4 Normal Diploid, Z5 CN Gain Low Diversity
- **H1/H3 驗證**：[報告](../research/zone_aware_validation/20260417_H1_H3_Zone_Validation_Report_01.md)
- **Z3 內部特徵探索（2026-04-18 NEGATIVE）**：[報告](experiments/in_progress/2026/04/20260418_Z3_Internal_Feature_Exploration_01.md)
  - Step 1 Z3 內 12 特徵 × 7 樣本 AUC：無特徵在 ≥3 樣本達 |AUC|≥0.60
  - Step 2.5 AF∈[0.4,0.6] × CN × NGroups 分層：僅 HCC1954 (1/7) 符合 germline pattern（TP rate=0.146）
  - Step 3 HCC1954 vs HCC1395 機制對比：HCC1954 Z3 FP 集中 chr5/8/17（HER2/MYC amplicon），FP NumReads=55 vs TP=37（p=4.7e-9）→ CNV amplicon artifact 驅動；HCC1395 均勻分佈
  - **結論**：Z3 內無跨樣本二階區分特徵；HCC1954 例外已由 F pilot canonical filter 覆蓋
  - H1 CONDITIONAL → Z1b 放寬後 TO 4.6% 覆蓋率、TP rate 0.965（7 變體 Pareto 最佳）
  - H3 PARTIAL：Paired 7/7 ≥ 89.1% 確認；TO 6/7 significant 但絕對值 ~72%
  - TO zone TP rate 範圍 0.61-0.94
- **QS 模擬 ❌ NEGATIVE**：[報告](../research/zone_aware_validation/20260417_QS_Simulation_Report_01.md)
  - 5 configs × 21 thresholds × 7 samples，max delta F1=+0.001
  - **根因：TO QS AUC=0.497 隨機，zone delta 無法修正**
  - QS 調整路線和 C++ 整合**暫停**
- **結論：Zone-Aware 價值確認僅在 characterization annotation，不在 F1 改進**

### 待完成項目

- ~~Zone Z1 定義放寬~~ — ✅ 完成，Z1b 為最佳（TO 4.6%, TP rate 0.965）
- ~~Zone-Aware QS 調整模擬~~ — ❌ NEGATIVE，max delta F1=+0.001
- haplotag + ISM 全量重跑（7 samples × paired + TO）— PON-only phasing 後
- Phase 2A Normal Methylation Reference baseline — 依賴重跑數據
- LOH Subclone 更精細 AF-bin 分析 — 現有數據即可
- 7 樣本全量 Phase 2A 驗證 — 依賴重跑數據
- Platform normalization（5kHz / DORADO / PAO / Google）

## 4. 阻塞與風險

1. `archive/deep/` 為 immutable 快照，保留歷史失效連結（不回寫）。
2. `purity-aware` / subsample 結果易受 tumor-normal 組織甲基差異混淆，暫不作主證據。
3. `pileup` 與 `full model` 的 TO 輸入不可混用，報告必須明確標示來源。
4. TO borderline rescue 必須使用 candidate-specific `lost_tp / removed_fp` run。

## 5. 每次研究啟動必查

1. `docs/CURRENT_FOCUS.md`（本檔案）
2. `docs/experiments/INDEX.md`
3. `docs/README.md`
4. `docs/concepts/2026/04/20260409_研究構想總索引_01.md`（研究大圖景、發展樹、理論基礎）
5. `docs/reports/research_landscape/00_INDEX.md`（需深度理解時）

## 6. 已完成研究概覽

所有已完成研究的詳細記錄已封存至 `docs/archive/2026/03/20260405_current_focus_completed_items.md`。

**關鍵里程碑**：
- **LOH 雙定義與特徵全面關閉**（166 圖表 × 16 判定）：SEQC2 Jaccard=0.928 驗證 LOH.bed 可信；10/10 filter FAIL；Non-LOH max AUC<0.58；ISM 定位轉向 characterization
- **O1-O10 系統性觀察**（82 圖表）：TO 無單一特徵 AUC>0.58；LOH penalty 是 QS 根因；Paired/TO 根本不同
- **O11-O13 甲基化假說**：三維度全 NEGATIVE（within-group heterogeneity / LOH scenario / cross-region correlation）
- **G1-G7 VCF 特徵**（48 圖表）：60+ 特徵全 AUC<0.64；TO germline FP post-hoc 識別正式關閉
- **Read-level 鑑別**：LOSO AUC=0.721（首次>0.70）但安全約束 FP removal=0%
- **ASM 定量**：POSITIVE — 32-66% SNV 位點有 ASM；ISM PERMANOVA 唯一捕捉 entropy imbalance
- **Phase 1A**：paired-pure delta F1=+0.0112 已鎖定；TO 模式負增益
- **Self-phasing 因果鏈**：CONFIRMED — 62% LOH 消失、somatic bias 17.3:1
- **SEQC2 CNV 分層觀察**（15 圖表 + 5 TSV）：Coverage_Multiple r=0.831 驗證可信；跨樣本 AUC ≤ 0.641；zone 排除全不可行；CNV 非特徵耗盡根因
