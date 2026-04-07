<!--
建立時間: 2026-01-12 00:00
更新時間: 2026-04-07 01:00
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
   - 目前狀態：規劃中，待 Self-phasing 修正後的重跑數據
   - 關鍵前置：PON-only phasing 驗證已完成（LOH.bed Jaccard=1.0, somatic bias 消除）

2. **Phase 2 方向 B**：Gene-level / mechanism-level evidence integration
   - 目前狀態：待 Phase 2A 完成後啟動

3. **Phase 2 方向 C**：CpG 功能分層與智慧選點
   - 目前狀態：待規劃

### Self-Phasing 修正後重跑

- PON-only phasing 已驗證：LOH.bed 不變、somatic bias 消除、N50 +99.7%、phased rate +23.6pp
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

### R1-R5 特徵設計研究（2026-04-07 確認）

- **R1**: CramersV 93% 為零 = 2×2 框架缺陷；HPMergedDelta 多群時反向；HPFineNGroups 已克服（AUC +0.125）
- **R2**: Excess groups 概念有效（跨子集統一 +0.059）但子集內無新信號，不需修改 C++
- **R3**: 結構清楚子集 AUC 反而下降 → **確認是 identifiability 問題而非特徵設計問題**
- **R4**: HPFineNGroups N≥4 + NR≥80 → TP rate 89.1%；低 AF (0.1-0.2) 信號最強（+50pp）
- **R5**: PairwiseMeanDist 與 HPFineN 正交（Spearman=0.07），微弱獨立增量
- **HPFineNGroups 確認為 somatic heterogeneity 標記** — 7/7 一致，residualized AUC=0.617，不能用於 filter 但有明確生物學價值

### 待完成項目

- haplotag + ISM 全量重跑（7 samples × paired + TO）— PON-only phasing 後
- Phase 2A Normal Methylation Reference baseline — 依賴重跑數據
- O9 FN 特徵觀察（需先跑 FN ISM 數據）
- TO-pure 獨立建模（與 paired-pure 分離）
- H2009 異質性深度診斷
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
4. `docs/reports/research_landscape/00_INDEX.md`（需深度理解時）

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
