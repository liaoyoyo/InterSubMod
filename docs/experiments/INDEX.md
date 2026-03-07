<!--
建立時間: 2026-03-05 10:00
目標: 實驗研究歷史主索引，聚合所有已試驗方向，供 AI 快速掌握研究脈絡
處理範圍: InterSubMod 專案 2025-11 至今的所有研究方向
關聯檔案:
  - docs/README.md
  - docs/CURRENT_FOCUS.md
  - docs/experiments/validated/
  - docs/experiments/finalized/
  - docs/ai_sessions/2026/02/20260209_InterSubMod專案完整分析報告_01.md
-->

# InterSubMod 實驗研究索引

> **AI 使用提示**：本索引為 Layer 2（歷史全貌）。從這裡了解哪些方向已探索、結論為何，再決定下一步。

## 狀態圖例

- ✅ 成功完成，結論明確
- ❌ 失敗或無效，不建議重試
- ⏳ 進行中，尚無定論
- 🔄 值得再探索，有改進空間

---

## Layer 1：研究方向總覽

### 01. 甲基化解析（MM/ML 標籤）

- **期間**: 2025-11
- **狀態**: ✅ 成功
- **關鍵結論**: 支援 ONT BAM 的 `MM`/`ML` 標籤解碼，含正反股 CIGAR 座標校正，CpG 定位正確
- **主要文件**: `docs/solutions/debugging/2025/11/20251130_methylation_parser_debug_01.md`
- **建議後續**: 補齊 MethylationParser 獨立單元測試（目前僅有整合測試）

### 02. CIGAR 座標映射

- **期間**: 2025-11
- **狀態**: ✅ 成功
- **關鍵結論**: 反向股需從 SEQ 末端往回迭代匹配 MM delta encoding；`ref[offset-1]='C' && ref[offset]='G'` 驗證邏輯正確但非顯而易見
- **主要文件**: `docs/solutions/debugging/2025/11/20251130_mm_ml_tags_analysis_01.md`
- **建議後續**: 建立 mock BAM record 測試 forward/reverse 各種 CIGAR 邊界情境

### 03. 距離計算與聚類分析

- **期間**: 2025-11 ~ 2025-12
- **狀態**: ✅ 成功
- **關鍵結論**: 實作 NHD / L1 / L2 / CORR / BERNOULLI / JACCARD 六種距離度量；UPGMA 為預設聚類方法
- **主要文件**: `docs/archive/2025/12/20251218_final_completion_report_01.md`
- **建議後續**: 考慮 nearest-neighbor chain 加速聚類（n>500 時 O(n³) 可能成瓶頸）

### 04. OpenMP 平行化

- **期間**: 2025-12
- **狀態**: ✅ 成功
- **關鍵結論**: Per-region 平行化，thread-local BamReader/FastaReader 避免鎖競爭；單 Region 平均 < 300ms
- **主要文件**: `docs/architecture/system_overview.md`
- **建議後續**: 相鄰 SNV window 重疊 reads 可考慮共享（可減少 30-50% BAM I/O）

### 05. 統計顯著性：Fisher / χ²

- **期間**: 2025-12
- **狀態**: ✅ 成功
- **關鍵結論**: Fisher-Freeman-Halton (RxC) 使用 Patefield 算法，Monte Carlo early stopping 有效節省計算；目前主流程的核心統計方法
- **主要文件**: `docs/archive/2025/12/20251202_phase3_report_01.md`
- **建議後續**: 加入多重檢驗校正（FDR/Benjamini-Hochberg），30K+ region 獨立測試存在 false discovery 風險

### 06. 統計顯著性：PERMANOVA

- **期間**: 2025-12 ~ 2026-01
- **狀態**: ✅ 已實作，⚠️ 預設關閉
- **關鍵結論**: 已實作 pseudo-F = SS_Between/SS_Within；目前主流程 `enable_permanova=false`（效能考量）；啟用時需完整距離矩陣（無 NaN）
- **主要文件**: `docs/experiments/validated/2025/12/20251230_bidirectional_verification_strategy_analysis_01.md`
- **建議後續**: 在 Label-First 驗證框架下（HP/Allele 分組）重新評估 PERMANOVA 的適用時機

### 07. Bernoulli 距離度量

- **期間**: 2025-12
- **狀態**: ✅ 成功
- **關鍵結論**: `w(p)=2×|p-0.5|`（信心權重），高信心位點才對距離有大貢獻；比 NHD 更適合 ONT 中間不確定值
- **主要文件**: `docs/plans/2025/12/20251222_Bernoulli_Implementation_Plan_01.md`
- **建議後續**: 在 purity 系列測試中比較 Bernoulli vs NHD 的效能差異

### 08. TP/FP 特徵富集分析

- **期間**: 2026-01
- **狀態**: ✅ 成功（見 Layer 2 深入摘要）
- **關鍵結論**: LOH_like 富集 TP（OR=7.33），HighCluster 富集 FP（OR=0.02），LowAF 富集 FP（OR=0.05）
- **主要文件**: `docs/reports/validated/2026/03/20260301_甲基化欄位對TPFP與subclone驗證比較_01.md`
- **建議後續**: 對全量 TP/FP（而非 448 筆抽樣）重算 OR/F1，補充統計穩健性

### 09. 雙向驗證（Label-First vs Cluster-First）

- **期間**: 2026-01
- **狀態**: ✅ 框架設計完成（見 Layer 2 深入摘要）
- **關鍵結論**: 三角互證框架：PERMANOVA（Label→Structure）+ Fisher（Cluster→Label）；4 類輸出：Strong / Subclone / Weak / Noise
- **主要文件**: `docs/experiments/validated/2025/12/20251230_bidirectional_verification_strategy_analysis_01.md`
- **建議後續**: 補齊 Bootstrap stability check 實作；LabelTest 多階段 HP 驗證需獨立單元測試

### 10. F1 最佳化

- **期間**: 2026-01
- **狀態**: ✅ 成功（見 Layer 2 深入摘要）
- **關鍵結論**: 最優策略 `AF≥0.3 OR VerificationClass=Subclone`，F1=0.8481（基線 0.8155）；LabelDelta>0.3 可微幅提升至 0.8158
- **主要文件**: `docs/experiments/validated/2026/01/20260107_F1_Optimization_Deep_Analysis_01.md`
- **建議後續**: 在 COLO829 + purity series 做交叉驗證；建立 ML 分類器整合多特徵

### 11. Subsample 混樣甲基化偏差

- **期間**: 2026-02 ~ 2026-03
- **狀態**: 🔄 初步驗證完成，待全基因組結果（見 Layer 2 深入摘要）
- **關鍵結論**: HCC1395 ONT subsample 無 MM/ML；DORADO 三路對照（chr19，737 regions）顯示 Tumor < Normal < Mixed 距離梯度，初步支持 H3（組織差異假設）；差異效應量中等（0.04~0.09）
- **主要文件**:
  - `docs/ai_sessions/2026/03/20260302_subsample混樣甲基化偏差_現況研究推論與驗證路線圖_01.md`
  - `docs/experiments/in_progress/2026/03/20260305_DORADO三路對照實驗_01.md`
- **建議後續**: 等待全基因組三路完成 → 加入 HP haplotagging → 計算 Label-First 統計量

### 13. 多樣本甲基化距離基準線（2026-03）

- **期間**: 2026-03
- **狀態**: ✅ chr19 完成，全基因組待排
- **關鍵結論**: 6 個 ONT 癌症樣本均有 MM/ML；HCC1395 DORADO 距離明顯高於其他（0.316 vs 0.097~0.215）；H2009（肺腺癌）顯著性率最高（3.8%，109/2843）；HCC1937 距離最低（0.097）
- **主要文件**: `docs/reports/validated/2026/03/20260305_多樣本自動化測試總覽_01.md`
- **建議後續**: 加入 haplotagging 後比較各樣本的 HP-based 顯著性，建立跨樣本 F1 評估

### 14. no-filter 模式 Segmentation Fault 修正（2026-03）

- **期間**: 2026-03
- **狀態**: ✅ 成功修正
- **關鍵結論**: `--no-filter` 模式跳過 FLAG 過濾，使 secondary 讀段（SEQ="*"）進入 MethylationParser，觸發非確定性 Segfault；修正為即使在 no-filter 模式下仍過濾 BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FDUP
- **影響**: H2009（含大量 secondary 讀段的 long-read BAM）現可穩定在 j=8 執行
- **主要文件**: `src/core/RegionProcessor.cpp`（直接修正）；`docs/experiments/in_progress/2026/03/20260305_Phase1_2_3執行結論總結_01.md`
- **建議後續**: 確認其他含 secondary reads 的 BAM（如混合樣本）是否有相同問題

### 15. 全樣本放寬標準拆解與 F1 提升驗證（2026-03）

- **期間**: 2026-03
- **狀態**: ✅ 成功
- **關鍵結論**: 全樣本單一標準最佳為 `AD>0.15 & CV<0.03 & VAF<0.08`（6/7 樣本提升或持平，mean F1 delta `+0.000061`）；單獨 `VAF` 或單獨 `AD` 無法穩定提升，需 `VAF + 甲基資訊(AD)` 組合才有正向增益（`AD+VAF` mean F1 delta `+0.000046`）
- **主要文件**:
  - `docs/experiments/validated/2026/03/20260305_全樣本放寬標準與拆解驗證_01.md`
  - `docs/reports/validated/2026/03/20260305_全樣本放寬標準拆解與TPFPF1比較報告_01.md`
  - `docs/reports/validated/2026/03/20260305_提升幅度偏小根因與LongPhaseTO驗證計畫_01.md`
- **建議後續**: 將 `AD+VAF` 作為規則核心，`CV` 由硬門檻改為輔助訊號；對 HCC1954 採樣本化門檻減少 TP 誤刪

### 16. 純樣本 round 自動化實際執行（2026-03-07）

- **期間**: 2026-03-07
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - 既有 `20260211` paired pure run 已成功接上新 round automation，並確認 `HCC1395 5kHz` 比 `DORADO` 更 noisy
  - 新版完整欄位正式 rerun 已完成：
    - `HCC1395 5kHz`：`ClairS F1=0.8443`、`LongPhase-S F1=0.8522`、`InterSubMod F1=0.8532`
    - `HCC1395_DORADO`：`ClairS F1=0.8565`、`LongPhase-S F1=0.8592`、`InterSubMod F1=0.8590`
  - `label-first` 在代表性 diagnostics 上有提供新訊息，不只是 cluster 包裝
  - `低 VAF + 偏高 AlleleDelta + 不穩定全域結構` 是目前最重要的候選特徵，但只在 `5kHz` 顯示正增益，尚不具跨平台可攜性
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260307_純樣本round執行與進度報告_01.md`
- **建議後續**:
  - 對 `5kHz` 規則觸發子集補做更多 region diagnostics，確認是否為平台特異性 artifact 特徵
  - 比較 `5kHz` 與 `DORADO` 的 `label_cluster_agreement.tsv`，找出只在 `5kHz` 大量出現的升級型 region
  - 暫時不要把 `AD+VAF` 規則升級為全域預設，先視為 `5kHz` 特化診斷候選

### 17. 低 VAF + 高 AlleleDelta 與升級型 shift 跨樣本分析（2026-03-07）

- **期間**: 2026-03-07
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - 跨 7 個 pure paired 樣本比較後，`低 VAF + 偏高 AlleleDelta` 類規則只有 `HCC1395 5kHz` 有明顯正增益（`+0.000990`）
  - `HCC1395_DORADO` 與其他樣本多數為負增益或接近無效，因此不適合作全域規則
  - `Weak->Strong` 雖然整體仍以 TP 為主，但 `HCC1395 5kHz` 的 FP rate 最高（`172/627 = 27.4%`）
  - `Noise->Strong` 幾乎只出現在 `HCC1395` 系列，`5kHz` 高於 `DORADO`
  - `label-first` 不是無效，但 `Strong` 類別可能需要再細分為較可信與較可疑的子類
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260307_低VAF高AlleleDelta與shift跨樣本分析_01.md`
- **建議後續**:
  - 將 `low VAF + high AlleleDelta` 暫時定位為 `5kHz` 特化候選規則或 `ArtifactSuspect` 訊號
  - 對 `Noise->Strong / Weak->Strong` 的 FP top regions 補做更多 diagnostics 與 samtools snapshot
  - 規劃 tumor-only benchmark bundle，驗證這組特徵是否也能抓 TO artifact

### 18. Strong 細分類與 haplotagged samtools 驗證（2026-03-07）

- **期間**: 2026-03-07
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - `Strong` 類別應先做細分，再決定是否納入規則，不宜直接將 `low VAF + high AlleleDelta` 寫成全域過濾
  - 在 `HCC1395 5kHz`，`Strong + low VAF + high AlleleDelta` 可移除 `68` 個 FP、僅誤傷 `1` 個 TP，`F1 +0.000806`
  - 在 `HCC1395_DORADO` 與其他 pure paired 樣本，此規則沒有跨平台正增益，多為無效或誤傷 TP
  - `Noise->Strong` 不能直接當刪除條件；在 `HCC1395 5kHz` 會移除 `1659` 個 TP、僅移除 `51` 個 FP
  - haplotagged BAM 的 samtools snapshot 顯示 suspect Strong 至少包含兩型：
    - 單側 HP 偏倚型
    - 高未分派 HP 型
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260307_Strong細分類與samtools驗證分析_01.md`
- **建議後續**:
  - 先把 `Strong-SuspectLowVAF` / `ArtifactSuspect` 做成 annotation，而不是直接改核心規則
  - 補上 `HP balance` 與 `NA HP ratio` 欄位，測試是否能減少 5kHz 唯一誤傷的 TP
  - 以 `HCC1395 5kHz` 建立 tumor-only pilot bundle，驗證 refined label 是否仍富集 FP

### 19. LongPhase 救回 TP 與甲基救援分析（2026-03-08）

- **期間**: 2026-03-08
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - `HCC1395 5kHz` 上，`LongPhase-S` 除了減少 FP，也誤刪了一批 caller TP
  - 在精確 key 對齊集合中，最佳 TP rescue 規則目前是 caller `GQ >= 20`，可救回 `45` 個 TP、`0` 個 FP，`F1 +0.000739`
  - 將 InterSubMod `VerificationClass / PairwiseDist / AlleleDelta` 納入 rescue 後，尚未超過 caller `GQ` 單獨使用
  - 目前甲基訊號對 `LongPhase-S` 的價值更偏向 FP triage 與 Strong 細分類，而非 TP rescue 主規則
  - `HCC1395 ONT ClairS-TO` 已有 benchmark tp/fp/fn，可作為後續 `LongPhase-TO + InterSubMod` 的 tumor-only pilot 起點
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260308_LongPhase救回TP與甲基救援分析_01.md`
- **建議後續**:
  - 回查 `LongPhase-S` 對高 GQ caller TP 的丟棄條件
  - 對 `GQ >= 20` 但被丟掉的 TP 補做 read-level 與 haplotagging 驗證
  - 先以 `HCC1395 ONT ClairS-TO` 建立 `LongPhase-TO + InterSubMod` tumor-only pilot

### 20. LongPhase-TO 空間需求與中間產物確認（2026-03-07）

- **期間**: 2026-03-07
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - 真正的 `LongPhase-TO + InterSubMod` 驗證仍需要 `haplotagged BAM`；目前 repo 的 `02_intersubmod.sh` 直接要求 `--tagged-bam`
  - 本次 `HCC1395 5kHz TO` 實跑後，`tagged BAM` 實測大小約 `260G`，`.bai` 約 `113M`
  - 實際驗證後，當 `/bip8_disk` 只剩約 `208G` 時，已不足以安全承接 `5kHz TO haplotag`；需改寫到 `/big8_disk`
  - `HCC1395 ONT` 雖已有 TO phased / benchmark 素材，但因原始 BAM 無 `MM/ML`，只適合作 caller/phasing 邏輯研究，不適合作 InterSubMod 主驗證
  - `HCC1395_DORADO` 已補齊 TO baseline benchmark，`ClairS-TO F1=0.7226`，目前是最合理的 methylation-capable TO pilot 候選
  - 既有 `tmp/phasing_output/merged.vcf.gz` 不是可直接拿來宣稱 LongPhase-TO 成效的最終 somatic callset；直接 benchmark 其 non-missing 子集會得到失真的極低 F1
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260307_LongPhaseTO_空間需求與中間產物確認_01.md`
- **建議後續**:
  - 使用者已改為先做 `HCC1395 5kHz tumor-only` 的 `LongPhase-TO + InterSubMod` pilot，`HCC1395_DORADO` 改為第二階段交叉驗證
  - `5kHz TO` 的真實起點不是 `LongPhase-TO`，而是先補 `ClairS-TO` caller 輸出
  - 啟動下一個 `200G+` 級 tagged BAM 任務前，先把既有 `5kHz TO tagged BAM` 移到 `Archive_pending_delete/`
  - tagged BAM 產出後，再正式把 `低 VAF + 高 AlleleDelta`、`GQ rescue`、`Weak->Strong / Noise->Strong` 搬到 tumor-only 驗證
  - `HCC1395 5kHz` 跑完並完成比較後，依實驗室規則將 tagged BAM 移至 `Archive_pending_delete/`，再切 `HCC1395_DORADO`

### 21. HCC1395 5kHz tumor-only LongPhase-TO pilot 啟動（2026-03-07）

- **期間**: 2026-03-07
- **狀態**: ⏳ 主流程完成，後續觀察進行中
- **關鍵結論**:
  - 使用者已將 TO 主線優先順序改為：`HCC1395 5kHz` 先跑，`HCC1395_DORADO` 第二階段驗證
  - 因 `5kHz` 沒有現成 `ClairS-TO` 輸出，真實起點是先補 caller，再接 `LongPhase-TO`
  - 已新增可重用流程：
    - `scripts/pipeline/utils/benchmark_split_snv_vcf.sh`
    - `scripts/analysis/run_longphase_to_intersubmod_pilot.sh`
  - 這輪新增的關鍵修正：
    - benchmark split 可自動處理 plain `.vcf`
    - pilot 腳本可用 `--tag-prefix` 將 `tagged BAM` 寫到 `/big8_disk`
  - 新 benchmark split 工具已先用 `HCC1395_DORADO ClairS-TO` 測通，計數與既有報告一致
  - 第一次正式 run 因 `docker` 未掛載 `/bip8_disk` 而在 `step01` 早期失敗；修正後第二次 run 已完整跑通 `ClairS-TO -> LongPhase-TO -> haplotag -> InterSubMod`
  - `HCC1395 5kHz TO` 本輪結果：
    - `ClairS-TO F1=0.7127`
    - `LongPhase-TO F1=0.7127`
    - `InterSubMod F1=0.7130`
    - `FP removed=36`
    - `TP removed=3`
    - `F1 delta=+0.0003`
  - `LongPhase-TO phase` 本身沒有改變 benchmark call set；這輪的正增益主要來自甲基過濾層
  - `label-first` 在 TO 下仍有訊號：
    - `label_upgrade=2756`
    - `conflict=3227`
    - `consistent_strong=1471`
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260307_HCC1395_5kHz_TO_pilot啟動與執行紀錄_01.md`
- **建議後續**:
  - 對 `5kHz TO` 的 `low VAF + high AlleleDelta`、`Weak->Strong / Noise->Strong` 做 TP/FP 細分與 diagnostics
  - 跑完後直接與：
    - `ClairS-TO baseline`
    - `LongPhase-TO`
    - `InterSubMod`
    - 既有 paired pure `HCC1395 5kHz`
    做同口徑比較
  - 比較完成後將 5kHz tagged BAM 移至 `Archive_pending_delete/`，再切 `HCC1395_DORADO`

### 22. TO 特徵回查與 label_upgrade / conflict 細分（2026-03-08）

- **期間**: 2026-03-08
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - `HCC1395 5kHz TO` 的 `label_upgrade=2756` 中，`TP=2526`、`FP=230`，`FP rate=8.35%`
  - `conflict=3227` 中，`TP=2465`、`FP=762`，`FP rate=23.61%`
  - `Weak->Strong`（`1849 TP / 176 FP`）與 `Noise->Strong`（`207 TP / 16 FP`）在 TO 下都以 `TP` 為主，**不能直接當刪除規則**
  - TO 下真正持續有效的 artifact 特徵仍是 `low VAF + high AlleleDelta`：
    - `TP=3`
    - `FP=36`
    - `FP rate=92.31%`
  - 再加 `low CramersV` 後更乾淨：
    - `TP=2`
    - `FP=36`
    - `F1 +0.000261`
  - 這批 TO artifact 幾乎**不落在** `label_upgrade/conflict`，而是主要落在：
    - `cluster_only`
    - `Strong->Weak`
    - `低 GQ / 低 QUAL`
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260308_TO特徵回查與label_upgrade_conflict分析_01.md`
- **建議後續**:
  - 下一輪 TO diagnostics 應優先看 `cluster_only + Strong->Weak + low VAF + high AlleleDelta`
  - 不要把 `label_upgrade`、`conflict`、`Weak->Strong`、`Noise->Strong` 整包升級成 artifact 規則
  - 直接用同一支 TO 回查腳本套到 `HCC1395_DORADO TO`，確認這個模式是否可跨平台重現

### 23. TO `cluster_only + Strong->Weak` diagnostics 與 verification scheme 調整（2026-03-08）

- **期間**: 2026-03-08
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - 本輪 `HCC1395 5kHz TO` diagnostics 使用的是 `ClairS-TO pileup` 路線，不是 full model；結論不能直接與其他 full-model run 混比
  - 原先鎖定的 TO top FP 中，`4/5` 聚集在 `chr11` 的 `1402 bp` 區段，顯示更像局部 artifact block，而不是彼此獨立的隨機 FP
  - `HP skew` 在 TO 下不夠區分 FP/TP，因為 TP 對照也可呈現極端單一 HP；不能單靠 `HP balance` 做刪除
  - 比較有辨識力的是：
    - `高 AlleleDelta`
    - `低 PairwiseMedianDist`
    - `很低 alt fraction`
    - `低 QUAL / GQ`
  - 現行 core 的 `Strong` 會混入「只有 allele sig、沒有 HP sig」的位點；這批 `current_strong_allele_only` 為 `5104 TP / 2042 FP`，不能整包視為 artifact
  - 新增 `cluster_plus_weak_label` 後，原先的 TO 可疑子集更精確地落在：
    - `cluster_plus_weak_label + Strong->Weak + lowVAF/highAD`
    - 其中 `lowVAF/highAD + low CramersV` 為 `25 FP / 1 TP`
  - verification scheme 的語意調整目前主要提升可解釋性；F1 仍未超過既有 `low VAF + high AlleleDelta + low CramersV`（`36 FP / 2 TP`, `F1 +0.000290`）
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260308_TO_cluster_only_StrongWeak_diagnostics與scheme調整分析_01.md`
- **建議後續**:
  - 後續 TO diagnostics 應將舊語意 `cluster_only + Strong->Weak` 修正為 `cluster_plus_weak_label + Strong->Weak`
  - 先把同一套 diagnostics 與 scheme evaluator 套到 `HCC1395_DORADO TO`
  - 若未來要改 core C++ 分類，優先考慮新增 annotation（例如 `Strong_AlleleOnly`），不要直接把整包 `Strong` 當過濾目標

### 24. ClairS 邊緣 TP rescue 與甲基輔助評估（2026-03-08）

- **期間**: 2026-03-08
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - 本輪固定只使用 **final VCF**，不使用 pileup 中間 candidate / tensor
  - `HCC1395 5kHz paired` 與 `HCC1395_DORADO paired` 都可用 caller-only 規則穩定救回一批 TP
  - 依計劃排序的最佳 safe rule 在 paired 兩條線都是 `candidate_gq_ge_15`：
    - `HCC1395 5kHz`: `106 TP / 75 FP`, `F1 +0.000825`
    - `HCC1395_DORADO`: `97 TP / 88 FP`, `F1 +0.000502`
  - 若單看最佳 `delta F1`，paired 兩條線都是 `gq>=20` 類規則更好：
    - `HCC1395 5kHz`: `59 TP / 8 FP`, `F1 +0.000871`
    - `HCC1395_DORADO`: `53 TP / 7 FP`, `F1 +0.000782`
  - `HCC1395 5kHz TO` 的 caller-only rescue 空間更大：
    - `candidate_qual_ge_10_or_gq_ge_20`: `491 TP / 118 FP`, `F1 +0.006824`
    - `candidate_gq_ge_15`: `365 TP / 71 FP`, `F1 +0.005233`
  - 目前甲基 rescue **尚未超過 caller-only**
  - `HCC1395 5kHz paired` 雖已有部分 candidate-specific InterSubMod coverage，但只覆蓋 `111/1052` lost TP、`802/12974` removed FP
  - `HCC1395_DORADO paired` 尚缺 candidate-specific InterSubMod summary，因此本輪不能把甲基 rescue 解讀成有效 negative result
  - `HCC1395 5kHz TO` 的 downstream lost/removed 與現有 baseline InterSubMod summary overlap 為 `0/0`；若要判斷 TO 下甲基是否能幫助 borderline TP rescue，必須先對這批 TO candidates 單獨跑 InterSubMod
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260308_ClairS邊緣TP_rescue與甲基輔助評估_01.md`
- **建議後續**:
  - 先跑 `HCC1395 5kHz TO` 的 candidate-specific `lost_tp / removed_fp` InterSubMod
  - 再補 `HCC1395 5kHz paired` 的 final VCF 全候選池 coverage，不只沿用舊的 PASS 子集
  - 對 `HCC1395_DORADO paired` 跑同型 candidate-specific InterSubMod，確認甲基 rescue 是否真為平台不適用
  - 若要進 pipeline，先加 `RescueCandidate` annotation，不先直接加 hard keep / hard filter

### 25. ClairS 邊緣 FN 探勘與甲基化救援（2026-03-08）

- **期間**: 2026-03-08
- **狀態**: ✅ 結論完整，不推薦升為正式規則
- **關鍵結論**:
  - `Pool B FN`（ClairS non-PASS 但 truth set 包含）：**840 個**，遠超預期（計劃估 30~500）
  - 策略 B 補跑 InterSubMod 後甲基覆蓋率達 **99.9%**（839/840）
  - **甲基唯一規則全數為負效益**：`pairwise_ge_020`（F1 delta=-0.006968），`class_strong_or_subclone`（-0.006265）
    - 根本原因：Pool B non-FN 含大量 germline 變異，Strong rate 反而比 FN 更高（44.7% vs 30.8%）
  - **最佳 Caller-only 規則**：`no_varcluster`（F1 delta=+0.003391）、`no_varcluster_and_gq15`（precision=71.9%, F1 delta=+0.001452）
  - **最佳 Combined 規則**：`gq15_and_allele_delta_low`（F1 delta=+0.002702），僅比最佳 Caller-only 好 +0.000386
  - **AlleleDelta 的 somatic 特異性**：Pool B FN 中位 AD=0.010（90% < 0.10），non-FN 中位 AD=0.065（40% ≥ 0.10）——可作為輔助訊號，但 precision 仍只 ~32%
  - **結論與實驗 19 一致**：甲基訊號對 TP rescue 貢獻受限，主要價值在 FP triage
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md`
  - `scripts/analysis/analyze_clairs_borderline_fn.py`
  - `scripts/analysis/run_clairs_borderline_fn_analysis.sh`
- **建議後續**:
  - `no_varcluster_and_gq15` 在 `HCC1395_DORADO TO` 交叉驗證（precision=71.9%，值得確認是否跨平台）
  - 研究重點維持 FP triage（`low VAF + high AlleleDelta` 已有明確正效益）
  - Pool A-light 不需再探索

### 12. Purity-Aware 過濾策略

- **期間**: 2026-02 ~ 2026-03
- **狀態**: 🔄 值得再探索（見 Layer 2 深入摘要）
- **關鍵結論**: 低 purity（<40%）時固定甲基化門檻嚴重傷害 Recall；需依 purity 分層的自適應策略
- **主要文件**: `docs/plans/2026/02/20260228_InterSubMod再驗證與再實驗執行計畫_01.md`
- **建議後續**: 低 purity 停用甲基化過濾，中/高 purity 保留；建立 purity-aware feature 融合模型

---

## Layer 2：重點方向深入摘要

> 以下 5 個方向為目前最具研究價值或待解問題，提供完整數據與失敗嘗試紀錄。

---

### 深入：TP/FP 特徵富集分析

**研究問題**: InterSubMod 的甲基化特徵是否能有效區分 Somatic SNV 的真陽性 (TP) 與假陽性 (FP)？

**方法**:
- 資料集：448 筆分層抽樣（TP 288、FP 160），HCC1395 + SEQC2 truth set
- 定義特徵：LOH_like, HighCluster（ClusterCount50kb>20）, LowAF（AF<0.3）, MethylSupport_SubcloneLike
- 統計：Fisher exact test（OR + p-value）+ 規則式 F1 評估

**數據結果**:
- `LOH_like` 富集 TP：OR=7.33, p=3.04e-16，TP rate=74.2%（非條件=28.1%）
- `HighCluster` 富集 FP：OR=0.0185, p=5.75e-26，TP rate=4.9%（非條件=73.6%）
- `LowAF` 富集 FP：OR=0.0458, p=4.44e-40，TP rate=34.1%（非條件=91.9%）
- 低 AF 子群下，`VerificationClass=Subclone`：OR=8.13, p=9.88e-05（顯著 TP 富集）

**最佳 F1 策略**:
| 策略 | F1 |
|---|---|
| 基線（全部保留）| 0.8155 |
| `AF≥0.3` | 0.8238 |
| `AF≥0.3 OR VerificationClass=Subclone` | **0.8481**（最優） |
| `SignificantOnly` | 0.1910（無效）|

**失敗嘗試**: 單獨使用 CramersV 過濾（AUC=0.52）無效；SignificantOnly 策略嚴重損失 Recall

**關鍵限制**: 核心比較集為分層抽樣（448筆），非全量 30K+ 位點；全量結論需再驗證

**建議下一步**:
1. 對全量 TP/FP 重算 LOH_like / HighCluster 的 OR 穩健性
2. 建立 logistic / tree 模型，組合 `VerificationClass + HPMergedDelta + ClusterCount50kb`
3. 在 COLO829 做 held-out 交叉驗證

**相關文件**:
- `docs/reports/validated/2026/03/20260301_甲基化欄位對TPFP與subclone驗證比較_01.md`
- `docs/experiments/validated/2026/01/20260107_F1_Optimization_Deep_Analysis_01.md`
- `docs/experiments/validated/2026/01/20260107_F1_and_Data_Optimization_Report_01.md`

---

### 深入：Purity-Aware 過濾策略

**研究問題**: 如何設計甲基化過濾策略，使其在不同腫瘤純度（purity）下都能保持穩定效能？

**方法**:
- 資料集：HCC1395 DORADO subsample（t7_n29 ~ t50_n00，有 MM/ML）
- 比較固定門檻 vs purity 分層策略的 F1 差異

**數據結果**:
- 低 purity（~19.4%，t7_n29）：固定甲基化門檻顯著傷害 Recall（F1 下降 0.05 級）
- 中高 purity（≥60%）：過濾效果趨近中性或微幅正向
- 結論：固定單一門檻跨 purity 直接套用不可行

**建議策略**:
```
< 40% purity  → 不做甲基化過濾
40-60% purity → 保守門檻或僅標註不刪除
≥ 60% purity  → 啟用甲基化過濾
```

**失敗嘗試**: 固定 `Significant=True` 在低 purity 下等同於強制移除大量 TP

**關鍵限制**: 目前缺乏跨 purity + 跨樣本的完整驗證矩陣；COLO829 與其他樣本尚未納入

**建議下一步**:
1. 停用低 purity 甲基化過濾，重跑 F1 曲線確認恢復
2. 建立統一 `purity_qc.tsv`：每個 subsample 記錄 MM_count, ML_count, analyzable
3. 規劃 Experiment-02（多特徵融合）作為固定門檻的替代方案

**相關文件**:
- `docs/plans/2026/02/20260228_InterSubMod再驗證與再實驗執行計畫_01.md`
- `docs/ai_sessions/2026/02/20260209_InterSubMod專案完整分析報告_01.md`

---

### 深入：Subsample 混樣甲基化偏差

**研究問題**: 為何 HCC1395 subsample 分析效果不佳？是 tumor/normal 組織來源差異（甲基化背景）造成的系統性偏差？

**方法**:
- 比較 HCC1395 ONT（無 MM/ML）與 HCC1395 DORADO（有 MM/ML）兩批資料
- 分析 tumor（乳腺癌細胞株）vs blood-normal（HCC1395BL，B-lymphoblast）的甲基化背景差異

**數據結果**:
- HCC1395 ONT subsample：無 MM/ML → `tp_regions=0, fp_regions=0`，InterSubMod 完全無有效輸出
- HCC1395 DORADO：有 MM/ML，但低 purity 時固定門檻傷害 Recall
- 外部文獻（Cell Reports Methods 2025）確認：HCC1395 vs HCC1395BL 甲基化 overlap 約五成，組織差異顯著

**已採取措施（2026-03-02）**:
- `scripts/pipeline/run_benchmark.sh`：新增 Methylation Guard（MM/ML 前置檢查）
- `scripts/analysis/run_purity_and_standard_verification.sh`：新增 source BAM MM/ML 檢查
- `scripts/pipeline/config.sh`：新增共用函式 `has_mm_ml_tags()`

**推論（依可信度）**:
1. H1（高）：MM/ML 缺失是第一層阻斷，直接造成無有效統計區域
2. H2（高）：低 purity 下固定門檻將大量 TP 錯誤過濾
3. H3（中高）：tumor/blood-normal 甲基化背景差異放大 composition effect
4. H4（中高）：InterSubMod 主程式未實際把 normal 納入核心矩陣，對混合 BAM 無來源分層

**關鍵限制**: 目前缺乏 tumor-only / normal-only / mixed 三路直接對照數據

**建議下一步**:
1. 使用 DORADO 重建 tumor-only / normal-only / mixed 三組，同一批位點比較 distance 分布
2. 在 InterSubMod 加入來源分層標籤輸出（RG tag 或 BAM 來源標記）
3. 對每個 purity 產出固定 `purity_qc.tsv`

**相關文件**:
- `docs/ai_sessions/2026/03/20260302_subsample混樣甲基化偏差_現況研究推論與驗證路線圖_01.md`
- `docs/ai_sessions/2026/02/20260213_HCC1395_subsample_purity_完整驗證報告_01.md`（Knowledge 路徑）

---

### 深入：雙向驗證（Label-First vs Cluster-First）

**研究問題**: 如何更準確地判斷甲基化距離矩陣中的結構是否有生物學意義？

**方法**:
- **Path A (Label-First)**：使用已知 HP/Allele 標籤分組，PERMANOVA 計算 Pseudo-F 與 p-value，Delta Distance = Between - Within
- **Path B (Cluster-First)**：層次聚類後，Bootstrap stability check（80% 重抽 20-50次），Fisher's Test 計算 Cramér's V

**4 類輸出判定矩陣**:

| 類別 | Label-First | Cluster-First | 解釋 |
|---|---|---|---|
| **Strong** | 顯著（高 R²）| 一致且穩定 | 標籤是甲基化變異主因 |
| **Subclone** | 不顯著 | 顯著且穩定 | 真實結構，但非 HP/Allele 主導 |
| **Weak** | 邊緣顯著 | 不穩定 | 效應存在但不明顯 |
| **Noise** | 不顯著 | 不穩定 | 隨機分群假象 |

**目前實作狀態**:
- PERMANOVA (`StructureTest.cpp`)：已實作，但主流程預設 `enable_permanova=false`
- Fisher (`GlobalTest.cpp`, `LocalTest.cpp`)：已啟用
- Bootstrap stability：已設計，實作待補齊
- LabelTest 多階段：已實作（Merged HP, Fine-grained HP, Allele, Unassigned affinity）

**關鍵限制**: 現有 PERMANOVA 呼叫基於 `cluster_labels`（非 HP/Allele 標籤），不符合 Label-First 定義；需修正

**建議下一步**:
1. 修正 PERMANOVA 呼叫，改用 HP/Allele 標籤作為分組依據
2. 補齊 Bootstrap stability check 的 C++ 實作
3. 建立 LabelTest 多階段 HP 驗證的獨立單元測試

**相關文件**:
- `docs/experiments/validated/2025/12/20251230_bidirectional_verification_strategy_analysis_01.md`
- `docs/architecture/system_overview.md`（S7: SignificanceAnalyzer）

---

### 深入：統計顯著性分析（PERMANOVA + Cramér's V）

**研究問題**: 哪些統計指標能有效度量甲基化分群的生物學顯著性？

**方法**:
- 5 層統計框架：GlobalTest（Fisher FFH）→ LocalTest（One-vs-Rest）→ StructureTest（PERMANOVA）→ LabelTest（HP/Allele）→ Bootstrap
- Cramér's V 作為效應量：`V = sqrt(χ²/n × 1/(min(r,c)-1))`

**關鍵數據**:
- Cramér's V 的 TP/FP 區分能力：AUC=0.52（幾乎無效）
- QUAL（VCF 原生）：AUC=0.9668；AF：AUC=0.9235
- 甲基化特徵整體 AUC 低，但在低 AF（<0.3）子群：NumReads AUC=0.84（有條件價值）
- 甲基化顯著位點中 94.8% 為 TP（高信心保留依據）

**失敗嘗試**: 使用 `Significant=True` 作為主要過濾器（F1=0.19）；使用 CramersV>0.1 過濾（損失過多 TP）

**核心發現**: 甲基化特徵不應作為獨立過濾器，而應作為 VCF 品質不足時的輔助依據（分層策略）

**建議下一步**:
1. 實作多重檢驗校正（FDR / Benjamini-Hochberg）
2. 建立分層過濾策略：QUAL≥0.8 直接接受 → QUAL<0.8 且甲基化顯著 → 保留（救援）→ 否則過濾
3. 在 LabelTest 加入 per-HP PERMANOVA（Label-First 框架）

**相關文件**:
- `docs/ai_sessions/2026/02/20260209_InterSubMod專案完整分析報告_01.md`（§4.2, §7.1）
- `docs/experiments/validated/2026/01/20260107_F1_and_Data_Optimization_Report_01.md`

---

## 附錄：待驗證方向（尚未正式啟動）

| 方向 | 期望收益 | 難度 | 依據 |
|---|---|---|---|
| 5hmC 雙通道距離矩陣 | 可能更有亞克隆特異性 | 高 | ONT 5kHz 同時提供 C+h |
| 跨 Region 亞克隆一致性 | 真正的亞克隆應跨多 SNV 一致 | 高 | 目前 per-SNV 獨立分析 |
| 機器學習組合特徵分類器 | 整合 15 個特徵的 ensemble | 中 | F1 研究發現特徵互補性 |
| PMD/ChromHMM Gating 啟用 | 降低甲基化背景噪聲 | 中 | 架構已有 `is_pmd` 欄位但未生效 |
| 多樣本交叉驗證（COLO829 等）| 泛化能力評估 | 中 | 目前僅 HCC1395 |
