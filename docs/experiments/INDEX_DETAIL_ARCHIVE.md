<!--
建立時間: 2026-04-05 22:00
目標: 實驗研究歷史詳情封存（Layer 2 完整描述）
處理範圍: InterSubMod 專案 2025-11 至 2026-04 的所有研究方向完整記錄
關聯檔案:
  - docs/experiments/INDEX.md（精簡活動版）
  - docs/CURRENT_FOCUS.md
-->

# InterSubMod 實驗研究詳情封存

> 本檔案為 INDEX.md 的完整歷史快照（2026-04-05 封存）。
> 精簡版活動索引見 `docs/experiments/INDEX.md`。

---



# InterSubMod 實驗研究索引

> **AI 使用提示**：本索引為 Layer 2（歷史全貌）。從這裡了解哪些方向已探索、結論為何，再決定下一步。
>
> **2026-03-24 後的使用方式**：live 主線與近期任務請先看 [CURRENT_FOCUS.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md) 與 [20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md)。本索引負責歷史收斂、失敗方向盤點與證據入口，不再直接充當近期任務清單。

## 研究文件 Agent 與 Skills 入口

若任務是整理研究週報、主題式證據鏈或收斂研究文件偏好，先看：

1. [主 Agent：InterSubMod 研究文件代理](/big8_disk/liaoyoyo2001/InterSubMod/.claude/agents/intersubmod-weekly-research-agent.md)
2. [研究文件使用手冊](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究報告Agent與Skills使用手冊_01.md)
3. [個人化研究寫作偏好設定規格](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_個人化研究寫作偏好設定規格_01.md)

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
  - `docs/provenance/ai_sessions/2026/03/20260302_subsample混樣甲基化偏差_現況研究推論與驗證路線圖_01.md`
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

### 25. HCC1395 5kHz TO candidate-specific 甲基 rescue 驗證（2026-03-08）

- **期間**: 2026-03-08
- **狀態**: ✅ 第一輪結論完整
- **關鍵結論**:
  - `HCC1395 5kHz TO` 的 candidate-specific InterSubMod 已正式補跑完成，不再是舊版 `0/0 overlap`
  - candidate export 完整成功：
    - `caller_lost_tp=773`，missing `0`
    - `caller_removed_fp=298`，missing `0`
  - candidate-specific InterSubMod 可分析率：
    - `lost_tp=675/773 (87.32%)`
    - `removed_fp=213/298 (71.48%)`
  - **最佳整體 rescue 仍是 caller-only**：
    - `qual>=10 or gq>=20`：`491 TP / 118 FP`，`F1 +0.006824`
  - **甲基資料不是全無效**：
    - `PairwiseMedianDist >= 0.20`：`300 TP / 68 FP`，`F1 +0.004219`
    - 表示 TO 下甲基資料可以幫 borderline TP rescue，但尚未超過最佳 caller-only
  - **label-first / cluster-first 的交叉訊號有實際價值**：
    - `agreement_positive`：`148 TP / 25 FP`，`F1 +0.002163`
    - 比單純 `Strong/Subclone`（`149 TP / 30 FP`）更乾淨，顯示 agreement 可改善 support precision
  - **artifact veto 不適合直接拿來救 TP**：
    - `lowVAF/highAlleleDelta/lowCramersV` 對 `gq>=15` pool 沒有額外效果
    - `combined_artifact_veto` 雖少 `6 FP`，但同時少 `57 TP`
  - **總結**：TO 下甲基訊號最合理的定位是 `support / ranking / annotation`，不是取代 caller 的主 rescue 規則
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260308_HCC1395_5kHz_TO_candidate_specific甲基rescue驗證_01.md`
  - `docs/reports/validated/2026/03/20260308_HCC1395_5kHz_TO_borderline_rescue特徵證據鏈整理_01.md`
  - `scripts/analysis/export_candidate_pool_vcfs.py`
  - `scripts/analysis/evaluate_rescue_with_methylation.py`
  - `scripts/analysis/validate_method_design.py`
- **建議後續**:
  - `HCC1395_DORADO TO` 的同型 candidate-specific InterSubMod 已於 `2026-03-11` 完成，確認 TO 下跨 `5kHz` 與 `DORADO` 都支持 `caller-first + methylation-support`
  - 對 `agreement_positive` 與 `PairwiseMedianDist>=0.20` rescue 到的 TP/FP 補做 diagnostics
  - 若要進 pipeline，先做 `RescueSupport` annotation，不先做 hard keep

### 26. ClairS 邊緣 FN 探勘與甲基化救援（2026-03-08）

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

### 27. HCC1395_DORADO paired candidate-specific 甲基 rescue 驗證（2026-03-09）

- **期間**: 2026-03-09
- **狀態**: ✅ 第一輪結論完整
- **關鍵結論**:
  - `HCC1395_DORADO paired` 的 candidate-specific InterSubMod 已正式補齊，不再是舊版「缺 summary、無法判讀」的狀態
  - 本輪沿用與實驗 24 相同的 caller / truth / baseline 定義：
    - caller final VCF: `ClairS v0.4.0 output.vcf.gz`
    - baseline TP/FP: `s-pure/HCC1395_DORADO/20260211`
  - candidate export 完整成功：
    - `caller_lost_tp = 1489`，`candidate_eligible = 1122`
    - `caller_removed_fp = 2533`，`candidate_eligible = 1658`
    - export missing 全部為 `0`
  - 因本機大磁碟空間不足，本輪將 `LongPhase-S` 重建寫到 root filesystem 的 scratch：
    - 重建 `HCC1395_DORADO_tagged.bam` 成功，大小 `223G`
    - `bam.bai` 已完成
    - LongPhase-S tagging fraction 約 `0.469462`
  - candidate-specific 可分析 coverage 明顯偏低：
    - `lost_tp = 93 / 1122 (8.29%)`
    - `removed_fp = 326 / 1658 (19.66%)`
  - **caller-only 仍穩定成立**：
    - `candidate_gq_ge_15`：`97 TP / 88 FP`，`F1 +0.000502`
    - `candidate_gq_ge_20`：`51 TP / 7 FP`，`F1 +0.000749`
    - `caller_any_gq_ge_20`：`53 TP / 7 FP`，`F1 +0.000782`
  - **甲基 support 規則在 DORADO paired 全部不成立**：
    - `gq>=10 + PairwiseMedianDist>=0.20`：`10 TP / 36 FP`，`F1 -0.000280`
    - `gq>=10 + agreement_positive`：`36 TP / 72 FP`，`F1 -0.000298`
    - `gq>=10 + Strong/Subclone`：`25 TP / 38 FP`，`F1 -0.000059`
  - 可分析子集的分布也支持這個負結果：
    - `lost_tp` 中 `Strong/Subclone = 29`
    - `removed_fp` 中 `Strong/Subclone = 56`
    - `agreement_positive` 在 `lost_tp` 為 `37`，在 `removed_fp` 為 `99`
    - `Pairwise>=0.20` 在 `lost_tp` 為 `13`，在 `removed_fp` 為 `58`
  - **總結**：
    - `5kHz TO` 的甲基 rescue 是真實存在的
    - 但目前**不能升級成跨樣本 / 跨模式穩定規則**
    - `DORADO paired` 新證據支持：真正穩定的是 `caller-first`，不是 `methylation-support`
- **主要文件**:
  - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260309_HCC1395_DORADO_paired_candidate_specific甲基rescue驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260309_HCC1395_DORADO_paired_candidate_specific甲基rescue驗證_01.md)
  - [/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/eval/rescue_rule_comparison.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/eval/rescue_rule_comparison.tsv)
  - [/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/design_validation/label_cluster_agreement.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/design_validation/label_cluster_agreement.tsv)
- **建議後續**:
  - `HCC1395_DORADO TO candidate-specific InterSubMod` 已於 `2026-03-11` 完成；對外口徑可更新為：
    - `5kHz TO` 與 `DORADO TO` 都支持 `caller-first` 穩定成立
    - 甲基單特徵可正向 rescue，但在同一 caller gate 下仍只作 support
    - `PairwiseMedianDist` 方向仍具樣本依賴，尚未證明可跨樣本直接共用
  - Pool A-light 不需再探索

### 28. 5kHz TO 與 DORADO paired 的甲基 rescue 特徵空間分析（2026-03-09）

- **期間**: 2026-03-09
- **狀態**: ✅ 第一輪完成
- **關鍵結論**:
  - 已用新腳本對兩個 candidate-specific `rescue_joined_features.tsv` 做同一套 feature-space 分析：
    - 單特徵分佈
    - 單特徵 rescue sweep
    - 雙特徵組合 sweep
    - `delta_f1_vs_baseline`
    - `delta_f1_vs_gate`
  - `HCC1395 5kHz TO`：
    - candidate-specific coverage 高（`lost_tp=87.32%`, `removed_fp=71.48%`）
    - `Quality_Score>=60`、`PairwiseMedianDist>=0.15/0.20`、`agreement_positive` 都能從 baseline 救回 TP
    - 但在 `gq>=10` gate 下，**沒有任何甲基或 label 特徵能進一步提高 F1**
    - `agreement_positive` 仍是較乾淨的窄 support，但屬 ranking / annotation，不是加成規則
    - `allele_assign_rate>=0.99` 與 `CramersV<=0.05` 看似數值高，但實際更像技術 proxy，不應升級為正式甲基 rescue 特徵
  - `HCC1395_DORADO paired`：
    - candidate-specific coverage 低（`lost_tp=8.29%`, `removed_fp=19.66%`）
    - 正向 support 方向與 `5kHz TO` 不同：不是看高 `Pairwise`，而是偏向 `高 hp_assign_rate + 低 Pairwise + 較高 Quality_Score`
    - 單特徵最佳增量為 `gq>=15 + hp_assign_rate>=0.99`：`50 TP / 15 FP`，相對 gate `F1 +0.000132`
    - 雙特徵最佳增量為 `gq>=15 + Pairwise<=0.20 + hp_assign>=0.99`：`46 TP / 12 FP`，相對 gate `F1 +0.000103`
    - 即使如此，整體仍沒有超過 `gq>=20` caller-only
  - **總結**：
    - `5kHz TO` 與 `DORADO paired` 的甲基 support 方向不一致
    - 目前仍找不到可跨樣本共用的單一 rescue 規則
    - 最合理定位仍是 `caller-first + mode-specific methylation support`
- **主要文件**:
  - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260309_5kHz_TO與DORADO_paired_甲基rescue特徵空間分析_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260309_5kHz_TO與DORADO_paired_甲基rescue特徵空間分析_01.md)
  - [/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_methylation_rescue_feature_space.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_methylation_rescue_feature_space.py)
  - [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/top_rules_by_gate.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/top_rules_by_gate.tsv)
  - [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/numeric_feature_distribution.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/numeric_feature_distribution.tsv)
- **建議後續**:
  - `5kHz TO`：針對 `Quality>=60`、`Pairwise>=0.15/0.20`、`agreement_positive` rescue 到的 top TP/FP 補 diagnostics
  - `DORADO TO`：優先驗證 `hp_assign>=0.99` 與 `低 Pairwise + 高 hp_assign` 是否在 TO 重現
  - 暫時不要再把 `Pairwise>=0.20` 當作全域規則
  - 也不要把 `allele_assign_rate` 或 `CramersV<=0.05` 升級成正式甲基 rescue 規則

### 29. GQ 與甲基 rescue 的系統性驗證（2026-03-11）

- **期間**: 2026-03-11
- **狀態**: ✅ 第一輪完成
- **關鍵結論**:
  - 已新增同一套系統性分析，直接比較：
    - `GQ only`
    - `甲基 only`
    - `GQ + 甲基`
    - `artifact role`
    - 區間性與正交特徵
  - 新腳本：
    - [/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_gq_methylation_rescue_matrix.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_gq_methylation_rescue_matrix.py)
  - `HCC1395 5kHz TO`：
    - `GQ` 仍是最穩定的 caller-first 訊號
    - `gq>=10`：`499 TP / 119 FP`，`F1 +0.006943`
    - `gq>=15`：`365 TP / 71 FP`，`F1 +0.005233`
    - `gq>=20`：`143 TP / 38 FP`，`F1 +0.001966`
    - `Pairwise>=0.20`、`Quality>=60`、`agreement_positive` 都能單獨從 baseline 救回 TP，但全部 **沒有超過 `gq>=10`**
    - `gq>=10 + Pairwise>=0.20`：`300 TP / 68 FP`，相對 gate `delta F1=-0.002724`
    - `gq>=10 + agreement_positive`：`148 TP / 25 FP`，相對 gate `delta F1=-0.004780`
    - `low VAF + high AlleleDelta (+ low CramersV)` 仍較適合後段 artifact triage，不適合提前升級為 TP rescue 主規則
  - `HCC1395 5kHz paired`：
    - `gq>=10` 為負增益（`183 TP / 616 FP`, `F1 -0.004459`）
    - paired 最穩定仍是 `gq>=15`（`+0.000825`）與 `gq>=20`（`+0.000871`）
    - 甲基單獨與 `gq>=10 + 甲基` 目前都沒有形成正向 rescue
  - `HCC1395_DORADO paired`：
    - `gq>=15`（`97 TP / 88 FP`, `F1 +0.000502`）與 `gq>=20`（`51 TP / 7 FP`, `F1 +0.000749`）仍是主軸
    - 甲基單獨最佳僅 `hp_assign_rate>=0.99`：`68 TP / 87 FP`, `F1 +0.000041`
    - `gq>=10 + quality>=60 + hp_assign>=0.99` 雖然相對 `gq>=10` gate 有正增量（`delta F1 vs gate=+0.001961`），但整體仍低於 `gq>=20`
  - `ClairS-TO` 不一開始就全域放寬 `GQ` 是合理的：
    - `HCC1395 5kHz TO` 的 `caller_removed_fp` 共有 `2,720,008` 筆
    - 真正 `candidate_eligible` 只有 `298`
    - 被排除的 FP 中，`98.85%` 的主因是 `NonSomatic`
    - 這說明 `gq>=10` 只能在 candidate-eligible pool 中討論，不能直接全域放寬
  - `PairwiseMedianDist` 呈現區間性與樣本依賴：
    - `5kHz TO` 在 `[0.15,0.25)` 為正向 support，但 `>=0.25` 不再更好
    - `DORADO paired` 方向相反，較低 pairwise 才較合理
  - `HCC1395_DORADO TO candidate-specific InterSubMod` 已補齊，且已納入同一套四資料集矩陣：
    - candidate-specific coverage：`247/247 lost_tp`、`94/94 removed_fp`
    - `gq>=10`：`40 TP / 11 FP`，`F1 +0.000540`
    - `Quality_Score>=60 only`：`208 TP / 77 FP`，`F1 +0.002620`
    - `gq>=10 + Quality>=60`：相對 gate `delta F1=-0.000064`
    - `gq>=10 + Pairwise<=0.20`：`30 TP / 9 FP`，相對 gate `delta F1=-0.000142`
    - `gq>=10 + Strong/Subclone`：`12 TP / 2 FP`，相對 gate `delta F1=-0.000366`
    - 結論：TO 下跨 `5kHz` 與 `DORADO` 都支持 `caller-first + methylation-support`
    - 甲基不是全負效果，但目前仍未超過同一 caller gate 的 `GQ only`
    - `PairwiseMedianDist` 方向未跨樣本穩定：
      - `5kHz TO` 偏高 pairwise 為正向 support
      - `DORADO TO` 偏低 pairwise（`<=0.20`）較合理
      - 因此不可升級成所有 TO 樣本的全域規則
  - **總結**：
    - 第一層：`caller-first`
    - 第二層：`methylation-support`
    - 獨立旁路：`artifact triage`
    - 目前仍**不支持**把 `Pairwise>=0.20` 或 `gq>=10 + 高 pairwise` 升級為跨樣本規則
- **主要文件**:
  - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_GQ與甲基rescue系統性驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_GQ與甲基rescue系統性驗證_01.md)
  - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_HCC1395_DORADO_TO_candidate_specific甲基rescue驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_HCC1395_DORADO_TO_candidate_specific甲基rescue驗證_01.md)
  - [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix/analysis_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix/analysis_summary.md)
  - [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix/gq_threshold_sweep.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix/gq_threshold_sweep.tsv)
  - [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix/gq_plus_methylation_rule_sweep.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix/gq_plus_methylation_rule_sweep.tsv)
  - [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix_with_dorado_to/analysis_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix_with_dorado_to/analysis_summary.md)
- **建議後續**:
  - 在 `5kHz TO` 上對 `Quality>=60 / Pairwise>=0.15 / agreement_positive` 的 top TP/FP 補做 diagnostics，確認是否對應不同 read-level 現象
  - 在 `DORADO TO` 上對 `Quality>=60 / Pairwise<=0.20 / hp_assign>=0.99` 的 top TP/FP 補做 diagnostics，確認 support 規則是否對應不同 read-level 現象
  - 若要進流程層，優先新增 annotation，而不是直接改成 hard filter / hard keep

### 30. TO support 特徵的 read-level diagnostics（2026-03-11）

- **期間**: 2026-03-11
- **狀態**: ✅ 第一輪完成
- **關鍵結論**:
  - 已新增批次 diagnostics 腳本：
    - [/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_to_support_feature_diagnostics.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_to_support_feature_diagnostics.py)
  - 本輪直接對 `HCC1395 5kHz TO` 與 `HCC1395_DORADO TO` 的代表性 `TP/FP` 補做：
    - matrix heatmap
    - distance heatmap
    - `samtools depth/mpileup/HP tag` snapshot
  - 代表性位點數：
    - `5kHz TO`: `18`（`caller_lost_tp=9`, `caller_removed_fp=9`）
    - `DORADO TO`: `16`（`caller_lost_tp=8`, `caller_removed_fp=8`）
  - `Quality_Score`：
    - 在兩個 TO 資料集中都不是負訊號
    - 但高 `Quality_Score` 的 `FP` 也很常見
    - `DORADO TO` 中，高 `Quality_Score` 的 `TP` 通常 alt fraction 中等、HP balance 較正常；高 `Quality_Score` 的 `FP` 則常有極高 alt fraction 與極端 haplotype skew
    - **定位**：`soft support / ranking / annotation`
  - `PairwiseMedianDist`：
    - `5kHz TO` 的代表性 support 偏向高 pairwise
    - `DORADO TO` 的代表性 support 偏向低 pairwise
    - 這種方向差異在 read-level diagnostics 中仍存在，因此不能升級成全域固定閾值
    - **定位**：`dataset-aware annotation`，不是跨樣本硬規則
  - `hp_assign_rate`：
    - 高 `hp_assign_rate` 的 `TP` 與 `FP` 都很多
    - 多個 `FP` 同時呈現 `hp_assign_rate=1.0 + 極端 alt fraction + 極端 haplotype skew`
    - **定位**：`phase/QC annotation`，不應作為 truth-level support 主規則
  - 綜合結論：
    - 第一層：`caller-first`
    - 第二層：`Quality_Score` 等甲基 support
    - annotation / QC：`PairwiseMedianDist`、`hp_assign_rate`、`NA HP fraction`、`haplotype skew`
- **主要文件**:
  - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md)
  - [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/selected_regions.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/selected_regions.tsv)
  - [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/diagnostic_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/diagnostic_summary.tsv)
- **建議後續**:
  - 對 `Quality_Score` 做更細的區間分析，確認是否存在比 `>=60` 更穩的 support 帶
  - 將 `hp_assign_rate` 從 support 候選正式降階為 `phase/QC annotation`
  - 若要進流程層，優先新增 annotation，不先加 hard keep / hard filter

### 31. HCC1395 四象限甲基 rescue 整合矩陣（2026-03-11）

- **期間**: 2026-03-11
- **狀態**: ✅ 第一輪完成
- **關鍵結論**:
  - 已新增整合腳本：
    - [/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_hcc1395_cross_quadrant_matrix.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_hcc1395_cross_quadrant_matrix.py)
  - 已將 `HCC1395 5kHz / HCC1395_DORADO × paired / TO` 的 benchmark、candidate-specific coverage、特徵單獨/組合、交集/正交性與 TO snapshot scope 差異放進同一個整合框架
  - 新整合輸出目錄：
    - [/bip8_disk/liaoyoyo2001/InterSubMod_out/output/research_rounds/20260311_hcc1395_cross_quadrant_matrix](/bip8_disk/liaoyoyo2001/InterSubMod_out/output/research_rounds/20260311_hcc1395_cross_quadrant_matrix)
  - 已新增的核心表：
    - `model_availability_matrix.tsv`
    - `comparability_matrix.tsv`
    - `layered_performance_matrix.tsv`
    - `candidate_pool_layer_matrix.tsv`
    - `feature_distribution_focus.tsv`
    - `selected_rule_matrix.tsv`
    - `feature_overlap_focus.tsv`
    - `snapshot_scope_bias_assessment.tsv`
    - `phase2_gap_tracker.tsv`
  - 目前最重要的高層結論：
    - `caller-first` 仍是四象限最穩定主規則
    - 甲基特徵不是全部負效果，但多半沒有超過同一 caller gate
    - `PairwiseMedianDist` 方向具有明顯 dataset-dependence，不可升級成全域規則
    - `5kHz TO full tagged BAM` 與 `DORADO TO subset tagged BAM` 的 read-level snapshot 不可做跨 dataset 絕對值硬比較
    - paired raw caller 已有 `pileup_filter.vcf` 與 `full_alignment_filter.vcf` 可供 phase 2 直接 benchmark；TO 目前仍只有 pileup-only availability
- **主要文件**:
  - [/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md)
  - [/bip8_disk/liaoyoyo2001/InterSubMod_out/output/research_rounds/20260311_hcc1395_cross_quadrant_matrix/final_integrated_report.md](/bip8_disk/liaoyoyo2001/InterSubMod_out/output/research_rounds/20260311_hcc1395_cross_quadrant_matrix/final_integrated_report.md)
- **建議後續**:
  - 在 paired 補同規格 support-feature read-level diagnostics
  - 做 paired `pileup_filter / full_alignment_filter / final output` 的 direct benchmark
  - 若空間允許，再評估是否補 `DORADO TO full baseline InterSubMod`

### 32. Phase 2：paired raw pileup/full 與 finer feature interval / orthogonality（2026-03-11）

- **期間**: 2026-03-11
- **狀態**: ✅ 第一輪完成
- **關鍵結論**:
  - 已新增 phase 2 腳本：
    - [/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase2_paired_model_feature_analysis.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase2_paired_model_feature_analysis.py)
  - paired raw `pileup_filter / full_alignment_filter / merged_output` 的 direct benchmark 已完成：
    - `HCC1395 5kHz paired`
      - raw `pileup_filter`：`TP=30595 / FP=7450 / FN=8852 / F1=0.7896`
      - raw `full_alignment_filter`：`TP=30759 / FP=12684 / FN=8688 / F1=0.7422`
      - merged output：`TP=29865 / FP=1430 / FN=9582 / F1=0.8443`
      - **結論**：raw full 雖多 `894 TP`，但多 `11254 FP`；raw recall 空間真實存在，但不值得直接放寬成最終輸出
    - `HCC1395_DORADO paired`
      - raw `pileup_filter`：`TP=30872 / FP=2232 / FN=8575 / F1=0.8510`
      - raw `full_alignment_filter`：`TP=31332 / FP=2559 / FN=8115 / F1=0.8545`
      - merged output：`TP=29986 / FP=588 / FN=9461 / F1=0.8565`
      - **結論**：raw full 優於 raw pileup，但仍低於 merged output；paired 最終仍需 merge + downstream precision control
  - finer feature interval 已確認：
    - `5kHz TO`
      - `GQ` 最佳 support 帶集中在 `[18,20)`
      - `Quality_Score` 最佳帶集中在 `[55,60)`
      - `PairwiseMedianDist` 正向帶集中在 `[0.18,0.20)`
    - `DORADO TO`
      - `GQ` 最佳帶集中在 `[20,25)`
      - `Quality_Score` 同樣偏 `[55,60)`
      - `PairwiseMedianDist` 最佳帶改為 `[0.12,0.15)`
    - **結論**：`Quality_Score` 較像穩定 support，`PairwiseMedianDist` 仍有強烈 dataset-dependence
  - orthogonality 已確認：
    - paired 目前沒有清楚的甲基正交補強，仍由 `GQ` 主導
    - `5kHz TO` 中 `GQ / Pairwise / hp_assign / agreement / Strong` 間有小幅正交增益
    - `DORADO TO` 也有弱正交增益，但幅度更小
    - **結論**：甲基特徵在 TO 下仍屬第二層 support，不足以取代 caller-first
- **主要文件**:
  - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_phase2_paired_raw與feature_interval_orthogonality分析_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_phase2_paired_raw與feature_interval_orthogonality分析_01.md)
  - [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/phase2_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/phase2_summary.md)
- **建議後續**:
  - 若要進一步擴 paired phase 2，優先補 paired support-feature 的 read-level diagnostics，而不是直接擴 hard rule
  - 若要把 TO support 再往前推，下一步應該做 annotation / ranking layer，而不是直接把區間規則寫成全域 filter
  - 若要再深化正交性，應先排除 coverage ceiling，再測 `GQ + Pairwise + Quality + Strong` 的 constrained union

### 33. Phase 2 finer interval 回接 annotation layer（2026-03-11）

- **期間**: 2026-03-11
- **狀態**: ✅ 第一輪完成
- **關鍵結論**:
  - 已新增 annotation builder：
    - [/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_phase2_annotation_layer.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_phase2_annotation_layer.py)
  - 已將 phase 2 的 top-bin interval 正式回接到四個 dataset 的 analysis-layer annotation：
    - `HCC1395 5kHz TO`
    - `HCC1395 5kHz paired`
    - `HCC1395_DORADO TO`
    - `HCC1395_DORADO paired`
  - annotation layer 已明確區分三類訊號：
    - `support`
    - `qc`
    - `artifact`
  - `Quality_Score` 已可穩定升級為 `support annotation`
  - `PairwiseMedianDist` 已可升級為 `dataset-aware support annotation`
  - `hp_assign_rate` 已明確降階為 `phase/QC annotation`
  - policy evaluation 已確認：
    - `HCC1395 5kHz TO`：最佳仍是 `caller_primary_only (gq>=10)`，annotation policy 全數略低於 caller-primary
    - `HCC1395 5kHz paired`：最佳仍是 `caller_primary_only (gq>=15)`，paired 目前沒有可升級的甲基 support
    - `HCC1395_DORADO TO`：最佳仍是 `caller_primary_only (gq>=10)`，annotation 幾乎追平但未超過
    - `HCC1395_DORADO paired`：最佳反而是更嚴格的 `caller_strict_only (gq>=20)`
  - **總結**：annotation layer 有價值，但現階段最合理的用途是排序、註記與解讀，不是直接改 hard keep / hard filter
- **主要文件**:
  - [/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_phase2_finer_interval回接annotation_layer驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_phase2_finer_interval回接annotation_layer驗證_01.md)
  - [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotation_layer_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotation_layer_summary.md)
  - [/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md)
- **建議後續**:
  - 先把 annotation 欄位接到整合矩陣、round summary 與 diagnostics dashboard
  - 若要再深化，優先做 annotation score / ranking，而不是直接改 pipeline 預設規則
  - paired 若要再往前走，先補 candidate-specific coverage ceiling，再決定是否要提升 paired support annotation 的權重

### 34. TO 雙模型可得性、snapshot scope 與 Pool B FN 收尾盤點（2026-03-12）

- **期間**: 2026-03-12
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - `TO pileup vs full model` 目前不是單純少一張 benchmark 表，而是現有 source tree 中尚未確認有 `full model` 產物
  - `v0_3_0` 與 `ss_v0_3_0` 目前較像不同 `pileup model family`，不應直接寫成 `pileup vs full model`
  - `5kHz TO` 的 read-level diagnostics 來自 `full tagged BAM`，`DORADO TO` 來自 `candidate-window subset tagged BAM`，因此只能做方向比較，不能做跨 dataset 絕對值硬比
  - `Pool B FN = 840` 與 caller-side rescue 空間已成立，但尚未接回四象限整合、phase 2 annotation 與主線週報 closing summary
- **主要文件**:
  - `docs/plans/2026/03/20260312_未完項closing_plan_01.md`
  - `docs/experiments/in_progress/2026/03/20260312_TO雙模型與snapshot_scope盤點_01.md`
- **建議後續**:
  - 先完成 `TO dual-model availability closeout`
  - 再對 `5kHz TO` 做 `same-scope subset snapshot`
  - 最後將 `Pool B FN` 正式接回四象限整合與 annotation 層

### 35. TO 雙模型可得性 closeout（2026-03-12）

- **期間**: 2026-03-12
- **狀態**: ✅ 已完成
- **關鍵結論**:
  - `HCC1395 5kHz TO` 與 `HCC1395_DORADO TO` 在目前 source tree 中尚未確認存在可直接做 benchmark 的 TO full-model 產物
  - 現有 `v0_3_0` 與 `ss_v0_3_0` 應視為不同 `pileup model family`，不應再描述成 `pileup vs full model`
  - 目前 TO caller 模型結論應統一標示為：`pileup-route only`
  - paired raw caller 的 `pileup/full` 對照已存在，但不能直接外推到 TO
- **主要文件**:
  - `docs/reports/validated/2026/03/20260312_TO雙模型可得性closeout_01.md`
  - `docs/experiments/in_progress/2026/03/20260312_TO雙模型與snapshot_scope盤點_01.md`
- **建議後續**:
  - 下一步直接進入 `snapshot scope same-scope control`
  - 若未來需要 TO 雙模型 benchmark，先補 full-model source 或重跑，不要直接假設資料已存在

### 36. TO snapshot scope same-scope control（2026-03-12）

- **期間**: 2026-03-12
- **狀態**: ✅ 已完成
- **關鍵結論**:
  - 已對 `HCC1395 5kHz TO` 的相同 18 個代表性 region，完成 `full tagged BAM` vs `candidate-window subset tagged BAM` 的 same-scope control
  - `18/18` 個 region 在 region-level comparison 中都完全一致
  - `read_count`、`target_depth`、`target_alt_fraction`、`na_hp_fraction`、`collapsed_hp_balance_delta`、`top_hp_tag` 等所有主要 snapshot 指標都為 `max_abs_delta = 0.0`
  - 這代表 `subset snapshot` 在同 dataset / 同 region 控制下，不會扭曲目前使用的 read-level diagnostics 指標
  - 先前 `5kHz TO full` vs `DORADO TO subset` 的限制，現在應更精確改寫為：
    - `subset` 技術本身未觀察到偏移
    - 仍不可做跨 dataset 絕對值硬比的主因是 dataset/platform 差異與 coverage ceiling，不是 subset 萃取本身
- **主要文件**:
  - `docs/reports/validated/2026/03/20260312_TO_snapshot_scope_same_scope_control_01.md`
  - `docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md`
  - `InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/same_scope_snapshot_bias_summary.md`
  - `InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/same_scope_metric_summary.tsv`
- **建議後續**:
  - 可將 `subset snapshot` 正式升級為 `safe_for_within_dataset_ranking`
  - 保留 `unsafe_for_cross_dataset_absolute_comparison` 的限制，但理由改成 dataset/platform 差異，而不是 subset 技術本身
  - 這條主線的下一步已在 **實驗 37** 完成 `Pool B FN integration closeout`

### 37. Pool B FN integration closeout（2026-03-12）

- **期間**: 2026-03-12
- **狀態**: ✅ 已完成
- **關鍵結論**:
  - `Pool B FN = 840` 再次證明 `ClairS-TO` 在 `HCC1395 5kHz tumor-only` 的 non-PASS 區存在明確 caller-side rescue 空間
  - 原始 caller-only 的 authoritative 最佳規則是 `no_varcluster`：`428 TP / 390 FP`, `F1 +0.003391`
  - 若追求較乾淨、較可實務採用的規則，則是 `no_varcluster_and_gq15`：`115 TP / 45 FP`, `precision=71.9%`, `F1 +0.001452`
  - 補跑甲基後的最佳 combined 規則是 `gq15_and_allele_delta_low`：`431 TP / 473 FP`, `F1 +0.002702`
  - 這個最佳 combined 只比相近 caller gate `gq15_and_af015` 多 `+0.000386`，但仍低於原始 caller-only 最佳 `no_varcluster`
  - `class_strong_or_subclone` 與 `pairwise_ge_020` 在 Pool B 中都為負效益，表示甲基訊號不適合作為獨立主 rescue 規則
  - `AlleleDelta` 在 Pool B 與 TO kept-set / artifact triage 的方向不同：Pool B 中 `AlleleDelta < 0.15` 較像 somatic-rich 訊號；不能直接把 `low VAF + high AlleleDelta` 的 kept-set 規則全域套用到 Pool B
  - 本次也正式釐清證據層級：
    - caller-side FILTER/GQ 規則必須使用原始 `pool_b_caller_rescue_rules.tsv`
    - `with_methyl` 目錄的 `passified VCF` 只適合拿來看 coverage、joined features 與 combined rescue，不能再當作原始 FILTER 語意的 caller-only 對照
- **主要文件**:
  - `docs/reports/validated/2026/03/20260312_PoolB_FN_integration_closeout_01.md`
  - `docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md`
  - `docs/reports/validated/2026/03/assets/20260312_pool_b_fn_integration_summary_01.tsv`
- **建議後續**:
  - Pool B 主線已可關閉，不需再當成未解 blocker
  - 若未來重跑腳本，應將原始 `filter_tags` 額外保存在 sidecar TSV 或 INFO，避免 passified 後失去 FILTER 語意
  - 研究重心維持在 `caller-first + methylation-support + artifact triage` 的主線整合，而不是再把 Pool B 當獨立探索方向

### 38. TO FP 來源分解與 paired 對照分析（2026-03-22）

- **期間**: 2026-03-22
- **狀態**: ✅ 成功
- **關鍵結論**:
  - big7 observation workspace 已完成，確認 `HCC1395` 與 `HCC1395_DORADO` 的 TO raw FP universe 中，`PON / NonSomatic` caller veto 佔比都超過 `99.48%`
  - `caller_pon_filtered` 在本次 workspace 中 100% 帶顯式 `PoN_1..PoN_4` tag，可直接視為 PON layer
  - `LongPhase-TO` 在這兩個 TO pilot 上對 truth-scoped raw FP 的額外移除為 `0`；真正進入 TO final residual 的 FP 約占 raw universe 的 `0.4246%`
  - TO final residual FP 中，paired 可解比例都超過 `99%`，且主因幾乎全是 `paired_raw_absent`，不是 `normal_alt_support`
  - 目前沒有任何新的 TO-only 規則可超過 current big7 TO final F1；`HCC1395` discovery 與 `DORADO` validation 都是 `0 trigger`
- **主要文件**:
  - `docs/reports/validated/2026/03/20260322_TO_FP來源分解與paired對照分析_01.md`
  - `docs/reports/validated/2026/03/20260322_TO_FP來源分解摘要_01.md`
  - `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260322_to_fp_provenance_analysis/analysis_report.md`
- **建議後續**:
  - 後續 hard-case 診斷應直接聚焦 `paired_persistent_final_fp`：`HCC1395=87`、`HCC1395_DORADO=77`
  - `paired_raw_absent` 應視為 missing-normal candidate universe expansion，不宜直接簡化成 germline
  - `PairwiseMedianDist / hp_assign_rate / agreement_type` 目前較適合留在 annotation / ranking，而不是升級成 TO hard filter

### 39. TO residual FP 深入分析：raw_absent 細分、recurrence 與 persistent diagnostics（2026-03-23）

- **期間**: 2026-03-23
- **狀態**: ✅ 成功
- **關鍵結論**:
  - `raw_absent` 可以再細分，但最有辨識力的切法不是 methylation hard filter，而是 `shared exact hotspot` vs `platform-specific tail`
  - `HCC1395` 與 `HCC1395_DORADO` 的 `raw_absent` 各自有 `11430 / 11424` 個，其中 `11220` 個為 cross-platform exact shared hotspot，兩邊 shared 比例都超過 `98.16%`
  - shared `raw_absent` 的中位 `AF/GQ/QUAL` 都偏高，不像低品質尾巴；這批更像 strict blacklist / PON candidate inventory，而不是可用 `Quality_Score / PairwiseMedianDist / VerificationClass` 直接硬切的群體
  - 針對整個 `raw_absent`、`shared raw_absent` 與 `persistent` 的同層級 rule sweep 全部沒有正向規則
  - 只有 `HCC1395 raw_absent_platform_specific` 出現一條極小正訊號：`ad_alt<=3`，可移除 `11 FP / 7 TP`，`delta F1=+0.000032`；`DORADO` 不支持，不能升級為穩定規則
  - `paired_persistent_final_fp` 的 cross-platform shared hard set 為 `45` 個，而且不只是 `Noise`，仍混有 `Strong / Weak / Subclone`；比較適合進一步做 read-level / IGV / local-context 診斷
- **主要文件**:
  - `docs/reports/validated/2026/03/20260323_TO_residual_FP_deep_dive_01.md`
  - `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260323_to_residual_fp_deep_dive/analysis_report.md`
  - `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260323_to_residual_fp_deep_dive/stage_best_rules.tsv`
  - `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260323_to_residual_fp_deep_dive/cross_platform_shared_hotspots.tsv`
- **建議後續**:
  - 把 `11220` 個 shared `raw_absent` 當成 strict blacklist / PON candidate inventory，到獨立樣本驗證 recurrence，而不是直接當 deployable 規則
  - `paired_persistent_final_fp` 後續直接聚焦 `45` 個 cross-platform shared 位點做 read-level 與 local-context 深查
  - 目前 `VerificationClass / agreement_type / PairwiseMedianDist / Quality_Score` 仍應維持 annotation / ranking 角色，不宜升級成 TO hard filter

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
- `docs/provenance/ai_sessions/2026/02/20260209_InterSubMod專案完整分析報告_01.md`

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
- `docs/provenance/ai_sessions/2026/03/20260302_subsample混樣甲基化偏差_現況研究推論與驗證路線圖_01.md`
- `docs/provenance/ai_sessions/2026/02/20260213_HCC1395_subsample_purity_完整驗證報告_01.md`（Knowledge 路徑）

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
- `docs/provenance/ai_sessions/2026/02/20260209_InterSubMod專案完整分析報告_01.md`（§4.2, §7.1）
- `docs/experiments/validated/2026/01/20260107_F1_and_Data_Optimization_Report_01.md`

---

### 40. 研究方法與相關文獻突破方向全域分析（2026-03-24）

- **期間**: 2026-03-24
- **狀態**: ✅ 成功
- **關鍵結論**:
  - 系統性對照 12 篇核心外部研究（Long-Read POG、TRACERx、EVOFLUx、DeepSomatic、LongPhase-S、t-nanoEM、NanoMethPhase、MethSig、MethPhaser、MethylBERT、PRISM、MHB pan-cancer）
  - 確認 ISM 的「以 somatic SNV 為錨點 + read-level methylation clustering + 雙向顯著性驗證」組合在文獻中具新穎性，尚無完整對標
  - ISM 最強學術定位：不是 variant filter，而是「以長讀序 somatic variant 為錨點，整合甲基化 pattern、haplotype、copy-number 資訊，提供 read-level epigenetic context for variant interpretation」
  - 突破方向優先序（研究者確認）：E（ML read classification）→ A+D（Normal ref + CN/Purity-aware）→ B（Gene-level integration）→ C（CpG 功能分層）
  - 學界共識：(1) Purity/CNV 校正必要；(2) 長讀序多模態價值；(3) LOH-SNV 天然分離標籤；(4) Promoter methylation + LOH = second hit；(5) ONT 甲基化/錯誤耦合
  - 學界分歧/未驗證：(1) 局部甲基化 cluster 生物意義；(2) 甲基化驅動定義多樣；(3) clustering vs ML-based 比較缺乏；(4) 跨平台穩定性；(5) germline vs somatic ASM 區分
- **主要文件**:
  - `docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md`
- **建議後續**:
  - 依 E → A+D → B → C 順序逐步推進突破方向
  - Phase 1：ML read classification，在現有框架內改進 read-level pattern recognition
  - Phase 2：引入 normal 甲基化基線 + purity/LOH/CNV 校正
  - Phase 3：Gene-level evidence integration，多點關聯 → biallelic event narrative
  - Phase 4：CpG 功能分層，篩選重要甲基位點

---

### 41. Phase 1 ML read classification 研究啟動與資料缺口盤點（2026-03-25）

- **期間**: 2026-03-25
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - `Phase 1` 的第一個實際阻塞點不是模型架構，而是缺少統一的 `per-read training/export layer`
  - 現有資產已足夠支撐第一版資料層：
    - `rescue_joined_features.tsv`
    - `reads.tsv`
    - `methylation.csv`
    - `distance/.../matrix.csv`
  - discovery / validation dataset 已明確收斂為：
    - `HCC1395 5kHz paired/TO`
    - `HCC1395_DORADO paired/TO`
  - 第一版最合理做法不是直接宣稱 tumor/normal deconvolution，而是先做 `read × region` 單位的 feature export 與 baseline-aligned supervision
  - 已新增兩條已驗證可行的資料底座：
    - candidate-specific pilot exporter（`77` read rows）
    - baseline full manifest v1（四象限共 `141,014` regions）
  - 四象限 smoke test 已確認：
    - `8/8` selected regions 可成功 resolve 回 `reads.tsv + methylation.csv`
    - 共匯出 `1,701` 筆 read rows
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260325_Phase1_ML_read_classification研究啟動與資料缺口盤點_01.md`
- **建議後續**:
  - 新增 Phase 1 專用 exporter，把 region outputs 展平成 `read × region` 訓練表
  - 依 `phase1_training_manifest.tsv` 做 sharded read export，而不是一次全量展開
  - 補 `5kHz / DORADO` harmonization 與 split 規則，避免 cross-platform 誤比
  - 視需要補 `DORADO candidate-specific` 路徑 provenance 回收，但它已不是 Phase 1 baseline 的 blocking issue

### 42. Phase 1 label schema 與 harmonization 規格（2026-03-25）

- **期間**: 2026-03-25
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - `Phase 1A` 已正式收斂為 `within-tumor alt-support read classification`
  - `Phase 1B tumor-vs-normal` 目前不能直接開始，因為 sampled shard `14,165` 筆 reads 的 `is_tumor` 全為 `1`
  - code-level 檢查已補強這個結論：
    - `normal_reader` 會初始化
    - 但 `process_single_region(...)` 目前只吃 `tumor_reader`
  - manifest-driven shard export 已打通：
    - smoke4：`8/8` resolve、`1,701` reads
    - sample80：`80/80` resolve、`14,165` reads
  - exporter schema 已加入：
    - `dataset_role`
    - `harmonization_group`
    - `phase1a_task`
    - `phase1a_region_label`
    - `phase1a_read_label`
    - `phase1b_ready`
    - `phase1b_blocker`
  - `5kHz` 與 `DORADO` 不應先當同分布 joint pool；第一版應保留 `platform|mode` 的 harmonization 邊界
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260325_Phase1_label_schema與harmonization規格_01.md`
- **建議後續**:
  - 建立 `Phase 1A` 的 head-to-head baseline table
  - 規劃 normal-read export layer，作為 `Phase 1B` 的資料前置
  - 補 `5kHz / DORADO` group-aware normalization
  - 視需要補 `DORADO candidate-specific` provenance 回收

### 43. Phase 1A split manifest 建立（2026-03-25）

- **期間**: 2026-03-25
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - `Phase 1A` 的 region-level split pool 已固定成可重跑輸出
  - discovery pool：
    - `HCC1395 5kHz paired`
    - `HCC1395 5kHz TO`
    - 共 `70,463` regions
  - external validation pool：
    - `HCC1395_DORADO paired`
    - `HCC1395_DORADO TO`
    - 共 `70,551` regions
  - split manifest 已保留：
    - `dataset_role`
    - `harmonization_group`
    - `phase1a_task`
    - `split_role`
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260325_Phase1_label_schema與harmonization規格_01.md`
- **主要輸出**:
  - `output/synthesis/research_rounds/20260325_phase1a_split_manifest_v1`
- **建議後續**:
  - 依 split manifest 建立 `Phase 1A` 的 head-to-head baseline table
  - 將 sampled shard 與 split manifest 串接成訓練/驗證 protocol

### 44. Phase 1A read classifier benchmark round 1（2026-03-25）

- **期間**: 2026-03-25
- **狀態**: ⏳ 進行中
- **關鍵結論**:
  - `Phase 1A` 已從資料層準備推進到第一版可重跑 benchmark
  - sample80 與 sample200 的結論方向一致：
    - `majority_ref` 無效
    - `methyl_mean_threshold` 與 `logistic_methyl_only` 只有弱訊號
    - `context-rich` 模型明顯最佳
  - sample200 最穩定結果：
    - `logistic_context_only`
      - discovery holdout：`F1=0.7816`
      - external validation：`F1=0.9010`
    - `logistic_methyl_context`
      - discovery holdout：`F1=0.7882`
      - external validation：`F1=0.8908`
  - `methyl + context` 相對於 `context-only` 的增益尚未被穩定證明
  - 真正的困難點不是 cross-platform 崩潰，而是 `to-pure`
    - `ONT_5kHz|to-pure` discovery holdout：`F1=0.6962`
    - `ONT_Dorado|to-pure` external validation：`F1=0.8311`
    - 初版 failure diagnosis 顯示錯誤主要集中在：
      - `FP + Subclone`
      - `FP + Weak`
      - 一部分 `FP + Strong`
  - `paired-pure` 顯著較容易：
    - `ONT_5kHz|paired-pure` discovery holdout：`F1=0.9157`
    - `ONT_Dorado|paired-pure` external validation：`F1=0.9848`
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260325_Phase1A_read_classifier_benchmark_round1_01.md`
- **主要輸出**:
  - `output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample80_v1`
  - `output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample200_v1`
- **建議後續**:
  - 做更嚴格的 methyl incremental test
  - 對 `to-pure` 做更大 shard 的 failure diagnosis 穩定性檢查
  - 擴到更大的 shard 或分批全量 read export

### 45. Phase 1A incremental test 與 multi-bio validation round 2（2026-03-25）

- **期間**: 2026-03-25
- **狀態**: ✅ 已完成
- **關鍵結論**:
  - 已完成更大 shard 的 paired comparison：
    - mode-mixed `sample400`
    - paired-only multi-bio `sample637`
  - `sample400` 顯示 `methyl + context` **不是全域更好**：
    - discovery holdout：`delta F1=-0.0088`
    - external validation：`delta F1=-0.0206`
    - bootstrap `95% CI` 皆小於 `0`
    - 負增益主要由 `HCC1395_DORADO TO` 驅動：`delta F1=-0.0328`
  - paired-only multi-bio external validation 顯示 `methyl + context` 有小幅但穩定整體增益：
    - `context-only F1=0.8611`
    - `methyl+context F1=0.8722`
    - `delta F1=+0.0112`
    - bootstrap `95% CI=[+0.0044, +0.0188]`
  - 這個增益不是只在單一 sample 出現：
    - `HCC1954`: `delta F1=+0.0496`
    - `COLO829`: `delta F1=+0.0491`
    - `HCC1395_DORADO paired`: `delta F1=+0.0055`
  - 但仍存在 sample heterogeneity：
    - `H1437` 正向但未收斂
    - `HCC1937` 幾乎無差異
    - `H2009` 輕微負向且未收斂
  - bucket-level diagnostics 已指出異質性來源：
    - `sample400 to-pure` 退步集中在 `TP REF` 誤判與 `FP Subclone / Strong`
    - `H2009` 退步主要集中在 `FP + Weak`
    - `HCC1954 / COLO829` 增益集中在 `Subclone / Weak / Strong` buckets
  - 正確結論應改為：
    - `methyl + context` 在 `paired-pure multi-bio` 有穩定小幅增益
    - `methyl + context` 在 `to-pure` 目前沒有增益，甚至會退步
- **主要文件**:
  - `docs/experiments/in_progress/2026/03/20260325_Phase1A_incremental_test與multi_bio_validation_round2_01.md`
- **主要輸出**:
  - `output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample400`
  - `output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample400_v1`
  - `output/synthesis/research_rounds/20260325_phase1a_incremental_test_sample400_v1`
  - `output/synthesis/research_rounds/20260325_phase1_training_manifest_paired_multibio_v1`
  - `output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_paired_multibio_sample637`
  - `output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_paired_multibio_sample637_v1`
  - `output/synthesis/research_rounds/20260325_phase1a_incremental_test_paired_multibio_sample637_v1`
- **建議後續**:
  - 將 `paired-pure multi-bio` 升級成 `Phase 1A` 的主要驗證軸
  - `to-pure` 與 `paired-pure` 分開建模與評估
  - 對 `H2009 / HCC1937 / H1437` 做 dataset-specific diagnostics
  - 補 `platform/modcall-aware normalization`

### 46. LOH Round 1 cross-sample audit 啟動與決策紀錄（2026-03-27）

- **期間**: 2026-03-27
- **狀態**: ✅ 已承接至項目 47（本項保留啟動決策紀錄）
- **關鍵結論**:
  - `LOH` 研究已正式從通用規格進入 Round 1 啟動階段
  - Round 1 不採最小閉環，而是直接納入 7 個具甲基資料的樣本做 cross-sample paired audit：
    - `HCC1395`（`ONT_5kHz`，即使用者口中的 `HCC1395_HKU_5kHZ`）
    - `HCC1395_DORADO`
    - `HCC1937`
    - `HCC1954`
    - `H1437`
    - `H2009`
    - `COLO829`
  - `LOH-like` 暫不先固定為單一最終閾值；Round 1 先以完整分箱與連續分佈為主，`<0.1 / >0.9` 只保留為 legacy reference
  - `PS` 缺口暫不先改 C++；若同位點或案例需要，先從 phased VCF 補查 `PS`
  - same-locus `paired-vs-TO` compare 已提升為 Round 1 必做項目
  - case panel 固定只做 `4–8` 個代表案例，並要求覆蓋 paired/TO、TP/FP、交集/非交集與高 `HP0/HP3` case
- **主要文件**:
  - `docs/plans/2026/03/20260326_LOH盤點執行規格_01.md`
  - `docs/experiments/in_progress/2026/03/20260327_LOH_round1_cross_sample_audit啟動與決策紀錄_01.md`
- **建議後續**:
  - 建立 Round 1 `input_manifest.tsv`，把 paired canonical 與 TO pilot availability 一次盤清楚
  - 先產出全樣本 `hp_ratio_core` / `effective_hp_reads` / `hp0_ratio` / `hp3_ratio` 分佈
  - 優先完成 same-locus overlap / non-overlap 對照與代表案例選點

### 47. LOH Round 1 cross-sample audit 完成報告（2026-03-27）

- **期間**: 2026-03-27
- **狀態**: ✅ 完成
- **關鍵結論**:
  - Round 1 已完成 `7 paired + 7 TO` cross-sample LOH audit，總計：
    - `748,391` 個 region rows
    - `459,782` 個 same-locus paired-vs-TO rows
    - `5` 個代表案例
  - `HP_Ratio` 不能單獨解讀：
    - `effective_hp_reads=0` 的 `69,807` 個位點，`HP_Ratio` 全部都是 `0.5`
    - 後續分析必須改用 `hp_ratio_core + effective_hp_reads`
  - paired 與 TO 的 `LOH-like` 角色不同：
    - paired 整體 `FP/TP enrichment = 1.194×`，可作為 mild risk signal
    - TO 整體 `FP/TP enrichment =` ~~`0.912×`~~ → **`0.805×`（TP 富集）**（[修正 2026-03-30] HP integer tag bug，舊值無效。詳見 `20260330_TO_LOH_enrichment_post_hp_fix_01.md`）。TO LOH-like 在 TP 更常見（與 paired 方向相反），不能當 FP 主分離器
  - TO 的主要 discordance 不是 `both_fp`，而是大量 `TO-only FP`：
    - `to_only_fp = 126,865`（佔 same-locus 全體 `27.59%`）
    - `paired_only_fp = 1,912`
    - `both_fp = 1,517`
  - `LOH` 適合作為 `annotation / diagnostics / evidence panel` 第一層，但目前不適合作為 TO 的單獨主判斷器
- **主要文件**:
  - `docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md`
  - `docs/plans/2026/03/20260326_LOH盤點執行規格_01.md`
  - `docs/experiments/in_progress/2026/03/20260327_LOH_round1_cross_sample_audit啟動與決策紀錄_01.md`
- **主要輸出**:
  - `output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/`
  - `output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/loh_enrichment_by_sample_mode.tsv`
  - `output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/same_locus_compare.tsv`
  - `output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/cases/`
- **建議後續**:
  - 先做 paired 高 enrichment 樣本的 `LOH-like FP vs TP` 深化 case study：`HCC1954 / H2009 / H1437 / HCC1937`
  - TO 先做 `support-aware stratification`：`effective_hp_reads / hp0_ratio / PS`
  - 下一輪再把 same-locus compare 與 `CNV/LOH interval`、`gene-level second hit` 接起來

---

### 29. 系統性大規模數據觀察 O1-O10（2026-03-31）

- **期間**: 2026-03-31
- **狀態**: ✅ 9/10 完成（O9 FN 待 ISM 數據）
- **關鍵結論**:
  - 748,391 rows × 116 cols 的 master dataset 系統性觀察，產出 82 張圖表 + 53 份數據表
  - **Top 3 發現**:
    1. TO 模式下無任何單一特徵 AUC > 0.58（O1,O4,O5,O8 四次獨立確認）
    2. LOH penalty 是 TO QS 失效根因（O2,O3），移除後 AUC +0.045，但 ceiling ~0.55
    3. Paired/TO 是根本不同問題空間（O1,O5,O7）：FP rate 1% vs 31%，HP r=0.001
  - **其他重要發現**: GQ paired AUC=0.811 最強; 甲基化 AUC 全部 < 0.55; VerificationClass 無分辨力(V=0.023); TO FP rate 跨樣本差 8.6×; Read-level 甲基化對 ALT/REF 無用
  - **交叉驗證通過**: 零矛盾，AUC 數字精確一致，既有 4 個 Memory 結論全部確認
- **主要文件**:
  - 全域索引: `big7_disk_output/synthesis/observation_workspaces/OBSERVATION_INDEX.md`
  - 交叉驗證: `big7_disk_output/synthesis/observation_workspaces/20260401_cross_validation_report.md`
  - Level A 整合: `docs/reports/validated/2026/04/20260401_systematic_observation_O1_O10_summary_01.md`
  - 9 個 workspace: `big7_disk_output/synthesis/observation_workspaces/20260401_O{01-10}_*/`
  - 10 個腳本: `scripts/analysis/build_observation_O{01-10}_*.py` + `observation_common.py`
- **建議後續**:
  - P0: 移除 TO LOH penalty + Paired/TO 分離模型
  - P1: Phase 1A ML 特徵集確定（GQ+DP+甲基化 5+effective_hp_reads）
  - P2: 執行 O9 FN 觀察（需先跑 FN ISM）
  - P3: 補充未覆蓋欄位觀察（ClusterPermanova, HPFine, Significant）

### 30. O11 Within-Group Methylation Heterogeneity（2026-03-31）

- **期間**: 2026-03-31
- **狀態**: ❌ NEGATIVE（假說否決）
- **關鍵結論**:
  - 假說：Germline ASM（mQTL 驅動）在 HP 組內比 cancer ASM 更一致 → 可用 heterogeneity 區分
  - Epipolymorphism AUC=0.845 **完全由 n_reads confound 驅動**（TP 有 1.87× 更多 reads）
  - 殘差化 n_reads + num_cpgs 後：所有 6 個特徵 AUC = 0.509-0.578（近隨機）
  - Read-count matched bin [81-120] 中 epipolymorphism AUC=0.560
  - **教訓**: 未來所有 region-level TP/FP 分析必須控制 n_reads
- **主要文件**:
  - Workspace: `big7_disk_output/synthesis/observation_workspaces/20260331_O11_methylation_heterogeneity/`
  - 腳本: `scripts/analysis/build_observation_O11_methylation_heterogeneity.py`
  - 文獻調查: `docs/references/20260331_甲基化區分germline_somatic_variant文獻調查_01.md`
- **建議後續**: 此方向正式關閉。Within-group heterogeneity 不區分 germline/somatic at ISM resolution。

### 31. O12 TO LOH Methylation Scenario Discrimination（2026-03-31）

- **期間**: 2026-03-31
- **狀態**: ❌ NEGATIVE（三場景不可區分）
- **關鍵結論**:
  - 三場景假說（LOH→sSNV / sSNV→LOH / germline LOH）甲基化區分
  - 數據量：175,542 TO LOH rows × 22 features（9 existing + 13 novel）× 7 samples × 4-level confound
  - **AlleleDelta**: Raw AUC 最高 0.66 (HCC1937)，但完全是 caller_af confound；AF-bin 內 AUC 全 < 0.59
  - **L2 Collider Bias 發現**: CramersV/HeuristicScore/HPMergedSig/PassedGating L2 AUC 高達 0.80，但為 collider bias（AF-bin = 0.50）。機制：近常數特徵 residualized on AF 產生虛假信號
  - **Novel spatial features**（transition_count, cpg_autocorrelation, block_length_mean）：全 < 0.55
  - **Novel global features**（region_methyl_mean, region_low_methyl_fraction）：2-4/7 samples 弱信號，不穩定
  - Sample-first 分析 → 無任何特徵在 ≥5/7 samples 有信號（corrected AUC > 0.58）
- **主要文件**:
  - Workspace: `big7_disk_output/synthesis/observation_workspaces/20260331_O12_loh_methylation_scenarios/`
  - 腳本: `scripts/analysis/build_observation_O12_loh_methylation_scenarios.py`
  - 8 主圖 + 7 per-sample summary + 7 data TSV/JSON + observation_report.md
- **建議後續**:
  - TO 甲基化區分正式關閉 → caller features only (AF, GQ, DP)
  - 甲基化區分唯一路徑 = Phase 2 Direction A+D: Normal Methylation Reference
  - L2 collider bias 偵測方法（L2 vs L3 交叉驗證）應納入所有後續分析

### 32. O13 Cross-Region Methylation Correlation（2026-03-31 ~ 2026-04-01）

- **期間**: 2026-03-31（O13v1 快速驗證）→ 2026-04-01（O13v2 嚴格 confound control）
- **狀態**: ❌ NEGATIVE（跨區域甲基化 correlation 無法區分 TP/FP）
- **關鍵結論**:
  - **假說**: 同 subclone 的 TP-TP pairs 跨區域甲基化 correlation > FP-FP pairs
  - **O13v1 數據**: 8,400 cross-region SNV pairs（50kb 內）× 7 samples
  - **O13v1 表面結果**: FP-FP median r=0.739 > TP-TP 0.623（p=8.38e-13）— **方向反轉**
  - **Confound 診斷**: FP-FP shared reads 系統性偏高（median 36 vs 21）→ inflated correlation
  - **O13v2 嚴格驗證**:
    - 雙重分層（distance × n_shared，6 strata）：3 個 TP>FP + 3 個 FP>TP → **方向不一致**
    - OLS residualized（controlled for n_shared + distance）：delta=0.028, p=0.464, d=-0.071
    - Propensity matched（n=500）：delta=-0.037, p=0.403
    - Per-sample：resid 3/7 TP>FP — 無一致方向
  - **生物學解讀**:
    - HP concordance = 1.0 → shared reads 天然來自同 haplotype（uninformative）
    - ALT concordance：TP-TP=0.700 vs FP-FP=0.994（variant calling 特性，非甲基化信號）
    - Distance decay 曲線 TP/FP 幾乎重疊（r=0.96 at <5kb → 0.35-0.47 at 40-50kb）
  - ISM C++ 架構審查：**單區域設計，無原生跨區域分析能力**
- **主要文件**:
  - O13v1 Workspace: `big7_disk_output/synthesis/observation_workspaces/20260331_O13_cross_region_correlation/`
  - O13v2 Workspace: `big7_disk_output/synthesis/observation_workspaces/20260401_O13v2_cross_region_confound_control/`
  - 腳本: `scripts/analysis/build_observation_O13_cross_region_correlation.py`（v1）、`build_observation_O13v2_cross_region_confound_control.py`（v2）
  - O13v2: 7 張圖 + 7 TSV + evidence_chain_report.md + logic_review_report.md + final_verdict_report.md
  - **補充驗證**（審查員 3 Gaps 全部關閉）:
    - Gap 1: 20+ 未測試 ISM 欄位 → all pooled AUC < 0.565（`gap1_untested_ism_columns_auc.tsv`）
    - Gap 2: ML 組合 → additive < +0.013 over caller（`gap2_ml_combination_auc.tsv`）
    - Gap 3: CpG-SNP proxy → AUC=0.449 方向反轉（`gap3_cpg_snp_proxy.tsv`）
  - 置信度：9/10
- **建議後續**:
  - 多點甲基化關聯性方向**正式關閉**
  - 與 O11（within-group）、O12（LOH scenario）共同構成完整甲基化可辨別性證據鏈（13 個觀察一致 NEGATIVE）
  - TO 甲基化在 within-region、LOH-specific、cross-region 三維度 + ML 組合 + 未測試欄位 + CpG-SNP 均 NEGATIVE
  - 確認 ISM C++ 不需擴展多區域功能
  - **TO 策略**: caller features only (AF, GQ, DP)；Phase 2 Normal Reference 為最具前景突破路徑

### 33. SNV-甲基化關聯性定量分析（2026-04-01）

- **期間**: 2026-04-01
- **狀態**: ✅ POSITIVE（ASM 廣泛存在，跨 7 samples 驗證，但無法區分 TP/FP — 與 O1-O13 不矛盾）
- **關鍵結論**:
  - **問題轉換**: 從「甲基化能否區分 TP/FP」（分類）→「多少 SNV 位點存在 ASM」（流行率）
  - **5 種方法交叉驗證（全 7 samples, 4,017 read-level regions）**:
    - M1 ISM PERMANOVA (HPMergedSig): 748,391 regions → 32.2% 有顯著 ASM
    - M2 ALT/REF Wilcoxon (7 samples pooled): FP 46.4%, TP 35.9%
    - M3 Per-CpG Fisher: FP 55.1%, TP 55.5%（M3 偵測更多 but 不分 TP/FP）
    - M4/M5: 與 M2 高度一致（Jaccard 0.78-0.83）
  - **方向一致性**: **7/7 samples** FP ASM > TP ASM
  - **Germline ASM > Somatic ASM（pooled 7 samples）**:
    - Stringent (|Δ|>0.1): FP 26.4% vs TP 7.6%（3.5× 差距）
    - Median |delta|: FP 0.038 vs TP 0.022（1.7× pooled；HCC1395 alone = 4.0×）
  - **樣本異質性**: HCC1395 FP ASM=65.6%（outlier），其餘 24-54%；effect size ratio 0.8-4.0×
  - **ISM PERMANOVA 獨特價值**: 文獻 96% ASM 是 entropy imbalance（Jenkinson 2020）→ mean-based 方法完全漏掉 → ISM 距離矩陣方法能部分偵測 → 唯一同時捕捉 location + dispersion 信號
  - **HP vs ALT/REF**: HP 分群等於或優於 ALT/REF（更多 reads, 更平衡, 捕捉 cis-mQTL）
  - **與 O1-O13 調和**: ASM 存在≠能區分 → TP 也有 36% ASM → 重疊 → AUC<0.60
- **主要文件**:
  - Workspace: `big7_disk_output/synthesis/observation_workspaces/20260401_snv_methylation_association/`
  - 報告: `20260401_SNV甲基化關聯性定量分析報告_01.md`（v2.0）
  - 文獻: `docs/references/20260401_ASM偵測方法與預期比例文獻調查_01.md`（16 篇文獻）
  - 腳本: `scripts/analysis/build_snv_methylation_association_analysis.py`（multiprocessing 46 cores）
  - 數據: 5 TSV（per_region_multi_method 4,017 rows + per_sample_independent_proportions）+ 9 PNG
- **建議後續**:
  - ISM 的**正面定位**：read-level epigenetic context characterization（not variant filter）
  - **Entropy imbalance ASM 的定量驗證**：可直接用 ISM dispersion test 結果量化 ISM+/M2- 的 regions 比例
  - Normal ASM baseline（Phase 2 Direction A+D）：比較 tumor vs normal 的 ASM pattern 可分離 germline vs cancer ASM
  - 論文可引用 Jenkinson 2020 定位 ISM PERMANOVA 的理論優勢

### 34. TO Germline FP 鑑別研究（2026-04-01）

- **期間**: 2026-04-01
- **狀態**: ❌ NO-GO（多角度 VCF 特徵組合無法識別 TO 殘餘 germline FP）
- **關鍵結論**:
  - **問題**: 能否利用 VCF 中未充分使用的特徵（strand bias, CpG context, GT/H flag, AF, PoN pattern）組合識別殘餘 germline FP？
  - **架構**: 7 模組（G1-G7），48 張圖表，516K variants × 40+ 特徵 × 7 samples
  - **Go/No-Go 判定**: AUC ≥ 0.70 & FP removal ≥ 10% at TP loss < 2% & ≥ 5/7 一致
  - **G1 VCF 特徵提取**: 從 7 samples 的 ClairS-TO raw VCF 提取 40+ 特徵（strand counts, SB, GT, H, PoN, trinuc context）
  - **G2 Germline Signatures (AUC=0.507)**: H_flag 0.502, is_hom_alt 0.527, Combined 0.507 — 幾乎隨機
  - **G3 AF Architecture (AUC=0.566)**: AF reversal 確認（高 AF = 更多 FP），但 TP/FP AF 高度重疊
  - **G4 Strand Bias (AUC=0.502)**: 預期 negative — germline FP 是真 variant → 正常 SB
  - **G5 CpG Context (AUC=0.616)**: is_transition AUC=0.604（最佳單一特徵），CpG+Ti combined 0.616
  - **G6 Combination**: Rule-based 0.563; Decision Tree CV 0.635±0.048; **LR LOSO 0.638** (0.491-0.725)
  - **G7 Validation**: Bootstrap AUC 0.503 [0.501, 0.504]; 0/7 samples ≥ 0.60; **FP removal at ≤2% TP loss = 0%**
  - **根因**: 殘餘 germline FP（逃脫 PoN 的罕見 germline variants）在 VCF 特徵上與 somatic TP 本質相似
  - **唯一 7/7 一致特徵**: `is_transition`（FP transition rate 高 15-29% than TP），但 AUC 僅 0.604
  - **與 O1-O13 一致**: 加上 G1-G7 的 40+ VCF 特徵，TO 模式下**已測試 60+ 特徵**，全部 AUC < 0.64
- **主要文件**:
  - Workspace: `big7_disk_output/synthesis/observation_workspaces/20260401_germline_fp_identification/`
  - README: `20260401_germline_fp_identification/README.md`（含完整裁決與數據）
  - Scripts: `scripts/analysis/build_germline_fp_G[1-7]_*.py`（7 scripts）
  - 48 figures across G2-G7
- **建議後續**:
  - **TO germline FP post-hoc 識別方向正式關閉** — 60+ 特徵 × 7 samples 充分驗證
  - 正確方向 = FN rescue（TO FP provenance 研究結論一致）
  - PoN 改進（更完整 population database）是最有效的 germline 過濾改善路徑
  - 完成 O1-O13 + G1-G7 證據鏈後，TO 模式下**已窮盡 single-sample post-hoc 特徵探索**

### 35. Read-Level Haplotype & Methylation Germline FP 鑑別（2026-04-02）

- **期間**: 2026-04-01 ~ 2026-04-02
- **狀態**: ⚠️ CONDITIONAL NO-GO（LOSO AUC 首次 > 0.70 但安全約束 FAIL；根因=高純度細胞株）
- **關鍵結論**:
  - **問題**: G1-G7 的 60+ site-level 特徵無法區分 TO residual germline FP。能否透過 read-level 特徵突破？
  - **四代理文獻研究**: smrest（Simpson 2024）啟發 haplotype balance ratio；CPEL（Jenkinson 2020）解釋 O11 per-site entropy 失敗（96% ASM = entropy imbalance）；pure co-occurrence 假說理論缺陷確認；spatial proximity NO-GO
  - **Direction A: within_dom_alt_frac AUC=0.679（inverted）**: Dominant HP 中 ALT fraction。FP med=1.0（germline het = pure ALT within HP），TP med=0.875。Per-sample: HCC1954=0.851, HCC1937=0.816, HCC1395=0.725, COLO829=0.522
  - **Direction B: NME delta AUC=0.504（NO-GO）**: Cross-HP entropy difference 需 CPEL 級精細度
  - **Direction D: region_low_frac AUC=0.583（弱）**: Low-methylation read proportion
  - **LOSO Combination**: Mean AUC=**0.721**（首次 > 0.70，對比 G1-G7: 0.638）
  - **安全約束**: TP loss ≤ 2% → FP removal = **0%**（FAIL）
  - **methyl_low_fraction TO 遷移**: Paired 0.737 → TO 0.576（退化 -0.161）；3/7 good, 2/7 random, 1/7 reversed
  - **根因**: 高純度細胞株（80-95%+）→ ~40% somatic TP 也有 within_dom=1.0（因所有 reads 都是 tumor）
  - **Purity 依賴性**: 預期 purity 30-70% 時 AUC 可達 0.85+（somatic TP within_dom 會下降到 0.5-0.7）
  - **完整證據鏈**: O1-O13（22+ 甲基化 AUC<0.58）+ G1-G7（40+ VCF AUC<0.64）+ Read-level（4 特徵 LOSO 0.721 但 FP removal=0%）→ TO single-sample post-hoc **正式關閉**
- **主要文件**:
  - 報告: `docs/reports/validated/2026/04/20260402_read_level_germline_fp_research_report_01.md`
  - 文獻: `docs/references/manual/20260401_read_level_linkage_germline_literature.md`
  - 文獻: `docs/references/manual/20260401_methylation_read_clustering_literature.md`
- **建議後續**:
  - **TO single-sample post-hoc germline FP 識別正式關閉**
  - 低純度臨床樣本驗證 within_dom_alt_frac（purity 30-70%）
  - smrest 完整 Bayesian 實作（> simple ratio）
  - CPEL-level multi-CpG joint entropy（非 ISM binarized NME）
  - FN rescue + PoN 改進 = 正確方向（三輪研究一致結論）

### 36. O14 LOH Haplotag Bias 觀察（2026-04-02）

- **期間**: 2026-04-02
- **狀態**: ✅ 成功
- **關鍵結論**:
  - TP LOH 區域 99.5% HP_Ratio 極端偏向（<0.1 或 >0.9）
  - 100% region 內部 HP 方向一致（753 regions 全 concordance=1.0）
  - TO vs Paired 方向一致性 52.9%（≈隨機）— phasing scaffold 獨立，HP label 任意
  - HP0 reads 比例在 LOH 區域較低（reads 被強制分配到單一 haplotype）
- **主要文件**: `research/loh_investigation/scripts/observe_loh_haplotag_bias.py`
- **建議後續**: 此觀察確認 LOH 區域 HP 結構正確但不可跨模式比較；O15 接續量化全面影響

### 37. O15 LOH 區域 ISM 指標全面量化（2026-04-02）

- **期間**: 2026-04-02
- **狀態**: ✅ 成功（32 張圖 + 2 報告 + 6 TSV，Phase 1 + Phase 2 完成）
- **關鍵結論**:
  - **Phase 1（HCC1395 + SEQC2 4-class, 16 圖）**：
    - 141,014 rows，SEQC2 TP/FP/FN/TN LOH zones
    - LOH 內 TOP AUC: caller_gq (paired 0.922), caller_ad_ref (TO 0.784), LabelAllelePermanovaF (TO 0.757)
    - LOH 內甲基化指標失效: CramersV AUC ~0.53, PairwiseMedianDist ~0.54
    - caller_af 在 LOH 內 Cohen's d = -1.204
    - VerificationClass: LOH TP zone ~54% Noise
  - **Phase 2（7 samples × LOH.bed binary, 16 圖）**：
    - 748,391 rows，跨 7 樣本泛化驗證
    - Phase 1↔2 Concordance: R=0.993 (inside), R=0.997 (outside)
    - 跨樣本一致: caller_af, caller_ad_ref 在所有 sample LOH 內最強；ISM 甲基化全部 ~0.50
    - LOH coverage vs AUC: 負趨勢 — 覆蓋越大退化越嚴重
  - **總結**: LOH 區域內 ISM 甲基化分析全面失效（7/7 samples），只有 caller 特徵保留區分力
- **主要文件**:
  - Phase 1: `research/loh_investigation/reports/20260402_O15_phase1_loh_zone_metrics_hcc1395.md`
  - Phase 2: `research/loh_investigation/reports/20260402_O15_phase2_loh_zone_metrics_cross_sample.md`
  - 腳本: `scripts/analysis/build_observation_O15_loh_zone_metrics_hcc1395.py`, `build_observation_O15b_loh_zone_metrics_cross_sample.py`
- **建議後續**:
  - TO QS 重設計: LOH 區域內甲基化權重歸零，改用 caller 特徵
  - LOH.bed 可直接替代 SEQC2 做 zone 分類（R > 0.99 concordance）
  - 為 Phase 2 Normal Methylation Reference 提供 baseline：LOH 區域是 normal ref 最需要的區域

---

## 附錄：待驗證方向（尚未正式啟動）

| 方向 | 期望收益 | 難度 | 依據 |
|---|---|---|---|
| 5hmC 雙通道距離矩陣 | 可能更有亞克隆特異性 | 高 | ONT 5kHz 同時提供 C+h |
| ~~跨 Region 亞克隆一致性~~ | ~~真正的亞克隆應跨多 SNV 一致~~ | ~~高~~ | **已驗證 NEGATIVE**（O13/O13v2: cross-region correlation confounded by shared read count）|
| 機器學習組合特徵分類器 | 整合 15 個特徵的 ensemble | 中 | F1 研究發現特徵互補性 |
| PMD/ChromHMM Gating 啟用 | 降低甲基化背景噪聲 | 中 | 架構已有 `is_pmd` 欄位但未生效 |
| paired-only multi-bio 全量擴充與 sample-specific diagnostics | 確認 `+0.0112 F1` 是否在更大 region universe 仍穩定，並釐清 H2009/H1437/HCC1937 異質性 | 中 | round 2 已完成 sample637 驗證，但仍是 sampled shard |
