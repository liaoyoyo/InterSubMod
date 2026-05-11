<!--
建立時間: 2026-04-23
用途: 口頭報告稿（逐 slide 對應段落）— 約 2500 字，可直接朗讀
目標時長: 20-25 min（PI 提問彈性時間外）
備註: 詳細佐證與 tier 分析見 source_weekly_report.md；本稿專注連貫敘事
-->

# 週報 0423 口頭報告稿

## 開場（Slide 1-2 · ~2 min）

各位好，今天要報告 2026-04-16 到 04-23 的研究進展。這一週最大的發現是：**NG=2 的真實機制從 "methylation bimodality" 正式修正為 "LOH-constrained phasing signatures"**，這個修正重新定位了我們 TO 層的論文主軸。整週有四條研究軌：CN 方法學校準、LOH×AF×CN filter 框架、`--germline-hp-only` flag 驗證，以及最後的機制揭露。它們不是並排獨立的四個主題，而是環環相扣的偵探故事：測量單位先校準，才能談切分框架；框架卡關引出特徵擴展；特徵擴展的 negative 結果反而成為機制揭露的線索。

## 承接上週（Slide 3 · ~1 min）

上週 0421 週報收斂了三件事：LOH.bed 與 self-phasing 經 Jaccard=1.0 確認完全獨立；R1-Global subset→global 訊號塌陷確認 F pilot 的高 AUC 來自 overfit；ClairS-TO 內建 Verdict 雖內部校準優秀但 F1 增益為 0，paired F1-filter 路線正式放棄。留下兩個 blocker：stale-binary artifact 讓整個 master dataset 共用 75× baseline、以及 `--germline-hp-only` flag 下 NG≥3 完全消失的意外現象。這兩個 blocker 正是本週四條軌道的直接入口。

## 本週提問（Slide 4 · ~1 min）

這一週的出發點是一個明確問題：**從 LOH × AF × CN 切分 TP/FP 出發，加入 HP 狀況、ISM 甲基、ISM 關聯與聚類特徵，能否進一步區分？** 這問題拆成四個子問題：Q1 CN 怎麼切才準確、Q2 LOH×AF×CN 能切到什麼程度、Q3 HP 雜訊是否在拖累、Q4 NG 到底是什麼。接下來的報告就是按這四個問題順序鋪開。

## Act I · 測量單位（Slide 5-9 · ~5 min）

**Q1 CN 校準**。舊 pipeline 用 `--expected-coverage 75.0` 作為所有樣本共用的 diploid baseline，從這個 baseline 再乘倍數切 CN tier。這個預設對 HCC1395 偏高 39%，對其他樣本更不一致，HCC1395 實際約 54×、COLO829 約 47×、H2009 約 65×。問題是：CN tier 切錯會讓後面 LOH×AF×CN 的 heatmap 變成 apples-to-oranges 比較。解方是本週的 KDE 雙 Pass 校準：Pass 1 並行計算各 region metrics，Mid 階段用 histogram + Gaussian smooth + 2nd-deriv peak detection 估 per-sample 2× diploid peak，Pass 2 重算 Coverage_Multiple。HCC1395 pilot 估到 53.0×，對 SEQC2 neutral median 54× 偏差僅 −1.9%，從 ±39% 降到 <2%。Commits 374fad4 和 12d9b3e 另外加上 WARN + fallback 的 audit column，9+1 unit tests 全通過。修正後 CovM median 從 0.880 變 1.245，恰好是 75 除以 53 的比例，Coverage_Category 大幅重分類，符合 CN gain 原本就該被正確識別的物理預期。

**Q2 LOH × AF × CN 切分**。在這個校準後的 CN tier 上做 HCC1395 TO 40,115 sites 的 heatmap，baseline TP 71.1% 被清晰切成雙極分佈：綠區是 LOH × Extreme AF（TP 88-96%，對應 deletion/cnLOH 下的單 hap 存活），以及 None × Near-half × CN 約 2 或 3（TP 93-96%，對應 canonical het somatic）；紅區是 None × Intermediate × CN≥3（TP 47-60%，是 FP hotspot）。我定義了七個 filter scheme S1 到 S7。**S3 Diploid Het 單模組 TP rate 95.5%**，FP reduction 99.85%，TP:FP fold-improvement 8.69 倍；S5 combo 把 S1、S2、S3 聯集後 TP rate 仍 91.8%，FP reduction 99.37%。這些是本週在 TO filter 方向上最強的正向結果。但有一個關鍵 bucket 是 S4 —— LOH=None 加 Extreme AF —— 這個 bucket 含了 75% 全部 TP 和 76% 全部 FP，完全混雜，用現有的 5 維座標（LOH × AF × CN × mode × sample）無法切分。S4 的存在是 Act II 的直接導火線。

## Act II · HP 特徵擴展 + tag 方案（Slide 10-17 · ~7 min）

**Q3 HP 雜訊**。S4 需要加入新特徵判別。用戶主軸明確：嘗試加入常用的 HP 狀況、ISM 甲基結果、ISM 關聯與聚類結論。但在嘗試之前要先解決一個已知雜訊：TO 模式下 LongPhase-TO 在 somatic block 階段會把 somatic 來源的 HP tag 混進 HP1-1 和 HP2-1 分群，這可能污染所有 HP-derived 特徵。上週已上游驗證 PON-only phasing 讓 LOH.bed 不變、somatic bias 消除，但 ReadParser 層仍保留 LongPhase 的所有 HP tag。如果 HPFineNGroups 的 N=4 subclone marker 訊號是在這種污染下測出來的 AUC 0.617，那它可能是 artifact 而非真 biology。

本週的方案是用戶提的設計：`--germline-hp-only` flag，用 PON 當 germline 錨點，只信任 PS tag 屬於 germline SNV block 的 HP 分配，somatic 來源的 HP 被 demote。Commits 775036c、a61779c、4dc2d73、2e2df22 實作並驗證。default=off 不影響既有流程，flag=on 作為研究工具。

Phase 1 我們在 HCC1395 TO 全量 40,115 sites × 4 runs 上做驗證。**機制驗證完全通過**：Audit 獨立性檢查顯示 NHP_Somatic11/21/33 三個全基因 sum 在 TP+FP 合併下 flag=off 和 flag=on identical，排除「flag 實作錯誤誤 remove reads」的替代解釋。這是機制正確性的關鍵證據。

**但 filter AUC Gate 未過**。18 個 HP-related 特徵沒有一個達到預設的 +0.02 Gate，最大正向只有 HP1FamilyN 的 +0.0084；而 HPFineNGroups、NHP3、HPMergedDelta、HPFine_NGroups_CF 四個主要 HP-derived 特徵反而下降了 −0.025 以上。AlleleDelta（唯一 HP-independent 的中等訊號特徵，AUC 0.6294）完全不動，符合機制預期 —— 這驗證了「flag 只影響 HP 依賴訊號」。

Phase 1 有一個意外觀察：flag=on 下，TP 22,732 regions 加 FP 8,148 regions 的 HPFineNGroups N=3 和 N=4 類別完全消失，97% TP regions 和類似比例 FP regions 全部集中在 NG=2。這是 Act III 最大的線索。

HP_Ratio 的 shift 也比 plan 預期小很多，plan 假設從 0.836 降到 0.55-0.65，實測從 0.549 降到 0.529。這說明 V3-Fixed TO BAM 的 HP_Ratio 本就較低，plan 引用的 0.836 可能來自舊版 haplotag。重要的是 Plan R3 「修正後仍偏 0 則 bug 在 hp_tag 解析」的條件沒被觸發，修正定位為「解析正確」，無需升級為 LongPhase-TO 上游 C++ 修復工作。

依照用戶決策 A-(b)「延後定論」，tag 方案的結論是：**機制正確但目前 filter 無增益，Phase 2 全樣本判定前結論懸掛**。default=off 保留，flag=on 作為 characterization 工具，不進 Phase 2 全量重跑。HPFineNGroups 的 marker 使用必須標註 pipeline dependency —— 因為 HCC1395 TO ClairS-TO raw split 上（TP rate 0.6944）無法重現 memory 記錄的 master dataset 89.1%（Fisher odds ratio 0.913 反向 p=3.5e-3）。CL-016 的穩定度從 ⭐4 降為 ⭐3。

## Act III · 機制揭露（Slide 18-24 · ~6 min）

**Q4 NG 到底是什麼**。這是本週的偵探時刻，也是整個故事的高潮。起因是 2026-04-22 下午 15:56 我的分析版本 V1 把 HPFineNGroups 按字面意思解讀為「HP-Fine-N-Groups」—— fine-grained methylation N groups，並基於此擬議了 "Haplotype-loss-dependent methylation bimodality" 的論文主軸。當天晚間用戶問了一句「NG=2 與甲基有關係嗎」，我回查 C++ 原始碼，到 20:40 完成版本 V2 的機制更正。4 小時 44 分的關鍵轉折。

真機制在 `src/core/LabelTest.cpp:265-302` 的 `hp_to_fine_labels()` 函數：它純粹把 reads 按 hp_tag 字串分到四個 bucket —— HP1（haplotype 1 ref）、HP1-1（haplotype 1 somatic）、HP2、HP2-1。HPFineNGroups 就是這四個 bucket 中被 populate 的數量。程式碼沒有任何 methylation 計算；`include/core/Stats.hpp:324` 的設計者註 `"Same haplotype, germline vs somatic"` 也明確說明原意是 phasing 而非 methylation。Methylation 只在 HPFineF 和 HPFineP 的 PERMANOVA 中作為**品質檢驗**參與，NGroups 本身不依賴 methylation。

理解真機制後，NG=2 的四種組成就有了生物學分類：same_HP1 是 HP1+HP1-1（單 hap 內部的 ref 子族 + somatic 子族）、same_HP2 是對稱版本、cross_het 是 HP1+HP2-1（兩個 hap 一個全 ref 一個全 somatic，也就是 canonical het phasing）、cross_het_inv 是對稱版本。同_hap 系列代表單 haplotype 存在（LOH 必然）；cross-het 系列代表雙 haplotype 保留，這裡的關鍵物理限制是 —— 真 somatic het 和 germline het 都產生完全相同的 cross-het phasing pattern，caller 在 TO 模式下無法區分。

obs18 腳本跨 6 TO 樣本拆分 NG=2 組成。**Inner 區（LOH 內）的 NG=2 全部 93% 以上是 same-haplotype**：HCC1395 93.2%、HCC1395_DORADO 99.0%、HCC1937 98.8%、HCC1954 96.5%、H2009 98.3%、H1437 97.0%，median 97%，6 個樣本完全一致。這個數字是物理必然的直接觀察 —— LOH 區只有單 hap 存在，somatic SNV 發生後必然產生 same-hap 分裂。

TP rate gap（Inner same_HP1 TP − Outer cross_het TP）median +0.37，範圍 +0.05 到 +0.52，6/6 正向，0 反向。這個 gap 就是 TO 模式 germline leak 的物理量化 —— Outer cross-het 因為 germline-somatic phasing 不可分混入了 germline FP，TP rate 才降下來。HCC1954 的 Outer cross-het TP 0.08 是 outlier，下週 P1 要做 Potential_LOH 可靠性專項分析。

這裡我要做一個重要限定：**H-D2 「物理必然」指的是 phasing pattern 必為 same-hap**（符合 LOH 單 hap 物理限制），**而非 TP rate 必為高**。HCC1954 的 Inner same-hap TP rate 0.43 是低的，這反映 TP rate 仍受 Potential_LOH annotation 可靠性、ClairS-TO caller 行為、CNV 異質性等外部因素影響；但 phasing composition 本身 —— 93-99% same-hap —— 不受這些因素干擾，跨 6 樣本穩定。

這裡還有 Wave 1 review 指出的另一個限定：Thread C 的 flag=on NG≥3=0 觀察**直接驗證的是 H-D1**（NG 本身是 phasing 而非 methylation）—— 因為 somatic HP bucket 被 demote 後，4-bucket 的 occupancy 只能落在 {0, 1, 2}。**H-D4 的完整定義「flag=on 應消除 Inner/Outer gap」需要下週 P0 paired 對照**才能完整驗證。所以本週的結論穩定度是 ⭐5 在 H-D1/D2/D3 三重證實之下的 TO 層結論，H-D4 gap-disappearance test 待下週。

基於以上，TO 層論文主軸從 "Haplotype-loss-dependent methylation bimodality" 正式 pivot 為 **"LOH-constrained phasing signatures distinguish somatic from germline-like variants in tumor-only sequencing"**。新主軸的優勢：機制純 phasing 不需 methylation 實驗驗證；可直接從 LongPhase output + LOH bed 重現；跨 basecall version（DORADO）、跨 pipeline（LongPhase-TO）穩定；直接連接 FACETS、Battenberg、TITAN 既有文獻模型。依用戶決策 C，paired mode AF × NGroups POSITIVE 的舊結論**保留不撤回**，僅加註「需獨立 phasing-vs-methylation 驗證」作為下週 P0 延伸項目。

## 結語（Slide 25-26 · ~3 min）

結論總表有六項變動：新增三項 CL-LCP-001 ⭐5 TO-only（論文主軸 pivot 核心）、CL-S3-001 ⭐4（filter framework）、CL-CN-KDE-001 ⭐4（CN 方法學）；降級一項 CL-016 從 ⭐4 到 ⭐3（HPFineN marker pipeline dependency）；加註一項 CL-LAF-001（paired AF × NGroups 保留不撤回 + 加註）；新增一項懸掛 CL-HP-ONLY-001 ⭐3（Phase 1 機制 ✓ filter 懸掛）。

下週 P0/P1 優先行動四項：P0 HCC1395 paired normal-pilot 的 obs18 對照，驗證 H-D3/D4 在 germline-排除條件下 gap 是否消失，1-2 天；P0 Archive TO 6 樣本重跑 ISM，補齊 KDE + LOH.bed 欄位做跨樣本 S1-S7 scheme 驗證，約 10 小時 parallel；P1 HCC1954 outlier 的 Potential_LOH 可靠性專項分析，1 天；P1 Wilcoxon signed-rank 正式統計加 bootstrap CI，半天。

里程碑時間軸：今天定稿 0423 週報；下週 04-30 完成 Paired 對照和 Archive 重跑；再下週 05-07 Phase 2B master × 兩 flag 重跑決定 HPFineN marker 生物學根基；約 05-14 做 S4 二級判別和 Wakhan pilot，以及 ISM 論文 Figure 3 和 Figure 4 重寫。

主要風險兩條：Paired 對照若 gap 未消失，Thread D 穩定度降為 ⭐4 TO-only；master × 兩 flag 若 HPFineN marker TP rate 暴跌，marker 生物學根基需重建。後者的緩解是 paired 層保留不撤回的決策提供了 fallback 空間。

以上就是本週進展。歡迎提問。
