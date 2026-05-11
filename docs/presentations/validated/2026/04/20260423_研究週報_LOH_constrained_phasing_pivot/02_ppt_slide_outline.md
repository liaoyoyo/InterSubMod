<!--
建立時間: 2026-04-23
用途: 每 slide 的結構化大綱（title + bullets + figure + notes），供快速導航
對應 build_pptx.py 的 slide_NN_* 函數
-->

# PPT Slide Outline（26 slides）

格式 per slide：`title_zh | title_en | 主視覺 | 核心 bullet | speaker notes 重點`

---

## Act I · 定錨與方法學（Slides 1-9）

### S1 · Cover
- **zh**: NG=2 LOH-constrained phasing 發現與 TO-層論文主軸 pivot
- **en**: Mechanism correction: from methylation bimodality to phasing signatures in tumor-only sequencing
- **fig**: obs18_NG2_composition_proportion.png (封面預告)
- **bullets**: 涵蓋區間 2026-04-16~04-23 · 報告人廖子游 · 2026-04-23
- **notes**: 用戶 04-22 晚間一句提問觸發 C++ 回查，推翻 methylation bimodality 誤解

### S2 · Roadmap
- **zh**: 研究脈絡
- **en**: Context · This-week deliverables · Next-week priorities
- **layout**: 三欄（過去已關閉 / 本週 4 Thread / 下週 P0-P2）
- **notes**: 4 Thread 非獨立 — A 是 B 的單位、B 的 S4 引出 C、C 的 NG≥3=0 是 D 的獨立 negative control

### S3 · Recap 0421
- **zh**: 上週結論回顧
- **en**: Settled / Outstanding from prior week
- **layout**: 3 已解卡（LOH.bed Jaccard=1.0 / R1-Global 0.527 / F1-filter 終結）+ 2 blocker（stale-binary / flag=on NG≥3=0）
- **notes**: 0421 核心訊息「ISM 定位 variant filter → epigenetic characterization」已三重確認

### S4 · Motivation 四子問題
- **zh**: 本週提問
- **en**: Extending the stratification beyond LOH × AF × CN
- **main**: 主問「加入 HP / ISM 甲基 / ISM 關聯聚類能否進一步區分 TP/FP？」
- **subs**: Q1 CN 準確? / Q2 LOH×AF×CN 切到什麼程度? / Q3 HP 雜訊拖累? / Q4 NG 是什麼?
- **notes**: 四問有因果：Q1→Q2→Q3→Q4 環環相扣，是 PPT 三幕結構骨架

### S5 · CN 老問題
- **zh**: ACT I · Thread A · CN 方法學 · 舊做法共用 75× 的問題
- **en**: Old pipeline: hardcoded --expected-coverage 75.0 for all samples
- **layout**: 左右對照（舊做法紅框 6 bullets / 實際情況綠框 6 bullets 列 7 樣本 per-sample coverage）
- **notes**: CN tier 是 B 的測量單位；測量單位不準→apples-to-oranges

### S6 · KDE 雙 Pass 方法
- **zh**: ACT I · Thread A · 解方 KDE 雙 Pass
- **en**: Two-pass KDE calibration · bias −1.9% vs SEQC2 54× ground truth
- **fig**: fig3_two_pass_architecture.png
- **side**: 4 key numbers（KDE 53.0× / SEQC2 54× / bias −1.9% / stale +39%）+ commits（374fad4, 12d9b3e）
- **notes**: Pass 1 並行 / Mid KDE / Pass 2 重算；5 種 fallback 路徑 TSV audit column；pilot 無 fallback 觸發

### S7 · KDE 量化驗收
- **zh**: CovM median 0.880 → 1.245 (×1.415 恰為 75/53) · Category 大幅重分類
- **en**: Distribution shift and category reclassification
- **figs**: fig1_covm_distribution_shift.png + fig2_category_reclassification.png
- **banner**: Normal −5,718 / CNV_Gain +2,956 / High_Copy +2,710
- **notes**: ×1.415 = 75/53 等比例拉伸；影響跨樣本 CovM 絕對值；不影響 percentile-based filter (scale-invariant)

### S8 · LOH×AF×CN TP 熱圖
- **zh**: ACT I · Thread B · TP 雙極分佈
- **en**: TP rate heatmap reveals bipolar distribution
- **fig**: fig_v2_1_to_tp_heatmap.png
- **side**: 綠區 TP 88-96% (LOH×Extreme) / 綠區 TP 93-96% (canonical het) / 紅區 TP 47-60% (FP hotspot)
- **notes**: 雙極分佈 = biology-informed filter 可行性視覺證據；綠區對應 FACETS/Battenberg 模型

### S9 · S1-S7 Scheme + 伏筆
- **zh**: S3 Diploid Het TP 95.5% · S5 combo FP ↓99.37% · S4 無辨別力
- **en**: Biology-informed stratified filter schemes
- **fig**: fig_v2_3_filter_scheme_bar.png
- **side**: 5 列 scheme 數字表 + **S4 伏筆卡** 75% TP+76% FP → 引出 Thread C/D
- **notes**: S3 TP:FP fold 8.69×；S4 = everything bucket 含 75% 資料；NG≥3 邊際貢獻 <1pp

---

## Act II · 卡關與方案（Slides 10-17）

### S10 · S4 卡關 → 擴展特徵
- **zh**: 5 維座標無法切 S4 · 需加入 HP / ISM
- **en**: Bridging to Thread C: HP / ISM / clustering features
- **layout**: 左 S4 bucket (TP 67% + FP 33%) + 箭頭 + 右 3 候選特徵卡（HP / ISM 甲基 / ISM 關聯）
- **notes**: 用戶主軸原話「加入常用的 HP 狀況、ISM 甲基結果、ISM 關聯與聚類結論」；Beyond-AUC 耗盡≤0.58 的背景

### S11 · Self-Phasing 根因
- **zh**: ACT II · Thread C · TO self-phasing 污染 HP-derived 特徵
- **en**: Self-phasing contaminates HP bucket assignment
- **flow**: LongPhase-TO step03 → tumor_tagged BAM → HP1/HP1-1/HP2/HP2-1 → HPFineN 可能 artifact
- **nums**: 17.3:1 somatic bias / 62% ISM LOH 消失 / 0.617 HPFineN AUC
- **notes**: PON-only 已上游驗證；下游 ReadParser 層仍未解決；0418 marker AUC 可能是 artifact → 需 orthogonal test

### S12 · `--germline-hp-only` 方案
- **zh**: PON 當 germline 錨點 · somatic demote
- **en**: Design: anchor HP assignment to germline phase-set blocks
- **layout**: 三欄（flag=off 現行 / flag=on 新設計 / 預期效果）
- **commits**: 775036c / a61779c / 4dc2d73 / 2e2df22
- **notes**: 用戶原話「PON 當 tag 錨點，somatic 分 HP1-1/HP2-1/HP3」；ReadParser.cpp:123 注入；default=off

### S13 · Phase 1 機制驗證 Audit
- **zh**: Audit 獨立性 · somatic tag sum 兩 flag identical
- **en**: Audit independence confirms mechanism correctness
- **fig**: fig_c4_audit_independence.png
- **side**: 3 ✓ check（NHP_Somatic11/21/33 都 identical）
- **notes**: 排除「flag 誤 remove reads」替代解釋；下游任何觀察可安全歸因機制

### S14 · Phase 1 AUC Gate FAIL
- **zh**: 18 個 HP-related 特徵無一達 +0.02 Gate
- **en**: None of 18 HP-related features reach the +0.02 gate
- **fig**: fig_c1_auc_before_after.png
- **side**: 最大 +0.0084 / 最大 −0.0350 / AlleleDelta 0.6294 不動對照
- **notes**: Gate 設 +0.02；實測全 <+0.02；4 HP-derived 反向 ≤−0.025；HP-independent AlleleDelta 不動符合機制預期

### S15 · NGroups 崩潰到 NG=2
- **zh**: flag=on 下 NG≥3 regions = 0 · marker 訊號源質疑
- **en**: NG≥3 categories → 0 in both TP and FP
- **fig**: fig_c2_ngroups_distribution.png
- **side**: 關鍵數字 TP 22,732 + FP 8,148 → 0
- **notes**: 雙重意義 (1) 機制正確 (2) HPFineN marker 根基問題；Fisher odds ratio 0.913 反向 p=3.5e-3

### S16 · HP_Ratio shift + R3 判定
- **zh**: shift 方向正確但幅度小 · Plan R3 未觸發
- **en**: HP_Ratio shift minor; upstream hp_tag parsing is correct
- **fig**: fig_c3_hpratio_shift.png
- **side**: 最大 shift −0.023 / Plan 預期 0.836→0.55-0.65 不成立 / R3 判定結論
- **notes**: Plan 0.836 baseline 可能出自舊版 haplotag；R3 未觸發→bug 不在 hp_tag 解析→無需 LongPhase-TO 上游修復

### S17 · Tag 方案定位
- **zh**: 機制 ✓ · filter 懸掛 · Phase 2 判定前結論懸掛
- **en**: Mechanism validated; filter hypothesis pending Phase 2
- **layout**: 3 卡（H-C1 ✓ / H-C2 ⏸ / H-C3 ⏸）+ 底部 banner
- **notes**: 按用戶決策 A·b「延後定論」；避免「tag 可採用」肯定語；default=off 保留

---

## Act III · 機制揭露與 pivot（Slides 18-24）

### S18 · 偵探時刻
- **zh**: 用戶一句「NG=2 與甲基有關係嗎」觸發 C++ 回查
- **en**: A single user question triggered the mechanism reinterpretation
- **timeline**: V1 2026-04-22 15:56 誤解 / User 晚間提問 / V2 20:40 更正（4h44m 轉折）
- **notes**: 偵探時刻；feature name 直覺解讀陷阱；memory feedback_feature_name_vs_definition_rule 起源

### S19 · HPFineNGroups 真機制
- **zh**: 4-bucket occupancy count · 與 methylation 無直接關係
- **en**: C++ source code proves: pure phasing × variant-presence count
- **code**: LabelTest.cpp:265-302 的 hp_to_fine_labels() 4 行 if-else
- **side**: Stats.hpp:324 註 "Same haplotype, germline vs somatic"
- **banner**: 綠色「HPFineNGroups = 純 phasing × variant-presence · 不是 methylation bimodality」
- **notes**: 原始碼 = 證據鏈終點；設計者 FineN 指 fine granularity in phasing，非 methylation

### S20 · NG=2 四種組成
- **zh**: 4 bucket 中 2 個被 populate · 組成決定生物意義
- **en**: NG=2 composition determines biological meaning
- **layout**: 2×2 卡（same_HP1 / same_HP2 / cross_het / cross_het_inv）
- **banner**: same-hap = somatic signature · cross-het = germline leak 風險區
- **notes**: obs18 腳本邏輯基礎；物理推論 LOH→same-hap / 非 LOH→cross-het

### S21 · obs18 6 樣本驗證
- **zh**: Inner × NG=2 全部 ≥93% same-hap · 6/6 一致
- **en**: 6/6 TO samples: Inner NG=2 ≥93% same-hap
- **fig**: obs18_NG2_composition_proportion.png
- **side**: 6 樣本 Inner same-hap% 列表 (93.2 / 99.0 / 98.8 / 96.5 / 98.3 / 97.0) + median 97%
- **notes**: 主旨在 Inner 的跨樣本一致；HCC1395_DORADO 的 Outer 高 same-hap 是 basecaller 特殊 pattern 不影響主結論

### S22 · Inner vs Outer TP gap +0.37
- **zh**: 6/6 正向 gap · median +0.37 (+0.05~+0.52)
- **en**: TP rate gap: Inner same-hap vs Outer cross-het
- **fig**: obs18_NG2_composition_heatmap.png
- **side**: 6 樣本 Inner/Outer/Δ 對照表 + median +0.37 大字 + HCC1954 outlier 註
- **notes**: 0 反向；H2009 +0.05 是 baseline 飽和；HCC1954 outlier 下週 P1 專項；Wilcoxon 下週

### S23 · LOH-constrained phasing 機制 + pivot
- **zh**: 兩種物理場景 · 論文主軸從 methylation bimodality → phasing
- **en**: Two physical scenarios explain somatic vs germline leak
- **layout**: 雙面板（左綠 LOH Inner somatic必然 / 右紅 非 LOH Outer germline leak）
- **banner**: 新論文主軸（TO 層）: LOH-constrained phasing signatures ...
- **notes**: TO 層論文主軸正式 pivot 瞬間；paired 層 AF×NGroups 保留不撤回（按決策 C）

### S24 · Thread C 是 Thread D 獨立 negative control
- **zh**: flag=on NG≥3=0 同時是 filter FAIL 也是機制驗證
- **en**: flag=on NG≥3=0 as independent negative control for Thread D
- **layout**: 左 Thread C 流程 / 右 Thread D 流程 / 中央雙向箭頭
- **banner**: somatic HP demote → NG≥3 消失 · NG 是 phasing 而非 methylation · 機制 ⭐5
- **notes**: Thread C filter FAIL 在加上 Thread D 後重新解讀；4 假說全 POSITIVE；⭐5 TO 層 / paired 層需獨立驗證

---

## 結語（Slides 25-26）

### S25 · 結論總表
- **zh**: 6 項變動 · 3 新增 / 1 降級 / 1 加註 / 1 懸掛
- **en**: Conclusion tier table delta
- **table**: CL-LCP-001 ⭐5 / CL-S3-001 ⭐4 / CL-CN-KDE-001 ⭐4 / CL-016 ⭐4→⭐3 / CL-LAF-001 加註 / CL-HP-ONLY-001 ⭐3
- **notes**: ⭐5 = 三重證實；⭐4 = pilot PASS 待全量；⭐3 = 單樣本或懸掛

### S26 · 下週 P0/P1 + 里程碑
- **zh**: 4 項主要行動 · 2 週內定案
- **en**: Next-week P0/P1 priorities
- **timeline**: 04-23 今天 / 04-30 Paired+Archive / 05-07 master×兩 flag / 05-14 S4+Wakhan
- **actions**: P0 Paired obs18 / P0 Archive TO 6 樣本 / P1 HCC1954 outlier / P1 Wilcoxon
- **risk**: Paired gap 若未消失 → D 降 ⭐4 / master × 兩 flag 若 marker 暴跌 → marker 根基重建
- **notes**: Paired 層保留不撤回提供 fallback；2 週內兩 milestone；3 週後論文 Figure 3/4 重寫
