---
title: 限制、泛化性質疑、未來工作
status: in_progress
last_updated: 2026-04-22
---

# 06. 限制、泛化性質疑、未來工作

---

## 6.1 五大主要限制

### 6.1.1 單樣本 testbed（H4 定量證偽，2026-04-22 更新）

**問題**：主要結論（Top-17 Pareto-dominates S3）僅在 HCC1395 TO 一個樣本上成立。

**證據（2026-04-22 新增定量）**：
- COLO829 paired_full 的 top cells 座標結構異於 HCC1395 TO（見 §4.7 + `05_biology_interpretation.md` §5.5）
- 5 個 paired_full 樣本 FP 稀少（<1,400 total），無法獨立做 5D cube 分析
- **`07_figure_layers.md` §7.5.1 定量**：以 "cell FP ≤ 0.5 × baseline FP" 為跨樣本一致性判準，Top-17 中只有 **1/17（5.9%）** 達到 `n_samples_high ≥ 5` 門檻（LOH_Noise + Extreme + T3 + NG=2 + high_NR，Top-17 rank=2）
- **L1/L2/L3 目視證據**：paired_full 所有樣本 baseline 已 0.94-0.99，LOH 各類別差異 <0.02，**完全無 HCC1395_TO 的 LOH 梯度**（見 `07_figure_layers.md` §7.2-7.4）

**影響**：原「Top-17 是可跨樣本部署的 filter」主張被**定量上駁斥**。跨樣本共識 cells 走向「Extreme AF + 單純 NG=2」而非 Top-17 強調的「NG=3 + methylation 複雜度」。

**緩解**：
- `04_comparison_narrative.md` §4.9 已標 H4 未證
- `07_figure_layers.md` §7.7 提出的替代方向：LOSO 聚焦 **7-cell 跨樣本共識子集**（n_samples_high≥5 者），而非整個 Top-17

### 6.1.2 Empirical bins 無平滑 + 過擬合風險

**問題**：5D cube 的 900 cells 中活躍 ~300-400，Top-17 從排序取得，可能對該 sample 的 noise 過擬合。

**證據**：
- Top-17 最後一列是 `LOH_Weak + Intermediate + T2 + NG=2 + high` n=24，TP rate 91.7%；n=24 Wilson CI [72.7%, 97.9%] — 置信區間下緣遠低於 Top-17 cum 的 96.1%
- 單純換 random seed 或重取 sample 可能讓這個 cell 掉出 Top-17

**緩解選項**（未做）：
- Bootstrap 重取 500 次，只保留在 ≥400 次中都 rank ≤ 17 的 cells 作「穩定 Top-17」
- Leave-one-chromosome-out 驗證 Top-17 是否跨 chromosome 穩定（延伸 obs04 空間分析）

### 6.1.3 無真 F1（無 FN 資料）

**問題**：本研究的「fold-improvement」是 caller TP:FP ratio 的比值，不是 F1 改進。無 FN (caller-missed truth) 無法算真 precision/recall。

**影響**：若 caller 漏報一大部分 truth (低 sensitivity)，Top-17 在「真 truth 集」上的 precision 可能遠低於 96.1%。

**緩解**：需要 SEQC2 FN 資料重新評估。這在 `01_research_question.md` §1.5 已明列為非本研究範疇。

### 6.1.4 NG 訊號依賴 flag=off

**問題**：HPFineNGroups 維度在 `--germline-hp-only` flag=off 下包含 somatic HP tag 的人為分群。2026-04-21 Phase 1 驗證確認 flag=on 下 NG≥3 幾乎消失。

**影響**：
- Group A 和 Group D 的高純度 cells 依賴 NG=3 — 若 flag=on 重跑，這些 cells 可能不存在
- Top-17 的 4 個 Group 中 Group B (None + Near-half，NG=2 為主) 是唯一不依賴 NG≥3 的

**緩解**：本研究用 flag=off 的 KDE baseline（與 v2 一致），結論在該 baseline 下自洽。flag=on 重驗證為獨立 HP-only 研究（見 `project_readparser_germline_hp_only_phase1_negative.md`）。

### 6.1.5 S3 定義 T1 包含 hemizygous 可能性

**問題**：S3 定義 `cn_tier_F ∈ {T1, T2}`，T1 對應 CovM 0.65-0.99 — 理論可對應 CN=1 hemizygous。

**影響**：
- CN=1 的 "Near-half AF" 在概念上不應存在（hemizygous 只有一個 allele）
- 若 T1 實際 CN=1，則 T1 Near-half 的 variant 是 mapping artifact 或 contamination

**證據**：
- v2 §3 顯示 T1 vs T2 的 TP rate 差異 <2%（均 >94%）
- obs02 CovM 分佈顯示 T1 主要為 0.85-0.95（偏 CN=2 而非 CN=1）

**緩解**：T1 的實際 CovM 中位數偏向 CN=2 邊界，合併入 S3 不嚴重。嚴格版本可定義 S3' = `T2 only`；n 降至 171 (rank 5 cell)，TP 不變高。

---

## 6.2 生物學假說的不完整驗證

`05_biology_interpretation.md` 提出 4 個生物學 Group (A-D) 的詮釋，但**四個假說都未直接驗證**：

| 假說 | 現況 | 驗證需求 |
|------|-----|---------|
| Group A (NG=3) 是 methylation-confirmed somatic LOH | 依賴 flag=off NG | flag=on 重跑；read-level methylation IGV 視覺檢視 |
| Group C (LOH_Noise + high NR) 是真實弱 LOH | 推測 | 疊加 LOH.bed 連續性、檢視 allele ratio 直接驗證 |
| Group D (LOH_Weak + Int + T1) 是亞克隆訊號 | 與既有結論一致 | cluster analysis 識別 subclone |
| S4 無法 rescue 是 germline leak | AlleleDelta 證據 | paired mode 比對 normal read |

**這些假說若部分被駁斥**，不影響本研究的「5D Pareto envelope 存在」主結論，但會影響對 Top-17 的 biological interpretation。

---

## 6.3 與已關閉結論的一致性檢查

本研究結論與 memory 中既存結論的相容性：

| 既存結論 | 本研究對照 | 相容性 |
|---------|----------|-------|
| `project_O12_loh_methylation_scenarios.md`: AlleleDelta=AF confound | Top-17 中 AF 是主軸；但 **已透過 5D cube 全聯合繞過此 confound** | ✅ 相容 |
| `project_hpfinengroups_subclone_marker.md`: NG≥4 + NR≥80 TP=89.1% | Top-17 中 NG=4 未出現（可能 n 太小被稀釋）；NG=3 主導 | ⚠ 需核對 |
| `feedback_L2_collider_bias.md`: residualize on AF 會產生虛假信號 | 本研究**無 residualize**，直接 categorical；避開此陷阱 | ✅ 相容 |
| `feedback_pooled_ols_residualization_trap.md`: 必須 within-group OLS | 本研究全用 categorical，無 OLS | ✅ 不適用 |
| `project_germline_fp_identification_nogo.md`: TO Germline FP G1-G7 全 AUC<0.64 | 本研究繞過 Germline FP 直接判斷方向，改看 cell-level purity | ✅ 不衝突 |

---

## 6.4 泛化性質疑（reviewer 會怎麼批）

預想 reviewer 會問的 5 個 hard questions：

### Q1: 單樣本結論能發表嗎？

**回答**：H4 未證，確實不足以作變異過濾器主張。但**方法學貢獻**（5D cube + Pareto envelope 框架）仍獨立；Top-17 作為 **HCC1395 TO 的 characterization marker** 有參考價值。

### Q2: Top-17 是否過擬合 HCC1395 TO？

**回答**：§6.1.2 已承認此風險。需 bootstrap 或 cross-validation 驗證，本研究未做。

### Q3: 換個 CN 分箱策略結論是否成立？

**回答**：v2 §3 測試了 A/B/C/F 四種 CN 策略，S3/S5 的 TP rate 相對穩定（±2%）；5D cube 的 Top-17 未重測其他策略下的結果。留待後續驗證。

### Q4: 為什麼只做 HCC1395 TO 而不做其他 TO？

**回答**：只有 HCC1395 有 TO mode 的 ISM 特徵輸出；其他 6 樣本 TO mode 需 C++ pipeline rerun（out-of-scope）。

### Q5: 5D cube 沒 smoothing，Top-17 裡 n=22 的 cell 能信嗎？

**回答**：不能單獨信。Top-17 的**集合**（n=1,099）Wilson CI [95.0%, 97.0%]（假設 binomial）仍可信；**單 cell 層級**確實有統計不確定性。

---

## 6.5 未來工作（優先順序）

### 6.5.1 Tier 1（必做，驗證 H4）

1. **LOSO 跨樣本驗證** — *2026-04-22 更新：已有初步定量替代驗證*
   - **初步結果（見 `07_figure_layers.md` §7.5）**：基於「cell FP ≤ 0.5 × baseline FP」判準，直接計算 Top-17 在 8 sample-mode 下的 high 次數
   - 結果：Top-17 僅 1/17 達 n_samples_high≥5 — H4 **目視 + 定量上已大致駁斥**
   - **正式 LOSO 仍建議做**：改用 Top-17 座標在 COLO829 pf 的精確 purity 計算（而非 binary halving）
   - 若 COLO829 pf 的 Top-17 purity < 94%（baseline），H4 明確 NO-GO
   - 預期時間：30 min（一個腳本）

2. **跨樣本共識 top cells**
   - 取 HCC1395 TO Top-20 ∩ COLO829 pf Top-20
   - 若交集 ≥5 cells 且 purity 在兩樣本都 ≥95% → 「共識 white-list」
   - 若交集 ≤2 cells → 「各樣本最優不同」結論
   - 預期時間：1 hr

### 6.5.2 Tier 2（驗證穩定性）

3. **Bootstrap Top-17 stability**
   - 500 次重抽，計算每個 cell 進 Top-17 的頻率
   - 保留頻率 ≥90% 的 cells 為「stable Top-17」
   - 預期時間：2 hr

4. **Leave-one-chromosome-out**
   - 移除某 chromosome，重跑 5D cube
   - 檢查 Top-17 是否轉移到該染色體以外的 cells
   - 驗證結論不依賴 HCC1395 特定 LOH hotspot
   - 預期時間：3 hr

### 6.5.3 Tier 3（若 PI 要求）

5. **真 F1 計算**
   - 取得 SEQC2 FN 資料
   - 在 FN+TP+FP 的 ground truth 集上重算 precision/recall
   - 比對 scheme filter 後 F1 變化
   - 預期時間：視 FN 資料可用性

6. **Logistic regression baseline**
   - 用 5D categorical features + 其他 (NumCpGs, AlleleDelta) 訓練 logistic
   - 比較 empirical Top-17 vs logistic predicted top-k
   - 若 logistic 明顯勝 → empirical 有過擬合；否則 empirical 可信
   - 預期時間：3 hr

### 6.5.4 Tier 4（延伸方向）

7. **flag=on 重跑驗證 NG 貢獻**
   - 整條 pipeline 在 `--germline-hp-only` flag=on 下重跑
   - 驗證 Top-17 去除 NG≥3 cells 後，剩餘 cells 是否仍 Pareto-dominate S3
   - 屬獨立 HP-only 研究範疇

8. **HCC1395 TO Group A/C biology 深度檢視**
   - 將 Top-17 中 Group A (LOH + Extreme + NG=3) 和 Group C (LOH_Noise + high NR) 的 variant 位置取出
   - IGV 視覺檢視：是否為真實 somatic pattern
   - 對照 HCC1395 known LOH hotspots
   - 預期時間：1 day

---

## 6.6 論文定位（暫論）

**若本研究寫 paper**：
- **主張**：方法學貢獻 = 5D cube + Pareto envelope 作為 ISM parameter space exploration 框架
- **不主張**：Top-17 as universal filter；需 H4 驗證後才能主張
- **biology story**：四個 Group (A-D) 各對應一個 epigenetic somatic event 類別
- **風險 reviewer 要求**：跨樣本 generalization（需 LOSO + 真 F1）

**替代路徑**（若 LOSO 失敗）：
- 改投 characterization-focused journal
- 將 5D cube 定位為 "HCC1395-specific cancer genome analysis"
- 放棄 universal filter 主張

---

## 6.7 本研究資料夾不做的事

為避免與其他研究方向重疊，**本資料夾明確不涵蓋**：

- ❌ normal BAM 整合（屬 Phase 2 read-level characterization）
- ❌ LOH.bed 重新產生（屬 LOH pipeline 研究）
- ❌ HP integer tag fix（已完成，見 `project_hp_integer_tag_fix.md`）
- ❌ Zone-Aware framework QS 整合（已 NEGATIVE，見 `project_zone_aware_framework.md`）
- ❌ TO mode 其他 6 樣本 rerun（C++ pipeline 範疇）
- ❌ Deep learning / XGBoost feature extraction
- ❌ Per-sample 另訓一個 Top-k list 的 auto-tuning pipeline

所有 ❌ 項目若 PI 認為必要，應開 **新的 research 資料夾** 而非擴展本資料夾。

---

## 6.8 結語

本研究在 HCC1395 TO 單樣本上實現了 E6（5D 綜觀）的方法學突破，證明 Top-17 Pareto-dominates 既有 S3/S5 biology-informed scheme。**主結論穩健**（大 cell 層級 n=1,099），但**泛化性未證**。

下一輪工作的明確動作：**LOSO（§6.5.1 item 1）— 30 min 可確認 H4 成立與否**。
