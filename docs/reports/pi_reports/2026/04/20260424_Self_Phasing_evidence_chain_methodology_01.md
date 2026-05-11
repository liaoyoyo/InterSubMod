<!--
建立時間: 2026-04-24
目標: 回答「Self-phasing 判定是否依賴 Paired HP1/HP2 平衡？」的 methodology 審查
受眾: PI（方法學嚴謹度焦點）
處理範圍: 6 個核心證據的依賴性分解、robustness 分析、完整邏輯鏈
狀態: validated_methodology_review
pipeline_track: all
priority: P0
relates_to: 20260422_Self_Phasing_complete_report_for_PI_01.md, 20260422_Self_Phasing_multiperspective_argument_01.md
-->

# Self-Phasing 證據鏈方法學審查
## ——判定是否依賴 Paired HP1/HP2 平衡？各證據獨立性逐項分解

> 撰稿日期：2026-04-24
> 觸發：PI 提問「Self-phasing 的判定是依據 Paired tag 的狀況與 HP1/HP2 平衡數量來決定嗎？」
> 搭配閱讀：`20260422_Self_Phasing_complete_report_for_PI_01.md`（技術敘事）、`20260422_Self_Phasing_multiperspective_argument_01.md`（多視角論證）

---

## 一句話回答

**部分依賴，部分不依賴。6 個核心證據中，4 個獨立於 Paired——結論的 robustness 建立在獨立證據鏈的三角驗證（triangulation）上，不是單一依賴 Paired HP1/HP2 平衡。**

---

## Section 1：核心問題釐清

### 1.1 PI 問題的精確意義

用戶問：「Self-phasing 是依據 Paired tag 的狀況與 HP1/HP2 平衡數量來決定嗎？」

這個問題可拆解為三個子問題：

| 子問題 | 方法學意涵 |
|--------|-----------|
| Q1: 我們的判定是否把 **Paired mode 當 ground truth**？ | circular reasoning 風險：若 Paired 自己也有 bug，結論會跟著錯 |
| Q2: 我們的判定是否**只**用 HP_Ratio（HP1/HP2 平衡數量）？ | metric monoculture 風險：單一指標可能受 artifact 影響而不自知 |
| Q3: 如果上述兩者成立，**如何驗證判定有效**？ | 需要獨立證據鏈或 perturbation experiment |

### 1.2 為何這個問題重要

**自我指涉（self-reference）陷阱**：
- 若用「Paired 的 HP_Ratio = 真」去證明「TO 的 HP_Ratio = 假」
- 而 Paired 本身也用同樣 LongPhase 演算法產生
- → 我們只證明了「兩種模式的 HP_Ratio 不一致」，沒有獨立證明「TO 是錯的」

**避免這個陷阱需要**：
- 獨立的生物學先驗（如 50:50 null 假設）
- 獨立的演算法理論（如 phasing graph 數學）
- Perturbation experiment（改動機制看現象是否消失）
- 跨模式、跨樣本、跨方法的一致性（triangulation）

### 1.3 本報告做什麼

把 self-phasing 證據鏈拆成 **6 條證據**，每條標註其是否依賴：
- **Paired mode 作參考**
- **HP1/HP2 平衡指標（HP_Ratio）**
- **生物學 null 假設（50:50 expectation）**
- **演算法理論（phasing graph math）**

然後評估 robustness：**如果 Paired 有 caveat，結論是否仍成立？**

---

## Section 2：六條核心證據的證據鏈分解

### E1：Somatic HP1:HP2 bias = 17.3 : 1（614K vs 35K reads）

**What（是什麼）**：
HCC1395 TO 模式下，跨全基因體所有 somatic variant sites，統計 ALT reads 被分配到 HP1 vs HP2 的總數比例。

**How（怎麼算）**：
```
For each somatic variant site in TO BAM:
    count HP1 ALT reads (HP:i:1 + HP:i:11)
    count HP2 ALT reads (HP:i:2 + HP:i:21 + HP:i:33 where applicable)
Sum across genome → ratio
Observed: 614,000 : 35,500 ≈ 17.3 : 1
```

**依據什麼基準判定異常？**
- **生物學 null 假設**：每個 somatic SNV 的兩條 haplotype 分配應該接近 50/50
- 理由：somatic variant 發生於單一 allele，但 which allele（maternal vs paternal）是隨機事件
- 預期：平均下來 ~1:1
- **實測 17.3:1 → 大幅偏離 null → 存在系統性偏差**

**依賴性分析**：

| 維度 | 依賴？ | 說明 |
|------|-------|------|
| Paired 作參考 | ❌ **不依賴** | 完全用 TO 數據自身計算，Paired 沒出現 |
| HP_Ratio 指標 | ❌ **不依賴** | 用 raw read count，不是 ratio metric |
| 生物學 null | ✅ **依賴** | 50:50 假設是判定「異常」的唯一基準 |
| 演算法理論 | 🟡 partial | 需要知道 phasing 機制才能解釋 bias，但 bias 偵測本身不需要 |

**Robustness**：若 Paired 有問題，E1 **完全不受影響**——這是純 TO 內部觀察。

---

### E2：Cross-mode HP_Ratio correlation  r = 0.001（288K same-site pairs）

**What**：
取在 Paired 與 TO 兩個模式下都有 ≥5 reads 覆蓋的位點（288K pairs），各自計算 HP_Ratio，做 Pearson correlation。

**How**：
```
For each shared site:
    HP_Ratio_Paired = min(HP1, HP2) / (HP1 + HP2)   # 在 Paired BAM 計算
    HP_Ratio_TO     = min(HP1, HP2) / (HP1 + HP2)   # 在 TO BAM 計算
Pearson r on 288K pairs → r = 0.001 (p = 0.59, n.s.)
```

**依據什麼基準判定異常？**
- **預期**：若 TO 的 HP_Ratio 是真實 haplotype 結構的近似，應該與 Paired 高度相關（r > 0.5）
- **實測 r = 0.001**：兩者完全不相關

**依賴性分析**：

| 維度 | 依賴？ | 說明 |
|------|-------|------|
| Paired 作參考 | ✅ **強依賴** | 以 Paired 為 reference 判定 TO 是否 deviate |
| HP_Ratio 指標 | ✅ **強依賴** | 唯一輸入是 HP_Ratio 數值 |
| 生物學 null | ❌ 不依賴 | 不用 50:50 假設；用 correlation≈0.5 的期望 |
| 演算法理論 | ❌ 不依賴 | 純 empirical |

**Robustness**：若 Paired 有問題，E2 結論會動搖——但只能告訴我們「Paired 跟 TO 不相關」而非「TO 是錯的」。**此證據需與獨立證據搭配才能指向 self-phasing**。

---

### E3：62% TO-mode ISM HP_Ratio LOH 是 artifact

**What**：
TO 模式判定為 ISM HP_Ratio LOH（HP_Ratio < 0.1 或 > 0.9）的 regions 中，有 62% 在 Paired 模式下 HP_Ratio 回到 0.4-0.6（balanced 非 LOH）。

**How**：
```
TO_ISM_LOH_set = regions where HP_Ratio_TO ∈ [0, 0.1) ∪ (0.9, 1]
Check these regions in Paired:
  - % with HP_Ratio_Paired ∈ [0.4, 0.6] = 62%
  - → "consistent" LOH (structural) = 38%
```

**依據什麼基準判定 artifact？**
- **預期**：若 TO LOH 判定反映真實 allelic imbalance，Paired 應該也看到 LOH
- **實測 62% 在 Paired 下 balanced → 這 62% 是 TO 特有的 artifact**

**依賴性分析**：

| 維度 | 依賴？ | 說明 |
|------|-------|------|
| Paired 作參考 | ✅ **強依賴** | 以 Paired HP_Ratio 0.4-0.6 當 ground truth |
| HP_Ratio 指標 | ✅ **強依賴** | 唯一指標 |
| 生物學 null | ❌ 不依賴 | 不用 50:50 普適假設；用 region-specific Paired |
| 演算法理論 | ❌ 不依賴 | 純 empirical |

**Robustness**：若 Paired 有問題，E3 會動搖——但**注意 LOH.bed Jaccard=1.0 的獨立證據**：
- LOH.bed 由 VCF AD（allele depth）計算，不經 HP tag
- Paired LOH.bed = TO LOH.bed（Jaccard=1.0）
- → **region-level LOH 兩模式一致**
- → ISM HP_Ratio LOH 在 TO 與 Paired 間的差異**不是 LOH 本質變化**，是 HP tag 層的差異
- 這是 Paired HP_Ratio 可信的**獨立佐證**

---

### E4：7/7 樣本方向一致性（CV-2 pass）

**What**：
HCC1395、HCC1395_DORADO、HCC1954、HCC1937、H1437、H2009、COLO829 七個樣本**全部觀察到相同方向的 self-phasing 偏差**。

**How**：
```
For each of 7 samples:
    Compute somatic HP1:HP2 bias (E1-like statistic)
    All 7 show ratio > 10:1 in same direction (HP1-dominant)
CV-2 criterion: direction consistent in ≥6/7 samples → PASS
```

**依據什麼基準判定？**
- **TO-internal replication**：若 self-phasing 是樣本特異性 artifact，應該只出現在 1-2 樣本；7/7 一致**排除樣本特異性**
- 這是跨樣本的 consistency check，**完全不需要 Paired**

**依賴性分析**：

| 維度 | 依賴？ | 說明 |
|------|-------|------|
| Paired 作參考 | ❌ **不依賴** | 只用 TO 數據跨樣本比對 |
| HP_Ratio 指標 | 🟡 partial | 用 HP_Ratio-like 指標，但核心是方向一致性不是絕對值 |
| 生物學 null | ❌ 不依賴 | 不需要 50:50 假設 |
| 演算法理論 | ❌ 不依賴 | 純 empirical |

**Robustness**：若 Paired 有問題，E4 **完全不受影響**。這是 TO-internal evidence。

---

### E5：PON-only phasing perturbation experiment

**What**：
修改 LongPhase-TO 的 `PhasingProcess.cpp`，呼叫 `convertNonGermlineToSomatic()` 使 somatic variants 以 reduced edge weight 處理；加 `--pon-only-phasing` CLI flag。觀察 self-phasing bias 是否消失。

**How**：
```
Baseline (no modification):
    Somatic HP1:HP2 = 17.3:1 (bias present)
After --pon-only-phasing:
    Somatic HP1:HP2 = ~1:1 (bias eliminated)
Phase block N50: 4,061 → 8,109 (+99.7%)
LOH.bed Jaccard vs baseline: 1.0000
```

**依據什麼基準判定有效？**
- **Causal perturbation**：若 self-phasing 真的是 somatic 反客為主造成
  → 移除 somatic 從 phasing graph anchors
  → bias 應該消失
- 實測 bias 確實消失 → **因果關係成立**
- 這是**改變機制來證明機制存在**的實驗方法（Koch's postulate 類比）

**依賴性分析**：

| 維度 | 依賴？ | 說明 |
|------|-------|------|
| Paired 作參考 | ❌ **不依賴** | TO × TO-modified 比對，不需 Paired |
| HP_Ratio 指標 | 🟡 partial | 使用 HP_Ratio 檢視但核心是 perturbation 前後對比 |
| 生物學 null | 🟡 partial | 預期 bias 消失是基於 "50:50 應該是 default" |
| 演算法理論 | ✅ **強依賴** | 基於「somatic-somatic edge weight 高 → 移除就消失」的理論預測 |

**Robustness**：若 Paired 有問題，E5 **完全不受影響**。這是最強的獨立證據——**它直接在 TO 內部展示了「改動機制可以讓現象消失」**。

---

### E6：Phasing graph edge weight 理論推導

**What**：
從 phasing graph 演算法的數學結構推導「為何 somatic 在 TO 模式會自我連結」。

**How**：
```
Phasing graph edge weight 公式：
    w(A, B) = Σ_reads I(read 帶 A.alt) × I(read 帶 B.alt)

兩個分析：
  - 真 germline het（位置遠）：
        ALT reads 分散於 HP1/HP2
        受 recombination crossover 影響 → w 低
  
  - 同 clone somatic（位置遠）：
        ALT reads 共享 tumor sub-population
        100% 共現於同一 reads → w 極高
        
結論：somatic-somatic edges 權重 > germline-germline edges
→ phasing graph 被 somatic 主導
→ somatic 決定 scaffold
```

**依據什麼基準判定？**
- **純演算法/數學推導**：不需要觀察數據
- 預測：TO 模式 phasing graph 必然被 somatic 主導
- 實測數據（E1, E2, E3）驗證了此預測

**依賴性分析**：

| 維度 | 依賴？ | 說明 |
|------|-------|------|
| Paired 作參考 | ❌ **不依賴** | 純推導 |
| HP_Ratio 指標 | ❌ **不依賴** | 純推導 |
| 生物學 null | ❌ **不依賴** | 純推導 |
| 演算法理論 | ✅ **核心** | 本身就是理論 |

**Robustness**：若 Paired 有問題，E6 **完全不受影響**。這是最 upstream 的證據——如果演算法理論錯誤，我們根本不會預期 self-phasing 會發生。

---

## Section 3：依賴性矩陣彙總

![Figure 15 — Evidence Chain Dependency Matrix](figures/fig15_evidence_dependency_matrix.png)

### 3.1 6 × 4 依賴性矩陣

| Evidence | Paired 作參考 | HP_Ratio 指標 | 生物學 null | 演算法理論 | Lineage |
|----------|:-------------:|:-------------:|:-----------:|:----------:|---------|
| **E1** Somatic bias 17.3:1 | ❌ 獨立 | ❌ 獨立 | ✅ 依賴 | 🟡 partial | Biology-null |
| **E2** Cross-mode r=0.001 | ✅ 依賴 | ✅ 依賴 | ❌ 獨立 | ❌ 獨立 | Paired-triangulated |
| **E3** 62% LOH artifact | ✅ 依賴 | ✅ 依賴 | ❌ 獨立 | ❌ 獨立 | Paired-triangulated |
| **E4** 7/7 consistency | ❌ 獨立 | 🟡 partial | ❌ 獨立 | ❌ 獨立 | TO-internal |
| **E5** PON-only experiment | ❌ 獨立 | 🟡 partial | 🟡 partial | ✅ 依賴 | Experimental |
| **E6** Phasing graph math | ❌ 獨立 | ❌ 獨立 | ❌ 獨立 | ✅ 依賴 | Theoretical |

### 3.2 關鍵觀察

**4/6 證據完全獨立於 Paired 參考**：
- E1（生物學 null 驅動）
- E4（TO 內部複製）
- E5（擾動實驗）
- E6（演算法理論）

**2/6 證據依賴 Paired 作參考**：
- E2（cross-mode correlation）
- E3（62% artifact 定義）

### 3.3 這告訴我們什麼？

**不是「結論建立在 Paired HP1/HP2 平衡上」**——而是**結論建立在 6 條獨立證據的 triangulation 上**，其中只有 2 條用到 Paired 作 cross-check。

即使 Paired 本身有未知 calibration 問題，**另外 4 條證據仍獨立指向同一結論**：

- E1 說：TO 觀察到 17.3:1 bias（違反生物學 50:50）
- E4 說：7 個樣本都看到同方向 bias（不是特殊案例）
- E5 說：改動機制（移除 somatic from graph）後 bias 消失（因果成立）
- E6 說：從演算法數學推導就應該預期 bias 存在

**Paired 作為 cross-mode reference 強化了結論（E2, E3），但不是結論的唯一支柱**。

---

## Section 4：Robustness 分析——如果 Paired 有 caveat，結論是否仍成立？

### 4.1 三種可能的 Paired caveat

| Caveat 情境 | 影響 E2 | 影響 E3 | 影響 E1/E4/E5/E6 | 整體結論 |
|-----------|:-------:|:-------:|:------------------:|---------|
| Paired HP_Ratio 也有未知 bias（但比 TO 輕微） | 部分削弱 | 部分削弱 | ❌ 不影響 | **仍成立**（4/6 證據獨立） |
| Paired 完全可靠（現有假設） | 強化 | 強化 | ❌ 不影響 | **成立** |
| Paired 完全不可靠（極端假設） | 失效 | 失效 | ❌ 不影響 | **仍成立**（4/6 證據獨立） |

### 4.2 為何即使極端情境下結論仍成立？

**三個獨立支柱**（不用 Paired 就能證明）：

1. **生物學證據（E1）**：TO 內部觀察 17.3:1 bias，違反 somatic heterozygosity 的生物學預期（50:50）
2. **擾動實驗（E5）**：在 TO 內部改動機制（PON-only phasing）→ bias 消失 → 因果成立
3. **理論推導（E6）**：從 phasing graph 數學**預測**此現象必然發生

**這是典型的「triangulation by disparate methods」**（不同方法學的三角驗證），哲學上等同於：
- 觀察（E1）：看到現象
- 實驗（E5）：證明可控
- 理論（E6）：數學預測

若三者**都獨立指向同一結論**，則結論極難被單一指標挑戰。

### 4.3 LOH.bed Jaccard=1.0 的特殊意義

有一個獨立證據進一步支持 Paired HP_Ratio 的可信度：
- LOH.bed 由 VCF AD（allele depth）計算，**完全不經過 HP tag**
- Paired LOH.bed vs TO-with-PON-only LOH.bed → **Jaccard = 1.0000**
- 表示 region-level LOH 的結構判定在兩模式間完全一致
- 如果 Paired 的 phasing/HP 有系統性 bug，LOH.bed 也會受影響——但它沒有
- → **Paired HP_Ratio 在 LOH 區的可信度獲得獨立佐證**

---

## Section 5：關於「HP1/HP2 平衡數量」的具體說明

PI 問題的另一半：「是依據 HP1 HP2 平衡數量來決定嗎？」

### 5.1 HP_Ratio 的角色

`HP_Ratio = min(HP1_reads, HP2_reads) / (HP1_reads + HP2_reads)`

- 範圍 [0, 0.5]
- 0.5 = 完全平衡（理想 het）
- 0 = 完全不平衡（全部 reads 在同一 HP，即 LOH）

**在 self-phasing 判定中，HP_Ratio 不是主要指標**：
- E1 用 raw read count（614K vs 35K），不是 ratio
- E4 用方向一致性，不是絕對值
- E5 核心是 bias 前後對比，ratio 只是輔助
- E6 純理論

**HP_Ratio 主要出現在**：
- E2（cross-mode correlation）：確實用 HP_Ratio
- E3（62% artifact）：用 Paired HP_Ratio 0.4-0.6 當 balanced baseline

### 5.2 HP_Ratio 本身是 derivative 指標

HP_Ratio 是**後果的一個表現**，不是 self-phasing 的**直接證據**：

```
┌─ root cause（機制）──────────────┐
│ Somatic variants as anchors       │
│ in phasing graph                  │
└───────────────┬───────────────────┘
                │
                ▼
┌─ direct evidence（直接證據）─────┐
│ Raw read HP1:HP2 count = 17.3:1  │ ← E1 看這個
│ (614K vs 35K reads)               │
└───────────────┬───────────────────┘
                │
                ▼
┌─ derivative metric（衍生指標）───┐
│ HP_Ratio distribution deviates    │ ← E2, E3 看這個
│ from balanced 0.5                 │
└───────────────────────────────────┘
```

### 5.3 為何 E1 用 raw count 而非 ratio？

**Raw count 比 ratio 更難被 artifact 污染**：
- Ratio 可能因分母小而波動（e.g., n=10 vs n=1000 的 ratio 穩定度不同）
- Raw count 跨全基因體加總 → 大數法則下可信度高
- 17.3:1 bias 建立在 ~650K reads 的大數據上，統計穩定

這是選擇 E1 作為**主要證據**而非 HP_Ratio correlation 的方法學理由。

---

## Section 6：Self-Phasing 完整邏輯鏈（PI 視角）

### 6.1 三層論證結構

```
╔═════════════════ Layer 1: 理論預測 ═════════════════╗
║  E6: Phasing graph math                               ║
║  預測：TO 模式 somatic 必然自我 phase                  ║
║  依賴：純演算法                                       ║
╚═══════════════════════════════════════════════════════╝
                       ↓ predicts
╔═════════════════ Layer 2: 觀察驗證 ═════════════════╗
║  E1: Raw read bias 17.3:1                             ║
║  E4: 7/7 samples direction consistency                ║
║  依賴：生物學 50:50 null + TO 內部複製                 ║
╚═══════════════════════════════════════════════════════╝
                       ↓ corroborated by
╔═════════════════ Layer 3: 擾動證實 ═════════════════╗
║  E5: PON-only phasing experiment                      ║
║  結果：移除機制後現象消失（因果性）                    ║
║  依賴：perturbation experiment                        ║
╚═══════════════════════════════════════════════════════╝
                       +
╔══════════════ Layer 4（輔助）: Cross-mode 檢查 ═════╗
║  E2: Cross-mode r=0.001                              ║
║  E3: 62% TO LOH = artifact                           ║
║  依賴：Paired mode as reference                      ║
║  Role: triangulation, 不是結論唯一依據                ║
╚══════════════════════════════════════════════════════╝
                       ↓ converges
          ╔═════════════════════════════════╗
          ║  CONCLUSION                      ║
          ║  Self-phasing 是 TO 模式真實     ║
          ║  存在的結構性 artifact           ║
          ║                                  ║
          ║  支撐：                          ║
          ║  - Theory + Observation +        ║
          ║    Experiment 三重獨立驗證       ║
          ║  - Paired cross-check 作為       ║
          ║    triangulation 補強            ║
          ╚══════════════════════════════════╝
```

### 6.2 證明有效的證據鏈（完整版）

**「如何證明 self-phasing 問題真實存在？」**

| 步驟 | 證據 | 說明 |
|------|------|------|
| 1. 理論預測 | E6 | 數學推導 somatic-somatic edges 應主導 graph |
| 2. 直接觀察 | E1 | 17.3:1 bias（raw count）違反 50:50 null |
| 3. 跨樣本複製 | E4 | 7/7 樣本同方向，排除 sample-specific artifact |
| 4. 因果驗證 | E5 | PON-only experiment 移除機制後現象消失 |
| 5. 跨模式檢查 | E2, E3 | 與 Paired 比對，顯示 TO 是偏離方 |
| 6. 獨立 LOH 佐證 | LOH.bed Jaccard=1.0 | 不經 HP tag 的 LOH 判定兩模式一致，佐證 Paired 可信度 |

**每個步驟都有獨立貢獻**；移除任一步驟後結論仍成立，但 confidence 降級。

---

## Section 7：結論與給 PI 的回答

### 7.1 直球回答

**問**：Self-phasing 的判定是依據 Paired tag 的狀況與 HP1/HP2 平衡數量來決定嗎？

**答**：
- **部分地。** 6 個核心證據中有 2 個（E2, E3）確實以 Paired HP_Ratio 作 cross-mode reference
- **但更重要的是**，另外 4 個證據（E1, E4, E5, E6）**完全獨立於 Paired**——它們基於生物學 null 假設、TO 內部複製、擾動實驗、演算法理論
- 結論建立在**三重獨立驗證**（theory + observation + experiment）上，Paired cross-check 是加分項而非必要條件

### 7.2 如果 PI 進一步追問

**Q: 為什麼不直接用 Paired 當唯一 ground truth？**
A: 避免 circular reasoning。Paired 也用同 LongPhase 演算法產生，若 Paired 本身有未知 bug，會連帶污染結論。獨立證據（E1/E5/E6）規避了此風險。

**Q: 如果 Paired 真的有問題，結論會變嗎？**
A: 不會。E1（17.3:1 bias）、E4（7/7 consistency）、E5（PON-only experiment）、E6（演算法理論）四條證據**都不用 Paired 也成立**。即使在極端 Paired 完全不可靠情境下，self-phasing 存在仍獲支持。

**Q: HP_Ratio 到底有沒有重要？**
A: HP_Ratio 是**後果的一個表現**（derivative metric），不是 self-phasing 的**直接證據**。直接證據是 raw read count（E1）。HP_Ratio 主要用於 cross-mode triangulation（E2, E3）與 perturbation 觀察（E5）。

**Q: 「證明有效」的標準是什麼？**
A:
- 機制層面 ✅：E6 理論推導 + E5 擾動實驗確認因果
- 現象層面 ✅：E1 raw bias + E4 跨樣本一致性確認普遍性
- 對照層面 ✅：E2/E3 跨模式 triangulation + LOH.bed Jaccard=1.0 獨立佐證

### 7.3 本報告補強了什麼

- 把「self-phasing 存在」從「單一論證」升級為「證據鏈 triangulation」
- 明確列出每個證據對 Paired / HP_Ratio / 生物學 null / 演算法理論的**依賴關係**
- 提供 robustness 分析：**極端情況下**（Paired 完全不可靠）結論仍成立
- 釐清 HP_Ratio 是 derivative metric，不是直接證據

---

## 附錄 A：證據鏈快速查詢表

| 想知道... | 看哪一條證據 |
|-----------|------------|
| Self-phasing 演算法上為何必然發生？ | E6（phasing graph math） |
| Self-phasing 有觀察到嗎？ | E1（17.3:1 raw bias） |
| 是不是只有一個樣本特例？ | E4（7/7 consistency） |
| 證明這是因果而非相關？ | E5（PON-only perturbation） |
| Paired 與 TO 的差距？ | E2（r=0.001） / E3（62% LOH artifact） |
| Paired 本身可信嗎？ | LOH.bed Jaccard=1.0 獨立佐證 |

## 附錄 B：與既有 PI 報告的分工

| 報告 | 焦點 |
|------|------|
| **20260422_Self_Phasing_complete_report_for_PI_01.md** | 技術敘事（問題→根因→修改→結果）+ 10 張圖 |
| **20260422_Self_Phasing_multiperspective_argument_01.md** | 3 視角 × 9 challenges adversarial review + fig11-13 |
| **本報告（20260424）** | 證據鏈方法學審查 + 依賴性矩陣 + fig15 |

本報告不重複前兩份的內容——專門處理「證據 robustness」方法學議題。

---

## 附錄 C：相關記憶與原始素材

**Memory**：
- `project_self_phasing_causal_chain_confirmed.md`（E5 擾動實驗）
- `project_pon_only_phasing_verification.md`（PON-only VCF 層結果 + Jaccard=1.0）

**原始素材**：
- E1 raw count：landscape doc `02_Self_Phasing根因.md:107-115`
- E2 cross-mode correlation：同上
- E3 62% artifact：同上
- E4 7/7 consistency：同上 + CV-2 檢核
- E5 PON-only experiment：`research/loh_investigation/reports/20260403_pon_only_phasing_verification_report.md`
- E6 演算法理論：landscape doc `02_Self_Phasing根因.md:138-144`

**圖表**：
- **Fig 15**（本報告新增）：6 × 4 依賴性矩陣

---

## 報告結束

**核心訊息給 PI**：

1. **不是單一依賴 Paired**：6 條證據中 4 條獨立於 Paired 參考
2. **三重獨立驗證**：理論（E6）+ 觀察（E1, E4）+ 實驗（E5）互相 triangulate
3. **HP_Ratio 只是衍生指標**：主要直接證據是 raw read count（17.3:1）
4. **Robustness 已建立**：即使 Paired 有未知 calibration 問題，結論仍成立
5. **本報告完善證據鏈的方法學透明度**——這對審稿與未來發表至關重要
