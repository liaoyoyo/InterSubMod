<!--
建立時間: 2026-03-24
目標: 深度追蹤 paired_persistent_final_fp 的機制，確認這些 FP 是否可被 ISM 改善
處理範圍: HCC1395(87個) + HCC1395_DORADO(77個) 持續性 FP，45個跨平台共享確認案例
關聯檔案:
  - big7_disk_output/synthesis/observation_workspaces/20260323_to_residual_fp_deep_dive/paired_persistent_candidate_regions.tsv
  - docs/reports/validated/2026/03/20260323_TO_residual_FP_deep_dive_01.md (前置分析)
  - docs/methodology/20260324_方法學審查全域結論報告_01.md (方法學背景)
-->

# Paired Persistent Final FP 深度追蹤報告

**狀態：分析完成 — 確認為不可降低 FP（Irreducible FP），研究方向關閉**

---

## 研究背景

在 TO residual FP 深度追蹤（2026-03-23）中，識別出兩類特殊 FP：
- **paired_persistent**（87 HCC1395 + 77 DORADO）：所有 rule sweep 均無法過濾
- 其中 **45 個跨平台共享確認案例**：在 5kHz 和 DORADO 兩個平台同時持續存在

本報告針對這些案例進行機制深度分析，確認是否存在任何尚未嘗試的 ISM 改進方向。

---

## 數據來源

| 來源文件 | 內容 |
|---------|------|
| `paired_persistent_candidate_regions.tsv` | 164 行（87 HCC1395 + 77 DORADO），含完整 ISM 特徵 |
| `stage_rule_scan.tsv` | 所有已測試 rule sweep 結果（全部 negative） |
| `analysis_report.md` | 前置深度追蹤報告 |

---

## 一、基本統計

### HCC1395 87 個持續性 FP

| VerificationClass | 數量 | 比例 |
|-----------------|------|------|
| Noise | 45 | 52% |
| Weak | 19 | 22% |
| Strong | 18 | 21% |
| Subclone | 5 | 6% |

| 特徵 | 中位數 | 範圍 |
|------|--------|------|
| AF | 0.20 | 0.07 ~ 0.58 |
| GQ | 14 | 4 ~ 26 |
| AlleleDelta | 0.009 | -0.065 ~ 0.136 |
| QualityScore | 75 | 25 ~ 100 |
| PotentialLOH | 67% True | — |

### 45 個跨平台共享確認案例

| VerificationClass | 數量 | 比例 |
|-----------------|------|------|
| Noise | 22 | 49% |
| Strong | 12 | 27% |
| Weak | 7 | 16% |
| Subclone | 4 | 9% |

**關鍵觀察：所有 45 個共享案例的 suggest_filter = False**（所有 rule sweep 均未觸發）

---

## 二、機制分類

### 機制 A：Strong/Weak FP — 真實甲基化信號但來源錯誤（27+16% = 43%）

這些 FP 有**真實的 ISM 甲基化差異信號**，但信號來源是胚系 ASM，而非體細胞變異：

```
機制：
  1. HCC1395 攜帶胚系 SNP 在特定位置
  2. 該位置有 LOH 或 germline ASM（等位基因特異性甲基化）
  3. 甲基化模式在兩個單倍型間確實有差異（ISM 正確偵測到）
  4. 體細胞 calling 錯誤地標記同位置的常見 SNP 為體細胞變異
  5. ISM 看到甲基化差異 → Strong/Weak（正確的 ISM 判斷，錯誤的前提）
```

**Strong FP（12/45）特徵：**
- AlleleDelta 中位數 = 0.067（相對較高，真實 HP 差異）
- QualityScore 中位數 = 75
- HP_assign_rate 中位數 = 0.548（中等 haplotagging 品質）
- PotentialLOH = 7/12（58% LOH 驅動）

**為何無法過濾：**
- ISM 偵測的甲基化差異是**真實存在的**（非假象）
- 問題在於體細胞 caller 的假陽性，而非 ISM 邏輯的錯誤
- 提高 gate 嚴格度會同等損失真 TP（TP 的 CramersV/AlleleDelta 分佈與這些 FP 高度重疊）

### 機制 B：Noise FP — 無甲基化信號但無法過濾（49%）

```
機制：
  1. 體細胞 caller 產生低信心 SNV（GQ 4-26，AF 7-58%）
  2. ISM 分析此位置：無統計顯著甲基化差異（Noise class）
  3. 沒有任何 ISM 特徵能識別這些 Noise FP
  4. 提高過濾閾值（如 Noise → filter）會同等損失 TP
      (33~67% 的 Noise class 是 TP，見方法學審查報告)
```

**Noise FP（22/45）特徵：**
- AlleleDelta 中位數 ≈ 0.001（幾乎無等位基因差異）
- QualityScore 中位數 = 75（與全局相同，無差別）
- PotentialLOH = 21/22（95% 有 LOH 信號）
- **無任何特徵能區分這些 Noise FP 與 Noise TP**

---

## 三、所有嘗試過的 rule sweep（全部 negative）

| 規則 | 觸發 shared FP | 問題 |
|------|-------------|------|
| QS < 50 | 1/45（2%） | 不足 |
| AlleleDelta > 0.15 | 0/45（0%） | 完全無效 |
| GQ < 5 | 1/45（2%） | 不足 |
| hp_assign_rate < 0.2 | 1/45（2%） | 不足 |
| Noise + LOH filter | 21/45（47%） | TP 損失率相似（已驗證 negative，見 N+1 假設） |
| AlleleDelta > 0.25 (hard filter) | 0/45（0%） | 跨樣本 AUROC 不一致 |

**結論：所有特徵組合均無法以可接受的 TP 損失率過濾這些 FP。**

---

## 四、為何 ISM 無法改善此問題

### 核心限制

```
ISM 的設計假設：
  - 如果 HP1 和 HP2 有甲基化差異（passed_gate=True）
  - 且有統計顯著的 label 差異（label_significant=True）
  - → 這個位點有「亞克隆甲基化特徵」

此假設在以下情況下失效：
  - 胚系 SNP 造成的 germline ASM（等位基因特異性甲基化）
  - 這類信號「完全符合」ISM 的亞克隆甲基化特徵定義
  - 但它來自胚系多態性，而非體細胞亞克隆分化
```

### 需要的額外信息

改善這些 FP 需要以下任一：

| 方法 | 所需資源 | 可行性 |
|------|---------|--------|
| 正常樣本甲基化參考（matched normal） | 正常 BAM + 甲基化數據 | TO pipeline 本質限制（無 normal） |
| 人群 germline ASM 資料庫 | 大型 cohort 甲基化 QTL 數據 | 研究性可行，實際部署複雜 |
| 非甲基化體細胞特徵（如 VAF 模式） | 需要額外 caller 信息 | Tier 3 改進方向 |

---

## 五、跨平台一致性分析

**45 個跨平台共享案例的特別重要性：**

```
HCC1395（5kHz）持續 FP ∩ HCC1395_DORADO（非 5kHz）持續 FP = 45 個

這表示這些 FP：
  1. 不是平台噪聲（5kHz vs DORADO 重現）
  2. 不是隨機甲基化噪音（跨平台一致）
  3. 是系統性胚系 ASM 信號
```

**這些 FP 的 ISM 分析一致性：**

| 特徵 | HCC1395 | DORADO | 一致性 |
|------|---------|--------|--------|
| VerificationClass 分佈 | Noise(49%)/Strong(27%)/Weak(16%) | 相似 | 高 |
| PotentialLOH | 75% | 相似 | 高 |
| paired_raw_filter | 100% PASS | 100% PASS | 完全一致 |

---

## 六、結論與研究方向評估

### 最終判決

| 問題 | 結論 |
|------|------|
| 這些 FP 是否可被 ISM 改善？ | **否** — 任何 ISM 內部規則均無法以可接受 TP 損失率過濾 |
| 機制是否已確認？ | **是** — germline ASM（Strong/Weak）+ caller FP 無 ISM 可判別特徵（Noise） |
| 是否需要繼續探索 ISM 規則？ | **否** — 所有規則空間已窮舉，研究潛力 = exhausted |
| 改善路徑是否存在？ | **條件性** — 需要 normal sample 甲基化參考或 germline ASM 資料庫 |

### 量化影響

這 87 個 HCC1395 持續性 FP 佔：
- 全量 FP（11843）的 **0.7%**
- 即使全部移除，F1 僅提升 ≈ 0.0005
- **不值得為 0.7% FP 投入 Tier 3 C++ 改動**

### 後續行動

| 行動 | 優先級 | 說明 |
|------|-------|------|
| 關閉此研究方向 | 立即 | 研究潛力 = exhausted，不再探索 |
| 記錄為已知限制 | 完成 | 見方法學審查全域結論報告 |
| 長期：germline ASM 過濾 | 研究性 | 需要外部資料庫，非短期目標 |

**研究方向狀態：CLOSED — paired_persistent_final_fp 為 ISM 架構下的不可降低 FP（irreducible FP），根本原因是 germline ASM 與體細胞 calling FP 的混合，無法在缺乏正常樣本甲基化參考的情況下改善。**
