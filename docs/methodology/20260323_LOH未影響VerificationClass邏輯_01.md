<!--
建立時間: 2026-03-23
問題類型: C++ 方法 | 分類邏輯
影響 track: TO | 兩者
狀態: pending_decision
-->

# LOH 未影響 VerificationClass 分類邏輯 審查文件

## 問題描述

ISM 計算 `potential_loh`（潛在雜合子缺失）標誌：

**程式碼位置**：`src/core/SignificanceAnalyzer.cpp:312-320`

```cpp
double hp_ratio = 0.5;
bool potential_loh = false;
if (result.hp_multistage.n_hp1_family > 0 || result.hp_multistage.n_hp2_family > 0) {
    int total_hp = result.hp_multistage.n_hp1_family + result.hp_multistage.n_hp2_family;
    hp_ratio = (n_hp1_family + 1) / (total_hp + 2);  // Laplace smoothing
    potential_loh = (hp_ratio < 0.1) || (hp_ratio > 0.9);  // 極端 HP 比例 → LOH
}
```

**問題**：LOH 被偵測到（potential_loh=true），但在 VerificationClass 決策中：

```cpp
if (label_sig && cluster_significant) {
    if (potential_loh) {
        result.verification_class = "Strong";  // 仍然是 Strong！
    } else {
        result.verification_class = "Strong";  // 同樣 Strong
    }
}
```

**LOH 偵測對 VerificationClass 沒有任何影響**。注解說「LOH_Subtype will track this」，但實際上只是寫入 `loh_subtype` 欄位作為參考。

---

## 量化影響

**LOH 在 label_first_metrics.tsv 中的比例**：

在 Case C（chr3:94578813，LOH case）：
- `PotentialLOH=true`
- `VerificationClass=Noise`（AF=0.978，passed_gating=False）
- LOH 直接導致 ClusterPassedGating=False（因為 global test 分組失敗）

**問題**：LOH 的問題不在 VerificationClass 決策，而在更早的 GlobalTest 環節：
- LOH 使 ALT/REF 完全偏斜（所有 reads 都是 ALT 或都是 REF）
- chi-square test 無法比較兩組
- 導致 passed_gate=False

從 label_first_metrics.tsv（HCC1395-5kHz）中分析 LOH 位點：
- 需從 `PotentialLOH` 欄位統計（待確認欄位名稱）

---

## 修改選項

### 方案 A：不修改（LOH 邏輯可接受）
- potential_loh=True 的位點在 TO mode 通常是 ClairS-TO 的 FP（LOH 使 VAF 接近 1）
- 這些位點的 VerificationClass 設為 Strong 並不會被最終過濾器採用
- **F1 影響**：0

### 方案 B：LOH 強制降階到 Weak 或 Subclone
- 修改位置：`src/core/SignificanceAnalyzer.cpp:323-330`
- 若 potential_loh=True → 強制 VerificationClass = "Subclone" 或 "Weak"
- **理由**：LOH 的甲基化信號不代表腫瘤亞克隆，降階可減少 Strong FP
- **風險**：LOH TP 也會被降階（若有 LOH TP 存在）
- **預估 F1 影響**：微小，因為 VerificationClass 未連接最終過濾

### 方案 C：先確認 LOH 在 TP/FP 中的分佈
- 從 label_first_metrics.tsv 統計 PotentialLOH=True 的 TP/FP 比例
- 再決定是否需要修改
- **成本**：純分析，無需改 C++

---

## 建議

**先執行方案 C**（數據確認），再決定 A 或 B。

特別需要確認：
1. HCC1954 的 LOH 比例（HCC1954 F1=0.378，最低，可能 LOH 比例高）
2. Strong FP 中有多少是 potential_loh=True

---

## 用戶決策

**選擇**：[ ] A / [ ] B / [ ] C（先分析）
**日期**：
**理由**：
