# 2026-01-07 顯著性篩選改進計畫

## 問題背景

根據 TP/FP 分析結果，發現以下核心問題：

| 問題 | 現況 | 影響 |
|:---|:---|:---|
| FP 顯著率 > TP | 40.8% vs 28.6% | 顯著性指標反向預測 |
| FP 高聚集 | chr9:41.7Mb 有 78 個位點/50kb | 聚集區域產生大量假顯著 |
| ROC AUC ≈ 0.5 | 所有指標無鑑別力 | 無法區分真假陽性 |
| FP Strong > TP | 37.1% vs 24.5% | VerificationClass 判定有誤 |

---

## 根本原因分析

1. **Gating 過於寬鬆**：僅依據 p-value ≤ 0.05，未加入效應量 (Cramér's V) 門檻
2. **缺少深度過濾**：低深度位點容易產生假顯著
3. **無聚集偵測**：高聚集區域被重複判定為顯著
4. **Significant 定義不嚴格**：`PassedGating AND GlobalP ≤ 0.05` 過於寬鬆

---

## 提議修改

### 第一部分：新增效應量閾值 (Cramér's V Gate)

#### [MODIFY] [GlobalTest.cpp](file:///big8_disk/liaoyoyo2001/InterSubMod/src/core/GlobalTest.cpp)

```diff
 void GlobalTest::apply_gating(GlobalTestResult& result) {
-    result.passed_gate = (result.fisher_ffh.p_value <= config_.gating_p_threshold);
+    // 雙重閾值：p-value 且 Cramér's V 效應量
+    bool p_pass = (result.fisher_ffh.p_value <= config_.gating_p_threshold);
+    bool v_pass = (result.cramers_v >= config_.min_cramers_v);
+    result.passed_gate = p_pass && v_pass;
 }
```

#### [MODIFY] [GlobalTest.hpp](file:///big8_disk/liaoyoyo2001/InterSubMod/include/core/GlobalTest.hpp)

```cpp
struct GlobalTestConfig {
    double gating_p_threshold = 0.05;
    double min_cramers_v = 0.1;  // NEW: 最低效應量閾值
    // ...
};
```

---

### 第二部分：新增最低深度過濾

#### [MODIFY] [RegionProcessor.cpp](file:///big8_disk/liaoyoyo2001/InterSubMod/src/core/RegionProcessor.cpp)

在 `write_significance_summary` 中修改顯著性判定：

```diff
-bool is_significant = r.passed_gating && (r.global_p_value <= 0.05);
+// 嚴格定義：Gating + P-value + 最低深度 + 效應量
+bool is_significant = r.passed_gating && 
+                      (r.global_p_value <= 0.05) && 
+                      (r.num_reads >= 20) &&           // 最低深度
+                      (r.cramers_v >= 0.1);            // 效應量門檻
```

---

### 第三部分：新增「無法判別」欄位

在輸出 CSV 中新增 `Undetermined` 欄位，標記以下狀況：

| 條件 | 原因 |
|:---|:---|
| `num_reads < 20` | 深度不足 |
| `permanova.valid == false` | 距離矩陣 NaN 過多 |
| `dispersion_warning == true` | PERMANOVA 被離散度驅動 |

#### [MODIFY] [RegionProcessor.cpp](file:///big8_disk/liaoyoyo2001/InterSubMod/src/core/RegionProcessor.cpp)

新增 `Undetermined` 和 `UndeterminedReason` 欄位到 CSV 輸出。

---

### 第四部分：強化 compare_vcf_results.py

#### [MODIFY] [compare_vcf_results.py](file:///big8_disk/liaoyoyo2001/InterSubMod/tools/compare_vcf_results.py)

新增以下分析功能：

1. **聚集偵測**：計算 50kb 窗口內顯著位點數，標記高聚集區域
2. **排除高聚集後重計 ROC**：驗證排除後是否改善鑑別力
3. **深度分層分析**：比較 <20, 20-50, >50 reads 各層的 TP/FP 分布

---

## 驗證計畫

### 自動化測試

```bash
# 編譯
cd /big8_disk/liaoyoyo2001/InterSubMod
mkdir -p build && cd build
cmake .. && make -j$(nproc)

# 運行批量分析
./scripts/run_batch_vcf_analysis.sh

# 運行增強版比較
python3 tools/compare_vcf_results.py \
    --output-dir output/.../analysis \
    --labels TP FP \
    --paths output/.../filtered_snv_tp output/.../filtered_snv_fp
```

### 預期改善

| 指標 | 修改前 | 預期修改後 |
|:---|:---|:---|
| TP 顯著率 | 28.6% | ~25% (過濾弱信號) |
| FP 顯著率 | 40.8% | <20% (排除假顯著) |
| ROC AUC (Cramér's V) | 0.519 | >0.6 |

---

## 檔案變更清單

| 動作 | 檔案 |
|:---|:---|
| MODIFY | `src/core/GlobalTest.cpp` |
| MODIFY | `include/core/GlobalTest.hpp` |
| MODIFY | `src/core/RegionProcessor.cpp` |
| MODIFY | `tools/compare_vcf_results.py` |
