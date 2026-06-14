# Verify: claim `threshold_use` — 距離矩陣用連續機率 vs 二值

> 對抗審查日期: 2026-06-13 · 審查員預設懷疑 · 只引實際檔案 + 行號 + 原始碼片段
> 待驗 claim:「距離矩陣用連續機率(非二值);methyl-high/low 閾值只在某些輸出二值化」
> 指派檔案: `src/core/DistanceMatrix.cpp` + Config + 解析路徑

---

## Verdict: **PARTIAL（核心成立，但需精修措辭）**

claim 的**結論在生產實況下成立**，但 claim 的**機制描述過度絕對化**。精確版：

- 距離計算**依 metric 而定**用連續或二值，**不是無條件用連續**。
  - `NHD` / `JACCARD` → 吃 `binary_matrix`（受 methyl-high/low 閾值二值化影響）
  - `L1` / `L2` / `CORR` / `BERNOULLI` → 吃 `raw_matrix`（連續，**不受** methyl-high/low 影響）
- **生產實際用 `BERNOULLI`（連續）** → 故「距離矩陣（生產）用連續機率」**為真**。
- 「methyl-high/low 閾值只在某些輸出二值化」**為真**：閾值唯一作用點是 `binary_matrix`，而 `binary_matrix` 只被 NHD/JACCARD 消費（生產不選），其餘只進 `methylation.csv` 連續輸出不經閾值。

---

## L1 證據鏈（檔:行 + 片段）

### 1. 距離 dispatcher 按 metric 硬選矩陣（非無條件連續）

`src/core/DistanceMatrix.cpp:309-345`（`calculate_distance_impl`）:
```cpp
case DistanceMetricType::NHD: {
    Eigen::VectorXi vec_i = mat.binary_matrix.row(row_i);   // 二值 (311)
    Eigen::VectorXi vec_j = mat.binary_matrix.row(row_j);
    return calculate_nhd(...);
}
case DistanceMetricType::L1: {
    Eigen::VectorXd vec_i = mat.raw_matrix.row(row_i);       // 連續 (317)
    ...
case DistanceMetricType::JACCARD: {
    Eigen::VectorXi vec_i = mat.binary_matrix.row(row_i);   // 二值 (335)
    ...
case DistanceMetricType::BERNOULLI: {
    Eigen::VectorXd vec_i = mat.raw_matrix.row(row_i);       // 連續 (342)
    return calculate_bernoulli(...);
}
```
→ **連續 vs 二值取決於 metric**。claim「用連續(非二值)」對 NHD/JACCARD 不成立；對 L1/L2/CORR/BERNOULLI 成立。

### 2. 生產實際 metric = BERNOULLI（連續）

- `scripts/run_vcf_all_snv.sh:65` `METRICS="BERNOULLI"`（64 行 `# METRICS="NHD"` 被註解掉）
- `scripts/analysis/run_pure_research_round.sh:17` `DISTANCE_METHOD="BERNOULLI"`
- 兩腳本透過 `--distance-metric ${m}` 顯式傳入（`run_vcf_all_snv.sh:315-316`）
- BERNOULLI 進 `calculate_bernoulli`（`DistanceMatrix.cpp:254-301`），全程用連續 `Eigen::VectorXd` raw 值，confidence weight `2·|p−0.5|`（262-264）、Bernoulli 期望差 `p_i(1-p_j)+(1-p_i)p_j`（281）——**完全不碰 binary_methyl_high/low**。

→ **生產距離矩陣 = 連續機率**，claim 結論在實況成立。

### 3. methyl-high/low 閾值「只」用於 binary_matrix（claim 後半成立）

唯一二值化點 `src/core/RegionProcessor.cpp:1410-1423`:
```cpp
double val = raw_matrix[i][j];
if (val < 0) { meth_mat.raw_matrix(i,j)=NAN; meth_mat.binary_matrix(i,j)=-1; }
else {
    meth_mat.raw_matrix(i,j)=val;                              // raw 永遠保連續
    if (val >= binary_methyl_high_)      binary_matrix(i,j)=1; // 0.8
    else if (val <= binary_methyl_low_)  binary_matrix(i,j)=0; // 0.2
    else                                 binary_matrix(i,j)=-1; // ambiguous
}
```
- 閾值預設 `binary_methyl_high=0.8 / binary_methyl_low=0.2`（`include/core/Config.hpp:33-34`）。
- 閾值只寫進 `binary_matrix`；`raw_matrix` 始終是連續值。
- `methylation.csv` 輸出走連續 raw（`RegionWriter.cpp` write_matrix_csv，4 位小數），**不經閾值**。
- 結論：閾值的下游消費者只有 NHD/JACCARD 距離 + `binary_matrix`-based 度量（如 NME/Epipoly per `PerCpgAsm`）。生產距離（BERNOULLI）與主 `methylation.csv` 均不受其影響 → claim「只在某些輸出二值化」**為真**。

---

## L2 關鍵陷阱：Config 預設 vs ArgParser 生效值（必須釐清，否則 claim 解讀錯）

存在**兩個互相矛盾的「預設 metric」**，需確認誰生效：

| 來源 | 預設 metric | 檔:行 |
|------|------------|-------|
| `Config.hpp` 結構預設 | `BERNOULLI`（連續）| `include/core/Config.hpp:40` |
| `ArgParser` CLI 字串預設 | `"NHD"`（二值）| `include/utils/ArgParser.hpp:86` |

**生效者 = ArgParser**。`include/utils/ArgParser.hpp:163-176`:
```cpp
config.distance_metrics.clear();                  // 163 無條件清空 Config 的 BERNOULLI
for (const auto& s : distance_metric_strs) {      // distance_metric_strs 預設 {"NHD"}
    ...
    config.distance_metrics.push_back(mit->second);
}
if (config.distance_metrics.empty())
    config.distance_metrics.push_back(DistanceMetricType::NHD);  // 175 fallback 也是 NHD
```
→ 經 CLI 啟動時，**Config.hpp 的 BERNOULLI 預設是 dead code**；若使用者**不傳** `--distance-metric`，生效 metric = **NHD（二值）**。

**主距離矩陣取 `distance_metrics_[0]`**（`RegionProcessor.cpp:854`「using first metric for main analysis」+ 419）。

> 🔴 重要結論：**「距離矩陣用連續」並非 binary 內建保證，而是因為生產腳本顯式傳 `BERNOULLI`**。
> - 若有人直接跑 binary **不帶** `--distance-metric` → 拿到的是 **NHD（二值，受 0.8/0.2 閾值影響）**，claim 在此情境**被推翻**。
> - 故 claim 對「ISM binary 預設行為」是 **REFUTED**；對「生產 pipeline 實際產物」是 **SUPPORTED**。合判 = PARTIAL。

---

## L3 correction（校正後正確版）

> **舊理解（claim 原文）**：「距離矩陣用連續機率(非二值);methyl-high/low 閾值只在某些輸出二值化」
>
> **校正版**：
> 1. ISM 距離計算**按所選 metric** 在連續(`raw_matrix`: L1/L2/CORR/BERNOULLI)與二值(`binary_matrix`: NHD/JACCARD)間切換（`DistanceMatrix.cpp:309-345`）。並非無條件連續。
> 2. **生產腳本顯式選 `BERNOULLI`（連續）**（`run_vcf_all_snv.sh:65`、`run_pure_research_round.sh:17`），故**實際論文用的距離矩陣是連續機率**——claim 結論在生產實況成立。
> 3. ⚠ ISM binary 的 **CLI 預設（不傳 `--distance-metric`）卻是 NHD（二值）**（`ArgParser.hpp:86,163-176` 覆蓋 Config.hpp:40 的 BERNOULLI）。論文/報告若寫「ISM 預設用連續距離」會與 binary 預設不符；正確措辭應為「ISM 生產配置選用 BERNOULLI 連續距離」。
> 4. `--methyl-high/--methyl-low`（0.8/0.2）唯一寫入 `binary_matrix`（`RegionProcessor.cpp:1417-1422`）；`raw_matrix` 與 `methylation.csv` 永遠連續，不經閾值——claim 後半「只在某些輸出二值化」**正確**。

---

## uncertain（誠實標註）

- 本輪未列舉**所有** `binary_matrix` 的下游消費者；已確認 NHD/JACCARD 距離 + PerCpgAsm 的 NME/Epipoly 用 binary（per extract `percpg.md`），但是否有其他模組讀 `binary_matrix` 未全 grep。不影響 claim 主結論（距離與主 CSV 連續成立）。
- `distance_use_binary`(`Config.hpp:50`)/`use_binary_matrix`(`DistanceMatrix.hpp`) 旗標在 dispatcher 中未見被讀取（dispatcher 純按 enum metric 選），疑為遺留 dead flag；不改變「按 metric 選矩陣」的事實。
- 未逐一審所有 handoff/HKU pipeline 是否都傳 BERNOULLI；只證實主 `run_vcf_all_snv.sh` 與 `run_pure_research_round.sh`。若存在不傳 `--distance-metric` 的路徑，該路徑會落回 NHD（二值）。
