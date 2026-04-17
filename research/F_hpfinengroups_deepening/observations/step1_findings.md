# F Pilot Step 1 — 觀察紀錄

**Date**: 2026-04-17
**Script**: `scripts/step1_baseline_and_param_sanity.py`
**Data**: `data/step1_{ng_scan,nr_scan,grid_2d,out_of_scope,per_sample}.tsv`

---

## Finding 1 ✅ Baseline 89.1% 成功重現

| Condition | n | TP rate | CI |
|---|---|---|---|
| TO NonLOH + NG≥4 + NR≥80 | 25,744 | **0.8912** | [0.8873, 0.8950] |
| 記載值 (memory) | — | 0.891 | — |

**結論**：與 memory `project_hpfinengroups_subclone_marker.md` 完全一致（差異 <0.001）。
**意涵**：master dataset 自 2026-04 以來未變動，後續分析可信。

---

## Finding 2 ⚠️ NGroups 非單調！（重大發現，之前未記錄）

固定 NR≥80，NonLOH：

| NGroups 值 | n | TP rate |
|---|---|---|
| ==1 | 3,655 | **0.7633** |
| ==2 | 51,559 | **0.6434** ⚠️ |
| ==3 | 64,653 | 0.7742 |
| ==4 | 25,744 | **0.8912** |

**關鍵觀察**：NGroups==2 的 TP rate (0.6434) **低於 NGroups==1** (0.7633)，也低於 overall baseline (0.6699)。
這是非單調關係，**先前所有分析（含 B.1-1/B.1-2/B.1-3）都隱含假設 monotone**。

**可能的生物學解釋**（需後續驗證）：
- NGroups=1：homozygous methylation pattern → 多為簡單 variant（clonal）
- **NGroups=2：biallelic simple pattern → germline ASM 富集（非 somatic！）**
- NGroups=3：intermediate → somatic 訊號開始
- NGroups=4：subclonal heterogeneity → somatic 最強

**意涵**：若 NGroups=2 主要是 germline ASM，則 "HPFineNGroups 單調 = somatic marker" 說法需修正為 "NGroups=4 才是 somatic marker"。

**對其他任務的影響**：
- **B.1-1/B.1-2/B.1-3 結論無變動**（因為他們比的都是 N=4 vs N<4 或 NR-bin matched）
- **但 memory 描述「NGroups 作為 monotone subclone marker」需要精確化**

---

## Finding 3 ⚠️ NR=80 不是最優 threshold

固定 NGroups≥4，NonLOH：

| NR threshold | n | TP rate |
|---|---|---|
| ≥20 | 34,506 | 0.8677 |
| ≥40 | 34,267 | 0.8688 |
| ≥60 | 32,555 | 0.8782 |
| ≥80 | 25,744 | 0.8912 |
| **≥100** | **17,687** | **0.8973** |
| ≥150 | 4,791 | 0.8814 (下降) |

**觀察**：
- TP rate peak 在 NR≥100（0.8973），而非 NR≥80
- NR≥150 反而下降（可能 high-coverage 區域本身有 caller bias）
- 從 NR=80 → 100 coverage loss = 8,057 regions (31%)；TP rate gain = +0.61pp

**Precision-recall trade-off**:
- NR≥80  (n=25,744) 是 **current recommended threshold**（平衡點）
- NR≥100 (n=17,687) 是 **precision-maximal** 選擇
- NR≥60  (n=32,555) 是 **recall-maximal** 可接受選擇（TP rate 仍 >0.87）

**推薦修訂**：
- 若研究目的 = 高置信驗證（如論文 figure）→ 用 NR≥100
- 若研究目的 = coverage 優先（如普查）→ 用 NR≥80（維持現況）

---

## Finding 4 🚨 per-sample TP rate 極度不均（89.1% 是 H2009 主導）

TO NonLOH + NG≥4 + NR≥80 per-sample 分佈：

| sample | n | TP rate | CI |
|---|---|---|---|
| H2009 | **19,979** (77.6%) | **0.9345** | [0.931, 0.938] |
| H1437 | 1,409 | 0.9212 | [0.906, 0.935] |
| HCC1395_DORADO | 637 | 0.9027 | [0.877, 0.925] |
| HCC1395 | 1,173 | 0.8099 | [0.786, 0.832] |
| HCC1937 | 890 | 0.7135 | [0.683, 0.743] |
| HCC1954 | 1,622 | **0.4969** ⚠️ | [0.472, 0.522] |
| COLO829 | 34 | **0.2353** ⚠️ | [0.107, 0.412] |

**危險訊號**：
- **89.1% overall 是 H2009 主導**（H2009 佔 77.6%；overall TP rate weighted by H2009 的 93.45%）
- **HCC1954 TP rate 0.497 幾乎隨機**（有效條件下）
- **COLO829 TP rate 0.235 比 base rate 還低**（但 n=34 太小，CI 寬）
- 若排除 H2009，overall TP rate = (22,943 - 18,671)/(25,744 - 19,979) = 4,272/5,765 = **0.741**（從 89.1% 掉到 74.1%）

**跨樣本可分類**：
- **Tier A**（有效）：H2009, H1437, HCC1395_DORADO（TP rate ≥90%）
- **Tier B**（中等）：HCC1395 (81%), HCC1937 (71%)
- **Tier C**（失效）：HCC1954 (50%), COLO829 (24%)

**對其他任務的影響**：
- **memory 記載「7/7 樣本 NR-bin weighted 後全 POS」在 effect direction 層級仍正確**
  （Δ=TP rate(N=4) - TP rate(N<4)，方向 7/7 正）
- **但 "89.1% TP rate" 無法外推到所有樣本**
- **HCC1954 n=34 COLO829 的 B.1-3 NEGATIVE finding 在這裡獲得進一步佐證**
- **建議任何未來 "TP rate ≥85% 保證" 類聲稱都必須加「僅對 H2009/H1437/HCC1395_DORADO 類高訊號樣本」限定**

---

## Finding 5 🌟 Paired + NG≥4 + NR≥80 TP rate = 99.85%（意外驚人）

| scope | n | TP rate |
|---|---|---|
| Paired (all) + NG≥4 + NR≥80 | **11,801** | **0.9985** |

**觀察**：Paired mode 在相同 NGroups/NR 條件下 TP rate 達 99.85%（n=11,801）。
- 相比 TO NonLOH 的 0.8912 高 10.7pp

**可能原因**：
1. Paired caller 本身 FP 很低（baseline TP rate 高）→ **需要驗證這假設**
2. Paired mode 不需要 LOH filter（因無 self-phasing bias）
3. NR + NGroups 組合對 Paired mode 是強 filter

**立即待驗證**：
- Paired 總體 TP rate（overall baseline，比較 filter gain）
- 若 Paired overall TP rate 本身就 >0.98 → NGroups filter 並無真正 gain

---

## Finding 6 ✅ LOH 區域 NGroups=4 極少（確認 memory 記載）

| scope | n | TP rate |
|---|---|---|
| TO LOH + NG≥4 + NR≥80 | **8** | 0.75 |

memory 記 n=11；實測 n=8（差異可能因不同 NR/truth 過濾）。
**確認**：LOH 區 NGroups=4 基本不存在 → 這個 filter 只在 NonLOH 有意義。

---

## 下一步決策

基於 5 個關鍵觀察：

### 已產生的新問題
1. **Finding 2**：NGroups=2 為何比 NGroups=1 差？是 germline ASM 富集嗎？→ Step 2 調查
2. **Finding 3**：NR=100 比 NR=80 更優的 precision，是否應改推薦？→ Step 2 檢討
3. **Finding 4**：HCC1954 n=1,622 但 TP rate 只 0.497，根因是什麼？→ Step 2 深度分析
4. **Finding 5**：Paired 99.85% 是 filter gain 還是 baseline 就高？→ Step 2 驗證

### Step 2 計畫
1. **NGroups 非單調調查**：NGroups=2 regions 的 AF 分佈、AlleleDelta、LOH、caller_af — 是否富集 germline ASM？
2. **Per-sample 異質性根因**：HCC1954 NG=4 regions 的 FP 為什麼仍多？
3. **Paired baseline 驗證**：Paired 99.85% 是否是 gain
4. **NR=80 vs 100 recommended threshold 檢討**

### 對其他任務的影響紀錄
- **B.1 HPFineNGroups 質疑點**：Finding 2 新增一個未被 B.1-1/B.1-2/B.1-3 覆蓋的面向（NGroups=2 具體性質）
- **memory 需更新**：精確化 "NGroups 作為 subclone marker" 為 "NGroups=4 vs <=3" 二分類而非 ordinal
- **B.2 LOH Subclone AF×Methylation**：不受影響（B.2 比較 AF bin，非 NGroups bin）
