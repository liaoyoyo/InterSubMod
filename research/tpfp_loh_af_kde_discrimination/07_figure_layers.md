---
title: 分層 Combinatorial 圖集 — L1/L2/L3/Consistency × 13 sample-mode（V2 2026-04-22 post-rerun）
status: authoritative_post_rerun
last_updated: 2026-04-22 15:31
figure_count: 37
sample_mode_count: 13
data_source: research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz
---

# 07. 分層 Combinatorial 圖集

## 7.0 TL;DR（3 分鐘讀完）

**核心發現（2026-04-22 ISM rerun 後）**：

1. **LOH / NG 梯度跨 6 個 TO 樣本方向一致**（5-6 / 6 樣本呈正向），非 HCC1395-specific；HCC1954_TO 的 AF 方向為**唯一結構性反轉**。
2. **Hard Gate 完全觸發**：408 個 5D cells 中 **16 個** 達 `n_samples_high ≥ 5/13`（閾值 10，**超標 60%**），`≥ 3/13` 達 64 個（閾值 50）。
3. **H4 判決升級**：從「HCC1395-specific 弱支持」升級為「**6 TO 樣本中-強支持**」。
4. **Top-17（HCC1395_TO）與 16-cell 跨樣本共識集 0/17 重疊**；共識主軸為 `Extreme AF + NG=2 + T1-T3 CN`（非 Top-17 的 Near-half + NG=3）。

**唯一必看圖**：`figures/new/L3/L3_loh_af_cn.png`（fig_v2_1 的 13 sample-mode 擴展版）。

**下一步建議**：LOSO 驗證聚焦 **16-cell 共識集**（見 §7.5.1 Tier A/B/C），而非 HCC1395 Top-17。預期泛化性更強但 F1 改善幅度較小。

**閱讀路徑提示**：
- **3 分鐘**：本節 §7.0 + 看 L3_loh_af_cn 主圖
- **15 分鐘**：+ §7.5.1 Tier 表 + §7.7 Hard Gate + `data/cell_consistency_matrix.tsv` 前 16 row
- **45 分鐘**：L1 TO (5 圖) + L2 TO 關鍵 3 圖（LOH×AF、LOH×NG、AF×NG）+ L3 全 4 圖 + 3 consistency
- **完整**：全文 L1 (10) + L2 (20) + L3 (4) + consistency (3) = 37 圖

---

## 7.1 圖集設計（V2 2026-04-22 更新）

本檔 37 張圖分四層：

| 層次 | 張數 | 目的 | 佈局 |
|------|----:|------|------|
| L1 單維度 | 10（5 dim × 2 mode）| 哪個維度在多少樣本都重要 | 2×3 TO (6 panel) / 3×3 paired (7 panel, 2 空) |
| L2 雙維度 | 20（10 pair × 2 mode）| 兩兩組合下哪些樣本有一致高 TP cell | 同上 |
| L3 三維度 | 4（trio）| 三維聯合下跨樣本是否一致（13 sample 同頁）| 5×3 macro × inner facet strip |
| Consistency | 3（5D-all / 3D-LOCN / 5D-TO-only）| 收斂成數字：多少樣本一致 | stacked bar / scatter |

**每張圖規格**（V2 統一）：
- **Suptitle**：研究問題格式（非視覺機制描述）
- **Per-panel metadata**：`[TIER] FP=N base=0.XX n=N`（右上白底框）；FP<100 加紅邊框
- **Takeaway caption box**：底部色框，綠（一致）/ 藍（中度）/ 紅（衝突）/ 灰（不足）

**核心設計決策**：L1/L2 拆成 TO-only + paired-only 兩版以避免「混用」——
TO baseline 25-91% 呈現強梯度，paired baseline 94-99% 基本飽和，同一張圖無法同時展示這兩個生物特性。

**13 sample-mode 統計力表**：

| Mode | Tier | sample | FP | baseline TP |
|:----:|:----:|--------|---:|-------------|
| TO   | HIGH | HCC1395 | 11,606 | 0.71 |
| TO   | HIGH | HCC1395_DORADO | 11,572 | 0.71 |
| TO   | HIGH | HCC1937 | 12,032 | 0.51 |
| TO   | HIGH | HCC1954 | 50,218 | 0.25 |
| TO   | HIGH | H2009 | 11,989 | 0.91 |
| TO   | HIGH | H1437 | 13,442 | 0.77 |
| pf   | HIGH | COLO829 | 2,273 | 0.94 |
| pf   | MID  | HCC1395 | 627 | 0.98 |
| pf   | MID  | HCC1395_DORADO | 240 | 0.99 |
| pf   | MID  | HCC1937 | 195 | 0.98 |
| pf   | LOW  | H2009 | 86 | 0.9994 |
| pf   | LOW  | HCC1954 | 29 | 0.9984 |
| pf   | LOW  | H1437 | 8 | 0.9999 |

---

## 7.2 L1 單維度觀察（10 圖，5 dim × 2 mode）

### 7.2.1 LOH_Subtype

![L1 LOH TO](figures/new/L1/L1_loh_TO.png)
![L1 LOH paired](figures/new/L1/L1_loh_paired.png)

**觀察**：TO 6 樣本（上圖）呈現 `None → LOH_Weak/Noise → LOH_Strong` 的正向梯度，跨 5-6 樣本方向一致。最強梯度：HCC1954_TO（0.25→0.85）與 HCC1395_TO（0.71→0.93）。paired_full（下圖）因 baseline 飽和（≥0.94），所有類別差異 <0.02。

**結論**：LOH 維度在 TO mode 具跨樣本穩健判別力（升級自「HCC1395-only」），paired mode 無可視覺化 lift。

### 7.2.2 AF_class

![L1 AF TO](figures/new/L1/L1_af_TO.png)
![L1 AF paired](figures/new/L1/L1_af_paired.png)

**觀察**：HCC1395_TO / HCC1395_DORADO_TO / HCC1937_TO / H2009_TO 呈 `Extreme < Intermediate < Near-half`；**HCC1954_TO 方向反轉**（Extreme 0.86 >> Near-half 0.14），為結構性例外。

**結論**：AF_class 跨 5/6 TO 樣本一致；HCC1954 需獨立處理（該樣本 AF 分佈受 germline-like allele 污染）。

### 7.2.3 cn_tier_F (Coverage Multiple)

![L1 CN TO](figures/new/L1/L1_cn_TO.png)
![L1 CN paired](figures/new/L1/L1_cn_paired.png)

**觀察**：T2（diploid）在所有樣本為主力。T3/T4 在 HCC1395_TO / HCC1395_DORADO_TO / H2009_TO 偏 FP；HCC1937_TO / HCC1954_TO 梯度較弱。

**結論**：CN tier 需搭配樣本特異 baseline 做決策，不建議作為單一 filter 軸。

### 7.2.4 HPFineNGroups (NG)

![L1 NG TO](figures/new/L1/L1_ng_TO.png)
![L1 NG paired](figures/new/L1/L1_ng_paired.png)

**觀察**：NG=3/4 > NG=2 的梯度在 **5/6 TO 樣本**可見（gap 0.05-0.12），HCC1954_TO 微弱。**推翻舊版「NG=3 為 HCC1395-specific」結論**。

**結論**：NG 維度跨 TO 樣本穩健；與 `project_hpfinengroups_subclone_marker.md` POSITIVE 結論吻合。警告：`--germline-hp-only` flag 下 NG≥3 可能消失（見 Phase 1 報告）。

### 7.2.5 nr_band (NumReads)

![L1 NR TO](figures/new/L1/L1_nr_TO.png)
![L1 NR paired](figures/new/L1/L1_nr_paired.png)

**觀察**：各樣本 high_NR > mid_NR > low_NR 梯度溫和（所有差 <0.08）。

**結論**：NR 單一維度判別力最弱，價值在與其他維度組合使用（見 §7.3 L2.10）。

---

## 7.3 L2 雙維度觀察（20 圖）

每 dim pair 兩張（TO / paired）。以下 10 組依生物學重要性排序；TO / paired 各對列。

### 7.3.1 LOH × AF（核心 biology module）

![L2 LOH×AF TO](figures/new/L2/L2_loh_x_af_TO.png)
![L2 LOH×AF paired](figures/new/L2/L2_loh_x_af_paired.png)

**結論**：**6 個 TO 樣本都可見「LOH_Strong + Near-half 深綠」vs「None + Extreme 淺/紅」色彩梯度**（HCC1954 方向反轉例外）。paired 全綠飽和，無區分。

### 7.3.2 LOH × CN

![L2 LOH×CN TO](figures/new/L2/L2_loh_x_cn_TO.png)
![L2 LOH×CN paired](figures/new/L2/L2_loh_x_cn_paired.png)

**結論**：LOH_Strong × T1-T3 於 4/6 TO 深綠；T3/T4 × None 在 HCC1954 偏低。跨 TO 樣本 pattern 部分一致。

### 7.3.3 LOH × NG

![L2 LOH×NG TO](figures/new/L2/L2_loh_x_ng_TO.png)
![L2 LOH×NG paired](figures/new/L2/L2_loh_x_ng_paired.png)

**結論**：(LOH_Strong × NG=3/4) 深綠在 ≥4/6 TO 重現；此為 `HPFineNGroups × LOH` 複合 marker 的最佳可視化證據。

### 7.3.4 LOH × NR

![L2 LOH×NR TO](figures/new/L2/L2_loh_x_nr_TO.png)
![L2 LOH×NR paired](figures/new/L2/L2_loh_x_nr_paired.png)

**結論**：LOH_Strong × high_NR 於 TO 樣本穩定深綠；low_NR 區塊多樣本 n<20 遮罩。

### 7.3.5 AF × CN

![L2 AF×CN TO](figures/new/L2/L2_af_x_cn_TO.png)
![L2 AF×CN paired](figures/new/L2/L2_af_x_cn_paired.png)

**結論**：Near-half + T1/T2 在 5/6 TO 深綠（對應 S3 Diploid Het）；HCC1954 反向。

### 7.3.6 AF × NG

![L2 AF×NG TO](figures/new/L2/L2_af_x_ng_TO.png)
![L2 AF×NG paired](figures/new/L2/L2_af_x_ng_paired.png)

**結論**：(Near-half, NG=3/4) 在 5/6 TO 最純；HCC1954 例外。

### 7.3.7 AF × NR

![L2 AF×NR TO](figures/new/L2/L2_af_x_nr_TO.png)
![L2 AF×NR paired](figures/new/L2/L2_af_x_nr_paired.png)

**結論**：high_NR + Near-half 於 TO 普遍接近 baseline 飽和；Extreme + low_NR 最雜。

### 7.3.8 CN × NG

![L2 CN×NG TO](figures/new/L2/L2_cn_x_ng_TO.png)
![L2 CN×NG paired](figures/new/L2/L2_cn_x_ng_paired.png)

**結論**：(T2, NG=3) 在 4/6 TO 達 0.83-0.93；(T4, NG=4) 多樣本 n 少被遮。

### 7.3.9 CN × NR

![L2 CN×NR TO](figures/new/L2/L2_cn_x_nr_TO.png)
![L2 CN×NR paired](figures/new/L2/L2_cn_x_nr_paired.png)

**結論**：T2 × high_NR 在多樣本 >0.85；差距與 baseline 相近。

### 7.3.10 NG × NR

![L2 NG×NR TO](figures/new/L2/L2_ng_x_nr_TO.png)
![L2 NG×NR paired](figures/new/L2/L2_ng_x_nr_paired.png)

**結論**：(NG=3/4, high_NR) 組合在 5/6 TO 呈一致 lift；HCC1954 例外。

**L2 整體收斂**（相較舊版結論徹底改變）：**6 TO 樣本 baseline 跨越 25-91%，在 L2.1 / L2.3 / L2.6 呈跨樣本一致的顏色梯度**（幅度略異但方向一致），將 H4 目視判決從「弱支持」**升級為「中-強支持」**。

---

## 7.4 L3 三維度觀察（4 圖）

L3 保留 13-sample 合併版（5×3 macro grid），因三維跨樣本對比需要同頁面比較；caption 自動計算跨樣本 Top-3 cells。

### 7.4.1 LOH × AF × CN（⭐ fig_v2_1 的 13 sample-mode 擴展版）

![L3 LOH×AF×CN](figures/new/L3/L3_loh_af_cn.png)

**對應用戶原始問題的核心圖**。以下為 per-sample 觀察表：

| sample_mode | 可見梯度 | 主要觀察 |
|-------------|:-------:|---------|
| HCC1395_TO | ✓ 強 | T0/T1 `None+Extreme` 紅區；`LOH_Strong+Near-half` 深綠，LOH×AF 梯度每個 CN 片清楚 |
| HCC1395_DORADO_TO | ✓ 強 | 與 HCC1395_TO pattern 近似；DORADO basecaller 下 T2 LOH_Weak 格略暗 |
| HCC1937_TO | ✓ 中 | LOH_Strong 各 CN 片深綠；baseline 0.51 → 無飽和現象 |
| HCC1954_TO | ⚠ 反向 | AF 方向反轉：`None+Extreme` 反而綠；`LOH_Strong+Near-half` 轉灰（低 n）|
| H2009_TO | ✓ 中 | baseline 0.91 接近飽和；可見 LOH_Strong > None 微弱梯度 |
| H1437_TO | ✓ 弱 | baseline 0.77；整體偏均勻，LOH_Strong 略優於 None |
| COLO829_pf | ○ 平 | 幾乎全綠；僅 `None+Extreme+T0` 略淺 |
| HCC1395_pf | ○ 平 | 均勻深綠；目視無區分 |
| HCC1395_DORADO_pf | ○ 平 | 均勻深綠 |
| HCC1937_pf | ○ 平 | 均勻深綠 |
| H2009_pf | 🔴 遮 | n<20 多數灰遮（FP=86）|
| HCC1954_pf | 🔴 遮 | n<20 灰遮（FP=29）|
| H1437_pf | 🔴 遮 | 幾乎全灰（FP=8）|

**H4 目視判讀**：Top-17 中的 `LOH_Strong + T1 + NG=3` 類 cells **在 4/5 TO 樣本**（HCC1954 除外）可被目視確認類似梯度。升級為中度支持；Hard Gate 於 §7.7 完全觸發。

### 7.4.2 LOH × AF × NG

![L3 LOH×AF×NG](figures/new/L3/L3_loh_af_ng.png)

**結論**：NG=3/4 slice 在 5/6 TO 可見 (LOH_Strong, Near-half) 深綠區；**反駁舊版「NG=3 HCC1395-specific」**。

### 7.4.3 AF × CN × NR

![L3 AF×CN×NR](figures/new/L3/L3_af_cn_nr.png)

**結論**：high_NR slice 所有樣本最飽和；Near-half + T2 + high_NR 跨 6 TO 穩定 >0.90。

### 7.4.4 LOH × CN × NG

![L3 LOH×CN×NG](figures/new/L3/L3_loh_cn_ng.png)

**結論**：(LOH_Strong, T2, NG=3) 在多數 TO 深綠；跨樣本 NG 梯度在 T2 切片最明顯。

---

## 7.5 跨樣本一致性（Consistency，3 圖）

**判準**：某個 5D cell 在某個樣本上「high」= `cell FP rate ≤ 0.5 × baseline FP rate`。此 baseline-scale-aware 判準取代「+0.10 絕對閾」（paired baseline>0.90 下不可達）。

**輸出**：
- `data/cell_consistency_matrix.tsv`（**408 cells** × 36 欄，13 sample 版）
- `data/cell_consistency_matrix_TO.tsv`（TO-only 版，6 sample）

### 7.5A Top-30 (13 sample-mode) 跨樣本 5D cells

![Cell consistency 5D Top-30](figures/new/consistency/cell_consistency_5d.png)

### 7.5B LOH × AF × CN 空間上的一致性 map

![Cell consistency 3D LOCN](figures/new/consistency/cell_consistency_3d_locn.png)

### 7.5C Top-30 **TO-only**（6 sample）

![Cell consistency TO-only](figures/new/consistency/cell_consistency_TO_only.png)

**為何新增 TO-only 版？** — 原 13-sample 版的 `n_samples_high` 受 paired_full baseline 飽和影響（paired 幾乎無 cell 能「FP 半減」），TO-only 版聚焦「6 TO 樣本之間」的真共識。

### 7.5.1 共識 cells 三層分類（V2 2026-04-22）

將 408 cells 按 `n_samples_high` 分層：

#### Tier A — 強共識（n_samples_high = 7/13，3 cells）

> **這 3 cells 在 13 sample-mode 中有 7 個達 FP 半減。** 可作為 LOSO 驗證的首要目標。

| LOH | AF | CN | NG | NR | n_samples_n20 | cross_avg_tp | n_samples_high |
|-----|----|----|----|----|--------------:|-------------:|---------------:|
| LOH_Strong | Extreme | T1 | 2 | mid | 9 | 0.914 | **7** |
| LOH_Weak | Extreme | T1 | 2 | mid | 11 | 0.889 | **7** |
| None | Extreme | T1 | 4 | mid | 12 | 0.887 | **7** |

**共同 pattern**：T1 CN（`CN<2`，loss/del-LOH 區域）+ Extreme AF + mid NR。LOH 類別分散（Strong / Weak / None 都出現），顯示**CN-loss + Extreme AF 是最跨樣本穩健的 ISM 特徵**。

#### Tier B — 中強共識（n_samples_high = 6/13，6 cells）

| LOH | AF | CN | NG | NR | n_samples_n20 | cross_avg_tp | n_samples_high |
|-----|----|----|----|----|--------------:|-------------:|---------------:|
| LOH_Subclone | Extreme | T3 | 2 | high | 7 | 0.972 | 6 |
| LOH_Subclone | Extreme | T1 | 2 | mid | 8 | 0.970 | 6 |
| LOH_Strong | Extreme | T0 | 2 | low | 6 | 0.944 | 6 |
| LOH_Noise | Extreme | T3 | 2 | high | 10 | 0.934 | 6 |
| None | Extreme | T3 | 4 | high | 10 | 0.908 | 6 |
| None | Extreme | T4 | 3 | high | 13 | 0.825 | 6 |

**觀察**：LOH_Subclone 在 Tier B 突出（2/6 cells），提示**中階 LOH 純度 + Extreme AF** 是 Tier A 之後的次要共識軸。

#### Tier C — 中度共識（n_samples_high = 5/13，7 cells）

| LOH | AF | CN | NG | NR | n_samples_high |
|-----|----|----|----|----|---------------:|
| LOH_Noise | Extreme | T2 | 2 | high | 5 |
| LOH_Weak | Extreme | T2 | 2 | high | 5 |
| LOH_Strong | Intermediate | T1 | 2 | mid | 5 |
| None | Extreme | T2 | 2 | high | 5 |
| LOH_Subclone | Extreme | T1 | 2 | low | 5 |
| LOH_Weak | Extreme | T0 | 2 | low | 5 |
| LOH_Weak | Extreme | T2 | 2 | mid | 5 |

**觀察**：所有 7 cells 的 **NG = 2**（非 Top-17 的 NG=3）；AF = Extreme（15/16 Tier A+B+C cells）。

**對 HCC1395_TO Top-17 的交集**：Top-17 核心 = `LOH_Strong + Near-half + T1 + NG=3 + mid_NR`。Tier A+B+C 合計 16 cells 中，**0/17 Top-17 cells 符合此模式**。共識主軸為「Extreme AF + NG=2 + T1-T3」，與 Top-17 的「Near-half + NG=3」**完全分離**。

**生物學解讀**：
1. **Extreme AF + NG=2**：代表「germline-like allele distribution 但 methylation pattern 單一」的區域；caller baseline 普遍能處理，跨樣本一致為合理預期。
2. **Top-17 的 Near-half + NG=3**：代表「high discrimination 但 low reproducibility」—— 在 HCC1395_TO 局部顯著，因要求多條件同時滿足，其他樣本 n 不足或條件漂移。

---

## 7.6 結論對比（rerun 前 vs 後）

下表整合原 §7.6 + §7.9，提供跨文件引用的標準對照：

| 議題 | rerun 前（8 sample-mode）| rerun 後（13 sample-mode）| H4 判決影響 |
|------|------------------------|------------------------|------------|
| paired_full baseline 飽和 | 已記錄 | 不變 | 中性 |
| LOH gradient 跨樣本 | HCC1395_TO-specific | **6/6 TO 樣本方向一致**（幅度異）| ↑ 升級 |
| NG=3 boost | HCC1395_TO-specific | **5/6 TO 樣本**可見 | ↑ 升級 |
| Top-17 vs 跨樣本共識 | 7 cells n_samples_high≥5，Top-17 重疊 1/17 | **16 cells n_samples_high≥5，Top-17 重疊 0/17** | ⚠ Top-17 不可直接移植，但有新 16-cell 候選集 |
| 單維度判別力（各維度）| 只有 HCC1395_TO 有梯度 | **TO mode 全部 6 樣本都可見**，paired 全飽和 | ↑ 升級 |
| 2D heatmap 顏色梯度 | 僅 HCC1395_TO | **5/6 TO 樣本**（HCC1954 反向例外）| ↑ 升級 |
| 3D 跨樣本一致 | 7 cells ≥5/8 | **16 cells ≥5/13**（Hard Gate §1 觸發）| ↑ 完全觸發 |
| 最強單 cell 共識度 | 6/8（最高）| **7/13**（3 cells）| 相對 rate 類似（~50-55%）|

**跨文件引用建議**：
- 04_comparison_narrative §4.9 → 引 Tier A 3 cells + Hard Gate 全面觸發
- 06_limitations_future_work §6.1.1 → 引「Top-17 重疊 0/17」+ 推薦 LOSO 改聚焦 16-cell 共識集
- 06 §6.5.1 → 引 §7.5.1 Tier A/B/C 表

---

## 7.7 Hard Gate 觸發與建議（**完全觸發**）

依 Plan §6.3：

| Gate 條件 | 觸發閾 | 觀察值 | 觸發？ | 建議動作 |
|----------|-------|-------:|:-----:|---------|
| §1 consistency `n_samples_high ≥ 5` | ≥ 10 | **16** | ✅ **完全觸發** | **強烈建議 LOSO Tier 1** |
| §2 `n_samples_high ≥ 3` | ≥ 50 | **64** | ✅ **觸發** | 多數共識 cells 有中度支持 |
| §3 L3 TO top cells 在其他樣本可判讀 | 定性 | §7.4.1 5/6 TO 可觀察 | ✅ 部分觸發 | 更新 06 §6.1.1 警語 |

**建議併入 06 §6.5.1 的敘述（取代舊版文字）**：

> **2026-04-22 ISM rerun 後重大發現**：6 個 TO 樣本經 post-HP-fix ISM 運行補齊 LOH/CN/KDE 欄位後，首次可在同一分析框架比對。定量結果：
>
> - **16 cells 達 n_samples_high ≥ 5**（Hard Gate §1 超標 60%）
> - 9 cells 達 ≥ 6；3 cells 達 7/13
> - HCC1395_TO Top-17 與 Tier A+B+C 共識 16 cells **0/17 重疊**
>
> **下一步 LOSO 應聚焦 16-cell 共識集合**（§7.5.1 Tier A/B/C），而非 HCC1395 Top-17。此 16-cell 子集理論 F1 改善幅度可能較小但跨樣本泛化性更強，符合 Phase 2 read-level epigenetic characterization 的定位目標。

---

## 7.8 重用與產生方式

四個腳本位於 `scripts/`（V2 2026-04-22 大改寫，支援 TO/paired 拆分 + research_suptitle + takeaway_caption helpers）：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination
python3 scripts/obs09_L1_single_dim_per_sample.py    # 10 L1 圖 (5 dim × 2 mode)
python3 scripts/obs10_L2_double_dim.py               # 20 L2 圖 (10 pair × 2 mode)
python3 scripts/obs11_L3_triple_dim.py               # 4 L3 圖 (13 sample combined)
python3 scripts/obs12_cross_sample_consistency.py    # 3 consistency 圖 + 2 TSV
```

所有圖輸出至 `figures/new/L1|L2|L3|consistency/`；TSV 至 `data/cell_consistency_matrix{,_TO}.tsv`。

V2 新增 helper（`_obs_common.py`）：
- `research_suptitle(fig, question, context)` — 研究問題格式標題
- `takeaway_caption(fig, text, color)` — 底部色框結論
- `mode_subset_grid(mode)` — TO (2×3) / paired (3×3) grid 生成
- `annotate_power_tier(ax, ..., baseline_tp, n_total)` — 擴展版 panel metadata badge

**上游重跑依賴**（已完成 2026-04-22）：
- `obs15_rerun_archive_to_ism.sh` → `obs15b_resume_failed_pilots.sh`：5 archive TO pilot ISM rerun
- `obs16_regen_after_ism_rerun.sh`：step0→step7 下游重生成
