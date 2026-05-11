---
id: ism-kb-05-data-formats-merged-dataset-pitfalls
name: "Merged 合成資料集已知陷阱 (AF 欄位與 HCC1395 phase1_new LOH)"
description: "記錄 2026-04-23 S5 PPT 圖製作發現的兩項 merged 合成資料集陷阱：`AF` 欄位非 caller_af、HCC1395 phase1_new LOH annotation 不完整。涵蓋根因、偵測方法、修正策略、下次使用注意事項。"
status: active
last_verified: 2026-04-23
content_nature: postmortem
doc_type: reference
verified_scope: "AF column mismatch + HCC1395 phase1_new LOH incompleteness pitfalls"
related_ids:
  - ism-kb-05-data-formats-master-dataset-schema
  - ism-kb-05-data-formats-per-region-files
  - ism-kb-02-samples-hcc1395
  - ism-kb-10-research-status-blockers-and-risks
entities:
  - {name: "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz", type: file}
  - {name: "caller_af", type: column}
  - {name: "Potential_LOH", type: column}
  - {name: "LOH_Subtype", type: column}
  - {name: "phase1_new", type: kde_status}
tags: [pitfall, postmortem, data-formats, merged-dataset, AF, caller_af, hcc1395, phase1_new, LOH, annotation-incomplete, 資料陷阱, 欄位誤用]
canonical_paths: [05_data_formats/06_merged_dataset_pitfalls.md]
alias_paths: []
---

# Merged 合成資料集已知陷阱

**一句結論**：`research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz` 合成檔有兩項關鍵陷阱 —— (1) `AF` 欄位**不是** caller 原始 VAF (caller_af) 而是一個殘缺 / 非標準的 delta-like 值；(2) HCC1395 的 `kde_status='phase1_new'` rows 的 LOH annotation 完全殘缺 (Inner 5.7% vs 正確值 58.8%)。使用前必須跳過這兩陷阱。

**適用對象**：任何用此合成檔畫 AF-axis 圖、跑 LOH-aware 分析、或比較跨樣本 HCC1395 vs 其他樣本。

**最後驗證日期**：2026-04-23（本 pitfall 首次在 S5 PPT 製作過程被發現並記錄）

---

## 陷阱 1 · `AF` 欄位不是 caller_af

### 症狀

以 merged 檔 `AF` 欄位畫 scatter plot（Y 軸）時，**非 HCC1395 的 5 個 TO 樣本**幾乎所有點都堆在 AF < 0.1 底部，Y 軸 > 0.4 區域幾乎無點。視覺上「不合理」。

### 量化差異

| 欄位 | 樣本 | mean | median | p75 | p90 | Extreme% |
|:-----|:-----|---:|---:|---:|---:|---:|
| `caller_af` (stale master archive) | HCC1395 TO | 0.545 | 0.479 | 0.813 | — | 合理 0-1 分散 |
| `caller_af` | HCC1395_DORADO TO | 0.552 | 0.485 | — | — | 合理 |
| `AF` (merged) | HCC1395_DORADO TO | **0.035** | **0.012** | **0.052** | 0.106 | **88.8%** |
| `AF` (merged) | HCC1937 TO | 0.007 | 0.003 | 0.027 | 0.059 | **97.2%** |
| `AF` (merged) | HCC1954 TO | 0.024 | 0.013 | 0.033 | 0.061 | **96.1%** |
| `AF` (merged) | H2009 TO | 0.011 | 0.004 | 0.015 | 0.033 | **98.0%** |
| `AF` (merged) | H1437 TO | 0.020 | 0.013 | 0.034 | 0.056 | **98.3%** |

### 根因（待進一步確認）

`AF` 欄位的數值分佈呈現「大部分極小 + 少數中小值」的特徵，**不是** VCF 的 `AF` (tumor allele fraction) 也不是 `caller_af`。可能情況：
- 是 `AlleleDelta`（HP1 vs HP2 的 methylation delta，通常 0-0.3 範圍）的誤拷貝
- 是 filter 後的 pre-selected subset（但如此仍不該分佈至 97% 皆 Extreme）
- 是上游 pipeline 的欄位名衝突（同一 key `AF` 被不同工具寫入）

**需要追查**：`research/ng_kde_rescaling/scripts/` 的合成步驟，確認 `AF` 欄位從哪個 source 取得。

### 正確做法

1. **所有 AF 分析統一使用 `caller_af`（stale master archive 有完整欄位）**：
   - 來源：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz`
   - 欄位：`caller_af` (所有 6 TO 樣本均 non-null，分佈 0-1 合理)
2. 若必須用 merged 檔，**不要用 `AF`**；改用 `caller_af`（merged 若有）或 merge with stale master 補 caller_af
3. 偵測腳本：
   ```python
   # 一鍵驗證：若 p75(AF) < 0.1 則 AF 欄位可疑
   suspect = df.groupby('sample')['AF'].quantile(0.75)
   print(suspect[suspect < 0.1])  # 列出可疑樣本
   ```

### 影響範圍

- `docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md`（v2）:
  使用 `AF_class` 作切分，`AF_class` 來源是否對也需檢查
- Phase 3 Synthesis `s1s7_per_sample.tsv`: 有使用 scheme 含 AF 條件，若 `AF_class` 推導自錯誤 `AF`，scheme 定義可能偏誤
- 任何用此合成檔做 `AF` 視覺化的腳本

---

## 陷阱 2 · HCC1395 `kde_status='phase1_new'` LOH annotation 不完整

### 症狀

HCC1395 在 merged 檔中的 `Potential_LOH=True` 比例為 **5.7%** (2,303 / 40,115)，遠低於同細胞株的 HCC1395_DORADO 60.0% (24,241 / 40,428) 與其他樣本 22-61%。

### 根因

`kde_status='phase1_new'` 標記的 run 是 2026-04-21 `--germline-hp-only` Phase 1 實驗的**輕量**執行，**未跑完整 LOH_Bed_Overlap / LOH_Subtype annotation pipeline**。LOH_Subtype 分佈對比：

| LOH_Subtype | phase1_new (HCC1395) | archive_rerun (HCC1395_DORADO) |
|:------------|:---:|:---:|
| None | 37,812 | 16,187 |
| LOH_Weak | 1,491 | 4,943 |
| **LOH_Noise** | **431** | **12,534** |
| LOH_Strong | 357 | 5,266 |
| LOH_Subclone | 24 | 1,498 |

`LOH_Noise` 差 29 倍 — 明顯是 annotation 殘缺而非 biology。

### 偵測方法

```python
# 對每樣本比較 Potential_LOH 比例；偏差大者警示
loh_pct = df.groupby('sample')['Potential_LOH'].mean() * 100
median_pct = loh_pct.median()
for s, pct in loh_pct.items():
    if abs(pct - median_pct) > 30:  # >30pp 偏差
        print(f"⚠ {s} LOH {pct:.1f}% deviates from median {median_pct:.1f}%")
```

### 正確做法

1. 遇到 `kde_status='phase1_new'`：**從 stale master archive 替換**（下面 stale archive 路徑）
2. 2026-04-30 (下週 P0) master × 兩 flag × 7 樣本重跑完成後可直接用新 master

### Stale Master Archive 路徑（HCC1395 完整 LOH annotation 來源）

```
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
    20260330_loh_round1_cross_sample_audit_post_to_hp_fix/
    all_region_rows.tsv.gz
```

- 欄位：`truth_label` (`TP`/`FP` 字串)、`caller_af`、`Potential_LOH`、`LOH_Subtype`、`NumReads`、`caller_af`、`mode` (`to`/`paired`)、`sample`
- 可用樣本（TO mode）：HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829 共 7 個
- Coverage_Multiple 是 stale 75× 校準的（需要 per-sample KDE baseline 重算）

---

## 修正範例：S5 PPT 圖的正確讀取模式

本 pitfall 發現過程的完整修正腳本：
`/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/20260423_s5_loh_af_cn_scatter.py`

核心邏輯：

```python
# 所有 6 TO 樣本統一從 stale master archive 讀
SAMPLE_KDE_BX = {
    "HCC1395": 54.0,          # SEQC2 neutral median
    "HCC1395_DORADO": 53.0,
    "HCC1937": 91.0,
    "HCC1954": 61.0,
    "H2009": 79.0,
    "H1437": 71.0,
}

arch = pd.read_csv(STALE_MASTER, sep="\t")
arch_to = arch[arch["mode"] == "to"]
arch_to["bx"] = arch_to["sample"].map(SAMPLE_KDE_BX)
arch_to["Coverage_Multiple"] = arch_to["NumReads"] / arch_to["bx"]  # KDE-equiv
arch_to["AF"] = arch_to["caller_af"]   # 真 VAF
arch_to["tp_label"] = (arch_to["truth_label"] == "TP").astype(int)
arch_to["loh_side"] = arch_to["Potential_LOH"].map({True: "Inner", False: "Outer"})
```

結果：6/6 樣本 Inner TP rate > baseline（Δ+0.6 至 +3.1pp），LOH 生物學預期成立。

---

## 時序紀錄

| 日期 | 事件 |
|:-----|:-----|
| 2026-03-30 | stale master `all_region_rows.tsv.gz` 產出，含正確 `caller_af` + 完整 LOH |
| 2026-04-21 | HCC1395 `--germline-hp-only` Phase 1 run：只為 flag 驗證，未跑完整 LOH annotation → 產生 `kde_status='phase1_new'` 殘缺資料 |
| 2026-04-22 | `merged_7samples_paired_full_plus_hcc1395_to.tsv.gz` 合成，納入 phase1_new HCC1395（LOH 殘缺）+ archive 其他樣本。合成過程中 `AF` 欄位來源不明確，與 `caller_af` 不一致 |
| **2026-04-23 (上午)** | S5 PPT 圖用 merged 檔畫出後發現 HCC1395 Inner n=2,303（vs DORADO 24,241） → 發現陷阱 2 |
| **2026-04-23 (下午)** | S5 圖用 archive 替換 HCC1395 後其他樣本 AF 仍集中 <0.1 → 發現陷阱 1 |
| **2026-04-23 (本記錄)** | 兩陷阱根因確認，全部 6 樣本統一從 archive 讀取修正 |

---

## 未來修正計畫

1. **下週 P0 · Archive TO 6 樣本 ISM 重跑**：重跑後新 master 將含一致的 `caller_af` + 完整 LOH + KDE baseline，取代 merged 合成檔
2. **合成腳本審計**：`research/ng_kde_rescaling/scripts/` 的欄位映射需確認 `AF` 應改名或修正來源
3. **validation script**：加入 `scripts/validate_merged_dataset.py` 自動檢查 `AF` 分佈合理性與 `Potential_LOH` 樣本差異

---

**相關檔案**：
- PPT 證據：`docs/presentations/validated/2026/04/20260423_研究週報_LOH_constrained_phasing_pivot/figures/fig_s5_loh_inner_outer_af_cn_per_sample.png`
- 修正腳本：`scripts/analysis/20260423_s5_loh_af_cn_scatter.py`
- S5 speaker notes：`build_pptx.py::s5_loh_af_cn_initial()`
