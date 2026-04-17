<!--
建立時間: 2026-04-17 23:20
目標: 定義 evidence_ledger.jsonl 結論的可信度分級 (tier 1-5) 與評分 rubric
處理範圍: research/autoresearch/evidence_ledger.jsonl、hypothesis_queue.json
關聯檔案:
  - docs/standards/research_terminology.json
  - scripts/analysis/score_evidence_tier.py
-->

# Evidence Tier 分級規範

本規範定義每筆 `evidence_ledger.jsonl` 條目的可信度 tier（1-5），供 landscape 報告、MEMORY.md、週報摘要統一引用。

## 設計原則

- **tier 越高，結論越不易翻盤**
- **tier ≤ 2 的 positive 結論不可寫入 MEMORY.md 的 Concluded 區段**
- **tier 升級需要補實驗**，不得以「再分析既有資料」直接升級

---

## 五級分級

| Tier | 名稱 | 準則（必須全部滿足） | 範例 |
|------|------|-------------------|------|
| **1** | `pilot_single_sample` | 單一樣本 pilot；未控 confound；無交叉驗證 | H011 QS≥50 rescue 首次 pilot on HCC1395_5kHz |
| **2** | `multi_sample_no_confound_control` | ≥2 樣本但未做 within-group OLS 或 AF-bin 交叉；或有多方法但無 permutation | 7 樣本觀察但只報 pooled 結果 |
| **3** | `confound_controlled_single_method` | 單一方法 + 至少一種 confound control（within-group OLS / AF-bin / nested CV） | O12 內做了 within-group OLS 並報 delta=0 |
| **4** | `multi_method_cross_validated` | ≥2 獨立方法（eg. RF + logistic regression）都給一致結論；含 permutation test p<0.05 | HPFineNGroups subclone marker (positive 在多 distance metric 下一致) |
| **5** | `independent_reproduction_or_literature_aligned` | 獨立樣本/獨立團隊重現；或與已發表文獻機制吻合 | Self-phasing circular dependency（與 LongPhase 文獻吻合） |

---

## 評分旗標（additive，填在 evidence_ledger 的 `tier_flags` 欄位）

以下每個 flag 代表一項已通過的 confound 控制或驗證：

| Flag | 意義 | Tier 最低要求 |
|------|------|--------------|
| `within_group_ols` | 用 within-group OLS 殘差化（而非 pooled OLS） | ≥3 |
| `af_bin_stratified` | 結果在 AF binning 下仍成立（抗 AF confound） | ≥3 |
| `permutation_tested` | permutation test p<0.05 | ≥3 |
| `nested_cv` | nested cross-validation | ≥3 |
| `multi_method_agreement` | ≥2 方法得一致結論 | ≥4 |
| `multi_sample_consistent` | ≥4 樣本方向一致（7 sample 時需 ≥5/7） | ≥4 |
| `literature_aligned` | 與已發表文獻機制吻合 | ≥5 |
| `independent_reproduction` | 獨立資料/獨立團隊重現 | 5 |

---

## Schema 擴充

新加入兩個欄位：

```json
{
  "tier": 3,                          // 1-5, 本條目最終 tier
  "tier_flags": [                      // array of strings, 對應上表
    "within_group_ols",
    "af_bin_stratified",
    "permutation_tested"
  ],
  "confidence_cap": 3,                 // tier_flags 推算出的最高可達 tier
  "tier_rationale": "單方法+三項 confound control，未跨方法"
}
```

**舊欄位 `tier_used` 保留**（記錄 pilot 執行時的運算 tier，與結論 tier 分開）。

---

## 自動化

- `scripts/analysis/score_evidence_tier.py` 讀 evidence_ledger 條目，根據 `tier_flags` 推算 `confidence_cap`
- 未來 hook 可在提交前警告「tier 宣告 > confidence_cap」

---

## MEMORY.md 使用規則

- **Active Research 區段**：允許 tier 1-2 條目（研究進行中）
- **Concluded 區段**：僅允許 tier ≥3 的結論
- **超過 1 年未升級 tier ≥3 的 Active 條目**：應移入 `Pending` 或標記「需升級驗證」

---

## 版本

- v1.0 (2026-04-17)：首版，引入 5 tier + 8 flag
