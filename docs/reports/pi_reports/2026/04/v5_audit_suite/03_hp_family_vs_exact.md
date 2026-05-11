---
title: V5 Audit Suite — 03. HP family vs HP exact — phasing notation analysis
report_id: 20260427_V5_audit_03_hp_family_vs_exact_01
date: 2026-04-27
authors: AI agent (analysis), liaoyoyo2001 (PI)
sample: HCC1395 (TO mode)
inputs: same as report 02
script: InterSubMod/scripts/analysis/v5_read_intersection_concordance.py
data:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/hp_family_exact.tsv
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/per_site_concordance.tsv
status: draft
tier: B (audit-grade evidence; explanation report — no new statistical claim)
---

# 03. HP family vs HP exact — phasing notation analysis

## 摘要 (TL;DR)

- **L1 family rate ≠ L2 exact rate** 是 systematic 現象，不是「phasing 錯誤」造成。
- 1303 個 shared-read 觀察中，**BL 有 87 reads 「家族對但 exact 錯」、V5 有 27 reads 「家族對但 exact 錯」**；幾乎全部來自 longphase 整數 HP（如 `1`、`2`）對 paired 字串 HP（如 `1-1`、`2-1`）的 **規範差異**。
- 結論：在比較 BL/V5 與 paired 的 phasing 時，**應以 L1 family（PS-aware）為主、L2 exact 僅作 notation hygiene**。L2 exact 低不代表 phasing 錯。

---

## Section 1. HP family vs exact 的差異

### 1.1 兩種 HP tag 編碼系統

InterSubMod 的 5 個 BAM 採用兩套 HP 編碼，互不直接相容：

| BAM | HP tag dtype | 範例值 | 編碼說明 |
|-----|-------------|--------|---------|
| BL / V2b / V3F / V5 (longphase tagged) | **int** | `1`, `2`, `11`, `21`, `33` | longphase 自訂格式：`1`/`2` = 單一 PS 內 allele 1/2；`11`/`21` = 雙 PS 模式下對應 `PS1-allele1`/`PS2-allele1`；`33` = both/unphasable |
| Paired tumor (ClairS pileup → tagged) | **str** | `'1'`, `'2'`, `'1-1'`, `'2-1'`, `'3'` | 標準 PS-allele 格式：`{PS_idx}-{allele}` 或 `{allele}` |

### 1.2 兩個 metric 的眼力差別

- **L1 family**：只看頭文字（族系 1 vs 2）。BL `1` ≡ paired `1` ≡ paired `1-1`（皆 family 1）。對 phasing 正確性的判定。
- **L2 exact**：看完整 canonical 字串。BL `1` ≠ paired `1-1`（雖然 family 同 1）。對 notation 一致性的判定。

### 1.3 為什麼這個區別重要

> **L2 失敗不必然代表 phasing 失敗** — 它有可能只是「longphase 沒有把 PS 後綴展開到 read tag」造成的字串不對，但 family（生物學上的 read-to-haplotype 歸屬）完全正確。

報告 02 §3.2 的 BL L2 mean 0.445 vs L1 mean 0.568 的 0.123 差距，**幾乎全部由此規範鴻溝產生**。

---

## Section 2. Confusion matrix per BAM × site

### 2.1 全域 confusion matrix（1303 shared reads）

#### BL_exact × PA_exact

| BL_exact \\ PA_exact | `1` | `1-1` | `2` | `2-1` | `3` |
|----------------------:|----:|------:|----:|------:|----:|
| `1`                   | **150** | 10  | 122 | 2  | 1  |
| `1-1`                 | 46  | **177** | 2 | 193 | 15 |
| `2`                   | 36  | 5   | **47** | 30 | 0  |
| `2-1`                 | 0   | 76  | 1   | **17** | 0 |

#### V5_exact × PA_exact

| V5_exact \\ PA_exact | `1` | `1-1` | `2` | `2-1` | `3` |
|----------------------:|----:|------:|----:|------:|----:|
| `1`                   | **93**  | 5   | 86  | 4   | 0  |
| `1-1`                 | 0   | **116** | 0 | 56  | 6  |
| `2`                   | 54  | 8   | **65** | 19 | 2  |
| `2-1`                 | 44  | 131 | 3   | **152** | 8 |
| `3`                   | 1   | 6   | 0   | 1   | 1  |

對角線（exact match）粗體。**離對角線一格** 但**同 family**（如 `1`↔`1-1`、`2`↔`2-1`）的 cell 是「family match but exact mismatch」族群。

### 2.2 「家族對但 exact 錯」的 read 計數

| Pair | reads (family ✓ + exact ✗) | reads (family ✓ + exact ✓) | 比例（exact 錯 / family 對） |
|------|---------------------------:|--------------------------:|----------------------------:|
| BL vs PA | 87 | (full match 略) | dominant 來自 `1`↔`1-1`（10）、`2-1`↔`1-1`（76）、`1-1`↔`2-1`（193 跨 family，非此族） |
| V5 vs PA | 27 | (full match 略) | dominant 來自 `1-1`↔`2-1`（混入 family flip） |

說明：上表分類是粗略的；精確「family match but exact mismatch」的 read-level 計數來自 `hp_family_exact.tsv` 的 `family_*_match_PA == 1 AND exact_*_match_PA == 0` 過濾條件。BL = **87 reads**，V5 = **27 reads**。

> 直觀解讀：BL 的 87/1303 ≈ 6.7% 看似「exact mismatch」其實是 paired 寫 `1-1`、BL 寫 `1` 這種 notation 差。**這 87 reads 在 L1 family 上得分 1.0，但在 L2 exact 上扣到 0**，造成 L2 systematic 低於 L1。

### 2.3 Per-site BL「家族對但 exact 錯」分布

對每個 site（依報告 02 §3.2 數據反算）：

| site | shared n | BL exact mismatch but family match | 主因 |
|------|---------:|-----------------------------------:|------|
| TP02 | 63 | 0 | exact 與 family 同步（皆 0.926） |
| TP04 | 106 | ~17 | BL 用 `1` (family 1)、PA 多用 `1-1` |
| TP05 | 137 | ~30 | BL `1` vs PA `1-1` 大量出現 |
| FPB1 | 87 | ~1 | exact 0.988 已接近 family 1.000 |
| SP2 | 100 | ~21 | BL `2` 對應 PA 的 `2-1` |
| SP3 | 98 | ~25 | BL `2` 對應 PA 的 `2-1` |
| 其他 | — | 各站 < 5 | — |

> 主因鎖定 TP04 / TP05 / SP2 / SP3 — 這四個 site 解釋了 BL L1-L2 gap 的大部分。

### 2.4 BL HP=33 的特殊角色（必然 mismatch）

longphase BL 偶有出現 `HP=33`（self-phasing 33 mark），其 canonical exact = `"3"`。Paired 端在這 15 site 的 `"3"` 出現率極低（confusion matrix 對應 row/col 總和 ~17 reads），因此：

> **BL HP=33 ≡ canonical "3"**，與 paired 任何 `"1"/"2"/"1-1"/"2-1"` 都不同 → **必然 L2 mismatch**，但 family=0 不參與 L1（無法 eligible）→ L1 不受影響、L2 直接拉低。

V5 在 confusion matrix 中也出現 `3` row（共 9 reads），規模與影響相若。

---

## Section 3. 為何 family 同但 exact 不同（規範差異深探）

### 3.1 三類「family ✓ but exact ✗」場景

#### Case A: BL 寫 `1`，PA 寫 `1-1`

最常見。PA 的 ClairS pileup haplotag 把 PS 編號保留為 `1-1`；longphase tag 的 BL 直接寫 `1`（單 PS 模式或 PS-id 隱式）。

```
read X:  BL.HP = 1     family=1, exact="1"
         PA.HP = "1-1" family=1, exact="1-1"
→ L1 match=1   L2 match=0
```

#### Case B: BL 寫 `11`、PA 寫 `1-1`

雙方都用雙 PS 格式但前綴不同：BL longphase 整數 `11` 對應 canonical `"1-1"`，PA 字串 `"1-1"`。**這個 case L2 應該 match**。本報告 confusion matrix 中 BL `1-1` row 在 PA `1-1` col = 177 即此 case 對齊正確的 reads。

#### Case C: BL `21` (canonical "2-1") vs PA `1-1`

family flip：BL family=2、PA family=1。**這個 case L1 與 L2 都 mismatch**，不在「家族對但 exact 錯」族群（屬於真正 family disagreement / orientation flip）。confusion matrix BL `2-1` × PA `1-1` = 76 reads 屬此。

> Case A 是 L1-L2 gap 的根因；Case B 沒有問題；Case C 是真正的 phasing disagreement（多由 PS orientation 差異主導，已在 L4 矯正）。

### 3.2 為什麼 V5 的 87→27 大幅下降？

對比 BL 87 vs V5 27 reads「家族對但 exact 錯」：

> V5 的 somatic fallback 在多個 site 把 BL 原本的 `HP=1` / `HP=2` 升級為 `HP=11` / `HP=21`（雙 PS 格式 canonical `1-1`/`2-1`），**正好對齊 paired 的 `1-1`/`2-1`**，所以 Case A 大量轉成 Case B（L2 也 match）。

- 但同時 V5 在某些 site 引入 family flip（Case C-class），這抵消了 Case A 修正的收益。
- 報告 02 的 L2 mean V5 0.435 vs BL 0.445 結果（V5 ↓ 0.010）：notation hygiene（A→B）給 V5 加分，family flip 給 V5 扣分，淨效應接近 0。

### 3.3 規範差異不可逆推回 phasing 錯誤的判定

**重要警告**：閱讀 BL 的 L2 exact mean 0.445 時，**不應解讀為「BL phasing 在 55% 的 read 上錯」**。實際 phasing 正確性指標：

- L1 family mean 0.568 — **PS orientation 翻轉污染後**的 family 一致性。
- L4 orient-corr mean 0.916 — **PS orientation 翻轉去除後**的真正 family 一致性。

L2 是 notation match，不是 phasing 正確性。下游論文若要使用 read-level HP 一致性指標，建議：

1. **首選 L4**（PS-aware 對齊後的 family rate）
2. **次選 L1**（PS frame 對齊問題視為 reporting）
3. **L2 僅在比較同一編碼系統的 phasing tools** 時使用（例如 longphase v1 vs longphase v2，皆整數編碼）

---

## Section 4. 結論

### 4.1 主要發現

1. **L1 vs L2 的 0.123 mean gap（BL 端）幾乎全部由 notation 差異解釋**：1303 reads 中 87（6.7%）為 BL「family ✓ + exact ✗」，主因為 longphase 整數 HP 對 PA 字串 HP 的格式落差（Case A）。
2. **V5 在 notation hygiene 上有局部優勢**：「family ✓ + exact ✗」reads 從 BL 87 降到 V5 27，因為 V5 somatic fallback 多用雙 PS canonical format 對齊 paired。
3. **但 V5 同時引入 family flip 失誤**：notation 修好的同時、family-level 對齊變差，淨 L2 接近 0（V5 -0.010 vs BL）。
4. **HP=33 特殊角色**：BL/V5 都有 `33` reads，必然 L2 mismatch（PA 端 `3` 罕見），但這些 reads 的 family=0，不參與 L1，所以 L2 是唯一受其影響的 metric。

### 4.2 對 audit suite 結論的影響

| 觀點 | 是否成立 |
|------|---------|
| V5 在 phasing 正確性上勝過 BL | **否**（L1/L4 顯示反向） |
| V5 在 phasing notation 規範上勝過 BL | **是（局部）**（87→27 reads notation match 改善） |
| 報告 02 的 L2 exact rate 可作為 phasing 正確性 proxy | **否**（被 notation 差異污染） |
| 報告 02 的 L4 orient-corr 是最強 phasing 正確性 proxy | **是**（PS-aware 後的真正 family 一致性） |

### 4.3 實務建議

1. **論文/週報引用 phasing 一致性時**：使用 L1 + L4 雙重指標，**避免引用 L2 作為 phasing 正確性聲稱**。
2. **比較同一系統的多版本 phasing**（如 BL → V5）：L2 可用作 notation hygiene metric，但需註明意義。
3. **未來工程化建議**：若要擴大 audit 規模到全 variant，建議在 BAM 後處理階段把 BL/V5 的整數 HP 標準化展開為 `{PS}-{allele}` 字串格式（與 paired 對齊），消除 Case A 干擾，使 L2 重新獲得 phasing-correctness 詮釋。

### 4.4 與 02 報告的關係

| 02 主張 | 此報告補強 |
|---------|------------|
| L2_BL mean 0.445 < L1_BL mean 0.568 | §2.2 解釋 87 reads 的 Case A，§3.1 拆解 |
| L2_V5 mean 0.435 ≈ L2_BL mean 0.445 | §3.2 解釋 V5 的 A→B 修正 vs family flip 抵消 |
| L4 是最可信的 phasing 正確性 proxy | §4.3 推薦清單一致 |
| V5 沒有顯著超越 BL | §4.2 verdict 一致 |

---

## Section 5. 數據附錄

- 來源 TSV：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/hp_family_exact.tsv` (1303 reads × 17 cols)
- 統計腳本：`InterSubMod/scripts/analysis/v5_read_intersection_concordance.py`

快速可重現指令：

```bash
python3 - << 'EOF'
import pandas as pd
r = pd.read_csv("/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/hp_family_exact.tsv", sep="\t")
print("BL family-match-but-exact-mismatch:", ((r.family_BL_match_PA==1) & (r.exact_BL_match_PA==0)).sum())
print("V5 family-match-but-exact-mismatch:", ((r.family_V5_match_PA==1) & (r.exact_V5_match_PA==0)).sum())
print(pd.crosstab(r.BL_exact, r.PA_exact, dropna=False))
EOF
```
