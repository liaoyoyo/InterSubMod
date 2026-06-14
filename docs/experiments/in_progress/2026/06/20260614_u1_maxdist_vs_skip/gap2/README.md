---
title: "gap2 — normal 對照判別真 somatic subclone 比例 + HP 顯著性判定方法審查"
date: 2026-06-14
type: experiment / methodology-audit
status: partial — germline 主導確認, 精確 somatic 比例待修方法
scope_flag: "⚠ 單樣本 HCC1395 paired; SKIP wg 11912 sig (全基因組)"
data_sources: gap2_subclone_discrimination.json
source_refs: "RegionProcessor.cpp:927-993 / LabelTest.cpp:149-166"
---

# gap2：真 somatic subclone 比例 + HP 顯著性判定方法審查

## L0 結論（謹慎）
U1 SKIP 的 11912 個顯著 region 中，**真 somatic subclone 上限 ~15%（1830）**，**主導為 germline-HP**
（強證據：normal_HP_valid 96.9% + tumor signed delta 0.032 < normal 0.043 + tumor/normal signed corr 0.01）。
**🔴 但 ISM 現有的 somatic HP 判別（`hp_residual_sig`）方法有缺陷**（顯著性借用 tumor 自己的 p、未對殘差檢定），
**精確 somatic 比例無法從現有欄位可信回答 → 需修方法**。gap2 確認方向（germline 主導），但完整判別屬第三部分方法審查待解。

## 數量觀察（← gap2_subclone_discrimination.json，SKIP wg 11912 sig）
| 指標 | 值 | 判讀 |
|------|---:|------|
| Normal_HP_Valid | 11539 (96.9%) | 絕大多數 sig 在 normal 也有 HP 分群 = germline |
| tumor_valid & normal_valid | 6934 (58%) | 兩邊都可測 HP |
| **normal_only**（tumor 無 HP / normal 有）| 4605 (39%) | tumor single-HP（somatic haplotag），只 normal 分群 |
| tumor_only | 117 (1%) | 極少 tumor-specific HP |
| **tumor vs normal signed-delta corr** | **0.01** | HP 甲基方向 tumor/normal 幾乎不相關 |
| \|tumor signed\| vs \|normal signed\| median | **0.032 < 0.043** | tumor HP 甲基差比 normal **還小** |
| somatic 候選（tumorValid & \|signed residual\|>0.1）| **1830 (15.4%)** | 真 somatic **上限**（均值維度粗估）|
| germline-dominant（\|signed residual\|≤0.05）| 3374 | 明確 germline |
| signed residual NaN（single-HP）| 3610 (30%) | tumor/normal single-HP 無法算殘差 |

## 🔴 方法質疑清單（源碼 L1 確認 — 對應你「判定方法正確性有質疑」）
1. **`hp_residual_sig = tumor_hp.significant`（RegionProcessor.cpp:978）** — 殘差 delta 算對了
   （tumor−normal :976），但**顯著性直接借 tumor 自己的 p-value，未對殘差做統計檢定**。
   germline 案例（tumor+normal 都顯著、delta 相近、殘差≈0）會被誤標 somatic sig → **44.9% HP_Residual_Sig 不可信**。
2. **HP test = distance permutation（LabelTest.cpp:149-166）非甲基均值差** — `delta = 群間距離 − 群內距離`，
   `delta>0 才 permutation`、`p≤alpha 顯著`。→ **CramersV（distance）高、signed delta（均值）小可並存**（兩個維度）；
   我 gap2 用 signed delta 判別**未對齊** ISM 的 distance-based HP test（故 15% 僅均值維度粗估）。
3. **殘差 = distance delta 相減**（tumor_hp.delta − normal_hp.delta）— distance delta 是 permutation-based
   非線性量，**相減的統計意義不明確**（非線性可加）。
4. **somatic 位點 tumor single-HP**（:971-973 註解，39% normal_only）— tumor HP test 常 invalid →
   `hp_residual_sig=false`（:992）；這類 region 的 somatic 訊號**無法用此路徑判**。
5. **99 permutation → p 下限 0.01** — 顯著性解析度粗（U1 大量 p=0.01 即此）。

## 結論與待解
- **germline-HP 主導 = 強證據確認**（96.9% normal 也分群 + tumor signed < normal + corr 0.01）→ 與「甲基=germline-haplotype 層級」一致。
- **真 somatic 精確比例待修方法**：需對「殘差（tumor−normal distance delta）」做**正規 permutation 檢定**（而非借 tumor p），才能可信判 somatic。屬第三部分方法審查產出。
- gap2 回答了「方向」（germline 主導、somatic ≤15%），但揭露 ISM somatic HP 判別需方法修正才能給精確數量。

## 產物
- `gap2_subclone_discrimination.json` — 數量分析
- 後續：第三部分（甲基分群 + 標籤顯著性判定方法）系統審查
