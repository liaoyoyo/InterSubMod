---
title: "CN/LOH-aware 結構分類體系 — 思路/場景/數據量化/分析/判別設計"
date: 2026-06-16
type: design (理解思路 → 場景 taxonomy + 數據量化 → 分析計算設計 → 判別設計；尚未實作)
trigger: 用戶提出更豐富的結構分類思路(Δβ保留 / 結構對標籤 / CN within-label多結構 / 整區CN / nested / LOH)
data_sources: /tmp/dbeta_wg_abc/significance_summary.csv (全基因組 30490), Coverage_Category/Potential_LOH/LOH_Subtype 欄
claim_levels: 場景=L1(用戶定義+源碼欄) / 量化=L1(全基因組 CSV) / 分析判別設計=L2(待決策)
status: AWAITING 理解確認 + 設計決策
---

# CN/LOH-aware 結構分類體系設計

## 1. 核心思路（一句話）

**germline-HP 標籤是「germline 相位」標籤；但 somatic 事件（CN 擴增 / LOH / subclone）改變了實際細胞群結構** → 真實甲基結構可能**比 2 個 germline-HP 多**（CN→多 copy / subclone→germline+carrier）或**少**（LOH→丟一 HP）。所以「結構不對上 germline-HP 標籤」常是**真生物訊號（CN/LOH/subclone），不是 noise**，應保留。

## 2. 場景 taxonomy + 實際數據量化

| # | 場景（你的描述）| 結構-標籤關係 | 成因 | 全基因組數據 | 應判 |
|---|------|------|------|------|------|
| **S1** | Δβ 有差異就保留 | mean β 差，**可無 distance 結構** | 真甲基均值差 | **~80% of 假陰性**（HP-AUC<0.7 但 Δβ sig）= **主因** | 保留(signal) |
| **S2** | 結構對上標籤 | cluster ↔ label | germline ASM/somatic | 原 Strong/Subclone | Strong |
| **S3** | Δβ + 結構對標籤都有 | both | clear | — | Strong(最高) |
| **S4** | 一標籤內 ≥2 結構 | structure > label, within-HP | **CN 擴增**(germline HP tag 但 somatic 複製) | CN-altered 佔假陰性 50%(1.2x)；明確場景4/5 ~5% | 保留(CN-substructure) |
| **S5** | 多標籤在一大結構各有小結構 | nested, labels in 1 big struct | **整區 CN** | CN 子集 | 保留(CN-region) |
| **S6** | 大結構內多小結構但標籤各自聚集 | nested, labels cluster within | 標籤結構相近 | 需 within-HP 分群測 | 保留 if labels cluster |
| **S7** | 一群甲基但標籤無效 | mono-structure, degenerate label | **LOH**(丟一 HP) | LOH 區 2605(9%), LOH-假陰性 183 | 特殊(LOH-mono) |
| **S8** | LOH 內甲基小結構 | sub-structure in LOH | LOH-subclone | LOH_Subclone 144 | 保留(LOH-subclone) |

**🔑 誠實主次**：S1（Δβ mean-shift，~80%）是假陰性主因；S4-S8（CN/LOH）真實但各 5-8% 少數。**修法優先序：S1(Δβ 保留) > S2(結構對標籤/HP-AUC) > S4/5(CN within-label) > S7/8(LOH)**。

## 3. 如何分析/計算（existing + 待補）

| 訊號 | 現有 | 待補 |
|------|------|------|
| S1 Δβ mean 差 | ✅ Δβ 模組(germline/subclone/alt) | — |
| S2 結構對標籤 | ✅ HP-AUC + GlobalTest | — |
| S4/S6 **within-HP 多結構** | 🔴 **無** | **新增：對單一 HP-family 的 read 再 clustering → NGroups_within_HP（測一個 HP 內是否 ≥2 甲基群）** |
| S4/S5 CN | ✅ Coverage_Category(read-count CN proxy) | **SEQC2 CN cross-check(外部真值)** + per-HP allelic coverage(allele-specific CN) |
| S7/S8 LOH | ✅ Potential_LOH/LOH_Subtype/LOH_Bed | LOH 區的 within-結構(同 within-HP 分群) |
| nested(S5/S6) | 🔴 部分(silhouette 偏小 k) | hierarchical k 多層 / 放寬 max_k |

**最關鍵待補 = within-HP 再分群**（測「一個 germline-HP 內有幾個甲基群」），直接抓 S4/S6（CN/subclone 造成的 within-label 結構），現有 HPFine 只分 germline-vs-carrier 不夠。

## 4. 如何判別（分層 decision，保留所有真結構）

> 原則：**任一真結構證據成立就保留**（Δβ / HP-AUC / within-HP 多群 / LOH 結構），並依成因標細類。

```
valid_signal = Δβ_sig (S1)
            OR cluster_matches_label (S2/S3)
            OR HP_AUC>=0.7 (S2 連續)
            OR within_HP_NGroups>=2 + CN-consistent (S4/S5/S6)
            OR LOH_internal_structure (S8)

判別細類:
  Δβ + 結構對標籤        → Strong
  結構對標籤(無Δβ)        → Strong/Subclone
  Δβ mean-shift(無結構)   → LabelShift (S1, 新類; 現被埋 Weak/Noise)
  within-HP 多群 + CN     → CN-Substructure (S4/5)
  LOH + 內結構           → LOH-Structure (S7/8)
  結構但正交標籤(非CN非LOH) → StructureNoLabel
  皆無                   → Noise
```

**CN/LOH cross-check 的角色**：不是判別門檻，而是**解釋成因 + 防偽**——within-HP 多群「是否 CN 造成」用 Coverage_Category/SEQC2 CN 確認；LOH-mono 用 Potential_LOH 確認。

## 5. 數據量化（2026-06-16）

**within-HP 再分群 pilot**（per-read mean β bimodality, 全基因組）：
- any-HP 內 ≥2 群（loose, 兩群各≥3）= **12416（40.7%）**；其中 CN-altered 41.4%（=基線）→ **非 CN 主導，多為 epiallele/subclone 層級**
- 🔴 40.7% 為**上界**（含「3 vs 57」小離群，同 subclone tiny-group 病）；refined（silhouette≥0.5+平衡+gap）收斂中 → 真乾淨多群 < 40%

**SEQC2 CN cross-check**（`/big8_disk/data/HCC1395/SEQC2/CNV/`, Masood 2024）：
- ISM region **93% 落 SEQC2-CN 區**（gain 56%/loh 37%/loss 1%）→ hyper-diploid，**CN/LOH 是常態**
- 假陰性落 CN 區 94% vs 全體 93%（**1.0x 未富集**）→ **CN-membership 無判別力（太普遍）= context 非門檻**
- ISM Coverage_Category vs SEQC2-gain: precision 79%/recall 47%；ISM Potential_LOH vs SEQC2-loh: **recall 僅 13%**（大量漏 LOH）

**🔑 數據定案**：① CN/LOH 普遍(93%)→germline-HP 幾乎總是潛在不足→**強力支持「預設保留、認結構超越標籤」哲學** ② 但 CN/LOH **不能當分類門檻**（太普遍），是成因解釋 ③ 判別靠**實際甲基結構**（乾淨多群/Δβ/HP-AUC）④ LOH-aware 用 **SEQC2 LOH 當真值**（ISM recall 13% 不可靠）

## 6. 精煉判別邏輯（用戶哲學：預設保留、只有真沒訊號才 Noise）

```
# 有結構證據（任一成立 → valid）— "innocent until proven noise"
valid = Δβ_sig                              # S1 level-shift (主因~80%FN)
     OR HP_AUC>=0.7                          # S2 distance 對 HP
     OR cluster_matches_label(GlobalTest)    # S2/S3
     OR clean_multigroup(sil>=0.5 & balanced)# 多群乾淨(不論對不對標籤; 用戶新加)
     OR within_HP_clean_multigroup           # S4/S6 (一HP內乾淨多群)
     OR LOH_structure(SEQC2-loh & 內結構)    # S7/S8

# 只有「真沒訊號」才 Noise（用戶定義，須全否上述 + 三子型）
Noise = NOT valid AND 標籤交錯無一致, 細分:
   Noise_Uniform   : 甲基均勻(low variance, NME 低) — 附近甲基一致
   Noise_Chaotic   : 高熵無乾淨群(NME 高 & silhouette 低) — 甲基變動混亂
   Noise_Uncorrelated: 有變動但不成群也不對標籤

# 細類(valid 內依成因)
Strong / LabelShift(Δβ無結構) / MultiGroupNoLabel(乾淨多群不對標籤) /
CN-Substructure(within-HP多群+SEQC2-gain) / LOH-Structure / StructureNoLabel
```

**關鍵區分（用戶最在意）**：乾淨多群(valid) vs 混亂高熵(noise) = **silhouette/ClusterPERMANOVA（有無乾淨群）× NME/Epipoly（熵）** 二維。

## 7. 決策點（待你確認）
- D1：§6 判別邏輯（預設保留 + clean-multigroup-vs-chaotic 二維 + CN 當 context 非門檻）符合你的哲學嗎？
- D2：LOH-aware 改用 **SEQC2 LOH 真值**（ISM recall 13%）— 同意嗎？
- D3：細類粒度（6 valid 類 + 3 noise 子型）OK 或先合併？
- D4：refined within-HP（乾淨多群真實規模）+ noise 三子型量化後 → pilot 新判別邏輯（如前 reclassify_pilot）→ 確認淨改善才 cpp-change。
