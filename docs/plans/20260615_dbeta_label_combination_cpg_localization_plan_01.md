---
title: "Δβ 標籤組合 + CpG 定位 + 盲點交叉驗證 — 查核發現與方案決策"
date: 2026-06-15
type: plan (查核 → 方案 → 決策；尚未實作)
trigger: 用戶問「Δβ 是否驗證所有標籤組合差異 / 高效計算 / 輸出哪些組合差異 / 定位差異 CpG / Δβ 盲點是否有他法驗證」
data_sources: /tmp/dbeta_chr1_v2 per-region reads.tsv (2624 region 量化), 6-method audit /tmp/w4ukyquky.output, dbeta_stage3_crosscheck.json
claim_levels: 查核發現=L1(源碼+本run量化) / 方案=L2(設計提案待決策)
status: AWAITING DECISION
---

# Δβ 標籤組合 + CpG 定位 + 盲點 — 查核與方案

## 1. 查核發現（現況 L1）

### 1.1 標籤其實有「兩個獨立軸」，現在 Δβ 只用一個
| 軸 | 來源 | 值 |
|----|------|----|
| **HP-tag** | somatic-haplotag **block 相位**（longphase）| 0/1/2/1-1/2-1/3 |
| **alt_support** | read 自己這顆 SNV 的**鹼基** vs VCF REF/ALT | ALT/REF/UNKNOWN |

**量化（chr1 2624 region, tumor reads）— 兩軸不對齊**：
- carrier-tag(1-1/2-1) 中 **13%（10696）alt=REF**（phased 到 carrier block 但此 SNV 不帶 somatic）
- germline-tag(1/2) **0% alt=ALT**（germline 群乾淨 = 純 REF）✅
- unphased(HP=0) 中 **237 reads alt=ALT**（帶 somatic 但無 HP tag → 被 Δβ 排除）

→ 現在 subclone Δβ「germline-tag vs carrier-tag」= 乾淨germline-REF vs（carrier-block: 87% ALT + 13% REF）。**carrier 群非純 ALT，且 237 個 unphased-ALT 完全沒進分析**。

### 1.2 現在 Δβ 算哪些組合（只 4 個 HP-tag 配對）
| Δβ | 群定義 | 軸 |
|----|--------|----|
| germline ASM | normal HP1-family vs HP2-family | HP-tag |
| subclone HP1 | tumor "1" vs "1-1" | HP-tag |
| subclone HP2 | tumor "2" vs "2-1" | HP-tag |
| somatic residual | (tumor HP1−HP2)−(normal HP1−HP2) | HP-tag |

**未做**：alt_support 軸的 subclone / cross-HP carrier 對（1-1 vs 2-1）/ 全組合枚舉 / 「哪些組合有差異」per-region 輸出。

### 1.3 CpG 定位現況（PerCpgAsm）
- `PerCpgAsm.cpp:244` 把 carrier **併入 family**（1-1→HP1-family）→ per-CpG 只定位 **germline HP1-vs-HP2 軸**。
- **無** subclone（germline vs carrier）或 alt 軸的 per-CpG 定位。

### 1.4 B 方向（structure-first）已覆蓋什麼
- `HPFine`（4-group HP1/HP1-1/HP2/HP2-1 PERMANOVA, `hp_fine_f/p`）+ `GlobalTest fine` → **有考慮 4-group 結構**，但 distance-multivariate、仍 HP-tag 軸、非 pairwise Δβ、不指方向、不定位 CpG。

---

## 2. 逐題回答（你的 5 問）

| # | 問題 | 答案（L1 事實）|
|---|------|------|
| 1 | Δβ 是否驗證所有標籤組合差異？合理嗎？ | **否**，只 4 個 HP-tag targeted 配對。合理（每個對應一個生物目標）但**漏 alt_support 軸 + unphased-ALT + cross 組合**。 |
| 2 | 能否高效算所有組合？ | **能**。per-read mean β **算一次**，任何重新分群（HP/alt/組合）都 O(reads)。組合枚舉幾乎零額外成本。 |
| 3 | 能否輸出「哪些組合有差異」？ | **能**（新增 per-region combo 表 / 欄）。需 FDR（多重檢定）。 |
| 4 | 能否定位差異的 CpG？ | **目前不能**（PerCpgAsm 只 germline 軸、併 carrier）。需擴成「任意群對 per-CpG Δβ + FDR」。 |
| 5 | 會影響目標分析嗎？ | **不會**（純加欄/加表，不改現有 subclone/germline Δβ）。targeted subclone 仍主，組合/定位是 enrichment。 |

---

## 3. Δβ 盲點 + 交叉驗證現況（你的第二組問）

**Δβ 三大盲點**：
1. **mean-only**：只看平均，**漏 bimodal/mixture**（一群內兩亞群平均相同）→ **B structure-first（clustering/PERMANOVA）catches**。
2. **單軸（HP-tag）**：忽略 alt_support 13% 不一致 + 237 unphased-ALT → **目前無方法交叉驗證此軸選擇**（盲點）。
3. **CpG-collapsed**：read mean β 抹平 per-CpG 模式 → **無 subclone per-CpG 定位**（盲點）。

**交叉驗證真的驗了嗎**（誠實）：
- ✅ stage 3（beta-reg/MWU）驗了**檢定本身**（perm 非 artifact，88-94% 一致）→ 但**用同一 HP-tag 分群**，**沒驗軸選擇盲點**。
- ✅ B 方向（HPFine PERMANOVA / GlobalTest）驗了**結構覆蓋**，且與 A 真的不一致 ~10%（gap#6: 185 region PERMANOVA-sig 但 Δβ 一尾 non-sig）→ **A/B 互補已實證**。
- 🔴 **未驗**：alt_support 軸的 subclone（盲點2）、subclone 的 per-CpG 定位（盲點3）= **本提案要補的交叉驗證**。

---

## 4. 方案（3 個可加元件，皆 additive 純加欄/表）

> 共用既算的 per-read mean β → 全部 O(reads) 高效；不改現有目標分析。

- **方案 A｜alt_support 軸 subclone Δβ**（targeted 交叉驗證盲點2）
  - tumor HP1 內 alt=ALT vs alt=REF；HP2 同理。= 更直接的 somatic 對比（read 自己的 allele，獨立於相位）。
  - 與現有 HP-tag subclone **並列輸出** → 兩定義一致則 robust，不一致則揭相位 vs allele 的差異。
  - 4 新欄（HP1/HP2 × Δβ/sig）。最小、最對齊你的問題。

- **方案 B｜全組合 Δβ 表**（exploratory）
  - 枚舉此位點實際出現的 (HP × alt × T/N) 群，兩兩 pairwise Δβ + perm + min-group guard + BH-FDR。
  - 輸出 per-region **combo 表**（哪些組合對有顯著差異 + Δβ + 群大小）→ 直接答「哪些組合有差異」。
  - 較重（per-region 副表檔 or 寬欄）；多重檢定需 FDR。

- **方案 C｜subclone per-CpG 定位**（補盲點3）
  - 擴 PerCpgAsm（或新「任意群對 per-CpG Δβ」）：對顯著的 subclone/組合，逐 CpG 算群差 + BH-FDR → 輸出**驅動差異的 CpG 清單**。
  - 答「定位差異 CpG」。中等成本（per-CpG 已有 Fisher 框架可複用）。

---

## 5. 決策點（待你選）

- **D1**：方案 A（alt 軸 subclone）做不做？（推薦做 — 最小、直接交叉驗證盲點2、對齊你的核心問題）
- **D2**：方案 B（全組合表）做不做？（exploratory，較重；或先 pilot 看有無新訊號再決定）
- **D3**：方案 C（subclone per-CpG 定位）做不做？
- **D4**：先 **pilot 驗證**（Python 小範圍跑，看 alt 軸/組合/CpG 有無新訊號）再決定 cpp-change，還是直接 cpp-change？（推薦 pilot-first — 避免加欄後發現無增益）

> 所有方案 **additive 不影響目標分析**（§2 Q5）。建議順序：**先 pilot（Python，重用 per-read β）驗證方案 A 的 alt 軸 subclone 與 HP-tag subclone 是否一致/互補 → 有價值才 cpp-change 內建**。
