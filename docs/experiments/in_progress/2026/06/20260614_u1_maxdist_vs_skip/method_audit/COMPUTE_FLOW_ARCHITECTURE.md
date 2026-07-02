---
title: "ISM 顯著性計算流程架構 — A/B 兩方向整合 + Δβ 檢定 + 效率設計"
date: 2026-06-14
type: methodology-architecture / 計算流程設計
status: DRAFT — 待用戶確認架構後逐步實作
claim_levels: 現有順序=L1(源碼) / 架構建議+效率=L2(提議,待驗證)
related: TWO_DIRECTIONS_FRAMEWORK.md, SIGNIFICANCE_INVENTORY.md
---

# ISM 顯著性計算流程架構

## 0. 先回答兩個問題

### Q1: somatic residual Δβ 是否應該內建到 ISM？→ **應該，且「值」已內建，只缺「檢定」**
- **值已在**：`HP_Signed_Residual`（RegionProcessor:1031）= `tumor_hp_signed_delta − normal_hp_signed_delta`
  = (tumor 的 mean HP1−mean HP2) − (normal 的 …) = **tumor-specific somatic HP 甲基差異**。tumor/normal signed
  delta 已各自算（:1019-1020），residual 只是相減 → **零額外計算成本**。
- **缺檢定**：目前無 p/sig（純描述）。應內建檢定（見 §4 #3 修法）。**檢定成本低**：Δβ = O(reads×CpG) 的 mean，
  permutation shuffle 比 distance permutation（O(reads²)）便宜一個量級 → 內建可行。
- **結論**：✅ 內建。把 `HP_Signed_Residual` 從「描述欄」升級為「somatic ASM 判定」（取代有缺陷的 `hp_residual_sig`）。

### Q2: 方向 A 對 subclone 的價值（修正你指出的）
A（label-first）**不是對 subclone 無用，而是「另一個角度的確認」**，且有 B 沒有的優勢：
| A 的優勢 | 說明 |
|---|---|
| **對齊外部方法** | DSS/methylKit/DAMEfinder 全是 label-first Δβ 檢定 → ISM 可直接對標、reviewer 熟悉、可比較 |
| **直接量化 + 可解釋** | Δβ 有方向(hypo/hyper)+大小，比 distance/CramersV 易解釋 |
| **明確標籤** | HP/sample 是已知輸入，不像 cluster 要先分群(無 overfit) |
| **somatic 確認** | sample 標籤(T/N) + somatic residual Δβ → 直接確認 somatic ASM(不需發現 subclone) |
| **per-CpG 定位** | Fisher per-CpG 知道「哪些位點」 |
→ **A 與 B 互補**：B 發現結構(subclone 候選) → A 確認「該結構對應標籤(tumor-specific)的甲基差異是否顯著」。

## 1. 計算依賴 + 共用中間量（效率的根基）
```
reads → [階段0] methylation matrix (raw β + binary) + labels(HP/sample/allele)   ← 算一次, 全下游共用
              ↓
        [階段1] distance matrix (read×read)                                        ← 算一次, B + A-distance 共用
              ↓
   ┌──────────┴──────────┐
[階段2 方向B structure]   [階段3 方向A label]
 用 distance              用 methylation raw + labels
```
**兩個一次性大中間量**：`distance matrix`（B 全用 + A 的 distance delta）、`methylation raw`（A 的 Δβ + PerCpgAsm）。
labels 算一次共用。**現有 ISM 已共用這兩個**（distance :843, raw :844）→ 效率基礎已對。

## 2. 計算流程架構（建議的清楚分階段）

| 階段 | 方向 | 方法 | 用什麼中間量 | 必要性 / 面向 |
|---|---|---|---|---|
| **0 基礎** | — | methylation matrix(raw+binary) + labels | reads | ✅ 必要 |
| **1 距離** | — | distance matrix | raw(binary/raw) | ✅ 必要 |
| **2 B-structure** | B | **HP-AUC**(距離直接 vs 標籤) | distance + labels | ✅ 結構-標籤關聯**最乾淨** |
| | B | clustering(UPGMA+cut) → cluster labels | distance | ✅ 結構發現(subclone 入口) |
| | B | ClusterPermanova(cluster 結構顯著) | distance + cluster | ✅ 結構顯著性 |
| | B | GlobalTest(cluster × label CramersV) | cluster + labels | ✅ 但與 HP-AUC 部分交集 |
| **3 A-label** | A | **Δβ 模組**: germline Δβ(normal HP1-HP2)+檢定 / **somatic residual Δβ(tumor-normal)+檢定** | raw + labels | ✅ **#3+#1 核心**, 對齊外部 |
| | A | PerCpgAsm Fisher(+NME/Epipoly) | raw + binary + HP | ✅ per-CpG 定位 + 異質性(不同面向) |
| **4 整合** | A+B | VerificationClass(綜合判定) | 全部 | ✅ |

## 3. 各方法優缺點 + 交集（你要的「是否有交集 / 可否合併」）

### 交集（測同樣事 = 可精簡）
| 交集組 | 各自做什麼 | 處置建議 |
|---|---|---|
| **HP-AUC ∩ PERMANOVA-label_hp ∩ CramersV(cluster×HP)** | 都測「甲基結構 vs HP 標籤」(pairwise / grouped / clustered 三角度) | **保留 HP-AUC(最乾淨)** 為主；PERMANOVA-label/CramersV 為補充角度(不冗餘,但要明示分工) |
| **LabelTest distance delta(HP/allele) ∩ PERMANOVA-label_hp/allele** | **完全冗餘**：都按標籤測距離結構(between-within delta vs pseudo-F) | 🟡 **擇一**(PERMANOVA pseudo-F 更標準) → 精簡掉 LabelTest distance delta |
| **hp_residual(distance) ∩ somatic residual Δβ** | 都測 somatic HP ASM(距離 vs 甲基率差) | 🔴 **棄 hp_residual(distance, 借 tumor p)**, 用 **somatic residual Δβ + 檢定**(可解釋+對齊外部+可檢定) |

### 不同面向（都要，不可省 — 互補非交集）
- **B HP-AUC**（結構是否對應標籤）vs **A Δβ**（標籤的甲基方向/大小）：一個問「結構像標籤嗎」、一個問「標籤甲基差多少」。
- **A Δβ region-level**（整體 HP 甲基差）vs **A PerCpgAsm per-CpG**（哪些位點）：巨觀 vs 定位。
- **A PerCpgAsm Fisher**（位點差異）vs **NME/Epipoly**（epiallele 熵）：差異 vs 異質性。

## 3.5 Δβ 配對矩陣設計（統一回答你問的所有 HP 配對：哪些必要、對目標有用）

**標籤層次**：family = HP1f(HP1+HP1-1) vs HP2f(HP2+HP2-1)；fine = HP1/HP1-1/HP2/HP2-1（"1-1"/"2-1" =
longphase 同 haplotype 內**帶 somatic 的 sub-haplotype**，ReadParser:128-138）；sample = tumor/normal。

🔴 **現狀**：ISM 有 fine pairwise **distance**（6 對 `HPFineD_*`，`Stats.hpp:324-326` 已註解語意），
**但無 Δβ pairwise + 無檢定**。你問的配對正是補這個 Δβ 線。

| Δβ 配對 | 層 | 意義（Stats.hpp 已註解）| 對目標 | 必要 | 檢定 |
|---|---|---|---|:---:|:---:|
| **HP1f − HP2f** | family | germline parental（納 somatic-carrier，**power 高**）| **germline ASM 主指標** | ✅ | ✅ |
| HP1 − HP2 | fine | germline parental（純，排 somatic）| germline ASM（純版，read 少）| ✅ | ✅ |
| **HP1 − HP1-1** | fine same-hap | **same haplotype, germline vs somatic** | 🔑 **HP1 內 subclone** | ✅✅ | ✅ |
| **HP2 − HP2-1** | fine same-hap | same haplotype, germline vs somatic | 🔑 **HP2 內 subclone** | ✅✅ | ✅ |
| HP1-1 − HP2-1 | fine | 兩 somatic-carrier 之間 | somatic 跨-hap 一致性 | 🟡 | 🟡 |
| HP1−HP2-1 / HP2−HP1-1 | fine cross | germline vs 異-hap somatic | 診斷 / confound 對照 | 🟢 | ✗ |
| **tumor residual − normal residual** | sample×HP | somatic residual Δβ（=`HP_Signed_Residual`，值已有）| somatic ASM（扣 germline baseline）| ✅ | ✅ |

**核心結論（回答「是否必要、對目標有用」）**：
- **family Δβ（HP1f−HP2f）= germline ASM 主指標**（對齊 GlobalTest family 主 gate，power 最高）→ 必要、對 germline ASM 目標有用。
- **fine same-hap Δβ（HP1−HP1-1 / HP2−HP2-1）= 方向 A 偵測 subclone 的核心**（haplotype 內 germline vs
  somatic 甲基差異 = P3b 的 label-first 版）→ **你問的配對裡最有目標價值的**（直接對 subclone）。
- somatic residual Δβ（tumor−normal）→ 必要（#3 修法，somatic ASM）。
- 跨配對 → 只當**診斷/confound 對照**（值存、不檢定）。

**效率（零成本配對）**：
1. **一次掃 raw matrix**，算 4 fine 組 ×{tumor,normal} 的 mean methyl（≤8 個 mean）。
2. **任意配對 Δβ = 兩 mean 相減 → 零額外成本**（所有配對都免費，包括跨配對診斷）。
3. **檢定只對 ✅ 必要配對做 permutation**（germline ASM family/fine + subclone same-hap + somatic residual），
   跨配對只存值。→ 計算集中在「算 8 個 mean + 對 ~5 個目標配對做 mean-permutation」，遠比現有
   6 對 distance pairwise（O(reads²)）便宜。

## 4. 效率改進 + 程式碼架構順序（建議）

### 效率改進點
1. ✅ distance/raw 已一次共用（基礎對）。
2. **somatic residual Δβ 值零成本**（已算 tumor/normal signed → 相減）。
3. **Δβ 檢定內建便宜**（mean-based permutation ≪ distance permutation）→ 加在 §3 Δβ 模組。
4. **精簡冗餘**：移除 LabelTest distance delta（HP/allele）+ hp_residual(distance) → 省 2 套 permutation(distance, 貴)。
5. **labels 算一次**（HP family / sample / allele）傳遞共用，避免各方法重複解析 hp_tag 字串。

### 建議程式碼順序（process_single_region 內，vs 現有）
| 現有順序(L1) | 建議 |
|---|---|
| :844 methylation → :869 **clustering+significance(B+A distance)** → :902 baseline → :927 HP delta(A distance) → :997 Δβ(A) → :1039 PerCpgAsm(A) | **階段化重排**：0 基礎 → 1 distance → 2 B(HP-AUC+clustering+PERMANOVA+GlobalTest) → 3 A(Δβ模組[germline+somatic residual+檢定] + PerCpgAsm) → 4 整合 |
| Δβ 散在 :997-1031, 無檢定 | **Δβ 模組集中**: `compute_dbeta_module()` 一次算 germline/somatic Δβ + permutation 檢定 + 輸出 dbeta/p/sig |
| hp_residual(distance) + HP_Signed_Residual(Δβ) **並存** | 棄 hp_residual(distance), 統一到 somatic residual Δβ |

## 5. 待用戶確認 → 逐步實作順序
1. **確認架構**（§2 分階段 + §3 交集處置 + §4 順序）對齊你的構想嗎？
2. 逐步實作（每步驗證）：
   - (a) **Δβ 檢定模組**（#3 修法核心）：somatic residual Δβ + germline Δβ + permutation 檢定 → 內建輸出欄。先 permutation(i)，beta-binomial(ii) 後續對照。
   - (b) 實際數據驗：somatic residual Δβ sig 比例 vs 舊 hp_residual_sig，與外部(DSS-like)對齊。
   - (c) 精簡冗餘（LabelTest distance delta / hp_residual distance）— 確認無下游依賴後移除。
   - (d) 階段化重排（效率）— 確認結果不變(雙守護 regression)。
> ⚠ 每步「對齊想法/目標 + 確認數據/實際狀況」。架構是工具，§2-§4 若有偏差先校正。
