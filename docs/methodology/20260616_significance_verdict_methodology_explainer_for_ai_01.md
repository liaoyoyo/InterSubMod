---
title: "ISM 顯著判定方法學 — 完整解釋敘述（供下游解釋 AI 引用）"
date: 2026-06-16
type: methodology explainer / handoff narrative（給「陪搭解釋 AI」整理詳細解釋用）
audience: 下游解釋 AI（產出對人類的詳細解釋） + 研究者複查
scope: HCC1395 paired_full, 全基因組 30490 TP SNV, single-sample single-pipeline (⭐3)
data_sources: >
  output/_wg_d1_unified/significance_summary.csv (D1, commit 990175d),
  output/_wg_d2_dispersion/significance_summary.csv (D2, commit 74eb2e2),
  src/core/SignificanceAnalyzer.cpp, src/core/StructureTest.cpp, src/core/RegionProcessor.cpp,
  scripts/regression/golden_chr1_w5000_bernoulli{,_skip}.tsv
related: 20260616_significance_verdict_too_strict_and_input_compliance_design_01 (decision/spec),
  20260616_strength_score_reclassify_cpp_change_spec_01 (Stage①②④ spec),
  20260616_structure_label_verdict_false_negative_audit_01 (原審查),
  memory project_ism_verdict_false_negative_audit_2026_06_16
commits: {Stage②+segfault: b6dde2b, D2 PERMDISP: 74eb2e2, D1 unified verdict: 990175d}
claim_levels: 架構/程式=L1(親讀 code) / 全數字=L1(第一手讀回 landed CSV) / caveat=L2(推論)
status: 完成（D2+D1 已 merge develop）；within-HP null / Stage③ pending
---

# ISM 顯著判定方法學 — 完整解釋敘述

> **這份文件的用途**：給下游「解釋 AI」一份**自含、可精準引用**的敘述，讓它能對人類（研究者 / 教授 / 論文讀者）產出詳細解釋。每個數字都標來源檔案；每個程式邏輯都標 file:line；每個結論都標證據層級。**禁止下游 AI 在此基礎上捏造或外推未列數字**。

---

## 0. 一句話總結（給解釋 AI 的 elevator pitch）

ISM 原本的「這個位點有沒有顯著甲基結構」判定**太嚴格**（全基因組只放行 39.1%），因為它**算了 PERMANOVA 卻丟掉、又把判別離散度混淆的 PERMDISP 檢定關掉了**。本次修正先修好 PERMDISP（D2），用它過濾掉「假性顯著」（87% 的候選其實是離散度差異而非位置差異），再把剩下**乾淨的位置型顯著**併入一個**統一結構判定**（D1），最終把位點乾淨分成三層：**可能有結構（98.5%）vs 真正無訊號（1.5%）**，且不誤傷任何原本的強訊號位點（FP=0）。

---

## 1. 背景：ISM 在做什麼 + 為何需要「顯著判定」

ISM（InterSubMod）對每個 somatic SNV 周邊 ±5000bp 視窗，取覆蓋該區的長 read，建 **read×read 甲基化距離矩陣**（BERNOULLI/NHD），做階層分群，再用一組統計檢定判斷「這個區域的甲基化是否呈現有意義的結構，且與標籤（HP 單倍型 / allele REF-ALT）相關」。最終每個區域得到一個 `VerificationClass`（結構分類）與一個 `Significant`（布林顯著）。

**為什麼這個判定重要**：它決定哪些位點「值得進一步看」（有結構）vs「沒訊號可丟」。判定太嚴 → 漏掉真結構（假陰性）；判定太鬆 → 把雜訊當訊號（假陽性）。本次處理的就是「太嚴」這一側。

---

## 2. 原本的判定架構（兩個常被混淆的東西）

ISM 其實同時輸出**兩個不同的判定**——解釋時務必分清：

| 判定 | 程式位置 | 原本公式 | 性質 |
|------|---------|---------|------|
| **`Significant`（布林欄）** | `RegionProcessor.cpp:1385`（原）| `passed_gating && global_p≤0.05 && CramersV≥0.1 && num_reads≥20` | 離散 cluster×label Fisher 的**嚴格 gate**（帶 F1-filter 歷史包袱）|
| **`VerificationClass`** | `SignificanceAnalyzer.cpp:326-340`（legacy 2×2）+ `RegionProcessor.cpp` Stage④ override | cluster_significant × label_sig 的 2×2 | 結構分類（Strong/Weak/Noise…）|

🔑 **關鍵問題 FN1（原審查坐實）**：ISM **有算** label-PERMANOVA（`SignificanceAnalyzer.cpp:301-304`，含 permutation p 與 dispersion），但這個結果**完全不進** `Significant` 也不進 verdict。等於做了正規的 permutation 顯著性檢定，卻把結果丟掉、改用較嚴的離散 Fisher gate。

---

## 3. 發現的三個問題（用戶觀察 → 程式坐實）

1. **顯著判定太嚴格**：全基因組 30490 個 TP 區域，舊 `Significant` 只放行 **39.1%（11912）**；在「被判 Noise/Weak」的 8180 個池子裡，舊 gate 只放行 **65 個**。
2. **被丟掉的 PERMANOVA 其實多半顯著**：那 8180 個池子裡，**7476（91.4%）有 valid 且 p<0.05 的 label-PERMANOVA**（HP 軸 6484 / allele 軸 6658，allele 略多 → 貼 somatic/cis）。這證明「太嚴」不是沒訊號，是訊號被丟。
3. **PERMDISP（離散度檢定）根本沒生效**：原本 dispersion warning 欄**全 0**。親讀程式找到兩個根因：
   - `RegionProcessor.cpp:2348` 明寫 `enable_dispersion=false`（pipeline 直接關掉，無註解理由）。
   - 即使開啟，`StructureTest.cpp:292` 的 dispersion p 是「F 統計量 3 段查表」(`F>4→0.01 / >2.5→0.05 / else 0.1`)，不是真 p。

> **為何 PERMDISP 重要（解釋 AI 必懂的核心概念）**：PERMANOVA 檢定「兩群的甲基化中心位置有沒有差」(location)，但它在「兩群的離散度（變異程度）不同」時**也會顯著**，即使中心位置一樣（Anderson 2006）。所以一個 PERMANOVA 顯著，可能是「真的位置不同（有結構）」，也可能只是「一群比另一群更雜（離散度差異）」。**PERMDISP 就是用來分辨這兩者的**。沒有它，你無法判斷 PERMANOVA 顯著到底是不是真結構。

---

## 4. 修法（分階段，每階段都驗證過才下一步）

> 鐵律：每個 C++ 改動都跑 ctest + 雙守護 regression（golden 比對），數字全部從 landed 檔案讀回才寫入。

### 4.0 過程中先修一個 bug — MAX_DIST clustering 越界 segfault（commit b6dde2b）
做 Stage② within-HP（見 4.1）時，regression 的 MAX_DIST 路徑 segfault。用 gdb 查到：崩在 jemalloc `free()`，faulting pointer = `0x4038000000000000`（= IEEE-754 double `24.0`）→ **堆積體被一個 double 越界寫壞**。根因：`TreeCutter::find_optimal_clusters` 回傳 `best_k`，但它的群標籤**不保證是連續的 [0,best_k)**（當合併高度大量並列時——MAX_DIST 的哨兵距離尤甚——會切出比 k 多的群）。within-HP 的 silhouette 迴圈用 `std::vector<double>(k)` 拿原始標籤直接索引 → 越界寫。修法：把標籤 remap 成連續 `[0,K)`、全程用實際群數 K（主分群路徑用 `std::map` 故免疫）。

### 4.1 Stage② within-HP 結構軸（pattern + level）（commit b6dde2b）
偵測「單一 germline-HP 內部的甲基亞結構」。關鍵方法學洞見：**距離分群只抓 pattern（每個 CpG 的圖樣差異），會漏掉「水平位移」（一群整體甲基較高/低，圖樣相似）**。chr1 解構：distance 只抓 1.2%，但 level（mean β 雙峰）有 26.1%，其中 25.6% 是 distance 完全漏掉的。所以加了 **level 軸**（per-read mean β 的最佳 2-split）。

### 4.2 D2 — 修 PERMDISP（commit 74eb2e2）
- 開啟 `enable_dispersion=true`；
- 加 `betai/betacf/f_dist_sf`（Numerical Recipes 正則化不完全 beta → F 分佈上尾 p），把 3 段查表換成**真 analytic F 分佈 p**（df1=群數-1, df2=觀測數-群數）。
- **驗證**：與 `scipy.stats.f.sf` 比對 8 組 (F,df1,df2)，最大絕對誤差 **1.1e-13**；regression 雙守護 **byte-identical**（dispersion 純加欄、不進 verdict，故 9 個 golden 欄完全不變）。
- **量化結果**：那 7476 個「顯著 PERMANOVA」中，**6507（87%）是 dispersion-warned（離散度混淆，dispersion p 常 1e-10 級，極端，非大-N 邊際）**，只有 **969（13%）是 clean location-shift（真位置型結構）**。

> **D2 的決定性意義**：若略過 D2 直接認列 7476，會把 6507 個離散度假象標成「真顯著」= 從「過嚴」直接換成「過鬆」。**「D2 先於 D1」的排序是必要的。**

### 4.3 D1 — 統一結構判定（commit 990175d）
**用戶定向**：F1-filter 已結論為 DEAD → 只要**一個統一輸出判定**；把「可能有結構」與「真無訊號」乾淨分開；保留所有 per-axis 原始數值欄 + StrengthScore。

- 在 Stage④ VC override 加 **`clean_location_permanova`** 軸（HP/allele PERMANOVA valid & p<0.05 & **非** dispersion-warned）→ 新細類 **`PermanovaLocation`**。
- **`dispersion_structure`**（同樣顯著但 dispersion-warned）→ 新獨立 tier **`DispersionStructure`**（承認「結構存在但離散度型、歸因模糊」，不混入 location 顯著）。
- `Significant`（col8）由舊 F1 gate **重定義 = VC ∈ valid-location 細類**（統一判定的布林投影）。
- per-axis 原始欄（所有 PERMANOVA/dispersion/Δβ/CramersV/HP-AUC）+ StrengthScore **全保留不動**。
- **驗證**：ctest 221/221；regression 前 7 欄 byte-identical、只 col8+col9 變（符設計）→ --update-golden + 雙守護 PASS。

---

## 5. 全基因組最終結果（30490 TP SNV，SKIP，第一手讀回 `output/_wg_d1_unified/significance_summary.csv`）

**三層分線（用戶要的「可能有結構 vs 真無訊號」）**：
| 層 | 類別 | 數量 | 佔比 | `Significant` |
|----|------|-----:|-----:|:---:|
| ① clean location structure | Strong 22310 / MultiGroupNoLabel 1899 / LabelShift 1809 / **PermanovaLocation 1756** / LOH-Structure 376 / StructureNoLabel 177 | **28327** | **92.9%** | true |
| ② dispersion-type | **DispersionStructure** | **1700** | 5.6% | false |
| ③ truly unhandleable | Noise_Uniform 311 / Noise_Uncorrelated 122 / Noise_Chaotic 30 | **463** | **1.5%** | false |

- **可能有結構（①+②）= 30027（98.5%） vs 真無訊號（③）= 463（1.5%）**。
- **`Significant` 舊 39.1%（11912）→ 新 92.9%（28327）**（delta +16415）。
- **FP=0**：legacy Strong 21015 個全部保留為 valid，**降級 0**（cluster_match→Strong 恆 valid，構造保證）。

---

## 6. 誠實 caveat（解釋 AI 必須一併傳達，不可只報好消息）

🟡 **92.9% valid-location 偏高**。兩個原因：
- (a) 合理生物：這些全是 TP（真 somatic）區域，本來就多半有甲基結構。
- (b) 方法限制：PERMANOVA 在 read 數多（大-N）時，對**微小**的 location 差異也會報顯著（過度敏感）。

→ **D1 達成了「分出那 1.5% 真無訊號」這個目標，但它是「結構特徵化的三層分層」，不是 TP 內部的細判別器**。而且 TP/FP 判別（用甲基當 F1 filter）本專案早已結論為 **DEAD**，所以這不是 D1 的目標。
→ 若未來要更利的分線：可加 **PERMANOVA effect-size 門檻**（不只 p<0.05，要求效應量夠大）/ 做 **TP-vs-FP noise-rate 對照**（看 FP 區域的 Noise 比例是否顯著更高）/ within-HP 的 permutation null（見 §8）。

---

## 7. 輸入合規護欄（重要原則，解釋 AI 須強調）

**分類 / 顯著判定的輸入只能用「資料本身可導出」的量**（BAM 衍生：甲基矩陣、HP tag、read coverage）。**外部真值（SEQC2 CN/LOH BED）只能在判定後做驗證/觀察，絕不可當分類輸入**——否則 (a) 循環/洩漏（用答案分類）、(b) 不可泛化（SEQC2 只有 HCC1395 有，方法就跑不了其他樣本）。對應：
- **LOH** ← **HP-tag 比例失衡**（`potential_loh = hp_ratio<0.1||>0.9`，`SignificanceAnalyzer.cpp:322`，已存在，合規）。
- **CN** ← **read coverage / depth peak 自估**（self-contained，Stage③ 待做）。
- **SEQC2** ← 降為**獨立 post-hoc 驗證腳本**。原 Stage③「載 SEQC2 BED 當輸入」**作廢重設計**。

---

## 8. 開放項（pending，解釋時標「尚未完成」）

1. **within-HP permutation null（記錄用）**：within-HP clean multigroup 目前是 heuristic 幾何判準（silhouette/gap/balance），**沒有統計 null**。要讓它能稱「顯著」需補 permutation null，但有 **select-then-test 循環性**陷阱（先選最佳分群再測會膨脹 p）：pattern 軸需 gap-statistic 式參考分佈、level 軸需 unimodal null（per-read mean β 對排序置換不變，單純 label 置換無效）。**設計尚未定案，未實作。**
2. **Stage③（coverage-peak 自估 CN）**：CN-Substructure / LOH-aware 細類，用自估 CN（非 SEQC2）。未實作。
3. **92.9% 過敏細化**（effort-size 門檻 / TP-FP 對照）。未實作。

---

## 9. 給解釋 AI 的「怎麼對人講」要點

- **核心類比**：PERMANOVA 像「測兩群人的平均身高有沒有差」；PERMDISP 像「測兩群人的身高離散度有沒有差」。一群人若平均一樣但一群高矮很不齊，PERMANOVA 也可能說「有差」——這時要靠 PERMDISP 拆穿。原本 ISM **沒開 PERMDISP**，所以分不出「真有結構」vs「只是比較雜」。
- **這次修正的精神**：不是「放寬標準把更多東西算顯著」，而是「先把被誤丟的真顯著（969 個乾淨的）找回來、同時把假性顯著（6507 個離散度型）誠實標成另一類」。**先嚴謹分辨，才談認列。**
- **必守的誠實**：修完 92.9% 算「有結構」是偏高的，要講清楚這是「分出真無訊號」而非「精準判別每個位點」。
- **必守的合規**：SEQC2 等外部真值只能驗證、不能當輸入——這是方法可信度與可泛化的底線。
- **不可外推**：單樣本（HCC1395）單 pipeline，⭐3。跨樣本復現未做。
