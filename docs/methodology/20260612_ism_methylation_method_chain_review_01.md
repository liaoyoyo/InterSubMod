---
title: ISM 甲基方法鏈 review — 取樣 / 位點保留判定 / 甲基差異與關聯驗證
date: 2026-06-12
type: methodology-review
status: in-discussion（與用戶逐點確認中）
branch: chore/ism-review-param-governance-202606
data_sources:
  - output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2/intersubmod_tp/significance_summary.csv
  - output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2/intersubmod_fp/significance_summary.csv
provenance:
  - 深讀 workflow wf_b87378c0-cbf（3 agents, 2026-06-12, read-only, file:line 對 develop HEAD 069cadb）
  - 數量已由主 agent 獨立 python 複核（2026-06-12，cascade 加總精確等於總數）
related:
  - InterSubMod/docs/method_comparison/20260609_ism_vs_external_methylation_tools/01_ism_method_spec_from_source.md
  - InterSubMod/docs/methodology/20260610_code_methodology_detail_audit_01.md
---

# ISM 甲基方法鏈 review（取樣 → 位點保留判定 → 甲基差異與關聯驗證）

> 標記：✓ 合理 / ⚠ 可討論 / 🔴 疑點。數字標 [複核] = 主 agent 親自 grep CSV 確認。

## L0 整體結論

> ISM 甲基方法鏈**整體 sound、底層數值教科書級正確、無翻轉結論方向的錯誤**。系統性議題集中在 **3 個會影響「解釋哪些位點 / 結論口徑」的設計選擇**（保留 gate 的 latent 張力、ASM 軸別 confound、跨基因組多重檢定缺失）+ 數個工程一致性疑點（C_min 雙預設、5hmC 丟棄、CramersV 雙口徑）。**沒有「方法錯到結論不能用」的情況，但有「方法可以更好 / 口徑要對齊」的明確改進點。**

---

## ① 甲基數據取樣的方式與範圍

**鏈路**：BAM `MM/ML` → MethylationParser → ReadParser/ReadAggregator 過濾 → MatrixBuilder（read×CpG 稀疏矩陣）→ binarize → DistanceMatrix。**取樣窗 = somatic SNV ± `window_size_bp`（預設 1000bp）**。

| 關鍵步驟 | file:line | 行為 |
|---|---|---|
| 甲基解析 | `MethylationParser.cpp:20-87` | **只解析 5mC（`C+m?`）；5hmC（`C+h?`）只 advance ml_offset 後丟棄** |
| ML 機率 | `:140,180` | prob = ML/255 |
| CIGAR seq→ref | `:232-281` | M/=/X 雙進；I/S seq-only；D/N ref-only |
| 雙鏈摺疊 | `:102-196` | forward C 正掃 / reverse G 反掃，都報 forward C 1-based，**摺疊到同一 CpG column** |
| read 過濾 | `ReadParser.cpp:26-69` | drop secondary/supp/dup/unmap；MAPQ≥20；len≥1000；須 MM+ML |
| tumor 錨定 | `ReadAggregator.cpp:24-71` | tumor 須 SNV 為 ALT/REF（UNKNOWN 丟）；normal 跳過 alt-filter；qname dedup |
| 建矩陣 | `MatrixBuilder.cpp:17-104` | 所有見過的 CpG 成 column，**無 per-site coverage filter**，missing=−1.0 |
| 二值化 | `RegionProcessor.cpp:1408-1426` | raw≥0.8→1，≤0.2→0，**中間 (0.2,0.8) 視為 ambiguous→missing** |
| 距離 C_min | `DistanceMatrix.cpp:31-55` | 共同 CpG ≥ C_min（runtime=3）否則 MAX_DIST=1.0 懲罰 |

**評價**：
- ✓ **雙鏈摺疊 forward CpG**（對稱 5mC 生物學正確；reverse 反向迭代正確）— T1 audit 曾標零測試，commit 6593f96 已補 golden test
- ⚠ **5mC-only，5hmC 丟棄**（partial）— 對 5mC-driven ASM OK，但**隱含未標註**，且與 Python tsg 層的 5mC/5hmC max-collapse 口徑不一致（T4 ASM-3）
- 🔴 **C_min 雙預設衝突**（concern, T5 CONFIRMED）— `Config.hpp:37`=3（runtime 生效）vs `DistanceMatrix.hpp:23`=5（被覆蓋）；**`01_spec` 文檔寫 5 但實跑用 3 → 文檔與碼不一致**。另 `min_site_coverage=5` 定義但**從未被引用 → 等於關閉**
- ⚠ 中間帶 (0.2,0.8) 當 missing，丟失 cell 比例未知（NEEDS_RUN，需在真實 region 上量）

**未來會調的取樣參數**：`window_size_bp`(1000) / `min_mapq`(20，MAPQ-60 blind spot) / `min_read_length`(1000) / `binary 0.8/0.2` / `C_min` / `use_full_read_span`(false)。

---

## ② 位點保留+解釋的判定 + 數量　③ 非保留原因

**唯一保留 gate**（`RegionProcessor.cpp:1143-1146`）：
> **`Significant = PassedGating AND GlobalP≤0.05 AND CramersV≥0.1 AND NumReads≥20`**（四條 AND）

判定前的多相分析：前置 gate（reads≥2/CpG≥1 才算距離；reads≥`clustering_min_reads`(10) 才分群+顯著性，否則根本不寫進 summary）→ 分群（階層+TreeCutter 找 k）→ Phase1 全域關聯（Fisher-FFH + CramersV，多層 HP 軸）→ Cochran 可靠性閘 → 雙閘 gating → Phase3 PERMANOVA → Phase5 VerificationClass。

### 數量 [複核 2026-06-12]

| | TP 集 | FP 集 |
|---|---|---|
| 計算位點（significance_computed）| **29,754** | 627 |
| **保留（Significant=true）** | **719（2.42%）** | **3（0.48%）** |
| 非保留 cascade：未過 gating | 26,925（90.5%）| 603 |
| 　　GlobalP>0.05 | 347 | 9 |
| 　　CramersV<0.1 | 1,763 | 12 |
| 　　reads<20 | 0 | 0 |
| VerificationClass | Weak 19435 / Noise 7840 / Strong 1940 / Subclone 539 | Weak 405 / Noise 207 / Strong 13 / Subclone 2 |

> cascade 加總精確 = 總數（719+26925+347+1763+0=29754 ✓）。**FP 保留率(0.48%) < TP(2.42%)，方向正確**。

### 🔴 核心張力：latent 位點（最大可改進點）

- **2,108** 個非保留位點其 PERMANOVA `valid` 且 `p≤0.05`（**有結構**），其中 **2,105（99.9%）是因 CramersV 被 Cochran gate 成 0** 才卡在保留 gate。[複核]
- 機制（`RegionProcessor.cpp:1592-1596`）：`result.cramers_v = max(各層 reliable? v : 0)` → 稀疏列聯表（小群 / HP-fine 多組）必然 Cochran-unreliable → V=0 → 卡 `CramersV≥0.1`。**等於用「卡方近似可靠性」否決了「permutation 已證實的結構」**。

### 評價

- ⚠ **保留 gate 四條 AND**（partial）：嚴格合理（顯著性×effect size×depth）、FP<TP 方向對，**但靈敏度極低（TP 僅 2.42% 保留）**
- 🔴 **CramersV gated 成 0 餵進保留 gate**（concern）：Cochran 規則本身 sound，問題在「把 unreliable V 設 0 再當保留條件」→ 系統性丟棄稀疏表真結構
- 🔴 **raw vs gated CramersV 雙口徑**（concern）：`apply_gating`(GlobalTest.cpp:141) 用 **raw** CramersV，最終保留用 **gated** CramersV → 邏輯雙重否定不透明
- ⚠ **reads≥20 hardcode**（`:1146`）與可計算門檻 10 不一致 → 10-19 reads 位點「算了但永不保留」灰區（此 TP 集覆蓋深，gate 不 binding；淺覆蓋樣本會 binding）
- ✓ **VerificationClass 四分類**（合理）：比單一 Significant 更能解釋非保留性質（Subclone=有結構無 label / Weak=有 label 無強結構）

**audit 背書**：`20260324 方法學審查` 已判定 Noise 含 33-67% 真 TP、Significant 捕獲率僅 1-12%，**改用 `VC!=Noise` 可把捕獲率提到 33-68% 且 precision 相近**。

---

## ④ 甲基差異與關聯的驗證/判定方法（重點）

**5 個正交統計引擎**：

| 引擎 | file:line | 量什麼 | 評價 |
|---|---|---|---|
| **per-CpG Fisher exact + BH-FDR** | `PerCpgAsm.cpp:274-338` | 率差存在性（HP1 vs HP2 meth 比例 2x2）| ⚠ partial |
| **LabelTest Δ=d_between−d_within + 999perm** | `LabelTest.cpp:147-217` | 距離空間群分離（**非率差**）| ✓ yes |
| **PERMANOVA pseudo-F + dispersion** | `StructureTest.cpp:93-201` | 結構顯著性（ISM 最 distinctive）| ✓ yes |
| **CramersV + FFH + Cochran gate** | `GlobalTest.cpp:73-143` | 關聯強度 | ⚠ partial |
| **normal-anchored cis-test** | `RegionProcessor.cpp:884-983` | somatic = tumor Δ − normal Δ（扣 germline 基線）| ✓ yes |

**逐引擎評價**：
1. **Fisher**（⚠）：底層教科書級正確（固定邊際超幾何 + log-sum-exp）。議題 = binary call 當獨立 Bernoulli 忽略同 read/clone over-dispersion。**但**：`08 audit` 用真資料 read-level permutation 反證 **FP 端不膨脹**（null FP 0.7-1.9% < 5%，合成模擬的 17-20 倍膨脹真資料不重現）；`11` 修正 **敏感度端漏 8.5-17%** 結構顯著但率測弱位點。→ **FP 端不需改，敏感度端值得加 effect-size 門檻或 beta-binomial-with-delta**
2. **LabelTest Δ**（✓）：distance-based permutation 恰當；Δ≤0 直接 p=1 避免測反方向。**警示**：Δ（距離差，無單位 NHD）≠ Δβ（率差），下游必須區分
3. **PERMANOVA**（✓）：ISM **最可防守新穎**（03 確認無外部工具同時有 read-read 距離矩陣+PERMANOVA）；數值工程漂亮（adj_term 防 catastrophic cancellation、完美分離回 1e9 sentinel）。**瑕疵**：dispersion p 用硬編碼階梯（f>4→0.01）非真 F CDF（僅 advisory）；**production 把 `enable_dispersion` 關了**（`:1575`）→ 失去「PERMANOVA 顯著是否只是 group 變異差」的內建 FP 偵測
4. **CramersV**（⚠）：方法 sound，但雙口徑 + 稀疏表 under-report（見 §② latent）；`chi_square_p` 是 placeholder（= fisher p，FISHER-2，下游無消費）
5. **normal-anchored cis-test**（✓）：**直擊 ALLELE-axis baseline confound 核心**；設計洞察正確（residual 距離平移不變→改 tumor/normal 子集相減；tumor 單 HP dominance 是 EXPECTED biology 非 bug）。**改進在口徑治理**（T4，見下）

---

## ⑤ 最值得逐點討論的爭議（排序）

| # | 爭議 | 類型 | 影響 |
|---|---|---|---|
| **A** | 保留 gate 改「**PERMANOVA 顯著 OR CramersV-reliable**」？現用純 CramersV AND 否決 2105 latent | 🔴 最大可改進 | 「能解釋多少位點」靈敏度 |
| **B** | **ASM 軸別 confound（T4）**：C++ 預設把 somatic 1-1 併入 germline family-1 → 預設 Δβ / per-CpG Fisher 量的是 **germline/ALLELE 軸（被 baseline confound）非 somatic-HP**；輸出欄名未明示軸別 | 🔴 最關乎結論正確性 | somatic ASM claim 是否成立 |
| **C** | Fisher over-dispersion 的 split 結論（FP 端維持/敏感度端改）是否接受 | ⚠ 方法選擇 | ASM 存在率口徑 |
| **D** | **跨 region/genome 多重檢定校正缺失**（per-CpG 層內有 BH，跨基因組無）；08 標這才是 ASM 存在率偏高真來源 | ⚠ 方法選擇 | 所有 ASM 存在率描述 |
| **E** | C_min 3 vs 5 統一 + min_site_coverage 未用 | ⚠ 工程一致性 | reproducibility |
| **F** | 5hmC 丟棄 vs Python max-collapse 口徑 | ⚠ 工程一致性 | 跨層口徑 |
| **G** | raw vs gated CramersV 雙口徑統一 | ⚠ 工程一致性 | 報告可解釋性 |
| **H** | reads≥20 hardcode 與可計算門檻 10 的灰區 | ⚠ 工程一致性 | 淺覆蓋樣本保留數 |
| **I** | production 關 dispersion（失去內建 FP 偵測）是否開回 | ⚠ 方法選擇 | PERMANOVA FP 把關 |

> 這些**多數是 C++ 改動（Hard Gate：commit 必編譯）**，須走 methodology-audit → cpp-change；本輪先逐點理解+確認方向，不動碼。
