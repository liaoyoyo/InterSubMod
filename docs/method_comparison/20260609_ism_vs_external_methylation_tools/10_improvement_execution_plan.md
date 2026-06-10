<!--
建立時間: 2026-06-11
狀態: PLAN (執行計畫 — 可改進處的清楚紀錄, 供後續執行; 尚未動 code)
報告類型: improvement_execution_plan
受眾: 執行者(自己/agent) · PI 核准 · cpp-change/methodology-audit 入口
framework: Verdict-Pyramid + 執行 checklist (每項: 差別 / 為何有用 / 怎麼改 file:line / 狀態+依據)
data_sources:
  - 05_improve_learn_recommendations.md (外部工具比較)
  - benchmark/08_external_program_audit.md (外部程式實機稽核, 修訂優先序)
  - docs/methodology/20260610_code_methodology_detail_audit_01.md (全碼方法學稽核 T1-T6)
  - 01_ism_method_spec_from_source.md (ISM 源碼 file:line)
provenance_note: 每項對應(a)一個外部工具/方法 +(b)一個 ISM 源碼 gap(file:line)。優先序已被 08 外部程式稽核「實機修訂」——標 DOWNGRADED 者為真資料證明非必要。本檔不含新跑數字, 引用值皆既有 validated/audit 來源。動 C++ 須走 /methodology-audit → /cpp-change(6步PDD)。
-->
<!-- provenance-verified: 引用數字(r=0.98/null FP 0.7-1.9%/BRCA2 等)皆 benchmark/08 + cis_scan_full validated; file:line 引 01 源碼確認; 本檔純規劃不產新數字。 -->

# 10 — ISM 可改進處執行計畫（差別 · 為何有用 · 怎麼改）

> **這份是什麼**：把所有「可用來改善 ISM 的部分」整理成**可後續執行**的清單。每項保留：**① 差別（現狀 vs 改進）② 為何可能有用 ③ 怎麼改（檔案/方向）④ 狀態（外部稽核後修訂）**。
> **三來源整合**：外部工具比較（05）+ 外部程式實機稽核（08）+ 全碼方法學稽核（T1-T6）。

---

## L0 — 一句話 + 最重要的修訂

**ISM 核心方法沒有「會翻轉結論」的錯誤**（外部程式稽核 08 + 全碼稽核一致）。可改進處分三類：**🔴 補洞/補測試（真高 ROI）· 🟡 統計校準（有用非錯誤）· 🟢 從外部工具學習採納（feature 增強）**。

> 🔴 **最重要的修訂（誠實）**：我先前列「Fisher→beta-binomial」為**最高 ROI #1**，但 **08 外部程式稽核用真資料推翻了它**——per-CpG Fisher 的 null 假陽性只有 **0.7–1.9%（<5%，沒偏樂觀）**，稽核引用的「53–68% 膨脹」是**合成模擬真資料不重現**。→ **per-CpG Fisher 不需改**；真正的統計 gap 是**跨 region/genome 的多重檢定校正**（不同東西）。**不要白做 beta-binomial 改寫。**

---

## L1 — 優先序總表（外部稽核修訂後）

| 優先 | 項目 | 類 | 來源 | 狀態 | 動哪 |
|:---:|------|----|------|------|------|
| **🔴 1** | 跨 region/genome **多重檢定校正** | 補洞 | 全碼稽核 SIGLABEL-2/REGION-3 | NEW（真 gap）| C++ + Python 匯總 |
| **🔴 2** | MM/ML **golden 單元測試** + 座標慣例測試 | 補測試 | 稽核 T1 + 08 Audit A | 邏輯已驗對(r=0.98)，須 pin | tests/（新增）|
| **🔴 3** | `fill_report.py` **null refuse** 補洞 | 補洞 | 稽核 T2 PYTOOL-1 | 反捏造基礎設施 | Python |
| **🟡 4** | PERMANOVA **99→999 perm** + 真 F 分布 | 校準 | 稽核 PERMANOVA-3/5 | 校準 | C++ |
| **🟡 5** | 刪/實作 `chi_square_p` placeholder + FisherExact p=0 修 | 校準 | 稽核 FISHER-2/4 | 誤導欄/inf | C++ |
| **🟢 6** | 距離矩陣用 **soft 機率**（非硬二值）| 學習 | cvlr / pycoMeth | feature 增強 | C++ |
| **🟢 7** | **CpG-pair tuple 共甲基 log-odds**（HP-stratified）| 學習 | DAMEfinder | 新軸(長讀優勢) | C++/Python |
| **🟢 8** | **Cohen's h** + M-value 交叉檢核 | 學習 | modkit / Du2010 | effect-size 穩健 | Python |
| **🟢 9** | **停 5mC/5hmC max-collapse** → 分軌 | 學習 | modkit / Du2010 | 真新穎軸 | Python |
| **🟢 10** | **qFDRP** 作 coverage-robust disorder 對照 | 學習 | Scherer 2020 | 操作化「結構非失序」| C++ |
| **🟢 11** | **imprinting ICR 正控** | 學習 | germline benchmark | 補正控 | 分析 |
| **🟢 12** | **modkit dmr** 當 genome-wide pre-filter | 學習 | modkit | 省 compute | pipeline |

> ❌ **不要做（外部稽核確認非必要 / 會白做）**：
> - per-CpG Fisher → beta-binomial **改寫**（08 Audit B' 證 per-CpG Fisher null FP<5%，不偏樂觀）。
> - **重寫 MM/ML 解析邏輯**（08 Audit A 證 r=0.98 正確；只補測試 pin，勿改邏輯）。

---

## L2 — 逐項執行細節（差別 / 為何有用 / 怎麼改）

### 🔴 1. 跨 region/genome 多重檢定校正 ← 真正的統計 gap
- **差別**：現狀 = 只有**單一 region 內** per-CpG BH-FDR（`PerCpgAsm`）；**跨 region / genome-wide 沒有任何 MT 校正**，且 `global_p = min(p_alt, p_hp, p_hp_family)` **取最小 p 不校正**（稽核 SIGLABEL-2）。改進 = region p 跨位點做 BH-FDR；min-p 配 Šidák/Bonferroni。
- **為何有用**：掃全基因組數萬位點時，每位點取 min-of-3-p 又不校正 → 「ASM 存在率」假陽性膨脹。**這才是 08 Audit B' 沒測到、但真實存在的偏樂觀來源**（per-CpG Fisher 本身 OK，問題在 aggregate 取 min + 跨位點未校正）。
- **怎麼改**：(a) `global_p` 組合處：min-p 改 `1−(1−min_p)^k` (Šidák, k=被 min 的檢定數) 或 Bonferroni；(b) Python region 匯總層：對「全位點 region p」加 BH-FDR 欄。檔案：`SignificanceAnalyzer.cpp` global_p 組合 + tsg Python 匯總腳本。**動 C++ 前走 /methodology-audit。**

### 🔴 2. MM/ML golden 單元測試 + 座標慣例
- **差別**：現狀 = `MethylationParser` 的 MM/ML 解析、反股**倒序**座標（碼自標 "CRITICAL" `MethylationParser.cpp:103-157`）、5mC/5hmC 分支 **零單元測試**。改進 = 加 hand-constructed golden test。
- **為何有用**：08 Audit A 證**邏輯正確**（vs modkit r=0.976/0.993），但**沒有測試鎖住** → 未來改動可能靜默改壞每個下游甲基數字。且 Audit A 匹配率只 **98/285（34%）** → 座標慣例（0/1-based、per-strand vs combine-strands）待釐清（稽核 IOREAD-6）。
- **怎麼改**：`tests/test_methylation_parser.cpp`（新增）：正股 + 反股 + 多 CpG + 邊界 read，手算答案比對；加座標慣例測試對齊 modkit pileup。**只補測試，勿改解析邏輯**（已驗對）。檔案：tests/ 新增 + 對照 `MethylationParser.cpp`。

### 🔴 3. fill_report.py null refuse 補洞
- **差別**：現狀 = `fill_report.py:38,54` 把 `None`/`NaN` 渲染成字面 `'None'`/`'nan'`，而非 refuse。改進 = None/NaN 視同 missing → refuse render。
- **為何有用**：§13-A「template+data 注入缺 key 必 refuse」**對 null 值破功** → 數據是 null 時報告靜默出現字串 `None`/`nan`（看起來像有值），繞過反捏造防線。
- **怎麼改**：`fill_report.py` 的 `resolve()` 把 `value is None or isnan` 一併放進 `missing`；`fmt()` 遇 None/NaN raise 或回 `{{待填}}`。檔案：`scripts/fill_report.py:38,54`。

### 🟡 4. PERMANOVA 99→999 perm + 真 F 分布
- **差別**：現狀 = 生產配置 **99 perm**（p 下限 0.01）；dispersion `anova_p` **硬寫死 3 段查表**（`StructureTest.cpp:292` F>4→0.01/F>2.5→0.05/else 0.1）。改進 = 999+ perm；anova_p 改真 F-CDF 或 permutation。
- **為何有用**：genome-wide 多重檢定需要 p<0.01 的解析度；硬查表不是真分布、會 0.5× 懲罰 heuristic_score。
- **怎麼改**：`StructureTestConfig.n_permutations` 99→999（或 9999）；`anova_f_test` 的 p 改 boost::math F-CDF 或 permutation。檔案：`StructureTest.cpp:292` + config。

### 🟡 5. chi_square_p placeholder + FisherExact p=0
- **差別**：現狀 = `chi_square_p = fisher_ffh.p_value; // Placeholder`（`GlobalTest.cpp:124`，欄名說卡方實為 Fisher p，下游無人用）；`FisherExact` RxC Monte-Carlo p = `n_extreme/n_samples` 強訊號時可得 **p=0**（→ −log10=inf）。改進 = 刪欄或實作真 χ² CDF；MC p 改 `(n_extreme+1)/(n_samples+1)`。
- **為何有用**：誤導欄位 + p=0 造成 inf。
- **怎麼改**：`GlobalTest.cpp:124` 刪 chi_square_p 或實作 χ² CDF；`FisherExact.cpp:285-314` 改 `(b+1)/(m+1)`（R `simulate.p.value` 慣例）。

### 🟢 6. 距離矩陣用 soft 機率（學 cvlr/pycoMeth）
- **差別**：現狀 = `DistanceMatrix.hpp:27-28` 先二值化 `>0.8 甲基/<0.2 未甲基`，丟掉 ML 信心；中間值算 ambiguous。改進 = 用 expected methylation（連續 p）算距離。
- **為何有用**：硬二值丟掉 ML 不確定度；低覆蓋（ISM 罰到 MAX_DIST=1.0 處）噪音大。cvlr Bernoulli-mixture + pycoMeth BernoulliPosterior 都保留 soft → 改善距離矩陣 + 降低覆蓋噪音。
- **怎麼改**：`DistanceMatrix` 加 `use_soft` 選項，NHD/Bernoulli 改用連續 p（如 |p_i − p_j| 期望距離）。檔案：`DistanceMatrix.cpp`。**動 C++ 走 /methodology-audit。**

### 🟢 7. CpG-pair tuple 共甲基 log-odds（學 DAMEfinder）
- **差別**：現狀 = per-CpG Fisher **丟掉 on-read cross-CpG pattern**（只在 NME/epipoly 用 unsupervised 版）。改進 = 加 HP-stratified tuple log-odds `(MM·UU)/(MU·UM)`。
- **為何有用**：DAMEfinder 核心洞見 = on-read MM/UU vs MU/UM 連動帶 allele-specificity，single-CpG 測不到。ISM 有理想基質（完整 read×CpG + **長讀可超短讀 150bp cap**，任意距離 CpG-pair）。
- **怎麼改**：`PerCpgAsm` 加 tuple 函數（HP1 vs HP1-1 的 CpG-pair concordance log-odds）。檔案：`PerCpgAsm.cpp` 新增。

### 🟢 8. Cohen's h + M-value 交叉檢核（學 modkit/Du2010）
- **差別**：現狀 = Δβ 只有 raw mean + Cohen's **d**。改進 = 加 **Cohen's h（arcsine 飽和-robust）** + headline loci 的 **M-value 交叉**。
- **為何有用**：raw Δβ 在 0/1 端 heteroscedastic（BRCA2/TBC1D16 都在飽和區）。**08 已驗證 modkit 原生輸出 Cohen's h**（BRCA2 −0.38[0.34–0.42]）→ 可直接照搬公式。Du2010：report Δβ 但 inference 用 M-value。
- **怎麼改**：`cis_asm_core.py` / `scripts/03` Δβ 計算加 Cohen's h（`2·asin√p1 − 2·asin√p2`）+ headline loci 算 Δ(M-value) 確認 sign/顯著存活。檔案：`pipeline/lib/cis_asm_core.py`, `scripts/03`。

### 🟢 9. 停 5mC/5hmC max-collapse → 分軌（學 modkit/Du2010）
- **差別**：現狀 = `collapse_modtype` 把每 read/CpG 的 5mC、5hmC 兩列**取最大值合成 any-modification**（`cis_asm_core.py:31-32`）。改進 = 保留兩 mod-type 分軌 / Dirichlet-Multinomial。
- **為何有用**：max-collapse 稀釋 5hmC（與已知 MSA Level1 dup-bug 同源）。**5mC/5hmC 分軌做 ASM 是文獻幾乎沒人做的真新穎軸**（05 KEEP #3）。
- **怎麼改**：`cis_asm_core.py:31-32` 改成保留 `m` 與 `h` 分別算 β，或 Dirichlet-Multinomial(h/m/C)。檔案：`pipeline/lib/cis_asm_core.py`。

### 🟢 10. qFDRP 作 coverage-robust disorder 對照（學 Scherer 2020）
- **差別**：現狀 = CORE4 用 coverage-fragile epipolymorphism 當失序對照。改進 = 加 **qFDRP**（與 ISM NHD **同 normalized-Hamming kernel**，最 coverage-robust）。
- **為何有用**：操作化「結構非失序」—— 做 2×2：高 qFDRP×非顯著 PERMANOVA=stochastic erosion；高 qFDRP×顯著 PERMANOVA=結構化 ASM。同時修 O11 暴露的 epipoly coverage artifact（AUC 0.845→0.530）。
- **怎麼改**：`PerCpgAsm` 對照槽加 qFDRP（read-pair 平均 normalized Hamming → per-CpG scalar）。檔案：`PerCpgAsm.cpp`。

### 🟢 11. imprinting ICR 正控
- **差別**：現狀 = ISM 只有 HP-shuffle null，**無 known-truth 正控**。改進 = 在 canonical imprinted ICR（H19/IGF2、SNRPN）跑 ISM 證明引擎偵得到已知 ASM。
- **為何有用**：cvlr/pycoMeth/NanoMethPhase 都這樣 validate；reviewer 會問「你的引擎在真有 ASM 時偵得到嗎」。便宜、高 credibility。
- **怎麼改**：分析層加 ICR 位點集當正控測試（不改 binary）。

### 🟢 12. modkit dmr 當 genome-wide pre-filter
- **差別**：現狀 = CORE 1-3 重（UPGMA O(N³) + 999 perm 每 region）。改進 = 全基因組先跑 modkit dmr（快 Rust）triage，再對通過率差/MHL screen 的 region 跑昂貴結構檢定。
- **為何有用**：降 compute 不改結論層；modkit 重疊 ISM 的 Δβ 層正好當快速 triage。
- **怎麼改**：pipeline 前置 `modkit dmr pair`（08 已驗指令）篩 region。檔案：pipeline 編排。

---

## L3 — 執行紀律與排程建議

1. **動 C++ 必走流程**：`/methodology-audit`（審查 + 量化影響 + 方案）→ 用戶選方案 → `/cpp-change`（6 步 PDD，每步 commit + 編譯 Hard Gate）。本檔的 🔴1/🟡4/🟡5/🟢6/🟢7/🟢10 都涉 C++。
2. **先做純加值、零風險的**：🔴2（補測試）、🔴3（Python 補洞）、🟢8/🟢9/🟢11（Python/分析層）—— 不動核心演算法、不會翻轉任何結論。
3. **🔴1 多重檢定校正**是唯一「會改變既有數字」的（會降低顯著位點數）→ 須先 /methodology-audit 量化「ASM 存在率」掉多少，並確認**不影響已 concluded 的方向結論**（filter-DEAD 等靠 LOSO/跨樣本，非單 CpG Fisher）。
4. **每項改完都要 regression**：對 BRCA2/TBC1D16 重跑，確認 Δβ/PERMANOVA/silhouette 不被改壞（08 的 audit JSON 可當 baseline）。
5. **誠實邊界**：所有「為何有用」是**機制論證 + 外部對照**，非已證 uplift。真正 uplift 要改完 A/B test。本檔是**方向**，不是保證。

## Provenance
- 優先序修訂依據：`benchmark/08_external_program_audit.md`（modkit Audit A r=0.98 / permutation Audit B' null FP 0.7–1.9%）。
- file:line：`01_ism_method_spec_from_source.md`（源碼確認）。
- 全碼稽核 T1-T6：`docs/methodology/20260610_code_methodology_detail_audit_01.md`。
- 外部工具比較：`05_improve_learn_recommendations.md`。
- ⚠ 範圍：改進方向基於單樣本 HCC1395 + 2 anchor loci 的稽核 + 文獻；genome-wide A/B 未做 → 是**規劃**非已證 uplift。
