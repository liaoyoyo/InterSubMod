<!--
建立時間: 2026-06-10
問題類型: C++ 方法 | 統計方法 | 特徵設計 | Python 分析層 | 測試覆蓋
影響 track: TO | paired | 兩者
狀態: pending_decision
scope: 全部程式碼（C++ 全 23 core+io+utils+tests+headers 窮盡；Python 策略性子集=論文數字腳本+反捏造工具，非全 181）— PARTIAL flag on Python
data_sources: docs/methodology/20260610_code_audit_assets/workflow_findings_parsed.json, docs/methodology/20260610_code_audit_assets/quantification.json
方法: 12-reviewer Dynamic Workflow (wf_ea727f9a-ef5, 28 agents) 靜態深讀 + CRITICAL/HIGH 對抗驗證 + 主迴圈 live 真值量化
敘述框架: Verdict-Pyramid + BLUF（結論先行）+ MECE 主題分群
-->

# InterSubMod 程式碼方法學與細節錯誤稽核（全碼窮盡）

> 用 **Verdict-Pyramid + BLUF**：結論先講，再分 6 主題降序揭露。**L0→L3 分層**——讀到夠就停。
> 行內標記：🔴 必修 / ⭐ 重要 / ◽ 細節｜CONFIRMED=對抗驗證通過 / REFUTED=假陽性已濾 / UNCERTAIN=機制真但量級存疑。

---

## L0 — 一句話結論

**程式碼沒有「會翻轉任一論文結論方向」的錯誤；但有 3 個必修的真問題（甲基化解析零測試、反捏造工具有洞、per-CpG Fisher 忽略 over-dispersion 使 ASM 顯著性偏樂觀），以及一批口徑/校正/magic-constant 的系統性議題。memory 記載的多數「已知 C++ bug」其實已修或本就合理。**

- 稽核規模：**105 findings**（16 HIGH → 對抗驗證後 **7 CONFIRMED / 3 REFUTED / 6 UNCERTAIN**）、42 MEDIUM、41 LOW、6 INFO。
- 共同主旨：**p 值系統性偏樂觀（偏小）+ ASM 軸別/口徑跨層不一致**。兩者都影響**量值與校準**，**不影響已發表 NEGATIVE 結論的方向**（filter-DEAD、單樣本 caveat、phasing Grade B+ 等都站得住）。
- ⭐ meta 結論：**現行碼比 memory/CURRENT_FOCUS 的「待修 bug 清單」乾淨**——引用已修、KDE 已 wired、CramersV gate 是合理的 Cochran 規則。舊筆記 stale。

---

## L1 — 重點邏輯（6 主題，降序；先修 T1-T3）

| # | 主題 | 代表 finding | 嚴重度 | 影響什麼 |
|---|------|------------|:---:|---------|
| **T1** 🔴 | **甲基化解析「最易錯的碼」零測試** | TESTS-1 (CONFIRMED, HIGH) | 🔴 高 | MM/ML strand-flip + 5mC 解析餵**所有**下游數字，卻無單元測試 pin 正確性 |
| **T2** 🔴 | **反捏造基礎設施有洞** | PYTOOL-1 (CONFIRMED, HIGH) | 🔴 高 | `fill_report.py` 把 null/NaN 渲染成字面 `'None'/'nan'` 而非 refuse → §13-A「構造上無法捏造」被繞過 |
| **T3** ⭐ | **統計校準：p 值系統性偏樂觀** | FISHER-1 (CONFIRMED, MEDIUM) | ⭐ 中 | per-CpG Fisher 忽略 read over-dispersion + region/genome 無 MT 校正 + placeholder p → ASM「存在性率」偏高 |
| **T4** ⭐ | **ASM 軸別與 Δβ 口徑跨層不一致** | ASM-2 (UNCERTAIN→MEDIUM) | ⭐ 中 | C++ 預設 Δβ/Fisher 走 germline/ALLELE 軸（`1-1` 併入 `1`）+ 5mC-only + cell-weighted；與論文要的 somatic-HP 軸 + any-mod 有口徑差 |
| **T5** ◽ | **距離/分群細節 + magic-constant 預設分歧** | DIST-3 / HEADERS-4 (CONFIRMED) | ◽ 低 | dead config 旗標(no-op)、metric 誤標(Jaccard=SMC)、`min_common_coverage` 3 vs 5 預設打架 |
| **T6** ◽ | **harness 工具 + 座標慣例無測試** | PYTOOL-3 / IOREAD-6 | ◽ 低 | queue 健康燈恆 blind、BAM/FASTA 1-based↔0-based 慣例無測試 pin |

**修正優先序**：T1（補解析測試，零方法學風險最高 ROI）→ T2（補 null refuse，防線補洞）→ T3（per-CpG 改 beta-binomial 或標 anti-conservative；這正是專案自己 memory 的 Fisher→beta-binomial 建議）→ T4（輸出欄名明示軸別 + 統一 collapse 口徑）→ T5/T6（清 magic constant + 補測試，逐步）。

---

## L2 — 細節（主題 × finding 表）

### T1 🔴 甲基化解析零測試（CONFIRMED HIGH）

| finding | file:line | 事實 | 量化/驗證 |
|---|---|---|---|
| **TESTS-1** | `tests/*`（grep `MethylationParser\|MatrixBuilder`=0 命中） | `MethylationParser.cpp`(290 行) 的 MM/ML 解析、反股 **倒序迭代** CpG 座標（碼內自標 "CRITICAL" 103-157）、5mC/5hmC 分支 **完全無單元測試**；現有測試只覆蓋「filter-to-zero」負路徑 | 對抗驗證三路嘗試反駁全失敗 → **CONFIRMED HIGH**。這段是專案 memory 自己標「歷史最易錯」的 strand-flip 邏輯 |

> **為何最該先修**：此碼是資料入口，錯了會**靜默污染每一個甲基化 call → 每一個 Δβ / 距離 / 分群**，且無號誌。補一組 hand-constructed MM/ML golden test（正股 + 反股 + 多 CpG + 邊界）是**零方法學爭議、最高 ROI** 的動作。

### T2 🔴 反捏造基礎設施有洞（CONFIRMED HIGH）

| finding | file:line | 事實（live 實測） | 量化 |
|---|---|---|---|
| **PYTOOL-1** | `scripts/fill_report.py:38,54` | `resolve({'auc':None},'auc')=(True,None)` → `found=True` 不進 `missing`；`fmt(None)='None'`、`fmt(nan)='nan'` | 🔴 §13-A「template+data 注入缺 key 必 refuse」**對 null 值失效**——數據是 null 時報告靜默出現字串 `None`/`nan` 而非拒絕 |
| PYTOOL-2 | `number_provenance.py` | validated/PI 路徑「零同層來源」→ exit 2，與自宣告 FAIL-OPEN 契約矛盾 | MEDIUM；行為其實偏嚴（擋），但與文件契約不符 |
| PYTOOL-4 | `number_provenance.py find_token` | 裸 substring 比對 → 捏造 metric 若為某真值子字串會誤判「有來源」 | MEDIUM 假陰性 |
| PYTOOL-5 | `number_provenance.py` | 抽 metric 時剝掉負號（unicode/ASCII minus）→ 方向敏感數字溯源失真 | LOW |

> 對一個把「數字誠信」當核心紀律（§13 三層）的專案，**B 層（gate）本身的漏洞**比一般 bug 更該補：PYTOOL-1 讓 A 層（構造防捏造）對 null 破功。修法：`fill_report` 把 `None`/NaN 視同 missing 一併 refuse。

### T3 ⭐ 統計校準：p 值系統性偏樂觀

| finding | file:line | 裁決 | 事實 + 量化 |
|---|---|---|---|
| **FISHER-1** | `PerCpgAsm.cpp:274-313` | CONFIRMED→MEDIUM | 每個 read×CpG binary call 當獨立 Bernoulli 餵 Fisher exact，**無 over-dispersion 處理**（全碼無 beta-binomial/quasi/GEE）。合成示意（ILLUSTRATIVE，非真資料）：ρ=0.3 read 內相依下，NULL 的 Fisher 名目 5% FP 膨脹到 **53–68%（~17–20×）**。真量級需 `methylation.csv` 跑 read-level permutation 對照（needs_rerun）。**與專案 memory 自己的 Fisher→beta-binomial+shrinkage 建議一致** |
| SIGLABEL-2 / FISHER-8 / REGION-3 / PYASM-2 | 多處 | UNCERTAIN/MEDIUM | region-level / genome-wide **無 multiple-testing 校正**（只有 PerCpgAsm 內部單 region 跨 CpG 有 BH-FDR）；`global_p=min(p_alt,p_hp,p_hp_family)` 取最小 p 未校正 |
| **FISHER-2 / SIGLABEL-8** | `GlobalTest.cpp:124` | CONFIRMED→LOW | `chi_square_p = fisher_ffh.p_value; // Placeholder`，欄名宣稱卡方 p 實為 Fisher MC p。**grep scripts/docs：0 下游消費** → 純 misleading 欄位，無實害（建議刪或實作卡方 CDF） |
| **PERMANOVA-3 / SIGLABEL-7** | `StructureTest.cpp:292` | CONFIRMED→LOW | PERMDISP 的 `anova_p = F>4?0.01:(F>2.5?0.05:0.1)` **硬寫死 3 段查表**，非 df-aware F-CDF/permutation（df 在 241-259 算了卻沒用）；此 p 會 0.5× 懲罰 heuristic_score |
| FISHER-3 | `LocalTest.cpp:285-321` | **REFUTED→INFO** | min-p selection 機制真，但「進 is_significant」的實害宣稱**經查為假**（`best_p` 未驅動 is_sig）→ 假陽性已濾 |
| FISHER-4 | `FisherExact.cpp:285-314` | （未列 HIGH）MEDIUM | RxC Monte-Carlo p 用 `n_extreme/n_samples`，強訊號時可得 **p=0**；標準應 `(b+1)/(m+1)`（R `simulate.p.value` 慣例）→ `-log10(p)=inf` |
| PERMANOVA-5 | 生產配置 | MEDIUM | PERMANOVA 僅 99 permutations → p 下限 0.01，與 genome-wide 不相容 |

> **共同主旨**：以上全部讓「顯著」偏多、p 偏小。**對描述性 ASM 存在性率（如 Fisher_Frac_Sig）影響量級可觀，對 TP/FP 方向、對已 concluded 的 NEGATIVE 結論影響小**（後者靠的是 read-level PERMANOVA + 跨樣本一致性，不靠單 CpG Fisher）。

### T4 ⭐ ASM 軸別與 Δβ 口徑（UNCERTAIN/MEDIUM；多為 known caveat 但碼層未明示）

| finding | file:line | 事實 |
|---|---|---|
| **ASM-2 / REGION-5** | `RegionProcessor.cpp:903-906`, `PerCpgAsm.cpp:244-249` | 預設把 `{"1","HP1","1-1"}→0`、`{"2","HP2","2-1"}→1` 合併 → signed Δβ 與 per-CpG Fisher 量的是 **HP1-family vs HP2-family = germline/ALLELE 軸**（somatic 子 haplotype `1-1` 被折進 family `1`）。**非 somatic-controlled HP 軸**（HP1 vs HP1-1）。與 memory `feedback_asm_allele_axis_baseline_confound` 記的 confound 一致——**碼忠實實作了會被 baseline allelic methylation confound 的軸**，但輸出欄名（`tumor_hp_signed_delta`）未明示軸別 → 下游易誤用 |
| ASM-3 | `MethylationParser.cpp:41-90` | C++ **5mC-only**（只解析 `C+m?`，`C+h?` 僅數 delta 推進 ml_offset）→ signed Δβ 非 any-mod max-collapse。**memory 說的「5mC+5hmC 雙列砍半 dup-bug」其實在 MSA Python 抽取層，不在 C++ binary** |
| ASM-1 | `RegionProcessor.cpp:987-1006` | signed Δβ 用 **cell-weighted pooled grand-mean**（每 read×CpG cell 等權）非 per-CpG 平均 → HP1/HP2 覆蓋組成不同時 Simpson-paradox 偏移。（UNCERTAIN：數學真但 provenance/impact 被部分反駁） |
| ASM-4 | `PerCpgAsm.cpp` | per-CpG Fisher 二值化 `raw>0.5` 硬切，與 binary_matrix 的 0.2/0.8（中間=missing）**口徑不一** |
| PYASM-3 / PYASM-7 | `03_step4` vs `33/34/49` | 餵同一篇論文的腳本 **Δβ口徑不一**（5mC-only vs 5mC/5hmC）；germline(HP1-vs-HP2) Δβ 被當 ASM 報告（baseline-confounded by construction） |
| PYASM-1 | `18_dual_axis_pivot.py:104`, `03_step4:157` | UNCERTAIN→MEDIUM：per-CpG Wilcoxon 把 ±600/1000bp 視窗內**空間自相關**的 CpG 當獨立配對 → p anti-conservative |

> **建議**：C++ 輸出欄改名/加註明示軸別（germline-family vs somatic-HP），並統一全 ASM 腳本的 5mC/5hmC collapse 口徑（擇一並文件化）。這是**口徑治理**問題，非單點 bug。

### T5 ◽ 距離/分群細節 + magic-constant 預設分歧

| finding | file:line | 事實（多數 live 確認） |
|---|---|---|
| **DIST-3** | `RegionProcessor.cpp:427` 寫入，`DistanceMatrix.cpp` **0 讀取** | `distance_use_binary` 旗標傳入 `distance_config_` 卻**從未被 metric 邏輯讀** → 設定 no-op（binary/raw 由 metric 型別 hardcode）。CONFIRMED |
| **HEADERS-4** | `Config.hpp:37`=3 vs `DistanceMatrix.hpp:23,49`=5 | `min_common_coverage`(C_min) 預設 **3 vs 5 分歧** → 算距離所需最少共同 CpG 數依「哪個 struct 預設生效」而不同。CONFIRMED |
| DIST-5 | `DistanceMatrix.cpp` | `include_unmeth=true` 的 "Jaccard" 實為 **Simple Matching Coefficient (Hamming)**，metric 誤標 |
| DIST-1 / DIST-2 / DIST-6 | `DistanceMatrix.cpp` | self-dist 硬寫 0（BERNOULLI 下違反 d(x,x)=0 對相異 read）；中間甲基化 CpG 對 NHD/JACCARD 靜默丟棄（與 no-coverage 混為 -1）；低重疊對 MAX_DIST fallback 可造假雙峰結構 |
| CLUSTER-1 | `HierarchicalClustering.cpp:511,533` | CONFIRMED→LOW：silhouette 選 k 時靜默丟單例群（`count_a==0`），系統性偏好較大 k（相對 scikit s=0 慣例） |
| CLUSTER-2 / CLUSTER-5 | `HierarchicalClustering.cpp` | Ward 用「平均成對平方距離」近似而非 Lance-Williams 遞推；`linkage_matrix.csv` 的 cluster id 存矩陣 slot index 非 scipy 累積 id（語意錯） |
| HEADERS-3 / HEADERS-5 | `Clustering.hpp` vs `Stats.hpp` | 重複 `ClusterStats` struct（ODR 地雷，目前 Clustering.hpp 無人 include=dead）；metric 預設 Config=BERNOULLI vs DistanceConfig=NHD 分歧 |

### T6 ◽ harness 工具 + 座標慣例

| finding | file:line | 事實（live 確認） |
|---|---|---|
| **PYTOOL-3** | `harness_health.py:310` | queue 燈只數 `status=='queued'`，真實 queue status=`rejected:7/closed:6/pending:5/adopted_annotation:3/adopted:2`、**queued=0** → 燈恆 blind，**5 個 pending 永遠不可見**。CONFIRMED |
| PYTOOL-6 / PYTOOL-7 / PYTOOL-9 | `A3_5features_TP_FP_analysis.py` | NG∈{0,1} 被誤標 'NG2'（`v<=2`）；sign-consistency 把 NaN 缺值當中性方向 0 懲罰；AUC bar y 軸 floor 0.3 視覺裁掉強 anti-discriminative(AUC<0.3) 特徵 |
| IOREAD-6 | `tests/*` | BAM 1-based fetch ↔ FASTA 0-based fetch ↔ `ref_offset` 的座標慣例**無測試** pin（IOREAD-1 的「+1bp shift」具體宣稱經查**不成立**：同一 `region_start` 同時餵 BAM+FASTA 且 FastaReader 已正確 0-based，內部一致） |
| IOREAD-2 | `LohBedAnnotator.cpp` | `find_overlap` 只查最近-by-start 的單一區間 → 漏 overlapping/nested LOH 區間 |
| IOREAD-7 | `RegionProcessor.cpp` | `require_mm_ml` 硬鎖 true（CLI 不可覆寫）→ 無 MM/ML 的 read 靜默全丟 |

---

## L3 — 溯源（reproducibility）

- **方法**：12-reviewer Dynamic Workflow（`wf_ea727f9a-ef5`，28 agents，2.99M subagent tokens，497 tool uses，~15min）跨模組靜態深讀 → 每個 CRITICAL/HIGH finding 由獨立 skeptic 對抗驗證（嘗試反駁，預設懷疑）→ 主迴圈用 **live 真值**量化。
- **覆蓋**：C++ **全窮盡**（23 `src/core/*.cpp` + io/utils/main + 15 `tests/*.cpp` + 30 `include/core/*.hpp`）；Python **策略性子集**（論文數字 ASM 腳本 + 反捏造工具，**非全 181** — ⚠ PARTIAL flag）。
- **溯源檔**：
  - `InterSubMod/docs/methodology/20260610_code_audit_assets/workflow_findings_parsed.json`（105 findings 全文 + 對抗驗證裁決）
  - `InterSubMod/docs/methodology/20260610_code_audit_assets/quantification.json`（本 session live 驗證真值）
- **數字誠信**：本報告每個量化數字皆本 session live 讀回（fill_report 實測、queue Counter、grep 計數、合成示意明標 ILLUSTRATIVE）；分析與撰寫**分批執行**（§13.0）。合成 over-dispersion 示意**非真資料 effect**，真量級標 needs_rerun。
- **HEAD**：`bbcc7a9`（2026-06-10）。

---

## 附錄 A — 已澄清的「非問題」（memory stale 修正；§13.7 不信任記憶自驗）

| memory 宣稱 | 現行碼實況 | 裁決 |
|---|---|---|
| `PerCpgAsm.cpp:6` 誤引「De Waele 2020」應為 Orjuela | 現行 = `DAMEfinder (Orjuela 2020)`；grep「De Waele」**0 命中** | ✅ **已修**（commit 891e04b），非 open finding |
| expected_coverage hardcoded 75.0，KDE 未啟用 | `Config.hpp:56` 預設 `0.0=auto KDE`；`RegionProcessor.cpp:687` Pass 2 已呼叫 `estimate_diploid_coverage()`；75.0 = Pass-1 placeholder「overwritten by Pass 2」 | ✅ **KDE 已 wired**（HEADERS-1 對抗驗證 REFUTED）；殘留只剩 `compute_coverage_multiple(...=75.0)` 的 default-arg latent foot-gun（LOW） |
| CramersV reliability gate 可疑 | `GlobalTest.cpp:112` 用 expected≥5 cell 比例的 **Cochran 規則**（標準統計實務） | ✅ **方法 sound**；唯一議題是 summary(gated) vs significance.json(raw) **雙口徑報告**（REGION-1 CONFIRMED→LOW；effect 已在 memory 量化=762 MISSED 中 487 latent gated-CV=0 但 PERMANOVA 顯著） |
| 5mC+5hmC 雙列砍半 dup-bug（C++） | C++ `MethylationParser` **5mC-only**（`C+m?`），dup-bug 在 MSA **Python** 抽取層 | ✅ 釐清：C++ 層無此 bug；真議題是 C++ **忽略 5hmC**（口徑 caveat，T4 ASM-3） |
| 3 個 HIGH 假陽性 | FISHER-3（best_p 未進 is_sig）/ SIGLABEL-6（OR 且他處有文件）/ HEADERS-1（Pass-1 placeholder） | **REFUTED**——對抗驗證有效濾掉 |

## 附錄 B — 修正建議優先序（含 SWOT 速覽）

| 序 | 動作 | 類型 | 為何先做 | 風險 |
|---|---|---|---|---|
| **1** 🔴 | TESTS-1：補 MM/ML 解析 golden test（正/反股 + 多 CpG + 邊界） | 純測試，非 C++ 邏輯改 | 入口碼錯=污染全鏈；零方法學爭議、最高 ROI | 幾乎無（只加 `tests/`） |
| **2** 🔴 | PYTOOL-1：`fill_report` 把 None/NaN 視同 missing → refuse | Python | 補 §13-A 防線洞，與專案核心紀律對齊 | 低；可能讓既有「合法 null」報告需顯式處理 |
| **3** ⭐ | FISHER-1：per-CpG 檢定改 beta-binomial / read-level permutation，或維持 Fisher 但輸出明標 "anti-conservative, independence-assumed" | C++ 方法（需 `/methodology-audit`→`/cpp-change`，Hard Gate） | ASM 存在性率校準；專案自己已建議 | 中；改檢定需重跑量化 + 跨樣本回歸驗證 |
| **4** ⭐ | T4：C++ ASM 輸出欄名明示軸別（germline-family vs somatic-HP）+ 統一 5mC/5hmC collapse 口徑 | C++ + Python 口徑治理 | 防下游把 germline-confounded Δβ 當 somatic ASM | 中；牽涉多腳本口徑對齊 |
| **5** ◽ | T5/T6 清理：DIST-3 dead flag、HEADERS-4 `min_common_coverage` 3↔5 統一、PYTOOL-3 queue 燈、FISHER-2/PERMANOVA-3 placeholder p、metric 誤標 | 雜項 | 降未來 drift；多為一行修 | 低；逐項可驗 |

> ⚠ **任何 T3/T4 的 C++ 改動**：先走 `/methodology-audit`（OPTIONS+SWOT）→ 用戶選方案 → `/cpp-change` 6 步（C++ commit 是 Hard Gate，必編譯）。**本報告不自行改 C++**。
> T1（測試）/T2（Python）可獨立先行，不阻塞研究主軸。

---

**一句話收尾**：程式碼方法學**沒有致命錯誤、無結論翻轉風險**；最該補的是「最易錯的甲基化解析竟無測試」與「反捏造工具自己有洞」這兩個基礎，加上把 per-CpG Fisher 的樂觀偏誤講清楚/改掉。其餘是口徑與 magic-constant 治理，可漸進。
