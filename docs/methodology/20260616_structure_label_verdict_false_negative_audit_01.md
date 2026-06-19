<!--
建立時間: 2026-06-16
問題類型: 統計方法 + 判定架構（VerificationClass verdict）
影響 track: paired（主軸）/ 全 ISM region-level
狀態: pending_decision
build_branch: research/subclonal-reconstruction-202606 @ b761336
data_sources: docs/methodology/_assets/20260616_fn_verdict_readback/HCC1395_paired_full.json
upstream_audit: workflow wf_fe623383-64b（23 agents，8 genuine FN / 15 機制，3 推翻）
-->

# ISM「結構↔標籤關聯」判定的大規模假陰性審查（VerificationClass verdict）

> **partial flag**：量化為**單樣本 HCC1395 paired_full（TP 29,754 / FP 627 region）、單 pipeline** → 結論封頂 **⭐3**；跨 6 樣本 scope-out 待做。

## 0. TL;DR（BLUF）

ISM 把 **91.7%** 的 region 判成 `Weak`/`Noise`，但其中 **Noise 74.3%、Weak 98.0%** 早已帶有一個 **valid 且顯著的 PERMANOVA**，只是最終 `VerificationClass` 的 2×2 判定**從不引用 PERMANOVA**。
**根因不是實作 bug，而是「結構↔標籤關聯」的定義不清（把 D1/D2/D3 三個不同問題混用）＋ verdict 圍繞的定義（D2）不符主軸目標（D1）。** 機械漏洞 FN1–FN15 是症狀。

## 1. 問題描述（親讀原始碼坐實）

### 1.1 verdict 邏輯（`src/core/SignificanceAnalyzer.cpp:309-339`，verified_in_file）

```cpp
cluster_significant = result.passed_gate &&
    (result.global_alt.fisher_ffh.p_value <= 0.05 ||
     result.global_hp_family.fisher_ffh.p_value <= 0.05);   // D2：離散列聯表 + gate
label_sig = result.label_significant;                        // = LabelTest.any_significant（D1′ 純量/單尾）
// 2×2：
//   label_sig &&  cluster_significant → "Strong"
//  !label_sig &&  cluster_significant → "Subclone"
//   label_sig && !cluster_significant → "Weak"
//  !label_sig && !cluster_significant → "Noise"
```

`result.label_hp_permanova`（多變量 label-PERMANOVA，第 **301** 行 `run_grouped_structure` 算好並寫入）、cluster PERMANOVA、校準 Δβ —— **均未被第 309-312 行引用**（Δβ 於本檔 0 次引用）。→ **算了卻丟棄**。

### 1.2 gate 邏輯（`src/core/GlobalTest.cpp:138-142`，verified_in_file）

```cpp
p_pass = (fisher_ffh.p_value <= gating_p_threshold);   // 預設 0.1
v_pass = (cramers_v        >= min_cramers_v);           // 預設 0.1
passed_gate = p_pass && v_pass;                         // 雙重 AND
```

`passed_gate` 是**串接前置**：cluster-first PERMANOVA 與 local test **只在過 gate 後才跑**（`SignificanceAnalyzer.cpp:94,99`）；`cluster_significant` 也要求 `passed_gate`。→ gate 的 type-II 被整條鏈繼承。

## 2. 量化影響（第一手讀回，非捏造）

**來源**：`docs/methodology/_assets/20260616_fn_verdict_readback/HCC1395_paired_full.json`（腳本 `scripts/analysis/fn_verdict_readback_audit.py`，純離線 group-by，零改碼／零重跑）。

### 2.1 階梯1 — 各關篩除（TP 29,754 region）

| 關卡 | 數 | 佔比 | 旗標 |
|---|---|---|---|
| 最終 Weak+Noise（非 Strong/Subclone） | 27,275 | **91.7%** | 🔴 異常大 |
| FN2 全域 gate 擋掉（PassedGating=false） | 26,925 | **90.5%** | 🔴 異常大 |
| └ 連帶 ClusterPERMANOVA 從未跑（=gate 數，完全相等 26,925） | 26,925 | 90.5% | 串接坐實 |
| CramersV 匯出欄被歸 0 | 29,003 | 97.5% | 呈現層（見 §3.4） |
| NumReads<10 / <20 / NumCpGs<10 | 0 / 0 / 0 | 0% | ✅ 證實 ±5000 無 CpG 飢餓（FN6 推翻） |
| ClusterDispersionWarn=true | 0 | 0% | ⚠ 欄位疑未填 → 無法排除離散度假象 |
| Potential_LOH=true | 14,116 | 47.4% | — |

VerificationClass：Weak 19,435 / Noise 7,840 / Strong 1,940 / Subclone 539。

### 2.2 階梯2 — verdict 忽略的顯著 PERMANOVA（FN1，決定性）

| 子集 | 帶 valid 顯著 PERMANOVA（被 verdict 忽略） | 佔該類 | 幾全在 label 側 |
|---|---|---|---|
| Noise 7,840 | **5,824** | **74.3%** | labelside 5,810 |
| Weak 19,435 | **19,037** | **98.0%** | labelside 19,034 |

FP 對照（627 region）：Noise 113/207 = **54.6%**、Weak 393/405 = **97.0%** → 此為**結構存在性**問題，**非 TP/FP 判別**（與 filter-DEAD 一致）。

### 2.3 跨 run 旁證（report_claim / 不同 run，標級引用）

- 全基因組 2,451 假陰性（94% p_fail + 6% v_fail，HP-AUC median 0.749）— `report_claim`（worktree commit 724aa51 + 可能已失效 `/tmp` CSV，未能重 grep）。
- reclassify pilot（Python）救回 **46.6–53.2%**，`orig_Strong_changed=0` — `report_claim`，**仍 Python pilot 未落 C++**。
- MAX_DIST→SKIP 全基因組 Significant 1,276→11,912 — `verified_in_file`（`docs/experiments/in_progress/2026/06/20260614_u1_maxdist_vs_skip/wg/comparison_summary.json`）；**現編譯預設已 SKIP（2026-06-14 起）→ FN4 僅遺留輸出殘留**。

### 2.4 T1 離線 A/B 結果（2026-06-16，第一手；`_assets/20260616_fn_verdict_readback/HCC1395_reclassify_t1.json`，腳本 `scripts/analysis/fn_verdict_reclassify_t1.py`）

實際重算 `verification_class`（`label_sig'=label_sig OR 顯著 label-PERMANOVA`，**排除 sample 軸**）：

| 指標 | TP | FP |
|---|---|---|
| Noise→Weak 救回（有真 label-PERMANOVA） | **5,810/7,840 = 74.1%** | 113/207 = 54.6% |
| **orig_Strong 被降級** | **0** ✅ | 0 ✅ |
| 救回主軸 | **allele 5,605** ≫ HP 1,119 | allele 108 ≫ HP 47 |
| 僅 sample 軸顯著（FN10 confound，已排除） | 1,404 | 82 |
| Subclone→Strong（label 補上） | 526 | 1 |
| 額外接 cluster-PERMANOVA：Noise→Subclone/Strong / Weak→Strong | 92 / 258 | 2 / 7 |

**判讀（GO/NO-GO）= FLAW CONFIRMED**：flip 74.1% 且 `orig_Strong_changed=0`（純加、不傷既有 Strong）→ Option B 安全。**救回主軸是 allele（ALT-vs-REF）非 germline-HP**，比預期更貼近 somatic/cis 訊號。
🔴 **未解邊界**：`Label*DispersionWarn` 欄全 false（疑未啟用）→「dispersion-clean 5,810」**無資訊**，無法從 CSV 排除離散度假象；allele 軸顯著需 **cis-control** 分 cis-ASM vs subclone。→ HTML 肉眼 / dispersion-aware 重算的落點。

## 3. 根因分析：方法學 bug／不符目標／定義不清？

### 3.1 「結構↔標籤關聯」被混用的三個定義

| 定義 | 問的問題 | 對應檢定 | verdict |
|---|---|---|---|
| **D1 omnibus** | 標籤能否解釋 read×CpG **距離結構** | label-PERMANOVA（多變量） | ❌ 算了棄用 |
| **D1′ 純量** | 同標籤距離是否較小（單尾） | `label_sig` | ✅ 用（弱） |
| **D2 重合**（最嚴） | **非監督幾何分群**離散化後與離散標籤對得上＋V≥0.1 | `cluster_significant`/gate | ✅ 當 gate |
| **D3 存在** | 是否存在 ≥2 離散亞群 | cluster PERMANOVA / silhouette | ❌ 算了棄用 |

### 3.2 為什麼 90.5% 過不了 gate / 91.7% 落 Weak-Noise

- **D2（gate）失敗 90.5%**：silhouette 最大化的是**幾何分離**（甲基資料常被覆蓋度/CpG pattern 主導），其主軸通常**不是標籤** → 離散列聯表對不上 → gate fail。
- **D1′（`label_sig`，單尾純量）** 比多變量 PERMANOVA 弱 → 漏掉多變量/離散度型關聯。
- 結果：真有 **D1 關聯**（label-PERMANOVA 顯著）但「幾何主軸≠標籤（D2 fail）＋純量單尾沒跨（D1′ fail）」的 region → 落 **Noise**。量到 **74.3% Noise** 即此。

### 3.3 判定

- **不是 (a) 實作 bug**：程式忠實照它寫的做；FN1–FN15 是症狀。
- **主因 (c) 定義不清（最深）**：「結構↔標籤關聯」從沒被釘成單一定義；verdict 默挑 D2（最強重合）當 gate、D1′（弱純量）當救援、棄用 D1（自然 omnibus）。**74–98% 落差＝這個定義模糊的量化大小。**
- **主因 (b) 方法不符目標**：主軸＝subclonal reconstruction（甲基被 characterize、有界），對應 **D1（germline ASM 表徵）**；verdict 卻圍繞 **D2**，對 ASM 表徵過嚴、對稀有 somatic subclone 又被 floor/peel-off 擋掉 → 兩目標都不貼。

### 3.4 心智模型更新（FN5/FN6 推翻；FN3 爭議未解）

- **FN3 CramersV 可靠性歸零 → DISPUTED（未解，勿當定論）**：workflow 親讀主張 gate 用 **raw 未歸零 V**（歸零只動匯出欄與次要 `is_significant` binary）；第一手旁證 **2,106 region 匯出 V=0 卻仍 PassedGating=true**（至少對這批 gate 不 key 匯出欄）。**但這與 `project_asm_locus_display_and_cramersv_reliability_gate` 的 script 92/93 發現衝突**（該處稱 2,538 個 raw CramersV≥0.3、median 0.81 的位點被可靠性閘控歸零、根因單一=reliability gate）。→ **乾淨定案需 emit 原始未歸零 V（T2 選配一行 C++）**；在此之前**不可斷言 CramersV gate 是/不是決定層 FN**。可確定者：487-latent 的 label-PERMANOVA 走 **ungated 路徑**（`SignificanceAnalyzer.cpp:301`），其「被忽略」屬 **FN1**——與 gate 歸零是**兩個可並存的獨立機制**，非互斥。
- **FN5 C_min/binarization → 推翻**：0.2/0.8→MISSING 只動 `binary_matrix`（NHD/Jaccard 非預設）；預設 BERNOULLI 用 `raw_matrix` 連續、只降權。
- **FN6 ±1000 window → 推翻為 footgun**：生產實跑 ±5000（`run_vcf_all_snv.sh` MODE=all-with-w5000）；§2.1 NumCpGs<10=0 坐實。

## 4. 修改選項 + SWOT

> ⚠ verdict / gate 邏輯在 `src/core/`，任何 B/C/D 落地需走 `/cpp-change` 6 步 PDD（pre-commit 編譯 Hard Gate）。但**決定性的離線 A/B 量化（T1）零改碼**先做。

### 方案 A：不修改（接受現狀 + 文件化）
補文件說明「`Noise`≠無結構」、PERMANOVA 欄供下游自取。

|  | Helpful | Harmful |
|---|---|---|
| **Internal** | S：零風險、不動 Strong | W：91.7% 假陰性留存，論文若引「ISM 判別力」會被低估 |
| **External** | O：資源投他處（normal-HP gap2） | T：reviewer 質疑「為何算了 PERMANOVA 卻不用」 |

### 方案 B：verdict 接入 label-PERMANOVA（rescue lane，對齊 D1）
`label_sig' = label_sig OR (label_hp_permanova.valid && p≤0.05 OR label_allele_permanova.valid && p≤0.05)`。對應 reclassify pilot（救回 46–53%、不動 Strong）。**先 Python 離線 A/B（T1）證實 verdict-changing 比例 + 確認非離散度/germline 假象，再 cpp-change。**

|  | Helpful | Harmful |
|---|---|---|
| **Internal** | S：救回大量假陰性、PERMANOVA 已算好（低成本）、orig_Strong 0 動 | W：救回多為 germline HP 軸，對 somatic subclone 淨增益未定（gap2）；新增類別需重訓下游消費者 |
| **External** | O：對齊「ISM=無監督距離 PERMANOVA」論文口徑 | T：若不先驗離散度/confound，恐把假陽性也救回（D1 多變量會被 dispersion 膨脹） |

### 方案 C：解耦 gate + 三軸並陳（D1/D2/D3 具名輸出）
不讓離散 contingency gate PERMANOVA；verdict 改報三個具名軸（label-associated / discrete-subclone / both），不塌縮成單一 categorical。

|  | Helpful | Harmful |
|---|---|---|
| **Internal** | S：定義透明、各軸可分別驗證、終結 D1/D2/D3 混用 | W：改動最大、需改 schema 與所有下游（SubcloneAnalyzer/CSV/報告） |
| **External** | O：方法學最乾淨、可寫成 paper 方法貢獻 | T：schema 變動衝擊既有 6 樣本輸出與既成結論的可比性 |

### 方案 D：對齊目標先做 normal-HP 對照（gap2），再定 verdict
為 subclonal reconstruction 目標，先以 normal-only HP 對照分離 germline vs somatic，**再**決定 verdict 報什麼。**研究 scope 決策，非純 code change。**

|  | Helpful | Harmful |
|---|---|---|
| **Internal** | S：直擊「真 somatic subclone vs germline ASM」的根本未決 | W：需 COLO829 normal 甲基（目前 5/6 ready，COLO829 缺）；週級工作 |
| **External** | O：解鎖 ⭐3→⭐4 跨樣本 somatic 主張 | T：normal 資料 SPOF（zhenyu112 帳號）、時程風險 |

## 5. 驗收標準

- [x] **T1 離線 A/B（零改碼）✅ 2026-06-16 完成**：Noise→Weak 救回 **74.1%**、`orig_Strong_demoted=0`、救回 **allele-dominated**；FN10 sample-only 1,404 已排除。見 §2.4。🔴 dispersion 欄未填 → 離散度假象未排除。
- [ ] **離散度/confound 把關**：對 flip 集確認非 `dispersion_warning` 假象、非純 tumor-vs-normal confound（sampleOrphan 要排除）。
- [ ] 若進 B/C：`test-quick` 通過 + `test-full` F1/分類分布不回退 + 跨 ≥3 樣本方向一致。
- [ ] 跨 6 樣本 scope-out：機制在 k/6 樣本復現（per-sample tally，**禁 pool**）。

## 6. 用戶決策（pending）

**選擇**：[ ] A 不改 / [ ] B 接入 label-PERMANOVA / [ ] C 三軸並陳 / [ ] D 先做 normal-HP / [x] **T1 已跑（2026-06-16）→ FLAW CONFIRMED（見 §2.4），B 證實安全（0 Strong 降級）**；改前建議先 HTML 肉眼／dispersion 重算排除離散度假象 + allele 軸 cis-control
**日期**：
**理由**：

---

## 7. cluster-first 循環修正：cluster×label 對應 typology + bootstrap-stability k（2026-06-18 收斂，方向 B+D）

> 本節是 §1–§6 的**延伸**：§3 指出舊 verdict 用循環的 cluster-first（D2）當門檻。這裡定義「如何把循環換成可信的對齊判定」，並把「標籤內裂成 subclone」明確納入。**用戶 2026-06-18 拍板方向 = B + D。**

### 7.1 cluster-first 為何循環、何時才有害

silhouette 先挑「最分得開的切法」→ PERMANOVA 再測「分得開嗎」→ 必然顯著，p 偏樂觀（正規 null 需每次置換**重新分群**，程式未做 → anti-conservative）。**但要分清哪題受害：**
- **對齊題**（標籤與甲基是否相關）：循環是「答錯題」——應改用 **label-first PERMANOVA**（外部標籤 HP=germline SNP 定相 / ALT-REF=變異等位；**免 k、非循環**）。
- **次結構題**（一個標籤內是否裂成多 subclone）：**這時才需要 k**；raw silhouette 的「過度切」偏誤會把雜訊切成假 subclone → **必須用穩定度選 k**（= D）。

### 7.2 三種 cluster×label 對應 = 三個生物結論（B）

| 對應形狀 | 白話 | 生物意義 | 現況 |
|---|---|---|---|
| 1 群 ↔ 1 標籤（雙射，ARI≈1） | 對齊 | 乾淨 ASM／標籤即主軸 | ≈ Strong |
| 1 標籤 → 多群 | 標籤內有次結構 | **subclone／第二表觀軸** | 平行 session Stage② within-HP 已部分做（WithinHP_CleanMultigroup） |
| 多標籤 → 1 群 | 標籤沒對上 | 主軸非標籤（覆蓋度等）／無對齊 | ❌ **未明確標記**（目前混入 Subclone/Noise） |

對齊強度用 **ARI（Adjusted Rand Index）或 NMI**（自動校正巧合）量化：ARI≈1 對齊、ARI≈0 獨立。

### 7.3 bootstrap-stability k：為什麼一定要（用戶指定，清楚說明原因）

**問題**：要回答「一個標籤內裂成幾個 subclone」就得選 k，但 raw silhouette **必定**回傳「某個 k 的最佳切」——**即使資料是純雜訊也會切出群**。直接拿這個 k 當 subclone 數 = 把循環偏誤帶進 subclone 判定。

**修法 = consensus / bootstrap-stability clustering**：
- 重抽資料 N 次（bootstrap reads 或 CpG），每次**重新分群**，算 co-clustering 矩陣（兩條 read 多常被分在同群）。
- 只保留「重抽下**穩定還在**」的群（co-cluster 比例 ≥ 閾值如 0.9）；雜訊切出的假群在重抽下會散掉、被濾除。
- k 取「**穩定群數**」而非 silhouette 最佳 k。

**為什麼非做不可（兩個原因）**：
1. **去循環**：穩定度把「優化挑最開」的樂觀偏誤過濾掉——真群在重抽下重現，假群不會。這是把循環偏誤擋在門外的唯一乾淨手段。
2. **subclone 是稀有低頻訊號、最易被假群淹沒**：不穩定度把關，等於直接決定「裂出來的群是真 subclone 還是雜訊過切」。

**對應碼**：`enable_bootstrap = false`（現況關閉，見 `RegionProcessor.cpp:1717`）→ 改 `true` + 接 consensus；**within-HP 子分群尤其需要**。

### 7.4 🔴 誠實可行性更正：B 非「零改碼純離線」

實查（2026-06-18）：`significance_summary.csv` 為逐 region 摘要、**無逐 read 群指派**；canonical 查無 `reads.tsv`/clustering csv；cluster-PERMANOVA 被 gate 擋到僅 ~10% region 有跑；只匯出 `p`+`CramersV` 摘要、**未匯出列聯表本身** → **真 cluster×label ARI 無法從現有 CSV 算**。
→ **更正**：realistic B = ISM 補**小 emit**（每 region 的 cluster×label 列聯表 + `optimal_k` + 每群標籤組成，**ungated**；D6-class provenance，**無演算法改動**）→ 之後 offline ARI/typology trivial。即「B 純離線零改碼」**更正為「小 emit（cpp）+ offline」**。

### 7.5 確認要修的清單（+ 原因）

1. **verdict 接 label-first**（§4 Option B；已部分落地 reclassify-v2）— 對齊題的**非循環**判定，免 k。
2. **cluster×label 對應 emit**（ungated 列聯表 + optimal_k + 群標籤組成）— **讓三 case 可算**；與既有「ISM provenance 輸出 pending」(FN15/FN7) 同 class。
3. **bootstrap-stability k**（`enable_bootstrap=true` + consensus，`RegionProcessor.cpp:1717`）— within-label subclone 裂群**去循環 + 可信**。
4. **ARI/NMI 對齊指標 + 三 case 具名 typology verdict** — 取代單一 p，並補上「多標籤→1 群（沒對齊）」這個目前沒被標的 case。

### 7.6 驗收 / next

- ② emit 落地後：offline 算 30,490 region 的 (aligned / 標籤內次結構 / 沒對齊) 三 case 佔比 + ARI 分布。
- ③ bootstrap 落地後：within-HP subclone 群數在重抽下穩定率 ≥ 閾值的比例；對照 raw-silhouette 過切量。
- 全程走 `/cpp-change` 6 步 PDD（pre-commit 編譯 Hard Gate）。

### 7.7 實作現況盤點（2026-06-18 親查 `feat/summary-nreadsvalid`，read-only git）

**③ bootstrap-stability k — 機具 90% 已在，沒接到該接處**

| 狀態 | 證據 |
|---|---|
| ✅ `Bootstrap` class + ARI(`adjusted_rand_index`) | `src/core/Bootstrap.cpp` |
| ✅ 已接**主分群**、每次重抽重新分群（正規 null） | `SignificanceAnalyzer.cpp:168 bootstrap_->compute(...,cluster_func)` |
| ❌ `enable_bootstrap=false`（從不跑） | `SignificanceAnalyzer.hpp:37` |
| ❌ within-HP 子分群**沒用 bootstrap**（raw silhouette + heuristic 閘 sil≥0.5+balanced+gap>0.15） | `RegionProcessor.cpp:1180-1188`、`compute_within_hp_substructure :2079` |

🔴 within-HP 過切鐵證：loose(sil>0.3)=**98.25%** vs clean=**14.2%**（`within_hp_clean_pilot.json`@4ec339b），silhouette 多落 0.6-1.0 → 98% 不可能全真 subclone → **非接 bootstrap 不可**。

**② cluster×label ARI / typology — 半成品**

| 狀態 | 證據 |
|---|---|
| ✅ `optimal_k` emit / within-HP 群數 emit / ARI 函式存在 | significance.json `:2523`、CSV `:1338`、Bootstrap.cpp |
| ❌ **ARI(幾何群,標籤) 從沒算**（ARI 只用於 bootstrap 比兩次分群） | grep 僅 `Bootstrap.cpp:154` |
| ❌ cluster×label 列聯表沒 emit；三 case typology（含「沒對齊」）只在 Python pilot | `reclassify_v2_pilot.py` |

**within-HP 三 case 初值（驗證自 committed pilot @4ec339b，n=30,490；🔴 /tmp 原始已失 = committed-summary 級）**

| case | 對應 | 數量/佔比 | caveat |
|---|---|---|---|
| 1↔1 對齊 | Strong | 22,311 (73.2%) | reclassify_v2 |
| 1 標籤→多群 | within-HP substructure | clean **14.2%**(4,322)/loose **40.7%**(12,416)/C++ clean_multigroup **26.5%** | 🔴 隨閾值漂 14→27→41% ＝循環現形；43% CN-altered |
| 多標籤→1 群（沒對齊） | StructureNoLabel120+MultiGroupNoLabel323 | ~443 (1.5%) | 目前未明確標 |
| noise | 三子型 | 4,369 (14.3%) | Chaotic/Uncorrelated/Uniform |

reclassify_v2 救回 **46.6%**(3,810/8,179)、Strong 誤動 **0**（對齊 T1 demoted=0）。

**最小剩餘 → cpp-change spec**：`InterSubMod/docs/plans/20260618_within_hp_bootstrap_and_cluster_label_ari_cpp_change_spec_01.md`。

---

### 數據溯源（§13）

| 數字 | 來源:鍵 | 狀態 |
|---|---|---|
| Weak 19435/Noise 7840/Strong 1940/Subclone 539；PassedGating=false 26925；CramersV=0 29003 | `_assets/20260616_fn_verdict_readback/HCC1395_paired_full.json:TP.stage1_per_gate_attrition` + `.verification_class_counts` | verified_in_file |
| Noise 5824/7840；Weak 19037/19435；labelside；CramersV=0 presentational 2106 | 同 JSON `TP.stage2_*` + `TP.T2_cramersv_zero_FN3` | verified_in_file |
| FP Noise 113/207；Weak 393/405 | 同 JSON `FP.*` | verified_in_file |
| verdict / gate 行號 | `src/core/SignificanceAnalyzer.cpp:301,309-339`、`src/core/GlobalTest.cpp:138-142` | verified_in_file（親讀） |
| 2451 假陰性 / reclassify 46-53% / MAX_DIST wg | workflow wf_fe623383-64b evidence_findings | report_claim / verified_in_file（見 §2.3 分級） |
