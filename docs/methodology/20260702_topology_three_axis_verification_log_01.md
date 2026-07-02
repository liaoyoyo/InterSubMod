---
title: 拓撲三軸/canonical 分析 — 迭代對抗驗證記錄（audit log）
date: 2026-07-02
status: verification_log
build_branch: research/subclonal-reconstruction-202606
method: 每輪 = 多 subagent 對抗質疑 + 獨立 python3 重算 + 主回合獨立覆核 + 修正/解釋 + 記錄;loop-until-dry
data_sources: docs/methodology/20260702_topology_cluster_shape_three_axis_analysis_01.md
---

# 拓撲三軸分析 — 迭代對抗驗證記錄

> 目的:對 `InterSubMod/docs/methodology/20260702_topology_cluster_shape_three_axis_analysis_01.md` 反覆「質疑→解釋/修正」直到 dry（連續一輪 0 個新真問題）。每輪記錄:質疑、覆核、我獨立重算、處置。

## 收斂狀態
| 輪 | 面向數 | agents | 真問題(覆核後) | major | 處置 | 狀態 |
|---|--:|--:|--:|--:|---|---|
| R1 | 8 | 33 | 21 | 3(同源) | 全修/解釋 | ✅ 已處置 |
| R2 | 4 | 17 | 7 | 1(殘留傳播) | 全修 | ✅ 已處置 |
| R3 | 傳播 sweep | — | 0(純傳播修正) | 0 | grep+HTML 驗 | ⏳ 待獨立 agent 確認 DRY |

**R2 裁決 NOT-DRY**:R1 核心修正全部確認落地正確(數字逐位重現、零捏造、核心裁決不變),但找到 **1 MAJOR + 6 MINOR 殘留** — 同一根因:R1 修對分析本體,未傳播到 display/provenance 4 層(caption/標題 + PI-facing HTML + provenance 引用 + JSON 命名)。R3 = 純傳播 sweep(無新計算)。

---

## Round 1（workflow wf_4106c85d-29a,8 面向 + 覆核,33 agents / 2.77M tok）

### 總裁決:NEEDS_FIX（範圍明確，非資料/生物問題）
所有 headline 數字獨立重算逐位吻合、**零捏造、零計算錯、核心裁決不變**（branched%=幾何上界、CN 未控、partly_artifact/L3、全基因組樹定不出來）。唯一實質缺陷 = 「一致性 100% 有效樹/可放心當定論」指標組。

### 🔴 MAJOR（3 項同源，主回合已獨立重算確認）
**M1 「成環 0 / 100% 有效樹 / 可放心當定論」= 同義反覆且與報告自身矛盾**
- 覆核確認:`cyclic_invalid` 算在**已 transitive-reduce 的 acyclic edges** 上 → AHU `(X)` 結構上不可能出現 → 對任何輸入恆 0。但報告自列 **1334 incompatible**（四配子/perfect-phylogeny 違反）。
- 🔬 **主回合獨立重算**（§13.0）:incompatible 合計 **1334/22676 = 5.9%**;逐樣本 **COLO829 0.4% / HCC1954 1.4% / DORADO 2.3% / HCC1937 2.7% / HCC1395 3.0% / H1437 5.5% / H2009 19.1%**;`has_cycle` 合計 **348**。→ 誠實違反率跨樣本差 0.4%↔19.1%,被「100%」反轉抹平。
- 處置:枚舉完整性（11 形狀·0 未分類=正確）與 perfect-phylogeny 有效性**分列**;移除「100% 有效/可放心當定論」;改報真實 incompatible% + has_cycle。

**M2 c≤k+1 off-by-one → 「c>k+1=0」vacuous**
- 覆核確認:c=含≥1 ALT 向量數（root 不計）;觀察到 germline root 時正確界為 **c≤k**。max(c−k)=1 → 「c>k+1」永不觸發。
- 🔬 主回合重算:c==k+1 且觀察到全-R root 的真違反 = **59 區**（HCC1395 30 為主）。「c>k+1=0」是 vacuous check。
- 處置:改用 c≤k（root 觀察時）;停止把「c>k+1=0」當 perfect-phylogeny 驗證。

**M3 valid_pct=100% 非獨立驗證（false confidence）**
- 覆核確認:兩分子（cyclic、c>k+1）結構性恆 0 → valid_pct 恆 100%,卻當 🔴 headline。處置同 M1:改用會變動的 incompatible/has_cycle 為分母。

> 三者=同一 consistency 指標的三面（tautology / off-by-one / 偽驗證）,屬**措辭+指標定義 overstatement,非錯誤數字/錯誤生物**。

### MINOR（18 項,全 confirmed_real;摘要,詳見 workflow output）
數字口徑:① sSNV 密度帶 2.09–2.41 僅 5/7 樣本（真全距至 H2009 17.41）② 「c≥3%」欄是 canonical 殘差非 count(c≥3)/total ③ Pearson 0.867 未標 n=7 leverage + 實為 br_sib% ④ CN「47.9 點」以 n=65 neutral 為錨脆弱,穩健錨=gain(1189)67.6 vs LOH(483)27.5。
定義/枚舉:⑤ §B topology_type 漏 incompatible(267×)列了 0 例 star/mixed ⑥ determinacy 漏 `other`(18.4%,4164;主回合確認)列了 0 例 A_noisy ⑦ tree_shape 漏 single_nested ⑧ chr8 範例同時計入 branched%+confirmed% 而 determinacy=incompatible,教學需標「局部有結構全域 incompatible」。
中間檔/程式:⑨ canonical_shapes.json 把 incompatible 區標 germline_only(應 incompatible/no-tree) ⑩ genotype_cap=8:H2009 17.7% 最受限 ⑪ HTML Math.round vs 報告一位小數 ⑫ HTML 硬編碼 CramérV/CN/「11種」有 live-vs-literal 漂移風險 ⑬ 2 處行號誤引(branched=107+111 非110;n_sSNV≥2 filter 位置)。
措辭/理論:⑭ TL;DR「c=2 唯一自由度」literally false（c≥3 也有,3.6%）⑮「幾何上界/下界/缺口」暗示夾住單一真值,但 confirmed% 非平行下界（4/7 缺口為負）⑯「11 種非 n^n」category mismatch（11=unlabeled AHU,n^n=labeled;真驅動=經驗淺薄 96.4% c≤2）。

### 已排除假警報（1）
dim3「canonical 非 byte-reproducible / 只 10 形狀 / HCC1395 single 1985 vs 1995」= **FALSE**。覆核:獨立 AHU 對 7 樣本 byte-for-byte 復現（11 形狀、(())=11490）。審計自己重建用錯模組的 has_cycle 欄排除了 acyclic 區才造成 1985+12=1995 差。number-reproduction 無缺陷。

### R1 處置清單（本輪落地）
1. 必修 M1/M2/M3 → 改寫報告 canonical 段 + 數據總表一致性欄 + canonical_shapes.json consistency + HTML badge。
2. minor → 定義表補 incompatible/other/single_nested;sSNV 密度全距;c≥3% relabel;CN gain-vs-LOH 錨;Pearson n=7 caveat;TL;DR/幾何措辭/n^n 改寫;HTML toFixed;行號修;canonical relabel。

---

## Round 2（workflow wf_966b37ef-46c,4 面向 + 覆核,17 agents / 1.4M tok）

### 裁決:NOT-DRY（找到 1 MAJOR + 6 MINOR 殘留,皆傳播層）
R2 從 7 份 raw JSON 各自獨立重算,**確認 R1 核心修正全部正確**:incompatible 1334/22676=5.9%、逐樣本 0.4–19.1%、has_cycle 348、c==k+1-root 真違反 59、other 桶 18.4%、六桶總和 100.0、11 形狀 0 未分類、headline branched%/confirmed% 全逐位吻合。**核心裁決不變且成立。**

### 經覆核成立的殘留（主回合 Read 二次確認 + R3 已修）
- **[MAJOR-1]** 報告 line 40 blockquote 殘留 forbidden「一致性欄全 0=可放心當定論(100% 有效樹)」+ line 21 標題仍「一致性」+ 「一致性欄」已不存在(現 incompatible% 欄,值非 0)→ 與 line 68/73 自相矛盾。**R3 修**:line 40 改「枚舉完整可放心;有效性非 100%(0.4–19.1%)」;line 21 標題→「incompatible%（有效性）」。
- **[MINOR-2]** HTML STAT_DICT n_clusters 仍「c≤k+1、實測中位 2」(HCC1395 中位=1)。**R3 修**:c≤k(root 觀察)、中位改「多為 1-2」。
- **[MINOR-3]** CN 錨未傳播:報告 caveat line 112 + HTML line 328 仍「27.5%↔75.4% 擺 48 點」(與 §③ 40 點矛盾)。**R3 修**:改 gain 67.6 vs LOH 27.5=40 點。
- **[MINOR-4]** §A filter 行號誤引 sm_region_integration.py:141(=nested_edges 輸出)。**R3 修**:改 topology_analysis.py:141 + sm_region_integration.py:200。
- **[MINOR-5]** canonical_shapes.json 兩個 incompatible(per_sample no-tree 1.0% vs consistency 3.0%)易混引。**R3 修**:報告 line 72 明指採 consistency.incompatible_pct、勿引 per_sample no-tree。
- **[MINOR-6]** HTML「11 種」「0.4–19.1%」「CramérV 0.227」仍 literal(現值正確,drift 風險)。**R3 處置**:接受(值正確,low-risk;incpct 已動態);註記待未來 build 期注入。
- **[MINOR-7]** incompatible%↔branched% 重疊未揭露(H2009 c=2-branched 17.5% 同時四配子違反)。**R3 修**:§② 加揭露句(分母不互扣但 H2009 branched≠乾淨平行;COLO829/HCC1954 重疊<0.2% 不影響)。

### 假警報（R2 正確排除）
- 「incompatible 應從 branched% 分母扣掉」= 排除正確(兩軸正交、各區各計一次,扣掉反而扭曲)。
- 報告 line 136「topology_type=incompatible ≡ has_cycle」細微不精確(242/267 為 has_cycle=False)= ~1% 未用丟棄類,load-bearing 的 determinacy incompatible 5.9% 計算正確 → 附帶 doc 瑕疵,不獨立列。

### R3 處置（純傳播 sweep,無重算;grep + HTML playwright 驗證）
全部殘留同根因=display/provenance 層未傳播。R3 全數修正並驗:報告 forbidden 斷言消失(只剩「不可寫」否定用法)、CN 舊錨/filter 舊行號清除、HTML c≤k+1/48點 清除、HTML 面板 CN 40點顯示正確、0 pageerror。
