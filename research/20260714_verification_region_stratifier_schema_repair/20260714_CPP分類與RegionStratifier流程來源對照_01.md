<!--
建立時間: 2026-07-15 00:04 +08:00
更新時間: 2026-07-15 00:38 +08:00
狀態: validated-source-map
目標: 建立 Verification schema v2 與 RegionStratification schema v1 的 C++ input→判斷→狀態變更→序列化→fail-closed→測試證據完整來源對照
處理範圍:
  - InterSubMod/include/core/RegionProcessor.hpp
  - InterSubMod/include/core/SubcloneAnalyzer.hpp
  - InterSubMod/src/core/RegionProcessor.cpp
  - InterSubMod/src/core/SubcloneAnalyzer.cpp
  - InterSubMod/tests/test_verification_schema_v2.cpp
  - InterSubMod/tests/test_region_stratifier.cpp
  - InterSubMod/tests/test_region_stratification_artifacts.cpp
關聯檔案:
  - InterSubMod/research/20260714_verification_region_stratifier_schema_repair/pre-decision-audit.md
  - InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md
  - /big7_disk/liaoyoyo2001/external_validation/subclone_reconstruction_method_comparison_20260714/20260714_InterSubMod分類與RegionStratifier修復任務規格_01.md
證據層級:
  - L1: 目前 source/test 可直接定位的 code fact
  - L2: 由多個 L1 facts 串接出的工程推論，必須明標
  - RISK: 尚無直接測試或仍需 consumer/runtime 約束的未解風險
排除範圍:
  - 不把任何 RegionStratum 或 legacy Subclone 類名宣稱為 biological clone
  - 不提供或推論 F1 改善
  - 不涵蓋 Python migration、downstream consumers 或大型 BAM runtime
-->

# C++ 分類與 RegionStratifier 流程來源對照

> **結論先行（L1）**：目前 C++ 主線把單一 region 的既有四態 class 與已完成的 evidence axes 送入純 `Verification schema-v2` ordered tree，產生 11-class current taxonomy、四個 evidence booleans、`EvidencePath`、`Significant` 與 legacy-derived LOH；跨 region 階段則使用固定四槽的 deterministic region stratification。兩條流程都在序列化前檢查 schema、identity、alias 與 count consistency。
>
> **使用限制（L1/L2）**：`region_stratification_status.tsv` 是 publication commit marker。三個 owned artifacts 各自 atomic rename，但不是跨檔單一 transaction；consumer 必須先讀 terminal status，不能只看某一個資料檔。本文件只解釋工程分類／stratification，不建立 biological clone 或效能改善 claim。

## 0. 30 秒摘要

| 主題 | 已由目前 C++ 直接證明的事實 | 邊界 |
|---|---|---|
| 11 classes | `classify_verification_v2()` 依固定 first-match precedence 回傳 11 個 canonical classes；前七個 clean branches 的 `Significant=true`，`DispersionStructure` 與三個 `Noise_*` 為 false | 這是工程 taxonomy，不等於 biological state |
| 四個 evidence fields | `LabelFirstSupport`、`ClusterFirstSupport` 來自保存的 legacy class；`WithinHPSupport`、`DispersionWarning` 來自 live evidence axes | 不可由 final class 反向猜測原始 evidence |
| Significant invariance | legacy `Strong` 與 `Subclone` 分別拆成 `Strong_Bidirectional` 與 `ClusterFirstOnly`，兩者仍為 true；CSV `Significant` 直接序列化 decision boolean | C++ truth-table test 證明 mapping；frozen corpus 逐列 invariance 不在這三支 C++ tests 的證明範圍 |
| LOH provenance | `LOH_Subtype_LegacyVC` 只用 `potential_loh + VerificationClass_Legacy`；舊 `LOH_Subtype` 是逐列 exact alias | 不存在 `LOH_Subtype_Current` |
| Eligibility | 共用 predicate：`success && num_reads>=10 && num_cpgs>=3 && significance_computed` | 這裡檢查的是「significance 已計算」，不是 `verification_significant=true` |
| 0/49/50 | 0 先被 empty guard 擋下；預設 minimum 50 下，49 不執行 assignment、50 全部執行 | 0 即使 `min_regions=0` 仍不是 valid |
| 四 fixed slots | storage 永遠四槽；occupied count 另用 distinct occupied IDs 計算；`n_subclones` 只是 deprecated mirror | 四槽是 region pattern codes，不是四個 cellular clones |
| Fail-closed publication | run 開始先寫 `FAILED/RUN_IN_PROGRESS`；所有 artifacts/CSV durable 後才最後寫 terminal status；non-valid run 覆寫舊 valid rows；RAII guard 在 temp publication 失敗時移除本輪 owned temp | 跨檔一致性仍依賴 status-first consumer；各檔不是單一跨檔 transaction |

## 1. 證據分級、snapshot 與非主張

### 1.1 證據標籤

| 標籤 | 本文件中的意思 | 可以怎麼引用 |
|---|---|---|
| **L1 code fact** | 可由列出的目前 C++ source 或 test 行號直接讀出 | 可寫「目前實作為……」 |
| **L2 engineering inference** | 由至少兩個 L1 facts 合併而成，但不是單一 assertion 直接測得 | 必須寫「依目前控制流推論……」 |
| **RISK** | 缺專門測試、需要 consumer discipline，或 compatibility removal 尚未落地 | 不可寫成「已完全證明」 |

### 1.2 行號綁定的 source snapshot

本文件行號以 2026-07-15T00:37:31+08:00、Git HEAD `6067568637088838a9f518955e41d222f057e4f1` 的工作樹內容為準。因 source 尚在 working tree，後續若 hash 改變，應重新產生行號對照。

| 檔案 | SHA-256 |
|---|---|
| `InterSubMod/include/core/RegionProcessor.hpp` | `cfbc9a82329293e56517981129b4c6d13479759cb4165af7da0bfeffdbfe3c4c` |
| `InterSubMod/include/core/SubcloneAnalyzer.hpp` | `4e2d4361e26fc8ef8e1c11f3c903caac6122595efa3a2c3b91087694385c7cf0` |
| `InterSubMod/src/core/RegionProcessor.cpp` | `024f8fb22b7fb58b9e38bc0e9463c8dac5a373fb5b03a0729bf80dce7fbeb9ba` |
| `InterSubMod/src/core/SubcloneAnalyzer.cpp` | `bd0bb14cc867cc079b7b0f753bd54c0078f8afd987aaa09220cabe7b857649b2` |
| `InterSubMod/tests/test_verification_schema_v2.cpp` | `72ba9683975ca5d93e7adfb34ab987bcaab55d5e0ee5d52b9d0fed1b93f17f3b` |
| `InterSubMod/tests/test_region_stratifier.cpp` | `ff0c4218e913863599c59e53ec14b562eebea66edc5029aa29bad3eddf8493bf` |
| `InterSubMod/tests/test_region_stratification_artifacts.cpp` | `aeeec85b74a4f28189f1c70b26bbf9f991064ac667d218b70ccea3a2a7b84abf` |

### 1.3 明確非主張

1. `RegionStratum::{BaselineLowAsm,HighHpAsm,LohFlagged,HighSampleAsm}` 是 deterministic region-level codes；header 已明文禁止解讀成 cellular subclones：`InterSubMod/include/core/SubcloneAnalyzer.hpp:18-29`。
2. `SubcloneAnalyzer`、`SubcloneSummary`、`n_subclones`、`Subclone_ID` 舊名只是 migration compatibility surface，不構成 biological interpretation。
3. 本文件不評估 performance metric，也不把「測試全綠」外推為生物學有效性。

## 2. 端到端控制流

### 2.1 單一 region：legacy result 到 Verification schema v2

```text
既有 clustering / label / Δβ / PERMANOVA / dispersion / epipoly evidence
    ↓
保存 Stage④ 前的 VerificationClass → VerificationClass_Legacy
    ↓
組 VerificationV2Input
    ├─ legacy class
    ├─ dbeta_sig / hp_auc_struct / potential_loh / within_hp
    ├─ clean_location_permanova / dispersion_structure / dispersion_warning
    └─ epipoly / pairwise mean
    ↓
classify_verification_v2()：first-match wins
    ↓
寫回 current class、v1 deprecated projection、四 evidence booleans、
EvidencePath、EvidenceDerivation、Significant
    ↓
determine_loh_subtype_legacy_vc(potential_loh, legacy class)
    ↓
LOH_Subtype_LegacyVC + deprecated exact alias LOH_Subtype
    ↓
write_significance_summary() 先 validate schema/aliases/sentinels，再 atomic publish CSV
```

主要來源：`InterSubMod/src/core/RegionProcessor.cpp:1794-1842`、`InterSubMod/src/core/RegionProcessor.cpp:1856-2205`。

### 2.2 全 run：status-first 到 terminal status-last

```text
atomic status = FAILED / RUN_IN_PROGRESS
    ↓
process all regions
    ↓
reset all RegionStratum fields to unassigned sentinels
    ↓
tumor-only? ─ yes → NOT_APPLICABLE_TUMOR_ONLY
    │
    no
    ↓
SubcloneAnalyzer compatibility façade
    → shared eligibility
    → 0/minimum guard
    → deterministic 4-way assignment
    → fixed-slot summaries + occupied distinct count
    ↓
validate ALL assignments
    ↓
commit ALL RegionResult mutations, or zero partial write-back
    ↓
construct stable-key artifact rows
    ↓
atomic owned artifacts
    ↓
atomic significance CSV/statistics
    ↓
atomic terminal status LAST
```

主要來源：`InterSubMod/src/core/RegionProcessor.cpp:911-1217`。

## 3. Verification schema v2 source map

### 3.1 逐 symbol ledger

| Symbol／來源 | Input | Predicate／precedence | State mutation | Serialized output | Fail-closed reason | Test evidence |
|---|---|---|---|---|---|---|
| `VerificationV2Input`／`VerificationV2Decision`，`InterSubMod/include/core/RegionProcessor.hpp:31-65` | legacy class、7 個結構/evidence flags、epipoly、pairwise mean | struct 本身不判斷；`schema_version` default 2、`evidence_derivation` default `LIVE`、`significant` default false | 無；是 pure classifier 的 value objects | decision 稍後映射到 `RegionResult` 與 CSV | 不適用；值域由 classifier 驗證 | `InterSubMod/tests/test_verification_schema_v2.cpp:13-25` 建 fixture |
| `is_valid_verification_legacy_class()`，`InterSubMod/src/core/RegionProcessor.cpp:44-46` | 任意 legacy 字串 | exact enum：`Strong`、`Subclone`、`Weak`、`Noise` | 無 | 不直接序列化 | false 會使 classifier／LOH helper throw `INVALID_VERIFICATION_CLASS_LEGACY` | `InterSubMod/tests/test_verification_schema_v2.cpp:125-139` |
| `classify_verification_v2()`，`InterSubMod/src/core/RegionProcessor.cpp:50-112` | `VerificationV2Input` | 先驗 legacy；再依 1a/1b→LOH→within-HP→Δβ→clean PERMANOVA→HP-AUC→dispersion→uniform→chaotic→uncorrelated first-match | 只組新 `VerificationV2Decision`；不改 caller object | current class、v1 projection、legacy、四 evidence booleans、path、derivation、Significant | invalid legacy 直接 throw；無 silent Noise fallback | 11-row truth table：`InterSubMod/tests/test_verification_schema_v2.cpp:30-87`；collision：`InterSubMod/tests/test_verification_schema_v2.cpp:89-123`；unknown：`InterSubMod/tests/test_verification_schema_v2.cpp:135-139` |
| Stage④ integration block，`InterSubMod/src/core/RegionProcessor.cpp:1794-1842` | `RegionResult` 中已完成的 legacy class 與 live axes | 聚合 Δβ flags；`hp_auc_all>=0.7`；PERMANOVA 必須 valid 且 `p<0.05`；clean 與 dispersion 依 warning 分流 | 先保存 `verification_class_legacy`，再寫回 schema/current/v1 projection/4 booleans/path/derivation/significant/LOH aliases | `significance_summary.csv` 的 current、legacy 與 tail provenance | classifier/helper throw 被 `process_single_region()` catch，該 region `success=false`、`error_message=e.what()`：`InterSubMod/src/core/RegionProcessor.cpp:1844-1848` | pure branch test 間接覆蓋 decision；沒有直接 fixture 驗證此 integration block 的每一個 `RegionResult` source field |
| `determine_loh_subtype_legacy_vc()`，`InterSubMod/src/core/RegionProcessor.cpp:115-126` | `potential_loh`、保存的 legacy class | 先驗 legacy；false→`None`；true 依四態→`LOH_Noise/Weak/Strong/Subclone` | pure return；Stage④ 寫入 canonical 與 alias | `LOH_Subtype_LegacyVC`；舊 `LOH_Subtype` | unknown legacy 即使 `potential_loh=false` 也 throw，不會降成 `None` | `InterSubMod/tests/test_verification_schema_v2.cpp:125-133` |
| Verification/LOH `RegionResult` fields，`InterSubMod/include/core/RegionProcessor.hpp:323-334`、`InterSubMod/include/core/RegionProcessor.hpp:369-370` | classifier/LOH decision | 欄位只承接 authoritative decision | 保存 current、v1 deprecated、legacy、4 booleans、path/derivation/significant、LOH canonical/alias | CSV current/legacy 在 `InterSubMod/src/core/RegionProcessor.cpp:2075-2079`；tail provenance 在 `InterSubMod/src/core/RegionProcessor.cpp:2181-2189` | CSV publication 前再驗 schema 與 aliases：`InterSubMod/src/core/RegionProcessor.cpp:1958-1995` | truth table 驗 decision；CSV row guard 沒有獨立 C++ test |
| `write_significance_summary()`，`InterSubMod/src/core/RegionProcessor.cpp:1856-2254` | 全部 `RegionResult` 與 `snvs_` | 只序列化 `success && significance_computed` rows；逐列驗 region ID、schema、LOH alias、RegionStratum alias/enum/sentinel | 不改 classification；累積統計；RAII-owned temp→flush→close→rename→mark published | `significance_summary.csv`、`significance_statistics.txt` | schema/alias/sentinel 與 I/O errors；未 mark published 的 owned temp 由 `TempFileGuard` 移除 | generic rename-failure cleanup：`InterSubMod/tests/test_region_stratification_artifacts.cpp:184-200`；private CSV/statistics writer仍無直接 failure-injection fixture |

### 3.2 Authoritative 11-class ordered decision table

以下是 `classify_verification_v2()` 的實際 first-match order；條件相撞時只採第一個命中 branch。來源：`InterSubMod/src/core/RegionProcessor.cpp:62-105`。

| 順序 | 實際 predicate | `VerificationClass` | `EvidencePath` | `Significant` |
|---:|---|---|---|---:|
| 1a | cluster-first support 且 legacy `Strong` | `Strong_Bidirectional` | `BIDIRECTIONAL` | true |
| 1b | cluster-first support 且 legacy `Subclone` | `ClusterFirstOnly` | `CLUSTER_FIRST_ONLY` | true |
| 2 | `potential_loh && (hp_auc_struct \|\| dbeta_sig)` | `LOH-Structure` | `LOH_STRUCTURE` | true |
| 3 | `within_hp` | `MultiGroupNoLabel` | `WITHIN_HP_MULTIGROUP` | true |
| 4 | `dbeta_sig` | `LabelShift` | `LABEL_SHIFT` | true |
| 5 | `clean_location_permanova` | `PermanovaLocation` | `PERMANOVA_LOCATION` | true |
| 6 | `hp_auc_struct` | `StructureNoLabel` | `HP_AUC_STRUCTURE_NO_LABEL` | true |
| 7 | `dispersion_structure` | `DispersionStructure` | `DISPERSION_STRUCTURE` | false |
| 8 | `pairwise_mean_distance < 0.15` 或 non-NaN `epipoly_hp1 < 0.2` | `Noise_Uniform` | `NOISE_UNIFORM` | false |
| 9 | non-NaN `epipoly_hp1 > 0.5` 或 `pairwise_mean_distance > 0.35` | `Noise_Chaotic` | `NOISE_CHAOTIC` | false |
| 10 | 以上皆未命中 | `Noise_Uncorrelated` | `NOISE_UNCORRELATED` | false |

精確邊界：

- `0.15` 與 `0.2` 本身不命中 uniform，因為使用 strict `<`。
- `0.5` 與 `0.35` 本身不命中 chaotic，因為使用 strict `>`。
- `epipoly_hp1` 是 NaN 時，兩個 epipoly predicate 都被明確跳過。
- precedence collisions 已由 `OrderedPrecedenceIsStableUnderCollisions` 驗證 cluster > LOH > within-HP > Δβ > clean PERMANOVA > HP-AUC：`InterSubMod/tests/test_verification_schema_v2.cpp:89-123`。
- collision test 尚未以單一 fixture 專門驗證「dispersion 不蓋過每一個 clean branch」；source order保證其位於所有 clean branches 之後，truth table 另驗 dispersion 本身為 false。

### 3.3 四個 evidence fields 的來源、寫回與輸出

| 欄位 | Live authoritative source | Decision assignment | `RegionResult` | CSV | Test |
|---|---|---|---|---|---|
| `LabelFirstSupport` | `legacy_class in {Strong,Weak}` | `InterSubMod/src/core/RegionProcessor.cpp:56` | `verification_label_first_support`，`InterSubMod/include/core/RegionProcessor.hpp:328` | header `InterSubMod/src/core/RegionProcessor.cpp:1941-1942`；value `InterSubMod/src/core/RegionProcessor.cpp:2182` | `InterSubMod/tests/test_verification_schema_v2.cpp:74-75` |
| `ClusterFirstSupport` | `legacy_class in {Strong,Subclone}` | `InterSubMod/src/core/RegionProcessor.cpp:57-58` | `verification_cluster_first_support`，`InterSubMod/include/core/RegionProcessor.hpp:329` | header `InterSubMod/src/core/RegionProcessor.cpp:1942`；value `InterSubMod/src/core/RegionProcessor.cpp:2183` | `InterSubMod/tests/test_verification_schema_v2.cpp:76-77` |
| `WithinHPSupport` | `result.within_hp_clean_multigroup` | integration `InterSubMod/src/core/RegionProcessor.cpp:1802`、`InterSubMod/src/core/RegionProcessor.cpp:1819`；decision copy `InterSubMod/src/core/RegionProcessor.cpp:59` | `verification_within_hp_support`，`InterSubMod/include/core/RegionProcessor.hpp:330` | header `InterSubMod/src/core/RegionProcessor.cpp:1942`；value `InterSubMod/src/core/RegionProcessor.cpp:2184` | `InterSubMod/tests/test_verification_schema_v2.cpp:78` |
| `DispersionWarning` | HP/allele dispersion warning 的 OR | integration `InterSubMod/src/core/RegionProcessor.cpp:1822-1823`；decision copy `InterSubMod/src/core/RegionProcessor.cpp:60` | `verification_dispersion_warning`，`InterSubMod/include/core/RegionProcessor.hpp:331` | header `InterSubMod/src/core/RegionProcessor.cpp:1942`；value `InterSubMod/src/core/RegionProcessor.cpp:2185` | `InterSubMod/tests/test_verification_schema_v2.cpp:79` |

另有兩個 provenance fields：

- `EvidencePath` 是實際 first-match branch，於 `InterSubMod/src/core/RegionProcessor.cpp:63-105` 設值。
- `EvidenceDerivation` live default 固定為 `LIVE`，定義於 `InterSubMod/include/core/RegionProcessor.hpp:59`，測試於 `InterSubMod/tests/test_verification_schema_v2.cpp:80`。

### 3.4 Significant invariance 的可證明範圍

**L1 code facts**

1. legacy `Strong` 命中 `Strong_Bidirectional` 且 `significant=true`；legacy `Subclone` 命中 `ClusterFirstOnly` 且 `significant=true`：`InterSubMod/src/core/RegionProcessor.cpp:63-70`。
2. 11-row test 明確 assert 每一 class 的 expected boolean：`InterSubMod/tests/test_verification_schema_v2.cpp:30-72`。
3. Stage④ 將 decision boolean 原樣寫到 `result.verification_significant`：`InterSubMod/src/core/RegionProcessor.cpp:1828-1837`。
4. CSV `Significant` 不再重新計算另一套 gate；直接取 `r.verification_significant`：`InterSubMod/src/core/RegionProcessor.cpp:2001-2007`，並在 `InterSubMod/src/core/RegionProcessor.cpp:2164-2169` 序列化。
5. `SuggestFilter` 仍是獨立舊欄位 `label_delta > 0.3`，不參與 `Significant`：`InterSubMod/src/core/RegionProcessor.cpp:2009-2010`。

**L2 inference**

因 split 前後的兩個 legacy positive branches 均維持 true，僅把它們拆成兩個 current labels，本次 class rename/split 本身不會改變這兩群的 binary `Significant`。這個推論由 pure truth table 支撐。

**證明邊界／RISK**

三支新 C++ test files（共 17 tests）沒有讀 frozen 14-file corpus，也沒有逐 stable key 比較修復前後 serialized `Significant` 與 raw metrics；因此「全 corpus 每列 invariance」必須引用獨立 golden migration evidence，不能只引用這份 C++ source map。

### 3.5 LOH legacy derivation 與 exact alias

| 步驟 | L1 fact | 來源 |
|---|---|---|
| 保存 provenance | Stage④ 在覆寫 current class 前，先把舊值放入 `verification_class_legacy` | `InterSubMod/src/core/RegionProcessor.cpp:1794-1798` |
| canonical derivation | helper 只接受四個 legacy enum；`potential_loh=false`→`None`，true 再依 legacy class 選 subtype | `InterSubMod/src/core/RegionProcessor.cpp:115-126` |
| state mutation | canonical 寫 `loh_subtype_legacy_vc`；deprecated `loh_subtype` 立即複製同值 | `InterSubMod/src/core/RegionProcessor.cpp:1839-1841` |
| schema declaration | 兩欄的 canonical/deprecated 意義寫在 `RegionResult` | `InterSubMod/include/core/RegionProcessor.hpp:369-370` |
| publication guard | 兩欄不相等即 `LOH_SUBTYPE_ALIAS_FAILED`，CSV 不發布 | `InterSubMod/src/core/RegionProcessor.cpp:1966-1968` |
| serialization | 舊 `LOH_Subtype` 位於既有欄序；canonical `LOH_Subtype_LegacyVC` append 在 tail | header `InterSubMod/src/core/RegionProcessor.cpp:1895-1897`、`InterSubMod/src/core/RegionProcessor.cpp:1940-1944`；values `InterSubMod/src/core/RegionProcessor.cpp:2091-2097`、`InterSubMod/src/core/RegionProcessor.cpp:2181-2189` |
| unit evidence | 五個合法 output 與 unknown true/false 兩種情境均驗證 | `InterSubMod/tests/test_verification_schema_v2.cpp:125-133` |

## 4. RegionStratifier source map

### 4.1 型別與固定契約

| Symbol／來源 | Input | Predicate／precedence | State mutation | Serialized output | Fail-closed reason | Test evidence |
|---|---|---|---|---|---|---|
| constants + `RegionStratum` + status enum，`InterSubMod/include/core/SubcloneAnalyzer.hpp:14-40` | 無 runtime input | schema=1、slot count=4、default minimum=50；IDs 固定 0..3；status 固定四態 | 無 | labels/reasons/status strings | invalid enum label helper 回 `Unknown/INVALID_STRATUM`；status fallback `FAILED` | IDs/slots/status 由全部 region/artifact tests 使用 |
| `RegionAsmProfile`，`InterSubMod/include/core/SubcloneAnalyzer.hpp:45-70` | eligible `RegionResult` 的 ASM/LOH/current/legacy/schema/QC | 無內部判斷 | 是 analyzer 的 immutable projection | 不直接作 canonical artifact row | current 版本在 summary counting 階段檢查 | index-preservation `InterSubMod/tests/test_region_stratifier.cpp:120-145` |
| `RegionStratumAssignment`，`InterSubMod/include/core/SubcloneAnalyzer.hpp:72-76` | original result index、region ID、typed stratum | typed stratum | 先暫存；commit 前不改 `RegionResult` | 之後轉 artifact row | complete-set validator 驗 identity/count | `InterSubMod/tests/test_region_stratifier.cpp:147-180` |
| `SubcloneSummary` compatibility struct，`InterSubMod/include/core/SubcloneAnalyzer.hpp:83-114` | 每一 fixed slot 的 profiles | current 11 leaves + unknown；legacy 4 leaves + unknown | 累積 count/means、dominant ties | canonical human summary只列 occupied slots | unknown 不 silent fallback；轉 warning count | `InterSubMod/tests/test_region_stratifier.cpp:182-240` |
| `SubcloneResult` compatibility struct，`InterSubMod/include/core/SubcloneAnalyzer.hpp:116-136` | analyzer 整體結果 | valid/status/reason + assigned/occupied counts | 保存四槽、assignments 與 deprecated mirrors | process lifecycle 轉成 status/artifacts | default `FAILED/RUN_IN_PROGRESS` | 0/49/50 與 occupied tests |

### 4.2 Eligibility：唯一共用 predicate

`is_region_stratification_eligible()` 的唯一條件位於 `InterSubMod/src/core/SubcloneAnalyzer.cpp:97-100`：

```cpp
result.success &&
result.num_reads >= 10 &&
result.num_cpgs >= 3 &&
result.significance_computed
```

| 邊界 | 結果 | Test evidence |
|---|---:|---|
| reads=9 | ineligible | `InterSubMod/tests/test_region_stratifier.cpp:120-135` |
| reads=10、CpGs=3、computed=true | eligible | `InterSubMod/tests/test_region_stratifier.cpp:16-28`、`InterSubMod/tests/test_region_stratifier.cpp:124`、`InterSubMod/tests/test_region_stratifier.cpp:127`、`InterSubMod/tests/test_region_stratifier.cpp:131-135` |
| CpGs=2 | ineligible | `InterSubMod/tests/test_region_stratifier.cpp:125-126`、`InterSubMod/tests/test_region_stratifier.cpp:133` |
| CpGs=3 | eligible | fixture default `InterSubMod/tests/test_region_stratifier.cpp:20-22` |
| `significance_computed=false` | ineligible | `InterSubMod/tests/test_region_stratifier.cpp:128-135` |

**L1 distinction**：eligibility 不要求 `verification_significant=true`。它只要求統計流程已執行，讓非-significant/noise regions 也能成為 cross-region stratification 的輸入。

### 4.3 RegionStratum ordered predicate

`assign_region_stratum()` 位於 `InterSubMod/src/core/SubcloneAnalyzer.cpp:102-110`；門檻常數位於 `InterSubMod/src/core/SubcloneAnalyzer.cpp:18-21`。

| 順序 | Predicate | ID | Label | Reason |
|---:|---|---:|---|---|
| 1 | `potential_loh \|\| loh_bed_overlap` | 2 | `LohFlagged` | `LOH_FLAGGED` |
| 2 | `abs(sample_asm_delta) > 0.10` | 3 | `HighSampleAsm` | `HIGH_SAMPLE_ASM` |
| 3 | `abs(hp_asm_delta) > 0.05` | 1 | `HighHpAsm` | `HIGH_HP_ASM` |
| 4 | otherwise | 0 | `BaselineLowAsm` | `BASELINE_LOW_ASM` |

精確邊界與碰撞 test：

- `0.05`、`0.10` 本身仍是 baseline，因為 strict `>`：`InterSubMod/tests/test_region_stratifier.cpp:106-110`。
- `0.050001`→HighHpAsm、`0.100001`→HighSampleAsm；再加 LOH→LohFlagged：`InterSubMod/tests/test_region_stratifier.cpp:112-117`。
- `region_stratum_label/reason()` 使用 fixed arrays，invalid typed value不會被當成某個合法 slot：`InterSubMod/src/core/SubcloneAnalyzer.cpp:23-26`、`InterSubMod/src/core/SubcloneAnalyzer.cpp:66-80`。

### 4.4 逐 symbol：projection、assignment、validation、commit

| Symbol／來源 | Input | Predicate／precedence | State mutation | Serialized output | Fail-closed reason | Test evidence |
|---|---|---|---|---|---|---|
| `extract_profiles()`，`InterSubMod/src/core/SubcloneAnalyzer.cpp:226-257` | full `results` vector | 只呼叫共用 eligibility | 建 profile 並保存原 `result_index`、`region_id`、current/legacy/schema | 間接進 summaries/artifact rows | ineligible row完全不投影 | `SharedEligibilityPreservesOriginalResultIndices`，`InterSubMod/tests/test_region_stratifier.cpp:120-145` |
| `assign_strata()`，`InterSubMod/src/core/SubcloneAnalyzer.cpp:259-268` | profiles，順序由 extraction 保存 | 每 profile 呼叫 ordered predicate | 建 `{result_index,region_id,stratum}`，不改 original results | 間接 | 無；完整一致性留給 validator | precedence test `InterSubMod/tests/test_region_stratifier.cpp:106-118` |
| `validate_region_stratification_assignments()`，`InterSubMod/src/core/SubcloneAnalyzer.cpp:112-147` | original results、完整 assignments、expected count | 先 count，再重新數 eligibility，再驗 index unique/range、eligible、region ID、typed stratum、無漏 assign | **零 mutation** | 無 | exact reasons：`ASSIGNMENT_COUNT_MISMATCH`、`ELIGIBLE_REGION_COUNT_MISMATCH`、`RESULT_INDEX_OUT_OF_RANGE`、`DUPLICATE_RESULT_INDEX`、`INELIGIBLE_ASSIGNMENT`、`REGION_ID_MISMATCH`、`INVALID_STRATUM`、`ELIGIBLE_REGION_UNASSIGNED` | duplicate/mismatch/out-of-range：`InterSubMod/tests/test_region_stratifier.cpp:147-170`；其餘 reason 尚無逐一專測 |
| `commit_region_stratification_assignments()`，`InterSubMod/src/core/SubcloneAnalyzer.cpp:149-168` | 同 validator | complete validation 必須先成功 | 才一次寫 schema、ID、label、reason、`subclone_id` alias | 後續 CSV/artifacts | validation 任一 false 即原樣 return，零 partial write-back | invalid cases保持 ID=-1，valid才全部寫回：`InterSubMod/tests/test_region_stratifier.cpp:147-180` |

### 4.5 0／49／50 與 fixed-slot state machine

`SubcloneAnalyzer::analyze()` 位於 `InterSubMod/src/core/SubcloneAnalyzer.cpp:170-224`。

| Eligible count | Guard／status | Assignment count | Occupied count | Summary storage | Test |
|---:|---|---:|---:|---:|---|
| 0 | `INSUFFICIENT_REGIONS / NO_ELIGIBLE_REGIONS`，在 minimum 比較前返回 | 0 | 0 | 4 | `InterSubMod/tests/test_region_stratifier.cpp:52-61` |
| 49（default minimum=50） | `INSUFFICIENT_REGIONS / BELOW_MIN_REGIONS` | 0 | 0 | 4 | `InterSubMod/tests/test_region_stratifier.cpp:63-70` |
| 50 | `VALID / OK`（無 unknown warning 時） | 50 | 依 distinct occupied IDs；fixture 為 1 | 4 | `InterSubMod/tests/test_region_stratifier.cpp:71-77` |

Fixed-slot／count facts：

1. 每次 analyze 一開始就以 `compute_summaries({}, {})` 建四槽，即使 early return 也保留四槽：`InterSubMod/src/core/SubcloneAnalyzer.cpp:170-181`。
2. occupied count 用 `std::array<bool,4>` 對實際 IDs 去重後計算：`InterSubMod/src/core/SubcloneAnalyzer.cpp:189-200`。
3. `n_subclones = n_occupied_region_strata` 只是 deprecated exact mirror：`InterSubMod/src/core/SubcloneAnalyzer.cpp:197-200`。
4. `compute_summaries()` 固定配置 4 筆並令 slot index 等於 compatibility `subclone_id`：`InterSubMod/src/core/SubcloneAnalyzer.cpp:270-275`。
5. artifact writer 要求 summary size 恰為 4，否則拒絕發布：`InterSubMod/src/core/RegionProcessor.cpp:275-278`。
6. only ID 3、IDs 0+3、only ID 2 均證明 storage 長度與 occupied count 已分離：`InterSubMod/tests/test_region_stratifier.cpp:79-104`。

### 4.6 Current／legacy leaf counters、unknown 與 ties

`compute_summaries()` 的 class counting 位於 `InterSubMod/src/core/SubcloneAnalyzer.cpp:297-335`。

| Taxonomy | Accepted leaves | Version gate | Unknown policy | Dominant tie |
|---|---|---|---|---|
| current v2 | 11 classes，固定 array order：`InterSubMod/src/core/SubcloneAnalyzer.cpp:28-31` | 只有 `verification_schema_version==2` 才查 current enum：`InterSubMod/src/core/SubcloneAnalyzer.cpp:297-305` | 進 `verification_class_v2_unknown`，不轉 Noise | 全部 max-count ties 依 array order以 `;` 串接：`InterSubMod/src/core/SubcloneAnalyzer.cpp:43-64`、`InterSubMod/src/core/SubcloneAnalyzer.cpp:329-331` |
| legacy | `Strong/Subclone/Weak/Noise`：`InterSubMod/src/core/SubcloneAnalyzer.cpp:32-33` | 無 current schema version gate；直接讀保存 legacy 欄 | 進 `verification_class_legacy_unknown`，不轉 Noise | 同上：`InterSubMod/src/core/SubcloneAnalyzer.cpp:332-335` |

守恆／warning：

- 11 個 current leaves各 1，加 1 unknown，current leaf sum=12；legacy leaf+unknown同樣=12：`InterSubMod/tests/test_region_stratifier.cpp:182-215`。
- current unknown 與 legacy unknown 各加一個 warning，reason 成為 `VALID_WITH_WARNINGS`：`InterSubMod/src/core/SubcloneAnalyzer.cpp:216-223`；test `InterSubMod/tests/test_region_stratifier.cpp:216-217`。
- v2 class 配 schema version 1 時不偷認為 canonical，明確落入 `UnknownCurrentClass` 並 warning：`InterSubMod/tests/test_region_stratifier.cpp:225-240`。

## 5. 序列化、stable key、status 與 stale overwrite

### 5.1 逐 symbol publication ledger

| Symbol／來源 | Input | Predicate／precedence | State mutation | Serialized output | Fail-closed reason | Test evidence |
|---|---|---|---|---|---|---|
| artifact/status records，`InterSubMod/include/core/RegionProcessor.hpp:67-102` | stable identity、stratum enum fields、status counts | schema defaults、status default `FAILED/RUN_IN_PROGRESS` | value records | status/assignment/summary/stub | writer再驗所有 invariants | 五個 artifact tests |
| `validate_region_status_record()`，`InterSubMod/src/core/RegionProcessor.cpp:149-184` | status record | schema、nonempty reason、nonnegative counts、occupied bounds；VALID/non-VALID consistency | 無 | 通過才允許任何 status/artifact output | `INVALID_REGION_STATUS_SCHEMA_VERSION`、`EMPTY_REGION_STATUS_REASON`、`NEGATIVE_REGION_STATUS_COUNT`、`INVALID_REGION_STATUS_OCCUPIED_COUNT`、`INCONSISTENT_VALID_REGION_STATUS`、`NONVALID_REGION_STATUS_HAS_ASSIGNMENTS` | contradictory status：`InterSubMod/tests/test_region_stratification_artifacts.cpp:160-181` |
| `TempFileGuard`／`atomic_temp_path()`，`InterSubMod/src/core/RegionProcessor.cpp:131-136`、`InterSubMod/src/core/RegionProcessor.cpp:186-198` | 本輪唯一 temp path | 未 `mark_published()` 即 destructor `std::remove()`；成功 rename 後才標 published | 只清除自己持有的 temp，不動 final path | generic artifacts、CSV、statistics 共用 cleanup policy | callback、flush、close、rename 或 scope exception 都不遺留 owned temp | forced rename failure：`InterSubMod/tests/test_region_stratification_artifacts.cpp:184-200` |
| `write_atomic_text()`，`InterSubMod/src/core/RegionProcessor.cpp:200-228` | final path + writer callback | 同目錄 temp，成功 flush/close 後 rename；guard 只在 rename 成功後 mark published | 只改 filesystem；成功時 final path 原子替換，失敗時移除 owned temp | 所有 canonical text artifacts | `CREATE_OUTPUT_DIRECTORY_FAILED`、`OPEN_TEMP_FAILED`、writer error、`FLUSH_TEMP_FAILED`、`CLOSE_TEMP_FAILED`、`ATOMIC_RENAME_FAILED` | success 後無 `.tmp.`：`InterSubMod/tests/test_region_stratification_artifacts.cpp:106-108`；rename failure也無 `.tmp.`：`InterSubMod/tests/test_region_stratification_artifacts.cpp:184-200` |
| `rfc3339_utc_now()`，`InterSubMod/src/core/RegionProcessor.cpp:232-239` | system UTC time | fixed `YYYY-MM-DDTHH:MM:SSZ` | 無 | status `generated_at` | 無顯式 clock failure branch | regex：`InterSubMod/tests/test_region_stratification_artifacts.cpp:39-54` |
| `write_region_stratification_status_atomic()`，`InterSubMod/src/core/RegionProcessor.cpp:241-258` | status record | 先完整 status validation；空 timestamp才現算 | atomic replace status marker | exact 9-column header + one value row | status validation/I/O failure回 false | exact header/defaults/timestamp：`InterSubMod/tests/test_region_stratification_artifacts.cpp:39-54` |
| `write_region_stratification_artifacts_atomic()`，`InterSubMod/src/core/RegionProcessor.cpp:260-391` | rows、恰四 summaries、status | 先驗 status、row count、row schema/enum/identity、stable key unique、summary slot/count、occupied count | 依序 atomic replace三個 owned artifacts | canonical TSV、summary、deprecated stub | `REGION_ARTIFACT_*` family；任一 validation fail在建立 output dir 前 return | valid：`InterSubMod/tests/test_region_stratification_artifacts.cpp:56-109`；stale overwrite：`InterSubMod/tests/test_region_stratification_artifacts.cpp:111-158`；contradiction：`InterSubMod/tests/test_region_stratification_artifacts.cpp:160-181`；rename cleanup：`InterSubMod/tests/test_region_stratification_artifacts.cpp:184-200` |
| `process_all_regions()` lifecycle，`InterSubMod/src/core/RegionProcessor.cpp:911-1217` | run inputs、all results、SNV table、output dir | status-first；paired/tumor-only；analyze→validate→commit；artifact→CSV→status-last | reset/commit RegionResult；filesystem publication | 全部 canonical artifacts + status marker | status/artifact/CSV failure皆 throw；artifact/CSV failure先嘗試 failed/empty overwrite | 三支 test files驗 pure/helper層；尚無完整 process lifecycle fault-injection test |

### 5.2 Canonical owned outputs

| 路徑 | VALID 內容 | non-VALID 內容 | Source |
|---|---|---|---|
| `region_stratification_status.tsv` | exact header + terminal status/counts/timestamp | exact header + non-VALID、assignment=0、occupied=0 | `InterSubMod/src/core/RegionProcessor.cpp:241-258` |
| `region_stratification.tsv` | header + 每個 assigned region 一列 | 只有 header、零 rows | `InterSubMod/src/core/RegionProcessor.cpp:325-340` |
| `region_stratification_summary.txt` | 標題/status/reason；只列 occupied slots及 current/legacy leaf counts | 只列標題/status/reason | `InterSubMod/src/core/RegionProcessor.cpp:342-381` |
| `subclone_structure.txt` | deprecation stub，指向 canonical summary | 同格式 deprecation stub，換本次 status/reason | `InterSubMod/src/core/RegionProcessor.cpp:383-390` |
| `significance_summary.csv` | 既有欄序 + tail-appended v2/RegionStratification provenance | 每個成功且 computed row仍含 `-1/Unassigned/<phase reason>` sentinel；失敗移除 owned temp | header `InterSubMod/src/core/RegionProcessor.cpp:1871-1944`；row `InterSubMod/src/core/RegionProcessor.cpp:1955-2189`；publish `InterSubMod/src/core/RegionProcessor.cpp:2191-2205` |
| `significance_statistics.txt` | 使用 `verification_significant` 統計，說明為 schema-v2 clean evidence projection；失敗移除 owned temp | 若寫入失敗，pipeline轉 terminal FAILED | `InterSubMod/src/core/RegionProcessor.cpp:2208-2254` |

### 5.3 Stable key 與 authoritative write-back identity

**L1 code facts**

1. assignment write-back authoritative identity 是 original `result_index`；`region_id` 只作一致性驗證：`InterSubMod/src/core/SubcloneAnalyzer.cpp:128-145`、`InterSubMod/src/core/SubcloneAnalyzer.cpp:157-166`。
2. artifact serialized stable key 是 `(RegionID,Chr,Pos,Ref,Alt)`，內部以 unit separator `\x1f` 連接；duplicate 立即 `DUPLICATE_REGION_ARTIFACT_STABLE_KEY`：`InterSubMod/src/core/RegionProcessor.cpp:280-306`。
3. process 建 row 時先驗 `result.region_id` 可索引 `snvs_`，再取 Chr/Pos/Ref/Alt；失敗就清空 rows、歸零 counts、把所有 state退回 failed sentinel：`InterSubMod/src/core/RegionProcessor.cpp:1128-1161`。
4. valid artifact test assert實際 serialization 含 `42 chr17 7579472 G A` stable identity：`InterSubMod/tests/test_region_stratification_artifacts.cpp:56-92`。

**RISK**

- 現有 artifact test沒有建兩列相同 stable key 以直接 assert `DUPLICATE_REGION_ARTIFACT_STABLE_KEY`；目前只有 source guard。
- `result_index` unique 與 stable key unique 是兩個不同 invariants，不能以其中一個取代另一個。

### 5.4 Status-first／status-last 與 stale overwrite

| 時序 | L1 行為 | Source |
|---:|---|---|
| 1 | 任何 region work 前，atomic 寫 `FAILED/RUN_IN_PROGRESS`；寫不成直接 throw | `InterSubMod/src/core/RegionProcessor.cpp:911-921` |
| 2 | 所有 results 先設 schema=1、ID=-1、label=Unassigned；依 terminal phase設定 reason | `InterSubMod/src/core/RegionProcessor.cpp:1060-1126` |
| 3 | assignment complete validation後才 commit；失敗則清空 assignment/count，所有列回 `-1/Unassigned/FAILED` | `InterSubMod/src/core/RegionProcessor.cpp:1087-1107` |
| 4 | owned artifacts先發布；若主寫入失敗，嘗試以 FAILED/empty representation覆寫全部 owned artifacts，再覆寫 failed status | `InterSubMod/src/core/RegionProcessor.cpp:1163-1192` |
| 5 | significance CSV/statistics再發布；失敗時同樣嘗試 failed artifacts/status | `InterSubMod/src/core/RegionProcessor.cpp:1194-1206` |
| 6 | terminal status 在全部 canonical artifacts/CSV durable 後最後發布 | `InterSubMod/src/core/RegionProcessor.cpp:1208-1213` |

Stale overwrite test先寫一個 valid `chr1` row，再依序以 `INSUFFICIENT_REGIONS`、`NOT_APPLICABLE_TUMOR_ONLY`、`FAILED` 重寫；每次都驗 assignment只剩 header、summary無舊 slot、stub帶本次 reason：`InterSubMod/tests/test_region_stratification_artifacts.cpp:111-158`。

**L2 operational inference**

若 process在 artifacts寫到一半時 crash，status marker仍是 `RUN_IN_PROGRESS`；遵守 status-first的 consumer會 fail closed，不會把混合世代 artifacts視為 terminal valid。這是由 status-first與status-last控制流推得，不是 crash-injection test。

## 6. Compatibility aliases 與移除邊界

| Deprecated surface | Canonical source／exact relation | Enforcement／serialization | Source |
|---|---|---|---|
| `VerificationClass_V1_Deprecated` | 兩個新 positive classes都投影成 `Strong`；其餘等於 current class | classifier直接產生；11-row test逐列驗 | `InterSubMod/src/core/RegionProcessor.cpp:107-111`；`InterSubMod/tests/test_verification_schema_v2.cpp:81-85` |
| `LOH_Subtype` | exact alias of `LOH_Subtype_LegacyVC` | Stage④同時寫；CSV前逐列 equality guard | `InterSubMod/src/core/RegionProcessor.cpp:1839-1841`、`InterSubMod/src/core/RegionProcessor.cpp:1966-1968` |
| `Subclone_ID`／`RegionResult::subclone_id` | exact alias of `RegionStratum_ID`，unassigned同為 -1 | commit同值；CSV前 equality guard | `InterSubMod/include/core/RegionProcessor.hpp:311-316`；`InterSubMod/src/core/SubcloneAnalyzer.cpp:157-166`；`InterSubMod/src/core/RegionProcessor.cpp:1970-1973` |
| `SubcloneSummary::subclone_id` | fixed slot ID，不是 cellular identity | summary index初始化成相同 ID；artifact writer驗 `summary[index].subclone_id==index` | `InterSubMod/include/core/SubcloneAnalyzer.hpp:83-89`；`InterSubMod/src/core/RegionProcessor.cpp:310-317` |
| `SubcloneResult::n_subclones` | exact mirror of distinct `n_occupied_region_strata` | analyze同時賦值 | `InterSubMod/include/core/SubcloneAnalyzer.hpp:119-125`；`InterSubMod/src/core/SubcloneAnalyzer.cpp:197-200` |
| `RegionAsmProfile::verification_class` | legacy-v1 compatibility alias，extraction讀 `verification_class_legacy` | canonical current/legacy另有明確欄位 | `InterSubMod/include/core/SubcloneAnalyzer.hpp:61-65`；`InterSubMod/src/core/SubcloneAnalyzer.cpp:247-250` |
| legacy summary counters／`dominant_verification_class` | mirror legacy leaf counts／legacy dominant ties | current v2另有完整11-leaf欄位 | `InterSubMod/include/core/SubcloneAnalyzer.hpp:101-113`；`InterSubMod/src/core/SubcloneAnalyzer.cpp:325-335` |
| class name `SubcloneAnalyzer` | compatibility façade，實作 deterministic region stratification | header明文標示 one compatibility release | `InterSubMod/include/core/SubcloneAnalyzer.hpp:157-179` |
| `subclone_structure.txt` | deprecation stub指向 `region_stratification_summary.txt` | valid/non-valid都原子覆寫同格式 stub | `InterSubMod/src/core/RegionProcessor.cpp:383-390` |
| `SubcloneAnalyzer::write_report()/write_assignments_tsv()` | historical public API compatibility writers | 非 canonical pipeline writer；open失敗只寫 stderr、無 bool error | `InterSubMod/src/core/SubcloneAnalyzer.cpp:340-376` |

**RISK**：source註解只說「one compatibility release」，但目前四個 core files沒有移除版本、deadline或 compile-time deprecation attribute；移除必須另有 migration issue/consumer inventory。

## 7. 三支新 C++ test files：17 tests 逐項 evidence

### 7.1 Verification schema v2（4 tests）

| Test／行號 | 直接證明 | 不直接證明 |
|---|---|---|
| `AuthoritativeElevenClassTruthTable`，`InterSubMod/tests/test_verification_schema_v2.cpp:30-87` | 11 classes、paths、Significant、schema=2、legacy保存、四 evidence fields、LIVE derivation、v1 projection | Stage④整合 block與CSV I/O |
| `OrderedPrecedenceIsStableUnderCollisions`，`InterSubMod/tests/test_verification_schema_v2.cpp:89-123` | cluster>LOH>within-HP>Δβ>PERMANOVA>HP-AUC 的代表 collisions | 每個可能 boolean組合的 exhaustive proof |
| `LegacyLohSubtypeAndDeprecatedAliasSourceAreExact`，`InterSubMod/tests/test_verification_schema_v2.cpp:125-133` | canonical LOH五值與 unknown fail closed | CSV alias equality guard |
| `UnknownLegacyClassFailsClosedBeforeClassification`，`InterSubMod/tests/test_verification_schema_v2.cpp:135-139` | invalid legacy不能被分類成 Noise或其他 class | process-level error publication |

### 7.2 RegionStratifier（8 tests）

| Test／行號 | 直接證明 | 不直接證明 |
|---|---|---|
| `EmptyGuardWinsEvenWhenMinimumIsZero`，`InterSubMod/tests/test_region_stratifier.cpp:52-61` | 0 profiles、min=0仍 insufficient；counts=0；storage=4 | process terminal artifact |
| `EligibilityBoundaryAtFortyNineAndFiftyIsExact`，`InterSubMod/tests/test_region_stratifier.cpp:63-77` | default 49不跑、50全 assign | 50筆多strata分布 |
| `OccupiedCountIsDistinctAndStorageAlwaysHasFourSlots`，`InterSubMod/tests/test_region_stratifier.cpp:79-104` | only 3、0+3、only 2 的 occupied distinct count與四槽 storage | 所有四槽同時 occupied |
| `PrecedenceAndStrictThresholdsRemainUnchanged`，`InterSubMod/tests/test_region_stratifier.cpp:106-118` | 0.05/0.10 strict boundary與 LOH>sample>HP | BED-only LOH碰撞 |
| `SharedEligibilityPreservesOriginalResultIndices`，`InterSubMod/tests/test_region_stratifier.cpp:120-145` | reads/CpG/computed eligibility與 original indices | `success=false` 專例 |
| `ValidateAllThenCommitAllPreventsPartialWriteback`，`InterSubMod/tests/test_region_stratifier.cpp:147-180` | duplicate、region mismatch、out-of-range零 partial mutation；valid全 commit；`Subclone_ID` alias | validator所有其他 error strings |
| `CurrentAndLegacyLeafCountersConserveUnknownsAndTies`，`InterSubMod/tests/test_region_stratifier.cpp:182-223` | 11 current leaves、4 legacy leaves、unknown buckets、守恆、warnings、deterministic ties | summary text serialization完整內容 |
| `CurrentClassRequiresExplicitSchemaVersionTwo`，`InterSubMod/tests/test_region_stratifier.cpp:225-240` | version!=2時 current進 unknown且 warning，legacy仍計數 | mixed versions across multiple rows 的獨立 fail |

### 7.3 Artifacts（5 tests）

| Test／行號 | 直接證明 | 不直接證明 |
|---|---|---|
| `StatusHasExactHeaderAndRfc3339UtcTimestamp`，`InterSubMod/tests/test_region_stratification_artifacts.cpp:39-54` | exact 9-column header、defaults、RFC3339 UTC | lifecycle status-first/last ordering |
| `ValidArtifactsUseStableKeyAndOnlyOccupiedSections`，`InterSubMod/tests/test_region_stratification_artifacts.cpp:56-109` | stable identity row、只列occupied summary、stub wording、success後無tmp | duplicate stable key rejection |
| `NonValidRunsOverwritePriorValidArtifactsWithoutStaleRows`，`InterSubMod/tests/test_region_stratification_artifacts.cpp:111-158` | 三個non-valid states都清除舊rows/sections並更新stub reason | mid-write crash/fault injection |
| `ContradictoryStatusAndRowsFailBeforePublication`，`InterSubMod/tests/test_region_stratification_artifacts.cpp:160-181` | contradictory valid rows/status在建立output dir前 fail | 所有 status/artifact validation reasons |
| `AtomicRenameFailureRemovesOwnedTempFile`，`InterSubMod/tests/test_region_stratification_artifacts.cpp:184-200` | final path被目錄占用時 rename明確失敗；guard移除本輪 `.tmp.`；既有 final directory不受影響 | private CSV/statistics writer各自的 failure injection |

## 8. Fail-closed reason catalog

### 8.1 Classification／LOH

| Reason／exception | Trigger | Source |
|---|---|---|
| `INVALID_VERIFICATION_CLASS_LEGACY: <value>` | classifier或LOH helper收到四態外字串 | `InterSubMod/src/core/RegionProcessor.cpp:50-53`、`InterSubMod/src/core/RegionProcessor.cpp:115-120` |
| `UNREACHABLE_VERIFICATION_CLASS_LEGACY` | 合法 enum卻未命中 helper branch；理論防衛分支 | `InterSubMod/src/core/RegionProcessor.cpp:121-126` |

### 8.2 Assignment validation

完整 exact reason set 位於 `InterSubMod/src/core/SubcloneAnalyzer.cpp:112-147`：

`ASSIGNMENT_COUNT_MISMATCH`、`ELIGIBLE_REGION_COUNT_MISMATCH`、`RESULT_INDEX_OUT_OF_RANGE`、`DUPLICATE_RESULT_INDEX`、`INELIGIBLE_ASSIGNMENT`、`REGION_ID_MISMATCH`、`INVALID_STRATUM`、`ELIGIBLE_REGION_UNASSIGNED`。任一 reason都在 commit mutation前返回。

### 8.3 Status validation

完整 exact reason set 位於 `InterSubMod/src/core/RegionProcessor.cpp:149-184`：

`INVALID_REGION_STATUS_SCHEMA_VERSION`、`EMPTY_REGION_STATUS_REASON`、`NEGATIVE_REGION_STATUS_COUNT`、`INVALID_REGION_STATUS_OCCUPIED_COUNT`、`INCONSISTENT_VALID_REGION_STATUS`、`NONVALID_REGION_STATUS_HAS_ASSIGNMENTS`。

### 8.4 Artifact validation

| Reason | Trigger | Source |
|---|---|---|
| `REGION_ARTIFACT_ROW_COUNT_MISMATCH` | valid count不等於 rows；或 non-valid仍帶 rows | `InterSubMod/src/core/RegionProcessor.cpp:270-274` |
| `REGION_ARTIFACT_SUMMARY_SLOT_COUNT_MISMATCH` | summaries不是四槽 | `InterSubMod/src/core/RegionProcessor.cpp:275-278` |
| `INVALID_REGION_ARTIFACT_ROW` | schema、ID、chrom或stratum越界 | `InterSubMod/src/core/RegionProcessor.cpp:283-289` |
| `REGION_ARTIFACT_ENUM_MISMATCH` | ID/label/reason不一致 | `InterSubMod/src/core/RegionProcessor.cpp:290-295` |
| `DUPLICATE_REGION_ARTIFACT_RESULT_INDEX` | 重複 original result index | `InterSubMod/src/core/RegionProcessor.cpp:296-299` |
| `DUPLICATE_REGION_ARTIFACT_STABLE_KEY` | 重複 `RegionID/Chr/Pos/Ref/Alt` | `InterSubMod/src/core/RegionProcessor.cpp:300-305` |
| `REGION_ARTIFACT_SUMMARY_COUNT_MISMATCH` | slot ID/count與 rows不一致 | `InterSubMod/src/core/RegionProcessor.cpp:310-317` |
| `REGION_ARTIFACT_OCCUPIED_COUNT_MISMATCH` | distinct occupied與status不一致 | `InterSubMod/src/core/RegionProcessor.cpp:318-323` |

### 8.5 Main CSV／lifecycle publication

1. CSV逐列 schema/alias/sentinel reasons：`InterSubMod/src/core/RegionProcessor.cpp:1958-1995`。
2. generic artifact temp path與RAII cleanup：`InterSubMod/src/core/RegionProcessor.cpp:131-136`、`InterSubMod/src/core/RegionProcessor.cpp:186-228`；CSV temp lifecycle：`InterSubMod/src/core/RegionProcessor.cpp:1860-1869`、`InterSubMod/src/core/RegionProcessor.cpp:2191-2205`；statistics temp lifecycle：`InterSubMod/src/core/RegionProcessor.cpp:2208-2251`。
3. artifact主寫入失敗轉 `ARTIFACT_WRITE_FAILED` 並嘗試 failed overwrite：`InterSubMod/src/core/RegionProcessor.cpp:1163-1191`。
4. CSV失敗轉 `SIGNIFICANCE_CSV_WRITE_FAILED` 並嘗試 failed overwrite/status：`InterSubMod/src/core/RegionProcessor.cpp:1194-1206`。
5. terminal status寫入失敗直接 throw `REGION_STRATIFICATION_TERMINAL_STATUS_WRITE_FAILED`；先前 marker仍不會被誤寫成 terminal valid：`InterSubMod/src/core/RegionProcessor.cpp:1208-1213`。

## 9. L2 工程推論

1. **分類 split 保留 binary retention**：兩個被拆開的 positive legacy branches都顯式 `Significant=true`，CSV又只投影該 boolean，所以 rename/split本身不改這兩群的 binary outcome。這不是 frozen corpus逐列證明。
2. **status是唯一跨檔 commit marker**：owned artifacts各自 atomic，但依序寫；只有最後的 terminal status能表示同一輪所有輸出已 durable。consumer若跳過status，無法獲得同樣的 fail-closed保證。
3. **unknown class是可觀察 warning，不是 fatal**：summary把 unknown放獨立leaf並增加warning；status可為 `VALID_WITH_WARNINGS`。這避免 silent Noise fallback，但下游仍需明確決定是否接受 warning。
4. **四槽與occupied count是兩個維度**：四槽確保ID-indexed storage安全；occupied count只表示本輪非空strata數。`n_subclones` 鏡像不能被解讀為推定cellular populations。
5. **write-back identity與serialized stable key分工不同**：vector index保證回寫正確row；stable key保證輸出可比與唯一。兩者都需通過，不能相互替代。

## 10. 未解風險與補強優先序

| 優先 | RISK | 現有保護 | 建議補強 |
|---:|---|---|---|
| P1 | 無完整 `process_all_regions()` fault-injection test證明 status-first/artifact/CSV/status-last在每個失敗點的實際filesystem狀態 | source控制流與helper tests | 加小型 fake writer／temporary output integration test，不需大型 BAM |
| P1 | 三個 owned artifacts不是跨檔transaction；直接 consumer可能讀到不同世代 | `RUN_IN_PROGRESS` commit marker | validator與consumer強制先讀且接受唯一 terminal `VALID` |
| P1 | duplicate stable key、summary count mismatch、部分 status reasons沒有逐reason單元測試 | writer在publication前驗證 | 參數化 negative tests覆蓋全部 reason strings |
| P2 | `write_significance_summary()` 若CSV已rename但statistics後續失敗，filesystem可暫存本輪CSV與failed artifacts混合；正確性依賴failed status | pipeline會改寫failed status/artifacts | consumer一律status-first；加此精確failure-order test |
| P2 | `SubcloneAnalyzer::write_report()/write_assignments_tsv()` 仍是public compatibility API，非atomic且open失敗只stderr | canonical pipeline不使用它們發布正式artifacts | call-site inventory後deprecate/remove，或改回傳bool |
| P2 | compatibility「one release」沒有機器可讀移除版本 | 註解與deprecation stub | 建立明確issue、target schema version與consumer completion gate |
| P3 | C++ tests驗 live decision，但不證 frozen corpus逐stable-key invariance | pure truth table 11/11 | 另引用 frozen golden migration validator；不可混寫為本文件證據 |

**已關閉風險（L1）**：generic text artifact、`significance_summary.csv` 與 `significance_statistics.txt` 現在都在建立 temp 後立即配置 `TempFileGuard`，只有成功 rename 才 `mark_published()`；forced atomic-rename failure test證明 owned `.tmp.` 被移除且 final path不被刪除。來源：`InterSubMod/src/core/RegionProcessor.cpp:186-228`、`InterSubMod/src/core/RegionProcessor.cpp:1860-1869`、`InterSubMod/src/core/RegionProcessor.cpp:2191-2251`；測試：`InterSubMod/tests/test_region_stratification_artifacts.cpp:184-200`。這不改變 remaining cross-file limitation：多個 final artifacts仍逐檔發布，必須以最後 terminal status作 commit marker。

## 11. Completeness matrix

| 指定項目 | Code source | Test source | 本文件節 | 完整度 |
|---|---|---|---|---|
| 11 classes | `InterSubMod/src/core/RegionProcessor.cpp:50-112` | `InterSubMod/tests/test_verification_schema_v2.cpp:30-87` | `3.2` | L1 complete |
| four evidence fields | `InterSubMod/src/core/RegionProcessor.cpp:56-60`、`InterSubMod/src/core/RegionProcessor.cpp:1822-1837` | `InterSubMod/tests/test_verification_schema_v2.cpp:74-80` | `3.3` | L1 complete |
| Significant invariance | `InterSubMod/src/core/RegionProcessor.cpp:63-105`、`InterSubMod/src/core/RegionProcessor.cpp:2001-2007` | `InterSubMod/tests/test_verification_schema_v2.cpp:30-72` | `3.4` | C++ mapping complete；corpus proof external |
| LOH legacy derivation | `InterSubMod/src/core/RegionProcessor.cpp:115-126`、`InterSubMod/src/core/RegionProcessor.cpp:1839-1841` | `InterSubMod/tests/test_verification_schema_v2.cpp:125-133` | `3.5` | L1 complete |
| eligibility | `InterSubMod/src/core/SubcloneAnalyzer.cpp:97-100`、`InterSubMod/src/core/SubcloneAnalyzer.cpp:226-254` | `InterSubMod/tests/test_region_stratifier.cpp:120-145` | `4.2` | L1 complete |
| 0/49/50 | `InterSubMod/src/core/SubcloneAnalyzer.cpp:170-200` | `InterSubMod/tests/test_region_stratifier.cpp:52-77` | `4.5` | L1 complete |
| four fixed slots | `InterSubMod/include/core/SubcloneAnalyzer.hpp:14-16`；`InterSubMod/src/core/SubcloneAnalyzer.cpp:170-201`、`InterSubMod/src/core/SubcloneAnalyzer.cpp:270-275` | `InterSubMod/tests/test_region_stratifier.cpp:79-104` | `4.5` | L1 complete |
| stable key | `InterSubMod/src/core/RegionProcessor.cpp:280-306`、`InterSubMod/src/core/RegionProcessor.cpp:1128-1161` | `InterSubMod/tests/test_region_stratification_artifacts.cpp:56-92` | `5.3` | serialization tested；duplicate path gap |
| status-first/last | `InterSubMod/src/core/RegionProcessor.cpp:911-921`、`InterSubMod/src/core/RegionProcessor.cpp:1163-1213` | helper tests only | `5.4` | L1 source complete；integration gap |
| stale overwrite | `InterSubMod/src/core/RegionProcessor.cpp:325-390`、`InterSubMod/src/core/RegionProcessor.cpp:1163-1206` | `InterSubMod/tests/test_region_stratification_artifacts.cpp:111-158` | `5.4` | helper-level complete |
| owned temp cleanup | `InterSubMod/src/core/RegionProcessor.cpp:186-228`、`InterSubMod/src/core/RegionProcessor.cpp:1860-1869`、`InterSubMod/src/core/RegionProcessor.cpp:2191-2251` | `InterSubMod/tests/test_region_stratification_artifacts.cpp:184-200` | `5.1`、`8.5`、`10` | rename-failure path L1 complete；cross-file transaction不在此 claim |
| compatibility aliases | `InterSubMod/include/core/RegionProcessor.hpp:311-370`；`InterSubMod/include/core/SubcloneAnalyzer.hpp:83-179` | verification/commit/artifact tests | `6` | current release complete；removal plan gap |
| fail-closed reasons | `InterSubMod/src/core/RegionProcessor.cpp:50-53`、`InterSubMod/src/core/RegionProcessor.cpp:149-228`、`InterSubMod/src/core/RegionProcessor.cpp:260-323`、`InterSubMod/src/core/RegionProcessor.cpp:1958-1995`；`InterSubMod/src/core/SubcloneAnalyzer.cpp:112-168` | selected negative paths | `8` | catalog complete；some reasons untested |
| no biological clone claim | `InterSubMod/include/core/SubcloneAnalyzer.hpp:18-23`、`InterSubMod/include/core/SubcloneAnalyzer.hpp:157-162`；stub `InterSubMod/src/core/RegionProcessor.cpp:383-390` | stub exact text `InterSubMod/tests/test_region_stratification_artifacts.cpp:100-104` | `1.3`、`6` | L1 complete |

## 12. 參考來源索引

### 12.1 四個 core files

1. `InterSubMod/include/core/RegionProcessor.hpp`
   - Verification structs：`InterSubMod/include/core/RegionProcessor.hpp:31-65`
   - artifact/status records：`InterSubMod/include/core/RegionProcessor.hpp:67-102`
   - RegionResult schema fields/defaults：`InterSubMod/include/core/RegionProcessor.hpp:311-370`、`InterSubMod/include/core/RegionProcessor.hpp:524-572`
2. `InterSubMod/include/core/SubcloneAnalyzer.hpp`
   - schema constants/enums：`InterSubMod/include/core/SubcloneAnalyzer.hpp:14-40`
   - profile/assignment/summary/result：`InterSubMod/include/core/SubcloneAnalyzer.hpp:45-136`
   - compatibility façade：`InterSubMod/include/core/SubcloneAnalyzer.hpp:138-179`
3. `InterSubMod/src/core/RegionProcessor.cpp`
   - classifier/LOH：`InterSubMod/src/core/RegionProcessor.cpp:44-126`
   - atomic temp guard／status／artifact writers：`InterSubMod/src/core/RegionProcessor.cpp:131-391`
   - run lifecycle：`InterSubMod/src/core/RegionProcessor.cpp:911-1217`
   - Stage④ integration：`InterSubMod/src/core/RegionProcessor.cpp:1794-1842`
   - significance serialization／temp publication：`InterSubMod/src/core/RegionProcessor.cpp:1856-2254`
4. `InterSubMod/src/core/SubcloneAnalyzer.cpp`
   - fixed enums/counters：`InterSubMod/src/core/SubcloneAnalyzer.cpp:18-95`
   - eligibility/precedence：`InterSubMod/src/core/SubcloneAnalyzer.cpp:97-110`
   - validate/commit：`InterSubMod/src/core/SubcloneAnalyzer.cpp:112-168`
   - analyze/profiles/summaries：`InterSubMod/src/core/SubcloneAnalyzer.cpp:170-338`
   - compatibility writers：`InterSubMod/src/core/SubcloneAnalyzer.cpp:340-376`

### 12.2 三支新 C++ test files（17 tests）

1. `InterSubMod/tests/test_verification_schema_v2.cpp:30-139` — 4 tests。
2. `InterSubMod/tests/test_region_stratifier.cpp:52-240` — 8 tests。
3. `InterSubMod/tests/test_region_stratification_artifacts.cpp:39-200` — 5 tests。

### 12.3 驗證狀態

本 source snapshot 的 Post-C++ verification 已執行：

- `cmake --build build -j2`：exit 0。
- `./build/bin/run_tests --gtest_filter=RegionStratificationArtifactsTest.AtomicRenameFailureRemovesOwnedTempFile`：1/1 passed。
- `./build/bin/run_tests`：251/251 passed（36 suites），其中上述三支新 C++ test files 共 17 tests。
- `ctest --test-dir build --output-on-failure`：251/251 passed，0 failed。

這些結果證明目前 source 可建置且既有 C++ test suite 全綠；它們不擴張本文件在 `10` 列出的 coverage 邊界。
