<!--
建立時間: 2026-07-15
目標: G4 多樣本一致性與 reproducibility；G5 業界級可外部驗證工程
處理範圍: Verification schema v2、RegionStratification schema v1、LOH provenance、74 個 consumer、frozen 14-file offline migration、測試與 fail-closed publication
關聯檔案:
  - /big7_disk/liaoyoyo2001/external_validation/subclone_reconstruction_method_comparison_20260714/20260714_InterSubMod分類與RegionStratifier修復任務規格_01.md
  - InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md
  - InterSubMod/research/20260714_verification_region_stratifier_schema_repair/20260714_CPP分類與RegionStratifier流程來源對照_01.md
  - InterSubMod/research/20260714_verification_region_stratifier_schema_repair/20260714_PythonConsumer契約與遷移對照_01.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/schema_migrations/20260714_verification_v2/migration_summary.json
-->

# Verification schema v2 與 RegionStratifier 修復完整流程說明

> **TL;DR — P0→P1→P2 全部完成，工程驗收 PASS：11-class schema、RegionStratifier、LOH provenance、74/74 consumers 與 frozen 14 檔 328,697 列均已 fail-closed 驗證並原子發布；這證明 schema/serialization migration 正確，不等於生物學 clone 或 caller accuracy 驗證。**（影響：高；信心：工程高、科學外推不適用）

**敘述框架：SCQA（Situation → Complication → Question → Answer）+ claim-to-evidence mapping。**

| SCQA | 本次答案 |
|---|---|
| Situation | Stage④ 已有 per-region 統計與 legacy provenance，但 final `Strong` 混合兩條 evidence path；Phase D 又用 `Subclone*` 命名固定規則的 region categories。 |
| Complication | current/legacy class、LOH subtype、region assignment、downstream cohort 與 frozen output 可能被混讀；stale artifact 或 destination race 會讓失敗看似成功。 |
| Question | 能否不改 raw metrics 與 binary `Significant`，把語意拆清楚、讓所有 consumer 明確選 schema，並以 frozen corpus 證明逐列不漂移？ |
| Answer | 能。canonical `VerificationSchemaVersion=2` 與 `RegionStratificationSchemaVersion=1` 已落地；正式 migration 為 `VALID 14/14`，三個 invariance flags 全 PASS，獨立 hash/row read-back 亦全 PASS。 |

## 0. 裁決、範圍與可信度

### 0.1 最終裁決

| Gate | 裁決 | 最強直接證據 | 證據層 |
|---|---|---|---|
| P0 schema / classification | PASS | 11-class truth table、collision precedence、typed evidence、binary retention tests | L1 ⭐⭐⭐⭐⭐ |
| P0 RegionStratifier | PASS | eligibility 單一 predicate、0/49/50、四 fixed slots、validate-all-before-commit | L1 ⭐⭐⭐⭐⭐ |
| P1 provenance / consumers | PASS | LOH exact alias；bounded inventory 74 unique、74 `MIGRATED`、0 unresolved | L1 ⭐⭐⭐⭐⭐ |
| P1 stale artifact / lifecycle | PASS | status-first、owned artifact overwrite、terminal status-last、temp cleanup failure regression | L1 ⭐⭐⭐⭐⭐ |
| P2 frozen migration | PASS | `VALID`、14/14、328,697→328,697、0 conflicts、三項 invariance PASS | L1 ⭐⭐⭐⭐⭐ |
| 生物學 clone／caller accuracy | **未主張** | 本任務未重跑 BAM、統計、truth benchmark 或 clone inference | scope ceiling |

### 0.2 任務分類與研究五問

- **Task type**：E（Hotfix/Bugfix）+ B（Comprehensive validation）。因規格使用「完整、驗證、全 14 檔」，subset 不可作終局證據。
- **服務目標**：G4 reproducibility、G5 外部可驗證工程；間接保護 G2/G3 的後續分析語意。
- **Thread D read-level epigenetic**：非直接研究假說；本次只修 per-region schema 與 downstream contract。
- **Thread B NO-GO**：不在撤回範圍；沒有重新主張 methylation F1 增益。
- **KDE-corrected**：本次 offline migration 不重算 coverage/KDE；所有既有 raw token 原樣保存。
- **VCF caller AF**：未使用，也沒有從 merged AF 推導新結論。
- **長計算／C++／搬移／gate**：有 C++ build、frozen 328,697-row migration 與 atomic directory publish；沒有 BAM/read-matrix/statistics rerun，也沒有刪除或原地覆寫舊結果。

### 0.3 證據分級

| 層 | 意義 | 本報告例子 |
|---|---|---|
| L1 | 可由 source、test、manifest、hash 或輸出直接重現 | 251/251 C++、14/14 SHA、59,910/6,931 |
| L2 | 由 L1 控制流推論，前提明列 | status-last 對遵守 status-first consumer 提供 fail-closed commit marker |
| L3 | 本次資料觀察但未有獨立 truth | 本報告刻意不新增 L3 生物學 claim |
| L4 | 文獻／外部 benchmark | 本任務不以外部工具 claim 取代本 repo 驗證 |
| L5 | 未確認假設 | repo 外 precompiled ABI consumer 是否存在 |

## 1. 端到端資料流：從 Stage④ 到 frozen provenance

```text
per-region raw metrics + legacy four-state class
                 │
                 ▼
      Verification schema v2 classifier
      ├─ 11-class current taxonomy
      ├─ v1 deprecated view
      ├─ four typed evidence booleans/path
      ├─ Significant boolean projection
      └─ LOH legacy provenance + exact alias
                 │
                 ▼
       RegionStratifier schema v1
      ├─ shared eligibility predicate
      ├─ four ordered region strata
      ├─ validate complete assignment set
      ├─ commit by result_index
      └─ current/legacy counters + warnings
                 │
                 ▼
 fail-closed publication boundary
      FAILED/RUN_IN_PROGRESS status first
      → atomic owned artifacts
      → atomic significance CSV/statistics
      → terminal status last
                 │
        ┌────────┴─────────┐
        ▼                  ▼
 shared Python contract    frozen offline migration
 C2/L4/E/R1/P/LOH-L/H1    byte-preserving append-only
        │                  │
        ▼                  ▼
 74 bounded consumers      no-replace directory publish
                           + independent validator
```

判讀順序不可反轉：**consumer 必須先接受 terminal status/provenance，才可讀 data artifact**。看到檔案存在，不等於該 generation 有效。

## 2. Verification schema v2：最基本判斷與 precedence

### 2.1 First-match-wins 決策表

下表是唯一 authoritative 順序。前一列命中後不再往下；本次不改任何原始 predicate 或統計門檻，只拆開原 final `Strong` 的來源並把 evidence 明確序列化。

| 優先序 | 最小條件 | canonical `VerificationClass` | `EvidencePath` | `Significant` |
|---:|---|---|---|---:|
| 1a | cluster-first match 且 legacy=`Strong` | `Strong_Bidirectional` | `BIDIRECTIONAL` | true |
| 1b | cluster-first match 且 legacy=`Subclone` | `ClusterFirstOnly` | `CLUSTER_FIRST_ONLY` | true |
| 2 | 尚未命中，且 `potential_loh && (hp_auc_struct || dbeta_sig)` | `LOH-Structure` | `LOH_STRUCTURE` | true |
| 3 | 尚未命中，且 within-HP clean multigroup | `MultiGroupNoLabel` | `WITHIN_HP_MULTIGROUP` | true |
| 4 | 尚未命中，且 Δβ label shift 顯著 | `LabelShift` | `LABEL_SHIFT` | true |
| 5 | 尚未命中，且 clean location PERMANOVA | `PermanovaLocation` | `PERMANOVA_LOCATION` | true |
| 6 | 尚未命中，且 HP-AUC structure | `StructureNoLabel` | `HP_AUC_STRUCTURE_NO_LABEL` | true |
| 7 | 沒有 clean branch，只有 dispersion structure | `DispersionStructure` | `DISPERSION_STRUCTURE` | false |
| 8 | 尚未命中，達 uniform noise threshold | `Noise_Uniform` | `NOISE_UNIFORM` | false |
| 9 | 尚未命中，達 chaotic noise threshold | `Noise_Chaotic` | `NOISE_CHAOTIC` | false |
| 10 | 其他 | `Noise_Uncorrelated` | `NOISE_UNCORRELATED` | false |

為何列為 11 class：第 1 有 1a/1b 兩個互斥 leaf，其後有九個 leaf；`LabelShift` 在 frozen corpus 本次觀察為 0，不代表它可從 enum 或報表拿掉。

### 2.2 欄位設定與不可混用的語意

| 欄位 | live 值域／來源 | offline frozen policy | Fail-closed 判斷 |
|---|---|---|---|
| `VerificationSchemaVersion` | integer `2` | 固定 `2` | missing/mixed/非 2 不可當 C2 |
| `VerificationClass` | 上表 11 enum | 由 v1 current + legacy pair 映射 | unknown 獨立 bucket 或拒絕，不可歸 Noise |
| `VerificationClass_V1_Deprecated` | v1 10-class view | 原 `VerificationClass` raw token | 只供一個 migration release |
| `VerificationClass_Legacy` | `Strong/Subclone/Weak/Noise` | 原欄，缺失即 unresolved | 不可從 current 猜回 legacy |
| `LabelFirstSupport` | legacy ∈ `{Strong,Weak}` | 可由 legacy 推得 | typed bool；missing/unknown 不得當 false |
| `ClusterFirstSupport` | legacy ∈ `{Strong,Subclone}` | 可由 legacy 推得 | typed bool；不是 cellular clone evidence |
| `WithinHPSupport` | live boolean | `NA` | offline 不從 class 猜 boolean |
| `DispersionWarning` | live HP/allele warning OR | `NA` | offline 不從 class 猜 boolean |
| `EvidencePath` | 命中的 ordered branch | 明列 pair mapping | 必須與 class、schema 相容 |
| `EvidenceDerivation` | `LIVE` | `LEGACY_CLASS` | 不可把推導值冒充 live evidence |

核心語意界線：

- `Strong_Bidirectional` 是 label-first 與 cluster-first 路徑同時支持的 schema class。
- `ClusterFirstOnly` 是 region-level cluster-first association，不是 cellular subclone。
- `Significant` 是既有 binary retention 的投影；兩個拆分 leaf 都維持 true。
- `VerificationClass` 不是 truth label；truth 必須來自明示的 truth 欄或外部 benchmark。

### 2.3 LOH provenance

`LOH_Subtype_LegacyVC` 固定值域為：

```text
None | LOH_Noise | LOH_Weak | LOH_Strong | LOH_Subclone
```

判斷最小化為：

1. `potential_loh=false` → `None`。
2. `potential_loh=true` → 只依 `VerificationClass_Legacy` 映射四個 `LOH_*`。
3. 其他 potential-LOH token 或 legacy class → fail closed。
4. deprecated `LOH_Subtype` 必須逐列 exact alias；本次不建立語意不清的 `LOH_Subtype_Current`。

## 3. RegionStratifier schema v1：設定、門檻與 write-back

### 3.1 Eligibility 是唯一共用 predicate

一列只有同時滿足下列四條才 eligible：

```text
result.success
&& num_reads >= 10
&& num_cpgs >= 3
&& significance_computed
```

Analyzer、assignment validator 與 caller write-back 共用此 predicate，避免「分析時一套、回寫時另一套」。

### 3.2 四個固定 region strata

| ID | Label | Reason | Ordered predicate |
|---:|---|---|---|
| 2 | `LohFlagged` | `LOH_FLAGGED` | `potential_loh || loh_bed_overlap` |
| 3 | `HighSampleAsm` | `HIGH_SAMPLE_ASM` | 非 LOH 且 `abs(sample_asm_delta) > 0.10` |
| 1 | `HighHpAsm` | `HIGH_HP_ASM` | 前兩者否且 `abs(hp_merged_delta) > 0.05` |
| 0 | `BaselineLowAsm` | `BASELINE_LOW_ASM` | 其他 eligible row |

precedence 固定為 `LOH → sample ASM → HP ASM → baseline`。門檻是嚴格大於；剛好 0.10 或 0.05 不會進 high branch。

### 3.3 0／49／50 與 sparse ID

| Eligible count（default min=50） | Status | Assignments | Occupied count |
|---:|---|---:|---:|
| 0 | `INSUFFICIENT_REGIONS / NO_ELIGIBLE_REGIONS` | 0 | 0 |
| 49 | `INSUFFICIENT_REGIONS / BELOW_MIN_REGIONS` | 0 | 0 |
| 50 | `VALID / OK` 或 `VALID_WITH_WARNINGS` | 50 | distinct non-empty IDs |

Storage 永遠配置四個 ID slots；`n_occupied_region_strata` 另計 distinct non-empty IDs。因此 only `{3}` 是 1，不是 `max_id+1=4`；`{0,3}` 是 2。deprecated `n_subclones` 只是一個 corrected mirror，不可解讀成 clone 數。

### 3.4 Validate-all-before-commit

每筆 assignment 帶 `{result_index, region_id, stratum}`：

1. 先驗整體 count 與重新計算的 eligible count。
2. 再驗 `result_index` range、index 唯一、該列 eligible、`region_id` 一致、stratum typed enum。
3. 再確認沒有 eligible row 漏 assignment。
4. **全部通過後**才依 `result_index` mutation；任何錯誤都零筆部分 write-back。

可能的 machine reasons 包括 `ASSIGNMENT_COUNT_MISMATCH`、`ELIGIBLE_REGION_COUNT_MISMATCH`、`RESULT_INDEX_OUT_OF_RANGE`、`DUPLICATE_RESULT_INDEX`、`INELIGIBLE_ASSIGNMENT`、`REGION_ID_MISMATCH`、`INVALID_STRATUM`、`ELIGIBLE_REGION_UNASSIGNED`。

### 3.5 Sentinel 與 aliases

未 assignment 的 canonical 值：

| 情境 | ID | Label | Reason |
|---|---:|---|---|
| valid phase 的單列不 eligible | -1 | `Unassigned` | `INELIGIBLE_REGION` |
| paired 但不足 | -1 | `Unassigned` | `INSUFFICIENT_REGIONS` |
| tumor-only | -1 | `Unassigned` | `NOT_APPLICABLE_TUMOR_ONLY` |
| phase failed | -1 | `Unassigned` | `FAILED` |

`Subclone_ID` 在相容期必須 exact 等於 `RegionStratum_ID`；未 assignment 時兩者都為 -1。這是 source compatibility，不是 biological alias。

## 4. Fail-closed publication 與 stale artifact 防護

### 4.1 Lifecycle

```text
[1] atomic status = FAILED / RUN_IN_PROGRESS
                    │
                    ▼
[2] compute + complete-set validation
                    │
          ┌─────────┴─────────┐
          │                   │
        failure             success
          │                   │
          ▼                   ▼
 overwrite owned       atomic assignment TSV
 artifacts as empty    atomic summary + deprecated stub
 + FAILED status       atomic significance CSV + statistics
          │                   │
          └─────────┬─────────┘
                    ▼
        [3] terminal status written last
```

`region_stratification_status.tsv` 是 commit marker；固定欄位為：

```text
status, reason, schema_version, eligible_region_count,
min_regions_required, assignment_count, n_occupied_region_strata,
warning_count, generated_at
```

只接受四種 status：`VALID`、`INSUFFICIENT_REGIONS`、`NOT_APPLICABLE_TUMOR_ONLY`、`FAILED`。非 `VALID` 的 assignment 與 occupied count 必須都是 0。

### 4.2 Canonical owned artifacts

| Artifact | `VALID` | 非 `VALID` |
|---|---|---|
| `region_stratification.tsv` | header + 一列／assignment | canonical header only |
| `region_stratification_summary.txt` | 標題、status/reason、occupied sections | 標題、status/reason，沒有 sections |
| `subclone_structure.txt` | deprecation stub | 同一 stub，帶本次 status/reason |

若 artifact validation 在碰檔案前即失敗，caller 仍會以四個空 slots 和 `FAILED/ARTIFACT_WRITE_FAILED` 覆寫三個 owned artifacts，再寫 FAILED status；上一輪 valid rows 不會留在原位冒充本輪資料。

### 4.3 Atomic temp cleanup

所有 owned text、significance CSV 與 statistics 都寫同目錄唯一 temp，再 rename。`TempFileGuard` 只有在 rename 成功後才解除；writer、flush、close 或 rename 任何失敗會呼叫 `std::remove()` 清理**本輪 owned temp**，不動 final path或其他檔案。新增的 rename-collision regression 已直接證明一般失敗路徑不殘留；destructor不拋出或回報 `std::remove()` 自身的失敗，外部同時改動目錄權限／namespace屬明列的 hostile-race範圍，不影響 terminal status仍為fail-closed。

### 4.4 仍需遵守的讀取規則

三個 artifact 不是跨檔案 filesystem transaction。其一致性來自 status-first protocol：reader 必須先讀 terminal `VALID`，讀完資料後若需要防 concurrent writer，再重讀 status／generation context。若未來允許同一 output directory 多 writer 或要求 immutable generation，應升級為 generation directory + manifest/symlink swap，而不是默默擴寫目前 fixed schema。

## 5. Downstream consumer migration

### 5.1 單一 shared contract helper

| Code | Consumer 應讀什麼 | 預設拒絕什麼 |
|---|---|---|
| C2 | schema v2 current 11-class | missing/mixed version、unknown current |
| L4 | legacy four-state class | 缺 legacy、由 current 猜 legacy |
| E | typed evidence booleans/path | string truthiness、missing 當 false |
| R1 | region status + assignments | 非 VALID 有 assignment、count/alias/sentinel 矛盾 |
| P | 完整 pass-through provenance | drop version/class/evidence 或混版 merge |
| LOH-L | canonical LOH legacy subtype + exact alias | missing/unknown 當 `None` |
| H1 | unversioned historical replication | production default 偷開 fallback |

production/default 不開 `allow_unversioned_*`。H1 只能用明確 CLI flag，並在輸出 metadata 留下 historical warning。

### 5.2 Inventory 結果

- 62 個 VerificationClass paths + 23 個 LOH paths，去重後 **74 unique consumers**。
- **74 `MIGRATED`，0 unresolved**；11 個 path 同時涉及兩類 contract。
- 74 列逐檔 ledger 位於 `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv`。
- 完整 old assumption、new field、guard 與 test seam 位於 `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_inventory.md`。

四支規格點名 consumer 的結果：

| Consumer | 修正後語意 |
|---|---|
| `tools/compare_vcf_results.py` | v2 顯示完整 11 classes + unknown；legacy view另選，不把零頻 class丟掉 |
| `tools/find_verification_candidates.py` | historical report明確讀 `VerificationClass_Legacy`，缺欄 fail |
| `scripts/analysis/analyze_methyl_cluster_allele_cooccurrence.py` | historical Strong/Noise只在 explicit H1 mode；記錄 selection field/value與 unknown count |
| `scripts/analysis/compare_subclone_validation.py` | support 改成 legacy/evidence contract，輸出不再稱 cellular clone |

另外三支既有 untracked historical consumer只有修改前 SHA prefix，沒有可誠實回推的完整 before digest；目前完整 after SHA 已保存。這是 provenance limitation，不以猜測補值。

## 6. P2 offline migration：byte preservation 與原子發布

### 6.1 輸入、命令、輸出

**輸入 manifest**：

```text
/big7_disk/liaoyoyo2001/external_validation/subclone_reconstruction_method_comparison_20260714/FROZEN_CORPUS_MANIFEST.tsv
SHA-256: 9ab4428e28d137b3c15e811563678d006e6258ff6d3f86e718cc4f56ead2e0ea
```

**輸入 root**：

```text
/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/20260620_allsample_subcluster_split/results
```

**唯一 formal 命令**：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 scripts/migration/migrate_verification_schema_v2.py \
  --manifest /big7_disk/liaoyoyo2001/external_validation/subclone_reconstruction_method_comparison_20260714/FROZEN_CORPUS_MANIFEST.tsv \
  --input-root /big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/20260620_allsample_subcluster_split/results \
  --output-dir output/schema_migrations/20260714_verification_v2 \
  --require-output-dir-absent
```

**實際 stdout / exit**：

```text
migration status=VALID valid_files=14 failed_files=0 input_rows=328697 output_rows=328697 output=output/schema_migrations/20260714_verification_v2
exit code: 0
```

**實際輸出**：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/schema_migrations/20260714_verification_v2
```

### 6.2 為何不是 dataframe round-trip

Migration 以 bytes 解析 CSV record/token：

1. 先核對 14 個 manifest SHA-256 與 row count；任一不符即停止。
2. 用 CSV grammar 辨識 record boundaries，但保留每個原 token bytes 與原 newline。
3. 只替換既有 `VerificationClass` token；原 token exact 複製到 `VerificationClass_V1_Deprecated`。
4. 其他既有 field token、quoting、NA 拼法與 row order不變；九個新欄 append 在尾端。
5. 逐 stable key `(RegionID,Chr,Pos,Ref,Alt)` 核對唯一、相同與 `Significant` raw token一致。
6. temp output fsync後再讀回做完整 audit，最後才進 staging canonical filename。

### 6.3 Directory-level fail closed

1. Authoritative output 在 formal run 前必須不存在；dangling symlink 也視為存在。
2. `mkdtemp` staging 綁定 resolved parent、prefix、EUID、device/inode與parent identity。
3. 全 14 檔、reports、manifest與 invariants 完成後，才呼叫 Linux `renameat2(RENAME_NOREPLACE)`。
4. destination race 的 `EEXIST/ENOTEMPTY` 直接拒絕；kernel/filesystem不支援也拒絕，沒有 `os.rename` fallback。
5. 例外 cleanup 只接受原 staging identity；不會刪 destination或不明路徑。

### 6.4 Formal frozen 結果

| 指標 | 預期 | 實際 | Verdict |
|---|---:|---:|---|
| Files | 14 | 14 valid / 0 failed | PASS |
| Rows | 328,697 | 328,697 input / 328,697 output | PASS |
| Input `Strong` | 66,841 | 66,841 | PASS |
| `Strong_Bidirectional` | 59,910 | 59,910 | PASS |
| `ClusterFirstOnly` | 6,931 | 6,931 | PASS |
| Split conservation | 66,841 | 59,910 + 6,931 = 66,841 | PASS |
| Unmapped/conflicts | 0 | 0 | PASS |
| HCC1395 input Strong | 9,228 | 9,228 | PASS |
| HCC1395 ClusterFirstOnly | 1,516 | 1,516 | PASS |
| Raw token preservation | PASS | PASS | PASS |
| Significant stable-key invariant | PASS | PASS | PASS |
| Stable-key uniqueness | PASS | PASS | PASS |
| Staging leftovers | 0 | 0 | PASS |

`migration_file_report.tsv` 有 14 data rows；`unmapped_conflicts.tsv` 與 `unresolved_files.tsv` 只有 header。

## 7. 驗證矩陣與實際命令

### 7.1 Baseline → post-change

| Check | Baseline | Post-change | Delta |
|---|---:|---:|---:|
| CMake build | exit 0 | exit 0 | 無 regression |
| GoogleTest | 234/234 | 251/251 | +17，0 failure |
| CTest | 234/234 | 251/251 | +17，0 failure |
| Python bounded contract/migration suites | n/a | 107/107 | 0 failure |
| Migration hardening suite | n/a | 24/24 | 0 failure |
| Formal frozen corpus | n/a | 14/14, 328,697 rows | VALID |
| Independent P2 validator | n/a | `VALID / 14 files / 328,697 rows` | PASS（exit 0） |

17 個新增 C++ cases = schema v2 4 + RegionStratifier 8 + artifact lifecycle 4 + root-review temp cleanup 1。

### 7.2 C++ commands

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
cmake --build build -j2
./build/bin/run_tests
ctest --test-dir build --output-on-failure
```

實際尾段：

```text
[==========] 251 tests from 36 test suites ran.
[  PASSED  ] 251 tests.

100% tests passed, 0 tests failed out of 251
```

### 7.3 Python bounded suites

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest -q \
  scripts.tests.test_migrate_verification_schema_v2 \
  scripts.tests.test_b0_clean1_verification_consumers \
  scripts.tests.test_b1_guards1_verification_consumers \
  scripts.tests.test_remaining_loh_critical_consumers \
  scripts.tests.test_remaining_pass_historical_consumers \
  tests.test_verification_schema_consumers \
  tests.test_verification_schema_historical_consumers \
  tests.test_verification_schema_b0_evidence \
  tests.test_verification_schema_b0_clean_consumers \
  tests.test_verification_schema_loh_round3 \
  tests.test_verification_schema_b1_guards2 \
  tests.test_verification_provenance_b2_consumers
```

實際結果：`Ran 107 tests ... OK`。H1 warning 是測試預期的 explicit historical warning，不是 failure。

### 7.4 第一輪獨立 output read-back

Formal tool之外，主 agent 重新讀 `migrated_outputs_manifest.tsv` 並從 14 個 output bytes重算：

```text
RECOMPUTED_FILES 14
HASH_PASS 14 ROW_PASS 14 TOTAL_ROWS 328697
ALL_PASS True
REPORT_ROWS 14 VALID 14 UNMAPPED 0
HCC_BEFORE_STRONG 9228
HCC_CLUSTER_FIRST_ONLY 1516
```

第二支獨立 validator 不 import migrator，另從 input/output逐列重驗；正式結果為：

```text
validator status=VALID files=14 rows=328697 input_strong=66841 strong_bidirectional=59910 cluster_first_only=6931 unmapped=0
exit code: 0
```

它的七個 synthetic negative-path tests亦為7/7 PASS，涵蓋缺 report、hash/count錯誤、duplicate stable key、Significant drift與既有 raw-token drift。

Machine-readable receipt：`InterSubMod/docs/reports/validated/2026/07/20260714_VerificationSchemaV2與RegionStratifier修復完整流程說明_01/qa/independent_p2_validator_receipt.json`。

## 8. 修改面、相容性與操作說明

### 8.1 Core 與 tests

| 類別 | 主要檔案 | 作用 |
|---|---|---|
| Build | `InterSubMod/CMakeLists.txt` | 註冊三支 schema/stratifier tests |
| C++ API | `InterSubMod/include/core/RegionProcessor.hpp`、`InterSubMod/include/core/SubcloneAnalyzer.hpp` | typed schema、status、assignment、aliases |
| C++ implementation | `InterSubMod/src/core/RegionProcessor.cpp`、`InterSubMod/src/core/SubcloneAnalyzer.cpp` | decision tree、LOH、stratification、publication |
| C++ tests | `InterSubMod/tests/test_verification_schema_v2.cpp`、`InterSubMod/tests/test_region_stratifier.cpp`、`InterSubMod/tests/test_region_stratification_artifacts.cpp` | truth table、boundaries、failure paths |
| Shared Python | `InterSubMod/scripts/lib/verification_schema_contract.py` | C2/L4/E/R1/P/LOH-L/H1 fail-closed loader |
| Migration | `InterSubMod/scripts/migration/migrate_verification_schema_v2.py` | frozen byte-preserving migration |
| Independent validation | `InterSubMod/scripts/validation/validate_verification_schema_v2_migration.py` | 唯讀重驗 formal output |

所有 74 consumer 的逐檔列表不在此重複貼 74 行，以 source map與 ledger為準；這避免報告摘要和 machine ledger各自漂移。

### 8.2 一個 migration release 的 aliases

| Deprecated surface | Canonical source | Exactness requirement | 移除 gate |
|---|---|---|---|
| `VerificationClass_V1_Deprecated` | v1 view | split leaf回 `Strong`，其他同 v2 | 外部 consumer inventory完成 |
| `LOH_Subtype` | `LOH_Subtype_LegacyVC` | 每列 exact equality | 下一 major schema |
| `Subclone_ID` | `RegionStratum_ID` | 每列 exact equality | region consumer migration完成 |
| `n_subclones` | `n_occupied_region_strata` | corrected exact mirror | source/ABI consumer gate |
| `subclone_structure.txt` | `region_stratification_summary.txt` | 只保留 deprecation stub | consumer不再讀舊檔 |
| C++ `SubcloneAnalyzer` name | RegionStratifier semantics | source compatibility only | 明確 release/ABI plan |

Repo 內 call sites已盤點；是否存在 repo 外 precompiled binary consumer無法由本 workspace直接證明。Header layout有新增欄，外部 binary若存在必須重新編譯。這是唯一 L5 open question，不影響 repo內 source/build驗收。

### 8.3 Production consumer 最小讀取流程

1. 先驗 schema/status/provenance；production 禁止 unversioned fallback。
2. 若要 current taxonomy，讀 C2；要 historical four-state，明確讀 L4。
3. 若問題是 support path，讀 E，不由 class名稱反推 boolean。
4. Region artifact 必須先接受 R1 terminal `VALID`；非 VALID不得讀 assignment當有效資料。
5. LOH 使用 `LOH_Subtype_LegacyVC`；若 alias存在，逐列驗 equality。
6. Pass-through 必須保留 P 的 version/current/legacy/evidence/derivation；不完整 provenance fail closed。
7. truth、clone、CCF/tree 等結論必須有各自明示來源；不能借用 Verification class或 RegionStratum名稱。

## 9. Definition of Done 對照

| 規格條件 | 結果 | 證據 |
|---|---|---|
| 保留 dirty worktree | PASS | 啟動時855 status lines；未 reset/刪除/覆寫無關檔 |
| canonical v2 兩個 split class + evidence | PASS | truth table + C++ tests |
| Significant/retention不變 | PASS | frozen逐 stable key raw-token audit |
| current/legacy分開計數 | PASS | fixed arrays、unknown buckets、conservation tests |
| occupied distinct IDs | PASS | sparse `{3}`、`{0,3}` tests |
| canonical不稱 cellular subclones | PASS | RegionStratifier naming + deprecated stub |
| 三 artifacts + 四 status | PASS | lifecycle/artifact tests |
| LOH canonical + exact alias | PASS | C++/Python contract + frozen audit |
| eligibility/write-back單一來源 | PASS | shared predicate、index identity、validate-all |
| tumor-only/insufficient/failed不留 stale valid | PASS | overwrite tests與 status protocol |
| 完整 consumer inventory | PASS | 74/74 migrated、0 unresolved |
| build/unit/integration/golden/hygiene | PASS | 251 C++、107 Python、24 migration + diff check |
| frozen 14 hashes/rows | PASS | preflight + read-back 14/14 |
| frozen locked counts | PASS | 59,910 / 6,931 / 0 conflicts |
| 最終 changed files/commands/risks | PASS | 本報告 §§7–10 |

## 10. 明確限制、風險與不應聲稱的結論

### 10.1 已接受的工程限制

- **跨檔 transaction**：owned files各自 atomic；整輪有效性依賴 status-first/last protocol。跳過 status的 consumer不受保證。
- **Same-UID hostile race**：P2 cleanup驗 EUID/path/dev/inode/parent identity，但不是 dirfd-based recursive deletion。受控 single-user staging threat model可接受；若威脅模型改變需再 harden。
- **C++ temp cleanup observability**：`TempFileGuard` 的 `std::remove()` 是 destructor best-effort且不回傳錯誤；forced rename test覆蓋正常權限下的cleanup，外部權限／namespace競爭不在本次 threat model。
- **External ABI**：repo外 precompiled C++ client未知；source compatibility不等於 binary ABI compatibility。
- **Compatibility deadline**：aliases設計為一個 migration release，但目前沒有機器可讀 removal version；應建立後續 issue，不可無限延長。

### 10.2 科學 claim ceiling

本任務**沒有**：

- 重跑 BAM、read matrices、distance、clustering、PERMANOVA或其他 statistics；
- 重算 F1/AUC或 caller benchmark；
- 證明 v1 upstream classification具生物學正確性；
- 證明 `ClusterFirstOnly` 是 subclone；
- 證明 `RegionStratum_ID` 是 cellular population、CCF cluster或演化樹節點；
- 驗證外部 caller/tool行為。

因此唯一可下的終局 claim 是：**在權威 frozen 14-file corpus與 repo code/test範圍內，schema mapping、serialization、provenance、consumer contract與 fail-closed publication符合鎖定規格。**

### 10.3 為何本次沒有啟動新 BAM 全量 run

權威規格明定先用 frozen corpus驗證，且不需要為此重跑 BAM/read matrices/distance/clustering/statistics。本次 formal「全量」指 manifest中的14個完整 frozen outputs，不是新的 BAM pipeline generation。若要啟動新的 production BAM run，應另立 run spec，鎖定 input BAM/VCF/reference hashes、validator版本、資源與成功條件；不得把本 migration output當作新 BAM run receipt。

## 11. 可追溯來源

### 11.1 權威規格與 provenance

- `/big7_disk/liaoyoyo2001/external_validation/subclone_reconstruction_method_comparison_20260714/20260714_InterSubMod分類與RegionStratifier修復任務規格_01.md`；SHA-256 `e713ad9dd0d716ff067e9523954dc9f429903730e4821ba24a37b676c89e38cd`。
- `/big7_disk/liaoyoyo2001/external_validation/subclone_reconstruction_method_comparison_20260714/FROZEN_CORPUS_MANIFEST.tsv`；SHA-256 `9ab4428e28d137b3c15e811563678d006e6258ff6d3f86e718cc4f56ead2e0ea`。
- `InterSubMod/scripts/migration/migrate_verification_schema_v2.py`；formal-run SHA-256 `589c52605e96d07af6bbfa6810dfdea1ca04c3875052c21e2a7c90100524d3cd`。

### 11.2 實作與 consumer source maps

- `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md`
- `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/20260714_CPP分類與RegionStratifier流程來源對照_01.md`
- `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/20260714_PythonConsumer契約與遷移對照_01.md`
- `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_inventory.md`
- `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv`

### 11.3 Formal output receipts

- `migration_summary.json` SHA-256 `454a259ec770cf1d8df7bc4f0194f9461e160efc2c88130e03e80a5c7d65532f`
- `migration_file_report.tsv` SHA-256 `c3ca00d8bd54af9e13856d625edf32040a0e07949dc6a10e2d0be19b924a3370`
- `migrated_outputs_manifest.tsv` SHA-256 `62433a28492b82783ed1c936a3aad9aeb1b30dbdd4a809ab2a1f3826d1d63c04`
- `migration_status.tsv`：`VALID / ALL_FILES_MIGRATED / 14 / 14 / 0 / 328697 / 328697 / 0 / PASS / PASS / PASS`。

## 12. 快速交接

若下一個 reviewer只做三件事：

1. 讀 `migration_status.tsv` 與 `migration_summary.json`，確認 status、14/14、328,697、三項 PASS。
2. 執行 independent validator，確認 output hash、stable key、Significant與 raw tokens全量重驗通過。
3. 讀兩份 C++/Python source map，確認任何新 consumer選擇的是 C2/L4/E/R1/P/LOH-L/H1之一，而不是裸字串推測。

## 13. Standalone HTML 瀏覽器 QA

Chromium 147 headless實測三次：第一次在390px抓到 mobile body overflow，第二次定位為single-column grid的639px min-content width；修成 `minmax(0, 1fr)` 並對 grid items加 `min-width:0` 後，第三次 PASS。

- Desktop 1440×1000：無 body overflow；11/11錨點；9/9表格有 wrapper；SVG 1041.41×483.50。
- Mobile 390×844：無 body overflow；11/11錨點；9/9表格可在內層捲動；SVG依比例縮為336.81×156.38。
- 兩者皆為0 console errors、0 page errors、0 external requests、0 placeholders；目視檢查 PASS。
- QA receipt：`InterSubMod/docs/reports/validated/2026/07/20260714_VerificationSchemaV2與RegionStratifier修復完整流程說明_01/qa/browser_qa_receipt.json`。

![Standalone HTML desktop QA](qa/desktop.png)

![Standalone HTML mobile QA](qa/mobile.png)

---

**Scope footer**：Comprehensive engineering validation；frozen 14-file full scope；非 BAM rerun、非生物學 truth、非 clone reconstruction validation。報告狀態：`VALIDATED`。
