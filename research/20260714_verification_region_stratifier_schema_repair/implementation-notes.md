<!--
建立時間: 2026-07-14 22:28 +08:00
目標: Verification schema v2 與 RegionStratifier 修復實作過程 living document
處理範圍: verification-region-stratifier-schema-repair
cycle_id: cycle_20260714-2216-verification-region-stratifier-schema-repair
spec_id: 20260714-verification-region-stratifier-schema-repair
status: in_progress
advisory: on
關聯檔案:
  - /big7_disk/liaoyoyo2001/external_validation/subclone_reconstruction_method_comparison_20260714/20260714_InterSubMod分類與RegionStratifier修復任務規格_01.md
  - InterSubMod/research/20260714_verification_region_stratifier_schema_repair/pre-decision-audit.md
  - InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_inventory.md
  - InterSubMod/.claude/skills/implementation-notes/SKILL.md
-->

# Implementation Notes: Verification schema v2 與 RegionStratifier 修復

> **Purpose**: 記錄實作期間的 protected decisions、相容性折衷、偏離與未決風險。
> **Current gate**: `COMPLETE_VALIDATED`；任一後續 Significant/raw-cell/predicate 變動必須重開 frozen gate。
> **Task type**: E Hotfix + B Comprehensive validation；服務 G4/G5。

## Frontmatter

- **Spec source**: `/big7_disk/liaoyoyo2001/external_validation/subclone_reconstruction_method_comparison_20260714/20260714_InterSubMod分類與RegionStratifier修復任務規格_01.md`
- **AI session**: Codex 2026-07-14
- **Last updated**: 2026-07-15（P0→P1→P2 完成，formal + independent validation PASS）
- **Pre-decision verdict**: `GO_WITH_FAIL_CLOSED`, 70/100

---

## 🔵 設計決定（Design Decisions）

### [2026-07-14 22:28] Machine schema 完全採 handed spec

- **Status**: Accepted / Protected
- **背景**: 較早方法比較文件曾使用 display phrase `ClusterOnlyAssociation`。
- **決定**: We will use only `Strong_Bidirectional` / `ClusterFirstOnly` and the exact 11-class/evidence enums from the handed repair spec.
- **理由**: spec §581 明定這些字串是 interoperability contract；display wording 不可回寫 machine data。
- **影響範圍**: C++ enum/string serializer、CSV、Python helper、migration、tests。
- **Revisit if**: 只有使用者發布新的 major schema spec。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

<!-- BEGIN USER-SPECIFIED -->
<!-- Decision: Verification v2 11 classes、Stage④ precedence、RegionStratum IDs/thresholds/status、frozen counts皆不可自行改字串或門檻。 -->
<!-- DO NOT change: 這些是 2026-07-14 handed spec 的鎖定互通契約。 -->
<!-- END USER-SPECIFIED -->

### [2026-07-14 22:28] RegionProcessor 擁有 stable-key artifact serialization

- **Status**: Accepted
- **背景**: `RegionAsmProfile.chrom/pos` 目前未初始化且沒有 Ref/Alt；`RegionProcessor` 才持有 `snvs_` 與 `chrom_index_`。
- **決定**: We will keep stratification logic in `SubcloneAnalyzer` compatibility surface, return typed `{result_index, region_id, stratum}`, and let `RegionProcessor` serialize `RegionID,Chr,Pos,Ref,Alt` after validate-all.
- **理由**: 避免把 metadata 複製進 profile 造成第二個 identity SoT；可在寫回前驗 index/region ID。
- **影響範圍**: `RegionProcessor.*`, `SubcloneAnalyzer.*`, artifact tests。
- **Revisit if**: 未來 RegionResult 自身成為完整 locus-owned record。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-07-14 22:28] Status 是跨 artifact commit marker

- **Status**: Accepted
- **背景**: per-file atomic rename 仍可能在 process crash 後留下混代檔。
- **決定**: We will atomically publish `status=FAILED, reason=RUN_IN_PROGRESS` before data publication, write/close/check all owned artifacts and the main CSV, then atomically publish the terminal status last.
- **理由**: 在固定四狀態契約內達到 single-writer crash fail-closed；任何中斷都不呈現 VALID。
- **影響範圍**: status writer、canonical TSV/summary/stub、significance CSV、reader contract。
- **Revisit if**: 需要 concurrent lock-free snapshot 或 cryptographic cross-file binding。
- **Evidence tier**: L2 ⭐⭐⭐⭐

### [2026-07-14 22:28] 保留 C++ legacy 命名作 source compatibility

- **Status**: Accepted
- **背景**: 正式 rename class/file 不是 correctness 必要條件，且 public headers 可能有 repo 外客戶。
- **決定**: We will keep `SubcloneAnalyzer`/deprecated mirrors for one migration release, but canonical methods, comments, artifacts and caller wording will use region stratification semantics only.
- **理由**: 最小化無關 build/ABI 破壞；canonical output 不再含 biological clone claim。
- **影響範圍**: public headers、logging、deprecated aliases/stub。
- **Revisit if**: 下一 major schema/API release。
- **Evidence tier**: L2 ⭐⭐⭐⭐

### [2026-07-14 22:28] Python consumer 使用單一 contract helper

- **Status**: Accepted
- **背景**: bounded inventory 有 62 支 VerificationClass、23 支 LOH readers，只有 2 支讀 legacy。
- **決定**: We will centralize C2/L4/E/R1/LOH-L/P/H1 validation and migrate consumers by behavior class; unknown never falls through to Noise.
- **理由**: 避免 62 套 enum/fallback 漂移，並讓 historical mode 明確留痕。
- **影響範圍**: shared Python helper、四支鎖定 consumer、B0 selectors、B1/B2 guards/provenance。
- **Revisit if**: schema package獨立發版。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-07-14 22:28] 固定 status reason 與時間格式

- **Status**: Accepted
- **決定**: We will use `OK`, `VALID_WITH_WARNINGS`, `NO_ELIGIBLE_REGIONS`, `BELOW_MIN_REGIONS`, `NORMAL_BAM_REQUIRED`, `RUN_IN_PROGRESS`, and machine-specific `*_FAILED` reason codes; `generated_at` uses RFC3339 UTC seconds (`YYYY-MM-DDTHH:MM:SSZ`).
- **理由**: status 欄固定但 reason 尚未列 literal；封閉字彙與 UTC 格式可測且跨時區。
- **影響範圍**: status artifact、tests、說明文件。
- **Revisit if**: status schema v2 引入 reason enum/generation ID。
- **Evidence tier**: L3 ⭐⭐⭐

### [2026-07-14] 非法 legacy provenance 在 live 邊界停止

- **Status**: Accepted / Root-review fix
- **背景**: 初版 `classify_verification_v2()` 與 LOH helper 對未知 legacy class 仍可能繼續分類或回傳 `None`。
- **決定**: legacy class 只接受 `Strong|Subclone|Weak|Noise`；其他值在 Verification 與 LOH helper 入口皆拋出 `INVALID_VERIFICATION_CLASS_LEGACY`。
- **理由**: 固定 enum 是 provenance contract；未知值不得被誤解為無 LOH 或 Noise。
- **驗證**: 新增 classifier/LOH unknown regression，targeted C++ 16/16 PASS。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-07-14] Artifact writer 自行驗 status、row、slot 與 stable key

- **Status**: Accepted / Root-review fix
- **背景**: 初版 writer 信任 caller，理論上可收到 `VALID` 但 assignment count、occupied count、四槽 summary 或 row enum 互相矛盾的資料。
- **決定**: 在任何檔案建立前驗 schema/count/status invariants、四個固定 slots、逐 stratum row count、result index 與 stable key uniqueness、label/reason exact mapping。
- **理由**: fail-closed 必須位於 publication boundary，而不能只依賴上游目前恰好正確。
- **驗證**: contradictory fixture 被拒且 output directory 未建立；完整 valid/stale fixtures PASS。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-07-14] P2 post-write audit 驗完整 appended provenance

- **Status**: Accepted / Root-review fix
- **決定**: byte preservation audit 除 V1 raw token 與 LOH raw alias外，也逐列核對 schema version、兩個 support booleans、兩個 `NA`、EvidencePath、EvidenceDerivation。
- **理由**: 只保證 raw cell 未改仍不足以證明新增 provenance 正確。
- **驗證**: table-driven tamper test 對七個欄位逐一竄改，全部 fail closed；migration tests 18/18 PASS。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-07-14] 三支 untracked historical consumer 採隔離式局部遷移

- **Status**: Accepted
- **背景**: `build_hcc1395_html.py`、`hcc1395_full_accounting.py`、`tp_fp_structure_label_association.py` 是任務啟動前即存在的 untracked user artifacts，但屬 B0 inventory。
- **決定**: 先記錄原 SHA-256，只局部加入 schema/status/legacy/LOH guard；不改資料路徑、門檻與統計方法。HTML 不再由 class 名稱重建 truth；association 明確使用 legacy view與 canonical within-HP row field。
- **驗證**: 三檔 py_compile、兩個 CLI help、historical selector regression 2/2 PASS；原始與修改後 SHA 均留存。
- **Evidence tier**: L2 ⭐⭐⭐⭐

### [2026-07-14] B0 consumer 全部以用途選擇 C2／L4／E，不從 class 名稱猜證據

- **Status**: Accepted / P1-B
- **決定**: live current panel 使用 `C2` 並驗 `VerificationSchemaVersion=2`；歷史 replication 使用 `L4`，unversioned v1 只在明示 H1 flag 下開啟；需要 label/cluster support 的流程直接讀 typed `E` 欄。truth label 不再由 VerificationClass 映射。
- **理由**: v2 把 evidence path 展開為 11 類後，四態字串不再同時代表 current taxonomy、historical cohort 與 evidence boolean。
- **驗證**: B0-Clean-1 targeted 12/12、B0-Clean-2 與 shared contract 20/20、evidence/shared group 21/21、historical group 2/2；主代理另做禁用 selector 精準掃描與 B0-Clean-2 20/20 重跑。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-07-14] Artifact publication 失敗也必須清除上一 run 的有效外觀

- **Status**: Accepted / Root-review fix
- **背景**: publication-boundary validation 可能在任何 canonical artifact 寫入前拒絕矛盾資料；只把 commit marker 改為 `FAILED` 仍會讓上一 run 的 assignment/summary 檔留在原位。
- **決定**: 主 artifact write 失敗時，以四個空 fixed slots 與 `FAILED/ARTIFACT_WRITE_FAILED` 再執行一次 owned-artifact 覆寫，最後寫 FAILED status；原始錯誤與 fallback 失敗皆保留在 exception。
- **理由**: spec P1-C 明定非 VALID 的 assignment/occupied 必須為 0，owned artifacts 不得沿用上一 run。
- **驗證**: incremental build exit 0；Verification v2、RegionStratifier、artifact lifecycle targeted 16/16 PASS。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-07-14] Consumer migration 以 inventory 為封閉母體完成 74/74

- **Status**: Accepted / P1 complete
- **決定**: 以既有 62 個 VerificationClass paths 與 23 個 LOH paths 合併成 74 個唯一 executable consumers；每一列保留 inventory contract、blocking tier、disposition、test/guard 與備註，不以新的全 repo 大掃描取代已審核清單。
- **結果**: `consumer_migration_status.tsv` 為 74 unique paths、74 `MIGRATED`、0 unresolved；production v2 default fail-closed，歷史 H1 只可 explicit opt-in。
- **驗證**: B0/B1/B2/LOH 分批 targeted regressions與跨批 suites全綠；X5 synthetic H1 default reject、explicit flag pass。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-07-14] P2 final publish 使用 kernel-level no-replace

- **Status**: Accepted / Red-team blocker closed
- **背景**: 初版先 `output_dir.exists()` 再 `os.rename(staging, output_dir)`；目的目錄若在兩者間出現，POSIX rename 可取代剛出現的空目錄。
- **決定**: 正式 migration 維持 NO-GO，改用 Linux atomic no-replace primitive；`EEXIST` 必須保留目的 sentinel/inode/content。任何不支援 no-replace 的 kernel/filesystem都明確拒絕，不得退回 `os.rename`。
- **附帶修正**: 例外時只清理由本程序 `mkdtemp` 建立且 identity/prefix/parent 全相符的 staging；不得碰 destination 或其他目錄。report 加 raw-token、Significant/stable-key invariance 明示 flags，時間改秒級 `Z`。
- **Gate result**: 24/24 migration unit tests、6/6獨立 hardening review、實際 target filesystem probe與 formal 14-file publish 全部通過；正式輸出在 staging 驗證完成前維持 absent，最後以 no-replace 原子發布。
- **Revisit if**: kernel/filesystem、publish primitive或 same-UID hostile-race threat model改變。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-07-14] Post-change C++ 完整驗證

- **Status**: Accepted / Verified
- **結果**: build exit 0；GoogleTest 250/250；CTest 250/250。相較 baseline 234/234 新增 16 cases，沒有新增失敗。
- **新增覆蓋**: VerificationSchemaV2 4、RegionStratifier 8、RegionStratificationArtifacts 4。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-07-14] Atomic temp file 失敗路徑改用 owned RAII cleanup

- **Status**: Accepted / Root-review hardening
- **背景**: source-map review 指出 `write_atomic_text()` 與 significance CSV/statistics 在 writer、flush、close 或 rename 失敗時，final path 雖不會被誤發布，但本輪 temporary file 可能殘留。
- **決定**: 每一個本輪建立的 temp path 都綁定 `TempFileGuard`；只有成功 rename 後才 `mark_published()`，其他 return/exception 路徑由 guard 移除該 temp，不觸碰 final path 或其他檔案。
- **驗證**: `cmake --build build -j2` exit 0；新增 rename-collision failure regression 1/1 PASS；`./build/bin/run_tests` 251/251；`ctest --test-dir build --output-on-failure` 251/251。
- **限制**: destructor內的 `std::remove()` 是 best-effort且不回報自身失敗；正常權限failure path已測，外部同時改動目錄權限／namespace屬 hostile-race threat model，terminal status仍維持fail-closed。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-07-15] P2 frozen 14-file formal migration 原子發布完成

- **Status**: Accepted / Formal output VALID
- **唯一執行命令**: `python3 scripts/migration/migrate_verification_schema_v2.py --manifest /big7_disk/liaoyoyo2001/external_validation/subclone_reconstruction_method_comparison_20260714/FROZEN_CORPUS_MANIFEST.tsv --input-root /big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/20260620_allsample_subcluster_split/results --output-dir output/schema_migrations/20260714_verification_v2 --require-output-dir-absent`
- **輸出**: `/big7_disk/liaoyoyo2001/big7_disk_output/schema_migrations/20260714_verification_v2`；publish 前 authoritative path absent，publish 後 staging leftovers=0。
- **結果**: status `VALID/ALL_FILES_MIGRATED`；14/14 valid、0 failed；328,697→328,697 rows；unmapped=0；raw-token preservation、Significant stable-key invariant、stable-key uniqueness皆 PASS。
- **鎖定計數**: input Strong=66,841；Strong_Bidirectional=59,910；ClusterFirstOnly=6,931；兩者合計66,841；HCC1395 input Strong=9,228、ClusterFirstOnly=1,516。
- **第一輪獨立 read-back**: 14/14 output SHA-256與row counts重算符合 `migrated_outputs_manifest.tsv`；14/14 file reports為 VALID；unmapped總和0；total rows 328,697。
- **第二輪獨立 byte-level validator**: 不 import migrator，以 read-only FD + mmap重驗14組input/output與七份 provenance reports；exit 0，輸出 `VALID files=14 rows=328697 input_strong=66841 strong_bidirectional=59910 cluster_first_only=6931 unmapped=0`。其七個 fail-closed synthetic tests為7/7 PASS。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

---

## 🟠 偏離之處（Deviations）

### [2026-07-14 22:28] 尚無已接受偏離

- **Status**: Closed
- **規範要求**: 依 P0→P1→P2 全量實作。
- **實作偏離**: 無。
- **風險評估**: 若後續發現無法安全遷移的 undocumented external interface，立即重開並依 spec §9 停止。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

---

## 🟡 折衷考量（Trade-offs）

### [2026-07-14 22:28] Status-last commit marker vs generation manifest

- **Status**: Accepted
- **方案 A**: We will use FAILED/RUN_IN_PROGRESS first and terminal status last under a single-writer output-directory contract.
- **方案 B**: 新增 generation_id + per-file hash manifest；更強，但會擴充 handed fixed-column schema。
- **方案 C**: 目錄世代 + symlink swap；跨多檔近似交易，但改變 canonical 路徑與現有 run layout。
- **採用 A 理由**: 能滿足本 spec 的 fail-closed/stale-artifact contract，且不自行升 schema。
- **Tier check**: process-crash 可否證；concurrent-reader TOCTOU 明列為限制。
- **Revisit if**: 同一 output directory 允許 concurrent writer/reader，或需要 immutable cryptographic provenance。
- **Evidence tier**: L2 ⭐⭐⭐⭐

### [2026-07-14 22:28] Frozen migration vs BAM/statistical rerun

- **Status**: Accepted / Protected
- **方案 A**: We will migrate only recorded v1 current+legacy provenance and verify byte preservation/Significant invariance.
- **方案 B**: 重跑 BAM/read matrix/clustering/statistics；成本高且超出 handed scope。
- **採用 A 理由**: spec §412 明確不需 BAM rerun；C++ truth table獨立接住新 code path。
- **Revisit if**: predicate、p-value 或 biological validity 進入 scope。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

---

## 🔴 未決問題（Open Questions）

### [2026-07-14 22:28] Repo 外已編譯 C++ client

- **Status**: open
- **Question**: 是否有未受 repo inventory 覆蓋、以舊 public struct layout連結的外部 binary client？
- **Context**: header新增欄位會改 layout；repo本身只建 static library、無 install/SOVERSION。
- **Default if no answer**: 保留 source-level deprecated mirrors並在最終相容性章節要求外部 client重新編譯。
- **Revisit if**: 使用者提供外部 ABI consumer。
- **Priority**: minor
- **Evidence tier**: L5 ⭐

---

## 📚 Lore（Prior Gotchas / Non-obvious Constraints）

### [2026-07-14] Dirty worktree 邊界

- **Constraint**: repo 啟動時有 855 筆 status；既有 `tests/*.cpp` 多為 100644→100755 mode-only changes。
- **Why it matters**: 新測試必須另建檔，只修改 content-clean CMake；不得 reset/覆蓋使用者模式變更。
- **Evidence**: baseline targeted `git status/diff` + core audit SHA readback。

### [2026-07-14] Frozen corpus 的證據邊界

- **Constraint**: 14/14 hash與rows可驗 migration，不可獨立證明 v1 upstream classification/biology。
- **Why it matters**: 終版只能 claim schema/serialization regression PASS。
- **Evidence**: handed spec §374–412；pre-decision red-team probe D。

---

## Provenance Footer

- **Commit hash**: `6067568637088838a9f518955e41d222f057e4f1`
- **Build baseline**: build exit 0；run_tests 234/234；ctest 234/234
- **Skill version**: `/implementation-notes` v0.1
- **Active cycle**: `InterSubMod/state/cycles/cycle_20260714-2216-verification-region-stratifier-schema-repair/`
- **Pre-decision audit**: `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/pre-decision-audit.md`
