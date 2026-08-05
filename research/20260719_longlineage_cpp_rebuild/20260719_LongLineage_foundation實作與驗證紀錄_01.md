<!--
建立時間: 2026-07-19
目標: 紀錄 LongLineage 新專案 foundation 實作、AI治理、資料契約與驗證證據
處理範圍: P0-P2/P6 foundation；P3-P8維持blocked
關聯檔案:
  - InterSubMod/research/20260719_longlineage_cpp_rebuild/00_INDEX.md
  - InterSubMod/research/20260719_longlineage_cpp_rebuild/implementation-notes.md
  - /big7_disk/liaoyoyo2001/LongLineage/docs/audits/20260719_foundation_verification.md
-->

# LongLineage foundation 實作與驗證紀錄

> **TL;DR**：新專案已建立可執行的AI治理、資料紀錄／格式／查詢契約及C++
> foundation；Debug、Release、ASan/UBSan各25/25 PASS，但strict release gate
> 正確地以12個缺證據項目阻擋，因此不宣稱P0-P8完成。

## 服務目標與scope

- Task type：B — Comprehensive validation。
- 服務G2/G3/G4/G5：把latest HP/PS、read-level methylation、sSNV
  co-occurrence與topology authority拆成可驗證C++ trust domains。
- 本輪只驗證repository foundation；沒有啟動真實資料長計算、沒有讀取truth、
  沒有執行Python科學程式。

## 輸入、命令、輸出

| 類別 | 輸入路徑 | 完整命令 | 輸出路徑／片段 | Exit |
|---|---|---|---|---:|
| Debug tests | `/tmp/LongLineage/build-verify-debug` | `ctest --test-dir build-verify-debug --output-on-failure` | `/tmp/LongLineage/build-verify-debug/Testing/Temporary/LastTest.log`；`25/25 PASS` | 0 |
| Release tests | `/tmp/LongLineage/build-verify-release` | `ctest --test-dir build-verify-release --output-on-failure` | `/tmp/LongLineage/build-verify-release/Testing/Temporary/LastTest.log`；`25/25 PASS` | 0 |
| Sanitizer tests | `/tmp/LongLineage/build-verify-sanitize` | `ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 ctest --test-dir build-verify-sanitize --output-on-failure` | `/tmp/LongLineage/build-verify-sanitize/Testing/Temporary/LastTest.log`；`25/25 PASS` | 0 |
| Debug foundation | `/tmp/LongLineage`＋Debug binaries | `scripts/ci/check_all.sh build-verify-debug` | `FOUNDATION_PASS`、`NO-NETWORK PASS` | 0 |
| Release foundation | `/tmp/LongLineage`＋Release binaries | `scripts/ci/check_all.sh build-verify-release` | `FOUNDATION_PASS`、`NO-NETWORK PASS` | 0 |
| Strict evidence | gate registry＋Release CTest | `scripts/ci/check_gate_test_coverage.sh build-verify-release --strict` | `FAIL missing=12`（預期fail closed） | 1 |
| AI readiness | `/tmp/LongLineage`＋Debug governance | `scripts/ai/check_readiness.sh build-verify-debug` | `PASS failures=0 warnings=1`；warning=Git尚未初始化 | 0 |

## AI工作環境與程式開發紀律

1. `AGENTS.md`定義cold-start Q1-Q5、task classification、goals、Step→Verify、
   truth isolation、Python boundary與禁止事項。
2. `.ai/README.md`與task template要求owner、parent/dependency DAG、typed
   write-set、lease、heartbeat、inputs、outputs、scope及驗證。
3. C++ governance負責驗證task registry、phase ledger、schema/catalog、status、
   gates與truth boundary；Markdown不能自行升級phase。
4. 開發流程要求short-lived branch、format/build/unit/negative/sanitizer、
   independent validation、review及clean merge；`main`是唯一永久branch。
5. `state/release_attestation.json=NOT_READY`；production kernel缺任何科學gate即
   非零退出，不能先做PASS receipt。

## 資料紀錄、格式與查詢

- `schema/id_registry.json`：22個offline schema ID＋physical SHA lock。
- `schema/catalog.json`：16種artifact的record schema、key、order與index contract。
- `contracts/v1/*.tsv`：87個status/reason與type/lifecycle/transform/query operator
  closed vocabulary。
- `provenance/source_to_target_manifest.json`：12項port逐一分離
  target presence、implementation與verification，不把檔案存在當作parity。
- nested JSON、Interval0 pair、lineage、membership、semantic digest與receipt chain
  皆有machine contract與positive/negative fixtures。
- Query只讀`VALIDATED_FROZEN` artifact；P6 query engine完成前fail closed。
- 查詢入口：
  - `/big7_disk/liaoyoyo2001/LongLineage/docs/data/RECORD_AND_QUERY_STANDARD.zh-TW.md`
  - `/big7_disk/liaoyoyo2001/LongLineage/docs/data/DATA_CONTRACTS.md`
  - `/big7_disk/liaoyoyo2001/LongLineage/docs/data/QUERY_GUIDE.md`

## 未完成與NO-GO邊界

- GitHub private remote尚未建立／驗證。
- Draft 2020-12 validator尚未pin；目前`jsonschema 3.2.0`只能算exercised
  compatible fixtures，不算宣告draft完整驗證。
- strict gate仍有12個fixture-only或缺executable fault-injection evidence項目。
- P2缺block reader、byte-bounded reorder sink與1/2/4/24/40 worker parity。
- P3-P5缺M1 RNG/logical vectors、formal full co-occurrence authority、完整topology
  parity與direct HiGHS。
- P6-P8缺獨立artifact replay、7-dataset全量雙跑與bilingual release QA。

## 證據鏈

新repo先在`/tmp/LongLineage`建立，以避免在未初始化Git前把staging誤認正式專案。
完成clean commit與immutable audit envelope後才搬入
`/big7_disk/liaoyoyo2001/LongLineage/`；GitHub remote若缺credential則保留
machine blocker，不會宣稱已完成。

### 正式本機路徑closeout

- 正式路徑：`/big7_disk/liaoyoyo2001/LongLineage/`
- Branch／HEAD：`main` /
  `3a789d3c8b384606dfad01ae0227834df01661ff`
- Git狀態：clean；199 tracked files；未搬入任何`build*`或`.cache`。
- Current machine audit：
  `state/audits/20260719-foundation-verification-002.json`
- Current audited source commit：
  `8b62261a384bd2dd2a469f5b2ad27df2e34f3c8d`
- Current audited tree SHA-256：
  `eb59c1f1856692569729742378d00ca396ad3a1ed125bc4aa0b395903a155bd3`
- `002` forward-supersedes `001`；兩者在正式路徑皆由
  `scripts/ci/check_audit_source_snapshot.sh`重播PASS。

正式路徑fresh smoke：

| Step | 完整命令 | 輸出路徑／片段 | Exit |
|---|---|---|---:|
| Configure | `/usr/bin/cmake -S /big7_disk/liaoyoyo2001/LongLineage -B /tmp/LongLineage-target-smoke -DCMAKE_BUILD_TYPE=Release -DLONGLINEAGE_WARNINGS_AS_ERRORS=ON` | `/tmp/LongLineage-target-smoke`；HTSlib 1.18、Jansson 2.13.1、OpenSSL 3.0.2 | 0 |
| Build | `/usr/bin/cmake --build /tmp/LongLineage-target-smoke --parallel 4` | 六個CLI＋tests；`Built target longlineage_governance` | 0 |
| CTest | `/usr/bin/ctest --test-dir /tmp/LongLineage-target-smoke --output-on-failure` | `100% tests passed, 0 tests failed out of 25` | 0 |
| Foundation gate | `scripts/ci/check_all.sh /tmp/LongLineage-target-smoke`（cwd為正式repo） | `FOUNDATION_PASS`、`NO-NETWORK RESULT: PASS` | 0 |

GitHub private remote仍未建立：本機GitHub credential無效且目前connector沒有
create-repository mutation；因此P0的`PRIVATE_REMOTE_NOT_YET_VERIFIED`保留。
