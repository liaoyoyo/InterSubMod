<!--
建立時間: 2026-08-13 11:30 +08:00
目標: 提供 InterSubMod／LongLineage research handoff snapshot 的單一入口
處理範圍: Task Type B Comprehensive validation + D External handoff；7 technical datasets／6 biological IDs／chr1-22
關聯檔案:
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv
  - InterSubMod/research/20260813_complete_research_handoff/release_freeze_manifest.json
驗證方式: authority SHA-256 replay、registry validator、clean build/test、synthetic E2E、reader acceptance test
證據等級: mixed；每一 claim 以 registry/文件內的 evidence status 為準
狀態: RESEARCH_HANDOFF_SNAPSHOT / RELEASE_BLOCKED
-->

# InterSubMod／LongLineage 完整研究交接索引

> **RESEARCH HANDOFF SNAPSHOT — 非 production release。** 本包整理的是截至 2026-08-13 的研究狀態、資料治理與可重現介面。InterSubMod 的 tag／GitHub Release，以及 LongLineage 的 public preview，必須等各自 gate 全數通過；目前不得稱為 release-ready。

用 SCQA + Claim–Evidence–Verdict：**現有 frozen science evidence 可重播且 19/19 bytes 相符，但公開文件、雙機、LongLineage license/source-origin 與 live publication gates 尚未全部閉合，所以可交接研究、不可發布 production。**

## 30 秒答案

| 問題 | 可安全回答 |
|---|---|
| 這是什麼專案？ | 以 ONT 長讀序列的 sSNV read linkage、候選 mutation-state reconstruction 與 read-level methylation association，研究癌症樣本內的分子狀態；不是已完成的 cellular lineage caller。 |
| 現在的科學結論？ | `confirmed cellular subclone = 0`、`confirmed linear ancestry = 0`。Methylation 僅為 association；CN/LOH 未整合進 frozen candidate reconstruction，未做 CN/LOH-corrected inference；InterSubMod optional LOH BED 只做 annotation／stratification。`88.2579%` 只是 frozen model 下 71,955 ranked units 的 rooted-unlabeled graph-shape 統計。 |
| 哪些數字可引用？ | 只從 2026-08-01 frozen authority manifest、denominator registry 及其 19/19 replay receipt引用，並同時保留 denominator、scope 與 claim boundary。 |
| 哪些資料是 final？ | 可引用的 science final 必須同時是 `AUTHORITY + FULL + FINAL_FOR_SCOPE` 且有 producer/hash；另有 1 筆 `VALIDATED_DERIVED` 只對 provenance adjudication final，不是 science final。Final 僅代表其 scope，並非整個研究 final。 |
| 軟體如何分工？ | LongPhase-S/TO 產生 phasing／HP/PS／recalibrated VCF／tagged reads；exact-PS／LongLineage 建候選 mutation-state families；InterSubMod 產 per-region methylation、read distance、read clustering與統計。Commit-pinned Python research solver 可是 science producer；publication builder／HTML 只呈現 validated data。 |
| bip7／bip8 怎麼跑？ | 先 bootstrap site profile，再跑 doctor、`run_benchmark.sh --plan`、clean build/test與 synthetic smoke。只有真正登入該主機取得的 receipt 才能宣稱該主機 PASS。 |
| 如何繼續研究？ | 先選 registry 中非 superseded artifact，建立新 cycle/pre-decision audit，固定 source commit、輸入 hash、scope、denominator與停止條件；active CNV/drilldown另走後續 PR。 |

## 科學權威與判讀順序

1. [2026-08-01 authority manifest（包內exact copy）](evidence/authority_manifest.json)
2. [2026-08-01 denominator registry（包內exact copy）](evidence/denominator_registry.tsv)
3. [本包 artifact registry](registries/artifact_registry.json)
4. [本包 authority／superseded crosswalk](registries/authority_superseded_crosswalk.json)
5. [2026-08-13 authority replay receipt（包內exact copy）](evidence/authority_replay_receipt.json)

禁止從目錄名稱、檔名中的 `final`、最新 mtime、目前 Git HEAD 或 HTML 視覺呈現自行推定權威。

## 文件地圖

| 讀者需求 | 文件／registry | 用途 |
|---|---|---|
| 現況、時間與結論 | [研究結論、時間與 finality](20260813_研究結論時間與Finality_01.md) | 每個核心 claim 的證據、限制與 verdict |
| 三套工具與資料流 | [軟體輸入輸出與研究流程](20260813_軟體輸入輸出與研究流程_01.md) | LongPhase、LongLineage/exact-PS、InterSubMod、Python/HTML 的責任邊界 |
| bip7／bip8 操作 | [雙機操作與驗證](20260813_bip7_bip8操作與驗證_01.md) | profile、doctor、plan、build/test、replay、fail-closed 規則 |
| AI cold start | [AI context bundle](ai_context/CONTEXT.md) | 無對話背景仍能正確回答六問；含禁止升格規則 |
| 資料集 | [dataset registry](registries/dataset_registry.json) | 7 technical datasets 與 6 biological IDs |
| Run | [run registry](registries/run_registry.json) | 18 logical rows合併至51 physical dirs（35 current＋16 pending archive）與取代關係 |
| Artifact | [artifact registry](registries/artifact_registry.json) | producer、commit、input、hash、scope、finality、限制與位置 |
| 公開 claims | [claim registry snapshot](registries/claim_registry.json) | hash-bind完整158-row registry、verdict分布與source/publication/release gates |
| Dataset alias | [dataset alias registry](registries/dataset_alias_registry.json) | site-profile operator key與canonical technical dataset ID的join規則 |
| 機器路徑 | [machine-path registry](registries/machine_path_registry.json) | bip7 本機路徑、遠端 mount view、missing source 與最後確認時間 |
| 中間資料根 | [storage-root manifest](registries/storage_root_manifest.json) | directory-level manifest、數量、大小量測狀態與重生方式 |
| Workflow／script治理 | [workflow registry](registries/workflow_registry.json) | Git tracked `scripts/**` exact inventory、portable support allowlist、legacy/archive界線、blob hash與absolute-path token盤點 |
| 大型tracked assets | [large-asset registry](registries/tracked_large_asset_registry.json) | >1 MiB逐檔SHA／大小／政策目標；不從副檔名推finality或license |
| 離線瀏覽 | [Standalone HTML overview](20260813_完整研究交接總覽_01.standalone.html) | 無外部 request 的單檔總覽；數字仍回指 registry |
| 實作紀錄 | [Implementation notes](implementation-notes.md) | 決策、偏離、折衷、驗證事件與發布阻塞 |
| 包內證據 | [Evidence manifest](evidence/EVIDENCE_MANIFEST.json) | authority、denominator、replay、preflight、claims、hygiene與LongLineage safety的SHA-256 |

### 歷史交接快照 crosswalk

| 歷史快照 | 現行判讀 |
|---|---|
| `InterSubMod/docs/handoff/20260805_系統交接與驗收_01/README.md` | `HISTORICAL / SUPERSEDED`；只供在完整 repository 中追蹤當日資料／程式／輸出判斷，不是本 standalone 包的必要依賴，也不可引用為 current gate。 |
| `InterSubMod/docs/handoff/20260806_LongLineage充分性稽核與路線裁決_01.md` | `HISTORICAL / SUPERSEDED`；現行以 candidate `b9aaa12`、private safety stack 與 P3/P4/P5/P7/P8 blocker set 為準。 |
| `InterSubMod/docs/handoff/20260806_兩repo端到端串接可行性盤點_01.md` | `HISTORICAL / SUPERSEDED`；現行 I/O 與 capability boundary 以本包「軟體輸入輸出與研究流程」為準。 |

## 固定狀態

| 項目 | 固定值／狀態 |
|---|---|
| InterSubMod baseline | `ddd8909a838318d8a77969313e9561c8ff9d01c2` |
| Excluded InterSubMod work | `73afaeac` 與 2026-08-13 drilldown/CNV dirty work：`IN_PROGRESS/PARTIAL` |
| LongLineage preview candidate | `b9aaa12`；research preview only |
| Excluded LongLineage work | `9ad976b`、`6ce62b2` 與未提交修改；另開科研 PR |
| Authority scope | 7 technical datasets／6 biological IDs／chr1-22 |
| Frozen byte replay | 19 MATCH／0 missing／0 mismatch |
| Public-claim audit | 初始58/158有問題、34 P0；完整working stack內33/33 local P0 guards ready，第34筆`C108` GitHub About已由live re-fetch取得bounded live `CONFIRMED_WITH_LIMITS`。完整registry現為69 `CONFIRMED`、69 `CONFIRMED_WITH_LIMITS`、20 `UNVERIFIED`；default branch、Wiki與Pages尚未發布並重抓，因此publication/release仍`BLOCKED`。 |
| LongLineage readiness | P3/P4/P5/P7/P8 `BLOCKED`；production `run` intentionally blocked |
| bip7 | local data preflight已取得 receipt；仍須合併後 fresh-clone全鏈驗收 |
| bip8 | 尚無在 bip8 hostname 上取得的 receipt；`BIP8_DATA_PREFLIGHT_BLOCKED` |
| Release/tag | `BLOCKED`；不得建立 `research-handoff-2026.08.1` tag／Release，直到所有 gate PASS |
| Large-asset data plane | 69,606,961-byte HTML已archive-first移至local prepublication archive並由Git改redirect；剩餘101個tracked >1 MiB資產尚待逐筆producer/finality/license與Release target裁定，`LARGE_ASSET_MIGRATION_BLOCKED` |

## Enum 快速判讀

- `evidence_status`: `AUTHORITY | VALIDATED_DERIVED | PARTIAL | HISTORICAL | INVALIDATED | IN_PROGRESS`
- `scope`: `FULL | PARTIAL | DEMO`
- `finality`: `FINAL_FOR_SCOPE | NON_FINAL | SUPERSEDED`
- `availability`: `GIT | GITHUB_RELEASE | LOCAL_CANONICAL | EXTERNAL_SOURCE | MISSING`
- workflow `classification`: `SUPPORTED | REPRODUCIBLE_LEGACY | ARCHIVED`

Workflow registry中的`REPRODUCIBLE_LEGACY`不代表腳本已壞；它只表示本交接不承諾該檔案具bip7／bip8 portable support。`scripts/pipeline/steps/04_cleanup.sh`因含直接刪除而固定不在`SUPPORTED` allowlist；portable pipeline須加`--skip-cleanup`，直到archive-first replacement完成驗證。

`FINAL_FOR_SCOPE` 只回答「在明示 scope 內是否固定」，不回答「是否具生物真實性」「是否 production-ready」「是否可以外推」。

Artifact registry共有20筆`FINAL_FOR_SCOPE`：其中19筆是frozen science/source byte authority（`AUTHORITY`），另1筆是append-only source adjudication（`VALIDATED_DERIVED`，只對provenance correction final）。後者不可被計作新增science authority或生物結果。

## 目前 gate verdict

| Gate | Verdict | 理由 |
|---|---|---|
| Frozen authority bytes | **PASS** | 19/19 SHA-256 MATCH |
| Structural repo hygiene | **PASS** | isolated HEAD 0 absolute/broken symlink；tracked local settings移除且可復原 |
| Large-asset data plane | **BLOCKED** | 原1個>50 MiB HTML已local archive驗hash並改redirect；目前仍有101個tracked >1 MiB、合計499,500,347 bytes待逐筆裁定／遷移。見[registry](registries/tracked_large_asset_registry.json)與[migration receipt](evidence/large_asset_migration_receipt.json) |
| Unified inventory | **PASS** | 51 physical runs、36 artifacts、16 machine paths、11 storage roots；core registry build receipt的schema checks PASS。Whole-package schema/test count由各驗收receipt動態產生，不把某個早期subset的數量當成整包總數。 |
| Clean build/C++/Python | **LOCAL PASS / HOSTED PENDING** | clean snapshot CTest與完整Python suite PASS；精確test/suite count由當次CTest/pytest receipt動態產生，GitHub-hosted receipt仍待PR執行 |
| Synthetic E2E | **DEMO PASS** | bip7由portable commit `fb806d9b` clean clone重播，驗199欄、12 read leaves與binary/output SHA-256；不得取代science |
| bip7 fresh-clone | **BLOCKED** | real-data paths存在但4個local checksum locator缺失；完整 fresh-clone receipt待補 |
| bip8 fresh-clone | **BLOCKED** | 尚未在 hostname=bip8執行 |
| 158 claims/source | **BLOCKED** | 完整working stack的33/33 local P0 guards PASS；第34筆`C108` About已bounded live `CONFIRMED_WITH_LIMITS`。Registry為69 confirmed、69 confirmed-with-limits、20 unverified；P1/P2 source closure與default branch／Wiki／Pages live驗證未閉合 |
| GitHub live surfaces | **BLOCKED** | About `C108`已bounded live確認；main／Wiki／Pages尚未綁定同一 release commit |
| LongLineage public preview | **BLOCKED** | private-first三層draft PR已建立；source/license/history blockers仍在，repo目前為private |
| Meeting notes / media curation | **PARTIAL / BLOCKED** | [storage-root manifest](registries/storage_root_manifest.json)只盤點Meeting root的10個 immediate children；尚未建立逐件決策、結論、figure provenance、finality與license curated Markdown index。原始PPTX/PDF/media未進Regular Git。 |
| Tag/GitHub Release | **NO-GO** | 上述任一 blocking gate存在即不得發布 |

交接包內嵌manifest、registry、receipt與來源裁定，並用relative links相互連結。Frozen authority的19個底層payload依資料分層政策留在local canonical/frozen data plane，不複製進Git；[authority replay receipt](evidence/authority_replay_receipt.json)保留其絕對路徑與實測hash。因此「包可獨立解讀」不等於「包含所有大型payload」。

Meeting 目前也只完成 root-level inventory，沒有足夠證據將任一舊簡報宣告為已被某份repo report逐件取代。新人應以本交接索引、2026-08-01 frozen authority與current registries為主；後續curation必須對每份PPTX/PDF/media記錄owner、決策/結論摘要、figure source、license、finality、supersedes與目標Git/Release/archive channel，才能關閉 `MEETING_NOTE_CURATION_PARTIAL`。

## 變更與發布政策

本交接包採 stacked、single-intent PR。讀者看到 branch/PR PASS，只代表該 PR 的 bounded gate；必須在整個 stack 合併後，從同一 commit重跑 aggregate acceptance，才可評估 tag／Release。任何真實 secret finding 都先撤銷憑證並另案處理 Git history；本輪不重寫歷史。

---

**Partial flag**：本索引的資料根盤點是 comprehensive inventory，但執行驗證仍有明列 blocked gate；不得把本檔單獨當 release receipt。
