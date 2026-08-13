<!--
建立時間: 2026-08-13 14:00 +08:00
目標: 逐步記錄完整研究交接包實作中的決策、偏離、折衷、驗證與未決事項
處理範圍: InterSubMod research handoff snapshot + LongLineage private research preview
關聯檔案:
  - InterSubMod/research/20260813_complete_research_handoff/implementation-notes.md
  - InterSubMod/docs/handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md
  - InterSubMod/docs/handoff/20260813_完整研究資料與軟體交接_01/registries/registry_build_receipt.json
驗證方式: Step→Verify；所有 PASS 必須對應 command、exit code、receipt或remote check
證據等級: engineering/research-governance record；不新增 science claim
狀態: IN_PROGRESS / RELEASE_BLOCKED
-->

# 完整研究交接 Implementation Notes

用 Claim–Evidence–Verdict：**交接包已閉合 frozen authority、structural repo hygiene、統一 registry 與 LongLineage private-first safety foundation；大型tracked asset migration、雙機、live publication與production release仍未閉合。**

## 設計決定

- [決策] 任務固定為 `(B) Comprehensive validation + (D) External handoff`；本輪交付是`research handoff snapshot`，不是production release。
- [決策] InterSubMod release基準固定`ddd8909a838318d8a77969313e9561c8ff9d01c2`；`73afaeac`與08-13 CNV/drilldown只登錄為`PARTIAL/IN_PROGRESS`。
- [決策] LongLineage candidate固定`b9aaa12a11fa00606bd174dabd0f172a5d112359`；post-candidate science不混入三個public-safety PR。
- [決策] 2026-08-01 authority manifest與denominator registry保持immutable；後續來源角色更正採append-only crosswalk。
- [決策] 所有final判定只由registry欄位與provenance作出；禁止依資料夾名稱、mtime或HTML視覺猜測。
- [決策] Tiny fixture固定為`DEMO`，只能驗clone/build/run/schema；不寫入scientific evidence ledger。
- [決策] LongLineage在source/license/history blockers歸零前維持`PRIVATE`，不建立tag或GitHub Release。
- [決策] `FINAL_FOR_SCOPE` artifact 的`regeneration_command`必須用`REPLAY_ONLY:`、`VERIFY_ONLY:`或`NOT_REGENERABLE_FROM_HANDOFF:`明示語意；不得把hash replay假稱science rerun。

## 已完成驗證

| 項目 | 輸入／命令 | 實際輸出 | Verdict |
|---|---|---|---|
| Authority replay | `scripts/handoff/replay_authority.py`，輸入2026-08-01 authority | `expected=19, match=19, missing=0, mismatch=0` | PASS；bytes only |
| InterSubMod clean build／CTest | `/tmp/ism-final-build.Bo3Oit`，Release configure/build＋dynamic CTest inventory | build exit 0；270/270 CTest PASS | engineering PASS；既有compiler warnings不升格為science證據 |
| InterSubMod official Python suite | `/tmp/ism-handoff-py310/bin/python -m pytest tests -q` | 326 passed、14 subtests、0 failed；43個legacy-schema warnings | official CI scope PASS |
| Structural repo hygiene | audit→archive-first cleanup→isolated HEAD scan | 0 broken/absolute symlink、0 tracked local settings | PASS；不涵蓋large-asset data plane |
| Large-asset data plane | `git ls-tree` baseline 102 files／569,107,308 bytes；archive-first移出唯一>50 MiB HTML | archived HTML 69,606,961 bytes SHA `36548d...c14`；現行registry 101 files／499,500,347 bytes | `LARGE_ASSET_MIGRATION_BLOCKED`；GitHub Release URL尚不存在 |
| Unified registry | `python3 scripts/handoff/build_registries.py ...` | 51 physical runs、36 artifacts、16 machine paths、11 storage roots；6/6 schemas | PASS |
| Registry tests | `/tmp/ism-handoff-py310/bin/python -m pytest -q tests/test_handoff_registries.py` | 12/12 PASS | PASS |
| Tiny synthetic E2E | clean clone `7c7fbd6f`＋isolated Release binary＋12-read fixture | exit 0；199 header/data columns、12 tree leaves、66/66 valid pairs；summary/full-receipt hashes固定 | DEMO PASS；不是real-data、machine acceptance或release gate |
| bip7 doctor | `scripts/site/doctor ... --mode real-preflight` | exit 5；4 local checksum locator missing | `BIP7_DATA_PREFLIGHT_BLOCKED` |
| bip8 | 尚無hostname=bip8 receipt | NFS view only | `BIP8_DATA_PREFLIGHT_BLOCKED` |
| Public claims | 158-row validator＋GitHub About live re-fetch | 69 CONFIRMED、69 CONFIRMED_WITH_LIMITS、20 UNVERIFIED；33/33 local P0 guards ready，`C108`為bounded live `CONFIRMED_WITH_LIMITS` | P0 source/About bounded gate PASS；P1/P2與publication/release仍blocked |
| Public HTML browser QA | Chromium，21頁×desktop/mobile/no-JS/print | 84/84 PASS；0 overflow、0 page/console error、0 external request | local source QA PASS；不是live Pages證據 |
| Algorithm／CLI crosswalk | 35-row immutable matrix＋commit/hash-bound source inspection | 35/35 rows；6 CONFIRMED、27 CONFIRMED_WITH_LIMITS、2 UNVERIFIED；validator 0 error | static asset-ready PASS；不是aggregate publication/release PASS |
| Reader acceptance contract | Draft 2020-12 schema＋15-file dual hash binding＋negative tests＋無前情 fresh agent | 6/6 questions、8/8 anti-promotion、15/15 current/tested-commit hashes，receipt validator PASS | reader comprehension PASS；不是science、host或publication驗證 |
| Package validator | `/tmp/ism-handoff-py310/bin/python scripts/handoff/validate_handoff_package.py <PACKAGE_ROOT>` | 15/15 PASS；reader與tiny E2E receipts均納入fail-closed checks；23 evidence records、0 missing/hash mismatch | local package structure PASS；release gates仍blocked |
| LongLineage clean foundation | clean clone `f60b5f3` configure/build/check_all | 49/49 CTest、no-network、FOUNDATION_PASS | foundation PASS；release blockers remain |
| LongLineage private docs closure | clean `0805bd5`；只含CHANGELOG、CURRENT_FOCUS與archived task receipt | fresh build exit 0、CTest 49/49、production `run`為預期`KernelBlocked` exit 6；worktree clean、未push/tag/release | engineering boundary PASS；`b9aaa12`仍是public-preview candidate，公開gate仍BLOCKED |
| LongLineage strict public gate | immutable audit base from safety receipt | expected exit 1；5 blocker classes | safe FAIL；keep private |

## 偏離與更正

- [偏離] 舊master manifest的「19」其實是1 header＋18 data rows；physical inventory為35 current＋16 pending archive。已改為logical metadata合併physical record，不發明第19個run。
- [偏離] 舊manifest的`tagged_bam_ready=false`與實體storage不一致。Targeted stat確認14 BAM／3,709,322,840,333 bytes；其中7個`paired_full`合計1,840,983,466,353 bytes，重播2026-07-11首／中／尾1 MiB sampled-identity receipt為7/7 MATCH。因未計full-file SHA-256、pre-fix/truth-aware歷史及未建立與現行sidecar的producer revision／內容等價，仍固定`PARTIAL/NON_FINAL`；294.2669×僅為跨世代footprint quotient，舊287×無效。
- [偏離] `InterSubMod_big7_runbook`廣泛搜尋仍找不到。Registry標`MISSING_SOURCE/HISTORICAL`，reconstructed index不可冒充原件。
- [偏離] LongLineage在稽核期間兩度被觀察為`PUBLIC`並收回`PRIVATE`。目前private不能證明暴露期間沒有clone/download；事件保留於safety receipt。
- [偏離] LongLineage PR #6最初以moving `origin/main`作history base，hosted CI 5 jobs fail。已改讀safety receipt的immutable `5daf50f…`，新版9/9 hosted checks PASS。
- [偏離] Unified machine-path registry一度仍記錄LongLineage preview `f60b5f3`且dirty=3。等docs-only closure `0805bd5`完成後，已在InterSubMod clean `15a6493b`時重建：兩個release worktree皆clean，LongLineage HEAD=`0805bd5`；這是`last_seen_at`操作快照，不把`0805bd5`升格為public candidate。
- [偏離] 一次full CTest誤用只build governance target的build dir，造成多個`Not Run`；該次明確作廢，不列PASS/FAIL證據。正式clean-clone full build後重跑49/49 PASS。
- [偏離] 直接執行repo-root bare `pytest -q`會收集`docs/`與`research/`內的frozen script/test素材，並因跨checkout import isolation產生11個collection errors；正式CI明確使用`pytest tests -q`，最終該scope為326 passed。後續可另案決定是否以`testpaths=tests`固定bare CLI語意，不能把collection錯誤冒充演算法回歸。
- [偏離] `EVIDENCE_MANIFEST.json`曾漏登3個tagged-BAM sampled-identity檔，且frozen-binary adjudication保留舊SHA。已補回baseline／replay manifest／replay receipt的`HISTORICAL|VALIDATED_DERIVED + PARTIAL + NON_FINAL`契約，並以現行exact-copy SHA重新綁定；package validator回到13/13 PASS。
- [偏離] Frozen binary adjudication曾把5/5 solver module byte identity與19/19 authority replay寫成「論文核心數字可從版本控制完整重建」。已改為binary/source identity與重編骨架；TiB/local data plane的science未全量重算。
- [偏離] 雙機操作文件的authority replay範例曾誤用registry builder的`--authority`參數；實際`replay_authority.py`介面是`--manifest`。首次命令exit 2且未產生receipt；修正後read-only replay為19/19 MATCH，並新增CLI文件regression test。
- [驗證事件] 21頁×4模式Chromium首輪為80/84 PASS；兩份handoff HTML只在mobile/no-JS出現table橫向overflow，desktop/print、HTML/XML/SVG parse、page/console error及external request均通過。該輪維持FAIL CLOSED，不能當最終browser receipt。

## 折衷

- [折衷] 本輪不讀取3.7 TB tagged BAM計完整SHA；保存exact path、bytes、mtime、生成責任與`NON_FINAL`，避免把storage inventory冒充byte authority。
- [折衷] 不重跑TiB級science；以19/19 frozen byte replay、clean build、synthetic E2E與雙機preflight作工程可重現性。這不等同scientific recomputation。
- [折衷] Meeting原始PPTX/PDF/media不直接塞入Regular Git；只有決策、結論與figure provenance整理稿才可進Git，其餘依finality進Release或local archive。

## 未決與發布阻塞

- [完成] InterSubMod clean build、動態CTest 270/270、official Python 326/326與commit `7c7fbd6f` fresh-clone tiny E2E均通過；eventual release commit仍須由外部／CI receipt自行綁定，避免receipt-in-commit自我參照。
- [完成] 無對話背景 fresh-agent reader receipt 已在固定15檔範圍通過：6問、8個禁止升格、15/15 current/tested-commit SHA一致；claim ceiling僅為reader comprehension。
- [未決] 101個tracked >1 MiB資產仍須逐筆producer、finality、license與Release target裁定；目前只有唯一>50 MiB HTML完成local archive-first保存，尚無Release upload/refetch SHA。
- [未決] bip7補齊4個data checksum locators並完成host-local fresh clone acceptance。
- [未決] bip8必須在`hostname=bip8`重新clone/build/test/smoke/preflight。
- [未決] 20個`UNVERIFIED` public claims、P1/P2 source closure、Wiki、Pages與default branch的live publication一致性；About `C108`已bounded live確認，但不能代表整體publication PASS。
- [未決] LongLineage 4 unresolved source rows、21 source-license reviews、11 dependency license determinations與7 history findings。
- [未決] InterSubMod branch protection／required CI、full-history secret scan與Release assets round-trip SHA。
- [未決] CN/LOH integration、cell barcode／orthogonal clone truth；未完成前`confirmed cellular subclone=0`、`confirmed linear ancestry=0`。

## 發布判定

`RELEASE_BLOCKED`。不得建立`research-handoff-2026.08.1` tag／GitHub Release，也不得把LongLineage切為public。只有所有阻塞項各自有immutable receipt且aggregate acceptance在同一commit重播PASS，才可重新評估。
