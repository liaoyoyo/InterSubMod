<!--
建立時間: 2026-07-16
目標: 獨立複核 M2 resource-gate、frozen release、resume 與 pilot/full verifier 的科學 provenance 一致性
處理範圍: Task Type B；程式、Runbook 與 synthetic/no-BAM fixtures；不執行 formal RUN_ROOT 或 154+154 真實計算
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2正式執行Runbook_02.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_extraction.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_ranking.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_m2_single_task_pilot.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_full_m2_receipts.py
-->

# M2 resource-gate 科學 provenance 一致性獨立 red-team

> **READER-TEST GO**：修補後targeted `117/117`與full regression `327/327`全部PASS；checkpoint aggregate/checks自洽重簽與ranking consolidated candidate內容＋metadata＋terminal自洽重簽的負向反例均被runner本身拒絕。P0=`0`、open P1=`0`、已修補P1=`5組`、P2=`2`。這代表可依Runbook嘗試formal execution，不代表`154 extraction + 154 ranking`真實數據已完成；formal actual-data仍為`WAIT`。（影響：高；reader-test：GO；formal actual-data：WAIT）

## 1. 任務分類、範圍與 claim boundary

- **Task Type**：B — Comprehensive validation。
- **服務目標**：G4 多資料集一致性與可重現；G5 外部可驗證證據鏈。
- **本輪輸入**：4 個 resource/provenance scripts、4 個 targeted test files、Runbook v2。
- **本輪不做**：不讀 BAM/VCF、不啟動 formal compute、不建立或修改 formal RUN_ROOT。
- **判讀邊界**：unit tests／synthetic fixtures 證明協定 fail-closed，不等於真實 `154 extraction + 154 ranking` 已完成。正式最終數字必須等兩份 terminal receipts、POST 與 independent final verifier 全部存在且通過。

## 2. 權威入口與數量口徑

原需求中的「`verify_m2_final_resource_gates.py`」目前不是實際檔名；對應權威入口為：

| 職責 | 實際權威檔案 |
|---|---|
| full extraction runner／gate capture | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_extraction.py` |
| full ranking runner | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_ranking.py` |
| direct-pilot independent verifier | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_m2_single_task_pilot.py` |
| final full independent verifier | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_full_m2_receipts.py` |

固定數量重新計算如下：

| 層級 | 計算 | 預期 |
|---|---:|---:|
| 每 stage tasks | 7 technical datasets × chr1–chr22 | 154 |
| completed cumulative checkpoints | 8, 24, 40, 56, 72, 88, 104, 120, 136, 152, 154 | 11 batches |
| 每 full stage resource gates | 1 session + 11 batch | 12 |
| extraction + ranking full gates | 12 + 12 | 24 |
| 每 direct pilot gates | extraction + ranking B0 + ranking B20 | 3 |
| 兩個 direct pilots | 2 × 3 | 6 |
| pilot GO receipts | 2 pilots ×（B0、B20） | 4 |

## 3. 已通過的靜態 provenance 鏈

### 3.1 Frozen release 與 canonical manifest

- extraction、ranking、full verifier 都固定 canonical manifest SHA-256：`16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45`。
- release loader要求 schema／authority／scope／完整 source-role map、canonical schema、canonical immutable copy、PRE copy、entrypoints、source path＋SHA、0444／single-link 與 adjacent frozen freezer deep verification。
- full extraction的 session gate直接綁完整 release-manifest identity與run-contract semantic SHA；ranking session另綁 parent extraction terminal identity。各 batch gate再綁 session SHA、canonical selected-task IDs與previous chain head，形成遞迴 provenance。
- pilot verifier要求 extraction gate、ranking gate、extraction child provenance三者使用同一 canonical manifest identity，並要求兩 gate使用完全相同的 frozen release-manifest identity。

### 3.2 Full orchestration／resume

- runner只接受合法 completed checkpoint prefix，不接受跳號或重排。
- 每個 completed batch都重驗 batch-start、resource gate、grants、completions、child receipts／declared outputs與checkpoint chain。
- open/orphan batch、resource gate重用／path swap、未引用 grant/completion、timeout marker、root inode交換與額外JSON gate都 fail-closed。
- final verifier獨立重建 11-batch chain，要求各 stage正好154 completions、12 authenticated gates、154 child receipts與輸出重雜湊。

### 3.3 Resource gate 本體

- schema與exact-key contract固定；gate identity包含 path、raw SHA、semantic SHA、gate ID、sidecar path與sidecar SHA。
- gate要求zero-conflict process snapshot、`f_bavail × f_frsize`可用空間、精確300 GiB reserve、producer source、target與host boot／monotonic timestamp。
- claim僅為launch-time snapshot；在non-hostile same-UID threat model下不是cryptographic process ancestry。

## 4. P0／P1／P2 findings

### P0

- **目前 0 項**。

### P1-1 — extraction terminal `REUSED_FINAL` 深驗不足（已修補）

Pre-fix `validated_existing_final()`只驗sidecar、schema、`all_pass`、154與run-contract equality；若terminal aggregate被自洽改寫並重新sidecar，而children仍原樣，早退路徑未直接從children重建terminal payload。

修補後先把154個terminal rows逐一綁回canonical reusable children，再用同一production builder重建完整terminal receipt並要求exact equality；release mode的orchestration chain已在此前由`load_release_orchestration_state()`重驗。重簽後的duplicate-for-missing result與aggregate tamper皆被拒絕。

- final cumulative runner SHA：`cf016b9a046c214bbefb6a4b2509955910710fce73d3186dce27b666d5c40fc4`
- final cumulative test SHA：`dbd6b0e44d80e9aea8babd3b31866a625b5817971dcbb99704f27ca078cd88e4`

### P1-2 — ranking terminal `REUSED_FINAL` 深驗不足（已修補）

同P1-1；pre-fix ranking早退只做terminal外層contract。修補後逐一綁定154個compact child rows、重算aggregate／conservation／runtime diagnostics、重雜湊且重算candidate-table semantic SHA，再用production terminal builder要求exact equality。duplicate-for-missing、aggregate tamper與candidate semantic-SHA tamper皆被拒絕。

- final cumulative runner SHA：`66bb175404c207ef320f213c650bb10c6d5fcf3c84cbc40b8ca25e68604da767`
- final cumulative test SHA：`ae8b324c6c18f845e1e49e65944c00110cb2887bb147f27e22502f81946d5103`

### P1-3 — direct evidence mode與Runbook immutable contract不一致（已修補）

三個同型deterministic blockers合併處理：

1. pre-fix pilot verifier的`write_receipt()`以一般`Path.write_text()`發布，JSON與sidecar在`umask 022`下是`0644`。
2. direct extraction／ranking繞過full-run wrapper，因此child receipt／sidecar也是`0644`。
3. PRE／POST producer以`mkstemp` hardlink發布，receipt／sidecar是`0600`。

Runbook的`authenticate_sidecar`要求`0444`、`nlink=1`；三者原本都會在內容PASS後deterministic fail。

**暫存重現**（未碰formal root）：

```text
{'receipt_mode': '0o644', 'receipt_nlink': 1,
 'sidecar_mode': '0o644', 'sidecar_nlink': 1,
 'runbook_requires_mode': '0o444'}

PRE/POST: receipt_mode=0o600, sidecar_mode=0o600
direct extraction/ranking: receipt_mode=0o644, sidecar_mode=0o644
```

**修補**：pilot verifier改成`O_CREAT|O_EXCL` fd writer、file與directory `fsync`、`fchmod(0444)`與single-link檢查；任一半preseed都拒絕，並移除`--overwrite`。Runbook另加入共用`seal_direct_receipt`：chmod前先驗regular／non-symlink、JSON與sidecar各`nlink=1`、same-directory basename、sidecar單行嚴格格式與實際SHA；之後才chmod `0444`並呼叫`authenticate_sidecar`。此helper統一用於PRE、POST、direct extraction與direct ranking；resume另重新authenticate child與GO receipts。

**修補後SHA**：

- verifier：`9d15ce2bf15af5cc2c4c690cd7718b131108fd8e3946f6a72da40487b06f1578`
- test：`0e622ed2230fb3cc5178008475dce9ebcdcb24d25c51878ab2afc3821c0adad9`
- final cumulative Runbook：`576caae20be596f2ba879a8c5b6f9ecce49d401da41335a705c691b35da9cba1`

**驗證**：shell helper實跑valid pair後兩檔皆`0444/nlink1`；hardlinked target在chmod前被拒絕；partial preseed、重複發布與overwrite入口皆被拒絕。

### P1-4 — completed checkpoint scientific payload可自洽重簽後通過resume（已修補）

`load_release_orchestration_state()`會深入驗session、batch-start、gate、grant、completion、result task set與chain head，但pre-fix版本未從已驗children重建每個completed checkpoint的完整aggregate／checks。獨立fixture將`checkpoint_008`的aggregate／checks改寫並重簽sidecar後，extraction與ranking loaders都接受。

這不一定直接改變terminal final數字，但completed checkpoint是合法resume權威；既然runner接受並宣告chain有效，checkpoint本身就必須與8個已驗children做deep equality，不能只把orchestration envelope當真。

**修補**：兩runner現在都從已驗證的child receipt與completion identity重建當前cumulative checkpoint。extraction使用deterministic compact child rows與production `build_extraction_receipt()`；ranking使用production `build_ranking_checkpoint()`。重建後要求完整receipt（aggregate、checks、canonical results、orchestration與integrity envelope）exact equality。負向fixture將`checkpoint_008`的aggregate與checks修改後重簽sidecar，extraction與ranking的loader均fail-closed。

- extraction runner/test SHA：`cf016b9a046c214bbefb6a4b2509955910710fce73d3186dce27b666d5c40fc4` / `dbd6b0e44d80e9aea8babd3b31866a625b5817971dcbb99704f27ca078cd88e4`
- ranking runner/test SHA：`66bb175404c207ef320f213c650bb10c6d5fcf3c84cbc40b8ca25e68604da767` / `ae8b324c6c18f845e1e49e65944c00110cb2887bb147f27e22502f81946d5103`

### P1-5 — ranking consolidated candidate可自洽改寫後被`REUSED_FINAL`接受（已修補）

獨立fixture把實體consolidated candidate的一列`profile_log_likelihood`由`-1`改為`-999`，同步更新size／raw SHA／semantic SHA、terminal metadata與sidecar。pre-fix `_verify_existing_candidate_table()`只驗consolidated table自身identity與schema，沒有從154個child candidate tables逐列重建，因此`validated_existing_final()`仍回傳有效，runner main會印出`REUSED_FINAL all_pass`。

後續independent final verifier預期會擋住此漂移，但terminal runner自身不應先宣告成功；修補必須把resume terminal candidate table與154 child sources做exact reconstruction或等價deep comparison。

**修補**：build與resume verifier共用child→canonical row iterator。resume先重驗154個child candidate artifact的path／size／SHA，再將實體consolidated gzip table與child-derived expected rows逐列lockstep比對，同時重算row count、unit count、max group size與semantic SHA。新fixture建立154個child artifacts並實際寫入1列candidate；將`profile_log_likelihood=-2.0`改成`-999.0`，再同步更新physical SHA、semantic SHA、terminal metadata與sidecar，`validated_existing_final()`仍拒絕，且`main()`不會輸出`REUSED_FINAL`。

- ranking runner SHA：`66bb175404c207ef320f213c650bb10c6d5fcf3c84cbc40b8ca25e68604da767`
- ranking test SHA：`ae8b324c6c18f845e1e49e65944c00110cb2887bb147f27e22502f81946d5103`

### P2-1 — pilot GNU-time原始檔可再加強publication-boundary重驗

Runbook會把`*.time.txt`設為0444且pilot verifier在GO receipt中保存path／size／SHA；但既有GO receipt resume路徑不會再次把目前time file重雜湊與receipt內SHA比對，也沒有獨立sidecar。這不會改變已sealed GO receipt當時使用的數值，但會降低後續原始resource evidence的drift可見性。建議未來將time files一併exclusive sidecar封存，或在每次resume重驗receipt內SHA。

### P2-2 — full verifier的extra-gate集合只枚舉`*.json`

所有必要gate sidecars都會被逐一驗證；但resource-gate namespace的exact-set檢查只枚舉`*.json`，孤立的額外`.json.sha256`不會被列為extra gate。它不能冒充有效gate，且會讓未來exclusive發布fail-closed；屬namespace hygiene／diagnostic hardening，不是科學證據誤接受。

## 5. 已執行命令與實際結果

### 5.1 輸入路徑

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/{run_full_m2_extraction.py,run_full_m2_ranking.py,verify_m2_single_task_pilot.py,verify_full_m2_receipts.py}`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_{full_m2_extraction,full_m2_ranking,verify_m2_single_task_pilot,verify_full_m2_receipts}.py`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2正式執行Runbook_02.md`

### 5.2 Targeted resource/provenance suite

```bash
python3 -m unittest -v \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_full_m2_extraction.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_full_m2_ranking.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_verify_full_m2_receipts.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_verify_m2_single_task_pilot.py
```

實際結果：`117/117 PASS`；test elapsed=`31.942 s`；`/usr/bin/time` wall=`32.88 s`、max RSS=`76,812 KiB`、exit=`0`。第一次重驗唯一failure是Runbook尚保留兩runner修補前SHA；更新這兩個allowlist SHA後117項全數PASS，科學反例沒有被隱藏。

### 5.3 Full research regression

```bash
python3 -m unittest discover -v \
  -s research/20260716_read_linked_hypercube_exact_likelihood_validation/tests \
  -p 'test_*.py'
```

實際結果：`327/327 PASS`；test elapsed=`161.383 s`；`/usr/bin/time` wall=`162.52 s`、max RSS=`111,492 KiB`、exit=`0`。測試前後重算四個修補檔SHA，4/4一致。

### 5.4 Runbook reader／CLI checks

| 檢查 | 結果 |
|---|---:|
| 從Runbook本身抽出exact-11 allowlist並`sha256sum -c` | 11/11 OK |
| Bash fences | 13 |
| 合併fences後`bash -n` | exit 0 |
| extraction／ranking／pilot／final verifier `--help` | 4/4 exit 0 |
| direct sealer valid fixture | PASS；0444／nlink1 |
| direct sealer hardlink fixture | PASS；chmod前拒絕 |

## 6. 最終已驗證 SHA-256

| 檔案 | SHA-256 |
|---|---|
| `run_full_m2_extraction.py` | `cf016b9a046c214bbefb6a4b2509955910710fce73d3186dce27b666d5c40fc4` |
| `run_full_m2_ranking.py` | `66bb175404c207ef320f213c650bb10c6d5fcf3c84cbc40b8ca25e68604da767` |
| `verify_m2_single_task_pilot.py` | `9d15ce2bf15af5cc2c4c690cd7718b131108fd8e3946f6a72da40487b06f1578` |
| `verify_full_m2_receipts.py` | `87fc3ddde1052cd32b64588aa3e483faa699474047bde52e1332a07f779e5558` |
| `test_full_m2_extraction.py` | `dbd6b0e44d80e9aea8babd3b31866a625b5817971dcbb99704f27ca078cd88e4` |
| `test_full_m2_ranking.py` | `ae8b324c6c18f845e1e49e65944c00110cb2887bb147f27e22502f81946d5103` |
| `test_verify_m2_single_task_pilot.py` | `0e622ed2230fb3cc5178008475dce9ebcdcb24d25c51878ab2afc3821c0adad9` |
| `test_verify_full_m2_receipts.py` | `7c72fe7e5b1f515839459890c51a5aef159addb6133c252b05d64f4dd7510e42` |
| `20260716_M2正式執行Runbook_02.md` | `576caae20be596f2ba879a8c5b6f9ecce49d401da41335a705c691b35da9cba1` |

## 7. 目前 verdict

| 維度 | 狀態 |
|---|---|
| resource/provenance contract | **GO；可依Runbook啟動，但每個live resource gate仍必須當場PASS** |
| P0／open P1 | **0／0** |
| resolved P1 | **5組；負向fixture與完整regression通過** |
| P2 hardening | **2；不造成科學證據誤接受** |
| reader-test | **GO；targeted 117/117、full 327/327、exact-11 11/11** |
| formal 154+154 actual-data evidence | **WAIT／不得宣稱完成** |
