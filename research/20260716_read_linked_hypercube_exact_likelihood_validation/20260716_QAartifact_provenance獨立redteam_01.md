<!--
建立時間: 2026-07-16 17:57:48 CST (+0800)
目標: 獨立 red-team 審核教授版 HTML browser-QA artifact 的 producer/content provenance、原子發布、失敗保留與 verify-only 契約。
處理範圍: Task Type B comprehensive-validation 的 QA 程式、單元測試與 presentation freezer 整合邊界；不讀取 BAM、不啟動正式 M2 RUN_ROOT、不產生生物數據。
服務目標: G4 多樣本一致性與 reproducibility；G5 可被外部重算的證據鏈。
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/qa_validated_html_report.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_qa_validated_html_report.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_presentation_snapshot.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_presentation_snapshot.py
-->

# QA artifact provenance 獨立 red-team

## TL;DR

**QA artifact provenance 修補後通過：P0=0、P1=0、P2=3；可進入正式 presentation freeze／HTML／QA 流程，但仍須由 frozen QA entrypoint 執行並在 QA 後跑 presentation `verify-only`（影響：高，信心：高）。**

用 PREP 敘述：

- **Point**：成功 QA 現在精確綁定被渲染的 HTML、實際執行的 QA producer、三份 rendered outputs，以及 receipt／sidecar 的固定五檔 layout。
- **Reason**：所有檔案均以 path、byte size、SHA-256、mode、`st_nlink`驗證；publish 使用 `renameat2(RENAME_NOREPLACE)`，receipt／sidecar 使用 `O_EXCL`。
- **Evidence**：審核時找到一個可重現的 sealing TOCTOU P1；修正後新增 HTML 與 producer 兩個 deterministic race tests，合併 presentation freezer 共 **24/24 tests PASS**；實際 Playwright transient smoke 的 frozen run／same-frozen verify 均 exit 0，而 live producer verify exit 2。
- **Point**：本稽核只判定 QA provenance 工程契約可用，不代表正式 154＋154 M2 數據、教授版 FINAL HTML 或任何生物結論已完成。

## 1. 任務分類、輸入與輸出

| 項目 | 定義 |
|---|---|
| Task type | **B — Comprehensive validation** |
| 服務目標 | **G4 / G5** |
| 審核單位 | QA producer、QA tests、presentation freezer integration boundary |
| 正式資料範圍 | **無**；未讀 BAM、未使用正式 M2 receipts |
| 本次可支持 | QA artifact path/content/mode/link provenance、fail-closed、no-replace、verify-only |
| 本次不可支持 | 正式數字、clone 數、拓撲、生物演化或教授版 FINAL report 已完成 |

### Step → Verify

1. 完整讀取 QA producer／tests → 驗證：指出 stable-read、publication、failure-preservation、receipt-validation、verify-only 的確切位置。
2. 攻擊 sealing boundary → 驗證：HTML 或 producer 在最後 capture 前改變時，不能產生有效 receipt。
3. 攻擊 publish boundary → 驗證：競爭者不能被覆寫；sidecar failure 必須使 canonical output 消失並保存失敗現場。
4. 驗證 frozen producer → 驗證：同一 frozen script `verify-only` 通過，live repo script 驗證同 artifact 必須失敗。
5. 列出殘餘限制 → 驗證：P0/P1 歸零後才給 GO，P2 明示且不得誤當正式數據。

## 2. 最終輸入 SHA-256

工作目錄：`/big7_disk/liaoyoyo2001/InterSubMod`

```bash
sha256sum \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/qa_validated_html_report.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_qa_validated_html_report.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_presentation_snapshot.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_presentation_snapshot.py
```

Actual output：

```text
cf66a8cf9b5fcd408ac2a487708617016b47bcd3551ac79881e1869e963d4e16  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/qa_validated_html_report.py
89d2059a09857155c7ebfab1b6101c0e57fd9c902de05e75d0275d4a1b63d3b2  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_qa_validated_html_report.py
639898a882b5d72f9c58ad95701a5afdfa22c0ad37f1aebc8748320d2efbd625  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_presentation_snapshot.py
af77886f98cf64e812b1879b49ebfe96e574396562dc01ebd0e8d8a6a36816b5  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_presentation_snapshot.py
```

本 red-team 依主 agent 授權，只修改前兩個 QA 檔案以修正 P1；未修改 presentation freezer 兩檔。

## 3. Provenance 契約逐項審核

### 3.1 Stable identity 與 alias 防護

`qa_validated_html_report.py:112-170` 使用 `lstat → O_NOFOLLOW open → fstat → read → fstat → lstat`，並核對：

- regular file、非 symlink；
- path 與 opened FD 的 `(st_dev, st_ino)` 相同；
- `st_nlink=1`；
- read 前後 `dev/inode/nlink/size/mtime/ctime` 不變；
- path、size、SHA-256、mode、nlink 形成 artifact identity。

這代表 receipt 不只記「檔名」，而是記錄當時實際讀到的完整 bytes 與檔案狀態。

### 3.2 HTML、producer 與三份 rendered outputs

成功 receipt 的 inputs／outputs 契約位於 `qa_validated_html_report.py:363-435`：

| Role | 綁定欄位 | 額外要求 |
|---|---|---|
| HTML | absolute lexical path、size、SHA、mode、nlink | mode=`0444`、nlink=1 |
| QA producer | absolute lexical path、size、SHA、mode、nlink | 必須等於執行 `verify-only` 的 producer path／identity |
| desktop PNG | path、size、SHA、mode、nlink | 固定檔名、0444、nlink=1 |
| mobile PNG | path、size、SHA、mode、nlink | 固定檔名、0444、nlink=1 |
| print PDF | path、size、SHA、mode、nlink | 固定檔名、0444、nlink=1 |

QA directory 本身必須是 `0555`，且**恰好五檔**：三份 rendered outputs、`browser_qa_receipt.json`、`browser_qa_receipt.json.sha256`。HTML 位於 QA directory 外，但其 identity 被 receipt 精確綁定。

### 3.3 No-replace 與 receipt authentication

- `qa_validated_html_report.py:217-232`：receipt／sidecar 以 `O_CREAT|O_EXCL` 新建，不 truncate 或 replace。
- `qa_validated_html_report.py:243-254`：directory publish 使用 Linux `renameat2(RENAME_NOREPLACE)`。
- `qa_validated_html_report.py:603-626`：publish 後才以 `O_EXCL`建立 sidecar、設定 0444／0555、fsync、執行完整 post-publication verify。
- `qa_validated_html_report.py:630-700`：`verify-only` 重查 exact five-file layout、receipt sidecar、HTML／producer／outputs identity、mode、nlink、symlink、extra files 與 cross-role inode alias。

### 3.4 失敗現場保存

- Pre-publication failure：`qa_validated_html_report.py:497-514` 將可能存在的 success filenames 改名為 `UNAUTHENTICATED_*`，加入 `QA_NOT_PUBLISHED.txt`，以 `failed-staging.*` 保存。
- Post-publication failure：`qa_validated_html_report.py:517-534` 同樣移除有效 success filenames，加入 `PUBLICATION_VERIFICATION_FAILED.txt`，以 `failed-publication.*` 保存。
- Production source 未使用 `unlink`、`remove`、`rmtree`、`rmdir` 或 replace 式覆寫。

## 4. 發現並修正的 P1：sealing identity TOCTOU

### 4.1 原始問題

原始流程先檢查 `pre_html_identity`／`pre_producer_identity`，稍後又重新 capture 最終 identity，卻未要求最終 capture 仍等於 pre-QA identity。可在兩次讀取間發生：

```text
browser 實際使用 OLD bytes
→ unchanged check 讀到 OLD bytes
→ 檔案改成 NEW bytes
→ receipt capture NEW bytes
→ post verify 對 NEW bytes 通過
```

因此 receipt 可能錯綁「未被瀏覽器渲染的新 HTML」或「未實際執行的新 producer」。這是 provenance 主張錯誤，列 **P1**。

### 4.2 修正

`qa_validated_html_report.py:570-578` 現在在 authoritative identity capture 後再次強制：

```text
final HTML (path,size,SHA) == pre_html_identity
final producer (path,size,SHA,mode,nlink) == pre_producer_identity
```

任何差異都在 receipt 建立前 fail-closed；`run_browser_qa()` 的 outer failure path 會保存 `failed-staging.*`，且不留下 canonical output、有效 receipt 或 sidecar。

### 4.3 Deterministic regression

`test_qa_validated_html_report.py:232-289` 精確在第四次 identity read 注入 HTML／producer 變更。兩個測試都要求：

- raise `QaArtifactError`；
- canonical QA output 不存在；
- 恰好一個 `failed-staging.*`；
- 有 `QA_NOT_PUBLISHED.txt`；
- 無有效 receipt／sidecar。

## 5. Targeted tests

### 輸入

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/qa_validated_html_report.py`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_qa_validated_html_report.py`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_presentation_snapshot.py`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_presentation_snapshot.py`

### 執行命令

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /bip7_disk/liaoyoyo2001/miniconda3/bin/python3 -m py_compile \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/qa_validated_html_report.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_qa_validated_html_report.py

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /bip7_disk/liaoyoyo2001/miniconda3/bin/python3 -m unittest -q \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_qa_validated_html_report.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_presentation_snapshot.py
```

### 輸出路徑

單元測試只建立 Python `TemporaryDirectory`；不建立正式 research-round artifact。程式修正輸出為前述 QA producer／test 兩檔。

### Actual output

```text
---------------------------------------------------------------------
Ran 24 tests in 0.366s

OK
```

其中 QA tests 13/13、presentation-freezer tests 11/11；包含成功綁定、HTML／rendered bytes tamper、receipt／sidecar tamper、mode、hardlink、extra file、symlink、resigned path swap、producer drift、preexisting output、failed staging、CLI parsing，以及新增的兩個 sealing races。

## 6. Transient adversarial smoke（非正式數據／非正式 artifact）

> 本節只記錄 code-path 攻擊測試。所有輸入都是臨時合成 HTML／dummy bytes，落在 `/tmp`；**不得當作正式 M2 evidence、正式教授版 QA receipt 或任何研究結果**。

### 6.1 實際 Playwright：same-frozen producer binding

執行介面：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
  /tmp/qa_provenance_redteam_postfix_ag9oebqh/frozen_qa_validated_html_report.py run \
  --html /tmp/qa_provenance_redteam_postfix_ag9oebqh/report.html \
  --output-dir /tmp/qa_provenance_redteam_postfix_ag9oebqh/qa \
  --expect-status final

/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
  /tmp/qa_provenance_redteam_postfix_ag9oebqh/frozen_qa_validated_html_report.py verify-only \
  --html /tmp/qa_provenance_redteam_postfix_ag9oebqh/report.html \
  --output-dir /tmp/qa_provenance_redteam_postfix_ag9oebqh/qa \
  --expect-status final
```

Actual output 摘要：

```text
SUCCESS_RUN_RC 0
VERIFY_FROZEN_RC 0
VERIFY_LIVE_RC 2
FILES ['browser_qa_receipt.json', 'browser_qa_receipt.json.sha256',
       'desktop_full.png', 'mobile_full.png', 'print_qa.pdf']
ROOT_MODE 0o555 HTML_MODE 0444 PRODUCER_MODE 0444
```

同一 frozen path／SHA 的 producer 通過；改用 live repo QA script 驗證同一 receipt 時，actual rejection 為：

```json
{"all_pass": false, "error": "receipt producer path differs from executing verifier"}
```

同一 transient smoke 另以錯誤 status 觸發真實 browser-QA failure：exit 2、canonical output 不存在、`failed-staging.*`為0555、含 diagnostic 與三份 rendered outputs，但無有效 receipt／sidecar。

### 6.2 No-replace race

臨時 monkeypatch 在 publication 前建立競爭者 directory／sentinel，再呼叫真正的 `_rename_noreplace()`。

Actual output：

```text
EXPECTED_FAILURE ... refusing to overwrite QA directory
FOREIGN_SENTINEL do not overwrite
PRESERVED_COUNT 1
VALID_RECEIPT False VALID_SIDECAR False
```

結論：競爭者未被搬移或覆寫；本次 staging 改名為 unauthenticated 並保存。

### 6.3 Sidecar publication failure

臨時 monkeypatch 只在 sidecar `O_EXCL` create 注入 `OSError`。

Actual output：

```text
EXPECTED_FAILURE QaArtifactError ... diagnostics preserved at ...failed-publication...
CANONICAL_EXISTS False STAGING_EXISTS False PRESERVED_COUNT 1
VALID_RECEIPT_EXISTS False VALID_SIDECAR_EXISTS False
UNAUTH_RECEIPT_EXISTS True FAIL_MARKER True
```

結論：普通 Python exception 的 post-rename failure 不會留下可驗證的 canonical success path。

## 7. Destructive-call 靜態掃描

### 命令

```bash
rg -n --pcre2 \
  '(os\.(unlink|remove|removedirs|rmdir|replace)|shutil\.rmtree|Path\.(unlink|rmdir|replace)|\.(unlink|rmdir)\s*\()' \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/qa_validated_html_report.py \
  || true
```

### Actual output

```text
(no matches)
```

字串正規化中的 `str.replace()` 不屬於 filesystem replacement；production QA source 沒有刪除式 API。

## 8. 最終 P0／P1／P2 判定

| Severity | 數量 | 判定 |
|---|---:|---|
| P0 | **0** | 無會產生錯誤正式 scientific result 的已知問題 |
| P1 | **0** | 原 sealing identity TOCTOU 已修正並有 deterministic regression |
| P2 | **3** | 下列非阻擋限制仍須明示／由 orchestration 控制 |

### P2-1：Hard-kill crash window

`qa_validated_html_report.py:606-615` 在 directory rename 成 canonical path 後，才建立 sidecar、設唯讀、fsync、post verify。普通 exception 會被封存；但 `SIGKILL`、power loss 或 kernel crash 可能留下**無法通過 verify-only**但占住 canonical path 的不完整 directory。

- 不會形成 false PASS，因 exact-five、mode、sidecar 或 post-verification 條件不成立。
- Runbook 必須將這種 canonical invalid directory搬至 archive／failed-publication，再以新 output path 或經審核的恢復流程重跑；不得直接刪除。

### P2-2：Frozen membership 由 presentation contract 證明

QA receipt 會綁定「實際執行的 producer path／SHA／mode／nlink」，但 QA script 本身不讀 `presentation_snapshot_manifest.json`，因此不單獨證明該 producer 是 presentation snapshot allowlist 的成員。

正式流程必須同時滿足：

1. QA 前 `freeze_presentation_snapshot.py verify-only` 通過；
2. 只執行 `$PRESENTATION/code_snapshot/qa_validated_html_report.py`；
3. QA artifact 以**同一 frozen entrypoint**執行 `verify-only`；
4. QA 後再次跑 presentation `verify-only`。

這是刻意的分層契約，不應把單一 QA receipt 說成完整 presentation-manifest attestation。

### P2-3：兩個 content-semantic check 可再收緊

這不影響 artifact path/SHA provenance，但會影響「內容 QA」的字面強度：

- `qa_validated_html_report.py:438-454` 的 `no_external_links` 只把 `http/https`列為 external；protocol-relative `//host` 或 `mailto:`目前不會列入 external。
- `qa_validated_html_report.py:761-765` 的 final status 使用 substring；理論上 `NOT FINAL`也含 `FINAL`。正式 builder 的固定 ribbon 不會使用此字串，但未來可改成 exact token／data attribute。

建議在正式 HTML builder 穩定後新增兩個 content-QA regression；它們不是本次 provenance GO 的 blocker。

## 9. Release gate

**Verdict：GO（QA provenance only）**，前提為：

1. presentation freeze 必須收錄 QA producer SHA `cf66a8cf...d4e16`；
2. formal runbook 不可只 `jq` receipt，必須使用同一 frozen QA entrypoint 跑 `verify-only`；
3. post-QA presentation snapshot 必須再 verify；
4. 正式 QA output path 先確認不存在，不得使用 overwrite；
5. 本文件的 `/tmp` smoke 不得被加入正式 evidence manifest 或引用為研究數據。

本稽核沒有產生、修改或宣稱任何全量 M2 數字。
