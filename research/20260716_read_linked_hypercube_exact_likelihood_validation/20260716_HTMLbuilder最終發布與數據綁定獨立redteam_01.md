<!--
建立時間: 2026-07-16 18:20 +08:00
目標: 獨立 red-team 最新教授版 HTML builder 的數據分母、S13/S16 provenance、FINAL/PARTIAL gate 與原子發布
處理範圍: Task Type B comprehensive validation；程式、合成 authenticated fixture 與 /tmp adversarial probes；未讀 BAM、未宣稱正式生物數字
服務目標: G4 多樣本一致性與 reproducibility；G5 可由外部重算的證據鏈
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_validated_html_report.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_validated_html_report.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/fixtures/validated_html_report_bundle.json
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_final_numeric_summary.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_presentation_snapshot.py
狀態: independent_red_team_complete; GO under documented non-hostile same-UID and same-frozen execution contract
-->

# HTML builder 最終發布與數據綁定獨立 red-team

## TL;DR

**GO：P0=0、P1=0、P2=3；最新 builder 的數據分母、8 overall＋56 dataset 分層、S13/S16 綁定、FINAL/PARTIAL 隔離與 no-replace 發布均通過（影響：高，信心：高）。**

用 PREP：

- **Point**：可進入正式 same-frozen presentation 執行與下一層整合驗證。
- **Reason**：37/37 builder tests 與 25/25 numeric/presentation-freezer tests PASS；另以十個隔離的 `/tmp` 案例驗證 denominator、tamper、hardlink、mode、parent、競爭發布、overwrite、partial 與 remote guard。
- **Evidence**：合成 FINAL 的 `T_EQ_1_V_EQ_1=0/49`、`single-only=7/35`、`exact-topology unique=21/49` 均以正確母體顯示；兩張 S16 ledger 都是 64 data rows＝8 overall＋56 dataset。
- **Point**：這是 builder/provenance 的 GO，不是正式 M2 生物數字的 GO；上述 0、7、21、35、49 都是 authenticated synthetic fixture 數值。

## 1. 任務邊界與 Step → Verify

| 項目 | 本次定義 |
|---|---|
| Task type | **B — Comprehensive validation** |
| 服務目標 | **G4 / G5** |
| 可支持 | Builder 的分母語意、來源綁定、fail-closed 與發布安全 |
| 不可支持 | 正式 7 datasets×chr1–22 的生物結果、真實 clone 數或唯一演化拓撲 |
| 威脅模型 | 現有 receipt 明示的 non-hostile same-UID；正式流程另須 presentation snapshot same-frozen `verify-only` |

1. 固定輸入 SHA → 驗證：測試前後 SHA 不變。
2. 重建合成 FINAL → 驗證：兩張表各有 8 overall＋56 dataset data rows。
3. 逐格驗 denominator → 驗證：0/49、7/35、21/49 不被錯寫為 0/0、7/49 或 21/21。
4. 攻擊 S16/S13 → 驗證：合法不同路徑同 bytes 通過；重簽 drift、hardlink、mode 與 parent drift 均拒絕。
5. 攻擊發布邊界 → 驗證：競爭者不被覆寫、0444 staging 保留、overwrite 舊 bytes 封存、PARTIAL 不建立 FINAL。
6. 攻擊輸出內容 → 驗證：HTML 特殊字元逃逸，remote URL payload 在 publish 前拒絕。

## 2. 審核輸入與 SHA-256

工作目錄：`/big7_disk/liaoyoyo2001/InterSubMod`

| 角色 | InterSubMod 路徑 | SHA-256 |
|---|---|---|
| HTML builder | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_validated_html_report.py` | `613d57a8ba956f625221cf1231625424e3b011c291858fd92a3bd6e7a8467e18` |
| Builder tests | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_validated_html_report.py` | `3891dfb5f7ac590e2655fe1461dd44f69e4541b4cf0ad984dbe80fe6e7af6cf4` |
| Builder fixture | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/fixtures/validated_html_report_bundle.json` | `56aa7c9f7eeb210e0bc2faed7ffcda22a5646389de3c4db5d4281a4d5b79298a` |
| S16 numeric producer | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_final_numeric_summary.py` | `8952ccb17a0e2514621ef99110b9bc589f084398855d95e19eee627f40d6a4cd` |
| Numeric fixture | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/fixtures/final_numeric_summary_bundle.json` | `044e48524c4a00311378e86b83c1a335b3240f6740b574b5d96a6a2e4980da45` |
| Presentation freezer | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_presentation_snapshot.py` | `639898a882b5d72f9c58ad95701a5afdfa22c0ad37f1aebc8748320d2efbd625` |
| Freezer tests | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_presentation_snapshot.py` | `af77886f98cf64e812b1879b49ebfe96e574396562dc01ebd0e8d8a6a36816b5` |

上述七個檔案在本稽核前後 SHA 相同；本 red-team 沒有修改 builder、tests、fixtures、numeric producer 或 freezer。

## 3. 數據與分母驗證

### 3.1 兩張 ledger 的完整分層

每張表使用相同 Cartesian scope：

```text
2 HP bases × 4 bridge thresholds = 8 OVERALL rows
7 datasets × 2 HP bases × 4 bridge thresholds = 56 dataset rows
合計 = 64 data rows（table header 不計）
```

獨立解析 rendered HTML 的結果：

```text
Extraction component 與 k=1/k>1 data_rows=64 overall=8 dataset=56
完整 T/V、effective route、h*、ranking 與 topology ledger data_rows=64 overall=8 dataset=56
```

另外逐 `HP×threshold` 驗證七 dataset 的 effective route、T/V 三分桶及四個 coarse-unique topology classes 加總，全部精確等於 overall cell；沒有 aggregate cancellation。

### 3.2 三個指定顯示值

| 顯示 | 正確 denominator | Rendered 結果 | 避免的錯誤 |
|---|---:|---:|---|
| `T_EQ_1_V_EQ_1` | 49 solver-complete units | `0 / 49 = 0.00%` | 0 不可因 numerator=0 被顯示為 `0/0 N/A` |
| `single-only` | 35 coarse-class-unique units | `7 / 35 = 20.00%` | 不可用全部 49 topology-evaluated units |
| exact-topology proven unique | 49 topology-evaluated units | `21 / 49 = 42.86%` | 不可自除成 `21/21 = 100%` |

這三個值是 synthetic fixture 的語意回歸測試，不是正式研究數字；它們證明的是「分母路由正確」。

### 3.3 程式證據位置

| 契約 | 程式位置 | 判讀 |
|---|---|---|
| Scope loop | `build_validated_html_report.py:4549-4622` | `(OVERALL, *7 datasets)×2 HP×4 thresholds` |
| T/V bucket denominator | `build_validated_html_report.py:4666-4679` | total 與 relative 都用 solver-complete denominator |
| Coarse class denominator | `build_validated_html_report.py:4572-4595, 4681-4694` | unique 與 ambiguous 各用自己的 partition denominator |
| Exact topology denominator | `build_validated_html_report.py:4709` | 用 topology-evaluated units，不用 exact-unique numerator |
| Per-stratum S16 reconciliation | `build_validated_html_report.py:2740-2910` | candidate stream、raw ranking cell 與 S16 cell 逐格核對 |

## 4. S16 numeric producer 綁定

正向案例刻意讓 receipt 記錄的 live producer 位於 fixture 目錄，而不是 HTML builder 旁的 sibling：

```text
recorded path != frozen sibling path
recorded size == sibling size
recorded SHA-256 == sibling SHA-256
recorded bytes == sibling bytes
結果：FINAL PASS
```

負向案例在 recorded live producer 追加 bytes，並同步更新 producer `size_bytes`、`sha256` 以及 numeric-summary sidecar，讓「receipt↔live」仍自洽。Builder 仍因「recorded live bytes≠sibling frozen bytes」拒絕且不建立輸出。

程式位置：`build_validated_html_report.py:3123-3178`。

**結論**：不同路徑是允許的；真正的 authority 是 size＋SHA＋bytes identity。只重簽 drift 無法穿透。

## 5. S13 實體 release contract

### 5.1 正向實體條件

- manifest、固定相鄰 sidecar、canonical copy 與 11 snapshots：regular non-symlink、`0444`、`st_nlink=1`。
- 從每個檔案 parent 到 release root：regular non-symlink directories、`0555`。
- manifest physical SHA、semantic SHA、sidecar exact line、snapshot role/path/size/SHA 均與 S13 binding 一致。

### 5.2 負向攻擊

| Tamper | 即使什麼仍一致 | 實際結果 |
|---|---|---|
| 對 `snapshot/extractor` 建外部 hardlink | bytes/SHA 不變 | 拒絕：`physical st_nlink 不是 1` |
| manifest 改為 `0644` | bytes/SHA 不變 | 拒絕：`physical mode 不是 0444` |
| snapshot parent 改為 `0755` | snapshot bytes/mode 不變 | 拒絕：`parent contract mode 不是 0555` |

程式位置：`build_validated_html_report.py:1174-1325`；presentation freezer 另由 11/11 tests 驗證建立與 verify-only 契約。

## 6. 原子發布、overwrite 與 PARTIAL

### 6.1 No-replace race

在 temp 已 fsync/chmod 後、正式 rename 前注入競爭者：

```text
destination = "competing writer\n"                    → 原 bytes 保留
renameat2(RENAME_NOREPLACE)                            → FileExistsError
staging = .race.html.3v3d8aa4.tmp                      → 保留
staging mode=0444, st_nlink=1                          → PASS
```

保留證據：`/tmp/intersubmod_html_builder_redteam_20260716_v3/07_publish_race/`。

### 6.2 Explicit overwrite

原檔不刪除、不截斷，先以 no-replace rename 封存：

```text
old bytes SHA-256 = 0306d981ed8e5eb35f7eec731f4dfcfbf337df06efd18598e97d4815b4e28388
archive count = 1
archive bytes = original bytes
new report mode = 0444
```

### 6.3 FINAL/PARTIAL 隔離

移除 full extraction/ranking sources 後，即使指定 FINAL 名稱並使用 `allow_partial=True`：

- 指定 FINAL 不存在；
- 實際輸出自動改名 `.partial-preview.html`；
- 內文固定含 `PARTIAL PREVIEW · NOT VALIDATION EVIDENCE`；
- CLI contract 的 `all_pass=false`、`final_ready=false` 已由正式 unit test 覆蓋。

### 6.4 No-delete 靜態掃描

對 builder 搜尋 `unlink/remove/rmtree/rmdir/replace` 等 destructive filesystem API，**無輸出**。發布只使用 exclusive temp、no-replace rename 與 preserving archive。

## 7. Escaping 與 standalone remote guard

獨立輸入：

```text
"><img src=x onerror=alert(1)>&
```

結果：

```text
HTML text = &quot;&gt;&lt;img src=x onerror=alert(1)&gt;&amp;
local href = %22%3E%3Cimg%20src%3Dx%20onerror%3Dalert%281%29%3E%26
```

另把已重簽 S16 的 `reason` 改為 `https://evil.invalid/...`；builder 在建立 temp 前拒絕 `standalone report unexpectedly contains remote URL or script`，destination 不存在。

結論：目前動態文字走 `html.escape(..., quote=True)`，local href 走 URL percent encoding，且整份 document 另有 standalone guard。

## 8. 分級發現

### P0 — 0

沒有會造成錯誤 FINAL 數字、錯誤 authority 或競爭覆寫的 critical 問題。

### P1 — 0

沒有阻擋本版進入正式 same-frozen integration 的 high-severity 問題。

### P2 — 3（不阻擋 GO）

1. **S13 stable read 仍是 pathname `lstat→read_bytes→lstat`**：沒有以 `open(O_NOFOLLOW)`＋同一 fd `fstat` 完成整段讀取，且 stability tuple 未含 `st_nlink`。在文件宣告的 non-hostile same-UID 模型內，持續 hardlink/mode tamper 都會被擋；若未來要防同 UID 主動競爭，建議改成 fd-bound read。
2. **`overwrite=True` 的舊 artifact 在新 document render/stage 前先封存，且封存沿用舊 mode**：若後續 render 失敗，指定 canonical path 暫時缺席但舊 bytes 仍在 `.superseded.*`；本 probe 的舊檔為 `0664`，封存後仍可寫。建議先完成新 temp 驗證，再封存／發布，並依 artifact policy 將 archive 封成 `0444` 後 fsync。
3. **發布邊界是 Linux `renameat2` 且只驗 immediate output parent**：不支援的系統會 fail-closed，不會不安全 fallback；但 portability 與 hostile ancestor-symlink 防禦尚未提供。若要跨平台或提升敵對模型，需加入安全 no-replace backend 與 output ancestor chain 驗證。

這三項不改變本次正確性結論；它們是 threat-model / availability / portability 的 defense-in-depth backlog。

## 9. 實際執行命令與輸出

### 9.1 Builder 全測試

輸入：builder、37-test suite、validated fixture。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 -m py_compile \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_validated_html_report.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_validated_html_report.py

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 -m unittest -v \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_validated_html_report.py
```

輸出：terminal；未建立 repo 內數據輸出。

```text
Ran 37 tests in 104.348s
OK
```

### 9.2 Numeric＋presentation freezer 回歸

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 -m unittest -v \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_final_numeric_summary.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_presentation_snapshot.py
```

```text
Ran 25 tests in 0.354s
OK
```

### 9.3 `/tmp` adversarial probe

輸入：同一 test materializer 與 fixture；每個案例使用獨立子目錄。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 - <<'PY'
# 建立 10 個隔離 synthetic authenticated bundles：
# metrics、S16 resign-drift、S13 baseline/hardlink/mode/parent、
# publication race、overwrite archive、PARTIAL、escaping/remote。
# 每一案例以獨立 assertions 比對 expected pass/fail 與 bytes/mode/nlink。
PY
```

輸出根：`/tmp/intersubmod_html_builder_redteam_20260716_v3/`。

```text
metrics: final_ready=true; rows=64; overall=8; dataset=56; target=0/49,7/35,21/49
s16: different_path_same_bytes_pass=true; drift_and_resign_rejected=true
s13: baseline_0444_nlink1_0555=true; hardlink/mode/parent rejected
race: destination_preserved=true; staging_mode=0444; staging_nlink=1
overwrite: archive bytes preserved; new_mode=0444
partial: final_absent=true; ribbon_present=true
escaping_remote: escaped/encoded=true; remote payload rejected; output_absent=true
```

合成 FINAL HTML SHA-256：

```text
dc1b07826f7980a41ff7e11209f130d29ba4c57c74780716c763e6f76e4e8084
```

## 10. 最終 verdict

**GO — 最新 HTML builder 可進入正式 presentation snapshot same-frozen 執行與整合。** 必須維持以下邊界：

- 正式 artifact 只能餵真實 full-scope S13/S16 及全部 final receipts；
- 先執行 presentation freezer，正式讀取前後走 same-frozen `verify-only`；
- fixture 的 0/49、7/35、21/49 只作 denominator regression，不可貼成研究結果；
- P2 三項列入 defense-in-depth backlog，不阻擋這版。
