<!--
建立時間: 2026-07-16 17:03:41 CST (+0800)
目標: 獨立 red-team 驗證 final numeric summary 的 T/V 三分桶、分層守恆與 terminal 綁定，以及 M2 release freezer 的 no-replace、失敗封存與 st_nlink=1 保證。
處理範圍: Task Type B comprehensive-validation 的程式與合成證據審核；未讀取 BAM，未啟動正式 M2 RUN_ROOT，未產生可作生物結論的正式數據。
服務目標: G4 多樣本一致性與 reproducibility；G5 可被外部驗證的業界級證據鏈。
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_final_numeric_summary.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_final_numeric_summary.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/fixtures/final_numeric_summary_bundle.json
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_m2_release_contract.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_m2_release_contract.py
-->

# Numeric summary 與 release freezer 修補獨立 red-team

## TL;DR

**獨立 red-team 通過：P0=0、P1=0；未發現會阻擋正式 M2 執行的正確性問題（影響：高，信心：高）。**

用 PREP 敘述：

- **Point**：T/V 三互斥桶確實由 authenticated candidate rows 重算，並在 overall、HP×threshold、dataset×HP×threshold 三層守恆；freezer 的 production path 無刪除式 API，且 no-replace 競爭失敗時不會搬走競爭者。
- **Reason**：每個 complete unit 使用 `V=len(candidate rows)` 與 `T=Σparent_choice_count`，強制 `T≥V≥1`，後續再逐分層精確對齊 terminal ranking 的 complete units/T/V。
- **Evidence**：34/34 targeted tests PASS；三種 numeric tamper 均 fail-closed；250 trials/9,873 units fuzz PASS；40 rounds×8 contenders no-replace race 每輪恰好一個 winner。
- **Point**：在本文明示的威脅模型下，這五個檔案可進入下一層整合驗證。

## 1. 任務分類與邊界

| 項目 | 定義 |
|---|---|
| Task type | **B — Comprehensive validation** |
| 服務目標 | **G4 / G5** |
| 審核單位 | 程式碼、測試、synthetic fixture、臨時合成 receipt bundle |
| 本次結論可支持 | Numeric/freezer 實作的守恆、fail-closed、no-replace 與可重現性 |
| 本次結論不支持 | 正式 7 datasets×chr1–22 的數值結果、生物 clone 數、真實拓撲或正式效能包絡 |

### Step → Verify

1. 固定輸入檔案與 SHA → 驗證：五個 SHA 與委派版本一致。
2. 追蹤 T/V 分桶與 terminal 綁定 → 驗證：指出 row→stratum→overall 的確切程式位置。
3. 建立最小 numeric tamper → 驗證：合法 bundle 通過；T/V 或分層交換必須被拒絕。
4. 建立隨機守恆測試 → 驗證：250 trials 的三層數據全部與獨立 Counter 一致。
5. 審查 freezer publish/failure path → 驗證：production source 無刪除 API，競爭者不被搬走，成功與封存檔案 `st_nlink=1`。

## 2. 輸入證據與 SHA-256

工作目錄：`/big7_disk/liaoyoyo2001/InterSubMod`

| 角色 | 輸入路徑 | SHA-256 |
|---|---|---|
| Numeric producer | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_final_numeric_summary.py` | `8952ccb17a0e2514621ef99110b9bc589f084398855d95e19eee627f40d6a4cd` |
| Numeric tests | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_final_numeric_summary.py` | `2ac2c7ecdab7bdd956aba13b644f2b396f1d1f67761365bb3122b9d333a37ce0` |
| Numeric fixture | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/fixtures/final_numeric_summary_bundle.json` | `044e48524c4a00311378e86b83c1a335b3240f6740b574b5d96a6a2e4980da45` |
| Release freezer | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_m2_release_contract.py` | `c734c1ed2142f6e2baed0155320cdaba8925e2d6b76965e60bd9e42e1bf4f7f1` |
| Freezer tests | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_m2_release_contract.py` | `83240880a2da104e7aba500fc788a5a25fc28fea6acb820ebc2eefda6429cbfc` |

執行前與執行後重算 SHA 一致；本 red-team 未修改上述五個檔案。

### SHA 命令

```bash
sha256sum \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_final_numeric_summary.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_final_numeric_summary.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/fixtures/final_numeric_summary_bundle.json \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_m2_release_contract.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_m2_release_contract.py
```

### Actual output

```text
8952ccb17a0e2514621ef99110b9bc589f084398855d95e19eee627f40d6a4cd  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_final_numeric_summary.py
2ac2c7ecdab7bdd956aba13b644f2b396f1d1f67761365bb3122b9d333a37ce0  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_final_numeric_summary.py
044e48524c4a00311378e86b83c1a335b3240f6740b574b5d96a6a2e4980da45  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/fixtures/final_numeric_summary_bundle.json
c734c1ed2142f6e2baed0155320cdaba8925e2d6b76965e60bd9e42e1bf4f7f1  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_m2_release_contract.py
83240880a2da104e7aba500fc788a5a25fc28fea6acb820ebc2eefda6429cbfc  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_m2_release_contract.py
```

## 3. Numeric T/V 實作審核

### 3.1 單位、公式與三分桶

對每個 solver-complete candidate-table unit：

\[
V_u = \lvert\text{candidate rows in unit }u\rvert
\]

\[
T_u = \sum_{r\in u}\texttt{parent\_choice\_count}_r
\]

程式強制每個 `parent_choice_count≥1`，因此 `T≥V≥1`，只可能落在三個互斥且穷盡的狀態：

| Bucket | 程式定義 | 解釋邊界 |
|---|---|---|
| `T_EQ_1_V_EQ_1` | `T=1,V=1` | 一個最佳 vertex set，且其 parent-edge assignment 唯一 |
| `T_GT_1_V_EQ_1` | `T>1,V=1` | vertex set 唯一，但 parent-edge assignment 不唯一 |
| `T_GT_1_V_GT_1` | `T>1,V>1` | 最佳 vertex set 不唯一；必然 `T≥V>1` |

這是「候選結構」分類，不是 clone 數的生物聲明。

### 3.2 程式證據鏈

| 證據點 | 程式位置 | 驗證內容 |
|---|---|---|
| Row-level T/V 重算 | `build_final_numeric_summary.py:790-810` | `V=len(rows)`、`T=sum(parent_choice_count)`、`T≥V≥1`、三分桶 |
| Unit 累積 | `build_final_numeric_summary.py:834-840` | 同時累積至 `(dataset,basis,threshold)` 與 overall |
| Bucket/分母/比例 | `build_final_numeric_summary.py:877-939` | exact keys、非負整數、加總守恆、ratio 重算 |
| Aggregate invariant | `build_final_numeric_summary.py:942-982` | aggregate `T≥V≥n_units` 與 outcome/h*/topology 守恆 |
| Dataset×HP×threshold | `build_final_numeric_summary.py:1574-1619` | candidate complete units/T/V/outcome/topology 精確對齊 authenticated ranking cells |
| Overall terminal binding | `build_final_numeric_summary.py:1620-1640` | candidate overall units/T/V 對齊 terminal child reaggregation |
| HP×threshold 整體輸出 | `build_final_numeric_summary.py:1714-1724` | 合併各 dataset 已累積 candidate strata 後再建 presentation cell |

### 3.3 Targeted test 命令

```bash
/usr/bin/time -f 'elapsed=%e sec maxrss=%M KiB exit=%x' \
  /bip7_disk/liaoyoyo2001/miniconda3/bin/python3 -m unittest -v \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_final_numeric_summary.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_m2_release_contract.py
```

### Actual output

```text
test_authenticated_bundle_derives_component_hstar_and_tied_topology ... ok
test_candidate_aggregate_T_less_than_V_fails_closed ... ok
test_candidate_byte_tamper_fails_closed ... ok
test_tree_vertex_partition_covers_all_three_mutually_exclusive_states ... ok
test_tree_vertex_partition_tampering_fails_closed ... ok
...
test_external_hardlink_to_snapshot_is_rejected ... ok
test_external_hardlink_to_source_is_rejected ... ok
test_publish_boundary_source_mutation_is_rejected ... ok

----------------------------------------------------------------------
Ran 34 tests in 2.951s

OK
elapsed=3.18 sec maxrss=25496 KiB exit=0
```

## 4. Numeric tamper 負測

所有 tamper 都在 Python `TemporaryDirectory` 內，使用
`test_build_final_numeric_summary.py::build_bundle()` 建立 synthetic authenticated bundle，修改後重算 synthetic receipt/sidecar SHA，以避免只測到最表層的 byte-tamper。最後直接呼叫：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 - <<'PY'
# 1. import test_build_final_numeric_summary.py and build_final_numeric_summary.py
# 2. build_bundle(TemporaryDirectory)
# 3. rewrite candidate TSV.GZ, recompute physical + semantic SHA
# 4. update synthetic terminal/final candidate identities and SHA sidecars
# 5. call build_summary(..., require_canonical_scope=False)
# 6. require SummaryError; print before/after T,V,bucket and exact rejection
PY
```

> 上述是本次實際 inline-Python 命令的執行介面；因臨時 tamper 程式不是生產程式，本次不另建立可執行檔。其輸入、變更、獨立預期與 actual exception 完整記錄如下。

### 4.1 Terminal T tamper

| 層級 | Tamper 前 | Tamper 後 |
|---|---:|---:|
| HP1 | units=1, T=1, V=1, `T=1,V=1`=1 | units=1, T=2, V=1, `T>1,V=1`=1 |
| HP2 | units=1, T=4, V=2 | units=1, T=4, V=2 |
| Overall | units=2, T=5, V=3 | units=2, T=6, V=3 |

Terminal ranking 的 HP1 T 故意保持 1；candidate table 為 2。

Actual output：

```text
REJECTED: candidate T mismatch: SYNTH/PS_HP1/3
```

### 4.2 HP stratum swap：overall T/V 不變

| 層級 | Swap 前 | Swap 後 |
|---|---:|---:|
| HP1 | units=1, T=1, V=1 | units=1, T=4, V=2 |
| HP2 | units=1, T=4, V=2 | units=1, T=1, V=1 |
| Overall | units=2, T=5, V=3 | units=2, T=5, V=3 |

Actual output：

```text
REJECTED: candidate V mismatch: SYNTH/PS_HP1/3
```

結論：不是只核對 overall；HP×threshold 的 cell 也有精確綁定。

### 4.3 Cross-dataset swap：global 不變

獨立將 fixture 扩成 `SYNTH` 與 `ALT` 兩 dataset：

1. 合法的兩 dataset bundle 先執行 `build_summary()` → **PASS**。
2. `PS_HP1/bridge=3` 中：SYNTH T=1，ALT T=2。
3. 將 candidate rows 的 dataset/unit-key 對調，使 global HP1 T=3/V=2 不變。
4. 重算 candidate physical/semantic SHA、terminal ranking SHA 與 final sidecar。
5. 再執行 `build_summary()`。

Actual output：

```json
{
  "valid_two_dataset_bundle": "PASS",
  "pre_dataset_HP1_T": {"ALT": 2, "SYNTH": 1},
  "tampered_candidate_dataset_HP1_T": {"ALT": 1, "SYNTH": 2},
  "overall_all_HP_T_V_after": [11, 6],
  "verdict": "REJECTED: candidate T mismatch: SYNTH/PS_HP1/3"
}
```

結論：dataset×HP×threshold 不會被全體總數相同所掩蓋。

## 5. Numeric 隨機守恆 fuzz

### 命令與參數

```bash
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 - <<'PY'
# deterministic seed
rng = random.Random(20260716)
# 250 trials; each trial randomly creates 1..80 units
# axes: datasets=(D0,D1,D2), bases=(PS_HP1,PS_HP2), thresholds=(3,5)
# V randomly 1..4; every candidate parent_choice_count randomly 1..4
# independently count n_units/T/V/three buckets with collections.Counter
# compare against:
#   finish_candidate_unit()
#   freeze_candidate_accumulator()
#   _merge_candidate_strata()
# assert equality at dataset×HP×threshold, HP×threshold, and overall
PY
```

### Actual output

```text
FUZZ_PASS trials=250 units=9873 strata_levels=dataset×HP×threshold,HP×threshold,overall seed=20260716
```

獨立 Counter 與 production summary 在 9,873 個隨機單位上完全一致；三桶每層的加總均等於 solver-complete unit 分母。

## 6. Freezer production-path 靜態審核

### 命令

```bash
rg -n --pcre2 \
  '(os\.(unlink|remove|removedirs|rmdir|replace)|shutil\.rmtree|Path\.(unlink|rmdir|replace)|\.(unlink|rmdir)\s*\()' \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_m2_release_contract.py \
  || true

rg -n '\.replace\(' \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_m2_release_contract.py
```

### Actual output

第一個 destructive-filesystem API 查詢無輸出。唯三個 `.replace()` 是字串日期格式轉換，不是 filesystem replace：

```text
1033: created_at = dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")
1256: created = dt.datetime.fromisoformat(str(run_manifest["created_at_utc"]).replace("Z", "+00:00"))
1470: "verified_at_utc": dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z"),
```

### 實作行為

| 行為 | 位置 | 審核結論 |
|---|---|---|
| 新檔建立 | `freeze_m2_release_contract.py:615-639` | `O_CREAT|O_EXCL`，新 inode 必須 `st_nlink=1` |
| 原子 publish | `freeze_m2_release_contract.py:659-675` | `renameat2(RENAME_NOREPLACE)`，不覆寫 destination |
| 失敗 tree 封存 | `freeze_m2_release_contract.py:678-694` | 使用 no-replace rename 搬到 `.failed-*`；封存失敗時保留原位 |
| Receipt half-pair 封存 | `freeze_m2_release_contract.py:697-721` | 只處理本程序的 published/temp remnants，不刪除 |
| Contract staging publish | `freeze_m2_release_contract.py:1133-1168` | 最後邊界再驗證後 no-replace publish；失敗則封存 staging |

## 7. Freezer competitor 負測

所有測試均使用同 filesystem 的 `TemporaryDirectory`、記錄 competitor inode 與內容，再以 mock 在 publish boundary 建立競爭目的檔。

| 情境 | Actual verdict | Competitor | 本程序證據 |
|---|---|---|---|
| `_rename_noreplace` destination 已存在 | REJECTED | source 與 destination inode/內容都保留 | 無覆寫 |
| Receipt 第一個 destination 競爭 | REJECTED | competitor 原 inode/內容保留 | 2 個 temp owned files 進 `.failed-publication.*`，均 `nlink=1` |
| Sidecar 第二個 destination 競爭 | REJECTED | sidecar competitor 原 inode/內容保留 | half-published receipt 與 temp sidecar 被封存，均 `nlink=1` |
| Contract root 競爭 | REJECTED | contract competitor sentinel 原 inode/內容保留 | 完整 staging 進 `.failed-staging.*`；manifest `nlink=1` |
| 成功 receipt | PASS | N/A | receipt/sidecar SHA 匹配、mode `0444`、`nlink=1` |

### Actual output

```json
{
  "contract_root_competitor": {
    "archive_count": 1,
    "archived_manifest_nlink": 1,
    "competitor_preserved": true,
    "staging_manifest_archived": true,
    "verdict": "REJECTED: refusing to overwrite release contract root"
  },
  "receipt_first_destination_competitor": {
    "archive_count": 1,
    "archived_all_nlink_1": true,
    "archived_owned_files": 2,
    "competitor_not_archived": true,
    "competitor_preserved": true,
    "sidecar_absent": true
  },
  "receipt_second_destination_competitor": {
    "archive_count": 1,
    "archived_all_nlink_1": true,
    "archived_owned_files": 2,
    "competitor_not_archived": true,
    "competitor_preserved": true,
    "owned_half_receipt_archived": true
  },
  "rename_noreplace_existing_destination": {
    "competitor_preserved": true,
    "source_preserved": true
  },
  "success_receipt_nlink": {
    "receipt_mode": "0444",
    "receipt_nlink": 1,
    "sha_matches": true,
    "sidecar_mode": "0444",
    "sidecar_nlink": 1
  }
}
```

## 8. No-replace concurrent race fuzz

### 完整命令

```bash
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 - <<'PY'
import concurrent.futures, importlib, sys, tempfile, threading
from pathlib import Path

scripts = Path('/big7_disk/liaoyoyo2001/InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts')
sys.path.insert(0, str(scripts))
freezer = importlib.import_module('freeze_m2_release_contract')

rounds = 40
contenders = 8
all_ok = True
for trial in range(rounds):
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        destination = root / 'winner'
        sources = []
        for index in range(contenders):
            source = root / f'source-{index}'
            source.mkdir()
            (source / 'marker').write_text(str(index))
            sources.append(source)
        barrier = threading.Barrier(contenders)

        def race(index):
            barrier.wait()
            try:
                freezer._rename_noreplace(sources[index], destination)
                return ('WIN', index)
            except freezer.ReleaseContractError:
                return ('LOSE', index)

        with concurrent.futures.ThreadPoolExecutor(max_workers=contenders) as executor:
            results = list(executor.map(race, range(contenders)))
        winners = [index for status, index in results if status == 'WIN']
        losers = [index for status, index in results if status == 'LOSE']
        ok = (
            len(winners) == 1
            and len(losers) == contenders - 1
            and (destination / 'marker').read_text() == str(winners[0])
            and all(
                sources[index].is_dir()
                and (sources[index] / 'marker').read_text() == str(index)
                for index in losers
            )
        all_ok &= ok
        if not ok:
            raise SystemExit(f'RACE_FAIL trial={trial} results={results}')

print(
    f'RENAME_NOREPLACE_RACE_PASS rounds={rounds} contenders={contenders} '
    f'exactly_one_winner_each={all_ok} loser_sources_preserved=true'
)
PY
```

### Actual output

```text
RENAME_NOREPLACE_RACE_PASS rounds=40 contenders=8 exactly_one_winner_each=True loser_sources_preserved=true
```

## 9. Severity findings

| Severity | 數量 | 結果 |
|---|---:|---|
| P0 | **0** | 未發現證據污染、覆寫競爭者、完全錯誤分桶或會讓正式結果無效的問題 |
| P1 | **0** | 未發現分層 T/V 可被 overall 掩蓋、terminal 不綁定或 failed archive 搬走競爭者的問題 |
| P2 | **2** | 以下為不阻擋正確性的 regression/diagnostic 改進點 |

### P2-1：持久化 numeric fixture 的 scope 較小

Repo 內現有 fixture 是 `SYNTH×chr1×bridge=3`，未將「兩 dataset 分層互換、global 不變」固化成 regression test。本次臨時 2-dataset red-team 已證明目前實作正確，因此這是未來防回歸建議，不是當前 correctness blocker。

### P2-2：no-replace 衝突訊息太通用

`_rename_noreplace()` 在 receipt/sidecar 目的路徑衝突時，仍回報 `refusing to overwrite release contract root`。行為是正確 fail-closed，但文字會讓診斷者誤以為衝突必然發生在 contract root。

## 10. 威脅模型與不可超譯邊界

本驗證採用正式協定的威脅模型：

- 假設沒有 hostile same-UID actor 在 syscall 之間刪除、置換或 hardlink 本程序已 publish 的 inode。
- 驗證重點是意外並行執行、正常競爭、部分 publish 失敗與 persisted evidence 漂移。
- SHA sidecar 與 receipt 是 persisted protocol authentication，不是對惡意 OS/user ancestry 的 cryptographic attestation。
- 當 destination 原本存在、或在 no-replace boundary 被另一正常執行者先取得時，本次負測已證明 competitor inode 與內容不變。

Numeric 的不可超譯邊界：

- T/V 是 solver-complete 候選結構數，不是 clone 數。
- `h*` 是 minimum extra-state objective，不是 hidden clone 數。
- 本次 red-team 不驗證正式 7-dataset 數值；它只驗證正式數值未來進入時的計算與凍結容器。

## 11. 未建立正式 RUN_ROOT 聲明

本 red-team：

- **未執行** BAM/VCF extraction。
- **未啟動** 7 datasets×chr1–22 正式 ranking。
- **未建立** `m2_frozen_release_v1` 或任何正式 `RUN_ROOT`。
- **未使用** `--ignore-resource-gate`。
- 所有 tamper/fuzz/race 輸出只存在 `TemporaryDirectory`，隨測試結束清理。
- 未將 synthetic fixture 當成 validation evidence，也未寫入 formal evidence ledger。

## 12. 最終判定

**Ready for the next integration-validation stage, with the two P2 items retained as non-blocking follow-up.**

可支持的精確結論：

1. Numeric summary 會從 authenticated candidate rows 重算 T/V 與三互斥桶。
2. Overall、HP×threshold、dataset×HP×threshold 均守恆，且逐 cell 綁定 terminal ranking T/V。
3. 即使 tamper 保持 overall 不變，HP 或 dataset 層的交換仍會 fail-closed。
4. Freezer production path 不使用刪除式 API，publish 使用 `RENAME_NOREPLACE`。
5. 在既定威脅模型下，競爭者不會被失敗封存搬走；成功與封存證據的 `st_nlink=1`。
