<!--
建立時間: 2026-07-22
目標: 保存 schema-recovery authority 簽署前四輪審查、拒絕理由、修補與最終三方判決
處理範圍: authority validator、Vrec、Rrec、Crec；不含科學資料重算
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/pre-decision-audit-recovery.md
-->

# Tumor-REF Schema-Recovery 預簽審查與修補紀錄

## 結論

前三組依序為 **REJECT / 主動撤回 / REJECT**，均未建立 authority、未退休 recovery
private key，也未寫入 19 個 recovery/downstream 輸出槽。第四組已取得兩個內部 agent 與
External Claude Opus 的獨立 **APPROVE**，三者皆為 High `0`、Medium `0`、unresolved `0`；
符合 `GO-to-sign` 門檻。科學資料、canonical receipt 與生物 claim ceiling 均未改動。

## 第一次 Frozen Source

| Role | SHA-256 | Reviewer 結果 |
|---|---|---|
| validator | `b93566dd49fb530f3ba48b9d710a9ed6013be83d194f1d16cbda432c743c3abd` | REJECT |
| Vrec | `77701c30f3929a844af28f95b2b8ea60cf0e87e930929d8fbd8afaff7d3afe18` | REJECT |
| Rrec | `42b6b0995d53001ca72e47c0d8a2b83645ff27d37e1249b7c84cf122e1bc9632` | REJECT |
| Crec | `369eb9462fc1f95315dc4732a1354210e2453a2e26a368ce8266a7e3510870ca` | REJECT |

Reviewer identities：Mendel `019f8ad7-5450-7ac3-a18f-14f312e2c6f9`；Nash
`019f8ad7-bf4e-7891-bd32-00aea77e3732`。舊 external session
`2cb4daf4-29f4-473f-b2c8-0df3b5ba9ab9` 在 source 被修補後已失去適用性並中止，不納入
正式 verdict。

## Findings 與修補

| Severity | Finding | 修補與驗證 |
|---|---|---|
| High | Vrec evidence 沒有 recovery validator artifact，Rrec 卻強制 relation 存在 | validator `_build_validation_evidence()` 新增 `authority_validator`；測試實際走 Rrec recursive relation discovery |
| High | Crec terminal verifier 把 legacy signed source/path 與 active recovery source/path 直接比較，必然在簽署後失敗 | 新增 `validate_signed_terminal_precommit_contract()`；legacy 只用 `authorized_*`，active source 另由 recovery authority `runtime_source` 綁定 |
| Medium | review validator 接受任意非空 probe command | 改為 strict equality `READONLY_PROBE_COMMAND`；測試證明 `/usr/bin/true` 被拒絕 |
| Medium | 原測試沒有覆蓋錯誤 terminal branch | 同一 regression test 驗證正確 old/active 分層可通過，將 legacy source 冒充 active source 必須 fail closed |

## 第二次 Frozen Source（主 Agent 自審後 Superseded）

| Role | Size | SHA-256 |
|---|---:|---|
| validator | 25,651 | `f79c6321eee7efb0748bb73e9689ea0a39bf7ac28345275c2114838f81191198` |
| Vrec | 125,470 | `ec4197cf58655fdbed0af363d5a6bb0cb8589c984b0943341d6044f234801b9c` |
| Rrec | 126,273 | `0fc1721029eac265112cf89875e5cd30680c7cae7500072322ea50fafab0532f` |
| Crec | 290,849 | `ab72c308777e16565dadd469c396e27c68bdc410fd9cb9dae7ae1854a7bbd139` |

第二組在 reviewer 形成 verdict 前由主 agent 主動撤回：validator 雖已要求 exact list，
但該 list 只含 Python argv，尚未把 `/usr/bin/env -i` 與 clean-room environment tokens
納入同一 contract。兩個內部 reviewer 均收到停止指令且未出具 verdict；external session
`17e187cd-2ccb-444a-bbbd-252c25c9c717` 同樣中止，不納入 authority。

## 第三次 Frozen Source（Formal REJECT）

| Role | Size | SHA-256 |
|---|---:|---|
| validator | 26,003 | `208e3e598ee9b2445b8aebd443a174fb5f18fc50ac09a7b6d46cb5e2bdc83351` |
| Vrec | 125,470 | `def82760ce80cc41c3c696397ac2eccfb468b6f3525d95dde82ccdcf62056b42` |
| Rrec | 126,273 | `6ad5e79e939ded1849e1bf59cf3db2e69c07b90ac5b8d1e928ae6c95b15c8035` |
| Crec | 290,849 | `e5c685aba84673eee0d5ee5682e40800b904e77d5fbe3fe6f2e3629b7e92d31d` |

第三組新增的唯一 production contract 變化，是將 review evidence 的 command 精確固定為
`/usr/bin/env -i`、全部環境 token、Python argv 與 probe path。第三組再次取得
`10 passed`、probe exit `0`、`no_output_writes=true`。

Mendel reviewer `019f8ae9-e987-7b30-a4c1-741b7d3bb89a` 判定 High `0`、Medium
`2`、REJECT：第一，Python container equality 仍允許 nested `False == 0`，JSON loader
也接受 duplicate key；第二，retired private key 只有 pathname `stat`，沒有 FD lease 與
return 前終局重檢。Nash 與 external session 在 source supersede 後停止，未出具 verdict。

## 第四次 Frozen Source（Final APPROVE）

| Role | Size | SHA-256 |
|---|---:|---|
| validator | 28,163 | `1dfb3f3efd62068535788a24bf711cabb1ff425a7b09458dcbd182193cb8ac89` |
| Vrec | 125,470 | `f1b29adf1908d58dec4520dba816036bdf27574ceab8dc3c5603c43b63e053c9` |
| Rrec | 126,273 | `dee756a0ca55884ca88260b3b5f1af71b397b26952a80ad1813855d03733ed3c` |
| Crec | 290,849 | `81175713f1a124d4d60347581723bd8ae8f79f3e4cf9539e14edd1354dbab94d` |

第四組加入遞迴 type-strict equality、duplicate-key/non-finite JSON 拒絕，以及
`O_PATH|O_NOFOLLOW` retired-key metadata FD lease。Regression test 以 monkeypatch 模擬 key
mode drift並確認終局 `_require_leases()` fail closed；沒有 chmod 或修改真實 key。實際
`10 passed in 0.34s`，probe exit `0`、19 slots absent，ceremony dry preflight亦
`pass=true/no_output_writes=true`。

### 三方正式判決

| Reviewer | Identity | Verdict | High | Medium | Unresolved |
|---|---|---|---:|---:|---:|
| Mendel | `019f8af7-757a-7c22-905c-b416760f7e34` | APPROVE | 0 | 0 | 0 |
| Nash | `019f8af7-a9bd-7fd0-b6e5-c88ec1fe41b7` | APPROVE | 0 | 0 | 0 |
| External Claude Opus 4.8 | `783f2548-4de3-4c68-a3f5-eec0d0fb0626` | APPROVE | 0 | 0 | 0 |

三位 reviewer 都重跑 exact clean-room probe，確認 exit `0`、`pass=true`、
`no_output_writes=true`、19 slots absent、10/10 regression tests、historical eight-key pilot
與 Vrec/Rrec/Crec recursive checks 全部通過。外部審查耗時 984,386 ms、37 turns，且沒有
permission denial。

External Claude Opus 記錄一個 **Low、non-exploitable** 的 defensive-depth observation：
Python `json.loads(..., parse_constant=...)` 不會攔截 `1e400` 這類 overflow decimal，可能先
解析成 `float('inf')`。此狀況不能產生錯誤 APPROVE，因所有 JSON-derived trusted numeric
不是要求 exact `int`，就是進入 recursive type-strict equality；`float` 無法冒充
`int`、`str` 或 `bool`。Reviewer 明確判定不需改變行為，無 unresolved condition。

第四組正式簽署門檻已滿足：三方判決均針對同一組 SHA，且簽署前最後一次主 agent probe
仍為 exit `0`、`pass=true`、`no_output_writes=true`。下一步可由 frozen ceremony builder
建立三份 review evidence、簽署 supplemental authority，並退休 one-time private key。

## Step → Verify

1. 執行 regression tests → 驗證：`10 passed in 0.33s`，exit `0`。
2. 執行 clean-room source probe → 驗證：exit `0`、`pass=true`、
   `no_output_writes=true`、19 slots absent、Crec 99 functions 自綁定通過。
3. 啟動 fresh 三方 review → 驗證：三方已對第四次 SHA 給出 High `0`、Medium `0`、
   unresolved `0`，全部 `APPROVE`。
4. 發布並驗證 recovery authority → 驗證：三份 review、authority JSON 與 detached signature
   原子建立；OpenSSL verify exit `0`；one-time private key mode 成為 `0000`。

Exact probe：

```bash
env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v1.py
```

本輪只修 release-engineering trust contract；未改動 469,849 位點資料、methyl/co-occurrence
計算、canonical tumor-REF receipt bytes 或生物學 claim ceiling。

## Post-Sign Execution Finding

Authority v1 於三方通過後正式簽署，OpenSSL 與三個 runtime-role validation 均 PASS；Vrec
也完成純驗證與 receipt 建立。Rrec v1 隨後在零 replay output 狀態揭露一個 formal review
未捕捉的必然矛盾：generic relation walker 將 signed historical `during_execution=0o664`
與 signed/live `current=0o444` 同時要求等於目前 inode mode。錯誤為
`Verification receipt relation drift: .../focal_alt_cluster_lib.py field=mode`。

Rrec v1 log/receipt 皆 absent，因此沒有 partial replay artifact。不得修改 v1 frozen source、
authority、reviews 或 retired key；後續採 append-only Rrec/Crec/authority v2，並重新取得
三方 review。此 finding 不改變前三輪與第四輪預簽判決的歷史事實，但證明 source-level
probe 尚未涵蓋 signed authorization payload 的 full runtime traversal；v2 必須新增該正例
與所有 context/mode/schema 負向測試。
