<!--
建立時間: 2026-07-22
目標: 在實作 signed-promotion verifier recovery 前完成 pre-decision evidence audit
處理範圍: historical during_execution 8-key schema 與 V/R/C recovery；不改寫任何已簽 artifact
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/verify_tumor_ref_receipt_promotion_v2.py
-->

# Signed Promotion Verifier Recovery Pre-Decision Audit

## §0 Cynefin Domain Gate

**Complex，強制 probe-first。** Promotion、兩份 detached signatures 與 canonical
receipt 已不可逆發布，但第一個 continuation verifier 在 historical metadata schema
邊界 fail-closed。相同 recovery action 尚未曾重複產生可預測結果，不能套用一般
one-line bugfix 流程。

## §1 Observation Completeness

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| V canonical invocation exit `1`，沒有 verification receipt | ✓ | L1 runtime | `verify_tumor_ref_receipt_promotion_v2.py --verify-and-record` stderr |
| Historical focal `during_execution` 只有 8 keys，沒有 `link_count` | ✓ | L1 artifact | `results/tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json:17` |
| P 生成 transition 時只沿用 historical 8-key record，沒有宣稱當時 link count | ✓ | L1 code | `audit_notes/promote_tumor_ref_source_receipt_v2.py:1858` |
| V 將該 historical record 誤套 `STAT_ARTIFACT_KEYS` 9-key live schema | ✓ | L1 code | `audit_notes/verify_tumor_ref_receipt_promotion_v2.py:1668` |
| 已簽 promotion/canonical artifacts 不可被改寫 | ✓ | L1 signature | prepare/promote detached signatures 均經 OpenSSL 驗證 |
| R/C 是否可只換 verifier 而不調整 source/command binding | □ | L2 code review | fresh agents Mill/Dalton + 本 session tracing |

## §2 Credibility Score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | historical observation 與 live full-stat 是不同 schema |
| 觀察支撐 | 20 | JSON、P 生成碼與 V error 三者一致 |
| 機制清晰度 | 20 | exact-key validator 在下游前直接拒絕 |
| 反例風險 | 0 | R/C transitive source/command binding 尚待完整盤點 |
| 所需資源 | 0 | 正式 supplemental authority + 三方 review 可能超過 6 小時 |
| **TOTAL** | **60/100** | 原始門檻可 GO，但 red-team gate 降為 PROBE |

**Falsifier observable**：若「窄幅 schema recovery 足夠」是假設錯誤，修正 historical
8-key validator 並區分 historical authorized source / executing recovery source 後，唯讀
`--verify` 仍會在其他 signed cross-link 失敗，或出現任何 formal/downstream 新產物。

## §3 Assumption Map

| importance | known | unknown |
|---|---|---|
| High | signed promotion artifacts 身分、兩把 promotion key 已退役、downstream 尚未開始 | R/C 所有 transitive path/command binding；recovery authority 是否能避免 circular hash |
| Low | 原始科學輸出與 source receipt bytes 未改變 | 最終 recovery 文件命名與展示細節 |

## §4 Quick Pilot

1. 建立**新路徑** recovery verifier，保留原 V mode-0444 bytes不變，只修 historical
   8-key type validator並把 signed authorization 的舊 verifier 身分與新執行 source 分開。
2. 以 clean environment 執行 recovery verifier `--verify`，禁止 record/replay/downstream。
3. 比對執行前後 formal/downstream slot inventory。

Checkpoint：

- `exit=0`、stdout `pass=true`、新增 artifact `0` → V pilot PASS，再完成 R/C 與
  supplemental signed recovery authority。
- 仍為 signed cross-link error → 保持 PROBE，依失敗點擴充最小修補集合。
- 任何 downstream/output 提前出現 → NO-GO，停止並封存 incident。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---:|
| R/C 對舊 verifier/replayer/continuation 的完整 hard binding 清單 | 決定最小修補集合 | 1-2h | P0 |
| Supplemental authority source/key/review schema | 新 source 無法可信取代舊 source | 2-4h | P0 |
| Recovery fault-injection + fresh external review | 不能升為正式 release | 2-4h | P0 |
| 科學輸出與 HTML | recovery 完成前不可啟動 | downstream runtime | P1 |

## §6 Evidence Conflict Scan

沒有與科學假說的既有 NEGATIVE 衝突；此問題是 release engineering。最相關反例是同一
topic 先前 mode-tightening failed formal attempt 與本輪 V exact-key failure，兩者都支持
「不可用寬鬆 bypass，必須保留 fail-closed 並建立窄幅、可追溯 supersession」。

## §7 Decision

**Verdict: GO（由 PROBE 經 pilot + 雙 agent red-team 升級）**

Pilot 實際結果：

- 新路徑 V pilot 先因 direct `--verify` 非 FD-bound 正確拒絕；零寫入。
- 改用 inherited Python/V descriptors 後，第一次精確揭露 signed continuation-gate
  仍混用新/舊 verifier record；修正分層後再次執行 `exit=0`、`pass=true`。
- 50 個 critical descriptors 保留至 process lifetime；兩份 Ed25519 signatures、
  canonical/source bytes、三方原始 reviews 與 v7 authority 全部通過。
- `tumor_ref_source_receipt_promotion_verification.v3.json` 仍不存在；下游零新增。

Fresh red-team：Mill 與 Dalton 均獨立確認只修 V 不可行，R/C 同時硬綁舊
source、command 與 receipt schema；兩者一致建議 append-only recovery authority +
Vrec/Rrec/Crec，且新 receipt/schema/path 不得冒充舊程式輸出。

最強反方 failure modes：

1. 只修 V，R 仍以 signed authorization 要求舊 V/R/C exact records，fresh verifier 無法通過。
2. Supplemental authority 把自身 SHA 寫回自身或讓未驗證 loader 先執行，形成 circular trust。
3. Recovery receipt 沿用舊 schema 卻沒有揭露 supersession，造成 provenance 誤導。

Next action：產生新 recovery authority key；建立正式 Vrec/Rrec/Crec 與共用 authority
validator；fault-inject 後凍結 source，再取得三方 fresh review、發布並簽署 supplemental
authority，最後依 Vrec → Rrec → Crec 順序執行。

Decision lock: `Y`。Supersession 僅限 historical metadata verifier recovery；不得改動
canonical bytes、科學 payload、原 promotion signatures 或其他 downstream gates。

## §8 Pre-Sign Review Audit（2026-07-22）

第一次 frozen source 在 authority 建立前被兩個獨立 reviewer 正確判為 REJECT：

1. Vrec evidence 缺 recovery validator artifact，無法滿足 Rrec recursive relation。
2. Crec terminal verifier 混用 legacy signed 與 active recovery source/path。
3. Review probe command 只要求非空 list，未 exact 綁定。
4. Regression test 沒有直接覆蓋 terminal precommit 分支。

上述缺陷修補後，第二組 source 又由主 agent 自我審查撤回：probe command 雖已 strict，
但尚未把 `/usr/bin/env -i` 和 clean-room environment tokens 納入完整 argv。第二組沒有
任何 reviewer verdict 可被沿用。

第三組 source 因 recursive bool/int equality、duplicate-key JSON 與 retired-key pathname-only
TOCTOU 風險被正式 REJECT。第四組 frozen source 為 validator `1dfb3f3e...`、Vrec
`f1b29adf...`、Rrec `dee756a0...`、Crec `81175713...`。Clean-room 結果為 `10 passed`、
probe exit `0`、`pass=true`、`no_output_writes=true`、19 個 output slots absent。

第四組 fresh review 已完成：Mendel `019f8af7-757a-7c22-905c-b416760f7e34`、Nash
`019f8af7-a9bd-7fd0-b6e5-c88ec1fe41b7`、External Claude Opus 4.8
`783f2548-4de3-4c68-a3f5-eec0d0fb0626` 全部 `APPROVE`，High `0`、Medium `0`、
unresolved `0`。External 僅記錄 `1e400 -> float('inf')` 的 non-exploitable Low finding；
trusted numeric 的 exact-type gates 會拒絕該值。簽署 gate 因而升為 **GO-to-sign**，但仍須
以 frozen ceremony builder 原子發布、驗簽並確認 key retirement，才能啟動 Vrec/Rrec/Crec。

資源 preflight：48 CPU cores、472 GiB available RAM；`/big7_disk` available 679 GiB。
正式 chain 可使用既定並行資源，但禁止為追求併發另複製大型中間資料。

## §9 Post-Sign Rrec v1 Incident 與第二層決策（2026-07-22）

第四組 authority 已經三方 APPROVE 並完成簽署；authority SHA-256
`af7c2cf68b7edbfe3883cd65b273213985bb9058e1799984860c3a09e1725fda`，detached signature
SHA-256 `ab2026b3ae7fce9c2e602e06dd95760f1d386a5ec90525471600f50ba5b60243`。OpenSSL 與三個
runtime role 均獨立驗證通過，one-time private key 已退休為 mode `0000`。第一次 builder
使用 cnvtools Python 時因該 interpreter 無 `os.memfd_create` 而在簽名前 fail closed；五個
output slots 全部仍 absent、key 仍為 `0400`。同一 frozen builder 改由具 Linux memfd API
的 `/usr/bin/python3 -I -B` 執行 preflight 與正式簽署後成功，builder bytes 未改。

Vrec inherited-FD `--verify` 與 direct `--verify-and-record` 均通過，建立 13,422-byte、
mode `0444` receipt，SHA-256
`b1f9b22d570ca68f91fd3039c2d2f9c1956c1bd769f4b30c983ddf6a122a13a5`。

Rrec v1 隨後在建立 log/receipt 前 fail closed：generic
`recursively_validate_artifact_relations()` 依序將 signed
`focal_source_identity_transition.during_execution`（historical mode `0o664`、精確八鍵）與
`current`（live mode `0o444`、精確九鍵）都當成同一 live artifact identity 要求。兩者的
bytes、path、device、inode、size、mtime 一致，但 mode 與 ctime 本來就被 signed transition
明示為歷史差異，因此同一 live inode 不可能同時符合 `0664` 與 `0444`。Rrec v1 output log
與 receipt 均未建立；不得用 chmod 將 current 降回 `0664`，否則會破壞 current signed
identity 並掩蓋 provenance。

**第二層 Verdict: GO-to-append-only-v2 / NO-GO-to-bypass。**

關鍵假設：

1. Vrec v1 receipt 已完整且可信，但只作為 immutable predecessor evidence；fresh formal
   design review 認為由 Crec 同時消費 v1/v2 authority 會增加不必要的翻譯風險，因此 active
   Vrec/Rrec/Crec 全部使用同一 v2 authority，並新建 Vrec `.recovery.v2` receipt。
2. 唯一允許的 Rrec semantic 變動，是在 exact context
   `evidence.focal_source_identity_transition.during_execution` 將精確八鍵 record 驗成
   historical observation，而非 rebound 為 live identity；同一 transition 的 current
   九鍵 record仍必須 descriptor-bind、mode `0444`、link count `1`。
3. v1 authority、reviews、source、Vrec receipt與 key retirement全部 immutable；Vrec/Rrec/Crec
   v2 及 supplemental authority v2 使用新 source path、新 key、新 review與新 authority。
4. Rrec/Crec v1 output 從未存在且永久保持 absent；v2 只可原子建立新的
   `.recovery.v2` receipt/log/attestation/signature/incident slots，不得回填 v1 空槽。
5. 科學 scripts、7 datasets、469,849 sites、32 data contracts、claim ceiling與 downstream
   output path均不變。

Step → Verify：

1. 建立 transition-aware Rrec v2 → 驗證：historical exact-8/current exact-9 正例通過；
   context、path、mode、ctime、unchanged fields、extra/missing key 任一漂移皆 fail closed。
2. 建立 Vrec/Crec/validator/probe/authority builder v2 → 驗證：source AST/self-binding、prior v1
   authority signature、Rrec v1 failure evidence、Vrec receipt SHA 與所有 output absence contract
   均被 hash-bind。
3. 執行 clean-room probe → 驗證：pytest 全通過、source hashes/modes/nlink 精確、只允許既有
   v1 authority/Vrec receipt present，Rrec/Crec/downstream outputs仍 absent。
4. fresh Mendel/Nash/External Claude review → 驗證：三者 High `0`、Medium `0`、unresolved `0`
   才簽 authority v2。
5. 執行 Rrec v2 → 驗證：exit `0`、receipt/log mode `0444`、runner lines `1-358` only、
   downstream output仍 absent。
6. 執行 Crec v2 → 驗證：完整 downstream terminal signature、all 7 datasets/469,849 sites、
   HTML/QA/ledger全部通過；任何 post-sign failure建立 incident並維持 release blocked。

Fresh design review 另提出並採納一項治理修正：builder v2 不得自行生成內容為 APPROVE 的
review JSON。三位 reviewer 必須先各自完成唯讀審查；主 agent 將其 verdict/severity/summary
保存成 mode `0444` review evidence。Builder 只能讀取既有三份 review，驗證 exact
source/prior/scope digests、High/Medium/unresolved皆空後，才可發布 authority/signature。

## §10 v25 C Incident：tumor-REF v1→v6 Audit Lineage（2026-07-24）

**Cynefin：Complicated recovery，Verdict：GO；不得修改 immutable tumor-REF manifest。**

### Observation completeness

| Observation | 狀態 | Tier | 證據 |
|---|---|---|---|
| v25 authority、V、R 皆 PASS，C 在 final builder fail-closed | ✓ | L1 runtime | `results/m2v5_downstream_continuation_incident.recovery.v25.json` + completion log |
| tumor-REF manifest exact 綁定 v1 pre-audit | ✓ | L1 artifact | tumor-REF `run_manifest.json:inputs.primary_artifact_audit_pre` |
| formal builder canonical input 為 source-authorized v6 pre-audit | ✓ | L1 code/runtime | `build_all_ssnv_final_report_dataset.py:232,4592` + failed command |
| v1/v6 審查同一 site/assignment identities、`102842/308526`、artifact-set SHA | ✓ | L1 independent artifacts | v1/v6 audit JSON |
| v25 Python contract 宣告 v6 output slots，composed shell prefix 仍實際使用 v5 slots | ✓ | L1 source/runtime | v25 continuation constants + completion log |
| 修復後 split-window 能否對所有 chronology 反例 fail-closed | □ | L2 test required | v26 regression + read-only probe |

### Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | append-only re-audit lineage 可將舊 producer epoch 與新 downstream epoch分開 |
| 觀察支撐 | 20 | 失敗 stack、兩份 audit 與 manifest 三方一致 |
| 機制清晰度 | 20 | exact-path check 先失敗，之後 chronology 亦需 split-window |
| 反例風險 | 10 | 全新 recovery builder/finalizer/report trust chain 需重審 |
| 所需資源 | 10 | 不重跑 102,842 個 tumor-REF/cooccurrence tasks，但 formal review 與 C 需數小時 |
| **TOTAL** | **80/100** | **GO，仍受 fresh review hard gate 約束** |

**Falsifier observable**：若 lineage 假設錯誤，v1/v6 的 audited input identity、counts、
artifact-set SHA、scope、pass semantics 或時間順序任一不一致；或 tumor-REF 不在
`v1 finished ≤ producer start ≤ producer finish ≤ v6 start`，或其他 downstream 不在
`v6 finished ≤ producer interval ≤ post start`。任一情形必須無 final output 地 fail-closed。

### Assumption map and decision lock

- **High / known**：v1/v6 審查同一 immutable data plane；v25 terminal key 尚未使用但必須退役。
- **High / unknown**：recovery builder、result finalizer、report builder 的 transitive source/command
  relations；必須由 regression、internal agents 與 external Claude 分別驗證。
- **Decision lock = Y**：不可改寫 tumor-REF manifest、v1/v6 audit、v7 producer sources或已有
  cooccurrence artifacts；不可將 lineage 簡化為忽略 path；不可沿用 v25 terminal key。

### Step → Verify

1. 封存 v25 已簽 failure chain 與所有 C outputs → 驗證：active v25 authority/V/R/C/review/output
   slots 全空，failure evidence + key quarantine 皆 mode `0444/0700`。
2. 建立 recovery builder split-window lineage → 驗證：正例通過；path/hash/count/time/schema/
   artifact-set 變異皆在 output publication 前被拒絕。
3. 對齊 composed shell 與 v7 producer canonical v5 output slots → 驗證：rendered bytes、Python
   constants、strict/matched/CN canonical commands 及 actual log 四方一致。
4. 建立 fresh v26 authority → 驗證：全 source mode `0444`、tests/probe PASS、Mendel/Nash/
   external Claude 均 APPROVE 後才允許簽署。
5. 新 terminal key 執行 V→R→C → 驗證：dataset/report detached signatures、desktop/mobile
   HTML QA、post-audit、immutability audit 與 terminal success witness 全 PASS。
