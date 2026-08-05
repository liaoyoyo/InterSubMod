<!--
建立時間: 2026-07-22
目標: 封存 v25 canonical 與未簽 v21 外審，保留 consumer rotation 缺口的完整證據
處理範圍: 7 datasets / 469,849 LongPhase-S recalibrated FILTER=PASS sSNV Task-B release chain
關聯檔案: canonical/, reviews/, source_snapshot/
-->

# v25 因 v7 consumer rotation 未完成而封存

## 判定

v25 canonical pytest 本身為 743 passed、0 failures/errors/skips，test-run manifest 亦已完成 Ed25519 簽章；但兩個正式 downstream runner 對 v7 authority 會 deterministic fail-closed，因此本批次不得用來啟動全量正式分析，也未簽署外部 review attestation、未鑄造 v7 source authority。

## 阻擋證據

- `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_cooccurrence_v6_source_locked.sh` 已指向 v7 authority/v10 public key，內部 `jq` 卻仍要求 v6 authority ID 與 v9 key。
- `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_m2v5_recovered_completion_chain.sh` 有相同矛盾。
- 外部 Claude Reviewer B 將此列為 `MEDIUM NB1`；主 agent 獨立重現同一靜態矛盾。

## 封存內容

- `canonical/`: v25 JUnit、stdout/stderr、raw evidence、signed test manifest、signature、bound test FD symlinks。
- `reviews/`: Reviewer A/B v21 request、command、stream transcript、stderr 與未簽 process attestation。
- `source_snapshot/`: v25 manifest 綁定的 protected、support、test source bytes，保留原 mode/timestamp。

## 後續門檻

修正兩個 consumer 並補正向 regression 後，必須使用新 test signing key 產生新的 canonical manifest，再重跑 A/B 外審；不得重用本批次的 source-set 或 approval claim。
