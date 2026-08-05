<!--
建立時間: 2026-07-22
目標: 紀錄 tumor-REF receipt promotion v3 的三方正式唯讀審查與 review artifact publication
處理範圍: frozen P/V/R/C、external read-only probe、Mendel/Nash/External Claude Opus verdict；不含科學結果確認
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/external_formal_readonly_probe_v1.py
-->

# Tumor-REF receipt promotion v3 三方正式審查紀錄

## 結論

三位 reviewer 對同一組 frozen source 均回傳無條件 `APPROVE`，沒有可重現的
HIGH/MEDIUM finding，`conditions=[]`。本 session 另重跑唯讀 probe，結果為
`pass=true`、`filesystem_artifacts_written=0`。因此三份 source-bound v3 review
artifacts 已可發布；這只解除 promotion 前的程式完整性 gate，不等於確認科學結論。

## 固定審查身分

| 角色 | size | SHA-256 | mode |
|---|---:|---|---:|
| P `promote_tumor_ref_source_receipt_v2.py` | 157,669 | `02fb9039b362fa261619b2236ddb764db278b23ae4467472fe2caa106770e06c` | `0444` |
| V `verify_tumor_ref_receipt_promotion_v2.py` | 120,163 | `03ff3b32368efffafa35491e04621508a46134d36407f060c3da12f90f2432a8` | `0444` |
| R `replay_m2v5_runner_only_gates_v1.py` | 118,987 | `10f1607aca3ef1a99da7fd77dcd6a207e0ba7003a6e3547a35b28926771fd694` | `0444` |
| C `continue_m2v5_after_tumor_ref_promotion_v1.py` | 277,381 | `f7b77bd16bd86ae1cbd97e85eebb38a882998b09bc9228fa5b045abfc0ffcfbd` | `0444` |
| Probe `external_formal_readonly_probe_v1.py` | 24,050 | `5d0c457d73547ab47810aabce647abc421ed73c2974a6257381094c45e907ed8` | `0444` |

## Reviewer 結果

| Reviewer | 執行身分 | Verdict | Finding / condition |
|---|---|---|---|
| Mendel | agent `019f8a9f-6da9-7d41-a45d-649466dc78fd` | `APPROVE` | 無；probe `pass=true`、零寫入 |
| Nash | agent `019f8a9f-730d-7a83-9f4c-b9e944e2a45a` | `APPROVE` | 無；`findings_closed=true`、`conditions=[]` |
| External Claude Opus | Claude Code session `9ea822ca-6737-4b44-bdf5-f87b6a765396`，`claude-opus-4-8` | `APPROVE` | 無；逐項覆核七個 failure surfaces |

外部 Claude 第一次 session `2f67f002-9468-4ecd-84c7-b27ee7fb405f` 使用
`--permission-mode plan`，只完成計畫與靜態定位，沒有執行 probe，因此明確判定為
**不構成正式 verdict**。第二次改用 `dontAsk`，停用 Edit/Write，Bash 自動允許範圍
限縮為 `stat`、`sha256sum` 與 sanctioned probe。主機缺 `socat`，Claude 內建 OS
sandbox 未啟用；此限制由唯讀 prompt、工具白名單、mode-0444 source 及審查後本
session 獨立 probe 補償，並保留在風險紀錄中。

## 本 Session 獨立 Probe

輸入：

`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/external_formal_readonly_probe_v1.py`

執行命令：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B \
  research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/external_formal_readonly_probe_v1.py
```

實際輸出摘要：

```text
pass=true
filesystem_artifacts_written=0
review_slots_absent=true
promotion_slots_absent=true
downstream_slots_absent=true
producer faults: exit0 accepted; exit1/SIGTERM/SIGKILL rejected
contract matrix: auth=17, completion=29, prepare=9, promote=12, runtime roles=11
```

## Review Artifact Publication

Builder：

`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/build_tumor_ref_promotion_reviews_v3.py`

- size: `10,296`
- SHA-256: `6c6ad303192f47d4aaf4fd34fe9c8805b9053071d3a02559e8b3464311f36f6b`
- 執行時 mode: `0664`

執行命令：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B \
  research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/build_tumor_ref_promotion_reviews_v3.py
```

實際輸出：`review_count=3`、`reviewed_source_count=15`、
`trusted_key_count=3`、`pass=true`。

| Artifact | size | SHA-256 | mode | inode |
|---|---:|---|---:|---:|
| `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/reviews/20260722_tumor_ref_receipt_promotion_mendel.v3.json` | 6,314 | `48012c89e3f1f4141fa7ed27d7ab6163b8ac2a06c798933e9b98e8a29b34269d` | `0444` | 615827522 |
| `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/reviews/20260722_tumor_ref_receipt_promotion_nash.v3.json` | 6,310 | `ac6b706da427aecef8b442984f930d5e4d3f4200b5dac7a230a9d3cbf2ed20e8` | `0444` | 615827523 |
| `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/reviews/20260722_tumor_ref_receipt_promotion_external_claude_opus.v3.json` | 6,342 | `4dd6dddb76ee1d3816409a334d5b35406fbd27065c606de25dbee7d7203a93c9` | `0444` | 615827524 |

三份 artifact 均為 link count 1、互異 inode，且各自精確綁定相同的 15 份
reviewed source 與 3 把 trusted public key。下一個 gate 是 P 的
`--prepare-authorization`；P 會再獨立重開並驗證這三份 artifact，而不是信任本文件。
