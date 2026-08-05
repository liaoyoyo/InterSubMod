<!--
建立時間: 2026-07-18
目標: 留存 positional-singleton 報告產生器的第二輪外部 Claude Code 唯讀複核
處理範圍: report builder、targeted tests、claim 分層、來源與簽章驗證
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/build_positional_singleton_report.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_build_positional_singleton_report.py
-->

# Positional singleton 報告產生器外部 Claude Code 複核

## 審查資訊

- Reviewer: 外部 Claude Code CLI
- 模式: `--permission-mode dontAsk --tools Read,Grep,Glob`
- 範圍: 唯讀；未修改檔案、未執行 pytest
- 審查時 targeted tests 狀態: `17 passed`
- 原始 verdict: **APPROVE_WITH_RESIDUAL_RISKS**

執行命令：

```bash
/bip7_disk/liaoyoyo2001/.local/bin/claude -p \
  --permission-mode dontAsk --tools 'Read,Grep,Glob' -- '<review prompt>'
```

## 已確認關閉的前次問題

1. Ed25519 receipt tamper、mode drift、public-key SHA drift 已有有效
   negative tests。
2. M2 分母已清楚區分：
   `30 PASS + 18 FAIL = 48 evaluable`、`5,913 NOT_EVALUABLE`、
   `44,471 NOT_RUN`。
3. `30/48` 明示為高度選擇的 conditional subset，不外推至全體。
4. Tumor-REF 的 `1,408/3,981` REF replication、`2,573/3,981`
   non-replication 與 M2 `7/7` non-replication 均未被寫成 subclone
   confirmation。
5. 原子 no-replace publication、identity checks、offline HTML 與主要
   accessibility 結構可接受。

## 第二輪 findings

| 等級 | Finding | 審查時狀態 |
|---|---|---|
| High | HTML claim ladder 把 `L1/L2` 同時用作已達成 evidence tier 與未達成 cellular tier，語意衝突 | 需修 |
| Medium | M1/M2 funnel 守恆依賴固定常數，缺獨立算術斷言 | 需修 |
| Medium | `status_census`、group distribution、breakdowns 等圖表來源未與正文 counts 機械對齊 | 需修 |
| Medium | Supplemental audit summary 本身不在 signed final release identity，衍生欄位需更強來源錨定 | residual |
| Medium | tumor-REF、candidate duplicate/count、XML errors 等部分 guard 缺 negative tests | 需補 |
| Low | Tumor-REF conservation 含恆真式、amber badge 對比不足、部分顯示值硬編碼 | 需清 |

## 本地修正紀錄

外部審查結束後完成：

1. 將未完成層改名為「細胞級正交確認」，不再重用 `L1/L2`。
2. 新增 M1 evaluability、M2-within-M1、M2-full-singleton 三組守恆。
3. 新增 M1/M2 status census、methyl-group distribution、
   dataset/truth/dataset-by-truth、50 kb positional recomputation、
   all-site branch 與 biological-sample identity 驗證。
4. 依 summary 記錄重新驗證 v4 source-authority Ed25519 signature、
   authority ID、source-set SHA、authorized HEAD、兩份 external
   `APPROVE` 與 supplemental auditor identity。
5. 補 candidate count/duplicate、XML errors/malformed、tumor-REF
   subset/missing/boolean/stability、technical overlap drift、HEAD drift、
   mode drift 等 negative tests。
6. 移除 tumor-REF 恆真 conservation guard；amber badge 改為較深色；
   動態化候選數與 illustrated dataset contract。

修正後 targeted tests：

```text
40 passed in 2.50s
```

真實正式 summary 新增 gate：

```text
summary_contract=PASS
source_authority_chain=PASS
checks=39
singleton_sites=50432
m1_flagged=5961
m2_pass=30
```

## 尚待終版後驗證

- 完成正式 signed final dataset/report 後，重跑全 canonical pytest。
- 產出 standalone HTML 後執行 desktop/mobile browser QA 與 screenshot
  檢視。
- 對最終 artifact 再做一次外部 Claude Code read-only review。
- Supplemental report receipt 需明示連結兩份 signed formal receipts；若另建
  detached signature，必須標示它是 supplemental integrity seal，不是原 v4
  release authority 的追溯核准。
