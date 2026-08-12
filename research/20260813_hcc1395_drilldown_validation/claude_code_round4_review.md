<!--
建立時間: 2026-08-13
目標: 保存 Claude Code 對 legacy methyl browser、final source、machine receipts、報告與 artifact 的最終獨立唯讀驗證
處理範圍: legacy method/crosswalk、drilldown source/tests、direct-generated browser receipt、主報告/source_notes/artifact；不修改 immutable bundles
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/legacy_browser_method_audit.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/legacy_browser_visual_audit.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/legacy_browser_crosswalk.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/legacy_browser_visual_audit.generated.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_HCC1395_drilldown完整驗證與多樣本改進_01.md
-->

# Claude Code Round 4 — legacy + final source adversarial review

## 執行契約

- **模式**：唯讀；允許 `Read/Grep/Bash`，禁用 `Edit/Write/WebSearch/WebFetch`，`--no-session-persistence`。
- **Claude Code**：`2.1.226`。
- **Repo HEAD**：`73afaeac8e61c767241fa59c1ca6043a1c95290c`；working tree source 另有本輪未提交修改，因此 reviewer 明確不把它當 clean release。
- **輸入路徑**：主報告、source notes、implementation notes、legacy 方法/視覺附錄、crosswalk/visual receipts、`InterSubMod/scripts/drilldown/`、`InterSubMod/tests/test_drilldown_{contract,multisample_contract}.py`。
- **輸出路徑**：本檔；Claude Code 原始程序只輸出 stdout，不寫 repo。
- **退出狀態**：`exit 0`。

## 執行命令

```bash
claude -p \
  --permission-mode bypassPermissions \
  --allowedTools Read Grep Bash \
  --disallowedTools Edit Write WebSearch WebFetch \
  --no-session-persistence \
  --output-format text \
  --effort high \
  '<Round 4 adversarial prompt：核對 A legacy taxonomy/linkage/selection/crosswalk；B observationScope/code 0·8·16/asset policy；C sample/provenance/LCA；D direct-generated receipt；E report/source_notes/artifact，並 fresh 重跑 pytest 與五個 node --check>'
```

## Claude Code 原始裁決摘要

```text
OVERALL = ACCEPT（工程與資料產品層；科學層維持 BLOCKED）
SCIENTIFIC_STATUS = BLOCKED for scientific / external claims
Multi-sample = PARTIAL
Claim ceiling = descriptive observation + internal QA only

43 tests passed；五個 JavaScript node --check 全 OK。
immutable HCC1395_v1、HCC1395_v3 與 2026-07-26 legacy standalone
不得因本輪 ACCEPT 而回溯升級。
```

## 獨立驗證結果

| 挑戰面 | 判定 | Claude Code fresh evidence |
|---|---|---|
| Legacy A/B taxonomy | **PASS for audit finding** | 30,077/30,077 可重建且 0 mismatch；但 B 採 priority-first，4,025 中 3,544（88.0%）也 allele-high，因此不能當純 HP negative control |
| Raw / powered linkage | **PASS for audit finding** | raw `coread≥2` 20,535、powered 17,749、powered+both-somatic 17,141 為不同 gate；headline 不可稱 confirmed linkage |
| 472→14 selection | **PASS for audit finding** | 472 中 68 含 FP、27 無 powered+both-somatic branch；18 個 A≥3 eligible 只展示 14，4 個 omitted 都有 heatmap，故是 selection-incomplete subset、不是 prevalence |
| Legacy/current crosswalk | **PASS** | 30,077 / 19,849 / 16,566，coordinate Jaccard 49.658%；兩側缺 REF/ALT，只能 coordinate-only。L4-available legacy A 555 的 ALT/HP 為 412/218；legacy B 2,149 的 ALT/HP 為 1,540/1,489，反證一對一 mapping |
| Observation encoding | **PASS** | immutable v1 axis code range 0–8；`AXIS_UNTESTED=16` 是後續 source contract，不回填舊 bundle。claim ceiling 與 missing/untested 分離誠實 |
| Light asset policy | **PASS with boundary** | staging 明示 `igv=skip, panels=0`；完整影像 bundle 維持 NOT EVALUATED，沒有偷渡成已驗證 |
| Sample fail-closed | **PASS** | topology/MLHP/ISM/LCA/lineage/strict-edge wrong-sample fixtures 都 fail closed；substring alias 誤配亦有反例 |
| Receipt / provenance | **PASS for source contract** | receipt v2 保存 commit、script hash、argv、claim ceiling 與 scientific/truth flags；不代表 legacy receipt 被補齊 |
| LCA causal gate | **PASS** | per-chrom input/threads mismatch 進 C11；15/22 same BAM、0/22 same threads，`causal_ab_gate=FAIL`，4.969×只可 descriptive |
| Direct-generated browser | **PASS for runtime light QA** | desktop 1440/1440、mobile 390/390、errors `[]`、4 tables、106 denominator labels、0 fake cells、mobile overlap false、2 FAIL/0 SKIP→BLOCKED |
| Report/source_notes/artifact | **PASS** | top-line metrics 均可回到 CSV/JSON；cluster 99.99–100% 仍標 invalid；多樣本與 truth 缺口未 over-claim |

## Reviewer 指出的 blocker 與本輪 closure

1. **Round 4 引用懸空**：review 當下本檔尚未存在。現在以本檔保存命令、版本、fresh checks 與裁決，已關閉內部引用鏈；這不改變 scientific `BLOCKED`。
2. **Truth-set 缺口**：仍開放。沒有 SEQC2 truth/HighConf BED/TP-FP-FN/F1 receipt，任何準確率 claim 禁止。
3. **Legacy bundle 不可完整重建**：仍開放。v1 ISM root 缺失；v3 只有 unverifiable directory；候選版必須 clean rebuild 到新路徑。

## Reviewer minor notes 與修正

| Minor note | 處理 |
|---|---|
| generated visual receipt 不獨立保存 10 PASS，只能直接驗 2 FAIL/0 SKIP | 主報告改為只引用 receipt 可驗的 2 FAIL/0 SKIP；PASS count 不作判定 |
| Crosswalk Jaccard 的 569 分母與主報告 L4-available 555 並存 | JSON 新增 `legacy_A_no_L4_axis_code_0=14`、`legacy_A_L4_available=555` 與 denominator note |
| threads 差異的因果語意可更清楚 | 主報告與 artifact 補註：threads 不是生物 confounder，但破壞「唯一變項只有 LCA」的 controlled-comparison 前提 |

## 實際輸出片段

```text
...........................................                              [100%]
43 passed in 1.58s

五個 scripts/drilldown/assets/*.js：node --check exit 0
OVERALL = ACCEPT（工程與資料產品層；科學層維持 BLOCKED）
```

## Claim boundary

本輪 `ACCEPT` 只代表 working-tree source、contract tests、crosswalk/report consistency 與 direct-generated light runtime QA 通過。**它不會回溯改變 immutable v1/v3 或 legacy standalone 的 `BLOCKED / DO-NOT-CITE / observation-only` 狀態，也不完成多樣本 ISM/lineage 或 truth-set validation。**

## Post-review follow-up：ISM-bearing Figure 18

主審後新增 machine-selected methyl-positive detail QA；因此再以 Claude Code `2.1.226`、唯讀 `Read/Grep/Bash`、禁用 Edit/Write/Web 重驗 runner、generated receipt、Figure 18、主報告、source notes 與 artifact。

```text
FOLLOWUP=ACCEPT
NO_HARDCODE_IN_SCRIPT
selectedLocus=chr1:1320793 / axisCode 2
read=105, CpG=371, ALT/REF raw p=0.0010, effect field=0.029
no-ISM fallback raises RuntimeError
igv=skip, panels=0; scientific BLOCKED and full assets NOT EVALUATED
Remaining issue: none
```

Reviewer 確認 locus 是從 `bootData.l1.axis` 依非循環 global-axis bit 機器選出，script 內無 `1320793` 字面；raw p 仍標 exploratory、非 FDR/因果/cohort-wide 結論。此 follow-up 只擴張 UI 正例 coverage，不改變任何 immutable claim boundary。
