<!--
建立時間: 2026-08-01 17:40 +0800
目標: Exact-PS 全 7 技術資料集最終 HTML observation report 的單一重生與驗證入口
處理範圍: authority manifest、denominator registry、Python builder、standalone HTML、responsive/no-JS/print QA、immutable release
狀態: validated
關聯檔案:
  - InterSubMod/scripts/analysis/build_exact_ps_observation_report.py
  - InterSubMod/research/20260801_exact_ps_observation_report/scripts/finalize_exact_ps_observation_report.py
  - InterSubMod/research/20260801_exact_ps_observation_report/scripts/qa_exact_ps_observation_report.py
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
-->

# Exact-PS 全資料最終 HTML Observation Report

> **結果：VALIDATED DERIVED OBSERVATION。** 已由 13 個 hash-verified authority
> artifacts 與 19 列 denominator registry 產生全 7 技術資料集 standalone HTML；
> 1440／1024／390／320 px、no-JS 與 A4 列印均通過。CN/LOH 固定標為
> `NOT_INTEGRATED`，甲基固定為 `association-only`。

## 1. 服務目標

- **G3**：把 read-level genetic／methyl evidence 以不越過辨識界線的方式呈現。
- **G4**：以同一契約整理 7 technical datasets／6 biological IDs。
- **G5**：提供可重生、可稽核、可交接的 report data、HTML、QA 與 receipt。

這是衍生觀察層，不修改 canonical JSON／receipts，也不建立第二份數值 SoT。

## 2. 輸入

```text
InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/
├── authority_manifest.json
└── denominator_registry.tsv
```

Builder 會在輸出前重新驗證：

1. authority manifest schema、13 個 exact artifact IDs 與 SHA-256；
2. strict READY 的 9 個 nested bundle identities；
3. denominator registry 的 19 個 metric、分母、百分比與 authority reference；
4. cohort、per-sample、topology、active-k、methyl 與 candidate-count 守恆；
5. CN/LOH、methyl 與 tree-decision 的 claim boundary。

Machine-readable 契約：

- [report_data.schema.json](schemas/report_data.schema.json)（JSON Schema draft-07）
- 22 個頂層必填欄位、7 個 exact sample IDs、23 個 named conservation checks。

## 3. 正式執行命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod

python3 research/20260801_exact_ps_observation_report/scripts/finalize_exact_ps_observation_report.py \
  --authority-manifest /big7_disk/liaoyoyo2001/InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json \
  --denominator-registry /big7_disk/liaoyoyo2001/InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv \
  --release-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report \
  --release-id all7_v1
```

這個 finalizer 是完整資料鏈最後的獨立 post-run stage：

```text
authority readiness
→ report_data.json
→ offline standalone HTML
→ responsive / no-JS / print QA
→ release_receipt.json（最後有效標記）
```

不直接修改 2026-07-24 frozen runner，因 topology census 與 methyl sidecar 是其後續
artifact；只有 runner 結束並不代表所有 HTML 依賴都已 ready。

## 4. 正式輸出

輸出根：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
└── 20260801_exact_ps_observation_report/all7_v1/
```

主要檔案：

- [Standalone HTML](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report/all7_v1/report.standalone.html)
- [Report data](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report/all7_v1/report_data.json)
- [Build receipt](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report/all7_v1/report_build_receipt.json)
- [Browser QA receipt](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report/all7_v1/qa/browser_qa_receipt.json)
- [A4 PDF](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report/all7_v1/qa/report_A4.pdf)
- [Release receipt](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report/all7_v1/release_receipt.json)

正式輸出片段：

```json
{
  "all_pass": true,
  "statuses": {
    "data_technical_status": "FROZEN_TECHNICAL_RESULTS_WITH_EXPLICIT_ABSTAIN",
    "scientific_completeness_status": "INCOMPLETE_WITH_ABSTAIN",
    "report_build_status": "PASS",
    "browser_qa_status": "PASS",
    "release_status": "VALIDATED_DERIVED_OBSERVATION"
  }
}
```

## 5. 報告內的主要數字

```text
98,955 final groups
├─ 13,014 no active ALT
└─ 85,941 mutation-bearing
   ├─ 75,224 complete minimum families
   │  ├─ 71,955 read-AF ranked
   │  │  ├─ 39,648 unique best tree
   │  │  ├─ 23,858 tied, same exact rooted-unlabeled topology
   │  │  └─  8,449 tied, cross exact topology
   │  ├─  3,224 zero denominator
   │  └─     45 AF recurrence-screen abstain
   └─ 10,717 resource abstain
```

另外：

- one exact rooted-unlabeled topology：63,506／71,955＝88.2579%。
- one coarse class：63,511；比 exact topology 多 5 個 cross-exact／same-coarse units。
- coarse geometry：Direct-only 36,267、Single-only 22,135、
  Sister+direct 2,890、Sister-only 2,219。
- methyl：1,045 formal、811 evaluable、3 robust、627 no robust、181 confounded、
  234 not evaluable。

## 6. 驗證

### Python 契約測試

```bash
python3 -m unittest -v \
  research/20260801_exact_ps_observation_report/tests/test_build_exact_ps_observation_report.py \
  research/20260801_exact_ps_observation_report/tests/test_finalize_exact_ps_observation_report.py \
  research/20260801_exact_ps_observation_report/tests/test_qa_exact_ps_observation_report.py
```

實際結果：

```text
Ran 12 tests in 9.735s
OK
```

涵蓋 authority hash 漂移、duplicate artifact、denominator row loss、non-empty output、
QA failure、偽造 build／QA receipt、HTML／JSON 數值不一致、schema drift 與既有
release overwrite。

### Browser／print QA

```text
desktop  1440×1000  PASS
laptop   1024×900   PASS
mobile    390×844   PASS
narrow    320×700   PASS
no-JS    1024×900   PASS
A4 PDF     12 pages PASS
console errors       0
page errors          0
external requests    0
```

每個 viewport 都會驗：

- HTML embedded JSON 與 `report_data.json` 完全相同；
- 7 個 exact dataset IDs 與每列 bound values；
- 4 個 headline metrics；
- 3 個 SVG chart 的 title／description；
- page-level horizontal overflow 為 0。

### Report-data schema

```bash
/usr/bin/jsonschema \
  -i /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report/all7_v1/report_data.json \
  research/20260801_exact_ps_observation_report/schemas/report_data.schema.json
```

結果：退出碼 `0`。重複／缺失 sample、把 CN/LOH 改成已整合、開啟 methyl
topology rescoring、偽稱 tree decision 已逐列物化或缺少 conservation check
都會被 schema 拒絕。

## 7. 科學邊界

- read-AF 是 candidate family 內、CN/LOH 未校正的 deterministic preference。
- 39,648 個 unique AF winner 只能稱 `RAW_AF_UNIQUE_REPRESENTATIVE`。
- 23,858 個同拓撲並列只能稱 `TOPOLOGY_REPRESENTATIVE_ONLY`。
- 8,449 個跨拓撲並列只能輸出 candidate set／consensus backbone。
- 10,717 個 family incomplete 必須 abstain。
- 不同 PS×HP local blocks 不可拼成一棵全樣本 biological phylogeny。
- graph geometry 不等於已確認 clone／subclone。
