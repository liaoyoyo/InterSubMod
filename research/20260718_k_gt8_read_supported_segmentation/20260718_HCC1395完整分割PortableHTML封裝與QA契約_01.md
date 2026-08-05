<!--
建立時間: 2026-07-18 07:05 +08:00
目標: 凍結 HCC1395 k>8 read-supported segmentation 最終 portable HTML 的封裝命令、發布閘門與量化 QA
處理範圍: artifact.json 建置後的 canonical validation、Data Analytics portable packaging、browser/structural receipt、交付 checksum
關聯檔案:
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/summarize_hcc1395_full.py
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/build_report_artifact.py
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/deliver_portable_artifact_scrollbar_safe.mjs
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/baseline-runtime-audit.md
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/pre-decision-audit.md
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/report_delivery/20260718T014238Z_full_v7/report.html
-->

# HCC1395 完整分割 Portable HTML 封裝與 QA 契約

> **用 SCQA：先完成 canonical artifact validation，再由唯一的 shared reader 封裝一次；只有 browser verification、精確內容數、來源互動、SHA-256 與 release marker 全部通過才可交付 `report.html`（影響：高，信心：高）。**

## 0. 任務分類、狀態與不可跨越的邊界

- **Task Type**：B — Comprehensive validation。
- **服務目標**：G4（chr1–22 完整可重現）與 G5（外部可攜、可稽核交付）。
- **讀者與表面**：technical report；唯一 delivery mode 是 self-contained portable HTML。
- **目前狀態**：**COMPLETE**。v7 的 canonical artifact、portable HTML、正式 Chromium
  verification、source interaction、desktop/mobile post-QA 與 checksum 均已通過。
- **最終交付**：
  `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/report_delivery/20260718T014238Z_full_v7/report.html`。
- **正式 scope**：HCC1395 autosomal chr1–chr22、PASS biallelic sSNV；不含 chrX/chrY/MT、
  indel，也不包含 candidate-tree/VAF 排名。
- **禁止事項**：
  - 不讀寫或搬移正在執行的 full output。
  - 不修改 locked runner、extractor、partitioner。
  - 不手寫第二套 HTML、SVG、JavaScript 或 chart runtime。
  - 不把 fixture、subset 或 `snapshot.status != "ready"` 的 artifact 當最終報告。
  - 不覆寫既有 `artifact.json`、`report.html`、receipt 或 checksum；失敗 attempt 原樣保留供稽核。

本輪只建立 packaging/QA SoT。正式封裝必須等下列 upstream gate 全部完成：

1. full runner `_SUCCESS` 與 `receipt.json` 已由獨立 summarizer 驗證；
2. summarizer 產出的 `summary.json`、`summary.tsv`、`component_all.tsv.gz` 與 SHA sidecars 完整；
3. 若最終主張包含 hard-span sensitivity，50/100/200 kb span-grid receipt 與 summary 必須先 all-pass；
4. `build_report_artifact.py` 產出 `snapshot.status="ready"` 的正式 artifact；
5. Data Analytics `validate_artifact` 對**完整 payload**通過，才執行一次 `report:deliver`。

上述五項現已全部完成。最終 identity：

| 項目 | 結果 |
|---|---|
| final attempt | `20260718T014238Z_full_v7` |
| artifact | `d43fef934ff0d24deca94e12be6e1f65f6801013cdd7bd98a26aa525b6438500` |
| HTML | `f719b534ba5a02313b3d28451d2ae6f535bc7dceb607e1b39c79b0fd0c4f7efa` |
| delivery receipt | `b26f692afd872584579cac64ebeac6a43beef6af23234ca8ca65fd11620522b4` |
| formal verifier | `validation/package/verification=passed`；viewports `[1440,390]` |
| source interaction | `sourceDialog=passed`；`keyboard_menu_semantic_click` |

## 1. 已檢查的 canonical packaging 實作

### 1.1 固定 plugin 與執行環境

```text
PLUGIN_ROOT=/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599
package name=datascience-mcp-widgets
package version=0.2.8
Node=v22.22.1
npm=10.9.4
Chromium headless-shell=/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium_headless_shell-1217/chrome-headless-shell-linux64/chrome-headless-shell
```

已回讀的固定檔案 SHA-256：

| 檔案 | SHA-256 |
|---|---|
| `package.json` | `856f7cc8a4737576167795ca11431a5b64721caae2137f5bc614d4e7ea6855fa` |
| `deliver_portable_artifact.mjs` | `09d7afe25f76429ba1e99d08a5793af44c27738482d91e9b0bd2351ffd516066` |
| `build_portable_artifact.mjs` | `0b86883810cf23c81f563a48a7708283a6a4cadf11228ad261efb64426ef728e` |
| `verify_portable_artifact.mjs` | `b495c4cc34113fb2918118eac302f4e7e2152c1b1b9b63e3646cea97ecbf9b3f` |

正式交付前應重算以上四個 SHA；任一漂移就停止，不沿用本文件的 QA 判讀。

### 1.2 `npm run report:deliver` 實際契約

`package.json` 的唯一命令為：

```text
report:deliver = node skills/build-report/scripts/deliver_portable_artifact.mjs
```

必要參數：

```text
--input <artifact.json>
--output <report.html>
```

可選的 bounded verifier 參數：

```text
--ready-timeout-ms <positive milliseconds>   # default 5000
--action-timeout-ms <positive milliseconds>  # default 2500
--timeout-ms <positive milliseconds>         # default 10000
--screenshot <failure-only PNG>
```

單一命令依序執行：

1. 呼叫 canonical `validate_artifact`；
2. 建立 shared reader 與同源 semantic fallback；
3. 若 artifact 有 chart，使用既有 Chromium 擷取 shared-renderer SVG，再由相同 artifact 重建；
4. 比較 HTML 內嵌 payload 與 canonical artifact 的 exact equality；
5. 驗證 desktop `1440×1000` 與 narrow `390×844`；
6. 驗證 block/chart/table/metric 數量、nonzero geometry、horizontal overflow、source menu/dialog、外部 network request 與 browser error；
7. 全部成功才將 temporary candidate 原子 rename 成指定 output；
8. stdout 回傳 compact JSON delivery receipt，CLI 失敗為 exit `1`。

`report:deliver` 本身在 POSIX 可取代同名 output，因此本研究的外層命令必須先要求全新 attempt 目錄與 output 不存在。

### 1.2.1 實測 classic-scrollbar 邊界與最小修正

標準命令在長報告的 sandbox iframe 中揭露 reader runtime 的既有 CSS 邊界：
`.analytics-top-bar` 使用 `width:100vw`；classic vertical scrollbar 佔 15 px 時，
名義 1440 px iframe 的 `documentElement.clientWidth=1425`，top bar 仍為 1440 px，
造成 `scrollWidth=1433` 與左右各 7.5 px overflow。

本研究沒有修改 plugin cache，也沒有手改成品 HTML。新增
`deliver_portable_artifact_scrollbar_safe.mjs`，透過官方
`buildPortableArtifact(..., {runtimeHtml})` 與 `deliverPortableArtifact(..., {build})`
擴充點，在 runtime head 最後注入：

```css
.analytics-top-bar {
  width: 100% !important;
  max-width: 100% !important;
  margin-right: 0 !important;
  margin-left: 0 !important;
}
```

其餘 canonical validation、static-chart extraction、artifact exact-equality、雙 viewport
browser verifier 與 source dialog 全部仍由 plugin 原函式執行。wrapper 單元測試 `3/3 PASS`；
patched-runtime 正式 verifier 與最終 v7 delivery 均 PASS。

### 1.3 環境測試結果

```bash
cd /bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599
node --test \
  tests/portable-artifact-delivery.test.mjs \
  tests/report-delivery-contract.test.mjs
```

實際結果：`17/17 PASS`、`0 fail`。

較廣的 `portable-artifact-verifier.test.mjs` 開發測試在本環境缺少 `playwright-core`，因此不能當成已通過的 plugin-wide CI。這不等於正式 delivery 可豁免 browser QA：目前 headless-shell 已存在，所以**本報告 final gate 要求真實 delivery receipt 的 `stages.verification="passed"`**；不得以 `structural_only` 取代。禁止為本任務下載或安裝 browser/dependency。

## 2. 正式輸入與 artifact 建置命令

### 2.1 輸入角色

| 角色 | 正式輸入 | 必要狀態 |
|---|---|---|
| full summary JSON | `<SUMMARY_DIR>/summary.json` | `all_pass=true`、`comprehensive_all_pass=true` |
| per-chrom summary | `<SUMMARY_DIR>/summary.tsv` | exact chr1–22 + `ALL`，與 JSON cross-check |
| component detail | `<SUMMARY_DIR>/component_all.tsv.gz` | exactly 408 k>8 components |
| baseline runtime | `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/baseline-runtime-audit.md` | frozen audit |
| hard-span sensitivity | `<SPAN_GRID_ROOT>/receipt.json` | 若主張 50/100/200 kb，必須 all-pass |
| artifact builder | `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/build_report_artifact.py` | source SHA 在建置前重讀 |

`summary.json` 會保留 machine-local evidence path 供本地稽核；`artifact.json` 的 portable source metadata 只可保留安全檔名、run identity 與 SHA-256，**不可嵌入 `/big7_disk/...`、`file://`、`../`、credential 或 credential-bearing URL**。

### 2.2 Step → Verify

1. 驗證 summarizer inputs  
   → 驗證：formal mode 拒絕 subset；exact `22` chromosomes、`408` components、`47,570` target sites、`79,687` sSNV sites。
2. 建立 canonical `artifact.json`  
   → 驗證：builder exit `0`；stdout receipt `ok=true`、`snapshot_status=ready`；output 不可預先存在。
3. 呼叫 Data Analytics `validate_artifact`  
   → 驗證：完整 manifest/snapshot/source payload 通過，沒有 validation issues；此步只驗證，不 render。
4. 執行一次 portable delivery  
   → 驗證：npm exit `0`；delivery receipt 三 stage PASS，browser verification 必須是 `passed`。
5. 建立 checksum 與 release marker  
   → 驗證：artifact/HTML/receipt/plugin 四類 identity 均寫入 SHA ledger；marker 最後建立。

### 2.3 Artifact builder 命令

先指定**不存在**的新 attempt 目錄；不要重用失敗 attempt：

```bash
set -euo pipefail

TOPIC=/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_k_gt8_read_supported_segmentation
SUMMARY_DIR='<ABSOLUTE_VALIDATED_SUMMARY_DIR>'
SPAN_GRID_ROOT='<ABSOLUTE_VALIDATED_SPAN_GRID_ROOT>'
ATTEMPT_ID='<YYYYMMDDTHHMMSSZ_or_frozen_run_id>'
DELIVERY_DIR="$TOPIC/report_delivery/$ATTEMPT_ID"
ARTIFACT="$DELIVERY_DIR/artifact.json"
ARTIFACT_BUILD_RECEIPT="$DELIVERY_DIR/artifact_build_receipt.json"
ARTIFACT_BUILD_STDERR="$DELIVERY_DIR/artifact_build.stderr.log"
GENERATED_AT='<FROZEN_ISO8601_TIMESTAMP>'

test -d "$SUMMARY_DIR"
test -f "$SUMMARY_DIR/summary.json"
test -f "$SUMMARY_DIR/summary.tsv"
test -f "$SUMMARY_DIR/component_all.tsv.gz"
test -f "$SPAN_GRID_ROOT/receipt.json"
test ! -e "$DELIVERY_DIR"
mkdir -p "$DELIVERY_DIR"

python3 "$TOPIC/scripts/build_report_artifact.py" \
  --summary-json "$SUMMARY_DIR/summary.json" \
  --summary-tsv "$SUMMARY_DIR/summary.tsv" \
  --component-all "$SUMMARY_DIR/component_all.tsv.gz" \
  --baseline-runtime-audit "$TOPIC/baseline-runtime-audit.md" \
  --span-grid-summary "$SPAN_GRID_ROOT" \
  --exact-log-audit '<ABSOLUTE_VALIDATED_EXACT_LOG_AUDIT_ROOT>' \
  --hp-ps-unit-audit '<ABSOLUTE_VALIDATED_HP_PS_AUDIT_ROOT>' \
  --top-components 25 \
  --generated-at "$GENERATED_AT" \
  --output "$ARTIFACT" \
  >"$ARTIFACT_BUILD_RECEIPT" \
  2>"$ARTIFACT_BUILD_STDERR"

jq -e \
  '.ok == true
   and .surface == "report"
   and .snapshot_status == "ready"
   and (.dataset_rows.headline_metrics == 1)
   and (.dataset_rows.per_chromosome_metrics == 22)
   and (.dataset_rows.retention_by_chrom_method == 42)
   and (.dataset_rows.runtime_by_chrom_stage == 88)
   and (.dataset_rows.genomic_components == 408)
   and (.dataset_rows.top_components == 25)
   and (.dataset_rows.hp_ps_summary_metrics == 1)
   and (.dataset_rows.hp_ps_unit_distribution == 22)
   and (.dataset_rows.hp_ps_worst_units == 25)
   and (.dataset_rows.span_sensitivity == 4)' \
  "$ARTIFACT_BUILD_RECEIPT"

sha256sum \
  "$ARTIFACT" \
  "$ARTIFACT_BUILD_RECEIPT" \
  "$ARTIFACT_BUILD_STDERR" \
  >"$DELIVERY_DIR/artifact_build.sha256"
```

預期 dataset profile（含 span sensitivity）：

| Dataset | rows |
|---|---:|
| `headline_metrics` | 1 |
| `per_chromosome_metrics` | 22 |
| `retention_by_chrom_method` | 42（chr21 無 k>8 denominator，不造假成 0%） |
| `runtime_by_chrom_stage` | 88（22 chromosomes × extraction、partition total、nested pattern load、nested component loop） |
| `genomic_components` | 408 |
| `top_components` | 25 |
| `span_sensitivity` | 4 |
| `hp_ps_summary_metrics` | 1 |
| `hp_ps_unit_distribution` | 22 |
| `hp_ps_worst_units` | 25 |
| **總列數** | **638** |

最終 manifest profile 為 `30 blocks / 12 cards / 8 charts / 3 tables / 10 datasets /
638 rows`；browser verifier 將一個 12-card metric strip 展開計數，因此 rendered counts 為
`41 blocks / 12 metrics / 8 charts / 3 tables / 0 HTML`。

## 3. Artifact pre-package hard gates

在呼叫 `report:deliver` 前，canonical validator 與人工 source map 至少要確認：

### 3.1 結構與 payload

- [ ] top-level `surface="report"`。
- [ ] `manifest.title` 為讀者可讀標題。
- [ ] 第一個 block 是與 title 完全相同的 `#` heading。
- [ ] `manifest.blocks` 有明確 top-to-bottom reading path；每個 major section 各自一個 markdown block。
- [ ] `snapshot.status="ready"`；不得用 `accessIssues` 裝 optional caveat。
- [ ] `snapshot.datasets` 是 object of plain row arrays，不是 `{columns, rows}`。
- [ ] datasets `≤50`；每 dataset `≤2,000` rows；snapshot `≤3 MiB`；inline source characters `≤200,000`。
- [ ] 所有數值 JSON finite；percent 使用 fractional rate（`0.98` 表 98%），或明確用 number + `%` unit。
- [ ] 沒有 secret、token、credential、absolute local path、`file://` 或 parent traversal。

### 3.2 來源與可稽核性

- [ ] 每張 card/chart/table 有 inline `source` 或可解析的 `sourceId`。
- [ ] quantitative markdown block 只有在全部數字同源時才設 block-wide `sourceId`。
- [ ] 每個 source 至少保存：資料角色、portable 檔名/run id、SHA-256、filters、metric definitions、generated/executed timestamp。
- [ ] 本地完整路徑只留在本 QA 文件、summary evidence 與 SHA ledger；不得嵌入 portable artifact。
- [ ] `manifest.sources[]` 與 top-level `sources` 一致。

### 3.3 Narrative 與方法主張

- [ ] 開頭 technical/executive summary 直接回答：工程 assignment、read retention、tree-ready/ABSTAIN、runtime 與 sensitivity 結果。
- [ ] 先定義 component、primary-active、read retention、TREE_READY_LOCAL、ABSTAIN、weight-stable，再依賴這些數字。
- [ ] 明確區分：
  - 100% site assignment；
  - observed read constraint retention；
  - 可進 candidate-tree inference 的 tree-ready local block；
  - 不可宣稱的唯一真實演化樹。
- [ ] VAF 僅可描述為 segmentation 後候選樹相對排序的推測證據，不可回頭證明切點或真實演化。
- [ ] hard-span/no-cap runtime denominator清楚；cached partition wall 不得冒充 full BAM end-to-end runtime。
- [ ] limitations、robustness、next steps 與 further questions 都有可見 section role。

### 3.4 圖表 QA

正式含 span profile 的 artifact 預期：

```json
{"blocks":41,"charts":8,"html":0,"metrics":12,"tables":3}
```

- [x] 8 charts 全部使用 native artifact chart，不使用手寫 SVG/canvas/remote library。
- [ ] 每張 chart 的前或後有相鄰解釋段落：takeaway、讀法、限制/含意。
- [ ] gene-region scatter 使用一致觀察 grain：一列一 component；x=genomic start Mb、y=chromosome、size=pre-cap k、color=evidence/weight status。
- [ ] chr6/chr16 高密度極端區仍可見；不能因 top-N table 取代完整 408-point genomic chart。
- [ ] 所有 chart 有 neutral descriptive title、單位/分母明確 subtitle。
- [ ] `by <dimension>` 的文字必須由 axis、color、series、facet 或 direct label 實際編碼。
- [ ] single-series chart 不顯示無意義 legend；grouped chart 的 group names 可見。
- [ ] 不以顏色作唯一區分；light/dark/grayscale 下仍可辨。
- [x] 三張 table 都有合法 `defaultSort`，且欄位存在於 `columns`。

## 4. 唯一 portable delivery 命令

此命令只能在 `validate_artifact` 已通過後對全新 attempt 執行一次。stdout 只有 JSON
receipt；所有 stderr 與 failure screenshot 都留在 attempt 目錄。由於 §1.2.1 的已重現
plugin runtime bug，最終命令使用 repo wrapper，但封裝與驗證仍呼叫 plugin 官方函式。

```bash
set -euo pipefail

PLUGIN_ROOT=/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599
TOPIC=/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_k_gt8_read_supported_segmentation
DELIVERY_DIR='<ABSOLUTE_FRESH_DELIVERY_DIR_FROM_SECTION_2>'
ARTIFACT="$DELIVERY_DIR/artifact.json"
HTML="$DELIVERY_DIR/report.html"
DELIVERY_PENDING="$DELIVERY_DIR/report.delivery_receipt.pending.json"
DELIVERY_RECEIPT="$DELIVERY_DIR/report.delivery_receipt.json"
DELIVERY_STDERR="$DELIVERY_DIR/report.delivery.stderr.log"
FAILURE_SCREENSHOT="$DELIVERY_DIR/report.verification-failure.png"
SHA_LEDGER="$DELIVERY_DIR/portable_delivery.sha256"
SUCCESS_MARKER="$DELIVERY_DIR/_PORTABLE_SUCCESS"
EXPECTED_COUNTS='{"blocks":41,"charts":8,"html":0,"metrics":12,"tables":3}'

test -f "$ARTIFACT"
test ! -e "$HTML"
test ! -e "$DELIVERY_PENDING"
test ! -e "$DELIVERY_RECEIPT"
test ! -e "$DELIVERY_STDERR"
test ! -e "$SHA_LEDGER"
test ! -e "$SUCCESS_MARKER"

node "$TOPIC/scripts/deliver_portable_artifact_scrollbar_safe.mjs" \
  --plugin-root "$PLUGIN_ROOT" \
  --input "$ARTIFACT" \
  --output "$HTML" \
  --ready-timeout-ms 15000 \
  --action-timeout-ms 10000 \
  --timeout-ms 60000 \
  --screenshot "$FAILURE_SCREENSHOT" \
  >"$DELIVERY_PENDING" \
  2>"$DELIVERY_STDERR"

jq -e \
  --arg expected_html "$(realpath -m "$HTML")" \
  --argjson expected_counts "$EXPECTED_COUNTS" \
  '.ok == true
   and .html == $expected_html
   and .stages.validation == "passed"
   and .stages.package == "passed"
   and .stages.verification == "passed"
   and (.browserWarning | not)
   and .counts == $expected_counts
   and .sourceDialog == "passed"
   and (.viewports == [1440, 390])' \
  "$DELIVERY_PENDING"

test -s "$HTML"
test ! -e "$FAILURE_SCREENSHOT"
mv -- "$DELIVERY_PENDING" "$DELIVERY_RECEIPT"

sha256sum \
  "$ARTIFACT" \
  "$HTML" \
  "$DELIVERY_RECEIPT" \
  "$DELIVERY_STDERR" \
  "$PLUGIN_ROOT/package.json" \
  "$PLUGIN_ROOT/skills/build-report/scripts/deliver_portable_artifact.mjs" \
  "$PLUGIN_ROOT/skills/build-report/scripts/build_portable_artifact.mjs" \
  "$PLUGIN_ROOT/skills/build-report/scripts/verify_portable_artifact.mjs" \
  >"$SHA_LEDGER"

jq -n \
  --arg created_at_utc "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  --arg artifact_sha256 "$(sha256sum "$ARTIFACT" | awk '{print $1}')" \
  --arg html_sha256 "$(sha256sum "$HTML" | awk '{print $1}')" \
  --arg delivery_receipt_sha256 "$(sha256sum "$DELIVERY_RECEIPT" | awk '{print $1}')" \
  --arg sha_ledger_sha256 "$(sha256sum "$SHA_LEDGER" | awk '{print $1}')" \
  '{
     schema_name: "intersubmod.portable_report_success",
     schema_version: "1.0.0",
     all_pass: true,
     created_at_utc: $created_at_utc,
     artifact_sha256: $artifact_sha256,
     html_sha256: $html_sha256,
     delivery_receipt_sha256: $delivery_receipt_sha256,
     sha_ledger_sha256: $sha_ledger_sha256
   }' >"$SUCCESS_MARKER"

jq -e '.all_pass == true' "$SUCCESS_MARKER"
```

### 失敗語意

- npm、validator、browser 或 `jq` 任一步非零：**沒有 `_PORTABLE_SUCCESS` 就不可交付**。
- 失敗 attempt 的 `artifact.json`、pending receipt、stderr、HTML 或 failure screenshot 全部保留；不得刪除或覆寫。
- 修正 artifact 後使用新的 `ATTEMPT_ID` 與新目錄重新開始。
- `stages.verification="structural_only"`：只能表示 payload equality + runtime/fallback roots 通過；本環境已有 browser，故 final 判定為 **HOLD**。
- `sourceDialog!="passed"`、counts 不一致、出現 `browserWarning` 或 failure screenshot：final 判定為 **HOLD**。

## 5. 最終可量化 QA 與預期輸出

### 5.1 成功 receipt 必須同時符合

| Gate | 期望值 |
|---|---|
| CLI exit | `0` |
| `ok` | `true` |
| `stages.validation` | `passed` |
| `stages.package` | `passed` |
| `stages.verification` | `passed` |
| browser warning | absent |
| blocks / metrics / charts / tables / html | `41 / 12 / 8 / 3 / 0` |
| source dialog | `passed` |
| viewports | `[1440, 390]` |
| failure screenshot | absent |
| external requests | `0`（由 verifier hard gate） |
| browser errors | `0`（由 verifier hard gate） |
| `_PORTABLE_SUCCESS` | present、`all_pass=true` |

### 5.2 最終 attempt 目錄

```text
report_delivery/<ATTEMPT_ID>/
├── artifact.json
├── artifact_build_receipt.json
├── artifact_build.time_v.txt
├── canonical_validation_receipt.json
├── report.html
├── report.delivery_receipt.json
├── report.delivery.stderr.log
├── postqa_1440.json / postqa_1440.png
├── postqa_390.json / postqa_390.png
├── portable_delivery.sha256
└── _PORTABLE_SUCCESS
```

`report.verification-failure.png` 只在失敗時存在；成功 package 不應有此檔。

### 5.3 交付前最後 readback

```bash
set -euo pipefail
DELIVERY_DIR='<ABSOLUTE_SUCCESSFUL_DELIVERY_DIR>'

jq -e '.all_pass == true' "$DELIVERY_DIR/_PORTABLE_SUCCESS"
jq -e \
  '.ok == true
   and .stages.verification == "passed"
   and .counts == {"blocks":41,"charts":8,"html":0,"metrics":12,"tables":3}
   and .sourceDialog == "passed"
   and .viewports == [1440,390]' \
  "$DELIVERY_DIR/report.delivery_receipt.json"

(cd "$DELIVERY_DIR" && sha256sum -c portable_delivery.sha256)
test -s "$DELIVERY_DIR/report.html"
test ! -e "$DELIVERY_DIR/report.verification-failure.png"
```

注意：目前 `portable_delivery.sha256` 同時含 attempt 內檔案與 plugin 絕對路徑；`sha256sum -c` 必須在原 workstation 與相同 plugin cache 上執行。對外分享的真正 portable artifact 只有單一 `report.html`；receipt、artifact 與 SHA ledger 是本地審核 companion，不是 HTML runtime sidecar。

## 6. 完成定義

只有以下全部成立，才可對使用者說「Portable HTML 已完成」：

1. 正式 full summary 與 optional span-grid 都符合本次主張範圍；
2. `artifact.json` 是 `ready`，不是 fixture/partial；
3. canonical `validate_artifact` 與 `report:deliver` 內部 revalidation 都通過；
4. 真實 Chromium browser verification 為 `passed`；
5. exact payload、內容 counts、source interaction、desktop/narrow layout 與 no-network gate 通過；
6. HTML、artifact、receipt、plugin hashes 已記錄；
7. `_PORTABLE_SUCCESS` 最後建立；
8. 對外 handoff 連結指向 `report.html`，並保留「局部 read-supported segmentation，不是唯一真實演化樹」的主張上限。
