<!--
建立時間: 2026-07-23
目標: 盤點 Data Analytics portable HTML artifact builder 的 schema、native chart、UPGMA/heatmap 呈現方式與 browser QA 合約
處理範圍: 只讀檢查 plugin 0.2.8-13ceeea1f599 與 InterSubMod 既有 delivery receipts；未修改最終 artifact
關聯檔案: InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/artifact.json
-->

# Portable HTML artifact contract 盤點

## TL;DR

**可用，但 UPGMA dendrogram 與帶 NA 的 read×CpG methylation heatmap 不能只靠 native heatmap。** `bar`、`funnel`、`heatmap` 都是原生 chart type；`image`、`dendrogram`不是。最穩妥的報告組成是：

1. 以 native `bar` 呈現 734 個 M1 位點的群數分布。
2. 以 native `funnel` 呈現 8,279 → 8,074 → 734 的 gate（若有完整檢核分母）。
3. 以 native `heatmap` 呈現無 NA 語意問題的 read×read distance matrix。
4. 以 reproducible PNG/SVG 將 UPGMA dendrogram + read×CpG methylation pattern 嵌入 `type: "html"` block。
5. 報告仍必須至少有一個 native chart block；只放 UPGMA 圖的 HTML block 會不合 report contract。

## 1. 可用性結論

| 需求 | Portable artifact 能力 | 本報告建議 |
|---|---|---|
| M1 群數分布 | native `bar` / `horizontalBar` | 使用 native `bar` |
| 8,279 → evaluable → 734 | native `funnel` | 可使用，但每個 stage 必須是真實包含關係 |
| read×read distance matrix | native `heatmap` | 可使用；用 tidy-long rows |
| read×CpG methylation（有 NA） | native heatmap 無法呈現 NA 為獨立灰色 | 用 reproducible static PNG/SVG |
| UPGMA dendrogram | 沒有 native dendrogram/tree chart | 與 methylation heatmap 一起產生靜態複合圖 |
| 外部 PNG/SVG 檔 | 沒有 native `image` block | 轉 data URI 嵌入 sandboxed `html` block |
| 離線交付 | self-contained HTML | 不可依賴 CDN、remote image、sidecar data |

直接證據：既有 singleton derivative v4 artifact 已同時放入 4 個 native charts（`horizontalStackedBar`、`funnel`、2×`heatmap`）與 2 個 data-URI PNG HTML blocks；經相容性 QA 後成功交付，輸出 HTML 為 1,372,176 bytes。

- Artifact: `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/artifact.json`
- Success receipt: `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/portable_builder_compat_receipt.json`

## 2. Canonical artifact 最小合法範例

下列範例已在本機以 JavaScript MCP server 的 `validate_artifact` 實際驗證，結果為：

```text
ok=true
surface=report
manifest_title=HCC1395 M1 穩定甲基多群
dataset_count=1
snapshot_status=ready
```

```json
{
  "surface": "report",
  "manifest": {
    "version": 1,
    "surface": "report",
    "title": "HCC1395 M1 穩定甲基多群",
    "description": "734 個 M1 穩定多群位點的觀察報告",
    "generatedAt": "2026-07-23T00:00:00Z",
    "charts": [
      {
        "id": "m1_group_counts",
        "title": "M1 位點依群數分布",
        "intent": "comparison",
        "question": "734 個 M1 位點主要分成幾群？",
        "rationale": "分類長條圖直接比較各群數對應的位點數。",
        "type": "bar",
        "dataset": "m1_group_counts",
        "source": {
          "id": "m1_counts_sql",
          "label": "M1 group-count staging query",
          "path": "queries/m1_group_counts.sql",
          "query": {
            "engine": "SQLite",
            "language": "SQL",
            "sql": "SELECT group_count, site_count FROM m1_group_counts ORDER BY group_count",
            "description": "Loads audited HCC1395 M1 site counts by methylation group count.",
            "tables_used": ["m1_group_counts"],
            "filters": ["dataset=HCC1395", "M1 flagged=true"]
          }
        },
        "encodings": {
          "x": {
            "field": "group_count",
            "type": "ordinal",
            "label": "甲基群數"
          },
          "y": {
            "field": "site_count",
            "type": "quantitative",
            "label": "M1 位點數"
          }
        },
        "valueFormat": "number",
        "palette": {"kind": "categorical"}
      }
    ],
    "blocks": [
      {
        "id": "title",
        "type": "markdown",
        "body": "# HCC1395 M1 穩定甲基多群"
      },
      {
        "id": "finding",
        "type": "markdown",
        "body": "## 734 個穩定多群位點主要為二群\n\n長條圖比較 M1 位點依甲基群數的分布。"
      },
      {
        "id": "m1_chart",
        "type": "chart",
        "chartId": "m1_group_counts"
      }
    ]
  },
  "snapshot": {
    "version": 1,
    "generatedAt": "2026-07-23T00:00:00Z",
    "status": "ready",
    "datasets": {
      "m1_group_counts": [
        {"group_count": "2", "site_count": 651},
        {"group_count": "3", "site_count": 73},
        {"group_count": "4", "site_count": 9},
        {"group_count": "5", "site_count": 1}
      ]
    }
  },
  "sources": []
}
```

### 不可省略的 contract

- Top level 至少要有 `surface`、`manifest`、`snapshot`；HTML packaging 建議顯式保留 `sources` 陣列。
- `manifest.version=1`、`snapshot.version=1`。
- Report 的 `manifest.blocks` 必須有可呈現內容，而且至少一個 `type: "chart"` block。
- 第一個 markdown block 應是與 `manifest.title` 完全相同的 `#` title。
- Native chart 必須使用 `encodings.x.field` 與 `encodings.y.field` 或 `encodings.y.fields`；舊的 `xField`/`series` 對 artifact chart 會被拒絕。
- 每個顯示中的 native card/chart/table 必須有 inline `source` 或可解析的 `sourceId`。
- 目前 validator 要求圖表來源附「實際可執行 SQL」或 `.sql` source path。對 Python/TSV 結果，建議依既有 InterSubMod 作法：先寫入 SQLite staging table，執行 `SELECT`，將回傳 rows 與 snapshot exact-compare，再把執行 SQL 存入 `source.query.sql`。
- `snapshot.datasets` 是 `{dataset_id: [row, ...]}`，不是 `{columns, rows}`。
- 上限：50 datasets、每 dataset 2,000 rows、artifact payload 3,000,000 bytes、單一 cell string 4,000 chars。

## 3. Native chart 支援範圍

Plugin 的 shared chart contract 列出 18 種：

```text
area, bar, boxPlot, funnel, heatmap, histogram,
horizontalBar, horizontalStackedBar, horizontalStackedBar100,
leaderboard, line, pie, scatter, sparkline, stackedArea,
stackedBar, stackedBar100, waterfall
```

### 本任務的具體 encoding

| 圖 | Dataset grain | Encoding | 限制 |
|---|---|---|---|
| `bar` | 每個 group count 一列 | `x=group_count`, `y=site_count` | 數量軸從 0 開始 |
| `funnel` | 每個 gate stage 一列 | `x=stage`, `y=count` | 只適用有序、單 series、真實包含的 stages |
| `heatmap` | tidy long：每個 row×column cell 一列 | `x=row_read`, `color=column_read`, `y=distance` | 多個 color categories 會在 renderer 內 pivot 成 columns |

Native heatmap 範例：

```json
{
  "id": "read_distance_heatmap",
  "title": "Representative ALT-read distance matrix",
  "intent": "relationship",
  "question": "UPGMA 分群在 read distance 上是否呈現 block structure？",
  "rationale": "二維距離矩陣可直接觀察群內與群間距離塊狀。",
  "type": "heatmap",
  "dataset": "read_distance_cells",
  "encodings": {
    "x": {"field": "row_read", "type": "nominal", "label": "ALT read"},
    "color": {"field": "column_read", "type": "nominal", "label": "ALT read"},
    "y": {"field": "distance", "type": "quantitative", "label": "Bernoulli distance"}
  },
  "palette": {"kind": "sequential"}
}
```

30 reads 的完整 pairwise matrix 是 900 rows，低於 2,000-row limit；45 reads 是 2,025 rows，會超限。因此建議「每群以全部 reads 計算統計，但圖只顯示 deterministic medoid-nearest representatives」，並在 caption 明確標記。

## 4. UPGMA dendrogram 與 methylation heatmap 要如何放

### 4.1 為何不能全部用 native heatmap

Native heatmap renderer 遇到 missing/non-numeric cell 時，會將它代入為整張圖的最小值再著色。對 read×CpG 甲基資料，「CpG 未呼叫」不等於低甲基，所以這會把 NA 誤當 0-like signal。

因此：

- read×read distance 若為完整 finite matrix，可用 native heatmap。
- read×CpG methylation 有 NA 時，使用可重現的 PNG/SVG，NA 固定為灰色。
- UPGMA dendrogram 沒有 native chart type，應與 read×CpG heatmap 用同一 row ordering 生成一張複合圖。

### 4.2 HTML block 容器範例

```json
{
  "id": "upgma_methylation_figure",
  "type": "html",
  "sourceId": "focal_alt_read_matrix",
  "body": "<figure style=\"margin:0;font-family:system-ui,sans-serif\"><img src=\"data:image/png;base64,BASE64_BYTES\" alt=\"UPGMA dendrogram and ALT-read by CpG methylation probability heatmap; rows are focal-ALT reads, columns are CpGs, gray cells are uncalled CpGs\" style=\"display:block;width:100%;height:auto;border:0\"><figcaption style=\"margin-top:10px;font-size:13px;line-height:1.5\">UPGMA is a distance-based visualization of focal-ALT read methylation patterns, not a mutation phylogeny. Gray denotes an uncalled CpG; values were not imputed.</figcaption></figure>"
}
```

Custom chart contract：

- 圖必須由 reviewed rows/matrix 用可重現程式產生；不手畫 CSS bars、canvas 或 JavaScript chart。
- HTML block 只是 title/caption/alt/source context 與圖片容器。
- 圖片用 data URI 嵌入，因為 final `report.html` 必須 self-contained；不可依賴 sibling PNG 或 remote URL。
- HTML block 會進 sandboxed iframe；CSP 是 `img-src data: blob:`、`connect-src none`、`script-src none`。
- Data URI 會計入 3 MB artifact payload（base64 約增加 33% bytes）；建議限制圖的 pixel size、使用 indexed/optimized PNG，並在 packaging 前驗證 artifact bytes。
- 報告仍要有一個 native chart，例如 734 個位點的 group-count bar。

### 4.3 圖上應顯示的證據

UPGMA 圖的視覺不應單獨被稱為「統計顯著」。每個案例應在圖旁同時列出：

- focal ALT raw/core read 數。
- 每群 read 數與最小群大小。
- between/within distance ratio。
- column-shuffle null percentile/gate 結果。
- 10 seeds 的 modal K support。
- minimum assignment ARI。
- 圖中使用的 representative rule，以及統計是否用全部 core reads。
- 限制語：UPGMA 是 methylation-distance visualization，不是 mutation ancestry tree，也不單獨證明 cellular subclone。

## 5. 實際 deliver 命令

### 輸入與輸出

- Plugin root: `/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599`
- 預期輸入: `InterSubMod/research/20260723_hcc1395_singleton_alt_methyl_multigroup_visual_report/report/artifact.json`
- 預期輸出: `InterSubMod/research/20260723_hcc1395_singleton_alt_methyl_multigroup_visual_report/report/report.html`
- Failure-only screenshot: `InterSubMod/research/20260723_hcc1395_singleton_alt_methyl_multigroup_visual_report/report/portable_delivery_failure.png`

從 plugin root 執行 canonical npm command：

```bash
cd /bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599
```

```bash
CHROMIUM_EXECUTABLE_PATH=/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium-1223/chrome-linux64/chrome npm run report:deliver -- --input /big7_disk/liaoyoyo2001/InterSubMod/research/20260723_hcc1395_singleton_alt_methyl_multigroup_visual_report/report/artifact.json --output /big7_disk/liaoyoyo2001/InterSubMod/research/20260723_hcc1395_singleton_alt_methyl_multigroup_visual_report/report/report.html --ready-timeout-ms 15000 --action-timeout-ms 5000 --timeout-ms 60000 --screenshot /big7_disk/liaoyoyo2001/InterSubMod/research/20260723_hcc1395_singleton_alt_methyl_multigroup_visual_report/report/portable_delivery_failure.png
```

也可不切換 working directory，直接呼叫入口：

```bash
CHROMIUM_EXECUTABLE_PATH=/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium-1223/chrome-linux64/chrome node /bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/deliver_portable_artifact.mjs --input /big7_disk/liaoyoyo2001/InterSubMod/research/20260723_hcc1395_singleton_alt_methyl_multigroup_visual_report/report/artifact.json --output /big7_disk/liaoyoyo2001/InterSubMod/research/20260723_hcc1395_singleton_alt_methyl_multigroup_visual_report/report/report.html --ready-timeout-ms 15000 --action-timeout-ms 5000 --timeout-ms 60000 --screenshot /big7_disk/liaoyoyo2001/InterSubMod/research/20260723_hcc1395_singleton_alt_methyl_multigroup_visual_report/report/portable_delivery_failure.png
```

Builder 的 `--help` 已在本機實際執行，支援上述四個 options。Node/npm 環境為 Node `v22.22.1`、npm `10.9.4`。

## 6. Browser QA 可用性與已知問題

### 6.1 可用性

本機 Chromium 存在且可執行：

```text
/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium-1223/chrome-linux64/chrome
```

因此 builder 可進行完整 browser QA，不只是 `structural_only`。成功 receipt 應有：

```json
{
  "ok": true,
  "stages": {
    "validation": "passed",
    "package": "passed",
    "verification": "passed"
  }
}
```

`verification=passed` 代表同一 artifact 已被驗證：

- canonical payload 與 HTML embedded payload exact equality。
- desktop 1440 px 與 narrow 390 px。
- 預期 block/chart/table 數量與可見性。
- chart non-zero geometry。
- horizontal overflow。
- representative source menu/dialog 流程。
- 無 external request、browser error。

Screenshot 只應在失敗時保留；成功不需另寫一套 Playwright QA。

### 6.2 Linux classic-scrollbar 已知 false-positive

目前 plugin reader 在 Linux non-overlay scrollbar 環境有一個已重現的交付問題：top bar 使用 `width:100vw` 與 negative viewport margins，當長報告產生垂直 scrollbar 時，Chromium 可回報約 8 px 的 desktop horizontal overflow。

既有 singleton derivative v4 的 canonical attempts 已留下三種失敗：

- default 5 s reader timeout。
- 延長 timeout 後，`horizontal_overflow` at 1440 px。
- 30 s total budget 下 static-chart extraction timeout。

但同一 canonical artifact 使用只調整 top-bar CSS 的 scoped compatibility wrapper 後，完整 QA 通過：4 charts、2 HTML image blocks、4 tables、viewports 1440/390、source dialog passed，總時間約 16.6 s。

處理原則：

1. 先執行 canonical builder，不預先關閉 overflow gate。
2. 若失敗只是已知 top-bar `100vw` classic-scrollbar 問題，應用新 topic-scoped、可審計、CSS-only compatibility wrapper 或上游修正。
3. 相容性 receipt 必須證明 canonical payload/scientific content 未改，並重跑 1440/390 QA。
4. 不可因為 official verifier 誤報就直接宣稱 pass，也不可單純刪掉 overflow check。

## 7. 建議的 Step → Verify

1. 建立 canonical `artifact.json`
   → 驗證：`validate_artifact` 回傳 `ok=true`，且 counts 與 source dataset 相符。
2. 用 native bar/funnel/finite-distance heatmap，UPGMA + CpG NA heatmap 用 embedded static image
   → 驗證：每張圖都有鄰接解釋段、alt text、caption 與 source metadata。
3. 執行 canonical `report:deliver`
   → 驗證：receipt 為 `validation/package/verification=passed`，輸出 HTML 存在且非空。
4. 若出現已知 classic-scrollbar false-positive，只做 audited CSS compatibility correction
   → 驗證：1440/390 的 overflow、source dialog、content counts 全部 pass，artifact SHA-256 不變。
5. 結果交付
   → 驗證：用戶打開的是自我包含 `report.html`，不是 screenshot 或依賴 sibling assets 的頁面。

## 8. 查核來源

- Portable HTML shared contract: `/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/src/analytics-app-core.md`
- Report skill: `/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/SKILL.md`
- Custom report chart rule: `/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/specifications/mcp-app-report.md`
- Visualization routing: `/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/visualize-data/SKILL.md`
- Artifact JSON schema and validator: `/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/mcp/server.cjs`
- Native chart type source: `/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/src/analytics-app/charting/chart-contract.ts`
- Native heatmap renderer: `/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/src/analytics-app/charting/ChartRenderer.tsx`
- Portable builder: `/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/deliver_portable_artifact.mjs`
- Portable builder tests/examples: `/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/tests/portable-artifact-builder.test.mjs`
- Prior successful embedded-image example: `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/artifact.json`

