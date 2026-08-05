<!--
建立時間: 2026-07-18 23:45 +08:00
目標: 紀錄 positional-singleton ALT 甲基子結構報告的進行中設計決定、偏離、折衷與未決
處理範圍: 7 datasets 完整母體 + HCC1395 兩個 M2 PASS loci 深入視覺化
關聯檔案:
  - InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/00_INDEX.md
  - InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/pre-decision-audit.md
-->

# Implementation notes

狀態：`completed`

## 設計決定

- [決策] Task Type B：singleton 統計採全 7 datasets / 50,432 sites；HCC1395 採全 8,279 sites，不做 subset。
- [決策] 將使用者的 `LAT` 解讀為 `ALT`，並在報告第一屏明示。
- [決策] 分成兩層候選：M1 stable multigroup（734/8,279）與 M2 measured-axis clear residual partition（2/8,279）。
- [決策] heatmap 的生物分析使用全部 core reads；視覺呈現每群 4 個 deterministic
  representative reads（每圖 8 reads），避免 100+ 欄不可讀。
- [決策] assignment labels 在報告中改名 `Group A/B`，避免與 LongPhase-S `HP=1-1/1-2` 混淆。
- [決策] VAF 只作位點 burden/context，不用於甲基群發現，也不推斷父子或 clone order。
- [決策] HTML 以 canonical Data Analytics `artifact.json` 為 source of truth，再由 portable artifact builder 產生 self-contained HTML。
- [決策] portable reader 的 fallback header 與 runtime `.analytics-top-bar` 在 classic-scrollbar Chromium
  使用 `100vw` 時造成 8 px 文件溢位；以 scoped compatibility wrapper 僅把兩者改成 container width，
  不變更 artifact、數據、圖表或分析值。
- [決策] canonical portable renderer 不會把 native heatmap 匯出為彩色圖，只會保留數字表格；
  因此 6 個矩陣改用 schema-supported `type: html` blocks，產生 script-free、逐格著色、
  含 A/B 邊界、色階、NA 與可存取標籤的 heatmaps。

## 偏離

- [偏離] 原問題希望標示特殊 clone/subclone 狀態；因正式 second-marker cooccurrence、matched-normal、CN/purity/CCF
  尚未完成，本輪只標 `M2 residual epigenetic partition candidate`，不得升級為 confirmed clone/subclone。

## 折衷

- [折衷] 8,279 sites 的全量數字與狀態表會完整進報告；read-level heatmap 只對兩個 M2 PASS loci 展開，
  因它們是 HCC1395 唯二通過全部 measured-axis guardrails 的例子。
- [折衷] 每個位點的 distance、shared-CpG、methyl 三張 heatmap 使用 representative reads
  以確保可讀性；所有 group count、distance summary 與 M2 判斷仍由全量 108/109 core reads 計算。

## 未決

- [未決] G1/G2/R1 genetic corroboration、matched-normal、tumor-REF specificity、CN/CCF 與跨技術同位點 replication
  需等待正式 downstream release，不能在這份 singleton report 內補推。

## Lore

- Positional singleton 是「同 dataset/chrom 的 ≤50 kb transitive component 只有一個 sSNV」，
  不是「全基因組只出現一顆變異」，也不是「該 read 上沒有其他遠端 marker」。
- M2 PASS 的 2/8,279 是 operational yield，不是 subclone prevalence。
- `confirmed clone/subclone=0` 的原因是必要驗證尚未執行，不代表已證明這些位點沒有 subclone。

## 完成結果

- 全量 HCC1395 site audit：8,279/8,279 rows；gzip integrity PASS。
- M1 stable multigroup：734/8,279（8.8658%；evaluable 分母為 734/8,074 = 9.0909%）。
- M2 measured-axis clear residual partition：2/8,279（0.02416%；每 100,000 singleton 24.16）。
- 正式 cellular clone/subclone：0 個已建立；語意是必要的 genetic/normal/CN-lineage gates 尚未執行。
- HTML：39 blocks、8 native charts、6 colored heatmap HTML blocks、8 tables、4 metrics；
  source dialog、1440 px desktop、390 px mobile、classic-scrollbar overflow 與獨立 verifier 全部 PASS。
- Heatmap DOM QA：6/6，每張 64 cells，桌面／手機皆有實際色階、A/B 標記且無裁切。
- HTML SHA-256：`13307dd73ee26856e28e699c764de7281cb304ac5fa022baa781da30b0fde3e9`。
- QA 收據：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/html_qa_receipt.json`。
- Heatmap rendering QA：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/heatmap_rendering_qa.json`。

## Independent derivative sidecar v4（2026-07-19）

> 本節是 append-only 的獨立 sidecar ledger；不取代上方另一個 parallel implementation。其正式根目錄為
> `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/`。

### 固定輸入與定義

- Audit root：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/positional_singleton_methyl_multigroup_audit_v1_source_attested/`。
- v10 root：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/`。
- Source hashes：summary `6d4128be...aedc`；site audit `2d5d2479...90ce`；PASS cases
  `2fe42d55...1dd0`；v10 sites `a8871af3...f74`；stable assignments `82af076d...a4ba`。
- Positional singleton 是 dataset/chrom 內相鄰 gap `<=50,000 bp` 連接的傳遞 component，
  `component_size=1`；不是 read-sharing graph degree zero，兩者不保證等價。
- 所有統計使用全部 core reads；視覺每群固定顯示最多 15 個 deterministic medoid-nearest
  representatives。群名只顯示 `Group A/B`，不是 HP。
- Native chart 用於全樣本狀態、HCC1395 funnel 與兩張 read-distance heatmaps；因 NA methyl
  state 不能被 native heatmap 忠實表示，兩張 read×CpG methyl heatmaps 使用可重現 static PNG。

### 正式命令

Focused tests：

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python3.11 -m pytest \
  -p no:cacheprovider \
  research/20260718_singleton_alt_methyl_substructure_validation/tests/test_build_derivative_validation.py -q
```

Builder：

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python3.11 \
  research/20260718_singleton_alt_methyl_substructure_validation/scripts/build_derivative_validation.py \
  --audit-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/positional_singleton_methyl_multigroup_audit_v1_source_attested \
  --v10-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full \
  --claim-contract research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md \
  --output-dir research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4
```

Canonical validator：呼叫 `mcp__dataAnalyticsWidgets__validate_artifact`，完整 input 為
`results/derivative_validation_v4/artifact.json`；receipt 為
`results/derivative_validation_v4/canonical_validation_receipt.json`。

Portable delivery：

```bash
node research/20260718_singleton_alt_methyl_substructure_validation/scripts/deliver_portable_artifact_scrollbar_compat.mjs \
  --input research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/artifact.json \
  --output research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/report.html \
  --receipt research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/portable_builder_compat_receipt.json \
  --failure-screenshot research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/portable_builder_compat_failure.png \
  --ready-timeout-ms 15000 --action-timeout-ms 5000 --timeout-ms 30000
```

Browser QA：

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  research/20260718_singleton_alt_methyl_substructure_validation/scripts/qa_portable_report.py \
  --html research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/report.html \
  --artifact research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/artifact.json \
  --output-dir research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/browser_qa_v3 \
  --chromium /bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium-1223/chrome-linux64/chrome
```

### 實際輸出片段

- 全量：`469,849` sites；singleton `50,432`；M1 evaluable `48,347`；M1 flagged
  `5,961`；M2 `NOT_RUN=44,471 / NOT_EVALUABLE=5,913 / FAIL=18 / PASS=30`。
- 比例：`50,432/469,849=10.733661240%`；`30/50,432=0.059486041%` 是 observed
  operational yield；`30/48=62.5%` 只是在 determinate subset 的 conditional fraction。
- HCC1395：singleton `8,279`；M1 evaluable `8,074`；flagged `734`；M2
  `NOT_RUN=7,545 / NOT_EVALUABLE=732 / FAIL=0 / PASS=2`。
- Exact joins：chr14 `108/108=86+22`；chr22 `109/109=88+21`。兩個 core distance
  matrices 均 finite、symmetric、diag=0；read/methyl/CpG identity 全通過。
- Focused tests：`22 passed in 1.30s`。
- Canonical artifact validator：`ok=true`，datasets `8`，snapshot `ready`。
- Portable verifier：`23 blocks / 4 charts / 4 tables / 2 HTML blocks`；viewports
  `[1440,390]`；source dialog 與 keyboard interaction PASS。
- Browser QA：desktop `scrollWidth/clientWidth=1440/1440`；mobile `390/390`；
  console/page/request errors 均為空。
- Artifact SHA-256：`32183456381e8e0ec6fae5dfe19581ec7730c1909d565e92b3365b15b9649f7d`。
- HTML SHA-256：`31981e8f3ff374e26b27ea7e497f2c2b98f637a8a7d1187dcb3ab2e53adcc371`。

### Fail-closed evidence 與 claim ceiling

- `derivative_validation_v1/`：nullable `unstable=NA` schema gate；修正後要求其 NA pattern
  精確等於 M1 non-evaluable。
- `derivative_validation_v2/`：`chr1-22` 被過寬 `1-2` substring guard 誤判；修正為 token-boundary
  source-label guard。
- `derivative_validation_v3/`：canonical validator 拒絕無 actual SQL 的 native item provenance；
  v4 新增實際執行、persist 且 exact round-trip 的 SQLite staging queries。
- v4 初始 portable attempts 分別保留 reader timeout、classic-scrollbar overflow 與 full-Chromium
  timeout receipts；最終 compatibility delivery 仍重用 canonical builder/extractor/verifier，只修正 top bar width。
- M2 measured axes 不含 PS。`confirmed cellular subclone=0`、`linear ancestry=0` 是未完成必要
  downstream gates 的未驗證數，不是真陰性。M2 PASS 是 read-level residual epigenetic partition
  operational result，不是 biological/subclone prevalence 或 lower bound。
