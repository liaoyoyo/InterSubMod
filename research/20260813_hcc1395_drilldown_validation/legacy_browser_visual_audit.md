<!--
建立時間: 2026-08-13
目標: 以真實 Chromium 比較 2026-07-26 standalone 與 source-complete generated drilldown staging 的桌機/行動視覺及互動契約
處理範圍: legacy standalone、generated HCC1395 v3 summary staging；不修改 source、結果 receipt 或外部 bundle
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/scripts/run_legacy_browser_visual_audit.py
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/legacy_browser_visual_audit.generated.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/legacy_browser_visual_audit.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/06_legacy_standalone_desktop.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/07_legacy_standalone_mobile.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/13_current_generated_desktop.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/14_current_generated_mobile.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/15_current_generated_cooccur_desktop.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/16_current_generated_selfcheck_mobile.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/17_current_generated_detail_desktop.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/18_current_generated_methyl_detail_desktop.png
-->

用 SCQA + Before/After：**TL;DR：source-complete generated staging 已取得 `DIRECT_GENERATED` 證據：desktop/mobile document 寬分別維持 1440/390、browser errors 皆為 `[]`，co-occurrence 為 4 tables／106 denominator labels／0 fake cells；mobile sticky chrome 由歷史 overlay 的 204.91px 降至 86.19px，但 selfcheck 仍為 2 FAIL／0 SKIP，且本次是 `igv=skip, panels=0` 的輕量 UI QA，不是完整影像 bundle 驗證（影響：中，信心：高）。**

# Legacy standalone × current drilldown 視覺與互動 QA

> [!NOTE]
> **這是 UI contract comparison，不是兩份科學結果的 A/B。** 舊頁只展示 14 個 curated regions；目前 dashboard 的母體是 19,849 個 sSNV。比較的是閱讀、互動、responsive、authority 與 metric semantics，不比較兩頁的 biological values。

> [!CAUTION]
> **兩頁的 `region` 不是同一統計單位。** Legacy 的 `7,415 → 518 → 472 → 14` 是高度選案 funnel：14/472 = **2.97%**，14/7,415 = **0.189%**；這只描述「符合準則後又被挑入展示」的比例，不能當全資料 prevalence，也不能與 current topology 的 11,590 regions 對接。

> [!IMPORTANT]
> **Primary evidence status：`DIRECT_GENERATED`。**「目前版」直接載入 source-complete generated/staging URL，沒有 response/DOM overlay。Machine receipt `legacy_browser_visual_audit.generated.json` 是本文件數值判定的 single source of truth。

> [!WARNING]
> **資產範圍仍有限。** staging 使用完整 HCC v3 summary，但 `asset_policy` 是 `igv=skip, panels=0`。因此本次可驗證 summary、axis、co-occurrence、selfcheck、detail missing-state 與 responsive geometry；不可外推成 IGV／panel PNG、image alt、fit/zoom 或完整影像 bundle integrity 已通過。

## 結論先行

1. **Responsive：generated page 在兩種 viewport 都無 body overflow。** Desktop 的 viewport/client/document/body width 均為 1440；mobile 均為 390。兩次 browser error arrays 都是 `[]`。相較之下，legacy mobile 的 document/body width 是 555。
2. **Sticky：mobile chrome 已實質縮短。** Generated mobile 的 authority 為 0–33.531px，levelbar 為 34–86.188px，無 overlap；相較歷史 overlay 的總高度 204.91px，已降至 **86.19px**。
3. **Authority：presentation truth 與科學狀態同步。** Authority badge 與 selfcheck body 都呈現 `BLOCKED · 2 FAIL · 0 SKIP`；這代表 UI 忠實揭露失敗，而不是 scientific validation pass。
4. **Co-occurrence：generated DOM 已完整出現四個表格。** Receipt 記錄 4 tables、106 denominator labels、0 fake filter cells；denominator 與 NCHS 小樣本 caveat 均存在。
5. **Detail：全量母體與 drilldown path 成立。** 初始 scope/listCount 都是 19,849；點第一個 sSNV 後 detail 載入 `chr1:98,311`。本次選中 locus 顯示 `無 ISM`，故只驗證 missing-state，不驗證 panel image rendering。
6. **Accessibility/semantics：目前可見圖形有基本文字語意。** Generated desktop/mobile 的 role-img 都是 3/3 有 label，visible SVG 都是 2/2 有 title/label；denominator、missing/untested 與 legend token 都存在。由於 `panels=0`，visible images 為 0，不能宣稱 raster alt contract 已被覆蓋。
7. **甲基 detail 有正例與缺失例。** Runner 不只點第一個無 ISM 位點，也從 L1 axis encoding 機器選出 `chr1:1,320,793`（axis code 2）：105 reads、371 CpG，ALT/REF raw p=0.0010、effect field=0.029；HP 不顯著，其他 unavailable/invalid axes 保留第三態與 circularity 警示。

## 實測矩陣

| Surface | Viewport | clientWidth | document scrollWidth | body scrollWidth | Horizontal overflow | Sticky geometry after scroll | Browser errors |
|---|---:|---:|---:|---:|---|---|---|
| Legacy desktop | 1440×1000 | 1440 | 1440 | 1440 | No | filter bar 0–55px；無 stacked overlap | `[]` |
| Legacy mobile | 390×844 | 390 | **555** | **555** | **Yes** | filter bar 0–127.25px；無 stacked overlap | `[]` |
| Current generated desktop | 1440×1000 | 1440 | 1440 | 1440 | No | authority 0–39.875px；levelbar 40–98.594px；無 overlap | `[]` |
| Current generated mobile | 390×844 | 390 | 390 | 390 | No | authority 0–33.531px；levelbar 34–86.188px；無 overlap | `[]` |

Primary receipt 明示 `current_mode=generated`、`evidence_status=DIRECT_GENERATED`；結果來源：`InterSubMod/research/20260813_hcc1395_drilldown_validation/results/legacy_browser_visual_audit.generated.json`。

## Legacy：桌機清楚、行動版存在真實水平溢位

桌機上，legacy 的 narrative path 很好：claim ceiling → selection funnel → filter bar → curated case。四張 KPI 是 `7,415 有甲基註記的串聯區域 → 518 含 allele 軸 → 472 分支＋甲基符合準則 → 14 本頁展示`，不是 prevalence chain；最末 14 只占 eligible 472 的 2.97%。filter 把 `allele 軸位點 ≥` 從預設 3 改為 4 後，實際從 14/14 降為 7/14；第一張 card 能展開，SVG、表格與 paired heatmap caption 共置，適合人工審閱少量 canonical examples。

![Legacy standalone desktop：篩選後展開第一個 region](figures/06_legacy_standalone_desktop.png)

行動版的 typography 與 KPI 兩欄仍可讀，但長 code/path、展開內容或固定寬 geometry 把 document 撐到 555px。這不是局部可橫捲的 table，而是 body-level overflow；右側文字可能被使用者忽略。

![Legacy standalone mobile：390px viewport 下 document 寬達 555px](figures/07_legacy_standalone_mobile.png)

### Legacy semantics 實測

| Check | Actual result | Interpretation |
|---|---|---|
| body text | 15.5px / line-height 26.35px | 主文可讀性佳 |
| filter/detail | `14/14 → 7/14`；`detailOpen=true` | 控制與結果回饋直接 |
| visible role-img | 7；7 有 label | SVG 基本語意完整 |
| visible images | 42；missing alt=0 | raster fallback 有文字替代 |
| visible figcaptions | 21 | image context 不只靠檔名 |
| explicit denominator token | No | 雖有 cohort/selection counts，圖表旁缺統一「分母」標示 |
| missing/legend text | Yes / Yes | `未對齊` 與軸別圖例有文字，不只靠色彩 |
| visible elements <12px | 112 / 991 | 多為 caption/metadata；行動版仍可能偏小 |

`<12px` 是 computed-style DOM element count，不是 unique text-node count；只作 density signal，不當 accessibility pass/fail。

## Current generated（DIRECT_GENERATED）：source-complete UI contract

Generated staging 直接呈現完整 HCC v3 summary。Desktop 初始畫面固定 observation population `19,849 / 19,849`、claim ceiling、KPI 與 selfcheck BLOCKED 狀態；document width 維持 1440px。

![Current generated desktop：完整 summary、claim ceiling 與 BLOCKED authority](figures/13_current_generated_desktop.png)

390px 下 authority 收斂成 sample + BLOCKED badge，levelbar 緊接其下；document/body width 都維持 390px，且兩層 sticky 不重疊。

![Current generated mobile：390px 無 body overflow，sticky chrome 降至 86.19px](figures/14_current_generated_mobile.png)

### Detail flow

Playwright 點擊第一個 `.vrow` 後，detail 載入 `chr1:98,311`。Receipt 驗證代表圖形的 role/title/label、denominator、missing-state 與 legend token；清單總數仍是 19,849。選中 locus 顯示 `無 ISM`，因此此圖證明的是 missing-state 可辨識，不是 panel asset rendering。

![Current generated detail：代表樹、locus 軸與無 ISM missing-state](figures/17_current_generated_detail_desktop.png)

### ISM-bearing locus：數字、第三態與 claim ceiling 同頁

為避免 visual QA 只證明 missing-state，runner 由 `bootData.l1.axis` 選第一個具有非循環全域訊號的 locus，再用 `window.__DD.select()` 走正式 detail path。Machine receipt 保存 index、axis code、染色體與座標，並拒絕 methyl section 回退到「此位點沒有 ISM 資料」。`chr1:1,320,793` 顯示 105 reads、371 CpG；ALT/REF raw p=0.0010、effect=0.029，HP p=1.0；HP-fine、cluster、lineage 的未檢定/invalid 狀態使用斜線底紋與理由文字。這能驗證資訊階層與誠實語意，但單一探索 locus 不能外推 prevalence、因果或 cohort-wide significance。

![Current generated methyl detail：有資料正例、raw p/effect、未檢定第三態與 circularity 警示](figures/18_current_generated_methyl_detail_desktop.png)

### Co-occurrence semantics

Generated DOM 實測為 **4 tables、106 個 denominator labels、0 fake filter cells**。Receipt 同時找到 k unique table、explicit denominator 與 NCHS `n<30` caveat；因此這次可把 co-occurrence presentation contract 升級為 direct-generated pass，但不把表格內容外推成 scientific validation pass。

![Current generated co-occurrence：四個表格、分母與小樣本 caveat](figures/15_current_generated_cooccur_desktop.png)

### Selfcheck authority

Generated page 在 authority badge 與 selfcheck body 一致呈現：

```text
heading = 自檢未通過：2 項不成立
status  = 自檢 BLOCKED · 2 FAIL · 0 SKIP。有 2 項守恆等式不成立；此產物不可宣稱驗證通過。另有 0 項無法檢查。
```

![Current generated selfcheck mobile：BLOCKED、2 FAIL、0 SKIP](figures/16_current_generated_selfcheck_mobile.png)

### Current semantics 實測

| Check | Actual result | Interpretation |
|---|---|---|
| evidence mode | `generated` / `DIRECT_GENERATED` | 直接載入 staging，沒有 overlay |
| body text | 15px / line-height 23.25px | 主文可讀；比 legacy 稍緊 |
| scope/detail | 19,849/19,849；listCount=19,849；detail=`chr1:98,311` | 全量母體與 detail flow 成立 |
| visible role-img | desktop/mobile 皆 3；3 有 label | 可見圖形有基本文字語意 |
| visible SVG | desktop/mobile 皆 2；2 有 title/label | tree/locus 不只靠視覺 |
| explicit denominator token | Yes；co-occurrence labels=106 | 分母長駐於 scope、charts 與 detail |
| missing/legend text | Yes / Yes | missing、untested 與圖例不是只靠色彩 |
| visible elements <12px | desktop 222/811；mobile 203/750 | metadata 密度仍高，適合作後續 typography audit |
| co-occurrence completeness | 4 tables；106 denominator labels | Generated presentation contract 成立 |
| co-occurrence pseudo interaction | 0 fake filter cells | no-fake-action gate 通過 |
| selfcheck | 2 FAIL；0 SKIP | UI 正確封鎖 validation claim |
| browser errors | desktop/mobile 皆 `[]` | 本次 navigation 未記錄 browser errors |

### 輕量 staging 的範圍限制

本次 staging 雖使用完整 HCC v3 summary，但資產政策是 `igv=skip, panels=0`，receipt 中 visible images 為 0。這是刻意縮小的 **light QA**：可確認 summary、axis table、ISM-bearing 與 no-ISM detail、co-occurrence、selfcheck、desktop/mobile layout 與 sticky geometry；不可確認 IGV snapshot、methylation panel PNG、image alt、local fit/zoom、lazy image loading 或完整 asset bundle 的檔案完整性。若要對外宣稱完整瀏覽器產物已驗證，必須另建含影像資產的 durable staging 再跑同一套 generated-mode QA。

### Historical overlay（PARTIAL；僅作診斷歷史）

較早的 `legacy_browser_visual_audit.json` 記錄 `current_mode=overlay`、`evidence_status=PARTIAL`：co-occurrence 只有 1 table／2 denominator labels，selfcheck 是 1 FAIL／1 SKIP，mobile sticky chrome bottom 為 204.906px。它已被本次 direct-generated 主判定取代；08–12 screenshots 僅保留為 source-overlay 診斷紀錄，不可拿來代表目前 generated page。

## 哪些 legacy 設計值得移植

### P0 — 應移植

1. **頁首 claim ceiling。** Legacy 的「本頁是描述性瀏覽，不是證據」比單一 BLOCKED badge 更直接。新版應讓 claim ceiling 與 validation status 同時存在：前者說能下什麼結論，後者說 gates 是否通過。
2. **Selection funnel + display count。** Legacy 把 7,415 → 518 → 472 → 14 與 selection recipe 寫在同一路徑。新版對 curated canonical/outlier gallery 也應保存 source population → axis eligible → branch eligible → displayed 四層分母，並明示 14/472 與 14/7,415，而不是暗示 prevalence。
3. **案例卡的證據共置。** 一個 legacy region 內將 locus SVG、per-sSNV axis table、linkage table、paired heatmap/caption 放在一起，適合作為新版 detail 的「evidence pack」子區塊；移植時必須重新綁定 current topology region ID，不能沿用舊 region 語意。
4. **每張 raster 都有 caption/alt。** 這個 contract 應保留到 panel lazy loading，不因資料量變大而退化。

### P1 — 有條件移植

1. **Compact KPI strip。** 可作 canonical/outlier gallery 的入口，但每張 KPI 必須附 numerator/denominator、scope、source 與 observation-only ceiling。
2. **簡潔 filter bar。** Legacy 的 2–4 個 select + 顯示 n/N 很適合 curated gallery；新版完整 dimensions 仍應保留，但 mobile 可先顯示 active/常用 filters，其餘折疊。
3. **Expand/collapse。** 適用於 10–20 個 curated cases，不適合 19,849 loci 全量頁；全量仍需 virtual list + detail pane。

## 哪些 legacy 設計不應移植

1. **13 MB base64 standalone 作全量交付。** 14 cases 已達 13 MB；擴到 thousands of loci 會破壞載入、cache 與 provenance。大型 panel 應 lazy/chunked，並保存 asset hash。
2. **固定寬 SVG/table 讓 body 溢位。** Mobile 555/390 是 hard failure；寬圖只能在明確 `.panel-scroll` 局部橫捲，document width 必須鎖在 viewport。
3. **只靠 footer 補 metric scope。** denominator、missing reason、selection criterion 要靠近 chart，不應要求讀者捲到底才知道母體。
4. **用 curated 14 cases 代表 cohort。** Curated examples 只能標 canonical/extreme/well-explained，不可代替 aggregate 或跨樣本分布。
5. **V 值只有欄名、沒有完整定義。** `V_allele` / `V_hp` 必須附 effect/test definition、unit、validity 與 multiplicity family；顏色/顯著 badge 不能代替 metric contract。

## Current UI 還可改善的地方

1. **把完整影像 bundle QA 補成下一個 hard gate。** 現在的 `igv=skip, panels=0` 不覆蓋任何 raster asset。應在新的 durable staging 產生 IGV/panels，再檢查 image load failure、alt、fit/zoom、局部 overflow、lazy loading 與 asset provenance；不可用本次 light QA 代替。
2. **降低小字與 metadata 重複。** Desktop 222/811、mobile 203/750 visible DOM elements 使用 <12px。優先把判讀必要文字升到 ≥12px；來源 badge/monospace ID 可保留較小，但需足夠 contrast。
3. **Mobile filter progressive disclosure。** 初始只顯示 scope、active filters 與 2–3 個高價值 dimensions；完整 dimensions 放在 collapsed advanced section，避免第一個 locus 清單遠離首屏。
4. **保留縮短後的 sticky contract。** Generated mobile 已從歷史 overlay 的 204.91px 降至 86.19px，且沒有 overlap；後續改版應把 `authorityBottom≈33.53`、`levelbarBottom≈86.19` 當 regression baseline，而不是再擴回多列 chrome。
5. **Detail 的 reading order 再壓縮。** Tree/locus 的 denominator 與 caveat 很完整，但 desktop 一次出現多層 legend；可把第一層保留判讀必需項，其餘放 `<details>`，不刪除 provenance。
6. **保留 scale/validation 架構，不退回 standalone。** Current 的全量 virtual list、capability/selfcheck authority 與 provenance affordance 才能支撐 19,849 sSNV；legacy 應只貢獻 funnel、折疊案例卡、paired heatmap 與簡潔 filter bar。影像 lazy-loading 的實際品質仍待 full-asset build 驗證。

## 執行與可重現證據

### 輸入路徑

- Legacy HTML：`InterSubMod/docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/20260726_ssnv_branch_x_methyl_browser.standalone.html`
- Current generated staging：`/tmp/dd-hcc1395-v3-light-final2.pUtxKy`（**ephemeral QA input；不是 durable deliverable**）
- Asset policy：完整 HCC v3 summary；`igv=skip, panels=0`
- Runner：`InterSubMod/research/20260813_hcc1395_drilldown_validation/scripts/run_legacy_browser_visual_audit.py`

### 執行命令

```bash
python3 /bip7_disk/liaoyoyo2001/.codex/skills/webapp-testing/scripts/with_server.py \
  --server "python3 -m http.server 8767 --directory /tmp/dd-hcc1395-v3-light-final2.pUtxKy" \
  --port 8767 -- \
  python3 research/20260813_hcc1395_drilldown_validation/scripts/run_legacy_browser_visual_audit.py \
    --legacy docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/20260726_ssnv_branch_x_methyl_browser.standalone.html \
    --current-mode generated \
    --current-url http://127.0.0.1:8767/index.html \
    --figures research/20260813_hcc1395_drilldown_validation/figures \
    --result research/20260813_hcc1395_drilldown_validation/results/legacy_browser_visual_audit.generated.json
```

上述 `/tmp/dd-hcc1395-v3-light-final2.pUtxKy` 是執行當下的 staging 絕對路徑，可能被系統清理，不應出現在對外交付清單或被視為可長期重跑的 artifact。可重現依據是 runner、generated receipt 與 13–18 screenshots；若要長期重跑，應先把等價 generated bundle 建到受控、具 provenance 的 durable output 路徑。

### 輸出路徑

- Primary machine-readable receipt：`InterSubMod/research/20260813_hcc1395_drilldown_validation/results/legacy_browser_visual_audit.generated.json`
- Current generated screenshots：`InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/13_current_generated_desktop.png` 至 `18_current_generated_methyl_detail_desktop.png`
- Historical overlay receipt：`InterSubMod/research/20260813_hcc1395_drilldown_validation/results/legacy_browser_visual_audit.json`（PARTIAL；非主判定）
- Legacy screenshots：`InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/06_legacy_standalone_desktop.png`、`07_legacy_standalone_mobile.png`
- 本文件：`InterSubMod/research/20260813_hcc1395_drilldown_validation/legacy_browser_visual_audit.md`

### 實際輸出片段（exit 0）

```text
legacy_desktop  viewport=1440x1000  scrollWidth=1440  overflow=false  errors=[]
legacy_mobile   viewport=390x844    scrollWidth=555   overflow=true   errors=[]
current_desktop viewport=1440x1000  scrollWidth=1440  overflow=false  errors=[]
current_mobile  viewport=390x844    scrollWidth=390   overflow=false  errors=[]

current_mode=generated, evidence_status=DIRECT_GENERATED
cooccur tables=4, denominatorLabels=106, fakeFilterCells=0
selfcheck fail=2, skip=0
methyl QA locus=chr1:1,320,793, axisCode=2, reads=105, CpG=371
ALT/REF raw p=0.0010, effect field=0.029; cluster/lineage remain untested or invalid

current mobile sticky:
authorityBottom=33.53125, levelbarTop=34, levelbarBottom=86.1875, overlap=false
```

## 最終判定

- **Legacy standalone**：desktop 的 curated evidence-card pattern 值得保留；mobile responsive 是 **REVISE**，不可原樣移植。
- **Current generated presentation/interaction**：**DIRECT_GENERATED PASS within light-QA scope**。Desktop/mobile geometry、sticky、authority、scope/detail、4-table co-occurrence semantics 與 no-fake-action gate 都有直接 browser receipt；mobile chrome 已由歷史 overlay 的 204.91px 降至 86.19px。
- **Full image bundle**：**NOT EVALUATED**。`igv=skip, panels=0` 使 IGV/panel PNG、image alt、fit/zoom、lazy-loading 與 bundle integrity 都不在本次 coverage；必須另建 full-asset durable staging 才能評估。
- **Scientific status**：**BLOCKED**。Selfcheck 是 2 FAIL／0 SKIP，頁面也明示不可宣稱驗證通過。`DIRECT_GENERATED` 只表示證據取得方式是直接 generated page，不等於研究結論通過。
- **Historical source overlay**：維持 **PARTIAL**，只作診斷歷史；不得覆蓋或降低上述 generated receipt 的主判定權重。
