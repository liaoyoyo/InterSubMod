<!--
建立時間: 2026-08-13
目標: 定義多 BAM 單頁 dashboard 的 UX、資訊階層、跨樣本 metric contract 與 browser QA acceptance
處理範圍: 7 份 topology datasets（6 biological samples + 1 HCC1395 technical replicate）；已實作 bounded offline snapshot，production/data science 仍維持 PARTIAL
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_HCC1395_drilldown完整驗證與多樣本改進_01.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/cohort_topology_metrics.csv
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/13_current_generated_desktop.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/14_current_generated_mobile.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/15_current_generated_cooccur_desktop.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/16_current_generated_selfcheck_mobile.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/17_current_generated_detail_desktop.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/18_current_generated_methyl_detail_desktop.png
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_artifact.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_multi_bam_analysis_overview.standalone.html
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/multi_bam_dashboard_browser_qa.json
-->

用 SCQA + decision-first hierarchy：**多 BAM dashboard 應預設呈現 6 個 biological samples／7 個 datasets 的 topology 可比範圍，5 秒內讓讀者看懂「technical identity/receipt 7/7 通過，但 dataset-level 全 family complete 與全 unit objective-certified 皆為 0/7」；technical replicate、PARTIAL/ABSENT、ISM/lineage 與 single-sample selfcheck 必須分層隔離，不能被平均或借用（影響：高，信心：高）。**

# 多 BAM 單頁 Dashboard UX 與 Browser QA Spec

> [!IMPORTANT]
> **狀態：IMPLEMENTED SNAPSHOT / DATA-PRODUCT QA PASS / SCIENCE PARTIAL。** Bounded artifact 與 standalone HTML 已完成；40/40 Playwright assertions、canonical packaging/verification 與 1440/1024/512/390/320 CSS-pixel 檢查通過。此狀態不代表 production route、truth benchmark 或多樣本 ISM/lineage bundles 已完成。

> [!WARNING]
> **Task Type B — Comprehensive validation。** Default scope 是 `cohort_topology_metrics.csv` 全 7 列，不得用 HCC1395 或少數 BAM 當 cohort 替身。Multi-sample science status 仍是 **PARTIAL**：topology 可逐 dataset 比較，等價的 7-dataset ISM/lineage bundles 不存在。

## 1. 目標、受眾與決策

### 1.1 服務目標

- **G3**：讓 read-level epigenetic 結果只在具有效 sample-specific evidence 時出現，並保留 untested／invalid／circularity 警示。
- **G4**：把 biological sample、dataset、technical replicate、capability 與 denominator 做成明確 contract，支援跨樣本一致性檢查。
- **G5**：讓每個畫面狀態都可由 receipt、資料表與 browser QA 外部稽核；錯誤時 fail closed。

### 1.2 Primary audience

| Audience | 需要在頁面完成的決策 | 不應被頁面誘導的決策 |
|---|---|---|
| PI／研究 lead | 哪個 biological sample 是 coverage/identifiability outlier；下一個該查哪個 sample/locus；目前結論 ceiling | 把 unique tree rate 當 clone prevalence；由單樣本 ISM 外推 cohort |
| 分析研究員 | 比較同口徑 topology metrics；確認 technical replicate；定位 PARTIAL/ABSENT 原因 | 把不同分母的百分比混在同圖；把 missing 當 0 |
| QA／外部 reviewer | 驗證 sample identity、hash、receipt、schema、denominator 與 selfcheck | 把 technical pipeline PASS 當 biological validation PASS |

### 1.3 頁面必須支援的三個問題

1. **整體狀況如何？** 六個 biological samples、七個 datasets；哪些 technical gates 通過，哪些 scientific/completeness gates 沒通過。
2. **哪個 dataset 值得 drill down？** 先看同口徑 tree coverage 與 unique rate 分布，再選 sample；不由 raw count 大小排序好壞。
3. **選定 sample 後能看什麼？** 只呈現該 sample 通過 identity/capability gate 的 topology、ISM、LCA/lineage 與 assets；其餘標 PARTIAL、ABSENT、UNTESTED 或 INVALID。

## 2. Default view：5 秒內要回答的事

Default 是 **All biological samples**，不是任一 BAM，也不是 HCC1395。首個 1440px viewport 不捲動時必須同時看見以下資訊：

1. `6 biological samples · 7 datasets · 1 technical replicate`。
2. `hash match 7/7 · sample identity 7/7 · topology receipt PASS 7/7`，旁邊直接寫「technical gate，不是 biological validation」。
3. `dataset-level all mutation-bearing families complete 0/7` 與 `all units objective-certified 0/7`；文案不可縮成「0 family complete」。
4. 兩個分開的 cohort plots：tree coverage 與 unique rate among tree。H2009 的 64.17%／35.29% 與 HCC1954 的 72.90%／90.58% 能被快速辨識，但不自動貼上 biological good/bad 標籤。
5. HCC1395_DORADO 必須顯示 `technical replicate of HCC1395`，且不納入 primary biological-sample macro summary。
6. 全頁 claim ceiling：`descriptive topology comparison; not clone prevalence, ancestry, lineage, or caller accuracy`。

### 2.1 主次資訊階層

| Layer | 內容 | Default | 收合規則 |
|---|---|---|---|
| L0 Authority | scope、selected sample/dataset、technical/scientific status、claim ceiling | Sticky；永遠可見 | 不可整段收合；mobile 只保留 selector + 最高嚴重度文字狀態 |
| L1 Decision | cohort gate strip、兩個 topology comparison plots、outlier annotation | 展開 | 不收合；這是 5 秒判讀核心 |
| L2 Coverage | sample × capability matrix、denominator-safe comparison table | 展開 | mobile 可收合 table，但 matrix 的 PARTIAL/ABSENT summary 常駐 |
| L3 Investigation | selected-sample funnel、axis summary、co-occurrence、locus list/detail | All-samples 時顯示引導；選 sample 後展開 | 依 module 收合；錯誤／BLOCKED 摘要不可藏 |
| L4 Audit | selfcheck details、receipt、hash、schema/model/AF basis、asset manifest | 摘要常駐，details 收合 | FAIL count、SKIP count、provenance dirty/mismatch 不可收合 |

## 3. Sample selector contract

### 3.1 資訊模型

- 第一層選擇是 **biological sample**：All、COLO829、H1437、H2009、HCC1395、HCC1937、HCC1954。
- 第二層只在有多 dataset 時出現。HCC1395 下列 `HCC1395 (primary)` 與 `HCC1395_DORADO (technical replicate)`；不得把兩者平均成一列，也不得讓 DORADO 增加 biological n。
- Selector label 同時顯示 `biological sample` 與 `dataset`。只寫「sample」不足以區分 HCC1395 與 DORADO。
- URL state 使用明確的 `biological_id` + `dataset_id`；未知、重複或不相容組合不得猜測。

### 3.2 行為與 fail-closed

| Interaction | 正常行為 | Fail-closed 行為 |
|---|---|---|
| 初次開頁 | 選 All biological samples；只顯示 cohort-comparable topology | 不以第一列或 HCC1395 作隱性預設 |
| 選 biological sample | 更新 capability、funnel、list、detail 與 audit context；清除前一 dataset locus | 任一 module 尚未完成 identity check 時顯示 `LOADING/VERIFYING`，不保留舊數值 |
| 選 HCC1395_DORADO | 顯示 technical replicate badge、paired-to HCC1395 關係與 macro exclusion | 若 paired identity 缺失，標 PARTIAL；不自動當獨立 biological sample |
| 切換 tab/view | 保留目前 biological/dataset context，更新 URL | module 不支援 All-samples 時禁用並附原因，不偷偷切 HCC1395 |
| 載入 PARTIAL capability | 顯示已載 numerator/expected denominator、缺項原因、可讀 claim | 不計入 cohort aggregate，不把可讀部分升格 AVAILABLE |
| 載入 ABSENT capability | 顯示 em dash + `ABSENT: no valid sample-specific source` | 不顯示 0，不借用其他 sample，不殘留上一個 sample 圖形 |
| invalid URL/sample | 回到 All-samples 並顯示可聚焦 error summary | 不猜 alias；不載入最接近名稱的 dataset |
| shard/network/parse error | module 進 `ERROR/BLOCKED`，保留 selector 與 retry；記錄 source path/token | 清空該 module 舊 DOM/data；其他 module 仍可用但不得把全頁標 PASS |
| 快速連續切 sample | 只接受最後一次 request token 的 response | 丟棄 late response，防止 cross-sample state contamination |

## 4. 建議 chart / table 順序

| Order | Module | 視覺形式 | Metric contract |
|---:|---|---|---|
| 1 | Authority + selector | 單列 sticky bar + 文字 status | sample/dataset/scope/claim ceiling；狀態不只靠顏色 |
| 2 | Cohort gate strip | 四張 compact cards | 7/7 hash、7/7 identity、7/7 receipt PASS、0/7 dataset-level complete/certified gates；technical 與 scientific 分列 |
| 3 | Topology comparison | 兩張上下對齊 dot plot，共用 sample row order | A=`tree_n / region_n`；B=`unique_tree_n / tree_n`。不可用 dual axis，不把 B 除以 all regions |
| 4 | Capability matrix | sample rows × topology/ISM/LCA/lineage/assets columns | AVAILABLE/PARTIAL/ABSENT/UNTESTED/INVALID；每格附 reason/denominator link |
| 5 | Denominator-safe table | 可排序 data table；預設 biological sample，technical replicate nested | 顯示 numerator、denominator、rate、schema/model/AF basis；排序只改視圖，不改 cohort |
| 6 | Selected-sample observation funnel | numerator/denominator funnel | 只有選 dataset 後出現；每一步標 source、scope、missingness；不可稱 prevalence |
| 7 | Axis summary | 分軸小 multiples + validity labels | ALT/REF、HP、HP-fine、cluster、lineage 分開；不做單一「有訊號率」總分 |
| 8 | Co-occurrence | 四表摘要，先 summary 再 details | denominator 貼近表；小樣本 suppress rule 常駐；格子非 filter 時不得做成可點 affordance |
| 9 | Locus list + detail | virtual list + detail pane | 切 sample 清空 selection；missing locus 與具 ISM locus 都有明確第三態 |
| 10 | Selfcheck + provenance | 常駐 fail summary + collapsible evidence | FAIL/SKIP/dirty/mismatch 置頂；hash/schema/path 放 details；失敗不能被 accordion 隱藏 |

### 4.1 Cohort chart 的 sample order

- Default 用固定 biological order，HCC1395_DORADO 緊鄰 HCC1395 並縮排；避免每次切 metric 後 row 跳動。
- 可讓使用者依 tree coverage 或 unique rate 排序，但 technical replicate 仍黏在 HCC1395 下方，且 macro exclusion badge 保留。
- 兩張 dot plot 使用同一 0–100% scale、同一 row order，但標題與 denominator 完整不同；不可疊成一個 composite score。
- 圖下必有等價 data table，不能要求讀者由像素估值。

## 5. 哪些指標可以／不可以跨樣本

### 5.1 可作逐 dataset 描述性比較

僅限 CSV 已證明 schema、model、AF basis、hash 與 sample identity 可比的 topology dataset：

- `tree_pct_all_regions = tree_n / region_n`。
- `unique_pct_among_tree = unique_tree_n / tree_n`。
- `region_n`、`distinct_ssnv_n`、`tree_n`、`unique_tree_n` 可並列作 workload/denominator context，但不能當品質排名。
- biological-sample macro summary 若要顯示，只能用六個 primary biological datasets 的 sample-level distribution；HCC1395_DORADO 另作 paired sensitivity，不加入 n。

### 5.2 禁止直接 pool、平均或排名

| 指標／狀態 | 為何不能直接跨樣本 | Dashboard 行為 |
|---|---|---|
| raw `region_n`、sSNV/read/CpG counts | exposure、coverage、sample biology 與 grouping 數量不同 | 只作 denominator/context；不畫「越多越好」rank |
| `unique_pct_all_regions` vs `unique_pct_among_tree` | 分母不同，數值不可互換 | UI 名稱包含 denominator；禁止同一 metric toggle 混用 |
| ISM linkage／availability／axis p/q/effect | 目前沒有 7 份等價 ISM bundle；tested set 與 missingness 不同 | 僅 selected-sample 顯示；All-samples 為 PARTIAL/ABSENT matrix，不算 cohort rate |
| raw p≤0.05 或 BH q≤0.05 counts | test family、eligible denominator、coverage 與 multiplicity scope 可能不同 | 不設跨樣本排行榜；必須先鎖定同一 test family 才能另開 cohort analysis |
| cluster-axis 或 self-clustering result | circularity；不能作獨立支持 | 標 INVALID/EXCLUDED，不納入 summary KPI |
| LCA/lineage numbers | source identity／pre-post scope 不等價，且 capability 不是全樣本可用 | sample-only；ABSENT 不補 0，不跨 sample pool |
| `receipt_all_pass`、hash/identity flags | technical integrity gate，不是效果量或 biological score | 分列 status；不得算「研究品質百分比」 |
| HCC1395 vs HCC1395_DORADO | 同 biological sample 的 technical datasets | paired sensitivity；不得當 n=2 biological evidence |
| legacy A/B、Cramér's V、raw coread | window、taxonomy、axis/denominator 與 current topology 不等價 | 不進 multi-BAM dashboard comparison layer |

## 6. Status semantics：文字先於色彩

| Status | 精確語意 | 顯示要求 | Aggregate 規則 |
|---|---|---|---|
| `AVAILABLE` | sample-specific source 存在、identity/schema gate 通過，module 可讀 | 寫出 AVAILABLE + source scope；不可寫「validated」 | 只在 metric contract 相同時入分母 |
| `PARTIAL` | 預期資料只有一部分可讀，或 identity/provenance/guard 有未閉合項 | 顯示 `loaded/expected`、缺項 reason、claim ceiling | 預設排除 cohort rate；需另列 partial count |
| `ABSENT` | 沒有有效的 sample-specific source | 顯示 em dash、ABSENT 與原因 | 不補 0；不入 numerator 或 denominator |
| `UNTESTED` | 有資料但此 metric 未執行合法檢定 | 顯示斜紋/文字 `未檢定` 與原因 | 不等於不顯著；不入 significant rate denominator |
| `INVALID` | 計算存在但因 circularity/identity/definition 不可作證 | 顯示 INVALID + 阻擋理由 | 永遠不入 conclusion KPI |
| `BLOCKED` | 一個或多個 hard gates 失敗，禁止宣稱驗證通過 | authority 常駐 FAIL/SKIP count；details 可展開 | cohort status 取最高嚴重度，不平均成黃色 |

所有狀態必須同時有：**完整文字 label + 非色彩視覺差異（border/pattern/shape）+ 可讀原因**。Color 只作冗餘提示；screen reader accessible name 必須包含 status、sample、module 與 reason。`PARTIAL` 不可用淡綠，`ABSENT` 不可顯示為數值 0，`UNTESTED` 不可顯示為「不顯著」。

## 7. Progressive disclosure 與 collapse details

### 永遠展開

- sample/dataset selector、cohort scope、claim ceiling。
- 最高嚴重度 status；FAIL、SKIP、dirty、identity mismatch 數量。
- Cohort gate strip、兩個 topology comparison plots。
- Selected-sample 的 numerator/denominator funnel 與 PARTIAL/ABSENT summary。

### 預設展開，但可由使用者收合

- Capability matrix、denominator-safe table。
- Co-occurrence summary 的四表標題、分母與 verdict；長表 body 可收合。
- Locus tree/locus spacing、read-state、methyl-axis summary。

### 預設收合

- 完整 receipt JSON、SHA-256、source path、build argv、schema/model details。
- 次要 legend、所有 per-row provenance 與完整 selfcheck equation trace。
- 大型 IGV/panel asset（未載入前顯示明確 load action、尺寸與 provenance）。

Collapse 不得改變資料或篩選狀態；展開區需用原生可達語意（例如 button + controlled region），支援 keyboard，並宣告 expanded/collapsed。錯誤摘要不得只存在於收合內容。

## 8. Mobile 與 accessibility contract

1. **Reflow**：320px、390px 均不得有 body-level horizontal overflow；寬表只允許具 label 的局部 scroll container。
2. **Sticky budget**：390px 下 authority + levelbar/selector 總高度目標不超過 96px；目前 single-sample baseline 是 86.19px。scroll 後只保留 sample/dataset、最高嚴重度 status 與 view control。
3. **Selector**：使用具 label 的 combobox/listbox；選項文字完整區分 biological sample、primary dataset 與 technical replicate。不可只靠縮排傳達關係。
4. **Reading order**：DOM 順序與視覺順序一致：authority → cohort status → comparisons → capability → investigation → audit。
5. **Keyboard/focus**：selector、tabs、sort、details、locus rows、load asset 都能鍵盤操作；切換 sample 後 focus 到更新後的 scope heading，不跳到頁首。
6. **State announcement**：sample change、load error、PARTIAL/ABSENT 與 filter result count 經 live region 簡短宣告；不逐列洗屏。
7. **Charts**：每張 chart 有 programmatic name、denominator 描述與等價 table；不以 hover 作唯一數值入口。
8. **Color independence**：status、significance、missing/untested 使用文字 + pattern/border；contrast 與灰階列印皆能分辨。
9. **Target/text**：主要 controls 的可點區至少 44×44 CSS px；主文目標 16px，判讀必要文字不得小於 12px。
10. **Zoom/motion**：200% zoom 仍保持 reading order 與功能；尊重 reduced-motion，不用 animation 表達唯一狀態變更。

## 9. 既有畫面的可保留與需調整部分

### 9.1 保留：authority、claim ceiling 與 mobile 收斂

Figures 13/14 已證明 current generated page 在 desktop/mobile 能讓 BLOCKED authority、scope 與 observation-first hierarchy共存；multi-BAM 版沿用其視覺語言，但 selector 必須升級為 biological sample + dataset 二層。

![Current generated desktop：可沿用 authority、claim ceiling 與 observation hierarchy](figures/13_current_generated_desktop.png)

![Current generated mobile：sticky chrome 86.19px，可作 multi-BAM regression baseline](figures/14_current_generated_mobile.png)

### 9.2 保留：分母鄰近資料；強化 sample context

Co-occurrence 四表與 selfcheck BLOCKED 都已有 direct-generated evidence。多 BAM 版需要在每個表格標題旁再加 biological/dataset context；`2 FAIL / 0 SKIP` 是當次 HCC1395 generated bundle 狀態，不得顯示成 cohort-wide count。

![Current generated co-occurrence：四表、分母與小樣本 caveat](figures/15_current_generated_cooccur_desktop.png)

![Current generated selfcheck：失敗狀態不可被折疊或平均](figures/16_current_generated_selfcheck_mobile.png)

### 9.3 同時保留 missing-state 與 evidence-present state

Figure 17 的 `無 ISM` 與 Figure 18 的 105 reads／371 CpG、global ALT/REF 探索關聯是兩個必要 QA state。多 BAM 版不得只測其中之一；Figure 18 的 cluster/lineage invalid/untested 文案也必須保留。

![Current generated detail：合法的無 ISM missing-state](figures/17_current_generated_detail_desktop.png)

![Current generated methyl detail：有資料、未檢定與 circularity 可同時辨識](figures/18_current_generated_methyl_detail_desktop.png)

## 10. Browser QA acceptance

### 10.1 Step → Verify

1. **載入 All-samples default** → 驗證：不點擊、不捲動即可讀到 6 biological／7 datasets／1 replicate、三個 7/7 technical gates、兩個 0/7 dataset-level completeness gates與 claim ceiling。
2. **驗證 cohort charts** → 驗證：7 dataset rows 全出現；HCC1395_DORADO 緊鄰 HCC1395 且標 technical replicate；tree coverage 與 unique-among-tree 使用不同標題/分母，沒有 composite score。
3. **逐一切換七個 datasets** → 驗證：authority、capability、funnel、list/detail、selfcheck 的 dataset ID 一致；每次切換先清舊 DOM，late response 不能回填。
4. **測 PARTIAL/ABSENT/UNTESTED/INVALID fixtures** → 驗證：文字、reason、denominator 與非色彩 pattern 正確；ABSENT 不產生 0；UNTESTED 不出現「不顯著」。
5. **測 technical replicate** → 驗證：HCC1395_DORADO 顯示 biological ID=HCC1395、technical replicate badge、macro excluded；不得使 biological n 從 6 變 7。
6. **測 detail 正負狀態** → 驗證：一個無 ISM locus 與一個有 ISM locus 都能載入；有資料 locus 保留 untested/invalid/circularity caveat。
7. **測 error recovery** → 驗證：invalid sample、404 shard、malformed row、identity mismatch、快速切換 race 都封鎖該 module且清除 stale values；修復後 retry 不需 reload 全頁。
8. **測 responsive/accessibility** → 驗證：1440×1000、1024×768、390×844、320×568、200% zoom 均無 body overflow/sticky overlap；keyboard flow、focus、live announcement 與 chart table fallback 可用。
9. **測 authority** → 驗證：任何 module FAIL/BLOCKED 時全頁不得顯示「validated/pass」總結；sample-level `2 FAIL / 0 SKIP` 不被誤算為 cohort count。
10. **保存 evidence** → 驗證：每個 acceptance state 有 screenshot、browser error array、selected IDs、fixture/source、viewport、DOM assertion 與 machine-readable receipt。

### 10.2 Data assertions

| Assertion | Expected |
|---|---|
| CSV dataset rows | 7 |
| unique biological IDs | 6 |
| technical replicate | exactly 1：HCC1395_DORADO → HCC1395 |
| hash match / sample identity / receipt all pass | 7/7 / 7/7 / 7/7 |
| dataset-level all mutation-bearing families complete | 0/7 |
| dataset-level all units objective-certified | 0/7 |
| H2009 comparison values | tree coverage 64.16958%；unique among tree 35.286233% |
| HCC1954 comparison values | tree coverage 72.902143%；unique among tree 90.579069% |
| All-samples ISM/lineage aggregate | must not exist until equivalent manifests pass |

### 10.3 Interaction fixtures 與預期結果

| Fixture/state | 必要 screenshot | 必要 machine assertion |
|---|---|---|
| All-samples default | desktop + mobile | scope counts、gate counts、replicate grouping、no body overflow |
| Primary biological dataset | HCC1395 primary | all visible modules carry HCC1395 dataset ID |
| Technical replicate | HCC1395_DORADO | `technical_replicate=true`、`biological_id=HCC1395`、macro excluded |
| Capability PARTIAL | module card + reason | loaded/expected present；aggregate exclusion true |
| Capability ABSENT | module + selector context | no zero value；no borrowed source；stale DOM count=0 |
| Untested/invalid metric | axis table | explicit status/reason；not counted significant |
| Missing-data locus | Figure-17-equivalent state | `無 ISM` text；no broken asset placeholder |
| Evidence-present locus | Figure-18-equivalent state | reads/CpG + axis results + circularity warning |
| Selfcheck BLOCKED | desktop + mobile | fail/skip counts match selected dataset receipt |
| Network/parse/identity error | error + recovery | blocked module, error captured, no prior-sample values |

### 10.4 Browser pass criteria

- `browserErrors=[]`、page errors 0、unexpected failed requests 0；預期故障 fixture 必須被 receipt 明列，不混入 success run。
- 1440/1024/390/320 的 `documentScrollWidth == clientWidth`；表格 overflow 只能發生在命名清楚的 local container。
- 390px sticky layers 不 overlap，總高度 ≤96px；sample/dataset 與最高嚴重度 status 全程可見。
- 所有互動 control 有 accessible name、keyboard path、visible focus；所有狀態有文字而非只靠色彩。
- Screenshot QA 必須同時包含 overview、technical replicate、PARTIAL、ABSENT、error、missing locus、evidence-present locus與 mobile；只有 happy path 不算通過。
- 若 asset policy 仍是 `igv=skip, panels=0`，receipt 必須寫 `full image bundle NOT EVALUATED`；不可把 light QA 升格成 full-asset PASS。

### 10.5 2026-08-13 bounded snapshot 實際驗收

| Spec 項目 | 實際狀態 | 證據／界線 |
|---|---|---|
| 7 datasets／6 biological／1 technical | **PASS** | Hero 與四層 evidence strip 明示；DORADO 選取時 macro `n=0`／`n/a` |
| hash／identity／receipt 三 gate | **PASS** | 三個獨立 badge 均 7/7；不可解讀為 truth validation |
| macro median + IQR | **PASS** | Tree 77.680% [73.997, 78.600]；unique/tree 62.386% [46.662, 84.173]；n=6 |
| 七個 dataset 狀態切換 | **PASS** | All + 7 dataset 共 8 states；數值逐列對 artifact；非 HCC 維持 `ABSENT_NO_EQUIVALENT_BUNDLE` |
| PARTIAL authority | **PASS** | Sticky top bar 有文字 `PARTIAL`；已知資料缺口使用「尚未產生」而非 query failure 文案 |
| 四層閱讀 | **PASS** | Aggregate／HCC1395 canonical／H2009 extreme observed／HCC1395_DORADO technical control |
| Responsive／keyboard／fallback | **PASS** | 1440、1024、512（200% zoom 等效）、390、320 無 body overflow；表格只在 local container 捲動；鍵盤選取、焦點回復、aria-live、6 份 chart fallback 通過 |
| Bounded BAM/reference identity | **PASS_BOUNDED** | 7/7 tumor/normal quickcheck、fixed chunks、BAI/FAI/dictionary；full BAM SHA 0/7、14 BAM 無 RG，不升格 full/biological identity |
| Producer tag denominator | **PASS_EXISTING_RECEIPT** | 7/7 HP/all、HP+PS/all、HP+PS/HP、duplicate/all exact N/D；分母是 captured alignment records，不是 primary reads |
| Invalid axis 與 exact rails | **PASS** | cluster 排除於 valid chart、保留 `INVALID · DOUBLE-DIPPING`；active-k 8 rows、axis 5 rows exact N/D |
| 兩層 cascading selector | **SNAPSHOT DEVIATION** | Canonical portable reader 本版採一個 dataset selector；biological ID 與 replicate 關係仍在 hero、chart tooltip、table、DORADO state 明示。Production route 才實作依 biological sample cascade |
| 兩張對齊 topology plots | **SNAPSHOT DEVIATION** | 本版採同一 0–100% 軸的 grouped bar + semantic table；兩個 numerator/denominator 未混成 composite score。Production route 可改 aligned panels |
| Network/shard retry、race、URL state | **NOT APPLICABLE TO OFFLINE SNAPSHOT** | Standalone HTML 無外部 HTTP(S) request；不得把未執行的 async fixture 宣稱 PASS |
| Multi-BAM locus detail 正負狀態 | **NOT INTEGRATED** | Figure 17/18 只驗證既有 HCC1395 direct-generated browser；尚未接入多樣本頁面 |
| 每 BAM QC、truth、等價 ISM/lineage | **PARTIAL / NOT COLLECTED** | bounded identity 與 producer alignment-tag receipt 已有；depth/N50/phase blocks/MM/ML/CpG/KDE/truth 保留 null，不以 0 代替 |

機器收據：`InterSubMod/research/20260813_hcc1395_drilldown_validation/results/multi_bam_dashboard_browser_qa.json`。完整執行紀錄：`InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_multi_bam_dashboard_validation_01.md`。

## 11. Production boundary 與未完成項

本文件初稿只定義 spec；後續依使用者同一任務的明確續作要求，已另外完成 bounded artifact／standalone HTML／QA receipt，沒有修改原始 HCC1395 bundle。若要升級成可持續載入多 BAM 的 production route，仍需鎖定：

1. canonical sample manifest（biological ID、dataset ID、technical replicate、aliases、truth scope）；
2. cohort metric registry（numerator、denominator、unit、validity、aggregate policy）；
3. 等價 ISM/lineage manifest 是否存在；不存在即維持 PARTIAL/ABSENT；
4. full-asset policy、durable output path 與 async browser QA fixture budget；
5. All-samples macro 目前已明定為六個 primary biological datasets 的 unweighted median/P25/P75；若要改成 pooled denominator 或 hierarchical model，必須另開統計決策，不可靜默替換；
6. biological → dataset cascading selector、URL state、retry/race/error fixtures 與 live BAM processing contract。

## 12. Evidence mapping

| Claim/design decision | Source |
|---|---|
| 7 datasets / 6 biological samples / 1 technical replicate | `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/cohort_topology_metrics.csv` |
| topology 可比、ISM/lineage multi-sample 仍 PARTIAL | `InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_HCC1395_drilldown完整驗證與多樣本改進_01.md` |
| generated desktop/mobile hierarchy 與 86.19px sticky baseline | `InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/13_current_generated_desktop.png`、`14_current_generated_mobile.png` |
| co-occurrence denominator 與 selected-sample selfcheck | `InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/15_current_generated_cooccur_desktop.png`、`16_current_generated_selfcheck_mobile.png` |
| missing-data 與 evidence-present detail 兩種狀態 | `InterSubMod/research/20260813_hcc1395_drilldown_validation/figures/17_current_generated_detail_desktop.png`、`18_current_generated_methyl_detail_desktop.png` |
| bounded snapshot artifact 與 standalone HTML | `InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_artifact.json`、`20260813_multi_bam_analysis_overview.standalone.html` |
| Bounded BAM/reference 與 producer tag contract | `InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_input_manifest.json`、`results/multi_bam_input_manifest_validation.json` |
| 8 selector states、responsive、accessibility、exact rails 與 screenshot QA | `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/multi_bam_dashboard_browser_qa.json`、`figures/19_multi_bam_dashboard_all_desktop.png` 至 `29_multi_bam_dashboard_denominator_rails.png` |

---

**PARTIAL footer — Task Type B scope disclosure**：bounded offline snapshot 已覆蓋全 7 列 topology dataset，並通過資料產品層 canonical verification 與 40/40 browser assertions；尚未有 7 份等價 ISM/lineage bundles、完整每 BAM QC/truth、production cascading selector 或 async error fixtures。因此只可標 **snapshot data-product QA PASS / science PARTIAL**，不可標 scientific validated 或 production complete。
