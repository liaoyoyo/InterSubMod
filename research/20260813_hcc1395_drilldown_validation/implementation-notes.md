<!--
建立時間: 2026-08-13 01:26
目標: 記錄 HCC1395 drilldown 全面稽核與 generator hardening 的設計決定、偏離、折衷與未決事項
處理範圍: HCC1395_v1 immutable audit、sample fail-closed、receipt/UI/test/report
cycle_id: cycle_20260813-0126-hcc1395-drilldown-validation
spec_id: 20260813_hcc1395_drilldown_validation
status: in_progress
advisory: on
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/pre-decision-audit.md
  - InterSubMod/scripts/build_drilldown_dashboard.py
  - InterSubMod/tests/test_drilldown_contract.py
-->

# Implementation Notes：HCC1395 drilldown validation

## Frontmatter

- **Spec source**: 使用者 2026-08-13 對話需求 + pre-decision audit
- **AI session**: 2026-08-13
- **Last updated**: 2026-08-13 04:05
- **Status**: in_progress

## 🔵 設計決定（Design Decisions）

### [2026-08-13 01:26] 保留 v1 為 immutable evidence

- **Status**: Accepted
- **背景**: v1 含可重現的 stale claim、1 FAIL/1 SKIP 與 receipt 缺口，原地重建會消除失敗證據。
- **決定**: We will not overwrite `/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/`; fixes land in source/tests and a separate audit report.
- **影響範圍**: 原始 1.7 GB bundle 唯讀；任何新 fixture 或重建落新路徑。
- **Revisit if**: 使用者明確要求封存/重建並指定新輸出目錄。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 01:26] claim ceiling 固定為 observation-only

- **Status**: Accepted
- **背景**: bundle 沒有 SEQC2 truth VCF、HighConf BED、som.py/hap.py 或 TP/FP/FN。
- **決定**: We will treat topology/methylation/lineage summaries as descriptive observations, not caller F1 validation or causal evidence.
- **影響範圍**: 報告、receipt status、UI 文案。
- **Revisit if**: truth、callset、benchmark region 與 metric contract 全部納入同一 provenance chain。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 01:26] 多樣本來源必須 fail-closed

- **Status**: Accepted
- **背景**: `--sample COLO829 --probe-only` 可讀入 HCC1395-only LCA/lineage/ISM roots。
- **決定**: We will validate sample identity before accepting a capability and add a low-linkage/mismatch regression fixture.
- **影響範圍**: source loaders、site config、CLI probe、tests。
- **Revisit if**: 建立 canonical sample alias registry；屆時以 manifest ID 取代 path-token fallback。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 03:10] 舊 standalone 只移植資訊架構，不移植分類結果

- **Status**: Accepted
- **背景**: 舊頁的 `A_ALLELE_clean` 與目前 dashboard 的 raw ALT-axis 並非同一 contract；30,077 舊 loci 與 19,849 current loci 只有 16,566 個座標交集。
- **決定**: We will reuse the selection funnel, progressive disclosure, compact curated-case filter and evidence co-location. We will not map legacy A/B classes, V thresholds, 472 candidates or 14 displayed cases into current ISM metrics.
- **影響範圍**: 新 observation scope、axis definition table、missing/untested third state、視覺比較與 report wording。
- **Revisit if**: 兩條 pipeline 共享 frozen input manifest、同一 locus universe、同一 metric definition 與 multiplicity family，且 crosswalk 驗證通過。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 03:10] 新視覺驗收必須使用 source-complete staging

- **Status**: Accepted
- **背景**: CSS/JS route overlay 能驗前端資產，但無法替換舊 `index.html` shell；第一次 machine receipt 仍看到舊 placeholder，與文字摘要不一致。
- **決定**: We will treat the overlay as a partial fixture only and rerun final screenshots against an index generated from the updated source before claiming the new information architecture is complete.
- **影響範圍**: visual QA receipt、screenshots、final report；legacy bundle 仍 immutable。
- **Revisit if**: source-complete fixture 已生成並通過 desktop/mobile/browser checks。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 04:05] 將 observation scope 與 detail availability 拆成兩條證據鏈

- **Status**: Accepted
- **背景**: legacy funnel 混合候選母體、統計可檢定性與案例展示；current 又可能把「有 detail」誤讀為「有顯著證據」。
- **決定**: We will persist `topology → ISM summary → non-circular testable → any global raw-p candidate` as the observation ledger, and expose locus-detail availability separately.
- **影響範圍**: payload、analysis UI、axis code 0/8/16、source-complete screenshots、report metric contract。
- **Revisit if**: cohort-wide multiplicity family 與 effect contract 完成；屆時另加 q-value decision layer，不覆寫 raw-p exploratory layer。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 04:05] Crosswalk 必須標示 immutable 舊 axis encoding

- **Status**: Accepted
- **背景**: current crosswalk 讀 HCC1395_v1 immutable index；該頁只含 axis code 0–8，而 source hardening 後才新增 `AXIS_UNTESTED=16`。
- **決定**: We will record the source encoding provenance in `legacy_browser_crosswalk.json` and prohibit interpreting absence of code 16 as evidence that every summary row was tested.
- **影響範圍**: crosswalk JSON、方法附錄、主報告。
- **Revisit if**: 以 clean source 重建新 bundle 後另做 versioned crosswalk；不得修改 legacy result。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 05:05] 多 BAM 總覽先做 canonical snapshot dashboard，不直接改 production route

- **Status**: Accepted
- **背景**: 使用者需要一頁式可選樣本總覽，但現有可比資料只有 7 份 topology；等價 ISM/lineage bundles 尚不存在。
- **決定**: We will build one canonical `surface=dashboard` artifact and package it as self-contained HTML. The page will use 7 complete topology datasets, biological macro summaries, a global sample selector, HCC-only extension panels, and explicit PARTIAL/ABSENT states.
- **影響範圍**: 新 metric/design specs、artifact builder、portable HTML 與 QA receipt；不修改 immutable bundles 或 production generator route。
- **Revisit if**: sample manifest、refresh ownership、全樣本 ISM/lineage 與 production output contract 凍結後，再把同一 metric model接入正式軟體。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

## 🟠 偏離之處（Deviations）

### [2026-08-13 01:26] 本輪不生成七資料集完整 dashboard

- **Status**: Accepted
- **規範要求**: Task Type B 的跨樣本結論需全樣本、全 scope。
- **實作偏離**: 只盤點可用 upstream 資料並設計 cohort contract；不宣稱已完成多樣本驗證。
- **理由**: `drilldown_out/` 只有同一 HCC1395 的 v1/v3，且 v3 也是 dirty/provenance 不完整；計算範圍需另行鎖定。
- **風險評估**: 不影響單 bundle 工程稽核；限制跨樣本結論。
- **回退路徑**: 新 cycle 以 sample-specific manifest 生成全七資料集並做 truth-aware cohort QA。
- **Revisit if**: 七份 inputs 與預算/輸出位置獲確認。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 03:10] Portable HTML report 未通過 shared renderer QA

- **Status**: Accepted
- **規範要求**: T3 複雜報告原預期同時交付 Markdown 與 portable HTML。
- **實作偏離**: Markdown、artifact JSON 與全部圖檔可交付；portable HTML 不列為通過成品。
- **理由**: shared renderer 的 `.analytics-top-bar` 使用 `100vw`，長頁出現 15px scrollbar 後造成 desktop `clientWidth=1425 / scrollWidth=1433`、mobile `375 / 383` 的水平溢位；不是 report data/content 問題。
- **風險評估**: 不影響 metric/report conclusions；若把失敗 HTML 交付會違反 visual QA contract。
- **回退路徑**: 以完整 Markdown 為 primary report，保留 artifact/debug evidence；等待 renderer 修正後再封裝。
- **Revisit if**: shared renderer 將 full-bleed top bar 改為不超出 client width，canonical deliver browser verification 通過。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 04:05] 完整圖像資產驗證延後，不把 light QA 冒充 full build

- **Status**: Accepted
- **規範要求**: 完整 bundle 驗收應同時檢查 panels、IGV、inventory 與 browser runtime。
- **實作偏離**: source-complete staging 讀完整 v3 summary，但採 `igv=skip, panels=0`；本輪只驗 IA、payload、authority、responsive runtime。
- **理由**: v3 IGV JS 約 4.34 GB；現有 upstream provenance/科學 gate 尚未通過，先重建全部靜態資產無法提升 claim tier。
- **風險評估**: 不影響 source/runtime contract；完整 asset delivery 仍標 NOT EVALUATED。
- **回退路徑**: clean commit + frozen manifest 後，在新 output 路徑執行 full asset build、hash inventory 與同一 browser suite。
- **Revisit if**: truth/provenance P0 關閉且使用者核定長計算輸出與容量預算。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

## 🟡 折衷考量（Trade-offs）

### [2026-08-13 01:26] 小 fixture 驗收優先於大型重建

- **Status**: Accepted
- **方案 A**: We will patch fail-closed/provenance/UI contracts and validate with unit/browser fixtures.
- **方案 B**: 立即重建 v3/v4；rejected，因來源本身含 dirty build 且會耗用數 GB，不能證明多樣本隔離。
- **方案 C**: 只寫報告不修 source；rejected，因 COLO829 contamination 是實際 defect。
- **採用 A 理由**: 可直接否證機制、可逆且不破壞原始 evidence。
- **Revisit if**: fixtures 全通過且有 clean source manifest，可啟動新 bundle build。
- **Evidence tier**: L2 ⭐⭐⭐⭐

## 🔴 未決問題（Open Questions）

### [2026-08-13 01:26] 七資料集 cohort 的正式輸入與真值範圍

- **Status**: open
- **Question**: 下一輪是否要生成完整七資料集 drilldown，及哪些樣本具 truth/HighConf 可進 validation subset？
- **Context**: upstream topology 有七資料集，但 dashboard bundle 只有 HCC1395；technical replicate 不可當獨立 biological n。
- **Default if no answer**: 本輪只交付 inventory、contract 與 no-claim 結論，不啟動長計算。
- **Revisit if**: 使用者確認 scope、輸出位置、truth availability 與計算預算。
- **Priority**: major
- **Evidence tier**: L5 ⭐

## 📚 Lore（Prior Gotchas / Non-obvious Constraints）

### [2026-08-13] Capability presence 不等於 validation

- **Constraint**: v1 同時顯示 `6/6 完整` 與 C10 FAIL/C8 SKIP。
- **Why it matters**: source 可讀、linkage coverage、自檢與科學有效性必須分開呈現。
- **Evidence**: v1 `receipt.json`、`SELFCHECK.md`、browser screenshot。

### [2026-08-13] v1 receipt 不能重建 ISM source

- **Constraint**: receipt 對 ISM directory 記為 size 0、無 hash；原路徑已不存在。
- **Why it matters**: 即使 topology/MLHP hashes 尚匹配，完整 bundle 仍不可重建。
- **Evidence**: v1 `receipt.json` + path existence/hash audit。

### [2026-08-13] panel bytes 只統計 base panel

- **Constraint**: v1 receipt 記 139,037,816 bytes，但全 PNG 實測 221,554,648 bytes；tumor-only panel 未加總。
- **Why it matters**: delivery size、storage budget 與 integrity summary 被低估約 82.5 MB。
- **Evidence**: `InterSubMod/scripts/drilldown/panels/bake.py` 與全量 PNG scan。

## Provenance Footer

- **Commit hash at start**: `73afaeac8e61c767241fa59c1ca6043a1c95290c`
- **Skill version**: `/implementation-notes` v0.1
- **Active cycle**: `InterSubMod/state/cycles/cycle_20260813-0126-hcc1395-drilldown-validation/`
- **Pre-decision audit**: `InterSubMod/research/20260813_hcc1395_drilldown_validation/pre-decision-audit.md`
