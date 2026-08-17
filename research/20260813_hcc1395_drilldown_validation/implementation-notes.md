<!--
建立時間: 2026-08-13 01:26
目標: 記錄 HCC1395 drilldown 全面稽核與 generator hardening 的設計決定、偏離、折衷與未決事項
處理範圍: HCC1395_v1 immutable audit、sample fail-closed、receipt/UI/test/report
cycle_id: cycle_20260813-0126-hcc1395-drilldown-validation
spec_id: 20260813_hcc1395_drilldown_validation
status: completed_bounded_snapshot
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
- **Last updated**: 2026-08-13 12:12
- **Status**: completed_bounded_snapshot（science PARTIAL；production NOT COMPLETE）

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

### [2026-08-13 session-continuation] Inventory zero 與 biological null 分開

- **Status**: Accepted
- **背景**: 非 HCC 樣本沒有 audited downstream bundle；若 hero card 顯示 n/a，使用者無法區分「尚未盤點」與「盤點後為 0」。但把 ISM/lineage 生物值寫 0 又會製造假陰性。
- **決定**: We will use `downstream_available_count=0` only for the reviewed bundle inventory, while sample-matched ISM/lineage numerators, rates and coverage stay null with `ABSENT_NO_EQUIVALENT_BUNDLE`.
- **影響範圍**: COLO829 等非 HCC sample card 顯示 0/1；HCC-only chart/table 為 explicit empty state，不借用資料。
- **Revisit if**: 新增 sample-matched bundle；屆時 count 與 capability status 由 manifest 自動更新。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 session-continuation] 主次資訊採 native dashboard + no-script evidence blocks

- **Status**: Accepted
- **背景**: 使用者要求主指標優先、細節摺疊；canonical artifact schema 沒有能包住 native filtered tables 的 collapse block。
- **決定**: We will keep filters/cards/charts/tables native and source-backed, then add three sandboxed no-script HTML blocks：four-layer strip、7-dataset fixed-label opportunity panel、以及含四個 closed `<details>` 的 scope/formula/missing/identity 區。
- **影響範圍**: 原生 sample filter 不受影響；fold 明示 selected-sample 值以上方動態 cards 為準。
- **Revisit if**: shared artifact schema 新增 native evidence strip/collapsible group；屆時移除 custom HTML blocks。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 12:05] v1.1 ingestion 契約固定 schema、receipt 與動態 identity

- **Status**: Accepted
- **背景**: 獨立對抗審查找出 permissive source schema、receipt summary 可篡改、dashboard 寫死 0/7 drift 與 implicit-input overwrite 缺口。
- **決定**: We pin canonical source schema `$id`+SHA, revalidate complete receipt/manifest coherence in the dashboard builder, dynamically aggregate all-role MATCH/drift, and protect every explicit/implicit input from atomic output replacement.
- **驗證**: 32/32 regressions PASS；all-MATCH=`7/0`，mixed All=`1/6`、HCC1395=`1/0`、H1437=`0/1`；permissive schema、receipt tampering、collision 均 fail-closed。
- **影響範圍**: input builder/schema/receipt、dashboard builder、QA runner、`.gitignore` trackability。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 12:08] LongPhase-S tag 分母與 duplicate guard 分開命名

- **Status**: Accepted
- **背景**: producer sidecar 統計的粒度是 captured mapped alignment records，不是 primary/unique reads；HCC1395 與 HCC1937 duplicate identity rates 偏高。
- **決定**: We expose HP/all, HP+PS/all, HP+PS/HP and duplicate/all with exact N/D; duplicate is a denominator-composition guard and never a biological/accuracy score.
- **驗證**: 7 datasets 四組 rates 與 source counts、42 chart rows、14 detail rows均一致。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 12:10] 有效甲基化軸與 circular cluster 分離

- **Status**: Accepted
- **決定**: We exclude cluster from the valid-axis chart, retain its exact N/D only in the validity rail, and label it `INVALID · DOUBLE-DIPPING`; active-k also receives an 8-row exact N/D rail.
- **驗證**: browser 40/40；cluster chart rows=0、validity rail=5 axes，N/D 與 source 一致。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

## 🟠 偏離之處（Deviations）

### [2026-08-13 session-continuation] Dashboard 封裝使用 bounded renderer compatibility patch

- **Status**: Accepted
- **規範要求**: portable dashboard 應由 canonical artifact renderer 產生並通過 desktop/mobile verifier。
- **實作偏離**: 仍使用 canonical builder/extractor/verifier，但在輸出加入 bounded compatibility layer：以 content width + gutter 取代 top-bar `100vw`、≤600px 限制 legend flex width，並補上 persistent PARTIAL、filter focus restoration/aria-live 與 known-gap 正確文案。
- **理由**: 共用 renderer 在 classic scrollbar iframe 將 `100vw` 算入 8–15px scrollbar gutter；甲基軸 legend 在 390px 又以 max-content 撐到 416.5px。
- **風險評估**: 不改 payload、資料、圖表 encodings 或 table scroll；JS 只增強既有 filter 的焦點/announcement 與靜態 status/copy，未加入資料載入 runtime，亦未用 global overflow clipping 隱藏錯誤。
- **驗證**: canonical delivery 35 rendered blocks/10 metrics/6 charts/10 tables/3 HTML blocks，1440/390/source dialog PASS；Playwright 40/40 assertions、0 console/page/external-request errors；1024/512/390/320 body width 全相等。
- **回退路徑**: shared renderer upstream 修正後，移除 wrapper CSS 並重跑同一 verifier/QA。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 session-continuation] Claude Code Round 5 因 session quota 未執行

- **Status**: Accepted as execution blocker
- **規範要求**: 最終 dashboard 與 Claude Code 再做獨立唯讀交叉審查。
- **實作偏離**: 預設 high 與 `--model haiku` 兩次都在任何 Read/Bash 前回傳 session limit、exit 1。
- **理由**: 最後一次 Claude Code 2.1.229 high-effort 嘗試在任何 Read/Bash 前回傳 session limit，當時顯示 14:20 Asia/Taipei 重置；不是 repo/test failure。
- **風險評估**: Round 4 只涵蓋 underlying source，不可冒充 final dashboard review。
- **回退路徑**: 先保存 `InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round5_review.md` execution receipt，使用 canonical validator、interaction QA 與獨立 Codex review；quota 重置後補跑原 prompt。
- **Evidence tier**: L1 execution receipt；review verdict absent

### [2026-08-13 session-continuation] 獨立 dashboard re-audit 限定接受 snapshot 層

- **Status**: Accepted
- **背景**: 第一輪獨立 reviewer 不接受 full design-spec product claim，指出 QA breadth、首屏 hierarchy、replicate/IQR、four-layer 與 accessibility 缺口。
- **決定**: We will close measurable snapshot blockers, rerun canonical/browser QA, and state the verdict only as `snapshot data-product QA PASS / science PARTIAL / production NOT COMPLETE`.
- **驗證**: 最終三路 re-audit ACCEPT/PASS；ingestion 32/32，source/value 682/682，browser 40/40；8 selector states、non-borrowing、1024/512/390/320、hash/payload/builder reconciliation、focus/aria-live 全一致。
- **未提升項**: 單層 selector、grouped chart 是明示 snapshot deviation；async/URL/retry/locus、per-BAM QC/truth/KDE、多樣本 ISM/lineage 仍未完成。
- **Evidence tier**: L1 independent reviewer receipt

### [2026-08-13 01:26；05:55 更新] 本輪不生成七資料集等價 ISM/lineage bundle

- **Status**: Accepted
- **規範要求**: Task Type B 的跨樣本結論需全樣本、全 scope。
- **實作偏離**: 後續已完成 7-dataset topology bounded snapshot dashboard，但沒有生成七份等價 ISM/lineage bundles，也不宣稱已完成多樣本科學驗證。
- **理由**: `drilldown_out/` 只有同一 HCC1395 的 v1/v3，且 v3 也是 dirty/provenance 不完整；計算範圍需另行鎖定。
- **風險評估**: 不影響單 bundle 工程稽核；限制跨樣本結論。
- **回退路徑**: 新 cycle 以 sample-specific manifest 生成全七資料集的等價 downstream bundle，補 per-BAM QC/truth 後做 truth-aware cohort QA。
- **Revisit if**: 七份 inputs 與預算/輸出位置獲確認。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 03:10] 長篇 Portable HTML report 未通過 shared renderer QA

- **Status**: Accepted
- **規範要求**: T3 複雜報告原預期同時交付 Markdown 與 portable HTML。
- **實作偏離**: Markdown、artifact JSON 與全部圖檔可交付；portable HTML 不列為通過成品。
- **理由**: shared renderer 的 `.analytics-top-bar` 使用 `100vw`，長頁出現 15px scrollbar 後造成 desktop `clientWidth=1425 / scrollWidth=1433`、mobile `375 / 383` 的水平溢位；不是 report data/content 問題。
- **風險評估**: 不影響 metric/report conclusions；若把失敗 HTML 交付會違反 visual QA contract。
- **回退路徑**: 以完整 Markdown 為 primary report，保留 artifact/debug evidence；等待 renderer 修正後再封裝。
- **Revisit if**: shared renderer 將 full-bleed top bar 改為不超出 client width，canonical deliver browser verification 通過。此項只指長篇 report；後續 multi-BAM dashboard 以獨立 artifact + bounded compatibility patch 通過，不回溯宣稱原 report HTML 已通過。
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
- **Evidence**: `/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/receipt.json`、`/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/SELFCHECK.md`、browser screenshot。

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
