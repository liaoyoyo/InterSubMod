<!--
建立時間: 2026-07-15 03:20
目標: Layered workstation GRCh38 分布與樣本全貌層實作過程 living document
處理範圍: layered-workstation-sample-overview
cycle_id: cycle_20260715-0320-layered-workstation-sample-overview
spec_id: layered-workstation-sample-overview
status: validated
advisory: on
關聯檔案:
  - InterSubMod/research/20260715_layered_workstation_sample_overview/pre-decision-audit.md
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_topology_workstation.py
-->

# Implementation Notes：Layered workstation 樣本全貌層

> **Purpose**: 記錄 GRCh38 coordinate view 與 sample-wide charts 的資料契約、偏離與 QA gate。
> **Task type**: B｜7 datasets × chr1–22 comprehensive validation。

## Frontmatter

- **Spec source**: 2026-07-15 使用者要求恢復舊頁的全基因分布與樣本重點觀察，但維持 current canonical 正確性。
- **AI session**: 2026-07-15。
- **Last updated**: 2026-07-15 04:22 +08:00。
- **Line count**: < 400。

## 🔵 設計決定（Design Decisions）

### 2026-07-15 03:20｜真正的 GRCh38 coordinate ideogram

- **Status**: Accepted。
- **背景**: current `chr1–22 全基因視圖` 是 22 格 count matrix，不是按座標落點的 genome distribution。
- **決定**: We will 由 current compact INDEX 的 `chrom/start/end` 與固定 GRCh38 autosome lengths 產生 22-row coordinate ideogram，並保留既有 chromosome count grid 作精確 lookup。
- **理由**: ideogram 回答「落在哪裡」；grid 回答「各 chromosome 有多少」，兩者任務不同。
- **影響範圍**: sample generator、7 generated HTML、Playwright assertions。
- **Revisit if**: scope 擴到 chrX/Y 或 reference build 改變。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐。

### 2026-07-15 03:20｜所有 overview 指標改綁 canonical v5

- **Status**: Accepted。
- **背景**: 舊 7,143 / 6,288 / 35,332 是 retired `is_somatic` universe。
- **決定**: We will 使用 current `W_tree/W_primary/complete`、`C_bins`、`hp_h3`、primary evidence 與 compact INDEX 重算；舊 count 一律不進 current page。
- **理由**: 避免 denominator shifting 與 archive snapshot 被誤認為現況。
- **影響範圍**: sample overview cards、圖例、說明文字、數值測試。
- **Revisit if**: canonical machine summary schema 變更。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐。

<!-- BEGIN USER-SPECIFIED -->
**Decision**: 保留 chr1–22 全基因視圖、加入清楚的樣本全部重點觀察、原始 JSON 維持隱藏連結。
**DO NOT change**: 不可用視覺精簡為理由移除全基因或把 JSON 放回主要閱讀流。
**Rationale**: 使用者 2026-07-15 明確要求。
<!-- END USER-SPECIFIED -->

### 2026-07-15 03:20｜Chart contract

- **Status**: Accepted（7/7 data audit + 4 viewport Chromium audit 通過）。
- **Analytical question**: regions 在 GRCh38 哪裡，以及該樣本的 candidate space、可辨識度、HP 組成與 region complexity 如何分布？
- **Takeaway ceiling**: 描述 current canonical reconstruction landscape；不推論 hotspot、enrichment、biological ancestry 或 clone truth。
- **Charts**:
  - coordinate distribution：22-row GRCh38 ideogram，region midpoint marks；mode=determinacy/evidence/primary HP/n_sSNV/CN sidecar。
  - composition：100% stacked bar + exact values，denominator 固定顯示。
  - count distributions：horizontal zero-baseline bars，依 count 排列或具語意順序。
- **Data sufficiency**: 每 dataset 3,612–9,674 W_tree rows；所有 marks 同一 region grain。
- **Palette**: 延用 current green/blue/magenta/slate；同時以文字、pattern 與排序提供非色彩辨識。
- **Final surface**: offline generated sample HTML；桌機與手機 Chromium。
- **Evidence tier**: L2 ⭐⭐⭐⭐。

### 2026-07-15 04:22｜舊 morphology vocabulary 只保留 boundary explanation

- **Status**: Accepted。
- **背景**: 舊 `single/linear/branched/star` 與 current `Topo_region` 並非相同 grain；stored current display tree 可能只保留前 32 棵，不能事後恢復 exhaustive morphology family。
- **決定**: 不報 current morphology count；在 closed explanation drawer 解釋為何改用 Topo / C / determinacy，並保留日後 producer versioned summary 的恢復條件。
- **理由**: 熟悉名詞可保留教學價值，但不可把不同 universe、candidate grain 與舊分類錯接成 current canonical。
- **恢復 gate**: producer 必須輸出 exhaustive、versioned、可閉合 W_primary 的 morphology-family summary，且明確處理 multi-primary-HP region。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐。

## 🟠 偏離之處（Deviations）

### 2026-07-15 03:20｜不直接復刻舊 single/linear/branched/star region counts

- **Status**: Accepted；依 7/7 audit 確認不可安全直接映射。
- **規範要求**: 使用者希望保留舊頁可一眼理解的 topology type 觀察。
- **實作偏離**: current page 不會直接搬入舊 region-level labels；先審核是否能在 current primary-unit grain 安全重算。
- **理由**: current `Topo_region` 是跨 primary HP units 的 shape-count product，與舊 cluster-first single tree ontology 不同。
- **風險評估**: 若只移除名稱，會降低熟悉度；若硬映射，會造成更嚴重科學誤導。
- **回退路徑**: 顯示 retired vocabulary 對照；不以 subset 取代 task type B 的 current sample-wide count。
- **Revisit if**: producer 提供 exhaustive versioned summary，且 7/7 candidate tree audit 證明分類閉合與 grain 清楚。
- **Evidence tier**: L2 ⭐⭐⭐⭐。

## 🟡 折衷考量（Trade-offs）

### 2026-07-15 03:20｜Ideogram 與 chromosome grid 並存

- **Status**: Accepted。
- **方案 A**: We will 保留 ideogram + count grid；前者看位置、後者查精確分母。
- **方案 B**: 只留 ideogram；精確 chromosome counts 難查。
- **方案 C**: 只留 grid；就是 current 缺失，無法回答 genomic distribution。
- **採用 A 理由**: 支援 overview → chromosome → region 的完整下鑽，不讓單一圖承擔兩種任務。
- **Revisit if**: mobile 頁長或效能超過驗收門檻。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐。

## 🔴 未決問題（Open Questions）

### 2026-07-15 03:20｜Current morphology 的安全 grain

- **Status**: resolved as non-reportable from current stored payload。
- **Question**: familiar `single/linear/branched/star` 是否只在 candidate-complete、Topo_unit=1 的 primary unit 可安全分類？
- **Context**: current region 可能同時有 HP1/HP2 兩個獨立 candidates，不能假裝為一棵合併樹。
- **Resolution**: 不報 morphology count，只顯示 current C/Topo/determinacy 與 retired vocabulary boundary；stored display tree 不足以支撐 exhaustive 事後分類。
- **Revisit if**: producer 新增 versioned exhaustive summary；不是靠 UI 猜測或舊 count 回填。
- **Priority**: critical。
- **Evidence tier**: L5 ⭐。

## 📚 Lore（Prior Gotchas / Non-obvious Constraints）

### 2026-07-15｜舊頁能畫圖是因為另接 ideogram_data

- **Constraint**: 舊 builder 讀 `ideogram_data.json`；HCC1395 還有專用 fallback 路徑。current generator 沒有這段接線。
- **Why it matters**: 缺圖不是 CSS 或 Chromium 問題。
- **Evidence**: `build_topology_workstation.py:62-86,279-353`。

### 2026-07-15｜舊 mobile baseline 本身有 overflow

- **Constraint**: 390px viewport 的 legacy document width 為 506px。
- **Why it matters**: 重用資訊結構時不可照抄固定寬度 SVG/卡片。
- **Evidence**: `/tmp/20260715_hg38_metric_baseline/metrics.json`。

### 2026-07-15｜Claude Code 與 Chromium 終驗

- **Pre-implementation Claude verdict**: `PASS_WITH_CHANGES`；三個 P0（分母固定顯示、n_sSNV 2–8 fail-closed、C null/closure）均已落地。
- **Post-implementation Claude verdict**: `PASS`、P0=0、P1=0、`COMMIT READY: YES`。
- **Post-review P2 也已處理**: 常駐 smoke 現在逐一驗證 5 種 mode；一次性 audit 直接比對 generated renderer SHA 與當前 renderer 檔案 SHA，而非只看 64 字元格式。
- **Chromium evidence**: 7 datasets × 4 viewports；32/32 page-runs PASS；0 console/page/request/external-network error；0 body overflow。
- **Aggregate（每 dataset 僅計一次）**: W_tree=51,815；W_primary=50,215；retained sSNV=194,149。
- **Evidence**: `InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_sample_overview_after/metrics.json`。

## Provenance Footer

- **Base commit**: `2bec873`
- **Build time**: 2026-07-15 04:22 +08:00
- **Skill version**: `/implementation-notes` v0.1
- **Canonical summary SHA-256**: `71da78b69f8afe5fb8e618179ab7b38a6940fdb17be6282f99a6ec4b720e5de7`
