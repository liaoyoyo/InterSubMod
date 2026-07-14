<!--
建立時間: 2026-07-15 Asia/Taipei
目標: 記錄 layered workstation 三維導讀實作的設計決定、偏離、折衷、未決與驗證
處理範圍: index + 7 canonical v5 sample pages；desktop/mobile；generator-first
cycle_id: 20260715-layered-workstation-concept-dimensions
spec_id: layered-workstation-concept-dimensions
status: validated
advisory: on
build_branch: research/subclonal-reconstruction-202606
build_commit: 6067568
worktree: /big7_disk/liaoyoyo2001/InterSubMod
data_sources: research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
關聯檔案:
  - InterSubMod/research/20260715_layered_workstation_concept_dimensions/pre-decision-audit.md
  - InterSubMod/research/20260714_layered_workstation_redesign/implementation-notes.md
  - InterSubMod/docs/methodology/_assets/topology_workstation/HCC1395.html
-->

# Implementation Notes: Layered Workstation Concept Dimensions

> **Purpose**: 把舊頁有用的概念教學形式，重新綁定到 canonical v5 的三個正確閱讀維度：拓樸型態、可辨識度、基因體位置。

## Frontmatter

- **Spec source**: 使用者 2026-07-15 指示 + current canonical v5 contract
- **AI session**: 2026-07-15
- **Last updated**: 2026-07-15
- **Line count**: <400

## 🔵 設計決定（Design Decisions）

### D-001 — 三維導讀綁定 current candidate-set contract

- **Status**: Accepted
- **背景**: Current page 有 C/Topo 與座標元件，但缺少「這維度在問什麼、如何讀、不能推論什麼」的集中入口。
- **決定**: We will add three explicit dimensions—拓樸型態、可辨識度、基因體位置—at overview and selected-region levels, with current data, reading rule, boundary, and direct interaction.
- **理由**: 讓讀者先建立心智模型，再讀數字與網路；不需回到 archived page 才懂名詞。
- **影響範圍**: v5 renderer、cohort index、7 generated pages、Chromium audit。
- **Revisit if**: 5 秒測試仍無法回答三個問題。
- **Evidence tier**: L1 user specification / L3 design choice

<!-- BEGIN USER-SPECIFIED -->
### D-002 — 舊頁僅作教學形式參考

- **Status**: Accepted
- **Decision**: 舊 `topology_workstation/HCC1395.html` 的概念標題與分面形式可參考；舊 `single/linear/branched/star`、A/B/C determinacy、clone-first 數字與分母不可移植。
- **DO NOT change**: current canonical v5 是唯一資料與語意主體。
- **Rationale**: 舊工作站已 deprecate；新頁回答的是 regional mutation-state candidate set，不是 confirmed biological clone/ancestry。
- **Evidence tier**: L1
<!-- END USER-SPECIFIED -->

<!-- BEGIN USER-SPECIFIED -->
### D-003 — 保留全基因視圖與隱藏 JSON 來源

- **Status**: Accepted
- **Decision**: We will strengthen, not replace, the chr1-22 view; raw `.json` links remain inside the default-collapsed evidence drawer.
- **DO NOT change**: 使用者前輪已明確指定兩項。
- **Evidence tier**: L1
<!-- END USER-SPECIFIED -->

### D-004 — 位置統計與落地列表分母必須同畫面明示

- **Status**: Accepted
- **背景**: 位置卡以 chromosome 的 W_primary 比較 count/rate，但點入 chromosome browser 時合理地保留全部 W_tree（含 no-primary）。
- **決定**: 卡片、CTA 與 genome status 同時寫出 W_primary 統計與 W_tree 落地筆數；CTA 不暗中排除 no-primary。
- **理由**: 保留全基因資料完整性，同時避免讀者把兩個 grain 混為同一分母。
- **Evidence tier**: L3 design / L5 full-data verification

## 🟠 偏離之處（Deviations）

### DV-001 — 不使用舊頁 emoji 作視覺資產

- **Status**: Accepted
- **規範要求**: 使用者以 🌳／🎯／📍 指認三個概念。
- **實作偏離**: We will use textual dimension labels and the current workstation design system, not emoji as icons.
- **理由**: emoji 在不同平台外觀不一致；概念名稱與數字才是必要訊息。
- **風險評估**: 低；保留完全相同的中文概念名稱。
- **回退路徑**: 若實測掃描性不足，可加入既有 icon library；不手工偽造 SVG。
- **Evidence tier**: L3

## 🟡 折衷考量（Trade-offs）

### T-001 — 三維導讀與既有 C/Topo 圖例的重複風險

- **Status**: Accepted
- **方案 A**: We will make the three cards task-oriented（問什麼／目前觀察／怎麼下鑽／限制），而既有 strip 保留 composition 詳細數字。
- **方案 B**: 只改現有三張 C/Topo 定義卡；rejected，因缺位置維度與 region-level reading model。
- **方案 C**: 回復舊 facet board；rejected，會重帶 archived schema。
- **採用 A 理由**: 分工清楚，避免同一段重複列相同數字。
- **Revisit if**: screenshot audit 顯示 overview 過長或資訊重複。
- **Evidence tier**: L3

### T-002 — 不以 `overflow-x:hidden` 掩蓋手機問題

- **Status**: Accepted
- **方案 A**: production body 保持 visible，runner 直接要求實測水平 overflow=0。
- **方案 B**: 固定 hidden/clip；rejected，因即使 descendant 溢出也可能被 CSS 裁掉而不被發現。
- **結果**: 32-context 專項與 24-context full audit 均為 0px body overflow；full runner 契約改為 visible 僅在 measured overflow=0 時通過。
- **Evidence tier**: L5 browser observation

## 🔴 未決問題（Open Questions）

### Q-001 — 最穩定的位置觀察句

- **Status**: answered (2026-07-15)
- **Question**: sample-level 位置卡應顯示最多 W_primary 的 chromosome，或顯示 incomplete rate 最高的 chromosome？
- **Context**: 兩者都只能作 descriptive observation，不可暗示 enrichment。
- **Options**: A 同時顯示 count leader + rate leader；B 只顯示目前 active chromosome；C 只教 coordinate/overlap search。
- **Answer**: 採 A；同時顯示 W_primary 數量最高 chromosome 與 incomplete rate 最高 chromosome，兩者都帶 numerator/denominator，並明示「描述性巡覽，未校正 chr 長度、sSNV density 或 region construction」。
- **Priority**: major
- **Evidence tier**: L5

## 📚 Lore

### 2026-07-15 — Determinacy 必須拆 exact 與 shape 兩層

- **Constraint**: `C=1` 是 exact candidate 唯一；`Topo=1` 是不帶節點標籤的 topology shape 唯一；`C>1,Topo=1` 不能被簡化成「不唯一」。
- **Why it matters**: 這正是 current 模型比舊 A/B/C determinacy 更精確的資訊。
- **Evidence**: `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py:142`

### 2026-07-15 — 位置只作描述與導覽

- **Constraint**: `chrom/start/end` 支持位置查找與分布描述，不自動支持 hotspot、mappability、telomere、centromere 或 artifact causal claim。
- **Why it matters**: 避免把 UI 排序誤寫成生物學結論。
- **Evidence**: current payload fields + pre-decision red-team failure mode C。

## 🧪 Verification Log

### V-001 — 改版前 reference/current 同 viewport capture

- **Input**: archived/current HCC1395 HTML。
- **Command**: `python3 InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_before/capture_concept_dimensions_before.py`
- **Output**: `InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_before/`
- **Observed result**: Chromium 147；desktop 1440×1000 + mobile 390×844；4 runs / 8 screenshots；0 page/console error；0 body overflow。Mobile current 頁需先越過約 1,496px 的 chromosome section 才到 C/Topo 定義。

完成狀態：7-page rebuild、Claude Code review、32-context focused QA、24-context full QA 與 hash receipts 全部完成。

### V-002 — Claude Code Sonnet read-only review

- **Input**: archived HCC1395、v5 renderer、index driver、canonical summary。
- **Command**: `claude -p '<三維／C-Topo／mobile／WCAG review prompt>' --model sonnet --effort high --allowedTools 'Read,Grep,Glob'`
- **Output**: stdout review；未修改檔案。
- **Observed result**: P0=0、P1=0。P2 的 index 位置層次、topology strip a11y、muted contrast 已納入實作；P3 zero-count evidence option 已改為動態 count + disabled。

### V-003 — 7/7 hash-bound rebuild

- **Input**: `InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json` + 7 canonical region-view JSONs。
- **Command**: `python3 InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py`
- **Output**: `InterSubMod/docs/methodology/_assets/layered_workstation/index.html` + 7 sample HTMLs。
- **Observed result**: 7/7 `BUILT hash-bound page`；summary SHA `71da78b69f8afe5fb8e618179ab7b38a6940fdb17be6282f99a6ec4b720e5de7`；HCC1395 inline JS `node --check` rc=0。

### V-004 — Three-dimension focused Chromium audit

- **Input**: index + 7 current sample HTMLs。
- **Command**: `python3 InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_after/capture_concept_dimensions_after.py`
- **Output**: `InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_after/metrics.json` + 34 screenshots。
- **Observed result**: Chromium 147.0.7727.15；8 pages × 4 viewports = 32 contexts；542/542 checks passed；0 console/page/request failures；0 horizontal overflow。

### V-005 — Generic layered + archived regression

- **Input**: layered/topology workstation indices + 7+7 sample pages。
- **Command**: `python3 InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/validate_workstation_ui.py --all --screenshot-mode off`
- **Output**: stdout JSON receipt。
- **Observed result**: 16 documents / 32 page runs；365/365 checks passed；0 failing coordinates。

### V-006 — Full layered visual audit

- **Input**: current layered index + 7 sample pages。
- **Command**: `python3 InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_after/capture_layered_workstation_after.py --input-dir InterSubMod/docs/methodology/_assets/layered_workstation --output-dir InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_full_after --metrics InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_full_after/metrics.json --before-metrics InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_before/metrics.json`
- **Output**: `InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_full_after/`。
- **Observed result**: 8 pages × 3 viewports；53 screenshots；669/669 checks passed；source HTML hashes unchanged。

### V-007 — Design QA

- **Input**: archived/current desktop + mobile clean combined comparison images。
- **Command**: visual inspection via local image viewer + required fidelity rubric。
- **Output**: `InterSubMod/design-qa.md` and `InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_after/20260715_三維概念改版後視覺稽核_01.md`。
- **Observed result**: no actionable P0/P1/P2；`final result: passed`。

### V-008 — Claude Code generated-artifact final gate

- **Input**: two builders、generated index、focused/full metrics、visual report、`design-qa.md`；刻意不讀 39 MB sample HTML，避免無效 context。
- **Command**: `claude -p '<final artifact gate>' --model sonnet --effort high --allowedTools 'Read,Grep,Glob'`
- **Output**: stdout read-only verdict。
- **Observed result**: P0=0、P1=0、P2=0，最終 `PASS`。Claude 獨立重算兩組 WCAG contrast 與 cohort counts 均吻合；唯一 cosmetic note 是 full metrics 的 `before_comparison.source_hashes_unchanged=null`（before schema 不同），但 current `source_html_hashes_unchanged` check 為 pass，不影響結論。

## Provenance Footer

- **Commit hash**: `6067568`
- **Build date**: 2026-07-15 Asia/Taipei
- **Skill version**: `/implementation-notes` v0.1
- **Canonical summary SHA-256**: `71da78b69f8afe5fb8e618179ab7b38a6940fdb17be6282f99a6ec4b720e5de7`
