<!--
建立時間: 2026-07-15 Asia/Taipei
目標: layered workstation 三維概念改版的設計對照與交付前 QA
處理範圍: archived HCC1395 reference vs canonical v5 index + 7 sample pages；desktop/tablet/mobile/narrow
關聯檔案:
  - InterSubMod/docs/methodology/_assets/topology_workstation/HCC1395.html
  - InterSubMod/docs/methodology/_assets/layered_workstation/index.html
  - InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_after/metrics.json
-->

# Design QA — Layered Workstation Concept Dimensions

## Comparison target

- **Source visual truth**: `InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_after/reference_desktop_three_facets_clean.png` and `reference_mobile_three_facets_clean.png`。
- **Implementation screenshots**: `InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_concept_dimensions_after/after_desktop_three_dimensions_clean.png` and `after_mobile_three_dimensions_clean.png`。
- **Combined evidence**: `comparison_reference_vs_after_desktop_clean.png` and `comparison_reference_vs_after_mobile_clean.png` in the same directory。
- **Viewport**: desktop 1440×1000；tablet 768×1024；mobile 390×844；narrow 320×720。
- **State**: HCC1395 initial sample overview；focused evidence additionally covers incomplete、multi-shape、no-primary region detail and candidate network。
- **Browser**: Chromium 147.0.7727.15, headless, local `file://` pages。

## Findings

No actionable P0/P1/P2 differences remain.

- **Information architecture**: 舊頁可取的「shape／determinacy／location」三軸形式已保留；實作改為 canonical v5 定義與數字，並置於 chromosome grid 前。舊 clone-first 類別、舊 A/B/C determinacy 與 archive 數字沒有移植。
- **Typography**: 既有 sans/mono 階層、數字 tabular rhythm 與中英輔助標籤一致；結果列三維摘要由 9.5px 提升至 11px。320px 未見截斷或不可讀的 CTA。
- **Spacing and layout**: desktop 三卡等寬；768/390/320 依序單欄，卡片資訊順序不靠 CSS `order` 改寫。Region detail 的三維卡固定在 assertion/evidence 前。
- **Colors and tokens**: 沿用 current blue/magenta/cyan semantic accents；index muted/panel contrast 經公式驗證為 5.56:1，sample muted/white 為 6.08:1。
- **Image/asset fidelity**: 此畫面是純資料工作站，參考頁的 emoji 與舊彩色 pills 不是 canonical asset；沒有以 CSS art、手工 SVG 或 placeholder 偽造資產。Candidate network 沿用既有資料生成 SVG，不是裝飾替代品。
- **Copy/content**: 每卡同時回答名詞、目前 dataset 觀察、下鑽入口與不可推論邊界。`Incomplete` 明示「尚未評估」，`no-primary` 明示「不適用」，位置明示非 hotspot/enrichment。
- **Interactions/states**: 三個 CTA 會清除舊 evidence/signal/query，更新 hash、焦點與結果；exact、shape-only、multi-shape、incomplete、no-primary 五種狀態均有正確 detail。
- **Accessibility**: 44px hit target、focus outline、reduced-motion、heading hierarchy、hidden JSON、keyboard selection、mobile detail reveal 皆通過；topology strip 設為 decorative，文字 legend 保留完整資訊。

## Comparison history

1. **Baseline finding — P1 information hierarchy**: current 頁的 C/Topo 解釋位於完整 22-chromosome grid 之後；390px 需先越過約 1,496px 才能建立概念。
   - **Fix**: 新增三維導讀並移到 genome grid 之前；region row 與 region detail 同步呈現三維。
   - **Post-fix evidence**: clean desktop/mobile combined comparisons；`metrics.json` 542/542 checks passed。
2. **Implementation review — P1 stale-filter interaction**: 初版 CTA 會保留 evidence/signal/query，可能得到 0 筆交集；位置卡的 W_primary 統計與 W_tree 落地未明示。
   - **Fix**: CTA 重置次要篩選與舊 detail；位置卡、status 與 CTA 明示兩個分母；endpoint distance 與 inclusive span 分欄。
   - **Post-fix evidence**: `topology_action_resets_filters`、`position_action_denominator` 與五類 detail assertions 全過。
3. **Harness review — false failure**: full runner 只接受 `overflow-x:hidden/clip`，即使實測 overflow=0 仍失敗，會鼓勵用 CSS 掩蓋問題。
   - **Fix**: 測試契約改為 visible 只在實測 0px 時允許；production 不隱藏 overflow。
   - **Post-fix evidence**: 8 pages × 3 viewports、53 screenshots、669/669 checks passed；source hashes unchanged。

## Open Questions

- 無阻擋問題。Screen-reader 實機與 200% browser zoom 未納入本次自動化，列為外部驗證缺口，不構成目前 P2。

## Implementation Checklist

- [x] Cohort index 三維導讀與 7 dataset chromosome leaders。
- [x] Sample overview 三維導讀、glossary 與 CTA。
- [x] Region list / detail 三維摘要。
- [x] 7/7 hash-bound rebuild + inline JS syntax check。
- [x] Chromium 1440/768/390/320 專項 QA。
- [x] Chromium 16-document generic regression。
- [x] Chromium 53-screenshot full regression。
- [x] Claude Code Sonnet read-only review；P0/P1=0。

## Follow-up Polish

- P3：若未來加入非 autosomal 或外部 annotation，可再評估真實 ideogram；目前不以不完整示意圖增加誤讀風險。

final result: passed
