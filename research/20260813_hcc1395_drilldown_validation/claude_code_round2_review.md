<!--
建立時間: 2026-08-13
目標: 保存 Claude Code 第二輪 pre-fix red-team challenge 與修補驗收條件
處理範圍: HCC1395 v1/v3 legacy bundles、當時 source state、sample contamination / provenance / UI / tests
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round1_review.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round3_review.md
-->

# Claude Code Round 2 — pre-fix challenge

> 時序註記：本輪讀到的是修補進行中的 source 與 legacy v1/v3；其 runtime 反證用來驅動後續修正，不代表 post-fix 狀態。最終是否關閉由 Round 3 與本地重跑決定。

## Verdict

**v1/v3 均為 `QUARANTINE / DO-NOT-CITE`。** 既有 bundle 早於 source 修正，不能作 truth validation、causal LCA、跨樣本或 external handoff。當時最嚴重反證是 COLO829 probe 仍可能沿用 HCC1395 的 ISM/LCA path。

## ACCEPT / REVISE

| 項目 | Round 2 裁決 | 驗收條件 |
|---|---|---|
| sample identity fail-closed | ACCEPT, P0 | topology/MLHP/ISM/LCA/lineage/strict edges 都不得借用其他 sample |
| linkage threshold | REVISE | 身分錯誤應由 identity gate 禁用，不能拿低 linkage 冒充生物性缺失 |
| receipt v2 | ACCEPT, P0 | commit/dirty/script hash/argv/input refs/inventory；status 必須機器推導 |
| LCA comparison | ACCEPT with revision | per-chrom shared set、input identity、non-LCA stats 分開；4.969× 不作 causal claim |
| dynamic selfcheck authority | ACCEPT, P0 | FAIL→BLOCKED；SKIP→INCOMPLETE；只有全 PASS 才能稱通過 |
| fake interaction / CN / mobile | ACCEPT | 無 filterMap 不可宣稱可點；移除硬編碼 CN；修正 390 px overlap/overflow |
| synthetic regression | ACCEPT, P0 | wrong-sample、LCA identity、receipt、panel bytes、C12 denominator 必須有反例 |
| immutable legacy bundles | ACCEPT | 不覆寫 v1/v3；新結果只能落新路徑 |

## 多樣本 aggregation 三條 fail-closed 契約

1. sample identity、ISM window、edge schema 任一不相容即拒絕該層聚合。
2. 逐樣本保留 numerator/denominator；跨樣本不 pool p-value，macro 使用 biological-sample median/IQR。
3. missing 是第三態，不等於 0 或不顯著；technical replicate 不增加 biological n。

## 對主流程的影響

本輪促成後續 source 加入：sample-aware routing、per-file identity、LCA causal gate、receipt v2、dynamic authority、panel byte totals、mobile contract 與 synthetic tests。Round 3 另行檢查這些項目是否真正 landed。
