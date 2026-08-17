<!--
建立時間: 2026-08-13 17:36 +08:00
目標: 預註冊 InterSubMod GitHub/Pages 最新資訊修正與 CCU 唯讀改進清單的完整驗證循環
處理範圍: 全部 24 個既定 InterSubMod P1/P2 問題 claim、17 個 Pages、全部 GitHub 公開入口；CCU 僅整理既有 32 findings 的重點清單
關聯檔案:
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/pre-decision-audit.md
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
  - InterSubMod/research/20260813_public_docs_p0_correction/00_INDEX.md
-->

# 2026-08-13 InterSubMod 公開介面資訊更新循環

> **Task B / comprehensive local-source correction**：修正 InterSubMod GitHub 與 Pages 的本地公開來源；CCU 只交付重點改進清單，不改動 CCU。

**最終狀態：`VALIDATED_LOCAL / PUBLICATION_PENDING`。** GitHub About 已在 live state 解決；default `main`、Wiki 與 Pages 仍待另行授權發布。

## 服務目標

- G2：清楚區分 longphase-s／longphase-to、caller-VAF、read-AF 與 CN/LOH 適用邊界。
- G3：把 read-level epigenetic 限定為 molecular evidence 與 association，不宣稱 cellular lineage／causality。
- G4：讓 EN/ZH、Wiki、Pages、source/generator 與版本數字一致。
- G5：以 claim registry、可重跑 QA、HTML 報告與 ledger 建立外部可驗證交付。

## Frozen scope

- 24 個問題 claim：C047、C048、C089、C090、C110、C111、C114–C123、C132、C133、C145–C147、C156–C158。
- GitHub source：`README.md`、`README.zh-TW.md`、`QUICKSTART.md`、`README_PROJECT_SUMMARY.md`、`docs/wiki/*.md`。
- Pages source：`docs/explain/index.standalone.html` 與 01–16 共 17 頁；若有 generator，generator 是必同步對象。
- CCU：只讀 `lab_tutorial_delta_reaudit.tsv` 與既有 receipt；不修改 CCU repo、live site 或 `ccu_source_correction_9eb1618.patch`。
- Remote：本輪不 push、merge、改 GitHub About、發布 Wiki 或部署 Pages。

## Pre-registration（confirmatory）

| 預測 H | 否證條件 | Decision threshold |
|---|---|---|
| H1：24 個剩餘問題 claim 全部可被閉合或誠實標為 external/unverifiable | 任一 claim 無 source mapping、修正文案、證據或 disposition | 24/24 claim registry 完整；0 orphan |
| H2：修正後 InterSubMod 公開來源不超過 frozen scientific ceiling | residual overclaim、錯誤分母、過時 test/code count 或無 receipt 性能數字仍存在 | claim guard 0 error；負向 mutation probe 必須 exit nonzero |
| H3：17 個 Pages 與所有 inline SVG 維持結構正確且可近用 | 任一 HTML FAIL、broken local link、SVG 無 title/desc、mobile/print/no-JS 失敗 | 17/17 page PASS；37/37 SVG title+desc；browser QA 0 error/overflow |
| H4：CCU 只產生清單且數量可閉合 | 修改任何 CCU material，或 32 findings 無法分成 prior-target 13 + remaining 16 + prior-resolved 3 | CCU source mutation=0；32=13+16+3；清單明示 16 remaining |

## Step → Verify

1. 凍結 24 claim 與 occurrence mapping → 驗證：TSV/registry 恰為 24 unique IDs，對應 inventory 母體 24/24。
2. 修正 GitHub 公開文件 → 驗證：EN/ZH、Wiki、Quickstart/Summary 的 claim guard、連結與命令檢查均 exit 0。
3. 修正 Pages 與 SVG → 驗證：17 頁結構 0 FAIL，37 個 SVG 均含非空 `title`/`desc`，generator-source parity PASS。
4. 建立 CCU 唯讀清單 → 驗證：16 remaining 逐項含問題、證據需求、最低修正與優先級；CCU 檔案 0 變動。
5. 整合 claim receipt 與 HTML 報告 → 驗證：24/24 disposition、P0 guard 未回歸、MD/HTML 核心數字 parity。
6. 執行 fresh-reader／run-evaluator → 驗證：scope、claim ceiling、local-vs-live 狀態與 CCU no-change 可由只讀者重建。
7. 寫入 evidence ledger → 驗證：entry 含命令、輸入、輸出、hash、verdict 與限制。

## 預期輸出

- `InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/github_correction_receipt.md`
- `InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/pages_correction_receipt.md`
- `InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/20260813_CCU教學站重點改進清單_01.md`
- `InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/scripts/p1_p2_claim_registry.json`
- `InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/scripts/validate_public_p1_p2_claims.py`
- `InterSubMod/docs/reports/validated/2026/08/20260813_InterSubMod公開資訊更新與CCU改進清單_01.md`
- `InterSubMod/docs/reports/validated/2026/08/20260813_InterSubMod公開資訊更新與CCU改進清單_01.standalone.html`

## 最終結果

| Gate | 結果 | 邊界 |
|---|---:|---|
| 公開 claim 問題閉合 | 58/58 | P0 34/34；P1/P2 24/24；21 修正、3 移除 |
| 全體 claim verdict | 158 | 69 confirmed、31 confirmed with limits、28 needs correction、26 contradicted、4 unverifiable |
| Pages 結構 | 17/17 pages；37/37 SVG | 0 structural error；8 組圖例 SVG/PNG 重生 PASS |
| Browser QA | 68/68 | desktop、390 px、no-JS、print；0 overflow/error |
| Core regression | 270/270 tests | 39 suites；只是版本綁定軟體回歸 |
| CCU 清單 | 32/32 分類 | 13 prior patch targets + 16 remaining + 3 prior resolved；0 CCU 改動 |
| 主交付 | official portable verifier PASS | 30 rendered blocks、3 charts、4 tables、6 metrics；1440/390 px PASS |

## 交付與收據

- 主 HTML：`InterSubMod/docs/reports/validated/2026/08/20260813_InterSubMod公開資訊更新與CCU改進清單_01.standalone.html`
- 輔助 Markdown：`InterSubMod/docs/reports/validated/2026/08/20260813_InterSubMod公開資訊更新與CCU改進清單_01.md`
- CCU 重點清單：`InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/20260813_CCU教學站重點改進清單_01.md`
- Claim guard：`InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/claim_guard_receipt.md`
- Pages QA：`InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/pages_structural_qa.json` 與 `pages_browser_qa.json`
- Report build：`InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/report_build_receipt.md`

## 結論邊界

- 可對外聲稱：同一 physical molecule 上的 allele/HP/PS/MM/ML called evidence；局部、允許 recurrence、模型條件下的 mutation-state candidates；genetic-pattern-conditioned regional methylation association。
- 不可升格：confirmed cellular clone/lineage、canonical CN/LOH-corrected CCF、causal methylation/function。
- 公開文件 QA 不是新的生物效能實驗；`validation_evidence_eligible=false`。
- 本輪未 commit、push、merge 或 deploy，也沒有改動 CCU source、patch、remote 或 live site。
