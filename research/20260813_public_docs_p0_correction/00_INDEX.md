<!--
建立時間: 2026-08-13
目標: 將 2026-08-12 的 158-row 公開文件稽核轉成可重跑、fail-closed 的修正與發布 Gate
處理範圍: Task B + D；34 個 P0 source corrections 與完整 158-claim remediation registry
release_baseline: ddd8909a838318d8a77969313e9561c8ff9d01c2
status: P0_SOURCE_READY__ABOUT_C108_LIVE_CONFIRMED__ALL_SOURCE_BLOCKED__PUBLICATION_BLOCKED__RELEASE_BLOCKED
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
  - InterSubMod/research/20260813_public_docs_p0_correction/claim_remediation_registry.json
  - InterSubMod/research/20260813_public_docs_p0_correction/claim_remediation_build_receipt.json
-->

# 2026-08-13 公開文件 claim 修正與發布 Gate

> **結論先行：**34 個 P0 中，33 個本機可控 claim 已通過 source guard；GitHub About（C108）已換成 bounded description、移除 `subclone`／`phylogeny` topics並re-fetch，receipt已保存。初始24個P1/P2中，C047/C048/C089/C114已以可追溯的有界證據收旂，其餘20個仍為醒目UNVERIFIED；default branch／Wiki／Pages尚未發布與抓回驗證。因此只允許說 **bounded source/evidence closure**，不可說全source-ready、publication-ready或release-ready。

## 任務與 claim ceiling

- Task type：B Comprehensive validation + D External handoff。
- 服務 G1、G3、G4、G5。
- 分子共現是 physical-molecule observation；不是 cellular co-membership truth。
- candidate topology 是 local、recurrence-allowed、model-conditional；確認的 cellular subclone = 0、linear ancestry = 0。
- methylation 僅為 genetic grouping 後的 association；CN/LOH 尚未整合。
- 63,506 / 71,955 = 88.2579% 是 model-conditional graph shape，不是生物正確率或 prevalence。

## 輸入 → 輸出 → 驗證

1. 輸入：InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
   → 輸出：InterSubMod/research/20260813_public_docs_p0_correction/claim_remediation_registry.json
   → 驗證：158 IDs 完全相等、無重複／遺漏。
2. 輸入：README EN/ZH、Quickstart、Summary、8 Wiki sources、17 Pages HTML、page-07 generator
   → 輸出：受控source corrections與可重跑guard
   → 驗證：P0 guard 34/34、27 target documents、77 target rules；42個cross-document guards在34個public source files執行500次document checks、285 required anchors、1,217 forbidden anchors、0 errors。
3. 輸入：2026-08-01 denominator_registry.tsv
   → 輸出：7 個固定 denominator assertions
   → 驗證：469,849、255,752、170,131、71,955、63,506、811、3 與精確百分比一致。
4. 輸入：branch/commit-specific capability statements
   → 驗證：LongLineage private baseline/main snapshot `5daf50f`、private public-preview candidate `b9aaa12`、revision-scoped `longlineage-tag-bam`、`NOT_READY`與P3/P4/P5/P7/P8 BLOCKED；InterSubMod 199-column schema anchors皆存在。
5. 輸入：`scripts/p0_claim_registry.json` 的跨文件 guard contract
   → 驗證：current schema 不得回退 59/193 欄、ISM 與 exact-PS solver 不得混用、66.5%／63,506／unique-tree 必帶 grain、README 不得手抄 test count、Python ≥3.10、fixture 必須是 Git stage-0 regular blobs、公開來源不得使用 `/develop/` raw URL。
6. 輸入：GitHub About API response
   → 輸出：`github_about_c108_receipt.json`
   → 驗證：description精確等於bounded wording、topics不含`subclone`／`phylogeny`、repository identity與re-fetch時間存在；C108為`CONFIRMED_WITH_LIMITS_AFTER_REFETCH`。
7. 發布後才可執行：re-fetch default branch／Wiki／Pages
   → 驗證：live bytes與指定commit/hash一致；本輪未執行，故這三面仍為UNVERIFIED_AFTER_20260812_LOCKED_SNAPSHOT。

## Registry 結果

| 維度 | 數量 |
|---|---:|
| claims | 158 |
| P0 / P1 / P2 / P3 | 34 / 20 / 35 / 69 |
| CONFIRMED | 69 |
| CONFIRMED_WITH_LIMITS | 69 |
| UNVERIFIED | 20 |
| 本機 SOURCE_READY P0 | 33 |
| SOURCE_EDITED_REVIEW_REQUIRED P1/P2 | 20 |
| EXTERNAL_LIVE_CONFIRMED_WITH_LIMITS | 1（C108） |

注意：current_verdict 是 evidence adjudication；source_status 是 checked-in source；live_status 是重新抓回的發布物。三者不可互相代替。

## Gate

| Gate | 狀態 | 原因 |
|---|---|---|
| P0_SOURCE_READY | PASS | 33個local P0 source claims通過semantic anchors；C108 About有bounded live receipt |
| SOURCE_READY | BLOCKED | 其餘20個P1/P2 source edits尚缺逐claim evidence closure |
| PUBLICATION_READY | BLOCKED | default branch、Wiki、Pages未發布／re-fetch；About C108單獨通過不會升格整體 |
| RELEASE_READY | BLOCKED | publication gate 與完整 handoff release gates 尚未通過 |

HTML companion的repository-owned Chromium runner綁定精確allowlist：17頁public explain + 4頁handoff，共21頁×desktop/mobile/no-JS/print=`84` cases，並檢查21份standalone SVG、inline SVG、local link、fragment target、水平overflow、page/console error與HTTP(S) runtime request。最終驗收只接受clean-source receipt；此receipt只驗呈現與runtime safety，不重算science。

## 重跑命令

~~~bash
python3 research/20260813_public_docs_p0_correction/scripts/run_claim_validation.py

# 上述 canonical runner 依序重建 registry，再執行以下分解步驟：
python3 research/20260813_public_docs_p0_correction/scripts/build_claim_remediation_registry.py
python3 research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py
python3 research/20260813_public_docs_p0_correction/scripts/validate_claim_remediation_registry.py
python3 -m unittest discover -v -s research/20260813_public_docs_p0_correction/tests -p 'test_*.py'
~~~

`p0_claim_registry.json` 是 guard 的唯一權威輸入；builder 必須先驗證它，並把其 SHA-256 寫入 generated registry 與 receipt。35個 fail-closed tests另覆蓋 duplicate／missing ID、illegal verdict、missing evidence、denominator drift、capability commit drift、guard-hash drift、59/193 current schema、solver 混用、裸分母、手抄 test count、Python 3.9、未追蹤 fixture、`/develop/` raw URL，以及把 local source 偷升 live confirmed 的錯誤。

## 尚未完成

- 初始 24 個 P1/P2 中尚有 20 個 source edits需依 registry 的 targets／inventory_minimum_rewrite補專屬guard、evidence closure與review；2026-08-12 inventory仍為不可變 audit input，當前 verdict以registry為準。
- main merge、Wiki push、Pages deploy，及四個 live surfaces 的 commit/hash re-fetch。
- default branch／Wiki／Pages live state未驗證前，任何local SOURCE_READY都不得改成published/resolved；C108僅依其獨立receipt例外。
