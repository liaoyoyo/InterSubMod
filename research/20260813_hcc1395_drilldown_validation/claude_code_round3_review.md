<!--
建立時間: 2026-08-13
目標: 保存 Claude Code 對 post-fix source 與回歸測試的獨立唯讀驗證
處理範圍: 指定 drilldown source、兩個 pytest suites；不讀／不改 legacy bundles
關聯檔案:
  - InterSubMod/tests/test_drilldown_multisample_contract.py
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round2_review.md
-->

# Claude Code Round 3 — post-fix verification

## 執行契約

- 模式：唯讀；`Read/Grep/Glob/Bash`，明確禁用 Edit/Write/Web。
- 測試：`python3 -m pytest tests/test_drilldown_contract.py tests/test_drilldown_multisample_contract.py -q`
- 實得：**37 passed**，無 skip/fail。

## 獨立裁決

| 項目 | 判定 | Claude Code 證據摘要 |
|---|---|---|
| sample fail-closed | **PASS** | ISM、LCA、lineage、strict edges 與 topology/MLHP identity checks 均落地；wrong-sample fixtures 全通過 |
| receipt validation / dirty / provenance | **PASS** | generator commit/dirty/script hash 與五道 validation gates 已機器化；`scientifically_validated=false`、`truth_set_validation=false` |
| LCA causal gate | **PASS** | 按染色體比較 input identity 與非 LCA stats；mismatch 會 PARTIAL，C11 會 FAIL，不可宣稱純 LCA effect |
| dynamic selfcheck | **PASS** | FAIL→blocked、SKIP 不等於 PASS、只有全 PASS 才稱通過；前端同步三態 |
| mobile source fix | **PASS（範圍內）** | 980/560 px breakpoint、單欄 layout、responsive SVG 已存在；此輪本身是 code review，實際 390 px screenshot 由主流程另驗 |

## 尚未關閉的技術債

1. `shell.py` 有舊的「待 P3 實作」placeholder 文案，與已渲染的 selfcheck/cooccur 不同步。
2. `build_drilldown_dashboard.py` 的 IGV path 組裝仍有 `if False` 死分支。
3. output inventory 只記 count/bytes，不對數萬個輸出檔逐一 hash；輸入層仍有 file hash。
4. `.locus-strip` / `.table-wrap` 的極窄螢幕體驗主要依賴既有 overflow，尚未作跨瀏覽器／touch 深度驗證。

## Claim boundary

本輪只驗證「source contract 與 tests 已修好」，**不會回溯改變 legacy v1/v3 的 BLOCKED / DO-NOT-CITE 判定**。
