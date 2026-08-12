<!--
建立時間: 2026-08-13 01:41
目標: 保存 Claude Code 對 HCC1395_v1/v3 的第一輪獨立 red-team 審查
處理範圍: 唯讀 code/output review；未授權任何修改
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/pre-decision-audit.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round2_review.md
-->

# Claude Code Round 1：獨立 red-team 審查

- **Reviewer**: Claude Code / `claude-opus-4-8`
- **Session ID**: `4574cff8-4840-43c6-a006-9dd8d73e6182`
- **Mode**: `plan`、Read/Grep/Glob/Bash only；Edit/Write/NotebookEdit disabled
- **Scope**: `/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/`，並以 `HCC1395_v3/` 對照
- **Result**: success；18 turns；無 permission denial

## Verdict

Claude Code 判定：

> v1 只適合作為「描述性、單樣本、拓撲骨幹內部自洽」的 observation UI，及有條件的 internal QA；不能作 caller F1、subclone、甲基因果或跨樣本研究結論。

## 與本地稽核一致的主要發現

1. **P0 provenance 不閉合**：receipt 只列 topology、MLHP、ISM directory 三項；strict edges、LCA receipts、lineage paths、annotations 未列 FileRef。ISM directory 又是 size=0、無 hash。
2. **C10 是科學有效性 gate，不是普通守恆式**：cluster 分組與檢定距離同源，20,903/20,904 raw p≤0.05 是 double-dipping by construction，不能當證據。
3. **C8 只比總數**：hidden count 相等是必要但不充分；不同 region 的集合交換仍可能假 PASS。
4. **C12 是 round-trip consistency**：component 本來由 passing endpoint edges 構造，因此 0/11,590 broken 主要證明回算一致，不是獨立生物學證據。active-only 394/11,590 = 3.40% 才是需顯示的橋接依賴指標。
5. **LCA A/B 交集策略正確，但 provenance 與檔數文案不足**：應記 pre/post/shared counts，不能把 shared 22 寫成「兩套各 22」；`lca_resolved` 不應與 net `lv_written` gain 混同。
6. **v3 成本不成比例**：有價值的新增主要是 C8 lineage_paths 接線；ISM linkage 僅 +109 loci，而整體體積增加約 2.9 GB，且 provenance 問題仍在。
7. **現有 test 假綠燈**：`test_receipt_has_input_hashes` 只檢查已在 receipt 的 items，且跳過 size=0 directory；完全看不到未列入的 4 個 capability sources。
8. **NCHS suppression fixture 不足**：cell n<30、margin n≥30 時仍可能顯示百分比；現有小 fixture 沒覆蓋這個反例。

## Claude Code 建議的最小新增指標

| 指標 | 明確公式 | 角色 |
|---|---|---|
| provenance 完整率 | 有 full SHA-256 的 consumed files / 全部 consumed files | provenance gate |
| shard 守恆 | Σ L2 regions = receipt regions；Σ L5 loci = panels made | mechanical invariant |
| 非自參軸顯著率 | p≤0.05 / valid tested，cluster 軸分開隔離 | heuristic warning |
| C8 集合一致率 | hidden vertex 集合相等 regions / 可比 regions | cross-source invariant |
| ISM 選擇偏誤 | 各 k 的 ISM availability + Cramér's V | missingness diagnostic |
| active-only bridge dependence | active-only broken blocks / original-sSNV blocks | structural diagnostic |
| multi-region consistency | HP/status 一致的 multi-region sSNV / 全 multi-region sSNV | duplicate-grain QA |

## 多樣本 contract

- 每樣本 receipt 必含 sample、capability states、linkage numerator/denominator、ISM window、input hashes、selfcheck。
- ISM window、endpoint-edge schema、region-linkage rule 任一不一致即拒絕聚合。
- 缺能力標 ABSENT，不插補、不把 missing 當 0。
- 先報逐樣本，再報 biological-sample macro median/IQR；HCC1395_DORADO 是 technical replicate，不能增加 biological n。
- 任一樣本的非自參軸出現極端顯著率或 sample identity mismatch，cohort status 應 BLOCKED。

## 第一輪 P0/P1 建議

1. 所有 consumed sources 記完整 FileRef；provenance completeness 必須 100%。
2. receipt 加 generator git SHA、dirty、CLI、site/sample identity、linkage contract。
3. 把 C10 從「守恆等式」語意分離，cluster 軸標示為 invalid-as-evidence。
4. C8 從 count equality 升級 per-region set equality。
5. 補負向 fixtures：壞樹、hidden 集合交換、small cell/large margin、缺 provenance。
6. 多樣本先過 manifest compatibility gate，再生成或聚合。

## 執行命令

```bash
claude -p '<red-team prompt>' \
  --session-id 4574cff8-4840-43c6-a006-9dd8d73e6182 \
  --output-format json --permission-mode plan \
  --tools Read,Grep,Glob,Bash \
  --disallowedTools Edit,Write,NotebookEdit \
  --max-budget-usd 5 --effort high \
  --add-dir /bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1 \
            /bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v3
```

實際輸出片段：`is_error=false`、`stop_reason=end_turn`、`permission_denials=[]`、`subtype=success`。

