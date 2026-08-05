<!--
建立時間: 2026-07-15T09:20:00+08:00
目標: 保存全 sSNV 下游分析正式執行前的多 agent 紅隊判定與修正閉環
處理範圍: screen / read-tag provenance / cooccurrence statistics / strict-normal controls / final claim ladder
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/implementation-notes.md；InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/
-->

# 2026-07-15 全 sSNV pre-execution multi-agent red-team

## 判定

**初始判定：NO-GO 作為最終科學報告資料源；GO 僅限全量 InterSubMod extraction 繼續執行。**

四位唯讀 reviewer 一致認為，469,849-site input/extraction 契約本身可保留，但原定 C3/C4 推論會高估證據。正式 focal-ALT screen 尚未開始，因此修正不需要重跑 C++ extraction；所有下游輸出必須在阻塞關閉後由新路徑產生。

## Reviewer 範圍

| Reviewer | Agent ID | 審查面 | 初始 verdict |
|---|---|---|---|
| Laplace | `019f6347-3c2a-76d1-8c9b-fd2a6643e42a` | 統計、permutation、strict | exploratory screen only |
| Confucius | `019f6347-46f9-70f2-a443-543af25574f2` | LongPhase-S input、tag、provenance | input PASS；terminal evidence pending |
| Archimedes | `019f6347-56a1-7cb0-9c18-ce925319e5d8` | 癌症演化、甲基與 claim ceiling | subclone-compatible NO-GO |
| Peirce | `019f6347-634d-7711-b85a-ac286402ccfd` | final integration/report contract | integration draft only |

## 阻塞與處理決定

| ID | 初始問題 | 風險 | 接受的修正 | 關閉驗證 |
|---|---|---|---|---|
| R1 | stable 只鎖 modal K，沒有鎖 read membership | 相同 K、不同 partition 被稱 stable | `modal_assignment_ari_min>=0.8` 納入 stable gate | 對相同 K/低 ARI fixture 必須 fail |
| R2 | 已知 HP/technical confound 未排除候選 | phase/geometry 可偽裝 methyl group | 候選要求 residual unexplained 且 HP/technical confound 均 false | confounded fixture 不得成 C3/C4 |
| R3 | analytic p 預選後才跑 conditional permutation | selection bias；稀疏 Kx2 p 不校準 | 所有 testable pair 先跑 HP-family x PS x strand permutation，再對 conditional p 做 global BH/BY | 稀疏表 spot-check、不可交換列 NOT_IDENTIFIABLE |
| R4 | `independent markers` 只代表相距 20 bp | 誤稱統計/生物獨立 | 正名 `spatially_separated_markers_20bp`；正式敘述不得使用 independent | schema/報告字串檢查 |
| R5 | 兩個 pairwise hits 沒有 same-read joint signature | 可由共享 reads/同一 artifact 重複計數 | coverage-preselected top markers 建同批 complete-read genotype signature，conditional p 後 global BY | joint signature 不通過不得升級 multi-marker candidate |
| R6 | 四狀態永遠容許至少 1 forbidden read | 低深度 25% violation 仍判 compatible | 固定 error ceiling 加 exact binomial/CI；低 power 為 NOT_IDENTIFIABLE | `RR=3,AR=3,RA=1,AA=3` 必須不 compatible |
| R7 | C4 可由不同 partner 拼接 gate | 沒有一條一致證據鏈 | 產生 same-pair BY + four-state witness keys；不得以不同 pair 拼接 | split-witness fixture 必須 fail |
| R8 | C4 對所有資料集強制 HCC1395 replication | 非 HCC 結構性不可評估卻被算失敗 | 技術重現改為獨立 endpoint；HCC pair要求雙平台 conditional-BY、方向/effect相容 | 非 HCC 不進 replication denominator |
| R9 | normal methyl-evaluable 被誤當 focal allele callable | all-UNKNOWN normal 可通過「無 ALT」 | normal REF-only methyl control；另要求 called depth與 REF depth；輸出 ALT/REF/UNKNOWN | all-UNKNOWN fixture 必須 fail |
| R10 | strict 小 permutation 也可輸出 robustness pass | 測試設定被誤當正式證據 | 正式需 permutations>=999、seeds>=10；低於只可 `EXPLORATORY_ONLY` 且 scientific pass=false | CLI/fixture 覆蓋兩種狀態 |

### R3 closure correction (2026-07-15)

R3 的初始建議「所有 testable pair 先跑 conditional permutation，再對 conditional p 做 global BH/BY」經可達解析度檢查後撤回。999 次 permutation 的最小非零 p 值為 0.001；在數萬至數十萬 pair 的 family 中，無法提供足夠解析度給全域 BY，且先依 marginal 訊號挑選再把 conditional p 當全域 family 亦會造成 selection bias。

正式 v2 契約改為：

1. 所有 testable Kx2 pair 使用固定邊際的 Freeman-Halton exact p，並在預先鎖定的 M2 family 內做 global BH/BY。
2. 只有 exact-BY discovery 再執行 HP-family x PS x strand conditional permutation，作為次級 sensitivity gate；此 conditional p 不宣稱具 global FDR。
3. formal pair 必須同時通過 exact-BY、效果量與 conditional sensitivity；不可交換 strata、超過 exact state-space ceiling 或低資訊列為 `NOT_IDENTIFIABLE`，不當成陰性。
4. same-read complete-read joint signature 另以 999 次 post-selection sensitivity 檢查，且不得替代 exact global family。

主 session 獨立驗證：`38 passed in 2.64s`，涵蓋 FFH golden table、M2-only family、state-space fail-closed、conditional gate、joint signature、20 bp spacing、four-state 低深度與 cross-platform fail-closed。
| R11 | strict/cooccurrence/final receipt SHA 鏈不完整 | key 相同、不同 assignment 可被混接 | 每階段 receipt 鎖 input SHA、family size、full scope與 code provenance | digest mismatch 必須在建 output 前 fail |
| R12 | `subclone-compatible` 缺 CN/purity/CCF | molecular mixture 被升級成 cellular clone | 最強 bulk endpoint改名 `multi-marker molecular-haplotype candidate`；CN/CCF與 cellular lineage另列 NOT_EVALUABLE/NOT_SUPPORTED | final schema 不得把此 endpoint命名為 subclone confirmed/compatible |

## Input 與 BAM 已確認事項

- 同次 LongPhase-S producer 的 recalibrated PASS records為 582,820；chr1-22 biallelic sSNV exact subset為 469,849。
- `469,849 = 335,296 TP + 7,745 FP + 126,808 UNASSESSED`，truth只可 post hoc 分層。
- producer VCF、sidecar、raw BAM realpath與 layered snapshot 7/7 綁定一致。
- 原始 BAM沒有被覆寫，也沒有持久化新的 tagged BAM；latest HP/PS保存於 sidecar。
- 既有 `latest_tag_projected_join_audit.json` 只能證明 7,745 FP sites 的 434,759 read rows HP 一致，不能冒充全 469,849-site HP/PS terminal evidence。

## Claim ceiling 重鎖

1. **C1**：focal-ALT read-level methyl heterogeneity。
2. **C1R**：已測量 confound 後仍存在、null/resampling-stable methyl partition。
3. **C2**：與 LongPhase-S PASS callset 的 same-molecule local co-segregation；不等同已確認 somatic marker。
4. **C3**：同批 complete reads 的 multi-marker molecular-haplotype candidate，另需 normal/tumor-REF/strict/four-state guardrails。
5. **C4/C5**：cellular subclone與 linear lineage 需要 allele-specific CN、purity、multiplicity/CCF與 single-cell/colony/spatial/multi-region 等正交證據；缺資料時必須 NOT_EVALUABLE 或 NOT_SUPPORTED。

HCC1395/DORADO只增加技術重現性，`biological n=1`。任何 read-level、HP-family 或 bulk topology 結果不得直接推定 clone數或 linear ancestry。

## Closure 要求

- 所有修正單元測試與完整 topic test suite通過。
- full screen receipt證明 469,849/469,849 有結果或明確不可評估原因。
- cooccurrence只對 stable subset，但 partner oracle必須與全 PASS VCF幾何 exact reconciliation。
- candidate、strict、tumor-REF、normal、CN/CCF annotation各自保留明確 denominator與不可評估原因。
- 修正後由新的多 agent reviewer重算分母與個案，再由外部 Claude Code唯讀審查；initial NO-GO不得從文件刪除。
