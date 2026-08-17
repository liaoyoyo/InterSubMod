<!--
建立時間: 2026-08-13
目標: 保存 multi-BAM dashboard 與 bounded ingestion contract 最終三路獨立唯讀複核
處理範圍: 7 datasets（6 biological + 1 technical）、artifact/HTML/browser QA、source/receipt/metric chain；不含 full BAM scan、truth/KDE 或 production release
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_input_manifest.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_artifact.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_multi_bam_analysis_overview.standalone.html
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/multi_bam_dashboard_browser_qa.json
-->

用 PREP：**ACCEPT — bounded offline snapshot data-product QA 通過；science 仍 `PARTIAL`、production `NOT COMPLETE`（影響：高，信心：高）。**

# Multi-BAM Dashboard 最終獨立唯讀 Audit Receipt

## 1. 最終交付 hash chain

| Artifact | SHA256 |
|---|---|
| Input manifest | `2e290c06be9e98be480ea913b62b3de30ad827eb9e6894a5a68b4e5dfa53f190` |
| Intake receipt | `0a18c63c35920d9c089ccf651850b1a4f3242caa681bcd8e4b87c447415d130b` |
| Dashboard artifact | `e8e9c2bf91d56de269a3d6868ac1e001e27c35dd93fd43a920e1d0d5f8a7ea0b` |
| Standalone HTML | `224122cfc41ba24acc886f7e04b7f6d3553887be1118125b1f541a9646d0476c` |
| Browser QA receipt | `0b72cb659f918913e2f25f6a6cb997354b031387cab9ca8a40c7fdfcf319b5ef` |

三位 Codex subagent 皆使用唯讀模式；沒有修改 repository。Claude Code Round 5 因 session quota 在任何 Read/Bash 前中止，狀態是 `NOT EXECUTED`，不列為 ACCEPT 或 REJECT 證據。

## 2. Ingestion 對抗審查：ACCEPT

實跑 `32/32` contract regressions，py_compile、`git diff --check` 與 canonical artifact contract 皆 exit 0/PASS。

| Gate | 獨立觀察 |
|---|---|
| Canonical source schema | `$id` 與 SHA256 固定；permissive schema 被拒，無 output |
| Receipt coherence | status、claim ceiling、完整 verification summary、schema/builder/source/topology/output SHA 任一篡改即拒絕 |
| Sidecar/sample binding | source、subject、producer receipt、sidecar/index/native validation 與 VCF inputs 依 dataset 綁定 |
| Chunk policy | 限定 first/middle/last，每段 1 MiB；錯 label/offset/size 拒絕 |
| I/O safety | intake output/receipt 與 dashboard output 皆不可覆寫 explicit 或 implicit inputs；atomic replacement |
| Dynamic identity | all-MATCH：`strict=7, drift=0`；mixed All=`1/6`、HCC1395=`1/0`、H1437=`0/1` |
| Dynamic narrative | all-MATCH artifact 的四個舊 drift-only 句子皆 0 hits；availability=`AVAILABLE_BOUNDED` |
| Trackability | 15/15 必要檔已被 `.gitignore` 精確放行；但新檔仍未 stage/commit |

## 3. 來源與數值審查：PASS

最終輕量重算 **682/682 PASS**，artifact contract 另外 PASS，未重讀大型 BAM payload。

- 8/8 unique artifact input pins 與 current bytes 一致，0 conflicting pins。
- Manifest/receipt status、summary、source schema、source manifest、topology、builder 與 output hash 閉合。
- 7 datasets／6 biological／1 technical 一致。
- 42-row tag chart、14-row producer detail rail、14-row identity rail、8-row BAM summary rail 與 exact counts 一致。
- HP/all、HP+PS/all、HP+PS/HP、duplicate/all 的 numerator、denominator、rate 在 All 與 selected views 無錯配。
- Canonical observation 仍是 all-role strict `0/7`、mount-device-only drift `7/7`、no RG `7/7`、cross-directory `1/7`、full BAM SHA `0/7`。

## 4. UI／metrics 審查：ACCEPT

- Browser receipt：**40/40 PASS**；console/page error=0，external HTTP(S) requests=0。
- All＋7 dataset selectors；All 為 7 BAM/7 tag rows，每個 selected state 只剩 1/1。
- 非 HCC1395 的 ISM/lineage rows 為 0，availability 是 `ABSENT_NO_EQUIVALENT_BUNDLE`，沒有借用或以 0 代 null。
- Cluster 從 valid-axis chart 排除，只在 validity rail 保留 `INVALID · DOUBLE-DIPPING` 與 exact N/D。
- Active-k 8-row exact N/D rail、duplicate guard 與七樣本 fixed-label opportunity panel 通過。
- 6 個 audit disclosures 的 `aria-controls`、Enter 開合、focus 與 `aria-live` 正常。
- 1024/512/390/320 px 皆無 page overflow；寬表只在 container 內局部橫捲。

## 5. 允許與禁止的 claim

**允許**：`bounded snapshot data-product QA PASS / science PARTIAL / production NOT COMPLETE`。

**禁止**：full BAM identity、biological identity、primary-read tag coverage、phasing accuracy、caller F1、KDE-corrected effect、platform effect、cluster biology 或 lineage causality。下一階段仍需 truth/HighConf/benchmark receipts、MM/ML/CpG/KDE provenance、逐 BAM depth/N50/phase blocks 與七份等價 ISM/lineage bundles。

---

**Handoff**：此 receipt 只證明當前工作樹的最終 bytes。新產物已 trackable 但未 stage/commit；未完成 commit 前，不宣稱 clean checkout 可重現。
