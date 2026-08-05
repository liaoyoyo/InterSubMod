<!--
建立時間: 2026-07-12
目標: 稽核 HCC1395 與 HCC1395_DORADO 的 VAF-supported 推測樹／shape 比對是否正確整合至主報告
處理範圍: chr1-22 exact-coordinate complete-both historical layered-v2 regions；artifact / Markdown / HTML 唯讀核對
關聯檔案:
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_vaf_pair_metrics.tsv
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_vaf_pair_summary.json
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_vaf_pair_checks.tsv
  - InterSubMod/docs/reports/in_progress/2026/07/20260712_HCC1395兩技術資料粗拓撲與癌症基因藥物一致性驗證_01/artifact.json
狀態: validated
-->

# HCC1395 pair VAF 推測樹比對與主報告最終稽核

用 Claim–Evidence–Boundary：先核對數字，再核對每個數字能代表什麼，最後檢查是否越界成 truth／accuracy／posterior。

## 最終 verdict：PASS

主報告的 artifact、Markdown 與 HTML 均與凍結的 VAF metrics／summary／checks 一致：

- **77.74%** 明確標成「兩側原始 `Topo>1`、各自由 VAF 第一順位集合縮至單一 **rooted-unlabeled shape** 後，允許 HP1↔HP2 交換的 shape agreement」；不是 exact-tree accuracy。
- **37.32%** 明確標成「兩側都實際使用 same-HP read-AF heuristic、且各自只剩一個 genomic-position＋candidate-ID **mutation-labeled exact forest** 後，允許 HP1↔HP2 交換的 agreement」；仍是推測，不是 biological truth。
- `truth`、`accuracy`、`posterior` 等詞在三種報告成品中只出現在明確的否定句或 claim ceiling，未發現正向誤稱。

## 1. 計算定義

### 1.1 VAF 第一順位不是機率真值

每個 primary HP unit 內，先用 reads 計算各 sSNV 位點的：

`read-AF = ALT reads / (REF reads + ALT reads)`

候選樹的 exact ordering score 為所有可比較祖先／後代位點的：

`Σ [read-AF(ancestor) − read-AF(descendant)]`

正式第一順位集合取 exact `Fraction` score 的 argmax；同分候選全部保留。temperature=0.05 softmax 只描述同一 candidate set 內的相對濃度，**不是 calibrated posterior，也不是候選樹為真的機率**。

### 1.2 Exact forest 與 shape 是不同 endpoint

- **Exact forest**：每個 HP component 使用 `(ordered genomic positions, exact-top candidate IDs)`；保留 hidden-support state pattern 與 mutation labels。HP swap sensitivity 只移除 HP1/HP2 名稱，不改 genomic positions。
- **Shape**：只保留 rooted-unlabeled branching skeleton；移除 mutation labels，因此不能識別哪個 sSNV 是祖先。
- 原始 `Topo=1` 的 shape 本來就由結構唯一決定，不需要 VAF。
- 原始 `Topo>1` 只有在 region `evaluable=True`、所有 ambiguous units 均可評估，且 `n_top_shapes_joint_exact=1` 時，才允許 materialize VAF-selected single shape；否則 fail closed。

## 2. 核心分母、分子與主報告對帳

所有 percentage 均由顯示的整數分子／分母重算。

| Endpoint | 分母 | Ordered | HP-swap tolerant | 主報告語意 |
|---|---:|---:|---:|---|
| 未排名完整 exact 候選樹集合 | 5,720 | 1,129 / 5,720 = 19.74% | 2,014 / 5,720 = 35.21% | VAF 前 candidate-set baseline |
| VAF exact-top 候選集合，保留 ties | 3,236 | 548 / 3,236 = 16.93% | 992 / 3,236 = 30.66% | 兩側都實際用 VAF；推測候選集合 |
| VAF 唯一 exact forest | 2,543 | 524 / 2,543 = 20.61% | **949 / 2,543 = 37.32%** | genomic-position＋candidate-ID 的 heuristic exact forest |
| 結構-first＋VAF single shape | 5,168 | 2,317 / 5,168 = 44.83% | 3,667 / 5,168 = 70.96% | 原 Topo=1 用結構；原 Topo>1 才需 VAF |
| 兩側皆 VAF-evaluable 且 single shape | 3,014 | 1,409 / 3,014 = 46.75% | 2,454 / 3,014 = 81.42% | rooted-unlabeled shape；分母有條件選擇 |
| 兩側原始皆 Topo>1，VAF 各縮至 single shape | 2,089 | 919 / 2,089 = 43.99% | **1,624 / 2,089 = 77.74%** | 純 VAF shape-rescue subset；不是 exact tree |

主報告 artifact snapshot 六列與 `hcc1395_vaf_pair_metrics.tsv` 的 denominator、ordered numerator、swap numerator 與未四捨五入比例逐列完全相同：**6/6 PASS**。

## 3. 官方 conservation 與資料 QA

| Check | Observed | Expected | 結果 |
|---|---:|---:|---|
| Exact-coordinate matches | 6,252 | 6,252 | PASS |
| Exact-coordinate complete-both | 5,720 | 5,720 | PASS |
| Pair ranked-unit candidate detail join | 10,446 | 10,446 | PASS |
| HCC1395 complete regions | 6,940 | 6,940 | PASS |
| HCC1395_DORADO complete regions | 6,750 | 6,750 | PASS |
| HCC1395 structural-or-VAF shape selectable | **6,798 = 3,438 structural Topo=1 + 3,360 VAF-selected** | 6,798 | PASS |
| HCC1395_DORADO structural-or-VAF shape selectable | **6,082 = 2,444 structural Topo=1 + 3,638 VAF-selected** | 6,082 | PASS |
| Pair both exact forest selectable | 5,347 | 5,347 | PASS |
| Pair both unique exact forest | 4,602 | 4,602 | PASS |
| Pair both actually VAF-evaluable | 3,236 | 3,236 | PASS |
| Pair both VAF unique exact forest | 2,543 | 2,543 | PASS |
| Pair both structural-or-VAF shape selectable | 5,168 | 5,168 | PASS |
| Pair both original Topo>1 and VAF single-shape | 2,089 | 2,089 | PASS |
| Fixed-T1 shape IDs converted to official TS namespace | 6,554 / 6,554 | 0 mismatch | PASS |

凍結 checks 為 **20/20 PASS**。主報告使用的 `data/` 複本與 agent 原始產物之 metrics、regions、checks 為 byte-identical；summary 移除輸出路徑欄位後 semantic-identical。

## 4. Artifact、Markdown、HTML 語意核對

| 成品 | 77.74% 標示 | 37.32% 標示 | 否定越界 claim | 結果 |
|---|---|---|---|---|
| `artifact.json` | `rooted-unlabeled selected shapes`、`branching skeleton`、`not exact-tree accuracy` | `HP-specific genomic-position tuples plus exact-score-maximizing candidate IDs`、heuristic inference | `not calibrated probabilities or biological truth` | PASS |
| Markdown | 「shape 不含 mutation direction」；77.74% 不是 exact-tree accuracy | 「mutation-labeled exact forests」；raw read-AF heuristic | 「不是 accuracy、posterior 或 biological truth」 | PASS |
| HTML | 與 artifact snapshot 同一組 2,089／1,624／77.74% 表格與 claim ceiling | 同一組 2,543／949／37.32% 與結論段 | 明示不能宣稱具體樹已穩定確認 | PASS |

注意：HTML 表格視覺格式會把 77.74% 顯示為 77.7% 的欄位格式，但同一 HTML 的 source snapshot、結論與 claim table 均保留 **77.74%**；底層數值仍為 `1624/2089 = 0.777405457...`，不是數據漂移。

詞彙掃描結果：`truth`／`accuracy`／`posterior` 的字面命中是刻意保留的否定警語，例如「不是 posterior 或 truth」「仍不是 biological tree accuracy」，不是正向主張。

## 5. 稽核輸入、命令與輸出片段

### 輸入

- `InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_vaf_pair_metrics.tsv`
- `InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_vaf_pair_summary.json`
- `InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_vaf_pair_checks.tsv`
- `InterSubMod/docs/reports/in_progress/2026/07/20260712_HCC1395兩技術資料粗拓撲與癌症基因藥物一致性驗證_01/artifact.json`
- 同目錄主 Markdown、HTML 與 `qa_receipt.json`

### 執行命令

```bash
python3 - <<'PY'
# Read metrics TSV, summary/checks JSON/TSV, artifact snapshot,
# rendered Markdown/HTML and qa_receipt hashes.
# Recompute every artifact row from integer numerators/denominators;
# verify required semantics and explicit claim negations.
PY
```

### 實際輸出片段

```text
ARTIFACT_SNAPSHOT_VS_METRICS 6/6 PASS
SUMMARY_CHECKS_AND_CONSERVATION 20/20 PASS
KEY_SEMANTICS shape=77.74%; heuristic_exact_forest=37.32% PASS
ARTIFACT_MARKDOWN_HTML SEMANTIC+NEGATION 3/3 PASS
DELIVERED_FILE_HASHES 3/3 PASS
FINAL VAF REPORT AUDIT PASS
```

### 被稽核交付檔 hash

| 檔案 | SHA-256 |
|---|---|
| artifact.json | `be2e021b5616a6e262420e218f132dbb9bfa9784bba54ac8506bc684dc6f9fe8` |
| 主 Markdown | `63ab7045eb31c14aa6f2e7af68d7392762b37ec67c049c057cb11d76443ba90f` |
| 主 HTML | `64811969fe3c63822b02d2f6d7288b436256f32d802c87eb3beb444faf7c2d16` |
| VAF metrics TSV | `19bfc37f4431546eeac99ccdfb014dccecd726f665e414100e7ec8f2acd031a0` |
| VAF summary JSON | `eb9206b23dc6908105d9900db24a5b40e84de2e66cfea769c74607bbe434f1e9` |
| VAF checks TSV | `a4a3bb9f25d8b082ac6b7950c076eedf8700385b242972b194dab106c1d2fa7d` |

## 6. 科學解讀邊界

這項稽核確認「主報告正確使用 VAF heuristic 結果比較兩個 HCC1395 technical datasets」，不確認推測樹為真。read-AF 排序與候選選擇使用同一批 reads，且未完整校正 CN、purity、mutation multiplicity 與 depth；兩個 datasets 也不是兩個獨立生物樣本。因此可以說：

> VAF 支援候選排序，且抽象 shape 在特定可選 subset 有較高技術一致性；具體 mutation-labeled exact forest 的 swap-tolerant agreement 只有 37.32%，仍不能當作 biological tree truth 或方法 accuracy。

**最終：PASS；主報告不需因本次 VAF 稽核修改。**
