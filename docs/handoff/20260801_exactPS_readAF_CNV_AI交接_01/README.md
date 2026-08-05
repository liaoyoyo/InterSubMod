<!--
建立時間: 2026-08-01
狀態: validated
目標: 提供 exact-PS×HP 區域 mutation-state tree 研究的單一 AI／研究者交接入口
處理範圍: 7 technical datasets / 6 biological IDs / chr1-22 / exact PS×primary HP / read-AF / topology / methyl auxiliary
關聯檔案:
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/20260801_研究狀態與AI交接手冊_01.md
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/20260801_readAF_CNV與單一樹決策規格_01.md
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
驗證命令:
  - python3 -m json.tool docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
  - python3 -m json.tool docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.schema.json
-->

# Exact-PS×HP 區域樹研究：AI／研究者交接入口

> **30 秒結論**：目前已完成的是 recurrence-allowed、read-compatible、
> minimum mutation-state candidate-tree 的全七 technical datasets 技術分析。
> read-AF 只提供未校正 CN/LOH 的條件式排序；不能稱為真實 cellular clone tree。
> 10,717 個 mutation-bearing units 因 search guard 明確 abstain。

## 1. 使用順序

新研究者或新 AI 應依序讀取：

1. [`authority_manifest.json`](authority_manifest.json)：primary machine-readable authority index。
2. [研究狀態與 AI 交接手冊](20260801_研究狀態與AI交接手冊_01.md)：背景、流程、數據與目前工作。
3. [read-AF、CNV 與單一樹決策規格](20260801_readAF_CNV與單一樹決策規格_01.md)：何時可輸出一棵樹、何時只能 representative 或 abstain。
4. [權威資料與過期來源登錄](20260801_權威資料與過期來源登錄_01.md)：可引用與不可回流的來源。
5. [`denominator_registry.tsv`](denominator_registry.tsv)：論文表格與百分比的固定分母。
6. [AI 新讀者測試與修正紀錄](20260801_AI新讀者測試與修正紀錄_01.md)：機器 QA 與無前情 reader test。
7. [全 7 資料集 HTML observation report](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report/all7_v1/report.standalone.html)：
   最終數據觀察、樣本表、SVG 圖、CN/LOH／methyl 狀態與 provenance。
8. [HTML builder 與重生說明](../../../research/20260801_exact_ps_observation_report/00_INDEX.md)：
   Python finalizer、輸入／輸出、negative tests、responsive／no-JS／print QA。

不得只讀舊論文、舊工作站或 `CURRENT_FOCUS.md` 的單一段落後自行推論最新狀態。

## 1.1 最新衍生觀察交付

正式 release：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
└── 20260801_exact_ps_observation_report/all7_v1/
```

狀態：

- report build：`PASS`
- browser QA：`PASS`
- release：`VALIDATED_DERIVED_OBSERVATION`
- scientific completeness：`INCOMPLETE_WITH_ABSTAIN`

HTML 是由本目錄的 authority manifest 與 denominator registry 重生的 derived view，
不取代 canonical artifacts。Release receipt 會把 builder、QA、report data 與 HTML
的 SHA-256 綁定；CN/LOH 仍為 `NOT_INTEGRATED`，甲基仍為 `association-only`。

## 2. 目前最重要的決策

### 2.1 read-AF 不是錯，但適用範圍比原先敘述窄

最終 clone／cell population size 會影響 read-AF 的差距與可辨識性；在理想
copy-neutral nested lineage 中，它不會破壞「祖先 carrier fraction 不低於後代」
的必要條件。但 read-AF 是終點 molecule snapshot，不是突變時間、CCF、likelihood
或 posterior。

### 2.2 CNV／LOH 未整合，不能把 AF winner 當作生物單一樹

allele-specific amplification、copy-neutral LOH、subclonal CNA、mutation-bearing
copy loss 與 multiplicity 都能使 molecule AF 與 cell prevalence 脫鉤。
exact HP 分組只能降低 homolog mixing，不能排除同一 HP 衍生 duplicated copies。

### 2.3 多候選可以顯示一棵 representative，但不能冒充唯一答案

- 完整候選集合本來只有一棵：可稱 model-conditional exact single candidate。
- AF 算術上唯一，但 CN/LOH 與 uncertainty 未通過：只能稱
  `RAW_AF_UNIQUE_REPRESENTATIVE`。
- 多棵但同一 rooted-unlabeled topology：顯示 topology representative，
  不宣稱 mutation-labeled tree 唯一。
- 跨 topology 並列：只能顯示 consensus backbone 或完整 candidate set。
- family 未列完：必須 abstain。

## 3. Canonical 技術結果

```text
98,955 final groups
├─ 13,014 no active ALT
└─ 85,941 mutation-bearing
   ├─ 75,224 complete minimum families
   │  ├─ 71,955 read-AF ranked
   │  │  ├─ 39,648 unique best tree
   │  │  ├─ 23,858 tied, same rooted-unlabeled topology signature
   │  │  └─  8,449 tied, cross topology
   │  ├─  3,224 zero denominator
   │  └─     45 AF_ABSTAIN_RECURRENCE_SCREEN
   └─ 10,717 family_status=ABSTAIN_RESOURCE_LIMIT
              → tree_decision=ABSTAIN_INCOMPLETE
              → abstain_reason=SEARCH_NODE_GUARD
```

在 71,955 個 ranked units 中，63,506（88.2579%）只有一種
root-preserving unlabeled topology。這是 ranked-only、model-conditional
graph-shape 統計，不是所有區域的 biological topology prevalence。

## 4. 專案狀態

| 層級 | 狀態 | 可主張 |
|---|---|---|
| 格式／receipt／hash | Technical PASS | 此 frozen output 可重現與守恆 |
| Minimum candidate family | 部分完成 | 75,224 complete；10,717 abstain |
| read-AF | 算式與 oracle parity 完成 | 未校正 CN/LOH 的 deterministic ordering |
| Exact topology census | ranked subset 完成 | AF-optimal mathematical tree shapes |
| Cellular clone lineage | 未驗證 | 不可主張 |
| CN/LOH-aware CCF | 未實作 | 不可主張 |
| Methyl auxiliary | association-only validated | 不可改 edge 或證明 ancestry |
| 論文 | 7/12 骨架、內容過期 | 不可直接送審 |

## 5. 新 AI 的固定啟動提示

可將下列提示與本目錄一起交給新 AI：

```text
先讀 authority_manifest.json，逐一確認 required artifact path 與 SHA-256。
本專案的分析單位是 local exact-PS×HP block，不是 cell。
Primary mutation model 是 recurrence-allowed upward Hamming-1 arborescence，
不是 strict perfect phylogeny。
read-AF 是 supported-pattern ALT/(REF+ALT) 的條件式排序，
不是 likelihood、posterior、CCF 或 CN/LOH-corrected prevalence。
任何 family incomplete、CN/LOH unknown、AF unstable 或 topology tie，
均依 tree_decision_policy 保留 representative／candidate set／ABSTAIN，
不可強迫選一棵 biological tree。
回答任何百分比前，先從 denominator_registry.tsv 指出 numerator、
denominator、excluded/unknown 與 authority。
```

## 6. 更新規則

1. 新結果不得直接改舊數字，先新增 artifact 與 SHA-256。
2. 只有通過 schema、receipt、denominator conservation 與 claim audit，才可更新
   `authority_manifest.json`。
3. `authority_manifest.json` 更新時同步更新本 README、分母登錄與過期來源登錄。
4. 不能以 Git HEAD 取代 artifact hash；目前 working tree 不是 clean release。
5. 不複製大型 BAM／JSONL 到本目錄，只保存 identity、path、scope 與 checksum。
6. `denominator_registry.tsv` 是 numerator／denominator／excluded 的機器表；
   manifest 是其 authority index，不以兩份文件建立互相競爭的事實來源。
