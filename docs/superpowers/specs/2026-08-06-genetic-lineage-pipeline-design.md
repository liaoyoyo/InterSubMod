---
title: 純遺傳 lineage 管線設計 — LongLineage → tagged BAM → InterSubMod → HTML
date: 2026-08-06
status: DESIGN — 待使用者 review
task_type: C production
decisions_locked:
  - 路線 C：read-linked component builder（新模組）
  - 拓撲目標：forest + 跨 unit 一致性標記（非單一整樹）
  - tagged BAM：post-freeze 獨立 export 工具
build_branch: research/subclonal-reconstruction-202606
---

# 純遺傳 lineage 管線設計

## 1. 目標與非目標

### 目標
建立一條**不依賴甲基資料**即可運行的 lineage 重建管線，四個模組各自獨立可用，
串接後產出完整可觀察的 HTML；甲基分析退到下游作為**觀察與對齊**，不進入重建。

### 非目標（明確排除，避免重蹈覆轍）
- ❌ **單一全域整樹** — ISM 已判定「跨 unit 整樹 ranking」為 CLOSED NEGATIVE（等機率、不可枚舉）
- ❌ **甲基進入 likelihood 或 gate** — 違反論文既有結論（cis-ASM 循環）
- ❌ **在 LongLineage run 內寫 BAM** — 架構層 const false 禁止（證據見 §4.1）
- ❌ **強制唯一的 read → lineage 標籤** — 模型固有非唯一性，見 §4 的 `ls` 狀態欄

---

## 2. 架構：四模組，各自獨立有意義

```
LongPhase-S tagged BAM（已 haplotag）+ PASS biallelic sSNV VCF + HP/PS sidecar + ref
        │
        ▼
┌───────────────────────────────────────────────────────────┐
│ M1. LongLineage + 新模組 src/linkage/                      │
│     用途（單獨）：純遺傳結構確認與再分群                    │
│     🔑 不需要甲基資料即可運行                               │
│     輸出：frozen artifacts + read_lineage_assignments      │
└───────────────────────────────────────────────────────────┘
        │  read_id = SHA-256(QNAME)
        ▼
┌───────────────────────────────────────────────────────────┐
│ M2. ll-bam-tag（新，獨立工具，post-freeze，run-root 外）    │
│     用途（單獨）：讓 IGV 能直接看到 lineage 結構            │
│     輸出：加標記的 BAM（samtools 可再 sort/index/view）    │
└───────────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────────┐
│ M3. InterSubMod                                            │
│     用途（單獨）：依 BAM 上「有哪些標籤」做甲基對齊分析     │
│     輸出：per-read 甲基 × 標籤統計                          │
└───────────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────────┐
│ M4. Python 整合層                                          │
│     用途（單獨）：依「實際存在的數據」產出觀察 HTML         │
│     🔑 缺件時產出誠實的部分結果，不失敗                     │
└───────────────────────────────────────────────────────────┘
```

---

## 3. M1 — LongLineage 新模組 `src/linkage/`

### 3.1 為什麼是新模組而非改 M1/M2 gate

實測（`docs/reports/20260720_..._parity報告_01.json`，已凍結）：
`79,687 sites → m1_stable 12,851 → m2_eligible 5 → joint_signature_pass 0 → topology_units 0`
甲基 gate 不是篩選而是**全滅**。且拆 gate 需重定義 partner 選擇規則、解除 MM/ML filter、
validator 同步（成本 2×）、與 P3 parity 目標衝突。

新模組**繞開整條 M1/M2 鏈**，輸入沿用既有 `site_reads` artifact（其 `read_id` 不受 gate 影響）。

### 3.2 演算法契約（照 ISM 現行規格）

參考實作：`InterSubMod/scripts/build_strict_ps_hp_regions.py`

| 概念 | 定義 |
|---|---|
| container | `chrom × exact non-missing PS × HP family` |
| node | 有 fixed R/A 的 sSNV |
| edge | ≥T 條 **distinct `read_id`** 兩端皆 fixed R/A |
| component | 上述圖的無向連通分量 |
| 距離 | **只報告，不建邊**（避免幾何假設） |

HP/PS 缺值採 **fail-closed**，規則必須與 ISM 逐字一致，否則數字不可比。

### 3.3 k>12 切分

移植 `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/cpp/exact_ps_k12_partition.cpp`
（39,737 bytes，C++17，**已有完整 Python↔C++ parity 測試套件** 7 檔）。
把 k>12 的 component 切成 k≤12 bounded block，帶 `retained / cut / unavoidable` disposition，
接上 `src/solver/topology_router.cpp:79-82,182-189` 目前直接 abstain 的路徑。

⚠ 需決定是否引入 `boost::multiprecision`（LongLineage 有 no-boost 的 CI 分支）。

### 3.4 solver 不動

`src/solver/` 3,021 行 grep `methyl|m1_|cpg|label` 零命中，
輸入 `ReadPatternObservation` 僅 `pattern_raox + base_qualities + multiplicity`。**保持原樣**。

### 3.5 forest + 一致性層（取代單一整樹）

1. 以 locus 為 key 建 unit 反向索引
2. 對共享 loci 的 unit pair 做 vertex set 投影比對
3. 標記 `COMPATIBLE / CONFLICT / NOT_COMPARABLE`
4. 輸出 **forest**，不合併成單一樹

### 3.6 新增 artifact

**`read_lineage_assignments`**（TSV_BGZF + site_index，沿用既有 artifact 機制）

| 欄位 | 型別 | 說明 |
|---|---|---|
| `dataset_order` | uint32 | |
| `dataset_id` | string | |
| `read_id` | opaque_sha256 | join key |
| `component_id` | string | read-linked 連通分量 |
| `unit_id` | string | topology unit（`.` 若無） |
| `vertex_label` | string | lineage 標籤（`.` 若非唯一） |
| `assignment_status` | enum | `UNIQUE` / `TIE_CLASS` / `PARTIAL` / `ABSTAIN` |
| `tie_class_size` | uint32 | >1 表示非唯一 |
| `abstain_reason` | enum | 沿用 topology_unit 的 reason 詞彙 |

**`topology_unit_membership`**：把已算好但只活在記憶體的 `TopologyMembershipPlan`
（`dataset_producer.cpp:1183-1255`）落檔，讓下游能把 unit 對回 `chrom:pos`。

⚠ 新增 artifact 必須同步：`schema/catalog.json` 的 `run_membership`、`producer_receipt`、
`semantic_digest`、`data_lineage`、`checksums.sha256`，以及
`artifact_validator.cpp:1294-1299, :4802-4805`（否則科學守恆 replay 會**安靜關閉**）。

---

## 4. M2 — `ll-bam-tag` 獨立工具

### 4.1 為什麼必須是獨立工具

| 障礙 | 證據 |
|---|---|
| JSON Schema const false | `production_input_authority.schema.json:69, :108` |
| 4 處 runtime 硬拒 | `longlineage_main.cpp:188`；`compat/regional_io.cpp:391,465`；`regional_compat_validator.cpp:1121` |
| FILE_CENSUS 拒絕多餘檔案 | `artifact_validator.cpp:3091-3101` |
| 禁止 in-place | `artifact_validator.cpp:1255,1614`（`input_snapshot_before == after`） |
| 自我不可重現 | aux tag 屬 `typed_aux` 摘要契約 → 讀回會改變 read identity 與 semantic SHA |
| 讀取路徑不適合 | `block_reader.cpp` 用 `sam_itr_queryi` + 5000bp halo → 跨 block 重複解碼 |

⇒ 獨立工具走 **linear 掃描**（`sam_read1` 全檔循序），不碰 run root。

### 4.2 BAM aux tag 規格

SAM 規格：**含小寫字母的 tag 保留給 local use**（`X?/Y?/Z?` 已被眾多工具佔用，不採用）。

| tag | 型別 | 內容 | 缺值 |
|---|---|---|---|
| `lc` | `Z` | component_id | 不寫該 tag |
| `lu` | `Z` | unit_id | 不寫 |
| `lv` | `Z` | vertex_label | 不寫 |
| `ls` | `A` | `U`=UNIQUE / `M`=TIE_CLASS / `P`=PARTIAL / `A`=ABSTAIN | 不寫 |

**`ls` 是防 overclaim 的機械保證**：稽核指出 read→vertex 在 partial-subcube read、
tie_class>1、family incomplete 三種情況下無唯一解。**任何 `lv` 出現時 `ls` 必須同時出現**，
且 `ls != U` 時 `lv` 僅為代表值，不得被下游當唯一真值。

### 4.3 CLI

```
ll-bam-tag --in-bam <tagged.bam> \
           --assignments <read_lineage_assignments.tsv.bgz> \
           --out-bam <out.bam> \
           [--regions <BED>]        # 限子集，省磁碟
           [--require-all-reads]    # fail-closed：任何 read 查無標籤即失敗
```

join：對每條 alignment 取 QNAME → 剝除 `/1` `/2` 後綴（`ReadParser.cpp:106-112` 的規則）
→ 計算 `SHA-256` → 查 assignments。

### 4.4 磁碟策略

`/big7_disk` 實測僅剩 **617 GB**（99% 已用）。
預設**必須**帶 `--regions`；全量輸出需明確 opt-in 並先確認餘量。

---

## 5. M3 — InterSubMod 甲基分析

### 5.1 新增能力：依 BAM 標籤分組

現況 ISM 只讀 HP tag（`ReadParser.cpp:126`），**PS tag 完全沒讀**（grep 零命中）。
新增：讀 `lc` / `lu` / `lv` / `ls` 四個 tag，作為甲基分析的分組軸。

新增 CLI 選項：
```
--group-by-tag <lc|lu|lv|HP>     # 甲基分析的分組依據，預設 HP（維持相容）
--require-tag-status <U|UM|any>  # 只用 ls=U 的 read / 也接受 M / 全收
```

### 5.2 必須同時修的既有缺陷

| 缺陷 | 證據 | 修法 |
|---|---|---|
| 甲基↔read 靠隱含列序 | `methylation.csv` 首欄是列號；綁定僅在 `ReadAggregator.cpp:63-64` | `methylation.csv` 首欄改為 `read_name`；加 assert |
| cluster id 不落檔 | 只餵給 SignificanceAnalyzer 後丟棄 | 落成 per-read TSV |
| 二值化三套閾值 | binary 0.8/0.2、PerCpgAsm 0.5、Python 128/255 | 統一由 CLI 參數驅動，並在輸出記錄實際使用值 |
| `linkage_matrix.csv` 其實是 TSV | 副檔名與分隔符不符 | 改副檔名為 `.tsv` |
| `significance_summary.csv` 無版本欄 | 源碼 193 欄 / 舊 run 59 欄 / kde 117 欄 | 加 `schema_version` 欄 |

⚠ 這些是 C++ 改動 → 撞 commit Hard Gate（必編譯 + 測試）。

---

## 6. M4 — Python 整合層與降級契約

### 6.1 降級契約（本設計的核心要求）

```
輸入缺件 → 對應面板標記為「不可用 + 原因」，其餘照常產出
```

| 缺什麼 | 影響 | HTML 行為 |
|---|---|---|
| 無 MM/ML 甲基資料 | ISM 甲基分析為空 | 顯示「無甲基資料，僅遺傳觀察」+ 列出被停用的面板 |
| 無 normal BAM | 無 SampleASM / HP_Residual | 該區塊標「需 normal BAM」 |
| topology abstain 比例高 | forest 稀疏 | 顯示 abstain 分佈與原因 census |
| `ls != U` 佔多數 | lineage 標籤非唯一 | **必須**顯示 tie_class 分佈，禁止只呈現代表值 |

### 6.2 實作方式

複用既有資產，不重寫：
- `.claude/skills/verify-workstation/tools/build_workstation.py`（317 行，缺必填值 `exit 3` 拒絕渲染）
- `InterSubMod/scripts/ism_heatmap_std.py`（492 行，甲基 × read 標籤熱圖 SoT，20+ renderer 使用）

新增 **spec builder**：`artifacts → workstation_spec.json`。
這是目前完全空缺的一層（`build_workstation.py` 只有 1 個複用者的根因）。

### 6.3 數字誠信

所有 metric 由 spec 注入，缺 key 即 refuse render（§13-A 由構造防捏造）。
Python **不得**自行計算科學數字，只呈現 artifact 內已凍結的值。

---

## 7. 實作順序

| 階段 | 內容 | 驗收 |
|---|---|---|
| **P1** | `src/linkage/` component builder + `read_lineage_assignments` artifact | 單樣本 HCC1395 跑出 component 數 > 0，與 `build_strict_ps_hp_regions.py` 逐筆 parity |
| **P2** | `exact_ps_k12_partition` 移植 + 接 topology_router | k>12 component 切分結果與既有 Python parity 測試一致 |
| **P3** | forest + 一致性層 + `topology_unit_membership` artifact | 一致性標記的 COMPATIBLE/CONFLICT 計數可重現 |
| **P4** | `ll-bam-tag` 獨立工具 | 輸出 BAM 通過 `samtools quickcheck`；IGV 可見 4 個 tag；`--require-all-reads` 能正確 fail |
| **P5** | ISM 讀 tag + 修 5 項既有缺陷 | 甲基分析可依 `lc/lu/lv` 分組；`methylation.csv` 首欄為 read_name |
| **P6** | spec builder + HTML | 降級契約 4 種缺件情境各有對應 HTML 產出 |
| **P7** | 端到端腳本 + 擴 7 樣本 | 單一命令跑通；7 樣本各自有 receipt |

**P1 完成即可獨立驗證整條路線是否成立。**

P1 go/no-go 初步標準（可調，見 §9.3）：
- **GO**：HCC1395 chr1-22 產出 ≥1,000 個含 ≥2 個 mutation-bearing loci 的 component，
  且與 `build_strict_ps_hp_regions.py` 的 component 邊界逐筆一致
- **NO-GO**：component 數為 0，或 ≥2 loci 的 component 少於 100 個
  → 停止，回頭檢視是 T 閾值、PS/HP 缺值率、或 sSNV 密度的問題

對照基準：ISM 既有 exact-PS 管線在 HCC1395 產出的量級記載於
`InterSubMod/docs/handoff/20260805_系統交接與驗收_01/`，P1 實跑後與其對帳。

---

## 8. 風險

| 風險 | 緩解 |
|---|---|
| 🔴 P1 若 component 數過低 → 整條路線不成立 | P1 設為 go/no-go 閘門，先跑單樣本 |
| 🔴 新增 artifact 漏改 validator → 守恆檢查安靜關閉 | 明列 `artifact_validator.cpp:1294-1299, :4802-4805` 為必改清單，並加測試 |
| 🔴 磁碟 617 GB | `ll-bam-tag` 預設必帶 `--regions` |
| 🟠 ISM C++ 改動撞 Hard Gate | 每步驟必編譯 + `run_tests`（現況 258 tests / 37 suites 全過，須維持） |
| 🟠 boost 依賴 | 移植時先評估能否避開；LongLineage 有 no-boost CI 分支 |
| 🟠 7 樣本資源未知 | HCC1937 BAM 472 GB 完全無實測；P7 前先單樣本測記憶體峰值 |
| 🟠 ADR 治理 | 新模組需寫 ADR 說明與 ADR-0007 claim boundary 的關係 |

---

## 9. 待 review 事項

1. `src/linkage/` 放在 LongLineage 內是否符合你對「模組分開」的期待？（也可獨立 repo）
2. BAM tag 用小寫 `lc/lu/lv/ls` 是否接受？（SAM 規格保留給 local use）
3. P1 的 go/no-go 閘門標準：component 數多少算「成立」？
4. ISM 的 5 項既有缺陷要在本次一併修，還是分開處理？
