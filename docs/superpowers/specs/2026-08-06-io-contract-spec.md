---
title: 純遺傳 lineage 管線 — I/O 契約規範（HCC1395 端到端）
date: 2026-08-06
status: SPEC — 實作依據
companion: 2026-08-06-genetic-lineage-pipeline-design.md
scope: HCC1395 單樣本端到端（P1-P6），7 樣本擴展見 P7
build_branch: research/subclonal-reconstruction-202606
---

# I/O 契約規範

本文定義每個模組的**精確**輸入輸出。所有路徑與欄位均經 2026-08-06 實測確認。
實作時以本文為準；任何偏離必須先改本文。

---

## 0. 實體輸入（HCC1395，已驗證存在）

| 角色 | 路徑 | 實測大小 (bytes) |
|---|---|---|
| **tagged BAM** | `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam` | 278,270,281,041 |
| tagged BAM index | 同上 `.bai` | 117,905,528 |
| PASS sSNV VCF | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/HCC1395.longphase_s.recalibrated.pass.vcf.gz` | 3,728,864 |
| VCF index | 同上 `.csi` | 344,074 |
| HP/PS sidecar | 同目錄 `HCC1395.read_tags.tsv.gz` | 1,433,423,193 |
| sidecar index | 同上 `.tbi` | 910,475 |
| reference | `/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta` | 3,144,230,986 |
| reference index | 同上 `.fai` | — |

### tagged BAM 的 tag 結構（實測 `samtools view chr1:1000000-1010000`）

```
NM:i ms:i AS:i nn:i tp:A cm:i s1:i s2:i de:f rl:i qs:i du:f ns:i ts:i
mx:i ch:i st:Z rn:i fn:Z sm:f sd:f sv:Z dx:i MN:i MM:Z ML:B RG:Z HP:Z PS:i PQ:i
```

關鍵：
- `HP:Z` **字串型**（非整數），值域為 LongPhase-S 九態詞彙（`1`/`2`/`1-1`/`2-1`/`3` 等）
- `PS:i` 整數 phase set
- `MM:Z` + `ML:B` + `MN:i` 甲基三件組齊全
- QNAME = 明文 UUID（如 `e728f803-b1b3-4209-845c-7b0e74e86950`）
- **`lc`/`lu`/`lv`/`ls` 未被任何現有 tag 佔用**（已逐一比對）

---

## 1. M1 — LongLineage `src/linkage/`

### 1.1 輸入 — 直接讀 BAM + VCF（2026-08-06 修正）

**原設計依賴既有 `site_reads.tsv.bgz`，已否決。**
理由：全磁碟搜尋（`find /big7_disk/liaoyoyo2001 -name "site_reads*.bgz"`）**零結果** —
該 artifact 從未落地，取得它需先跑 `dataset-gate`（受甲基 gate 影響、且硬編碼 HCC1395）。
繞一圈反而把依賴引回甲基鏈。

**M1 直接讀原始輸入**（與既有 `build_strict_ps_hp_regions.py` 同源）：

| # | 輸入 | 用途 |
|---|---|---|
| 1 | LongPhase-S tagged BAM + `.bai` | 取 QNAME / FLAG / MAPQ / CIGAR / `HP:Z` / `PS:i` / 各 sSNV 位點的鹼基 |
| 2 | PASS biallelic sSNV VCF + index | 定義 node（sSNV 位點集合） |
| 3 | reference FASTA + `.fai` | 位點鹼基比對 |

read filter：**僅 FLAG + MAPQ ≥ 20**（不要求 MM/ML）
→ 這是「無甲基資料也能跑」降級契約的實現點。

⚠ HP/PS 直接讀 BAM tag（實測 tagged BAM 有 `HP:Z` + `PS:i`），
**不經 sidecar** —— sidecar 是 LongLineage production 路徑為了 truth-isolation 的設計，
M1 作為獨立分群模組不需要那層間接。

### 1.2 演算法參數（manifest 可配置，非硬編碼）

| 參數 | 型別 | 預設 | 說明 |
|---|---|---|---|
| `min_shared_reads` (T) | uint32 | **3** | 建 edge 所需的最少 distinct `read_id` 數。**P1 需做敏感度掃描** |
| `require_exact_ps` | bool | true | PS 必須 non-missing 且精確相等 |
| `require_hp_family` | bool | true | HP family 必須 non-missing |
| `allele_states_for_edge` | enum | `R_AND_A` | 兩端皆須為 fixed R 或 A（`O`/`X` 不算） |
| `max_component_loci` | uint32 | 12 | 超過此值走 §1.5 切分 |

### 1.3 演算法

```
container := (chrom, PS, hp_family)         # PS/HP 任一 missing → 該 read 排除（fail-closed）
node      := sSNV site 在該 container 內有 ≥1 條 fixed R/A 的 read
edge(u,v) := |{ read_id : read 在 u 與 v 皆為 fixed R/A }| >= T
component := 上述無向圖的連通分量
```
**距離只記錄不建邊** — 不引入任何幾何/距離假設。

### 1.4 輸出 artifact A：`linkage_components`

格式 `TSV_BGZF` + site_index，沿用既有 artifact 機制。

| 欄位 | 型別 | 說明 |
|---|---|---|
| `dataset_order` | uint32 | |
| `dataset_id` | string | |
| `component_id` | string | `{chrom}_{ps}_{hp_family}_{ordinal}` |
| `chrom` | string | |
| `phase_set` | uint64 | |
| `hp_family` | string | |
| `n_loci` | uint32 | component 內 sSNV 數 |
| `n_mutation_bearing_loci` | uint32 | 至少有 1 條 ALT read 的 loci 數 |
| `n_reads` | uint32 | distinct read_id 數 |
| `span_start1` / `span_end1` | Position1 | component 基因組跨度 |
| `min_edge_support` | uint32 | 該 component 內最小 edge 的 read 支持數 |
| `partition_disposition` | enum | `INTACT` / `PARTITIONED` / `UNAVOIDABLE_CUT` |

### 1.5 k>12 切分

移植 `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/cpp/exact_ps_k12_partition.cpp`
（39,737 bytes，C++17）。既有 parity 測試 7 檔必須全數移植並通過。

輸出 artifact B：`linkage_blocks`

| 欄位 | 型別 | 說明 |
|---|---|---|
| `component_id` | string | 父 component |
| `block_id` | string | |
| `loci` | string | 逗號分隔的 site_order |
| `k` | uint32 | ≤ `max_component_loci` |
| `disposition` | enum | `RETAINED` / `CUT` / `UNAVOIDABLE` |
| `cut_edges_lost` | uint32 | 切分損失的 edge 數 |

### 1.6 輸出 artifact C：`read_lineage_assignments`

**這是與下游的唯一交界。**

| 欄位 | 型別 | 說明 |
|---|---|---|
| `dataset_order` | uint32 | |
| `dataset_id` | string | |
| `read_id` | opaque_sha256 | = `lowercase SHA-256(QNAME)` |
| `component_id` | string | `.` 若未分派 |
| `block_id` | string | `.` 若無 |
| `unit_id` | string | `.` 若無拓撲解 |
| `vertex_label` | string | `.` 若非唯一 |
| `assignment_status` | enum | `UNIQUE` / `TIE_CLASS` / `PARTIAL` / `ABSTAIN` / `UNASSIGNED` |
| `tie_class_size` | uint32 | >1 表示非唯一；`UNIQUE` 時必為 1 |
| `abstain_reason` | enum | 沿用 `topology_unit` 的 reason 詞彙；`NONE` 若不適用 |

**不變式（必須在 validator 中檢查）**
1. `assignment_status == UNIQUE` ⟺ `tie_class_size == 1` 且 `vertex_label != "."`
2. `vertex_label != "."` ⟹ `assignment_status ∈ {UNIQUE, TIE_CLASS}`
3. 同一 `read_id` 在同一 `component_id` 內只能出現一次
4. `component_id == "."` ⟹ `assignment_status == UNASSIGNED`

### 1.7 輸出 artifact D：`unit_consistency`（forest 一致性層）

| 欄位 | 型別 | 說明 |
|---|---|---|
| `unit_id_a` / `unit_id_b` | string | 共享 ≥1 locus 的 unit 對 |
| `shared_loci` | string | 逗號分隔 |
| `relation` | enum | `COMPATIBLE` / `CONFLICT` / `NOT_COMPARABLE` |
| `reason` | enum | 判定依據 |

⚠ **明確不輸出單一全域樹**（ISM 已判 CLOSED NEGATIVE）。

### 1.8 必須同步的契約（漏改會讓守恆檢查安靜關閉）

- `schema/catalog.json` 的 `run_membership` 四個陣列
- `schema/records/` 新增 4 個 record schema
- `producer_receipt` / `semantic_digest` / `data_lineage` / `checksums.sha256`
- **`src/validation/artifact_validator.cpp:1294-1299` 與 `:4802-4805`**
  （`:4805-4810` 在 catalog 集合不等於硬編碼 9 元集合時 return false，
  `:5224-5230` 只記成 "not claimed" 而 `all_pass` 仍為 true → **fail-silent**）

---

## 2. M2 — `ll-bam-tag`

### 2.1 定位

**獨立 executable，post-freeze，輸出到 run-root 外。** 不是 LongLineage run 的產物。
理由見 design doc §4.1（6 條架構層禁止）。

### 2.2 CLI

```
ll-bam-tag --in-bam <tagged.bam>
           --in-bam-index <tagged.bam.bai>
           --assignments <read_lineage_assignments.tsv.bgz>
           --out-bam <out.bam>
           [--regions <BED>]          # 強烈建議：磁碟僅剩 617 GB
           [--require-all-reads]      # fail-closed
           [--threads N]
```

### 2.3 join 規則

```
QNAME → 剝除 /1 /2 後綴（僅當 paired flag 要求時）→ lowercase SHA-256 → 查 assignments
```
ONT 長 read 為單端，不受後綴規則影響（實測 QNAME 為純 UUID）。

### 2.4 寫入的 aux tag

| tag | SAM 型別 | 來源欄位 | 缺值行為 |
|---|---|---|---|
| `lc` | `Z` | `component_id` | 不寫該 tag |
| `lu` | `Z` | `unit_id` | 不寫 |
| `lv` | `Z` | `vertex_label` | 不寫 |
| `ls` | `A` | `U`/`M`/`P`/`A`（對應 UNIQUE/TIE_CLASS/PARTIAL/ABSTAIN） | 不寫 |

**寫入不變式**
1. `lv` 存在 ⟹ `ls` 必存在（防 overclaim 的機械保證）
2. `ls != 'U'` 時 `lv` 僅為代表值，下游**不得**當唯一真值
3. 既有 tag 一律保留不動（含 `HP`/`PS`/`MM`/`ML`/`RG`）
4. **絕不 in-place 修改輸入 BAM**

### 2.5 實作方式

**linear 掃描**（`sam_read1` 全檔循序 + `bam_aux_append`），
不用 `sam_itr_queryi`（`block_reader.cpp` 的 halo 設計會重複解碼、遺漏遠端 read）。

### 2.6 驗收

- `samtools quickcheck <out.bam>` exit 0
- `samtools view <out.bam> | head` 可見四個 tag
- IGV 載入後可依 tag 分組/著色
- read 總數與輸入完全一致（`samtools view -c`）

---

## 3. M3 — InterSubMod 甲基分析

### 3.1 新增輸入能力

現況 ISM 只讀 `HP`（`ReadParser.cpp:126`），**PS 完全沒讀**（grep 零命中）。

新增讀取：`lc` / `lu` / `lv` / `ls` / `PS`

新增 CLI：
```
--group-by-tag <lc|lu|lv|HP>        # 甲基分析的分組軸，預設 HP（維持向後相容）
--require-tag-status <U|UM|any>     # U=只用唯一解；UM=含 tie；any=全收。預設 U
--emit-per-read-methylation         # 輸出 per-read 甲基 × 標籤長表
```

### 3.2 新增輸出：`per_read_methylation.tsv`

| 欄位 | 型別 | 說明 |
|---|---|---|
| `read_name` | string | **明文 QNAME**（非 hash，ISM 端保持可讀） |
| `read_id` | string | `SHA-256(QNAME)`，供與 LongLineage join |
| `region_key` | string | `chrom:pos:ref:alt` |
| `hp` | string | BAM 的 `HP:Z` |
| `phase_set` | uint64 | BAM 的 `PS:i`，`.` 若無 |
| `lc` / `lu` / `lv` | string | lineage 標籤，`.` 若無 |
| `ls` | char | 狀態，`.` 若無 |
| `allele_call` | enum | `R`/`A`/`O`/`X` |
| `n_cpg_covered` | uint32 | |
| `n_cpg_methylated` | uint32 | 依 `--methyl-high` 判定 |
| `n_cpg_unmethylated` | uint32 | 依 `--methyl-low` 判定 |
| `n_cpg_ambiguous` | uint32 | 落在兩閾值之間 |
| `mean_beta` | float | 平均甲基化程度 |

### 3.3 必須同時修正的既有缺陷

| # | 缺陷 | 證據 | 修法 |
|---|---|---|---|
| 1 | 甲基↔read 靠隱含列序 | `methylation.csv` 首欄是矩陣列號；唯一綁定在 `ReadAggregator.cpp:63-64` | 首欄改為 `read_name`；加 assert 檢查長度一致 |
| 2 | cluster id 完全不落檔 | 只餵 SignificanceAnalyzer 後丟棄 | 併入 `per_read_methylation.tsv` |
| 3 | 二值化三套閾值 | binary 0.8/0.2、PerCpgAsm 0.5、Python 128/255 | 統一由 `--methyl-high/--methyl-low` 驅動；輸出檔記錄實際使用值 |
| 4 | `linkage_matrix.csv` 實為 TSV | 副檔名與分隔符不符 | 改名 `.tsv` |
| 5 | `significance_summary.csv` 無版本欄 | 源碼 193 欄 / 舊 run 59 欄 / kde 117 欄 | 加 `schema_version` 首欄 |

⚠ 全為 C++ 改動 → 每步必編譯 + `run_tests`（現況 **258 tests / 37 suites 全過**，必須維持）。

---

## 4. M4 — Python 整合與 HTML

### 4.1 輸入（全部可缺）

| 輸入 | 缺件時 |
|---|---|
| `linkage_components.tsv.bgz` | 結構面板不可用 |
| `read_lineage_assignments.tsv.bgz` | 標籤面板不可用 |
| `unit_consistency.tsv.bgz` | 一致性面板不可用 |
| `per_read_methylation.tsv` | **甲基面板不可用 → 顯示「無甲基資料，僅遺傳觀察」** |
| `significance_summary.csv` | 統計面板不可用 |

### 4.2 兩階段

```
階段 A: build_lineage_spec.py   artifacts → workstation_spec.json
階段 B: build_workstation.py    spec → standalone HTML
```
階段 B **複用既有** `.claude/skills/verify-workstation/tools/build_workstation.py`（317 行，
缺必填值 `exit 3` 拒絕渲染）。階段 A 是目前完全空缺的一層。

### 4.3 降級契約（本規範的核心要求）

`workstation_spec.json` 必含：
```json
{
  "availability": {
    "<panel_id>": {
      "available": true|false,
      "reason": "<缺哪個檔案 / 哪個欄位>",
      "input_path": "<檢查過的路徑>"
    }
  }
}
```
HTML **必須**渲染不可用面板的佔位區塊與原因，**不得靜默省略**。

### 4.4 必要面板

| panel_id | 內容 | 依賴 |
|---|---|---|
| `component_overview` | component 數 / loci 分佈 / span 分佈 | A |
| `assignment_census` | `assignment_status` 五態分佈、**`tie_class_size` 分佈** | C |
| `consistency_matrix` | COMPATIBLE/CONFLICT/NOT_COMPARABLE 計數 | D |
| `methyl_by_label` | read×CpG β 熱圖，側欄依 `--group-by-tag` | ISM |
| `methyl_availability` | 甲基覆蓋率、缺件說明 | ISM（可缺） |
| `provenance` | 所有輸入檔的路徑 + SHA-256 + 大小 | 必需 |

⚠ `assignment_census` **必須**顯示 tie_class 分佈 —— 只呈現代表值即為 overclaim。

### 4.5 數字誠信

Python **不得**自行計算科學數字，只呈現 artifact 內已凍結的值。
所有 metric 由 spec 注入，缺 key 即 refuse render（§13-A）。

---

## 5. 端到端執行順序（HCC1395）

```bash
# 步驟 0：LongLineage 產出 site_reads（既有能力，dataset-gate 路徑）
#   ⚠ 現況 dataset-gate 硬編碼 HCC1395，正好是本次目標樣本

# 步驟 1：M1 linkage
longlineage-linkage --site-reads <run>/site_reads.tsv.bgz \
                    --min-shared-reads 3 \
                    --out-root <linkage_out>

# 步驟 2：M2 寫回 BAM（限區域，磁碟只剩 617 GB）
ll-bam-tag --in-bam HCC1395_tagged.bam \
           --assignments <linkage_out>/read_lineage_assignments.tsv.bgz \
           --regions <target.bed> \
           --out-bam HCC1395_lineage_tagged.bam
samtools index HCC1395_lineage_tagged.bam

# 步驟 3：M3 ISM 甲基分析
inter_sub_mod -t HCC1395_lineage_tagged.bam \
              -r GRCh38_no_alt_analysis_set.fasta \
              -v HCC1395.longphase_s.recalibrated.pass.vcf.gz \
              --group-by-tag lc --require-tag-status U \
              --emit-per-read-methylation \
              -o <ism_out>

# 步驟 4：M4 HTML
python3 build_lineage_spec.py --linkage <linkage_out> --ism <ism_out> \
                              -o spec.json
python3 build_workstation.py spec.json -o HCC1395_lineage_observation.html
```

---

## 6. 逐步驗證點（每步必過才進下一步）

### P1 對照基準 —— 既有 exact-PS 的實測 canonical funnel

來源（已凍結，`all_pass: true`，8 項守恆檢查全過）：
`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/summary/all7_summary.json`

**HCC1395 實測值（P1 必須合理一致）**

| 指標 | 值 |
|---|---|
| `groups_total` | **11,590** |
| `mutation_bearing_units` | 9,624 |
| `mutation_family_complete_units` | 9,357 |
| `mutation_family_abstain_units` | 267 |
| `objective_certified_units` | 11,323 |
| `objective_abstain_units` | 267 |
| `no_active_alt_units` | 1,966 |
| `ranked_units` | 9,130 |
| `best_tree_unique_units` | 7,047 |
| `best_tree_tied_units` | 2,083 |
| `recurrence_required_units` | 0 |

`active_k` 分佈：`0→1966, 1→3259, 2→4250, 3→1279, 4→405, 5→189, 6→89, 7→49, 8→26, 9→27, 10→13, 11→15, 12→23`

既有 guards：`max_active_bits=12`, `max_family_size=100000`, `max_search_nodes=1000`
（⚠ 與 LongLineage solver 的 default「無上限」不同，跨系統比對 incomplete 率會失真）

**P1 通過標準（取代原「≥1,000」臆測值）**

1. `linkage_components` 的 group 總數與 **11,590** 同量級（容差待 P1 實跑後與使用者共同定）
2. `active_k` 分佈形狀與上表一致（特別是 k=0 的 1,966 與 k=2 的 4,250 兩個峰）
3. 與 `InterSubMod/scripts/build_strict_ps_hp_regions.py` 的 component 邊界逐筆一致
4. `no_active_alt_units`（無 ALT 的純參考 group）比例接近 1,966/11,590 ≈ 17.0%

**NO-GO**：group 總數與 11,590 相差一個數量級以上，或 `active_k` 分佈形狀明顯不同
→ 停止，檢視 T 閾值 / PS-HP 缺值處理 / sSNV 輸入是否同源

---

| 步 | 驗證 | 通過標準 |
|---|---|---|
| P1 | `linkage_components` 產出 | 見上表四項對照 |
| P1b | T 敏感度 | T ∈ {2,3,4,5} 各跑一次，component 數變化曲線合理（非斷崖） |
| P2 | k>12 切分 | 7 檔 parity 測試全過；`cut_edges_lost` 總和可重現 |
| P3 | forest + 一致性 | 四個 artifact 的守恆等式成立；validator 13/13 PASS |
| P4 | `ll-bam-tag` | `samtools quickcheck` exit 0；read 數與輸入一致；IGV 可見四 tag |
| P5 | ISM | `run_tests` 258/258 維持；`per_read_methylation.tsv` 的 read_id 能與 assignments 100% join |
| P6 | HTML | 4 種缺件情境各產出一份 HTML，缺件面板顯示原因 |

**P1 是 go/no-go 閘門**：若 ≥2 loci 的 component < 100 個，停止並重新評估 T / PS-HP 缺值率 / sSNV 密度。

---

## 7. 未決事項（實作前需確認）

1. `src/linkage/` 置於 LongLineage 內，或獨立 repo？
2. `min_shared_reads` 預設 3 是否合理？（P1b 會給實測依據）
3. `boost::multiprecision` 是否引入（LongLineage 有 no-boost CI 分支）
4. `--regions` 的目標 BED 用哪些區域？（建議先用已驗證的 chr2:18M 等位點）
