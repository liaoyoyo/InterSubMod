---
title: LongPhase-S → tagged BAM → 甲基分析 → 觀察 HTML 的介面契約
date: 2026-08-08
status: 實測驗證（HCC1395 全基因組）
scope: chr1–22；單樣本 HCC1395；其餘 6 樣本僅拓撲層可用
data_sources: /bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/receipt.json
---

# 介面契約

三個獨立專案串成一條鏈。每一段的輸入／輸出都是檔案，可各自單獨執行、單獨驗證。

```
LongPhase-S tagged BAM ─┐
somatic VCF ────────────┼─→ ① exact-PS partition ─→ topology.jsonl + MLHP + endpoint_edges
reference FASTA ────────┘                              │
                                                       ├─→ ② LongLineage tag_bam ─→ lineage_tagged BAM
                                                       │                              │
                                                       │                              ├─→ ③ ISM ─→ 逐位點甲基產物
                                                       │                              │
                                                       └──────────────────────────────┴─→ ④ drilldown ─→ HTML
```

---

## ① exact-PS partition（上游，本文件不涵蓋實作）

| | |
|---|---|
| 輸入 | LongPhase-S tagged BAM、somatic VCF、reference FASTA |
| 輸出根 | `<ROUND>/samples/<SAMPLE>/` |
| 關鍵產物 | `<S>.topology.jsonl`（**硬核心**）、`<S>.exact_ps_mlhp.json`、`chromosomes/<chr>/strict_regions/*.endpoint_edges.tsv.gz` |

**連鎖規則**（`region_linkage_rule`）：`strict fixed-R/A endpoint-pair support connected component`，門檻 `MINREAD=3`、`MAX_SNV=12`。

⚠ **三份產物的位點集合不同，混用會對錯位**：

| 檔案 | 欄位 | 內容 |
|---|---|---|
| `topology.jsonl` | `active_positions` | **只有 active 位點**（該 HP 家族內有 ALT 的） |
| `exact_ps_mlhp.json` | `positions` | **全部 sSNV**（含非 active） |
| — | `active_original_indices` | 前者在後者中的索引 |

實測 **60.8%（7,049/11,590）的 region 兩者長度不同**。`pattern_vector` / `populations_by_hp` 的字元數對應的是 **MLHP 的完整 positions**，不是 `active_positions`。

---

## ② LongLineage tag_bam → lineage_tagged BAM

| | |
|---|---|
| 輸入 | 原始 BAM（含 `MM`/`ML`）、sidecar（HP/PS）、`topology.jsonl`、`unit_lineage_paths.tsv.gz`、`read_lineage_assignments.tsv.gz` |
| 輸出 | `<S>.<chr>.lineage_tagged.bam` ×22 + `.bai` + `.tag_bam.receipt.json` |
| 實測 | 259.6 GB / 22 檔；39,617,373 alignment；1,251,092 帶 lineage 標籤 |

### BAM aux tag 契約

| tag | 型別 | 內容 | 不變式 |
|---|---|---|---|
| `HP` | **`Z`（字串）** | 九態：`1` `2` `1-1` `2-1` `3` | 由 sidecar 注入 |
| `PS` | `i` | phase set | 同上 |
| `lc` | `Z` | lineage component（unit_id） | — |
| `lu` | `Z` | block id | — |
| `lv` | `Z` | 階層路徑 `HP2-1-1`；`+` 後綴 = **此節點或其後代** | **存在必有 `ls`** |
| `lp` | `Z` | 觀察 pattern（R/A/X） | 事實，非推論 |
| `lo` | `Z` | 突變順序（`pos>pos>pos`） | **`ls=A` 時不寫** |
| `ls` | `A` | `U` 唯一 / `M` 同分並列 / `P` 部分覆蓋 / `A` 拓撲未定 | `ls≠U` 時 `lv` 僅為代表值 |
| `lg` | `Z` | **合併觀察標籤** = `lv` 去 `+`；跨 block 取最深那段 | 與 `lv` 同時寫；實測 13 個相異值 |
| `ln` | `i` | **相容頂點數**（LCA 候選集大小；命中頂點為 1；跨 block 取 max） | 與 `lv` 同時寫；`ls=U` ⟹ `ln=1` |

**原有 `MM`/`ML` 甲基標記完整保留** —— 同一份 BAM 既可看譜系也可做甲基分析。

### `lg` 與 `ln` 為什麼要存在（2026-08-13 新增）

**`lg`**：`lv` 把 `HP2-1` 與 `HP2-1+` 分成兩個字串，而 IGV 的 color-by-tag
**按字串精確配色**，於是「同一支譜系，不論確定或推論」在畫面上永遠是兩色，看不出結構。
`lg` 把兩者併成同一值。不照抄逗號串的理由：多 block 的 `lv` 形如 `HP2+,HP2-1`，
直接配色會產生一堆只出現一兩次的組合值（chr21 實測原值 **30 種**）；
取最深那段 = 「這條分子最確定能到達的位置」，值收斂到 **13 種**，IGV 才可讀。

```
Color alignments by → tag → lg     ← 譜系（HP2-1 與 HP2-1+ 同色）
Group alignments by → tag → ls     ← 可信度（U / M / P / A 分開）
```

**`ln`**：原本 `compatible_labels()` 的候選數算完就丟，只進 receipt 的
`lca_candidates_sum` 總和 —— 相容 2 個頂點與相容 16 個頂點的 read 在 BAM 裡完全一樣。

🔑 **chr21 實測揭出一件重要的事**：

| `ls` | `ln==1` | `ln>1` |
|---|---:|---:|
| `U` | 614 | 0 |
| `M` | 491 | 0 |
| `P` | **1,840** | 1,866 |

`ls=P` 中有 **49.6% 是 `ln==1`** —— 帶著 `+`，但相容集只有一個頂點，
表示後代都因某位點讀到 `R` 而被排除，**read 其實就在該節點上**，`+` 保守到過頭。

⇒ 下游若要「已確認」的集合，用 **`ls=U OR ln==1`** 比只用 `ls=U` 好：
chr21 從 614 條擴大到 2,454 條（**4×**），而且是把原本丟掉的資訊撿回來，**不是放寬標準**。
全基因組 `ls_U` 僅佔有標籤 read 的 7.7%，這個擴大對樣本量很關鍵。

### 🔴 IGV 可載入性（實測結論）

**BAM 本身完全合法，沒有缺任何資訊**：

- `samtools quickcheck` ✓
- `@HD SO:coordinate` ✓（座標排序）
- `.bai` 索引存在 ✓
- `@SQ` 195 條，22 個檔的 header **完全一致**（可 `samtools merge`）
- 每個檔只含自己那條染色體的 read

**不能「直接」進 IGV 的兩個原因，都不是資訊缺失**：

1. **檔案被切成 22 個**（每染色體一個，合計 259.6 GB）。
   IGV 要嘛載 22 條 track，要嘛先合併：
   ```bash
   samtools merge -@ 16 -o HCC1395.lineage_tagged.bam bam/HCC1395.chr*.lineage_tagged.bam
   samtools index -@ 16 HCC1395.lineage_tagged.bam
   ```
   （合併後仍是 259.6 GB，確認磁碟餘量）

2. **`HP` 是字串型別 `HP:Z`，不是 IGV 內建 haplotype 分組認的 `HP:i` 整數。**
   這是刻意的 —— 九態 HP（`1-1` / `2-1` / `3`）無法用整數表達。
   IGV 仍可分組，只是要走「**Group alignments by → tag → 輸入 `HP`**」，
   而不是內建的 haplotype 模式。同樣方式可用 `lv`（lineage 階層）、`ls`（可信度）分組。

**IGV 建議操作**：
- Color/Group by tag `lv` → 看 read 的譜系歸屬
- Group by tag `ls` → 分開唯一解與部分覆蓋
- 開啟 base modification 顯示 → 直接看 `MM`/`ML` 甲基

---

## ③ ISM（inter_sub_mod）→ 逐位點甲基產物

| | |
|---|---|
| 輸入 | lineage_tagged BAM、normal BAM、somatic VCF、reference FASTA |
| 輸出 | `<out>/<chrom>/<chrom>/<chrom>/<chrom>_<pos>/<chrom>_<start>_<end>/` |
| 實測 | 29,754 位點 / 26.8 分鐘（24 threads） |

### 產物

```
metadata.txt
methylation/methylation.csv        read × CpG 的 β；第一欄 read_id，其餘欄名即 CpG 座標
distance/<METRIC>/matrix.csv       read × read 距離
clustering/leaf_order.txt          UPGMA 葉序（存 read **name**）
clustering/linkage_matrix.csv      cluster_i, cluster_j, distance, new_cluster_id, size
clustering/significance.json       含 optimal_k
reads/reads.tsv                    read_id, read_name, chr, start, end, mapq, hp,
                                   alt_support, is_tumor, strand + 7 個 lineage 欄
<out>/<chrom>/significance_summary.csv   逐位點 199 欄的多軸統計
```

### 🔴 三個必須知道的 ID 空間差異

| 檔案 | 用哪個 ID | 後果 |
|---|---|---|
| `methylation.csv` / `distance/matrix.csv` | **read_id（窗內整數索引）** | 跨窗不通用 |
| `clustering/leaf_order.txt` | **read_name（UUID）** | 與上者交集為 **0**，必須經 `reads.tsv` 轉換 |
| `read_lineage_assignments.tsv.gz` | **`qname_sha256`** | `= sha256(read_name)`；且 **同一 read 會落多個 block，join 必須用 `(qname_sha256, region_id)` 複合鍵** |

### 分組軸

`--group-by-tag HP,ALT,lc,lu,lv` + `--require-tag-status {U|UM|any}`

**預設用 `UM`**（實測 chr21 370 位點）：

| gate | 可檢定 | 顯著率 | 語意 |
|---|---:|---:|---|
| `U` | 36（9.7%） | 66.7% | 只收唯一解 |
| **`UM`** | **127（34.3%）** | **56.7%** | U + M，**都是點斷言** |
| `any` | 152（41.1%） | **90.8%** | 含 `+` 子樹斷言 |

`any` 的 90.8% 是 artifact：`+` 只出現在 `ls=P`（實測 2,092/2,094，U 與 M 皆為 0），而 `+` read 中位跨度 **19,125 bp vs 非 `+` 的 28,394 bp（1.48×）** —— 用含 `+` 的標籤分群等同部分依 read 長度分群，與甲基構成 confound。

---

## ④ drilldown → 觀察 HTML

```bash
python3 InterSubMod/scripts/build_drilldown_dashboard.py \
    --sample HCC1395 --out <DIR> \
    --bake-panels all \
    --lineage-assign <...>/read_lineage_assignments.tsv.gz \
    --lineage-paths  <...>/unit_lineage_paths.tsv.gz
```

### 能力探測（capability detection）

**硬核心只有 `topology.jsonl`**（缺則 `exit 3`）。其餘六層各自獨立探測，缺了就把對應面板降級並在頁面上寫明缺什麼、探測停在哪一段 —— **不會整份拒繪**。

| 能力 | 缺了會怎樣 |
|---|---|
| `topology`（CORE） | REFUSE，`exit 3` |
| `mlhp` | read state matrix 不可用 |
| `ism_dirs` | 甲基面板、軸軌道、IGV 圖全部降級 |
| `strict_edges` | 邊的 read 支持與門檻 what-if 不可用；自檢 C12 變「無法檢查」 |
| `lca_ab` | LCA 增益 KPI 與自檢 C11 不可用 |
| `annotations` | 無 drop-in 篩選維度 |

### 註釋 drop-in

把檔案丟進 `<out>/annotations/`，**零程式碼改動**就變成篩選維度：

| 格式 | 需要的欄位 |
|---|---|
| `.bed` | `chrom start end [name]` |
| **SAVANA / CN segment `.tsv`** | `chromosome start end copyNumber [minorAlleleCopyNumber]` → 自動分 gain(CN≥3) / loss(CN≤1) / loh(minor<0.5) / neutral，**門檻與 `scripts/analysis/savana_to_smcnbed.py` 一致** |
| 位點表 `.tsv`/`.csv` | 含 chrom + pos（容忍多種欄名寫法） |

檔名去副檔名即維度名。解析失敗不靜默略過 —— 能力矩陣逐檔說明讀到什麼、命中幾個、或為什麼解析不了。

### 輸出結構（可整包搬走）

```
index.html                  826 KB   shell + L1 columnar（19,849 sSNV）
data/L{2,4}.<chr>.js         32 MB   region + 樹 + ISM 統計（22 片，最大 chr7 2.78 MB）
data/L5.<chr>.js                     甲基面板 manifest
panels/<chr>/<pos>[.T].png  133 MB   16,302 位點 × 全部/tumor-only 兩版
igv/<region_id>.js                   IGV 式 read 對齊圖（點擊才載入）
annotations/                         drop-in 資料夾
receipt.json                         所有輸入的 sha256 + 能力矩陣
SELFCHECK.md                         12 條守恆等式
```

⚠ **分片用 `<script src>` 注入而非 `fetch`** —— Chrome 對 `file://` 的 fetch 一律以 origin `null` 阻擋，但 subresource script 可正常載入。載入失敗顯示錯誤卡片而非留白。

### 自檢（12 條守恆等式，實測 10 通過 / 1 不成立 / 1 無法檢查）

唯一的 ❌ 是 **C10 循環論證**：「甲基自身分群」軸 20,903/20,904 = 100.0% 顯著 —— 分組標籤（cluster_labels）與檢定用的距離矩陣同源（`SignificanceAnalyzer.cpp:123`）。該軸已在三處隔離：表格標紅籤、結論句排除、軸軌道不含它。

---

## 已知限制（對外使用前必讀）

1. **單樣本**。7 個樣本都有 topology，但 lineage tag 與 ISM 只有 HCC1395。比較區塊會標「從缺」。
2. **ISM 覆蓋 81.6%**（16,195/19,849 sSNV）。缺的**不是隨機** —— 是被 ISM 自己的 coverage/CpG gate 濾掉，與 CpG 密度、覆蓋深度相關。任何「拓撲 × 甲基」的結論都必須聲明母體。
3. **ISM 窗只有 ±1 kb**。寬 region 的甲基是斷續的（實測 chr4 span 10,140 只有 53.0% 有 CpG）。
4. **k≥7 的 region 沒有樹**。153 個全是 `family_incomplete` / `ABSTAIN_RESOURCE_LIMIT` / `SEARCH_NODE_LIMIT_REACHED`，solver 主動放棄。那些 region 的 IGV 圖只能用等位 pattern 分組。
5. **lineage vertex ≠ subclone**。節點由 read 共現的簡約解推得；稱它為 subclone 需要 single-cell 或 multi-region 的獨立證據，本管線沒有。
6. **read 比例 ≠ CCF**。本樣本共現區 94% 為 CN-altered。
7. **甲基不參與拓撲推論**，只在事後 characterize。

---

## 每一段都可單獨驗證

| 段 | 驗證方式 |
|---|---|
| ① | `strict_regions/*.endpoint_edges.tsv.gz` 的 `passes_primary_threshold` 重算連通性 → 實測 11,590/11,590 = 100% |
| ② | `bam/*.tag_bam.receipt.json` 的 `ls_P+ls_M+ls_A+ls_U == lineage_written` → 實測 1,251,092 分毫不差 |
| ③ | `run_summary.json` 的 `regions.succeeded` vs `significance_summary.csv` 列數 |
| ④ | `--probe-only` 印能力矩陣；`SELFCHECK.md` 的 12 條等式；`receipt.json` 的輸入 sha256 |
