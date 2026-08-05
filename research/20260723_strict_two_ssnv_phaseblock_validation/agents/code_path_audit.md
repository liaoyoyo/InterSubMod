<!--
建立時間: 2026-07-23
目標: 完整稽核最新版 strict endpoint region 從 sSNV/molecule 到 block/tree input 的實際程式語意
處理範圍: strict graph、region builder、production runner、k<=12 partition、tree adapter、solver 入口、HCC1395 chr1-22 receipts
服務目標: G1、G4、G5
任務類型: B（Comprehensive validation）
關聯檔案:
  - InterSubMod/tools/strict_endpoint_graph.py
  - InterSubMod/scripts/build_strict_ps_hp_regions.py
  - InterSubMod/scripts/run_layered_v4_strict.py
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/extract_lossless_read_linkage_collapsing.py
  - InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/scripts/exact_ps_k12_partition.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/exact_ps_partition_to_mlhp.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/layered_tree_reconstruction.py
-->

# Strict two-sSNV phase-block 程式路徑稽核

> **框架：輸入契約 → edge → container → component → bounded block → tree input。**
>
> **TL;DR**：strict source-region 建構的核心語意正確：同一 canonical molecule 必須在兩個不同 sSNV 各有 fixed `R/A` 才能對該 pair 投票；edge 不跨 chromosome、exact nonmissing PS 或 primary HP family；A–B 與 B–C 可依連通性把 A/B/C 放在同一個 source region，但不會虛構 A–C direct edge，也不能據此宣稱同一 molecule、cell 或 clone。需要修正的是下游用詞與防呆：graph registry 含 `k=1` singleton，只有 `k>1` 才是 tree-eligible region；`k>12` partition 後的 block 並未強制誘導子圖連通；adapter 也會從全部原始 molecule 重新投影，而非只使用 `retained` constraints。HCC1395 實測在 adapter gate 後所有 11,590 個 tree-input blocks 仍 strict-connected，但其中 75 個 block 的 partition retained weight 為 0，另有 383 個只剩單一 fixed-site 的 MINREAD-supported patterns，必須明示為計算輸入而非共同 clone 證據。

## 1. 五項問題的直接結論

| 問題 | 結論 | 信心 | 需要保留的限制 |
|---|---|---:|---|
| 1. 同一 canonical molecule 是否必須對至少兩個不同 sSNV 各有 fixed `R/A` 才投 edge？ | **是。** 少於兩個 fixed calls 時可成為 active singleton node，但不產生 edge vote。 | 高 | `support_total` 包含 `RR/RA/AR/AA`，不是只算 ALT-bearing states。 |
| 2. 每個 edge 的兩端是否確實是 sSNV？ | **本次 production 資料中是；reusable producer/builder 尚未完整 fail-closed 驗證。** | 高 | Extractor 只強制 PASS、biallelic、single-base；本輪 catalog audit 另證明 `A/C/G/T` 且 REF≠ALT。graph core 把 `site_index` 視為 opaque node。 |
| 3. container 是否嚴格 chromosome × exact nonmissing PS × primary HP？ | **是，但 HP 是 primary HP family 1/2，不是原始 HP subtype。** | 高 | upstream 會把 `1/1-1/1-2` 合為 family 1，把 `2/2-1/2-2` 合為 family 2。 |
| 4. region 是否 `k>=2 connected component`？ | **要分兩層。** registry 的 component 包含 `k=1` singleton；tree-eligible source region 才是 `k>1` connected component。 | 高 | 不可把 39,846 個 all components 全稱為可建樹 regions；HCC1395 tree denominator 是 11,462。 |
| 5. 轉成 block/tree input 後有無語意偏移？ | **有受控轉換，也有兩處需補契約。** source `component_id`、PS、HP 都保留，但 `k>12` 會變成多個 bounded blocks；adapter 重新投影所有 molecule 並做 exact-vector `MINREAD`，不是直接傳遞 source edges/retained constraints。 | 高 | block 不是新的 biological region；tree 是 regional mutation-state candidate，不是 confirmed clone tree。 |

## 2. Step → Verify 稽核鏈

1. 追溯 sSNV catalog 來源  
   → 驗證：VCF loader 只保留單一 ALT、`len(REF)=len(ALT)=1`、`FILTER=PASS`。
2. 追溯 molecule→edge 投票  
   → 驗證：只有 fixed `R/A` 進入 pair combinations；molecule ID 與 site index 唯一。
3. 追溯 exact container  
   → 驗證：chrom scope fail-closed，container key 為 `(HP family, raw exact PS)`，不同 PS/HP 不共用 accumulator。
4. 追溯 component 與 tree gate  
   → 驗證：threshold-qualified edges 才 union；multisite component 重新做 connectivity assertion；singleton 永不進 tree。
5. 追溯 partition/adapter/solver  
   → 驗證：source identity preserved、`k_B<=12`、cross-PS/HP=0；盤點 block connectivity、pattern projection 與 claim ceiling。

## 3. 問題 1：同一 molecule 何時能投 edge？

### 3.1 程式語意

`InterSubMod/tools/strict_endpoint_graph.py`：

- `:31`：fixed call vocabulary 僅為 `{"R","A"}`。
- `:43-76`：`MoleculeCalls` 要求 site index 嚴格遞增且唯一，call code 只能是 `R/A`。
- `:162-171`：duplicate `molecule_id` fail closed；同一 molecule 的 fixed calls 經 `combinations(..., 2)` 產生所有 pair votes。
- `:94-111`：pair endpoints 必須 `site_i < site_j`，且 `total = RR+RA+AR+AA > 0`。

`InterSubMod/scripts/build_strict_ps_hp_regions.py`：

- `:279-291`：每列 site indices 必須排序且唯一，並逐一綁定 site catalog position。
- `:301-311`：先把 `X/L/O/D/S` 全部排除，只把 fixed `R/A` 傳給 `EndpointPairAccumulator`。

因此：

```text
read r1: A=R, B=A     → 對 edge A—B 投 1 票（RA）
read r2: A=R          → A 是 active node，但不產生 edge
read r3: A=R, B=X     → A 是 active node，B 不因 X 建 edge
```

一條 molecule 有三個 fixed calls A/B/C 時，會各投 A–B、A–C、B–C 一票；同一 molecule 對同一 pair 至多一票。

### 3.2 A–B、B–C 串聯

`InterSubMod/tests/test_strict_endpoint_graph.py:26-33` 明確測試：

```text
3 molecules support A—B
3 different molecules support B—C
0 molecules support A—C
```

結果是 component `(A,B,C)`，但 edge table **沒有** `(A,C)`。

所以 A–B/B–C chaining 對「connected analysis domain」是標準 graph transitivity；它只證明存在一條 evidence path：

```text
A ↔ B ↔ C
```

不證明：

- 有同一條 read 同時觀察 A/B/C；
- A 與 C direct co-observed；
- 三者在同一 tumor cell；
- 三者一定屬於同一 clone；
- A→B→C 是唯一時間順序。

## 4. 問題 2：edge endpoints 是否真的都是 sSNV？

### 4.1 本次 production data：是

`InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/extract_lossless_read_linkage_collapsing.py`：

- `:244-257`：VCF loader 只保留 `len(REF)=1`、exactly one ALT、`len(ALT)=1`；非 PASS record 直接報錯。
- `:258-262`：依 position 排序並拒絕同 chromosome 重複 positional sSNV。
- `:707-720`：site catalog 直接由上述 `Variant(ref, alt)` 寫出。

`InterSubMod/scripts/run_layered_v4_strict.py`：

- `:55-60`：production default extractor 與 strict region builder 綁定為上述版本。
- `:916-927`：每條 chromosome 先 materialize extraction，再取得同一 extraction 的 site catalog/molecule calls。
- `:929-949`：strict builder 直接使用該 catalog，完成後驗 receipt。

HCC1395 production provenance 也綁定：

```text
extractor sha256 = 2ca7ccb67c89e816fae9284f4e2ba21b186378105086c6b0128ed5445a133490
tree VCF = HCC1395.longphase_s.recalibrated.pass.vcf.gz
```

本次獨立重算 HCC1395 chr1-22：

| 檢查 | 結果 |
|---|---:|
| site catalog rows | 79,687 |
| 非單鹼基 REF/ALT | 0 |
| REF/ALT 非 A/C/G/T | 0 |
| REF=ALT | 0 |
| edge endpoint 不在同 chromosome catalog | 0 |
| duplicate container-edge key | 0 |

### 4.2 Standalone 契約缺口

`InterSubMod/scripts/build_strict_ps_hp_regions.py:217-244` 的 `load_site_catalog()` 只驗：

- header 有 `ref/alt` 欄；
- chromosome、site index、position、locus uniqueness；
- position 嚴格遞增。

它**沒有使用** `ref/alt` 值，也沒有驗：

- `len(ref)=len(alt)=1`；
- `ref/alt ∈ {A,C,G,T}`；
- `ref != alt`；
- catalog 是否確實來自指定 somatic PASS VCF。

`InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/scripts/exact_ps_k12_partition.py:288-314` 同樣只讀入 `ref/alt`，不驗 SNV。

結論要寫成：

> 在本次 production runner 路徑中，每個 endpoint 都由 PASS、biallelic、single-base somatic VCF catalog 產生；全量 catalog audit 另確認 REF/ALT 為互異的 canonical DNA bases，因此本批 endpoints 是 sSNV。但 extractor 尚未自行拒絕非 A/C/G/T 與 REF=ALT，`build_strict_ps_hp_regions.py` standalone 使用時也未自行驗證此條件。

建議 P1：抽出共用 `validate_ssnv_catalog_row()`，在 extractor、strict builder、partitioner 三處 fail closed。

## 5. 問題 3：container 是否確實不跨 chromosome／PS／HP？

### 5.1 Python strict builder

`InterSubMod/scripts/build_strict_ps_hp_regions.py`：

- `:272-273`：每個 molecule row 必須等於 CLI 指定的 dataset/chrom。
- `:293-300`：只接受 HP family `1/2` 且 PS 非空。
- `:309-311`：container key 是 `(hp, phase_set)`；chrom 已由整個 function invocation 固定。
- `:359-369`：每個 `(hp, exact PS)` accumulator 獨立建 graph/components。

因此完整 container 是：

```text
dataset × chromosome × exact nonmissing PS string × primary HP family
```

### 5.2 Cross-language verifier

`InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/cpp/strict_endpoint_graph.hpp`：

- `:28-40`：`ContainerKey = dataset, chrom, phase_set, hp_family`。
- `:224-249`：全域 molecule ID 去重，非 primary HP/missing PS 排除，再建立 container。
- `:271-285`：只在該 container 內收集 fixed nodes 與 endpoint pairs。

HCC1395 Python/C++ frozen production parity：

```text
22/22 chromosome parity PASS
endpoint/component mismatch total = 0
```

### 5.3 兩個需要說清楚的細節

1. **HP 是 family，不是 exact raw subtype。**  
   `extract_lossless_read_linkage_collapsing.py:167-177` 把 raw `1/1-1/1-2` 合為 family 1，把 `2/2-1/2-2` 合為 family 2。論文應寫「primary HP family 1/2」，不要寫「每個 raw HP tag 完全分離」。

2. **Python/C++ 對 missing-PS sentinel 的 standalone 規則不完全相同。**  
   Python builder `:298-300` 只把空字串視為 missing；C++ `strict_endpoint_graph.hpp:141-143` 也排除 `. / NA / N/A / none / None / null`。預設 extractor `:180-186, :448-451` 會把 `.` 正規化成空字串，所以 production parity 不受影響；但外部 catalog/molecule input 仍應共用同一 normalization。

## 6. 問題 4：region 是否必須 `k>=2`？

### 6.1 Graph registry 與 tree region 是不同集合

`InterSubMod/tools/strict_endpoint_graph.py:274-295` 對所有 active nodes 做 union-find，因此無 qualifying edge 的 node 仍會形成 singleton component。

`InterSubMod/scripts/build_strict_ps_hp_regions.py`：

- `:436-467`：`k>1` → `READ_LINKED_MULTISITE_REGION`；`k=1` → `UNLINKED_SINGLETON_COMPONENT`。
- `:478-506`：只有 `k>1` 是 `PRIMARY_PS_AWARE` 並送往 partition；singleton 是 `ABSTAIN_SINGLETON_UNLINKED`。
- `:526-531`：每個 multisite component 以 retained edges 重新做一次 connectivity assertion。
- `:573-585`：receipt 另驗 singleton 永不 primary、primary 必須 multisite 且至少一條 retained edge。

`InterSubMod/scripts/run_layered_v4_strict.py:637-667` 再次驗證同一 gate。

精確用詞：

```text
component registry = k>=1（含 singleton，以守恆）
source tree-eligible region W = connected component 且 k_W>=2
```

### 6.2 HCC1395 frozen production 實數

threshold = 3、chr1-22：

| 層級 | 數量 |
|---|---:|
| candidate sSNV S | 79,687 |
| active node memberships | 62,651 |
| all graph components | 39,846 |
| singleton abstentions (`k=1`) | 28,384 |
| tree-eligible source regions W (`k>1`) | 11,462 |
| retained endpoint edges | 76,202 |
| `k_W>12` source regions | 90 |
| max `k_W` | 153 |

本次獨立重算 39,846 個 threshold-3 components：

```text
endpoint/catalog violations = 0
mixed PS/HP components = 0
source multisite disconnected components = 0
singleton/tree-role mismatches = 0
```

## 7. 問題 5：source region → block → tree input 的語意

### 7.1 Partition 保留什麼？

`InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/scripts/exact_ps_k12_partition.py`：

- `:317-454`：只接受 threshold 3、`PRIMARY_PS_AWARE`、known exact PS、`PS_HP1/2`；每個 source `component_id` 形成一個 unit。
- `:492-563`：molecule fixed calls 依 `(HP family, exact PS, site_index)` route；同 unit 至少兩個 fixed calls才成 constraint。
- `:634-750`：每個 unit 個別送入 ordered-hypergraph partition，block row 保留 `component_id/PS/HP`。
- `:753-818`：驗 site assignment、`k_B<=12`、small unit one block、cross-PS=0、cross-HP=0、constraint disposition mass conservation。
- `:821-827`：正式執行強制 threshold=3、`max_block_size=12`。

所以 `k>12` 的正確敘述是：

```text
完整 source region W 保留不變
        ↓
切成一個或多個 k_B<=12 computational blocks B
        ↓
每個 B 做 local mutation-state tree
```

block 不是新的 biological region，也不能把多個 block local trees 無條件接成唯一 global clone tree。

### 7.2 Partition 沒有保證什麼？

`_validate_outputs()` 的 checks（`:753-818`）**沒有**：

- child block induced strict-edge graph connectivity；
- child block 至少一個 retained constraint；
- child block 至少一個 MINREAD-supported multi-site pattern。

HCC1395 實測：

| 指標 | 數量 |
|---|---:|
| source W | 11,462 |
| bounded blocks | 11,712 |
| `k_B=1` child blocks | 12 |
| `k_B>=2` 但 retained weight=0 | 134 |
| `k_B>=2` induced threshold-3 graph disconnected | 4 |

4 個 disconnected blocks 全都沒有 adapter-supported pattern，依現行 gate 不會進 solver；因此目前 frozen HCC1395 tree input 沒有 disconnected block，但這是 `MINREAD` filter 的結果，不是 partition connectivity invariant。

### 7.3 Adapter 做的是「重新投影」，不是單純轉檔

`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/exact_ps_partition_to_mlhp.py`：

- `:213-296`：production mode 要求 strict receipt PASS，並用 SHA-256 證明 strict membership 正是 partition input。
- `:331-342`：block 不得缺 PS/primary HP，且 positions/k 一致。
- `:377-402`：讀 dispositions，只對 `retained` constraint 驗 block/PS/HP 與累計 retained weight。
- `:404-462`：接著重新讀取**全部 authoritative molecule rows**，按 block route；`O/D/S/L/X` 與 block 中未觀察位置都轉成 `X`。
- `:464-480`：排除 `k=1`；每個 exact vector count 必須 `>=MINREAD`，否則整個 block 不進 tree input。
- `:481-539`：無 `X` 的 vector 成 full population，有 `X` 的成 partial pattern；輸出明示 `vaf_eligible=false`。

這造成一個文件／實作落差：

- module docstring `:4-7` 寫「Only constraints classified as retained are projected」；
- 實作 `:404-462` 則從全部 molecule rows 重新投影，包含原 constraint 被標成 cut/unavoidable 的 molecule 在各 child block 的局部資訊。

後者可作為合理的 local marginalization policy，但必須明講，不能說 tree 只用了 retained constraints。

### 7.4 HCC1395 依現行 adapter 規則的獨立投影

以 `MINREAD=3`、11,712 blocks 與 3,201,802 canonical molecule rows 重算：

| funnel | 數量 |
|---|---:|
| all blocks | 11,712 |
| `k_B=1` excluded | 12 |
| `k_B>=2` | 11,700 |
| no supported exact vector, excluded | 110 |
| predicted tree-input blocks | 11,590 |
| tree-input blocks strict-induced connected | 11,590 |
| tree-input blocks strict-induced disconnected | 0 |
| retained weight=0 但仍因重新投影進 tree input | 75 |
| 只有 single-fixed-site supported patterns | 383 |
| 至少一個 multi-fixed-site supported pattern | 11,207 |

解讀：

- **正面**：現行 HCC1395 實際送 solver 的 blocks 都保留 source strict connectivity。
- **待補 gate**：75 個 zero-retained blocks 仍能由 raw projection 進 tree，證明「partition retained evidence」與「solver-visible evidence」是兩個不同 denominator。
- **待補標記**：383 個 block 雖有 strict graph 背景，但 MINREAD 後 solver 看見的每個 exact vector 最多只有一個 fixed site；應標成 `partial-only / no MINREAD-supported multi-site pattern`，不可用來宣稱組合共現。

### 7.5 Solver 最後能宣稱什麼？

`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/layered_tree_reconstruction.py`：

- `:115-138`：每個 block/HP 的 full/partial patterns 個別送入 minimal-tree enumeration。
- `:181-224`：保留所有等價最小候選樹、ambiguity/capped/verification 狀態。
- `:444-468`：analysis contract 明示 exact PS、HP、bounded block 與 no cross-PS join。
- `:458`：claim scope 是 regional mutation-state trees，不是 confirmed cell clones。

`InterSubMod/scripts/run_layered_v4_strict.py:1183, :1318-1321` 也設定同一 claim ceiling。

因此 A–B/B–C chain 最多支援：

> A/B/C 屬於同一 exact-PS/HP strict read-linkage source region，且 downstream 可在 bounded block 內枚舉與現有 full/partial mutation-state patterns 相容的最小候選樹。

不能升級成：

> A/B/C 已證實是同一 clone，或 A→B→C 是唯一腫瘤演化路徑。

## 8. 需要修正的項目

### P0：正式論文／口試用詞

1. 把「region 都是 `k>=2`」改成：
   - graph component registry 含 `k=1`；
   - tree-eligible source region 才要求 `k>1`。
2. 把「A–B、B–C 串起來表示同一 clone」改成：
   - 表示同一 evidence-connected analysis region；
   - clone identity 仍未知。
3. 把「`k>12` region 被切掉」改成：
   - source region 保留；
   - 只切 bounded computational blocks。
4. 把「tree input = strict multisite components」改成：
   - `MINREAD-supported k>=2 bounded blocks derived from strict multisite components`。

### P1：程式 fail-closed

1. strict builder 與 partitioner 都要驗 `ref/alt` 是合法 single-base SNV。
2. Python/C++ 共用 missing-PS normalization。
3. partition receipt 增加：
   - `block_induced_strict_graph_connected`；
   - `block_retained_weight_positive`；
   - `block_has_minread_multisite_pattern`。
4. adapter receipt 分開報：
   - source strict edge evidence；
   - retained partition constraint evidence；
   - raw-projected solver-visible pattern evidence。

### P2：文件與 runner receipt

1. 修正 adapter docstring「only retained constraints projected」與實作不一致。
2. `InterSubMod/scripts/run_layered_v4_strict.py:1059, :1309` 的 `tree_input` 契約要明示 bounded block、`MINREAD` 與 pattern-supported gate。
3. `build_strict_ps_hp_regions.py:570` 的 connectivity receipt flag 現為常數 `True`；雖然 `:526-531` 已 fail-closed assertion，仍建議改為由實際計數產生。
4. adapter strict-receipt production gate目前沒有獨立整合測試；應新增 tampered receipt/membership SHA、cross-PS、disconnected block、zero-pattern block fixtures。

## 9. 測試與 production 證據

### 9.1 執行命令

輸入路徑：

```text
InterSubMod/tests/test_strict_endpoint_graph.py
InterSubMod/tests/test_build_strict_ps_hp_regions.py
InterSubMod/tests/test_run_layered_v4_strict.py
InterSubMod/research/20260718_k_gt8_read_supported_segmentation/tests/test_extract_lossless_read_linkage_collapsing.py
InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_exact_ps_k12_partition.py
InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_exact_ps_partition_to_mlhp.py
InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_cpp_exact_ps_k12.py
```

命令：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m pytest -q \
  tests/test_strict_endpoint_graph.py \
  tests/test_build_strict_ps_hp_regions.py \
  tests/test_run_layered_v4_strict.py \
  research/20260718_k_gt8_read_supported_segmentation/tests/test_extract_lossless_read_linkage_collapsing.py \
  research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_exact_ps_k12_partition.py \
  research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_exact_ps_partition_to_mlhp.py \
  research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_cpp_exact_ps_k12.py
```

實際輸出：

```text
.................................................................. [100%]
66 passed in 16.45s
```

### 9.2 Frozen production receipts

Strict source regions：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
20260723_production_exact_ps_strict_read_linkage/
hcc1395_strict_regions_v2/
```

Partition：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
20260723_production_exact_ps_strict_read_linkage/
hcc1395_partition_v2/
```

實際摘要：

```text
strict chr receipts: 22/22 PASS
strict Python/C++ parity: 22/22 PASS, mismatch=0
partition receipts: 22/22 PASS
partition Python/C++ comparisons: 22/22 PASS, mismatch=0
partition constraints: 57,629 rows / molecule weight 285,596
retained/cut/unavoidable weight: 281,685 / 1,242 / 2,669
```

### 9.3 All-technical-dataset production 狀態

截至 2026-07-23 13:24（Asia/Taipei）的只讀 snapshot：

```text
HCC1395             strict 22/22 + summary PASS（既有 frozen）
HCC1395_DORADO      extraction 22/22, strict 22/22, summary PASS
COLO829             extraction 22/22, strict 22/22, summary PASS
H1437               extraction 22/22, strict 22/22, summary PASS
H2009               extraction 22/22, strict 22/22, summary PASS
HCC1937             extraction 5/22 completed receipts, strict 0/22
HCC1954             extraction 0/22, strict 0/22
```

因此本報告的 code-path verdict、HCC1395 chr1-22 parity 與五組 completed datasets 的 strict audit 是完整的；七 technical datasets 的跨樣本 production 結論仍是 `PARTIAL`，不可寫成 7/7 完成。

## 10. 建議口試的一句話

> 我們先在同一條染色體、同一個 exact phase set 與 primary haplotype family 內，要求同一條 canonical read 在兩個不同 sSNV 都有明確 REF 或 ALT call，才對該位點對提供一票；至少三條不同 reads 支援的 pair 才形成 edge。A–B 與 B–C 可以用 graph connectivity 定義同一個候選分析區域，但這只代表有證據路徑，不代表 A/B/C 已證明存在於同一細胞或同一 clone。當區域超過 12 個位點時，我們保留完整 source region，只為計算把它切成最多 12 位點的 local blocks，再用 full/partial read patterns 枚舉候選 mutation-state trees。
