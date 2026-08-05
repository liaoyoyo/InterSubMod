<!--
建立時間: 2026-07-23 13:15 +08:00
最後復核: 2026-07-23 13:24 +08:00
目標: 驗證 strict read-linkage 是否真的以同一 canonical molecule 的兩個不同 sSNV fixed calls 建邊、是否限制於同一 exact phase block，並核對最新 production 數據與 downstream claim boundary
處理範圍: current workspace 程式、HCC1395 strict/k12 artifacts、截至快照完成的 5/7 technical datasets
關聯檔案: InterSubMod/tools/strict_endpoint_graph.py；InterSubMod/scripts/build_strict_ps_hp_regions.py；InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/scripts/exact_ps_k12_partition.py；InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/exact_ps_partition_to_mlhp.py
-->

# 兩個 sSNV 與 exact phase block：程式及最新數據驗證

> **PARTIAL／Needs revision：strict region 建構層為 PASS；但尚不能宣稱 end-to-end tree input 全部保有「同一 read 至少兩個 sSNV」證據。**  
> 影響：高；region-level 信心：高；end-to-end 信心：中。  
> 任務類型：B Comprehensive validation；服務 G1／G4／G5。

## 1. 先講結論

### 1.1 可以確認的部分

1. **本次五組 production 的候選位點確實是 sSNV。** Extractor 程式強制 `PASS`、exactly one ALT 與 REF/ALT 長度皆為 1；本輪再獨立確認 428,759 個 catalog loci 的 REF/ALT 全為互異的 `A/C/G/T`，indel、MNV、symbolic allele、`N` 與 REF=ALT 皆為 0。要注意 producer 本身尚未 fail closed 驗 `A/C/G/T` 與 REF≠ALT。
2. **一條 direct edge 確實需要同一 canonical molecule 在兩個不同 sSNV 都有 fixed `R/A` call。** `X/L/O/D/S` 不會建 edge；同一 molecule ID 不可重複計票。
3. **production primary edge 不是一條 read 即成立，而是至少 3 個 distinct canonical molecules。**
4. **每條 edge 與每個 source region 都限制在同一 chromosome、同一 exact nonmissing `PS`、同一 primary `HP1/HP2` container。** 完成五組共 1,133,073 條 primary edges、72,742 個 tree-eligible source regions，獨立重算跨 PS、跨 HP、不連通與 hash/mass 違規皆為 0。
5. **source region 的 `k>=2`。** `k=1` singleton 被保留作 funnel conservation，但不進 tree inference。
6. **k≤12 partition 的 structural constraint 也要求同一 molecule 在同一 source unit 至少兩個 fixed `R/A` calls。**

### 1.2 現在不能直接說的部分

1. **source region 是 connected component，不代表一條 read 穿過 region 的全部 sSNV。** A–B 與 B–C 可由不同 molecules 支援而形成同一 source region；這只證明圖上有 transitive read linkage，不證明 A–C direct co-observation，也不證明同一 clone。
2. **`PS` 一致是 operational phase-block boundary，不等於已獨立驗證 phase accuracy。** Exact `PS` 可防止已知 phase blocks 被混合；switch error、錯誤 haplotag 或 HP orientation 仍需另外驗證。
3. **current adapter 可能把只有一個 fixed 位點的 supported pattern 納入 tree input。** 因此不能把 source-region PASS 改寫成「所有 solver-visible patterns 都跨至少兩個 sSNV」。
4. **current strict partition → adapter 的 production binding 尚未 PASS。** HCC1395 partition root 缺 `run_receipt.json`；partition receipt 的 identity 欄位使用 `bytes`，adapter 卻讀 `size_bytes`。
5. **目前只完成 5/7 technical datasets。** HCC1937、HCC1954 尚未形成同層級完整 receipt chain。
6. **這批 strict 程式仍是 untracked current-workspace files。** 結論綁定本文件列出的 SHA-256，不是已提交的 release commit。

## 2. 建議採用的正式定義

> 在同一 `chromosome × exact nonmissing PS × primary HP family` 容器內，若至少 3 個 distinct canonical molecules 在兩個不同、VCF-PASS、biallelic sSNV endpoints 上皆具有 `MAPQ≥20` 且 `base quality≥20` 的 fixed `REF/ALT` call，則建立一條 direct endpoint edge。由達門檻 edges 形成的 `k≥2` connected component，定義為 read-linked source region `W`；`k=1` 為 singleton abstention，不進 tree inference。

這個定義比「一條 read 穿過兩個位點就算一個區域」更精確，因為它同時說清楚：

- 判定單位是 **distinct canonical molecule**；
- edge 有 **兩個不同 sSNV endpoints**；
- primary threshold 是 **3**；
- phase boundary 是 **exact `PS × HP`**；
- `W` 是 **connected component**，不是「一條 read 跨完整區域」；
- R/R、R/A、A/R、A/A 都是 callable state evidence；不是只有含 ALT 才能建 edge。

## 3. 程式路徑逐層驗證

| 層級 | 實際規則 | 程式證據 | 判定 |
|---|---|---|---|
| Candidate sSNV | 單一 REF、單一 ALT、各長度 1；tree VCF 非 PASS 即 fail；本輪資料另驗證 `A/C/G/T` 且 REF≠ALT | `extract_lossless_read_linkage_collapsing.py:244-262`；全量 catalog audit | Current data PASS；producer 防呆待補 |
| Read call | MAPQ 20；base quality 20；REF→R、ALT→A，其餘→O/L/D/S/X | 同檔 `call_alignment()` 與 extraction runner | PASS |
| Canonical identity | primary alignment；exact sidecar join；duplicate identity 只可在 observation 完全相同時 collapse；molecule ID 重複 fail | 同檔 extraction loop；`strict_endpoint_graph.py:162-171` | PASS |
| Two-site edge | 只保留 fixed R/A；`combinations(molecule.calls, 2)` 對兩個不同 site 建一票 | `build_strict_ps_hp_regions.py:301-311`；`strict_endpoint_graph.py:162-171` | PASS |
| Exact phase container | key 為 `(hp, phase_set)`，只接受 HP1/HP2 與 nonempty PS | `build_strict_ps_hp_regions.py:293-311,359-394` | PASS；schema robustness 待補 |
| Primary threshold | `support_total>=3`；state mass = RR+RA+AR+AA | `build_strict_ps_hp_regions.py:370-394,537-547` | PASS |
| Region | retained edges 的 connected component；只有 `k>1` tree eligible | `build_strict_ps_hp_regions.py:398-531` | PASS |
| k≤12 constraint | 同 exact unit；一個 molecule 在該 unit routed fixed calls `<2` 不形成 structural constraint | `exact_ps_k12_partition.py:502-563` | PASS |
| Adapter projection | 從 authoritative raw molecules 對 block 重投影；只要 vector 不全 X 且相同 vector count≥MINREAD 即 supported | `exact_ps_partition_to_mlhp.py:431-482` | **需補 two-fixed-site gate** |
| Formal binding | strict receipt、partition identity、root run receipt | `exact_ps_partition_to_mlhp.py:210-296,564-568` | **BLOCKED** |

補充：production runner 的 upstream extractor 有 `PASS + biallelic + single-base` gate，但尚未明確拒絕非 `A/C/G/T` 或 REF=ALT；standalone `build_strict_ps_hp_regions.py:217-244` 與 `exact_ps_k12_partition.py:288-314` 更只讀 `ref/alt` 欄，沒有自行驗 allele 長度、字元與 REF≠ALT。Current production data 經全量 catalog audit 為正確，但三個 reusable components 都尚未完全 fail closed。

### 3.1 為什麼這確實是「兩個 sSNV」

`MoleculeCalls` 只接受：

- unique、strictly increasing 的 `site_index`；
- call code 只能是 `R` 或 `A`。

`EndpointPairAccumulator.add()` 對同一 molecule 的 calls 執行兩兩組合；因此每一票必然有兩個不同 endpoints。`X/L/O/D/S` 在進 accumulator 前已被移除，無法被當成「穿過」證據。

### 3.2 為什麼這確實在同一 phase block

Builder 先以 `(hp_family, phase_set)` 分 container，再各自建立 accumulator；不同 `PS` 或不同 `HP` 從資料結構上就不會進入同一 pair counter。LongPhase 官方文件亦明確說明：haplotag 將 read 指派為 HP1/HP2，並把該 read 的 haplotype block 存在 `PS` tag。因此以 exact `PS × HP` 作為 read-level phase boundary 是合理的最小防混合規則。

但「same PS」只驗證 identifier 一致，不驗證 phase call 本身正確。正式論文應寫「同一 exact PS-tagged phase block 內」，不要寫成「已證明位於真實且無 switch error 的同一 haplotype block」。

### 3.3 A–B 與 B–C 串接的正確解釋

如果：

- 3 條 molecules 直接支援 A–B；
- 另外 3 條 molecules 直接支援 B–C；
- 全部具有同一 exact `PS × HP`；

則 A、B、C 可成為同一 connected component。這在 graph region 定義上合理，但只代表：

> A、B、C 位於同一 read-linked analysis domain。

它不代表：

- 有 read 直接同時觀察 A 與 C；
- A、B、C 一定在同一細胞；
- A、B、C 一定屬於同一 biological clone；
- tree 的 parent-child 關係已唯一確定。

## 4. 最新 production 數據

完整資料快照以 2026-07-23 13:18:04（Asia/Taipei）已完成的五組資料為準；13:24:18 重新執行全量獨立稽核，數值不變且 `violations_total=0`：

| technical dataset | candidate S | primary edges | all components | singleton abstain | tree-eligible W | k>12 W | max k |
|---|---:|---:|---:|---:|---:|---:|---:|
| HCC1395 | 79,687 | 76,202 | 39,846 | 28,384 | 11,462 | 90 | 153 |
| COLO829 | 37,788 | 41,736 | 42,189 | 28,256 | 13,933 | 30 | 40 |
| H1437 | 77,080 | 140,566 | 31,588 | 15,262 | 16,326 | 904 | 138 |
| HCC1395_DORADO | 79,739 | 19,453 | 43,310 | 36,482 | 6,828 | 34 | 106 |
| H2009 | 154,465 | 855,116 | 49,392 | 25,199 | 24,193 | 4,206 | 567 |
| **合計** | **428,759** | **1,133,073** | **206,325** | **133,583** | **72,742** | **5,264** | — |

全量重算結果：

- 110/110 chromosome strict receipts：`all_pass=true`，所有 checks 為 true。
- 1,133,073/1,133,073 primary edges：兩端不同、catalog match、`support_total≥3`、same exact PS/HP。
- 72,742/72,742 tree-eligible `W`：`k≥2`，induced primary graph connected。
- 133,583/133,583 singletons：`k=1` 且不進 tree。
- 880 個 receipt input/output identities、110 個 summary→child identities：size 與 SHA-256 全吻合。
- 0 cross-PS、0 cross-HP、0 disconnected、0 state-mass、0 hash violations。

### 4.1 `phase_set` 與 sSNV catalog 實際資料品質

- 18,141,508 raw molecule rows 中，4,389,652 筆 `PS` 為空字串，已被排除。
- 非空的 `. / NA / N/A / None / null / UNKNOWN / MISSING` 等 sentinel：0。
- 非空但非正整數 `PS`：0。
- strict membership 2,169,816 rows 與 strict edges 1,564,789 rows：空／sentinel／非正整數 `PS` 均為 0。
- 428,759 candidate loci：REF/ALT 長度不為 1、非 A/C/G/T、REF=ALT 均為 0。

所以 current completed data 不受 `if not phase_set` 只擋空字串的缺口影響；但程式仍應把「positive integer PS」寫成 schema gate，而不是依賴目前輸入剛好乾淨。

## 5. 實際發現的不合理或未完成狀況

### P0：adapter 會讓 two-site structural evidence 退化

Partitioner 的 constraint 很嚴格：同一 molecule 在同一 unit 少於兩個 fixed R/A 時不形成 constraint。可是 adapter 不只讀 retained constraints，而是重新讀所有 authoritative molecule rows，將每個 block 的 raw projection 聚合；只要 pattern 有一個 R/A 且 count≥3 就能成為 supported partial pattern。

因此 `AXX...`、`RXX...` 可能進入 solver-visible input。單點 read 可以作 coverage／likelihood evidence，但不應單獨建立 topology constraint。

HCC1395 依 current adapter eligibility contract、`MINREAD=3` 的完整重算：

| adapter／pattern gate | blocks |
|---|---:|
| bounded blocks | 11,712 |
| `k=1`，排除 | 12 |
| `k≥2` | 11,700 |
| 無任何 exact vector 達 MINREAD，排除 | 110 |
| **would-be tree-input groups** | **11,590** |
| 至少一個 supported pattern fixed-observe ≥2 sites | 11,207 |
| **只有 single-fixed-site supported patterns** | **383（3.30%）** |
| supported-pattern graph 覆蓋且連通整個 block | 9,717（83.84%） |
| **supported-pattern graph 不連通** | **1,873（16.16%）** |
| pattern-connected 且 ALT-supported loci≥2 | **5,009（43.22%）** |

這些是對現行規則的 read-only 重現；因 formal adapter 尚未通過 root receipt gate，11,590 是 **would-be inputs**，不是已完成 solver run。

Strict edge 與 solver-visible pattern pair 也不是同一分母：

| HCC1395 strict edge disposition | edges | 占 76,202 |
|---|---:|---:|
| partition 後跨 block | 35,008 | 45.94% |
| 同 block但整個 block pattern unsupported | 3,170 | 4.16% |
| 同 tree block但 exact-vector 分桶後 pair 不再達 MINREAD | 7,680 | 10.08% |
| **仍由至少一個 supported vector 同時 fixed-observe兩端** | **30,344** | **39.82%** |

這不是 pair edge 計算錯誤，而是 threshold 的 aggregation grain 從「endpoint pair」換成「完整 R/A/X vector」。若要做 block-wide topology，建議明確分流：

- `fixed_count>=2`：可作 structural pattern；
- `fixed_count=1`：只作 likelihood／coverage，不可讓 block 取得 tree eligibility；
- primary topology 至少要求 supported-pattern graph 連通全部 block sites；
- 若要稱 multi-mutation topology，再要求 ALT-supported loci≥2；
- 其他 block 分成 single-mutation control、reference-only control 或 abstain。

### P0：strict partition → adapter binding 尚未可執行

實際命令：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/exact_ps_partition_to_mlhp.py \
  --partition-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_partition_v2 \
  --output /tmp/hcc1395_strict_adapter_validation.json \
  --sample HCC1395 --chroms chr1 --min-read 3 \
  --require-strict-endpoint-receipt
```

實際結果：exit 1；找不到：

```text
.../hcc1395_partition_v2/run_receipt.json
```

即使補上 root receipt，partition receipt 的 membership identity 使用 `bytes`，adapter 比較的是 `size_bytes`；需要統一 schema 並加 integration regression test。

### P1：RR-only edge 的語意要說清楚

Current edge support 計算 `RR+RA+AR+AA`。因此一條 edge 可以只由 RR observations 建立。這對「兩個 sSNV endpoints 在同一 molecule 上都可判定」是合理的；但 RR-only 不是 mutation co-occurrence evidence。

HCC1395 的 primary edges 中，18,116/76,202（23.77%）為 RR-only。若論文要說「突變共同出現」，需改用含 ALT 的 state support 或另報 ALT-bearing sensitivity；若只說「read-linked mutation-state region」，保留 RR 是合理的。

### P1：程式版本尚未 release 化

本輪主要 current-workspace SHA-256：

| 檔案 | SHA-256 |
|---|---|
| `tools/strict_endpoint_graph.py` | `df3d6f37615b0bfa5d382a082152cdb16ce206eaa075cd3c0dd418b5c320e37e` |
| `scripts/build_strict_ps_hp_regions.py` | `7260a7631b30cbb5e4878159583b8a5b27a153de07e8d001303417dd2f29aedd` |
| `scripts/run_layered_v4_strict.py` | `2c2a113be341a7e7d4a57226f9390ccd1e2f46c47b20657f54d6aa4f3d9edcdd` |
| `exact_ps_partition_to_mlhp.py` | `7044f42af8599476545a8a519363c00e95785a8c00b13070038253cc25c5e06a` |
| lossless extractor | `2ca7ccb67c89e816fae9284f4e2ba21b186378105086c6b0128ed5445a133490` |

前四份與其主要 tests 在 `git status` 都是 `??`。在 commit、frozen manifest 與 full receipt binding 前，不宜稱為正式 release。

### P1：runner 與 adapter 的文字契約需同步實作

- Adapter docstring 寫「only retained constraints are projected」，實作卻重投影全部 raw molecule rows。
- Runner receipt 把 tree input 簡寫成 multisite strict component；實際上是「由 strict source component 衍生、經 k≤12 partition、且至少一個 exact R/A/X vector 達 MINREAD 的 bounded block」。
- `determined` 只代表在目前 mutation universe 與 solver model 下候選唯一，不代表有 read 同時跨兩個突變位點。Synthetic `AX, k=2` 也可得到 determined 單一候選。

## 6. 測試與交叉實作

### 6.1 Python regression tests

輸入：

- `tests/test_strict_endpoint_graph.py`
- `tests/test_build_strict_ps_hp_regions.py`
- `research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_exact_ps_partition_to_mlhp.py`
- `research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_verify_big7_input_binding.py`

命令：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m pytest -q \
  tests/test_strict_endpoint_graph.py \
  tests/test_build_strict_ps_hp_regions.py \
  research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_exact_ps_partition_to_mlhp.py \
  research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_verify_big7_input_binding.py
```

實際輸出：

```text
27 passed in 1.38s
```

測試已包含 A–B/B–C transitive component、X 不建 edge、duplicate molecule fail closed、>50 kb 不作硬切與 cross-PS rejection。

較廣的 extractor／runner／partition／C++ suite 另執行：

```text
66 passed in 16.45s
```

涵蓋 strict core、builder、runner、lossless extractor、k≤12 partition、adapter 與 C++ partition tests。現有 suite 尚未把 only-single-fixed pattern、pattern graph disconnected 與 strict-edge→solver-visible-pair conservation 設為 fail gate。

### 6.2 Python／C++ independent parity

HCC1395 chr1–22 的獨立 C++ strict graph 重建為 all-PASS：

- observed edge rows：117,760；
- components：39,846；
- edge mismatch：0；
- component／role mismatch：0。

這比只重讀 Python receipt 更強，因為第二個實作重新建 edge 與 component，再逐列比對。

## 7. 建議修正順序（Step → Verify）

1. **統一 partition identity schema：`size_bytes`。**  
   → 驗證：strict-bound adapter 不再因欄位名稱失敗，並對錯誤 size/SHA 各有 fail-closed test。
2. **建立完整 `run_receipt.json` 並綁定 22 chromosome partition／comparison receipts。**  
   → 驗證：HCC1395 `--require-strict-endpoint-receipt` exit 0；output receipt `all_pass=true`。
3. **adapter 加 structural／likelihood evidence 分流。**  
   → 驗證：`AXX×3` 只進 likelihood；`AAX×3` 可建立 pattern edge。
4. **每個 primary topology block 加 supported-pattern connectivity gate。**  
   → 驗證：pattern edges 可連通全部 block sites；否則標 `ABSTAIN_PATTERN_DISCONNECTED`。若稱 multi-mutation tree，再要求 ALT-supported loci≥2。
5. **把 strict pair constraints 與 exact-vector patterns 分開保存。**  
   → 驗證：receipt 可守恆報告 strict edges 的 internal／cross-block／solver-visible disposition，不再把「molecule 還在」誤當「pair relation 還在」。
6. **extractor／strict builder／partitioner 共用並強制驗證 sSNV catalog。**  
   → 驗證：indel、MNV、symbolic allele、N、REF=ALT fixtures 全部 fail closed。
7. **`phase_set` 改為 positive-integer schema validation。**  
   → 驗證：空、`.`、`NA`、`None`、0、負數與非數字都 fail/exclude；現有 completed data 數字不變。
8. **凍結程式 SHA、commit、重跑 integration tests。**  
   → 驗證：receipt 內 producer hashes 等於 frozen commit files。
9. **等 7/7 production 完成後重跑全量 independent audit。**  
   → 驗證：154/154 chromosome receipts、7 genome summaries、全 invariant 0 violations。

## 8. 論文與口試可用說法

### 可直接使用

> 本研究先將 reads 依 chromosome、exact phase-set identifier 與 primary haplotype family 分開。在同一容器內，只有當至少三個不同 canonical molecules 都能在兩個不同 sSNV 位點作出高品質 REF/ALT 判定時，才建立位點間的 direct edge；再以這些 edges 的 connected components 定義 read-linked source regions。這個 region 表示局部 mutation-state evidence 可透過 reads 串接，不代表所有位點由同一條 read 覆蓋，也不直接等同於 biological clone。

### 修正後才能使用

> 每一個送入 tree solver 的 block，都至少具有同一 molecule 跨兩個 sSNV 的內部結構證據。

Current adapter 尚未滿足這句話，完成 P0 gate 後才可使用。

更嚴格的 primary claim 應寫成：

> 每個納入 multi-mutation topology 分析的 block，其達 MINREAD patterns 所形成的 fixed-site graph可連通全部 block sites，且至少兩個 loci 有 ALT support；不足者分流為 control 或 abstain。

## 9. 參考與可追溯來源

- LongPhase 官方文件：haplotag 將 reads 指派至 HP1/HP2，並把 read 的 haplotype block 存在 `PS` tag：<https://github.com/twolinin/longphase>
- HTS specifications canonical repository：<https://github.com/samtools/hts-specs>
- 內部 KB：`/big8_disk/Knowledge/06_workflows/phasing-workflow.md`
- 最新獨立資料稽核：`InterSubMod/research/20260723_strict_two_ssnv_phaseblock_validation/agents/latest_data_audit.md`
- 程式路徑稽核：`InterSubMod/research/20260723_strict_two_ssnv_phaseblock_validation/agents/code_path_audit.md`
- Tree-input 語意稽核：`InterSubMod/research/20260723_strict_two_ssnv_phaseblock_validation/agents/tree_input_semantics_audit.md`

## Claim ceiling

本輪可證明：完成五組資料的 strict source-region engineering contract 正確，且每條 primary edge 都是 same exact PS/HP container 內至少三個 canonical molecules 的兩個 fixed sSNV endpoint co-calls。  
本輪不可證明：七組全量完成、每個 final solver block 都保有兩位點 structural evidence、唯一 tree topology、真實 clone/subclone 數、跨 HP cell-level pairing，或 phase truth 無 switch error。
