<!--
建立時間: 2026-07-10 (Asia/Taipei)
目標: 稽核 layered reconstruction v2 從 LongPhase-S 產物、multi-locus read census、per-HP solver、CN/read-AF/methyl auxiliary 到 HTML workstation 的實際程式鏈、判斷、參數、I/O、例外與驗證缺口
處理範圍: 全 7 dataset / 6 biological samples 的 current-working-tree code path；未重跑全量 BAM，執行驗證限於靜態追蹤、既有 manifest lock、solver golden、合成 unit tests 與最小 probe
關聯檔案: InterSubMod/research/20260710_layered_reconstruction_v2/{00_INDEX.md,pre-decision-audit.md,implementation-notes.md,input_manifest_lock.json}; InterSubMod/docs/CURRENT_FOCUS.md; InterSubMod/docs/experiments/INDEX.md
-->

# Layered reconstruction v2：完整程式鏈與判斷稽核

> **PARTIAL EVIDENCE — Task type B（Comprehensive validation）的程式碼全鏈稽核已完成，但依任務限制未重讀全量 BAM，因此不是全資料 execution validation。**

用 SCQA + stage/decision ledger：

- **Situation**：目前已有 7-dataset manifest、per-HP-family tree solver、V1–V7、CN/read-AF auxiliary 與 workstation 程式。
- **Complication**：這些元件並未形成單一、契約一致的 production chain；legacy 與 v2 入口同時存在，部分 user-facing HTML 仍使用舊欄位、舊數字或禁止的排序邏輯。
- **Question**：是否已經能無模糊地宣稱「從輸入到 HTML 的完整流程皆合理且全面驗證」？
- **Answer**：**目前 NO-GO。核心 v2 JSON solver 可進入修正後的全量驗證，但完整 end-to-end/HTML canonical claim 被 8 個 P0 阻擋。**

## 0. 一句結論與裁決

**完成 current-working-tree 全鏈稽核 — solver golden 與 7 個合成測試通過，但 upstream `somatic_pass` 生命週期、group→family 守恆、region capped 表示、input-lock enforcement、獨立 verifier、HTML 候選全集/欄位契約與 stress 計數仍有 P0；修正前不可把現有 workstation 或「4,000/4,000」當完整驗證證據（影響：高，信心：高）。**

| 範圍 | 裁決 | 信心 | 理由 |
|---|---|---:|---|
| v2 solver 的小型手算案例 | PASS | 高 | 8/8 golden，V1–V7 皆通過。 |
| current v2 Python regression | PASS with warning | 高 | 7/7 `unittest` 通過；有未關檔 `ResourceWarning`。 |
| 7-dataset input lock | PASS as snapshot | 高 | 7 dataset / 6 biological samples、7/7 preflight pass，manifest SHA 相符。 |
| v2 JSON core full run | **未執行** | — | 本任務明確禁止全量 BAM；不能用靜態稽核代替 full execution。 |
| end-to-end provenance | **NO-GO** | 高 | upstream 同一路徑覆寫/刪除 `somatic_pass.vcf.gz`，runner 又未強制 lock。 |
| canonical HTML/workstation | **NO-GO** | 高 | 現有 builder/已生成 index 混用舊 schema、硬編 scientific claims，且 client 端用顯示前綴做 read-AF 排序/forced-edge。 |
| methyl L3 | **not evaluated in v2** | 高 | v2 driver 正確標 `not_evaluated`；現存 methyl script 仍讀舊 topology 並會排序，不能接入 v2。 |

### 證據標記

- `[F-L1]`：直接由 current-working-tree 程式碼或 immutable input 讀得的事實。
- `[O-L1]`：本次實際執行的小型命令/probe 結果。
- `[F-L2]`：repo/KB 明文契約，但尚未以本次全資料重算。
- `[I-L2]`：由兩個以上 L1 事實推得的工程/科學風險。
- `[U-L?]`：無法在不讀全量 BAM/不做外部驗證下確定。

## 1. 任務 gate、研究目標與邊界

| 問題 | 本次答案 |
|---|---|
| Task type | **B Comprehensive validation**；scope 是全 7 dataset 的完整程式鏈，不允許用 subset 冒充 full validation。 |
| 服務目標 | G3（read-level epigenetic 價值邊界）、G4（多樣本一致/可重現）、G5（外部可驗證工程證據）。 |
| Thread D | 有關：methyl 只能是 bounded auxiliary，不能 rank/confirm tree。 |
| Thread B 撤回範圍 | 非核心；舊 pooled topology/`is_somatic` 路徑只作 legacy 比對，不能回流 canonical。 |
| KDE-corrected | 本流程不直接使用 KDE；methyl GMM 也不是 KDE。 |
| VCF caller AF | 核心建樹不使用 caller AF；auxiliary 使用 family-specific read ALT fraction，明確不是 purity/CN-corrected CCF。 |
| 長計算/BAM | 有；本次依限制只做 code audit 與最小合成驗證，**未跑 full BAM**。 |

已讀 SoT：`InterSubMod/docs/CURRENT_FOCUS.md`、`InterSubMod/docs/experiments/INDEX.md`、`InterSubMod/docs/README.md`、`InterSubMod/research/20260710_layered_reconstruction_v2/{00_INDEX.md,pre-decision-audit.md,implementation-notes.md}`。外部工具契約透過 knowledge KB 查閱 `05_tools/longphase-s.md` 與 `06_workflows/phasing-workflow.md`；KB 的 HP 完整詞彙包含 `./1/2/3/4/1-1/2-1/1-2/2-2`，因此 HP4 不是可默認為 unphased 的未知字串。

## 2. 稽核 snapshot 與可重現性

### 2.1 Git/worktree

- Git HEAD：`4fb9e742482b63a660de19a1f1bd07d49d713111`。
- 本次依 parent 指示稽核 **current working-tree snapshot**，不是 HEAD-only。
- 核心檔當時為 dirty/untracked；下列 SHA 才是本報告實際讀到的版本：

| 檔案 | SHA-256 |
|---|---|
| `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_linkage_genomewide.py` | `a5b1938bd197b3c350a8341b13182e6beb52af57fca040d6b753d04e207881af` |
| `.../sm_multilocus_combinations.py` | `f28c3b4ec5211350f90fbbc91e0c30a806b6f1900cc71d17b5021cab48d319ca` |
| `.../tree_enumeration_solver.py` | `2c21ba497c2c74e7612facd988627de0d4dc35d3a27b487e87b76ded1c025d3c` |
| `.../layered_tree_reconstruction.py` | `4eaafa0ba7d4396d08615757ac0b8fb1d2ef3de923c5c727e4d31b99f70fb150` |
| `.../build_region_view.py` | `35e816e8b364c0f272925f4de0ef2e3b5f725b151b0cc71e012797db8e65cba6` |
| `.../run_layered_7samples_newbb.sh` | `fa6a050a4cd5e5ae4a2bbe1d2277626ef11f9f601164d024727226cc7111130c` |
| `.../verify_layered_v2.py` | `50f49ba5756414162668381f1920583798012b3307c34b2cea28ee6bf581479b` |

### 2.2 Input lock

`InterSubMod/research/20260710_layered_reconstruction_v2/input_manifest_lock.json`：

- schema `2.0`；7 datasets；6 biological samples；`all_pass=true`。
- manifest SHA `fd1a9aec...15dc`，與目前 `layered_v2_input_manifest.json` 相符。
- CN：HCC1395/HCC1395_DORADO=SEQC2；H1437/H2009/HCC1954=SAVANA；COLO829/HCC1937=unavailable。
- VCF autosomal counts：HCC1395 80,234；DORADO 79,120；COLO829 36,585；H1437 73,243；H2009 150,370；HCC1937 15,915；HCC1954 19,743。

這只證明 **lock 建立當時** 的檔案、header/hash/index/CN label 基本檢查通過；runner 並未消費這份 lock，見 P0-04。

## 3. 真正的 end-to-end call graph

```text
paired tumor/normal BAM + ClairS VCF + germline phased VCF + reference
  │
  ├─ scripts/pipeline/steps/01_longphase_s.sh
  │    ├─ bcftools PASS prefilter
  │    ├─ LongPhase-S somatic_haplotag (-q 20, --tagSupplementary)
  │    ├─ <sample>_tagged.bam + .bai
  │    └─ <sample>_tagged_sc.vcf → intended PASS VCF / TP / FP
  │
  ├─ validate_layered_v2_inputs.py        [目前是手動 preflight，不在 runner 內]
  │    └─ input_manifest_lock.json
  │
  └─ run_layered_7samples_newbb.sh        [canonical v2 JSON runner]
       ├─ 每 sample 五個 chromosome shard 並行
       │    └─ sm_multilocus_combinations.py
       │         ├─ import sm_linkage_genomewide.make_groups
       │         └─ import sm_linkage_genomewide.per_read_calls
       │              → mlhp_part_1..5.json
       ├─ layered_tree_reconstruction.py
       │    └─ tree_enumeration_solver.py
       │         → layered_reconstruction_<sample>.json
       ├─ build_region_view.py
       │         → layered_region_view_<sample>.json
       └─ verify_layered_v2.py
                 → verification_summary.{json,tsv}

手動/未串接 auxiliary：
  read_af_tree_ordering_multisample.py → read_af_tree_ordering.json
  methyl_auxiliary_annotation.py       → 舊 topology 輸入；目前不相容 v2
  build_layered_workstation.py         → standalone HTML；目前不相容完整候選/新欄位
  build_layered_per_sample.py          → 讀舊固定路徑並產 layered_workstation/
```

**關鍵事實 `[F-L1]`**：canonical runner 目前停在 JSON verifier；它沒有呼叫 read-AF、methyl 或 HTML builder。因此 repo 中沒有一個單一命令可重建「輸入到 canonical HTML」的全鏈。

## 4. Canonical、helper、legacy 的邊界

| 元件 | 實際角色 | 是否可進 canonical claim |
|---|---|---|
| `scripts/pipeline/steps/01_longphase_s.sh` | upstream tagged BAM/VCF 產生器 | **需要 P0-01 修正後才可** |
| `validate_layered_v2_inputs.py` | preflight/lock | 可用，但需串入 runner |
| `run_layered_7samples_newbb.sh` | v2 core orchestrator | 是；current snapshot 已有 top-level 與 child `set -euo pipefail`（line 3, 174） |
| `sm_linkage_genomewide.py::make_groups/per_read_calls` | v2 helper library | 是 |
| `sm_linkage_genomewide.py::main` | TP/FP + normal `is_somatic` 舊流程 | **否，legacy** |
| `sm_multilocus_combinations.py` | canonical Stage 1 | 是，待 P0/P1 |
| `tree_enumeration_solver.py` + `layered_tree_reconstruction.py` | canonical Stage 2 | 是，待 P1 |
| `build_region_view.py` + `verify_layered_v2.py` | canonical Stage 3/4 | 是，待 P0/P1 |
| `sm_region_integration.py` | 舊 pairwise/integration | 否；目前只可能供 HCC1395 `U3_linkage_regions_full_span` 額外計數 |
| `run_layered_6samples.sh` | 舊六樣本 runner；未設 `SM_SOMATIC_VCF`、`VERIFY_EVERY=20`、無 fail-fast | **否，必須明示 deprecated** |
| `funnel_census_{HCC1395,7samples}.py` | 2026-07-10 report-layer hotfix | 不應取代 v2 embedded funnel；7-sample CN 表仍寫僅 HCC1395 有 CN |
| `full_census_hierarchy.py` | 舊 HCC report | 否；`isolated=S-retained` 再度混合 out-of-scope/cap/read-unsupported |
| `read_af_tree_ordering_multisample.py` | v2 exploratory auxiliary | 可保留為 auxiliary，不是 validation/CCF |
| `ccf_and_cn_multisample.py` | deprecated wrapper | 否；名稱不得回流 UI/claim |
| `methyl_auxiliary_annotation.py` | 舊 topology + TP/FP methyl script | **不相容 v2** |
| `build_layered_workstation.py` / per-sample driver | transition-era HTML | **目前 NO-GO** |

## 5. Stage-by-stage 程式與方法細節

### S-1. LongPhase-S upstream contract

**入口**：`InterSubMod/scripts/pipeline/steps/01_longphase_s.sh`。

**輸入**：tumor BAM、normal BAM、caller somatic VCF、germline phased VCF、reference、truth VCF/BED。truth 只用來切 TP/FP；v2 操作宇宙是 PASS VCF，不是 truth label 子集。

**實際判斷/設定**：

- line 124–150：先把 caller VCF 過濾至 PASS。
- line 155–168：LongPhase-S `somatic_haplotag`，`-t <threads>`、`--tagSupplementary`、`-q 20`、`--output-somatic-vcf`。
- line 192–209：tagged BAM 建 index，另跑 haplotag QC。
- line 220–246：再由 `<sample>_tagged_sc.vcf` 取 PASS。
- line 248–280：以 truth 做 `bcftools isec` 產 TP/FP，**這不是 v2 backbone**。

**HP 輸入契約 `[F-L2]`**：KB 列出 `./1/2/3/4/1-1/2-1/1-2/2-2`；HP3 是 unresolved somatic-integrated，不是第三條已確認 germline lineage。repo invocation 另明確要求 tag supplementary。

**重大例外**：同一 `somatic_pass.vcf.gz` 同時被用作 caller PASS 暫存（line 126）與 LongPhase `_sc.vcf` PASS 結果（line 235），最後又於 line 291 刪除。詳 P0-01。

### S0. Manifest/preflight

**入口**：`validate_layered_v2_inputs.py`。

**已檢查**：file/index 存在、VCF 僅 PASS、ClairS source header、VCF/VCF-index hash、BAM header/reference-dictionary/index hash、CN labels 是否在 gain/loss/loh/neutral。

**未檢查**：

- sample name 唯一/非空、恰 7 dataset/6 biological 的 policy；
- tagged BAM 實際 HP vocabulary、HP4/extended tag 分布、MM/ML tag；
- BAM content hash（只有 metadata/header/index；可接受但要明示）；
- VCF record-key 與 raw caller PASS / LongPhase `_sc` 的轉換一致性；
- CN interval coordinate convention、排序/重疊/完整性與 categorical/integer BED 彼此對齊；
- manifest `cn_source` 與實際 BED provenance 的外部來源證據。
- upstream canonical bundle 的 `manifest/run_context.json`、`manifest/upstream_dependency.tsv` 與 `metrics/completeness.tsv`；依 `InterSubMod/docs/standards/20260314_big7_canonical輸出與延續研究規範_01.md`，`completeness.tsv` 才是 canonical 完整性入口，單看 BAM/VCF 存在不足以取代它。

### S1. PASS sSNV → multi-locus read census

**入口**：`sm_multilocus_combinations.py main(chroms,out_path)`；runner 對 chr1–22 切 5 shards。

#### 5.1.1 位置宇宙與 grouping

1. `load_somatic_universe` 透過 Tabix 讀指定 chromosome，只保留 `len(REF)==len(ALT)==1`；不再以 normal `is_somatic` 二次篩選（lines 28–43）。
2. 任何 `ValueError/OSError` 直接回空 list，錯誤原因不進 JSON（lines 34–42）。
3. `make_groups` 以相鄰位置差 `<=TIER_R` 做 connected component；預設 `TIER_R=50,000`。**總 span 不封頂**（`sm_linkage_genomewide.py:145–156`）。
4. group `<2` sSNV → positional singleton。
5. group `>MAX_SNV=8` → 取 span 最小的連續 8 個；tie 由 Python `min` 取最左 window（`sm_multilocus_combinations.py:108–114`）。被移除的 sSNV 計 cap-excluded。

#### 5.1.2 BAM read 判斷

`per_read_calls` 對整個 group span 做一次 pileup：

- 預設 `MAPQ_MIN=20`、`BASEQ_MIN=0`（`sm_linkage_genomewide.py:32–35`）。
- `stepper="samtools"`，pysam 0.23.3 本機 doc 顯示 default flag filter = unmapped/secondary/QC-fail/duplicate；**supplementary 未在 default filter 中**，orphans 預設忽略。
- deletion/refskip/query-position missing 跳過；base=REF → `R`，base=ALT → `A`，其他/未覆蓋 → `X`。
- dict key 是 `query_name`；同 QNAME 多 alignment/同位置時後寫覆蓋先寫，沒有 conflict counter。
- HP/PS 也是同 QNAME 後寫覆蓋；PS 在 canonical Stage 1 之後未使用。

#### 5.1.3 Family、support 與保留

- HP `1*→family1`、`2*→family2`、恰為 `3→family3`，其他/缺值→`none`（lines 69–83）。因此 `1-2/2-2` 能歸 1/2，但 HP4 被歸 `none`。
- 每 QNAME 只要覆蓋至少一位點便進 `reads_by_hp`。
- full vector：沒有 X；`combo[_fam]` 計數。
- `subread[_fam]`：**目前所有 vector 都計入，包括 full vector**（lines 129–145）。
- `MINREAD=3`：full population 與 subread genotype 各自保留 count≥3。
- group retention：`pooled n_full_cov_reads>=3` **或** pooled 任一 subread genotype≥3（lines 148–154）。
- family unit 只由 per-family genotype count≥3 產生（lines 155–161）。pooled 通過不保證任一 family 通過。
- categorical CN 在 `positions[len(positions)//2]`（偶數 k 是右側中央 sSNV）取值；CN BED 存在但沒有 segment hit→neutral，BED 不存在→unavailable。

#### 5.1.4 Stage 1 JSON schema（ad hoc）

```text
schema_version: "2.0"
params: MINREAD/MAX_SNV/TIER_R/MAPQ_MIN/BASEQ_MIN/backbone/CN source
input_funnel: scope_input/singleton/pre-cap/cap-excluded/read-unsupported/retained + conservation
n_groups_analyzed
aggregate_all / aggregate_clean_loh_neutral / aggregate_strict_neutral
groups[]:
  chrom/start/end/span/n_sSNV/pre_cap_n_sSNV/positions
  n_full_cov_reads/populations/subread_groups/col_coverage
  populations_by_hp/subread_groups_by_hp/col_coverage_by_hp/reads_by_hp
  cn/cn_data_available/dominant_hp_count/same_hp
```

沒有 JSON Schema；runner 只用 `jq` 檢查 `schema_version` 與一個 funnel boolean。

### S2. Per-family minimal-tree enumeration

#### 5.2.1 數學模型

對每個 region×family：

- full genotype `g∈{R,A}^k` → observed node `altset(g)`；所有 full node 必須在樹中。
- partial `p∈{R,A,X}^k` → subcube `C(p)`；最終 node set 至少與每個 subcube 相交。
- root=`∅`；每個非根 node 必須有一個只少 1 個 ALT bit 的 predecessor（unit-flip closure）。
- 目標：最小化 `|N \ ({root} ∪ full)|`，也就是 hidden/partial-supported node 數。
- 對所有最小 feasible node sets，枚舉每個 node 的所有合法 parent choice；candidate tree 數為各 node predecessor 數乘積、再跨 node set 加總。
- read count 在 `MINREAD` threshold 後只決定 presence；solver 不用 count 作 likelihood/weight。

來源：`tree_enumeration_solver.py:98–205`。

#### 5.2.2 Search/cap/fallback

- `extra_cap=4` hidden nodes；每一 hidden level 的組合數 budget=`150,000`；任一 level 超出即 capped（lines 141–161）。
- capped 時 `_greedy_closure` 產一個覆蓋用 fallback node set；它不是完整 candidate set。
- solver 函式預設 `tree_cap=32`；v2 driver 以 `ANALYSIS_TREE_CAP=0` 覆寫成全候選，另 `DISPLAY_TREE_CAP=32` 只保留顯示前綴。
- 無 RNG；但 `_parent_choice_trees` 對 set/frozenset 未顯式排序，輸出順序不是語言層保證。

#### 5.2.3 Class 與 V1–V7

| 分類 | 精確條件 |
|---|---|
| `determined` | non-capped、存在 recurrence-free minimal solution、exact `n_trees=1` |
| `ambiguous` | non-capped、存在 recurrence-free minimal solution、`n_trees>1` |
| `recurrence_required` | 所有 minimal solution 都需同一 mutation bit 在不同分支重複取得 |
| `capped` | search budget/hidden cap 內未完成；不得宣稱唯一/完整 |

V1 unit-flip arborescence；V2 full/partial coverage；V3 full nodes consistency；V4 獨立檢查 `e_min-1` 不可行；V5 獨立重算 **candidate count**；V6 recurrence flag/class consistency；V7 determined 不過度宣稱。V5 不比較 candidate edge-set identity/digest。

#### 5.2.4 Driver role/denominator

`layered_tree_reconstruction.py:115–210`：

- `mutation_bearing`：任一 MINREAD-supported full/partial genotype 含 A。
- primary lineage：mutation-bearing HP1/HP2。
- reference-only control：family 1/2/3 且不 mutation-bearing；排除 primary denominator。
- HP3 mutation-bearing：`unresolved_H3_auxiliary`，排除 primary denominator。
- `none`：unphased auxiliary。
- driver 把 `subread_groups_by_hp` 所有 key 傳作 partial，因此 full vector 被重複當 singleton subcube；不改可行解，但 `n_partial`/UI 語義錯。
- V1–V7 預設 `VERIFY_EVERY=1` 全 non-capped unit；current runner line 123 明示 full V1–V7。

### S3. L2 CN、region view 與 funnel

#### 5.3.1 CN post-tree 判斷

只對 mutation-bearing、non-capped、`recurrence_required` unit 執行：

| 條件 | verdict |
|---|---|
| CN missing/unknown/unavailable | unavailable |
| integer gain CN≥3 | artifact_drop |
| integer loss CN=0 | artifact_drop |
| categorical gain、無 integer | gain_confounded |
| neutral | candidate_keep，記 cn_int=2 |
| loss | loss_confounded |
| LOH 且 max family read-AF>0.7 | `likely_artifact(高VAF)`（仍是 LOH_unresolved） |
| LOH 且 ≤0.7 | `likely_recurrence(低VAF)`（仍是 LOH_unresolved） |

來源：`layered_tree_reconstruction.py:43–62`。0.7 沒有 prereg calibration/evidence；read-AF 也不是 caller AF/CCF。

CN lookup point 是 `(first_position+last_position)//2`（line 119），不同於 Stage 1 的右中央 sSNV。既有 HCC1395 7,928-group snapshot 的只讀 probe 發現 12/7,928（0.151%）categorical state 不一致；這是舊 snapshot 的風險量化，不是 v2 full-run 結果。

#### 5.3.2 Region determinacy

`build_region_view.py:35–45` 只看 primary lineages：

1. 無 primary → `no_primary_lineage`；
2. 全 determined → `all_determined`；
3. 任一 recurrence → `has_recurrence`；
4. 否則任一 capped → `has_capped`；
5. 其餘 → `has_ambiguous`。

因此同時含 recurrence 與 capped family 的 region 只標 recurrence，會隱藏 incomplete enumeration（P0-03）。

#### 5.3.3 Funnel

autosomal sSNV 應滿足：

```text
L1 all PASS universe
 = L2 non-autosomal out-of-scope
 + L3 autosomal positional singleton
 + L5 cap-excluded sSNV
 + L5 read-unsupported sSNV
 + L6 retained sSNV
```

Stage 1 shard 內先檢 `scope_input = singleton + cap_excluded + read_unsupported + retained`；Stage 2 合併；region view 再從 VCF 重數全體/autosomal/non-autosomal。這是 sSNV-level 守恆，**沒有檢查 retained group 必然映射到 ≥1 family unit/region view**。

#### 5.3.4 Region-view join

- raw mlhp 以 `region=chrom:start-end` join；glob 未排序、duplicate key 後寫覆蓋。
- malformed shard JSON broad-except 後靜默跳過。
- region 只由 layered `detail` 建；Stage 1 retained 但沒有 family unit 的 group 不會出現在 `regions[]`。
- `is_primary_lineage` 缺值時 fallback 舊 `is_lineage`，即使上游是舊 schema，也會輸出 schema 2.0。
- U5 tree numerator只算 non-capped primary；U6 hidden 卻含 capped fallback。
- JSON 非 atomic write，部分寫入在 crash 時可能留下可見檔。

### S4. Runner 與 final verifier

#### 5.4.1 Runner

`run_layered_7samples_newbb.sh` current snapshot：

- top-level line 3 與 xargs child line 174 都有 `set -euo pipefail`；本報告**不沿用舊版「child 未 fail-fast」結論**。
- default `PARALLEL_SAMPLES=2`；每 sample 同時啟 5 BAM shard，峰值最多約 10 個 Stage 1 reader。
- 5 shard 分割互斥覆蓋 chr1–22。
- input 只檢 BAM/BAI/VCF/CN 檔存在；沒有要求 VCF index，也未呼叫 `validate_layered_v2_inputs.py`。
- 建 run-root copy 的 `input_manifest.json`，但 worker `manifest_value` 與 final verifier 仍讀原始 `$MANIFEST`。
- input hash含 VCF、BAI、VCF index（若存在）、CN；不 hash 大 BAM content。
- `code.sha256` 不含 runner 自己、input validator、read-AF、methyl、HTML builders。
- `PROFILE` 宣告但未使用；status 多 process append 無 file lock；無 trap 寫 final failed summary。
- core 成功後只跑 final JSON verifier，不產 HTML/auxiliary。

#### 5.4.2 Verifier

`verify_layered_v2.py` 檢：schema字串、role counts、zero skipped/fail、candidate generated count、exact-shape 非空、funnel booleans、missing-CN unavailable、region multiplicity。

它**沒有**：

- 由 VCF/mlhp 獨立重算 funnel；
- 重新執行 solver 或比對所有 candidate edge digest；
- 驗證 manifest/lock/input hashes；
- 檢 retained group vs region view；
- 檢 CN available sample 的 categorical/integer一致；
- 檢 duplicate regions/shards；
- 驗證 HTML/auxiliary。

所以名稱雖為 final verifier，實際是輸出內部 consistency checker，不能單獨支持「獨立全面驗證」。

### S5. read-AF ordering auxiliary

`read_af_tree_ordering_multisample.py`：

- 僅 non-capped、primary、ambiguous unit。
- 從 mlhp 重建 `full/partial`，`tree_cap=0` 重枚舉全部 candidate；只比 `n_trees`，不比 tree digest。
- family per-locus read-AF=`ALT/(REF+ALT)`，含 partial reads。
- 每條 parent→child 邊：對 child 新取得 mutation 與 parent 已有 ancestor mutations累加 `AF_ancestor-AF_child`。
- softmax temperature grid 0.025/0.05/0.1；threshold 0.6/0.75/0.9；violation margin 0.03/0.05/0.1；default 0.05/0.6/0.05。
- 無 RNG、無 purity/CN correction、無 read-depth uncertainty；輸出已正確標 exploratory/not CCF/not independent validation。
- `strict_neutral_reach_fraction` 分母是**所有 prepared units**，不是 neutral prepared units；名稱不足以讓讀者立即分辨兩種 reach。

### S6. methyl auxiliary

current `methyl_auxiliary_annotation.py` **不是 v2 模組**：

- 讀 `topology_per_region.json` 與舊 TP/FP VCF，不讀 layered v2 JSON/manifest。
- 會對 ambiguous region 產生 `ordering`，直接違反 v2 `L3.prohibited_uses=[tree_ranking,...]`。
- 參數：MAPQ20、每 read≥3 CpG、group≥12 reads、GMM peak separation≥0.2、小 component≥4、READ_CAP=500、HP dominance 0.7。
- GMM `random_state=0,n_init=1`；可重現 seed，但對第一 500 alignment 的 BAM order 敏感。
- `tb.fetch(chrom,s,e+1)` 把 1-based `s` 直接當 0-based start，左界 off-by-one。
- 只接受覆蓋所有 sSNV 的 read；QNAME 後寫覆蓋，但 `n_read` 仍逐 alignment 加，可能未達 500 unique read 就停止。
- CN 只把 `gain` 當 confound；missing/LOH/loss 不足以被排除。
- broad exception 靜默；chunk outputs 沒有 merge/completeness gate；JSON 無 schema/provenance/sample/input hash。

v2 driver 把 L3 固定為 `not_evaluated(bounded auxiliary;禁止 rank/confirm)` 是目前唯一合理 canonical 狀態。

### S7. HTML/workstation

#### 5.7.1 Per-sample driver

`build_layered_per_sample.py`：

- 固定讀舊 HCC pilot與 `/big7_disk/.../multisample_subclone`，不讀 immutable v2 run-root。
- builder fail 時把 sample 標 pending 後繼續，整體仍可 exit 0。
- 用 `python3` 而不是當前 interpreter。
- index ribbon 硬編 `主分母=region(7100)`；欄位/說明仍寫 HCC1395 23,810 等舊值。
- 已生成 `InterSubMod/docs/methodology/_assets/layered_workstation/index.html` 仍含 pre-P0 counts（如 HCC1395 multi-HP 4,443、COLO829 77.6%、HCC1954 90.8%）與互相矛盾的 sSNV 說明，不可視為 v2 canonical output。
- README 仍指示 `run_layered_6samples.sh` 重生，該 runner 是舊 backbone。

#### 5.7.2 Single-sample builder

`build_layered_workstation.py` 的主要問題：

- line 117 表頭硬編 `HCC1395 計數`；line 118 舊敘述；line 129 宣稱「過半區」；line 164 硬編「94% CN-altered / 98% VAF / 77–86% methyl」等非 JSON 注入的 scientific claims，與檔頭「所有數字從 JSON、無 hardcode」矛盾。
- `regTag` 仍認 `no_germline_lineage`，新 view 是 `no_primary_lineage`；dashboard 讀 `rd.no_germline_lineage`，會把真實 no-primary 顯為 0。
- fail 欄讀 `L1.n_verify_fail`，新欄為 `n_verification_fail`。
- missing CN 的 `unavailable` 被放進 `other`，但表只顯 `unknown`，所以 COLO829/HCC1937 CN table 不守恆。
- U5 numerator=non-capped tree sum，denominator卻用所有 primary lineages；capped sample平均值偏低。
- U6 把 capped greedy fallback 的 hidden 數納入精確總數。
- client 重新以 raw ALT 判 mutation-bearing/root-only/c，定義與 core MINREAD-supported per-family role 不同；`c` 又跨 family pooled，重新引入 allelic/clonal 混淆。
- `obs_subreads` 含 full vector，但 UI 一律標 `partial(X=未覆蓋)`。
- `ccfBlock` 對 `trees[]`（最多 display 32）做 softmax排序；不是全候選，且稱 CCF。
- stable/forced edge 也只由 stored/display prefix 交集計算；`n_trees>display` 時可高估 forced edge。
- capped/partial verification 一律顯 `V✗`，而不是 capped=N/A、partial=未完整。

結論：HTML 層不是單純樣式問題，而是會改變科學判斷；必須在 core JSON 契約修好後重建。

## 6. 單位、分母與公式 registry

| ID | 單位/分母 | 精確定義 | 不可混用 |
|---|---|---|---|
| U1-all | VCF sSNV | final operational PASS VCF 的 1bp REF/ALT records，含非 autosome | 不是 truth TP；不是生物真值 |
| U1-auto | autosomal sSNV | U1-all 中 chr1–22 | runner/solver scope |
| G-pre | positional group | 相鄰 gap≤50kb connected component，cap 前 | total span 可>50kb |
| G-ret | retained group | cap 後且 pooled full≥3 或 pooled subread genotype≥3 | 目前不保證有 family unit |
| R-call | read/QNAME | pileup 中至少一位點有 R/A/OTHER call 的 QNAME | 不是 alignment count；supplementary policy未鎖 |
| Pop-full | genotype population | 同 family、完整 R/A vector count≥MINREAD | support magnitude過 threshold 後不進 solver |
| Subcube | partial observation | 應為含 X 的 vector count≥MINREAD | current JSON 錯含 full vector |
| F-unit | family unit | region×family，有 full 或 subread group | 包含 controls/H3/none |
| P-unit | primary lineage | mutation-bearing HP1/HP2 family unit | primary determinacy denominator |
| Region | region | 至少一個 layered detail unit 的 `chrom:start-end` | current view可能漏 retained zero-family group |
| Tree | candidate tree | non-capped unit 的所有 minimal node-set/parent-choice arborescences | display prefix不是全集 |
| Hidden | inferred node | `N \ ({root}∪full)` | capped fallback不得和 exact hidden混算 |
| Read-AF reach | prepared ambiguous primary unit | default softmax top≥0.6 | strict-neutral目前分母仍是 all prepared |

必守公式：

1. `U1_auto = singleton + cap_excluded + read_unsupported + retained_sSNV`。
2. `retained_groups = groups_with_≥1_family_unit + explicit_zero_family_branch`；目前缺此式。
3. `primary units = mutation-bearing HP1 + mutation-bearing HP2`；reference/H3/none 不進分母。
4. `region determinacy` 必須能同時保留 completeness risk（capped）與 biology/model flag（recurrence），不能單一互斥 label 隱藏其中一項。
5. `avg trees = Σ n_trees(non-capped primary) / n_noncapped_primary`，不能除以全 primary。
6. `exact hidden sum` 只對 analysis-complete non-capped units；capped fallback需另欄。

## 7. 最基本 decision registry

| ID | 判斷 | Default/threshold | 來源 |
|---|---|---:|---|
| D01 | 分析 chromosome | chr1–22 | `run_layered_7samples_newbb.sh:19–20` |
| D02 | group split | adjacent gap>50,000 才切 | `sm_linkage_genomewide.py:145–156` |
| D03 | singleton | group sSNV<2 | `sm_multilocus_combinations.py:101–104` |
| D04 | densest cap | 最小 span 的連續 8 位點，tie 最左 | `sm_multilocus_combinations.py:108–114` |
| D05 | mapping quality | MAPQ≥20 | `sm_linkage_genomewide.py:34,104–106` |
| D06 | base quality | BASEQ≥0 | `sm_linkage_genomewide.py:35,96` |
| D07 | pileup flags | pysam samtools default；supplementary 未明排 | `sm_linkage_genomewide.py:96` + local pysam doc |
| D08 | allele | exact REF→R、exact ALT→A、其餘→X | `sm_linkage_genomewide.py:108–109`; `sm_multilocus...:129–145` |
| D09 | family | `1*→1,2*→2,3→3,else→none` | `sm_multilocus...:69–83` |
| D10 | support | genotype count≥3 | `sm_multilocus...:19,149,155–160` |
| D11 | group retention | pooled full≥3 OR pooled subread-group存在 | `sm_multilocus...:148–154` |
| D12 | CN Stage1 anchor | 右中央 sSNV | `sm_multilocus...:165` |
| D13 | primary lineage | mutation-bearing HP1/HP2 | `layered_tree_reconstruction.py:137–151` |
| D14 | solver objective | minimal hidden nodes | `tree_enumeration_solver.py:141–156` |
| D15 | solver search cap | hidden≤4、每 level≤150,000 combinations | `tree_enumeration_solver.py:98,141–164` |
| D16 | analytical candidate cap | 0=all | `layered_tree_reconstruction.py:78,132` |
| D17 | display candidate cap | 32 | `layered_tree_reconstruction.py:79,176` |
| D18 | verification cadence | 1=每 non-capped unit full V1–V7 | `layered_tree_reconstruction.py:80,134–173` |
| D19 | recurrence | mutation bit在≥2 edge獨立取得 | `tree_enumeration_solver.py:88–94` |
| D20 | CN Stage2 anchor | first/last coordinate arithmetic midpoint | `layered_tree_reconstruction.py:119` |
| D21 | LOH heuristic | max family read-AF>0.7 | `layered_tree_reconstruction.py:58–61` |
| D22 | L3 v2 | not evaluated；禁止 rank/confirm | `layered_tree_reconstruction.py:343–346` |
| D23 | sample concurrency | 2 samples × 5 shards | `run_layered_7samples_newbb.sh:10,98–110,173–174` |
| D24 | read-AF default | temp .05 / top .60 / margin .05 | `read_af_tree_ordering_multisample.py:146,166–180` |
| D25 | methyl legacy | MAPQ20/MINGRP12/SEP.2/MINSUB4/cap500 | `methyl_auxiliary_annotation.py:27,85–100,116` |

## 8. P0 — 阻擋 comprehensive/canonical claim

### P0-01. Upstream `somatic_pass.vcf.gz` 同名暫存、覆寫、刪除

- **證據 `[F-L1]`**：`01_longphase_s.sh:124–150` 先寫/index `somatic_pass.vcf.gz`；`:235–246` 再把 LongPhase `_sc.vcf` PASS 寫回同一路徑；`:289–292` 最後刪除它與 index。
- **影響**：同名檔同時代表 caller input PASS 與 LongPhase output PASS；index 可能 stale/overwrite衝突；current script rerun後又刪除 layered manifest 所需的 canonical backbone。現存 manifest 檔如何保留下來無法由 current script單獨重建。
- **可驗收修正**：
  1. 暫存改為 `caller_pass_input.vcf.gz`；final 改為唯一 `longphase_somatic_pass.vcf.gz`（或明確 canonical 名）。
  2. final VCF/index不得刪；TP/FP cleanup只刪 isec temp。
  3. lock記 raw caller、prefilter、`_sc.vcf`、final PASS 四個 SHA 與 record-key diff。
  4. 小型 fixture 重跑兩次，兩次都 exit 0、index對應 final hash、manifest target存在。

### P0-02. Stage 1 retained group 可完全沒有 per-family unit，且 verifier 不會發現

- **證據 `[F-L1]`**：retention 用 pooled support（`sm_multilocus...:148–154`），family output另各自 threshold（`:155–161`）；layered detail只遍歷有 per-family key 的 family（`layered_tree_reconstruction.py:124–130`）；region view只從 detail建 region（`build_region_view.py:47–52`）。
- **最小反例 `[I-L2]`**：同一 pooled genotype 3 reads分散 HP1=2、HP2=1；pooled通過、每 family均<3，group計 retained sSNV但沒有任何 unit/region。
- **影響**：sSNV funnel仍可全部 true，tree/region denominator卻漏 group；目前 final verifier不比 retained group count。
- **可驗收修正**：retention 改成「至少一個 per-family supported genotype」或新增 `zero_family_support` 明確分支；新增 invariant `retained_group_count = region/detail group count + explicit_zero_family_count`；合成 2+1 split test必不得靜默消失。

### P0-03. Region 單一 label 會用 recurrence 隱藏 capped family

- **證據 `[F-L1]`**：`build_region_view.py:39–45` priority 是 R 在 C 前。
- **影響**：含一個 recurrence family及一個 capped family的 region顯示 `has_recurrence`，`has_capped` census/UI低估，使用者看不到候選集未完整。
- **可驗收修正**：輸出獨立 booleans `any_capped/any_recurrence/any_ambiguous/all_determined`；若仍需 main label，completeness風險 capped優先；mixed R+C fixture必同時保留兩 flag。

### P0-04. Input lock 未被 runner 強制消費，worker/final verifier仍讀 mutable manifest

- **證據 `[F-L1]`**：runner line 30 copy manifest，但 `manifest_value` line 41與 final verifier line 177仍用 `$MANIFEST`；沒有呼叫 validator/lock。
- **影響**：preflight證據與真正執行輸入可分離；run進行中 manifest改動會讓不同 sample/verify讀不同版本。
- **可驗收修正**：runner啟動先跑 validator或指定 lock；核對 source manifest SHA；此後只讀 `$RUN_ROOT/input_manifest.json`；每 sample執行前比對 lock中的 VCF/index/BAI/CN hash；verification summary嵌 manifest/lock SHA。

### P0-05. Final verifier 多數相信 upstream self-report，不能稱 independent comprehensive verification

- **證據 `[F-L1]`**：`verify_layered_v2.py:63–86` 讀 `verification_status` 與 funnel boolean；`:70–79` 只比 generated count/exact shape欄；沒有重跑 solver或重算 funnel。
- **影響**：同一 implementation若共同算錯，producer與verifier仍可一致 PASS；P0-02就是實例。
- **可驗收修正**：至少獨立重算 VCF→funnel、mlhp→role/region count、all-candidate canonical edge digest；抽樣/全量 solver oracle policy明確；任何缺欄不得因 `.get()` falsey巧合通過；輸出 verifier code hash。

### P0-06. HTML 用 display 前綴做 tree ordering 與 forced-edge，直接違反全集/禁止排序契約

- **證據 `[F-L1]`**：core只存 `DISPLAY_TREE_CAP=32`（layered driver lines 175–201）；workstation `ccfBlock` 對 `L.trees` 排序（builder lines 312–332）；stable edge也對 stored trees交集（lines 380–395）。
- **影響**：n_trees>32 時 winner posterior與 forced edge都可能因漏候選而改變；UI又把 read-AF稱 CCF。即使 core analysis完整，HTML會產生另一個不完整科學結果。
- **可驗收修正**：HTML不得自行排序；只顯 core/auxiliary已驗證結果。若要 stable edge，core需基於全集輸出 `all_candidate_stable_edges` + digest；read-AF只能消費 `tree_cap=0` auxiliary且 candidate digest相同；不滿足即顯 N/A。

### P0-07. 現有 workstation/schema/數字與 v2 current contract 不一致

- **證據 `[F-L1]`**：per-sample builder fixed legacy paths與 hardcoded 7100；single builder讀舊欄位、漏 unavailable CN、硬編 scientific percentages；已生成 index保留 pre-P0 multi-HP counts與舊 somatic解釋。
- **影響**：即使 JSON修正，使用者打開現有 HTML仍會看到錯誤 denominator/claim；這正是本任務要求排除的「模糊狀況」。
- **可驗收修正**：builder只接受 immutable v2 run-root與 schema contract；所有 scientific number/claim由 versioned JSON欄位注入；missing CN列需守恆；7 sample fixture做字段/總和 snapshot test；舊 HTML加顯著 `DEPRECATED/PRE-P0` ribbon或移至 archive（依 repo archive規範，不直接刪）。

### P0-08. 「seeded stress 4,000 cases → 0 mismatch」把 capped 未驗案例算成 PASS

- **證據 `[F-L1]`**：`full_v4v5_verification.py:91–98` 若 capped，直接 `stress_pass += 1`；summary line 103–107不另報 capped skipped。
- **影響**：4,000是 generated cases，不是 4,000個經V4/V5檢查的 cases；`CURRENT_FOCUS` 的 4,000/0 mismatch claim證據口徑不成立。
- **可驗收修正**：輸出 `generated/checked/capped_skipped/pass/fail`；`all_checked_pass`只對checked，另設最小coverage gate；重新執行並把 CURRENT_FOCUS/報告改成精確口徑。

## 9. P1 — 下一次 full run 前應修

| ID | 問題與證據 | 影響 | 驗收條件 |
|---|---|---|---|
| P1-01 | `[F-L1]` HP4 被 `germ_family(4)→none`；KB明列 HP4 | 把「somatic on both germline families」與真正 unphased混在一起 | 定義 H4 role；先做7 sample tag census；HP4 fixture與denominator policy |
| P1-02 | `[F-L1]` `subread_groups` 收所有 vector，包括 full | `n_partial`、UI、trace錯；full constraint重複 | partial欄只容許含X；另保留 all-read vector欄；schema/test拒絕無X partial |
| P1-03 | `[O-L1]` capped fallback若 full含 root，top-level `n_hidden=len(N)-len(full)-1` 少1；probe回 top=0/tree=1 | capped U6/trace不一致 | 統一 `len(N-({ROOT}|full))`；root+capped fixture；verifier檢 top/tree一致 |
| P1-04 | `[O-L1]` Stage1右中央sSNV vs Stage2座標中點；舊HCC snapshot 12/7928 CN不同 | 同一region L1與L2可能用不同CN | manifest宣告唯一anchor；兩層共用；boundary fixture |
| P1-05 | `[F-L1]` CN `.bed` 實際像1-based inclusive，但validator未宣告/驗證 | off-by-one與跨來源不可審核 | manifest加 coordinate_system；validator檢sorted/nonoverlap；轉檔 provenance |
| P1-06 | `[F-L1]` supplementary未明排、QNAME conflict後寫覆蓋；LongPhase又 `--tagSupplementary` | chimeric/supplementary可能橋接或覆蓋call | 明定include/exclude；輸出flag/QNAME conflict census；合成primary+supp fixture |
| P1-07 | `[F-L1]` BASEQ=0、support=3、MAX_SNV=8、50kb與LOH .7缺 prereg sensitivity | 閾值可能驅動結論 | 在params+spec逐項寫來源；至少MINREAD/MAX_SNV/TIER_R/MAPQ/BASEQ敏感度 |
| P1-08 | `[F-L1]` `_parent_choice_trees` 未排序，digest又不sort tree list；V5只比count | candidate digest跨環境不穩、相同數量不同集合可漏 | canonical sort node/edge/tree；V5比edge-set digest；重跑hash一致 |
| P1-09 | `[F-L1]` positive ANALYSIS_TREE_CAP仍可能 `verification_status=full_pass`，只有final H6另擋 | 中間JSON `all_V1V7_pass`會過度表述 | full_pass需同时 analysis_complete；或拆 verification/exhaustiveness兩軸 |
| P1-10 | `[F-L1]` region view fallback legacy `is_lineage`，raw JSON broad-except/duplicate overwrite | 舊輸入可偽裝schema2；read evidence靜默缺 | 嚴格schema與required fields；duplicate/malformed即fail |
| P1-11 | `[F-L1]` SAVANA只有categorical CN；SEQC2另有integer CN | L2 artifact resolution深度跨樣本不等 | 每verdict帶 evidence_level/source；跨樣本表分層，不混 raw count |
| P1-12 | `[F-L1]` read-AF只比n_trees不比digest；softmax稱posterior但未校準；strict-neutral分母模糊 | winner/reach易被誤讀為機率/CCF | 比candidate digest；改稱relative weight或校準；同時報all/neutral兩分母 |
| P1-13 | `[F-L1]` methyl讀舊topology、off-by-one、first500 alignment、chunk無merge gate | 若誤接v2會產錯誤ranking/residual | 明示DEPRECATED；另寫v2-only bounded negative-screen module與schema/merge gate |
| P1-14 | `[O-L1]` unit tests有 open file `ResourceWarning`；多處JSON非atomic | FD累積/partial output | `with open`、temp+fsync+rename；warning-as-error regression |
| P1-15 | `[F-L1]` runner status append無lock、無trap、PROFILE未用；VCF index不要求 | 失敗run審計不完整 | flock/per-sample status後merge；trap final FAIL；index required；profile落地或移除 |
| P1-16 | `[F-L1]` HTML avg trees numerator non-capped、denominator all primary；U6含capped fallback | 基本分母錯 | 加 exact denominator欄；HTML只讀已算好 ratio，fixture守恆 |
| P1-17 | `[F-L2]` v2 preflight只驗artifact，未讀canonical `run_context/upstream_dependency/completeness.tsv` | 上游run可能是partial或provenance未閉合，卻仍進v2 | manifest納入canonical bundle root；validator要求/鎖定三份canonical provenance與`completeness_state` policy |

## 10. Error/fallback registry

| 位置 | 現在的行為 | 是否足夠 |
|---|---|---|
| Stage1 VCF fetch error | 回空 list | **否**；可能把整chrom當0且只靠後面funnel碰巧攔截 |
| BAM pileup/開檔 error | 未捕捉，process fail | 是；runner可攔，但需status detail |
| CN BED missing | unavailable | 是；這一點 current v2 正確 |
| CN segment no hit | neutral | 需 source-specific contract；不可默認所有 BED皆如此 |
| Solver budget/hidden cap | greedy fallback + capped | 合理作展示/診斷；不得進exact統計/forced edge |
| Region raw mlhp malformed | 靜默 continue | **否**；會丟read evidence |
| Optional integration missing | 靜默不加欄 | 可接受，但 schema需標 `not_available` |
| Runner shard fail | current child/top都 `set -euo`; xargs nonzero終止 | 基本合理；缺 final FAIL summary/trap |
| HTML sample build fail | 標 pending、繼續、可能exit 0 | **否**，Task B需整體fail |
| Methyl exception | 多處 broad-except→no data | **否**，缺error census |

## 11. Seeds、截斷與平行性

| 模組 | RNG/順序 | 截斷/上限 | 平行 |
|---|---|---|---|
| Stage1 | 無 RNG；pileup/BAM順序與QNAME後寫可影響 conflict | MAX_SNV=8；MINREAD=3 | 5 chromosome shards/sample |
| Core solver | 無 RNG；set iteration順序未canonical | hidden≤4；level budget150k；analysis cap0；display32 | sample內serial units |
| Runner | 無 RNG | 無 timeout | 2 samples × 5 shard readers |
| read-AF | 無 RNG | 只 ambiguous/noncapped；全candidate re-enumeration | samples serial |
| methyl legacy | GMM seed0,n_init1；BAM順序 | first500 alignments/region | 可strided chunk但無merge gate |
| stress | seeds 20260710/1/42/777/12345，各800 | capped被錯算pass | serial |

## 12. 本次實際驗證命令與輸出片段

### V-A. Input lock

**輸入**：`InterSubMod/research/20260710_layered_reconstruction_v2/input_manifest_lock.json`。

**命令**：

```bash
jq -r '[.schema_version,.dataset_count,.biological_sample_count,.all_pass,.manifest_sha256] | @tsv' \
  research/20260710_layered_reconstruction_v2/input_manifest_lock.json
```

**實際輸出**：

```text
2.0  7  6  true  fd1a9aec8514e602e7ae1407b6e388735488a2679e80a0ac17b41527a59415dc
```

### V-B. Runner syntax/current fail-fast snapshot

**輸入**：`InterSubMod/docs/.../scripts/run_layered_7samples_newbb.sh`。

```bash
bash -n docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/run_layered_7samples_newbb.sh
```

**實際輸出**：exit code 0；current line 3與174皆有 `set -euo pipefail`。

### V-C. Solver golden

```bash
PYTHONDONTWRITEBYTECODE=1 /bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
  docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py
```

**實際輸出摘要**：8個手算案例全部 PASS；每案 V1–V7 PASS；exit 0。

### V-D. v2 regression tests

```bash
PYTHONDONTWRITEBYTECODE=1 /bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
  docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/test_layered_reconstruction_v2.py
```

**實際輸出摘要**：`Ran 7 tests ... OK`；另有 `layered_tree_reconstruction.py:238/348 ResourceWarning: unclosed file`。

環境未安裝 pytest；直接執行 `unittest` 檔完成同一測試集。

### V-E. HP family probe

```text
{1:'1',2:'2',3:'3',4:'none','1-1':'1','2-1':'2','1-2':'1','2-2':'2',None:'none'}
```

這直接確認 HP4 conflation。

### V-F. pysam pileup flag contract

本機 pysam=`0.23.3`；doc 實際輸出 default `flag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP`，supplementary 不在其中。

### V-G. capped `n_hidden` probe

合成 full=`RR,AA`，強制 `per_level_budget=0`：

```text
top-level n_hidden=0
stored tree n_hidden=1
capped=true
```

確認 P1-03 是實際程式錯誤，不只是理論。

### V-H. CN anchor只讀風險量化

**輸入**：既有 `InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/mlhp_part_*.json` + `/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed`。

**輸出**：7,928 groups 中 12（0.151%）在「右中央 sSNV」與「first/last midpoint」得到不同 CN；例 `chr2:10178748-10256938 loss→loh`。這只是既有 snapshot 的風險量化，不能當 v2 full dataset結果。

## 13. 尚未驗證的邊界

1. `[U-L?]` 7 tagged BAM 中 HP4 的實際 read/region 數；現有 haplotag QC只聚合 `hp_other`，例如 H1437有 1,638,979 `hp_other`，無法拆 HP3/4/extended。
2. `[U-L?]` P0-02 zero-family retained group在 7 dataset current v2實際發生幾次；full run前要先用 dry census/probe算。
3. `[U-L?]` raw ClairS PASS、LongPhase `_sc.vcf`與 manifest `somatic_pass.vcf.gz` 的 record-key差異；header同為ClairS不等於內容identity。
4. `[U-L?]` 7 dataset full runtime/RSS/I/O contention；本次沒有跑10-reader峰值。
5. `[U-L?]` SAVANA與SEQC2 CN coordinate/transformation的外部原始規格；檔案形態看似1-based inclusive，但repo manifest未聲明。
6. `[U-L?]` read-AF heuristic的 calibration、depth uncertainty與truth discrimination；self-consistency不是獨立驗證。
7. `[U-L?]` methyl v2的任何 biological claim；current canonical狀態就是 not evaluated。
8. `[U-L?]` canonical HTML render/runtime/accessibility；因資料契約本身先NO-GO，本次沒有把舊HTML視為可驗證終點。

## 14. 建議修正順序（Step → Verify）

1. 修 P0-01 upstream artifact命名/保存 → 驗證：小型重跑兩次，final PASS VCF/index/hash皆存在且一致。
2. 修 P0-02/P0-03 group/region守恆與多軸flags → 驗證：pooled 2+1 family、mixed recurrence+capped fixtures都精確命中。
3. 把 validator/lock串入runner並只用run-root manifest → 驗證：執行中篡改source manifest不影響run，hash mismatch會在preflight fail。
4. 強化 independent verifier與canonical candidate digest → 驗證：故意改一棵edge、一個funnel count、一個region，三者皆fail。
5. 修 HP4/subread/capped hidden/CN anchor/explicit pileup policy → 驗證：新增對應unit tests，warning-as-error 7/7 pass。
6. 重跑7 dataset core → 驗證：7/7 lock/hash、funnel、group→unit、candidate completeness、V1–V7、no silent branch全部pass。
7. 重跑 stress並分開 checked/capped → 驗證：報告精確列 denominator，不再把skip算pass。
8. 重寫/重接 auxiliary與HTML → 驗證：HTML只讀immutable v2 JSON；missing CN守恆；候選不完整時不rank/不標forced；跨7 sample field snapshot+browser smoke pass。

## 15. 最終可對外方法敘述（修正後應採的口徑）

本流程重建的是 **bulk ONT 中、局部 autosomal region、依 germline HP1/HP2 分開的 mutation-state trees**。操作宇宙是 version-locked 的 LongPhase-S/ClairS PASS sSNV；read在每個局部region被編碼為R/A/X向量，只有達MINREAD的family-specific evidence進入solver。solver枚舉所有最小hidden-node、unit-flip相容樹；candidate多於一棵即保留不可辨識性，capped則明示未完整。CN只在樹完成後註記recurrence confound；family read-AF只作exploratory relative ordering，不是CCF或獨立驗證；methyl在沒有正交per-unit驗證前只可作bounded negative screen/residual flag，不能rank或confirm tree。region tree不等於完整cell clone tree，也不等於single-cell clone confirmation。

在 P0/P1 修正與7-dataset full execution evidence完成前，對外只能說「**current-working-tree code audit完成；core solver small-case驗證通過；full end-to-end validation pending**」。
