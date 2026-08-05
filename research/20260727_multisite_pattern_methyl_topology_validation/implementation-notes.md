<!--
建立時間: 2026-07-27 09:37
更新時間: 2026-07-27（Task B 全範圍驗證、read 聚集圖與交付完成）
目標: 記錄多 sSNV pattern × methyl association spec 實作中的決策、偏離與折衷
處理範圍: extraction / statistics / topology overlay / HTML
spec_id: 20260727_multisite_pattern_methyl_topology_validation
status: validated
advisory: on
關聯檔案:
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/00_INDEX.md
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/pre-decision-audit.md
-->

# Implementation Notes：多 sSNV pattern × 甲基關聯

## 🔵 設計決定

### [2026-07-27 09:37] Primary grain固定為 exact PS × exact raw HP

- **Status**: Accepted
- **Decision**: We will keep raw HP tags (`1`, `1-1`, `1-2`, `2`, `2-1`,
  `2-2`, `3`, `4`, `.`) separate in the primary analysis.
- **Consequence**: 現有 topology的 HP family只提供定位；同一 family row可掛多個
  raw-HP evidence bundles，但不可合併計數。
- **Revisit if**: 有獨立 phaser證明特定 raw tag可交換。
- **Evidence tier**: L1

<!-- BEGIN USER-SPECIFIED -->
**Decision**: 甲基只作 association overlay，永不重算 topology方向、edge或selected tree。
**DO NOT change**: 使用者已於本 spec 前確認。
<!-- END USER-SPECIFIED -->

### [2026-07-27 09:37] Partial R/A/X不作硬填補

- **Status**: Accepted
- **Decision**: We will retain X-containing signatures as subcube/terminal
  interval evidence and never impute a complete topology state.
- **Revisit if**: 有獨立 base-quality-aware posterior且事先註冊。
- **Evidence tier**: L2

### [2026-07-27 09:37] Methyl artifact由 qname SHA重新綁定

- **Status**: Accepted
- **Decision**: We will treat sparse `hp_raw`/R-A calls as label authority and
  methylation.csv as methyl authority; legacy/canonical reads.tsv is QA only.
- **Consequence**: 既有 BERNOULLI matrix不得跨 tile拼貼；從 long-form CpG union
  後重算 distance。
- **Revisit if**: 新 canonical artifact同時帶 exact raw HP且 receipt綁定。
- **Evidence tier**: L1

### [2026-07-27 10:31] Primary artifact改用 July-universe all-sSNV v2

- **Status**: Accepted；取代下方「Canonical優先、attested fallback」作為
  primary 路徑，但保留後者為 reference sensitivity。
- **Decision**: We will use the receipt-validated 20260715 all-sSNV v2 artifacts
  aligned to the same 7-dataset, chr1–22, LongPhase-S recalibrated PASS
  biallelic-sSNV universe as the July topology.
- **Authority root**:
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/intersubmod_all_ssnv_v2_verification_fix/`
- **Coverage**: 469,849/469,849 marker loci contain reads.tsv,
  methylation.csv and BERNOULLI matrix according to upstream receipts；本輪仍對
  frozen candidate marker逐檔計算 SHA-256。
- **Consequence**: artifact內 `hp=0` 是刻意設計，不可當標籤；所有 HP/PS/pattern
  仍以 strict sparse qname digest重新綁定。
- **Evidence tier**: L1

### [2026-07-27 10:58] Multiple-testing family與 permutation解析度先驗分層

- **Status**: Accepted before full run
- **Decision**: confirmatory family固定為 pair full4與 k≥3 formal units；其餘
  k=2 any-two-state完整掃描，但放在 secondary family，不升格 robust claim。
- **Reason**: 999 permutations的最小可報 p為0.001；若把數萬 secondary pairs
  與數百 confirmatory units混為一個 BY family，離散解析度會在數學上禁止任何
  confirmatory判定。兩階段 adaptive refinement只追加 pre-gated confirmatory
  floor cases，並保留 requested/realized次數。
- **Consequence**: comprehensive denominator仍涵蓋所有 formal n5 units；
  reader-facing robust association只來自 confirmatory family。
- **Evidence tier**: L2

### [2026-07-27 10:49] Adaptive permutation上限凍結

- **Status**: Accepted before full-run results
- **Decision**: confirmatory先跑999次；僅 `p=1/1000` floor且預先通過
  effect、PERMDISP、geometry、n8、equal-N與CpG-rarefaction gates的單元，
  以同一 deterministic seed提高到49,999次。secondary不精煉。
- **Interpretation**: 49,999次仍無法解析BY門檻者標
  `resolution-limited`，不把 permutation floor外推成更小p值。
- **Evidence tier**: L1

### [2026-07-27 10:55] BERNOULLI artifact parity門檻凍結

- **Status**: Accepted before formal catalog results
- **Decision**: 每個 formal marker以 SHA-256 deterministic抽至多16 reads，
  使用所有 CpGs重算 BERNOULLI；invalid mask零差異，finite cell最大絕對誤差
  `≤5×10⁻⁴`。
- **Reason**: frozen matrix輸出6位小數、methylation.csv輸出4位小數；此門檻
  容納序列化誤差但仍足以攔截公式或索引漂移。
- **Smoke evidence**: HCC1954 chr22的5/5 formal markers、600個pair cells
  PASS，invalid mismatch=0，最大誤差 `4.19×10⁻⁵`。
- **Evidence tier**: L1 engineering

### [2026-07-27 11:10] Restricted structure-test API完成

- **Status**: Implemented and verified
- **Files**: `include/core/Stats.hpp`、`include/core/StructureTest.hpp`、
  `src/core/StructureTest.cpp`、`tests/test_structure_bootstrap.cpp`。
- **Result**: 新增 within-strata permutation contract、PERMANOVA R²、
  requested/realized permutation audit、restricted PERMDISP centroid重算與
  NOT_EVALUABLE fail-closed；legacy unrestricted overload保留。
- **Verification**: targeted 18/18；完整 C++ `run_tests` 258/258 PASS。
- **Evidence tier**: L1

### [2026-07-27 11:10] HCC1954 chr22端到端 smoke完成

- **Status**: PASS（DEMO/engineering smoke，非 validation evidence）
- **Input**: 15,671 sparse rows、HCC1954 chr22 topology與同-universe
  20260715 all-sSNV artifacts。
- **Output**:
  `/tmp/intersubmod_pattern_methyl_hcc1954_chr22_v2_20260727/`。
- **Observed**: 3,610 read projections、73 exact raw-HP strata、2 formal n5
  units、5 unique markers全數 catalog PASS；1 unit可評估、1 unit因
  post-exchangeability support不足而 NOT_EVALUABLE。
- **Why it matters**: 驗證 census→formal selection→artifact hash→multi-tile
  CpG union→restricted statistics的 schema可連通；不取代全7×22正式長跑。
- **Evidence tier**: L1 engineering / L5 scientific

### [2026-07-27 12:26] 獨立紅隊後改為完整 fail-closed data binding

- **Status**: Implemented before formal methyl run
- **Decision**: analyzer在消費時逐 unit重算 candidate
  identity、complete/partial counts與 state counters，並重驗 reads、
  methylation及 topology JSONL的 size/SHA-256。
- **Reason**: 只驗 manifest或四欄 join仍可能讓兩個各自 hash-valid、但來自不同
  census run的檔案被混用。
- **Verification**: drift/mismatched-identity對抗測試均 fail closed。
- **Evidence tier**: L1 engineering

### [2026-07-27 12:26] Common basis改為保守 monotone fixed point

- **Status**: Implemented before formal methyl run
- **Decision**: methylation前凍結所有 `n≥5`正式 complete states；
  exchangeability、read coverage與 state-CpG coverage只允許刪 row/CpG並
  共同收斂。任何正式 state降到 `n<5`、總 N<40或 basis<10即
  `NOT_EVALUABLE`。fixed point後仍有 invalid pair distance則整個 unit
  fail closed，不作資料驅動 read peel。
- **Reason**: 紅隊重現舊邏輯可將 AA n=5刪成 n=4後繼續檢定，或保留 final
  state coverage 7/9=0.7778的 CpG。
- **Verification**: 新增 n5 boundary、7/9 coverage與 monomorphic-stratum
  regression；analyzer + merge 36/36、研究工具全套72/72 PASS。
- **Evidence tier**: L1 engineering

### [2026-07-27 12:26] Adaptive輸出只能作 provisional merge input

- **Status**: Implemented before formal methyl run
- **Decision**: `--unit-key-file`輸出不計 subset-only BY/Holm，assessment標
  `PROVISIONAL_REFINEMENT`；summary固定標
  `PROVISIONAL_SUBSET_REFINEMENT_REQUIRES_FULL_FAMILY_MERGE`。merge必須
  核對兩個完整 family的 source/config契約，並在完整 family內重算BY/Holm。
- **Reason**: 合成100-hypothesis例中，同一 raw p=0.001的 subset q=0.001，
  但完整BY q=0.51873775。
- **Evidence tier**: L1 engineering

### [2026-07-27 12:26] Topology關係一律改為無方向 association overlay

- **Status**: Implemented before final HTML
- **Decision**: 報告只顯示 state pair是否屬於 frozen unanimous
  Hamming-1 edge；不使用 `first→second`箭頭。report與sidecar另逐欄核對
  census/evidence/detail/catalog marker set及所有 unordered pair effects。
- **Reason**: best pair排序不是 parent/child方向，固定用字串順序畫箭頭會反轉
  authoritative topology並暗示甲基決定 lineage。
- **Evidence tier**: L1 engineering

### [2026-07-27] 全 7-dataset formal run 完成

- **Status**: VALIDATED
- **Input scope**: 7 technical datasets / 6 biological samples / chr1–22 /
  exact PS × exact raw HP。
- **Census**: 21,601,606 sparse rows、8,893,098 candidate read projections、
  154,132 pattern strata、1,045 formal n5 units。
- **Formal result**: 811 evaluable；3 robust、627 evaluable no robust、
  181 confound-gate failures、234 not evaluable。
- **Direct answer**: formal full-four = 0；三個 robust 均為 k≥3；exact
  二位點 RA 對 AA/AR/RR 的 robust 均為 0。
- **Robust loci**: H1437 chr22 `AARR/RRAR/RRRA`、H2009 chr3
  `AAA/RAA`、HCC1937 chr10 `AAR/RRA`。
- **Evidence tier**: L2 statistical association；lineage/clone claim維持 L5。

### [2026-07-27] Standalone HTML 與 sidecar 完成 hash-bound 交付

- **Status**: VALIDATED
- **HTML**:
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/20260727_multisite_pattern_methyl_association_06.html`
- **Sidecar receipt**:
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/sidecar/pattern_methyl_sidecar.v1.receipt.json`
- **Report behavior**: aggregate / canonical / extreme / explained 四層、
  exact raw-HP filter、state profile heatmap、Bernoulli pair matrix、無方向
  topology relation、兩個指定 sentinel。
- **Independent audit fixes**: active-marker 以 genomic position 對 CpG
  ordinal 軸分段內插；最強 secondary 固定卡明示 no robust；COLO829
  zero-formal 可篩選；固定母體指標全部標「全域」。
- **Final closeout audit**: P0=0、P1=0、P2=0；全域指標歧義已修正並加入
  browser regression。
- **Evidence tier**: L1 engineering。

### [2026-07-27] 舊 BERNOULLI 樣式的 read 聚集圖加入正式 HTML

- **Status**: VALIDATED
- **Decision**: 對 source 保存合法原始 read × read Bernoulli matrix 的
  557 個 formal units，以只使用 distance 的 deterministic UPGMA
  average-linkage 排序；完整 R/A pattern、strand 與 RG 只在排序完成後疊加。
- **Coverage**: 29,316 reads、1,648,104 matrix cells。488 個不可畫單元
  分別為 234 個未通過正式可評估條件、253 個未達 detail trigger，以及
  1 個 N=169 超過 source `N≤160` 保存上限。
- **Encoding / privacy**: 矩陣在固定排序後量化為 uint8 row-major base64；
  JSON/HTML 不輸出 qname digest 或原始 read identifier，互動只顯示
  `r001` 類匿名 ordinal；RG 僅保留 `RG 1 / RG 2 / … / RG .` 類別。
- **Tie policy**: exact-distance ties 以 source leaf ordinal deterministic
  決勝；因此排序不使用 pattern/strand/RG，但不宣稱對任意 row permutation
  完全 invariant。
- **Claim boundary**: dendrogram 與 heatmap 是描述性 read-distance
  視覺；不設定 cluster cut，不把塊狀關係稱為 clone、祖源或 lineage。
- **Verification**: report unit tests 18/18；Chromium desktop/mobile/print
  PASS，console/page errors 0；canvas bytes/pixels、匿名 hover、鍵盤讀值、
  手機局部橫向捲動與 print visibility 均通過。
- **Output**:
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/20260727_multisite_pattern_methyl_association_06.html`
- **Evidence tier**: L1 engineering；L2 descriptive association。

## 🟠 偏離之處

### [2026-07-27 09:37] Canonical workstation完整重建暫不執行

- **Status**: Accepted
- **規範要求**: builder SoT重建 index + 7 sample pages後才算正式整合。
- **實作偏離**: 先產生 versioned pattern-methyl sidecar/report與 builder-ready
  contract，不繞過 authority。
- **理由**: candidate receipt綁定 LongLineage header 1054 bytes/SHA `1d19…`，
  現檔為1976 bytes/SHA `79d4…`；現有 unit test fail closed。
- **風險評估**: canonical sample-page node/edge halo需等 upstream authority重新
  簽定；本輪仍可完成資料、統計、standalone HTML與安全 adapter。
- **回退路徑**: upstream solver凍結並重簽 candidate factorization後，跑完整
  builder + 16-page browser QA。
- **Revisit if**: 新 candidate authority receipt可用。
- **Evidence tier**: L1

## 🟡 折衷考量

### [2026-07-27 09:37] Association census GO、lineage upgrade NO-GO

- **Status**: Accepted
- **方案 A**: We will publish unit-level association with strict confound gates.
- **方案 B**: 以顯著甲基差直接宣稱 clone/lineage；拒絕，證據層級不足。
- **方案 C**: 等第二 phaser再做所有分析；延後作 causal/lineage升級。
- **採用 A 理由**: 能回答使用者「是否有明顯位點」，且不越過現有負結果。
- **Revisit if**: marker-held-out/orthogonal phasing或cellular validation完成。
- **Evidence tier**: L2

### [2026-07-27 09:37] Canonical優先、attested fallback

- **Status**: Superseded for primary；reference sensitivity only
- **方案 A**: We will resolve each marker tile through a hash-bound source catalog,
  preferring canonical but allowing verified archive artifacts.
- **方案 B**: hardcode canonical-only；拒絕，已知 H2009 full4 marker缺 artifact。
- **方案 C**: 重掃 BAM；延後，現有 artifact足以先完成且避免大型重算。
- **Revisit if**: artifact catalog coverage低於 frozen candidate需求。
- **Evidence tier**: L1

## 🔴 未決問題

### [2026-07-27 09:37] LongLineage authority何時可重簽

- **Status**: open
- **Question**: 現行 solver header變動是否已定稿並可建立新 candidate authority？
- **Default if no answer**: 不重簽、不繞過；交付 standalone report +
  builder-ready sidecar。
- **Priority**: major
- **Evidence tier**: L5

## 📚 Lore

### [2026-07-27] 舊 chr22圖的 methyl matrix正確但 HP label過期

- **Constraint**: archive與 canonical methylation.csv SHA相同，但 reads.tsv
  對同 reads標成 `2-1` vs strict sparse `1-1`。
- **Why it matters**: 圖像樣式可沿用，HP/pattern標籤必須重新綁定。
- **Evidence**: 2026-07-27 input-map audit。

### [2026-07-27] HCC1395與DORADO不是 biological n=2

- **Constraint**: 兩者是同一 biological cell line的technical datasets。
- **Why it matters**: cross-sample consistency分母固定6，不得灌水。
- **Evidence**: exact-PS workstation cohort contract。

## ✅ Final verification

- Artifact catalog：2,313/2,313 markers PASS。
- Bernoulli parity：2,313/2,313 markers、277,560 pair cells PASS；
  invalid mismatch 0；max error `1.0538×10⁻⁴ ≤ 5×10⁻⁴`。
- Sidecar：1,045 bundles、811 evaluable、513 edge halos、288 pair bands、
  24 matrices；receipt PASS。
- Research Python tests：88/88 PASS。
- Combined-run merge tests：21/21 PASS。
- C++ full `./build/bin/run_tests`：258/258 PASS。
- Report tests：18/18 PASS。
- Browser QA：desktop/mobile/print PASS；console/page errors 0；
  embedded/external JSON equal；COLO829 zero filter下全域指標標籤 PASS；
  read-cluster bytes/pixels、匿名 hover、鍵盤與局部橫向捲動 PASS。
- Final HTML SHA-256：
  `cdc2e55165f05c430177d02993ee52535919ba4bda99a7036c23e31961361b1d`。
- Final data JSON SHA-256：
  `98a980c9fa28ed9859fe9656f3ab089ccfddd7eb907bd3580947ae375d56e3fc`。
- Final sidecar SHA-256：
  `45488a5fd2bdbbe5155a38aeaea3faa19ed02fa49a15be3fdb6139b85bf243f1`。
- Formal report：
  `InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/20260727_多sSNV_pattern與甲基關聯全面驗證_01.md`。

## Provenance

- Commit: `387a101e6a3292e0d7f230ba8d20271c7434972a`
- Created: 2026-07-27 09:37 +0800
- Completed: 2026-07-27
- Status: validated
