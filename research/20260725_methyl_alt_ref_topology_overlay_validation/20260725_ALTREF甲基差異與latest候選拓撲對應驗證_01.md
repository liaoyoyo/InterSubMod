<!--
建立時間: 2026-07-25
目標: 驗證 focal ALT/REF、linked-partner REF/ALT 甲基差異能否對應最新 exact-PS 候選拓撲
處理範圍: 7 個 formal G1 pair、7 個 focal ALT/REF controls、10 個 exact-PS×HP lanes、HCC1395/DORADO 同位點比較
關聯檔案: artifact.json; results/validation_receipt.json; results/html_delivery_receipt.json; 20260725_ALTREF甲基差異與latest候選拓撲對應驗證_01.html
-->

# ALT/REF 甲基差異 × 最新 exact-PS 候選拓撲對應驗證

**Task type B — comprehensive validation；服務 G3、G4。**

## TL;DR

可以把 ALT/REF 證據標到最新結構，但要分成兩種問題：

1. **Linked-partner REF/ALT 與甲基群的共分離：7/7 formal G1 pairs 成立。** 7/7 都落在同一 strict W，10/10 HP lanes 有 threshold-qualified direct read-linkage；638 reads 來自 7 條 formal-signal lanes，另有 152 reads 來自 3 條配對 RR-only background lanes，合計 790。
2. **Focal 位點本身 ALT vs REF 的甲基差異：只有 3/7 可檢定、1/7 通過 allele-axis gate。** HCC1395 focal 為 ALT=98、REF=2，不能檢定 focal ALT/REF methyl difference。
3. **最新候選結構：整棵 candidate T 唯一為 2/7；局部 focal→partner relation 可解為 3/7。** 第三個是 H2009 chr18：整樹有 6 個 AF-best ties，但 6/6 都在 focal-bearing state 後取得 partner。

因此樹上可標：

> `partner-ALT-bearing derived molecular state ↔ methyl-group enrichment`

但不可升級成「突變造成甲基變化」、「MG 就是 clone/subclone」或「真實 cellular lineage 已證明」。

## 兩種 ALT/REF 定義不能混用

| 問題 | Read population | 本次結果 | 可標記內容 |
|---|---|---:|---|
| G1 partner REF/ALT | focal-ALT methyl-core 且 partner callable | 7/7 formal positives | partner state 與 MG enrichment 關聯 |
| Focal ALT/REF control | focal ALT+REF 重新 joint clustering | testable 3/7；aligned 1/7 | focal allele-associated methyl badge |
| Strict edge | exact-PS×HP fixed endpoint molecules | 10/10 lanes；signal 638 + RR-only background 152 = 790 | 無向 read-linkage 與 RR/RA/AR/AA |
| C++ candidate T | whole-block supported patterns，weight≥3 | unique 2/7；tied 1/7；abstain 4/7 | model-conditional acquisition relation |

Focal ALT+REF control 是獨立重分群，沒有沿用 ALT-only M1 centroid，因此其群標不能直接當成原 MG 1-1／1-2。

## 七個 formal G1 pair 對位結果

| Dataset / focal→partner | Partner MG×allele | Focal ALT/REF | Full candidate T | Pair relation |
|---|---:|---|---|---|
| H2009 chr12:4,413,414→4,414,974 | V=0.929 | below effect / non-significant | resource abstain | unresolved |
| H2009 chr13:28,815,939→28,817,498 | V=0.847 | REF=0，不可檢定 | resource abstain | unresolved |
| H2009 chr18:567,920→563,687 | V=0.894 | REF=0，不可檢定 | 6 AF-best ties | **focal→partner，6/6** |
| H2009 chr4:2,307,521→2,304,921 | V=0.964 | **aligned；V=0.618，p=0.002** | resource abstain | unresolved |
| H2009 chr5:44,615,693→44,612,223 | V=0.614 | joint clustering unstable | resource abstain | unresolved |
| HCC1395 chr3:127,986,757→127,981,978 | V=0.645 | ALT=98、REF=2，不可檢定 | **unique；Sister＋direct** | **focal→partner，1/1** |
| HCC1954 chr8:100,071,382→100,070,832 | V=0.868 | V=0.245，未達效應門檻 | **unique；Direct-only** | **focal→partner，1/1** |

這是 407,738 個 evaluated rows 中已通過 formal gate 的 7 個陽性 pair；**7/7 不代表方法準確率 100%**。

## HCC1395 的具體對應

### 1. Focal ALT vs REF

- focal chr3:127,986,757 G>A：ALT=98、REF=2。
- REF 少於 joint gate 與 REF-only gate 所需數量，判定 `NOT_EVALUABLE`。
- 所以不能說 HCC1395 此 focal mutation 本身具有已確認的 ALT/REF 甲基差異。

### 2. Focal-ALT 內的 methyl group × linked-partner allele

Partner chr3:127,981,978 C>G：

| Methyl group | Partner REF | Partner ALT |
|---|---:|---:|
| MG 1-1 | 21 | 11 |
| MG 1-2 | 0 | 19 |

共 51 個 focal-ALT methyl-core、partner-callable reads；V=0.645、Δ ALT fraction=0.656、global BY q=1.08×10⁻⁴、conditional permutation p=0.001。

合理解釋是：**MG 1-2 對 partner-ALT 明顯 enrichment**。這是 read-level 共分離，不是因果。

### 3. Strict edge 與 candidate T

- exact PS=126975958、HP2、W k=5。
- focal-first strict states RR/RA/AR/AA = `2/0/21/32`，direct support=55。
- C++ 有 44 個 minimum vertex sets／44 棵 trees；supported-pattern read-AF 選出 1 個 winner。
- 唯一 candidate 包含：

```text
ROOT
└─ focal chr3:127,986,757 G>A
   ├─ sister branch
   └─ partner chr3:127,981,978 C>G
      └─ downstream branch
```

因此可以把 `partner-ALT enrichment` 疊在 focal→partner candidate branch，但它只是 pairwise projection；51 個 methyl reads 不能當成完整五位點 node membership。

### 4. HCC1395 與 DORADO

| 指標 | HCC1395 | HCC1395_DORADO |
|---|---:|---:|
| Strict W k | 5 | 2 |
| RR/RA/AR/AA | 2/21/0/32 | 1/23/0/26 |
| P(partner ALT \| focal ALT) | 60.38% | 53.06% |
| Supported-pattern partner AF | 54.43% | 54.08% |
| Candidate pair relation | focal→partner | focal→partner |
| Focal M1 methyl multigroup | K=2 stable | K=1，未重現 |

Strict conditional fraction 差 7.32 percentage points，Fisher p=0.549；未偵測差異，但不是等效性證明。可說 **genetic candidate substructure 跨技術重現**，不可說 methyl association 已重現。

## 為何 H2009 chr18 可算第三個「局部已解」

2026-07-25 exact factorization 的目標 row：

- active positions = `[563687, 567920, 581932, 632412]`
- AF-best trees = 6
- global-best edge incidence 包含 `14→15: 6/6`
- vertex 14 已含 focal bit（567,920），14→15 才取得 partner bit（563,687）

所以：

- Full T：不唯一，仍有 6 棵。
- Coarse topology：6 棵皆 Direct-only。
- Focal/partner 局部 relation：6/6 都是 focal before partner，故可標記。

## 執行與驗證

### 輸入

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity/methyl_ssnv_pair_results.tsv.gz`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel/all_ssnv_tumor_ref_control_site_results.tsv.gz`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260725_exact_ps_candidate_factorization/all7_v2/`

### 分析命令

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  research/20260725_methyl_alt_ref_topology_overlay_validation/scripts/build_alt_ref_topology_report.py
```

實際輸出：

```json
{"all_pass":true,"checks":32,"formal_pairs":7,"strict_lanes":10,
 "unique_candidate_pairs":2,"tied_candidate_pairs":1,
 "resource_abstain_pairs":4,"resolved_pair_relations":3,
 "focal_axis_aligned":1,"signal_direct_support":638,
 "background_direct_support":152,"hcc_dorado_order":"FOCAL_BEFORE_PARTNER"}
```

### HTML 封裝命令

```bash
CHROMIUM_EXECUTABLE_PATH=/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium_headless_shell-1223/chrome-headless-shell-linux64/chrome-headless-shell \
node research/20260725_methyl_alt_ref_topology_overlay_validation/scripts/deliver_portable_with_topbar_fix.mjs
```

實際驗證：

```text
blocks=32, charts=6, tables=6, metrics=5
viewports=1440,390
sourceDialog=passed
sourceInteraction=keyboard_menu_semantic_click
```

## Step → Verify

1. 鎖定 7 個 formal G1 rows → 驗證：formal pair count=7。
2. 以 coordinate 接 focal ALT/REF control → 驗證：7/7 有 control row；testable=3、aligned=1。
3. 以 exact PS×HP×component 接 strict W → 驗證：7/7 same W；10 lanes；signal support=638、RR-only background support=152、total=790。
4. 以 sample×region_id 接 C++ topology 與 exact factorization → 驗證：full unique=2、AF tied=1、resource abstain=4。
5. 對每棵 global-best tree 重算 focal/partner relation → 驗證：resolved=3；三者皆 focal→partner；H2009 chr18=6/6。
6. 產生 bounded SQLite 與 portable HTML → 驗證：12/12 SQL datasets、32/32 analysis checks、desktop/mobile browser QA PASS。

## 限制與 claim ceiling

- Read-AF 是 whole-pattern support-filtered ALT/(REF+ALT)，不是 caller VAF、CCF、posterior 或獨立驗證。
- 四個 formal pair 因 1,000-node guard fail-closed abstain；cohort `validation_evidence_eligible=false`。
- 未正式補 matched-normal、CN/LOH、purity、multiplicity、CCF 或 single-cell/colony truth。
- HCC1395 與 DORADO 是同一 cell line 的技術來源，不是獨立生物複本。
- Factorization current `all7_v2` receipt 已綁定目前 workspace source；current SHA 與 receipt-bound SHA 均為 `c52231…`，byte-identical source identity 成立。後續 renderer 仍須逐次 fail-closed 驗證 receipt、TSV 與 source identity。

## 輸出

- HTML：`InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/20260725_ALTREF甲基差異與latest候選拓撲對應驗證_01.html`
- Pair join：`InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/data/formal_pair_alt_ref_topology_join.tsv`
- HP lane join：`InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/data/strict_hp_lane_cpp_topology_join.tsv`
- Focal control：`InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/data/focal_alt_ref_control_join.tsv`
- HCC cross-source：`InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/data/hcc1395_crosssource_pair_substructure.tsv`
- Analysis receipt：`InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/results/validation_receipt.json`
- HTML receipt：`InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/results/html_delivery_receipt.json`
