<!--
建立時間: 2026-07-17 04:43 +08:00
目標: 記錄 M2 frozen v2 正式 HCC1395_DORADO chr6 pilot 的可驗證終端結果，釐清 partial-read group-Steiner 實際處理，並界定哪些數字可用、哪些數字禁止當成最終研究結果
處理範圍: Task Type B Comprehensive validation 的第一個 single-chromosome gate；HCC1395_DORADO chr6 extraction 與 bootstrap=0 ranking；不外推 chr1-22 × 7 technical datasets
關聯檔案: InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/implementation-notes.md；InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/partial-read-method-audit.md；InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2正式執行Runbook_02.md
-->

# M2 frozen v2 正式 pilot：NO-GO 驗證紀錄

用 SCQA（Situation → Complication → Question → Answer）：先分開「方法是否合理」與「這次正式數據是否完成」，再逐層列出可用證據、不可用數字與下一版 exact-first 修復。

> **TL;DR：partial read 的正式 group-constraint＋likelihood 方法合理，且比逐 completion 建樹更正確、更有效率；但第一個正式 HCC1395_DORADO chr6 ranking pilot 在 8 小時上限被停止，缺少完成 receipt，fail-closed verifier 判定 NO-GO。因此本次只有抽取漏斗可作有效 pilot 證據；截斷 ranking 只能診斷效能，不能提供最終候選樹、拓撲、clone/subclone 或七資料集比例。**（影響：高；方法判讀信心：高；全量數據結論：尚不可用）

> **PARTIAL／NO-GO ribbon**：這不是教授版 FINAL 全量結果。正式範圍原訂 chr1–22 × 7 technical datasets／6 biological samples；本文件只記錄第一個 gate `HCC1395_DORADO × chr6`。HCC1395 與 HCC1395_DORADO 是同一 biological sample 的兩個 technical datasets。

## 1. 本次回答的四個問題

1. **正式流程有沒有照 frozen 計畫執行？** 有；PRE、freeze、verify-only、resource gate與 extraction皆通過，ranking依凍結參數執行至8小時停止。
2. **最後 gate 結果是什麼？** `NO_GO`；原因是 ranking沒有完成 `receipt.json`，不是把部分輸出視為PASS。
3. **partial read 是否會展開成多種樹，任一成功就保留？** 不會；實際方法是一條partial pattern形成一個subcube group constraint，全部groups由同一個候選state set聯合滿足，再保留全部全域最小的vertex sets。
4. **現在能否列出最終拓撲與七資料集數字？** 不能；B0未完成，B20、第二pilot、154-task full與final numeric summary均未啟動。

本任務服務：G3（read-level mutation-state evidence）、G4（多資料集可重現性）、G5（外部可驗證的fail-closed證據鏈）。

## 2. 正式執行的輸入、命令、輸出與終端結果

### 2.1 Frozen release 與來源身分

| 項目 | Frozen SHA-256 |
|---|---|
| M2 extractor | `5f910c836fd5ebaa3d5c933cfa6f0a97e36c7e4b4c38f0b206a4343e5e2d913a` |
| M2 ranker | `c82210f25870d1405cc070c18096fb7d1c2b908fb4d8a7858aece7ac4b151d28` |
| Hypercube exact solver | `9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95` |
| Runbook v2 | `8a347a67d1a6f59b25b6ea8d56d812ec2a084cbe4f5816d64d0431f03b47ee96` |

**輸入根目錄**：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
20260716_read_linked_hypercube_exact_likelihood_validation/m2_frozen_release_v2/
```

**Pilot單位**：`HCC1395_DORADO × chr6`。

### 2.2 Extraction command

```bash
timeout --signal=TERM --kill-after=5m 8h nice -n 10 ionice -c2 -n7 \
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
.../scripts/extract_lossless_read_linkage.py \
  --manifest .../release_contract/input_contract/canonical_manifest.json \
  --sample HCC1395_DORADO \
  --chrom chr6 \
  --output-dir .../pilots/HCC1395_DORADO_chr6/extraction \
  --mapq-min 20 \
  --baseq-min 20 \
  --bridge-thresholds 1,2,3,5 \
  --samtools-threads 1 \
  --require-existing-empty-output-dir
```

**輸出**：`.../pilots/HCC1395_DORADO_chr6/extraction/`。

**實際結果**：exit=`0`、wall=`7:04.43`、CPU=`135%`、max RSS=`754,028 KiB`；`receipt.json` 的 `all_pass=true`，receipt SHA-256=`bcfc8c1346711ca6f524350848a1c7c4a008aa462cb821a67b44ba84654a2f7c`，外部sidecar重驗為`OK`。

### 2.3 Ranking B0 command

```bash
timeout --signal=TERM --kill-after=5m 8h nice -n 10 ionice -c2 -n7 \
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
.../scripts/build_m2_patterns_and_rank.py \
  --input-dir .../pilots/HCC1395_DORADO_chr6/extraction \
  --output-dir .../pilots/HCC1395_DORADO_chr6/ranking_bootstrap0 \
  --threshold 1 --threshold 2 --threshold 3 --threshold 5 \
  --component-basis PS_HP1 --component-basis PS_HP2 \
  --hp-family 1 --hp-family 2 \
  --structural-exact-pattern-minread-grid 1,2,3,5 \
  --primary-structural-exact-pattern-minread 3 \
  --exact-k-max 12 \
  --max-vertex-sets 256 \
  --solver-time-limit-seconds 30 \
  --fixed-error-grid 0.005,0.01,0.02,0.05 \
  --minimum-bq-error-rate 0.000001 \
  --maximum-bq-error-rate 0.25 \
  --conditional-candidate-ranking-bootstrap-replicates 0 \
  --conditional-candidate-ranking-bootstrap-seed 20260716 \
  --tie-tolerance 0.000001 \
  --require-existing-empty-output-dir
```

**輸出**：`.../pilots/HCC1395_DORADO_chr6/ranking_bootstrap0/`。

**實際終端輸出**：

```text
Command exited with non-zero status 124
Elapsed (wall clock) time: 8:00:00
User time: 29267.25 s
System time: 899.07 s
CPU: 104%
Maximum resident set size: 510208 KiB
Exit status: 124
```

沒有產生 `ranking_bootstrap0/receipt.json`；四個 `.tsv.gz` 均為截斷串流，`gzip -t` 回傳`1`與`unexpected end of file`。

### 2.4 Fail-closed verifier

**輸出**：

```text
.../pilots/HCC1395_DORADO_chr6/pilot_gate_verification_receipt.json
```

**實際結果**：

```json
{
  "all_pass": false,
  "release_gate": {"verdict": "NO_GO", "all_metrics_go": false},
  "failure": "PilotVerificationError: missing JSON: .../ranking_bootstrap0/receipt.json",
  "verification_independence": {
    "reads_bam": false,
    "reads_vcf": false,
    "launches_pilot_or_samtools": false
  }
}
```

Verifier receipt SHA-256=`dc1381a517ba70b0a3f4a0a90b8dbf143ff4b515758a631e3a4ede7fbe456926`；sidecar在receipt所在目錄執行`sha256sum -c`後為`OK`。

## 3. 哪些數據有效、哪些只能診斷

| 證據層 | 狀態 | 可否引用為研究結果 | 原因 |
|---|---|---|---|
| Frozen來源與PRE／verify-only | PASS | 可以，作provenance | source SHA、input identity與release contract已綁定 |
| HCC1395_DORADO chr6 extraction | PASS | 可以，但只代表單一pilot | 完整receipt、守恆、sidecar與exit 0 |
| Ranking B0四個gzip | TRUNCATED | 不可以；只可作效能診斷 | 8 h timeout、gzip不完整、無receipt |
| B20 bootstrap | NOT STARTED | 不可以 | B0 gate未過 |
| H2009 chr2 stress pilot | NOT STARTED | 不可以 | 第一pilot已NO-GO |
| 154-task full extraction/ranking | NOT STARTED | 不可以 | pilot gate阻擋 |
| NO-GO pilot HTML | PASS | 可以，但只能報告本次gate與有效extraction | official artifact validation／package／desktop＋mobile browser verification皆通過 |
| final numeric summary／FINAL topology HTML | NOT CREATED | 不可以 | 沒有全量verified receipts |

## 4. 有效的 HCC1395_DORADO chr6 extraction 漏斗

### 4.1 Alignment／molecule funnel

| 步驟 | 數量 | 相對前層 | 相對raw | 守恆解釋 |
|---|---:|---:|---:|---|
| Raw overlapping alignments | 570,915 | 100.0000% | 100.0000% | 分母 |
| Primary alignments | 521,714 | 91.3821% of raw | 91.3821% | raw＝primary＋secondary＋supplementary |
| Secondary | 5,522 | 0.9672% of raw | 0.9672% | 與supplementary共同排除 |
| Supplementary | 43,679 | 7.6507% of raw | 7.6507% | 與secondary共同排除 |
| Excluded by flag | 49,201 | 8.6179% of raw | 8.6179% | `5,522＋43,679＝49,201` |
| MAPQ rejected after flag | 451 | 0.0864% of primary | 0.0790% | primary後再篩 |
| Canonical eligible alignments | 521,263 | 99.9136% of primary | 91.3031% | `521,714−451` |
| Exact sidecar matches | 521,263 | 100.0000% of eligible | 91.3031% | missing sidecar＝0 |
| Molecule sparse rows | 521,263 | 100.0000% of eligible | 91.3031% | unique molecule IDs＝521,263 |

這個漏斗只回答「抽取與join是否守恆」，不回答拓撲是否唯一。

### 4.2 Site-call evidence

| 單位 | 數量 | 分母與比例 | 意義 |
|---|---:|---:|---|
| Sparse site-call rows | 2,528,172 | 100.0000% | molecule×target-site呼叫總列數 |
| Fixed R/A calls | 2,211,241 | 87.4640% of site-call rows | 可判定REF或ALT |
| 其他O/D/S/L/X類 | 316,931 | 12.5360% of site-call rows | 缺失、非R/A或不可用資訊，保留原因而非硬補值 |
| ALT calls | 1,499,656 | 67.8197% of fixed R/A | 位點投影中的ALT呼叫；不是獨立clone數 |
| REF calls | 711,585 | 32.1803% of fixed R/A | `fixed R/A−ALT` |

### 4.3 HP1／HP2 exact-PS contract

| 單位 | Known PS | Missing PS | Known比例 | Missing比例 |
|---|---:|---:|---:|---:|
| HP1/HP2 molecule rows | 74,573 | 6,182 | 92.3447% | 7.6553% |
| HP1/HP2 fixed R/A calls | 89,925 | 15,395 | 85.3826% | 14.6174% |

Primary inference只使用known exact PS；missing PS保留為diagnostic，不混入primary component。

### 4.4 Read-linkage component units

`bridge threshold`是相鄰位點切口至少需要幾個unique molecules同時在兩側有固定R/A資訊；它與後續`structural exact-pattern minread`是兩個不同門檻。

| Bridge threshold | PS_HP1 components | PS_HP2 components | 合計 | k=1 | k>1 |
|---:|---:|---:|---:|---:|---:|
| 1 | 798 | 765 | 1,563 | 1,223 | 340 |
| 2 | 874 | 864 | 1,738 | 1,447 | 291 |
| 3 | 913 | 911 | 1,824 | 1,560 | 264 |
| 5 | 966 | 967 | 1,933 | 1,702 | 231 |
| Grid合計 | 3,551 | 3,507 | 7,058 | 5,932 | 1,126 |

注意：`7,058`是`component basis × bridge threshold × component`的**參數格點unit instances**，不是7,058個互斥生物區域。以bridge threshold=3單獨看，1,824個component instances中，k=1為1,560（85.5263%），k>1為264（14.4737%）。

## 5. Ranking 截斷進度：只可解釋「為何超時」

每個structural-minread設定都要處理上述7,058個unit instances；四個設定共規劃28,232次unit-evaluations。

| Structural exact-pattern minread | 可讀runtime rows | 規劃 | 可讀進度 | 科學資格 |
|---:|---:|---:|---:|---|
| 1 | 7,058 | 7,058 | 100.0000% | 截斷檔內診斷，不是sealed result |
| 2 | 7,058 | 7,058 | 100.0000% | 同上 |
| 3（primary） | 5,154 | 7,058 | 73.0235% | 未完成 |
| 5 | 0 | 7,058 | 0.0000% | 未開始 |
| 合計 | 19,270 | 28,232 | 68.2559% | 不可作拓撲分母 |

截斷串流另可讀：symbolic pattern rows=`433,302`、candidate rows=`20,233`、winner-responsibility rows=`73,426`。這些數字只表示程式停止前寫入多少列；不代表完整候選樹數、完整winner數或區域比例。

### 5.1 確認的長尾瓶頸

minread=1的7,058個runtime rows中：

| 狀態 | Units | Candidate generation time | Likelihood time |
|---|---:|---:|---:|
| Candidate＋likelihood皆invoked | 4,852 | 588.655131 s | 338.985172 s |
| Candidate invoked、likelihood未invoked | 33 | 14,715.729568 s | 0 s |
| Solver未invoked | 2,173 | 0 s | 0 s |

這33個長尾units占minread=1 candidate-generation時間的`96.153683%`。它們都有partial structural patterns，但正確因果不是「程式把所有partial completions各跑一棵樹」；而是partial ambiguity可間接產生很多同分的最優vertex sets，現行程式為每個新set重新建立MILP並再呼叫solver。

`solver-time-limit-seconds=30`是**每次MILP solve**的上限，不是每個unit總上限。一個unit可先求`h*`，再反覆加入no-good cut、重建模型、重新求下一個最優set，最多保存256 sets，因此單unit可累積遠超過30秒。

## 6. Partial read 的正確數學與實際程式流程

### 6.1 一條partial pattern有多少種可能？

對k個sSNV，完整mutation state為：

\[
v\in\{R,A\}^{k}.
\]

一條partial pattern為：

\[
p\in\{R,A,X\}^{k},
\]

其中`X`表示該位置沒有可用R/A觀測，不表示50% R＋50% A，也不等於另一條read。

若有`u`個X，在完整k維立方體中概念上有：

\[
2^u
\]

個相容完整states。現行active observed-ALT universe中，真正進solver的有效數量是：

\[
2^{u_{eff}},\qquad
u_{eff}=\operatorname{popcount}(free\_mask\;\&\;active\_mask).
\]

### 6.2 實際不建立 `2^u` 棵樹

每個達`structural_exact_pattern_minread`的distinct exact R/A/X pattern形成一個群組：

\[
G(p)=\{v:((v\oplus alt\_mask)\;\&\;fixed\_mask)=0\}.
\]

候選vertex set `N`必須滿足：

\[
N\cap G(p)\neq\varnothing.
\]

所有retained groups都在**同一個MILP中同時滿足**。程式會把一個group相容的active vertex indices寫成一條sparse hit row，但不把每個completion當成獨立tree world，也不建立跨reads的completion Cartesian product。

### 6.3 Structural objective

令`F`為完整R/A patterns強制的observed states，root為全REF state。現行objective為：

\[
\min_N\;|N\setminus(F\cup\{root\})|.
\]

限制包括：

1. root與全部full-observed states必須在`N`；
2. 每個被選非根vertex至少有一個Hamming-distance=1的已選predecessor；
3. 每個partial group至少由`N`中的一個vertex命中。

先求全域最小extra-state數`h*`，再固定objective=`h*`並加入no-good cuts，列出所有不同的全域最小vertex sets。欄名`minimum_hidden_nodes`更精確的意義是minimum-extra-state count；它包含partial representative與connector，不等於真實hidden clone數。

### 6.4 為何「任一completion成功就保留」不正確？

例：full patterns=`RRA, AAA`，partial=`RAX`。

```text
G(RAX) = {RAR, RAA}
```

若強迫`RAR`，最少需要2個non-full extra vertices；若選`RAA`，只需要1個。若剛好先測到`RAR`並在「可行」時停止，就會保留非全域最小解。

如果要用explicit completion方式得到與group constraint等價的結果，必須：

1. 對全部retained partial groups聯合枚舉completion assignment；
2. 每一個joint assignment都求minimum-extra tree；
3. 跨全部worlds比較同一個全域`h*`；
4. 只保留所有達全域`h*`的vertex sets；
5. 對重複vertex sets去重。

若第i個group有`u_i`個X，最壞需要：

\[
2^{\sum_i u_i}
\]

個joint worlds。因此現行symbolic group constraint比explicit world enumeration更有效率，也避免first-success偏差。

## 7. Likelihood、VAF與read數量各自扮演什麼角色？

### 7.1 Read-pattern likelihood

對候選vertex set `N`，每個state有混合比例`π_v`。對每一條pattern，只在fixed R/A位置計算BQ/error-aware emission，X作missing marginalization：

\[
\ell(N)=\max_{\pi\in\Delta}
\sum_p n_p\log\left(\sum_{v\in N}\pi_vP(p\mid v,BQ,error)\right).
\]

這個式子同時使用：

- 每種read pattern的數量`n_p`；
- 哪些位點在同一條molecule上共同被觀測；
- 每個固定R/A call的BQ／錯誤率；
- partial X對可能latent states的marginalization。

低於structural minread的informative molecules不形成hard structural constraint，但仍可進likelihood。

### 7.2 VAF不是第二份獨立證據

由同一批reads計算的VAF，本質上已包含在pattern counts與估計`π`中。若再把同read VAF加成另一個分數，會double counting。VAF適合：

- 把估計state mixture轉成教授可理解的比例摘要；
- 與獨立caller AF、CN/LOH或matched-normal資訊做敏感度分析；
- 檢查祖先高於後代的敘述是否與外部模型一致。

但它不應在同read likelihood後再重複加權。

### 7.3 為何不直接用read數量修改Hypercube edge weight？

一條snapshot read支持的是「某個mutation state或partial subcube」，不是「某次祖先→後代transition真的發生」。把read數量分攤到流入邊會：

1. 把state abundance誤當成edge transition evidence；
2. 對同一read在結構與edge score重複計數；
3. 任意分攤時產生模型依賴的pseudo-evidence；
4. 仍無法區分相同vertex set下的不同parent edges。

所以較合理順序是：先用unit Hamming edges定義mutation相容性與minimum-extra候選，再用read-pattern likelihood排序不同vertex sets。相同vertex set的不同edges目前應標`EDGE_NONIDENTIFIABLE`，除非加入真正edge-informative的獨立時間、單細胞、跨樣本lineage或其他資料。

## 8. 方法結論：合理，但現行枚舉器不夠可擴充

### 8.1 已合理的部分

- Partial read是subcube group，不是硬補值read。
- 所有groups由同一個candidate set聯合滿足。
- 先求minimum-extra，再保留全部同分optimal vertex sets。
- Likelihood對X做marginalization並使用BQ/error，而非選一個completion。
- Read counts與同read VAF只使用一次。
- 無法辨識的edge不假裝唯一。

### 8.2 真正需要修正的部分

| 優先級 | Exact-first修復 | 是否改變primary estimand | 驗證方式 |
|---|---|---|---|
| P0 | Canonical structural-problem cache：相同full/group signature跨threshold/minread重用 | 否 | candidate set、`h*`、complete flag與semantic SHA逐案相同 |
| P0 | Persistent solver＋incremental no-good cuts，不為每個set重建MILP | 否 | k≤6 exhaustive oracle；33個長尾unit frozen regression |
| P0 | 真正per-unit總deadline；超時輸出`complete=false`與零candidate | 不改候選定義，但增加abstain | fault test確保超時不產生winner |
| P1 | Atomic per-unit checkpoint/resume＋deterministic merge | 否 | 中斷續跑的final semantic SHA等於不中斷run |
| P1 | Shared emission matrix與fit cache | 否 | LL、KKT gap、winner/tie逐列一致 |
| P1 | B20只讀已驗B0候選並增加bootstrap | 否 | B20不重跑candidate generation；綁定B0 receipt |
| P1 | Deterministic unit sharding | 否 | 不同worker數輸出semantic SHA一致 |
| Contract change | 只對primary minread=3做完整ranking；其他minread只做structural funnel | Primary不變、sensitivity deliverable改變 | 必須新release與明示`STRUCTURAL_ONLY_NOT_RANKED_BY_DESIGN` |

不可以把completion抽樣、first-success、beam top-N、硬補X、任意k-cut、read/VAF edge cost或未證明heuristic pruning包裝成exact結果；這些會改變candidate universe或estimand。

## 9. 現在能下的結論與不能下的結論

### 可以下的結論

1. HCC1395_DORADO chr6 extraction流程在本frozen contract下守恆且receipt PASS。
2. Partial-read symbolic group formulation與error-aware likelihood的科學語意合理。
3. 目前正式ranking實作不具足夠可擴充性；第一pilot在8小時gate失敗。
4. 長尾主要發生在exact candidate enumeration反覆重建MILP，不是explicit completion-world enumeration。
5. 下一版應先做exact-preserving solver/cache/checkpoint修復，再重跑pilot。

### 不能下的結論

- 不能報告七資料集最終區域數、T、Topo、唯一第一／並列第一比例。
- 不能用20,233個截斷candidate rows當作完整T或樹數。
- 不能由minimum-extra count宣稱真實clone數。
- 不能由同一vertex set的likelihood宣稱唯一parent-edge topology。
- 不能由單一chr6 pilot外推所有樣本一致性。
- 不能把舊的`98,941／48,963／39,885／29,053`或其他舊口徑數字混入本次M2 final。

## 10. 下一個release的 Step → Verify 計畫

1. **建立長尾frozen regression panel**：保存33個slow/incomplete units與complete controls。  
   → 驗證：每個unit有輸入semantic hash、舊`h*`／candidate status與預期abstain界線。
2. **實作canonical structural cache＋persistent exact enumerator**。  
   → 驗證：k≤6 brute-force完整candidate set 0 mismatch；33個長尾units不丟失任何已證明optimal set。
3. **加入per-unit總deadline與原子checkpoint**。  
   → 驗證：deadline後`complete=false`、candidate rows=0；中斷resume與uninterrupted final semantic SHA一致。
4. **拆分B0與B20、共用likelihood fit**。  
   → 驗證：B20只增加bootstrap；primary LL／winner／tie與未拆分參考逐列一致。
5. **重新pre-decision audit、freeze與HCC1395_DORADO chr6 pilot**。  
   → 驗證：wall≤4 h才為GO；incomplete rate≤1%；exact-limit coverage≥90%；receipt與sidecar PASS。
6. **再執行H2009 chr2 tail-stress pilot**。  
   → 驗證：同一gate全部PASS。
7. **才啟動154 extraction＋154 ranking tasks與final verifier**。  
   → 驗證：154/154 child receipts、POST=PRE、cross-stage binding與final numeric summary全部PASS。
8. **最後才生成教授版FINAL HTML**。  
   → 驗證：每個數字具有資料範圍、unit、總分母、相對分母、比例、abstain原因與source receipt。

## 11. 已交付的NO-GO HTML與官方QA

教授版短報告不是FINAL topology報告，而是本次正式gate的blocked／NO-GO證據摘要：

- HTML：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2_pilot_NO_GO驗證報告_01.html`
- Canonical artifact：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/html_no_go_report_v1/artifact.json`
- Delivery receipt：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/html_no_go_report_v1/delivery_receipt.json`

| HTML區塊 | 視覺／內容 | 資料來源 | 讀者要理解的結論 |
|---|---|---|---|
| 正式結論 | Answer-first文字 | 本驗證紀錄＋verifier receipt | Extraction可用；ranking與T／Topo不可用 |
| 有效抽取漏斗 | Native funnel chart | extraction `receipt.json` | 570,915 raw → 521,263 exact-matched canonical molecules |
| Partial-read方法 | 結論區內的精簡公式 | method audit＋frozen source | `2^u`是概念states；實際用joint group constraint，不做first-success |
| 下一步 | Exact-first修復摘要 | scaling audit＋source inspection | cache／persistent solver／deadline／checkpoint後重跑雙pilot |

Official delivery receipt：`validation=passed`、`package=passed`、`verification=passed`；source dialog與鍵盤互動通過，viewports=`1440, 390`，rendered counts=`2 blocks／1 chart`。完整命令、比例、截斷診斷與主張上限以本Markdown為權威細節來源。

## 12. 證據路徑

- Frozen run root：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_frozen_release_v2/`
- Extraction receipt：`.../pilots/HCC1395_DORADO_chr6/extraction/receipt.json`
- Extraction time：`.../pilots/HCC1395_DORADO_chr6/extraction.time.txt`
- Ranking time：`.../pilots/HCC1395_DORADO_chr6/ranking_bootstrap0.time.txt`
- Fail-closed verifier receipt：`.../pilots/HCC1395_DORADO_chr6/pilot_gate_verification_receipt.json`
- Attempt logs：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_frozen_release_v2_attempts/20260716T122154.841147978Z-2098709-14255/`
- Method audit：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/partial-read-method-audit.md`
- Living notes：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/implementation-notes.md`
- NO-GO professor HTML：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2_pilot_NO_GO驗證報告_01.html`
- HTML QA receipt：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/html_no_go_report_v1/delivery_receipt.json`

## 13. 最終判斷

本次「最後數據結果」是**流程gate的終端結果**，不是生物拓撲結果：

```text
Extraction: PASS
Ranking B0: TIMEOUT at 8:00:00, exit 124
Independent pilot verifier: NO_GO
B20 / second pilot / full 154 / final summary: NOT STARTED
Professor FINAL topology report: NOT ELIGIBLE
```

方法方面，partial read應保留為joint subcube groups，先求全部minimum-extra candidate vertex sets，再用不重複計數VAF的read-pattern likelihood判斷哪些候選能被資料區分；不應逐completion first-success，也不應把read abundance直接當edge transition weight。工程方面，下一版必須把「相同結構重算＋逐set重建MILP」改成cache、persistent solver、per-unit deadline與checkpoint，通過雙pilot後才值得展開全量。
