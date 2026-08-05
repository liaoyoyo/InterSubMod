<!--
建立時間: 2026-07-15T12:00:00+08:00
目標: 全量驗證 focal-ALT 甲基多群是否與同一 read 的 sSNV 共現狀態相關，並界定單位點可推估的證據上限
處理範圍: 7 datasets / 6 biological samples；chr1-22 最新 LongPhase-S recalibrated FILTER=PASS biallelic sSNV
cycle_id: cycle_20260715-all-ssnv-focal-alt-cooccurrence
topic: all_ssnv_focal_alt_multigroup_cooccurrence_validation
status: verdict_GO_validation_claim_PROBE
audit_version: 0.1
關聯檔案:
  - InterSubMod/docs/CURRENT_FOCUS.md
  - InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/pre-decision-audit.md
  - InterSubMod/docs/experiments/in_progress/2026/07/20260706_methyl_vs_ssnv_lineage_01.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md
-->

# Pre-Decision Audit: 全 sSNV focal-ALT 甲基多群與共現關聯

> **雙判定**：執行全量、可否證的 association validation 為 **GO（70/100）**；把單位點甲基多群直接稱為 clone/subclone 仍為 **PROBE**。
> **研究目標**：G3 / G4 / G5；Task Type B（469,849 sites、7 datasets、chr1-22，不以 subset 代替 final）。

## 啟動 5 問

1. **Thread D**：是，核心是 read-level epigenetic latent structure。
2. **Thread B 撤回範圍**：否；本輪不把甲基當全域 FP filter，也不宣稱改善 caller F1。
3. **KDE-corrected**：不適用；本輪不使用 coverage multiple 或 KDE-normalized score 作主要 endpoint。
4. **VCF caller AF**：只作 post-hoc 描述，固定讀 LongPhase-S recalibrated PASS VCF 保留的 caller `AF/AD`，不使用 merged AF。
5. **高影響操作**：是，屬長計算與大輸出；不改寫或覆蓋 BAM，不寫 canonical，不刪除既有檔案。

## §0 Cynefin Domain Gate

- **Domain**: **Complex**。
- **理由**：共同 ancestral ALT 下可能混有 descendant genotype states，但 methyl group 也可由 cis-ASM、HP、CN/LOH、read quality、strand、CpG coverage 或 stochastic epigenetic state 生成。
- **策略**：先 truth-blind 建立 methyl groups，再以未參與分群的 sSNV read-level R/A pattern 作 post-hoc external label；禁止用共現標籤預先挑 K 或形成群。

## §1 Observation Completeness

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| 最新 LongPhase-S PASS autosomal biallelic sSNV 共 469,849 | ✓ | L1 | `InterSubMod/docs/CURRENT_FOCUS.md` |
| 最新 layered read-tag exact matches 11,513,224 / 11,513,224 | ✓ | L1 | `InterSubMod/research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md` |
| FP-only focal-ALT 全量已證明 stable methyl multigroup 存在，但無 subclone confirmation | ✓ | L1/L3 | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/` |
| 早期 region-lineage pilot 曾觀察 methyl 與 multi-sSNV genotype 差異 | ✓ | L3 | `InterSubMod/docs/experiments/in_progress/2026/07/20260706_methyl_vs_ssnv_lineage_01.md` |
| 全 469,849 sites focal-ALT 結果與明確不可評估原因 | □ | 本輪 P0 | 尚待執行 |
| methyl group x 同 read 鄰近 sSNV R/A/genotype association + global FDR | □ | 本輪 P1 | 尚待執行 |
| focal ancestral relation 在 complete candidate trees 的一致性 | □ | 本輪 P1 | 尚待執行 |
| single-cell / clone isolate / multi-region orthogonal cellular truth | ✗ | 不可用 | 現有資料限制 |

## §2 Credibility Score

| 維度 | 評分 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | ancestral S1-ALT reads 可同時包含 S2-REF 與 S2-ALT descendant states，具明確可觀察預測 |
| 觀察支撐 | 10 | FP-only 與早期 lineage pilot 有部分支持，但尚非全 sSNV / 最新 tag / 同一 endpoint |
| 機制清晰度 | 20 | DAG 可分 methyl group、genetic co-occurrence、HP/PS、technical、CN/LOH 與 truth class |
| 反例風險 | 0 | 既有研究顯示多群常見且可由多種非 clone 機制生成，必須強控制 |
| 所需資源 | 20 | 方法 pilot 已由 7,745 FP + 7,745 matched TP 完成；全量本身超過 6 小時但不阻擋可行性 |
| **TOTAL** | **70 / 100** | **GO full validation；claim 維持 PROBE** |

**Falsifier observable**：若 methyl group 不能輔助辨識共同 ancestral ALT 下的 latent genetic substructure，則在有足夠 partner coverage 的 stable multigroup sites 中，group x partner R/A 或 multi-marker genotype 的 FDR 顯著率不高於 permutation expectation、效應不跨 HP-family 存活、同一 HCC1395 位點不跨平台重現，且 complete trees 不支持一致 ancestral-to-descendant relation。

**Reality-test 三反例**：

1. methyl groups 與 partner sSNV 看似對齊，但完全由 HP family / phase set 分層解釋。
2. group x genotype 顯著只來自少數 reads、重複 readset、dense/capped component 或多重檢定。
3. 遺傳 association 存在，但 candidate trees 對 focal/partner 順序不唯一，只能稱 genetic-state association，不能稱 linear ancestry。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant / action |
|---|---|---|---|
| focal ALT 必須由 raw read 在焦點 base 判定，不以 HP tag 代替 | HIGH | KNOWN | schema hard gate |
| 最新 HP/PS 可由 sidecar 對 output read exact join | HIGH | contract known / 本輪尚待重驗 | MUST validate 100% |
| methyl grouping 不讀取 partner label | HIGH | KNOWN by design | leakage hard gate |
| selected-group partner 位點足以代表可驗 topology scope | HIGH | KNOWN for retained；FALSE for capped/singleton | 分支明示，不外推 |
| group x partner association 可等價 cellular clone | HIGH | UNKNOWN/FALSE as default | claim blocker |
| 單位點無 partner 時可由 calibration 推估個案 clone | HIGH | UNKNOWN | 只能報 empirical context，不確認個案 |

## §4 Quick Pilot / Entry Gate

1. 以既有最新-tag FP 7,745 sites 驗證 focal ALT、methyl matrix、BERNOULLI distance 與 stable grouping 可重現。
2. 驗證 raw BAM 產生的 methyl/distance 不依賴 HP tag；最新 HP/PS 僅作 post-hoc annotation/confound test。
3. 對已 retained 且有 selected group 的 stable cases，試算 group x partner R/A permutation association。

**Checkpoint**：只有 `site count exact + read-tag join 100% + no label leakage + null/FDR implementation tests PASS` 才啟動全量；否則 fail closed。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| 全 sSNV manifest、truth post-hoc label 與 selected-group partner contract | HIGH | 1-2 h | P0 |
| 469,849 per-site raw MM/ML + exact BERNOULLI artifacts | HIGH | 6-18 h | P0 |
| 最新 HP/PS sidecar exact read join + PS 欄位 | HIGH | 1-3 h | P0 |
| group x partner/genotype association、BH-FDR、within-HP sensitivity | HIGH | 3-8 h | P1 |
| topology relation parser與手算個案 | HIGH | 2-4 h | P1 |
| cellular clone orthogonal truth | HIGH | external | P2 / unavailable |

## §6 Evidence Conflict Scan

Repo root 無 `MEMORY.md`（2026-07-15 實查）；改查 `CURRENT_FOCUS`、experiments index、validated/in-progress methodology 與前一 cycle。

| Prior conclusion | Tier | Relation | Source |
|---|---|---|---|
| FP stable focal-ALT multigroup 583/4,967，但 orthogonal subclone confirmed=0 | L1/L3 | conflict with direct-clone claim | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/` |
| 甲基可有 CN-independent within-lineage bimodality，但 bulk 不能證 subclone | L1/L3 | support phenomenon / conflict claim | `InterSubMod/docs/experiments/in_progress/2026/07/20260706_methyl_vs_ssnv_lineage_01.md` |
| 甲基不能替代 sSNV co-occurrence backbone，只能 bounded-auxiliary | L1/L3 | dependent guardrail | `InterSubMod/docs/methodology/20260628_subclone_reconstruction_master_spec_01.md` |
| 單一 sSNV 無法識別 linear/branching，需多 marker relation | mathematical | conflict with single-site topology claim | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/` |

## §7 Red-Team Gate 與 Decision Path

**最強反方**：

1. 對齊可能是 allele-specific/cis methylation 或 phase-set linkage，不是 clone identity；即使 q-value 顯著也只證 read-level association。
2. 只在 genetically testable subset 評估會產生 selection bias；其 association rate 不能直接套到 positional singleton。
3. 同一 pipeline 產生 LongPhase tag 與 layered topology，不能當作完全獨立 orthogonal confirmation。

- **Validation verdict**: **GO**，因為 endpoint 可否證、既有 pilot 已通過可行性、資料與資源可用。
- **Biological claim verdict**: **PROBE**；只有 methyl-genetic association 可升為 `genetically anchored latent-substructure candidate`。
- **Decision lock**：沒有 partner genetic association / cross-platform replication / orthogonal cellular truth時，不得把單位點 methyl multigroup稱為 confirmed/high-probability subclone；沒有 candidate-tree一致順序時，不得稱 linear ancestry。

## Step -> Verify

1. 建立 469,849-site manifest與共現 contract -> 驗證：7/7 VCF SHA256、每 dataset count、總數與 latest SoT exact match。
2. 全量產生 focal read/methyl/distance artifacts -> 驗證：每個 VCF site恰有一組輸出或明確 fail-closed receipt，exit code 0。
3. 最新 sidecar join與 truth-blind grouping -> 驗證：所有分析 reads 100% exact join；partner/truth欄位不進入 clustering function。
4. 共現 association與 topology -> 驗證：permutation、BH-FDR、within-HP sensitivity、人工個案重算一致。
5. 報告與交接 -> 驗證：TSV/JSON/figures/HTML、browser QA、ledger/CURRENT_FOCUS/INDEX位置全部可讀。

## 2026-07-16 Runtime Rescue Addendum

> **Verdict: GO_WITH_EQUIVALENCE_GATES (75/100)**。這不是改變 M1 統計方法，而是把彼此獨立的
> 10 個 coarse seed runs 與 1 個 fine run 分派到不同 process；每個 run 內的 RNG seed、40 次
> column-null 呼叫順序、遞迴判定與輸出欄位維持不變。只有 serial/parallel exact-equivalence、
> 7-dataset merge reconciliation 與 input immutability 全 PASS 才可取代正在執行的 serial tail。

### 觸發觀察與 Cynefin gate

- 2026-07-16 只讀診斷確認 `HCC1954 chr17:39851484` 在 `4,610` 個 peeled focal-ALT reads下，
  每個 seed最多有 `82` 個可執行 null 的不平衡節點、最大深度 `74`，平方成本約為
  `65.5588` 個完整 root equivalents。
- 原 runner 以46個outer workers跑到tail後，只剩一個worker串行處理該ordered chunk；
  worker CPU time持續增加，沒有deadlock/OOM，但預估尚需數十小時。
- Domain為 **Complicated**：統計分解可證明獨立，但process nesting、RNG與merge provenance
  需要工程驗證，不可直接假定等價。

### Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | 11個 `phylo_label` runs只共享immutable distance/methylation arrays，沒有跨run狀態 |
| 觀察支撐 | 20 | serial straggler、82-node depth與CPU持續前進已由實際矩陣重算 |
| 機制清晰度 | 20 | seed-level parallelism不改每個seed內RNG stream；merge按canonical dataset/site order |
| 反例風險 | 5 | nested fork、浮點/RNG drift、prefix/replacement重複或遺漏皆可能 |
| 所需資源 | 10 | 實作與測試1-6h；replacement run仍屬長計算 |
| **TOTAL** | **75/100** | **GO_WITH_EQUIVALENCE_GATES** |

**Falsifier observable**：任一固定synthetic/real fixture的serial與parallel labels、split traces、
ARI或完整result payload不相等；replacement與merged output任一sample/site/assignment key缺失、重複；
或frozen input SHA改變。任一發生即不得採用rescue output。

### Assumption map與red team

| 假設 | 重要性 | 已知度 | 動作 |
|---|---|---|---|
| coarse/fine runs互相獨立 | HIGH | KNOWN from code | unit test exact payload |
| fork children只讀共享arrays | HIGH | KNOWN on Linux / runtime待驗 | nested-process smoke + RSS/exit audit |
| HCC1954 replacement可與前6 datasets合併 | HIGH | UNKNOWN until run | exact per-sample/key/assignment reconciliation |
| 目前partial gzip可當完整前6 datasets來源 | HIGH | KNOWN only for prefix | 只接受manifest固定前6 counts；忽略所有HCC1954 partial rows |

最強反方：一是平行執行可能改變RNG stream或代表run選擇；二是recovery merge可能把partial HCC1954
混入replacement；三是多層processes可能過度配置CPU/RAM。對應gate分別是exact-equivalence、
replace-sample排他key audit、outer+inner worker上限與terminal receipts。

### Step -> Verify

1. 新增seed-level parallel dispatch -> 驗證：serial/parallel完整payload exact equal。
2. 以新目錄重跑HCC1954 -> 驗證：`22,400/22,400`、summary/manifest PASS，不覆寫既有輸出。
3. 前6 datasets + replacement HCC1954合併 -> 驗證：469,849 unique keys、7個sample count exact、
   stable row與assignment key set exact、terminal tag audit PASS。
4. rescue前後input audit -> 驗證：77/77 frozen artifacts hash不變。

## 2026-07-17 Terminal Recovery + M2 Semantics Addendum

> **Verdict: GO_WITH_SOURCE_ATTESTATION (82/100)**。全量M1結果不重算；M2只修正文句與加入
> assignment-derived category proof。Tumor-REF只接受目前完整gzip中的canonical ordered prefix，缺失suffix
> 以已通過exact-equivalence的seed-parallel implementation重算，最後重新產生dedup欄位、summary與manifest。

### 新觀察與決策

- Ordered writer目前停在`102,000 / 102,842`，不是全流程沒有計算；後續chunk可能已在worker記憶中完成，但
  不能在前方straggler完成前寫出，因此不能當成可稽核結果。
- 實際長尾包含HCC1954 chr17高深度位點，例如`chr17:39520424 G>T`（ALT=987、REF=3,562、
  joint=4,549）與`chr17:39560316 C>G`（ALT=1,710、REF=2,006、joint=3,716）。Serial REF+joint各需
  11次RNULL40 recursive phylo runs，單worker已累積超過11小時CPU。
- 既有`phylo_parallel_exact_equivalence.pinned_390228ce.v1.json`已在synthetic及real fixture證明
  serial/parallel完整payload SHA相等。Recovery仍必須另做prefix exact-key、suffix exact-key、final
  102,842 unique key與source hash前後一致gate。
- M2 reviewer重算得到`eligible=919`、`evaluable_ineligible=948`；其中121個evaluable rows含已顯著aligned
  axis但該軸低於陰性判定所需planning power。v4文字將其一律寫成NOT_EVALUABLE，屬合約敘述錯誤；採用
  非對稱gate：陽性混雜可保守排除，未對齊才要求80% power。

### Assumption map與falsifiers

| 假設 | 重要性 | 驗證 |
|---|---|---|
| open gzip只取完整newline rows可形成canonical prefix | HIGH | 每row key逐一等於task order prefix；EOF/truncated tail不採用 |
| suffix seed-parallel與serial統計完全相同 | HIGH | pinned real+synthetic exact payload receipt；新run鎖定dependency SHA |
| 合併後dedup/summary可由完整rows重建 | HIGH | 重新按canonical order跑DedupAnnotator與SummaryBook；counts exact |
| categorical statistic缺值等於constant axis | HIGH | stable assignment直接計算HP exact/family/strand level count=1；若>=2則hard fail |
| 低power陽性混雜可用來排除M2 | HIGH | effect與permutation p均已實際觀察；僅作保守FAIL，不宣稱陰性power充分 |

**Falsifier observable**：prefix任一key不等於canonical task order、suffix與prefix重疊或缺漏、final key數非
102,842、source/input SHA在run中改變、assignment group counts不等於screen coarse cluster_sizes、或任一
categorical statistic缺值但observed levels>=2。任一發生即NO-GO，不得進入source-attested release。

### Step -> Verify

1. 鎖定M2 v5與observed category proof -> 驗證：logic-independent auditor不import production gate，
   102,842 assignment keys exact，919/948及121重算一致。
2. 建立tumor-REF recovery producer -> 驗證：prefix/suffix/final key conservation、102,842 rows、source/input
   SHA前後一致、serial-parallel dependency receipt exact。
3. 以fresh paths重跑cooccurrence與terminal downstream -> 驗證：G1/G2/R1/B1/C1全為實際完成或structured
   NOT_APPLICABLE，不沿用舊M2 output。
4. 最終多agent+Claude review -> 驗證：每項finding有修正或明示residual limitation，HTML desktop/mobile QA PASS。
