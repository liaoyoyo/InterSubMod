<!--
建立時間: 2026-07-27 10:42
目標: 凍結 multi-sSNV pattern × methylation 正式統計與分類門檻
處理範圍: 7 technical datasets / chr1-22 / exact PS × exact raw HP
status: frozen_before_full_run
關聯檔案:
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/pre-decision-audit.md
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/implementation-notes.md
-->

# Analysis Contract v1

## 1. Frozen universe

- July LongPhase-S recalibrated PASS autosomal biallelic sSNVs：469,849 sites。
- July exact-PS topology：98,955 final groups。
- Primary methyl source：receipt-validated 20260715 all-sSNV v2，與 site universe
  逐 `(dataset, chromosome, position, REF, ALT)` 對齊。
- Primary grain：`dataset × chromosome × exact PS × exact raw HP × topology region`。
- Candidate gate：active marker `k≥2`、complete R/A reads `N≥40`，且至少兩個
  complete states各 `n≥5`。pair full4及 `k≥3` 另列 strata，但不是唯一候選。

## 2. Read and CpG gates

1. `SHA256(raw read_name)` 必須唯一對上 strict sparse `qname_sha256`。
2. artifact HP/PS不具權威性；只採 sparse `hp_raw`、`phase_set`與 pattern。
   `pattern_counts` 與 candidate shard 必須逐 unit 重算
   `unit/PS/HP/active positions/N/complete/partial/state counts` 並完全相等。
   reads/methylation 與 topology JSONL 必須在消費時重驗 size + SHA-256。
3. 含 X pattern只列 partial subcube，不進 complete-state正式檢定。
4. 每個 marker取 ±5 kb，跨 marker以 `(qname_sha256, CpG position)` union；
   overlap值衝突即 fail closed。
5. 遮蔽所有 marker ±100 bp；distal sensitivity另只保留距所有 marker
   `>1,000 bp` 的 CpG。
6. common basis要求每個正式 state在該 CpG finite coverage ≥80%；每個 read
   在 basis finite coverage ≥80%；basis至少10個 distinct CpG。正式 state
   集合在 methylation 前凍結；exchangeability、read coverage與 state-CpG
   coverage以只會刪 row/CpG 的保守 fixed point共同收斂。任一輪
   `N<40`、任一凍結 state `n<5`或 basis<10即 `NOT_EVALUABLE`，不得刪掉
   該 state後繼續檢定。
7. Bernoulli distance完全對齊 C++：
   `delta=p(1-q)+(1-p)q`、`w=2|p-.5|×2|q-.5|`、共同 CpG至少3；
   invalid保持 NA，禁止改成1；fixed point後若仍有 invalid read pair，整個
   unit `NOT_EVALUABLE`，不以資料驅動 peel刪 read。
8. 每個 formal marker以 SHA-256 deterministic抽取最多16 reads、保留全部
   CpGs，對照 frozen BERNOULLI matrix；invalid mask必須零差異，因
   methylation.csv為小數4位序列化，finite cell絕對誤差容許上限預先固定為
   `5×10⁻⁴`。

## 3. Primary statistics

- 位置檢定：restricted PERMANOVA，固定 seed `20260727`，primary 999
  permutations，回報 pseudo-F、R²、raw p、realized permutations。
- Exchangeability：先只保留同 `read_group × strand` strata內至少有兩個
  pattern的 reads；若重新檢查後 `N<40` 或任一正式 state `n<5`，
  `NOT_EVALUABLE`，不得退回 unrestricted。
- Dispersion：同一 restricted permutation contract重算 centroid的 PERMDISP。
- Multiplicity預先分兩個 family：`pair_full4 ∪ k≥3` 是 confirmatory family；
  其餘 `k=2` any-two-state 是 comprehensive secondary pair census。各 family
  分別做 Benjamini–Yekutieli，Holm作 family-wise sensitivity；secondary
  finding只能列 exploratory ranking，不升格 robust claim。
- permutation採兩階段：confirmatory family先以999次 screen；secondary
  pair census以199次作 exploratory ranking。只有落在離散 p-value floor且
  已過 effect/confound gate的 confirmatory units，才以固定 seed增加
  permutation到49,999次上限。adaptive資格還必須先通過 `R²≥0.10`、
  最佳 pair contrast `≥0.10`、standardized effect `≥0.50`、
  PERMDISP `p≥0.05`、geometry SMD `<0.50`、全 state `n≥8`、等 N與
  rarefaction retention皆 `≥0.50`。若上限仍不足以判斷 BY門檻，標
  `resolution-limited`；不得把 `1/(B+1)` 當作更小的 p-value。
- adaptive 49,999次輸出一律標
  `PROVISIONAL_SUBSET_REFINEMENT_REQUIRES_FULL_FAMILY_MERGE`，其中的
  subset-only BY/Holm與 assessment不得解讀；只有重新併回完整 frozen
  family並重算 multiplicity後的 combined evidence可作正式結果。
- Pair relation：回報 within-state、between-state Bernoulli mean及
  `between − pooled within`；standardized effect固定為該 contrast除以
  concatenated between與within distance cells的樣本標準差。Hamming>1
  永不投影成 topology edge。

## 4. Robustness gates

- 明顯 effect門檻：primary `R²≥0.10`、最佳 state pair
  `between − pooled within≥0.10`。
- 顯著門檻：`BY q≤0.05`。
- dispersion gate：PERMDISP `p≥0.05`。
- read geometry gate：MAPQ、read length、start及end的跨 state最大
  pairwise standardized mean difference `<0.50`。
- support sensitivity：所有正式 state `n≥8`；`n≥10`另報，不作 primary硬門。
- equal-N sensitivity：固定 seed下每 state降採到最小 n，R²方向為正且
  至少保留 primary R²的50%。
- CpG rarefaction：10次固定 seed的80% CpG subsample，中位 R²至少保留
  primary R²的50%。

## 5. Reader-facing assessment

- `ROBUST_ASSOCIATION`：所有 primary與 robustness gates通過。
- `LOCAL_CIS_COMPATIBLE`：primary通過，但 distal basis不可評估或 distal
  R²低於 primary 50%；只能說與局部 cis相容。
- `TAG_DEPENDENT`：保留為未來正式 tag×pattern interaction/equivalence
  檢定的分類；本版不得以「另一 raw HP不顯著」推論 tag間有差異，因此不產生
  此 assessment。
- `CONFOUNDED`：有 raw association/effect，但 dispersion、geometry或
  已量測 geometry gate失敗。support/equal-N/rarefaction失敗不稱為
  confounded，而歸入未達 robust。
- `EVALUABLE_NO_ROBUST_ASSOCIATION`：可評估但未過顯著/effect門檻。
- `NOT_EVALUABLE`：資料、common-CpG、support或 exchangeability不足。

以上分類均為 association；禁止改寫 topology、ancestry、clone或因果方向。
