<!--
建立時間: 2026-07-17T00:00:00+08:00
目標: 固定 source-attested release 的 M2 群數分母與 B1 背景重現 predicate
處理範圍: 7 datasets / 6 biological samples / chr1-22 / 469,849 LongPhase-S recalibrated FILTER=PASS dataset-sites
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v2.md
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v3.md
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/m2_screen_gate.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/build_all_ssnv_final_report_dataset.py
-->

# all-sSNV focal-ALT source-attested release claim contract v4

`terminal-claim-contract-v4`

## 與 source-locked v2 及 terminal v3 的關係

- `claim-contract-v2.md` 是 screen 啟動時的 source-locked screen contract v2；SHA-256 為
  `904abd6c9fbcf152770be72a3cd2b12f38d0058b4d668b2c6ced9005b86afc22`，不得回寫。
- `claim-contract-v3.md` 是 terminal downstream 的 integration contract；SHA-256 為
  `da94a50d0717174ff007b75f2edad2de79bf3aebf6b15df179eb736e8d8f526e`，不得回寫。
- 本 v4 是完成全量 screen 後、正式 source-attested release 前的語意封閉契約。它不更改 methyl group
  assignment、M1 screen output、M2 gate 數值或任何已執行統計；只明確鎖定 v3 已設定的 2-10 群
  planning model 分母，以及實際 B1 background-control replication predicate。
- intermediate terminal report 可保留 v3 identity；只有通過 tumor-REF bounded retrospective source
  identity gate 的 Task Type B release 才可使用本 v4。正式 release 必須同時保存 v2、v3 與 v4 的用途與 identity。

## 統計單位

- Primary key：`(dataset, CHROM, POS, REF, ALT)`。
- Biological-site key：`(biological_id, CHROM, POS, REF, ALT)`；HCC1395_DORADO 不增加 biological n。
- Pair key：primary key加上 `(partner_CHROM, partner_POS, partner_REF, partner_ALT)`。
- `469,849` 是 dataset-sites，不是 469,849 個獨立病人、truth-confirmed somatic sites或 unique loci。
- 每個 ratio 必須保存 numerator、denominator、denominator definition、not-evaluable 與 not-run。

## Evidence ladder

| ID | 正式名稱 | 必要條件 | 可支持 | 不可支持 |
|---|---|---|---|---|
| M1 | focal-ALT operational stable methyl-multigroup screen flag | focal ALT reads通過預定資料量與聚類 screen；coarse K>=2；modal K與membership穩定 | 被此 operational screen標記的 ALT-carrier read-level methyl heterogeneity | 未標記site為同質；生物 prevalence；genetic clone |
| M2 | residual robust epigenetic partition | M1；2-10個methyl groups；`m2-measured-axis-v3_effect-permutation-and-power-evaluability`；所有八個 measured axes可判定且無對齊 | 已測量 HP/strand/geometry/read-quality axes後仍存在的epigenetic partition | 已排除所有confound；somatic clone |
| G1 | LongPhase-S callset-anchored local co-segregation | M2；完整 eligible pair family 的 fixed-margins exact test通過 global BY；V>=0.30；delta AF>=0.50；conditional與callability gates通過 | methyl group與同次 LongPhase-S PASS call在相同 molecules局部共分離 | partner為正交確認somatic；cellular clone |
| G2 | multi-marker molecular-haplotype base candidate | G1；同一complete-read set至少2個effect-supported markers且相距>=20 bp；joint test可執行；`joint_signature_q_global_by<=0.05`且`joint_signature_global_by_discovery=true` | 共同 ancestral ALT之下存在可重現的 multi-marker molecular substructure candidate | markers獨立；subclone-compatible；cellular clone；linear ancestry |
| R1 | strict methyl-partition robustness | G2；999 permutations x 10 seeds；預定null/resampling gates通過 | partition對預定擾動穩健 | post-selection genetic FDR confirmation；clone |
| B1 | background-controlled molecular-haplotype candidate | G2+R1；tumor-REF與matched-normal可評估且未通過lenient背景重現predicate；只用預先指定的單一 G1 pair通過four-state gate | 經目前bulk guardrails仍相容的molecular-haplotype mixture | cellular subclone、clone數、lineage |
| C1 | CN/CCF-conditioned candidate | B1；authority-reviewed exact-locus CN、purity與fit-local CCF可用且相容 | 指定CN/CCF模型下的候選 | model-independent clone truth |
| L1 | cellular subclone supported | C1；single-cell/colony/spatial/multi-region等正交資料對應同一cell population | cellular subclone support | linear ancestry，除非另有order evidence |
| L2 | lineage/order supported | L1；>=3 mutation states、CCF/order與獨立資料一致，且排除branching/recurrence/CN替代解 | 指定模型下的lineage/order support | 唯一真樹，除非可識別性成立 |

## M2 measured-axis 與群數 gate

- 正式 contract：`m2-measured-axis-v3_effect-permutation-and-power-evaluability`。
- 八軸：HP exact、HP family、strand、read start、read end、read length、MAPQ、called CpGs。
- Categorical effect threshold為 Cramer's V>=0.30；continuous threshold為 eta-squared>=0.14；permutation
  `p<0.05`，499 permutations使用 add-one p-grid。
- Planning-power目標為80%、alpha=0.05、每群至少5 reads、只支援2-10個methyl groups。低於80% power、
  缺統計量，或 effect達門檻但 permutation p未過者均為 NOT_EVALUABLE。
- M1 FLAGGED但觀察到超過10群的位點保留M1 PASS，M2必須標為
  `NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM`；不得加入M2 PASS/FAIL denominator，
  且G1/G2/B1為NOT_RUN。禁止為了納入觀察結果而事後擴張10群上限。
- Effect與permutation p來自 frozen source-locked screen producer；terminal gate會重算分類、n/grid/power，
  但不宣稱從原始 reads重新估計這些 axis statistics。

## G1/G2 multiple testing

- G1 pair family包含所有 M2-eligible、exact-testable pairs；正式 discovery 使用 global BY，BH只作較不保守對照。
- G2 joint family為 `ELIGIBLE_M2_JOINT_SIGNATURE_GLOBAL_FDR_FAMILY`；禁止只在已顯著或已選候選上做FDR。
- G2正式欄位必含 `joint_signature_global_fdr_family_status`、`joint_signature_q_global_bh`、
  `joint_signature_q_global_by`、`joint_signature_global_bh_discovery`、
  `joint_signature_global_by_discovery`。
- Raw joint p或單一 marker p不得單獨形成 G2。

## B1 background-control replication gate

- Tumor-REF與matched-normal REF使用同一個lenient predicate：`coarse_ng>=2` 且 `unstable=false`；
  `unstable=false`包含modal fraction gate，但不要求membership `ARI>=0.8`。
- 在同一個background payload上，lenient predicate是加入ARI條件之predicate的superset。B1要求背景不重現，
  因此lenient gate不可能增加B1 PASS；只可能把K穩定但membership不穩定的背景結構保守列為重現而減少候選。
- Superset敘述不可跨payload套用：不得宣稱observed focal ALT、tumor-REF與normal-REF實際flag集合互為superset。
- B1 pair在查看four-state結果前，固定以 G1 formal pairs 中 `endpoint_a_n_informative`最大者選取；再以
  absolute distance與partner identity固定tie-break。正式字串為
  `among_G1_formal_pairs_max_endpoint_a_n_informative_then_abs_distance_then_partner_identity_without_four_state_result`。
- RR/AR/RA/AA 同時評估三個relation models，使用 `bonferroni_three_relation_models`。每模型使用98.3333%
  單側upper bound維持familywise 95% confidence；2% error ceiling下的 minimum zero-violation depth = 203。
- Four-state只表示fixed error model下的mutation-order compatibility，不能確認祖先、直系或唯一樹。

## Denominator與敘述鎖

1. M1主結果是469,849 dataset-sites中的operational flag yield；不得稱biological prevalence。
2. M2分母是M1 PASS、觀察2-10個methyl groups、八軸皆可判定且達power要求者；群數超過10、低power、
   缺統計量或高effect但p未過者為NOT_EVALUABLE，不得併入FAIL。
3. G1分母是M2 eligible且至少一個exact-testable partner者；正式pair discovery採完整family global BY。
4. G2分母是具有預定multi-marker opportunity、joint complete reads與可執行joint test的M2 sites；正式結果另需
   joint global BY discovery。
5. R1/B1/C1只在上游base candidates內計算；零候選必須是structured NOT_APPLICABLE。
6. B1只有在tumor-REF與matched-normal兩個lenient背景gate均可評估且不重現時才可PASS。
7. HCC1395與HCC1395_DORADO只構成biological n=1的technical comparison，不得作跨樣本重現推論。
8. 在L1前只能稱molecular-haplotype或read-level epigenetic evidence；禁止稱subclone-compatible、confirmed
   subclone、clone數或linear evolution。
