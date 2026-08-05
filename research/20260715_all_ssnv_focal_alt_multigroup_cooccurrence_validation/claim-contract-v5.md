<!--
建立時間: 2026-07-17T10:09:08+08:00
目標: 修正 M2 功效與陽性混雜證據的非對稱決策語意，並鎖定 observed categorical-level proof
處理範圍: 7 datasets / 6 biological samples / chr1-22 / 469,849 LongPhase-S recalibrated FILTER=PASS dataset-sites
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v4.md
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/m2_screen_gate.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/audit_independent_m2_gate_recount.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/build_all_ssnv_final_report_dataset.py
-->

# all-sSNV focal-ALT source-attested release claim contract v5

`terminal-claim-contract-v5`

## 版本關係與修正理由

- `claim-contract-v2.md` 是 immutable source-locked screen contract v2；SHA-256 為
  `904abd6c9fbcf152770be72a3cd2b12f38d0058b4d668b2c6ced9005b86afc22`。
- `claim-contract-v3.md` 是 terminal integration contract；SHA-256 為
  `da94a50d0717174ff007b75f2edad2de79bf3aebf6b15df179eb736e8d8f526e`。
- `claim-contract-v4.md` 鎖定 2-10 群與背景重現 predicate；SHA-256 為
  `3bcc2d042708e34ec7c75db7fe9c65baef9b80494f668f911f9551ef0e588515`，保持 immutable。
- v4 把「低於 80% planning power」一律寫成 NOT_EVALUABLE，與已執行 gate 的非對稱決策不完全一致。
  本 v5 明確修正：已觀察到 effect 達門檻且 permutation `p<0.05` 是陽性混雜證據，可在 planning
  power 低於 80% 時保守排除；只有要把未對齊結果解讀為「未觀察到指定大小混雜」時，才要求足夠
  planning power。此修正不改 methyl assignment、effect、p-value、M1 或既有 screen 數值。

## 統計單位

- Primary key：`(dataset, CHROM, POS, REF, ALT)`。
- Biological-site key：`(biological_id, CHROM, POS, REF, ALT)`；HCC1395_DORADO 不增加 biological n。
- Pair key：primary key加 `(partner_CHROM, partner_POS, partner_REF, partner_ALT)`。
- 469,849 是 dataset-sites，不是獨立病人數、truth-confirmed somatic sites或 globally unique loci。
- 所有比率必須保存 numerator、denominator、denominator definition、NOT_EVALUABLE 與 NOT_RUN。

## Evidence ladder

| ID | 正式名稱 | 必要條件 | 可支持 | 不可支持 |
|---|---|---|---|---|
| M1 | focal-ALT operational stable methyl-multigroup screen flag | focal ALT reads通過預定資料量與聚類screen；coarse K>=2；modal K與membership穩定 | operational screen標記的ALT-carrier read-level methyl heterogeneity | biological prevalence、genetic clone |
| M2 | residual robust epigenetic partition | M1；2-10群；`m2-measured-axis-v4_asymmetric-positive-confound-and-observed-categorical-levels`；八軸無已檢出對齊且所有未對齊軸可作陰性判定 | 已測HP/strand/geometry/read-quality後仍存在的epigenetic partition | 排除所有confound、somatic clone |
| G1 | LongPhase-S callset-anchored local co-segregation | M2；完整eligible pair family fixed-margins exact test通過global BY；V>=0.30；delta AF>=0.50；conditional與callability gates通過 | methyl group與同次LongPhase-S PASS call在相同molecules局部共分離 | partner為正交somatic truth、cellular clone |
| G2 | multi-marker molecular-haplotype base candidate | G1；同一complete-read set至少2個effect-supported markers且相距>=20 bp；joint test可執行；`joint_signature_q_global_by<=0.05`且`joint_signature_global_by_discovery=true` | 共同ancestral ALT下的multi-marker molecular substructure candidate | marker獨立、cellular clone、linear ancestry |
| R1 | strict methyl-partition robustness | G2；999 permutations x 10 seeds；預定null/resampling gates通過 | partition對預定擾動穩健 | post-selection genetic FDR、clone |
| B1 | background-controlled molecular-haplotype candidate | G2+R1；tumor-REF與matched-normal可評估且未通過lenient背景重現predicate；單一預選G1 pair通過four-state gate | 經目前bulk guardrails仍相容的molecular-haplotype mixture | cellular subclone、clone數、lineage |
| C1 | CN/CCF-conditioned candidate | B1；authority-reviewed exact-locus CN、purity與fit-local CCF可用且相容 | 指定CN/CCF模型下候選 | model-independent clone truth |
| L1 | cellular subclone supported | C1；single-cell/colony/spatial/multi-region正交資料對應同一cell population | cellular subclone support | linear ancestry，除非另有order evidence |
| L2 | lineage/order supported | L1；>=3 mutation states、CCF/order與獨立資料一致，且排除branching/recurrence/CN替代解 | 指定模型下lineage/order support | 唯一真樹，除非可識別性成立 |

## M2 asymmetric measured-axis gate

- 正式 contract：`m2-measured-axis-v4_asymmetric-positive-confound-and-observed-categorical-levels`。
- 八軸：HP exact、HP family、strand、read start、read end、read length、MAPQ、called CpGs。
- Categorical threshold為 Cramer's V>=0.30；continuous threshold為 eta-squared>=0.14；499次
  permutations使用add-one grid且要求 `p<0.05`。
- **陽性判定**：effect達門檻且permutation `p<0.05` 即為 observed aligned confound evidence；即使
  planning power低於80%，仍可保守列為M2 FAIL。低power不會推翻已觀察到的陽性證據。
- **陰性判定**：effect低於門檻時，只有 planning power>=80%、alpha=0.05且每群至少5 reads，才能把該軸
  列為 determinate not-aligned；否則為 NOT_EVALUABLE。effect達門檻但permutation p未過亦為
  NOT_EVALUABLE，不可當成陰性。
- Categorical effect/p同時缺值時，必須由 stable assignment 的 `core_reads.latest_hp`、HP-family mapping、
  `core_reads.strand` 實算 observed level count=1，才可列為 determinate constant axis；觀察>=2類卻缺統計量
  必須hard fail。production gate與logic-independent auditor均需驗證 coarse `labels` count等於screen
  `cluster_sizes`。
- Planning model只支援2-10個methyl groups。超過10群保留M1 PASS，但M2標為
  `NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM`，且G1/G2/B1為NOT_RUN。
- Effect與permutation p來自frozen source-locked screen；terminal只重算分類、n/grid/power及assignment
  category proof，不宣稱重新由raw reads估計axis statistics。

## G1/G2 multiple testing

- G1 family包含所有M2-eligible、exact-testable pairs；正式discovery使用global BY，BH只作對照。
- G2 joint family為 `ELIGIBLE_M2_JOINT_SIGNATURE_GLOBAL_FDR_FAMILY`；禁止只在顯著候選做FDR。
- G2正式欄位必含 `joint_signature_global_fdr_family_status`、`joint_signature_q_global_bh`、
  `joint_signature_q_global_by`、`joint_signature_global_bh_discovery`、
  `joint_signature_global_by_discovery`。Raw p或單marker p不得單獨形成G2。

## B1 background-control replication gate

- Tumor-REF與matched-normal REF共用lenient predicate：`coarse_ng>=2` 且 `unstable=false`；不要求
  membership ARI>=0.8。
- 同一background payload下，此predicate是加入ARI條件之predicate的superset；但不得宣稱observed focal
  ALT、tumor-REF與normal-REF實際flag集合互為superset。
- B1 pair在看four-state結果前固定為
  `among_G1_formal_pairs_max_endpoint_a_n_informative_then_abs_distance_then_partner_identity_without_four_state_result`。
- RR/AR/RA/AA使用 `bonferroni_three_relation_models`；每模型98.3333%單側upper bound維持familywise 95%
  confidence；2% error ceiling下的 minimum zero-violation depth = 203。
- Four-state只支持fixed error model下的mutation-order compatibility，不確認祖先、直系或唯一樹。

## Denominator與敘述鎖

1. M1主結果是469,849 dataset-sites中的operational yield，不是biological prevalence。
2. M2 PASS/FAIL分母包含：M1 PASS、2-10群，且每個軸不是已觀察到陽性對齊，就是具有足夠陰性資訊，或由
   observed category count=1證明constant。M2 PASS另要求沒有任何aligned confound；低功效且未對齊、缺統計
   proof、effect高但p未過及K>10均為NOT_EVALUABLE。
3. 陽性混雜證據即使低於80% planning power仍可進入M2 FAIL；這是保守排除，不得被描述為M2可確認陰性。
4. G1分母是M2 PASS且至少一個exact-testable partner者；正式pair discovery採完整family global BY。
5. G2分母是具有預定multi-marker opportunity、joint complete reads與可執行joint test的M2 sites；正式結果
   另需joint global BY discovery。
6. R1/B1/C1只在上游base candidates計算；零候選必須是structured NOT_APPLICABLE。
7. B1只有tumor-REF與matched-normal兩個lenient背景gate均可評估且不重現時才可PASS。
8. HCC1395與HCC1395_DORADO是biological n=1的technical comparison，不得作跨樣本重現推論。
9. L1前只能稱molecular-haplotype或read-level epigenetic evidence；禁止稱confirmed subclone、clone數或
   linear evolution。單一focal ALT的methyl groups不能識別linear與branching ancestry。
