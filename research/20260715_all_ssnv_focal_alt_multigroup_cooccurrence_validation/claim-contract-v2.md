<!--
建立時間: 2026-07-15T09:30:00+08:00
目標: 固定全 sSNV focal-ALT methyl multigroup 下游整合與報告的 v2 claim contract
處理範圍: 469,849 LongPhase-S PASS sSNV；stable subset cooccurrence；candidate controls；CN/CCF/cellular ceiling
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/reviews/20260715_preexecution_multiagent_red_team.md
-->

# all-sSNV focal-ALT claim contract v2

## 統計單位

- Primary key：`(dataset, CHROM, POS, REF, ALT)`。
- Biological-site key：`(biological_id, CHROM, POS, REF, ALT)`；HCC1395_DORADO不增加 biological n。
- Pair key：primary key加上`(partner_CHROM, partner_POS, partner_REF, partner_ALT)`。
- 所有 ratio必須同時保存 numerator、denominator、denominator definition與 not-evaluable count。

## Evidence ladder

| ID | 正式名稱 | 必要條件 | 可支持 | 不可支持 |
|---|---|---|---|---|
| M1 | focal-ALT stable methyl multigroup | evaluable；coarse K>=2；modal K穩定；modal read membership ARI min>=0.8 | ALT-carrier read-level methyl heterogeneity | genetic clone |
| M2 | residual robust epigenetic partition | M1；HP/technical measured axes不對齊；預定 read/group/phase門檻 | 已測量 confound後仍存在的 epigenetic partition | somatic clone |
| G1 | LongPhase-S callset-anchored local co-segregation | M2；partner Kx2 fixed-margins exact p通過global BY；V>=0.30；delta AF>=0.50；HP-family×PS×strand conditional sensitivity pass；differential-callability gate pass | methyl group與同次 LongPhase-S PASS call在相同 molecules局部共分離，且未偵測到預定義的強 callability 對齊 | partner已被正交確認為somatic；cellular clone；已排除所有 read-geometry confound |
| G2 | multi-marker molecular-haplotype base candidate | G1 pair中至少2個相隔>=20 bp的spatially separated markers；至少2個在同一 coverage-preselected complete-read set 各自 R/A-testable 且保留 V>=0.30、delta AF>=0.50；同批complete reads的joint signature sensitivity pass | 同一 focal ALT內有與多個PASS calls一致的latent molecular substructure | marker統計獨立；subclone confirmed/compatible |
| R1 | strict methyl-partition robustness | G2；999 permutations×10 seeds；column與row-circular null；exact partition；fixed-K LOO；artifact hash lock | partition對預定null/resampling穩健 | FDR-confirmed genetic association；clone |
| B1 | background-controlled molecular-haplotype candidate | G2+R1；tumor-REF控制可評估且不重現；normal focal called>=10、REF>=5、ALT=0；normal REF-only methyl control可評估且不重現；至少1個同一BY pair通過four-state error-model witness | 經目前bulk guardrails仍相容的molecular-haplotype mixture | cellular subclone、clone數、lineage |
| C1 | CN/CCF-conditioned candidate | B1；authority-reviewed allele-specific CN、purity、mutation multiplicity/CCF可用且相容 | 在指定CN/CCF模型下與subclonal mixture相容 | model-independent clone truth |
| L1 | cellular subclone supported | C1；single-cell/colony/spatial/multi-region等正交資料對應同一cell population | cellular subclone supported | linear ancestry，除非另有order evidence |
| L2 | lineage/order supported | L1；>=3 mutation states、CCF/order與獨立資料一致，排除branching/recurrence/CN替代解 | 指定模型下的lineage/order support | 唯一真樹，除非可識別性成立 |

## Required machine fields

### Screen

- `stable_null_multigroup`
- `modal_assignment_ari_min`
- `hp_axis_confound`
- `technical_axis_confound`
- `residual_unexplained_multigroup`
- `phase_anchored_robust_epigenetic_candidate`
- latest HP/PS join status/counts與 primary artifact SHA。

### Pair inference

- descriptive：analytic test/p。
- primary：fixed-margins exact p、global BH q、global BY q、exact state-space status。
- sensitivity：HP-family×PS×strand conditional permutation p、permutations、permutable/status。
- effect：Cramer's V、delta ALT fraction、informative reads。
- four-state：RR/AR/RA/AA/O/X、called depth、violation rate、exact-binomial p/CI、error ceiling、所有相容 relation models 與 status；禁止用優先順序隱藏多模型相容。
- technical replication：雙平台同一 exact pair的BY/sensitivity/effect-direction/four-state狀態；`biological_n=1`。

### Site inference

- `n_pair_by_confirmed`
- `pair_by_confirmed_positions`
- `n_spatially_separated_pair_by_20bp`
- `top_marker_positions`
- `n_top_marker_pair_by_confirmed`
- joint complete-read signature table/effect/sensitivity p/status。
- `multi_marker_molecular_haplotype_base_candidate`
- 任何 legacy `independent_marker`欄位只能作deprecated alias，不得出現在正式敘述。

### Controls and integration

- strict formal/exploratory status與`pass_semantics`。
- tumor-REF evaluability、methyl partition與 not-evaluable reason。
- normal ALT/REF/UNKNOWN/called counts、focal-callability gate、REF-only methyl evaluability、not-evaluable reason。
- same-pair witness key list；不同pair不得拼接gate。
- CN/purity/CCF source、exact-locus coverage、model status與blocked reason。
- technical replication與biological replication分開。

## Denominator rules

1. 469,849是M1 screen母體，不是所有下游endpoint的共同分母。
2. cooccurrence母體是M1 stable sites；沒有PASS partner也要保留0-partner狀態。
3. G1 denominator是至少1個exact-testable partner的M2 sites，另報全部M1 denominator。
4. G2 denominator是具有至少2個spatially separated testable top markers且joint complete reads足夠的sites。
5. R1、B1、C1只在其上游base candidate集合內計算；不可評估不能合併為fail。
6. technical replication denominator只限HCC1395 pair具有exact shared focal+partner opportunity的rows；非HCC不進分母。
7. CN/CCF denominator只限authority-reviewed source覆蓋的候選；COLO829/HCC1937不得假設CN=2。

## Formal language lock

- 使用`spatially separated markers`，禁止稱`independent markers`。
- 使用`LongPhase-S PASS callset-anchored`，除非truth/normal/CN另證，禁止簡寫成`confirmed somatic anchor`。
- 使用`molecular-haplotype candidate`，在L1前禁止稱`subclone-compatible`或`confirmed subclone`。
- four-state只稱`mutation-order compatibility under a fixed error model`，禁止稱祖先已證明；零 violation 的 95% 上界需由相應 opportunity denominator（預定 2% ceiling 時至少 149）支持。
- strict只稱`methyl-partition robustness audit`，禁止稱post-selection FDR confirmation。
