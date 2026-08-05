<!--
建立時間: 2026-07-15T12:00:00+08:00
目標: 保存全 sSNV focal-ALT 甲基多群與 sSNV 共現關聯驗證的設計決定、偏離、折衷與未決事項
處理範圍: 7 datasets / 6 biological samples；469,849 latest LongPhase-S PASS chr1-22 biallelic sSNV
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/pre-decision-audit.md；InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_input_manifest.json
-->

# Implementation Notes: all-sSNV focal-ALT methyl multigroup x co-occurrence

## 設計決定

- 主單位為 `(dataset, CHROM, POS, REF, ALT)`；所有 469,849 位點都必須得到結果或明確不可評估原因。
- focal ALT 只由 raw BAM 在焦點 base、baseQ>=20 的 `alt_support=ALT` 判定；LongPhase-S HP/PS 不參與 focal allele 定義。
- methyl groups 使用 ±5 kb MM/ML、BERNOULLI `C_min=3`、complete-link peel、phylo-v4.1 column-null、10 seeds；partner sSNV label 在群建立後才加入。
- 最新 LongPhase-S HP/PS 從 `20260711_longphase_s_raw_all_production_sidecars_v2` sidecar逐 alignment join；原始 BAM不覆寫，BAM內舊 HP不能作最新結果。
- 共現分兩層：site-level `ssnv_branch/component/selected_group/family/topology` 註記；read-level 則對 ±5 kb 內**所有 frozen LongPhase-S PASS partner**建立 R/A/X genotype。selected-group 不限制 partner universe，避免只觀察 layered retained subset。
- 主要 genetic endpoint為 focal-ALT core reads 的 `methyl group x partner R/A` 與 `methyl group x multi-marker genotype` association，做 global BH-FDR、effect-size gate與 within-HP-family/strand conditional sensitivity。
- ancestry compatibility 必須另納入 focal-REF reads，計算 focal/partner `00,10,01,11` 四狀態；只看 focal ALT 的 `10/11` 無法驗證 focal 是否為共同祖先。
- topology只對 retained、candidate-complete unit 判 ancestral/reverse/coincident/ambiguous；capped/singleton不得補推 topology。
- truth TP/FP/UNASSESSED只在群建立後作 post-hoc分層，不是輸入或選群條件。

## 偏離

- 不建立或覆寫全基因組 latest-tag BAM。甲基與 allele直接讀 raw MM/ML BAM，HP/PS由 frozen sidecar在分析階段 exact join，避免約 2 TB BAM複製。
- 早期 `methyl_vs_ssnv_lineage.py`使用整 read平均甲基與舊 tagged BAM；本輪僅作歷史敏感度，不作主要 endpoint。
- layered `max_snv_excluded` 與 `positional_singleton` 沒有完整 selected-group topology；仍分析 focal methyl，但共現/topology欄位必須明示不可用原因。

## 折衷

- 全量 per-site artifacts使用既有 InterSubMod binary與相同參數，保留 exact BERNOULLI matrix語義；代價是較長執行時間與約 0.1 TB級輸出。
- HCC1395與DORADO先分 dataset報告，再以同 biological site作跨平台 replication；不可把兩者當獨立 biological replicates計顯著性。
- [歷史數值，2026-07-18 已 supersede，不得引用] HCC1395/HCC1395_DORADO最新LongPhase-S PASS exact shared records=`76,787`，分別占兩callset的`96.3608%`與`96.2979%`，site Jaccard=`0.929186`。
- [更正 2026-07-18] 以 `(CHROM,POS,REF,ALT)` variant key 重算：intersection=`76,721`、union=`82,705`、Jaccard=`92.764645%`；intersection分別占HCC1395與HCC1395_DORADO callset的`96.277937%`與`96.215152%`。另有1個同座標但ALT不同的位點：`chr1:70439946 C>A`對`C>T`。技術重現只在exact-key shared universe評估，不把platform-specific sites當replication failure。
- genetic-testable subset的 association rate只描述「可觀察 partner 的位點」；不得直接當 positional singleton個案的 posterior probability。
- 全量 phylo-v4.1 `RNULL=40` 是與既有 FP/TP 可比較的 discovery screen，不是 genome-wide FDR；genetic association 同時報全域 BH 與較保守 BY。strict methyl null只作候選後穩健性稽核，因候選已使用甲基資訊選出，其 p/BH/BY 不具正式 post-selection FDR 校準，不得稱「確認」。

## 未決

- 全量 stable multigroup數量未知；完成主分析後再以預先固定規則啟動共現 pileup，不能先看結果改 threshold。
- candidate-tree relation parser需以至少三個人工 state-tree個案與現有 V4/V5 solver輸出交叉驗算。
- 是否能對單位點建立可校準推估模型，取決於跨樣本/cross-platform holdout discrimination；若校準不佳，只報 population-level context，不輸出個案 clone probability。

## 2026-07-15 執行證據

- 全量 manifest：`results/all_ssnv_input_manifest.json`；`pass=true`，`469,849 = 335,296 TP + 7,745 FP + 126,808 UNASSESSED`，7/7 樣本 VCF key 與 latest PASS ledger exact match。
- latest-tag join synthetic tests：`7 passed in 0.13s`。真實回歸：`7,745/7,745` 既有 FP sites、`434,759/434,759` output read rows 與同 run sidecar HP 一致；eligible projected full identity multiplicity全部為 1；結果在 `results/latest_tag_projected_join_audit.json`。
- 全PASS幾何partner上限（不讀甲基/群標籤）：`238,664/469,849=50.7959%` focal sites在±5 kb內至少有一個其他PASS sSNV，共`864,255` unordered site pairs。逐dataset focal opportunity/pairs為HCC1395 `49,341/79,687; 287,727`、DORADO `48,957/79,739; 285,217`、COLO829 `8,351/37,788; 6,975`、H1437 `38,233/77,080; 142,185`、H2009 `86,763/154,465; 122,004`、HCC1937 `3,148/18,690; 10,572`、HCC1954 `3,871/22,400; 9,575`。這是高度相依的幾何上限，不是joint-spanning/testable或獨立檢定分母。
- 發現並處理的邊界：sidecar 可有同投影但異 tag的 supplementary records；join 明確複製 InterSubMod FLAG exclusion（secondary/supplementary/duplicate/unmapped），eligible primary 若仍衝突則 hard fail。
- v1 執行發現 verification schema 邊界：current class `Noise_Uncorrelated` 被錯誤寫入只接受 legacy 四分類的欄位。v1 已中止並保留在 `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/intersubmod_all_ssnv_v1/`，receipt `results/v1_aborted_run_incident.json` 明示 `pass=false`，不得納入分析。
- 修正於 `RegionProcessor::resolve_verification_legacy_input`：current class若不是合法 legacy值，保留既有合法 legacy class；兩者皆不合法才 fail。isolated build binary SHA256=`8d082edfea79b1e75600e377d6f00f90dabdf071b739054d65fdfa4806a10eb7`；相關 GTest 5/5、全 CTest 252/252通過。
- 修正 smoke輸入包含 HCC1395 `chr1:98311` 與曾失敗的 `chr1:13201492`；2/2 success、0 region failure，兩版本 reads/methylation/BERNOULLI 共6組 artifact byte-identical，receipt在 `results/verification_fix_smoke.json`。
- 正式 v2 command：`python3 research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_all_ssnv_intersubmod.py --binary build_all_ssnv_verification_fix/bin/inter_sub_mod --output-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/intersubmod_all_ssnv_v2_verification_fix --workers 7 --threads-per-sample 6 --minimum-free-gib 300 --summary research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_intersubmod_batch.v2_verification_fix.json`。
- v2 使用42 threads（主機48 logical CPUs），BAM唯讀且沒有 materialize/overwrite。完成後必以 `audit_all_ssnv_output_completeness.py` 對 `(sample,chrom,pos)` 集合做 exact reconciliation，不能只依 batch receipt。
- frozen input pre-run audit：`results/frozen_input_immutability.pre_v2.json`，7/7 samples與77/77 input artifacts PASS。大型BAM/sidecar依frozen size+nanosecond mtime，其他manifest artifacts另做SHA-256；v2結束後須另產生post-run receipt，不能覆寫pre-run檔。
- 全量 screen analyzer：`scripts/analyze_all_ssnv_focal_alt_multigroup.py`。它在 ALT選取前回接 latest HP/PS，truth/ledger/cooccurrence只在 clustering結束後加入，並以 bounded chunk futures輸出 gzip TSV/JSONL；primary screen的 strict status固定為 `NOT_RUN`，避免把 RNULL=40 screen誤稱確認。
- analyzer real-artifact smoke：HCC1395 `chr1:98311` latest-tag join 13/13、projection multimatch=0、舊HP被取代8、ALT=7、`coarse_ng=1`；曾觸發v1 schema failure的`chr1:13201492` join 9/9、multimatch=0、舊HP被取代9、ALT=9、`coarse_ng=2`且screen stable。兩者皆`strict_confirm_status=NOT_RUN`，證明修正後artifacts可由新analyzer讀取，但不把screen當strict結果。
- 2026-07-15 07:43 +08:00 topic Python tests為 `26 passed`，涵蓋 manifest/run、latest-tag join、exact reconciliation、all-sSNV screen、cooccurrence primitives與normal BAM audit。

## Matched-normal 控制輸入

- 稽核 command：`/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/audit_matched_normal_methyl_tags.py --target-reads 1000 --max-scanned 100000 --output research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/matched_normal_methyl_tag_audit.v2.json`。v2 receipt含完整command與script SHA-256，且已驗證重跑同一路徑會fail closed拒絕覆寫。
- 結果：7/7 BAM與index可讀，7/7各抽樣1,000 primary reads皆為 MM+ML pair fraction=`1.0`。HCC1395與COLO829同時含`C+m?`/`C+h?`；其餘為`C+m?`。這更新了舊文件中「COLO829無normal甲基」的過期狀態。
- 所有抽樣normal reads均無HP tag；因此normal只能作局部甲基背景與REF-read控制，不能直接作latest LongPhase-S haplotype或clone標籤。
- paired pilot使用H2009 `chr3:193395128 C>T`的LongPhase-S PASS單位點VCF、tumor BAM與H2009BL BAM；command為`build_all_ssnv_verification_fix/bin/inter_sub_mod -t <H2009 tumor BAM> -n <H2009BL BAM> -r <GRCh38> -v <single-site VCF> -o <pilot output> -w 5000 -j 2 --distance-metric BERNOULLI --min-common-coverage 3 --nan-distance-strategy SKIP --methyl-low 0.2 --methyl-high 0.8`，exit=0、1/1 region success。
- paired pilot output為`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/matched_normal_pilot/H2009/output`；reads.tsv實際為tumor `77 ALT / 3 REF`、normal `2 ALT / 15 REF / 14 UNKNOWN`。這是normal genetic background可影響FP解釋的正例，不是全體比例。

## Claim Ladder（red-team lock）

1. stable focal-ALT methyl multigroup：只支持 A-carrier read-level epigenetic heterogeneity。
2. methyl group 與至少一個 partner genotype 經 FDR/confound存活：支持 local methyl-genetic co-segregation。
3. 至少兩個獨立 markers、完整四狀態/CN/CCF相容、跨平台或樣本重現：可稱 bulk subclone-compatible candidate。
4. single-cell/colony/spatial/longitudinal同時確認：才可稱 cellular subclone或 lineage-supported。

現有 bulk long-read資料即使通過第 2 層，也不得把 methyl group數當 clone數、把 HP tag當 clone label，或把 local linear-compatible state直接寫成已證明 linear evolution。

## 2026-07-15 strict robustness audit 方法稽核

- 外部 subagent 初審判定 strict 舊設計為 `NO-GO`：候選後 BH 不能宣稱 FDR、modal-K representative seed 造成事後選擇、`ARI>=0.8` 不足以鎖定 frozen genetic association 使用的分群、row-circular null 排除 identity shift、leave-one-out 可能以替代 read 維持群數，以及 candidate/artifact 身分契約不足。
- 修正版預設只接 `genetically_anchored_multi_marker_candidate_by_sensitivity`（上游 global BY-sensitive genetic gate），strict 的 BH/BY 改名為 `postselection_*_descriptive` 並明示 `fdr_calibrated=false`；科學結果欄位為 `strict_null_robustness_pass`，不再輸出 `strict_final` 或 FDR gate。
- 每種 null 以單一預先固定 seed做推論與 empirical p；另跑多 seed只作穩定性稽核。screen vs strict、column vs row-circular、以及全部 multi-seed結果均要求 label-invariant 的**完全相同 partition**，ARI只留描述。
- leave-one-out 每次移除一條 read後，原始每群仍須 `>=3`，且剩餘 reads必須得到完全相同 partition；因此原 size=3群被移除任一成員時必然 fail，不會由其他 read補位。
- row-circular null從完整循環群 `shift=0..m-1`抽樣，納入 identity shift；不再從`1..m-1`條件式排除零位移。
- stable assignment 現保存 `reads.tsv`、`distance/BERNOULLI/matrix.csv`、`methylation/methylation.csv` 的 absolute path、size與SHA-256；strict另要求 schema、screen contract、REF/ALT、ALT readset digest完全一致，任一不符在建立輸出目錄前 hard fail。
- execution與scientific status已分離：summary `pass=true`只表示 execution integrity，`pass_semantics=execution_integrity_only_not_scientific_confirmation`；零上游候選是合法、可稽核結果。
- 修正後首輪 topic tests=`81 passed`，shared focal null library=`8 passed`；strict專屬=`12 passed`。外部 reviewer再審已判 `blocking=0/high=0, GO for robustness audit`，但再次鎖定 FDR/subclone confirmation為NO-GO。
- reviewer留下的3個medium均已補：CLI不再允許覆寫BY-sensitive selection；synthetic regression固定要求正向`ROBUSTNESS_PASS_NOT_FDR_CALIBRATED`；run manifest加入strict script與shared library SHA-256。另將ALT membership/index改為O(n)查表。正式執行仍須使用鎖定的default selection，且最終只報robustness。
- closure-only再審結果：指定strict範圍 `blocking=0/high=0/medium=0/low=0`，strict `13/13`與shared null `8/8`；正式run判定GO。這個GO只代表方法可執行，FDR-confirmed subclone與lineage claim仍為NO-GO。

## 2026-07-15 v2 red-team 修正與資源決策

- 本節與 `claim-contract-v2.md` **取代**上方早期 Claim Ladder 第 3 層的「獨立 markers／subclone-compatible」用語。正式欄位只使用 `spatially_separated_markers_20bp`、`multi-marker molecular-haplotype base candidate`；在 L1 正交 cellular evidence 前禁止稱 subclone-compatible。
- M1 gate 已固定為 `coarse_ng>=2 AND not unstable AND modal_assignment_ari_min>=0.8`。只穩定重現群數、但 read membership 不穩定者不再算 M1。site row、stable assignment detail與 legacy alias均套相同 gate；legacy `strict_confirm_candidate` 明示不是 R1 claim。
- output completeness改以完整 `(dataset,chrom,pos,ref,alt)` key核對；artifact未保存 allele時，只容許 frozen manifest 中該位置恰有一個 allele，否則 fail closed。receipt同時驗 sample、exit/pass、VCF path/hash、binary path/hash、reference、output directory與 validation。
- existing-output reuse在讀取任何舊結果前先驗原始 receipt provenance，成功時直接回傳原 receipt，不重寫或重新簽發。COLO829真實 receipt測試前後 SHA-256皆為 `3592f8667b1b2adf2f3b02f50a816dca169ff0b752eda8b2e57132ca8885aeec`。
- `latest_tag_projected_join_audit.json` 的既有真實回歸只驗 **HP**；它不能支持全 sSNV PS completeness。正式 all-site analyzer仍須在 ALT selection前對每一 output read做 latest HP/PS exact join，並由 full-scope run manifest保存 join counts/status。
- 上述 screen/provenance修補的指定回歸測試為 `29 passed in 1.07s`。這只驗程式契約；469,849-site正式結果仍以 full run與post-run reconciliation為準。
- 主機有48 logical CPUs與約475 GiB available memory；v2 extraction以7 datasets並行、每dataset 6 threads，峰值42 worker threads。2026-07-15 09:34共享load約87，剩H2009單dataset執行時不重啟為更高thread數，避免丟失已完成進度及加劇共享BAM I/O。後段Python分析預設最多32 processes，並固定 `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`；啟動前保留`/big7_disk`至少300 GiB free gate。
- tree/input獨立 receipt：`results/latest_tree_input_contract_audit.json`，SHA-256=`526f89fc1d5b87c052c6a48fb70820aa87349ee8508ba4f938bd864328678d42`。7/7 `tree_vcf`與同sample `longphase_recalibrated_pass_vcf`的path、size、recorded/observed SHA均相同；top-level 7 gates與逐sample全部gates PASS。469,849 focal universe明確是同一 PASS backbone 的chr1-22 biallelic-sSNV frozen subset。重跑同一路徑已驗證 `FileExistsError` fail closed。

## 2026-07-15 真實資料整合契約補強

- 首次 HCC1937 terminal-tag smoke 在啟動前 fail closed：focal VCF/truth/output皆為 `18,690`，但完整 LongPhase-S PASS ledger為 `52,115`。這不是漏位點，而是 ledger 合法包含非 focal scope context。失敗輸出未刪除，已移至 `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/archive/20260715/screen_smoke_HCC1937_terminal_tag_v2_failed_superset_contract/`；incident receipt為 `results/screen_smoke_v2_failed_superset_contract.json`，明示不得作科學結果。
- focal key contract改為：VCF keys必須與truth keys及output positions exact equal，且必須是 LongPhase-S PASS ledger keys的subset；ledger extra context允許但要計數。缺少任一 focal key仍 hard fail。對應 analyzer schema升為`1.2.0`，指定回歸測試`20/20` PASS。
- terminal latest-tag audit以 **site-read observation** 為計數單位，同一read跨多焦點會重複計數，不能稱全域unique reads。每列強制`fetched >= eligible >= joined == reads.tsv rows`、PS/replaced/multimatch不超過joined、projection multimatch必為0；summary與run manifest都保存`latest_hp_ps_terminal_join_audit`。正式469,849-site輸出完成前，不以FP-only HP audit代替全site HP/PS completeness。
- 修正版真實 smoke輸出落點為 `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/screen_smoke_HCC1937_terminal_tag_v3/`。它只作工程驗證；只有全18,690位點完成且summary/run manifest terminal audit PASS後才可引用其契約，不得把中途stable count當正式比例。
- CN/CCF annotator對非零候選輸出`PASS/EXECUTION_PASS`；零候選輸出`NOT_APPLICABLE`、`reason=ZERO_SELECTED_CANDIDATES`、`is_negative_result=false`、`c1_formed=false`與header-only TSV。獨立測試`12/12` PASS。focal-only CN註記本身不能建立joint witness pair multiplicity，因此不得自動形成C1。
- matched-normal runner對明確selection欄位的零命中只建立`not_applicable_receipt.json`，含Task Type B、六項SHA-256 provenance與`pass=true` execution semantics；不建立sample output、不啟動C++、不寫一般`run_receipt.json`。缺selection欄仍 hard fail。主agent獨立重跑`19/19` tests PASS。
- HCC1937 terminal-tag v3 smoke已完整結束，scope明示`full_469849=false`：`18,690/18,690` sites、TSV `18,691` lines（含header）、stable assignment JSONL `1,938` lines，summary與run manifest皆`EXECUTION_PASS`。這些比例僅作工程驗證，不能替代7-dataset正式結果。
- smoke terminal audit實際為`2,058,400/2,058,400` reads.tsv site-read observations exact HP/PS join，PS present=`1,855,466`、source HP replaced=`1,945,028`、sidecar fetched/eligible=`3,672,484/3,219,976`、projection multimatch=`0`、zero-read sites=`0`。計數單位明示不是globally unique read。
- 修正後專題整體測試：topic tests加root cooccurrence tests合計`170 passed`；僅2個SciPy呼叫NumPy deprecated API warning，無功能失敗。final integration builder另以跨producer/zero-path組合`18/18`獨立重跑PASS並通過`py_compile`。

## 2026-07-15 final-builder 零候選語意 red-team closure

- 唯讀 explorer 對 `build_all_ssnv_final_report_dataset.py` 做負向 fixture 稽核，發現 strict 零候選 receipt 可被竄改成 scientific negative semantics、matched-normal 的 not-evaluable 可被竄改成 negative result而仍通過；另發現 strict/CN 零候選整合狀態可能落為`NOT_RUN`。初始判定為兩項High blocker。
- strict producer的`NOT_APPLICABLE` receipt現同時固定`pass_semantics=execution_integrity_only_not_scientific_confirmation`、top-level與nested `is_negative_result=false`；final builder逐欄強制驗證，並回傳`NOT_APPLICABLE_VALIDATED_RECEIPT`，不再降級成`NOT_RUN`。
- matched-normal summary與receipt現都固定`not_evaluable_is_negative_result=false`，builder同時驗兩份。CN/CCF在零G2且未提供目錄時改為`NOT_APPLICABLE_ZERO_SELECTED_CANDIDATES_OMITTED`；正式流程仍會執行native annotator，保留header-only TSV與receipt作較強證據。
- topic-local `scripts/ssnv_cooccurrence_lib.py`已改為repo canonical `InterSubMod/scripts/ssnv_cooccurrence_lib.py`的compatibility bridge，正式producer與局部測試不再各自執行不同four-state enum。
- regression command：`/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m pytest -q research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_build_all_ssnv_final_report_dataset.py research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_strict_methyl_candidate_confirmation.py research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_matched_normal_candidate_controls.py research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_ssnv_cooccurrence_lib.py`；結果`68 passed, 2 warnings in 3.23s`。
- four-state個案驗算明確分層：`RR=6, AR=6, RA=0, AA=6`雖符合使用者描述的幾何，但在2% error ceiling/95% confidence下為`NOT_IDENTIFIABLE_FIXED_ERROR_CEILING`；零violation分母達預定`149`後才可列`FOCAL_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL`，仍不是cellular ancestry證明。
- 修補後完整方法測試分兩個pytest invocation執行，避免topic/root同名module collection衝突：topic `176 passed, 2 warnings in 9.25s`，repo canonical cooccurrence `16 passed in 0.95s`，合計`192 passed`。report artifact builder另重跑`16 passed in 1.08s`。

## 2026-07-15 v2 全量抽取完成與 reconciliation 進行中

- v2 runner於七樣本全部結束後寫出`results/all_ssnv_intersubmod_batch.v2_verification_fix.json`；實際總計為`expected_vcf_sites=469,849`、`reads_files=469,849`、`bernoulli_matrix_files=469,849`、`methylation_files=469,849`、`failures=[]`、`pass=true`。每個dataset各有獨立`run_receipt.json`與exit=0。
- `significance_summary.csv`總列數為`469,109`，比frozen位點少`740`；分布為HCC1395 `49`、DORADO `83`、COLO829 `271`、H1437 `32`、H2009 `150`、HCC1937 `119`、HCC1954 `36`。runner將此summary視為非權威衍生表，只要求`0 < rows <= expected`；正式screen不使用它作位點母體。
- 740-row差額在獨立exact reconciliation完成前維持`PENDING_INTERPRETATION`。正在執行的command為`python3 research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/audit_all_ssnv_output_completeness.py --manifest research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_input_manifest.json --output-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/intersubmod_all_ssnv_v2_verification_fix --output research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_output_reconciliation.v2_verification_fix.json`。
- exact audit的hard gate為三類per-site artifacts各自完整五欄key reconciliation，且missing/extra/duplicate/empty均為0；未通過前不得啟動正式全量screen或發布科學比例。

## 2026-07-15 正式 screen 與 portable report QA

- exact reconciliation已完成：`results/all_ssnv_output_reconciliation.v2_verification_fix.json`對469,849個五欄dataset-site keys逐一核對，reads/methylation/BERNOULLI的missing、extra、duplicate、empty皆為0，`pass=true`。`significance_summary.csv`少740列只屬非權威衍生表，沒有造成正式screen輸入缺位點。
- extraction後 frozen input audit為`results/frozen_input_immutability.post_v2.json`：7/7 samples、77/77 frozen artifacts PASS。matched-normal pre-candidate audit為`results/matched_normal_methyl_tag_audit.v3_pre_candidate.json`：7/7 BAM可讀且每個樣本1,000/1,000 primary reads具MM+ML；所有抽樣normal HP仍為空，符合normal不使用HP的契約。
- 正式screen輸入為`intersubmod_all_ssnv_v2_verification_fix`，輸出為`all_ssnv_focal_alt_multigroup_v3_terminal_tag`；command固定46 workers、chunk-size 32、max-pending 92，並將OMP/OpenBLAS/MKL/NumExpr各鎖為1 thread。主機48 logical CPUs，保留2 cores供OS與稽核；約450 MB RSS/worker，瓶頸為BAM/I/O，故不重啟或再超額增加worker。
- 2026-07-15約19:05 +08:00的中途進度為`207,000/469,849`、stable=`32,711`。這些只是progress log，不是正式分母或科學比例；只有`all_ssnv_summary.json`與`run_manifest.json`同時完成、469,849 exact且terminal HP/PS audit PASS後才可引用。
- report builder的完整Markdown與portable artifact刻意分層：`report.md`保留完整30段敘述、全欄位表格與hash；artifact/HTML只顯示3-block可視化摘要，但仍內含13個bounded datasets、完整table/chart definitions及12個source mappings。`report_build_receipt.json`以`artifact_presentation_scope`明示此差別，避免把HTML短頁誤解成資料刪減。
- 官方Data Analytics synthetic QA：`validate_artifact ok=true`（13 datasets、12 sources）；`deliver_portable_artifact.mjs`的validation/package/verification皆PASS，source dialog鍵盤互動PASS，1440/390 viewports皆無overflow。先前長頁失敗已定位為外部portable runtime在垂直scrollbar時sticky top bar的`100vw`產生8 px overflow；未修改外部renderer，而以正式visual-summary/full-Markdown雙層交付解決。
- report/analysis topic tests最新為`190 passed, 2 warnings in 15.91s`；兩個warning均為SciPy呼叫NumPy deprecated API，無功能失敗。外部Claude Code版本固定為`2.1.202`，待正式結果完成後以read-only工具做結果級終審。

## 2026-07-15 extraction reference identity 證據缺口封閉

- extraction receipts的`reference_sha256=null`不再由路徑相等默認補足。獨立稽核`results/extraction_reference_identity_audit.v1.json`已對`/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta`與`.fai`做完整SHA-256；reference=`9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7`、FAI=`502a1b8fb73ccd53285c28a0f12df90c818b4fe3de1e862ef47c593ef1a0a4b4`。
- 該receipt同時綁定7/7 extraction receipts的reference path/command path、既有frozen size+sampled chunks及FAI hash；所有checks PASS。限制仍明示：這是post-extraction full hash加pre-existing frozen chunk identity，不是倒填每個immutable extraction receipt。
- `build_all_ssnv_report_artifact.py`現將此receipt列為必填`--reference-identity-audit`與第5個release-blocking audit；schema、Task Type B、pass semantics、469,849 scope、7樣本membership、receipt hash/exit/binding、所有checks與limitations任一不符即在建立輸出目錄前fail closed。
- 驗證輸入：`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/extraction_reference_identity_audit.v1.json`。驗證command：`/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m pytest -q research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests`。實際結果：`191 passed, 2 warnings in 16.16s`；warning仍僅為SciPy/NumPy deprecation。真實audit另經同一validator直接載入，輸出`datasets=7, full_reference_sha256=true, pass=true`。

## 2026-07-15 portable report 整體與個案單圖 QA

- 增加第二張chart或metric strip會讓官方portable reader在產生垂直scrollbar後，因其sticky header `100vw`/viewport寬度處理於1440或390觸發`horizontal_overflow`；full/half、4-claim與2-card synthetic均fail，因此未納入正式交付，也未修改外部renderer。
- 正式artifact改用單一共同百分比尺度的native bar chart：左側是M1/M2/G1/G2的PASS/evaluable比例（dataset-site與biological-site兩個series），右側是actual case每個methyl group內partner R/A read composition。tooltip保留panel、原始分子、分母與status；這不把read count或甲基群升級為cellular clone。
- 新bounded dataset為`overview_case_chart`，正向fixture共12 rows（8 claim + 4 case）；來源`src-overview-case-view`明確綁定final schema 2.0.0 dataset與cooccurrence pair table。HTML仍是3個visible blocks，但現在同一張圖同時回答整體與實際個案；完整九層、strata、four-state、joint signature與hash仍保存在`report.md`及artifact datasets/tables/sources。
- 三種presentation path均經官方`deliver_portable_artifact.mjs`通過`validation/package/verification`、source keyboard interaction及1440/390 viewports：G2>0、G2=0但有non-confirming witness、以及無eligible pair的structured N/A。N/A分支不再保留count全null的不可渲染舊case chart definition，summary明示`actual case=N/A`，但N/A rows與來源仍保留。

## 2026-07-15 下游 preflight 與 data-plane SHA closure

- 獨立 explorer `019f65bb-a395-79c2-bbe7-9d53d592cbad` 在正式下游啟動前找到一項 blocking contract drift：cooccurrence producer 的 G2 要求同一 complete-read set 上至少2個 effect-supported、相距至少20 bp的top markers，但 final builder舊重算式漏掉此gate。若合法site其producer selection為false，builder會錯誤要求true並中止。
- `build_all_ssnv_final_report_dataset.py`現把`joint_signature_complete_marker_effect_supported_positions`與producer衍生的same-complete-read count/positions列為必填，逐欄重算intersection、20 bp spacing與G2 selection；新增正反例證明只有1個complete-read marker時合法維持G2=false。targeted final-builder tests=`22 passed`。
- 新增`scripts/audit_stable_primary_artifacts.py`：對完整M1 stable assignment key set與screen site set做exact reconciliation，並對每一site的`reads.tsv`、BERNOULLI matrix、methylation matrix逐檔重算resolved path、size與SHA-256。正式流程會在所有data-plane consumers前後各跑一次，final builder要求兩次`artifact_set_sha256`相同。
- final builder另要求pre audit完成時間早於cooccurrence/tumor-REF開始，post audit開始時間晚於兩者結束；cooccurrence receipt新增`started_at_utc/finished_at_utc`。因此正式`VERIFIED_SHA256`不再只代表producer記錄，而是相同artifact set包住實際data-plane讀取區間。
- matched-normal analyzer現重新掃描runner輸出，依runner相同relative-path+file-SHA演算法重算每sample `artifact_set_sha256`並寫入summary/receipt；final builder逐sample與paired runner receipt比對。竄改任一paired artifact會在analysis前hard fail。matched-normal tests=`20 passed`。
- reviewer指出G2=0時runner只產生`not_applicable_receipt.json`，analyzer不能讀此路徑。正式orchestration已固定條件分支：零G2跳過analyzer，runner目錄直接交final builder；非零才啟動analyzer。這不是未處理錯誤。
- reviewer另指出假設性M1=0分支尚未由tumor-REF/cooccurrence producer原生閉合；本次真實full-scope run在中途已產生超過4萬個M1 stable assignments，因此不會進入該分支。此限制不得改寫成生物負結果，並留作workflow generalization項目，不影響本次469,849-site實際資料完成性。
- 修補後topic tests=`196 passed, 2 warnings`；repo canonical cooccurrence tests=`17 passed`。兩個warning仍僅為SciPy呼叫NumPy deprecated API。4支修改/新增producer均通過`py_compile`。
- 正式screen在H2009前置154,465-site目錄盤點期間曾停在`271,000/469,849`；PID與46 workers全程存活，parent `read_bytes`持續增加。盤點完成後恢復，最新progress為`277,000/469,849`、stable assignments=`42,338`。此為進度紀錄，不是發布比例。

## 2026-07-15 第二輪 data-plane / G2 red-team closure

- 第二輪唯讀紅隊發現前一版primary-artifact audit使用了錯誤的自洽`SCREEN_CONTRACT`字串；這會讓測試與audit彼此一致，卻無法讀取正式screen輸出。現已改為producer canonical `phylo-v4.1_column_null95_modal_K10_RNULL40_min_group3`，並新增AST cross-producer constant test，禁止兩邊一起漂移。
- G2 final integration不再信任`joint_signature_sensitivity_pass`或effect-supported positions等producer衍生布林。builder會由`testable + permutable + p<=0.05`重算joint pass，逐top marker由同一complete-read support的`testable + Cramer's V>=0.30 + delta ALT fraction>=0.50`重算effect gate，並把site `top_marker_positions`與pair-table `top_coverage_marker=true`集合精確對帳。
- pre/post primary-artifact audits均保存`started_at_utc/finished_at_utc/created_at_utc`，要求`started<=finished==created`。pre audit的`consumer_receipts=[]`；cooccurrence與tumor-REF producer各以path/size/SHA-256綁定同一pre audit；post audit再以path/size/SHA-256綁定兩份producer receipt。final builder同時驗證時間包圍與雙向雜湊鏈。
- audit輸出採`os.path.lexists`加exclusive `open("x")`，dangling symlink與競態建立均fail closed；site/assignment/consumer inputs在verification前後重算artifact identity，執行中改動不得產生receipt。
- 指定故障注入測試：primary audit `11 passed`；final integration `27 passed`；tumor-REF `13 passed`；cooccurrence root/topic分開收集為`9 + 26 passed`。涵蓋same-size SHA tamper、worker exception、joint p drift、marker effect drift、top-marker跨表漂移、audit時間倒置與matched-normal per-sample digest mismatch。
- 2026-07-15約21:20 +08:00正式screen進度為`323,000/469,849`，仍位於H2009區段。48 logical CPUs上46 workers實測約90% user CPU；資料碟瞬時100% util。此時增加process會造成CPU/I/O oversubscription，且重啟會失去八小時進度，因此維持46 workers是本輪可稽核的最大合理資源配置。
- 上述數字仍是工程進度，不是科學結果。正式M1/M2/G1/G2比例、個案與結論只可由完成後的summary/run manifest、pre/post audit、cooccurrence/control outputs與final integration dataset計算。

## 2026-07-15 producer BAM-output 報告 gate 補強

- 最終report builder現在會從frozen layered manifest取出7份`producer_capture_receipt_v2.json`，逐份重算path、size與SHA-256，並要求schema、version及sample identity完全一致。
- 每份receipt都必須同時滿足`transport=named_fifo`、`persisted_bam=false`、`regular_bam_count=0`與`is_fifo_at_closeout=true`；任一欄漂移均在建立report output directory前fail closed。
- 因此正式報告可明確寫：LongPhase-S的BAM-named output只是即時消費的FIFO，沒有持久化regular tagged BAM，也沒有以新檔覆寫原始或既有tagged BAM。大型raw BAM不宣稱full-byte hash；輸入未變仍受frozen metadata+sampled-chunk assurance與post-run audit限制。
- Report flow、release-blocking audit、source inventory與artifact snapshot均加入此契約；7份producer receipt各自保留可追溯SHA-256。新增故障注入把重新簽署但`persisted_bam=true`的receipt判為fail。
- Report專屬測試為`30 passed`；本專題全套為`212 passed, 2 warnings`。兩個warning均為SciPy呼叫NumPy deprecated API，沒有功能測試失敗。

## 2026-07-16 高深度 tail 診斷與獨立統計審查

- 正式screen已完成至少`467,000/469,849`的ordered輸出，但HCC1954尾端一個32-site chunk長時間未回傳。parent與46 workers持續存活；慢worker CPU time逐秒增加、無OOM或exception，因此不是deadlock。
- 只讀重算`HCC1954 chr17:39851484`：原矩陣4,752 reads、focal ALT 4,673、complete-distance peel後4,610；其樹在只套`SEP_MIN`的上限路徑可有82個null-tested nodes、深度74、平方成本65.5588 root-equivalents。這解釋單worker數十小時尾端，而非資料遺失。
- `pre-decision-audit.md`新增runtime rescue addendum，判定`GO_WITH_EQUIVALENCE_GATES`：只平行11個彼此獨立的coarse/fine runs，不更動每run seed、null順序或遞迴。先做exact-equivalence，再以新目錄重跑HCC1954，最後用前6 datasets完整prefix加replacement HCC1954做469,849-key merge；原partial輸出不覆寫、不刪除。
- 第一位獨立統計agent判定下游發布前有4項High gate：strict低effect應列FAIL而非N/A、final builder須完整重算G1、four-state須由RR/AR/RA/AA重算95% upper bound與2% ceiling、HCC1395 biological rate不得採兩平台any-PASS union；另需收緊G2 evaluability。
- M1 `insufficient_valid_null`與真正screen negative混合是conditional blocker。正式報告若未補齊全site null-calibratability，不得把legacy `stable_fraction_evaluable`稱prevalence；至少報全469,849的operational flag yield，並把conditional null-calibrated比例標為not estimated。是否需要重算由完成output可取得的failure census與方法claim決定。

## 2026-07-16 正式 recovery、thread pinning 與 terminal chain

- HCC1954 nested-parallel 與原 serial producer 已完成 real high-depth site及synthetic payload的exact-equivalence。正式receipt為`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/phylo_parallel_exact_equivalence.pinned_390228ce.v1.json`，`pass=true`；real site為`HCC1954 chr5:832126`、`n_reads=368`、ALT peel後230 reads，serial/parallel完整row SHA-256皆為`09fd9a7b2a19dc39318ce525858f0e3313f604a68b18df7ef760d201ccb3667b`、detail SHA-256皆為`3c0db58b172d4bbcd5fd10ba7d1a68ec95923cd40a05b0076a83f0be48a83e4a`。正式 analyzer已固定為`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/analyze_all_ssnv_focal_alt_multigroup.pinned_390228ce.py`，完整SHA-256=`390228ce1fb98b59409f2517341334ab94329ef6c73a4838555fca85878027b2`且mode=`0444`。
- 原全量tail run與效能探索run均保留但不得發布：舊full partial、v4約68,820、v5 threshold-400，以及first-six v7的13,000-site partial皆屬nonformal。v7被中止的工程原因是44個外層worker各自啟動48條BLAS threads，約2,112 threads造成load >1,100與高system CPU；其前13,000位點stable=`1,071`只作重現檢查，不可作比例。
- 正式first-six替代run為`all_ssnv_focal_alt_multigroup_v9_first_six_thread_pinned_source_locked`：40外層workers、phylo inner=1，並固定`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=NUMEXPR_NUM_THREADS=BLIS_NUM_THREADS=1`。逐worker實測`NLWP=1`、約86-97% CPU；first 13,000位點stable同為`1,071`，與非正式v7前綴一致。正式HCC1954替代run為`all_ssnv_focal_alt_multigroup_v6_hcc1954_seed_parallel_200`：3外層workers x 11 nested workers，只在`n_reads>=200`啟用exact-equivalent路徑。
- 2026-07-16本段更新時，first-six為`19,000/447,449`、stable=`1,755`；HCC1954為`9,000/22,400`、stable=`1,982`。這些仍是progress counters，不是科學分母。資料碟瞬時util約96.8%，40個first-six workers已近滿CPU；增加worker會擴大random-I/O競爭，因此維持現配置是本機48 logical CPUs的最大合理用量。
- 完成後唯一正式screen將由`merge_all_ssnv_screen_recovery.py`合併v9 first-six與v6 HCC1954，落在`all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full`。merger要求prefix source-lock receipt、replacement exact-equivalence receipt、469,849五欄key exact、stable assignment exact及所有producer provenance；不會覆寫任一partial run。
- 新增可審核terminal runner：`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_terminal_validation_chain.sh`。固定順序為recovery merge -> frozen audit -> primary artifact pre-audit -> tumor-REF -> all-PASS sSNV cooccurrence -> strict robustness -> matched-normal -> CN/CCF annotation -> consumer/post audit -> terminal frozen audit -> final machine dataset -> Markdown/artifact -> portable HTML/desktop/mobile QA。matched-normal只讀`multi_marker_molecular_haplotype_base_candidate=true`；零G2時保留native N/A receipt並跳過不適用的normal analyzer。
- terminal runner所有新output path均要求不存在，BLAS相關thread環境全固定為1，post primary audit的consumer順序固定為cooccurrence、tumor-REF、strict、matched-normal paired、matched-normal analysis（若適用）。report使用downstream全部完成後的`frozen_input_immutability.post_terminal_downstream_v1.json`，避免只證明screen完成時未變。
- portable report現為完整36 blocks、4 charts、15 tables的canonical artifact，不再沿用舊3-block摘要限制。官方reader在non-overlay scrollbar下的`100vw` top bar問題，以專題wrapper只注入兩個header selector的`width:100%`修正，再呼叫官方builder/chart extractor/verifier；synthetic fixture desktop `1425/1425`與mobile `375/375` document width、console errors=0。
- `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/qa_portable_report_layout.py`新增`.analytics-layout-row` pairwise overlap檢查；synthetic desktop/mobile均`overlapCount=0`、`pass=true`。正式terminal gate同時要求官方delivery receipt `pass=true`、desktop/mobile page overflow/console gate PASS與`overlapCount=0`。
- 正式數值結論仍保持未發布：只有v10 screen、pre/post audits、cooccurrence/strict/normal/CN outputs、final dataset與result-level multiagent/Claude稽核全部完成後，才可填入M1/M2/G1/G2比例與個案敘述。

## 2026-07-16 M2/B1 獨立稽核修正

- 第一輪獨立唯讀方法學review指出兩個發布阻塞點：舊M2只信任`*_aligned`布林，因此效應已達門檻但permutation p未達0.05的低資訊軸可能被當成「未對齊」；舊B1則可在多個G1正式pair中事後挑任一four-state compatible pair，且把`NOT_IDENTIFIABLE_FIXED_ERROR_CEILING`錯列FAIL。
- 新共用gate為`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/m2_screen_gate.py`。八軸逐一重算並保存四態：`ALIGNED_EFFECT_AND_PERMUTATION_P_PASS`、`NOT_ALIGNED_EFFECT_BELOW_THRESHOLD`、`NOT_ALIGNED_AXIS_HAS_NO_OBSERVED_VARIATION`、`INDETERMINATE_EFFECT_ABOVE_THRESHOLD_WITHOUT_P_LT_0_05`。最後一態令M2=`NOT_EVALUABLE`，不得進入G1/G2 family；無變異軸因不可能解釋當前分群，可保留為determinately not aligned。
- cooccurrence site TSV新增gate contract、evaluable/eligible、每軸完整證據與indeterminate/aligned/constant axis清單。final builder由原screen數值獨立重算，逐欄比對cooccurrence輸出；receipt另鎖定M2 gate程式SHA，final整合再次驗證目前source identity。
- B1現在只使用G1 formal pairs中`endpoint_a_n_informative`最大、再以abs distance與partner identity固定tie-break的單一pair；排序完全不讀four-state結果。只有該pair compatible才PASS；insufficient depth或fixed-error ceiling均`NOT_EVALUABLE`；明確`INCOMPATIBLE_OR_COMPLEX`才FAIL。所有compatible pairs仍可作descriptive opportunities，但不參與B1 gate。
- M1正式名稱改為operational stable methyl-multigroup screen flag；非stable的FAIL只代表未被此screen標記，不能當作同質性證據。final dataset另保存`global_null_validity_exported_for_nonstable_sites=false`、biological prevalence=`null`與不可區分true negative/null-invalid的限制。
- HCC1395/DORADO aggregate技術重現不再輸出`PASS/FAIL`，改為`ANY_CONCORDANT_EXACT_PAIR_OBSERVED`等描述狀態；固定`replication_claim_status=NOT_EVALUABLE_BIOLOGICAL_N1`、無inferential CI且pair independence=false。site metric的truth明確只指focal；另新增9-cell focal x partner truth matrix、table與stacked chart。
- 關鍵回歸包括：高效應但p不顯著軸必為M2 NOT_EVALUABLE；constant軸可明確判為不可能對齊；有另一個compatible pair時，若預先指定pair為fixed-error ceiling，B1仍須NOT_EVALUABLE，禁止事後換pair。final integration=`36 passed`（新增source identity測試後待全套重跑）；report artifact=`32 passed`；terminal chain=`5 passed`。

## 2026-07-16 terminal contract v3與machine TSV provenance closure

- fresh-context integration red-team `019f68c8-88bf-71b3-9a50-0c2aa3d416af` 找到2項High與2項Medium：舊`claim-contract-v2.md`缺新版M1名稱且會令report builder中止；v2仍描述舊G2/B1/149-depth語意；兩份machine TSV漏新版M2/G2/four-state provenance；terminal runner任一步驟失敗後無法保留已完成階段續跑。
- screen-time source-lock未被修改。修補前後四個SHA-256均固定：pinned analyzer=`390228ce1fb98b59409f2517341334ab94329ef6c73a4838555fca85878027b2`、focal library=`ed5f3c99461248248b20a9f49597ab5de7340a1e2055d77a1d83dbcc2799b72a`、latest-tag join=`7ddadb05481bb61f2dc3a489762944ec347ab9f5265c69fb17538efd44b543cd`、`claim-contract-v2.md`=`904abd6c9fbcf152770be72a3cd2b12f38d0058b4d668b2c6ced9005b86afc22`。
- 新增`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v3.md`作post-screen、pre-terminal integration contract。它明示v2仍是screen source-lock歷史契約，v3只收緊terminal統計與敘述：M1 operational flag、M2 80% planning-power evaluability、G1完整pair family global BY、G2 joint完整family global BY、B1 outcome-blind單一預選pair、3 relation Bonferroni、familywise 95%與zero-violation depth 203。
- report builder不再只grep claim IDs/names；它同時固定完整v3 SHA-256=`da94a50d0717174ff007b75f2edad2de79bf3aebf6b15df179eb736e8d8f526e`，任何門檻、denominator或claim ceiling文字漂移都會在輸出建立前hard fail。真實v3直接驗證輸出為`PASS da94a50d...8f526e`。terminal runner只在report階段讀v3；screen recovery仍經prefix source-lock receipt保存v2 identity。
- `candidate_catalog.tsv`補`m2_low_power_axes`與5個joint global FDR欄位；`candidate_witness_pairs.tsv`補`four_state_familywise_confidence`、`four_state_relation_family_size`與`four_state_multiplicity_method`。report builder現在要求final receipt逐一SHA綁定兩份TSV，並逐key/逐required-field與final JSON交叉核對；刪除、未重新簽署竄改、重新簽署但與JSON不一致三類故障注入皆fail closed。
- 第二位fresh-context reviewer證明只檢查`pass=true`的初版`--resume`可接受偽receipt，因此該功能已完全撤回。terminal runner維持one-shot、所有output預先不存在與不覆寫；若真實鏈失敗，保留既有輸出後以新tag或另經完整schema/SHA驗證的續跑設計處理，不以便利性降低證據標準。
- 新增root `pytest.ini`固定`--import-mode=importlib`，同一invocation收集topic與root兩組同名cooccurrence tests不再發生module mismatch。最終相關聯合命令為`/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m pytest -q tests/test_ssnv_cooccurrence_lib.py tests/test_analyze_methyl_ssnv_cooccurrence.py research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests`，實際`272 passed, 2 warnings in 15.44s`；warning僅SciPy呼叫NumPy deprecated API。全repo Python collection另受該env缺`seaborn`與`statsmodels`的5個既有dependency errors阻擋，已明確分開，不誤報成方法測試失敗。`bash -n`與`git diff --check`均exit 0。
- 最後一位fresh-context closure auditor `019f68ff-0a3a-7280-9039-fa66d93b8c96` 對A-D做限定故障注入後判定全部CLOSED、無可重現繞過：one-shot overwrite gate、v3 full-SHA gate、兩份TSV未簽署/重新簽署drift與root/topic同名測試同批收集皆PASS。此closure只涵蓋方法與integration契約；正式469,849-site final dataset仍須等active screen與terminal chain完成，不能先行標科學結果CLOSED。
- 本段完成時source-locked first-six progress=`156,000/447,449`；HCC1954仍在高read-depth nested chunks。兩者均保持active，這些progress/stable counters仍不可作科學比例。

## 2026-07-16 M1 全 dataset-site 分母契約修正

- 服務研究目標：G3（read-level epigenetic 的可解釋性）與 G5（外部可驗證的統計口徑）。Task Type 維持 B comprehensive validation，正式母體固定為 7 datasets、chr1-22、469,849 個 LongPhase-S recalibrated FILTER=PASS biallelic focal dataset-sites。
- 在 terminal chain 啟動前重新對照 `claim-contract-v3.md`，發現 final dataset builder 舊實作仍以 screen-evaluable 子集作 M1 主分母；這會把 operational screen yield 與技術可評估比例混為一談，屬 publication blocker。執行中的 source-locked screen、輸入與輸出均未修改。
- `scripts/build_all_ssnv_final_report_dataset.py` 現令每個 dataset-site 的 M1 僅有 operational `PASS=FLAGGED` 或 `FAIL=NOT_FLAGGED`，主分母必須為 469,849；另以 `m1_operational_screen` 保存 `n_screen_evaluable`、`n_screen_not_evaluable`、全位點 flag yield 與可評估子集次要比例。not-flagged 不得解釋為生物陰性，`biological_prevalence_estimate` 固定為 null。
- `scripts/build_all_ssnv_report_artifact.py` 在建立任何報告輸出前，會獨立驗證 M1 pooled population/denominator=469,849、`NOT_EVALUABLE=0`、flagged/not-flagged 與技術 evaluability 兩組計數各自守恆、兩個比例可重算、正式 denominator definition 完全一致。報告文字會同時顯示全位點主比例及可評估子集次要比例，並明示 M2 以後才保留 claim-level `NOT_EVALUABLE/NOT_RUN`，避免用通則誤描述 M1。
- 驗證輸入：上述兩支 builder、`claim-contract-v3.md` 與兩份對應 test。驗證命令：`/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m pytest -q tests/test_ssnv_cooccurrence_lib.py tests/test_analyze_methyl_ssnv_cooccurrence.py research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests`。實際輸出：`273 passed, 2 warnings in 18.99s`；warning 僅為 SciPy 呼叫 NumPy deprecated API。
- screen-time source locks 重算後未漂移：pinned analyzer=`390228ce1fb98b59409f2517341334ab94329ef6c73a4838555fca85878027b2`、focal library=`ed5f3c99461248248b20a9f49597ab5de7340a1e2055d77a1d83dbcc2799b72a`、latest-tag join=`7ddadb05481bb61f2dc3a489762944ec347ab9f5265c69fb17538efd44b543cd`、screen-time claim contract v2=`904abd6c9fbcf152770be72a3cd2b12f38d0058b4d668b2c6ced9005b86afc22`；terminal v3=`da94a50d0717174ff007b75f2edad2de79bf3aebf6b15df179eb736e8d8f526e`。
- 獨立唯讀 reviewer `019f696a-f290-7f73-8f6d-b94c9ca09536` 確認 machine dataset 的 M1 主分母與 visible Markdown 已正確，但找到 canonical artifact 的中度發布 blocker：overview/claim chart 與 source metadata 仍把 M1-M2-G1-G2 統稱 `PASS/evaluable`，且 M1 biological-site bar 未在圖內明示只屬描述性 aggregation。這會讓正確數值搭配錯誤語意，故初審 verdict 為 NO-GO。
- 修正後 M1 的 chart series、subtitle、tooltip、source description 與 metric definitions 全部固定為 `FLAGGED / all dataset-sites`；M1 biological-site ratio 固定標為 descriptive aggregation/non-prevalence；只有 M2-L2 使用 `PASS / claim-specific evaluable denominator`。原本要求舊字串的測試已反轉，並直接驗證 M1 兩種 series 與 source metric definitions。專項結果 `79 passed in 9.96s`，完整相關聯合結果 `273 passed, 2 warnings in 18.29s`。
- 原 reviewer 完成 bounded re-review：無 High/Medium/Low finding，synthetic artifact 的 M1 dataset series=`10/469,849` 且定義為 `FLAGGED / all dataset-sites`，biological series 明示 descriptive/non-prevalence，M2 才使用 claim-specific evaluable denominator；其獨立 targeted tests=`79 passed in 9.09s`，v3 hash未漂移。最終判定為 GO，僅涵蓋本次 M1 report-layer blocker closure。

## 2026-07-17 tumor-REF 長尾與背景重現 gate 語意封閉

- terminal chain 的 tumor-REF stage 已完成 ordered output `102,000/102,842`，最後一個 HCC1954 chr14 16-site chunk 由單一 worker 持續計算。worker CPU time 等量增加、RSS 約 3.3–3.5 GiB、無 OOM/exception；其餘 39 workers 已完成 pending chunks但受 ordered future 的 head-of-line blocking。`phylo_label` 每個 site 需 REF 與 joint 兩個subset，各做10 coarse+1 fine、每個可測node 40個column null；因此判定為有限高乘法成本，不是資料遺失或無限遞迴。不得 timeout 或排除長尾位點。
- 獨立唯讀 reviewer `019f6d29-475a-7403-bc35-9343efcba173` 找到背景控制與primary M1的穩定性語意不完全對稱：tumor-REF與matched-normal REF producer以`coarse_ng>=2 && !unstable`判重現，沒有primary M1額外的`modal_assignment_ari_min>=0.8`。在**同一個背景payload**上，較寬鬆predicate是ARI-qualified predicate的superset；這不是宣稱ALT、tumor-REF與normal-REF三組實際flag集合互為superset。B1只在較寬鬆背景predicate仍不重現時通過，因此相對於同一payload的ARI-qualified判定不會增加B1 PASS，只可能把membership不穩定的背景K結構保守地當成重現而減少候選。
- active tumor-REF producer、輸入與輸出均未修改。final machine dataset新增`background_control_replication_gate`，契約固定為`lenient_coarse_modal_K2_without_membership_ARI_requirement_v1`，保存scope、no-ARI disclosure、同一背景payload上的predicate包含關係及false-positive/false-negative方向。B1 reason改用`LENIENT_COARSE_MODAL`字串；Markdown/HTML方法與限制明示這不是完全對稱的primary-M1 replay。
- 執行期source identity補充證據：tumor-REF analyzer SHA-256=`95bf7cdca5b636eb0905693ee3c35f1bab699cf4872370fe6f22157cac4c8b87`，mtime/ctime=`2026-07-15 21:55:20 +08:00`；focal library SHA-256=`ed5f3c99461248248b20a9f49597ab5de7340a1e2055d77a1d83dbcc2799b72a`，mtime/ctime=`2026-07-15 08:26:55 +08:00`。兩者ctime均早於terminal launch `2026-07-16 20:03:56 +08:00`，完成後仍須與producer receipt hash再次核對。這是本輪posthoc source-identity佐證，不取代未來應新增的prelaunch source-lock receipt。
- 新增`audit_retrospective_running_source_identity.py`後，PID 2096992存活期間建立不可覆寫snapshot；除SHA-256外也鎖定resolved path、device、inode、size、mode、mtime_ns、ctime_ns、live cmdline與process start。完成後必再與producer `run_manifest.json`的source artifacts、command與時間區間逐項核對，receipt並hash-bind正式manifest。其限制明示為bounded retrospective source-file identity，不是prelaunch lock或完整環境attestation；正向/漂移拒絕測試=`2 passed`。
- 驗證命令：`/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m pytest -q tests/test_ssnv_cooccurrence_lib.py tests/test_analyze_methyl_ssnv_cooccurrence.py research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests`。實際結果：`278 passed, 2 warnings in 15.43s`；warning仍只來自SciPy呼叫NumPy deprecated API。專項final/report builders為`82 passed in 8.36s`。

## 2026-07-17 M2 2-10群邊界、contract v4與source-attested release

- 全量v10 screen已正式完成469,849 dataset-sites，M1 operational stable multigroup共102,842。獨立M2 gate重算對所有102,842個M1位點逐列重建結果：eligible=`919`、evaluable-ineligible=`948`、axis-indeterminate NOT_EVALUABLE=`100,974`、group-count>10 NOT_EVALUABLE=`1`，四類合計精確等於102,842。正式receipt為`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/independent_m2_gate_recount.v1.json`，SHA-256=`a4a3231ea00c684d90b6a83bed5c505ebfcf1c42a51d714978c95148ecb4f9e7`、全部5項checks PASS。這些是M2 screen-gate狀態，不是G1/G2 discovery或cellular subclone比例。
- M1群數分布為K=2 `76,402`、K=3 `19,003`、K=4 `5,430`、K=5 `1,471`、K=6 `398`、K=7 `93`、K=8 `33`、K=9 `10`、K=10 `1`、K=11 `1`。唯一超過planning model上限的位點是`HCC1954 chr5:751076 G>A`，`n_alt_after_peel=119`、11個核心群；正確machine decision為M1 PASS、M2 NOT_EVALUABLE、G1/G2/B1 NOT_RUN，不得擴張模型或納入M2 PASS/FAIL分母。
- `m2_screen_gate.py`現在可一致接收JSON字串或原生Mapping形式的`cluster_sizes`；11群在任何axis power計算前結構化返回`NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM:observed=11:maximum=10`。final dataset另保存M2 NOT_EVALUABLE reason census、2-10群分母定義、超群案例及下游claim行為；report machine validator逐項守恆並在Markdown/portable HTML顯示專表。
- `claim-contract-v3.md`保持SHA-256=`da94a50d0717174ff007b75f2edad2de79bf3aebf6b15df179eb736e8d8f526e`不變。新增不可回寫的`claim-contract-v4.md`，SHA-256=`3bcc2d042708e34ec7c75db7fe9c65baef9b80494f668f911f9551ef0e588515`；它只封閉v3既有2-10群上限與實際background lenient predicate，不改變已執行screen或統計。Intermediate terminal build只能使用exact v3；tumor-REF source-attestation PASS的Task Type B release必須使用exact v4，交叉使用會fail closed。
- Tumor-REF執行期間的權威v3 source snapshot為`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tumor_ref_source_identity_audit_v3/observed_during_execution.snapshot.json`，SHA-256=`8e71cf2e0236f82017186b662149095eb8315a85a6a7212d72233387af59f159`、mode=`0444`。完成後`run_source_attested_release_chain.sh`會先建立post-run receipt，再重建v3 source-attested final dataset、v4 report、portable HTML及desktop/mobile QA；原terminal report維持intermediate，不覆寫。
- 最新完整相關測試為`291 passed, 2 warnings in 16.48s`；新增v4 release與11群路徑的定向測試為`96 passed in 8.62s`。兩個warning仍是SciPy呼叫NumPy deprecated API；Python `py_compile`與terminal/release兩支shell的`bash -n`均exit 0。正式G1/G2/R1/B1/C1/L1/L2數值仍須等待active tumor-REF與terminal downstream全部完成後才能發布。
- 原方法reviewer `019f6d8c-f6f3-7c51-a640-d5ea590d48f8` 完成bounded re-review：前輪4項findings全部CLOSED，High/Medium/Low均為0，最終判定GO；其唯讀定向重跑為`104 passed in 8.49s`。此GO只封閉contract v4 SHA gate、intermediate v3相容、M2 reason conservation、11-group regression、Mapping API與B1同payload敘述，不代表尚未完成的source-attested release或科學結果已發布。

## 2026-07-17 M2 v5語意修正與tumor-REF suffix recovery決策

- 新一輪logic review確認既有程式重算`eligible=919`、`evaluable-ineligible=948`本身一致，但
  `claim-contract-v4.md`把所有低於80% planning power軸一律描述成NOT_EVALUABLE，與實際gate「已觀察到
  effect達門檻且permutation p<0.05時可保守排除」的非對稱語意衝突。121個M2 evaluable rows含至少一個
  aligned-but-below-negative-power軸；若照v4字面錯改為NOT_EVALUABLE，分母會由1,867變1,746且eligible率
  會被事後抬高到52.6346%。此替代不採用。
- 新`claim-contract-v5.md`明示：陽性混雜證據不因事前power不足而消失；只有未對齊軸要形成陰性判定時，
  才要求80% planning power與每群至少5 reads。正式M2仍為919/1,867=49.2234%；這是screen eligibility，
  不是subclone比例。
- M2 categorical缺值不再只靠production library行為推定constant。cooccurrence producer將直接從stable
  assignment的`core_reads.latest_hp`、HP-family mapping與`core_reads.strand`計算observed level counts，並
  驗證`core_reads.label` counts與screen coarse `cluster_sizes`完全一致。observed levels>=2但effect/p缺值、
  或levels=1卻有statistic均hard fail。
- 既有`independent_m2_gate_recount.v1.json`仍保留，但其script曾import production M2 gate，不能稱
  logic-independent。v2 auditor將自行實作常數、power、status與flag reconciliation，不import production gate，
  並exact join全部102,842 stable assignments；正式報告只把v2當獨立驗算證據。
- 先前把tumor-REF長尾記成HCC1954 chr14屬runtime診斷誤判。依canonical suffix與實際矩陣重查，主要極端
  位點位於chr17，包括`chr17:39520424 G>T`（joint 4,549 reads）與`chr17:39560316 C>G`（joint 3,716）。
  active serial worker仍持續前進，但ordered writer只持久化到102,000/102,842。
- 採fresh recovery output：只接受open gzip中完整且依canonical task order的prefix；缺失suffix以pinned
  seed-parallel helper重算；最後重新建立dedup annotations、summary、manifest與source attestation。舊serial
  output不覆寫、不刪除；recovery通過後才停止並封存其未完成chain。

## 2026-07-17 cooccurrence Mapping incident 與 release identity 補強

- M2v5共現全量第一次啟動以12 workers、chunk-size 8、999 permutations執行，在worker第一次呼叫
  `m2_categorical_level_counts()`時以`NameError: name 'Mapping' is not defined`結束。要求的輸出目錄
  `methyl_ssnv_cooccurrence_v2_m2v5_source_locked`並未建立，沒有formal receipt，也沒有修改任何上游輸入；
  因此此執行標為`ABORTED_EXCLUDED`，不可作科學結果。
- 完整incident、輸入、命令、exit code與修正證據位於
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_m2v5_mapping_nameerror_incident.v1.json`。
  修正只新增`typing.Mapping`匯入，並增加直接runtime regression；定向測試`115 passed`，完整topic測試
  `277 passed, 2 known deprecation warnings`。
- fresh rerun固定輸出到
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v2_m2v5_mapping_fix_source_locked`，
  禁止覆寫或復用前次要求路徑。此第二次啟動進一步暴露`cluster_sizes`雖存在完整screen TSV，卻未列入
  candidate-screen required fields，也未投影到worker task；因此以`KeyError: cluster_sizes`結束，仍未建立
  output或receipt。修正後缺欄會在worker啟動前fail，完整欄位會進入task；定向`29 passed`、完整topic
  `281 passed, 2 known deprecation warnings`。第三次正式路徑改為
  `methyl_ssnv_cooccurrence_v3_m2v5_task_contract_fix_source_locked`，啟動前另跑全部102,842 tasks真資料preflight。
- `independent_m2_gate_recount.v2.json`建立後，production gate source於10:29另有程式身份變更，導致v2記錄
  SHA-256=`05ad...`與目前`m2_screen_gate.py` SHA-256=`f3d...`不一致。數值不因此自動失效，但v2不得作
  current-release source identity；正式收尾改用fresh `independent_m2_gate_recount.v3.json`全量重算並綁定
  現行gate source。
- 新增`run_m2v5_recovered_completion_chain.sh`，只在cooccurrence、recovered tumor-REF、source snapshot與
  independent M2 v3全部PASS時才執行strict、matched-normal、CN/CCF、post audits、final dataset、v5報告、
  portable HTML與desktop/mobile零overlap QA。腳本`bash -n`與10項runner contract tests均PASS。

## 2026-07-17 全task預檢、M2 release接縫修正與正式v3共現啟動

- 正式真資料preflight命令逐筆重建102,842個reduced worker tasks，輸入為frozen manifest、v10 site/assignment、
  v2 extraction root與`independent_m2_gate_recount.v3.json`；輸出
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_task_contract_preflight.v1.json`。
  實際exit=0、`pass=true`：M2 eligible=919、evaluable=1,867、evaluable-ineligible=948、axis-indeterminate=100,974、
  K>10=1均exact；raw assignment constant counts為HP-exact/family/strand=93,534/93,534/39，gate-evaluated為
  93,533/93,533/39。唯一差值由`HCC1954 chr5:751076 G>A`的K=11先行group-limit gate完整解釋。
- Preflight綁定analyzer SHA-256=`29d1754d07a5d71523c4ab40fc4fcc8bfb1fb731d243f7407d9caddecab9c46c`、
  M2 gate SHA-256=`f3d932388bf8ddc90f8fa1a136f6203dfe8d9f5aa9979ae32e9d936c2a413b79`；analyzer、gate、
  `ssnv_cooccurrence_lib.py`、`latest_tag_join.py`、preflight script與receipt均設為mode 0444。preflight只證明
  execution/task contract，不是科學結果。
- 獨立M2 method agent指出final dataset已產生6項`downstream_checks`，但report validator仍只接受舊4項；
  同時final denominator/report methods仍殘留「所有低power軸皆NOT_EVALUABLE」的對稱舊語意。兩項High blocker
  已修正：report精確要求6項checks，final/report均採非對稱denominator。另結構化保存categorical planning
  ceilings HP-exact/family/strand=7/5/2，並明示assignment observed levels只作constant-axis proof、不取代planning model。
- 新隔離回歸以兩群各70 reads驗證HP-exact最低80% power n=152；n=140時若effect/p已顯著對齊、其餘七軸可判，
  正確結果仍為`evaluable=true / INELIGIBLE_SCREEN_HP_CONFOUND`。少任一非對稱downstream check或planning ceiling
  漂移均在建立report output前fail closed。最新完整相關測試為`303 passed, 2 known deprecation warnings in 17.07s`。
- 第三次正式cooccurrence於source freeze與preflight PASS後啟動，session=`50515`；命令固定28 workers、chunk-size=8、
  max-pending=56、999 permutations、top-markers=3、exact-state-space-ceiling=250,000，輸出fresh路徑
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v3_m2v5_task_contract_fix_source_locked`。
  啟動時tumor-REF仍有11個高深度inner workers，故未直接使用40個共現workers；待前者自然完成後CPU會釋放給現行
  process pool，不重啟或覆寫任何producer。

## 2026-07-17 SciPy runtime incident、v2 preflight與fresh v4共現

- 第三次共現session `50515`在第一個worker batch的3x2 categorical association立即以
  `AttributeError: 'tuple' object has no attribute 'pvalue'` fail-closed。固定cnvtools環境為Python 3.11 / SciPy 1.9.2，
  `chi2_contingency`實際回傳四元素tuple；analyzer錯用較新版named result的`.pvalue`。要求的v3 output path不存在，
  formal receipt未建立，上游輸入未修改；此執行加入原incident receipt並標`ABORTED_EXCLUDED`。
- 修正只把Pearson p-value讀取改為跨版本穩定的tuple index `[1]`；2x2 Fisher原本即使用index 1。新增實際runtime
  3x2 Pearson與2x2 Fisher probes，並掃描本專題其他SciPy named-result依賴。新analyzer SHA-256=
  `31414080f899ec36cfbdfe6fd407def1a00c7e3bbb7ea0c7d1564d9b68ad2731`。
- 不沿用舊preflight。fresh `cooccurrence_task_contract_preflight.v2_runtime_compat.json`重新掃描全部102,842 tasks，
  同時要求Pearson/Fisher五個runtime checks；實際exit=0、schema 1.1.0、所有原task/M2守恆與新runtime checks皆true。
  Preflight script SHA-256=`0ef01ece371533e738b277be096b4ab82a36a5e79a61f280076712f70fedcfa4`，receipt綁定current
  analyzer/gate後設為mode 0444。修正後完整相關測試為`305 passed, 6 known deprecation warnings in 16.95s`。
- completion runner改為只接受v4 cooccurrence路徑，並新增v2 runtime preflight `pass=true`硬條件；`bash -n`與
  10項runner tests全PASS。fresh v4正式共現session=`20378`，28 workers、chunk-size=8、max-pending=56、
  999 permutations、top-markers=3、exact-state-space-ceiling=250,000；輸出
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v4_m2v5_scipy_runtime_fix_source_locked`。

## 2026-07-17 canonical library runtime closure、preflight v3與fresh v5共現

- 第四次共現session `20378`在第一批worker進入repository-canonical
  `InterSubMod/scripts/ssnv_cooccurrence_lib.py`的3x2 Pearson路徑時，以相同SciPy 1.9.2
  `AttributeError: 'tuple' object has no attribute 'pvalue'` fail-closed。v4輸出目錄不存在、formal receipt未建立、
  上游輸入未修改；已納入incident並標`ABORTED_EXCLUDED`。根因是v2 preflight只執行analyzer probe，未綁定或
  執行另一個實際worker會載入的canonical library。
- canonical library改用`chi2_contingency(...)[1]`，preflight schema升至1.2.0並同時執行analyzer 3x2、
  canonical-library 3x2與Fisher 2x2真實runtime路徑。完整相關測試實際`306 passed, 10 known deprecation warnings`
 ；warning均為SciPy呼叫NumPy deprecated API。
- fresh formal preflight輸出為
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_task_contract_preflight.v3_full_runtime.json`，
  exit=0、`pass=true`。102,842 tasks/assignments/sites exact；M2 eligible=919、evaluable=1,867、
  evaluable-ineligible=948、axis-indeterminate=100,974、K>10=1。runtime probes全部true。
- source readback與receipt一致：analyzer=`31414080f899ec36cfbdfe6fd407def1a00c7e3bbb7ea0c7d1564d9b68ad2731`、
  preflight=`4260cc48d8420412496f992cbabf709dc7c6cfc66737c0414076bd1f613cc84c`、canonical library=
  `3d41023b03f1f047312011f8513b685bf7ffe25f5d0795d253bc9013b9ec88ff`、M2 gate=
  `f3d932388bf8ddc90f8fa1a136f6203dfe8d9f5aa9979ae32e9d936c2a413b79`、latest-tag join=
  `7ddadb05481bb61f2dc3a489762944ec347ab9f5265c69fb17538efd44b543cd`；相關source與receipt均mode 0444。
- fresh v5於session `96734`啟動，參數28 workers、chunk-size=8、max-pending=56、999 permutations、
  top-markers=3、exact-state-space-ceiling=250,000。輸出fresh路徑為
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v5_m2v5_full_runtime_fix_source_locked`；
  producer完成並產生`pass=true` receipt前不得引用為科學結果。
- 生物學獨立reviewer另指出G1 global family、G2=0 fallback、C1=0、ancestral ALT與background superset的
  5項報告語意風險。報告層已改為：所有M2 exact-testable pairs先進global BY；G1=0不虛構局部共分離；
  C1=0標為缺joint model的結構性NOT_EVALUABLE；ancestral ALT只稱模型預測相容；superset只限同一background
  payload predicates。定向report tests=`47 passed`，未修改凍結producer或claim-contract-v5。

## 2026-07-17 recovered tumor-REF全量完成

- recovery producer exit=0、`pass=true`，正式輸出102,842/102,842個M1位點；保存open-gzip中102,307個
  canonical完整prefix rows，僅重算缺失suffix 535 rows。prefix/suffix key sets disjoint、canonical task order exact、
  source screen fields與assignment digest exact；原live prefix檔案before/after SHA-256同為
  `7debf12d5db46439054b4bbfda5513d33061e4f720591161504e97634ab3e5ae`，沒有修改或覆寫。
- 正式site TSV SHA-256=`c85e8c41753a2563765287d4b2e1f5cf335cf55140ad16a6ea883cd7ba2a9312`；
  `gzip -t` exit=0，行數102,843=header+102,842 data。run manifest鎖定recovery analyzer、serial analyzer、
  focal library與pinned seed-parallel helper，且serial/parallel exact-equivalence receipt為PASS。
- site-weighted tumor-REF描述性結果：REF evaluable=83,688/102,842，lenient stable multigroup=
  33,206/83,688=39.6783%；joint evaluable=85,777，joint stable=53,191/85,777=62.0108%。這是背景
  partition的描述與B1保守veto來源，不是subclone比例；同一background payload的lenient predicate缺membership
  ARI，不能與primary M1 ALT實際flag集合直接互稱superset。
- report wording修正後的完整相關suite重跑為`307 passed, 10 known deprecation warnings`；completion runner
  `bash -n` exit=0。fresh cooccurrence v5仍active，G1/G2/R1/B1/C1仍不得發布。

## 2026-07-17 source-attestation relative command blocker與schema 1.2 closure

- recovery完成後先以正式verifier邏輯做`/tmp` pre-probe，舊v1 verifier exit=1：snapshot live cmdline在移除
  Python launcher後與manifest逐token完全相同，但manifest的script token是repo-relative，v1另要求該raw token
  等於absolute analyzer identity，故必然fail。正式receipt尚未建立，下游未啟動，沒有混合release輸出。
- v1 snapshot creator檔不可修改，因snapshot已鎖定其SHA-256=
  `7927a86de1082860959b5a51ecb74ee8d4be49c55974a4edb783e757f2acb9f0`。新增post-run verifier v2，
  receipt schema 1.2.0分開保存：執行中snapshot creator before/after identity，以及完成後verifier identity。v2只允許
  absolute exact path或無空段/`.`/`..`的relative suffix綁定；Python launcher後所有其他tokens仍須exact。
- v2 verifier SHA-256=`f02461f61fdda41333f065750c13d45dcebf06b304a5325a6df21912129436b9`、mode 0444。
  真實snapshot/manifest pre-probe exit=0、`pass=true`，command mode=`relative_suffix`，全部16個checks通過。
  `other.py`、`../producer.py`、`./producer.py`、manifest extra token與producer source drift均有fail-closed regression。
- final dataset loader只接受schema 1.2.0與新增verifier/token checks，並在建立dataset時重驗post-run verifier當前
  path/size/SHA。current `run_m2v5_recovered_completion_chain.sh`改用v2；已mode 0444的歷史
  `run_source_attested_release_chain.sh`不修改且不得用於本次release。完整相關suite=`313 passed, 10 known warnings`，
  定向source/final-builder/runner=`69 passed`，Python compile與current runner `bash -n`均exit 0。

## 2026-07-17 source-attestation獨立review與strict repo-relative closure

- 有效唯讀reviewer `019f6eac-6ea8-7fd1-8271-87beb582dd72`指出：v2 verifier雖已驗完整
  `command_binding`，final builder仍只接受receipt的boolean checks與宣告verifier hash；另外relative token採任意suffix，
  basename-only在一般化情境仍可能誤配。原判定無High、2 Medium、2 Low；正式receipt尚未建立屬預期pending，
  不能由pre-probe代替。
- verifier現只接受absolute exact或`source.resolve().relative_to(repo_root)`的完整精確token；relative mode固定為
  `repo_relative_exact`，basename-only、`scripts/<name>`、`.`與`..`全部拒絕。current verifier SHA-256=
  `dc785d129b3eb3c34bf8bce12907d57f7069be5cf71fccbf6ae9b0307c9733ea`、mode=`0444`。
- final builder現固定可信v2 verifier的current path/size/SHA，讀回hash-bound snapshot後獨立重驗：snapshot與receipt
  source/creator identities、`live_cmdline[1:] == manifest.command`、完整analyzer token、expected command fragment與
  全部command-binding fields。錯verifier或自行改寫binding的替代receipt會在任何output建立前被拒絕。
- 定向verifier/final-builder測試實際`51 passed in 7.41s`。真實snapshot+recovered manifest fresh pre-probe輸出為
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json`，
  `pass=true`、schema1.2.0、mode=`repo_relative_exact`、verifier mode=`0444`。正式
  `post_run_source_identity.receipt.json`仍保持不存在，只能由current completion runner建立。
- 完整review與逐項closure：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/reviews/20260717_tumor_REF來源attestation與release接縫獨立稽核_01.md`。
- 原reviewer完成bounded re-review：High/Medium/Low=`0/0/0`、verdict=`GO`。其唯讀runtime replay確認
  complete repo-relative/absolute exact通過，basename/短suffix/dot/parent拒絕；fresh pre-probe經builder為
  `release_gate_pass=true`；wrong verifier、receipt binding drift、manifest command drift與snapshot SHA drift皆被拒絕。
  修補後完整相關suite為`314 passed, 10 known deprecation warnings in 19.80s`，`py_compile`與runner `bash -n`
  exit 0。此GO只封閉source-attestation接縫，不代表active cooccurrence或final scientific release完成。

## 2026-07-17 cooccurrence v5 raw-BAM identity incident、嚴格RG-only契約與全量preflight v4

- 第五次cooccurrence session `96734`以28 workers執行約2小時23分後fail-closed；錯誤位點為
  `HCC1395 chr1:13200585 A>C`（UNASSESSED，不是G>A），6個預期focal-ALT projections中有2個各對到2筆
  C++-eligible raw BAM records。要求的v5路徑
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v5_m2v5_full_runtime_fix_source_locked`
  不存在，formal receipt=0、上游修改=0；attempts 1-5全部`ABORTED_EXCLUDED`。完整incident v2為
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_execution_incident.v2_through_attempt5.json`。
- 兩個QNAME為`45a0a6af-df78-4fb4-8d94-689b30bd7efa`與
  `d9fcc979-8e0e-483a-863b-888763e58a15`。逐record比較確認FLAG、reference start/end、MAPQ、CIGAR、
  mate fields、template length、sequence、qualities、MM、ML與其他auxiliary tags完全相同；唯一差異是`RG`由
  原ID變為同ID加`-683557E`。LongPhase-S sidecar對兩個QNAME也各含2筆完全相同的full-identity/HP/PS rows，
  原sidecar契約以set折疊後distinct full identity count=1；因此是同一分析payload的read-group duplicate，
  不是兩條不同alignment或不同molecule的證據。
- 新raw duplicate契約固定為`sam_core_and_all_aux_tags_except_RG_exact_v1`：只有SAM全部core fields與除`RG`
  外所有typed auxiliary tags exact時，才能以BAM fetch順序選第一筆等價代表，並輸出折疊audit。sequence、quality、
  CIGAR、FLAG、MM、ML、HP、AS等8類定向衝突測試均hard fail；exact duplicate與RG-only duplicate各有正向測試。
  真實位點完整重播為expected=6、matched records=8、RG-only duplicate projections=2、extra records=2、
  focal ALT readset digest exact。
- v5 analyzer與v3 preflight source已另存唯讀快照：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/aborted_sources/analyze_methyl_ssnv_cooccurrence.v5.sha31414080.py`
  SHA-256=`31414080f899ec36cfbdfe6fd407def1a00c7e3bbb7ea0c7d1564d9b68ad2731`；current v6 analyzer
  SHA-256=`2b0df7e68e583ff88e409bb1d117455d66595399c469db69585c461c23a4dd2d`、preflight v4 script
  SHA-256=`cc0edd01f3148605ef03fc852f05a7f97d691d64bb628741b80adae0bad85a97`，兩者均mode 0444。
- v3 preflight的缺口是只建立102,842 tasks與執行統計API probes，未對每個task執行raw BAM recovery。v4新增
  全task `recover_site_reads`重播、focal call與ALT readset digest重驗、逐樣本/site-weighted折疊計數與digest、
  missing/conflicting=0 gate、v5 incident位點必須exact replay，以及preflight期間source before/after identity不變。
  輸入為frozen manifest、v10 assignments/sites、primary artifact audit、independent M2 v3與v2 extraction root；命令
  使用40 workers、chunk-size=8、max-pending=80、BLAS/OpenMP各1 thread。active session=`58967`，預期輸出
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_task_contract_preflight.v4_raw_identity_full_runtime.json`；只有`pass=true`才允許執行
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_cooccurrence_v6_source_locked.sh`。
- fresh v6正式輸出固定為
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v6_m2v5_raw_identity_contract_source_locked`。
  launch與completion runners會比較preflight/analyzer SHA、raw identity counts與site-weighted digest；任何差異阻擋下游。
  最新完整相關測試為`325 passed, 10 known deprecation warnings in 16.32s`，新增runner定向為
  `52 passed, 4 known deprecation warnings`，`py_compile`、`bash -n`、targeted `git diff --check`均exit 0。

## 2026-07-17 partial preflight v4排除與dependency-attested v5修正

- `[偏離]` session `58967`在24,000/102,842 tasks後主動停止，exit 130；原因不是資料位點失敗，而是獨立review發現release契約仍會遺失SAM `B` array subtype、允許site audit欄位互相矛盾，且以硬編碼`missing/conflicting=0`造成錯誤保證。v4沒有JSON receipt，partial log SHA-256=`e1a2b760c313c2166dff0569cb255fc20ac215e12613ac091e0756f6992090d3`，狀態固定為`PRECHECK_ABORTED_EXCLUDED`。
- `[決策]` v4a analyzer/preflight/launcher已保留mode 0444快照於`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/aborted_sources/`，SHA分別為`2b0df7e68e583ff88e409bb1d117455d66595399c469db69585c461c23a4dd2d`、`cc0edd01f3148605ef03fc852f05a7f97d691d64bb628741b80adae0bad85a97`、`effdecac19bd2d408ef1f16651323940c66a4b2a47caeafb89cc6b29bb8ee5bd`；不可刪除或替代。
- `[決策]` current analyzer改以`(SAM type, array.typecode, length, raw bytes)`比較typed array；`PS/HP/MM/ML`及所有非RG aux差異均hard fail。新增逐位點invariants：`matched=expected+extras`、`exact+RG-only=duplicate projections`、`duplicates<=expected`、`extras>=duplicates`、digest/empty semantics與完整sparse duplicate record reconciliation。
- `[決策]` 不再宣稱未實際計數的`missing/conflicting=0`；receipt只記錄`hard_fail_before_site_result` policy、所有已產出site rows通過invariant validation、failure counts未物化。新增`raw_identity_duplicate_audit.tsv.gz`保留所有site-weighted duplicate projections、record count、classification與差異tags。
- `[決策]` preflight source lock擴為preflight、analyzer、`ssnv_cooccurrence_lib.py`、`latest_tag_join.py`、`m2_screen_gate.py`五檔，before/after SHA/size/path/mode 0444全部一致；launcher與completion runner逐檔重驗。唯一下一個receipt為`cooccurrence_task_contract_preflight.v5_dependency_attested_raw_identity_full_runtime.json`，不得沿用v4名稱。
- `[驗證]` 聚焦回歸`53 passed`；完整相關suite`336 passed, 10 known SciPy/NumPy deprecation warnings in 15.91s`；`py_compile`、`bash -n`與targeted `git diff --check`exit 0。`ruff`未安裝，未將此工具列為通過證據。current frozen SHA-256：analyzer=`bdf6644293c9744cf43a8fa925039b9e4e261a7cf0499034dc5814cbcc2c2908`；preflight=`fabbfd3e993e757cc1295266026807bab06d4915159a65053a17bb22c630ac1f`。

## 2026-07-18 attempt 6 summary-totality incident與v5 source freeze

- `[偏離]` 正式cooccurrence attempt 6於`2026-07-18 10:30:42+08:00`啟動，40 workers、
  chunk-size 8、999 permutations；約3小時54分後所有worker futures已進入`build_summary`，
  但舊summary直接讀取每一列`m2_axis_statuses.hp_exact.status`，對唯一K=11列拋出
  `KeyError: hp_exact`並exit 1。要求的v6 output directory不存在、run/release receipt均未建立；
  只有5403-byte log保留，SHA-256=`7858039ffb1c41b70b81849f85cbf070d9fe3e6d504a21a5c8981547f07cc0ad`。
- `[根因]` `HCC1954 chr5:751076 G>A`有11個M1 methyl groups；M2 planning model固定只接受K<=10，
  因此gate正確在axis evaluation前回傳`NOT_EVALUABLE_M2_GROUP_COUNT_GT10`與空axis object。
  錯誤只在summary totality，不是M2判定、read identity、worker計算或資料缺失。全102,842個M1列的K分布為
  2:`76,402`、3:`19,003`、4:`5,430`、5:`1,471`、6:`398`、7:`93`、8:`33`、9:`10`、
  10:`1`、11:`1`。
- `[決策]` 新summary對八個固定axis逐一保存完整102,842-row census。pre-axis停止列只寫
  `NOT_EVALUATED_M2_GATE_PRE_AXIS:<原gate status>`；不補effect、不補p-value、不改寫成aligned、
  unaligned或negative evidence。空axis但不是明確`NOT_EVALUABLE_M2_*`、partial/extra axis、
  非object evidence、缺status或任一axis不守恆都hard fail。
- `[驗證]` K=11正向與partial-axis負向回歸已加入；targeted=`96 passed`。固定隔離環境的完整canonical
  suite=`448 passed, 0 failed/errors/skipped`，JUnit SHA-256=
  `957af7558bb1c0d2d1d2daaae03221449abb22da2e5e98a0db4a201ba39349f1`。
  23個protected source全部mode 0444，source-set SHA-256=
  `e42db240ddfffe560f9c1e22eede01294f6c1dc7d4ea9ede79ab8540cd3c4066`。
- `[未決]` attempt 6不得作科學輸出或部分完成資料。下一次只允許fresh v5 source authority、
  `stable_primary_artifact_audit.v4_source_authorized_pre_downstream.json`、
  `cooccurrence_task_contract_preflight.v8_summary_hotfix_full_runtime.json`與新
  `methyl_ssnv_cooccurrence_v7_m2v5_raw_identity_contract_source_locked_summary_hotfix`路徑；
  v6 log不得覆寫。完整incident為
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_execution_incident.v3_attempt6_summary_totality.json`。
- `[review closure]` 獨立agent指出空axis map若只用`NOT_EVALUABLE_M2_*`前綴，理論上也會接受
  post-axis `NOT_EVALUABLE_M2_AXIS_INDETERMINATE`。雖然current producer不可達，source freeze前仍改為
  三個真實pre-axis status明確allowlist，並新增post-axis空map必須fail的負向測試；runner成功訊息同步改為
  `Cooccurrence v7 summary-hotfix PASS`。targeted=`71 passed`、canonical=`449 passed`，
  JUnit SHA-256=`2a6489fb4c2386331e7e3db0684c524b148c9e4d2cd9b40f83822b7041a79551`；
  final protected source-set SHA-256=`088ccd47ecee4462876c0ce4bb7a4f054f50d28dc52ce995434a0b0cd9221295`。
  先前v8 reviews只適用舊digest並標記superseded，不能用於v5 authority；fresh v9 review為下一個gate。
- `[review closure]` fresh v9 A/B均判定pre-axis allowlist與source-authority機制`APPROVE`，但B發現completion
  shell已宣告v5 final dataset/report，Python builder/finalizer仍硬鎖v4 canonical output，正式鏈會fail-closed而無法
  完成；A另指出extra-axis雖由與partial-axis相同guard拒絕，沒有獨立負向測試。source authority尚未建立或簽署，
  因此先修補而不消耗one-time key。兩個Python canonical路徑已對齊v5，新增跨producer路徑一致性與extra-axis
  專屬回歸；targeted新增`2 passed`，canonical=`451 passed, 0 failed/errors/skipped`，JUnit
  SHA-256=`a2d4e94ac9cf2cc411602d28ae4c8734904351340b923ff63a29bf2e1c308dab`。23個protected source仍全為
  mode 0444，新source-set SHA-256=`4685ff85bbad8c3174dcfcd3b60dcd7d20a940a6caf11d43558ca079e22029b1`。
  v9 reviews只適用被取代的`088ccd47...` source set；fresh v10 review為簽署前最後gate。
- `[review closure]` fresh v10 A/B再次獨立指出：v5 output雖已對齊，但builder的7個canonical input仍指向
  舊v2/v6，會在production `validate_canonical_task_b_paths` fail-closed；B判定`REQUEST_CHANGES`。另外
  matched-normal正式runner在有候選時產生analysis dir，零候選時只產生controls dir的明示NOT_APPLICABLE receipt，
  因此不能只硬鎖單一路徑。authority仍未建立或簽署。
- `[決策]` cooccurrence、preflight、strict、primary pre/post、matched-normal與CN/CCF canonical inputs已完整
  對齊current v7/v8/v3/v4；builder/finalizer只接受兩條matched-normal canonical commands：
  analysis-v3（非零候選）或controls-v3（零候選明示NOT_APPLICABLE）。定向tests=`2 passed`；獨立shell parser
  對11個固定inputs、2個matched branches、2條builder commands及v5 outputs全部exact PASS。
- `[驗證]` canonical=`452 passed, 0 failed/errors/skipped`，JUnit SHA-256=
  `511366ce5fe27a8382820b99fd1ddbf6d199c32f1ec58f1865e637ab0bd427fd`；23個protected sources仍全mode
  0444，新source-set SHA-256=`20232fc5433fdc5a5210c731bfd3b2cbf8ed0c1c04abc2781196ea67b08aad12`。
  v10 reviews只適用被取代的`4685ff85...`；fresh v11全鏈path review為簽署前gate。
- `[review closure]` fresh v11 A錯把builder層對齊當成全鏈完成；B正確發現matched-normal runner/analyzer
  自身仍硬鎖cooccurrence-v6與matched-v2，因此formal authority下會在`resolve_release_command`直接fail。
  全23-source stale-token掃描再多發現cooccurrence release finalizer仍鎖primary-pre-v2/preflight-v6/
  producer-v6，CN/CCF annotator仍鎖cooccurrence-v6/output-v2。authority仍未建立或簽署。
- `[決策]` 四個leaf producers均已對齊：cooccurrence finalizer=`primary-pre v4 + preflight v8 + producer v7`；
  matched runner/analyzer=`cooccurrence v7 + controls/analysis v3`；CN/CCF=`cooccurrence v7 + output v3`。
  新增5個回歸：三個producer常數、completion shell option/value跨source、以及遍歷
  `EXPECTED_SOURCE_PATHS`全部23角色的superseded exact-token zero gate。定向tests與跨模組映射全部PASS。
- `[驗證]` canonical=`457 passed, 0 failed/errors/skipped`，JUnit SHA-256=
  `b671f9b92d081a318e5ffe60beb60ee60ced010ff27d8424ab9ee9ab2c3ff94e`；23個protected source全mode
  0444，新source-set SHA-256=`db0c33f89c31e6b1a3e4cd68a2943df3c109f0e76ddaff40e9e0d2bf78faa3ed`。
  v11 reviews只適用被取代的`20232fc5...`；fresh v12 all-producer review為簽署前gate。

## 2026-07-18 v12 外部審查 session-limit 事件與 supplemental preformal 驗證

- `[偏離]` fresh v12 reviewer A/B 於`2026-07-18 16:36:59+08:00`同時啟動，但Claude Code在讀取
  review request後立即以session limit結束，exit=1，訊息為
  `You've hit your session limit · resets 7:40pm (Asia/Taipei)`。兩次均未讀取producer、未回傳verdict，
  不可作release approval。
- `[證據]` 兩個61-byte raw outputs已各自保留為
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/20260718_external_claude_reviewer_{A,B}_v12_all_producers_attempt1_session_limit_20260718T1637+0800.raw.txt`，
  mode=`0444`、SHA-256皆為
  `75b9b20a8cfcf38bc07923da604ce44e467079e153816aa4c41565926afc3e68`。對應結構化failure records明示
  `FAILED_NO_VERDICT`與`verdict_usable_for_release=false`。
- `[未決]` 兩個fresh v12唯讀review已排程於`2026-07-18 19:42:00+08:00`重新執行，fresh success raw/JSON
  必須是兩個不同reviewer UUID、`APPROVE`、blocking findings空、source digest=`db0c33...a3ed`，才能組裝及
  一次性簽署v5 source authority；失敗attempt不得被assembler讀取。
- `[驗證]` 等待期間獨立重算23個protected roles仍全mode=`0444`且source-set digest精確為
  `db0c33f89c31e6b1a3e4cd68a2943df3c109f0e76ddaff40e9e0d2bf78faa3ed`；canonical XML仍為
  457 tests、0 failures/errors/skips，SHA-256=`b671f9...f94e`。
- `[決策]` 既有positional audit由v4 authority簽章，supplemental report另要求未來v5 final
  dataset/report雙簽章；這是刻意的雙鏈來源，不得把v4 audit chain直接改寫成v5。
- `[決策]` M2八軸包含HP exact/family，但不包含PS。從30個M2 PASS frozen assignments重算：
  3,505條core reads中131條missing PS，10/30 sites至少一條missing PS，僅1/30跨多個非空PS。
  `HCC1937 chr5:43849776 T>C`為主要caveat：88/109 missing PS，其餘21條跨PS
  43,668,888與43,913,176；methyl group `1-2`的5/5 reads皆missing PS。因此報告不得宣稱phase-set
  confounding已排除。
- `[決策]` 真實read-by-CpG圖由每個positive dataset最大core-read案例，加上固定
  `HCC1395_DORADO chr1:20467811 G>C`、`chr22:47518662 A>G`與
  `HCC1937 chr5:43849776 T>C`，去重後共8個案例；每panel顯示PS missing與非空PS數。
- `[歷史驗證]` 新PS census、必要個案與fail-closed tests初次加入後，positional supplemental audit、report
  builder及release finalizer三組聚焦測試為`72 passed in 3.43s`。該次JUnit位於
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_positional_singleton_delivery_ps_caveat_v2.xml`，
  mode=`0444`、SHA-256=`0e392ca2c398596d1756953ed0642a91cb5114eeab16a315ba2911afa7726c59`。
  pre-PS v1與72-test v2 JUnit均保留為歷史preformal證據，不作current report-builder驗證。
- `[決策]` supplemental report builder不再只淺層讀取final dataset/report receipts。它現在要求canonical
  v5 receipt/signature/public-key路徑，直接呼叫final release validator重驗兩條Ed25519簽章、source
  authority、dataset/report互鏈、全部正式輸出identity與一次性private key mode=`000`退休狀態。
- `[決策]` supplemental report在staging Markdown、HTML、receipt與`_SUCCESS`完成後、atomic
  no-replace publish之前，再次重驗全部輸入、repo/source identity、current JUnit及正式雙簽章鏈，並與
  起始snapshot逐欄比較；任何晚期輸入或signature drift均fail closed，不發布部分結果。
- `[歷史驗證]` v16完整canonical suite納入既有producer tests、PS census與案例選取、正式release-chain
  validators及pre-publish drift guards，共`471 passed, 0 failures/errors/skipped`。JUnit位於
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_supplemental_v16_ps_release_publish_guard.xml`，
  mode=`0444`、SHA-256=`813217316082cac97114aa0b33e81f4b22ee850ffb93d75b54e4f5c106d46283`。
  執行後23個protected source仍全mode=`0444`，source-set digest仍精確為
  `db0c33f89c31e6b1a3e4cd68a2943df3c109f0e76ddaff40e9e0d2bf78faa3ed`。
- `[質疑與修補]` internal reviewer指出v16仍可由任意471-test同名synthetic XML冒充，且late guard漏掉
  v4 positional authority manifest/approval/signature/public key與supplemental auditor。F4正式
  dataset/report fd-bound雙簽章重驗已判CLOSED；F3/F5因此再修補，未改科學門檻、分母或claim ceiling。
- `[決策]` current supplemental test evidence改為一次性Ed25519簽署：formal builder只接受固定v17 XML
  路徑，signed manifest同時綁定XML、report builder source與test source的path/size/SHA-256/mode，
  並要求private key已退休mode=`000`。Manifest與signature分別為
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/supplemental_report_test_evidence.v1.json`
  與同路徑`.ed25519.sig`；signature SHA-256=
  `3aeba864e0ff22c271c258d89fca8b93aae4c97718f2cfd98ce7cdfdac39af8d`。
- `[決策]` v4 transitive authority四檔與supplemental auditor現納入起始`input_paths` identity snapshot、
  report receipt與pre-publish snapshot；late guard另再次執行完整`validate_source_authority_chain()`。
  任一late mutation會在atomic rename前fail closed。
- `[驗證]` current v17 suite為`476 passed, 0 failures/errors/skipped`、476 testcase、17/17 required tests；
  XML mode=`0444`、SHA-256=
  `203d923faa17ddb6309dbb28a4985dda52e5739646f9264eba8fe22a6b1d150a`。Frozen builder/test source
  SHA-256分別為`7e75997e6ca50a38a01465516e6868403e17a33e31a712f79d623b27122a8a23`與
  `032d445375914a90ffb1aa3fdcdd5c17f4065e5273f3a2c3d12b0cd00993126e`，皆mode=`0444`。
  真實FD-bound validator回傳`signed_evidence_verified=true`；23個protected source digest仍為
  `db0c33f89c31e6b1a3e4cd68a2943df3c109f0e76ddaff40e9e0d2bf78faa3ed`。
- `[決策]` fresh v12 Claude raw output不以人工或寬鬆`json.loads`直接升格。新增
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/normalize_external_source_review_v12.py`，
  逐層拒絕duplicate JSON keys、trailing prose/code fence、缺漏或額外top-level keys、非canonical UUIDv4、
  A/B reviewer identity錯置、錯誤HEAD/source digest及與findings矛盾的`APPROVE`；normalized output以
  `O_EXCL`完整write loop建立並固定mode=`0444`。
- `[驗證]` normalizer targeted tests用
  `/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -m pytest -q
  research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_normalize_external_source_review_v12.py`
  執行，實際為`11 passed in 0.08s`。第一次誤用不存在的
  `/bip7_disk/liaoyoyo2001/InterSubMod_envs/cnvtools/bin/python`僅得到exit=`127`、未執行tests，已明確
  排除不作驗證證據。Frozen normalizer/test SHA-256分別為
  `46e68a2964d3b3b6488458b12f6e2004b39a0d9c7ae6a88c8d523428323282b7`與
  `2cbf6ea10931e3fb955d0b32ada5b40963e733cb981949310e0fcd50786777aa`，皆mode=`0444`；獨立重算
  23個protected roles仍全mode=`0444`，source-set digest仍精確為
  `db0c33f89c31e6b1a3e4cd68a2943df3c109f0e76ddaff40e9e0d2bf78faa3ed`。
- `[質疑與修補]` internal reviewer第三輪前複核判定F3 CLOSED，但F5仍有post-validation publish
  window：`validate_pre_publish_state()`返回後仍執行staging directory fsync與chmod，才做no-replace
  rename。這兩個staging mutation已全部移到final guard之前；新`validate_and_publish()`函式體只允許
  `validate_pre_publish_state(...)`後立即`rename_no_replace(staging, output_dir)`，guard不再先返回main。
- `[驗證]` 新
  `test_validate_and_publish_contains_only_guard_then_atomic_rename`同時以AST要求函式內expression calls恰為
  `validate_pre_publish_state, rename_no_replace`，並以monkeypatch事件序列重驗`guard, rename`；main ordering
  test另要求`success write < staging fsync < staging chmod < validate_and_publish`且main不得直接呼叫rename。
  專項為`64 passed in 3.39s`。一次版本矩陣尚留476時曾得到`1 failed, 63 passed`，已只將canonical
  contract參數同步為488/487後全綠，不隱藏該中間失敗。
- `[驗證]` v18 canonical scope為root兩組共現tests加topic全tests；collect精確為488。dry run為
  `488 passed, 38 warnings in 39.32s`，final JUnit run為
  `488 passed, 38 warnings in 38.61s`；warnings全為既知SciPy/NumPy deprecation。XML位於
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_supplemental_v18_atomic_guard_publish.xml`，
  mode=`0444`、SHA-256=
  `b7396deb743cc2bc50a2419ecf19556396aab8041df8f8244add4e72e2d09268`；builder/test mode皆`0444`、
  SHA-256分別為`4c02533e9ca9783aba6dada9cb4f36a1d3fe58e8dbe8a749869cca7c59570bf2`與
  `776eb8ecabdb1ce1b0b9a09aa6934e5dafd89ca970624f40efaa6e4b6efab295`。
- `[決策]` v17/v1 signed evidence保留為歷史pre-fix證據，不覆寫。Current v18改用v2 manifest
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/supplemental_report_test_evidence.v2.json`
  （SHA-256=`ac8d7a0f4ea7a262403323884cec1b985c1bdebad1fc4d6a5aaf829ab3d923cb`）及detached
  signature（SHA-256=`5a7ee1a73afcbb1dcd8e3aa14706542cefbd6aa09d9532df72112d2e70f03ae7`），皆mode=`0444`。
  新public key SHA-256=`a7d7933deff31339f7f7af92ac1cf4db280c12b49b1c36b5db1a955e30220b3a`，
  private key簽後已退休mode=`000`。正式FD-bound validator實際回傳
  `signed_evidence_verified=true`、488 tests、488 testcase、18/18 required、0 failure/error/skip；
  23-source digest仍為`db0c33...a3ed`。
- `[決策]` v5 source-authority assembler不再由一般`json.loads`單獨接收external reviews。
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/assemble_release_source_authority_v5.py`
  固定normalizer SHA-256=`46e68a...282b7`與mode=`0444`，在組裝時再次執行duplicate-key、exact schema、
  UUIDv4、A/B identity、HEAD、source digest及clean-APPROVE驗證，再交給base assembler檢查兩個reviewer
  ID獨立性。舊v11 review實測被`Strict review validation failed`拒絕；frozen assembler mode=`0444`、
  SHA-256=`29c75a59ac1f15af166e5297c52a3c55969e0b810823f182b8a06a3f76a72361`，且23-source digest不變。
- `[review closure]` internal reviewer以第三次fresh唯讀複核判定`APPROVE`：F3仍CLOSED、F5現已CLOSED，
  未找到可重現blocker。它獨立確認`validate_and_publish()` AST body只有兩個Expr calls
  `validate_pre_publish_state, rename_no_replace`，main所有staging fsync/chmod均在其前，v18/v2全部
  frozen identities、488/18 test contract與FD-bound signature驗證一致。Claim ceiling未提高，仍為
  `M2_read_level_residual_epigenetic_partition`並排除cellular subclone/lineage confirmation。此
  APPROVE只關閉supplemental F3/F5；external A/B與formal producers/signers仍是後續F1/F2 gate。
- `[狀態同步]` `InterSubMod/docs/CURRENT_FOCUS.md`頂端已更新為本線最新Task-B狀態，明示469,849
  LongPhase-S PASS scope、50,432 positional-singleton分母、M1/M2比例、BAM未覆寫、claim ceiling、
  20個formal outputs仍缺席及剩餘release gates；7/17段保留並標為歷史快照。這是其他session的live入口，
  不代表正式release已完成。
- `[資源]` `2026-07-18 18:22+08:00`發現`/dev/sdc`即時`%util=100`不是本線producer，而是四個其他
  session遺留的全workspace唯讀搜尋：34小時`find .../envs/.../python`、24小時N50廣域`rg`、19小時舊
  authority token `grep -rl`及6小時`Path.rglob`文字掃描。逐一確認cmdline、cwd、parent與fd後只對這四個
  搜尋程序送`SIGTERM`，未終止任何producer或四個one-time signers；四PID全退出，下一個`iostat`樣本為
  `r/s=0, aqu-sz=0, %util=0`。正式全量執行仍維持凍結的40/48 workers與每個BLAS一thread，不擅自提高。
- `[未決]` fresh Claude v12 Reviewer A/B既有排程仍在等待`2026-07-18 19:42:00+08:00`session reset；
  source/result/report/supplemental四個signer sessions均存活且尚未收到SIGN。只有兩份strict-normalized
  clean APPROVE JSON與v5 source-authority簽章完成後，才能依序啟動primary-v4、preflight-v8、
  cooccurrence-v7及completion chain。

## 2026-07-18 positional-singleton claim-language v19 修補

- `[目標]` 服務G2/G3；依獨立科學review修補singleton supplemental report的claim ladder與M2分母敘述，
  不改任何frozen數據、screen分類或23-role正式producer source。
- `[質疑與修補]` `InterSubMod/docs/CURRENT_FOCUS.md`不再把48個M2-determinate位點寫成預先選取；
  現明示M2套用於全部5,961個M1 flags，只有48個取得determinate PASS/FAIL，5,913個NOT_EVALUABLE，
  44,471個M1未標記位點為M2 NOT_RUN。`30/50,432=0.059486%`只稱observed operational M2-PASS
  yield，不稱biological/subclone prevalence下限；`30/48`只稱M2-determinate conditional rate。
- `[質疑與修補]` supplemental builder移除自創`L1 technical/L2 residual/L4 hypothesis`標籤，只用
  claim-contract-v5正式`M1/M2/G1/G2/R1/B1/C1/L1/L2`階梯；Markdown與HTML均新增
  candidate-level locus matched-normal未執行、exact-locus CN/purity/CCF尚未註記，因此
  B1/C1/L1/L2未達。新增reader-facing字串回歸測試，canonical最低測試數升為489。
- `[驗證]` 語意與pytest-contract專項命令：
  `env -u PYTHONPATH PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -m pytest -q research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_build_positional_singleton_report.py::test_reader_facing_report_uses_official_claim_ladder_and_determinate_language research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_build_positional_singleton_report.py::test_validate_pytest_xml_contract`；
  實際`7 passed in 2.23s`。兩支Python source compile為`syntax_ok=2/2`，`git diff --check` exit 0。
- `[驗證]` 第一次full dry run漏設numeric-library single-thread環境，觸發正式runtime fail-closed，
  實際`54 failed, 435 passed, 38 warnings`；所有54項首因均為
  `Strict producer requires ... single-threaded numeric-library settings`，不作成功證據且不修改測試。
  以正式環境
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1`
  重跑root兩組共現tests加topic全tests，實際`489 passed, 38 warnings in 39.80s`；warnings仍只為既知
  SciPy/NumPy deprecation。
- `[驗證]` 變動後獨立重算23個protected roles仍為23/23 mode=`0444`，source-set digest精確維持
  `db0c33f89c31e6b1a3e4cd68a2943df3c109f0e76ddaff40e9e0d2bf78faa3ed`。
- `[review closure]` science reviewer第二輪發現30/50,432的non-prevalence限制只被source-level測試錯誤
  綁到7/50,432 joint tumor-REF yield；修補後Markdown與HTML均把`30/50,432=0.059486%`直接稱為
  observed operational M2-PASS yield且明示不是biological/subclone prevalence或其lower bound。
  測試改為分別檢查`build_markdown()`與`build_html()`，並移除「選出的子集／高度選擇」用語；中間
  targeted run如預期先得到`1 failed, 6 passed`並捕捉換行造成的source binding缺口，修正後為
  `7 passed in 2.31s`。同一science reviewer最終回傳`Verdict=APPROVE`、無剩餘blocking語意。
- `[驗證]` final dry full suite為`489 passed, 38 warnings in 39.33s`；fresh canonical JUnit run為
  `489 passed, 38 warnings in 39.36s`。XML：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_supplemental_v19_claim_language.xml`，
  mode=`0444`、SHA-256=`db1036b43bed6e6ff0a1afee28d320ef5cc6e9ade567db49803ea1ae95457067`。
- `[決策]` v19 builder/test已凍結為mode=`0444`，SHA-256分別為
  `a5c4d736491209e4d65f163bbdd0052ab19837b0aafec56fde0ec7f01845ca07`與
  `73bf91c6246e504a4499d258da25c678e33bfef52910839fb4484037880c81e2`；同hash no-clobber frozen
  copies位於`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/frozen_sources/`。
- `[驗證]` v3 signed evidence manifest
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/supplemental_report_test_evidence.v3.json`
  SHA-256=`08569ca3ae8fc44dcc8d66797408f5471743cb8ea410f6df41b99e252e962b57`；detached signature
  SHA-256=`56cdb5da3758b9b530cd7eb6439a8a1b59f2dbff783a2a101c007f4c3dae1199`，兩者mode=`0444`。
  Public key SHA-256=`1a87955b536cf2945e8fbd5e5d433a0e851e530cd607a58f0210632a37f76ce1`，
  one-time encrypted private key簽後mode=`000`。正式FD-bound validator回傳
  `signed_evidence_verified=true`、489 tests/testcases、19/19 required、0 failure/error/skip。
- `[release readiness]` 逐一依current source-authority、cooccurrence runner與completion runner列出fresh
  targets，實測`checked=21, absent=21, existing=0`；即20個formal artifact slots加1個completion log，
  沒有其他session產生半成品。四個正式source/result/report/supplemental signer sessions均仍存活。
- `[資源]` `/big7`可用716G、`/bip7`可用1.8T、RAM available 472Gi、CPU 48；fresh第二個iostat sample
  為`sdc r/s=0,w/s=0,aqu-sz=0,%util=0`。正式參數維持frozen 40 workers及五個numeric-library
  thread env均為1，不臨時提高或改變canonical argv。
- `[驗證]` primary artifact audit canonical argv順序固定為`site-results,stable-assignments,output,
  workers=40,chunk-size=64,max-pending-chunks=80,progress-every=1000`；preflight v8固定為manifest、
  assignments、sites、InterSubMod root、independent M2 audit、primary-pre、output及`40/8/80/1000`。
  兩者均要求direct CLI `/proc/self/cmdline` exact match。
- `[review closure]` v19仍驗證historical singleton summary內的v4 source authority，與current final
  dataset/report v5 authority是刻意雙provenance。Raman獨立唯讀review回傳`APPROVE`：真實v4 chain
  `validate_source_authority_chain()` PASS；v5 dataset/report verification只要求彼此相同，不要求等於
  historical v4；initial snapshot、output receipt與prepublish guard均同時綁兩條鏈。Metadata
  `data_version=v4`正確描述singleton producer provenance，但最終reader文件應另註formal release為v5。

## 2026-07-18 v13 source-key rotation 與 authority ceremony 恢復

- `[external review transport]` v12 Reviewer A/B 使用`claude --json-schema`後仍在JSON前加入prose；
  原`normalize_external_source_review_v12.py`正確fail-closed。新增
  `audit_notes/extract_external_source_review_transport.py`與10個測試，只允許「無brace前綴＋唯一
  EOF JSON object」、拒絕duplicate/non-finite/trailing/multiple object，O_EXCL輸出mode `0444`。
  A經extractor與原normalizer後為clean APPROVE。B第一次payload的`evidence`錯為13-item array，
  normalizer拒絕；續接同一Reviewer B session只做格式更正，獨立`jq -e`確認其他14欄完全相等、
  `evidence.checks`逐項13/13相等後，原normalizer才接受。v12 A/B normalized SHA-256分別為
  `4e8caa8522c02fe437c58b959fe204af13ceb1fcaa2ded0d019c04c232aafd40`與
  `7fe6e1ee28472379dde2bff7902f185fcce00842f4a96fc9bbf02d6c860889e1`。
- `[偏離 / hard NO-GO]` 第一次v5 authority組裝成功（23 roles、source digest
  `db0c33f89c31e6b1a3e4cd68a2943df3c109f0e76ddaff40e9e0d2bf78faa3ed`），但預先啟動的source
  signer在`openssl pkeyutl -sign -rawin`回報unknown option後process退出；沒有`.sig`或`.pending`，
  隨機PASS未落盤且遺失，因此該encrypted private key不可再用，formal chain未啟動。
- `[archive]` unsigned authority SHA-256=`99d743da27779721c75bfe015b1e0d094f69345518279fe564928e9cad1df1b9`、
  approval SHA-256=`60e15e5ce92c36d025de24199e5c71b768abcd3b9d6b421b2c480cf202095aa9`、
  舊public key SHA-256=`cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920`
  已用hash/mode/absence-gated script搬至
  `InterSubMod/docs/archive/2026/07/20260718_unsigned_source_authority_v5_signer_failure_01/`；
  key directory完整搬至
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/archive/20260718_all_ssnv_v5_summary_hotfix_failed_signer_01/`。
  狀態固定為`UNSIGNED / NOT_AUTHORIZED / NEVER_USED_FOR_FORMAL_EXECUTION`，無刪除。
- `[signer verification]` 以正式`run_one_time_ed25519_signer_v2.sh`、乾淨OpenSSL環境完成獨立smoke
  ceremony：signature verify PASS、SHA-256=
  `419b87a33abb0b91d8afcda852bc2d86739dc861b4fefe07e958e5677d1c7991`、private key mode=`000`。
  新source signer session=`67803`正在等待canonical v5 approval；新public key SHA-256=
  `e91076b68495d9d879b378f3af8f16c00b3aa7bd3711943bf1dc4bff66d85cbe`、public/private modes
  `0444/0400`。result/report/supplemental三個既有v2 signer未終止。
- `[source repair]` 只在`release_source_authority.py`、cooccurrence runner與completion runner替換
  public-key hash；另關閉v12 Reviewer B的唯一hygiene finding：matched-normal runner的
  `--normal-audit` default改為`CANONICAL_NORMAL_AUDIT`，測試加入舊
  `matched_normal_methyl_tag_audit.v2.json`的superseded-token zero gate。targeted第一次因測試變數
  誤名得到`1 failed, 27 passed`，修正為`RUNNER.CANONICAL_NORMAL_AUDIT`後`28 passed in 1.42s`。
- `[new source authority anchor]` 23/23 protected roles全mode=`0444`；舊public hash命中0、舊normal
  audit token命中0；新source-set SHA-256=
  `8dac1eaab5f624568b11a36eeb8b8698d371b253a368194f22f9ca73ba5cb99b`。
  舊v12 reviews只綁`db0c33...a3ed`，不可授權新source set。
- `[tests]` key rotation後v16完整scope為`499 passed, 38 warnings in 39.47s`。為不改寫v12
  normalizer，新增v13 normalizer/schema與11個strict tests；v13 targeted=`11 passed in 0.09s`，
  最終v17 canonical=`510 passed, 38 warnings in 39.72s`。XML：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v17_key_rotation_v13_normalizer_canonical.xml`，
  mode=`0444`、SHA-256=`f11bd1427719e3c24ae0e6792df2d004736d7b3093ec44345edfcc4352f62067`、
  510 testcase、0 fail/error/skip。
- `[current gate]` fresh v13 external Reviewer A/B sessions=`43199/45811`正在獨立唯讀審查新digest、
  新key、archive isolation、stale-default closure與完整科學/command non-regression。
  `assemble_release_source_authority_v6.py`已凍結mode=`0444`、compile PASS、SHA-256=
  `6f072b6726c4bc7e7ce72c9cc2d2a534ae12277351a0b27782768abd4c52de17`，但未執行。
  只有兩份v13 normalized clean APPROVE及新detached source signature完成後才可啟動primary-v4。

## 2026-07-18 v3 FD-bound signer / v7 assembler / v14 review gate

- `[質疑 / hard stop]` Raman對v2 ceremony做fresh engineering review後回傳`REQUEST_CHANGES`，提出3個
  HIGH blocker：v2 shell signer只有`set -u -o pipefail`，`mv/chmod`等post-verify錯誤可能仍繼續並
  exit 0；key directory實際`0775`且sign/verify由path重開，沒有same-FD/no-clobber保證；v6 assembler先
  `artifact(review)`後又`path.read_text()`，review identity與payload存在TOCTOU。另有MEDIUM：
  v13 normalizer接受non-standard `NaN/Infinity`。因此v13 reviews雖為clean APPROVE，仍只保留歷史，
  未組裝/簽署authority，也未啟動任何formal producer。
- `[實作 / signer v3]`
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_one_time_ed25519_signer_v3.py`
  先以FD讀取並pin `/usr/bin/openssl`（SHA-256=`38c064...ecd08`），key dir用`0700`、key files用
  exclusive create；target/private/public/signature皆由同一FD交給OpenSSL。Signature先寫入
  `O_TMPFILE`匿名inode，再由`/proc/self/fd/<fd>`配合`linkat(..., AT_SYMLINK_FOLLOW)`以no-replace
  發布；第二次FD-bound verify、target/key/path binding與directory fsync全通過後才將private key改為
  `000`。任何例外回傳非零且不印`SIGNED`。Source mode=`0444`、SHA-256=
  `8243ced80242ace2ceae4ae71dad7c685a11f8ebe754dec56e7a1baa509c69d1`。
- `[驗證 / signer v3]` 9個tests涵蓋成功簽署、既有/racing signature collision、target replacement、
  retirement failure、signature chmod failure、post-publication directory fsync failure、錯誤命令EOF、
  existing key dir；實際`9 passed`。初版`AT_EMPTY_PATH`在無CAP環境得到`ENOENT`，正確fail-closed；
  改用Linux man-page建議的proc-FD hardlink方式後全綠，未觸及formal slots。
- `[archive / v2 signers]` 四個v2 signer全部只送EOF，均輸出
  `SIGNER_INPUT_CLOSED_WITHOUT_SIGNATURE`、exit=`7`；12個target/sig/pending slots全空、private仍
  `0400`。第一次archive script因system Python不支援`zip(strict=True)`在任何status write/move前
  exit 1；唯讀v1 ledger（SHA-256=`c1d830...e5240`）保留。相容修補後v2 ledger
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/20260718_unused_v2_signer_key_archive.v2.jsonl`
  mode=`0444`、10 events、SHA-256=
  `49860985d5c20c3801f2bb92d7124ea97d079f64f6a500a418244d697b05b94b`；4/4舊key dirs以
  `renameat2(RENAME_NOREPLACE)`完整移入各自`archive/`、mode收緊`0700`，每組有
  `UNUSED_KEY_ARCHIVE_RECORD.v1.json`，狀態`ARCHIVED_UNUSED_NEVER_SIGNED`，沒有刪除。
- `[key rotation v3]` 四個新signer sessions仍等待精確target：
  source=`18112`，public SHA-256=`e4a09d9e43f76efed408330f43cee66ac1fde457f5b4ee6889d11338935e5b6c`；
  result=`27662`，`2d37e58d837c216061ac50c6b123442593a2f87b28db71d162090c98a463d318`；
  report=`35264`，`3fb508f62ee2f476fe217f9ade8d3ee7ffcf05d431cbddb7e00e0ec92f79e585`；
  supplemental=`36796`，`a67d41ff36af6c3ff92c776ed04c2256cc0547b59d1db50013bab66f5920aca3`。
  四組key-dir/public/private modes皆`0700/0444/0400`，尚未送出任何`SIGN`。
- `[source anchor]` 新hash已回填source validator、兩個formal runners、final result/report consumers與
  current singleton supplemental consumers；23/23 protected roles全mode=`0444`，四個舊hash命中0，
  source-set SHA-256=
  `1d7166f0a192848a9b6ad812e93dac4404b65caaffaa6509fb53156c2ca8eab4`。
  相較v13 anchor，protected差異限於source/result/report trust-anchor消費者；M1/M2與科學producer未改。
- `[實作 / v14 normalizer]` 新normalizer除v13 exact-schema/duplicate/trailing/UUID/digest/HEAD/
  approval-consistency gates外，以`parse_constant`拒絕nested `NaN/Infinity/-Infinity`。Source
  mode=`0444`、SHA-256=
  `a4592b19a2cf56c9b50882f413f41cf1bd20d7596e919d7d831d9612aa73f89f`；
  normalizer+signer targeted合計`23 passed`。
- `[實作 / assembler v7]` standalone v7不載入v1/v6，也不用`importlib`、`Path.read_text()`或
  `module.artifact()`重讀review。Release module、v14 normalizer由已綁定bytes直接`compile/exec`；
  reviews的identity與JSON payload取自同一bytes；23 sources、JUnit、claim contract、Git binary、
  key files與output parent directory全保留FD到final check；JUnit attributes與testcase elements獨立
  重算；authority/approval以bound parent-dir FD + `O_EXCL`建立。Source mode=`0444`、SHA-256=
  `113345e9574725787c67f3cec888aa67328188746b2d2e5daa4ebdf5a4bc6573`；
  11個FD/race/JUnit/Git/source-digest tests實際全通過。
- `[canonical tests]` v18 dry與final JUnit皆`544 passed, 38 warnings`，warnings仍只為既知
  SciPy/NumPy deprecation。XML
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v18_v3_signer_v14_normalizer_v7_assembler_canonical.xml`
  mode=`0444`、120008 bytes、SHA-256=
  `2bce0db29e92864a740df7b2fa1492d49a624fe5fe5ca5ba3f6f55a0e998dff7`；
  attributes/elements皆544 tests、0 failure/error/skip。
- `[transport-only supplemental]` v18凍結後新增Claude CLI envelope extractor，強制command/envelope
  session binding、success state、structured object、reviewer token、duplicate/non-finite/trailing rejection與
  O_EXCL output；其12個post-v18 supplemental tests全通過。這12個不宣稱在544-test JUnit內，且不是
  protected producer；v14 normalizer與v7 assembler仍會再次完整驗證抽取結果。
- `[external v14]` schema與A/B requests皆mode=`0444`，分別綁新digest、544-test XML、v3/v7/v14
  source hashes、四把v3 keys與v2 archive ledger。第一次CLI attempt因錯誤空MCP config `{}`在review
  開始前exit 1；A/B envelope皆空SHA=`e3b0c4...b855`、stderr明示
  `expected record, received undefined`，完整保留且不授權。attempt2改用
  `{"mcpServers":{}}`、獨立O_EXCL prefix並鏈回attempt1 hashes；external A/B與internal
  Raman/Pascal fresh reviews現正進行。
- `[current hard gate]` live v5 authority/approval/signature與所有formal downstream slots仍不存在；
  本線仍是scientific release **NO-GO**。只有fresh v14 A/B clean APPROVE、Raman blocker closure、
  v7 assembly與v3 source detached signature全部通過後，才可啟動fresh primary/preflight與40-worker
  full-scope chain。

## 2026-07-19 v6 signer two-stage authority protocol / canonical JUnit gate

- `[archive / unused v5 keys]` fresh Raman review發現v5 signer仍可在terminal snapshot之後、
  Python `finally`真正return之前替換canonical target/signature並回報`SIGNED`。四個v5 signer只送EOF，
  全部exit=`7`且輸出`SIGNER_INPUT_CLOSED_WITHOUT_SIGNATURE`；沒有任何formal target被簽署。稽核腳本
  `audit_notes/archive_unused_v5_signer_keys.py` mode=`0444`、SHA-256=
  `4397c303c40e091bc6e63423bf64e846b97ab26b156b6f106b0ef02c10b14674`；ledger
  `audit_notes/20260719_unused_v5_signer_key_archive.v1.jsonl` mode=`0444`、10 events、SHA-256=
  `58aabdd0417159b42a01bda2286fd84c8fc21d839cecbd460e94717169b886f8`。四組key directory均以
  no-replace方式搬至各自`archive/`，狀態`ARCHIVED_UNUSED_NEVER_SIGNED`，沒有刪除。
- `[決策 / signer responsibility boundary]` v6不再聲稱檔案系統在caller取得結果前後具有不可改變的
  terminal state，也禁止回傳`SIGNED`。新成功事件為
  `CEREMONY_OUTPUT_AVAILABLE_REQUIRES_INDEPENDENT_VERIFICATION`，固定帶
  `consumer_verification_required=true`、`release_authority_granted=false`與
  `point_in_time_checks=PASS_BEFORE_DESCRIPTOR_TEARDOWN`。Signer的責任只限no-replace publication與
  one-time private-key retirement；正式authority只能由signer process結束後、獨立consumer重新開啟
  live target/signature/public/private/key-dir並執行crypto、identity、mode與protected-source驗證後授予。
- `[驗證 / adversarial teardown race]` 新test在private mode=`000`後，精確攔截terminal helper關閉
  reopened signature inode的`os.close`，隨即替換live signature。Signer只能回上述provisional事件；
  獨立OpenSSL consumer對live path必須非零拒絕，而被移開的原signature仍可驗證。第一次hook條件過寬，
  在OpenSSL pipe close即注入，結果`21 passed, 1 failed`且signer正確fail-closed；縮限至same-inode
  close後，專項結果`22 passed in 1.79s`。
- `[驗證 / fresh internal gate]` Raman對修補後v6做fresh adversarial read-only review，結論
  `APPROVE`、無HIGH/MEDIUM；source/test SHA-256分別為
  `e032a5c56dcb40b81b4c963413d7190a97894c7c5258c67f29700bc64503dd49`與
  `431f8412a3e86200e3d30fbabe489a514d429886487471b6bc6d187ccbce4897`。唯一LOW為production
  independent authority仍須fresh canonical與實際consumer chain另證。完整closure matrix：
  `audit_notes/20260719_internal_raman_v6_signer_two_stage_protocol_review_01.md`。
- `[偏離 / canonical JUnit circularity]` 舊v20 JUnit早於最終v9 assembler bytes 49秒，因此不能作
  最終證據。新方案不在assembler內硬編自己的hash或JUnit hash：先凍結v10 assembler bytes，再跑
  canonical tests；A/B external review payload各自帶exact `reviewed_assembler`與`canonical_junit`
  identity/counts。Assembler執行時以bound FD計算自身與live XML identity，要求兩份review完全一致、
  XML attributes與elements一致且fail/error/skip皆0。這消除self-reference，也確保JUnit晚於最終
  assembler source。
- `[current hard gate]` source/result/report/supplemental authority、approval與signature仍全部不存在；
  primary/preflight/cooccurrence及完整downstream尚未啟動。正式狀態維持
  `NO-GO / UNSIGNED / NOT_AUTHORIZED`。

## 2026-07-22 v7 signed source authority 與正式 producer 啟動

- `[canonical evidence]` final v27 canonical JUnit為`743 passed`、0 failure/error/skip；XML位於
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v27_git_env_alignment_attested_canonical.xml`，
  SHA-256=`2b191fd358d68098abdeab0664466f56061071326a60c0a28c84a4d96dce2a58`。其signed test-run manifest
  SHA-256=`b38395054c1d548e39feaabcc9a0f863a360cb429118639ed91a58be50490436`；source-set
  SHA-256=`f5ba8a9e786971c4261b51e283a0f9df6e807a8aca695bf5041271fee5420f58`。
- `[external review]` 第三輪Reviewer A/B皆使用`claude-opus-4-8`、max effort與唯讀工具，對同一HEAD、
  source-set與canonical evidence回傳`APPROVE`、blocking findings空、前輪F1/F2已解決。兩份process
  attestation SHA-256分別為`87616675b7165e196dcac7407c0415eabe5c0b813471f55c69eb63e1a84f6e98`
  與`5333fe461463dad40ed91adfa9c9a89fdd61bd5dfeadc4b64c04addceef41e5b`，detached signatures均經
  `/usr/bin/openssl`獨立驗證。剩餘LOW/INFO只涵蓋同UID/root/ptrace、bootstrap trust及尚未由source
  chain建立的生物統計claims，不是formal producer blocker。
- `[source authority v7]` signed authority為
  `InterSubMod/docs/provenance/source_authorities/20260722_all_ssnv_focal_alt_release_source_authority.v7.json`，
  SHA-256=`ab68a9d3e8578545e6f930e8f8c79912a132f02692adec3b19542276605b6a0c`；approval
  SHA-256=`a4f2ac6a7605842122a8a59e47e4867fbef0e1fbfb1aa668eaa17475c504ca99`；signature
  SHA-256=`7e0b39abca65451de1e68b6f756f98b75f675137e63dde212935d8aeae38aa44`。獨立consumer
  驗證`pass=true`、23/23 roles、2 signed reviews、exact source-set與147個FD lease descriptors；release
  private key已退役為mode`0000`。
- `[archived superseded evidence]` v25因兩個runner仍殘留v6/v9 authority literals而封存至
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260722_v25_superseded_by_v7_consumer_rotation/`；
  v26因live Git subprocess未明確帶clean env而封存至同層
  `20260722_v26_superseded_by_git_env_context_fix/`。兩者均保留完整source snapshot與SUMMARY，未刪除、
  未作current authority。
- `[primary artifact audit]` 正式absolute-entrypoint命令以40 workers重算`102,842` stable sites、
  `102,842` assignment objects與`308,526/308,526` artifacts，receipt
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream.json`
  為`pass=true`、`formal_task_b_release_eligible=true`，artifact-set SHA-256=
  `195e77d571576f37debf1627edb57f9dc7edb935849bd0ae6e29b323b380b2ca`。第一次repo-relative
  entrypoint被direct-CLI canonical-process gate正確拒絕，無JSON產物，cache保留於
  `audit_notes/failed_formal_runs/20260722_primary_v6_relative_entrypoint_rejected/`。
- `[current runtime]` full raw-identity preflight v9正以40 workers、chunk-size 8、max pending 80執行全部
  `102,842` tasks；session=`7921`，預定輸出
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_task_contract_preflight.v9_command_parity_full_runtime.json`。
  它必須證明每個result row的SAM core、除RG外所有typed auxiliary tags、latest HP/PS、assignment label與
  partner calls符合identity invariants，且8個protected consumers的source identity before/after完全相同。
  preflight PASS前不得啟動v8 cooccurrence；v8完成前不得啟動terminal completion chain。

## 2026-07-22 positional-singleton supplemental v4 test evidence 修復

- `[pre-decision / isolated hotfix]` 唯讀agent發現singleton report builder仍綁定已封存v25 JUnit，會使
  supplemental test-evidence v4無法組裝。修正範圍只包含
  `scripts/build_positional_singleton_report.py`的canonical XML常數與
  `tests/test_build_positional_singleton_report.py`的一項反漂移測試；兩檔修後均恢復mode`0444`。
  23個v7 main producer sources、active preflight與result/report keys均未修改。
- `[targeted verification]` singleton report專項測試實際`69 passed in 7.46s`；完整51-file dry suite實際
  `744 passed, 34 warnings in 77.46s`。warnings只為既知SciPy/NumPy deprecation。
- `[canonical evidence v28]`
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_supplemental_v28_v7_authority_alignment_canonical.xml`
  為mode`0444`、size=`79,647`、SHA-256=
  `ce6de8a04371e80b39e9520312a3346fc2ab99ce29a04a56f42413b6e16a2450`；獨立parser重算
  tests/testcases=`744/744`、failures/errors/skipped=`0/0/0`、required supplemental tests=`19`。
- `[signed supplemental test evidence]`
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/supplemental_report_test_evidence.v4.json`
  為mode`0444`、SHA-256=`d84c69cf97068656ca19b73b87ae851eaa3126bf541fc9c0d14d9bc3a2db95c6`；
  detached signature SHA-256=`73a6183e88281f152e61fec6306e53e7e62ef2d24977e59819a7a8f977db4bec`。
  一次性signer成功後private key mode=`0000`；獨立OpenSSL與report consumer重開live paths均PASS，
  `signed_evidence_verified=true`。這只授權supplemental report的測試證據，不代表主Task-B或生物claims已完成。

## 2026-07-22 positional-singleton supplemental v5 v7-inventory recovery（權威修正）

- `[修正 / v4 invalidated]` 上一節的v4/744-test證據是歷史紀錄，不再具有release authority。新增的
  regression test造成v7 `test_sources` identity drift，因此已恢復v7 manifest所列的精確50個test source，
  並將v28 XML、v4 manifest與signature完整移至
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260722_supplemental_v4_invalidated_by_v7_test_inventory_restore/`；
  未刪除任何檔案，舊一次性key維持退役且不得重用。先前紀錄的「51-file」也是摘要誤植，權威manifest
  實際為50個test source。
- `[canonical evidence v29]` 以v7 manifest的50個test source、隔離Python環境與新的pycache目錄執行，
  正式結果為`743 passed, 34 warnings in 74.99s`；warnings僅為SciPy/NumPy deprecation。XML
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_supplemental_v29_v7_test_inventory_recovery_canonical.xml`
  為mode=`0444`、size=`79,536`、SHA-256=
  `20fa8cc630c8dde9566433197dd0341f15df8aa9cff1da62f3e1b5de0eca8e27`；builder validator重算
  tests/testcases=`743/743`、failures/errors/skipped=`0/0/0`、required supplemental tests=`19`。
- `[signed supplemental test evidence v5]`
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/supplemental_report_test_evidence.v5.json`
  為mode=`0444`、size=`4,386`、SHA-256=
  `d1b18adc900d9d5d80c5b86f1c82b004a1629cd05dd9ae1aa9aeb5fe7630fe67`；detached signature為
  64 bytes、mode=`0444`、SHA-256=
  `2a7e20266c3f16bca05fd82f61e1fe2b8051eb5c056bc3fad4bead25d4d016f4`。一次性private key已退役
  為mode=`0000`；獨立`/usr/bin/openssl pkeyutl -verify`與report consumer均PASS，回報
  `signed_evidence_verified=true`。此v5是目前唯一有效的singleton supplemental test authority；主Task-B
  與生物claims仍須等待full raw-identity preflight、cooccurrence及terminal completion chain完成。

## 2026-07-22 fresh positional-singleton schema-2 audit與identity-schema修補

- `[fail-closed / root cause]` fresh audit前兩次都在atomic publication前拒絕、正式output directory不存在。
  v7 canonical authority identity只定義`path/size_bytes/sha256/mode`，但singleton auditor本地artifact另保存
  `mtime_ns`並以whole-dict equality比較，因此cryptographic identity完全相同仍被誤判drift；第二個
  bound-FD comparison也有同一問題。修正為單一`canonical_identity()`只比較四個canonical fields，
  canonical expected record仍須exact四欄schema；其他同來源before/after的`mtime_ns` race detection保留。
  修後source mode=`0444`、SHA-256=
  `f865344811c42ec7c86990ebee6191f20aab3de08de15963638e0651a61564df`，既有audit tests兩輪皆
  `16/16 passed`，沒有新增或修改v7 test source；主v7 consumer再驗證23 roles、147 FD leases與
  source-set `f5ba8a9e786971c4261b51e283a0f9df6e807a8aca695bf5041271fee5420f58`仍`pass=true`。
- `[fresh atomic audit]` 修補後使用v10 screen、primary v6、tree audit、independent M2 v3、signed v7
  authority與claim-contract-v5重跑全部資料，原子發布至
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/positional_singleton_methyl_multigroup_audit_v2_source_authority_v7/`。
  四個輸出皆mode=`0444`：summary SHA-256=
  `755e1df53afeff0f0d5504dc5f763aecdeffebfabeb6064582ba818c21cd6818`、site audit=
  `f8f21190e0b03cf6354523d6669e897aacf49f5944ed7f47ee5a4f499ada10eb`、30-case table=
  `2fe42d55924d4d73f8a0d7b436cc8c517e9d06834da9fd7a1e882ce24c4f1dd0`、`_SUCCESS`=
  `999dbae73b5d92985c67af1c2babab527c1eb7e5d35c62e1ecfb1384903f4081`；marker內嵌hash與live files一致。
- `[recomputed results]` singleton=`50,432/469,849=10.733661%`；M1 evaluable=`48,347`、flagged=`5,961`
  （all-singleton yield=`11.819876%`、evaluable conditional=`12.329617%`）；M2 PASS/FAIL/NOT_EVALUABLE/
  NOT_RUN=`30/18/5,913/44,471`。`30/48=62.5%`只是在selected M2-determinate subset中的conditional
  rate；`30/50,432=0.059486%`才是all-singleton operational yield，兩者都不是cellular subclone
  prevalence。39項checks全true，confirmed cellular subclone與linear ancestry仍皆為0。
- `[independent review]` 唯讀agent Meitner對修補與release transitive binding判定`APPROVE`，無
  blocking/high/medium finding。v5 signed test evidence只綁canonical JUnit、report builder與其test，
  不宣稱凍結singleton auditor；fresh audit summary會記錄修後auditor identity，report builder會核對
  live identity，最終supplemental receipt再transitively seal，因此不需v5 key/test-evidence rotation。
- `[independent data recount]` 唯讀agent Faraday由raw TSV獨立重算全部分母、key uniqueness、gzip與
  `_SUCCESS` hashes，完整性blocker為0。30個M2 PASS全為TP；per-dataset依序為COLO829/H1437/H2009/
  HCC1395/HCC1395_DORADO/HCC1937/HCC1954=`0/2/9/2/2/14/1`，其中HCC1937+H2009占`23/30`；
  HCC1395與DORADO四個PASS沒有同一biological-site overlap。29/30為2 groups，唯一3-group案例為
  `H2009 chr6:52060111 C>A`；30/30六個technical axes均LOW_EFFECT，HP aligned為0。M2 determinate僅
  `48/5,961=0.805234%`，即99.1948%的M1 flags無determinate M2，因此`30/48=62.5%`高度selection-
  conditional。這些是必須保留的發表限制；只要報告維持read-level residual epigenetic partition ceiling，
  就不是release data blocker。

## 2026-07-22 positional-singleton standalone QA schema 2.4 frozen-source closure

- `[independent review blocker]` 唯讀agent Hooke先判定`NOT APPROVE`：含border的box ratio會誤殺
  正常`2237x883`寬圖；固定5% non-dominant pixel門檻會誤殺正常稀疏科學圖；預期失敗receipt重建
  payload時會遺失input、executed copy、QA script、browser與runtime的after identity。這三項分別為
  1 HIGH、2 MEDIUM，未完成前禁止使用formal HTML QA。
- `[修正]` `scripts/qa_positional_singleton_standalone.py`升至schema `2.4.0`：只以扣除CSS border的
  `clientWidth/clientHeight`判原始比例，box ratio只保留描述；像素門檻改成dense直接通過，或
  sparse圖必須同時達到pixel、bbox、8x8 grid與unique-color門檻；failure path保留既有payload並在
  早期失敗時補抓全部after identities。`scripts/test_qa_positional_singleton_negative_mutations.py`升至
  schema `1.6.0`並逐child強制after identity相等；fixture加入`2237x883`多色稀疏寬圖canary。
- `[frozen sources]` 三檔均已mode=`0444`：QA script size=`97,547`、SHA-256=
  `5fe5af3048de09a344dfefefed4df3fa185e39763d6bb2ed47852a4ebcbea0a7`；negative harness
  size=`26,499`、SHA-256=`f120913ec750776eec43822483800d52d483d73d4786167c6eb971d89ac9be0f`；
  fixture size=`5,980`、SHA-256=`f6aebb0f7ab086bc08b4fffa88c810496ae71a6ff816a9803eaf17726e3d1122`。
- `[positive frozen command/output]` command為
  `/bip7_disk/liaoyoyo2001/miniconda3/bin/python3.9 -I -X pycache_prefix=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/.python_cache_qa_hotfix_v1 InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/qa_positional_singleton_standalone.py --html InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/fixtures/qa_positional_singleton_standalone_valid.html --executable-path /bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium_headless_shell-1217/chrome-headless-shell-linux64/chrome-headless-shell --output-dir InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/qa_positional_singleton_standalone_selftest_v12_frozen_source --evidence-kind self_test --wait-ms 50`。
  實際5/5 profiles PASS；320px profile的寬圖content error=`0.4%`、border-box error=`1.1%`、
  non-dominant=`3.5%`、bbox=`67.4%`、grid=`40.6%`、unique colors=`644`，以
  `sparse_spatially_distributed`通過。receipt mode=`0444`、size=`1,535,607`、SHA-256=
  `d155d91ad2ff594cb6d7f8b975fcae06bddff8198db3f1a8fee6794e08d2abc1`。
- `[negative frozen command/output]` 同一Python執行
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/test_qa_positional_singleton_negative_mutations.py`，輸入上述frozen QA/fixture/positive receipt與固定Chromium，
  輸出`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/qa_positional_singleton_standalone_negative_mutations_v9_frozen_source/`。
  實際10/10 mutations拒絕、10/10 complete after identity、2/2 publication failures quarantine、
  1/1 unsigned formal context拒絕；97% white sparse mutation為non-dominant=`3%`、bbox=`88%`、
  grid=`75%`但只有2 colors，故仍正確拒絕。receipt mode=`0444`、size=`40,807`、SHA-256=
  `d0f24918be2982ffb7cfbc339a93b5b7169426b256a638babfd41582c3c985e8`。
- `[independent closure]` Hooke以修後來源、正向v11與負向v8先行證據重驗，判定`APPROVE`、沒有
  新blocking/HIGH/MEDIUM；主agent再以frozen v12/v9獨立重算同一契約PASS。這仍是Task-F harness
  self-test，不能取代實際supplemental PNG/HTML/PDF與簽章後formal QA；低密度精確雙色圖會刻意被
  拒絕，適用範圍限定本報告的多色Matplotlib圖。

## 2026-07-24 Task-B signed release bounded stop

- `[決策]` 科學資料與分析維持完成，但停止無限延伸的 schema-recovery/signing 修補。允許一次 bounded
  attempt 關閉第二輪 C/finalizer key split-brain；若第三輪仍有 blocking finding，就不進第四輪。
- `[已修]` C runner 組合時 exact-once 替換3個舊key paths與2個舊public-key digests；新增direct v30
  finalizer綁定result/report-v9 signer targets。重寫後prefix實際shell exit=`0`、stderr空白。
- `[驗證]` Frozen validator SHA-256=`74c3e445fec567cf81184611b983d260b1b60344b8def63fd2feac099004daf1`；
  11-source set SHA-256=`83bf28b8c9455df8bc74390841ff235de39e3e78b486537a6decc6877ddf8c70`。
  Exact read-only probe exit=`0`、`pass=true`、`no_output_writes=true`、canonical=`733 passed`；獨立
  sidecar=`775 passed`。
- `[第三輪 blocker]` Mendel發現`build_all_ssnv_report_artifact_schema_recovery_v29.py`仍遞迴呼叫
  v29 finalizer及已撤離的result/report-v6 keys；1 HIGH + 1 MEDIUM，verdict=`REQUEST_CHANGES`。
  既有test只執行C prefix，未覆蓋nested report-builder dependency。因此不建立authority、不啟動
  Nash或external Claude正式審查、不執行V/R/C。
- `[安全停機]` result/report signer均以KeyboardInterrupt結束；formal reviews、authority bundle、兩個
  release receipts與signatures皆不存在。兩把未使用private keys仍為mode=`0400`、link count=`1`、
  size=`290`，未簽署、未消耗。
- `[狀態上限]` 7 datasets / 469,849 same-run LongPhase-S recalibrated PASS sSNV 的資料分析結論不變；
  可稱data analysis validated，不能稱signed Task-B release完成。稽核紀錄：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/rejected_pre_authority_reviews/20260724_v30_round3_transitive_report_builder_finalizer_split_brain/SUMMARY.md`。
