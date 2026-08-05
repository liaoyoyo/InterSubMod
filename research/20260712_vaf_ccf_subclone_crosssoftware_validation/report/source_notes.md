<!--
建立時間: 2026-07-12
目標: 定義 VAF/CCF/外部 subclone 驗證 HTML 的報告工作、證據 spine、圖表契約與 omitted-claim 理由
處理範圍: 技術讀者；chr1-22；7 dataset rows／6 biological cell lines；HCC1395 cross-source pair
關聯檔案:
  - InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/pre-decision-audit.md
  - InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/
狀態: in_progress
-->

# HTML report source notes

## Reporting job

- **User question**：七個 dataset 的 VAF／subclone 分佈是否相似；兩個 HCC1395 技術來源是否反映相同 clone/subclone 現象；外部知名工具是否支持；這能否驗證 InterSubMod 方法。
- **Audience**：technical（方法、統計與驗證為主要價值）。
- **Delivery mode**：單一 self-contained portable HTML，由 canonical `artifact.json` 與 Data Analytics portable builder 生成。
- **Scope**：chr1–22；latest LongPhase-S recalibrated PASS tree-backbone VCF 為 primary；canonical caller-PASS baseline 為 sensitivity；historical layered-v2 topology 僅作有標籤的對照；clean layered-v3 topology release 未完成時不可取代 historical 結果。
- **Comparison basis**：HCC1395 vs HCC1395_DORADO 是 same-cell-line cross-source technical comparison，不是 biological replicate。其餘 cell lines 只比較 profile summary，不比較 mutation-cluster identity。
- **Decision-useful success**：raw shared-locus effect sizes與CI、CN/purity-aware external clusters、mutation→region bridge、cause-decomposition sensitivity與資料 gate 指向一致；明確區分 technical reproducibility、computational concordance、identifiability drivers、biological truth。

## Answer-first report spine

1. HCC pair raw VAF spectrum與逐位點結果是否一致，差異相對二項抽樣底噪多大。
2. 七個 datasets 的 truth-confirmed raw VAF distribution 如何分群；HCC pair 是否互為最近鄰。
3. PyClone-VI joint/separate fits 是否得到相近 cluster prevalence／assignment；CN rounding敏感度是否改變結論。
4. 固定 5,720 regions 中，site/allele identity 與 k=2、k=3、k>=4 的 coarse/read/VAF/shape/exact endpoints 如何分離；fixed 與 conditional 分母是否一致揭露。
5. Independent PyClone labels 能否接回 regional mutations/partitions/read directions；perfect/high agreement 是否由 selection 或 single-cluster dominance 造成。
6. Site sampling、k、candidate multiplicity、HP mapping、depth、|delta VAF| 與 shared CN/LOH 哪些只相關、哪些在 multivariable diagnostic 後仍穩定。
7. InterSubMod 的 coarse shape、mutation-labeled forest與外部 CCF clustering在什麼解析度一致或不一致。
8. Fang/SEQC2 已知 HCC1395 S1–S10 heterogeneous architecture 與本輪結果是否「相容」，而非 exact validation。
9. 哪些缺口阻止 tree truth／method accuracy claim，以及最短補強路徑。

## Required technical-report structure mapping

| Required role | Visible report section |
|---|---|
| Title | `HCC1395 跨來源 VAF／CCF 與 subclone 外部驗證` |
| Technical summary | `結論先講：頻率結構高度再現，真實演化樹尚未被證明` |
| Key findings + visuals | raw VAF relationship、7-sample distributions、PyClone cluster results、k-stratified resolution chart、clone-region bridge tables、cause audit table、topology endpoint table |
| Scope/data/definitions | raw VAF、CCF、C、T、Topo、shared-locus、truth-confirmed、primary/sensitivity |
| Methodology/model | exact allele join、count reconstruction、block bootstrap、SAVANA CN rounding、PyClone beta-binomial |
| Limitations/robustness | shared CN、union recount gap、purity conflicts、historical topology、no single-cell truth |
| Recommended next steps | DORADO source-specific CN、BAM union recount、within-source split-half、orthogonal truth |
| Further questions | passage/source identity、multiplicity、cluster-to-regional-edge mapping |

## Chart map

| Segment | Question | Family / native type | Data fields | Claim supported | Palette policy |
|---|---|---|---|---|---|
| 7-sample distribution | HCC pair是否互為最近鄰，其他cell lines是否不同 | Distribution / `line`（50-bin density） | sample, VAF_bin, density, source, truth_scope | aggregate raw-VAF spectrum | relaxed multi-category，≤5 roots + line style/open tones；HCC pair直接標示 |
| VAF composition | 各dataset落在相同VAF bands的比例 | Composition / `stackedBar100` | sample, band, fraction, n | distribution shape差異，不當clone labels | relaxed multi-category；bands有序 |
| HCC locus concordance | 逐位點是否沿 y=x | Relationship / `scatter`（deterministic bounded sample） | locus, VAF_A, VAF_B, DP_A, DP_B, truth | exact-locus technical concordance | single blue root + neutral y=x |
| HCC difference | 誤差與VAF大小是否有關 | Relationship / `scatter`（Bland–Altman bounded sample） | mean_VAF, delta, DP | bias與heteroscedasticity | single orange root + neutral/bias references |
| PyClone joint clusters | 兩來源每個共同cluster的CCF是否接近 | Comparison / grouped `bar` | cluster, sample, CCF, mutation_n, assignment_prob | conditional external-cluster prevalence reproducibility | hard two-root cap |
| All-ready PyClone | 不同cell lines的cluster count/clonal fraction/entropy | Comparison / horizontal `bar` + exact table | dataset, active_clusters, clonal_fraction, entropy, status | profile summary only | categorical≤5 roots；blocked neutral |
| Regional k resolution | 位點數增加時，coarse與read/VAF細粒度endpoint是否同向 | Comparison / grouped `bar` + full audit table | k_stratum, endpoint, value, numerator, stratum_fixed_denominator, evaluable_denominator, denominator_basis | site-set/allele高度一致但exact/substructure隨k下降；k>=4 coarse受unresolved↔unresolved 888/1,093堆高 | relaxed four-category；denominator以tooltip與鄰接文字明示 |
| Endpoint evidence ladder | raw/coarse/exact/CCF各自證到哪裡 | Table | endpoint, denominator, agreement, independence, claim ceiling | resolution-dependent reproducibility | neutral table；不以顏色代替語意 |

## Omitted or bounded outputs

- PNG/SVG figures由分析程式產生作 QA/supporting evidence；正式 HTML 使用 canonical native charts，不嵌入第二套 chart runtime。
- 不畫「唯一真實演化樹」：現有資料與外部工具都不支持此輸出。
- 不把 C 當 clone count：C 是 primary HP 中 read-supported mutation-state groups。
- 不把 raw VAF bands標成 clonal/subclonal；只有 PyClone cellular prevalence可在模型條件下使用 CCF 語意。
- COLO829/HCC1937 不產 CN-aware PyClone 結果：可信 allele-specific CN 缺失，顯示為 blocked。
- Pairtree 若未有 within-source ceiling、source-specific CN與穩定 clusters，不執行或只列 future conditional step；避免最高分單樹誤導。
- Clone-region bridge 與 cause-decomposition 不新增 chart：join/selection/partition/edge 與 unadjusted/adjusted estimates 的分母不同，使用三張小型精確 audit tables，避免柱高掩蓋 vacuous selection 或 confounding。
- `assignment>=0.8` regional Jaccard=1 明示為 subclonal union=1 的退化結果；partition exact 96.9% 明示由 5,007 both-single regions 主導。

## Visual design direction

採 **scientific field notebook / restrained editorial**：高密度但單欄閱讀、深墨文字、藍／赭雙根標示 HCC pair、其他樣本退為低飽和分類色；數字以 compact metric cards，方法與 caveat 緊貼證據。由 portable reader處理 light/dark、responsive、來源 affordance 與 semantic fallback，不另寫平行 HTML/CSS runtime。
