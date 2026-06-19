# 研究任務樹（WBS） — TASKS.md

> **自動產生，請勿手改。** 唯一真值 = `state/tasks/graph.json`；改後重跑 `python3 scripts/tasks/task_graph.py --graph state/tasks/graph.json render`。
> 生成 2026-06-20 · build `c660179` @ `research/subclonal-reconstruction-202606` · 節點 83（✅34 ◐20 ⛔1 ☐27 ✗1） · focus = `T-SL`

## ⓪ 驗證（validate + check）
- validate: ✅ PASS
- check: 1 findings
  - DRIFT `T-C1`: 指向的 cycle 20260602-1521-g1-zar1l-brca2-asm-survey 在 active.json 已 stale 18 天

## 聚焦路徑　論文（碩論：Subclonal reconstruction） ＞ 完成 ISM 程式 ＞ 修正「結構驅動切區塊 × 標籤驗證」

## 任務樹（含括 parent；縮排=層級）
- ◐ `T-PAPER` 論文（碩論：Subclonal reconstruction） [里程碑/進行中]　↳project_thesis_writing_architecture
  - ◐ `T-ASM` ASM characterization [里程碑/進行中]
    - ☐ `T-C-CROSS` 跨樣本 ASM 正式圖表（R6） [分析/待辦]　in:●6 樣本 excess …　缺:待 G-A 跨樣本統計　↳project_cross_sample_asm_reproducibility
    - ✅ `T-C1` ZAR1L/BRCA2 ASM 驗證 ⭐3 POSITIVE [程式/已完成]　in:●tagged BAM,●somatic VCF　↳20260602-1521-g1-zar1l-brca2-asm-survey
    - ✅ `T-C2` 跨 6 樣本 ASM 復現 ⭐3 [程式/已完成]　in:●ISM TP/FP re…　↳project_cross_sample_asm_reproducibility
    - ✅ `T-C3` ISM 完整 TP/FP/FN + cis 存在性 ⭐4 [程式/已完成]　in:●tagged BAM,●somatic VCF　↳project_ism_complete_tpfpfn_existence_cis
    - ✅ `T-C4` ASM×CN confound pilot ⭐3 [分析/已完成]　in:●HP-axis 甲基矩陣,●CN call　↳project_asm_cn_confound_pilot
    - ✅ `T-C5` ASM 位點顯示頁 + CramersV 閘控 [分析/已完成]　in:●significance…　↳project_asm_locus_display_and_cramersv_reliability_gate
    - ☐ `T-C6` COLO829 單樣本 ASM 補強（解單樣本上限） [程式/待辦]　in:●normal 甲基 ba…,●COLO829 tagg…　缺:COLO829 normal 甲基缺　↳project_zar1l_brca2_asm_verification
  - ◐ `T-ASSET` 資料資產 [里程碑/進行中]
    - ✅ `T-A1` 6 樣本 tagged BAM + somatic VCF + ISM TP/FP 齊備 [資料/已完成]　in:●HKU/DORADO O…　↳project_subclonal_reconstruction_paper_focus
    - ◐ `T-A2` 6 樣本 normal 甲基（5/6 ready，COLO829 缺） [資料/進行中]　in:●tagged BAM,●normal BAM (…　缺:COLO829 normal 甲基缺；6 normal 全 zhenyu112 帳號 = SPOF　↳project_subclonal_reconstruction_paper_focus
    - ☐ `T-A3` 6 normal 甲基帳號 SPOF 備份 [資料/待辦]　缺:6 normal 全 zhenyu112 帳號 = 單點失效　↳project_subclonal_reconstruction_paper_focus
    - ⛔ `T-M1` 建立 6 樣本 normal 甲基 baseline [程式/阻塞]　in:●normal 甲基 5m…　缺:COLO829 normal 甲基缺 → 無法全 6 樣本　↳project_subclonal_reconstruction_paper_focus
  - ◐ `T-EXT` 外部驗證 [里程碑/進行中]
    - ✅ `T-E1` 外部文獻驗證庫（61 源親讀，稽核 CLEAN） [資料/已完成]　in:●文獻 PDF/repo …　↳project_external_validation_library
    - ◐ `T-E2` 6 樣本 clone/subclone 外部真值盤點 [分析/進行中]　in:●6 樣本清單,●外部 truth (Fa…　缺:H1437/H2009/HCC1937/HCC1954 缺專門解答（需自證）　↳project_six_sample_clone_subclone_external_truth
    - ☐ `T-E3` citation-verification（投稿前） [撰寫/待辦]　↳project_external_validation_library
  - ◐ `T-GATE` 決策 gates（OPEN） [里程碑/進行中]
    - ☐ `T-GATE-GA` G-A：跨 6 樣本 → 定單樣本 ⭐3 vs ⭐4 [分析/待辦]　缺:待 COLO829 normal 甲基；G-A 統計未跑　↳project_methyl_phasing_assist_line
    - ☐ `T-GATE-GB` G-B：within-hap somatic null → 定甲基-subclone 故事 [分析/待辦]　缺:未跑前甲基-subclone 只能寫存在性窄+負　↳project_apriori_subclone_classification_model
    - ☐ `T-GATE-HD1` HD-1：R-SELFREF → 定論文 Grade [分析/待辦]　缺:未決；影響論文主張強度　↳project_subclonal_reconstruction_paper_focus
  - ◐ `T-ISM` 完成 ISM 程式 [里程碑/進行中]　↳project_ism_method_soundness_validation
    - ☐ `T-CSC` 驗證判讀確認 → clone / subclone [分析/待辦]　缺:ISM 完成才開跑；六層框架 ⑥ subclone 下游　↳project_apriori_subclone_classification_model
    - ◐ `T-S1` a-priori 4-pop subclone 分類模型 ADOPT_WITH_CORRECTIONS [分析/進行中]　in:●TP/FP/FN+cis…,●haplotag 1-1…　缺:B 排序 illustrative / somatic 未定待 G-B　↳20260617_keep_remove_classification_conditioned_axes
    - ✗ `T-S2` tumor-only 非監督軸 NEGATIVE（勿再開） [分析/已放棄]　in:●a-priori 分類　↳20260617_tumor_only_unsupervised_axis_negative
    - ◐ `T-S5` 甲基『幾群』判定 — 收斂 Path B model-based [分析/進行中]　in:●a-priori 軸　缺:A 路相關感知 null 真實資料失敗 81% → 需 model-based 大改　↳project_subcluster_cluster_count_determination
    - ✅ `T-S6` cluster×label 對齊 = paired 非 tumor-only（盤點） [分析/已完成]　in:●a-priori 分類　↳project_cluster_label_alignment_readset_paired
    - ★ `T-SL` 修正「結構驅動切區塊 × 標籤驗證」 [里程碑/進行中]　缺:救回主軸 allele≫HP 需 cis-control 分 cis-ASM vs subclone　↳20260617_keep_remove_classification_conditioned_axes
      - ✅ `T-SL-ARI` 觀察：C2 ARI（幾何群×標籤）genome-wide [里程碑/已完成]
        - ✅ `T-SL-ARI-F` 修：emit ARI_Cluster_HP / ARI_Cluster_Allele C++ [程式/已完成]　in:●幾何群 labels +…,●Bootstrap::a…　↳project_ism_verdict_false_negative_audit_2026_06_16
        - ✅ `T-SL-ARI-O` 觀察：ARI median HP0.145/allele0.005；對齊23.9%/沒對齊49.3% [分析/已完成]　↳project_ism_verdict_false_negative_audit_2026_06_16
      - ◐ `T-SL-C1` ① within-HP bootstrap-stability 選 k [里程碑/進行中]
        - ☐ `T-SL-C1-F` 修：enable_bootstrap=true + within-HP 用 bootstrap [程式/待辦]　in:●RegionProces…,●Bootstrap cl…　缺:C1 spec 20260618_within_hp_bootstrap_and_cluster_label_ari_cpp_change_spec_01.md pending_approval　↳project_ism_verdict_false_negative_audit_2026_06_16
        - ✅ `T-SL-C1-O` 觀察：raw silhouette 過切（within-HP substructure 30.9%） [分析/已完成]　↳project_ism_verdict_false_negative_audit_2026_06_16
        - ☐ `T-SL-C1-V` 驗：consensus 重抽只留穩定群 [分析/待辦]　缺:待 C1-F 落地後跑　↳project_subcluster_cluster_count_determination
      - ☐ `T-SL-C3` ④ C3 HP-fine 組合對齊測試 emit [程式/待辦]　in:●HP-fine sub-…　缺:C3 spec pending；落地序 C1→C3　↳project_ism_verdict_false_negative_audit_2026_06_16
      - ✅ `T-SL-PERM` ② PERMANOVA + PERMDISP 確認對齊哪標籤 [里程碑/已完成]　↳project_ism_verdict_false_negative_audit_2026_06_16
        - ✅ `T-SL-PERM-F1` 修：enable_dispersion=true + analytic F-p [程式/已完成]　in:●ISM C++ Regi…,●距離矩陣/cluster…　↳project_ism_verdict_false_negative_audit_2026_06_16
        - ✅ `T-SL-PERM-F2` 修：D1 統一結構判定（棄 F1-Significant） [程式/已完成]　in:●clean_locati…　↳project_ism_verdict_false_negative_audit_2026_06_16
        - ✅ `T-SL-PERM-O` 觀察：8180 legacy 91.4% 顯著但 87% dispersion-混淆 [分析/已完成]　↳project_ism_verdict_false_negative_audit_2026_06_16
      - ✅ `T-SL-RECLS` ③ reclassify-v2（Option-B 家族） [里程碑/已完成]
        - ✅ `T-SL-RECLS-1` Stage① 連續強度分數（取代飽和 HeuristicScore） [程式/已完成]　in:●region 結構統計　↳project_ism_verdict_false_negative_audit_2026_06_16
        - ✅ `T-SL-RECLS-2` Stage② within-HP 結構軸（pattern + level bimodality） [程式/已完成]　in:●within-HP pa…　↳project_ism_verdict_false_negative_audit_2026_06_16
        - ✅ `T-SL-RECLS-4` Stage④ 新 VC（任一結構證據→保留） [程式/已完成]　in:●Δβ_sig/HP_AU…　↳project_ism_verdict_false_negative_audit_2026_06_16
        - ✅ `T-SL-RECLS-BUG` 🐛 修 MAX_DIST clustering heap-OOB segfault [程式/已完成]　in:●TreeCutter::…　↳project_ism_sampling_layer_review
      - ☐ `T-SL-RO` ⑤ tumor-only 獨立重跑 + 兩層事實表 [里程碑/待辦]
        - ☐ `T-SL-RO-1` P-readset：--normal-bam 是否 optional → tumor-only rerun [程式/待辦]　in:●tumor reads …,●--normal-bam　缺:先驗 --normal-bam 是否 optional（可能零改碼）
        - ☐ `T-SL-RO-2` P-facttable：tumor+normal 兩層事實表（只記不判） [分析/待辦]　缺:待 P-readset
      - ☐ `T-SL-S3` Stage③ coverage-peak CN 自估 + SEQC2 驗證腳本 [程式/待辦]　in:●coverage pea…,●SEQC2 BED（只驗…　缺:🔴 SEQC2 不可當分類輸入（用戶定）；原載 SEQC2 BED 當輸入作廢重設計　↳project_ism_verdict_false_negative_audit_2026_06_16
  - ◐ `T-LOCUS` 位點層驗證 [里程碑/進行中]
    - ✅ `T-L1` chr2:18M subclone 位點驗證 L2 [分析/已完成]　in:●tagged BAM (…,●6 sSNV 候選　↳20260615_chr2_18M_subclone_independent_validation
    - ✅ `T-L2` chr2:18M × SEQC2 外部驗證 + AI 解釋指南 [分析/已完成]　in:●chr2:18M L2 …,●SEQC2 truth　↳project_hcc1395_chr2_18M_subclone_external_validation
  - ◐ `T-METHOD` 方法健全性 [里程碑/進行中]
    - ✅ `T-SUP-CLUST` clusterability vs coverage/CN（支撐 ⭐3） [分析/已完成]　↳project_subclonal_reconstruction_paper_focus
    - ✅ `T-SUP-COPY` 甲基×copy confound（SEQC2，支撐 ⭐3） [分析/已完成]　↳project_asm_cn_confound_pilot
    - ✅ `T-SUP-LOCUS` 單位點甲基組合詮釋窮舉（支撐 ⭐3） [分析/已完成]　↳project_O12_loh_methylation_scenarios
    - ✅ `T-V1` ISM 方法健全性驗證 SOUND_WITH_CAVEATS [分析/已完成]　in:●ISM C++ 源碼　↳project_ism_method_soundness_validation
    - ✅ `T-V2` 距離矩陣分群驗證方法選單 [分析/已完成]　in:●7 guardrails…　↳project_distance_matrix_cluster_validation_methods
    - ☐ `T-V3` Fisher over-dispersion 必修（C++） [程式/待辦]　in:●7 guardrails…,●ISM C++ (Per…　↳project_code_methodology_audit_2026_06_10
    - ◐ `T-V4` ISM 取樣層 review（MAX_DIST vs SKIP 收斂） [程式/進行中]　in:●paired BAM　缺:U1-U5 不可寫定論；待 ISM 加取樣 provenance 輸出(cpp-change)　↳project_ism_sampling_layer_review
    - ☐ `T-V5` ISM vs 外部甲基法比較 Phase B [分析/待辦]　in:●ISM 方法定位,●外部工具 (modkit…　缺:Phase B 待核准　↳project_ism_vs_external_methylation_tools_comparison
  - ◐ `T-NEG` filter-NEGATIVE 四道脊柱 [里程碑/進行中]
    - ✅ `T-NEG-D1` D1 乾淨 somatic-cis 稀有 [分析/已完成]　↳project_paper_claim_audit_consensus_base_2026_06_12
    - ✅ `T-NEG-D2` D2 甲基=germline-allelic 多重共線 [分析/已完成]　↳project_methyl_phasing_assist_line
    - ✅ `T-NEG-D3` D3 甲基分群=germline 非 somatic 驅動 [分析/已完成]　↳project_tumor_only_axis_negative_subclone_classification
    - ✅ `T-NEG-D4` D4 NGroups=phasing 非甲基訊號 [分析/已完成]　↳project_hpfinengroups_subclone_marker
  - ◐ `T-PHASE` 甲基-assist phasing [里程碑/進行中]
    - ✅ `T-P1` 甲基救 unphase / haplotag assist V1-V12 [程式/已完成]　in:●tagged BAM,●longphase-S …　↳project_methyl_phasing_assist_line
    - ☐ `T-P2` 修正 T2 OVERSTATED 口徑（1-1/2-1 可分歸 H3 未證） [撰寫/待辦]　in:●V1-V12 結果　↳project_methyl_phasing_assist_line
  - ◐ `T-WRITE` 論文撰寫 [里程碑/進行中]
    - ☐ `T-W-ABS` 摘要（首段點明分工） [撰寫/待辦]　缺:🔴 須點明分工；待 G-A 數字　↳project_thesis_writing_architecture
    - ☐ `T-W-CH1` Ch1 緒論 [撰寫/待辦]　↳project_thesis_writing_architecture
    - ☐ `T-W-CH2` Ch2 文獻探討 [撰寫/待辦]　↳project_ism_vs_external_methylation_tools_comparison
    - ◐ `T-W-CH3` Ch3 材料與方法（骨架已有） [撰寫/進行中]　↳project_thesis_writing_architecture
    - ◐ `T-W-CH4` Ch4 結果（骨架已有/填數中） [撰寫/進行中]　缺:R6 待 G-A；數字待真值　↳project_thesis_writing_architecture
    - ☐ `T-W-CH5` Ch5 討論 [撰寫/待辦]　↳project_g6_paper_framing_external_corroboration
    - ☐ `T-W-CH6` Ch6 結論 [撰寫/待辦]
    - ☐ `T-W-FIG1` Fig1 方法總覽 schematic ⭐最高優先 [撰寫/待辦]　↳project_thesis_writing_architecture
    - ☐ `T-W-FIGS` Fig2-6 其餘圖（含跨樣本 R6） [撰寫/待辦]　缺:5/6 Fig 物理不存在
    - ☐ `T-W-TABS` Table1-3（樣本統計/工具對照/NEGATIVE 摘要） [撰寫/待辦]　缺:Table1/3 缺
    - ☐ `T-W3` 整合篇章（ASM-char + 4 道 NEGATIVE + LOH 脊柱） [撰寫/待辦]　in:●COLO829 pair…,●Stage③ 結果　缺:待上游 COLO829 / Stage③ / 外部真值補齊　↳project_subclonal_reconstruction_paper_focus
    - ✅ `T-W4` 論文敘述對抗稽核共識底座（51 agents / 606 claim / F=0） [撰寫/已完成]　in:●論文 spec,●606 claim 集　缺:Hard-Gate 待修正 ledger:95 / CURRENT_FOCUS:137　↳project_paper_claim_audit_consensus_base_2026_06_12

## Ready — 可立即開跑（19）
- `T-A3` 6 normal 甲基帳號 SPOF 備份　owner: user
- `T-C-CROSS` 跨樣本 ASM 正式圖表（R6）　owner: claude
- `T-E3` citation-verification（投稿前）　owner: claude
- `T-GATE-GB` G-B：within-hap somatic null → 定甲基-subclone 故事　owner: claude+user
- `T-GATE-HD1` HD-1：R-SELFREF → 定論文 Grade　owner: claude+user
- `T-P2` 修正 T2 OVERSTATED 口徑（1-1/2-1 可分歸 H3 未證）　owner: claude
- `T-SL-C1-F` 修：enable_bootstrap=true + within-HP 用 bootstrap　owner: claude
- `T-SL-C3` ④ C3 HP-fine 組合對齊測試 emit　owner: claude
- `T-SL-RO-1` P-readset：--normal-bam 是否 optional → tumor-only rerun　owner: claude
- `T-SL-S3` Stage③ coverage-peak CN 自估 + SEQC2 驗證腳本　owner: claude
- `T-V3` Fisher over-dispersion 必修（C++）　owner: claude
- `T-V5` ISM vs 外部甲基法比較 Phase B　owner: claude
- `T-W-ABS` 摘要（首段點明分工）　owner: claude+user
- `T-W-CH1` Ch1 緒論　owner: claude+user
- `T-W-CH2` Ch2 文獻探討　owner: claude+user
- `T-W-CH5` Ch5 討論　owner: claude+user
- `T-W-CH6` Ch6 結論　owner: claude+user
- `T-W-FIG1` Fig1 方法總覽 schematic ⭐最高優先　owner: claude+user
- `T-W-TABS` Table1-3（樣本統計/工具對照/NEGATIVE 摘要）　owner: claude+user

## 細節不清楚 / 缺資訊
- ✅ 無結構性細節缺口
**各任務 missing_info：**
- `T-SL`: 救回主軸 allele≫HP 需 cis-control 分 cis-ASM vs subclone
- `T-SL-C1-F`: C1 spec 20260618_within_hp_bootstrap_and_cluster_label_ari_cpp_change_spec_01.md pending_approval
- `T-SL-C1-V`: 待 C1-F 落地後跑
- `T-SL-C3`: C3 spec pending；落地序 C1→C3
- `T-SL-RO-1`: 先驗 --normal-bam 是否 optional（可能零改碼）
- `T-SL-RO-2`: 待 P-readset
- `T-SL-S3`: 🔴 SEQC2 不可當分類輸入（用戶定）；原載 SEQC2 BED 當輸入作廢重設計
- `T-CSC`: ISM 完成才開跑；六層框架 ⑥ subclone 下游
- `T-A2`: COLO829 normal 甲基缺；6 normal 全 zhenyu112 帳號 = SPOF
- `T-M1`: COLO829 normal 甲基缺 → 無法全 6 樣本
- `T-C6`: COLO829 normal 甲基缺
- `T-S1`: B 排序 illustrative / somatic 未定待 G-B
- `T-S5`: A 路相關感知 null 真實資料失敗 81% → 需 model-based 大改
- `T-V4`: U1-U5 不可寫定論；待 ISM 加取樣 provenance 輸出(cpp-change)
- `T-V5`: Phase B 待核准
- `T-E2`: H1437/H2009/HCC1937/HCC1954 缺專門解答（需自證）
- `T-A3`: 6 normal 全 zhenyu112 帳號 = 單點失效
- `T-C-CROSS`: 待 G-A 跨樣本統計
- `T-GATE-HD1`: 未決；影響論文主張強度
- `T-GATE-GA`: 待 COLO829 normal 甲基；G-A 統計未跑
- `T-GATE-GB`: 未跑前甲基-subclone 只能寫存在性窄+負
- `T-W-CH4`: R6 待 G-A；數字待真值
- `T-W-ABS`: 🔴 須點明分工；待 G-A 數字
- `T-W-FIGS`: 5/6 Fig 物理不存在
- `T-W-TABS`: Table1/3 缺
- `T-W3`: 待上游 COLO829 / Stage③ / 外部真值補齊
- `T-W4`: Hard-Gate 待修正 ledger:95 / CURRENT_FOCUS:137

---
單一所有權：本層只擁含括/依賴/分派/狀態/缺漏/io；verdict・證據・敘述歸 cycle/ledger/memory。
