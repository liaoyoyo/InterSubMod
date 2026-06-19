# 研究任務樹（WBS） — TASKS_ism.md

> **自動產生，請勿手改。** 唯一真值 = `state/tasks/graph_ism.json`；改後重跑 `python3 scripts/tasks/task_graph.py --graph state/tasks/graph_ism.json render`。
> 生成 2026-06-20 · build `78339d2` @ `research/subclonal-reconstruction-202606` · 節點 26（✅13 ◐4 ⛔0 ☐9 ✗0） · focus = `T-SL`

## ⓪ 驗證（validate + check）
- validate: ✅ PASS
- check: ✅ 0 findings

## 聚焦路徑　論文（碩論：Subclonal reconstruction） ＞ 完成 ISM 程式 ＞ 修正「結構驅動切區塊 × 標籤驗證」

## 任務樹（含括 parent；縮排=層級）
- ◐ `T-PAPER` 論文（碩論：Subclonal reconstruction） [里程碑/進行中]　↳project_thesis_writing_architecture
  - ☐ `T-CSC` 驗證判讀確認 → clone / subclone [分析/待辦]　缺:ISM 完成才開跑；六層框架 ⑥ subclone 下游　↳project_apriori_subclone_classification_model
  - ◐ `T-ISM` 完成 ISM 程式 [里程碑/進行中]　↳project_ism_method_soundness_validation
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
  - ☐ `T-WRITE` 整理論文 [撰寫/待辦]　缺:🔴 摘要首段須點明分工（reconstruction骨幹/甲基characterize有界）　↳project_thesis_writing_architecture

## Ready — 可立即開跑（4）
- `T-SL-C1-F` 修：enable_bootstrap=true + within-HP 用 bootstrap　owner: claude
- `T-SL-C3` ④ C3 HP-fine 組合對齊測試 emit　owner: claude
- `T-SL-RO-1` P-readset：--normal-bam 是否 optional → tumor-only rerun　owner: claude
- `T-SL-S3` Stage③ coverage-peak CN 自估 + SEQC2 驗證腳本　owner: claude

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
- `T-WRITE`: 🔴 摘要首段須點明分工（reconstruction骨幹/甲基characterize有界）

---
單一所有權：本層只擁含括/依賴/分派/狀態/缺漏/io；verdict・證據・敘述歸 cycle/ledger/memory。
