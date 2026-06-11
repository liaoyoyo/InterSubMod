<!--
建立時間: 2026-06-09
狀態: handoff (研究與觀察任務清單 — 給其他 AI session 獨立領取執行)
報告類型: paper_focus_session_handoff_tasklist
受眾: 其他 AI session（執行者）· 廖子游
provenance_note: 任務由 03 方向卡裁決 + 06-09 決策 + 05 任務樹導出；數字沿用已驗證集合。
-->
<!-- provenance-verified: 任務由 docs/paper_focus 既有裁決導出；本檔為 handoff spec，非新分析。 -->

# 研究與觀察任務清單 — 給其他 session 執行（G6 論文）

> **L0 一眼結論**：9 個可獨立分派的任務（CODE / OBSERVATION / VALIDATE / PROVENANCE），每個**自包含**（目標/輸入/步驟/輸出/驗收/依賴/caveat），其他 session 直接領。**最該先做＝T-CODE-1 catalog 骨架**（不卡 HD-1、高 ROI）。
>
> **L1 執行鐵則（每個任務都套）**：
> ① **守 characterization≠filter**：catalog 是描述目錄，TAG-D（FP-prone）只標註不當 filter。
> ② **§13 序列**：分析→落檔(.json/.tsv)→Read 讀回驗證非 error→**才**寫報告；產數字與寫報告**不同 batch**。
> ③ **§8 長 compute**：用背景 `Bash(run_in_background)` 跑，**不放 workflow agent step**；≥3 樣本實跑用 parallel-benchmark 但先 resource preflight。
> ④ **Scope（task-type B/C）**：catalog / cross-sample 預設**全 6 樣本 + 全基因組**；single-pipeline ClairS-TO → tier 封頂 ⭐3。
> ⑤ **數字帶 tier**：🟢 原檔對賬 / 🟡 待對賬（投稿前補）。
>
> **權威輸入**：方向裁決 `01_focus_notes/03` · 驗證機制 `04` · catalog schema `02_paper_framework/位點甲基分群catalog_schema提案.md` · provenance `02` · session 任務 #15–#20。

---

## T-CODE-1 ｜位點甲基分群 catalog 骨架（ISM-4，⭐ 先做，session 任務 #15）

- **目標**：依定稿 schema（16 欄+7 TAG）組「一位點一列」catalog；先跑**骨架**（現有欄），Δβ/cis 之後補。
- **輸入**：6-sample ISM native run（ledger `20260604_ISM_complete_TPFPFN`，294723 TP/3177 FP/34805 FN）+ `genome_survey_v2/` + schema 提案檔。
- **步驟**：(1) 抽現有欄 `clustering_reliability`(CramersV+PERMANOVA)/`clustering_origin`(ARI ruler)/`somatic_status`/`axis`/`LOH_status`/`coverage`/`n_CpG`/`gene` per locus → (2) 套 7-TAG 規則 → (3) 輸出 `catalog_skeleton.tsv` + per-TAG 計數。
- **輸出**：`research/.../catalog_skeleton.tsv`（每欄帶實際數值）+ `catalog_tag_counts.json`。
- **驗收**：每位點有 1 個 TAG；TAG-A（乾淨 cis）能列出（chr17 應在內）；per-TAG 計數合理（TAG-C/G 佔多數）。
- **依賴**：閾值校準（用 imprinting 正控 GNAS/RB1=1.000 + chr17 校 PERMANOVA p）。**不卡 HD-1**。
- **caveat**：Δβ 欄（依 T-CODE-3）+ cis 欄（依 T-VAL-1）骨架階段先留空；§13 落檔→Read 驗。

## T-CODE-2 ｜CramersV reliability gate 過嚴審查（ISM-2，#16）

- **目標**：審 CramersV gate（`CramersV=reliable?v:0`，Cochran 最小期望格≥5，RegionProcessor.cpp:1592）是否過嚴。
- **輸入**：762 MISSED（487 latent = CV gated=0 但 PERMANOVA 顯著）。
- **步驟**：/methodology-audit → 量 487 latent 的 TP/FP 組成 + minN 分佈 → 評估二元 0/1 改三態（reliable/latent/unreliable）。
- **輸出**：methodology-audit 決策文件（不直接改碼）。
- **驗收**：量化「過嚴」影響 + 給三態提案 + **明確：納回只提升 characterization，不可當 filter**（minN<15 FP-leaning）。
- **依賴**：無。**🔴 若決定改 C++ → 走 /cpp-change + 編譯 Hard Gate**。
- **caveat**：不可放寬成 filter。

## T-CODE-3 ｜ISM 合併 Δβ 互補發現通道（ISM-3，#17）

- **目標**：ISM 輸出加 per-position `dbeta_max` 欄 + union 發現標記（discrete-clustered **OR** large-Δβ-且有-clustering）。
- **輸入**：ISM read×CpG + genome_survey Δβ。
- **步驟**：算 per-position allele Δβ → 加欄 → union 標記（**不 naive 加**：高-Δβ regime FP-enriched）→ 輸出兩通道 Venn。
- **輸出**：catalog `dbeta_max` 欄 + `discovery_venn.tsv`（ISM-only / Δβ-only / 交集 + 各 TP/FP）。
- **驗收**：BRCA2（clustering-only Δβ~0.007）在 ISM-only；高-Δβ FP 在 Δβ-only 且被標 TAG-D。
- **依賴**：feeds T-CODE-1 的 `dbeta_magnitude` 欄。**🔴 涉 C++/Python 改動。**
- **caveat**：守 characterization；Δβ 不取代 CramersV。

## T-OBS-1 ｜catalog 6 張觀察分佈圖（P1–P6，見樹也見林）

- **目標**：產 schema C 節的 6 張分佈圖（forest view）。
- **輸入**：catalog_skeleton.tsv（+ Δβ/cis 補完後）。
- **步驟**：P1 Δβ 分佈 / P2 CramersV 分佈(reliable 閾值+latent 區) / P3 ARI 分佈(somatic vs het-null vs imprinting) / P4 within_d vs HP_d scatter / **P5 per-TAG 計數長條(=Results R6 主圖)** / P6 coverage vs reliability。
- **輸出**：`04_figures/P{1-6}_*.png`（matplotlib + `scripts/lib/plot_setup.py` CJK 字型）。
- **驗收**：P5 能當 R6 主圖；P6 確認 reliability 非純被覆蓋驅動。
- **依賴**：T-CODE-1（P2/P3/P6 可先跑）；P1/P4 依 T-CODE-3/T-VAL-1。

## T-OBS-2 ｜每 TAG canonical 例 + chr17 第二樣本重現

- **目標**：每個 TAG 挑 1 個 canonical 位點做 read×CpG 視覺化（見樹）；chr17/TBC1D16 在第二樣本重現觀察。
- **輸入**：catalog + ISM display_v2（read×CpG + ISM 距離矩陣，需 figs/ 重生）。
- **步驟**：per-TAG 抽代表位點畫 read×CpG heatmap + ISM 距離矩陣；chr17 在另一乳腺樣本(HCC1937/1954)查是否同向。
- **輸出**：`04_figures/TAG_{A-G}_example_*.png` + chr17 第二樣本 .md。
- **驗收**：TAG-A 例（chr17）清楚；chr17 第二樣本有/無重現都誠實記。
- **caveat**：chr17 單 locus 單樣本 nominal p → 第二樣本是升信心關鍵。

## T-VAL-1 ｜cis-candidate 大規模驗證協議（D3，#18）

- **目標**：設計+跑可規模化 normal-anchored cis-test，補 catalog `cis_status` 欄。
- **輸入**：816 HP-axis loci + normal BAM（`HCC1395BL_...5mCG_5hmCG_tagged.bam`）+ SEQC2 CN。
- **步驟**：normal HP1/tumor HP1/tumor HP1-1 三路 cis-test + copy-partition + **anchor SEQC2 CN** + 對 pure-ALT/CGI-desert untestable 標註。
- **輸出**：`cis_status` 欄 + `cis_protocol_results.json`。
- **驗收**：copy-partition 接真值；untestable（6 個）明確標；cis（chr17）保留。
- **依賴**：COLO829（T-VAL-3）解 untestable。feeds T-CODE-1。
- **caveat**：🟡 chr17 perm p / BRCA2 ~80%copy 須對賬原檔（見 T-PROV）。

## T-VAL-2 ｜HD-4 closeout（AF→NGroups partial-corr knockout，N9，#19）

- **目標**：ISM re-run 匯出 `fine_group_counts{HP1,HP1-1,HP2,HP2-1}/variant` + `n_clusters`(甲基-only count) 做乾淨負控，封閉 HD-4=phasing。
- **輸入**：ISM + HCC1395 baseline（n=36854）。
- **步驟**：re-run 匯這兩欄 → partial-corr knockout（控 n_clusters 後 AF→NGroups r 是否仍在）。
- **輸出**：`hd4_knockout.json`（補刀 definitional+mechanism+分布三證）。
- **驗收**：控 n_clusters 後 r 仍高 → 確認 phasing 非甲基。**off 關鍵路徑（HD-4 已 RESOLVED，此為精修）**。

## T-VAL-3 ｜umtag yardstick + COLO829 ⭐4（D1/N8/N10，#19）

- **目標**：(a) COLO829 升 ⭐4（解 ASM 天花板 + 給 phase-truth）；(b) umtag switch-error/N50 vs MethPhaser/HapBridge + 真-unphase 救援實測。
- **輸入**：COLO829 BAM + selfphasing_v6 + methyl-phasing pilot。
- **步驟**：COLO829 跑 ISM ASM → ⭐4；umtag 算 switch-error/N50 + 真救援（非 held-out）+ V10 非-copy 在 COLO829 重現。
- **輸出**：COLO829 ASM results + umtag yardstick.json。
- **驗收**：⭐4 解鎖；umtag 有 yardstick 才可寫 Discussion tooling。
- **caveat**：🟡 held-out ≠ 真救援；單樣本→COLO829 是關鍵。**off 關鍵路徑（資源到位再做）**。

## T-PROV ｜投稿前 provenance 對賬（~10 個 🟡 數字，#20）

- **目標**：把 `02 文件庫` 🟡 清單定位正確源檔再 grep，升 P-level。
- **清單**：condition_fp_consensus.json 正確檔名→OR 8.63/4.09/5.84/0.194 ｜cycle2 transfer→−0.377 ｜上游 TSV→null ARI 0.177 ｜tsg 另一檔→chr17 perm p+BRCA2 ~80%copy ｜*_gwasm.json→6/6 excess ｜補 E[overlap]=0.16 ｜same-hap 93% 源 ｜umtag 數字 ｜/citation-verification 核外部 PMID/DOI。
- **驗收**：每數字 grep 得到 → 升 🟢；找不到 → 標 retract/restate。
- **依賴**：tsg 專案 06-07 仍在寫 → 以其定稿為準。**投稿前必做**。

---

## L1 — 建議分派序（不卡 HD-1）
1. **T-CODE-1 catalog 骨架**（最高 ROI）+ **T-CODE-2 CramersV audit**（並行，純 audit）
2. **T-CODE-3 merge Δβ** → 補 catalog dbeta 欄 → **T-OBS-1 分佈圖**（P5=R6 主圖）+ **T-OBS-2 例**
3. **T-PROV provenance 對賬**（投稿前）
4. **T-VAL-1 cis 協議** + **T-VAL-3 COLO829**（資源到位）→ **T-VAL-2 HD-4 closeout**（精修）

> 🔴 **HD-1（phasing R-SELFREF）獨立、用戶暫緩** → phasing 維持 Grade B+，**這批任務都不依賴 HD-1**。
