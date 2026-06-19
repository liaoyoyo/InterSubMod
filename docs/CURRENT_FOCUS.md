<!--
建立時間: 2026-01-12 00:00
更新時間: 2026-06-11 (🎯 主軸轉向：Subclonal reconstruction 取代 G6；G6/G1 降支撐；見固定焦點區塊 + foundation doc) | 2026-06-02 (固定焦點區塊 + G6/G1 cycle); 5/30: SoT reconcile
狀態: validated
資料來源:
  - docs/standards/20260228_文件命名與狀態管理規範_01.md
  - docs/standards/20260228_output軟連結與版本控管規範_01.md
  - scripts/analysis/check_ai_agent_readiness.sh
  - docs/reports/validated/2026/04/20260423_研究週報_20260416_20260423_NG2_LOH_constrained_phasing與TO_pivot_01.md (latest)
  - docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md
-->

# 當前目標

## 2026-06-19 — 📚 外部文獻驗證庫就緒（69 源親讀，稽核 CLEAN）⭐ 最新

> **新 session / 論文撰寫先讀**：外部論文/程式碼實體驗證庫已建於 **repo 外** `/big7_disk/liaoyoyo2001/external_validation/`（**69 源** CONTEXT 卡 + 39 repo + 36 PDF；tier A×13·B×34·C×22；軸 18/12/11/14/8/6；61/69 親讀 + 8 cited_secondary〔mcf7/deamination + 6 why-hard L3-abstract〕；2026-06-19 +8 why-hard axis1 問題陳述卡=SNV-only/CNV 困難一手源，doc `docs/method_comparison/20260619_subclone_analysis_interpretation_full_framework_01.md`）。**整體準確性對抗稽核（20 卡）= CLEAN**（0 MAJOR）；2026-06-15 5 cis-basis 卡親讀 PDF 升 fulltext_verified（並修出 Onuchic「NPD」捏造方法名 / Min heritability 二手數字 / Turcan citation 年份錯三處瑕疵）。repo 內橋接索引 = `InterSubMod/docs/method_comparison/20260613_external_validation_library_index_01.md`（含論文 Ch2 逐軸對應 + paper-critical 結論）；整體論文文獻地圖 = `external_validation/_landscape/08_paper_literature_map.md`。memory: `project_external_validation_library`。
>
> **3 個投稿必守口徑**：① 與 69 源 **0 真 CONFLICT**（最尖對比物皆 regime 差，EVOFLUx 自承 subclone 罕見反佐證弱-subclone）；② **ISM 創新點口徑 = 無監督 read×read 距離矩陣結構 PERMANOVA + normal-baseline cis-test + somatic-subclone 目標**（🔴**禁**用「對手二代定序」或「對手缺顯著性檢定」當差異 — cvlr/ASMS/MethylBERT 都 ONT-capable 且都有 randomization 檢定，會被 reviewer 打臉）；③ cancer 甲基-phasing 白地 source+全文雙證，**LongHap 2026=germline-only 不威脅**（Fig3C 強不對稱甲基反佐證 R2）。

---

## 2026-06-11 — 🎯 主軸轉向（用戶確認）：Subclonal reconstruction 取代 G6 ⭐

> **新 session 先讀**：`InterSubMod/docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md`（新主軸基礎：6 樣本×3 癌種資產盤點 + V1-V12 成果 map + 誠實 gap + 交接清單）。
>
> **決策**：新論文主軸 = **Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing**，**取代 G6 phasing**；G6 LOH-phasing / G1 ASM **降為支撐材料 park**。本 session 已「整理放緩」所有現有任務（甲基救 unphase/tag 矯正線 V1-V12 SoT 完整），論文細節由後續 session 完成（見 foundation §6 交接清單）。
>
> **資產（緩解單樣本限制）**：6 cell line × 3 癌種（HCC1395/1937/1954 乳腺 · H1437/H2009 肺 · COLO829 黑色素瘤）皆有 somatic-haplotag BAM + somatic VCF + 甲基 → 可衝 ⭐4。
>
> **誠實天花板（對外必守）**：甲基=germline-haplotype 層級；T3 存在性窄翻、可用性 NEGATIVE；T2 只證 1-1/2-1 可分非歸 H3；非 variant filter(DEAD)。詳見 foundation §4。
>
> ⬇ 以下 6/08 段為轉向前的 G6 主軸紀錄，保留供 audit；G6/G1 成果併入新主軸為支撐。

---

## 2026-06-08 — 論文就緒收斂（跨 7 線獨立稽核）〔轉向前 G6 主軸，已降支撐〕

> **新 session 先讀**：`InterSubMod/docs/concepts/2026/06/20260608_研究現況地圖_整體目標與流程_給其他AI_01.md`（單頁現況地圖 + HTML）+ `InterSubMod/knowledge/11_external_literature/10_paper_readiness_convergence.md`（論文就緒收斂，最完整）。

**論文主軸定稿**：read-level LOH/haplotype + 甲基 **characterization + tooling** 論文（非 variant filter）。主體 = **三~四道防彈 NEGATIVE** + phasing 脊柱(**Grade B+ 非 A**) + ASM copy-confounded 支撐。**一篇能過 review 的論文今天就存在**，卡關只在 HD-1。

**🔴 別踩的雷（必讀）**：
1. 甲基→TP/FP filter 已 concluded NEGATIVE/**死四道** — 勿重開。
2. characterization-set（628/15391/strong-ASM）**絕不寫成 filter**（15391 TP/FP=1.0）。
3. **BRCA2 ≠ 乾淨 cis 錨點**（06-07 重分析 = **subclone/copy-confounded 主導**，HP1-1 是 longphase-S 的 somatic subclone tag〔germline-H1+somatic-ALT〕**非 copy、非 CN-dosage**〔dosage 已 REFUTED〕；focal cis 殘餘 d_within=−0.023 邊際·單樣本不可分離；**% split 不 robust**）；乾淨 cis 改 **chr17/TBC1D16**；hypo≠canonical hyper，勿當 TSG-silencing 證據。已 amend lit-07 + master_draft §2.2（06-09 R1/R2/R3 統一口徑落地）。
4. **「Grade A」是 Grade B+**（R-SELFREF circularity 對照未跑；p=0.0078 是方向一致性非 effect；n=7=bio-n=6）。
5. 所有 ASM 單樣本 HCC1395 ⭐3 / 單 pipeline；跨樣本 6/6 是 phenomenon-level；private 0/38 是 underpowered。

**開放決策**：🔴 **HD-1（用戶決定·hold）** phasing by-construction 循環 → 跑 R-SELFREF(~25-50hr C++) or 降為 characterization｜HD-3 ✅ BRCA2 已 amend｜HD-4 ✅ AF→NGroups=**phasing 非甲基**（NGroups=HP-tag count；ledger 20260608_HD4）→ T6 無正向甲基發現｜HD-5 ✅ T7 umtag 已註冊 ledger（FUTURE-WORK）｜HD-2/6 投稿前。

**並行**：tsg ASM/cis 線（另一 session，06-07 仍在寫，T2/T3 數字可能再動）。active cycle = G6 LOH-phasing (P4 ⭐3) + G1 ASM (P3 ⭐3)。

---

<!-- ┌─ 固定焦點區塊（2026-06-02 起；SessionStart 注入最先看到；每次主軸/背景變動時更新此區，不新增日期段；機械態 SoT = state/active.json+cycles，本區為人讀鏡像）─┐ -->

## ★ 當前焦點（pinned；2026-06-12 更新 — 資料驗證 + MEMORY 收斂；2026-06-11 主軸轉向）

> **🎯 新主軸（取代 G6）**：**Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing**。
> **兩個互補 SoT 面**：① 甲基-phasing-assist `InterSubMod/docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md`（V1-V12 + 6 樣本資產 + G-A~E gap）② ASM-characterization + 四道 NEGATIVE + LOH-phasing 脊柱 + HD-1 gate + cross-line reconcile `InterSubMod/docs/concepts/2026/06/20260611_Subclonal_Reconstruction_Paper_Focus_整合篇章_01.md`。
> **2026-06-11 `/pivot-direction` 已正式記錄轉向**（research_direction.md）。論文細節由後續 session 完成。

> **🟢 2026-06-12 資料驗證 + 整理收斂（本 session）**：
> - **G-A 資料 ready 升級**：6 樣本 tagged BAM+somatic VCF+ISM TP/FP 實測齊；**matched-normal 甲基實測 5/6 有**（HCC1395 5mC+5hmC · HCC1937/1954/H1437/H2009 5mC；**只 COLO829 缺**）→ V10 跨樣本對 **5 樣本（乳腺3+肺2）可直接跑衝 ⭐4**，G-A 非「全卡」。契約 `InterSubMod/docs/data_specs/20260612_external_data_dependencies_01.md`（6 normal 全 zhenyu112 帳號=SPOF）。
> - **數據準確度/可尋性 7 改進**已 commit（3380681/805d66a/069cadb）：provenance stamp · P-17(盤點別憑記憶) · redirect banner · harness 10 燈 · 外部依賴契約。稽核 `InterSubMod/docs/data_specs/20260612_data_accuracy_findability_improvement_audit_01.md`。
> - **MEMORY.md 收斂** 34KB→13.5KB、重組到新主軸（topic 檔全保留）。

| 主軸 | tier | 狀態 | 下一步 |
|------|------|------|--------|
| **🆕 Subclonal reconstruction（somatic haplotag + methylation）**（新論文主軸）| ⭐3→⭐4 候選 | 兩面 foundation 已立；6 樣本資產齊 + **normal 甲基 5/6 ready**；🔴 gate = **HD-1**（R-SELFREF 跑 or 降 characterization）+ **G-B**（subclone 甲基 somatic-specific vs germline-allelic）| ①定 HD-1 ②先寫四道 NEGATIVE methods（HD-1 獨立、今天防彈）③**G-A 跨 5 樣本重現**（資料已 ready）④G-B 對照 ⑤COLO829 補甲基 normal → 論文 outline |
| ~~G6 LOH-constrained phasing~~（**降為支撐材料**）| ⭐3 | ✅ park（=新主軸 phasing 脊柱，補充面詳述）| — |
| ~~G1 ZAR1L/BRCA2 ASM~~（**降為支撐材料**）| ⭐3 | ✅ park（=新主軸 ASM characterization 層；本 session 工作站證據基座）| — |

**甲基救 unphase / tag 矯正研究線（本 session 主體，已整理）**：V1-V12 SoT 完整（`20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md`）。誠實定論：甲基=germline-haplotype 層級；T1 unphase 救援 SUPPORTED+caveats(0.885，僅~6% unphase 可嘗試)；T2 OVERSTATED(只證1-1/2-1可分)；**T3 存在性窄翻案 + 可用性 NEGATIVE**（local-allele 亞群甲基可分 farCpG AUC 0.85，但救 ambiguous read lean<0.5）。memory `project_methyl_phasing_assist_line`。

> ❌ **DEAD（勿再開）**：甲基化當 FP filter（⭐2 L4）；T3 subclone「可用性」（救 ambiguous read 偏向反了）。
> 🔗 新主軸基礎 `docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md`｜真值 SoT `...A0_assets/VERIFIED_RESULTS.md`｜機械態 `/cycle-state`

<!-- └─ 固定焦點區塊結束 ─┘ -->

---

## 2026-05-31 — Harness 全方位 audit + 自我稽核儀表板 + 任務 reconcile ⭐

> **當前主軸（不變）**：兩條 LIVE 高價值方向 — **① LOH-constrained phasing ⭐3（Grade B+→A，論文主軸候選）** + **② ZAR1L/BRCA2 ASM ⭐3（characterization-only）**。甲基化當 FP filter ❌ DEAD（⭐2 L4，勿再開）。詳見下方 5/30 section。

**本 session（harness 工程，非研究）**：業界 AI 工作流研究 + 全 harness 盤點 + restraint 對抗驗證（11-agent workflow）→ **0 個新外部工具**，價值全在修 latent bug + drift。完整報告：`InterSubMod/docs/references/migration/20260531_harness_audit_dashboard_design_01.md`。

**已落地**：(1) 修 doc/agent drift（CLAUDE.md §4 Hard Gate 5→4、release co-author 4.8、research-orchestrator 路由、researcher MCP、cache_telemetry flat-rate 註解）；(2) agent 權限硬化（headless-research +isolation:worktree、reviewer 去 Write、literature-reviewer 標 deprecated、paper-miner Bash 收窄）；(3) **wire `pre_tier_upgrade_check` 成 Hard Gate**（state.json exit-2 / 散文 advisory）；(4) 新 `/harness-health` skill + `scripts/harness_health.py` 6-燈儀表板（45→46 skills）。

**任務 reconcile（解決儀表板紅燈 4&5）**：
- **cycle3 = 已 concluded**（active.json `recently_concluded` = Phase2 Cycle3 P6_COMMIT NEGATIVE @ 2026-05-30T07:40Z）→ 5/30 §D3「下一步跑 /conclude-research」**已完成**，非待辦。
- **5/31 ASM tier 語意收斂回寫**（6 條 ledger entry，CURRENT_FOCUS 先前未反映）：ZAR1L/BRCA2 ASM **real but non-directional + non-discriminative + coverage-modulated**（連續 |Δβ| AUC=0.505 NEUTRAL；strong-ASM 5× FP enrichment = LOH→single-hap→low-cov→extreme baseline 的 regression-to-extreme artifact，ONE regime）；B2 clustering = germline allelic (Layer A)，somatic NEGATIVE 經 imprinting 正控（NORMAL GNAS/RB1 ARI=1.0）validated 為真 biology。**對外引用 ASM 必用此收斂口徑（勿用「方向 POSITIVE」over-claim）**。tier ⭐3 數字不變，語意收窄。

**待用戶決定（本 session 新發現，未動）**：
- ✅ **[已解決，2026-06-12 fresh 驗證]** `pre_commit_compile_check` + `kb_schema_check` neuter **早已修復**（settings 現為 `1>&2`，2026-05-31 移除 `|| exit 0` mask）+ stale marker `/tmp/ism_cpp_pending_compile.txt` **已不存在**。harness_health 燈#2 HARD-GATE-TRUTH 持續監看 = GREEN。**勿再當 live 問題處理**（曾被就緒度稽核誤當 live 傳播）。
- 🟡 **stale queue H013-018**（filter/caller-f1 已 DEAD 仍 queued）→ 待 `/pivot-direction` 降權。

---

## 2026-05-30 — 研究狀態 audit（驗證後）+ SoT reconcile + ASM magnitude rerun 啟動 ⭐

> **當前主軸（2026-05-30 reconcile，用戶確認）**：兩條 LIVE 高價值方向 —
> **① LOH-constrained phasing ⭐3（Grade B+→A，論文主軸候選）** + **② ZAR1L/BRCA2 ASM ⭐3**。
> 此二者為**當前最有發表價值 + 最能成「說明報告」**的方向。
> 甲基化當 FP filter ❌ DEAD（⭐2 L4，direction exhausted，勿再開）。HKU handoff 已交付降級（見下方 [CLOSED]）。

**16-agent verified status audit 結論（見樹也見林）**：研究是「一條死路（FP filter）+ 兩條活路（LOH-phasing + ASM）」格局。

### 三大結論

| 主軸 | verdict | tier | 一句話 |
|------|---------|------|--------|
| **LOH-constrained phasing**（論文主軸候選）| ✅ LIVE / 最有發表潛力 | ⭐3 **Grade B+ 接近 A** | NG=2 Inner same-hap，TP gap **n=7 全 7/7 正向 Wilcoxon W=28 p=0.0078**（2026-05-30 加 COLO829 達成 Grade A 要求 #1）；剩 R-SELFREF 全 7-sample flag-on 負控（~25-50hr C++ 重跑）升 full A |
| **甲基化當 FP filter** | ❌ DEAD | ⭐2 L4 | Phase 2 LOSO 證實 +0.02236 = 100% circularity；mean −0.00004；methylation 5th-rank vestigial；filter direction EXHAUSTED |
| **ZAR1L/BRCA2 ASM** | ✅ 方向 POSITIVE + magnitude RESOLVED | ⭐3 (⭐4 需 COLO829) | BRCA2 HP-axis Δβ=**−0.122** / ALLELE −0.099（全 TP 重跑定案，與 script03 收斂）；全基因組 strong-ASM 172/51,171；**B-discrimination anti-discriminative**（strong-ASM 在 FP enriched 5×）= 真實但不能 filter |

### 本 session 動作（2026-05-30）

1. **SoT reconcile D1-D5**（修 5 處漂移）：
   - ✅ D1 MEMORY.md Phase2 Cycle1 ⭐3→⭐2 L4 DOWNGRADED
   - ✅ D2 ZAR1L magnitude 標 PROVISIONAL/superseded（memory + ledger correction entry 20260530）
   - ✅ D3 cycle3 **已 concluded**（2026-05-31 reconcile 確認：active.json `recently_concluded` = Phase2 Cycle3 P6_COMMIT `NEGATIVE_filter_direction_failed` @ 2026-05-30T07:40Z）。先前標「下一步跑 /conclude-research」已完成，非待辦。
   - ✅ D4 INDEX.md 補 5/21 PI signoff + 5/24 HKU + 5/29 ZAR1L + 修 header 日期
   - ✅ D5 本段（CURRENT_FOCUS refresh，解 6.1 天 staleness）
2. **ASM magnitude rerun — ✅ 完成定案（5/30）**：全 39,447 TP（51,171 records, 22 chr）max-collapse 修正口徑重跑完成。**BRCA2 HP-axis Δβ=−0.122 / ALLELE −0.099**（已驗證 v2 TSV，與 script03 收斂；buggy −0.054 = MSA Level1 5mC+5hmC 雙列砍半 artifact，修法在 Python `18_dual_axis_pivot.py` MAX-collapse，C++ 根因外部 repo 未動）。新發現 **B-discrimination anti-discriminative**（strong-ASM 在 FP enriched 5×, OR=0.194 p=1.8e-28）。MEMORY + INDEX + report banner 全對齊。

### 勿再開（concluded NEGATIVE guardrail）

pure-methylation TP/FP 判別 / TO germline-FP filter / TO QS post-hoc patch / ASM-as-discriminator / LR methylation-augmented filter — 全已 concluded dead，C1/C2/C3 productive-failure 條件未滿足前不 reopen。

---

## ✅ [CLOSED] 2026-05-22 ClairS-TO HKU Luo Lab 甲基 subclone handoff（per D4 已入 INDEX，5/24 產出）

> **2026-05-30 降級**：此 P0 已過期，且 per 本日 audit D4「INDEX.md 補 5/24 HKU」= 已交付 → **不再是 active 主軸**（當前主軸見頂部 5/30 section = LOH-phasing + ASM）。以下保留原始 plan 供 audit / 若需後續重啟。

**[原始 5/22 指派紀錄]**：5/24 結束前產出 standalone HTML 報告給香港 ClairS-TO Luo Lab，敘述 InterSubMod 既有「甲基 / phasing 相關 subclone-aware 功能」+ 列出兩 section 問題清單（對外 collaboration ask + 內部 audit list）。

### Scope + deliverable

- **目標受眾**：ClairS-TO Luo Lab (HKU) — caller 開發者
- **目的**：(1) 説明既有功能 (2) 評估能否幫他們提升 ClairS-TO caller F1 (3) 釐清需要對方協助/補資訊的 ask
- **格式**：standalone HTML (`html-report-build` standalone mode，sticky TOC + 折疊 cards + claim card tier badge)
- **路徑**：`InterSubMod/docs/reports/pi_reports/2026/05/20260524_ClairS_TO_HKU_methyl_subclone_handoff_01.{md,standalone.html}`

### 3 段時程 (D1 → D3)

| 日 | 階段 | 產出 |
|---|------|------|
| 5/22 (今天) | D1 PLAN + outline | `docs/plans/20260522_ClairS-TO_HKU_methyl_subclone_handoff_01_plan.md` (本 plan) + narrative-frame N1-N6 + 大綱 ack |
| 5/23 | D2 內容寫作 + evidence 引用驗證 | .md 完稿 (引用每條 claim 標 tier + commit hash) + 初版 standalone HTML |
| 5/24 | D3 self-review + 5 秒測試 + PI 終版 | scientific-rigor §2-§7 對齊 + design_principles 12 條 checklist + 最終 HTML + 兩 section 問題清單定稿 |

### 敘述地雷區（必避）

1. **HPFineNGroups ≠ methylation 訊號** — 是 `{HP1, HP1-1, HP2, HP2-1}` 4-bucket phasing × variant occupancy count（5/22 記憶 `feedback_feature_name_vs_definition_rule` 已警示）
2. **Phase 2 Cycle 1 +0.02236 是 in-distribution 結果** — 5/20 LOSO 證明 cross-sample ΔF1=-0.00004 (mean)，**不能對外宣稱 cross-sample 有效**
3. **paired_full LOH × AF × NGroups ΔNG=+0.787** — 是 paired mode（normal+tumor）結論，TO mode (tumor-only) 已 pivot 為 LOH-constrained phasing
4. **methylation 5 features 在 global LR 是 5th-rank** — 非 dominant signal；不可對 ClairS-TO 説「我們靠甲基」
5. **外部 caller F1 數字必查 KB** — known-pitfalls P-14 outside-claim 規則；不可從本專案 report 推論 ClairS-TO 行為

### 可 claim 的（已驗證 tier ≥ ⭐3）

- **LOH-constrained phasing discovery (TO mode)** ⭐3：NG=2 Inner ≥93% same-hap, TP gap +0.37, 6/6 sample 一致 (`project_loh_constrained_phasing_discovery`)
- **HPFineNGroups canonical filter** ⭐3 pipeline-dependent：NG=4+AF<0.4+NR≥80 NonLOH TP rate 92.8% (5/7 sample ≥0.85), master dataset+flag=off 條件
- **LOH × AF × Methylation Paired POSITIVE** ⭐3 (paired mode only)：ΔNG=+0.787, 7/7 p<1e-65
- **Phase 2 V3F→V5→V6 binary evolution**：V6 marker coverage +9.0% vs V3F, hp=33 ambiguous handling, caller F1 0.7166 不變

### Skill 啟動序列

- 立即：`/narrative-frame` (Tier 3 必跑 N1-N6 + 框架推薦 = Audience-Scenario-Pitch + Verdict-Pyramid 混合)
- D2 寫作前：`/scientific-rigor` §2-§7 (claim tier 標記 + L1-L5)
- D2 evidence 引用時：`/known-pitfalls` P-14 (外部 claim 必查 KB)
- D2 HTML 生成：`/html-report-build` standalone mode
- D3 self-review：`/scientific-rigor` §8.4 provenance footer + commit hash + 12 條 design_principles checklist

### Plan 文件

`InterSubMod/docs/plans/20260522_ClairS-TO_HKU_methyl_subclone_handoff_01_plan.md`

---

## 2026-05-20 (深夜) — H_NEW_2 FAIL + H_NEW_4 SANITY VIOLATED (Unexpected HCC1395 +0.00699) ⭐

**Session anchor**: 用戶 5/20 立刻啟動 H_NEW_2 + H_NEW_4 LOSO 驗證 observation-driven hypotheses → 兩個結果，一 FAIL 一 partial positive

### 3-Way LOSO Comparison

| Sample | Baseline LOSO (10f) | H_NEW_2 (2f loh+HPFineF) | H_NEW_4 (9f drop caller_af) |
|---|--:|--:|--:|
| HCC1395 | -0.00012 | -0.00012 | **+0.00699** ★ |
| HCC1937 | +0.00000 | +0.00000 | +0.00000 |
| HCC1954 | -0.00008 | -0.00000 | -0.00008 |
| H1437 | -0.00001 | -0.00000 | -0.00001 |
| H2009 | -0.00001 | +0.00000 | -0.00001 |
| **Mean** | -0.00004 | -0.00002 | **+0.00138** |

### Pre-Reg 結果對照（HARKing 防護）

| H | Prior | Observed | Match? |
|---|---:|---|---|
| H_NEW_2 ≥2/5 ΔF1 > +0.002 | 25% | 0/5 above threshold | ✅ consistent (conservative prior) |
| H_NEW_4 HCC1395 ≈ 0 | 80% | HCC1395 = +0.00699 | ❌ **VIOLATED** (post-hoc unexpected) |

### Mechanism 解讀

- caller_af 在 cross-sample LR train 中是 **confusing signal** (HCC1395 d=+1.60 vs HCC1937 d=-1.41)
- Drop caller_af 後 LR 純粹用 LOH + Coverage + NG + 5 methyl + chr8 train → 4 sample 的 weak coherent signal 對 HCC1395 marginal effective
- 其他 4 樣本仍 caller-F1-ceiling 卡住 (best τ=0.10 keep all)

### Verdict + Cycle 5+ Path

| Confirmed | Unknown |
|---|---|
| caller_af direction-inconsistent 是 LOSO 災難主因 | HCC1395 +0.00699 是 single-sample artifact 還是真 cross-sample? |
| 4 高 caller F1 sample 受 ceiling 限制 | 其他 low-F1 sample 在 H_NEW_4 設定下表現? |
| 2-feature (loh + HPFineF) LR 不夠強 | Per-zone / RF / interaction terms 是否能改善? |
| LR 在 universal production 設定下不 deployable | Non-LR framework (zone rules / boolean filter) 是否更適合? |

**建議 Cycle 5 路徑** (按 ROI):
1. 🟡 **Path A (推薦)**: Pivot phase_block_3d / thread_d — LR direction 已 exhausted, ROI 偏低
2. 🟡 Path B: H_NEW_3 chr8-specific zone gate — 不依賴 LR threshold
3. 🟡 Path C: 找 low-F1 panel 驗 HCC1395 +0.00699 generalize — 仍 caller_af-bounded

### Artifacts (2026-05-20 深夜)

- Findings: `InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/loso_hnew_findings.md`
- Scripts: `InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/scripts/run_loso_hnew{2,4}.py` + `loso_3way_comparison.py`
- Data: `InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/data/loso_hnew{2,4}_results.tsv`
- Figures: `InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/figures/loso_3way_comparison.png` + `loso_hnew{2,4}_5sample.png`
- Ledger: `cycle_id: 20260520_loso_hnew2_hnew4_observation_validation` (line 51)

---

## 2026-05-20 (晚) — 🔴 LOSO Sample-Level CV: Cycle 1 +0.02236 ≡ 100% Sample-Level Circularity Bias ⭐⭐

**Session anchor**: 用戶 2026-05-20 直接質疑「LR filter 用 HCC1395 數據訓練,又用 HCC1395 數據驗證,這敘述合理嗎」 → 觸發 LOSO 真正 sample-level cross-validation → **Phase 2 結論需 major reframe**

### LOSO 結果 (5 sample, 30 sec wall clock, reuse cycle 2 functions)

| Test sample | LOSO ΔF1 | Best τ |
|---|--:|--:|
| HCC1395 | **-0.00012** (vs in-dist +0.02236) | 0.10 |
| HCC1937 | +0.00000 | 0.10 |
| HCC1954 | -0.00008 | 0.10 |
| H1437 | -0.00001 | 0.10 |
| H2009 | -0.00001 | 0.10 |
| **Mean** | **-0.00004** | — |

- HCC1395 in-distribution → LOSO drop = **+0.02248 (100% effect from sample-level circularity)**
- 5/5 best τ → 0.10 (=keep everything) = LR 對任何 held-out sample 找不到 useful threshold
- Wilcoxon p=0.125, **0/5 pos, 4/5 neg → DIRECTION_NEGATIVE**

### Tier 重評（核心 reframe）

| Claim | Before | After LOSO |
|---|---|---|
| HCC1395 +0.02236 | ⭐⭐⭐⭐ L2 | **⭐⭐ L4** (HCC1395 in-distribution case study only) |
| BAM-invariant V3F/V5/V6 | ⭐⭐⭐⭐ L2 | ⭐⭐⭐⭐ L2 (unchanged) |
| Cross-sample transfer NEGATIVE | ⭐⭐⭐ L3 | **⭐⭐⭐⭐ L2** (LOSO 強化) |
| ISM vestigial in LR | ⭐⭐⭐⭐ L2 | ⭐⭐⭐ L3 (scope 變窄) |
| Caller-F1-headroom mechanism | ⭐⭐⭐ L3 | **⭐⭐ L4** (filter 整體失效 → mechanism moot) |
| **LOSO sample-level negative** | (new) | **⭐⭐⭐⭐ L2** (new core claim) |

### Phase 2 Decision Reframe

| 動作 | 狀態 |
|---|---|
| Cycle 1 +0.02236 ⭐⭐⭐⭐ → ⭐⭐ | ✅ 落地 ledger + memory + INDEX + PI Trust HTML |
| Cycle 4 Trial A/B/C prior PASS 顯著下降 (75/30/35% → 15/10/10%) | 建議 deprioritize |
| Paper §3 完全撤回 "ISM-augmented filter" 宣稱 | 待 paper draft 啟動 |
| 改為 "ISM characterization study + LR sample-level negative finding" | framing 已 reframe |
| **Pivot phase_block_3d 或 thread_d** | 🟠 等用戶確認方向 |
| ISM characterization (v0.3 cycle ⭐3 不變) | ✅ 保留 (LOSO 不影響 mechanistic understanding) |

### LOSO Artifacts

- Findings: `InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/loso_findings.md`
- Script: `InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/scripts/run_loso_cv.py`
- Data TSV: `InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/data/loso_cv_results.tsv`
- Figure: `InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/figures/loso_5sample_dF1.png`
- Updated PI Trust HTML: `InterSubMod/research/methyl_augmented_filter_phase2/phase2_pi_verification/phase2_pi_trust_framework.standalone.html` (Section 1.5)
- Ledger entry: `cycle_id: 20260520_loso_sample_level_circularity_revealed` (line 50)

---

## 2026-05-20 — Cycle 3 Step 1.5 Ablation: ISM Vestigial Confirmed + Filter Reframe ⭐

**Session anchor**：5/19 cycle 3 Step 1 PASS (qualifying mean +0.01499) 但 74.6% by HCC1395 alone → 用戶質疑 ISM 是否真貢獻 → plan v2.1 pre-reg 4 H ablation → 最小評估 (~25 sec on cycle 2 既有資料) → **H_M1a FAIL + H_A1 MARGINAL**

### Ablation Verdicts (n=5, refit OOF + transfer coef shrinkage)

| H | Pre-reg threshold | Computed | Verdict | 含義 |
|---|---|---|---|---|
| **H_M1a** drop ISM → HCC1395 refit ΔF1 drop | PASS ≥ +0.003 / FAIL < +0.001 | **+0.00065** | **FAIL** | ISM 為 vestigial covariate |
| **H_A1** caller_af shrink → HCC1954 ΔF1 改善 | PASS > +0.30 / FAIL ≤ +0.10 | **+0.25114** | **MARGINAL** | caller_af = 67% disaster confound |

**5-sample refit pivot 關鍵發現**：
- full mean +0.00619 vs **no-methyl mean +0.00630** (no-methyl 略勝)
- 4/5 樣本 no-methyl ΔF1 ≥ full (HCC1395 -0.00065 是唯一 marginal drop)
- HPFineF coef +0.75 (cycle 1 LR rank 5) 為 caller_af L2 ridge-split 而非 incremental signal
- 與 cycle 1 主報告 Step 5c「methylation 訊號實為 caller_af proxy」線索完全一致

### Cycle 3 Reframe Action（5/20）

| 動作 | 狀態 |
|---|---|
| Filter 命名: "methyl-augmented" → **"caller-F1-headroom-gated 4-feature filter"** (caller_af + LOH_inner + Coverage_Multiple + NG) | ✅ 落地 |
| `cycle3_caller_f1_gate.json` features_production = 4 | ✅ 更新 |
| Paper §3 撤回 "methylation-augmented filter" 宣稱 | 待 paper draft |
| ISM characterization (v0.3 cycle ⭐3) 保留作 mechanistic understanding | ✅ 不變 |
| Cycle 3 Step 2 panel survey 改驗 4-feature filter cross-sample | 待啟動 |
| 大規模 ablation (M2/M3/M0/LOFO) 降級 paper supplementary | ✅ deferred |
| A2 caller_af KS test (~30 min) 補強 H_A1 殘餘 33% | 建議跑 |

### 跨 session 結論交叉確認

5/19 V6 vs baseline 報告 §13 Day 3 5-Goal validation 已測 H_ABL_1 methylation contribution +0.0005-0.0007 — **與本 cycle 3 Step 1.5 H_M1a FAIL (+0.00065) 同方向同 magnitude** — 兩條獨立路徑收斂同結論強化信度。

### Artifacts

- Findings: `InterSubMod/research/methyl_augmented_filter_phase2/cycle3/ablation/cycle3_step1_5_min_ablation_findings.md`
- Script: `InterSubMod/research/methyl_augmented_filter_phase2/cycle3/ablation/scripts/run_ablation_variants.py`
- Data TSV: `InterSubMod/research/methyl_augmented_filter_phase2/cycle3/ablation/data/cycle3_step1_5_min_ablation.tsv`
- Figure: `InterSubMod/research/methyl_augmented_filter_phase2/cycle3/ablation/figures/cycle3_step1_5_min_ablation.png`
- Evidence ledger: `cycle_id: 20260520_cycle3_step1_5_ism_ablation_vestigial` (line 48)
- Gate config: `InterSubMod/research/methyl_augmented_filter_phase2/cycle3/cycle3_caller_f1_gate.json` (features_production 4-feature)

### Pre-reg prediction 命中

plan v2.1 事前預測 55% A1 PASS + M1a FAIL → 實際 A1 MARGINAL (between PASS-MARGINAL) + M1a FAIL → 命中方向，避免 confirmation bias 解讀。

---

## 2026-05-19 — Phase 2 Cycle 2 結束 (⭐3 保持 + caveat) + Cycle 3 啟動

**Session anchor**：2026-05-18 cycle 1 ⭐3 strong → 2026-05-19 cycle 2 cross-sample DIRECTION_NEGATIVE + cross-binary PASS → user 5/19 三選決議

### Cycle 2 結論

| Hypothesis | Verdict | 細節 |
|---|---|---|
| **H_C1_5** cross-sample n=5 Wilcoxon | **DIRECTION_NEGATIVE** | Transfer 1+/4- p=0.1875 ΔRecall p=0.0625；re-fit 3+/0-/2≈0 MIXED；僅 HCC1937 (F1=0.37) re-fit +0.00761 |
| **H_C1_6** V3F/V5/V6 HCC1395 cross-binary | **PASS** | max_var transfer 0.00073 / re-fit 0.00055；V6 re-fit bit-exact reproduce cycle 1 drift 0 |

**Mechanism**: caller-F1-headroom-bounded — 3/4 新樣本 caller F1 > 0.83 + FP density < 4% 不留 filtering room；HCC1954 transfer 災難 -0.377 = caller_af coef overfit HCC1395 AF 非 methylation 失敗。

### User 5/19 三選決議

| 決策 | 選項 |
|---|---|
| Cycle 1 tier | **保持 ⭐3 + caveat** — internal valid + cross-sample bounded |
| HTML PI preview | **等 cycle 3 結果後再做** — 避免在不完整故事上 PI 決策 |
| Cycle 3 方向 | **啟動 Caller-F1-headroom-aware redesign** — IF caller F1<0.80 + FP density>0.10 apply per-sample re-fit ELSE skip |

### 新增四軌狀態（patch v2）

| 軌道 | 5/18 狀態 | 5/19 狀態 |
|------|----------|----------|
| methyl_filter_phase2 | cycle 1 ⭐3 strong | **cycle 2 ⭐3 + caveat (cross-sample bounded)，cycle 3 caller-F1-headroom-aware 啟動** |
| thread_d_paper | 5/22 Tier 2 啟動 | unchanged |
| selfphasing_v6_production | 5/19-22 4-day | **5/19 整合 H_C1_6 PASS evidence (V6 zero F1 regression)** |
| phase_block_3d | 5/23 init-research | unchanged |

### Cycle 3 Plan 概要

- **H_C3_1**: Caller-F1-headroom rule (caller F1 < 0.80) qualifying subset (HCC1395 + HCC1937) re-fit mean ΔF1 ≥ +0.01
- **H_C3_2**: 擴 low-F1 panel n≥4 Wilcoxon p<0.05
- **H_C3_3**: High-F1 (>0.83) 樣本 skip filter caller F1 preserved drift < 0.001

### Cycle 2 artifacts

- Coordinator synthesis: `InterSubMod/research/methyl_augmented_filter_phase2/cycle2/cycle2_findings.md`
- H_C1_5 detail: `InterSubMod/research/methyl_augmented_filter_phase2/cycle2/cycle2_step_b3_b4_findings.md`
- H_C1_6 detail: `InterSubMod/research/methyl_augmented_filter_phase2/cycle2/cycle2_step_c1_h_c1_6_sanity.md`
- Evidence ledger entry: `cycle_id: 20260519_phase2_cycle2_cross_sample_negative`

### 🔍 5/19 (晚) Multi-Agent Audit 揭露 + Cycle 3 啟動延後

**4-agent Explore audit（平行 × 4）校對 cycle 2 claim 與原始資料**：

- ✅ **數字 fidelity 100% verified** — line 23-26 所有數字（transfer/refit/sanity max_var）在原始報告精確找到
- ⚠️ **Mechanism 量化證據不足** — HCC1954 caller_af overfit + caller-headroom-bounded 定性 OK 但**缺視覺化**（AF 分佈 / coef inspection / scatter plot）；paper draft 前需補
- 🔴 **Plan v2 R-MENTAL-DRIFT 0% compliance** — 48hr cooling-off + NEGATIVE postmortem 完全缺；同日 cycle 2 NEGATIVE → cycle 3 啟動違反 plan v2 自定紀律
- 🔴 **cycle_id narrative bias** — `_negative` 預定 verdict（append-only ledger 精神違反）
- 🔴 **Cycle 3 H_C3_2 low-F1 panel 樣本名單不存在** — 僅 HCC1395+HCC1937，缺 2-3 個；無 BAM 可用驗證
- 🔴 **Provenance 缺漏** — binary_commit / dataset_id / pre-reg link 全缺

### Cycle 3 啟動延後（plan v2 R-MENTAL-DRIFT 紀律執行）

| 動作 | 狀態 |
|------|------|
| 48hr cooling-off period | ✅ 2026-05-19 → 2026-05-21 |
| NEGATIVE postmortem .md | ✅ `InterSubMod/research/methyl_augmented_filter_phase2/cycle2/cycle2_NEGATIVE_postmortem.md` |
| evidence_ledger postmortem entry | ✅ cycle_id `20260519_phase2_cycle2_negative_postmortem` |
| Cycle 3 H_C3_2 low-F1 panel 樣本名單 | ⏳ 待 explore archive 或 user 提供 |
| Cycle 3 啟動時間 | 🟡 **2026-05-21 後**（cooling-off 完成 + panel 名單就緒 + 4 must-fix 解除） |

### Cycle 3 啟動前 4 個 Must-Fix

1. 🔴 H_C3_2 low-F1 panel 具體 4 樣本名單（含 BAM 可用性驗證）
2. 🟠 cycle_id 重命名規範：未來 cycle 用 `_validation` 不 `_negative`（verdict 在 entry 不在 name）
3. 🟠 ledger entry 補 binary_commit hash / dataset_id / pre-reg link
4. 🟡 H_C3_1 target 計算澄清：cycle 2 re-fit mean (+0.015) vs cycle 3 重新 re-fit

**Plan v2 R-MENTAL-DRIFT 紀律檢驗**：本次為 plan v2 commitment 的真實測試 — cycle 2 NEGATIVE 後若立即啟動 cycle 3 屬於 confirmation loop。經 4-agent audit 揭露 + 用戶 explicit 暫停決策後**正式延後執行**，紀律保持 commitment 而非表演。

---

## 2026-05-18 (晚晚) — Plan v2 patch（cycle 1 strong 後四軌定案）

**Session anchor**：7 輪 Socratic 燒烤 (5/15+5/18) 累積 27 燒烤點 → plan v2 寫入 `~/.claude/plans/tender-pondering-blossom.md` → cycle 1 strong 結果 patch

### v2 核心變更（vs v1）

- 雙軌 → **四軌平行**（thread_d + V6 + phase_block_3d + methyl_filter_phase2 cycle 2）
- 5 目標**階段化 active**（A→G4+G5, B→G2, C→G1+G3）取代 v1「降權」二元語言
- 三階段 paper（A framework W3-W8 / B clone preprint M2-M5 / C ASM+two-hit M6-M12）
- B/C 從 sub-analysis 抽取（不新建 project）
- methyl_filter_phase2 升 4th 軌（不 merge phase_block_3d，因 cycle 1 ΔF1=+0.02236 strong）

### 四軌定案

| 軌道 | 狀態 |
|------|------|
| thread_d_paper | 5/15 init, Tier 2 啟動 5/22+ |
| selfphasing_v6_production | 5/15 init, 4-day workflow 5/19-22 |
| phase_block_3d | 5/18-19 inject 3 H, 5/23 init-research |
| methyl_filter_phase2 | 5/15 init, **5/18 晚 cycle 1 ⭐3 strong**, cycle 2 cross-sample 待 V6 4 樣本 ISM |

### Phase A/B/C timeline（四軌延長為 6-8 週）

| Phase | 時程 | Deliverable | G 對齊 |
|-------|------|-------------|--------|
| A | W3-W8 (5/19-6/30, 6-8 wk) | framework draft + V6 + phase_block_3d + methyl cycle 2 + G1/G2/G3 skeleton | G4 + G5 |
| B | M2-M5 (7/1-9/30) | framework 投稿 + clone preprint（thread_d §4 抽取） | + G2 |
| C | M6-M12 (10/1+) | ASM short + two-hit notes（phase_block_3d X1/X2 抽取） | + G1 + G3 |

### 3 Hypothesis 並託 phase_block_3d（待 5/18-19 inject）

- **H013 (X1)**: Phase block 邊界內 CN-stratified methylation 高一致 → TP enrichment (cis-effect)
- **H014 (X2)**: Phase block 跨越大 CN segment 時 methylation 不連續 → FP marker (boundary artifact)
- **H015 (X3)**: Joint multi-feature zone score: NG × CN tier × methylation rank-2 — **methylation 為 weak axis (rank-5 per methyl cycle 1)，X3 zone-stratified 探索 zone 內 methyl 是否強化**

### 啟動時序

- **5/18-19**: inject H013/H014/H015 進 hypothesis_queue.json（純 metadata, 5 min）
- **5/19-22**: V6 production 4-day workflow 獨享 focus
- **5/22 後**: thread_d_paper Tier 2 + phase_block_3d init-research scaffolding + methyl_filter_phase2 cycle 2

### 5 Risk Items

- R-CONFIRMATION-LOOP: cycle 1 結果決定（5/18 晚 strong → (a) 路徑通過）
- R-A-SCOPE: A 階段 G1/G2/G3 skeleton only, no deep dive
- R-PAPER-THIN: B/C 接受 thin（preprint 不是 main journal）
- R-OPS-CREEP: ops <2 hr/week
- R-MENTAL-DRIFT: 下次 reactivate 前 48 hr cooling-off + NEGATIVE postmortem

### Mental Model Shift 日誌

| Date | Shift | Trigger |
|------|-------|---------|
| 5/15 | G3/G5 降權（plan v1）| 4 輪 Socratic 燒烤 |
| 5/18 上午 | G3/G5 reactivate（plan v2 階段化 active）| methyl filter v1.0 ΔF1=+0.00242 marginal |
| **5/18 晚** | **methyl_filter_phase2 升 4th 軌** | Phase 2 Cycle 1 ΔF1=+0.02236 strong (9.24× v1.0, 2.24× Cohen ribbon) |

---

## 2026-05-18 (晚) — Phase 2 Cycle 1 Global FP Filter 完成 (⭐3 strong, ΔF1 +0.02236 = 9.24× v1.0)

**Session anchor**: v1.0 cycle ⭐3 marginal → pivot global FP exploration → ⭐3 strong (HCC1395-internal validated)

### 完成項目（Phase 2 Cycle 1, ~30 min wall clock, Track A only）

- ✅ **Step 0 Global FP audit** (Agent A1, ~0.7 min): D1 94.22% top 10 cells / D2 ΔF1 +0.02637 Strategy B / D3 heterogeneous FAIL / D4 high-AF FAIL → **Path B selected** (pure global LR)
- ✅ **Step 1 Filter design** (Agent A2): VIF=217 → drop NumReads_master + L2 C=1.0；NaN MNAR confirmed → impute correct；final filter 10 features
- ✅ **Step 2 HCC1395 ΔF1 verdict**: **ΔF1 +0.02236** (9.24× v1.0)，τ*=0.39，Step 5c lost TP 81% rescued
- 🔵 **Track B (cross-sample) DEFERRED to cycle 2**: 4 樣本 V3F/V5 BAM 物理不存在
- 主報告: [InterSubMod/docs/experiments/in_progress/2026/05/20260518_Phase2_Cycle1_Global_FP_Filter_01.md](experiments/in_progress/2026/05/20260518_Phase2_Cycle1_Global_FP_Filter_01.md)
- 研究目錄: `research/methyl_augmented_filter_phase2/cycle1/`

### 關鍵發現

**🎯 ΔF1 +0.02236 從 v1.0 marginal 升至 substantial** (2.24× Cohen 小 effect ribbon)
- caller_af (+3.44) > LOH_inner (+1.46) > Cov (+1.27) > NG (+1.07) > **HPFineF (+0.75)** + 4 個更小 methyl coef
- TP loss 只 1.56% / FP removed 70.20% → Precision 0.9541 / Recall 0.6030

**🔍 Step 5c lost TP 17/21 (81%) rescued by cycle 1 global LR** vs v1.0 cell-gated 0% — 證明 v1.0 over-restrictive

**⚠️ Methylation 是 5th-rank covariate** — filter 主導為非-methyl 軸；reframe "multi-axis filter incl. methylation"

**Multi-seed std = 5e-5** (20× below threshold) — high intra-sample stability

**⚠️ R-Step0 5 caveats**: 4 RESOLVED (NaN MNAR / VIF collinearity / lost TP overlap / methyl marginal) + 1 OPEN (HCC1395-only)

### 評估與後續

- **Tier ⭐3 strong** (HCC1395-internal validated): 3/4 H PASS + multi-seed stable + lost TP 81% rescue
- ⭐3 而非 ⭐4 因 cross-sample DEFERRED (V3F/V5 4 樣本 BAM 不存在)
- **⭐4 升級必要條件 (cycle 2)**: V6 4 樣本 ISM rerun (~3.2 hr, BAM 已存在 v6_5sample_extension/) + Wilcoxon n=5 + HCC1395 phaseC V3F/V5/V6 三向 cycle 1 filter apply
- **Cycle 2 可與 V6 production 4-day workflow 共用 ISM data** (Day 2 6-sample marker coverage)

---

## 2026-05-18 — Methylation-Augmented FP Filter Pilot 完成 (⭐3 PARTIAL POSITIVE marginal) + Phase 2 project init

**Session anchor**: 跨越 characterization → filter F1 evaluation；plan v1.0 cycle 完成 + phase 2 project scaffolded

### 完成項目（v1.0 cycle, ~80 min Step -1 ISM rerun + ~1 hr 後分析 multi-agent）

- ✅ Step -1: phaseC ISM 12 runs 重跑 with significance (移除 `--no-distance-matrix` flag) — 80 min, 12/12 全成功
- ✅ Step 0: augmented master TSV (35,332 × 202 cols, 13 methylation features × V3F+V5+V6 × off/on)
- ✅ Step 1+2: 138 augmented LR + LRT + 12 FP-rich τ sweep — H1 16/30 cells q<0.05 / H5 V5≈V6 Δβ=1.87e-5
- ✅ Step 3: ΔF1 vs **ClairS-TO ssrs caller F1=0.7166** (per 09_V6_caller_F1_verification.md, tumor-only mode, ClairS-TO v0.3.0 ssrs model + LongPhase-TO + ISM pipeline; ⚠️ 早期 reports 誤標 "paired-pileup" 已 5/20 framing patch — canonical alias = `clairs_to_ssrs`, see `InterSubMod/docs/data_specs/20260411_工作區命名與目錄結構_01.md` §1.5) — **+0.00242 @ τ*=0.52 marginal** (POSITIVE-but-marginal)
- ✅ Step 4: 13 mechanism candidates + 14 PubMed refs (H4 POSITIVE relaxed gate)
- ✅ Step 5c TP rescue **NEGATIVE** — 95.2% lost TP 是 low-AF subclone，methylation 訊號實為 caller_af proxy
- ✅ Step 5d robustness GREEN with caveats — ΔF1 std 2e-5 stable, 4 unique LRT cells, NaN borderline
- 主報告: `InterSubMod/docs/experiments/in_progress/2026/05/20260518_V6_Methyl_Filter_Pilot_01.md`
- 研究目錄: `research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/`
- **Phase 2 project init**: `research/methyl_augmented_filter_phase2/` (4 H, 3-5 cycles 預期)

### 關鍵發現

**⭐3 Verdict (ALL 5 H POSITIVE but ΔF1 marginal)**:
- H1 LRT q<0.05 在 16/30 cells (top p=1.8e-58)
- H2 ΔF1 +0.00242 < +0.005 Cohen ribbon "marginal"
- H3 FP_removal 98.3% > TP_loss 35.0%
- H4 13 mechanism × 14 PubMed refs (cis-mQTL/cancer ASM/allele-imbalance/repeat/replication timing)
- H5 V5/V6 LR β median diff = 1.87e-5 (V6 重用 V5 phased VCF)

**Prior art (TumorLens/ROCIT/SGZ/Wakhan/SAVANA) 全無同口徑** — read-level multi-axis methylation-augmented filter 為新貢獻

**TP rescue NEGATIVE 機制解釋**: lost TP 與 removed FP 在 caller_af 空間嚴重重疊；methylation 沒有獨立 axis 訊號 → rescue rule 必然 reimport FP

**Plan v1.0 → v2.0 pivot**:
- v1.0 cycle 完成 ⭐3 marginal
- v2.0 plan (2026-05-18 批准): Phase 2 Cycle 1 pivot 為 global FP exploration + heterogeneous threshold filter (Track A) → cross-sample 4 樣本 ISM rerun (Track B) — 目標 ΔF1 ≥ +0.01

### 評估與後續

- **Tier ⭐3 PARTIAL POSITIVE (marginal)**: 5 H POSITIVE 但 effect size < +0.005 → 需 cross-sample 才能升 ⭐4
- **與 V6 production 4-day workflow 整合**: phase 2 cycle 1 Track B (4 樣本 V3F+V5 ISM ~1 hr) 與 V6 production COLO829 V6 ISM (~2 hr Day 1) 可平行；7-sample marker coverage data 共用
- **Pre-registration 5 H 已寫入 Phase 2 manifest** (research/methyl_augmented_filter_phase2/manifest.yaml)

---

## 2026-05-18 — Agent Harness Audit P0-P4 完整收尾 (11 commits)

**Session anchor**: 2026-05-18 single-day sprint 完成 InterSubMod agent harness 7-phase audit (P1-P7) + 5-tier fix (P0-P4) 共 29 fix items × 11 commits × ~10 hr。

### 完成項目（依 commit chain）

| Commit | 範圍 |
|--------|------|
| `ee648fb` | G1+G4 hooks（SessionStart + skill_change_audit）+ verification_guide |
| `61055b8` | 7-phase audit deliverables (P1-P7 HTML/JSON × 5) |
| `d5db8dc` | P0 critical (5 YAML invalid + 16 silent failure + query template) |
| `a7c1495` | P1 (13 D2 + 5 D3 cross-ref + 5W2H + SMART + D7-1) |
| `1a9379e` | P2 (TL;DR / SWOT / allow_audit / cache_telemetry / ledger v2 / manual) |
| `d64c0e7` | P3 (worktree / Eisenhower / subagent_logger / rules paths / compact_test + H4 no-op) |
| `43d8b28` | P4 Industry Deep Audit (Anthropic + OpenAI + Walking Labs × architect/researcher cross-validate) |
| `3ff0980` | P4 Top 3 (Evaluator agent + Evidence gate + Watch dashboard) |
| `93659a8` | P4 剩 7 (E5 5→3 入口 / E7 claim hook / E9 recall / E10 injection / cold-start / cleanup / spec) |
| `6998470` | wrap-up (gitignore + log baseline) |
| `56e00c7` | 00_COMPLETION_REPORT.md |

### Verification 結果

- **Cache hit rate 96.8%**（業界 Anthropic claim 90% 超越）
- **30 hooks 跨 6 events**（silent failure 0/30）
- **42/42 skills YAML valid**
- **Cold-start test Q1-Q5 全 ✅**（Walking Labs L03）
- 4 audit HTML reports（P1 / P2 / P3 / P4 industry / P4 final）
- 11 new hook scripts + 3 new templates + 1 evaluator agent + 1 SKILL spec

### 完整 deliverables 索引

`InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/00_COMPLETION_REPORT.md` — 282 行完整 §1-§10 結構化文件

### 業界對齊

| 維度 | InterSubMod 達標 |
|------|---------------|
| Anthropic 3-agent harness | ✅ evaluator agent 已建 |
| Walking Labs 12 lectures | ✅ 12/12 對齊 (L03 cold-start / L06 init / L11 observability / L12 clean state) |
| OpenAI 3 pillars | ✅ Context / Constraint / GC (recall_logger + skill_audit) |
| cwc-long-running-agents | ✅ evaluator + verify_gate + evidence_tracker |
| Prompt Caching 90% claim | ✅ 96.8% 超越 |

### 下一階段

V6 production 4-day workflow (W3 deadline 5/22) — `InterSubMod/research/selfphasing_v6_production/4day_compressed_workflow.md`

---

## 2026-05-17 — Tier 1-4 序列化執行 + T1.1/T1.3 完成

**Session anchor**：4 輪 Socratic 燒烤對話收斂 9 條決策 + plan `~/.claude/plans/tender-pondering-blossom.md`

### Plan tender-pondering-blossom Tier 1-4 進度追蹤

```
Tier 1 (W3 2026-05-15~22) 必須前置
  T1.1 Thread D 主軸正名         ✅ DONE 2026-05-16
       → 加 banner + §2.5 paradigm reframe + 主標題改名
       → InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_*.md (338→381 行)
       → 新名稱「TP-enriched phasing signatures (LOH × cross_het)」
  T1.2 V6 production tag finalize 🔴 Hard Gate (待執行)
       → 詳見 InterSubMod/research/selfphasing_v6_production/00_PLAN.md §Tier 1.2 workflow
  T1.3 init-research scaffolding ✅ DONE 2026-05-16
       → InterSubMod/research/thread_d_paper/         (165 行 00_PLAN + 108 行 manifest)
       → InterSubMod/research/selfphasing_v6_production/ (154 行 00_PLAN + 122 行 manifest)

Tier 2 (W3-W4) 證據強化
  T2.1 Z-AUTO KDE 跨 4 樣本擴展     ⏳ ⭐3 → ⭐4 升級必要條件
  T2.2 HCC1395 primary discovery 章節骨架  ⏳ Strategy A §3
  T2.3 6-sample replication cohort 章節骨架 ⏳ Strategy A §4（含 HCC1954/HCC1937）

Tier 3 (W4-W6) Paper draft + 工程平行
  T3.1 Paper full outline (abstract + 6 章 + 6 主圖)  ⏳ Tool paper (Bioinformatics / NAR GB)
  T3.2 GitHub repo 整理 + public 化      ⏳ Reproducible release
  T3.3 Docker image build + tutorial      ⏳ 1-hour install + run
  T3.4 Benchmark suite 公開化             ⏳ HG002 subset or HCC1395 SEQC2 公開部分
  T3.5 Discussion 章節                    ⏳ 63% gap + cancer-only + Normal BAM future

Tier 4 (W6+) reactive 擴展
  T4.1 Phase 2A Normal BAM cross-sample   ⏳ G4 characterization (45% → 70%)
  T4.2 GC/mappability/repeat 新軸 pilot   ⏳ reactive (if reviewer 質疑 framework gap)
  T4.3 PI Report 4-29 errata + V6 sign-off ⏳ T1.2 完成後一併打包（written email）
  T4.4 HG002 non-cancer pilot             ⏳ reactive (if reviewer 強質疑 cancer-only)
```

### T1.2 V6 Production Tag — Hard Gate 5-day workflow

詳見 `InterSubMod/research/selfphasing_v6_production/00_PLAN.md` §"Tier 1.2 V6 Production Tag Workflow"：

| Day | Action | Gate level |
|-----|--------|------------|
| 1-2 | COLO829 V6 ISM 補完（Archive TO rerun + KDE-corrected） | 🟢 normal |
| 3   | 7-sample marker coverage + caller F1 比較 (V3F vs V5 vs V6 vs SEQC2) | 🟢 normal |
| 4   | Binary commit hash 寫 `manifest.yaml` | 🟡 review |
| 4   | `git tag v6-prod-{YYYYMMDD}` | 🔴 **Hard Gate**（不可逆） |
| 5   | PI errata 5 條 + V6 sign-off written email draft | 🟡 review |
| 5   | User review email → send | 🔴 **Hard Gate**（送出後不可逆） |

**T1.2 Gate 通過後解鎖**：
- thread_d_paper Tier 2 Archive TO 7-sample rerun (V6 binary, ~10 hr parallel)
- T4.3 PI errata package 同期完成

### 已收斂的 9 條決策（plan §"4 輪燒烤確認的 9 條決策"）

主軸與目標：(1) 主軸雙軌平行；(2) G1/G2/G4 → characterization-only / G3 暫緩 / G5 降權
論文：(3) Tool paper（Bioinformatics / NAR GB）；(4) 核心 contribution = read-level framework；(5) 主軸正名「TP-enriched phasing signatures (LOH × cross_het)」；(6) Strategy A 樣本階層 HCC1395 primary + 6 replication；(7) 接受 37% framework + 63% gap discussion；(8) Cancer-only 接受 limitation
釋出：(9) 完全 reproducible（binary + Docker + benchmark suite + GitHub）

### 與下方 2026-05-13/2026-05-15 區塊的關係

- **2026-05-13**：原 3 週 V6 序列化估計（已 deprecated；保留 historical reference）
- **2026-05-15**：multi-agent fan-out paradigm reframe（evidence base，本區塊 T1.1 + T1.3 的執行依據）
- **2026-05-17（本區塊）**：plan tender-pondering-blossom 細化 + 序列化執行追蹤（live progress）

---

## 2026-05-15 — V3F/V5/V6 ISM 三向 × LOH × HP × CN characterization 完成（⭐3 PARTIAL POSITIVE）

**Session anchor**：`"驗證ISM"` (multi-agent fan-out A/B/C/D/E + Coordinator)

### 完成項目（HCC1395 pilot + 4 樣本擴展，~3.5 hr 全執行）

- ✅ phaseC 12 個 ISM run 整合 → step1_master_three_way.tsv (35,332 rows × 64 cols)
- ✅ 3 軸 50-cell grid + power gate (46% powered) + LR + 7 道 confound guard
- ✅ 4 FP zone deep dive (Z-OCH / Z-CHR8 / Z-GL / Z-AUTO)
- ✅ 4 樣本擴展 V6 ISM (H1437/H2009/HCC1954/HCC1937; COLO829 deferred)
- ✅ Prior art (TumorLens/ROCIT/SGZ/Wakhan/SAVANA) 全無同口徑等效物
- 主報告：[InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md](experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md)
- 研究目錄：`research/v6_bam_tpfp_hp_loh_cn/` (00_PLAN + 01_data + 02_methodology + 02_prior_art + step1-4)

### 關鍵發現

**🤯 Paradigm reframe**: 2/3 預設「FP-rich」zone 實際 **TP-pure signatures**（不是 FP markers）：
- Z-OCH (Outer cross_het): FP rate 0.017 << global 0.137 (Fisher p=3.8e-62 for **TP**-enrichment) — cross_het = somatic-evidence marker
- Z-GL (Inner gain+LOH): FP rate 0.003 (0.022× global) — gain on somatic hap = somatic signature

**✅ H4 POSITIVE — chr8 hotspot CN+AF 主導**:
- LR deviance: caller_af 0.393 > **CN 0.211** > HP 0.063 > LOH 0.038
- (LOH+CN) − HP = +0.186 (3.7× threshold 0.05)
- chr8 FP enrichment 2.31× (highest of 23 chr)，捕獲 20% 全 FP，但 sample-specific

**🔍 V5 over-promote 直接證據**: Inner LOH NG=2 region V5=8,136 (+60% over V3F 5,064)，TP rate 沒升；V6 修補回 V3F 水準 (5,353)。V5/V3F top cell ratio 達 5.95× (Inner|cross_het_inv|cov_normal)，**集中 cross_het bucket**（same_HP* 全 ~1.0×）→ Layer 1.5 機制只在 somatic-fallback heterozygous reads 作用

**⚠️ Framework coverage gap**:
- Z-CHR8 + Z-AUTO 共捕獲 ~37% 全 FP
- **剩 ~63% FP 不被此 framework 解釋** → 需新軸（mappability/repeat/GC/SV）

**Cross-sample (n=5)**: 唯一 signature candidate `Outer|other|cov_high_gain` 5/5 同方向 (Wilcoxon p=0.0625, Δ=+0.0069)；受 caller 飽和 (≥0.998) 限制。HCC1937 BRCA1 outlier chr15/17/14 driver；chr17 FP 不專屬（HCC1395 共有）

### 評估與後續

- **Tier ⭐3 PARTIAL POSITIVE**：H4 POSITIVE + paradigm reframe + 1 cross-sample candidate ✅；但 H7 confound guard 通過<5 + Z-CHR8 sample-specific + 63% FP unexplained → 未到 ⭐4
- **升 ⭐4 需做**：(a) Z-AUTO KDE 跨 4 樣本各自做（驗證機制是否 recur）、(b) 加新軸測 63% unexplained FP
- **V6 production tag finalize**（CURRENT_FOCUS W3 phase 1 gate）：本 cycle 確認 V6 marker coverage +9% over V3F + caller F1 不變 → 可作 V5→V6 production 升級依據

---

## 2026-05-13 — 主軸切換：V6 Production → Thread D Paper（4-6 週序列化雙軌）

**Session anchor**：本次大盤點對話 `"未來目標方向的釐清與探索"`

### 決策摘要

- **主軸切換**：Self-Phasing V6 production 化（Track B 先）→ Thread D Paper 主軸（Track A 後）
- **G5 處理**：保留但降權；F1 路線實質關閉，等 V6 + Phase 2 rerun 後再評估一次
- **PI 報告 4-29 errata**：延後到 V6 production 完成後（W3 末）一併打包；目前不主動推進
- **論文期刊**：暫不鎖定，先寫 framework paper 骨架
- **Scaffolding**：建立 `research/thread_d_paper/` + `research/selfphasing_v6_production/` 兩個 init-research 專案目錄
- **Archive TO 6 樣本 rerun binary**：等 V6 production tag 完成後一次到位（避免重跑）

### 4-6 週序列化雙軌時程

```
Phase 1 (W1-3): Track B — Self-Phasing V6 production 化
  W1: F-paired-D3 V5 Layer 1.5 ISM 影響量化 (chr19+全基因組)
  W2: V6 7-sample expansion (1-2 day workload)
  W3: V6 production tag finalize + binary commit hash 寫 manifest
      Gate: V6 marker coverage > V3F/V5 且 0 critical regression
      Deliverable: PI 報告 errata package 一併出（含 4-29 5 條 errata + V6 sign-off）

Phase 2 (W3-5): Track A — Thread D evidence 升 grade A
  W3-4: Archive TO 6 樣本 ISM rerun (V6 binary, ~10 hr parallel)
  W4:   Wilcoxon n=7 重算（同方向預期 p<0.0156）
  W4-5: HPFineNGroups marker × 兩 flag 跨樣本驗證
        HCC1954 outlier deep panel + Normal BAM pilot (R17)

Phase 3 (W5-6+): Thread D framework paper draft
  W5-6: paper framework 骨架（暫不鎖期刊）
        - 主軸：read-level epigenetic characterization
        - G1 ASM 全域章節 + G2 cross-sample 4-group subclone
        - HCC1954 專章（outlier + Thread B 撤回脈絡）
```

### 五大目標成功標準切換（2026-05-13 正式生效）

| 目標 | 舊標準 | 新標準 |
|------|--------|--------|
| G1 per-CpG ASM | F1 提升 + 可解釋性 | **僅可解釋性**（characterization-only） |
| G2 clone 結構 | F1 提升 + 可解釋性 | **僅可解釋性**（4-group subclone characterization） |
| G3 二次打擊順序 | 推論成功 | **依賴 G1 突破，暫緩** |
| G4 TO normal 補強 | F1 提升 + 可解釋性 | **僅可解釋性**（R1-Global 已 NEGATIVE for F1） |
| G5 evidence panel F1 | F1 提升 | **保留 ceiling +0.0112；不投入新資源** |

### 已關閉死路（不再投入）

- 純 methylation TP/FP filter (CL-008 ⭐5)
- LOH/CNV zone-aware filter (CL-024 ⭐5)
- FN methylation rescue (O9 ⭐5)
- S3 cross-sample whitelist (Thread B 撤回 2026-04-26)
- TO single-sample post-hoc germline FP filter
- G3 事件順序（暫緩等 G1 突破）

### Immediate next actions（執行序）

1. ✅ 2026-05-13：本區塊已寫入 CURRENT_FOCUS.md
2. ⏳ 待執行：`init-research` 建立 `research/thread_d_paper/` + `research/selfphasing_v6_production/` scaffolding
3. ⏳ 待執行：W1 F-paired-D3 V5 Layer 1.5 ISM 影響量化（chr19+全基因組）

### 為何序列化而非完全平行

Archive TO 6 樣本 rerun 直接吃 V6 binary 可避免重跑；V6 production 化包含 ISM 影響量化（F-paired-D3）本身就需要先量化 V5→V6 改變才能 sign-off — 兩者實質有 binary 依賴關係。序列化讓 W3 後 Phase 2 一次到位，整體節省 ~10 hr × 1 次重跑成本。

---

## 2026-05-11 — html-preview skill shipped (Phase 2)

✅ `InterSubMod/.claude/skills/html-preview/` shipped. End-to-end demo at
`InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01/index.html`
(topic-folder mode, 281-line source).

Spec: `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md` (D1-D20)
SOP: `InterSubMod/docs/references/manual/20260511_HTML_MD_PPTX輸出格式SOP_01.md`
Phase 2 plan: `InterSubMod/docs/plans/2026/05/20260511_Phase2_html_preview_implementation_01.md`

Ready to plan Phase 3 (Tier A 6 skill 接入 html-preview as companion).

> **主軸（2026-04-26 切換）**：**Thread D LOH-constrained phasing signatures**（TO 層論文主軸）。主軸報告：[InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md](reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md)。Thread B（LOH × AF × CN 跨樣本 whitelist filter）已於 2026-04-26 正式撤回 filter 用途宣稱，HCC1395 single-sample case study 與 per-sample characterization 保留。撤回宣告：[InterSubMod/docs/reports/validated/2026/04/20260426_Thread_B_Retraction_Whitelist_Cross_Sample_01.md](reports/validated/2026/04/20260426_Thread_B_Retraction_Whitelist_Cross_Sample_01.md)。

## 1. 目前狀態

1. docs 重整已完成核心落地：
   - 命名：`YYYYMMDD_主題_流水號.md`
   - 報告分層：`reports/validated|finalized`
   - 實驗分層：`experiments/in_progress|validated|finalized`
   - 資訊分層：Active / Recent / Archive（詳見 AGENTS.md）
2. `output/` 入口已固定為軟連結：
   - `output -> /big7_disk/liaoyoyo2001/big7_disk_output`
3. Knowledge MCP 已接入：
   - `.mcp.json` 指向 `/big8_disk/liaoyoyo2001/knowledge/scripts/mcp/knowledge_server.py`

## 2. AI Agent 主要入口

1. 啟動壓縮上下文：`docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md`
2. docs 導航：`docs/README.md`
3. 研究歷史索引：`docs/experiments/INDEX.md`
4. 研究全景：`docs/reports/research_landscape/00_INDEX.md`
5. Agent 手冊：`docs/references/manual/20260301_AI_Agent_快速操作手冊_01.md`
6. 健康檢查：`scripts/analysis/check_ai_agent_readiness.sh`
7. 文件規範：`docs/standards/README.md`

### Agent 上下文控制面分工（2026-04-27）

| 入口 | 角色 | 不應承擔 |
|------|------|----------|
| `AGENTS.md` | repo 內硬規則、資料/輸出位置、Knowledge Base 查閱義務 | 每週研究結論細節 |
| `.claude/CLAUDE.md` | Claude Code 執行策略、hooks、確認矩陣、壓縮保留規則 | 正式研究結論來源 |
| `docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md` | 啟動壓縮上下文、重要數據、任務順序、待決策矩陣 | 取代 validated report 或 Knowledge Base |
| `docs/CURRENT_FOCUS.md` | live 主軸、阻塞、入口索引 | 長期歷史總覽 |
| `research/autoresearch/research_direction.md` | 候選 queue 與 research-loop 定向 | 自動執行觸發 |

## 3. 當前進行中

### 2026-05-10 Self-Phasing 整合主軸（5 個 commit 鏈完成）

**主軸 commit 鏈（5/8 → 5/10）**：
- `951e7c9` 5/8 Self-Phasing 完整觀察整合報告（10 段 + 5 figures + 1202 行）
- `f17754f` 5/9 PI 報告 4-29 errata companion + 原報告 banner（4 條 errata）
- `6ed8a0d` 5/9 paired audit Step A+C — paired 沒 priority bug
- `766ec5f` 5/9 paired audit Step D — V5 Layer 1.5 設計缺陷揭露
- `df5137e` 5/10 整合 5/9 paired 發現至 5/8 主報告 §8.6
- `2553e96` + `71d21bd` 5/10 補強 E5 errata + F6+F7 figures

**主結論**：
- self-phasing 機制因果確立（17.3:1 + 34,855 read-level victims + 100% V3F/V5 修正）
- V5 整體可作 production tag baseline
- **5/9 新發現 V5 Layer 1.5 在 germline-absent 區域繼承 priority bug 偏移（4.19:1 偏 HP1，與 baseline 完全相同）；V3F 標 hp=33 反而更穩健** — 設計缺陷待 F-paired-D3 ISM 影響量化
- PI 報告 4-29 5 處 errata 已 patch（companion + banner）

**待 follow-up（5/9 新加 4 條）**：
- F-paired-D1：全基因組 germline-absent 擴展（chr19 → 全 chr，~150K events）
- F-paired-D2：phase block 內 axis-aligned 分析
- F-paired-D3：V5 Layer 1.5 改回 V3F 的 ISM 影響量化
- F-paired-D4：E5 PI errata 補強 ✅ DONE 5/10 (commit 2553e96)

主入口：[InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md](reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md)

### 2026-04-23 週報後 P0/P1 行動（之前主軸，部分仍 active）

| 優先 | 行動 | 預期產出 | 估時 | 依據 |
|:---:|------|---------|:---:|------|
| **P0** | **Paired normal-pilot obs18 對照**：HCC1395 paired mode 做 NG=2 same-hap vs cross-het 拆分，驗證 H-D4 gap 是否在 germline-排除條件下消失 | `step7_hcc1395_normal_paired_obs18.md` + Wilcoxon signed-rank 初步統計 | 1-2 天 | 本週 Wave 1 review 指出 Thread C 僅驗證 H-D1；H-D4 gap-disappearance test 待 paired 對照 |
| **P0** | **Archive TO 6 樣本重跑 ISM**（含新 KDE + LOH.bed + germline-hp-only default=off）| 擴充 `master.tsv.gz`；跨樣本 S1-S7 scheme 驗證；HPFineN marker master × 兩 flag 對照 | ~10 hr parallel | 本週 Archive 文件 `20260422_Archive_TO_Rerun_ISM_Requirement_01.md` |
| **P1** | **HCC1954 outlier 專項分析**：Outer cross-het TP 0.08 根因（Potential_LOH 可靠性 / AF / CovM / IGV） | HCC1954 focused report | 1 天 | obs18 6 樣本中唯一 outlier |
| **P1** | **Formal stats on NG=2 gap**：Wilcoxon signed-rank on 6/6 samples + bootstrap CI | Thread D 證據卡 tier 升級 | 0.5 天 | Wave 1 指出目前僅有 descriptive stats |
| **P2** | ~~**S4 內部二級判別 pilot**：HCC1395 TO S4 subset (n=30,432, 8,780 FP) 加入 ReadParser 特徵跑 logistic regression~~ `[RETRACTED 2026-04-26]` Thread B 跨樣本 whitelist filter 已撤回 → 詳見 [撤回宣告](reports/validated/2026/04/20260426_Thread_B_Retraction_Whitelist_Cross_Sample_01.md) | ~~S4 子 filter scheme~~ characterization-only 保留 | 1-2 天 | 本週 Thread B 伏筆 |

### Phase 2 研究方向（優先序）

1. **Phase 2 方向 A+D**：Normal Methylation Reference + CN/Purity-aware correction
   - **Phase A-D 程式碼已完成（2026-04-13）**：Normal BAM integration, Sample ASM, Somatic HP ASM, LOH BED annotation, Cross-region subclone analysis
   - HCC1395 全基因體驗證通過：Sample ASM 97.3% sig, Normal Baseline 100% valid, LOH concordance 94.1%, 4-group subclone
   - [驗證報告](experiments/validated/2026/04/20260413_Phase_BCD_Dual_BAM_Validation_01.md)
   - 待執行：7 samples 全量驗證、haplotag 重跑後的 Phase 2A 正式分析

2. **Phase 2 方向 B**：Gene-level / mechanism-level evidence integration
   - 目前狀態：待 Phase 2A 全量分析完成後啟動

3. **Phase 2 方向 C**：CpG 功能分層與智慧選點
   - 目前狀態：待規劃

### Self-Phasing 修正後重跑

- PON-only phasing 已驗證：LOH.bed 不變、somatic bias 消除、N50 +99.7%、phased rate +23.6pp
- **P0-3 LOH.bed 生成機制已確認（2026-04-11）**：LOH.bed 使用 VCF AF/VAF（`PhasingGraph.cpp:1817`，VAF >= 0.8 → HOM），ISM hp_ratio 使用 BAM HP tags（`ReadParser.cpp:123`）。兩套系統使用不同數據源，解釋 Jaccard=1.0 與 62% ISM LOH 消失的矛盾。
- 待執行：haplotag + ISM 全量重跑（7 samples × paired + TO）
- 重跑後可啟動：Phase 2A normal methylation reference baseline

### 🚨 Self-Phasing V5 Provenance Audit（2026-05-05 新發現，**P0**）

**CRITICAL**：[V5 Data Provenance Audit](reports/validated/2026/05/20260505_self_phasing_V5_data_provenance_audit_01.md)

- PI 報告（4-29）所有 V5 數值（sanity 15/15、+13.3pp、AMB 17.5→8%）= **Pass 1 only**（ploidy bug 讓 Pass 2 從未觸發）
- d0bcd8c (4-30) 已修 ploidy bug；938f0df (4-30 cherry-pick) threshold 0.95 → 0.9
- 4-30/5-01 已重跑：threshold_compare/{baseline_09 (purity=0.977), v5_flag (purity=0.984)} + v5_flag_force_path2only (purity=0.984)
- **P0 待辦 T1**：對 4-30 BAM 跑 ISM benchmark + sanity check 對齊 PI §6.4/§6.5 數字 → 三種情境（A 結論強化 / B Pass 2 冗餘 / C 需修訂）
- **P1 待辦 T2-T6**：trace 4-12 BAM 重產判定、7 樣本擴展、PI 報告 errata、manifest 加 binary_commit_hash、0.6 purity simulation 重做

### Phase 1A 已鎖定

- 最優方案：`methyl+context`，paired-pure multi-bio external validation delta F1=+0.0112, CI=[+0.0044,+0.0188]
- TO 模式下甲基化增益為負（delta F1=-0.0206）
- LOH feature 對 read-level 分類無用（LOH 是 region-level）
- 已正式關閉 TO single-sample post-hoc germline FP 識別

### LOH 雙定義與特徵探索全面關閉（2026-04-06 確定性結論）

- SEQC2 外部驗證 Jaccard=0.928：LOH.bed 與 FDA 金標準幾乎完全吻合
- LOH 區域 10/10 filter 策略全失敗（FP rate=0.239 < Non-LOH 0.338，LOH 是 TP-enriched）
- Non-LOH max AUC<0.58，多特徵 Voting AUC=0.577
- cnLOH 表面 0.587 是 Simpson's Paradox（per-sample mean=0.50）
- QS mode-aware 已實作：TO 模式停用 LOH penalty 與 verify bonus
- **結論：ISM 價值在 read-level epigenetic characterization，非 variant filter**

### Coverage_Multiple GC 校正與甲基化-CN 驗證 — 全 NO-GO（2026-04-11 確定性結論）

- **TO Pipeline 數據**：ClairS-TO → LongPhase-TO → ISM（TP=28,383, FP=11,830, F1=0.7127）
- **GC-Content 校正**：delta-r = -0.0002（Go/No-Go 門檻 ≥0.03）；ONT 5kHz GC bias 極小（98.7% regions 變化<5%）
- **甲基化-CN 相關**：所有 HP-free 特徵 residualized |r| < 0.07（甲基化對 CN 完全無感）
- **HPFineNGroups-CN**：raw r=0.495 → residualized r=0.160（68% 是 NumReads confound）
- **KDE auto-estimation 確認正確**：CN 分類準確度 6.2%→43.8%（已實作，`--expected-coverage` CLI 可用）
- **Coverage_Multiple 作為 CN proxy 已足夠**：TO r=0.827，Paired r=0.831
- 腳本：`scripts/analysis/gc_correction_to_validation.py`, `scripts/analysis/methylation_cn_validation.py`
- 圖表：Figures 28-33（`research/seqc2_cnv_stratification/figures/`）
- **結論：GC 校正不需實作；甲基化無法驗證 CN；Coverage_Multiple 現有精度已足夠**

### SEQC2 CNV 分層觀察 — CNV zone-aware filter 關閉（2026-04-10 確定性結論）

- **SEQC2 正交驗證 CNV truth set**（6 callers × 21 replicates × 3 technologies）用於 HCC1395 分層
- **Coverage_Multiple 驗證**：與 SEQC2 真實 CN Pearson r=0.831，ISM read-depth 代理可信
- **HCC1395 單樣本假象**：Gain+LOH zone AUC=0.782（AlleleDelta），但跨 7 樣本驗證 mean AUC ≤ 0.641
- **Simpson's Paradox 排除**：CNV 非 Quality_Score pooling 問題根源（分層後 QS diff=-0.042，pooling 反而有利）
- **Zone 排除策略全不可行**：所有排除策略 trade-off 低於 break-even（如排除 CN_Loss 移除 45% FP 但損失 11% TP）
- **Gain+LOH 根因**：CN=3+LOH 環境造成最高 FP rate 12.9%（allele imbalance 誤導 caller），FP rate 隨 CN 增加反而下降
- **15 張圖表 + 5 TSV 數據檔**：報告見 `docs/experiments/in_progress/2026/04/20260409_SEQC2_CNV分層觀察_01.md`
- **結論：CNV 不是特徵空間耗盡的根因；zone-aware filter 無可行策略，正式關閉**

### R1-R5 特徵設計研究（2026-04-07 確認）

- **R1**: CramersV 93% 為零 = 2×2 框架缺陷；HPMergedDelta 多群時反向；HPFineNGroups 已克服（AUC +0.125）
- **R2**: Excess groups 概念有效（跨子集統一 +0.059）但子集內無新信號，不需修改 C++
- **R3**: 結構清楚子集 AUC 反而下降 → **確認是 identifiability 問題而非特徵設計問題**
- **R4**: HPFineNGroups 新 canonical filter **NG=4 + AF<0.4 + NR≥80** → TP rate **92.81%**（舊 filter N≥4+NR≥80 為 89.12%，ΔTP +3.7pp；HCC1954 +21pp 挽救），F pilot 2026-04-18 驗證；低 AF (0.1-0.2) 信號最強（+50pp）
- **R5**: PairwiseMeanDist 與 HPFineN 正交（Spearman=0.07），微弱獨立增量
- **HPFineNGroups 確認為 somatic heterogeneity 標記** — 7/7 一致，residualized AUC=0.617，不能用於 filter 但有明確生物學價值

### Option C 雙路測試 NEGATIVE（2026-04-07 確定性結論）

- **架構確認**：cluster_labels 已是 HP-free（甲基化 hierarchical clustering），ClusterPermanova 被 passed_gate 擋住（TO 22% 有效）
- ClusterPermanovaF AUC=0.512（n=90,572）— 純甲基化 cluster 品質完全隨機
- HP-free 5 features combo AUC=0.564 vs HP-dependent combo AUC=0.598 vs 全部 AUC=0.607
- HP-free 特徵僅增加 +0.009 AUC — 無實用價值
- **結論：所有區分力來自 HP tags，純甲基化 clustering 無法突破 identifiability problem**
- C++ 代碼修改取消，確認正確策略為 PON-only phasing 上游修正

### O9 FN 特徵觀察 NO-GO（2026-04-08 確定性結論）

- 7 samples × 2 modes (Paired+TO)，122,790 FN regions 完整 ISM 執行
- **HP-free 甲基化特徵全 AUC < 0.53** — 純甲基化空間 FN 與 TP 不可區分
- 最強信號 LabelAllelePermanovaF=0.664（Paired）/ 0.636（TO）是 AF 代理非甲基化
- **TO Quality_Score AUC=0.338（嚴重反轉）** — FN 的 QS 比 TP 高，Cohen's d=-0.671
- FN VerificationClass: 55% Noise, 23% Weak, 18% Strong — 與 TP 分布相似
- **NO-GO：ISM 無法 rescue FN，甲基化空間 FN≡TP**

### TO-pure 獨立建模 NEGATIVE（2026-04-08 確定性結論）

- LOSO 7-fold：HP-free AUC=0.53、All-ISM=0.60-0.64、Caller-only=0.63、ISM+Caller=0.66
- caller_af (AUC=0.654) 單獨超越全部 ISM 特徵組合
- ISM 在 Caller 之上僅增加 +0.003（LR）~ +0.030（RF）
- HP-free ISM 完全隨機（AUC~0.53），HP-dependent +0.07~0.10 但可能循環
- **TO 模式 ISM 在 HCC1395 LOSO 上 caller_af 為主要判別器**（2026-04-22 修正：此結論限 HCC1395 TO 單樣本；跨樣本 TO 泛化性尚未驗證 — archive 5 樣本 TO 於 `research/tpfp_loh_af_kde_discrimination/08_archive_to_crosssample.md` 彙整觀察，HPFineNGroups≥4 於 5/6 樣本 TP%≥69%，但 LOH/CN 欄位缺失，需重跑 ISM 才能完整驗證 filter 跨樣本）

### 獨立分析完成（2026-04-11）

- **PON 跨樣本移除率驗證 ✓**：7 樣本 raw PON rate 95.16-98.81%（mean 97.77%），refined FP-level 全 > 98%。結論穩定度 3→4。H2009 最低（95.16%）因高突變負荷非 PON 失效。
- **H2009 負向根因 ✓**：Paired FP rate 僅 0.06%（86/132,994），改進天花板 +0.0004 F1。76.7% FP 在 LOH、89.5% Noise class。根因=caller 已近乎完美，ISM 只能誤傷 TP。甲基化特徵區分力反而優於平均。

### LOH / CN / AF 結論速查

完整三維度統合見 [07_LOH_CN_AF_研究總整理](reports/research_landscape/07_LOH_CN_AF_研究總整理.md)。

| 維度 | Filter 方向 | Characterization 方向 |
|------|-----------|---------------------|
| **LOH** | ❌ 全面關閉（10/10 策略失敗） | ✅ Subclone AF×Methylation 雙模式確認 |
| **CN** | ❌ Zone 排除全不可行 | ✅ Coverage_Multiple r=0.831 代理可信 |
| **AF** | ⚠️ 唯一有效信號但來自 caller | ✅ AF gradient 預測 NGroups |

### Zone-Aware Confidence Framework（2026-04-17 完整驗證）

- **構想文件**：[Zone-Aware Confidence Framework](concepts/2026/04/20260417_Zone_Aware_Confidence_Framework_01.md)
- **5 Zone 定義**：Z1 LOH Subclonal Active, Z2 High Somatic Hetero, Z3 Complete LOH, Z4 Normal Diploid, Z5 CN Gain Low Diversity
- **H1/H3 驗證**：[報告](../research/zone_aware_validation/20260417_H1_H3_Zone_Validation_Report_01.md)
- **Z3 內部特徵探索（2026-04-18 NEGATIVE）**：[報告](experiments/in_progress/2026/04/20260418_Z3_Internal_Feature_Exploration_01.md)
  - Step 1 Z3 內 12 特徵 × 7 樣本 AUC：無特徵在 ≥3 樣本達 |AUC|≥0.60
  - Step 2.5 AF∈[0.4,0.6] × CN × NGroups 分層：僅 HCC1954 (1/7) 符合 germline pattern（TP rate=0.146）
  - Step 3 HCC1954 vs HCC1395 機制對比：HCC1954 Z3 FP 集中 chr5/8/17（HER2/MYC amplicon），FP NumReads=55 vs TP=37（p=4.7e-9）→ CNV amplicon artifact 驅動；HCC1395 均勻分佈
  - **結論**：Z3 內無跨樣本二階區分特徵；HCC1954 例外已由 F pilot canonical filter 覆蓋
- **Z3 × HCC1954 Amplicon Blacklist Pilot（2026-04-19 CONDITIONAL）**：[設計](experiments/in_progress/2026/04/20260419_Z3_amplicon_blacklist_pilot_design_01.md) · [結果](experiments/in_progress/2026/04/20260419_Z3_amplicon_blacklist_pilot_result_01.md)
  - S2 whole-chr5/8/17 ∩ Z3：HCC1954 ΔF1=+0.0065（ceiling +0.0075 的 87%）
  - 其他 6 樣本 mean ΔF1=−0.0044（5/6 hurt）→ 非 global canonical filter
  - Circularity guard（S3 CovM 95%ile non-Z3 baseline）信號過弱 Δ≈+0.0002
  - **結論**：HCC1954-local CONDITIONAL；不納入預設 filter；Zone-Aware Framework 定位不變
- **CovM Baseline 準確度驗證（2026-04-20 MIXED；2026-04-19 H-CN1 重審降為 CONDITIONAL）**：[報告](experiments/in_progress/2026/04/20260420_CovM_Baseline_Accuracy_Validation_01.md) · [方法審計](methodology/20260419_KDE_expected_coverage_audit_01.md)
  - **H-CN1 CONDITIONAL**：master dataset（2026-03-30）由 KDE commit 8d0a0c8（2026-04-13）**之前**的 binary 產出 → 全部 75.0 是 stale-binary artifact，非「KDE 未啟用」亦非 KDE 邏輯缺陷；per-sample baseline bias 需以現行 binary 重跑才能確立
  - **H-CN2 POSITIVE（觀察仍成立但需以新 master 重跑量化）**：HCC1395 SEQC2 benchmark — Gain recall 14.6%（CN=3 僅 0.15% 標為 Gain）vs Loss recall 86.9%
  - **H-CN3 NEGATIVE**：HCC1395 oracle truth_cn≥3.5 ∩ Z3 ΔF1=+0.0011；Variant A 數學證明 percentile filter 對 per-sample re-centering 免疫
  - **結論**：「cov 需要修正」仍為真（stale-binary 導致），但**修正後不會救 Z3 pilot 跨樣本失敗**；後者是真實生物學 CNV 異質性
  - **2026-04-19 C++ 改動**（commits `374fad4` + `12d9b3e`）：KDE fallback WARN + `Diploid_Coverage_Used` TSV column + DiploidEstimate struct（避免 value==75.0 sentinel 誤判）；evidence ledger H_KDE_001
  - **後續動作（Blocker）**：以現行 binary 重跑 master dataset（7 樣本 × 2 modes，~4-6 hr），取得 per-sample KDE baseline 後重量化 H-CN1/H-CN2；所有 2026-04-19 前的 cross-sample CovM 分析需標註此前提
  - **KDE 新算法覆蓋狀態（2026-04-23 澄清更新）**：
    - **paired_full 7/7 全部用新 KDE binary 直接重跑完成**（含 HCC1954：`kde_rerun_B_14combos/HCC1954_paired_full_tp/` 2026-04-21 產出，Diploid_Coverage_Used=61× 全 17,909 rows 一致；更正先前「post-hoc 除法」誤記）
    - **TO 1/7**（僅 HCC1395 Phase 1）；其他 6 TO 樣本需等 Archive TO rerun（P0 項 ~10 hr parallel）與 COLO829 TO 背景 pipeline 完成
    - 詳見 `docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md` §5.0.1
  - H1 CONDITIONAL → Z1b 放寬後 TO 4.6% 覆蓋率、TP rate 0.965（7 變體 Pareto 最佳）
  - H3 PARTIAL：Paired 7/7 ≥ 89.1% 確認；TO 6/7 significant 但絕對值 ~72%
  - TO zone TP rate 範圍 0.61-0.94
- **QS 模擬 ❌ NEGATIVE**：[報告](../research/zone_aware_validation/20260417_QS_Simulation_Report_01.md)
  - 5 configs × 21 thresholds × 7 samples，max delta F1=+0.001
  - **根因：TO QS AUC=0.497 隨機，zone delta 無法修正**
  - QS 調整路線和 C++ 整合**暫停**
- **ClairS-TO Verdict Characterization Pilot（2026-04-20 NEGATIVE on F1 / POSITIVE on calibration）**：[報告](experiments/in_progress/2026/04/20260420_ClairS_TO_Verdict_Characterization_Pilot_01.md) · [外部 CN 工具 survey](references/2026/04/20260420_external_CN_tools_survey_01.md)
  - scope：HCC1395 subsample t20_n30（purity ≈ 0.40），因其他 6 樣本缺 CNA loci 無 Verdict 標籤
  - H-V1 POSITIVE：Verdict_Germline FP rate = 96.1%（3,463/3,602, Wilson 95% CI 0.955–0.967）
  - H-V2 POSITIVE：Verdict_Somatic TP rate on PASS = 91.8%；Verdict_SubclonalSomatic TP rate = 94.9%
  - **F1 直接增益 = 0**：Verdict_Germline 100% 落在 LowQual（從未出現在 PASS）→ S1 「remove Verdict_Germline from PASS」ΔF1 = +0.0000
  - 最激進 S2「只保留 PASS ∩ Verdict_Somatic∪Subclonal」 ΔF1 = −0.2007（recall 崩塌）
  - 根因：Verdict 與 LowQual 共用同一套 ASCAT / 二項式資訊；PASS 95.4% FP 落在 no_Verdict 子集，Verdict 無決策權
  - **結論**：不升級 ClairS-TO 內建 Verdict 為 production filter；Verdict 標籤保留作 per-variant annotation；主升級路徑改為 **Wakhan（haplotype-specific CN）+ SAVANA（SV/CNA）external CN pilot**
- **結論：Zone-Aware 價值確認僅在 characterization annotation，不在 F1 改進**

### 待完成項目

- ~~Zone Z1 定義放寬~~ — ✅ 完成，Z1b 為最佳（TO 4.6%, TP rate 0.965）
- ~~Zone-Aware QS 調整模擬~~ — ❌ NEGATIVE，max delta F1=+0.001
- haplotag + ISM 全量重跑（7 samples × paired + TO）— PON-only phasing 後
- Phase 2A Normal Methylation Reference baseline — 依賴重跑數據
- LOH Subclone 更精細 AF-bin 分析 — 現有數據即可
- 7 樣本全量 Phase 2A 驗證 — 依賴重跑數據
- Platform normalization（5kHz / DORADO / PAO / Google）

## 4. 阻塞與風險

1. `archive/deep/` 為 immutable 快照，保留歷史失效連結（不回寫）。
2. `purity-aware` / subsample 結果易受 tumor-normal 組織甲基差異混淆，暫不作主證據。
3. `pileup` 與 `full model` 的 TO 輸入不可混用，報告必須明確標示來源。
4. TO borderline rescue 必須使用 candidate-specific `lost_tp / removed_fp` run。

## 5. 每次研究啟動必查

1. `docs/CURRENT_FOCUS.md`（本檔案）
2. `docs/experiments/INDEX.md`
3. `docs/README.md`
4. `docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md`（高密度壓縮上下文、重要數據、任務順序、待決策矩陣）
5. `docs/concepts/2026/04/20260409_研究構想總索引_01.md`（研究大圖景、發展樹、理論基礎）
6. `docs/reports/research_landscape/00_INDEX.md`（需深度理解時）

## 2026-05-10 — image-gen + image-vision-check skills shipped (Phase 1)

Two new skills shipped:
- `InterSubMod/.claude/skills/image-gen/` — dual-track image generation (AI via codex `$imagegen` + cairo program rendering)
- `InterSubMod/.claude/skills/image-vision-check/` — 6-dim vision audit via Claude Read

End-to-end demo validated at `InterSubMod/examples/phase1_demo/demo_topic/`. Ready to plan Phase 2 (html-preview skill).

Spec: `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md`
Phase 1 plan: `InterSubMod/docs/plans/2026/05/20260510_Phase1_image_gen_vision_check_implementation_01.md`

## 6. 已完成研究概覽

所有已完成研究的詳細記錄已封存至 `docs/archive/2026/03/20260405_current_focus_completed_items.md`。

**關鍵里程碑**：
- **LOH 雙定義與特徵全面關閉**（166 圖表 × 16 判定）：SEQC2 Jaccard=0.928 驗證 LOH.bed 可信；10/10 filter FAIL；Non-LOH max AUC<0.58；ISM 定位轉向 characterization
- **O1-O10 系統性觀察**（82 圖表）：TO 無單一特徵 AUC>0.58；LOH penalty 是 QS 根因；Paired/TO 根本不同
- **O11-O13 甲基化假說**：三維度全 NEGATIVE（within-group heterogeneity / LOH scenario / cross-region correlation）
- **G1-G7 VCF 特徵**（48 圖表）：60+ 特徵全 AUC<0.64；TO germline FP post-hoc 識別正式關閉
- **Read-level 鑑別**：LOSO AUC=0.721（首次>0.70）但安全約束 FP removal=0%
- **ASM 定量**：POSITIVE — 32-66% SNV 位點有 ASM；ISM PERMANOVA 唯一捕捉 entropy imbalance
- **Phase 1A**：paired-pure delta F1=+0.0112 已鎖定；TO 模式負增益
- **Self-phasing 因果鏈**：CONFIRMED — 62% ISM HP_Ratio LOH 消失（非 LOH.bed；LOH.bed Jaccard=1.0 不變）、somatic bias 17.3:1
- **SEQC2 CNV 分層觀察**（15 圖表 + 5 TSV）：Coverage_Multiple r=0.831 驗證可信；跨樣本 AUC ≤ 0.641；zone 排除全不可行；CNV 非特徵耗盡根因
