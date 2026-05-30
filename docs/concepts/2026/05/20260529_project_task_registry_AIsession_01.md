---
title: InterSubMod 全專案目標與任務登錄表（AI-session handoff）
generated: 2026-05-29
source: 7-thread multi-agent survey (run wf_ab45fe26-785), verifier verdict PASS
status: in_progress
scope: 全專案所有研究主軸（comprehensive, per task-type D handoff）
companion_html: InterSubMod/docs/reports/in_progress/2026/05/20260529_project_goal_landscape_dashboard_01.standalone.html
machine_json: InterSubMod/docs/concepts/2026/05/20260529_project_tasks_01.json
---

# InterSubMod — 全專案目標與任務登錄表

> **給接手的 AI session**：本檔是整個 InterSubMod 專案（不只 HKU）的目標樹 + 任務 DAG 快照，由 7 條主軸的 fresh-context 並行盤點 + 合成 + 對抗式驗證產生。讀完即可掌握「現在在哪、什麼緊急、下一步、什麼不能碰」。互動版見 companion HTML；機器可讀版見 machine_json。

## 0. 如何使用本檔

1. 先讀 §1 即時狀態 + §2 大局（6 目標 + 戰略 pivot）。
2. 要動手 → 看 §3 critical path + §4 優先序 Top-15。
3. 找特定任務 → §5 全任務（依主軸分組）；每個 `task_id` 與 machine_json 一一對應。
4. 動之前 → 看 §6 blockers + §8 護欄（勿重啟已 concluded NEGATIVE）。
5. governance：C++ 改動走 `/methodology-audit`+compile Hard Gate；研究 NO-GO 判定 / 刪檔需用戶確認；cycle 由用戶觸發不自動啟。

## 1. 即時狀態（2026-05-29 親自 df / ps / git 驗證）

| 項目 | 狀態 |
|------|------|
| 磁碟 /big7_disk | 🔴 **97%** used（剩 1.4T/42T）|
| 磁碟 /big8_disk | 🔴 **98%** used（剩 689G/27T）|
| 磁碟 / (incl /tmp) | 🟢 **28%** used（剩 703G/1008G）|
| ASM survey | ⏳ PID 411460 — scripts/19_full_tp_asm_batch.sh (ASM genome survey, actively writing) |
| HKU 交付 git | 📦 untracked (?? — 3 deliverables + 22 figs + 10 scripts not yet committed) |
| 分支 | refactor/phase1-safety @ 274f152 |

> ⚠️ `/big7_disk` 97% + `/big8_disk` 98% — 啟動任何長計算前先 `df -h` 並留 margin；ASM survey (PID 411460) 完成前勿並行其他大 pipeline。

## 2. 大局 — 6 個最終目標 + 戰略 pivot

**關鍵轉折（2026-05-13 reframe + 2026-05-20 LOSO）**：

- G1/G2/G4 → **characterization-only**（不再以 F1 為成功標準）。
- G3 → **blocked on G1**（region-level 已 CONDITIONAL NEGATIVE；需 G1 per-read epigenotype 才續）。
- **G5（甲基增強 F1 filter）路線實質關閉** — LOSO 證實 Cycle1 ΔF1=+0.02236 是 100% sample-level circularity（LOSO mean −0.00004，0/5 positive）。
- **G6 論文主軸 pivot 到 LOH-constrained phasing**，把 LOSO NEGATIVE 當成嚴謹 methods finding，ASM characterization 當支撐證據。

### G1 — 目標 1 — per-CpG 甲基位點多標籤關聯性評分：分析每個 CpG 位點甲基化與 ALT/REF、HP1/HP2 的統計關聯強度（per-CpG methylation multi-label association scoring vs ALT/REF + HP1/HP2）。2026-05-13 重新定位為僅做特徵描述（characterization-only）。為 G3 的前置條件。
- **狀態**：ACTIVE，僅做特徵描述 — 單一位點 BRCA2/ZAR1L ASM 已驗證 ⭐3（POSITIVE 但幅度有限，單樣本 HCC1395）；全基因組 ASM 普查進行中（22 條中 chr1-4 的 TP 已完成）；待修正 5mC/5hmC collapse、跨樣本提升與 Level1 dup-bug 修正後才能超越 tier 3。
- **貢獻主軸**：ASM, Vision/Queue, Phasing

### G2 — 目標 2 — clone 結構分析：單位點的 sub-clone 分析 + 跨位點的 clone 共演化分析（clone structure analysis: per-site subclone + cross-site clone co-evolution）。2026-05-13 重新定位為僅做特徵描述（4 群 subclone 描述）；排入 Phase B clone preprint。
- **狀態**：DEFERRED 至 Phase B — HPFineNGroups subclone marker 維持 ⭐3 且依賴 pipeline；phase_block_3d（H013-H015）為目標專案，尚未建立骨架。
- **貢獻主軸**：Phasing, Vision/Queue, Guardrails

### G3 — 目標 3 — 二次打擊與事件順序推論：依賴目標 1/2，結合 LOH/CNV/HP 推論打擊順序（second-hit & event-order inference; depends on G1/G2, combines LOH/CNV/HP）。
- **狀態**：BLOCKED 於 G1 — region-level 方法已證實為 CONDITIONAL NEGATIVE（P4 pilot 2026-04-17，|Δ|=0.043<<0.15）。正確路徑 = 等待 G1 的 per-read epigenotype，再做 phased haplotype × methylation sub-pattern。在沒有 C3（G1 per-read epigenotype 資料）前不要重啟。
- **貢獻主軸**：Vision/Queue, Guardrails, ASM

### G4 — 目標 4 — TO normal 資訊補強：在無配對 normal 的情況下，利用甲基/germline SNP/HP 估計 normal 背景（TO normal-info augmentation without paired normal）。2026-05-13 重新定位為僅做特徵描述（R1-Global 在 F1 上已為 NEGATIVE）。
- **狀態**：ACTIVE，僅做特徵描述 — Normal BAM 跨樣本列為 Tier 4.1 reactive（45%->70%）；HCC1395 normal-BAM Sample-ASM vs Normal-ASM pilot 尚未開始（BAM 複製中）。
- **貢獻主軸**：Infra/Data, Vision/Queue, Guardrails

### G5 — 目標 5 — 整合 evidence panel 提升 F1：補充 caller 輸出，保留更多 TP、過濾 FP（integrated read-level evidence panel to raise somatic-SNV F1 above ClairS-TO baseline F1=0.7166）。
- **狀態**：F1-FILTER 路線實質關閉 — 2026-05-20 LOSO 證實 Cycle 1 的 +0.02236 完全是 sample-level circularity（LOSO 平均 -0.00004，0/5 為正）；ISM 為殘留的 caller_af proxy。僅能作為嚴謹的 NEGATIVE 方法學發現留存。重啟需 C2（非 LR 框架）或 C1（真正的 low-F1 跨樣本 panel）。
- **貢獻主軸**：Phase2 FP-filter, Vision/Queue, Guardrails, Phasing

### G6-PAPER — 論文定位目標 — 將 InterSubMod 定位為一項『read-level 表觀遺傳 + haplotype 脈絡用於 somatic 變異判讀』的研究（characterization，而非 variant filter），以 LOH-constrained phasing signature 為候選主軸，以記錄完整的 LOSO sample-level-circularity NEGATIVE 發現作為嚴謹的方法學貢獻，並以 ASM 特徵描述（BRCA2/ZAR1L + 全基因組普查）作為佐證。策略順序 E->A+D->B->C（2026-03-24 確認）。
- **狀態**：ACTIVE PRIMARY — 主軸已轉向 phasing（LOH-constrained，可由 LongPhase+LOH bed 重現，不需 methylation）；發表前需正式統計 + paired phasing-vs-methylation 拆解 + 全基因組延伸；Phase A 框架論文（thread_d）受 V6 production tag 限制。
- **貢獻主軸**：Phasing, ASM, Phase2 FP-filter, HKU handoff, Vision/Queue, Infra/Data, Guardrails

## 3. Critical Path（通往最近高價值目標的最短鏈）

1. **維護操作護欄：>100MB 的 pipeline 須設 TMPDIR=/big7_disk + 明確指定 -o + df 預檢；讀取 merged file 時須 archive caller_af + 過濾 phase1_new；新伺服器須 mount|grep tmp** (`infra-maintain-guardrails`)
2. **對 6/6 TP-gap（中位數 +0.37）做正式的 Wilcoxon signed-rank，為 LOH-constrained phasing 發現提供可發表的顯著性檢定** (`phasing-formal-wilcoxon-tpgap`)
3. **執行 --germline-hp-only=true negative control，確認它能消除 Inner-vs-Outer 的 NG=2 TP gap（如預測；Phase 1 實驗可重用）** (`phasing-negative-control-germline-hp-only`)
4. **以校正後的 MAX-collapse 5mC+5hmC 進行完整 39,447-TP 重跑（取代 PROVISIONAL 折半的 beta；幅度約 2x 更大）— genome_survey_v2 批次 PID 411460 目前跑到 chr5** (`asm-corrected-collapse-rerun`)
5. **完成 genome_survey_v2 雙軸 TP-vs-FP ALLELE-axis 調查直到 chr22（FP 已完成 5,149；TP chr1-4 已完成，chr5-22 待處理）** (`asm-genome-survey-v2-complete`)
6. **從源頭修正 Level1 重複輸出 bug（MSA C++ 將 p 與 1-p 輸出為每個 read-CpG 2 列 -> 2.6x 膨脹），或確認 MAX-collapse 完全抵銷，才能回報任何量化的 beta/n/p** (`asm-fix-level1-dup-bug`)
7. **用校正後的 pipeline 重建 ASM negative control，採用 n>=50-100 個配對對照（TVAF + 局部 CpG 密度 + coverage 配對）（目前 n=5，僅 3 個有效）** (`asm-rebuild-negative-control`)
8. **待 corrected-collapse 量值定案後，將 ASM 實驗文件 + evidence_ledger 更新為最終版（去除 PROVISIONAL 標記）；透過 conclude-research 把 in_progress -> validated** (`asm-finalize-experiment-doc`)
9. **正式 /conclude-research NEGATIVE 撰寫 + 封存 methyl_augmented_filter_phase2（已驗證 manifest 仍為 status:initiated、conclusion:null；pre-reg 指出 all-H-NEGATIVE -> 結論 NEGATIVE）。Hard Gate：NO-GO 需使用者確認** (`methylfilter-conclude-negative-archive`)
10. **草擬統一的論文大綱，將 phasing 主軸 + ASM 特徵刻畫 + LOSO sample-level-NEGATIVE 方法學發現整合成一條 read-level-epigenetic-context 敘事** (`paper-unify-three-axes-outline`)

## 4. 推薦優先序 Top-15

> heuristic：value × confidence ÷ effort，但尊重依賴與 blocker。你可自行重排 —— 這是建議不是命令。

| # | 任務 | 主軸 | 價值 | 狀態 | 工時 | 為何這個排名 |
|---|------|------|------|------|------|-------------|
| 1 | 維護操作護欄：>100MB 的 pipeline 須設 TMPDIR=/big7_disk + 明確指定 -o + df 預檢；讀取 merged file 時須 archive caller_af + 過濾 phase1_new；新伺服器須 mount|grep tmp<br>`infra-maintain-guardrails` | Infra/Data | 🟥高 | 🔵進行 | 持續進行 | 零成本的持續護欄；/tmp 事故曾讓全部 33 個帳號當機，故它必須先於每一項已在跑（ASM batch）或即將跑（Archive TO）的 long compute。已 in_progress — 維持強制執行。 |
| 2 | 以校正後的 MAX-collapse 5mC+5hmC 進行完整 39,447-TP 重跑（取代 PROVISIONAL 折半的 beta；幅度約 2x 更大）— genome_survey_v2 批次 PID 411460 目前跑到 chr5<br>`asm-corrected-collapse-rerun` | ASM | 🟥高 | 🔵進行 | 執行中（數小時） | 已確認正在執行中（PID 514459 chr5）。對 G1/G6 高價值，且已在進行中 — 讓它跑完並監控；它能解鎖確定性的 ASM 量級並移除錯誤的暫時歸因。 |
| 3 | git commit v4 HKU handoff（research/hku_collaboration/ + docs/concepts/2026/05 + docs/references/2026/05 + docs/reports/in_progress）— 非 C++，非 Hard Gate<br>`hku-commit-v4-deliverables` | HKU handoff | 🟥高 | ⚪todo | 30 分鐘 | 30min、高價值、零依賴、非 Hard-Gate。已驗證完全未被追蹤 — 立即 commit 可保護約 24 張圖 + 12 個腳本 + 3 份文件不致遺失。 |
| 4 | 讓 HKU standalone HTML 適合 email 寄送 — 將 18 個外部相對路徑的 <img> 轉為內嵌 base64（已驗證 0 個 base64 / 18 個外部 img）<br>`hku-html-email-safe-base64` | HKU handoff | 🟥高 | ✅done | 1-2 小時 | 1-2h、高價值、無依賴。已驗證 0 個 base64 / 18 張外部 img — 此「standalone」HTML 一旦寄出就會壞掉；這是 D-handoff 關鍵路徑上成本最低的待清障礙。 |
| 5 | 對 6/6 TP-gap（中位數 +0.37）做正式的 Wilcoxon signed-rank，為 LOH-constrained phasing 發現提供可發表的顯著性檢定<br>`phasing-formal-wilcoxon-tpgap` | Phasing | 🟥高 | ✅done | 1 小時 | 1h、高信心、高價值、無依賴。將 6/6 的描述性 TP-gap 轉成候選論文 MAIN 主軸（G6）可發表的顯著性檢定。 |
| 6 | 執行 --germline-hp-only=true negative control，確認它能消除 Inner-vs-Outer 的 NG=2 TP gap（如預測；Phase 1 實驗可重用）<br>`phasing-negative-control-germline-hp-only` | Phasing | 🟥高 | ⚪todo | 2-4 小時 | 2-4h、高信心、高價值。此為預測中的 negative control，可確認 phasing 主軸為真（非 artifact）；屬可重用的 Phase-1 實驗。 |
| 7 | 正式 /conclude-research NEGATIVE 撰寫 + 封存 methyl_augmented_filter_phase2（已驗證 manifest 仍為 status:initiated、conclusion:null；pre-reg 指出 all-H-NEGATIVE -> 結論 NEGATIVE）。Hard Gate：NO-GO 需使用者確認<br>`methylfilter-conclude-negative-archive` | Phase2 FP-filter | 🟥高 | ⚪todo | 2-4h | 已驗證 manifest 仍為 status:initiated / conclusion:null。將其收尾（Hard Gate、需 user 確認）可把嚴謹的 LOSO 循環性發現轉為可引用的 NEGATIVE 方法學貢獻，並解鎖 Cycle 5 pivot。 |
| 8 | 對外交付前，將 F4 HCC1395 baseline HP counts（863K/266K/79K/23K/630）與 evidence_ledger 交叉核對（README 標註為憑記憶引用、未驗證）<br>`hku-verify-f4-hp-counts` | HKU handoff | 🟧中 | ✅done | 30 分鐘 | 對外交付前花 30min 驗證 memory 引用的數字；existing-artifacts-must-verify 規則所要求，成本低，位於 HKU 關鍵路徑上。 |
| 9 | 對外交付前，將 F5 合成 rng.beta ISM heatmap placeholder 替換為真實 HCC1395 ISM region（候選 chr2:18,086,020 +-32kb）<br>`hku-replace-f5-synthetic-heatmap` | HKU handoff | 🟧中 | ✅done | 2-4 小時 | 把合成的 placeholder 資料交付給外部合作者是誠信風險；PI sign-off 前必須換成真實的 HCC1395 ISM。 |
| 10 | 完成 genome_survey_v2 雙軸 TP-vs-FP ALLELE-axis 調查直到 chr22（FP 已完成 5,149；TP chr1-4 已完成，chr5-22 待處理）<br>`asm-genome-survey-v2-complete` | ASM | 🟥高 | 🔵進行 | 執行中（數小時） | 已在執行中；完成 chr5-22 的 TP 即為當前 TP-vs-FP 可區分性檢定（G5 probe），並解鎖 strong-site 排名 + Level1 bug-fix 驗證。 |
| 11 | 對全部 25+ 項已收尾的 NEGATIVE/NO-GO 結果強制執行 REOPEN GATE：任何 reopen 都必須引用具體的 C1（新數據）/ C2（新方法）/ C3（新前置條件），並在開始計算前記錄於 research/<topic>/00_INDEX.md 的 Pre-registration<br>`guardrails-reopen-gate-enforce` | Guardrails | 🟥高 | 🔵進行 | 持續進行 | 持續性 governance，保護價值高；確保 25+ 個已收尾的 NEGATIVE 不會在缺乏指名 C1/C2/C3 的情況下被重啟，保護後續所有 cycle。 |
| 12 | 從源頭修正 Level1 重複輸出 bug（MSA C++ 將 p 與 1-p 輸出為每個 read-CpG 2 列 -> 2.6x 膨脹），或確認 MAX-collapse 完全抵銷，才能回報任何量化的 beta/n/p<br>`asm-fix-level1-dup-bug` | ASM | 🟥高 | ⚪todo | 2-4h（C++ 變更：methodology-audit + 編譯 gate + 使用者確認） | 把關每一次量化的 ASM tier 升級；屬 C++ Hard Gate（methodology-audit + compile + user 確認），故應排在執行中的 survey 確認 MAX-collapse 是否已將其抵銷之後。 |
| 13 | Cycle 5 轉向決策，於 Path A（phase_block_3d/thread_d，建議）、Path B（chr8 zone gate）、Path C（low-F1 panel）之間抉擇。Hard Gate：宣告 LR filter NO-GO 需使用者確認；已驗證無 cycle5/ 目錄<br>`methylfilter-cycle5-pivot-decision` | Phase2 FP-filter | 🟥高 | 🔴blocked | decision（使用者觸發） | 由 user 觸發的 Hard Gate；推薦 Path A（phase_block_3d），可使已降權的 filter 路線對齊論文主軸 — 解鎖 phaseblock3d scaffolding。 |
| 14 | phase_block_3d 專案 init-research 鷹架建置（規劃於 5/23 V6 截止後）— H013/H014/H015 的目的地（G1/G2/G3 特性描述，Phase C ASM + 二次打擊）<br>`phaseblock3d-init-research` | Vision/Queue | 🟥高 | ⚪todo | 2-4 小時 | 一旦選定 pivot，即搭建 Path-A characterization 專案（G1/G2/G3）；中等工作量，對論文價值高。 |
| 15 | 候選 #4 / Archive TO 6 樣本 ISM 重跑，套用新 KDE + LOH.bed + germline-hp-only=off，將 LOH_Subtype/CovM_used/baseline_x 加入 master_extended（Path 2 約 4.5h 或全量約 10h）。USER 觸發的長計算；export TMPDIR + -o<br>`infra-archive-to-ism-rerun` | Infra/Data | 🟥高 | ⚪todo | 4.5-10h 平行 | 由 USER 觸發的長計算，是 unified master、paired phasing-vs-methylation 拆解、以及 HPFineN two-flag 的資料前置；但受 V6 tag 把關，故位於 production-decision 鏈之後。 |

## 5. 全任務清單（依主軸分組）

欄位：★=necessary（在通往某最終目標的關鍵路徑上）◉=observation-first（須先觀察/pilot 才能 commit）。

### Vision/Queue （15 任務）

| # | task_id | 任務 | 狀態 | 價值 | 信心 | 難度 | 工時 | 標記 | 依賴 |
|---|---------|------|------|------|------|------|------|------|------|
| 14 | `phaseblock3d-init-research` | phase_block_3d 專案 init-research 鷹架建置（規劃於 5/23 V6 截止後）— H013/H014/H015 的目的地（G1/G2/G3 特性描述，Phase C ASM + 二次打擊）<br><sub>Path A pivot 目標；為服務 G1/G2/G3 與論文主軸的特性描述專案建立鷹架。</sub> | ⚪todo | 🟥高 | medium | medium | 2-4 小時 | ★ | `methylfilter-cycle5-pivot-decision` |
|  | `paper-unify-three-axes-outline` | 草擬統一的論文大綱，將 phasing 主軸 + ASM 特徵刻畫 + LOSO sample-level-NEGATIVE 方法學發現整合成一條 read-level-epigenetic-context 敘事<br><sub>論文定位這個目標需要在三條可挽救的脈絡各自有可辯護的結論後，明確地將它們整合起來。</sub> | 💡idea | 🟥高 | medium | medium | 半天 | ★ | `phasing-formal-wilcoxon-tpgap`<br>`asm-finalize-experiment-doc`<br>`methylfilter-conclude-negative-archive` |
|  | `phaseblock3d-h013-cis-effect` | H013 X1 cis 效應：phase-block 內 CN 分層的甲基化高一致性 -> TP 富集（HCC1395_V6_TO，Wilcoxon p<0.05 + Cohen d>0.5）<br><sub>測試 cis 甲基化一致性是否為 phase block 內的 TP 訊號（G1/G2）。</sub> | ⚪todo | 🟧中 | medium | medium | Phase A W3-W8 | ◉ | `phaseblock3d-init-research` |
|  | `phaseblock3d-h014-boundary-artifact` | H014 X2 邊界假影：phase block 跨越大型 CN 區段 -> 甲基化不連續 -> FP 標記（HCC1395_V6_TO，Fisher exact p<0.05）<br><sub>測試 CN 邊界甲基化不連續是否為 FP 標記（G5 特性描述）。</sub> | ⚪todo | 🟧中 | medium | medium | Phase A W3-W8 | ◉ | `phaseblock3d-init-research` |
|  | `phaseblock3d-h015-zone-joint-score` | H015 X3 zone 分層聯合評分：每個 zone 套用 NG×CN tier×甲基化 rank-2 異質性 tau 過濾；pre-reg ΔF1>=+0.005 並做 residualized-AF 檢查（無 caller_af proxy）<br><sub>逐 zone 測試甲基化訊號是否在 Z-OCH/Z-GL 中增強；residualized-AF 檢查可防範 caller_af proxy guardrail。</sub> | ⚪todo | 🟧中 | medium | hard | Phase A W3-W8 | ◉ | `phaseblock3d-init-research` |
|  | `vision-candidate-caller-benchmark` | 候選 #3：caller 層級 benchmark（DeepVariant/Strelka2 對比 ClairS-TO）以驗證 HCC1954 caller 上限（待用戶核可，1-2 天）。必須先查 KB 取得外部工具 F1 口徑<br><sub>驗證 HCC1954 的低 F1 是否為 caller 上限所致；依 outside-claim 規則，外部工具 metrics 必須先經 KB 查核。</sub> | ⚪todo | 🟧中 | medium | medium | 1-2 天 | ◉ | — |
|  | `vision-candidate-external-cn-sv-pilot` | 候選 #2：外部 CN/SV pilot（在 HCC1395 上跑 Wakhan + SAVANA）以強化 thread_d §3（待用戶核准，6-8h）<br><sub>外部 CN/SV calls 強化 framework 論文的 CN 論點（論文目標）。</sub> | ⚪todo | 🟧中 | medium | medium | 6-8 小時 | ◉ | — |
|  | `vision-candidate-hpfinen-reverify-phase2b` | 候選 #1：Phase 2B HPFineN marker 重新驗證（master + flag=on）以釋出 R-SELFREF（待用戶核准，2-4h）<br><sub>釋出 subclone marker（G2）上的自我參照 blocker；以用戶核准為前提。</sub> | ⚪todo | 🟧中 | medium | easy | 2-4 小時 |  | `infra-archive-to-ism-rerun` |
|  | `vision-cycle3-h018-negative-control` | H018 H_C3_3 negative control：在 high-F1（>0.83）樣本上套用 caller-F1-headroom skip（H1437/H2009/HCC1954 V6_TO），pre-reg \|ΔF1\|<0.001。<br><sub>重用 cycle-2 函式的低成本配套 negative control；只有在 cycle 3 啟動時才有意義。</sub> | ⚪todo | 🟩低 | high | easy | 1-2 小時 |  | `vision-cycle3-write-planjson` |
|  | `vision-cycle3-resolve-mustfix` | 解決 cycle 3（H016）剩餘的 3 個 must-fix 項目（cycle_id 用 _validation 而非 _negative；ledger entry 缺 binary_commit/dataset_id/pre-reg link；H_C3_1 target-calc 釐清）— 已驗證該 cycle 目錄只有 state.json。<br><sub>冷卻期已過但 cycle 仍停擺；這些修起來成本低，但整個 cycle 被 panel block，因此除非找到 panel 否則優先度低。</sub> | 🔴blocked | 🟩低 | high | easy | 1-2h |  | — |
|  | `vision-cycle3-write-planjson` | Cycle 3（H016）P1 PLAN：撰寫 plan.json 以推進過 P0_REGISTER（針對 caller F1<0.80 子集的 caller-F1-headroom filter；target mean ΔF1>=+0.01）。USER 觸發的 cycle 工作。<br><sub>Cycle 自 5/19 起卡在 P0；被 H017 panel 與使用者 GO 所 block；鑑於 G5 F1 路線已大致被 LOSO 關閉，價值偏低。</sub> | 🔴blocked | 🟩低 | medium | medium | 2-4h |  | `vision-cycle3-resolve-mustfix`<br>`vision-low-f1-panel-survey-h017` |
|  | `vision-low-f1-panel-survey-h017` | H017 H_C3_2 panel 擴充：普查 TCGA/SEQC2 truth-pairs，尋找 >=4 個 low-F1（caller F1<0.80）且具備 truth set + V6 BAM 的樣本（5 個中只有 HCC1937 F1=0.37 符合；COLO829 truth-set 權限待核准）。<br><sub>整條 caller-F1-headroom 路線的關鍵阻塞點；真正的 panel 會是 C1 新數據，但這條路線本身在 LOSO 後屬於 low-ROI。</sub> | 🔴blocked | 🟩低 | high | hard | data-availability blocked | ◉ | — |
|  | `vision-register-fp-unexplained-axis` | 為 LOH/CN/cross_het 無法解釋的約 63% FP 註冊一個假說（Z-CHR8+Z-AUTO 僅捕捉約 37%）；需要新軸向（mappability/repeat/GC/SV）— 屬 Tier 4.2 reactive，尚未成為假說<br><sub>多數 FP 仍無法解釋；新軸向 pilot（mappability/repeat/SV）是唯一尚未耗盡的 FP 特徵刻畫方向（註：GC 已依護欄判為 NO-GO，故排除 GC）。</sub> | 💡idea | 🟧中 | low | hard | pilot 2-4h 後再註冊 | ◉ | — |
|  | `vision-register-zauto-kde-extension-t21` | 將 Z-AUTO KDE 跨 4 樣本擴展（T2.1，明訂為 V6 特徵刻畫 cycle 的 ⭐3->⭐4 升級需求）註冊進 hypothesis_queue.json<br><sub>這個已具名的 tier 升級需求目前未在佇列中追蹤；註冊後可讓 V6 特徵刻畫 cycle 得以完成。</sub> | 💡idea | 🟧中 | medium | medium | 註冊 + 跑 2-4h |  | `infra-archive-to-ism-rerun` |
|  | `vision-v6-production-tag-t12` | T1.2 V6 production tag 定案（Hard Gate，git tag v6-prod-{date}）：解鎖 Archive TO 7 樣本重跑 + T4.3 PI errata 套件。需要 COLO829 V6 ISM + 7 樣本 marker coverage + binary_commit_hash 寫入 manifest<br><sub>此 Hard Gate 把關 framework 論文 Tier 2 交付物與 Archive TO 7 樣本重跑。</sub> | ⚪todo | 🟥高 | medium | medium | 5 天工作流程（Hard Gate tag） | ★ | `phasing-v5-v6-production-decision`<br>`infra-colo829-v6-ism` |

### Phasing （12 任務）

| # | task_id | 任務 | 狀態 | 價值 | 信心 | 難度 | 工時 | 標記 | 依賴 |
|---|---------|------|------|------|------|------|------|------|------|
| 5 | `phasing-formal-wilcoxon-tpgap` | 對 6/6 TP-gap（中位數 +0.37）做正式的 Wilcoxon signed-rank，為 LOH-constrained phasing 發現提供可發表的顯著性檢定<br><sub>此發現目前僅靠描述性的 6/6 方向一致性；候選論文主軸需要正式檢定。</sub> | ✅done | 🟥高 | high | easy | 1 小時 | ★ | — |
| 6 | `phasing-negative-control-germline-hp-only` | 執行 --germline-hp-only=true negative control，確認它能消除 Inner-vs-Outer 的 NG=2 TP gap（如預測；Phase 1 實驗可重用）<br><sub>phasing 主軸宣稱需要其預測的 negative control，確認這個 gap 是 phasing 驅動而非 artifact。</sub> | ⚪todo | 🟥高 | high | medium | 2-4 小時 | ★◉ | — |
|  | `phasing-genome-wide-extension-4.19bias` | F-paired-D1：將 germline-absent 4.19:1 priority-bug bias（約 150K events）做全基因組延伸，確認超出 chr19（5,789 events）的跨染色體一致性<br><sub>把僅限 chr19 的 V5 Layer-1.5 caveat 量化推廣，這是 production 決策所需。</sub> | ⚪todo | 🟧中 | medium | hard | 半天 |  | — |
|  | `phasing-hcc1954-outlier-rootcause` | HCC1954 離群值根本原因：Potential_LOH 的可靠性（gap +0.35 但 Inner TP 絕對值最低 0.43；HER2+ 受 CNV 驅動的 germline-het AF 漂移）<br><sub>為論文的 caveats 章節解釋 6/6 發現中最弱的樣本。</sub> | ⚪todo | 🟩低 | medium | medium | 2-4h |  | — |
|  | `phasing-hpfinen-master-twoflag-phaseC` | Phase 2B/C：在 master dataset × 兩種 --germline-hp-only flag 上重新驗證 HPFineNGroups（chr19 + 4 樣本已完成；7 樣本 V5 BAM 全基因組 × 2-flag 待辦）— 解除 R-SELFREF<br><sub>決定 subclone marker（G2）的生物學立足點；需要 TO ISM 重跑後的 master。</sub> | ⚪todo | 🟧中 | medium | hard | multi-day | ★ | `infra-archive-to-ism-rerun` |
|  | `phasing-inner-samehap-flag-idea` | Filter 設計候選：將 Inner + same-hap NG=2 作為 TO somatic 高可信度 FLAG（非二元 FP filter；LOH 作為 filter 在 F1 上已 NEGATIVE）<br><sub>高可信度 flag（而非 filter）符合 characterization-only 的定位，並可避開 LOH-filter NEGATIVE 的護欄。</sub> | 💡idea | 🟧中 | medium | medium | 2-4h | ◉ | `phasing-negative-control-germline-hp-only`<br>`phasing-formal-wilcoxon-tpgap` |
|  | `phasing-ism-impact-v5-to-v3f` | F-paired-D3：量化把 V5 Layer 1.5 還原成 V3F（germline-absent -> hp=33）對下游 ISM 的影響；V6 patch 已部分處理<br><sub>為 production 升級決策量化 V5/V6 的 trade-off。</sub> | 🔵進行 | 🟧中 | medium | hard | 半天（需要 binary patch + ISM benchmark） |  | — |
|  | `phasing-paired-vs-methylation-decompose` | 在 PAIRED 模式下做獨立的 phasing-vs-methylation 驗證（等同 obs18 的 Inner/Outer NG=2 組成），以判定 paired AF×NGroups POSITIVE 是同一個 phasing 現象，還是真正的 methylation×AF coupling — 論文 KEY gate<br><sub>論文的核心宣稱（機制可單靠 phasing 重現）需要把 paired 層拆解；這需要現行 master 缺少的 paired_loh_bed_hit 欄位。</sub> | ⚪todo | 🟥高 | medium | hard | 多天 | ★◉ | `infra-unified-postkde-master` |
|  | `phasing-paired-vs-to-axis-align` | F-paired-D2：在 phase-block 內做 axis-aligned 分析，解決 paired-vs-TO 的 HP1/HP2 axis-swap caveat<br><sub>解決 phasing 敘事中的一個詮釋 caveat。</sub> | ⚪todo | 🟩低 | medium | medium | 2-4 小時 |  | — |
|  | `phasing-sp123-readlevel-audit` | SP1/2/3 read-level 稽核以釐清 H3（ClairS-TO PASS NonSomatic + GT=0\|0 + AF~0.5；相同投票卻給不同 tag）— 機制目前仍處於「待稽核」<br><sub>釐清 V5 HP2 與 V6 HP1 相同投票卻不同 tag 的機制；對主要目標價值不高。</sub> | ⚪todo | 🟩低 | low | hard | half-day |  | — |
|  | `phasing-top17-consensus-integration` | Top-17 / 16-cell 共識整合：Tier A「Extreme + T1 CN-loss + NG=2」= same-hap + LOH 雙重訊號，併入 phasing-signature 敘事。<br><sub>藉由結合雙重訊號強化敘事；重要性低。</sub> | ⚪todo | 🟩低 | low | medium | 2-4h |  | `phasing-formal-wilcoxon-tpgap` |
|  | `phasing-v5-v6-production-decision` | V5->V6 production 升級決策（V6 marker coverage 最佳 23,980 + germline-absent 保守採 hp=33，caller F1 相同；V5 在 SP1/2/3 phasing-weak loci 勝出）<br><sub>餵入 V6 production tag（T1.2 Hard Gate），而該 tag 把關 framework 論文與 Archive TO 重跑。</sub> | ⚪todo | 🟥高 | medium | medium | decision（已由使用者確認） | ★ | `phasing-ism-impact-v5-to-v3f`<br>`phasing-genome-wide-extension-4.19bias` |

### ASM （11 任務）

| # | task_id | 任務 | 狀態 | 價值 | 信心 | 難度 | 工時 | 標記 | 依賴 |
|---|---------|------|------|------|------|------|------|------|------|
| 2 | `asm-corrected-collapse-rerun` | 以校正後的 MAX-collapse 5mC+5hmC 進行完整 39,447-TP 重跑（取代 PROVISIONAL 折半的 beta；幅度約 2x 更大）— genome_survey_v2 批次 PID 411460 目前跑到 chr5<br><sub>已確認正在執行（PID 514459 chr5 TP）；解除對確定性 BRCA2 幅度的阻塞，並移除 WRONG 暫定的 CpG-set 歸因。</sub> | 🔵進行 | 🟥高 | high | medium | 執行中（數小時） | ★◉ | — |
| 10 | `asm-genome-survey-v2-complete` | 完成 genome_survey_v2 雙軸 TP-vs-FP ALLELE-axis 調查直到 chr22（FP 已完成 5,149；TP chr1-4 已完成，chr5-22 待處理）<br><sub>已確認執行中；這是即時測試任何 ASM 子集是否能區分 TP 與 FP（直接探測 G5 + 強化非區分性的先驗）。</sub> | 🔵進行 | 🟥高 | high | medium | 執行中（數小時） | ★◉ | `asm-corrected-collapse-rerun` |
| 12 | `asm-fix-level1-dup-bug` | 從源頭修正 Level1 重複輸出 bug（MSA C++ 將 p 與 1-p 輸出為每個 read-CpG 2 列 -> 2.6x 膨脹），或確認 MAX-collapse 完全抵銷，才能回報任何量化的 beta/n/p<br><sub>Agent A 高嚴重度 bug 使 n 膨脹約 2.4x，導致絕對 beta/p 不可信；卡住任何超過 3 的 tier 升級。C++ 變更需要 methodology-audit + 編譯 + 使用者確認。</sub> | ⚪todo | 🟥高 | high | medium | 2-4h（C++ 變更：methodology-audit + 編譯 gate + 使用者確認） | ★ | `asm-genome-survey-v2-complete` |
|  | `asm-audit-prior-msa-beta-numbers` | 稽核 Level1 dup-emission 是否影響專案內所有先前以 MSA 為基礎的 beta 數字（跨 thread 稽核需求）<br><sub>若此 bug 出在共用的 MSA export，其他 thread 的 methylation 數字可能需要相同修正。</sub> | ⚪todo | 🟧中 | medium | medium | 半天 |  | `asm-fix-level1-dup-bug` |
|  | `asm-confound-guard-tpfp` | 若要用 ASM 來區分 TP/FP，必須通過 /auc-confound-guard（within-group OLS + AF-bin + permutation）；既有結論指出 ASM 不具區分力 AUC<0.60，而 v2 的 TP-vs-FP run 正是這個即時測試<br><sub>任何 TP/FP 區分力宣稱前的強制關卡；先前證據使 positive 結果機率偏低（C1 = 更細尺度的 v2 資料足以正當化重新測試）。</sub> | ⚪todo | 🟩低 | low | hard | 半天 |  | `asm-genome-survey-v2-complete` |
|  | `asm-cross-sample-validation` | 跨樣本 ASM 驗證，涵蓋 7 個樣本中的 >=4 個（目前僅單樣本 HCC1395 paired_full），以將 tier 提升至 3 以上<br><sub>在 cross-sample 方向/效應一致性獲得確認前，tier 上限為 3；論文中任何可泛化的 ASM 主張都需要此驗證。</sub> | ⚪todo | 🟥高 | medium | hard | 多天 | ★ | `asm-fix-level1-dup-bug` |
|  | `asm-finalize-experiment-doc` | 待 corrected-collapse 量值定案後，將 ASM 實驗文件 + evidence_ledger 更新為最終版（去除 PROVISIONAL 標記）；透過 conclude-research 把 in_progress -> validated<br><sub>PROVISIONAL 數字必須替換、ledger 須定案，才能用於論文。</sub> | ⚪todo | 🟧中 | high | easy | 1-2 小時 | ★ | `asm-fix-level1-dup-bug`<br>`asm-rebuild-negative-control` |
|  | `asm-powerup-tsg-enrichment` | 強化 TSG 啟動子 enrichment 檢定（目前只有 3/14,222 個 variant 落在 TSG 視窗 -> power 不足），在宣稱「no enrichment」前先擴大視窗或合併樣本以提升 power<br><sub>目前 enrichment=0.00x 是 power 不足而非證據；擴大視窗可避免過度宣稱 negative 結果。</sub> | 💡idea | 🟩低 | medium | medium | 2-4 小時 | ◉ | `asm-cross-sample-validation` |
|  | `asm-rank-strong-sites` | 排序/特徵化 19 個 Bonferroni+\|Delta-beta\|>=0.1 強訊號位點 + ALLELE-axis 高 \|Delta-beta\| 位點（chr1:200029019 -0.350，ALT_vs_REF 達 -0.61），以找出超越 BRCA2 那個 modest -0.063 的效應<br><sub>效應較大的位點比 BRCA2 那個 modest 效應更適合作為論文 showcase。</sub> | 🔵進行 | 🟧中 | medium | medium | 2-4h | ◉ | `asm-genome-survey-v2-complete` |
|  | `asm-rebuild-negative-control` | 用校正後的 pipeline 重建 ASM negative control，採用 n>=50-100 個配對對照（TVAF + 局部 CpG 密度 + coverage 配對）（目前 n=5，僅 3 個有效）<br><sub>目前的 empirical null 在統計上太弱，無法支撐 BRCA2 p-value 主張。</sub> | ⚪todo | 🟧中 | high | medium | 2-4h | ★ | `asm-fix-level1-dup-bug` |
|  | `asm-validate-allele-axis-confounds` | 驗證 ALLELE-axis 可用性：在 non-LOH 區證明 HP-axis 與 ALLELE-axis 的 sign 一致性，並加入 MAPQ/read-length 過濾以控制 ALT-read mapping bias / collider bias<br><sub>ALLELE-axis 的大 effect（-0.61）唯有在控制 collider/mapping-bias 後才可信。</sub> | ⚪todo | 🟧中 | medium | hard | 半天 |  | `asm-genome-survey-v2-complete` |

### Phase2 FP-filter （5 任務）

| # | task_id | 任務 | 狀態 | 價值 | 信心 | 難度 | 工時 | 標記 | 依賴 |
|---|---------|------|------|------|------|------|------|------|------|
| 7 | `methylfilter-conclude-negative-archive` | 正式 /conclude-research NEGATIVE 撰寫 + 封存 methyl_augmented_filter_phase2（已驗證 manifest 仍為 status:initiated、conclusion:null；pre-reg 指出 all-H-NEGATIVE -> 結論 NEGATIVE）。Hard Gate：NO-GO 需使用者確認<br><sub>已驗證 manifest 尚未收尾；將其關閉可把 LOSO 發現轉為可引用的嚴謹 NEGATIVE，且屬於需使用者確認的 Hard Gate。</sub> | ⚪todo | 🟥高 | medium | easy | 2-4h | ★ | — |
| 13 | `methylfilter-cycle5-pivot-decision` | Cycle 5 轉向決策，於 Path A（phase_block_3d/thread_d，建議）、Path B（chr8 zone gate）、Path C（low-F1 panel）之間抉擇。Hard Gate：宣告 LR filter NO-GO 需使用者確認；已驗證無 cycle5/ 目錄<br><sub>已驗證沒有 cycle5 目錄；整條 filter 路線都停在這個使用者 pivot 決策上；Path A 與論文主軸一致。</sub> | 🔴blocked | 🟥高 | medium | medium | decision（使用者觸發） | ★ | `methylfilter-conclude-negative-archive` |
|  | `methylfilter-hnew4-low-f1-panel-followup` | 調查 HCC1395 H_NEW_4 +0.00699（drop-caller_af, tau=0.95）：這是真實的跨樣本效果，還是單樣本 artifact？需要 low-F1 panel 普查（Path C）。標記為 post-hoc／非 confirmatory。<br><sub>Post-hoc 且 SANITY-VIOLATED 的發現；只有真正的 low-F1 panel（C1 新數據）才能區分是 artifact 還是 signal — 受 caller_af 限制，先驗機率低。</sub> | ⚪todo | 🟩低 | low | hard | multi-day | ◉ | `methylfilter-cycle5-pivot-decision`<br>`vision-low-f1-panel-survey-h017` |
|  | `methylfilter-track-b-v3f-v5-rerun` | 延後的 Track B 跨樣本：為 H_PHASE2_2 重跑 4 樣本 V3F+V5 ISM（約 10hr 平行）— BLOCKED，V3F/V5 的 4 樣本 BAM 實際上不存在（只有 V6 phaseD_v6_5sample）。<br><sub>因 BAM 不存在而被實體 block；只有在 filter 方向復活時才有意義，但 LOSO 已實質上把它關閉。</sub> | 🔴blocked | 🟩低 | high | hard | multi-day（資料不可用） |  | — |
|  | `methylfilter-update-sot-after-loso` | 在 LOSO reframe 後完成 SoT 更新（PI Trust HTML + Cycle1 報告 banner + 論文 Section 3；ledger/memory/CURRENT_FOCUS 已確認完成）。<br><sub>其餘 SoT 介面（PI Trust HTML、報告 banner、論文 Section 3）必須反映 +0.02236 的降級，以避免過時的標題式 claim。</sub> | 🔵進行 | 🟧中 | medium | easy | 1-2h | ★ | — |

### HKU handoff （13 任務）

| # | task_id | 任務 | 狀態 | 價值 | 信心 | 難度 | 工時 | 標記 | 依賴 |
|---|---------|------|------|------|------|------|------|------|------|
| 3 | `hku-commit-v4-deliverables` | git commit v4 HKU handoff（research/hku_collaboration/ + docs/concepts/2026/05 + docs/references/2026/05 + docs/reports/in_progress）— 非 C++，非 Hard Gate<br><sub>整個 D-handoff 都未納入版控；提交可保護約 24 張圖 + 12 個腳本 + 3 份文件不致遺失，並消除 orphan-audit 風險。</sub> | ⚪todo | 🟥高 | high | easy | 30 分鐘 | ★ | — |
| 4 | `hku-html-email-safe-base64` | 讓 HKU standalone HTML 適合 email 寄送 — 將 18 個外部相對路徑的 <img> 轉為內嵌 base64（已驗證 0 個 base64 / 18 個外部 img）<br><sub>已驗證該『standalone』HTML 有 0 張內嵌圖片與 18 個移動後即失效的外部參照，因此無法照現狀 email 給 HKU。</sub> | ✅done | 🟥高 | high | easy | 1-2 小時 | ★ | — |
| 8 | `hku-verify-f4-hp-counts` | 對外交付前，將 F4 HCC1395 baseline HP counts（863K/266K/79K/23K/630）與 evidence_ledger 交叉核對（README 標註為憑記憶引用、未驗證）<br><sub>依『既有 artifact 必先驗證』規則，外部交付物中憑記憶引用的數字必須經 ledger 確認。</sub> | ✅done | 🟧中 | medium | easy | 30 分鐘 | ★ | — |
| 9 | `hku-replace-f5-synthetic-heatmap` | 對外交付前，將 F5 合成 rng.beta ISM heatmap placeholder 替換為真實 HCC1395 ISM region（候選 chr2:18,086,020 +-32kb）<br><sub>將合成的 placeholder 資料交付給外部合作者是科學誠信風險；必須使用真實 ISM 輸出。</sub> | ✅done | 🟧中 | medium | medium | 2-4 小時 | ★ | — |
|  | `hku-align-open-questions` | 與 HKU 對齊 4 個開放問題（linear-depth 命名採用、ISM cluster cutoff 度量 + tau、X-ratio 上限、HP1-vs-HP2 不平衡作為 LOH proxy）<br><sub>這些屬雙邊討論點，僅能在 HKU 收到並回覆報告後才能解決。</sub> | 🔴blocked | 🟧中 | high | medium | 交付後討論 |  | `hku-pi-signoff-email` |
|  | `hku-cross-sample-extension` | 將 HKU 發現跨樣本延伸至 COLO829/H1437/H2009（設計上延後至 HKU 提出需求為止；目前僅 HCC1395）<br><sub>依 plan v4 延後；chr8 99% LOH 為 HCC1395 特有，僅在 HKU 提出要求時才動作（需 BAM 可用）。</sub> | 🔴blocked | 🟧中 | high | hard | 多日 | ◉ | `hku-pi-signoff-email` |
|  | `hku-index-registration` | 為 HKU handoff 新增 experiments INDEX.md（或適當的 SoT）條目，使其不會在 data-audit/provenance-tier-audit 中被標為 orphan<br><sub>grep 確認 HKU handoff 沒有任何 INDEX 條目；SoT 登錄可補上此 provenance 缺口。</sub> | ⚪todo | 🟧中 | medium | easy | 30 分鐘 | ★ | `hku-commit-v4-deliverables` |
|  | `hku-memory-landing` | 記憶落地：建立 project_hku_collaboration_handoff.md + 於 MEMORY.md 新增索引行（plan 列為完成後義務）<br><sub>記憶條目為必要的完成後義務，可避免未來 orphan/recall 失敗。</sub> | ⚪todo | 🟧中 | high | easy | 30 分鐘 | ★ | `hku-commit-v4-deliverables` |
|  | `hku-pi-signoff-email` | 對外寄送 HKU 報告前先取得 PI 簽核（依文檔治理規範經由 pi_reports/ 流程；不可將 in_progress 草稿交付外部單位）<br><sub>對外交付需 PI 核准；in_progress 草稿不可擅自送交外部協作方。</sub> | ⚪todo | 🟥高 | medium | easy | 30min + PI 回覆時間 | ★ | `hku-translate-html-en`<br>`hku-replace-f5-synthetic-heatmap`<br>`hku-verify-f4-hp-counts` |
|  | `hku-reconcile-naming-convention` | 釐清命名的待決事項（sibling-vs-child）；確認沒有殘留的 sibling 用語（交付物已採用 linear-depth child 慣例）<br><sub>plan v4 Open Decisions #3 仍標註存在模糊；對交付物 grep 殘留的 sibling 標記。</sub> | ⚪todo | 🟩低 | low | easy | 30 分鐘 |  | — |
|  | `hku-relinearize-axes-verify` | 驗證 A2_1/A2_2 軸確實為線性（review #2 log->linear；圖於 5/23 19:02 重新產生但尚未獨立確認）<br><sub>開啟 PNG 確認某個 user review 點確實已關閉。</sub> | ⚪todo | 🟩低 | medium | easy | 15 分鐘 |  | — |
|  | `hku-translate-bc-qcards-en` | 將 HKU 交付物 B + C 翻成英文（或至少 Q1-Q4 這 4 張討論點卡片）<br><sub>Q1-Q4 卡片必須讓 HKU 能回答；完整 B/C 翻譯則屬於 nice-to-have。</sub> | ⚪todo | 🟧中 | low | medium | 半天 |  | `hku-translate-html-en` |
|  | `hku-translate-html-en` | 將 HKU 交付物 A（HTML）翻成英文，供說英語的 HKU somatic team 閱讀（目前為 zh-TW）<br><sub>Task type D handoff 要求交付物須能被接收方閱讀；HKU 為說英語的外部單位。</sub> | ⚪todo | 🟥高 | medium | medium | 半天 | ★ | `hku-html-email-safe-base64` |

### Infra/Data （17 任務）

| # | task_id | 任務 | 狀態 | 價值 | 信心 | 難度 | 工時 | 標記 | 依賴 |
|---|---------|------|------|------|------|------|------|------|------|
| 1 | `infra-maintain-guardrails` | 維護操作護欄：>100MB 的 pipeline 須設 TMPDIR=/big7_disk + 明確指定 -o + df 預檢；讀取 merged file 時須 archive caller_af + 過濾 phase1_new；新伺服器須 mount\|grep tmp<br><sub>永久性基礎設施護欄；/tmp 磁碟寫滿事故曾讓全部 33 個帳號當機，故每條 long pipeline 都必須遵守（由 pipeline_block_check.sh Hard Gate 把關）。</sub> | 🔵進行 | 🟥高 | high | easy | 持續進行 | ★ | — |
| 15 | `infra-archive-to-ism-rerun` | 候選 #4 / Archive TO 6 樣本 ISM 重跑，套用新 KDE + LOH.bed + germline-hp-only=off，將 LOH_Subtype/CovM_used/baseline_x 加入 master_extended（Path 2 約 4.5h 或全量約 10h）。USER 觸發的長計算；export TMPDIR + -o<br><sub>已驗證 master_extended 為 pre-KDE（缺 LOH/CN）；此重跑是跨樣本 TO 驗證、HPFineN two-flag 與統一 master 的前置條件。必須遵守 TMPDIR guardrail。</sub> | ⚪todo | 🟥高 | high | medium | 4.5-10h 平行 | ★ | `vision-v6-production-tag-t12` |
|  | `infra-chr8-hotspot-deep-analysis` | 深入分析 HCC1395 chr8 LOH+ASM hotspot（7.4x FP enrichment，sample-specific）：拉取 chr8 LOH bed + methylation heatmap 以刻畫 ASM 機制；候選的 sample-aware feature。<br><sub>刻畫 chr8 hotspot 機制（H4 POSITIVE）並餵入 methylfilter Path B 的 chr8 zone gate；屬 sample-specific，故非全域 marker。</sub> | ⚪todo | 🟧中 | medium | medium | 2-4h | ◉ | — |
|  | `infra-colo829-bam-decision-a-vs-c` | COLO829 BAM 決策：選項 A（沿用現有）對比 C（約 12h tumor BAM 重建 + KDE 校正 QS 重跑）。USER 決策，gating P1/P2/P3<br><sub>未決的 user A-vs-C 決策阻擋 COLO829 QS 重跑及其下游（P1/P2/P3）。</sub> | 🔴blocked | 🟧中 | medium | easy | user 決策 | ★ | — |
|  | `infra-colo829-melanoma-ism-p2` | P2：在 KDE 後 QS 上重新測試「ISM 對黑色素瘤（COLO829）無效」假設<br><sub>相依的 characterization 觀察；鑑於 COLO829 無 methylation basecall，價值偏低。</sub> | 🔴blocked | 🟩低 | low | medium | 2-4h | ◉ | `infra-colo829-qs-rerun-p1` |
|  | `infra-colo829-qs-rerun-p1` | P1：以 KDE 校正的 CNV_Loss penalty 重跑 COLO829 QS（QS 可能回升，約 12h）<br><sub>QS 僅供 characterization（依 guardrail，TO QS AUC=0.497 為隨機）；重跑只解決 stale-binary caveat，不影響 F1。</sub> | 🔴blocked | 🟩低 | medium | hard | 約 12h |  | `infra-colo829-bam-decision-a-vs-c` |
|  | `infra-colo829-v6-ism` | COLO829 V6 ISM 完成（Archive TO 重跑 + KDE 校正）— 7 樣本 marker coverage 所需，供 V6 production tag 使用<br><sub>COLO829 step05_intersubmod 為空；其 V6 ISM 為 7 樣本 coverage 與 production tag 所必需（注意：COLO829 ONT_R10 無 methylation，故 methylation 相關 claim 排除它）。</sub> | 🔴blocked | 🟥高 | medium | hard | 約 12h | ★ | `infra-colo829-bam-decision-a-vs-c` |
|  | `infra-copy-hcc1395-normal-bam` | 將 HCC1395 normal tagged BAM 從 /big8_disk 複製到 /big7_disk（優先使用 samtools view -L region-subset 約 2-3k 個 NonLOH 區域，以避免完整 136GB 傳輸；保護接近滿載的 /big8）。<br><sub>為 Normal-BAM Sample-ASM pilot（G4）的前置條件；region-subset 可避免塞爆接近滿載的 /big8。</sub> | ⚪todo | 🟧中 | high | easy | transfer-bound | ★ | — |
|  | `infra-cpp-normalbaseline-bugfix` | 候選 #5：透過 /cpp-change 修正 NormalBaseline.cpp writer bug（R-DATA-GAP）（待使用者核准，2-4h）。C++ Hard Gate：編譯 + commit。<br><sub>修正一個 data-gap writer bug；需使用者核准 + C++ gate；唯有在它阻擋下游 consumer 時才有必要。</sub> | ⚪todo | 🟩低 | medium | medium | 2-4h（C++ Hard Gate） |  | — |
|  | `infra-cpp-normalreader-phase1b` | C++ Phase 1B：擴展 RegionProcessor.cpp 的 process_single_region() 以接受 normal_reader（目前僅 tumor_reader）。C++ Hard Gate：methodology-audit + 編譯 + 使用者確認。<br><sub>讓 pipeline 內可使用 normal methylation reference（G4）；唯有 pilot 為正向時才值得跑 C++ Hard Gate。</sub> | ⚪todo | 🟧中 | medium | hard | half-day（C++ Hard Gate） |  | `infra-normal-asm-pilot` |
|  | `infra-cross-sample-qs-p3` | P3：KDE 後重跑的跨樣本 QS 觀察<br><sub>QS 僅供 characterization；不在任何 F1 或主軸關鍵路徑上。</sub> | 🔴blocked | 🟩低 | low | medium | 2-4h | ◉ | `infra-colo829-qs-rerun-p1` |
|  | `infra-normal-asm-pilot` | 執行 HCC1395 normal-BAM Sample-ASM vs Normal-ASM pilot；僅在 Delta_ASM AUC>=0.60 且方向明確時，才擴展至其他 5 個 normal BAM（先以 ls 驗證路徑）。<br><sub>對持續性 FP（G4）以 normal-reference ASM 進行 pilot-gated（AUC>=0.60）測試；有防護的擴展可保護磁碟。</sub> | ⚪todo | 🟧中 | medium | medium | 2-4h | ◉ | `infra-copy-hcc1395-normal-bam` |
|  | `infra-per-sample-heterogeneity-diag` | 針對 H2009/HCC1937/H1437 的 per-sample 異質性診斷（paired multi-bio 增益不穩定；H2009 有 68% FP 落在 LOH 區）<br><sub>解釋跨樣本不穩定性；僅作為輔助性描述特徵刻畫。</sub> | ⚪todo | 🟩低 | low | medium | 2-4h | ◉ | — |
|  | `infra-platform-methyl-normalization` | 平台/分組感知的 methylation 特徵正規化（5kHz/DORADO/PAO/Google 在 read 品質與 methylation baseline 上各有差異）— 目前尚未註冊<br><sub>任何跨 basecaller 泛化宣稱都需要跨平台正規化（已證實 AF filters 無法跨 basecaller 移植）。</sub> | 💡idea | 🟧中 | low | hard | 多天 | ◉ | `infra-unified-postkde-master` |
|  | `infra-requantify-pre0419-covm` | 用新的 per-sample KDE baseline 重新量化 H-CN1/H-CN2 及所有 2026-04-19 前的跨樣本 CovM 分析；為 stale-master 結論標註 75.0 caveat<br><sub>從先前的 CN 結論中移除過時的 75.0 hardcoded-coverage artifact。</sub> | 🔵進行 | 🟧中 | medium | medium | 2-4h |  | — |
|  | `infra-to-pure-diagnostics-phase1a` | Phase 1A TO-pure 診斷：追查 methyl-delta F1 regression（-0.0206）的根本原因，集中於 TP-REF 誤分類 + FP Subclone/Strong；可能拆分 TO/paired 建模。<br><sub>對已知 TO regression 的診斷；鑑於 TO ISM 依 guardrail 近乎無用，價值偏低。</sub> | ⚪todo | 🟩低 | low | medium | 2-4h | ◉ | — |
|  | `infra-unified-postkde-master` | 產出新的統一 post-KDE master（7 樣本 × 2 模式，一致的 caller_af + 完整 LOH + KDE baseline），以汰除過時的 2026-03-30 master 以及 merged_7samples 陷阱檔案。<br><sub>單一可信賴的 master 是 paired phasing-vs-methylation 分解與所有跨樣本工作的資料基礎；汰除兩個已知的資料陷阱。</sub> | ⚪todo | 🟥高 | medium | medium | 2-4h | ★ | `infra-archive-to-ism-rerun` |

### Guardrails （3 任務）

| # | task_id | 任務 | 狀態 | 價值 | 信心 | 難度 | 工時 | 標記 | 依賴 |
|---|---------|------|------|------|------|------|------|------|------|
| 11 | `guardrails-reopen-gate-enforce` | 對全部 25+ 項已收尾的 NEGATIVE/NO-GO 結果強制執行 REOPEN GATE：任何 reopen 都必須引用具體的 C1（新數據）/ C2（新方法）/ C3（新前置條件），並在開始計算前記錄於 research/<topic>/00_INDEX.md 的 Pre-registration<br><sub>範圍保護治理；防止浪費 cycle 並避免五大目標中各處的 HARKing（依 feedback_productive_failure_reopen_threshold）。</sub> | 🔵進行 | 🟥高 | high | easy | 持續進行 | ★ | — |
|  | `guardrails-memory-index-reconcile` | Memory-consolidation 巡檢：確認每一列 INDEX.md NO-GO 都有對應的 MEMORY.md Concluded 條目（P3 window、P4 second-hit、Zone-Aware F1、SEQC2 CNV、GC/CN、L1-L4 文獻有重疊但非 1:1）<br><sub>弭平 INDEX 與 MEMORY 的 drift，讓未來的 agent 能可靠地看到每一道護欄。</sub> | ⚪todo | 🟩低 | high | easy | 1-2h |  | — |
|  | `guardrails-reverify-cpp-fileline` | 在任何新決策中將其當作 live fact 宣稱前，先對照目前原始碼重新驗證已收尾機制內承載關鍵的 C++ file:line 引用（例如 LabelTest.cpp:265-302）<br><sub>Provenance 新鮮度；file:line drift 可能讓某個用於新決策的機制 claim 失效。</sub> | ⚪todo | 🟩低 | medium | easy | 每次使用 30min |  | — |

## 6. 跨主軸 Blockers

| 嚴重度 | Blocker | 阻擋的任務 |
|--------|---------|-----------|
| HIGH | V6 production tag（T1.2 Hard Gate）尚未定案 — 需要把 COLO829 V6 ISM + 7-sample marker coverage + binary_commit_hash 寫入 manifest，再由 user 核准 git tag。 | `vision-v6-production-tag-t12` · `infra-archive-to-ism-rerun` · `infra-unified-postkde-master` · `phasing-paired-vs-methylation-decompose` · `phasing-hpfinen-master-twoflag-phaseC` · `vision-candidate-hpfinen-reverify-phase2b` |
| HIGH | COLO829 BAM A-vs-C 的 user 決策 + 空的 step05_intersubmod — 把關著餵入 V6 tag 的 7-sample coverage 以及所有 QS 重跑。 | `infra-colo829-bam-decision-a-vs-c` · `infra-colo829-v6-ism` · `infra-colo829-qs-rerun-p1` · `infra-colo829-melanoma-ism-p2` · `infra-cross-sample-qs-p3` |
| HIGH | Cycle 5 pivot 方向（A/B/C）需要 USER 決策；宣告 LR filter 為 NO-GO 屬 Hard Gate。目前不存在 cycle5/ 目錄（已驗證）。 | `methylfilter-cycle5-pivot-decision` · `methylfilter-hnew4-low-f1-panel-followup` · `phaseblock3d-init-research` · `phaseblock3d-h013-cis-effect` · `phaseblock3d-h014-boundary-artifact` · `phaseblock3d-h015-zone-joint-score` |
| HIGH | Low-F1（caller F1<0.80）跨樣本 panel 並不存在 — 5 個樣本中僅 HCC1937 符合；n>=4 沒有 BAM（COLO829 truth-set 權限待批）。對整條 caller-F1-headroom filter 路線（cycle 3）是硬性 blocker。 | `vision-low-f1-panel-survey-h017` · `vision-cycle3-write-planjson` · `vision-cycle3-h018-negative-control` · `methylfilter-hnew4-low-f1-panel-followup` |
| HIGH | ASM PROVISIONAL 量級（5mC+5hmC 合併使 beta 減半）+ Level1 重複輸出 C++ bug（2.6x 行數膨脹），使絕對 beta/n/p 不可信；只有 Delta-beta 的正負號是安全的。Tier 上限為 3。 | `asm-fix-level1-dup-bug` · `asm-rebuild-negative-control` · `asm-cross-sample-validation` · `asm-finalize-experiment-doc` · `asm-confound-guard-tpfp` |
| HIGH | 外部 HKU 交付需要 PI sign-off + 英文翻譯 + email-safe HTML + 真實（非合成）圖；目前 artifacts 為 zh-TW、in_progress、未 commit，外部 img 引用壞掉（已驗證），且含一個 F5 合成 placeholder。 | `hku-pi-signoff-email` · `hku-cross-sample-extension` · `hku-align-open-questions` |
| MEDIUM | Paired phasing-vs-methylation 拆解（論文 KEY gate）需要當前 master 缺少的 paired_loh_bed_hit 欄位；需重跑 cross_sample_audit / unified master。 | `phasing-paired-vs-methylation-decompose` |
| MEDIUM | master_extended 為 pre-KDE（已驗證缺 LOH_Subtype/CovM_used/baseline_x）且 COLO829 step05 為空 — 僅 AF×NG×NumReads 可觀察；把關著跨樣本 TO 驗證與 HPFineN two-flag。 | `infra-archive-to-ism-rerun` · `phasing-hpfinen-master-twoflag-phaseC` · `vision-candidate-hpfinen-reverify-phase2b` · `infra-unified-postkde-master` |
| HIGH | 25+ 個已收尾 NEGATIVE/NO-GO guardrails 的 REGRESSION 風險 — 任何缺乏指名 C1/C2/C3 的重啟都會浪費 cycle 並有 HARKing 風險。 | `guardrails-reopen-gate-enforce` · `vision-register-fp-unexplained-axis` · `asm-confound-guard-tpfp` · `methylfilter-hnew4-low-f1-panel-followup` |
| MEDIUM | 所有 .cpp 變更都受 C++ Hard Gate（methodology-audit + compile + user 確認）— 適用於 ASM Level1 dedup、normal_reader 延伸、NormalBaseline writer 修正。 | `asm-fix-level1-dup-bug` · `infra-cpp-normalreader-phase1b` · `infra-cpp-normalbaseline-bugfix` |
| HIGH | /tmp 磁碟寫滿的 infra 危害 — 任何把 GB 級中間檔預設寫到 root /tmp 的 pipeline 都可能拖垮共用的 33 帳號伺服器（由 pipeline_block_check.sh Hard Gate 把關）。 | `infra-archive-to-ism-rerun` · `infra-colo829-v6-ism` · `asm-corrected-collapse-rerun` · `asm-genome-survey-v2-complete` |

## 7. 新任務點子（合成提出，待用戶決定是否納入）

- **[🟥高] 將三條可挽救的主線（LOH-constrained phasing 主軸 + ASM characterization + LOSO sample-level-NEGATIVE 方法學發現）整合成單一篇「read-level epigenetic & haplotype context for somatic variant interpretation」論文，把 NEGATIVE 當作嚴謹的方法學章節，而非被丟棄的結果。**
  - 為何：三條主線各自獨立收斂到同一定位（characterization，而非 filter）；統一敘事可把已降權的 F1 失敗轉為方法學貢獻，並讓 phasing 發現有所歸屬。已記錄為 task paper-unify-three-axes-outline。
- **[🟥高] 建立全專案的可重現性 manifest（透過 /pipeline-manifest），橫跨 ASM、phasing、methyl-filter、HKU，對應每個 script -> inputs -> outputs -> figures -> reports，在任何投稿前偵測 orphan（例如 F5 合成 placeholder、未 commit 的 HKU 腳本）。**
  - 為何：多條主線都有 provenance 缺口（HKU 未 commit/未編入索引、F5 合成資料、過時的 2026-03-30 master 與 merged-trap 並存）；單一 manifest 可把這些全部攤出，且為 supplementary methods 附錄所必需。
- **[🟥高] 設計一套明確的跨樣本泛化協議（LOSO-as-circularity-test + 各 feature 方向一致性 + caller-F1-headroom 稽核）作為可重用的方法學資產，並一致地套用到 ASM 跨樣本驗證與 phasing 全基因組延伸。**
  - 為何：Phase-2 主線已產出此 playbook（VIF 稽核 + L2 ridge + MNAR NaN + multi-seed + LOSO）；將其形式化可避免目前仍是單樣本/僅 chr19 的 ASM 與 phasing 主線重蹈同樣的 single-sample 循環性陷阱。
- **[🟩低] 純粹為了裁決 HCC1395 H_NEW_4 +0.00699（drop-caller_af）這個 post-hoc 發現而執行一個受控的 low-F1 panel survey：它是真實的跨樣本效應，還是單樣本 artifact？預先登記約為 0 的 prior，使其維持 confirmatory 而非 HARKing。**
  - 為何：這是原本已收尾的 G5 路線中（CURRENT_FOCUS 5/20）唯一真正未決的問題；一個預先登記的 panel test（C1 新數據）可乾淨地解決它，而不重啟 filter 敘事。已記錄為 methylfilter-hnew4-low-f1-panel-followup。
- **[🟧中] 為 LOH/CN/cross_het 無法解釋的那約 63% FP，登記一個「new FP axis」pilot（mappability / repeat / SV — 明確排除已 NO-GO 的 GC），作為投入假設前的 observation-first pilot。**
  - 為何：大多數 FP 仍無法解釋，而既有的每條軸都已耗盡（guardrails）；一條全新的軸是唯一尚未關閉的 FP-characterization 方向，且可能實質推進論文。已記錄為 vision-register-fp-unexplained-axis。
- **[🟧中] 跨主線 MSA beta 稽核：判定 Level1 重複輸出（duplicate-emission）bug 是否汙染了 ASM 主線以外的任何先前專案 methylation 數字（例如 methyl-filter pilot、snv_methylation_association）。**
  - 為何：若此 bug 存在於共用的 MSA C++ export，多項已收尾的 methylation 結果可能帶有膨脹的 n / 減半的 beta；一次性稽核可避免把錯誤數字傳播進論文。已記錄為 asm-audit-prior-msa-beta-numbers。
- **[🟧中] 做一次 CURRENT_FOCUS 的 live 刷新，移除如今已過期的 5/24 HKU deadline 框架，並把文件重新錨定在實際的當前關鍵路徑上（論文主軸 + ASM survey 完成），因為 SessionStart hook 每個 session 都會注入它。**
  - 為何：CURRENT_FOCUS 會被自動注入，目前卻被一個過期 deadline 主導；context 頂端文件一旦過時，就會誤導每個未來 session。低工作量、高槓桿。

## 8. 護欄 — 已 concluded NEGATIVE / NO-GO（勿重啟，除非具名 C1/C2/C3）

> C1=新數據 / C2=新方法 / C3=新前置條件。任何想重啟以下方向的任務，rationale 必須具名其中一項。

- REOPEN GATE（governance）：任何已結案的 NEGATIVE/NO-GO 需滿足以下至少一項：C1 new data（new sample/platform/truth set/scale）／C2 new method（繞過已知限制，例如 pooled->within-group OLS）／C3 new precondition（upstream dependency 已修復，例如 KDE bug fix、V5->V6 binary、HP tag fix）。Spaced recall（30d/90d）並不授權 reopen。Reopen 必須記錄於 research/<topic>/00_INDEX.md 的 Pre-registration 中。
- ReadParser --germline-hp-only Phase 1 CONDITIONAL NEGATIVE (2026-04-21)：機制正確，但 HCC1395 TO 全基因組（40,115 sites）顯示沒有任何 TSV 特徵的 AUC 改善 >=0.02；flag=on 使 NG>=3 區域 collapse 為 0；不進入 Phase 2。HPFineNGroups marker 結論需要 flag=on 重新驗證。
- R1-Global Sample ASM NEGATIVE (2026-04-21)：HCC1395 全部 40,237 variants；SampleASM_Delta residualized AUC 0.527 [0.520, 0.533]（隨機）；CL-025a 降級 star3->star2；CL-008 <=0.58 上限提升至 star5。
- ClairS-TO Verdict Pilot 在 F1 上 NEGATIVE：內部校準 POSITIVE（Germline FP 96.1%、Somatic TP 91.8%），但 Verdict_Germline 100% 落在 LowQual；S1 deltaF1=0；主要升級路徑改導向 Wakhan/SAVANA。
- Z3 Internal Exploration NEGATIVE：Z3 內部 12 個特徵全部 AUC<0.61；AF germline 分層僅在 1/7 樣本有效；HCC1954 amplicon 為已知。
- O11 Within-Group Methylation Heterogeneity NEGATIVE (2026-03-31)：epipolymorphism AUC=0.845 完全是 n_reads confound（residualized 後 -> 0.530）。
- O12 TO LOH Methylation Scenario Discrimination NEGATIVE (2026-03-31)：AlleleDelta 是 AF confound；CramersV L2 AUC=0.80 為 collider bias。
- O13 Cross-Region Methylation Correlation NEGATIVE (2026-03-31~04-01)：FP>TP correlation 是 shared-read-count confound；在 stratify+residualize+matching 後消失。
- G1-G7 TO Germline FP identification NO-GO (2026-04-01)：48 張圖 x 60+ 特徵全部 AUC<0.64；在 TP loss <=2% 的條件下，FP removal = 0%。
- Read-Level Haplotype & Methylation germline-FP CONDITIONAL NO-GO (2026-04-02)：LOSO AUC=0.721（首次 >0.70），但在安全約束下 FP removal=0%；根本原因為 high-purity 細胞株。low-purity 樣本仍保有潛力（唯一的 C1 式 reopen 槓桿）。
- Option C LabelTest 雙路徑（HP-free vs HP-dependent）NEGATIVE（2026-04-07）：cluster_labels 已經是 HP-free；ClusterPermanovaF AUC=0.512（隨機）；HP-free 5 特徵組合 AUC=0.564 vs HP-dependent 0.598；純 methylation 分群無區辨力；C++ 變更取消。
- O9 FN 特徵觀察 NO-GO（2026-04-08）：7 樣本 x 2 模式，122,790 個 FN 區域；HP-free methylation 特徵全部 AUC<0.53；最強訊號 LabelAllelePermanovaF=0.664 是 AF proxy 而非 methylation；TO QS AUC=0.338（FN>TP 反轉）。ISM 無法挽救 FN；在 methylation 空間中 FN==TP。
- TO-pure 獨立建模 NEGATIVE（2026-04-08）：7 樣本 LOSO；HP-free AUC=0.53、All-ISM 0.60-0.64、Caller-only 0.63、ISM+Caller 0.66；ISM 僅增加 +0.003~+0.030。單靠 caller_af（0.654）即超越所有 ISM 特徵；ISM 在 TO 模式下近乎無用。
- Fine-Pairwise Distance 分析 NEGATIVE（2026-04-08）：HP 4 群 6 個 pairwise 距離，748,391 個區域；Paired 全部 AUC<0.50（反轉：germline ASM > somatic ASM）；TO 最高 0.579（<0.58 閾值）。ISM methylation 特徵空間耗盡。
- Beyond-AUC 7 方法綜合驗證 EXHAUSTED CONFIRMED（2026-04-09）：7 種統計方法 x 分層 x 25 特徵 x 748K 區域；純 methylation 特徵全部 AUC<=0.58；pooled OLS residualization 被揭露為 data snooping；唯一 positive = HPFineNGroups TO non-LOH low-AF AUC=0.72（somatic 異質性 marker）。ISM methylation 特徵空間正式 CLOSED。
- SEQC2 CNV zone-aware 過濾器 CLOSED（2026-04-10）：FP-rate 模式跨樣本不一致（CN_HighGain>CN_Normal 僅 4/7）；各樣本 zone-specific 平均 AUC<=0.641；Simpson's Paradox 被否決。CNV 並非特徵空間耗盡的根因；所有 zone-exclusion 取捨皆不可行。
- Coverage_Multiple GC 校正 + methylation-CN 驗證全部 NO-GO（2026-04-11）：GC delta-r=-0.0002（閾值 >=0.03）；ONT 5kHz GC bias 極小；所有 HP-free 特徵 residualized |r|<0.07（CN-blind）；HPFineNGroups 68% 是 NumReads confound。不需 GC 校正；methylation 無法驗證 CN。（附註：KDE 確認正確，CN 準確度經 --expected-coverage 由 6.2%->43.8%。）
- 文獻假說交叉驗證（L1-L4）多數 NEGATIVE（2026-04-12）：60+ 篇論文、340K 區域 x 7 樣本。L1 directional ASM（epiTRACERx）：TP signed_delta 隨機 p=0.854；L2 PMD：TP/FP AUC<=0.622；L4 fCpG（EVOFLUx）：TP/FP CpG 變異相同 p=0.77。文獻 positive 源自任務差異（tumor-vs-normal != TP-vs-FP variant）。（L3 Normal baseline 可行 -> Phase 2A，並非 NEGATIVE。）
- Zone-Aware Confidence Framework：Characterization POSITIVE 但 F1 NEGATIVE（2026-04-17）：QS 模擬 5 種 delta 配置 x 21 閾值 x 7 樣本，最大 deltaF1=+0.001；根因為 TO QS AUC=0.497（隨機）。價值僅在 characterization 標註，而非 F1 改善。
- P3 Window Aggregation Pilot NEGATIVE（2026-04-17）：1Mb window 聚合；naive 4/7 樣本看似突破（H2009 deltaAUC +0.342），但在 mid-TP-rate confound 檢查下反轉（H2009 -0.346）；僅 1/8 TO 測試保留增益。AUC 幾乎完全由 TP/FP 基因組空間 auto-correlation 驅動。重啟需 shuffle-within-chr null + mid-TP-rate delta>+0.03 雙重 gate。
- P4 Second-Hit Order Pilot CONDITIONAL NEGATIVE（2026-04-17）：region-level AF x methylation Bimodality Coefficient 平均 |delta|=0.043（<<0.15 閾值）；HPFineNGroups delta=-0.637 是 LOH 副作用而非順序訊號。單一區域特徵摘要無法推論 two-hit 順序；需要 per-read epigenotype（ISM goal 1）。Goal 3 依賴 Goal 1；延後至 Goal 1 有進展。
- Subsample/Purity-Aware 過濾 NEGATIVE（item 11/12，2025-12~2026-02）：purity-aware 受 tumor-normal 組織差異混淆；subsample 無法可靠模擬 purity 效應。
- Phase 2 LOSO Sample-Level CV：Production 過濾器方向 FAILED（2026-05-20）：Cycle 1 +0.02236 in-distribution vs LOSO held-out -0.00012 = 100% 效應來自 sample-level circularity bias；5 樣本 LOSO 0/5 positive，Wilcoxon p=0.125，最佳 tau 退化為 0.10（保留全部）。論文 section3 撤回「ISM-augmented filter」claim -> 「ISM characterization + LR sample-level negative」。ISM characterization v0.3 cycle star3 維持不變。
- Phase 2 Cycle 3 ablation：ISM 殘留 + caller_af 主導（2026-05-20）：移除 5 個 methyl 特徵使 HCC1395 refit deltaF1 由 +0.02236 -> +0.02171（僅佔提升的 3%）；4/5 樣本 no-methyl >= full；ISM 是殘留 covariate（methylation = caller_af proxy）。過濾器更名為 caller-F1-headroom-gated 4-feature（caller_af + LOH_inner + Coverage_Multiple + NG）。
- Cross-sample / cross-Region subclone-consistency ML-ensemble 分類器關閉（backlog 附錄）：O13/O13v2 shared-read-count confound；Wave 3 J13 voting AUC=0.577<0.58；Beyond-AUC 2026-04-09 特徵空間耗盡。15 特徵的 ML ensemble 不可行。
- PMD/ChromHMM gating + CNV zone-aware + GC 校正全部列為 INDEX 附錄中已關閉、受 reopen 保護的方向（刪除線標示）：PMD CpG 變異 TP=FP AUC<=0.622 fCpG p=0.77；CNV 跨樣本不一致 AUC<=0.641；GC delta-r=-0.0002。
- [RETRACTED 2026-04-26] LOH x AF x CN 生物學知識導向分層過濾器：cross-sample FILTER 用途已 RETRACTED（X6 cross-sample：S3 TP>=0.85 僅 1/6 樣本，Wilcoxon S3>baseline p=1，原始 95.5% 為 stale-binary CN_tier artifact）。HCC1395 single-sample characterization 保留。

## 9. Provenance & Governance

- **來源**：7-thread fresh-context multi-agent survey（gather → synthesize → adversarial-verify），run `wf_ab45fe26-785`。
- **驗證**：verifier verdict **PASS**（filesystem-grounded，12 項獨立確認，0 fabrication，DAG 無環，全 7 thread 覆蓋）。
- **任務數**：76；**最終目標**：6；**blockers**：11；**護欄**：27。
- **Hard Gate**：C++ 改動（compile + commit）/ 研究 NO-GO / 刪檔 / evidence_ledger 覆寫 → 永遠需用戶確認。
- **快照性質**：合成時點快照；live ASM survey 完成後 ASM 相關數字（beta/n/rank）會自我更新，屆時刷新本檔。
