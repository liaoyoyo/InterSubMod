<!--
建立時間: 2026-06-01
目標: Cleanup 稽核 Area [2/6] — output 現役區（canonical baseline + synthesis 非 archive + 頂層 standalone）
處理範圍（皆在 /big7_disk/liaoyoyo2001/big7_disk_output/ 下）:
  - canonical/（per-sample du，不進 per-region；含 longphase_s/*_tagged.bam 大檔）
  - synthesis/（concluded / observation_workspaces / research_rounds(active) / kde_rerun_* / manifests / master_report）
  - 頂層 standalone（multilayer_hp_benchmark / hcc1395_normal_pilot{,_global} / kde_smoke_test / v5_provenance_followup / InterSubMod_big7_runbook）
方法: 唯讀 du -sh(per-dir timeout) / ls -la / scoped find(-maxdepth + timeout + ALLOW_FULL_SCAN) / stat / 讀 README/INDEX/manifest。絕不刪除/搬移。
排除（屬其他 area）: bip8_output_archive / big8_output_archive（= area 3 archive 區）。
-->

# Area [2/6]：output 現役區（canonical + synthesis 非 archive + 頂層 standalone）

> **scan_status**: PARTIAL_COMPLETE — 所有 artifact 已定位 + 信任度/結論/verdict 判定完成。**精確 du 在 canonical 與 2 個 synthesis 大目錄（observation_workspaces / research_rounds）整樹會 timeout**（per-region ISM CSV 子目錄數十萬個 + 3.37 TB tagged BAM），故大檔改用 targeted `stat`/`find -maxdepth` 取真實 byte，per-region 分析樹標 `timeout(large)` bytes=-1。
>
> **本區頭條（且 area 6 漏掉的重大發現）**：`canonical/` 內藏 **14 個 `longphase_s/*_tagged.bam` = 3.37 TB**（complete_matrix runs），是整個 big7_disk_output 最大單一磁碟消費。area 6 報告（line 73/101）誤稱「canonical 全 ISM metadata、無 BAM >100MB、磁碟真實 tagged BAM 只有 synthesis 3 個」— **實測錯誤**：canonical paired BAM 才是大宗。這 14 個是 paper LIVE baseline → **KEEP**（per-sample 一份）；舊/重複/smokecheck run（同 sample 多版）標 NEEDS_USER_DECISION。

---

## 0. 關鍵發現摘要（先給結論）

| # | artifact | 大小 | verdict | 一句理由 |
|---|----------|------|---------|---------|
| 1 | `canonical/*/paired_{full,pileup}/*_complete_matrix/longphase_s/*_tagged.bam`（14 檔）| **3.37 TB** | **KEEP**（baseline 本體）| paper LIVE baseline 的 paired tagged BAM；canonical KEEP 主軸；重生需重跑 longphase-s + source BAM。**area 6 漏計** |
| 2 | `canonical/HCC1395/paired_full/` 等同 sample 多餘 run 目錄（20260211/20260307/_1/_2/no-suffix/20260420×3 等，多為 metadata-only 或空 longphase_s）| 小（KB-MB 級/目錄）| NEEDS_USER_DECISION | 同 sample-mode 多次重跑，僅 `_complete_matrix` 帶 BAM；其餘是早期/重複 ISM run 殘留，可清但需確認非引用 |
| 3 | `synthesis/research_rounds/legacy_partials/.../tumor_tagged.bam`（260G du）| 278.3 GB BAM | **→ area 6 owns** | TO legacy partial tagged BAM；**已由 area 6 列管**（本區只計其非 BAM = README 894 B）|
| 4 | `synthesis/research_rounds/20260315_hcc1395_to_pilot`（BAM→area6）| 非 BAM 996.7 MB | KEEP | HCC1395 TO 完整 run 分析產物（legacy README 指此為權威）；BAM(278.3 GB)由 area 6 列管 |
| 5 | `synthesis/research_rounds/20260423_colo829_to_pilot`（BAM→area6）| 非 BAM 914.4 MB | KEEP | COLO829 TO pilot 分析產物；撐 LIVE LOH-phasing n7 Grade A；BAM(95.4 GB)由 area 6 列管 |
| 6 | `synthesis/kde_rerun_B_14combos/`（per-region 大）| timeout(large) | ARCHIVE | KDE 全量重跑 14 組合，撐 validated 「KDE Fix Acceptance」結論；aggregate_summary + all_region_rows_kde_B 為證據本體 |
| 7 | `synthesis/observation_workspaces/`（57 子目錄 + 索引）| ~3.9 GB（已知部分）| KEEP/ARCHIVE 混 | O1-O13 系統性觀察證據；3× loh_round1 ~620M 各（含 master dataset）KEEP；NEGATIVE 結案 O 系列 ARCHIVE |
| 8 | `synthesis/concluded/`（6 子目錄）| 202M | KEEP | 已 concluded 一次性分析（chr19 verify / clairsto / loh）；撐 self-phasing + LOH 結論，收尾證據 |
| 9 | `kde_smoke_test/`（chr19 smoke，per-region 大）| timeout(large) | NEEDS_USER_DECISION | KDE chr19 smoke（BATCH COMPLETE）；已被 kde_rerun_B 全量取代 + 餵 validated acceptance；smoke 屬 TRANSIENT 但撐過渡證據 |
| 10 | `synthesis/kde_rerun_pilot/`（per-region 大）| timeout(large) | NEEDS_USER_DECISION | KDE pilot（單樣本 HCC1395_to_tp），已被 kde_rerun_B 14-combo 取代；過程檔 |
| 11 | `synthesis/research_rounds/20260423_colo829_to_kde_full/`（pipeline.log 13.7 KB only）| 13.7 KB | SAFE_DELETE | stub log-only 目錄；實際輸出寫進 sibling `20260423_colo829_to_pilot`，此處只剩孤兒 log |
| 12 | `hcc1395_normal_pilot{,_global}/`（BAM→area6）| 非 BAM ~5 MB | NEEDS_USER_DECISION | normal-read FP pilot（CONDITIONAL NO-GO）；subset BAM 1.84 GB 由 area 6 列管；VCF/regions 小 |

> **本區可回收量（扣除 area 6 已列管 BAM 後）≈ 0**：14 個 canonical BAM（3.37 TB）= KEEP（baseline）；synthesis 分析產物全是低成本但撐結論的證據（KEEP/ARCHIVE）。唯一乾淨 SAFE_DELETE = `colo829_to_kde_full` stub log（13.7 KB，可忽略）。真正釋放空間的決策（canonical 多餘 run BAM / synthesis TO BAM）全 NEEDS_USER_DECISION，且 3 個 synthesis 大 BAM 已歸 area 6 計。

---

## 1. 完整 artifact 表

> 欄位：path | 人類可讀 | bytes(timeout=-1) | purpose | trust_tier | conclusion_status | verdict | reclaimable_bytes | referenced_by | rationale
> 註：標「→area6」者 BAM byte 不在本區 reclaimable 計（避免雙重計）；本區只計其非 BAM 部分。

### 1.1 canonical/（baseline — 原則 KEEP / PRE-FIX 標記）

| path | 大小 | bytes | purpose | trust_tier | conclusion | verdict | reclaimable | referenced_by | rationale |
|------|------|-------|---------|-----------|-----------|---------|------------|---------------|-----------|
| `canonical/`（整體 19 canonical runs index）| timeout(large) | -1 | 7 樣本 × 3 模式 ISM baseline；CURRENT_FOCUS 主軸資料 | PRE-FIX (`missing_tagged_bam` per manifest) | LIVE | **KEEP** | 0 | OUTPUT_INDEX §4、master_run_manifest.tsv、`docs/standards/20260314_big7_canonical輸出規範`、近乎所有 analysis script | 重生成本=重跑 C++ ISM pipeline；paper baseline 必保 |
| `canonical/*/paired_{full,pileup}/*_complete_matrix/longphase_s/*_tagged.bam`（**14 檔**）| **3.37 TB** | 3709322840333 | longphase-s paired-mode tagged BAM（HP tag 來源）；HPRatio/HPSig 等 HP 欄位的原始載體 | PRE-FIX (HP；self-phasing 修正前) | LIVE | **KEEP** | 0 | self-phasing 驗證、LOH-phasing 主軸、HP 欄位全部分析 | **area 6 漏計**；單檔 95-453 GB（HCC1937 最大 452.6+447.2 GB）；canonical 主軸本體，per-sample 一份必保；可重生但重跑成本極高 |
| `canonical/{sample}/paired_*/` 內**非 complete_matrix** run 目錄（smokecheck / parallel_test / _1 / _2 / no-suffix / 20260211 / 20260307 / 20260420×N / 20260421）| 多為 KB-MB（多數 longphase_s 空、無 BAM；如 smokecheck=32K）| -1（per-dir 小但數量多）| 同 sample-mode 早期/重複/平行測試 ISM run | PRE-FIX | TRANSIENT / SUPERSEDED（被 complete_matrix 取代）| **NEEDS_USER_DECISION** | 視 user 確認（多數 metadata-only，BAM 只在 complete_matrix）| OUTPUT_INDEX 只列 complete_matrix 為 canonical run | 同 sample 多餘 run（HCC1395 paired_full 有 10 個 run 目錄，僅 1 帶 BAM）；保守標 NEEDS_DECISION，可清但需確認無 script 硬路徑引用 |
| `canonical/README.md` | 2 KB | 2023 | canonical 目錄説明 | CURRENT | LIVE | KEEP | 0 | — | 索引 |

> ⚠ **canonical BAM trust_tier 修正點**：manifest 標 `tagged_bam_ready=false` + master_report.md 稱「canonical bundles only materialize manifests and small summary files」— 但實測 14 個 paired tagged BAM 物理存在（3.37 TB）。manifest/master_report 與磁碟現實不一致（見 §3 anomalies）。

### 1.2 synthesis/ — manifests & 全域檔（KEEP）

| path | 大小 | bytes | purpose | trust_tier | conclusion | verdict | reclaimable | referenced_by | rationale |
|------|------|-------|---------|-----------|-----------|---------|------------|---------------|-----------|
| `synthesis/master_run_manifest.tsv` | 7.7 KB | 7856 | 19 canonical run 狀態追蹤 SoT | CURRENT | LIVE | KEEP | 0 | OUTPUT_INDEX、master_report | run 登記 SoT |
| `synthesis/master_experiment_matrix.tsv` | 6.2 KB | 6387 | 樣本×特徵×結論矩陣 | CURRENT | LIVE | KEEP | 0 | OUTPUT_INDEX | 實驗矩陣 |
| `synthesis/archive_synthesis_manifest.tsv` | 15 KB | 15376 | bip8/big8→big7 遷移對照 | CURRENT | LIVE | KEEP | 0 | OUTPUT_INDEX | 遷移對照表 |
| `synthesis/master_report.md` | 2.7 KB | 2734 | synthesis 總結（含 canonical 計數）| CURRENT（內容有過時宣稱）| LIVE | KEEP | 0 | OUTPUT_INDEX | 總結；但「metadata only」宣稱與 §1.1 BAM 現實衝突 |

### 1.3 synthesis/concluded/（一次性已收尾分析 — KEEP，202M 總）

| path | 大小 | bytes | purpose | trust_tier | conclusion | verdict | reclaimable | referenced_by | rationale |
|------|------|-------|---------|-----------|-----------|---------|------------|---------------|-----------|
| `concluded/chr19_full_refactor_verify/` | 187M | 196083712 | chr19 重構驗證完整版（self-phasing 修正驗證）| CURRENT | CONCLUDED_POSITIVE | KEEP | 0 | self-phasing 根因報告 | 收尾證據；撐 self-phasing 結論 |
| `concluded/chr19_refactor_verify/` | 20K | 20480 | chr19 重構驗證初版 | CURRENT | CONCLUDED_POSITIVE | KEEP | 0 | 同上 | 小；收尾 |
| `concluded/chr19_refactor_verify_new/` | 232K | 237568 | chr19 重構驗證新版 | CURRENT | CONCLUDED_POSITIVE | KEEP | 0 | 同上 | 小；收尾 |
| `concluded/clairsto_correction_analysis/` | 6.8M | 7130112 | ClairS-TO VCF 矯正分析 | CURRENT | CONCLUDED | KEEP | 0 | experiments #180、TO 特徵研究 | TO 校正證據 |
| `concluded/loh_verification_analysis/` | 872K | 892928 | LOH AlleleDelta 深度驗證 | CURRENT | CONCLUDED_NEGATIVE | KEEP | 0 | experiments #175、LOH 研究 | LOH 證據（O12 collider 教訓）；LIVE 主軸鄰接 |
| `concluded/to_loh_dual_definition_report/` | 7.7M | 8074035 | TO LOH 雙重定義報告 | CURRENT | CONCLUDED | KEEP | 0 | experiments #176-178、LOH 研究 | LOH 定義證據；LIVE 主軸鄰接 |

### 1.4 synthesis/observation_workspaces/（O1-O13 系統性觀察 — KEEP/ARCHIVE 混）

> 57 個分析子目錄 + OBSERVATION_INDEX.md + README。多數 KB-MB；大者為含 master dataset 的 loh_round1。NEGATIVE 結案 O 系列降 ARCHIVE，主軸鄰接保 KEEP。

| path | 大小 | bytes | purpose | trust_tier | conclusion | verdict | reclaimable | referenced_by | rationale |
|------|------|-------|---------|-----------|-----------|---------|------------|---------------|-----------|
| `20260327_loh_round1_cross_sample_audit/` | 621M | 651165696 | LOH round1 跨樣本 audit（含 `all_region_rows.tsv.gz` master dataset）| CURRENT | LIVE | **KEEP** | 0 | README「Master Dataset」、LOH 主軸全部分析 | **master dataset 入口**（post-fix）；LIVE 主軸命脈 |
| `20260330_loh_round1_cross_sample_audit_post_to_hp_fix/` | 621M | 651165696 | 同上 TO HP fix 後 | CURRENT | LIVE | KEEP | 0 | LOH 主軸 | post-fix master dataset |
| `20260327_loh_round1_cross_sample_audit_before_hp_fix/` | 620M | 650117120 | 同上 HP fix 前對照 | PRE-FIX | SUPERSEDED（被 post-fix 取代）| **NEEDS_USER_DECISION** | ~620M（若 user 確認 post-fix 足夠）| — | before-fix 對照；可清但屬 fix 前後對比基準，保守標 |
| `20260401_germline_fp_identification/` | 200M | 209715200 | Germline FP 識別 G1-G7 | CURRENT | CONCLUDED_NEGATIVE (NO-GO) | ARCHIVE | 0（撐 NO-GO 結論）| memory `germline_fp_identification_nogo` | TO Germline FP NO-GO 證據本體 |
| `20260322_to_fp_provenance_analysis_before_hp_fix/` | 113M | 118489088 | TO FP provenance（HP fix 前）| PRE-FIX | CONCLUDED_NEGATIVE | NEEDS_USER_DECISION | ~113M | — | before-fix 版，被 post 版取代 |
| `20260330_to_fp_provenance_analysis_post_hp_fix/` | 113M | 118489088 | TO FP provenance（HP fix 後）| CURRENT | CONCLUDED_NEGATIVE | ARCHIVE | 0 | memory `to_fp_provenance` | TO FP 無法 ISM 過濾證據 |
| `20260319_to_fp_normal_solved_analysis/` | 106M | 111149056 | TO FP normal-solved 分析 | CURRENT | CONCLUDED | ARCHIVE | 0 | — | 收尾證據 |
| `20260327_loh_round1_cross_sample_audit_smoketest_h1437/` | 104M | 109051904 | loh_round1 smoke（單樣本 H1437）| PRE-FIX | TRANSIENT | NEEDS_USER_DECISION | ~104M | — | smoke test，被全量取代 |
| `20260317_cross_sample_methylation_observation_workspace/` | 47M | 49283072 | 跨樣本甲基化觀察 | CURRENT | CONCLUDED | ARCHIVE | 0 | O5 系列 | 甲基化觀察證據 |
| 其餘 ~48 個 O 系列子目錄（O01-O13、allele、phasing_causal_chain、beyond_auc 等）| 各 KB-30M | -1（逐一小）| O 系列系統性觀察 / 因果鏈 / ablation | CURRENT/PRE-FIX 混 | 多 CONCLUDED_POSITIVE/NEGATIVE | KEEP（主軸）/ ARCHIVE（NEGATIVE 結案）| 0 | OBSERVATION_INDEX.md、experiments INDEX | 82 圖 36 報告證據庫；單檔小，整體撐多條已發佈結論 |
| `20260408_fine_pairwise_analysis/` | timeout(large) | -1 | Fine-pairwise 距離分析（per-region 樹）| CURRENT | CONCLUDED_NEGATIVE | ARCHIVE | 0 | memory `fine_pairwise_negative` | 6 pairwise 距離全無效；特徵空間耗盡證據 |
| `20260408_O09_fn_characterization/` | timeout(large) | -1 | FN 特徵化（per-region 樹；maxdepth2=9.9M 但深層大）| CURRENT | CONCLUDED_NEGATIVE (NO-GO) | ARCHIVE | 0 | memory `O9_fn_characterization_nogo` | FN≡TP in methylation space；NO-GO 證據 |
| `OBSERVATION_INDEX.md` + `README.md` | 小 | -1 | observation 索引 | CURRENT | LIVE | KEEP | 0 | OUTPUT_INDEX | 索引 |

### 1.5 synthesis/research_rounds/（active 研究 round — KEEP，含 area6 BAM）

| path | 大小 | bytes | purpose | trust_tier | conclusion | verdict | reclaimable | referenced_by | rationale |
|------|------|-------|---------|-----------|-----------|---------|------------|---------------|-----------|
| `legacy_partials/20260307_hcc1395_to_pilot_1_tagged_bam_only/` | 260G（du）| **非 BAM 894 B** | 保留舊 20260307 TO pilot 的 `tumor_tagged.bam` 供歷史回查 | SUPERSEDED | CONCLUDED (legacy partial) | **→area6**（BAM）/ KEEP（README）| 0（BAM 計 area6）| area6 報告、`docs/.../20260307_*TO_pilot*`、`docs/standards/20260314_*` | README 自述「非完整 bundle，不應當主結果入口」；278.3 GB BAM 由 area 6 列管 NEEDS_DECISION |
| `20260315_hcc1395_to_pilot/` | timeout(large)（含 278.3 GB BAM→area6）| 非 BAM 1045165530 (996.7 MB) | HCC1395 TO「完整」run（legacy README 指此為權威來源）| PRE-FIX | CONCLUDED (foundational TO run) | **KEEP**（分析）/ →area6（BAM）| 0 | legacy README、INDEX | TO 主 pilot；分析產物保留，BAM area 6 列管 |
| `20260423_colo829_to_pilot/` | timeout(large)（含 95.4 GB BAM→area6）| 非 BAM 958814353 (914.4 MB) | COLO829 TO pilot（ΔF1=+0.0001 NEGATIVE，但撐 LIVE phasing）| PRE-FIX | LIVE-adjacent | **KEEP**（分析）/ →area6（BAM）| 0 | `20260530_LOH_phasing_n7_Grade_A`(LIVE)、`20260423_COLO829_TO_Append_Plan`、INDEX | COLO829 撐 LIVE LOH-phasing n7 Grade A + ASM ⭐4 |
| `20260423_colo829_to_kde_full/` | 20K（只 pipeline.log 13.7 KB）| 14080 | COLO829 TO KDE full run **stub** | NA | TRANSIENT | **SAFE_DELETE** | 14080 | — | 孤兒 log：實際輸出寫進 sibling `colo829_to_pilot`，此目錄只剩 log；無 docs 引用 |
| `20260325_phase1a_*_sample637_*`（3 個）+ `20260325_phase1_*manifest*`（2 個）| 20K-96M | -1（各）| Phase 1A 最終 benchmark / 訓練 manifest / shard export（637 樣本）| CURRENT | LIVE | KEEP | 0 | experiments INDEX #44-45 | Phase 1A 最終版；訓練數據 manifest |
| `20260328_phase1a_round3_loh_feature_v1/` | 340K | 348160 | LOH 特徵 Phase 1A 第三輪 | CURRENT | LIVE | KEEP | 0 | — | LOH 特徵分析 |
| `20260330_loh_round2_ps_export_and_to_block_audit/` | 232M | 243269632 | LOH round2 PS export + TO block audit | PRE-FIX | LIVE | KEEP | 0 | — | LOH round2 證據 |
| `20260330_loh_round2_ps_export_..._post_to_hp_fix/` | 232M | 243269632 | 同上 TO HP fix 後 | CURRENT | LIVE | KEEP | 0 | — | post-fix LOH round2 |
| `20260330_to_hp_fix_reanalysis/` | 199M | 208666624 | TO HP fix 重新分析 | CURRENT | LIVE | KEEP | 0 | — | HP fix 重分析證據 |
| `20260330_to_loh_corrections/` + `20260330_to_loh_enrichment_post_hp_fix/` | 468K + 416K | -1 | TO LOH 校正 / enrichment | CURRENT | LIVE | KEEP | 0 | — | LOH 校正證據 |
| `20260314_hcc1395_to_pilot/`（32K）+ `20260315_hcc1395_to_pilot`(見上) | 32K | 32768 | 原始 TO baseline pilot | PRE-FIX | SUPERSEDED | NEEDS_USER_DECISION | ~32K | — | 早期 TO pilot；被後續取代 |
| `archive/`（34 dirs：duplicates/early_pilots/phase1a_iter/smoke）| timeout(large) | -1 | 被取代的歷史 research round | SUPERSEDED | DEAD/SUPERSEDED | **→ area 3/archive 區**（屬 research_rounds/archive，非本區現役）| 0 | archive/INDEX.md | INDEX 標 duplicates/ = SAFE_DELETE；建議交 archive area 處理 |
| `hp_fix_rerun_all.log`（70K）+ `README.md` | 小 | -1 | HP fix 批次重跑日誌 | CURRENT | LIVE | KEEP | 0 | — | 重跑日誌；索引 |

### 1.6 synthesis/kde_rerun_*（KDE 修正驗證 — ARCHIVE/NEEDS_DECISION）

| path | 大小 | bytes | purpose | trust_tier | conclusion | verdict | reclaimable | referenced_by | rationale |
|------|------|-------|---------|-----------|-----------|---------|------------|---------------|-----------|
| `kde_rerun_B_14combos/` | timeout(large)（per-region；top files 含 aggregate + all_region_rows_kde_B 10.5M）| -1 | KDE 全量重跑 14 樣本×模式組合（expected_coverage 75.0 hardcoded bug 修正後）| CURRENT | CONCLUDED_POSITIVE | **ARCHIVE** | 0 | `20260420_KDE_Fix_Acceptance_Validation`(validated)、`20260424_X5_CrossSample_obs18_KDE_Verified` | 撐 validated「KDE Fix Acceptance」結論的全量證據本體；aggregate_summary 是 SoT |
| `kde_rerun_pilot/`（HCC1395_to_tp）| timeout(large) | -1 | KDE pilot（單樣本，全量前的試跑）| CURRENT | TRANSIENT（被 14combos 取代）| **NEEDS_USER_DECISION** | per-region 大（待 user 確認可重生）| `20260419_KDE_*acceptance` | 單樣本 pilot；已被 kde_rerun_B 全量取代；過程檔 |

### 1.7 頂層 standalone

| path | 大小 | bytes | purpose | trust_tier | conclusion | verdict | reclaimable | referenced_by | rationale |
|------|------|-------|---------|-----------|-----------|---------|------------|---------------|-----------|
| `multilayer_hp_benchmark/` | timeout(large) | -1 | Self-phasing 修正前後 HP tag 多層比較 benchmark | PRE-FIX（benchmark 目的即比較 fix 前後）| LIVE（仍被引用）| **KEEP** | 0 | `docs/reports/research_landscape/02_Self_Phasing根因.md`、OUTPUT_INDEX §7 | self-phasing 驗證仍引用；benchmark 本體 |
| `hcc1395_normal_pilot/` | 3.0G（含 1.84 GB subset BAM→area6）| 非 BAM 4900642 (4.7 MB) | normal-read FP pilot（subset BAM + paired/TO filter VCF）| PRE-FIX | CONDITIONAL NO-GO | **NEEDS_USER_DECISION** | 0（BAM area6；VCF 小）| memory `read_level_germline_fp`、`normal_bam_progressive_copy` | normal pilot；CONDITIONAL NO-GO；subset BAM 由 area 6 列管 |
| `hcc1395_normal_pilot_global/` | timeout(large)（TO global VCF + logs，無 >100M BAM per area6）| -1 | normal pilot global 版（TO only）| PRE-FIX | CONDITIONAL NO-GO | NEEDS_USER_DECISION | 待 user 確認 | 同上 memory | global 版；與 hcc1395_normal_pilot 關係待 user 釐清（重複/取代？）|
| `kde_smoke_test/` | timeout(large)（batch + x1_archive_to_rerun，per-region）| -1 | KDE chr19 smoke test（BATCH COMPLETE）| CURRENT | TRANSIENT（被 kde_rerun_B 取代）| **NEEDS_USER_DECISION** | per-region 大（待 user 確認）| `20260423_KDE_Smoke_Test_Cross_Sample_Validation` | chr19 smoke；餵 validated acceptance 後被全量 14combos 取代；TRANSIENT 但撐過渡證據 |
| `v5_provenance_followup/chr19_pilot_20260506/` | 328K | 335872 | V5 provenance audit chr19 單點 pilot | CURRENT | LIVE | KEEP | 0 | memory `v5_somatic_fallback_verification`、commit 388a437 | V5 audit 證據；小 |
| `InterSubMod_big7_runbook/` | timeout(large) | -1 | big7 遷移執行紀錄（遷移手冊 + 狀態 + scripts + manifests）| CURRENT（標籤）/ 1.5 月未動 | CONCLUDED（遷移完成）| **NEEDS_USER_DECISION** | 待 user 確認（遷移若已 finalize 可降 SUPERSEDED）| OUTPUT_INDEX §7 | 遷移紀錄；歷史價值高，但遷移已完成 → user 決定是否降級封存 |
| `20260307_hcc1395_to_pilot_1`（symlink → legacy_partials）| 0（symlink）| 0 | 舊路徑相容 symlink | NA | NA | KEEP | 0 | 舊 docs/scripts | symlink 不佔空間；維持舊引用相容 |
| `OUTPUT_INDEX.md` + `README.md`（頂層）| 15.7K + 1.8K | 17546 | output 唯一入口索引 + 導航 | CURRENT（但 §5.6 列的 final_closeout 已不存在）| LIVE | KEEP（建議更新）| 0 | README、所有 area | output SoT 索引；需更新（見 §3 anomaly）|

---

## 2. Totals（本區，已扣除 area 6 列管 BAM 避免雙重計）

| 類別 | bytes | 人類可讀 | 説明 |
|------|-------|---------|------|
| **KEEP** | 3,713,000,000,000 (約) | ~3.37 TB | 主體=14 canonical tagged BAM 3.37 TB（baseline）+ synthesis 分析產物（concluded 202M / loh_round1 master ~1.86G / research_rounds active ~1.1G / observation O 系列 ~1.5G）|
| **ARCHIVE** | ~530,000,000 | ~530 MB（可量部分）| germline_fp 200M + to_fp_provenance_post 113M + cross_sample_methyl 47M + kde_rerun_B(per-region timeout) + NEGATIVE 結案 O 系列（fine_pairwise/O09 timeout）|
| **SAFE_DELETE** | 14,080 | 13.7 KB | 僅 `colo829_to_kde_full` stub log |
| **NEEDS_USER_DECISION** | ~2,200,000,000+（可量部分；多項 per-region timeout 無法精量）| ~2.2 GB+ | before_hp_fix 對照(620M+113M)、smoke(104M)、normal_pilot、kde_smoke/pilot(per-region 大)、runbook、canonical 多餘 run、各 timeout 大目錄 |
| **total_scanned（本區可量部分）** | ~3,716,000,000,000 | ~3.38 TB | 不含 area6 列管的 651.9 GB synthesis BAM + 1.84 GB normal subset；不含 per-region timeout 樹精確值 |

> ⚠ **bytes 精度限制**：canonical per-region ISM 樹、observation_workspaces 深層、kde_rerun per-region、archive 樹 du 全 timeout → 該等 reclaimable 標 -1/「per-region 大」。可精量的大檔（14 canonical BAM + 各 standalone 非 BAM）已用 `stat`/`find -maxdepth` 取真實 byte。

---

## 3. Anomalies（與先驗/索引/manifest 衝突處）

1. **🔴 area 6 漏計 3.37 TB canonical tagged BAM**：area_6_bam_transient.md line 73 稱「canonical/（19 runs，全 ISM metadata，`tagged_bam_ready=false`）...均無 BAM >100MB」、line 101「磁碟上真實 tagged BAM 只有 synthesis 這 3 個」。**實測錯誤**：`canonical/*/paired_*/*_complete_matrix/longphase_s/*_tagged.bam` 共 14 檔 = **3.37 TB**（單檔 95-453 GB），是 big7_disk_output 最大磁碟消費。area 6 淺掃（maxdepth 4-6）疑似未觸及 depth-6 的 longphase_s 或 find 提前 timeout → 假陰性「0 BAM」。**主 agent 加總時務必以本區 14 BAM 為準，勿信 area 6 的 canonical=0-BAM 宣稱。**
2. **manifest / master_report 與磁碟現實不符**：`master_run_manifest.tsv` 全 19 run 標 `tagged_bam_ready=false`；`master_report.md` 稱「canonical bundles only materialize manifests and small summary files; large BAM/VCF payloads remain at archive source paths」。但 14 個 paired tagged BAM（3.37 TB）物理存在於 canonical。manifest 欄位與實況脱鈎（可能 BAM 是後來 longphase-s 跑出但未回寫 manifest）。
3. **OUTPUT_INDEX §5.6 列的 `final_closeout/` 與 `final_closeout_debug/` 已不存在**：實測 `synthesis/` 下無此二目錄（2026-05-10 inventory proposal §2.2 仍列為「review/封存」候選）。已被清理或搬移但 OUTPUT_INDEX 未更新 → 索引 stale。
4. **OUTPUT_INDEX §5.4「batch_runs/（10 dirs）」已不存在於 synthesis/**：實測 `synthesis/` 下無 `batch_runs/`（任務指引也提到「batch_runs(10 PRE-FIX)」）。可能已搬入 archive 或刪除；OUTPUT_INDEX/任務先驗 stale。本區未見此目錄。
5. **2 個 colo829 新 round（20260423）未進 OUTPUT_INDEX**：`20260423_colo829_to_pilot` + `20260423_colo829_to_kde_full` 為 Apr 23 新增，OUTPUT_INDEX §5.2（2026-04-14 建立）未列。其中 to_kde_full 是 stub log。
6. **kde_rerun_B / kde_rerun_pilot / kde_smoke_test 未進 OUTPUT_INDEX**：三者皆 Apr 19-24 新增於 OUTPUT_INDEX（4-14）之後；2026-05-10 inventory proposal 有列為待 review。本區判 kde_rerun_B=ARCHIVE（撐 validated 結論）、pilot/smoke=NEEDS_DECISION（被取代的過程檔）。
7. **同 sample-mode 多版 canonical run**：HCC1395/paired_full 有 10 個 run 目錄（20260211/20260307/20260314×4/20260420×3），僅 `_complete_matrix` 帶 3.37TB 內的 BAM；其餘多 metadata-only/空 longphase_s。屬可清的重複 run，但保守標 NEEDS_USER_DECISION。
8. **legacy_partials BAM 與 20260315 BAM 非 byte-identical**（area 6 已註）：278306777912 vs 278305302826，差 1.47M → 兩次獨立 longphase run，非純副本；故 legacy 標 NEEDS_DECISION 而非 duplicate→SAFE_DELETE（保守）。

---

## 4. 跨 area 邊界聲明（避免雙重計）

- **area 6 (BAM/transient) 已列管**：synthesis 3 個 `tumor_tagged.bam`（278.3+278.3+95.4 = 651.9 GB）+ `hcc1395_normal_pilot` subset BAM（1.84 GB）。本區只計其**非 BAM** 部分（legacy README 894B / 20260315 996.7M / colo829_to_pilot 914.4M / normal_pilot 4.7M）。
- **本區獨家（area 6 漏）**：14 個 **canonical** `longphase_s/*_tagged.bam` = 3.37 TB。建議主 agent 把這 14 BAM 歸入本區計，且**修正 area 6 的 canonical=0-BAM 結論**。
- **archive 區（area 3?）**：`synthesis/research_rounds/archive/`（34 dirs）+ `bip8/big8_output_archive` 不在本區現役 scope，未深掃。

---

*稽核者: Area 2 subagent | 2026-06-01 | 唯讀，未刪除/搬移任何檔案*
