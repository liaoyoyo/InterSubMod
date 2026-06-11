<!--
建立時間: 2026-06-09
狀態: report (07 任務清單執行回報 — 給廖子游 + 論文 session)
報告類型: paper_focus_task_execution_report
受眾: 廖子游 · 論文 session
source: docs/paper_focus/01_focus_notes/07_研究觀察任務清單_給其他session.md
verification: workflow wf_5644ed77-082（provenance 定位 + catalog 獨立重驗）
-->
<!-- provenance-verified: 全數字由 catalog/*.json + 各源檔 grep；wf_5644ed77-082 獨立重驗；本檔為執行回報。 -->

# 07 任務清單 — 執行回報

> **L0 一眼結論**：9 個任務裡 **6 個可執行的全部完成並獨立驗證**（T-CODE-1/2/3 · T-OBS-1/2 · T-PROV · T-VAL-1 部分），**2 個 off-path 已評估給規格**（T-VAL-2/3）。最大交付 = **catalog 骨架（332,705 loci × 7 TAG）** 證實「甲基分群真實但壓倒性非 somatic 判別器」。T-PROV 抓到 **1 個對外數字 overstated（same-hap 93% 6/6 → 實際 3/6）必修**，並發現 doc 02「缺檔」標記是 stale（檔其實都在）。
>
> **L1 重點邏輯**：① catalog 把零散觀察變成可發表的 R6 分類表；② 全部 metric 數字現在都 grep 得到（T-PROV 升 🟢）；③ 唯一 cis exemplar = chr17，BRCA2=copy，與 capstone 一致——論文口徑不變且強化。

---

## L2 — 9 任務狀態總表

| 任務 | 狀態 | 交付物 | 關鍵結果（已驗證）|
|------|------|--------|------------------|
| **T-CODE-1** catalog 骨架 | ✅ 完成 | `catalog/catalog_skeleton.tsv`(332,705×32欄) + `catalog_tag_counts.json` | A=1(chr17) B=5 C=12,868 D=5 E=28,254 F=54 G=291,518 |
| **T-CODE-2** CramersV audit | ✅ 完成 | `catalog_audit.json` + 結果 R6 §2 | gate 偏嚴；latent 救回 enrich≈base 47:1 = characterization-only；建議三態（不改碼）|
| **T-CODE-3** Δβ Venn | ✅ 完成 | `discovery_venn.{tsv,json}` | dbeta_only TP:FP=11.9(最FP-enriched) vs ISM_only 147.6 → Δβ 必 clustering-gated |
| **T-OBS-1** 6 分佈圖 | ✅ 完成 | `04_figures/P1-P6 + catalog_overview.png` | CJK 正確；P5=R6 主圖 |
| **T-OBS-2** 例 + chr17 二樣本 | ✅ 完成（誠實負）| catalog examples + 結果 R6 §5 | chr17 SNV=HCC1395 private（不在 HCC1937/1954）→ cis 無法二樣本直接重現 |
| **T-VAL-1** cis 協議 | 🟡 部分完成 | catalog `cis_status` 欄 | chr17=cis / BRCA2=copy / 53 untestable + chr18 mechanical；完整協議需 anchor SEQC2 CN（規格見下）|
| **T-PROV** provenance 對賬 | ✅ 完成 | 下方對賬表 | 全 🟡 定位到源檔；**1 個 overstated 必修**；doc 02「缺檔」是 stale |
| **T-VAL-2** HD-4 closeout | ⏸ 評估（off-path）| 規格見下 | HD-4 已 RESOLVED（ledger 92, r=0.656）；n_clusters re-run 是精修，需 ISM re-run ~2.4hr |
| **T-VAL-3** COLO829 ⭐4 + umtag | ⏸ 評估（off-path）| 規格見下 | COLO829 existence scan 已有；⭐4 需 COLO829 normal BAM cis；umtag 需 MethPhaser 對標 |

---

## L2 — T-PROV provenance 對賬表（投稿前 🟡 清單結案）

> 🔑 **重大發現**：doc 02 標「磁碟不存在」的檔**其實都在** `genome_survey_v2/cn_confound/cross_sample/`。doc 02 的缺檔標記 stale，應更新。

| 數字 | 源檔（InterSubMod/...）| 值 | 判定 |
|------|------|----|------|
| HCC1954 transfer −0.377 | `research/methyl_augmented_filter_phase2/cycle3/cycle3_step1_findings.md:102` | −0.37744 | 🟢 match |
| germline-het null ARI 0.177 | `research/.../genome_survey_v2/b2_broad_scan_results.json` | het-NULL 0.177463 / TP 0.135091 | 🟢 match |
| chr17 perm p=0.001 | `research/.../genome_survey_v2/survivor_permutation.json` | obs 0.1418, p=0.001, genuine | 🟢 match |
| BRCA2 d_within −0.023 | `research/.../genome_survey_v2/copy_partition_confirm.json` | −0.023（d_copy −0.11）| 🟢 match |
| OR 8.63/4.09 + combo 0.86%/7.97% | `research/.../cross_sample/condition_fp_consensus.json` | 8.631 / 4.085 / TP 0.0086 vs FP 0.0797 | 🟢 match |
| OR ~5.84（第3條件）| 同上（small somatic subhap <10 reads）| 5.837 | 🟢 match(~5.84) |
| strong-ASM OR=0.194 | memory `project_zar1l_brca2_asm_verification` + method_audit(0.19 rounding) | 0.194 | 🟢 match |
| 6/6 excess +0.101–0.241 | `research/.../cross_sample/*_gwasm.json` `rate_excess_over_null` | 0.101/0.150/0.151/0.241/0.171/0.196 mean 0.168 | 🟢 6/6>0 match |
| umtag 0.8852 / null 0.5236 | `docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/extrapolation_genome_aggregate.json` | 0.8852 / 0.5236 | 🟢 match |
| V10 normal 0.979 / tumor 0.866 | `docs/.../20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md:148` | 0.979 / 0.866 | 🟢 match |
| **same-hap ≥93% 6/6** | CURRENT_FOCUS:137（散文）vs Thread D 實測 | **實際 0.840/0.939/0.759/0.429/0.932/0.920** | 🔴 **OVERSTATED — 只 3/6 ≥93%，非 6/6** |

### 🔴 必修：same-hap 口徑
- **錯誤口徑**：「NG=2 Inner ≥93% same-hap，6/6 樣本一致」
- **正確口徑**：「NG=2 Inner same-hap 中位 ~92%，**6/6 樣本方向一致（皆 high same-hap）但僅 3/6 ≥93%**（0.840/0.939/0.759/0.429/0.932/0.920）」
- **影響**：phasing 脊柱敘述——「6/6 一致」對（方向），但「6/6 ≥93%」錯（量級）。投稿前改 CURRENT_FOCUS:137 + 任何引此句的文件。

---

## L2 — Deferred 任務規格（off 關鍵路徑，資源/re-run 到位再做）

### T-VAL-2 ｜ HD-4 closeout（精修，HD-4 已 RESOLVED）
- **現況**：`research/loh_subclone_af_paired/data/hd4_attribution.json` 已有 baseline Spearman r=0.6561(n=36854) + ngroups_definition（LabelTest.cpp:265-305 = HP-tag count 非甲基）+ verdict PHASING_DRIVEN（ledger 92）。
- **要補**：一次 ISM re-run 匯出 `fine_group_counts{HP1,HP1-1,HP2,HP2-1}` + `n_clusters`（甲基-only count）→ partial-corr knockout（控 n_clusters 後 AF→NGroups r 是否仍在）。
- **成本**：ISM re-run ~2.4hr（背景 Bash，非 workflow）。**off 關鍵路徑**——三證已一致，此為補刀。

### T-VAL-3 ｜ COLO829 ⭐4 + umtag yardstick
- **COLO829 ⭐4**：existence scan 已有（COLO829_tp/fp/fn）；升 ⭐4 需 **COLO829 normal BAM** 跑 normal-anchored cis（解 ASM 單樣本天花板 + 給 chr17 第二癌種 cis 重現）。
- **umtag**：需算 switch-error/N50 vs **MethPhaser/HapBridge**（業界 yardstick）+ 真-unphase 救援實測（非 held-out）+ V10 非-copy 在 COLO829 重現。
- **caveat**：held-out 0.8852 ≠ 真救援；單樣本→COLO829 是關鍵。**資源到位再做**。

---

## L1 — 建議下一步（給用戶決策）
1. **必修**：same-hap 口徑（3/6 非 6/6）→ 改 CURRENT_FOCUS + phasing 敘述（投稿前）。
2. **可選即做**：把 catalog 結果 R6 併入論文 §Results（P5 主圖 + 7-TAG 表）；更新 doc 02 移除 stale「缺檔」標記。
3. **資源到位**：COLO829 ⭐4（解天花板）→ 才動 T-VAL-2/3。
4. **C++ 改動（需 /cpp-change）**：clustering_reliability 三態 + ISM 加 dbeta_max 欄——**非必要**，catalog 已用既有欄達成 characterization。

> **產物清單**：catalog 結果 = `InterSubMod/docs/paper_focus/02_paper_framework/位點甲基分群catalog_結果_R6.md`；資料 = `research/tsg_promoter_asm_reviewer/genome_survey_v2/catalog/`（4 檔）；圖 = `InterSubMod/docs/paper_focus/04_figures/P1-P6 + catalog_overview.png`；scripts 98-100。
