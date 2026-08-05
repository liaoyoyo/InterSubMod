<!--
建立時間: 2026-06-09
狀態: library (各文件庫 — 權威來源地圖 + 三色 provenance 對賬)
報告類型: paper_focus_document_library
受眾: 廖子游 · PI · 其他 AI session
provenance_note: provenance 分級來自 wf_a8ccbb34-3f7（8 組原檔 grep-back）；本檔不複製內容，只索引 + 標 tier。
-->
<!-- provenance-verified: 三色對賬來自 wf_a8ccbb34-3f7；本檔為來源索引。 -->

# 各文件庫 — 權威來源 + 三色 provenance 對賬

> **L0 一眼結論**：論文 headline 數字裡 **11 條已原檔對賬（🟢 可直接引用）**、**~10 條結論在但指定源檔本輪 grep 不到（🟡 投稿前必補）**。引用任何數字前，先在這查它是 🟢 還是 🟀。
>
> **L1 重點邏輯**：① 權威真值有順序（衝突時 ledger 為準）；② tsg 專案還在寫 → T2/T3/T5 是 moving target；③ 🟡 清單就是投稿前的 provenance 待辦（→ `05 任務樹 P-provenance`）。

---

## L2 — 權威來源地圖（要查 X 去哪）

| 要查 | 去哪 | 性質 |
|------|------|------|
| 論文就緒/可寫/矛盾/GO-NO-GO（最完整）| `InterSubMod/knowledge/11_external_literature/10_paper_readiness_convergence.md` | 收斂稽核 |
| 跨 cycle 歷史真理（append-only，衝突以此為準）| `InterSubMod/research/autoresearch/evidence_ledger.jsonl` | SoT |
| active cycle 快照 | `InterSubMod/state/active.json` | 機械態 |
| 本週 validated 數字 | `InterSubMod/docs/reports/validated/2026/06/20260602_…/master_draft.md` | validated |
| filter-DEAD 證據 | `InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/` | 🟢 stable |
| \|Δβ\| AUC | `InterSubMod/docs/experiments/in_progress/2026/05/20260531_ISM_aux_tag_observation_funnel_01.md` | 🟢 SoT-doc |
| copy/cis（dosage falsify）| `InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/fast_cnv_validation.json` | 🟡 live-tsg |
| AF↔NGroups（HD-4）| `InterSubMod/research/loh_subclone_af_paired/data/hd4_attribution.json` | 🟢 |
| 外部文獻地景（7 角度）| `InterSubMod/knowledge/11_external_literature/00-09` | L3 |

---

## L2 — 🟢 已原檔對賬（11 條，可直接引用）

filter 死四道（LOSO +0.02236→−0.00012 / mean −0.00004 p=0.125 / \|Δβ\| AUC=0.505 / COLO829 TP≈FP）· copy-dosage REFUTED（MW p=0.6183, ρ=−0.0829）· ASM TP 3.95%>FP 1.07% · phasing 7/7 W=28 p=0.0078125 · ARI 0.135<null + GNAS/RB1=1.000 · AF↔NGroups r=0.656=phasing · 同-locus 0/38 · paired 負控 NS · filter 補強三道。
> 完整逐條 + 檔:行 → 視覺版 `…20260609_G6_研究構想三層架構…standalone.html` §2 🟢 區。

## L2 — 🟡 原檔待對賬（投稿前必補 — 這就是 provenance 待辦）

> ✅ **2026-06-09 T-PROV 對賬結案**（wf_5644ed77-082 + 本地 grep；詳見 `01_focus_notes/08_任務執行回報_catalog建構與provenance.md`）：
> - **下表 🟡 已全部定位到源檔 → 升 🟢**。⚠ **doc 缺檔標記 stale**：`condition_fp_consensus.json` / `*_gwasm.json` **其實都在** `research/.../genome_survey_v2/cn_confound/cross_sample/`（OR 8.631/4.085/5.837、6/6 excess 0.101–0.241）。
> - **chr17 perm p=0.001** → `survivor_permutation.json`；**germline-het null ARI 0.177** → `b2_broad_scan_results.json`(0.177463)；**umtag 0.8852/V10 0.979-0.866** → `20260531_methyl_phasing_A0_assets/`。
> - 🟢 **same-hap occupancy 校正（原「3/6 OVERSTATED」誤判已撤回 — 2026-06-20 raw 逐樣本重算）**：「same-hap ≥93% **6/6**」**正確**。occupancy = (n_sameHP1 + n_sameHP2) / Inner 總 = **0.932/0.990/0.988/0.965/0.983/0.970**（源 `research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv`，6/6 ≥93%）。先前「3/6 修正」誤把 **same_HP1 tp_rate**（0.959/0.939/0.759/0.429/0.932/0.920，HCC1954 0.429 = TP outlier）當 occupancy；ledger:95 的 HCC1395=0.840 來自另一 run，混 run 又混度量、連自身來源都不自洽。🔴 **度量必明標 occupancy(phasing-purity) ≠ tp_rate**，勿再混。⚠ ledger:95 撤回 + CURRENT_FOCUS same-hap 引用修正屬 **Hard-Gate（owner session 處理）**；引用行號已 drift（稽核標 :158，今實測在 :166），owner 須重新 grep 定位，**勿照抄 :137**。

| 數字 | 問題 | 要去哪補 |
|------|------|---------|
| HCC1954 transfer −0.377 | 在 cycle4 LOSO 檔 grep 不到（LOSO 值是 −8.24e-05）| cycle2 cross-sample transfer 輸出 |
| regression-to-extreme OR 8.63/4.09/5.84 + combo | `condition_fp_consensus.json` 磁碟不存在 | 定位正確 emitted 檔名 |
| strong-ASM OR=0.194 | 同上缺檔 | 同上 |
| chr17 perm p=0.001 | 不在 fast_cnv_validation.json | tsg survey 另一檔 |
| BRCA2 d_within −0.023（**subclone/copy 主導，% 不 robust**，勿寫 ~80%）| ✅ 在 `copy_partition_confirm.json`（d_copy −0.11/d_within −0.023）+ `survivor_permutation.json`（perm p=0.024, 19/19）| ✅ 已對賬 🟢 |
| germline-het null ARI 0.177 | master_draft 只寫 p=0.922 | 上游 TSV |
| 6/6 excess +0.101–0.241 | 不在 cross_sample_synthesis.json | *_gwasm.json |
| E[overlap]=0.16 | derived，未 file-stored | 補計算說明 |
| same-hap occupancy ≥93% 6/6（0.932/0.990/0.988/0.965/0.983/0.970）| ✅ raw 重算對賬（occupancy ≠ same_HP1 tp_rate）| ✅ `obs18_NG2_composition_by_sample.tsv` 🟢 |
| umtag 0.8852 / V10 0.979-0.866 / 45.84% | 不在已驗證集合 | selfphasing_v6 / methyl-phasing pilot |

---

## L1 — 提醒
- 引用 🟡 數字前**必標「待原檔對賬」**，不可當 🟢 L1 用（否則 PI/reviewer 對賬會抓）。
- tsg 專案 06-07 仍在寫 → 投稿前以該專案定稿為準。
