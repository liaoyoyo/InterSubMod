<!--
建立時間: 2026-06-03
報告類型: 開發紀錄 (dev log) — 甲基救 unphase 研究 pilot
任務類型: A pilot — partial: HCC1395 單樣本
-->

# 甲基救 unphase read — 開發紀錄

> 本目錄是「用甲基分群救 longphase-S unphase read」研究 pilot 的全部腳本 + 真值 + 圖。
> **所有數字的單一真值來源**：`VERIFIED_RESULTS.md`（V1-V6）。引用前對照 `DATA_PROVENANCE_LEDGER.md` 分級。

## 研究問題
longphase-S 把 45.84% read 標成 unphase（無相位）。能否用 per-read 甲基模式把它們救回正確 haplotype？

## 證據鏈（已建立，全程 HCC1395 paired 單樣本）
1. **機會空間**（V1）：全基因組 45.84% unphase（1382 萬 read，守恆 PASS）。
2. **甲基能分 HP**（V2-V5）：anchor AUC，過 shuffle / CN-confound(SEQC2) / CpG-SNP 三道控制。
3. **能真的救 unphase**（V6）：PS-block held-out 外推，全基因組 183 窗 accuracy median 88.5%（null 52.4%）。
4. **救援原型**（本檔）：對真 unphase read（無 HP tag）救援，chr20 smoke 87.5% 指派、守恆通過。

## 腳本清單（皆唯讀 BAM、可重跑）
| 腳本 | 用途 | 輸出 |
|------|------|------|
| `a0a_hp_distribution.sh` | 全基因組 HP 分布 | per_chr/*.tsv |
| `extract_per_read_methyl.py` | per-read×CpG 甲基矩陣 | *_matrix/meta.tsv |
| `heatmap_and_separation.py` | 甲基熱圖 + anchor AUC + null | *_separation.json + heatmap.png |
| `germline_het_null.py` | germline-het null AUC 分布 | germline_het_null_results.json |
| `shuffle_control.py` | shuffle-label 控制 | shuffle_control_results.json |
| `seqc2_cn_methyl.py` / `_v2.py` | SEQC2 CN/LOH × 甲基（v2 加 CpG-SNP 排除）| seqc2_cn_methyl(_v2)_chr*.json |
| `aggregate_*.py` | 彙總 | seqc2_*aggregate.json |
| `extrapolation_validation.py` | ★ PS-block held-out 外推救援驗證 | extrapolation_chr*.json |
| `unphase_rescue_prototype.py` | ★ 真 unphase read 救援原型 | unphase_rescue_chr*.json |
| `plot_*.py` | 圖 | fig_*.png |

## 開發決策紀錄
- **2026-05-31~06-03**：採 Python 旁路原型（用戶確認），**不改 ISM C++**（避 Hard Gate / 不可逆 / 重編譯）。方向驗證後再評估回灌 C++。
- **救援原型設計**（design_02 對齊）：不覆寫原 tag，輸出 `original_hp` + `reassigned_hp` + `reassign_source` + `reassign_conf`，可回溯；雙口徑（全 target vs 有甲基覆蓋）；守恆強制（assigned + 無法救 = 總 target）。
- **confidence 公式**：`margin × cov_factor`，margin=|d1-d2|/(d1+d2)，cov_factor=min(1, n_cpg/8)。

## 數據誠信
本 session 曾捏造數字 3 次（已全更正 + 立防線）；主動發現修正取樣 bug（重疊 window，421→317 去重）。完整事件：`InterSubMod/docs/postmortems/20260601_fabricated_metric_in_html_preview_postmortem.md`。所有報告數字 grep-verified。

## V9 開發紀錄（2026-06-03）— aDMR × LOH 富集對照文獻 79%
- **動機**：V8 引文獻「79% aDMR 落 CNV/LOH」推論 LOH 是 ASM 富集區。本輪不靠引用，**直接在本資料算**。
- **腳本**：`admr_loh_enrichment.py`（±2kb 窗算 HP1 vs HP2 |Δβ| + Mann-Whitney → 窗級 aDMR）+ `admr_aggregate.py`（Fisher OR）+ `plot_admr.py`（雙圖）+ `run_admr_genome.sh`（12 染色體平行）。
- **產出**：`admr_aggregate.json` / `admr_chr*.json` / `fig_admr_enrichment.png`。
- **關鍵自我修正**：原想用「絕對 % 對照 79%」→ 發現 HCC1395 背景 CNV/LOH 窗級覆蓋達 **96.5%**（chr 級 68-99%）→ 絕對 % 必然高、無區辨力 → **改用 Fisher OR**。
- **真值**（admr_aggregate.json）：1012 窗 / 816 aDMR(80.6%) / aDMR 落 CNV/LOH 96.8% vs 背景 96.5% vs 非 aDMR 95.4% / **OR=1.462 Fisher p=0.382 不顯著** / aDMR maxΔβ median 0.802 vs 非 0.269。
- **結論（誠實）**：(1) 數字對得上文獻 79%（甚至更高）；(2) 但為背景假象非真富集（OR≈1, p NS, 每染色體 aDMR≈背景）；(3) aDMR 本身是真甲基訊號（maxΔ 0.80 vs 0.27），只是不偏好 LOH 區。
- **方法學教訓**：「對照文獻數字相符 ≠ 證實機制」。要真驗證 ASM×LOH 富集需低背景樣本（近二倍體）或 matched 對照。
- **報告**：`20260603_loh_admr_enrichment_literature_check_01.standalone.html`（6-taboo PASS / 圖路徑驗 / 數字全 grep-verified）。

## V10 開發紀錄（2026-06-08）— 決定性 matched normal 對照：甲基 allele 差異「不是 copy」
- **動機**：用戶問「HP 群內外甲基差異是 copy 還是真 haplotype ASM，如何驗證」。最乾淨 = matched normal（copy-clean 二倍體）對照。
- **設計**：單 SNP REF-vs-ALT allele 標籤（不靠 longphase HP tag、不需相位 → tumor/normal 對等獨立）；normal HCC1395BL 二倍體無 CNV/LOH → 若 copy 造成分離 normal 應低。
- **腳本**：`allele_asm_auc.py`（tumor+normal 同區域配對，內含 depth-matched + CN status）+ `aggregate_allele_asm.py` + `imprint_bimodality.py`（正控）+ `plot_allele_asm.py`。背景跑 `run_allele_asm.sh`（6 chr × 2 樣本，無 error）。
- **產出**：`allele_asm_{tumor,normal}_chr*.json` + `allele_asm_aggregate.json` + `imprint_bimodality.json` + `fig_allele_asm_tumor_vs_normal.png`。
- **真值**（6 chr；tumor 638/normal 720/136 配對）：整體 tumor 0.866 vs **normal 0.979**（6/6 染色體 normal>tumor）；depth-matched 0.859/0.982（P-06 否證）；by status tumor neutral 0.854 / gain 0.871 / **loh 0.782 最低**，normal@同位點 0.786/0.989/0.932；配對 delta −0.046。imprinting GNAS 雙峰 Δ=0.49 n=464。
- **結論**：(1) ✅ copy 決定性排除（normal copy-clean 卻 ≥ tumor；tumor copy 事件反而降低分離）；(2) ✅ depth/P-06 否證；(3) ✅ 真 haplotype-linked copy-independent 訊號（V6 88.5% 根基）；(4) ⚠ 絕對 AUC 含方法樂觀（normal ~0.98 全基因組不合理）→ 誠實 effect 用相對 null（real vs shuffle tumor 0.87/0.67、normal 0.98/0.79）+ V6 88.5%；(5) ✅ imprinting 正控過關。
- **報告**：`20260608_allele_asm_normal_control_not_copy_01.standalone.html`（6-taboo PASS / 圖驗 / 數字全 grep-verified）。

## V11 開發紀錄（2026-06-09）— 甲基輔助 longphase-S tag 矯正三 target（T1/T2 ✅、T3 ❌）
- **動機**：用戶 — 用甲基弱數據輔助矯正 longphase-S：T1 unphase→H1/H2、T2 H3→H1-1/H2-1、T3 拆 plain H1 是否亞群(H1-1)。
- **先讀原始碼確認 cascade**（known-pitfall P-10/P-14，不猜 tag 名）：`HaplotagStrategy.cpp:452-600`(judgeSomaticReadHap)+`SomaticHaplotagProcess.cpp:461-527`(inheritHaplotype)，門檻 0.6。6-state：unTag / H1·H2(germline-only) / H1-1·H2-1(somatic 歸 germline 分支) / H3(somatic germline 未知)。plain H1=純germline沒踩somatic位點；H1-1=踩到且歸H1分支。
- **腳本**：T1=`assist_tag_separable.py`（held-out，train 量分離度/test 量 assist，**無循環**）；T2/T3=`methyl_assist_targets.py`（somatic 位點窗按既有 HP tag 分群算 LOO-AUC+shuffle）+`aggregate_assist_3targets.py`+`plot`。背景跑 run_assist_tag.sh / run_methyl_assist.sh。
- **真值**（assist_3targets_aggregate.json）：T1 可分位點 held-out **0.926**(conf_yield 0.794/conf_acc 0.974, 不可分 0.682, n=141)；T2 **0.900**(null 0.732, sig 0.729, n=436)；T3-H1 **0.652**(null 0.707, sig 0.359, n=1139)；T3-H2 **0.649**(null 0.702, n=1262)。
- **結論**：T1✅(分不同 haplotype) + T2✅(分 germline 分支) 可行；**T3❌ NEGATIVE**(AUC 0.65 **低於** shuffle null 0.71，甲基無法區分同 haplotype 內亞群vs母本)。
- **機制（一致）**：甲基訊號=germline-haplotype 層級(V10)→分不同 haplotype 強、分同 haplotype 內亞群弱；subclone 層級 ASM ≪ haplotype 層級。
- **報告**：`20260609_methyl_assist_longphase_tag_3targets_01.standalone.html`（6-taboo PASS / 圖驗 / 18 數字全 grep-verified）。

## V11b 開發紀錄（2026-06-09）— 三 target 嚴格反例驗證 + 4-agent 對抗審查（修正 V11）
- **動機**：用戶質疑 T3 反駁是否成立 + 要求 T1/T2 信心/反例/邊緣處理 + 多 agent 詳細檢視。
- **腳本**：`rigor_t1.py`(失敗特徵化)/`rigor_t2.py`(H3 germline 邊緣分類)/`rigor_t3.py`(Δβ+噪音地板+正控)+`aggregate_rigor.py`+`plot_rigor_corrected.py`。背景跑 run_rigor.sh(t1/t2/t3 × chr7/15/20/22)。
- **4-agent workflow** wf_b5b391cf-ea4（3 審查+1 綜整，462k tokens）→ **所有反駁獨立從 per-chr JSON 重現**，我再第三次獨立重現（§13.7）。
- **修正結論（取代 V11 T2/T3 措辭，見 VERIFIED_RESULTS V11b）**：
  - **T1 = SUPPORTED+caveats(0 blocking)**：held-out 真無循環；誠實 headline=V6 全基因組 0.885 非 0.926(僅可分子集)；🔴 失敗=低 train_sep 非低覆蓋；237/240 窗 gain、LOH 基本未測；abstain gate chicken-egg。⭐3。
  - **T2 = OVERSTATED(2 blocking)**：0.90 只在有真值 1-1/2-1；用甲基歸 H3 未驗證(H3 無 AUC/null、僅 15-18% 可指派、56% margin 無 null)；inconsistent germline 建議 unTag。⭐2-3。
  - **T3 = NEGATIVE 維持**：我的「+0.022 excess」是 unmatched-window artifact(已撤)；正確 matched 同窗配對 亞群−噪音 median −0.0043 **p=0.924**(=噪音)、正控存活 +0.0425 **p=1.77e-49**。用戶必要條件不滿足。
  - **機制天花板**：甲基=germline-haplotype 層級(V10)→分不同 haplotype 強(T1/T2)、within-haplotype 弱(T3)；LOH/subclone 救援情境恰是機制衰退處。
- **報告**：`20260609_methyl_assist_3target_adversarial_review_01.standalone.html`（含 matched 檢定箱型圖 + verdict + 共通風險）。

## V11c+V12 開發紀錄（2026-06-11）— T3 local-allele 窄翻案 + unphase/HP3 量清點
- **V12 unphase/HP3 量**（`unphase_inventory.py`，chr15/19/20/22 實掃）：unphase 13,821,877 中**僅 ~6%(~80萬)可嘗試救援**（有甲基+本地HP1/HP2錨點），**94% 無本地錨點無法救**（88.5% 是可嘗試池內正確率）。HP3 76,738：90% no_germline、可嘗試指派~1萬但無真值。
- **V11c T3 local-allele**（用戶質疑驅動）：改用該位點 ALT/REF 鹼基當亞群/母本(限同HP1)。`t3_local_allele.py`(GATE held-out AUC+噪音地板+farCpG+lean)+`t3_reconcile.py`(matched Δβ+覆蓋confound)。對抗 wf_5553cd83-ee9(2攻擊+1裁決,408k,雙重重現)→ **REVERSAL_OVERSTATED**：
  - **存在性窄翻案**：亞群vs母本 farCpG AUC 0.85 vs 噪音0.50、matched Δβ+0.04(3/4chr sig)、無覆蓋confound、非somatic-CpG → 先前 HP-tag NEGATIVE 部分是標籤污染(亞群read多不蓋稀疏突變被標「1」)。
  - **可用性NEGATIVE維持**：ambiguous lean 8/8 全<0.5(偏母本)、無ground truth。
  - 未解：farCpG±100bp太窄(cohesion≠cis)、單樣本⭐3未達reopen；Attack2 germline_het_null「fatal」被裁決長否決(跨-hap對照軸錯)。
- **housekeeping**：rigor_aggregate.json 的已撤 subclone_excess_over_noise 0.0224 標 RETRACTED。
- 複驗(§13.7)：lean 8/8<0.5 + germline_het_null=HP1/HP2 + housekeeping 全自查確認。

## 待辦
- LOH 純-LOH 分得開異常 IGV 深究（76% 分得開的 LOH 是 pureLOH，非 cnLOH）
- 🔬 T3 reopen 前置：≥1 樣本(COLO829) + within-haplotype somatic-vs-baseline 對照 + 更寬空間控制
- SEQC2 大染色體速度（預載甲基索引）
- 跨樣本 COLO829（真驗證 ASM×LOH 富集需低背景樣本）
- 評估回灌 longphase-S / ISM C++
