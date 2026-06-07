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

## 待辦
- LOH 純-LOH 分得開異常 IGV 深究（76% 分得開的 LOH 是 pureLOH，非 cnLOH）
- SEQC2 大染色體速度（預載甲基索引）
- 跨樣本 COLO829（真驗證 ASM×LOH 富集需低背景樣本）
- 評估回灌 longphase-S / ISM C++
