# 20260618 subcluster pilot — 資產與重現

> 標籤內子分群 / cluster×label「幾群」判定方法學 pilot。**HCC1395 單樣本 · 全基因組 30,490 TP SNV · tumor reads · ⭐2 偵測非驗證**。
> 主文件：`../20260618_within_label_subcluster_methodology_validation_01.md`。核心 memory：`project_subcluster_cluster_count_determination`。

> **2026-06-19~20 延伸（k-sweep 對齊式 k-selection + k-profile 切法品質 + tumor/merged + dual-panel）**：實驗紀錄 = `../../../experiments/in_progress/2026/06/20260618_ksweep_alignment_k_selection_wg_01.md`（觀察 A-D + 總結論）。腳本：`ksweep_wg.py`(全基因組 k-sweep)→`ksweep_analyze.py`(gate+FDR)→`ksweep_wg_v2.py`(tumor/merged dual-mode)→`ksweep_kprofile.py`(margin+三態)→`ksweep_fdr_dump.py`/`select_kprofile_examples.py`(mini-VCF)→`plot_{fdr,kprofile}_heatmaps.py`(**SoT dual-panel 甲基+距離**, import `scripts/ism_heatmap_std.py`)→`build_{ksweep,kprofile}_explainer.py`/`build_fdr_workstation.py`(HTML)。HTML：`20260618_{ksweep,kprofile}_explainer_01.standalone.html` + `..._workstation_01_lite.standalone.html`(+`figs_fdr/`)。大檔(records/loci/figs/log/vcf)gitignore，從 binary 5c39051 重生。

## HTML（自含，base64 圖，可直接開）
- `20260618_cluster_label_correspondence_wg.standalone.html` — 8 段：對應 Venn / 沒訊號拆解 / 為何 0.4 / Fisher paired-vs-tumor 兩軸。
- `20260618_subcluster_confirmation_spectrum_decision_01.standalone.html` — ≥2 群確認光譜 + 對齊軸 + 三態決策 + 完整切群總帳（**收尾 SoT**）。
- `20260618_within_label_subcluster_methodology_validation_01.standalone.html` — 合併方法學+三項驗證確認。
- `20260618_subcluster_workstation_chr1_01.standalone.html` — chr1 45 案逐項判讀工作站（嵌熱圖）。

## 鎖檔數據（注入 HTML 的真值，data-lock §13）
`spectrum_decision / split_accounting / contingency_summary / nosignal_breakdown / sensitivity(_ext) / section8 / case_validation / tumor_layers / tumor_vs_paired / fisher_cramersv.json`、`cases.json`、`spec.json`。

## 不入版控（.gitignore）
- `figs/`（52 PNG，已 base64 內嵌進 HTML）
- `records.json` / `records_wg.json` / `records_wg2.json`（中間數據，~18MB，可重生）

## 重現
```bash
# 1. 全 WG cluster×label 列聯表（逐 chr disk-safe，~31min，產 records_wg2.json）
python3 scripts/wg_contingency.py
# 2. 分類 / 拆解 / 門檻 / 對照（讀 records_wg2.json）
python3 scripts/classify_contingency.py 0.5 0.4 3
python3 scripts/classify_nosignal.py
python3 scripts/section8_data.py ; python3 scripts/tumor_vs_paired.py
# 3. 建 HTML（讀鎖檔 json）
python3 scripts/build_spectrum_html.py   # 等
```
> binary：`feat/summary-nreadsvalid @ 5c39051`；tumor=ClairS_pileup_v040、normal=5khz_simplex（分析端 filter 至 tumor）。
