# self_phasing_synthesis_PI/figures/ — 復原狀態（2026-05-11）

此目錄因 2026-05-11 git history rewrite 操作意外影響，原 13 PNG 全部消失。
以下為復原狀態：

## ✅ 已復原（從 v3/figures/ cp 過來，3 個）
- igv/D_SP1_chr19_17565944.png
- igv/D_SP2_chr19_12452332.png
- igv/D_SP3_chr19_12467180.png

## 🔄 需用 image-gen skill 重生（4 個，prompts/G*.txt 還在）
- G1_player_as_referee.png — 用 `prompts/G1_player_as_referee.txt`
- G2_pass_two_flow.png — 用 `prompts/G2_pass_two_flow.txt`
- G3_getVote_three_layer.png — 用 `prompts/G3_getVote_three_layer.txt`
- G4_germline_absent_three_versions.png — 用 `prompts/G4_germline_absent_three_versions.txt`

執行命令：呼叫 `/image-gen` skill，指向 prompts/ 資料夾。

## ⚠️ 永久遺失（6 個 master 數據圖，需手動重產或從 PI 拿）
- master/F1_priority_bug_mechanism.png
- master/F2_priority_bug_per_chr_enrichment.png
- master/F3_binary_commit_timeline.png
- master/F4_chr19_752_victims_scatter.png
- master/F5_layer15_zero_sum_4quadrant.png
- master/F6_paired_vs_TO_HP_distribution.png

來源候選：
- research/paired_priority_bug_audit/ 下的 .py scripts（需找對應 generation script）
- PI Drive / Google Drive 的副本
- 從 phaseB/C/D run 數據重產（dataset 還在 .gitignore 路徑下）

## 引用方
`build_html.py` + `build_pptx.py` 引用這些 PNG 於 5/11 HTML preview slides。
