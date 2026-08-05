<!--
建立時間: 2026-07-27 09:37
更新時間: 2026-07-27（全面驗證與 read 聚集圖完成）
目標: 提供多 sSNV pattern × 甲基關聯研究專案的人類入口
處理範圍: project navigation
關聯檔案:
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/00_INDEX.md
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/20260727_多sSNV_pattern與甲基關聯全面驗證_01.md
-->

# 多 sSNV pattern × 甲基關聯驗證

本專案以 7 個 technical datasets、chr1–22 的 frozen exact-PS strict
read-linkage資料，重新在 **exact raw HP tag** 內建立 RR/RA/AR/AA與較長
R/A/X read pattern，並測試 pattern-conditioned regional methylation
association。

**Status：validated。** 1,045 個 formal units 中找到 3 個較長 signature 的
robust regional methylation associations；formal full-four 與 exact
二位點 RA robust 都是 0。

入口：

- `00_INDEX.md`：pre-registration與方法契約。
- `20260727_多sSNV_pattern與甲基關聯全面驗證_01.md`：正式結果、三個位點、
  指定位點回查、限制與完整證據鏈。
- `pre-decision-audit.md`：GO/NO-GO邊界與紅隊結果。
- `implementation-notes.md`：設計決定、偏離、紅隊修正與驗證紀錄。
- `scripts/`：可重跑 producer/validator。
- 大型輸出：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/`
- 互動 HTML：
  `../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/20260727_multisite_pattern_methyl_association_06.html`

最新版 HTML 為 557 個有合法原始 matrix 的 formal units 加入
label-independent UPGMA dendrogram、對稱 read-distance heatmap 與
R/A／strand／RG strips；488 個不可畫的單元會顯示具體原因。

任何結果均限於 read-level association，不等同 subclone或 ancestry。
