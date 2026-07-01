---
title: 多樣本拓樸工作站 — 產生程式與 HTML 檔案地圖 / 建置文件
date: 2026-07-01
status: reference
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_topology_workstation.py, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json
---

# 多樣本拓樸工作站 — 建置地圖

## 🎯 你要的兩個路徑（repo 相對 + 絕對）

| 檔 | repo 相對路徑（用 InterSubMod/ 前綴開） | 絕對路徑（瀏覽器 file:// / 執行用） |
|---|---|---|
| **產生 HTML 的程式** | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_topology_workstation.py` | `/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_topology_workstation.py` |
| **產出的 HTML（主結果·33MB）** | `InterSubMod/docs/methodology/_assets/20260629_multisample_topology_workstation.standalone.html` | `/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260629_multisample_topology_workstation.standalone.html` |

> HTML 是 **standalone**（所有資料 inline 嵌入，單檔可直接用瀏覽器開，無需 server）。33MB 可 git 追蹤但**可從 build script 重生**。

## 如何重建 HTML

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts
python3 build_topology_workstation.py
# → 寫出 ../../20260629_multisample_topology_workstation.standalone.html
```

**環境變數（可選，皆有預設）**：
- `SM_OUT=<path>`：改輸出 HTML 路徑（預設 `20260629_multisample_topology_workstation.standalone.html`）。
- `SM_SAMPLES="name:dir,name:dir"`：指定樣本與資料目錄（預設自動納入 HCC1395 凍結 + `multisample_subclone/` 下有 topology 的樣本）。

## 資料流（上游 pipeline → JSON → build script → HTML）

```
BAM(ONT tumor+normal)
  │  sm_linkage_genomewide.py (O2 region-pileup, commit 37ddd73)  ← 單分子共現骨幹
  ▼
sm_region_integration.py (C2 四配子計數/incompatible, 44ed19c)
  ▼
topology_analysis.py (C3 去噪 eps floor + determinacy, 979c54d)  → topology_per_region.json
  ├→ candidate_scoring.py (D4 移除甲基排序 overclaim)            → candidate_scoring.json
  ├→ backbone_resolution.py (parsimony 機率回填)
  └→ enumerate_candidate_trees.py (B1 誠實候選樹)               → candidate_trees.json
        ＋ region_gene_annotation.json / single_snv_accounting.json / rd_perregion.json / ideogram_data.json
        ＋ longphase-TO tumor_phased_LOH.bed (LOH 底帶)
  ▼
build_topology_workstation.py  ← 讀上述 JSON，前端純手刻 SVG/JS，數字 inline 注入(§13-A)
  ▼
20260629_multisample_topology_workstation.standalone.html  ← 主結果(7 分頁樣本)
```

## 每樣本輸入檔（`_load_sample()` 讀取；build script L34-67）

| 檔名（每樣本目錄下） | 用途 | 缺檔行為 |
|---|---|---|
| `topology_per_region.json` | **核心**：每區 genotype 向量→拓樸/determinacy/HP/node_hp/pos_vaf/pos_nested/cn/span | 必需（無則該樣本不載入） |
| `candidate_scoring.json` | 評分佇列 + situation_dist（🔴 前端 regSit 須逐字複製其 situation()） | 空佇列 |
| `region_gene_annotation.json` | 基因註釋（GENCODE+DGIdb+COSMIC 癌症基因；癌基因命中篩選來源） | 空 |
| `single_snv_accounting.json` | 全 sSNV 宇宙帳本（linked/underpowered/isolated 三桶） | 宇宙帳本卡留白 |
| `candidate_trees.json` | 上游 enumerate 的等機率候選樹（carousel） | 前端不顯示 carousel |
| `rd_perregion.json` | read-driven 交叉確認（22-way） | detail 顯 — |
| `ideogram_data.json` | per-sample HG38 ideogram（染色體長度/樹位點/密度軌；census 衍生免 BAM） | ideogram 留白 |
| `{S}_clp_cn.bed` | （fallback）CLP CN caller 的 loh 段 | 無則走 longphase-TO |

**目錄位置**：
- **HCC1395（凍結主樣本）**：`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/`（+ ideogram 從 `.../20260630_perregion_workstation/data/ideogram_data.json` 補；+ `chr17_tree_data.json` 供固定教學圖鑑）。
- **6 樣本**：`/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone/{COLO829,H1437,H2009,HCC1937,HCC1954,HCC1395_DORADO}/`。

**LOH 底帶（全 7 樣本用 longphase-TO `tumor_phased_LOH.bed`，SEQC2-validated；build script L79 `LOH_BEDS`）**：
- HCC1395 / HCC1395_DORADO：`/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_phased_LOH.bed`（1094 段；DORADO 同細胞株借用）
- COLO829：`.../big7_disk_output/synthesis/research_rounds/20260423_colo829_to_pilot/step03_longphase_to/tumor_phased_LOH.bed`（191）
- H1437 / H2009 / HCC1937 / HCC1954：`/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/{S}/tumor_phased_LOH.bed`（34/134/538/18）

## HTML 功能總覽（使用者看得到什麼）

**頂部整體觀察區**（隨分頁樣本變）：① 每樣本記分卡（已確定/pairwise/單群欠定/跨HP/順序/衝突/截斷/LOH/癌基因）② 全 sSNV 宇宙帳本 ③ **HG38 全基因組分布 ideogram（可摺疊）** ④ 4 張統計卡（拓樸型態/群數 c/determinacy/HP 根數，點看放大+解釋）。

**HG38 ideogram — 4 個 toggle 模式**（皆疊 LOH 紫背景帶 + CN gain紅/loss藍 細帶）：
1. `結果+LOH`：每區 tick 依 situation 上色（綠已確定/藍pairwise/灰單群欠定/橙跨HP/黃順序/紅衝突）。
2. `結果/每sSNV`：每個已定相 sSNV 一點·同區連底線=同一棵樹。
3. `HP單倍型`：每區依 germline HP 家族上色（藍H1/紫H2/橙H1H2跨HP/灰未定相）。
4. `樹形`：樹位點依形狀（綠full_tree/藍結構/紅成環）+ 密度軌（橙underpowered/灰isolated）。

**逐區檢視**：左側 region list（篩選：chr/排序/最少群數/TP-FP/僅LOH區/僅無法定義/**僅癌基因命中**/搜尋；tag：拓樸/ctx/🟣LOH/🧬癌）→ 點區 → 右側 detail（locus 突變排列 lollipop〔VAF Y軸+距離比例尺〕 / 克隆樹或跨HP兩棵HP樹 / 全位點分群樹〔截斷區〕 / 可能固定樹 carousel〔成環〕 / 替代整樹候選 / 甲基裁決 / resolution_path / read-driven / 細胞群表 / 基因註釋）。

**底部**：候選評分確認佇列（左右判讀 + localStorage + 匯出 JSON）。

## 資料誠信 caveat（用前必讀）
- 🔴 **前端 `regSit()` 必逐字複製 `candidate_scoring.py` 的 `situation()`**，否則 ideogram/記分卡/佇列計數漂移。
- **CN gain/loss 只有 HCC1395**（cn 欄完整 2331 gain/65 loss）；6 樣本 cn=unknown → CN 帶留空（只有 LOH bed）。
- **determinacy 為 C2/C3 修正後值**（A_determined 1741 / incompatible 118；非修前 1812/12）。
- **complete-case**：HCC1395 有 40.4%(2887/7143) n_sSNV≥2 區因未全覆蓋被丟；細節見演算法稽核 memo。
- 建構有效性依賴 infinite-sites，癌症 LOH/CNV 系統違反。

## 變更歷史（本輪，2026-07-01；`git log` 於 build_branch）
`89d0a0b` 10 項整合 → `9ded20b` LOH全7樣本+節點美化+per-sSNV → `4b71166` locus Y軸/比例尺 → `9a7cc81` LOH標籤+正名+透明度+稽核memo → `5fbb819` 刷新C2/C3/D4資料+癌基因篩選 → `d571ef9` ideogram 可摺疊 → `f3d7bc5` HP模式+CN帶 → `d99ae7d` locus 未定相(H?)不顯示 bug 修正。

## 相關文件（SoT）
- 演算法稽核（Q1-Q4 + 上游待辦 ②④）：`InterSubMod/docs/methodology/20260701_topology_algorithm_audit_findings_01.md`
- sSNV 骨幹方法規格 + C2/C3/D4：`InterSubMod/docs/methodology/20260701_ssnv_backbone_method_spec_and_correctness_audit_01.md`
- 主 spec：`InterSubMod/docs/methodology/20260628_subclone_reconstruction_master_spec_01.md`
- 記憶索引：memory `project_topology_workstation_features_state`
