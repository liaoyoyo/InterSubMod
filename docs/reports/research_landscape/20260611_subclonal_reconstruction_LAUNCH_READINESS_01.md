---
title: Subclonal Reconstruction 研究啟動就緒架構（單一 launch 入口）
date: 2026-06-11
type: research landscape / launch-index
status: launch-ready (HD-1 為唯一研究啟動 gate，用戶決定)
owner: liaotzuyu000@gmail.com
purpose: 把新主軸的「軸 / 數據 / 知識 / 流程 / 決策」一頁串起，供未來 session 乾淨啟動論文研究
---

# 🚀 研究啟動就緒架構 — Subclonal Reconstruction（未來任務從這裡開始）

> **本檔職責 = 單一 launch 入口**（thin navigation + 就緒驗證 + 啟動決策樹）。內容不重造，指向各權威 SoT。
> **2026-06-11 整理+整合完成後產出**：所有散落工作已收斂到 `develop`(=活躍線 54b02ab)；本檔是「知識+數據+架構整理清楚、可開始研究」的確認頁。

---

## L0 一眼結論

**主軸**：Subclonal reconstruction using somatic haplotagging + methylation profiles, Nanopore（取代 G6；G6/G1 降支撐）。
**就緒度**：✅ 數據齊（6 樣本×3 癌種）· ✅ 知識整理（V1-V12 + 4 道 NEGATIVE + 紅線 + gaps）· ✅ 流程（master-workflow + cycle skills）· ✅ git 整合（develop=canonical）。**唯一啟動 gate = HD-1（你決定）**。

---

## L1 數據資產（已 spot-check ✓）

6 cell line × 3 癌種，三要素齊（tagged BAM + somatic VCF + 5mCG）：HCC1395/1937/1954（乳腺）· H1437/H2009（肺）· COLO829（黑色素瘤）+ HCC1395_DORADO（basecaller 對照）。
- 路徑型式（**完整路徑，勿縮寫**；2026-06-12 ls 驗證）：`/big7_disk/liaoyoyo2001/big7_disk_output/canonical/{sample}/paired_full/{date}_*_complete_matrix/longphase_s/{sample}_tagged.bam` (+ 同目錄 `somatic_pass.vcf.gz`)。⚠ 根 = `big7_disk_output/`（非 repo 內 `InterSubMod/canonical/`）。{date}：HCC1395=20260314、其餘=20260315。ISM 輸出：`intersubmod_tp/filtered_snv_tp/` + `intersubmod_fp/`。
- normal 甲基（V10 用）：`/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam`。
- **意涵**：單樣本 ⭐3 上限可被打破 → 跨 6 樣本可達 ⭐4（但 single-pipeline longphase-S 自我參照 = tier 風險，需正交對照 G-E）。
- 真值 SoT（V1-V12 全數字）：`InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md`

## L2 知識地圖（過去結論整理 → 論文元件）

| 論文元件 | 結論 | 級 | SoT |
|---|---|---|---|
| **脊柱 §R1** LOH-constrained 子克隆重建 | POSITIVE（Grade B+，需 HD-1 caveat）| ⭐3 | 整合篇章 §2-3 |
| **§Methods-Neg** 甲基不能 filter/判別（死四道）| 最 reviewer-proof，今天防彈 | — | 整合篇章 §3 NEGATIVE |
| **§R2** 子克隆內 ASM characterization | 存在·小·非方向·跨 3 癌種現象復現 | ⭐3 | foundation §3 + ASM 工作站 |
| **甲基貢獻** | characterization + 誠實負目錄（**非重建驅動**）| L2 | 兩 foundation reconcile |
| 🔴 **必守紅線** | 甲基=germline-haplotype 層級·非 filter(DEAD)·T3 可用性 NEGATIVE·T2 勿宣稱歸 H3·BRCA2 非乾淨 cis(改 chr17/TBC1D16)·cohesion≠cis | — | foundation §4 + 整合篇章 §5 |
| **開放 gap** | G-A 跨 6 樣本(tier 解鎖)·**G-B within-hap somatic-vs-baseline null(gate 甲基貢獻)**·G-C cis 足跡·G-D 重建 demo·G-E 正交 pipeline | — | foundation §5 |

**權威 SoT 兩份**（互補面）：
- ① 甲基-assist 面 + V1-V12 + 資產：`InterSubMod/docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md`
- ② ASM-characterization + 四道 NEGATIVE + 脊柱 + HD-1：`InterSubMod/docs/concepts/2026/06/20260611_Subclonal_Reconstruction_Paper_Focus_整合篇章_01.md`
- KB 結論層：`knowledge/09_conclusions/` + `knowledge/11_external_literature/10_paper_readiness_convergence.md`（7-thread 稽核）
- ❌ **DEAD 勿再開**：甲基→TP/FP filter（死四道）· T3 subclone 可用性（救 ambiguous read 偏向反了）

## L3 流程（怎麼執行）

- **主工作流地圖**：`InterSubMod/docs/references/20260611_master_workflow_architecture_01.md`（agentic-loop 四層 + git 四決策表 + 鉤子）。
- **研究 cycle**：P0 `/cycle-init` → P1 `/research-loop`(plan.json OSF 預註冊+stop_criteria) → P3 `/feature-layered-observation` → P4 `/multi-sample-consistency` → P5 `/run-evaluator` → P6 `/conclude-research`。
- **論文撰寫**：`/narrative-frame` 或 `/structured-tech-report`（title/abstract/figure 骨架）；對外宣稱必繼承 `/scientific-rigor` §2-§7（evidence tier + Cohen + DAG + pre-reg）。
- **並行 session**：各開 `bash scripts/new_session_worktree.sh <topic>`（勿共用主 dir+HEAD）。
- **數字誠信**：§13（產數字→落檔→Read→才寫報告；number_provenance gate 已 live）。

## L4 啟動決策樹（唯一 gate = HD-1）

```
HD-1（你 GO/NO-GO，AI 不代決）：phasing 重建脊柱 by-construction 部分循環
├─ 選 1：跑 R-SELFREF 全 7 樣本 flag-on（~25-50hr C++ 長計算，主回合背景 Bash 落檔）→ 脊柱寫成 positive 主結果
└─ 選 2：phasing 降 characterization observation → 四道 NEGATIVE 扛論文（methods/negative 定位，最快防彈）

撰寫關鍵路徑（整合篇章 §7）：
① 定 HD-1 → ② 先寫四道 NEGATIVE methods（今天防彈，HD-1 獨立）→ ③ 依 HD-1 grade 寫脊柱
→ ④ ASM characterization + chr17 exemplar + 跨樣本現象 → ⑤ 改 stale 文件（BRCA2 已 amend）
平行可先做：G-A 跨 6 樣本重現（tier 解鎖）、G-B null 對照（gate 甲基貢獻）
```

## 待清項（非啟動阻擋；研究啟動後或投稿前處理）

- 🟡 KB 5 個 snapshot 仍 G6-era 過時（`knowledge/0{1,2}_*.md` 等 2026-05-18）→ 用 `knowledge/scripts/refresh_last_verified.py` 或人工更新到本軸（freshness hook 持續警告中）。
- 🟡 `main` 仍舊 unrelated 快照（develop 已=活躍線）；redundant branch cis-asm/g6-paper-focus 可 archive。
- 🟡 shared-state（ledger/queue/active.json/CURRENT_FOCUS/postmortem logs）+ pre-existing 他 session 改 + settings hook wiring：未 commit，依 §C 留用戶。
- 🟡 deferred：cross_sample_benchmark fan-out cap + per-sample checkpoint（跑 benchmark 前套）。

> **使用方式**：未來啟動論文研究 → 先讀本檔 L0-L4 → 定 HD-1 → 依 L4 關鍵路徑用 L3 流程執行 → 數據 L1、知識/紅線 L2 隨時對照。
