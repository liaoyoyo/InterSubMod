<!--
建立時間: 2026-07-06
類型: Handoff SoT (D — 交付/跨 session 理解)
主軸: Subclonal reconstruction — 甲基/CN 線本 session 統整
build_branch: research/subclonal-reconstruction-202606
關聯: docs/methodology/20260628_subclone_reconstruction_master_spec_01.md (數字 SoT), memory MEMORY.md
-->

# 甲基/CN 線 Handoff SoT（2026-07-06）— 新 session 先讀

> **一句話**：本 session 補了 **SAVANA CN（4/7 可靠）** 並用它做**甲基 vs sSNV-lineage 因果驗證** → 結論**甲基群內雙峰真實且 CN-獨立（非 CN-gain 假象）但屬 within-lineage 隱藏子結構 L3、bulk 不能證 subclone** → **用 CN 定量佐證「甲基=有界輔助」主軸**。數字 SoT 仍是 master_spec（並行 session 正重算，見 §5）。

## §1 主軸（穩定，不變）

**論文**：Subclonal reconstruction using somatic haplotagging and methylation profiles（中文碩論）
- 🔴 遺傳 sSNV 單分子共現 = **唯一非循環重建骨幹**
- HP tag = allelic/clonal **鑑別器**（非確認器）
- 🔴 甲基 = **bounded-auxiliary**（非 subclone 偵測器）
- 單樣本 **⭐3** 封頂；多區域整合「**定不出來即答案**」（field-endorsed）
- 真天花板 = **single-cell / 正交 ground truth**（非只等資料；⭐3→⭐4）

## §2 本 session 甲基/CN 線（完成，L1 CN + L3 甲基）

### 2a. SAVANA CN（補 cn=unknown）— commit 8aedc86 / 97240f1
- **cna-only**（跳過 SV `run` 的 488G OOM 殺手）；raw big8 tumor BAM（tagged 的自訂 HP tag 撞 SAVANA）；`-v` matched-normal het。
- **4/7 可靠 CN**：HCC1395(SEQC2 既有) + **H1437 0.95 / H2009 0.95 / HCC1954 0.66(自洽)**。
- **3/7 cn=unknown**：COLO829(0.53 mis-fit,BAF~1.0 矛盾)、HCC1937(0.62 mis-fit)、DORADO(共用 HCC1395)。
- 工具：`scripts/analysis/savana_to_smcnbed.py`(CN→SM_CNBED) + `rerun_cn_integration_from_savana.sh`。
- 記憶 `project_savana_cna_only_pipeline_and_colo829_purity`。

### 2b. 甲基 vs sSNV-lineage 因果驗證 — commit 2689087 / 1a23d74
- 用並行 session 既有 `methyl_auxiliary_annotation.py`（每 sSNV-genotype 群測甲基雙峰,GMM,扣 HP/CN）+ 我的 SAVANA CN。
- 🔴 **因果證明（雙峰率 by CN,`methyl_bimodality_cn_rate.py` 8-chunk）**：**率跨 CN 平坦** — H2009 gain **0.170** ≈ neutral **0.163**（**1.04x**,最有 power）;HCC1954 1.25x;H1437 2.41x(neutral n=49 太小不可靠)。→ **CN-gain 未提高群內甲基雙峰率**。
- 🔴 **自我修正**：先前「77-86% residual 在 gain = CN 假象」是 **base-rate 假象**（gain 區佔基因組大→分母效應）。**教訓：count 富集≠因果,必看 rate(帶分母)。**
- **裁決**：群內甲基雙峰 ~15% = **CN-獨立的真實現象**（甲基有超越 CN 的獨立內容）;**但 within-lineage 隱藏子結構 L3**(cis-ASM/HP-殘留/技術未排除),**bulk 不能證 subclone**。
- 記憶 `project_methyl_lineage_cn_gain_confound`。

### 2c. Q5 難建樹區甲基輔助旗標（用戶想法）— commit 1a23d74
- `C_underdetermined`(sSNV 定不出樹) 區 + ALT-cluster + cn-clean(neutral) + 扣 HP 甲基雙峰 = **7 候選**(H2009 6+HCC1954 1,Δβ 0.2-0.3)= L3 難建樹區甲基輔助（候選 subclone/branching）。**稀少 + 需 normal-cis/single-cell 確認。**

## §3 觀察 HTML 位置（verify-workstation 格式）

- `InterSubMod/docs/experiments/in_progress/2026/07/_assets/methyl_lineage/methyl_lineage_H2009_workstation.html`（H2009 pilot 10 樹區,備查）
- `InterSubMod/docs/experiments/in_progress/2026/07/_assets/methyl_lineage/q5_methyl_hardtree_workstation.html`（Q5 7 候選,逐項判讀）

## §4 遺傳骨幹線（並行 session,活躍;數字在流動）

- incompatible 重分類 118→18（隱藏祖先樹）;gap#1 partial-subcube 救回全 7 樣本 +15,979 區;basecaller 重現性(5khz vs DORADO branched% 29 vs 30);7 樣本 per-sample 拓撲工作站;三軸報告重算到 07-05/06。
- 🔴 **數字 SoT = `docs/methodology/20260628_subclone_reconstruction_master_spec_01.md`**（並行正重算,舊 35,332 sSNV/7,143 區部分過時,如 HCC1395 3885→6288）→ **統一數字待重算穩定。**

## §5 之後方向 / 目標 / 任務

| 優先 | 任務 | 依賴 |
|---|---|---|
| 🔴 真天花板 | single-cell / 正交 ground truth（⭐3→⭐4）| 缺資料 |
| ⭐4 blocker | COLO829 normal 甲基（bounded-auxiliary 負向重現）| 缺 |
| 本線延伸 | 甲基-lineage cn-clean + Q5 7 候選 → **normal-cis control**（分辨 cis-ASM vs subclone）| tumor+normal BAM 已在 |
| 收尾 | 統一數字：待並行重算穩定 → 更新 master_spec + handoff | 並行 session |
| 治理 | 主 branch ahead origin ~194 未 push | 需用戶確認 |

## §6 commit 索引（本 session 甲基/CN 線）
8aedc86(轉檔器+driver) → 97240f1(SAVANA CN 結果 doc) → 2689087(甲基-lineage CN cnclean) → 1a23d74(CN-aware 率修正+Q5)。docs：`docs/experiments/in_progress/2026/07/20260705_savana_cna_6sample_cn_results_01.md`、`20260705_methyl_clustering_vs_cn_crosssample_01.md`、`20260706_methyl_vs_ssnv_lineage_01.md`。
