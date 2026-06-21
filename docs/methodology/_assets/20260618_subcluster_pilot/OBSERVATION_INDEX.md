# 觀察索引 — 位點甲基切群 / per-CpG 歸屬 / ISM vs modkit（2026-06-20）

> **這份是什麼**：把「位點甲基能否切群 → per-CpG 定位 → ISM vs modkit 對照」這串研究的**所有可肉眼觀察 HTML + 圖 + 資料 + doc** 集中索引，之後要觀察直接從這裡找。
> **HTML 開法**：主 repo 同層的 `*.standalone.html` 自含（base64 圖），瀏覽器直接開。**圖在 `figs_*/`**。**資料 JSON 在 worktree（gitignore 可重生，見各 doc §重生）**。

## 觀察 HTML（主 repo `docs/methodology/_assets/20260618_subcluster_pilot/`）

| HTML | 看什麼 | 對應 doc | 資料源 |
|---|---|---|---|
| `20260620_decisionflow_5state_classification_01.standalone.html` | 判別流程 5 態組成 + 門檻依據 + 切不出 mean-shift（tumor/merged）| decisionflow doc | decisionflow_summary.json |
| `20260620_brca2_percpg_observation_01.standalone.html` | BRCA2 per-CpG 多軸歸屬：定位 track / 三軸 / dual-panel / 多軸表 / 共線（逐項判讀）| brca2 pilot doc | brca2_{idea1,multiaxis,percpg_table}.json |
| `20260620_ism_vs_modkit_percpg_comparison_01.standalone.html` | ISM joint vs modkit-equiv per-CpG 全資料 Venn（both/ISM-only/modkit-only）+ by-state | ism_vs_modkit doc | percpg_compare_summary.json |
| `20260620_s3_percpg0_joint_observation_01.standalone.html` | 7 案「modkit per-CpG=0 但 ISM joint 顯著」dual-panel（location-clean vs dispersion，逐項判讀）| ism_vs_modkit doc §3b | s3_index.json |

## 觀察圖資料夾（`figs_*/`）

- `figs_decisionflow/` — 5 態 4 圖（組成 / minority 分佈 / 門檻敏感度 / mean-shift 分解）
- `figs_brca2/` — BRCA2：track_carrier_dbeta / track_3axis / dualpanel_brca2 / dualpanel_brca2_tumor
- `figs_cmp/` — cmp_venn_state（Venn + by-state）
- `figs_s3/` — 7 案 dual-panel（chr*_pos.png）

## doc（主 repo `docs/experiments/in_progress/2026/06/`）

- `20260620_decisionflow_5state_classification_wg_01.md` — 判別流程 5 態（N=30,490，tumor/merged）+ 驗證表
- `20260620_per_cpg_multiaxis_attribution_brca2_pilot_01.md` — BRCA2 per-CpG 多軸 + ISM 輸出方向
- `20260620_ism_vs_modkit_percpg_comparison_wg_01.md` — ISM vs modkit-equiv 全資料對照（含 §0 modkit 能力 / §3b S3 dual-panel 判讀）

## 🔗 姊妹工作（真 modkit/DSS binary，另一線）

- `docs/experiments/in_progress/2026/06/20260620_subclone_label_methods_pilot/` — 三法（想法1 per-CpG / **想法2 真 modkit/DSS** / 想法3 per-label Δβ）BRCA2 pilot；INDEX 線 ≈「三個甲基↔標籤關聯方法」。本線用 **modkit-equiv（Fisher）**，姊妹線有 **真 modkit binary** head-to-head（22305 Jaccard 0.76 / 07 區域 modkit −0.159 vs ISM −0.122）。兩線互補。
- 推論 SoT：`docs/plans/20260620_subclone_situation_verification_reasoning_01.md`

## 一句話結論（已驗證）

modkit/DSS 式逐-CpG 定位**夠用**（per-CpG 非 ISM 創新）；**ISM 差異化 = read-level 結構偵測**（diffuse location + heterogeneity/dispersion + 無監督群發現，三者 modkit 皆無對應）。訊號弱（F 2.4-7.2≪BRCA2 15.4）、單樣本 ⭐2-3、結構≠subclone。memory：`project_ism_vs_external_methylation_tools_comparison`、`project_decisionflow_5state_classification_wg`。
