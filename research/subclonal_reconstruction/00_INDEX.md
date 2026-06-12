<!--
建立時間: 2026-06-12
目標: Subclonal Reconstruction 新主軸的輸出索引（最終主輸出 + 每次比較測試輸出）— 後續 AI 找輸出的單一入口
狀態: skeleton（cycles 隨研究累積；final 隨論文圖表累積）
build_branch: research/subclonal-reconstruction-202606
關聯: docs/data_specs/20260612_output_organization_guide_subclonal_01.md（組織規則）
-->

# Subclonal Reconstruction — 輸出索引（00_INDEX）

> **後續 AI 找輸出先讀這裡。** 組織規則見 `InterSubMod/docs/data_specs/20260612_output_organization_guide_subclonal_01.md`。

## L0 放哪（決策表）

| 類型 | 放這裡 | 命名 |
|---|---|---|
| **最終主輸出**（論文 figure/table/result）| `final/{figures,tables,results}/` | `Fig{N}_{desc}` / `Tab{N}_{desc}` |
| **每次比較測試輸出**（per-cycle）| `cycles/{YYYYMMDD}_{topic}/` | 內 `00_README.md`+`data/`+`figures/`+`results.json` |
| **正典輸入**（ISM run，唯讀勿動）| `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/{sample}/paired_full/{date}_*_complete_matrix/` | — |
| **真值/可重跑腳本** | `docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/` | VERIFIED_RESULTS=SoT |

## L1 規則（避免混亂）
- 一 cycle 一資料夾，必含 `00_README.md`（目的/輸入/輸出/verdict/tier）+ `results.json`（真值）。
- 升上論文的圖表才進 `final/`；探索/比較留 `cycles/`。
- 每 cycle 完成**回填下方表一行**。
- 大中間檔（BAM/大 TSV）寫 `cycles/.../data/` 並 gitignore；數字落 `results.json`。

## Cycles（每次測試/比較輸出 — 完成回填）

| date | topic | verdict | tier | path |
|---|---|---|---|---|
| _(尚無；第一個預期 = G-A 跨 5 樣本 V10 重現)_ | | | | |

## Final（論文最終主輸出 — 升級後回填）

| Fig/Tab | 主張 | 資料來源 | path |
|---|---|---|---|
| _(尚無)_ | | | |

## 重用既有輸出（provenance）
- 7 樣本 ISM 正典輸入：見組織指南 §1 + §5。
- A0_assets V1-V12 真值：`VERIFIED_RESULTS.md`。
- normal BAM 契約：`docs/data_specs/20260612_external_data_dependencies_01.md`。
