# 案例模板 — case templates（依目標組合小物件）

> case = 一個「解釋目標」對應的物件組合 + 佈局。**先問「要讓人懂什麼」→ 再選物件**。
> case = `{primitives 列, 各自參數綁定, data_source 路徑, synthetic_flag}`。

## goal → objects 對照（怎麼挑物件）

| 解釋目標（要讓人懂的） | 選哪些物件 | case |
|---|---|---|
| 「一個 read-level 方法怎麼把分子收斂成一個方向性判斷」 | P1 read×CpG + P2 分軌 + P3 β bar + P4 Δβ step + P5 normal-cis | **C1** |
| 「為何舊方法抓不到、新方法抓到」（互補非取代） | 左 P6 舊法 V box ‖ 右 P3+P4+P5 新法，共用座標/色標 | **C2** |
| 「這方法的系統性盲點 + 怎麼補救」 | P7 交叉型（synthetic）+ 對照 max|Δ|/ARI | **C3** |
| 「一個統計量的逐步定義」（純流程） | P3→P4→P5 串接，省略 P1 示意 | Cx_custom |

挑物件三原則：**一圖一概念**（一個 case 一個敘事單元，別塞）/ **先直覺後公式**（每物件前一句白話）/ **示意與真值分軌**（P1/P7 示意，P3/P4/P5/P6 真值）。

---

## C1 — 方向型 anchor 對照（directional anchor）
- **目標**：用單一真實 anchor 位點，展示「同 germline 單倍型依 somatic allele 分兩群 → 多數 CpG 一致地差 → Δβ 抓到 + normal 基線拆 cis/drift + copy-partition 誠實標真 cis 小」。
- **物件**：P1（read×CpG，HP1+HP1-1 兩組）→ P3（β bar 0.228 vs 0.108）→ P4（Δβ=−0.122 主鍵 + max|Δ| 副鍵 + Wilcoxon）→ P5（normal-cis 三條 + d_cis/d_drift/d_within）。
- **真值源**：`control_cohesion_cistest.json`（β/Δβ/d_cis/d_drift/n_shared）+ `brca2_pertag_cohesion.json`（silhouette）+ `copy_partition_confirm.json`（d_within/d_copy/perm p）。
- **誠實標註**：單樣本 ★3；真 cis 小（d_within=−0.023 marginal）；label-flip（HP1-1≡HP2-1）。
- **pilot**：`examples/c1_brca2_dbeta/`（已驗證可跑）。

## C2 — 新舊方法並列（old vs new）— *待擴 P6*
- **目標**：同一 tag、同一 locus，差別在**指標設計**：舊 ISM Cramér's V（無監督分群類別關聯，混 germline-allelic+copy、tumor-only 無 normal、無方向 → BRCA2 淹沒在數千高-V 中）vs 新 Δβ（方向性 + somatic-controlled + normal-anchored + copy-partitioned + 可排序）。
- **物件**：左 P6（cluster↔hp 列聯 + V hp=0.80/V alt=0.51 + 56 reads/229 CpG/gating）‖ 右 P3+P4+P5（新法 3 步），中間箭頭，共用座標/色標。
- **互補口徑（重要）**：呈現為**互補非取代**——舊法給「**存在性 YES**」、新法給「**cis 判別 + 可排序**」。
- **真值源**：⚠ 舊 ISM V/56/229/gating 在離線 `bip8_output_archive significance.json`（**轉 validated 前須落檔到可 grep 源**，否則 number_provenance_check.sh 對 validated 路徑 exit 2）。
- **擴充**：實作 `p_cluster_v_box()`。

## C3 — 盲點合成模擬（blindspot, synthetic）— *待擴 P7*
- **目標**：誠實展示 Δβ-mean 主鍵的系統性盲點——前後反向（交叉型）pattern → Δβ mean≈0/Wilcoxon NS（主鍵漏）；但 max|Δ| 大/盲分群 ARI 抓得到（副鍵救）→ 雙鍵動機（工程任務 #19）。
- **物件**：P7（交叉型 read pattern，**全 synthetic**）+ 對照 max|Δ|/ARI 數字（合成，標 synthetic，不進 verified）。
- **鐵則**：**只准合成不准混真值**；獨立 namespace + watermark；標題/figcap/footer 三處明標「合成模擬，非真實位點」。
- **擴充**：實作 `p_blindspot_scenario()`。

---

## 如何加一個新 case（通用化）

1. 先寫「解釋目標」一句話（要讓人懂什麼）。
2. 查上表 goal→objects 挑物件（或定義新 primitive，見 `object_library.md`）。
3. 在 `figure_spec.schema.json` 的 `case_template_id` enum 加值（或用 `Cx_custom`）。
4. 建 `examples/<case>/data.json`（verified 附 src + schematic）+ `figure_spec.json`。
5. 跑 renderer → 對 `explainer_checklist.md` 驗證 → 迭代圖例細節。
