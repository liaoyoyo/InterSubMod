---
title: 切割穩定性驗證 — 位點甲基切群「同一刀是否重現 / 夠好或更好」
date: 2026-06-22
status: in_progress
tier: 3
sample: HCC1395 tumor-only（單樣本 ⭐2-3 characterization）
scope: kprofile_examples 29 代表位點（splittable subset；全基因組 sweep 未跑）
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/kprofile_stability.json
build_commit: 578763c
observation_standard: true
---

# 切割穩定性驗證 — 同一刀是否重現、是否夠好/更好

> **核心問題**：位點附近甲基被切成 k 群後，(1) 這刀在擾動下能重現嗎？(2) 結果夠好或更好嗎？(3) 哪些位點圖能肉眼確認？
> **HTML（肉眼確認）**：`InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/20260622_clustering_stability_validation_01.standalone.html`
> **🔴 一句話**：raw clusterboot Jaccard **全 0.74–1.0**（幾乎全「看似很穩」）→ raw 穩定性**零鑑別力**；真正鑑別 = **excess-over-null + 同一 k coherence + Cochran e≥5 可靠性**。

## 0. 為什麼不能只看穩定性（先講陷阱）

直覺：bootstrap-Jaccard 高 = 切得好。**錯**。本資料已實證「可重現 ≠ 離散」：read-內甲基相關 / gradient 會**穩定地回傳同一刀**（A-path 相關感知 null 在真實資料失敗 81%；Şenbabaoğlu 2014 consensus 在 null 也假穩）。本輪實測直接證實：**29 個代表位點 raw Jaccard 全 0.74–1.0**，若只看 raw 全部「穩」。

## 1. 怎麼驗證 — 3 軸 + 雙 null（chance-corrected）

| 軸 | 工具 | null（chance 校正） |
|---|---|---|
| **① 穩定性** | clusterboot Jaccard（80% subsample × 300，Hennig） | **within-1-group**：逐 CpG 欄內非 NaN 重排 → 重算 BERNOULLI（40 draws）→ **excess = real − null(95pct)** |
| **② 對齊（非循環）** | CramérV(cut vs HP/carrier/allele 軸) | **label-shuffle**（500 次）+ **Cochran e≥5 可靠性閘** |

- **判準**：穩定 = `real>null95 AND excess≥0.10 AND jac≥0.6`；對齊 = `V≥0.3 AND shuffle_p<0.05 AND e≥5`。
- **GOOD = 同一 k 同時過**（coherent）；headline k = coherent 中 max-excess（無則整體 max-excess）。
- **雙 null 理由**（user 06-22 指定由資料特性選）：read-內相關是「假穩」根因 → within-1-group 擋假穩、label-shuffle 擋假對齊。
- **🟢 BERNOULLI 重算自檢**：重算距離 vs 既有 matrix.csv `max|Δ| = 1.1e-4`（僅來自存檔 3 位小數捨入）→ null 距離重算正確（`bernoulli_selfcheck_max_absdiff`）。

## 2. 結果：夠好或更好

### 2a. 「更好」（vs raw-only / 舊版）— 是
- **raw Jaccard 0.74–1.0 全高 → 零鑑別力**；**excess-over-null 跨 0.06–1.0** → excess 才有鑑別力（呼應 cross-sample ASM「看 excess 非 raw rate」）。
- 3 軸 gate **正確降級 13/29**（8 UNALIGNED + 5 DIFF_K）= raw-only 會誤判為「穩」者 → 嚴格優於只看穩定性。

### 2b. 「夠好」的邊界 — 單樣本天花板
| verdict | n | 意義 |
|---|---:|---|
| **GOOD_CUT** | 16 | 同一 k 穩定(excess≥0.10)+ 可靠對齊(e≥5) |
| **STABLE+ALIGNED_DIFF_K** | 5 | 細 k 穩定(exc 0.7–1.0,V=1.0)但該 k Cochran e<5 不可靠；可靠對齊的 k=2 又 excess 低 |
| **STABLE_BUT_UNALIGNED** | 8 | cut 可重現(raw 甚至 0.998)但不對齊任何軸 / 對齊不可靠（**陷阱**） |

**by group**（GOOD 集中高覆蓋）：
- multi-resolution 10/10 GOOD、single-k-forced 3/3 GOOD（n=350-370 高覆蓋撐得起細 k 對齊）
- confident-unique 3 GOOD / 4 DIFF_K / 1 UNALIGNED（「silhouette-confident」≠ coherent good）
- ambiguous-near-tie 7 UNALIGNED / 1 DIFF_K（穩但不對齊 = 陷阱主集中地）

🔴 **單樣本最關鍵發現**：「最穩的刀」（細 k、高 excess）與「最可靠對齊的刀」（k=2、e≥5）**常不是同一個 k**。要把 STABLE+ALIGNED_DIFF_K 升為「確認」需更多樣本/覆蓋（COLO829、cross-sample）。

## 3. 代表位點圖（6 個，肉眼確認；見 HTML）
| 位點 | verdict | 讀法 |
|---|---|---|
| chr11:61508413 | GOOD | k=4 real 0.998 / null 0 / V=1.0 e=5.3 → 穩+可靠對齊同 k |
| chr1:119365298 | GOOD（multi-res） | 多 k 皆過，穩定且對齊跨解析度 |
| chr15:94640645 | STABLE+ALIGNED_DIFF_K | k=4 exc 0.967 V=1.0 但 e=2.4<5；可靠對齊在別的 k = 天花板 |
| chr7:157702979 | STABLE_BUT_UNALIGNED | k=3 real 1.0 null 0（真結構!）但對齊 e<5 全不可靠 = 陷阱 |
| chr1:218115208 | STABLE_BUT_UNALIGNED | raw 0.998、V=1.0 但 e=1.1（over-fragment）= 假對齊 |
| chr4:190112507 | GOOD（single-k） | n=370 高覆蓋乾淨 2 群，exc 0.30 V=0.80 e=10.5 |

每張：左=穩定性(real vs within-1-group null) · 右=對齊(綠=可靠 e≥5；紅=不可靠) · 下=既有 dual-panel 實際結構。

## 4. 驗證表（headline 數字全 L1，從 kprofile_stability.json）
| 數字 | 值 | 來源 key | L |
|---|---|---|---|
| 代表位點數 | 29 | n_loci | L1 |
| BERNOULLI 自檢 max\|Δ\| | 1.1e-4 | bernoulli_selfcheck_max_absdiff | L1 |
| verdict 分布 | GOOD 16 / DIFF_K 5 / UNALIGNED 8 | verdict_counts | L1 |
| raw Jaccard 範圍 | 0.74–1.0 | loci[].per_k[].jac_real_mean | L1 |
| excess 範圍 | 0.06–1.0 | loci[].headline.stab_excess | L1 |
| multi-res GOOD | 10/10 | by_group | L1 |
| ambiguous UNALIGNED | 7/8 | by_group | L1 |

## 5. 限制
- **代表位點集**：取自 kprofile_examples（皆已 splittable）→ 無 state② 純無訊號位點；excess 下界偏高。**全基因組穩定性 sweep 未跑**（本輪只小規模代表）。
- **單樣本 HCC1395 ⭐2-3**：穩定/對齊皆 cis-ASM characterization，**「好切法」≠ subclone**（仍需 normal cis-control）。
- 門檻（excess≥0.10 / V≥0.3 / e≥5 / subsample 80% / null 40 / shuffle 500 / seed 20260622）為約定；per-k 全數字在 json 可 re-threshold。

## 6. 重生
```bash
cd InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot
python3 scripts/kprofile_stability.py          # → kprofile_stability.json（含自檢）
python3 scripts/plot_kprofile_stability.py     # → overview + 6 per-locus figs
python3 scripts/build_stability_observation.py # → standalone HTML（§13-A 注入）
```

## 7. 下一步（候選，未決）
- **全基因組穩定性 sweep**：對全 2215 (tumor) / 9194 (merged) k-choice 位點跑 chance-corrected stability gate，量化 GOOD/DIFF_K/UNALIGNED 全域比例。
- **跨樣本破天花板**：STABLE+ALIGNED_DIFF_K 在 COLO829 / 6 樣本是否升 GOOD（覆蓋疊加 → e≥5）。
- 關聯：[[project_distance_matrix_cluster_validation_methods]]、[[project_subcluster_cluster_count_determination]]、[[project_cross_sample_asm_reproducibility]]（excess-over-null）。
