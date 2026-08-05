<!--
建立時間：2026-07-25
目標：在建立七資料集跨樣本比對矩陣與 HTML 前，界定可驗證假設、反例與主張上限
處理範圍：chr1–22；7 technical datasets；6 biological IDs；21 technical pairs
關聯檔案：research/20260725_crosssample_topology_matrix_report/implementation-notes.md
-->

# 跨樣本 sSNV／區域拓撲矩陣 pre-decision audit

## §0 Cynefin domain

**Complex → probe-first。** 相似度會隨比較層（精確位點、VAF 邊際分布、區域邊、粗拓撲形狀）改變；不能用單一最佳實務直接宣稱 clone tree 相同。本輪將各層平行量化，再依 denominator 與 biological-ID 去重結果收斂。

## §1 Observation completeness

| Observation | 狀態 | Tier | 來源 |
|---|---:|---:|---|
| 7 technical datasets × chr1–22 strict site/edge sets 已有 21-pair receipt | ✓ | L1 | `research/20260723_hcc1395_crosssource_topology_resolution/strict_pair_validation/results/strict_pairwise_metrics_01.json` |
| 7 datasets exact-PS topology census 為 71,955 ranked units 且 receipt checks pass | ✓ | L1 | `output/synthesis/observation_workspaces/20260724_exact_ps_cpp_topology_signature_census/all7_v1/receipt.v2.json` |
| HCC1395／DORADO coordinate-matched topology comparator 已有 machine checks | ✓ | L1 | `output/synthesis/observation_workspaces/20260725_exact_ps_cohort_similarity/all7_v1/cohort_similarity.v1.json` |
| 同 biological ID positive pair | ✓，但僅 1 組 | L2 | HCC1395 × HCC1395_DORADO |
| 同癌種不同 ID | ✓，僅 4 個去重 biological-ID pairs | L2 | breast 3 pairs + lung 1 pair |
| 外部 clone truth 可確認每條 parent→child edge | ✗ | L5 | 本輪無單細胞／多時間點 clone truth |

## §2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20/20 | 同 ID 應共享精確變異背景；同癌種可能只共享較粗分布 |
| 觀察支撐 | 10/20 | 全 7 datasets，但 same-ID positive 僅 1 組 |
| 機制清晰度 | 10/20 | phase、region、partial-read、edge identifiability 合理但尚無完整因果 ablation |
| 反例風險 | 10/20 | topology-shape common-mode 與 technical duplicate pseudoreplication 皆可能 |
| 所需資源 | 20/20 | 主要輸入已凍結，可在本輪完成矩陣與 HTML |
| **TOTAL** | **70/100** | **GO：建立分層比較報告；不得升格為全域 clone tree validation** |

**Falsifier observable：**

1. 若方法只產生共同形狀偏好，HCC same-ID pair 在精確位點／邊不會顯著領先不同 ID。
2. 若粗拓撲具強 biological-ID specificity，HCC pair 應在 chromosome-block bootstrap 穩定排名第一且與次佳差的 95% interval 全為正。
3. 若癌種造成可重現結構相似，同癌種不同 ID 在 biological-ID 去重後應持續高於跨癌種，而非只由 HCC technical duplicate 重複計數造成。

## §3 Assumption map

| 重要性 | 已知 | 未知／本輪先驗證 |
|---|---|---|
| 高 | sample/cancer mapping；chr1–22 scope；strict set metric definitions | 同癌種效應是否跨 metric；粗拓撲領先是否 bootstrap-stable |
| 低 | 視覺排序與色盤 | 不同 profile 權重的替代選擇 |

## §4 Probe 與 checkpoint

1. 讀 21-pair strict set metrics → 驗證：pair count=21、HCC target ranks 可重算。
2. 讀 7-sample VAF 與 exact-PS profile → 驗證：每矩陣 7×7、對稱、對角線固定。
3. 以 6 biological IDs 做同癌種探索性 label-permutation → 驗證：technical 與 biological-ID 去重結果分開。

Checkpoint：若同一主張跨至少兩個不同證據層且 receipt checks pass，才列為支持；只有 marginal distribution／coarse profile 支持時，維持「描述性相似」。

## §5 Gaps

| Missing | Impact | Effort | Priority |
|---|---|---:|---:|
| 第二組以上 same-ID technical replicate | 無法估計一般化 sensitivity | 高 | P0 |
| 單細胞／多區域 clone truth | 無法確認 global parent→child tree | 高 | P0 |
| matched depth／PS fragmentation ablation | 無法量化技術差異的因果占比 | 中 | P1 |
| 更多癌種與 cell lines | 同癌種檢定 power 很低 | 高 | P1 |

## §6 Evidence conflict scan

| Prior conclusion | Tier | Relation | Source |
|---|---:|---|---|
| production strict directed topology、clone count、fusion tree 皆未完成 | L1 | 限制本輪主張上限 | `docs/experiments/INDEX.md` |
| endpoint edge 是無向 read linkage，不是 evolutionary edge | L1 | 防止把 edge overlap 寫成 ancestry | `docs/experiments/INDEX.md` |
| bulk／regional topology 不可直接寫成 confirmed clone ancestry | L2 | 限制 biologic interpretation | `docs/CURRENT_FOCUS.md` |

## §7 Verdict

**GO（報告與分層驗證）；PROBE（同癌種效應）；NO-GO（宣稱已確認一棵全域 clone/subclone tree）。**

Red-team failure modes：

1. HCC technical duplicate 造成 same-cancer 組重複計權。
2. 所有樣本共享 k=1/2、Direct-only 等常見形狀，使粗拓撲相似但位點完全不同。
3. phase block fragmentation／coverage 改變 edge opportunity，使 Jaccard 低估 conditional containment。

Decision lock：只有新增獨立 same-ID replicate 或正交 clone truth，才能把本輪結論升格為一般化 clone-tree validation。
