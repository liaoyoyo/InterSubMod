<!--
建立時間: 2026-07-23 00:00 +08:00
目標: 在 production pipeline 導入 exact-PS×HP strict endpoint read-linkage region 前完成可否證的決策前審計
處理範圍: Python/C++ region grouping、production runner、tree adapter、HCC1395 chr1-22 與後續 7-dataset promotion gate
關聯檔案:
  - InterSubMod/research/20260722_hcc1395_ps_block_hp_orientation_audit/20260722_HCC1395跨PS區域與HP方向敏感度驗證_01.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/lossless_read_contract.py
  - InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/scripts/exact_ps_k12_partition.py
-->

# Production exact-PS strict read-linkage：Pre-decision audit

## §0 Cynefin front-gate

- Domain：**Complicated → Complex boundary**。Exact PS 隔離與 endpoint co-call graph 是明確工程契約；但 `MINREAD=3`、partial-read 證據分散與 k>12 切割對拓撲的影響需以實際資料 probe。
- Task type：**C — Production deployment**。
- 服務目標：**G1、G4、G5**；G3 僅作後續 likelihood／甲基輔助介面，不在本次宣稱範圍。

## 啟動研究任務五問

1. Thread D read-level epigenetic：否；本輪是 genetic linkage 基礎層。
2. Thread B 已撤回方向：否。
3. KDE-corrected：不適用。
4. Caller AF：region grouping 不使用；VAF/CCF 留在候選後驗證層。
5. 長計算／C++／搬移／gate：**是**；涉及 Python/C++、chr1-22 BAM-derived artifacts 與 production promotion gate。

## §1 Observation completeness

| Observation | 狀態 | Tier | 證據 |
|---|---|---|---|
| 舊 mixed-PS grouping 對 block orientation 敏感 | ✓ | L1 | HCC1395 865 mixed-PS regions；flip probe 可使 T 1→3、Topo 1→2 |
| exact-PS pilot 無跨 PS／HP | ✓ | L1 | HCC1395 chr1-22：cross_PS=0、cross_HP=0、Python/C++ mismatch=0 |
| cut-span threshold=3 不等於 pairwise endpoint threshold=3 | ✓ | L1 | 205/11,435 multisite components、127/9,600 primary blocks 會裂開 |
| 所有現行 multisite components 至少存在 threshold=1 read chain | ✓ | L1 | HCC1395 chr1-22 strict pairwise census |
| 7 datasets 皆具備可直接重算的 exact-PS molecule artifacts | □ | L5 | 本輪先盤點，未確認前不得宣稱 full promotion |
| strict graph 重跑後 topology 影響 | □ | L5 | 必須重建 region、partition、candidate trees 後才可計算 |

## §2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | PS 是局部 phase block；endpoint co-observation 是直接可觀測 linkage |
| 觀察支撐 | 10 | HCC1395 全 chr1-22 已重算，但尚未跨樣本 |
| 機制清晰度 | 20 | cut-span 過度連接與 adapter MINREAD 不一致已有可重現反例 |
| 反例風險 | 10 | strict pair threshold 可能丟失分散 partial evidence，需 sensitivity |
| 所需資源 | 0 | production + 7-dataset full validation 預期 >6 小時 |
| **TOTAL** | **60/100** | **GO implementation；跨樣本 promotion 另設 gate** |

Falsifier observable：若 strict implementation 後仍出現 cross-PS／cross-HP、任一輸出 region 在 endpoint graph 不連通、Python/C++ partition 不一致，或 read/membership 質量不守恆，則 production promotion 失敗。

## §3 Assumption map

| 重要性 | 已知 | 未知 |
|---|---|---|
| 高 | exact PS×HP 不可跨容器；canonical molecule 必須去重；endpoint 必須 fixed R/A | 7 datasets exact artifacts 可用性；strict threshold 對 topology 的跨樣本影響 |
| 低 | threshold 1/2/3/5 可作 sensitivity | 是否需要未來以統計 edge confidence 取代固定 threshold |

## §4 Quick pilot

1. HCC1395 chr1 synthetic＋chr1 smoke → 驗證：strict graph、PS/HP、mass conservation、Python/C++ parity 全 PASS。
2. HCC1395 chr1-22 → 驗證：每個 region endpoint graph connected；cross PS/HP=0；正式 census 完整。
3. 7-dataset inventory → 驗證：輸入身份、LongPhase-S PASS、latest HP/PS sidecar 與必要 molecule rows 均可追溯；缺任一項則停止 promotion。

## §5 Gap diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| Production runner 尚走 legacy distance/HP-only path | 正式結果仍可能跨 PS 或無 read linkage | 中 | P0 |
| strict pairwise Python/C++ 共用契約 | 無法保證 parity | 中 | P0 |
| 7-dataset exact molecule inventory | 無法發布全樣本數量 | 中 | P1 |
| likelihood／CN/CCF ranking | 不能確認真實樹，只能列候選 | 高 | P1，非本次 grouping blocker |

## §6 Conflict scan／red team

- 反方一：pairwise `>=3` 可能比 cut-span 過度保守，將分散於不同 endpoint pairs 的真實 partial evidence拆掉。回應：保存 threshold 1/2/3/5 sensitivity；primary 固定 3，但輸出 edge weights 與裂解原因，後續由 molecule likelihood 檢驗。
- 反方二：A–B、B–C 的 transitive chain 可能由單一弱 articulation edge造成長距離串接。回應：輸出 minimum edge support、bridge/articulation flags 與 leave-one-molecule sensitivity。
- 反方三：k>12 bounded blocks 被誤當生物 region。回應：永久保留 `source_component_id`；region 與 computational block 分層命名。
- 既有 NEGATIVE 結論主要針對 methylation FP filter，與本次 genetic linkage correctness 無直接衝突。

## §7 Verdict

- Verdict：**GO implementation；7-dataset publication promotion = conditional GO**。
- Decision lock：production scientific region 固定為 `chromosome × exact PS × HP × strict endpoint connected component`；未經 signed cross-PS orientation bridge，不得跨 PS 合併。
- 下一步：建立 living implementation notes → 修改 production route → synthetic/parity → HCC1395 chr1-22 → 7-dataset gate。
