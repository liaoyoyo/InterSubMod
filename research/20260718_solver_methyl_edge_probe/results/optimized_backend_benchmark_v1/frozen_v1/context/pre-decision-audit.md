<!--
建立時間: 2026-07-18 02:24 +08:00
目標: 小步驗證 persistent MILP、uncovered-group exact B&B、vertex-set-level read-AF cache 與 methylation edge tie-breaker 是否值得進入正式實作
處理範圍: synthetic exact cases + 至多兩個既有 real capped units + HCC1395_DORADO chr6 methylation feasibility scan；PARTIAL exploratory pilot
cycle_id: cycle_20260718-solver-methyl-edge-probe
topic: 20260718_solver_methyl_edge_probe
status: verdict_PROBE
audit_version: 1.0
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/pre-decision-audit.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py
  - InterSubMod/docs/CURRENT_FOCUS.md
-->

# Pre-Decision Audit：solver 與甲基 edge 小步驗證

> **Verdict：PROBE（50/100）**。允許在全新隔離目錄做 synthetic 與極小 real pilot；不允許替換 canonical solver、全樣本執行，或把甲基距離寫成祖先／真實 edge 證明。
>
> **服務目標**：G3（read-level epigenetic evidence）與 G5（可重現、可外部檢查的 exact solver）。

## §0 任務分類與範圍

- **Task Type**：A — Exploratory pilot。
- **PARTIAL scope**：synthetic k≤5 exact oracle、`AAAA,k=4` 24-optima stress、最多兩個既有 capped real units，以及 HCC1395_DORADO chr6 methylation feasibility。
- **Domain**：Complex；採 safe-to-fail probe。
- **禁止推廣**：本輪結果不得代表 chr1-22、全樣本 runtime、edge accuracy 或 clone/subclone accuracy。

### 研究啟動五問

1. **Thread D？** 是；甲基只作同一 vertex set 內的 exploratory edge regularizer。
2. **Thread B 撤回範圍？** 不直接涉及 FP filter；但既有「methylation cluster ≠ biological subclone」結論仍是硬邊界。
3. **KDE-corrected？** 不適用於結構 solver；甲基資料直接使用 MM/ML，不用 KDE 結論。
4. **VCF caller AF？** 不需要；現行 primary score 是 read-pattern mixture likelihood。不得用相同 reads 再重複加入 VAF。
5. **長計算／C++／搬移／NO-GO？** 不改 C++、不搬檔；每個 pilot 設短時限，只寫本研究目錄與 `/tmp` 隔離依賴。

## §1 已知觀察

| 觀察 | 狀態 | 證據 |
|---|---:|---|
| 現行 SciPy `milp` 每找到一個 candidate 會重建整個模型 | ✓ | `hypercube_exact.py:285,412,466` |
| `AAAA,k=4` 有 24 個最佳 vertex sets，因此現行路徑 build/solve 25 次 | ✓ | 獨立 read-only audit |
| SciPy 1.13.1 wrapper 沒有可增量 `addRow` 的 persistent handle | ✓ | 環境與 API audit |
| 本環境目前沒有 `highspy` | ✓ | `ModuleNotFoundError` |
| persistent model 只能減少 structural rebuild，不能消除 V+1 次求解與 V 個輸出 | ✓ | no-good enumeration contract |
| 現行 `rank_unit` 已先按 `vertex_sets` 做一次 mixture fit，再解析計數 parent assignments | ✓ | `build_m2_patterns_and_rank.py:886+` |
| 同 vertex set、不同 parent edges 的 read likelihood 必須同分 | ✓ | 現行模型註解、既有 regression |
| HCC1395_DORADO chr6 在 MINREAD 1/2/3/5 都沒有同一 component×HP 同時出現 10、01、11 exact states | ✓ | 本輪 lossless sidecar feasibility scan |
| HCC1395 是 hyper-diploid（ploidy 2.85），大量 gain/LOH；甲基 edge 正例需 CN gate | ✓ | KB `02_samples/hcc1395.md`，last verified 2026-07-11 |
| read×CpG cluster 是 methylation similarity，不是已驗證 biological subclone | ✓ | KB `05_tools/intersubmod.md`，last verified 2026-07-11 |

## §2 Credibility Score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | MILP incremental cuts、exact branching completeness、vertex-set likelihood invariance皆可形式化 |
| 觀察支撐 | 10 | 已知 build bottleneck、real MILP optimum pilot與現行 cache路徑，但尚無新 backend實測 |
| 機制清晰度 | 10 | solver三層可分；甲基 edge score仍依賴可逆性較弱的相似性假設 |
| 反例風險 | 0 | 舊 full pilot曾8小時 timeout；tie數量仍可爆炸；CNA/LOH與甲基收斂可使 edge誤導 |
| 所需資源 | 10 | synthetic短；新 backend與real enumeration仍有版本、尾端runtime及輸出量風險 |
| **TOTAL** | **50 / 100** | **PROBE** |

## §3 Falsifiers 與停止條件

1. persistent 與現行 SciPy/exhaustive oracle 的完整 vertex-set digest 任一不一致 → **FAIL**。
2. 未先固定 `sum(extra)=h*` 就加入 sparse no-good，或把高一層合法解錯刪 → **FAIL**。
3. B&B 在 seeded k≤5 任一 case 的 optimum／完整最佳 set 不等於 MILP/exhaustive oracle → **FAIL**。
4. B&B node count/runtime在強 partial 小案例沒有改善，或記憶體因 memo/output迅速放大 → 僅保留替代 backend，不升 canonical。
5. primary read-AF fit呼叫次數不等於 distinct vertex sets，或不同 edge得到不同 read-AF → **FAIL**。
6. synthetic methylation「11接近10／接近01／等距」不能分別回復 T10／T01／TIE → **FAIL**。
7. 真實資料缺 MM/ML、三個完整 state、共同 CpG、HP×PS一致或 CN-clean gate → **ABSTAIN**，不得補值。

## §4 Assumption Map

| 假設 | 重要性 | 已知？ | Action |
|---|---:|---:|---|
| 新 `highspy` 與 SciPy內嵌 HiGHS 給相同整數解集合 | High | No | pin version；完整 set digest雙 backend對照 |
| sparse no-good在固定h*後只排除一個 set | High | Yes | regression含高一層反例 |
| uncovered-group branching涵蓋每個可行延伸 | High | Yes | formal proof + k≤5 exhaustive oracle |
| hidden node少、partial constraints強是常見實例 | High | Partial | small real capped units只作描述 |
| snapshot read likelihood可分parent edge | High | Yes=False | 同V只算一次，edge維持 unresolved |
| methylation沿真edge保持相似且不收斂／反轉 | High | No | synthetic sensitivity；real資料只稱 regularizer |
| HP/PS已排除copy-number混淆 | High | No | CNA/LOH/allele-specific amplification gate |

## §5 Quick Pilot

1. **Persistent solver**  
   一次建模求 `h*`，固定 cardinality，再增量加入 sparse no-good rows。  
   → 驗證：`AAAA,k=4` 完整24 sets；model builds=1；sets digest等於現行25-build結果。
2. **Exact uncovered-group B&B**  
   對 orphan predecessor 或最小 uncovered group 做完備分枝，使用安全 lower bound與 memo。  
   → 驗證：seeded k≤5 optimum與完整sets 0 mismatch；回報 visited/pruned/runtime。
3. **Vertex-set-level read-AF**  
   對相同V的兩個 parent assignments插入instrumentation。  
   → 驗證：fit calls=distinct V=1、parent choices=2、score完全相同。
4. **Methyl edge tie-breaker**  
   `Δ = D_M(01,11) − D_M(10,11)`；`Δ>0`只表示模型偏好 `10→11`。  
   → 驗證：正向、鏡像、等距、uninformative與missing-CpG synthetic controls。
5. **Real feasibility gate**  
   先查同HP×PS的10/01/11、MM/ML、共同CpG與CN/LOH。  
   → 驗證：不合格即明確 `ABSTAIN_NO_EVALUABLE_REAL_UNIT`。

## §6 Conflict Scan 與 Red Team

| 既有結論／反例 | 關係 | 處理 |
|---|---|---|
| HCC1395 chr6 formal M2排名曾8小時timeout | 直接風險 | 本輪只跑短 synthetic／單元，不升full |
| methylation-only subclone usability為negative | claim conflict | 只測edge-similarity scorer，不稱clone |
| same-V/different-edge對snapshot sSNV不可辨識 | 方法硬限制 | read-AF只算一次；甲基另列exploratory score |
| allele-specific amplification可在同cell/同HP產生10與11 | 生物反例 | CN/LOH confounded即ABSTAIN |
| methylation可逆、收斂或受strand/CpG coverage驅動 | 模型反例 | coverage matching、bootstrap/permutation與TIE |

獨立 red team 的最強反方：

- **Persistent**：仍需 V+1 次 `run()`；候選全集很大時，重用模型不能解決 output explosion。
- **B&B**：若 group subcubes高度重疊且下界弱，仍會接近 `O(2^M)`。
- **甲基**：即使 `Δ` 顯著，也可能是 CNA、LOH、甲基狀態收斂或 technical geometry，不是真parent。

## §7 決策

- **允許**：全新 research scripts/tests、`/tmp` pinned dependency、synthetic與至多兩個small real units。
- **不允許**：修改 canonical C++/Python、全資料執行、publication claim、真實edge／clone／subclone宣稱。
- **升級條件**：正確性全PASS，且至少一項 solver在預註冊case有實質改善；甲基需另找到滿足三state+共同CpG+CN-clean的real unit才可升下一個PROBE。
- **Claim ceiling**：solver-certified regional mutation-state candidate sets；甲基最多是 methyl-similarity regularized candidate-edge preference。

## Provenance

- Git commit：`0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Current exact solver SHA256：`9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95`
- Current ranker SHA256：`c82210f25870d1405cc070c18096fb7d1c2b908fb4d8a7858aece7ac4b151d28`
- Runtime baseline：Python / SciPy 1.13.1；`highspy` absent at audit。
- Resource：`/big7_disk` 99% used、716G available；pilot outputs必須保持小型。

