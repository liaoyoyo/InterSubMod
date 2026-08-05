<!--
建立時間: 2026-07-18 00:00 +08:00
目標: 在修改 Hypercube exact solver 與 ranking cache 前，先固定 exact-preserving 邊界、反證條件與驗證門檻
處理範圍: prepared base MILP reuse、complete-only structural enumeration cache、diagnostics；不含 coverage/estimand 變更
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2正式pilot_NO_GO驗證紀錄_01.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/partial-read-method-audit.md
-->

# Hypercube exact-preserving solver/cache remediation — Pre-decision audit

> **任務類型**：E Hotfix / performance remediation。  
> **Task scope**：只修重複建模與完整結果重算，不改 group constraints、minimum-extra objective、candidate universe、likelihood、T/V/Topo 或任何生物分母。  
> **狀態**：`complete`；三方唯讀審查已完成，鎖定為 bounded `PROBE`。

## §0 Cynefin domain gate

- **Domain**：Complicated。
- **理由**：相同結構問題的候選集合與 `h*` 有小 k exhaustive oracle；可用逐案相等檢查判斷 refactor 是否正確。
- **禁止擴張**：不把目前 SciPy stateless MILP 包裝成真正 persistent solver；本輪只承諾「base matrix prepared once」。

## §1 Observation completeness

| Observation | 狀態 | Tier | 來源 |
|---|---|---:|---|
| 每次 `solve_min_hidden()` 都呼叫 `_build_problem()` | ✓ | L1 | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py:412-430` |
| optimal-set enumeration 每找到一個集合就再呼叫一次 `solve_min_hidden()` | ✓ | L1 | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py:466-529` |
| base predecessor/group rows在每次 enumeration solve 都重建 | ✓ | L1 | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py:285-408` |
| frozen v2第一pilot B0 於8小時timeout，正式verdict NO-GO | ✓ | L1 | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2正式pilot_NO_GO驗證紀錄_01.md` |
| minread=1 diagnostic中33個長尾units占candidate-generation時間96.153683% | ✓ | L2 diagnostic | 同上；只解釋效能，不是生物比例 |
| 跨units／minread的實際cache hit率 | □ | L4 | 本輪只建立diagnostic；效能收益須由新pilot量測 |

## §2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20/20 | 固定base constraints後只附加objective equality/no-good rows，數學模型不變 |
| 觀察支撐 | 20/20 | source trace與正式timeout evidence直接指出重複build/solve長尾 |
| 機制清晰度 | 20/20 | prepare→solve first→append exclusions→solve next，可逐constraint比對 |
| 反例風險 | 10/20 | sparse no-good、forced groups、timeout/max-set與mutable cache copy皆可能出錯 |
| 所需資源 | 20/20 | targeted oracle與synthetic benchmark可在30分鐘內完成 |
| **exact-preserving patch可信度** | **90/100** | 可進入受限實作與等價性測試 |

> **Scale-level獨立紅隊分數：40/100（PROBE）**。原因不是數學 refactor
> 不可測，而是尚無新 source 的 33-tail／chr6／H2009 實測；不得把單元測試
> 或 Python-side build 減少誤寫成正式 scaling 已解鎖。

**Falsifier observable**：任一案例出現 `h*`、`complete`、stop status、sorted vertex-set IDs、partial-group coverage或likelihood結果與baseline不同，即否決本patch；cache若重用`complete=false`結果也立即否決。

## §3 Assumption map

| Importance × knowledge | Assumption | Action |
|---|---|---|
| High × known | base predecessor/group rows在同一structural problem內固定 | source＋matrix equality test |
| High × known | complete result與time limit無關，但本輪key仍納入solver contract參數 | canonical tuple key＋payload equality |
| High × unknown | 真實pilot跨unit/minread有足夠cache hit率 | 只記錄hit/miss；不先承諾加速倍數 |
| High × unknown | SciPy/HiGHS能否保留solver instance | 本輪不宣稱；只reuse prepared sparse base |
| Low × known | cache entry需deep copy避免caller mutation | mutation-isolation test |

## §4 Quick pilot

1. Baseline：重跑`test_hypercube_exact.py`與`test_build_m2_patterns_and_rank.py`。  
   → 實際結果：`Ran 51 tests in 2.638s`，`OK`，exit 0。
2. Probe：對具有24個optimal sets的固定案例計數`_build_problem()`呼叫次數。  
   → 修改前實際結果：`AAAA, k=4, all_loci`完整取得24個候選，
   `_build_problem()`呼叫25次，wall 0.266119秒；修補後prepared-base build須固定為1。
3. Equivalence：k≤6 independent exhaustive oracle＋固定RRA/AAA/RAX、AX/XA、timeout/max-set負向測試。  
   → Checkpoint：candidate catalog 0 mismatch、failure semantics不變。

## §5 Gap diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---:|
| 真實長尾33 units在新source的wall-time | 無法宣稱正式pilot已解鎖 | >8h風險 | 後續新release雙pilot |
| cache hit率與節省時間 | 無法承諾cache收益 | <30min synthetic；真實需pilot | P1 diagnostic |
| 真正persistent HiGHS model | 無法消除每次solver presolve | 需highspy/API spike | P2 research |
| per-unit deadline | coverage會改變 | 中 | 另立contract，不在本輪 |

## §6 Evidence conflict scan

- `MEMORY.md`在repo root不存在；以`InterSubMod/docs/CURRENT_FOCUS.md`、cycle `00_INDEX.md`、NO-GO紀錄與validated NEGATIVE索引替代。
- 唯一相關NO-GO是frozen v2 ranking timeout，屬**支持修補必要性**，不是反對exact-preserving refactor。
- `InterSubMod/docs/reports/validated/2026/04/20260401_LOH_weekly_review/06_methylation_hypothesis_negative.md`與本工程效能patch無直接依賴。
- 不重新開啟已否決的methylation hard-filter方向，也不把HP/VAF/甲基加入edge cost。

## §7 Decision path

- **最終 verdict**：`PROBE`。
- **局部工程判定**：
  - Python-side prepared-base refactor：`GO_WITH_EXACT_EQUIVALENCE_GATES`。
  - process-local complete-only structural cache：`GO_WITH_CACHE_SAFETY_GATES`。
  - 宣稱 persistent HiGHS model、正式 chr6 scaling已解鎖或啟動full 154：`NO-GO`。
- **Decision lock**：只允許prepared base reuse、complete-only cache與diagnostics。
- **硬中止條件**：
  1. 任一semantic output mismatch；
  2. cache保存或回傳incomplete結果；
  3. timeout/max-set被誤標complete；
  4. 需要更動candidate universe、objective或ranking分母；
  5. cold/cache的ordered candidate IDs、likelihood或semantic hash出現差異。
- **紅隊追加 gates**：
  1. cache dictionary以完整frozen tuple判等；SHA僅供診斷；
  2. `complete=false`、limit、numerical failure與`MAX_SETS_REACHED`一律不存；
  3. cache store與hit各做一次deep copy；
  4. likelihood、BQ、fixed-error與bootstrap結果永不快取；
  5. timeout後不得把partial candidates排成winner。
- **下一步**：source patch → targeted/full regression → 33-tail regression → 新release
  chr6＋H2009雙pilot；本輪不得碰frozen v2失敗root。

## §8 Independent review record

| 視角 | Verdict | 關鍵限制 |
|---|---|---|
| prepared-base solver audit | GO（局部） | `_build_problem`可降25→1，但SciPy `milp()`與HiGHS model仍每輪重建 |
| structural cache audit | GO_WITH_GATES | process-local、完整tuple、complete-only、deep-copy、likelihood不快取 |
| scale/red-team audit | PROBE | 只允許新release bounded implementation；正式雙pilot前不得升級GO |

三方均為唯讀審查，未修改source。最強共同反證條件是：任何`h*`、`complete`、
ordered `vertex_set_ids`、candidate likelihood、winner/tie或scientific semantic hash差異。
