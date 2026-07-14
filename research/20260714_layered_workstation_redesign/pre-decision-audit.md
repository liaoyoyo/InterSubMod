<!--
建立時間: 2026-07-14 Asia/Taipei
目標: 在 layered_workstation 全站重構前，鎖定 canonical 資料、展示邊界、失敗條件與驗收路徑
處理範圍: InterSubMod/docs/methodology/_assets/layered_workstation/index.html 與 7 個 per-sample 工作站；chr1-22 全基因視圖
關聯檔案:
  - InterSubMod/docs/CURRENT_FOCUS.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json
  - InterSubMod/research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation.py
cycle_id: 20260714-layered-workstation-redesign
topic: 20260714_layered_workstation_redesign
status: verdict_GO
audit_version: 0.1
task_type: B comprehensive validation
-->

# Pre-Decision Audit: Layered Workstation 全站重構

> **Verdict: GO（80/100）**。現行介面可運作，但資料來源仍是 2026-07-10 historical snapshot；2026-07-14 LongPhase-S PASS canonical v5 已有 7/7 完整 region-view，可用 generator-first 方式安全升級。研究結論只允許 regional mutation-state candidate tree，不升級為 cellular clone/ancestry。

## §0 Cynefin Domain Gate

- **Domain**: Complicated；視覺層以可逆的 Chromium probe 反覆收斂。
- **Test**: 同一份 canonical region-view 經 deterministic builder 可重複得到相同 census 與 region payload；版面與互動則需桌機／手機實測。
- **禁止事項**: 不手改上游 JSON，不把 ClairS PASS sensitivity v6 取代 canonical v5，不把 L3 缺測當 neutral，不把 read AF 寫成 CCF 或 tree posterior。

## §1 Observation Completeness

| Observation | 狀態 | Tier | 來源 |
|---|---:|---|---|
| 7/7 canonical v5 layered_region_view JSON 存在 | ✓ | L1 | canonical root samples/*/layered_region_view_*.json |
| canonical aggregate W_tree=51,815、W_primary=50,215 | ✓ | L1 | current_layered_topology_v3_raw_all_v1.json |
| complete/incomplete=42,240/7,975；三類 C/Topo=11,582/10,737/19,921 | ✓ | L1 | 2026-07-14 latest analysis |
| 現行 builder SAMPLES 仍指向舊 pilot／multisample 路徑 | ✓ | L1 | build_layered_per_sample.py |
| 現行首頁顯示 data snapshot 2026-07-10 與舊 census | ✓ | L1 | layered_workstation/index.html |
| 8 頁桌機／手機 current-run 視覺 baseline | □ | — | Playwright capture 進行中 |
| Claude Code 與獨立 IA／資料紅隊 | □ | — | 本輪完成後寫入 audit notes |

## §2 Credibility Score

| 維度 | 評分 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | assertion–evidence、progressive disclosure、WCAG 與 canonical data contract 明確 |
| 觀察支撐 | 20 | 7/7 v5 region-view 與 machine summary 已存在 |
| 機制清晰度 | 20 | stale source、錯誤術語、頁面層級與 responsive 問題均可在 generator 定位 |
| 反例風險 | 10 | 138MB 離線頁面、舊欄位相容、candidate set/capped precedence 需 browser regression |
| 所需資源 | 10 | 估計 1–6 小時的 generator、全頁重生與 QA |
| **TOTAL** | **80/100** | **GO** |

**Falsifier observable**：任一頁仍顯示舊 ClairS PASS 主幹、2026-07-10 snapshot、與 canonical v5 不一致的 7-sample 數字；任一 sample 缺 chr1-22 全基因入口；C/Topo/complete/incomplete 守恆失敗；network 把推測祖先畫成實測、把 capped 當唯一；390px/1440px overflow、console error、broken link 或 keyboard dead-end，則此次重構失敗。

**Reality-test 三反例**：

1. 只更新 index 會讓首頁與 7 個 sample 頁互相矛盾，故需 8 頁同版資料契約。
2. 增加漂亮 network 圖可能弱化 evidence state，故 observed/inferred、forced/variable、complete/incomplete 必須同時以文字與非色彩 encoding 顯示。
3. 全基因 overview 若只做聚合圖會遮蔽可追溯 region，故 chromosome overview 必須能直接篩選並回到 region detail。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant / action |
|---|---|---|---|
| canonical v5 是唯一主結果；v6 只是 sensitivity | HIGH | KNOWN | 鎖定 |
| 7 個 region-view 的 census 與 regions 同一 snapshot | HIGH | KNOWN | hash／manifest 驗證 |
| C_region 與 Topo_region 可由 analysis-complete primary HP1/HP2 units 重算 | HIGH | KNOWN | generator contract test |
| 全基因 view 指 chr1-22，不含 chrX/Y primary claim | HIGH | KNOWN | scope badge + chromosome navigator |
| 138MB 頁面在手機載入與互動仍可接受 | HIGH | UNKNOWN | Chromium 實測；必要時減少重複 DOM |
| L3 可稱 pending | HIGH | KNOWN FALSE | 依 payload顯示 not_evaluated / bounded auxiliary |

## §4 Quick Pilot / Checkpoint

1. 先用 HCC1395 canonical v5 生成單頁 smoke。
   → 驗證：summary、chr1-22、C/Topo 與 payload 重算一致；console error=0。
2. 在 1440×1000 與 390×844 檢查 index + HCC1395。
   → 驗證：body overflow=0；導航、篩選、region deep link、tree switch、keyboard focus 全通。
3. 通過後重生 7/7 sample pages。
   → 驗證：8 頁 links 7/7、embedded JSON 7/7 parse、canonical hash/provenance 7/7 顯示。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| builder 尚未指向 canonical v5 | HIGH | <1h | P0 |
| index 缺最新 C/Topo 與 full-genome cohort map | HIGH | 1–2h | P0 |
| sample 頁缺 chromosome-level genome overview | HIGH | 1–2h | P0 |
| no_primary_lineage、L3、backbone 等新 schema 顯示相容 | HIGH | 1h | P0 |
| CCF/posterior、clone-c 等過度語意需移除或降級 | HIGH | 1h | P0 |
| mobile comparison、focus、large-DOM performance QA | MED | 1–2h | P1 |

## §6 Evidence Conflict Scan

| Prior conclusion | Tier | Relation | Source |
|---|---|---|---|
| 7/11 historical W=48,959 不再代表 current canonical | L1 | conflict with current UI | 2026-07-14 latest analysis |
| bulk regional candidate tree 不確認 cellular clone K/ancestry | L1/L3 boundary | constrains wording | CURRENT_FOCUS.md |
| L3 methylation 未 join；禁止 rank/confirm tree | L1 contract | constrains UI | canonical region-view L3_methyl |
| PS 只作 phase-block QC，不作 topology edge | L1 contract | constrains network | canonical analysis_contract |
| sensitivity verdict=backbone_sensitive | L1 | requires explicit comparison caveat | current machine summary |

## §7 Decision Path

- **Verdict**: GO。
- **Next action**: baseline screenshot → independent red team → generator-first implementation → HCC1395 smoke → 7/7 build → full Chromium/data QA。
- **Decision lock**: Y。
- **Revisit if**: canonical root 改版、region-view schema 變更、或取得正交 single-cell／multi-region truth。

### Independent red-team gate

- Failure mode 1：index 更新而 sample HTML 沒重生。Gate：8 頁 generated-at/backbone/snapshot contract。
- Failure mode 2：C/Topo、determined、candidate-complete 混成同一指標。Gate：分開顯示 candidate count、shape count、completeness。
- Failure mode 3：手機只看見圖的一小角或列表把 detail 推到 30k px 後。Gate：chromosome overview 可捲動／折行、selected detail 在 mobile 可直接聚焦。
- Failure mode 4：為美觀隱藏重要限制。Gate：claim ceiling、L3/PS/CN/partial-read boundary 在首頁與 sample 首屏都可見。
