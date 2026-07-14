<!--
建立時間: 2026-07-15 Asia/Taipei
目標: 在 layered workstation 新增拓樸型態、可辨識度與基因體位置三維導讀前完成 pre-decision audit
處理範圍: 1 個 cohort index + 7 個 canonical v5 dataset pages；chr1-22；desktop/mobile Chromium
cycle_id: 20260715-layered-workstation-concept-dimensions
topic: layered-workstation-concept-dimensions
status: verdict_GO
audit_version: 0.1
build_branch: research/subclonal-reconstruction-202606
build_commit: 6067568
worktree: /big7_disk/liaoyoyo2001/InterSubMod
data_sources: research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py,docs/methodology/_assets/topology_workstation/HCC1395.html
關聯檔案:
  - InterSubMod/.claude/skills/pre-decision-audit/SKILL.md
  - InterSubMod/research/20260714_layered_workstation_redesign/implementation-notes.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/12_layered_workstation_postimplementation_redteam.md
-->

# Pre-Decision Audit: Layered Workstation Concept Dimensions

> **Verdict: GO（80/100）** — 沿用舊頁「先解釋概念、再用數字與篩選下鑽」的教學形式；禁止沿用舊 clone-first 分類、舊分母與舊資料。

## Frontmatter

- **Triggered by**: new-spec / Task B comprehensive validation
- **服務目標**: G3（read-level evidence 語意）、G4（7 datasets 一致）、G5（可重生與可稽核）
- **Last updated**: 2026-07-15
- **Cycle ref**: `InterSubMod/research/20260715_layered_workstation_concept_dimensions/`

## §0 Cynefin Domain Gate

- **Domain**: Complicated
- **Test**: 同一 generator-first、canonical-v5、Playwright QA 流程已於 2026-07-14 重複產生可預測結果。
- **Rationale**: 資料契約已固定；風險主要是資訊架構、語意邊界與 responsive 呈現，需要專家檢視但不需研究型 probe。

## §1 Observation Completeness Checklist

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| 使用者點名三個舊頁維度在新頁不夠清楚 | ✓ | L1 user specification | `InterSubMod/docs/methodology/_assets/topology_workstation/HCC1395.html:109` |
| Current renderer 已正確計算 C、Topo、complete/incomplete/no-primary | ✓ | L1 | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py:90` |
| Current page 已保留 chr1-22 與 coordinate overlap search | ✓ | L1 | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py:648` |
| 舊 topology/determinacy 定義與 current candidate-set contract 相容 | ✗ | L1 counterexample | 舊頁 `single/linear/branched/star` 與 A/B/C determinacy 是 archived clone-first universe，不可直接移植 |
| Canonical v5 7/7 sample pages 可由單一 renderer 重生 | ✓ | L1 | `InterSubMod/research/20260714_layered_workstation_redesign/implementation-notes.md` V-001/V-007 |

## §2 Credibility Score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | 三維導讀符合 overview → filter → detail-on-demand，且沿用既有視覺語言 |
| 觀察支撐 | 20 | 使用者明確缺口 + current renderer/source + 7/14 全頁 QA |
| 機制清晰度 | 20 | 每個維度可綁定 canonical 欄位與既有互動入口 |
| 反例風險 | 10 | 舊術語會把 regional candidate tree 誤讀為 biological clone/ancestry |
| 所需資源 | 10 | 全頁重生與 24-context QA 預估 1-6 小時 |
| **TOTAL** | **80 / 100** | **GO** |

**Falsifier observable**：若三維導讀無法由 canonical 欄位重生，或新增後仍無法在 5 秒內回答「形狀是否唯一、辨識到哪一層、位於哪裡」，或任一 7-sample page 出現數值／overflow／互動 regression，則此設計失敗並回退。

**Reality-test 三個反例觀察**：

1. 若任何文案將 `Topo=1` 寫成 exact tree 唯一或 biological ancestry 已確認，語意 gate 失敗。
2. 若「位置觀察」宣稱 hotspot／偽影／富集而沒有外部註解與統計，科學邊界失敗。
3. 若 mobile 只能看到卡片標題而看不到數字、座標或下鑽入口，資訊架構失敗。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant |
|---|---|---|---|
| 舊頁只作教學形式參考，不作資料真值 | HIGH | KNOWN | verify quickly |
| `C/Topo/topology_class` 可作 current determinacy 的完整分層 | HIGH | KNOWN | verify quickly |
| start/end 是可顯示的 inclusive genomic coordinate | HIGH | KNOWN | verify formula |
| chromosome counts 可描述分布但不能推論 biological enrichment | HIGH | KNOWN | document boundary |
| emoji 是理解所必需 | LOW | UNKNOWN | defer；沿用既有無 emoji 視覺系統 |

## §4 Quick Pilot Guide

GO 已成立；正式實作仍先以 HCC1395 desktop/mobile capture 作視覺基準，再 generator-first 擴至全部 7 pages。

1. 擷取 old/new HCC1395 同 viewport → 驗證：兩張 accepted screenshot 均可見三維入口或其缺口。
2. 將三維卡綁到 `topology_class`、`C/Topo`、`chrom/start/end` → 驗證：每欄皆可追溯 source field。
3. 全量重生 → 驗證：7/7 hash-bound pages、無手打 sample 數字。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| Current-run old/new same-state screenshot comparison | HIGH | 0.5h | P0 |
| Region-level 三維摘要與名詞欄位對照 | HIGH | 1-2h | P0 |
| 7-sample desktop/mobile regression | HIGH | 1-2h | P0 |
| Screen-reader 實機 | MED | external | P2 |

## §6 Evidence Conflict Scan

- Repo root 沒有 `MEMORY.md`，因此強制 grep 以 **source absent** 記錄；未用記憶補猜。
- `InterSubMod/docs/methodology/_assets/topology_workstation/` 已由 commit `82bd02e` deprecate，與「直接移植舊定義／舊數字」衝突。
- `InterSubMod/docs/CURRENT_FOCUS.md` 明示 canonical v5 regional mutation-state tree 不是 cellular clone/ancestry truth，支持 current 邊界。
- 既有 post-implementation red-team 支持 C/Topo 與 primary/auxiliary contract，但不表示三維導讀已清楚。

## §7 Decision Threshold + Path

- **TOTAL**: 80 / 100
- **Verdict**: **GO**
- **Decision lock**: Yes — canonical v5 only；舊頁只保留教學結構；全基因視圖不可移除；JSON 連結保持收合。
- **Next action**: capture baseline → implementation notes → renderer/index → 7-page rebuild → Claude Code read-only review → Chromium Task B QA。

### Independent red-team

1. **Failure mode A**: 用舊 `single/linear/branched/star` 讓讀者以為已重建真實 clone lineage。防線：只使用 exact candidate / unlabeled topology shape / incomplete / no-primary。
2. **Failure mode B**: 三張卡只是重複既有數字，沒有連到篩選與 region detail。防線：卡片必提供 current observation、how-to-read、限制與操作入口。
3. **Failure mode C**: 把 chromosome 描述升格為 hotspot 或 artifact risk。防線：只報 descriptive count/rate，明示非 enrichment test。

紅隊未觸發降級；verdict 維持 GO。

## Provenance Footer

- **Commit hash**: `6067568`
- **Build date**: 2026-07-15 Asia/Taipei
- **Skill version**: `/pre-decision-audit` v0.1
- **Canonical summary SHA-256**: `71da78b69f8afe5fb8e618179ab7b38a6940fdb17be6282f99a6ec4b720e5de7`
