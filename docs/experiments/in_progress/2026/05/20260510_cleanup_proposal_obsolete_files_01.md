---
title: InterSubMod 過時檔案 cleanup proposal — 待用戶確認後手動執行
date: 2026-05-10
status: pending_user_review
type: cleanup_proposal
classification: housekeeping
generated_by: AI scan (no deletion executed)
related: InterSubMod/docs/experiments/in_progress/2026/05/20260510_v1.8_implementation_retro_01.md
---

# InterSubMod 過時檔案掃描 — Cleanup Proposal

> **Bottom line**：掃描 `InterSubMod/` 與相關輸出後，找到 **3 類候選**：(1) 確認過時可封存 (Tier 1，~600 KB，2 項) / (2) 可能過時建議 review (Tier 2，~26 MB，3 項) / (3) 已自封存或屬正常使用 (Tier 3，不需動)。**未執行任何刪除/搬移**；下列每項附建議行動 + 完整 path（用戶手動確認後可執行 §6 helper script）。**CJK font 已修復**（plot_setup.py 加 `Noto Sans CJK JP` 到 chain + 改用直接 family list 設定法）。

---

## §1 CJK font 安裝確認

**結果**：✅ **CJK 渲染正常**（修補 `plot_setup.py` 後）。

| 項 | 結果 |
|---|---|
| 系統 fc-list | `Noto Sans CJK JP`（Debian/Ubuntu fonts-noto-cjk 套件 unified font，含全部繁簡日韓 glyphs） |
| 第一次 plot_setup 測試 | ❌ glyph missing warnings（matplotlib 3.6.2 `font.sans-serif` chain fallback 不可靠） |
| 修補：加 JP 到 chain + `font.family` 直接設 actual list | ✅ `applied: ['DejaVu Sans', 'Noto Sans CJK JP', 'Droid Sans Fallback', 'sans-serif']` |
| 最終 render 測試 | ✅ 13.5 KB output / 0 warning / 繁體中文渲染正常（測試文字「測試：中文 + ISM + COLO829\n甲基化分析 ΔF1=+0.0112 樣本基因」）|

**修補內容（已修但未 commit）**：
- `InterSubMod/scripts/lib/plot_setup.py`：DEFAULT_FONT_CHAIN 加 `"Noto Sans CJK JP"`
- 改 `font.family = applied` 而非 `font.family = ['sans-serif']`（matplotlib 3.6.2 caveat workaround）

**建議**：v1.8 後續 commit 「fix(infra): plot_setup CJK fallback uses direct family list (matplotlib 3.6.2 caveat)」。

---

## §2 Tier 1 — 確認過時（建議封存）

### 2.1 `.agents/skills/` 整個目錄（untracked，~600 KB）

**狀態**：repo 內 untracked 但物理存在，git 從未 add。

**過時根據**：
- 與 `.claude/skills/` 部分重複（14 個目錄都在 `.claude/skills/` 找得到）
- 缺 **11 個 newer skill**：`citation-verification` / `cpp-change` / `data-audit` / `doc-standards` / `fast-learning-coach` / `memory-consolidation` / `methodology-audit` / `myPPT` / `observation-analysis` / `pptx-build` / `problem-framing-ideation` / `report` / `results-analysis` / `structured-tech-report` / `verification-loop` / `weekly-report`
- 結論：v1.0 / v1.5 plan 期間建立的舊版 skill 結構，已被 `.claude/skills/` 完全取代

**完整 path**：`InterSubMod/.agents/skills/`（含 `auc-confound-guard / check-staleness / conclude-research / confirmation-protocol / cycle-init / cycle-state / feature-layered-observation / grill-me / init-research / inject-hypothesis / known-pitfalls / multi-sample-consistency / pivot-direction / provenance-tier-audit / research-dashboard / research-loop / results-report / review-evidence / run-evaluator / skills` + `README.md`）

**建議行動**：搬到 `InterSubMod/docs/archive/2026/05/20260510_agents_skills_obsolete/`

### 2.2 `.codex` 空檔案（untracked，0 bytes）

**狀態**：`InterSubMod/.codex`，0 bytes 空檔案，2026-05-08 17:05 建立

**過時根據**：
- codex tool 殘留 marker（已 cross-ref 到 `docs/archive/2026/05/20260509_codex_nested_skills_duplicate/` 的封存案例 — 證實 codex 整合早已被 deprecated）
- 大小 0 bytes，無功能

**完整 path**：`InterSubMod/.codex`

**建議行動**：搬到 `InterSubMod/docs/archive/2026/05/20260510_codex_marker_obsolete/`（保留 audit 軌跡）

---

## §3 Tier 2 — 可能過時（建議用戶 review 後決定）

### 3.1 D2-A placeholder cycle（state/cycles/，已從 active.json 移除）

**狀態**：cycle 完成（synthetic placeholder data，非真跑），已從 `active.json` 移除但目錄仍在 `state/cycles/`

**完整 path**：`InterSubMod/state/cycles/20260507-2112-d2a-colo829-kde-rerun/`（6 件 artifacts: state / plan / precheck / pilot / generalize / evaluation）

**Plan §10 易錯點 #2 提示**：「Path A 與 Path B hand-off」期間，placeholder cycle 應視為 reference fixture

**評估**：
- ✅ Drill 2 retro 報告（`InterSubMod/docs/experiments/in_progress/2026/05/20260508_Drill2_End_to_End_Cycle_Walkthrough_01.md`）已引用此 cycle
- ⚠️ 但 cycle 仍在 active 目錄可能造成下次 `/cycle-state` dashboard 誤判

**建議行動**：搬到 `InterSubMod/state/cycles_archived/` (該目錄目前空，等用)

### 3.2 NG_KDE_Rescaled retracted 報告

**狀態**：含「已被取代 / superseded」字眼

**完整 path**：`InterSubMod/docs/experiments/in_progress/2026/04/20260421_NG_KDE_Rescaled_Multi_CN_Analysis_01.md`

**評估**：
- 內文標 superseded → 應搬到 archive 而非留在 in_progress
- in_progress/ 應只保留 active 或 validated 報告

**建議行動**：搬到 `InterSubMod/docs/archive/2026/04/20260421_NG_KDE_Rescaled_superseded.md`

### 3.3 research/v5_provenance_followup vote_dump 中間檔（5 個 .tsv.gz，~26 MB）

**狀態**：v5 audit 的中間 vote dump 檔

**完整 path**：
```
InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/
├── vote_dump_baseline_genome.tsv.gz    (大檔)
├── vote_dump_baseline_chr19.tsv.gz
├── vote_dump_v3f_chr19.tsv.gz
├── vote_dump_v5_chr19.tsv.gz
└── vote_dump_v5_genome.tsv.gz
```

**評估**：
- 大檔案，佔空間
- 是 v5 provenance audit 的研究中間檔
- 若 v5 audit 已 finalize（用戶最近 commit `388a437 docs(self-phasing): sync E5/D4 狀態 DONE` 暗示 v5 已 done）→ 可封存
- 若仍需 reproduce vote analysis → 保留

**建議行動**：**請用戶確認 v5 audit 是否已 finalize**。若已 done → 移到 `InterSubMod/research/v5_provenance_followup/_archive/`；若仍 active → 保留

---

## §4 Tier 3 — 已自封存或屬正常使用（不需動）

| 路徑 | 狀態 | 為何不動 |
|---|---|---|
| `InterSubMod/docs/plans/2025/ARCHIVED.md` | 已自我標記 ARCHIVED | 自我封存完備 |
| `InterSubMod/docs/archive/2025/` `2026/` `deep/` | 既有封存目錄 | 正常使用中 |
| `InterSubMod/state/cycles_archived/` (空 .gitkeep) | 等待 §3.1 D2-A 進駐 | 預留接收 archive |
| `InterSubMod/state/retro_cycles/` (8 cycles) | Drill 1 regression test fixture | 每次 pitfalls/evaluator 改動後跑 regression — **不可動** |
| `InterSubMod/build/bin/*` (5 binaries) | C++ 編譯產物 | gitignored；保留作執行 |
| `InterSubMod/docs/plans/2026/01/*` (10+ files) | 1 月歷史 plan | plans/ 是計畫存檔；歷史價值 |
| `InterSubMod/.claude/scheduled_tasks.lock` | scheduled task lock | 系統管理檔 |

---

## §5 預期空間節省

| Tier | 大小 | 數量 |
|---|---|---|
| Tier 1 (.agents/skills + .codex) | ~600 KB | 2 項 |
| Tier 2 D2-A cycle | ~10 KB | 1 cycle / 6 artifacts |
| Tier 2 NG_KDE retracted | ~30 KB | 1 file |
| Tier 2 vote_dump (若用戶確認封存) | ~26 MB | 5 files |
| **總計（若全執行）** | **~26.5 MB** | **9 項** |

---

## §6 Helper Script（用戶確認後手動執行）

> **未執行任何操作**。用戶 review 後決定是否跑下列 commands。

### 6.1 Tier 1 一定可以做（低風險）

```bash
# 建立 archive 目錄
mkdir -p InterSubMod/docs/archive/2026/05/20260510_agents_skills_obsolete
mkdir -p InterSubMod/docs/archive/2026/05/20260510_codex_marker_obsolete

# 6.1.a — .agents/skills/ → archive
mv InterSubMod/.agents/skills InterSubMod/docs/archive/2026/05/20260510_agents_skills_obsolete/
# 確認 .agents/ 整個目錄是否還有其他內容；若無則一起搬：
ls InterSubMod/.agents/  # 若只剩空目錄 → rmdir InterSubMod/.agents

# 6.1.b — .codex marker → archive
mv InterSubMod/.codex InterSubMod/docs/archive/2026/05/20260510_codex_marker_obsolete/
```

### 6.2 Tier 2 確認後執行

```bash
# 6.2.a — D2-A placeholder cycle → cycles_archived
mv InterSubMod/state/cycles/20260507-2112-d2a-colo829-kde-rerun \
   InterSubMod/state/cycles_archived/

# 6.2.b — NG_KDE retracted → docs/archive
mkdir -p InterSubMod/docs/archive/2026/04
mv InterSubMod/docs/experiments/in_progress/2026/04/20260421_NG_KDE_Rescaled_Multi_CN_Analysis_01.md \
   InterSubMod/docs/archive/2026/04/

# 6.2.c — v5 vote_dump（**僅當 v5 audit 已 finalize**）
mkdir -p InterSubMod/research/v5_provenance_followup/_archive
mv InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/vote_dump_*.tsv.gz \
   InterSubMod/research/v5_provenance_followup/_archive/
```

### 6.3 Commit 建議

執行 §6.1 / §6.2 後：

```bash
git add InterSubMod/docs/archive/2026/05/20260510_agents_skills_obsolete \
        InterSubMod/docs/archive/2026/05/20260510_codex_marker_obsolete

# Tier 2 額外：
# git add -u InterSubMod/state/cycles/ InterSubMod/state/cycles_archived/
# git add InterSubMod/docs/archive/2026/04/

git commit -m "$(cat <<'EOF'
chore(housekeeping): archive obsolete .agents/skills + .codex marker (cleanup proposal §6)

- .agents/skills/ (600K, 14 dirs) → docs/archive/2026/05/20260510_agents_skills_obsolete/
  Reason: superseded by .claude/skills/ (which has 11 newer skills not in .agents).
  Detected during sanity check post v1.8.

- .codex (0 bytes empty marker) → docs/archive/2026/05/20260510_codex_marker_obsolete/
  Reason: codex tool integration deprecated; cross-ref docs/archive/2026/05/20260509_codex_nested_skills_duplicate/

Cleanup proposal: InterSubMod/docs/experiments/in_progress/2026/05/20260510_cleanup_proposal_obsolete_files_01.md

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## §7 待用戶決定的問題

| # | 問題 | 預設選項 |
|---|---|---|
| 1 | Tier 1（.agents/skills + .codex）是否執行 §6.1？ | ✅ 建議直接做（兩者都明確過時） |
| 2 | D2-A cycle 是否搬到 cycles_archived？ | ✅ 建議做（cycle 已 done；不影響 reference 引用）|
| 3 | NG_KDE_Rescaled 是否搬到 docs/archive/2026/04/？ | ✅ 建議做（已標 superseded） |
| 4 | v5 vote_dump 是否封存？ | ⚠️ **需用戶確認 v5 audit 已 finalize 才執行** |
| 5 | plot_setup.py CJK fix 是否 commit？ | ✅ 建議做（fix infra small follow-up） |

---

## §8 References

- v1.8 retro 報告：`InterSubMod/docs/experiments/in_progress/2026/05/20260510_v1.8_implementation_retro_01.md`
- 既有封存目錄結構：`InterSubMod/docs/archive/`
- Drill 2 retro（引用 D2-A cycle）：`InterSubMod/docs/experiments/in_progress/2026/05/20260508_Drill2_End_to_End_Cycle_Walkthrough_01.md`
- 已封存 codex 案例（cross-ref）：`InterSubMod/docs/archive/2026/05/20260509_codex_nested_skills_duplicate/`
