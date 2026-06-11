<!--
build_date: 2026-05-22
agent: weekly-report handoff D (auto, default)
status: in_progress
report_class: next-week-plan
parent: master_draft.md (same dir)
last_verified: 2026-05-22
-->

# 下週計畫 (2026-05-23 ~ 2026-05-29) — V6 後續 + Pivot 決策

> 本檔配合 master_draft.md 一同產出（weekly-report handoff option D default）；待 user 在週一上午 PI 1-on-1 後依 PI Q3 回答決定 pivot 主軸 (a/b/c)。

## 0. Hard Gate Pending (本週交接)

| Item | Status | Owner |
|---|---|---|
| Send V6 sign-off email to PI (`20260521_PI_V6_signoff_email_draft_5goal.md`) | 🔴 user 親自 copy 到 mail client | user |
| `git push origin fix/pon-only-phasing v6-prod-20260520` | 🔴 等 user 明示授權 | user |

## 1. PI 回答 → Pivot 決策分叉 (Plan A / B / C)

### Plan A：read-level pivot (default 推薦)

| Step | Skill | Cost |
|---|---|---|
| 1 | `/init-research phase_block_3d` (建立 `research/phase_block_3d/` scaffolding) | 30 min |
| 2 | `/inject-hypothesis` H_PB3D_001 (read-level epigenotype + methyl-phasing) | 15 min |
| 3 | `/cycle-init` for cycle 1 (read-level pilot, HCC1395 chr19 subset) | 30 min |
| 4 | Pilot: 跑 read-level pipeline (~2 hr ISM rerun with read-level annotation) | 2 hr |
| 5 | `/feature-layered-observation` Step 1-6 (read-level features × LOH×AF×CN 32-cell) | 1 day |

### Plan B：caller-F1-headroom redesign (cycle 5+)

| Step | Skill | Cost |
|---|---|---|
| 1 | `/cycle-init` for cycle 5 (caller-F1-headroom-aware) | 30 min |
| 2 | `/research-loop` P1 PLAN (plan v3 — caller F1 < 0.80 + FP density > 0.10 gate per-sample re-fit) | 1 hr |
| 3 | low-F1 panel sample identify (search ONT WGS public cohort) | 1 day |
| 4 | apply Plan B + LOSO | 2 day |

### Plan C：low-F1 panel sample expansion

| Step | Skill | Cost |
|---|---|---|
| 1 | identify HCC1937-like (caller F1 < 0.5) sample 3-5 個 | 1 day |
| 2 | ISM pipeline rerun on new samples (each ~3 hr) | 1 week |
| 3 | apply A4 multi-algo on extended cohort | 1 day |

### 推薦並行: A + C

Plan A 為主軸 (read-level pivot, breakthrough potential); Plan C 為強 negative control (低 F1 panel 是否提供 generalize signal)。Plan B deprioritize 因 cycle 3 已試類似邏輯仍 sample-level circular。

## 2. Operational Maintenance

- weekly-report C4 confirmation：user review master_draft.md + auto_decisions log; 決定 handoff (A/B/C/D)
- PI 4-goal HTML build：`/html-report-build standalone` mode for `docs/presentations/in_progress/20260522_PI_V6_signoff_4goal/` (use A0 audit verified source list)
- `evidence_ledger.jsonl` append entry #52 (Phase A 7-task + Phase B weekly-report)
- MEMORY.md update：`project_phase2_a_completeness_audit.md` 新增
- CURRENT_FOCUS.md update：5/22 區塊加 V6 sign-off complete + Phase A 7-task complete + Phase B 週報 produced

## 3. Risk

| Risk | Likelihood | Mitigation |
|---|---|---|
| PI 對 LR LOSO 失敗結論不同意 (e.g. "你還沒試 X algo") | low | A4-Ext 6 algo inventory + A4 4 algo 100-fold = 10 algo space 已覆蓋 |
| PI 要求 V3F 重啟 Goal 1 dedicated rerun | medium | ~2 hr ISM rerun 可接受；preserve V3F BAM 不刪 |
| Plan A read-level pilot 仍 NEGATIVE | medium | Plan C 並行作 backup |
| pivot 拖到 2 週才有結果 | medium | 設 cycle 5 strict probe-first (per Cynefin Complex 域) — 2 week budget cap |
