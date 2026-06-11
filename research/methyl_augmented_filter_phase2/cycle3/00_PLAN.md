<!--
建立時間: 2026-05-19
agent: main session Coordinator
status: P0_REGISTER → P1_PLAN draft (待 research-loop refine)
report_class: cycle plan
audience: cycle 3 executor / coordinator / panel-search agent
scope: methyl_augmented_filter_phase2 cycle 3
parent_cycle: state/cycles/20260519-0031-cycle3-caller-f1-headroom/state.json
predecessor_cycle: 20260519_phase2_cycle2_cross_sample_negative (cycle 2 ledger entry)
last_verified: 2026-05-19
-->

# Phase 2 Cycle 3 — Caller-F1-headroom-aware Filter Redesign

> **Cycle ID**: `20260519-0031-cycle3-caller-f1-headroom`
> **Primary hypothesis**: H016 (H_C3_1) + companion H017 (H_C3_2) + H018 (H_C3_3)
> **Predecessor**: Cycle 2 H_C1_5 DIRECTION_NEGATIVE + H_C1_6 PASS — filter caller-F1-headroom-bounded
> **Goal**: 設計 production-deployable per-sample filter rule that gates on caller F1 + FP density，並透過 n≥4 low-F1 panel 達 Wilcoxon p<0.05

---

## 1. Context

### 1.1 Cycle 2 結論 (caller-F1-headroom-bounded)

| Sample | Caller F1 | FP density | Transfer ΔF1 | Re-fit ΔF1 |
|--------|--:|--:|--:|--:|
| HCC1395 | 0.7166 | 14% | **+0.02232** | **+0.02236** |
| H1437 | 0.8670 | 1.1% | -0.03597 | ~0 |
| H2009 | 0.8863 | 1.0% | -0.00450 | ~0 |
| HCC1954 | 0.8385 | 3.4% | **-0.37744** | +0.00095 |
| HCC1937 | 0.3692 | 16% | -0.07068 | **+0.00761** |

**Pattern**: Filter only works when caller F1 < 0.80 AND FP density > ~14%. 3/4 新樣本 caller F1 > 0.83 留無 filtering headroom；HCC1937 為 only generalize-positive case。

### 1.2 與已關閉方向差異化 (Research Direction Guard)

| 已關閉 | 本 cycle 差異化 |
|---|---|
| [NEGATIVE] Methyl Filter Pilot v1.0 cell-gated (2026-05-18) | Cycle 3 不重做 cell-gated；continue cycle 1 global LR 框架 |
| [DIRECTION_NEGATIVE] Cycle 2 cross-sample (2026-05-19) | **明確**新增 caller-F1 + FP density gate，**不再宣稱 universal transfer** |
| [NO-GO] TO Germline FP G1-G7 (AUC<0.64) | Cycle 3 仍在 paired-pileup + multi-axis LR；不是 TO germline FP framework |
| [NO-GO] Pure methylation Beyond-AUC (2026-04-09) | Methyl 仍是 5th-rank covariate；filter 設計重點是 gating rule，非 methylation 強化 |

---

## 2. Pre-Registration (Hard Gate per scientific-rigor §7.1)

| ID | Hypothesis | Prediction | Falsification | Decision threshold |
|----|-----------|-----------|---------------|-------------------|
| **H_C3_1** (H016) | Caller-F1-headroom rule | n=2 qualifying (HCC1395 + HCC1937) re-fit mean ΔF1 ≥ +0.01 | mean ΔF1 < +0.005 | Pre-reg 2 樣本 strict mean |
| **H_C3_2** (H017) | Panel expansion | 擴 n≥4 low-F1 panel 後 Wilcoxon paired p < 0.05 | n<4 OR p ≥ 0.05 | Panel survey + ISM rerun |
| **H_C3_3** (H018) | High-F1 skip preserves caller F1 | n=3 high-F1 (H1437/H2009/HCC1954) skip filter, caller F1 \|drift\| < 0.001 | any sample \|drift\| ≥ 0.001 | Sanity check |

**NO-GO 條件**:
- H_C3_1 PASS + H_C3_2 FAIL (n<4 panel) → cycle 3 結論 "qualifying subset only POSITIVE，⭐3 strict subset"，⭐4 升級擱置
- H_C3_1 FAIL → 整 cycle NEGATIVE，cycle 4 considering distributional adaptation alternative
- H_C3_3 FAIL → mechanism contradiction，需 root-cause why high-F1 skip 仍引入 drift

---

## 3. Plan Structure — 3 Step

### Step 1 — Caller-F1-headroom rule formalization (~30 min)

**Goal**: 把 "IF caller F1 < 0.80 AND FP density > 0.10 apply per-sample re-fit ELSE skip" 寫成 deployable rule + JSON config。

| Activity | Detail |
|---|---|
| 1a | 寫 `cycle3_caller_f1_gate.json` — gate threshold + per-sample re-fit fallback config |
| 1b | 寫 `scripts/apply_gated_filter.py` — input: 5 sample master TSV + caller F1 + FP density, output: gated ΔF1 |
| 1c | Apply on cycle 2 5 samples reproduce → expected H_C3_1 (HCC1395+HCC1937) PASS / H_C3_3 (3 high-F1) PASS |
| 1d | 輸出 `cycle3_step1_findings.md` |

**Expected**: HCC1395 +0.02236 + HCC1937 +0.00761 mean = +0.0150 → PASS H_C3_1 (1.5× threshold) | 3 high-F1 skip → drift 0 (trivial PASS H_C3_3)

### Step 2 — Low-F1 panel survey (~2-3 hr)

**Goal**: 找 n ≥ 4 low-F1 (caller F1 < 0.80) samples with truth set + V6 BAM 已 tagged + ISM 可跑。

| Activity | Detail |
|---|---|
| 2a | Inventory all InterSubMod archive samples with caller F1 < 0.80 + truth pair |
| 2b | 對 candidate list 驗 V6 BAM 是否 tagged (longphase-to-mod V6) |
| 2c | 如 candidate <4: 列 COLO829 truth set status (chenhan112 權限)、SEQC2 derivatives、其他 |
| 2d | 輸出 `cycle3_step2_panel_inventory.md` |

**Risk**: HCC1937 是唯一已知 low-F1 sample；如無 ≥3 新 low-F1 candidate → H_C3_2 stuck (panel 不夠)。

**降級 plan**:
- 若 n=2 (HCC1395+HCC1937 only): H_C3_2 標 "panel insufficient, n<4"，cycle 3 結論限 qualifying subset
- 若 n=3-4: 跑 Wilcoxon exact，min p=0.0625 (n=4 全同方向) 至 0.0312 (n=5 全同方向)

### Step 3 — Cycle 3 verdict synthesis (~1 hr)

| Activity | Detail |
|---|---|
| 3a | 整合 Step 1+2 → H_C3_1/2/3 verdicts |
| 3b | Cycle 3 tier 判定（per run-evaluator P5） |
| 3c | Cycle 3 結論 → Cycle 4 入口（若 H_C3_2 FAIL panel）或 production deployment (若 PASS) |
| 3d | 輸出 `cycle3_findings.md` + ledger entry + INDEX/CURRENT_FOCUS update |

---

## 4. Expected timeline

| Step | 估時 | 依賴 |
|------|------|------|
| Step 1 Gate rule formalization | 30 min | cycle 2 master TSV (existing) |
| Step 2 Panel survey | 2-3 hr | InterSubMod archive metadata |
| Step 2 ISM rerun (if new sample found) | +2-4 hr | V6 BAM existence check |
| Step 3 Verdict synthesis | 1 hr | Step 1 + 2 |
| **Total active** | **4-8 hr** | (depends on panel hit) |

---

## 5. Out-of-scope (避免 scope creep)

- ❌ Distributional adaptation (robust scaler / IQR / per-sample mean-shift) → cycle 4+ if cycle 3 stuck
- ❌ Methylation-specific gate redesign (methyl 5th-rank, 不是 root cause) → cycle 5+
- ❌ Read-level methylation embedding (T5+ multi-month)
- ❌ COLO829 truth set 取得（chenhan112 權限 pending，不阻塞 cycle 3）

---

## 6. Verification (end-to-end)

### 6.1 Step verification gates

1. Step 1 通過: Gate JSON 寫好 + 5 樣本 reproduce ΔF1 + H_C3_1/H_C3_3 verdicts 明確
2. Step 2 通過: Panel inventory 列 ≥1 candidate (或明確 n=2 only) + 各 candidate ISM status
3. Step 3 通過: Cycle 3 ledger entry append + state.json phase P6_COMMIT + tier 寫入

### 6.2 命令模板

```bash
# Step 1
python3 research/methyl_augmented_filter_phase2/cycle3/scripts/apply_gated_filter.py \
    --gate-config cycle3/cycle3_caller_f1_gate.json \
    --samples HCC1395,H1437,H2009,HCC1954,HCC1937 \
    --output cycle3/data/cycle3_step1_gated_delta_f1.tsv

# Step 2 (panel survey 主要是手動 + inventory script)
python3 scripts/inventory/list_caller_f1_lt_080.py \
    --archive-root /big7_disk/liaoyoyo2001/InterSubMod \
    --output cycle3/data/cycle3_step2_panel_candidates.tsv
```

---

## 7. Multi-agent fan-out 設計（optional）

| Agent | 任務 | Mode |
|-------|------|------|
| Agent 1 | Step 1 gate rule formalization | foreground |
| Agent 2 | Step 2 panel survey | background (independent) |
| Agent 3 | (if Step 2 found candidate) ISM rerun | background |
| Coordinator (main) | Step 3 synthesis | — |

---

## 8. Decision rule (post-cycle 3)

1. **H_C3_1 PASS + H_C3_2 PASS (n≥4 panel)** → ⭐4 cross-sample validated for qualifying subset
2. **H_C3_1 PASS + H_C3_2 FAIL (panel <4)** → ⭐3 strict subset only，filter limited to caller F1 < 0.80
3. **H_C3_1 FAIL** → cycle 4 distributional adaptation or methyl filter direction abandon
4. **H_C3_3 FAIL** (gate rule 漏 high-F1 sample) → root-cause investigation 必要

---

## 9. SoT links

- Cycle state: `InterSubMod/state/cycles/20260519-0031-cycle3-caller-f1-headroom/state.json`
- Hypothesis_queue entries: H016 (priority 88) / H017 (86) / H018 (84) all queued
- Cycle 2 predecessor ledger: `cycle_id: 20260519_phase2_cycle2_cross_sample_negative`
- Cycle 1 deployable filter: `cycle1/cycle1_track_a_filter.json` (reuse for gated mode)
- CURRENT_FOCUS: 2026-05-19 section

---

## 10. Risk + Mitigation

| Risk | Likelihood | Mitigation |
|------|---|---|
| Panel survey yields n=2 only | HIGH (HCC1937 was unique) | Pre-reg n<4 → 限 qualifying subset 結論 ⭐3，不退 NEGATIVE |
| ISM rerun on new candidate fails (BAM issues) | MEDIUM | Step 2c 列 ≥3 backup candidates |
| User pivot 到 phase_block_3d / thread_d before cycle 3 done | MEDIUM | Plan keep cycle 3 ≤8 hr，不重於四軌平行 |
| Cycle 3 結論 "panel insufficient" 不滿足升 tier | LOW | 接受 ⭐3 strict subset 為 final，方向 archive |

---

**End of Cycle 3 Plan v0.1** — 待 `/research-loop` 進入 P1 PLAN 後寫入 `plan.json` (OSF-style 4-section)。
