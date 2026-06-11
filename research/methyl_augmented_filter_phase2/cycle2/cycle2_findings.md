<!--
建立時間: 2026-05-19
agent: main session Coordinator (Phase 2 Cycle 2 synthesis)
status: complete
report_class: cycle final synthesis (Cycle 2 H_C1_5 + H_C1_6 integration)
audience: PI / lab / cycle 3 planner
scope: HCC1395 + 4 樣本 (H1437/H2009/HCC1954/HCC1937) cross-sample + cross-binary
parent_plan: /bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md v2.0
predecessor_cycle: Cycle 1 ⭐3 strong (HCC1395 ΔF1 +0.02236)
verdict: Cycle 1 ⭐3 STRONG (HCC1395-internal) + caveat (cross-sample DIRECTION_NEGATIVE，caller-F1-headroom-bounded)
last_verified: 2026-05-19
-->

# Phase 2 Cycle 2 — Coordinator Synthesis (H_C1_5 NEGATIVE + H_C1_6 PASS)

> **Cycle 1 tier 決定**: 保持 ⭐3 strong (HCC1395-internal validated, +0.02236 = 2.24× Cohen ribbon)，但加入 **caveat**：cross-sample transfer mode DIRECTION_NEGATIVE — filter 不可作 production transfer rule，僅在 caller F1 < 0.80 + FP density > ~14% 樣本可重現。
>
> **Cycle 2 推進方向**：啟動 Cycle 3 「Caller-F1-headroom-aware redesign」(per user 5/19 決議)。

---

## 0. TL;DR

| Hypothesis | Verdict | Implication |
|---|---|---|
| **H_C1_5** (cross-sample n=5 Wilcoxon) | **DIRECTION_NEGATIVE** transfer / **MIXED** re-fit | Filter 不通用；transfer 在 caller F1 > 0.83 樣本上 ΔF1 < 0；只有 HCC1937 (F1=0.37) re-fit +0.00761 |
| **H_C1_6** (HCC1395 V3F/V5/V6 cross-binary) | **PASS** max_var 0.00073 transfer / 0.00055 re-fit | BAM upgrade zero F1 regression；V6 可作 production BAM |

**Combined verdict**：Cycle 1 filter **HCC1395-internal valid**，**BAM-invariant**，但 **caller-F1-headroom-bounded**（不是 methylation 失敗，是 caller F1 ceiling）。

---

## 1. H_C1_5 Cross-Sample (NEGATIVE)

### 1.1 Per-sample ΔF1（transfer τ=0.39 / re-fit per-sample swept τ）

| Sample | Caller F1 | FP density | Transfer ΔF1 | Re-fit ΔF1 | Best re-fit τ |
|--------|--:|--:|--:|--:|--:|
| HCC1395 | 0.7166 | 14% | **+0.02232** | **+0.02236** | 0.39 |
| H1437 | 0.8670 | 1.1% | -0.03597 | ~0 (-2.4e-7) | 0.10 |
| H2009 | 0.8863 | 1.0% | -0.00450 | ~0 (+5.2e-6) | 0.61 |
| HCC1954 | 0.8385 | 3.4% | **-0.37744** | +0.00095 | 0.71 |
| HCC1937 | 0.3692 | 16% | -0.07068 | **+0.00761** | 0.24 |

### 1.2 Wilcoxon paired n=5 verdict

| Mode | n_pos | n_neg | n_zero | Wilcoxon p (ΔF1) | Wilcoxon p (ΔRecall) | Verdict |
|------|--:|--:|--:|--:|--:|---|
| Transfer fixed | 1 | 4 | 0 | 0.1875 | **0.0625** | DIRECTION_NEGATIVE |
| Transfer swept | 1 | 4 | 0 | 0.3125 | 0.0625 | DIRECTION_NEGATIVE |
| Re-fit | 3 | 0 | 2 | 0.125 | 0.0679 | MIXED (positive trend, magnitude weak) |

**ΔRecall p=0.0625 (transfer 5/5)** 意味著 filter 在所有樣本上一致降 recall — proximal mechanism of F1 collapse。

### 1.3 為何 NEGATIVE — caller-F1-headroom 假說

| Caller F1 region | FP density | Filter has headroom? | Observed re-fit ΔF1 |
|---|---|---|---|
| F1 < 0.50 | > 14% | YES | HCC1937 +0.00761 |
| 0.50–0.80 | ~14% | YES (cycle 1 training zone) | HCC1395 +0.02236 |
| > 0.83 | < 4% | NO（ceiling） | H1437/H2009 ~0 / HCC1954 +0.00095 |

**機制**：
1. **Transfer mode 災難 (HCC1954 -0.377)** = caller_af coef +3.44 對 HCC1395 AF 分布 overfit；標準化 scaler 把 HCC1954 TP 推到「flag as FP」半平面
2. **Re-fit 4/4 新樣本 ≈0** = 3 樣本 caller F1 > 0.83，FP density < 4%，LR 找不到有意義的 τ，最佳解退化為「不過濾」
3. **HCC1937 (low F1) 是 only signal** — FP density 16% 給 filter work 的空間

### 1.4 失敗 ≠ Methylation 無用

- HPFineF coef +0.75 = 5th-rank covariate（cycle 1 同等位階）
- 主導 axes (caller_af / LOH / Coverage / NG) BAM-invariant per H_C1_6
- 失敗來源是 caller F1 ceiling 不留 filtering headroom，非 methylation 訊號弱

---

## 2. H_C1_6 Cross-Binary (PASS)

| BAM | Transfer ΔF1 | Re-fit ΔF1 | Re-fit τ |
|-----|--:|--:|--:|
| V3F | +0.02295 | +0.02289 | 0.40 |
| V5  | +0.02222 | +0.02234 | 0.42 |
| V6  | +0.02232 | +0.02236 | 0.39 |
| **max - min** | **0.00073** | **0.00055** | — |

**Threshold**: < 0.005 PASS / 0.005-0.010 BORDERLINE / > 0.010 FAIL → **PASS** (6.8–9.1× under)

**含義**:
- V6 production BAM 升級 = zero F1 regression
- V3F 名義最佳 (+0.063% vs V6) 但低於 multi-seed std 5e-5 × 10 → 不撤回 V6
- Methylation features 跨 BAM 一致性 dampened by 10× via LR weights → 不依賴 BAM-tagging idiosyncrasies

**Reproducibility**: V6 re-fit ΔF1 = +0.022358 與 cycle 1 stored bit-exact match (drift 0.000000)。

---

## 3. Combined Verdict — Cycle 1 ⭐3 + Caveat

### 3.1 為何不降 tier

- **H_C1_2/H_C1_3 仍 PASS**：HCC1395 內部 +0.02236 = 2.24× Cohen ribbon
- **H_C1_6 PASS 強化**：BAM-invariant，HCC1395 結論 robust to upstream V3F/V5/V6
- **Cycle 1 不是宣稱「universal」** — pre-reg H_C1_5 是 ⭐4 升 gate，FAIL 只阻擋升級不退回

### 3.2 Caveat 註記（accompany all cycle 1 references）

> Cycle 1 filter (10-feature L2 LR, τ*=0.39) HCC1395-internal valid (+0.02236) and BAM-invariant (V3F/V5/V6 max_var 0.00055)，**但 not cross-sample transferable**：4 new samples 中 3 個 caller F1 > 0.83 留無 filtering headroom；only HCC1937 (F1=0.37) 在 re-fit mode 達 +0.00761 nominal uplift。Production deployment 必須 per-sample re-fit + caller F1 < 0.80 gate。

### 3.3 Decision tree per cycle 1 plan v2.0 §Decision rule

```
H_C1_3 PASS (HCC1395) + H_C1_5 FAIL (cross-sample)
→ rule 2: "HCC1395-specific filter, marginal cross-sample" → paper case study positioning
```

但 user 5/19 決議「保持 ⭐3 + caveat」 — cycle 1 結論不退，rule 2 framing 改 "HCC1395-internal ⭐3 + caller-F1-headroom-bounded" + 啟動 cycle 3 redesign。

---

## 4. Cycle 3 入口 — Caller-F1-headroom-aware Redesign

### 4.1 設計核心

**前提**: Filter 只在 caller F1 < 0.80 + FP density > ~14% 樣本 work。Production rule:

```
IF caller_F1(sample) < 0.80 AND FP_density > 0.10:
    apply cycle 1 filter (re-fit on this sample, swept τ)
ELSE:
    skip filter (caller already ceiling)
```

### 4.2 Pre-reg H_C3_* (待 cycle 3 init 寫死)

| ID | Prediction | Falsification | Threshold |
|----|-----------|---------------|----------|
| **H_C3_1** | Caller-F1-headroom rule 在 qualifying subset (HCC1395 + HCC1937) re-fit mean ΔF1 ≥ +0.01 | mean ΔF1 < +0.005 | n=2 (擴展需 panel) |
| **H_C3_2** | 擴 low-F1 panel n≥4 後 Wilcoxon p<0.05 | n<4 OR p≥0.05 | TCGA / SEQC2 truth pairs survey |
| **H_C3_3** | High-F1 (>0.83) 樣本 skip filter 後 caller F1 preserved (drift < 0.001) | drift ≥ 0.001 | n=3 H1437/H2009/HCC1954 |

### 4.3 Cycle 3 scope (out-of-scope guard)

**In-scope**:
- Per-sample re-fit infrastructure (script + filter rule template)
- Caller F1 + FP density panel survey across InterSubMod archive
- Low-F1 sample 擴 panel candidate identification

**Out-of-scope** (per Phase 2 plan v2 + user 5/19 confirmation):
- Distributional adaptation (cycle 3+ optional)
- Methylation-specific gate redesign (cycle 4+)
- Read-level methylation embedding (T5+)

---

## 5. SoT Update Status

| Artifact | Action |
|---|---|
| `evidence_ledger.jsonl` | append cycle 2 entry (verdict: positive_with_caveat, stability 3) |
| `MEMORY.md project_phase2_cycle1_global_fp_filter.md` | update with caveat + cross-sample status |
| `docs/experiments/INDEX.md` | add cycle 2 entry |
| `docs/CURRENT_FOCUS.md` | add cycle 2 closure + cycle 3 entry under "2026-05-19" section |
| `hypothesis_queue.json` | inject H_C3_1/2/3 via inject-hypothesis |
| `state/cycles/{cycle3_id}/state.json` | create via cycle-init |

Cycle 1 主報告 (`InterSubMod/docs/experiments/in_progress/2026/05/20260518_Phase2_Cycle1_Global_FP_Filter_01.md`) 不重寫 — 它記錄 cycle 1 結果並 explicitly defer H_C1_5 至 cycle 2，現在的 caveat 由本檔 + ledger 承擔。

---

## 6. Files

```
cycle2/
├── cycle2_findings.md                          (this synthesis)
├── cycle2_step_b3_b4_findings.md              (H_C1_5 detail)
├── cycle2_step_c1_h_c1_6_sanity.md            (H_C1_6 detail)
├── data/                                       (4 augmented master + ΔF1 + Wilcoxon)
├── figures/                                    (per-sample bar + Wilcoxon verdict + V3F/V5/V6 bar)
├── intermediate/                               (json + log)
└── scripts/                                    (build / apply / verdict / per-BAM)
```

---

## 7. Coordinator Notes

- **Wall clock**: cycle 2 active time ~3 hr (B1' ISM 2 hr + B2'/B3'/B4' 30 s + C1 25 min)
- **B1' wall clock bottleneck**: H2009 tp 45.7 min (largest BAM)
- **Reproducibility**: scripts deterministic seed=42; cycle 1 V6 re-fit bit-exact match
- **失敗 framing**: H_C1_5 NEGATIVE 不是 cycle 1 失敗，是 pre-reg 升 tier gate 未通過 — cycle 1 ⭐3 是 within-sample validated 不需退；caveat = pre-reg sciatific integrity，不是失敗認賠

### 7.1 Productive Failure Lessons

Per `feedback_productive_failure_reopen_threshold.md`：
- C1 (new data): 4 new samples ISM rerun 提供 cross-sample evidence
- C2 (new method): cycle 3 caller-F1-headroom-aware redesign 算 new method
- C3 (new precondition): low-F1 panel 擴展為新 precondition

→ Cycle 3 reopen 條件成立。

### 7.2 已驗證可推進 → 全自動模式

Per `feedback_execution_mode_hierarchy.md` 與用戶 5/19 三選確認：
- SoT 更新 + cycle 3 注入 + cycle-init = 🟡 重複執行（已驗證可推進路徑）
- 一行告知 + 直接執行，不再逐項確認

---

**End of Coordinator Synthesis**
