<!--
建立時間: 2026-05-19
類型: NEGATIVE Postmortem (SRE blameless 5-段)
觸發: plan v2 R-MENTAL-DRIFT 紀律 + 4-agent multi-agent audit (5/19)
Parent cycle: 20260519_phase2_cycle2_cross_sample_negative
規範來源: ~/.claude/plans/tender-pondering-blossom.md §6 Risk Items
-->

# Cycle 2 H_C1_5 DIRECTION_NEGATIVE Postmortem

> **Blameless 5-段格式（SRE 範式）**
> **目的**：plan v2 R-MENTAL-DRIFT 規範要求 — 每次 cycle NEGATIVE 後 48hr cooling-off + 寫 postmortem 才能 reactivate
> **5/19 4-agent audit 揭露**：本 postmortem 為紀律補正

---

## 1. Summary（事實摘要，無歸咎）

2026-05-19 完成 methyl_filter_phase2 cycle 2 cross-sample validation。Cycle 1 HCC1395-trained filter (10-feature L2 LR, τ*=0.39) 在 4 新樣本（H1437/H2009/HCC1954/HCC1937）V6 ISM transfer + re-fit 兩 mode 驗證：

- **Transfer fixed τ=0.39**：1+/4− (Wilcoxon p=0.1875)；ΔRecall p=0.0625 達 n=5 exact 最小值（5/5 ΔRecall<0）
- **Re-fit per-sample**：3+/0−/2≈0 MIXED (p=0.125)；唯一非邊緣 POSITIVE 為 HCC1937 (caller F1=0.37) re-fit +0.00761

**Verdict**：H_C1_5 DIRECTION_NEGATIVE — 非 underpowered，而是 **mechanism-bounded NEGATIVE**（caller-F1-headroom-bounded）。

**Caveat**：cycle 1 HCC1395-internal ⭐3 結果（ΔF1=+0.02236）不受影響；H_C1_6 V3F/V5/V6 cross-binary PASS (max_var transfer 0.00073 / re-fit 0.00055；V6 re-fit drift=0 bit-exact)。

---

## 2. Impact（影響範圍）

### Direct Impact
- methyl_filter_phase2 cycle 1 ⭐3 **不退**（HCC1395-internal valid）
- Cross-sample F1 提升路徑 **bounded but not closed** — caller F1<0.80 + FP density>0.10 仍有 conditional uplift potential
- G5 (F1 提升) plan v2 path (a) NOT activated，NOT closed（在 plan v2 R-CONFIRMATION-LOOP 邊界灰色地帶）

### Indirect Impact
- Plan v2 三 papers (A framework / B clone / C ASM+two-hit) **不直接受影響**（A framework 主軸是 V6 production + thread_d + phase_block_3d，不依賴 methyl filter）
- H_C1_6 PASS **強化 V6 production tag** — BAM upgrade zero F1 regression 已驗證
- HCC1937 +0.00761 valuable per-sample signal — 但 effect size < Cohen small (0.2)，需更多 low-F1 樣本擴充

### No Impact
- thread_d_paper 主軸 (TP-enriched phasing signatures)
- selfphasing_v6_production 4-day workflow
- phase_block_3d (5/23 init)
- Plan v2 §3 Phase A/B/C timeline

---

## 3. Root Cause（根因分析，技術而非人）

### 3.1 Primary mechanism: caller-F1-headroom-bounded

| Caller F1 region | n samples | FP density | Filtering room | Re-fit ΔF1 |
|------------------|-----------|------------|----------------|------------|
| < 0.50 (HCC1937) | 1 | 16% | YES | +0.00761 |
| 0.50–0.80 (HCC1395) | 1 | 14% | YES | +0.02236 |
| > 0.83 (H1437/H2009/HCC1954) | 3 | 1.1–3.4% | **NO** | ~0 |

**機制鏈**：高 caller F1 → 少 FP candidates → 沒有 filtering room → filter τ→0.10 (keep everything) → ΔF1≈0。

### 3.2 Secondary: HCC1954 catastrophic transfer (-0.377)

- Cycle 1 caller_af coefficient +3.44 (top rank)
- HCC1395 TP caller_af 中心 ~0.45；HCC1954 TP 中心 ~0.35
- Cycle 1 標準化 scaler 應用到 HCC1954 → TP 集中到 "filter out" 半平面 → 災難性誤刪

**根因**：filter coef HCC1395 AF 分布 overfit，非 methylation 機制失敗。

### 3.3 Tertiary: methylation 5th-rank in cycle 1 LR

- cycle 1 coef: caller_af (+3.44) > LOH_inner (+1.46) > Cov (+1.27) > NG (+1.07) > HPFineF (+0.75) > [4 其他 methyl]
- "Multi-axis filter incl. methylation" reframing 在 cycle 2 cross-sample 仍 hold — methyl 是 weak axis 但非完全無效

---

## 4. Detection / Response（偵測與回應時序）

| Date / Time | Event |
|-------------|-------|
| 2026-05-18 17:00 | v1.0 cycle marginal +0.00242 |
| 2026-05-18 18:30 | Phase 2 Cycle 1 strong ΔF1=+0.02236（9.24×） |
| 2026-05-19 00:30 | Cycle 2 H_C1_5 cross-sample DIRECTION_NEGATIVE + H_C1_6 cross-binary PASS |
| 2026-05-19 (早) | Cycle 3 redesign proposed in `cycle2_findings.md §4` |
| 2026-05-19 (本 audit) | **4-agent multi-agent audit (Explore × 4) 揭露 critical issues**：<br>• R-MENTAL-DRIFT 0% compliance（無 cooling-off, 無 postmortem）<br>• cycle_id `_negative` 命名 narrative bias<br>• Provenance 缺漏（binary_commit, dataset_id）<br>• Cycle 3 H_C3_2 low-F1 panel 樣本名單不存在 |
| 2026-05-19 (本 postmortem) | **Cycle 3 啟動暫停 + 補 P0 3 項** |

### Audit 揭露重要 finding

**4-agent 平行 audit 在 cycle 3 即將啟動前介入** — 若未介入，cycle 3 會在 cycle 2 NEGATIVE 不到 24 hr 內啟動，違反 plan v2 自定 48hr cooling-off 紀律。這是 plan v2 R-MENTAL-DRIFT 紀律的**第一次真實測試**：紀律是 commitment 還是表演？

**Response**：用戶 explicit 選擇「暫停 cycle 3，補 P0 3 項」→ 紀律保持 commitment。

---

## 5. Lessons & Action Items

### 已執行（本 postmortem 階段）

- ✅ Postmortem .md 寫入（本檔）
- ✅ evidence_ledger 加 postmortem entry (cycle_id `20260519_phase2_cycle2_negative_postmortem`)
- ✅ CURRENT_FOCUS amendment（cycle 3 延後 + audit findings 紀錄）
- ✅ 48hr cooling-off period: 2026-05-19 → 2026-05-21

### Cycle 3 啟動前 must-fix（解除暫停的前置條件）

| Pri | 項目 | 對應 audit issue |
|-----|------|-----------------|
| 🔴 P0 | H_C3_2 low-F1 panel n≥4 具體樣本名單（含 BAM 可用驗證）| Agent C Q3 |
| 🟠 P1 | cycle_id 重命名規範：`_validation` 不 `_negative`（verdict 在 entry 不在 name） | Agent D Q1 |
| 🟠 P1 | ledger entry 補 binary_commit hash / dataset_id / pre-reg link | Agent D Q5 |
| 🟡 P2 | H_C3_1 target 計算澄清：cycle 2 re-fit mean (+0.015) vs cycle 3 重新 re-fit | Agent C Q2 |
| 🟡 P2 | HCC1954 overfit + caller-headroom-bounded mechanism 視覺化 figures（AF histogram / coef inspection / scatter plot）| Agent A Q3/Q4 |

### Plan v2 紀律加強（avoid future drift）

- **R-MENTAL-DRIFT enforcement**：每次 cycle NEGATIVE 自動 trigger postmortem template + 48hr cooling-off ledger entry（建議 hook 自動化）
- **cycle naming convention**：`{date}_{project}_{phase}_{cycle_n}_{scope}_validation` — verdict 不入 cycle_id
- **ledger schema 強化**：必填欄位 `binary_commit` + `dataset_id` + `pre_registration_link`
- **Pre-cycle-init audit**：cycle 啟動前自動跑 multi-agent rationality check（hook）

### 對 plan v2 R-CONFIRMATION-LOOP 條款的反思

Plan v2 §6 R-CONFIRMATION-LOOP 寫「cycle 1 結果決定（NEGATIVE → G5 真正關閉）」。但實際情況更複雜：

- Cycle 1 strong (NOT NEGATIVE) → path (a) 通過 ✓
- Cycle 2 NEGATIVE cross-sample → 不是「cycle 1 NEGATIVE」 → R-CONFIRMATION-LOOP 條款**未明確涵蓋此情況**
- Cycle 3 「caller-F1-headroom-aware」屬於 **productive failure escalation** 還是 **confirmation loop 修補**？

**判定**：因 cycle 2 確認**機制邊界**（caller-F1-headroom）而非「重複跑找 marginal」，屬 productive failure escalation。但需設**硬性 stop rule**：

> Plan v2 amendment proposal: 若 cycle 3 H_C3_1 target ≥+0.01 未達 + cycle 3-4 仍 cross-sample NEGATIVE → methyl_filter_phase2 正式關閉，G5 path (a) 確定走向 close。

---

## 6. References

- Cycle 2 main findings: `InterSubMod/research/methyl_augmented_filter_phase2/cycle2/cycle2_findings.md`
- Cycle 2 H_C1_5 detail: `cycle2_step_b3_b4_findings.md`
- Cycle 2 H_C1_6 detail: `cycle2_step_c1_h_c1_6_sanity.md`
- Cycle 1 findings: `InterSubMod/research/methyl_augmented_filter_phase2/cycle1/cycle1_findings.md`
- Plan v2: `~/.claude/plans/tender-pondering-blossom.md` (§6 Risk Items)
- 4-agent audit conversation: 5/19 session log

---

**Postmortem 結論**：cycle 2 NEGATIVE 是 plan v2 紀律的第一次真實測試。經 4-agent audit + 用戶 explicit 暫停決策後，紀律保持 commitment。Cycle 3 啟動延後到 2026-05-21 後（48hr cooling-off 完成 + P0 must-fix 解除）。
