# Example: PI 報告新發現 — SCQA + STAR per item

> 場景: PI 1-on-1 (5min), Who=PI, Why=報告進度, What=新發現, How=口頭 + 簡單 .md companion

```markdown
---
framework: SCQA + STAR per item
audience: PI (1-on-1)
when: 5min lab meeting
applied_by: narrative-frame v0.1
cycle_id: 20260522-XXXX
---

# V6 修正 chr19 priority bug 後 HCC1395 F1 ↑0.022（cross-sample LOSO confirmed）

## §S — Situation

過去 V3F 在 chr19 報 752 reads 異常標 hp=33（非 1/2），HCC1395 paired pileup F1 0.7153。

(source: `InterSubMod/research/paired_priority_bug_audit/09_V6_caller_F1_verification.md:23`)

## §C — Complication

V3F 保守 hp=33 處理對齊 paired，但 V5 在 germline-absent 區（5,789 chr19 events）與 baseline 4.19:1 偏 HP1 完全相同 — Layer 1.5 設計缺陷繼承 priority bug。

(source: `InterSubMod/research/paired_priority_bug_audit/08_self_phasing_step_D_analysis.md:67`)

## §Q — Question

V6 能否 reproduce V5 的 paired alignment 同時不繼承 priority bug？

## §A — Answer

**Verdict POSITIVE on HCC1395**：V6 修 100% chr19 victims (752/752 reads)；HCC1395 paired-pileup ΔF1 = +0.022（vs V3F baseline，9.24× v1.0）。

### Supporting evidence（STAR per finding）

**Finding 1: chr19 victim rescue**
- S: 過去 chr19 752 reads hp=33
- T: V6 修正 hp tag priority
- A: 套用 V6 binary 重跑 paired pileup
- R: 100% rescued (752/752, source: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/INDEX.md:42`)

**Finding 2: F1 cross-sample uplift**
- S: HCC1395 LOSO baseline F1=0.7153
- T: 用 V6 重跑 + cross-sample LR
- A: 用 9 feature LR（drop caller_af, H_NEW_4）
- R: HCC1395 ΔF1=+0.00699 (single-sample marginal positive, source: `InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/loso_hnew_findings.md:51`)

### ⚠ Caveats

- HCC1395 +0.00699 是 single-sample 結果；其他 4 樣本仍 caller-F1-ceiling 卡住（best τ=0.10 keep all）
- LR direction 已 exhausted；下個 cycle 應 pivot 至 zone gate / boolean filter

---

## 預期 PI 追問（PREP per question）

**Q1: 為什麼 4 樣本沒效果？**
- P: caller_af direction-inconsistent 是 LOSO 災難主因
- R: HCC1395 d=+1.60 vs HCC1937 d=-1.41（cross-sample 矛盾）
- E: drop caller_af 後 LR 只能 train on weak coherent signal
- P: 應改 non-LR framework（zone rules / boolean filter）

**Q2: HCC1395 +0.00699 是 single-sample artifact 嗎？**
- P: 目前 inconclusive
- R: 4 sample 在 H_NEW_4 設定下被 caller-F1-ceiling 限制
- E: 找 low-F1 panel 才能驗 generalize
- P: 建議 Cycle 5 Path C — low-F1 panel 驗證

---

Framework: SCQA (Minto《Pyramid Principle》2020) + STAR per item (Behavioral Interviewing 1970s) + PREP for Q&A (Toastmasters)

Provenance footer:
- Cycle: 20260522-XXXX
- Sources: 3 .md files, 7 line citations
- Tier: 3 (full N1-N6)
- 5 秒測試: PASS（首段傳達 verdict + 數字）
```

---

## 為什麼選 SCQA + STAR + PREP（N2 推薦理由）

**5W 識別**:
- Who: PI（formality 高）
- Why: 報告進度 + 説服繼續方向
- What: 新發現 + 量化結果 + 預期問答
- When: 5min（lab quick）
- How: 口頭 + 簡單 .md companion

**主框架 SCQA** — 命中 PI / 報告進度 / 結論先行需求
**Sub-framework STAR per finding** — 每個 finding 用 behavioral case 強化 credibility
**Q&A 預備用 PREP** — 預測追問 + 30s 答結構

**不選擇**:
- Pixar Spine: PI 1-on-1 不需戲劇性故事弧
- ELI5: 受眾是 PI，過簡會 condescending
- Pure Pyramid: 缺敘事張力，PI 易失去 attention
