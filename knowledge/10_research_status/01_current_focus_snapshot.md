---
id: ism-kb-10-research-status-current-focus-snapshot
name: "CURRENT_FOCUS Snapshot"
description: "docs/CURRENT_FOCUS.md 的結構化鏡像；Thread D paper + V6 production tag 雙軌序列化（Tier 1-4 2026-05-17~W6+）。⚠ 2 週有效。"
status: active
last_verified: 2026-05-18
content_nature: runtime-fact
doc_type: reference
verified_scope: "mirror of docs/CURRENT_FOCUS.md §2026-05-17 Tier 1-4 serialized execution"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-active-hypotheses
  - ism-kb-10-research-status-blockers-and-risks
  - ism-kb-10-research-status-next-milestones
  - ism-kb-01-project-overview-current-phase
  - ism-kb-03-pipelines-f1-baseline-canonical
  - ism-kb-09-conclusions-hypothesis-queue-snapshot
tags: [status, snapshot, current-focus, tier-1-4, thread-d, v6-production]
canonical_paths: [10_research_status/01_current_focus_snapshot.md]
alias_paths: []
---

# CURRENT_FOCUS Snapshot

> ⚠️ **此為 2026-05-18 快照（基於 2026-05-17 plan tender-pondering-blossom 收斂結果），2 週有效至 2026-06-01**
> 最新權威：[docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md) §2026-05-17

- 一句結論：Tier 1-4 序列化執行中；T1.1 Thread D 主軸正名 ✅ + T1.3 init-research scaffolding ✅ 完成；T1.2 V6 production tag 🔴 Hard Gate 待 5-day workflow 執行；Tier 2-4 排隊
- 適用對象：決策前快速了解現況、安排當前 priority、判斷阻塞解除順序
- 可直接執行命令：
  ```bash
  sed -n '15,120p' /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md
  ```

---

## 當前 Phase（2026-05-18）

**主軸雙軌平行**（決策 #1，2026-05-13 切換）：
- **Track A**：Thread D paper（Bioinformatics / NAR GB tool paper）— 4-6 週序列化
- **Track B**：Self-Phasing V6 production tag finalize — Hard Gate 5-day workflow

> **主軸正名**（決策 #5，2026-05-17 T1.1 完成）：
> 從「Self-Phasing causal chain validation」改為 **「TP-enriched phasing signatures (LOH × cross_het)」**

---

## Tier 1-4 進度追蹤（plan tender-pondering-blossom）

### Tier 1 (W3 2026-05-15~22) — 必須前置

| ID | 項目 | 狀態 | 詳細 |
|----|------|------|------|
| T1.1 | Thread D 主軸正名 | ✅ DONE 2026-05-16 | `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_*.md` 338→381 行（加 banner + §2.5 paradigm reframe） |
| T1.2 | V6 production tag finalize | 🔴 **Hard Gate** 待執行 | 5-day workflow（COLO829 V6 ISM → 7-sample marker → manifest → `git tag` → PI errata email） |
| T1.3 | init-research scaffolding | ✅ DONE 2026-05-16 | `research/thread_d_paper/`（165+108 行）+ `research/selfphasing_v6_production/`（154+122 行） |

### Tier 2 (W3-W4) — 證據強化

| ID | 項目 | 狀態 | 升級條件 |
|----|------|------|---------|
| T2.1 | Z-AUTO KDE 跨 4 樣本擴展 | ⏳ | ⭐3 → ⭐4 必要條件 |
| T2.2 | HCC1395 primary discovery 章節骨架 | ⏳ | Strategy A §3 |
| T2.3 | 6-sample replication cohort 章節骨架 | ⏳ | Strategy A §4（含 HCC1954/HCC1937） |

### Tier 3 (W4-W6) — Paper draft + 工程平行

| ID | 項目 | 狀態 | 備註 |
|----|------|------|------|
| T3.1 | Paper full outline | ⏳ | Abstract + 6 章 + 6 主圖（Tool paper 格式） |
| T3.2 | GitHub repo 整理 + public 化 | ⏳ | Reproducible release |
| T3.3 | Docker image build + tutorial | ⏳ | 1-hour install + run |
| T3.4 | Benchmark suite 公開化 | ⏳ | HG002 subset 或 HCC1395 SEQC2 公開部分 |
| T3.5 | Discussion 章節 | ⏳ | 63% gap + cancer-only + Normal BAM future |

### Tier 4 (W6+) — Reactive 擴展

| ID | 項目 | 觸發條件 |
|----|------|---------|
| T4.1 | Phase 2A Normal BAM cross-sample | G4 characterization 45% → 70% |
| T4.2 | GC/mappability/repeat 新軸 pilot | Reviewer 質疑 framework gap |
| T4.3 | PI Report 4-29 errata + V6 sign-off | T1.2 完成後一併打包 |
| T4.4 | HG002 non-cancer pilot | Reviewer 強質疑 cancer-only |

---

## T1.2 V6 Production Tag Hard Gate Workflow

| Day | Action | Gate Level |
|-----|--------|-----------|
| 1-2 | COLO829 V6 ISM 補完（Archive TO rerun + KDE-corrected） | 🟢 normal |
| 3 | 7-sample marker coverage + caller F1 比較（V3F vs V5 vs V6 vs SEQC2） | 🟢 normal |
| 4 | Binary commit hash 寫 `manifest.yaml` | 🟡 review |
| 4 | `git tag v6-prod-{YYYYMMDD}` | 🔴 **Hard Gate**（不可逆） |
| 5 | PI errata 5 條 + V6 sign-off written email draft | 🟡 review |
| 5 | User review email → send | 🔴 **Hard Gate**（送出後不可逆） |

**T1.2 通過後解鎖**：thread_d_paper Tier 2 Archive TO 7-sample rerun（V6 binary, ~10 hr parallel）+ T4.3 PI errata package。

---

## 已收斂的 9 條決策（2026-05-17 4 輪 Socratic）

| # | 領域 | 決策 |
|---|------|------|
| 1 | 主軸 | 雙軌平行（Thread D paper + V6 production） |
| 2 | G1-G5 | G1/G2/G4 characterization-only / G3 暫緩 / G5 降權 |
| 3 | 論文類型 | Tool paper（Bioinformatics / NAR GB） |
| 4 | 核心 contribution | read-level framework（不是新 ML 方法） |
| 5 | 主軸命名 | 「TP-enriched phasing signatures (LOH × cross_het)」 |
| 6 | 樣本階層 | Strategy A（HCC1395 primary + 6 replication） |
| 7 | Framework coverage | 接受 37% framework + 63% gap discussion |
| 8 | 限制 | Cancer-only 接受 limitation |
| 9 | 釋出 | 完全 reproducible（binary + Docker + benchmark + GitHub） |

---

## 2026-05-15 V6 ⭐3 PARTIAL POSITIVE（key evidence base）

**multi-agent fan-out A/B/C/D/E + Coordinator**（~3.5 hr 全執行）：
- phaseC 12 個 ISM run 整合 → step1_master_three_way.tsv (35,332 rows × 64 cols)
- 3 軸 50-cell grid + power gate (46% powered) + LR + 7 道 confound guard
- 4 FP zone deep dive (Z-OCH / Z-CHR8 / Z-GL / Z-AUTO)
- 4 樣本擴展 V6 ISM (H1437/H2009/HCC1954/HCC1937; COLO829 deferred)

**🤯 Paradigm reframe**：2/3 預設「FP-rich」zone 實際是 **TP-pure signatures**：
- Z-OCH: FP rate 0.017 << global 0.137（Fisher p=3.8e-62 for TP-enrichment）
- Z-GL: FP rate 0.003（0.022× global）→ gain on somatic hap = somatic signature

**✅ H4 POSITIVE — chr8 hotspot CN+AF 主導**：
- LR deviance：caller_af 0.393 > **CN 0.211** > HP 0.063 > LOH 0.038
- (LOH+CN) − HP = +0.186（3.7× threshold 0.05）

**🔍 V5 over-promote 直接證據**：
- Inner LOH NG=2: V5=8,136（+60% over V3F 5,064），TP rate 沒升
- V6 修補回 V3F 水準（5,353）
- V5/V3F top cell ratio 達 5.95×（Inner|cross_het_inv|cov_normal）

**⚠️ Framework coverage gap**：~63% FP 不被此 framework 解釋 → 需新軸（GC/mappability/repeat/SV）

主報告：[InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md](../../docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md)

---

## 與下方歷史區塊的關係

| Region | 角色 | 是否 deprecated |
|--------|------|---------------|
| **2026-05-17（本快照基礎）** | Plan tender-pondering-blossom 細化 + 序列化執行追蹤（live progress） | 當前權威 |
| **2026-05-15** | multi-agent fan-out paradigm reframe（evidence base） | 仍有效，T1.1/T1.3 執行依據 |
| **2026-05-13** | 原 3 週 V6 序列化估計 | ⛔ deprecated（保留 historical reference） |
| **2026-04-22**（前快照） | Phase 2 A+D 進行中、HCC1395 pilot POSITIVE | ⛔ 歷史鏡像 |

---

## Phase 1A Locked（不再改動）

- paired-pure ΔF1 = **+0.0112**（locked，~Cohen's small 0.2）
- TO ΔF1 = **-0.0206**（歷史觀察值；非定論，雙軌平行架構下仍可作為 Track A 替代探索基礎 — 依 `/scientific-rigor §2` L2 證據級）
- **完整 provenance（CI、方法、樣本、限制）** → [../03_pipelines/05_f1_baseline_canonical.md](../03_pipelines/05_f1_baseline_canonical.md)（SoT）

---

## 使用情境

若使用者問「現在在做什麼？」：
1. 讀本文件 §「Tier 1-4 進度追蹤」（最快概覽）
2. 若深度需求 → 讀 [docs/CURRENT_FOCUS.md §2026-05-17](../../docs/CURRENT_FOCUS.md)（含 plan 細節）
3. 需假說細節 → [02_active_hypotheses.md](02_active_hypotheses.md)
4. 需阻塞 → [04_blockers_and_risks.md](04_blockers_and_risks.md)
5. 需里程碑 → [05_next_milestones.md](05_next_milestones.md)
6. 需證據細節 → [03_evidence_ledger_format.md](03_evidence_ledger_format.md) + `research/autoresearch/evidence_ledger.jsonl`

---

## Governance 變更（2026-05-17）

- **v3 D2 分流** (commit `696c7c1`): CLAUDE.md 161 行 + AGENTS.md 286 行重構 — Claude Code 特定 vs 跨 agent governance
- **新元 skill** `/scientific-rigor` (commit `42217cf`): 14 章元方法論層，整合啟發式學習工作流映射
- **8 skill cross-reference** (commits `c1dde00` + `dce837f`): known-pitfalls / methodology-audit / verification-loop / validation-protocol / fast-learning-coach / memory-consolidation / check-staleness / auc-confound-guard 全雙向引用
- **Q5 erratum** (commit `f3611a7`): biorxiv/ensembl MCP 實測非僵屍（researcher claim 必實測升 L1）
- **5 KB snapshot refresh** (commit `d11b270`)：本快照在 2026-05-18 進入深度更新 phase

---

## 相關

- 最新權威：[../../docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)
- 活躍假說：[02_active_hypotheses.md](02_active_hypotheses.md)
- 阻塞：[04_blockers_and_risks.md](04_blockers_and_risks.md)
- Next milestones：[05_next_milestones.md](05_next_milestones.md)
- Hypothesis queue：[../09_conclusions/05_hypothesis_queue_snapshot.md](../09_conclusions/05_hypothesis_queue_snapshot.md)
- Plan：`~/.claude/plans/tender-pondering-blossom.md`
- Thread D paper plan：`InterSubMod/research/thread_d_paper/00_PLAN.md`
- V6 production plan：`InterSubMod/research/selfphasing_v6_production/00_PLAN.md`
