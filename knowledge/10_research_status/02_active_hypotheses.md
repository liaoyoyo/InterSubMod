---
id: ism-kb-10-research-status-active-hypotheses
name: "Active Hypotheses"
description: "當前 active 假說：Thread D paper TP-enriched signatures (H4 chr8 ⭐3 POSITIVE)、V6 production tag pre-registration、Z-AUTO KDE 4-sample 擴展、Pre-registration 機制（依 /scientific-rigor §7.1）。⚠ 2 週有效。"
status: active
last_verified: 2026-05-18
content_nature: runtime-fact
doc_type: reference
verified_scope: "active hypotheses against docs/CURRENT_FOCUS.md §2026-05-17 + §2026-05-15 Tier 1-4 plan tender-pondering-blossom"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-current-focus-snapshot
  - ism-kb-10-research-status-evidence-ledger-format
  - ism-kb-09-conclusions-hypothesis-queue-snapshot
tags: [status, hypotheses, active, thread-d, v6-production, pre-registration]
canonical_paths: [10_research_status/02_active_hypotheses.md]
alias_paths: []
---

<!-- STALE-REDIRECT-BANNER (scripts/stale_redirect_banner.sh) -->
> ⚠ **此檔為 2026-05-18 前後快照，可能已過時** — 現役主軸/狀態以 `InterSubMod/docs/CURRENT_FOCUS.md` 為準（主軸已於 2026-06-11 pivot 至 Subclonal reconstruction（取代 G6））。本檔僅供歷史對照，勿據此判斷現況。


# Active Hypotheses

> ⚠️ **此為 2026-05-18 快照**，最新權威：`research/autoresearch/hypothesis_queue.json` + `docs/CURRENT_FOCUS.md §2026-05-17`
>
> **本快照已從 2026-04-22 H011/H012 priority filter 鏡像深度更新到 Tier 1-4 Thread D paper × V6 production 雙軌假說階層。**

- 一句結論：當前 active 假說集中於 Thread D paper TP-enriched phasing signatures（H4 chr8 ⭐3 PARTIAL POSITIVE）+ V6 production tag pre-registration；早期 H011/H012 仍 adopted 但 priority 已降至 baseline supporting layer
- 適用對象：研究迴圈下一輪選擇、Tier 2 證據強化排序、Pre-registration 條件確認
- 可直接執行命令：
  ```bash
  # Tier 1-4 進度
  sed -n '15,80p' /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md
  # 假說 queue
  jq '.hypotheses[] | select(.status == "adopted")' \
    /big7_disk/liaoyoyo2001/InterSubMod/research/autoresearch/hypothesis_queue.json
  ```

---

## §1 當前主軸假說（Thread D paper / V6 production）

### H_THREAD_D_MAIN：TP-enriched phasing signatures (LOH × cross_het)
- **主軸命名**（決策 #5，2026-05-17 T1.1）
- **Priority**：99（最高，雙軌核心）
- **Status**：⭐3 PARTIAL POSITIVE → 升 ⭐4 條件 T2.1 + T2.3
- **Evidence Tier**：L2（依 `/scientific-rigor §2`）— 主要數據支持 + 機制部分清楚
- **Track**：Thread D paper
- **Pre-registration 條件**（依 `/scientific-rigor §7.1`）：
  - H 預測：cross_het bucket 在 ≥4 樣本同方向 TP-enrichment
  - 否證條件：Wilcoxon p > 0.1 across samples OR Δ < +0.005
  - Decision threshold：⭐3 → ⭐4 需 (a) Z-AUTO KDE 跨 4 樣本 recur (b) 加新軸測 63% gap
- **參考**：[../../docs/CURRENT_FOCUS.md §2026-05-15](../../docs/CURRENT_FOCUS.md)

### H4：HCC1395 chr8 hotspot CN+AF dominant
- **Status**：✅ POSITIVE（2026-05-15 multi-agent fan-out）
- **Evidence Tier**：L2（caller_af LR deviance 0.393 > CN 0.211 > HP 0.063 > LOH 0.038）
- **Effect Size**（依 `/scientific-rigor §3`）：(LOH+CN) − HP = +0.186（3.7× threshold 0.05；中等）
- **Sample-specific 限制**：chr8 FP enrichment 2.31× highest of 23 chr，但無法跨樣本 generalize
- **後續**：併入 H_THREAD_D_MAIN §3 HCC1395 primary 章節

### H_V6_PROD：V6 binary production tag finalize
- **Track**：V6 production
- **Priority**：95（T1.2 Hard Gate）
- **Pre-registration 條件**：
  - H 預測：V6 marker coverage +9% over V3F + caller F1 不變
  - 否證條件：V6 F1 < V5 F1 OR marker coverage < V3F baseline
  - Decision threshold：通過 → `git tag v6-prod-{YYYYMMDD}` + PI errata email
- **Hard Gate Workflow**：5-day（Day 4 tag + Day 5 PI email 不可逆）
- **參考**：`InterSubMod/research/selfphasing_v6_production/00_PLAN.md`

### H_Z_AUTO_RECUR：Z-AUTO KDE mechanism cross-sample replication
- **Track**：Thread D paper Tier 2
- **Status**：⏳ 未驗證（T2.1）
- **Pre-registration 條件**：
  - H 預測：Z-AUTO mechanism 在 H1437/H2009/HCC1954/HCC1937 重複出現
  - 否證條件：4 樣本中 ≥2 樣本 KDE pattern 不同方向
  - Decision threshold：通過 → ⭐4 升級

---

## §2 Paradigm Reframe（2026-05-15 critical update）

| Zone | 舊框架預設 | 新框架實證 | Mechanism |
|------|-----------|----------|----------|
| Z-OCH (Outer cross_het) | FP-rich | **TP-pure signature** (FP rate 0.017 vs global 0.137) | cross_het = somatic-evidence marker |
| Z-GL (Inner gain+LOH) | FP-rich | **TP-pure signature** (FP rate 0.003) | gain on somatic hap = somatic signature |
| Z-CHR8 | (新) | FP-rich (2.31× enrichment) | chr8 sample-specific hotspot |
| Z-AUTO | (新) | KDE-mechanism (待 cross-sample) | Coverage-aware correction |

**意義**：先前以為「FP markers」的 zones 實為 TP signatures；framework 方向重定向。

---

## §3 V5 over-promote 直接證據

依 `/scientific-rigor §4` DAG（confounder/collider/mediator 分流）：
- Inner LOH NG=2 region：V5=8,136 vs V3F=5,064（+60% over-promote）
- V5 TP rate **沒升**（mechanism 失敗證據）
- V6 修補回 V3F 水準（5,353）
- V5/V3F top cell ratio 達 5.95×（Inner|cross_het_inv|cov_normal），**集中 cross_het bucket**

→ Layer 1.5 機制只在 somatic-fallback heterozygous reads 作用（特定 trade-off，非全域 regression）

→ 對應 memory `project_v5_v6_tradeoff_sp123` + `project_v5_layer15_design_caveat`

---

## §4 Adopted 假說（baseline supporting layer）

> 自 2026-04-22 以來保留但 priority 已降至 baseline；雙軌主軸下不再 active 推進。

| ID | 名稱 | Delta F1 | 角色 |
|----|------|---------|------|
| H011 | QS≥50 TO rescue | +0.008556 | adopted（TO track，主軸已 pivot） |
| H012 | GQ≥3 | +0.009365 | adopted（baseline filter） |
| H_COMBO | 組合 filter | — | adopted（infra） |
| H_KDE_001 | KDE audit logging | — | adopted（infra） |

---

## §5 Conditional Positive / Annotation Only

| ID | 狀態 | 備註 |
|----|------|------|
| LOH-constrained phasing discovery (NG=2) | conditional_positive | 2026-04-22 新發現，Inner 93-99% same-hap（6/6 樣本）；TP gap +0.37（更正 NG 非 methylation 訊號）|
| LOH × AF × Methylation Paired POSITIVE | conditional_positive | 7/7 樣本，TO 層 pivot 為 phasing；Paired 層保留需獨立 phasing-vs-methylation 驗證 |

---

## §6 Rejected 假說（防重複調查）

> 跨參考至 [09_conclusions/05_hypothesis_queue_snapshot.md](../09_conclusions/05_hypothesis_queue_snapshot.md) §「狀態分類」+ memory「Concluded」區。

| ID | 名稱 | 理由 |
|----|------|------|
| H001 | HPP bias | rejected |
| H003/H007 | (TO haplotagging) | TO haplotagging 不可靠 |
| H005 | VAF/AlleleDelta in TO | rejected |
| H008/H009 | (CramersV track-dep) | rejected |
| Read-Level Germline FP Phase 1 | CONDITIONAL NO-GO | LOSO AUC 0.721 但 FP removal=0% |
| Q5 biorxiv/ensembl MCP「僵屍」誤判 | NEGATIVE→修正 | commit `f3611a7` erratum；實測非僵屍 |

---

## §7 最近 Cycle 結果（2026-04-11 → 2026-05-17）

從 `research/autoresearch/evidence_ledger.jsonl` + multi-agent runs：

| Cycle / Run | 日期 | 主題 | Verdict |
|-------------|------|------|---------|
| cycle_011-019 | 2026-04-11~21 | PON / Phase B/C/D / KDE / germline-hp-only | ✅/⚠/❌ |
| multi-agent A-E + coord | 2026-05-15 | V6 TPFP HP LOH CN characterization | ⭐3 PARTIAL POSITIVE |
| T1.1 Thread D 主軸正名 | 2026-05-16 | banner + §2.5 paradigm reframe | ✅ DONE |
| T1.3 init-research scaffolding | 2026-05-16 | thread_d_paper + selfphasing_v6_production | ✅ DONE |
| governance v3 + scientific-rigor | 2026-05-17 | D2 分流 + 元 skill + 8 cross-ref | ✅ commit f3611a7 / df76e45 |
| T1.2 V6 production tag | ⏳ pending | 5-day Hard Gate workflow | 🔴 待執行 |

---

## §8 下一候選（Tier 2-3 排隊）

依 Tier 1-4 序列化執行：
1. **T1.2 完成優先**（5-day workflow + PI errata email）
2. T2.1 Z-AUTO KDE 4 樣本擴展（驗證 mechanism 是否 recur）
3. T2.2 HCC1395 primary discovery 章節骨架
4. T2.3 6-sample replication cohort（HCC1954 / HCC1937 / H1437 / H2009 / COLO829 / DORADO）
5. T3.1 Paper full outline + 6 主圖
6. T3.2-T3.4 GitHub / Docker / Benchmark public release

---

## §9 Pre-registration 機制（依 /scientific-rigor §7.1，2026-05-17 新增）

所有 Tier 2-4 新研究方向開跑前須在 `research/<topic>/00_INDEX.md` 強制 3 欄：
- **H 預測**（confirmatory 假設）
- **否證條件**（NO-GO threshold）
- **Decision threshold**（升 ⭐4 / 廢棄 / 推進 paper）

達 NO-GO 條件 → verdict 不可事後改寫（Hard Gate）

範本：`InterSubMod/templates/research_index.md`（commit `a031d21`）

---

## §10 相關

- 當前主軸：[01_current_focus_snapshot.md](01_current_focus_snapshot.md)
- 阻塞：[04_blockers_and_risks.md](04_blockers_and_risks.md)
- Milestones：[05_next_milestones.md](05_next_milestones.md)
- Queue snapshot：[../09_conclusions/05_hypothesis_queue_snapshot.md](../09_conclusions/05_hypothesis_queue_snapshot.md)
- Evidence ledger：[03_evidence_ledger_format.md](03_evidence_ledger_format.md)
- Pre-registration 範本：`InterSubMod/templates/research_index.md`
- /scientific-rigor §7.1：`InterSubMod/.claude/skills/scientific-rigor/SKILL.md`
