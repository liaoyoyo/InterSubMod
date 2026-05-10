<!--
build_date: 2026-05-10
agent: pptx-build P2 outline
status: in_progress (awaiting C2 user ack)
report_class: presentation_outline
audience: PI / lab meeting (領域專家但不熟 C++)
target_duration_min: 20
report_type: improvement_report (主) + academic_defense (深度) + concept_walkthrough (機制段) + data_showcase (量化段)
main_thesis: "Self-Phasing 修補成熟，V5 Layer 1.5 germline-absent 區仍待補強"
parent_synthesis: InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md
parent_errata: InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md
narrative_order: observation-first (對齊母稿)
-->

# Self-Phasing 整合 PPT — Outline (P2)

## Main Thesis (≤30 字)

**Self-Phasing 修補成熟，V5 Layer 1.5 germline-absent 區仍待補強**

雙焦點平衡敘事：
- **正面焦點**：V3F + V5 兩層修補在 read-level 100% 修正全基因組 34,855 victims；20 指標 0 regression；caller F1 三版完全相同（V5 不改 caller）
- **caveat 焦點**：5/9 paired cross-ref 揭露 V5 Layer 1.5 在 germline-absent 區與 baseline 4.19:1 偏移完全相同（priority bug 的 feature 化非修補）；V3F 標 hp=33 反而更穩健

## 6 Section + 18 Main Slide + 3 Backup = 21 slide

> **C2 後更新（2026-05-10）**：S7 拆 2（17 errata / 18 follow-up）；conclusion ratio 從 11.8% → 17.6% ✓

### S0. Cover + TL;DR (2 slide)

**Section thesis**：修補主線 100% 但 5/9 新發現有 germline-absent 區 V5 設計缺陷

| Slide | Title (assertion-evidence) | Focal point (≤20 字) |
|-------|---------------------------|---------------------|
| 01 | Self-Phasing 整合觀察與 V5 Layer 1.5 設計缺陷揭露 | (cover) |
| 02 | 修補成熟 V3F+V5 100% read-level；V5 Layer 1.5 germline-absent 區仍待補強 | TL;DR：5 個關鍵數字 + 1 caveat |

### S1. 觀察起點 — 17.3:1 全基因組偏移 (2 slide)

**Section thesis**：LongPhase-TO 全基因組 94.6% somatic→HP1 是 systematic bias 非樣本性質

| Slide | Title | Focal point | 母稿來源 |
|-------|-------|-------------|---------|
| 03 | 全基因組 HP1:HP2 = 17.3:1 偏移；3 條獨立論證真實 systematic bias | 生物學 + 跨 23 chr 一致 + paired 對照 1:1 | §2.1 |
| 04 | chr19 SP1/2/3 三個近 100% 失衡位點 — IGV 6-BAM 並列鐵證 | 113:0 / 109:1 / 108:0 + V5 翻回 paired | §2.2 + IGV |

### S2. 機制 — phasing + tagging 兩層 bug (3 slide)

**Section thesis**：phasing 層球員兼裁判 + tagging 層 getVote priority bug，1 票 somatic 蓋過 5 票 germline

| Slide | Title | Focal point | 母稿來源 |
|-------|-------|-------------|---------|
| 05 | phasing 層 — somatic 進 graph 球員兼裁判，自我增強迴圈 | germline 50/50 vs somatic 100/0 共現 | §3.1-3.2 |
| 06 | tagging 層 — `getVote()` vector 順序錯 + break early；1 票 somatic 觸發 priority bug | F1 mechanism + 5-vote 範例 read | §3.3 + F1 |
| 07 | 兩層 bug 兩層修補對應表 — 任一層單獨修不夠 | phasing → PON-only flag / tagging → V3F + V5 | §3.4 |

### S3. 量化鐵證 — read-level 34,855 victims (2 slide)

**Section thesis**：chr19 752 + 全基因組 34,855 victims，V3F+V5 修正率 100% 單向

| Slide | Title | Focal point | 母稿來源 |
|-------|-------|-------------|---------|
| 08 | chr19 read-level 752 victims — 4-path 驗證 3.5/4 PASS；100% 單向 baseline=11→V3F=21 | F4 chr19 1Mb hotspot + 統一指紋表 | §4.1, F4 |
| 09 | 全基因組擴展 34,855 victims（46×）；priority bug 主要分佈 chr7/chr2/chr1 — chr19 占比僅 2.16% | F2 per-chr 雙 panel | §4.2-4.3, F2 |

### S4. 修補設計演進 — 5 commits 兩層三版 (2 slide)

**Section thesis**：baseline → V3F (兩層 + INDEL guard) → V5 (Layer 1.5 + ploidy fix + threshold)，5 commit stacking 才達完整

| Slide | Title | Focal point | 母稿來源 |
|-------|-------|-------------|---------|
| 10 | 5 commits 時間軸 — phasing layer 藍 / tagging layer 綠 / 跨層紫 | F3 timeline + commit hash | §5.1-5.4, F3 |
| 11 | `getVote()` 三版程式碼對照 — V3F 兩層 + V5 Layer 1.5 fallback | code 三版 side-by-side | §5.6 |

### S5. 驗證 + no-regression — 20 指標全綠 (3 slide)

**Section thesis**：read-level 100% 修正 + paired GT +13.3 pp + 20 指標 0 regression + caller F1 三版完全相同

| Slide | Title | Focal point | 母稿來源 |
|-------|-------|-------------|---------|
| 12 | SP1/2/3 修正後對齊 paired 3/3；HP1:HP2 17.3:1 → ~1:1 | IGV V5 vs paired 翻轉示意 + 結構表 | §6.3-6.4 |
| 13 | 20 指標 no regression — ISM/HP_Ratio/methylation/paired GT/HP 結構/LOH | metric table + +13.3 pp / +99.7% 高亮 | §8.5.1 |
| 14 | Caller F1 vs SEQC2 三版完全相同（V5 不改 caller）；purity 0.6 完整對照 0 critical regression | 0.93 = 0.7166 / 0.6 = 0.6273 + 因果鏈 + **「→ 但 5/9 paired cross-ref 揭露另一面...」cliffhanger** | §8.5.2 |

### S6. 5/9 新發現 — V5 Layer 1.5 設計缺陷 (2 slide)

**Section thesis**：germline-absent 區 V5 與 baseline 4.19:1 完全相同 — priority bug 的 feature 化非修補；V3F 標 hp=33 反而更穩健

| Slide | Title | Focal point | 母稿來源 |
|-------|-------|-------------|---------|
| 15 | paired mode 整體無 priority bug — som_ratio 0.462 跨 windows 全範圍分布 | F6 paired vs TO chr19 | §8.6.2-3 |
| 16 | germline-absent 區 V5 = baseline 4.19:1；V3F 全標 hp=33 保守不選邊更穩健 | F7 三版 cross-tab + V5 Layer 1.5 機制詮釋 | §8.6.4-5 |

### S7. 結論 + 5 errata + 後續 (2 slide — C2 確認拆 2)

**Section thesis**：V5 仍可作 production baseline；5 條 PI errata 已 patch；3 個 paired follow-up + T3 跨樣本 + T1.3 ablation 待跑

| Slide | Title | Focal point | 母稿來源 |
|-------|-------|-------------|---------|
| 17 | 5 條 PI errata 已 patch（E1-E5） | E1 chr19 hotspot 降級 / E2 V5 commit 完成 / E3 證據鏈升級 / E4 V5=Pass1only / E5 Layer 1.5 設計缺陷 | §9.2 |
| 18 | 整體成熟度 + 5 項 follow-up | 12 維度成熟度（10 ✅ / 1 ⚠️ / 1 待跑）+ F-paired-D1/2/3 + T3 + T1.3 | §9.1, §9.3 |

### Q&A Backup (3 slide)

| Slide | Title | 觸發 PI 問題 |
|-------|-------|------------|
| B1 | Pass 2 second round 機制 — 只重跑 2-point edgeConnectResult | 「Pass 2 為什麼只跑 2-point？」 |
| B2 | purity 0.6 樣本完整對照表 — 6 caller F1 + 9 結構指標 | 「低純度樣本是否變差？」 |
| B3 | 7 樣本擴展 + cnLOH 雙親同源待開放方向 | 「跨樣本一致性？」 |

---

## 自我審查（架構師視角）

### 跳躍合理性

| 跳躍 | 風險 | 緩解 |
|------|------|------|
| S1 → S2 | 低 | 母稿 §2.3 已 funnel 4 問，自然進機制 |
| S3 → S4 | 低 | 量化 confirm 機制真實 → 才解釋修補 |
| **S5 → S6** | **中** | **S5 結束「修對了」氣氛，S6 突然「沒修對」需 cliffhanger transition** |
| S6 → S7 | 低 | 母稿已採「成熟 + caveat」雙重結論 |

**S5 → S6 transition 處理**：S5 第 14 slide 結尾用「**但 5/9 paired cross-ref 揭露另一面**」一句 hook，銜接 S6（手法：inductive ordering 對齊母稿 frontmatter）。

### 17 主 slide ratio check（playbook §20.D）

| 類別 | slide 數 | 占比 | playbook 規格 |
|------|---------|------|--------------|
| Definition + Prerequisite | S2 (3) | 17.6% | ≤ 25% ✓ |
| Body | S1+S3+S4+S5+S6 (12) | 70.6% | ≥ 60% ✓ |
| Conclusion | S0+S7 (2) | 11.8% | ≥ 15% ⚠ 略低 |

→ Conclusion 略低於 15%，建議 S7 拆 2 張（一張 errata + 一張 follow-up）→ 總 18 主 slide。

### 個人風格規則套用清單

| 規則 | 應用 slide |
|------|----------|
| R-G2 術語密度 | 高密度: S2 (05/06/07), S4 (10/11), S6 (15/16) |
| R-G3 Agent-N 數字驗證 | S0/S2 (TL;DR 數字), S3 (752/34,855), S5 (各指標), S6 (4.19:1) |
| Metric scope 明示 | S5 (chr19 vs 全基因組 / Pass 1 only vs Pass 1+2 / 0.93 vs 0.6 purity) |

---

## C2 confirmation

請審 outline 6 段 + thesis sentence + slide ratio：
1. 6 段切割是否合理？S5/S6 是否需要更多 slide？
2. S5 → S6 transition 的 cliffhanger 設計是否接受？
3. S7 是否拆 2 張（提升 conclusion ratio 到 ≥ 15%）？
4. Q&A backup 3 張選題（Pass 2 / purity 0.6 / 跨樣本）是否符合 PI 預期問題？
