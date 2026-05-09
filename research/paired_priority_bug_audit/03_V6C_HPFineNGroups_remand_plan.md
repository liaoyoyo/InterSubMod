<!--
build_date: 2026-05-10
agent: V6-C HPFineNGroups subclone marker 重評估計劃
status: validated (plan)
report_class: research-plan
audience: PI / lab member / 自己未來
parent_v6_eval: InterSubMod/research/paired_priority_bug_audit/02_V6_proposal_evaluation.md
parent_synthesis: InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md
inputs:
  - HPFineNGroups marker memory (⭐3 pipeline-dependent, 2026-04-23 降級)
  - LOH-constrained phasing discovery memory (NG=2 Inner same-hap 6/7 一致, 2026-04-22)
  - F pilot NG=4+AF<0.4+NR≥80 NonLOH TP rate 92.8% (2026-04-18)
  - Phase 1 audit (2026-04-21) flag=on NGroups ≥3=0
  - obs18 6 TO 樣本 same-hap composition data
  - paired_priority_bug_audit Step D V5 Layer 1.5 設計缺陷
outputs:
  - 本檔（V6-C 三 phase 計劃）
verdict: 三 phase 階梯式驗證 — Phase A (理論推導, 1 hr) → Phase B (chr19 子集 flag=on/off 對比, 半天) → Phase C (7 樣本 master × flag 重驗, 1-2 day)；先做 Phase A 看是否值得 B/C 投入
last_verified: 2026-05-10
-->

# V6-C HPFineNGroups Subclone Marker 重評估計劃

## 0. TL;DR

> **V6 提案評估發現** ISM `germline_hp_only` flag 已等效；2026-04-21 Phase 1 audit 顯示 flag=on 時 NGroups ≥3 完全消失（22,732 → 0 TP regions）。這對 **HPFineNGroups subclone marker** 有重大意義 — 原 ⭐4 結論可能是 priority bug 的 phasing artifact，不是真實生物學 sub-clone signal。本計劃用三 phase 階梯式驗證確認 marker 機制是否真實，影響 LOH-constrained phasing 主軸論文 thesis。

---

## 1. 問題定義

### 1.1 待驗證假說

**H-NG-1**：HPFineNGroups（特別是 NG=2 Inner same-hap pattern）是真實 sub-clone heterogeneity signal，非 priority bug 副產物
**H-NG-2**（對立）：HPFineNGroups NG≥3 部分 / NG=2 部分訊號是 priority bug 在 BAM 內 hp=11/21 額外分類造成的 artifact

### 1.2 已知事實

| 觀察 | 數值 / 結論 | 來源 |
|---|---|---|
| F pilot NG=4+AF<0.4+NR≥80 NonLOH | TP rate 0.9281（n=14,197）跨 5/7 樣本 ≥0.85 | F pilot 4-18 主報告 |
| obs18 NG=2 Inner × same_HP composition | 93-99% same-hap（HP1+HP1-1 或 HP2+HP2-1）跨 6/6 樣本一致 | LOH-constrained discovery 4-22 |
| obs18 TP gap (Inner same_HP1 − Outer cross_het) | +0.05 ~ +0.52, median +0.37, 6/6 正向 | 同上 |
| Phase 1 audit flag=on NGroups ≥3 TP regions | 22,732 → **0**（完全消失）| 2026-04-21 Phase 1 |
| Phase 1 audit HPFineNGroups AUC delta | flag=on 衰退 -0.026（filter 角度 NEGATIVE）| 同上 |
| HCC1395 TO ClairS-TO raw split TP rate | 0.694 ≈ baseline 0.699（Fisher odds 0.913 反向 p=3.5e-3）| 4-23 週報降級 |
| Step D paired germline-absent | 5,789 events, V5/baseline 4.19:1 偏 HP1 | 2026-05-09 |

### 1.3 機制兩個可能解讀

| 解讀 A：marker 真實 | 解讀 B：marker artifact |
|---|---|
| NG=2 Inner same-hap = LOH 區域真實 sub-clone phasing 訊號 | NG=2 Inner same-hap = LOH 區域 + priority bug 把 somatic 都聚到同 haplotype 的 artifact |
| obs18 6/6 一致 → 強生物學 signal | obs18 6/6 一致 → 6/6 都被 priority bug 影響 |
| F pilot 5/7 樣本 ≥0.85 → 跨樣本驗證 | F pilot 5/7 樣本 ≥0.85 → 5/7 樣本 priority bug 表現一致 |
| flag=on NG≥3 消失 = flag 過度 demote 真實訊號 | flag=on NG≥3 消失 = flag 正確移除 priority bug 副產物 |

→ 兩個解讀在現有資料上都 **partial fit**；需要 **flag=on master dataset × 7 樣本** 才能定論。

---

## 2. 三 phase 階梯式計劃

### Phase A — 理論推導（1-2 hr，不需重跑）

**目標**：用既有 5/8 §8.6 paired audit + obs18 + F pilot 數據做機制推導，估計 marker artifact rate 上界 / 下界。

**步驟**：

A1. 對 obs18 NG=2 Inner same-hap 數據做 priority bug 受影響度估計：
   - 5,789 chr19 germline-absent events × 全基因組擴展約 ~150K events
   - 對比 obs18 NG=2 Inner total events (per sample)
   - 計算 「priority bug 影響的 NG=2 events 占 obs18 same-hap 觀察的比例」

A2. 從 F pilot 6/7 樣本 NG=4+AF<0.4+NR≥80 NonLOH 14,197 regions 推算：
   - 多少 region 在 germline-absent area（priority bug 主要影響範圍）
   - 多少 region 純 germline_vote>0（priority bug 不影響的「乾淨」signal）

A3. 結論輸出：
   - 若 priority bug 影響 < 10% F pilot regions → marker 大部分真實
   - 若影響 > 50% → marker 可能 artifact 主導
   - 若 10-50% → 需 Phase B/C 量化才能定論

**Deliverable**：`InterSubMod/research/paired_priority_bug_audit/04_V6C_phaseA_theory_estimate.md`

**Decision gate**：依 A3 比例決定是否進 Phase B。

### Phase B — chr19 子集 ISM run flag=on / flag=off 對比（~半天，需重跑 ISM）

**目標**：對 HCC1395 TO chr19 子集（速度優先）跑 ISM with flag=on，與既有 flag=off 對比 NG distribution + same-hap composition 變化。

**前提**：Phase A 結論顯示 priority bug 影響度 ≥10%（足夠值得驗證）。

**步驟**：

B1. ISM chr19 重跑（flag=off baseline 已有；flag=on 需新跑）：
   ```bash
   ./scripts/run_vcf_all_snv.sh --mode chr19-verification --germline-hp-only=true
   ```
   估時：~5 min（chr19 子集）

B2. 比對 NG distribution：
   - flag=off NG=1/2/3/4 reigon counts
   - flag=on NG=1/2/3/4 region counts
   - 計算 same-hap fraction in NG=2 Inner under flag=on

B3. 若 flag=on 下：
   - NG=2 Inner same-hap fraction 仍 ≥80% → marker 真實（保留 ⭐3 / 升級 ⭐4）
   - NG=2 Inner same-hap fraction < 50% → marker artifact 主導（降級 ⭐2）
   - 50-80% → 中間區，需 Phase C 7 樣本確認

**Deliverable**：`InterSubMod/research/paired_priority_bug_audit/05_V6C_phaseB_chr19_remand.md`

**Decision gate**：依 B3 結果決定是否進 Phase C。

### Phase C — 7 樣本 master × flag 全量重驗（1-2 day，需重跑 ISM × 7 樣本 × 2 flag）

**目標**：對 7 樣本 master_extended × flag=on/off 完整重跑，重算 F pilot canonical filter + obs18 6 樣本 NG=2 Inner same-hap 結論。

**前提**：Phase B 結論顯示中間區或顯示部分樣本 flag=on 有保留 marker（需跨樣本驗證）。

**步驟**：

C1. 7 樣本 ISM 全量 flag=on 跑：
   - HCC1395 / HCC1395_DORADO / HCC1937 / HCC1954 / H1437 / H2009 / COLO829
   - 估時 ~10 hr parallel（依 Archive 文件 `20260422_Archive_TO_Rerun_ISM_Requirement_01.md`）
   - 條件：每 sample 對齊 KDE-fixed binary + 兩 flag

C2. 重算指標：
   - F pilot canonical filter NG=4+AF<0.4+NR≥80 NonLOH **flag=on** 下 TP rate per-sample
   - obs18 NG=2 Inner same-hap composition **flag=on** 下 per-sample
   - Per-sample Cohen's h flag=on 變化

C3. 結論決策：
   | flag=on 下保留信號 | 解讀 | marker 結論 |
   |---|---|---|
   | TP rate ≥0.85 in 5/7 + same-hap ≥80% in 6/7 | 真實生物學 marker | ⭐4 升級 |
   | TP rate <0.70 in ≥4/7 + same-hap <50% | 主要 priority bug artifact | ⭐2 降級 |
   | 中間區 | 部分真實 + 部分 artifact | ⭐3 pipeline-dependent 維持 |

**Deliverable**：`InterSubMod/research/paired_priority_bug_audit/06_V6C_phaseC_7sample_remand.md`

**對下游影響**：
- **LOH-constrained phasing 主軸論文 thesis 修訂**（NG=2 same-hap 結論可能需精確化）
- **MEMORY `project_hpfinengroups_subclone_marker.md` 升級降級**
- **MEMORY `project_loh_constrained_phasing_discovery.md` 機制重新詮釋**

---

## 3. 風險與已知 caveat

| 風險 ID | 描述 | Mitigation |
|---|---|---|
| R1 | Phase 1 audit 中間檔在 /tmp/ 已被災情清掉 | 不能 reuse；Phase B/C 需從頭跑 |
| R2 | 7 樣本 master_extended 是 flag=off 跑出 → 不能直接對比 flag=on | 必須跑新一輪 flag=on |
| R3 | Phase B chr19 子集可能 underpower（germline-absent 比例小）| 跑全 7 sample chr19 + 看 cross-sample 一致性 |
| R4 | flag=on 下 NGroups distribution 重算可能與 obs18 模式不同（demote 後 region 統計改變）| Phase B 多 metrics 對比；不單看 NG count |
| R5 | LOH-constrained discovery 機制可能需修訂 | 接受並記錄；論文 thesis 需修訂 |

---

## 4. 預期時程

```
Phase A   1-2 hr            ← 理論推導，不需重跑
   │
   ├── 結論 < 10% 影響 → marker 真實，停止（不需 B/C）
   ├── 結論 10-50% 影響 → 進 Phase B
   └── 結論 > 50% 影響 → 直接進 Phase C
   
Phase B   半天 (~5 hr)       ← chr19 子集 ISM 重跑 flag=on
   │
   ├── same-hap ≥80% → marker 真實，停止
   ├── same-hap <50% → 直接進 Phase C 確認
   └── 50-80% → 進 Phase C 7 樣本 cross-validation

Phase C   1-2 day             ← 7 樣本 master × flag=on 全量重驗
   └── 結論決策（升級 ⭐4 / 維持 ⭐3 / 降級 ⭐2）
```

**最快路徑**：Phase A 顯示影響微小 → 1-2 hr 結束。
**完整路徑**：A → B → C → 1-2 day total。

---

## 5. 啟動條件

User 確認後啟動 Phase A（不需新 BAM、純分析）。

Phase B 啟動條件：
- Phase A 結論影響度 ≥10%
- 用戶確認 chr19 ISM 重跑（5 min）
- KDE-fixed binary 仍可用

Phase C 啟動條件：
- Phase B 結論「中間區」或「部分樣本保留」
- 用戶確認 7 樣本全量重跑（~10 hr parallel）
- 7 樣本 BAM + germline VCF 都齊全

---

## 6. 對主軸論文 thesis 的潛在影響

### 6.1 LOH-constrained phasing discovery (4-22 主軸)

原 thesis：「**TO mode 下 NG=2 在 LOH Inner 區 93-99% same-haplotype 分裂，TP rate ~0.93；Outer 區以 cross-het 為主**」

可能修訂方向：
- **若 V6-C 確認 marker 真實** → thesis 不變，但加 caveat「signal 不依賴 priority bug, 在 flag=on 下仍 ≥80%」
- **若 V6-C 揭露 artifact 主導** → thesis 大幅修訂為「在 priority bug 主導機制下 NG=2 Inner 呈 same-hap，去除 artifact 後 signal 弱化或消失」
- **混合（最可能）** → thesis 加精確化「signal 包含真實 LOH-constrained phasing + priority bug 副產物，比例量化於 V6-C」

### 6.2 HPFineNGroups marker 結論層級

| V6-C 結論 | marker 等級 | F pilot canonical filter |
|---|---|---|
| 真實 | ⭐4（升級回原階）| 保留 NG=4+AF<0.4+NR≥80 NonLOH 92.8% |
| 中間 | ⭐3 pipeline-dependent（維持）| 加 「flag=off only」caveat |
| Artifact 主導 | ⭐2 重大降級 | filter 撤回；改為 sanity check 而非 marker |

---

## 7. 引用

- F pilot 主報告：`InterSubMod/docs/experiments/in_progress/2026/04/20260418_F_HPFineNGroups_deepening_POSITIVE_01.md`
- F pilot 子專案：`InterSubMod/research/F_hpfinengroups_deepening/`
- LOH-constrained discovery：`InterSubMod/research/tpfp_loh_af_kde_discrimination/09_TO_sample_af_lohside_ng.md`
- obs18 data：`InterSubMod/research/tpfp_loh_af_kde_discrimination/data/obs18_*`
- Phase 1 audit：`InterSubMod/docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_HCC1395_01.md`
- 5/8 整合報告 §8.6 + V6 evaluation：`InterSubMod/research/paired_priority_bug_audit/02_V6_proposal_evaluation.md`
