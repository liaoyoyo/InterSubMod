<!--
build_date: 2026-05-10
agent: V3F vs V5 BAM head-to-head evaluation (chr19)
status: validated
report_class: comparative-empirical
audience: PI / lab member / 自己未來
parent_phase_B: InterSubMod/research/paired_priority_bug_audit/05_V6C_phaseB_findings.md
inputs:
  - InterSubMod/research/paired_priority_bug_audit/v6c_phaseB_runs/ (V3F BAM)
  - InterSubMod/research/paired_priority_bug_audit/v6c_phaseB_v5bam_runs/ (V5 BAM)
  - 5/8 主報告 §8.4 / §8.5 (V5 vs baseline 全指標)
  - 5/9 paired audit Step D (V5 Layer 1.5 4.19:1 偏移)
outputs:
  - 本檔（V3F vs V5 統合 evaluation）
  - InterSubMod/research/paired_priority_bug_audit/v6c_phaseB_v5bam_runs/v3f_vs_v5_summary.tsv
verdict: chr19 head-to-head + 5/8 + 5/9 證據鏈合計給出**清楚分區結論** — germline-existent 區域兩者並列、germline-absent 區域 V3F 勝、caller F1 平手、region-level marker coverage V3F 略勝 (+23.5% TP regions 抓到, marker rate +0.023~0.030 pp)；V5 仍可作 production baseline (caller F1 不變, sanity 全通過)，但若 ISM downstream marker 為主軸，V3F BAM 是略佳選擇
last_verified: 2026-05-10
decision: V5 不撤回 production baseline；germline-absent 區域改回 V3F 設計（V6-D binary patch）為高 ROI follow-up；7 樣本 V5 BAM × 完整 marker filter 為下一驗證步驟
-->

# V3F vs V5 BAM Head-to-Head Evaluation — 完整清楚評估

## 0. TL;DR

> chr19 head-to-head 用同樣 ISM × 2 flag protocol 跑 V3F BAM (`pononly_v3_fixed`) 與 V5 BAM (`threshold_compare/v5_flag`)，揭露：(1) **V5 Layer 1.5 把 reads 集中到 hp=1-1** — V5 hp=1-1 比 V3F 多 10,225 reads (+190%)，hp=3 (somatic ambiguous) 從 V3F 2,625 → V5 325 (-88%)；(2) **V5 hp=1-1 vs hp=2-1 = 1.86:1 偏 HP1**（priority bug feature 化）vs V3F = 0.59:1 中性；(3) **region-level marker filter (NG≥3) V3F 略勝 V5**：V3F 抓到 489 regions (463 TP) vs V5 396 regions (366 TP)，coverage +23.5% TP regions；marker rate V3F 0.947 vs V5 0.924 (+0.023 pp)；(4) **51.7% region (397/768) 兩 BAM 給的 NG_off 不同**，主要是 V5 把 V3F NG=5 region 塌成 NG=3 或 NG=4。**結論：V3F BAM 在 region-level downstream marker engineering 略勝 V5 BAM；V5 仍可作 production baseline（caller F1 = V3F）但 germline-absent 區改回 V3F 設計值得 follow-up**。

## 1. 完整 V3F vs V5 evaluation matrix

| 維度 | V3F 行為 | V5 行為 | 哪個穩健 | 證據來源 |
|---|---|---|---|---|
| **caller F1 vs SEQC2 truth (HCC1395 0.93/0.6)** | 0.7166 / 0.6273 | 0.7166 / 0.6273 | **平手**（V5 不改 caller）| 5/8 §8.5.2 |
| **17.3:1 baseline 偏移修正（germline-existent 區）** | 100% (chr19 752 / 全基因組 34,855 victims 修正) | 100% (與 V3F 同) | **平手** | 5/8 §5/§6 |
| **germline-absent 區域 hp tag 行為（5,789 chr19 events）** | 全標 hp=33（保守 ambiguous）| 4.19:1 偏 HP1（priority bug 繼承）| **V3F 勝** ⭐ | 5/9 Step D paired cross-ref |
| **chr19 read-level hp=1-1 集中度** | 5,380 reads (12.4% / non-zero buckets) | 15,605 reads (37.0% / non-zero) | **V3F 中性，V5 集中** | V5 BAM head-to-head (本報告) |
| **chr19 read-level hp=1-1:hp=2-1 ratio** | 0.59:1（中性） | **1.86:1（偏 HP1）** | V3F 中性 | 本報告 §2.1 |
| **chr19 read-level hp=3 (somatic ambiguous)** | 2,625 reads (保留 ambiguity) | 325 reads (Layer 1.5 改派為 hp=11/21) | V3F 保守 | 本報告 §2.1 |
| **region-level NG distribution** | 平均 NG 較高（NG=5: 122 regions）| 平均 NG 較低（NG=5: 32 regions, -74%）| V3F 多樣性高 | 本報告 §2.2 |
| **marker filter coverage (NG_off ≥ 3)** | 489 regions (463 TP) | 396 regions (366 TP) — coverage **-19%** | **V3F 勝** ⭐ | 本報告 §3 |
| **marker filter TP rate (flag=off)** | 0.947 | 0.924 (-0.023) | **V3F 略勝** | 本報告 §3 |
| **marker filter TP rate (flag=on, NG_on=2)** | 0.915 | 0.885 (-0.030) | **V3F 略勝** | 本報告 §3 |
| **marker decision gate (≥0.85)** | 兩者都通過 ✓ | 兩者都通過 ✓ | **平手**（marker 真實性確認）| 本報告 §3 |
| **per-region NG_off agreement (V3F vs V5)** | n/a | 397/768 (51.7%) regions 兩 BAM 給不同 NG | 結構性差異存在 | 本報告 §4 |
| **AMB% reduction（V5 設計目標）** | 17.5% | 8.0% (-54%) | **V5 設計目標達成** | 5/8 §8.5.1 |
| **HP:i:33 reduction（V5 設計目標）** | 239,679 | 110,197 (-54%) | **V5 設計目標達成** | 5/8 §8.5.1 |
| **20 No-regression indicators** | baseline ≈ V3F | V5 vs baseline = 0 critical regression | **V5 通過 sanity** | 5/8 §8.5 |

## 2. Read-level 詳細數據（chr19 head-to-head）

### 2.1 hp tag distribution by BAM × flag

| hp value | V3F off | V3F on | V5 off | V5 on | V5-V3F off Δ | 詮釋 |
|---|---|---|---|---|---|---|
| `0` (unphased) | 6,477 | 23,644 | 3,680 | 27,987 | -2,797 | V5 把更多 reads 進 phased buckets |
| `1` (germline HP1) | 11,216 | 11,216 | 12,181 | 12,181 | +965 | V5 略多 germline HP1 tag |
| `1-1` (somatic on HP1) | 5,380 | 0 | **15,605** | 0 | **+10,225 (+190%)** | **V5 Layer 1.5 集中化** ⚠ |
| `2` (germline HP2) | 13,303 | 13,303 | 7,995 | 7,995 | -5,308 | V5 較少 germline HP2 tag (?) |
| `2-1` (somatic on HP2) | 9,162 | 0 | 8,377 | 0 | -785 | 略少 |
| `3` (somatic ambiguous) | 2,625 | 0 | **325** | 0 | **-2,300 (-88%)** | **V5 Layer 1.5 改派 hp=33** |
| **TOTAL** | 48,163 | 48,163 | 48,163 | 48,163 | 0 | 守恆 ✓ |

**機制詮釋**：
- V5 Layer 1.5 在 germline 票數平手或缺席時用 somatic vote 決方向
- 結果：V3F 標 hp=33 (保守) 的 reads，V5 大多改派 hp=11 (somatic HP1) — 因為 self-phasing 機制下 sub-clone somatic 100% 共現偏向同邊
- **chr19 hp=1-1:hp=2-1 ratio**：
  - V3F: 5,380 / 9,162 = **0.59:1**（中性偏 HP2）
  - V5:  15,605 / 8,377 = **1.86:1**（偏 HP1，priority bug feature 化）
  - 對比 baseline 4.19:1（5/9 Step D 量化的 priority bug 主峰）

→ **V5 Layer 1.5 的 read-level 偏移在全 chr19 (含 germline-existent 區) 也存在**，不只 germline-absent 區。

### 2.2 Region-level NG distribution by BAM

| NG_off | V3F TP+FP regions | V5 TP+FP regions | V5/V3F ratio | 變化 |
|---|---|---|---|---|
| 0 | 1 | 0 | 0 | -1 |
| 1 | 54 | 96 | 1.78× | +42 |
| 2 | 154 | 281 | 1.82× | +127 |
| 3 | 228 | 184 | 0.81× | -44 |
| 4 | 113 | 75 | 0.66× | -38 |
| 5 | **122** | **32** | **0.26×** | **-90** |

→ **V5 BAM 把 NG distribution 整體往低 NG 移**（NG=5 急減 -74%）— Layer 1.5 把原本 V3F 上散在 hp=33 的 reads 集中到 hp=1-1，導致 bucket 多樣性下降。

## 3. Marker filter (HPFineNGroups) — V3F vs V5 直接比較

### 3.1 Marker (NG_off ≥ 3) coverage 與 TP rate

| BAM | regions | TP | FP | TP rate | flag=on NG_on=2 對應 cell |
|---|---|---|---|---|---|
| **V3F** | 489 | 463 | 26 | **0.947** | TP=367 / FP=34 / rate=**0.915** |
| **V5** | 396 | 366 | 30 | **0.924** | TP=270 / FP=35 / rate=**0.885** |
| V3F - V5 | +93 (+23.5%) | +97 | -4 | +0.023 | +97 / -1 / +0.030 |

**結論**：
- **V3F BAM marker coverage 比 V5 多 23.5% TP regions**（+97 TP regions 抓到）
- V3F marker rate 0.947 略勝 V5 0.924（+0.023 pp）
- 兩者都 ≥ 0.85 gate（marker 真實性都確認）
- V5 BAM 在 region-level marker engineering 上**輕微劣於 V3F BAM**

### 3.2 為什麼 V5 marker coverage 較低

V5 Layer 1.5 把 reads 集中 hp=1-1 → 同 region 內 bucket 多樣性下降 → 較多 region 從 NG=5/4 塌成 NG=3/2/1：
- V3F NG=5 (122) → V5 NG=5 (20) + NG=4 (56) + NG=3 (39) + NG=2 (8)
- V3F NG=4 (113) → V5 NG=4 (51) + NG=3 (51) + NG=2 (9)

→ V5 把 V3F 高 NG region 拆散，導致 marker filter (NG≥3) 漏抓部分 TP 案例。

### 3.3 Per-cell TP rate（哪個 BAM 信號更純？）

| Cell | V3F TP/FP/rate | V5 TP/FP/rate |
|---|---|---|
| NG_off=5 → NG_on=2（最強）| 122/1/**0.992** | 31/1/**0.969** |
| NG_off=4 → NG_on=2 | 93/6/0.939 | 103/5/**0.954** ⬆ |
| NG_off=4 → NG_on=1 | 20/0/1.000 | 27/0/1.000 |
| NG_off=3 → NG_on=2 | 110/15/0.880 | 130/22/0.855 |
| NG_off=3 → NG_on=1 | 115/4/0.966 | 73/2/**0.973** ⬆ |
| NG_off=2 → NG_on=1 | 96/26/0.787 | 215/36/**0.857** ⬆ |
| NG_off=2 → NG_on=2 | 42/12/0.778 | 6/7/0.462 ⬇ |
| NG_off=1 → NG_on=1 | 45/32/0.584 | 9/23/0.281 ⬇ |

→ V5 在低 NG cell（NG_off=2 → NG_on=2 / NG_off=1 → NG_on=1）TP rate 較低（更多 FP 流入）；高 NG cell 與 V3F 差不多。整體可解讀為：**V5 把高 NG TP 集中到中 NG（NG=4）的同時也讓低 NG cell 更 noisy**。

## 4. 51.7% region NG_off disagreement — 實質意義

768 regions 中 397 個（51.7%）兩 BAM 的 NG_off 不同：

| 主要 disagreement pattern | count | 詮釋 |
|---|---|---|
| V3F=3 → V5=2 | 105 | V5 集中化讓 NG 降 |
| V3F=4 → V5=3 | 51 | 同上 |
| V3F=5 → V5=4 | 56 | 同上 |
| V3F=2 → V5=3 | 27 | 反向 — V5 在某些區 NG 增 (Layer 1.5 把 hp=0 reads 派到 hp=11/21) |
| V3F=2 → V5=1 | 28 | V5 把 reads 集中到單 bucket |
| V3F=3 → V5=4 | 25 | 反向 |

→ **V3F vs V5 給的 region-level NG signal 結構性不同**，不只是 read-level 偏移；51.7% 比例不算邊緣。

→ **論文 implication**：若 HPFineNGroups marker 用於發表生物學結論，必須標註 BAM 來源（V3F vs V5），因為兩者給的 NG signal 在過半數 region 上不同。

## 5. caller F1 vs region-level marker 的 decoupling

關鍵發現（與 5/8 §8.5.2 互補）：

| 指標 | V3F | V5 | 差異 |
|---|---|---|---|
| caller F1 (HCC1395 0.93) | 0.7166 | 0.7166 | 0（完全相同）|
| caller F1 (HCC1395 0.6) | 0.6273 | 0.6273 | 0（完全相同）|
| ISM HPFineNGroups marker rate (chr19) | 0.947 | 0.924 | -0.023 |
| ISM marker coverage (chr19) | 489 regions | 396 regions | -19% |

→ **V3F vs V5 在 caller 層級無差異，但在 region-level downstream feature engineering 有量化差異**。
→ 對 production callers 而言 V5 = V3F；對 ISM downstream 研究而言 V3F BAM 略佳。

## 6. 結論與建議

### 6.1 「V3F 比 V5 有效」的清楚評估結論

**分區回答**：

| 區域 / 用途 | 結論 |
|---|---|
| **germline-existent 區 priority bug 修正** | V3F = V5（兩者都 100%）✓ |
| **germline-absent 區 hp tag 穩健性** | **V3F 勝**（V3F hp=33 保守 vs V5 4.19:1 priority bug 繼承）|
| **caller F1 vs truth** | V3F = V5（相同 0.7166/0.6273）✓ |
| **chr19 read-level hp=1-1 集中度** | V5 Layer 1.5 集中（+190% reads）— 設計選擇而非 bug |
| **chr19 region-level marker coverage** | **V3F 勝**（+23.5% TP regions 抓到, +0.023 marker rate）|
| **chr19 region-level marker decision gate ≥0.85** | V3F = V5（兩者都通過）|
| **AMB% / HP:i:33 reduction** | V5 達成設計目標（V3F 沒此設計）|
| **20 no-regression indicators** | V5 vs baseline 0 critical regression |

### 6.2 整體判斷

**「V3F 嚴格比 V5 更有效」= 否**：
- caller F1 平手
- germline-existent 區 priority bug 修正都 100%
- region-level marker 兩者都通過 0.85 gate

**「V3F 在某些區域比 V5 穩健」= 是**（已量化）：
- germline-absent 區 V5 繼承 priority bug 4.19:1 偏移
- region-level marker coverage V3F 多 23.5% TP regions
- region-level marker rate V3F 略高 +0.023 pp

**「V3F 應取代 V5 作 production」= 否**：
- V5 達成設計目標（AMB% -54% / HP:i:33 -54%）
- V5 通過 sanity 15/15 + no-regression 20 indicators
- V5 caller F1 與 V3F 完全相同

### 6.3 行動建議

**短期（不阻擋 production）**：
- V5 仍作 production tag baseline（caller F1 不變、sanity 通過）
- HPFineNGroups marker 發表時若用 ISM downstream feature，**標註 BAM 來源**（V3F vs V5）以避免 51.7% region NG signal 結構性差異被忽略

**中期（建議 follow-up）**：
- **V6-D binary patch**：V5 主結構 + germline-absent 區改回 V3F 行為（hp=33 conservative）→ 預期 read-level 4.19:1 偏移消除 + region-level marker coverage 接近 V3F + 保留 V5 AMB% reduction 設計目標
- 啟動條件：用戶批准 longphase-to-mod 第 6 commit
- ROI：read-level audit + region-level marker engineering 雙贏；估時 ~1 day binary patch + ~3 day 7 樣本驗證

**長期**：
- Phase C 7 樣本 V5 BAM × 完整 marker filter（NG ∧ AF<0.4 ∧ NR≥80 ∧ NonLOH）
- 比較 V3F BAM vs V5 BAM 7 樣本級別 marker robustness（本 chr19 結論的擴展）

## 7. Caveat

| Caveat | 影響 | 補強 |
|---|---|---|
| chr19 only（占 2.16%）| 其他 chr 行為未驗 | Phase C 7 樣本擴展 |
| V3F BAM 沒 caller VCF re-call | 用 V5 caller VCF 跑 V3F BAM 可能不完全公平 | 兩 BAM 都用同一 ClairS-TO VCF（caller 在 BAM 之前），所以公平 |
| HPFineNGroups marker 只測 NG 維度 | 完整 marker = NG ∧ AF<0.4 ∧ NR≥80 ∧ NonLOH | Phase C joint filter |
| BAM 存儲路徑跨 disk（V3F big7 / V5 big7） | I/O 時間可能影響 | 兩者都從同一 disk 跑，無 bias |

## 8. 引用 / Reproducibility

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
bash research/paired_priority_bug_audit/scripts/run_v6c_phaseB_chr19.sh         # V3F BAM
bash research/paired_priority_bug_audit/scripts/run_v6c_phaseB_chr19_v5bam.sh   # V5 BAM
python3 research/paired_priority_bug_audit/scripts/v6c_phaseB_v3f_vs_v5_compare.py
```

完整數據：
- V3F BAM 結果：`InterSubMod/research/paired_priority_bug_audit/v6c_phaseB_runs/`
- V5 BAM 結果：`InterSubMod/research/paired_priority_bug_audit/v6c_phaseB_v5bam_runs/`
- 比較 summary：`InterSubMod/research/paired_priority_bug_audit/v6c_phaseB_v5bam_runs/v3f_vs_v5_summary.tsv`
- 完整 log：`InterSubMod/research/paired_priority_bug_audit/v6c_phaseB_v5bam_runs/v3f_vs_v5_output.log`

關聯文件：
- 5/8 主報告：`InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`
- 5/9 PI 報告 errata：`InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md`
- 5/10 V6-C Phase B chr19 (V3F BAM)：`InterSubMod/research/paired_priority_bug_audit/05_V6C_phaseB_findings.md`
- 5/9 paired audit Step D：`InterSubMod/research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md`
