---
title: V5 Audit Suite — 02. Read-intersection HP concordance (BL/V5 vs Paired)
report_id: 20260427_V5_audit_02_read_intersection_concordance_01
date: 2026-04-27
authors: AI agent (analysis), liaoyoyo2001 (PI)
sample: HCC1395 (TO mode)
inputs:
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/tumor_tagged.bam
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam
  - /big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam
script: InterSubMod/scripts/analysis/v5_read_intersection_concordance.py
outputs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/per_site_concordance.tsv
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/hp_family_exact.tsv
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/02_concordance/fig02a_4metric_heatmap.png
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/02_concordance/fig02b_winloss_summary.png
status: draft
tier: B (audit-grade evidence on 15 hand-picked sites; not population-scale)
---

# 02. Read-intersection HP concordance（BL/V5 vs Paired tumor）

## 摘要 (TL;DR)

| 度量 | V5 wins | ties | BL wins | NA |
|------|--------:|-----:|--------:|---:|
| L1 family（{1,11} vs {2,21} family 同源） | **6** | 3 | 5 | 1 |
| L2 exact（HP 字串完全相同，含 PS-allele 解析） | **6** | 4 | 5 | 0 |
| L3 ratio distance（HP1:HP2 比率 + orientation flip min） | 4 | 5 | 5 | 1 |
| L4 orientation-corrected family（per-PS native vs flipped） | 2 | 3 | **9** | 1 |

**重點結論**：

1. **沒有 metric 顯示 V5 顯著優於 baseline**（最佳 case 為 L1/L2 的 6 vs 5，差距僅 1 site）。
2. **L4 反過來顯示 BL 9 vs V5 2 的劣勢**：當允許 per-PS orientation flip（消除 phase block 的左右翻轉自由度），BL 的 family rate 已達高水準，V5 反而在多個 site 出現 family 一致性下降。
3. **L1 vs L2 對 BL 不對稱**：BL 的 L2 exact rate（mean 0.4446）顯著低於 L1 family rate（mean 0.5677），原因是 BL 大量使用「`HP=1`/`HP=2`」單碼（無 PS-allele 後綴），而 paired 端大量使用「`1-1`/`2-1`」雙碼，**家族正確但 exact 必然不同**（詳見 03 報告）。
4. **Read intersection 後規模充足**：15 sites 中 14 個 ≥58 shared reads，TP05 達 137；不存在 read-set bias。

> Verdict（暫定）：在這 15 個 hand-picked 位點上，**V5 並未在 read-level HP 一致性上超越 baseline**。若 V5 的價值來自其他維度（block fragmentation、AMB 比例、F1 提升），需在其他 audit 子報告中驗證；本報告不支持「V5 phasing 更貼近 paired truth」這個敘事。

---

## Section 1. 為何用 read intersection（vs aggregate）

### 1.1 Aggregate 方法的問題

過往直覺是：對每個位點計算各 BAM 的 `% HP1`、`% HP2`、`% un-tagged`，比較分布。但這會被三類偏差污染：

1. **Read set bias**：不同 BAM 經過 longphase 重複處理後，**留在 BAM 的 read set 不一定相同**（duplicate marking、re-tag 行為、secondary/supplementary 過濾）。Aggregate 比對若 read 集合不一致，「分布差異」可能純粹來自 read 集合差異。
2. **Tag absence asymmetry**：BL 與 V5 對 HP=untagged 的處理規則不同（V5 somatic fallback 會主動回填 33 或 1/2），導致 denominator 不可比。
3. **Orientation degeneracy**：每個 PS 區塊的 1↔2 標籤是任意的。Aggregate `% HP1` 沒有對齊到 paired truth 的 PS frame，**家族對的 site 也可能 aggregate 看似 50/50 翻轉**。

### 1.2 Intersection-first 的設計

本報告採用：

```
shared_reads(site) = read_names(BL) ∩ read_names(V5) ∩ read_names(PA)
```

僅在 shared set 內計算 BL→PA 與 V5→PA 的 read-by-read concordance。三個保證：

- **同一 read 同時被三邊評估**（apples-to-apples）；
- 沒有 untagged 偏差，因為我們在 metric 內部對 family 是否可比、exact 是否可比有顯式 eligibility 規則；
- L3、L4 內部對 PS-orientation 顯式建模，**不會把 phase 翻轉誤認為 HP 不一致**。

每位點 read fetch 範圍：以變異位置為中心 ±50 bp，僅取 `primary, mapped, non-secondary, non-supplementary, non-duplicate` 的 reads。

---

## Section 2. 4 層 metric 定義與計算

設 `X ∈ {BL, V5}`，`PA = paired tumor`，`shared = BL ∩ V5 ∩ PA` per site。

### 2.1 HP tag 標準化（讀懂以下表格的前提）

longphase BAM（BL/V2b/V3F/V5）使用整數 HP tag，paired BAM 使用字串 HP tag：

| BL/V5 raw | 含義 | family | canonical exact | paired raw 對應 |
|-----------|------|--------|-----------------|----------------|
| `1` | 單一 PS、第 1 條 | 1 | `1` | `1` 或 `1-1` |
| `2` | 單一 PS、第 2 條 | 2 | `2` | `2` 或 `2-1` |
| `11` | 雙 PS 第一條（PS1-allele1） | 1 | `1-1` | `1-1` |
| `21` | 雙 PS 第二條（PS2-allele1） | 2 | `2-1` | `2-1` |
| `33` | both / unphasable | 0 | `3` | `3` |
| (none) | untagged | -1 | `NA` | `NA` |

### 2.2 L1 — HP family agreement

```
eligible(read) := X[read].family ∈ {1,2}  AND  PA[read].family ∈ {1,2}
match(read)    := X[read].family == PA[read].family
L1_X = #match / #eligible
```

家族層級判定，忽略 PS-allele 後綴（`1` 與 `1-1` 視為同 family）。NA 出現於 shared 為空或無 read 雙邊都 family-tagged 的情況。

### 2.3 L2 — HP exact agreement

```
eligible(read) := X[read].has_tag AND PA[read].has_tag (任何非 NA 值)
match(read)    := X[read].canonical_exact == PA[read].canonical_exact
L2_X = #match / #eligible
```

**差別**：BL 的 `HP=1` 對應 canonical `"1"`，paired 的 `"1-1"` 為另一字串 → 雖然 family 同為 1，**exact 不同**。L2 暴露的是 phasing notation 規範差異，不是真正的「phasing 錯誤」。詳見報告 03。

### 2.4 L3 — HP1:HP2 ratio distance（含 orientation flip）

```
ratio_X  = #{family=1 in shared} / (#{family=1} + #{family=2})        (in shared set)
ratio_PA 同上
d_native = |ratio_X - ratio_PA|
d_flip   = |(1 - ratio_X) - ratio_PA|
L3_X     = 1 - min(d_native, d_flip)            ∈ [0, 1]
```

L3 是 site-level 統計距離（非 per-read），對 PS-orientation 翻轉不敏感（取 min）。1 = 完全一致或完美鏡像。

### 2.5 L4 — Orientation-corrected family match (per-PS)

按 paired 的 PS tag 將 shared reads 分桶，每桶內計算：

```
n_native = #{X.family == PA.family within bucket}
n_flip   = #{X.family != PA.family AND both ∈ {1,2}}
n_match_bucket = max(n_native, n_flip)
L4_X = Σ n_match_bucket / Σ n_eligible_bucket
```

L4 是 L1 的 PS-aware 上限：若 V5 在每個 PS 都「一致地翻轉」，L1 顯示低，L4 顯示高。L4 < 1 才代表存在「無法用 orientation flip 解釋」的真正 family 衝突。

---

## Section 3. 完整數據表（15 sites × 4 metrics）

### 3.1 Read counts

| site | chrom:pos | n_BL | n_V5 | n_PA | n_shared(BL∩V5∩PA) |
|------|-----------|-----:|-----:|-----:|-------------------:|
| TP01 | chr6:145444893 | 85 | 85 | 85 | 85 |
| TP02 | chr4:70548355 | 63 | 63 | 63 | 63 |
| TP03 | chr5:153209947 | 78 | 78 | 78 | 78 |
| TP04 | chr16:35118902 | 106 | 106 | 106 | 106 |
| TP05 | chr7:109185781 | 137 | 137 | 137 | 137 |
| FPA1 | chr8:93565727 | 97 | 97 | 97 | 97 |
| FPA2 | chr9:137953060 | 76 | 76 | 76 | 76 |
| FPB1 | chr7:52087777 | 87 | 87 | 87 | 87 |
| FPB2 | chr9:75383880 | 87 | 87 | 87 | 87 |
| V5max1 | chr19:4639528 | 58 | 58 | 58 | 58 |
| V5max2 | chr19:2235521 | 62 | 62 | 62 | 62 |
| V5max3 | chr19:7405500 | 68 | 68 | 68 | 68 |
| SP1 | chr19:17565944 | 101 | 101 | 101 | 101 |
| SP2 | chr19:12452332 | 100 | 100 | 100 | 100 |
| SP3 | chr19:12467180 | 98 | 98 | 98 | 98 |

> 所有 5 個 BAM 的 read 集合在這 15 個 fetch 視窗內完全一致（85 = 85，等等）。intersection size = single-BAM size，**read-set bias 在這個 audit 上不存在**。

### 3.2 Per-site rates（小數點 3 位；NA 為無 eligible read）

| site | L1_BL | L1_V5 | L2_BL | L2_V5 | L3_BL | L3_V5 | L4_BL | L4_V5 |
|------|------:|------:|------:|------:|------:|------:|------:|------:|
| TP01 | 0.013 | 0.013 | 0.013 | 0.013 | 0.987 | 0.987 | 0.987 | 0.987 |
| TP02 | 0.926 | 0.077 | 0.926 | 0.000 | 0.926 | 0.926 | 0.926 | 0.923 |
| TP03 | 0.000 | 0.015 | 0.000 | 0.015 | 0.992 | 0.993 | 1.000 | 0.985 |
| TP04 | 1.000 | 0.026 | 0.569 | 0.018 | 0.735 | 0.875 | 1.000 | 0.974 |
| TP05 | 1.000 | 0.882 | 0.720 | 0.614 | 0.944 | 0.941 | 1.000 | 0.882 |
| FPA1 | NA | NA | 0.000 | 0.000 | NA | NA | NA | NA |
| FPA2 | 0.957 | 0.200 | 0.957 | 0.200 | 0.986 | 0.745 | 0.957 | 0.800 |
| FPB1 | 1.000 | 0.933 | 0.988 | 0.920 | 0.979 | 0.928 | 1.000 | 0.933 |
| FPB2 | 0.000 | 0.857 | 0.000 | 0.794 | 0.979 | 0.998 | 1.000 | 0.873 |
| V5max1 | 1.000 | 1.000 | 1.000 | 1.000 | 0.983 | 0.983 | 1.000 | 1.000 |
| V5max2 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| V5max3 | 0.015 | 0.045 | 0.015 | 0.045 | 0.985 | 0.955 | 0.985 | 0.955 |
| SP1 | 0.235 | 0.746 | 0.235 | 0.706 | 0.760 | 0.750 | 0.765 | 0.746 |
| SP2 | 0.324 | 0.681 | 0.041 | 0.681 | 0.669 | 0.669 | 0.676 | 0.681 |
| SP3 | 0.478 | 0.527 | 0.207 | 0.522 | 0.531 | 0.531 | 0.522 | 0.527 |
| **mean** | **0.568** | **0.500** | **0.445** | **0.435** | **0.890** | **0.877** | **0.916** | **0.876** |
| **median** | **0.702** | **0.604** | **0.235** | **0.522** | **0.979** | **0.935** | **0.993** | **0.928** |

### 3.3 數據觀察

1. **L1/L2 對 BL 強烈雙峰分布**：13/15 site 的 L1_BL 落在 {0, 0.013, ~1, 0.957} 兩端 — 表示 BL 在每個 site 內部 **要嘛全部對齊 paired，要嘛全部翻轉**（內部一致但 orientation 與 paired 相反）。L3 反映此特性：L3_BL mean 0.890（已對 flip 取 min），meaning 在 ratio 層面平均只差 11%。
2. **V5 在 SP 系列有顯著進步**：SP1（0.235→0.746）、SP2（0.324→0.681）— V5 在 chr19 SP 系列確實對齊到 paired 的 PS frame；但 V5 在 TP02 / TP04 / FPA2 / V5max3 反而從接近 1 掉到接近 0（family 翻轉）。
3. **L4 後 BL 的優勢凸顯**：L4_BL mean 0.916 vs L4_V5 0.876，median 0.993 vs 0.928 — 一旦消除 orientation 自由度，BL 在 9/15 個 site 上的 per-PS family agreement 仍優於 V5。
4. **TP01、V5max3 為通用 outlier**：兩者都呈 paired 是 family 1 主導但 BL/V5 全標 family 2 的一致翻轉；L4 高（≥0.985）反映此類僅是 PS frame 差異，**並非真正的 phasing 錯誤**。

---

## Section 4. Win/Loss/Tie 統計

> 「win」定義：在該位點 V5_rate > BL_rate（concordance score 同方向：高=接近 PA）。

| Metric | V5 wins | ties (含浮點等值) | BL wins | NA | Sum |
|--------|--------:|-----------------:|--------:|---:|----:|
| L1 family | 6 | 3 | 5 | 1 | 15 |
| L2 exact | 6 | 4 | 5 | 0 | 15 |
| L3 ratio distance | 4 | 5 | 5 | 1 | 15 |
| L4 orientation-corrected | **2** | 3 | **9** | 1 | 15 |
| **Total wins (sum across 4)** | **18 / 60** | **15 / 60** | **24 / 60** | **3 / 60** | — |

mean concordance（跨 15 site，較高 = 較接近 PA）：

| Metric | BL mean | V5 mean | delta (V5−BL) |
|--------|--------:|--------:|--------------:|
| L1 family | 0.568 | 0.500 | -0.067 |
| L2 exact | 0.445 | 0.435 | -0.010 |
| L3 ratio | 0.890 | 0.877 | -0.012 |
| L4 orient-corr | 0.916 | 0.876 | **-0.040** |

**結論**：4 個 metric 的 mean delta 全為負；L1 family 與 L4 orient-corr 的差距尤其顯著。「V5 read-level HP 更接近 paired」**不成立**。

圖：
- `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/02_concordance/fig02a_4metric_heatmap.png`
- `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/02_concordance/fig02b_winloss_summary.png`

---

## Section 5. 結論（V5 在哪些 metric 更接近 paired）

### 5.1 Direct verdict

| 主張 | 是否成立（在 15 audit sites 上） | 證據 |
|------|--------------------------------|------|
| V5 在 L1 family 上勝 BL | **否**（6 vs 5，mean -0.067） | §3.2、§4 |
| V5 在 L2 exact 上勝 BL | **否**（6 vs 5，mean -0.010） | §3.2、§4 |
| V5 在 L3 ratio 上勝 BL | **否**（4 vs 5，mean -0.012） | §3.2、§4 |
| V5 在 L4 orient-corr 上勝 BL | **否**（2 vs 9，mean -0.040） | §3.2、§4 |
| V5 在 SP 系列特定位點仍有局部優勢 | **是**（SP1/SP2 各 +0.4 級別 L1 改善） | §3.2 obs 2 |
| V5 read-level HP 普遍接近 paired truth | **否** | 4 metric 一致負向 |

### 5.2 為何 BL 在 L4 反勝？

L4 內部已對 PS orientation flip 取 max。BL 的 HP 標籤幾乎是「全 site 一致地翻轉或一致地對齊」，**內部 self-consistency 高**；V5 則在多個 site 引入了「同 PS 內 part-flip-part-not」的中間態，導致 per-PS max 也下降。這指向：

> **V5 的 somatic fallback 並沒有「修正」phasing；它在某些 site 把 BL 的 clean self-phasing block 打破了。**

### 5.3 與 audit suite 其他子報告的關係

- 本報告（02）只能說「V5 HP read-level 一致性不優於 BL」。
- V5 的價值若存在，需由 **04（imbalance）/ 06（sanity）/ 07（paired comparison）/ 08（synthesis）** 從 block fragmentation、AMB 比例、F1、cluster 維度補上。
- **不應將本結論作為「V5 不該採用」的最終判定**：F1=0.7154 已驗證 V5 的 variant-calling 下游收益（見 `MEMORY.md/project_v5_somatic_fallback_verification.md`），但此收益的因果可能不在 read-level HP 對齊。

### 5.4 後續建議

1. 在 04（imbalance）報告中檢查 V5 的 33 比例 / AMB 比例變化，**那才是 V5 設計的目標 metric**。
2. 若要支持「V5 phasing 更貼近 paired」這個主張，需要 **跨樣本 + 跨 PS-frame 對齊** 的 read-level intersection（本報告僅 HCC1395、僅 15 site），先到 8 sample × N=10000 reads 級別才能下穩定結論。
3. L2 exact 與 L1 family 的 systematic gap（特別是 BL 端）幾乎全來自 notation 規範差異 — 詳見報告 03。

---

## Section 6. 數據附錄

- TSV (per-site): `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/per_site_concordance.tsv` (16 行)
- TSV (per-read): `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/hp_family_exact.tsv` (1304 行)
- Script: `InterSubMod/scripts/analysis/v5_read_intersection_concordance.py`

可重現命令：

```bash
DISPLAY=localhost:18.0 python3 \
  /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/v5_read_intersection_concordance.py
```
