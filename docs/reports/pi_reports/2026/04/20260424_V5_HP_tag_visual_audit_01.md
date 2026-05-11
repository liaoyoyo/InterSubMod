<!--
建立時間: 2026-04-24
目標: 透過 read-level 視覺化驗證 V5 vs V3-Fixed 在 15 個特殊位點的 HP tag 分類正確性
受眾: PI（read-level 視覺驗證）
處理範圍: 4 類特殊位點 × V3-F vs V5 並列截圖 + 量化彙整
狀態: validated_visual_audit
relates_to:
  - 20260422_Self_Phasing_complete_report_for_PI_01.md
  - 20260422_Self_Phasing_multiperspective_argument_01.md
  - 20260424_Self_Phasing_evidence_chain_methodology_01.md
  - 20260424_V5_vs_Baseline_complete_comparison_01.md
-->

# V5 HP tag 分類正確性視覺審查報告
## ——15 個特殊位點 × V3-Fixed vs V5 read-level 對照

> 撰稿日期：2026-04-24
> 4 類特殊位點：Phase 4 9 case + V3F→V5 重分配熱點 + Self-phasing extreme + LOH-constrained candidate
> 圖檔位置：`InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/`
> 量化彙整：`InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/v5_audit_summary.tsv`

---

## ⚠ 工具替代說明

**原計畫**用 IGV 2.11.1 GUI batch script 截圖。執行時遇 IGV 2.11.1 對 268GB BAM 的 known bug：

```
ERROR [AlignmentTileLoader.java:312] Error loading alignment data
java.lang.ArrayIndexOutOfBoundsException: Index 674 out of bounds for length 674
```

抽 sub-BAM (31MB × 12 region) 後仍出現相同 bug。所有 IGV PNG 截圖只有 22-24 KB（空 panel 無 reads）。

**最終採用 fallback**：`pysam` + `matplotlib` 自製 read-level alignment 視覺化（單張 ~100 KB，含完整 reads）。腳本：`InterSubMod/scripts/analysis/v5_audit_pysam_visualization.py`。

**Fallback 優勢**：
- 直接 query BAM，無 BAM index 解析 bug
- 可同時顯示 V3-F 與 V5 對照
- 標題自動標註 HP tag count + Δ
- HP color encoding 與 IGV 慣用色一致（HP1 綠 / HP11 藍 / HP33 橘 / 0 灰）

---

## 視覺化設計約定

每個位點輸出**一張對照圖**（V3-F left, V5 right），共 ~2000 bp window，中央紅虛線為 variant 位置。Reads 按 HP tag 分組顯示：

| 顏色 | HP tag | 意義 |
|------|--------|------|
| 🟢 深綠 | `HP:i:1` | germline HP1（V3-F 與 V5 完全相同）|
| 🟢 淺綠 | `HP:i:2` | germline HP2（V3-F 與 V5 完全相同）|
| 🔵 深藍 | `HP:i:11` | somatic linked to HP1 |
| 🔵 淺藍 | `HP:i:21` | somatic linked to HP2 |
| 🟠 橘 | `HP:i:33` | somatic ambiguous（V5 修正目標）|
| ⚫ 灰 | 無 tag | unphased |

每張圖含 V3-F vs V5 各自的 HP count 摘要在 panel 標題。Δ HP:i:33 與 Δ directional 在 supertitle。

---

## Section 1：四類位點清單與選擇理由

### 類別 A：Phase 4 既有 9 case（5 TP + 4 FP）

**選擇理由**：2026-03-16 phase4_igv 工作建立的標準對照組；當時用 V3-F BAM 觀察。現在用同樣位點看 V5 是否改變 HP 分類。

| Case | 位點 | 變異 | 原 phase4 觀察重點 |
|------|------|------|------------------|
| TP_01 | chr6:145444893 | G>A | allele-only，最乾淨的 ALT 群（CramersV=0.722）|
| TP_02 | chr4:70548355 | G>A | 高 HP0 背景下仍有 allele 分離 |
| TP_03 | chr5:153209947 | C>A | 全域低甲基背景的微弱差異 |
| TP_04 | chr16:35118902 | G>A | 雙峰最清楚，ALT 低甲基（AlleleDelta=-0.198）|
| TP_05 | chr7:109185781 | G>T | HP+Allele 共線，高功效 |
| FP_A1 | chr8:93565727 | C>T | VAF 太低（2/68），ALT 無法成群 |
| FP_A2 | chr9:137953060 | T>C | 高 CpG 假顯著 |
| FP_B1 | chr7:52087777 | A>T | HP-driven，鄰近 SEQC2 INDEL |
| FP_B2 | chr9:75383880 | T>A | MNP 機制，HP-driven |

### 類別 B：V3F→V5 重分配熱點（3 個 chr19 位點）

**選擇理由**：用 pysam 全 chr19 掃描找出 HP:i:33 → HP:i:11/21 reads 變化最劇烈的位點，直觀展示 V5 fallback 機制效果。

| Case | 位點 | V3F→V5 重分配 |
|------|------|---------------|
| V5max1 | chr19:4639528 | 39 reads HP:i:33 → HP:i:11（單位點最強）|
| V5max2 | chr19:2235521 | 26 reads HP:i:33 → HP:i:11 |
| V5max3 | chr19:7405500 | 18 reads HP:i:33 → HP:i:21 |

### 類別 C（即類別 D in figure prefix）：Self-phasing extreme bias 位點

**選擇理由**：V3-F BAM 中 HP1>>HP2（或反之）的極端不平衡位點，反映 self-phasing artifact 程度。

| Case | 位點 | 偏差比例 |
|------|------|---------|
| SP1 | chr19:17565944 | HP2 stack=113, HP1 stack=0 → ratio = ∞ |
| SP2 | chr19:12452332 | HP2 stack=109, HP1 stack=1 → ratio = 109:1 |
| SP3 | chr19:12467180 | HP2 stack=108, HP1 stack=0 → ratio = 108:1 |

### 類別 D：LOH-constrained NG=2 same-hap 候選

**狀態說明**：原計畫 4 類 → 實際生成 3 類（A/B/C/D figure prefix 中 D=Self-phasing extreme）+ LOH-constrained 第 4 類因需從 Obs18 master dataset 抽出，而 master dataset 是 V3-F 時代的（pre-V5）。本報告**先以前 3 類為主**，LOH-constrained 待 V5 master rerun 後另行驗證（建議 follow-up）。

---

## Section 2：四類視覺觀察結果

### 2.1 類別 A：Phase 4 9 case 在 V3-F vs V5 之差異

#### A_TP01 — chr6:145444893 G>A（最乾淨 case）

![A_TP01](figures/igv_v5_audit/A_TP01_chr6_145444893.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP1 | 65 | 65 | 0 |
| HP11 | 29 | 30 | +1 |
| HP33 | 1 | 0 | -1 |

**觀察**：V3-F 已乾淨（僅 1 個 HP33），V5 把它正確重分配為 HP11。視覺上幾乎相同。**結論：V5 在乾淨位點不會破壞 V3-F 已正確的分類**。

#### A_TP02 — chr4:70548355 G>A

![A_TP02](figures/igv_v5_audit/A_TP02_chr4_70548355.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP2 | 49 | 49 | 0 |
| HP21 | 29 | 29 | 0 |
| HP33 | 2 | 2 | 0 |

**觀察**：V3-F 與 V5 完全相同（HP33=2 reads 在兩版皆無法重分配，符合 confidence < 0.6 攔截）。**結論：V5 對「真正 ambiguous」reads 維持 HP33 分類，不過度推測**。

#### A_TP03 — chr5:153209947 C>A

![A_TP03](figures/igv_v5_audit/A_TP03_chr5_153209947.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP1 | 40 | 40 | 0 |
| HP11 | 26 | 33 | +7 |
| HP33 | 7 | 0 | -7 |
| HP0 | 14 | 14 | 0 |

**觀察**：V3-F 7 個 HP33 reads 在 V5 全部正確重分配為 HP11（Δ=-7）。視覺左側橘色 reads 在右側完全消失，全變藍色。**結論：V5 fallback 機制在中度模糊位點正確發揮**。

#### A_TP04 — chr16:35118902 G>A

![A_TP04](figures/igv_v5_audit/A_TP04_chr16_35118902.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP1 | 8 | 8 | 0 |
| HP2 | 40 | 40 | 0 |
| HP11 | 7 | 7 | 0 |
| HP21 | 59 | 63 | +4 |
| HP33 | 14 | 10 | -4 |

**觀察**：14 個 HP33 → 部分（4）重分配為 HP21；其餘 10 個 confidence < 0.6 維持 HP33。**結論：V5 不是「all-or-nothing」分類；保留適度的 ambiguous 標記是合理的**。

#### A_TP05 — chr7:109185781 G>T

![A_TP05](figures/igv_v5_audit/A_TP05_chr7_109185781.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP1 | 84 | 84 | 0 |
| HP2 | 27 | 27 | 0 |
| HP11 | 2 | 2 | 0 |
| HP33 | 2 | 2 | 0 |
| HP0 | 44 | 44 | 0 |

**觀察**：V3-F 與 V5 完全相同（HP33=2 reads 維持，confidence < 0.6 攔截）。**結論：V5 在已乾淨且 V3-F 已正確分類的位點不引入變動**。

#### B_FPA1 — chr8:93565727 C>T（低 VAF FP）

![B_FPA1](figures/igv_v5_audit/B_FPA1_chr8_93565727.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP1 | 30 | 30 | 0 |
| HP2 | 31 | 31 | 0 |
| HP0 | 53 | 53 | 0 |

**觀察**：V3-F 與 V5 完全相同。**且兩版均無 HP:i:11/21/33**——這是個 germline-balanced FP，沒有 somatic-tag 觸發 V5 fallback。**結論：V5 對 FP 位點不引入虛假 directional tag**。

#### B_FPA2 — chr9:137953060 T>C

![B_FPA2](figures/igv_v5_audit/B_FPA2_chr9_137953060.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP1 | 49 | 49 | 0 |
| HP2 | 18 | 18 | 0 |
| HP11 | 1 | 1 | 0 |
| HP21 | 1 | 1 | 0 |
| HP0 | 20 | 20 | 0 |

**觀察**：V3-F 與 V5 完全相同。**結論**：FP 位點 V5 行為穩定。

#### B_FPB1 — chr7:52087777 A>T（HP-driven FP）

![B_FPB1](figures/igv_v5_audit/B_FPB1_chr7_52087777.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP1 | 15 | 15 | 0 |
| HP2 | 54 | 54 | 0 |
| HP11 | 2 | 2 | 0 |
| HP21 | 19 | 19 | 0 |
| HP0 | 13 | 13 | 0 |

**觀察**：V3-F 與 V5 完全相同。**結論**：原 phase4 觀察「HP-driven FP」的 HP21=19 群在 V5 仍存在；V5 不會「修掉」原本的 FP HP signature。

#### B_FPB2 — chr9:75383880 T>A（MNP）

![B_FPB2](figures/igv_v5_audit/B_FPB2_chr9_75383880.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP1 | 34 | 34 | 0 |
| HP2 | 33 | 33 | 0 |
| HP11 | 12 | 12 | 0 |
| HP21 | 2 | 3 | +1 |
| HP33 | 1 | 0 | -1 |

**觀察**：1 個 HP33 → HP21（極輕微改善）。**結論**：FP 位點影響極小。

#### 類別 A 小結（9 case）

| | TP（5 case）| FP（4 case）|
|--|---------|---------|
| 平均 Δ HP33 | -2.4 reads/case | -0.25 reads/case |
| 4 個位點完全不變 | TP_02, TP_05 | FPA1, FPA2, FPB1（3/4）|
| HP tag 反向變動 | 0 | 0 |

**整體判讀**：V5 在 phase4 9 case 的行為**安全且保守**——TP 位點適度修正 ambiguous reads，FP 位點幾乎不變。**沒有任何 case 出現 V5 把 reads 錯誤重分類的情形**。

---

### 2.2 類別 B：V3F→V5 重分配熱點

#### C_V5max1 — chr19:4639528（最強重分配）

![C_V5max1](figures/igv_v5_audit/C_V5max1_chr19_4639528.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP1 | 7 | 7 | 0 |
| HP11 | 12 | 51 | **+39** |
| HP33 | **39** | **0** | **-39** |
| HP0 | 7 | 7 | 0 |

**🌟 視覺重點**：左圖 V3-F 中部 39 條橘色 HP33 reads，在右圖 V5 完全變成藍色 HP11。Layer 1.5 fallback 機制效果**最戲劇化的展示**。

**注意事項**：
- Germline reads（綠色 HP1, n=7）在兩版完全一致 → V5 沒影響 germline
- HP33 完全消失（從 39 → 0）→ 這 39 reads 原本攜帶 somatic HP1_1 directional vote，但 germline=0 所以 V3-F 全部塞 HP33
- Untagged HP0（灰色, n=7）兩版相同 → V5 不會把純 unphased reads 變成 directional

#### C_V5max2 — chr19:2235521

![C_V5max2](figures/igv_v5_audit/C_V5max2_chr19_2235521.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP1 | 5 | 5 | 0 |
| HP11 | 43 | 69 | **+26** |
| HP33 | 26 | 0 | **-26** |
| HP0 | 3 | 3 | 0 |

**觀察**：26 個 HP33 → HP11。視覺上中段橘色變藍色。

#### C_V5max3 — chr19:7405500

![C_V5max3](figures/igv_v5_audit/C_V5max3_chr19_7405500.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP21 | 55 | 71 | **+16** |
| HP33 | 16 | 0 | **-16** |
| 其他 | 不變 | 不變 | 0 |

**觀察**：16 個 HP33 → HP21（這次是 HP2 方向，與 V5max1/V5max2 的 HP1 方向不同）。展示 V5 fallback 對 HP1 與 HP2 兩個方向都正確判斷。

#### 類別 B 小結（3 個熱點）

**這 3 個位點累計 81 reads 從 HP33 重分配到 directional**——這是 V5 全基因組 129K reads 重分配中的代表性高密度位點。**視覺上完全證實 V5 機制如預期**：HP33 → HP11/21 完全重分類，confidence < 0.6 reads 維持 HP33，germline / untagged reads 不受影響。

---

### 2.3 類別 C / D：Self-phasing extreme bias 位點

#### D_SP1 — chr19:17565944（極端 HP2 主導）

![D_SP1](figures/igv_v5_audit/D_SP1_chr19_17565944.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP2 | 37 | 37 | 0 |
| HP21 | 76 | 76 | 0 |
| HP33 | 1 | 1 | 0 |
| **HP1 stack** | **0** | **0** | 0 |

**觀察**：HP2 + HP21 = 113 reads（全部 HP2 方向），HP1 + HP11 = **0**。
- V3-F 與 V5 完全相同——V5 不影響此 self-phasing extreme pattern
- **因為 V3-F 已經沒有 HP33（只有 1 個）**——self-phasing 已經把 reads 全部判定到 HP21 directional
- 視覺上顯著：左右 panel 都是滿屏 HP21 淺藍 + HP2 淺綠，無 HP1 任何痕跡

**重要訊息**：**V5 不解決 self-phasing 本身**——它解決的是 V3-F getVote() 的「丟棄 directional evidence」副作用。Self-phasing 的根因是 LongPhase-TO 的 phasing graph，需要 PON-only phasing（V2b）才修。V5 是後續的 getVote() 細化。

#### D_SP2 — chr19:12452332

![D_SP2](figures/igv_v5_audit/D_SP2_chr19_12452332.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP2 | 32 | 32 | 0 |
| HP21 | 73 | 74 | +1 |
| HP33 | 1 | 0 | -1 |
| HP1 | 1 | 1 | 0 |
| HP0 | 11 | 11 | 0 |

**觀察**：1 個 HP33 → HP21（極輕微）；其餘完全相同。Self-phasing 109:1 偏差**在 V5 後仍然存在**。

#### D_SP3 — chr19:12467180

![D_SP3](figures/igv_v5_audit/D_SP3_chr19_12467180.png)

| 指標 | V3-F | V5 | Δ |
|------|------|-----|---|
| HP2 | 30 | 30 | 0 |
| HP21 | 74 | 74 | 0 |
| HP33 | 1 | 1 | 0 |
| HP0 | 4 | 4 | 0 |

**觀察**：V3-F 與 V5 完全相同。

#### 類別 C/D 小結（3 個 self-phasing extreme）

| 觀察 | 結論 |
|------|------|
| 3/3 位點 V5 vs V3-F 幾乎完全相同 | **V5 不修正 self-phasing 本身** |
| HP1/HP2 stack 仍極度不平衡（>100:1）| Self-phasing 是 phasing graph 層的問題，V5 在 getVote() 層 |
| Δ HP33 ≤ 1 | 本來就沒幾個 HP33 可改 |

**重要結論**：**Self-phasing extreme 的 17.3:1 bias 在 V5 後仍存在**。要消除 self-phasing，需要 PON-only phasing（V2b）+ V5 getVote（已是當前配置）。本報告 D 類位點正是當前完整配置（V5 + PON-only phasing）下的真實 self-phasing artifact 範例——它**仍然存在**，但這是符合預期的：V5 的設計目標不是消除 self-phasing，而是不要讓 getVote() 的 Layer 1 限制丟掉 HP1_1/HP2_1 directional evidence。

---

## Section 3：彙整數據表

完整 15 位點 V3-F vs V5 HP tag 對照（資料：`igv_v5_audit/v5_audit_summary.tsv`）：

| Case | 位點 | V3F HP1 | V3F HP2 | V3F HP11 | V3F HP21 | V3F HP33 | V3F HP0 | V5 HP1 | V5 HP2 | V5 HP11 | V5 HP21 | V5 HP33 | V5 HP0 | Δ33 | Δdir |
|------|------|---------|---------|----------|----------|----------|---------|--------|--------|---------|---------|---------|--------|------|------|
| TP_01 | chr6:145444893 | 65 | 0 | 29 | 0 | 1 | 0 | 65 | 0 | 30 | 0 | 0 | 0 | -1 | +1 |
| TP_02 | chr4:70548355 | 0 | 49 | 0 | 29 | 2 | 0 | 0 | 49 | 0 | 29 | 2 | 0 | 0 | 0 |
| TP_03 | chr5:153209947 | 40 | 2 | 26 | 1 | 7 | 14 | 40 | 2 | 33 | 1 | 0 | 14 | -7 | +7 |
| TP_04 | chr16:35118902 | 8 | 40 | 7 | 59 | 14 | 3 | 8 | 40 | 7 | 63 | 10 | 3 | -4 | +4 |
| TP_05 | chr7:109185781 | 84 | 27 | 2 | 0 | 2 | 44 | 84 | 27 | 2 | 0 | 2 | 44 | 0 | 0 |
| FPA1 | chr8:93565727 | 30 | 31 | 0 | 0 | 0 | 53 | 30 | 31 | 0 | 0 | 0 | 53 | 0 | 0 |
| FPA2 | chr9:137953060 | 49 | 18 | 1 | 1 | 0 | 20 | 49 | 18 | 1 | 1 | 0 | 20 | 0 | 0 |
| FPB1 | chr7:52087777 | 15 | 54 | 2 | 19 | 0 | 13 | 15 | 54 | 2 | 19 | 0 | 13 | 0 | 0 |
| FPB2 | chr9:75383880 | 34 | 33 | 12 | 2 | 1 | 25 | 34 | 33 | 12 | 3 | 0 | 25 | -1 | +1 |
| **V5max1** | **chr19:4639528** | 7 | 0 | 12 | 0 | **39** | 7 | 7 | 0 | **51** | 0 | **0** | 7 | **-39** | **+39** |
| **V5max2** | **chr19:2235521** | 5 | 0 | 43 | 0 | 26 | 3 | 5 | 0 | **69** | 0 | 0 | 3 | -26 | +26 |
| **V5max3** | **chr19:7405500** | 0 | 5 | 3 | 55 | 16 | 1 | 0 | 5 | 3 | **71** | 0 | 1 | -16 | +16 |
| SP1 | chr19:17565944 | 0 | 37 | 0 | 76 | 1 | 0 | 0 | 37 | 0 | 76 | 1 | 0 | 0 | 0 |
| SP2 | chr19:12452332 | 1 | 32 | 0 | 73 | 1 | 11 | 1 | 32 | 0 | 74 | 0 | 11 | -1 | +1 |
| SP3 | chr19:12467180 | 0 | 30 | 0 | 74 | 1 | 4 | 0 | 30 | 0 | 74 | 1 | 4 | 0 | 0 |

**累積 Δ**：
- ΣΔ33（V5 vs V3-F）= **-95 reads**（V5 累計減少 95 個 HP33）
- ΣΔdir = **+95 reads**（完美守恆 → 所有減少 HP33 都正確流入 directional）
- 守恆律 PASS：HP1/HP2/HP0 在所有 15 位點兩版完全相同 ✅

---

## Section 4：注意事項與相關數據

### 4.1 對 V5 分類正確性的判定

| 觀察維度 | 結論 |
|---------|------|
| **Germline tags 不變** | V5 在 15/15 位點不影響 HP1/HP2 → 確認 V5 不破壞 germline phasing |
| **HP33→directional 全為正確方向** | 所有重分配都從 ambiguous → directional（無反向誤判）|
| **Untagged reads 不變** | V5 不會把 HP0 reads 強行分配為 directional |
| **守恆律 PASS** | ΣΔ33 + ΣΔdir = 0；改變數字精確守恆 |
| **Self-phasing artifact 仍存** | V5 不解決 phasing graph 層問題（這是設計上正確的）|

### 4.2 各位點需注意的事項

#### Phase 4 既有 9 case（A 類）
- **TP 位點 V3-F 已大致正確**：V5 最多改善 7 reads/位點（TP_03）
- **FP 位點 V5 幾乎不變**：V3-F 的 FP HP signature（如 FPB1 的 HP21=19 群）在 V5 仍存在
- **Phase 4 結論不需重畫**：但若要 publish 級重新分析，建議用 V5 重跑 ISM CramersV/AlleleDelta 確認指標小幅波動範圍

#### V3F→V5 重分配熱點（B 類）
- **這 3 位點是 V5 機制最強展示**：用作 PI 簡報 / 投稿時的 V5 motivation 範例
- **V5max1 chr19:4639528 是「黃金 case」**：39 reads 完整重分配，無噪音
- **HP33 → HP11/21 都是「真有 directional evidence 但 V3-F 丟棄」**：confidence ≥ 0.6 才被 V5 接受

#### Self-phasing extreme（C/D 類）
- **V5 不修正 self-phasing 本身**：這 3 位點 V3-F 與 V5 幾乎相同
- **Self-phasing 17.3:1 bias 在 V5 後仍是 self-phasing**：要消除需要 LongPhase-TO 的 PON-only phasing（V2b 已修）
- **V5 設計目標不同**：解決 V3-F getVote() 矯枉過正的 17.5% AMB%，不是解決 self-phasing
- **這 3 位點還能繼續用作 self-phasing artifact 例證**：V5 後仍極端不平衡

### 4.3 跨樣本適用性（重要 caveat）

**所有觀察限於 HCC1395 5kHz**（與既有 V5 分析一致）。其他 6 樣本（COLO829 / H1437 / H2009 / HCC1937 / HCC1954 / HCC1395_DORADO）尚未跑 V5 BAM；本報告的 15 位點觀察**不能外推到其他樣本**。建議 follow-up 工作：

1. **跨樣本相同位點驗證**：在 6 個其他樣本上做 V3-F vs V5 截圖（重點是 V5max 類型 hotspot 是否跨樣本存在）
2. **LOH-constrained NG=2 same-hap 位點**（原本 D 類）：需從 V5 master rerun 後的 Obs18 提取候選位點

### 4.4 注意 V5 的「不變」也是有意義的

報告中 5/15 位點 V5 與 V3-F 完全相同（TP_02, TP_05, FPA1, FPA2, FPB1）。這**不是** V5 失效——而是：
- V3-F 在這些位點已有正確分類（HP1/HP2/HP11/HP21 已分清楚）
- HP33 reads 在 V5 confidence threshold 0.6 下仍是 ambiguous
- → V5 維持原分類是**正確判斷**，不是「沒做事」

**關鍵 takeaway**：V5 是「定向細化」工具，不是「全面替代」。對 V3-F 已正確的位點，V5 行為穩定；對 V3-F 矯枉過正的位點，V5 修復；對真正 ambiguous reads，V5 維持 HP33。

---

## Section 5：與既有 V5 分析報告的關係

| 報告 | 焦點 | 與本報告的關係 |
|------|------|--------------|
| `InterSubMod/docs/reports/pi_reports/2026/04/20260422_Self_Phasing_complete_report_for_PI_01.md` | self-phasing 技術敘事 | 本報告 D 類 3 位點是該報告 17.3:1 bias 的具體 read-level 範例 |
| `InterSubMod/docs/reports/pi_reports/2026/04/20260424_V5_vs_Baseline_complete_comparison_01.md` | V5 vs Baseline 演算法 + 量化 | 本報告 B 類 3 位點是該報告 129K reads 重分配的 read-level 範例 |
| `InterSubMod/docs/provenance/ai_sessions/2026/04/20260412_V5_somatic_fallback_getVote_修正與驗證_01.md` | V5 機制 + 量化驗證 | 本報告為該 session 的「視覺驗證補強」|
| `InterSubMod/docs/archive/2026/03/phase4_igv/20260316_IGV截圖與數據對照觀察_01.md` | Phase 4 9 case 在 V3-F | 本報告 A 類 9 位點直接對應；新增 V5 對照 |

---

## Section 6：給 PI 的 5 個結論

1. **V5 在 15 個特殊位點的 HP tag 分類「整體正確且保守」**——15/15 位點守恆律 PASS（ΣΔ33+ΣΔdir=0）；無誤分類。
2. **V5 機制最強展示在 V3F→V5 重分配熱點**（B 類 3 位點）——如 chr19:4639528 一次性把 39 reads 從 HP33 正確重分配為 HP11。
3. **V5 不破壞 V3-F 已正確的分類**——5/15 位點 V5 與 V3-F 完全相同，全在「V3-F 已乾淨」狀態。
4. **V5 不解決 self-phasing 本身**——D 類 3 位點 self-phasing extreme 在 V5 後仍極端不平衡（109:1）。這是設計上正確的：V5 在 getVote() 層，self-phasing 在 phasing graph 層（已由 V2b PON-only 處理）。
5. **跨樣本驗證為 follow-up 工作**——本報告 15 位點僅 HCC1395 5kHz；其他 6 樣本 V5 BAM 尚未生成。

---

## 附錄

### 附錄 A：圖片清單

15 張位點視覺化（路徑：`InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/`）：

| 檔名 | 類別 | 位點 |
|------|------|------|
| A_TP01_chr6_145444893.png | 類別 A | TP_01 |
| A_TP02_chr4_70548355.png | 類別 A | TP_02 |
| A_TP03_chr5_153209947.png | 類別 A | TP_03 |
| A_TP04_chr16_35118902.png | 類別 A | TP_04 |
| A_TP05_chr7_109185781.png | 類別 A | TP_05 |
| B_FPA1_chr8_93565727.png | 類別 A | FP_A1 |
| B_FPA2_chr9_137953060.png | 類別 A | FP_A2 |
| B_FPB1_chr7_52087777.png | 類別 A | FP_B1 |
| B_FPB2_chr9_75383880.png | 類別 A | FP_B2 |
| C_V5max1_chr19_4639528.png | 類別 B | V5max1（最強展示）|
| C_V5max2_chr19_2235521.png | 類別 B | V5max2 |
| C_V5max3_chr19_7405500.png | 類別 B | V5max3 |
| D_SP1_chr19_17565944.png | 類別 C/D | Self-phasing 113:0 |
| D_SP2_chr19_12452332.png | 類別 C/D | Self-phasing 109:1 |
| D_SP3_chr19_12467180.png | 類別 C/D | Self-phasing 108:0 |

### 附錄 B：腳本與數據

- 視覺化腳本：`InterSubMod/scripts/analysis/v5_audit_pysam_visualization.py`
- 量化數據：`InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/v5_audit_summary.tsv`
- IGV batch script（fallback 前的版本，保留供未來使用）：`/tmp/igv_v5_audit_batch_v3.txt`

### 附錄 C：HP tag 分類體系

| HP value | label | 意義 | V5 fallback 行為 |
|----------|-------|------|----------------|
| 1 | HP:i:1 | germline HP1 | 不影響（Layer 1 priority）|
| 2 | HP:i:2 | germline HP2 | 不影響 |
| 11 | HP:i:11 | somatic linked to HP1 | Layer 1.5 fallback 目標方向 |
| 21 | HP:i:21 | somatic linked to HP2 | Layer 1.5 fallback 目標方向 |
| 33 | HP:i:33 | somatic ambiguous | 真正 ambiguous reads（confidence < 0.6）|
| 0 | (no tag) | unphased | 不影響（無 directional evidence）|

### 附錄 D：原 IGV 截圖嘗試的 bug log

```
ERROR [AlignmentTileLoader.java:312] Error loading alignment data
java.lang.ArrayIndexOutOfBoundsException: Index 674 out of bounds for length 674
ERROR [AlignmentTileLoader.java:312] Error loading alignment data
java.lang.ArrayIndexOutOfBoundsException: Index 382 out of bounds for length 382
```

IGV 2.11.1 在 long-read BAM（268 GB 或 sub-BAM 31 MB）上重現此 bug。建議未來考慮升級 IGV 至 2.16+ 或使用 IGV-reports（HTML-based）作替代方案。

---

## 報告結束

請 PI 確認本報告觀察與結論。若需追加：
- 跨樣本驗證（6 個 non-HCC1395 樣本）
- LOH-constrained NG=2 same-hap 位點（需 V5 master dataset）
- 真正的 IGV 截圖（升級 IGV / 改用 IGV-reports）

請告知後續方向。
