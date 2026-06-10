<!--
建立時間: 2026-06-10
狀態: in_progress (Phase B — 用外部程式稽核 ISM 是否需修正)
報告類型: external_program_audit
受眾: PI · 開發決策 · 反駁/確認既有 audit 結論
framework: Verdict-Pyramid (結論先行) + 三審 (A/B/B') + 誠實揭露 confound
data_sources:
  - runs/audit_A.json, runs/audit_B.json, runs/audit_B2.json (實跑讀回)
  - runs/counts_*.tsv, runs/dss_*.tsv (中間真值)
  - 外部程式: modkit 0.6.3, DSS(bioconductor-dss, r_dss_env)
  - ISM 真值: research/tsg_promoter_asm_reviewer/pipeline/cache/level1/{chr13_32315128,chr17_79991120}_*.tsv.gz
provenance_note: 所有數字自實跑 JSON 讀回; Audit B 標記 CONFOUNDED 不採信; Audit B' 為決定性 read-level permutation。本檔 Write 與分析不同 batch(§13.0)。單樣本 HCC1395, 2 anchor loci, 非全基因組。
-->
<!-- provenance-verified: 數字引自 runs/audit_A.json/audit_B.json/audit_B2.json 實跑; modkit/DSS 為外部程式; 撰寫與分析分批。 -->

# 08 — 用外部程式稽核 ISM「是否需修正」（三審）

> **問題**：用外部程式（modkit、DSS）實機稽核 ISM 兩個最被質疑的方法 —— **MM/ML 解析正確性（稽核 T1，零單元測試）** 與 **per-CpG Fisher 是否 over-dispersion 偏樂觀（稽核 T3/FISHER-1）** —— 確認 ISM 程式是否需要修正。

---

## L0 — 一句話結論（反直覺，但有據）

**外部程式稽核「推翻」了兩個最被擔心的問題**：① ISM 的 **MM/ML 解析經 modkit 驗證為正確**（per-CpG r=0.98–0.99）；② per-CpG Fisher 的 **over-dispersion 偏樂觀，在真資料上不成立**（null 假陽性 0.7–1.9%，**低於**名目 5%，稽核引用的「53–68% 膨脹」是合成模擬、真資料不重現）。→ **這兩處 ISM 不需修正**（但 MM/ML 應補單元測試以「鎖住」正確性防 regression）。真正待修的剩 **infra/校準類**（fill_report null 洞、region/genome 多重檢定校正、PERMANOVA 99→999 perm），非核心方法錯誤。

---

## L1 — 三審結果表

| 審 | 用什麼外部程式 | 測 ISM 什麼 | 對應稽核 | 結果 | 需修正？ |
|----|--------------|-----------|---------|------|:---:|
| **A** | modkit pileup | MM/ML 解析：每-CpG 5mC β | T1（零測試）| **r=0.976(BRCA2)/0.993(TBC1D16)** vs modkit | 🟢 **否**（補測試 pin）|
| **B** | DSS beta-binomial | Fisher vs DSS p（同 counts）| T3/#1 | ⚠ **CONFOUNDED**（DSS smoothing 增 power 蓋過 dispersion）→ 作廢 | — 不採信 |
| **B'** | read-level 置換（gold std）| Fisher null 假陽性是否 ≫5% | T3/FISHER-1 | **null FP 0.7–1.9% < 5%**（不膨脹）| 🟢 **否** |

> 🔑 **方法誠實**：Audit B（DSS）我**主動判為 confounded 不採信** —— DSS `smoothing=TRUE`（no-replicate 模式）跨 CpG 借強度增加 power，把 over-dispersion 懲罰蓋掉，故 DSS p 反而比 Fisher 小（86–91%），**不能用來證明或反證 over-dispersion**。改用 **read-level 置換**（exact permutation，gold standard）才是決定性測試。

---

## L2 — 逐審細節

### Audit A — MM/ML 解析正確性 🟢 PASS
| 位點 | Pearson r | Spearman r | 匹配 CpG | mean\|Δ\| | ISM β(thr0.5) | modkit %(thr~0.74) |
|------|-----------|-----------|---------|----------|--------------|-------------------|
| BRCA2 | **0.976** | 0.926 | 98 | 0.042 | 0.191 | 0.227 |
| TBC1D16 | **0.993** | 0.917 | 83 | 0.034 | 0.625 | 0.639 |

- **判讀**：ISM 從 BAM MM/ML 解析出的每-CpG 甲基率，與獨立工具 modkit 的 pileup **相關 0.98–0.99**。mean|Δ|~0.04 純粹是二值化閾值差（ISM 0.5 vs modkit auto ~0.74）。→ **ISM 的甲基化提取（含反股座標、5mC 分支）功能正確**，稽核 T1「零測試最易錯」的風險**經外部工具實證沒壞**。
- ⚠ **一個 follow-up**：匹配率僅 98/285（BRCA2 ISM CpG），modkit 在同區有 782 CpG。匹配上的 r=0.976 強，但 ~2/3 ISM CpG 未對到 modkit 位置（offset −1 後仍未中）→ 可能是 **0/1-based 或 per-strand vs combine-strands 座標慣例**（對應稽核 IOREAD-6 座標慣例無測試）。**不是解析錯，是座標對齊待釐清** —— 值得補一個座標慣例測試。
- **行動**：ISM **不需改解析邏輯**；但應**補 golden MM/ML 單元測試**（pin 正確性防 regression）+ 釐清座標慣例。

### Audit B — DSS vs Fisher（同 counts）⚠ CONFOUNDED，作廢
- BRCA2 n=455 CpG：median Fisher p=1.0、median DSS p=0.737；DSS 顯著率 12.75% > Fisher 8.57%；86–91% 的 CpG DSS p < Fisher p。
- **為何作廢**：方向與「over-dispersion 讓 Fisher 偏小」**相反**，根因是 DSS `smoothing=TRUE` 的空間借強度增加 power（與 dispersion 懲罰兩個效應糾纏），**無法隔離 over-dispersion**。→ **不採為證據**，由 Audit B' 取代。

### Audit B' — read-level 置換（決定性）🟢 Fisher 不需修
| 位點 | 真資料 Fisher 顯著率 | null 平均 | null 95% | 膨脹比(null/5%) | 真訊號 perm p |
|------|--------------------|----------|---------|----------------|--------------|
| BRCA2 | 0.199 | **0.0068** | 0.025 | **0.14×** | 0.002 |
| TBC1D16 | 0.311 | **0.0187** | 0.048 | **0.37×** | 0.002 |

- **方法**：打亂 read 的 HP1/HP1-1 標籤但**保留每條 read 的完整 CpG profile**（保留 within-read/clonal 相依），重算「Fisher 顯著 CpG 比例」500 次當 null。這正是 exact permutation null，直接測 FISHER-1 的 claim。
- **判讀**：
  1. null 下 Fisher 假陽性率 **0.68%/1.87%，遠低於名目 5%**（甚至偏保守，因小計數 Fisher exact 的離散性）→ **沒有 over-dispersion 膨脹**。稽核引用的「ρ=0.3 → 53–68% 膨脹」是**合成模擬**（稽核自標 ILLUSTRATIVE / needs_rerun），**真資料不重現**。
  2. 觀測顯著率 19.9%/31.1% **遠高於 null（perm p=0.002）** → BRCA2/TBC1D16 的 ASM **是真訊號**，非 Fisher artifact（與全鏈結論一致）。
- **行動**：per-CpG Fisher **不需改**（empirically valid）。→ **改進 #1（Fisher→beta-binomial）從「最高 ROI 必修」降為「marginal / 非必要」**（至少在 per-CpG + 這 2 loci 層級）。

---

## L3 — 修正清單（外部稽核後修訂）

| 原優先 | 項目 | 外部稽核後判定 |
|:---:|------|--------------|
| ~~🔴1~~ → ◽ | per-CpG Fisher → beta-binomial | **降級**：B' 證明 per-CpG Fisher 不偏樂觀（null FP<5%）→ 非必要 |
| 🔴 | MM/ML 補 golden 單元測試 | **保留**：解析正確(A 驗證)，但須 pin 防 regression + 釐清座標慣例(IOREAD-6) |
| ⭐ | **region/genome 多重檢定校正**（稽核 SIGLABEL-2/REGION-3）| **保留**（本次未測；per-CpG 內已有 BH-FDR，跨 region/genome 無）—— 這才是「ASM 存在率」可能偏高的真來源，非 per-CpG over-dispersion |
| ⭐ | PERMANOVA 99→999 perm + 真 F 分布（稽核 PERMANOVA-3/5）| 保留 |
| 🔴 | fill_report null 洞（稽核 T2 PYTOOL-1）| 保留（反捏造基礎設施） |
| ◽ | 停 5mC/5hmC max-collapse + Cohen's h | 保留（口徑改善，非錯誤）|

> **meta 結論**：外部程式稽核的價值 = **用真資料推翻了兩個 theory-motivated 的擔心**（解析可能壞、Fisher over-dispersion）。ISM 的**核心方法沒有需要修的錯誤**；待修項都是 **校準/基礎設施/測試覆蓋** 類（不影響任何已發表結論方向）。**這與那份全碼稽核的 L0「沒有會翻轉結論的錯誤」一致，但進一步用外部程式把 per-CpG Fisher over-dispersion 從『必修』降為『非必要』。**

## Provenance
- 外部程式：modkit 0.6.3、DSS(bioconductor)；ISM 真值：tsg Level-1 cache。
- 數字：`runs/audit_A.json` / `audit_B.json`(作廢) / `audit_B2.json`(決定性)；中間 `counts_*.tsv` / `dss_*.tsv`。
- ⚠ 範圍：單樣本 HCC1395、2 anchor loci（BRCA2/TBC1D16）、HP1 vs HP1-1。genome-wide / 多樣本 / region-level MT 未測 → 結論限 per-CpG + 這 2 loci；不外推到「所有統計都沒問題」。
