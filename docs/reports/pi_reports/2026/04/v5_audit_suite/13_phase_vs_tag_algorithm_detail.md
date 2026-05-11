---
title: Phase 與 Tag 演算法細節對比 — 程式碼層級拆解
date: 2026-04-28
author: liaoyoyo2001
tags: [audit, longphase-to, phasing, haplotag, algorithm, baseline, v5]
status: validated_complete
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/11_code_issues_inventory.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/12_gt_distribution_audit.md
---

# Phase 與 Tag 演算法細節對比

## 一句話摘要

> **Phase 在 VCF 上分配 GT/PS（決定 site 的 hap 歸屬）**；
> **Tag 在 BAM 上分配 HP:i:（決定 read 的 hap 歸屬）**。
> Baseline 的兩階段都用 somatic site 作 anchor，V5 在 phase 階段限制只用 PON germline、tag 階段加 somatic fallback，徹底切斷 self-phasing 連鎖。

---

## Section 1 — 兩階段功能定位

```
ClairS-TO VCF + BAM
        │
        ▼
┌──────────────────────────────────┐
│ longphase-to phase                │
│ (PhasingProcess.cpp + Graph.cpp) │
│                                   │
│ 對每個 variant 決定 GT/PS:        │
│   0|1, 1|0, .|0, 0|., 0|0...     │
└──────────────────────────────────┘
        │
        ▼  phased.vcf
        │
        ▼
┌──────────────────────────────────┐
│ longphase-to haplotag             │
│ (HaplotagProcess.cpp)            │
│                                   │
│ 對每條 read 決定 HP:i:           │
│   0, 1, 2, 11, 21, 33            │
└──────────────────────────────────┘
        │
        ▼  tagged.bam
```

| 維度 | Phase | Tag |
|------|-------|-----|
| **作用對象** | variant (site) | read |
| **資料結構** | phasing graph (variants 為節點) | read 為單位的 vote |
| **演算法核心** | edge weight + greedy graph cut | 三層 vote 多數決 |
| **輸出位置** | VCF FORMAT 欄位 (GT/PS) | BAM auxiliary (HP:i:) |
| **時間複雜度** | O(V × E) | O(R × V) |
| **是否使用 PON** | ✅ 決定哪些 variant 可當 anchor | ❌ 不直接用，間接受 phase 結果影響 |

---

## Section 2 — Phase 演算法細節

### 2.1 流程

![Figure A — Phase 演算法流程](figures/13_phase_vs_tag_algo/figA_phase_algorithm_flow.png)

| Step | 函式 | 動作 |
|:----:|------|------|
| 1 | `PhasingProcess::process()` | 讀 VCF + BAM，對每個 chunk 建 graph |
| 2 | variant origin 標記 | PON / SOMATIC / ORIGIN_UNDEFINED 三類 |
| 2.5 | **`convertNonGermlineToSomatic()` [V5 only]** | 當 `--pon-only-phasing` 為 true，**把所有非 PON 標為 SOMATIC**，使其不參與 phasing graph 主決策 |
| 3 | `addEdge()` | 對每對相鄰 variants：若 reads 同時 cover → 計算 edge weight `(refRefCount, altAltCount, refAltCount)` |
| 4 | `edgeConnectResult()` | 解 graph → 分配 hap1/hap2 → 寫 GT, PS |
| 5 | `somaticCalling()` | 標 LOH 區的 sub-genotype (GT2/GT3, PS2/PS3) |

### 2.2 V5 Phase 階段關鍵修復

**Before (Baseline)**:
```cpp
// PhasingProcess.cpp（baseline 沒有此分支）
// → 所有 ClairS-TO PASS 的 somatic 都可能成為 anchor
//   只要某 read 同時涵蓋 2 個 somatic site，就會建 edge
//   → 同 clone 的 reads 共現 ALT pattern → 強連結
//   → graph 把 3 個 somatic 全綁同一 phase block
```

**After (V5)**:
```cpp
// PhasingProcess.cpp:154-157
if(params.ponOnlyPhasing){
    vGraph->convertNonGermlineToSomatic();
}

// PhasingGraph.cpp:1139-1145
void VairiantGraph::convertNonGermlineToSomatic() {
    for(auto variantIter = variantPosType->begin();
        variantIter != variantPosType->end(); variantIter++) {
        if(variantIter->second.origin != PON) {
            variantIter->second.origin = SOMATIC;
        }
    }
}
```

→ **關鍵效果**：phasing graph 只把 PON-confirmed germline 當作 anchor。其他 variants（包含 ClairS-TO PASS 的 somatic）標為 SOMATIC origin → 進入 sub-genotype track（GT2/GT3）但不主導主 phase block 方向。

### 2.3 Phase 輸出 GT 分類（呼應 §12_gt_distribution_audit）

| GT | 含義 | Origin 來源 |
|----|------|------------|
| `0\|1`, `1\|0` | Germline het | PON |
| `.\|1`, `1\|.` | Germline hom / LOH | PON |
| `0\|0` | Somatic not in LOH | non-PON |
| `.\|0`, `0\|.` | Somatic in LOH | non-PON, has GT2/GT3 |
| `0/1`, `1/1` | Unphased | edge 不足或無 anchor |

---

## Section 3 — Tag 演算法細節

### 3.1 流程

![Figure B — Tag 演算法 getVote() 三層決策](figures/13_phase_vs_tag_algo/figB_tag_getVote_decision.png)

對每條 read：
1. 走訪 read 涵蓋的所有 variants
2. 對每個 variant：
   - `countSNPHaplotype()` (HaplotagProcess.cpp:484-495)：依 base 比對結果累計 `countMap[HAPLOTYPE1/2/HP1_1/HP2_1/HP3]`
3. `getVote()` (HaplotagProcess.cpp:512-563)：三層決策
4. confidence threshold 0.6 過濾

### 3.2 `countMap` 五個 bucket 含義

| Bucket | 來源 | 含義 |
|--------|------|------|
| `HAPLOTYPE1` | germline het ref-side / hom alt | read 應屬 HP1 |
| `HAPLOTYPE2` | germline het alt-side | read 應屬 HP2 |
| `HAPLOTYPE1_1` | somatic 在 LOH，hap1 邊 | somatic-traceable HP1 |
| `HAPLOTYPE2_1` | somatic 在 LOH，hap2 邊 | somatic-traceable HP2 |
| `HAPLOTYPE3` | somatic 模糊 (GT2=1\|1) | hp3 ambiguous |

### 3.3 Baseline (Layer 1 only) vs V5 (三層)

**Baseline `getVote()`** （概念）：
```cpp
if (germlineHP1 > 0 || germlineHP2 > 0) {
    // 多數決 → HP:i:1 或 HP:i:2
} else {
    hpResult = 0;  // ❌ 直接無 tag (没 fallback)
}
```

**V5 `getVote()`** (HaplotagProcess.cpp:512-563)：
```cpp
// Layer 1: Germline (highest priority)
if (germlineHP1 > 0 || germlineHP2 > 0) {
    germlineResult = (germlineHP1 >= germlineHP2) ? 1 : 2;
}
// Layer 1.5 [V5 NEW]: Somatic fallback
else if (somaticHP1 > 0 || somaticHP2 > 0) {
    germlineResult = (somaticHP1 >= somaticHP2) ? 1 : 2;
}
// else: undetermined (germlineResult = 0)

// Layer 2: somatic encoding
if (somaticTotal > 0) {
    if      (germlineResult == 1) hpResult = 11;  // HP:i:11
    else if (germlineResult == 2) hpResult = 21;  // HP:i:21
    else                          hpResult = 33;  // HP:i:33
} else {
    hpResult = germlineResult;  // 0/1/2
}
```

### 3.4 V5 Tag 階段三大改動

| Bug ID | 行為改動 |
|:----:|---------|
| **Bug 1-1** | 加 Layer 1.5 somatic fallback（HaplotagProcess.cpp:537-548） |
| **Bug 1-2** | HP:i:33 改用 enum 比對而非 integer literal |
| **Bug 1-3** | `countSNPHaplotype` / `countINDELHaplotype` 加 UNDEFINED guard，避免無 hap 資訊的 site 污染 vote |

---

## Section 4 — 具體例子（圖示）

![Figure C — 具體例子：3 reads × 5 variants in LOH region](figures/13_phase_vs_tag_algo/figC_concrete_example.png)

### 4.1 場景設計

LOH 區 5 個 variants：
- `g1`, `g2`：germline het（PON 已知）
- `s1`, `s2`, `s3`：somatic（ClairS-TO PASS）

3 條 reads：
- Read A：完整 cover g1 + s1 + s2 + s3 + g2（長 read）
- Read B：只 cover s1 + s2 + s3（短 read，落在 LOH 中段）
- Read C：只 cover s1 + s2（更短，無 germline）

### 4.2 Baseline 結果

**Phase 階段**：
- Read A 同時涵蓋 g1, s1, s2, s3, g2 → 建邊 `(g1, s1), (s1, s2), (s2, s3), (s3, g2)`
- s1/s2/s3 被當作 anchor 加入 phasing graph
- 3 個 somatic 因 read 共現被綁同一 phase block
- s1/s2/s3 GT = `0|0`（hap1=ref, hap2=ref）

**Tag 階段（Layer 1 only）**：
- Read A：g1 投 HP1, g2 投 HP1 → countMap[HP1]=2 → HP:i:1
- Read B：g1=g2=0，無 germline anchor → countMap[HP1]=HP[2]=0 → HP:i:0（**lost**）
- Read C：同 Read B → HP:i:0（**lost**）

→ 大量 reads 失去 tag；剩下的 Read A 全堆 HP1 → **17.3:1 bias**

### 4.3 V5 結果

**Phase 階段（`--pon-only-phasing`）**：
- `convertNonGermlineToSomatic()` → s1/s2/s3 origin = SOMATIC
- phasing graph 只用 g1, g2 當 anchor
- s1/s2/s3 進入 sub-genotype track，GT 標為 `.|0` 或 `0|.`（LOH-aware）
- Phase block 方向由 g1, g2 決定，不被 somatic 共現綁架

**Tag 階段（三層）**：
- Read A：Layer 1 (g1+g2 → HP1) + somaticTotal>0 → HP:i:11
- Read B：Layer 1=0 → Layer 1.5 (HP1_1=3, HP2_1=0) → germlineResult=1 → HP:i:11
- Read C：Layer 1=0 → Layer 1.5 (HP1_1=2, HP2_1=0) → germlineResult=1 → HP:i:11

→ 3 條 reads 全保留 + 可追溯。**全基因組平均後 HP1 ≈ HP2**，因為不同 LOH 區 / 不同 clone 有不同的 "majority direction"，自然平衡。

### 4.4 為什麼 V5 跨基因組能達到 1:1？

Baseline 的 17.3:1 是 **同 clone reads 的 self-phasing 連鎖效應**：
- 同 clone read 共享 ALT pattern → graph 強連結 → 全推向同一 phase
- 全基因組多 LOH 區可能因「同基準 reference allele 配對」傾向同一 phase
- 跨樣本固定 bias

V5 修復後：
- somatic 不主導 phase 方向 → 每個 LOH 區的 phase 由 PON germline 獨立決定
- 不同 LOH 區的 g1/g2 對應 hap1/hap2 是隨機的（由 PS block 內部選擇）
- 全基因組統計 → HP1 ≈ HP2

---

## Section 5 — 數據驗證（對應證據）

| 假說 | 證據（已驗證） | 來源 |
|------|---------------|------|
| **Phase 階段 baseline 與 V5 對 somatic GT 的決策幾乎一致** | PASS somatic GT 21,304 / 11,983 兩版本 100% 一致 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/12_gt_distribution_audit.md` §3.1 |
| **PS block orientation 在 V5 翻向（germline het 改邊）** | `1\|0` 508K → 35K, `0\|1` 559K → 1131K | `12_gt_distribution_audit.md` §3.2 |
| **Tag 階段是 self-phasing bias 的根因** | baseline.bam 17.3:1 → V5.bam ≈ 1:1，BAM 差異無法由 VCF GT 差異解釋 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/04_self_phasing_circular.md` |
| **Layer 1.5 拯救 LOH 區 reads** | AMB% 17.5% → 8.0%，clean blocks 70% → 95% | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/06_v5_sanity_bug_check.md` |
| **F1 不被影響** | ClairS-TO calling 固定 → F1 = 0.7157 vs 0.7154（噪音） | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md` |

---

## Section 6 — 關鍵程式碼引用

| 檔案 | 行號 | 函式 | 角色 |
|------|------|------|------|
| `PhasingProcess.cpp` | 154-157 | `if(params.ponOnlyPhasing)` | V5 入口：開啟 PON-only mode |
| `PhasingGraph.cpp` | 1139-1145 | `convertNonGermlineToSomatic()` | V5 核心：把非 PON 標為 SOMATIC |
| `PhasingGraph.cpp` | 615-830 | `VairiantGraph::addEdge()` | Phasing graph 建邊（共用） |
| `PhasingGraph.cpp` | 380 | `edgeConnectResult()` | Graph 解 + 寫 GT/PS |
| `HaplotagProcess.cpp` | 484-495 | `countSNPHaplotype()` | 累計 vote（V5 加 UNDEFINED guard） |
| `HaplotagProcess.cpp` | 512-563 | `getVote()` | V5 核心：三層決策 |
| `HaplotagProcess.cpp` | 565-725 | `judgeHaplotype()` | 整合 vote + confidence |

---

## Section 7 — 一頁速查

```
┌─────────────────────────────────────────────────────────────────┐
│   階段     │ Baseline                │ V5                       │
├────────────┼─────────────────────────┼──────────────────────────┤
│ Phase      │ ClairS-TO PASS somatic │ --pon-only-phasing       │
│ (VCF)      │ 可當 anchor             │ → 只用 PON germline       │
│            │ → graph 連結 somatic    │ → s1/s2/s3 標為 SOMATIC  │
│            │ → s1-s3 GT=0|0         │ → s1-s3 GT=.|0/0|.       │
│            │ (somatic NoLOH)         │ (somatic in LOH)          │
├────────────┼─────────────────────────┼──────────────────────────┤
│ Tag        │ Layer 1 only            │ Layer 1 + 1.5 + 2        │
│ (BAM)      │ → 無 germline 直接 lost │ → somatic fallback       │
│            │ → reads 大量丟          │ → HP:i:11/21/33 編碼     │
│            │ → 剩下的全堆 HP1        │ → 跨基因組自然平衡        │
├────────────┼─────────────────────────┼──────────────────────────┤
│ 結果       │ HP1:HP2 = 17.3:1 ❌    │ HP1:HP2 ≈ 1:1 ✓         │
│            │ AMB% = 17.5%            │ AMB% = 8.0%               │
│            │ Clean block ~70%        │ Clean block 95%           │
│            │ F1 = 0.7157             │ F1 = 0.7154 (Δ 噪音)     │
└─────────────────────────────────────────────────────────────────┘
```

---

## 圖檔索引

| 圖檔 | 內容 |
|------|------|
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/13_phase_vs_tag_algo/figA_phase_algorithm_flow.png` | Phase 演算法流程 baseline vs V5 |
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/13_phase_vs_tag_algo/figB_tag_getVote_decision.png` | Tag getVote() 三層決策樹 |
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/13_phase_vs_tag_algo/figC_concrete_example.png` | 3 reads × 5 variants 具體例子 |

## 跨檔交叉引用

- 程式碼 diff 詳細：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md`
- Bug 編號清單：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/11_code_issues_inventory.md`
- GT 分布證據：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/12_gt_distribution_audit.md`
- Self-phasing 因果鏈：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/04_self_phasing_circular.md`
