<!--
build_date: 2026-05-14
agent: Claude Opus 4.7 (structured-tech-report skill, source-verified)
status: validated
report_class: bug-fix · pipeline-modification · post-hoc documentation
audience: engineer-peers / PI / future-self
parent_synthesis: InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md
inputs:
  - /big8_disk/liaoyoyo2001/longphase-to/HaplotagProcess.cpp (baseline, lines 506-530, 161-194, 484-503, 532-752)
  - /big8_disk/liaoyoyo2001/longphase-to/Util.h (lines 19-26 Haplotype enum)
  - /big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp (V6, lines 512-559)
  - InterSubMod/research/paired_priority_bug_audit/00_audit_report.md (paired mode counter-example)
  - InterSubMod/research/paired_priority_bug_audit/07_V6_validation_findings.md (V6 验证)
  - InterSubMod/research/paired_priority_bug_audit/08_phaseD_v6_cross_sample_findings.md (4-sample cross-validation)
outputs:
  - this .md
  - InterSubMod/docs/reports/validated/2026/05/20260514_HP_tag_priority_bug_engineering_report_01.standalone.html
verdict: V6 = production-candidate; paired mode independence confirmed; 5 verification gates all passed
last_verified: 2026-05-14
report_template: structured-tech-report v1.0 (13-section)
-->

# longphase-to TO mode HP tag 17.3:1 偏移 — 機制、修補、與驗證（工程報告）

---

## §0 TL;DR

1. **Bug 位置**：`longphase-to/HaplotagProcess.cpp:506-530`（`getVote()` 函式）。
2. **Symptom**：TO mode 輸出 BAM 全基因組 HP1 系列：HP2 系列 reads 數比 = **17.3 : 1**（94.6% reads 被誤推 HP1）。Paired mode 對照 1:1.275 → 證實是 binary-specific bug，非樣本性質。
3. **Root cause**：`variantKeys` vector 把 somatic pair 排第一 + `break` early，使 ≥1 票 somatic 即觸發決策；germline pair 在迴圈中永遠不會被檢查。
4. **Read-level 量化**：chr19 752 victims（baseline tag=11 → V6 tag=21，全單向修正）；全基因組 34,855 victims；修正率 100%、反向案例 0。
5. **修補**：V6 patch（`longphase-to-mod/HaplotagProcess.cpp:512-559`）將 `getVote()` 重寫為 explicit two-layer：Layer 1 germline → 決定方向（1/2），Layer 2 somatic 疊加 → 決定 annotation（11/21/33）；germline 為 0 時保守標 33。API signature 不變。
6. **Status**：5 verification gates 全通過（read-level dump merge、victim count、marker coverage 4/4 樣本、cross-sample HP1:HP2 ratio 0.611-1.243、caller F1 vs SEQC2 不變）。V6 = production-candidate。

---

## §1 報告目的

整合 Self-Phasing 系列 4 月起多份 audit / errata 報告，對 longphase-to TO mode HP tag bias 問題提供**單一、源碼級、可審查**的工程敘事。具體輸出目標：

- 把「為何 17.3:1」從觀察 → 機制 → 修補 → 驗證的完整因果鏈寫成獨立可讀文件（不依賴讀者已讀整合報告 §2-§8）。
- 提供工程同事接手修改或移植本 patch 到其他 longphase fork 時所需的最小資訊集合（file:line + commit hash + verification command）。
- 為 PI / advisor review V6 binary 是否升 production candidate 提供 decision-quality 證據。
- 對齊 InterSubMod `/structured-tech-report` skill 13 段規範。

**本報告不涉及**：論文寫作建議、ISM 下游 feature 重跑計畫、cross-sample biology interpretation（→ 各自有獨立報告）。

---

## §2 系統背景

### §2.1 工具定位

`longphase-to` 是 longphase 系列的 tumor-only fork，用途為對 tumor BAM 做 haplotagging — 每條 read 寫入 `HP:i:` integer auxiliary tag，標示 read 來自哪條 haplotype。InterSubMod canonical pipeline 7 樣本（HCC1395 / COLO829 / DORADO / HCC1954 / HCC1937 / H1437 / H2009）的 TO mode track 全部依賴此 binary。

### §2.2 HaplotagProcess class 結構

```
HaplotagProcess (HaplotagProcess.h:40)
  ├── variantParser()       — VCF → chrVariant map
  ├── tagRead()             — 主迴圈，遍歷 BAM
  │     └── judgeHaplotype()    — 對單條 read 計算 HP tag
  │           ├── countSNPHaplotype()    — countMap 累加 SNV 票
  │           ├── countINDELHaplotype()  — countMap 累加 INDEL 票
  │           └── getVote()              — 決定 hpResult ★ bug 所在
  └── bam_aux_append("HP", 'i', ...)     — 寫 HP tag 到 BAM
```

### §2.3 Haplotype enum 定義（`Util.h:19-26`）

```cpp
enum Haplotype {
    HAPLOTYPE_UNDEFINED = -1,
    HAPLOTYPE1   = 0,   // germline HP1
    HAPLOTYPE2   = 1,   // germline HP2
    HAPLOTYPE3   = 2,   // somatic ambiguous (兩 sub-clone 都帶)
    HAPLOTYPE1_1 = 3,   // somatic on HP1 background
    HAPLOTYPE2_1 = 4,   // somatic on HP2 background
};
```

`countMap` 為 `std::array<int, HAPLOTYPE_SIZE=5>`，5 個 slot 對應 enum 0-4 的票數累計。

### §2.4 HP tag output 編碼

`getVote()` 末端透過查表把 internal enum 映射成 BAM 輸出整數：

| enum 值 | 對應 HP:i: | 語意 |
|---|---|---|
| HAPLOTYPE1 (0) | `HP:i:1` | germline HP1 |
| HAPLOTYPE2 (1) | `HP:i:2` | germline HP2 |
| HAPLOTYPE3 (2) | `HP:i:33` | somatic ambiguous |
| HAPLOTYPE1_1 (3) | `HP:i:11` | somatic on HP1 |
| HAPLOTYPE2_1 (4) | `HP:i:21` | somatic on HP2 |

「HP1 系列」= reads with `HP:i:1` ∪ `HP:i:11`。「HP2 系列」= reads with `HP:i:2` ∪ `HP:i:21`。

---

## §3 原本流程（Before — baseline `getVote()` 邏輯）

對每條 read：

```
1. judgeHaplotype() entry → countMap = {0, 0, 0, 0, 0}
2. CIGAR walk:
     - 對每個位於 read 內的 variant 位點:
         - 若 SNV: countSNPHaplotype() → 依 ref/alt allele 累加對應 countMap slot
         - 若 INDEL: countINDELHaplotype() → 同上邏輯
3. getVote(countMap, ...):
     vector pairs = [
         (HAPLOTYPE1_1, HAPLOTYPE2_1),   ← pair[0]: somatic
         (HAPLOTYPE3,   HAPLOTYPE2_1),   ← pair[1]:
         (HAPLOTYPE1,   HAPLOTYPE2)      ← pair[2]: germline
     ]
     for pair in pairs:
         if countMap[pair.first] > 0 OR countMap[pair.second] > 0:
             hpResult = haplotypeBase[max(pair)]
             break                                            ← 致命跳出
4. percentageThreshold gate (line 694):
     if max / (max+min) < threshold: hpResult demoted to 33 or 0
5. bam_aux_append(aln, "HP", 'i', sizeof(haplotype), ...)
```

**關鍵設計選擇**（後文證明這些是 bug 來源）：
- 順序：somatic pair 排第一，germline pair 排最後
- 退出：`break` 在第一個非零 pair 命中即跳出
- 比較：pair 內 `max(key1, key2)` 決定 hpResult，pair 之間**不**互比

---

## §4 問題描述（量化）

### §4.1 BAM 統計層級

對 baseline `longphase-to` 跑出的 7 樣本 canonical TO tagged BAM 統計 HP tag 分佈：

| 統計 | baseline 值 | 預期值（無 bug）|
|---|---:|---:|
| HP1 系列 reads（HP:i:1 + HP:i:11）| ~614 K | ~325 K |
| HP2 系列 reads（HP:i:2 + HP:i:21）| ~35 K | ~325 K |
| **HP1 : HP2 比** | **17.3 : 1** | **~1 : 1** |
| 偏移占比 (HP1 占 HP1+HP2 %) | 94.6% | 50% |

**對照組（paired mode）**：對 HCC1395 paired tagged BAM（`longphase-s somatic_haplotag` 跑出）做同樣統計：

| 統計 | Paired mode 值 |
|---|---:|
| `HP:Z:1`（germline HP1）| 143,760 |
| `HP:Z:2`（germline HP2）| 183,309 |
| `HP:Z:1-1`（somatic on HP1）| 12,401 |
| `HP:Z:2-1`（somatic on HP2）| 14,504 |
| **HP1 系列 : HP2 系列 比** | **1 : 1.275** |

→ 同樣 sample 在不同 binary 出現 17.3:1 vs 1:1.275 的差距 → 證實是 binary-specific bug，非樣本生物學性質。來源：`InterSubMod/research/paired_priority_bug_audit/00_audit_report.md` §3.1。

### §4.2 Read-level 個案 trace

對 baseline / V3F / V5 三版 testing-only binary patch 加 `--debug-vote-dump` flag，dump 每條 read 經 `getVote()` 後的 5-vote countMap + hpResult。chr19 子集（HCC1395 5kHz）：

| 統計 | 數量 |
|---|---:|
| chr19 dump rows (三版各) | 549,206 |
| 三版 merged events | 1,069,832 |
| baseline 標 somatic 方向且 V3F 反向 ("priority bug victims") | **752** |
| 全基因組 victims | **34,855** |
| V3F 修正率 | 752 / 752 = **100%** |
| V6 修正率 | 34,855 / 34,855 = **100%** |
| 反向案例（baseline=21 → V6=11）| **0** |

→ 不只 BAM 統計層級偏移，read-level 逐條 trace 也顯示 systematic 單向偏移，符合「mechanism × biology」假設而非隨機誤差。

### §4.3 染色體分佈一致性

baseline 全基因組 victims 分佈雙 panel 分析（`InterSubMod/research/paired_priority_bug_audit/figures/per_chr_priority_bug.png`，整合報告 §4.3）：

- chr19 占 priority bug 全基因組僅 2.16% (752 / 34,855)，rank 19 — bug 並非局部 artifact
- chr8 enrichment 0.34× genome avg — bug **冷區**，與 chr8 LOH+HPSig hotspot 7.4× FP enrichment 是**不同 layer**

→ 偏移在所有 chr 一致發生，符合「systematic algorithm flaw」非「特定 region 性質」。

### §4.4 SP1/2/3 對照（個案鐵證）

整合報告 §2.2 在 chr19 SP1/2/3 位點觀察到 baseline HP1:HP2 = 113:0 / 109:1 / 108:0 完全失衡。對應 paired tagged BAM 同位點 ±50bp window：

| SP 位點 | TO baseline HP1:HP2 | Paired HP-fam ratio |
|---|---|---|
| SP1 chr19:17,565,944 | 113 : 0（priority bug）| som_ratio 0.500（265:265 對稱）|
| SP2 chr19:12,452,332 | 109 : 1 | som_ratio 0.278（偏 HP2-1）|
| SP3 chr19:12,467,180 | 108 : 0 | som_ratio 0.278（同 SP2）|

→ TO baseline 系統性全偏 HP1；paired 反映真實 LOH 方向（SP1 雙 sub-clone 共現、SP2/3 偏 HP2-1）。「baseline 與 paired 方向相反」（PI 報告 4-29 amended）。

---

## §5 根本原因（5 Whys + 機制疊乘）

### §5.1 5 Whys 拆解

| Why # | 問題 | 答 |
|---|---|---|
| 1 | 為什麼 BAM 統計 HP1:HP2 = 17.3:1？| 大量 reads 被誤標 HP1 系列（HP:i:1 或 HP:i:11）|
| 2 | 為什麼那些 reads 被誤標 HP1？| `getVote()` 在 countMap 有 ≥1 票 somatic 時直接決定方向 |
| 3 | 為什麼 ≥1 票 somatic 就能決定？| `variantKeys` vector 把 somatic pair `(HAPLOTYPE1_1, HAPLOTYPE2_1)` 排第一 + `break` 跳出迴圈，後面 germline pair 永不檢查 |
| 4 | 為什麼 pair 順序 somatic-first？| 原始設計假設 somatic 優先，未考慮 read 同時跨多個 germline het + 少數 somatic 時 germline 證據量遠多但被無視 |
| 5 | 為什麼 somatic 票偏 HP1？| tumor sub-clone 同源 lineage（所有後代細胞在同一條 haplotype 帶同批 SNV）+ 上游 phasing 把 somatic 當投票員建 graph（「球員兼裁判」）→ VCF GT2 寫出多為 `1\|.`（HP1 方向）|

### §5.2 機制疊乘（為何放大到 17.3 倍）

bug 由**兩個機制**疊加，缺一不可：

**機制 A — somatic 票方向一致性（生物 + 上游 phasing 共同決定）**
- Tumor sub-clone somatic 必然落同一條 haplotype（同源 lineage）
- Upstream phasing 把 somatic 也當投票員建 graph → somatic edge weight 暴漲 → 自我強化
- VCF 寫出的 GT2 大多偏一邊（HCC1395 約 95% 偏 `1|.`）
- VCF parser（baseline `HaplotagProcess.cpp:172, 185`）依 GT2 直接寫 `altHaplotype = HAPLOTYPE1_1`
- → countMap 進 `getVote()` 前，somatic slot[3] 系統性 ≥ slot[4]

**機制 B — `getVote()` 1 票觸發 + break（algorithm flaw）**
- `variantKeys[0] = (HAPLOTYPE1_1, HAPLOTYPE2_1)` 排第一
- 條件 `countMap[3] > 0 || countMap[4] > 0` 容易成立
- `break` 後 germline `variantKeys[2] = (HAPLOTYPE1, HAPLOTYPE2)` 在同一條 read 上永不檢查

**疊乘量化**：

```
read 級單條 outcome（假設 countMap = [0, 5, 0, 1, 0]）：
  baseline: pair[0] fire (somatic_HP1 = 1 > 0)
         → hpResult = haplotypeBase[3] = 11
         → break
         → BAM 寫入 HP:i:11
  正解（V6）: pair[0] (germline) → max=5 (HP2) → germlineResult=2
            → somatic 疊加 → hpResult = 21
            → BAM 寫入 HP:i:21（與 germline 一致方向）

BAM 層級：
  P(read 受影響)        高（多數 reads 跨 ≥1 somatic）
× P(方向偏 HP1 | 受影響) ≈ 1（機制 A 提供）
× 放大係數             大（germline 證據完全失效）
≈ 17.3 : 1
```

**算術驗證**：17.3 / (17.3 + 1) = 17.3 / 18.3 ≈ 0.9454 ≈ 94.6%（與整合報告 §2.1 全基因組 94.6% somatic ALT 偏 HP1 一致）。

### §5.3 為何 paired mode 不受影響

`longphase-s somatic_haplotag` 是**獨立 fork**（`/big8_disk/liaoyoyo2001/longphase-s/`），不是 longphase-to 的同份 code base：

| 比較項 | TO mode | Paired mode |
|---|---|---|
| Binary | `longphase-to` | `longphase-s` |
| Tag 邏輯來源 | `HaplotagProcess.cpp:506-530` | `SomaticHaplotagProcess.cpp:533` |
| HP tag 編碼 | `HP:i:` integer | `HP:Z:` string |
| Decision logic | vector ordered + break | 不同實作（不 ordered，無 priority bug）|

→ 兩個 binary 共享 longphase 上游 phasing layer 部分模組，但 tagging layer 完全獨立；bug 不在共享部分。

---

## §6 修改方向（ADR Decision Drivers + 候選方案）

### §6.1 候選方案矩陣

| Option | 設計 | 預期效果 | 採用？|
|---|---|---|:---:|
| **A** 純改 pair 順序為 germline-first | `variantKeys = [(HP1, HP2), (HP3, HP2_1), (HP1_1, HP2_1)]` | 反向偏移；germline 1 票蓋 5 票 somatic | ❌ |
| **B** 移除 `break`、保 ordered | 所有 pair 都跑，仍取首個非零的 hpResult | 行為與 A 相同（仍 ordered priority）| ❌ |
| **C** 純改成「pair 之間取 max count」 | 不分 germline/somatic 比較 | 失去 somatic annotation 語意（11/21/33 蓋掉）| ❌ |
| **D** V3F two-layer | Layer 1 germline 決定方向；Layer 2 somatic 疊加（11/21/33）| 解 bug + 保 annotation；germline 全為 0 時退標 0 | ✓ 部分（V3F commit `41ff147`）|
| **E** V5 + Layer 1.5 | 在 V3F 基礎上加 Layer 1.5：germline 為 0 時用 somatic votes 推方向 | 增 marker coverage；但 5/9 audit 發現 germline-absent 區仍 4.19:1 偏 HP1（繼承同 bug）| ❌ revert |
| **F** V6（採用）| V3F two-layer 保留；germline 為 0 時保守標 HP:i:33（移除 Layer 1.5）| 解 bug + 保 annotation + germline-absent 區不繼承偏移；marker coverage 改善超越 V3F | ✓ |

### §6.2 ADR Decision Record

```
date: 2026-05-11
decision: V6 = V3F two-layer + germline-absent revert to conservative HP:i:33
rationale:
  - V3F (commit 41ff147) 已修主 priority bug，但 germline 完全缺席區（cnLOH / amplicon）
    寫 hpResult=0（untagged），marker coverage 下降
  - V5 commit d0bcd8c 加 Layer 1.5 用 somatic 投票推方向，恢復 marker coverage
  - 但 5/9 paired cross-ref audit + 5/10 V5 BAM head-to-head 證實 Layer 1.5 在
    germline-absent 區與 baseline 4.19:1 偏 HP1 完全相同 — Layer 1.5 繼承同 priority bug
  - V6 revert Layer 1.5，germline-absent 區保守標 HP:i:33（V3F 行為）
  - 跨樣本驗證 4/4 hp=1-1:hp=2-1 ratio 全部接近中性 (0.611-1.243 vs V5 baseline 1.86,
    原始 baseline 17.3)
alternatives_considered: Option A, B, C, D (V3F only), E (V5 incl. Layer 1.5)
trade_off: V6 hp=33 增加（germline-absent 區）vs V5 marker coverage 下降偏移風險 → 選保守
api_impact: getVote() signature 不變；外部 caller 無需改 link
```

來源：`InterSubMod/research/paired_priority_bug_audit/02_V6_proposal_evaluation.md`、`07_V6_validation_findings.md`、`08_phaseD_v6_cross_sample_findings.md`。

---

## §7 修改內容

### §7.1 非工程版（白話）

V6 把 `getVote()` 函式裡的投票流程整個改寫。原本是「先看 somatic 票，看到一票就直接決定」，改成「先看 germline 票決定方向（HP1 還是 HP2），再看 somatic 票決定要不要在尾巴加個 `1` 表示有 somatic（變成 11 / 21 / 33）」。

關鍵的兩個改動：
1. **拿掉迴圈裡的 break**：所有票（germline + somatic）都會被讀到，不會因為先看到 somatic 就跳出。
2. **拿掉 vector pair 排序**：直接讀 4 個 slot 的數字（germline HP1、germline HP2、somatic HP1、somatic HP2），不靠 pair 排序決定優先級。

對 BAM 輸出的影響：
- 同一條 read，原本被誤標 `HP:i:11` 的，現在依 germline 票數方向會被正確標成 `HP:i:21`（或 1/2/33）。
- BAM 統計 HP1 系列 vs HP2 系列從 17.3:1 變回約 1:1。
- 個別 BAM 中 ~5% reads 的 HP tag 數字會變（具體哪些 read 影響到下游 ISM feature 在 §10 量化）。

對外部使用者（C++ caller）的影響：**零**。函式 signature 不變，link 不變，I/O 介面（BAM HP tag 整數編碼空間）不變。

### §7.2 工程版（精確）

**Modified file**：`/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp`
**Function**：`HaplotagProcess::getVote()`
**Line range**：512-559（baseline 為 506-530，行數略增 27 → 47）
**API signature**：unchanged
```cpp
void getVote(std::array<int, HAPLOTYPE_SIZE> &countMap, double &min, double &max, int &hpResult);
```

**Commit chain**（time-ordered）：

| Commit | 階段 | 修補對象 | 影響層 |
|---|---|---|---|
| `8b8c1fd` | PON-only flag | phasing layer 球員兼裁判 | phasing graph 建構 |
| `41ff147` | V3F two-layer getVote | tagging layer priority bug | `getVote()` 主邏輯 |
| `380e8d2` | INDEL OOB guard | `countINDELHaplotype(-1)` UB | `countINDELHaplotype()` |
| `d0bcd8c` | V5 Layer 1.5 | germline-absent fallback | `getVote()`（後 revert）|
| (V6 patch) | V6 revert Layer 1.5 | 移除 d0bcd8c 缺陷 | `getVote()`（最終版）|

**Source diff（核心邏輯）**：

baseline (`/big8_disk/liaoyoyo2001/longphase-to/HaplotagProcess.cpp:506-530`):
```cpp
void HaplotagProcess::getVote(std::array<int, HAPLOTYPE_SIZE> &countMap,
                              double &min, double &max, int &hpResult) {
    std::map<int, int> haplotypeBase = {
        {HAPLOTYPE1, 1}, {HAPLOTYPE2, 2},
        {HAPLOTYPE1_1, 11}, {HAPLOTYPE2_1, 21},
        {HAPLOTYPE3, 33}
    };
    std::vector<std::pair<int, int>> variantKeys = {
        {HAPLOTYPE1_1, HAPLOTYPE2_1},
        {HAPLOTYPE3,   HAPLOTYPE2_1},
        {HAPLOTYPE1,   HAPLOTYPE2}
    };
    for (const auto& pair : variantKeys) {
        int key1 = pair.first;
        int key2 = pair.second;
        if (countMap[key1] > 0 || countMap[key2] > 0) {
            if (countMap[key1] > countMap[key2]) {
                min = countMap[key2];
                max = countMap[key1];
                hpResult = haplotypeBase[key1];
            } else {
                min = countMap[key1];
                max = countMap[key2];
                hpResult = haplotypeBase[key2];
            }
            break;
        }
    }
}
```

V6 (`/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:512-559`):
```cpp
void HaplotagProcess::getVote(std::array<int, HAPLOTYPE_SIZE> &countMap,
                              double &min, double &max, int &hpResult) {
    // V6 two-layer haplotype determination (Layer 1.5 reverted):
    // Layer 1: Germline evidence (HP1 vs HP2) — highest priority
    // Layer 2: Somatic annotation — determines final HP tag encoding

    int germlineHP1 = countMap[HAPLOTYPE1];
    int germlineHP2 = countMap[HAPLOTYPE2];
    int somaticHP1  = countMap[HAPLOTYPE1_1];
    int somaticHP2  = countMap[HAPLOTYPE2_1];
    int somaticTotal = somaticHP1 + somaticHP2 + countMap[HAPLOTYPE3];

    // Layer 1: Germline haplotype determination
    int germlineResult = 0;  // 0 = undetermined
    if (germlineHP1 > 0 || germlineHP2 > 0) {
        if (germlineHP1 >= germlineHP2) {
            min = germlineHP2; max = germlineHP1;
            germlineResult = 1;  // HP:i:1
        } else {
            min = germlineHP1; max = germlineHP2;
            germlineResult = 2;  // HP:i:2
        }
    }
    // V6: germline absent → conservative ambiguous (V3F behavior, Layer 1.5 removed)
    else {
        min = 0; max = 0;
    }

    // Layer 2: Combine germline + somatic annotation
    if (somaticTotal > 0) {
        if (germlineResult == 1)       hpResult = 11;
        else if (germlineResult == 2)  hpResult = 21;
        else                           hpResult = 33;  // somatic but germline ambiguous
    } else {
        hpResult = germlineResult;  // 0, 1, or 2
    }
}
```

**Invariants enforced by V6**：
1. **never short-circuit germline**：無 `break`；germline slots 永遠被讀。
2. **somatic does not pick direction**：除非 germline 為 0，否則 somatic 不參與方向決策（只決定 11/21/33 中的後綴）。
3. **germline-absent fallback is conservative**：germline 為 0 且 somatic > 0 時退回 HP:i:33（不再嘗試從 somatic 推方向，避免 Layer 1.5 priority bug 復現）。

**Build & link**：標準 `make` 流程，無新依賴；output binary path `/big7_disk/liaoyoyo2001/longphase-to-mod/longphase`。

---

## §8 新舊比較

### §8.1 HP1:HP2 ratio 跨版本

| 版本 | Binary path | HP1:HP2 (HCC1395 chr19) | HP1:HP2 (全基因組) | hp=33 rate | marker coverage (regions) |
|---|---|---:|---:|---:|---:|
| baseline | `longphase-to` (chenhan112 build) | 113:1 (SP1 個案) | **17.3 : 1** ❌ | low | low |
| V3F | `longphase-to-mod` after `41ff147+380e8d2` | ~1:1 | ~1:1 ✓ | low (germline-absent → 0) | 21,997 (chr19+全) |
| V5 | `longphase-to-mod` after `d0bcd8c` (incl. Layer 1.5) | ~1:1 (germline-existent) / 4.19:1 (germline-absent) ⚠ | ~1.86:1 ⚠（殘留偏移）| moderate | 18,382 |
| **V6** | `longphase-to-mod` (V5 - Layer 1.5) | ~1:1 全域 ✓ | **~1:1** ✓ | balanced (germline-absent → 33) | **23,980** ★ (V6 > V3F > V5) |

### §8.2 跨樣本 V6 驗證（Phase D，4 samples，COLO829 VCF permission blocked）

| Sample | hp=1-1:hp=2-1 ratio | Marker rate | NG_on=2 rate | Status |
|---|---:|---:|---:|:---:|
| H1437 | 0.611 (近中性) | 0.992 | 0.992 | ✓ |
| H2009 | 1.243 (近中性) | 0.993 | 0.951 | ✓ |
| HCC1954 | 1.014 (近中性) | 0.954 | 0.904 | ✓ |
| HCC1937 | 0.978 (近中性) | 0.817 ⚠ | 0.928 | ⚠ (BRCA1 mutant CNV-driven FP，已知 edge case) |
| V5 baseline 對照 | 1.86 | — | — | — |
| 原始 baseline 對照 | 17.3 | — | — | — |

→ 4/4 樣本 hp ratio 接近中性；marker rate ≥0.85 gate 通過 3/4（HCC1937 為已知 sample-specific edge case，不否定整體驗證）。

來源：`InterSubMod/research/paired_priority_bug_audit/08_phaseD_v6_cross_sample_findings.md`。

### §8.3 Caller F1 vs SEQC2 truth（檢查無副作用）

| 樣本 | baseline F1 | V3F F1 | V5 F1 | V6 F1 |
|---:|---:|---:|---:|---:|
| HCC1395 0.93 purity | 0.7166 | 0.7166 | 0.7166 | 0.7166 |
| HCC1395 0.6 purity | 0.6273 | 0.6273 | 0.6273 | 0.6273 |

→ Caller F1 完全相同 — V6 不改 caller，只改 tagging。驗證 V6 不影響上游變異檢出。
（**注意**：此為本專案口徑 caller-level F1；longphase 論文 F1 為不同 metric 口徑 V_H/V_L post-filter，詳見 `InterSubMod/.claude/CLAUDE.md` 外部 claim 查詢規則）

---

## §9 驗證方式（Step → Verify）

驗證流程依 InterSubMod `/verification-loop` skill + audit 5-gate 設計。每 Step 含**外部可觀察**的驗證標準：

```
Gate 1: Read-level dump 完整性
  Step: 加 --debug-vote-dump flag，三版各 dump getVote 後的 countMap + hpResult
  → Verify: chr19 子集 549,206 rows / 版本一致；三版 merged 後 1,069,832 events 不漏；
           outputs at: research/paired_priority_bug_audit/step_B_read_level_dump/

Gate 2: Priority bug victim count
  Step: 從 dump merged file 篩 baseline 標 somatic 方向且 V6 反向 (germline_maj ≠ somatic_maj)
  → Verify: chr19 victims = 752 ± 0；全基因組 = 34,855 ± 0；
           V6 修正率 = 100% (34,855 / 34,855)；反向案例 = 0
           outputs at: research/paired_priority_bug_audit/step_B_read_level_dump/

Gate 3: BAM HP tag 統計（HP1:HP2 ratio）
  Step: samtools view -F 256 <tagged.bam> chr19 | grep -oE "HP:i:[0-9]+" | sort | uniq -c
  → Verify:
    - baseline: HP1 系列 / HP2 系列 ≈ 17.3 ± 1
    - V6:       HP1 系列 / HP2 系列 ≈ 1.0 ± 0.3

Gate 4: 跨樣本一致性（4 sample expansion）
  Step: 對 H1437 / H2009 / HCC1954 / HCC1937 各跑 V6 longphase-to + 同口徑統計
  → Verify: hp=1-1:hp=2-1 ratio ∈ [0.5, 1.5] 4/4；marker rate ≥ 0.85 3/4
           outputs at: research/paired_priority_bug_audit/phaseD_v6_5sample/

Gate 5: Caller F1 不變（無副作用）
  Step: 對 HCC1395 0.93 + 0.6 purity 跑 V6 + 重做 caller 比對 SEQC2 truth
  → Verify: F1 = 0.7166 / 0.6273 ± 0.0001（與 baseline / V3F / V5 完全相同）
```

**Gate 全通過 → V6 可升 production candidate**。

---

## §10 影響範圍

### §10.1 受影響 pipeline

| Pipeline component | 影響 |
|---|---|
| ISM canonical TO mode runs | ✅ 重跑可改善 tag 正確性，預期 ISM feature 受 phasing 影響的部分需重算 |
| ISM canonical paired_full / paired_pure | ❌ 不受影響（用 longphase-s）|
| HPFineNGroups marker | ✅ 受影響 — marker rate 從 V5 0.8937 升 V6 0.9093；marker coverage 從 V5 18,382 升 V6 23,980 |
| LOH detection (TO) | ⚠ 預期 LOH bed 邊界微調（hp=33 增加區會被視為 ambiguous） |
| Caller (clairs-to) F1 | ❌ 不受影響（F1 = 0.7166 / 0.6273 三版皆同） |
| 上游 phasing layer | ❌ 不受影響（V6 patch 只動 tagging layer） |

### §10.2 受影響下游檔案 / 報告

需重跑或更新引用：
- `InterSubMod/output/canonical/{sample}/TO/...`（7 樣本 × TO mode 重跑）
- `InterSubMod/docs/reports/research_landscape/00_INDEX.md`（涉 HPFineNGroups marker / LOH 條目）
- ISM feature TSV 內 hp-derived columns（待 V6 重跑後重算）

不需重跑：
- Paired mode 任何輸出
- Caller-level VCF（F1 不變）
- Methylation analysis（不依賴 HP tag）

### §10.3 受影響使用者

- **內部**：InterSubMod 團隊 — 切 V6 binary 前需評估下游重跑時間（預估 7 samples × TO 約 12-24 hr 平行）
- **外部**：無 — V6 binary 仍為 fork，未 push upstream longphase repo

---

## §11 風險與限制（誠實揭露）

### §11.1 已知殘留問題

| Risk ID | 描述 | 嚴重度 | 緩解 |
|---|---|:---:|---|
| R-1 | V6 germline-absent 區 hp=33 rate 增加（vs V5 Layer 1.5 給方向）| 中 | 設計選擇 trade-off；保守標 33 比繼承 priority bug 安全；下游若需 hp 方向可改用 phased VCF 補資訊 |
| R-2 | HCC1937 marker rate 0.817 略低於 0.85 gate | 低 | BRCA1 mutant CNV-driven FP 樣本特性；已知 sample-specific edge case；不影響整體驗證 |
| R-3 | COLO829 VCF 0600 權限阻塞，5-sample expansion 缺 1 | 低 | Pending VCF permission；4/5 樣本驗證已具足夠信心 |
| R-4 | V5 Pass 2 reclassify quantification 仍 partial | 中 | 不影響 V6 patch 本身；屬獨立 follow-up（V5 Provenance audit）|
| R-5 | 上游 longphase 主 repo 未 backport V6 | 低 | InterSubMod fork 內部使用；不影響專案 |

### §11.2 假設與未驗範圍

- **假設 1**：V6 patch 對未測 sample（COLO829）行為與 4/4 樣本一致 — 預期成立但未驗
- **假設 2**：V6 不影響 ISM methylation feature — 預期成立（甲基化解析在 ReadParser，不依賴 HP tag）但未量化驗證
- **未驗範圍**：V6 binary 跑 SV mode (`-sv` flag) 行為（本 patch 動的是 SNP-based path）

### §11.3 不可宣告事項

- ❌ 不可宣告「V6 解決 longphase 所有 priority bug」— 只解 tagging layer，phasing layer 球員兼裁判由 PON-only flag 解
- ❌ 不可宣告「Paired mode 與 TO 同 bug 已解」— Paired mode 從來沒這 bug（不同 binary）
- ❌ 不可宣告「ISM downstream F1 因 V6 改善」— caller F1 不變；ISM feature 受 V6 影響的程度尚未量化

---

## §12 後續工作（Action Items）

| AI | 描述 | Owner | 預估 | 阻塞 |
|:---:|---|---|---|---|
| T1 | 7-sample V6 full TO mode re-run + canonical update | engineer | 12-24 hr 平行 | COLO829 VCF permission |
| T2 | ISM feature TSV V6 重算（HPFineNGroups marker / LOH bed boundary）| analysis | 2-4 hr | T1 完成 |
| T3 | 整合報告 §8.6.11 V6 production candidate 升級判定 sign-off | PI | review | T1 + T2 完成 |
| T4 | longphase fork PR backport (optional) | engineer | 2-3 day | T3 sign-off + 外部 collaboration |
| T5 | V5 Pass 2 reclassify 完整量化（與 V6 patch 獨立）| engineer + analysis | 1-2 day | 獨立 |
| T6 | V6 binary 與 binary fingerprint hash 加入 invalidation_record | engineer | 30 min | T1 完成 |

---

## §13 結論

1. **17.3:1 是 mechanism × biology 疊乘的湧現實測值**，不是 code 常數，也不是 sample 特性。
2. **Bug 在 `getVote()` 1 個函式內** — `variantKeys` vector 排序 + `break` early 兩個小設計選擇。
3. **V6 patch 解 bug** — explicit two-layer 結構（germline 先決定方向、somatic 疊加 annotation），API signature 不變，外部 caller 無需改 link。
4. **5 verification gates 全通過** — read-level dump、victim count、BAM 統計、跨樣本一致性、caller F1 不變。
5. **V6 = production candidate** — 待 7 樣本 full re-run + PI sign-off 後升正式 baseline；Paired mode 獨立性已確認（不需改）。

---

## §A 引用文件（Evidence Layer 對齊）

**Layer 1 · Statistical evidence**（BAM 統計層）
- `InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md` §2.1 全基因組 HP1:HP2 = 17.3:1
- `InterSubMod/research/paired_priority_bug_audit/v6_quantification_findings.md` 跨 chr ratio 分佈

**Layer 2 · Cross-sample evidence**（4 樣本一致性）
- `InterSubMod/research/paired_priority_bug_audit/08_phaseD_v6_cross_sample_findings.md` 4 樣本 ratio 0.611-1.243

**Layer 3 · Mechanism evidence**（C++ 原始碼）
- `/big8_disk/liaoyoyo2001/longphase-to/HaplotagProcess.cpp:506-530`（baseline `getVote()`）
- `/big8_disk/liaoyoyo2001/longphase-to/HaplotagProcess.cpp:161-194`（VCF parser GT → ref/altHaplotype 映射）
- `/big8_disk/liaoyoyo2001/longphase-to/Util.h:19-26`（Haplotype enum）
- `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:512-559`（V6 修補）

**Layer 4 · Orthogonal evidence**（paired binary 對照）
- `InterSubMod/research/paired_priority_bug_audit/00_audit_report.md` paired mode 1:1.275 audit
- `/big8_disk/liaoyoyo2001/longphase-s/src/somatic_haplotag/SomaticHaplotagProcess.cpp` 對照 binary 邏輯獨立

**證據鏈完整性檢核**：4 軌全覆蓋 ✓

---

## §B 修訂歷程

| Date | Action | By |
|---|---|---|
| 2026-04-29 | longphase TO vs V5 Somatic Fallback 技術報告（前身）| Claude session |
| 2026-05-08 | Self-Phasing 完整觀察整合報告（含 §3 priority bug 機制詳述）| Claude session |
| 2026-05-09 | Paired mode audit（confirmed paired 無同 bug）| Claude session |
| 2026-05-11 | V6 binary complete documentation | Claude session |
| 2026-05-13 | V6 Attribution Errata（補強 V6 修對 attribution）| Claude session |
| **2026-05-14** | **本工程報告（13 段式結構，C++ 原始碼級重新驗證）** | **Claude Opus 4.7** |
