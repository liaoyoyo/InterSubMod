---
title: Baseline GT2/GT3 sub-genotype 機制再確認 — 既有報告的精確化更正
date: 2026-04-28
author: liaoyoyo2001
tags: [erratum, baseline, V5, sub-genotype, GT2, GT3, getVote]
status: validated_complete
audience: developer + PI
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/11_code_issues_inventory.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/14_user_report_integrated.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/15_software_engineering_perspective.md
---

# Baseline GT2/GT3 sub-genotype 機制再確認

## §0 一句話結論

> **Baseline 有完整的 GT2/GT3 sub-genotype 機制 + 在 `getVote()` 內以 SOMATIC-first 順序投票**，並非「沒有 fallback」。V5 的真實修法是**順序反轉**（somatic-first → germline-first）+ **encoding 拆分**（直接判 11/21/33 → 先決定 germlineResult 再加 somatic encoding）。17.3:1 bias 與 0.6 simulation 結論**仍然成立**，但根因描述需要精確化。

---

## §1 問題起源與用戶提問

> **用戶提問**：「確認當時 baseline 是否將判斷 somatic 放到其他 GT2 GT3 去判讀，例如要判 52 萬 read 是否是 HP1-1 需要多個標籤確認，這會影響到 baseline phase 與 tag 的結論和之前觀察的問題結論嗎」

→ 這個問題揭露了之前報告的描述偏差。實際查閱程式碼發現：**baseline 確實有 sub-genotype 機制**，本文做完整更正。

---

## §2 Baseline 真實程式碼（git show 8b8c1fd^）

### 2.1 Baseline `getVote()` 實際內容

```cpp
// Baseline (8b8c1fd^ 之前版本) HaplotagProcess.cpp::getVote()
void HaplotagProcess::getVote(std::array<int, HAPLOTYPE_SIZE> &countMap,
                              double &min, double &max, int &hpResult) {
    std::map<int, int> haplotypeBase = {
        {HAPLOTYPE1, 1}, {HAPLOTYPE2, 2},
        {HAPLOTYPE1_1, 11}, {HAPLOTYPE2_1, 21},
        {HAPLOTYPE3, 33}
    };
    std::vector<std::pair<int, int>> variantKeys = {
        {HAPLOTYPE1_1, HAPLOTYPE2_1},   // ← 第 1 優先（somatic-traceable）
        {HAPLOTYPE3, HAPLOTYPE2_1},      // ← 第 2 優先（HP3 vs HP2_1）
        {HAPLOTYPE1, HAPLOTYPE2}         // ← 第 3 優先（germline）
    };
    for (const auto& pair : variantKeys) {
        int key1 = pair.first;
        int key2 = pair.second;
        if (countMap[key1] > 0 || countMap[key2] > 0) {
            if (countMap[key1] > countMap[key2]) {
                hpResult = haplotypeBase[key1];  // 直接判 11/21/33
            } else {
                hpResult = haplotypeBase[key2];
            }
            break;  // ← 找到第一組就跳出
        }
    }
}
```

### 2.2 關鍵觀察

1. **Baseline 已有 5 個 HAPLOTYPE buckets**：
   ```cpp
   // Util.h:24 (baseline 與 V5 共用)
   enum {
       HAPLOTYPE1 = 0,    // germline HP1
       HAPLOTYPE2 = 1,    // germline HP2
       HAPLOTYPE3 = 2,    // ambiguous somatic
       HAPLOTYPE1_1 = 3,  // somatic-traceable HP1
       HAPLOTYPE2_1 = 4,  // somatic-traceable HP2
   };
   ```

2. **Baseline 已會解析 GT2/GT3 寫入 sub-genotype haplotype**：
   ```cpp
   // HaplotagProcess.cpp:160-189 (parser logic, baseline 與 V5 共用)
   else if (gtValue[0] == '0' && gtValue[2] == '0') {  // GT = 0|0 (somatic NoLOH)
       std::string gt2Value = formatSample.getValue("GT2");
       if (gt2Value[0] == '1' && gt2Value[2] == '.') {
           chrVariant[chr][pos].altHaplotype = HAPLOTYPE1_1;  // ← sub-genotype HP1
       }
       else if (gt2Value[0] == '.' && gt2Value[2] == '1') {
           chrVariant[chr][pos].altHaplotype = HAPLOTYPE2_1;  // ← sub-genotype HP2
       }
       else if (gt2Value[0] == '1' && gt2Value[2] == '1') {
           chrVariant[chr][pos].altHaplotype = HAPLOTYPE3;    // ← ambiguous
       }
   }
   ```

3. **`getVote()` 在 baseline 已用 sub-genotype 投票**，且**順序為 somatic-first**：
   - 第 1 優先：`(HAPLOTYPE1_1, HAPLOTYPE2_1)` → somatic-traceable
   - 第 2 優先：`(HAPLOTYPE3, HAPLOTYPE2_1)` → HP3 模糊
   - 第 3 優先：`(HAPLOTYPE1, HAPLOTYPE2)` → germline

### 2.3 V5 commit 41ff147 的修法明確揭露根因

V5 commit message 直接寫：

```
fix(haplotag): two-layer getVote — germline first, somatic second

Root cause: getVote() checked somatic pairs (HP1_1/HP2_1) BEFORE germline
pair (HP1/HP2). A single somatic vote overrode all germline votes, causing
99.9% of reads to get HP:i:21 in PON-only mode.
```

→ 完全確認：baseline `getVote()` 是 **somatic-first**，V5 修為 **germline-first**。

---

## §3 V5 修法的精確版本

### 3.1 V5 `getVote()`（HaplotagProcess.cpp:512-563）

```cpp
void HaplotagProcess::getVote(...) {
    int germlineHP1 = countMap[HAPLOTYPE1];
    int germlineHP2 = countMap[HAPLOTYPE2];
    int somaticHP1 = countMap[HAPLOTYPE1_1];
    int somaticHP2 = countMap[HAPLOTYPE2_1];
    int somaticTotal = somaticHP1 + somaticHP2 + countMap[HAPLOTYPE3];

    // Layer 1: Germline first（順序反轉）
    int germlineResult = 0;
    if (germlineHP1 > 0 || germlineHP2 > 0) {
        germlineResult = (germlineHP1 >= germlineHP2) ? 1 : 2;
    }
    // Layer 1.5: Somatic fallback（baseline 原是第 1 優先，V5 改為 fallback）
    else if (somaticHP1 > 0 || somaticHP2 > 0) {
        germlineResult = (somaticHP1 >= somaticHP2) ? 1 : 2;
    }

    // Layer 2: Encoding（baseline 直接判 11/21/33，V5 從 germlineResult 推）
    if (somaticTotal > 0) {
        if      (germlineResult == 1) hpResult = 11;
        else if (germlineResult == 2) hpResult = 21;
        else                          hpResult = 33;
    } else {
        hpResult = germlineResult;
    }
}
```

### 3.2 Baseline vs V5 真實差異

| 項目 | Baseline | V5 |
|------|----------|----|
| 是否有 sub-genotype 機制 | ✅ 有（GT2/GT3 + HAPLOTYPE1_1/HAPLOTYPE2_1/HAPLOTYPE3）| ✅ 有（與 baseline 同） |
| `getVote()` 投票順序 | **SOMATIC first** | **GERMLINE first** |
| Encoding 邏輯 | 直接判 11/21/33（vote winner 直接給 hpResult） | 先決 germlineResult，再加 somatic encoding |
| 投票結果可逆推 germline | ❌ 不行（直接給 11/21） | ✅ 可（germlineResult 是 1/2） |

---

## §4 17.3:1 bias 根因再釐清

### 4.1 兩個獨立 bug 同時存在

baseline 17.3:1 bias 來自**兩個層面的 bug 疊加**，不是單一機制：

| Bug | 階段 | 影響 | V5 修法 |
|-----|------|------|--------|
| **Bug-A** | Phase | ClairS-TO PASS somatic 進 phasing graph anchor → graph 把 somatic 全綁同一 phase block | `--pon-only-phasing` flag (commit 8b8c1fd) |
| **Bug-B** | Tag (getVote) | SOMATIC-first 投票 → 一個 somatic vote 蓋過所有 germline votes | two-layer 順序反轉 (commit 41ff147) |

### 4.2 為什麼之前報告誤把 Bug-B 描述成「無 fallback」

**用戶提問**揭露了報告描述偏差。錯誤敘述源於：
1. 看 V5 commit message 寫「Layer 1.5 somatic fallback」誤以為 baseline「無 fallback」
2. 沒實際 `git show 8b8c1fd^` 看 baseline 真實 `getVote()` 內容
3. baseline 的 SOMATIC-first 機制其實**就是 fallback**，但**順序錯**

**正確說法**：baseline 有 fallback（even more aggressive than V5 — somatic 是第一優先），但 fallback 順序錯了，造成 99.9% reads 被 somatic 標記覆蓋 germline 標記。

### 4.3 Phase vs Tag 階段如何同時造成 17.3:1

```
ClairS-TO somatic
       │
       ▼
[Bug-A] Phase 階段: somatic 進 graph anchor
       │
       │  → 同 clone reads 共現 ALT
       │  → graph 把 somatic 全推到同一 phase block (HP1)
       │  → phased VCF 中 somatic site 被標 GT=0|0 (with GT2 偏 HP1)
       ▼
[Bug-B] Tag 階段 (baseline getVote): SOMATIC-first 投票
       │
       │  → read 涵蓋 somatic site → countMap[HAPLOTYPE1_1] > 0
       │  → 第 1 優先 (HAPLOTYPE1_1, HAPLOTYPE2_1) 直接決定 hpResult = 11
       │  → 即使 read 也涵蓋 germline 證據（HP2 votes），也被 somatic 蓋掉
       ▼
17.3:1 bias 形成
```

→ 兩個 bug 互相強化：phase 偏 HP1 → tag 也用 somatic 投票偏 HP1 → bias 雙重放大。

---

## §5 對既有 audit suite 結論的影響

### 5.1 結論不變（17.3:1 bias 是真實 bug，V5 修復有效）

| 結論 | 是否仍成立 |
|------|:----:|
| baseline 17.3:1 self-phasing bias | ✅ |
| V5 修復 bias 至 ≈ 1:1 | ✅ |
| V5 不傷 ClairS-TO calling | ✅ |
| V5 sanity check 全 PASS | ✅ |
| 0.6 simulation：baseline self-phasing 衰減 | ✅ |
| 0.6 simulation：V5 conservative tagging | ✅ |

### 5.2 描述精確化（既有報告需修正的點）

| 報告 | 原描述 | 精確化更正 |
|------|--------|----------|
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/11_code_issues_inventory.md` Bug 1-1 | 「baseline 沒有 Layer 1.5 somatic fallback」 | baseline 有 fallback（甚至 somatic-first），V5 是**順序反轉**為 germline-first |
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md` §3.3 | 「Baseline `getVote()` 只有 Layer 1 (germline only)」 | baseline 是三優先序（somatic, HP3, germline），但**順序為 somatic-first**；V5 改為 germline-first 兩層 |
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/14_user_report_integrated.md` §6 | 「baseline 的 phasing graph 把 somatic 也當 anchor」 | 正確（這是 Bug-A），但需補充 Bug-B：baseline tag 階段也是 somatic-first 投票（雙重 bug） |
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/15_software_engineering_perspective.md` §2.1 | 「baseline 單體 if-else，無 fallback」 | baseline 有完整 fallback 機制（用 `vector<pair>` for-loop iteration），V5 是**重排優先序 + encoding 拆分**而非「加 fallback」 |

### 5.3 Bug-A vs Bug-B 對 17.3:1 的相對貢獻

從 commit log 可推得：
- **8b8c1fd** 加 `--pon-only-phasing` 後（只修 Bug-A），出現 99.9% HP:i:21 issue → 暗示 Bug-B 的影響在 PON-only mode 下被放大暴露
- **41ff147** 修 Bug-B 後（順序反轉），99.9% HP:i:21 消失

→ 真 baseline-only mode 下：**Bug-A 主導 17.3:1**（phase 階段偏向已決定 reads 走向）；Bug-B 在 PON-only 開啟後才成為主要問題。
→ V5 必須**雙修**才能修復完整。

---

## §6 對 0.6 simulation 結論的影響

### 6.1 結論不變

`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` 的結論完全不受影響：

- ✅ baseline @ 0.6 self-phasing 強度減弱（chr19-22 ratio 1.33:1 → 1:1.14）
- ✅ V5 在 0.6 仍 conservative tagging（HP33 比例 12.4%）
- ✅ V5 vs baseline gap 在 0.6 縮小

### 6.2 補充解讀

V5 的 **conservative tagging（HP33 多）** 在新理解下更合理：
- baseline 在 LOH 區若 read 同時覆蓋 germline + somatic site：
  - **somatic vote** 強行覆蓋 germline → HP:i:11 或 21（不見得正確）
- V5 在同樣情況：
  - **germline vote first** → 若 germline 投票一致 → HP:i:11/21 (確定)
  - **若 germline=0, somatic≠0** → fallback → HP:i:11/21 (仍打但保守)
  - **若 somatic 投票分散** (HP1_1=HP2_1) → HP:i:33 (ambiguous)

→ V5 把 ambiguous 顯式標 HP33，這是 baseline 沒有的設計益處。

---

## §7 軟工視角的精確化更新

### 7.1 Baseline `getVote()` 真實設計評估

修正 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/15_software_engineering_perspective.md` §9.1：

| 維度 | Baseline 真實狀態 | V5 |
|------|:----:|:--:|
| Sub-genotype 機制存在 | ✅ 有（GT2/GT3 解析 + HAPLOTYPE1_1/HAPLOTYPE2_1）| ✅ 有 |
| Fallback 設計存在 | ✅ 有（3 層 priority via for-loop） | ✅ 有 |
| 投票順序合理 | ❌ SOMATIC-first（造成 99.9% reads HP:i:21） | ✅ Germline-first（正確） |
| Encoding 與 voting 分離 | ❌ 直接給 hpResult = 11/21/33 | ✅ 先 germlineResult 再 encoding |
| 可推回 germline 來源 | ❌ 不行 | ✅ germlineResult 保留 |

### 7.2 修正：V5 是「順序反轉 + encoding 拆分」而非「新增 Layer」

更精準的 SE 描述：

```
Baseline:
  for (priority in [(HP1_1, HP2_1), (HP3, HP2_1), (HP1, HP2)]):  # somatic first
    if has_vote(priority):
      hpResult = direct_encode(winner)  # 直接給 11/21/33
      break
   ↓ Bug: somatic-first 順序 + encoding 與 voting 耦合

V5:
  # Phase 1: germline-first voting
  if has_germline_vote:
    germlineResult = winner_of(HP1, HP2)
  elif has_somatic_vote:  # fallback
    germlineResult = winner_of(HP1_1, HP2_1)

  # Phase 2: encoding (separated)
  if has_somatic:
    hpResult = encode(germlineResult)  # 11/21/33
  else:
    hpResult = germlineResult  # 0/1/2
   ✓ 順序反轉 + 兩階段拆分
```

### 7.3 SE 原則對應的精確化

| SE 原則 | Baseline 真實狀態 | V5 改進 |
|---------|------|--------|
| **Single Responsibility** | for-loop 把 priority+encoding 混在一起 | 兩階段：voting 一階段，encoding 一階段 |
| **Strategy Pattern** | priority list 是 strategies | 改為 layered chain (Layer 1 → 1.5 → 2) |
| **Default Behavior** | 所有 reads 預設走 somatic-first（aggressive） | 預設走 germline-first（conservative） |
| **Code Smell** | encoding logic 散落在 priority winner determination | encoding logic 獨立 |

---

## §8 給後續開發者的 takeaways（更新）

1. **「沒有 fallback」是錯誤推論**：baseline 有 fallback，問題是**順序錯**
2. **Code review 必查 priority order**：當 priority list 涉及 confidence levels，必須確認最可信的在第一順位
3. **Encoding 與 voting 分離有用**：V5 的 germlineResult 是 intermediate state，可在 encoding 前 inspect
4. **commit message 是黃金**：commit `41ff147` message 已明確記載 root cause，但需要實際讀程式碼才能完整理解
5. **用戶問題揭露盲點**：本文件由用戶問題觸發，提醒「程式碼層面的精確描述比抽象總結重要」

---

## §9 跨檔索引

### 程式碼引用（baseline 真實內容）

| 檔案 | 行號 | 函式 | 位置查法 |
|------|------|------|--------|
| `HaplotagProcess.cpp` (baseline 8b8c1fd^) | getVote 全函式 | priority for-loop | `git show 8b8c1fd^:HaplotagProcess.cpp` |
| `HaplotagProcess.cpp` (V5 working tree) | 512-563 | layered getVote | 直接 Read 即可 |
| `HaplotagProcess.cpp` | 160-189 | GT2/GT3 解析 (baseline + V5 共用) | 直接 Read |
| `Util.h` | 24-28 | HAPLOTYPE enum (baseline + V5 共用) | 直接 Read |

### 受影響的 audit suite 文件

| # | 文件 | 影響 |
|---|------|------|
| 11 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/11_code_issues_inventory.md` | Bug 1-1 描述精確化 |
| 13 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md` | §3.3 baseline `getVote()` 描述更正 |
| 14 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/14_user_report_integrated.md` | §6 補充 Bug-A + Bug-B 雙重 bug 說明 |
| 15 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/15_software_engineering_perspective.md` | §2 Baseline 真實程式碼補充 |

---

## §10 結論：對用戶問題的明確答覆

> **「確認當時 baseline 是否將判斷 somatic 放到其他 GT2 GT3 去判讀？」**

✅ **是**，baseline 完整支援 GT2/GT3 sub-genotype 判讀：
- VCF 輸出 GT2/GT3 欄位（baseline 與 V5 都有）
- `HaplotagProcess.cpp:160-189` 解析 GT2/GT3 寫入 `HAPLOTYPE1_1`/`HAPLOTYPE2_1`/`HAPLOTYPE3`
- `getVote()` 用 sub-genotype haps 投票

> **「會影響 baseline phase 與 tag 結論嗎？」**

⚠ **部分影響 — 結論方向不變但根因更精確**：
- 17.3:1 bias 仍然是真實 bug（不變）
- 但 bug 是「**雙重**」：Phase 階段（Bug-A: somatic 當 anchor）+ Tag 階段（Bug-B: SOMATIC-first 投票）
- V5 的修法是「**順序反轉 + encoding 拆分**」而非「加新 fallback layer」
- 既有 audit suite 結論（17.3:1, V5 修復, 0.6 conservative tagging）**全部保持成立**

> **「會影響之前觀察的問題結論嗎？」**

✅ **不會** — 觀察到的 phenomenon（17.3:1, V5 修復, 0.6 self-phasing decay）都是真實的，數據與圖示無誤。只是程式碼層級描述需要精確化（已在本文件 §5.2 列出更正項）。
