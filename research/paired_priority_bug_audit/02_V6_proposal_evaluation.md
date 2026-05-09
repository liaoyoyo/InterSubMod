<!--
build_date: 2026-05-10
agent: V6 proposal evaluation (V5 + V3F hybrid for germline-absent region)
status: validated
report_class: design-decision-evaluation
audience: PI / lab member / 自己未來
parent_synthesis: InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md
parent_step_d: InterSubMod/research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md
inputs:
  - longphase-to-mod HaplotagProcess.cpp:512-563 (V5 getVote 三層)
  - InterSubMod/src/core/ReadParser.cpp:121-153 (HP tag demotion + germline_hp_only flag)
  - InterSubMod/src/core/PerCpgAsm.cpp:244-247 / SignificanceAnalyzer.cpp:288-290 / LabelTest.cpp:230-282 (HP1/HP2 family 計算)
  - 2026-04-21 ReadParser germline-hp-only Phase 1 audit (cycle 20260423_weekly_report)
outputs:
  - 本檔（V6 提案評估）
verdict: NOT NEEDED — ISM 端已有 `germline_hp_only` flag 提供 V5 Layer 1.5 在 germline-absent 區域的等效 demotion；Phase 1 audit 已驗證對 filter F1 無實質改善；HPFineNGroups ≥3 marker 需用 flag=on 重評估
last_verified: 2026-05-10
report_template: design-decision v1.0
decision: 不做 V6 binary patch；改用既有 ISM `--germline-hp-only` flag 在 characterization 層配置；HPFineNGroups marker 重評估列為 follow-up
-->

# V6 提案（V5 + V3F hybrid for germline-absent）評估 — NOT NEEDED

## 0. TL;DR

> **V6 binary 不需做** — ISM 端 `ReadParser.cpp:150` 已有 `germline_hp_only` config flag，當 enabled 時把 hp=1-1/2-1/3 reads demote 為 unphased，**等效於 V5 Layer 1.5 移除後的 V3F 行為**。Phase 1 audit (2026-04-21) 對 HCC1395 TO 全量驗證該 flag 對 filter AUC **無 ≥+0.02 delta 改善**；機制正確但 filter 方向 CONDITIONAL NEGATIVE。**改 longphase-to 端做 V6 binary = 重複勞動**。但 5/9 Step D 揭露的 V5 Layer 1.5 設計問題對 **HPFineNGroups subclone marker** 結論有實質影響（需重評估）。

---

## 1. 用戶提問（5/10）

> 「V5 主結構 + germline-absent 區域回 V3F 設計，這樣不會造成其他部分的結果變遭嗎？還是說只會改善有缺陷的部分？」
> 「合理的邏輯 tag HP33 才可以清楚的讓 ISM 更明確知道 read 預設狀況 是否合理」

評估這個 V6 提案是否值得做。

## 2. ISM 端已有等效機制（關鍵發現）

### 2.1 `ReadParser.cpp:121-153` HP tag demotion 邏輯

```cpp
// Extract HP tag from BAM (HP:i: 整數 or HP:Z: 字串)
std::string hp_raw = "0";
uint8_t* hp_aux = bam_aux_get(b, "HP");
if (hp_aux) {
    char type = hp_aux[0];
    if (type == 'Z' || type == 'H') {
        hp_raw = bam_aux2Z(hp_aux);  // longphase-s 字串
    } else if (type == 'c' || ...) {
        // longphase-to 整數 → 字串
        switch (hp_int) {
            case 1:  hp_raw = "1";   break;
            case 2:  hp_raw = "2";   break;
            case 11: hp_raw = "1-1"; break;
            case 21: hp_raw = "2-1"; break;
            case 33: hp_raw = "3";   break;
        }
    }
}
info.hp_tag_raw = hp_raw;  // 永遠保留原值供 audit

// Self-phasing fallback: treat somatic HP tags as unphased
if (config_.germline_hp_only && (hp_raw == "1-1" || hp_raw == "2-1" || hp_raw == "3")) {
    info.hp_tag = "0";  // ← demote 為 unphased
} else {
    info.hp_tag = hp_raw;  // default: 保留 1-1/2-1/3
}
```

→ **`--germline-hp-only` flag 等效於 V6 提案的 ISM 後端版本**：
- V5 BAM 內 hp=11/21（含 4.19:1 偏移污染） → ISM `info.hp_tag = "0"` → 不算進 HP1/HP2 family
- 達成「V3F 標 hp=33 → ISM 不算進 family」的同樣效果

### 2.2 ISM 內 HP1/HP2 family 計算邏輯

從多處 source code 確認 ISM 對 hp 值的分類：

| Source | line | HP1_family | HP2_family | 不算 family |
|---|---|---|---|---|
| `PerCpgAsm.cpp:244-247` | hp ∈ {"1","HP1","1-1"} | hp ∈ {"2","HP2","2-1"} | hp = "3" / "0" |
| `SignificanceAnalyzer.cpp:288-290` | 同上 | 同上 | 同上 |
| `LabelTest.cpp:230-282` | 同上 | 同上 | 同上 |

→ **default mode（`germline_hp_only=false`）** 下：
- V5 BAM 的 hp=11 reads → ISM hp_tag="1-1" → **算進 HP1_family**（含 4.19:1 偏移污染）
- V5 BAM 的 hp=21 reads → ISM hp_tag="2-1" → **算進 HP2_family**（同樣污染）

→ **flag=on 模式（`germline_hp_only=true`）** 下：
- V5 BAM 的 hp=11/21/33 reads → ISM hp_tag="0" → **不算進任何 family**
- 等同於 V3F 標 hp=33 對 ISM 的處理效果

## 3. V6 binary patch vs ISM flag 的等價性矩陣

| 處理方式 | longphase-to 端 | ISM 端 | 結果 |
|---|---|---|---|
| V5 default + flag=off | hp=11/21（4.19:1 偏移）| 算進 HP1/HP2 family | ISM 受偏移污染 |
| **V6 binary**（V5 移除 Layer 1.5）+ flag=off | **hp=33** | hp="3" → 不算 family | ISM 乾淨 ✅ |
| V5 default + **flag=on** | hp=11/21（BAM 上仍偏移）| **hp_tag="0" → 不算 family** | ISM 乾淨 ✅ **與 V6 等效** |

→ **V6 binary 與 ISM flag=on 在 ISM features 計算上等價**。

差別只在：
- **BAM-level**：V6 BAM 內 hp=33 vs V5 BAM 內 hp=11/21（外部下游工具如 IGV、samtools stats 行為不同）
- **ISM-internal**：完全相同（兩者都不算進 HP1/HP2 family）

## 4. Phase 1 Audit 已實證 ISM flag 影響（2026-04-21）

對 HCC1395 TO 全量（28,509 TP + 11,606 FP）跑 `germline_hp_only` flag=on vs flag=off 對比：

| 指標 | flag=off | flag=on | Δ |
|---|---:|---:|---:|
| HP-related features 最大 AUC delta | — | — | **無 ≥ +0.02 改善** |
| HPFineNGroups AUC | — | — | **−0.026**（衰退）|
| HPMergedDelta AUC | — | — | **−0.025**（衰退）|
| NHP3 AUC | — | — | **−0.035**（衰退）|
| NGroups ≥3 TP regions | 22,732 + 8,148 | **0 / 28,509** | NGroups ≥3 完全消失 |
| 整體 filter F1 改善 | — | — | **無顯著改善** |

→ **filter 方向 CONDITIONAL NEGATIVE**：機制正確（demote 邏輯 12 unit tests pass、audit 獨立、數學守恆），但對 HCC1395 TO 整體 filter performance 無實質改善。

來源：cycle `20260423_weekly_report` (PROV-V5-001 系列前的探索)；[Phase 1 audit 報告](../../docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_HCC1395_01.md)。

## 5. 為什麼 ISM filter 沒改善？

### 5.1 機制推導

V5 Layer 1.5 reads（hp=11/21）約占 tagged reads ~3% (~560K / 18.9M)。即使這 3% 全部偏 HP1，對整體 ISM filter 的 AUC 影響有限：

- HP_Ratio AUC baseline 0.526（PI §6.2，已接近 random）
- 即使移除 priority bug 污染 reads，HP_Ratio 仍接近 random
- **HP_Ratio 本身不是強 filter feature**（PI 報告 §6.2 結論）

### 5.2 實證對齊

Phase 1 audit 數據：HP-related 18 features 最大 delta AUC 僅 +0.0084（HP1FamilyN）。即使把所有 priority bug 污染移除，filter 性能上界很低。

→ **ISM filter 對 V5 Layer 1.5 偏移不敏感**。

## 6. V6 提案的真實價值（不在 filter，在 characterization）

雖然 filter 方向 NEGATIVE，但 Phase 1 audit 揭露**重要 characterization 議題**：

> Phase 1 顯示 flag=on 時 NGroups ≥3 **完全消失**（0/28,509 TP regions）
> 強烈暗示原 NGroups≥4 訊號來自 HP:i:11/21/33 的額外 label 分類，而非真實 subclone

**HPFineNGroups subclone marker 結論可能是 priority bug artifact**：
- 原 5/8 整合報告（藉由 LOH-constrained phasing discovery memory）將 NG=2 same-hap 視為 sub-clone marker
- 若 NGroups ≥3 訊號**僅在 flag=off (V5 hp=11/21 偏移污染) 時存在**，則該 marker 不是真實生物學訊號
- → **HPFineNGroups subclone marker 需重評估**

## 7. 整體判決

| 問題 | 答案 |
|---|---|
| V6 binary patch 是否需要做？ | **否** — ISM `germline_hp_only` flag 已等效 |
| 移除 V5 Layer 1.5 對 ISM filter 是否有改善？ | **否（CONDITIONAL NEGATIVE）** — Phase 1 audit 已實證無 ≥+0.02 AUC delta |
| ISM 預設應 flag=on 嗎？ | **default=off 保留**（filter 不改善；characterization 用戶可選 on）|
| HPFineNGroups subclone marker 是否需重評估？ | **是 ⭐** — flag=on 下 NG≥3 消失，原 ⭐4 結論可能降級至 ⭐3 (pipeline-dependent) |
| 「合理的邏輯 tag hp=33」對 ISM 更明確嗎？ | **理論上是，實證上不顯著** — 對 filter 無影響；對 characterization 可由 flag 控制 |

## 8. 推薦行動

| ID | 動作 | 優先 | 估時 |
|---|---|---|---|
| **V6-A** | **不做 V6 binary**（既有 ISM flag 等效）| — | 0 |
| **V6-B** | 把 V5 Layer 1.5 設計缺陷 + ISM flag 對應關係寫入 PI errata E5（已 5/10 commit 2553e96）| ✅ DONE | 0 |
| **V6-C** | HPFineNGroups subclone marker 用 flag=on 重評估（priority HIGH，影響 LOH-constrained phasing 主軸論文 thesis）| **HIGH** | 1-2 day |
| **V6-D** | （條件）若 V6-C 確認 marker 是 artifact → MEMORY `project_hpfinengroups_subclone_marker.md` 升級降級 | depends V6-C | 30 min |
| **V6-E** | （條件）若 lab member 主導 paired pipeline 也想要相同 demotion → 確認 longphase-s 也走 ISM 端 flag | LOW | n/a |

## 9. 為什麼 5/8 整合報告 §8.6 沒提 ISM flag 等效？

5/8 整合報告 §8.6 把 V5 Layer 1.5 設計缺陷描述為「需要 V6 binary 才能修」是不精確的 — 5/10 本評估補上釐清：**ISM 端已有等效 user-facing flag**，binary 改動不必要。

§8.6.4 + 5/9 Step D 結論不變（V5 Layer 1.5 在 germline-absent 區域 4.19:1 偏移確實存在），但「修補方式」有兩條路徑：
1. ✅ **ISM `germline_hp_only=true`（既有 flag, 預設 off）**
2. ❌ **V6 binary patch**（多餘，已有等效）

下次主報告 amend 時應補此釐清；本檔可作為 5/8 §8.6 的 V6 evaluation companion。

## 10. 引用文件

- [longphase-to-mod HaplotagProcess.cpp:512-563](file:///big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp) — V5 getVote 三層
- [InterSubMod ReadParser.cpp:121-153](../../src/core/ReadParser.cpp) — HP demotion 邏輯
- [Phase 1 audit memory](.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/project_readparser_germline_hp_only_phase1_negative.md) — flag CONDITIONAL NEGATIVE
- [HPFineNGroups marker memory](.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/project_hpfinengroups_subclone_marker.md) — ⭐4→⭐3 降級依據
- [5/8 Self-Phasing 整合報告 §8.6](../../docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md)
- [5/9 Step D germline-absent finding](01_step_D_germline_absent_finding.md)
- [5/9 PI errata companion E5](../../docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md)
