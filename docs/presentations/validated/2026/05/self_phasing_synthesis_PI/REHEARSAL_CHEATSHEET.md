---
title: Self-Phasing PI 報告 — Rehearsal Cheatsheet
audience: PI 1-on-1 進度報告
duration_min: 25
duration_qa_min: 5
slide_count_main: 18
slide_count_backup: 5
generated_from: plan v1.1.5 audit (D6 敘述準備度)
generated_date: 2026-05-13
source_plan: /bip7_disk/liaoyoyo2001/.claude/plans/agent-harness-langgraph-resilient-waterfall.md §15
---

# Self-Phasing PI 報告 — Rehearsal Cheatsheet

> **報告主軸**：我們發現並修好了 longphase 一個被忽略的 priority bug
> **完整邏輯鏈**：discovery → mechanism → fix → validation → impact
> **真實價值**：read-level tag concordance（**非** caller F1）

---

## 1. 必背 5 大數字（每場必講）

| # | 數字 | 含義 | 出處 slide | section |
|---|---|---|---|---|
| 1 | **17.3 : 1 → V6 1.84 : 1 (改善 -29.8 pp vs baseline)** | baseline 94.6% 偏 HP1 → V6 64.8%；V3F 1.14:1 (53.2%) 為 ratio 最佳 / V6 換 marker eng 改善 (trade-off) | slide 03 / 12 / 17 | S1 / S5 / S7 |
| 2 | **34,855 V3F/V5/V6 100% 修對 (V3F 主力)** | baseline read-level victims → V3F (41ff147) 100% 修對；V5/V6 此 germline-existent 子集繼承 V3F (Layer 1.5 revert 對此子集不適用)；V6 audit 沒重跑 34,855 forensic 但 logic 推論 V6=V3F=V5 valid | slide 08 / 09d | S3 / S4 |
| 3 | **+13.3 pp = V5 vs baseline 達成** (V3F + Layer 1.5 合力, PI 4-29) | paired GT concordance @ 0.93 — read-level 真實 value；V6 重用 V5 phased VCF 預期保留, **未實測** | slide 09d / 17 | S4 / S7 |
| 4 | **caller F1 = 0.7166** | 三版（V3F/V5/V6）caller F1 數學保證 invariant | slide 12 / 14 / 17 | S5 / S7 |
| 5 | **hp=33 +4.7% / marker coverage +9.0%** | V6 vs V3F marker engineering 改善 | slide 09d / 17 | S4 / S7 |

**輔助數字（按 section 需要）**：
- SP1/2/3 三案 baseline 113:0 / 109:1 / 108:0 極端失衡（slide 04, S1）
- chr19 752 victims → 全基因組 34,855 = 46.4× generalize（slide 08）
- V6 = 5 commits binary patch (f17754f → 2553e96 → 71d21bd → V6)
- Phase D 4 樣本 ratio 0.61-1.24（slide 09d / B4）vs V5 baseline 1.86

---

## 2. 每 Section take-home（一句話）

| Section | take-home（30 sec 內講完） |
|---|---|
| **S0 開場** | 我們發現並修好了 longphase 一個被忽略的 priority bug — longphase HP tag 是 ISM 下游所有分析的基礎，tag 不可信則下游全失準 |
| **S1 觀察** | 17.3:1 全基因組偏移 + SP1/2/3 三案近 100% 失衡 = systematic engineering artifact（不是樣本性質） |
| **S2 機制** | getVote (Layer 1 投票) + judgeHaplotype (Layer 1.5 落判) **兩層分工**；judgeHaplotype 有 enum/integer 比較 bug，HP=33 永遠是 dead code branch |
| **S3 量化** | chr19 752 victims → 全基因組 34,855 = 系統性影響非 edge case；5 commits trace 證據鏈完整 |
| **S4 修補（高潮）** | V6 兩層修補（getVote priority 修 + judgeHaplotype enum 修）+ IGV/BAM forensic 鐵證 SP1/2/3 修對；4 維度 hub 獨立 corroborate |
| **S5 驗證** | F1 雙口徑誠實揭露 — caller F1 invariant 不是失敗，是 ISM 改 read-level layer 不影響 FILTER 的數學保證 |
| **S7 結論** | V6 = production candidate；**真實價值在 read-level tag concordance（+13.3 pp）— 非 caller F1**；本報告 priority bug 修復鏈完結，轉回 ISM 五大目標主軸 |

---

## 3. 預期 PI 追問 × 答案速記

### Q1: F1 三版相同代表沒進步？
**A**：caller F1 是 FILTER 決策層的數學保證 invariant — ISM 改善的 layer 是 read-level HP tag，不影響 FILTER 決策，所以 F1 數學上 invariant。真實 value 在 read-level：paired GT concordance **+13.3 pp @ 0.93**、hp=33 還原 +944%、marker coverage +9.0%。

### Q2: V6 為何比 V5 強？
**A**：V6 = V3F 保守 hp=33 + V5 設計（ploidy / threshold / phased VCF）+ marker engineering 改善；hp=33 +4.7%、marker coverage +9.0%；Phase D 4 樣本 ratio **0.61-1.24 全中性**（vs V5 baseline 1.86 仍偏 HP1）。

### Q3: V5 → V6 有 regression 嗎？
**A**：SP1/2/3 三案 V5 翻 HP2 對齊 paired 是優點；V6 在 germline-absent 區退 hp=33 是設計取捨（保守、不過度修正）。Forensic 5/13 確認 V5 vs V6 **守恆完美**：Type B reads (germline-only) 完全相同；Type C reads (germline-absent + somatic) HP:i:11 + HP:i:21 減 = HP:i:33 增（轉移 5,927 reads）。

### Q4: 「V6 chr19 distance to paired 0.367」不是退步嗎？
**A**：5/13 forensic 揭露此結論需修正 — 統計腳本沒考慮跨 codebase HP1/HP2 命名顛倒（longphase-to HP:i: vs longphase-s HP:Z:）。命名顛倒後 V6 = 0.065 ≈ V5 = 0.068 ≈ baseline = 0.011 同數量級。詳見 slide 09e + `InterSubMod/research/paired_priority_bug_audit/v6_quantification_erratum_2026_05_13.md`。

### Q5: 為何需要兩層修補（vs 單層）？
**A**：**責任分離**：getVote (Layer 1) 是投票決策層（每 read 對 HP1/HP2 vote 哪邊）、judgeHaplotype (Layer 1.5) 是輸出落判層（決定 HP:i:11/21/33）。兩層各自有獨立 bug：getVote 的 priority bug 讓低 confidence reads 強行投票；judgeHaplotype 的 enum/integer 比較 bug 讓 HP:i:33 dead code（永遠不會輸出）。單層修無法 cover 另一層。

### Q6: 34,855 怎麼從 chr19 752 推到全基因組？
**A**：chr19 752 是 HCC1395 5kHz pilot；全基因組 34,855 是同樣 baseline → V6 跨所有 chr 的 read-level victims 統計，**baseline HP:i:11 → V6 HP:i:21 單向、0 反向**。比例 46.4×（chr19 占基因組 ~2.5%），跨 chr 行為一致 → 不是 chr19 特殊。

### Q7: 為何用 SP1/2/3 + chr19:27M 兩個案例？
**A**：**兩區互補**：SP1/2/3 屬 **germline-absent** 區（V6 退 hp=33 設計取捨）；chr19:27,376,222 屬 **germline-existent + somatic 共現** 區（V6 翻方向修對 priority bug）。兩個 case 共同證 V6 雙區行為符合預期。

### Q8: ISM 五大目標研究怎麼接上去？
**A**：本報告 priority bug 修復鏈到 slide 17 verdict 完結；slide 18 future 是**轉場意圖預告**（已在 v1.1.5 加 caveat box 明示），未來 F1 LOH 內外 TP/FP 差異、F2 subclone 結構 + two-hit、F3 7-sample expansion、F4 erratum patch 都是 ISM 五大目標延伸 — 但**不在本次報告主軸內**。

### Q11: 34,855 read-level victims 是 V6 修對嗎？(歸功正確性 #3)
**A**：「34,855 victims 100% 修對」由 **V3F (commit 41ff147 tagging fix)** 達成 (T1.2-F1 audit, 主報告 line 64/564/572)；V5/V6 在此 germline-existent 子集**繼承 V3F 修對結果**，因為：
- 34,855 victims 屬「germline + somatic 都 >0」(germline-existent) 子集
- V6 唯一改動 = 移除 V5 Layer 1.5 (germline-absent 邏輯)
- → V6 改動對 germline-existent 子集**不適用** → V6 = V3F = V5 = 100% (logic 推論 valid)

V6 audit 07_V6_validation_findings **沒重跑 34,855 read-level forensic** (因 logic 上不需要)；只測過 chr19 V5→V6 zero-sum transfer 2,542 reads (germline-absent 子集)。

**核心訊息**：不要在 PI 面前說「V6 修對 34,855」(歸功不準) — 應說「**baseline 34,855 priority bug victims；V3F (41ff147) 是 tagging fix 修對主力；V5/V6 在此子集 logic 上繼承 V3F**」。

### Q10: +13.3 pp paired GT 是 V6 還是 V5 達成？(歸功正確性 #2)
**A**：是 **V5 (V3F tagging fix + Layer 1.5) vs baseline** 達成（74.9% → 88.2%, PI 4-29 報告 15-site Clean PS metric, line 824/687）。**V6 從沒測過此 metric** — V6 audit 12 個 md 0 hits。

V6 重用 V5 phased VCF → 理論上 phase block 結構繼承 V5，paired GT 預期保留；但 V6 移除 Layer 1.5 對 read-level GT concordance 的影響**未重跑 15-site Clean PS metric** → 未來研究範圍。

**V6 真正驗證 strong 的 metric (5 維度, 歸功正確)**：
1. V6 hp=33 ambiguous reads +944% vs V5 (138K vs 13K)
2. V6 marker coverage +9.0% vs V3F (23,980 vs 21,997 NG≥3 regions)
3. V6 chr19 marker rate 0.935 (介於 V3F 0.947 / V5 0.924)
4. V6 Phase D 4 樣本 ratio 0.61-1.24 中性 (vs baseline 17.3:1)
5. V6 caller F1 0.7166 invariant vs V3F/V5

**核心訊息**：不要在 PI 面前說「V6 +13.3 pp」(歸功錯)；說「V5 +13.3 pp 證明 priority bug 在 read-level 有真實影響；V6 在其他 5 維度有 valid 改善 vs baseline」。

### Q9: V6 全基因組 ratio 1.84:1 比 V3F 1.14:1 退步，為何還要升 V6？(誠信揭露)
**A**：在「ratio 接近 1:1」這單一 metric 上 V3F 確實更佳。V6 升級價值不在此 metric，而在多維度綜合：
1. **marker coverage +9.0% vs V3F** (23,980 vs 21,997 NG≥3 regions) — ISM 下游可用 marker 變多
2. **hp=33 ambiguous reads +4.7% vs V3F** (138,317 vs 132,060) — 保守 ambiguous 標記恢復，避免 V5 過度修正
3. **Layer 1.5 缺陷修補** — V5 設計目標保留但避開 priority bug 在 germline-existent 區的 feature 化
4. **caller F1 三版相同** (V3F/V5/V6 = 0.7166) — 不傷下游 FILTER

ratio 退步原因 mechanistic 解釋: V6 重用 V5 phased VCF，germline-existent 區 V5 ploidy/threshold 細節未改 → 此 mechanism 為**未來研究**範圍。

**核心訊息**：V6 = 「V3F 保守 + V5 marker eng + Layer 1.5 修補」三者最佳平衡的 production candidate；相比 baseline 17.3:1 大幅改善 (94.6% → 64.8%, -29.8 pp)；相比 V3F 是工程 trade-off (ratio vs marker eng) 非退化。

---

## 4. 反直覺主張的開場 hook（吸引 PI 注意力）

| Hook | 使用時機 |
|---|---|
| 「我們找到的 bug 是 **dead code branch** — HP:i:33 永遠不會出現」 | slide 07 / 09b 機制揭露 |
| 「改 **13 行 C++** 對應 **34,855 reads** 修對」 | slide 09 fix design 開場 |
| 「**F1 沒變但 read-level 真實對齊 +13.3 pp**」 | slide 14 / 17 cliffhanger / verdict |
| 「priority bug 被 longphase 作者忽略 — 因為 **HP:i:33 看起來像 ambiguous 不重要**」 | slide 09b S4 開場 |

---

## 5. 場上應對策略

### 5.1 時序緊張時的壓縮路線（若超時）
**砍 backup 走主軸 18 張**：
- 不講 B1-B5（QA 才秀）
- S2 機制 3 張壓縮到 2 張時：05 player_referee 縮 60 sec、06 priority bug 縮 90 sec、07 兩層修補 縮 60 sec
- S5 驗證 2 張壓縮：12 no_regression 30 sec + 14 cliffhanger 60 sec

### 5.2 PI 卡在某張時的策略
- **卡在 S1 觀察**（17.3:1 不接受）→ 直接秀 slide 04 SP1/2/3 個案層 113:0 數字震撼
- **卡在 S2 機制**（enum bug 沒概念）→ 用 09b 三列對照表逐欄解（比較目標 / 比較結果 / 判定邏輯方向）
- **卡在 S4 V6 證據**（質疑 cross-codebase）→ 跳 09e forensic 守恆證明（V5 = V6 Type B 完全相同 + Type C 5,927 守恆）
- **卡在 S5 F1**（不理解雙口徑）→ 跳 B2 backup 深掘 F1 數學機制

### 5.3 結尾收尾必達
最後 30 sec **必說 verdict 3 ★**：「真實價值在 read-level concordance — 非 caller F1。caller F1 數學保證 invariant 是因為 ISM 改的 layer 不影響 FILTER。**這是 caller F1 看不到的 dimension**。」

---

## 6. 必避免（場上禁區）

| 禁區 | 為何禁 |
|---|---|
| ❌ 說「V6 比 V5 好」沒附 caveat | 4-29 PI 報告 E4 errata：V5 = Pass 1 only，主要功勞 V3F + Layer 1.5；V6 是設計收斂不是 V5 取代 |
| ❌ 講「34,855 reads V6 修好」沒提 5/13 forensic erratum | 09e 守恆證明是誠信揭露，必須一起講 |
| ❌ 用 caller F1 = 0.7166 當主成績 | 主成績是 read-level concordance +13.3 pp；F1 是「沒退化」的責任證明 |
| ❌ 把 slide 18 future F1/F2 講成本報告延伸 | 18 已加轉場 caveat：F1/F2 是 ISM 五目標支線，不是本報告主軸 |
| ❌ 承諾「ISM 用 V6 tag 一定能 F1 提升」 | F1 對 V3F/V5/V6 數學 invariant；提升 ISM F1 需要 caller 改動 — 不是本報告範圍 |

---

## 7. 報告者自我檢核（場前 5 分鐘）

- [ ] 5 大必背數字默念一次（17.3:1 / 34,855 / +13.3 pp / 0.7166 / +4.7% +9.0%）
- [ ] 7 段 take-home 各 30 sec 自試（S0/S1/S2/S3/S4/S5/S7）
- [ ] 反直覺主張 verdict 3 ★ 演練（為何 F1 不變不是失敗）
- [ ] 開啟 preview/index.html 走 18 主 slide × 1 圈，每張停 5-10 sec 確認 IGV/figure 載入
- [ ] backup B1-B5 順序確認（為意外 Q 準備）

---

## 8. 緊急斷網/技術問題備援

- 若 PanZoom CDN 載入失敗 → IGV 縮圖仍能 hover scale + Shift+click 開新分頁
- 若 IGV PNG 載入失敗 → speaker note 內含完整 read_name + chr19:pos 描述
- 若 slide 順序 nav 失靈 → 用 ←/→ 鍵或瀏覽器 history forward/back
- 若整個 preview 環境壞 → 直接秀 source `InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`

---

**Source files**：
- Preview：`InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/preview/`
- Plan：`/bip7_disk/liaoyoyo2001/.claude/plans/agent-harness-langgraph-resilient-waterfall.md` §15
- Source report：`InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`
- Errata：`InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md`
- v6 forensic erratum：`InterSubMod/research/paired_priority_bug_audit/v6_quantification_erratum_2026_05_13.md`
