<!--
build_date: 2026-05-13
agent: V6 attribution errata (從 PPT v1.6/v1.7 誠信檢查產生; b/c/d 三輪 attribution audit)
status: validated
report_class: errata-companion (修訂 5/8 整合報告 V6 attribution 模糊處)
audience: PI / lab member / 自己未來
parent_report: InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md
inputs:
  - InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md (5/8 整合報告原文)
  - InterSubMod/research/paired_priority_bug_audit/07_V6_validation_findings.md (V6 audit 實測)
  - InterSubMod/research/paired_priority_bug_audit/06_V3F_vs_V5_evaluation.md (V3F vs V5 比較)
  - PPT v1.7 修正歷史 (preview/slide_08, slide_09d, slide_12, slide_17, REHEARSAL_CHEATSHEET)
outputs:
  - 本檔（新 erratum companion, 補強 V6 attribution）
  - 5/8 整合報告頂部加 erratum banner 引用本檔
  - PPT v1.6/v1.7 改動已套用至 preview/*.html (見 §5)
verdict: 3 個 V6 attribution errata；5/8 報告**主要結論不撤回**（priority bug 修對機制因果確立、V6 production candidate）；attribution 模糊處補強，禁止「V6 修對 17.3:1 / +13.3 pp / 34,855」單獨 framing
last_verified: 2026-05-13
report_template: errata-companion v1.0
-->

# Errata for 5/8 整合報告 — V6 Attribution 補強（b/c/d 三輪歸功 audit）

## 0. TL;DR

5/8 整合報告 [`InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`](20260508_Self_Phasing_完整觀察整合報告_01.md) **主要結論不撤回**：

- ✅ priority bug 修對機制因果確立（17.3:1 全基因組 + 34,855 read-level victims）
- ✅ V6 = production candidate (融合 V3F 保守 + V5 marker eng + Layer 1.5 修補)
- ✅ 20 指標 no regression / caller F1 invariant

本 erratum 補強 **3 處 V6 attribution 模糊**（避免讀者誤以為 V6 是修對 17.3:1 / +13.3 pp / 34,855 的單獨功勞）：

| Errata | 議題 | 5/8 報告現狀 | 補強後 |
|---|---|---|---|
| **E1** | 17.3:1 → 1:1 修對 | §2.1 / §9.x 表沒指版本 | V3F 1.14:1 是 ratio 最佳；V5 退步至 2.00；V6 改善至 1.84 (vs baseline 17.3) 但仍偏 V3F；V6 換取 marker eng 改善是 trade-off |
| **E2** | +13.3 pp paired GT @ 0.93 | TL;DR ② 已正確標 V5；但易被讀為 V6 成就 | V5 (V3F + Layer 1.5) vs baseline 達成；V6 重用 V5 phased VCF 預期保留但**未重跑 15-site Clean PS metric** |
| **E3** | 34,855 100% 修對 | §0/§4/§6 已正確標 V3F/V5；V6 未提 | V3F (commit 41ff147 tagging fix) 是修對主力；V5/V6 在此 germline-existent 子集 **logic 繼承 V3F**（V6 唯一改動 Layer 1.5 revert 對 germline-absent 子集，與 34,855 不重疊）；V6 audit 沒重跑 forensic 但 logic 推論 V6 = V3F = V5 valid |
| **E4** ⚠ | judgeHaplotype enum 數值錯誤 | PPT slide 07 line 345 寫「`HAPLOTYPE1_1=2` enum 與 HP tag int `11` 不一致」(2026-05-13 v1.6-E audit 發現)| 實際 `HAPLOTYPE1_1 = 3` (Util.h:24 確認, 不是 2)。dead code 機制不變: enum=3 被當 int 與 hpResult int 11 比較永遠 ≠ → if 永遠 false → HP:i:33 dead code 是**長期未經 read-level tag audit 因而未生效**（非 longphase 作者疏失，是該 codebase 著重 caller F1 / phase block 等 caller-level metric，較少 audit 個別 HP tag 行為）。**修正範圍**: PPT slide 07 (v1.6-E HTML 改完) + 不影響 5/8 報告主結論 (報告 §3.4 機制描述未誤用此數值) |

證據來源：
- V6 audit 12 個 md 0 hits for {`13.3`, `Clean PS`, `15-site`, `88.2`}（V6 未實測 +13.3 pp）
- V6 audit `07_V6_validation_findings.md` §3.1 「全基因組 hp=1-1:hp=2-1 ratio V6=1.838 vs V3F=1.138」
- 5/8 報告 line 674：「17.3:1 → ~1:1 修正功勞**主要在 V3F**」（已寫對但 §2.1 等地方 attribution 模糊）
- 5/8 報告 line 64：「**V3F/V5 修正率均 100%**」（V3F/V5 明確；V6 attribution 未提）

---

## 1. E1 — §2.1 / §9.x「HP1:HP2 17.3:1 → ~1:1 修對」攻略歸 V3F 主力

### 1.1 5/8 報告中 attribution 模糊處

- **line 88-96 §2.1 表**：「| HP1:HP2 | 17.3:1 | 1:1 | 17.3× |」(僅列 baseline vs 隨機，沒列哪版修對到 1:1)
- **line 604 §7.1 表**：「| HP1:HP2 ratio | 17.3:1 | ~1:1 | **消除偏移** |」(沒指版本)
- **line 826 §9.x 表**：「| HP 結構 | HP1:HP2 ratio | 17.3:1 | ~1:1 | 消除 |」(沒指版本)
- **line 684 §7.2 表**：「| HP1:HP2 17.3:1 → ~1:1 | V3F tagging fix（`41ff147`）| 應持平 |」(✅ 已對)

### 1.2 補強表（精確版本歸功）

**全基因組 hp=1-1:hp=2-1 ratio**（priority bug 主指標, HCC1395 5kHz purity 0.93）：

| 版本 | ratio | HP1 占比 | vs baseline 改善 | 歸功來源 |
|---|---:|---:|---:|---|
| baseline | 17.3 : 1 | 94.6% ❌ | — | PI 4-29 報告 |
| **V3F** ⭐ | **1.138 : 1** | **53.2%** ✅ ≈ 1:1 | **大幅修對** | commit 41ff147 tagging fix |
| V5 | 2.003 : 1 | 66.7% ⚠ | 中等 (Layer 1.5 反引偏移) | V5 commit chain |
| V6 | **1.838 : 1** | **64.8%** | -29.8 pp vs baseline (改善但比 V3F 退步) | V6 重用 V5 phased VCF + Layer 1.5 revert |

### 1.3 補強 attribution

- 「**17.3:1 → ~1:1 修對主力為 V3F**（commit 41ff147 tagging fix）」
- V5 加 Layer 1.5 後 ratio 退步至 2.00（germline-absent 區重新引入偏移）
- **V6** 移除 Layer 1.5 後改善至 1.84:1（相比 baseline -29.8 pp；相比 V3F 仍略偏 HP1）
- V6 在 ratio 上是 trade-off：以 ratio 部分退步換取 marker coverage +9.0% / hp=33 +4.7% / Layer 1.5 缺陷修補
- **不能單獨說「V6 修對 17.3:1」** — 應標明「V3F 主力修對；V6 在多維度價值 (marker eng) 改善 vs V3F」

---

## 2. E2 — TL;DR ② / §X.X「+13.3 pp paired GT @ 0.93」歸 V5；V6 未實測

### 2.1 5/8 報告中（已對但易誤讀）

- **line 65 TL;DR ②**：「V5 vs baseline 全表...purity 0.93 sample **+13.3 pp paired GT**」(✅ V5 attribution)
- **line 687 §7.2 表**：「+13.3 pp paired GT concordance = V3F (tagging fix) + V5 Layer 1.5¹」(✅ V3F + Layer 1.5)
- **line 824 §9.x 表**：「15-site Clean PS (11) | 74.9% | 88.2% | +13.3 pp |」(只列數字無版本)
- **line 893 §X.X**：「V5 真實價值在 read-level tag concordance（**+13.3 pp paired GT @ 0.93**）」(✅ V5)

→ 5/8 報告 **+13.3 pp = V5 達成是明確的**。但 PPT 早期版本 (09d 鐵證 1 / 17 verdict ③) 把它寫成 V6 真實價值 — 此 errata 對齊 PPT 修正。

### 2.2 V6 對此 metric 的真實 status

V6 audit `07_V6_validation_findings.md` **0 hits** for {`13.3`, `Clean PS`, `15-site`, `88.2`} — V6 **沒重跑 15-site Clean PS paired GT** metric。

V6 重用 V5 phased VCF（07_V6 §1.3）：
- phase block 結構繼承 V5 → paired GT concordance **理論上保留** ≈ V5 88.2%
- 但 V6 移除 Layer 1.5 → 部分 reads 從 hp=11/21 移到 hp=33（hp=33 在 paired GT 算 ambiguous，可能不算 correct 也不算 wrong）
- → V6 paired GT 確切值 **未知**

### 2.3 補強 attribution

- 「**+13.3 pp paired GT @ 0.93 由 V5 (V3F + Layer 1.5) 達成**」
- V6 重用 V5 phased VCF → 預期保留（但**未實測**）
- Layer 1.5 移除對 read-level GT 的影響為**未來研究**範圍
- **不能說「V6 +13.3 pp」** — 應說「V5 vs baseline +13.3 pp；V6 預期保留但未實測」

---

## 3. E3 — §4 / §6「34,855 victims 100% 修對」歸 V3F；V5/V6 logic 繼承

### 3.1 5/8 報告中（已對但 V6 attribution 未提）

- **line 64 TL;DR ①**：「34,855 victims，**V3F/V5 修正率均 100%**」(✅ V3F/V5)
- **line 262 §4.1**：「全 752 條都是 baseline=11 → **V3F=21 → V5=21**」(✅ V3F=V5)
- **line 564 §6.2**：「**100% V3F+V5 修正**」(✅ V3F=V5)
- **line 572 §6.2**：「**V3F+V5 修對在 read-level 鐵證確立**」(✅ V3F+V5)
- **line 701 §X.X**：「機制因果不變（priority bug 仍是真，**V3F/V5 仍 100% 修對**）」(✅ V3F/V5)
- **line 1240 verdict**：「**100% V3F/V5 修正**」(✅ V3F/V5)

→ 5/8 報告 **34,855 = V3F/V5 100% 修對是明確的**。但 PPT 早期版本 (08 title / 09d 鐵證 2 / 17 verdict ① / cheatsheet 必背 #2) 把它寫成「V6 100% 修對」— 此 errata 對齊 PPT 修正並補 V6 logic 推論。

### 3.2 V6 對此 metric 的真實 status

V6 audit `07_V6_validation_findings.md` 沒重跑 34,855 read-level forensic；只測過 chr19 V5 → V6 zero-sum transfer 2,542 reads（07_V6 §2.1）。

**Logic 推論**（V6 = V3F = V5 valid）：
- 34,855 victims 屬「germline + somatic 都 >0」（germline-existent）子集
- V6 唯一改動 = 移除 V5 Layer 1.5（germline-absent 邏輯）
- → V6 改動對 germline-existent 子集**不適用** → V6 在此子集繼承 V5 phased VCF → V5 在此子集 = V3F 修對結果（Layer 1.5 對 germline-existent 也不改變 hp）
- → V6 = V3F = V5 = 100% 修對在 34,855 子集 **logic valid 推論**

### 3.3 補強 attribution

- 「**34,855 victims 100% 修對由 V3F (commit 41ff147 tagging fix) 主力達成**」
- V5/V6 在 germline-existent 34,855 子集**logic 繼承 V3F 修對結果**
- V6 audit 沒重跑 34,855 forensic（因 logic 上不需要 — V6 改動對此子集不適用）
- **可以說「V3F/V5/V6 三版在此子集 100% 修對」** — 但**必須標明 V3F 是主力**
- **不能說「V6 100% 修對 34,855」單獨 framing** — 給 PI 錯印象 V6 是修對版本

---

## 4. 對 PI 報告 (4-29) errata 的關聯

本 erratum (5/13) **不取代** [`20260509_PI_Report_4_29_Errata_01.md`](20260509_PI_Report_4_29_Errata_01.md)（5/9 PI 報告 4-29 errata, 5 條 E1-E5）。兩者並存：
- **5/9 errata (E1-E5)**：修訂 PI 報告 **4-29** 5 處表述（chr19 hotspot / V5 commit / 證據鏈 / V5=Pass1 only / Layer 1.5 設計缺陷）
- **5/13 errata (E1-E3)**：補強 5/8 整合報告 V6 attribution 模糊（17.3:1 / +13.3 pp / 34,855 三處）

兩者交集 = V5 attribution；本 5/13 erratum 進一步區分 V3F vs V5 vs V6 三版功勞。

---

## 5. 對應 PPT (`InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/`) 改動清單

| Slide / 檔 | E1 (17.3:1) | E2 (+13.3 pp) | E3 (34,855) |
|---|---|---|---|
| `preview/slide_03_genome_173.html` | — (S1 純揭露問題，S5 才收尾) | — | — |
| `preview/slide_08_quant_evidence.html` | — | — | ✅ title + 表 + arrow + speaker note 4 處 |
| `preview/slide_09d_v6_evidence_summary.html` | — | ✅ 鐵證 1 row + speaker note | ✅ 鐵證 2 row + speaker note |
| `preview/slide_12_no_regression.html` | ✅ 加「呼應 S1 + V3F/V6 trade-off」段 | — | — |
| `preview/slide_17_main_verdict.html` | ✅ verdict ① 表格 + speaker note | ✅ verdict ③ 表格 + speaker note | ✅ verdict ① 表格 + speaker note |
| `REHEARSAL_CHEATSHEET.md` | ✅ 必背 #1 + Q9 | ✅ 必背 #3 + Q10 | ✅ 必背 #2 + Q11 |

實際改動已 apply 至 PPT preview files（v1.6 / v1.7 兩輪 commit）。

---

## 6. 對未來工作影響

### 6.1 短期（PI 報告當下）

- 用戶報告 PI 時使用 [`InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/REHEARSAL_CHEATSHEET.md`](../../presentations/validated/2026/05/self_phasing_synthesis_PI/REHEARSAL_CHEATSHEET.md) 已含完整 attribution 速答 (Q9/Q10/Q11)
- 5 大必背數字已標明 attribution

### 6.2 中期（V6 跨樣本擴展 + Phase D 7 樣本）

V6 vs V3F vs V5 三版在 ratio metric 上的差異 (V6 1.84 vs V3F 1.14) 需要 mechanistic 解釋：
- germline-existent 區 V5 ploidy / threshold 細節為何讓 V6 ratio 比 V3F 偏 HP1？
- 此差異是否 sample-specific？(Phase D 4 樣本 ratio 0.61-1.24 中性，HCC1395 是 outlier?)
- 未來研究方向：重跑 V6 15-site Clean PS metric 確認 paired GT 是否保留 ≥V5

### 6.3 長期（report 命名規範）

未來所有「V vs V」對比表格 / verdict 數字 / 必背 anchor，必須**明示版本**（V3F / V5 / V6 / 三版同）。
- ❌ 禁止「~1:1 修對」「100% 修正」「+13.3 pp」無版本 anchor
- ✅ 必寫「**V3F** 1.14:1 ratio 最佳」「**V5** +13.3 pp」「**V3F/V5/V6** 在 34,855 子集 100% 修對 (V3F 主力)」

---

## 7. References

- 5/8 整合報告原文：[`20260508_Self_Phasing_完整觀察整合報告_01.md`](20260508_Self_Phasing_完整觀察整合報告_01.md)
- 5/9 PI 報告 4-29 errata：[`20260509_PI_Report_4_29_Errata_01.md`](20260509_PI_Report_4_29_Errata_01.md)
- V6 audit findings：`InterSubMod/research/paired_priority_bug_audit/07_V6_validation_findings.md`
- V3F vs V5 比較：`InterSubMod/research/paired_priority_bug_audit/06_V3F_vs_V5_evaluation.md`
- PPT preview：`InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/preview/`
- Rehearsal cheatsheet：`InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/REHEARSAL_CHEATSHEET.md`

---

## 8. Decision Log

| Date | Decision | 理由 |
|---|---|---|
| 2026-05-13 | 不撤回 5/8 報告主結論 | priority bug 修對機制 + V6 production candidate 兩主結論在 V6 audit 多維度實測下成立 |
| 2026-05-13 | 不重做 5/8 報告（不動原文）| validated 報告原則 + 5/8 報告 attribution 多數已對（line 64/262/564/572/684/701/1240）；模糊處透過本 erratum 補強 |
| 2026-05-13 | 補 V6 對 34,855 的 logic 推論而非要求重跑 forensic | V6 audit 邏輯上不需重跑（V6 改動對此子集不適用）；重跑為 future work 但非當下必要 |
| 2026-05-13 | +13.3 pp 標 V6 預期保留但未實測 | V6 重用 V5 phased VCF → 推論保留；但精確值 V6 audit 沒測 → 誠信標「未實測」|
| 2026-05-13 | E1 ratio attribution 用 trade-off framing | 避免貶低 V6；保留 V6 production candidate 結論；明確標明 V3F ratio 最佳 / V6 marker eng 改善 |
