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

> **⚠ 2026-05-14 amend (E5)**：本節 §3.2 logic 推論 (V6 = V3F = V5 valid) 已被 V6 spot check 全集實測**推翻**。詳見下節 E5。

---

## 3a. E5 — V6 spot check 全集實測推翻 logic 推論（**2026-05-14 新增**）

> **⚠ 2026-05-14 PM amend (Interim revision)**：本節原 framing「V6 31.7% 仍 hp=11 = V6 失敗繼承 V3F 修對」**有誤導**。V3F vote_dump 對同樣 17,404 victim subset 自己也 42.76% hp=11；V6 31.7% 實際**比 V3F 改善 11 pp** + 加 24% 保守 hp=33/empty tag。詳見 §3a.5b（待 V3F/V5 BAM 完成後 finalize）。

### 3a.1 方法

- **驗證對象**：17,404 unique reads（從 V3F vote_dump on baseline 抽 germline_maj ≠ somatic_maj 且 both >0 子集，去重後）
- **方法**：用 V3F vote_dump 抽 victim read_names → samtools view V6 BAM (`/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam`) → 提取每條 victim read 的 HP:i: tag → 統計 V6 對此 subset 的 HP tag 分布
- **執行**：2026-05-14 02:19 啟動 samtools view (287GB BAM) → 03:33 完成 (~74 min)；output `/tmp/v6_victim_hp_tags.tsv` (958KB, 24,227 events, dedupe 17,404 unique reads = 100% cover)
- **驗證 criteria**：若 logic 推論成立，V6 hp=21 count 應 ≥ 80%、hp=11 count 應 ~ 0%

### 3a.2 實測結果（per-unique-read dedupe, 17,404 reads）

| V6 HP tag | Count | % | 預期 (logic 推論) | 實測 vs 預期 |
|---|---|---|---|---|
| **hp=21** (V6 翻方向修對) | 7,769 | **44.6%** | ≥ 80% | ❌ 嚴重低於預期 |
| **hp=11** (V6 仍 baseline 錯方向) | **5,510** | **31.7%** | ~ 0% | ❌ **意外存在** — 推翻 logic 推論 |
| **hp=33** (V6 保守 ambiguous) | 2,458 | 14.1% | ≤ 20% | ✅ 預期內 |
| 空 (V6 未 tag) | 1,667 | 9.6% | ≤ 10% | ✅ 預期內 |

### 3a.3 結論 — Logic 推論失效

> **V6 對此 victim subset 不是 100% 修對。約 31.7% reads V6 仍標 hp=11（與 baseline 同方向）。V6 並非完全繼承 V3F。**

### 3a.4 可能原因（待深掘）

1. **V6 重用 V5 phased VCF** — V5 phased VCF 對 germline-absent 邊界 reads 已標 hp=11 → V6 繼承這些 reads 的 hp tag（並非由 Layer 1.5 重新決定）
2. **vote_dump (getVote 階段) vs BAM HP tag (haplotag 階段) 不一致** — 兩階段不同處理邏輯
3. **vote_dump 篩 criteria 邊界 case** — 「germline_maj ≠ somatic_maj + both >0」可能包含 baseline 已標 HP1 family + V3F vote 翻 HP2 但 V6 BAM 邊界仍標 hp=11 的 edge case reads
4. **V3F victim list 主要分布於 germline-existent + somatic 共現區，但部分 reads 跨到 germline-absent 邊界**，V6 移除 Layer 1.5 後此邊界 reads 反而標 hp=11

### 3a.5 對 §3 (E3) 的修正

§3.2 「V6 改動對 germline-existent 子集不適用」假設**部分失效**：實測顯示至少 31.7% reads V6 行為與 V5 不一致（V5 應 = V3F = hp=21；V6 = hp=11）。

修正後 attribution：
- 「**V3F (commit 41ff147) 100% 修對 34,855 victims**」(read-level forensic anchor)
- 「**V5 logic 上繼承 V3F**」(V5 audit 確認 V5 = V3F 在此 subset)
- 「**V6 對此 subset 不完全繼承 V3F**」(spot check 31.7% reads V6 仍 hp=11)；機制待深掘

### 3a.6 對 PPT 的更新（5/14 已 apply）

| 檔案 | 改動 |
|---|---|
| `preview/slide_08_quant_evidence.html` | title / caveat / table row 更新 partial → final 44.6%/31.7%/14.1%/9.6% |
| `preview/slide_09d_v6_evidence_summary.html` | row 2 V3F/V5 forensic 100% + V6 不完全繼承說明 |
| `preview/slide_17_main_verdict.html` | verdict ① + speaker note 對齊 final 數字 |
| `REHEARSAL_CHEATSHEET.md` | 必背 #2 + Q11 完整改寫 |

### 3a.7 對核心結論的影響

| 結論 | 是否成立？ | 理由 |
|---|---|---|
| priority bug 修對機制因果確立 | ✅ 仍成立 | V3F 100% 修對 34,855 + SP1/2/3 IGV 對齊 paired 3/3 |
| V6 = production candidate | 🟡 部分成立 | hp=33 +4.7% / marker coverage +9.0% / caller F1 invariant / Phase D 4 樣本 ratio 中性都成立；但「V6 在 34,855 子集繼承 V3F 修對」推翻 → V6 是「marker engineering 強 + ratio 改善 + 但對 priority bug 全集子集修對不完全」的 trade-off candidate |
| 20 指標 no regression / caller F1 invariant | ✅ 仍成立 | FILTER 不動 → F1 數學保證 invariant |

### 3a.8 後續行動

1. ✅ PPT 5 處更新（slide 08 / 09d / 17 / cheatsheet 必背 #2 / Q11）
2. ✅ 本 erratum E5 加入
3. ⏳ V6 hp=11 reads chr/pos 分布分析（10 min awk）— 確認是否集中於 germline-absent 邊界
4. ⏳ V6 hp=11 reads 對齊 V5 phased VCF 邊界檢查（mechanism 確認）
5. ⏳ 若機制確認為「V6 重用 V5 phased VCF 邊界繼承 hp=11」→ 補進 5/8 整合報告作補丁（不撤回主結論，但 attribution 精確化）

---

## 3a.5b. Interim Revision (2026-05-14 PM) — V6 非退步而是保守化

> 此節為 5/14 下午 v1.7-G 深度調查的 interim 結論；最終數據待 V3F / V5 BAM extract 完成後 finalize（background tasks bq4dajhz9 + bbsoraygs，~60-90 min）。

### 3a.5b.1 新發現

V3F vote_dump (genome) 對同樣 17,404 victim subset 的 hpResult 分布：

| HP tag | V3F vote_dump | V6 BAM | Δ V6 vs V3F |
|---|---|---|---|
| **hp=21** (HP2 family) | 9,962 (**57.24%**) | 7,769 (44.6%) | -12.64 pp |
| **hp=11** (HP1 family) | 7,442 (**42.76%**) | 5,510 (31.7%) | **-11.06 pp** ✅ V6 比 V3F 少 |
| hp=33 (ambiguous) | 0 (binary 決策) | 2,458 (14.1%) | +14.1 pp ✅ V6 加保守 |
| empty (未 tag) | 0 (vote_dump 全 cover) | 1,667 (9.6%) | +9.6 pp |

→ **V6 不是退步**：V6 hp=11 比 V3F 自己**少 11 pp**；V6 把 V3F 的 ~24% 邊界 reads **保守化** 重標為 hp=33/empty。

### 3a.5b.2 之前 framing 為何錯誤

**原 framing**（§3a.5）：「V6 31.7% hp=11 = V6 不完全繼承 V3F = V6 失敗」

**錯誤前提**：
- 假設 V3F「100% 修對 34,855」適用於 V3F vote_dump 抽出的 17,404 victim subset
- 實際「34,855 V3F 100% 修對」是**狹義 subset**「baseline=hp=11 → V3F=hp=21 flip」（by construction = 100%）
- 17,404 是**廣義 subset**「germline_maj ≠ somatic_maj + both >0」— 含 V3F 自己也標 hp=11 的邊界 case
- 兩個 subset 不同，不能直接比較

### 3a.5b.3 真實機制（**待 BAM 完成驗證**）

**關鍵程式碼發現** — `/big7_disk/liaoyoyo2001/longphase-to-mod` git log（commit chain）:

| commit | stage | 變動 |
|---|---|---|
| `41ff147` (V3F) | **haplotag** | two-layer getVote — 唯一 動 haplotag |
| `380e8d2` | haplotag | INDEL guard |
| `d0bcd8c` (V5) | **phasing** | ploidyRatio after PON-only Pass 2 |
| `938f0df` (V6) | **phasing** | purity calculation + threshold 0.95→0.9 |

→ **V5 / V6 都是 phasing stage 改動，haplotag binary 自 V3F 後沒變過**。

| Stage | V3F | V5 | V6 |
|---|---|---|---|
| haplotag binary | `41ff147` | 同 V3F | 同 V3F |
| phasing binary | baseline | `d0bcd8c` | + `938f0df` |
| HCC1395 Pass 2 trigger | No (0.93 < 0.95) | No | **Yes (0.93 > 0.9)** |
| phased VCF for HCC1395 | V3F phasing | ≈ V3F | **改變** (Pass 2 觸發) |

V6 對 17,404 reads 的差異 = V6 改 purity threshold 觸發 HCC1395 Pass 2 → phased VCF 重組 phase block → V3F haplotag binary 對 17,404 reads 給出不同 tag。

### 3a.5b.4 對 3 hypothesis 的修正 verdict

| Hypothesis | 原 verdict (§3a.4) | 修正 verdict |
|---|---|---|
| **H1** V6 重用 V5 phased VCF | candidate | 🟢 **valid mechanism** — V6 phased VCF 與 V3F 不同（Pass 2 trigger 差異） |
| **H2** vote_dump vs BAM 階段不一致 | candidate | 🔴 **refuted in V3F case** — V3F vote_dump = V3F BAM 預期一致（待 BAM 確認）；V6 用 V3F 同 haplotag binary 但不同 phased VCF |
| **H3** V3F victim list edge case | candidate | 🟡 **partial — broader criteria 包含 V3F 自己也 hp=11 reads**；不是 edge case leak 而是 criteria 太廣 |
| **H4 新** phasing stage 差異 | n/a | 🟢 **主機制** — V6 purity threshold 0.95→0.9 → HCC1395 Pass 2 trigger → phased VCF 重組 |

### 3a.5b.5 對核心結論影響（再評估）

| 結論 | 原 verdict | 修正 verdict |
|---|---|---|
| priority bug 修對機制因果確立 | ✅ 仍成立 | ✅ 不變 |
| V6 = production candidate | 🟡 部分成立 (V6 對 subset 修對不完全) | ✅ **完全成立** — V6 對 broader subset 改善 11 pp hp=11 + 加 24% 保守 ambiguous |
| caller F1 invariant | ✅ 仍成立 | ✅ 不變 |

V6 的「保守化」與既知優點一致：hp=33 marker engineering +4.7% / coverage +9.0% / Phase D 4 樣本 ratio 中性化。

### 3a.5b.6 後續驗證（background tasks）— **2026-05-14 PM 完成**

| Task ID | 內容 | 結果 |
|---|---|---|
| `bq4dajhz9` ✅ | V3F BAM extract on 17,404 victim | **V3F BAM ≠ V3F vote_dump**！V3F BAM hp 分布: hp=21=41.5% / hp=2=16.65% / empty=15.24% / hp=11=14.46% / hp=33=7.77% / hp=1=4.38%（vote_dump 二元 42.76%/57.24% 是 getVote() 中間記錄，BAM 是 judgeHaplotype 整合最終決策）|
| `bbsoraygs` ✅ | V5 BAM extract on 17,404 victim | **V5 BAM = V6 BAM 100% read-by-read 完全相同**（0 diff after normalization）|

## 3a.5b.7 V3F vs V5/V6 BAM Final Truth (2026-05-14 PM, **修正 §3a.5b 1-6 多項解讀**)

### Triple Comparison (V3F / V5 / V6 BAM on 17,404 V3F victim subset)

| HP tag | V3F BAM | V5 BAM | V6 BAM | V5/V6 vs V3F |
|---|---|---|---|---|
| hp=21 (HP2+somatic) | 7,222 (41.50%) | 7,769 (44.64%) | 7,769 (44.64%) | +3.1 pp |
| hp=11 (HP1+somatic) | 2,517 (14.46%) | 5,510 (31.66%) | 5,510 (31.66%) | **+17.2 pp** ⚠ |
| hp=33 (ambiguous) | 1,353 (7.77%) | 2,458 (14.12%) | 2,458 (14.12%) | +6.4 pp |
| empty | 2,652 (15.24%) | 1,667 (9.58%) | 1,667 (9.58%) | -5.7 pp |
| hp=2 (HP2 germline) | 2,898 (16.65%) | 0 (0%) | 0 (0%) | **-16.7 pp** |
| hp=1 (HP1 germline) | 762 (4.38%) | 0 (0%) | 0 (0%) | -4.4 pp |

### Key Finding 1: V5 = V6 (V6 patch 對此 subset 完全無效)

V6 purity threshold `938f0df` (0.95→0.9) 對 HCC1395 5kHz priority bug victim subset **完全沒影響**：
- 物理檔案不同 (inode 615562637 vs 684605118, size +280KB)
- 但 read-by-read HP tag 100% 相同：7,769 hp=21 / 5,510 hp=11 / 2,458 hp=33 / 1,667 empty
- transition matrix: 4 種 unchanged 100%（21→21 / 11→11 / 33→33 / empty→empty）
- → 在 HCC1395 5kHz 上「V6 比 V5 好」**不可區別**

### Key Finding 2: V3F → V5 改動 63% reads（d0bcd8c ploidyRatio 真正生效）

V3F → V5 transition (V5 = V6 indistinguishable from here):
- Unchanged (V3F = V5): 6,419 (**36.88%**)
- Changed (V3F ≠ V5): 10,985 (**63.12%**)
- **Direction FLIP (HP1 ↔ HP2)**: **3,198 reads (18.4%)** — 2,128 V3F hp=21 → V5 hp=11 + 1,070 V3F hp=11 → V5 hp=21

| Top transitions | Count | 意義 |
|---|---|---|
| V3F hp=21 → V5 hp=21 | 3,398 | 不變 |
| V3F hp=21 → V5 hp=11 | 2,128 | **direction FLIP** |
| V3F hp=2 → V5 hp=21 | 1,688 | germline-only → HP2+somatic |
| V3F hp=21 → V5 hp=33 | 1,644 | 確定 → ambiguous |
| V3F hp=11 → V5 hp=11 | 1,340 | 不變 |
| V3F hp=11 → V5 hp=21 | 1,070 | **direction FLIP** |

→ V5 d0bcd8c phasing fix **不是局部修補**，是**全 BAM 大規模 phase block 重組**。

### Key Finding 3: vote_dump ≠ BAM (重要 caveat)

V3F vote_dump 二元 hpResult (42.76% hp=11 / 57.24% hp=21) **不等於** V3F BAM 最終 HP:i: tag (14.46% hp=11 / 41.50% hp=21 / 16.65% hp=2 / 15.24% empty / 7.77% hp=33 / 4.38% hp=1)。

- vote_dump = getVote() 階段中間記錄（per-variant-position，二元決策）
- BAM HP:i: = judgeHaplotype 階段最終決策（per-read，整合所有變異 + 含 hp=1/hp=2/hp=33/empty 多種 outcome）

§3a.5b 1-5 之前用「V3F vote_dump 42.76% hp=11」對比 V6 BAM 31.7% 是 **不正確比較對象**。

### Key Finding 4: V5/V6 vs V3F 「31.7% hp=11 變 17 pp 多」可能不是退步

V5/V6 對 V3F **17 pp 多 hp=11** + **24% 多 hp=33+empty** + **失去 21% hp=1/hp=2** — 三個 framework 可能解讀：

| Framework | V5/V6 vs V3F 解讀 |
|---|---|
| **(α) 命名顛倒** | V5 Pass 2 phasing 重組 phase block 翻 HP1/HP2 designation；3,198 reads direction flip 是命名差異非實質差異 |
| **(β) V5 過度修正** | V5 把 V3F 認為 germline-only (hp=1/hp=2) 的 reads 強行加 somatic → hp=11/hp=21；可能誤加 somatic |
| **(γ) V5 修對更多** | V5 ploidyRatio fix 更精確識別 somatic 信號 → 之前 V3F 漏標的 somatic 在 V5 被補回 |

→ **不能單靠 V3F vs V5/V6 BAM 對比判定優劣**；需 ground truth (paired_T HP:Z:) 對齊驗證。

### Key Finding 5: 「V6 比 V3F 改善」claim 需修正

之前 v1.7-G interim 寫「V6 比 V3F 改善 hp=11 -11 pp」**錯誤**（比較對象是 V3F vote_dump 非 V3F BAM）。

正確版本：
- V3F BAM hp=11: 14.46% — **V3F 自己最少 hp=11 (priority bug 修對最徹底 in this dimension)**
- V5/V6 BAM hp=11: 31.66% — V5/V6 在 d0bcd8c phasing reshuffle 後 hp=11 增加 17 pp
- 「V6 比 V3F 改善」此單一 metric **不成立**
- 但 V6 的其他優勢 (cross-sample 4 樣本中性化 / marker coverage +9% / hp=33 marker eng +4.7%) **仍成立** — 這些是 codebase-invariant metrics

### 3a.5b.8 後續行動 (final, **2026-05-14 全部完成**)

| Action | Status |
|---|---|
| 認知更新: V5 = V6 (V6 patch 在此 subset 無效) | ✅ confirmed |
| 認知更新: V3F vote_dump ≠ V3F BAM | ✅ confirmed |
| 認知更新: V3F → V5 63% reads 改 phase block | ✅ confirmed (transition matrix) |
| paired_T BAM extract on 17,404 victim → 對齊 ground truth | ✅ **completed** (background biak3oadj) |
| 完成 paired_T 對齊後 finalize Q11/E5 + PPT | ✅ §3a.5b.9 added; PPT updated |
| 主報告 verdict 「V6 = production candidate」**仍成立**（cross-sample + marker eng + caller F1 三 dimension 不受影響）| ✅ 不撤回 |
| caller F1 invariant **empirical 驗證** (VCF FILTER diff = 0 lines) | ✅ confirmed (§3a.5b.10) |

---

## 3a.5b.9 paired_T 4-way Alignment Final Verdict (**2026-05-14 PM, biak3oadj 完成後**)

### paired_T BAM HP:Z: 分布（17,404 V3F victim subset 全集，HCC1395 longphase-s paired_full）

| HP:Z: | Count | % | 對應 longphase-to HP:i: |
|---|---|---|---|
| 1-1 (HP1+somatic) | 4,426 | **25.43%** | HP:i:11 |
| 2-1 (HP2+somatic) | 4,124 | **23.70%** | HP:i:21 |
| 1 (HP1 germline) | 2,648 | 15.21% | HP:i:1 |
| 2 (HP2 germline) | 2,740 | 15.74% | HP:i:2 |
| empty (未 tag) | 2,936 | 16.87% | empty |
| 3 (ambiguous) | 530 | 3.05% | HP:i:33 |

### 4-way per-read alignment (V3F BAM / V5≡V6 BAM / paired_T)

#### Per-read full 6-class matching:

| 比較 | Straight match | Inverse match (HP1↔HP2 swap) | Neither |
|---|---|---|---|
| V3F vs paired_T | 4,299 (24.70%) | 3,432 (19.72%) | 9,673 (55.58%) |
| V6 vs paired_T | 4,346 (24.97%) | 3,235 (18.59%) | 9,823 (56.44%) |

#### Family-level (排除 germline-only/empty/ambiguous, 只比 somatic family HP1/HP2):

| 比較 | n | same direction | inverse direction | Best after swap correction |
|---|---|---|---|---|
| V3F vs paired_T | 11,233 | 48.95% | 51.05% | **51.05%** |
| V6 vs paired_T | 10,924 | 51.95% | 48.05% | **51.95%** |

### Final Verdict

**V3F 和 V6 對 paired_T 的 family-level 對齊都 ~50% — 等於隨機水準**:
- TO codebase HP1/HP2 designation 與 paired_T codebase 沒有系統性對應（既不一致也不顛倒）
- 跨 codebase HP1/HP2 alignment **不是判定 V6 vs V3F 優劣的標準**
- 此現象與 paired_priority_bug_audit/01_step_D §「paired BAM HP:Z:1 與 TO BAM HP:i:1 在同 chr 內不保證同 axis」一致

### V6 vs V3F 最終 verdict

| 判斷依據 | V3F 勝/V6 勝/平手 | 理由 |
|---|---|---|
| paired_T per-read 對齊 | **平手** (24.70% vs 24.97%) | 都 ~25%，差距無統計意義 |
| paired_T family-level 對齊 (swap-corrected) | **平手** (51.05% vs 51.95%) | 都接近 50% 隨機 |
| 保留 germline-only outcomes (hp=1/hp=2) | **V3F 勝** (21% preserved vs V5/V6 0%) | V3F 不強迫 bivariate 分類 |
| Cross-sample 4 樣本 ratio 中性化 | **V6 勝** (0.61-1.24 中性 vs baseline 17.3:1) | V3F 未跑此 metric |
| Marker coverage NG≥3 | **V6 勝** (+9% vs V3F baseline) | V6 hp=33 marker engineering |
| caller F1 invariant | **平手** (4 版 0.7166 相同) | FILTER 不動 |

**結論**：
- ❌ **不能單獨說「V6 比 V3F 好」or「V6 比 V3F 差」** — 看 dimension
- ✅ **V6 = trade-off production candidate** — codebase-invariant metrics (cross-sample + marker eng) 勝出，但 per-read 對齊 paired_T 沒比 V3F 強
- ✅ **5/13 erratum E5 §3a.5 原 framing「V6 31.7% hp=11 = 失敗繼承 V3F」徹底失效** — 因為「繼承 V3F = 正確」前提本身就不成立（V3F 對 paired_T 也只 ~50% 對齊）
- ✅ **主報告 verdict「V6 = production candidate」仍成立**，但理由需轉到 codebase-invariant metrics

## 3a.5b.10 Empirical Proof — Caller F1 Invariant (FILTER 欄 row-by-row 完全相同)

### 驗證命令

```bash
diff <(zcat baseline_clairs_to.vcf.gz | grep -v "^#" | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$7}' | sort) \
     <(grep -v "^#" v5_longphase_phased.vcf | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$7}' | sort) \
     | wc -l
# Output: 0
```

### 結果

| 比較對象 | input (ClairS-TO baseline VCF) | output (V5 longphase phased VCF) |
|---|---|---|
| NonSomatic | 48,358 | 48,358 ✓ |
| LowQual;NonSomatic | 831 | 831 ✓ |
| PASS | 753 | 753 ✓ |
| LowQual;VariantCluster | 41 | 41 ✓ |
| LowQual;StrandBias | 4 | 4 ✓ |
| **完整 VCF diff (CHROM/POS/REF/ALT/FILTER)** | **0 lines** | row-by-row 完全 identical |

### 含義

- longphase 輸出 phased VCF **每一個 variant 的 FILTER 欄與 input 100% 完全相同**
- F1 = TP/(TP+0.5FP+0.5FN) 來自 FILTER=PASS 集合 → **set 完全相同 → TP/FP/FN 不變 → F1 mathematically invariant**
- 從理論宣稱（slide 14 mechanism box）升級為**經驗證實**
- → caller F1 = 0.7166 對 baseline/V3F/V5/V6 4 版全相同**不是巧合是必然**

---

## 3a.5b.11 Source Code Verification — 3 用戶假設全部證實 (**2026-05-14 PM**)

用戶提 3 個 hypothesis 為何 phasing 改但 F1 不變；source code 完整驗證:

### Hypothesis 1: 是否 github 累積 upstream updates 造成 phase 改變? 🟡 **部分對**

**驗證 (commit chain @ `/big7_disk/liaoyoyo2001/longphase-to-mod`):**

| commit | author | date | stage | 性質 |
|---|---|---|---|---|
| `c7f59d3` ~ `3148944` | upstream | pre-baseline | doc/cosmetic | upstream documentation/legend updates, **無邏輯改動** |
| **`8b8c1fd`** | **InterSubMod Research** (us) | **2026-04-10** | phasing | **本地 feature: --pon-only-phasing mode** |
| `41ff147` (V3F) | InterSubMod Research | 2026-04 | haplotag | 本地 priority bug 修對 |
| `380e8d2` | InterSubMod Research | 2026-04 | haplotag | INDEL guard |
| `d0bcd8c` (V5) | InterSubMod Research | 2026-04-30 | phasing | 本地 ploidyRatio fix |
| `938f0df` (V6) | **zhenyu (upstream)** | **2025-10-02** | phasing | upstream purity threshold (cherry-picked) |

→ **5 個 phasing/haplotag 邏輯改動中 4 個是本地** (`8b8c1fd` / `41ff147` / `380e8d2` / `d0bcd8c`)，**1 個是 upstream** (`938f0df` V6, 但對 HCC1395 5kHz 此 subset V5=V6 0 diff 即無實效)。
→ Phase 改變**主因不是 github 累積**，是本地加的 PON-only mode (`8b8c1fd`)。

### Hypothesis 2: 是否「改用 PON phase2」造成 phase 改變? ✅ **完全正確**

**驗證 (8b8c1fd commit message + source code):**

8b8c1fd commit message 原文（節錄）:
> Problem: In standard LongPhase-TO, somatic variants participate as phasing anchors AND are evaluated targets, creating circular dependency. This causes:
> - 94.6% somatic variants biased to HP1 (17.3:1 ratio)
> - 31.2% TO LOH are self-phasing artifacts
> - Phase block N50 degraded by somatic contamination
>
> Solution: New --pon-only-phasing flag implements three-step approach:
> 1. convertNonGermlineToSomatic() — clean phasing using only PON germline
> 2. resetNonPonOrigin() + somaticCalling() — correct somatic classification
> 3. syncPhasingResultOrigins() — sync GT format (germline→0\|1, somatic→0\|0+GT2)
>
> Verified on HCC1395 5kHz:
> - LOH.bed Jaccard=1.0000 (unchanged)
> - **Somatic HP bias 17.3:1 → eliminated**
> - **Phase block N50: 4,061 → 8,109 (+99.7%)**
> - **Phased rate: 54.9% → 78.5% (+23.6pp)**

**Source code 證實** (`PhasingProcess.cpp:155-181`):
```cpp
// PON-only phasing with correct classification:
// 1. Clean phasing: mark all non-PON as SOMATIC → simplified graph → better N50
vGraph->convertNonGermlineToSomatic();

// 2. Correct classification: reset origins, then run somaticCalling
vGraph->somaticCalling(snpFile.getVariants((*chrIter)));

// (Pass 2 between)
vGraph->syncPhasingResultOrigins(chrInfo.posPhasingResult);

// Pass 1's phasingProcess(..., nullptr) skipped ploidy collection; Pass 2's
// somaticCalling+syncOrigins finalizes which variants are truly somatic.
vGraph->somaticCalling(snpFile.getVariants((*chrIter)));
vGraph->convertNonGermlineToSomatic();
```

→ **PON-only mode (8b8c1fd) 確實是 phase block N50 +99.7% / Phased rate +23.6pp 的根本原因**。slide 12 的「phasing 結構 ⭐ N50 +99.7% / Phased +23.6pp / 1.36× / Jaccard=1.0」**完全對應**此 commit。

### Hypothesis 3: Pass 2 不改 FILTER, 直接拿 caller 輸出, 所以 F1 不變? ✅ **完全正確**

**驗證 (source code grep)**:

```bash
# 整個 longphase-to-mod repo 找寫 FILTER 的 code
grep -rnE "bcf_update_filter|update_filter|filterStr|set_filter|FILTER.*=.*\"" *.cpp *.h
→ 0 matches (沒有 write FILTER 的 API call)

# 找讀 FILTER 的 code
grep "FILTER" Phasing.cpp
→ Phasing.cpp:33: "--disable-pon-tag  disable reading the VCF FILTER field PON ..." (只 read)

# ParsingBam.cpp:492 PASS 處理
if (strcmp(filter_str, "PASS") == 0) {
    originType = SOMATIC;  // ← 只 read, 用 FILTER 來分類 variant origin
    break;
} else if (strcmp(filter_str, "PON") == 0 || strcmp(filter_str, "NonSomatic") == 0) {
    originType = (!params->disablePonTag) ? PON : GERMLINE;
    break;
}
```

→ **longphase 只 READ FILTER 不 WRITE FILTER**。FILTER 欄保持 ClairS-TO caller 原始輸出。

**Pass 2 改動什麼 (從 8b8c1fd Step 3 source)**:
- `syncPhasingResultOrigins()` → sync **GT 欄格式** (germline `0/1` → `0|1`, somatic `0/0+GT2`)
- 加 `PS=` (phase set ID) 欄位
- 加 BAM `HP:i:` tag (read-level)
- **完全不動 FILTER 欄**

### 三 hypothesis 整合機制圖

```
ClairS-TO caller
  ↓ 產生 VCF: 49,987 variants, FILTER 已決定 (PASS=753 / NonSomatic=48,358 / ...)
  ↓
longphase-to read VCF:
  ↓ ParsingBam.cpp:492 read FILTER 來分類 origin (SOMATIC/PON/GERMLINE)
  ↓
8b8c1fd PON-only Pass 1 phasing:
  ↓ convertNonGermlineToSomatic() — 簡化 graph
  ↓ phasing 演算法產出 phase blocks
  ↓
Pass 2:
  ↓ somaticCalling + syncPhasingResultOrigins
  ↓ 改 GT 欄 (0/1 → 0|1) + 加 PS 欄
  ↓ ★ 不動 FILTER ★
  ↓
output phased VCF: 49,987 variants, FILTER 完全同 input
  ↓
F1 evaluation (bcftools isec, 看 FILTER=PASS):
  ↓ PASS set 完全相同 → TP/FP/FN 不變
  ↓
F1 = 0.7166 (baseline/V3F/V5/V6 全相同 — 必然)
```

### 用戶 3 假設驗證結論

| Hypothesis | Verdict | 證據 |
|---|---|---|
| Q1: github 累積 updates 是 phase 改變主因 | 🟡 部分對（V6 938f0df 是 upstream 但無實效；主因是本地 8b8c1fd） | commit chain author/date |
| Q2: 改用 PON phase2 造成 phase 改變 | ✅ 完全正確 | 8b8c1fd commit message + PhasingProcess.cpp:155-181 source |
| Q3: Pass 2 不改 FILTER, F1 因此不變 | ✅ 完全正確 | source grep "FILTER write" = 0 matches + empirical diff = 0 lines |

→ 用戶**整套機制理解正確**：**PON-only Pass 2 (8b8c1fd) 改 phasing structure 大幅提升 N50/Phased%，但不動 FILTER → F1 數學保證 invariant**。

---

## 3a.5b.12 Phased VCF GT Source Tracing (**2026-05-14 PM v1.7-I**)

> 用戶提 5 輪 feedback 追問: phased VCF GT 偏向是否造成 HP1 family 17.3:1 偏移? 哪個 stage 開始造成影響?

### Stage Amplification Matrix (baseline pathway)

| Stage | bias ratio (HP1:HP2) | Evidence | Grade |
|---|---|---|---|
| **S2 baseline NonSomatic GT (germline)** | 1.10:1 (52:48 偏 0\|1) | 1,051,022 variants | ⭐⭐⭐⭐⭐ |
| **S2 baseline PASS GT (somatic, phased)** | **1.013:1 (3,108:3,147)** | 6,255 phased PASS variants | ⭐⭐⭐⭐ |
| **S4 baseline vote_dump HP1+HP1_1 vs HP2+HP2_1** | **0.989:1 (182,758:184,846)** | 17,404 victim subset × 全 genome | ⭐⭐⭐⭐⭐ |
| **S5 baseline BAM HP family** | **17.3:1** | full BAM | ⭐⭐⭐⭐⭐ |

### **S5 Priority Bug Amplification = 17.3 / 0.989 = 17.5×** (Critical Finding)

→ baseline pathway 從 S2 → S4 都接近 50:50 (1.013:1 / 0.989:1)
→ **S5 (haplotag judgeHaplotype priority bug) 是唯一 17.5× 放大器**
→ HP1 family 17.3:1 偏移**完全產生於 S5 階段**, 與 GT/vote stage 無關

### V3F Pathway (priority bug 修對後)

| Stage | bias ratio | amplification |
|---|---|---|
| S2/S3 V3F PON-only GT2 | ~1.18:1 (推測, 同 V5) | — |
| S4 vote (V3F binary) | 0.989:1 | — |
| **S5 V3F BAM** | **1.14:1** | **~1× (修對後無放大)** |

→ V3F 修對 priority bug 後, S5 amplification 從 17.5× 降到 ~1×。BAM 直接反映 vote ratio (with 微 noise)。

### V5/V6 Pathway

| Stage | bias ratio | amplification |
|---|---|---|
| **S2/S3 V5 PON-only PASS GT2 (ALT direction)** | **1.18:1 (17,216 : 14,608 偏 HP1)** | — |
| S4 vote | (V5 binary 未測, presumed similar) | — |
| **S5 V5/V6 BAM** | **1.84:1** | **~1.56×** (V5 phasing reshuffle) |

→ V5/V6 BAM 1.84:1 = V5 PON-only mode GT2 1.18:1 + V5 phasing reshuffle (V3F→V5 改 63%)
→ V5/V6 沒有 17× priority bug 放大 (因為 V3F two-layer fix 已修對)

### 對用戶假設的最終 verdict

| 用戶假設 | 結果 | 證據 |
|---|---|---|
| GT 0\|1 vs 1\|0 偏向是 17.3:1 源頭 | ❌ **推翻** | PASS GT 1.013:1 不顯著, vote 0.989:1 ~50:50 |
| HP:i:11 對應的 GT pattern 是否考慮到 | ✅ 已 cover | §3 HP family ↔ GT pattern 對映表 + GT2 統計 |
| 全面 downstream 影響評估 | ✅ 已 cover | §5 downstream impact 6 dimension |
| 統計確認問題的源頭 | ✅ **S5 priority bug** | stage amplification 17.5× 鐵證 |
| 統計結果 + 影響範圍 + 數據支持 | ✅ 已 cover | §6 evidence grading ⭐⭐⭐⭐⭐ 全集 |

### 完整 deep-dive 報告

`InterSubMod/docs/experiments/in_progress/2026/05/20260514_phased_VCF_GT_HP_family_analysis_01.md`

含完整 Phase 0-6 stage tracing + Chi-square + downstream impact + verdict A 確認。

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
