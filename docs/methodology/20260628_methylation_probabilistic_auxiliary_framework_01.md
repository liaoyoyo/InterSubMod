<!--
建立時間: 2026-06-28
類型: methodology — 甲基作為機率輔助層的猜想裁決 + 完整 3-tier 證據框架
狀態: in_progress(框架設計 sound;甲基層 validation 待 T-GATE-GB cis-control)
build_branch: research/subclonal-reconstruction-202606
provenance: chr17 per-read 實測(cis-dominated) + external_validation axis4/5(Sgootr/Gaiti/CAMDAC) + 凍結資料
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data
-->

# 甲基機率輔助層 — 猜想裁決 + 完整 3-tier 拓樸證據框架

> 框架：scientific-rigor §2(每裁決標 L1-L5)+ §4(DAG)。重述猜想 → 文獻+資料查證 → 方法學裁決 → 完整框架 → 條件與邊界。

## §0 猜想（用戶，重述確認）

在資訊缺乏（pairwise 缺口、欠定）的區，甲基可否當**機率輔助**：
1. 用「與 normal **同 HP** 的甲基差異距離」判哪個分群較可能是**子突變**（較發散=較晚/衍生）。
2. 用「甲基相似性 read 分群」當 **sSNV ALT 之外的第二軸**，輔助驗證串接 + read 距離，以**機率方式**定可能拓樸。
3. 證據層級：**直接突變共現（最可信）> HP 分哪棵樹 > 甲基機率提升**。

## §1 查證（文獻 + 本資料）

| 查證點 | 結果 | tier |
|---|---|---|
| 甲基能否重建 lineage（前研究）| ✅ **Sgootr / Gaiti：單細胞甲基可重建 lineage**（external_validation axis4）→ 概念有先例，但**單細胞**；bulk read-level 未建立 | L1(文獻) |
| normal-baseline 處理 cis 是否正確 | ✅ **CAMDAC（TRACERx 2025）用 matched-normal 做 cancer ASM**（axis3）→ normal-same-HP 基線是**已建立的正確 cis 處理**（= 本專案 T-GATE-GB）| L1(文獻) |
| 甲基「發散=較晚」是否有分子鐘依據 | 🔴 **無已驗證的 bulk 甲基鐘**；只有 fluctuating-CpG **單細胞**鐘（Gabbutt 等，待搜尋）→ 「離 normal 越遠=越衍生」是**啟發式非定律** | L3 |
| 甲基 read 分群是否獨立於 genotype（本資料）| 🔴 **chr17 實測：分群對齊 cis-genotype 軸（α 23 ≫ lineage 6）**、cis-control 0/740 未跑 → 甲基距離**被突變 cis 效應主導**（double-dip）| L1(本資料) |

## §2 方法學裁決

**框架（3-tier 分層證據）= ✅ 方法學上適當**（這是標準的分層/貝氏式證據整合，排序正確）：
- Tier 1 直接 sSNV 共現 = **likelihood，決定性、非循環**（單分子物理連鎖）。
- Tier 2 HP 定根 = **partition，決定性**（germline phasing 分哪棵樹）。
- Tier 3 甲基 = **prior，機率性、非決定性**（只在 Tier 1/2 欠定時微調候選機率）。

**甲基 Tier-3 層 = 🟡 設計合理但目前未驗證（validation-blocked）**：
- ✅ **正確處**：(a) 排序正確（甲基永遠最低、不override genetic）；(b) **用 normal-same-HP 基線是正確的 cis 處理**（CAMDAC 先例）；(c) 機率而非硬判（誠實表達不確定）。
- 🔴 **未驗證處**：(a) 「離 normal 越遠=越衍生」**無 bulk 甲基鐘依據**（啟發式）；(b) 本資料甲基距離被 cis-genotype 主導（chr17）→ **未扣 cis 前，甲基的「機率提升」是循環的**；(c) 跨分子比甲基（缺口的兩 read 不同分子）更混淆。

**→ 裁決：框架適當、可寫進方法當「設計的未來機率層」；但甲基層的可信度 BLOCKED 在 T-GATE-GB（matched-normal cis-control）——cis-control 跑通且證實 normal-baseline 甲基差異獨立於突變 cis 之後，Tier-3 才從「設計」升「可用」。** 在此之前甲基只標註、不參與機率裁決。

## §3 完整框架（形式化）

對每個區的拓樸 T，後驗機率：

```
P(T | data) ∝ P(sSNV共現 | T)              [Tier 1 likelihood，決定性]
              × 1[HP-root partition 相容]    [Tier 2，硬約束:跨HP拆獨立樹]
              × P(甲基 | T, normal-baseline)^w  [Tier 3 prior，w↓、僅欠定區、cis-control-gated]
```

- **Tier 1**：只認同分子共現（AA 格定 vertical/horizontal、零格定方向、ε=2%）。determined 區（~11%）由此唯一定。
- **Tier 2**：HP tag 分根（1850 兩根/3697 一根）；跨-HP 強制拆獨立樹（非機率，硬約束）。
- **Tier 3**：僅當 Tier 1 留多棵相容樹（欠定 78.5%）時啟用：
  - 用 **normal-same-HP 甲基距離** 當基線（扣 germline-ASM cis）；
  - read×read 甲基相似度當第二軸，給每棵候選樹一個 weight；
  - **w 必須 ↓（弱權）+ 標 cis-control-pending + 永不 override Tier 1/2**。
  - 輸出 = 候選樹**機率排序**（非單一確定樹）。

## §4 條件與紅線（對外必守）

1. 🔴 **甲基永遠 Tier 3、永不 override genetic/HP**；genetic 已 determined 處**不用**甲基（冗餘且 cis 循環）。
2. 🔴 **必須 normal-same-HP 基線**（扣 cis-ASM）；未扣 cis 的甲基機率**不可用**（chr17 證 cis 主導）。
3. 🔴 **「發散=衍生」是啟發式 prior 非定律**（無 bulk 甲基鐘）；只給機率傾向、不給確定方向。
4. 🔴 **validation gate = T-GATE-GB**：cis-control 跑通 + 證 normal-baseline 甲基獨立於突變 cis，Tier-3 才升「可用」。
5. 🔴 **對外用詞**：單細胞甲基-lineage 已有（Sgootr/Gaiti）→ 本方法定位 = 「bulk long-read 上 sSNV-骨幹為主、甲基為**有界機率輔助**」，**不宣稱甲基重建 lineage**。
6. 🔴 全程 ⭐3 單樣本上限；甲基層升級需 ≥5/7 樣本 + COLO829。

## §5 與已 commit 文件的關係

- Tier 1/2 的判準 + 標籤 → `20260628_lineage_label_definition_01.md`（已 commit 8339605）
- ε=2% 門檻 → `20260628_sSNV_linkage_threshold_decision_eps2_01.md`
- 甲基 bounded-auxiliary 證據（chr17 cis）→ `20260628_reconstruction_model_verification_01.md` §甲基
- 本檔 = 把甲基從「bounded-auxiliary（descriptive）」進一步形式化為「**Tier-3 機率輔助層（設計 sound、validation 待 cis-control）**」。

**一句話**：你的猜想**框架方法學上適當**（3-tier 排序正確、normal-baseline 是正確 cis 處理、機率而非硬判），但甲基 Tier-3 層**目前未驗證、blocked 在 T-GATE-GB cis-control**；在 cis-control 證實 normal-baseline 甲基獨立於突變 cis 之前，甲基只標註不參與機率裁決。
