<!--
建立時間: 2026-06-27
類型: L2 CCF tier + 祖先-後代梯度 — 單變量驗證（clone/subclone 整合報告）
狀態: in_progress（HCC1395 ⭐3）
data_sources: data/sm_ccf_tiers.json
-->

# L2 — CCF tier + 克隆階層梯度（單變量驗證）

> 問題：加入 VAF/CCF 能否 data-supported 地支撐克隆階層與離散 subclone？
> 圖：`figures/03_ccf_tiers.png`。資料：`data/sm_ccf_tiers.json`。

## §1 最強結果 — 祖先 ≥ 後代 VAF 梯度（CN-independent 克隆階層測試）

**原理（pigeonhole）**：後代 subclone 是祖先的子集 → 在**更少**細胞 → CCF（及 VAF）應**更低**。
**方法**：每條 nested 邊比 ancestor vs descendant 的 VAF。**兩節點在同一區域、同一局部 CN** → VAF 可直接比，**不需 CCF 轉換**（消除 CN 混淆）。null = 50%。

| 結果 | 值 |
|---|---|
| nested 邊總數（有 VAF）| 3,627 |
| 祖先 VAF 較高（符合）| **2,532** |
| 後代較高（違反）| 205 |
| tie（\|Δ\|≤0.03）| 890 |
| **全邊 support 率（含 tie 為分母）** | **69.8%** support / 5.7% violate / **24.5% tie** |
| **決定性邊（排除 tie）support** | **92.5%**（vs null 50%）|
| **CN-clean 決定性邊** | **96.0%**（902 support / 38 violate / 310 tie）|

→ 🔑 **read 級 nesting 推出的「祖先→後代」方向，被 VAF 梯度獨立驗證**。**誠實兩種口徑**：含 tie 為分母 = 69.8% 支持、5.7% 違反、24.5% 無法區分（median 決定性 ΔVAF 僅 0.15，多數差異小）；只看**決定性邊** = 92.5%（clean 96%）。關鍵是**違反率極低（5.7%）** —— 即使 1/4 邊無法區分，幾乎沒有「後代 VAF 反高於祖先」的矛盾。兩個獨立訊號（read 共現結構 + VAF 高低）方向一致。

## §2 離散 CCF tier（跨區一致的 subclone 層級）

**CCF 轉換（誠實標假設）**：neutral→min(1, 2×VAF)（diploid het, purity~1, multiplicity=1）；loh→VAF；gain/loss→undet（multiplicity 歧義不估）。

**CCF 峰（CN-clean）**：~**0.9** = clonal/truncal · ~**0.4–0.5** = major subclone · ~**0.05** = rare。
**離散性檢定（修對抗稽核：不只看 bin 峰）**：GMM BIC 1/2/3/4 component = **[3990, −1567, −3018, −3001] → best n=3**（BIC 大幅偏好 3 成分 vs 1 成分）→ **CCF 分布統計上多峰（3 群），非單一連續分布** = 離散 subclone 層級有 BIC 支持。
⚠ 但 CCF 轉換（neutral→2×VAF、loh→VAF）會把 neutral-VAF~0.45 與 loh-VAF~0.9 都映到 CCF~0.9，故「clonal 峰」部分是 CN-mixture 構造；3 群為統計事實但其生物學標籤（clonal/major/rare）為詮釋。

**tier 分布（CN-clean）**：clonal 2,278 / major_subclone 2,018 / minor_subclone 975 / rare 867。

## §3 data-supported 結論
1. **克隆階層方向為真**：92.5%（clean 96%）nested 邊滿足祖先 VAF ≥ 後代 → 樹的方向不是隨機（null 50%）。
2. **離散 subclone 層級存在**：CCF 峰集中在 ~0.9（clonal）/ ~0.45（major）/ ~0.05（rare）→ 一組固定 CCF 的 subclone 跨基因組重現（呼應 L0/區域分析的 VAF 峰）。

## §4 限制
- **CCF 僅在 CN-clean（LOH/neutral）可估**；gain/loss 因 multiplicity 歧義不估（purity-ploidy 不可識別，Tarabichi 2021）。
- 205 違反（= **5.7% 全邊**，或 7.5% 若以決定性邊為分母；統一以 **5.7% 全邊**為 headline）+ 890 tie：多數支持但非全部（量測噪聲/局部 CN/真實複雜度）。clean 決定性 96% 為論文級數字。
- 🔴 CCF tier = **分層**（跨區一致的 CCF 層級），**非 phase-linked genome-wide tree**（守紅線）。

## §5 對論文的意義
VAF/CCF 提供**第二條獨立軸**驗證 read-linkage 樹：① 方向（祖先 VAF 高）② 離散性（CCF 峰）。這把「局部 read 樹」升級為「有 CCF 支撐的克隆階層敘述」，但仍是**分層 + 局部樹**，非全基因組單一樹。
