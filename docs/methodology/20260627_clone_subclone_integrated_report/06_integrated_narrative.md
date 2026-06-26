<!--
建立時間: 2026-06-27
類型: L5 整合敘述 + 對抗稽核 verdict（clone/subclone 整合報告）
狀態: in_progress（HCC1395 ⭐3；對抗稽核 NEEDS_WORK → 修正已套用）
-->

# L6（整合敘述）— 用 sSNV + ONT read 建構 clone/subclone：data-supported 綜述

> 敘述框架：Verdict-Pyramid。已過 fresh-context Workflow 對抗稽核（7 agent，verdict **NEEDS_WORK**）→ **8 項修正全數套用**（見 §4）。每 claim 標證據層級。

## §1 結論金字塔（頂層）
**HCC1395 的局部克隆/亞克隆結構，可由 ONT long read 上的 sSNV 單分子共現非循環地重建為 7,143 個局部克隆樹；HP tag 提供「sibling vs allelic」的決定性鑑別、VAF 提供獨立的階層方向驗證；甲基為佐證（待補）。所有推論限 ⭐3 單樣本、局部（非 genome-wide tree）。**

## §2 各層 data-supported 貢獻（修正後）

| 層 | 貢獻 | 關鍵數據（verified） | 證據層 |
|---|---|---|---|
| **骨幹 sSNV 連鎖** | 非循環克隆共現 → 局部樹 | 7,143 區域；full_tree 677；chr17 自動重建 | L1 |
| **L1 HP** | **sibling vs allelic 鑑別器**（非「確認」）| mutual_excl same-HP **0.86× DEPLETED**（baseline 0.500）→ HP 移除 5,238 allelic（57%）| L1 |
| **L2 CCF** | 階層方向**獨立驗證** + 離散層級 | 祖先≥後代：違反僅 **5.7%**（決定性 92.5%/clean 96%）；GMM BIC **n=3** CCF 群 | L1/L2 |
| **L3 PS** | **reliability 旗標** | 92.7% 區域單一 PS（HP 可信）；Tier-PS 不延伸克隆（同 PS CCF 僅 42% 一致）| L1 |
| **L4 甲基** | corroborate（待重抽）| 既有輸出覆蓋 **0.19%**（9/4678）；小樣本 4/8 corroborate；🔴 非 detect | L2 |

## §3 整合敘述（中層）
1. **骨幹**：ONT read 夠長 → 同分子讀多 sSNV → 共現 2×2 → classify（互斥/巢狀/共連）→ 局部克隆樹。這是**唯一非循環**的克隆共現證據。
2. **HP 加值（修正框架）**：same-HP 高在「正共現」關係（nested/co_linked/independent 1.7–1.87×）是**區域背景**（能共讀分析的 sSNV 對本就同單倍型），**非克隆特異證據**；HP 的真價值在**互斥**：互斥 same-HP **低於背景（0.86×）**，因為它混了 sibling（同 HP）與 allelic（異 HP）→ **HP 移除 57% allelic**，把 read 級互斥轉成細胞層 sibling。
3. **CCF 加值**：祖先 VAF ≥ 後代（違反僅 5.7%）= 用 VAF 這條**獨立軸**驗證 read-樹方向；CCF 統計多峰（n=3）= 離散 subclone 層級。
4. **PS**：標記 phase-uncertain 區（排除以提升可信度）；不升級成克隆連鎖（germline phasing ≠ clonal）。
5. **甲基**：corroborate 非 detect；既有輸出幾乎不覆蓋（0.19%）→ 需從 BAM 重抽（後補）。
6. **precedence**：sSNV 連鎖 > HP > 甲基；衝突 genetic 勝。

## §4 對抗稽核 verdict + 修正（誠實留底）
**Workflow（7 agent）verdict = NEEDS_WORK**；以下 8 項修正**全部已套用**：
1. ✅ **L1 chance baseline 0.443→0.500**（只用 gateable HP）→ enrichment 修正、mutual_excl 由「≈chance」改「**DEPLETED 0.86×**」。
2. ✅ **L1 框架**：same-HP 高 = 背景非克隆證據（independent 最高佐證）；HP = 鑑別器非確認器。
3. ✅ **L2 祖先≥後代**：誠實報含 tie（69.8% support / 5.7% violate / 24.5% tie）+ 決定性 92.5%。
4. ✅ **L2 離散性**：補 GMM BIC（n=3）取代「目視峰」；標 CN-mixture 構造警告。
5. ✅ **L3 分母**：含 71 no_ps 區 → 92.7%（非 94%）。
6. ✅ **L4 分母**：覆蓋 0.19%（/4678）非 3%；n=8 小樣本明標。
7. ✅ **00b**：problem PS block **34** 塊（非 37-40）；HP 85-90% 標「跨版本 concordance 非本資料量測」；「96%」= V5-改變 sites 非誤差集中。
8. ✅ **L0**：somatic=None 30% 分母透明化；hp=克隆改「必要非充分」；「linked」操作定義明述。

**對抗稽核認定 solid（未被推翻）**：HP code/BERNOULLI 公式源碼逐字驗證；HP 對 mutual_excl 的 allelic 移除（57%）為真；VAF 梯度違反率極低（5.7%）；甲基 corroborate-not-detect + cis-ASM double-dip 框架源碼驗證；單樣本紀律無跨樣本 overclaim。

## §5 🔴 限制（無可被質疑的未佐證推論）
⭐3 單樣本；**regional（≤read-span）非 genome-wide tree**；HP ~85-90% phasing 誤差 + 34 problem PS block（已標排除）；CCF 僅 CN-clean 可估（gain multiplicity 歧義）+ CN-mixture 構造；甲基 corroborate 非 detect + 既有覆蓋 0.19%（待重抽）；64% sSNV 在 CN-gain 混淆（乾淨集 = LOH/neutral）；precedence「genetic 勝」僅 1 個實證衝突 + n=8 甲基測試（未大規模壓測）；分子證據非 single-cell confirmation。**對外勿稱「甲基偵測 subclone / genome-wide tree / 對手缺檢定」。**
