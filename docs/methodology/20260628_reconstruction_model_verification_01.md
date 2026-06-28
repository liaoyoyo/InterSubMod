<!--
建立時間: 2026-06-28
類型: methodology — 用戶 sSNV 重建證據階層 + 6 主張 + 拓樸/分群/重建邏輯的大規模對抗驗證裁決
狀態: in_progress（15-agent workflow wf_f2b070ea-64c 驗證收斂；待用戶逐項確認）
build_branch: research/subclonal-reconstruction-202606
provenance: workflow wf_f2b070ea-64c（15 agent：5 碼 + 6 階層/主張 + 3 對抗 + 1 綜合）；原始碼/資料源 branch feat/summary-nreadsvalid@5308d9e（pending-merge）凍結於 _assets/20260627_subclone_4axis_teaching/{data,scripts}
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data, docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts
-->

# 重建模型大規模驗證裁決（拓樸 / 分群 / 重建 / 整合驗證）

> 框架：Verdict-Pyramid + scientific-rigor §2（每裁決標 L1-L5 + verdict）。重點＝確認**計算 sSNV 的實作細節**是否符合資料/生物/資訊定義 + 逐項裁決用戶模型。

## §0 headline

用戶提出的「同-read 共現骨幹 → 跨-read 交集 → 巢狀深度 → VAF → HP/LOH → CNV → 甲基」階層**方向大致正確、骨幹（單分子共現）是對的非循環脊柱**，但驗證出 **3 處結構性修正 + 1 主張 INVALID + 2 實作碼 bug**：
1. **CNV/LOH 不是證據 rung，是 confound precondition**（先扣除才能讀 VAF；影響 53–69% 骨幹）。
2. **VAF 是 consistency-check 非獨立證據**；**HP 是鑑別器（只在互斥有診斷力）非確認器**。
3. **主張 3（用甲基判突變先後時序）= INVALID**（同基因型 read 零 ordering 訊號 + 甲基非分子鐘 + double-dip）。
4. 實作碼 2 個 CODE_DISCREPANCY：**somatic 定義不一致（==0 vs <5%）** + **classify() 噪聲容忍不對稱（AA 容忍 1、off-diagonal 嚴格 ==0）**。

## §1 階層重構（原 7-rung → 修正結構）

| 原 rung | verdict | 修正 |
|---|---|---|
| 1 同-read 多點 sSNV | ✅ SOUND | 唯一非循環骨幹，但 read=cell 在 chimera(ONT 1-3%)/segdup-paralog/CN-gain 失效；骨幹 **68.9% 落 CN-gain** → 非循環性在 ≥53% 骨幹 **UNDETERMINED**；限 CN-clean∩non-segdup 才宣稱方向 |
| 2 跨-read 交集 | 🟡 NEEDS_CARE | 順序對、弱於 #1；cross-PS(>50kb)=GAP；within-read 對 phase-switch 免疫是關鍵 asymmetry |
| 3 巢狀深度 | 🟡 NEEDS_CARE | **非獨立證據**，是 #1-2 的**重建輸出 topology**；計數 L1 / 生物樹詮釋 L2 |
| 4 VAF tiebreaker | 🔴 RISKY | **consistency-check 非獨立**；必先 condition on CNV、只 CN-clean 可估；永遠連報 69.8/5.7/24.5（禁只報 5.7%）|
| 5 HP（LOH）| 🔴 RISKY | **HP≠LOH 必拆**：HP=鑑別器（只在互斥 0.86× DEPLETED）；**LOH 不是 HP 子類，是 confound**（HP 塌陷 + imprinting-unmask 82-91%）|
| 6 CNV | 🔴 INVALID（當證據）| **改列 condition-on confound（precondition）**；standalone subclone 證據已 REFUTED；僅 SV-anchor 才有 L3 |
| 7 甲基 | 🟡 NEEDS_CARE | 位置對但應 **off-ladder bounded-auxiliary**；6.6% 弱 corroborator；cis-control 0/740 → subclone-specificity UNDETERMINED |

**修正後結構（非平鋪 7-rung）**：
- **A 非循環骨幹（genetic 共現）**：A1[L1] 同-read 多點 sSNV（最強，限 read-span/CN-clean∩non-segdup）> A2[L1] 跨-read 交集（cross-PS GAP）> A3[L1計數/L2詮釋] nesting topology（是輸出非獨立證據）
- **B 條件式一致性檢查（非獨立）**：VAF magnitude[L2]（先 condition on CNV、只 CN-clean、連報三數）
- **C 鑑別/rule-out（非確認）**：HP tag[L1]（只在互斥有診斷力）
- **precondition（confound，非證據）**：CNV[L1 confound] + LOH[confound]（與 CNV 同類）
- **D off-ladder 佐證（非偵測）**：甲基[bounded-auxiliary；subclone-signal UNDETERMINED]
- precedence「衝突時 genetic 勝甲基」= **設計立場非已驗階層**（壓測僅 n=1 衝突 + n=8 甲基）

## §2 6 主張逐項裁決

| 主張 | verdict | 核心理由 | tier |
|---|---|---|---|
| **1 二代含一代（nesting=parent-containment）** | ✅ SOUND | = infinite-sites→perfect-phylogeny（所有 bulk 重建地基）。⚠ 修詞：het→hom 主機制是 **LOH 不是第二次點突變** → AA 應稱「二代**事件**」；記號 1-1/1-1-1 是 lineage 代數，與 HP haplotag 1/2 **正交須標軸別** | L1/L3 |
| **2 互斥解析（diff-HP=allelic / same-HP=sibling）** | 🟡 NEEDS_CARE | 鑑別方向對。但 (a) same-HP 本身不是 sibling 確認（85-94% 是背景）；(b) 🔴 **把 LOH 當「乾淨定義 subclone」方向相反**——LOH 是假陽最高發處（corroborated 41/49 正落 LOH〔來源 06_narrative §3，非本 bundle 凍結 JSON〕）；(c) 🔴 **用甲基距離定 1-1 = cis-ASM/double-dip**（甲基 recover partition=0/740）| L1-L2/L3 |
| **3 甲基判突變先後時序（RRR/AAA）** | 🔴 **INVALID** | (i) 同基因型 read **零 ordering 訊號**（先後需 nested 基因型）；(ii) 甲基**非分子鐘**無時序映射；(iii) double-dip。a fortiori：有遺傳訊號的 740 區甲基都 recover 0 | L1/L5 |
| **4 更多位點→更多群→更完整演化** | 🟡 NEEDS_CARE | 「更多群」power 內成立；「更完整」被 read-span ≤50kb **封頂為 regional**；**中位 sSNV 間距 51.4kb〔來源 20260622 spacing doc，非本 bundle 凍結〕、P(單 read 跨≥2)≈0.74%、94.2% read 一個都碰不到**；計算量非瓶頸、**sSNV 稀疏才是** | L1/L2/L3 |
| **5 甲基定義 clone→外推單位點** | 🔴 RISKY | 前提「甲基能定義 subclone」**未成立而非已成立**（recover 0、cis-control 0/740 未測）；+ 層級轉移謬誤。標 RISKY 非 INVALID 因前提 UNDETERMINED（待 T-GATE-GB），是「高風險待驗證未來方向」 | L1/pending |
| **6 還有哪些衝突/算法細節** | 🔴 有實質須修 | chr17 headline 硬編 4-sSNV vs 資料 3-sSNV；somatic 雙定義；classify aa<2 vs 文件宣稱 AA=0；孤兒 chr17_*.json；γ FP-source 被畫成確定 subclone | L1/L2 |

## §3 實作碼計算細節（你最想查的）

| # | 項目 | verdict | 發現 |
|---|---|---|---|
| A1 | classify() 2×2→config | CODE_OK（邏輯）+ ⚠ | 方向判定正確；但見 A1-bug |
| **A1-bug** | **噪聲容忍不對稱** | 🔴 CODE_DISCREPANCY | `aa<2`（容忍 1 個 AA）vs off-diagonal `==0`（嚴格）→ ONT 數% 錯誤率下單一雜訊/嵌合 read 把真 co_linked 降 nested、把 nested 降 independent → **系統性低估 lineage 桶**；反向單一假 AA 即推翻 mutual_excl。建議改 Fisher/binomial over-dispersion（呼應專案「Fisher over-dispersion 必修」）|
| **A2** | **somatic 定義一致性** | 🔴 CODE_DISCREPANCY | 全基因組 `is_somatic` = nref≥4 且 **normal VAF<5%**；chr17 per-read 檔用嚴格 **==0** → 同一 β1（normal 3.4%）一處 True 一處 false；**全基因組差 2,230 sSNV（9.05%）**；嚴格==0 在高覆蓋 normal **系統性 false-negative**。建議統一 + depth-adaptive 統計檢定 |
| A3 | 方向指派 | ✅ CODE_OK | **方向純由零格拓樸決定，完全不碰 VAF/read 數**（證實非循環）；chr17 VAF 同向是巧合佐證非依據。脆弱性：coread≥6 即判零格（absence≠evidence）|
| A4 | component/成環 | ✅ CODE_OK | Kahn 拓樸排序偵測；23 環/6146(0.37%)保留+標 inconsistent+排除 full_tree（**正確保守**）；19/23 落 CN-gain、10/23 為 2-node 純反對稱 |
| A5 | HP-gate 非循環 | ✅ SOUND_WITH_CAVEATS | HP tag 來自 BAM 既有（longphase-S，**非從受測 sSNV 衍生→非循環**）。⚠ HP3('3') 未排除 → ('3','3')誤判 same_hp → 污染 sibling_sameHP **6.5%**；same_hp 無一致性門檻；JSON 無 n_hp3 不可稽核 |

## §4 成環能否用 HP/甲基 resolve？（你直接問的）

**答案：不能。** 有碼根據：
1. **HP 已用盡**——建圖只用 same_hp 邊（`sm_region_integration.py:152` filter `p['same_hp']`），環是在 same-HP 條件下「仍」發生 → HP 無剩餘鑑別力；且 HP 二元、與祖先/後代方向正交、無排序力。
2. **甲基非分子鐘**——專案紅線定為 characterize 非偵測、cis-control 0/740、無時序映射。
3. **VAF 也救不了**——8/10 小環 |ΔVAF|<0.05（tie）、19/23 落 CN-gain（multiplicity 扭曲 VAF）。

環的數學意義 = **nested 偏序自相矛盾（four-gamete/ISA 違反）= 「資料不支持單一克隆樹」的誠實訊號**；強解成樹 = 捏造 lineage order。現行（保留+標 inconsistent+排除 full_tree）正確。真正的「解」= 補 mappability/segdup mask + condition-on CN 把 19/23 CN-gain 環分離為候選 artifact；10 個 2-node 小環可考慮 collapse 成 co_linked 單節點。

## §5 資料輸出完整性（你要確認的「有清楚定義與輸出供日後驗證」）

**✅ 可驗（grep-able）**：cn_gain 52.8%、HP 0.86×/移除 57%、甲基 recover 0/740·6.6%、CCF 69.8/5.7/24.5·GMM n=3、phaseset 7143·full_tree 677、linked CN-gain 68.9%、Fisher≈permutation 97.7%、by_shape 加總、sSNV 間距 51.4kb〔來源 20260622 spacing doc，非本 bundle 凍結〕。

**🔴 缺口（待補才完整可驗）**：
1. `no_confirmed_structure=2443` + `inconsistent=22` 無 JSON key（只能 7143−4678−22 反推）
2. CLEAN full_tree=205/763/408/357 **硬編** build script，無 JSON，grep '205' 三處錯誤碰撞
3. chr17 headline（4 sSNV、7/10/15）硬編 SVG **與孤兒檔 chr17_subclone_data.json（3 sSNV、9/20/19）矛盾**
4. `sm_region_integration.json` + `regions.tsv` **只在 gitignored worktree** → 7143 區無法 drill 到個別 shape/members/CN/edges/reads；FF-merge 或清 worktree 後**主張 4 證據全失**
5. 無 CN-clean per-region flag、無 6-bucket shape_distribution、無 segdup/mappability mask 輸出
6. 無 HP-coverage/n_hp3 → same_hp 桶污染（6.5%）不可稽核

## §6 待用戶拍板的決策（排序，最 load-bearing 在前）

1. **[最重要] CNV/LOH 改列 condition-on confound（非證據 rung）？** 影響 53-69% 骨幹；連帶 VAF 降條件式、LOH 從 HP 拆出。
2. **somatic 定義統一為 <5%（+nref≥4）還是 depth-adaptive 檢定？** 決定 chr17 canonical 是 2/3/4 sSNV、全基因組 2,230 sSNV 口徑——是 teaching dossier 內部矛盾的根。
3. **classify() off-diagonal 改 Fisher/binomial over-dispersion（對稱容忍）？** 先跑 ε=1/coread sensitivity 量化變動。
4. **甲基定案** off-ladder BOUNDED_AUXILIARY + 封殺主張 3 + 凍結主張 5（待 T-GATE-GB）+ headline 用 6.6% 非 50%。
5. **§13-A 合規修補**：chr17 SVG 改 data-injected、CLEAN 205 落 clean_subset.json、刪/回填孤兒檔。
6. **output completeness**：sm_region_integration.json aggregate + per_region.tsv 落版控 + 補 no_confirmed_structure/inconsistent/6-bucket JSON key。
7. **對外 novelty**：丟「first read-level methyl lineage」（Foltz 2024 已有 single-molecule co-occurrence primitive）；novelty 收窄為 native ONT + read×read PERMANOVA + normal-anchored cis-test + 有界甲基；單樣本封頂 ⭐3。
8. **HP gate 收尾**：same_hp 排除 '3'/None（sibling_sameHP −6.5%）+ 加 HP-coverage 欄位 + phasing provenance stamp。
