<!--
建立時間: 2026-06-27
類型: L5 整合敘述 + 對抗稽核 verdict（clone/subclone 整合報告）
狀態: in_progress（HCC1395 ⭐3；對抗稽核 NEEDS_WORK → 修正已套用）
-->

# L6（整合敘述）— 用 sSNV + ONT read 建構 clone/subclone：data-supported 綜述

> 敘述框架：Verdict-Pyramid。已過 fresh-context Workflow 對抗稽核（7 agent，verdict **NEEDS_WORK**）→ **8 項修正全數套用**（見 §4）。每 claim 標證據層級。

## §1 結論金字塔（頂層）
**HCC1395 的局部克隆/亞克隆結構，可由 ONT long read 上的 sSNV 單分子共現非循環地重建：35,332 sSNV → 7,143 個局部區域（其中 677 完整分支樹 full_tree / 4,678 有確認結構 / 其餘為 sparse/single-pair）；HP tag 提供「sibling vs allelic」的決定性鑑別、VAF 提供獨立的階層方向驗證；甲基為弱佐證（6.6% 區域）。所有推論限 ⭐3 單樣本、局部（非 genome-wide tree）。**

## §2 各層 data-supported 貢獻（修正後）

| 層 | 貢獻 | 關鍵數據（verified） | 證據層 |
|---|---|---|---|
| **骨幹 sSNV 連鎖** | 非循環克隆共現 → 局部樹 | 7,143 區域；full_tree 677；chr17 自動重建 | L1 |
| **L1 HP** | **sibling vs allelic 鑑別器**（非「確認」）| mutual_excl same-HP **0.86× DEPLETED**（baseline 0.500）→ HP 移除 5,238 allelic（57%）| L1 |
| **L2 CCF** | 階層方向**獨立驗證** + 離散層級 | 祖先≥後代：違反僅 **5.7%**（決定性 92.5%/clean 96%）；GMM BIC **n=3** CCF 群 | L1/L2 |
| **L3 PS** | **reliability 旗標** | 92.7% 區域單一 PS（HP 可信）；Tier-PS 不延伸克隆（同 PS CCF 僅 42% 一致）| L1 |
| **L4 甲基** | **有界弱** corroborate（已重抽 + 有效性審查）| 重抽 740 區：**49（6.6%）corroborate**；🔑 power-gated（高功率區 54.9%，僅 51 區達）；🔴 cis-control **14/49 可評估、11 germline-cis**（str/int bug 修正）→ 甲基多為 allele-specific **非** subclone；0 新 partition；非 detect（見 `08`）| L1/L2 |

## §3 整合敘述（中層）
1. **骨幹**：ONT read 夠長 → 同分子讀多 sSNV → 共現 2×2 → classify（互斥/巢狀/共連）→ 局部克隆樹。這是**唯一非循環**的克隆共現證據。
2. **HP 加值（修正框架）**：same-HP 高在「正共現」關係（nested/co_linked/independent 1.7–1.87×）是**區域背景**（能共讀分析的 sSNV 對本就同單倍型），**非克隆特異證據**；HP 的真價值在**互斥**：互斥 same-HP **低於背景（0.86×）**，因為它混了 sibling（同 HP）與 allelic（異 HP）→ **HP 移除 57% allelic**，把 read 級互斥轉成細胞層 sibling。
3. **CCF 加值**：祖先 VAF ≥ 後代（違反僅 5.7%）= 用 VAF 這條**獨立軸**驗證 read-樹方向；CCF 統計多峰（n=3）= 離散 subclone 層級。
4. **PS**：標記 phase-uncertain 區（排除以提升可信度）；不升級成克隆連鎖（germline phasing ≠ clonal）。
5. **甲基（已重抽 genome-wide + L8 有效性審查）**：corroborate 非 detect；從 BAM 重抽 740 區，僅 **6.6%（49）corroborate**；🔑 訊號**存在但 power-gated**（popB_n≥20 達 54.9%，僅 51 區達）；🔴 **cis-control 14/49 可評估、11 germline-cis（甲基非 subclone）、3 候選**（2026-06-29 str/int bug 修正；先前「0/740 不可評估」為型別 bug 假象）；CN 分層（使用者校正）：neutral 8/8 可評估・7 germline-cis；LOH 35/41 單倍型不可評估；0 新 partition → **有界弱 corroborator（where testable 多為 allele-specific 非 subclone）**，詳見 `08_methylation_sufficiency_audit.md`。
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

**對抗稽核 round 2（re-audit）= NEEDS_WORK → 再修**（root cause：修正只在 prose、JSON/HTML 未同步 = §13-A 違規）：
- ✅ **F-B**：重算 JSON（`sm_phaseset_extension.json` reliable_rate→0.927 含 71 no_ps；`sm_ccf_tiers.json` 加 `gmm_bic_1to4 [3990,−1567,−3018,−3001] best_n=3`）→ HTML 重新注入正確值（非手打）。
- ✅ **F-A**：「7,143 個局部克隆樹」改「7,143 區域（677 full_tree / 4,678 結構 / 其餘 sparse）」。
- ✅ **#5/#6**：04 殘留「94%」→92.7%；05「~3%」→0.19%（既有輸出）/ 6.6%（重抽結果）。
- ✅ **A1 甲基重抽**：n=8 anecdotal → genome-wide 740 區、6.6% corroborate、驗證 corr=1.000。

**對抗稽核 round 3（final confirm）= NEEDS_WORK → 修**（fix #2 框架層未同步 = 與 round-2 同根因，只是發生在 HP 框架而非數字）：
- ✅ **F1**：HTML L1 box + `sm_hp_contribution.json` verdict 仍留 round-1 已駁回的「same-HP 高 = 真同單倍型克隆連鎖」→ 改正為「same-HP 高（含 independent 最高）= 區域背景屬性、非克隆特異；HP 診斷力在 mutual_excl DEPLETED = sibling vs allelic 鑑別器」；HTML 從修正後重生。
- ✅ **F2**：HTML 骨幹 box `full_tree —`（未注入的 placeholder，因 generator 引用不存在的 `sm_region_stats.json`）→ 改從已載入 `sm_phaseset_extension.json` 的 `by_shape_reliable_vs_uncertain` 衍生 full_tree **677**（596+81）+ structured **4,678**。
- ✅ **F3**：partition 用詞改「4,678 有確認結構（含 677 full_tree）/ 其餘 2,465 sparse」杜絕 677+4,678 雙重計數誤讀。
- ✅ **F4**：L4「0% cis」headline 加「HP-axis proxy cis-control，完整 normal-baseline NACT 為精煉項」限定詞。

**三輪共通根因（已內化）**：修正若只落在 prose 而 JSON/HTML 未同步 = §13-A 違規；每輪修正必 prose+JSON+HTML 三處一致（grep 三方驗證）。

**對抗稽核 round 4（final confirm）= NEEDS_WORK → 修**（F1-F4 全 confirmed CLOSED；新發現非層同步類，而是 own-data 矛盾）：
- ✅ **限制段「64% sSNV 在 CN-gain」無法從本報告 data/ grep**（為骨幹 commit 0a8658d 不同 CN 分段的 stale 值，prose+HTML 彼此一致但都與 master TSV 矛盾 ~11pp）→ 從 `data/sm_locus_master.tsv` 重算（獨立 verify）並落 `sm_locus_master_summary.json` cn_somatic_pct：**somatic gain 52.8%（12,569/23,810）/ CN-confounded 53.2% / 乾淨集 LOH+neutral 46.8%**；HTML 注入 + 06 §5 改正。⚠ 04 §3 的 64% 是另一個合法值（Tier-PS same-HP 80,980/127,183），不動。
- ✅ **（minor）L2 HTML headline 0.925** 加「決定性邊；含 tie 全邊 69.8% 支持/違反 5.7%」對齊 prose 誠實口徑。

**對抗稽核 round 5（收斂確認）= PASS** ✅：64%→52.8% 修正三層閉合 + 獨立從 TSV 重算 verify（12,569/23,810=52.8%）；全報告 own-data-contradiction 掃描乾淨；round 1-4 修正無回歸。唯一 surface 的 prose-only 數字「clean 310 tie」（`03` §1）已補進 `sm_ccf_tiers.json` clean_only.tie（重算 verify=310，902+38+310=1,250 clean edges）→ **每個數字皆可從 data/ grep**。**5 輪 fresh-context 稽核收斂，報告 finalize-ready。**

**對抗稽核認定 solid（兩輪未被推翻）**：HP code/BERNOULLI 公式源碼逐字驗證；HP 對 mutual_excl 的 allelic 移除（57%）為真（5238/9187）；VAF 梯度違反率極低（5.7%）；chr17 canonical 對 TSV 逐行驗證；35,332=30,490+4,842；甲基 corroborate-not-detect 兩輪一致；單樣本紀律無跨樣本/genome-wide-tree overclaim。

## §5 🔴 限制（無可被質疑的未佐證推論）
⭐3 單樣本；**regional（≤read-span）非 genome-wide tree**；HP ~85-90% phasing 誤差 + 34 problem PS block（已標排除）；CCF 僅 CN-clean 可估（gain multiplicity 歧義）+ CN-mixture 構造；甲基 corroborate 非 detect + 既有覆蓋 0.19%（待重抽）；**52.8% somatic sSNV 在 CN-gain（CN-confounded 53.2%；乾淨集 LOH/neutral 46.8%；來源 `data/sm_locus_master_summary.json` cn_somatic_pct）**；precedence「genetic 勝」僅 1 個實證衝突 + n=8 甲基測試（未大規模壓測）；分子證據非 single-cell confirmation。**對外勿稱「甲基偵測 subclone / genome-wide tree / 對手缺檢定」。**
