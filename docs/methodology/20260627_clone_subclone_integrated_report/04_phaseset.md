<!--
建立時間: 2026-06-27
類型: L3 phase-set — 單變量（clone/subclone 整合報告）
狀態: in_progress（HCC1395 ⭐3）
data_sources: data/sm_phaseset_extension.json, data/region_ps_flag.json
-->

# L3 — phase-set（PS）：reliability 旗標 + Tier-PS 延伸評估（單變量）

> 兩個用途（依 `00b` grounding）：① **PS-reliability 旗標**（多 PS region = HP 判別不可信）② **Tier-PS 是否延伸克隆連鎖**。
> 圖：`figures/04_phaseset.png`。資料：`data/sm_phaseset_extension.json` + `data/region_ps_flag.json`（逐區域旗標）。
> PS 抽取：23,099/23,810 somatic sSNV 有 PS（97%）。

## §1 PS-reliability 旗標（credibility 精煉）

每區域的 somatic sSNV 跨幾個 PS？單一 PS = 同 phase block（HP 一致可信）；跨 >1 PS = phase-switch 風險（HP same/diff 判別不可信）。

| | 數量 | 比例 |
|---|---|---|
| **single_ps_reliable** | **6,623** | 92.7%（/7,143 全區）|
| multi_ps_uncertain | 449 | 6.3% |
| no_ps（無 PS 資訊，同樣不可 gate）| 71 | 1.0% |

（修對抗稽核：分母用全 7,143 區；含 71 個 no_ps 區 → single-PS 可信 = 6,623/7,143 = **92.7%**；不可信或未知 = 520/7,143 = 7.3%。）

**結構區（reliable / uncertain）**：full_tree 596/81 · sibling_only 1,130/105 · linear_nested 1,784/124 · co_linked 812/46。

→ 🔑 **94% 區域在單一 phase block，HP 判別可信**；**6%（449）跨 PS = HP 不可信，最乾淨的 sibling subclone 結論應排除之**（逐區旗標在 `region_ps_flag.json`）。這正是 `00b` grounding 要求的「problem-PS 排除」credibility 保護。

## §2 Tier-PS 是否延伸克隆連鎖？— 不延伸（data-supported）

same-PS 但 >50kb apart 的 somatic sSNV 對（跨 Tier-R region，**非同分子**）= Tier-PS 候選。問：它們是否屬同克隆？用 CCF 一致性檢驗。

| 指標 | 值 |
|---|---|
| same-PS >50kb 候選對（capped sample）| 127,183 |
| 其中 same-HP | 80,980（64%）|
| CCF 一致（\|ΔVAF\|≤0.1）| **53,939 / 127,183 = 42.4%** |

→ 🔑 **同 PS（同 germline 單倍型 block）的遠距 sSNV 只 42.4% CCF 一致**（58% 不一致）= **它們在同一 germline 單倍型上、卻屬不同 subclone**。
→ **PS（germline phasing）不建立克隆共現**：same-PS ≠ same clone。Tier-PS 只給「germline 單倍型 context」，**不延伸克隆連鎖**（符合 `00b` precedence：phasing 是 germline 層級，非 clonal）。

## §3 data-supported 結論
1. **PS 的價值 = reliability 旗標**：94% region HP 可信、6% 需排除 → 提升 sibling subclone 結論的可信度。
2. **PS 不延伸克隆連鎖**：CCF 一致率僅 42.4% → 同單倍型 ≠ 同克隆。故克隆連鎖**只認 Tier-R（same-read 同分子）**，PS 不升級成 clonal linkage。

## §4 限制
- Tier-PS 候選用 capped sample（giant PS block 取最密 40，hard cap 300k）→ 比例為估計；方向（CCF 多不一致）穩健。
- PS 本身受 longphase phasing 誤差（problem block）影響；本層正是用來**標記**這類風險。
