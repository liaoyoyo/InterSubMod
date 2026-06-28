<!--
建立時間: 2026-06-29
類型: methodology final — 骨幹欠定問題解決層 + 甲基註釋角色窮舉 + 結構穩定性 + 待定義 gap
狀態: concluded(全問題已有可執行方法/界定;結構穩定性誠實量化;9 gap 列出)
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/backbone_resolution.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/backbone_stability_audit.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/h3_methyl_phasing.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/h3_unresolved_grouping.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/eps_minread_sensitivity.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/sm_phaseset_extension.json
provenance: HCC1395 ⭐3 單樣本;凍結 topology @ feat/summary-nreadsvalid@5308d9e;4-角度對抗完整性審查 workflow wf_857b82d9;每數字 grep-able
-->

# 骨幹欠定問題解決層 + 甲基註釋角色 + 結構穩定性 — 最終定義

> 框架：Verdict-Pyramid。HCC1395 ⭐3 單樣本。回答「6 問題是否都有可執行解法 + 甲基還有哪些角色 + 最終結構穩定性 + 哪些細節需定義」。每數字 grep-able（§13-A）。經 4-角度對抗完整性審查（workflow）。

## §0 TL;DR（誠實 headline）

1. **6 個問題全有可執行方法或明確界定**（成環/欠定/未定/H3/ε門檻/Tier-PS）——但多數是「界定為資料極限 + 誠實標不確定」，非「全部解決」。
2. **甲基註釋角色窮舉完畢（13 個），無遺漏**：可用 = HP 定相(germline 軸) + 負向篩選 + 有界 characterization；不可用 = 定群/specificity/排序/外推。與 bounded-auxiliary 一致。
3. **🔴 結構穩定性關鍵修正**：「47% determined」是**軟上界**；同時 **非循環+TP-backed+多位點(≥3)+非CN-gain** 的**真正穩固核心 = 57 區（1.5%）**。CN-gain multiplicity 是最大威脅。
4. **9 個待定義 gap**（denominator/incompatible 重算/isolated 計數/截斷/覆蓋功率/輸出格式…）已列優先序。

## §1 六個問題 — 成因 + 可執行解法（驗證數字）

| # | 問題 | 成因（數據）| 解法 / 裁決 |
|---|---|---|---|
| **1** | 衝突/成環 | 12 區;**verified incompatible = 0**（10/12 npop≤2 不可能真違反、8/12 genotype 截斷、9/12 CN-gain）| ✅ **無真樹衝突訊號**（reassuring）→ 標 artifact 不建樹 + 修截斷/分類邏輯（`backbone_stability_audit.json`）|
| **2** | ε/min-read 減資訊? | MINREAD=3 每區丟 0.48 ALT 群;singleton 1230 中 **76% 與大群差 1-Hamming = 定序錯誤**;coherent 933(可能真稀有)| ✅ 合理 specificity/sensitivity 取捨（移除多數噪聲）;933 coherent 為已知邊界（`eps_minread_sensitivity.json`）|
| **3** | 多樹相容欠定 有無方法? | **550 區全只有單一 ALT 群**（非「多樹歧義」是「只觀測到一群」）;parsimony resolvable=0 | ⚠️ **需深覆蓋**觀測中間群（非演算法可解）;甲基串接=ordering=L3 弱不可靠 |
| **4** | 未定/缺中間群 用最小突變定序? | A_ambiguous 76;parsimony 第一順位中位機率 **0.5、0 高信心、55 缺中間群** | ✅ parsimony(最大簡約)是**正確原則**但**信心低**（跳步太不確定）→ 給第一順位+標 L3（`backbone_resolution.json`）|
| **5** | 相容但欠定 如何定義? | 同 #3，單群為主 | ✅ **輸出=相容樹集合+第一順位(機率)**，非硬塞單樹（誠實輸出不確定性）;深覆蓋為唯一增益 |
| **6** | H3-unphased 用甲基定 HP? | 93 區;**可測 10→10/10 指派既有 HP、0 真第三群**;45 區 HP1≈HP2 無 germline-ASM 無法定相;38 缺一 HP | ✅ **支持你的理論**（H3=定相失敗非新群）;但僅 germline-ASM 存在處可定（`h3_methyl_phasing.json`）|

### §1.1 Q6 延伸（45 無-ASM + 38 缺-HP 能否確認「就是一群」）
`h3_unresolved_grouping.json`（83 區）：
- **63/83 ALT reads 甲基 unimodal**（無次結構）→ 支持「同一群」，但 **= absence of evidence ≠ proof**。
- **mut_methyl_coherent = 0/80**：這些無-ASM 區，ALT vs REF 甲基根本不可分 → **甲基在此完全靜默**（連 corroborate 都做不到）。
- **13 區 ALT 集中單一 HP** → 「誰帶突變」分群**基因型已定**，甲基只需確認無次結構（負向篩選）。
- 結論：對 45+38 區，甲基能做的 = **負向篩選支持單群**（弱、非證明）;不能定相、不能 corroborate。

### §1.2 Tier-PS（遠距 same-PS partner 補連接）
`sm_phaseset_extension.json`：same-PS>50kb 對 80,980 同-HP 候選，但 **CCF 僅 42.4% 一致** → **PS = germline 單倍型 context、非克隆連鎖**。
→ ⚠️ **部分不合理**：能排除跨-HP + 放置 isolated，**不能**補 cell-level 巢狀。

## §2 甲基註釋角色窮舉（13 個，無 sound-但-未測）

| 角色 | 裁決 | 依據 |
|---|---|---|
| R1 定義群數/指派 read | ❌ 不可用（循環 double-dip）| `project_subcluster_cluster_count_determination` |
| R2/① 群間甲基距離 | ⚠️ 有界（純描述、cis-confounded）| concordance R²≈0.17 |
| R3/② 突變 read 群聚 corroborate | ⚠️ 有界-弱（6.6% corroborate、cis 0/740）| reextract;本輪 mut_coherent 0/80 |
| R5 cis-control 立 specificity | ❌ 不可用（structural;needs_methyl∩乾淨≈0）| cis_control_verdict |
| R6/③ 表觀距離排 root→leaf | ⚠️ L3 弱未顯著（ρ0.18 p0.06）| ordering_pilot |
| R8/⑤ 負向篩選(無次結構) | ✅ **可用**（unimodal 95.8%）/ L3 候選 / 確認不可 | aux_annotation |
| R9/⑥ cluster-count sanity | ⚠️ 有界（=R3+R8）| 同上 |
| R10/④ germline HP 定相輔助 | ✅ **可用-有界**（10/10、headline 0.885）| h3_phasing + V1-V11 |
| ③' 甲基外推到未觀測位點 | ❌ **foreclosed**（已觀測都 recover 0/弱→無可外推訊號）| concordance |
| 獨立估 CCF / 跨區表觀連結 subclone / fCpG 絕對時序 | ❌/reopen（單-bulk 無功率，需 C1 多樣本/C3 single-cell）| — |

**窮舉結論**：**沒有任何 sound 角色未測**。可用 = ④HP 定相 + ⑤負向篩選 + 有界 characterization;其餘界定或 foreclosed。對外口徑勿主張「甲基驅動重建/定群/定樹/cis-control 已驗 specificity」。

## §3 結構穩定性（🔴 誠實量化 — robustness ladder）

`backbone_stability_audit.json`：

| 層 | 條件 | n | %(of 3885) |
|---|---|--:|--:|
| L1 | A_determined（單分子向量）| 1,812 | 47% |
| L2 | +TP-backed（fp=0,tp>0）| 931 | 24% |
| L3 | +多位點（n_sSNV≥3）| 260 | 6.7% |
| **L4** | **+非 CN-gain（真正穩固核心）** | **57** | **1.5%** |

🔴 **「47% determined」是軟上界**。determined core 組成脆弱：**79% 只有 2 位點**（非多節點樹，只說 2 突變共享/分裂一個 lineage）、**69.3% CN-gain**（multiplicity artifact 區）、**43.5% 無 truth label**。

**CN-gain inflation = 最大威脅**：determinacy rate 在 CN-gain **最高（53.9%** vs loh 36% / neutral 44.5% / loss 10.8%）→ multiplicity 可能製造假共現→假 determined。任何 tree-level claim 前應 mask/控制 CN-gain。

**denominator**：47% = 1812/3885（有向量區）= **25%/7143（全區）**。引用必標分母。

**穩 vs 上界**：
- **穩（reproducible）**：incompatible=0(真衝突)、B1 max-prob 0.5/0 高信心、B2=550 覆蓋限制、topology byte-可重現（已修 tie）、A-framing 質性。
- **上界（會縮）**：47% determined、把 2-位點區當「subclone tree」、CN-gain determined 計數、任何用 3885 而非 7143 的 %。

## §4 待定義 gap（9 個，優先序）

| Gap | 級 | 問題 | 必做 |
|---|---|---|---|
| G1 denominator | 🔴 | stats(7143) vs detail(3885) 兩套矛盾數 | 定**單一 canonical 分母**;對齊 incompatible 22vs12 等 |
| G2 incompatible 重算 | 🔴 | 12 全 artifact（verified 0）| 改判準 npop≥3 + 有違反對;報 verified=0 |
| G3 genotype 截斷 | 🟠 | cap=8 截斷 **42 區**（31 帶 TP、8 incompatible）| 提高 cap 或標 `truncated` 排除 |
| G4 isolated 計數 | 🟠 | 8320/11520/266 三套、不在 topology output | 定單一定義 + 加呈現列 |
| G5 覆蓋功率 | 🟡 | 「深覆蓋可解」未量化 | 算 per-region 目標深度 + 2x/4x 預期 resolved |
| G6 輸出格式 | 🟡 | 無最終呈現 contract | 定各類(determined/ambiguous/…/H3/isolated)圖表呈現 |
| G7 inconsistent 計數 | 🟡 | 22/12/5 跨檔不一 + CN-clean subset 未宣告 | 宣告報告 subset |
| G8 B1 set-vs-tree | 🟢 | 76 全 coin-flip 無呈現規則 | 定「相容集合 vs 單樹」呈現 |
| G9 provenance stamp | 🟢 | 計數未標 (ε,MINREAD) | output 加 ε=0.02/MINREAD=3 + 敏感度帶 |

## §5 最終輸出結構定義（誠實分層機率）

```
每區 →
├ A_determined(1812,47%軟上界;真正穩固 L4=57)：單一樹
├ A_ambiguous(76)：parsimony 第一順位 + 機率(中位 0.5,L3)
├ B_pairwise(958)：拼接結構,標非單分子
├ C_underdetermined(550,全單群)：相容集合 + 標「需深覆蓋」
├ incompatible(0 verified)：artifact,不建樹
├ isolated(待定義單一計數)：Tier-PS 放上單倍型(非巢狀)
└ H3-unphased(93)：germline-ASM 在→甲基定相(10/10);否→標 unresolved
附 CN-gain mask 警示 + 甲基軟標(負篩/HP定相,L3)
```

## §6 數字溯源（§13-C）
| 數字 | 值 | 來源 |
|---|---|---|
| robustness ladder L1/L2/L3/L4 | 1812/931/260/57 | `backbone_stability_audit.json` |
| core 組成 n2/gain/notruth | 79%/69.3%/43.5% | 同上 |
| verified incompatible | 0 | 同上 |
| determinacy rate by CN | gain53.9/loh36/neu44.5/loss10.8 | 同上 |
| 截斷 | 42(31 TP) | 同上 |
| 成環/npop≤2 | 12/10 | `backbone_resolution.json` |
| B1 中位機率/缺中間群 | 0.5/55 | 同上 |
| B2 單群 | 550/550 | 同上 |
| H3 可測/指派 | 10/10 | `h3_methyl_phasing.json` |
| Q6-ext unimodal/coherent/單HP | 63/0-of-80/13 | `h3_unresolved_grouping.json` |
| eps 多救/1-Hamming% | 1863/76% | `eps_minread_sensitivity.json` |
| Tier-PS CCF 一致 | 42.4% | `sm_phaseset_extension.json` |
