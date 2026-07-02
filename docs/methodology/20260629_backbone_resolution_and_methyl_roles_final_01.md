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

> **📍 2026-07-02 correction（此檔原為 06-29 快照，本輪已回填 canonical）**：robustness ladder 已更新為 **07-02 重生值 L1-L4 = 1741/889/218/34**（源 `backbone_stability_audit.json` 07-02，C2/C3/D4 後）。**HCC1395 canonical determinacy = A_determined 1741 / incompatible 118 / A_ambiguous 62 / B_pairwise 943 / C_underdetermined 544**（`topology_per_region.json` 07-01；權威 `20260701_ssnv_backbone_method_spec_and_correctness_audit_01.md` + `20260628_..._master_spec_01.md §6`）。⚠ 內文其餘散在的 06-29 診斷數字（upstream has_cycle 22、§1 表的 incompatible/A_ambiguous）如與上列衝突，**以 canonical 為準**。**定性結論（軟上界 / CN-gain multiplicity 為最大威脅 / 甲基 bounded-auxiliary）不變**。

## §0 TL;DR（誠實 headline）

1. **6 個問題全有可執行方法或明確界定**（成環/欠定/未定/H3/ε門檻/Tier-PS）——但多數是「界定為資料極限 + 誠實標不確定」，非「全部解決」。
2. **甲基註釋角色窮舉完畢（13 個），無遺漏**：可用 = HP 定相(germline 軸) + 負向篩選 + 有界 characterization；不可用 = 定群/specificity/排序/外推。與 bounded-auxiliary 一致。
3. **🔴 結構穩定性關鍵修正**：「44.8% determined」是**軟上界**；同時 **非循環+TP-backed+多位點(≥3)+非CN-gain** 的**真正穩固核心 = 34 區（0.9%）**（07-02 canonical；修前 57/1.5%）。CN-gain multiplicity 是最大威脅。
4. **9 個待定義 gap**（denominator/incompatible 重算/isolated 計數/截斷/覆蓋功率/輸出格式…）已列優先序。

## §1 六個問題 — 成因 + 可執行解法（驗證數字）

| # | 問題 | 成因（數據）| 解法 / 裁決 |
|---|---|---|---|
| **1** | 衝突/成環 | **上游 has_cycle = 22 真 cycle**（全 sSNV pairwise 圖,如雙向 nesting;12 有向量+10 無向量;**77% CN-gain**）;stored 截斷(cap=8)查不到（≠「0 真 cycle」,前述「verified=0」是測錯截斷資料的誤判,已修 G2）| ✅ cycle 真實但 **77% = CN multiplicity artifact 非真演化衝突** → 標 likely-CN-multiplicity-artifact 不建樹;真驗證須提高 cap + CN-aware 模型（`backbone_stability_audit.json`）|
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
| L1 | A_determined（單分子向量）| 1,741 | 44.8% |
| L2 | +TP-backed（fp=0,tp>0）| 889 | 22.9% |
| L3 | +多位點（n_sSNV≥3）| 218 | 5.6% |
| **L4** | **+非 CN-gain（真正穩固核心）** | **34** | **0.9%** |

〔07-02 重生 canonical；修前 06-29：1,812/931/260/57 = 47/24/6.7/1.5%〕

🔴 **「44.8% determined」是軟上界**。determined core 組成脆弱：**81.2% 只有 2 位點**（非多節點樹，只說 2 突變共享/分裂一個 lineage）、**70.9% CN-gain**（multiplicity artifact 區）、**44.1% 無 truth label**。

**CN-gain inflation = 最大威脅**：determinacy rate 在 CN-gain **最高（53.9%** vs loh 36% / neutral 44.5% / loss 10.8%）→ multiplicity 可能製造假共現→假 determined。任何 tree-level claim 前應 mask/控制 CN-gain。

**denominator（G1 已修 canonical）**：determinacy 一律以 **3885（有向量區）**為分母 → 44.8%=1741/3885 = **24.4%/7143（全區覆蓋脈絡）**；全 7143 區覆蓋改看 `region_coverage`（with_vector 3885 / germline_only 371 / no_vector 2887）。stats.determinacy 現 == detail（sum 一致），不再有 7143-混算的 incompatible 22/B 2760/C 1658。

**穩 vs 上界**：
- **穩（reproducible）**：incompatible=0(真衝突)、B1 max-prob 0.5/0 高信心、B2=550 覆蓋限制、topology byte-可重現（已修 tie）、A-framing 質性。
- **上界（會縮）**：44.8% determined、把 2-位點區當「subclone tree」、CN-gain determined 計數、任何用 3885 而非 7143 的 %。

## §4 待定義 gap（G1-G9 全處理完畢 06-29）

| Gap | 級 | 問題 | 狀態 |
|---|---|---|---|
| G1 denominator | ✅ **已修** | stats(7143) vs detail(3885) 兩套矛盾數 | **determinacy canonical=3885**（stats==detail，sum 一致）;全區覆蓋改 region_coverage（3885/371/2887）;master spec §6 已更（棄用舊 B2760/C1658/incomp22）|
| G2 incompatible 重算 | ✅ **已修** | 「12 全 artifact/verified 0」誤判 | 真相 = 上游 has_cycle **22**（12 有向量+10 無）真 cycle、**77% CN-gain**=multiplicity artifact;stored 截斷查不到非「0 真 cycle」;detail 加 `cycle_cause` 註記 |
| G3 genotype 截斷 | ✅ **已處理** | cap=8 截斷 42 區 | detail 加 `truncated`/`genotype_len` 旗標;42 區(31 TP/8 incompat/僅 3 在 core)→ core 標 determined-on-subset,主結論穩（§7）|
| G4 isolated 計數 | ✅ **已收斂** | 8320/11520/266/13778 四套 | 測不同東西,定 **canonical=8320**（n_partners=0）+ 各定義（§7,`gaps_g4_g5_resolution.json`）|
| G5 覆蓋功率 | ✅ **已量化** | 「深覆蓋可解」未量化 | underpowered 5458;2x 救 39.4%/4x 56.7%;但 43%(2364)幾何需長 read（§7）|
| G6 輸出格式 | ✅ **已定義** | 無呈現 contract | 各類呈現表（§5+§7 G6 表）|
| G7 count 宣告 | ✅ **已宣告** | 計數跨檔不一 | incompatible 22(全)/12(有向量);骨幹錨 CN-clean full_tree 205/677;determinacy 錨 3885（§7）|
| G8 B1 set-vs-tree | ✅ **已定規則** | 76 coin-flip 呈現 | parsimony 第一順位 + 標「≥2 相容樹、機率≈0.5」（§7）|
| G9 provenance stamp | ✅ **已加** | 計數未標參數 | topology output 加 `provenance`（cap8/eps0.02/min3/coread6/分母3885/byte-repro）（§7）|

## §5 最終輸出結構定義（誠實分層機率）

```
每區 →
├ A_determined(1741,44.8%軟上界;真正穩固 L4=34)：單一樹
├ A_ambiguous(76)：parsimony 第一順位 + 機率(中位 0.5,L3)
├ B_pairwise(958)：拼接結構,標非單分子
├ C_underdetermined(550,全單群)：相容集合 + 標「需深覆蓋」
├ incompatible(has_cycle 22,77% CN-gain)：likely CN-multiplicity artifact,不建樹
├ isolated(待定義單一計數)：Tier-PS 放上單倍型(非巢狀)
└ H3-unphased(93)：germline-ASM 在→甲基定相(10/10);否→標 unresolved
附 CN-gain mask 警示 + 甲基軟標(負篩/HP定相,L3)
```

## §7 Gap 解決（G3-G9，2026-06-29 全處理）

### G3 genotype 截斷（cap=8）✅ 已加旗標
- detail 每區加 `truncated`（n_sSNV>genotype_len）+ `genotype_len`;**42 區截斷**（31 帶 TP、8 為 incompatible、**僅 3 在 A_determined core**）。
- 處置：core 的 3 區標「determined-on-subset(8 位點)」;incompatible 8 區的 cycle 落在 >8 位點（與 G2 一致）。**主結論穩**（core 僅 3/1741 受影響）。真解須提高 upstream cap（需重跑 pipeline,非本層）。

### G4 isolated 計數收斂 ✅（`gaps_g4_g5_resolution.json`）
四套數字測**不同東西**（非矛盾,是混用未標基底）：

| 數字 | 基底 | 角色 |
|--:|---|---|
| **8,320** | per-sSNV;`n_partners_le50k==0`（read-span 內無 partner）| 🔵 **CANONICAL isolated**（骨幹連鎖本義;TP 7555/FP 765）|
| 11,520 | 全 loci tree_role（含 germline,廣基底）| 不同分母,非 isolation 本義 |
| 266 | 只 somatic-confirmed loci tree_role | somatic 子集 |
| 13,778 | 只含 1 sSNV 的**區**（區層級）| 與 isolation 正交 |

→ 對外一律用 **8,320**;isolated loci **不在 topology detail**（無 tree）→ 需獨立「isolated 表」（caller VAF + 可能 Tier-PS 放單倍型,非巢狀）。

### G5 覆蓋功率（量化「深覆蓋可解」）✅
- underpowered（有 partner 無 co-read link）= **5,458**。
- 線性估救回：**2x = 2,151（39.4%）/ 4x = 3,094（56.7%）**。
- 🔴 **但 2,364（43%）max_coread=0 = 幾何限制（span>read）→ 覆蓋救不了、需更長 read**（非深度問題）。
- → 「深覆蓋可解」精確版：**~57% 區 4x 可救（depth-limited）;~43% 是幾何（需長 read）**。C_underdetermined 550 單群亦受此限。

### G6 最終輸出格式 contract ✅（定義見 §5;每類呈現）
| 類別 | 呈現 | 標註 |
|---|---|---|
| A_determined 1741（穩固 34）| 單一樹 + lineage 標籤 | 標 robustness tier(L1-L4)+truncated flag |
| A_ambiguous 76 | 單樹 + parsimony 第一順位 | 標機率(中位 0.5,L3)|
| B_pairwise 958 | 拼接結構 | 標非單分子 |
| C_underdetermined 550 | 「需深覆蓋」+ 救回估 | 標 depth vs 幾何 |
| incompatible 12(has_cycle 22) | 不建樹 | 標 likely-CN-multiplicity-artifact |
| isolated 8320 | 獨立表(無 tree)| caller VAF + Tier-PS 單倍型 |
| H3 93 | germline-ASM 在→甲基定相 | 否則標 unresolved |

### G7 count 宣告 ✅
- **incompatible**：has_cycle 全 22 / 有向量 12（report 對外用「22 真 cycle、12 在 tree-detail」,勿混）。
- **structure / CN-clean subset**：對外引「reconstruction 骨幹」錨定 **CN-clean full_tree 205 / 全 full_tree 677**（master spec §6）;determinacy % 錨 3885（G1）。

### G8 B1 set-vs-single-tree 規則 ✅
A_ambiguous 76 全 first_rank_prob≤0.5（0 高信心）→ **呈現規則：給 parsimony 第一順位樹但標「≥2 相容樹、機率≈0.5」**，不呈現為單一定論。

### G9 provenance stamp ✅
topology_per_region.json 加 `provenance`：genotype_cap=8 / eps=0.02 / min_read=3 / coread≥6 / canonical 分母 3885 / byte_reproducible=true / tiebreak 紀錄。→ determinacy 計數現可重現且綁定參數。

## §6 數字溯源（§13-C）
| 數字 | 值 | 來源 |
|---|---|---|
| robustness ladder L1/L2/L3/L4 | 1741/889/218/34 | `backbone_stability_audit.json`（07-02 重生）|
| core 組成 n2/gain/notruth | 81.2%/70.9%/44.1% | 同上（07-02 重生）|
| incompatible has_cycle(全/有向量/CN-gain%) | 22 / 118 / 77% | 同上 incompatible_recompute（07-02；有向量 12→118 隨 C2）|
| determinacy canonical sum(==detail) | 3885 | 同上 denominators |
| region_coverage(向量/germline/無向量) | 3885/371/2887 | 同上 |
| determinacy rate by CN | gain53.9/loh36/neu44.5/loss10.8 | 同上 |
| 截斷 | 42(31 TP) | 同上 |
| 成環/npop≤2 | 12/10 | `backbone_resolution.json` |
| B1 中位機率/缺中間群 | 0.5/55 | 同上 |
| B2 單群 | 550/550 | 同上 |
| H3 可測/指派 | 10/10 | `h3_methyl_phasing.json` |
| Q6-ext unimodal/coherent/單HP | 63/0-of-80/13 | `h3_unresolved_grouping.json` |
| eps 多救/1-Hamming% | 1863/76% | `eps_minread_sensitivity.json` |
| Tier-PS CCF 一致 | 42.4% | `sm_phaseset_extension.json` |
| G3 截斷區/在core/帶TP | 42/3/31 | topology detail truncated |
| G4 isolated canonical | 8320(TP7555/FP765) | `gaps_g4_g5_resolution.json` |
| G5 underpowered/2x救/4x救/幾何 | 5458/2151/3094/2364 | 同上 G5 |

## §8 補充交付（Tier-R 刻畫 + 4-gamete AA 歸屬 + 佇列 3 欄；06-29）

### §8.1 Tier-R 刻畫（#3 underpowered + #4 isolated）`tier_r_characterization.json`
- **underpowered 5,458**：100% 有 CCF（census VAF）→ clonal 2,992 / mid 2,174 / low 292；depth 2x 救 2,151(39.4%)/4x 3,094(56.7%)；**2,364(43%) max_coread=0=幾何(span>read)覆蓋救不了需長 read**。
- **isolated 8,320**：census 無共讀 VAF → 補 **caller VCF AF（FORMAT.AF）→ 100% 可刻畫**：clonal 4,702 / mid 3,218 / low 400（TP 7,555/FP 765）。
- 結論：兩類**非無法處理**，有 CCF/caller VAF 可放 clonal 譜；差別僅能否建樹。

### §8.2 4-gamete AA 歸屬（用戶提議：AA 屬 AR 還是 RA 譜系？）`four_gamete_aa_resolution.json`
- 51 incompatible/noise 區 → 37 個 4-gamete pair。**same-HP 前提**：14 成立 / **23(62%) 跨HP失敗**（AA=doublet/CN artifact 非譜系後代）→ **前提檢查必要**。
- same-HP ∩ CN-clean = **0**（14 個 same-HP **全 CN-gain**，13/14 AA 膨脹=CN multiplicity 訊號）→ VAF 不可靠。
- 裁決：方法在 same-HP+CN-clean 時可行，但 **4-gamete 本身壓倒性是 CN-gain artifact**（呼應 G2），此單樣本 0 例乾淨可用。

### §8.3 確認佇列 3 欄（#1）`candidate_scoring.json` + 工作站顯示
- **why_conflict**：成環=cycle_cause(9 CN-gain/3 other)、550 單群缺連接、76 缺中間群、99 跨HP、903 pairwise。
- **parsimony_first_rank_prob**：76 ambiguous 區。**methyl_applicability**：負篩可用/排序 L3 弱/specificity 不可；H3 區標 germline-ASM 定相。
- 工作站佇列顯示這 3 欄（commit 7d95436）。

### §8.4 樹節點歸類驗證（#2）
✅ **1,851/1,867（99.1%）樹形區所有 ALT 群正確歸節點**；16 區 drop 噪聲群（17 群）；HP 兩樹 20 區正確拆；甲基分支定位 L3 弱。
> commits: 7d95436（#1/#3/#4）+ 7958406（AA 歸屬）。
