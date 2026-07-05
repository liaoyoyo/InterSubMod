<!--
建立時間: 2026-07-04
更新: 2026-07-04 v2 — 已實作+7樣本重跑+HTML更新;更正 v1 的 over-claim(「117可救/118→1」錯,基於有損 population 向量;正確用 pairwise 四型=30可救)
類型: 方法正確性修正 (topology_analysis.py solve_topology 的 altset-row-laminar 檢定過嚴)
狀態: LANDED — solve_topology+candidate_scoring+workstation 已改, 7 樣本重跑, 回歸乾淨
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json + sm_region_integration.json (本輪 Read-back+回歸 diff), scripts/topology_analysis.py, candidate_scoring.py, build_topology_workstation.py
-->

# incompatible 重分類：altset-row-laminar 過嚴 → pairwise 四型分流(救30 + 70 recurrence + 18 真artifact)

> 2026-07-04 從純數學形式化(unit-flip / perfect-phylogeny / Steiner tree)反查程式碼揪出並修正。**已實作落地 + 7 樣本重跑 + 回歸乾淨 + HTML 更新。**

## 🔴 誠實更正紀錄(v1→v2)
v1 曾宣稱「117 可救 / incompatible 118→1」— **錯,是 over-claim**。根因:v1 的多路字元四型跑在**有損的 population 向量**(distinct 多位點 genotype)上,漏掉只覆蓋 2 位點的 partial-read 的第 4 型。**正確信號 = pairwise co-read 四型**(`n_independent_clean`,已在 codebase)。用它重算:83/118 有真四型衝突。**真可救只有 30。**

## 核心裁決(L1 已驗)
- `solve_topology` 用 `_laminar()` 檢查 genotype **ROWS** 的 altset 是否巢狀/互斥(L58-61)。
- 🔴 這對 perfect-phylogeny 是**錯的檢定**:正確 = **字元(欄/位置)兩兩四型相容**(Gusfield)。sibling clone 共享祖先突變、各有私有 → altset 天生**非 row-laminar**,即使祖先已觀測(chr5:72534548 鐵證:祖先{2}已觀測 RRA 仍被誤殺)。
- **`n_independent_clean`(pairwise co-read `classify()` 的 `independent`=aa≥2&ra>0&ar>0=全四型) = 正確完整四型信號**;`_laminar`(row) 是過嚴誤判源。

## 118 的正確拆解(pairwise 四型)+ 修正動作
| 類 | 數 | 判準 | 修正後 determinacy | 動作 |
|---|---|---|---|---|
| ✅ 真可救 | **30** | nic==0 & no cycle | A_determined(23)/A_ambiguous(7) | `build_hidden_node_tree` 建隱藏祖先深分支樹 |
| 🟡 非-gain 真四型 | **70** | nic>0 & cn≠gain & no cycle | `recurrence_required(Model A候選;需m通道)` | 重標(Model A;CN=m通道) |
| 🔴 真 artifact | 18 | has_cycle 或 (nic>0 & cn==gain) | incompatible | 維持 |

**HCC1395 headline**: incompatible **118 → 18**;+30 救(23 determined + 7 ambiguous);+70 recurrence-required。
**18 的精確 cross-tab**(對抗驗證校正,原「14 gain+4 cycle」的 gain-first 標法會混淆 gain∩cycle 重疊):**4 cycle-only(非gain) + 6 gain-only(無cycle) + 8 gain∩cycle** = 18;cn 分布 = 15 gain + 3 loh。

## 🔴 對抗驗證揪出並修復的 bug(2026-07-04,workflow wf_520ddd16)
6-agent workflow(5 驗證 + 1 skeptic)獨立重查:surgical-safety / 四型判準 / 分流互斥 / 下游一致性 **全 PASS**;但 algo-correctness + skeptic 獨立揪出 **1 個真 bug**:
- **H1437 chr20:1899981-1922732**:救成 A_determined,但 population `{RAR:20,ARA:12,AAR:3}` 在欄(0,1)**本身四型衝突**(rr=0 無 ancestral 觀測)。
- **根因 = gate/builder 資料源脫鉤**:rescue gate 用 read-level `n_independent_clean`(pairwise co-read,此區=0),但 `build_hidden_node_tree` 吃 **population 向量矩陣**(此區有衝突)→ Gusfield spine emit 同位元翻兩次(homoplasy)的**非法樹**卻標 determined。頻率 1/150 rescued(**HCC1395 的 30 個全乾淨,headline 不受影響**)。
- **修**:`build_hidden_node_tree` 開頭加 **population 欄 rooted 三型自檢**((1,0)(0,1)(1,1) 同現 → 回 incompatible 不建),確保 gate 與 builder 用同一份資料判可建性。
- **修後全 7 樣本 rescued 149,0 個非法樹**(homoplasy/pop-conflict);H1437 chr20 → incompatible(65→66);其餘不變;回歸仍乾淨(僅舊 incompatible 區變動)。

## 逐區鐵證(Read-back + 跑真 solve_topology/builder)
- **chr8:61402948**(`ARAR`,`RRRA`,`AARR`):nic=1 → **非可救**(有真 pairwise 四型衝突,v1 手算漏)→ 落 70 recurrence(loh)。**v1 誤判此為乾淨可救。**
- **chr5:72534548**(`RRA{2}`,`RAA{1,2}`,`ARA{0,2}`):nic=0 → **可救**;祖先{2}已觀測仍被 row-laminar 誤殺 → build_hidden_node_tree 建 determined。
- **chr8:52128513**(gain,`ARRA`,`RRAA`,`AARA`):nic=0(raw n_independent=1 歸 gain-multiplicity)→ 可救 → 建 root→隱藏{3}→分岔。

## 實作(3 檔,已落地)
1. **topology_analysis.py**:+`build_hidden_node_tree(pops)`(Gusfield,補隱藏祖先,label `H_<geno>`);主迴圈 `ttype=="incompatible"` 時三向分流(nic/cycle/cn)。**僅動 incompatible 區→不碰既有 laminar 建樹區(含 107 個 gain+四型建樹區,回歸驗證 0 誤改)。**
2. **candidate_scoring.py**:BASE+recurrence(35);situation/resolution/why_conflict/ncand 加 recurrence 分支。
3. **build_topology_workstation.py**:regSit 鏡射 situation();situation 圖例/卡/位置圖補 recurrence;topology_type 映射補「深分支(推斷祖先)」;更正過時「row-laminar 畫不出深分支/1112 只根分支/incompatible 0.4-19.1%」claim。

## 7 樣本重跑 + 回歸(全乾淨:僅舊 incompatible 區變動,0 誤改非-incompatible)
| 樣本 | inc | 救(det/amb) | recurrence |
|---|---|---|---|
| HCC1395 | 118→18 | 30 | 70 |
| COLO829 | 15→1 | 1 | 13 |
| H1437 | 263→**66** | **26** | 171 |
| H2009 | 812→265 | 66 | 481 |
| HCC1395_DORADO | 55→4 | 11 | 40 |
| HCC1937 | 44→0 | 9 | 35 |
| HCC1954 | 27→1 | 6 | 20 |

〔H1437 為 bug 修復後值:chr20:1899981 從 A_determined 改回 incompatible。全 7 樣本 rescued 合計 **149**(修前 150,扣掉該 1 個 homoplasy 假樹)。〕跨樣本 incompatible% : 0.0%–6.2%(原 0.4%–19.1% 含 row-laminar 誤殺)。

## 科學立場(誠實標註)
- **30 救**=SOLID(perfect phylogeny 證明存在;pairwise 全四型相容);代價=**posits 未觀測祖先 clone**(標準 perfect-phylogeny + 對齊 Steiner 形式化,但須在論文標「inferred ancestral clone」非觀測)。
- **70 recurrence-required**=Model A 候選,**弱主張**:非-gain 四型 → 需信任 CN(SEQC2)為獨立 m-通道拆 artifact(m>1) vs 真recurrence(m≤1);未證真 recurrence,只是「非乾淨 conflict、非乾淨 determined」的中間態。
- **未變**:HCC1395 ⭐3 封頂結論、甲基 bounded-auxiliary、「定不出來即答案」不受影響(此修正只把「假 incompatible」正名,不新增乾淨 subclone 證據)。

關聯:形式化 spec `20260704_formal_problem_statement_topology_from_cooccurrence_01.md`(Model A §4/§9);`20260628_subclone_reconstruction_master_spec_01.md`(§6/§9 determinacy census,需同步)。

---

## 🆕 gap#3(2026-07-04 已落地）— Model A m-通道拆分 70 recurrence 區
形式化 §4/§9 的 Model A m-通道由「粗 `cn!='gain'` categorical 代理」升為**整數 CN**（SEQC2 `ngs_benchmark_cnv_gain_cn.bed` CN 3/4/5 + `_loss_cn.bed`）。HCC1395 70 recurrence → **三分流**（自算=實作輸出，逐區吻合）:
| m-通道判決 | 區數 | 判準 | 處置 |
|---|---|---|---|
| `recurrence_artifact(m>1;CN-amp)` | **11** | 整數 CN≥3（gain_cn）| 棄（粗 gate 漏掉的 GAINLOH，如 chr8:61402948 CN3 / chr6:89339042 CN4）|
| `recurrence_candidate(m=1)` | **9** | diploid CN2 / loss CN≥1 | 留（真 recurrence 候選；含 chr14:16093807 原「真四型衝突」）|
| `recurrence_LOH_unresolved` | **50** | copy-neutral LOH total=2、無 allele-specific CN | 未解（VAF L3：35 likely_artifact / 15 likely_recurrence，**不硬判**）|
**最大 win**：11 GAINLOH 是粗 categorical 漏掉的多重度 artifact，接整數 CN 後正確判棄 → 推翻「70 全候選」。**誠實上限**：50/70 copy-neutral LOH 因 SEQC2 無 allele-specific CN 無法定 m∈{1,2}，只給 VAF L3 軟旗標；6 非-HCC1395 樣本 cn=unknown → `recurrence_required(m通道不可用)`（優雅降級）。實作:`topology_analysis.py:m_channel_split` + candidate_scoring 4 子標籤 + workstation m-通道判決框。**additive:不動 A_determined 1764 / 3885 分母。**

## 🆕 gap#2(2026-07-04 已落地）— 全枚舉 + identifiability
`enumerate_candidate_trees.py` 舊只枚舉 per-edge 順序、對「無 multi-flip 邊」的 parent-choice 歧義 **silently `return None`（32 區憑空消失）**、`B1_equiprobable=0` **假稱無歧義**。修:`+parent_choice_enum`（枚舉等大觀測父）+ isomorphism dedup → `n_minimal_trees` + `identifiability`。HCC1395:
- **69 A_ambiguous 全枚舉**（舊只 37，+32 捕捉）→ **44 枚舉後收斂 determined + 25 真歧義**（修正假「0 歧義」）。
- **完整 identifiability_table 分區 3885**：determined 1808 / ambiguous 25 / recurrence 70 / conflict 18 / nonenumerable(B_pairwise+C) 1487 / other 477。
**additive:此為 candidate_trees.json 上層,不改 topology 的 determinacy census（determinacy 仍 A_ambiguous 69）。** 🔴 排序絕不用甲基。

## 🆕 gap#1(partial-read subcube；2026-07-04 已落地 canonical，controlled 重跑）
**Provenance 先確認**：當前 C2-fixed `sm_region_integration.py` 跑 20260618 census **完美重現 canonical nic**（全 7143 區 0 不一致）→ 20260618 census = canonical 源,安全。
**Stage 1**：`sm_multilocus_combinations.py` 加 `subread_groups`（'X'=未覆蓋=cube face）+ `col_coverage` + 不再因無 full-cov 就丟整組;BAM 重跑（chr1-22）→ 7100 groups,**2885 空-population 但有 subread**。`sm_region_integration.py` 傳遞 subread_groups。
**Stage 2**：`topology_analysis.py` — 無 full-cov population 但有 ≥2-位點 partial read(co-occurrence substrate)= **subcube-recovered**,以 tree_shape/nic 分類為 `E_subcube_recovered(pairwise樹/欠定/含衝突)`,🔴 **不僭稱 A_determined(無 full-cov single-molecule)**。
**HCC1395 結果**（controlled 重跑,20260618 workspace 驗證後 promote 4axis canonical）:
- **救回 2403 區**（先前 empty 掉出重建）:pairwise樹 **1561** + 欠定 676 + 含衝突 166。
- **denominator 3885 → 6288**;**A_determined 1764 完全不變**（full-cov census 0 regression）。
- 這是誠實的 headline 成長:救回區透明標 E_subcube_recovered、與 full-cov 3885 分開,不混入 A_determined。
🔴 **僅 HCC1395**（6 樣本 cn/subread 未重跑,其 topology subcube 空 → 不受影響,graceful）。備份 `scratchpad/CANON_pregap1_*`。
