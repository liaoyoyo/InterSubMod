<!--
建立時間: 2026-07-01
類型: methodology 定稿 — sSNV 單分子共現重建骨幹「方法規格 + 邊界正確性稽核」
狀態: concluded(方法規格已逐階段對源碼;13 邊界確認;3 must-fix 已獨立驗證)
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/sm_region_integration.json,.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot/sm_linkage_genomewide.json
provenance: HCC1395 canonical。數字皆本輪 Read-back 獨立驗證(非 subagent 轉述):C3(topology_per_region.json drop_noise_frac)、C2(sm_linkage pairs rel==independent)、D2、c群分布、pairwise 救回量。source 檔為 gitignore artifact,腳本已 commit 可重現。
-->

# sSNV 單分子共現重建骨幹 — 方法規格 + 邊界正確性稽核

> 框架：Verdict-Pyramid（裁決先行 → 規格 → 證據 → 邊界）。每個數字標來源檔（§13-C）+ 證據級（scientific-rigor §2）。
> 對象：`docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/` 五支腳本。目的：改效能前先確認方法正確性。

---

## §1 TL;DR（裁決）

| 問題 | 裁決 | 證據級 |
|---|---|---|
| 方法整體是否 sound？ | ✅ **是**（骨幹概念+實作站得住；多數邊界=誠實解析度降級非 bug；誤差方向多保守） | L2 |
| c 群為何有 0 群？ | 3258 區(46%) c=0 = **無 read 全覆蓋(2887)或全覆蓋卻只有全R(371)** | L1 實測 |
| 不需單 read 全覆蓋能否串樹？ | ✅ **能，已實作 = B_pairwise_structure 層**（A連B+B連C 傳遞式）**救回 2140 區** | L1 實測 |
| 3 個 must-fix（往過度自信方向） | **C3**（去噪丟 60.9% 仍標已定）/ **C2**（4-gamete 靜默丟）/ **D4**（甲基排序 overclaim 字串） | L1 實測 |

**一句話**：核心骨幹 sound；**c=0 與 spread 區靠 pairwise 傳遞式串接救回（B 層，非全覆蓋）**；唯二往「過度自信」方向的洞（C3/C2）須改效能前先補，其餘邊界誤差保守（與「定不出來即答案」同向）。

---

## §2 方法規格 — BAM + VCF → 演化樹（5 階段）

```
VCF(filtered_snv_tp/fp) + tumor.bam + normal.bam
 ├─[S0-S2] sm_linkage_genomewide.py → pairs + census
 ├─[ML]   sm_multilocus_combinations.py → populations(genotype 向量)
 └────────→ sm_region_integration.py → tree_shape + CN
                 → topology_analysis.py → determinacy + edges
                 → candidate_scoring.py → 佇列
```

| 階段 | 職責 | 核心判斷 | file:line |
|---|---|---|---|
| **S0 宇宙** | union TP∪FP per-chrom,TP 優先 | 每 sSNV 標 source | `sm_linkage:70-79` |
| **S1 分組** | 相鄰 gap>TIER_R(50kb)切 linkage group;每 sSNV 只 pileup 一次 | group-based 效率核心 | `sm_linkage:134-145,82-103` |
| **S2 共現** | 每對 co-read≥6 的 2×2(RR/RA/AR/AA)→ classify 6 類 | `aa<2`→mutual_excl/sparse;`ra==0&ar==0`→co_linked;`ar==0`→nested_a_in_b;`ra==0`→nested_b_in_a;else→independent | `sm_linkage:122-131` |
| **gate** | powered(coread≥6)+ normal=REF(is_somatic)+ same_hp | 三道閘才入 tally | `sm_linkage:40-44,219-223,235` |
| **ML** | 覆蓋全 sSNV 的 read → R/A genotype 向量;MINREAD≥3;MAX_SNV=8 | populations={向量:count} | `sm_multilocus:62-74` |
| **整合** | co_linked→union-find 併;nested→有向邊;mutual_excl→sibling;transitive_reduction 去冗餘+偵環 | tree_shape 7 類 | `sm_region_integration:93-150` |
| **拓樸** | population perfect-phylogeny laminar solve + 貪婪去噪 + determinacy | A_determined/B_pairwise/C/incompatible | `topology_analysis:63-109,163-168` |
| **評分** | confidence + resolution + 佇列 | needs_methyl 等 | `candidate_scoring:19-100` |

🔑 **非循環骨幹**：樹由 sSNV 單分子共現(genotype 軸)建,**不靠 HP tag、不靠甲基**;HP tag 只當鑑別器,normal=REF 確認 somatic。

---

## §3 c 群（n_clusters）的意義 + 為何有 0 群

**定義**（`topology_analysis:153-154`）：`n_clusters = len([g for g in populations if "A" in g])` = 該區觀察到的**帶突變的相異 genotype 向量數** = read-level 偵測到的相異突變組合(候選 subclone 節點)數。理論：觀察 c 群 → perfect-phylogeny 把拓樸限縮到 c 節點樹。

**c 群分布（全 7143 區,n_sSNV≥2）**：`{0:3258, 1:1995, 2:1744, 3:141, 4:3, 5:1, 6:1}`

**c=0（3258 區,46%）的兩個機制（實測）**：
| 機制 | 區數 | 意義 |
|---|--:|---|
| populations 空 | **2887** | 無 read 全覆蓋所有 sSNV(spread group;A2 全覆蓋約束)→ 拿不到聯合向量 |
| 只有全 R 群 | **371** | 有全覆蓋 read 但全 REF(低 VAF 突變未被全覆蓋 read 抓到) |

→ **c=0 = 「read-level 觀察不到任何帶突變的聯合組合」**,故 population 路徑畫不出節點。**但不代表無樹**——見 §4。

---

## §4 不需單 read 全覆蓋的串接 = pairwise 傳遞式（B 層,已實作）

**用戶問**：A連B、B連C、無 A連C 的 read,能否用 A-B + B-C 共同建樹?
**答：能,而且現行 `sm_region_integration.build_tree` 正是這樣做**（不需單 read 全覆蓋）：
- 每對(co-read≥6)給關係 → nested 建有向邊 → **A→B、B→C 傳遞式串成 A→B→C**,`transitive_reduction`(`:60-90`)去掉多餘 A→C。
- 這條路徑**只需成對共讀,不需單分子整跨**。

**實測救回量**（全 7143 區）：
| 層 | 定義 | 區數 |
|---|---|--:|
| **A_determined** | 有單分子全覆蓋向量（最強證據） | 3885 |
| **B_pairwise（c=0 但 pairwise 串出樹）** | **A連B+B連C 傳遞式,不需全覆蓋** | **2140** |
| C_underdetermined | c=0 且無 pairwise 結構 | 1118 |

B 層救回的 tree_shape：`linear_nested 760 / sibling_only 546 / full_tree 496 / co_linked_lineage 338`。
**範例**：`chr2:189984375-190158792`,11 sSNV / 174kb / **depth-4** linear_nested — 純靠 pairwise nested 邊串成 4 層樹,**無單一 read 跨全長**。

**🔴 盲點（為何 B 弱於 A,對齊稽核 C2/C5）**：pairwise 傳遞式可**推論出從未在單分子上直接觀察到的組合**。若 A-B、B-C 各自 pairwise-相容,但三方聯合 A-B-C 其實是 4-gamete 違反(只有 3-way 資料看得到),pairwise 串接會**拼出一棵 pairwise-fit 但實際不存在的 lineage**。故 B 層標為較弱證據;且此盲點被 C2（4-gamete 靜默丟）放大。

---

## §5 邊界正確性稽核（20 邊界 / 13 確認;3 must-fix 已獨立驗證）

### 確認合理（no-action,給安心）
A1(MAX_SNV cap 該保留:去 cap 會讓大 group 靜默 drop 非爆炸)、A2(全覆蓋→降B 是預期降級)、A4(向量長度無 bug)、B3(TIER_R=50kb 涵蓋 read 長尾)、C1(co_linked 併節點=identifiability 極限正解)、C4(父選擇)、D1(incompatible 是 dead code、非灌向 A;20萬壓測+真實3885區 0 觸發)。

### 需注意（確認為真問題,獨立驗證數字）
| ID | 問題 | 嚴重 | 影響正確性 | file:line |
|---|---|:--:|:--:|---|
| **C3** | 貪婪去噪無 eps floor、drop_frac 不 gate determinacy。**實測:39 區 drop>0.05、max 0.609、38/39 標 A_determined、0 incompatible、全 CN gain(22)/loh(17)** → 真 multiplicity 被當噪聲丟後仍高信心 | **HIGH** | 是 | `topology_analysis:70-79,159` |
| **C2** | `independent`(4-gamete)在 build_tree 靜默丟、不標 inconsistent。**實測:5198 對被丟、1137 區(15%)含違反、306 個 A_determined 內藏 4-gamete、63% min(RA,AR)≤1、多 CN-gain(774)** | **MED** | 是 | `sm_region_integration:108-117` |
| **D4** | resolution 字串「VAF tie→甲基輔助」違反核心結論『排序絕不用甲基』且 compute 中甲基缺席 = overclaim | LOW(口徑) | 否 | `candidate_scoring:36-40` |
| B2 | classify `aa<2` 絕對門檻:低深度真 nested(aa=1)誤判 sibling | MED | 是 | `sm_linkage:122-131` |
| B5 | same_hp 用 most_common 無純度門檻、None→diffHP;是 build_tree 硬過濾 | MED | 是 | `sm_linkage:219-223` |
| B1 | classify 精確==0 對雜訊脆弱(方向保守,不製造假 confident tree) | LOW | 是(aggregate 小偏) | `sm_linkage:122-131` |
| B4 | is_somatic NORM_REF_MIN=4 偏低(危險方向 bounded 且小) | LOW | 是(bounded) | `sm_linkage:40-44` |
| D1b | determinacy B 漏 single_nested/co_linked_lineage→other(分布低估) | LOW | 是(低估) | `topology_analysis:165` |
| D5 | incompatible resolution 一律「likely-artifact」,CN-gain multiplicity 誤標假影 | LOW | 是(處置字串) | `candidate_scoring:32,60-64` |
| C5 | 兩建樹路徑不 reconcile(真殘留=drop_frac 不 gate + 無 path-disagreement flag) | LOW | 是(<1%高估) | `topology_analysis:163-185` |
| A3 | MINREAD=3 對單-ALT 低複雜度假群敏感 | LOW | 是(n_pop敏感) | `sm_multilocus:74` |
| **D2** | A_determined 的 `sum(pops)>=6` 是 **no-op 死條件**(**實測 1812 A_determined、0 個 sum<6**,因 MINREAD=3×≥2向量恆真) | LOW | 否 | `topology_analysis:164` |
| D3 | dewey 姊妹編號用子樹 read 數當 CCF proxy,CN-gain 失真(純顯示) | LOW | 否 | `topology_analysis:111-134` |

---

## §6 改效能前必修的 3 處（正確性前提）

1. **C3(HIGH)**：加 eps noise floor（victim.count/total<0.02 才移除,否則回 incompatible）+ drop_frac>閾值降級 determinacy。現況把 40-61% multiplicity 群當噪聲丟後 collapse 成乾淨樹仍給 BASE=80。**最嚴重的正確性洞。**
2. **C2(MED)**：build_tree 計數 `n_independent`,4-gamete 經「噪聲容錯(排 min≤1)+ CN 分層(CN-gain 標 multiplicity 非錯誤)」後當品質旗標寫入 rec,非 CN-gain 的乾淨 4-gamete 觸發 inconsistent。
3. **D4(口徑)**：移除「VAF tie→甲基輔助」字串,改標 L3-weak / 甲基不作 resolver（純字串,不需重跑）。

> C3/C2 修正需重跑全基因組 pipeline(長 compute);D4 純字串低成本。

---

## §7 數字溯源表（§13-C,皆本輪 Read-back 獨立驗證）

| 數字 | 值 | 來源 |
|---|---|---|
| c 群分布 | {0:3258,1:1995,2:1744,...} | topology_per_region.json stats |
| c=0 機制 | 空 2887 / 全R 371 | sm_region_integration.json populations |
| A/B_rescue/C 層 | 3885 / 2140 / 1118 | 同上 × tree_shape |
| C3 drop | 39 區>0.05、max 0.609、38 A_det、0 incompat | topology_per_region.json drop_noise_frac |
| C2 4-gamete | 5198 對、1137 區(15%)、306 A_det、63% min≤1 | sm_linkage pairs rel==independent |
| D2 no-op | 1812 A_det、0 個 sum<6 | topology_per_region.json |

---

## §8 效能改進（workflow wf_8537f50f 已產計畫,待測）

效能審計主瓶頸定論 = **大 BAM 隨機讀,且 tumor BAM 被 pileup 兩次**（sm_linkage + sm_multilocus 對同批 sSNV 各一次全 pass）。前 3 優化：O1(合併雙 pass→單次 BAM 讀)、O3(normal pileup per-group 批次)、O4(by-chrom→by-group 平行吃滿 44 核)。**須測試 byte-identical 輸出後才採納**（正確性不變量:classify 6 類/populations/determinacy/byte-reproducible）。詳見效能 workflow 結果,下一階段實作+驗證。

---

## §9 REFLECTION

**警示指標**：C3/C2 未修前,「A_determined 乾淨定樹比例」被高估(306 A_det 內藏 4-gamete、38 A_det 丟>5% reads);論文若引用「已定樹比例」須先修或標 caveat。
**根因（double-loop）**：兩處洞都源於「衝突偵測不完整」——貪婪去噪(C3)掩蓋 population 層衝突、build_tree 漏接 independent(C2)掩蓋 pairwise 層衝突,兩層安全網同時漏 → incompatible 計數虛低為近 0。
**Reopen/待辦**：C3/C2/D4 修正 + 效能測試（§8）。
**Spaced recall**：2026-07-31。
