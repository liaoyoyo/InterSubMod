---
title: ISM 計算複雜度定位 — hypercube-Steiner 形式化 + 數學島 + 最近鄰優勢（論文 Ch2/Discussion）
date: 2026-07-05
status: in_progress
task_type: positioning-note（理論定位，非新實驗；合併本 session 3 個 workflow/agent 的對抗驗證產出）
build_branch: research/subclonal-reconstruction-202606
tier_of_claims: >
  形式化裁決 = L1（推理）；外部 citation = L3（本 session web/一手核，投稿前一律過 /citation-verification）；
  內部數字 = 既有 ⭐3 prior work（已 grep 回真值）。
data_sources: >
  本 session workflows: wf_8655ec3f-69e(steiner-formalization) / wf_33fc5661-204(nearest-neighbor) / agent(math-solved);
  docs/method_comparison/20260630_ism_positioning_vs_prior_work_01.md,
  docs/method_comparison/20260703_cn_axis_perread_feasibility_and_sSNV_methyl_role_01.md,
  docs/methodology/20260628_subclone_reconstruction_master_spec_01.md
related:
  - docs/method_comparison/20260630_ism_positioning_vs_prior_work_01.md
  - docs/method_comparison/20260703_cn_axis_perread_feasibility_and_sSNV_methyl_role_01.md
  - docs/methodology/20260621_baseline_discipline_conceptual_definitions_01.md
---

# ISM 計算複雜度定位：hypercube-Steiner 形式化 + 數學島 + 最近鄰優勢

> **一句話**：ISM 的可解階段（骨幹 + partial-read 補全）踩在**已被證明多項式可解**的數學島（Gusfield perfect-phylogeny + Pe'er 2004 Incomplete Directed Perfect Phylogeny）；它拒絕/enumerate 的階段全部是**已被證明 NP-hard** 的島。→「⭐3 characterization + 定不出來即答案」= **對齊真實複雜度分界的原則性設計**，非保守選擇。
>
> **🔴 原創性界線**：把**腫瘤 subclone 重建** casting 成「布爾超立方體上有向有根 Steiner tree」是**本研究的原創綜合**，非既有 named result（通用系統發生版有 Mahapatra 2025，腫瘤版無人做）。

---

## §0 目的、範圍、證據紀律

**目的**：給 ISM 一個計算複雜度層的理論定位，供論文 Ch2（Related Work）+ Discussion 直接引用，並作為 reviewer 攻防口徑。
**範圍**：形式化 + 數學已解問題 + 最近鄰競品優勢——**非新實驗**；合併本 session 三個對抗驗證過的 workflow/agent 產出。
**證據紀律**：
- **所有外部 citation = L3**（本 session WebSearch/WebFetch 或一手 CONTEXT 卡核；投稿前一律過 `/citation-verification`，見 §5）。
- 內部數字（partial-read 40.4%、k_ISM⊥k_CN 等）= 既有 ⭐3 work，已 grep 回真值（源見 related docs）。
- **誠實標註**：同實驗室（LongPhase-S）、attribution UNCERTAIN（min-flip「Chen-Kolmogorov」查無 Kolmogorov）、arXiv-only（Dollo-1 poly）、原創綜合（tumor-hypercube-Steiner）— 均在對應處標。

---

## §1 形式化：布爾超立方體上的有向有根（group）Steiner tree

### §1.1 物件映射

| 形式化元素 | ISM 對應 | 源碼 |
|---|---|---|
| 超立方體頂點 {0,1}ⁿ | n 個 somatic sSNV 的二元基因型向量（R=0/A=1，n≤8=`MAX_SNV`） | `sm_multilocus_combinations.py:70` |
| 根 = 0ⁿ | germline all-REF（RR），parsimony 生根 | `topology_analysis.py:305` |
| 有向邊 0→1（單調） | 突變獲得，infinite-sites 下遞增 Hamming weight | `topology_analysis.py:126`（父=maximal proper subset） |
| four-gamete test | perfect phylogeny 相容性判準（Gusfield）：linear=nested / branched=mutual-excl / incompatible=四配子全滿 | `sm_linkage_genomewide.py:132-141` |
| Steiner 節點 | 未觀測的祖先 clone（`H_*` 隱藏節點，只在 perfect phylogeny 可證存在處插入） | `topology_analysis.py:152-207` |

### §1.2 三層誠實裁決（對抗驗證後）

| 層 | 宣稱 | 裁決 | 理由 |
|---|---|---|---|
| ①**hypercube + directed + rooted-at-RR** | ✅ **準確** | max-parsimony over binary ≡ Hamming hypercube 樹；infinite-sites 假設 explicit，違反是 first-class 輸出（`recurrence_required`），非默默 forced |
| ②**"Steiner tree"** | ⚠ **物件準/計算不準** | H_* = 真 Steiner 點（物件對）；但 ISM **只算多項式 perfect-phylogeny 特例**，**不解** NP-hard Steiner 最佳化（無 min-cost 目標、無 branch-length 搜尋、無 enumerate-and-score） |
| ③**"GROUP" Steiner** | 🔴 **過度宣稱(as solved)** | 有真結構類比（partial-read `X`→subcube 群，程式註解literally「供 group-Steiner 覆蓋條件」；multiplicity 不確定集）但 ISM 不做 cover-and-minimize；**白地**（無腫瘤 phylo 用 group-Steiner）→ 當**原創 motivation/future-work** |

🔴 **casting 腫瘤重建成 hypercube-Steiner = 原創綜合**：通用系統發生版有 Mahapatra 2025（MSA-DH；*Acta Inf* 62(1):6，arXiv:2110.02830）；腫瘤文獻用 perfect phylogeny（ISA）+ ancestry-DAG 的 **spanning arborescence**（El-Kebir 2015，ISA 下無 latent 節點故非 Steiner）。→ 論文寫成「**橋接 El-Kebir spanning-arborescence 與經典 parsimony-Steiner-on-hypercube 的新形式化**」，不可寫成引用既有結果。

### §1.3 一句話可辯護形式化（直接進論文）

> 「每個區域內，ISM 以多項式 four-gamete/pairwise-compatibility 檢定（Gusfield 1991）重建 somatic-SNV 基因型向量的**有向有根 perfect phylogeny**——布爾超立方體 {0,1}ⁿ 頂點、根在 germline all-REF 角——只在 perfect phylogeny 可證存在處插入未觀測祖先基因型作為 Steiner 點，並在 infinite-sites 被違反處拒絕或枚舉等可能樹。」

**明確不宣稱**：不解 NP-hard (group-)Steiner 最佳化；只算多項式 perfect-phylogeny 特例，對 incompatible/multiplicity 不確定情形 flag/枚舉，不強擬最小成本樹。

---

## §2 數學島：誰解過類似問題 + ISM 各階段坐落

### §2.1 數學島全圖

| ISM 階段 | 數學問題 | 複雜度 | 已解? | Citation（L3） |
|---|---|---|---|---|
| 骨幹建構（完整 read） | directed perfect phylogeny / four-gamete | poly | ✅ SOLVED | Gusfield 1991 *Networks* 21(1):19–28, DOI 10.1002/net.3230210104 |
| 🎯**partial-read 補全（X=缺失）** | **Incomplete Directed Perfect Phylogeny (IDP)** | poly→linear | ✅ **SOLVED** | **Pe'er, Pupko, Shamir, Sharan 2004** *SIAM J Comput* 33(3):590–607, DOI 10.1137/S0097539702406510；linear: SPIRE 2021 DOI 10.1007/978-3-030-83508-8_13 |
| 有向-root vs 無向 | directed 版更易 | poly | ✅ | 同上（有向 rooted = 可解分支） |
| incompatible 選丟哪些 loci | Min Character Removal ≡ **Vertex Cover** | NP-hard, FPT-in-k | flag/enumerate | Day & Sankoff 1986 *Syst Zool* 35:224–229, DOI 10.2307/2413304 |
| 混合 bulk row 拆成 clone | Min Conflict-Free Row Split | NP-hard + **inapproximable**, FPT | flag | Hujdurović et al. 2018 *IEEE/ACM TCBB* 15(1):96–108, DOI 10.1109/TCBB.2016.2606620 |
| 錯誤翻轉去衝突 | Min-flip | NP-hard, FPT | flag | Chen et al. 2006 *IEEE/ACM TCBB* 3(2):165–173, DOI 10.1109/TCBB.2006.26 ⚠「Chen-Kolmogorov」查無 Kolmogorov |
| Steiner 補全 hypercube | Steiner in {0,1}ⁿ | NP-/MAX-SNP-hard, FPT-only | 不解 | Foulds & Graham 1982 *Adv Appl Math* 3:43–49；FPT: Mahapatra et al. 2025 *Acta Inf* 62(1):6, DOI 10.1007/s00236-024-00474-8 |
| 允許甲基 loss（1→0） | Dollo-k | **k=1 poly / k≥2 NP-complete** | 不入（保持單調） | Bonizzoni et al.（Dollo-1 arXiv:1611.01017=L3；constrained PPP PMC4240218） |
| multiplicity/CCF 唯一性 | non-identifiable | **impossibility（證明）** | 定不出來 | DeCiFer（Satas 2021）*Cell Syst* 12(10):1004, DOI 10.1016/j.cels.2021.07.006 |

### §2.2 🎯 Pe'er 2004 IDP = partial-read 補全的嚴格基礎

「**有缺失資料的有向 perfect phylogeny 多項式（→線性）可解**」正是 ISM partial-read（`X`=未覆蓋位點=「?」entry）的數學對應。→ **「救 40.4% 空向量」不是工程 trick，而是有嚴格 Incomplete Directed Perfect Phylogeny 基礎的補全**，前提：單調獲得 + all-zeros 生根、X 當「?」非強制 0。

**🔴 caveat（誠實）**：Pe'er 2004 保證**可行性（realizability）+ 一棵 canonical 完成樹**（L1-solid），但**不唯一決定**多完成並存的 entry——那部分 heuristic/non-identifiable，與 DeCiFer multiplicity 不可識別一致。
→ **宣稱口徑**：partial-read 補全的**可行性/一致性**有嚴格數學基礎；**選哪個完成**則 bounded/非唯一。

### §2.3 複雜度邊界正當化（最強設計論述）

> **ISM 的演算法邊界 = 真實計算複雜度邊界**：多項式可解者（perfect phylogeny + IDP 補全）就算出來、確定；NP-hard 者（挑 loci ≡ Vertex Cover、拆 row inapproximable、min-flip、Steiner 補全、Dollo-k≥2）就 flag/enumerate。→「⭐3 characterization + 定不出來即答案」= **對齊已證明複雜度分層的原則性設計**，reviewer 難攻擊。

---

## §3 限制與加速：別人怎麼做 + ISM 借用 lever

### §3.1 限制（constrain — 生物先驗把搜尋縮到多項式）

| 限制 | 代表 | ISM 現況 | 可借? |
|---|---|---|---|
| ISA → perfect phylogeny | Gusfield；El-Kebir 2015（AncesTree）；Canopy | ✅ 已用（four-gamete） | 核心 |
| Sum-condition / VAF pigeonhole（父 CCF≥Σ子） | AncesTree, CITUP, Pairtree | ⚠ 只用在 sibling 排序 | ✅ 可借剪枝 |
| k-Dollo（有界丟失/復發） | SPhyR, SASC | ❌ 全 flag incompatible/recurrence | ✅ 可借處理有界 recurrence（70 recurrence_required 或救一部分） |
| Multi-state CN / multiplicity | SPRUCE, DeCiFer | ⚠ CN 走獨立 m-channel | ⚠ 借需防 CN confound |
| Multi-sample / longitudinal | CITUP, CALDER | ❌ 單樣本 | future work（serial ctDNA） |

### §3.2 加速（accelerate — 演算法）

| 加速 | 代表 | ISM 現況 | 對應 |
|---|---|---|---|
| **Pairwise 關係/tensor 預算 O(n²)** | Pairtree pairs tensor | ✅ **已用（O2 four-gamete 普查）** | **直接對應** |
| 區域分解 divide-and-conquer | (ISM 特色) | ✅ 7,143 區 | ISM 原生 |
| ILP/CSP 精確解 | PhISCS, El-Kebir | ❌ | 解 incompatible 才需 |
| MCMC 後驗抽樣 | PhyloWGS, Canopy, Pairtree | ❌（改 enumerate） | 替代路線 |
| Beam search | Orchard | ❌ | ✅ 可借做候選枚舉 |
| 近似演算法 | Charikar 1999(directed); GKR 2000(group) | ❌（停多項式島） | 一般 Steiner 才需 |

### §3.3 🎯 ISM 繞過 NP-hard VAF-factorization（觀測 vs 推斷）

ISM 用**單分子 read-level 共現普查**=**直接觀測**超立方體頂點；bulk VAF 方法（PhyloWGS/Pairtree）從 marginal **推斷**頂點——那個推斷 = **VAF factorization = NP-hard（El-Kebir VAFFP 證）**，且需**多樣本**才可識別。
→ **ISM 用 ONT 單分子共現「繞過」了 bulk short-read 必須解的 NP-hard 反卷積**：把「哪些突變共存」直接讀出來，而非從邊際頻率反推。這是「單分子共現=命脈」的複雜度版論述（verified、可辯護）。**代價**：n≤8 cap、單樣本、hard case 拒解 → ⭐3。

---

## §4 最近鄰競品：運作 / 限制 / 結果 / 我方優勢

### §4.1 最近鄰 4 名 + 唯一關鍵差

| 名次 | 最近鄰 | 為何最近 | 唯一關鍵差 vs ISM |
|---|---|---|---|
| 1 | **Foltz 2024**（SomaticHaplotype） | 同「somatic 共現→subclone」邏輯+HP；就是 chr2:18M 互斥招 | 載體=10X linked-read barcode（非原生 ONT 單分子）；零甲基；只計數無結構檢定/樹 |
| 2 | **sgootr 2023** | distance-based 甲基譜系樹=read×read 方法類 | 單細胞；甲基是 DRIVER；無 normal-cis/PERMANOVA/somatic-haplotag |
| 3 | **TumorLens 2026** | 同 bulk-ONT + matched-normal + allele/甲基生態 | 偵測/註記器非重建器；甲基=10kb 窗 DMR 非 read-level |
| 4 | **LongPhase-S 2025** | 是 ISM 自己的骨幹工具 | 做標記非重建；⚠ 同實驗室（非獨立佐證） |

### §4.2 誠實優勢表（skeptic 模式）

| 最近鄰 | 他們做什麼 | 限制 | ISM 優勢 — 或「無優勢」 |
|---|---|---|---|
| **Foltz 2024**（PMID 39149342） | 單分子代理 somatic 共現→subclone；HP；NRAS G13R+Q61K 從不同時 ALT→獨立 subclone | linked-read 非 ONT；零甲基；無結構檢定；MM-only；無 cis | 🔴 **共現本身非獨有**（Foltz 有）。edge=原生 ONT 單分子 + 同分子甲基 + normal-cis + PERMANOVA。差異化打**軸**非「觀測到共現」 |
| **sgootr / MethylTree / Gaiti** | 甲基譜系樹 ~100%（MethylTree Q≈1；Gaiti SF3B1 Fisher 7.4e-9） | 單細胞 + 監督式 | 🔴 **他們 lineage 確認比 ISM 好 — 老實講**。ISM 無解析度/確認優勢；edge=bulk-ONT 實用性。**絕不宣稱「distance-methylation reconstruction」為新**（sgootr 擁有此方法類） |
| **TumorLens**（PMID 41891012） | bulk-ONT 多組學偵測 | 非重建；窗式 DMR；無 cis；無 read×read | ✅ **真優勢**：無監督 subclone 結構檢定 + read-level 甲基 + somatic-cis；三軸不重疊 |
| **LongPhase-S**（bioRxiv 2025.11.20.689492） | somatic haplotag（零甲基） | 只標記，不重建；⚠ 同實驗室 | ✅ 在其上加重建層；是依賴非對手 |
| **cvlr / ASMS**（PMC9887406 / bioRxiv 2024.12.18.629129） | read×CpG 甲基分群 + permutation | germline/imprinting；無 normal-cis；測甲基差非距離結構 | ⚠ 窄但真=normal-cis + distance pseudo-F。🔴 **絕不說「他們沒顯著檢定」**（都有 permutation）——打統計量+somatic 目標 |
| **Pairtree / PhyloWGS**（PMID 36129821 / 25786235） | VAF 推斷樹 ~30 subclone，多樣本 | 零 read 共現；需多樣本；無甲基/HP | ✅ **互補非對手**（大規模重建贏 ISM）。edge=觀測 vs 推斷、繞 NP-hard；**不宣稱勝過** |

### §4.3 三可宣稱 + 兩絕不可

**✅ 可宣稱**：
1. **實用性/成本**：單一 bulk ONT run 取 read-level subclone 訊號——免單細胞製備/多樣本 panel/linked-read library（每個更高解析度鄰居都需其一）。
2. **原生同分子甲基**：甲基與 sSNV 讀自同一條物理 ONT 分子——無需另一 assay、無 barcode/imputation 拼接。
3. **非循環設計**：HP 來自零甲基 LongPhase-S；甲基嚴格 label-first/bounded-auxiliary，永不進 likelihood（MethPhaser cancer-循環警告要求的防火牆）。

**🔴 絕不可宣稱**：
1. **解析度/確認優越**——單細胞法（MethylTree/sgootr Q≈1、Gaiti 7.4e-9）+ 多樣本 Pairtree 都重建更好；ISM=⭐3 characterization、ensemble-only。
2. **單分子 somatic 共現為獨有**（Foltz 已發表）；且**不可從共現單獨宣稱已確認 subclone**（necessary-not-sufficient；只有互斥→分開可靠）。

---

## §5 citation queue（投稿前 `/citation-verification`）+ 誠實標註

> ✅ **2026-07-05 全 35 cite 已過 `/citation-verification`（workflow wf_630cf93a-ed8，WebSearch+Scholar 4-field）→ verified `.bib`：`InterSubMod/docs/method_comparison/_assets/ism_positioning_refs.bib`（0 CITATION_NEEDED；修正 log 見 .bib 檔尾）。以下標註已更新為 verified 後狀態。**

**理論骨幹（metadata 已核）**：Gusfield 1991 · Foulds-Graham 1982 · Charikar 1999 · Garg-Konjevod-Ravi 2000 · Halperin-Krauthgamer 2003 · **Mahapatra 2025**（*Acta Inf* 62(1):6，非 2024/vol61）· El-Kebir 2015 · **Pe'er 2004 IDP**（partial-read 基礎，poly-time 逐字確認） · Day-Sankoff 1986 · Hujdurović 2018 · **Chen 2006**（⚠ Kolmogorov 確認非作者） · **Bonizzoni Dollo-1**（arXiv-only；作者 Soto Gomez/Trucco 非 Dondi） · DeCiFer 2021。

**最近鄰**：Foltz 2024 · sgootr 2023 · TumorLens 2026 · LongPhase-S 2025(⚠ 同實驗室) · cvlr 2023 · ASMS 2024 · CpelNano 2021 · Gaiti 2019 · MethylTree 2025 · Epiclomal 2020 · MethPhaser 2024(cancer deferred=白地實證) · ROCIT 2026(🔴 PacBio 非 ONT) · MethylBERT 2025(🔴 code 吃 ONT，勿以「ONT vs bisulfite」差異化)。

**🔴 已解決/剩餘**：min-flip「Chen-Kolmogorov」歸屬 → **已核 Kolmogorov 非作者**（Chen/Eulenstein/Fernández-Baca/Sanderson）；PhISCS PMID 31628257→**31628256**（原為 SiCloneFit）；cvlr 描述「call variants」**係錯誤**（實為 methylation read clustering）；Mahapatra 2025 後無其他 hypercube-Steiner exact/parameterized 結果（UNCERTAIN=none）；Dollo-1 poly 目前 arXiv-only。**剩餘投稿前工作**：`.bib` 中 author 含 "others" 的 7 條補全完整作者名單（Foltz/Wakhan/Severus/MethPhaser/MethylTree/ROCIT/Epiclomal）+ ASMS 全作者（bioRxiv 403）。

---

## 關聯文件

- ISM 定位主文：`InterSubMod/docs/method_comparison/20260630_ism_positioning_vs_prior_work_01.md`
- CN 軸 per-read + 角色分工：`InterSubMod/docs/method_comparison/20260703_cn_axis_perread_feasibility_and_sSNV_methyl_role_01.md`
- baseline 詞典（confirm vs describe、cis-ASM vs lineage）：`InterSubMod/docs/methodology/20260621_baseline_discipline_conceptual_definitions_01.md`
- 重建 master spec：`InterSubMod/docs/methodology/20260628_subclone_reconstruction_master_spec_01.md`
