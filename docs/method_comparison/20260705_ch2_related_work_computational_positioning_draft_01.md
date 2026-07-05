---
title: Ch2 相關研究與計算定位（可投稿草稿）— 從 20260705 定位 doc 轉寫
date: 2026-07-05
status: draft
task_type: paper-prose-draft（步驟 3/3；由定位 doc + verified .bib 轉寫，非新分析）
build_branch: research/subclonal-reconstruction-202606
tier_of_claims: 形式化=L1 推理；外部 cite = 已過 /citation-verification（.bib）；內部數字=既有 ⭐3
data_sources: >
  docs/method_comparison/20260705_ism_computational_complexity_positioning_01.md,
  docs/method_comparison/_assets/ism_positioning_refs.bib
note: >
  \cite{key} 對應 _assets/ism_positioning_refs.bib；投稿前 .bib 中 author=others 的 7 條需補全作者。
  誠實邊界（原創綜合 / ⭐3 / necessary-not-sufficient / 不宣稱解析度優越）已內建，勿於編修時移除。
---

# 第二章（節選）相關研究與計算定位

> 撰寫紀律：本草稿把 `20260705_ism_computational_complexity_positioning_01.md` 的 bullet 定位轉為可投稿散文。所有引用鍵對應 verified `.bib`。誠實 hedge（原創綜合、⭐3、必要非充分、不宣稱勝過單細胞/多樣本）為載重口徑，編修時不可稀釋。

## 2.1 亞克隆重建的計算形式化：從完美系統發生到斯坦納樹

腫瘤亞克隆重建在計算上可被視為一個定義於**布爾超立方體** $\{0,1\}^n$ 上的樹重建問題：以 $n$ 個體細胞單核苷酸變異（somatic SNV, sSNV）為維度，每一個克隆對應一個二元基因型向量（帶有=1、未帶有=0），即超立方體的一個頂點；生殖系（germline）為全 0 的根頂點；在無限位點假設（infinite-sites assumption, ISA）下，突變只獲得一次且不丟失，因此演化邊為單調的 $0\to1$ 轉移，整棵樹為一個**有向、以全 0 為根的完美系統發生（perfect phylogeny）**。對二元字元而言，Gusfield 證明完美系統發生存在的充要條件為任兩字元皆相容——亦即不存在任一對位點同時出現四種配子組合 $\{00,01,10,11\}$（四配子檢定），且此判準可於線性時間內判定 \cite{gusfield1991efficient}。此「成對相容即全域相容」性質為二元字元所特有，對三態以上並不成立 \cite{gusfield1991efficient}。

當無限位點假設被放寬——允許復發（recurrence）、回復突變（back-mutation）或丟失（Dollo 型），樹上需要引入未觀測的中間基因型作為**斯坦納點（Steiner points）**，問題即退化為超立方體 Hamming 度量下的最小斯坦納樹，Foulds 與 Graham 證明此「系統發生斯坦納問題」為 NP-complete \cite{foulds1982steiner}。純理論層面，此類有向有根的超立方體斯坦納樹（Steiner arborescence）近年已有參數化（FPT）演算法 \cite{mahapatra2025parameterized}；而在近似層面，有向斯坦納樹有 $O(k^\epsilon)$ 多項式時間近似 \cite{charikar1999approximation}，其群組（group）變體有多對數近似 \cite{garg2000group} 與相應的不可近似下界 \cite{halperin2003polylog}。

值得強調的是，將**腫瘤**亞克隆重建明確 casting 為「布爾超立方體上的有向有根斯坦納樹」是本研究的形式化綜合：通用系統發生層面雖有超立方體斯坦納 arborescence 的處理 \cite{mahapatra2025parameterized}，但既有腫瘤系統發生文獻採用的是無限位點下的完美系統發生，並將解表述為祖先圖（ancestry graph）上滿足和條件（sum-condition）的**生成 arborescence（spanning arborescence）**——在 ISA 成立時並無隱藏節點，故非斯坦納樹 \cite{elkebir2015reconstruction}。因此本研究的定位為橋接 El-Kebir 的生成-arborescence 觀點 \cite{elkebir2015reconstruction} 與經典的超立方體-parsimony-斯坦納觀點 \cite{foulds1982steiner,mahapatra2025parameterized}，而非引用既有的具名結果。

## 2.2 複雜度分層與本方法的定位

本方法（ISM）在此形式化下的關鍵設計，是**只在多項式可解的完美系統發生特例中運算，並在證明為 NP-hard 的情形明確拒絕求解或列舉等可能解**，而非強擬一個最小成本斯坦納樹。具體而言：

**可解階段坐落於已被證明多項式可解的數學島。** 完整讀取下的骨幹建構即 Gusfield 的完美系統發生判定，多項式可解 \cite{gusfield1991efficient}。更關鍵地，本方法對**部分讀取（partial read）**的處理——長讀分子未涵蓋的位點被視為缺失項（"?"）而非強制填 0——恰好對應 Pe'er 等人所解的**不完全有向完美系統發生（Incomplete Directed Perfect Phylogeny, IDP）**問題：以圖夾擠（graph-sandwich）形式化，於多項式（$\tilde{O}(nm)$）時間內判定並補全一個有向、以全 0 為根、單調的完美系統發生 \cite{peer2004incomplete}（後續改良至線性時間）。此結果使本方法「補救約 40% 空基因型向量」的步驟具有嚴格的計算理論基礎，而非僅為工程啟發式。惟須誠實界定：IDP 保證的是**可行性（realizability）與一棵典範完成樹**（此為穩固結論），但在多個完成並存時**不唯一決定**具體填補——與多重性/CCF 的不可識別性一致 \cite{satas2021decifer}。

**被拒絕或列舉的階段，恰為已被證明 NP-hard 的數學島。** 在不相容（四配子違反）的區域，選擇丟棄最少位點以回復系統發生等價於最小頂點覆蓋，NP-hard \cite{day1986compatibility}；將混合（bulk）基因型列拆解為克隆列為最小無衝突列拆分，NP-hard 且不可近似 \cite{hujdurovic2018mixed}；以最少翻轉消除衝突為 min-flip，NP-hard \cite{chen2006minimumflip}；在超立方體上補全斯坦納結構為 NP-/MAX-SNP-hard \cite{foulds1982steiner,mahapatra2025parameterized}；允許甲基丟失（$1\to0$）則進入 Dollo-$k$，$k\ge2$ 時 NP-complete \cite{bonizzoni2016dollo1}。本方法在這些情形一律標記為疑似假象或列舉等可能候選，而不求最小成本解。由此，「單樣本 characterization（⭐3）＋定不出來即答案」並非保守的研究取捨，而是**對齊已證明複雜度分層的原則性設計**。

**繞過 NP-hard 的 VAF 反卷積。** 本方法以單分子讀取層級的共現普查**直接觀測**超立方體頂點；相對地，主流 bulk 方法（如 \cite{deshwar2015phylowgs,wintersinger2022pairtree}）從變異等位頻率（VAF）的邊際分布**推斷**頂點，而此推斷即 El-Kebir 證明為 NP-complete 的 VAF 分解問題（VAFFP），且通常需要多樣本方能識別 \cite{elkebir2015reconstruction}。因此，長讀單分子共現使本方法繞過了短讀方法必須求解或近似的反卷積難題——代價為超立方體維度上界（$n\le8$）、單樣本、以及對真正困難情形的拒解。

## 2.3 最近鄰方法與本方法的區別

在「以體細胞共現為遺傳骨幹、輔以單倍型與甲基」的設計軸上，最近鄰方法與本方法的差異可歸為四類。

**單分子體細胞共現（遺傳骨幹）。** 最近鄰為 Foltz 等人的工作，其以 10x 連鎖讀取（linked-read）的條碼共享判定兩個已定相的體細胞突變屬「同單倍型序列」抑或「分屬不同亞克隆」，並以 NRAS 的 G13R 與 Q61K 從不同時出現於同一分子推得兩者屬獨立亞克隆 \cite{foltz2024somatichaplotype}——此與本方法的互斥判定邏輯一致。然而其載體為連鎖讀取的液滴代理而非原生長讀單分子，且無甲基軸、無結構性檢定、不建樹。此案例亦同時提醒：共現→同譜系為**必要非充分**（互斥→分屬可靠，共現→同源仍需 VAF/CCF），本方法繼承此保留 \cite{foltz2024somatichaplotype}。本方法的骨幹工具為長讀體細胞單倍型標記，其體細胞路徑不使用任何甲基資訊，此為本方法非循環性的來源與保證 \cite{ho2025longphases}（惟屬同一實驗室產出，非獨立佐證）。

**VAF 推斷的系統發生學派。** PhyloWGS \cite{deshwar2015phylowgs} 與 Pairtree \cite{wintersinger2022pairtree} 等以 VAF 反推克隆樹，能重建至數十個亞克隆並隨樣本數增加而更準 \cite{wintersinger2022pairtree}；此類方法在規模上優於本方法，屬互補而非競爭。特別地，Pairtree 的「成對關係張量（pairs tensor）」在概念上對應本方法以 $O(n^2)$ 成對四配子普查加速骨幹的做法 \cite{wintersinger2022pairtree}。組合學派（AncesTree \cite{elkebir2015reconstruction}、SPRUCE \cite{elkebir2016spruce}、SPhyR \cite{elkebir2018sphyr}、PhISCS \cite{malikic2019phiscs}、CITUP \cite{malikic2015citup}）則以 ILP/枚舉求解，其中 SPRUCE 的「枚舉所有相容樹」\cite{elkebir2016spruce} 與 MACHINA/MACH2 對非唯一解空間的刻畫 \cite{elkebir2018machina,roddur2025mach2} 為本方法「定不出來即答案」的領域先例；SPhyR 的 $k$-Dollo 有界丟失 \cite{elkebir2018sphyr} 則為本方法可借鑑以處理有界復發的方向。

**甲基驅動的譜系重建。** 距離型甲基譜系樹（sgootr \cite{liu2023sgootr}）、單細胞甲基譜系（MethylTree \cite{chen2025methyltree}、Epiclomal \cite{desouza2020epiclomal}）與 CLL 的單細胞甲基-遺傳對應（Gaiti \cite{gaiti2019epigenetic}）皆以甲基為驅動訊號，且達到接近完全的譜系解析——但皆為**單細胞、且多為監督式**。此界定了本方法在 bulk 下僅能以甲基弱佐證的天花板：sgootr 已擁有「距離型甲基重建」此一方法類 \cite{liu2023sgootr}，故本方法不得宣稱該方法類為新穎，其甲基定位嚴格為有界輔助。讀取層級甲基分群（cvlr \cite{raineri2023cvlr}、ASMS \cite{raineri2024asms}、CpelNano \cite{abante2021cpelnano}）與本方法共用基質但主從相反，且皆具置換檢定，故本方法**不得以「對手缺顯著性檢定」或「對手非 ONT」作為差異點**（三者皆 ONT-capable 且具隨機化檢定）\cite{raineri2023cvlr,abante2021cpelnano}；監督式讀取來源分類（ROCIT 為 PacBio 長讀 \cite{baker2026rocit}、MethylBERT 之程式碼可讀 ONT \cite{jeong2025methylbert}）則與本方法的無監督結構檢定任務正交。甲基延長相鄰定相區塊之工作（MethPhaser）明確將癌症留為未來工作 \cite{fu2024methphaser}，此正證本方法的白地為真，並印證「以甲基建單倍型將導致循環」的設計顧慮。

**同生態的 bulk-ONT 多組學。** TumorLens 與本方法共享 bulk-ONT、配對正常樣本與等位/甲基感知的生態，但其為偵測/註記管線而非重建器，甲基以 10 kb 窗差異甲基區（DMR）呈現而非讀取層級結構，且無體細胞順式檢定 \cite{paulin2026tumorlens}；此為本方法三軸（無監督讀取結構檢定、正常錨定體細胞順式、體細胞亞克隆目標）不重疊之處。

## 2.4 定位口徑：可宣稱與不可宣稱

綜上，本方法可辯護的貢獻為三點：其一，**實用性與成本**——自單一 bulk ONT 定序即取得亞克隆相關的讀取層級訊號，無需單細胞製備、多樣本 panel 或連鎖讀取文庫，而上述每個更高解析度的鄰近方法皆需其一；其二，**同分子原生甲基**——甲基與 sSNV 讀自同一條物理長讀分子，無需另一實驗或條碼/插補拼接；其三，**非循環設計**——單倍型標記源自零甲基的體細胞定相 \cite{ho2025longphases}，甲基嚴格為標籤優先（label-first）之有界輔助而從不進入似然，正符 MethPhaser 對癌症甲基定相循環風險的告誡 \cite{fu2024methphaser}。

反之，本方法**不主張**下列二者：其一，解析度或確認的優越性——單細胞甲基方法 \cite{chen2025methyltree,liu2023sgootr,gaiti2019epigenetic} 與多樣本 VAF 方法 \cite{wintersinger2022pairtree} 在重建上均更強，本方法定位為 ⭐3 characterization、僅輸出集成而非單一確認樹；其二，單分子體細胞共現為本方法所獨有——Foltz 以連鎖讀取先行發表 \cite{foltz2024somatichaplotype}，且不可僅憑共現宣稱已確認亞克隆（必要非充分，僅互斥→分屬方向可靠）。

---

## 待辦（投稿前）
1. `.bib` 中 author 含 "others" 的 7 條補全完整作者（Foltz/Wakhan/Severus/MethPhaser/MethylTree/ROCIT/Epiclomal）+ ASMS 全作者。
2. 全文對照定位 doc §5「絕不可宣稱」清單，確保編修未引入 overclaim。
3. 依實際投稿格式將 \cite{} 轉為對應 citation 樣式。
