<!--
建立時間: 2026-03-01 00:00
目標: InterSubMod 專案的研究主軸、方法論與可突破方向的深度分析報告
處理範圍: 長讀序腫瘤/正常資料中體細胞變異與讀序甲基化的連結方法、顯著性分析、LOH/CNV 品質分層
關聯檔案:
  - docs/README.md
  - docs/reports/research_landscape/00_INDEX.md
-->

# InterSubMod：以長讀序腫瘤/正常資料連結體細胞變異與讀序甲基化的研究主軸、方法與可突破方向

## 問題重新定義與研究定位

你想做的核心其實是把「同一條長讀序（single molecule/read）」上同時存在的三種訊號串起來：體細胞突變（somatic SNV/INDEL/SV）、相位/單倍型（phasing/haplotype；HP/PS 或等效資訊）、以及 CpG 甲基化（methylation）圖樣，並進一步回答三類更難的問題：一是突變真偽與技術偽陽性（TP/FP）辨識；二是等位/單倍型特異性甲基化（ASM/allele-specific methylation）與「第二擊（second hit）」/雙等位失活（biallelic inactivation）的機制整合；三是腫瘤內亞克隆（subclone）在表觀遺傳層面的結構與演化順序。InterSubMod 在 repo 的系統總覽文件中把這件事定義得很清楚：以 ONT 長讀序（含 MM/ML 的甲基化標籤、HP/PS 的相位標籤）為前提，針對每個 somatic SNV 形成一個視窗，建立 read×CpG 甲基化矩陣、計算 read–read 距離、做階層式聚類，然後檢定聚類與標籤（ALT/REF、HP1/HP2、Tumor/Normal）的關聯性，作為「變異驗證」「ASM」「亞克隆結構」的基礎。fileciteturn34file0L20-L35

此定位與近年大型癌症長讀序資源型研究高度相符：例如 Cell Genomics 的 Long-Read POG（advanced cancer cohort）明確展示「長距離相位可用來確認腫瘤抑制基因雙等位失活、找 aDMR（allelically differentially methylated regions）、與 allele-specific expression」，並指出像 BRCA1、RAD51C 等修復基因的啟動子甲基化，可能在缺乏 coding driver 時成為 HRD（同源重組缺陷）的驅動因素，且常需搭配第二擊（例如 LOH）才形成臨床相關表型。citeturn13view0turn13view3turn12view0

## InterSubMod 目前的研究方法與系統流程

InterSubMod 的做法可拆成「以變異為錨點（SNV-anchored）」「以單分子甲基化圖樣作為距離/聚類依據」「用標籤關聯檢定把聚類結果轉成可解釋訊號」三條主線。repo 的專案摘要與架構文件，把資料流與假設寫得相當具體：只處理 biallelic SNV（single ALT），每個 SNV 建立 ±window 的 ROI，從 BAM 擷取 reads，解析 MM/ML 得到 CpG 機率或二值化狀態，形成 read×CpG 矩陣，並可對 PMD/異染色質等高噪訊區域做 gating（排除或降權），最後計算距離矩陣與做聚類、輸出熱圖/樹/summary。fileciteturn5file0L12-L45

在「甲基化訊號抽取」層，InterSubMod 直接吃 ONT BAM 的 MM/ML 標籤：以 CIGAR 建立 read-seq index ↔ reference coordinate 的對應，解碼 MM 的 delta-encoding 位置，依 ML 取出對應機率並確認 CpG context（含 strand/反向互補處理）。這部分的實作在 `MethylationParser.cpp` 中可看到對 MM/ML、CIGAR、正反股 CpG 的細節處理，是整個系統正確性的地基。fileciteturn6file0L1-L8

在「矩陣與距離」層，InterSubMod 同時保留 raw 機率矩陣與二值矩陣（missing 以 NaN 或 -1 表示），並提供多種距離度量：NHD（normalized Hamming，基於二值矩陣）、L1/L2（基於 raw 機率）、相關距離、Jaccard，以及特別針對 ONT 機率不確定性設計的 Bernoulli 距離（用「期望不一致」與「信心加權」降低 p≈0.5 的雜訊位點支配距離）。fileciteturn29file0L1-L7

在「聚類」層，InterSubMod 以階層式聚類（UPGMA、Ward、single/complete linkage）建立樹並支援切樹得到 cluster labels；文件中也提到可用 silhouette 找最佳群數，並輸出 Newick 等格式。這些能力與你的目標吻合，因為「先形成可重現的甲基化結構，再討論與突變/單倍型的對應」能避免一開始就把模型綁死在某個生物假設上。fileciteturn41file0L60-L90

## 顯著性分析、雙向驗證與 LOH/CNV 品質分層

InterSubMod 最具辨識度的地方在於它不是只畫熱圖，而是把「聚類→關聯檢定→結構檢定→再回頭用標籤驗證」做成一套可量化的多層顯著性架構。`SignificanceAnalyzer` 把分析拆成多個 phase：先做全域列聯檢定（Fisher–Freeman–Halton；同時給 Cramér’s V 效應量），通過 gate 才做 local tests 與結構檢定（PERMANOVA、dispersion），最後再做「Label-First」多階段驗證（直接依 HP 或 allele 分組，看距離矩陣是否能分出組間差異），並將 cluster-first 與 label-first 的一致性映射成 `VerificationClass`（Strong/Subclone/Weak/Noise）。fileciteturn33file0L67-L110

更貼近你問的 LOH/CNV/second hit：InterSubMod 在 2026-01 的多層驗證實作報告中，新增了 HP ratio（HP1/(HP1+HP2)；含 Laplace smoothing）、Potential LOH（<0.1 或 >0.9）、coverage multiple（NumReads/75 的近似倍數）、coverage category（CNV_Loss/CNV_Gain 等）、以及 LOH subtype（LOH_Noise/Weak/Strong/Subclone）與 quality score（0–100）分層。這些設計本質上是在承認：在只有 tumor/normal（甚至 tumor-only）長讀序的情境下，LOH 與 copy-number 不只是「背景」而是會直接扭曲你對甲基化–突變關聯的判讀，因此必須先把「可疑 LOH、coverage 異常、CpG/reads 太少」標成低可信度。fileciteturn39file0L12-L33

repo 內部的全面分析報告也提供了很重要的現況診斷：InterSubMod 嘗試用 TP/FP（以 HCC1395 類資料的既有結果）比較，發現僅用原本的 Cramér’s V、p-value 等指標，TP/FP 的 AUC 接近 0.5，鑑別力不足；之後加入 LOH/coverage/quality score 後，quality score 的 AUC 提升到約 0.62，但單獨用作過濾器仍難以超越既有 baseline（例如 VCF QUAL 的策略），更適合當「風險分層與診斷輔助」。這個結論非常關鍵：它暗示 InterSubMod 的突破點，可能不是再加一個更複雜的閾值，而是要改變「特徵/證據整合方式」或「選擇哪些 CpG/區域才值得讓甲基化訊號參與決策」。fileciteturn38file0L266-L353 fileciteturn39file0L118-L176

## 與核心外部研究的對照：哪些方法可直接借鑑

Cell Genomics 的 Long-Read POG 研究非常像你想做的「臨床級第二擊/相位/甲基化整合」範例：他們用長讀序相位去確認腫瘤抑制基因雙等位失活，並指出多數雙突變是在 trans、也會出現 cis；更重要的是，他們把 ASE 的主要驅動因素指向 CNA/LOH（ASE genes 傾向落在 LOH 與 copy number imbalance 區域），同時也觀察到在 copy number 平衡區域，promoter 的 allelic methylation 反而更常見，且常在 trans 對應 major expressed allele（暗示被甲基化的等位基因表現較低）。對 InterSubMod 來說，這是兩個很實務的提醒：其一，「甲基化與表現/突變的關聯」必須顯式引入 CNA/LOH 作為解釋變項；其二，「allelic promoter methylation」本身可成為因果候選（對應表現下降），但你需要能在等位層級做判讀，而不是只看 bulk methylation rate。citeturn13view0turn13view2turn12view0

Long-Read POG 在方法面也提供可複製的工程路線：他們以 MBASED 做 ASE 的 beta-binomial 檢定與 haplotype 聚合，並在後處理階段把 CNV（Ploidetect）、allelic methylation（NanoMethPhase）與 somatic calls 一起丟進 IMPALA 之類的整合流程，用 bedtools intersect 等把不同證據對齊到同一基因/區段，產生可解釋 summary。InterSubMod 目前是以「每個 SNV 視窗」為單位；POG 的思路則更像「以基因/機制（TSG inactivation, HRD, promoter methylation）」為單位，把多種訊號融合。這提供一個很具體的突破方向：InterSubMod 可以保留 per-SNV 的細粒度，但在上層加一個「gene-level/segment-level evidence integrator」，讓 second hit 的推論不必被侷限在單一 SNV 視窗是否顯著。citeturn13view2turn13view3

Nature Genetics（TRACERx NSCLC）那篇（你特別點名要深讀）則直接命中 InterSubMod 的 LOH/CNV 痛點：作者指出 bulk 腫瘤樣本的 methylation rate 會被 purity 與 CN instability 混淆，因此用 CAMDAC（copy number-aware methylation deconvolution）把純腫瘤細胞的 methylation rate 從 bulk 中解卷積出來；他們甚至把同樣的理念延伸到 PDR（proportion of discordant reads；一種甲基化「隨機性/雜訊」指標），並用「LOH 區域的 SNV 當作天然標籤」來驗證 deconvolution 後的 PDR：在 LOH 區域，帶 SNV 的 reads 幾乎可視為腫瘤來源，wild-type reads 近似正常污染來源，因此可以用 SNV 分組估計 PDR，並觀察到 CAMDAC 的 PDR 與 SNV-based PDR 有顯著高相關（R>0.8）。這對 InterSubMod 幾乎是「可以直接搬來用」的驗證框架：你現在已經在做 per-SNV read 標記（ALT/REF/UNKNOWN），只差把這個標記從「關聯檢定」升級成「用來估計與校正腫瘤/正常甲基化背景」的工具。citeturn20view0turn20view2

同一篇 Nature Genetics 也提出兩個更上層的概念，能幫 InterSubMod 重新定義「什麼叫甲基化驅動」：其一是 intratumoral methylation distance（量化腫瘤內甲基化異質性），其二是把 promoter 的 CpG 分成 regulatory 與 nonregulatory，建立類似 dN/dS 的 \(M_R/M_N\) 指標，用來區分「偏向功能性調控的高甲基化」（可能受正向選汰）與「偏向非調控位點的高甲基化」（較像乘客事件）。InterSubMod 若要把「甲基化現象」推到「驅動突變/驅動機制」層級，可能必須從這類「選汰/功能關聯」的框架汲取想法，而不只是依賴 cluster–label 的顯著性。citeturn0view2turn20view2

Nature 的 EVOFLUx（fluctuating CpGs）又提供另一個互補觀點：某些 CpG 位點的甲基化會以近似獨立、持續波動的方式變動，形成可追蹤系譜的「演化條碼（evolving barcode）」；他們用 bulk methylation array 找出 978 個 fCpGs，並用 β-mixture model 把每個樣本的 fCpG methylation 分佈離散化成 0/1/2（0%、50%、100% 等位甲基化狀態），再用分佈形狀推回族群演化史；也用 nanopore 長讀序確認 fCpG 的變動不是由底層 somatic mutation 造成，並建立 fCpG methylation haplotypes。InterSubMod 目前是「以 SNV 視窗取 CpGs 來聚類」；EVOFLUx 暗示你也可以反過來：先找「最像條碼」的 CpG 集合（高異質、近似中間甲基化、且跨樣本/區域可比較），再把它們用於亞克隆/演化推論，而 SNV 只作為校準與生物解釋的輔助。citeturn19view3

最後，t-nanoEM（targeted long-read methylation + hybridization capture）雖然用的是轉換式（EM-seq）長讀序而非 MM/ML，但它在「臨床低 input」「目標區域高深度」「提供 haplotype 與 mutated allele-specific methylation workflow」上很有啟發：作者明確寫到因為 converted reads 的相位分析缺乏現成工具，所以他們基於 WhatsHap 自建 pipeline；並展示 t-nanoEM 在目標 panel 的 fold-enrichment、bait coverage、與 short-read EM-seq 的高相關，且可以在臨床乳癌/肺癌切片的多區域取樣中重現腫瘤/非腫瘤差異 DMR。這件事對 InterSubMod 很重要，因為它指向一條可行的臨床落地路線：若全基因體長讀序太貴或 coverage 不穩定，改用 targeted panel + 高深度，反而更適合做「second hit + promoter methylation + LOH」的臨床問題。citeturn17view0turn17view3

## 共識、分歧與尚未被清楚驗證的關鍵點

在「共識」上，至少有三點已相當明確。第一，癌症甲基化訊號強烈受 purity 與 CNV/LOH 影響；不做校正容易把混合訊號誤判成腫瘤特異甲基化，TRACERx 的 CAMDAC 與其用 LOH-SNV 驗證 PDR，是非常直接的證據鏈。citeturn20view0 第二，長讀序的價值不只在 methylation，而是「把相位、突變、結構變異、甲基化放在同一分子上」讓你能確認 biallelic event、allelic methylation、與第二擊（例如 promoter methylation + LOH），Long-Read POG 的多個案例（TSG biallelic、BRCA1/RAD51C promoter methylation + LOH）就是臨床導向的示範。citeturn13view0turn13view3turn12view0 第三，Nanopore 的 methylation/錯誤型態會反過來影響 variant calling 與 phasing；例如 LongPhase 的方法段落直接討論 Nanopore 在甲基化位點容易出現系統性 SNP 誤呼叫，並設計 heuristic 降低影響，顯示「甲基化 ↔ 變異」之間存在技術層的耦合，做 somatic 解析時必須正視。citeturn1view1

在「分歧/未清楚驗證」上，InterSubMod 目前遇到的 TP/FP 鑑別力瓶頸，其實正好反映了學界仍在爭論或尚難充分驗證的點。其一，局部 CpG pattern 的 cluster–label 顯著性，究竟多大程度代表「生物學上有意義的亞克隆/調控事件」，而不是 PMD/異染色質、重複序列、mapping bias、或 copy-number 造成的偽結構；InterSubMod 已用 PMD gating、dispersion/permanova、LOH/coverage quality 分層試圖處理，但從報告看仍不足以把鑑別力推到超越 baseline 的程度。fileciteturn39file0L118-L176 其二，「甲基化驅動」在癌症中的定義其實很多樣：TRACERx 用 \(M_R/M_N\) 與表現關聯去靠近「受選汰的功能性高甲基化」；MethSig（Cancer Discovery）則是用統計模型校正「背景隨機高甲基化率差異」來偵測 driver methylation；EVOFLUx 又指出有些 CpG 的波動可能是近乎中性的條碼訊號，重點在 lineage tracking 而非功能調控。這三種觀點彼此互補，但也意味著 InterSubMod 若把目標同時設定成「TP/FP」「driver methylation」「subclone lineage」，可能需要在方法上明確拆成不同任務與不同證據標準。citeturn20view2turn5view0turn19view3

## 對 InterSubMod 的推論：目前想法、差一個關鍵什麼、以及可往哪裡突破

綜合 repo 的架構文件、顯著性分析設計、以及 TP/FP 的驗證報告來看，我會把 InterSubMod 的「目前想法與目標」歸納成一個很具體的假設鏈：

你希望在「每個 somatic SNV 的局部視窗」內，讀序甲基化圖樣能自然形成穩定群集；若該 SNV 是真實腫瘤事件或與亞克隆/等位機制相關，則群集應能與 ALT/REF、HP1/HP2、Tumor/Normal 的標籤呈現統計顯著關聯；反之若是 FP 或純噪訊背景，群集與標籤的關聯應較弱或不可重現。這條假設鏈的工程實作已相當完整（距離、聚類、global/local/permanova/label-first、LOH/coverage/quality 分層）。fileciteturn34file0L20-L35 fileciteturn33file0L67-L110

但驗證顯示目前最大的缺口是：**「以單一 SNV 視窗為單位」的顯著性訊號，太容易被 LOH/CNV、區域性技術偏差、或低資訊量（CpG 太少、共同覆蓋不足）稀釋或誤導**。你已經觀察到 FP 的 LOH 率更高、且 LOH subtype 的分布差異顯著，這其實暗示了一個更核心的突破方向：**不要把 LOH/CNV 當成事後的「品質扣分」，而要把它提升成「生成模型的一部分」**。fileciteturn39file0L118-L151

下面是我認為最可能帶來「質變」的幾條路（每條都盡量對齊現有研究可借鑑的方法），並標註哪裡已有先例、哪裡仍是空缺值得你做出差異化。

第一條突破是做「copy-number / purity aware 的 read-label 建模」，把 CAMDAC 那種思想引入 InterSubMod 的 per-SNV 框架：你已經有 tumor/normal、也有 per-read ALT/REF/UNKNOWN，遇到 LOH 區域時甚至擁有「把 reads 分成更接近腫瘤來源 vs 正常來源」的天然標籤；下一步應是把這個標籤拿來估計「腫瘤純度/正常污染對 methylation 的混合比例」，並把校正後的腫瘤 methylation（或 PDR/heterogeneity 指標）作為聚類與顯著性檢定的輸入，而不是直接在混合訊號上聚類。TRACERx 已示範可以用 LOH-SNV 來驗證 deconvolution 的正確性（PDR correlation >0.8），因此這條路有很強的外部證據支持。citeturn20view0

第二條突破是把「SNV 視窗顯著性」升級成「基因/機制層級的證據整合」，類似 Long-Read POG 的做法：POG 把 CNV、allelic methylation（NanoMethPhase）、ASE（MBASED）與 somatic calls 整合（甚至提供 IMPALA 後處理流程），並能回答像「BRCA1 promoter methylation + LOH → HRD」這樣具臨床意義的 second-hit 敘事。InterSubMod 若把輸出從「region significance」再往上匯總成「gene evidence panel」（例如 TSG 的 biallelic event、DNA repair gene 的 promoter methylation + LOH、或 ecDNA 上的 methylation 異常），你在 TP/FP 的問題上可能也會更接近臨床可用的判讀維度。citeturn13view2turn13view3

第三條突破是重新思考「哪些 CpG 才該被當成亞克隆線索」，把 EVOFLUx 與 MHB（methylation haplotype block）概念納進來：EVOFLUx 指出 fCpGs 可以作為中性演化條碼，適合 lineage tracking；HGG Advances 的 MHB/LD-R² 研究則顯示癌症樣本的 methylation haplotype block 更小、methylation LD 更低，並能在變異附近找到 haplotype-specific methylated regions 且與 ATAC-seq 的 allele-specific accessibility 一致。這兩篇共同暗示：如果你的目標是「亞克隆/演化」，你可能要優先挑選能最大化 lineage 訊號的 CpG 集合（高異質、低相關、近似獨立波動），而不是把視窗內所有 CpG 平等看待；如果目標是「調控/驅動」，你可能要把 CpG 依 regulatory/nonregulatory 或 TF motif/ChIP peak proximity 做功能分層，再做類似 \(M_R/M_N\) 的選汰推論。這會讓 InterSubMod 從「聚類顯著」走到「可解釋的驅動機制」。citeturn19view3turn21view0turn20view2

第四條突破是把「somatic haplotyping」與「tumor-only」慣用的偵錯策略納入你的長讀序 pipeline，讓 second hit 的先後順序推論更可靠。LongPhase-S 的核心主張是：腫瘤不是乾淨的二倍體，subclonal heterogeneity、aneuploidy、contamination 會破壞一般 diploid phasing 的假設；因此它主打用 somatic haplotypes 來做 purity estimation，並用 purity-aware 的方式 recalibrate somatic variants，據稱能提升 ClairS 與 DeepSomatic 的 F1。即使你拿不到完整 preprint 內容，光這個問題定義就非常契合 InterSubMod：你現在已經看到 LOH/coverage 是 FP 指標，代表「相位/純度/拷貝數」的建模不足會直接反映在 variant truthfulness 上。對 tumor-only 的情境，Nature Communications 的 ClairS-TO 則提供了另一個可借鑑的工程路徑：用雙網路（affirmative/negational）輸出 posterior probability，並結合 panel of normals、以及使用 purity/ploidy/copy-number profile 做統計分類（germline/somatic/subclonal somatic）。這些策略都可以成為 InterSubMod 的上游品質控制，減少你在下游用甲基化「救火」的壓力。citeturn15search0turn16search0turn16search2

最後，若你把「second hit 的先後順序」當成主要科研亮點，InterSubMod 目前其實已經具備一個關鍵優勢：它能在 read-level 同時看到（1）突變/等位標記，（2）HP/PS，（3）局部 methylation pattern。Long-Read POG 已展示「雙突變 cis/trans」「突變與 LOH/小缺失的組合」的解析價值；TRACERx 展示「用 LOH-SNV 來切割 tumor vs normal 訊號」的驗證方法。InterSubMod 若能在方法上更明確定義幾種 second-hit archetype（例如：突變→LOH、LOH→突變、promoter methylation→LOH、LOH→promoter methylation），並設計可在 tumor/normal（或 tumor-only）資料下辨識的「機率式」判別規則（不是硬閾值），加上上層 gene-level 證據整合，你就會非常接近「已有人做但尚未真正產品化」的空缺地帶：可在臨床/研究實務中，面向長讀序做 biallelic mechanism reconstruction。citeturn13view0turn20view0turn39file0