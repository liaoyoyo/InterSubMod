<!--
建立時間: 2026-06-12
報告類型: 碩士論文正文草稿 — 第三章 材料與方法
狀態: draft v1（依架構 spec 方案1 / 標題定位 B characterization）；待用戶 review 校準語氣與深度
data_sources:
  - docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md (V1-V12 真值)
  - docs/method_comparison/20260609_ism_vs_external_methylation_tools/01_ism_method_spec_from_source.md (ISM 6 核心 + file:line)
  - docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md (資產盤點)
provenance_note: 方法參數(NHD/C_min=5/UPGMA/999perm/threshold 0.6/binary 0.2-0.8/cis-test)引自源碼 file:line(經 01_ism_method_spec 確認)；樣本級數字(unphase 45.84%/30.15M read)引自 VERIFIED_RESULTS V1。未實測之樣本統計(purity/depth/ploidy)、軟體版本、repo 連結以 {{待填}} 留白，不填預期值(§13)。
-->

# 第三章　材料與方法

本章描述本研究所使用之資料、方法與對照設計。需先界定本研究之定位：腫瘤亞克隆結構之重建（subclonal reconstruction）以既有之 somatic haplotagging 流程（longphase-S）所建立之單倍型骨幹（haplotype backbone）為基礎，而本研究之核心貢獻在於**特徵化（characterize）** 此骨幹之上 ONT 原生 per-read 甲基化訊號所能提供之**額外資訊及其解析度上限**——亦即量化甲基化在 germline-haplotype 層級、subclone 層級、以及做為 somatic 變異真偽過濾器等不同用途下之能力與邊界。換言之，重建之骨幹來自既有之 haplotagging，甲基化之加值與其界限則為本研究所欲釐清者。因此，本章將依序説明：所用之細胞株與定序資料（3.1）、做為重建骨幹與分析脈絡之 somatic haplotagging（3.2）、per-read 甲基化之矩陣表示（3.3）、用以特徵化甲基結構之 ISM（InterSubMod）六項核心分析（3.4）、以甲基化進行單倍型指派與相位救援之方法（3.5）、以及貫穿全研究、用以界定（bound）各項主張之對照與虛無模型設計（3.6）。最後説明統計分析原則（3.7）與軟體、資料及程式碼之可得性（3.8）。

全章所引用之具體數值，其唯一真值來源為已驗證真值紀錄（`VERIFIED_RESULTS.md`，文中以 V1–V12 標示）；方法參數與閾值則直接引自 C++/Python 原始碼，並標註檔案與行號。

---

## 3.1 樣本與定序資料

本研究使用六株具有配對正常樣本（matched normal）之癌症細胞株，涵蓋三種癌別，皆以 Oxford Nanopore Technologies（ONT）R10 long-read 平台定序，並於 basecalling 階段同時取得 5-methylcytosine（5mCG）之原生甲基化訊號：

- **乳腺癌**：HCC1395、HCC1937、HCC1954
- **肺癌**：H1437、H2009
- **黑色素瘤**：COLO829

每一細胞株皆具備三項要素：經 somatic haplotagging 之 tumor BAM（`{sample}_tagged.bam`）、somatic 變異呼叫結果（`somatic_pass.vcf.gz`，含 TP/FP 標籤），以及 per-read 甲基化（MM/ML tags）。各樣本之純度（purity）、覆蓋深度（depth）、倍體（ploidy）、體細胞突變數與拷貝數變異（CNV）涵蓋率等基本統計，整理於 Table 1（**{{待填：六樣本統計實測值；尚未彙整}}**）。

需特別指出本研究證據範圍之分層：第四章中以甲基化進行相位救援與單倍型輔助之深度驗證（V1–V12）目前係建立於**單一細胞株 HCC1395** 之上；跨六樣本之復現（將 V10 配對正常對照與 V11c 之亞群存在性測試推廣至全部樣本）為一獨立之開放工作，其結果不在本章之單樣本方法範圍內。配對正常樣本之原生甲基化目前於六樣本中有五株可得（HCC1395 同時具 5mC 與 5hmC，HCC1937、HCC1954、H1437、H2009 具 5mC），COLO829 之 R10 正常樣本缺 MM tag。

**參考真值（ground truth）之來源**：HCC1395 之拷貝數與雜合性缺失（loss of heterozygosity, LOH）區段採用 SEQC2 官方標註（`ngs_benchmark_cnvs_gain_loss_loh.bed`，660 segments；該樣本 ploidy 2.85）；germline 雜合位點（heterozygous SNP）取自配對正常樣本 HCC1395BL 之 Clair3 呼叫結果；somatic 變異之真偽（TP/FP）以 somatic VCF 之標籤為準。

---

## 3.2 Somatic haplotagging 骨幹（longphase-S）

本研究所特徵化之單倍型骨幹由 longphase-S 之 somatic haplotagging 流程產生。該流程先以配對正常樣本之 germline 雜合位點建立 germline 單倍型（HP1/HP2），再將 tumor reads 依其攜帶之 germline 與 somatic 證據指派至對應之單倍型或子單倍型。其 tag 指派之決策串接（cascade）見於原始碼 `HaplotagStrategy.cpp:452`（`judgeSomaticReadHap`）與 `SomaticHaplotagProcess.cpp:461`（`inheritHaplotype`），判定門檻 `percentageThreshold = 0.6`。

由此產生之 read 標籤具有六種狀態，其語意如下：

| 標籤 | 語意 |
|------|------|
| `unTag`（unphase）| 無足夠 germline 證據可定相之 read |
| `H1` / `H2` | 僅具 germline 證據、歸入 germline 單倍型 HP1 或 HP2 之 read |
| `H1-1` / `H2-1` | 同時攜帶 somatic 變異、且該 somatic 變異被歸入 germline 分支 HP1 或 HP2 之 read（即 somatic 子單倍型）|
| `H3` | 攜帶 somatic 變異、但其所屬 germline 分支未知之 read |

在 HCC1395 配對模式之全基因組統計中（V1，共 30,149,552 條 primary read），各標籤之分布為：unphase 45.84%（13,821,877）、HP1 25.35%、HP2 24.63%、HP1-1 2.00%、HP2-1 1.93%、HP3 0.25%。其中 unphase reads 構成本研究相位救援分析（3.5、第四章 R3）之主要對象池；HP3 在配對模式下極為稀少。

需強調此一單倍型骨幹係**既有流程之產物、為本研究之分析對象（脈絡）而非貢獻**；同時，由於本研究後續以甲基化「矯正」之 read 標籤本身即為 longphase-S 之輸出，故部分分析存在單一流程之自我參照（self-reference）特性，此一限制將於第五章討論並列為證據層級之上限依據。

---

## 3.3 Per-read 甲基化之矩陣表示（MethylationMatrix）

對每一分析區域（region），本研究自 BAM 之 `MM/ML` tags 取出每條 read 在各 CpG 位點之 5mC（及 5hmC）修飾機率，建構一以 read 為列、CpG 為行、元素為原始修飾機率之矩陣（MethylationMatrix）。

在需要二元化（binarization）之分析中，採用雙閾值：修飾機率 `≥ 0.8` 視為甲基化、`≤ 0.2` 視為未甲基化，介於兩者之間者視為不確定（ambiguous）（`DistanceMatrix.hpp:27-28`）。在 Python 層之 Δβ 計算中，以 ML 值 `≥ 200`（約 0.78）為甲基、`≤ 50`（約 0.20）為未甲基（`03_step4_ism_methylation_diff.py:30-31`），或於 `cis_asm_core.py` 採 `THR = 0.5`。當一條 read 在同一 CpG 同時具 5mC 與 5hmC 兩列訊號時，以取最大值方式塌陷（max-collapse）為單一「any-modification」訊號（`cis_asm_core.py:31-32`）；本研究主要結果由 5mC 驅動，5hmC 之邊際貢獻於討論中説明。

依分析需求，read 之 HP 標籤可整併為兩種分群方式（`LabelTest.cpp:244-305`）：

- **二分群（HP-merged）**：HP1-family（`1` 與 `1-1`）對 HP2-family（`2` 與 `2-1`）。
- **細四分群（fine 4-group）**：HP1 / HP1-1 / HP2 / HP2-1；HP3、HP0 與 unphased 於此分析中排除。

需註明一項目前之方法限制：MM/ML tag 之解析模組（`MethylationParser`，含正/反股訊號之處理）尚缺單元測試，本研究已於程式碼方法學稽核中列為待補之 golden test（涵蓋正反股、多 CpG 與邊界情形），相關不確定性於討論中誠實標註。

---

## 3.4 甲基結構之特徵化：ISM 六項核心分析

本研究以自行開發之 C++17 引擎 ISM（InterSubMod）量化 read 之間的甲基結構（between-molecule structure），亦即「哪些 read 屬於同一甲基亞群」之連結關係，而非單一 read 內部之失序（within-read disorder，如 epipolymorphism / PDR）。該引擎包含六項核心分析，整理如下（參數為原始碼預設值）：

**核心一　read–read 距離矩陣**（`DistanceMatrix.hpp:21-57`）。對區域內每一對 read 計算其甲基模式之距離，預設度量為正規化漢明距離（normalized Hamming distance, NHD），另提供 L1、Bernoulli、Pearson（mean-centered）與 Jaccard。每一對 read 需共同覆蓋至少五個 CpG（`min_common_coverage = 5`）方為有效，不足者其距離以懲罰值 `MAX_DIST = 1.0` 記錄而非捨棄；必要時可分正反股各自計算。此「顯式建構 N×N read↔read 距離矩陣」係本方法最具區辨性之一步，下游之 clustering 與結構檢定皆建立於此矩陣之上。

**核心二　階層式 subclonal clustering**（`HierarchicalClustering.hpp:23-207`）。以距離矩陣建構樹狀結構，預設採 UPGMA（隱含分子鐘假設），另提供 Ward、Single、Complete linkage。切樹方式可依距離閾值、固定群數 k，或以 silhouette score 於 k = 2…10 之間自動選取最佳群數（`find_optimal_clusters`）。

**核心三　PERMANOVA 結構檢定**（`StructureTest.cpp:141-205`）。檢定 HP 或 somatic 標籤所定義之分群在距離空間中是否呈現結構分離。其平方和分解為 `SS_total = (Σ_{i<j} d_ij²)/N`、`SS_within = Σ_group (Σ_{i<j∈group} d_ij²)/n_k`、`SS_between = SS_total − SS_within`；pseudo-F = `(SS_between/(k−1)) / (SS_within/(N−k))`。以 **999 次標籤排列（label permutation）** 建立虛無分布，p = (#{F_perm ≥ F_obs} + 1)/(999 + 1)。並另行進行離散度同質性檢定（各群至質心平均距離之單因子 ANOVA-F），用以辨別群間差異係源於位置差或離散度差。其衍生之 LabelTest 統計量 Δ = d_between − d_within（`LabelTest.cpp:213`）僅在 Δ > 0（群間較群內更遠、符合生物方向）時進行排列檢定。須特別區分：此 Δ 為**距離**之差（建於 read–read 距離矩陣），與 3.4 後段之 Δβ（甲基**率** β 之差）為不同量綱。

**核心四　per-CpG ASM 與失序度量**（`PerCpgAsm.hpp:3-98`）。同時計算三組文獻既有之度量供對照：per-CpG Fisher 精確檢定（搭配 BH-FDR 校正）、normalized methylation entropy（NME），以及四-CpG 滑窗之 epipolymorphism。此三度量係刻意實作既有文獻方法（Fisher per-CpG 對應 DAMEfinder/pycoMeth、NME 對應 CPEL/Jenkinson 2020、epipolymorphism 對應 methclone/Metheor）以供本方法之結構軸與之對比。

**核心五　Cramér's V 與可靠性閘**（`GlobalTest.cpp:128-141`）。計算 HP × 甲基狀態之全域關聯強度 Cramér's V，並以 Cochran 條件（列聯表最小期望格 ≥ 5）標記其可靠性；於稀疏列聯表時 `cramers_v_reliable = false`，以避免回報假性偏高之關聯。

**核心六　NormalBaseline 與 normal-anchored cis-test**（`NormalBaseline.hpp:38-67`；`34_control_loci_cohesion_cistest.py:108-152`）。以配對正常樣本之 per-CpG 平均甲基率扣除 germline 基線（residual = raw − normal mean），並以三組 per-CpG β——A = normal HP1（突變前基線）、B = tumor HP1（germline 單倍型、無 somatic allele）、C = tumor HP1-1（攜 somatic allele）——構成三角差：d_somatic = C − B（tumor 內 ASM）、d_cis = C − A（相對正常基線）、d_drift = B − A（germline 單倍型自身漂移）。判定為真 cis-driven 之準則為 `p_cis < 0.05` 且 `|d_cis| > 1.8×|d_drift| + 0.02`（`34_*.py:149`）。此設計於同一單倍型內比較，使 CN、ploidy 與比對特性維持恆定。

此外，於 Python/tsg 層計算 per-CpG 甲基率差 Δβ（`03_step4_ism_methylation_diff.py`、`cis_asm_core.py`）：將兩分群之 per-CpG β 依 CpG 配對後，計算平均差 `mean(β_som − β_germ)`、最大絕對差、Cohen's d（pooled SD）與 Wilcoxon signed-rank 配對檢定 p（並以 Mann-Whitney 非配對檢定為 sanity check），每分群需 ≥ 5 個 CpG、每 CpG 每分群需 ≥ 3 條 read。Δβ 沿兩條軸計算：HP-axis（HP1 對 HP1-1，同一 germline 母單倍型、germline-tag 對 somatic-重建-tag，屬 somatic-controlled）與 ALLELE-axis（ALT 對 REF）。

---

## 3.5 以甲基化進行單倍型指派與相位救援

為量化甲基化在「將無相位 read 指派回正確單倍型」此一任務上之能力，本研究採用以下方法，並一律以對照（3.6）界定其效應之真偽與強度：

**Anchor AUC（leave-one-out 質心法）**。於具有已知 HP1/HP2 標籤之 read（anchor）所在區域，遮去待測 read 之標籤，僅以其甲基模式相對兩單倍型質心之距離預測其歸屬，並與真實 HP 比對計算 AUC。此法用於量測「甲基模式能否區分兩單倍型」（V2、V3）。需特別説明：此類絕對 AUC 帶有方法樂觀成分（見 3.6 與第五章），故報告時一律改引相對於虛無模型之效應量，而非絕對 AUC 值。

**Held-out PS-block 外推（相位救援之直接模擬）**。由於 unphase read 本身無真值，「甲基能分 anchor」是否等同「甲基能救 unphase」屬外推，需獨立驗證。本研究以 germline VCF 之 phase-block（PS）做為獨立真值，在每一局部窗（±2 kb）中以一部分 read 學習 HP1/HP2 之甲基 profile，再對另一部分 read **遮去其 germline 證據（模擬 unphase）** 後僅以甲基預測其單倍型，並對照 PS 真值計算救援正確率與 shuffle 虛無模型（V6 設計，`extrapolation_validation.py`）。此設定如實反映真實救援情境：預測模型自其他 read 學習、被預測之 read 自身無 germline 證據。

**Allele-based REF/ALT AUC（獨立於 HP tag）**。為排除 longphase 跨 SNP 相位之依賴，於每一 germline 雜合 SNP 位點，依該位點之實際鹼基（REF/ALT）將 read 分為兩群計算甲基 anchor AUC（`allele_asm_auc.py`）。由於此標籤為單一位點之局部資訊、不需跨 SNP 相位，tumor 與 normal 之比較得以完全對等獨立，構成 3.6 中「非 copy」對照之基礎（V10）。

**Tag 輔助之三項目標（T1/T2/T3）**。在 somatic haplotagging 之 6-state 架構下，定義三項甲基輔助任務：T1 為將 unphase read 指派至 H1/H2；T2 為將攜 somatic 變異之 read 歸入 germline 分支 H1-1/H2-1；T3 為判別一條「未踩到 somatic 位點而被標為 H1/H2 之 read」究竟屬亞群（如 H1-1）抑或母本（H1）。各任務之 ground truth 來源不同（T1 為 held-out germline-masked、T2 為既有 germline 分支標籤、T3 於對抗驗證中改採該位點實際 ALT/REF 鹼基為亞群/母本標籤以求更乾淨之標籤），其驗證結果與限制詳見第四章 R4。

**可救援量之清點（unphase inventory）**。為界定相位救援之實際適用範圍，以單次掃描之 2 kb 分箱漏斗（`unphase_inventory.py`）統計 unphase read 中具足夠甲基覆蓋、且其鄰近具備 HP1/HP2 本地錨點（可建立本地參考）之比例（V12）。

---

## 3.6 對照與虛無模型設計

本研究之核心方法學原則為：**任何「甲基帶有某層級資訊」之主張，皆須以對應之對照或虛無模型界定其效應之真偽、來源與強度**；且因絕對 AUC 帶方法樂觀成分，效應量一律以「相對虛無模型之超出量」表述。所採用之對照如下：

- **germline-het 虛無模型**（V3，`germline_het_null.py`）：以與目標分析完全相同之方法，於隨機選取、非 ASM、非 imprinting 之 germline 雜合區計算 HP1-vs-HP2 甲基 anchor AUC 分布，做為「方法本身在一般 het 區可達之分離度」之基線。
- **shuffle-label 對照**（V4，`shuffle_control.py`）：於同一距離矩陣上，比較真實 HP 標籤與打亂標籤之 leave-one-out AUC，以 per-region 排列檢定（real 是否超過 shuffle 之 p95）判定訊號是否為標籤相依（label-dependent）而非方法洩漏（leak）。
- **配對正常樣本之 copy-clean 對照**（V10，`allele_asm_auc.py` + `imprint_bimodality.py`）：配對正常樣本 HCC1395BL 為 copy-clean 二倍體；若甲基分離度係由拷貝數造成，則正常樣本應較低。此對照用以決定性地判別觀察到之甲基差異是否為 copy artifact，並以 depth-matched 分析排除 read 數混淆。
- **SEQC2 CN/LOH 混淆控制**（V5，`seqc2_cn_methyl.py`）：以 SEQC2 之 CN 與 LOH 標註，檢定甲基分離度是否與拷貝數相關（Spearman）、以及 gain/loss/neutral/LOH 各狀態之分離度差異。
- **imprinting 正控**（V10）：於已知 imprinted DMR（GNAS、SNRPN、H19）以高斯混合模型檢定甲基之雙峰結構，做為「方法確能撈到已知真實 ASM」之陽性對照。

此外，本研究亦明確標示一項目前尚未完全排除之系統性混淆：CpG-SNP pseudo-ASM（germline 雜合 SNP 本身為 C/T 時改變 CpG 而製造假性 ASM）。shuffle-label 對照無法捕捉此一與真標籤相關之混淆，故其影響於第五章誠實討論並列為效應量解讀之限制。

---

## 3.7 統計分析

本研究之統計報告遵循計算方法評比之既有準則（Weber et al., Genome Biology 2019）：以效應量（effect size）搭配信賴區間或相對虛無模型之超出量為主要報告依據，而非單以 p 值判定；涉及分群結構與單倍型指派之檢定皆採排列檢定（permutation test，999 次或以上）建立虛無分布；多位點之檢定於模組內以 BH-FDR 校正，跨分析之多重比較於對應處標註；模型外推一律於獨立之 held-out / PS-block 真值上評估，避免純粹之留一法（leave-one-out）所造成之循環膨脹。所有絕對 AUC 之解讀皆對照其相對虛無模型之效應量。

---

## 3.8 軟體、資料與程式碼之可得性

**軟體版本**：longphase-S {{待填：版本}}、basecaller 與甲基化呼叫模型 {{待填}}、Clair3 {{待填}}、ISM C++ 引擎 commit {{待填：hash}}。

**程式碼可得性**：ISM C++ 引擎與全部分析腳本（3.4–3.6 所引之 Python/Shell）將公開於 {{待填：GitHub repository + commit hash}}，並存放於具永久識別碼之典藏（{{待填：Zenodo/Figshare DOI}}）。

**資料可得性**：六細胞株之 tumor/normal BAM、somatic VCF 與甲基化矩陣之存放位置與取得方式 {{待填：deposition statement}}。

> 註：本節之 {{待填}} 項為尚未彙整之事實性資料，依本研究之資料誠實原則留白，不以預期值填補；將於投稿/送審前補齊並各標來源。
