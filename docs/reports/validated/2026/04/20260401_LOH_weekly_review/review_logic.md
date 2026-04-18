<!--
建立時間: 2026-04-01 18:00
目標: 對 LOH 週報審閱系列（00-09）10 份文件的邏輯推論合理性進行嚴格科學審查
處理範圍: 假設合理性、替代解釋、統計推論、因果推論、跨文件一致性、結論強度
關聯檔案:
  - docs/reports/validated/2026/04/20260401_LOH_weekly_review/00_background.md ~ 09_conclusions_and_actions.md
-->

# LOH 週報審閱系列——邏輯推論審查報告

**審查日期**：2026-04-01
**審查範圍**：`00_background.md` ~ `09_conclusions_and_actions.md`，共 10 份文件
**審查原則**：對每個重要推論逐一挑戰其前提、替代解釋、統計解讀、因果方向與結論強度

---

## 質疑清單

---

### Q01

- **文件**: `04_mechanism_investigation.md` (Section 3.2)
- **推論**: 「TO somatic allele 造成 HP imbalance」——TP 的 somatic variant 天然偏向特定 haplotype，導致 HP_Ratio 偏離 0.5，因此 TO LOH 呈現 TP-enriched。
- **質疑**: 此推論是**間接推導**，缺乏直接 per-read 級別驗證。所呈現的證據鏈為：(a) TO TP extreme LOH 比例 44.6% > FP 35.9%；(b) TO TP LOH-like AF=0.580 > non-LOH AF=0.433；(c) 同位點 TO_only_LOH 是 paired 的 16-52 倍。但這三項證據都是相關性觀察，並非因果證明。替代解釋至少包括：(1) TO phasing 演算法本身對高 AF variant 位點有系統性 bias（與 somatic origin 無關）；(2) ClairS-TO 的 PASS filter 選擇性地保留高 AF TP（survivorship bias），而非 somatic allele 本身造成 HP imbalance；(3) 同位點 concordance 差異可能源自 TO 與 paired 使用不同 phasing block 結構，而非 somatic allele 驅動。文件第 7 節「不能宣稱的」中有部分承認，但推論主體（Section 3.2-3.5 的「機制總結」）仍以因果語氣呈現。
- **嚴重度**: 高
- **建議**: (1) 將 Section 3 標題從「機制假說：逐步推導」改為「機制假說：間接證據推導」；(2) 選取若干同位點，逐 read 比較 TO 與 paired 的 HP assignment，直接驗證 somatic allele 是否是 HP 偏斜的原因；(3) 檢驗 TO phasing 演算法是否對高 AF variant 有獨立的 bias（例如用 germline-only 位點做控制組）。

---

### Q02

- **文件**: `06_methylation_hypothesis_negative.md` (Section 2.4)
- **推論**: 「L2 Collider Bias」——對近常數特徵做 OLS residualization on AF 時，殘差 inherit AF 的信號方向，產生虛假 AUC。
- **質疑**: 數學描述基本正確（residual ≈ constant - beta * X），但**嚴謹度不足**。(1) 文件僅提供「數學直覺」而非嚴格推導。collider bias 在流行病學中有精確定義（conditioning on a collider variable introduces spurious associations between its parents），此處的現象更接近「induced confounding by adjustment」或「suppression effect」，而非經典 collider bias。術語使用不精確可能誤導後續研究者。(2) 文件提出的偵測標準（inflation = AUC_L2 - max_L3 > 0.10）缺乏理論推導。0.10 門檻是經驗值還是有統計基礎？(3) 文件聲稱此問題僅在「近乎常數」特徵時出現，但未提供 variance 門檻的量化標準。什麼程度的 variance 算「近乎常數」？(4) 文件未討論 L3 AF-bin stratification 本身的局限性：bin 內部仍可能有 AF gradient，且分 4 個 bins 會大幅降低統計力，導致真實但微弱的信號也被判定為 negative。
- **嚴重度**: 中
- **建議**: (1) 使用更精確的統計術語（如「adjustment-induced bias」或「regression artifact on near-constant features」），並附上正式的數學推導（variance decomposition）；(2) 為 inflation > 0.10 門檻提供 simulation-based 支持或 theoretical bound；(3) 定義「近乎常數」的 operational criterion（例如 CV < 0.05 或 ICC < 0.01）；(4) 討論 L3 的 power loss 問題並報告各 bin 的 sample size。

---

### Q03

- **文件**: `02_loh_evidence_panel_rounds1_4.md` (Round 4)；`05_systematic_observation_O1_O10.md` (O8)
- **推論**: 「LOH+HPMergedSig 87.5% 來自 HCC1395 chr8」，因此被歸結為「樣本特異性」現象。
- **質疑**: 此結論的方向正確（確實是 HCC1395 主導），但推論的完整性不足。(1) 「樣本特異性」不等於「完全是 artifact」。HCC1395 chr8 的 LOH + ASM 共現是真實的生物學現象（已知的 BRCA1 缺陷細胞系在 chr8 有 recurrent LOH）。真正的問題不是它「是否是 artifact」，而是「此生物學機制在其他腫瘤中的普遍性如何」。文件過早地否定了這個特徵的價值，卻未系統性調查此機制在其他高 LOH 腫瘤類型中是否重現。(2) 排除 HCC1395+HCC1954 後 enrichment 從 7.4x 降至 1.3x，但 1.3x 仍 > 1.0。這個殘餘 enrichment 是否統計顯著？未報告 p-value。(3) 80 個 FP 的樣本量極小，在 7 個樣本的分解下（0-70 per sample），per-sample 推論的統計力幾乎為零。「H1437 = 0 FP」不代表 H1437 不存在此現象，可能僅是因為 n 太小。
- **嚴重度**: 中
- **建議**: (1) 將結論修正為「HCC1395 chr8 LOH+ASM 驅動了本資料集中的大部分信號，但此生物學機制在高 LOH 腫瘤中可能具有普遍性」；(2) 報告排除 HCC1395 後的 1.3x enrichment 的 p-value 和 confidence interval；(3) 進行 power analysis 評估 80 FP 的 per-sample 分解是否有足夠統計力下結論；(4) 在更大的樣本集（如 TCGA LOH-enriched tumors）上驗證。

---

### Q04

- **文件**: `05_systematic_observation_O1_O10.md`；`09_conclusions_and_actions.md` (C4)
- **推論**: 「TO 無單一有效特徵」——所有 TO 特徵 AUC < 0.58，因此被判定為「無效」。
- **質疑**: (1) AUC < 0.58 作為「無效」的閾值缺乏明確依據。文件 `00_background.md` 6.8 節定義 AUC 0.6-0.7 為「弱區分力」，但 0.58 既不是 0.6 也沒有獨立的統計論證。為什麼 0.58 是臨界值而非 0.60 或 0.55？(2) 在 TO FP rate = 30.6% 的場景下（prevalence = 0.306），即使單一特徵 AUC = 0.57，若作為多特徵 ensemble 的一員，其增量貢獻可能是有意義的。文件 O5 已提到「甲基化與 caller 特徵近乎正交 (|r| < 0.26)，可提供獨立增量，但個別信號太弱」，但未量化 ensemble 上限。(3) 以 7 個樣本全部 < 0.58 來判定「無效」，但其中部分樣本（如 H1437 L3 max AUC=0.594, H2009 L3 max=0.597）已接近 0.60。0.58 的 cutoff 在此恰好將這些「borderline signal」劃為無效。(4) AUC 是 rank-based metric，對不平衡 prevalence 不敏感。在 TO FP rate 30.6% 時，更適合以 precision-recall AUC (PR-AUC) 或 calibrated likelihood ratio 評估。
- **嚴重度**: 中
- **建議**: (1) 明確定義「無效」的 operational threshold 並提供統計理由（如 bootstrap 95% CI 包含 0.5）；(2) 量化 2-3 個正交特徵的 ensemble AUC 上限（即使每個 < 0.58，組合可能 > 0.60）；(3) 報告 PR-AUC 作為補充指標；(4) 區分「no single useful feature」和「no useful feature」——前者不排除 ensemble 有效性。

---

### Q05

- **文件**: `06_methylation_hypothesis_negative.md` (Section 1)
- **推論**: O11 epipolymorphism AUC 從 0.845 降至 0.530（residualized），結論為「完全被 n_reads confound 驅動」。
- **質疑**: (1) Residualization 是否是唯一正確的去 confound 方式？OLS residualization（`feature ~ n_reads + num_cpgs`）假設線性關係且殘差同質。如果 epipolymorphism 與 n_reads 的關係是非線性的（如 log 或飽和曲線），OLS 殘差可能 over-correct 或 under-correct。文件未報告 OLS 的 R-squared 或殘差診斷。(2) Read-count-matched bin [81-120] 的 AUC = 0.560 雖確認了 residualization 結論，但此 bin 的 N 很小（TP=42, FP=117），且 TP:FP 比例嚴重不平衡（1:2.8 vs 全域 1.5:1），可能導致 AUC 估計不穩定。(3) 更重要的是：n_reads 與 truth label 的相關性（AUC=0.926）本身需要解釋。文件將其歸因於「paired manifest 的 sample selection artifact」，但如果此 artifact 也存在於實際應用場景中，n_reads 就不是需要「去除」的 confound，而是一個有效的 proxy feature。**是否應該去 confound 取決於應用場景：如果在實際使用中 TP 總是有更多 reads，那 n_reads（和因此 inflated 的 epipolymorphism）就是有用的 predictor。**
- **嚴重度**: 高
- **建議**: (1) 報告 OLS residualization 的 R-squared 和殘差 Q-Q plot，確認線性假設合理；(2) 使用非線性方法（如 LOESS 或 quantile matching）作為 robustness check；(3) 增大 read-count-matched bin 的範圍或使用多個 bins 以改善統計力；(4) 明確討論 n_reads confound 的因果結構——如果 n_reads 是 TP 的 causal consequence（高品質 variant calling 需要高 coverage），那麼 conditioning on n_reads 可能是不當的（conditioning on a mediator）。提供 DAG 圖示。

---

### Q06

- **文件**: `03_post_hp_fix_loh_enrichment.md` (Table 4)；`04_mechanism_investigation.md`
- **推論**: TO LOH enrichment 在所有 7 個樣本中一致 < 1.0，因此結論為「全面呈現 TP 富集」，強度標記為「確定」。
- **質疑**: (1) Table 4 顯示 COLO829 Tier A p=0.189 和 H2009 Tier A p=0.270，**不顯著**。全域 All-Tier 中雖全部 < 1.0，但 7 個 enrichment 值的範圍是 0.852-0.956。0.956 接近中性。文件將「方向一致但幅度差異大」的結果標記為「確定」，但未進行 meta-analysis 或 random-effects model 來正式評估 cross-sample 的一致性。(2) 所有 7 樣本來自癌症細胞系，不是獨立的隨機抽樣。細胞系之間可能共享培養條件、passage artifacts 等系統性偏差。7/7 一致的結果可能反映的是細胞系的共同特性而非 TO phasing 的普遍行為。(3) 文件明確指出 COLO829 是「低信心」，但結論 C2 (09_conclusions_and_actions.md) 仍以「7 個樣本一致 0.852-0.956x」呈現，未排除 COLO829。結論強度應反映最弱環節。
- **嚴重度**: 中
- **建議**: (1) 進行 random-effects meta-analysis（7 studies, effect = log(enrichment), variance from Fisher test），報告 pooled estimate 和 I-squared 異質性指標；(2) 將結論從「確定」調整為「強支持（6/7 樣本顯著，1 樣本方向一致但不顯著）」；(3) 在 C2 結論中標注 COLO829 的 caveat。

---

### Q07

- **文件**: `02_loh_evidence_panel_rounds1_4.md` (Round 4, Section 4)
- **推論**: Tier A(30-49) enrichment = 0.43x（TP 富集）與 Tier A+(>=50) enrichment = 2.02x（FP 富集），方向相反。
- **質疑**: (1) 此結果的生物學解釋（「30-49 reads 正確偵測 somatic LOH」vs「>=50 reads 是古老 LOH 捕捉 germline SNP」）是 post-hoc rationalization，缺乏先驗假說。相同的數據可以有其他解釋：例如 eff_hp 30-49 的 region 可能是 medium-depth 區域，caller 在此深度的表現本身不同。(2) **Tier A combined (>=30) enrichment = 1.169x 是加權平均結果**，但 0.43x 和 2.02x 的「混合」不一定等於 1.169x，需確認是否為 Simpson's paradox 的另一個案例。(3) Tier A 的 FP 數僅 1,065（其中 LOH-like 227），Tier A+ 的 FP 數更少。在如此小的 FP 樣本量下，enrichment 2.02x 的 confidence interval 是多少？(4) 此發現是 Round 4 的新發現，但 Round 2/3 報告中 Tier A(>=30) 被統一處理。如果 A 和 A+ 方向相反是真實的，那 Round 2/3 的所有 Tier A(>=30) 結論都需要重新解讀。
- **嚴重度**: 高
- **建議**: (1) 報告 Tier A(30-49) 和 A+(>=50) enrichment 的 95% CI；(2) 進行 per-sample 的 Tier A vs A+ 分解，確認方向翻轉是否在所有樣本中一致（還是被少數樣本驅動）；(3) 在 Round 2/3 相關結論處加入勘誤，提示 A(>=30) 是混合效應；(4) 探索更細的 eff_hp bins（如 30-39, 40-49, 50-59, 60+）來確認 enrichment 的轉折點。

---

### Q08

- **文件**: `07_qs_mode_aware_change.md` (Section 3.1)
- **推論**: 移除 TO LOH penalty 和 verify bonus 後，TO QS AUC 從 0.497 升至 ~0.546（+0.049）。
- **質疑**: (1) AUC 改善從 0.497 到 0.546 雖然是正向的，但 **0.546 仍然極低**，距離任何實用門檻（0.7+）仍有巨大差距。文件 Section 3.2 承認「效果有限」，但 Section 07 標題（「QS Mode-Aware 程式碼修改」）和 09_conclusions 中 C10 的措辭（「QS Mode-Aware 修改已完成」）可能給讀者留下過度正面的印象。(2) +0.049 的預估值來源不明。是 simulation 還是 component decomposition 推算？Section 3.1 列出 0.546 但標記「~」（近似），未報告此估計的方法和誤差範圍。(3) AUC 0.497 -> 0.546 的改善是否統計顯著？在 n=419,692 的資料量下，DeLong test 可能顯著，但效應量極小（Cohen's d 差異可能 < 0.05）。
- **嚴重度**: 低
- **建議**: (1) 在結論中明確強調「此修改消除了反向效果，但 TO QS 仍然無效（AUC ~0.55），需要全新 feature engineering」；(2) 報告 AUC 估計的方法（analytical decomposition vs empirical rerun）和 confidence interval；(3) 進行修改後的實際 benchmark rerun 以獲得精確的 post-fix AUC。

---

### Q09

- **文件**: `00_background.md` (Section 2.2)；`03_post_hp_fix_loh_enrichment.md` (Section 6)
- **推論**: 「LOH 的 enrichment 方向在 Paired 和 TO 模式下完全相反」，Paired 為 FP-enriched，TO 為 TP-enriched。
- **質疑**: (1) 「完全相反」的表述可能過於強烈。Paired enrichment 範圍 1.02-3.18x，其中 HCC1395 paired 僅 1.02x（幾乎中性）。TO enrichment 範圍 0.852-0.956x，其中 COLO829 TO 0.956x（接近中性）。在兩端各有接近 1.0 的樣本，用「完全相反」描述可能誇大了效應的一致性。(2) Paired 的 enrichment 受極端小樣本影響：H1437 僅 8 FP、HCC1954 僅 29 FP。這些樣本的 enrichment（1.79x、3.18x）統計不穩定，confidence interval 極寬。如果只看 FP > 100 的樣本（COLO829=2244, HCC1395=627, HCC1395_DORADO=240, HCC1937=195），enrichment 範圍縮窄為 1.02-1.50x，遠不如表面數字暗示的那麼極端。(3) 未對 Paired enrichment 進行全域 Fisher test 或 meta-analysis。
- **嚴重度**: 中
- **建議**: (1) 將「完全相反」修改為「方向一致地相反（Paired 全部 >=1.0, TO 全部 <=1.0），但效應幅度差異大」；(2) 報告 Paired enrichment 的 per-sample Fisher p-value 和 sample-size-weighted average；(3) 在呈現 per-sample 數據時加入 n_FP 和 confidence interval。

---

### Q10

- **文件**: `05_systematic_observation_O1_O10.md` (O7)
- **推論**: 「TO/paired HP_Ratio 完全不相關：Pearson r=0.001——兩套完全獨立的 haplotype assignment 系統」。
- **質疑**: (1) Pearson r=0.001 不代表兩套系統「完全獨立」。HP_Ratio 本身的分佈是高度偏態的（多峰結構，大量值在 0 或 1 附近），Pearson correlation 在此類分佈下嚴重低估真實關聯。應以 Spearman rank correlation 或 concordance（同位點 LOH/non-LOH 分類一致率）來衡量。(2) 事實上，04_mechanism_investigation.md Table 1 報告同位點 concordance 為 81.8-89.8%，這與「完全不相關」的說法矛盾。85%+ concordance 意味著兩套系統在 LOH 判定上有相當程度的一致性。(3) r=0.001 與 85% concordance 的矛盾可能源自 HP_Ratio 的雙峰分佈——大量位點在兩端（0 或 1）高度一致，但連續數值層面的 correlation 被中間區域的噪聲淹沒。
- **嚴重度**: 高
- **建議**: (1) 報告 Spearman rank correlation 作為替代；(2) 將結論從「完全不相關」修正為「連續數值層面低相關（Pearson r=0.001），但二元分類（LOH vs non-LOH）有 85%+ concordance」；(3) 解釋兩個指標的矛盾原因（雙峰分佈效應）。

---

### Q11

- **文件**: `06_methylation_hypothesis_negative.md` (Section 2.9)
- **推論**: 「此方向正式關閉」——甲基化無法區分 TO LOH 區域中的 germline 和 somatic variants。
- **質疑**: (1) O12 的資料規模為 175,542 rows，特徵為 22 個，使用 4 層 confound 控制。但 negative result 的結論強度取決於 **power to detect a true effect**。文件未報告 power analysis：在 175K rows 中，AUC 0.55 和 0.58 之間的差異（delta AUC = 0.03）需要多少樣本才能偵測？以此資料量，power 可能足夠偵測 AUC > 0.52 的信號，因此 negative result 在此水準上是可信的。但文件的斷語「正式關閉」比「在目前的特徵集和解析度下無有效信號」更強。(2) 22 個特徵是否窮盡了所有可能？文件指出 ISM 有 116 欄位，其中 42 個已深度分析。剩餘的未分析欄位中是否有 LOH-specific 的可能？(3) 關閉一個研究方向是重大決策。「corrected AUC 全 < 0.58」是否排除了非線性組合有效的可能性（如 decision tree 或 neural network 在 feature space 中找到非線性 boundary）？
- **嚴重度**: 中
- **建議**: (1) 將「正式關閉」改為「在當前線性特徵框架下正式排除，非線性組合或 read-level 方法仍為開放問題」；(2) 補充 power analysis（minimum detectable AUC given n=175K and alpha=0.05）；(3) 以 random forest 或 gradient boosting 在 22 features + AF 上做一次 quick pilot，確認非線性組合也無效後再下最終結論。

---

### Q12

- **文件**: `04_mechanism_investigation.md` (Section 5, Table 3)
- **推論**: TO FP LOH-like 的 AF 中位數 = 0.969（接近 homozygous），被解釋為「germline homozygous variant」。
- **質疑**: (1) AF=0.969 不一定代表 germline homozygous。在高 purity 腫瘤中，somatic variant 也可以有接近 1.0 的 AF（clonal variant in LOH region）。真正的區分需要 genotype information（GT field）而非僅靠 AF。(2) 文件將 TO FP 等同於 germline variant，但 TO FP 的組成可能更複雜——可能包含 sequencing artifacts、mapping errors、或 mosaic germline variants。將所有 FP 統一歸因於「germline」過度簡化了問題。(3) TO TP LOH-like AF=0.580 被解釋為「somatic clonality 的自然結果」，但 0.580 也可能反映 subclonal TP（非 LOH 環境下的 partial allele fraction）。將 AF 差異直接歸因於 LOH 機制，未充分考慮 tumor purity 和 clonality 的 confound。
- **嚴重度**: 中
- **建議**: (1) 結合 GT field（如果可用）來區分 germline homozygous（1/1）和 somatic high-AF；(2) 進行 TO FP provenance 分類（germline vs artifact vs other），不要假設 FP = germline；(3) 在 AF 差異分析中控制 tumor purity。

---

### Q13

- **文件**: `05_systematic_observation_O1_O10.md` (Top 10, #1)
- **推論**: 「TO 模式下無任何單一特徵有效區分 TP/FP」，結論等級為 Level A（跨觀察交叉驗證通過）。
- **質疑**: (1) Level A 的判定標準為「至少 2 個獨立觀察確認，且零矛盾」。但 O1、O4、O5、O8 的分析都是基於同一份 master dataset（748,391 rows），因此不是統計意義上的「獨立觀察」——它們共享相同的 sampling bias 和 measurement error。真正的獨立驗證需要外部數據集。(2) O10 的 read-level AUC = 0.547 for methyl_bimodal_score 被宣告為「無效」，但 O10 的 read-level 分析僅用了 86,521 reads（paired-pure only），未在 TO 模式下測試。因此「TO 無有效特徵」的結論排除了一個尚未被完整測試的方向。
- **嚴重度**: 低
- **建議**: (1) 將 Level A 的定義補充為「至少 2 個獨立 analysis pipeline 在同一 dataset 上確認」，並承認 external validation 尚缺；(2) 在 Top 10 #1 結論後附加「read-level 特徵在 TO 模式下尚未被完整測試（O10 僅覆蓋 paired-pure）」。

---

### Q14

- **文件**: `02_loh_evidence_panel_rounds1_4.md` (Round 3, Section 1)
- **推論**: HP0 filter 假說被否定——「High HP0 TP% = 76.7% > Low HP0 74.8%」。
- **質疑**: (1) TP% 差異僅 1.9pp（76.7% vs 74.8%），在 148K 的資料量下可能統計顯著但效應極小。未報告 p-value 和 effect size。(2) 更重要的是，此分析是在 TO Tier A 的 LOH-like region 中進行的。HP0 filter 的目標不是改善全域 TP%，而是識別「phasing quality 差的 LOH-like region」。即使 High HP0 的 TP% 略高，HP0 可能仍然是 phasing quality 的有用指標（例如在 HP0 threshold 與其他特徵的交互作用中）。(3) 因果方向不清：High HP0 → TP% 高可能是因為某些 TP-enriched 的 genomic regions（如 pericentromeric, low mappability）天然 HP0 更高，而非 HP0 本身是 TP 的指標。
- **嚴重度**: 低
- **建議**: (1) 報告 1.9pp 差異的 p-value 和 Cohen's h；(2) 測試 HP0 與其他特徵的交互作用（而非僅作為獨立 filter）；(3) 控制 genomic context（如 mappability）後重新評估。

---

### Q15

- **文件**: `08_literature_survey.md` (Section 3.2)
- **推論**: 文獻「部分支持」ISM 甲基化區分力弱的結論，理由是 ROCIT/MethylBERT 在特定條件下有效但 ISM 條件不滿足。
- **質疑**: (1) 文獻的引用是選擇性的。文件列出了支持 negative result 的文獻（Do et al. 癌症 ASM 增加、PMD stochastic drift），但未系統性搜索是否有文獻在類似條件下（ONT、region-level、tumor-only）成功使用甲基化區分 germline/somatic。如果沒有，應明確指出「無先例」；如果有，需解釋與本研究結果的差異。(2) ROCIT（2026 年 bioRxiv）尚未經過 peer review，引用時應標注 preprint 狀態。(3) Section 3.3 認為文獻「不矛盾但也不支持」O12 結論，這個雙重否定的措辭模糊。應明確表述為「O12 是此假說在 ISM 框架下的首次實證檢驗」。
- **嚴重度**: 低
- **建議**: (1) 補充系統性文獻搜索（PubMed/Google Scholar: "tumor-only" AND "methylation" AND "germline filter"），確認是否遺漏相關文獻；(2) 標注 ROCIT 的 preprint 狀態；(3) 將「不矛盾但也不支持」改為更明確的表述。

---

### Q16

- **文件**: `01_hp_integer_tag_fix.md` (Section 6.1)
- **推論**: 修正後 TO LOH enrichment 全面 < 1.0（0.85-0.96，TP 富集），修正前接近 1.0 或略 > 1.0。
- **質疑**: 跨文件數字一致性問題。(1) Section 6.1 報告修正後 enrichment 為「0.85-0.96」，但 03_post_hp_fix_loh_enrichment.md Table 3 報告的 All-Tier enrichment 範圍為 0.852-0.956x。兩者吻合。(2) 然而，01_hp_integer_tag_fix.md Section 6.1 表中列出 TO Tier A+ enrichment 修正前為「接近中性（~1.03）」，修正後為「0.70」。但 03_post_hp_fix_loh_enrichment.md Table 2 顯示 Tier A+ enrichment = 0.766x（非 0.70）。**0.70 與 0.766 的差異**需要解釋——是不同計算口徑（全域 vs per-tier）還是數據錯誤？
- **嚴重度**: 中
- **建議**: 核實 01_hp_integer_tag_fix.md Section 6.1 中「0.70」的來源數據，若與 03 中的 0.766 計算口徑不同，需標注一致的口徑；若為錯誤則修正。

---

### Q17

- **文件**: `05_systematic_observation_O1_O10.md` (O10)
- **推論**: 「FP variants 傾向高甲基化（非活性）區域：FP mean methyl=0.679 vs TP=0.463」。
- **質疑**: (1) 此觀察是在 paired-pure 的 86,521 reads 上得出的，但 O10 也提到「受 clustering 膨脹」。methyl_low_fraction AUC=0.737 在 region-level 上可能因為同一 region 的多條 reads 高度相關（within-region correlation）而被人為膨脹。文件承認「有效獨立 N 僅 620 regions」，但仍以 0.737 作為報告值。這個數字應以 region-level bootstrap 重新估算後再呈現。(2) FP mean methyl=0.679 > TP 0.463 的因果解釋是什麼？是 FP 傾向出現在高甲基化區域（genomic context confound），還是 FP 的 reads 本身有更高的甲基化水平（biological signal）？如果是前者，則 methyl_mean 只是 genomic position 的 proxy，而非獨立的 epigenetic evidence。
- **嚴重度**: 中
- **建議**: (1) 以 region-level bootstrap（N=620 regions）重新估算 AUC 並報告 95% CI；(2) 控制 CpG island / shore / shelf / open sea 的 genomic annotation 後重新評估 methyl_mean 的區分力；(3) 在「行動建議」中將此項從 P2 提升（或保持 P2 但標注需要 confound control）。

---

### Q18

- **文件**: `09_conclusions_and_actions.md` (Section 3, N4)
- **推論**: 「LOH 作為 binary FP filter」被否決，理由為所有 filter 的 F1 delta < 0。
- **質疑**: (1) F1 delta < 0 的根本原因是 paired FP 絕對數量太少（3,429 / 328,699 = 1.04%）。在此 prevalence 下，任何 binary filter 都幾乎必然 F1 < 0，因為 TP 被移除的損失（denominator 減少）遠大於 FP 被移除的收益（numerator 增加）。**這不是 LOH 特有的問題，而是 low-prevalence FP filter 的通用限制**。文件應明確區分「LOH 不能作為 FP filter（方向問題）」和「任何特徵都不能在 1% FP rate 下做 binary filter（base-rate 問題）」。(2) 在 TO 模式下（FP rate = 30.6%），binary filter 的 prevalence 條件更有利。但文件未報告 TO 的 F1 delta——可能因為 TO LOH 是 TP-enriched（方向反轉），移除 LOH 會損害更多。如果補充 TO 的 F1 delta 數據，可以區分「方向不對」和「base-rate 問題」。
- **嚴重度**: 低
- **建議**: (1) 在 N4 結論中明確指出 F1 delta < 0 同時受「方向問題」和「base-rate 問題」影響；(2) 補充 TO mode 的 F1 delta 作為對照；(3) 討論在不同 FP rate 下（如 5%, 10%, 20%）LOH filter 是否可能有效。

---

### Q19

- **文件**: `02_loh_evidence_panel_rounds1_4.md` (Round 2, Section 2)
- **推論**: Tier B 的 enrichment 0.901x 是 Simpson's paradox，由 COLO829 主導。
- **質疑**: (1) 文件對 Simpson's paradox 的診斷方向正確，但解釋不夠精確。COLO829 「貢獻了 Tier B 中 93% 的 FP（1,206/1,299）」，而 COLO829 本身 Tier B enrichment = 1.093x（正向）。那麼為什麼全域 Tier B = 0.901x（負向）？文件歸因於「非 COLO829 樣本的 TP LOH-like fraction 很高」但未提供數學驗證。Simpson's paradox 的標準呈現應是分組後每組 enrichment 方向一致但合併後翻轉，但此處 COLO829 = 1.093x（正向）而全域 = 0.901x（翻轉）。其他樣本各自的 Tier B enrichment 是多少？如果部分樣本 enrichment 確實 < 1，那就不是 Simpson's paradox，而是 sample heterogeneity。(2) 稱「Tier B 整體反轉不是真實的生物訊號」的措辭過強——它可能是真實但被 heterogeneous 的信號，而非純 artifact。
- **嚴重度**: 低
- **建議**: (1) 列出所有 7 個樣本在 Tier B 的 enrichment；(2) 正式驗證是否符合 Simpson's paradox 的定義（每個子群方向一致但合併後翻轉）；(3) 若部分樣本 Tier B 也 < 1.0，則改為「sample composition effect」而非 Simpson's paradox。

---

### Q20

- **文件**: `00_background.md` (Section 5)；跨文件
- **推論**: 整個研究基於 748,391 rows x 116 columns 的 master dataset，涵蓋 7 個細胞系。
- **質疑**: (1) **External validity 風險**：7 個細胞系全部是已建立的商業化癌症細胞系（HCC1395、COLO829 等），經過長期體外培養。這些細胞系的 LOH pattern、甲基化 landscape、基因組結構可能與原始腫瘤或臨床樣本有系統性差異（如 culture-induced methylation drift、selection for growth advantage）。文件 00_background.md 待驗證問題 #3 提到「cell line vs clinical sample」，但未在任何結論中加入此 caveat。所有結論都可能在臨床樣本上不成立。(2) **樣本多樣性不足**：7 個樣本中 3 個是乳腺癌（HCC1395、HCC1937、HCC1954），2 個肺腺癌（H1437、H2009），1 個黑色素瘤（COLO829），加上 1 個平台重複（HCC1395_DORADO）。有效獨立樣本可能只有 5-6 個。乳腺癌在 LOH 行為上可能有系統性偏向。(3) 所有統計推論（enrichment、AUC、concordance）都在 pooled dataset 上進行，但池化了不同腫瘤類型。cancer-type-specific 效應可能被平均掉。
- **嚴重度**: 高
- **建議**: (1) 在所有核心結論後加入「本結論基於 7 個癌症細胞系，推廣至臨床樣本需要外部驗證」的 disclaimer；(2) 在 O8 基礎上明確標注有效獨立樣本數（排除 HCC1395_DORADO 後 = 6）；(3) 在 09_conclusions_and_actions.md 中將「外部驗證計劃」列為 P1 行動項。

---

## 跨文件一致性檢查摘要

| 項目 | 文件 A | 文件 B | 數字 A | 數字 B | 一致性 |
|------|--------|--------|--------|--------|--------|
| TO LOH enrichment (All) | 01 Sec 6.1 | 03 Table 2 | 0.85-0.96 | 0.852-0.956x | 一致 |
| TO Tier A+ enrichment | 01 Sec 6.1 | 03 Table 2 | ~~0.70~~ → 0.766 (已修正) | 0.766x | **已解決**：01 原誤標 Tier A(30-49)=0.706 為 A+；已修正為 A(30-49)=0.706, A+(≥50)=0.766 |
| TO FP rate | 00 Sec 4.2 | 03 Sec 6.2 | 30.6% | 30.6% | 一致 |
| Paired FP enrichment | 00 Sec 2.2 | 02 Round 1 | 1.194x | 1.194x | 一致 |
| LOH penalty trigger TP | 07 Sec 1.2 | 05 O2 | 44.5% | 44.5% | 一致 |
| epipolymorphism AUC drop | 06 Sec 1.3 | 09 N1 | 0.845->0.530 | 0.845->0.530 | 一致 |
| HP concordance | 05 O3 | 04 Table 1 | 85.5% | 81.8-89.8% | 一致（05 為平均，04 為 per-sample 範圍） |
| TO/paired HP_Ratio r | 05 O7 | 04 concordance | r=0.001 | 85%+ concordance | **已解決**：haplotype swap 導致散點圖呈 X 形，正負 r 抵消；LOH binary 不受 swap 影響 |

---

## 總體評估

### 優點

1. **研究設計的系統性**：四輪 LOH Evidence Panel + 10 個系統性觀察單元，從基線盤點到假說檢驗再到 confound 控制，邏輯鏈清晰。
2. **Negative results 的嚴謹對待**：O11、O12 的否決過程包含多層 confound 控制，避免了過早放棄或過度宣稱。
3. **數據透明度**：每項結論都附有具體數字、來源檔案路徑、分析腳本，可追溯性高。
4. **自我批判意識**：文件中多處自我標記「待驗證」、「不能宣稱」、「低信心」，反映了謹慎的科學態度。

### 核心風險

1. **因果推論不足**：最重要的機制推論（TO somatic allele -> HP imbalance）缺乏直接 per-read 驗證，仍停留在相關性推導階段（Q01）。
2. **External validity 風險**：所有結論基於 7 個細胞系，無臨床樣本驗證（Q20）。
3. ~~**跨文件數字不一致**~~ **已解決**：Tier A+ enrichment 差異為 tier 標錯（01 已修正為 A=0.706, A+=0.766）；HP_Ratio r=0.001 與 concordance 85% 為 haplotype swap 效應（05 已補充解釋）。
4. **AUC < 0.58 = 「無效」的閾值缺乏理論基礎**：多個重要結論依賴此閾值，但其合理性未被論證（Q04）。

---

*本審查報告為獨立邏輯審查，不替代數據重分析。所有高嚴重度質疑建議在下一輪研究中優先處理。*
