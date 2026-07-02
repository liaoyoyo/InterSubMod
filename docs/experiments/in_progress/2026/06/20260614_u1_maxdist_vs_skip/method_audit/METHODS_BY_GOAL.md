---
title: "ISM 方法 — 以目標/子問題為先的敘述（目標→方法→敘述→解釋→結果→有效→替代）"
date: 2026-06-14
type: methodology / goal-oriented inventory
status: DRAFT — 逐個驗證的清單；待用戶確認後一個個驗
format: 每方法 7 欄 [目標/子問題 | 方法 | 敘述 | 簡單解釋 | 結果(實際數據) | 是否有效 | 替代方法]
related: TWO_DIRECTIONS_FRAMEWORK / COMPUTE_FLOW_ARCHITECTURE / SIGNIFICANCE_INVENTORY
claim: 敘述/結果=L1(源碼+本session數據) / 有效·替代=L2(待逐一驗證)
---

# ISM 方法（以目標為先）

> 7 欄模板：**目標/子問題 → 方法 → 敘述 → 簡單解釋 → 結果 → 是否有效 → 替代方法**
> 方向標記：A=label-first（標籤→甲基差異）/ B=structure-first（甲基結構→對應標籤）

---

## 目標 1：germline ASM 是否存在（haplotype 甲基差異）
**子問題**：HP1 vs HP2（parental haplotype）的甲基有沒有顯著差異？

- **方法 1a｜Δβ family（A，建議新增）**
  - 敘述：mean(HP1-family 甲基率) − mean(HP2-family) + permutation 檢定。
  - 簡單解釋：兩條 parental 染色體甲基平均差多少，差多少算顯著。
  - 結果：目前**只有 distance 版**（HPMergedDelta/P），Δβ 版未實作。
  - 有效？：✅ 預期有效（germline ASM 是真生物；HP-AUC normal 0.788 證距離抓得到 HP）。
  - 替代：DSS/methylKit beta-binomial（Δβ 顯著性的金標準）。

- **方法 1b｜LabelTest HP distance delta（A distance，現有）**
  - 敘述：按 HP 分組，between−within 距離 delta + 999-perm。
  - 簡單解釋：HP1/HP2 兩群 read 在甲基「距離」上分不分得開。
  - 結果：HPMergedDelta/P/Sig 已輸出。
  - 有效？：🟡 可用但**與 PERMANOVA-label_hp 冗餘**（#6）。
  - 替代：PERMANOVA pseudo-F（更標準）。

- **方法 1c｜GlobalTest HP（B，現有）**
  - 敘述：cluster_label × HP contingency → Fisher-FH + Cramér's V（Cochran gate）。
  - 簡單解釋：甲基自己分的群，跟 HP 標籤對不對得上。
  - 結果：GlobalP/CramersV 已輸出。
  - 有效？：✅（B 視角；但 cluster 為中介，HP-AUC 更直接）。
  - 替代：HP-AUC（B 直接版）、ARI。

---

## 目標 2：somatic ASM（tumor-specific 甲基差異，扣 germline baseline）
**子問題**：tumor 比 normal 多出來的 HP 甲基差異（不是 germline 本來就有的）？

- **方法 2a｜somatic residual Δβ + 檢定（A，#3 修法核心，建議）**
  - 敘述：(tumor: mean HP1−mean HP2) − (normal: mean HP1−mean HP2) + permutation（shuffle T/N）。**值已在 `HP_Signed_Residual`:1031**，缺檢定。
  - 簡單解釋：tumor 的 HP 甲基差，扣掉 normal 本來就有的那份，剩下的才是 somatic。
  - 結果：值已輸出（描述），**無 p/sig**。
  - 有效？：✅ 預期最對（somatic ASM 正解 + 對齊外部 DSS interaction）。
  - 替代：beta-binomial GLM `methyl~HP+sample+HP:sample`（interaction = somatic）。

- **方法 2b｜hp_residual（A distance，現有，缺陷）**
  - 敘述：tumor_hp.delta − normal_hp.delta（距離 delta 相減），sig 借 tumor p。
  - 簡單解釋：同 2a 但用距離、且顯著性是抄 tumor 的（沒真檢定殘差）。
  - 結果：HP_Residual_Delta/P/Sig 已輸出，**44.9% sig 不可信**（gap2）。
  - 有效？：🔴 **缺陷**（借 tumor p，距離相減意義不明）→ 應由 2a 取代。
  - 替代：2a（Δβ 殘差 + 真檢定）。

- **方法 2c｜SampleASM（A，現有）**
  - 敘述：tumor vs normal read 整體甲基距離 delta + perm。
  - 簡單解釋：tumor 跟 normal 的甲基整體像不像（不分 HP）。
  - 結果：SampleASM_* 已輸出（U1: 99.8% sig）。
  - 有效？：🟡 太敏感（99.8% sig，可能含 batch/depth confound）；非 HP-controlled。
  - 替代：2a（HP-controlled somatic）。

---

## 目標 3：subclone（haplotype 內 germline vs somatic 結構）
**子問題**：同一個 HP 內，germline read vs 帶 somatic 的 read（HP1-1/HP2-1），甲基是否不同？

- **方法 3a｜fine same-hap Δβ + 檢定（A，建議，最直接）**
  - 敘述：mean(HP1) − mean(HP1-1) 與 mean(HP2) − mean(HP2-1) + perm。
  - 簡單解釋：同一條染色體上，帶 somatic 的 read 甲基有沒有跟 germline read 不一樣。
  - 結果：**未實作 Δβ 版**（只有 distance pairwise HPFineD_HP1_HP1S）。
  - 有效？：✅✅ 預期是 label-first 偵測 subclone 的核心（直接、可解釋）。
  - 替代：fine PERMANOVA（B 版，已有 HPFineF/P）、beta-binomial（germline vs somatic carrier）。

- **方法 3b｜HP fine PERMANOVA（B，現有）**
  - 敘述：4 群（HP1/HP1-1/HP2/HP2-1）pseudo-F + perm。
  - 簡單解釋：四群 read 在甲基距離空間整體分不分得開。
  - 結果：HPFineF/P/Sig + NGroups 已輸出。
  - 有效？：🟡 整體顯著但不指方向（哪兩群差）；需配 3a 的 pairwise 才知是不是 same-hap subclone。
  - 替代：3a（pairwise Δβ 指方向）。

- **方法 3c｜clustering 子分群（B，P3b）**
  - 敘述：對單一 HP 的 read 再 clustering，看是否再分兩甲基群。
  - 簡單解釋：HP1 的 read 自己會不會再分成兩群（= HP1 內 subclone）。
  - 結果：**未實作**（plan P3b）。
  - 有效？：✅ B 視角偵測 subclone；但缺 ground truth（gap2 困境）。
  - 替代：3a（label-first，有 somatic-carrier 標籤當 ground truth）。

---

## 目標 4：甲基結構有效性（距離是否抓到真結構，非雜訊）
**子問題**：甲基距離能不能還原已知標籤（sanity floor）？

- **方法 4a｜HP-AUC（B，已實作 P1a，最乾淨）**
  - 敘述：P(dist(異HP對) > dist(同HP對))，NaN 跳過。
  - 簡單解釋：距離有沒有把同 HP 的 read 拉近、異 HP 推遠。
  - 結果：✅ 全基因組 normal **median 0.788, 76% strong**（地基確認）。
  - 有效？：✅ 已驗（C++/python 交叉一致）；B 方向最乾淨（不經 clustering）。
  - 替代：silhouette(HP label)、ARI。

- **方法 4b｜ClusterPermanova + silhouette（B，現有）**
  - 敘述：cluster 在距離空間 pseudo-F + perm；silhouette 選 k。
  - 簡單解釋：甲基自己分的群是不是真的有結構（非硬切）。
  - 結果：ClusterPermanovaF/P 已輸出（U1: SKIP 大量 sig）。
  - 有效？：✅ 但 99-perm（p floor 0.01）+ overfit 風險；HP-AUC 為前置 sanity。
  - 替代：gap statistic、HP-AUC（前置）。

---

## 目標 5：CpG 定位（哪些位點驅動差異）
**子問題**：是哪些單一 CpG 在 HP1 vs HP2 顯著？

- **方法 5a｜PerCpgAsm Fisher + BH-FDR（A，現有）**
  - 敘述：per-CpG 2×2(HP×meth) Fisher exact + BH-FDR。
  - 簡單解釋：逐個 CpG 看 HP1/HP2 甲基比例差不差，多重校正。
  - 結果：Fisher_N_Sig/Frac_Sig 已輸出。
  - 有效？：🟡 合理但 Fisher 忽略 over-dispersion（#4）。
  - 替代：DSS/methylKit beta-binomial per-CpG（更 powerful）。

---

## 目標 6：epiallele 異質性（甲基模式多樣性）
**子問題**：每個 haplotype 的 read 甲基模式有多亂（subclone 混雜的代理）？

- **方法 6a｜NME / Epipoly（A，現有）**
  - 敘述：NME=Shannon(read pattern)/Hmax；Epipoly=1−Σpᵢ²(4-CpG window)。
  - 簡單解釋：一條 haplotype 的 read 甲基模式越多種、熵越高 = 可能混了 subclone。
  - 結果：NME_HP1/2, Epipoly_*, Entropy_Imbalance 已輸出。
  - 有效？：✅ 標準 epiallele 指標（Jenkinson2020/Li2014）；異質性面向（與差異/定位互補）。
  - 替代：PDR、MHL。

---

## 逐個驗證順序（每個「對齊目標 + 確認數據」）
1. **目標 2/3 的 Δβ 模組**（2a somatic residual + 3a fine same-hap）— 最高槓桿（#3 修法 + subclone 核心），先做。
2. 實際數據驗：Δβ sig vs 舊 hp_residual_sig；fine same-hap Δβ 在 tumor both-HP region 的分布。
3. 對齊外部：與 DSS-like beta-binomial 結果比（(i)perm vs (ii)beta-binomial）。
4. 精簡冗餘（1b/2b distance）；目標 4/5/6 現有方法的細節驗（perm 數、over-dispersion）。
