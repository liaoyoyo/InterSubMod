---
title: "ISM 甲基-標籤驗證的兩個方向框架（label-first vs structure-first）+ #3 hp_residual 修法"
date: 2026-06-14
type: methodology-framework / 突破想法整理
status: DRAFT — 待用戶確認對齊後逐個驗證
claim_levels: 二分法=概念框架(用戶構想) / ISM 歸類=L1(源碼) / 修法+合理性評估=L2(提議, 待驗證)
related: SIGNIFICANCE_INVENTORY.md（6 方法盤點）, gap2/README.md（hp_residual 缺陷）
---

# 兩個驗證方向框架 — 對齊用戶構想

## 0. 核心二分法（你的突破想法，我的複述）
ISM 對「甲基 vs 標籤」其實有**兩種根本不同的運算與統計確認方向**，目前混在一起，必須分清：

| | **方向 A：label-first（標籤驅動）** | **方向 B：structure-first（結構驅動）** |
|---|---|---|
| **流程** | 先按**標籤**(HP/allele/T-N)分組 → 算組間甲基差異 → 檢定 | 先按**甲基結構**(read 距離)分群 → 再驗證「群」與標籤的關聯 |
| **問的問題** | 「這個**已知標籤**，在甲基上有沒有差異？」 | 「甲基**自己**會分成幾群？這些群對應哪個標籤 / 是不是新結構？」 |
| **性質** | 確認性 / supervised（標籤是輸入） | 探索性 / unsupervised（標籤是事後對照） |
| **核心量** | Δβ（甲基率差）或 distance delta（組間-組內） | cluster labels → 與標籤的關聯(CramersV) / 距離對標籤的還原(HP-AUC) |
| **對 subclone** | ❌ 不行（沒有 subclone 標籤可分組） | ✅ **正確方向**（從甲基結構發現群，再看是否 tumor-specific） |

**為什麼這是重要觀念**：A 量化「已知標籤的甲基效應」（如 germline ASM = HP 標籤的甲基差異）；
B 探索「甲基結構本身」（subclone 必須走 B，因為沒有 subclone 標籤）。**兩者回答不同問題、各有 confound、
各有驗證方式** — 混用會誤判（把 B 的 cluster 結果當 A 的標籤差異，或反之）。

## 1. ISM 現有方法歸類（L1 源碼，← SIGNIFICANCE_INVENTORY）

| 方向 | ISM 方法 | 算什麼 |
|---|---|---|
| **A label-first** | LabelTest distance delta（HP/allele/sample）| 按標籤分組，組間-組內距離 + permutation |
| | **Δβ signed delta** | 按 HP 分組，mean(HP1)-mean(HP2) 甲基率（**無檢定**）|
| | PerCpgAsm Fisher | 按 HP 分組，per-CpG 2×2 Fisher+FDR |
| | SampleASM / hp_residual | 按 sample(T/N) 分組 + 扣 normal baseline |
| **B structure-first** | clustering(甲基距離 UPGMA + cut) | 不看標籤，純甲基距離分群 |
| | GlobalTest(cluster × label) | cluster 與標籤的 contingency(Fisher+CramersV) |
| | ClusterPermanova | cluster 在距離空間是否真有結構 |
| | **HP-AUC**(我加的, P1a) | 距離**直接** vs 標籤(不經 clustering) — B 方向最乾淨的關聯度量 |

**🔑 關鍵釐清**：HP-AUC（P1a）其實是**方向 B 最乾淨的驗證** — 它直接量「甲基距離能否還原標籤」，
少了 clustering→cut→CramersV 那層 overfit 風險。U1/gap2 的發現（SKIP 去污染後 cluster 對應 HP）
就是 B 方向的結果（甲基結構 → 對應 germline-HP，非 subclone）。

## 2. 各方向的 confound + 驗證方式（是否合理有效，L2）

### 方向 A（label-first）
- **適用**：標籤可信時，量化標籤的甲基效應（germline ASM、somatic ASM）。
- **confound**：①標籤本身錯（HP tag 錯誤）②需 baseline（germline HP 差異 ≠ somatic；故要扣 normal）③read 內 CpG 相關 / over-dispersion。
- **驗證方式 + 合理性**：
  - Δβ + **檢定**（beta-binomial GLM / arcsine-t / permutation）— ISM **只有 Δβ 描述、缺檢定**（gap #1）。
  - distance delta + permutation — ISM 有（LabelTest），合理但與 PERMANOVA 冗餘（#6）。
  - Fisher per-CpG — ISM 有，但 over-dispersion（#4，beta-binomial 更佳）。

### 方向 B（structure-first）
- **適用**：探索未知結構（subclone）；驗證「甲基結構是否對應某標籤」。
- **confound**：①clustering overfit（強行分 k 群）②cluster-label 關聯 artifact ③k 選擇敏感。
- **驗證方式 + 合理性**：
  - **HP-AUC**（距離直接 vs 標籤）— 最乾淨，不經 clustering（✅ 已加，全基因組 normal 0.788）。
  - silhouette + PERMANOVA（cluster 結構顯著）— ISM 有，合理（Anderson）。
  - cluster × label CramersV（GlobalTest）— ISM 有，合理（Cochran gate）。
  - ⚠ **subclone 缺 ground truth**：B 對 subclone 只能用 HP/sample 標籤**間接**對照（gap2 困境）。

## 3. #3 hp_residual 殘差檢定怎麼修（你要先確認的）

**現狀缺陷**（RegionProcessor:976-993）：
- `hp_residual_delta = tumor_hp.delta − normal_hp.delta`（**距離 delta 相減**，非線性可加量，統計意義不明）
- `hp_residual_sig = tumor_hp.significant`（**借 tumor 自己的 p，根本沒檢定殘差**）

**修法（推薦走方向 A 的 Δβ 線，一次統一 #1+#3+Δβ內建）**：
- somatic residual Δβ = `(tumor: mean HP1 − mean HP2) − (normal: mean HP1 − mean HP2)`
  = tumor ASM 扣掉 germline ASM baseline = **tumor-specific(somatic) HP 甲基差異**（這才是 somatic ASM 的正解）。
- **檢定方式（兩個候選，待驗證哪個合理）**：
  - **(i) permutation**：shuffle 「sample 標籤(T/N)」within 各 HP，重算 residual Δβ，看 observed 是否極端 → p。
    （直接、無分布假設、對齊 ISM 既有 permutation 風格）
  - **(ii) beta-binomial GLM**：`methyl ~ HP + sample + HP:sample`，**interaction term HP:sample 的顯著性** = somatic ASM。
    （DSS/methylKit 標準、考慮 over-dispersion、更 powerful，但需引入 GLM）
- **這同時回答你的另一問**：是的，ISM **需要內建 Δβ 運算 + 對 Δβ 做檢定**（目前只描述）。
  somatic residual Δβ 是「方向 A 對 somatic ASM 的正確驗證」。

## 4. 對齊確認 + 逐個驗證計劃（你要「一個個驗證方法/細節/理論依據/實際狀況」）
**先確認對齊（§0-§1 二分法 + 歸類對嗎），再逐個驗證**：
1. **#3 修法選型**：permutation(i) vs beta-binomial(ii) — 各自理論依據 + 對 ISM 既有架構的契合 + 實際數據試算。
2. **Δβ 內建**：把 Δβ + 檢定做成方向 A 的標準輸出（per-region somatic residual Δβ + p + sig）。
3. **方向 A 各驗證方式逐一驗**：Δβ 檢定 / distance delta / Fisher 在實際數據上的一致性 + 差異。
4. **方向 B 各驗證方式逐一驗**：HP-AUC / silhouette / PERMANOVA / CramersV 在實際數據上是否一致指向同結論。
5. **A vs B 交叉**：同一 region，A（標籤甲基差異）與 B（結構對應標籤）結論是否一致？不一致代表什麼？
6. 每步「對齊想法/目標 + 確認數據/實際狀況」。

> ⚠ 本框架是**對齊工具**，不是定論。§0-§1 若與你構想有偏差，先校正再往下。
