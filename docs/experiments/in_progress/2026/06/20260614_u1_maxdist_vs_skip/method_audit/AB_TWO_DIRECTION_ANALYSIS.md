---
title: "A label-first vs B structure-first 兩方向方法分析 — 觀察/數據/圖表/方法細節驗證"
date: 2026-06-15
type: methodology analysis (兩方向方法理解 + 全基因組數據 + 圖表)
scope: B comprehensive — 全基因組 30490 region (HCC1395 paired, develop e4f4bbc)
data_sources: /tmp/dbeta_wg_abc/significance_summary.csv, gen_ab_direction_figure.py, 6-method audit /tmp/w4ukyquky.output
claim_levels: 觀察數據=L1(全基因組 CSV) / 方法性質判定=L1(源碼 audit) / 詮釋=L2
related: TWO_DIRECTIONS_FRAMEWORK / SIGNIFICANCE_INVENTORY / PART2_METHOD_VERIFICATION_REPORT
---

# A label-first vs B structure-first 兩方向方法分析

## 1. 概念定義（根本差異）

| | **A — label-first（標籤驅動）** | **B — structure-first（結構驅動）** |
|---|---|---|
| 起點 | **已知標籤**（HP/allele/sample/germline-carrier）| **甲基距離/結構**（無監督）|
| 問句 | 「這個標籤的 read，甲基**差不差**？」| 「甲基結構**對不對得上**某標籤？」|
| 方向 | 標籤 → 甲基差異 | 結構 → 標籤關聯 |
| 代表 | Δβ family（germline/somatic/subclone/alt/combo）、PerCpgAsm、LabelTest distance | HP-AUC、clustering+PERMANOVA、GlobalTest、HPFine |
| 輸出 | 方向+大小（Δβ）+ p | 結構強度（AUC/pseudo-F）+ 關聯（CramersV）|

---

## 2. A 方向 — label-first（觀察/數據/驗證）

| 方法 | 做什麼 | 全基因組數據 | 性質 |
|------|--------|------|------|
| **Δβ germline ASM** | normal HP1 vs HP2 mean β 差 + perm | **25% sig**（7733/30490）| ✅ 校準保守（stage3 驗 perm≈beta-reg≈MWU 88-94%）|
| **Δβ somatic residual** | (tumor−normal) HP 差殘差 + perm | ~10% sig（min-group guard 後）| ✅ #3 修法，對齊 gap2 上限 |
| **Δβ subclone** | tumor germline vs carrier + perm | HP1/HP2 各 ~8% | ✅ tiny-group guard |
| **Δβ alt 軸（A 新）** | tumor HP-family ALT vs REF | HP1 sig 2302/HP2 2180 | ✅ 互補 HP-tag |
| **PerCpgAsm Fisher** | per-CpG HP1 vs HP2 | median Frac_Sig 0.19 | 🟡 over-dispersion 膨脹（characterization-only）|

**A 性質**：**校準、保守、指方向**。每個 sig 都是「該標籤真有甲基差異」的嚴格證據。Δβ 已三方獨立交叉驗證（stage3）。

---

## 3. B 方向 — structure-first（觀察/數據/驗證）

> **🔑 B 有兩個子類，性質天差地別**：

### 3a. B-clean（HP-AUC）— 唯一乾淨校準的 B
- **做什麼**：P(dist(異HP對) > dist(同HP對))，NaN 跳。**無 clustering 中介**，直接測「距離有沒有還原 HP」。
- **數據**：normal **median 0.788 / 76% strong**；tumor **0.641 / 40% strong**。
- **性質**：✅ **乾淨**（不經 clustering 最佳化，無循環）。是**地基 sanity floor** + somatic 線索（tumor < normal）。

### 3b. B-clustering（PERMANOVA/GlobalTest/HPFine）— 過敏/循環，慎用
- **做什麼**：silhouette 最佳化 k 切群 → 在**同一距離資料**上測群是否顯著。
- **數據**：ClusterPERMANOVA **valid 中 100% sig** / HPFine 4-group **95% sig** / GlobalTest **74% sig**。
- **性質**：🔴 **過敏 + 循環**（先最佳化分離再測分離 = 幾乎必顯著）。**「有結構」幾乎恆真，不判別哪些 region 特殊**。CramersV median 0.561（effect 中等）但 p 不可信為發現門檻。

**B 性質**：HP-AUC 乾淨可信；clustering 系列過敏（確認結構存在，非發現工具）。

---

## 4. A vs B 關係（圖 `figures/ab_direction_comparison.png`）

**不是同一嚴格度**（全基因組）：
- A germline Δβ sig **25%** vs B GlobalTest sig **74%** → 一致率僅 **42%**，**B-only 16295（53%）= B 大幅 over-call**（GlobalP sig 但 Δβ 不 sig）。
- HPMergedSig（A-distance）vs ClusterPERMANOVA（B）**25% 不一致**（gap#6）。

**互補但不可互換**：
- **A 計算「哪個標籤真有差異」**（嚴格發現）。
- **B-clean(HP-AUC) 測「距離有沒有抓到標籤」**（地基）。
- **B-clustering 測「有沒有任何結構」**（過敏，慎當發現）。
- 個案實證（`figures/case_subclone_robust_chr1_58061351.png`）：germ32/car13 robust subclone Δβ **sig（A）但 CramersV=0（B-clustering 漏）** → **A 抓到 mean-shift，B-clustering 因非雙峰漏掉** = 互補的直接視覺證據。

---

## 5. 方法細節驗證 verdict（你要的「合理沒問題嗎」）

| 方向/方法 | verdict | 用途定位 |
|------|:---:|------|
| A Δβ family | ✅ 校準可信（stage3 三方驗）| **發現工具**：哪個標籤真有甲基差異 |
| A PerCpgAsm Fisher | 🟡 over-dispersion 膨脹 | characterization-only，勿單獨升 tier |
| B HP-AUC | ✅ 乾淨可信 | **地基 sanity** + somatic 線索（tumor<normal）|
| B clustering/PERMANOVA/GlobalTest/HPFine | 🟡 過敏/循環 | 確認「有結構」，**非發現門檻**；CramersV 可看 effect，p 慎用 |

**淨結論**：
1. **發現用 A（Δβ，校準）**；**地基/sanity 用 B-clean（HP-AUC）**；**B-clustering 當「結構存在性」確認，不當判別門檻**。
2. A 與 B **互補**（A 抓 mean-shift、B 抓多元結構/雙峰），但**嚴格度不同**，禁直接比 p 或混用門檻。
3. 全基因組地基穩固（HP-AUC normal 0.788）；tumor 弱化（0.641）是 somatic 線索但須 confound 排除（L2，非定論）。
