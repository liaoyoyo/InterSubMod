---
title: "最終判別假陰性 — 有結構卻沒被分類：理解 + 查驗計劃"
date: 2026-06-15
type: plan (理解 → 確認觀察 → 定位問題 → 查驗計劃；尚未實作/修正)
trigger: 用戶手動觀察『許多認為應有結構、應進有效分類的位點，卻沒被分類』
data_sources: /tmp/dbeta_wg_abc/significance_summary.csv (全基因組 30490), SignificanceAnalyzer.cpp:307-340, GlobalTest.cpp:142
claim_levels: pipeline=L1(源碼) / 量化=L1(全基因組 CSV) / 根因=L2(待個案查驗)
status: AWAITING UNDERSTANDING-CONFIRM + PLAN-APPROVAL
---

# 最終判別假陰性 — 有結構卻沒被分類

## 1. 最終判別 pipeline（一步步）

每個 region 最後得一個 `VerificationClass` ∈ {Strong, Subclone, Weak, Noise}，決策 = **2×2 of (label_sig × cluster_sig)**（`SignificanceAnalyzer.cpp:326-339`）：

| | cluster_sig=true | cluster_sig=false |
|---|---|---|
| **label_sig=true** | **Strong** | **Weak** |
| **label_sig=false** | **Subclone** | **Noise** |

- `label_sig` = LabelTest 距離 delta 顯著（A label-first 的距離版）
- `cluster_sig = passed_gate AND (GlobalTest_alt p≤0.05 OR GlobalTest_hpfamily p≤0.05)`
  - = 「**無監督 cluster 是否對得上 allele/HP 標籤**」（用 GlobalTest，**非** ClusterPERMANOVA 純結構）
  - `passed_gate = p_pass && v_pass`（`GlobalTest.cpp:142`）= GlobalTest p 顯著 **且** CramersV **reliable**（Cochran：<20% 格期望<5）

**全基因組分佈**：Strong 21015 / Weak 6122 / Noise 2058 / Subclone 1295；**PassedGating 75%**（25% 沒過 gate）。

## 2. 確認你的觀察（數據證實 = 真的）

「A 抓到真甲基差異(Δβ sig) 但最終沒進有效結構分類」：
- A-signal region **11736**，其中 **VC=Weak/Noise 2451（21%）**（Weak 2040 + Noise 411）
- 其中 **PassedGating=false 2162（18%）**

→ **你看到結構、A 也證實有甲基差異，但最終判別放進 Weak/Noise 的，確實約 1/5**。觀察成立。

## 3. 問題在哪（3 個機制，待個案查驗定位）

1. **🔴 gate 的 CramersV-reliable 要求（最可疑，~18%）**：`passed_gate = p_pass && v_pass`。v_pass = CramersV reliable（Cochran：格子期望count 不足→unreliable→**CramersV 被歸 0 + gate fail**）。**reads 少 / HP 不平衡 / cluster 小 的 region，即使有真結構，CramersV 不可信 → gate fail → cluster_sig=false → 最多 Weak**。
2. **🟡 label-anchoring（標籤錨定）**：cluster_sig 要求 cluster **對得上 allele/HP**。與標籤**正交的甲基結構**（如 mean-shift、或非 HP/allele 軸的 subclone）→ GlobalTest 不 sig → Weak/Noise。對應 A/B 互補（58061351：Δβ sig 但 cluster 不分）。
3. **🟡 cluster k 欠分辨 + 距離不分 mean-shift**：silhouette tie-break 偏小 k（audit 已揭）；distance-clustering 對 **mean-shift（非雙峰）不分群** → 你肉眼看到的群差異，clustering 沒形成 cluster → cluster_sig=false。

## 4. 如何查驗（計劃 — 尚未執行）

> 目標：對「A-sig 但 Weak/Noise」的假陰性，**逐案定位是哪個機制**，確認符合你的觀察，再決定修不修。

- **步驟1｜分層抽樣假陰性**：從 2451 個「A-sig 但 Weak/Noise」中，按根因候選分層抽代表案例：
  - (a) PassedGating=false 且 CramersV-unreliable（gate 卡）
  - (b) PassedGating=true 但 GlobalTest 不 sig（label-anchoring / cluster 不對標籤）
  - (c) cluster k=2 但 Δβ/heatmap 顯示 ≥3 群（k 欠分辨）
- **步驟2｜產假陰性個案圖**：每案甲基熱圖+距離熱圖+cluster 標註 + 標明「為何沒被分類」（gate fail 原因 / GlobalTest p / CramersV reliable? / NGroups）→ **你肉眼確認是否符合觀察**。
- **步驟3｜量化各機制佔比**：2451 假陰性中 gate-卡 vs label-anchoring vs k-欠分辨 各多少 → 定位主因。
- **步驟4｜決定修法**（理解後）：
  - 若主因=gate CramersV-reliable 過嚴 → 放寬 / 改用 ClusterPERMANOVA(純結構) 補一條「structure-present 但 label-orthogonal」分類
  - 若主因=label-anchoring → 加一個「label-free 結構」判別軸（純 ClusterPERMANOVA / HP-AUC 不經 GlobalTest）
  - 若主因=k 欠分辨 → 調 silhouette tie-break / max_k

## 5. 機制量化結果（2026-06-15 完成）

全基因組「A-sig 但 Weak/Noise」假陰性 **2451**：

| 機制 | 數量 | 佔比 | 本質 |
|------|---|---|------|
| **p_fail**（GlobalP>0.05，cluster 對不上 allele/HP）| 2307 | **94%** | A 抓到甲基差異但離散 cluster-cut 對不上標籤；HP-AUC median **0.749**（結構其實在）|
| **v_fail**（GlobalP≤0.05 但被 CramersV-reliability gate 殺）| 144 | 6% | 關聯顯著但 Cochran gate 過嚴（NumReads med 130 不小）|

**🔴 核心根因（源碼確認）**：
1. `label_sig = label_result.any_significant`（`SignificanceAnalyzer.cpp:282`）= **LabelTest 距離版**；**Δβ 模組在 SignificanceAnalyzer 0 個引用 → 校準的 Δβ 完全不參與最終 verdict**。
2. `cluster_sig` 錨定**離散 cluster × 標籤的 GlobalTest 配對**（最脆弱一步）；即使 HP-AUC=0.969（距離極強抓到 HP），離散切群對不上 HP → GlobalTest 失敗 → 判 Weak。
3. 個案確認（圖 `figures/case_pfail_*`/`case_vfail_*`）：6202793 HP-AUC 0.969 / 248732962 0.952 / 182241036 GlobalP=0 被 reliability gate 殺 — 距離熱圖 HP block 清晰但 VC=Weak。

**淨結論**：判別**沒用上更好的訊號**（連續 HP-AUC + 校準 Δβ），只認脆弱的離散 cluster×標籤 → 你看到的結構被丟。

## 6. 修法方案（待決策，皆改 VerificationClass=golden 欄，需更新 golden + 驗證）

- **方案 A｜HP-AUC 接進 cluster_sig**（解 p_fail 主體）
  - `cluster_sig = 原邏輯 OR (HP_AUC_Normal>=0.7 或 _All>=0.7)`：距離有 HP 結構就認，不必離散 cluster 對上標籤。
  - 影響：~94% p_fail（高 HP-AUC 者）改判 Strong/Weak→更正確。最直接對應根因2。

- **方案 B｜Δβ 接進 label_sig**（解「Δβ 不參與判別」）
  - `label_sig = 原邏輯 OR (germline ASM Δβ sig 或 subclone Δβ sig)`：校準的 mean-shift 進 verdict。
  - 影響：A-only mean-shift 不再被埋進 Noise；對應根因1。

- **方案 C｜放寬 v_fail 的 reliability gate**（解 6%）
  - GlobalP 顯著但 CramersV unreliable 時，改用 raw 關聯（不歸 0）或軟化 Cochran 門檻。

- **方案 D｜新增 VerificationClass「StructureNoLabel」**（保守，不改既有 4 類語意）
  - 對「HP-AUC/PERMANOVA 有結構但 GlobalTest 對不上標籤」獨立標一類，不混進 Weak/Noise。

> **建議**：A+B（接 HP-AUC + Δβ）最對根因，但改 VerificationClass → 須 `--update-golden` + 全基因組重驗分佈。建議**先 pilot**（Python 重算新判別邏輯，看 2451 假陰性有多少被救 + 是否引入假陽性）→ 確認淨改善才 cpp-change。
> **決策**：選方案（A/B/C/D 組合）+ pilot-first vs 直接 cpp-change。
