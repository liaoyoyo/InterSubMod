<!--
建立時間: 2026-06-18
類型: 分類驗證 spec（design；落地分件走 /cpp-change 或 offline 分析）
狀態: approved_design（用戶 2026-06-18 鎖定確認）
build_branch: research/subclonal-reconstruction-202606（doc）；C++ 件目標 = feat/summary-nreadsvalid（含 within-HP/Bootstrap/PERMDISP/C2-ARI）
depends: C2（已落地 commit 64e68d9 feat/cluster-label-ari）、C1（待做，spec=20260618_within_hp_bootstrap_and_cluster_label_ari_cpp_change_spec_01.md §7.3）
data_sources: docs/methodology/_assets/20260616_fn_verdict_readback/HCC1395_C2_ari_typology.json
-->

# 分類驗證 spec — cluster×label 正確分類（subclone 為下游）

> **核心紀律（用戶定）**：先建立**可驗證、正確的分類**；subclone 判讀是**下游**；本階段**只檢驗與分類，不做 germline/somatic 詮釋**。

## 0. 目的 / 範圍 / 非目標

- **目的**：對每個 region，**清楚地檢驗並分類**「甲基無監督結構 ↔ 原始標籤」的對應關係，且分成 **tumor-only** 與 **tumor+normal** 兩層、保留**標籤判讀**與**結構確認**兩種讀法。
- **非目標（明確不做）**：subclone calling、germline/somatic 判定、CCF/tree 重建。這些建在本分類**之上**，屬下游。
- **scope**：先單樣本 HCC1395（⭐3），spec 通過後同腳本 scope-out 6 樣本。

## 1. 六層模型 + 三鎖定決定

```
①甲基無監督 UPGMA 建距離+分群   ← 切 k 由 C1 bootstrap-stability 給可信答案
②PERMANOVA + PERMDISP 確認對齊哪種標籤（位置差 vs 離散度假象）  ← 已在 binary（feat D2）
③標籤判讀(label-first) + 結構確認(structure-first) 兩讀法並存
④HP-fine sub-label 對應 typology（含「一群含多標籤」flag + 標籤組合對齊測試）  ← 新增 C3
⑤tumor-only 分類 + tumor+normal 分類 + 兩層差異「事實表」（不判別）  ← 新增（read-set 切換 + offline）
⑥subclone 分析 / 其他目的  ← 下游，建其上（不在本 spec）
```

**三鎖定決定（用戶 2026-06-18）**：
1. **tumor-only 獨立重跑一套**：只用 tumor reads 重建距離矩陣 + UPGMA + PERMANOVA/PERMDISP（與 paired 各一套輸出）。
2. **sub-label（HP-fine）對應 + 組合測試 + flag**（見 §4），要能**正確分析+標記**「一大結構含兩標籤」與「兩標籤組合是否對齊結構」。
3. **tumor-only vs tumor+normal 只記事實、不判別**：對每軸記兩 read-set 各自顯著性，直接寫「tumor 顯著 / (加 normal) 不顯著」等事實，**不下 germline/somatic 結論**，意義留下游。

## 2. 與現況對照 + 實作分件

| 要件 | 狀態 | 落地分件 |
|---|---|---|
| ②PERMANOVA + PERMDISP | ✅ binary（feat D2 commit 74eb2e2） | 消費既有 |
| HP-merged / allele ARI | ✅ C2（commit 64e68d9，paired） | 消費既有 |
| ①切 k bootstrap 穩定 | ❌ | **C1**（前置依賴；spec §7.3） |
| ④HP-fine sub-label 對應 + 組合 + flag | ❌ | **C3**（C++ emit，本 spec §4） |
| ⑤tumor-only 獨立重跑 | ❌ | **P-readset**（執行層，§5） |
| ⑤兩層事實表 + 差異 | ❌ | **P-facttable**（offline，§6） |

> ⚠ **依賴順序**：④的「per-cluster 子標籤組成」建在 cluster 指派上 → **cluster 須先由 C1 穩定**，否則組成是「過切假群」的組成。**建議 C1 先於 C3**（或 C3 先輸出、但 k 標 provisional）。

## 3. 共用前置：PERMANOVA + PERMDISP 雙讀法（③，消費既有）

每個 region、每個候選標籤分組，輸出兩讀法：
- **label-first**：`run_grouped_structure(partition)` → PERMANOVA F/p（位置）+ **PERMDISP p（離散度）** → 「位置差顯著且離散度乾淨」才算真對齊（非假象）。
- **structure-first**：UPGMA cluster（k 由 C1）+ **ARI(cluster, partition)**（C2/C3）。
- 「對齊」= 兩讀法一致（ARI 高 **且** PERMANOVA 位置顯著 **且** PERMDISP 不警告）。

## 4. C3 — HP-fine sub-label 對應分析（C++ emit，新增）

HP-fine 原始子標籤 = `1` / `1-1` / `2` / `2-1`（germline 與其 somatic-carrier 子型）。

### 4.1 候選標籤分組（對每個都算 ARI + PERMANOVA + PERMDISP）
| 分組 | 內容 | 問的問題 |
|---|---|---|
| **Fine-4** | {1},{1-1},{2},{2-1} | 甲基結構對齊到完整 4 子標籤？ |
| **HP-merged** | {1,1-1},{2,2-1} | 對齊 germline 單倍型？（= C2 的 HP，可重用） |
| **Somatic-vs-germline** | {1-1,2-1} vs {1,2} | somatic-carrier 身分是否對齊甲基結構？ |
| Allele | ALT vs REF | （= C2 的 allele，可重用） |

→ emit 每個分組的 `ARI_*`、`Permanova*_P`、`Dispersion*_Warn`；**標出哪個分組對齊最好**（最高 ARI ∩ 位置顯著 ∩ 離散度乾淨）。

### 4.2 per-cluster 子標籤組成 + 兩個 flag（2.1 的核心）
對每個 region 的幾何 cluster：算每 cluster 內各 HP-fine 子標籤的 read 數 →
- 🚩 **`ClusterMultiSublabel`**：存在一個 cluster 內含 ≥2 個 HP-fine 子標籤（且每個 ≥ minN）→「**一大結構含多標籤**」（一結構→多子標籤）。
- 🚩 **`SublabelMultiCluster`**：存在一個 HP-fine 子標籤散落 ≥2 個 cluster → 「**多小結構→一子標籤**」。
- emit：`NCluster`(=optimal_k)、`ClusterPurity_Mean`（各 cluster 主子標籤佔比平均）、上兩 flag。

### 4.3 「兩標籤組合是否對齊」的判讀（2.1）
由 §4.1 的分組 ARI 比較直接得：
- 若 **Somatic-vs-germline** 分組 ARI 高 → 甲基結構是按「有無 somatic-carrier」分（兩標籤組合 {1-1,2-1} 對齊）。
- 若 **HP-merged** 高但 Fine-4 不再更高 → 結構按 germline 單倍型分，**somatic 子型在甲基上與母單倍型不可分**（→ `ClusterMultiSublabel` 會 true）。
- 全部**只標記分組對齊狀態，不詮釋 subclone**。

### 4.4 C3 emit 欄（significance_summary.csv 純加）
`ARI_Cluster_HPFine, Permanova_HPFine_P, Dispersion_HPFine_Warn, ARI_Cluster_SomaticVsGermline, Permanova_SomaticVsGermline_P, NCluster, ClusterPurity_Mean, ClusterMultiSublabel, SublabelMultiCluster`（observation-only、不進 verdict、不改行為；沿用 C2 模式 + sentinel −2.0）。

### 4.5 C3 驗收
- ctest 全綠 + 新增單元測試（純對齊→ARI=1 / 一群含兩子標籤→`ClusterMultiSublabel`=true / somatic 對齊→Somatic-vs-germline ARI 高）。
- chr2:18M HCC1395（已知 6 sSNV）：`ClusterMultiSublabel`/組合對齊應與該位點實測一致。

## 5. P-readset — tumor-only 獨立重跑（執行層）

- **做法（先驗證可行）**：跑 binary **不帶 `--normal-bam`**（只 `--tumor-bam`）→ 距離/分群/HP/allele/HP-fine **全在 tumor reads**；normal 依賴欄（SampleASM/Normal_HP）自然 NA。**若 binary 允許省略 `--normal-bam` 即零改碼**；否則加 `--tumor-only-clustering` 小旗標（restrict 距離建構的 read-set）。← spec checklist：先確認 `--normal-bam` 是否 optional。
- 輸出：`<sample>_tumor_only/significance_summary.csv`（與 paired 同 schema，含 C2/C3 欄）。
- 與 paired 對照：同一批 region、同 join key（anchor）。

## 6. P-facttable — tumor-only vs tumor+normal 事實表（offline，不判別）

offline 腳本 join `tumor_only` 與 `paired` 兩 CSV（on anchor），對每軸（HP/allele/HP-fine/somatic-vs-germline）輸出**事實**：
| 欄 | 意義（純事實） |
|---|---|
| `<axis>_TumorOnly_Sig` | tumor-only PERMANOVA 位置顯著 ∩ PERMDISP 乾淨？ |
| `<axis>_Paired_Sig` | tumor+normal 同上？ |
| `<axis>_ARI_TumorOnly` / `_Paired` | 兩 read-set 的 ARI |
| `<axis>_Change` | 事實標籤：`tumor_sig_paired_notsig` / `both_sig` / `both_notsig` / `paired_sig_tumor_notsig` |
- 🔴 **`_Change` 只是事實字串，不映射 germline/somatic**（用戶定；意義留下游 ⑥）。
- 產彙總：各 `_Change` 類別佔比 + ARI 分布（tumor-only vs paired 對照；呼應 C2 paired median ARI 0.15）。

## 7. 落地順序（建議）

1. **C1**（within-HP bootstrap k 穩定）— 讓 cluster 可信，是 C3 per-cluster 組成的前置。
2. **C3**（HP-fine 對應 + flag + 組合，C++ emit）。
3. **P-readset**（tumor-only 重跑；先驗 `--normal-bam` optional）。
4. **P-facttable**（offline 兩層事實表）。
- 每個 C++ 件走 `/cpp-change` 6 步（pre-commit 編譯 Hard Gate）；offline 件 §13.0 run→寫→讀→報告。

## 8. 總驗收（分類正確性，非 subclone）

- [ ] 每軸 ARI + PERMANOVA + PERMDISP 三者一致定義「對齊」（§3）。
- [ ] HP-fine 4 分組 + per-cluster 組成 + 兩 flag 正確（C3 單元測試 + chr2:18M 對照）。
- [ ] tumor-only 與 paired 各一套輸出、同 schema、可 join。
- [ ] 事實表 `_Change` 四類齊全且**無 germline/somatic 詮釋**。
- [ ] 全程 observation-only：**不改 verdict、F1/分類分布不回退**（C2 已立此先例）。
- [ ] 6 樣本 scope-out（per-sample tally，禁 pool）。

## 9. 用戶決策

**先做哪件？** [ ] C1（k 穩定，C3 前置） / [ ] C3（HP-fine 對應，可先輸出 k 標 provisional） / [ ] P-readset（tumor-only 重跑先看差多少） / [ ] 全部依 §7 順序
**備註 / 調整**：
