<!--
建立時間: 2026-06-22
問題類型: 方法規格（phylo-v3.1 切群標籤 — C++ 內建前的權威 spec）
影響 track: subclonal reconstruction pilot（tumor-only 切群層）
狀態: spec_pending_user_confirm（定義+方法確認後才動 C++）
data_sources: _assets/20260618_subcluster_pilot/phylo_pilot_stats.json,_assets/20260618_subcluster_pilot/phylo_sensitivity.json,_assets/20260618_subcluster_pilot/phylo_v3_validate.json,_assets/20260618_subcluster_pilot/phylo_v3_wg_summary.json
build_branch: feat/summary-nreadsvalid
-->

# phylo-v3.1 切群標籤 — 方法規格 + 完整統計（C++ 內建前定稿）

> **此文件職責**：在改 C++ 前，把「判別 + 驗證的邏輯順序、門檻定義（數據依據）、完整統計、邊界案例」一次定清楚。修正歷程（double-dip → FM1）見 `20260622_phylo_subcluster_labeling_doubledip_correction_01.md`。

## §1 判別 + 驗證 邏輯順序（決策序列）

每個位點（一筆資料，不分 TP/FP）依序：

```
[輸入] read×CpG 甲基矩陣 M + read×read BERNOULLI 距離 D（C++ binary 已輸出）
  │
① peel        : 移除 NaN-距離無法完整配對的 read（filter_reads_for_complete_matrix）
  │            → 不足 2·MINSZ(=6) read → 終止「不可分群」
② UPGMA 樹    : 對 D 建 average-linkage 樹（HierarchicalClustering，已存在）
  │
③ 遞迴判別（rec，從根開始）對每個 node：
  │   3a. descend : 沿樹剝離 (tiny,big) caterpillar 單離群，直到「兩子群皆 ≥MINSZ」的平衡節點
  │                （修 FM1：避免單一離群讓整群判 1 群）
  │   3b. 分離比 r = between(C1,C2)/within  ── r < SEP_MIN(1.3) → 此 node 不切（pre-filter）
  │   3c. null 檢定 : 對該子群 read 做「逐 CpG 欄內置換 → 重新 UPGMA → 取其 root 分裂比」×RNULL(40)
  │                  r > null 的 95 百分位 → 分裂顯著（specificity 核心閘）
  │   3d. 物有所值 : 僅在 3c 通過時才把 descend 剝離的 read 標 outlier；不通過 → 整 node 保留不侵蝕
  │   3e. 階層標籤 : 通過 → 較多 read 子群編號較小（1 / 1-1 / 1-2 / 2…），遞迴兩子群
  │
④ 後過濾      : 終端標籤 read 數 <MINSZ → 標 outlier
⑤ 對齊歸因    : 各群 dominant hp + allele(REF/ALT)；群間 allele 不同 → allele 軸(cis-ASM)；
  │             hp 不同 → hp 軸；同 allele+hp → within-germline(subclone 候選)
⑥ 三態輸出    : (a)單群=無離散結構 (b)多群+對齊 germline 軸=cis-ASM (c)多群+within-germline=subclone 候選
```

**驗證一個切割是否真實（§回答「如何驗證判別」）= 四者一致**：
clusterboot Jaccard 穩定（Hennig >0.75）+ 對齊 CramérV(hp/allele) + 平衡 r ≥ 噪音 null95 + 全域純噪音 FP ≈0%。

## §2 門檻定義 + 數據依據（敏感度掃描 `phylo_sensitivity.py`）

| 門檻 | 值 | 類型 | 數據依據 |
|---|---|---|---|
| **null95（per-subgroup 重分群）** | 95pct | ★ data-derived 核心 | **唯一控 FP 的閘**（SEP_MIN 1.0-1.5 噪音 FP 全 0% → FP 來自此 null）；百分位控敏感度（90→12 / 95→10 / 99→8 多群，都 0% FP）|
| SEP_MIN | 1.3 | 約定·安全平台 | 1.0-1.3 全給 10 多群+0%FP（冗餘於 null95）；1.4↑ 才砍真群 → 可降無害 |
| MINSZ | 3 | 約定·穩健 | 3/4/5 結果不變（都 10 多群）→ 非關鍵；=最小可估群統計量 |
| RNULL | 40 | data-informed | 穩定位點 12 已收斂；borderline 位點任何 RNULL 都翻 → 40 對穩定足夠，borderline 另以多-seed 標 ambiguous |
| Jaccard >0.75 | — | 文獻 Hennig 2007 | cluster 穩定性標準 |
| CramérV≥0.3 / Cochran e≥5 | — | 文獻 | 效應量小門檻 / χ² 可靠性 |

## §3 完整統計（pilot 30 位點，HCC1395 tumor-only chr20-22，**curated 非隨機**）

> 全基因組（30490 位點隨機真值）統計見 §5；本節為 curated pilot（過取樣有結構者）。

**A. 結構/無結構**：有結構(≥2群) **10 (33.3%)** / 無結構(1群) **20 (66.7%)**
**B. 群數分佈**：1群 20(66.7%) · 2群 7(23.3%) · 3群 2(6.7%) · 4群 1(3.3%)
**C. 分層深度**（多群 10）：深度1(如 1-1) 5 位點 · 深度2(如 1-1-1) 5 位點
**D. 對齊型態**（多群 10）：**allele 軸 REF/ALT (cis-ASM) 10 (100%)** · hp 軸 0 · within-germline(subclone 候選) **0**
**E. 4-gate fine_class × v3 群數 交叉**：

| 4-gate class | 位點 | v3 多群 | v3 單群 | 解讀 |
|---|---|---|---|---|
| CONFIRMED | 6 | **5** | 1 | v3 多切（real+aligned）|
| NEAR_CONFIRMED | 6 | 2 | 4 | 部分 |
| REAL_NOVEL | 6 | 3 | 3 | 部分 |
| REAL_DIFFUSE | 6 | **0** | 6 | **v3 全收單群** — diffuse/gradient 非離散群 |
| NO_CLEAR | 6 | 0 | 6 | 一致無結構 |

→ **v3(10) vs 4-gate(24) 差距 = REAL_DIFFUSE**：4-gate「有無可重現結構」把 gradient 算結構；v3「有無離散群」正確拒 gradient。

## §4 邊界案例（需特別觀察）

| ID | 案例 | 位點 | 為何要看 |
|---|---|---|---|
| **F1** | verdict 不穩 | chr20_42981498 | RNULL=80 仍 1/2/3 翻動 = **資料在決策邊界** → 該標 ambiguous，非強判 |
| **F3** | 深層階層 depth≥2 | chr20_56191455 / chr22_30454004 / chr20_21855867 / chr20_30274614 / chr20_24575172（5）| within-allele 子結構（如 chr22_30454004 三個 ALT 子群）= subclone 候選 vs cis-ASM 漸層，需肉眼+normal cis-control |
| **F4** | 低 CpG（≤40，FP 略升區）| chr21_40743336 (C=38) | 校準顯示 C≤40 FP 略升；此例 clean REF/ALT 對齊支持為真但標 borderline |
| **F5** | 大 n 高信心多群 | chr22_44981696 / chr20_56191455 / chr22_30454004 / chr22_26939195 / chr20_24575172（5）| n≥60 + 強對齊 = 最可信，可當 canonical |

> F2（locus 層 within-germline）= 0：所有多群位點都跨 REF/ALT（cis-ASM）；within-germline 子結構只在 F3 深層階層內。

**全基因組邊界盤點**（`phylo_edge_analysis.py`，TP 30077/FP 4659）：

| ID | 案例 | 量 | 判定 |
|---|---|---|---|
| **E1** | 0群退化（全 read 變離群）| TP31/FP26（n 中位 9）| 小 n caterpillar 無凝聚群 → 合理；**應報「無結構」非 0-state**（報告修）|
| **E2** | 高群數 ≥5 | TP180（n 中位 99）| **aligned 98.9%** → 真 cis-ASM 豐富分層，**非過切** ✓ |
| **E3** | 深層階層 depth≥4 | TP129（n 中位 90）| 大 n 合理；depth-6(1) 須看 |
| **E4** | FP 帶 unaligned 結構 | FP402(8.6%)，V_allele=0 | **FP/TP=2.0×** → 噪音非 somatic（specificity 確認，非 subclone）✓ |
| **E5** | 低 CpG ≤20 帶結構 | TP 僅 3 | genome-wide 極罕見（±5000 窗）→ 非問題 ✓ |
| **E6** | 大 n ≥150 | TP52（aligned 90%）| 最可信 canonical ✓ |

> **待處理 2 項**（C++ 內建時納入）：① E1 0群 → 報「無結構」；② F1 instability flag（多 seed 一致性 → unstable 標 ambiguous）。其餘邊界合理。代表圖見邊界工作站 `phylo_edge_verify_workstation.standalone.html`（repo 外 `/big7_disk/liaoyoyo2001/ism_html_review/20260622_phylo_edge_verify_workstation.html`）。

## §5 vs 之前「~30% 有結構」的比較（全基因組）

**「30%」來源**：`cluster_redesign_wg`(4-gate) 全基因組 CONFIRMED+REAL_NOVEL = (6198+3206)/30077 = **31.3%**；與 decisionflow `split_align_rate 31.34%` 一致 = 「有清楚離散/對齊結構」比例。

**v3 全基因組（`phylo_v3_wg.py`，RNULL=25，elapsed 1962s）**：

| | TP (n=30077) | FP (n=4659) | TP/FP 富集 |
|---|---|---|---|
| **有結構(≥2群)** | **24.3%** | 24.06% | 1.0×（總量不判別）|
| ├ aligned (cis-ASM) | **19.98%** | 15.43% | **1.30×（TP 偏多）** |
| └ unaligned (subclone 候選) | **4.32%** | **8.63%** | **0.50×（FP 偏多！）** |
| 無結構(1群) | 75.7% | 75.94% | — |

群數分佈 TP：1群 22738 · 2群 4204 · 3群 2249 · 4群 675 · 5群 135 · 6群 38。深度：depth1 3455 · depth2 2910 · depth3 815 · depth4 113。

**差別與原因（全基因組定案）**：
1. **v3 24.3% < 4-gate 31.3%（低 7pt）**：v3 較嚴（per-subgroup null + 離散群要求 + descend），把 4-gate「REAL_NOVEL(10.7%) 含 double-dip 假象」與部分 diffuse 拒掉；FM1 修救回一些但淨值較嚴。**「30%」鬆在它把 gradient/big_gap 都算結構**。
2. **🔴 組成才是關鍵（修正後的乾淨結論）**：
   - **aligned cis-ASM 20.0%（TP 富集 1.30× over FP）** = 唯一 somatic-偏向的結構訊號。
   - **unaligned「subclone 候選」4.32%，但 FP 更多(8.63%)=0.50× → 非 somatic 特異** → **單樣本無法當 subclone**（與舊結論一致，但 v3 用對齊軸更乾淨地證實）。
3. **雙向修正抵銷**：double-dip 修移除 v1 假 subclone（降）+ FM1 修救回單離群遮蔽真切割（升）→ 淨 24.3% 比 4-gate 31.3% 低但**每一個都過 null 閘、可信**。
4. **與舊「30%」最大概念差**：舊 31.3% 是「有清楚/對齊結構」混合；v3 拆成 **20% cis-ASM(偏 TP) + 4.3% 未對齊(偏 FP 非 subclone)**，把「有結構」誠實分解為「絕大多數是 cis-ASM、subclone 候選非特異」。

## §6 C++ 內建規格（**待 §1-§5 用戶確認後才實作**）

> 此節定義 C++ 該做什麼，**非實作**。C++ 是 Hard Gate，須走 `/methodology-audit` → `/cpp-change` 6 步。

**程式流程**（§1 邏輯順序的 C++ 對應；~70% 零件已存在）：
- 既有：`DistanceMatrix`(BERNOULLI) / `HierarchicalClustering`(UPGMA) / `Tree`+`MergeRecord` / `StructureTest::filter_reads_for_complete_matrix`(=peel) / `run_permanova`(置換範本) / `Bootstrap`(clusterboot 槽)。
- 新增：① descend-quarantine（樹遞迴剝單離群至平衡節點）② per-CpG 欄內置換 null（重算 BERNOULLI，reuse DistanceMatrix）③ 遞迴分裂 + null95 閘 + 物有所值 quarantine ④ 階層標籤產生 ⑤ 對齊歸因 + 三態。
- 取代：現 `TreeCutter::find_optimal_clusters`(silhouette best_k = 證實 75-91% 噪音 FP) → phylo-v3 遞迴。

**v4.1 改善（用戶回饋 + 對抗再驗證落地，`phylo_v4.py`）**：
- **① other 群**：殘留 outlier ≥MINSZ → 標 `other`（residual 紀錄）+ hidden_het 警告（n_other/n>0.3）+ merge-back（other_recovered_fine）。chr8 40 殘離群→記錄。
- **② instability redesign（modal）**：K=10 seeds → **modal_ng（多數決 robust verdict，取代脆弱單 seed）+ modal_frac（信心）**；`unstable = modal_frac<0.7`（只真分裂，非單 seed 偏離）。修對抗驗證的 K-敏感 + 用戶 chr20:42981498 = **modal 3@0.8**（=「應 3-4 群」）。
- **③ coarse/fine 兩層**：coarse=null95 modal / fine=null90（候選低信心，標 `_confidence:low`）。chr4 modal 4·5（=「大群可分更多」）。

**對抗再驗證（workflow，PROBE→修）+ 全基因組 v4.1 驗證（`phylo_v4_wg.py`，RNULL=40 K=5）**：
- 5 必修 #1(modal)#3(merge-back/warn)#4(誠實標籤)#2(全基因組) **全 PASS**；#5 C++ RNG harness = cpp-change 驗收。
- 全基因組 TP n=30077：**structure 23.83%（確認 v3.1 24.3% 不漂移）· aligned 19.59%(FP 15.15%=1.29×) · unaligned 4.24%(FP 8.24%=0.51×非subclone) · fine_multi 5.29% · other TP14%/FP26.6% · unstable 5.61% · hidden_het TP0.52%/FP2.38%(4.6×)**。FP 更噪音（other/hidden_het 偏多）再證 specificity。
- 驗證儀表板 `phylo_v41_final_dashboard.standalone.html`（repo 外 `/big7_disk/liaoyoyo2001/ism_html_review/20260622_phylo_v41_final_dashboard.html`）。

**輸出格式**（讓 Python 只讀畫圖，含 v4 欄）：
- `phylo_groups.tsv`：`read_id / coarse_label / fine_label / is_other / is_outlier / group_n / seed / sep_min / rnull`
- `phylo_groups_summary.json`：per-group `{n, dominant_hp, dominant_allele, align_type, parent_label}` + locus `{coarse_ng, fine_ng, n_other, maxdepth, structure_state, unstable, ng_range}`
- **不塞進 reads.tsv**（穩定 schema vs 衍生結果分層）。

> C++ 內建定稿 = v4（coarse/fine + other + instability），非僅 v3。重驗工作站 `phylo_v4_verify_workstation.standalone.html`（repo 外 `/big7_disk/liaoyoyo2001/ism_html_review/20260622_phylo_v4_verify_workstation.html`）。

## §7 驗證表

| # | 數字 | 來源 | L 級 | 狀態 |
|---|---|---|---|---|
| 1 | pilot 多群 10/30 (33.3%) | `phylo_pilot_stats.json` | L1 | ✓ |
| 2 | 多群 100% allele 軸 cis-ASM | `phylo_pilot_stats.json` atype_dist | L1 | ✓ |
| 3 | REAL_DIFFUSE v3 全收單群 (0/6) | `phylo_pilot_stats.json` 交叉 | L1 | ✓ |
| 4 | SEP_MIN 1.0-1.5 噪音 FP 全 0% | `phylo_sensitivity.json` | L2(模擬) | ✓ |
| 5 | chr20_42981498 RNULL=80 仍翻 | `phylo_sensitivity.json` rnull | L2 | ✓ |
| 6 | 「30%」=4-gate CONFIRMED+REAL_NOVEL 31.3% | `cluster_redesign_wg_summary.json` | L1 | ✓ |
| 7 | chr22 v3 structure 27.8%/aligned 26.96% | `phylo_v3_wg`(chr22 smoke) | L1 | ✓ |
| 8 | 全基因組 v3 統計 | `phylo_v3_wg_summary.json` | L1 | {{背景跑中}} |

**限制**：單樣本 HCC1395 tumor-only、pilot 30 為 curated（非隨機，比例不代表全基因組，看 §5）、⭐2-3 characterization 非 subclone 確認。
