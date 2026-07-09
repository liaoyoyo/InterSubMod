<!--
建立時間: 2026-07-09
類型: 缺口量化 + 小規模驗證 — E_subcube_recovered(partial-read 救回)區未真正建/枚舉 Steiner 樹
狀態: finding(使用者直覺 catch → 數據確認);可一鍵重驗
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/candidate_trees.json
provenance: 數字由 scripts/subcube_recovery_gap_check.py 於 2026-07-09 對 topology_per_region.json 重算(edges/n_sSNV/subread_groups 直接欄位)。partial flag: 單樣本 HCC1395。
-->

# 缺口量化:partial-read 救回區未真正建/枚舉 Steiner 樹(HCC1395)

> **使用者直覺 catch**:沒有 read 穿 A-B-C 時,正解 = 把整區當「一個 n-cube + 子面(pairwise)group 約束 + 隱藏 ABC 節點」,枚舉所有相容 Steiner 樹;而非分開拼 pairwise。**懷疑現行碼沒真做。** → 本文以數據確認:**懷疑成立。**

## §1 假設
E_subcube_recovered(gap#1 partial-read 救回)區只被「辨識 + 標 tree_shape」,**未跑 sub-face group-Steiner 建樹/枚舉**(edges 應全空)。

## §2 量化(全 = 直接欄位讀,verify 腳本重算)〔L1〕
| 指標 | 值 |
|---|---|
| E_subcube_recovered 區總數 | **2,403** |
| **有非空 edges(真建了樹)** | **0 / 2,403** 🔴 |
| n=2(單對·無需枚舉) | 51(2.1%) |
| **n≥3(需真 Steiner)** | **2,352(97.9%)** |
| 子標 pairwise樹 / 欠定 / 含衝突 | 1,561 / 676 / 166 |
| 乾淨 n≥3(veclen==n) | 2,148 |
| **其中含未觀測對(該枚舉多候選/需隱藏節點)** | **2,148 = 100.0%** 🔴 |

- **enumerate_candidate_trees 只枚舉 `ambiguous` 區;E_subcube_recovered 只被計數不被枚舉**(`enumerate_candidate_trees.py:331`)→ 這 2,403 區從頭到尾沒建/枚舉過樹。

## §3 小規模驗證(5 個乾淨 n=3 recovered 區逐個追)〔L1〕
| 區 | tree_shape 標籤 | 觀測到的對 | 應輸出 | 現況 |
|---|---|---|---|---|
| chr4:57969733-57996887 | full_tree | **1/3**(只對(1,2)) | 需枚舉(缺 2 對→多候選/隱藏節點) | edges=[] |
| chr1:24463672-24515927 | linear_nested | 2/3 | 需枚舉(缺 1 對) | edges=[] |
| chr1:38736335-38763132 | linear_nested | 2/3 | 需枚舉(缺 1 對) | edges=[] |
| chr1:53135760-53164656 | sibling_only | 2/3(皆稀疏) | 需枚舉 | edges=[] |
| chr1:53414848-53440268 | full_tree | **1/3** | 需枚舉(缺 2 對) | edges=[] |

🔴 **關鍵**:chr4 只觀測到 1/3 對就標 `full_tree` → **tree_shape 標籤在資料不完整時就掛上,且未真的建/枚舉樹**。

## §4 判斷狀況(嚴重度)
- **完備性:MEDIUM-HIGH 缺口** — 幾乎全部 recovered 區(乾淨 n≥3 的 100%)該輸出 co-optimal 樹集(含隱藏節點),卻只拿到一個 tree_shape 標籤 + edges=[]。使用者「算法有問題」的直覺**數據上成立**。
- **對 headline 結論:LOW 影響** — recovered 區**不計入 A_determined(1,764)**,主「已定樹」結論不受污染;誠實性(不僭稱 A_determined)有保住。
- 🔴 **需更正的措辭**:「partial-read 救回 2,403 區」的精確語意 = **『辨識出 2,403 區有 partial 共現底料』,不是『建出/枚舉出 2,403 棵樹』**。之前 explainer §10.2 說「新增部分樹結構(上界)」仍偏樂觀 → 更正為「**辨識底料,樹未建/未枚舉**」。

## §5 結果是否合理且對應假設驗證
✅ **合理且強力對應假設**:假設「recovered 未真建/枚舉」→ 數據 edges=0/2403、100% 乾淨 n≥3 含未觀測對 → 完全吻合。這**不是**演算法輸出「錯的樹」(它根本沒輸出樹),而是**「辨識層 ≠ 求解層」的實作缺口**:gap#1 做到「認出這些區 + 標記」,**沒做到「用子面 group 建/枚舉 Steiner 樹」**。

## §6 修法(使用者指出的正解)
把 recovered 區(尤其乾淨 n≥3 的 2,148)**路由進真正的 sub-face group-Steiner 枚舉**:以觀測到的 pairwise 子面為 group 約束,建隱藏 ABC 節點,**枚舉所有相容 Steiner 樹**(未觀測對 → 保持自由 → 輸出多候選),取代現在的「標 tree_shape + edges=[]」。⚠ 需先解 §另述的 Gap B(CN-gain 共享位點多拷貝)+ 那 126 個 veclen>n 的資料不一致。

## §7 一鍵重驗
```bash
cd InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts
python3 subcube_recovery_gap_check.py
```
只讀不寫,重算 edges=0/2403、n≥3=97.9%、乾淨 n≥3 含未觀測對=100%。2026-07-09 實跑確認。

**關聯**:explainer `20260709_pipeline_explainer_and_verification_01.md`(§10.2 已更正);形式化 `20260704_formal_problem_statement_topology_from_cooccurrence_01.md`;commit `68ced6f`(gap#1 實作);review `20260705_topology_tractability_review_and_gaps_01.md`。
