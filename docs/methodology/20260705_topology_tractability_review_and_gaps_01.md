<!--
建立時間: 2026-07-05
類型: review 互補紀錄 — 對 68ced6f(Track A group-Steiner gap#1/2/3)的獨立驗證 + 運算量/剪枝可解性論證 + 查證後引用修正 + 真正剩餘缺口
狀態: review 紀錄(不重述形式化;形式化 SoT = 20260704_formal_problem_statement_topology_from_cooccurrence_01.md)
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/candidate_trees.json, docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md
provenance: 數字全引用 candidate_trees.json(commit 68ced6f 後)+ 68ced6f commit message;本檔為 review 互補,不產生新分析數字。運算量論證為源碼結構推論(標 L1/L4)。
-->

# 拓樸重建 — 68ced6f 獨立驗證 + 運算量可解性論證 + 引用修正 + 剩餘缺口

> **本檔定位**:形式化主體在 `20260704_formal_problem_statement_topology_from_cooccurrence_01.md`(並行 session,v3 Model A)。本檔**不重述**,只補三件它沒有的:① 對其實作(68ced6f)的獨立驗證;② 「read-span/grouping ⇒ 運算量與樹數受資料界定」的可解性論證;③ workflow 查證後的引用修正 + 真正剩餘缺口。

## §1 68ced6f 已完成(誠實盤點 + 獨立驗證)

並行 session commit `68ced6f`「Track A group-Steiner 實作補完 gap#1/2/3 + 甲基排序 Bernoulli 否決」已落地並經 2 workflow 對抗驗證。本輪獨立複核:

| 項 | 內容 | 獨立驗證 |
|---|---|---|
| **gap#3 m-通道(點 1 CN 分流)** | 70 recurrence 接整數 CN(SEQC2 gain_cn)→ **11 artifact(m>1)棄 + 9 candidate(m=1)留 + 50 copy-neutral LOH 未解(VAF L3)** | identifiability_table `recurrence=70` ✓ |
| **gap#2 全枚舉 + identifiability** | 舊碼 silently drop 32 parent-choice 歧義 + 假稱「0 歧義」→ +parent_choice_enum + dedup;**69 A_ambiguous → 44 determined + 25 真歧義** | B1 `69=25+44` ✓;`conflict=18` |
| **gap#1 partial subcube 救回(點: partial-read as group)** | 「無 full-cov 但 ≥2 位點 partial read」建 `E_subcube_recovered`(不僭稱 A_determined)。**救回 2403 區、分母 3885→6288、A_determined 1764 零 regression** | 桶加總 **3885** ✓;+2403 = **6288** ✓ |
| **甲基排序 Bernoulli 否決** | read-level Bernoulli 取代 mean-Δβ,實測 distal ρ=0.068<0.144 更差、全 n.s. → 「換 Bernoulli 救 ordering」**已測否決** | 對齊本研究「甲基不排序」紅線 |

🟢 **數字自洽性(本輪 grep-verify)**:`determined 1808 + ambiguous 25 + recurrence 70 + conflict 18 + B_pairwise 943 + C_underdetermined 544 + other 477 = 3885`(with-vector);`+ subcube_recovered 2403 = 6288`(總分母)。全部精確加總 ✓。〔L1〕

## §2 運算量/樹數可解性論證(形式化文件缺此,補上)

> 回應「read 有多長、過多少 sSNV 就處理多長,運算量與可能樹結構是否其實影響不大」——**是,且這正是逐區精確可解的底層機制。**〔源碼✓ L1 + 推論 L4〕

**運算量不由 n(或 2ⁿ)決定,由資料的實際連結結構決定**,兩個都小:
1. **|R| = 相異觀測基因型數**(非 2ⁿ)。`populations` 每區常 2–4 個。MSA-DH 精確法 **Õ(3^|R|)** 對 |R| 指數、對 n 不指數 → |R| 小即便宜。〔文獻✓〕
2. **共讀 grouping 剪枝後的殘餘自由排列數**。`resolve_and_order`:co_linked→union-find 併 block(消排列)、nested→方向定死(消排列)、**只有「無共讀」對留自由排列**;候選樹 = groups 偏序的線性延伸,再 hard-cap `CAP=24`。determined 區=1。

🔑 **read 越長 → 對可解性越有利**:一條 read 跨 k 個 sSNV = 直接觀測那 k 位點聯合基因型 → 塌成一觀測頂點/co_linked block → **一次消掉 k! 排列**。組合爆炸只在「從未被任何 read 共觀測」的對上殘存 → 而那裡誠實輸出「等機率 / 非可枚舉」(B_pairwise 943 + C_underdetermined 544),**不是硬算**。

**結論**:運算量與樹數受資料連結結構自然界定;會爆的地方本來就標欠定而非硬算。兩道硬界(每區 densest-8 → n≤8;CAP=24)保證常數上界。實證:全基因組 O2 平行 ~6 分鐘。
🟢 **CAP 實測佐證(本輪 candidate_trees.json)**:69 個 B1 枚舉區中 **capped 區 = 0、最大 true_count = 2**(所有區候選樹 ≤2 棵)→ CAP=24 遠不 binding,共讀 grouping 把樹空間壓到連 cap 邊都碰不到。〔L1〕 → §2 論證經驗成立:無合法候選被截。
⚠ **綁定條件**:grouping 剪枝正確性依賴共讀計數可信 → 這就是「≥2-read 承重、單讀輔助、衝突非對稱門檻」與本論證綁在一起的原因(見 §4 缺口)。

## §3 查證後引用修正(refine 形式化文件 §13,workflow wf_272a702f 已驗)

- 🔴 **directed rooted 大簡約 NP-complete 出處 = Day–Johnson–Sankoff 1986**(Math. Biosci. 81:33-42),**非** Foulds–Graham 1982(後者是無根 Wagner)。〔文獻✓〕
- **MSA-DH**(arXiv:2110.02830)= 有向超立方體最小斯坦納樹形圖的**精確抽象問題** + 複雜度骨架,精確法 **Õ(3^|R|)**(對相異終端數,非 O(1))。〔文獻✓〕
- 🔴 **Scelestial(PLOS CB 2022)不是 group-Steiner-on-subcube**:它是 neighbor-joining 推廣,解**幾何(歐氏)斯坦納樹近似** + lineage imputation 補缺失。**不能引用為「最近形式匹配」**——反而**強化白地**(無工具做 group-Steiner-on-subcube)。〔文獻✓〕
- **partial read = 子面(subface)= IDPP**(Pe'er et al. 2004)的標準表示,**非新貢獻**;novelty 收窄為「單分子 long-read BULK + germline-有向根 + 支持度門檻 + 全等成本枚舉」四者組合。〔文獻✓〕
- 複雜度正名:GST 有**多項式**時間 polylog 近似(GKR 2000);有向 Steiner 2024 亦有(arXiv:2412.10744)。但 n≤8 全無關,暴力精確更強。

## §4 真正剩餘缺口(68ced6f 之外,誠實列)

| 缺口 | 現況 | 規模 | 是否該做 |
|---|---|---|---|
| **VAF/CCF pigeonhole PRUNE** | 只有 `vaf_ok` 交叉核對旗標(enumerate line 163-168),**未拿來 prune 候選集** | 中 | ⚠ 與形式化 §8「通道 A/M 只精修不重定義」設計相關,須先確認不與其 non-circular 紀律衝突 |
| **single-read 輔助旗標**(點 2) | grep 無;現只 ≥2-read 承重 | 小 | 可加(純 annotation,不驅動結構) |
| ~~**CAP=24 是否截掉合法候選**~~ ✅ **已驗** | 實測 **capped=0、max true_count=2** → 未截任何合法候選 | — | ✅ 完成(本輪) |
| **HTML 每位點狀況更細** | 已有 per-region 明細(genotype 向量 + 樹 SVG + regSit 上色 + rd_perregion 交叉);逐 **sSNV** 的 2×2/VAF/CN 是否入表待確認 | 中 | 你最後那句的重點,待確認現況後決定增補幅度 |

## §5 建議(給使用者)

核心已完成 + 驗證,**不重做**。剩餘為精修 + HTML 細化。因**目前 5 個並行 session 活躍於共用檔**,任何對 `enumerate/topology/build_topology_workstation` 的編輯建議**開 git worktree** 隔離後做,避免 clobber。優先序:① CAP 截斷量化(零風險驗證)② HTML 逐 sSNV 狀況入表(你的明確需求)③ 視需要 VAF pigeonhole prune / single-read aux(須先對齊形式化 §8)。

**關聯**:形式化 SoT `20260704_formal_problem_statement_topology_from_cooccurrence_01.md`;架構裁決 `20260703_clone_read_probabilistic_labeling_architecture_01.md`;commit `68ced6f`;memory `project_clone_read_labeling_probabilistic_architecture` / `project_incompatible_reclassification_hidden_node_pairwise_fourgamete`。
