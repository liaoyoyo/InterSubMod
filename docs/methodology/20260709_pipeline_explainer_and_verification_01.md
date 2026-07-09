<!--
建立時間: 2026-07-09
類型: 端到端 explainer + 可獨立驗證 — sSNV 共現 → 布爾超立方體群組斯坦納樹 → 輔助驗證 全流程 + 每數字 provenance + 一鍵驗證腳本
狀態: explainer(理解)+ verification(驗證);形式化 SoT = 20260704_formal_problem_statement_topology_from_cooccurrence_01.md
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/sm_completeness_ledger.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/sm_region_integration.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/candidate_trees.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/partial_read_recoverable_fraction.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/four_gamete_aa_methyl_demo.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json
provenance: 所有數字由 scripts/verify_pipeline_numbers.py 於 2026-07-09 對回原始檔重算 + 一致性檢查(全 ✓ 通過)。跑該腳本即可獨立驗證。
-->

# 端到端 explainer + 可獨立驗證 — sSNV 共現 → 群組斯坦納樹 → 輔助驗證(HCC1395)

> **兩個目的**:① **理解**(下方 7 站流程)② **驗證**(§8 一鍵重算,每數字對回原始檔)。

## 核心一句話
把「癌症演化樹重建」變成「在布爾超立方體 `{0,1}ⁿ` 上、以 germline(`0ⁿ`)為根、每邊只翻一個 0→1、每條 read(=一個子面/group)至少相容一節點的**最小有向斯坦納樹形圖枚舉**」。**輸出不是一棵樹,是「所有最小相容樹 + 可辨識性分類」。**〔形式化 SoT:`20260704_formal_problem_statement...`〕

---

## §0 問題怎麼擺（數學 ↔ 生物）
| 數學物件 | 生物意義 |
|---|---|
| 超立方體頂點 | 一個基因型(哪些 sSNV 已獲得)|
| 根 `0ⁿ` | germline |
| 有向邊(翻一位) | 一次 infinite-sites 獲得突變 |
| Steiner 內部節點 | 推斷的隱藏祖先 clone(標 `inferred`)|
| 一條 partial read | 超立方體的一個**子面/group**(未覆蓋位點=自由座標)|
| 最小成本樹 | parsimony |
逐區 **n≤8** → 超立方體 ≤256 頂點 → **常數大小 → 可精確窮舉**。

## §1–§5 主流程（每站一個數字）

| 站 | 處理 | 量化(HCC1395) |
|---|---|---|
| **1 輸入→sSNV** | BAM(腫瘤+normal)+VCF → somatic sSNV | **35,332 sSNV**(TP 30,490 / FP 4,842)|
| **2 sm_linkage 單分子共現** | 逐 read 讀各位點 REF/ALT(**天生 0/1**)+ 每對 2×2 共現 | **linked 21,554** / underpowered 5,458 / isolated 8,320 |
| **3 region 整合** | ≤50kb 共現串成區;能被單 read 蓋滿的組合 = `populations` | **7,143 區**(n≥2)= 有向量 4,256 + **空向量 2,887(40.4%)** |
| **4 topology 可辨識性** | 2×2 哪格為零定關係:co_linked(併block)/ nested(定向)/ 四配子(送 recurrence-衝突) | 見下表(3,885 有向量進樹)|
| **5 enumerate co-optimal** | 共讀 grouping 剪枝 → 枚舉全等成本樹 + softmax + equiprobable | **B1 69 → 44 determined + 25 真歧義**;CAP **0 截斷**、最大 **2** 棵 |
| **5b partial-read 救回** | 空向量區用「蓋≥2 位點」的 partial read 建 E_subcube_recovered | 救回 **2,403 區** → 分母 3,885→**6,288** |

**站 4 可辨識性(3,885 有向量區,桶加總=3,885 ✓)**:

| 類別 | 數 | 白話 |
|---|---:|---|
| **determined** | **1,808** | 唯一樹(= 單分子直接 A_determined **1,764** + 枚舉後收斂 **44**)|
| ambiguous(真歧義) | 25 | 多棵等成本樹 |
| recurrence-required | 70 | 位點翻≥2次(送 m-通道)|
| conflict(真 artifact) | 18 | 允許 recurrence 仍無樹(成環)|
| **非可枚舉 B_pairwise** | 943 | 只 pairwise 拼、非單分子 |
| **非可枚舉 C_underdetermined** | 544 | 單 ALT 群、缺連接 read |
| other | 477 | — |

**站 5b partial-read**:空向量 2,887 中 **2,316(80.2%)有 pairwise 可救** / **1,692(58.6%)連≥3位點可拼樹** / **571(19.8%)真單位點救不回**(80% 是 n=2 稀疏)。

---

## §6 輔助驗證 — 3 個有界、非循環通道 🔑

**都只精修/篩選/拆分,絕不重定義樹、絕不排序。**

1. **CN 多重度通道**（拆 recurrence 的 artifact vs 真;接 SEQC2 整數拷貝數）
   70 recurrence → **11 artifact(m>1 棄)** + **9 candidate(m=1 留)** + **50 copy-neutral LOH(m 未定,VAF L3 軟旗標)**。（11+9+50=70 ✓）最大 win = 11 是粗 categorical gate 漏掉的。
2. **VAF 交叉核對**：較早突變應較高 VAF → `vaf_ok` 旗標核對 nested 定向（目前核對,尚未 prune 候選集）。
3. **甲基通道（實測否決,不進目標函數）**：
   - 四配子拓撲方向:20 對中唯二 distal 顯著**全在 cross-HP CN-gain**(抓 confound 非 lineage),6 個乾淨案例 distal **全不顯著**。
   - Bernoulli 換 mean-Δβ 救 ordering:distal ρ=0.068 更差、全 n.s. → **否決**。
   - → 甲基 = 突變 cis 足跡,**bounded-auxiliary**:只負篩 + L3 註記 + 上游 germline phasing;**絕不進 likelihood/排序/破 tie**。

---

## §7 輸出與誠實邊界

| | 意思 | 現況 |
|---|---|---|
| **SOLVE-1 枚舉精確解集** | 每區精確窮舉全等成本樹 + 分類 | ✅ 可以(n≤8 常數規模、暴力精確、CAP 0 截斷)|
| **SOLVE-2 指派唯一全域樹** | 每區/全域一棵確定樹 | ❌ 理論不可(資訊極限非計算極限)|

- **非可辨識 = 943 + 544 = 1,487 區**:單 bulk 資訊不足 →「定不出來」就是精確答案,需 multi-region/single-cell。
- **封頂 ⭐3**:單 bulk 只能 characterize、不能 confirm subclone;Steiner 節點標 `inferred`。

---

## §8 一鍵驗證(讓你獨立確認每個數字)

```bash
cd InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts
python3 verify_pipeline_numbers.py
```
該腳本**只讀不寫**,把上方每個 headline 數字對回原始 JSON 重算 + 一致性檢查(三桶加總、桶=3885、+subcube=6288、B1=44+25、CAP=0、m-通道=70…),**並含 §10 的 Q1 n_sSNV 分佈 / Q2 pairwise 假關聯 / Q3 跨樣本 determinacy**(後兩者 optional,需 sm_linkage pairs / MSROOT)。**全 ✓ = 本文數字可信;任何 ✗ = 該數字與原始檔漂移需查。** 2026-07-09 實跑全 ✓。

### Provenance 表（數字 → 來源檔:key）
| 數字 | 值 | 來源 |
|---|---|---|
| sSNV total / TP / FP | 35,332 / 30,490 / 4,842 | `sm_completeness_ledger.json:universe_*` |
| 三桶 linked/under/isolated | 21,554 / 5,458 / 8,320 | `sm_completeness_ledger.json:buckets` |
| 區數 / 有向量 / 空向量 | 7,143 / 4,256 / 2,887 | `sm_region_integration.json:regions`(populations 空與否) |
| determinacy 7 桶 | 1808/25/70/18/943/544/477 | `candidate_trees.json:summary.identifiability_table` |
| subcube_recovered / 分母 | 2,403 / 6,288 | 同上 + 3885 |
| B1 / determined / ambiguous / CAP | 69 / 44 / 25 / 0截斷 max2 | `candidate_trees.json:summary` + `candidate_trees[]` |
| partial-read 2316/1692/571 | 80.2%/58.6%/19.8% | `partial_read_recoverable_fraction.json:summary` |
| CN m-通道 11/9/50 | =70 | `topology_per_region.json:detail[].determinacy`(recurrence_*) |
| 甲基四配子 20/12/6/2 | — | `four_gamete_aa_methyl_demo.json:summary` |

---

## §9 一圖串起來
```
35,332 sSNV
  → sm_linkage: linked 21,554 / underpowered 5,458 / isolated 8,320
  → 7,143 區 = 有向量 4,256 + 空向量 2,887(40.4%)
  → topology(3,885): determined 1,808 / ambiguous 25 / recurrence 70
                     / conflict 18 / 非可枚舉 1,487 / other 477
  → enumerate: 69 → 44 determined + 25 真歧義(CAP 0·max 2 棵)
  → partial-read 救回 2,403 → 分母 6,288
  → 輔助驗證: CN m-通道(11棄/9留/50未定) · VAF 核對 · 甲基 Bernoulli 否決
  → SOLVE-1 精確枚舉 ✓ / SOLVE-2 唯一樹 ✗(非可辨識 1,487) / ⭐3
```

---

## §10 補充:方法細節合理性(數據支持 + 反思)

### §10.1 n≤8 cap 是否合理〔L1〕
n_sSNV 分佈極度集中低端:**n≤8 覆蓋 97.51% 區**(僅 178=2.49% 超過,max=150 長尾極稀),截斷只丟 **6.73% sSNV**;且 commit `151a86d` 已對 >8 區加**全 pairwise 建樹**(去 MAX_SNV=8)→ 大區走 pairwise 不遺失。累積:n≤2 51.9% / ≤4 85.1% / ≤5 91.0% / ≤8 97.5%。→ **cap 是資料驅動的合理工程界**(精確窮舉路徑上限),非任意;僅 178 個 >8 區的**精確全向量樹**降級(有 pairwise 替代)。

### §10.2 partial-read 救回反思(5 邊界)〔L1 機制 + L4 反思〕
✅ 非循環(用真共現)/ ✅ 不 overclaim(標獨立 `E_subcube_recovered` 弱層,不當 A_determined)/ ⚠ 「救回 2,403 區」= 新增部分樹結構(**上界**,非完整還原;實誠 tier = chain≥3 的 58.6%)/ ⚠ 帶錯 partial read 可造假衝突(緩解:≥MINREAD + 對抗驗證 wf_526f6c8f,非零風險)/ ⚠ 已知 nic/population 不一致 edge(源碼自標,rr=0 時)/ ⚠ 僅 HCC1395。

### §10.3 「靠共享 B 串」是否合理完備 + 能否用超立方體〔L2〕
- **範圍澄清**:B_pairwise 943 中 **n=2(單對·無串接)614=65.1%**、**n≥3(真靠共享 B 串)329=34.9%**。65% 根本不涉傳遞。
- **✅ 合理**:Gusfield 定理(所有對兩兩四型相容 ⟺ perfect phylogeny 存在);實證 **nic=0**(觀測對全相容)+ **98.7% 全對直接觀測** → 單拷貝下串接數學上保證正確、非假關聯。
- **🔴 完備? 不完全,兩缺口**:
  - **Gap A 傳遞假設(小)**:12 區(1.3%)含未觀測對 → 超立方體枚舉對未觀測對保持自由 → 轉多候選,已緩解。
  - **🔴 Gap B 共享 B 可能多拷貝(大,主缺口)**:B_pairwise **51.6%(487)是 CN-gain、僅 2.9%(27)CN-neutral**。CN-gain 下共享位點 B 可能 2+ 拷貝 →「A-B on read1」+「B-C on read2」可能跨不同拷貝 → 誤併譜系。最該盯 = **n≥3 且 CN-gain 的 187 區(19.8%)**。漏網機制:A-C 未觀測時 / 或「pairwise-相容但實為多拷貝結構」(pairwise-相容 ≠ 單拷貝譜系);**m-通道目前只拆 four-gamete recurrence,未覆蓋 pairwise-相容 CN-gain**。
- **⚠ 更正(§13.7)**:先前「pairwise 假關聯風險只 1.3%」**只算 Gap A、漏 Gap B**;完整風險面 = 1.3% 傳遞 + CN-gain 主導的多拷貝串接不確定(51.6% gain / 最壞 187 區)。
- **能否用超立方體?** ✅ 可以(pairwise co-read = 2-位點子面 = group;gap#1 已做),對 Gap A 更安全(枚舉全相容樹非硬串);但**對 Gap B 無自動偵測** → 需 CN 資訊介入(擴 m-通道到 pairwise-相容 CN-gain 區)。

### §10.4 跨樣本拓撲分佈 + 誠實邊界〔L1 分類 / L3 比較〕
| 樣本 | n | determined% | conflict | B_pairwise | C_underdet |
|---|--:|--:|--:|--:|--:|
| HCC1937 | 1,886 | **52.7%** | 0 | 285 | 132 |
| HCC1954 | 2,837 | 40.2% | 1 | 343 | 222 |
| HCC1395_DORADO | 3,418 | 34.5% | 4 | 412 | 307 |
| H1437 | 8,169 | 25.3% | 66 | 1,027 | 382 |
| H2009 | 9,444 | 18.6% | 265 | 827 | 246 |
| COLO829 | 6,613 | 18.5% | 1 | 660 | 607 |

✅ 拓撲可乾淨分 6-7 類 + 算比例(每樣本一張)。🔴 **但 determined% 差異主由資料品質(coread 深度/read 長度)驅動,非生物**(COLO829 18.5% = 低 coread artifact,非 clone 少);pseudoreplication(真 n=7,p 值無意義)、CN 混淆、不同 truth set → **跨樣本比較 = L3 描述性,非驗證整體 clone 的生物結論**;⭐3 封頂(缺正交真值)。**可合法用**:樣本內「結構可解析度」誠實陳述;不可用:「A 比 B clone 複雜」。

---

**關聯**:形式化 SoT `20260704_formal_problem_statement_topology_from_cooccurrence_01.md`;review `20260705_topology_tractability_review_and_gaps_01.md`;架構 `20260703_clone_read_probabilistic_labeling_architecture_01.md`;驗證腳本 `_assets/20260627_subclone_4axis_teaching/scripts/verify_pipeline_numbers.py`。
