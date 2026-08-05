<!--
建立時間: 2026-07-24
目標: 實作前稽核 exact-PS C++ topology/read-AF HCC1395 gate 與七資料集推廣
處理範圍: Task Type B / comprehensive validation / HCC1395 gate then all seven datasets
關聯檔案: InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/implementation-notes.md
-->

# Exact-PS C++ topology + read-AF 實作前稽核

## 1. Verdict

**PROBE → conditional GO**。先允許 HCC1395 C++ 實作與逐 unit oracle 驗證；
七資料集執行維持鎖定，直到 HCC1395 同時通過 family completeness、Python/C++
語意一致、AF exact-rational parity、shape parity、fail-closed 與加速門檻。

服務目標：G1、G3、G4、G5。

## 2. Task type 與範圍

- Task Type：**B — Comprehensive validation**。
- Stage 1：HCC1395 chr1–22，使用 7/22 exact-PS、k<=12 authority。
- Stage 2：相同 frozen contract 的七資料集 chr1–22；不可用 subset 冒充完成。
- 不重新複製 BAM／VCF；只讀 authority，輸出壓縮 unit result、aggregate 與 receipt。

## 3. 已知事實

1. 7/22 HCC1395 分區已通過 22/22 chromosome、5,342 checks，
   `cross_PS=0`、`cross_HP=0`、Python/C++ partition mismatch=0。
2. 7/23 Python topology adapter 的 9,600 primary mutation-bearing route-blocks 中，
   9,180 complete、420 因舊 `extra_cap/per_level_budget` incomplete。
3. 舊 read-AF score 對每條取得 mutation `j` 的邊 `p→c` 加總
   `sum_{i in p}(AF_i-AF_j)`，並以 Python `Fraction` 比較。
4. 目前 production all-seven authority只證明 exact-PS read linkage；尚未證明
   topology 或 AF ranking 已完成。

## 4. 關鍵假設、反證與防線

| 假設 | 主要風險 | 反證／防線 |
|---|---|---|
| 固定 vertex set 時 AF objective 可依 child 分解 | 若 score 含全樹 normalization、CN/LOH 或跨邊項，parent 選擇會耦合 | 逐行核對舊 Python；以 exhaustive tiny-tree oracle 比較最佳分數、並列數與 winner digest |
| C++ 可完整找出全部最小 vertex sets | family 數長尾仍可能爆炸；只找到 objective 不等於 family complete | objective certificate 與 family certificate 分離；cap/deadline 一律 incomplete/ABSTAIN |
| recurrence-allowed parent choices彼此獨立 | strict infinite-sites 會加入跨邊 character constraint | 凍結 recurrence-allowed contract；strict perfect phylogeny 不在本輪暗中加入 |
| 七資料集具有同一 exact-PS input schema | 樣本路徑、sidecar、call alphabet 或 chromosome completeness 不一致 | Stage 0 manifest 稽核；缺任一 chr/authority/schema 即不啟動該 cohort claim |
| AF 分母可由 exact-PS block evidence重現 | 舊 `col_coverage_by_hp` 與新 molecule evidence口徑可能不同 | HCC 建立雙路 oracle；先做 exact legacy parity，再把新口徑列為獨立 sensitivity |

## 5. AF factorization 要證明的命題

對固定可行 vertex set `V`，每個非根 child `c` 的候選 parent 集合為
`P(c)={p in V: p<c, d_H(p,c)=1}`。若總分只有 edge-additive：

`S(T)=sum_{c != root} s(parent(c),c)`，

且每個 child 選 parent 不會改變其他 child 的可行集合，則：

- `S*(V)=sum_c max_{p in P(c)} s(p,c)`；
- 最佳 parent-tree 並列數為
  `product_c |argmax_{p in P(c)} s(p,c)|`；
- 當且僅當每個 child 都只有一個最佳 parent 時，固定 `V` 的最佳樹唯一。

Directed Boolean hypercube 依 Hamming weight 嚴格遞增，因此任選一個
Hamming-1 predecessor 都不會成環；每個非根一個 parent 且 root 可達，即形成
`|E|=|V|-1` 的 rooted arborescence。

此命題**不**消除 optimal vertex-set family 的列舉成本；跨 vertex sets 仍須比較
`S*(V)`，並對全域最高分的 sets 加總最佳樹並列數。

## 6. 預先固定成功門檻

### HCC1395 gate

1. 所有舊 Python complete units 的 mutation universe、minimum hidden count、
   optimal vertex-set digest、tree count、shape census一致；mismatch=0。
2. AF 有效 units 的 exact rational site AF、最佳 score、最佳樹並列數、
   unique/tied status 與 winner edge digest一致；mismatch=0。
3. 420 個歷史 capped units不得被誤報 winner；若新 C++ 完整解出，須另列
   `newly_completed` 並以獨立 exhaustive/synthetic certificate 驗證。
4. corrupted input、缺 AF、zero denominator、cap、deadline、overflow fixture
   全部 fail closed。
5. 在同一凍結輸入、單 worker、warm/cold 分開的 controlled benchmark，
   C++ 核心 wall time 至少比 Python oracle 快 2×；否則不宣稱加速。

### 七資料集 gate

1. 七資料集 × chr1–22 authority/input checks全通過。
2. 每個 requested chromosome有 terminal receipt；任何缺失使 cohort completion
   claim 失敗。
3. aggregate 分母只納入 `complete=true`，另列 incomplete reason census。
4. 輸出空間預估加 20% 安全係數後仍低於可用空間；否則先停，不寫部分 cohort。

## 7. 主要失敗模式

1. 固定 parent factorization 成立，但 optimal vertex-set family 本身太大。
2. 舊 Python AF 依賴 legacy HP 聚合分母，不能直接移植到 exact-PS block。
3. LongLineage C++ solver 與舊 recurrence-allowed Python 的 group-coverage語意不同。
4. k<=12 只限制 `2^k` vertex universe，不能給定固定 wall-time 上限。
5. allele-specific CNA、LOH、purity 與 mapping/phasing error 使 molecule AF
   不等於 cell prevalence；即使計算正確也不能升格成真實 clone 結論。

## 8. Step → Verify

1. 凍結 source/input manifests與 Python oracle。  
   → 驗證：絕對路徑、size、SHA-256、schema、sample/chromosome coverage有 receipt。
2. 建立 tiny exhaustive oracle與 factorization property tests。  
   → 驗證：所有 k<=5 exhaustively generated fixed-vertex cases分數/並列數一致。
3. 實作 C++ exact solver、AF factorization、shape與 compressed output。  
   → 驗證：build/test exit 0；arborescence、group coverage、minimum objective、
   family completeness、score守恆 validator全 PASS。
4. HCC1395 chr1–22 paired run。  
   → 驗證：逐 unit parity=0 mismatch；controlled timing輸出 p50/p95/max/RSS。
5. HCC gate通過後執行七資料集。  
   → 驗證：all requested chromosomes terminal；per-sample與aggregate數據及時間可重算。
6. 獨立 red-team稽核。  
   → 驗證：以 receipt重新計算關鍵 totals/digests，且未讀取 producer aggregate。

## 9. 決策鎖

- HCC implementation：**unlocked**。
- 七資料集 full run：**locked until HCC gate PASS**。
- biological clone/subclone claim：**locked；本輪只允許 graph geometry + AF ranking**。

