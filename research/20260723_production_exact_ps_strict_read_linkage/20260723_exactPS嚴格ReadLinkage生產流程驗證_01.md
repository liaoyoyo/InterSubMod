<!--
建立時間: 2026-07-23 11:30 +08:00
目標: 紀錄 chromosome × exact PS × HP strict endpoint graph 的 production 實作、全7資料集結果與完成層級
處理範圍: 7 technical datasets × chr1-22 L1正式完成；production strict directed topology尚未執行
關聯檔案:
  - InterSubMod/scripts/build_strict_ps_hp_regions.py
  - InterSubMod/scripts/summarize_strict_ps_hp_regions.py
  - InterSubMod/tools/strict_endpoint_graph.py
  - InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/implementation-notes.md
-->

# Exact PS × HP strict read-linkage production 驗證

> 狀態：`validated_L1_only`。7 technical datasets × chr1–22 的 strict read-linkage regions 已完成；current production strict directed topology 為 0/7，clone／parent–child／fusion 為 0/7 驗證。

## 1. 先講結論

舊 production grouping 不符合研究問題：它先用 50 kb 幾何距離分區，再跨 exact PS 聚合同一數字 HP。新版把 primary region 改為：

1. 容器：`dataset × chromosome × exact nonmissing PS × primary HP family`。
2. linkage graph 節點：容器中至少被一條 canonical molecule 以 fixed `R/A` 觀察到的 sSNV 位點。
3. linkage graph 邊：同一 canonical molecule 在兩個 endpoint 都是 fixed `R/A`；distinct molecule support ≥3。
4. tree-eligible region：上述邊形成的 `k>1` connected component。
5. `k=1`：只記為 `ABSTAIN_SINGLETON_UNLINKED`，表示該 membership 在此 exact container、threshold=3 下沒有合格 edge；它仍保留在 graph component audit registry，但不進下游 k≤12 computational partition、Steiner tree、topology 或 clone/subclone 分母。

HCC1395 chr1–22 全量結果：`S=79,687`，其中形成 62,651 個 exact-container node memberships；threshold=3 有 39,846 個 graph components，但 28,384 是 singleton abstention，真正可進建樹流程的 read-linked candidate regions 是 **11,462**。

## 2. 兩種「節點」不可混用

| 層 | 圖 | 節點 | 邊 | 用途 |
|---|---|---|---|---|
| L-linkage | sSNV read-linkage graph | 單一 sSNV 位點 | 同 molecule 對兩 endpoint 的 fixed R/A 支持 | 決定哪些位點有資格共同分析 |
| L-tree | Boolean hypercube / mutation-state graph | k 位點的 R/A 狀態向量，例如 `000、100、110` | Hamming distance=1 的一次突變 | 枚舉最小相容 mutation-state trees |

因此，「同一 PS 內的所有 sSNV 直接作一棵 Steiner tree」並不充分。PS 僅提供 LongPhase-S 演算法所定義的局部 HP orientation 座標系；如果兩組 sSNV 間沒有 read chain，資料無法辨識其共現順序，把它們硬放入同一個 k 維超立方體只會創造沒有觀測根據的組合。

若同一 exact PS×primary HP family 內所有 sSNV 經 A–B、B–C 等 threshold-qualified endpoint edges 連成同一 component，則可共同描述一個局部候選演化結構。這仍只支持候選 mutation-state tree，不等於唯一細胞 clone tree。

### 2.1 Canonical region 定義與記號

- `P`：unique `(dataset, chromosome, exact nonmissing PS)`；相同 PS 數字只在同一 dataset/染色體內有意義。
- `H`：`(dataset, chromosome, exact nonmissing PS, primary HP family)`，是建 linkage graph 的容器。raw `1`、`1-1`、`1-2` 合併為 HP1，raw `2`、`2-1`、`2-2` 合併為 HP2；HP3、HP4、unphased 與 missing PS 不進 primary graph。
- `W_i`：容器內由 support≥3 的 direct endpoint edges 形成的 maximal connected component，且 `k_i=|W_i|≥2`；這是本流程唯一稱為 read-linked candidate region 的單位。
- `start_i=min(position in W_i)`、`end_i=max(position in W_i)`、`span_i=end_i−start_i`；這只是顯示用 coordinate envelope，`W_i` 實際是位點集合，不是從 start 到 end 間全部位點的連續區間。
- `B_ij`：k>12 的 `W_i` 經 k≤12 限制切出的 computational block；`B` 不是新 region，必須用 source component ID 回聚到 `W_i`。
- `T`、`Topo`、`C` 都屬下一階段的候選樹、拓撲與 mutation-state group 數，不參與 `W_i` 的區域定義。

PS 是 LongPhase-S 的演算法 phase-block 標籤，只提供局部 HP orientation 座標系；它不是 phasing 必然正確、真實 clone boundary 或生物演化邊界的保證。不同 PS 間 HP1/HP2 的相對 orientation 在本層未定義。

### 2.2 Canonical molecule 與 fixed call

- canonical molecule ID 定義為 `SHA256(dataset, RG-or-".", QNAME)`。
- 只納入 mapped primary alignment，排除 unmapped、secondary、supplementary、QC-fail 與 duplicate，且要求 MAPQ≥20。
- 相同 alignment identity 只有在 calls、qualities、MAPQ、HP、PS 全部一致時才能 collapse；衝突則 fail closed。
- fixed `R/A` 表示 endpoint 的 BQ≥20，且鹼基等於 REF 或 ALT。
- `X/L/O/D/S` 分別表示未覆蓋、低品質、其他鹼基、deletion、ref-skip；它們都不替該 endpoint 建 edge。
- 每個 canonical molecule 對同一 endpoint pair 最多投一票；按 genomic site order 分為 `RR/RA/AR/AA`，edge total support 必須等於四種 state support 總和，且 total≥3 才保留。

## 3. Partial read 的規則

- 對排序位點 `s1<s2<s3`，read pattern `R-X-A` 只對 `(s1,s3)` endpoint pair 投一票 `RA`。
- `X/L/O/D/S` 不會建立 `(s1,s2)` 或 `(s2,s3)`，也不會因為 read 的座標 span 而被動補邊。
- `(s1,s2)≥3`、`(s2,s3)≥3` 時，即使 `(s1,s3)=0`，三點仍可由 graph transitivity 合成同一 component；兩條 edge 可以由不同 molecule 集合支持。
- 因此 W 不表示有一條 read、同一細胞或同一 clone 橫跨全部 k 個位點；`RA/AR` 依 genomic site order 編碼，也不是祖先—後代方向。
- 距離不參與連邊；若三條 molecule 都 fixed-observe 相距 >50 kb 的兩 endpoint，該邊仍合法。

## 4. HCC1395 chr1–22 正式漏斗

來源：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_strict_regions_v2/summary/summary.json`，`all_pass=true`。

| 層 | 單位 | 數量 | 比例／說明 |
|---|---|---:|---|
| S | LongPhase-S PASS autosomal sSNV loci | 79,687 | 100% |
| Active loci | 至少進入一個 exact PS×HP 容器的 physical loci | 36,384 | 45.6586% of S |
| Containers | 有 active node 的 exact PS×HP | 7,118 | 3,971 unique chromosome×PS |
| Node memberships | `(locus, exact PS, HP)` memberships | 62,651 | 同一 physical locus 可重複屬於不同容器 |
| Qualifying edges | direct endpoint support ≥3 | 76,202 | 每列含 RR/RA/AR/AA read counts |
| All graph components | 含 singleton 的守恆 registry | 39,846 | 不可直接當 tree region 分母 |
| Singleton abstain | k=1、沒有合格連接 | 28,384 | 71.2343% of all components；排除建樹 |
| **Tree-eligible regions** | **k>1 read-linked candidate components** | **11,462** | **28.7657% of all components；tree-eligible analysis 主分母，不是 clone prevalence 分母** |
| Large regions | k>12 | 90 | 0.7852% of tree-eligible；max k=153 |

單點排除的集合守恆（threshold=3）：

- exact-container memberships：`62,651 = 28,384 singleton abstain + 34,267 tree-eligible`；建樹使用 54.6951%，排除 45.3049%。
- 去重 physical loci：36,384 個 active loci 中，22,689 個至少在一個 PS×HP 容器進入 read-linked region；13,695 個只出現為 singleton。
- 另有 5,717 個 physical loci 在某一 PS×HP membership 是 singleton，但在另一 HP/PS membership 又有合格連接。因此排除必須在 membership 層執行，不能把該 genomic locus 全域刪除。

HP 分層：

| HP | all components | singleton abstain | tree-eligible regions | node memberships |
|---|---:|---:|---:|---:|
| HP1 | 20,059 | 14,261 | 5,798 | 31,696 |
| HP2 | 19,787 | 14,123 | 5,664 | 30,955 |

tree-eligible region 的 k 分布前段：k=2 6,840（59.6754%）；k=3 2,578（22.4917%）；k=4 1,041（9.0822%）；k=5 455（3.9696%）。完整表在 `.../summary/k_distribution.tsv`。

Region 尺度（W=11,462）：平均 k=2.9896，median=2，p90=4，p95=5，p99=10，max=153。Coordinate span 的 min/median/p90/p95/p99/max 為 1/17,837/48,364/63,537/103,429/272,774 bp。

三種 50 kb 指標必須分開，因為事件與分母不同：

| 指標 | 精確定義 | 數量 | 分母 |
|---|---|---:|---:|
| direct edge >50 kb | qualifying endpoint pair 的兩端座標距離 >50 kb | 47 | 76,202 retained direct edges |
| W span >50 kb | `max(position in W)−min(position in W)>50 kb` | 1,064 | 11,462 W |
| W adjacent gap >50 kb | 將 W 成員依座標排序後，至少一個相鄰 membership gap >50 kb | 4 | 11,462 W |

47 條長 direct edges 中，support=3 有 35 條、support=4 有 12 條；24 條為 RR-only callability linkage，23 條含 RA/AR/AA ALT-informative pattern。1,064 個長 span W 中只有 4 個存在 adjacent membership gap >50 kb，顯示多數長 W 是由多段合格 read chain 連成，不是由 50 kb 幾何區間定義。三種距離統計不可互相代稱或合併成單一比例。

50 kb hard-edge-cutoff 反事實：47 條長邊分布在 22 個 W；移除後只有 4 個原始 W 的 linkage partition 改變，4 個都仍保留一個 k≥2 W，但各掉出一個 singleton membership。原始與反事實 W 總數同為 11,462，卻有 4 個 W 的 k/成員集合已改變；因此只比較 W 總數會漏掉 cutoff 對後續候選樹輸入的影響。

QC 提醒：chr6 與 chr16 合計占 upstream S 的 57.6229%，卻只占 active loci 的 9.2073%。這不破壞 region 守恆，但是論文引用 S→active 漏斗前，必須另外稽核 upstream VCF、exact PS 可用性與 fixed R/A callable 率，不可直接解釋為生物差異。

## 5. Edge threshold sensitivity

| minimum distinct molecules | qualifying edges | singleton abstain | tree-eligible k>1 regions |
|---:|---:|---:|---:|
| 1 | 117,760 | 19,238 | 13,210 |
| 2 | 91,090 | 24,931 | 12,160 |
| **3（primary）** | **76,202** | **28,384** | **11,462** |
| 5 | 58,214 | 33,385 | 10,290 |

提高 threshold 時 retained edge 必然單調不增、graph partition 只會維持或細分、singleton 必然不減；但 `k>1` region 數沒有一般性的單調保證，因為一個大 component 可能裂成兩個仍各自 `k>1` 的 component。本次 HCC1395 的 multisite region 數恰好隨 threshold 下降，這是觀測結果而非定理。Primary threshold=3 是可稽核的工程設定，不等同已完成最佳統計校準；論文需保留 sensitivity。

驗證層級需分開：threshold=3 已由 component、membership 與 edge TSV 逐列獨立重算；threshold 1/2/5 則由目前 22 份 chromosome receipts 重新加總並與 dataset summary 逐欄對帳，尚未在本報告層第二次重讀所有非 primary threshold rows。

## 6. k>12 的正確語意

90 個 k>12 component 仍是 90 個 source read-linked candidate regions；後續 k≤12 DP 只產生 computational blocks。這個名稱描述分析單位，不宣稱每個 component 就是一個真實 biological clone/subclone region。禁止把 block 數當 region 數。每個 block 必須保留 `source_component_id`，才能在 topology 與報告層回聚。

新版 partition 已把 28,384 singleton 全數排除：

- eligible source units：11,462，等於 tree-eligible k>1 regions。
- eligible unit memberships：34,267。
- partition 實際輸入的 `min(k)=2`；`k=1` unit=0；singleton component ID 洩漏=0。
- strict eligible component ID 與 partition unit component ID 逐一集合相等：11,462/11,462，missing=0，unexpected=0。
- bounded blocks：11,712。
- 其中 12 個 `k=1` child blocks 來自 k=13–153 的大型 source regions 切割，不是 source singleton 洩漏；12/12 的 retained weight/pattern count 為 0，adapter 必須並已依 `len(positions) < 2` 排除建樹。
- 因此 `W=11,462 source regions` 與 `B=11,712 computational blocks` 是不同單位；不可將 block 數作為 region 數或直接作為 topology 分母。
- structural constraint rows：57,629。
- constraint molecule weight：total 285,596 = retained 281,685 + cut 1,242 + unavoidable 2,669。
- Python/C++ k≤12 partition：chr1–22 共 22/22 PASS，mismatch=0。

這一層只完成 segmentation/constraint preservation；是否形成 tree input，還要通過 block 內 pattern MINREAD 與後續 solver。

## 7. 驗證證據

### 7.1 Python synthetic tests

輸入：`InterSubMod/tests/test_strict_endpoint_graph.py`、`InterSubMod/tests/test_build_strict_ps_hp_regions.py`。

命令：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  -m pytest -q tests/test_strict_endpoint_graph.py \
  tests/test_build_strict_ps_hp_regions.py
```

實際結果：`15 passed`。涵蓋 transitivity、`R-X-A` partial-read endpoint、不跨 PS、threshold 2/3、>50 kb、duplicate molecule fail-closed、singleton abstention 與 deterministic digest。

### 7.2 C++ tests

輸入：`InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/cpp/`。

命令：

```bash
bash research/20260723_production_exact_ps_strict_read_linkage/cpp/tests/run_tests.sh
```

實際結果：`PASS`；使用 `-std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror`，另經 ASan/UBSan probe 通過。

### 7.3 真實資料獨立驗證

- Python genome summarizer 重新讀 edge/membership TSV，對每個 component 再做 connected-component 驗證；39,846/39,846 全部通過。
- cross-PS violations=0；cross-HP violations=0。
- retained edge minimum support=3；RR+RA+AR+AA=total 全部守恆。
- chr22 strict graph Python/C++ semantic parity：1,600 observed edge rows、532 components，mismatch=0。
- chr1–22 strict graph Python/C++ semantic parity：22/22 PASS；117,760 observed pair rows、39,846 components，edge mismatch=0、component/role mismatch=0。

## 8. 七 technical datasets readiness

Canonical manifest 中 7/7 BAM、VCF、read-tag sidecar 與索引可讀，7 BAM `samtools quickcheck` exit=0。本輪開始時只有 HCC1395 22/22 chromosome 具最新 molecule artifacts；其餘 132 tasks 隨後已由 BAM 重抽完成。正式 full extraction 落點：

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/all7_production_v1/`

最終 154/154 extraction receipts、154/154 strict receipts 與 7/7 summaries 均 PASS；全資料集數字、圖表及完成層級請以 `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/` 為準。

## 9. Claim ceiling

目前可以主張：primary analysis unit 限制在 exact PS×primary HP family，且每個 tree-eligible candidate region 由 threshold-qualified direct read endpoint edges 連通；singleton 未進建樹；50 kb 不參與 membership 判定。

目前不可以主張：唯一 mutation-state topology、真實 clone tree、clone/subclone 數、祖先順序、跨 HP 的同細胞配對、跨 PS 的 HP orientation、單一 molecule 橫跨整個 W、VAF/CCF/CN/甲基支持或 phase accuracy。`k=1` 只表示該 exact-container membership 在指定 threshold 下沒有合格 edge，不表示該 locus 全域沒有 read support。RR-only edge 只證明共同 callable 的 reference-state read，不證明 somatic ALT 共現；RA/AR/AA 雖為 ALT-informative，仍不是 clone truth。所有 endpoint edges 都是無向 linkage，不是 evolutionary parent→child。七 technical datasets 實際代表六個 biological samples；HCC1395 與 HCC1395_DORADO 不得視為獨立 biological replicates。最新 completion audit 為 L1=7/7、production strict topology=0/7、clone/parent/fusion=0/7。
