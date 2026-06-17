<!--
建立: 2026-06-18
報告類型: 方法學參考 — read×read 甲基距離矩陣的「無監督分群驗證」方法選單 + ISM 推薦 stack + 單細胞甲基先例對照
任務類型: D handoff（供論文 Methods/Related Works 撰寫；給其他 AI/人查方法）
狀態: 文獻映射（L3）— method/citation 來自 workflow wdzxkbpwf（3 角度 WebSearch+WebFetch）+ 本地 ISM source 交叉驗證；落地 ISM 前須 pilot
data_sources:
  - workflow wdzxkbpwf（distance-matrix cluster-validation survey, 3 角度 researcher agents）
  - ISM source: src/core/{StructureTest,SignificanceAnalyzer,HierarchicalClustering,DistanceMatrix}.cpp（PERMANOVA pseudo-F+perm 已實測）
  - external_validation 卡: cvlr / asms / qfdrp / cpelnano / epiclomal / smallwood + 新增 pdclust / scmelody
provenance_note: 本檔不產新數字；method/citation 為文獻層（L3），引用前過 /citation-verification 補一手。ISM 現況 claim 有 source file 佐證。
-->

# read×read 甲基距離矩陣的無監督分群驗證方法 — 選單 + ISM 推薦 stack + 單細胞先例對照

> **要回答的問題**：找到 read×read 甲基距離後，怎麼**無監督**確認 (1) 分群真實存在、(2) 結構穩固、(3) 分幾群？並擴展到一般 distance-matrix 分群驗證。
> **一句話**：ISM 目前只有**存在性檢定（PERMANOVA，已是甲基領域 best-in-class）**，缺 **k-selection + 穩定性**；缺口**不能借甲基工具的 k-selection（全在 model/特徵空間，不可移植到距離矩陣）**，須用**距離矩陣原生**的通用驗證工具補（單細胞甲基領域 PDclust/scMelody/SINCEF 已例行於 cell-cell 距離上這麼做）。

## 0. ISM 現況（src 實測）
read×read 甲基距離（NHD/L1/Bernoulli/Jaccard，**precomputed pairwise，非 feature 空間**）→ 階層分群（UPGMA/Ward）→ **PERMANOVA pseudo-F + N-perm**（`StructureTest.cpp`）+ per-CpG Fisher。**缺**：(a) 無監督決定 k、(b) cluster 穩定性、(c) 內部 validity index。

## 1. 方法選單（★=可直接吃 ISM precomputed 距離矩陣）

| 方法 | 驗證什麼 | 吃距離矩陣? | 工具 | canonical citation |
|---|---|:--:|---|---|
| **PERMANOVA / adonis2** | 結構存在性 | ★ YES | R vegan::adonis2 / py skbio | Anderson 2001 |
| **betadisper (PERMDISP)** | 排除「離散度差異」confound | ★ YES | R vegan::betadisper | Anderson 2006 |
| **ANOSIM · MRPP** | 存在性 + 效應量(MRPP A) | ★ YES | R vegan::anosim/mrpp | Clarke 1993 / Mielke 2001 |
| **Mantel · partial Mantel** | 兩距離矩陣相關 | ★ YES | R vegan::mantel | Mantel 1967 |
| **Silhouette width** | k-selection + validity | ★ YES | cluster::silhouette(dist) | Rousseeuw 1987 |
| **Dunn index** | 內部 validity（無需質心）| ★ YES | clValid::dunn / fpc | Dunn 1974 |
| **clusterboot（Jaccard bootstrap）** | 穩定性（逐群）| ★ YES | fpc::clusterboot(distances=TRUE, disthclustCBI) | Hennig 2007 |
| **nselectboot** | k-selection（bootstrap）| ★ YES | fpc::nselectboot | Fang & Wang 2012 |
| **Consensus clustering + PAC** | k + 穩定性 | ★ YES | ConsensusClusterPlus / diceR::PAC | Monti 2003 / Şenbabaoğlu 2014 |
| **Eigengap (spectral)** | k-selection | ★ YES（距離→affinity）| sklearn SpectralClustering | von Luxburg 2007 |
| Gap statistic | k-selection | ⚠ ADAPT（null 預設特徵空間，須自建）| cluster::clusGap | Tibshirani 2001 |
| Prediction strength | k + 穩定性 | ⚠ ADAPT（須 nearest-medoid 分類器）| fpc::prediction.strength | Tibshirani-Walther 2005 |
| SigClust / SigClust2-SHC | 存在性 + 逐節點 split 顯著性 | ⚠ ADAPT（餵 raw read×CpG + Gaussian null）| sigclust / pkimes/sigclust2 | Liu 2008 / Kimes 2017 |
| Dip test | 單峰 vs 雙峰 | ⚠ ADAPT（投影 1-D：PCoA 軸或單 CpG β）| diptest::dip.test | Hartigan 1985 |
| Calinski-Harabasz · Davies-Bouldin | 內部 validity | ❌ NO（需質心，距離矩陣無向量空間）| sklearn | Caliński 1974 / Davies 1979 |

## 2. 推薦 ISM validation stack（次序鐵則：existence 先 gate → k → stability）
```
1. EXISTENCE (gate, 已有)：PERMANOVA pseudo-F + N-perm（續用）
   ＋ betadisper/PERMDISP 排除 dispersion confound
   ＋ 可選 ANOSIM/MRPP 三角驗證（MRPP 的 A 補效應量）
2. k-SELECTION (僅在第 1 層顯著後跑)：average silhouette over k（precomputed dist, 對角=0）
   次選 fpc::nselectboot（與 stability 共用 resampling）
3. STABILITY (對選定 k)：fpc::clusterboot 逐群 Jaccard
   （≥0.75 stable / ≥0.85 highly stable；小樣本改 continuous-weights/subsetting）
   可選 ConsensusClusterPlus + PAC（⚠ chance-partition，見 §4）
4. (可選 raw-matrix 補強)：SigClust2-SHC 逐節點 split 顯著性 / dip on PCoA 第一軸
```
> **為何這樣排**：existence gate 先判「這 region 有沒有結構」——ISM 自己 verdict 顯示 ~91.7% region 為 Weak/Noise（多數無結構），不先 gate 就在無結構區硬選 k = 製造假群。

## 3. 單細胞甲基先例對照（ISM 的「細胞層表親」）
這 3 個工具做「pairwise 甲基差異 → 距離矩陣 → 分群 → 矩陣上驗證」，跟 ISM 同骨架，差在**單位是細胞不是 read**。

| 軸 | PDclust | scMelody | SINCEF | **ISM** |
|---|---|---|---|---|
| 單位 | 細胞 | 細胞 | 細胞 | **read（區域內）** |
| regime | scBS 全基因組 | scBS 全基因組 | scBS 全基因組 | **ONT bulk 長讀 · 單 locus** |
| 距離/相似度 | PD(平均絕對差) | Cosine/Hamming/Pearson 整合 | 多距離融合 | NHD/L1/Bernoulli/Jaccard |
| 分群 | HC | spectral+HC | spectral(融合) | UPGMA/Ward HC |
| **存在性檢定** | ✗ | ✗ | ✗ | ✅ **PERMANOVA**（獨有）|
| k-selection | HC 切(手動) | ✅ silhouette 掃 k | ✅ silhouette/CH | ❌ 無 |
| 穩定性 | ✗ | ✅ consensus/NMI | △ | ❌ 無 |
| citation | Hui 2018 Stem Cell Reports (PMC6093082) | Tian 2022 Front Bioeng Biotechnol 10:842019 | IEEE BIBM 2021 (doc 9562895) | 本專案 |

**5 個關鍵差異**：
1. **read vs cell = 不同 regime**，不可直接等同（細胞層能重建 lineage ≠ read 層；守 §13 紅線）。
2. **🟢 ISM 獨有 existence-gate**（這 3 個都預設結構存在）→ 論文可寫「我方比單細胞甲基分群多一道存在性檢定」。
3. **🔴 ISM 缺 k+stability，這 3 個有** → 借 scMelody（silhouette-k+consensus）/ SINCEF（融合+eigengap），且皆距離原生可直接套。
4. **多距離融合**：ISM 算 4 個距離只用一個；**SINCEF 是融合先例**（ISM 可考慮融合，須先 pilot）。
5. **PDclust = ISM 直系結構雙胞胎**（PD 平均絕對差 → HC），只是 cell-level。

## 4. 🔴 誠實陷阱（落地前必守）
1. **假結構警告（對 ISM 最關鍵）**：Şenbabaoğlu 2014（Sci Rep 4:6207）證 consensus clustering 會把**無結構單峰隨機資料**切成「看似穩定」的群；silhouette 偏好凸/緊緻群；eigengap 僅群分明時可靠 → **任何 k>1 必先過 existence gate（PERMANOVA perm-null），兩者一致才採信**。
2. **小樣本**（每 region 數十~百 reads）：gap-null/bootstrap/split-half 在 n 小方差大 → bootstrap 改 continuous-weights（Möller-Radke 2006）；絕對閾值降級為「相對 + permutation 顯著性」。
3. **無單一決定性方法** → 必組合（existence→k→stability 三層各答不同問題）。
4. **甲基領域 k-selection 不可借**：cvlr 固定 k=2 / Epiclomal DIC / Melissa VB-ELBO 全在 model/特徵空間操作 raw read×CpG 矩陣，**不可移植到 precomputed 距離矩陣** → ISM 的 k 須用距離原生法（silhouette/consensus）。
5. **紅線一致**：此 stack 屬 structure **characterization**（確認甲基分群真實/幾群/多穩），**非 variant filter**（filter direction=DEAD）；甲基訊號處 germline-haplotype 層、非重建驅動。
6. **L3 文獻映射**：method/citation 已查（agent WebFetch），但 ISM 落地前須在真實 region pilot（尤其小樣本 + chance-partition 對照）；引用前過 `/citation-verification` 補一手。

## 5. 指引
- 完整方法調查原始：workflow `wdzxkbpwf`（3 角度）。上一層方法選單見本檔 §1；ISM 現況 source 見 §0。
- 單細胞先例卡（repo 外）：`external_validation/axis5_methylation_clustering_distance/{pdclust,scmelody}/`（2026-06-18 新增）。
- 同軸甲基對手既有卡：cvlr / asms / qfdrp / cpelnano（見 `external_validation/_landscape/05`）。
- memory: `project_distance_matrix_cluster_validation_methods`。
