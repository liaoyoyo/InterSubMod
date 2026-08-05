<!--
建立時間: 2026-06-21
類型: 觀察分析結果 (per-read 整合)
狀態: in_progress
主軸: Subclonal reconstruction — 甲基分群 × CN/SV 整合
plan: /bip7_disk/liaoyoyo2001/.claude/plans/ism-k-ism-zany-tiger.md
data_sources: _assets/aggregate_tp.json, _assets/aggregate_fp.json, _assets/region_table_tp.tsv, /big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna_normalhet/HCC1395_segmented_absolute_copy_number.tsv
build_branch: research/subclonal-reconstruction-202606
-->

# ISM 甲基分群 k_ISM × CN/SV per-read 整合分析（HCC1395 全基因組）

> **版本** v1.0 (2026-06-21) · 單樣本 **⭐3** · 證據 tier：**L1**=第一手磁碟數據重現
> 腳本 `InterSubMod/scripts/analysis/analyze_kism_vs_cn_perread.py`；數據 `_assets/`

## 0. 一句話結論

**ISM 甲基分群與 CN 大致無關**：k_ISM 與 k_CN **無相關**（Spearman ρ=−0.038）；全基因組只 **9.5%** somatic 區有顯著甲基 cluster 結構，且其中**僅 0.64% 對齊 germline-HP（CN/allele-dosage 可解釋）**、het 區 **85% unaligned**。→ **甲基分群不是 CN 鏡子**（有獨立內容），但存在的結構也**大多不對齊任何遺傳軸**（HP/ALT/SV）= 來源未定的 epigenetic 異質（bulk 下無法形式判 subclone）。per-read 群指派**可行**（345,714 read，n_clusters=optimal_k 98.2%）。

## 1. 方法（零重跑 ISM）

- k_ISM = per-region `clustering/significance.json::optimal_k`（UPGMA silhouette，floor=2）。**k_eff** = PERMANOVA gate 後真實群數：`cluster_structure.permanova_valid & permanova_p<0.05 & dispersion_p≥0.05` → optimal_k，否則 1。
- per-read 群標籤：從 `distance/BERNOULLI/matrix.csv`（read×read）**scipy 重算 UPGMA** 後 `fcluster(maxclust=optimal_k)`。⚠ C++ `linkage_matrix.csv` 的 cluster-id 方案與 scipy fcluster 不相容（`"uses the same cluster more than once"`）→ 改用同距離重算 UPGMA（同演算法，n_clusters==optimal_k 一致率 **98.2%**）。
- k_CN = SAVANA `cna_normalhet` 該 locus segment：minorCN≥0.5(het)→2 / <0.5(LOH)→1。
- 對齊 = cluster × {germline hp_family / somatic ALT-REF / SV-support} 列聯表 CramérV+Fisher（門檻 V≥0.7+p<0.05）。
- TP 29,754 區 + FP 627 區分跑（負對照）。

## 2. 四問結果（全基因組 TP，L1）

### Q1 — k_ISM vs k_CN：**無相關**
- Spearman **ρ=−0.038**（p=7.7e-11，僅因 n 大而「顯著」，效應量≈0）。
- k_eff 中位數：LOH=1、het=1（CN 狀態不改變甲基群數分佈）。
- crosstab：k_eff=1 在 k_CN=1（12,628）與 k_CN=2（14,296）皆主導，分佈近乎相同。

### Q3 — 多少甲基分群受 CN 影響：**極少（<1.4%）**
- 有顯著 cluster 結構：**9.5%（2,827/29,754）**（90% 無結構，呼應既有 91.7% Weak/Noise）。
- 結構區可用軸：hp 32.5% / alt 98.2%。
- 結構區對齊分類：

| class | n | % |
|---|---|---|
| unaligned | 1,119 | **39.6%** |
| candidate_beyond_CN（LOH 區 split 無 germline allele 可解釋）| 1,083 | 38.3% |
| aligned_somatic_allele_in_LOH | 401 | 14.2% |
| ambiguous | 139 | 4.9% |
| no_usable_allele_label | 50 | 1.8% |
| **CN_explained_HP（對齊 germline 單倍型）** | **18** | **0.64%** |
| aligned_somatic_allele（het）| 17 | 0.6% |

- **het-retained 結構區（1,317）**：CN_explained_HP **僅 1.37%**、aligned_somatic_allele 1.29%、**unaligned 85.0%**。
- **LOH 結構區（1,510）**：candidate_beyond_CN **71.7%**。
- → **甲基分群幾乎不由 CN/germline-allele 解釋**（非 CN 鏡子）；het 區結構絕大多數獨立於遺傳軸。

### Q2 — SV 斷點當 read 標籤：**稀疏 + 弱**
- SV-informative 區（≥3 tumour-support read）= 1,292。
- SV-support 對齊甲基 cluster **僅 5.1%**（同 subset HP 對齊 4.5%）→ SV 軸即使存在，對齊率也低。

### Q4 — per-read 群指派：**可行**
- 2,827 結構區產 **345,714 read 指派**（`perread_table_tp.tsv.gz`，欄：region_id/read_name/cluster/hp/hp_family/alt_support/is_tumor/sv_support/k_eff/k_CN/alignment_class）。
- n_clusters==optimal_k **98.2%**（驗證重建忠實）；cluster_purity wrt hp_family 0.867；median ARI(cluster,hp)=**0.0**（群不對應 germline 單倍型分區）。

### context（label-level，非 cluster-level）
- **88.8%** 區 allele-PERMANOVA 顯著（ALT-vs-REF 甲基差 = ASM 訊號）— 但這是 label 層的 allele-甲基差，**與 cluster 結構對齊是兩件事**（cluster 多 unaligned）。

## 3. FP 負對照（627 區）

- 結構率 **3.8%**（< TP 9.5% → **TP 區甲基結構約 2.5× 多於 FP**）。
- 同樣 **CN_explained_HP=0%**、candidate_beyond_CN 58%、no_usable 33%。→「不由 CN 解釋」對 TP/FP 皆成立；差別在結構「量」非「對齊型態」。

## 4. 詮釋與紅線

- **回答原始疑慮**：先前擔心「若甲基對齊 CN = 甲基只是 CN 鏡子（壞消息）」—— 實測**否定**此擔憂（<1.4% CN-explained）→ 甲基分群有獨立於 CN 的內容。
- **但翻面**：存在的結構**大多不對齊任何遺傳軸**（HP/ALT/SV）= 來源未定 epigenetic 異質。bulk 下**無法形式判定是 subclone / cis-ASM / 技術**（candidate_beyond_CN 僅「候選」，需 normal-cis control / multi-sSNV CCF / single-cell；I2 confound 天花板）。
- 🔴 **對外勿稱「對齊=subclone」或「甲基偵測 subclone」**；甲基在此框架是 characterize（有界）。
- 單樣本封頂 **⭐3**；跨樣本需 COLO829。

## 4b. 觀察圖（肉眼可看）

![k_ISM vs CN/SV 6-panel](figs/kism_vs_cn_observation.png)

`figs/kism_vs_cn_observation.png`（6-panel，`scripts/analysis/plot_kism_vs_cn_observation.py` 產）：① 結構率 TP 9.5% vs FP 3.8% ② k_eff 分佈不隨 k_CN（LOH/het 重疊）③ 對齊分類（CN_explained_HP 0.64%）④ cluster×HP vs ×ALT CramérV 散點（多落低-低）⑤ het 85% unaligned / LOH 72% candidate_beyond_CN 堆疊 ⑥ SV 軸 5.1% vs HP 4.5%。
觀察入口同步登錄於 `InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/OBSERVATION_INDEX.md`（甲基分群×CN/SV 整合觀察段）。

## 4c. 條件分析：結構區內，分群是否被 read 數 / CN 驅動（2026-06-22）

> 限定在「甲基判定有結構」的 2,827 區，問分群狀況是否其實是 coverage/CN 假象。script `structure_vs_readcount_cn.py`、data `_assets/structure_vs_readcount_cn.json`（L1）。

**結構區內（n=2,827）分群數 vs read數/CN — 全部無相關**：
| 對 | Spearman ρ |
|---|---|
| optimal_k vs n_reads（read 數）| **−0.047** |
| optimal_k vs SEQC2 CN | **−0.031**（NS）|
| optimal_k vs SAVANA copyNumber | **−0.025**（NS）|
| n_clusters vs n_reads / SEQC2 CN | −0.048 / −0.032 |
| 結構強度(-log10 PERMANOVA p) vs n_reads / SEQC2 CN | 0.007 / −0.001 |

**結構『存在與否』也不被 coverage/CN 驅動（全 29,754）**：結構區 n_reads median **117 < 無結構 124**（不增反略減）；SEQC2 CN median 3.00 = 3.00；結構率隨 SEQC2 CN **不升反略降**（CN≤3 ~10-12% → CN≥5 4-7%）。

→ **甲基「有結構」的分群狀況既非 coverage 假象、也非 CN 假象**（全 |ρ|<0.05）。CN 數量**不預測**甲基分群數 → 用 CN 偵測軟體「驗證數量」**無法當判別甲基結構的 discriminator**；其價值在 confound 排除（證明結構非 CN 造）+ LOH 區 mask + allele context（characterization 非 discrimination）。

## 5. Provenance / caveat
- **L1**：所有 % 由 `_assets/aggregate_{tp,fp}.json` + `region_table_tp.tsv` 第一手算（腳本 `analyze_kism_vs_cn_perread.py`，env cnvtools）。
- **caveat**：①k_eff 用 PERMANOVA gate（optimal_k floor=2 不可當「有結構」）②cluster 重建用 scipy 重算 UPGMA（C++ linkage 不相容 fcluster；98.2% 一致）③germline-HP 軸在 somatic 區僅 32.5% 可用（read 多 hp=3/0）④CN 參考=SAVANA cna_normalhet（LOH Jaccard 0.962 vs SEQC2）⑤candidate_beyond_CN 受 cis-ASM/技術 confound 未排除。
