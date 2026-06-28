<!--
建立時間: 2026-06-28
類型: master spec — sSNV 克隆/亞克隆重建：完整資訊帳本 + 標籤格式 + 算法 + 證據階層(整合本 session 全部)
狀態: in_progress(方法定稿;甲基/HP3 待 T-GATE-GB cis-control;⭐3 單樣本)
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data
provenance: 凍結資料 @ branch feat/summary-nreadsvalid@5308d9e;本 session commits 8339605/b51e071/58e01c7
-->

# sSNV 克隆/亞克隆重建 — Master Spec（完整資訊、格式、算法、證據階層）

> 框架：Diátaxis(reference)。本檔為整個研究的單一索引；每數字 grep-able。HCC1395 ⭐3 單樣本。

## §1 全 sSNV 宇宙帳本（無遺漏，整體各情況比例）

**35,332 sSNV = TP 30,490 + FP 4,842**（SEQC2 標；TP/FP 只觀察、不進前處理）。三桶（per_sSNV_census）：

| 桶 | n | % | 狀態 / 是否有資訊 |
|---|--:|--:|---|
| **linked 可建樹** | 21,554 | 61.0% | 進拓樸分析（§3）|
| **underpowered** | 5,458 | 15.4% | 有 partner 無共讀 link → **加深覆蓋可救**；**有 census VAF→CCF 可刻畫**（clonal 2992 / mid 2174 / low 292）|
| **isolated** | 8,320 | 23.5% | read-span(≤50kb)內無 partner（Tier-R 樹外）；**有 caller VAF 可刻畫 + 可能 same-PS(>50kb) partner（Tier-PS 未做）** |

🔑 **單位點非「全無法處理」**（回答研究問題）：
- underpowered（15.4%）= 有 CCF + 深覆蓋可救 → **有資訊可用**。
- isolated（23.5%）= 此共現 census 無 VAF，但 caller(ClairS) VCF 仍有 VAF（可放 clonal 譜）+ 可能 Tier-PS 連鎖 → **部分可救、非永久 dead**。
- **真正拓樸-dead** = isolated 中無 same-PS partner 者（待 Tier-PS 細分量化）。
- by_source：TP|linked 17983 / TP|underpowered 4952 / TP|isolated 7555；FP|linked 3571 / FP|underpowered 506 / FP|isolated 765。

## §2 標籤與格式定義（一致命名）

| 物件 | 標籤 | 定義 |
|---|---|---|
| somatic sSNV | **S1, S2, …, Sk** | 區內依座標排序的 sSNV（genotype 向量第 i 位 = Si）|
| read 群 | **r / population** | 同基因型向量的 read 集（worked example 用 r1.. per-read）|
| 甲基位點 | **m1, m2, …** | 顯著差異 CpG（chr17: m1-m16）|
| haplotype 根 | **H1 / H2 / H3?(unphased)** | germline HP tag（"1-1"→H1、"2-1"→H2、"3"→H3?）|
| lineage 標籤 | **HP{h}-{b1}-{b2}…** | Dewey 路徑（root→分支；VAF 遞減編號）|
| 基因型向量 | **R/A 串**（如 RAR）| 位置=Si，A=帶突變/R=不帶 |
| 關係（V/H）| **vertical(直系/nested) / horizontal(姊妹/sibling) / co_linked(同節點)** | 由 2×2 AA 格定（有 AA=直系、AA=0 same-HP=姊妹）|
| situation tier | **A 單分子整跨 / B 可整跨pairwise / C 必鏈接 / HP-mixed / 不相容** | span vs read 長 + multi-support |
| determinacy | **A_determined / A_ambiguous(缺中間群) / B_pairwise / C_underdetermined / incompatible** | 連通+相容+單分子 |
| genome_ctx | **telomere / centromere / arm** | region vs hg38 chrlen+centromere（±3Mb 近似）|
| 來源 | **TP / FP**（SEQC2）| 僅觀察、不進前處理 |

## §3 算法（cluster-count-first topology）

```
每區: ① n_clusters(genotype 向量) → ② HP 定根(1850 兩根/3433 一根) →
③ 突變集 laminar 解樹(perfect-phylogeny;noise 過濾貪婪移除最小衝突群;ambig 偵測=跳>1 突變/缺中間群) →
④ V/H 分型 → ⑤ determinacy + HP{h}-path 標籤
約束: n_pop ≤ k+1(99.9%),拓樸 ∈ c-節點樹 ≪ 2^k(中位 c=2)
```
理論：二元字元 perfect-phylogeny → **pairwise 相容即足以定樹**（不需單分子整跨）。實作見 `topology_analysis.py`。

## §4 證據階層（3-tier，排序固定）

```
P(拓樸 T) ∝ P(sSNV共現|T) × 1[HP-root 相容] × P(甲基|T,normal-baseline)^w
            Tier1 決定性,非循環   Tier2 硬約束        Tier3 prior,弱權,cis-control-gated
```
- **Tier1** sSNV 單分子共現（AA 格定 V/H、零格定方向、ε=2%）= 唯一非循環骨幹。
- **Tier2** HP tag 定根（跨-HP 拆獨立樹）。
- **Tier3** 甲基（normal-baseline）= 機率輔助，**validation 待 T-GATE-GB**。

## §5 甲基裁決（chr17 實測）

甲基 read 分群對齊 **cis-genotype 軸（α 23 ≫ lineage 6）** → **不能當 HP3 定相 / 分叉的獨立驗證器**（cis-confounded）→ 維持 bounded-auxiliary / off-ladder。最佳角色 = cluster-count sanity（甲基群數 > genetic → cis 過切警示）。詳 `20260628_methylation_probabilistic_auxiliary_framework_01.md`。

## §6 驗證狀態 + 統計

- **拓樸型態**（7143 區）：single 2018 / branched 1113 / linear 754 / germline 371 / no-vec 2887。
- **determinacy**：A_determined 1812 / A_ambiguous 76 / B_pairwise 2760 / C_underdetermined 1658 / incompatible 22。
- **可辨識度**（Q1/Q2）：穿越充分 12.3% / 拓樸可辨識 10.9%；欠定根因 跨HP 36%/幾何 32%/功率 26%/HP3 6%。
- **chr17 worked**：S2(α,VAF0.82,TP,H1) 祖先 → S1+S3(β,co_linked) 後代；linear；dropped 2 noise；16 sig CpG。
- **驗證抓修 2 bug**：population 噪聲過濾（chr17 假 incompatible）+ ambiguous-parentage（chr14 缺中間群）。

## §7 誠實邊界 / open gates

1. ⭐3 單樣本 single-pipeline 封頂；升 ⭐4 需 ≥5/7 + COLO829。
2. regional(≤read-span) 非 genome-wide tree；分子共現 ≠ single-cell。
3. 拓樸可信僅 ~11%（determined）；其餘相容但欠定。
4. **#1 load-bearing：T-GATE-GB matched-normal cis-control** → 解鎖甲基 Tier-3 + HP3 定相 + subclone-specificity。
5. 其他 gate：T-METHYL-REEXTRACT 深覆蓋（救 underpowered）/ Tier-PS（救 isolated same-PS）/ mappability mask（清偽影）。

## §8 檔案 / 格式索引

**docs/methodology/**（commit 8339605/b51e071/58e01c7 + 本檔）:
- 20260628_subclone_reconstruction_master_spec_01.md（本檔，總索引）
- 20260628_lineage_label_definition_01.md（HP{h}-path 定義）
- 20260628_sSNV_linkage_threshold_decision_eps2_01.md（ε=2% ADR）
- 20260628_methylation_probabilistic_auxiliary_framework_01.md（甲基 Tier-3）
- 20260628_reconstruction_model_verification_01.md（模型驗證）
- 20260627_subclone_4axis_teaching.standalone.html（4 軸教學）
- 20260628_topology_workstation.standalone.html（互動工作站：篩選/排序/樹/chr17/宇宙帳本）

**_assets/20260627_subclone_4axis_teaching/**:
- `data/`：sm_*.json(凍結真值) · regions.tsv · per_sSNV_census.tsv · lists/*.tsv(per-pair) · topology_per_region.json · single_snv_accounting.json · 各分析 JSON
- `scripts/`：topology_analysis.py · build_topology_workstation.py · single_snv_accounting.py · apply_eps2_canonical.py · build_lineage_labels.py · 等（皆可重跑複算）
- `PROVENANCE.md`：每數字→來源 + branch 凍結紀錄

## §9 候選評分 + 確認佇列（衝突/2-3 順位/需輔助驗證）

對每區算 **confidence_score(0-100)** = base(determinacy) + 覆蓋 + CN-clean + 單HP − FP主導 − 順序未定 − 偽影區(centromere/telomere)；並標 **resolution_path**（需什麼資訊定先後/驗證）：

| situation | 區數 | 評分傾向 | resolution（需哪些資訊）|
|---|--:|---|---|
| 已確定 | 2,245 | 高 | genetic 足夠（單分子向量）|
| pairwise 拼接 | 903 | 中 | 需更長 read / 單分子整跨 |
| 多樹相容(欠定) | 550 | 低 | 加深覆蓋 / Tier-PS 連遠 partner |
| 跨HP(兩棵樹) | 99 | 中 | HP 已拆，各樹獨立 |
| **順序 2-3 順位待定** | 76 | 中低 | **VAF/CCF(CN-clean) 定先後；VAF tie→甲基輔助(cis-control 後)** |
| 衝突(成環) | 12 | 極低 | likely-artifact（補 mappability mask、不強建樹）|

- **需確認佇列 2,118 區**（非已確定或 <70 分），最低分為 chr9:41.8M/chr14:16M/chr16（已知 dense/centromere 偽影，自動排前）。
- **需甲基輔助 624 區**（VAF tie 或欠定）→ 這正是甲基 Tier-3 機率層的適用處，但**待 T-GATE-GB cis-control 才可信**。
- 互動確認：`topology_workstation` 下方「確認佇列」可**左右選項判讀**（✓同意rank1 / ⇄偏好其他 / ?需更多資訊）+ 觀察評分，存 localStorage、可匯出 JSON。
- 產物：`scripts/candidate_scoring.py` + `data/candidate_scoring.json`。

> 一句話：本研究 = 在 HCC1395 ONT 上，用 **sSNV 單分子共現（非循環骨幹）+ HP 定根 + 甲基有界輔助** 重建區域級克隆樹；35,332 sSNV 中 61% 可建樹、~11% 拓樸可辨識、2,118 區待確認（624 區可能需甲基輔助、待 cis-control）；單位點非全無法處理（有 CCF/Tier-PS）；甲基 cis-confounded 待 cis-control；⭐3 單樣本。工作站自帶 28 詞名詞解釋 + 左右確認評分。
