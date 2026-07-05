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

> ✅ **CLOSE-LOOP：tumor BAM 1.18× 重複對三桶無影響（2026-07-02 讀碼驗證）**：census 用的 HCC1395 tumor BAM 含 ~18% read 兩份 primary（memory `reference_hcc1395_ont5khz_normal_bam_doubled`；chr19:1-2M 實測 distinct QNAME 6,402 中 1,140 出現兩次）。**但 `sm_linkage_genomewide.py` 的 co-read 計數是 QNAME-keyed**（L108 `calls[rn][pos]` + L221-223 `multi={rn:c for rn,c in calls.items()...}; coread=len(multi)`）→ 重複 QNAME 自動覆寫、coread 數 distinct QNAME，**1.18× 重複被程式自動去重**。三桶（linked=coread-based / isolated=position-based n_partners / underpowered）**全不受灌水**；VAF=alt/(ref+alt) 對稱加倍亦不變（CCF tier 穩）。→ 三桶 61/15.4/23.5% 為 canonical，**無需 dedup 重跑**。（唯一名目受影響 = 絕對 tumor depth ~1.18×，但不用於分桶。）

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

## §5 甲基裁決（chr17 實測 + cis-control pilot 06-28）

甲基 read 分群對齊 **cis-genotype 軸（α 23 ≫ lineage 6）** → **不能當 HP3 定相 / 分叉的獨立驗證器**（cis-confounded）→ 維持 bounded-auxiliary / off-ladder。最佳角色 = cluster-count sanity（甲基群數 > genetic → cis 過切警示）。詳 `20260628_methylation_probabilistic_auxiliary_framework_01.md`。

**🔴 cis-control pilot 裁決（06-28，`20260628_cis_control_scope_pilot_verdict_01.md`）**：matched-normal HP cis-control **僅對 CROSS-HP 區有效**（全資料 663/1874 = 35.4%；其中乾淨 CN-neutral+loss 僅 43 區）。SAME-HP 59%（subclone 在同一 germline HP 內，somatic-cis 主型）→ normal-HP 軸**正交、無法 control**（germline-ASM 因共用單倍型自動抵消、本非 confound）。**對接 624 needs_methyl：可分類 72、CROSS-HP 12（全 CN-gain multiplicity 混淆）、乾淨 CN-neutral = 0** → 「cis-control 解鎖甲基 Tier-3 給 624 區」**原假設不成立**。subclone-specificity 對 SAME-HP 多數區 = 結構性 UNDETERMINED（normal 無對應 within-HP 軸 = structural zero）。4-角度對抗驗證（workflow，全 high-conf）證實 + 新發現 **tumor/normal HP 標籤指向同一實體染色體**（同 phased VCF、240 het 98.85% 一致）。

**🔴 甲基拓樸「排序」pilot（06-28，verdict doc §11）**：分群❌/specificity❌ 後唯一可能用途 = 對基因型已定群提供排序/距離。實測（A_determined+SAME-HP+CN-clean，326→55 區）：甲基離 root 距離 vs 譜系深度 Spearman ρ≈0.18、permutation **p≈0.06-0.08 未達顯著**；distal(0.636)≥near(0.588) 非純 cis 痕跡但太弱；功率不足（淺樹為主）。→ **甲基不能當可信拓樸 resolver（L3 弱提示）**，僅可標軟提示。ISM = 存在性偵測器非排序器（read 無 genotype 標籤），但有 `compute_group_distances` + read-level BERNOULLI 距離；接骨幹 genotype 標籤可產群間距離原料，惟底層訊號弱、上限有限。

**🔴 甲基有界軟標記 pilot（06-29，verdict doc §12）**：(b) 隱藏次結構旗標——3,416 群測甲基雙峰+扣 HP/CN：unimodal 95.8%（可信負向篩選）、residual-candidate 46（1.35%），再排「區域內在雙穩態」12 + 只留 ALT 群 → **22 候選（0.64%）**；但這 22 = 無遺傳佐證的無監督甲基切分 = double-dip 路徑 → **L3 候選非確認**。(a) A_ambiguous 76 區→35 個 L3 軟提示。**裁決：甲基軟標記可當保守負向篩選+候選旗標（L3），不能確認群數/subclone**。產物 `methyl_auxiliary_annotation.json`。

## §6 驗證狀態 + 統計

> **📍 2026-07-01 correction（C2/C3/D4 修正落地；權威裁決 `20260701_ssnv_backbone_method_spec_and_correctness_audit_01.md`）**：下方 **topology_type 與 determinacy 已回填 07-01 canonical**（源 `_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json`，07-01 14:36 重生，本輪 Read-back 驗證）。✅ **robustness ladder（L2-L4）、及 7 樣本表** 的源檔（`backbone_stability_audit.json` / `multisample_summary.json`）**已於 07-02 重生**（`gen_multisample_summary.py` + `backbone_stability_audit.py`）→ 下方全部 07-02 canonical，本輪 Read-back 驗證。

- **拓樸型態 / structure**（7143 全區，pairwise；**07-04 canonical**）：single 1995 / branched(直系+姊妹) 1110 / **branched(需隱藏祖先) 23 + branched(順序未定,需隱藏祖先) 7〔07-04 新:深分支推斷祖先樹〕** / linear 741 / germline 371 / no-vec 2887 / **incompatible 9**〔07-04:39→9,30 區 row-laminar 誤殺→改建隱藏祖先樹〕。〔07-01：…/incompatible 39；修前 06-29：single 2018 / branched 1113 / linear 754〕
- **🔴 determinacy（canonical 分母 = 有 genotype 向量 3885 區；**07-04 gap#3 後**）**：A_determined **1764** / A_ambiguous **69** / B_pairwise 943 / C_underdetermined 544 / **〔recurrence 70 拆為 gap#3 m-通道:artifact 11 + candidate 9 + LOH_unresolved 50〕** / **incompatible 18** / other 477（sum=3885）。**region_coverage**（7143）= with_vector 3885 / germline_only 371 / no_vector 2887。determinacy % 分母 3885（45.4%=1764/3885）。
  - **gap#3 m-通道（Model A，權威 `20260704_incompatible_reclassification_hidden_node_finding_01.md`）**：70 recurrence 接整數 CN → **11 artifact(m>1;CN≥3)棄 + 9 candidate(m=1)留 + 50 copy-neutral LOH 未解(VAF L3)**；6 非-HCC1395 樣本 cn=unknown → m通道不可用。
  - **gap#2 identifiability（candidate_trees.json 上層，不改 determinacy census）**：69 A_ambiguous 全枚舉 → **44 收斂 determined + 25 真歧義**；完整 identifiability_table 分區 3885（determined 1808 / ambiguous 25 / recurrence 70 / conflict 18 / nonenumerable 1487 / other 477）。修正舊 `B1_equiprobable=0` 假「0 歧義」。
  - 〔07-01 C2/C3：A_determined 1741 / incompatible 118；修前 06-29：A_determined 1812 / incompatible 12〕
  - 🆕 **gap#1（partial-read subcube）2026-07-04 已落地 canonical（controlled 重跑，provenance 先確認 nic 0 不一致）**：`sm_multilocus_combinations` 加 subread_groups + BAM 重跑 → `topology_analysis` 對「無 full-cov 但 ≥2-位點 partial read」建 subcube-recovered。**HCC1395 救回 2403 區**（E_subcube_recovered：pairwise樹 1561 + 欠定 676 + 含衝突 166）→ **determinacy denominator 3885 → 6288**；🔴 **A_determined 1764 完全不變**（full-cov census 0 regression，救回區透明標 E_ 不僭稱 A_determined）。僅 HCC1395（6 樣本 subcube 未重跑，graceful 不受影響）。權威 `20260704_incompatible_reclassification_hidden_node_finding_01.md`。
- **🔴 incompatible 重分類（07-04；altset-row-laminar 過嚴 → pairwise 四型 `n_independent_clean` 分流；權威 `20260704_incompatible_reclassification_hidden_node_finding_01.md`）**：舊 determinacy incompatible **118 → 三分流**：① **30 救**（nic==0：pairwise 四型乾淨=perfect phylogeny 存在,row-laminar 誤殺）→ `build_hidden_node_tree` 建隱藏祖先深分支樹 → A_determined 23 + A_ambiguous 7；② **70 recurrence_required**（nic>0 & cn≠gain：非-gain 真四型 → Model A recurrence 候選,需獨立 m(CN)通道拆 artifact vs 真recurrence,弱主張）；③ **18 維持 incompatible**（精確 cross-tab：4 cycle-only〔非gain〕+ 6 gain-only〔無cycle〕+ 8 gain∩cycle；cn 分布 15 gain+3 loh；= CN-multiplicity artifact / pairwise 成環）。回歸乾淨（僅舊 118 區變動,0 誤改）。⚠ **對抗驗證(wf_520ddd16)揪出並修復 1 個 builder bug**：H1437 chr20 因 rescue-gate(read-level nic) 與 builder(population 向量) 資料源脫鉤 → 建出 homoplasy 假樹;已加 population 欄四型自檢修復（HCC1395 30 救全乾淨不受影響）。〔舊 07-01：C2 +106 判 118；has_cycle=22 為不同指標〕。處置：incompatible 標 likely-artifact 不建樹;recurrence 標「需 m 通道」;救回區為推斷祖先樹(標 inferred ancestral clone)。
- **真正穩固核心**（robustness ladder，`backbone_stability_audit.json` **07-04 重生**）：L1 A_determined **1764**(45.4%) → +TP-backed **910**(23.4%) → +多位點 **239**(6.2%) → **+非CN-gain = 44 區(1.1%)**〔07-01：1741/889/218/34；修前 06-29：1812/931/260/57〕；core 80.1% 僅2位點 / 70.7% CN-gain / 43.6% 無truth → 「45.4%」是軟上界。⚠ 救回的 23 區含**推斷祖先(inferred ancestral clone)**,其中 +10 進真正穩固核心 L4(34→44)→ 論文須標此 caveat。詳 `20260629_backbone_resolution_and_methyl_roles_final_01.md`。
- **可辨識度**（Q1/Q2）：穿越充分 12.3% / 拓樸可辨識 10.9%；欠定根因 跨HP 36%/幾何 32%/功率 26%/HP3 6%。
- **chr17 worked**：S2(α,VAF0.82,TP,H1) 祖先 → S1+S3(β,co_linked) 後代；linear；dropped 2 noise；16 sig CpG。
- **驗證抓修 2 bug**：population 噪聲過濾（chr17 假 incompatible）+ ambiguous-parentage（chr14 缺中間群）。

## §7 誠實邊界 / open gates

1. ⭐3 單樣本 single-pipeline 封頂；升 ⭐4 需 ≥5/7 + COLO829。
2. regional(≤read-span) 非 genome-wide tree；分子共現 ≠ single-cell。
3. 拓樸可信僅 ~11%（determined）；其餘相容但欠定。
4. ~~**#1 load-bearing：T-GATE-GB matched-normal cis-control** → 解鎖甲基 Tier-3 + HP3 定相 + subclone-specificity。~~ **✅ 06-28 已測 → 結論：條件式適用（僅 CROSS-HP 35.4%），對 624 needs_methyl 乾淨可用 ≈ 0**（§5 + `20260628_cis_control_scope_pilot_verdict_01.md`）。甲基維持 bounded-auxiliary；subclone-specificity 對 SAME-HP 59% 為結構性 UNDETERMINED（需 single-cell/multi-region）。**殘餘可做 = 43 個乾淨 CROSS-HP 區的 PoC**（非解 needs_methyl）。
5. 其他 gate：T-METHYL-REEXTRACT 深覆蓋（救 underpowered）/ Tier-PS（救 isolated same-PS）/ mappability mask（清偽影）。

## §8 檔案 / 格式索引

**docs/methodology/**（commit 8339605/b51e071/58e01c7 + 本檔）:
- 20260628_subclone_reconstruction_master_spec_01.md（本檔，總索引）
- 20260628_lineage_label_definition_01.md（HP{h}-path 定義）
- 20260628_sSNV_linkage_threshold_decision_eps2_01.md（ε=2% ADR）
- 20260628_methylation_probabilistic_auxiliary_framework_01.md（甲基 Tier-3）
- 20260628_cis_control_scope_pilot_verdict_01.md ⭐（cis-control pilot 結果 + 適用 scope 裁決；甲基 bounded-auxiliary 機制證據）
- 20260628_reconstruction_model_verification_01.md（模型驗證）
- 20260627_subclone_4axis_teaching.standalone.html（4 軸教學）
- _assets/20260629_multisample_topology_workstation.standalone.html（**多樣本互動工作站＝主結果**：7 樣本分頁/篩選/排序/克隆樹/chr17/宇宙帳本/四配子衝突框/逐區 SAME-CROSS-HP 甲基判定。**取代**舊單樣本 20260628_topology_workstation.standalone.html〔已 deprecated，build script 不再產生〕）

**_assets/20260627_subclone_4axis_teaching/**:
- `data/`：sm_*.json(凍結真值) · regions.tsv · per_sSNV_census.tsv · lists/*.tsv(per-pair) · topology_per_region.json · single_snv_accounting.json · 各分析 JSON
- `scripts/`：topology_analysis.py · build_topology_workstation.py · single_snv_accounting.py · apply_eps2_canonical.py · build_lineage_labels.py · 等（皆可重跑複算）
- `PROVENANCE.md`：每數字→來源 + branch 凍結紀錄

## §9 候選評分 + 確認佇列（衝突/2-3 順位/需輔助驗證）

對每區算 **confidence_score(0-100)** = base(determinacy) + 覆蓋 + CN-clean + 單HP − FP主導 − 順序未定 − 偽影區(centromere/telomere)；並標 **resolution_path**（需什麼資訊定先後/驗證）：

| situation | 區數 | 評分傾向 | resolution（需哪些資訊）|
|---|--:|---|---|
| 已確定 | 2,198 | 高 | genetic 足夠（單分子向量；含 07-04 救回隱藏祖先深分支樹 +20）|
| pairwise 拼接 | 892 | 中 | 需更長 read / 單分子整跨 |
| 多樹相容(欠定) | 544 | 低 | 加深覆蓋 / Tier-PS 連遠 partner |
| 跨HP(兩棵樹) | 94 | 中 | HP 已拆，各樹獨立 |
| **順序 2-3 順位待定** | 69 | 中低 | **VAF/CCF(CN-clean) 定先後；🔴 甲基不用於排序**（含 07-04 救回的 7 個順序未定深分支樹；D4 已移除甲基 tie-breaker；06-30 政策「排序絕不用甲基」）|
| **recurrence(Model A m-通道已拆)** | **70** | 低 | **07-04 gap#3**：非-gain 四型 → 整數 CN m-通道拆 → **11 artifact(m>1;CN≥3)棄 + 9 candidate(m=1)留 + 50 copy-neutral LOH 未解(VAF L3,不硬判)**；6 樣本 cn=unknown 不可用 |
| 衝突(成環) | **18** | 極低 | likely-artifact（CN-gain multiplicity / pairwise 成環；補 mappability mask、不強建樹）|

> 〔**07-04 canonical**，源 `candidate_scoring.json`（07-04 重生,含 incompatible reclassification）；07-01：已確定 2178 / 跨HP 91 / 順序 62 / 成環 118；修前 06-29：已確定 2245 / 成環 12〕

- **需確認佇列 2,185 區**（07-01 canonical；非已確定或 <70 分），最低分為 chr9:41.8M/chr14:16M/chr16（已知 dense/centromere 偽影，自動排前）。
- **需甲基輔助 605 區**（VAF tie 或欠定）→ 原以為是甲基 Tier-3 適用處，但 **06-28 cis-control pilot 證實乾淨可用 ≈ 0**（72 可分類中 CROSS-HP 12 全 CN-gain、neutral 0）→ 甲基**無法**作這些區的 resolver（§5）。
- 互動確認：`topology_workstation` 下方「確認佇列」可**左右選項判讀**（✓同意rank1 / ⇄偏好其他 / ?需更多資訊）+ 觀察評分，存 localStorage、可匯出 JSON。
- 產物：`scripts/candidate_scoring.py` + `data/candidate_scoring.json`。

> 一句話：本研究 = 在 HCC1395 ONT 上，用 **sSNV 單分子共現（非循環骨幹）+ HP 定根 + 甲基有界輔助** 重建區域級克隆樹；35,332 sSNV 中 61% 可建樹、~11% 拓樸可辨識、2,185 區待確認（605 區可能需甲基輔助、待 cis-control）；單位點非全無法處理（有 CCF/Tier-PS）；甲基 cis-confounded 待 cis-control；⭐3 單樣本。工作站自帶 28 詞名詞解釋 + 左右確認評分。
