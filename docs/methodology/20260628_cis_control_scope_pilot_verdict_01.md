<!--
建立時間: 2026-06-28
類型: methodology verdict — T-GATE-GB matched-normal cis-control pilot 結果 + 全資料集適用 scope 裁決
狀態: concluded(pilot 設計驗證完成;結論 = 條件式適用 + needs_methyl ∩ 乾淨可用 ≈ 0)
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/cis_control_scope_summary.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/pilot_cis_control.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/hp_alignment_fullscan.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/methyl_ordering_pilot.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/methyl_auxiliary_annotation.json
provenance: HCC1395 tumor.bam+normal.bam(big8) + 凍結 topology_per_region.json @ branch feat/summary-nreadsvalid@5308d9e;每數字 grep-able 於 _assets/.../data/*.json
-->

# T-GATE-GB matched-normal cis-control — Pilot 結果 + 適用 Scope 裁決

> 框架：Verdict-Pyramid（結論先行 → 證據 → 邊界）。HCC1395 ⭐3 單樣本。**每個數字標來源檔**（§13-A）。
> 一句話：matched-normal HP cis-control **不是壞掉，是只對 CROSS-HP 區（全資料 35.4%）有效**；而真正需要甲基輔助的 624 區裡，乾淨可用的 ≈ **0**。甲基維持 bounded-auxiliary，subclone-specificity 對 SAME-HP 多數區（59%）為**結構性 UNDETERMINED**。

---

## §1 TL;DR（裁決）

| 問題 | 裁決 | 證據級 |
|---|---|---|
| cis-control 設計可跑嗎？ | ✅ 可跑（8 區 pilot、2,073 CpG 對） | L1 |
| corr≈0 是「無 cis-confound」好消息嗎？ | ❌ **不是**——是 tumor genotype 軸與 normal HP 軸**正交**的必然 | L2 |
| matched-normal HP 基線何時有效？ | **僅 CROSS-HP 區**（兩 subclone 群跨不同 germline HP）= 全資料 663/1874 = **35.4%** | L2 |
| 多數區（SAME-HP 59%）能做嗎？ | ❌ 結構性不行（normal 無對應 within-HP 軸）；但 germline-ASM 因共用單倍型**自動抵消**故也非 confound | L2 |
| 能解鎖原本 624 個 needs_methyl 區嗎？ | ❌ **幾乎不能**——624 中可分類 72 區、CROSS-HP 僅 12（全 CN-gain）、乾淨 CN-neutral = **0** | L2 |
| 對研究敘述的意涵 | 甲基維持 **bounded-auxiliary / characterization-only**；「cis-control 解鎖甲基 Tier-3」的原假設**不成立** | L2 |
| 甲基能否當**拓樸排序 resolver**（§11 延伸 pilot）| ⚠️ **弱提示、未達顯著、功率不足**（rho≈0.18, perm p≈0.06-0.08）→ 不能當可信 resolver | L3 |

---

## §2 設計回顧 — cis-control 原本要做什麼

**動機**（master spec §7 #1 load-bearing gate）：原假設「matched-normal cis-control 可解鎖甲基 Tier-3 + subclone-specificity」，用於 624 個 needs_methyl 區（VAF tie / 欠定）的拓樸輔助。

**Pilot 設計**：每區
- **tumor** 依 genotype cluster（somatic sSNV 向量）分群 → 取最大兩群的 per-CpG **Δβ_tumor**（subclone 候選訊號）；
- **normal** 依 germline **HP tag** 分群（H1 vs H2）→ 同 CpG 的 per-CpG **Δβ_normal**（germline-ASM 基線）；
- 比較：`Δβ_tumor 高 + Δβ_normal 高` = germline-ASM（非 subclone）；`Δβ_tumor 高 + Δβ_normal≈0` = 候選 subclone-specific。

腳本：`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/pilot_cis_control.py`

---

## §3 Pilot 觀察（分布，先不設門檻）

8 區 pilot、6 區有 CpG、**2,073 CpG 對**（`pilot_cis_control.json`）：

| 量 | tumor (genotype 軸) | normal (HP 軸) |
|---|---|---|
| median \|Δβ\| | 0.047 | 0.044 |
| p90 \|Δβ\| | 0.229 | **0.501** |
| max \|Δβ\| | 0.966 | 0.979 |

聯合分布（`cis_control_scope_summary.json` → `joint_distribution`）：
- tumor `\|Δβ\|≥0.2`（有訊號）= 256 / 2073 = **12.3%**
- 其中 normal 也高（germline-ASM）= 68（27%）
- normal 低 `<0.1`（若軸對齊則為 subclone 候選）= 138（54%）
- **corr(Δβ_tumor, Δβ_normal) = −0.026**（近零）

散點圖：`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/pilot_cis_control_scatter.png`（點雲對稱、無對角趨勢）。

⚠ **這 6 區全是 LOH**（見 §4）→ 此散點圖主要展示「正交軸」而非通則，不可單憑它宣稱「54% 是 subclone-specific」。

---

## §4 關鍵機制 — 為什麼 corr≈0：軸對齊問題

corr≈0 有兩種相反解讀。直接查 tumor genotype 群的 **HP tag 組成**（`pilot_axis_alignment.json`）判定：

> **解析陷阱（已修）**：tumor BAM 的 HP tag 是**字串**（`"1-1"`/`"1"`/`"2"`），normal 是**整數**。初版用整數 key 比對字串 → 全抓不到 → 假的「MIXED」。修正（取 `-` 前數字正規化）後結論才正確。

**Pilot 6 區（全 LOH）→ 6/6 SAME-HP**：每區兩 subclone 群都落在**同一條 germline 單倍型**。但 LOH 只剩一條 HP，SAME-HP 是必然 → 選區偏差。

**補測 CN-neutral 區（兩條 HP 都在，`axis_alignment_neutral.json`）→ 10 區中 5 SAME-HP / 5 CROSS-HP**：證實 SAME-HP 不是純 LOH artifact，但也非全部。

**機制定義**：

| 類型 | 機制 | normal-HP cis-control |
|---|---|---|
| **CROSS-HP**（如 chr21: AR在H1、RA在H2）| 兩突變在**不同**親代染色體（獨立 somatic 事件）| ✅ **有效**——subclone Δβ 含 germline-HP 成分，normal 可扣除 |
| **SAME-HP**（如 chr2:216540437: AA與RA都在H1）| 突變在**同一**染色體的 within-HP lineage（nested/sibling）| ⚠ normal-HP **正交、無法 control**；但兩群共用 germline 單倍型 → germline-ASM **自動抵消**（本來就非 differential confound）；殘差 = somatic-cis（normal 結構上無對應軸）|

🔑 **生物學合理性**：somatic 突變是單一親代染色體上的 cis 事件 → within-HP lineage（SAME-HP）是預期主型；CROSS-HP 發生於兩條染色體各自獨立突變。

---

## §5 全資料集適用 Scope（權威交叉表）

對全部 **1,874 個 ≥2-ALT-群的區**（1890 目標，skip 16）分類（`hp_alignment_fullscan.json`）：

| CN ＼ alignment | CROSS-HP | SAME-HP | MIXED | 小計 |
|---|--:|--:|--:|--:|
| gain | 610 | 615 | 75 | 1,300 |
| loh | 10 | 467 | 27 | 504 |
| neutral | 42 | 17 | 4 | 63 |
| loss | 1 | 6 | 0 | 7 |
| **總計** | **663 (35.4%)** | **1,105 (59.0%)** | **106 (5.7%)** | **1,874** |

**讀法**：
- **CROSS-HP 663（35.4%）= cis-control 有效**——但其中 610 是 **CN-gain**（受 multiplicity 混淆：一個「群」可能是一份 copy 非一個 cell），真乾淨僅 neutral 42 + loss 1。
- **SAME-HP 1,105（59.0%）= 正交**——含 467 LOH（93% 的 LOH 區）。
- **最乾淨 cis-control 適用集**（`cn∈{neutral,loss} ∧ CROSS-HP`）= **43 區**（`cis_control_scope_summary.json` → `cleanest_applicable_set`）。

---

## §6 對接 needs_methyl — 嚴峻的 scope 發現

cis-control 原本要解鎖的是 **624 個 needs_methyl 區**（candidate_scoring：VAF tie / 欠定）。對接結果（`needs_methyl_intersection`）：

| 項 | 數 |
|---|--:|
| needs_methyl 總數 | 624 |
| 其中有 ≥2 ALT 群、可分類 | **72**（其餘 552 區無清楚 ≥2 ALT 群，是因單群/欠連結而需助）|
| ↳ SAME-HP | 54（75%）|
| ↳ CROSS-HP | 12（17%）|
| ↳ MIXED | 6 |
| ↳ CN 分布 | gain 48 / loh 24 / **neutral 0** |
| **CROSS-HP 且 cis-control 可助** | 12（全 CN-gain，multiplicity 混淆）|
| **CROSS-HP 且 CN-neutral（乾淨可用）** | **0** |

🔴 **核心發現**：`{需要甲基輔助}` ∩ `{cis-control 乾淨可用}` ≈ **空集**。
- 需要甲基幫忙的區，75% 是 SAME-HP → 甲基與 genotype **共分離（co-segregate）** → 對拓樸**不提供獨立資訊**（genotype 已定 partition）。
- 少數 CROSS-HP（12）全是 CN-gain → multiplicity 混淆。
- 反之，43 個最乾淨可用區，是**已定、不需甲基幫忙**的區。

**直白結論**：原假設「cis-control 解鎖甲基 Tier-3 給 624 區」**不成立**。不是執行失敗，是結構性質——需要幫忙的區，正好是甲基無法獨立幫忙的區。

---

## §7 裁決 + 對論文敘述的意涵

1. **甲基維持 bounded-auxiliary / characterization-only**（與既有敘述一致，但現有機制證據）：
   - SAME-HP 區：甲基 Δβ 與 somatic genotype 共分離 → 描述 subclone 甲基狀態 OK，但**不增加重建力**（partition 已由 genotype 決定）。
   - subclone-specificity 對 SAME-HP 多數區 = **結構性 UNDETERMINED**（不是「已驗陰性」，是 normal 無對應 within-HP 對比軸 = structural zero）。這正是 memory `project_subclone_snv_linkage_verification_pipeline` 記「740 structural-zero → UNDETERMINED」的**機制原因**。
2. **論文 A-framing 獲強化**（reconstruction 歸 sSNV 骨幹 / 甲基 characterize 有界）：pilot 直接數據顯示甲基沿 cis-genotype 軸、normal 無法獨立驗證 → 不可寫成「甲基驅動重建」或「cis-control 已驗證 subclone-specific 甲基」。
3. **silver lining**：SAME-HP 的共用單倍型使 germline-cis-ASM **自動抵消**（兩 subclone 群同 germline 序列）→ subclone Δβ **不是** germline-ASM artifact；殘差純為 somatic-cis（突變的局部足跡）。這對「甲基差異真實存在」是正面的，但無助於證明「subclone-specific lineage marker」。

---

## §8 誠實邊界 + 殘餘可做

**做不到的（結構性，非覆蓋問題）**：
- 用 matched-normal 證明 SAME-HP 區甲基差異是 subclone lineage marker vs somatic-mutation-cis 後果 → 兩者在單一 bulk 樣本**不可識別**（normal 無 somatic 突變，無法重建 within-HP subclone 軸）。需 single-cell / multi-region（Tarabichi LEARN：單 bulk 只能 characterization）。

**可做的（小但乾淨）**：
- **43 個 CN-neutral/loss CROSS-HP 區**可做 cis-control-validated subclone-specific 甲基的 **proof-of-concept**（扣除 normal-HP germline-ASM 後殘差 = somatic 候選）。範例：chr2:204644837-204648404、chr21:22306211-22315240 等（`cleanest_applicable_set.examples`）。⚠ 這是**展示方法可行性**，非解決 needs_methyl 區。
- 甲基最佳角色維持 master spec §5 的 **cluster-count sanity**（甲基群數 > genetic 群數 → cis 過切警示），不作 resolver。

**門檻問題（用戶原問「先觀察分布再定」）**：→ **無需定門檻**。因為決定性瓶頸不是「Δβ 門檻多少」，而是「軸是否對齊」。在不對齊（SAME-HP）的多數區，任何 Δβ 門檻都無意義；在對齊（CROSS-HP）的乾淨少數區，再談門檻（建議跟既有 ASM 口徑 `|Δβ|≥0.2` + permutation null，但樣本太少先不固化）。

---

## §9 REFLECTION（給下次同方向 agent）

**警示指標**：若再提「用 matched-normal cis-control 解鎖甲基給欠定區做拓樸 resolver」→ 停。本輪已證 `{需助} ∩ {可用} ≈ 0`。

**根因**（double-loop）：原假設的隱含前提是「subclone 分群 ⟂ germline 甲基 confound，可用 normal 扣除」。但 subclone 是 somatic-cis 事件 → 多在**同一條 germline HP 內**（SAME-HP 59%）→ normal 的對比軸（HP1 vs HP2）與 subclone 軸正交，扣不到；且本就無 germline-ASM 可扣。前提錯，非方法錯。

**改進方向（若重試）**：需 (a) single-cell 或 multi-region 真值錨；或 (b) 把目標從「subclone-specificity 驗證」改為「CROSS-HP 乾淨集的 PoC 展示」（43 區）。

**Reopen 條件**（scientific-rigor §8.3.1）：
- C1 新數據：COLO829 / 第二樣本 multi-region；
- C2 新方法：per-read CN-aware 分群解 multiplicity（解鎖 CN-gain 610 區的 CROSS-HP）；
- C3 新前置：single-cell 甲基真值可用。

**Spaced recall**：2026-07-28（30d）。

---

## §10.5 4-角度對抗驗證（workflow wf_251715fe，4 agents/362K token，全 high-confidence）

對「同-HP 條件化後 within-HP 甲基差異 = germline-free 的突變狀態差異」做 4 獨立角度查證：

| 角度 | verdict | 重點 |
|---|---|---|
| 邏輯/代數 | PARTIAL | 加法模型 `β=μ+g(HP)+c(somatic-cis)+t(lineage)+κ·CN+ε`：同-HP 下 germline 主效應 `g(h)−g(h)=0` **抵消（需真同-HP+相位正確+CN 相符）**；但 `c(somatic-cis)` **不可移除**且與 `t(lineage)` 混淆。**關鍵釐清**：`A−B` 是「癌 vs 癌」（兩 subclone 都帶 somatic 突變），**不是**「癌 vs 正常原始」（後者是 `A−N`）。matched-normal 對 `A−B` 在 SAME-HP **貢獻為零**（`(A−N)−(B−N)=A−B`），其唯一非冗餘角色 = **per-arm 極性錨 `A−N`/`B−N`**（定 somatic-vs-ancestral 方向）。|
| **實證 HP 一致性** | **SUPPORTS** | 🆕 **tumor 與 normal 的 HP1/HP2 指向同一條實體染色體**：① provenance — 兩 BAM 用**同一份 germline phased VCF**（md5 `c48e9509…` byte-identical）haplotag、tumor `--somaticMode` anchor 在 normal 單倍型；② 讀層 — 240 germline het SNV，tumor-HPk 與 normal-HPk 同 allele **98.85%（259/262）**、19 CONSISTENT/0 FLIPPED。→ 「normal 同-HP 基線可與 tumor-subclone 對齊」**結構前提成立**。|
| 殘餘 confound | PARTIAL | 嚴重度排序（針對「獨立 lineage 證據」）：①somatic-cis footprint（CRITICAL/結構）>②同分子共分離（HIGH/循環）>④CN/multiplicity（HIGH/經驗最大，CN-gain 53-69%）>⑤ASE/imprinting（多自動抵消）>③隨機雜訊>⑥覆蓋/解析。|
| 文獻/memory 一致性 | SUPPORTS | 與 verdict §7.3 + master spec §5 + baseline-dependence 鐵則 + Do 2020（somatic→cis 甲基真實但稀有）一致。修正：「germline-free」僅 SAME-HP 成立；「代表 somatic 軸」精確說法 = 「與 somatic-cis COMPATIBLE、單一 bulk 不可識別 footprint vs lineage」。|

**綜合**：用戶「同-HP 條件化 → 群內 germline-free」**成立（L2，4 角度一致）**；但「代表突變狀態分群差異」只在**關聯/共分離 + germline-clean** 意義成立，**不構成 subclone-lineage 的獨立證據**（殘差 = somatic-cis，單一 bulk 不可識別）。

## §11 延伸 pilot — 甲基能否當「拓樸排序」resolver（非分群、非 specificity）

> 動機：分群 ❌（循環）、subclone-specificity ❌（單樣本不可識別）後，唯一還可能的甲基用途 = 對**基因型已定的群**提供**排序/距離**（nested 深淺、sibling 平行）。此用法**非循環**（群已由基因型定），且 same-HP 下 germline-free。本 pilot 實證測試。

**設計**（`methyl_ordering_pilot.py`）：在「基因型已能定先後（A_determined）」+ SAME-HP + CN-clean 的區當 **ground truth**，測「甲基離 root(RR 祖先群)距離」是否隨基因型譜系深度單調增加（Spearman + permutation null）。🔑 **CpG 分 near(<1kb 任一 somatic)/distal(>1kb)**：只有 near 隨深度漲 = somatic-cis 痕跡累積（假時鐘）；distal 也漲 = 真譜系時鐘。

**結果**（326 候選 → 55 區可用 → 76 對深度-距離；271 區因無覆蓋足夠的 RR 祖先群被排除）：

| 測試 | all CpG | distal（真時鐘）| near（cis 痕跡）|
|---|--:|--:|--:|
| Pooled Spearman ρ | 0.182 | 0.177 | 0.164 |
| permutation p | 0.080 | 0.060 | 0.070 |
| nested 對「深群距離>淺群」 | 0.794 (27/34) | 0.636 (21/33) | 0.588 (20/34) |

**判讀（L3，弱訊號、未達顯著、功率不足）**：
1. **方向對但不顯著**：ρ≈0.18（深群傾向離 root 較遠），但 permutation p 全 > 0.05（0.06-0.08）。
2. **nested 比例 0.794 是假象**：nested 對不獨立（同區多對共用群）→ binomial 高估；正確口徑 = permutation p≈0.06-0.08。
3. **非純 cis 痕跡（小好消息）**：near(0.588) ≤ distal(0.636)，若純突變局部痕跡則 near 應最強 → 有一絲遠端訊號但太弱。
4. **功率嚴重不足**：A_determined 區多為淺樹（55 區幾乎全 2 層深度，僅 2 區 3 層）；深 nested 鏈（最需排序處）本就稀有 → 單樣本測不出單調梯度；per-region ρ 僅 1 區可算。

**裁決**：甲基**不能**當可信的拓樸 resolver（弱、未顯著、功率不足）；**可**當「軟提示」標於欠定區（須標 L3、附 ρ/p、非定論）。reopen 條件：COLO829/第二樣本（C1 新數據增功率）或 ISM read-level BERNOULLI 距離替代平均 Δβ（但訊號弱，上限有限）。

**架構備註（ISM 連接，源碼驗證）**：ISM 是**存在性偵測器**非拓樸排序器（per-region 對 HP/allele/tumor-normal 軸測甲基結構，read 無 genotype 標籤）。但 ISM 已有 `LabelTest.compute_group_distances`（輸出 6 條 HP-fine 群間距離 `HPFineD_*`，RegionProcessor.cpp:1098-1203），且 read-level BERNOULLI 距離（DistanceMatrix.cpp:254-302）比平均 Δβ 更有力。**若**外掛骨幹的 multi-SNV 共現分群給 read 貼 genotype 群標籤，ISM **能**產出 subclone 群間甲基距離原料——但本 pilot 顯示底層排序訊號本身太弱，換更強度量上限有限。⚠ ISM 內 `verification_class=="Subclone"` = `cluster_significant && !label_sig`（非監督 cluster ≠ HP/allele 軸），是 double-dip confounded 軸，**非基因型驗證的 subclone**，對外勿引。

## §12 甲基有界軟標記 — (b) 隱藏次結構旗標 + (a) ambiguous 軟提示

> 用途：分群❌/specificity❌/排序⚠️L3 後，把甲基當**有界軟標記**——(b) 標記「某基因型群內可能藏 ≥2 次群」、(a) 給欠定區排序軟提示。腳本 `methyl_auxiliary_annotation.py`（4 路平行）。

### (b) 隱藏次結構旗標（每群測甲基雙峰 + 扣 HP/CN，保守判準 GMM BIC + 峰差≥0.2 + 小組件≥4）

**漏斗**（`methyl_auxiliary_annotation.json` → b_hidden_substructure）：

| 階段 | 數 | 說明 |
|---|--:|---|
| 測試群數（≥12 reads + 甲基）| 3,416 | — |
| **unimodal（無次結構）** | 3,272（95.8%）| ✅ 信心說「無隱藏次結構」|
| bimodal | 144（4.2%）| 進去混淆 |
| ↳ cn-flagged（CN-gain multiplicity）| 96 | 排除（一群可能是拷貝）|
| ↳ hp-explained（雙峰沿 HP=germline）| 2 | 排除 |
| ↳ **residual-candidate（扣完仍雙峰）** | **46（1.35%）** | 全 CN-clean |

**residual 再嚴格過濾（critical sanity）**：
- **12/46 = 區域內在雙穩態**：同區 ≥2 群皆雙峰且**峰均值相近**（如 chr2:207617973 的 AR(HP1) 與 RR(HP2) 都 ~0.59/0.88）→ 與基因型/HP 皆無關 = 區域 epigenetic 雙穩，**非 subclone** → 排除。
- 剩 34 solo → 30 same-HP → **22 為 ALT 群**（只 ALT 群可能是 somatic subclone；RR 祖先群雙峰 definitionally 非 somatic）= **0.64% of tested**。

**🔴 誠實邊界（22 候選的真義）**：這 22 個 = 「基因型看不到、只有甲基看到的次分群」→ **正是已知 NEGATIVE 的『無監督甲基分群』路徑（double-dip）**：群內無第二個 sSNV 可佐證該甲基切分 → **無遺傳 corroboration** → 不能確認是 subclone，可能是 intrinsic epi-heterogeneity。故 = **L3 候選旗標（待 single-cell/orthogonal 驗證）**，非確認。

### (a) A_ambiguous 軟提示

76 個 A_ambiguous 區（缺中間群、nested vs sibling 待定）：**35 產生 L3 軟提示**（甲基離 root 距離排序）、41 無資料（缺 RR root 或 <2 群覆蓋）。⚠ 全域 ordering ρ≈0.18（§11）→ 這 35 個提示**低信心**，僅供人工複核起點，不可自動定案。

### 裁決：資訊是否有效可用？

| 用途 | 有效性 | 說明 |
|---|---|---|
| **負向篩選**（確認「無隱藏次結構」）| ✅ **可用** | 95.8% 群乾淨 unimodal，可信地排除 |
| **候選旗標**（標「可能藏次群」供跟進）| ✅ 有限可用（L3）| 22 區（0.64%）候選，但需 orthogonal 驗證、多數可能 intrinsic |
| **確認 subclone/群數** | ❌ **不可** | 無遺傳佐證 = double-dip 路徑；甲基自己定的群數不可信 |
| **ambiguous 排序** | ⚠️ L3 軟提示 | 低信心人工複核起點 |

**一句話**：甲基有界軟標記**能當保守的負向篩選 + 候選旗標（L3）**，但**不能確認**——與全研究一致：甲基 = bounded-auxiliary。22 個 ALT 候選 + 35 個 ambiguous 提示已落 `methyl_auxiliary_annotation.json`，供未來多樣本/單細胞優先檢視。

## §10 數字溯源表（§13-C）

| 數字 | 值 | 來源檔（grep-able）|
|---|---|---|
| pilot 區 / 有 CpG / CpG 對 | 8 / 6 / 2073 | `pilot_cis_control.json` → summary |
| tumor/normal Δβ p90 | 0.229 / 0.501 | 同上 |
| tumor 有訊號 CpG | 256 (12.3%) | `cis_control_scope_summary.json` → joint_distribution |
| 其中 germline-ASM / 候選 | 68(27%) / 138(54%) | 同上 |
| corr | −0.026 | `pilot_axis_alignment.json` |
| neutral pilot SAME/CROSS | 5 / 5 | `axis_alignment_neutral.json` |
| 全掃描分類數 | 1874 (skip 16) | `hp_alignment_fullscan.json` |
| CROSS / SAME / MIXED | 663(35.4%) / 1105(59.0%) / 106(5.7%) | 同上 + scope_summary |
| 交叉表 by_cn | gain 610/615/75; loh 10/467/27; neutral 42/17/4; loss 1/6/0 | `hp_alignment_fullscan.json` → by_cn |
| 最乾淨適用集 | 43 | `cis_control_scope_summary.json` → cleanest_applicable_set |
| needs_methyl 可分類 | 72（SAME54/CROSS12/MIXED6）| 同上 → needs_methyl_intersection |
| needs_methyl 乾淨可用 | 0 | 同上 → n_cross_hp_AND_neutral_clean |
| ordering pilot 候選/可用 | 326 / 55 | `methyl_ordering_pilot.json` → summary |
| ordering ρ (all/distal/near) | 0.182 / 0.177 / 0.164 | 同上 → pooled |
| ordering perm p (all/distal/near) | 0.080 / 0.060 / 0.070 | 同上 → pooled |
| nested 深>淺 (all/distal/near) | 0.794 / 0.636 / 0.588 | 同上 → nested_pair_frac |
| aux 測試群數 / unimodal | 3416 / 3272 | `methyl_auxiliary_annotation.json` → b_hidden_substructure |
| bimodal / residual-candidate | 144 / 46 | 同上 |
| residual 漏斗 (內在雙穩態/ALT候選) | 12 / 22 | 同上 hidden_candidates(過濾) |
| A_ambiguous 區/產生提示 | 76 / 35 | 同上 → a_ambiguous_hint |

**腳本**（皆可重跑複算）：`scripts/pilot_cis_control.py`、`analyze_pilot_distribution.py`、`axis_alignment_neutral.py`、`classify_hp_alignment_all.py`、`finalize_cis_control_scope.py`（皆於 `_assets/20260627_subclone_4axis_teaching/scripts/`）。
