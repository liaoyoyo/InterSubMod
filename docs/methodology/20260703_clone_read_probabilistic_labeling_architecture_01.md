<!--
建立時間: 2026-07-03
類型: 方法規格 — clone read 標記的機率式三層架構 + 甲基在遺傳欠定處為何不能定案（四配子拓撲 tie-break / 單位點非監督補結構）
狀態: 方法裁決（conceptual/methodology）；所有數字為引用既有 verified 文件之結果（非新分析）
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/20260701_ssnv_backbone_method_spec_and_correctness_audit_01.md, docs/methodology/20260701_asm_vs_clonal_definition_and_loh_methyl_verification_01.md, docs/methodology/20260701_topology_algorithm_audit_findings_01.md, docs/methodology/20260628_subclone_reconstruction_master_spec_01.md, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/four_gamete_aa_methyl_demo.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/sm_completeness_ledger.json
provenance: 概念/方法裁決文件。每個數字皆引用上列既有 verified 文件（§8 溯源表逐項對應），本文件不產生新分析數字。裁決由 2026-07-03 workflow（6 agents：文獻對照 + 循環對抗 + 逐特徵 + 自我一致，對抗 verify=sound_with_corrections）產出並落地。
-->

# Clone read 標記 — 機率式三層架構 + 甲基在欠定處為何不能定案

> **TL;DR**：clone read 標記可以做，且該做成 **Cardelino 式 beta-binomial 機率指派**（遺傳驅動、帶不確定性）；但**甲基（與 region 級 CN）絕不能進分群目標函數**——在單 bulk、cis-ASM 主導的資料上，用甲基去「排序四配子拓撲」或「非監督補單位點結構」都是循環，只能做**事後 L3 描述性註記 + 負向篩選（僅『不要 over-split』）**。

> **摘要（3 行）**：① **為何做這件事** — 使用者提案把 [somatic SNV + HP + SV + methylation + CN] 融成一個機率式 read 分群；核心風險是把已證實 NEGATIVE 的 cis-ASM 循環重新引入。② **怎麼做** — 三層角色分離：LAYER-1 目標函數只含非循環 genetic（sSNV 共現 + SV），LAYER-2 錨（HP-germline + CN），LAYER-3 甲基事後。③ **結果** — 四配子拓撲方向與單位點 subclone 結構在單 bulk **不可由甲基定案**（cis-ASM 對稱 + CNV multiplicity + 無真值 / tumor-only 無監督 NEGATIVE）；甲基的合法角色 = 描述性註記 + 保守負篩，正確加值來自遺傳側的 beta-binomial + partial-read 救援。

---

## §0 這份文件回答什麼

回答三個具體提問（後兩者為 2026-07-03 使用者精確提問）：

| # | 提問 | 一句話裁決 | tier |
|---|------|-----------|------|
| **基本** | 能否用 read 多特徵做 clone 標記？ | ✅ 能，用 Cardelino 式機率指派，**遺傳驅動、甲基排除在 likelihood 外** | L2 |
| **Q1** | 四配子（環）位置，能否用甲基判 `AR→AA` vs `RA→AA` 各自機率、解多環？ | ❌ 排序不行（循環）；✅ 標記/貼標籤/枚舉等機率候選樹可行且已做 | L2 |
| **Q2** | 單 sSNV 稀疏（read 只過一位點）位置，能否用甲基非監督結構補標、確認 ALT 內 CNV/subclone？ | ❌ 確認 subclone 不行（tumor-only NEGATIVE）；✅ L3 描述性註記 + HP-split + CN 偵測（非甲基）可行 | L2 |

**核心一句**：三題想用甲基做的動作都是「**在遺傳欠定處，用甲基把答案定出來**」——這是被 cis-ASM 循環擋死的那條線。甲基能加**旁邊的描述層**，不能加**定案力**。

---

## §1 三層架構（哪些進 likelihood、哪些是錨、哪些事後）

**不是「五特徵 concatenate 進一個 objective」，而是三層。**

```
LAYER 1  分群目標 (likelihood，非循環 genetic only)   = { somatic SNV 單分子共現(主),  SV breakpoint(共錨) }
LAYER 2  conditioning 錨 (固定 prior/offset，不當距離軸) = { HP-germline(1/2),  CN / local coverage }
LAYER 3  bounded-auxiliary 事後 (群凍結後才動)          = { methylation }
```

### 逐特徵角色表

| 特徵 | 循環性 | 角色 | 一句話 | tier |
|---|---|---|---|---|
| **somatic SNV 單分子共現** | 非循環 | **目標（主）** | 唯一遺傳真值軸；read×sSNV genotype(REF/ALT/OTHER) 共現 = likelihood 本體，非「五欄之一」 | L2 |
| **SV breakpoint support** | 非循環 | **目標（共錨）** | 唯一第二非循環 anchor（獨立於點突變與甲基）；⚠ 目前 pipeline 尚無 per-read SV support，屬 PROBE ⭐3 未實作 | L3 |
| **HP — germline (1/2)** | 非循環（allelic 軸） | **conditioning 錨（不進距離）** | germline haplotype family，做 cis-control / cis-trans；**非 clonal 分群維度** | L3 |
| **HP — somatic-integrated (1-1/2-1/3)** | 條件循環（衍生自 sSNV） | **剔除（雙重計數）** | 與 SNV 骨幹同源 → 進 objective = pseudo-replication；🔴 HP1-1 ≠ 已確認 subclone | L2 |
| **local coverage / CN state** | 條件循環（region-level） | **covariate/offset 錨（不進距離）** | 可做 VAF→CCF multiplicity 校正 + LOH flag；但 region 級**無 read 內判別力** → 當距離軸會讓 copy-number/mappability 主導分群（standalone REFUTED） | L2 |
| **methylation profile** | 循環（cis-ASM） | **bounded-auxiliary 事後** | 群凍結後才用：保守負篩 + L3-weak 旗標 + 上游 germline-HP phasing assist；**絕不進 likelihood、絕不 rank 樹** | L2 |

---

## §2 為什麼甲基不得進分群目標函數（cis-ASM 循環）

### DAG（scientific-rigor §4）

- 令 subclone label `L ≡ f(G_som)`（由 somatic SNV genotype 定義）。
- read 甲基 `M` 的父節點 = { `A`（germline allele → germline-ASM）, `G_som`（somatic-cis passenger 足跡）, `CN/LOH`（dosage unmask）}。
- 因此 `M` 是 `L` 之父 `G_som` 的**後裔**，也是 confounder `A`、`CN` 的後裔。

兩個獨立失效機制指向同一結論：
1. **循環/雙重計數**：把後裔 `M` 放進 `P(read|cluster)` 去推 `L` = 把 label 的下游後果餵回 label 推斷；`M` 攜帶的 lineage 資訊已在 `G_som` 內。
2. **collider/selection**：`M ← A` 且 `M ← G_som`，讓 `M` 影響 assignment = 條件於共同效果，開啟 `A ↔ G_som` 的 spurious 路徑（germline-ASM 軸洩漏進 somatic-lineage 指派）。

### 直接經驗證據 — 跨 CN 態 near-Δβ 幾乎相同（cis 主導）〔L2〕

`loh_genotype_methyl_verify.json`（2048 配對，near ±1kb / distal + 置換 + BH-FDR）：

| CN 態（RR_vs_1ALT） | n | Δβ_near 中位 | Δβ_distal 中位 | distal 顯著 (FDR<.05) |
|---|---:|---:|---:|---:|
| loh | 197 | **0.040** | 0.0254 | 6.1% |
| neutral | 49 | **0.040** | 0.0227 | 10.2% |
| gain | 1135 | **0.041** | 0.0247 | 13.0% |

- **near-Δβ ≈ 0.040 跨三種 CN 態幾乎相同** → genotype–甲基差**主要是突變的局部 somatic-cis 足跡**（`M ← G_som(cis)`），而 `G_som ≡` 定義 label。LOH 全集 near/distal = 1.46，去 cis（distal）後**僅 6–13% CpG 顯著** → lineage 殘餘極小。
- **甲基連 impute genotype 都做不到**：founder-imputation LOOCV balanced_acc 中位 **0.566 ≈ 隨機**；弱 ASM 區（<0.1，佔 **95%**, n=1513）**0.560 ≈ 隨機**；預測力由 ASM 強度驅動（corr=0.617）而非 lineage。

→ 任何非零權重把 `M` 放進 likelihood = **主動把 partition 拉向 cis-ASM/germline 軸**，且只在遺傳曖昧區（無真值可抓錯處）製造假確定性。

### 文獻防火牆（外部佐證）〔L3〕

**無任何已發表方法把甲基灌進 genetic-clone 指派 likelihood。** 三種防火牆：
- **anchor-to-independent-genetics**：clonealign 用 CN 定 clone、把第二模態（表現）當 dependent（CN→expression），從不當 clone-definer；Cardelino 樹取自獨立 bulk exome，scRNA 變異只做**指派**。
- **separate-axis / 改名**：Epiclomal 用純甲基分群但刻意叫「**epiclone**」，明說與 CN-defined clone concordant-**或**-discordant，從不宣稱甲基群 = genetic clone。
- **先除 confounder / 用 clock 非 function**：CAMDAC 先解 CN+purity 才詮釋腫瘤內甲基異質；MethylTree 用 stochastic **epimutation** 當 lineage clock 並 regress out cell-type methylation（非用 functional/allele-specific 甲基）。
- **read-level 甲基唯一既有角色** = germline-haplotype **phasing aid**（MethPhaser/NanoMethPhase），germline-only、從不跨進 somatic —— 對映本專案 LongPhase-S somatic 軸用零甲基。

---

## §3 Q1 — 四配子（環）位置能否用甲基判拓撲方向

**情境**：loci 1,2 四配子俱全（RR/RA/AR/AA），infinite-sites 違反 → sSNV 單獨無法定序。使用者提案：用甲基判「AA 離 AR 還是 RA 近」當機率，決定 `AR→AA` vs `RA→AA`。

### ✅ 可做且已做
- **清楚標記只用 sSNV 的環** — determinacy `incompatible = 118`（修後）+ upstream 環偵測。〔L1〕
- **給 sSNV 多位點關聯的 read 群貼標籤** — read×sSNV genotype 向量。〔L1〕
- **枚舉所有合理拓撲並列等機率** — `enumerate_candidate_trees`；單 bulk 多樹 posterior（Pairtree 式）不可行 → **維持等機率候選集**。〔L1〕

### ❌ 不能做：用甲基距離當機率排序拓撲（循環）〔L2〕

**理由 1 — cis-ASM 讓甲基距離「天生對稱」（最關鍵）。** 四配子下 AA 同時：
- 與 **AR 共享 locus-1 ALT** → locus-1 附近 cis-甲基 `AA≈AR`（cis-ASM 保證 = 重讀 genotype）
- 與 **RA 共享 locus-2 ALT** → locus-2 附近 cis-甲基 `AA≈RA`（同樣是重讀 genotype）

→ 局部 cis-甲基**在構造上對稱**，給不出方向。任何不對稱只能來自全基因組 lineage epimutation，而該殘差實測**去 cis 後僅 6–13% CpG 顯著、LOOCV 0.566≈隨機** → 不是拓撲機率，是 cis 足跡（對稱、無資訊）+ 不可用噪音。

**理由 2 — 四配子最常見成因是 CNV multiplicity，不是二次突變。** 稽核實測：被丟的 `independent`（4-gamete）對中 **CN-gain 佔大宗（774 對）**；C3 的 39 個高丟棄區 **22 個 CN-gain + 17 LOH**。此時「AA」是**同一分子上兩拷貝**（一帶 mut-1、一帶 mut-2）的 multiplicity 假象，根本不是 lineage 節點。而**同一 allele 多拷貝甲基完全相同** → 甲基原理上無法區分 multiplicity。（118 incompatible 是**非-CN-gain 的殘餘乾淨違反**，屬 LOH/mappability，CN-aware 建模也修不掉。）

**理由 3 — 等機率區無遺傳真值可校驗。** 甲基排序只能拿甲基自己驗 = 循環；且「甲基距離→拓撲機率」需 calibrated 的 methylation→lineage-time 時鐘（不存在且會被 confound）→ 連「機率」都未校準。〔L3〕

### 唯一理論非循環路（目前不可用）〔L3〕
**隨機漲落 CpG 時鐘（fCpG / stochastic epimutation）** 非 cis、不受 cis-ASM 汙染，single-cell 有成功案例（MethylTree）。但需 **single-cell 解析度 + 校準時鐘**；bulk read-level ONT 沒有，且本資料殘差訊號實測太弱。→ **reopen 條件 C3（single-cell）**，非現況能力。

### Q1 正解
環的**標記 / read 群標籤 / 枚舉等機率候選樹**可做且已做；**用甲基排序拓撲、給各拓撲機率**不行。區分「真 homoplasy vs CNV multiplicity」是 **CN 資料**的事（不是甲基）；定不出方向的區 = **「定不出來即答案」**（保留高熵，不用甲基抹平）。

> ⚠ **§9 實證修正**：本節原「cis 讓 near-band 天生對稱」的**機制**經 2026-07-03 demo（§9）**未被乾淨證實**（near-band 甲基噪音、非完美對稱，因 somatic-cis ASM 多數位點弱）。但「甲基無法判四配子拓撲方向」的**結論經更直接路徑成立** —— 見 §9。

---

## §4 Q2 — 單 sSNV 稀疏位置能否用甲基非監督補結構

**情境**：read 只覆蓋單一 somatic locus（無共現可用），現況無法 linkage。使用者提案：用 ISM 甲基非監督結構數 / 標籤補標，確認 ALT 內是否有幾群、是否有 CNV/subclone。

### ❌ 不能做：甲基非監督分群「確認 ALT 內有 subclone」（循環）〔L2〕
單位點**沒有第二個遺傳 marker 在同一批 read 上** → ALT reads 間甲基差異**無法歸因**：可能是 (a) germline cis-ASM、(b) LOH-unmask、(c) CNV dosage、(d) 真 subclone lineage，四者在單位點**不可拆**。這正是已判 **NEGATIVE 的 tumor-only 無監督甲基分群（double-dip）**：先用甲基挑最分離 partition 再宣稱 subclone = 循環。ISM 輸出的「幾群 / 幾種標籤」是 **epi-cluster 數，不是 subclone 數**（被 germline-ASM + CNV dosage confound 主導）。

### ✅ 可做且有真加值
1. **HP tag 把單位點 ALT reads 依 haplotype 拆開**（多數 read 已 HP-tagged）—— 遺傳、非循環；拆後可 characterize 每 haplotype 甲基，但仍只到「有無 ASM」，**不升級成 confirmed subclone**。〔L2〕
2. **CNV 用 coverage/CN 偵測，不是用甲基** —— 「ALT 內有無 CNV」是 read-depth/allele-ratio 的事，orthogonal 非循環；甲基是被 CN confound 的下游，不能當 CNV 偵測器。〔L2〕
3. **甲基當 L3-weak annotation** —— 對單位點標「是否 bimodal / 有無 ASM」，**分開軸的描述性註記**（此層確實比現況多資訊），但是 characterization 非 confirmation。〔L3〕
4. **label-first 一致性檢定（state③）** —— 固定距離矩陣、置換**外部** HP/germline-het 標籤、群內分層 permutation（非循環結構）；測出「真關聯但主軸 = germline-haplotype **非 subclone**」。〔L2〕

### Q2 正解
甲基**可替單位點加一層描述性註記**（有無 ASM/bimodal、對照 HP），比現況多資訊；但**不能**把單位點「確認成有 subclone 結構」。想在單位點加**遺傳**資訊，正解 = (a) HP-split、(b) CN/coverage 偵 CNV、(c) partial-read soft-likelihood 救回共現（見 §5），三者皆非循環，甲基只做旁邊註記欄。

---

## §5 機率框架 — 真加值 vs 只是換記法

### ✅ 真加值（全在 genetic 通道，甲基零介入）
1. **beta-binomial / ONT-error-aware read-count likelihood** 取代硬 REF/ALT/OTHER → 直補稽核兩洞：C3（貪婪去噪無 eps floor，39 區 drop>5%、max 0.609 仍標 A_determined）+ C2（4-gamete 在 build_tree 靜默丟）。硬閾值 conflict → **機率 conflict score**。〔L4 前瞻〕
2. **per-read soft likelihood 吃 partial-coverage reads** —— read-level 框架的**真獨有優勢**：現況 complete-case genotyping 丟掉 **40.4%（2887/7143）空向量區**；一條覆蓋 sSNV 1&3 未覆蓋 2 的 read 仍供真實 1–3 共現證據 → **真多出 identifiability**（非換記法）。⚠ 可救比例須先算（`pairs_eps2 ∩ 空區`），**不可把 40.4% 當可救目標**。〔L3〕
3. **不確定性量化本身 = 交付成果**：determinacy 從 4 個離散桶 → per-region posterior entropy / CCF credible interval。〔L3〕

### 🔁 只是換記法（不加 identifiability）
- **樹層軟指派已跑完**：候選樹**維持等機率**；likelihood 對候選樹平 → posterior 必吐 P=1/N。連續 posterior 在欠定區**只是換記法**，變不出 reads 裡不存在的資訊。〔L2〕
- ⚠ **最大陷阱**：若模型預設輸出 MAP 單一樹（PyClone/PhyloWGS 風格），會把平的 posterior 藏成**假點估計** = 撞我們最想避免的 overclaim。**高熵 posterior 本身就是「定不出來」的量化陳述，要保留而非抹除。**

---

## §6 最小可行 pilot（不違反循環）

**Scope（pilot 級，對稽核發現做 falsification）**：只在 HCC1395 稽核標 HIGH/MED 那批（C3 高丟棄區 + C2 內藏 4-gamete + incompatible）；不全基因組、不先跨樣本。

**只換 likelihood 層，重用既有 topology 骨架**：
1. beta-binomial per-sSNV read-count model（ONT per-base error 從 BAM QS）→ soft genotype posterior 取代硬 REF/ALT/OTHER
2. 4-gamete/perfect-phylogeny → 機率 conflict score（取代 C3 eps=2% 硬 floor + C2 硬丟棄）
3. CN-aware allele-count 項處理 gain-multiplicity；SV breakpoint 作非循環共現 anchor（若 per-read SV 可得）
4. 餵進 `enumerate_candidate_trees` 既有 softmax → per-region 候選樹 posterior + CCF credible interval

**🔴 硬約束（缺一即退回循環）**：
- (i) likelihood **只含 genetic**（sSNV counts + CN + SV）；**甲基進入零項**，只在群凍結後做 downstream L3-weak annotation。
- (ii) HP = **固定 conditioning 軸**（germline 兩棵 HP 樹分開），非 mixture component。
- (iii) **不做全基因組 mixture**；跨區整合明標 non-identifiable，把「定不出來」量化成 posterior entropy **保留而非抹除**。

**甲基事後（群凍結後）只能**：
- (a) **保守負篩** —— 僅「甲基 unimodal → **不要 over-split**」（從不移除/降級遺傳 read）。🔴 **絕不**「甲基矛盾時降信心/剔除遺傳 read」（那是把 label 的 cis 下游餵回指派 = 循環）。
- (b) 加**分開標示**的 annotation 欄（read-only，不改遺傳群信心/成員）。
- (c) 上游 **germline-HP phasing assist**（作用在 somatic-label 之上游的 germline-haplotype 層、phasing 本體用零甲基，正交合法）。
- **絕不**：確認群 / 創群 / rank 拓樸 / impute genotype / 破 equiprobable tie。

**可驗證成功判準（對稽核真值 falsify，非自證）**：conflict score 能否把「A_determined 內藏 4-gamete」+「丟>5% reads 仍標 A_det」正確重分類為 incompatible/降信心，且把 CN-gain 假 4-gamete 溶回可定樹。能 → 真修補；不能、仍判平 → 證實那些區本就欠定，**高熵 posterior 就是最終答案**。**兩種結果都 publishable、都不觸循環。**

**流程紀律**：`/methodology-audit` → Python pilot（動 C++ 呼叫層才 `/cpp-change`）；§13.0 產數字先落 `.json` → Read 驗證非 error → 才寫報告。

---

## §7 誠實邊界 + evidence tier + reopen 條件

- **本裁決整體 tier = L2**：循環的邏輯核心靠兩條腿 —— (a) 單 bulk 下「somatic-cis passenger 與 lineage epimutation **不可區分**」（最硬，近演繹）+ (b) 跨-CN near-Δβ 近乎相同（0.040/0.040/0.041）+ LOOCV≈隨機的**經驗**支持（L2）。**不是「幅度無關的演繹鐵律」**，而是「單 bulk 非可識別 + 經驗上 cis 主導」。
- **未達 L1** 之因：confound 幅度仍以 **HCC1395 單細胞株**為主，缺 single-cell / multi-region 正交真值錨。**此不確定性不影響「甲基不得進目標函數」的結論方向**，只影響「未來是否存在可用的殘餘 lineage 訊號」。
- **reopen 條件**：**C1 COLO829 多樣本 / C3 single-cell**。若正交真值顯示存在「非-allelic ∩ 非-dosage ∩ excess-over-null」的 lineage epimutation，可重審甲基「僅確認、capped 權重」形式 —— 但仍須先過 label-first exchangeability（固定距離矩陣、置換外部標籤、群內分層 permutation），且僅限「確認已有遺傳支持的 split」，**永不破 equiprobable tie**。

---

## §8 Provenance（數字 → 來源，§13-A/§8.4）

| 數字 | 值 | 來源檔案:定位 |
|---|---|---|
| determinacy（修後）A_determined / A_ambiguous / B_pairwise / incompatible | 1741 / 62 / 943 / 118 | `20260701_ssnv_backbone_method_spec_and_correctness_audit_01.md:19` |
| 修前診斷數（勿用）1812 / 12 / 306 | — | 同上:19,108,154 |
| C3 高丟棄區 | 39 區 drop>5%、max 0.609、22 CN-gain + 17 LOH | 同上:107 |
| C2 靜默丟 independent 對 / CN-gain | 5198 對、1137 區(15%)、774 CN-gain | 同上:108 |
| 跨-CN near-Δβ（loh/neutral/gain） | 0.040 / 0.040 / 0.041 | `20260701_asm_vs_clonal_definition_and_loh_methyl_verification_01.md:117,119,121` |
| LOH near/distal 比 | 1.46（僅 6–13% distal 顯著） | 同上:124,26 |
| founder-imputation LOOCV balanced_acc | 0.566/0.560≈隨機；弱ASM 95%(n=1513) 0.560；corr=0.617 | 同上:89,90,91 |
| 空向量丟棄 | 40.4%（2887/7143） | `20260701_topology_algorithm_audit_findings_01.md:25` |
| 候選樹策略 | 單 bulk 多樹 posterior 不可行 → 維持等機率；絕不用甲基排序 | 同上:35；`project_candidate_tree_ranking_impossible` |
| 甲基 bounded-auxiliary / germline-HP phasing assist | 4 用途窮盡；phasing assist 屬 germline-haplotype 層 | `project_methylation_use_exhausted_bounded_auxiliary`, `project_methyl_phasing_assist_line` |

**關聯 memory**：`project_clone_read_labeling_probabilistic_architecture`（本裁決 hook）、`project_methylation_use_exhausted_bounded_auxiliary`、`project_candidate_tree_ranking_impossible`、`reference_hp_tag_definition_and_subclone_caveat`、`project_ssnv_backbone_o2_and_correctness_fixes`、`project_tumor_only_axis_negative_subclone_classification`。

---

## §9 Q1 實證 demo — 四配子甲基能否判拓撲方向（2026-07-03，HCC1395）

> **partial flag**：single-sample **HCC1395** / `incompatible`（四配子）子集 / 小 n（20 pairs）—— 非全樣本、非全基因組（demo/falsification scope）。
> script `_assets/20260627_subclone_4axis_teaching/scripts/four_gamete_aa_methyl_demo.py` → `four_gamete_aa_methyl_demo.json`；5:22 wall / 234MB RSS；seed 20260703。

### §9.1 設計
118 個 incompatible 區 → 找齊 RR/RA/AR/AA 四配子的 sSNV pair（各群 ≥4 reads）→ 每 read 甲基拆 **near-i（±1kb of locus i）/ near-j / distal（>1kb of both）** 三帶；**same-HP 閘**（AR/RA/AA 同 germline HP，≥0.6 dominant）+ **CN 分層** + **distal 置換 null**（2000×）。

### §9.2 Funnel（20 對四配子 pair）
| 關卡 | 數 | 意義 |
|---|---|---|
| 四配子 pair | 20 | RR/RA/AR/AA 俱全 |
| same-HP 前提成立 | 12 | 可問「AR→AA vs RA→AA」 |
| **cross/mixed-HP（前提不成立）** | **8** | 兩突變在不同單體型 → **AA = doublet/CN artifact 非譜系**，問題本身不成立（40%） |
| same-HP 且 CN-clean(neutral/loh) | 6 | 全 loh |

### §9.3 🔑 關鍵發現（比原理論預測更直接）
1. **唯二 distal 顯著（p<0.05）的 pair 全在 cross-HP CN-gain 區**（chr9:41777788：S=−0.168 **p=0.0015** / S=0.119 **p=0.0025**）→ 甲基「顯著」時抓的是 **cross-HP / CN-gain doublet 的 confound**（兩條不同單體型甲基不同），**不是 lineage**。〔L2〕
2. **6 個 well-posed clean 案例（same-HP + loh）distal 置換 p 全 > 0.05**（min 0.066）→ 問題問得成立時，甲基距離**從不顯著**。〔L2〕
3. **parent-call 不一致**：clean 案例 near-based vs distal-based「較近母系」多數相反（~4/6）、與 VAF 也常相反 → **無穩定方向**。〔L2〕
4. **cis 標籤驗證 = 混雜**：near-band 甲基距離小（~0.08–0.12）且**不乾淨 track 共享 allele**（near_i 共享(AR) 0.105 vs 非共享(RA) 0.076；cis 僅 ~2–3/6 案例成立）→ 因 somatic-cis ASM 在多數位點弱（呼應「僅 ~5% 強 ASM 區」），near-band 是**噪音而非乾淨 cis readout**。〔L2〕

### §9.4 誠實修正（§13.7）
原 §3「cis 讓 near-band 天生對稱」的**機制未被乾淨證實**（near-band 噪音、非完美對稱）。但「甲基無法判四配子拓撲方向」的**結論經更直接路徑成立**：(a) 40% pair 問題不成立（cross-HP）；(b) well-posed 乾淨案例 distal **從不顯著**、方向不一致；(c) **唯二「顯著」是 cross-HP CN-gain confound**。→ **Q1 用甲基排序四配子拓撲 = 經驗否證。** 正解仍是「定不出來即答案」+ **CN 資料**區分 multiplicity vs homoplasy（demo 中 chr9/chr22/chr14/chr2 大 AA + 微 AR/RA + 多 cross-HP = 典型 CN-multiplicity artifact）。

### §9.5 §6 pilot 確認（demo 落地驗證）
- **pilot target 118 incompatible 確認**（CN：loh 77 / gain 31 / neutral 9 / loss 1）。
- **demo 驗證 §6 硬約束（甲基 zero-in-likelihood）正確**：最難的四配子區甲基也加不了可靠解析 → pilot **只用 genetic**（beta-binomial + CN-aware conflict score），不碰甲基。
- **§6 falsification target 具體化**：CN-gain 四配子 pair（chr9/chr22/chr14/chr2：大 AA + 微 AR/RA + cross-HP）= multiplicity artifact = §6 CN-aware conflict score 應重分類的對象。
- **Q2 單位點 scope 量化**（`sm_completeness_ledger.json`，sSNV 層）：isolated_singleton **8320** + underpowered **5458**（= 單位點稀疏）vs linked **21554**（有共現）。
- ⚠ **§6 開放項（誠實邊界）**：partial-read 可救比例（40.4% 空向量 ∩ pairs_eps2）**仍未算** → **不可把 40.4% 當可救目標**（沿用 audit doc 邊界，需另跑 `pairs_eps2 ∩ 空區`）。

### §9.6 Provenance（§9 數字 → 來源）
| 數字 | 值 | 來源 |
|---|---|---|
| funnel / distal 顯著 / cis 帶 | 20 pair、same-HP 12、cross-HP 8、clean 6、distal-sig 2(皆 chr9 cross-HP gain)、near_i 0.105/0.076 | `four_gamete_aa_methyl_demo.json`（本輪 §13.0：先 compute→Read 驗→才寫） |
| chr9 顯著 pair | S=−0.168 p=0.0015 / S=0.119 p=0.0025 | 同上 `cases[]` |
| Q2 單位點 scope | isolated 8320 / underpowered 5458 / linked 21554 | `sm_completeness_ledger.json:buckets` |
| pilot target CN 分層 | loh 77 / gain 31 / neutral 9 / loss 1 | `topology_per_region.json`（det where determinacy==incompatible） |
