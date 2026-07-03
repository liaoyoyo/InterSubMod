---
title: CN 軸 per-read 可行性 + SAVANA 能力邊界 + sSNV/甲基角色分工（方法定位）
date: 2026-07-03
status: in_progress
task_type: positioning-note（方法定位，非新實驗）
build_branch: research/subclonal-reconstruction-202606
tier_of_claims: 概念裁決=L1（推理）/ SAVANA 能力=primary-source 對抗驗證，投稿前 L3 待 /citation-verification / 內部數字=既有 ⭐3 prior work（已 grep 回真值）
data_sources: >
  docs/experiments/in_progress/2026/06/20260621_kism_vs_cn_perread/20260621_kism_vs_cn_perread_result_01.md,
  docs/plans/20260620_ont_cnv_sv_subclone_verification_feasibility_01.md,
  docs/methodology/20260629_clp_cn_crossval_integration_01.md,
  SAVANA Nat Methods 2025 (DOI 10.1038/s41592-025-02708-0; PMC12240814; GitHub README v1.3.8 2026-06-30)
related:
  - docs/method_comparison/20260630_ism_positioning_vs_prior_work_01.md
  - docs/methodology/20260621_baseline_discipline_conceptual_definitions_01.md
---

# CN 軸 per-read 可行性 + SAVANA 能力邊界 + sSNV/甲基角色分工

> **一句話定位**：CNV 軟體能告訴你「**某區域幾倍**」（segment-level），但**無法**把倍數落到單條 read、無法指認「哪些 read 是被複製的 copy」、無法用 CN 把 read 分群。read-level 判別靠的是 **somatic sSNV 共現骨幹（＝演化標籤唯一來源，非循環）**；**甲基是有界的註記/資訊層（characterize，不 confirm、不產標籤）**。單樣本天花板 ⭐3。
>
> **這份文件回答**：(1) CN 能不能做到 per-read？(2) SAVANA 精確能到哪、到不了哪？(3) 那 sSNV + 甲基 要怎麼「輔助分群、標演化標籤、提供資訊」才不越界？

---

## §0 起因與範圍

沿論文 slide 3-4 方法圖審查延伸的一串釐清：能不能用 CNV 工具（尤其 **SAVANA**）標記單條 read、或估完 CNV 再把 reads 分群？延伸到「sSNV + 甲基 read-距離/群聚」能否當輔助分群、掛演化標籤、提供資訊給後續分析。本文把結論釘死，作為投稿時回應 reviewer「你怎麼不用 CN 分 read / 甲基到底做什麼」的口徑依據。**scope＝方法定位（概念 + 工具邊界 + 角色分工）**，非新實驗；引用的內部數字均為既有 ⭐3 work 且已 grep 回真值。

---

## §1 概念裁決 — CN 本質不是 per-read

Copy number 是從**聚合訊號**（binned depth log-R + het-SNP 的 BAF）+ purity/ploidy 擬合出的**片段（segment）屬性**。單一分子告訴不了你「這個位點有幾個 copy」——那是族群計數量。把聚合拿掉，訊號就消失。

| 問題 | 能不能知道？ |
|---|---|
| 「**這個區域**是幾倍」（segment-level 整數 CN，total + minor allele） | ✅ **能**（SAVANA/Wakhan 的正職輸出） |
| 把那個倍數**落到「哪些 read」**（per-read CN 標記） | ❌ 不能 |
| 指認「**哪些 read 是多出來的那份 copy**」 | ❌ 不能 |
| 用 CN 把 reads **分群**（read→subclone 指派） | ❌ 不能 |

**⚠ 常見誤解校正**：「無法知道此區域有幾倍」是錯的——**區域幾倍是知道的**（這正是 CN caller 的核心輸出）。不能知道的是把倍數**落到單 read / 單 copy**。

**為何「哪些 read 是被複製的」原理上不可能**：一段純 CN gain（如 3 copies）裡，多出來的 copy 是**序列完全相同**的——深度變高，但**沒有任何 per-molecule 標記**能區分「這條 read 來自原本那份 vs 多出來那份」。深度只給你「更多 reads」，給不了「哪條屬哪份 copy」。

> **唯一能把同區 copy 分開的分子途徑，是它們攜帶不同的 somatic mutation**——而那一刻，做分辨工作的是 **somatic 共現/sSNV 骨幹**，CN 只是 context。這反而**強化 ISM 設計**：read-level 判別器是共現骨幹，不是 CN。

---

## §2 兩個唯一的 per-read「CN-linked」hook

read 能攜帶 CN 相關訊號只有兩條路（其餘皆不行）：

| Hook | 什麼標記了 read | 對甲基 β 獨立？ |
|---|---|---|
| **(a) SV/CNA 斷點** — split/supplementary alignment 實際跨越 subclonal junction 的 read | junction-supporting read-ID（SAVANA/Severus `*_read_support.tsv`；來自 CIGAR/supplementary，與 β 正交） | 🟢 **真正獨立** — read 有沒有跨斷點與其 CpG β 無因果耦合（feasibility §「唯一非循環 anchor / GO」） |
| **(b) HP-phased read 在 AI/LOH 區** — read 靠 HP tag 繼承該段 allele-CN | HP tag → segment allele-CN（Wakhan HP1/HP2） | 🔴 **不獨立＝循環** — CN=1（LOH）本身決定哪個 allele 存活（LOH-unmask，W2/I2）；「甲基群對齊此軸」可由單一 CN 事件人為製造 |

**⚠ HP 標籤的歸屬**：haplotype 標籤來自**上游 phaser**（LongPhase/haplotag，靠 germline SNP），**不是** CN caller 產的。SAVANA 只是**清點** SV-supporting reads 各屬哪個 HP（VCF INFO `TUMOUR_ALT_HP`/`NORMAL_ALT_HP` 計數），它不製造 HP tag、也不把 CN phase 到 HP。

---

## §3 SAVANA 精確能力邊界（本 session 對抗驗證）

一手源：Nat Methods 2025（DOI 10.1038/s41592-025-02708-0；PMC12240814；bioRxiv 2024.07.25.604944）+ GitHub README v1.3.8（2026-06-30，Apache-2.0，Cortés-Ciriano lab）。經 4-agent workflow 的對抗驗證 phase 逐條 REFUTE/CONFIRM。

**SAVANA 有的**：somatic SV（RF 分類器 → `.classified.somatic.vcf`/`.bedpe`）、**allele-specific CN（total + minor，segment 層）**（`*_segmented_absolute_copy_number.tsv`）、purity/ploidy（單一 grid-search fit）、ONT + PacBio、tumour-normal + tumour-only（`--g1000_vcf`）、**per-read SV 斷點支持**（`*_sv_breakpoints_read_support.tsv`，逐 junction 列支持 read-ID）。

**SAVANA 沒有的**（對抗驗證 REFUTED）：
- ❌ **per-read CN 標記**（CN 是 depth-based bin/segment 層，任何 read 上無 CN 欄位）。
- ❌ **用 CN 分群 reads**（README/paper 無 read→subpopulation 分群）。
- ❌ **subclonal CN**（只有單一全基因組 purity/ploidy fit，無多-clone 分解）。
- ❌ **genome-wide HP1-vs-HP2 per-segment CN 標籤**（是 major/minor-allele 意義的 allele-specific，非 HP1/HP2 標籤；且 SV 在 reads 有 phase 時達 single-haplotype resolution，但那是 SV 不是 CN）。

**Bottom line**：SAVANA 不能標單 read 的 CN（NO）、不能用 CN 分群 reads（NO）。唯一 per-read 輸出是 SV 斷點 read-ID 支持 + SV-supporting reads 的 HP 計數。

### 工具全景（bulk 長讀）

| 工具 | HP/allele-specific CN | subclonal CN | 分群 reads？ | tier |
|---|---|---|---|---|
| **SAVANA** | allele total+minor，無 HP1/HP2 split | ❌ segment 層 | 🔴 否 | A（Nat Methods 2025） |
| **Wakhan**（Kolmogorov） | ✅ **唯一** HP1/HP2 分離（`*_segments_HP_1/2.bed`） | ✅ `--copynumbers-subclonal-enable`→分數 CN+CCF，**仍 segment 層** | 🔴 否（per-read 無） | B（preprint；自建 self-benchmark，非獨立第三方） |
| **Severus** | SV 工具、haplotype-specific SV + phased 斷點圖；`--output-read-ids` 給 per-read junction 支持 | SV 層 | 🟡 只能靠 junction-read（hook a） | A（Nat Biotech 2025） |
| Spectre / HiFiCNV | 深度型、germline 導向 | 🔴 | 🔴 | C（somatic 用途 L3） |
| ClairS-TO | 不產 CN 圖；verdict 模組**消費** CNA+purity 標變異 germline/somatic/subclonal-somatic | 標變異非建 CN | 🔴 | B（內部 verdict = filter NEG） |

> **單細胞對比**：Ginkgo/SCOPE/**CHISEL** 能 per-**cell** 全基因組 CN 並分群 cells——但那是 per-**細胞**（每細胞＝自己一個 library，聚合發生在單細胞基因組內），**不是 bulk 的 per-read**。這是單細胞/多區域仍是確認金標準、單 bulk 封頂 ⭐3 的結構性原因。

---

## §4 接回 ISM 既有裁決（勿重推）

（源：`20260620_ont_cnv_sv_subclone_verification_feasibility_01.md`、`20260621_kism_vs_cn_perread_result_01.md`、`20260629_clp_cn_crossval_integration_01.md`；已 grep 回真值。）

1. **「CN 分群 ↔ 甲基群 對齊 → 驗 subclone」= 循環，standalone REFUTED**：CN 軸與甲基軸經 LOH-unmask-ASM 隱性耦合，「對齊」可由單一 CN 事件人為製造 → 對齊 ≠ subclone ≠ 獨立確認。
2. **k_ISM ⊥ k_CN（實測 ⭐3）**：Spearman **ρ=−0.038**（p=7.7e-11，僅因 n 大而顯著，效應量≈0）；345,714 read 指派，median **ARI(cluster,hp)=0.0**；optimal_k vs read 數 ρ=−0.047、vs SEQC2 CN ρ=−0.031(NS)。→ 甲基分群**不是 CN 鏡子**（有獨立內容），但存在的結構**大多不對齊任何遺傳軸**（bulk 下無法形式判 subclone）。
3. **SV 斷點軸 = 唯一非循環 anchor（GO），但實測稀疏弱**（SV 對齊 **5.1%** vs HP 4.5%）→ 原料在（per-read SV 支持 TSV），但**尚未 operationalize 成 production read-level anchor**。
4. **per-read CN 指派：全領域無、ISM pipeline 也無**（只有一次 one-off 分析用 region-overlap 假指派，ARI≈0）。
5. **CN 在 ISM 的正確用途 = 固定的 segment 層參考層**（算 k_ISM 相對 k_CN 的 excess、LOH masking、allele context 做 characterization/confound 排除），**非** per-read 分群器。
6. **SAVANA-as-fixed-reference 缺口**：6 樣本中 **5/6 無等價 WGS CN 參考**（CLP fill rate 0%/0%/7.4%，HCC1395 的 100% 來自 SEQC2 非 CLP）→ `20260629_clp_cn_crossval_integration_01.md` 明講「需 SAVANA」；`sm_region_integration.py` 已加 `SM_CN_SPARSE` + CN BED 介面留好接口，**尚未接進 production topology pipeline**。
   - ⚠ **LOH Jaccard（SAVANA vs SEQC2）口徑歧異未定**：內部記 0.847（C13 卡）vs 0.928（memory）——region-level vs bp-level 基準差異，**引用前必釘**，本文不寫死單一數字。

---

## §5 sSNV / 甲基 角色分工（核心）

### 5.1 兩者地位不同，不可並列為「輔助分群」

| 資訊 | 角色 | 循環？ |
|---|---|---|
| **somatic sSNV 共現** | 🟢 **主骨幹 = 演化標籤的唯一來源**（不是輔助） | 非循環（HP tag 來自零甲基 phasing） |
| **甲基 read-距離/群聚** | 🟡 **輔助 / 註記層**（描述、負篩、L3 flag） | 用來「確認/分群/貼標籤」就循環 |

sSNV 共現不是「輔助確認」——它**就是**判別器與演化標籤的發源地；甲基才是輔助。

### 5.2 方向決定合法性（label-first vs unsupervised-then-label）

- ✅ **label-first（合法，非循環）**：先用 **sSNV 共現 + HP** 給每群 read 演化標籤（HP-family + nested/sibling，four-gamete）→ **再**問「甲基群聚**是否對齊**這個既有標籤」（PERMANOVA HP-label 檢定，非循環，因 HP 來自零甲基 phasing）→ 對齊程度當 **characterization** 報告。
- ❌ **unsupervised-then-label（越界）**：先用**甲基群聚**分群 → 再貼演化標籤 → 宣稱 subclone。這是已證實的 **tumor-only double-dip NEGATIVE**，也是 slide 3「用甲基確認分群」那句的本質。

### 5.3 三條護欄

1. **標籤來源護欄**：演化標籤**只能**來自遺傳骨幹（sSNV 共現 / HP），**永不**來自甲基。甲基只能**掛在**已被遺傳標記的 read 上當註記。
2. **方向護欄**：label-first（遺傳標籤 → 查甲基對齊），非 cluster-first-then-label。
3. **動詞護欄**：甲基「**描述 / corroborate / L3-flag**」，不寫「**confirm / discriminate subclone**」。甲基相似 ≠ 譜系相似（可能是 cis-ASM/CN；見 baseline-discipline 詞典）。

### 5.4「輔助之後的分析」——合法用途 vs 禁止

- ✅ **負篩**：甲基說「這群 read 甲基同質」→ 降優先。
- ✅ **L3 假設 flag**：「遺傳群 X 內部**還有**甲基次結構」→ 標**候選**，建議 single-cell/multi-region 跟進（生假設，不下結論）。
- ✅ **characterization 註記**：把每個遺傳定義的 lineage 掛上其甲基 β-profile 做描述。
- ❌ **禁**：「甲基群 = subclone」「甲基確認了譜系」「演化標籤由甲基產生」。

---

## §6 投稿口徑（reviewer 防守）

- 「你怎麼不用 CN 把 read 分群？」→ **CN 是 segment 屬性，per-read CN / read→subclone 全領域無工具**（SAVANA/Wakhan 全 segment 層）；同區 copy 序列相同無 per-molecule 標記，唯一分辨靠不同 somatic mutation＝我方共現骨幹。CN 在本法定位為固定參考層（excess-over-CN、LOH mask、confound 排除）。
- 「甲基到底做什麼？」→ **bounded-auxiliary 資訊/註記層**：label-first 下查甲基是否對齊遺傳標籤（characterize），可負篩、可生 L3 假設，**不 confirm subclone、不產演化標籤**。天花板 ⭐3；升確認需 single-cell/multi-region。
- 🔴 **不可宣稱**：「甲基驗證 subclone」「用 CN 分 read」「演化標籤來自甲基」。

---

## §7 待辦 / citation-verification queue（投稿前）

1. SAVANA Nat Methods 2025（DOI 10.1038/s41592-025-02708-0 / PMC12240814）逐字核（能力邊界句）。
2. LOH Jaccard（SAVANA vs SEQC2）口徑歧異 0.847 vs 0.928 釘死（region- vs bp-level）。
3. Wakhan「subclonal CN」為自建 CASTLE self-benchmark，非獨立第三方——引用時標明。
4. SAVANA 接進 production topology pipeline 補 5/6 樣本 WGS CN 缺口（`sm_region_integration.py` `SM_CN_SPARSE` 介面）＝工程待辦，非本文 claim。

---

## §8 接進 slide 3-4 敘述（paste-ready）

> 把 §1-§7 定位落到論文方法圖 slide 3（共現骨幹→區域局部拓撲）與 slide 4（甲基結構+雙軌驗證）的講稿與 caption。直接可貼進 PPT speaker notes / 圖說。

### §8.1 一句話串接 slide 3 → 4

> **slide 3 的 somatic 共現骨幹「給出並命名」譜系（＝演化標籤）；slide 4 的甲基與 CN「描述並註記」這些已命名的譜系（characterize）。分工：骨幹判別、甲基/CN 輔助；天花板 ⭐3。**

### §8.2 Slide 3 敘述 — somatic 共現骨幹 → 區域局部拓撲

**主敘述（speaker note）**：我們用 ONT 單分子上的 **somatic sSNV 共現（four-gamete / census）當非循環骨幹**重建區域內譜系拓撲——這是**演化標籤的唯一來源**（非循環，因 HP tag 來自零甲基 phasing）。

**三態講法**：
- **linear（nested，缺一 single-ALT 配子）**：直系；從 RR 最小突變（parsimony 生根）。
- **branched（mutual-excl，缺 AA）**：姊妹；依**子群比例（VAF / read-count）**顯示（heuristic，**非時序**）。
- **incompatible（four-gamete 全滿）**：**不強建單樹**；標 likely **CN-multiplicity / mapping artifact + FP**，加 mappability mask，列舉等可能候選（虛擬節點、等機率）。

**🔴 caption 修訂（必改）**：
- ❌ before：「incompatible 衝突(成環) → **使用甲基或其他資訊 確認 分群與歸類**」
- ✅ after：「incompatible（four-gamete 違反）→ 標 likely **CN-multiplicity/mapping artifact** + mappability mask + **不建單樹**；列舉等可能候選。**甲基不用於分群/排序**」
- 依據：程式 `candidate_scoring.py:32,34-37`（D4 fix 2026-07-01）本就是此行為（甲基已降 L3-weak non-resolver、`n_candidates="0(不支持單樹)"`）；舊 caption 是圖落後於碼。

**兩個 provenance 註記（方法章講清楚，免混淆）**：
1. slide 3 的**骨幹樹不是 UPGMA**——是基因型向量 **perfect-phylogeny / laminar**（`topology_analysis.py`）；UPGMA 只用在 slide 4 的甲基分群。**別讓 reviewer 以為拿甲基距離建演化樹（那就循環）。**
2. 深標籤 `HP1-1-1 / HP1-2 / HP2-2` 是**樹重建的 Dewey 路徑**，不是 BAM haplotag（BAM 只有 `1/2/1-1/2-1/3`）。圖上區分「淺=tag、深=重建輸出」。

**CN 在 slide 3 的唯一角色**：只出現在 incompatible 的 **artifact 解釋**——gain 造成的 multiplicity 會偽造第 4 配子；CN 用來**排除混淆**，不建樹、不分 read。

### §8.3 Slide 4 敘述 — 甲基 read×read 結構 + 雙軌驗證（characterization 層）

**主敘述（speaker note）**：甲基是**有界的資訊/註記層**——在 **label-first** 下檢查甲基結構是否**對齊**既有遺傳標籤（characterize），**不產生演化標籤、不確認 subclone**。

**角色分工講法**：
- **sSNV 共現（slide 3）＝標籤源；甲基（slide 4）＝註記 / 對齊檢查。**
- **label-first**（遺傳標籤 → 查甲基對齊；非循環，因 HP 來自零甲基 phasing）≠ cluster-first-then-label（＝tumor-only **double-dip NEGATIVE**）。
- 動詞紀律：甲基「**describe / corroborate / L3-flag**」，不「**confirm / discriminate subclone**」。

**CN 在 slide 4 的角色**：**固定參考層**——算 k_ISM 相對 k_CN 的 **excess-over-CN**、LOH masking、confound 排除。**k_ISM ⊥ k_CN（ρ=−0.038）**證甲基非 CN 鏡子（有獨立內容），但 CN **不是 per-read discriminator**。

**參數 / 口徑註記**：
- 「±5000bp」標成 **production 參數（預設 1000）**。
- **顯著 PERMANOVA ≠ subclone**（可為 cis-ASM / CN）；甲基作為 hard filter 為 DEAD——只做 characterize。
- 「共同分類與解釋」的 `?/?` 格＝誠實的 **undetermined**（定不出來即答案），**保留**。

### §8.4 三種合法「輔助後續分析」用途（slide 4 收尾可點）

- ✅ **負篩**：甲基同質 → 降優先。
- ✅ **L3 假設 flag**：遺傳群內部還有甲基次結構 → 標候選，建議 single-cell / multi-region 跟進。
- ✅ **characterization 註記**：每個遺傳定義 lineage 掛甲基 β-profile 做描述。
- ❌ 禁：「甲基群＝subclone」「甲基確認譜系」「演化標籤來自甲基」。

---

## 關聯文件

- 定位主文：`InterSubMod/docs/method_comparison/20260630_ism_positioning_vs_prior_work_01.md`
- baseline 詞典（confirm vs describe、cis-ASM vs lineage）：`InterSubMod/docs/methodology/20260621_baseline_discipline_conceptual_definitions_01.md`
- CN/SV 可行性 SoT：`InterSubMod/docs/plans/20260620_ont_cnv_sv_subclone_verification_feasibility_01.md`
- k_ISM×CN per-read 實測：`InterSubMod/docs/experiments/in_progress/2026/06/20260621_kism_vs_cn_perread/20260621_kism_vs_cn_perread_result_01.md`
- CLP CN 缺口 + 需 SAVANA：`InterSubMod/docs/methodology/20260629_clp_cn_crossval_integration_01.md`
