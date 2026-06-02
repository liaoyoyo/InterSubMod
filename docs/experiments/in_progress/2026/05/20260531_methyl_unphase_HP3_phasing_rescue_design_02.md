<!--
建立時間: 2026-05-31
作者: Claude (Opus 4.8) + InterSubMod Research
報告類型: design spec (brainstorming → writing-plans 前置設計文件)
任務類型: A (exploratory pilot) — partial flag: HCC1395 單樣本 pilot，未跑全樣本
phase: design (brainstorming HARD-GATE 已過：設計呈現+用戶 ack)
framework: ADR (Architecture Decision Record) + Pre-mortem
provenance:
  - scoping: InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_unphase_HP3_phasing_rescue_scoping_01.md
  - workflow1: wf_41624212 (6 Explore + 1 綜整, 卡關根因 + 可行性矩陣)
  - workflow2: wf_7cbe71fc (5 深讀 + 2 對抗驗證 + 1 綜整, longphase-S inheritance 逐位元組 + method decision ledger)
  - 用戶拍板: 2 輪 AskUserQuestion (success metric / 改哪邊 / 首版範圍 / scope / MD1 / Phase A 工具)
partial_flag: true
-->

# 甲基分群拯救 unphase/HP3 與更遠 read 的相位 tag — Design Spec v2

> **狀態**：design 定案，待用戶 review → 通過後 invoke writing-plans 拆實作計畫。
> **本文件職責**：鎖定設計決策 + 完整 method decision ledger（可檢核）+ Phase A 詳細規劃。探索證據在 scoping_01。

---

## 1. 任務定位（鎖死，不可漂移）

**一句話**：在 longphase-S 既有 phasing/tag 已盡力後，對它**放棄或救不動**的 read（unphase / 卡住的 H3），用**甲基分群**嘗試再救一層相位歸屬，並量化這件事能不能做、做得多好。

**success metric（鎖死 — 這是「新方向」vs「滑回已死路」的唯一分水嶺）**：
- ✅ **phasing rescue accuracy**：被救回的 read 與 ground-truth/anchor haplotype 一致率、被救回的 phaseable read 數/比例。
- ✅ **對 LOH-constrained phasing 活軸（NG=2 same-hap）的增益**（嚴格分離觀察）。
- 🔴 **絕不可回頭包裝成 variant TP/FP F1 改善** — 一旦這樣 framing 立即落入已 concluded 死路（甲基化當 FP filter，⭐2 L4），需 Productive-Failure C1/C2/C3 reopen 條件。

**兩條紅線**（pre-mortem 防呆）：
1. metric 鎖在 read tag / phasing 層級，報告與 spec 任何地方不得宣稱改善 variant F1。
2. **ALLELE-axis confound**：HP-axis（HP1 vs HP1-1，somatic-controlled）為主證據軸；ALLELE-axis（ALT vs REF read 甲基差）受 germline-het 基線 ASM confound（memory: TP 11.1% < het null 15.2%），只作描述且強制 germline-het negative control。

---

## 2. 已坐實的關鍵事實（逐位元組 + 雙重對抗驗證 confirmed）

### 2.1 longphase-S 軟體身分（更正前一輪 agent 誤判）
- **= `/big8_disk/liaoyoyo2001/longphase-s`**（LongPhase-S，擴展 LongPhase v2.0.0 的 somatic 變體），**不是** longphase-to-mod。
- ISM paired pipeline `01_longphase_s.sh` 用它 `somatic_haplotag` 模式產 `{SAMPLE}_tagged.bam` → ISM step02 `inter_sub_mod` 消費。
- ⚠ pipeline 實際 binary = `config LONGPHASE_S_BIN`（/big8_disk/.../Knowledge/codebase/longphase-s/，Feb 7 build）；源碼在 `/big8_disk/liaoyoyo2001/longphase-s/src/`（模組化）。

### 2.2 longphase-S 既有 read→HP 決策（HaplotagStrategy.cpp:452-602, percentageThreshold=0.6）
- 每條 read 由 `judgeSomaticReadHap` 一次決定，吐 unTag/H1/H2/H3/H4/H1_1/H2_1/H1_2/H2_2。
- 用整數計數 hpCount[1..4]（germline H1/H2 + somatic H3/H4）算 `tumHPsimilarity` / `norHPsimilarity` = max/(max+min)，與 0.6 比較。
- **H3 成因**：有 somatic 變異（tumHPsim≥0.6）但 norHPsim<0.6（germline 方向不明）。
- **unphase(unTag) 三來源**：tumHPsim<0.6 / 跨 phase block（norCountPS.size()>1）/ 無任何變異。
- unTag 與 H3 **都仍寫進 BAM**（HP:Z 字串："." / "3" / "1-1" / "2-1"）。

### 2.3 inheritance 只救 H3、完全不碰 unphase（SomaticHaplotagProcess.cpp:461-527）
- `inheritHaplotype` 僅當第一段==H3 才呼叫；用 per-position `somaticReadDeriveByHP` 投 deriveByH1 vs deriveByH2，`deriveByHpSimilarity`=max/(max+min)，≥0.6 升 H1_1/H2_1，否則（含平手0.5、無證據）留 H3。
- **SSRS 實測基線（HCC1395, 120,834 somatic 位點, read-instance 口徑, 守恆已驗證）**：

| read-instance | before inheritance | after | 變化 |
|---|---:|---:|---:|
| H1-1 | 1,296,335 | 2,595,184 | +1,298,849 |
| H2-1 | 1,007,638 | 1,531,667 | +524,029 |
| H3（卡住）| 4,343,202 | 2,520,324 | **−1,822,878（救 42%）** |
| **untag(unphase)** | 1,777,511 | 1,777,511 | **0（完全不變）** |

- 守恆 ✓：H1-1+H2-1 增量 = H3 減量 = 1,822,878。
- 位點級 DeriveHP：H0（最難）65,033 (53.8%) / H1 33,947 / H2 21,854；H0 子集 derive_by_H1_H2 25,653。

### 2.4 甲基缺口（雙重 confirmed）
- longphase-S 的 haplotag/somatic voting/inheritance **零甲基**（全 body + 全目錄 grep `bam_parse_basemod`/`methWeight`/`WithMeth` 零命中）。
- 唯一解 MM/ML 的是**解耦的 modcall 子命令**（functional production-grade htslib `bam_parse_basemod`，5mC qual≥204=mod / ≤51=unmod，但對外 `//` 隱藏、輸出走 VCF 不回 per-read→per-HP）。
- 把甲基接進 longphase-S inheritance = **從零接 O(300-500) 行 C++**，無現成 hook。
- 源 BAM = `ONT_5khz_simplex_5mCG_5hmCG` → **甲基已在 BAM MM/ML，無需重跑 basecalling**。

### 2.5 RR/RA/AR/AA 原料現況
- per-read 有序 (pos, baseHP) 序列存在 `ReadVarHpCount.posHpPairs`（SomaticVarCaller.cpp:434-446）。
- ⚠ baseHP 把 somatic-ALT 編成 SOMATIC_H3，**無顯式 somatic-REF token**；無跨「一對 somatic 位點」的 4-way joint 結構 → 需 ISM/Python 端 derive。`_alleleCount.out`(39,447 變異) 給 per-position 基礎。

---

## 3. 設計決策（ADR）

### ADR-1：救援目標分三層 + 優先序（取代「救全部」）
基於 SSRS 鐵證，longphase-S 留下三類待救援族群，優先序：
1. **(P1) unphase read**（longphase-S 主動放棄、最乾淨、無既有方法競爭）— 最大機會空間。
2. **(P2) inheritance 後仍卡 H3 的 read**（deriveByHpSimilarity<0.6）。
3. **(P3) H0 難定相位點**（65,033 個 / derive_by_H1_H2 25,653 子集）。
- distant read（跨中介傳遞）：列**後期獨立階段**，首版不做重指派（最高紅旗：定義未定 + transitive 累積誤差 + 對抗 longphase 跨 PS-block 保守設計）。

### ADR-2：介入層 = C（先 ISM/Python 旁路，過 GO 門檻才開 C++）【用戶拍板 MD1=C】
- **Phase A = 純 Python 旁路**（連 ISM C++ 都不改）：讀 paired tagged BAM 的 MM/ML + HP tag，直接算甲基 affinity 判別力。
- 只有 Phase A **GO**（甲基 affinity AUC>0.58 且過 germline-het null + /auc-confound-guard 三關）才升級到 C++ 增強（B 路徑，另開 spec）。
- 理由：用戶明示「先 ISM 原型低風險可逆」+「先 paired」；C++ 是 Hard Gate（必編譯）+ 278GB BAM 重產不可逆；unphase 是乾淨空間 post-hoc 即可驗證。

### ADR-3：信心量化 = 現成 longphase-S 訊號當 baseline + 自定甲基信心（不重新發明）
見 §5 confidence metrics shortlist。

---

## 4. Method Decision Ledger（用戶要求：每個方法分岔記錄觀察數據 + 量化判準 + 決定，供檢核驗證）

> 此表是本任務的**可檢核抉擇核心**。每跑完一步觀察，回填「實測值」欄並標決定是否符合判準。

| ID | 方法分岔 | 選項 | 要觀察的數據 | 量化判準 | 當前決定 | 理由 | 實測值（待回填）|
|----|---------|------|------------|---------|---------|------|------------|
| **MD1** | 介入層 | A ISM旁路 / B longphase-S C++ / C 先A後B | unphase+卡H3 子集甲基 affinity AUC；C++ 工程量(O(300-500)行+Hard Gate) | AUC>0.58 過 het null 才升 B；≤0.58 兩路皆 NO-GO | **C（已拍板）** | 用戶先 ISM 低風險；C++ 不可逆 | — |
| **MD2** | 甲基證據軸 | A 純HP-axis / B 純ALLELE / C HP主+ALLELE描述 | HP1 vs HP1-1 Δβ(somatic-ctrl) vs germline-het ALT-REF null Δβ | HP-axis OR/AUC 須>het null(memory OR=1.79)；ALLELE 不過 null 只描述 | **C** | feedback_asm_allele_axis_baseline_confound | — |
| **MD3** | per-read 信心指標 | A deriveByHpSim / B PQ / C germlineVarSim / D 組合 | 三訊號在 unphase/H3 子集分布與 dynamic range | 須覆蓋 unphase(deriveByHpSim 對 unphase 無定義) | **A 作H3 baseline + 自定甲基信心覆蓋 unphase** | deriveByHpSim 僅 H3 適用，unphase 須新信心 | — |
| **MD4** | 甲基特徵表示 | A binary / B raw ML matrix / C binary先+raw備 | binary vs raw 判別 AUC 差；noise band(51-204) read 佔比 | raw 比 binary 高>0.02(Cohen) 才用 raw | **C（binary 先 pilot）** | modcall binary 閾值 production-validated；ISM 已有 Read×CpG matrix | — |
| **MD5** | 救援範圍 | A 只救卡H3 / B 只救unphase / C 兩者(unphase優先) | unphase vs 卡H3 各自 unique read 數、甲基覆蓋率、affinity | 各子集 affinity AUC 分別評 | **C（unphase 優先）** | unphase 是乾淨無競爭空間；機制不同須分開 | — |
| **MD6** | 重指派接受門檻 | A 沿用0.6 / B ROC Youden / C 雙門檻 | 甲基信心 vs ground-truth(germline_phased HP) ROC | cutoff 須使 precision≥已知 HP read baseline 一致率 | **B（ROC 決定）+ 並列 0.6 對照** | 需可量化抉擇依據，避免任意門檻 | — |
| **MD7** | LOH/CN read gate | A 不gate / B 標記但仍救 / C 排除高LOH | LOH vs non-LOH affinity 差；LOH 是否 cnLOH 單套致 axis 退化 | LOH 顯著異於 non-LOH → 分層報告不排除 | **B（標記分層、仍救、觀察 LOH-phasing 增益）** | success metric 含 LOH-phasing 增益，排除會抹掉觀察目標 | — |
| **MD8** | 無甲基覆蓋 read | A 保留原tag標記 / B fallback / C 排除分母 | unphase/H3 read 甲基 CpG 覆蓋數分布；零覆蓋佔比 | 雙口徑報告（全 vs 有甲基覆蓋子集） | **A（保留+標 no-meth-evidence+雙口徑）** | 嚴禁分母操弄虛報 rescue 率 | — |
| **MD9** | 口徑/版本對齊 | A 用SSRS Tmode dump / B paired重產 / C SSRS基線+paired真口徑 | SSRS(Tmode HP:i:3) vs paired(HP:i:11/21/33) inheritance 規則一致性 | 兩 build 規則須驗證一致才可混用；最終 effect-size 須來自 paired | **C** | SSRS 是另一 build Tmode 非 paired；混用前須對齊 | — |
| **MD10** | HP-tag 細類統計 | A 擴充haplotag_qc.py / B 自寫 / C ISM RegionProcessor counter | genome-wide unique-read 9 類 count+frac | 須分 unphase/HP3/HP1-1/HP2-1（hp_other 混合不可用）| **A 或 C（擴充現成）** | haplotag_qc.py hp_other 混所有 somatic+HP3 不能當 HP3 數 | — |

---

## 5. Confidence Metrics Shortlist（現成 longphase-S 訊號優先）

| 指標 | 來源 file:line | 意義 + 適用範圍 |
|------|--------------|--------------|
| **deriveByHpSimilarity** | SomaticHaplotagProcess.cpp:503（max/(max+min), 門檻0.6@515）| ★per-read H3 重指派信心 baseline；值域[0.5,1]∪{0}。⚠ **僅 H3 適用**，unphase 無 derive 證據用不了 → 正是甲基須補的子集 |
| **PQ** | HaplotagStrategy.cpp:287（−10·log10(min/(max+min)), 純→40）| 現成 Phred-scale germline 信心；可覆蓋 unphase（deriveByHpSim 不能）|
| **somaticReadDeriveByHP + existDeriveByH1andH2** | SomaticVarCaller.cpp:1480-1505 | per-variant 對偶；existDeriveByH1andH2=true = 同一 SNP 的 H3 read 同時 derive 自 H1+H2 的歸屬混亂位點 = 甲基最該拆解對象 |
| **germlineVarSimilarity** | `_totalRead.out`（per-read 欄, 1.07GB）| 現成 dump 可零改碼當 baseline 信心對照 |
| **pure_H3_readRatio / Mixed_HP_readRatio** | SomaticVarCaller.cpp:548-551（[0,1]）| per-variant：pure_H3 高=純H3無derive=救援目標；Mixed 高=混亂 |
| **before/after inheritance 轉換率** | `_readDistri_before/afterInheritance.out`（21欄, 120,834 site, 守恆已驗證）| longphase-S 既有救援成效權威基線（救 42% H3、unphase 不變）；甲基成效 = 比這份的增量 |

---

## 6. Phase A 詳細設計（純 Python 旁路，零 C++，TDD 拆最小單元）

> **目標**：在零改碼下回答「甲基對 unphase/卡住H3 read 有沒有足夠 HP 判別力」這個 GO/NO-GO 問題。
> **資料**：paired HCC1395_tagged.bam（278GB, /big7.../paired_full/.../longphase_s/）+ germline_phased_merged.vcf.gz(22MB) + filtered_snv_tp/fp.vcf.gz。

### A0（前置，可逆）：HP-tag 9 類細分統計
- 對 paired tagged BAM 算 genome-wide unique-read 級 9 類 count+frac（不可用 haplotag_qc.py 的 hp_other）。
- **驗收**：unphase / HP3 / HP1-1 / HP2-1 各自數字 + 守恆檢查（總和=總 read 數）。
- **見樹也見林四層**：aggregate（全基因組）/ canonical（典型 region）/ extreme outlier（如 chrX 高 HP3）/ well-explained。

### A1：三類救援族群刻畫
- 對 unphase / 卡住H3 / H0 位點，量化 unique read 數、caller_af 分布、LOH/CN 位置、甲基 CpG 覆蓋率。
- **confound 前置檢查**：unphase/H3 是否同 Phase2 般 low-AF 主導（若是，甲基可能仍是 caller_af proxy）。

### A2（核心 GO/NO-GO gate）：甲基 affinity 判別力
- 對 unphase('.') + 卡住H3('3') 子集：以「已知 HP1/HP2 read（germline_phased het 投票或 SEQC2 truth）」為 anchor，算甲基 cluster affinity → HP1 vs HP2 的判別 AUC。
- **HP-axis 主軸**（somatic-controlled）；ALLELE-axis 僅描述。
- **強制 germline-het negative control**（HP1 vs HP1 same-allele null）。
- **過 /auc-confound-guard 三關**（within-group OLS + AF-bin + permutation）。

### A3：RR/RA/AR/AA 兩點 read 觀察（描述性）
- 從 posHpPairs + alleleCount derive 跨兩 somatic 位點 read 的 Ref-Alt/Alt-Ref 組合 + 各組甲基分群。
- 強制 germline-het negative control（ALLELE-axis confound 防線）。

### GO/NO-GO 判準（量化）
- **GO**（升級到 B 路徑評估，另開 spec）：unphase 與/或卡住H3 子集甲基 affinity **AUC>0.58 且顯著超過 germline-het null**（過 confound 三關）。
- **NO-GO**（兩路皆停，不投入任何 C++）：AUC≤0.58，或 AUC>0.58 但被 germline-het null 解釋。
- **PROBE**：AUC>0.58 僅特定子集（如 non-LOH）成立 → 限定 scope 重評。

### Hard Gate 標記
- Phase A 全程**零 C++ commit、無刪檔、無 NO-GO 研究方向判定** → **非 Hard Gate，可直接執行**。
- ⚠ AUC 顯著性宣告前**必過 /auc-confound-guard**（重蹈 ISM 多次 methylation AUC confound 教訓的防線）。
- ⚠ MD9：Phase A 若用 SSRS dump 當概念基線，最終 effect-size 口徑須來自 paired BAM 重產（需確認 paired 模式產 readDistri）。

---

## 6+. 第二輪用戶要求（2026-05-31 round 2）— R1-R7 納入追蹤

> 用戶補 7 項要求。R4/R7 框架本輪寫完整（方法學內部、不依賴外部）；R1/R2/R5/R6 寫方法、待文獻 workflow（wf_8a3d2e1c）+ real-data preview 補實；R3 由文獻 workflow 產出。

| # | 要求 | 狀態 | 落點 |
|---|------|------|------|
| R1 | unphase 也用甲基分 HP1/HP2，**嚴格驗證合理有效** | 方法定，驗證協議見 §6A2 強化 | §6A2 + §6+R1 |
| R2 | **同 tag 甲基空間分密集多群**（copy/alignment）→ 偵測+紀錄 | 新增偵測設計 | §6+R2 |
| R3 | 現狀 + **論文預敘述 + 差驗證** | 文獻 workflow 進行中 | 待 wf_8a3d2e1c |
| R4 | 資料**優先級+影響+可信度+轉 tag 可信度數值** | ✅ 本輪寫完整 | §6+R4 |
| R5 | 更多視覺化 + **數據限制與定義** | viz 清單定，待數據 | §6+R5 |
| R6 | **IGV + tag/甲基分群熱圖** + 驗證合理 | 工具✓，待 BAM 定位 | §6+R6 |
| R7 | 各步驟 vs **baseline delta + 可驗證** | ✅ 本輪寫完整 | §6+R7 |

### §6+R1 — unphase 重指派的嚴格驗證協議（用戶強調「小心驗證」）

unphase read 無任何 germline derive 證據（deriveByHpSimilarity 不適用），是甲基唯一可救但**最易過度宣稱**的子集。驗證三道關，全過才宣稱「救得合理」：
1. **Anchor 一致率**：對「有 germline 證據的已知 HP1/HP2 read」做 leave-out — 把它們當未知，用甲基分群預測，比對真 HP。一致率須顯著 > 0.5（隨機）+ 過 permutation。
2. **內部凝聚度**：unphase read 甲基分群的 silhouette / Davies-Bouldin，須顯示「真的有兩群」而非硬切連續分布（避免把 noise 切成假群）。
3. **負控**：germline-het null（同 allele 不應有 somatic 甲基差）+ 打亂甲基標籤的 shuffle null，重指派一致率須掉回隨機。
> 量化門檻：anchor 一致率 AUC>0.58 過 null + silhouette>0.25（弱結構下限，待文獻校準）才 GO。

### §6+R2 — 同-tag 甲基多群偵測（copy / alignment confound）⭐ 新洞察

用戶觀察：同一 HP tag 的 read 可能在甲基空間分成**密集清楚的多群** → 可能是 (a) copy-number gain / 多拷貝 paralog、(b) multi-mapping / repeat alignment artifact、(c) cnLOH / subclone 異質性。若不處理會污染重指派（把「拷貝差異」誤當「haplotype 甲基差異」）。
- **偵測**：對每個 (region × HP-tag) 群，跑甲基 modality test（Hartigan dip test / GMM BIC 比 1-cluster vs k-cluster）。多模態 flag。
- **分流歸因**：多模態群交叉檢查 — 局部 coverage 異常（>2× 周邊 = 多拷貝嫌疑）/ MAPQ 分布（低 MAPQ 富集 = mismapping）/ LOH-CN track（落多拷貝區）/ supplementary alignment 比例。
- **處理（用戶要「避免影響或清楚紀錄延伸處理」）**：(i) copy/alignment 來源 → 標記 `multimodal_artifact` 排除出重指派分母但**保留紀錄**；(ii) subclone 真實生物多群 → 標記 `multimodal_biological` 分層分析（連結 LOH-phasing 活軸）。雙口徑報告。

### §6+R4 — 資料可信度框架 + tag 可信度數值轉換（用戶要求）✅

**(a) 資料/資訊優先級 + 可信度分級**（高→低，重指派證據權重）：

| 級 | 證據 | 可信度 | 影響權重 | 驗證方式 |
|----|------|:---:|:---:|------|
| T1 | germline het SNP 直接投票（norHPsim≥0.6）| 最高（直接 phase）| 1.0 | longphase-S 既驗 |
| T2 | somatic derive-by-HP（deriveByHpSim≥0.6）| 高（跨read共識）| 0.8 | inheritance 守恆 |
| T3 | **甲基 affinity（本任務新增）** | **中-低（待 Phase A 證）** | **0.3-0.6 待定** | anchor 一致率 + null |
| T4 | 鄰近 read / PS-block 傳遞 | 低（transitive 累積誤差）| ≤0.2 | distant read 後期 |
| X | copy/alignment multimodal | 不可信（confound）| 0（排除）| §6+R2 偵測 |

**(b) tag 可信度數值轉換公式**（把異質證據統一成 [0,1] confidence）：
- 有 germline/derive 證據：`conf = max/(max+min)`（沿用 longphase-S，T1/T2）。
- 純甲基重指派（T3，unphase）：`conf_methyl = w · affinity_margin · coverage_factor · (1 − multimodal_penalty)`，其中 affinity_margin = |dist_HP1 − dist_HP2|/(dist_HP1+dist_HP2)、coverage_factor = min(1, n_CpG/n_CpG_min)、w 由 anchor ROC 校準。
- **兩套 confidence 分開存**（`phase_conf` vs `methyl_reassign_conf`），不混用；最終 tag 標來源 tier（T1-T4）+ 數值。
> 量化判準：methyl_reassign_conf 的校準須使「同 conf bin 的 anchor 實際一致率」單調對齊（calibration curve，reliability diagram 驗證）。

### §6+R7 — Delta-tracking 框架（每步 vs baseline 可驗證）✅

用戶要求：每個動作對「原始數據 / 現況 baseline」的差異都要可驗證、有數據支持。

**Baseline 定義（三層，凍結）**：
- **B0 原始**：longphase-S tagged BAM 的原始 HP tag 分布（unphase/HP3/HP1-1/...）。
- **B1 既有 inheritance**：longphase-S 自己救完的 after-inheritance 分布（救42% H3、unphase不變 — 已驗證守恆）。
- **B2 甲基重指派後**：本任務每步輸出。

**每步必記錄的 delta 卡片**（schema）：
```
step_id / 動作 / 輸入(B?) / 輸出 /
  Δ各HP類count / 守恆檢查(救出=減少?) /
  正確性樣本(anchor一致率) / 對照(無甲基baseline) /
  異常數(multimodal/no-coverage) / verdict(合理?有數據支持?)
```
- **守恆強制**：任何「unphase→HP1/HP2」的數量，必等於 unphase 減少量（如 inheritance 守恆驗證 1,822,878 那樣）。不守恆 = bug，紅旗停。
- **可回溯**：每步 reassigned_hp_tag 不覆寫原 tag，保留 `original_hp` + `reassign_source` + `reassign_conf` 三欄，任何人可重算 diff。
- **每步 commit**（TDD），delta 卡片進報告 + evidence_ledger。

### §6+R5 — 視覺化清單（待數據補實）+ 每圖限制與定義

| 圖 | 類型 | 數據需求 | 限制/定義須標 |
|----|------|---------|------------|
| 既有 SVG 決策樹 / before-after / GO流程 | 概念 SVG | ✅ 已做 | read-instance≠unique read |
| **tag × 甲基 read 分群熱圖** | seaborn clustermap | per-read×CpG 矩陣（待BAM）| CpG 覆蓋稀疏定義、-1 缺值處理 |
| **IGV per-read 甲基 view** | igv-reports standalone | BAM region（待定位）| 顯示的 loci 選擇理由、coverage |
| anchor 一致率 calibration | reliability diagram | Phase A 跑完 | bin 數、CI |
| multimodal 偵測 | dip-test p / GMM BIC 圖 | per-group 甲基 | modality 判準閾值 |
| 三類救援族群 caller_af 分布 | violin | Phase A | low-AF confound 警示 |

### §6+R6 — IGV + 熱圖（工具就緒，待 BAM 定位）

工具已確認（2026-05-31）：`pysam 0.23`（解 MM/ML → per-read×CpG 矩陣）+ `seaborn clustermap`（甲基熱圖含 HP tag 側欄）+ `igv-reports`（standalone IGV.js HTML，含 modified-base track，無需 X display）。
- preview loci 選擇原則（待確認）：(1) 已知 ASM 位點（BRCA2 chr13:32,315,128）做 positive control、(2) 高 unphase 比例 region、(3) 已知多拷貝 region 做 §6+R2 multimodal demo。
- **阻塞**：真實 paired tagged BAM 路徑需重新定位（見 §6 ⚠）。

---

## 7. 後續 Phase（條件化，Phase A GO 才啟動）

- **Phase B（條件式 Hard Gate）**：longphase-S C++ inheritance 甲基增強（inheritHaplotype:503-515 對 deriveByHpSim<0.6 注入甲基 tie-break）O(300-500) 行 → 走 /methodology-audit → /cpp-change，另開 spec。
- **Phase C**：distant read 跨中介傳遞（先探索性可視化，定義「distant」語意後再評重指派）。
- **Phase D**：跨 7 樣本 LOSO 泛化（防 Phase2 Cycle1 in-distribution overfit 覆轍）。

---

## 8. Pre-mortem（這個設計可能怎麼失敗）

1. **甲基 AUC 撞 ≤0.58 ceiling**（pure-methylation 系統根因）→ Phase A 就 NO-GO，但這正是最便宜的擋點，沒浪費 C++。
2. **unphase read 甲基覆蓋太稀疏** → A1 會先量化，若零覆蓋佔比過高則 rescue 上限受限，MD8 雙口徑誠實報告。
3. **LOH 區 cnLOH 單套讓甲基 axis 退化** → MD7 分層處理，不混報。
4. **口徑混用**（SSRS Tmode vs paired）→ MD9 強制 paired 真口徑。
5. **滑回 F1 死路** → §1 紅線 + metric 鎖死防呆。

---

*探索證據：scoping_01 + workflow wf_41624212 + wf_7cbe71fc。本 spec 待用戶 review → 通過後 invoke writing-plans 拆 Phase A 實作計畫（TDD 最小單元 + 每步 commit + 驗證）。*
