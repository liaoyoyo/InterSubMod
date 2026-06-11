<!--
建立時間: 2026-06-04
狀態: in_progress (方法解釋例子的檢核標準 + skill spec 參考；workflow wf_abcc0746-388 產出整理)
報告類型: method_explanation_criteria_and_skill_spec
受眾: 廖子游 + 未來建 skill (/skill-creator) + 論文 methods 圖製作
data_sources:
  - research/tsg_promoter_asm_reviewer/genome_survey_v2/control_cohesion_cistest.json   # β 0.25/0.228/0.108, d_cis −0.142, d_drift −0.022, n_shared_cpg 193, ARI 0.79, silhouette
  - research/tsg_promoter_asm_reviewer/genome_survey_v2/copy_partition_confirm.json      # d_within −0.023 (p=0.022), d_copy −0.11
  - /big7_disk/.../bip8_output_archive/20260119_all-with-w5000_1/.../chr13_32315128/clustering/significance.json  # Cramér's V 0.80/0.51, 56 reads, 229 CpG
  - research/tsg_promoter_asm_reviewer/scripts/15_genome_asm_pivot.py / 34_control_loci_cohesion_cistest.py / 37_copy_partition_confirm.py
  - src/core/ReadParser.cpp:128-138,150 / src/core/NormalBaseline.cpp / SignificanceAnalyzer.cpp
provenance_note: 12 criteria + 12 example requirements + skill spec 由 workflow wf_abcc0746-388（15 agents：web研究→定標準→三段draft→review→fix→re-review→skill）產出。所有 metric 為已驗證真值（grep 自上列 JSON/significance）；⚠ 部分 JSON 為 gitignored 工作樹輸出（in_progress 可 grep 即可；轉 validated/PI 前須確認落檔 + frontmatter data_sources 供 number_provenance_check.sh）。
-->

# 方法解釋例子：檢核標準 + 繪製要求 + skill spec
## Δβ somatic-cis（新法）vs 無監督分群 + Cramér's V（舊 ISM）— 以 BRCA2 chr13:32,315,128 貫穿

> **這份是什麼**：把「怎麼產生並驗證一個好的方法解釋例子」整理成可重複的檢核標準 + 例子要求 + skill spec。三段例子的**已渲染版**＝同層 `method_explainer_dbeta_vs_ism_v1.html`（2 SVG，已通過先前驗證）。
> **產製**：workflow `wf_abcc0746-388`（web 研究好例子特徵 → 定 5 維檢核標準 → 三段 draft→review→fix→re-review 迭代 → skill spec）。
> **迭代誠實結果**：A（Δβ法）PASS · C（新舊對照）PASS · **B（舊 ISM 段）FAIL** — fixer 過度矯正，把已驗證真值 blank 成 `{{待填}}`（誤把「不在白名單」當「必須刪真值」）。**§13 要的是「數字可 grep 到來源」，不是「未授權就刪真值」** → 修法＝回填真值 + 標來源檔。

---

## 1. 三段內容（已渲染例子見 HTML explainer）

| 段 | 內容 | 迭代 verdict |
|---|---|---|
| **A** | Δβ somatic-cis：read×CpG 如何整合成一個判斷（per-read 甲基→per-CpG 群率→跨 CpG 配對 Wilcoxon + normal cis）| ✅ PASS |
| **B** | 舊 ISM（w5000 binary）做法：10kb region + k=4 無監督分群 + Cramér's V + gating；C++ 解析/丟棄行為 | ⚠ FAIL（真值被過度 blank，需回填）|
| **C** | 新舊差異與**互補**：舊擅長任意 pattern（含交叉型）、新擅長方向性 somatic-cis + de-confound + 可排序 | ✅ PASS |

> 已渲染例子：`InterSubMod/docs/presentations/in_progress/20260603_methyl_haplotype_story_deck_v2/method_explainer_dbeta_vs_ism_v1.html`（圖1 Δβ生物示意 + 圖2 舊 vs 新；數字已驗證、read 標記為示意）。

---

## 2. 檢核標準 rubric（12 條 / 5 維度）

### 生物忠實 (BIO)
- **BIO-1** read→CpG→haplotype 三層關係正確；HP1-1 是 **HP1 germline 下的 somatic-allele 子標**（非獨立第三 haplotype）；同帶 label flip：**HP1-1（分析側）≡ HP2-1（舊 ISM/圖側）** 須一致或註明。FAIL：CpG 當 read 屬性 / HP1-1 當第三 haplotype / 甲基方向標反。
- **BIO-2** β＝某 haplotype 群在某 CpG 的甲基比例 (0-1, per-CpG 群率)；Δβ=β_som−β_germ 有方向。引數字須與真值一致（normalHP1 **0.25** / tumorHP1 **0.228** / tumorHP1-1 **0.108** / Δβ **−0.122**）。

### 工程科學嚴謹 (ENG)
- **ENG-1** Δβ 五步齊全且順序對：(1) per-CpG 兩群各算 β (2) 取 shared CpG (3) per-CpG Δ (4) mean=Δβ 主鍵 + max_abs_delta 副鍵 (5) **paired Wilcoxon over CpG**；說明對 read 數不均穩健。FAIL：漏 shared-CpG / Wilcoxon 講成對 read / 漏 max|Δ| 副鍵。
- **ENG-2** 舊 ISM 描述正確：10kb region + k=4 無監督分群 + Cramér's V(cluster↔hp / cluster↔alt) + gating；真值 **V hp=0.80 (p=0) / V alt=0.51 (p=0.0005) / passed_gating / 56 reads / 229 CpG**；指出 `germline_hp_only` 丟 somatic 子標的後果。
- **ENG-3** 每個數字可追溯（L1 實測），無新編造：因果分解 **d_cis=−0.142 / d_drift=−0.022 NS / d_somatic=−0.12 / d_within=−0.023 (perm p=0.022) / d_copy=−0.11**、n_shared_cpg **193**、silhouette **0.267**、blind ARI **0.79**。

### 解釋清晰 (CLR)
- **CLR-1** 一圖一概念、先直覺後公式、5 秒讀懂（直覺：「比同一人 germline 那條 vs 長了 somatic 突變那條，在同批 CpG 上甲基差多少」；公式 Δβ 後置）。
- **CLR-2** 新舊差異用對照結構呈現**互補**（非取代）：舊=無方向類別關聯（混 germline-allelic+copy、tumor-only、高 V 在常見 germline-ASM 也成立→BRCA2 淹沒）；新=方向性+somatic-controlled+normal-anchored+copy-partition+可排序。

### 誠實 (HON)
- **HON-1** 示意 vs 真值明標：read 甲基標記＝「示意，方向已放大」、聚合數字＝「實測真值」，視覺與文字可區分。
- **HON-2** 不過度宣稱：明標 **BRCA2 真 cis 小（d_within=−0.023 marginal, p=0.022）/ 單樣本 ★3 / 不宣稱 strong cis-driven**；不可把 d_cis=−0.142 當強因果而隱去 d_within。
- **HON-3** 誠實標 **Δβ-mean 交叉型盲點**：前後反向 pattern → Δβ mean≈0、paired Wilcoxon p=0.438 → 主排序漏；但 max|Δ|=0.583 + blind ARI=0.544 → 分群層抓得到 → 補 ARI/max|Δ|/PERMANOVA 雙鍵（＝工程任務 #19 動機）。

### 完整 (CMP)
- **CMP-1** A/B/C 三段皆有實質且自洽（read×CpG 整合 / ISM 做法 / 差異互補）。
- **CMP-2** 同一 **BRCA2 位點 chr13:32,315,128 G>A 貫穿新舊兩法**對照（舊法 56 reads/229 CpG/V=0.80 淹沒 vs 新法 n=193/Δβ=−0.122/可排序突出）；交叉型盲點以**合成模擬例**補充（與 BRCA2 真實例區分）。

---

## 3. 例子繪製要求（12 條）

1. 畫 **read×CpG 矩陣**：行=reads（依 HP 分組，至少 HP1=germline + HP1-1=somatic 子標兩組）、列=CpG 座標、格子=甲基狀態（實心=甲基/空心=未甲基）→ 讓人看出 β＝某列某組 reads 的甲基比例。
2. read 標記必標「示意：方向已放大」；聚合數字必標「實測真值 + 來源」，兩者視覺可區分（示意=淺色/虛線框）。
3. β/Δβ **一圖一步**：①per-CpG 兩群算 β ②取 shared CpG ③per-CpG Δ ④mean=Δβ + max|Δ| ⑤paired Wilcoxon(對 CpG)，每步單一概念。
4. **同一 BRCA2 真實位點** chr13:32,315,128 貫穿；真值一致：normalHP1 0.25 / tumorHP1 0.228 / tumorHP1-1 0.108 / Δβ −0.122 / n_shared_cpg 193 / blind ARI 0.79。
5. 舊 ISM 畫 k=4 無監督分群 → cluster↔hp/alt 列聯 → Cramér's V + gating；標 BRCA2 Region 真值（56 reads/229 CpG/V hp=0.80 p=0/V alt=0.51 p=0.0005/passed_gating）。
6. 新舊用**並列對照**呈互補（舊無方向混淆 vs 新方向性+somatic-controlled+normal-anchored+copy-partition+可排序）。
7. 畫/標 **交叉型盲點**（合成模擬，明標非 BRCA2）：HP1 前低後高/HP1-1 前高後低 → Δβ mean≈0、p=0.438；對照 max|Δ|=0.583/ARI=0.544。
8. HP1-1 圖例明標「HP1 germline 下 somatic-allele 子標，非第三 haplotype」；沿用 HP2-1 命名須註明 **HP1-1≡HP2-1（label flip）**。
9. 誠實區塊：明寫「真 cis 小 d_within=−0.023 marginal (p=0.022, n=19/19)、單樣本 ★3、不宣稱 strong cis-driven」。
10. 每圖開頭一句白話直覺（先直覺後公式），過 5 秒測試。
11. A/B/C 三段各有對應視覺/段落，缺一不合格。
12. 引 C++ 行為須正確標：`ReadParser.cpp:128-138` 解析 HP1-1/2-1/3、`:150` `germline_hp_only` 把 1-1/2-1/3 當 unphased 丟棄（＝舊管線丟 somatic 子標、需任務 #19 接回的根因）、`NormalBaseline.cpp` 提供 normal 基線供 d_drift。

---

## 4. 好例子特徵（web 研究）+ 來源 + pitfalls

**好例子特徵（節選，完整 12 條見 workflow 輸出）**：先直覺後公式 · 真實分子/read 當主角（非只有聚合條）· haplotype 分軌共用 x 軸與色標 · 一圖一概念 · 軸/單位標清楚 · 色盲友善 + 形狀冗餘（實心/空心 + 紅/藍）· 誠實分示意 vs 真值 · 顯示不確定性/read 數 · 漸進抽象（raw→modified-base→統計→聚合）· read 排序有語意 · 尊重 IGV/modbam 慣例 · caption 可獨立讀懂。

**來源（6, web 驗證）**：
1. **Methylartist**（Bioinformatics, PMC9154218）— read-level 甲基圖範本：基因模型→每 read CpG 實心/空心→reduced modified-base space→log-likelihood→sliding-window；WhatsHap HP/PS 分 haplotype。
2. **Modbamtools**（bioRxiv 2022.07.07.499188）— 每 read 一條 bar、紅=甲基/藍=未甲基；HP tag 分群 + read clustering 顯亞群。
3. **Ten Simple Rules for Better Figures**（Rougier et al., PLOS Comput Biol 2014）— 一圖一論點、色盲友善、避 chartjunk 的 canonical checklist。
4. **The Science of Visual Data Communication**（Franconeri et al., PSPI 2021）— 比較子集是慢任務→haplotype 差異靠空間對齊 + preattentive 編碼。
5. **MethPhaser**（Nat Commun 2024, s41467-024-49588-0）— 把方法步驟 schematic 與真實結果分開呈現的對照範本。
6. **Computational analysis of DNA methylation from long-read sequencing**（Nat Rev Genet 2025, s41576-025-00822-5）— single-read plot / haplotype-specific count 標準視覺化指引。

**常見爛例子 pitfalls（節選）**：只給聚合藏單分子（→看不出低覆蓋 artifact，對應 median|Δβ|<0.10/看 excess-over-null 教訓）· 示意當量化（→§13 捏造紅線）· haplotype 不對齊/色標不一致 · rainbow/紅綠無冗餘 · 一圖塞太多 · 軸未標/截斷誤導 · 先公式後直覺 · 忽略不確定性/單樣本當通則 · 不同口徑（read-instance vs unique / Tmode vs paired）並列。

---

## 5. Skill spec（可交 /skill-creator 建正式 skill）

**name**: `methods-example-generate-and-verify`
**when_to_use**: 為「分析方法」產生/檢核解釋例子或圖（尤其新舊方法互補對 PI 講、read-level haplotype 圖）。SKIP 純資料圖、PPTX 微調、無關 diagram。
**6 步流程**：
- **STEP 0 LOCK-AND-GATHER**（§13.0）：鎖定方法 + 一個真實 anchor 位點（BRCA2 chr13:32315128）貫穿新舊；每個數字問「現在哪個檔 grep 得到」；分析跑完→落 json/tsv→Read 回；無來源數字標 TODO 不憑記憶；標 C++ ReadParser:128-138 + :150 germline_hp_only 丟子標。
- **STEP 1 GENERATE DRAFT**：A read→β/Δβ、B 舊 ISM、C 新舊互補；一圖一概念、先直覺後公式、HP1 vs HP1-1 共用軸/色標 + 形狀冗餘；read=示意、數字=實測或 TODO。
- **STEP 2 VERIFY 12 criteria**（BIO/ENG/CLR/HON/CMP）：FAIL 條件＝CpG 當 read 屬性 / HP1-1 當第三 haplotype / Δβ 無方向 / Wilcoxon 對 read / 數字 grep 不到 / 示意當量化 / d_cis 隱去 d_within / 漏交叉型盲點。
- **STEP 3 FIX minimally**：⚠ **「不在白名單」≠「必須 blank 已交叉確認的真值」**；§13 要 traceable 不是 delete → 回填真值 + 標 src + frontmatter data_sources，否則保留 TODO；**分析 Bash 與 例子 Write 永不同 batch**。
- **STEP 4 RE-VERIFY 全 12**：先前 FAIL 已解或誠實降為示意；無新不可追溯數字；A/B/C 自洽。
- **STEP 5 FINALIZE/RECONCILE**：定義不共用一張圖；label 撞號標來源（如 max|Δ|=0.583 撞 chr7:88569202）；過 5 秒測試；frontmatter data_sources 完整供 `number_provenance_check.sh`。
**anti-patterns**：只給聚合藏單分子 · 示意未標當量化 · **過度矯正 blank 真值（本次 B FAIL）** · rainbow/未對齊 haplotype · 軸超載/未標 · 先公式後直覺 · 單樣本當跨樣本 · HP1-1 當第三 haplotype / CpG 當 read 屬性 / Wilcoxon 對 read · d_cis 隱 d_within · 多定義擠一圖。
**overall_verdict**: A PASS / C PASS / B FAIL（過度矯正把 9 組交叉確認真值 blank 成 TODO；§13 反捏造 ≠ 反真值）。值得建 skill：可重複的 generate→verify→fix→re-verify→finalize loop（先直覺後公式 + 示意/真值分離 + 互補性 + 每數字可 grep）。

---

## 6. 本次迭代的誠實發現（lesson + 待辦）

1. **B 段 FAIL 根因＝過度矯正**：fixer 把已驗證真值（0.25/0.228/0.108/−0.122/193/ARI 0.79/V 0.80/56/229）blank 成 `{{待填}}`。**修法**：回填真值 + 標來源檔路徑（這些值我已直接 grep 自 JSON）。已渲染的 HTML explainer **沒有此問題**（真值已在）。
2. **Provenance 細節**：上列真值多在 **gitignored 研究 JSON**（工作樹可 grep，但未 commit）。in_progress 可接受；**轉 validated/PI 前須確認 JSON 落檔 + frontmatter data_sources**（否則 `number_provenance_check.sh` 對 validated 路徑會 exit 2）。
3. **口徑差異（重要）**：cis-test 的 `mean_beta` = **0.25/0.228/0.108**（control_cohesion_cistest.json，給 d_cis=−0.142 的那組）；deepdive 另有 **5mC-only fraction germline 0.627 / somatic 0.592**（不同量、不同正規化）→ **兩組不可並列**；explainer/報告用 cis-test 的 0.25/0.228/0.108 並標明口徑。
4. **label 撞號**：合成交叉型 max|Δ|=0.583 與真實位點 chr7:88569202 dbeta_HP=0.583 同值 → 標來源避免誤讀。

---

## Provenance footer
- 12 criteria + 12 requirements + skill spec：workflow `wf_abcc0746-388`（15 agents）；web 來源經 researcher agent 檢索。
- 全 metric 為已驗證真值，grep 自 frontmatter data_sources 列之 JSON/significance/C++ 行號；read 級為示意。
- 三段已渲染例子：`method_explainer_dbeta_vs_ism_v1.html`（先前驗證 PASS）。
- in_progress 草稿；建正式 skill → `/skill-creator`；接回 binary → 工程任務 #19。
