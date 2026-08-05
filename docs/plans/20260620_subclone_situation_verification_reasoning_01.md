<!--
建立時間: 2026-06-20
狀態: 推論整理 + 驗證設計（pre-verification reasoning；用戶明示「先整理推論，再分別驗證各種狀況」）
報告類型: reasoning_consolidation + verification_design
受眾: 廖子游 · 碩論 subclone 主軸 · 後續「分別驗證各狀況」階段的 SoT
framework: Verdict-Pyramid + Assertion-Evidence (Tier 4 paper-scope)
task_type: B comprehensive validation 的 pre-reg 推論層（scope=全 30,490 位點盤點 + chr2/chr8 驗證單元）
data_sources:
  - docs/methodology/20260617_structure_label_situation_inventory_01.md (三狀況 + 六格 + §0.5-0.7，verify_breakdown 25/25 PASS)
  - docs/methodology/_assets/20260618_subcluster_pilot/{split_accounting,spectrum_decision,contingency_summary,fisher_cramersv,nosignal_breakdown}.json (切群總帳/三態，加總已驗)
  - docs/paper_focus/02_paper_framework/20260615_chr2_subclone_case_and_method_concept_01.md (chr2 9 步驗證模板 + ASM-control)
  - /big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed (chr8 LOH 範圍，本輪 awk 計算)
  - /tmp/perlocus_flags.csv (reclassify-v2 全 30,490 dump；⚠ 暫存未固化)
  - memory: project_tumor_only_axis_negative_subclone_classification / project_apriori_subclone_classification_model / project_chr2_18m_subclone_locus_verification / project_subcluster_cluster_count_determination
provenance_note: 所有數字引自上列 verified 文件 + 本輪 grounding workflow（wf_3b296dc8-143，加總驗算全 PASS）。本檔=推論整理，不產新數字；§13.0 撰寫與分析分離。判別率 7.6%/16.4% 等屬 L3（memory + /tmp JSON，未 committed），引用入投稿前須重算固化。
-->
<!-- provenance-verified: 切群總帳/三態/六格/四格加總全驗（workflow wf_3b296dc8-143）；chr8 LOH 範圍 awk over CNV bed；判別率/Jaccard 標 L3 memory 來源待固化。 -->

# Subclone 狀況分類 — 推論整理 + 驗證設計（pre-verification SoT）

> **這份是什麼**：把「以甲基距離、依標籤分割 30,490 個 somatic 位點」的整體結果，整理成 **(1) 資料現況確認 → (2) 能否清楚分類定位 → (3) 各狀況 ↔ clone/subclone 對應 → (4) 各狀況驗證設計（chr2 模板 + chr8 全-LOH 想法）→ (5) 推論合理性誠實裁決**。
> **用途**：作為下一階段「分別驗證各種狀況」的 SoT；本檔只整理推論，**不下已驗證結論**。
> **整體 tier**：盤點數字 = L2（verified 文件，加總已驗）；狀況→subclone 對應 = **L3 推論**；chr8 驗證 = **L4 假設+可行性已證**。

---

## L0 — 一頁裁決（先看這裡）

**問：能不能清楚分類、定位各狀況，並對應 clone/subclone、設計驗證？**

| 子問題 | 裁決 | 一句話 |
|---|---|---|
| 資料整體可信？ | ✅ **是**（L2）| 30,490 加總自洽，多重 partition 機械真值對得上 |
| 能「清楚分類定位」？ | 🟡 **有條件成立** | 只在**互斥 partition（§2 六格 / §0.7 四格）**成立；§0「三狀況」彼此**重疊**，不能直接定位 |
| 狀況② = subclone 最相關？ | 🟡 **方向對、口徑要收斂** | 對的是 tumor-only **方向**；但須收斂到 PATTERN(354) 非 LEVEL(7910)，且 85% 清楚多群落 germline/carrier 軸 = 可能 cis-ASM 非 subclone |
| chr2 模板可逐一驗證？ | 🟡 **可，但非「逐一」** | chr2 = 手做 L2 單位點；30,490 不可能逐一手做 → 只能 **canonical/outlier 抽樣** + 自動化（未證）|
| chr8 全-LOH 也好用？ | 🟡 **PARTIALLY_SOUND** | 資料齊（chr8 96% LOH、284 個 ≥3-sSNV 簇、normal 覆蓋齊）；但 LOH 是**必要鋪墊非充分**，且 chr8 = 7.4× FP 富集區（雙面刃）→ 鎖**高密度 40kb 子窗 + normal-control**，非「整段」|
| 推論整體合理？ | 🟡 **PARTIALLY_SOUND（合理但須 5 處修正）** | 見 §5 |

**一句總結**：**推論方向正確、資料健全，但「清楚分類」「逐一驗證」「狀況②=subclone」「chr8 整段」四個措辭都偏樂觀，須各自加界線**（詳 §5）。

---

## §1 資料整體情況確認（L2）

### 1.1 輸入與方法
- **輸入** = caller 輸出 somatic 位點總數（TP+FP）= **30,490 region**（單一 paired run，HCC1395）。
- **方法**：每位點 ±5000bp 窗 → read×CpG 甲基矩陣 → BERNOULLI 距離 → UPGMA 階層聚類（**完全無監督，標籤不進入**）→ 再與 a-priori 標籤（haplotag HP1/HP2/HP1-1/HP2-1、ALT/REF）做 cluster×label 列聯表檢定（CramérV + Fisher）。

### 1.2 三套數字系統（🔴 不可混算）
盤點用了**三種不同切面**，各自加總自洽，但**框架不同不可並列**：

| 系統 | 框架 | 互斥？ | 用途 |
|---|---|---|---|
| **§0 三狀況** | 純度軸關係（對齊/一群多標籤/多群一標籤）| ❌ **重疊** | 描述「結構↔標籤」三種關係 |
| **§2 六格 / §0.7 四格** | 純度軸 × within-HP 軸 / T×N 軸 | ✅ **互斥（Σ=30,490）** | **定位某位點屬哪格** |
| **切群總帳** | silhouette best_k≥2 光譜 | ✅ 互斥 | 「能不能切 ≥2 群、多清楚」|

> ⚠ 還有第四套（Venn 互斥區 `contingency_summary.json`，doc §3 用，th_within=0.5）與切群總帳（silhouette best_k）數字不同（如「within-carrier split 351」vs「clear≥2 2941」），**各自自洽不可混**。

### 1.3 機械真值（加總已驗，workflow PASS）
- **三狀況**（重疊）：①對齊 22,497（73.78%）/ ②一標籤多群 8,077（26.49%）/ ③一群多標籤 510（1.67%）
- **§2 六格**（Σ=30,490）：對齊×單群 16,317 / 對齊×多群 6,180 / 一群多標籤×單群 371 / ×多群 139 / 無群×單群 5,725 / 無群×多群 1,758
- **§0.7 四格 T×N**（Σ=30,490）：①tumor 結構未歸因 4,698 / ②both 7,066 / ③僅 normal 9,182 / ④無結構無歸因 9,544
- **切群總帳**（Σ=30,490）：**切不出來 24,493（80.3%）**＝不足 412＋一團同質 24,081；**切得出來 5,997（19.7%）**＝clear≥2 **2,941（9.6%）**＋中等 765＋弱 506＋極弱 1,249＋退化 536
- **clear≥2 對齊軸**（分母 2,941）：germline/carrier **2,500（85%）** / HP1-HP2 441（15%）/ 隨機 **0**
- **三態決策**（Σ=30,490）：不能問 12,039（39.5%）/ **確認 1 群 13,855（45.4%）** / **多群候選 4,596（15.1%）**

### 1.4 誠實底線（§7，必守）
1. **「是訊號」≠「是 subclone」** — 每狀況只代表「過了哪些檢定→字面意義」，生物身分需逐一驗。
2. **對齊基礎是 paired**（混 normal reads）；真正 tumor-only 的只有狀況②與 Subclone Δβ 軸（C++ `!is_tumor continue` 源碼證實）。
3. **單一 paired run、單樣本 HCC1395**；跨樣本未做。

---

## §2 能否「清楚分類與定位」？（核心方法學點）

**裁決：能 — 但必須用互斥 partition 語言，不能用「三狀況」語言。**

- **§0 三狀況彼此重疊**（如「對齊 ∩ 一標籤多群」= 6,180）→ 問「某位點屬哪狀況」會**歧義**（一個位點可同時對齊又一標籤多群）。
- **§2 六格 / §0.7 四格是乾淨 partition** → **每位點剛好落一格、可加總、可定位**。
- 因此「逐一定位驗證」**必須先 collapse 到六格或四格**；用三狀況措辭會踩重疊歧義。

> **可定位的最強分類 = 切群總帳 + 三態決策**：先問「能不能切 ≥2 群」（80.3% 不能 / 19.7% 能），能切的再問「對齊哪條 a-priori 軸」（clear≥2 = 2,941）。**這是目前最接近「清楚定位」的框架**。

---

## §3 各狀況 ↔ clone/subclone 對應（L3 推論）

> 原則（§7 + chr2 L0）：**列全「可能是」，除非證據排除才不列；確認 ≥2 群 ≠ 確認 subclone**。下表的「subclone 相關度」是**方向性**，非已證。

### 3.1 主對應表

| 狀況 / 格 | n | clone/subclone 相關度 | 最可能身分（窮舉，未排除全列） | 驗證錨點 | 信心 |
|---|--:|---|---|---|---|
| **①對齊×單群** | 16,317 | 低-中 | germline ASM / cis-ASM / somatic residual / **true subclone(需 linkage)** / technical | 需 normal-anchored cis-control + multi-sSNV linkage | 多數非 subclone |
| **①對齊×多群** | 6,180 | 中 | 同上 + within-HP 子結構 | 同上 | 混合 |
| **②一標籤多群（tumor-only 軸）** | 8,077 | **最高方向** | **true subclone(HP 內 carrier vs germline 甲基不同)** / **epiallele** / **LEVEL artifact** / cis-ASM / LOH 偏移 | Subclone Δβ + PATTERN silhouette | 🔴 **PATTERN(354) 可信 / LEVEL(7910,97.9%) 低信心** |
| **③一群多標籤** | 510 | 低 | LOH 單親消失 / 非 HP-allele 軸驅動 / 覆蓋稀 / 真 epiallele | cis-control + entropy | 稀少、多 LOH/覆蓋 |
| **無群** | 7,483 | 極低 | 真平坦（≥20 reads）| — | 非 subclone |

### 3.2 三個 load-bearing 修正（投稿口徑，🔴 必守）

1. **狀況② 不等於 subclone 候選整包**：8,077 中大宗是 **LEVEL（mean-β 雙峰，7,910，低信心）**；嚴格 silhouette 驗的 **PATTERN 僅 354（1.16%）**；Subclone/AltSubclone Δβ-sig = 1,786。**subclone 候選須收斂到 PATTERN ∪ Δβ-sig，不可把 LEVEL 大宗算進去**。

2. **「確認 ≥2 群」85% 落 germline/carrier 軸 = 可能 cis-ASM**：clear≥2（2,941）對齊 germline/carrier **2,500（85%）** / HP **441** / 隨機 0。→ 切得很清楚 ≠ subclone；**多數清楚多群可能是 cis-ASM（變異誘發等位甲基），需 normal-anchored cis-control 才分 cis vs subclone**。

3. **subclone-specific 上界 = 5,473（KEEP 的 19.3%，30,490 的 18%）**，**仍含 cis-ASM 未排除** → 是**上界非真值**；最 load-bearing 缺口 = **G-B within-haplotype somatic-vs-baseline 控制未跑**。

### 3.3 判別性的誠實口徑（🔴 禁用詞）

| ❌ 禁用 | ✅ 正確替代 | 數據 |
|---|---|---|
| 「a-priori 在 Noise 0%」 | 「對 legacy VC Noise 非循環判別率 = **7.6%**（Strong **16.4%**，~2× 富集）」 | 0% 是 reclassify cascade 的 **tautology**（用 a-priori-sig 定義 valid）|
| 「more discriminative」 | 「**double-dip removed** / 合法 permutation null（免 collider）/ 可解釋」 | unsup-clean 6.8%/15.3% 與 a-priori **統計打平** → 優勢非偵測力 |
| 「獨有顯著性檢定」 | 「成熟 PERMANOVA/permutation 框架（Anderson 2001/McArdle/Legendre）」 | cvlr/ASMS 等都有 randomization 檢定 |

> **within_clean ≠ subclone**（前一份分析已更正）：unsupervised「乾淨多群」(4,323) 與 a-priori subclone(4,596) **Jaccard 僅 0.123**；within_clean **55% 無 anchor = epiallele**。→ 無監督切出的多群多是 epiallele，不可標 subclone。

---

## §4 各狀況的驗證設計（chr2 模板 + chr8）

### 4.1 chr2:18M = 驗證模板（已 L2 verified，可複用）

chr2:18M 是**唯一已跑通的旗艦驗證單元**，其 **9 步 pipeline** = 各狀況驗證的標準工具箱：

| 步 | 工具 | 對應「可得資訊」 | 在 chr8 是否可得 |
|:--:|---|---|---|
| 1 萃取 | pysam read×sSNV×5mC 矩陣 | tagged BAM | ✅ |
| 2 LOH context | SEQC2 CNV/LOH bed + HP imbalance + normal ALT=0 | bed ✓ | ✅ |
| 3 lineage by linkage | 00/10/01/11 互斥（0 違反=分叉）| multi-sSNV | ✅（chr8 284 個 ≥3-sSNV 簇）|
| 4 artifact guard | homopolymer context（pos4 20bp poly-T 教訓）| GRCh38 ref | ✅ |
| 5 甲基判別 | MWU + BH-FDR + **cross-basecaller**（HKU+DORADO）| MM/ML in BAM | ✅ |
| 6 🔴 **normal-HP ASM-control** | normal HP1-vs-HP2 Δβ → 排既存 germline-ASM | **normal BAM** | ✅（chr8 normal 1.4M reads）|
| 7 tumor/normal differential | all-read MWU FDR | paired | ✅ |
| 8 建可防守樹 | parsimony；nesting 才排序；**VAF 不可排 siblings** | linkage | ✅ |
| 9 誠實 tier | regional/operational；禁「confirmed N-subclone」| — | ✅ |

> **核心增值**（步 6）= chr2 把「tumor-acquired subclone 甲基」與「既存 germline-ASM 被 LOH 保留」拆開：clean 子集 {2.1,2.2,3.1,3.2,5.1,5.2}（normal 無 HP 差）vs confounded {3.3,3.4,3.5,4.1}（normal HP1-vs-HP2 顯著）。**這是 ISM 對 subclone 案例的方法學關鍵**。

### 4.2 chr8 全-LOH 想法 — 可行性 ✅ + 裁決 PARTIALLY_SOUND

**資料盤點（四項全 [verified-exists]）：**
- chr8 **~96% LOH**（139.3 Mb / 26 段；最大單段 31.68 Mb）— 遠超 chr2 單一 5.95 Mb 區
- chr8 **2,605 個 HC sSNV**（99.2% 落 LOH 內）；**284 個 ≥3-sSNV / 40kb 簇**、**10 個 ≥6-sSNV**（= chr2:18M 同等規模）
- chr8 normal tagged BAM 覆蓋齊（1,398,803 reads）→ normal-anchored control 可做
- chr8 ISM 已分類 2,094 loci（Strong 1610 / LOH-Structure 139 / …）

**🟡 PARTIALLY_SOUND（對抗裁決，4 個 confound）：**

1. **LOH 是必要鋪墊、非充分條件**：LOH collapse 掉 live germline HP1/HP2 軸（讓甲基異質性更可詮釋為 lineage）— 但**真正分開 lineage 的是 sSNV linkage，不是 LOH**。純 LOH 段（無 sSNV）的甲基分群會退化成**無監督路線（已 NEGATIVE，double-dip，勿再開）**。

2. **雙面刃（最強反證）**：chr8 正是 **LOH+HPSig FP 7.4× 富集**的同一塊地（80 個 FP 中 66 落 chr8，根因 = germline SNP 在 LOH+ASM 區被誤判 somatic）。→ 你用來建 lineage 骨幹的 sSNV linkage **本身可能含 germline-SNP 污染** = 假骨幹風險。chr2 靠 SEQC2 truth + normal 6 點全 REF 排掉；**chr8 整段未做同等 per-locus normal-REF 驗證**。

3. **retained germline-ASM confound 在 chr8 更嚴重**：chr8 是已知大型 ASM block → 「既存 germline-ASM 被 LOH 保留再被 somatic 切分」的 confound 通道**更寬**，normal-anchored control 工作量與失敗風險都比 chr2 大。

4. **單樣本天花板**：chr8 若仍只有 HCC1395，封頂同 chr2 = ⭐3；**加 chr8 不會跳到跨樣本層**（需 COLO829）。

**→ 正確用法（必要前置）：**
- **不挑「整段 chr8」**，鎖定 **1-2 個高 HC-sSNV 密度 40kb 子窗**（候選：**chr8:133.68M=7 sSNV / 144.84M=6 / 132.52M=6**）。
- 每個建骨幹的 sSNV 必須：(a) SEQC2 HC 可評估或經 normal 0% ALT + 雙股 + mapQ + 非 homopolymer 驗真；(b) **normal 該位點全 REF 排 germline 污染**（直擊 chr8 7.4× FP 根因）。
- 跑 normal-anchored ASM-control（複製 chr2 步 6）。
- cross-basecaller（HKU+DORADO）。
- 口徑限定 **characterization**（LOH-constrained, somatic-haplotag-conditioned subclonal structure），**甲基 = 被刻畫對象，非重建驅動**。

> **chr8 的真正增量價值** = **破單位點的第二 LOH 驗證點**（升 chr2 結論 tier 的關鍵前置），前提是「選對子窗 + 配 normal-control」，不是「整段染色體當單一驗證標的」。

### 4.3 「逐一驗證」不可行 → 抽樣策略

chr2 是**手做 L2 單位點**（pysam + 手標 6 sSNV + 4-agent 複核）；30,490（甚至狀況② 8,077）**不可能逐一手做**。
- **可做**：canonical（chr2:18M）+ outlier（chr8 高密度子窗、最強 silhouette 位點如 chr8 sil=0.812）+ well-explained 例的**抽樣示範**（「見樹也見林」四層）。
- **要 scale**：須先把 chr2 九步**自動化成 ISM 模組**並標準化判準（互斥閾值 / FDR / cross-basecaller / homopolymer mask / ASM-gate），**且自動化模組本身需先驗** — 這是 future work（handoff B2）。

---

## §5 推論合理性誠實裁決（§scientific-rigor §2）

**整體：PARTIALLY_SOUND — 合理，但須 5 處修正。**

| # | 原推論措辭 | 問題 | 修正 |
|---|---|---|---|
| 1 | 「清楚分類定位各狀況」 | §0 三狀況**重疊非互斥** | 改用 §2 六格 / §0.7 四格互斥 partition 語言 |
| 2 | 「狀況② 最與 subclone 相關」 | LEVEL(7910) 低信心混入；85% clear≥2 落 germline/carrier 軸 | 收斂到 **PATTERN(354) ∪ Δβ-sig(1786)**；85% 標可能 cis-ASM |
| 3 | 「可用 chr2 模板**逐一**驗證」 | chr2 手做 L2，無法 scale | 改「**抽樣** canonical/outlier + 自動化(未證)」 |
| 4 | 「chr8 **整段** LOH 也好用」 | LOH 非充分 + 7.4× FP 雙面刃 | 改「高密度 **40kb 子窗** + normal-control + 排 germline」 |
| 5 | 任何「confirmed/偵測力/0% noise/獨有檢定」 | reviewer 會重現打臉 | 單樣本只到 **regional/operational L2~⭐3**；口徑見 §3.3 |

**什麼站得住（保留）：**
- 盤點/校準層**健全自洽**（多 partition 加總 = 30,490 機械真值）。
- 文件層**誠實邊界已守**（§7「是訊號≠subclone」、PATTERN vs LEVEL、chr2 L0 禁用詞）。
- chr2 九步模板**方法學成熟**（4-agent + 第二獨立 audit byte-identical），normal-anchored ASM-control 是**真實增值**。
- 狀況②=tumor-only **方向**有 C++ 源碼錨點（非空想）。

---

## §6 下一階段「分別驗證各狀況」計劃（待用戶確認）

> 本檔到此為止只整理推論。以下為**驗證階段提案**，待確認 scope/優先序後執行。

**建議驗證單元（由可行→困難）：**

| 優先 | 驗證單元 | 對應狀況 | 方法 | 產出 tier |
|:--:|---|---|---|---|
| P1 | **chr8 高密度子窗 ×2-3**（133.68M/144.84M/132.52M）| 狀況②+① in LOH | chr2 九步（含 normal-REF 排 germline + ASM-control）| 第二 LOH 點，破單位點 → 仍 ⭐3 |
| P2 | **狀況② PATTERN 子集（354）canonical 抽樣** | 狀況② | Subclone Δβ + normal-anchored cis-control 分 cis vs subclone | 量化「真 subclone vs cis-ASM」比例上界 |
| P3 | **clear≥2 中 HP-aligned 441 vs gc-aligned 2500 對比** | 切群總帳 | normal cis-control：HP-aligned 是否較像真 subclone | 守「85% 可能 cis-ASM」界線 |
| P4 | **跨樣本（COLO829）** | 全部 | 重跑 chr8/狀況② 於第二樣本 | 破單樣本 → ⭐4（需 COLO829 normal 甲基，目前缺）|

**每單元的「可得驗證資訊」**（皆 ✅ available for HCC1395）：tagged tumor+normal BAM（含 MM/ML 5mC）、SEQC2 truth VCF + HC bed + CNV/LOH bed、DORADO cross-basecaller、GRCh38 ref、ISM perlocus 分類。

**🔴 待用戶定奪**：(a) 先驗哪個單元？(b) chr8 子窗數量（2 vs 3 vs 更多）？(c) 是否現在就接 COLO829 跨樣本（需先確認其 normal 甲基資料）？

---

## §7 三個驗證方法想法（2026-06-20 用戶提出；待設計確認後寫子任務）

> 三者皆為「**甲基 ↔ 標籤關聯**」的偵測/輸出方法，可互相做**集合比較**（containment / intersection / 差集）。本節只記錄重點與設計方向，實作細節與小步驟待 grounding 後確認、用戶 ack 後才寫子任務列表。

### 想法 1 — cluster 切分結果 → 回頭做「附近 CpG 的篩選與定位 + 軸歸屬」
- **目標**：在 sSNV 位點用「read×CpG 距離矩陣→切分→標籤對齊」拿到 cluster-level 對齊結果後，**回頭逐一檢視窗內每個 CpG**，找出哪些 CpG 符合這個切法 = **與標籤高度相關**，並標出**與哪條軸**（HP-family / HP-fine / allele）相關。
- **輸出**：per-CpG → 「與哪個 sSNV 的哪條軸相關」的篩選+定位表，供下游分析用。
- **驗證標的**：用 **BRCA2 位點**（已 verified ASM，HP-axis Δβ=−0.122）測這想法的結果是否成立。

### 想法 2 — 我方（甲基距離壓縮→切分→判讀）vs modkit 式（逐 CpG 按軸驗差異）
- **目標**：比較兩種範式 —— 我方 = read×read 距離壓縮資訊後切分；modkit/DSS 式 = 找所有 CpG，逐一按標籤軸分類測該 CpG 是否有明顯差異。
- **問題**：兩法各有何優缺點？逐 CpG 驗軸相關是否更好/更可行？
- **做法**：**小規模 pilot** 比較兩法結果差異。（對齊既有 `project_ism_vs_external_methylation_tools_comparison` Phase B；ISM 站「軸 C read-read 距離+clustering」）
- ⚠ **與已關閉方向區隔**：此非 Fine-Pairwise distance（那是 TP/FP 判別 NEGATIVE）；此為 label-association characterization，目標不同。

### 想法 3 — 直接 per-label 平均甲基 Δβ（過 β 閾值）+ 與距離切分結果做集合比較
- **目標**：另一種更簡單的方法 —— 區域內 read 甲基先平均、再在每個標籤組內平均，比較標籤組間 **Δβ 是否超過某 β 閾值**。
- **問題**：(a) 這法偵測到的結果是否更多？分佈如何？β 閾值設多少合適？(b) 是否**完全包含**「距離圖依標籤切顯著」的結果，還是有**交集/差集**？數量如何？差集有哪些案例可觀察？(c) 與想法 1、想法 2 結果差異如何？
- ⚠ 可能與既有 **SubcloneDbeta / AltSubcloneDbeta / LEVEL bimodal**（C++ 已實作）高度重疊 → grounding 須確認「新做 vs 已有」。

### 三想法的共同評估軸
- **集合關係**：三法各自的「顯著位點集合」vs 距離切分顯著集合的 Venn（containment / ∩ / 差集 + 差集案例）。
- **軸歸屬 confound**：ALLELE 軸需 germline-het null（`feedback_asm_allele_axis_baseline_confound`）；HP 軸才 somatic-controlled。
- **scope**：用戶提「小規模」→ pilot（Type A）；但「確認整個狀況/分佈」屬全 scope → 須標 partial flag 並確認是否擴全。

---

## §8 三想法實作設計（2026-06-20 grounding workflow wf_d14add83-9e2；待用戶確認後寫子任務）

> **跨想法關鍵發現**：三法**都能用既有 per-region 矩陣以 Python 跑，不需重跑 C++**（per-region `methylation/methylation.csv` + `reads/reads.tsv` + `distance/BERNOULLI/matrix.csv` + `clustering/` 全 [verified-exists]，HCC1395 TP=29,754 region / FP=627）。⚠ 30,490 universe vs 29,754 落檔矩陣（~736 位點無原始矩陣，只在上游分類）。BRCA2 pilot（5 region 22305-22309，canonical=22305 chr13:32,315,128）**零重跑可立即跑**。

### 想法 1 — per-CpG 軸歸屬（部分已實作）
- **已有**：C++ `PerCpgAsm.cpp:274-319` 已逐 CpG Fisher + BH-FDR，**但只測 HP-family 一軸 + per-CpG 向量未落檔**（只 aggregate 成 region-level 3 數）。GlobalTest 是 read-cluster 層（3 軸）非 per-CpG。
- **新做**：Python per-CpG × **三軸**（HP-family / HP-fine / allele）MWU/Kruskal + 各軸獨立 BH-FDR → `per_cpg_axis_attribution.tsv`（cpg→過FDR?+dominant_axis+屬哪sSNV+normal-ASM-confound?）。
- **可行性**：BRCA2 RegionID 22305（56 reads×229 CpG）✓；20260420 run 有 45 normal reads → normal-control 可做；真值 cross-check = `research/tsg_promoter_asm_reviewer/output/step4_ism_results.json`（Δβ=−0.122）。
- **小步驟**：S1[gate]確認 run/軸/裁決規則 → S2 join+QC → S3 HP-family 單軸對真值 → S4 三軸 → S5 normal-ASM-control 欄（複現 chr2 §A5 clean/confounded）→ S6[gate]allele germline-het null + dominant_axis → S7 報告。
- **confound**：ALLELE 軸需 germline-het null；**characterization 非 TP/FP 判別**（ASM B-discrimination NEGATIVE）；229×3 多重檢定；per-CpG↔cluster 循環性（用 a-priori 標籤不用 cluster_label）；normal hp 多 unphased(0) 可能不足分 HP。

### 想法 2 — modkit/DSS 逐 CpG vs 我方 read-distance（基礎齊、延伸 Phase B）
- **已有**：modkit v0.6.3 binary ✓（`benchmark/tools/dist_modkit_v0.6.3_26c3f9e/`，絕對路徑呼叫）+ DSS R env ✓ + 我方 PerCpgAsm 即等價逐 CpG 元件 ✓。**Phase B 已跑 region 級**（07：modkit BRCA2 −0.159 vs ISM −0.122 收斂）；816-loci 已有 Python proxy（11：結構顯著但率測弱 8.5%）。retag BAM `runs/HCC1395_somatic_hp.bam` ✓。
- **新做**：**同位點同軸的逐-CpG 集合 head-to-head**（07 只到 region 級）→ 位點級 2×2 concordance（結構顯著 × 有≥1顯著CpG）+ CpG 級 Jaccard/Venn + discordant 案例熱圖。
- **小步驟**：S1[gate]位點集+軸+粒度 → S2 modkit per-CpG → S3 DSS per-CpG → S4 我方 PERMANOVA+PerCpgAsm → S5[gate]三表 join 成 concordance → S6 discordant 熱圖+報告。
- **confound**：軸需先 retag（modkit 原生只 HP1vsHP2）；**壓縮維度本質差**（marginal per-CpG vs haplotype-coherence，非 bug）；閾值敏感需 sweep；**與 Fine-Pairwise NEGATIVE 嚴格區隔**（那是 TP/FP 判別 DEAD，此為 label-association 描述對照）；三法皆無 cis-control。

### 想法 3 — per-label Δβ + 集合比較（C++ 大部分已實作）
- **已有**：**想法3 字面即 `feat/summary-nreadsvalid:RegionProcessor.cpp:2018 compute_group_dbeta_test()`**（per-read mean β→組平均→Δβ）；顯著性=**permutation test（999 洗牌 p≤0.05）+ min-group≥3，無固定 Δβ 大小閾值**；5 軸已接（germline_asm / subclone_hp1/hp2 / alt_subclone_hp1/hp2）；per-CpG driver 有硬編 |Δβ|>0.2；LEVEL bimodal（gap>0.15）已實作。Python 原型 `dbeta_combo_pilot.py` ✓。
- **新做**：**Δβ 閾值 sweep**（純後處理 τ∈{0.10,0.15,0.20,0.25}）+ **集合論 Venn**（set A=Δβ法命中 vs set B=距離依標籤顯著）+ 差集案例挑選 + confound 量化（driver_cpg_n 對照 LEVEL-shift；allele null）。
- **資料**：可 Python 從既有矩陣算 Δβ（不需 C++ 重跑）；若要原生 dbeta 欄全基因組/6 樣本則需 feat 分支重跑（長計算，非必需）。
- **小步驟**：S1[gate]set B 定義（22,497 對齊 vs 2,941 clear≥2）+ scope → S2/S3 取/算 Δβ → S4 Venn+sweep → S5[verify]加總自洽 → S6 差集案例 → S7[gate]選閾值 → S8 confound → S9 報告。
- **confound**：region-mean 洗掉 per-CpG 結構（LEVEL artifact，用 driver_cpg_n 對照）；ALLELE baseline confound；**與 tumor-only 非監督 NEGATIVE 區隔**（a-priori conditioned 非無監督）；閾值非中性（ONT β 雜訊大，0.1 可能假陽偏高）；「偵測更多」是方向性（Δβ 對 same-level-different-pattern 天生盲 = 差集賣點）。

### 三想法統一比較框架
| Evaluation | 比誰 | 集合定義 | 指標 |
|---|---|---|---|
| (a) 位點層 Venn | 想法3 vs 距離切分 | region→{顯著/否}（significance_summary.csv 直接 join，秒級）| Jaccard + 4-cell |
| (b) CpG 層 | 想法1 vs 想法2 | (region,CpG)→{與軸相關/否}（讀 methylation.csv）| per-CpG Jaccard + Δβ 同號率 |
| (c) 對齊軸一致 | 三法各自 dominant 軸 | region→對齊哪條 a-priori 軸（CramérV≥0.7）| 軸一致率 + Cohen's κ |

### 共通鐵則
1. **scope**：pilot（5 BRCA2 region，零重跑）vs full-scope（讀 29,754 矩陣 + 集合運算，**仍非重跑 C++**，I/O-bound）。
2. **三法皆 characterization 非判別** — 全程禁「TP/FP 判別 / filter / more discriminative / 0% noise / 獨有檢定」。
3. **ALLELE 軸**一律配 germline-het null；**HP 軸**才 somatic-controlled。
4. §13.0：先落 .tsv→Read→才寫報告；撰寫與分析不同批。

---

## §9 三方法 ↔ clone/subclone 驗證（2026-06-20 workflow wf_6ac241ae-a79；對抗裁決）

> **裁決 = CHARACTERIZATION_ONLY**：三方法對 subclone「**存在**」的獨立驗證淨貢獻 ≈ 0；它們是「刻畫已由基因型定義之 subclone 的甲基」的工具。對外口徑：**「characterize methylation of genotype-defined subclones」，不可寫「verify/discover subclones」**。

### 9.1 為什麼（循環性）
三方法**都 conditioned on a-priori 標籤軸**（haplotag germline-tag/carrier-tag、ALT/REF）—— 軸是**輸入不是輸出**。它們回答「哪 CpG/read 跟著已給定的軸走」，**不檢驗該軸是否=真 subclone**。若軸本身是 cis-ASM/germline，方法照樣輸出「符合」（**無否證力**）。對齊既有定論：tumor-only 非監督 NEGATIVE（double-dip）、確認≥2群≠subclone（85% 對齊 germline/carrier=可能 cis-ASM）、within_clean≠subclone（Jaccard 0.123）。

### 9.2 外部論文定位（grounding）
- **subclone-methylation 的「確立」只在 single-cell 有先例**：Gaiti2019（SF3B1 subclone 精準映射純-epimutation clade，Fisher P=7.4e-9）/ Epiclomal / MethylTree / EPI-Clone —— 它們有**細胞身分**。
- **bulk read-level 甲基 subclone 驗證 = GAP（未建立）**：bulk 平均遮蔽 subpopulation；跨 phase-block 無 read/cell 連結（MethPhaser PMID38909018）。我方三方法落 bulk regime → 對應 germline-ASM 發現工具（cvlr/ASMS/CpelNano）或 supervised deconv（MethylBERT），**非 subclone 驗證**。
- **白地**：「bulk 甲基 + 基因型 anchor 聯合驗 subclone」**無 established 先例** → 機會，但須 framing 為 characterization。
- 🔴 **LOH-unmask-ASM**（Martin-Trujillo PMID28883545：tumor imprinted-DMR 甲基變化 **82–91% 由 CN/LOH 解釋非 epimutation**）+ cis-ASM（Do2020 罕見）= bulk 把既存 ASM 誤認 subclone 的主通道；外部控法 = normal cis-control + imprinting blocklist（Rosenski PMID40069157）。

### 9.3 實際位點接法
| 位點 | 有 sSNV linkage？ | 三方法角色 | 裁決 |
|---|---|---|---|
| **chr2:18M** | ✅ 6 sSNV | 想法1 = **系統化 §A5 手做的 clean/confounded 拆解**；但承重骨幹（00/10/01/11 互斥 linkage）來自 pysam 非三方法 | characterization 層**疊在**既有 genetic anchor 上 = 「方法有幫」的範例 |
| **BRCA2** | ❌ 單 TSG 無 linkage | 想法1 抓 18 carrier-軸 CpG（方向符 −0.122）；但無 linkage + normal unphased | **只能 ASM characterization，不能驗 subclone** |
| **chr8 純-LOH 段** | ❌（純 LOH 無 sSNV）| 想法3 差集 A−B（chr8:118671269 LOH+無driver）**正是 LOH-unmask-ASM 偽 subclone** | **放大 confound 非助力**；助力僅在高密度 40kb 子窗 + genetic anchor + normal-control |

### 9.4 承重柱三方法給不出
| verification anchor | 來源 | 三方法產出？ |
|---|---|---|
| sSNV linkage（互斥/巢狀）| pysam read×sSNV | ❌ 只吃 read×CpG + 標籤 |
| CCF 梯度/lineage 排序 | nested genotype | ❌（VAF 不可排 siblings）|
| normal-anchored cis/ASM-control | normal-HP Δβ | ⚠ 概念含但**實測 normal unphased 跑不了** |

### 9.5 但是 — 建設性淨貢獻（在 pipeline 內的位置）
雖**非 verifier**，三方法在「先建 genetic anchor → 再刻畫」的 pipeline 中有用：
1. **characterization 層**（anchor 下游）：刻畫已驗證 subclone 的 epigenotype（論文 Ch4 chr2 demo）。
2. **triage/優先序**：想法1 normal-control-clean 子集 = 候選真 subclone-甲基；想法3 差集 A−B = LOH-unmask confound flag → **把昂貴的 genetic-anchor 驗證集中到有希望的位點**。
3. **系統化 chr2 手做 ASM-control**（想法1 per-CpG clean/confounded，**前提 normal phasing 可得**）。
4. **白地新穎性**：bulk 甲基 + 基因型 anchor 聯合 characterization 無先例 → 可為新框架的 bulk 刻畫側（framing=characterization）。

### 9.6 要真推進 subclone 驗證須先補
- **#1 blocker = normal phasing**（解 unphased → 啟用 clean/confounded cis-control）。
- genetic anchor：chr2:18M（已有）systematize；chr8 高密度子窗先建 anchor 再刻畫。
- 獨立確立：single-cell 或跨樣本外部真值（Fang2021/COLO829）。

---

## §10 想法1 兩角色定位 × 論文題目對齊（2026-06-21 方向建議；workflow wf_34d1adba-f62）

> **方向建議**（非已驗證結論）：把想法1（per-CpG ASM 偵測/軸歸屬）定位成**兩個角色**，並評估其與論文題目 *Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing* 的貼合度。

### 10.1 題目本意（A-framing，thesis spec L19/L38 SoT）
- **reconstruction = 動詞，承重者 = somatic haplotagging 骨幹**（sSNV linkage/CCF/lineage 由 pysam read×sSNV）。
- **"using methylation profiles" = 被 characterize 的有界附加層，非 reconstruction driver**。
- 🔴 **題目 integrity 條件**：摘要首段必須立即點明分工，否則「reconstruction using methylation profiles」字面會被讀成「甲基化重建了亞克隆」（spec D2 候選 C = 高 overclaim 風險，hedge 全靠摘要+動詞紀律外掛）。

### 10.2 兩角色 × 對齊裁決
| 角色 | 對題目 | 定位 / 放哪 |
|---|---|---|
| **角色1 = subclone-甲基宣稱的 confound 分離引擎** | 🟢 **STRONG_FIT** | 讓「using methylation profiles」在 subclone 題目裡 legitimate 的**承重柱**；= carrying-claim 的 **de-confound floor（robust/fallback-safe）**。放 Ch4 chr2 demo methods-foundation + Ch3 方法。疊在 genetic anchor 下游，非 headline driver |
| **角色2 = germline/cis-ASM landscape 刻畫** | 🟡 **ORTHOGONAL** | established（cvlr/ASMS 線）但 **germline ASM 本質非 subclonal** → methods-foundation / 獨立 secondary aim，放 Ch3 + Ch5。**標籤必須「ASM landscape」非「subclone evidence」**（否則踩「確認≥2群≠subclone」紅線） |

### 10.3 鎖定 framing（對抗裁決定稿）
> 「以 somatic haplotagging 重建 subclone 骨幹；在該骨幹上用 per-read 甲基化做兩件有界 characterization —— (角色2) 刻畫 germline/cis-ASM landscape，並以此 (角色1) 把『subclone-acquired 甲基』與既存 ASM 分離，使任何 subclone-甲基陳述都是扣除 confound 後、conditioned on 基因型 anchor 的誠實刻畫，而非甲基化獨立偵測或重建 subclone。」

### 10.4 4 守法條件
1. **動詞紀律**：甲基 = characterize/corroborate/quantify-bounded；**禁** enable/improve/reconstruct（指甲基時）；reconstruction 只描述 haplotagging 骨幹。
2. **摘要首段點明分工**（題目 integrity 條件）。
3. **角色2 標「ASM landscape」非 subclone 證據**；混入即 overclaim；allele 軸需 germline-het null。
4. **圖文標清「甲基層疊在 genetic anchor 上」**（chr2 骨幹 00/10/01/11 linkage 來自 pysam 非甲基）。

### 10.5 2 個誠實 caveat（升級/投稿前）
- 🔴 **角色1 目前只能做半邊**：confound 分離的「扣 germline-ASM」需 **phased normal**，pilot 實測 **normal 全 unphased**（BRCA2 5/5；chr2 可做、BRCA2 跑不了）→ 對外只能寫「**在 normal-phasable 子集示範**」，**不可寫「已系統化扣除 ASM」**。**#1 blocker = normal phasing**。
- **量化數字多 L3**（判別率/Jaccard 0.123/clean-confounded 計數在 memory + /tmp JSON 未 committed）→ 投稿前須從 /tmp 重算固化。

### 10.6 方向建議（收斂）
- **押角色1（de-confound floor）為甲基層的承重論述** —— robust、fallback-safe、直接服務題目。
- **角色2 當 methods foundation / 獨立 aim**，標籤分開、別混進 subclone headline。
- **先解 #1 blocker（normal phasing）** 才能讓角色1 從半邊變完整；可先在 **chr2 normal-phasable 子集**實跑 germline/cis/carrier 三軸分離出實際數字（驗證角色1 可行性）。

---

## §11 🔴 散度（PERMDISP）未計算 — 「標籤顯著」過於樂觀（2026-06-21 實證）

> **觸發**：用戶質疑「標籤顯著是否用散度驗證、是不是過於樂觀」。查證 = **20260420 run 根本沒算散度**。

### 11.1 鐵證：散度欄是 stub
- `LabelHPDispersionP / LabelAlleleDispersionP / ClusterDispersionP` 在全 29,754 region **全 == 1.0**（DispersionWarn 全 False）→ 統計上不可能是真跑（真 PERMDISP 給分佈）= **stub 預設值（沒跑）**。
- 原因：**20260420 run 用舊 binary，早於 `check_dispersion` 實作**（C++ `src/core/StructureTest.cpp:262` 已有，但該 run 在其前）。**非 C++ 缺功能**。
- ⚠ **之前已 flag**：`docs/methodology/20260616_..._false_negative_audit_01.md:164`「dispersion 欄未填 → 離散度假象未排除」= 同一問題 06-16 已記。

### 11.2 後處理實證（skbio PERMDISP on 既有距離矩陣；PERMANOVA 重現 CSV→對齊可信）
| Region | HP軸 PERMANOVA(重現CSV) | HP軸 PERMDISP | HP 判讀 | allele PERMANOVA | allele PERMDISP | allele 判讀 |
|---|--:|--:|---|--:|--:|---|
| 22305 | 0.015(CSV0.01) | **0.03** | 🔴 confounded | 0.005 | 0.16 | 🟢 clean |
| 22306 | 0.265(CSV0.34) | 0.07 | n.s. | 0.005 | 0.905 | 🟢 clean |
| 22307 | 0.045(CSV0.05) | **0.01** | 🔴 confounded | 0.365 | 0.75 | n.s. |
| 22308 | —(HP2<3) | — | — | 0.035 | 0.945 | 🟢 clean |
| 22309 | —(HP2<3) | — | — | 0.785 | 0.16 | n.s. |

### 11.3 關鍵：dispersion confounding 是**軸特異**的
- **HP 軸（germline）**：可測 2/3（22305/22307）**dispersion-confounded** → 連旗艦 BRCA2 22305（Δβ=−0.122）的 HP-ASM 都中招 → **germline-HP「標籤顯著」過於樂觀**。
- **allele 軸（somatic/ALT-REF）**：3/3 可測 **location-clean**（PERMDISP 0.16/0.905/0.945）→ **somatic 軸是真群心位移、沒被散度污染**（對 subclone 主軸是好消息）。
- → **「標籤顯著 91.8%」必須分軸 + 控散度重評**；單一百分比無意義。前述 2×2（§討論用）的「標籤」軸標註「**未控 dispersion，HP 軸 over-optimistic**」。

### 11.4 誠實 caveat
- 199 perm 較粗 + n=3 可測 HP 極小 = pilot 級（機制真實、magnitude 待全 scope）。
- 🔴 **BERNOULLI 非歐**（PCoA 負特徵值 −0.216≈正 +0.613）→ **PERMDISP/betadisper 在非歐空間是近似**，**C++ analytic-F 與 skbio permutation 同樣中招** → 「非歐距離的正確散度檢定」= 獨立方法學待辦（Cailliez/Lingoes 校正或距離原生法）。
- committed「Dispersion 1,700(5.6%)」vs L3「72%」vs pilot（HP 軸 2/3）三者不一致 → magnitude 未定，須**全 scope 分軸**定論。

### 11.5 路線決策
- **不修 C++**（已有 check_dispersion）；**觀察先（後處理）**最快又最正確（同距離矩陣 + permutation/analytic）。
- **全 scope 分軸 PERMDISP**（✅ 完成，見 §11.6）：用 **analytic betadisper**（ANOVA-F on 到群心距離，免 permutation 迴圈才跑得完）對 HP/allele 軸各算 location-clean vs dispersion-confounded 比例。
- 數據 SoT：`_assets/data/permdisp_check_brca2.json` + `fullscope_permdisp_summary.json` + `fullscope_permdisp.tsv`。

### 11.6 ✅ 全 scope 分軸 PERMDISP 定論（2026-06-21；27,304 label-sig region，analytic betadisper）
| 軸 | loc-sig | 可測(group≥3) | **location-clean(真差異)** | **dispersion-confounded(假樂觀)** |
|---|--:|--:|--:|--:|
| **HP（germline）** | 11,959 | 8,360 | **5,482 (65.6%)** | **2,878 (34.4%)** |
| **allele（somatic）** | 26,720 | 21,990 | **14,049 (63.9%)** | **7,941 (36.1%)** |

- **🎯 magnitude 定論**：「標籤顯著」中 **~34-36% 是 dispersion-confounded（過於樂觀）/ ~64% 是真 location 差異**。→ over-optimism **約 1/3，非 72%**。三方落差釐清：committed「1,700(5.6%)」= 更嚴 REMOVE 組合 gate；L3「72%」= **高估**；BRCA2 pilot「HP 2/3」= 小樣本巧合。
- 🔴 **修正 §11.3 的軸特異說法**：BRCA2 pilot 看似「HP 髒 / allele 乾淨」，**全 scope 推翻** —— **兩軸 confound 率相近（HP 34.4% / allele 36.1%）**。「allele 乾淨對 subclone 是好消息」**不成立（n=3 巧合）**，誠實降級。
- ⚠ caveat：(1) analytic ANOVA-F（非 permutation）→ 近似，但對 BRCA2 與 skbio binary 判讀一致（22305/22306/22307）；(2) **BERNOULLI 非歐** → PERMDISP 本身近似（C++ analytic-F 與 skbio 同病）；(3) testable < loc-sig（部分 region group<3 不可測）。
- **對下游意義**：任何用「label-PERMANOVA 顯著」的結論（含 2×2 §討論、想法3 全 scope Venn 的 set B）**須扣掉 ~1/3 dispersion-confounded**；真正「有 location 差異」的 label-sig ≈ **64% of testable**。投稿任何 ASM/結構顯著數字**必附 PERMDISP location-clean 比例**。

---

## Provenance
- 盤點數字：`docs/methodology/20260617_structure_label_situation_inventory_01.md`（verify_breakdown 25/25 PASS）+ `_assets/20260618_subcluster_pilot/{split_accounting,spectrum_decision}.json`（本輪 workflow 加總驗算 PASS）。
- chr2 模板：`docs/paper_focus/02_paper_framework/20260615_chr2_subclone_case_and_method_concept_01.md`（verdict_02 定案）。
- chr8 可行性：`/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed`（awk）+ `significance_summary.csv` + `/tmp/perlocus_flags.csv`（⚠ 暫存未固化，引用前須重算落檔）。
- 判別率 7.6%/16.4%、Jaccard 0.123 = **L3**（memory + /tmp JSON，未 committed）；投稿前從 `/tmp/{apriori_vs_unsup,set_final}.json` 重算固化。
- 對抗裁決：workflow `wf_3b296dc8-143`（chr8_soundness + overall_reasonableness 皆 PARTIALLY_SOUND）。
