<!--
建立時間: 2026-06-12
報告類型: 中文碩論「章節↔真實脊柱↔素材↔待補」對齊表（取代架構 spec 的 Results 結構）
任務類型: D handoff — 對齊既有 paper_focus 框架後的中文碩論撰寫藍圖
狀態: alignment table / 待用戶 confirm（決策：標題保留 + de-confound 防彈脊柱 + phasing-assist 移 Discussion）
data_sources:
  - docs/paper_focus/00_共識證據台帳_20260612_01.md (45 條硬共識 + 13 矛盾 + 4 校正)
  - docs/paper_focus/02_paper_framework/論文架構_正式學術版_Slide2Thesis格式.md (v2 Results 脊柱)
  - docs/paper_focus/00_論文章節材料對應總表_20260612_01.md (master 對應表 §0/§A)
  - docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md (phasing-assist→Discussion future)
provenance_note: 所有數字引自共識台帳 §1（45 條 grep-confirmed）或 v2 架構（🟢 T-PROV）；已套用 06-12 稽核 4 校正；本檔不產新數字。撰寫每節前回查對應 primary。
-->

# 中文碩論章節對齊表 — Subclonal reconstruction（de-confound 防彈脊柱版）

> **取代**：架構 spec（`20260612_thesis_architecture_spec_01.md`）的 §2 章節結構 + §3 Results 脊柱（那版建在 A 面 phasing-assist，已校正）。spec 的 Slide2Thesis 轉換法、10 紅線、格式慣例仍有效。
> **本檔職責**：把中文碩論六章對齊到「已對抗稽核的真實脊柱（de-confound floor）」，供逐章撰寫。

---

## §0 兩決策落地（2026-06-12 用戶拍板）

1. **標題保留**：維持 **「Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing」**。
2. **🎯 主軸（用戶親述構想後定，2026-06-13 update）= 個案驅動 subclonal reconstruction（proof-of-concept）**：
   - 「**ISM 程式能找到**『甲基差異 + 多點 somatic 突變可**共同輔助區分/驗證** subclone 與演化』的位點（exemplars）」= 論文 climax/主軸。甲基 = **輔助/corroborate**（與多點突變並用），**非單獨判別器**。
   - **de-confound catalog = 脈絡**（全基因組甲基分群多 germline、乾淨 somatic-cis 罕見 → **凸顯個案的特別**）。Rule(de-confound 脈絡)+Exception(exemplar 重建 climax)。
   - **exemplar**：chr2:18M 區=主重建 demo（Fig1，假想樹）/ chr17/TBC1D16=乾淨 cis / BRCA2=illustrative（copy-confound caveat）。
   - phasing-assist V1-V12 → 第五章 future tooling；不押 HD-1。
3. **🔴 誠實護欄（用戶接受，全文強制）**：假想樹標 **hypothetical/consistent-with-data**；甲基寫 **corroborate 非 validated**；大規模驗證+判準=**future**；single-sample/single-pipeline **⭐3** 明標；subclone-甲基 somatic-specificity=**undetermined**。動詞：甲基 characterize/corroborate；reconstruction 描述「個案展示的途徑」非已驗證主張；禁「甲基已驗證重建亞克隆」。

**4 處稽核校正（全文強制，取代 stale 值）**：
| # | ✅ 正確口徑 |
|---|---|
| 1 | same-hap NG=2 Inner = **occupancy 6/6 ≥93%**（0.932/0.990/0.988/0.965/0.983/0.970）；標明度量=same-hap occupancy（非 TP-rate）|
| 2 | 乾淨 somatic cis = **HCC1395 816 可測 HP-axis loci 中僅 1**（chr17, ~0.12%）；**勿用 1/332,705**（332,705 未全 cis-tested）|
| 3 | chr17/TBC1D16 = **最強 cis 候選（copy 軸已控）非「乾淨」**：殘餘 allele-axis d=0.183>d_within、未過 Cochran gate、單樣本 n=30 Bonferroni-marginal |
| 4 | subclone-甲基 somatic-specificity = **undetermined**（押 G-B within-hap null，未跑）；germline-het null 是錯對照、p=0.922 是 fail-to-reject。**勿寫成已 NEGATIVE 也勿寫成已 POSITIVE** |

---

## §1 中文碩論六章對齊總表

| 章 | 對應 v2 | 核心 | 主要素材（共識台帳 §1）| tier |
|---|---|---|---|---|
| **第一章 緒論** | Intro | long-read 單分子同讀 somatic+haplotype+5mC → read-level 問「甲基化加什麼、在哪失效」；hook=de-confound floor + rule/exception；reconstruction 歸骨幹 | framing + foundation §1 | 🔵 |
| **第二章 文獻探討** | Intro 1.2 + background | 業界甲基 haplotype 分析=軸A per-CpG 率差（germline 軸）；ISM 站軸C read-read 距離+PERMANOVA+normal-cis；癌症 ASM（Do2020/Martin-Trujillo）；epipolymorphism（Landan2012）| method_comparison 02/03/12 + KB 11 | 🔵 + L3 |
| **第三章 材料與方法** | Methods §4 | ISM read-level pipeline（read×CpG→距離→reliability-gated 分群→HP/LOH tag→normal-cis）+ catalog 7-TAG + 9 驗證機制 + FDR 校準 | 共識台帳 §1.4/§1.6；既有 ch3_methods.md 大半可重用 | 🟢 |
| **第四章 結果** | Results §2 | de-confound 防彈脊柱 R1-R7（見 §2 細表）| 共識台帳 §1.1-1.5 全 🟢 | 🟢/🟡 |
| **第五章 討論** | Discussion §3 | 規則意義（de-confound 貢獻）+ 例外（chr17）+ 邊界機制 + **future tooling（phasing-assist V1-V12）** + limitations | §3.1-3.4 + VERIFIED_RESULTS（future）| 🔵+🟢 |
| **第六章 結論** | Conclusion | 三訊號 same-pipeline cross-check 同一 haplotype 重塑；甲基能 co-validate、不能 filter；交付 catalog + 負結果目錄 + tooling | framing | 🔵 |

---

## §2 第四章 結果 — de-confound 防彈脊柱（R1-R7）

> **排序（build-up to exemplar climax）**：先建立 de-confound 脈絡（R1-R4：ASM 存在→不能 filter→catalog RULE→copy 非 driver）→ **climax 主軸 R5：個案展示 ISM 找到「甲基+多點突變共同輔助區分/驗證 subclone」的位點（chr2 demo + chr17 cis）** → R6 LOH phasing 骨幹 → R7 跨樣本（個案 sample-specific）。每節一圖、先報不解釋；R5 為論文亮點。

| # | 節題（會說話）| 真值（共識台帳 §1 / v2，🟢 除非標註）| 圖 | 待補 |
|---|---|---|---|---|
| **R1** | ASM 存在且 somatic-enriched | TP significant **3.95% > FP 1.07%**（subhap-matched 3.77 vs 1.09 ~3.5×）；但 ~4% 低 sensitivity【🟢 ledger:85】| fig2 ✅ | — |
| **R2（最硬）** | 甲基**不能** filter somatic TP/FP（四道）| (a) LOSO 100% circularity（+0.02236→−0.00012, mean **−0.00004 p=0.125**）；(b) \|Δβ\| **AUC=0.505** 落 null；(c) ~4% sens + COLO829 TP≈FP（0.0089≈0.0103）；(d) regression-to-extreme OR **8.63/4.09/5.84**【🟢 全 grep】| fig4 ✅ | OR provenance tier 標註 |
| **R3 ⭐ RULE** | read-level catalog：甲基分群是 germline-allelic 非 somatic | 332,705 位點→**7 TAG**：reliable **12,868（3.87%, TAG-C）全 germline-allelic**；latent 28,254（TAG-E，enrich≈base 無判別）；none 291,518（87.6%, TAG-G）；**乾淨 somatic cis = 816 可測中僅 1（chr17）**【🟢 catalog 6/6 驗，⚠校正#2 分母】| fig3（P5 per-TAG）✅ 主圖 | 口徑註：12868=TAG-C / 11689=reliable(TP+FP) / 12876=all（矛盾#3）|
| **R4** | copy-dosage **非** driver；cis vs copy partition | MW neutral-vs-gain **p=0.6183**；signed ρ(\|Δβ\|,CN)=**−0.0829** → REFUTED；BRCA2=copy-confounded（d_within −0.023 ≪ d_copy −0.11）【🟢】| fig6 ✅ | — |
| **R5 ⭐主軸 climax**<br>EXEMPLAR 重建 | 個案展示：ISM 找到「甲基狀態 + 多點 somatic 突變**共同輔助區分/驗證** subclone 與演化」的位點（proof-of-concept）| **chr2:18M 區=主 demo**（HCC1395，假想樹 2-1/2-1-1/2-2-1/2-2-2/2-2-3；ClairS+DeepVariant 一致）⚠圖數字未驗證、寫前溯源；**chr17/TBC1D16=乾淨 cis**（within 0.142≫copy 0.029, perm p=0.001【🟢】，⚠校正#3「最強候選非乾淨」）；**BRCA2=illustrative**（標 copy-confound）| chr2 假想樹圖（Fig1 候選）+ fig5 | 🔴 護欄：假想樹標 hypothetical、甲基 corroborate、大規模驗證 future；FDR #24 |
| **R6 脊柱（reconstruction 元素）** | LOH-constrained phasing 軸：子克隆重塑 | NG=2 Inner same-HP1 > Outer，**方向 7/7（W=28, p=0.0078125）**；same-hap **occupancy 6/6 ≥93%**（⚠校正#1）；paired 負控 NS；AF→NGroups=**phasing 非甲基（HD-4, r=0.656）**【🟢】| fig5 ✅ | by-construction caveat（HD-1 R-SELFREF 未跑）；single-pipeline ⭐3；**不押 HD-1，寫 Grade B+ 支撐** |
| **R7** | 跨樣本復現（現象層）| **6/6 excess-over-null +0.101–0.241（mean 0.168，3 癌種）**；同-locus **0/38 private**（UNDERPOWERED E[overlap]=0.16）【🟢】| —（待製彙整圖）| 6/6 在 *_gwasm.json 對賬 |

> 🔴 **subclone-甲基 somatic-specificity 不在 R 系列當已定案**：standalone「甲基能/不能拆 subclone」一律標 **undetermined（押 G-B null）**，放第五章討論為 open-question，**不寫 Methods-Neg-4 為已 NEGATIVE**（校正#4）。

---

## §3 既有 docs/thesis 草稿處置

| 檔案 | 處置 |
|---|---|
| `chapters/ch3_methods.md` | ✅ **大半保留**；補：catalog 7-TAG schema（§4.3）、9 驗證機制清單、FDR #24 校準、ISM reliability-gate 三態、與 epipolymorphism 區隔（O11 AUC 0.845→0.530）。phasing-assist 方法（3.5）部分移述為 Discussion future 的方法。|
| `chapters/ch4_results.md` | 🔁 **重構**為本 §2 的 R1-R7 de-confound 脊柱；現有 V1-V12 內容（V10/V6/V11）移到第五章 Discussion future tooling。|
| `chapters/00_abstract.md` | 🔁 **重寫**：headline 改 de-confound floor（甲基 germline 非 somatic、不能 filter、catalog）；reconstruction 歸骨幹；subclude 改 undetermined。|
| `chapters/ch1_introduction.md` | 🔁 **重寫 hook**：rule（de-confound）+ exception（chr17 + phasing co-validation）+ 誠實 catalog；gap 三項（read-level de-confound 未做 / 甲基-as-filter 未公開證偽 / somatic-controlled ASM 未界定）。|
| `20260612_thesis_architecture_spec_01.md` | 📌 §2/§3 被本對齊表取代；其餘（Slide2Thesis 法、10 紅線、格式）保留。|

---

## §4 待補 checklist（影響哪節強度）

| 項 | 內容 | 影響 | 狀態 |
|---|---|---|---|
| **G-A** | V10/V11c + catalog 跨 6 樣本（normal 甲基 5/6 ready）| R3/R7 ⭐3→⭐4 | P0 未跑 |
| **G-B** 🔴 | within-haplotype somatic-vs-baseline 正確 null | subclone undetermined→定案 | P0 未跑（reopen 前置）|
| **HD-1** | R-SELFREF（~25-50hr C++）破 by-construction 循環 | R6 grade（不押→保守可不跑）| 暫緩 |
| **FDR #24** 🔴 | reliable 12,868 + chr17 perm p=0.001 跨 816 須 null/FDR 校準（含 n_reads 校正，O11 教訓）| R3/R5 嚴謹度 | 投稿前必補 |
| **Fig 1** | rule+exception 整合 schematic | Intro/Methods 首圖 | 待製 |
| **catalog 口徑註** | 12868/11689/12876 三值口徑（矛盾#3）、none 291,518（矛盾#6）| R3 正確性 | 撰寫時標 |
| **citation** | `{{cite:…}}` 過 /citation-verification（守 3 不對齊）| 全文 | 投稿前 |
| **Table 1** | 6 樣本統計（purity/depth/ploidy）| Methods | {{待填}} |
| **軟體版本 / Data&Code Availability** | longphase-S/ClairS-TO/ISM commit + repo + deposition | Methods | {{待填}} |

---

## §5 誠實紅線 + 校正速查（撰寫每節對照）

**6 紅線**（master 表 L1）：① 甲基=germline-haplotype 層級；② 甲基絕非 filter（死四道）；③ T3 usability NEGATIVE（可分≠可救）；④ T2 勿宣稱歸 H3；⑤ BRCA2=copy-confounded 非乾淨 cis 錨（chr17 最強候選非「乾淨」）；⑥ cohesion≠cis、single-pipeline ⭐3。

**4 校正**：same-hap 6/6 occupancy ｜ cis 1/816（非 1/332,705）｜ chr17 最強候選非乾淨 ｜ subclone undetermined。

**framing 守則**：三訊號用「**same-pipeline cross-check**」非「orthogonal co-validation」（矛盾#10）；reconstruction 只描述 haplotagging 骨幹；甲基一律 characterize。
