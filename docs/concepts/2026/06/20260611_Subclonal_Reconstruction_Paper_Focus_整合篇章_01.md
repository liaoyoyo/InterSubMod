# 重點論文聚焦篇章（補充面）— Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing

> **⚠ 本檔不是新主軸 SoT**。新主軸 SoT = `InterSubMod/docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md`（另一 session 2026-06-11 建立，聚焦**甲基-phasing-assist 線 V1-V12** + 6 樣本資產盤點 + G-A~G-E gap）。
> **本檔職責 = 補充另一個互補面**：ASM-characterization + **四道 NEGATIVE methods** + LOH-phasing 脊柱（即 foundation 文件所稱「G6/G1 降為支撐材料」的詳細證據映射）。數字 consolidate 自 `InterSubMod/knowledge/11_external_literature/10_paper_readiness_convergence.md`（7-thread 稽核 2026-06-08）+ 本 session ASM 工作站/驗證/SEQC2-CN（commit 8dbb931→5e69a99）。
> **生成 2026-06-11｜Tier 4 paper-scope｜誠實評估，非歡呼就緒。**

## 🔗 兩研究線 reconcile（同一論文，兩個面）

| 面 | SoT 文件 | 涵蓋 | 對論文 |
|---|---|---|---|
| **A. 甲基-phasing-assist**（foundation）| `..._paper_foundation_01.md` | V1-V12：甲基幫忙定相/救 unphase read、T1/T2/T3、subclone 甲基存在性 V11c、6 樣本資產 | 甲基**能為重建加什麼**（assist + subclone 甲基存在）|
| **B. ASM-characterization + NEGATIVE**（本檔）| 本檔 + convergence | LOH-phasing 脊柱(T1)、四道 filter-NEGATIVE、copy-confound、chr17 cis、BRCA2 退役、HD-1 | 甲基**不能做什麼**（誠實目錄）+ 重建脊柱 grade |

**🔴 兩線共享的一個開放問題（必先解，gate 論文甲基貢獻）**：subclone-discriminating 甲基（foundation V11c「亞群 ALT vs 母本 REF 可分」）到底是 **somatic-specific** 還是 **germline-allelic baseline**？
- foundation 框：存在性窄翻（V11c ⭐3 核心候選），但 §5 **G-B 正確 null 未跑**（承認 germline-het null 是錯對照）。
- convergence 框：用 blind-ARI/germline-het null → **germline-allelic 非 somatic-specific**（Methods-Neg-4）。
- **honest 收斂**：= 同一個未決問題（foundation G-B）。**論文甲基-subclone 貢獻的強度押在 G-B 對照**（within-haplotype somatic-vs-baseline，非跨-hap germline-het）。**未跑 G-B 前，甲基-subclone 寫「存在性窄訊號 + 可用性負」，勿寫「甲基可重建子克隆」。**

---

## TL;DR（一句話）

**這篇論文今天就能開寫**：用 ONT somatic haplotagging 重建 **LOH-constrained 的 subclonal haplotype 結構**（genetic 引擎，Grade B+ ⭐3），並以 read-level 甲基化 **characterize 子克隆內的 allele-specific methylation + 交付一份機制解明的「甲基能做/不能做什麼」誠實目錄**。**唯一決定性卡關 = HD-1（R-SELFREF 循環對照跑或不跑），由你決定。**

---

## 1. 標題拆解 → 誠實證據映射（核心：不能把甲基當重建驅動力）

| 標題成分 | 是什麼 | 既有證據支撐 | 角色 | 證據級 |
|---|---|---|---|---|
| **Somatic haplotagging** | LongPhase-S：germline HP1/HP2 + somatic 子單倍型 HP1-1/HP2-1（HP1-1 = germline-H1 + somatic-ALT，`HaplotagStrategy.cpp:505-516`）| 強 infra；V5/V6 self-phasing fork 已建 | **重建引擎（genetic）** | L1-L2 |
| **Subclonal reconstruction** | somatic SNV 在 LOH 區形成同-haplotype 子家族 = 子克隆譜系 | NG=2 Inner same-HP1 > Outer **7/7**（W=28, p=0.0078, gap 0.345）；Inner×NG=2 ≥93% same-hap 6/6；AF→NGroups 7/7（HD-4 確認 = phasing/occupancy 非甲基）| **主產品（Grade B+ ⭐3）** | L2（circularity caveat）|
| **Methylation profiles** | read-level 5mCG/5hmCG | ASM 真但 ~4% sens、非方向、copy-dosage-independent、**不判別 TP/FP**、subclone-甲基 = germline-allelic 非 somatic-specific | **characterization + 誠實負結果（非重建驅動）** | L2 |
| **Nanopore sequencing** | ONT 長讀（跨距離 phasing + 直接甲基）| 全程 substrate | 平台 | L1 |

> 🔴 **最高優先框架紀律**：subclonal reconstruction 由 **haplotagging 驅動**；**甲基不獨立重建子克隆**（甲基側 subclone 訊號是 germline-allelic、不 somatic-specific、不能 filter）。把甲基寫成「重建子克隆的判別力」= 與三~四道 NEGATIVE 直接矛盾，必撤。

---

## 2. 既有成果 → 論文章節映射（7 線收斂）

| 既有研究線 | → 論文章節 | 可寫狀態 |
|---|---|---|
| **T1 LOH-constrained phasing** | **§Result 1 — 脊柱：somatic haplotagging 重建 LOH-constrained 子克隆結構** | POSITIVE（Grade B+，需 HD-1 caveat）|
| **T6 AF→NGroups（HD-4 = phasing）** | §Result 1 支撐 — haplotype-occupancy signature（隨 AF 變）| POSITIVE（改述為 phasing 非甲基）|
| **T2 ASM 存在 + T5 跨樣本** | §Result 2 — 子克隆內 ASM characterization（存在、小、非方向、跨 3 癌種現象復現）| READY-AS-CHARACTERIZATION |
| **T3 cis-test / copy-partition** | §Result 2 支撐 — chr17/TBC1D16 唯一乾淨 cis exemplar；BRCA2 = 小 copy/subclone-confounded 例 | READY-AS-METHODS |
| **T4 methyl→TP/FP filter** | **§Methods-Neg — 甲基不能 filter/判別子克隆（死四道）** | READY-AS-NEGATIVE（最強）|
| **本 session: ASM 工作站/驗證/SEQC2-CN** | §Methods-Neg 證據基座 + §Confound-check | ✅ 已落地（見 §4）|
| **T7 methyl-assist phasing (umtag)** | §Discussion / future-work + 1 張 V10 copy-control 圖 | NOT-READY（roadmap 非 result）|

---

## 3. 現在就能寫的（可寫清單，caveat 必隨行）

### POSITIVE（脊柱）
- **§R1 主**：NG=2 Inner same-HP1 TP-rate > Outer cross-het，7/7，W=28 p=0.0078125, median gap 0.345 CI[0.104,0.459]；機制 `LabelTest.cpp:265-302`。**caveat**：① Grade **B+** 非 A；② p 量方向一致性非 effect（magnitude 跨樣本 11×）；③ n=7 生物上 n=6（HCC1395 算兩次）；④ single-pipeline ClairS-TO → ⭐3；⑤ **R-SELFREF circularity 首段明說**。
- **§R1 負控**：paired gap median≈−0.0002, p=0.578 NS → 訊號 TO-specific、乾淨。

### NEGATIVE（methods 主體 — 最 reviewer-proof）
- **甲基不能 filter 變異（死四道）**：LOSO 100% circularity（in-dist +0.0224→held-out −0.0001）｜\|Δβ\| **AUC=0.505**（null 內）｜powered 6-sample ~4% sens、COLO829 TP≈FP｜Zone-Aware QS NEG + CNV-zone filter（4 月關閉）。
- **copy-dosage 非 driver**：MW neutral-vs-gain p=0.6183、signed ρ(\|Δβ\|,CN)=−0.083 反向 → dosage REFUTED（全程式最乾淨單一結果）。
- **high-Δβ→FP 因 regression-to-extreme**：no-clustering OR=8.6 + LOH OR=4.1 + small-subhap OR=5.8（非低覆蓋 OR=0.87）；combo TP 0.86% vs FP 7.97%（~9×）。
- **subclone-甲基 = germline-allelic 非 somatic-specific**：blind-ARI（imprinting 正控 GNAS/RB1=1.000）somatic TP 0.135 < germline-het null 0.177, p=0.922。

### CHARACTERIZATION
- ASM 真、6/6 excess-over-null +0.101–0.241（3 癌種，現象層）｜chr17/TBC1D16 唯一乾淨 cis（perm p=0.001）｜BRCA2 = 小 copy/subclone-confounded 例（勿寫精確 %、勿當 cis 錨）。

---

## 4. 本 session 對此論文的具體貢獻（已落地，commit 8dbb931→5e69a99）

把 §Methods-Neg + §Confound-check 從「結論」升級為「有全位點證據 + 互動可檢視 + 獨立驗證 + 真值對照」：
- **全 30,350 位點（chr1-22，TP 29,723/FP 627）逐位點雙圖**（甲基 read×CpG + ISM 距離矩陣）+ 判斷/匯出工作站 → 子克隆內 ASM 的 read-level 證據基座。
- **TP/FP-free 5 層獨立驗證**（HP-permutation + split-half + 雙股 + bootstrap，非循環）：Tier A/B 85-86% 通過 → ASM 結構真實可複製（characterization 站得住）。
- **CramersV 可靠性閘控機制解明**（Cochran 期望<5 歸零；2,788 over-strict 91% 為稀疏歸零）→ 「過嚴格」是統計安全閥非 bug。
- **SEQC2 CN 真值 × ISM 分類關連**：ASM 在 SEQC2-LOH 最低（富集 0.46×）= **需兩條單倍型才可能有 ASM**（生物學方向，非 CN 偽造）→ 直接 back-stop copy-confound 警告，強化 §Methods-Neg 2。
- 教授級說明書（文字版 + 圖解版）+ provenance 表（樣本/突變/範圍）。

> 這些**全部支撐 ASM-characterization 軸 + 負結果軸**，**不**新增任何「甲基正向重建子克隆」claim（與框架紀律一致）。

---

## 5. 🔴 決定性 gate + 邊界（投稿前必守）

- **HD-1（GO/NO-GO，你決定，AI 不代決）**：subclonal reconstruction 的 headline（Inner same-HP1）**by construction 部分循環**（Inner 用 somatic-attributed HP1-1 定義）。只有兩條路：
  1. **跑 R-SELFREF 全 7 樣本 flag-on 對照（~25-50hr C++）** → 讓重建脊柱防彈、可寫成 positive 主結果；或
  2. **把 phasing/重建降為 characterization observation**，讓三~四道 NEGATIVE 扛論文（接受 methods/negative 論文定位）。
- **邊界（每節隨行）**：① 單樣本 HCC1395（+ 配對正常）→ ⭐3；② single-pipeline ClairS-TO；③ chr1-22（無 X/Y）、SNV only；④ 甲基是 characterization 非重建驅動；⑤ 「重建」= LOH-constrained 同-haplotype 子家族分割，**非**完整 clonal phylogeny/CCF tree（勿 overclaim 成 full subclonal tree tool）。
- **5 大 framing trap**：characterization-set 當 filter｜cis-anchor 當 subclone/copy-confound（BRCA2）｜private 當 underpowered（0/38）｜phenomenon-recurrence 當 cis-recurrence（6/6）｜held-out 當 true-rescue（umtag）。

---

## 6. 整理放緩的任務（consolidate 後狀態）

| 任務線 | 狀態 | 處置 |
|---|---|---|
| ASM 顯示/工作站/驗證/SEQC2-CN（本 session）| ✅ 完成、已 commit | **凍結為 §Methods-Neg 證據基座**，不再加功能（除非寫論文需要特定圖）|
| 甲基→TP/FP filter | ❌ DEAD（死四道）| **勿重開**；寫成最強 NEGATIVE |
| G6 三軸 framing | → reframe 為本聚焦論文 | 標籤改：Grade B+、ASM=characterization、BRCA2 退役、umtag=future-work |
| G1 ASM survey（chr5-22 dup-bug）| 暫緩 | 非關鍵路徑；characterization 已足夠寫 |
| H019-H023（queue）| ASM funnel/CN/LOH 探索 | 已併入 characterization；非新方向 |
| umtag (T7) | future-work | Discussion + V10 圖；ledger 註冊待 HD-5 |

---

## 7. 就緒度誠實評估（可以開啟論文目標了嗎？）

**可以開寫 = YES，但有一個你必須先做的決定（HD-1）。**

| 問題 | 答案 |
|---|---|
| 有沒有可寫的 positive 主結果？ | ✅ 有（LOH-constrained 子克隆重建，Grade B+）— 但需 HD-1 決定寫成 positive-spine 還是 characterization |
| 有沒有可寫的 methods 主體？ | ✅ 有（死四道 NEGATIVE，最 reviewer-proof，今天就防彈）|
| 甲基的角色清楚了嗎？ | ✅ characterization + 誠實負結果（非重建驅動）|
| 證據基座扎實嗎？ | ✅ 全位點 + 獨立驗證 + 真值 confound-check（本 session 落地）|
| 還缺什麼？ | 🔴 **HD-1 決定**（跑 R-SELFREF or 降 characterization）；⭐4 需 COLO829（非必要）；clonal-tree 若要當 tool 需新建（建議不做、維持 characterization）|

**建議撰寫關鍵路徑**：① 先定 HD-1 → ② 先寫死四道 NEGATIVE methods 段（今天防彈）→ ③ 依 HD-1 grade 寫重建脊柱 → ④ 寫 ASM characterization + chr17 exemplar + 跨樣本現象 → ⑤ 改 stale 文件（HD-3 BRCA2）。

---

## See Also
- 證據 SoT：`InterSubMod/knowledge/11_external_literature/10_paper_readiness_convergence.md`
- 現況地圖：`InterSubMod/docs/concepts/2026/06/20260608_研究現況地圖_整體目標與流程_給其他AI_01.md`
- 本 session 工作站：`InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/display_v2/20260610_judgment_workstation_01.html`
- ASM-CN 機制：memory `project_asm_locus_display_and_cramersv_reliability_gate`
