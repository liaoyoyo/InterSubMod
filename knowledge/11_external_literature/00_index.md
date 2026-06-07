---
id: ism-kb-11-external-literature-index
name: "External Literature Landscape — 7 角度交叉驗證索引"
description: "ISM 癌症長讀甲基化研究 7 大問題的外部文獻地景：每問題的內部發現(L1-L2) × 外部佐證(L3) × 衝突度，加跨問題依賴鏈(dead/live/open dreams)。經 workflow wf_37b2cc97-663 實際 WebFetch 查證。"
status: active
last_verified: 2026-06-05
content_nature: paper-derived
doc_type: explanation
verified_scope: "external papers verified via WebFetch in workflow wf_37b2cc97-663 (2026-06-05); internal claims cross-referenced against research story §7 + master_draft (validated SoT)"
decision_refs: []
related_ids: []
tags: [literature, external-corroboration, asm, methylation, long-read, subclone, phasing, tp-fp, cancer, landscape, index]
canonical_paths: [11_external_literature/00_index.md]
alias_paths: []
---

# External Literature Landscape — 7 角度交叉驗證索引

- **一句結論**：ISM 的整體方向（read-level haplotype-linked 甲基 characterization + tooling、誠實負結果）在文獻中**大方向被證實、無真正反駁**；4 個「夢」已死、3 個 LIVE（最強 = LOH 下甲基輔助 phasing 的白地）、3 個 OPEN/未證；所有 ASM 結論受**單樣本 HCC1395 ⭐3 / 單一 pipeline** 天花板限制。
- **適用對象**：論文 Intro/Discussion 骨架、PI/協作者全景說明、reviewer 防守。
- **這份怎麼來的**：7 個研究問題各派 1 個 scout（5-7 次 WebSearch）+ 1 個 verify agent（逐篇 WebFetch 實證 PMID/DOI + 分類 CONFIRM/REFUTE/DIFFERENT-VIEW vs 我們的內部 claim），再跨綜合 + 完整性 critic。**證據分級：內部 L1-L2（本專案實測）/ 外部 L3（文獻，附真實識別碼）。內外不對齊一律誠實標，不假裝外部證實內部。**

> ⚠ **本文件性質**：外部文獻 paper-derived 索引；數字為「文獻報告值」與「本專案已 validated 內部值」的對照，非本輪新跑分析。引用前查每篇 credibility 欄與 `09_references_and_provenance.md` 的識別碼核對註記。

---

## 7 問題 → 子文件

| # | 問題 | 子文件 | 一句裁決 |
|---|------|--------|----------|
| Q1 | 目標對齊：長/短讀方法、時間軸、paired vs TO | [01_goal_alignment_methods_timeline.md](01_goal_alignment_methods_timeline.md) | **CONFIRM** — 長讀 read-level 前提 + array→WGBS→ONT 時間軸 + paired 為標準，全被背書 |
| Q2 | 癌症甲基差異位點的做法與困難 | [02_finding_methyl_diff_loci.md](02_finding_methyl_diff_loci.md) | **CONFIRM 困難 + 部分 REFUTE「找不到很多」** — 失序公認，但 Lin/POG 在 cohort 找到 ~23.6K aDMR/腫瘤 |
| Q3 | ISM 方法是否合理有效×相似/差異工具 | [03_ism_method_comparison.md](03_ism_method_comparison.md) | **CONFIRM framing，niche 收窄** — structure≠disorder 成立；cis-test + LOH 耦合零競爭，read-read 距離單獨非新 |
| Q4 | subclone 研究對印我們觀察 | [04_subclone_landscape.md](04_subclone_landscape.md) | **CONFIRM 概念，「reconstruction」過度** — 缺 phylogeny/CN-deconv/單細胞三大金標準 |
| Q5 | 甲基區分 TP/FP 或 germline/somatic | [05_methyl_tp_fp_germline_somatic.md](05_methyl_tp_fp_germline_somatic.md) | **DEAD 但被背書** — 循環性是公認 ML 陷阱(Kapoor/Soneson)；無 caller 用甲基當 standalone FP 訊號 |
| Q6 | 癌症甲基輔助 phase/tag | [06_methyl_assisted_phasing.md](06_methyl_assisted_phasing.md) | **LIVE 白地** — germline 機制已證(MethPhaser/HapBridge)，**tumor/LOH 沒人做過(含我們)** |
| Q7 | ASM/cis 與癌症影響 | [07_asm_cis_cancer_impact.md](07_asm_cis_cancer_impact.md) | **方法 CONFIRM、方向部分** — Martin-Trujillo CN 背書最強；Do2020 是口徑差非矛盾；BRCA2 hypo≠canonical hyper |
| ★ | 跨問題依賴鏈 + 可信度衝突矩陣 | [08_cross_question_correlation.md](08_cross_question_correlation.md) | 理想鏈逐環 REAL/BLOCKED/POSSIBLE |
| ★ | 完整參考文獻 + provenance + 未覆蓋工具 | [09_references_and_provenance.md](09_references_and_provenance.md) | 識別碼核對 + critic gap |
| 🔴★ | **論文就緒收斂**（usable/failed/needs-thought + GO-NO-GO）| [10_paper_readiness_convergence.md](10_paper_readiness_convergence.md) | **2026-06-08 跨 7 線獨立稽核**：可寫清單 + 2 真矛盾 + R-SELFREF 決定性 GO/NO-GO + BRCA2→chr17 demotion |

---

## 跨問題依賴鏈 — 哪些夢死了、哪些活著、哪些開放

> 你的理想情境鏈：**找甲基差異位點(Q2) → 助 phase/tag(Q6) → ISM 改造成可做這些(Q3) → 解 subclone(Q4) → 分 TP/FP(Q5)；底層 ASM/cis 生物(Q7)、長讀使能(Q1)。** 逐環裁決（完整見 08）：

| 環 | 內容 | 裁決 |
|----|------|------|
| L0 | 長讀使能 read-level haplotype-linked 甲基（Q1）| 🟢 **REAL**（field-validated，無爭議）|
| L1 | 找「大量、乾淨、可統計」的甲基差異位點（Q2）| 🔴 **BLOCKED**（單樣本群體統計）／🟢 REAL（稀有但真實、IGV 可見、cohort 可定位）|
| L2 | 甲基差異 → 助 phase/tag（Q6）| 🟡 **POSSIBLE（最強 LIVE 白地）** — 機制 germline 已證、tumor/LOH 無人做 |
| L3 | ISM 改造成做這些（Q3）| 🟢 REAL（distinct 方法）／🟡 POSSIBLE（是否「更好」未證，無 head-to-head）|
| L4 | 解 subclone（Q4）| 🔴 BLOCKED（「reconstruction」）／🟡 POSSIBLE（「characterization」）|
| L5 | 分 TP/FP（Q5）| 🔴 **BLOCKED — 已 concluded NEGATIVE，鏈在此終止** |
| U1 | 底層 ASM/cis 生物（Q7）| 🟢 REAL（存在+cis+CN-controlled）／🟡 方向部分／🔴 recurrence BLOCKED |

**夢的狀態**：
- 🔴 **DEAD**：甲基當 TP/FP filter（L5）；跨 region 二次打擊**順序**（L4）；「單樣本群體統計找出大量乾淨位點」（L1）；BRCA2 當 TSG-silencing 證據（Q7 方向）。
- 🟢 **LIVE（最強）**：LOH 下甲基輔助 phasing/tagging（L2，白地）；LOH-constrained phasing 主軸（Grade A）；ASM 存在+cis+CN-controlled（U1）。
- 🟡 **OPEN/未證**：ISM 方法**優於** pycoMeth/cvlr（L3，無 head-to-head）；subclone **characterization**（非 reconstruction）（L4）；5mC/5hmC 分軌新穎性（文獻 silent）。

---

## 最可防守的論文主張（給 reviewer 看）

> *「在 ONT 長讀癌症基因組中，allele/haplotype-specific 甲基化（ASM）真實、cis-driven、且設計上獨立於 copy number——但它是**結構/連結的 characterization 訊號，不是 somatic variant filter**。我們貢獻：(a) normal-anchored somatic cis-test，用 HP1-vs-HP1-1 軸 held-constant CN/ploidy/alignment 來分辨真 cis-driven ASM 與 drift；(b) 把 haplotype-linked 甲基用於 **LOH 區（germline-het anchor 消失處）的 phasing/tagging**——這是現有 germline methyl-phaser（MethPhaser/NanoMethPhase/HapBridge）明確留為 future work 的白地。」*

此主張能撐過每個內部負結果：定位為 characterization+tooling（對齊 filter 的 NEGATIVE）、主打 CN-confound 擊破（唯一鐵證外部背書 = Martin-Trujillo）、押 LOH-phasing 白地（真正未填、可行性已證），不過度宣稱 subclone reconstruction 或 recurrence。

**3 大 reviewer 脆弱點 + 文獻先發防守**（完整見 08 §4）：
1. **單樣本/單 pipeline** → 用 Sakamoto22（private hap-DMR 證實我們 private 0/38 = cohort-一致而非失敗）+ POG（recurrence 是 cohort-scale future work）防守。
2. **read-read 距離不新** → 自己先窄化（cvlr/pycoMeth/Metheor 已佔 clustering），只押**未佔的組合**（somatic cis-test + LOH/CN 耦合 + 5mC/5hmC 分軌）。
3. **BRCA2 hypo 與 canonical hyper 矛盾** → 自己先標方向，澄清 HP-axis Δβ 是**等位結構差（含 allele+copy），非 promoter-mean silencing**，明說 BRCA2 不當 TSG-silencing 證據。

---

## 與既有專案文件的關係

- 研究故事母稿（內外證據 §7 表 C1-C8）：[../../docs/concepts/2026/06/20260603_研究故事與困難敘述_甲基haplotype論文_內外證據_01.md](../../docs/concepts/2026/06/20260603_研究故事與困難敘述_甲基haplotype論文_內外證據_01.md) — 本地景是其**大規模外部擴展**（Q1/Q4/Q6 為新增深度）。
- 內部數字 SoT：[../../docs/reports/validated/2026/06/20260602_研究週報_20260526_20260602_phasing_GradeA與ASM收斂/master_draft.md](../../docs/reports/validated/2026/06/20260602_研究週報_20260526_20260602_phasing_GradeA與ASM收斂/master_draft.md)
- 五大目標：[../01_project_overview/02_five_research_goals.md](../01_project_overview/02_five_research_goals.md)
- 既有結論索引：[../09_conclusions/00_index.md](../09_conclusions/00_index.md)
- HTML 視覺版：`00_index.standalone.html`（同目錄 companion）

---

## Provenance

- 外部文獻：workflow **wf_37b2cc97-663**（16 agents，7 scout + 7 verify + cross + critic，210 tool uses，~8 min，2026-06-05）實際 WebSearch/WebFetch 檢索；每篇附 relation + 真實 PMID/DOI；無法 fetch 確認者標 UNVERIFIED。
- 內部數字：沿用 `master_draft.md` / 研究故事 §7（validated SoT），本輪**未跑新分析**。
- 撰寫此 KB 文件的 Write 與外部查證的 WebFetch **在不同 batch**（§13.0 物理隔離）。
