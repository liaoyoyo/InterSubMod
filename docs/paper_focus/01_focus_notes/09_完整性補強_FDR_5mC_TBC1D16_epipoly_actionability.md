<!--
建立時間: 2026-06-10
狀態: focus (完整性補強 — 5 個遺漏的詳細內容 + 佐證方法 + 圖片 spec)
報告類型: paper_focus_completeness_supplement
受眾: 廖子游 · PI · 執行 FDR/5mC/figures 的其他 session
provenance_note: FDR/5mC/actionability 為方法設計(🔵)+待跑任務；TBC1D16 為外部文獻(L3，web verified)；epipoly 區隔用內部 O11(🟢 NEGATIVE)。
-->
<!-- provenance-verified: TBC1D16 = WebSearch 2026-06-10 verified(Nat Med 2015/Cancer Discovery 2015/Clin Epigenetics 2019)；O11 epipoly NEGATIVE = memory project_O11_heterogeneity_negative(🟢)；FDR/5mC = 待跑任務 spec。 -->

# 完整性補強 — FDR / 5mC / TBC1D16 / epipolymorphism / actionability

> **L0 一眼結論**：5 個 reviewer 級遺漏的補強。**最該先做 = FDR/null 校準（#24）**——它若失敗會同時砍掉「規則 12,868」和「例外 chr17」；其餘 4 個是**誠實口徑 + 內容補充**（5mC 範圍縮、TBC1D16 生物動機已查到、epipoly 區隔用自家 O11、actionability + 圖 triage）。
>
> **L1 重點邏輯**：① FDR 是唯一可能改結論的（必先做）；② 5mC 縮範圍即解（5mC 主結論 + 5hmC 觀察/future）；③ TBC1D16 是真癌症甲基基因（web verified）→ 例外有生物動機，但跨癌+bulk-vs-ASM caveat，COLO829 是理想驗證；④ epipoly 區隔我們有自家 O11 NEGATIVE 當證據；⑤ actionability + 圖 triage 是投稿包裝。

---

## §1 FDR / null 校準（#24，必補，最該先做）

**問題**：(a) 規則的 12,868 reliable = 3.87% 通過 gate，但**低於 p≤0.05 的 5% 噪音地板** → 沒 null 校準無法宣稱是真訊號；(b) 例外 chr17 perm p=0.001 跨 816 cis-test → **Bonferroni 0.001×816≈0.82 (NS)**，可能是多重檢定 artifact。

**要做的方法（給 #24）**：
| 對象 | 方法 |
|------|------|
| 12,868 reliable | **label-shuffle null**：permute HP/somatic 標籤重算 reliable count → 經驗 FDR = null_reliable/observed；**或** per-locus PERMANOVA p → BH-FDR across 332,705 → q<0.05 計數 |
| chr17 cis | **Bonferroni + BH-FDR** across 816；**genome-wide permutation**：shuffle 後整條 pipeline 產幾個「乾淨 cis」→ null 期望 vs observed=1 |
| 共同 | 🔴 **必含 n_reads/coverage 校正**（O11 教訓：epipolymorphism AUC 0.845→0.530 after n_reads correction；不校正會假陽性）|

**🟢 可以佐證的其他方法（用戶問的）**：
1. **統計**：Storey q-value / local-FDR；empirical-Bayes（limma-style）on methylation difference。
2. **外部佐證**：chr17/TBC1D16 是已知癌症甲基基因（§3）→ 文獻 orthogonal 支持（但 bulk 非 ASM）。
3. **生物複製 > 統計**：COLO829（黑色素瘤，TBC1D16 主場）重現 chr17 = 最強佐證（取代單純 p 校正）。
4. **bulk WGBS/array ASM** 對這些細胞株若有 → cross-platform 一致性。

> **誠實後果**：若 chr17 不過 FDR → 例外降為「suggestive，靠生物動機(TBC1D16)+COLO829 複製撐」非統計唯一；若 12,868 不過 null → 規則改述「reliable 分群率不顯著高於噪音，更強化甲基非 somatic 判別器」（**反而強化 de-confound**）。**任一結果 de-confound floor 都不垮。**

---

## §2 5mC vs 5hmC 範圍（#25，用戶指示：5mC 主結論 + 5hmC 觀察）

**決策（用戶 2026-06-10）**：論文**主結論以 5mC 為準**（well-established、dup-bug 不影響純 5mC）；**5hmC 只說「亦可輸入觀察、但未詳細驗證」**（不納主結論）。
- **要確認**：catalog/de-confound 結論在**純 5mC** 上重跑確認（剝離 5hmC）→ 確保「甲基跟 germline」是 5mC 結論。
- **5hmC**：Methods 一句「ISM 亦接受 5hmC 輸入；5hmC 生物意義不同（常 active gene-body），本研究未詳細驗證，列為 future observation」。
- **依賴**：MSA Level1 5mC/5hmC dup-bug 修（否則 5mC 純度受污染，重演 BRCA2 −0.054 artifact）。

---

## §3 chr17/TBC1D16 生物機制（#3，web verified 2026-06-10）

**外部文獻（L3，已查證）**：
- TBC1D16 = **Rab GTPase-activating protein**；**hypomethylation 重活化 cryptic 轉錄本 TBC1D16-47KD** → 經 **RAB5C → EGFR** 驅動**黑色素瘤生長/轉移**（Vizoso/Esteller et al., *Nature Medicine* 2015）。
- 「**Hypomethylation of TBC1D16 Drives Melanoma Metastasis**」（*Cancer Discovery* 2015）。
- hypomethylation 關聯**較短 PFS/OS**，但賦予 **BRAF/MEK 抑制劑敏感性**。
- 多癌種：EBF3 + TBC1D16 甲基變化關聯 tumour progression/metastasis（*Clinical Epigenetics* 2019, PMID 31383000）。

**對論文的意義**：
- ✅ **例外有生物動機**：chr17/TBC1D16 不是隨機位點，是**有名的癌症甲基化驅動基因**；方向 **hypo**，與我們 BRCA2「hypo≠canonical hyper」口徑一致（cancer methylation 不全是 hyper-silencing）。
- 🔴 **兩個 caveat（必寫）**：(a) 發表是**黑色素瘤 + bulk hypomethylation**，我們是**乳腺 HCC1395 + allele-specific cis-ASM** → 跨癌種 + bulk-vs-ASM 不可直接等同；(b) 須說「我們的 read-level allele-specific cis 觀察**與已知 TBC1D16 甲基-驅動生物學一致**」而非「證實」。
- ⭐ **COLO829 = 理想驗證樣本**：COLO829 是**黑色素瘤**（TBC1D16 主場）→ deferred 的 COLO829 工作現在有**強生物動機**：在 TBC1D16 的原生癌種驗 chr17 cis-ASM = §1 的最強生物佐證。

**Sources**（投稿前 /citation-verification 核 PMID/DOI）：
- Vizoso M, et al. *Nat Med* 2015, nm.3863 / PMC4968631.
- *Cancer Discovery* 2015 「Hypomethylation of TBC1D16 Drives Melanoma Metastasis」.
- *Clinical Epigenetics* 2019, PMID 31383000 / PMC6683458.

---

## §4 catalog 當資源的 actionability（別人拿來做什麼）

論文 catalog 不能只是「我們看了一眼」，要是**可用資源**：
1. **De-confounding 參考表**：別人做 read-level ASM 時，可查某位點是 reliable-but-germline (TAG-C) 還 untestable (TAG-F)，避免把 germline/copy 當 cis（直接服務文獻 Martin-Trujillo confound 警告）。
2. **負面 hotspot 清單**：TAG-D（high-Δβ FP-prone）+ TAG-F（untestable）= 「別在這些位點宣稱 cis-ASM」的 caution list。
3. **方法 sanity-check**：別的工具可拿 TAG-A（chr17）+ TAG-C 當已知 positive/germline 校準點。
4. **釋出形式**：catalog_skeleton.tsv（332,705×32 欄）當 **supplementary table + 可查詢**（建議附 schema 說明）。

---

## §5 Figure 主/補 triage（GB ~4-6 主圖）

| Fig | 內容 | 主/補 |
|-----|------|-------|
| **F1** rule+exception 整合圖（LOH→germline 主背景 + 罕見 somatic co-opt chr17；證據強度標）| **主** ⭐ 待製 |
| **F2** catalog P5 per-TAG 計數（de-confound 規則 主結果）| **主** ✅ 有 |
| **F3** filter 死四道對照 | **主** 待製 |
| **F4** copy-dosage falsification + within_d vs HP_d（P4）| **主** ✅ 有 |
| **F5** chr17/TBC1D16 cis-test 三路 + 生物機制 inset | **主** ⭐ 待製 |
| F6 NG=2 phasing forest | 補（contingent #21）|
| F7 ISM vs 工具 能力對照（#23）| 補 |
| P1/P2/P3/P6 分佈圖 | 補（supplementary）|

> GB 主圖鎖 **5 張**（F1-F5）：rule+exception / catalog / filter-dead / dosage-falsify / chr17。其餘 supplementary。

---

## §6 與 epipolymorphism 方法的概念區隔（用自家 O11）

**🟢 我們有自家證據**：O11 epipolymorphism heterogeneity 假說**已 NEGATIVE**（AUC **0.845→0.530** after n_reads correction；`memory/project_O11_heterogeneity_negative`）。

**區隔論述（Methods/Discussion）**：
1. epipolymorphism/methylation-entropy（Landan2012 等）測 read-level 甲基**異質性/熵** → 我們證實其判別力是 **n_reads-confounded artifact**（O11）。
2. **ISM 不同**：不靠熵/異質性判別，而是 **reliability-gated 分群（CramersV Cochran + PERMANOVA，控稀疏表）+ normal-anchored cis-test（控 germline/copy）**——正是 epipolymorphism 漏掉的 confound 控制。
3. → ISM 的貢獻是「**把 read-level 甲基異質性從 confounded 描述，升級為 de-confounded characterization**」。這條把 O11 的 NEGATIVE 變成方法論賣點。
> （與 #23 工具 benchmark 互補：#23 是 modkit/NanoMethPhase 功能對標；本節是 epipolymorphism 概念區隔。）

---

## §7 對映 + 任務
- **#24 FDR/null 校準**（必先做）→ §1；威脅 12,868 + chr17；含 n_reads 校正（O11）。
- **#25 5mC 純度驗證 + 5hmC observational**（依 dup-bug 修）→ §2。
- **#26 解釋圖**（F1 rule+exception / F5 chr17 機制 / epipoly 區隔示意）→ §5/§6，用 `/methods-example`。
- TBC1D16（§3）→ 併入論文架構 Discussion §3.2 + Results §2.5；COLO829 驗證動機 → T-VAL-3。
- actionability（§4）→ 論文「Data availability / Resource」+ Discussion。
