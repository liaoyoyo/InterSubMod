<!--
建立時間: 2026-06-10
狀態: PLAN (method-comparison study, part 6/6 — Phase B 實機 benchmark 方案；尚未執行)
報告類型: benchmark_plan
受眾: 執行者(自己/agent) · PI 核准
framework: 分階段執行計劃(B0-B4) + 工具安裝矩陣 + 對照設計
data_sources:
  - 03_method_comparison_matrix.md, 05_improve_learn_recommendations.md
  - _assets/survey_digest.md (GAPS: 無工具在 PATH; ASMS 最像; 無 runtime 比較)
provenance_note: 本檔為「規劃」,不含實測數字。實機執行屬 Phase B(下一步), 須背景 Bash 跑 + 落檔 + Read 驗 + 才寫結果報告(§13.0)。工具版本/可裝性以執行時實測為準。
-->

# 06 — Phase B 實機 benchmark 方案（下載 + 運行；⚠ 尚未執行）

> **這份是什麼**：第 6 部分 —— 把「先整理分析」之後的「**實際下載運行分析結果**」規劃成**可執行、分階段、資源安全**的方案。**本檔不含任何實測數字**；執行屬下一步（須背景 Bash 跑 → 落檔 → Read 驗 → 才寫結果報告，§13.0）。
> **核心問題（synth gap 指出）**：紙上比較已知 ISM 的組合無人佔，但**「組合無人佔 ≠ 更好/有用」**。Phase B 要回答：**當 modkit dmr 說「無率差」但 ISM PERMANOVA 說「有結構」時，ISM 的 call 是真訊號還是噪音？**

---

## L0 — 一句話方案

在主樣本 **HCC1395**（有 SEQC2 truth + 既有 haplotagged BAM）上，對一組定義好的位點（BRCA2、TBC1D16、imprinting ICR 正控 + panel），**並排跑 ISM vs 3 個最相關對手**（modkit dmr＝率差基準 / cvlr＝read-clustering 近鄰 / ASMS＝最像的 no-phasing ASM），產出 **concordance 表 + discordant-case 深析 + runtime**，判定 ISM 的 PERMANOVA/cis-test 是否帶來對手沒有的真訊號。

---

## L1 — 工具安裝矩陣（難易 + 來源 + 角色）

| 工具 | 角色（對哪個軸）| 安裝難度 | 來源 | 備註 |
|------|---------------|---------|------|------|
| **modkit** ⭐ | 軸A 率差基準（ONT 官方）| 🟢 易（單 Rust binary / conda）| github.com/nanoporetech/modkit releases | 先裝這個；`modkit pileup --phased` + `modkit dmr pair` |
| **ASMS** ⭐ | 軸C 最像 ISM（no-phasing read-clustering ASM）| 🟡 中（bioRxiv，須找 repo）| Raineri 2024 bioRxiv 2024.12.18.629129（repo 待查）| **最關鍵對照**；若 repo 不可得則降為紙上 |
| **cvlr** | 軸C read-clustering EM | 🔴 難（**GitHub 404**，作者 0 repo @2026-06-10）| 待查 CNAG/聯絡作者 | 可能無法裝 → 紙上 + 自實作 Bernoulli-EM 替代 |
| **pycoMeth** | 軸A/F haplotype-aware（read-level）| 🟢 易（pip / conda）| github.com/PMBio/pycoMeth | 需 Nanopolish/MetH5 前置（較重）|
| **NanoMethPhase + DSS** | 軸F HP1/HP2 + 軸A DSS 統計 | 🟡 中（pip + R/Bioconductor DSS）| github.com/vahidAK/NanoMethPhase; Bioc DSS | DSS 同時是 #1 改進的參考實作 |
| **DAMEfinder** | 軸B tuple（短讀，**僅概念對照**）| 🟡 中（R/Bioconductor）| Bioc DAMEfinder | 短讀工具，HCC1395 ONT 不直接可跑 → 概念對照 |
| **Metheor / WSHPackage(qFDRP)** | 軸E disorder 對照 | 🟢🟡 | github.com/dohlee/metheor; MPIIComputationalEpigenetics/WSHPackage | qFDRP 是 #7 改進參考 |

> **可行性現況（synth GAP 實測）**：執行時 `which modkit`、`~/.cargo`、R/Bioc DSS/cvlr/DAMEfinder/pycoMeth **皆不在 PATH/磁碟** → Phase B 必含**安裝階段**。**modkit 與 ASMS 是優先且最直接 ONT 可比**；短讀工具（DSS/DAMEfinder/Metheor）需 HCC1395 的 bisulfite 資料（無）或僅作方法藍圖，**不強跑**。

---

## L2 — 分階段執行計劃（B0–B4）

### 既有可用輸入（已確認路徑）
- **HCC1395 haplotagged BAM**（含 MM/ML + HP）：`/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam`
- **MSA Level-1 raw methylation**（已抽）：`research/tsg_promoter_asm_reviewer/genome_survey_v2/.../level1_raw_methylation_details.tsv.gz`
- **ISM binary**：`/big8_disk/liaoyoyo2001/InterSubMod/build/`（read-read 距離/clustering/PERMANOVA）+ MSA `build/bin/msa`
- **位點集**：BRCA2 chr13:32,315,128（~80% copy）、**TBC1D16 chr17:79,991,120（唯一乾淨 cis）**、+ imprinting ICR 正控（如 H19/IGF2、SNRPN）、+ 10-20 control loci（TP/FP/het-NULL，manifest 已有）。

### B0 — 安裝 + smoke test（先做，輕）
- 裝 modkit（binary）；試裝 ASMS / pycoMeth；確認 ISM/MSA 可跑。
- preflight：`mount | grep tmp`、`export TMPDIR=/big7_disk/...`（避免 /tmp 寫滿，memory 教訓）；磁碟 + CPU/mem 檢查（`/infra-ops`）。
- smoke：modkit pileup 單 region 跑通；ISM 單 region 跑通。

### B1 — modkit 率差基準（最易，先出結果）
```bash
# 1) phased pileup (HP1/HP2 bedMethyl)
modkit pileup --phased --cpg --ref <ref.fa> HCC1395_tagged.bam pileup_out/
# 2) HP1 vs HP2 DMR (Bayesian LR + MAP p + Cohen's h)
modkit dmr pair -a pileup_out/hp1.bedmethyl.gz -b pileup_out/hp2.bedmethyl.gz \
  --ref <ref.fa> --regions <loci.bed> -o modkit_dmr.tsv
# 3) entropy (disorder 對照)
modkit entropy --regions <loci.bed> HCC1395_tagged.bam > modkit_entropy.bedgraph
```
- 產：每位點 modkit Δβ(率差) + Cohen's h + LR score；對照 ISM Δβ（HP-axis）。
- ⚠ modkit dmr 是 HP1 vs HP2（germline 軸），ISM 主軸是 HP1 vs HP1-1（somatic）—— 對照時須明標軸差。

### B2 — read-clustering 近鄰（ISM 主場對照）
- **ASMS**（若 repo 可得）：對同位點集跑 no-phasing read-clustering ASM；比 ISM 的 cluster 指派 + PERMANOVA 是否一致；ASMS 無 normal-anchor/PERMANOVA → 看 ISM 額外能力是否在 discordant 位點顯現。
- **cvlr**（若可裝）：EM Bernoulli-mixture k=2；比 cluster 指派 ARI vs ISM UPGMA 切群。
- 若兩者皆不可裝 → 在 ISM 內以「無顯著性檢定的純 clustering」模擬對照（自實作 baseline）。

### B3 — 嚴謹統計對照（#1/#7 改進的實證）
- **NanoMethPhase + DSS**：HP1 vs HP2 用 DSS beta-binomial（dispersion shrinkage）算 region DMR；對照 ISM 的 Fisher/Wilcoxon Δβ —— 量化「Fisher 是否在低 HP 覆蓋膨脹顯著性」。
- **qFDRP**（WSHPackage / 自實作）：算各位點 qFDRP（coverage-robust disorder），與 ISM PERMANOVA 做 **2×2**（高 qFDRP×非顯著 PERMANOVA=stochastic erosion；高 qFDRP×顯著=結構化 ASM）—— **直接操作化「結構非失序」**。

### B4 — 綜合 + discordant-case 深析
- **concordance 表**：每位點 × {ISM PERMANOVA p / ISM Δβ / modkit Δβ+Cohen h / ASMS cluster / DSS DMR / qFDRP}。
- **discordant cases（核心交付）**：找「modkit 說無率差 but ISM PERMANOVA 說有結構」與反向案例，**人眼 + IGV + read×CpG 圖**判定 ISM call 是真訊號還是噪音 → 回答 synth gap。
- **runtime/scalability**：各工具單位 region 與全基因組外推 wall-clock/mem（ISM O(N³) UPGMA + 999perm vs modkit Rust）。
- **正控**：imprinting ICR 上 ISM 應偵得已知 ASM（#8）。

---

## L3 — 執行紀律（§8 + §13.0，務必遵守）
1. **長 compute 不放 workflow agent step** —— modkit pileup / ISM 全跑 / DSS 等用 **背景 `Bash(run_in_background)`** 落檔，**不**包進 Dynamic Workflow（workflow 只用於「讀已落檔結果」的匯總）。
2. **資源 preflight** —— N 個 BAM 並行恐撞 CPU/mem → 預設**循序背景跑**；`export TMPDIR=/big7_disk/...`（/tmp 災情 memory 教訓）。
3. **§13.0 序列** —— 跑分析 → 落 .json/.tsv → **Read 讀回確認非 error/未完成** → **才**寫 benchmark 結果報告（`07_benchmark_results.md`，Phase B 完成後新增）；產數字與寫報告**絕不同 batch**。
4. **NO COMPLETION CLAIMS WITHOUT FRESH VERIFICATION**（§13.7）—— 宣稱「modkit 比 ISM X」前必有讀回的真值支持。
5. **資料夾**：所有 Phase B 產物落 `benchmark/`（見下），與 Phase A 文件分離。

### 預計資料夾佈局（Phase B）
```
20260609_ism_vs_external_methylation_tools/
├── 00_INDEX.md ... 06_benchmark_plan_phaseB.md   ← Phase A（本批，已完成）
├── _assets/                                       ← workflow 原始結果 + digest
└── benchmark/                                     ← Phase B（執行時建）
    ├── tools/        (安裝的 modkit/ASMS/...)
    ├── loci.bed      (位點集)
    ├── runs/         (各工具 raw 輸出, 落檔)
    ├── compare/      (concordance 表 + discordant cases)
    └── 07_benchmark_results.md  (讀回真值後才寫)
```

---

## 決策點（需用戶確認才進 Phase B）
- **B 範圍**：先只跑 **modkit + ISM**（最易、最直接 ONT 率差 vs 結構對照），還是含 ASMS/cvlr/DSS 全套？
- **位點集大小**：先 ~5 位點（BRCA2/TBC1D16/2 ICR/1 control）pilot，還是直接 panel（20+）？
- **是否現在啟動 B0 安裝**（modkit 單 binary，~10 min，可背景）？

> 我**未自動啟動 Phase B**（涉及下載外部工具 + 長 compute；§1 需用戶當輪明示或核准）。Phase A（文獻+方法+交叉+建議）已可獨立交付。
