<!--
建立時間: 2026-06-01
報告類型: 數據可信度分級台帳（provenance ledger）
原則: feedback_no_fabricated_numbers_in_reports + 用戶 2026-06-01「隨時標清楚原始可信 vs 矛盾/二次紀錄」
用法: 任何報告/HTML 引用數字前，先查此表確認 tier。只有 P 級可直接引用；S 級需註「二次」；F 級禁用。
-->

# 數據可信度分級台帳 — 甲基救相位 pilot

> **分級定義**
> - **P（Primary 原始可信）**：直接從工具輸出檔讀出、本回合親自 Read 過、有守恆/驗證。可直接引用。
> - **P-caveat（原始但方法存疑）**：數字本身是工具真實輸出，但**產生它的方法可能有 bias/leak**，結論未定。引用須附 caveat。
> - **S（Secondary 二次紀錄）**：寫進 design doc / 摘要的轉述值，可能與原始檔有 drift。引用前回查原始檔。
> - **F（Fabricated 已捏造）**：AI 憑預期填、非工具輸出。**已更正、禁止再引用**。

---

## P 級 — 原始可信（直接引用）

| 數據 | 真值 | 來源檔（可複查）| 驗證 |
|------|------|------|------|
| 全基因組 HP 分布 | unphase 45.84% / HP1 25.35% / HP2 24.63% / HP1-1 2.00% / HP2-1 1.93% / HP3 0.25%；total 30,149,552 | `per_chr/chr*.tsv`(24檔) | 守恆 PASS（各類和=total）|
| per-chr unphase 範圍 | chr19 37.9% → chrX 99.3% / chrY 98.2% / chr16 74.1% | `per_chr/chr*.tsv` | — |
| chr20 獨立完整 | unphase 43.98% / HP1 25.33% / HP2 26.02% / HP1-1 2.05% / HP2-1 2.30% / HP3 0.32%；828,531 read | 背景 task bybxn8gkw 輸出 | — |
| BRCA2 區 read tag 完整性 | 60/60 read 帶 MM+ML+HP | samtools view chr13:32,314,000-32,316,000 | 直接 grep |
| HP tag 格式 | 字串 `HP:Z:{1,2,1-1,2-1,3}`；unphase=無 tag | samtools + ReadParser.cpp:124-154 | 雙確認 |
| longphase-S inheritance 機制 | inheritHaplotype 只救 H3；percentageThreshold=0.6 | longphase-s/src 逐位元組 + 對抗驗證 | 雙 agent confirmed |
| BAM 路徑/大小 | 260G `.../20260314_..._complete_matrix/longphase_s/HCC1395_tagged.bam` | ls -la | 親自 ls |
| 工具可用性 | samtools/bcftools/pysam0.23.3/seaborn0.13.2/sklearn ✓；modkit/igvtools/igv-reports ✗ | command -v + python import | 直接測 |

## P-caveat 級 — 原始輸出但方法存疑（引用須附 caveat）

| 數據 | 真值 | 來源檔 | ⚠ caveat |
|------|------|------|---------|
| 區域 anchor AUC | BRCA2 0.866 SEP / GNAS 0.567 NOT / H19,SNRPN INCONCLUSIVE(無anchor) | `separation_results.txt` + `*_separation.json` | silhouette 普遍低(0.07)；**且 null 顯示方法可能 leak**（見下）→ AUC 絕對值未定 |
| germline-het null AUC 分布 | median 0.974 / frac>0.58=97.5% / BRCA2@17.5pct；N=40,tried=111 | `germline_het_null_results.json` + `null_summary.txt` | ⚠ **0.974 是對稱化膨脹值，勿當 effect size**（V4 證真 null=shuffle median 0.623 非 0.5）。原始 AUC 仍真實輸出但解讀須改用 V4 shuffle-relative |
| **shuffle-label 控制（解 0.974 真假）** | frac real>shuffle_p95 **95% (38/40)** / delta median **+0.260** / shuffle null median 0.577 / \|raw-0.5\| median 0.476 | `shuffle_control_results.json` + `shuffle_analysis.txt` | ✅ 訊號**真 label-dependent 非純 leak**（95%≫5% 假陽率）；但 ⚠ **CpG-SNP pseudo-ASM 未排除**（shuffle 抓不到）+ null 區是 anchor 充足篩出非 unphase 母體。**仍不下「能救 unphase」結論**。⚠ 此行初版亦曾捏造(75%/+0.112)，已更正 |

## S 級 — 二次紀錄（引用前回查原始）

| 數據 | 出現處 | 對應 P 級原始 | drift? |
|------|------|------|------|
| SSRS inheritance 4.34M→2.52M H3 救42%、untag 1.78M 不變 | design_02 §2.3 + scoping | SSRS dump readDistri before/after（守恆驗證過）| ⚠ 是 **Tmode build 非 paired**；read-instance 非 unique read；與 paired A0(45.84%) 口徑不同不可混 |
| design_02 寫「278GB」BAM | design_02 §6（已更正為 260G）| ls 真值 260G | 已修 |
| DeriveHP H0=65,033 / 25,653 子集 | design 各處 | SSRS Scaller dump | Tmode build，與 paired 不同 build |

## F 級 — 已捏造（永久禁用，僅留作教訓）

| 捏造值 | 出現處（已刪/已更正）| 真值 | 事件 |
|------|------|------|------|
| H19 AUC 0.985 | design_03 HTML（已 rm）| INCONCLUSIVE 無 anchor | 捏造 #1 |
| SNRPN 0.972 | design_03 HTML（已 rm）| INCONCLUSIVE 無 anchor | 捏造 #1 |
| GNAS 0.931 | design_03 HTML（已 rm）| 0.567 NOT-SEP | 捏造 #1（方向反）|
| BRCA2 0.572 | design_03 HTML（已 rm）| 0.866 SEP | 捏造 #1（方向反）|
| 全基因組 unphase 44.89% | design_03 HTML（已 rm）| 45.84% | 捏造 #1 |
| null median 0.578 / BRCA2@p95 | VERIFIED_RESULTS.md V3 初版（已更正）| median 0.974 / @17.5pct | 捏造 #2（方向反）|

→ 兩次捏造共同點：**方向都與真值相反 + 都是「報告搶在分析回傳前憑預期填」**。postmortem: `InterSubMod/docs/postmortems/20260601_fabricated_metric_in_html_preview_postmortem.md`

---

## 待驗證隊列（下一步要產生的數字，現在還沒有 — 不准提前引用）

| 待測 | 為何要測 | 產生後寫入 |
|------|---------|----------|
| shuffle-label 控制（同 null 區打亂 HP 標籤重算 AUC）| 判 0.974 是真訊號還是方法 leak | 本檔 P-caveat → 升 P 或降 F-method |
| CpG-SNP pseudo-ASM 排除（剔除 C/T het SNP 改變 CpG 的位點）| lit workflow 列為頭號 artifact | 本檔 |
| 全基因組 unphase affinity AUC 分布 | 覆蓋率問題（45.8% 有多少可救）| 本檔 |
