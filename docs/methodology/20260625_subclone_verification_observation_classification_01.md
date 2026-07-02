<!--
建立時間: 2026-06-25
類型: methodology — Subclone 驗證研究 Phase 1：清楚觀察 + 詳細分類與定義（tumor-only / tumor-normal 雙軌）
狀態: in_progress（單樣本 HCC1395 ⭐2-3 觀察層；cis-clean 重標待 Phase 1.3）
build_branch: docs/method-comparison-ism-external-202606
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/corrected_tree.json, docs/methodology/_assets/20260618_subcluster_pilot/contingency_summary.json, docs/methodology/_assets/20260618_subcluster_pilot/percpg_summary.json, docs/methodology/_assets/20260618_subcluster_pilot/decisionflow_summary.json, docs/methodology/_assets/20260618_subcluster_pilot/tumor_vs_paired.json, docs/methodology/_assets/20260618_subcluster_pilot/fptp_attribution.json, docs/methodology/_assets/20260618_subcluster_pilot/percpg_cpg_classification.json, docs/methodology/_assets/20260618_subcluster_pilot/percpg_distance.json, docs/methodology/_assets/20260618_subcluster_pilot/phylo_cpp_wg_full_summary.json
-->

# Subclone 驗證 Phase 1 — 清楚觀察 + 詳細分類與定義（HCC1395，tumor-only / tumor-normal 雙軌）

> **任務 P1.1**：把 HCC1395 全位點的既有分類整理成一套**清楚、精確、grounded** 的觀察 + 定義，作為後續「驗證甲基狀況 + 長 read 多 sSNV」的乾淨基座。**數字全來自上列 data_sources JSON（這一輪 Read 回真值）。** 單樣本 HCC1395 ⭐2-3。

## §0 目的與執行原則

- **先觀察/分類定義 → 後驗證**（用戶 2026-06-25 定）。本文件 = 觀察 + 定義交付；驗證在 Phase 2。
- **🔴 tumor-only 與 tumor/normal 分開處理**：既有 `decisionflow_wg.py` 已實作雙軌（`-t tumor.bam -n normal.bam` 跑，每位點建 `tum`=只 is_tumor reads / `mer`=T+N 全 reads，各做 5-態分類）→ **可重用**。§4 直接用此雙軌展示 T/N confound。
- **兩個 universe**（勿混）：
  - **records_v6 universe = 34,736**（= TP 30,077 + FP 4,659；含 FP，per `phylo_cpp_wg_full_summary.json`）→ §1-§3, §5-§7。
  - **decisionflow universe = 30,490**（TP-only，跑 `filtered_snv_tp_*.vcf.gz`；precond tumor 30,077）→ §4。

## §1 主分類定義表（8 類，N=34,736，corrected_tree.json）

| 類 | 精確定義 | n | %總 | 主排除/歸因理由 |
|---|---|---|---|---|
| **D** 可測無結構 | 有軸(REF/ALT 或 HP)但**單群**（coarse=1） | 21,963 | 63.23% | 可測但無子結構 → 無 subclone 可言 |
| **C-S1** 有結構+有軸 | coarse≥2 **且對齊標籤**（cis-ASM-like，最乾淨子層） | 6,857 | 19.74% | 結構被既有標籤解釋 → cis-ASM 嫌疑 |
| **B1** 真單群-LOH | LOH 致生物性同質單群 | 1,978 | 5.69% | LOH unmask，非技術缺陷 → 不可分 |
| **C-S2** 有結構+有軸 | coarse≥2 對齊（次乾淨） | 1,569 | 4.52% | 同 C-S1 |
| **C-S3** 有結構+有軸 | coarse≥2 但 modal_frac<0.7 **不穩定** | 1,407 | 4.05% | 切法相依 → 存疑 |
| **🎯 A** subclone 候選 | **單標籤 + 多結構**（coarse≥2 但只一個 haplotag 標籤） | **723** | **2.08%** | 結構不被單一標籤解釋 → **候選（必要非充分）** |
| **B3** 真單群-非LOH足覆蓋 | 足覆蓋、非 LOH、生物同質單群 | 193 | 0.56% | 生物性單群 |
| **B2** 無法區分-低覆蓋 | 覆蓋不足無法切 | 46 | 0.13% | 技術性 |

**驗算 sum_check=true**（8 類加總 = 34,736）。**LOH+somatic 軸+結構**交集 highlight = 3,553（10.23%）= 好的「LOH 內出現 HP1/HP1-1 結構」說明例。

## §2 結構×標籤關係（contingency，n_multi_structured=10,556）

| 型 | 定義 | n | %多結構 | germline 雙家族 | somatic 單家族 | **單標籤-somatic** |
|---|---|---|---|---|---|---|
| **mixed** | 兼有多對一/一對多 | 3,980 | 37.7% | 2,691 | 1,289 | — |
| **many1** 結構>標籤 | 一標籤內多結構（**最像 subclone**） | 3,187 | 30.19% | 1,453 | 642 | **1,043** |
| **one2one** | 1 標籤 ↔ 1 結構（cis-ASM） | 1,739 | 16.47% | 1,024 | 715 | — |
| **1many** 跨標籤 | 一結構跨多標籤（無 ASM） | 1,650 | 15.63% | 1,552 | 98 | — |

🔑 **最乾淨 subclone 候選核 = many1 × single_label_somatic = 1,043**（一個 somatic 子單倍型標籤內出現多結構；single_label_germline 僅 49）。

## §3 per-CpG verdict + FP/TP 歸因（N=34,736）

**verdict 分佈**（percpg_summary.json）：no_cluster 24,180 / **cis-ASM(cluster≈標籤) 7,098** / **subclone_novel 523** / cluster_no_sigCpG 2,284 / weak_scattered 651。cis-ASM 細分：REF/ALT 3,804・tumor/normal 1,406・HP1/HP2 1,239・HP1/HP1-1 339・HP2/HP2-1 310。

**FP/TP 歸因（≥3 差異 CpG，fptp_attribution.json）— 哪些軸 TP 專一可做研究**：

| 軸/狀況 | %位點 | TP | FP | FP/TP | 機制 | 研究可用 |
|---|---|---|---|---|---|---|
| tumor/normal | 95.8% | 28,875 | 4,396 | 0.98 | 普遍 somatic 偏移，不專一 | ⚠ 描述可 |
| HP1/HP2 germline | 77.2% | 23,298 | 3,523 | 0.98 | germline cis-ASM，不專一 | ⚠ |
| **somatic (HP-1)** | **70.8%** | **23,414** | **1,192** | **0.33** | somatic 變異定義子單倍型 → **TP 專一（最佳 somatic 定位）** | ✅ |
| REF/ALT | 43.3% | 13,954 | 1,070 | 0.50 | TP 真 ALT；FP=germline/artifact | ✅ |
| cluster ≥3 | 23.4% | 7,011 | 1,118 | 1.03 | 子群結構略 FP | ✗ |
| **cluster ≥5** | 22.6% | 6,835 | 1,032 | 0.97 | **subclone-like cis/CN 驅動 → FP 富集** | ✗ 單樣本不可宣稱 |

🔴 **subclone（cluster）軸 FP 富集（≥5 ratio 0.97，越嚴越富 FP）→ 單樣本不可當 somatic subclone 偵測基礎**；唯一 TP 專一定位 = **somatic HP-1 軸（0.33×）**。

## §4 🔴 tumor-only vs merged 雙軌（5-態，decisionflow，TP-only N=30,490）

**5 態定義**：S1 不可驗證(n_complete<6)・S2 1群/無訊號・S3 監督可分 mean-shift(PERMANOVA-FDR 顯著；location vs dispersion)・S4 可分但**未對齊**(epiallele?雜訊?)・S5 **確認真結構**(離散切得出且對齊生物軸)。

| 態 | tumor-only (只 tumor reads) | % | merged (T+N reads) | % |
|---|---|---|---|---|
| S1 不可驗證 | 413 | 1.35 | 377 | 1.24 |
| S2 1群/無訊號 | 5,051 | 16.57 | 715 | 2.35 |
| S3 監督可分 | 5,525 | 18.12 | 2,061 | 6.76 |
| S4 可分未對齊 | 13,390 | 43.92 | 16,748 | 54.93 |
| **S5 確認真結構** | **6,111** | **20.04** | **10,589** | **34.73** |
| split_align_rate | — | 31.34% | — | 38.74% |

🔑 **merged（T+N）的 S4+S5（結構）= 89.66% 遠高於 tumor-only 的 63.96%** → **加入 normal reads 把結構灌大**（T-vs-N 組成本身造結構）= **直接證實混 T+N 的 confound**。**tumor-only S5=20.04% 才是 tumor 自身內在結構的乾淨視角**。佐證 `tumor_vs_paired.json`：tumor-only 顯著 4,185 vs paired-allele 顯著 22,632（混入 normal 大幅膨脹）。→ 用戶「分開處理」的直覺正確，後續候選一律以 **tumor-only 軌** 為觀察基座、normal 只作獨立 cis baseline。

## §5 差異 CpG 方向 / 強度 / 距 SNV 分佈（percpg_cpg_classification + percpg_distance）

- **方向（hyper% = group1 增甲基佔比）**：tumor/normal 90.1%・HP1/HP1-1 85.3%・HP2/HP2-1 85.4%（**somatic 軸增甲基 gain 偏向**）；germline/cluster 平衡（HP1/HP2 50.5%・REF/ALT 46.7%・cluster 42.9%）。
- **總差異 CpG instance = 2,819,893**（跨軸）；tumor/normal 軸最多（1,240,666）。
- **距 SNV 分佈**：背景 CpG 0-500bp 佔 10.4%；各軸差異 CpG **near_snv_enrich ≈ 1.0-1.15**（REF/ALT 1.12、cluster 1.15 略高）→ **差異甲基非點狀富集於 SNV，呈區域性分佈**（與「regional partition」一致，非單點 ASM）。

## §6 生物檢核（corrected_tree.json）

**HP1-1 ⟹ ALT 一致性 = 99.36%**（34,372 / 34,595 帶 somatic 子單倍型的位點同時有 ALT 差異）→ 證實「存在 HP1/HP1-1 必伴 REF/ALT 差異」的知識背景，分類自洽。

## §7 研究可用性裁決（觀察層）

- ✅ **可做研究的 TP 專一定位 = somatic HP-1 軸（70.8%，0.33×）+ REF/ALT（43.3%，0.50×）** → 標記「哪些 read 帶 somatic 變異 + 其甲基」。
- 🎯 **subclone 候選核 = A 類 723 / many1×single_label_somatic 1,043 / subclone_novel 523** → 但 cluster 軸整體 FP 富集，**須 normal cis-test 清洗**才知真候選數。
- ✗ **cluster/subclone 軸（23.4%）FP 富集** → 單樣本不可當 somatic subclone 偵測。
- ⚠ tumor/normal（95.8%）+ germline HP（77.2%）普遍但不專一 → 描述可、判別不可。

## §8 已知限制 + 下一步

1. **🔴 cis-clean 未做**：subclone 候選（A/subclone_novel/many1）尚未扣 normal cis → 不知多少是 germline cis-residual。**→ Phase 1.3 normal-anchored cis-attribution（HCC1395 normal ready）**，重標 {cis-ASM / 候選-subclone(residual) / 無法判定}。
2. **per_cpg_diff 混 T+N**：除 tumor_vs_normal 軸外，per-CpG 分組用全 reads（confound）→ **Phase 1.2 Track A 須加 is_tumor filter** 重算 tumor-only per-CpG（decisionflow 結構層已是 tumor-only，per-CpG 層待補）。
3. **多 sSNV 連鎖未做**：唯一非循環錨（Phase 2.2/2.3）。
4. strand 未分（5mC/5hmC 混）；單樣本無法宣稱 somatic subclone（需 normal/multi-region/single-cell）。

## §9 Provenance + 自驗

| 數字群 | 來源檔（grep 可得） |
|---|---|
| 8 類 + 生物檢核 | `corrected_tree.json` |
| contingency 型 | `contingency_summary.json` |
| verdict 分佈 | `percpg_summary.json` |
| 5-態雙軌 | `decisionflow_summary.json` |
| tumor vs paired | `tumor_vs_paired.json` |
| FP/TP 歸因 | `fptp_attribution.json` |
| 差異 CpG 方向/強度 | `percpg_cpg_classification.json` |
| 距 SNV 分佈 | `percpg_distance.json` |
| TP/FP 總數 | `phylo_cpp_wg_full_summary.json` |

> 紅線：⭐3 單樣本 HCC1395；本文件為**觀察 + 分類定義**（pre-cis-clean），subclone 候選為「必要非充分」；cluster 軸 FP 富集非 somatic-specific；唯一 TP 專一 = somatic HP-1 軸。cis-clean 重標 + 多 sSNV 連鎖驗證見後續 Phase。§13 數字不手打，全 JSON 溯源。
