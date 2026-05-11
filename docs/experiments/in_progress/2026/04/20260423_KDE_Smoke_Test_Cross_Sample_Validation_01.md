---
title: KDE 自動估算 smoke test 跨樣本驗證（chr19 subset, 13 runs）
date: 2026-04-23
status: SMOKE_VALIDATED（chr19 subset 通過；全基因體重跑待決策）
hypothesis_id: H_KDE_001 (extension)
verdict: KDE 邏輯跨 7 樣本 × 2 modes 一致啟用；master dataset 75.0 default 對所有樣本偏差顯著
priority: P0
estimated_effort: chr19 13 runs = 5.3 min（完成）；全基因體估 60-120 min（待決策）
related:
  - docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md
  - docs/experiments/in_progress/2026/04/20260420_CovM_Baseline_Accuracy_Validation_01.md
  - docs/methodology/20260419_KDE_expected_coverage_audit_01.md
  - docs/CURRENT_FOCUS.md
batch_output: /big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/batch_20260422_234255/
cpp_binary:
  path: /big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod
  mtime: 2026-04-21 20:01
  commits_included: ["374fad4", "12d9b3e", "5abc659"]
rollback_anchor: "pre-kde-logging-audit (commit ec0608b)"
---

# KDE 跨樣本 smoke test 驗證報告

> **結論**：13 / 13 runs KDE `auto_kde` 正常啟用，跨 7 樣本 × 2 modes diploid coverage range **29x–111x（3.8×）**；master dataset 75.0 default 對所有樣本均偏差，COLO829 高估 159%、HCC1937 低估 32%；Paired 與 TO 在同樣本 KDE 結果 gap ≤3%（6/6 有 TO 的樣本），驗證 tag.bam 深度在兩模式間等價。chr19 subset 不等於全基因體 baseline（HCC1395 兩次 smoke 123x vs 65x 差 1.9×），**最終 per-sample baseline 需全基因體重跑**。

---

## 1. 背景（為何跑這份驗證）

### 1.1 前置狀態

- 2026-04-19 commit `374fad4+12d9b3e+5abc659`：C++ 新增 `Diploid_Coverage_Used` TSV 欄位 + KDE fallback WARN + DiploidEstimate struct
- 2026-04-20 HCC1395 TO pilot 於單樣本驗收通過（`20260420_KDE_Fix_Acceptance_Validation_01.md`，dc=53.0x）
- 2026-04-20 CovM baseline 驗證揭露 master dataset 問題：全 7 樣本共用 75.0 default（`20260420_CovM_Baseline_Accuracy_Validation_01.md`）
- 2026-04-22 盤點確認：`scripts/run_vcf_all_snv.sh:301` 硬編碼指向 big8 舊 binary（2026-03-06，無 KDE）→ 這是 master dataset 75.0 的真正根因

### 1.2 本次驗證目的

回答三個問題：
1. 新 binary 是否跨所有 7 樣本 × 2 modes 都能正確啟用 KDE？
2. 不同樣本的真實 per-sample baseline range 為何？Master 75.0 default 的偏差有多大？
3. Paired 與 TO 模式在同樣本的 KDE 結果是否一致？

---

## 2. 方法

### 2.1 Binary

```
/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod
mtime: 2026-04-21 20:01, size 6.8 MB
```

CLI：`--expected-coverage 0` 觸發 KDE auto-estimation。

### 2.2 執行範圍

| 維度 | 選擇 | 理由 |
|------|------|------|
| 樣本 | 7 細胞株 | COLO829 / H1437 / H2009 / HCC1395 / HCC1395_DORADO / HCC1937 / HCC1954 |
| 模式 | Paired + TO | 各樣本同時測兩模式 |
| 區域 | chr19 only | 快速 smoke，非全基因體 |
| VCF 子集 | filtered_snv_tp 的 chr19 variants | 每 run 322-2518 variants |
| Normal BAM | 不提供 | 統一 conditions，僅驗證 KDE 邏輯；SampleASM 欄位非本次重點 |
| 其他參數 | `--window-size 5000 --threads 16 --distance-metric BERNOULLI` | 與既有 master dataset 配置一致 |

### 2.3 Tag.bam 來源

- Paired：`canonical/{sample}/paired_full/20260314-20260315_*/longphase_s/{sample}_tagged.bam`
- TO：`synthesis/research_rounds/archive/202603_early_pilots/{date}_{sample}_to_pilot*/step03_longphase_to/tumor_tagged.bam`
- COLO829 TO：**缺**（step03 空，需獨立 longphase-to 重跑）

### 2.4 批次執行

Script：`/big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/batch_20260422_234255/run_batch.sh`

對每個 (sample, mode) 執行：
1. `bcftools view -r chr19` 產 chr19 subset VCF
2. 呼叫 binary `--expected-coverage 0`
3. 解析 `run.log` 取 `diploid_coverage` 與 `source=...`
4. 寫入 `kde_summary.tsv`

---

## 3. 結果

### 3.1 Per-sample KDE 結果表

| 樣本 | chr19 variants | Paired DC | TO DC | gap | vs master 75.0 |
|------|:--------------:|:---------:|:-----:|:---:|:--------------:|
| COLO829 | 750 | **29.0x** | N/A | — | **+159% 高估** 🔴 |
| H1437 | 503 / 480 | **43.0x** | **43.0x** | 0 | **+74% 高估** 🔴 |
| HCC1954 | 294 / 284 | **57.0x** | **57.0x** | 0 | +32% 高估 🟡 |
| HCC1395 | 650 / 615 | **65.0x** | **63.0x** | +3.2% | +15% 高估 |
| HCC1395_DORADO | 645 / 623 | **75.0x** | **75.0x** | 0 | 0%（巧合）|
| H2009 | 2,518 / 2,393 | **83.0x** | **83.0x** | 0 | -10% 低估 |
| HCC1937 | 322 / 332 | **111.0x** | **111.0x** | 0 | **-32% 低估** 🔴 |

### 3.2 核心發現

**發現 1：KDE 跨樣本全部成功啟用**
- 13/13 runs `source=auto_kde`（無 fallback WARN）
- 無 `histogram range too narrow` / `insufficient valid regions` 觸發
- 5.3 min 完成全部（單 run 10-91 秒）

**發現 2：Paired vs TO 高度一致**
- 6 / 6 有 TO 資料的樣本，gap ≤3%（5 / 6 完全一致 dc=相同，HCC1395 僅差 +3.2%）
- 此結果支持「tag.bam 深度在 paired（longphase-s）與 TO（longphase-to）處理後差異極小」
- 意味：現有 master dataset 全 75.0 同時錯在兩個模式（而非 paired 對、TO 錯或反之）

**發現 3：master 75.0 default 偏差跨樣本巨大**
- Range：29x（COLO829）→ 111x（HCC1937），**3.8× fold**
- COLO829 被高估 159%（真實 29 → master 75）
- HCC1937 被低估 32%（真實 111 → master 75）
- 7 樣本中 6 個 bias >10%，僅 HCC1395_DORADO 恰好 75x（純巧合）

**發現 4：chr19 subset ≠ 全基因體 baseline**
- 首次 smoke test（2026-04-22）HCC1395 paired dc = **123x**
  - VCF = `/big8/.../filtered_snv_tp_chr19.vcf.gz`（672 variants）
  - 有給 `--normal-bam`
- 本次 batch HCC1395 paired dc = **65x**
  - VCF = `canonical/.../filtered_snv_tp.vcf.gz` 經 bcftools chr19 subset（650 variants）
  - 未給 normal BAM
- 差異 1.9× → KDE 對 VCF 子集與 normal BAM 設定敏感

### 3.3 執行效能

| 指標 | 值 |
|------|-----|
| 總 wall time | 317 秒（5.3 min） |
| 單 run range | 10-91 秒 |
| 最慢 run | H2009 paired（91 秒，2,518 variants） |
| 最快 run | H1437 paired（10 秒，503 variants） |
| Throughput | ~35 regions/秒（16 threads） |

---

## 4. 與 20260420 pilot 的差異

| 面向 | 20260420 pilot | 20260423 本次 |
|------|----------------|--------------|
| 範圍 | 單樣本（HCC1395 TO） | 跨 7 樣本 × 2 modes |
| 區域 | 全基因體 28,495 regions | chr19 subset（322-2,518 regions/run）|
| HCC1395 TO dc | 53.0x（全基因體）| 63.0x（chr19 subset）|
| 差異成因 | chr19 子集 TP 富集於高覆蓋區 | — |

**補充**：HCC1395 TO 全基因體 dc=53，chr19 subset dc=63，差異 +19% — 說明 chr19 的 TP 區域平均覆蓋度略高於 genome-wide neutral。

---

## 5. 限制與可能誤判

### 5.1 已知限制

- **L1**：本次 chr19 subset 結果**不可直接用作 master dataset CovM baseline**
- **L2**：COLO829 TO 缺（step03_longphase_to 空），TO 跨樣本對照不完整
- **L3**：Normal BAM 未提供 → SampleASM / Normal baseline 欄位本次無法驗證
- **L4**：不同 VCF subset（TP only vs TP+FP vs 全 VCF）可能產生不同 KDE 估值，本次僅測 TP

### 5.2 可能反例（需後續驗證）

- HCC1395_DORADO 剛好 dc=75 是否為 KDE 特殊偏誤（而非真實 baseline）？建議全基因體獨立驗證
- HCC1937 dc=111 顯著高於其他樣本，需核對 tag.bam 覆蓋度是否異常（read length / mapq filter 影響）

---

## 6. 結論

1. ✅ **H_KDE_001 擴展驗證通過**：KDE 邏輯在 7 樣本 × 2 modes 全部正常；原 20260420 的 HCC1395 單樣本結論可擴展至跨樣本
2. ✅ **Paired/TO KDE 等價假設成立**：tag.bam 在 longphase-s 與 longphase-to 處理後 read density 分布幾乎相同（gap ≤3%）
3. 🔴 **Master dataset 75.0 default 必須重跑**：跨樣本 baseline range 29x-111x，6/7 樣本 bias >10%；**所有使用 master CovM 的跨樣本分析（H-CN1/H-CN2、HPFineNGroups CN tiers、Zone-Aware 覆蓋率、LOH×AF×TP:FP 等）結論需在重跑後重算**
4. ⚠️ **chr19 subset 不可取代全基因體**：HCC1395 paired chr19 不同 VCF 差 1.9×；最終 per-sample baseline 需全基因體重跑

---

## 7. 後續動作（依優先級）

### P0（阻塞 P1 分析）

- [ ] **H1**：13 runs 全基因體重跑（使用新 binary + canonical VCF 全 autosomes），~60-120 min
- [ ] 重建 master dataset，驗證反推 14 rows `ec` 各自獨立（不再共用 75）
- [ ] 修改 `scripts/run_vcf_all_snv.sh:301` EXECUTABLE 指向 big7 新 binary（獨立 `/cpp-change` commit）
- [ ] COLO829 TO 補跑 longphase-to（取得 step03 tagged BAM，~30 min），再跑 ISM

### P1（KDE fix 連動重量化）

- [ ] H-CN1 / H-CN2 per-sample baseline 重算（`20260420_CovM_Baseline_Accuracy_Validation_01.md`）
- [ ] HPFineNGroups F pilot CN 跨 tiers 0.90-0.94 穩定性重驗證
- [ ] Zone-Aware Z1-Z5 覆蓋率重算
- [ ] LOH × AF × TP:FP 比例與生物意義重分析（配合新 per-sample baseline）

### P2（診斷已知 anomaly）

- [ ] HCC1395 paired chr19 兩次 smoke 1.9× 差異根因（VCF subset size / normal BAM 影響分析）
- [ ] HCC1395_DORADO dc=75 巧合是否為 KDE 偏誤
- [ ] HCC1937 dc=111 高覆蓋是否合理（核對 tag.bam QC）

---

## 8. 相關檔案索引

### 輸出資料

- 批次結果：`/big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/batch_20260422_234255/`
  - `run_configs.tsv` — 13 runs 設定
  - `run_batch.sh` — 批次執行 script
  - `kde_summary.tsv` — 關鍵結果彙總
  - `batch.log` — 每 run 狀態
  - `{sample}_{mode}/` — 各 run 獨立輸出目錄

### 先前驗證

- `docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md` — 單樣本 pilot
- `docs/experiments/in_progress/2026/04/20260420_CovM_Baseline_Accuracy_Validation_01.md` — H-CN1 / H-CN2 MIXED 結論
- `docs/methodology/20260419_KDE_expected_coverage_audit_01.md` — 方法論審計
- `research/autoresearch/evidence_ledger.jsonl` — H_KDE_001 entry

### C++ 程式碼

- `src/core/RegionProcessor.cpp:681-691` — KDE Pass 2 分支邏輯
- `src/core/RegionProcessor.cpp:78-160` — `estimate_diploid_coverage()` 實作

---

## 9. Evidence ledger 更新建議

在 `research/autoresearch/evidence_ledger.jsonl` 新增條目：

```json
{
  "id": "H_KDE_001.v2",
  "date": "2026-04-23",
  "status": "SMOKE_VALIDATED_CROSS_SAMPLE",
  "parent": "H_KDE_001",
  "summary": "KDE auto_kde 跨 7 樣本 × 2 modes 全部啟用；chr19 subset range 29x-111x；master 75.0 default 6/7 樣本 bias >10%；全基因體重跑待決策",
  "evidence_path": "docs/experiments/in_progress/2026/04/20260423_KDE_Smoke_Test_Cross_Sample_Validation_01.md",
  "batch_data": "/big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/batch_20260422_234255/"
}
```
