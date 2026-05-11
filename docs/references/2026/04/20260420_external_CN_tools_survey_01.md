<!--
建立時間: 2026-04-20
更新時間: 2026-04-20
狀態: active
資料來源:
  - https://github.com/cortes-ciriano-lab/savana
  - https://github.com/KolmogorovLab/Wakhan
  - research/clairs_to_verdict_pilot（內部關聯）
-->

# External CN Tools Survey — Wakhan / SAVANA

## 動機

ClairS-TO Verdict Characterization Pilot（2026-04-20）結論：**Verdict 在 HCC1395 subsample t20_n30 內部校準正確（H-V1/H-V2 POSITIVE）**，但 Verdict_Germline 僅在 LowQual 變異上觸發（從未出現在 PASS）→ **啟用 Verdict 對 ClairS-TO 現有 PASS 輸出的 F1 直接增益 = 0**。原因：Verdict 與 LowQual 都依賴同一組 ASCAT 內部特徵，資訊重疊。

**解法候選**：改用**獨立於 ClairS-TO 的外部 CN 工具**生成 sample-level CN segmentation，再在 variant annotation 層面交叉標註，讓 CN-aware filter 從 **ClairS-TO 正交訊號** 進入決策。兩個主要候選：**Wakhan**（Kolmogorov Lab · medRxiv 2025）與 **SAVANA**（Cortes-Ciriano Lab · Nature Methods 2025）。

## 工具對照

| 特性 | Wakhan | SAVANA |
|------|--------|--------|
| **來源** | KolmogorovLab/Wakhan (v0.4.2, MIT) | cortes-ciriano-lab/savana |
| **發表** | medRxiv preprint (2025) | Nature publication reference (2025) |
| **Tumor-only 支援** | ✅（需 tumor phased VCF） | ✅（`savana to` 子命令） |
| **Paired mode** | ✅（tumor BAM + normal phased VCF） | ✅（tumor + normal BAM/CRAM，推薦使用） |
| **輸出** | haplotype-specific CN 段，purity，ploidy，interactive HTML，BED + VCF | SV (VCF/BEDPE)，CNA (CBS 段)，absolute CN，purity，ploidy |
| **Haplotype-specific CN** | ✅（核心能力） | ❌（single-scale CN） |
| **SV 支援** | 可選（讀取外部 SV） | ✅（內建 SV caller） |
| **ONT long-read** | ✅（建議 Clair3 + LongPhase phasing） | ✅（測試 ONT + PacBio HiFi + minimap2/winnowmap） |
| **Runtime（參考）** | 未明列 | COLO829 subset (50k reads) 約 49 秒（8 CPU, 16GB RAM）|
| **HCC1395 benchmark（README）** | 未明列（指向 medRxiv 預印本） | 未明列（指向 Nature 論文） |

## 評估：是否適合接入 InterSubMod

### Wakhan

**優勢**：
- **haplotype-specific CN** 正好與 InterSubMod 的 HP-tagged read-level ISM 結構對應，可整合進 Zone-Aware Framework 作 per-haplotype zone 定義
- tumor-only 模式成熟，可直接套用於 6 個缺 matched normal 的 in-house 樣本
- 需要 Clair3 phased VCF → 本專案已用 LongPhase + Clair3 產出，輸入準備成本低

**風險**：
- HCC1395 benchmark 未在 README 列明 → 需自行驗證（~2 day pilot）
- 純度/ploidy 收斂度未知，可能需手動調參

**預估 pilot 成本**：1 人 2-3 天（單樣本 HCC1395 全量跑 + 與 SEQC2 CNV truth 交叉）

### SAVANA

**優勢**：
- 同時產出 SV + CNA 雙輸出 → 可直接取代目前缺失的 SV caller 層
- Random Forest 內建 ONT 分類器，reproducibility 佳
- Nature Methods 2025 有明確 peer-reviewed 驗證

**風險**：
- 不提供 haplotype-specific CN → **對 HP-tagged ISM 整合價值低於 Wakhan**
- `savana to` 為 tumor-only 但官方「strongly recommend」paired — 本專案 6 樣本 TO 場景可能落在次優區

**預估 pilot 成本**：1 人 3-4 天（需先處理 SV/CN 雙路並行並與既有 pipeline 對接）

## 建議路徑（若 Verdict pilot 判定 NEGATIVE）

1. **優先 Wakhan**：haplotype-specific CN 與 InterSubMod HP-aware ISM 架構高度互補
2. **HCC1395 單樣本 pilot**（2-3 天）：
   - 輸入：tumor BAM + Clair3+LongPhase phased VCF
   - 產出：per-haplotype CN segments + purity + ploidy
   - 驗證口徑：Wakhan CN ∩ SEQC2 CNV gain/loss/LOH truth → 計算 region-level concordance
3. **若 HCC1395 通過**：橫擴 6 個 in-house 樣本
4. **可選 parallel**：SAVANA 作為 SV 互補層，與 Wakhan CN 合併
5. **整合層**：Wakhan CN 輸出轉成 BED → 加入 RegionProcessor 的外部 annotation 輸入；不修改 Verdict path

## 不建議的做法

- **改繞 ClairS-TO purity safety gate + 補 CNA loci**：成本高（Docker image 重建 + 1000G loci 下載 + ASCAT 高純度行為調查），收益已被本 pilot Step 2 量化為 0
- **直接用 Wakhan/SAVANA output 作 hard filter**：先做 characterization annotation；符合本專案「characterization before filter」慣例

## 參考連結

- **Wakhan GitHub**：<https://github.com/KolmogorovLab/Wakhan>
- **Wakhan medRxiv**（從 README 所指）：待本次 pilot 後再檢索精確連結
- **SAVANA GitHub**：<https://github.com/cortes-ciriano-lab/savana>
- **SAVANA Nature Methods 2025**：同上
- **ClairS-TO Verdict spec**：Knowledge base `03_file_formats/vcf-clairs-to.md`

## 下一步

本 survey 是 trigger document。若 Verdict pilot POSITIVE → 收入 archive 作為備選。若 Verdict pilot NEGATIVE（如實際結果）→ 建立獨立 pilot 任務 `research/wakhan_pilot/`。
