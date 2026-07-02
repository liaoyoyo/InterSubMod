<!--
建立時間: 2026-06-27
類型: L6 跨樣本能力（參數化 + 文件化，不執行）— clone/subclone 整合報告
狀態: in_progress（HCC1395 ⭐3）
-->

# L6 — 跨樣本能力（參數化 + 資料就緒盤點）

> 用戶分階段裁決：① 先確 HCC1395 → ② 驗其他 HCC1395 資料 → ③ 同架構驗其他 ONT tag BAM。
> 本層**只參數化 + 文件化資料就緒**，不執行跨樣本（多數樣本資料缺）。

## §1 pipeline 參數化（任何 ONT tag BAM + VCF 可跑）
整條 pipeline 由 3 路徑驅動（`sm_linkage_genomewide.py` 頂部）：
- `VD` = per-chrom `filtered_snv_{tp,fp}_{chr}.vcf.gz` 目錄
- `TBAM` = tumor BAM（**需 longphase-S HP/PS tag**）
- `NBAM` = normal BAM
→ 換樣本只需換這 3 個 + CN/truth（可選）。**架構與方法不變**。

## §2 資料就緒盤點（實測 ls 確認）

| 樣本 | tumor BAM(tagged) | normal BAM(ONT) | TP/FP VCF | truth/CN | 可跑？ |
|---|---|---|---|---|---|
| **HCC1395（現行）** | ✅ Tmode_tagged_ClairS_v040 | ✅ HCC1395BL 5khz 5mC | ✅ pileup/ | ✅ SEQC2 | **已完成** |
| HCC1395 其他化學 | ⚠ Dorado / 5khz 原始存在（`/big8_disk/data/HCC1395/ONT_Dorado`、`ONT_5khz...`），**需 longphase-tag + 對應 VCF** | （同上 normal）| 需對應 | ✅ | **②within-sample 驗證：需先 tag Dorado BAM** |
| COLO829 | ✅ ONT R10/PAO | ❌ 只 Illumina（NYGC，非 ONT 甲基）| ❌ 缺 TP/FP VCF | ⚠ | ✗ blocked（VCF + ONT normal 甲基）|
| H1437/H2009/HCC1937/HCC1954 | ❌ 無 BAM（只 caller 輸出）| ❌ | ❌ | ❌ | ✗ blocked（無 BAM）|

## §3 分階段路徑（用戶裁決）
1. **① HCC1395 確認**：✅ 完成（骨幹 commit 0a8658d + L0-L4 整合，本報告夾）。
2. **② 其他 HCC1395 資料 within-sample 復現**：Dorado 化學存在 → 對 Dorado 跑 longphase-tag + 同 VCF → 同 pipeline 驗證骨幹/HP/CCF 是否化學間一致（**下一步可做，需 tag Dorado BAM**）。
3. **③ 其他 ONT tag BAM**：COLO829 待 TP/FP VCF（跑 caller）+ ONT normal 甲基；4 細胞株待 BAM 取得。皆**另案**。

## §4 結論
- **方法/架構已跨樣本可複用**（參數化），非 HCC1395-hardcoded。
- **現可執行的跨驗證 = HCC1395 Dorado 化學**（②）；其餘樣本（③）待資料齊。
- ⭐3 單樣本定論不變；跨樣本 ⭐4 需 ②/③ 補齊後另案。
