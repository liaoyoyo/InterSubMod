# C18: HCC1954 CNV-driven reversal 機制（補充結論 18）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 18（補充，2026-04-19 新增）
- **來源文件**: `06_結論穩定性審查.md:388-403` / `20260419_Batch2_cross_track_synthesis_01.md`
- **穩定度評分**: ⭐4 穩固
- **狀態**: CONFIRMED（雙獨立路徑）
- **所屬證據鏈**: 🔗10（HCC1954 HER2+/pseudo-tetraploid CNV 新鏈候選）
- **所屬樣本**: HCC1954-specific

## D1 內部一致性
- ✅ frac_intermediate=0.727 / 其他 0.118-0.383 一致
- ✅ Z3 chr5/8/17 arm-level 85% FP 一致
- ✅ HCC1954-only blacklist ΔF1=+0.0065 一致

## D2 方法論健康度
- ✅ 雙獨立路徑（B.2-1 AF 分佈 + Z3 chromosome enrichment）
- ✅ 生物機制對接（HER2+/chr17q amp / chr8p loss / TP53 LOH / chr5 rearrangement）
- ⚠️ **單樣本機制泛化問題**：僅 HCC1954 ONT 數據，其他 HER2+ cell lines 無獨立重現
- ⚠️ Coverage_Multiple 受 CovM bug 影響 → 但 per-sample 特徵仍存在
- ❌ 無 Normal BAM 驗證（R17 批次 3 候選）

## D3 證據鏈
- **依賴結論**: C02, C17, C19, C20
- **被依賴結論**: C22（Zone Characterization 樣本特異性）
- **鏈完整度**: ⚠️ 單樣本雙獨立路徑 → 機制合理但泛化待驗

## D4 數據信任度
- **dataset 版本**: Haplotag v5 + B.2 批次 1 + Z3 pilot
- **CovM bug 影響**: 🟡 HCC1954 Z3 FP CovM=0.733 per-sample 層面明顯，但 bug 修正後絕對值會變（相對排序應穩定）
- **重跑必要性**: CovM bug 修後重算 per-sample z-score

## D5 統計嚴謹度
- **Effect size**: 單樣本效應大（frac_intermediate 外 2σ）
- **CI 覆蓋率**: ⚠️ 單樣本無 CI
- **Power 評估**: N/A（機制歸因）
- **Multiple testing**: N/A

## D6 知識庫交叉驗證
- **支持文獻**: HER2+ breast cancer molecular pathology（knowledge 02）；chr8p loss + chr17q amp 為典型
- **挑戰文獻**: 無（HCC1954 pathology 標準文獻支持）
- **缺口**: patient-derived HER2+ ONT 資料 / 其他 HER2+ cell lines

## 修正建議
- **P1**: R17 批次 3 HCC1954 Normal BAM pilot
- **P1**: 其他 HER2+ cell lines 獨立驗證
- **P1**: CovM bug 修後 z-score 重算確認 per-sample 特徵穩定
- **對應 R-id**: R-01
- **對應 Q-id**: Q-07（單樣本外推）、Q-10（臨床 HER2+ cohort）

## 整體評分
**✅ 雙路徑機制穩固 — 需跨 HER2+ 樣本泛化驗證，對 R-01 CovM 修正相對 robust。**
