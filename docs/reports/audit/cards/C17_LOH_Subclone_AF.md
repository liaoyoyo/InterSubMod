# C17: LOH Subclone AF×Methylation 雙證據鏈（補充結論 17）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 17（補充）
- **來源文件**: `07_LOH_CN_AF_研究總整理.md` / `06_結論穩定性審查.md:345-385`
- **穩定度評分**: ⭐4 穩固 ⬆️（2026-04-19 批次 2 升級）
- **狀態**: POSITIVE（雙模式確認）
- **所屬證據鏈**: 🔗9（Part B）+ 新鏈候選 🔗10
- **所屬樣本**: 7/7 step2；5/7 step3 明確 + 2/7 功效不足

## D1 內部一致性
- ✅ Δ_NG TO +0.705 / Paired +0.787 一致
- ✅ step3 meta-ρ TO +0.367 / Paired +0.357 一致
- ❌ **`20260417_Zone_Aware_Confidence_Framework_01.md:49,73` 誤標 Coverage_Multiple r=0.997**（實為 C20 r=0.831 OR C21 Jaccard=1.0000 OR LOH.bed SEQC2 0.962）→ **P1-A**
- ⚠️ HCC1954 n=34 step3 「反向」與「功效不足」兩種解讀並存，B.2-1 已確認為噪音

## D2 方法論健康度
- ❌ **Pooled OLS residualization 疑慮**（C-STAT-1）：Δ_NG Inter−Ext 跨樣本聚合是 pooled OLS 還是 within-sample？→ **P0-B**
- ✅ chr-shuffle null 已 pass（spatial 無 confound）
- ✅ NR-matched 7/7 POS（B.2-3 駁回 NR confound）
- ✅ cnLOH vs deletion-LOH 7/7 / 6/6 方向一致（B.2-5 駁回混合擔憂）
- ⚠️ HCC1954 n=34 underpowered 已標，但未傳遞至其他文件 → **P3**
- ❌ step3 meta-ρ 缺 per-sample bootstrap CI
- ❌ 無 FDR

## D3 證據鏈
- **依賴結論**: C01, C02, C09, C12, C15, C16, C21
- **被依賴結論**: C18（HCC1954 reversal 機制）、主軸 characterization
- **鏈完整度**: ✅ 完整（雙證據鏈：step2 AF + step3 segment ρ）

## D4 數據信任度
- **dataset 版本**: Haplotag v5 + loh_subclone_af_paired run
- **CovM bug 影響**: 🔴 **B.2-2 Coverage_Multiple proxy + CN1 only 限定依賴 CovM baseline** → hardcoded 75.0 bug 直接影響 CN1 identification → **R-01 已決定立即修+重跑**
- **重跑必要性**: **必需**

## D5 統計嚴謹度
- **Effect size**: Δ_NG +0.705（大）；meta-ρ +0.367（中）
- **CI 覆蓋率**: ✅ step3 leave-one-out sensitivity 已做；❌ step2 per-sample bootstrap 缺
- **Power 評估**: ⚠️ HCC1954/H2009 TO n<50 underpowered
- **Multiple testing**: ❌ 7 樣本 × step2/step3 × CN tier 無 FDR

## D6 知識庫交叉驗證
- **支持文獻**: epiTRACERx / EVOFLUx 等 subclone marker 文獻（knowledge 03, L2-L4）
- **挑戰文獻**: epiTRACERx 用不同 subclone marker；需比較 Δ_NG AF×Methylation 是否文獻空白
- **缺口**: 文獻空白（user memory 標為「文獻空白」）；Phase 2 A+D 補齊

## 修正建議
- **P0-A**: CovM bug 修 + 重跑（R-01 已決定）
- **P0-B**: within-group OLS 重算 Δ_NG
- **P1-A**: 修正 20260417_Zone_Aware_Confidence_Framework_01.md:49,73 r=0.997 數值錯標
- **P1-C**: step2/step3 per-sample bootstrap CI
- **P3**: HCC1954 underpowered 傳遞至 CURRENT_FOCUS / INDEX
- **對應 R-id**: R-01, R-02, R-05
- **對應 Q-id**: Q-10（臨床 downstream cohort）

## 整體評分
**✅ POSITIVE 雙證據鏈穩固，但 CovM bug 影響下游 CN1 分層 → R-01 修復後需 within-group OLS 重算 → 主軸 characterization 核心 Subclone Marker 功能。**
