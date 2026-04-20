# C11: Phase 1A paired F1+0.011（結論 11）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 11
- **來源文件**: `04_暫停判定與重評估.md` / `06_結論穩定性審查.md:234-263`
- **穩定度評分**: ⭐3 需注意
- **狀態**: POSITIVE（弱，paired-only）
- **所屬證據鏈**: 🔗6（Phase 1A ML）
- **所屬樣本**: 7/7（6/7 正向, H2009 負向）

## D1 內部一致性
- ✅ External +0.0112 [+0.0044, +0.0188] / Discovery +0.0073 CI 含零 一致
- ✅ Mixed mode -0.0206 / H2009 -0.0039 一致
- ⚠️ McNemar p=1.44e-51 vs effect size 0.011 巨大反差（樣本量驅動）

## D2 方法論健康度
- ✅ External validation split（非 data snooping）
- ✅ McNemar 合適於 paired classifier 比較
- ❌ **Effect size 過小**（0.011 F1）- 是否臨床/實務顯著？C-STAT-3
- ❌ Discovery CI 含零 → 統計顯著但實務邊緣
- ⚠️ Mixed mode -0.0206 有害 → 適用範圍極窄
- ⚠️ H2009 負向根因未確認
- ❌ 無 per-sample CI overlap 揭露

## D3 證據鏈
- **依賴結論**: C01, C03, C10
- **被依賴結論**: 論文 variant filter 故事線（C-STRAT-1 認為過弱）
- **鏈完整度**: ⚠️ Mixed mode 破壞跨模式可攜性

## D4 數據信任度
- **dataset 版本**: Haplotag v5
- **CovM bug 影響**: 🟢 無（Phase 1A 不用 CovM baseline）
- **重跑必要性**: 無

## D5 統計嚴謹度
- **Effect size**: ΔF1 +0.0112（極小）
- **CI 覆蓋率**: ✅ External validation；❌ per-sample CI overlap 缺
- **Power 評估**: ✅ n 足夠偵測小效應
- **Multiple testing**: ⚠️ 特徵組合 × 超參掃描無 FDR

## D6 知識庫交叉驗證
- **支持文獻**: ML variant refinement 工具（Strelka2 etc.）F1 改善典型 0.01-0.03 → 0.0112 在範圍內
- **挑戰文獻**: ΔF1<0.02 通常不被 major journal 視為 significant gain
- **缺口**: 無

## 修正建議
- **P1**: 揭露 per-sample CI overlap + H2009 負向根因
- **P2**: 下修定位為「proof-of-concept，需更多外部樣本驗證」
- **對應 R-id**: R-02（論文定位已決定，但 C11 故事線需下調）
- **對應 Q-id**: Q-04（negative result 包裝可能性）

## 整體評分
**⚠️ POSITIVE 但弱 — 論文定位層級需重新思考；對應用戶新決策「Subclone Marker 為功能之一」時，C11 應降為次要支點。**
