Development Spec

Scope: 目前不實作 calibration/ML；但需保留未來可用的輸出欄位與格式，確保後續能直接接 ML/校準流程。

Language & Ownership: 所有計算與統計（Fisher/Chi-square/PERMANOVA/Dispersion/Bootstrap/BH/效應量/彙整判定）必須由 C++ Core 完成；Python 僅允許視覺化或報表格式化，不得改變統計結論。

Parallelization: 可平行則平行（region-level、距離矩陣、permutation/MC）；需明確規範執行緒層級，避免巢狀 OMP 過度訂閱；資料結構生命週期必須清楚「建立→使用→釋放」，並保證記憶體用量合理。

OOP & Modularity: 統計與顯著性分析需拆分為獨立模組（Statistics/Significance/StructureTest/Dispersion/Bootstrap/Report），以單一入口整合到 RegionProcessor；檔案與函式分工清楚、可測試。

Comments: 重要演算法與統計定義必須有完整英文註解（two-sided Fisher 定義、MC 停止規則、PERMANOVA/Dispersion、邊界情況與數值穩定性）。

Processing Modes: 提供兩種模式且可切換：Tumor-only（只分析 tumor reads）與 Mixed（Tumor/Normal 共同分析）；Mixed 模式需保留 sample label，並允許檢驗「是否位點顯著與 Tumor/Normal 無關」。

Method Toggles: 每個顯著性方法可獨立啟用或同時啟用；必須有一個統一的 SignificanceAggregator 整合所有方法輸出與最終判定。

Reproducibility: 所有隨機流程（MC Fisher、permutation、bootstrap）必須可設定 seed，並寫入輸出。

Label Definitions: 標籤至少包含 Allele(ALT/REF/UNKNOWN), HP(HP1/HP2/HP1-1/HP2-1/HP3/unphase), Strand(F/R), Sample(Tumor/Normal)；HP 以 longphase-s 規則為準。

Dynamic Encoding: 以「類別相乘」方式形成組合碼：HP * Allele * Strand * Sample；需輸出 code 對照表與版本號，確保跨版本可追溯。

Composability: 必須能從組合碼投影回單一維度（例如只看 ALT/REF 或只看 HP）以支援不同檢定與解釋。

PERMANOVA/Dispersion Preprocessing: 若要把 PERMANOVA/Dispersion 作為嚴格證據，必須先過濾 reads 使距離矩陣完整；SKIP/MAX_DIST 僅可用於探索或視覺化，不得作為正式顯著性依據。

Insufficient N: 若子集合 N 太小（例如 < 5）或有效重疊不足，輸出 valid=false 並 reason=insufficient_overlap；可只保留簡單效應量或直接跳過檢定。

Decision Rule (Default): 最終「顯著/不顯著」需結合 p-value、q-value 與整合分數；預設以常見門檻（例如 p≤0.05、q≤0.1）並設定可調的 score_threshold，所有門檻必須可由參數調整以利後續驗證。

Score Definition: 整合分數至少由 -log10(p_global) 與結構性證據（PERMANOVA/Dispersion 通過時才納入）組成；若方法未啟用則自動移除該項並重新正規化權重。

Per-region Output: 每個位點資料夾下新增 significance/，至少包含：summary.(tsv|json)（顯著判定、p/q/score、valid/reason、方法啟用狀態）、tables/（contingency/cluster 統計）、metadata.(tsv|json)（seed、門檻、矩陣有效比例）。

Run-level Output: 輸出 significance_report.csv（每個 VCF 位點統計值與顯著結果）、significance_features.csv（保留給未來校準/ML 的特徵欄位）、label_codebook.tsv（組合碼對照）。

Report Summary: 產出整體統計說明文件，包含每染色體與全域的顯著比例、無效原因分布、方法覆蓋率與效能摘要。

Validation: 每個主要功能完成即寫單元測試，編譯並通過自動測試；測試結果記錄於 YYYYMMDD_測試內容.md。

Integration Test: 先跑單點驗證 run_full_vcf_test.sh --mode chr19-verification --no-plots，確認流程可完整執行並產生合理顯著結果。

Performance Target: 30,000 點在 120 執行緒下約 10 分鐘完成；若未達標需紀錄瓶頸與調整策略（如 gating/減少 permutation 次數）。