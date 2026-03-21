<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# 甲基顯著計算：下一步目標與觀察計劃

## 目標與驗收
- 提升甲基顯著判斷對整體 F1 的貢獻，短期目標：在現有 F1 0.8155 基礎上，嘗試 ≥0.02 提升，同時 Recall 下滑不超過 0.005。
- 兩條路徑並行：**嚴格降 FP**（偏向 Precision）與 **寬鬆救 FN**（偏向 Recall），觀察 FP/TP 與 FN/FP 交換率是否符合「FP:TP < 1:0.69、FN:FP < 1:1.449」的有效區間。
- 所有程式碼實驗在獨立分支進行，紀錄變更點與指令；分析輸出統一放入 `docs/plans/2026/01/observations/`（新建可用）。

## 資料對齊與參考來源
- 基線資料：purity 1.0 的 HCC1395 tumor/normal，ClairS_pileup → longphase-s 強篩選 → HP tagged BAM/VCF。
- 交叉來源：
  - Gold truth：`data/answer/SEQC2/`
  - Phase 定位：`data/answer/phase_pos/`
  - ClairS 中間/最終輸出：`data/answer/ClairS_pileup/`、`data/answer/ClairS_output/`
  - 篩選後輸入：`data/vcf/HCC1395/pileup/`

## 觀察面向（TP/FP 與甲基特徵）
- 甲基化差異：TP vs FP 的 CpG methylation delta、entropy、missing rate；分 haplotype/strand 分層。
- 覆蓋/等位：TP vs FP 的 depth、AF、HP 一致性（HP 衝突率、phase block 長度）。
- 區域結構：距離矩陣稀疏度、DominantLabel 分布、Cramér’s V/LabelDelta 分位數。
- 品質/雜訊：MAPQ、base quality、indel 負載、strand bias 與顯著標記的聯動。

## 小範圍核對（不改碼）
- 單點檢視：`samtools view -h data/bam/HCC1395/tumor.bam chr:pos-window` + `bcftools view -r chr:pos-window data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz`，比對 MM/ML/HP/PS 與 VCF 取樣資訊。
- Python 快速檢查：用 pysam/pandas 計算指定位點的甲基化比例、AF、HP 分層後的 TP/FP 分布，輸出簡表存 `observations/`。

## 實驗步驟（可迭代）
1) **資料切片**：選取 TP/FP/FN 代表位點各 20–50 個；同步拉取 SEQC2 標籤與 phase_pos 定位。
2) **特徵切換試驗**（無程式改動）：在分析腳本增加篩選條件的虛擬開關（離線統計），例如調高/降低 Cramér’s V、LabelDelta 門檻，觀察 TP/FP 交換率。
3) **門檻搜尋**：用 Python 網格搜尋（如 V、p-value、最少 CpG/reads）計算 Precision/Recall/F1，確認是否有 ≥0.01 的可行提升區域。
4) **程式碼調整試驗**（獨立分支）：若門檻或統計方法需改動，在新分支實作、`ctest --test-dir build` + `./scripts/verify_output.sh`，並以 `compare_vcf_results.py` 或自製對照表量化變化。
5) **回歸驗證**：確保 baseline 全流程（`run_vcf_all_snv.sh --mode all-with-w1000 --plot-type no`）不退化；新增報告存 `observations/`。

## 開放觀察與後續
- 若嚴格路徑仍留 FP，考慮引入「phase-block 長度 + methylation entropy」複合閾值；若寬鬆路徑 Recall 不足，考慮放寬 low-depth 區段但加上 strand/HP 一致性補償。
- 若發現新可疑特徵（如特定 context 的錯報高發），在 `observations/` 加上案例與建議過濾規則，下一輪迭代再評估程式碼更動。
