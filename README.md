# InterSubMod

InterSubMod 是一個專為長讀取 (Long-read) 測序設計的表觀遺傳變異分析工具。它能夠整合甲基化 (Methylation)、體細胞變異 (Somatic SNVs) 與單倍體型 (Haplotypes)，藉此解析腫瘤內的亞克隆結構 (Subclonal Structure)。

---

## 📚 專案文件導引 (Documentation)

請根據您的需求參考以下文件：

| 文件 | 描述 |
| :--- | :--- |
| **[🚀 快速開始 (Quick Start)](QUICKSTART.md)** | **推薦優先閱讀**。包含環境配置、編譯步驟，以及如何使用自動化腳本執行分析。 |
| **[📖 專案全貌 (Project Summary)](README_PROJECT_SUMMARY.md)** | 詳細的專案技術總結。包含完整分析流程、核心模組架構、演算法說明與技術亮點。 |

---

## ✨ 核心功能 (Key Features)

*   **亞克隆解析 (Subclonal Resolution)**: 利用 Read-level 甲基化模式，區分不同的細胞群體。
*   **多樣化距離度量**: 支援 L1、NHD 等多種距離算法，精確量化表觀遺傳異質性。
*   **視覺化整合**: 自動生成距離熱圖 (Distance Heatmap) 與分群熱圖 (Cluster Heatmap)。
*   **高效能運算**: 基於 C++17 與 OpenMP 平行化架構，快速處理大規模測序數據。
*   **精準位點映射**: 修正了 Indel 對甲基化座標的影響，確保每個 CpG 位點的精確對齊。

---

## 🛠️ 快速執行範例

編譯完成後，即可使用以下指令進行標準全流程測試：

```bash
./scripts/run_full_vcf_test.sh --mode all-with-w1000
```

詳細參數與使用方式請參閱 [Quick Start](QUICKSTART.md)。
