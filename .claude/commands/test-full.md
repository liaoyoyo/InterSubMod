---
allowed-tools: Bash(./scripts/run_batch_vcf_analysis.sh:*)
description: 執行完整流程測試 (約 5 分鐘)
---

執行完整流程測試，包含所有輸出和圖表生成。

```bash
./scripts/run_batch_vcf_analysis.sh
```

測試完成後：
1. 確認所有輸出無錯誤
2. 檢查 TP/FP 分析結果
3. 驗證圖表生成正常
4. 確認適合發布
