---
allowed-tools: Bash(./scripts/run_vcf_all_snv.sh:*)
description: 執行完整數據測試 (約 1 分鐘)
---

執行完整數據測試，驗證主要數據正確性。

```bash
./scripts/run_vcf_all_snv.sh --mode all-with-w5000 --plot-type no
```

測試完成後：
1. 比對與上次測試的差異
2. 報告處理區域數、顯著區域數
3. 檢查是否有新的錯誤或警告
