---
allowed-tools: Bash(./scripts/run_vcf_all_snv.sh:*)
description: 執行快速單點測試 (chr19 驗證，< 30 秒)
---

執行快速單點測試，驗證主流程未被破壞。

```bash
./scripts/run_vcf_all_snv.sh --mode chr19-verification
```

測試完成後：
1. 確認執行無錯誤
2. 報告關鍵指標（處理區域數、成功率）
3. 如有異常立即報告
