---
allowed-tools: Bash(./scripts/validation/*), Bash(python3:*), Read, Glob, Grep
description: 執行自動化驗證框架 — 建置、跑樣本、比較 baseline、產出報告
---

# 驗證框架 (validate.sh)

根據使用者意圖，選擇合適的驗證模式：

## 使用者要求 Quick 驗證（或未指定模式）

```bash
./scripts/validation/validate.sh \
  --hypothesis "$ARGUMENTS" \
  --mode quick
```

## 使用者要求 Full 驗證

```bash
./scripts/validation/validate.sh \
  --hypothesis "$ARGUMENTS" \
  --mode full
```

## 使用者要求 TO (Tumor-Only) 模式

```bash
./scripts/validation/validate.sh \
  --hypothesis "$ARGUMENTS" \
  --track to \
  --mode quick
```

## 使用者要求平行多實驗

```bash
python3 scripts/validation/agents/orchestrator.py \
  --experiments scripts/validation/agents/experiments.yaml
```

## 執行後

1. 讀取實驗報告 (`experiment_report.md`)
2. 讀取比較結果 (`comparison.json`)
3. 向使用者報告 VERDICT 和關鍵指標
4. 如果 verdict = positive_pilot，建議跑 full 驗證確認
5. 如果 verdict = negative，分析回歸原因
