# 修改層級與規模梯度

## 修改層級 (Modification Level, ML)

> **注意**：本檔案的「ML1/ML2/ML3」指**程式碼修改層級**。
> 週報 `.claude/skills/weekly-report/references/LAYER_STRUCTURE.md` 的「Tier 1/2/3」指**證據卡報告詳細度**，兩者無關。
> `.claude/skills/validation-protocol/SKILL.md` 的「L1/L2/L3/L4」指**信號信任層級**，也是不同系統。

### ML1：Python 篩選腳本（預設）

修改目標：
- `scripts/evaluate_rescue_with_methylation.py`
- `scripts/analyze_gq_methylation_rescue_matrix.py`
- `scripts/run_batch_vcf_analysis.sh`（規則參數部分）

修改前必須：
1. Read 目標腳本相關片段
2. 確認不會破壞其他功能
3. 沙盒測試語法
4. Git commit baseline

### ML2：Python 特徵提取腳本

目標：`scripts/research_common.py` 或新建分析腳本

**升級條件**：ML1 確認有效的特徵，需提取新欄位。

### ML3：C++ 核心（最後手段）

目標：`src/core/*.cpp`, `include/core/*.hpp`

**升級條件（全部滿足）**：
- (a) ML1/ML2 在 S2 Medium 驗證 delta_f1 > 0.002
- (b) 需全基因組或多樣本高效能執行
- (c) 人類明確確認「請修改 C++」

未全部滿足時，呈現建議等待確認：
```
[ML3 建議] 建議修改 src/core/[file].cpp:
  位置: [函數名:行號區間]
  改動: [說明]
⚠ 需要您確認後才會執行（輸入「確認修改 C++」）
```

---

## Scale Ladder 規模梯度

| Scale | 條件 | 資料集 | 預期時間 |
|-------|------|--------|---------|
| S1 pilot | 新假設預設 | HCC1395 only | ~2-5 min |
| S2 medium | pilot 通過（delta > 0） | HCC1395 + 1-2 樣本 | ~10 min |
| S3 full | medium 一致 | 全部 6 樣本 | ~30-60 min |

---

## 平行執行（S2/S3 自動啟用）

**可平行**：
- 多 dataset × 同一分析腳本（讀取獨立、輸出獨立）
- 多個獨立 Python 分析腳本
- C++ build + Python 分析
- 同 dataset × 不同參數（輸出目錄不同）

**不可平行**：
- 後續腳本依賴前一腳本輸出（串行）

**S2 平行模板**（Agent tool 平行調度）：
```
Agent(prompt="執行 HCC1395_5kHz_TO benchmark...",
      subagent_type="tester", run_in_background=True)
Agent(prompt="執行 HCC1395_DORADO_TO benchmark...",
      subagent_type="tester", run_in_background=True)
```

**S3 平行模板**（parallel-benchmark agent）：
```
Agent(prompt="平行執行 S3 full benchmark:
  TASKS: [6 datasets]
  OUTPUT_DIR: cycles/${CYCLE_ID}/",
  subagent_type="parallel-benchmark")
```

---

## Pipeline Track 基線表

| Track | 資料集 | 基線 F1 | FP | FN | 命令 |
|-------|--------|---------|----|----|------|
| paired_full | HCC1395_5kHz | 0.8532 | 少 | 少 | `--mode HCC1395_5kHz_paired` |
| paired_pileup | HCC1395_DORADO | 0.8590 | 少 | 少 | `--mode HCC1395_DORADO_paired` |
| TO | HCC1395_5kHz | 0.7127 | 多 | 多 | `--mode HCC1395_5kHz_TO` |
| TO | HCC1395_DORADO | 0.7226 | 多 | 多 | `--mode HCC1395_DORADO_TO` |

**TO 優先規則**：paired track 連續 3 輪 `|delta| < 0.001` → 建議切換 TO。

---

## 擴展指南

### 新增 Pipeline Track

1. 在「Pipeline Track 基線表」新增行：Track 名稱、資料集、基線 F1、FP/FN 特性、`--mode` 命令
2. 基線 F1 必須從未修改的 baseline 程式跑出（不可估算）
3. 同步更新 `research-loop/SKILL.md` 的 frontmatter description（如增加新 track 名稱）

### 新增 Dataset（擴展 S2/S3 範圍）

1. 在「Scale Ladder 規模梯度」表中更新 S2/S3 的資料集欄位
2. 在 `.claude/skills/validation-protocol/SKILL.md` L3 固定順序清單中新增樣本名稱
3. 更新 L3 通過閾值的分母（如 7→8）
4. 確認 `scripts/run_batch_vcf_analysis.sh` 支援新 dataset 的 `--mode`
