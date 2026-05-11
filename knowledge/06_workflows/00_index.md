---
id: ism-kb-06-workflows-index
name: "Workflows 目錄索引"
description: "ISM 執行流程 howto 索引：build/test、full VCF 分析、batch TP/FP、F1 benchmark、debug、新增樣本。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "workflows directory"
related_ids:
  - ism-kb-06-workflows-build-and-test
  - ism-kb-06-workflows-full-vcf-analysis
  - ism-kb-06-workflows-batch-tp-fp-analysis
  - ism-kb-06-workflows-f1-benchmark
  - ism-kb-06-workflows-debug-and-logging
  - ism-kb-06-workflows-adding-new-sample
  - ism-kb-06-workflows-cpp-change-pdd
  - ism-kb-06-workflows-analysis-scripts-index
  - ism-kb-06-workflows-pptx-and-weekly-report
tags: [workflow, index, howto]
canonical_paths: [06_workflows/00_index.md]
alias_paths: []
---

# Workflows 目錄索引

- 一句結論：6 個工作流程文件，從 build/test 到 F1 benchmark；每份含可直接執行命令
- 適用對象：第一次用 ISM、執行 benchmark、debug
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/InterSubMod/scripts/
  ```

---

## 文件列表

| 檔案 | 主題 | 時長 |
|------|------|------|
| [01_build_and_test.md](01_build_and_test.md) | cmake/make/ctest + 3 種 test 指令 | 2-5 min |
| [02_full_vcf_analysis.md](02_full_vcf_analysis.md) | `run_vcf_all_snv.sh` 三個 mode | 30s-5min |
| [03_batch_tp_fp_analysis.md](03_batch_tp_fp_analysis.md) | `run_batch_vcf_analysis.sh` TP+FP 並行 | 5-10 min |
| [04_f1_benchmark.md](04_f1_benchmark.md) | som.py / hap.py / pipeline/run_benchmark.sh | 10-20 min |
| [05_debug_and_logging.md](05_debug_and_logging.md) | log-level debug + output-filtered-reads | — |
| [06_adding_new_sample.md](06_adding_new_sample.md) | 新樣本入庫（truth set 登記、VCF symlink） | — |
| [07_cpp_change_pdd.md](07_cpp_change_pdd.md) | ★ **C++ 修改 PDD 6 步驟**（Hard rule） | 前置 methodology-audit |
| [08_analysis_scripts_index.md](08_analysis_scripts_index.md) | 155 個分析腳本按類別索引 + 15 高頻詳細表 | 找可用分析工具、寫新腳本前 |
| [09_pptx_and_weekly_report.md](09_pptx_and_weekly_report.md) | PPTX 與週報 3 skill 對照 + manual 索引 | 每週週報、實驗報告前 |

---

## 常用 /slash-commands

這些是專案定義的 slash commands（`.claude/commands/`）：

| Command | 對應 workflow | 時長 |
|---------|--------------|------|
| `/build` | [01_build_and_test.md](01_build_and_test.md) | ~1 min |
| `/test-quick` | chr19 verification | <30 sec |
| `/test-data` | light data check | ~1 min |
| `/test-full` | full batch TP+FP | ~5 min |
| `/validate` | 自動化驗證框架 | ~5 min |

---

## 執行 ISM 的最短路徑

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod

# 1. Build
cd build && cmake .. && make -j$(nproc) && cd ..

# 2. Quick verification
./scripts/run_vcf_all_snv.sh --mode chr19-verification

# 3. Full benchmark
./scripts/run_batch_vcf_analysis.sh
```

---

## 相關

- Pipelines：[../03_pipelines/](../03_pipelines/)
- 參數：[../04_parameters/](../04_parameters/)
