---
id: ism-kb-00-governance-query-protocol
name: "查詢流程 SOP（5 種情境）"
description: "★ AI agent 或協作者查詢 KB 時的標準流程：情境 A-E 對應典型問題，每情境有步驟、範例、可執行命令。避免猜測文件位置或跳過 index。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: howto
verified_scope: "query protocol SOP with 5 scenarios"
related_ids:
  - ism-kb-00-governance-cross-reference-rules
  - ism-kb-00-governance-frontmatter-schema
  - ism-kb-00-governance-new-info-protocol
  - ism-kb-06-workflows-cpp-change-pdd
  - ism-kb-00-governance-confirmation-protocol
  - ism-kb-00-governance-think-before-code
tags: [governance, query, sop, workflow, agent, protocol]
canonical_paths: [00_governance/05_query_protocol.md]
alias_paths: []
---

# 查詢流程 SOP（5 種情境）

- 一句結論：收到查詢時先分類為情境 A-E，再按對應 SOP 導航；**不要猜檔名、不要跳過 00_index**
- 適用對象：AI agent（包括 Claude）、新進研究者、自己
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 示範：查「--min-mapq 預設是多少？」
  cat /big7_disk/liaoyoyo2001/InterSubMod/knowledge/04_parameters/00_index.md | head -30
  ```

---

## 📐 通用 SOP（所有情境共用）

1. **分類問題** → 對照情境 A-E
2. **進對應目錄的 `00_index.md`** → 讓 index 告訴你該讀哪份
3. **讀具體文件的 frontmatter** → 確認 `last_verified` 與 `source_type`
4. **讀 body 的「可直接執行命令」** → 必要時驗證現況
5. **引用回使用者** → 給 KB 路徑 + 必要時跳 docs/ 或 file:line

---

## 🅰️ 情境 A — 「某參數的意義 / 預設值 / 取值範圍？」

### 觸發關鍵字
"--min-mapq、"距離度量、"methyl_high"、"filter 閾值"、"OpenMP 幾 threads"

### SOP
```
Step 1: 讀 knowledge/04_parameters/00_index.md 判斷分組
  ├── CLI 參數 → 01_cli_arguments.md
  ├── 距離度量 → 02_distance_metrics.md
  ├── 統計方法 → 03_statistical_methods.md
  ├── 過濾規則 → 04_filtering_rules.md
  └── 內部 Config 常數 → 05_config_defaults.md

Step 2: 讀對應檔案的參數表，找目標參數
Step 3: 若需實作細節 → 跳原始碼 file:line（frontmatter 或表格內會列）
Step 4: 用「可直接執行命令」驗證：
  ./build/bin/inter_sub_mod --help 2>&1 | grep -A 1 "<param-name>"
```

### 範例對話
```
User: --min-mapq 預設是什麼？
Agent:
  1. 情境 A → 04_parameters/00_index.md
  2. Index 指向 01_cli_arguments.md 表格「1.3 讀取過濾參數」
  3. 答案：預設 20，取值 [0, 60]，出處 include/utils/ArgParser.hpp:67-68
  4. 用途：最小映射品質，低於此值的 read 會被過濾
  5. （可選）跑 `./build/bin/inter_sub_mod --help` 驗證
```

---

## 🅱️ 情境 B — 「某 pipeline 怎麼跑 / 輸入輸出是什麼？」

### 觸發關鍵字
"paired_full 怎麼跑"、"TO pipeline"、"執行 benchmark"、"VCF 來源"

### SOP
```
Step 1: 讀 knowledge/03_pipelines/00_index.md 選 pipeline
  ├── paired_full（canonical benchmark） → 01_paired_full.md
  ├── paired_pileup（pileup 模型）       → 02_paired_pileup.md
  └── tumor_only（有 self-phasing bias）  → 03_tumor_only.md
  （對比表 → 04_pipeline_comparison.md）

Step 2: 讀對應 pipeline 文件的「輸入鏈路 → 參數 → 命令 → 輸出結構」
Step 3: 若要實際跑 → 跳 06_workflows/02_full_vcf_analysis.md 看批次 script
Step 4: 若要 F1 benchmark → 跳 08_truth_and_benchmark/02_f1_calculation.md
```

### 範例對話
```
User: paired_full pipeline 的典型命令是什麼？
Agent:
  1. 情境 B → 03_pipelines/00_index.md
  2. 選 01_paired_full.md
  3. 文件 body 有 canonical command:
     ./scripts/run_vcf_all_snv.sh --caller-mode paired --mode all-with-w5000 \
       --metrics BERNOULLI --threads 120 -o output/paired_full_run
  4. 輸入：ClairS paired full VCF + tumor+normal BAM + LongPhase phased.bed
  5. 輸出：output/canonical/<sample>/paired_full/<timestamp>/
  6. F1 計算：跳 08_truth_and_benchmark/02_f1_calculation.md
```

---

## 🅲 情境 C — 「某主題的研究結論？」

### 觸發關鍵字
"HPFineNGroups 能不能當 filter"、"LOH 與甲基化的關係"、"Zone-Aware 結論"、"xxx 還在嗎？"

### SOP
```
Step 1: 讀 knowledge/09_conclusions/00_index.md 看主題狀態
  ├── positive（已驗證為真實信號） → 01_positive_findings.md
  ├── characterization（能描述但不能做 filter） → 02_characterization_only.md
  ├── NEGATIVE / NO-GO → 03_concluded_negative.md
  └── 仍進行中 → 10_research_status/

Step 2a [positive]: 讀 07_derived_features/ 的具體 feature 文件
  → 跳 docs/reports/research_landscape/XX.md 看完整證據

Step 2b [characterization]: 同上；但在回答時明確警告「不可做 variant filter」

Step 2c [NEGATIVE]: 讀 03_concluded_negative.md 看失敗原因
  → **絕對不要建議重做**，除非使用者明確表示要挑戰既有結論

Step 2d [ongoing]: 讀 10_research_status/02_active_hypotheses.md
  → 跳 research/autoresearch/hypothesis_queue.json 看詳細假說

Step 3: 引用時同時給 KB 路徑 + docs/reports/research_landscape/ 路徑
```

### 範例對話
```
User: HPFineNGroups 能當 variant filter 嗎？
Agent:
  1. 情境 C → 09_conclusions/00_index.md
  2. 狀態：characterization only（能描述，不能 filter）
  3. 細節：07_derived_features/01_hpfinengroups.md
     → N≥4+NR≥80 canonical filter，TP rate 92.81%
     → 但若用作 variant filter 會產生負 F1 增益
  4. 權威報告：docs/reports/research_landscape/09_Part_B.md
  5. ⚠️ 警告：2026-04-21 發現 flag=on 下 N≥3 可能是 somatic HP tag artifact
```

---

## 🅳 情境 D — 「樣本資訊 / truth set / F1 計算？」

### 觸發關鍵字
"HCC1395 truth set 在哪"、"COLO829 幾個 truth variants"、"F1 怎麼算"、"哪個樣本最適合 benchmark"

### SOP
```
Step 1: 分辨要查樣本還是 truth set
  ├── 單一樣本資訊 → 02_samples/<sample>.md
  ├── 7 樣本 truth set 全表 → 08_truth_and_benchmark/01_truth_set_registry.md
  ├── F1 計算公式 → 08_truth_and_benchmark/02_f1_calculation.md
  └── benchmark 工具（som.py / hap.py）→ 08_truth_and_benchmark/03_benchmark_protocols.md

Step 2: 進具體文件
Step 3: 跳權威文件 docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md
Step 4: 驗證路徑存在：
  ls /big8_disk/data/<SAMPLE>/
```

### 範例對話
```
User: HCC1395 的 truth set 在哪？多少 variants？
Agent:
  1. 情境 D → 08_truth_and_benchmark/01_truth_set_registry.md
  2. HCC1395 使用 SEQC2 high-confidence set：
     /big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
     TRUTH_TOTAL: 39,447
     有 HC BED 過濾區
  3. 權威協議：docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md
  4. 驗證：bcftools view <path> | grep -v "^#" | wc -l
```

---

## 🅴 情境 E — 「現在在做什麼 / 下一步？」

### 觸發關鍵字
"current focus"、"最高優先"、"下一步"、"阻塞在哪"、"Phase 2 進度"

### SOP
```
Step 1: 讀 knowledge/10_research_status/00_index.md
  ├── 當前焦點快照 → 01_current_focus_snapshot.md
  ├── active 假說 → 02_active_hypotheses.md
  ├── evidence ledger 規範 → 03_evidence_ledger_format.md
  ├── 阻塞/風險 → 04_blockers_and_risks.md
  └── 下一里程碑 → 05_next_milestones.md

Step 2: 讀快照文件；注意開頭「⚠ 有效期 2 週」警示
Step 3: 若快照 > 2 週或有疑慮 → 跳 docs/CURRENT_FOCUS.md 核對
Step 4: 若需最新 evidence → 跳 research/autoresearch/evidence_ledger.jsonl
```

### 範例對話
```
User: 現在在做什麼？下一步是什麼？
Agent:
  1. 情境 E → 10_research_status/01_current_focus_snapshot.md
  2. 快照：Phase 2 方向 A+D（Normal Methylation Reference + CN/Purity-aware）
  3. 阻塞：haplotag + ISM 全量重跑待執行（需 PON-only phasing）
  4. 下一里程碑：7 樣本 × paired_full 重跑完成後正式 Phase 2A 分析
  5. ⚠️ 快照日期：2026-04-22；若已超過 2 週請讀 docs/CURRENT_FOCUS.md 核對
```

---

## 📌 遇模糊查詢的處理原則

使用者問題橫跨多情境或不清楚時：

1. **Clarifying question 優先**：問使用者想達成什麼目標
2. **默認順序**：A（參數） → B（執行） → C（結論） → D（樣本/truth） → E（現況）
3. **廣覆蓋查詢**：從 [README.md](../README.md) 的「快速導航」表開始，或讀相關多份 index

### 範例：模糊查詢
```
User: LOH 相關的東西
Agent:
  不清楚你想知道：
  - (A) LOH 相關的 CLI 參數 (--loh-bed)
  - (B) 如何產生 LOH BED 檔案（workflow）
  - (C) LOH 對甲基化訊號的研究結論
  - (D) 各樣本的 LOH 狀況
  請問具體是？
```

---

## ✅ 反模式（禁止）

- ❌ 收到查詢直接 `grep -r` 搜整個 repo（浪費；應先進 KB）
- ❌ 跳過 `00_index.md` 直接猜檔名（可能遺漏相關資訊）
- ❌ 引用 frontmatter 未列的「相關」文件而不核對
- ❌ 看到 `last_verified > 90d` 仍直接引用不驗證
- ❌ 複製 docs/ 內容到 KB 回答（應跳轉到 docs/ 權威）

---

**相關**：
- AGENT 指引：[../AGENT.md](../AGENT.md)
- 交叉引用：[03_cross_reference_rules.md](03_cross_reference_rules.md)
- Freshness：[04_freshness_policy.md](04_freshness_policy.md)
