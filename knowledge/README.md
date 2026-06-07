# InterSubMod Knowledge Base

**目標讀者**：AI Agent、協作研究者、未來的自己
**設計目的**：單一入口查詢 ISM 的 pipeline / 參數 / 樣本 / truth set / F1 / 結論，避免在 208 份 `docs/*.md` 中迷失
**建立日期**：2026-04-22
**對齊規範**：參考 `/big8_disk/liaoyoyo2001/knowledge/` 格式，但本庫不接 MCP（由 `scripts/` 靜態驗證）

---

## 🎯 5 秒查詢決策樹

問自己：**「我要查的是什麼？」**

```
              ┌──────────────────────────────────────────┐
              │        我要查什麼？                       │
              └──────────────────────────────────────────┘
                   │
        ┌──────────┼──────────┬──────────┬──────────┬──────────┐
        ▼          ▼          ▼          ▼          ▼          ▼
     參數/        pipeline   樣本       結論/       現在在     研究
     欄位格式     怎麼跑     資訊       歷史       做什麼      方法
        │          │          │          │          │          │
        ▼          ▼          ▼          ▼          ▼          ▼
   04_parameters 03_pipelines 02_samples 09_conclusions 10_research 07_derived
   05_data_formats           08_truth_                  _status    _features
                            and_benchmark
```

**記住**：每個頂層目錄有 `00_index.md`，先讀它再進具體文件。

---

## 📚 快速導航

### 最常見需求速查表

| 需求 | 去哪裡 |
|------|--------|
| **ISM 有哪三條 pipeline？差在哪？** | [03_pipelines/00_index.md](03_pipelines/00_index.md) |
| **某個 CLI 參數的預設值與用途** | [04_parameters/01_cli_arguments.md](04_parameters/01_cli_arguments.md) |
| **距離度量（NHD/L1/L2/BERNOULLI）的數學定義** | [04_parameters/02_distance_metrics.md](04_parameters/02_distance_metrics.md) |
| **significance_summary.csv 的 59 欄意義** | [05_data_formats/01_significance_summary_schema.md](05_data_formats/01_significance_summary_schema.md) |
| **7 個樣本的 truth set VCF 在哪** | [08_truth_and_benchmark/01_truth_set_registry.md](08_truth_and_benchmark/01_truth_set_registry.md) |
| **F1 怎麼算** | [08_truth_and_benchmark/02_f1_calculation.md](08_truth_and_benchmark/02_f1_calculation.md) |
| **怎麼跑一次完整 benchmark** | [06_workflows/02_full_vcf_analysis.md](06_workflows/02_full_vcf_analysis.md) |
| **編譯與快速測試** | [06_workflows/01_build_and_test.md](06_workflows/01_build_and_test.md) |
| **C++ 修改 PDD 6 步驟（Hard rule）** | [06_workflows/07_cpp_change_pdd.md](06_workflows/07_cpp_change_pdd.md) |
| **ΔF1 locked 數字（SoT，含完整 provenance）** | [03_pipelines/05_f1_baseline_canonical.md](03_pipelines/05_f1_baseline_canonical.md) |
| **新資訊如何補充與驗證** | [00_governance/07_new_info_protocol.md](00_governance/07_new_info_protocol.md) |
| **AI 執行模式切換（互動 / 全自動） + 暫停級別** | [00_governance/09_confirmation_protocol.md](00_governance/09_confirmation_protocol.md) |
| **`git commit` 被擋 / Hooks 配置** | [00_governance/08_hooks_and_automation.md](00_governance/08_hooks_and_automation.md) |
| **Opus 4.7 特性 + subagent 觸發語** | [AGENT.md](AGENT.md) §「Opus 4.7 Prompt 策略」 |
| **實作前準則（假設陳述 + 暫停矩陣）** | [00_governance/10_think_before_code.md](00_governance/10_think_before_code.md) |
| **可用分析腳本索引（155 個）** | [06_workflows/08_analysis_scripts_index.md](06_workflows/08_analysis_scripts_index.md) |
| **週報 / PPTX 怎麼做** | [06_workflows/09_pptx_and_weekly_report.md](06_workflows/09_pptx_and_weekly_report.md) |
| **HPFineNGroups 是什麼？結論如何？** | [07_derived_features/01_hpfinengroups.md](07_derived_features/01_hpfinengroups.md) |
| **哪些方向已被證實為 NEGATIVE（別重做）** | [09_conclusions/03_concluded_negative.md](09_conclusions/03_concluded_negative.md) |
| **外部文獻地景（7 角度：我們做法/困難/觀察的內外證據對照 + 跨問題依賴鏈）** | [11_external_literature/00_index.md](11_external_literature/00_index.md) |
| **目前最高優先研究方向** | [10_research_status/01_current_focus_snapshot.md](10_research_status/01_current_focus_snapshot.md) |
| **我想寫新的 KB 文件** | [00_governance/06_update_workflow.md](00_governance/06_update_workflow.md) |

---

## 🗂️ 目錄總覽

| 目錄 | 內容 | 文件數 | 典型 doc_type |
|------|------|--------|---------------|
| [00_governance/](00_governance/) | 本 KB 本身的規範（frontmatter schema、命名、查詢 SOP、更新流程） | 7 | reference |
| [01_project_overview/](01_project_overview/) | 專案定位、五大研究目標、當前 phase、突破策略 | 5 | explanation |
| [02_samples/](02_samples/) | 7 個主樣本各一份（HCC1395、COLO829、H1437、H2009、HCC1937、HCC1954、HCC1395_DORADO） | 8 | reference |
| [03_pipelines/](03_pipelines/) | ★ 三條 pipeline（paired_full / paired_pileup / TO）完整規格 | 5 | reference + howto |
| [04_parameters/](04_parameters/) | C++ CLI 參數、距離度量、統計方法、過濾規則 | 6 | reference |
| [05_data_formats/](05_data_formats/) | 輸出檔案欄位字典（significance_summary、master_dataset、per-region） | 6 | reference |
| [06_workflows/](06_workflows/) | build/test、full VCF 分析、batch TP/FP、F1 benchmark、debug | 7 | howto |
| [07_derived_features/](07_derived_features/) | HPFineNGroups、Coverage_Multiple、Zone-Aware 等衍生特徵 | 6 | explanation |
| [08_truth_and_benchmark/](08_truth_and_benchmark/) | ★ 7 樣本 truth set + F1 計算 + benchmark 協議 | 4 | reference |
| [09_conclusions/](09_conclusions/) | 研究結論索引（positive / characterization / NEGATIVE） | 6 | reference-summary |
| [10_research_status/](10_research_status/) | 時間敏感：CURRENT_FOCUS 快照、active 假說、阻塞、里程碑 | 6 | reference |
| [scripts/](scripts/) | KB 維運腳本（frontmatter 驗證、related_ids 對稱性、等） | 4 | — |

**總計**：約 **70 份文件**（含 scripts），對齊既有 `/big8_disk/knowledge/` 規模。

---

## 🔑 關鍵術語對照

| 術語 | 英文 | 在本專案的意義 |
|------|------|---------------|
| ISM | Inter-Subclonal Methylation Analysis | 本專案名稱 |
| paired_full | — | ClairS paired full VCF + LongPhase-S haplotag pipeline（canonical benchmark） |
| paired_pileup | — | ClairS pileup VCF + LongPhase-S pipeline |
| TO / tumor_only | — | ClairS-TO + LongPhase-TO pipeline（有 self-phasing bias） |
| HP tag | Haplotype tag | `HP:i:1/2` germline、`HP:i:11/21/33` somatic（LongPhase-S 衍生） |
| HPFineNGroups | — | 四分群（HP1-G/HP1-S/HP2-G/HP2-S）中有效群組數；positive characterization marker |
| Coverage_Multiple | — | NumReads / 期望覆蓋度；CN proxy（r≈0.83） |
| Zone-Aware | — | 基於 AF/LOH/HP 的 5 分區框架（Z1-Z5）；characterization only |
| Self-Phasing | — | TO 模式無 normal 參照導致 HP tag 偏差（94.6% → HP1） |

---

## 🛠️ KB 使用守則（必讀）

1. **先讀 00_index.md 再進具體文件** — 每個目錄的 index 會告訴你哪份是你要的
2. **信任 canonical_paths 不信任搜尋** — 每份文件 frontmatter 的 `canonical_paths` 是官方路徑；`alias_paths` 只是歷史別名
3. **看 last_verified** — 如果 `last_verified` 超過 90 天，執行前請用文件內「可直接執行命令」驗證現況
4. **研究結論類文件只做索引** — 真正詳盡報告在 `docs/reports/research_landscape/`，KB 只給狀態與跳轉連結
5. **別改變既有 docs/** — KB 建立後 `docs/` 仍是最完整權威來源，KB 是精煉與導航層
6. **遇模糊主題** — 到 [00_governance/05_query_protocol.md](00_governance/05_query_protocol.md) 看 5 種典型查詢 SOP

---

## 🔄 更新頻率

| 目錄 | 更新頻率 | 觸發條件 |
|------|----------|----------|
| 00_governance/ | 穩定 | KB 架構變更 |
| 01-05（規格類） | 低頻 | 參數/格式異動、新 pipeline 加入 |
| 06_workflows/ | 中頻 | 腳本異動 |
| 07_derived_features/ | 中頻 | 新 feature 或 feature 結論變化 |
| 08_truth_and_benchmark/ | 低頻 | 新樣本、truth set 更新 |
| 09_conclusions/ | 中頻 | 新結論產生時追加索引 |
| 10_research_status/ | **高頻（每週）** | 每週研究迴圈、CURRENT_FOCUS 變動 |

---

## 🚀 進入點（按身份）

- **🤖 AI Agent**：先讀 [AGENT.md](AGENT.md) 再進 [00_governance/05_query_protocol.md](00_governance/05_query_protocol.md)
- **👤 新進研究者**：按順序讀 [01_project_overview/](01_project_overview/) → [03_pipelines/00_index.md](03_pipelines/00_index.md) → [02_samples/00_index.md](02_samples/00_index.md)
- **🔧 工程師（想跑 pipeline）**：直接讀 [06_workflows/01_build_and_test.md](06_workflows/01_build_and_test.md) → [06_workflows/02_full_vcf_analysis.md](06_workflows/02_full_vcf_analysis.md)
- **🔬 研究者（想知道結論）**：讀 [09_conclusions/00_index.md](09_conclusions/00_index.md)

---

## 📎 延伸資源（既有專案文件，KB 不重寫）

- **研究景觀總索引**：`docs/reports/research_landscape/00_INDEX.md`（11 份分題報告）
- **當前進度**：`docs/CURRENT_FOCUS.md`
- **實驗歷史**：`docs/experiments/INDEX.md`
- **專案規則**：`.claude/CLAUDE.md`
- **AI memory**：`/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/MEMORY.md`

---

**KB 版本**：v0.1（Phase 1 骨架）
**上次整體驗證**：2026-04-22
**維護者**：InterSubMod Research Team
