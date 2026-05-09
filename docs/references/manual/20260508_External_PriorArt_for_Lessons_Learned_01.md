---
title: 開發流程複盤 — 外部 Prior Art 掃描報告
date: 2026-05-08
status: validated
type: external_research
classification: prior_art_review
related: docs/reports/validated/2026/05/20260508_Development_Process_Lessons_Learned_01.md
audience: liaotzuyu000 / InterSubMod harness designer
---

> **Disclaimer**: 以下方案均來自網路搜尋（2026-05-08 wall-time ~25 min），**非親身驗證使用**。引用閾值與設計片段以官方文件描述為準，採納前需在 InterSubMod 上下文 evaluate。所有 URL 經搜尋結果引用，未深入訪問每篇全文；若採納應再讀官方原文確認版本與當前 API。

---

## §1 Bottom Line

InterSubMod 在 2026-04 ~ 2026-05 暴露的 5 大類問題（前置條件 / 結論健壯性 / 認知流程 / 基礎設施 / 報告呈現），**前 4 類在外部生態都有成熟 prior art 可直接 cherry-pick**，第 5 類 prior art 偏分散但有局部可用方案。最具直接借鏡價值的三個方案：

1. **dbt source freshness + model contracts**（Category A）— `warn_after`/`error_after` 二段門檻、`source_status:fresher+` 級聯重跑語法可直接對應 InterSubMod 的 `/check-staleness` skill 與 `post_cpp_commit_invalidate.sh`，把 KB 14d freshness 與 dataset schema 一起納入單一框架
2. **Anthropic evaluator-optimizer + Devin plan checkpoint**（Category C）— Devin 的「plan as the most valuable checkpoint」對應本專案 anchor #6 + L4 mandatory；`anthropic-cookbook/patterns/agents/evaluator_optimizer.ipynb` 是直接可抄的 reference 實作
3. **Prometheus node_exporter textfile collector + systemd timer disk monitor**（Category D）— 解 800 GB 災情的最低成本方案；textfile collector + 10-min 週期 + 90% 警 / 95% 阻擋，可在無 Prometheus server 時退化為純 logfile + cron，仍保留可觀測性

整體投入產出比評估：**前 4 類採納 = high ROI**（可直接抄設計，1-2 天落地）；第 5 類 = medium ROI（需自寫 lint，外部僅給靈感）。

---

## §2 Category A — 前置條件問題

對應 InterSubMod 痛點：dataset schema 變動下游未 invalidate / merged dataset AF 欄位漂移 / symlink 指向錯誤 caller / KB 14d freshness 警告無強制重認證。已有 `/check-staleness` skill + `post_cpp_commit_invalidate.sh` hook。

### A-1 dbt source freshness + model contracts
- URL: <https://docs.getdbt.com/docs/deploy/source-freshness> / <https://docs.getdbt.com/reference/resource-properties/freshness> / <https://medium.com/@faouaghzene/dbt-model-contracts-beyond-tests-toward-robust-governance-6814c1eda642>
- 適用範圍：data warehouse SQL 變數模型，但其 freshness/contract 語意架構可直接借到 file-based research pipeline
- **可借鏡點**：
  - `freshness: { warn_after: {count: 12, period: hour}, error_after: {count: 24, period: hour} }` 二段門檻語法 — 可直接套用到 KB 14d / 30d 強制重認證
  - `loaded_at_field` + 自動 SQL 查 max(timestamp) 的 schema → 對應到 dataset 寫入 mtime 自動掃描
  - `dbt build --select source_status:fresher+` 級聯重跑語法 — 對應 InterSubMod「下游 stage 自動失效」的 invalidation 需求
  - **model contracts**：YAML 宣告 column types / constraints，break 即 fail — 可用來防 `caller_af` 欄位語意漂移
- **採納建議**：B 部分採。InterSubMod 不需要整套 dbt，但 freshness YAML schema + `fresher+` selector 概念值得抄到 `/check-staleness` skill
- **既有 cross-ref**：`/check-staleness` skill、`post_cpp_commit_invalidate.sh` hook、`feedback_merged_dataset_af_and_loh_pitfalls`

### A-2 DVC pipelines（content-hash invalidation）
- URL: <https://doc.dvc.org/user-guide/pipelines> / <https://doc.dvc.org/command-reference/repro> / <https://doc.dvc.org/start/data-pipelines/data-pipelines>
- 適用範圍：ML/research file-based pipeline，git-friendly
- **可借鏡點**：
  - DVC 對每個 stage 的 dependency 計算 **content hash**（非 mtime），content 沒變不重跑 — 對應 InterSubMod 重跑判斷
  - `dvc repro` 自動偵測 dirty stages 並級聯 — 對應 binary commit 後下游 stale
  - `dvc.yaml` 宣告式 stage / deps / outs / params — 比純 shell hook 更穩定
- **採納建議**：A 直接採（小範圍 pilot）。建議先在 `output/synthesis/` 一個 round 試裝 DVC stage，driver 跑通再評估全 repo 採用
- **既有 cross-ref**：`output/OUTPUT_INDEX.md` 信任度規範、`docs/data_specs/20260414_output資料信任度與生命週期_01.md`

### A-3 Great Expectations data validation
- URL: <https://greatexpectations.io/> / <https://docs.greatexpectations.io/docs/reference/learn/data_quality_use_cases/schema/>
- 適用範圍：Python pandas/Spark 資料 batch 驗證；schema drift / distribution drift
- **可借鏡點**：
  - **Expectation Suite** = YAML 宣告期望（column 存在 / type / 值域 / null %）
  - `Anomaly Detection` UI flag — 自動偵測「columns don't diverge over time」
  - 對接 CI/CD：validation fail → block pipeline
  - 對 InterSubMod 最直接的用法：在 merged dataset 載入處 attach Expectation，AF 欄位若 p75 < 0.06 直接 fail（對應 S5 教訓）
- **採納建議**：A 直接採。可寫 5-10 條 expectation 給 master_extended / merged datasets，當作 dataset-level pitfalls 自動化執行
- **既有 cross-ref**：pitfalls_table、`feedback_merged_dataset_af_and_loh_pitfalls`、`auc-confound-guard` skill

### A-4 lakeFS / Pachyderm（重型方案，僅參考）
- URL: <https://github.com/pachyderm/pachyderm> / <https://lakefs.io/blog/data-lineage-for-data-lakes/>
- 適用範圍：企業級 data lake 版本控制 + lineage
- **可借鏡點**：containerized stage + git-like commit on data；對 single-researcher 太重
- **採納建議**：C 不採（理由：基礎設施成本過高，single-researcher 無 ROI）。但「git-like commit on data」概念值得記在心裡；若未來開源/多人合作可再評估
- **既有 cross-ref**：N/A

### A-5 Data Contract CLI
- URL: <http://cli.datacontract.com/>
- 適用範圍：跨團隊 data contract 宣告（schema + SLA + quality rules + tests）
- **可借鏡點**：YAML 宣告 contract → CLI lint + test；single source of truth
- **採納建議**：B 部分採。InterSubMod 的 `caller_af` / `LOH_flag` / `master_extended` 可寫 mini contract YAML，當作 schema 文件 + 自動測試入口
- **既有 cross-ref**：`docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md`

---

## §3 Category B — 結論健壯性問題

對應 InterSubMod 痛點：cross-sample 不一致 / effect size noise / confound 未 sweep / pitfall 未編入規則。已有 6-component composite + per-component override + auc-confound-guard。

### B-1 nf-core / Snakemake 多樣本 QC pattern
- URL: <https://nf-co.re/rnaseq/3.14.0/> / <https://github.com/nf-core/rnaseq> / <https://www.biorxiv.org/content/10.1101/610741v3.full>
- 適用範圍：bioinformatics 多樣本 pipeline 標準化
- **可借鏡點**：
  - nf-core 強制 7 大 QC 檢查（raw read / alignment / gene biotype / sample similarity / strand-specificity / contamination / duplication）
  - 「peer-reviewed」community pattern — 任何新 pipeline 要進主線需通過 review
  - **portability/reproducibility/scalability/parallelism** 四維度 inherent — 是 pipeline 設計時的硬性約束，不是事後補
- **採納建議**：B 部分採。InterSubMod 不會切到 Nextflow，但 nf-core 的「7 維度 QC + peer review gate」可對應到 `multi-sample-consistency` skill 的 6-component composite
- **既有 cross-ref**：`multi-sample-consistency` skill、6-component composite、`auc-confound-guard`

### B-2 MultiQC 跨樣本 aggregation
- URL: <https://github.com/MultiQC/MultiQC> / <https://academic.oup.com/bioinformatics/article/32/19/3047/2196507>
- 適用範圍：bioinformatics 跨工具跨樣本 QC 報告聚合
- **可借鏡點**：
  - 自動掃 log 檔 → 單一 HTML 報告 → side-by-side 表
  - 「shared plots allow accurate comparison between samples, allowing detection of subtle differences」— 直接對應 InterSubMod 的「肉眼檢視 2x2 × 7 samples 固定順序」需求
  - Built-in 150+ tools — 對應 plug-in based 設計
- **採納建議**：A 直接採。InterSubMod 的 weekly report 可仿 MultiQC pattern，建立「7 樣本 × N feature × M run」的單一 HTML aggregate，作為 cross-sample consistency 視覺驗證入口
- **既有 cross-ref**：`feedback_figure_layout_standard`、`feedback_visual_inspection_requirement`

### B-3 Anthropic evaluator-optimizer pattern
- URL: <https://www.anthropic.com/research/building-effective-agents> / <https://github.com/anthropics/anthropic-cookbook/blob/main/patterns/agents/evaluator_optimizer.ipynb> / <https://www.anthropic.com/engineering/demystifying-evals-for-ai-agents>
- 適用範圍：LLM agent loop（generator → evaluator → feedback → regenerate）
- **可借鏡點**：
  - 「two signs of good fit: LLM responses can be demonstrably improved when a human articulates feedback, AND the LLM can provide such feedback」— InterSubMod 的 hypothesis 驗證與 conclude-research 流程符合
  - cookbook ipynb 是直接可抄的 reference 實作
  - **避免使用情境**：first-attempt OK / 主觀標準 / time-cost 限制 — 對應 InterSubMod 不要對所有 cycle 都加 evaluator，只在「結論可能撤回」高風險時段啟用
- **採納建議**：A 直接採。InterSubMod 的 `run-evaluator` skill 可仿照 cookbook 結構升級
- **既有 cross-ref**：`run-evaluator` skill、`review-evidence` skill、L4 mandatory

### B-4 GATK Best Practices（VQSR + hard filter）
- URL: <https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices> / <https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows>
- 適用範圍：variant calling 標準
- **可借鏡點**：
  - 兩階段 filter：VQSR（machine-learning, 需大樣本）+ hard filter（rule-based, 小樣本）— 對應 InterSubMod 對「composite score」與「per-component override」的二段判斷
  - 階段化 pipeline（pre-process → variant discovery → callset evaluation → filter）— 每階段有獨立驗收標準
- **採納建議**：B 部分採。GATK 的「callset evaluation」階段可借鏡到 InterSubMod 的 conclude-research 前 evidence-tier-audit
- **既有 cross-ref**：`provenance-tier-audit` skill、6-component composite

### B-5 OSF Preregistration
- URL: <https://www.cos.io/initiatives/prereg> / <https://help.osf.io/article/330-welcome-to-registrations>
- 適用範圍：學術 hypothesis 公開預註冊
- **可借鏡點**：
  - 「conditional analysis plan」— 預先寫「if normality 不符，改用 non-parametric」— 直接對應 InterSubMod 的 hypothesis branching 與 sensitivity 6/6 / specificity 2/2 預設值
  - **time-stamped read-only** — 防止 post-hoc 修改假說
  - 「separates hypothesis-generating from hypothesis-testing」— 對應 InterSubMod 多次撤回時 exploratory vs confirmatory 模糊化的問題
- **採納建議**：A 直接採（內部版）。InterSubMod 不需公開 OSF，但可在 `cycle-init` skill 強制寫 conditional plan + sensitivity threshold + 預期 fail mode，git commit 鎖定後才開跑
- **既有 cross-ref**：`cycle-init` skill、`inject-hypothesis` skill、Drill 1 retrospective

---

## §4 Category C — 認知/流程/風格問題

對應 InterSubMod 痛點：pivot 過頻 / 過早結論 / 模糊指令 / 執行模式分級 / anchor 風格硬化 / weekly review ritual fatigue。已有暫停判斷矩陣 + execution mode hierarchy + L4 mandatory + 5 段報告骨架。

### C-1 Devin plan-as-checkpoint pattern
- URL: <https://sureprompts.com/blog/devin-ai-prompting-guide> / <https://www.builder.io/blog/devin-vs-cursor> / <https://battleaitools.com/ai-agents/cursor-agent-vs-devin-vs-claude-code/>
- 適用範圍：long-running autonomous coding agent
- **可借鏡點**：
  - 「Devin proposes a plan before executing — treat as the most valuable checkpoint. A bad plan runs for an hour before you notice; a reviewed plan catches the drift in a minute.」— 直接對應 InterSubMod 的「strategy 同意後逐項實作確認」memory 規則
  - 「Keep agent sessions under 2 hours or add re-indexing checkpoints」— 對應 anchor #6 + L4 mandatory
  - 「.cursorrules / CLAUDE.md / .mdc 是 mandatory — 自主 agent 需要嚴格規則」— 已實踐
- **採納建議**：A 直接採。anchor #6 PARKED 可解鎖實作；明確把「plan review = 最高價值 checkpoint」寫進 `cycle-init` skill 與 confirmation-protocol
- **既有 cross-ref**：`cycle-init` skill、`confirmation-protocol`、anchor #6、`feedback_strategy_then_per_item_confirmation`、`feedback_execution_mode_hierarchy`

### C-2 Anthropic 6 composable agent patterns
- URL: <https://resources.anthropic.com/building-effective-ai-agents> / <https://aimultiple.com/building-ai-agents>
- 適用範圍：agent 架構選型
- **可借鏡點**：
  - 6 patterns：prompt chaining / routing / parallelization / orchestrator-workers / evaluator-optimizer / autonomous agents
  - 「start simple, add complexity only when demonstrably improves performance」— 對 InterSubMod 7-phase × M1+M2+M4 governance harness 的設計理念吻合
  - **明確 anti-pattern**：複雜框架 / specialized libraries 反而拖慢成功率
- **採納建議**：B 部分採。已實踐大部分；可作為新 skill 設計時的 sanity check 對照表
- **既有 cross-ref**：`research-loop` skill、harness plan v1.7

### C-3 Cognitive forcing strategies + DECLARE
- URL: <https://bmcmededuc.biomedcentral.com/articles/10.1186/s12909-018-1444-3> / <https://pmc.ncbi.nlm.nih.gov/articles/PMC10149772/> / <https://pmc.ncbi.nlm.nih.gov/articles/PMC3786644/>
- 適用範圍：醫學診斷 cognitive bias 對抗
- **可借鏡點**：
  - 「SLOW」助記符 — 字面提醒 slow down + 每字母對應一個 specific bias counterforce
  - **DECLARE**（Differential / Evidence / Confounders / Likelihood / Alternatives / Reassess / Engage）— 7 字 cognitive checklist，對 InterSubMod 的 conclude-research 5 段骨架可作為對照升級
  - **三層強度**：universal / generic / specific — 對應 confirmation-protocol 的 Hard Gate / Gate / Review / FYI 四級
  - **負面證據**：CFS 在新手 randomized trial 不顯效（PMID 24423999）— 提醒「checklist 不可取代 evidence-driven validation」
- **採納建議**：B 部分採。DECLARE 7 字可抄入 `conclude-research` skill 與 `review-evidence` skill 作為 explicit checklist
- **既有 cross-ref**：`conclude-research`、`review-evidence`、`confirmation-protocol`、暫停判斷矩陣

### C-4 Quarto + Sumatra（reproducible research workflow）
- URL: <https://quarto.org/about.html> / <https://sumatra.readthedocs.io/en/latest/> / <https://github.com/open-research/sumatra> / <https://www.fharrell.com/talk/rflow/>
- 適用範圍：computational reproducibility 工具鏈
- **可借鏡點**：
  - **Quarto**：literate programming + 文檔即可執行；對應 InterSubMod 的「驗證 .md 內 numbers 對得上」
  - **Sumatra**：自動 capture 每次 run 的 git SHA / params / output / dependency — 對應 evidence_ledger 的自動化
  - `smt run` 命令 — 直接捕獲 context；可仿 InterSubMod 的 `provenance-tier-audit` skill
- **採納建議**：B 部分採。Quarto 可考慮做 weekly report 渲染替代手寫 .md（中長期）；Sumatra 短期可仿其 metadata schema
- **既有 cross-ref**：`provenance-tier-audit` skill、evidence_ledger、weekly-report v2

### C-5 Turing Way reproducibility checklist
- URL: <https://book.the-turing-way.org/reproducible-research/reproducible-research/> / <https://the-turing-way.netlify.app/reproducible-research/rdm/rdm-checklist.html> / <https://github.com/the-turing-way/the-turing-way>
- 適用範圍：data science 可複現性 community handbook
- **可借鏡點**：
  - 「don't touch raw data, keep read-only backup」— 對應 InterSubMod canonical 19-runs 的不可改原則
  - 「another researcher from outside the group should reproduce」— 對 conclude-research 的 sensitivity 6/6 可作為「外部複現性 proxy」
  - checklist 條列式格式 — 可直接抄入 `validation-protocol` skill
- **採納建議**：B 部分採。checklist 結構好抄；hypothesis registry 部分搜尋未明確找到，需另查
- **既有 cross-ref**：`validation-protocol`、canonical 19-runs、`feedback_evidence_driven_iteration_workflow`

---

## §5 Category D — 基礎設施問題（最高 P1）

對應 InterSubMod 痛點：disk monitoring 無 / TMPDIR 未強制 export / KB freshness 14d 警告但無強制重認證 / expected_coverage hardcoded。**2026-05-08 已造成 800 GB disk 災情**。

### D-1 Prometheus node_exporter textfile collector
- URL: <https://github.com/prometheus/node_exporter> / <https://prometheus.github.io/client_python/exporting/textfile/> / <https://omarghader.github.io/monitor-disk-space-node-exporter-prometheus-grafana/>
- 適用範圍：通用 host metrics 暴露
- **可借鏡點**：
  - **textfile collector** = `--collector.textfile.directory` 指向資料夾，cron / script 寫 `*.prom` 檔，node_exporter 自動讀取暴露
  - **atomic write** pattern：先寫 temp file → mv 到目標 — 防止 partial read（也是 InterSubMod /tmp 災情可借鏡的 atomic write 設計）
  - 適合監控 cron job / 自定義 metric / 不適合長 daemon
  - **直接設計片段**：
    ```bash
    # /etc/cron.d/disk_monitor → 每 10 分鐘
    df --output=source,pcent,target | awk 'NR>1 {gsub("%","",$2); print "disk_use_pct{mount=\""$3"\"} " $2}' > /tmp/disk.prom.$$
    mv /tmp/disk.prom.$$ /var/lib/node_exporter/textfile/disk.prom
    ```
- **採納建議**：A 直接採。即使無 Prometheus server，textfile 模式仍可退化為「純 logfile + grep + alert」，不損失可觀測性
- **既有 cross-ref**：`feedback_tmp_disk_full_pipeline_pitfall`、`scripts/analysis/check_ai_agent_readiness.sh`

### D-2 Bash + cron / systemd timer 磁碟告警
- URL: <https://www.cyberciti.biz/tips/shell-script-to-watch-the-disk-space.html> / <https://tecadmin.net/shell-script-to-check-disk-space-and-send-alert/> / <https://computingforgeeks.com/systemd-timers-linux/>
- 適用範圍：簡易 host-level 告警，無監控基建
- **可借鏡點**：
  - 二段門檻：warn 80% / critical 90% — 對應 InterSubMod 800 GB 災情前可早 12-24 小時發現
  - **systemd timer 比 cron 強的點**：journal log / 補跑 missed runs / `RandomizedDelaySec` 防 thundering herd / `systemctl start` 立即測試
  - cron `@daily` minimal env 陷阱 — 必用絕對路徑（對應 InterSubMod hook 的注意事項）
  - **直接設計片段**：
    ```bash
    # /usr/local/bin/disk_guard.sh
    USAGE=$(df / | awk 'NR==2 {gsub("%","",$5); print $5}')
    [ "$USAGE" -ge 90 ] && logger -t disk_guard "CRITICAL: / at ${USAGE}%" && \
        echo "DISK CRITICAL ${USAGE}%" >> /big7_disk/.disk_alert.log
    [ "$USAGE" -ge 95 ] && touch /big7_disk/.PIPELINE_BLOCKED  # 給 hook 讀
    ```
- **採納建議**：A 直接採。配 PreToolUse hook 讀 `.PIPELINE_BLOCKED` 即可硬阻擋長 pipeline
- **既有 cross-ref**：`feedback_tmp_disk_full_pipeline_pitfall`、PreToolUse hook 模式、`.claude/settings.local.json`

### D-3 Snakemake `resources.tmpdir` + wrapper 注入
- URL: <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html> / <https://github.com/snakemake/snakemake/issues/1059>
- 適用範圍：HPC pipeline TMPDIR 強制注入
- **可借鏡點**：
  - `resources: tmpdir="/big7_disk/tmp"` 自動 export `$TMPDIR` 到 shell / wrapper / notebook
  - delayed evaluation — cluster job 執行時才確定 path
  - **教訓**：Snakemake issue #1059 顯示 wrapper 不一定遵守 — InterSubMod 必須**自己驗證** longphase-to / clairs-to 等 caller wrapper 是否真的吃 `$TMPDIR`，不能信賴 wrapper 的 doc
  - **直接設計片段**：所有 pipeline launcher 第一行
    ```bash
    export TMPDIR="${TMPDIR:-/big7_disk/tmp}"
    mkdir -p "$TMPDIR"
    [ "$(df / --output=avail | tail -1)" -lt 10000000 ] && { echo "ABORT: / has <10GB"; exit 1; }
    ```
- **採納建議**：A 直接採。可寫成 `scripts/lib/tmpdir_guard.sh` 給所有 launcher source
- **既有 cross-ref**：`feedback_tmp_disk_full_pipeline_pitfall`、`/big7_disk/liaoyoyo2001/InterSubMod/scripts/`

### D-4 FAIR data principles（freshness 概念補強）
- URL: <https://www.go-fair.org/fair-principles/> / <https://www.nature.com/articles/sdata201618> / <https://book.the-turing-way.org/reproducible-research/rdm/rdm-fair/>
- 適用範圍：scientific data stewardship 通用原則
- **可借鏡點**：
  - 「machine-actionability」— metadata 必須機器可讀，對應 InterSubMod KB 的 `freshness_class` / `last_verified` field
  - **R = Reusable** 強調「sufficiently annotated」— 對應 KB 14d freshness 警告但無強制重認證的缺口：應補 `last_verified` field 過期時自動重認證
  - **注意**：FAIR 文獻**未明確定義 freshness/reverification 機制** — 屬於 InterSubMod 自己要定義的領域，FAIR 只給原則
- **採納建議**：B 部分採。原則層引用，實作仍需 InterSubMod 自定
- **既有 cross-ref**：KB `last_verified` field、`/check-staleness` skill

### D-5 Pachyderm container-stage versioning
- URL: <https://github.com/pachyderm/pachyderm> / <https://atlan.com/pachyderm-data-lineage/>
- 適用範圍：containerized pipeline + auto data versioning
- **可借鏡點**：每 stage 有獨立 container + git-like commit on data；config drift 不會發生
- **採納建議**：C 不採（理由：對 single-researcher 過重 + InterSubMod C++ 核心已有 build hook，無需 container 層）。但「config drift impossible by design」概念可記住
- **既有 cross-ref**：N/A

---

## §6 Category E — 報告呈現問題

對應 InterSubMod 痛點：matplotlib CJK 字型亂碼 / PPTX font fallback / 圖片擠壓 / caption checklist 5 層尚未實作（anchor #6 PARKED）。

### E-1 Matplotlib font fallback（>=3.6）
- URL: <https://matplotlib.org/stable/users/explain/text/fonts.html> / <https://github.com/matplotlib/matplotlib/pull/20740> / <https://matplotlib.org/stable/gallery/text_labels_and_annotations/font_family_rc.html>
- 適用範圍：matplotlib 多語言字型
- **可借鏡點**：
  - **3.6+ 已支援** Agg/SVG/PDF/PS backend per-glyph fallback — 不需手寫 fallback 邏輯
  - 推薦 rcParams：
    ```python
    import matplotlib
    matplotlib.rcParams['font.family'] = ['sans-serif']
    matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Noto Sans CJK TC', 'Source Han Sans TW', 'Droid Sans Fallback']
    ```
  - 「good practice: append generic-family as last resort」
  - **已知 bug**（issue #29173）：math formula + CJK 混排仍有問題 — 提醒 InterSubMod 圖表盡量避免 LaTeX math + 中文同時出現
- **採納建議**：A 直接採。可在 `scripts/lib/plot_setup.py` 集中設定 rcParams，所有 plot 強制 import
- **既有 cross-ref**：`feedback_matplotlib_cjk_font_rule`、`feedback_pptx_screenshot_rendering_rules`

### E-2 Quarto + pandoc-crossref figure caption
- URL: <https://quarto.org/docs/authoring/figures.html> / <https://quarto.org/docs/authoring/cross-references.html> / <http://lierdakil.github.io/pandoc-crossref/>
- 適用範圍：scientific 文件 figure 自動編號與 cross-ref
- **可借鏡點**：
  - `#fig-xxx` label 強制命名規範 — 對應 anchor #6 caption checklist 的「label 必須有結構性前綴」
  - cross-reference 必須以 type 開頭（fig- / tbl- / eq-）— 直接 lint friendly
  - 編號自動 — 不會錯
- **採納建議**：B 部分採。InterSubMod 短期不轉 Quarto，但 caption checklist 5 層可仿 Quarto 的「title / type / label / units / source」結構
- **既有 cross-ref**：anchor #6 PARKED、5 層圖表 caption

### E-3 Cell STAR Methods caption 規範
- URL: <https://www.cell.com/information-for-authors/star-authors-guide> / <https://www.cell.com/pb-assets/journals/research/cell/methods/Methods_Guide_Cell-1724873926857.pdf>
- 適用範圍：Cell Press 期刊強制統一格式
- **可借鏡點**：每個 figure caption 必含：(1) declarative title (2) brief methods (3) statistical details: tests / n / center definition / dispersion measures (4) p-value asterisks 規範
- **採納建議**：A 直接採。InterSubMod 的 5 層 caption checklist 可對齊：title / methods / n+samples / stat test + p / source ledger ID
- **既有 cross-ref**：anchor #6、`structured-tech-report` skill、weekly-report v2

### E-4 Caltech writing center figure caption guide
- URL: <https://writing.caltech.edu/documents/27629/HWC-FigureCaptionHandout.1-2024.pdf> / <https://www.aje.com/arc/writing-effective-figure-legend>
- 適用範圍：通用 scientific 寫作指南
- **可借鏡點**：caption「title-then-detail」結構；caption 必須 stand alone（不依賴 main text 也能懂）
- **採納建議**：A 直接採。可作為 `weekly-report` skill 的 figure 寫作 reference 文檔
- **既有 cross-ref**：weekly-report v2、`feedback_pptx_term_explanation_rule`

### E-5 Diátaxis documentation framework
- URL: <https://diataxis.fr/> / <https://diataxis.fr/start-here/>
- 適用範圍：技術文件分類框架
- **可借鏡點**：
  - 4 種文件類型：Tutorial / How-to / Reference / Explanation — 各自有不同需求結構
  - 「action vs knowledge × study vs work」矩陣 — 對應 InterSubMod 的 docs/ 子目錄分類（references / guides / reports / experiments）
- **採納建議**：B 部分採。InterSubMod docs/ 已部分對齊；可在 `docs/README.md` 顯式標註每子目錄屬於哪 1-2 個 Diátaxis 類型，避免分類混亂
- **既有 cross-ref**：`docs/README.md`、`/doc-standards` skill

---

## §7 整合建議：對 5 個 P1 提案的最直接 prior art 對應

複盤報告列出的 5 個 P1 提案，逐一給「最接近 prior art + 可直接抄的設計片段」：

### P1-1 pitfalls table review（每週 / 每 cycle）
- **最近 prior art**：Google SRE blameless postmortem trend analysis（<https://sre.google/sre-book/postmortem-culture/>）+ DECLARE checklist（<https://pmc.ncbi.nlm.nih.gov/articles/PMC10149772/>）
- **可抄設計**：
  - SRE：「standard postmortem template to consistently capture root cause and trigger data, enabling trend analysis to target improvements addressing systemic root-cause types」— 對應 pitfalls table 應分類「systemic root-cause type」（不是平鋪 N 條）
  - DECLARE 7 字 checklist（Differential / Evidence / Confounders / Likelihood / Alternatives / Reassess / Engage）— 直接 7 行 review 表單
  - **action**：每週週報加「pitfalls trend」一節，標出本週新增 pitfall + 已知 pitfall 觸發次數
- **既有**：pitfalls_table 6 條、weekly-report v2

### P1-2 plan template（cycle-init 強制 conditional plan）
- **最近 prior art**：OSF Preregistration（<https://www.cos.io/initiatives/prereg>）+ Anthropic 6 patterns（<https://resources.anthropic.com/building-effective-ai-agents>）+ Devin plan checkpoint（<https://sureprompts.com/blog/devin-ai-prompting-guide>）
- **可抄設計**：
  - OSF conditional plan：「if normality 不符 → 改用 non-parametric」結構 — 直接抄 cycle-init plan template 的 fail-mode 預設
  - Devin：「plan as the most valuable checkpoint, reviewed plan catches drift in 1 min」— 把 plan review 設為 L4 mandatory
  - **action**：`cycle-init` skill 強制 4 段：(1) hypothesis + sensitivity threshold (2) confound list + sweep plan (3) conditional analysis branches (4) expected fail mode + fallback
- **既有**：`cycle-init` skill、Drill 1 retrospective sensitivity 6/6 / specificity 2/2

### P1-3 pitfalls 擴充（merge dataset / TMPDIR / caller_af 等）
- **最近 prior art**：Great Expectations Expectation Suite（<https://greatexpectations.io/>）+ dbt source freshness（<https://docs.getdbt.com/docs/deploy/source-freshness>）
- **可抄設計**：
  - Great Expectations：每個 pitfall 寫成可執行 expectation（如 `expect_column_quantile_values_to_be_between(caller_af, 0.75, [0.05, 0.95])`）— 自動化驗證
  - dbt freshness：`warn_after` / `error_after` 二段門檻語法
  - **action**：把 6 條 pitfalls 改寫成 YAML / JSON spec，driver 載入時自動跑；fail 即 abort
- **既有**：`feedback_merged_dataset_af_and_loh_pitfalls`、`feedback_tmp_disk_full_pipeline_pitfall`、`feedback_feature_name_vs_definition_rule`

### P1-4 disk health monitor
- **最近 prior art**：Prometheus textfile collector（<https://prometheus.github.io/client_python/exporting/textfile/>）+ bash df cron（<https://www.cyberciti.biz/tips/shell-script-to-watch-the-disk-space.html>）+ systemd timer（<https://computingforgeeks.com/systemd-timers-linux/>）
- **可抄設計**：
  ```bash
  # /usr/local/bin/disk_guard.sh （每 10 min systemd timer）
  ROOT_USAGE=$(df / | awk 'NR==2 {gsub("%","",$5); print $5}')
  TMP_USAGE=$(df /tmp | awk 'NR==2 {gsub("%","",$5); print $5}')
  TS=$(date -Iseconds)
  echo "${TS} root=${ROOT_USAGE} tmp=${TMP_USAGE}" >> /var/log/disk_guard.log

  if [ "$ROOT_USAGE" -ge 90 ] || [ "$TMP_USAGE" -ge 90 ]; then
      logger -t disk_guard "WARN: root=${ROOT_USAGE}% tmp=${TMP_USAGE}%"
      touch /big7_disk/.DISK_WARN
  fi
  if [ "$ROOT_USAGE" -ge 95 ] || [ "$TMP_USAGE" -ge 95 ]; then
      touch /big7_disk/.PIPELINE_BLOCKED  # PreToolUse hook 讀取硬阻擋
  fi
  ```
- **既有**：`feedback_tmp_disk_full_pipeline_pitfall`、`scripts/analysis/check_ai_agent_readiness.sh`

### P1-5 TMPDIR wrapper
- **最近 prior art**：Snakemake `resources.tmpdir` delayed evaluation（<https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html>）+ Snakemake issue #1059 教訓（<https://github.com/snakemake/snakemake/issues/1059>）
- **可抄設計**：
  ```bash
  # scripts/lib/tmpdir_guard.sh （所有 pipeline launcher 第一行 source）
  set -e
  : "${TMPDIR_BASE:=/big7_disk/tmp}"
  export TMPDIR="${TMPDIR_BASE}/$(whoami)/$$"
  mkdir -p "$TMPDIR"
  trap 'rm -rf "$TMPDIR"' EXIT

  # Pre-flight check
  AVAIL_KB=$(df --output=avail "$TMPDIR_BASE" | tail -1)
  [ "$AVAIL_KB" -lt 10000000 ] && { echo "ABORT: $TMPDIR_BASE has <10GB free"; exit 1; }

  # 驗證 wrapper 真的吃 TMPDIR（issue #1059 教訓）
  echo "TMPDIR=$TMPDIR" >&2
  ```
- **教訓**：Snakemake issue #1059 顯示「設了 tmpdir resource 但 wrapper 不一定吃」— InterSubMod 必須對 longphase-to / clairs-to / 其他 caller 各別驗證一次
- **既有**：`feedback_tmp_disk_full_pipeline_pitfall`

---

## §8 References

### Category A
1. [DVC Pipelines documentation](https://doc.dvc.org/user-guide/pipelines)
2. [DVC repro command reference](https://doc.dvc.org/command-reference/repro)
3. [dbt source freshness](https://docs.getdbt.com/docs/deploy/source-freshness)
4. [dbt freshness reference](https://docs.getdbt.com/reference/resource-properties/freshness)
5. [dbt model contracts (Medium)](https://medium.com/@faouaghzene/dbt-model-contracts-beyond-tests-toward-robust-governance-6814c1eda642)
6. [Great Expectations](https://greatexpectations.io/)
7. [Great Expectations schema validation](https://docs.greatexpectations.io/docs/reference/learn/data_quality_use_cases/schema/)
8. [Pachyderm GitHub](https://github.com/pachyderm/pachyderm)
9. [lakeFS data lineage blog](https://lakefs.io/blog/data-lineage-for-data-lakes/)
10. [Data Contract CLI](http://cli.datacontract.com/)

### Category B
11. [nf-core rnaseq pipeline](https://github.com/nf-core/rnaseq)
12. [nf-core community paper (bioRxiv)](https://www.biorxiv.org/content/10.1101/610741v3.full)
13. [MultiQC GitHub](https://github.com/MultiQC/MultiQC)
14. [MultiQC paper (Bioinformatics 2016)](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)
15. [Anthropic Building Effective Agents](https://www.anthropic.com/research/building-effective-agents)
16. [Anthropic cookbook evaluator-optimizer](https://github.com/anthropics/anthropic-cookbook/blob/main/patterns/agents/evaluator_optimizer.ipynb)
17. [Anthropic demystifying evals](https://www.anthropic.com/engineering/demystifying-evals-for-ai-agents)
18. [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices)
19. [OSF Preregistration](https://www.cos.io/initiatives/prereg)
20. [OSF preregistration help](https://help.osf.io/article/330-welcome-to-registrations)

### Category C
21. [Devin AI Prompting Guide](https://sureprompts.com/blog/devin-ai-prompting-guide)
22. [Devin vs Cursor (builder.io)](https://www.builder.io/blog/devin-vs-cursor)
23. [Cursor vs Devin vs Claude Code (BattleAITools)](https://battleaitools.com/ai-agents/cursor-agent-vs-devin-vs-claude-code/)
24. [Anthropic 6 composable agent patterns (aimultiple)](https://aimultiple.com/building-ai-agents)
25. [Cognitive forcing tool RCT (BMC Med Educ 2018)](https://bmcmededuc.biomedcentral.com/articles/10.1186/s12909-018-1444-3)
26. [DECLARE strategy (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10149772/)
27. [Cognitive debiasing strategies (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3786644/)
28. [Quarto](https://quarto.org/about.html)
29. [Sumatra docs](https://sumatra.readthedocs.io/en/latest/)
30. [Sumatra GitHub](https://github.com/open-research/sumatra)
31. [Turing Way reproducible research](https://book.the-turing-way.org/reproducible-research/reproducible-research/)

### Category D
32. [Prometheus node_exporter GitHub](https://github.com/prometheus/node_exporter)
33. [Node exporter textfile collector](https://prometheus.github.io/client_python/exporting/textfile/)
34. [Disk monitor with node_exporter (Omar Ghader)](https://omarghader.github.io/monitor-disk-space-node-exporter-prometheus-grafana/)
35. [Bash df disk alert script (cyberciti)](https://www.cyberciti.biz/tips/shell-script-to-watch-the-disk-space.html)
36. [Systemd timers guide](https://computingforgeeks.com/systemd-timers-linux/)
37. [Snakemake rules + tmpdir](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)
38. [Snakemake tmpdir wrapper bug #1059](https://github.com/snakemake/snakemake/issues/1059)
39. [FAIR principles (Nature 2016)](https://www.nature.com/articles/sdata201618)
40. [GO FAIR principles](https://www.go-fair.org/fair-principles/)

### Category E
41. [Matplotlib fonts explain](https://matplotlib.org/stable/users/explain/text/fonts.html)
42. [Matplotlib font fallback PR #20740](https://github.com/matplotlib/matplotlib/pull/20740)
43. [Matplotlib CJK math bug #29173](https://github.com/matplotlib/matplotlib/issues/29173)
44. [Quarto figures](https://quarto.org/docs/authoring/figures.html)
45. [Quarto cross-references](https://quarto.org/docs/authoring/cross-references.html)
46. [pandoc-crossref](http://lierdakil.github.io/pandoc-crossref/)
47. [Cell STAR Methods authors guide](https://www.cell.com/information-for-authors/star-authors-guide)
48. [Caltech figure caption handout](https://writing.caltech.edu/documents/27629/HWC-FigureCaptionHandout.1-2024.pdf)
49. [AJE figure legend guide](https://www.aje.com/arc/writing-effective-figure-legend)
50. [Diátaxis framework](https://diataxis.fr/)

### Cross-cutting
51. [Google SRE postmortem culture](https://sre.google/sre-book/postmortem-culture/)
52. [Google SRE workbook postmortem analysis](https://sre.google/workbook/postmortem-analysis/)

---

*報告完。所有 prior art 採納前需在 InterSubMod 上下文 evaluate；建議下一步：先驗證 Category D 的 disk_guard.sh + tmpdir_guard.sh（最高 P1 + 最低成本），再評估 Category A 的 dbt freshness / Great Expectations YAML 是否值得引入。*
