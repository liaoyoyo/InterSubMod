---
title: 數據誠信三層防捏造系統
date: 2026-06-01
type: migration / harness-design
task_type: E hotfix (流程缺陷修補) → 機制落地
data_sources: docs/postmortems/20260601_fabricated_metric_in_html_preview_postmortem.md
status: landed
---

# 數據誠信 — 三層防捏造系統（2026-06-01）

> **一句話**：把「報告數字必有檔案來源」從失效的純文字規則，升級成 **由構造防止 + 寫入時 gate + 任務結束溯源** 三層機械防線，平衡 真實 × 消耗 × 錯誤率。
> **觸發**：`20260601_fabricated_metric_in_html_preview_postmortem.md`（AI 把預期 AUC 當真值寫進 PI HTML，方向還相反）。
> **設計約束（用戶）**：在「抓到捏造(真實) / token-latency(消耗) / 誤報漏報(錯誤率)」三軸平衡。

---

## 1. 為什麼純規則不夠（先講結論）

postmortem §9 是整件事的核心：寫完 memory + postmortem + pitfall 卡之後 **30 分鐘內又捏造一次**。

| 防線 | 對「憑預期打數字」有效？ |
|------|:---:|
| 文字規則（memory / CLAUDE.md / pitfall）| ❌ 實測失效兩次 |
| 模型自律（Opus 4.8 literal）| ❌ literal 只保證不泛化指令，不涵蓋「不憑預期打數字」的衝動 |
| 同回合自我比對 | 🟡 有效但太晚（已寫進檔案）|
| **機械 gate / 由構造防止** | ✅ 唯一可靠 |

→ 方向鎖定：要在「打字進檔案」這步機械攔截，或乾脆讓手打不可能。

## 2. 業界方法對三軸的打分（2026 web 調查）

| 方法 | 代表 | 真實 | 消耗 | 錯誤率(對數字) | 採用 |
|------|------|:---:|:---:|:---:|:---:|
| **由構造防止** | Quarto / RMarkdown inline code（數字 render 時算出，不手打）| ★★★ | ★★★ | ★★★ | ✅ Layer A |
| **來源歸因** | evidence_ledger / `[verified:path:line]` / CiteGuard / PaperTrail(CHI 2026) | ★★ | ★★ | ★ 可調 | ✅ Layer B/C |
| **LLM-judge 偵測** | HalluMat / RAGAS / NLI fact-check | ★ | ✗ 每數字一次 call | ✗ EMNLP 2025：對數字/量詞最不準 | ❌ 不採用 |

關鍵：LLM-judge 偵測（業界主流）正好三軸最差，且對「數字」這類 claim 天生不準 → **不裝**，符合 restraint。

## 3. 三層防線（已落地）

### Layer A — 由構造防止（最優先）
- 工具：`scripts/fill_report.py <template> <data.json> -o <out>`
- 機制：template 用 `{{dotted.key}}` 佔位，render 時從 data.json 注入；**任一 key 缺失 → refuse 不寫出**（exit 1）。
- 保證：報告每個數字 = `data.json[key]`，物理上無法手打一個不在資料檔的值。
- 證據：`scripts/harness_health.py` 就是這模式（grep 磁碟→render），從沒捏造過。

### Layer B — 寫入時 gate（backstop）
- Hook：`scripts/hooks/number_provenance_check.sh`（PreToolUse Edit|Write `.md`/`.html`）
- 抽「metric 形數字」(`AUC=`/`p=`/`%`/`Δ`/`≥2 位小數`) → 去 bounded 來源（frontmatter `data_sources:` + 同層 `_assets/` + 同目錄 `.json|.tsv|.txt|.csv`）grep。
- **分級**：`validated/`、`pi_reports/` 無來源 → **exit 2 阻擋**；其他路徑 → advisory exit 0 提醒。
- **fail-OPEN**：python 缺/解析錯/無來源 → exit 0（broken hook 絕不擋所有 Write；與 `|| exit 0` neutering bug 本質不同 — 真偵測在 strict 路徑仍 exit 2）。
- override：內文 `<!-- provenance-verified: 理由 -->`。

### Layer C — 任務結束溯源表（紀錄依據）
- 工具：`python3 scripts/number_provenance.py audit <report> [--sources P ...]`
- 產 markdown 表：`| metric | 檔案:行 | 狀態 |`，同時就是「每個數字的依據」紀錄。
- 時機：validated / PI / handoff 收尾前跑一次（低頻 → 成本攤平）。

## 4. 三軸如何平衡

| 軸 | 怎麼顧 |
|----|--------|
| **真實** | A 由構造杜絕；B 抓手打捏造（validated 強制擋）；C 邊界紀錄 |
| **消耗** | 全確定性 grep，flat-rate≈0 token；bounded 來源集（§12 不掃 repo-root）；hook 非 .md/.html 早退 |
| **錯誤率** | 只抓 metric 形（跳年份/日期/章節號/裸整數）降假陽性；advisory 預設不擋降摩擦；衍生數字假陰性靠 A 補；substring 比對容忍尾數精度（report `0.866` 命中 source `0.8662`）|

## 5. 驗收

- `number_provenance_check.sh`：11/11 整合測試通過（非 doc 路徑/in_progress advisory/validated block/override/empty/fail-open/Edit/Read 忽略）。
- `fill_report.py`：注入 OK；缺 key 拒寫 OK。
- `number_provenance.py audit`：sourced/unsourced 分類 OK。
- hook 數 38→39；CLAUDE.md §4 + §13 + memory `feedback_no_fabricated_numbers_in_reports` 已同步。

## 6. 落地的 postmortem action items

A4（grep 自審）→ Layer B；A5（hook 從待議升必做）→ `number_provenance_check.sh`；A7（read-back-then-fill）→ Layer A `fill_report.py` + §13 鐵則。

## 7. 殘餘限制（誠實揭露）

- B 層假陰性：心算的衍生數字（和/均值）grep 不到 → 只能靠 A 層（衍生值也由 render 算出）。
- B 層假陽性：sign 不符（report `-0.122` vs source `0.122`）會 flag，advisory surface，可接受。
- C 層非強制自動跑（避免 Stop hook 噪音）→ 靠收尾紀律 + 文件提醒。

> 來源：事件與真值引自 `data_sources` 宣告的 postmortem（H19=0.985 捏造 / BRCA2 真值 0.866 / GNAS 0.567 / unphase 45.84%）。業界方法引自 2026 web 調查（Quarto literate programming、CiteGuard、PaperTrail CHI 2026、EMNLP 2025 hallucination-detection 再評估）。
