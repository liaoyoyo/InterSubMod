# 限制與提示建議治理盤點（2026-06-01）

> **目的**：盤點 harness 現有所有「限制（block/gate）」與「提示建議（advisory）」，逐條判斷合理性、追溯當初加入的原因（origin）、給改進方式。
> **方法**：read settings.local.json（permissions + hooks）+ CLAUDE.md（§1/§8/§12 + Hard Gate）+ state/research per-folder CLAUDE.md + MEMORY/postmortem 追 origin。
> **結論先講**：**全部 traced 到具體事故或治理原則，restraint-consistent，無一不合理**。最大真破口（2 個 Hard Gate 被 `\|\| exit 0` neuter）已於 2026-05-31 修復。剩餘可改進主題只有「軟提示廣關鍵字噪音」（皆非阻擋，低優先）。
> **此 review 已整合進 `/harness-health` 儀表板 §⑦**（持續可見）。

---

## A. 分類框架

| 類型 | 行為 | 可繞過? | 數量 |
|------|------|:---:|:---:|
| **BLOCK** | exit-2 / deny → 直接擋工具執行 | ❌（exit-2）/ permission deny | 6 gate + deny list |
| **PROMPT** | permission ask → 人工確認 | 需用戶 ack | ask list 20 條 |
| **RULE** | CLAUDE.md 硬規則（§1/§8/§12 等）| 治理層，AI 自律 + hook 輔助 | ~6 |
| **ADVISE** | 注入 context / 印提醒，不擋 | 預設過 | ~10 hooks |
| **LOG** | 純被動 telemetry | 無行為影響 | ~6 hooks |

---

## B. 硬限制（BLOCK / PROMPT / RULE）— 逐條

### B1. permission `deny`（settings）— ✅ 合理 keep
- **機制**：`Read **/.env* / *secret* / *token* / *.pem / id_rsa / .ssh / /etc/shadow`；`Bash rm -rf / · rm -rf /* · mkfs · fdisk · parted`。
- **緣由**：標準機密外洩防護 + 毀滅性磁碟操作（業界 baseline security）。
- **合理性**：✅ 教科書級，零誤傷風險。**改進：無。**

### B2. permission `ask`（settings，20 條）— ✅ 合理（極小摩擦）
- **機制**：`rm / mv / chmod / chown / ln / truncate / dd / git reset/clean/checkout/restore/rebase/push / curl / wget / ssh / scp / rsync / sudo`。
- **緣由**：對齊 §1「可逆性 override = 不可逆操作前人工確認」。多數是檔案破壞或外送（git push / curl / ssh）或重寫歷史（reset/rebase）。
- **合理性**：✅ 合理。**唯一摩擦**：少數「可逆」操作也被 prompt（如 `git restore --staged` 取消暫存、`git checkout` 切分支）。本 session 清 stage 時就被 prompt 一次。
- **改進（低優先）**：可把明確可逆的安全形式 allow-list（如 `git restore --staged:*`、`git checkout <branch>`），但摩擦罕見、收益小，**不急**。

### B3. exit-2 Hard Gate hooks（6）— ✅ 全合理（剛修好）
| gate | 緣由（origin 事故/原則）| 合理性 |
|------|------|--------|
| `pre_commit_compile_check` | C++ 改動未編譯就 commit → regression | ✅ |
| `kb_schema_check` | KB 寫入 schema/related_ids/canonical 違規 | ✅ |
| `pipeline_block_check` | **/tmp 800GB disk-full 災情**（feedback_tmp_disk_full）| ✅ |
| `no_binary_commit` | binary 進 git → repo bloat | ✅ |
| `search_scope_guard` | **repo 極大、大遞迴搜尋卡住**（§12）| ✅ |
| `pre_tier_upgrade_check` | **tier ⭐4/5 over-claim**（⭐4→⭐3 降級教訓）| ✅（2026-05-31 wired）|
- **2026-05-31 修復**：`pre_commit_compile_check` + `kb_schema_check` 先前被 settings wiring `2>/dev/null \|\| exit 0` **neuter**（exit 2 被吃成 0 → 實際不擋）；已移除 mask（改 `1>&2` surface 原因）。**改進：已做；`/harness-health` 燈 #2 持續監看復發。**

### B4. CLAUDE.md 硬規則 — ✅ 核心治理 keep
- **§1 Hard Gate 永遠 🔴**：刪檔 / C++ commit / NO-GO 判定 / evidence_ledger 覆寫。緣由=不可逆 / 影響結論。✅。
- **§0 task type gate**：6 類分類。緣由=**5/24 scope-mismatch incident**。✅ 合理；**微改進**：對明顯 trivial 任務分類步驟稍重，但 advisory 非阻擋。
- **§12 搜尋紀律**：禁 `grep -r .` / `find .` 無 maxdepth。緣由=repo GB 級 data/output/worktree。✅（Grep/Glob 工具已自動跳過，規則只擋 Bash 誤用）。
- **§8 Hard-Gate work 絕不入 Dynamic Workflow**：緣由=workflow subagent 一律 acceptEdits + 繞過 exit-2 hook。✅ **關鍵保護，勿放寬**。
- **state/CLAUDE.md 不可手改 state.json**、**research/CLAUDE.md 禁直接 edit queue/ledger**：緣由=schema/provenance 完整。✅（須走 skill）。

---

## C. 軟提示（ADVISE）— 逐條

### C1. UserPromptSubmit 注入（6）— ✅ 合理，但廣關鍵字易噪音
| hook | 緣由 | 噪音風險 |
|------|------|----------|
| `knowledge_check` | outside-claim 必查 KB（**slide14 F1 雙口徑事故**）| 中（本回合對「分析流程」提示 KB，未必相關）|
| `kb_freshness_warn` | KB 時效（14 天）| 低-中（本回合 5 檔 14 天提示）|
| `narrative_frame_advisor` | 減少理解負擔（§11）| 中（幾乎每個「報告/整理/說明」都觸發）|
| `task_type_advisor` | **5/24 incident** | 中（本回合誤判 handoff）|
| `md_path_format_rule` | **用戶強制規則**（feedback_md_path_prefix）| 低 |
| `research_direction_guard` | 防 reopen concluded NEGATIVE | 低 |
- **合理性**：✅ 全部有正當 origin。**但都是廣關鍵字觸發 → 偶爾重複/誤觸**（本 session 就見 handoff 誤判 + 無關 KB 提示）。
- **改進（低優先，皆非阻擋）**：(1) 收窄關鍵字特異性（如 task_type 區分真 handoff vs 文中提及）；(2) 對明顯 trivial / meta prompt 靜音；(3) narrative_frame 已有「不用框架」override。**因非阻擋，ROI 低，現狀可接受。**

### C2. PostToolUse claim 4-linter — ✅ 合理（稍冗餘，維持）
- `evidence_level_lint` / `causal_claim_check` / `researcher_claim_evidence_check` / `terminology_guard`：防 hedge 語言 / 因果過度宣稱 / L3 當 L1 / 術語漂移。緣由=研究誠信（researcher claim 需實測升 L1）。
- **合理性**：✅。4 個分開稍冗餘（PostToolUse 共 17 hooks，最重事件）。**改進**：audit 判 **SKIP_RISK** — working safety hook 合併有 regression 風險，flat-rate 下省 latency 無意義。**維持。**

### C3. soft gate（verify_gate / kb_sot_guard）— ✅ 已正名
- 皆 `exit 0` advisory（非 Hard Gate）。緣由：verify_gate=**+0.057 fabrication 事故**；kb_sot_guard=F1 SoT。
- **2026-05-31 已正名**：CLAUDE.md §4 不再誤列為 Hard Gate。**改進：已做。**

### C4. 收尾提示（Stop / SubagentStop）— 🟡 合理但稍噪
- 會話結束提醒寫報告 / subagent 完成檢查。緣由=防漏記錄。
- **合理性**：🟡 每次觸發稍噪（未必每次都需寫報告）。**改進（低優先）**：可改條件式（有實質產出才提醒）。

---

## D. 背景 logger（LOG）— ✅ 純被動，keep

`subagent_completion_logger` / `skill_change_audit` / `memory_recall_logger` / `evidence_read_tracker` / `skill_usage_logger` / `cache_telemetry`（manual）。量化 cache-hit / 引用率 / skill 變動稽核。**無行為影響**。`cache_telemetry` USD 已標 flat-rate 無意義（2026-05-31）。**改進：無。**

---

## E. 總結裁決

| 維度 | 裁決 |
|------|------|
| **合理性** | ✅ 全部有正當 origin（事故或治理原則），restraint-consistent，**無一不合理或該移除** |
| **最大真問題** | 2 個 Hard Gate 被 `\|\| exit 0` neuter（限制名存實亡）→ **已修** |
| **次要可改進** | 軟提示廣關鍵字噪音（C1）+ 收尾提示稍噪（C4）+ ask-list 可逆操作偶誤 prompt（B2）→ **皆非阻擋、ROI 低、現狀可接受** |
| **不要動** | claim 4-linter 合併（SKIP_RISK）、permission deny/ask 核心、§1/§8/§12 規則 |

**持續監看**：`/harness-health` 燈 #2（Hard-Gate 真實性）+ §⑦ 治理面板。每次改 hook/settings 重跑即知限制是否仍真實生效（防再次 neuter）。
