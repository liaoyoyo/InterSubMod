<!--
建立時間: 2026-06-15
報告類型: 工作流稽核 + 改進清單（最近研究流程/資料整理/運行）供 AI session 理解確認
任務類型: A pilot（稽核盤點，read-only + 改進建議；scope = 最近 arc + harness）
狀態: ✅ 盤點查驗完成；改進待用戶勾選才執行
data_sources: 主回合 2026-06-15 git log/status/check-ignore + ls -d（非 du）+ harness_health.py + hook 源碼讀取 + workflow 6-agent grounded 稽核（每條附指令證據）
build_branch: research/subclonal-reconstruction-202606
build_commit: 1cf5116
worktree: /big7_disk/liaoyoyo2001/InterSubMod （⚠ 6 worktree 並行，stamp 僅相對本 branch）
驗證方式: 6-agent workflow 每條發現附「實際指令+輸出」；主回合獨立複核 7 個高嚴重度爭議點（SC1-SC7），抓到 1 個 agent 建議錯誤（D2 chr2 歸屬）
-->

# 工作流稽核 + 改進清單 — 2026-06-15

> 觸發：用戶「確認最近研究流程/資料整理/運行是否有可改進，特別留意常見錯誤、高重複、可複利加速、高效率、減少風險」。
> 方法：Ultracode workflow 6 維度 grounded 稽核（每條附指令證據）+ 主回合獨立複核（§13.7 不採信單方回報）。

## L0 結論

- **整體健康良好**：harness 8 GREEN / 2 YELLOW / 0 RED；計數零漂移（47 skills / 18 agents / 42 hooks）；資料組織如實落地（archive 16 run、canonical 每樣本 2 run）；commit message 品質高、無 `git add -A` sweep。
- **三大機械錯誤防線運作中**：數字捏造（§13 三層）、跨 branch commit（trunk guard）、injection（sanitizer）皆已守。
- **🔴 最高價值缺口（1 件，本 session 親身踩到）**：`search_scope_guard.sh` 對 `big7_disk_output` 路徑的 `du/find/ls -R/wc -l` **全部放行** → 「卡住的動作」#1 主因零機械覆蓋。修復便宜（擴充既有 hook ~15 行）。
- **🔴 立即風險（不是我能修的）**：並行 session 的 chr2-18M 2291 行交付物未提交、共用 HEAD → loss risk。**我依 §C 不碰**，標示給你/原 session。
- **可複利加速**：4 個 anti-drift 腳本（provenance_stamp / number_provenance audit / fill_report / safe-ls）已存在但 manual-only / 零部署；archive 流程手做 2 次、還有 ≥3 批待做 → 值得腳本化。

## 主回合獨立複核對照（SC1-SC7）

| # | 複核項 | 指令 | 結果 | 對 agent 回報 |
|---|---|---|---|---|
| SC1 | D1-1 shared-state 被 commit？ | `git log -8 -- CURRENT_FOCUS.md` + `check-ignore` | 7/8 commit 含 CURRENT_FOCUS；check-ignore exit=1（未 ignore）| ✅ 屬實，但**規則矛盾**非活躍錯誤（見下）|
| SC2 | working tree 規模 | `git status --short \| uniq -c` | 36 ?? / 5 D / 26 M = 67 | ✅ 確認 |
| SC3 | D2-1 chr2 歸屬 | `ls -lt` mtime + `git status` | ?? 未提交、mtime 今天 05-17 點、**並行 session 作品** | 🔴 **修正 agent 建議**：不可由我 commit（§C）|
| SC4 | D5 harness | `harness_health.py` | 8G/2Y/0R、47/18/42/2 | ✅ 確認 |
| SC5 | D4 archive 落地 | `ls -d ... \| wc -l`（非 du）| 16 archived；HCC1395/COLO829 各 2 run | ✅ 確認 |
| SC6 | D3-1 hang guard | 灌 du/find/ls-R 進 hook + 讀源碼 | 全 exit=0 放行；hook 只擋 `find .`/`grep -r .`/`du */` | ✅ 確認**真缺口** |
| SC7 | D6-2 stamp 缺漏 | `grep -c build_branch:` audit 文件 | audit 文件 0、僅 3 manifest 有 | ✅ 確認 |

## 6 維度發現摘要

### D1 Git 衛生 — DONE_WITH_CONCERNS
- **D1-1/D1-2（高）規則矛盾**：§C「絕不 commit shared-state」vs §6「CURRENT_FOCUS+evidence_ledger 是權威 SoT（應版控）」。實況：shared-state 被版控+commit（多由別 session）。`git_branch_commit_guard.sh` 只擋 trunk，對 shared-state 零覆蓋。→ **reconcile 規則文字**（傾向承認 shared-state 該版控、但禁 `-A` sweep），非加 exit-2。
- **D1-3（中）**：`main` 停在 2026-03-14、落後 feature branch 479 commit；`develop` 才是實質整合 trunk（最新 06-14）。→ 文件講清楚 main=release-snapshot，否則被誤判 divergence。
- **D1-4（中）**：feature branch 比 develop 多 11 commit + 67 條 working tree → 大混合樹提高誤 `-A` sweep 風險（連動 D1-1）。
- **D1-5（低，優點）**：commit 嚴格 scoped、conventional-commit、無 sweep。

### D2 未提交風險 — DONE_WITH_CONCERNS
- **D2-1/D2-4（高，已修正歸屬）**：chr2-18M 2291 行（3 報告+concept+6 repro .py）全未提交、共用 HEAD（70 transcript session）→ loss risk。**但屬並行 session**，我依 §C「other-session 留用戶」**不 commit**；建議原 session/你處理（或開 worktree 隔離）。
- **D2-6（低，正確）**：shared-state（active.json/ledger/postmortem log）+ other-session churn（tsg_promoter / flagship_chr2 ism_out）已正確識別為「留著不動」。⚠ `research/flagship_chr2_.../ism_out` 未 gitignore + 體積大 → 若有人 `git add research/` 會誤提交。

### D3 常見錯誤（postmortem 挖掘）— DONE
- **D3-1（高）🔴 真缺口**：output-dir `du/find/ls -R/wc -l` 可 hang，`search_scope_guard` 完全沒擋（實測+源碼雙證）。→ **擴充既有 hook**：擋打在 `big7_disk_output`/`/output/` 且無 `-maxdepth` 的 `find|du|ls -R|wc -l`。
- **D3-2（高）已守**：數字捏造 → §13 三層（fill_report 構造層 + number_provenance_check exit-2 + audit 表）皆 wired。維護：light#2 盯 `|| exit 0` neutering 復發。
- **D3-3（中）半守**：output-token session 丟失 = friction #1；logger 有捕但**從不告警**（CLAUDE.md §8 承諾的 advisory 未 wire）。實測 7 個 session out>3M、單 session 最高 10.4M。→ wire 一行 threshold echo。
- **D3-4/D3-5（低）已守**：跨 branch commit（trunk guard + concurrent advisor）、injection（sanitizer 4 次皆 benign）。
- **D3-6（中）半守**：validated 報告 stale `last_read` 編輯（evidence_gate_bypass.log 有記但 advisory 不擋）。→ 可讓 number_provenance_check 順帶 flag stale last_read。

### D4 資料組織落地查驗 — DONE（全如實）
- **D4-1**：archive 16 run（5+11）磁碟=manifest 數學一致；mv 非刪、可還原。
- **D4-2**：canonical 每樣本 2 run（complete_matrix + B3-KEEP），2 個 KEEP 皆真被 `20260423_B3_paired_obs18.py:76/82` 硬引用。
- **D4-3**：data_specs guide+2 manifest+external-dependency 契約皆存在。
- **D4-4（中）**：新主軸輸出家 `research/subclonal_reconstruction/` 仍空骨架（只 00_INDEX.md）→ 第一個 cycle（G-A 5 樣本 V10）跑時要 mkdir + backfill，否則 AI 找不到會散落回 ad-hoc。
- **D4-5（低）**：archive 偏 HCC1395=7 已解釋（歷史 test/retry 最多、全 0 BAM）。

### D5 Harness 漂移 — DONE_WITH_CONCERNS
- **D5-1/D5-6（低，優點）**：計數零漂移；Hard-Gate truth/tier-gate/hook-wiring 全 GREEN（歷史 neutering 未復發）。
- **D5-2（高）**：KB 研究狀態 5 檔 last_verified 凍在 2026-05-18 = 28 天（2× 14 天 gate）→ 每次 UserPromptSubmit 都噴警示。**必須先重新核對內容再 refresh**（直接 bump 日期 = §13 捏造）。
- **D5-3（中）**：snapshot 有人工 banner 但 last_verified 仍 stale → 散文掩蓋、機械 gate 沒清。→ 考慮 auto-derive from CURRENT_FOCUS。
- **D5-4（低）**：ledger 領先 CURRENT_FOCUS 1 天（append-lag，非 over-claim，但是 canary）。
- **D5-5（中）**：1 個 orphan memory（`feedback_html_for_conclusions_md_for_final_check.md`）未進 MEMORY.md 索引 → 不會被 recall；疑與 `feedback_review_format_html_ask_md` 近重複 → 走 /memory-consolidation 判合併 or 索引。

### D6 高重複 → 可複利自動化 — DONE
- **D6-1（高）**：archive 流程（grep-check→safe-ls mtime→大BAM檢→manifest row→mv）全手做、已 2 次、≥3 批待做 → `scripts/infra/archive_candidate.sh` 自動產 manifest row + 印 mv（不自動跑，保 Hard-Gate 刪除確認）。
- **D6-2（高）**：`provenance_stamp.sh` 存在但 manual-only；連定義此規則的 audit 文件自己都沒蓋（本報告已蓋以正視聽）→ wire 成 PostToolUse Write advisory（盤點/audit/status 類 .md 缺 build_branch → 印可貼 stamp）。
- **D6-3（中）**：`number_provenance.py audit`（§13-C）已建但零部署到任何 conclude/report 收尾 → 在 conclude-research/results-report/html-report-build SKILL 收尾清單加一行。
- **D6-4（中）**：safe-ls 方法散文重述 6 份文件、無 helper → `scripts/infra/safe_ls.sh` 把 §12 規則編碼一次（用 verify_output.sh 既有 `ls -d ... | wc -l` idiom 推廣，不重造）。
- **D6-5（中）**：`fill_report.py`（§13-A 構造層，最強 anti-fabrication）零真實部署 → 挑 1 個資料密集 artifact 轉 template+data.json 當範例 + SKILL 指標。

## 改進優先序（payoff / effort）

| P | 改進 | 維度 | effort | 類別 | 我的 scope？ |
|---|---|---|:---:|---|:---:|
| **P1** | search_scope_guard 補 output-dir hang 防護 | D3-1 | S | 風險↓×常見錯誤 | ✅ |
| **P2** | chr2 並行 session 未提交 at-risk | D2-1 | — | 風險↓ | ❌ 留你/原 session |
| **P3** | reconcile §C shared-state 規則矛盾 | D1-1/2 | S | 風險↓ | ⚠ 需你定調 |
| **P4** | archive_candidate.sh + safe_ls.sh 腳本化 | D6-1/4 | M | 可複利×高重複 | ✅ |
| **P5** | provenance_stamp wire advisory | D6-2 | S | 防 P-17 | ✅ |
| **P6** | KB 5 檔 28d stale 重核+refresh | D5-2/3 | M | 風險↓ | ⚠ 需重核內容 |
| **P7** | output-token >3M advisory wire | D3-3 | S | 風險↓ | ✅ |
| **P8** | fill_report/number_provenance audit 部署範例 | D6-3/5 | M | 可複利 | ✅ |
| **P9** | MEMORY.md orphan 索引/合併 | D5-5 | S | 風險↓ | ✅ |
| **P10** | 新主軸輸出 final/ tree 建骨架 | D4-4 | S | 高效率 | ✅ |

## 維護紀錄
- 本報告數字全來自 2026-06-15 真實指令輸出（git/ls -d/harness_health/hook 源碼）；無 du/find 遞迴 ISM 輸出目錄（§12）。
- 改進 P1-P10 待用戶勾選才執行（「先盤點 → 讓我知道 → 然後改進」）。chr2（P2）絕不由本 session commit。
