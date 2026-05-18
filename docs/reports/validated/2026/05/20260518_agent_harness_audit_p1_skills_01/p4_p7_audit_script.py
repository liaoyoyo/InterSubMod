#!/usr/bin/env python3
"""P4-P7 + Final Consolidated Fix Plan."""
import json
from datetime import datetime
from pathlib import Path

OUTDIR = Path("/tmp/skill_audit_20260518")

# =============================================================================
# P4 Hooks Deep Audit (補 M4 未深入點)
# =============================================================================
P4_HOOKS = {
    "multi_hook_per_event": [
        ("UserPromptSubmit", "—", 4, "並行執行各獲獨立 stdin 副本；hook 之間無 explicit ordering 保證"),
        ("PreToolUse", "Bash", 5, "5 hooks 含 rm-rf guard / pre_commit_compile / kb_schema / pipeline_block / no_binary — 可整合為單一 router_pre_bash.sh"),
        ("PostToolUse", "Edit|Write", 8, "8 hooks 含 cpp_edit / trigger_routing / evidence_ledger / terminology / md_link / kb_sot / standalone / skill_change_audit — 高耦合"),
        ("PostToolUse", "Bash", 3, "compile_success + bash_logger + post_cpp_commit"),
    ],
    "stdin_consumers": 17,  # of 21 total hooks
    "hard_gate_count": 5,
    "silent_failure_count": 6,
    "issues": [
        ("H4-1", "5 個 PreToolUse Bash hooks 整合為單一 router_pre_bash.sh → 降 maintenance cost / 統一 stdin 處理", "P3"),
        ("H4-2", "8 個 PostToolUse Edit|Write hooks 是否全必要？評估 deprecate 候選（terminology_guard / md_link_check 與 kb_sot_guard 是否重疊）", "P3"),
        ("H4-3", "6 silent failure hooks 改 stderr → InterSubMod/docs/postmortems/hook_failures.log", "P0"),
        ("H4-4", "Hook 並行執行假設未驗證 — 應跑 hook_concurrency_test.sh 確認", "P3"),
        ("H4-5", "缺 hook latency telemetry — 22 hooks × 每 turn 可能加 0.5-2s", "P4"),
    ],
}

# =============================================================================
# P5 Agents Audit
# =============================================================================
P5_AGENTS = {
    "active": [
        "architect", "developer", "headless-research", "literature-reviewer",
        "methodology-reviewer", "optimizer", "paper-miner", "parallel-analysis",
        "parallel-benchmark", "release", "researcher", "research-orchestrator",
        "reviewer", "tester",
    ],
    "archived": ["intersubmod-weekly-research-agent"],
    "plugin_overlaps": [
        ("architect", "feature-dev:code-architect", "high",
         "兩者都做 design — InterSubMod architect 更研究導向 (docs/plans/) / feature-dev 更 codebase 導向"),
        ("reviewer", "feature-dev:code-reviewer + pr-review-toolkit:*",
         "high", "InterSubMod reviewer = 科學家驗證數據 / plugin = 程式碼審查 — scope 不同但命名衝突"),
        ("researcher", "Explore (builtin)", "medium",
         "InterSubMod researcher = web search 論文 / Explore = codebase 搜尋 — scope 不同"),
        ("optimizer", "code-simplifier:code-simplifier (plugin)", "medium",
         "InterSubMod optimizer 更廣 (品質+效能+可讀性) / plugin 只簡化"),
        ("developer", "—", "n/a", "無顯著 plugin 對應；InterSubMod 專屬 C++ PDD"),
    ],
    "issues": [
        ("A5-1", "archive/ 只有 1 個 .md + README — 是否該真刪 intersubmod-weekly-research-agent？", "P3"),
        ("A5-2", "reviewer 命名與 plugin code-reviewer 衝突 — 評估改名為 scientist-reviewer", "P4"),
        ("A5-3", "缺 cache_hit / cost / completion telemetry — 不知道哪些 agent 真有用", "P3"),
        ("A5-4", "parallel-benchmark + parallel-analysis 命名一致 ✅", "—"),
        ("A5-5", "research-orchestrator 與 myPPT 風格類似（dispatcher）— 但 scope 不同（meta 路由 vs PPT 場景）", "—"),
    ],
}

# =============================================================================
# P6 Commands / Templates / 5 Entry Points
# =============================================================================
P6_COMPONENTS = {
    "commands": [
        ("/build", "編譯專案 Release", "OK"),
        ("/git-commit", "提交變更含 Co-Author", "OK"),
        ("/git-finish", "合併分支回 develop", "OK"),
        ("/git-start", "建立功能分支", "OK"),
        ("/test-quick", "快速 chr19 驗證 <30s", "OK"),
        ("/test-data", "完整數據測試 ~1min", "OK"),
        ("/test-full", "完整流程 ~5min", "OK"),
        ("/validate", "研究假說 benchmark", "OK"),
    ],
    "templates": [
        ("templates/postmortem.md", "SRE Blameless 8 段", "OK"),
        ("templates/research_index.md", "Pre-reg 3 欄 + G1-G5", "OK"),
    ],
    "entry_points": [
        ("AGENTS.md", "0d", "0d, 12.5KB", "✅"),
        (".claude/CLAUDE.md", "0d", "0d, 7.8KB", "✅"),
        ("docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md", "9d", "9d, 16KB", "⚠ stale (>7d cap)"),
        ("docs/CURRENT_FOCUS.md", "1d", "1d, 35KB (hook 截 3KB)", "✅"),
        ("research/autoresearch/research_direction.md", "9d", "9d, 3.7KB", "⚠ stale (>7d cap)"),
    ],
    "issues": [
        ("E6-1", "20260424_AI啟動壓縮上下文 9d 未更新 — 屬 SoT 文件需週級維護", "P2"),
        ("E6-2", "research_direction.md 9d 未更新 — autoresearch queue 是否仍 active？", "P3"),
        ("E6-3", "8 commands 描述短 (各 <30 char) ✅ 業界推薦長度", "—"),
        ("E6-4", "缺 templates 第 3 個：M15 提問模板 (P0 fix item)", "P0"),
        ("E6-5", "缺 templates 第 4 個：M10 四象限 task triage (P3 fix item)", "P3"),
    ],
}

# =============================================================================
# P7 Documentation Drift
# =============================================================================
P7_DRIFT = {
    "broken_refs_real": [
        ("InterSubMod/.claude/skills/weekly-master-draft/outlines/UPGRADE_PLAN_FOR_WEEKLY_REPORT.md",
         "Referenced in: InterSubMod/.claude/skills/pptx-build/SKILL.md",
         "weekly-master-draft skill 已被 /weekly-report 取代但 pptx-build 仍引用舊路徑"),
    ],
    "placeholder_refs": [
        "InterSubMod/docs/concepts/DAG/<topic>.md (placeholder, OK)",
        "InterSubMod/docs/postmortems/<YYYYMMDD>_<topic>.md (template syntax, OK)",
        "InterSubMod/docs/reports/validated/YYYY/MM/... (template syntax, OK)",
    ],
    "issues": [
        ("D7-1", "pptx-build SKILL.md 第 X 行引用 weekly-master-draft → 改為 /weekly-report 對應 master 路徑", "P1"),
        ("D7-2", "其他 placeholder 路徑（<topic>, {YYYY}）保持作為 template syntax — 無 fix", "—"),
        ("D7-3", "缺自動 broken-ref scanner — 評估加 PostToolUse markdown_link_check (已存 md_link_check.sh 看實際效果)", "P3"),
    ],
}

# =============================================================================
# CONSOLIDATED FIX PLAN
# =============================================================================
CONSOLIDATED = {
    "P0": {  # Critical (1-2 hr)
        "label": "🔴 P0 Critical (1.5 hr)",
        "color": "#cf222e",
        "items": [
            ("M15", "新建 InterSubMod/templates/user_query_template.md 6 欄提問模板", "30 min"),
            ("H4-3", "改 6 silent failure hooks → 失敗寫 InterSubMod/docs/postmortems/hook_failures.log", "45 min"),
            ("P1-YAML", "修 5 YAML invalid skill description (feature-layered-observation / problem-framing-ideation / provenance-tier-audit / research-loop / verification-loop)", "30 min"),
        ],
    },
    "P1": {  # This week (5 hr)
        "label": "🟠 P1 This Week (5 hr)",
        "color": "#fb8500",
        "items": [
            ("P1-D2", "13 個非工具類 skill 補 D2 診斷段 (Phase chain / Dependencies / Failure Mode)", "3.5 hr"),
            ("P1-D3", "5 個非工具類 skill 加 /scientific-rigor cross-ref", "30 min"),
            ("M2", "/problem-framing-ideation §1 加 How Much (5W2H 升級)", "30 min"),
            ("M3", "templates/research_index.md Pre-reg 加 SMART 5-字 checklist", "30 min"),
            ("D7-1", "pptx-build 引用改 weekly-master-draft → /weekly-report 路徑", "5 min"),
        ],
    },
    "P2": {  # Next 2 weeks (3 hr)
        "label": "🟡 P2 Next 2 weeks (3 hr)",
        "color": "#bf8700",
        "items": [
            ("M3-Permission-Audit", "allow_list_audit.sh 偵測 158 entries 未用 / 過度具體", "30 min"),
            ("M7-Cache-Telemetry", "cache_telemetry.sh — Stop hook 解析 cost report 抽 cache_read / cache_creation tokens", "1 hr"),
            ("M4-SWOT", "/methodology-audit Step 3 OPTIONS 加 SWOT 2x2 matrix template", "45 min"),
            ("M5-TLDR", "AGENTS.md §15 加「Tier 3+ 首句 TL;DR」規則", "15 min"),
            ("M13-Ledger", "evidence_ledger.jsonl schema 加 next_action + identified_issues 欄", "1 hr"),
            ("E6-1", "更新 20260424_AI啟動壓縮上下文與研究索引_01.md (9d stale)", "30 min"),
        ],
    },
    "P3": {  # Next month (4 hr)
        "label": "🔵 P3 Next Month (4 hr)",
        "color": "#0969da",
        "items": [
            ("H4-1", "整合 5 個 PreToolUse Bash hooks → 單一 router_pre_bash.sh", "1 hr"),
            ("H4-2", "evaluate 8 個 PostToolUse Edit|Write 重疊 → deprecate 候選 2-3 個", "45 min"),
            ("Compact-Test", "compact_test.sh — /compact 後保留指令自動回歸測試", "45 min"),
            ("Subagent-Logger", "subagent_completion_logger.sh — SubagentStop 記錄 cost/cache/status", "30 min"),
            ("M6-Worktree", "parallel-benchmark agent 預設 isolation: worktree", "15 min"),
            ("Rules-Paths", "4 個 .claude/rules/ 加 paths frontmatter (省 5KB always-loaded)", "30 min"),
            ("M2-Desc-Shorten", "5 個 max desc > 600 chars skill 縮短到 ≤400 chars", "30 min"),
            ("M10-Matrix", "新建 templates/task_priority_matrix.md (Eisenhower)", "20 min"),
        ],
    },
    "P4": {  # Future / opt-in (3 hr)
        "label": "⚫ P4 Future / Opt-in (3 hr)",
        "color": "#6e7781",
        "items": [
            ("Injection-Guard", "external_input_sanitizer.sh — WebFetch / 大 Read 偵測 prompt injection", "1 hr"),
            ("Memory-Archive", "memory archive 機制 (>1 年 stale → memory/archive/)", "30 min"),
            ("Memory-Recall-Log", "memory_recall_logger.sh — 量化 feedback memory 引用率", "45 min"),
            ("Modify-Hook-Demo", "dangerous_path_redirect.sh — PreToolUse modify input 示範", "30 min"),
            ("Cross-Version-Baseline", "Q1-Q13 跨 Claude version baseline JSON (每 major 升版前跑)", "—"),
        ],
    },
}

# =============================================================================
# Render HTML
# =============================================================================
html = f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>P4-P7 Audit + Final Consolidated Fix Plan</title>
<style>
  body {{ font-family: -apple-system, "Segoe UI", "Microsoft JhengHei", "Noto Sans CJK TC", sans-serif;
         max-width: 1500px; margin: 2em auto; padding: 0 1em; line-height: 1.55; color: #1f2328; background: #fafbfc; }}
  h1 {{ border-bottom: 2px solid #d0d7de; padding-bottom: 0.3em; }}
  h2 {{ border-bottom: 1px solid #d0d7de; padding-bottom: 0.2em; margin-top: 2em; }}
  h3 {{ margin-top: 1.5em; }}
  .meta {{ color: #57606a; font-size: 0.9em; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.88em; margin: 1em 0; background: #fff; }}
  th, td {{ border: 1px solid #d0d7de; padding: 6px 10px; text-align: left; vertical-align: top; }}
  th {{ background: #f6f8fa; font-weight: 600; }}
  code {{ background: #f6f8fa; padding: 1px 4px; border-radius: 3px; font-family: ui-monospace,monospace; font-size: 0.9em; }}
  .phase-box {{ background: #fff; border: 1px solid #d0d7de; border-left: 5px solid; border-radius: 6px; padding: 1em 1.5em; margin: 1em 0; }}
  .priority {{ display: inline-block; padding: 2px 8px; border-radius: 3px; color: white; font-weight: bold; font-size: 0.85em; }}
  .status {{ display: inline-block; padding: 2px 6px; border-radius: 3px; font-weight: bold; font-size: 0.85em; }}
  .stat-ok {{ background: #dafbe1; color: #1a7f37; }}
  .stat-warn {{ background: #fff8c5; color: #bf8700; }}
  .stat-bad {{ background: #ffeef0; color: #cf222e; }}
  details {{ background: #fff; border: 1px solid #d0d7de; border-radius: 4px; padding: 0.5em 1em; margin: 0.5em 0; }}
  summary {{ cursor: pointer; font-weight: bold; }}
  .item-id {{ font-family: ui-monospace,monospace; font-size: 0.9em; color: #0969da; font-weight: bold; }}
</style></head>
<body>
<h1>P4-P7 Audit + Final Consolidated Fix Plan</h1>
<p class="meta">Audit date: {datetime.now().strftime('%Y-%m-%d %H:%M')} · This report consolidates P1 (Skills) + P2 (Industry) + P3 (Methodology) + P4 (Hooks deep) + P5 (Agents) + P6 (Commands/Templates/Entry) + P7 (Doc drift) findings into a single execution plan.</p>

<h2>📊 Overall Harness Health</h2>
<table>
<thead><tr><th>Phase</th><th>Scope</th><th>Status</th><th>Critical findings count</th></tr></thead>
<tbody>
<tr><td><strong>P1 Skills</strong></td><td>42 skills × 6 維度</td><td><span class="status stat-warn">⚠ 18 issues</span></td><td>5 YAML invalid + 13 D2 缺診斷 + 5 D3 缺 cross-ref</td></tr>
<tr><td><strong>P2 Industry</strong></td><td>7 modules + 5 challenges</td><td><span class="status stat-warn">⚠ Mixed</span></td><td>M4 28% silent failure / M3 allow list bloat / M7 no telemetry</td></tr>
<tr><td><strong>P3 Methodology</strong></td><td>17 framework 對照</td><td><span class="status stat-ok">✅ 41% strong</span></td><td>1 P0 (M15 提問模板)</td></tr>
<tr><td><strong>P4 Hooks Deep</strong></td><td>22 hooks routing</td><td><span class="status stat-warn">⚠ 5 issues</span></td><td>20 hooks in multi-matcher; stdin handling 17/21</td></tr>
<tr><td><strong>P5 Agents</strong></td><td>14 active + 2 archive</td><td><span class="status stat-ok">✅ OK</span></td><td>1 archive cleanup candidate; plugin naming overlap (low impact)</td></tr>
<tr><td><strong>P6 Commands/Templates/Entry</strong></td><td>8 + 2 + 5</td><td><span class="status stat-ok">✅ OK</span></td><td>2 entry 9d stale (manual + research_direction)</td></tr>
<tr><td><strong>P7 Doc drift</strong></td><td>broken refs scan</td><td><span class="status stat-ok">✅ OK</span></td><td>1 真 broken ref (pptx-build → weekly-master-draft)</td></tr>
</tbody>
</table>

<h2>P4 — Hooks Deep Dive (補 M4 未深入)</h2>
<h3>Multi-hook Per Event 分析</h3>
<table>
<thead><tr><th>Event</th><th>Matcher</th><th>Hooks #</th><th>備註</th></tr></thead>
<tbody>
"""
for evt, mt, cnt, note in P4_HOOKS["multi_hook_per_event"]:
    html += f"<tr><td>{evt}</td><td><code>{mt}</code></td><td>{cnt}</td><td>{note}</td></tr>\n"
html += f"""</tbody></table>

<p><strong>Stdin consumers</strong>: {P4_HOOKS['stdin_consumers']}/21 hooks 讀 stdin (jq -r / cat)</p>
<p><strong>Hard gate (exit 2)</strong>: {P4_HOOKS['hard_gate_count']} hooks ✅</p>
<p><strong>Silent failure</strong>: {P4_HOOKS['silent_failure_count']}/21 (28%) <span class="status stat-bad">需修</span></p>

<h3>P4 Issues</h3>
<table>
<thead><tr><th>ID</th><th>Issue</th><th>Priority</th></tr></thead>
<tbody>
"""
for iid, desc, pri in P4_HOOKS["issues"]:
    pc = {"P0":"#cf222e","P1":"#fb8500","P2":"#bf8700","P3":"#0969da","P4":"#6e7781"}.get(pri,"#6e7781")
    html += f'<tr><td><span class="item-id">{iid}</span></td><td>{desc}</td><td><span class="priority" style="background:{pc}">{pri}</span></td></tr>\n'
html += "</tbody></table>\n"

html += f"""
<h2>P5 — Agents Audit (14 active + 1 archived)</h2>
<h3>Active agents</h3>
<p>{', '.join('<code>'+a+'</code>' for a in P5_AGENTS['active'])}</p>

<h3>Plugin overlap 分析</h3>
<table>
<thead><tr><th>InterSubMod agent</th><th>Plugin equivalent</th><th>重疊度</th><th>備註</th></tr></thead>
<tbody>
"""
for ism, plg, deg, note in P5_AGENTS["plugin_overlaps"]:
    badge = {"high":"#fb8500","medium":"#bf8700","n/a":"#6e7781"}.get(deg, "#6e7781")
    html += f'<tr><td><code>{ism}</code></td><td><code>{plg}</code></td><td><span class="status" style="background:{badge};color:white">{deg}</span></td><td>{note}</td></tr>\n'
html += "</tbody></table>\n"

html += "<h3>P5 Issues</h3>\n<table><thead><tr><th>ID</th><th>Issue</th><th>Priority</th></tr></thead><tbody>\n"
for iid, desc, pri in P5_AGENTS["issues"]:
    if pri == "—": continue
    pc = {"P0":"#cf222e","P1":"#fb8500","P2":"#bf8700","P3":"#0969da","P4":"#6e7781"}.get(pri,"#6e7781")
    html += f'<tr><td><span class="item-id">{iid}</span></td><td>{desc}</td><td><span class="priority" style="background:{pc}">{pri}</span></td></tr>\n'
html += "</tbody></table>\n"

html += "<h2>P6 — Commands / Templates / 5 Entry Points</h2>\n"
html += "<h3>5 Entry Points (CLAUDE.md §9)</h3>\n"
html += "<table><thead><tr><th>Path</th><th>Status</th></tr></thead><tbody>\n"
for path, age, info, st in P6_COMPONENTS["entry_points"]:
    badge = "stat-ok" if "✅" in st else "stat-warn"
    html += f'<tr><td><code>InterSubMod/{path}</code></td><td><span class="status {badge}">{st}</span> {info}</td></tr>\n'
html += "</tbody></table>\n"

html += "<h3>P6 Issues</h3>\n<table><thead><tr><th>ID</th><th>Issue</th><th>Priority</th></tr></thead><tbody>\n"
for iid, desc, pri in P6_COMPONENTS["issues"]:
    if pri == "—": continue
    pc = {"P0":"#cf222e","P1":"#fb8500","P2":"#bf8700","P3":"#0969da","P4":"#6e7781"}.get(pri,"#6e7781")
    html += f'<tr><td><span class="item-id">{iid}</span></td><td>{desc}</td><td><span class="priority" style="background:{pc}">{pri}</span></td></tr>\n'
html += "</tbody></table>\n"

html += "<h2>P7 — Documentation Drift</h2>\n"
html += "<h3>真 Broken References (排除 placeholder)</h3>\n"
html += "<table><thead><tr><th>Broken target</th><th>Source</th><th>Reason</th></tr></thead><tbody>\n"
for tgt, src, reason in P7_DRIFT["broken_refs_real"]:
    html += f'<tr><td><code>{tgt}</code></td><td>{src}</td><td>{reason}</td></tr>\n'
html += "</tbody></table>\n"
html += f"<p>Placeholder refs (NOT broken — template syntax): {len(P7_DRIFT['placeholder_refs'])} entries</p>\n"

# =============================================================================
# Final Consolidated Fix Plan
# =============================================================================
html += """<h2>🎯 Final Consolidated Fix Plan</h2>
<p>合併 P1-P7 所有 fix items，依優先級分 5 階段。每階段含完整實作路徑 + 工時。</p>
"""
total_hours = 0
for phase, info in CONSOLIDATED.items():
    phase_hours = 0
    html += f"""<div class="phase-box" style="border-left-color:{info['color']}">
<h3 style="margin:0 0 0.5em 0;color:{info['color']}">{info['label']}</h3>
<table><thead><tr><th>Item ID</th><th>Action</th><th>Est. Time</th></tr></thead><tbody>
"""
    for iid, action, hrs in info["items"]:
        html += f'<tr><td><span class="item-id">{iid}</span></td><td>{action}</td><td>{hrs}</td></tr>\n'
        # Parse hours from string for total
        try:
            if "hr" in hrs:
                phase_hours += float(hrs.split()[0])
            elif "min" in hrs:
                phase_hours += float(hrs.split()[0]) / 60
        except Exception:
            pass
    html += f"</tbody></table>\n<p><strong>Phase total: ~{phase_hours:.1f} hr</strong></p>\n</div>\n"
    total_hours += phase_hours

html += f"""<h2>📈 Total Effort Estimate</h2>
<table>
<thead><tr><th>Phase</th><th>Hours</th><th>Cumulative</th><th>Recommended Window</th></tr></thead>
<tbody>
<tr><td>P0 Critical</td><td>1.5</td><td>1.5</td><td>今天 / 立刻</td></tr>
<tr><td>P1 This week</td><td>5.0</td><td>6.5</td><td>本週剩餘工作日</td></tr>
<tr><td>P2 Next 2 weeks</td><td>3.0</td><td>9.5</td><td>W4-W5</td></tr>
<tr><td>P3 Next month</td><td>4.0</td><td>13.5</td><td>W6-W7</td></tr>
<tr><td>P4 Future/opt-in</td><td>3.0</td><td>16.5</td><td>W8+ / 評估後</td></tr>
<tr><td><strong>Total</strong></td><td><strong>{total_hours:.1f} hr</strong></td><td>—</td><td>~4 週可全清</td></tr>
</tbody>
</table>

<h2>⚖️ Dependency Map（為什麼依此順序）</h2>
<table>
<thead><tr><th>Dependency</th><th>順序理由</th></tr></thead>
<tbody>
<tr><td>P0 YAML fix → P1 D2 diagnostic 段補</td><td>D2 補時若 YAML 無效，YAML parser 仍寬鬆 trigger，但補完 YAML 後 D2 diff 更可靠</td></tr>
<tr><td>P0 hook silent failure 修 → P2 cache_telemetry</td><td>cache telemetry 屬新 hook，先解 silent failure 模式才能保新 hook 正確 fail-loud</td></tr>
<tr><td>P1 D2 補完 → P3 H4-1 router 整合</td><td>整合 5 PreToolUse hooks 前需要每個 hook 有 Failure Mode 段（D2）才能評估 deprecate 候選</td></tr>
<tr><td>P2 ledger 7 欄擴充 → P4 memory archive</td><td>archive 需依新 schema 標記 stale entries</td></tr>
<tr><td>M15 templates 模板（P0）→ M2/M3/M10 其他模板（P1/P3）</td><td>建立 InterSubMod/templates/ 統一風格規範後其他 template 才一致</td></tr>
</tbody>
</table>

<h2>🏗️ Cross-cutting Architecture Improvements (NOT in fix list)</h2>
<details>
<summary>未列入 P0-P4 但值得評估的 architecture-level 改進</summary>
<ul>
<li><strong>Skill ↔ Hook ↔ Agent 三層責任邊界文件化</strong>: CLAUDE.md §X 加一張表釐清 (skill = behavior recipe / hook = guardrail / agent = autonomous worker)</li>
<li><strong>5 入口聲明降到 3</strong>: CURRENT_FOCUS 已被 SessionStart hook 注入，是否仍需 AGENTS.md §9 列為「主動 Read」？評估合併</li>
<li><strong>Permission mode 顯式化</strong>: 在 settings.local.json 加 <code>"permissionMode": "default"</code> 等 (CLI default 但顯式化便於 audit)</li>
<li><strong>Hook chain visualization</strong>: 22 hooks 跨 5 events 用 mermaid graph 自動產生（git-pre-push）</li>
<li><strong>Memory feedback hit-rate 量化</strong>: 加 grep counter — 哪些 feedback memory 從未被引用？候選 archive</li>
</ul>
</details>

<h2>📦 Deliverables 索引</h2>
<table>
<thead><tr><th>檔案</th><th>內容</th></tr></thead>
<tbody>
<tr><td><code>InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/audit_report.standalone.html</code></td><td>P1 Skills audit (42 × 6)</td></tr>
<tr><td><code>InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/p2_industry_audit.standalone.html</code></td><td>P2 業界對照 (7 modules + 5 challenges)</td></tr>
<tr><td><code>InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/p3_methodology_alignment.standalone.html</code></td><td>P3 17 方法論對照</td></tr>
<tr><td><code>InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/p4_to_p7_final.standalone.html</code></td><td>🎯 本報告 (P4-P7 + Final Fix Plan)</td></tr>
</tbody></table>

<hr>
<p class="meta">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} · Total audit effort: ~3 hr spent / 16.5 hr fix backlog</p>
</body></html>
"""
(OUTDIR / "p4_to_p7_final.standalone.html").write_text(html, encoding="utf-8")

# Save consolidated JSON
final_json = {
    "audit_date": datetime.now().strftime("%Y-%m-%d %H:%M"),
    "p4_hooks": P4_HOOKS,
    "p5_agents": P5_AGENTS,
    "p6_components": P6_COMPONENTS,
    "p7_drift": P7_DRIFT,
    "consolidated_fix_plan": CONSOLIDATED,
    "total_hours_estimate": total_hours,
}
(OUTDIR / "p4_to_p7_final.json").write_text(json.dumps(final_json, ensure_ascii=False, indent=2))
print(f"HTML: {OUTDIR}/p4_to_p7_final.standalone.html ({len(html)} chars)")
print(f"JSON: {OUTDIR}/p4_to_p7_final.json")
print(f"\nTotal fix effort: {total_hours:.1f} hr")
print(f"P0 items: {len(CONSOLIDATED['P0']['items'])} ({sum(1 for _,_,h in CONSOLIDATED['P0']['items'])} fixes)")
print(f"P1 items: {len(CONSOLIDATED['P1']['items'])}")
print(f"P2 items: {len(CONSOLIDATED['P2']['items'])}")
print(f"P3 items: {len(CONSOLIDATED['P3']['items'])}")
print(f"P4 items: {len(CONSOLIDATED['P4']['items'])}")
