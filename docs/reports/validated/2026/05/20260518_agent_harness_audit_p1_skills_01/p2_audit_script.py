#!/usr/bin/env python3
"""P2 Industry Best Practice Alignment — InterSubMod vs 7 Harness Engineering modules."""
import json
import re
from pathlib import Path
from datetime import datetime

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
OUTDIR = Path("/tmp/skill_audit_20260518")

# Reuse P1 audit JSON
p1 = json.loads((OUTDIR / "audit_results.json").read_text(encoding="utf-8"))

# Collect M3 Permission stats
settings = json.loads((ROOT / ".claude/settings.local.json").read_text())
allow_list = settings.get("permissions", {}).get("allow", [])
deny_list = settings.get("permissions", {}).get("deny", [])
hooks_cfg = settings.get("hooks", {})

# Collect M4 Hook stats
hook_scripts = sorted((ROOT / "scripts/hooks").glob("*.sh"))
hook_metrics = {
    "total_scripts": len(hook_scripts),
    "silent_failure_count": 0,
    "hard_gate_count": 0,
    "stdin_consume_count": 0,
    "scripts": [],
}
for hs in hook_scripts:
    text = hs.read_text()
    has_silent = "2>/dev/null || true" in text or "|| exit 0" in text
    has_hard = "exit 2" in text
    has_stdin = "tool_input" in text or 'jq -r' in text
    hook_metrics["silent_failure_count"] += int(has_silent)
    hook_metrics["hard_gate_count"] += int(has_hard)
    hook_metrics["stdin_consume_count"] += int(has_stdin)
    hook_metrics["scripts"].append({
        "name": hs.name,
        "silent_failure": has_silent,
        "hard_gate": has_hard,
        "consumes_stdin": has_stdin,
        "line_count": len(text.splitlines()),
    })

# Collect M5 Memory stats
mem_dir = Path("/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory")
feedback_files = sorted(mem_dir.glob("feedback_*.md"))
project_files = sorted(mem_dir.glob("project_*.md"))
reference_files = sorted(mem_dir.glob("reference_*.md"))
memory_metrics = {
    "feedback_count": len(feedback_files),
    "project_count": len(project_files),
    "reference_count": len(reference_files),
    "MEMORY_md_lines": len((mem_dir / "MEMORY.md").read_text().splitlines()) if (mem_dir / "MEMORY.md").exists() else 0,
}

# Collect M6 Agent stats
agent_dir = ROOT / ".claude/agents"
active_agents = [f.name for f in agent_dir.glob("*.md") if f.is_file()]
archived_agents = [f.name for f in (agent_dir / "archive").glob("*.md") if f.is_file()] if (agent_dir / "archive").exists() else []

# Check if any agent uses worktree default in frontmatter
worktree_default = 0
for af in (agent_dir / a for a in active_agents):
    try:
        if "isolation: worktree" in af.read_text() or "isolation: \"worktree\"" in af.read_text():
            worktree_default += 1
    except Exception:
        pass

# Skill metrics from P1
skill_metrics = {
    "total": p1["summary"]["total_skills"],
    "user_invocable_real": 33,  # corrected manually (PyYAML bool bug in P1)
    "paths_scoped": p1["summary"]["paths_scoped_count"],
    "yaml_invalid": p1["summary"]["yaml_invalid_count"],
    "avg_desc_len": p1["summary"]["avg_desc_len"],
    "max_desc_len": max(r["desc_len"] for r in p1["results"]),
}

# CLAUDE.md / AGENTS.md / rules sizes
static_sizes = {
    "CLAUDE.md": (ROOT / ".claude/CLAUDE.md").stat().st_size,
    "AGENTS.md": (ROOT / "AGENTS.md").stat().st_size if (ROOT / "AGENTS.md").exists() else 0,
    "rules/total": sum(f.stat().st_size for f in (ROOT / ".claude/rules").glob("*.md")) if (ROOT / ".claude/rules").exists() else 0,
    "CURRENT_FOCUS.md_full": (ROOT / "docs/CURRENT_FOCUS.md").stat().st_size,
    "CURRENT_FOCUS_hook_inject_target": 3000,
}

# 7 modules verdict
verdicts = {
    "M1_PromptCaching": {
        "score": "🟢 Mostly aligned",
        "wins": [
            f"靜態 prefix 25KB（CLAUDE.md {static_sizes['CLAUDE.md']//1024}KB + AGENTS.md {static_sizes['AGENTS.md']//1024}KB + rules {static_sizes['rules/total']//1024}KB）= cache-friendly",
            f"Memory MEMORY.md 116 行 (~10KB) always-loaded prefix",
            "CURRENT_FOCUS 完整 35KB → hook 截 3KB 注入（避免破 cache）",
        ],
        "gaps": [
            "**無 cache hit telemetry** — 不知道實際命中率，無法驗證 cache 策略效果",
            "**rules/ 4 個 always-loaded（5KB）** — CLAUDE.md §5 已標規劃加 paths frontmatter 變條件式，未執行",
            "Hook 注入位置在 UserPromptSubmit（每次都加 = 每 prompt 後 invalidate 部分 cache）— 改 SessionStart-only 或加 cache_control 可省",
        ],
        "fixes": [
            "P3 加 cache_hit_telemetry hook（記錄 turn 後 cache_read_input_tokens / cache_creation_input_tokens）",
            "P4 補 rules/ paths frontmatter（4 個 rules → 條件式，省 5KB always-loaded）",
            "P5 評估改 SessionStart-only inject CURRENT_FOCUS（不 UserPromptSubmit 重複）",
        ],
    },
    "M2_ToolDesign": {
        "score": "🟡 Mixed",
        "wins": [
            "命名統一（slash + kebab-case）",
            "user-invocable 33/42 (78%)",
            "USE WHEN / SKIP WHEN trigger 詞清楚（D4 ✅42/42）",
        ],
        "gaps": [
            f"**avg desc len {skill_metrics['avg_desc_len']} chars / max {skill_metrics['max_desc_len']}** — 業界推薦 150-300 chars 純 trigger，偏長稀釋訊號",
            f"**paths-scoped 僅 {skill_metrics['paths_scoped']}/{skill_metrics['total']}** — 大多 skill 描述 always-discoverable，違反 Deferred Loading 精神",
            f"**5 YAML invalid description** — 跨 parser 行為不一致 latent bug（P1 critical fix）",
            "skill description 含大量中英混 + markdown ** + 引號 — LLM-friendly 但 token 浪費",
        ],
        "fixes": [
            "P0 修 5 YAML invalid（30 min）",
            "P3 評估 max-1411-char desc 縮短到 ≤400 chars（重 trigger 詞、移移 markdown）",
            "P4 補 paths frontmatter 至 ≥20 個 file-scoped skills",
        ],
    },
    "M3_Permission": {
        "score": "🟡 Mixed (allow list bloat)",
        "wins": [
            f"deny list {len(deny_list)} entries（顯式 Deny-First）",
            "Docker sandbox enabled",
            f"Hooks 縱深防禦（hard gate {hook_metrics['hard_gate_count']} 個 exit 2）",
        ],
        "gaps": [
            f"**allow list {len(allow_list)} entries** — 累積無 audit，含許多一次性命令（如 specific bash invocations）",
            "未顯式設 default permission mode（CLI default = default）",
            "缺季度 allow list review 機制",
            "缺 PreToolUse modify input 範例（業界推薦的「透明攔截」未使用）",
        ],
        "fixes": [
            "P3 寫 allow_list_audit.sh — 偵測 > 30 天未用 entry + 過度具體 entry",
            "P4 評估明確設 permissionMode（plan / default / acceptEdits 選一）",
            "P5 設計一個 PreToolUse modify hook 示範（如 production DB → test DB redirect）",
        ],
    },
    "M4_HookSystem": {
        "score": "🟠 Improvement needed (silent failure anti-pattern)",
        "wins": [
            f"{hook_metrics['total_scripts']} hook scripts 覆蓋 5 events",
            f"{hook_metrics['hard_gate_count']} 個 hard gate (exit 2) ✅",
            f"{hook_metrics['stdin_consume_count']} 個 hooks 處理 stdin tool_input",
            "新加 SessionStart + skill_change_audit ✅",
        ],
        "gaps": [
            f"**{hook_metrics['silent_failure_count']}/{hook_metrics['total_scripts']} ({hook_metrics['silent_failure_count']*100//hook_metrics['total_scripts']}%) 用 `2>/dev/null || true`** — 業界明確的 anti-pattern",
            "無統一 hook log 目的地（debug 時很難 trace）",
            "缺 PreToolUse modify input 示範 — Anthropic 官方提的「透明攔截」未實作",
            "PreToolUse Bash matcher 多次（5 個）— 可整合為單一 router hook",
            "缺 hook failure alert 機制（失敗就 silent）",
        ],
        "fixes": [
            "P0 改 6 silent failure hooks → 失敗寫 `docs/postmortems/hook_failures.log`（不 silent）",
            "P3 整合 5 個 PreToolUse Bash hooks → 單一 router_pre_bash.sh",
            "P4 設計示範 modify hook（如 dangerous_path_redirect.sh）",
            "P5 建 hook_health_check.sh weekly 跑",
        ],
    },
    "M5_Memory": {
        "score": "🟢 Best in class",
        "wins": [
            f"{memory_metrics['feedback_count']} feedback + {memory_metrics['project_count']} project + {memory_metrics['reference_count']} reference memory（業界級數量）",
            f"MEMORY.md {memory_metrics['MEMORY_md_lines']} 行索引（< 200 line cap）",
            "Concluded 區明確分隔（防 re-investigation）— 創新做法",
            "Auto memory 系統（4 types user/feedback/project/reference）符合 Cookbook 推薦",
            "今天新加 Cynefin + PF reopen 2 條 — 主動 evolve",
        ],
        "gaps": [
            "**無 Memory Tool（工具呼叫式）** — Anthropic 尚未釋出獨立 Memory Tool，N/A",
            "Concluded 區條目達 30 行（接近 200 line cap 30%） — 未來需 archive 機制",
            "缺 memory recall 命中監測（哪些 memory 真被引用了？）",
        ],
        "fixes": [
            "P4 設計 memory_recall_logger.sh — 偵測 AI response 引用 memory feedback file 並計數",
            "P5 評估 memory archive 機制（>1 年 stale conclusions → memory/archive/ 子目錄）",
        ],
    },
    "M6_SubAgents": {
        "score": "🟡 Functional but no telemetry",
        "wins": [
            f"{len(active_agents)} active agents + {len(archived_agents)} archived（合理分類）",
            "specialist agents (parallel-benchmark / methodology-reviewer / headless-research) 覆蓋場景",
            "CLAUDE.md §8 明確「預設不 spawn」+ 明示觸發語 — 符合 Opus 4.7 literal",
        ],
        "gaps": [
            f"**worktree isolation 預設 OFF**（active agents 中 {worktree_default} 個明示 isolation: worktree） — 業界 Cursor 2.0 / Devin 推薦預設 ON for parallel safety",
            "缺 sub-agent cache hit telemetry — 平行 agent 是否真的省 cache？無資料",
            "缺 sub-agent 失敗回報協議標準化（agent 跑掛了主 agent 不知道）",
        ],
        "fixes": [
            "P4 parallel-benchmark agent 加 default isolation: worktree（安全並行）",
            "P5 加 subagent_completion_logger.sh — SubagentStop hook 記錄 cost / cache / status",
        ],
    },
    "M7_PromptCache": {
        "score": "🟠 Unknown effectiveness (no telemetry)",
        "wins": [
            "CLAUDE.md + AGENTS.md + rules 結構符合 cache prefix 模式（靜態前端）",
            "Skill description progressive disclosure（按需載入）",
        ],
        "gaps": [
            "**無 cache_read_input_tokens / cache_creation_input_tokens telemetry** — 無法量化 90% 降本是否實際達成",
            "CLI 自動處理 cache_control — 無人為控制空間（trade-off）",
            "Sub-agent 獨立 cache 但無命中率資料",
            "UserPromptSubmit hooks 4 個都加 dynamic context — 每次 prompt 後 cache 部分失效",
        ],
        "fixes": [
            "P3 加 cache_telemetry.sh — 在 Stop hook 解析 cost report 抽 cache stats",
            "P5 評估改 SessionStart-only inject（不 UserPromptSubmit 重複）省 cache invalidation",
        ],
    },
}

challenges = {
    "C1_PromptInjection": {
        "status": "🟠 Partial",
        "current": "PreToolUse Bash hooks 攔截 rm -rf / git push --force / pipeline_block_check",
        "gap": "缺 external data (WebFetch / Read of foreign files) PreToolUse sanitization",
        "fix": "P5 加 external_input_sanitizer.sh — WebFetch / Read 大檔時偵測「Ignore previous instructions」pattern",
    },
    "C2_ContextDrift": {
        "status": "🟢 Good",
        "current": "CLAUDE.md §7 /compact 保留指令明確 + auto memory 系統補長期記憶",
        "gap": "無自動測試 compact 後保留度",
        "fix": "P5 加 compact_test.sh — 跑 /compact 後檢核保留指令是否仍在",
    },
    "C3_Evaluation": {
        "status": "🟢 Best in class",
        "current": "verification_guide.md 13 query × 6 維度 fresh CLI 套組 + audit_report.html 機械掃描",
        "gap": "缺跨版本 regression baseline（Q1-Q13 在 Opus 4.7 vs 4.8 diff）",
        "fix": "P6 每 Claude major version 升級時跑 13 query，存 baseline JSON",
    },
    "C4_CostVsLatency": {
        "status": "🟠 No dashboard",
        "current": "skill_validation 跑時記錄 cost ($1.21 Q9+Q10) + sub-agent isolated cache",
        "gap": "缺日常 cost dashboard / 缺 turn-by-turn 比較",
        "fix": "P5 加 weekly_cost_summary.sh — 解析 ~/.claude transcripts 統計",
    },
    "C5_Reproducibility": {
        "status": "🟢 Best in class",
        "current": "evidence_ledger.jsonl + commit hash binding + ai_sessions 歸檔 + manifest.yaml",
        "gap": "skill_validation rerun 偶爾 cost variance ±20%（cache hit/miss）",
        "fix": "P6 加 noise floor 文件化（接受 ±20% variance 為 normal）",
    },
}

# Save analysis JSON
out_data = {
    "metadata": {
        "audit_date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "reference": "Lessons from Building Claude Code: Prompt Caching is Everything + Anthropic Cookbook",
        "scope": "7 Harness Engineering modules + 5 real-world challenges",
    },
    "metrics": {
        "static_prefix_sizes": static_sizes,
        "skill_metrics": skill_metrics,
        "permission_metrics": {"allow_count": len(allow_list), "deny_count": len(deny_list)},
        "hook_metrics": hook_metrics,
        "memory_metrics": memory_metrics,
        "agent_metrics": {"active": len(active_agents), "archived": len(archived_agents), "worktree_default": worktree_default},
    },
    "verdicts": verdicts,
    "challenges": challenges,
}
(OUTDIR / "p2_industry_audit.json").write_text(json.dumps(out_data, ensure_ascii=False, indent=2))
print(f"JSON: {OUTDIR}/p2_industry_audit.json")

# Render HTML
def score_color(score):
    if "🟢" in score: return "#1a7f37"
    if "🟡" in score: return "#bf8700"
    if "🟠" in score: return "#fb8500"
    if "🔴" in score: return "#cf222e"
    return "#6e7781"

html = f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>P2 Industry Best Practice Audit — InterSubMod</title>
<style>
  body {{ font-family: -apple-system, "Segoe UI", "Microsoft JhengHei", "Noto Sans CJK TC", sans-serif;
         max-width: 1400px; margin: 2em auto; padding: 0 1em; line-height: 1.55; color: #1f2328; background: #fafbfc; }}
  h1 {{ border-bottom: 2px solid #d0d7de; padding-bottom: 0.3em; }}
  h2 {{ border-bottom: 1px solid #d0d7de; padding-bottom: 0.2em; margin-top: 2em; }}
  h3 {{ margin-top: 1.5em; }}
  .meta {{ color: #57606a; font-size: 0.9em; }}
  .module {{ background: #fff; border: 1px solid #d0d7de; border-radius: 6px; padding: 1em 1.5em; margin: 1.5em 0; }}
  .module-header {{ display: flex; justify-content: space-between; align-items: center; }}
  .module-score {{ font-weight: bold; padding: 0.3em 0.8em; border-radius: 4px; color: white; }}
  .wins {{ background: #dafbe1; padding: 0.5em 1em; border-radius: 4px; margin: 0.5em 0; }}
  .gaps {{ background: #fff8c5; padding: 0.5em 1em; border-radius: 4px; margin: 0.5em 0; }}
  .fixes {{ background: #ddf4ff; padding: 0.5em 1em; border-radius: 4px; margin: 0.5em 0; }}
  ul {{ margin: 0.3em 0; padding-left: 1.5em; }}
  li {{ margin: 0.3em 0; }}
  code {{ background: #f6f8fa; padding: 1px 4px; border-radius: 3px; font-family: ui-monospace,monospace; font-size: 0.9em; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.9em; margin: 1em 0; background: #fff; }}
  th, td {{ border: 1px solid #d0d7de; padding: 6px 10px; text-align: left; }}
  th {{ background: #f6f8fa; font-weight: 600; }}
  .metric-grid {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 0.5em; margin: 1em 0; }}
  .metric-cell {{ background: #fff; border: 1px solid #d0d7de; padding: 0.6em; border-radius: 4px; }}
  .metric-cell .label {{ font-size: 0.85em; color: #57606a; }}
  .metric-cell .value {{ font-size: 1.2em; font-weight: bold; margin-top: 0.2em; }}
  .priority-fix {{ background: #fff; border-left: 4px solid; padding: 0.8em; margin: 0.5em 0; }}
  .p0 {{ border-color: #cf222e; }}
  .p1 {{ border-color: #fb8500; }}
  .p2 {{ border-color: #bf8700; }}
  .p3 {{ border-color: #0969da; }}
  .p4 {{ border-color: #6e7781; }}
  .badge {{ display: inline-block; padding: 1px 6px; border-radius: 3px; font-size: 0.85em; color: white; }}
</style></head>
<body>
<h1>P2 Industry Best Practice Alignment Audit</h1>
<p class="meta">Audit date: {out_data['metadata']['audit_date']} · Reference: <em>{out_data['metadata']['reference']}</em></p>

<h2>Executive Summary</h2>
<table>
<thead><tr><th>Module</th><th>Score</th><th>Top Gap</th></tr></thead>
<tbody>
"""
for mid, v in verdicts.items():
    color = score_color(v["score"])
    top_gap = v["gaps"][0][:120] + ("..." if len(v["gaps"][0]) > 120 else "") if v["gaps"] else "—"
    html += f'<tr><td><strong>{mid}</strong></td><td><span class="module-score" style="background:{color}">{v["score"]}</span></td><td>{top_gap}</td></tr>\n'

html += """</tbody></table>

<h2>Quick Metrics</h2>
<div class="metric-grid">
"""
metric_cells = [
    ("CLAUDE.md", f"{static_sizes['CLAUDE.md']//1024} KB"),
    ("AGENTS.md", f"{static_sizes['AGENTS.md']//1024} KB"),
    ("rules/ total", f"{static_sizes['rules/total']//1024} KB always-loaded"),
    ("Skills total", f"{skill_metrics['total']}"),
    ("avg skill desc", f"{skill_metrics['avg_desc_len']} chars"),
    ("max skill desc", f"{skill_metrics['max_desc_len']} chars"),
    ("allow list", f"{len(allow_list)} entries"),
    ("deny list", f"{len(deny_list)} entries"),
    ("Hook silent failure", f"{hook_metrics['silent_failure_count']} / {hook_metrics['total_scripts']} ({hook_metrics['silent_failure_count']*100//hook_metrics['total_scripts']}%)"),
    ("Hook hard gate", f"{hook_metrics['hard_gate_count']} (exit 2)"),
    ("Memory feedbacks", f"{memory_metrics['feedback_count']}"),
    ("Active agents", f"{len(active_agents)}"),
]
for label, val in metric_cells:
    html += f'<div class="metric-cell"><div class="label">{label}</div><div class="value">{val}</div></div>\n'
html += "</div>\n"

html += "<h2>7 Modules Detail</h2>\n"
for mid, v in verdicts.items():
    color = score_color(v["score"])
    html += f"""<div class="module">
<div class="module-header">
  <h3 style="margin:0">{mid}</h3>
  <span class="module-score" style="background:{color}">{v["score"]}</span>
</div>
<div class="wins"><strong>✅ Wins:</strong><ul>"""
    for w in v["wins"]: html += f"<li>{w}</li>"
    html += "</ul></div>\n"
    html += '<div class="gaps"><strong>⚠ Gaps:</strong><ul>'
    for g in v["gaps"]: html += f"<li>{g}</li>"
    html += "</ul></div>\n"
    html += '<div class="fixes"><strong>🔧 Fixes:</strong><ul>'
    for f in v["fixes"]: html += f"<li>{f}</li>"
    html += "</ul></div>\n</div>\n"

html += "<h2>Real-world Challenges (5)</h2>\n"
for cid, c in challenges.items():
    color = score_color(c["status"])
    html += f"""<div class="module">
<div class="module-header">
  <h3 style="margin:0">{cid}</h3>
  <span class="module-score" style="background:{color}">{c["status"]}</span>
</div>
<p><strong>當前:</strong> {c["current"]}</p>
<p><strong>Gap:</strong> {c["gap"]}</p>
<p><strong>Fix:</strong> {c["fix"]}</p>
</div>
"""

html += """<h2>Consolidated Fix Priority List (P0-P5)</h2>

<div class="priority-fix p0">
<span class="badge" style="background:#cf222e">P0 — Critical (30 min)</span>
<ul>
<li>修 5 YAML invalid skill description（feature-layered-observation / problem-framing-ideation / provenance-tier-audit / research-loop / verification-loop）</li>
<li>改 6 silent failure hooks → 失敗寫 <code>InterSubMod/docs/postmortems/hook_failures.log</code></li>
</ul>
</div>

<div class="priority-fix p1">
<span class="badge" style="background:#fb8500">P1 — This week (3-5 hr)</span>
<ul>
<li>P1 audit 13 個 D2 缺診斷段補完（每 skill 15 min × 13 = 3.5hr）</li>
<li>P1 audit 5 個 D3 cross-ref backlink 補（30 min）</li>
<li>allow list audit script — 偵測未用 / 過度具體 entry</li>
</ul>
</div>

<div class="priority-fix p2">
<span class="badge" style="background:#bf8700">P2 — Next 2 weeks</span>
<ul>
<li>cache_telemetry.sh — Stop hook 解析 cost report 抽 cache_read / cache_creation tokens（量化 90% 降本）</li>
<li>compact_test.sh — /compact 後保留指令自動回歸測試</li>
<li>subagent_completion_logger.sh — SubagentStop 記錄 cost / cache / status</li>
</ul>
</div>

<div class="priority-fix p3">
<span class="badge" style="background:#0969da">P3 — Next month</span>
<ul>
<li>整合 5 個 PreToolUse Bash hooks → 單一 router_pre_bash.sh（降 hook 維護成本）</li>
<li>parallel-benchmark agent 加 default isolation: worktree（業界推薦）</li>
<li>補 rules/ paths frontmatter（4 個 rules → 條件式，省 5KB always-loaded）</li>
<li>skill description max 1411 chars → ≤400 chars 縮短評估</li>
<li>weekly_cost_summary.sh — 每週 cost dashboard</li>
</ul>
</div>

<div class="priority-fix p4">
<span class="badge" style="background:#6e7781">P4 — Future / opt-in</span>
<ul>
<li>external_input_sanitizer.sh — Prompt Injection 防禦</li>
<li>memory archive 機制（&gt;1 年 stale conclusions）</li>
<li>memory_recall_logger.sh — 量化 memory 引用率</li>
<li>PreToolUse modify input 示範（dangerous_path_redirect.sh）</li>
<li>Q1-Q13 跨 Claude version baseline JSON（每 major 升版前跑）</li>
</ul>
</div>

<h2>後續 audit 階段</h2>
<table>
<thead><tr><th>Phase</th><th>原規劃</th><th>狀態</th></tr></thead>
<tbody>
<tr><td>P1 Skills audit</td><td>42 skills × 6 維度</td><td>✅ 完成 — <code>InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/</code></td></tr>
<tr><td>P2 Industry alignment</td><td>新增（用戶要求）</td><td>✅ 本報告</td></tr>
<tr><td>P3 Hooks audit</td><td>22 hooks 衝突 / silent failure</td><td>部分覆蓋於 M4；剩 hook-order 衝突待查</td></tr>
<tr><td>P4 Agents audit</td><td>16 agents role 清楚度</td><td>部分覆蓋於 M6；剩 plugin agent 重疊待查</td></tr>
<tr><td>P5 Prompt 庫 / commands / templates</td><td>8 commands + 2 templates + 5 入口</td><td>未開始</td></tr>
<tr><td>P6 Documentation drift</td><td>broken refs / stale</td><td>未開始</td></tr>
</tbody></table>

<hr>
<p class="meta">Source: <code>/tmp/p2_industry_audit.py</code> · Data: <code>/tmp/skill_audit_20260518/p2_industry_audit.json</code></p>
</body></html>
"""
(OUTDIR / "p2_industry_audit.html").write_text(html, encoding="utf-8")
print(f"HTML: {OUTDIR}/p2_industry_audit.html ({len(html)} chars)")
