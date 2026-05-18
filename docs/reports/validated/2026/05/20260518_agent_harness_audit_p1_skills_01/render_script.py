#!/usr/bin/env python3
"""Render audit JSON to Markdown table + engineering HTML report."""
import json
from pathlib import Path
from datetime import datetime

OUTDIR = Path("/tmp/skill_audit_20260518")
data = json.loads((OUTDIR / "audit_results.json").read_text(encoding="utf-8"))
results = data["results"]
summary = data["summary"]
yaml_invalid = set(summary["yaml_invalid_skills"])

# Tool-utility skills (D3 cross-ref optional)
TOOL_SKILLS = {"image-gen", "image-vision-check", "html-preview", "html-report-build",
               "pptx-build", "myPPT", "report", "results-report", "weekly-report",
               "structured-tech-report", "infra-ops", "data-audit", "doc-standards",
               "citation-verification", "init-research", "research-dashboard",
               "research-context-loader"}

# ---- Markdown ----
md = []
md.append("# P1 Skills Audit — InterSubMod Agent Harness")
md.append("")
md.append(f"**Audit date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
md.append(f"**Skills scanned**: {summary['total_skills']}")
md.append(f"**Avg line count**: {summary['avg_line_count']} / **Avg desc len**: {summary['avg_desc_len']} chars")
md.append("")
md.append("## Executive Summary")
md.append("")
md.append(f"- 🔴 **CRITICAL — YAML invalid**: {summary['yaml_invalid_count']} skills（description 以 `**` 開頭被 YAML 解讀為 alias reference）")
md.append(f"- 🟠 **D2 內容診斷段缺**: {summary['scores']['d2']['❌']} / {summary['total_skills']} skills 缺 Phase chain / Dependencies / Failure Mode 段")
md.append(f"- 🟡 **D3 cross-ref 缺**: {summary['scores']['d3']['❌']} skills 無與 /scientific-rigor 連結（含 {len([r for r in results if r['d3_score']=='❌' and r['name'] not in TOOL_SKILLS])} 個非工具類真問題）")
md.append(f"- ✅ **D1/D4/D5 全綠**: frontmatter / trigger words / staleness 都 OK")
md.append("")
md.append("## 6 維度評分")
md.append("")
md.append("| Dim | ✅ | ⚠ | ❌ | 維度說明 |")
md.append("|-----|----|----|----|---------|")
for dim, label in [("d1","frontmatter 完整度"), ("d2","內容診斷段"),
                    ("d3","cross-ref to /scientific-rigor"),
                    ("d4","trigger 詞 USE WHEN/SKIP WHEN"),
                    ("d5","staleness mtime")]:
    s = summary["scores"][dim]
    md.append(f"| **{dim.upper()}** | {s['✅']} | {s['⚠']} | {s['❌']} | {label} |")
md.append("")
md.append("## 詳細表 — 42 skills × 6 維度")
md.append("")
md.append("| Skill | Group | D1 frontmatter | D2 診斷 | D3 cross-ref | D4 trigger | D5 mtime | YAML | Lines | Age (d) |")
md.append("|-------|-------|----------------|---------|--------------|------------|----------|------|-------|---------|")

# Sort by group then name
group_order = {"元方法論": 1, "7-Phase Waterfall": 2, "程式開發": 3,
               "文件管理": 4, "報告生成": 5, "研究專用": 6, "未分類": 7}
for r in sorted(results, key=lambda x: (group_order.get(x["group"], 99), x["name"])):
    yv = "🔴" if r["name"] in yaml_invalid else "✅"
    md.append(f"| `{r['name']}` | {r['group']} | {r['d1_score']} | {r['d2_score']} | {r['d3_score']} | {r['d4_score']} | {r['d5_score']} | {yv} | {r['line_count']} | {r['mtime_age_days']} |")
md.append("")
md.append("## 🔴 Critical Fix List (5 YAML invalid)")
md.append("")
md.append("這 5 個 SKILL.md description 以 `**` 開頭，違反 YAML 規範（`**` 視為 alias reference）。Anthropic skill loader 寬鬆度未知，是 latent bug。")
md.append("")
md.append("| Skill | 修法 |")
md.append("|-------|------|")
for sk in summary["yaml_invalid_skills"]:
    md.append(f"| `{sk}` | description 用 `\"...\"` quoted string 包起來，或移除開頭 `**` markdown |")
md.append("")
md.append("**修法範例**：")
md.append("```yaml")
md.append("# ❌ Before")
md.append("description: **P3 PILOT 主要分析 skill** — ...")
md.append("# ✅ After (quoted)")
md.append('description: "**P3 PILOT 主要分析 skill** — ..."')
md.append("# ✅ After (no markdown)")
md.append("description: P3 PILOT 主要分析 skill — ...")
md.append("```")
md.append("")
md.append("## 🟠 D2 內容診斷段缺 (Priority by group)")
md.append("")
md.append("依 user feedback `feedback_skill_md_must_state_dependencies_and_diagnostics`，所有 SKILL.md 應含 3 段：")
md.append("- `Phase & Chain Position`（在工作流哪一階段）")
md.append("- `Dependencies` (Uses / Used-by / Reads / Writes)")
md.append("- `Failure Mode & Diagnostics`（失敗訊號 + 排查路徑）")
md.append("")
md.append("**非工具類 13 skill 缺診斷段（高優先）**：")
md.append("")
non_tool_d2 = [r for r in results if r["d2_score"] == "❌" and r["name"] not in TOOL_SKILLS]
md.append("| Skill | Group | Lines | 缺項 |")
md.append("|-------|-------|-------|------|")
for r in sorted(non_tool_d2, key=lambda x: (group_order.get(x["group"],99), x["name"])):
    missing = []
    if not r["d2_phase"]: missing.append("Phase chain")
    if not r["d2_dependencies"]: missing.append("Deps")
    if not r["d2_failure"]: missing.append("Failure")
    md.append(f"| `{r['name']}` | {r['group']} | {r['line_count']} | {', '.join(missing)} |")
md.append("")
md.append("## 🟡 D3 cross-ref 缺 — 5 個應補（非工具類）")
md.append("")
non_tool_d3 = [r for r in results if r["d3_score"] == "❌" and r["name"] not in TOOL_SKILLS]
md.append("| Skill | Group | 建議 cross-ref hub |")
md.append("|-------|-------|---------------------|")
for r in non_tool_d3:
    hint = {"cpp-change": "scientific-rigor §6 消融原則",
            "inject-hypothesis": "scientific-rigor §7.1 Pre-registration",
            "observation-analysis": "scientific-rigor §5 + auc-confound-guard",
            "results-analysis": "scientific-rigor §3 Effect Size",
            "review-evidence": "scientific-rigor §2 Evidence Tier + §8.4"}
    md.append(f"| `{r['name']}` | {r['group']} | {hint.get(r['name'], 'TBD')} |")
md.append("")
md.append("## Recommendations")
md.append("")
md.append("1. **即刻修 (P0)**: 5 YAML invalid → 30 min 機械式 description quote")
md.append("2. **本週修 (P1)**: 5 個非工具類 D3 cross-ref → 30 min 加 backlink")
md.append("3. **2 週修 (P2)**: 13 個非工具類 D2 診斷段補完 → 每個 15 min × 13 = 3.5hr")
md.append("4. **可選 (P3)**: 12 個工具類 D2/D3 → 工具類本就低 cross-ref 需求，pending 評估")
md.append("")
md.append("## 後續審查階段")
md.append("")
md.append("- P2 Hooks audit")
md.append("- P3 Agents audit")
md.append("- P4 Prompt 庫 / templates audit")
md.append("- P5 Hooks × Skills × Settings 互動 audit")
md.append("- P6 Documentation drift audit")

(OUTDIR / "audit_report.md").write_text("\n".join(md), encoding="utf-8")
print(f"Markdown: {OUTDIR}/audit_report.md ({len(md)} lines)")

# ---- HTML ----
def cell_color(score):
    return {"✅": "#1a7f37", "⚠": "#bf8700", "❌": "#cf222e"}.get(score, "#6e7781")

html = f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>P1 Skills Audit — InterSubMod</title>
<style>
  * {{ box-sizing: border-box; }}
  body {{ font-family: -apple-system, "Segoe UI", "Microsoft JhengHei", "Noto Sans CJK TC", monospace, sans-serif;
         max-width: 1400px; margin: 2em auto; padding: 0 1em; line-height: 1.5; color: #1f2328; background: #fafbfc; }}
  h1 {{ border-bottom: 2px solid #d0d7de; padding-bottom: 0.3em; }}
  h2 {{ border-bottom: 1px solid #d0d7de; padding-bottom: 0.2em; margin-top: 2em; }}
  .meta {{ color: #57606a; font-size: 0.9em; }}
  .summary {{ background: #fff; border: 1px solid #d0d7de; border-radius: 6px; padding: 1em; margin: 1em 0; }}
  .summary-grid {{ display: grid; grid-template-columns: repeat(5, 1fr); gap: 0.5em; margin-top: 0.5em; }}
  .summary-cell {{ background: #f6f8fa; padding: 0.5em; border-radius: 4px; text-align: center; }}
  .summary-cell .num {{ font-size: 1.5em; font-weight: bold; }}
  .green {{ color: #1a7f37; }}
  .yellow {{ color: #bf8700; }}
  .red {{ color: #cf222e; }}
  .gray {{ color: #6e7781; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.85em; margin: 1em 0; background: #fff; }}
  th, td {{ border: 1px solid #d0d7de; padding: 4px 8px; text-align: left; }}
  th {{ background: #f6f8fa; font-weight: 600; position: sticky; top: 0; }}
  tr:hover {{ background: #f6f8fa; }}
  .score {{ text-align: center; font-weight: bold; font-size: 1.1em; }}
  .crit {{ background: #ffeef0; }}
  code {{ background: #f6f8fa; padding: 1px 4px; border-radius: 3px; font-family: ui-monospace,SFMono-Regular,monospace; font-size: 0.9em; }}
  pre {{ background: #f6f8fa; padding: 1em; border-radius: 6px; overflow-x: auto; }}
  .group-元方法論 {{ color: #8250df; }}
  .group-7-Phase {{ color: #0969da; }}
  .group-程式開發 {{ color: #cf222e; }}
  .group-文件管理 {{ color: #1a7f37; }}
  .group-報告生成 {{ color: #bf8700; }}
  .group-研究專用 {{ color: #6639ba; }}
  details {{ margin: 1em 0; }}
  summary {{ cursor: pointer; font-weight: bold; color: #0969da; }}
  .priority-badge {{ display: inline-block; padding: 2px 6px; border-radius: 3px; font-size: 0.8em; font-weight: bold; color: white; }}
  .p0 {{ background: #cf222e; }}
  .p1 {{ background: #bf8700; }}
  .p2 {{ background: #0969da; }}
  .p3 {{ background: #6e7781; }}
</style>
</head>
<body>

<h1>P1 Skills Audit — InterSubMod Agent Harness</h1>
<p class="meta">Audit date: {datetime.now().strftime('%Y-%m-%d %H:%M')} · Skills scanned: <strong>{summary['total_skills']}</strong> · Source: <code>InterSubMod/.claude/skills/*/SKILL.md</code></p>

<div class="summary">
  <h2 style="margin-top:0;border:none">Executive Summary</h2>
  <ul>
    <li>🔴 <strong class="red">CRITICAL — YAML invalid</strong>: {summary['yaml_invalid_count']} skills（description 以 <code>**</code> 開頭被 YAML 解讀為 alias reference）</li>
    <li>🟠 <strong class="yellow">D2 內容診斷段缺</strong>: {summary['scores']['d2']['❌']} / {summary['total_skills']} skills（含 {len([r for r in results if r['d2_score']=='❌' and r['name'] not in TOOL_SKILLS])} 個非工具類真問題）</li>
    <li>🟡 <strong class="yellow">D3 cross-ref 缺</strong>: {summary['scores']['d3']['❌']} skills（含 {len([r for r in results if r['d3_score']=='❌' and r['name'] not in TOOL_SKILLS])} 個非工具類真問題）</li>
    <li>✅ <strong class="green">D1 / D4 / D5 全綠</strong>: frontmatter / trigger words / staleness 都 OK</li>
  </ul>
  <div class="summary-grid">
"""
for dim, label in [("d1","D1 frontmatter"), ("d2","D2 內容診斷"), ("d3","D3 cross-ref"),
                    ("d4","D4 trigger"), ("d5","D5 mtime")]:
    s = summary["scores"][dim]
    html += f"""    <div class="summary-cell">
      <div style="font-size:0.8em; color:#57606a">{label}</div>
      <div class="num"><span class="green">{s['✅']}</span> / <span class="yellow">{s['⚠']}</span> / <span class="red">{s['❌']}</span></div>
    </div>
"""
html += """  </div>
</div>

<h2>🔴 Critical Fix List (P0) — 5 YAML invalid</h2>
<p>這 5 個 SKILL.md description 以 <code>**</code> 開頭，違反 YAML 規範（<code>**</code> 視為 alias reference）。Anthropic skill loader 寬鬆度未知 — 是 <strong>latent bug</strong>，跨 platform / 工具可能行為不一致。</p>
<table>
<thead><tr><th>Skill</th><th>修法</th></tr></thead>
<tbody>
"""
for sk in summary["yaml_invalid_skills"]:
    html += f'<tr class="crit"><td><code>{sk}</code></td><td>description 用 <code>"..."</code> quoted string 包起來，或移除開頭 <code>**</code> markdown</td></tr>\n'
html += """</tbody></table>

<pre><code># ❌ Before
description: **P3 PILOT 主要分析 skill** — ...

# ✅ After (quoted)
description: "**P3 PILOT 主要分析 skill** — ..."

# ✅ After (no markdown)
description: P3 PILOT 主要分析 skill — ...
</code></pre>

<h2>詳細表 — 42 skills × 6 維度</h2>
<table>
<thead>
<tr>
  <th>Skill</th><th>Group</th>
  <th>D1<br>結構</th><th>D2<br>診斷段</th><th>D3<br>cross-ref</th><th>D4<br>trigger</th><th>D5<br>mtime</th>
  <th>YAML</th><th>Lines</th><th>Age (d)</th>
</tr>
</thead>
<tbody>
"""
for r in sorted(results, key=lambda x: (group_order.get(x["group"], 99), x["name"])):
    yv_color = "red" if r["name"] in yaml_invalid else "green"
    yv_sym = "🔴" if r["name"] in yaml_invalid else "✅"
    row_class = ' class="crit"' if r["name"] in yaml_invalid else ''
    html += f'<tr{row_class}>'
    html += f'<td><code>{r["name"]}</code></td>'
    html += f'<td><span class="group-{r["group"].split()[0]}">{r["group"]}</span></td>'
    for dim in ["d1","d2","d3","d4","d5"]:
        score = r[f"{dim}_score"]
        color_class = {"✅":"green","⚠":"yellow","❌":"red"}.get(score, "gray")
        html += f'<td class="score {color_class}">{score}</td>'
    html += f'<td class="score {yv_color}">{yv_sym}</td>'
    html += f'<td>{r["line_count"]}</td><td>{r["mtime_age_days"]}</td>'
    html += "</tr>\n"
html += """</tbody></table>

<h2>🟠 P1 — D2 內容診斷段缺 (13 非工具類)</h2>
<p>依 user feedback <code>feedback_skill_md_must_state_dependencies_and_diagnostics</code>，所有 SKILL.md 應含 3 段：
<strong>Phase &amp; Chain Position</strong> / <strong>Dependencies</strong> (Uses+Used-by+Reads+Writes) / <strong>Failure Mode &amp; Diagnostics</strong>。</p>

<table>
<thead><tr><th>Skill</th><th>Group</th><th>Lines</th><th>缺項</th></tr></thead>
<tbody>
"""
non_tool_d2 = [r for r in results if r["d2_score"] == "❌" and r["name"] not in TOOL_SKILLS]
for r in sorted(non_tool_d2, key=lambda x: (group_order.get(x["group"],99), x["name"])):
    missing = []
    if not r["d2_phase"]: missing.append("Phase chain")
    if not r["d2_dependencies"]: missing.append("Deps")
    if not r["d2_failure"]: missing.append("Failure")
    html += f'<tr><td><code>{r["name"]}</code></td><td>{r["group"]}</td><td>{r["line_count"]}</td><td>{", ".join(missing)}</td></tr>\n'
html += """</tbody></table>

<h2>🟡 P1 — D3 cross-ref 缺 (5 非工具類)</h2>
<table>
<thead><tr><th>Skill</th><th>Group</th><th>建議 cross-ref hub</th></tr></thead>
<tbody>
"""
non_tool_d3 = [r for r in results if r["d3_score"] == "❌" and r["name"] not in TOOL_SKILLS]
hint_map = {"cpp-change": "/scientific-rigor §6 消融原則",
            "inject-hypothesis": "/scientific-rigor §7.1 Pre-registration + §8.3.1 Reopen Threshold",
            "observation-analysis": "/scientific-rigor §5 對照組 + /auc-confound-guard 3-gate",
            "results-analysis": "/scientific-rigor §3 Effect Size + §2 Evidence Tier",
            "review-evidence": "/scientific-rigor §2 Evidence Tier + §8.4 Provenance audit"}
for r in non_tool_d3:
    html += f'<tr><td><code>{r["name"]}</code></td><td>{r["group"]}</td><td>{hint_map.get(r["name"], "TBD")}</td></tr>\n'

html += """</tbody></table>

<h2>📋 Recommendations</h2>
<table>
<thead><tr><th>優先級</th><th>動作</th><th>影響</th><th>工時</th></tr></thead>
<tbody>
<tr><td><span class="priority-badge p0">P0</span></td><td>修 5 YAML invalid description</td><td>消除 latent skill load 失敗風險（跨 platform 一致性）</td><td>30 min</td></tr>
<tr><td><span class="priority-badge p1">P1</span></td><td>5 非工具類 D3 cross-ref 加 backlink</td><td>discoverability（用戶找得到 hub）</td><td>30 min</td></tr>
<tr><td><span class="priority-badge p1">P1</span></td><td>13 非工具類 D2 補診斷段</td><td>失敗時可診斷 + 跨 session 維護</td><td>3.5 hr (每 skill ~15 min)</td></tr>
<tr><td><span class="priority-badge p2">P2</span></td><td>12 工具類 D2/D3 補（optional）</td><td>completeness，影響低</td><td>2 hr</td></tr>
</tbody>
</table>

<h2>📊 後續審查階段（P2-P6）</h2>
<ol>
<li><strong>P2 Hooks audit</strong>: 22 hooks 衝突 / silent failure / Hard Gate 真實性</li>
<li><strong>P3 Agents audit</strong>: 16 agents role 清楚度 / 與 plugin agents 重疊</li>
<li><strong>P4 Prompt 庫 / commands / templates audit</strong>: 8 commands + 2 templates + 5 入口同步</li>
<li><strong>P5 Hooks × Skills × Settings 互動 audit</strong>: order 依賴 / silent failure / log 一致性</li>
<li><strong>P6 Documentation drift audit</strong>: 引用 broken / stale 偵測</li>
</ol>

<details>
<summary>方法論說明（click to expand）</summary>
<p>本 audit 使用機械掃描 + LLM 細查混合：</p>
<ol>
<li>Bash + Python (PyYAML) 掃 42 個 SKILL.md frontmatter / body</li>
<li>YAML invalid 觸發 regex fallback parser（避免漏判）</li>
<li>D1-D5 自動評分，D6 重複偵測留待後續人工 review</li>
<li>結果以 JSON / Markdown / HTML 三格式輸出</li>
</ol>
<p>原始數據：<code>/tmp/skill_audit_20260518/audit_results.json</code></p>
<p>掃描腳本：<code>/tmp/skill_audit_p1.py</code></p>
</details>

<hr>
<p class="meta">Generated by <code>/tmp/skill_audit_p1_render.py</code> · InterSubMod Agent Harness Phase 5 Audit</p>

</body>
</html>
"""

(OUTDIR / "audit_report.html").write_text(html, encoding="utf-8")
print(f"HTML: {OUTDIR}/audit_report.html ({len(html)} chars)")
