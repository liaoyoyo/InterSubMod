#!/usr/bin/env python3
"""P3 Methodology Alignment Audit — 17 業界方法論 vs InterSubMod 體現程度。"""
import json
from datetime import datetime
from pathlib import Path

OUTDIR = Path("/tmp/skill_audit_20260518")

# Methodology coverage matrix
# Format: (id, name, evidence_locations, coverage_level, gap_description, recommendation, priority, rationale)
METHODOLOGIES = [
    # id, name, locations, level, gap, recommendation, priority, ROI rationale
    ("M01", "PDCA 循環",
     ["/scientific-rigor §9 持續改進迴圈 (Plan-Do-Check-Act)"],
     "STRONG",
     "—",
     "已落地。Plan → Pre-registration / Do → Step→Verify execution / Check → /validation-protocol L1-L4 / Act → §9.2 SRE Postmortem",
     "—", "重複既有，不必新建"),

    ("M02", "5W2H 提問",
     ["/problem-framing-ideation §1 (5W1H 框架: What/Why/Who/When/Where/How)",
      "references/5w1h-framework.md"],
     "MEDIUM",
     "缺第 7 字 H = How Much（成本/時間估算）",
     "升 5W1H → 5W2H：在 /problem-framing-ideation Step 1 加「How Much」欄（成本估 + 工時估）",
     "P1", "現有 5W1H 缺 cost dimension；對研究 cycle pre-reg 有實質幫助（避免低估投入時間）"),

    ("M03", "SMART 法則",
     ["/problem-framing-ideation references/research-question-design.md (mentioned)",
      "templates/research_index.md Pre-reg 3 欄 (隱含 measurable + time-bound)"],
     "MEDIUM",
     "未在 Pre-registration template 強制 5 字 (Specific/Measurable/Achievable/Relevant/Time-bound) 對齊",
     "templates/research_index.md Pre-reg 3 欄旁加 SMART 5-character mapping checkbox",
     "P1", "Pre-reg 已有 measurable threshold + time deadline；補 SMART 對齊使顯式 + 自我審查"),

    ("M04", "SWOT 分析",
     ["/problem-framing-ideation (mentioned in references)",
      "/methodology-audit Step 3 OPTIONS (隱含 trade-off 比較)"],
     "WEAK",
     "/methodology-audit Step 3 沒明示 SWOT 4 象限結構",
     "/methodology-audit Step 3 OPTIONS template 加 SWOT 2x2 matrix（S 內優 / W 內劣 / O 外機 / T 外威）",
     "P2", "Step 3 OPTIONS 目前只列方案不結構化比較；SWOT 強制 4 象限提升決策品質"),

    ("M05", "30 秒法則 (TL;DR Lead)",
     ["AGENTS.md §15 回應分級（隱含）",
      "CLAUDE.md 無顯式"],
     "WEAK",
     "複雜任務 AI 回應開頭未強制 1-sentence summary",
     "AGENTS.md §15 回應分級加「Tier 3+ 任務強制 TL;DR 第一句」規則",
     "P2", "降低用戶 cognitive load；對 long-form report 特別有效"),

    ("M06", "二八原理 (Pareto)",
     ["/fast-learning-coach Step 2 帕雷托分析",
      "/scientific-rigor §0.5 最小可用子集 (高影響 20% 章節)",
      "feedback_small_scale_validation_first.md (個人 anchor)"],
     "STRONG",
     "—",
     "已落地",
     "—", "重複既有"),

    ("M07", "6S 現場管理",
     ["/data-audit (整理 raw data 階段性 cleanup)",
      "/doc-standards (文件命名與目錄結構)"],
     "MEDIUM",
     "6S 名詞未顯式（Seiri/Seiton/Seiso/Seiketsu/Shitsuke/Safety）",
     "不建議獨立 skill；可在 /infra-ops 加 §6S checklist 月度跑（評估後）",
     "P3", "重複既有 skill；6S 是製造業術語，研究專案場景中價值有限"),

    ("M08", "管理 4R 法則",
     ["evidence_ledger (Result tracking)",
      "/run-evaluator (Review)"],
     "WEAK",
     "4R 全字 (Result/Responsibility/Routine/Review) 未明示",
     "不建議獨立 skill；可在 weekly-report 強制 4R 4 段（已有類似）",
     "P3", "與 weekly-report 重疊；4R 對個人工作管理有效，AI agent flow 部分重複"),

    ("M09", "高效人士七個習慣",
     ["—（範圍 mismatch）"],
     "N/A",
     "領導力/個人成長 framework，非 AI agent 適用範圍",
     "不採納",
     "—", "Stephen Covey framework 為個人領導力，scope 與 agent harness 不符"),

    ("M10", "時間管理四象限 (Eisenhower)",
     ["—（完全缺）"],
     "MISSING",
     "無任何 skill / template 提及四象限",
     "評估新 skill `/task-prioritize` 或加入 myPPT 風格的入口 skill：4 象限 (Urgent×Important) triage",
     "P3", "對研究 backlog 管理有實際價值；但每日 triage 應由用戶手動做，AI 提供 template 即可"),

    ("M11", "七階段問題解決 (提問→拆解→執行→紀錄→判讀→改進→固化)",
     ["7-Phase Waterfall: P0 /cycle-init → P1 /research-loop → P2 /check-staleness → "
      "P3 /feature-layered-observation → P4 /multi-sample-consistency → P5 /run-evaluator → P6 /conclude-research"],
     "STRONG",
     "—",
     "已落地（7-Phase = 七階段對應）",
     "—", "重複既有"),

    ("M12", "核心能力 6 步 (模糊→清楚→拆小→記錄→判讀→改善→流程)",
     ["/scientific-rigor §11 12 步協作圖",
      "/problem-framing-ideation (模糊→清楚)",
      "/methodology-audit (拆小)",
      "evidence_ledger (記錄)",
      "/run-evaluator (判讀)",
      "/scientific-rigor §9 (改善+流程化)"],
     "STRONG",
     "—",
     "已落地",
     "—", "重複既有"),

    ("M13", "紀錄格式 7 欄 (日期/目標/做法/結果/問題/原因/下一步)",
     ["evidence_ledger.jsonl (5 欄: cycle_id/hypothesis/delta_f1/decision/key_observations)",
      "templates/postmortem.md (8 段含 timeline/root cause/action items)"],
     "MEDIUM",
     "evidence_ledger 缺 2 欄: 「問題」「下一步」 (action items)",
     "evidence_ledger.jsonl schema 加 next_action + identified_issues 2 欄",
     "P2", "現有結構是 ledger 風格（過去式 verdict），加 action 欄使可 trace 後續決策"),

    ("M14", "解決問題 7 步 (描述→影響→原因→排除→驗證→行動→檢查)",
     ["/scientific-rigor §4 DAG (描述+原因)",
      "/methodology-audit Step 1-3 (排除+驗證)",
      "/verification-loop (行動+檢查)"],
     "STRONG",
     "—",
     "已落地（散在 3 個 skill）",
     "—", "重複既有"),

    ("M15", "提問格式 6 欄 (目標/已做/問題/懷疑/限制/期望答案)",
     ["—（無顯式 prompt template）",
      "用戶日常以對話形式提問，AI 主動 5W1H 反問補齊"],
     "WEAK",
     "用戶提問時無強制模板；只在 AI 反問時提供",
     "建立 InterSubMod/templates/user_query_template.md 6 欄模板（可選用，不強制）",
     "P0", "對用戶端最大價值（一次提問減少 3-5 輪 AI 反問），最低修改成本"),

    ("M16", "完整做事 9 步",
     ["7-Phase Waterfall + /problem-framing-ideation + /scientific-rigor 完整協作"],
     "STRONG",
     "—",
     "已落地",
     "—", "重複既有"),

    ("M17", "9 個 do 原則 (DoD/拆小/MVP/80-20/timebox/checklist/PDCA/復盤/分類)",
     ["DoD: AGENTS.md §6 Step→Verify",
      "拆小: /cpp-change 6 steps",
      "MVP: known-pitfalls + cpp-change + memory-consolidation 等 6 skill",
      "80-20: /fast-learning-coach + /scientific-rigor §0.5",
      "timebox: feedback_small_scale_validation_first (<2hr)",
      "checklist: 8 skill 含 quality checklist",
      "PDCA: /scientific-rigor §9",
      "復盤: /scientific-rigor §9.2 + templates/postmortem.md (5 skill)",
      "分類: /data-audit + /doc-standards"],
     "STRONG",
     "—",
     "已落地（9 條全覆蓋）",
     "—", "重複既有"),
]

# Score summary
score_count = {"STRONG": 0, "MEDIUM": 0, "WEAK": 0, "MISSING": 0, "N/A": 0}
for m in METHODOLOGIES:
    score_count[m[3]] += 1

# Priority count
priority_count = {"P0": 0, "P1": 0, "P2": 0, "P3": 0, "—": 0}
for m in METHODOLOGIES:
    priority_count[m[6]] = priority_count.get(m[6], 0) + 1

# Color mapping
COLOR = {
    "STRONG": "#1a7f37",
    "MEDIUM": "#bf8700",
    "WEAK": "#fb8500",
    "MISSING": "#cf222e",
    "N/A": "#6e7781",
}
PRI_COLOR = {"P0": "#cf222e", "P1": "#fb8500", "P2": "#bf8700", "P3": "#0969da", "—": "#6e7781"}

# Render HTML
html = f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>P3 Methodology Alignment Audit — InterSubMod</title>
<style>
  body {{ font-family: -apple-system, "Segoe UI", "Microsoft JhengHei", "Noto Sans CJK TC", sans-serif;
         max-width: 1500px; margin: 2em auto; padding: 0 1em; line-height: 1.55; color: #1f2328; background: #fafbfc; }}
  h1 {{ border-bottom: 2px solid #d0d7de; padding-bottom: 0.3em; }}
  h2 {{ border-bottom: 1px solid #d0d7de; padding-bottom: 0.2em; margin-top: 2em; }}
  .meta {{ color: #57606a; font-size: 0.9em; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.85em; margin: 1em 0; background: #fff; }}
  th, td {{ border: 1px solid #d0d7de; padding: 6px 10px; text-align: left; vertical-align: top; }}
  th {{ background: #f6f8fa; font-weight: 600; }}
  .level {{ display: inline-block; padding: 2px 8px; border-radius: 3px; color: white; font-weight: bold; font-size: 0.85em; }}
  .priority {{ display: inline-block; padding: 2px 6px; border-radius: 3px; color: white; font-weight: bold; font-size: 0.8em; }}
  code {{ background: #f6f8fa; padding: 1px 4px; border-radius: 3px; font-family: ui-monospace,monospace; font-size: 0.9em; }}
  .summary-grid {{ display: grid; grid-template-columns: repeat(5, 1fr); gap: 0.5em; margin: 1em 0; }}
  .summary-cell {{ background: #fff; border: 1px solid #d0d7de; padding: 0.7em; border-radius: 4px; text-align: center; }}
  .summary-cell .num {{ font-size: 1.5em; font-weight: bold; }}
  .recommendation {{ background: #ddf4ff; padding: 0.5em 0.8em; border-radius: 3px; }}
  .gap {{ background: #fff8c5; padding: 0.5em 0.8em; border-radius: 3px; }}
  .rationale {{ font-size: 0.85em; color: #57606a; font-style: italic; }}
  ul, ol {{ margin: 0.3em 0; padding-left: 1.5em; }}
  details {{ background: #fff; border: 1px solid #d0d7de; border-radius: 4px; padding: 0.5em 1em; margin: 0.5em 0; }}
  summary {{ cursor: pointer; font-weight: bold; }}
</style></head>
<body>
<h1>P3 Methodology Alignment Audit — 17 業界方法論對照</h1>
<p class="meta">Audit date: {datetime.now().strftime('%Y-%m-%d %H:%M')} ·
Scope: <em>提問/拆解/執行/紀錄/判讀/改進/固化</em> 7 階段元方法 +
17 specific frameworks (PDCA, 5W2H, SMART, SWOT, Pareto, 6S, 4R, Eisenhower, 7-habits...)</p>

<h2>1. Coverage Level Summary</h2>
<div class="summary-grid">
"""
for lvl, cnt in score_count.items():
    color = COLOR[lvl]
    html += f'<div class="summary-cell"><div style="font-size:0.85em;color:#57606a">{lvl}</div>'
    html += f'<div class="num" style="color:{color}">{cnt}</div>'
    html += f'<div style="font-size:0.8em;color:#57606a">/{len(METHODOLOGIES)}</div></div>\n'
html += """</div>

<p><strong>解讀</strong>：</p>
<ul>
<li><strong style="color:#1a7f37">STRONG (7/17 = 41%)</strong>: 已完整落地，重複既有，不必新建</li>
<li><strong style="color:#bf8700">MEDIUM (3/17 = 18%)</strong>: 部分落地，可補強深度</li>
<li><strong style="color:#fb8500">WEAK (4/17 = 24%)</strong>: 隱含或淺，需明文化</li>
<li><strong style="color:#cf222e">MISSING (1/17 = 6%)</strong>: 完全缺（四象限）</li>
<li><strong style="color:#6e7781">N/A (1/17 = 6%)</strong>: 範圍 mismatch（七習慣）</li>
</ul>

<h2>2. 完整對照表 (17 methodologies × 6 columns)</h2>
<table>
<thead>
<tr>
  <th>ID</th><th>方法論</th><th>體現位置</th><th>強度</th><th>Gap</th><th>建議落地</th><th>優先級</th>
</tr>
</thead>
<tbody>
"""
for m in METHODOLOGIES:
    mid, name, locs, lvl, gap, rec, pri, rat = m
    color = COLOR[lvl]
    pri_color = PRI_COLOR.get(pri, "#6e7781")
    loc_html = "<ul style='margin:0;padding-left:1.2em'>" + "".join(f"<li>{l}</li>" for l in locs) + "</ul>"
    html += f"""<tr>
  <td><strong>{mid}</strong></td>
  <td><strong>{name}</strong></td>
  <td>{loc_html}</td>
  <td><span class="level" style="background:{color}">{lvl}</span></td>
  <td>{gap}</td>
  <td>{rec}</td>
  <td><span class="priority" style="background:{pri_color}">{pri}</span></td>
</tr>
"""
html += """</tbody></table>

<h2>3. 落地決策 — P0 必做 (1 項, 30 min)</h2>
<table>
<thead><tr><th>ID</th><th>方法論</th><th>實作位置</th><th>實作內容</th><th>工時</th></tr></thead>
<tbody>
<tr>
<td><strong>M15</strong></td>
<td>提問格式 6 欄</td>
<td><code>InterSubMod/templates/user_query_template.md</code> (新建)</td>
<td>
6 欄模板：目標 / 已做 / 問題 / 懷疑原因 / 限制 / 期望答案形式。
不強制（用戶可選用）；在 CLAUDE.md §6 或 README 提供連結供用戶複製。
</td>
<td>30 min</td>
</tr>
</tbody>
</table>
<p class="rationale">為什麼 P0？對用戶端 ROI 最高（一次提問減少 3-5 輪 AI 反問），實作最簡單（純文件模板）。</p>

<h2>4. 落地決策 — P1 本週做 (2 項, 1hr)</h2>
<table>
<thead><tr><th>ID</th><th>方法論</th><th>實作位置</th><th>實作內容</th><th>工時</th></tr></thead>
<tbody>
<tr>
<td><strong>M02</strong></td>
<td>5W2H 升級</td>
<td><code>InterSubMod/.claude/skills/problem-framing-ideation/SKILL.md</code> §1 Step 1</td>
<td>5W1H 表加第 7 欄「How Much (cost/time)」— 估算 USD 成本 + 工時上限</td>
<td>30 min</td>
</tr>
<tr>
<td><strong>M03</strong></td>
<td>SMART 強制對齊</td>
<td><code>InterSubMod/templates/research_index.md</code> Pre-reg 3 欄</td>
<td>Pre-reg 旁加 SMART 5-字 self-checklist：[ ] Specific [ ] Measurable [ ] Achievable [ ] Relevant [ ] Time-bound</td>
<td>30 min</td>
</tr>
</tbody>
</table>

<h2>5. 落地決策 — P2 下兩週做 (3 項, 2hr)</h2>
<table>
<thead><tr><th>ID</th><th>方法論</th><th>實作位置</th><th>實作內容</th><th>工時</th></tr></thead>
<tbody>
<tr>
<td><strong>M04</strong></td>
<td>SWOT 4 象限</td>
<td><code>InterSubMod/.claude/skills/methodology-audit/SKILL.md</code> Step 3 OPTIONS</td>
<td>方案比較表加 SWOT 2x2 matrix template (Strengths/Weaknesses/Opportunities/Threats)</td>
<td>45 min</td>
</tr>
<tr>
<td><strong>M05</strong></td>
<td>30 秒法則 (TL;DR)</td>
<td><code>InterSubMod/AGENTS.md</code> §15 回應分級</td>
<td>加「Tier 3+ 複雜任務回應首句強制 TL;DR (≤30 sec read)」規則</td>
<td>15 min</td>
</tr>
<tr>
<td><strong>M13</strong></td>
<td>紀錄 7 欄補完</td>
<td><code>InterSubMod/research/autoresearch/evidence_ledger_schema.md</code> + scripts</td>
<td>evidence_ledger.jsonl schema 加 next_action + identified_issues 2 欄</td>
<td>1hr (schema migration + 既有 entry backfill)</td>
</tr>
</tbody>
</table>

<h2>6. 落地決策 — P3 評估後再做 (3 項)</h2>
<details>
<summary>M07 — 6S 現場管理</summary>
<p><strong>建議</strong>：不獨立 skill；在 <code>InterSubMod/.claude/skills/infra-ops/SKILL.md</code> 加 §6S 月度 checklist (整理 / 整頓 / 清掃 / 清潔 / 教養 / 安全)。</p>
<p class="rationale"><strong>理由</strong>：與 /data-audit + /doc-standards 重疊；6S 名詞為製造業術語，研究專案 ROI 邊際。</p>
</details>
<details>
<summary>M08 — 4R 法則</summary>
<p><strong>建議</strong>：不獨立 skill；在 <code>InterSubMod/.claude/skills/weekly-report/SKILL.md</code> 加 4R 章節對齊 (Result / Responsibility / Routine / Review)。</p>
<p class="rationale"><strong>理由</strong>：weekly-report 已隱含類似結構；4R 強調個人責任，AI 場景部分轉化為 cycle_id 追蹤。</p>
</details>
<details>
<summary>M10 — 時間管理四象限 (Eisenhower)</summary>
<p><strong>建議</strong>：新建 <code>InterSubMod/templates/task_priority_matrix.md</code> 模板（不建 skill），供用戶日常 backlog triage。</p>
<p class="rationale"><strong>理由</strong>：每日 triage 是用戶端工作，AI 提供 template 即可；建獨立 skill 過 over-engineering。</p>
</details>

<h2>7. 不採納 (1 項)</h2>
<details>
<summary>M09 — 高效人士七個習慣 (Stephen Covey)</summary>
<p><strong>不採納理由</strong>：Covey 框架為個人領導力 / 人際關係 framework（主動積極 / 以終為始 / 要事第一 / 雙贏 / 知彼解己 / 統合綜效 / 不斷更新）。Scope 與 InterSubMod agent harness (technical agent flow) 不符；強行對齊會稀釋既有方法論。</p>
<p class="rationale">業界對照：Anthropic Cookbook / Claude Code 文件均未引用此框架；OpenAI / Cursor / Devin 等 agent product 也不採納。</p>
</details>

<h2>8. 業界格式 vs 教學對話格式 — 元層說明</h2>
<table>
<thead><tr><th>場景</th><th>用業界嚴謹格式</th><th>用教學對話格式</th></tr></thead>
<tbody>
<tr><td>研究結論宣告 / 證據評估</td><td>✅ default</td><td>—</td></tr>
<tr><td>方案比較 / 決策表</td><td>✅ default</td><td>—</td></tr>
<tr><td>程式碼修改 / verification</td><td>✅ default</td><td>—</td></tr>
<tr><td>狀態盤點 / progress report</td><td>✅ default</td><td>—</td></tr>
<tr><td>新概念解釋（用戶不熟）</td><td>—</td><td>✅ 觸發 /fast-learning-coach</td></tr>
<tr><td>用戶明示「不懂」「解釋一下」</td><td>—</td><td>✅</td></tr>
<tr><td>方法論討論（如本報告）</td><td>✅</td><td>—</td></tr>
</tbody></table>

<h2>9. Implementation Schedule</h2>
<p>合併 P1 Skills audit + P2 Industry audit + 本 P3 Methodology audit 的所有 fix items：</p>
<table>
<thead><tr><th>Phase</th><th>項目</th><th>工時</th><th>累計</th></tr></thead>
<tbody>
<tr><td><span class="priority" style="background:#cf222e">P0</span></td><td>5 YAML invalid + 6 silent failure hooks + M15 提問模板</td><td>1.5 hr</td><td>1.5 hr</td></tr>
<tr><td><span class="priority" style="background:#fb8500">P1</span></td><td>13 D2 診斷段 + 5 D3 cross-ref + allow audit + M02 5W2H + M03 SMART</td><td>5 hr</td><td>6.5 hr</td></tr>
<tr><td><span class="priority" style="background:#bf8700">P2</span></td><td>cache_telemetry + compact_test + subagent_logger + M04 SWOT + M05 TL;DR + M13 ledger 7 欄</td><td>3 hr</td><td>9.5 hr</td></tr>
<tr><td><span class="priority" style="background:#0969da">P3</span></td><td>router_pre_bash + worktree default + rules paths + skill desc 縮短 + cost_summary + M07/08/10 評估</td><td>4 hr</td><td>13.5 hr</td></tr>
<tr><td><span class="priority" style="background:#6e7781">P4</span></td><td>injection + memory archive + cross-version baseline + modify hook 示範</td><td>3 hr</td><td>16.5 hr</td></tr>
</tbody></table>

<hr>
<p class="meta">Generated by <code>/tmp/p3_methodology_alignment.py</code> · Data: <code>/tmp/skill_audit_20260518/p3_methodology_alignment.json</code></p>
</body></html>
"""

(OUTDIR / "p3_methodology_alignment.html").write_text(html, encoding="utf-8")

# Also save JSON
out_data = {
    "audit_date": datetime.now().strftime("%Y-%m-%d %H:%M"),
    "methodologies": [
        {"id": m[0], "name": m[1], "locations": m[2], "level": m[3],
         "gap": m[4], "recommendation": m[5], "priority": m[6], "rationale": m[7]}
        for m in METHODOLOGIES
    ],
    "score_summary": score_count,
    "priority_summary": priority_count,
}
(OUTDIR / "p3_methodology_alignment.json").write_text(json.dumps(out_data, ensure_ascii=False, indent=2))

print(f"HTML: {OUTDIR}/p3_methodology_alignment.html ({len(html)} chars)")
print(f"JSON: {OUTDIR}/p3_methodology_alignment.json")
print(f"\nCoverage summary:")
for lvl, cnt in score_count.items():
    print(f"  {lvl:<10} {cnt:>3}/{len(METHODOLOGIES)} ({cnt*100//len(METHODOLOGIES)}%)")
print(f"\nPriority distribution:")
for pri, cnt in priority_count.items():
    print(f"  {pri:<6} {cnt}")
