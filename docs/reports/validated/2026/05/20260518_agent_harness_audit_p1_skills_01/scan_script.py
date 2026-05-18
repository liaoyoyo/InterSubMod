#!/usr/bin/env python3
"""P1 Skills Audit — scan all 42 SKILL.md across 6 dimensions, output JSON + Markdown + HTML."""
import json
import os
import re
import time
import yaml
from pathlib import Path

SKILLS_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/.claude/skills")
OUTDIR = Path("/tmp/skill_audit_20260518")
OUTDIR.mkdir(parents=True, exist_ok=True)

# Skill grouping per CLAUDE.md §3
GROUPS = {
    "元方法論": ["confirmation-protocol", "known-pitfalls", "cycle-state", "grill-me",
              "research-context-loader", "fast-learning-coach", "scientific-rigor"],
    "7-Phase Waterfall": ["cycle-init", "research-loop", "check-staleness",
                          "feature-layered-observation", "multi-sample-consistency",
                          "run-evaluator", "conclude-research"],
    "程式開發": ["cpp-change", "methodology-audit", "infra-ops", "verification-loop"],
    "文件管理": ["doc-standards", "data-audit", "memory-consolidation", "citation-verification"],
    "報告生成": ["weekly-report", "pptx-build", "html-report-build", "results-report",
              "structured-tech-report", "myPPT", "report", "html-preview"],
    "研究專用": ["auc-confound-guard", "pivot-direction", "inject-hypothesis",
              "init-research", "observation-analysis", "results-analysis",
              "review-evidence", "problem-framing-ideation", "validation-protocol",
              "provenance-tier-audit", "research-dashboard", "image-gen",
              "image-vision-check"],
}

def find_group(skill_name):
    for g, sks in GROUPS.items():
        if skill_name in sks:
            return g
    return "未分類"

def parse_frontmatter(text):
    """Extract YAML frontmatter as dict + remaining body. Handles ** markdown in description (YAML-invalid)."""
    m = re.match(r"^---\s*\n(.*?)\n---\s*\n(.*)", text, re.DOTALL)
    if not m:
        return {}, text, False
    fm_text, body = m.group(1), m.group(2)
    yaml_valid = True
    try:
        fm = yaml.safe_load(fm_text) or {}
        if not isinstance(fm, dict):
            fm = {}
            yaml_valid = False
    except Exception:
        fm = {}
        yaml_valid = False

    # Fallback regex parser if YAML failed
    if not fm or not yaml_valid:
        for line in fm_text.split("\n"):
            mm = re.match(r"^([\w\-]+):\s*(.*)$", line)
            if mm:
                fm[mm.group(1)] = mm.group(2)

    fm = {k: (str(v) if not isinstance(v, str) else v) for k, v in fm.items()}
    return fm, body, yaml_valid

def audit_skill(skill_dir: Path):
    """Audit one skill across 6 dimensions."""
    name = skill_dir.name
    skill_md = skill_dir / "SKILL.md"

    if not skill_md.exists():
        return {"name": name, "error": "missing SKILL.md"}

    text = skill_md.read_text(encoding="utf-8")
    fm, body, yaml_valid = parse_frontmatter(text)
    line_count = len(text.splitlines())
    byte_size = len(text.encode("utf-8"))
    mtime_age_days = int((time.time() - skill_md.stat().st_mtime) / 86400)

    # D1 結構 — frontmatter 完整度
    d1_required = ["name", "description"]
    d1_optional = ["user-invocable", "paths", "allowed-tools"]
    d1_missing_req = [k for k in d1_required if k not in fm]
    d1_missing_opt = [k for k in d1_optional if k not in fm]
    d1_score = "✅" if not d1_missing_req else ("⚠" if len(d1_missing_req) == 1 else "❌")

    # D2 內容 — 必含段落（per feedback_skill_md_must_state_dependencies_and_diagnostics）
    d2_dependencies = bool(re.search(r"##.*Dependencies|##.*依賴|Phase.*Chain Position|##.*前置條件", body))
    d2_failure = bool(re.search(r"##.*Failure Mode|##.*Diagnostics|##.*失敗|##.*診斷", body))
    d2_phase = bool(re.search(r"Phase Chain|流程位置|工作流位置|Phase &amp; Chain", body))
    d2_present = sum([d2_dependencies, d2_failure, d2_phase])
    d2_score = "✅" if d2_present >= 2 else ("⚠" if d2_present == 1 else "❌")

    # D3 cross-ref to /scientific-rigor or 元方法論
    d3_scirig = len(re.findall(r"/scientific-rigor|元方法論", body))
    d3_other_skills = len(re.findall(r"/(known-pitfalls|methodology-audit|verification-loop|"
                                      r"validation-protocol|memory-consolidation|cycle-init|"
                                      r"research-loop|conclude-research|auc-confound-guard|"
                                      r"confirmation-protocol|fast-learning-coach|"
                                      r"problem-framing-ideation|inject-hypothesis|"
                                      r"check-staleness|run-evaluator|pivot-direction|"
                                      r"init-research|pptx-build|weekly-report)", body))
    d3_score = "✅" if d3_scirig >= 1 else ("⚠" if d3_other_skills >= 2 else "❌")

    # D4 trigger 詞
    desc = fm.get("description", "")
    d4_use_when = bool(re.search(r"USE WHEN|觸發", desc))
    d4_skip_when = bool(re.search(r"SKIP WHEN|DO NOT USE WHEN|不適用", desc))
    d4_score = "✅" if (d4_use_when and d4_skip_when) else ("⚠" if d4_use_when else "❌")

    # D5 staleness
    if mtime_age_days < 30:
        d5_score = "✅"
    elif mtime_age_days < 90:
        d5_score = "⚠"
    else:
        d5_score = "❌"

    # D6 重複/重疊（heuristic — 用 description keyword overlap，後續人工 audit）
    # 此處只記錄 description 第一個關鍵動詞短語供後續人工比對
    d6_keyword = re.search(r"^([^。.]+?)(?:[。.]|USE WHEN)", desc)
    d6_first_clause = d6_keyword.group(1).strip()[:60] if d6_keyword else desc[:60]

    return {
        "name": name,
        "group": find_group(name),
        "line_count": line_count,
        "byte_size": byte_size,
        "mtime_age_days": mtime_age_days,
        "desc_len": len(desc),
        "d1_score": d1_score,
        "d1_missing_req": d1_missing_req,
        "d1_missing_opt": d1_missing_opt,
        "d2_score": d2_score,
        "d2_dependencies": d2_dependencies,
        "d2_failure": d2_failure,
        "d2_phase": d2_phase,
        "d3_score": d3_score,
        "d3_scirig_count": d3_scirig,
        "d3_other_skill_count": d3_other_skills,
        "d4_score": d4_score,
        "d4_use_when": d4_use_when,
        "d4_skip_when": d4_skip_when,
        "d5_score": d5_score,
        "d6_first_clause": d6_first_clause,
        "user_invocable": fm.get("user-invocable", "—"),
        "paths_scoped": "paths" in fm,
        "yaml_valid": yaml_valid,
    }

# Run audit
results = []
for skill_dir in sorted(SKILLS_DIR.iterdir()):
    if skill_dir.is_dir():
        results.append(audit_skill(skill_dir))

# Summary stats
def count_score(key, sym):
    return sum(1 for r in results if r.get(key) == sym)

summary = {
    "total_skills": len(results),
    "scores": {
        d: {"✅": count_score(f"{d}_score", "✅"),
            "⚠": count_score(f"{d}_score", "⚠"),
            "❌": count_score(f"{d}_score", "❌")}
        for d in ["d1", "d2", "d3", "d4", "d5"]
    },
    "user_invocable_count": sum(1 for r in results if r.get("user_invocable") == "true"),
    "paths_scoped_count": sum(1 for r in results if r.get("paths_scoped")),
    "yaml_invalid_count": sum(1 for r in results if r.get("yaml_valid") is False),
    "yaml_invalid_skills": [r["name"] for r in results if r.get("yaml_valid") is False],
    "avg_line_count": int(sum(r.get("line_count", 0) for r in results) / max(len(results), 1)),
    "avg_desc_len": int(sum(r.get("desc_len", 0) for r in results) / max(len(results), 1)),
}

# Priority fix list — sort by total ❌+⚠ count
def total_issues(r):
    return sum(1 for d in ["d1","d2","d3","d4","d5"] if r.get(f"{d}_score") == "❌") * 2 \
         + sum(1 for d in ["d1","d2","d3","d4","d5"] if r.get(f"{d}_score") == "⚠")

ranked = sorted(results, key=total_issues, reverse=True)
priority_fixes = [r for r in ranked if total_issues(r) >= 3][:15]

# Save JSON
with open(OUTDIR / "audit_results.json", "w", encoding="utf-8") as f:
    json.dump({"summary": summary, "results": results, "priority_fixes": priority_fixes},
              f, ensure_ascii=False, indent=2)

print(f"Saved JSON: {OUTDIR}/audit_results.json")
print(f"Total skills audited: {len(results)}")
print(f"Summary scores:")
for d in ["d1","d2","d3","d4","d5"]:
    s = summary["scores"][d]
    print(f"  {d}: ✅{s['✅']} / ⚠{s['⚠']} / ❌{s['❌']}")
print(f"Priority fixes (top 15): {len(priority_fixes)}")
