#!/usr/bin/env python3
"""harness_health.py — InterSubMod harness self-audit + architecture dashboard (read-only).

Generates a single standalone HTML dashboard that lets a solo researcher, on every
harness change, confirm BOTH the moving parts AND the whole architecture:

  1. L0 health lights (8)  — count drift / Hard-Gate truth / tier-gate / state<->doc / ledger / queue / memory-drift / doc-path-currency
  2. Architecture map      — 3 entrypoints · 7-Phase Waterfall · meta + 3-storey assurance ·
                             18 agents by role · 38 hooks by event · state machine · evidence chain · workflows
  3. Component detail       — skills by category · agent table (tools/model/isolation) · hooks by event
  4. Hard-Gate panel        — real exit-2 vs advisory vs neutered (|| exit 0)
  5. Recent issues          — postmortem tails + recent commits
  6. Industry alignment     — which top patterns are ALREADY_COVERED / SKIP (restraint verdicts)
  7. How to update          — regenerate cmd + localStorage diff cadence

NO agents, NO web, NO writes outside state/health_snapshots/. Pure grep/parse of disk truth,
so the dashboard is always-current: re-run after any change to see what moved (vs last snapshot).

Run: python3 scripts/harness_health.py     Design: docs/references/migration/20260531_harness_audit_dashboard_design_01.md
"""
import json
import os
import re
import glob
import html as _html
import datetime
import subprocess

REPO = "/big7_disk/liaoyoyo2001/InterSubMod"
SNAP_DIR = os.path.join(REPO, "state", "health_snapshots")
CLAUDE_MD = os.path.join(REPO, ".claude", "CLAUDE.md")
SETTINGS = os.path.join(REPO, ".claude", "settings.local.json")
CURRENT_FOCUS = os.path.join(REPO, "docs", "CURRENT_FOCUS.md")
LEDGER = os.path.join(REPO, "research", "autoresearch", "evidence_ledger.jsonl")
ACTIVE = os.path.join(REPO, "state", "active.json")
QUEUE = os.path.join(REPO, "research", "autoresearch", "hypothesis_queue.json")
MEM_DIR = "/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory"
COMPILE_MARKER = "/tmp/ism_cpp_pending_compile.txt"
POSTMORTEMS = ["evidence_gate_bypass.log", "hook_failures.log", "injection_attempts.log",
               "recall_log_202605.log", "skill_audit_202605.log", "subagent_completion_202605.log"]

COLOR = {"GREEN": "#1a9850", "YELLOW": "#e0a800", "RED": "#d73027", "GREY": "#8a8f98", "BLUE": "#2c6fbb"}

# ── curated STABLE structure (verified against disk where possible) ───────────
ENTRYPOINTS = [
    ("AGENTS.md", "跨 agent governance（語言/結構/build/目標/KB/輸出/回應分級）", "concatenate"),
    (".claude/CLAUDE.md", "Claude Code 模式（確認矩陣/Skills/Hooks/Rules/subagent）", "auto-load"),
    ("docs/CURRENT_FOCUS.md", "live 主軸/阻塞/active cycle", "SessionStart hook 注入"),
]
WATERFALL = [
    ("P0", "cycle-init", "REGISTER"), ("P1", "research-loop", "PLAN"), ("P2", "check-staleness", "PRECHECK"),
    ("P3", "feature-layered-observation", "PILOT"), ("P4", "multi-sample-consistency", "GENERALIZE"),
    ("P5", "run-evaluator", "EVALUATE"), ("P6", "conclude-research", "CONCLUDE"),
]
ASSURANCE = [("pre ≤30min", "pre-decision-audit", "GO/PROBE/NO-GO"),
             ("process live", "implementation-notes", "設計決定/偏離/折衷"),
             ("post P5", "run-evaluator", "6-component retraction risk")]
AGENT_BUCKET = {
    "architect": "規劃/設計", "research-orchestrator": "規劃/設計",
    "developer": "生成(C++/發布)", "optimizer": "生成(C++/發布)", "release": "生成(C++/發布)",
    "evaluator": "驗證 (read-only)", "methodology-reviewer": "驗證 (read-only)",
    "reproducibility-audit": "驗證 (read-only)", "security-reviewer": "驗證 (read-only)", "reviewer": "驗證 (read-only)",
    "tester": "平行/執行", "parallel-benchmark": "平行/執行", "parallel-analysis": "平行/執行", "headless-research": "平行/執行",
    "researcher": "研究/文獻", "literature-reviewer": "研究/文獻", "paper-miner": "研究/文獻",
    "narrative-organizer": "敘述萃取",
}
BUCKET_ORDER = ["規劃/設計", "生成(C++/發布)", "驗證 (read-only)", "平行/執行", "研究/文獻", "敘述萃取"]
EVENT_ORDER = ["SessionStart", "UserPromptSubmit", "PreToolUse", "PostToolUse", "SubagentStop", "Stop", "PreCompact"]
INDUSTRY = [
    ("理解 AI 產生 (observability)", "ALREADY_COVERED",
     "cache_telemetry + subagent_completion_logger + evidence_ledger L1-L5 + Step→Verify",
     "SKIP: Langfuse/LangSmith/OTel SaaS（self-host=disk-full；flat-rate 下 cost=0）"),
    ("驗證 AI 正確 (verification)", "ALREADY_COVERED",
     "4 read-only verifier + deterministic 6-component run-evaluator + kb_schema exit-2 + provenance-tier-audit",
     "SKIP: 通用 LLM-as-judge（self-preference bias / prose brittle keyword match）"),
    ("評估與計劃 (planning)", "ALREADY_COVERED",
     "research-loop plan.json + pre-decision-audit(Pre-Mortem) + confirmation-protocol(Cynefin) + cpp-change PDD",
     "SKIP: Spec Kit/Kiro（dual-SoT 衝突）/ self-reflection loop（業界證實會退化）"),
    ("orchestration 編排", "ALREADY_COVERED",
     "7-Phase chaining + task_type 路由 + parallel-* + 2 Dynamic Workflow + evaluator-optimizer fleet",
     "KEEP: §8『Hard-Gate work 絕不入 workflow』是最重要的業界級保護，勿放寬"),
]

# 限制與提示治理盤點（2026-06-01；curated 分析 — 每條：機制 / 當初原因 / 合理性 / 改進）
GOVERNANCE = [
    ("硬·deny", "BLOCK", ".env/.ssh/secret/*.pem 讀取 · rm -rf // · mkfs/fdisk",
     "標準機密外洩 + 毀滅性操作防護", "✅ 合理 keep", "無"),
    ("硬·ask", "PROMPT", "rm/mv/chmod/ln/dd · git push/reset/restore/rebase · curl/wget/ssh/sudo",
     "不可逆操作前人工確認（§1 可逆性 override）", "✅ 合理", "git restore --staged 等可逆操作偶誤 prompt；可細分 allow 安全形式（低優先）"),
    ("硬·exit-2 gate", "BLOCK", "pre_commit_compile · kb_schema · pipeline_block · no_binary · search_scope · pre_tier_upgrade",
     "C++ regression / KB 完整 / tmp 800GB disk-full / repo bloat / 搜尋浪費 / tier over-claim 六類事故",
     "✅ 全合理（5-31 修好 neutering 後全生效）", "/harness-health 燈#2 持續監看復發"),
    ("硬規則·§1", "RULE", "刪檔 / C++ commit / NO-GO / ledger 覆寫 永遠 🔴",
     "不可逆 / 影響結論", "✅ 核心治理 keep", "無"),
    ("硬規則·§12", "BLOCK", "禁 grep -r . / find . 無 maxdepth",
     "repo 極大（GB data/output/worktree）大遞迴=卡住主因", "✅ 合理", "Grep/Glob 工具已自動跳過重目錄，規則只擋 Bash 誤用"),
    ("硬規則·§8", "RULE", "Hard-Gate work 絕不入 Dynamic Workflow",
     "workflow subagent 一律 acceptEdits + 繞過 exit-2", "✅ 關鍵保護 keep", "無"),
    ("軟提示·UserPrompt", "ADVISE", "knowledge_check · kb_freshness · narrative_frame · task_type · md_path · research_direction（6）",
     "outside-claim 查KB(slide14) / KB 時效 / 理解負擔(§11) / 5-24 scope事故 / 路徑規則 / 防 reopen concluded",
     "✅ 合理但廣關鍵字易重複/誤觸（本回合 handoff 誤判 + KB 提示無關）", "收窄關鍵字特異性 / trivial prompt 靜音（非阻擋，低優先）"),
    ("軟提示·claim 4-linter", "ADVISE", "evidence_level_lint · causal_claim_check · researcher_claim_evidence · terminology_guard",
     "防 hedge 語言 / 因果過度宣稱 / L3 當 L1 / 術語漂移", "✅ 合理（稍冗餘）", "audit 判 SKIP_RISK：working safety 重構風險 > 省 latency，維持"),
    ("軟提示·soft gate", "ADVISE", "verify_gate · kb_sot_guard（皆 exit 0 advisory）",
     "+0.057 fabrication 事故 / F1 雙口徑(slide14)", "✅ 合理（已正名 advisory 非 Hard Gate）", "無"),
    ("背景 logger", "LOG", "subagent_completion · skill_change_audit · memory_recall · evidence_read · skill_usage · cache_telemetry",
     "量化 cache-hit / 引用率 / skill 變動稽核", "✅ 純被動無行為影響 keep", "cache_telemetry USD 已標 flat-rate 無意義"),
    ("收尾提示·Stop", "ADVISE", "會話結束提醒寫報告 · SubagentStop 完成檢查",
     "防漏記錄", "🟡 合理但每次觸發稍噪", "可改條件式（有實質產出才提醒）低優先"),
]


def _read(path):
    try:
        with open(path, encoding="utf-8") as f:
            return f.read()
    except Exception:
        return ""


def _claim(pattern, text, default=None):
    nums = [int(m) for m in re.findall(pattern, text)]
    return max(nums) if nums else default


# ── disk ground truth ────────────────────────────────────────────────────────
def disk_counts():
    return {
        "skills": len(glob.glob(os.path.join(REPO, ".claude/skills/*/SKILL.md"))),
        "agents": len(glob.glob(os.path.join(REPO, ".claude/agents/*.md"))),
        "hooks": len(glob.glob(os.path.join(REPO, "scripts/hooks/*.sh"))),
        "workflows": len(glob.glob(os.path.join(REPO, ".claude/workflows/*.js"))),
    }


def parse_agents():
    out = []
    for f in sorted(glob.glob(os.path.join(REPO, ".claude/agents/*.md"))):
        txt = _read(f)
        m = re.search(r"^---\s*\n(.*?)\n---\s*\n", txt, re.S)
        fm = m.group(1) if m else txt[:1500]
        name = (re.search(r"^name:\s*(.+)$", fm, re.M) or [None, os.path.basename(f)[:-3]])[1].strip()
        model = (re.search(r"^model:\s*(.+)$", fm, re.M) or [None, "inherit"])[1].strip()
        iso = "worktree" if re.search(r"^isolation:\s*worktree", fm, re.M) else "—"
        tline = re.search(r"^tools:\s*(.+)$", fm, re.M)
        tools = tline.group(1).strip() if tline else "?"
        if tools.startswith("["):
            try:
                tools = ", ".join(json.loads(tools))
            except Exception:
                tools = tools.strip("[]")
        desc = (re.search(r"^description:\s*(.+)$", fm, re.M) or [None, ""])[1].strip().strip('"')
        role = re.split(r"USE WHEN|。|\. ", desc)[0][:62]
        deprecated = "DEPRECATED" in desc
        out.append({"name": name, "role": role, "tools": tools, "model": model,
                    "iso": iso, "bucket": AGENT_BUCKET.get(name, "其他"), "deprecated": deprecated})
    return out


def parse_skill_categories(cmd):
    cats = []
    for m in re.finditer(r"^- \*\*(.+?)（(\d+)）\*\*[:：]\s*(.+)$", cmd, re.M):
        cat, _n, body = m.group(1), m.group(2), m.group(3)
        skills = re.findall(r"`(/[a-z0-9-]+)`", body)
        if skills:
            cats.append((cat, skills))
    return cats


def parse_hooks_by_event(settings):
    by_event = {}
    for event, blocks in (settings.get("hooks", {}) or {}).items():
        rows = []
        for b in blocks:
            for h in b.get("hooks", []):
                command = h.get("command", "")
                mm = re.search(r"(/[^\s]+\.sh)", command)
                if not mm:
                    inline = re.search(r"(jq -r[^\"]+)", command)
                    rows.append({"name": "(inline)", "exit2": False, "masked": False, "soft": False})
                    continue
                spath = mm.group(1)
                body = _read(spath)
                exit2 = ("exit 2" in body) or ("gate_exit" in body)
                masked = bool(re.search(r"\|\|\s*exit 0", command) or re.search(r"\|\|\s*true", command))
                soft = bool(re.search(r"SOFT gate|advisory only|advisory\)", body))
                rows.append({"name": os.path.basename(spath), "exit2": exit2, "masked": masked, "soft": soft})
        by_event[event] = rows
    return by_event


def recent_commits(n=6):
    try:
        r = subprocess.run(["git", "-C", REPO, "log", "--oneline", f"-{n}"],
                           capture_output=True, text=True, timeout=8)
        return [l for l in r.stdout.splitlines() if l.strip()]
    except Exception:
        return []


def postmortem_tails(n=2):
    out = []
    for fn in POSTMORTEMS:
        p = os.path.join(REPO, "docs", "postmortems", fn)
        if not os.path.exists(p):
            continue
        try:
            lines = [l.rstrip() for l in open(p, encoding="utf-8") if l.strip()]
        except Exception:
            lines = []
        if lines:
            out.append((fn, lines[-n:]))
    return out


# ── L0 lights (unchanged logic) ───────────────────────────────────────────────
def light_count(disk, cmd):
    sk = _claim(r"(\d+)\s*個?\s*SKILL\.md", cmd)
    ag = _claim(r"(\d+)\s*個\s*project agent", cmd) or 18
    hk = _claim(r"\*\*(\d+)\s*hook scripts", cmd)
    rows, worst = [], "GREEN"
    for name, dval, cval in [("skills", disk["skills"], sk), ("agents", disk["agents"], ag), ("hooks", disk["hooks"], hk)]:
        if cval is None:
            rows.append(f"{name}: disk={dval} doc=?")
            continue
        diff = abs(dval - cval)
        st = "GREEN" if diff == 0 else ("YELLOW" if diff == 1 else "RED")
        if st == "RED":
            worst = "RED"
        elif st == "YELLOW" and worst == "GREEN":
            worst = "YELLOW"
        rows.append(f"{name}: disk={dval} doc={cval} {'OK' if diff==0 else f'DRIFT {diff}'}")
    return worst, rows


def light_hardgate(settings, cmd):
    real, neutered, soft_int = [], [], []
    for event, blocks in (settings.get("hooks", {}) or {}).items():
        for b in blocks:
            for h in b.get("hooks", []):
                command = h.get("command", "")
                m = re.search(r"(/[^\s]+\.sh)", command)
                if not m:
                    continue
                body = _read(m.group(1))
                has_exit2 = ("exit 2" in body) or ("gate_exit" in body)
                masked = bool(re.search(r"\|\|\s*exit 0", command) or re.search(r"\|\|\s*true", command))
                soft = bool(re.search(r"SOFT gate|advisory only|advisory\)", body))
                base = os.path.basename(m.group(1))
                if has_exit2 and not masked:
                    real.append(base)
                elif has_exit2 and masked and soft:
                    soft_int.append(base)
                elif has_exit2 and masked:
                    neutered.append((base, event))
    claim = _claim(r"校正.{0,8}=\s*\*\*(\d+)", cmd)
    rows = [f"real exit-2 (unmasked): {len(set(real))} → {sorted(set(real))}"]
    if neutered:
        rows.append("⚠ NEUTERED BUG (exit-2 但 wiring '|| exit 0' 吃掉): " + ", ".join(f"{n}({e})" for n, e in neutered))
    if soft_int:
        rows.append("intentionally-soft: " + ", ".join(sorted(set(soft_int))))
    return ("RED" if neutered else "GREEN"), rows, [n for n, _ in neutered]


def light_tiergate(raw):
    if "pre_tier_upgrade_check" in raw:
        return "GREEN", ["pre_tier_upgrade_check.sh WIRED (state.json exit-2 / 散文 advisory)"]
    if os.path.exists(os.path.join(REPO, "scripts/hooks/pre_tier_upgrade_check.sh")):
        return "RED", ["pre_tier_upgrade_check.sh 存在但 NOT wired（4 SKILL.md 宣稱它強制 ⭐4/5）"]
    return "YELLOW", ["pre_tier_upgrade_check.sh 缺"]


def light_state(active, focus):
    rows, worst, concluded = [], "GREEN", []
    for c in active.get("recently_concluded", []) if isinstance(active, dict) else []:
        concluded.append(c.get("cycle_id", ""))
    for cid in concluded:
        if cid and re.search(r"(待\s*conclude|需\s*/?conclude-research|尚未\s*conclude)", focus):
            short = cid.split("-")[-1]
            if "cycle3" in focus or short in focus:
                rows.append(f"⚠ {cid} 已 concluded 但 CURRENT_FOCUS 仍描述待 conclude")
                worst = "RED"
    if not rows:
        rows.append(f"recently_concluded={concluded or 'none'}；無明顯衝突")
    return worst, rows


def light_ledger(focus):
    last = ""
    try:
        for line in open(LEDGER, encoding="utf-8"):
            if line.strip():
                last = line
    except Exception:
        pass
    m = re.search(r'"?(?:date|timestamp|ts|concluded_at|ran_at)"?\s*[:=]\s*"?(\d{4}-\d{2}-\d{2})', last)
    led = m.group(1) if m else None
    fm = re.search(r"##\s+(?:\S+\s+)?(\d{4}-\d{2}-\d{2})", focus)
    foc = fm.group(1) if fm else None
    rows = [f"ledger last={led}；CURRENT_FOCUS header={foc}"]
    worst = "GREEN"
    if led and foc and led > foc:
        worst = "YELLOW"
        rows.append("⚠ ledger 領先 → 最新證據未回寫（tier over-claim 風險）")
    return worst, rows


def light_queue(queue):
    items = queue if isinstance(queue, list) else (queue.get("hypotheses") or queue.get("queue") or [])
    queued = [h.get("id", "?") for h in items if isinstance(h, dict) and str(h.get("status", "")).lower() == "queued"]
    rawq = json.dumps(queue, ensure_ascii=False).lower()
    dead = any(t in rawq for t in ["filter", "caller-f1", "caller_f1", "headroom"])
    rows = [f"queued={len(queued)} {queued}"]
    worst = "GREEN"
    if len(queued) >= 3:
        worst = "RED" if dead else "YELLOW"
        rows.append("⚠ queued 含 concluded-DEAD 方向 → /pivot-direction")
    return worst, rows


def light_memory_drift():
    """#7 — memory index↔files drift（升 /memory-consolidation Step-4 手動 grep 為持續偵測；read-only）。"""
    idx = _read(os.path.join(MEM_DIR, "MEMORY.md"))
    files = [os.path.basename(f) for f in glob.glob(os.path.join(MEM_DIR, "*.md"))
             if os.path.basename(f) != "MEMORY.md"]
    if not idx or not files:
        return "GREY", ["MEMORY.md / memory/ 不可讀（外部路徑）"]
    fileset = set(files)
    linked = set(re.findall(r"\(([a-z0-9_][a-z0-9_-]*\.md)\)", idx))   # [title](slug.md) bare-name links
    orphan_files = sorted(f for f in files if f not in idx)            # exists on disk but not in index
    dead_links = sorted(l for l in linked if l not in fileset)        # index points to non-existent file
    rows = [f"memory={len(files)} 檔；index links={len(linked)}"]
    worst = "GREEN"
    if dead_links:
        worst = "YELLOW"
        rows.append(f"⚠ {len(dead_links)} 索引連結指向不存在檔: {dead_links[:5]}")
    if orphan_files:
        worst = "YELLOW"
        rows.append(f"⚠ {len(orphan_files)} 檔未進 MEMORY.md 索引（不會被 recall）: {orphan_files[:5]}")
    if worst == "GREEN":
        rows.append("✓ index↔files 一致")
    return worst, rows


def light_doc_paths(cmd):
    """#8 — DOC PATH CURRENCY: CLAUDE.md 引用的 explicit scripts/*.sh|.py 路徑是否存在（references-to-deleted）。
    只查 explicit scripts/ 路徑，跳過散文裸名以避免假陽性（critique 建議）。"""
    paths = sorted(set(re.findall(r"(scripts/[A-Za-z0-9_./-]+\.(?:sh|py))", cmd)))
    missing = [p for p in paths if not os.path.exists(os.path.join(REPO, p))]
    rows = [f"CLAUDE.md 引用 {len(paths)} 個 explicit scripts/ 路徑"]
    worst = "GREEN"
    if missing:
        worst = "RED"
        rows.append(f"⚠ {len(missing)} 個不存在（references-to-deleted）: {missing[:6]}")
    else:
        rows.append("✓ 全部存在於磁碟")
    return worst, rows


def light_hook_wiring(settings, cmd):
    """#9 — HOOK WIRING: CLAUDE.md §4「Hard Gate hooks」宣稱不可繞過的 hook，是否真的在 settings
    出現。抓「宣稱 wired 但 settings 缺」的盲點 — #2 HARD-GATE TRUTH 只查 exit-2 被 '|| exit 0'
    neuter，不查 non-wiring；number_provenance 曾因此靜默漏 9 天（2026-06-10 補）。read-only grep。"""
    wired = set()
    for _ev, blocks in (settings.get("hooks", {}) or {}).items():
        for b in blocks:
            for h in b.get("hooks", []):
                m = re.search(r"hooks/(\w+\.sh)", h.get("command", ""))
                if m:
                    wired.add(m.group(1))
    sec = re.search(r"Hard Gate hooks.*?(?=\n> |\n## |\Z)", cmd, re.S)
    claimed = set(re.findall(r"`(\w+\.sh)`", sec.group(0))) if sec else set()
    if not claimed:
        return "GREY", ["CLAUDE.md §4 'Hard Gate hooks' 段不可解析"]
    missing = sorted(h for h in claimed if h not in wired)
    rows = [f"§4 宣稱 Hard Gate hooks={len(claimed)}；settings 實際 wired={len(wired)}"]
    if missing:
        rows.append(f"⚠ 宣稱 wired 但 settings 缺（false-capability，如 2026-06-10 number_provenance）: {missing}")
        return "YELLOW", rows
    rows.append("✓ 所有宣稱的 Hard Gate hook 都實際 wired 於 settings")
    return "GREEN", rows


def compile_marker_status():
    if os.path.exists(COMPILE_MARKER) and os.path.getsize(COMPILE_MARKER) > 0:
        try:
            mt = datetime.datetime.fromtimestamp(os.path.getmtime(COMPILE_MARKER)).date().isoformat()
            n = sum(1 for _ in open(COMPILE_MARKER))
        except Exception:
            mt, n = "?", "?"
        return {"present": True, "files": n, "mtime": mt}
    return {"present": False}


# ── HTML render ───────────────────────────────────────────────────────────────
def esc(s):
    return _html.escape(str(s))


def render(snap, arch):
    data = json.dumps(snap, ensure_ascii=False)
    lights = snap["lights"]
    ov = snap["overall"]
    oc = COLOR["RED"] if ov["red"] else (COLOR["YELLOW"] if ov["yellow"] else COLOR["GREEN"])

    def light_card(i, L):
        c = COLOR.get(L["status"], "#888")
        det = "<br>".join(esc(r) for r in L["rows"])
        return (f'<div class="light" style="border-left:5px solid {c}"><div class="lh">'
                f'<span class="dot" style="background:{c}"></span><b>{i} {esc(L["name"])}</b>'
                f'<span class="st" style="color:{c}">{L["status"]}</span></div>'
                f'<div class="ld">{det}</div></div>')
    lights_html = "".join(light_card(i + 1, L) for i, L in enumerate(lights))

    # architecture: entrypoints
    ep = "".join(f'<div class="box"><b>{esc(a)}</b><div class="sub">{esc(b)}</div>'
                 f'<div class="tag">{esc(c)}</div></div>' for a, b, c in ENTRYPOINTS)
    # waterfall
    wf = "".join(f'<div class="phase"><b>{esc(p)}</b><div class="sk">/{esc(s)}</div>'
                 f'<div class="sub">{esc(g)}</div></div>' +
                 ('<div class="arr">→</div>' if p != "P6" else "")
                 for p, s, g in WATERFALL)
    asr = "".join(f'<div class="box"><b>{esc(a)}</b><div class="sk">/{esc(s)}</div>'
                  f'<div class="sub">{esc(d)}</div></div>' for a, s, d in ASSURANCE)
    # agents by bucket
    agg = {}
    for a in arch["agents"]:
        agg.setdefault(a["bucket"], []).append(a)
    ag_cols = ""
    for bucket in BUCKET_ORDER + [k for k in agg if k not in BUCKET_ORDER]:
        if bucket not in agg:
            continue
        chips = ""
        for a in agg[bucket]:
            cls = "chip dep" if a["deprecated"] else "chip"
            badge = ""
            if a["iso"] == "worktree":
                badge += '<span class="b b-iso">iso</span>'
            if a["model"] == "haiku":
                badge += '<span class="b b-hk">haiku</span>'
            if "Write" not in a["tools"] and "Edit" not in a["tools"] and a["bucket"] == "驗證 (read-only)":
                badge += '<span class="b b-ro">RO</span>'
            chips += f'<span class="{cls}" title="{esc(a["tools"])}">{esc(a["name"])}{badge}</span>'
        ag_cols += f'<div class="acol"><div class="ahd">{esc(bucket)} ({len(agg[bucket])})</div>{chips}</div>'
    # hooks by event
    hk_rows = ""
    for ev in EVENT_ORDER + [k for k in arch["hooks"] if k not in EVENT_ORDER]:
        rows = arch["hooks"].get(ev)
        if not rows:
            continue
        chips = ""
        for h in rows:
            if h["name"] == "(inline)":
                chips += '<span class="chip sm">(inline)</span>'
                continue
            cls, t = "chip sm", h["name"]
            if h["exit2"] and not h["masked"]:
                cls += " gate"; t += " ⛔exit2"
            elif h["exit2"] and h["masked"] and h["soft"]:
                cls += " soft"; t += " ~soft"
            elif h["exit2"] and h["masked"]:
                cls += " bad"; t += " ⚠masked"
            chips += f'<span class="{cls}">{esc(t)}</span>'
        hk_rows += f'<tr><td class="ev">{esc(ev)} ({len(rows)})</td><td>{chips}</td></tr>'
    # skills by category
    sk_html = ""
    for cat, sks in arch["skills"]:
        chips = "".join(f'<span class="chip sm">{esc(s)}</span>' for s in sks)
        sk_html += f'<div class="acol"><div class="ahd">{esc(cat)} ({len(sks)})</div>{chips}</div>'
    # hard gate panel
    hg = snap["hardgate"]
    hg_html = (f'<p><b style="color:{COLOR["GREEN"]}">真 exit-2 Hard Gate ({len(hg["real"])})</b>: '
               + ", ".join(esc(x) for x in hg["real"]) + "</p>")
    if hg["neutered"]:
        hg_html += (f'<p><b style="color:{COLOR["RED"]}">⚠ NEUTERED ({len(hg["neutered"])})</b>: '
                    + ", ".join(esc(x) for x in hg["neutered"]) + " — wiring '|| exit 0' 吃掉 exit2</p>")
    else:
        hg_html += f'<p style="color:{COLOR["GREEN"]}">✓ 無 neutered gate</p>'
    if hg["soft"]:
        hg_html += '<p><b style="color:%s">intentionally-soft</b>: %s</p>' % (COLOR["YELLOW"], ", ".join(esc(x) for x in hg["soft"]))
    mk = snap["compile_marker"]
    if mk.get("present"):
        hg_html += f'<p style="color:{COLOR["YELLOW"]}">stale compile marker: {mk["files"]} 檔 @ {mk["mtime"]}</p>'
    # recent
    commits = "".join(f'<li><code>{esc(c[:9])}</code> {esc(c[10:])}</li>' for c in arch["commits"])
    pm = ""
    for fn, tail in arch["postmortems"]:
        pm += f'<div class="pm"><b>{esc(fn)}</b><pre>{esc(chr(10).join(tail))}</pre></div>'
    # industry
    ind = "".join(
        f'<div class="box"><b>{esc(axis)}</b> <span class="tag tg-ok">{esc(v)}</span>'
        f'<div class="sub">✓ {esc(cov)}</div><div class="sub sk2">✗ {esc(sk)}</div></div>'
        for axis, v, cov, sk in INDUSTRY)
    # governance (limits & advisories)
    _gc = {"BLOCK": COLOR["RED"], "PROMPT": COLOR["YELLOW"], "RULE": COLOR["BLUE"], "ADVISE": COLOR["GREY"], "LOG": "#aab"}
    gov_rows = "".join(
        f'<tr><td class="ev">{esc(cat)}<br><span class="b" style="background:{_gc.get(typ,"#ccc")};color:#fff">{typ}</span></td>'
        f'<td><div>{esc(mech)}</div><div class="sub">緣由：{esc(origin)}</div>'
        f'<div class="sub" style="color:#1a7f37">{esc(verdict)}</div>'
        f'<div class="sub sk2">改進：{esc(impr)}</div></td></tr>'
        for cat, typ, mech, origin, verdict, impr in GOVERNANCE)

    return f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>InterSubMod Harness Dashboard</title>
<style>
:root{{--fg:#1f2328;--mut:#636c76;--bd:#d8dee4;--bg:#f6f8fa;--card:#fff}}
*{{box-sizing:border-box}}
body{{font-family:-apple-system,"Segoe UI",Roboto,"Droid Sans Fallback",sans-serif;color:var(--fg);background:var(--bg);margin:0;line-height:1.5}}
.wrap{{max-width:1180px;margin:0 auto;padding:20px}}
h1{{font-size:21px;margin:0 0 2px}} h2{{font-size:16px;margin:26px 0 10px;padding-bottom:5px;border-bottom:2px solid var(--bd)}}
.sub{{color:var(--mut);font-size:12px;margin-top:3px}} .sk2{{color:#a04}}
nav{{position:sticky;top:0;background:rgba(246,248,250,.96);backdrop-filter:blur(4px);padding:8px 0;margin:0 -20px 4px;border-bottom:1px solid var(--bd);z-index:9}}
nav a{{font-size:12px;color:var(--mut);text-decoration:none;margin:0 8px;white-space:nowrap}} nav a:hover{{color:var(--fg)}}
.bar{{background:var(--card);border-radius:8px;padding:12px 16px;margin:10px 0;box-shadow:0 1px 2px rgba(0,0,0,.06);border-left:6px solid {oc}}}
.diff{{font-size:12px;color:var(--mut);margin-top:5px}}
.grid{{display:grid;grid-template-columns:1fr 1fr 1fr;gap:10px}} @media(max-width:780px){{.grid{{grid-template-columns:1fr}}}}
.light{{background:var(--card);border-radius:8px;padding:11px 13px;box-shadow:0 1px 2px rgba(0,0,0,.06)}}
.lh{{display:flex;align-items:center;gap:7px;font-size:13px}} .lh .st{{margin-left:auto;font-weight:700;font-size:11px}}
.dot{{width:11px;height:11px;border-radius:50%;display:inline-block}}
.ld{{margin-top:7px;font-size:11.5px;color:#444;word-break:break-word}}
.layer{{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px 12px;margin:8px 0}}
.layer .lt{{font-size:11px;color:var(--mut);text-transform:uppercase;letter-spacing:.04em;margin-bottom:7px}}
.row{{display:flex;flex-wrap:wrap;gap:8px;align-items:stretch}}
.box{{background:#fbfcfd;border:1px solid var(--bd);border-radius:6px;padding:8px 10px;flex:1;min-width:150px}}
.box b{{font-size:12.5px}} .tag{{display:inline-block;font-size:10px;background:#eaeef2;color:#555;border-radius:4px;padding:1px 5px;margin-top:4px}}
.tg-ok{{background:#e6f4ea;color:#1a7f37}}
.phase{{background:#eef3fb;border:1px solid #cfe0f5;border-radius:6px;padding:6px 9px;text-align:center;min-width:64px}}
.phase b{{font-size:13px;color:{COLOR['BLUE']}}} .phase .sk{{font-size:10.5px;color:#333}} .phase .sub{{font-size:9.5px}}
.arr{{align-self:center;color:#9aa4ae;font-weight:700}}
.acol{{background:#fbfcfd;border:1px solid var(--bd);border-radius:6px;padding:8px 10px;flex:1;min-width:200px}}
.ahd{{font-size:11.5px;font-weight:700;color:#333;margin-bottom:6px;border-bottom:1px dashed var(--bd);padding-bottom:3px}}
.chip{{display:inline-block;font-size:11px;background:#eef1f4;border:1px solid #e1e6ea;border-radius:4px;padding:2px 6px;margin:2px}}
.chip.sm{{font-size:10.5px}} .chip.dep{{background:#fbe9e9;color:#a33;text-decoration:line-through}}
.chip.gate{{background:#e6f4ea;border-color:#bfe3c8;color:#1a7f37;font-weight:600}}
.chip.bad{{background:#fdecea;border-color:#f3c0bb;color:#c0392b;font-weight:600}}
.chip.soft{{background:#fdf6e3;border-color:#f0e3b8;color:#9a7d23}}
.b{{font-size:8.5px;border-radius:3px;padding:0 3px;margin-left:3px;vertical-align:middle}}
.b-iso{{background:#dceefc;color:#1d5fa8}} .b-hk{{background:#efe3fb;color:#7a3ec0}} .b-ro{{background:#e6f4ea;color:#1a7f37}}
table{{width:100%;border-collapse:collapse;font-size:12px;background:var(--card);border-radius:8px;overflow:hidden}}
td{{padding:6px 9px;border-bottom:1px solid var(--bd);vertical-align:top}} td.ev{{font-weight:600;white-space:nowrap;color:#333;width:160px}}
details{{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:6px 12px;margin:8px 0}}
summary{{cursor:pointer;font-weight:600;font-size:13px;padding:4px 0}}
.pm{{margin:6px 0}} .pm pre{{background:#0d1117;color:#c9d1d9;font-size:10.5px;padding:8px;border-radius:5px;overflow-x:auto;margin:3px 0}}
code{{background:#eef1f4;padding:1px 4px;border-radius:3px;font-size:11.5px}} ul{{margin:4px 0;padding-left:18px;font-size:12px}}
.foot{{color:var(--mut);font-size:11.5px;margin-top:18px;border-top:1px solid var(--bd);padding-top:10px}}
</style></head><body><div class="wrap">
<h1>InterSubMod HARNESS DASHBOARD</h1>
<div class="sub">snapshot {snap['ts']} · model {esc(snap['settings'].get('model','?'))} · effort {esc(snap['settings'].get('effort','?'))}
 · {snap['disk']['skills']} skills · {snap['disk']['agents']} agents · {snap['disk']['hooks']} hooks · {snap['disk']['workflows']} workflows</div>
<nav><a href="#health">健康燈</a><a href="#arch">整體架構</a><a href="#gate">Hard Gate</a><a href="#comp">各部分明細</a><a href="#issues">最近問題</a><a href="#tech">業界對齊</a><a href="#gov">限制與提示</a><a href="#update">如何更新</a></nav>

<div class="bar"><b>OVERALL</b> {ov['green']} GREEN · {ov['yellow']} YELLOW · {ov['red']} RED
<div class="diff" id="diff">vs last: (computing…)</div></div>

<h2 id="health">① L0 健康燈（8）</h2><div class="grid">{lights_html}</div>

<h2 id="arch">② 整體架構（一眼掌握）</h2>
<div class="layer"><div class="lt">3 入口（context 載入）</div><div class="row">{ep}</div></div>
<div class="layer"><div class="lt">7-Phase Waterfall 研究 cycle 主流程</div><div class="row">{wf}</div></div>
<div class="layer"><div class="lt">假說驗證三層樓 assurance（pre → process → post）</div><div class="row">{asr}</div></div>
<div class="layer"><div class="lt">{snap['disk']['agents']} Agents（依角色；iso=worktree隔離 · RO=read-only verifier · haiku=降智力省token）</div><div class="row">{ag_cols}</div></div>
<div class="layer"><div class="lt">{snap['disk']['hooks']} Hooks × 7 events（<span style="color:#1a7f37">綠=真 exit-2 gate</span> · <span style="color:#c0392b">紅=neutered</span> · 黃=soft）</div><table>{hk_rows}</table></div>
<div class="layer"><div class="lt">State 機 + Evidence 鏈（SoT）</div><div class="row">
<div class="box"><b>state/active.json</b><div class="sub">active cycles ≤5（skill 寫，不可手改）</div></div>
<div class="box"><b>evidence_ledger.jsonl</b><div class="sub">L1-L5 + ⭐tier provenance（append-only SoT）</div></div>
<div class="box"><b>hypothesis_queue.json</b><div class="sub">假說佇列（inject/pivot 寫）</div></div>
<div class="box"><b>{snap['disk']['workflows']} Dynamic Workflows</b><div class="sub">cross_sample_benchmark · harness_audit_2026</div></div>
</div></div>

<h2 id="gate">③ Hard Gate 真實性</h2><div class="layer">{hg_html}</div>

<h2 id="comp">④ 各部分明細</h2>
<details open><summary>{snap['disk']['skills']} Skills（依分類）</summary><div class="row" style="margin-top:8px">{sk_html}</div></details>
<details><summary>{snap['disk']['agents']} Agents（tools / model / isolation）</summary>
<table><tr><td class="ev">agent</td><td>role · tools · model · iso</td></tr>
{''.join(f'<tr><td class="ev">{esc(a["name"])}{" 🚫" if a["deprecated"] else ""}</td><td><span class="sub">{esc(a["role"])}</span><br><code>{esc(a["tools"])}</code> · {esc(a["model"])} · {esc(a["iso"])}</td></tr>' for a in arch["agents"])}
</table></details>

<h2 id="issues">⑤ 最近問題 / 修復</h2>
<div class="layer"><div class="lt">最近 commits</div><ul>{commits}</ul></div>
<details><summary>Postmortem 尾端（最新 entry）</summary>{pm or '<p class="sub">無</p>'}</details>

<h2 id="tech">⑥ 新技術 / 業界對齊（restraint 裁決）</h2>
<div class="layer"><div class="lt">頂尖模式 → 你的覆蓋（0 新外部工具；完整 19 條見 audit 報告）</div><div class="row">{ind}</div></div>

<h2 id="gov">⑦ 限制與提示治理（盤點 · 緣由 · 合理性 · 改進）</h2>
<div class="layer"><div class="lt">每條皆 traced 到事故或治理原則 · restraint-consistent · <span style="color:#c0392b">BLOCK</span>=擋 · <span style="color:#9a7d23">PROMPT</span>=人工確認 · <span style="color:{COLOR['BLUE']}">RULE</span>=規範 · ADVISE=軟提示 · LOG=被動</div>
<table>{gov_rows}</table>
<p class="sub" style="margin-top:8px">主要可改進主題：軟提示廣關鍵字噪音收窄（皆非阻擋，低優先）；最大真破口（neutered Hard Gate）已於 2026-05-31 修復。詳見 <code>InterSubMod/docs/references/migration/20260601_limits_advisories_governance_review_01.md</code></p></div>

<h2 id="update">⑧ 如何持續更新</h2>
<div class="layer">
<p>改了 skill/agent/hook/settings 後，重跑：<code>python3 scripts/harness_health.py</code> → 本檔自動更新，頂部「vs last」顯示哪些燈變了。</p>
<p class="sub">read-only · 無 agent · 無 web · localStorage 存上次 snapshot。完整 audit + 設計：<code>InterSubMod/docs/references/migration/20260531_harness_audit_dashboard_design_01.md</code></p>
</div>
<div class="foot">InterSubMod harness self-audit dashboard · /harness-health · snapshot 留存 state/health_snapshots/</div>
</div>
<script>
const SNAP={data};
try{{
 const prev=JSON.parse(localStorage.getItem('ism_harness_health')||'null');
 const el=document.getElementById('diff');
 if(prev&&prev.ts!==SNAP.ts){{
  const pm={{}};(prev.lights||[]).forEach(l=>pm[l.name]=l.status);
  const ch=(SNAP.lights||[]).filter(l=>pm[l.name]&&pm[l.name]!==l.status).map(l=>`${{l.name}}: ${{pm[l.name]}}→${{l.status}}`);
  const dd=(prev.disk&&SNAP.disk)?['skills','agents','hooks','workflows'].filter(k=>prev.disk[k]!==SNAP.disk[k]).map(k=>`${{k}} ${{prev.disk[k]}}→${{SNAP.disk[k]}}`):[];
  el.textContent='vs last ('+prev.ts+'): '+([...ch,...dd].length?[...ch,...dd].join(' · '):'無變化');
 }}else{{el.textContent='vs last: (首次 baseline 已存)';}}
 localStorage.setItem('ism_harness_health',JSON.stringify(SNAP));
}}catch(e){{}}
</script></body></html>"""


def main():
    os.makedirs(SNAP_DIR, exist_ok=True)
    cmd = _read(CLAUDE_MD)
    raw = _read(SETTINGS)
    try:
        settings = json.loads(raw)
    except Exception:
        settings = {}
    focus = _read(CURRENT_FOCUS)
    try:
        active = json.loads(_read(ACTIVE))
    except Exception:
        active = {}
    try:
        queue = json.loads(_read(QUEUE))
    except Exception:
        queue = {}

    disk = disk_counts()
    l1 = light_count(disk, cmd)
    l2 = light_hardgate(settings, cmd)
    l3 = light_tiergate(raw)
    l4 = light_state(active, focus)
    l5 = light_ledger(focus)
    l6 = light_queue(queue)
    l7 = light_memory_drift()
    l8 = light_doc_paths(cmd)
    l9 = light_hook_wiring(settings, cmd)
    lights = [
        {"name": "COUNT DRIFT", "status": l1[0], "rows": l1[1]},
        {"name": "HARD-GATE TRUTH", "status": l2[0], "rows": l2[1]},
        {"name": "TIER-GATE WIRED", "status": l3[0], "rows": l3[1]},
        {"name": "STATE<->DOC SYNC", "status": l4[0], "rows": l4[1]},
        {"name": "LEDGER FRESH", "status": l5[0], "rows": l5[1]},
        {"name": "QUEUE HYGIENE", "status": l6[0], "rows": l6[1]},
        {"name": "MEMORY DRIFT", "status": l7[0], "rows": l7[1]},
        {"name": "DOC PATH CURRENCY", "status": l8[0], "rows": l8[1]},
        {"name": "HOOK WIRING", "status": l9[0], "rows": l9[1]},
    ]
    mk = compile_marker_status()
    if mk.get("present"):
        lights[1]["rows"].append(f"stale compile marker: {mk['files']} 檔 @ {mk['mtime']}")

    # hard gate detail
    hg_real, hg_neut, hg_soft = [], [], []
    for ev, blocks in (settings.get("hooks", {}) or {}).items():
        for b in blocks:
            for h in b.get("hooks", []):
                command = h.get("command", "")
                m = re.search(r"(/[^\s]+\.sh)", command)
                if not m:
                    continue
                body = _read(m.group(1))
                e2 = ("exit 2" in body) or ("gate_exit" in body)
                masked = bool(re.search(r"\|\|\s*exit 0", command) or re.search(r"\|\|\s*true", command))
                soft = bool(re.search(r"SOFT gate|advisory only|advisory\)", body))
                base = os.path.basename(m.group(1))
                if e2 and not masked:
                    hg_real.append(base)
                elif e2 and masked and soft:
                    hg_soft.append(base)
                elif e2 and masked:
                    hg_neut.append(base)

    arch = {
        "agents": parse_agents(),
        "skills": parse_skill_categories(cmd),
        "hooks": parse_hooks_by_event(settings),
        "commits": recent_commits(6),
        "postmortems": postmortem_tails(2),
    }
    sts = {s: sum(1 for L in lights if L["status"] == s) for s in ("GREEN", "YELLOW", "RED")}
    snap = {
        "ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "disk": disk,
        "settings": {"model": settings.get("model", "?"), "effort": settings.get("effortLevel", "?")},
        "lights": lights,
        "overall": {"green": sts["GREEN"], "yellow": sts["YELLOW"], "red": sts["RED"]},
        "hardgate": {"real": sorted(set(hg_real)), "neutered": sorted(set(hg_neut)), "soft": sorted(set(hg_soft))},
        "compile_marker": mk,
    }
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M")
    with open(os.path.join(SNAP_DIR, f"{stamp}.json"), "w", encoding="utf-8") as f:
        json.dump(snap, f, ensure_ascii=False, indent=1)
    html_path = os.path.join(SNAP_DIR, "dashboard.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(render(snap, arch))

    sym = "●"
    print(f"# InterSubMod Harness Dashboard — {snap['ts']}")
    print(f"  {disk['skills']} skills / {disk['agents']} agents / {disk['hooks']} hooks / {disk['workflows']} workflows")
    for i, L in enumerate(lights, 1):
        print(f"  {sym} {i} {L['name']:<18} [{L['status']}]")
    print(f"  OVERALL: {sts['GREEN']} GREEN · {sts['YELLOW']} YELLOW · {sts['RED']} RED")
    print(f"  dashboard → {html_path}")


if __name__ == "__main__":
    main()
