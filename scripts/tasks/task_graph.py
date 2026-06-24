#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
task_graph.py — research task DAG/WBS layer driver (state/tasks/graph*.json).

A THIN, grep-able, file-based task layer ABOVE cycle/ledger/memory. A *map for
human-gated dispatch* — shows which tasks can run / are blocked; never auto-executes
(full-auto research loops are a project tombstone).

Schema v2 = flat component grouping (state/tasks/graph.json, 30 nodes).
Schema v3 = nested WBS tree (`parent`) + `focus` (the task currently being solved);
            UI auto-expands the path to focus and collapses other branches.
Both supported; the script branches on is_v3(). Two relations coexist:
  parent     = containment (X is a subtask of Y)        → the WBS tree / fold-expand
  depends_on = sequence (Y can't start until X is done) → ready / blocked / flowchart edges
CWL-style io.inputs[{name,required,ref}] marks required (solid) vs optional (dashed).

SoT = the graph JSON. TASKS*.md + *board*.html are AUTO-GENERATED mirrors (never hand-edit).
HTML reuses .claude/skills/verify-workstation/tools/build_workstation.py (§13-A refuse-on-missing
+ per-node human verdict + export). Flowchart/tree = pure hand-coded SVG/HTML + tiny inline JS
(focus switching) — zero external dependency, self-contained, grep-able. No external jsonschema dep.

Usage: python3 task_graph.py [--graph PATH] {validate|check|ready|render|render-html|all}
"""
import json, os, sys, argparse, datetime, subprocess, importlib.util, html, re

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
LEDGER = os.path.join(ROOT, "research/autoresearch/evidence_ledger.jsonl")
ACTIVE = os.path.join(ROOT, "state/active.json")
CYCLES_DIR = os.path.join(ROOT, "state/cycles")
ARCHIVED_DIR = os.path.join(ROOT, "state/cycles_archived")
BUILD_WS = os.path.join(ROOT, ".claude/skills/verify-workstation/tools/build_workstation.py")
MEM_DIR = "/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory"

# set by main() from --graph
GRAPH = os.path.join(ROOT, "state/tasks/graph.json")
TASKS_MD = os.path.join(ROOT, "state/tasks/TASKS.md")
BOARD_HTML = os.path.join(ROOT, "state/tasks/tasks_board.html")

VALID_STATUS = ["todo", "in_progress", "blocked", "done", "dropped"]
VALID_KIND = ["data", "compute", "analysis", "writing", "milestone"]
VALID_OWNER = ["claude", "user", "claude+user", None]
STATUS_LABEL = {"todo": "待辦", "in_progress": "進行中", "blocked": "阻塞", "done": "已完成", "dropped": "已放棄"}
STATUS_COLOR = {"todo": "#94a3b8", "in_progress": "#38bdf8", "blocked": "#ef4444", "done": "#22c55e", "dropped": "#64748b"}
STATUS_BADGE = {"todo": "mut", "in_progress": "info", "blocked": "no", "done": "ok", "dropped": "mut"}
STATUS_MARK = {"todo": "☐", "in_progress": "◐", "blocked": "⛔", "done": "✅", "dropped": "✗"}
KIND_LABEL = {"data": "資料", "compute": "程式", "analysis": "分析", "writing": "撰寫", "milestone": "里程碑"}
ID_RE = re.compile(r"^T-[A-Z0-9-]+$")


# ----------------------------------------------------------------- io / git
def load_graph():
    with open(GRAPH, encoding="utf-8") as fh:
        return json.load(fh)


def out_paths(graph_path):
    d = os.path.dirname(graph_path)
    base = os.path.splitext(os.path.basename(graph_path))[0]
    if base == "graph":
        return os.path.join(d, "TASKS.md"), os.path.join(d, "tasks_board.html")
    suf = base[len("graph_"):] if base.startswith("graph_") else base
    return os.path.join(d, f"TASKS_{suf}.md"), os.path.join(d, f"board_{suf}.html")


def git_info():
    def g(args):
        try:
            return subprocess.check_output(["git"] + args, cwd=ROOT, stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return "?"
    return g(["rev-parse", "--short", "HEAD"]), g(["rev-parse", "--abbrev-ref", "HEAD"])


def git_panel():
    """Snapshot of git state at render time: current branch, last 5 commits, worktrees (each AI's branch)."""
    def g(args):
        try:
            return subprocess.check_output(["git"] + args, cwd=ROOT, stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return ""
    return (g(["rev-parse", "--abbrev-ref", "HEAD"]) or "?",
            g(["log", "--oneline", "-5"]),
            g(["worktree", "list"]))


def known_cycle_ids():
    ids = set()
    try:
        a = json.load(open(ACTIVE, encoding="utf-8"))
        for grp in ("cycles", "recently_concluded"):
            for c in a.get(grp, []):
                if c.get("cycle_id"):
                    ids.add(c["cycle_id"])
    except Exception:
        pass
    for d in (CYCLES_DIR, ARCHIVED_DIR):
        if os.path.isdir(d):
            ids.update(os.listdir(d))
    try:
        for ln in open(LEDGER, encoding="utf-8"):
            ln = ln.strip()
            if ln:
                try:
                    cid = json.loads(ln).get("cycle_id")
                    if cid:
                        ids.add(cid)
                except Exception:
                    continue
    except Exception:
        pass
    return ids


def stale_active_cycle_ids(stale_days=7):
    out = {}
    try:
        a = json.load(open(ACTIVE, encoding="utf-8"))
        today = datetime.date.today()
        for c in a.get("cycles", []):
            la = (c.get("last_advanced_at") or c.get("started_at") or "")[:10]
            try:
                age = (today - datetime.date.fromisoformat(la)).days
            except Exception:
                age = None
            if age is not None and age > stale_days:
                out[c["cycle_id"]] = age
    except Exception:
        pass
    return out


def is_v3(data):
    return bool(data.get("focus")) or data.get("schema_version") == "3.0" or any(n.get("parent") for n in data.get("nodes", []))


def trunc(s, n=12):
    s = s or ""
    return s if len(s) <= n else s[:n] + "…"


# ----------------------------------------------------------------- validate
def validate(data):
    nodes = data.get("nodes", [])
    comps = {c["id"] for c in data.get("components", [])}
    errors = []
    idset = {n.get("id") for n in nodes}
    seen = set()
    for n in nodes:
        nid = n.get("id", "<no-id>")
        if not n.get("id") or not ID_RE.match(n.get("id", "")):
            errors.append(f"bad id pattern: {nid}")
        if nid in seen:
            errors.append(f"duplicate node id: {nid}")
        seen.add(nid)
        if not n.get("title"):
            errors.append(f"{nid}: missing title")
        if comps and n.get("component") and n["component"] not in comps:
            errors.append(f"{nid}: component {n.get('component')!r} not in components[]")
        if n.get("kind") not in VALID_KIND:
            errors.append(f"{nid}: bad kind {n.get('kind')!r}")
        if n.get("status") not in VALID_STATUS:
            errors.append(f"{nid}: bad status {n.get('status')!r}")
        if n.get("owner") not in VALID_OWNER:
            errors.append(f"{nid}: bad owner {n.get('owner')!r}")
        p = n.get("parent")
        if p is not None and p not in idset:
            errors.append(f"{nid}: parent {p!r} not a known node")
        for dep in n.get("depends_on", []):
            if dep not in idset:
                errors.append(f"{nid}: depends_on missing node {dep!r}")
        io = n.get("io")
        # compute 節點必須宣告非空 io.inputs+outputs（2026-06-25 架構稽核 D1-2：從 check() advisory 升 validate error）
        if n.get("kind") == "compute" and not (io and io.get("inputs") and io.get("outputs")):
            errors.append(f"{nid}: compute 節點必須宣告非空 io.inputs+outputs")
        if io is not None:
            if "inputs" not in io or "outputs" not in io:
                errors.append(f"{nid}: io must have inputs+outputs")
            deps = set(n.get("depends_on", []))
            for inp in io.get("inputs", []):
                ref = inp.get("ref")
                if ref is not None and ref not in idset:
                    errors.append(f"{nid}: io.input ref {ref!r} not a known node")
                if ref is not None and inp.get("required") and ref not in deps:
                    errors.append(f"{nid}: required io.input ref {ref!r} not in depends_on")
        prob = n.get("problem")
        if prob is not None and not isinstance(prob, dict):
            errors.append(f"{nid}: problem must be object")
        if isinstance(prob, dict) and prob.get("fix_child") and prob["fix_child"] not in idset:
            errors.append(f"{nid}: problem.fix_child {prob['fix_child']!r} not a known node")
        if n.get("verify") is not None and not isinstance(n.get("verify"), dict):
            errors.append(f"{nid}: verify must be object")
        if n.get("observe") is not None and not isinstance(n.get("observe"), list):
            errors.append(f"{nid}: observe must be array")
    # depends_on cycle
    adj = {n["id"]: [d for d in n.get("depends_on", []) if d in idset] for n in nodes if n.get("id")}
    _cycle_check(adj, errors, "dependency")
    # parent cycle
    padj = {n["id"]: ([n["parent"]] if n.get("parent") in idset else []) for n in nodes if n.get("id")}
    _cycle_check(padj, errors, "parent")
    return errors


def _cycle_check(adj, errors, label):
    color = {}

    def dfs(u, stack):
        color[u] = 1
        for v in adj.get(u, []):
            if color.get(v) == 1:
                errors.append(f"{label} cycle: {' -> '.join(stack + [u, v])}")
                return
            if color.get(v, 0) == 0:
                dfs(v, stack + [u])
        color[u] = 2

    for nid in adj:
        if color.get(nid, 0) == 0:
            dfs(nid, [])


# ----------------------------------------------------------------- check + detail gaps
def detail_gaps(nodes):
    gaps = []
    for n in nodes:
        io = n.get("io")
        if io is not None:
            names = [i.get("name") for i in io.get("inputs", [])] + list(io.get("outputs", []))
            if "TBD" in names or any((not x) for x in names):
                gaps.append((n["id"], "I/O 含 TBD/空白"))
        if n.get("kind") == "compute" and not io:
            gaps.append((n["id"], "compute 任務缺 I/O 規格"))
        if not n.get("owner"):
            gaps.append((n["id"], "未指派 owner"))
        # 缺說明：milestone 容器不需；headline/summary/notes 任一即算有說明
        has_desc = n.get("notes") or n.get("headline") or n.get("summary")
        if n.get("kind") != "milestone" and not has_desc and not (n.get("links") or {}):
            gaps.append((n["id"], "缺說明(headline/notes)與連結(links)"))
    return gaps


def check(data):
    nodes = data["nodes"]
    findings = []
    kc, stale = known_cycle_ids(), stale_active_cycle_ids()
    for n in nodes:
        nid, links = n.get("id"), (n.get("links") or {})
        has_prov = any([links.get("cycle_id"), links.get("ledger_entries"), links.get("reports"), links.get("memory")])
        if n.get("status") == "done" and n.get("kind") != "milestone" and not has_prov:
            findings.append(("WARN", nid, "status=done 但無 links 證據"))
        cid = links.get("cycle_id")
        if cid and cid not in kc:
            findings.append(("WARN", nid, f"links.cycle_id={cid} 找不到 (dangling/archived)"))
        if cid and cid in stale:
            findings.append(("DRIFT", nid, f"指向的 cycle {cid} 在 active.json 已 stale {stale[cid]} 天"))
    # F1: branch collision (two in_progress tasks on the same branch = risk of clobber)
    for br, tasks in branch_task_map(nodes).items():
        ips = [t for t in tasks if t.get("status") == "in_progress"]
        if len(ips) >= 2:
            findings.append(("WARN", br, f"分支碰撞：{', '.join(t['id'] for t in ips)} 同時 in_progress 同分支（各開 worktree）"))
    # F4/F5: verify failed but no fix task scaffolded
    for n in nodes:
        if (n.get("verify") or {}).get("result") == "fail" and not (n.get("problem") or {}).get("fix_child"):
            findings.append(("WARN", n["id"], "verify=fail 但無修任務（跑 new-problem 建子修任務）"))
    for nid, why in detail_gaps(nodes):
        findings.append(("DETAIL", nid, why))
    return findings, stale


# ----------------------------------------------------------------- graph helpers (depends_on)
def ready(nodes):
    by = {n["id"]: n for n in nodes}
    return [n for n in nodes if n.get("status") == "todo" and n.get("kind") != "milestone"
            and all(by.get(d, {}).get("status") == "done" for d in n.get("depends_on", []))]


def blocked_view(nodes):
    by = {n["id"]: n for n in nodes}
    out = []
    for n in nodes:
        if n.get("status") not in ("blocked", "todo"):
            continue
        unmet = [d for d in n.get("depends_on", []) if by.get(d, {}).get("status") != "done"]
        if n.get("status") == "blocked" or unmet or n.get("missing_info"):
            out.append((n, unmet))
    return out


def levels(nodes, subset_ids=None):
    by = {n["id"]: n for n in nodes}
    allow = subset_ids if subset_ids is not None else set(by)
    memo = {}

    def lvl(nid, seen):
        if nid in memo:
            return memo[nid]
        if nid in seen:
            return 0
        deps = [d for d in by.get(nid, {}).get("depends_on", []) if d in by and d in allow]
        v = 0 if not deps else 1 + max(lvl(d, seen | {nid}) for d in deps)
        memo[nid] = v
        return v

    return {nid: lvl(nid, set()) for nid in (by if subset_ids is None else subset_ids)}


def degrees(nodes):
    indeg = {n["id"]: len(n.get("depends_on", [])) for n in nodes}
    outdeg = {n["id"]: 0 for n in nodes}
    for n in nodes:
        for d in n.get("depends_on", []):
            if d in outdeg:
                outdeg[d] += 1
    return indeg, outdeg


def edge_required(node, dep):
    for inp in (node.get("io") or {}).get("inputs", []):
        if inp.get("ref") == dep:
            return bool(inp.get("required", True))
    return True


def input_readiness(node, by):
    """Per-input state: ready (ref task done) / wait (ref task not done) / external (ref null)."""
    out = []
    for inp in (node.get("io") or {}).get("inputs", []):
        ref = inp.get("ref")
        if ref is None:
            state = "external"
        elif by.get(ref, {}).get("status") == "done":
            state = "ready"
        else:
            state = "wait"
        out.append({"name": inp.get("name", ""), "required": bool(inp.get("required", True)), "ref": ref, "state": state})
    return out


def impact(node_id, nodes):
    """Who completing this task unlocks = nodes that depend_on it."""
    return [n["id"] for n in nodes if node_id in n.get("depends_on", [])]


READY_MARK = {"ready": "✓就緒", "wait": "⧖待產出", "external": "⊘需外部備"}


# ----------------------------------------------------------------- parent-tree helpers (v3)
def tree_maps(nodes):
    by = {n["id"]: n for n in nodes}
    kids, roots = {}, []
    for n in nodes:
        p = n.get("parent")
        if p:
            kids.setdefault(p, []).append(n["id"])
        else:
            roots.append(n["id"])
    for k in kids:
        kids[k] = sorted(kids[k])
    return by, kids, sorted(roots)


def ancestors(nid, by):
    out, cur = [], by.get(nid, {}).get("parent")
    while cur:
        out.append(cur)
        cur = by.get(cur, {}).get("parent")
    return list(reversed(out))


def subtree(nid, kids):
    acc, stack = {nid}, [nid]
    while stack:
        c = stack.pop()
        for ch in kids.get(c, []):
            if ch not in acc:
                acc.add(ch)
                stack.append(ch)
    return acc


# ----------------------------------------------------------------- SVG dependency DAG (all)
def svg_dag(nodes):
    by = {n["id"]: n for n in nodes}
    lv = levels(nodes)
    cols = {}
    for n in sorted(nodes, key=lambda x: x["id"]):
        cols.setdefault(lv[n["id"]], []).append(n["id"])
    maxlevel = max(lv.values()) if lv else 0
    maxcount = max((len(c) for c in cols.values()), default=1)
    W, H, GX, GY = 200, 50, 250, 66
    pos = {}
    for level, members in cols.items():
        for i, nid in enumerate(members):
            pos[nid] = (24 + level * GX, 36 + i * GY)
    width, height = 48 + (maxlevel + 1) * GX, 56 + maxcount * GY
    edges = []
    for n in nodes:
        x2, y2 = pos[n["id"]]
        for d in n.get("depends_on", []):
            if d not in pos:
                continue
            x1, y1 = pos[d]
            sx, sy, ex, ey = x1 + W, y1 + H / 2, x2, y2 + H / 2
            mx = (sx + ex) / 2
            col = STATUS_COLOR.get(by.get(d, {}).get("status", "todo"), "#94a3b8")
            dash = "" if edge_required(n, d) else ' stroke-dasharray="5 4"'
            edges.append(f'<path d="M{sx},{sy} C{mx},{sy} {mx},{ey} {ex-6},{ey}" fill="none" stroke="{col}" stroke-width="1.7" opacity="0.85"{dash} marker-end="url(#arr)"/>')
    boxes = []
    for n in nodes:
        x, y = pos[n["id"]]
        st = n.get("status", "todo")
        col = STATUS_COLOR.get(st, "#94a3b8")
        bd = ' stroke-dasharray="4 3"' if st == "dropped" else ""
        boxes.append(
            f'<g style="cursor:pointer" onclick="openM(\'{n["id"]}\')"><rect x="{x}" y="{y}" rx="8" width="{W}" height="{H}" fill="#1e293b" stroke="{col}" stroke-width="2.2"{bd}/>'
            f'<circle cx="{x+13}" cy="{y+15}" r="5" fill="{col}"/>'
            f'<text x="{x+25}" y="{y+19}" fill="#e2e8f0" font-size="12" font-weight="700">{html.escape(n["id"])}</text>'
            f'<text x="{x+W-10}" y="{y+19}" fill="#94a3b8" font-size="10" text-anchor="end">{html.escape(KIND_LABEL.get(n.get("kind"),""))}·{html.escape(STATUS_LABEL.get(st,st))}</text>'
            f'<text x="{x+10}" y="{y+38}" fill="#cbd5e1" font-size="10.5">{html.escape(trunc(n["title"],14))}</text></g>')
    return (f'<svg viewBox="0 0 {width} {height}" style="width:100%;min-width:700px;font-family:\'Noto Sans CJK TC\',system-ui,sans-serif">'
            f'<defs><marker id="arr" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#64748b"/></marker></defs>'
            + "".join(edges) + "".join(boxes) + "</svg>"
            + '<div class="sub" style="margin-top:6px">點節點看細節 ｜ A→B = A 須先完成才能跑 B（<b>實線=必需</b>／<b>虛線=可選</b>；邊色=上游狀態）</div>')


# ----------------------------------------------------------------- SVG program flowchart (data+compute)
def svg_pipeline(nodes):
    prog = [n for n in nodes if n.get("kind") in ("data", "compute")]
    if not prog:
        return '<div class="sub">此範圍無 data/compute（程式/資料）節點</div>'
    ids = {n["id"] for n in prog}
    lv = levels(nodes, subset_ids=ids)
    cols = {}
    for n in sorted(prog, key=lambda x: x["id"]):
        cols.setdefault(lv[n["id"]], []).append(n)
    maxlevel = max(lv.values()) if lv else 0
    maxcount = max((len(c) for c in cols.values()), default=1)
    W, H, GX, GY = 234, 116, 304, 146
    pos = {}
    for level, members in cols.items():
        for i, n in enumerate(members):
            pos[n["id"]] = (28 + level * GX, 30 + i * GY)
    width, height = 56 + (maxlevel + 1) * GX, 40 + maxcount * GY
    by = {n["id"]: n for n in nodes}
    edges, boxes = [], []
    for n in prog:
        x2, y2 = pos[n["id"]]
        for d in n.get("depends_on", []):
            if d not in pos:
                continue
            x1, y1 = pos[d]
            sx, sy, ex, ey = x1 + W, y1 + H / 2, x2, y2 + H / 2
            mx = (sx + ex) / 2
            col = STATUS_COLOR.get(by.get(d, {}).get("status", "todo"), "#94a3b8")
            dash = "" if edge_required(n, d) else ' stroke-dasharray="5 4"'
            edges.append(f'<path d="M{sx},{sy} C{mx},{sy} {mx},{ey} {ex-6},{ey}" fill="none" stroke="{col}" stroke-width="1.8" opacity="0.85"{dash} marker-end="url(#parr)"/>')
    for n in prog:
        x, y = pos[n["id"]]
        st = n.get("status", "todo")
        col = STATUS_COLOR.get(st, "#94a3b8")
        io = n.get("io") or {}
        ins = io.get("inputs", [])
        in_lines = []
        for inp in ins[:3]:
            mark = "●" if inp.get("required") else "○"
            ext = "" if inp.get("ref") else " (外部)"
            in_lines.append(f'{mark} {trunc(inp.get("name",""),18)}{ext}')
        if len(ins) > 3:
            in_lines.append(f'… +{len(ins)-3}')
        outs = io.get("outputs", [])
        out_line = "▶ " + "；".join(trunc(o, 16) for o in outs[:2]) + (f" +{len(outs)-2}" if len(outs) > 2 else "")
        shape = (f'<rect x="{x}" y="{y}" rx="6" width="{W}" height="{H}" fill="#11202f" stroke="{col}" stroke-width="2.4"/>'
                 if n.get("kind") == "compute"
                 else f'<rect x="{x}" y="{y}" rx="18" width="{W}" height="{H}" fill="#0b1a14" stroke="{col}" stroke-width="2.4"/>')
        txt = (f'<text x="{x+12}" y="{y+19}" fill="#e2e8f0" font-size="12" font-weight="700">{html.escape(n["id"])} '
               f'<tspan fill="{col}" font-size="10">{html.escape(KIND_LABEL.get(n.get("kind"),""))}·{html.escape(STATUS_LABEL.get(st,st))}</tspan></text>'
               f'<text x="{x+12}" y="{y+35}" fill="#cbd5e1" font-size="10.5">{html.escape(trunc(n["title"],24))}</text>')
        for j, ln in enumerate(in_lines):
            txt += f'<text x="{x+12}" y="{y+52+j*14}" fill="#93c5fd" font-size="9.5">{html.escape(ln)}</text>'
        txt += f'<text x="{x+12}" y="{y+H-8}" fill="#86efac" font-size="9.5">{html.escape(out_line)}</text>'
        boxes.append(f'<g style="cursor:pointer" onclick="openM(\'{n["id"]}\')">{shape}{txt}</g>')
    return (f'<svg viewBox="0 0 {width} {height}" style="width:100%;min-width:760px;font-family:\'Noto Sans CJK TC\',system-ui,sans-serif">'
            f'<defs><marker id="parr" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#64748b"/></marker></defs>'
            + "".join(edges) + "".join(boxes) + "</svg>"
            + '<div class="sub" style="margin-top:6px">只含<b>程式/資料</b>節點（圓角=資料、方角=程式）｜ in: <b>●必需 ○可選</b>（標「外部」=原始輸入）｜ ▶ out: 產物 ｜ 實線/虛線=必需/可選上游</div>')


def node_mini_svg(n, indeg, outdeg):
    st = n.get("status", "todo")
    col = STATUS_COLOR.get(st, "#94a3b8")
    io = n.get("io") or {}
    req = sum(1 for i in io.get("inputs", []) if i.get("required"))
    opt = len(io.get("inputs", [])) - req
    return (f'<svg viewBox="0 0 320 50" style="width:100%">'
            f'<rect x="2" y="2" rx="6" width="316" height="20" fill="{col}" opacity="0.22"/>'
            f'<rect x="2" y="2" rx="6" width="316" height="20" fill="none" stroke="{col}" stroke-width="1.5"/>'
            f'<text x="10" y="16" fill="{col}" font-size="12" font-weight="700">{html.escape(KIND_LABEL.get(n.get("kind"),""))}·{html.escape(STATUS_LABEL.get(st,st))}</text>'
            f'<text x="312" y="16" fill="#94a3b8" font-size="11" text-anchor="end">{html.escape(n.get("parent") or n.get("component") or "")}</text>'
            f'<text x="10" y="42" fill="#cbd5e1" font-size="11">in ●{req} ○{opt} ｜ out {len(io.get("outputs",[]))} ｜ 被依賴 {outdeg} ｜ 缺 {len(n.get("missing_info") or [])}</text></svg>')


# ----------------------------------------------------------------- v3 focus hierarchy
def _has_flow(nid, by, kids):
    return any(by[x].get("kind") in ("data", "compute") for x in subtree(nid, kids))


def _tree_left(data):
    """Left column: breadcrumb + focus select + nav tree (click any node → focusFlow updates right panel)."""
    nodes = data["nodes"]
    by, kids, roots = tree_maps(nodes)
    focus = data.get("focus")
    open_set = (set(ancestors(focus, by)) | subtree(focus, kids)) if focus else set()

    def sub_counts(nid):
        desc = [x for x in subtree(nid, kids) if x != nid]
        return sum(1 for x in desc if by[x].get("status") == "done"), len(desc)

    def render(nid):
        n = by[nid]
        st = n.get("status", "todo")
        col = STATUS_COLOR.get(st, "#94a3b8")
        ch = kids.get(nid, [])
        hl = ' style="box-shadow:0 0 0 2px #f2c94c"' if nid == focus else ''
        probmark = ' <span class="badge b-no" title="有問題">⚠</span>' if n.get("problem") else ''
        if ch:
            done, tot = sub_counts(nid)
            opened = " open" if nid in open_set else ""
            star = "★ " if nid == focus else ""
            summary = (f'<summary class="tgsum" style="cursor:pointer">{STATUS_MARK.get(st,"·")} '
                       f'<span class="tgfocus" onclick="event.preventDefault();focusFlow(\'{nid}\')"><b style="color:{col}">{star}{html.escape(n["id"])}</b> {html.escape(trunc(n["title"],20))}</span> '
                       f'<span class="badge b-info">{done}/{tot}</span>{probmark}</summary>')
            inner = "".join(render(c) for c in ch)
            return (f'<details id="tgn-{nid}" data-tid="{nid}"{opened}{hl} style="margin:3px 0 3px 10px;border-left:2px solid #334155;padding:2px 0 2px 8px">'
                    f'{summary}<div>{inner}</div></details>')
        rdy_s = "".join(("✓" if i["state"] == "ready" else "⧖" if i["state"] == "wait" else "⊘") for i in input_readiness(n, by)) or "—"
        return (f'<div id="tgn-{nid}" data-tid="{nid}"{hl} style="margin:2px 0 2px 22px;cursor:pointer;padding:2px 6px;border-left:2px solid #334155" '
                f'onclick="focusFlow(\'{nid}\')"><span style="color:{col}">{STATUS_MARK.get(st,"·")} {html.escape(n["id"])} {html.escape(trunc(n["title"],18))}</span>{probmark} '
                f'<span class="sub">{rdy_s}</span></div>')

    tree = "".join(render(r) for r in roots)
    bc_ids = (ancestors(focus, by) + [focus]) if focus else []
    bc = " ＞ ".join(f'<a style="cursor:pointer;color:{"#f2c94c" if i==focus else "#60a5fa"}" onclick="focusFlow(\'{i}\')">{html.escape(by[i]["title"])}</a>' for i in bc_ids)
    opts = "".join(f'<option value="{n["id"]}"{" selected" if n["id"]==focus else ""}>{"　"*len(ancestors(n["id"],by))}{html.escape(n["id"])} {html.escape(trunc(n["title"],20))}</option>' for n in nodes)
    return (f'<div id="tg-bc" style="margin:0 0 6px;font-size:13px">{bc}</div>'
            f'<div class="row" style="margin-bottom:6px"><span class="sub">聚焦：</span><select id="tg-focus" onchange="focusFlow(this.value)">{opts}</select></div>'
            f'<div class="sub" style="margin-bottom:4px">▸ 三角形＝展開/收合；點<b>任務名稱</b> → 右欄顯示「輸入 → 任務 → 後續」＋細節＋文件連結。葉節點末標輸入就緒 ✓/⧖/⊘；⚠=有問題。</div>'
            f'<div>{tree}</div>')


def _abs_href(ref):
    if not ref:
        return None
    r = str(ref).strip()
    if r.startswith("InterSubMod/"):
        return "file://" + os.path.join(ROOT, r[len("InterSubMod/"):])
    if r.endswith(".md") or r.endswith(".html"):
        return "file://" + (r if os.path.isabs(r) else os.path.join(ROOT, r))
    if r.startswith(("project_", "feedback_", "reference_")):
        return "file://" + os.path.join(MEM_DIR, r + ".md")
    if re.match(r"^[0-9]{8}[-_]", r) and os.path.isdir(os.path.join(ROOT, "state/cycles", r)):
        return "file://" + os.path.join(ROOT, "state/cycles", r, "state.json")
    return None


def linkify(ref, label=None):
    lab = html.escape(label if label is not None else str(ref))
    href = _abs_href(ref)
    return f'<a href="{html.escape(href)}" target="_blank" style="color:#60a5fa;text-decoration:underline">{lab}</a>' if href else lab


def neighborhood_html(fid, data, by, kids):
    """Local focus view per user request: needed inputs → this/these task(s) → subsequent tasks only."""
    sub = subtree(fid, kids)
    is_c = bool(kids.get(fid))
    center = [by[c] for c in kids.get(fid, [])] if is_c else [by[fid]]
    seen, inputs = set(), []
    for nn in ([by[fid]] + center if is_c else [by[fid]]):
        for inp in (nn.get("io") or {}).get("inputs", []):
            if inp.get("ref") in sub:
                continue
            key = (inp.get("name"), inp.get("ref"))
            if key not in seen:
                seen.add(key)
                inputs.append(inp)
    subseq = [n for n in data["nodes"] if n["id"] not in sub and any(d in sub for d in n.get("depends_on", []))]

    def in_pill(inp):
        ref = inp.get("ref")
        sc, sl = (("#64748b", "外部") if ref is None
                  else ("#22c55e", "就緒") if by.get(ref, {}).get("status") == "done" else ("#f59e0b", "待產出"))
        return (f'<div class="tg-pill" style="border-color:{sc}"><span style="color:{sc}">{"●" if inp.get("required") else "○"}</span> '
                f'{html.escape(inp.get("name",""))}<span class="tg-pill-s" style="color:{sc}">{sl}</span></div>')

    def task_box(n):
        col = STATUS_COLOR.get(n.get("status"), "#94a3b8")
        cls = " tg-focusbox" if n["id"] == fid else ""
        click = "" if n["id"] == fid else ' onclick="focusFlow(\'' + n["id"] + '\')" style="cursor:pointer;border-color:' + col + '"'
        style = f' style="border-color:{col}"' if n["id"] == fid else ""
        return (f'<div class="tg-taskbox{cls}"{style}{click}>'
                f'<div style="color:{col};font-weight:700;font-size:12px">{STATUS_MARK.get(n.get("status"),"·")} {html.escape(n["id"])}</div>'
                f'<div style="font-size:11px;color:#e2e8f0">{html.escape(trunc(n["title"],24))}</div>'
                f'<div class="sub" style="font-size:9.5px">{KIND_LABEL.get(n.get("kind"),"")}·{STATUS_LABEL.get(n.get("status"),"")}</div></div>')

    def seq_pill(n):
        col = STATUS_COLOR.get(n.get("status"), "#94a3b8")
        return (f'<div class="tg-pill" style="border-color:{col};cursor:pointer" onclick="focusFlow(\'{n["id"]}\')">'
                f'<span style="color:{col}">{STATUS_MARK.get(n.get("status"),"·")}</span> {html.escape(n["id"])} {html.escape(trunc(n["title"],14))}</div>')

    in_h = "".join(in_pill(i) for i in inputs) or '<div class="sub">（無外部輸入）</div>'
    task_h = "".join(task_box(n) for n in center)
    seq_h = "".join(seq_pill(n) for n in subseq) or '<div class="sub">（無後續任務）</div>'
    return ('<div class="tg-nb">'
            f'<div class="tg-nbcol"><div class="tg-nbh">需要的輸入</div>{in_h}</div><div class="tg-arrow">➜</div>'
            f'<div class="tg-nbcol"><div class="tg-nbh">{"區塊任務" if is_c else "此任務"}</div>{task_h}</div><div class="tg-arrow">➜</div>'
            f'<div class="tg-nbcol"><div class="tg-nbh">後續任務</div>{seq_h}</div></div>')


def detail_panel_html(nid, data, by, kids):
    """Right panel, layered: L0 headline + L1 summary (always) + neighborhood + key facts + L2 細節(folded) + verdict."""
    n = by[nid]
    st = n.get("status", "todo")
    col = STATUS_COLOR.get(st, "#94a3b8")
    headline = n.get("headline") or n.get("title") or ""
    summary = n.get("summary") or n.get("notes") or ""
    P = [f'<div class="row" style="justify-content:space-between;border-bottom:2px solid {col};padding-bottom:6px;margin-bottom:8px">'
         f'<b style="font-size:14px;color:{col}">{STATUS_MARK.get(st,"·")} {html.escape(n["id"])} {html.escape(n["title"])}</b>'
         f'<span class="badge b-{STATUS_BADGE.get(st,"mut")}">{KIND_LABEL.get(n.get("kind"),"")}·{STATUS_LABEL.get(st,st)}</span></div>']
    # L0 重點 (always) + L1 摘要 (always)
    P.append(f'<div class="tg-l0">▍{html.escape(headline)}</div>')
    if summary and summary != headline:
        P.append(f'<div class="tg-l1">{html.escape(summary)}</div>')
    goal = n.get("goal")
    goal_v = html.escape(goal) if goal else '<span style="color:#64748b">（待確認填入）</span>'
    P.append(f'<div class="kv"><span class="k">🎯 目標/完成判準</span><span class="v">{goal_v}</span></div>')
    # core visual: neighborhood (always)
    P.append('<div class="tg-h">聚焦流程（輸入 → 任務 → 後續）</div>' + neighborhood_html(nid, data, by, kids))
    # key facts (always = 重點)
    out = (n.get("io") or {}).get("outputs", [])
    if out:
        P.append(f'<div class="kv"><span class="k">產出 / 交接（其他 AI 讀）</span><span class="v">{html.escape("；".join(out))}</span></div>')
    if n.get("missing_info"):
        P.append(f'<div class="kv"><span class="k" style="color:#f59e0b">缺資訊</span><span class="v" style="color:#f59e0b">{html.escape("；".join(n["missing_info"]))}</span></div>')
    lk0 = n.get("links") or {}
    if lk0.get("branch") or lk0.get("worktree"):
        P.append(f'<div class="kv"><span class="k">分支 / 隔離(worktree)</span><span class="v">{html.escape(lk0.get("branch") or "—")} / {html.escape(lk0.get("worktree") or "—")}</span></div>')
    # L2 細節解釋 (folded)
    L2 = []
    if n.get("observe"):
        L2.append(f'<div class="kv"><span class="k">觀察材料</span><span class="v">{html.escape("；".join(n["observe"]))}</span></div>')
    lk = n.get("links") or {}
    rows = []
    if lk.get("cycle_id"):
        rows.append("cycle " + linkify(lk["cycle_id"]))
    for mslug in lk.get("memory", []):
        rows.append("memory " + linkify(mslug))
    for rp in lk.get("reports", []):
        rows.append("report " + linkify(rp))
    if rows:
        L2.append(f'<div class="kv"><span class="k">文件連結</span><span class="v">{" ｜ ".join(rows)}</span></div>')
    hf = n.get("handoff")
    hp = os.path.join(ROOT, "state/tasks/handoff", nid + ".md")
    if hf or os.path.exists(hp):
        stc = (hf or {}).get("status") or _status_to_contract(n.get("status"))
        hlink = (f' ↳<a href="file://{html.escape(hp)}" target="_blank" style="color:#60a5fa">交接包.md</a>'
                 if os.path.exists(hp) else '（跑 handoff 產生）')
        L2.append(f'<div class="kv"><span class="k">AI 交接</span><span class="v">contract <b>{html.escape(str(stc))}</b>{hlink}</span></div>')
    p = n.get("problem")
    if p:
        body = "<b>⚠ 問題：</b>" + html.escape(p.get("summary", ""))
        if p.get("cause"):
            body += "<br><b>根因：</b>" + html.escape(p["cause"])
        if p.get("options"):
            body += "<br><b>選項：</b>" + html.escape("；".join(p["options"]))
        if p.get("ref"):
            body += "<br><b>出處：</b>" + linkify(p["ref"])
        if p.get("fix_child"):
            fc = p["fix_child"]
            body += (f'<br><b>修任務：</b><a href="#" onclick="event.preventDefault();focusFlow(\'{fc}\')" style="color:#60a5fa">'
                     f'{html.escape(fc)} {html.escape(by.get(fc,{}).get("title",""))}</a>')
        L2.append(f'<div class="note">{body}</div>')
    vf = n.get("verify")
    if vf:
        res = vf.get("result")
        rescol = {"pass": "#22c55e", "fail": "#ef4444", "pending": "#f59e0b"}.get(res, "#94a3b8")
        resstr = f' · <b style="color:{rescol}">result={res}</b>' if res else ""
        ev = (" ↳" + linkify(vf.get("evidence"))) if vf.get("evidence") else ""
        L2.append(f'<div class="kv"><span class="k">驗證方式</span><span class="v">{html.escape(vf.get("how",""))} ↳{linkify(vf.get("ref")) if vf.get("ref") else "—"}{resstr}{ev}</span></div>')
    if n.get("on_fail"):
        L2.append(f'<div class="kv"><span class="k">驗證失敗回饋</span><span class="v">{html.escape(n["on_fail"])}</span></div>')
    if n.get("notes") and n.get("notes") != summary:
        L2.append(f'<div class="sub" style="margin-top:6px">{html.escape(n["notes"])}</div>')
    if L2:
        P.append('<details class="tg-l2"><summary style="cursor:pointer;font-size:11px;color:#7dd3fc;margin-top:8px">▸ 細節解釋（觀察 / 問題 / 驗證 / 連結 / notes）</summary>'
                 '<div style="padding:6px 0">' + "".join(L2) + '</div></details>')
    P.append('<div class="jrow" style="margin-top:12px" data-vfor="' + nid + '">'
             + "".join(f'<button class="jb" data-v="{k}" onclick="tgVerdict(\'{nid}\',\'{k}\')">{lab}</button>' for k, lab in [("in", "確認無誤"), ("doubt", "待查"), ("ex", "有問題"), ("un", "清除")])
             + '</div>')
    return f'<div class="tg-detail" id="detail-{nid}" style="display:none">' + "".join(P) + '</div>'


def status_header_html(data):
    """Top of console: current task + progress + completed-in-order + git/branch snapshot."""
    nodes = data["nodes"]
    by = {n["id"]: n for n in nodes}
    focus = data.get("focus")
    cnt = {s: sum(1 for n in nodes if n.get("status") == s) for s in VALID_STATUS}
    lv = levels(nodes)
    done = sorted([n for n in nodes if n.get("status") == "done" and n.get("kind") != "milestone"], key=lambda n: (lv[n["id"]], n["id"]))
    fn = by.get(focus, {})
    fcol = STATUS_COLOR.get(fn.get("status"), "#94a3b8")
    fhead = (fn.get("headline") or fn.get("summary") or fn.get("notes") or "") if fn else ""
    chips = "".join(f'<span class="tg-chip" onclick="focusFlow(\'{n["id"]}\')">✅ {html.escape(n["id"])} {html.escape(trunc(n["title"],9))}</span>' for n in done) or '<span class="sub">（尚無）</span>'
    branch, log, wts = git_panel()
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    loglines = "".join(f'<div class="tg-gitline">{html.escape(l)}</div>' for l in log.splitlines()) or '<div class="sub">（無）</div>'
    wtlines = "".join(f'<div class="tg-gitline">{html.escape(w)}</div>' for w in wts.splitlines()) or '<div class="sub">（單一工作樹）</div>'
    # 分支 ↔ 任務 對照 + 碰撞旗標 (F1)
    bm, wb = branch_task_map(nodes), worktree_branches()
    btrows = ""
    for br in sorted(bm):
        tks = bm[br]
        ips = [t for t in tks if t.get("status") == "in_progress"]
        flag = ' <span style="color:#ef4444">⚠碰撞</span>' if len(ips) >= 2 else ""
        iso = "" if br in wb else ' <span style="color:#f59e0b">⚠未隔離(不在worktree)</span>'
        btrows += f'<div class="tg-gitline">{html.escape(br)} ← {html.escape(", ".join(t["id"] for t in tks))}{flag}{iso}</div>'
    btrows = btrows or '<div class="sub">（尚無任務標分支；用 task_graph.py claim --task T-x --branch ...）</div>'
    return (
        '<div class="tg-status">'
        f'<div style="font-size:14px"><b>📍 現在處理：</b><span style="color:{fcol}">★ {html.escape(focus or "—")} {html.escape(fn.get("title","") if fn else "")} [{STATUS_LABEL.get(fn.get("status"),"")}]</span>'
        f'　｜　進度 <b style="color:#22c55e">✅{cnt["done"]}</b> ◐{cnt["in_progress"]} ⛔{cnt["blocked"]} ☐{cnt["todo"]} ✗{cnt["dropped"]}（共 {len(nodes)}）</div>'
        + (f'<div style="margin-top:3px;font-size:12px;color:#cbd5e1">▍{html.escape(fhead)}</div>' if fhead else "")
        + f'<div style="margin-top:7px"><span class="sub">已完成（依依賴順序，可點跳）：</span> {chips}</div>'
        f'<details style="margin-top:8px"><summary style="cursor:pointer;font-size:12px;color:#7dd3fc">🔀 Git 狀況 / 各 AI 分支（快照 {ts}；重跑 render-html 刷新）</summary>'
        f'<div style="padding:6px 0;font-size:11px"><div>當前分支：<code>{html.escape(branch)}</code></div>'
        f'<div style="margin-top:5px;color:#94a3b8">最近 5 commit：</div>{loglines}'
        f'<div style="margin-top:5px;color:#94a3b8">worktrees（各 AI 分支）：</div>{wtlines}'
        f'<div style="margin-top:5px;color:#94a3b8">分支 ↔ 任務（碰撞=同分支多 in_progress）：</div>{btrows}</div></details>'
        '</div>')


def focus_console_html(data):
    """Two-column focus console: left sticky nav tree | right per-node detail panel."""
    nodes = data["nodes"]
    by, kids, _ = tree_maps(nodes)
    right = "".join(detail_panel_html(n["id"], data, by, kids) for n in nodes)
    # hide build_workstation's flat judgment grid (verdict engine kept; we drive it from the right panel)
    css = ("<style>"
           ".panel:has(#grid){display:none!important}h2:has(+ .panel #grid){display:none!important}"
           "details[data-tid]>summary.tgsum:hover{background:#16233a;border-radius:5px}"
           ".tgfocus{cursor:pointer;padding:1px 4px;border-radius:4px;transition:background .1s}.tgfocus:hover{background:#1d4ed8;color:#fff}"
           ".tg-h{font-size:11px;color:#7dd3fc;margin:13px 0 5px;font-weight:700;letter-spacing:.4px;border-left:3px solid #38bdf8;padding-left:7px}"
           ".tg-detail .kv{padding:3px 0;border-bottom:1px solid #182234;font-size:12px}.tg-detail .jb{flex:1;font-size:11px;padding:6px 2px}"
           ".tg-detail a{word-break:break-all}"
           ".tg-nb{display:flex;gap:10px;align-items:stretch;background:#0b1222;border:1px solid #334155;border-radius:9px;padding:11px;overflow-x:auto}"
           ".tg-nbcol{flex:1;min-width:118px;display:flex;flex-direction:column;gap:6px}"
           ".tg-nbh{font-size:10px;color:#94a3b8;font-weight:700;text-align:center;border-bottom:1px solid #334155;padding-bottom:4px;margin-bottom:2px}"
           ".tg-arrow{display:flex;align-items:center;color:#64748b;font-size:18px}"
           ".tg-pill{background:#1e293b;border:1px solid #334155;border-left-width:3px;border-radius:6px;padding:4px 8px;font-size:10.5px;color:#cbd5e1;line-height:1.5}"
           ".tg-pill-s{float:right;font-size:9px;font-weight:700}"
           ".tg-taskbox{background:#11202f;border:2px solid #334155;border-radius:8px;padding:7px 9px;text-align:center;line-height:1.45}"
           ".tg-focusbox{box-shadow:0 0 0 2px #f2c94c;background:#17243f}"
           ".tg-l0{font-size:13.5px;font-weight:700;color:#f2c94c;margin:2px 0 4px;line-height:1.45}"
           ".tg-l1{font-size:12px;color:#cbd5e1;margin-bottom:4px;line-height:1.55}"
           ".tg-status{background:#0b1222;border:1px solid #3b4a63;border-radius:10px;padding:12px 14px;margin-bottom:14px}"
           ".tg-chip{display:inline-block;background:#0b2417;border:1px solid #22c55e;border-radius:12px;padding:1px 8px;margin:2px 3px 2px 0;font-size:10.5px;color:#86efac;cursor:pointer}"
           ".tg-gitline{font-family:monospace;font-size:10.5px;color:#cbd5e1;padding:1px 0;word-break:break-all}"
           "</style>")
    return (css + status_header_html(data)
            + '<div style="display:flex;gap:18px;align-items:flex-start;flex-wrap:wrap">'
            f'<div style="flex:1 1 350px;min-width:300px;position:sticky;top:14px;max-height:calc(100vh - 36px);overflow:auto;background:#0b1222;border:1px solid #334155;border-radius:10px;padding:12px">{_tree_left(data)}</div>'
            f'<div style="flex:1 1 500px;min-width:340px;background:#1e293b;border:1px solid #3b4a63;border-radius:10px;padding:16px;box-shadow:0 2px 12px rgba(0,0,0,.3)" id="tg-focus-panel">{right}</div>'
            '</div>')


def focus_js(data):
    nodes = data["nodes"]
    focus = data.get("focus")
    tparent = {n["id"]: n.get("parent") for n in nodes}
    tgt = {n["id"]: n["title"] for n in nodes}
    init = f"focusFlow({json.dumps(focus)});" if focus else ""
    # NOTE: do NOT declare LSKEY here — build_workstation's main script already declares a
    # top-level `const LSKEY`; classic scripts share the global lexical scope, so re-declaring
    # throws "Identifier already declared" and silently kills this whole script. Reuse the global.
    return ("const TPARENT=" + json.dumps(tparent, ensure_ascii=False) + ",TGT=" + json.dumps(tgt, ensure_ascii=False) + ";"
            "function tgDesc(f){let s=new Set([f]),ch=true;while(ch){ch=false;for(const c in TPARENT){if(s.has(TPARENT[c])&&!s.has(c)){s.add(c);ch=true;}}}return s;}"
            "function tgPanelHi(id){var J={};try{J=JSON.parse(localStorage.getItem(LSKEY)||'{}')}catch(e){}var c=(J[id]||{}).c||'un';"
            "var box=document.querySelector('[data-vfor=\"'+id+'\"]');if(box)box.querySelectorAll('.jb').forEach(function(b){b.className='jb'+(b.getAttribute('data-v')===c?(' on-'+c):'');});}"
            "function tgVerdict(id,c){if(window.setJ)window.setJ(id,c);tgPanelHi(id);}"
            "function focusFlow(f){"
            "document.querySelectorAll('.tg-detail').forEach(e=>e.style.display='none');var dt=document.getElementById('detail-'+f);if(dt){dt.style.display='block';tgPanelHi(f);}"
            "var c=TPARENT[f];while(c){var ae=document.querySelector('details[data-tid=\"'+c+'\"]');if(ae)ae.open=true;c=TPARENT[c];}"
            "document.querySelectorAll('[data-tid]').forEach(e=>e.style.boxShadow='');var fe=document.querySelector('[data-tid=\"'+f+'\"]');if(fe)fe.style.boxShadow='0 0 0 2px #f2c94c';"
            "const sel=document.getElementById('tg-focus');if(sel)sel.value=f;"
            "const ch=[];let x=f;while(x){ch.unshift(x);x=TPARENT[x];}const bc=document.getElementById('tg-bc');"
            "if(bc)bc.innerHTML=ch.map((id,i)=>'<a style=\"cursor:pointer;color:'+(i===ch.length-1?'#f2c94c':'#60a5fa')+'\" onclick=\"focusFlow(\\''+id+'\\')\">'+(TGT[id]||id)+'</a>').join(' ＞ ');}"
            + init)


# ----------------------------------------------------------------- render TASKS md
def render_md(data):
    nodes = data["nodes"]
    by = {n["id"]: n for n in nodes}
    errs = validate(data)
    findings, stale = check(data)
    sha, branch = git_info()
    cnt = {s: sum(1 for n in nodes if n.get("status") == s) for s in VALID_STATUS}
    gaps = detail_gaps(nodes)
    v3 = is_v3(data)
    L = ["# 研究任務" + ("樹（WBS）" if v3 else "有向圖") + " — " + os.path.basename(TASKS_MD), "",
         "> **自動產生，請勿手改。** 唯一真值 = `" + os.path.relpath(GRAPH, ROOT) + "`；改後重跑 `python3 scripts/tasks/task_graph.py --graph " + os.path.relpath(GRAPH, ROOT) + " render`。",
         f"> 生成 {datetime.date.today().isoformat()} · build `{sha}` @ `{branch}` · 節點 {len(nodes)}（✅{cnt['done']} ◐{cnt['in_progress']} ⛔{cnt['blocked']} ☐{cnt['todo']} ✗{cnt['dropped']}）"
         + (f" · focus = `{data.get('focus')}`" if v3 and data.get("focus") else ""), "",
         "## ⓪ 驗證（validate + check）",
         f"- validate: {'✅ PASS' if not errs else '❌ FAIL'}"]
    L += [f"  - ❌ {e}" for e in errs]
    L.append(f"- check: {'✅ 0 findings' if not findings else str(len(findings))+' findings'}")
    L += [f"  - {k} `{nid}`: {m}" for k, nid, m in findings]
    L.append("")
    if v3:
        by_, kids, roots = tree_maps(nodes)
        focus = data.get("focus")
        if focus:
            L.append("## 聚焦路徑　" + " ＞ ".join(by_[i]["title"] for i in ancestors(focus, by_) + [focus]))
            L.append("")
        L.append("## 任務樹（含括 parent；縮排=層級）")

        def emit(nid, depth):
            n = by_[nid]
            io = n.get("io") or {}
            req = [i["name"] for i in io.get("inputs", []) if i.get("required")]
            ioabbr = ("　in:" + ",".join(("●" + trunc(i.get("name", ""), 12)) for i in io.get("inputs", [])[:2])) if io.get("inputs") else ""
            mk = "★" if nid == focus else STATUS_MARK.get(n.get("status"), "·")
            miss = ("　缺:" + "；".join(n["missing_info"])) if n.get("missing_info") else ""
            lk = n.get("links") or {}
            ref = lk.get("cycle_id") or (lk.get("memory") or [None])[0]
            L.append("  " * depth + f"- {mk} `{nid}` {n['title']} [{KIND_LABEL.get(n.get('kind'),'')}/{STATUS_LABEL.get(n.get('status'),'')}]" + ioabbr + miss + (f"　↳{ref}" if ref else ""))
            for c in kids.get(nid, []):
                emit(c, depth + 1)

        for r in roots:
            emit(r, 0)
        L.append("")
    else:
        L.append("## 各 component × 子任務")
        for c in data.get("components", []):
            members = sorted([n for n in nodes if n.get("component") == c["id"]], key=lambda x: x['id'])
            d = sum(1 for m in members if m.get("status") == "done")
            L.append(f"### {'★ ' if c.get('is_main') else ''}{c['title']} [{d}/{len(members)} done]")
            for m in members:
                L.append(f"- {STATUS_MARK.get(m.get('status'),'·')} `{m['id']}` {m['title']}")
            L.append("")
    rdy = ready(nodes)
    L.append(f"## Ready — 可立即開跑（{len(rdy)}）")
    L += [f"- `{n['id']}` {n['title']}　owner: {n.get('owner') or '未指派'}" for n in sorted(rdy, key=lambda x: x['id'])] or ["- （無）"]
    L.append("")
    L.append("## 細節不清楚 / 缺資訊")
    L += ([f"- `{nid}`: {why}" for nid, why in gaps] if gaps else ["- ✅ 無結構性細節缺口"])
    mi = [(n['id'], n['missing_info']) for n in nodes if n.get('missing_info')]
    if mi:
        L.append("**各任務 missing_info：**")
        L += [f"- `{nid}`: {'；'.join(info)}" for nid, info in mi]
    L.append("\n---\n單一所有權：本層只擁含括/依賴/分派/狀態/缺漏/io；verdict・證據・敘述歸 cycle/ledger/memory。")
    with open(TASKS_MD, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return errs, findings


# ----------------------------------------------------------------- render HTML
def load_build_ws():
    spec = importlib.util.spec_from_file_location("build_workstation", BUILD_WS)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _item_cards(nodes, indeg, outdeg):
    items = []
    by = {x["id"]: x for x in nodes}
    for n in nodes:
        nid, st = n["id"], n.get("status", "todo")
        lk, io = n.get("links") or {}, n.get("io") or {}
        badges = [{"label": STATUS_LABEL.get(st, st), "kind": STATUS_BADGE.get(st, "mut")},
                  {"label": KIND_LABEL.get(n.get("kind"), ""), "kind": "mut"}]
        if n.get("owner"):
            badges.append({"label": n["owner"], "kind": "info"})
        metrics = [{"k": "狀態", "v": STATUS_LABEL.get(st, st)}, {"k": "類型", "v": KIND_LABEL.get(n.get("kind"), "")}]
        if n.get("parent"):
            metrics.append({"k": "上層", "v": n["parent"], "src": n["parent"]})
        for inp in io.get("inputs", []):
            metrics.append({"k": ("輸入●必需" if inp.get("required") else "輸入○可選"),
                            "v": inp.get("name", "") + ("" if inp.get("ref") else "（外部）"), "src": inp.get("ref")})
        for o in io.get("outputs", []):
            metrics.append({"k": "輸出▶", "v": o})
        if n.get("depends_on"):
            metrics.append({"k": "依賴(先後)", "v": ", ".join(n["depends_on"])})
        metrics.append({"k": "被依賴", "v": str(outdeg.get(nid, 0))})
        metrics.append({"k": "owner", "v": n.get("owner") or "未指派"})
        if n.get("missing_info"):
            metrics.append({"k": "缺", "v": "；".join(n["missing_info"])})
        if lk.get("cycle_id"):
            metrics.append({"k": "cycle", "v": lk["cycle_id"], "src": lk["cycle_id"]})
        if lk.get("memory"):
            metrics.append({"k": "memory", "v": ", ".join(lk["memory"]), "src": "memory/" + lk["memory"][0]})
        if lk.get("reports"):
            metrics.append({"k": "report", "v": lk["reports"][0], "src": lk["reports"][0]})
        # --- B: per-task readiness / impact / observe / handoff ---
        rd = input_readiness(n, by)
        if rd:
            metrics.append({"k": "輸入就緒", "v": " ｜ ".join(READY_MARK[i["state"]] + " " + i["name"] for i in rd)})
        imp = impact(nid, nodes)
        if imp:
            metrics.append({"k": "完成後解鎖", "v": ", ".join(imp)})
        if n.get("observe"):
            metrics.append({"k": "觀察材料", "v": "；".join(n["observe"])})
        handoff = "；".join(io.get("outputs", [])) or "（見 cycle/memory/report links）"
        metrics.append({"k": "交接(其他AI讀)", "v": handoff})
        # --- C: verify / on_fail / problem ---
        if n.get("verify"):
            metrics.append({"k": "驗證方式", "v": n["verify"].get("how", ""), "src": n["verify"].get("ref")})
        if n.get("on_fail"):
            metrics.append({"k": "驗證失敗回饋", "v": n["on_fail"]})
        if n.get("problem"):
            metrics.append({"k": "⚠問題", "v": n["problem"].get("summary", ""), "src": n["problem"].get("ref")})
        guide = html.escape(n.get("notes", "") or "")
        p = n.get("problem")
        if p:
            guide += "<br><b style='color:#fca5a5'>⚠ 問題：</b>" + html.escape(p.get("summary", ""))
            if p.get("cause"):
                guide += "<br><b>根因：</b>" + html.escape(p["cause"])
            if p.get("options"):
                guide += "<br><b>選項：</b>" + html.escape("；".join(p["options"]))
            if p.get("ref"):
                guide += "<br><b>出處：</b><code>" + html.escape(p["ref"]) + "</code>"
            if p.get("fix_child"):
                guide += "<br><b>修任務：</b><code>" + html.escape(p["fix_child"]) + "</code>"
        items.append({"id": nid, "title": n["title"], "subtitle": f'{n.get("parent") or n.get("component") or ""} · {KIND_LABEL.get(n.get("kind"),"")} · {nid}',
                      "badges": badges, "metrics": metrics, "meta_metrics": ["狀態", "類型"],
                      "figures": {"mode": "svg", "svg": node_mini_svg(n, indeg.get(nid, 0), outdeg.get(nid, 0))},
                      "modal_stats": metrics, "reading_guide": guide})
    return items


def build_spec(data):
    nodes = data["nodes"]
    by = {n["id"]: n for n in nodes}
    errs = validate(data)
    findings, stale = check(data)
    indeg, outdeg = degrees(nodes)
    sha, branch = git_info()
    cnt = {s: sum(1 for n in nodes if n.get("status") == s) for s in VALID_STATUS}
    v3 = is_v3(data)

    fcol = {"WARN": "#f59e0b", "DRIFT": "#f59e0b", "DETAIL": "#a78bfa"}
    vrows = (f'<div class="kv"><span class="k">validate</span><span class="v" style="color:{"#22c55e" if not errs else "#ef4444"}">{"✅ PASS" if not errs else f"❌ FAIL ({len(errs)})"}</span></div>')
    for e in errs:
        vrows += f'<div class="kv"><span class="k">❌</span><span class="v">{html.escape(e)}</span></div>'
    vrows += f'<div class="kv"><span class="k">check</span><span class="v" style="color:{"#22c55e" if not findings else "#f59e0b"}">{"✅ 0" if not findings else str(len(findings))+" findings"}</span></div>'
    for k, nid, m in findings:
        vrows += f'<div class="kv"><span class="k" style="color:{fcol.get(k,"#94a3b8")}">{k} {html.escape(nid or "")}</span><span class="v">{html.escape(m)}</span></div>'

    gaps = detail_gaps(nodes)
    dg = "<b>結構性細節缺口：</b>" + (("<ul>" + "".join(f"<li><code>{html.escape(nid)}</code> {html.escape(why)}</li>" for nid, why in gaps) + "</ul>") if gaps else " <span style='color:#22c55e'>✅ 無</span><br>")
    mi = [(n['id'], n['missing_info']) for n in nodes if n.get('missing_info')]
    dg += "<b>各任務缺的資訊（missing_info）：</b><ul>" + "".join(f"<li><code>{html.escape(nid)}</code> <span style='color:#f59e0b'>{html.escape('；'.join(info))}</span></li>" for nid, info in mi) + "</ul>"

    rdy = ready(nodes)
    rb = "<b>Ready（可立即開跑）：</b><ul>" + ("".join(f"<li><code>{html.escape(n['id'])}</code> {html.escape(n['title'])} <span class='sub'>· {html.escape(n.get('owner') or '未指派')}</span></li>" for n in sorted(rdy, key=lambda x: x['id'])) or "<li class='sub'>（無）</li>") + "</ul>"

    if v3:
        probs = [n for n in nodes if n.get("problem")]
        ph = ("<div class='sub'>問題只連既有 postmortem/ledger/cycle（不重造驗證機制）；流程＝問題→建 child 修任務→verify→on_fail 回饋環。</div>")
        if probs:
            for n in probs:
                p = n["problem"]
                ph += (f'<div class="panel" style="border-color:#7f1d1d"><b style="color:#fca5a5">⚠ {html.escape(n["id"])} {html.escape(n["title"])}</b>'
                       f'<div class="kv"><span class="k">現象</span><span class="v">{html.escape(p.get("summary",""))}</span></div>')
                if p.get("cause"):
                    ph += f'<div class="kv"><span class="k">根因</span><span class="v">{html.escape(p["cause"])}</span></div>'
                if p.get("options"):
                    ph += f'<div class="kv"><span class="k">選項</span><span class="v">{html.escape("；".join(p["options"]))}</span></div>'
                if p.get("ref"):
                    ph += f'<div class="kv"><span class="k">出處(連既有)</span><span class="v"><code>{html.escape(p["ref"])}</code></span></div>'
                if p.get("fix_child"):
                    fc = p["fix_child"]
                    ph += f'<div class="kv"><span class="k">修任務(child)</span><span class="v"><code>{html.escape(fc)}</code> {html.escape(by.get(fc,{}).get("title",""))}</span></div>'
                if n.get("verify"):
                    ph += f'<div class="kv"><span class="k">驗證</span><span class="v">{html.escape(n["verify"].get("how",""))} ↳<code>{html.escape(n["verify"].get("ref") or "")}</code></span></div>'
                if n.get("on_fail"):
                    ph += f'<div class="kv"><span class="k">失敗回饋</span><span class="v">{html.escape(n["on_fail"])}</span></div>'
                ph += "</div>"
        else:
            ph += "<div class='sub'>✅ 無標記問題</div>"
        # verify-only nodes (no problem) — show their verification spec too
        vonly = [n for n in nodes if n.get("verify") and not n.get("problem")]
        if vonly:
            ph += "<b>驗證規格（無問題的任務）：</b><ul>" + "".join(
                f'<li><code>{html.escape(n["id"])}</code> {html.escape(n["verify"].get("how",""))} ↳<code>{html.escape(n["verify"].get("ref") or "")}</code>'
                + (f'　<span class="sub">失敗回饋：{html.escape(n["on_fail"])}</span>' if n.get("on_fail") else "") + "</li>" for n in vonly) + "</ul>"
        gate_nodes = [n for n in nodes if n.get("is_gate")]
        gh = '<div class="sub">OPEN 決策 gate 跨切彙整（掛在各領域下，但影響論文 Grade/tier/故事/精確度）：</div><table><tr><th>gate</th><th>狀態</th><th>領域</th><th>目標/完成判準</th></tr>'
        for gn in sorted(gate_nodes, key=lambda x: x["id"]):
            gh += (f'<tr><td><code style="cursor:pointer;color:#60a5fa" onclick="focusFlow(\'{gn["id"]}\')">{html.escape(gn["id"])}</code> {html.escape(gn["title"][:28])}</td>'
                   f'<td style="color:{STATUS_COLOR.get(gn.get("status"),"#94a3b8")}">{STATUS_MARK.get(gn.get("status"),"·")} {STATUS_LABEL.get(gn.get("status"),"")}</td>'
                   f'<td>{html.escape(gn.get("parent",""))}</td><td class="sub">{html.escape((gn.get("goal") or "")[:50])}</td></tr>')
        gh += "</table>"
        folded = (
            (f'<details open><summary style="cursor:pointer;color:#f2c94c">🚩 決策 gates 彙整（{len(gate_nodes)} 個 OPEN 決策，跨領域）</summary><div style="padding:8px 0">{gh}</div></details>' if gate_nodes else "")
            + f'<details><summary style="cursor:pointer">驗證結果（validate + check）</summary><div style="padding:8px 0">{vrows}</div></details>'
            f'<details><summary style="cursor:pointer">Ready — 可立即開跑</summary><div style="padding:8px 0">{rb}</div></details>'
            f'<details><summary style="cursor:pointer">細節不清楚 / 缺資訊</summary><div style="padding:8px 0">{dg}</div></details>'
            f'<details><summary style="cursor:pointer">問題與驗證回饋彙整（連既有 postmortem/run-evaluator/ledger）</summary><div style="padding:8px 0">{ph}</div></details>'
            f'<details><summary style="cursor:pointer">完整依賴 DAG（跨任務先後；實線必需／虛線可選）</summary><div style="padding:8px 0">{svg_dag(nodes)}</div></details>')
        sections = [
            {"title": "① 焦點控制台（左樹導航 ｜ 右焦點面板：流程圖 + I/O就緒 + 影響 + 觀察 + 問題/驗證 + 判讀）", "html": focus_console_html(data)},
            {"title": "② 其他檢視（點開）", "html": folded},
        ]
        title = "研究任務樹 — 聚焦控制台（ISM 小樹）"
        subtitle = "左樹導航｜右焦點面板（點任一節點即更新）· 互動流程圖 + I/O就緒/影響/觀察 + 問題/驗證回饋 · 非自動執行器 · focus=" + (data.get("focus") or "")
        lskey = "task_graph_v3_" + os.path.splitext(os.path.basename(GRAPH))[0]
    else:
        sections = [
            {"title": "① 程式流程圖（●必需 ○可選輸入）", "html": svg_pipeline(nodes)},
            {"title": "② 完整依賴 DAG（實線必需／虛線可選）", "html": svg_dag(nodes)},
            {"title": "③ 驗證結果", "html": f"<div>{vrows}</div>"},
            {"title": "④ 細節不清楚 / 缺資訊", "html": dg},
            {"title": "⑤ Ready 可開跑", "html": rb},
        ]
        title, subtitle, lskey = "研究任務有向圖 — 確認與驗證看板", "Subclonal reconstruction · 任務 DAG", "task_graph_v2"

    changelog = {"data_status": {"summary": "節點 status 由 memory 種子化(last-known)；請逐項確認/校正"},
                 "phase": {"summary": f"{len(nodes)} 節點 · ✅{cnt['done']} ◐{cnt['in_progress']} ⛔{cnt['blocked']} ☐{cnt['todo']} ✗{cnt['dropped']}"},
                 "audit_conclusion": {"summary": "人工分派地圖 + 驗證顯示，非自動執行器；任何開跑由使用者決定。"}}
    if stale:
        changelog["corrections"] = [{"id": cid, "what": f"active.json 此 cycle stale {age} 天", "status": "superseded", "effect": "需 /conclude-research 退役", "src": "state/active.json"} for cid, age in stale.items()]

    return {
        "meta": {"title": title, "subtitle": subtitle, "build_commit": sha, "build_branch": branch, "tier": "harness", "scope": "task-graph-layer", "lskey": lskey},
        "banner": {"text": "節點 status 為 memory 種子化 last-known；請逐項『確認無誤/待查/有問題』校正；判讀存 localStorage 可匯出。"},
        "item_config": {"verdict_states": [{"key": "in", "label": "確認無誤"}, {"key": "doubt", "label": "待查"}, {"key": "ex", "label": "有問題"}],
                        "reason_options": [{"key": "", "label": "— 原因（可選）—"}, {"key": "status", "label": "狀態過時"}, {"key": "deps", "label": "依賴不對"}, {"key": "io", "label": "I/O 不對"}, {"key": "missing", "label": "缺漏未列"}, {"key": "other", "label": "其他"}],
                        "per_page": 40, "required_metrics": ["狀態"]},
        "sections": sections, "changelog": changelog, "items": _item_cards(nodes, indeg, outdeg),
        "extra_js": (focus_js(data) if v3 else None),
        "provenance": {"sources": [os.path.relpath(GRAPH, ROOT), "state/schemas/task.schema.json", "scripts/tasks/task_graph.py", "research/autoresearch/evidence_ledger.jsonl", "state/active.json"],
                       "note": "樹/流程圖/I/O/狀態由 graph JSON 注入（§13-A 缺值即 refuse）；純手刻 SVG + 原生 <details> + focus-switch JS（build_workstation extra_js 一次注入，不進序列化），零外部依賴。"},
    }


def render_html(data):
    m = load_build_ws()
    spec = build_spec(data)
    try:
        doc = m.build(spec)
    except m.Refuse as e:
        sys.stderr.write(f"[REFUSE §13-A] {e}\n")
        sys.exit(3)
    with open(BOARD_HTML, "w", encoding="utf-8") as fh:
        fh.write(doc)
    print(f"[ok] wrote {BOARD_HTML} ({os.path.getsize(BOARD_HTML)//1024} KB, {len(spec['items'])} nodes)")


# ----------------------------------------------------------------- v5 lifecycle helpers
def dump_graph(data, path):
    """Write graph JSON one-node-per-line (clean git diffs)."""
    parts = []
    for k in ("schema_version", "updated_at", "focus", "note"):
        if k in data:
            parts.append(f'  {json.dumps(k, ensure_ascii=False)}: {json.dumps(data[k], ensure_ascii=False)}')
    if "components" in data:
        parts.append('  "components": [\n' + ",\n".join("    " + json.dumps(c, ensure_ascii=False) for c in data["components"]) + '\n  ]')
    parts.append('  "nodes": [\n' + ",\n".join("    " + json.dumps(n, ensure_ascii=False) for n in data["nodes"]) + '\n  ]')
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("{\n" + ",\n".join(parts) + "\n}\n")


def _find(data, tid):
    for n in data["nodes"]:
        if n["id"] == tid:
            return n
    sys.stderr.write(f"[error] task {tid!r} not found in {os.path.relpath(GRAPH, ROOT)}\n")
    sys.exit(1)


def _status_to_contract(st):
    return {"done": "DONE", "blocked": "BLOCKED", "dropped": "BLOCKED"}.get(st, "NEEDS_CONTEXT")


def worktree_branches():
    out = {}
    try:
        txt = subprocess.check_output(["git", "worktree", "list"], cwd=ROOT, stderr=subprocess.DEVNULL).decode()
        for ln in txt.splitlines():
            m = re.search(r"\[([^\]]+)\]", ln)
            if m:
                out[m.group(1)] = ln.split()[0]
    except Exception:
        pass
    return out


def branch_task_map(nodes):
    bm = {}
    for n in nodes:
        b = (n.get("links") or {}).get("branch")
        if b:
            bm.setdefault(b, []).append(n)
    return bm


def cmd_claim(data, a):
    n = _find(data, a.task)
    n.setdefault("links", {})["branch"] = a.branch
    if a.worktree:
        n["links"]["worktree"] = a.worktree
    dump_graph(data, GRAPH)
    print(f"[ok] {a.task} → branch '{a.branch}'" + (f" / worktree '{a.worktree}'" if a.worktree else "")
          + "（在該 worktree 起 session；勿與他人同分支同時 in_progress；commit 前 git_branch_commit_guard 仍把關）")


def cmd_handoff(data, a):
    n = _find(data, a.task)
    by = {x["id"]: x for x in data["nodes"]}
    io, lk = n.get("io") or {}, n.get("links") or {}
    contract = n.get("handoff") or {"status": _status_to_contract(n.get("status")), "key_metrics": {}, "anomalies": [], "path_to_landed_doc": (lk.get("reports") or [None])[0]}
    pkg = {"task_id": n["id"], "title": n["title"], "headline": n.get("headline") or n["title"],
           "goal": n.get("goal"),
           "status": n.get("status"), "kind": n.get("kind"), "owner": n.get("owner"), "contract": contract,
           "inputs": [{"name": i.get("name"), "required": i.get("required"), "ref": i.get("ref"),
                       "ready": (by.get(i["ref"], {}).get("status") == "done") if i.get("ref") else None} for i in io.get("inputs", [])],
           "outputs": io.get("outputs", []), "verify": n.get("verify"), "links": lk,
           "branch": lk.get("branch"), "worktree": lk.get("worktree"),
           "next_unlocked": impact(n["id"], data["nodes"]), "missing_info": n.get("missing_info", []),
           "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
           "note": "AI handoff (CLAUDE.md §8 return-contract)；唯一真值在 graph，本檔為交接快照。"}
    hdir = os.path.join(ROOT, "state/tasks/handoff")
    os.makedirs(hdir, exist_ok=True)
    jp, mp = os.path.join(hdir, n["id"] + ".json"), os.path.join(hdir, n["id"] + ".md")
    with open(jp, "w", encoding="utf-8") as fh:
        json.dump(pkg, fh, ensure_ascii=False, indent=2)
    v = n.get("verify") or {}
    md = [f"# 交接包 — {n['id']} {n['title']}", "",
          f"- **狀態**: {n.get('status')} / contract `{contract.get('status')}`",
          f"- **重點(headline)**: {n.get('headline') or n['title']}",
          f"- **目標(goal/完成判準)**: {n.get('goal') or '—'}",
          f"- **分支/worktree**: {lk.get('branch') or '—'} / {lk.get('worktree') or '—'}",
          f"- **已驗證輸出(交接)**: {'；'.join(io.get('outputs', [])) or '—'}",
          f"- **驗證**: {v.get('how','—')} → result `{v.get('result','—')}` ↳ {v.get('evidence','—')}",
          f"- **landed doc**: {contract.get('path_to_landed_doc') or '—'}",
          f"- **後續解鎖**: {', '.join(pkg['next_unlocked']) or '—'}",
          f"- **缺資訊**: {'；'.join(n.get('missing_info', [])) or '—'}", "", "## 輸入就緒"]
    for i in pkg["inputs"]:
        rdy = "就緒" if i["ready"] else ("外部需備" if i["ready"] is None else "待產出")
        md.append(f"- {'●必需' if i['required'] else '○可選'} {i['name']} — {rdy}" + (f" ↳{i['ref']}" if i["ref"] else ""))
    md.append("\n## 關聯（連既有 SoT，可 grep / 開檔）")
    for kk in ("cycle_id", "memory", "reports", "ledger_entries"):
        if lk.get(kk):
            md.append(f"- {kk}: {lk[kk]}")
    md.append(f"\n_generated {pkg['generated_at']} · 重跑 handoff 更新_")
    with open(mp, "w", encoding="utf-8") as fh:
        fh.write("\n".join(md) + "\n")
    print(f"[ok] handoff → {os.path.relpath(jp, ROOT)}\n            {os.path.relpath(mp, ROOT)}")


def _parse_io_in(s):
    out = []
    for item in (s or "").split(","):
        item = item.strip()
        if not item:
            continue
        bits = item.split(":")
        ref = bits[2] if len(bits) > 2 and bits[2] else None
        out.append({"name": bits[0], "required": (len(bits) < 2 or bits[1].lower() in ("req", "required", "true", "1", "")), "ref": ref})
    return out


def cmd_add_task(data, a):
    if not a.id or not ID_RE.match(a.id):
        sys.stderr.write(f"[error] need --id matching ^T-[A-Z0-9-]+$ (got {a.id!r})\n"); sys.exit(1)
    if any(n["id"] == a.id for n in data["nodes"]):
        sys.stderr.write(f"[error] id {a.id} already exists\n"); sys.exit(1)
    if not a.title or not a.kind or not a.owner:
        sys.stderr.write("[error] need --title --kind --owner\n"); sys.exit(1)
    node = {"id": a.id, "title": a.title, "parent": a.parent, "kind": a.kind,
            "status": a.status or "todo", "depends_on": [d for d in (a.depends_on or "").split(",") if d], "owner": a.owner}
    if a.io_in or a.io_out:
        node["io"] = {"inputs": _parse_io_in(a.io_in), "outputs": [o for o in (a.io_out or "").split(";") if o]}
    if a.headline:
        node["headline"] = a.headline
    if a.summary:
        node["summary"] = a.summary
    if a.observe:
        node["observe"] = [o for o in a.observe.split(";") if o]
    node["added_at"] = datetime.date.today().isoformat()
    errs = validate({**data, "nodes": data["nodes"] + [node]})
    print("提議新節點：\n" + json.dumps(node, ensure_ascii=False, indent=2))
    if errs:
        print("\n❌ validate 失敗（未寫入）：")
        for e in errs:
            print("  -", e)
        sys.exit(1)
    if a.dry_run:
        print("\n（--dry-run：validate ✅ 但未寫入；去掉 --dry-run 才寫入。）")
        return
    data["nodes"].append(node)
    dump_graph(data, GRAPH)
    render_md(data)
    print(f"\n[ok] 已加入 {a.id}；跑 render-html 更新看板。")


def cmd_verify_result(data, a):
    n = _find(data, a.task)
    if not (a.passed or a.fail):
        sys.stderr.write("[error] need --pass or --fail\n"); sys.exit(1)
    n.setdefault("verify", {})
    n["verify"]["result"] = "pass" if a.passed else "fail"
    n["verify"]["verified_at"] = datetime.datetime.now().isoformat(timespec="seconds")
    if a.ref:
        n["verify"]["evidence"] = a.ref
    if a.passed:
        n["status"] = "done"
        dump_graph(data, GRAPH); render_md(data)
        print(f"[ok] {a.task} 驗證 PASS → status=done, evidence={a.ref or '—'}。跑 render-html 更新。")
    else:
        dump_graph(data, GRAPH); render_md(data)
        print(f"[記錄] {a.task} 驗證 FAIL（verify.result=fail）。")
        print("→ 半自動下一步（建子修任務並標 parent）：")
        print(f"   python3 scripts/tasks/task_graph.py --graph {os.path.relpath(GRAPH, ROOT)} new-problem --task {a.task} --summary '<失敗現象>' [--cause ... --options 'a;b' --postmortem InterSubMod/docs/... --reopen '<可重試條件>'] --dry-run")
        print("   實際驗證 RUN 走 /run-evaluator 或 /verification-loop；postmortem 走 /structured-tech-report。")


def cmd_new_problem(data, a):
    if not a.summary:
        sys.stderr.write("[error] need --summary\n"); sys.exit(1)
    parent = _find(data, a.task)
    fixid = parent["id"] + "-FIX"
    if any(n["id"] == fixid for n in data["nodes"]):
        fixid += "2"
    child = {"id": fixid, "title": "修：" + a.summary[:28], "parent": parent["id"], "kind": "compute",
             "status": "todo", "depends_on": [], "owner": parent.get("owner"),
             "io": {"inputs": [{"name": "失敗證據/觀察", "required": True, "ref": None}], "outputs": ["修正後重驗結果"]},
             "added_at": datetime.date.today().isoformat()}
    prob = {"summary": a.summary}
    if a.cause:
        prob["cause"] = a.cause
    if a.options:
        prob["options"] = [o for o in a.options.split(";") if o]
    if a.postmortem:
        prob["ref"] = a.postmortem
    if a.reopen:
        prob["reopen_threshold"] = a.reopen
    prob["fix_child"] = fixid
    print("提議：")
    print(f"  parent {parent['id']} → status=blocked, problem=" + json.dumps(prob, ensure_ascii=False))
    print("  新增子修任務：" + json.dumps(child, ensure_ascii=False))
    if a.dry_run:
        print("\n（--dry-run：未寫入；去掉 --dry-run 才寫入。postmortem 內容走 /structured-tech-report。）")
        return
    parent["problem"] = prob
    parent["status"] = "blocked"
    if not parent.get("on_fail"):
        parent["on_fail"] = f"見子修任務 {fixid}；postmortem 走 /structured-tech-report（13 段 SRE）"
    data["nodes"].append(child)
    dump_graph(data, GRAPH); render_md(data)
    print(f"\n[ok] parent 標 blocked + 子修任務 {fixid} 已建；跑 render-html 更新看板。")


# ----------------------------------------------------------------- main
def main():
    global GRAPH, TASKS_MD, BOARD_HTML
    ap = argparse.ArgumentParser(description="research task DAG/WBS layer driver")
    ap.add_argument("--graph", default=GRAPH, help="path to graph JSON (default state/tasks/graph.json)")
    ap.add_argument("cmd", choices=["validate", "check", "ready", "render", "render-html", "all",
                                    "claim", "handoff", "add-task", "verify-result", "new-problem"])
    ap.add_argument("--task")
    ap.add_argument("--branch")
    ap.add_argument("--worktree")
    ap.add_argument("--pass", dest="passed", action="store_true")
    ap.add_argument("--fail", action="store_true")
    ap.add_argument("--ref")
    ap.add_argument("--note")
    ap.add_argument("--dry-run", dest="dry_run", action="store_true")
    ap.add_argument("--id")
    ap.add_argument("--title")
    ap.add_argument("--parent")
    ap.add_argument("--kind", choices=VALID_KIND)
    ap.add_argument("--owner", choices=["claude", "user", "claude+user"])
    ap.add_argument("--status", choices=VALID_STATUS)
    ap.add_argument("--depends-on", dest="depends_on")
    ap.add_argument("--headline")
    ap.add_argument("--summary")
    ap.add_argument("--observe")
    ap.add_argument("--io-in", dest="io_in")
    ap.add_argument("--io-out", dest="io_out")
    ap.add_argument("--cause")
    ap.add_argument("--options")
    ap.add_argument("--postmortem")
    ap.add_argument("--reopen")
    a = ap.parse_args()
    GRAPH = a.graph if os.path.isabs(a.graph) else os.path.join(ROOT, a.graph)
    TASKS_MD, BOARD_HTML = out_paths(GRAPH)
    data = load_graph()
    nodes = data["nodes"]
    lifecycle = {"claim": cmd_claim, "handoff": cmd_handoff, "add-task": cmd_add_task,
                 "verify-result": cmd_verify_result, "new-problem": cmd_new_problem}
    if a.cmd in lifecycle:
        lifecycle[a.cmd](data, a)
        return
    if a.cmd in ("validate", "all"):
        errs = validate(data)
        if errs:
            print("validate: ❌ FAIL")
            for e in errs:
                print("  -", e)
            if a.cmd == "validate":
                sys.exit(1)
        else:
            print(f"validate: ✅ PASS ({len(nodes)} nodes{' / v3 tree' if is_v3(data) else ' / v2'}, no error)")
    if a.cmd in ("check", "all"):
        findings, _ = check(data)
        print("check: ✅ 0 findings" if not findings else f"check: {len(findings)} findings")
        for k, nid, m in findings:
            print(f"  {k} {nid}: {m}")
    if a.cmd in ("ready", "all"):
        rdy = ready(nodes)
        print(f"ready: {len(rdy)} dispatchable")
        for n in sorted(rdy, key=lambda x: x['id']):
            print(f"  {n['id']}  {n['title']}")
    if a.cmd in ("render", "all"):
        errs, findings = render_md(data)
        print(f"[ok] wrote {TASKS_MD} (validate {'PASS' if not errs else 'FAIL'}, {len(findings)} findings)")
    if a.cmd in ("render-html", "all"):
        render_html(data)


if __name__ == "__main__":
    main()
