#!/usr/bin/env python3
"""harness_health.py — InterSubMod harness self-audit health check (read-only).

Computes 6 L0 health lights by reading disk ground-truth and comparing against
doc claims. Emits a snapshot JSON to state/health_snapshots/ and regenerates a
standalone HTML dashboard (data inlined, localStorage diff vs previous snapshot).

NO agents, NO web, NO writes outside state/health_snapshots/. Pure grep/parse.
Run on demand or monthly: `python3 scripts/harness_health.py`.

Design: docs/references/migration/20260531_harness_audit_dashboard_design_01.md
Lights deliberately surface known drift so first run proves the dashboard works.
"""
import json
import os
import re
import glob
import datetime
import html as _html

REPO = "/big7_disk/liaoyoyo2001/InterSubMod"
SNAP_DIR = os.path.join(REPO, "state", "health_snapshots")
CLAUDE_MD = os.path.join(REPO, ".claude", "CLAUDE.md")
SETTINGS = os.path.join(REPO, ".claude", "settings.local.json")
CURRENT_FOCUS = os.path.join(REPO, "docs", "CURRENT_FOCUS.md")
LEDGER = os.path.join(REPO, "research", "autoresearch", "evidence_ledger.jsonl")
ACTIVE = os.path.join(REPO, "state", "active.json")
QUEUE = os.path.join(REPO, "research", "autoresearch", "hypothesis_queue.json")
COMPILE_MARKER = "/tmp/ism_cpp_pending_compile.txt"


def _read(path):
    try:
        with open(path, encoding="utf-8") as f:
            return f.read()
    except Exception:
        return ""


def _claim(pattern, text, default=None):
    """Largest integer matching a pattern group in text (handles drift-prone counts)."""
    nums = [int(m) for m in re.findall(pattern, text)]
    return max(nums) if nums else default


# ── disk ground truth ──────────────────────────────────────────────
def disk_counts():
    skills = len(glob.glob(os.path.join(REPO, ".claude/skills/*/SKILL.md")))
    agents = len(glob.glob(os.path.join(REPO, ".claude/agents/*.md")))
    hooks = len(glob.glob(os.path.join(REPO, "scripts/hooks/*.sh")))
    workflows = len(glob.glob(os.path.join(REPO, ".claude/workflows/*.js")))
    return {"skills": skills, "agents": agents, "hooks": hooks, "workflows": workflows}


# ── Light 1: COUNT DRIFT ───────────────────────────────────────────
def light_count(disk, cmd):
    sk = _claim(r"(\d+)\s*個?\s*SKILL\.md", cmd) or _claim(r"SKILL\.md.{0,12}?(\d{2})", cmd)
    ag = _claim(r"(\d+)\s*個\s*project agent", cmd) or _claim(r"18", cmd)
    hk = _claim(r"\*\*(\d+)\s*hook scripts", cmd)
    rows, worst = [], "GREEN"
    for name, dval, cval in [("skills", disk["skills"], sk), ("agents", disk["agents"], ag), ("hooks", disk["hooks"], hk)]:
        if cval is None:
            rows.append(f"{name}: disk={dval} doc=?(unparsed)")
            continue
        diff = abs(dval - cval)
        st = "GREEN" if diff == 0 else ("YELLOW" if diff == 1 else "RED")
        if st == "RED":
            worst = "RED"
        elif st == "YELLOW" and worst == "GREEN":
            worst = "YELLOW"
        rows.append(f"{name}: disk={dval} doc={cval} {'OK' if diff==0 else f'DRIFT {diff}'}")
    return worst, rows


# ── Light 2: HARD-GATE TRUTH (detects exit-2 hooks neutered by || exit 0) ──
def light_hardgate(settings, cmd):
    """Scan every wired hook command: a script that contains `exit 2` but whose
    wiring appends `|| exit 0` / `|| true` is NEUTERED (claims to block, doesn't)."""
    real, neutered, soft_intentional = [], [], []
    try:
        hooks = settings.get("hooks", {})
    except Exception:
        hooks = {}
    for event, blocks in hooks.items():
        for b in blocks:
            for h in b.get("hooks", []):
                command = h.get("command", "")
                m = re.search(r"(/[^\s]+\.sh)", command)
                if not m:
                    continue
                spath = m.group(1)
                body = _read(spath)
                has_exit2 = bool(re.search(r"^\s*exit 2\b", body, re.M)) or "exit 2" in body or "gate_exit" in body
                masked = bool(re.search(r"\|\|\s*exit 0", command) or re.search(r"\|\|\s*true", command))
                # intentionally-soft = script self-declares SOFT/advisory in its header
                declares_soft = bool(re.search(r"SOFT gate|advisory only|advisory\)", body))
                base = os.path.basename(spath)
                if has_exit2 and not masked:
                    real.append(base)
                elif has_exit2 and masked and declares_soft:
                    soft_intentional.append(base)
                elif has_exit2 and masked:
                    neutered.append((base, event))
    # claim parse scoped to the Hard Gate section's "校正 = **N 個**"
    claim = _claim(r"校正\s*=\s*\*\*(\d+)", cmd)
    rows = [f"real exit-2 (unmasked): {len(set(real))} → {sorted(set(real))}"]
    if neutered:
        rows.append("⚠ NEUTERED BUG (script says exit-2 but wiring '|| exit 0' eats it): "
                    + ", ".join(f"{n}({e})" for n, e in neutered))
    if soft_intentional:
        rows.append("intentionally-soft (self-declared advisory): " + ", ".join(sorted(set(soft_intentional))))
    if claim is not None:
        rows.append(f"doc claims {claim} core Hard Gates")
    worst = "RED" if neutered else "GREEN"
    return worst, rows, [n for n, _ in neutered]


# ── Light 3: TIER-GATE WIRED ───────────────────────────────────────
def light_tiergate(settings_raw):
    wired = "pre_tier_upgrade_check" in settings_raw
    marker_referenced = os.path.exists(os.path.join(REPO, "scripts/hooks/pre_tier_upgrade_check.sh"))
    if wired:
        return "GREEN", ["pre_tier_upgrade_check.sh WIRED (state.json exit-2, prose advisory)"]
    if marker_referenced:
        return "RED", ["pre_tier_upgrade_check.sh exists but NOT wired (4 SKILL.md assume it enforces ⭐4/5)"]
    return "YELLOW", ["pre_tier_upgrade_check.sh missing"]


# ── Light 4: STATE <-> DOC SYNC ────────────────────────────────────
def light_state(active, focus_text):
    rows, worst = [], "GREEN"
    concluded = []
    try:
        for c in active.get("recently_concluded", []):
            concluded.append(c.get("cycle_id", ""))
    except Exception:
        pass
    for cid in concluded:
        if not cid:
            continue
        # crude: cycle short id (e.g. cycle3) referenced near a 'todo/conclude/待' marker
        short = cid.split("-")[-1]
        if re.search(r"(待\s*conclude|需\s*/?conclude-research|P0[^\n]{0,40}cycle|尚未\s*conclude)", focus_text):
            if "cycle3" in focus_text or short in focus_text:
                rows.append(f"⚠ {cid} in active.recently_concluded but CURRENT_FOCUS still describes a pending conclude")
                worst = "RED"
    if not rows:
        rows.append(f"recently_concluded={concluded or 'none'}; no obvious state↔doc conflict")
    return worst, rows


# ── Light 5: LEDGER FRESHNESS ──────────────────────────────────────
def light_ledger(focus_text):
    last = ""
    try:
        with open(LEDGER, encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    last = line
    except Exception:
        pass
    led_date = None
    m = re.search(r'"?(?:date|timestamp|ts|concluded_at|ran_at)"?\s*[:=]\s*"?(\d{4}-\d{2}-\d{2})', last)
    if m:
        led_date = m.group(1)
    foc = re.search(r"##\s+(?:\S+\s+)?(\d{4}-\d{2}-\d{2})", focus_text)
    foc_date = foc.group(1) if foc else None
    rows = [f"ledger last entry date={led_date}; CURRENT_FOCUS header date={foc_date}"]
    worst = "GREEN"
    if led_date and foc_date:
        if led_date > foc_date:
            worst = "YELLOW"
            rows.append("⚠ ledger leads CURRENT_FOCUS → newest evidence not yet written back (tier over-claim risk)")
    return worst, rows


# ── Light 6: QUEUE HYGIENE ─────────────────────────────────────────
def light_queue(queue, active):
    queued = []
    try:
        items = queue if isinstance(queue, list) else queue.get("hypotheses", queue.get("queue", []))
        for h in items:
            if isinstance(h, dict) and str(h.get("status", "")).lower() == "queued":
                queued.append(h.get("id", h.get("hypothesis_id", "?")))
    except Exception:
        pass
    dead_terms = ["filter", "caller-f1", "caller_f1", "headroom", "methyl_filter"]
    rawq = json.dumps(queue, ensure_ascii=False).lower()
    suspect = [q for q in queued if any(t in rawq for t in dead_terms)]
    rows = [f"queued={len(queued)} {queued}"]
    worst = "GREEN"
    if len(queued) >= 3:
        worst = "RED" if suspect else "YELLOW"
        rows.append(f"⚠ {len(queued)} queued; some reference concluded-DEAD directions (filter/caller-f1) → run /pivot-direction")
    return worst, rows


# ── extra: stale compile marker (feeds Hard-Gate truth context) ────
def compile_marker_status():
    if os.path.exists(COMPILE_MARKER) and os.path.getsize(COMPILE_MARKER) > 0:
        try:
            mtime = datetime.datetime.fromtimestamp(os.path.getmtime(COMPILE_MARKER)).date().isoformat()
            n = sum(1 for _ in open(COMPILE_MARKER))
        except Exception:
            mtime, n = "?", "?"
        return {"present": True, "files": n, "mtime": mtime,
                "note": "stale marker → if compile gate were unmasked it would block ALL commits"}
    return {"present": False}


# ── HTML renderer (standalone, data inlined, localStorage diff) ────
def render_html(snap):
    data = json.dumps(snap, ensure_ascii=False)
    lights = snap["lights"]
    color = {"GREEN": "#1a9850", "YELLOW": "#e0a800", "RED": "#d73027"}
    cards = []
    for i, L in enumerate(lights, 1):
        c = color.get(L["status"], "#888")
        detail = "<br>".join(_html.escape(r) for r in L["rows"])
        cards.append(f"""
      <div class="light" style="border-left:6px solid {c}">
        <div class="lh"><span class="dot" style="background:{c}"></span>
          <b>{i} {_html.escape(L['name'])}</b><span class="st" style="color:{c}">{L['status']}</span></div>
        <div class="ld">{detail}</div>
      </div>""")
    overall = snap["overall"]
    oc = color.get("RED" if overall["red"] else ("YELLOW" if overall["yellow"] else "GREEN"), "#888")
    return f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>InterSubMod Harness Health</title>
<style>
 body{{font-family:-apple-system,"Segoe UI",Roboto,"Droid Sans Fallback",sans-serif;background:#f5f6f8;color:#222;margin:0;padding:24px;max-width:1000px;margin:0 auto}}
 h1{{font-size:20px;margin:0 0 2px}} .sub{{color:#666;font-size:13px;margin-bottom:16px}}
 .grid{{display:grid;grid-template-columns:1fr 1fr 1fr;gap:12px}}
 @media(max-width:760px){{.grid{{grid-template-columns:1fr}}}}
 .light{{background:#fff;border-radius:8px;padding:12px 14px;box-shadow:0 1px 3px rgba(0,0,0,.08)}}
 .lh{{display:flex;align-items:center;gap:8px;font-size:14px}} .lh .st{{margin-left:auto;font-weight:700;font-size:12px}}
 .dot{{width:12px;height:12px;border-radius:50%;display:inline-block}}
 .ld{{margin-top:8px;font-size:12px;color:#444;line-height:1.5;word-break:break-word}}
 .bar{{background:#fff;border-radius:8px;padding:12px 16px;margin:14px 0;box-shadow:0 1px 3px rgba(0,0,0,.08);border-left:6px solid {oc}}}
 .diff{{font-size:12px;color:#666;margin-top:6px}} code{{background:#eee;padding:1px 4px;border-radius:3px}}
</style></head><body>
<h1>InterSubMod HARNESS HEALTH</h1>
<div class="sub">snapshot {snap['ts']} · model {snap['settings'].get('model','?')} · effort {snap['settings'].get('effort','?')}
 · disk {snap['disk']['skills']} sk / {snap['disk']['agents']} ag / {snap['disk']['hooks']} hk / {snap['disk']['workflows']} wf</div>
<div class="bar"><b>OVERALL</b>: {overall['green']} GREEN · {overall['yellow']} YELLOW · {overall['red']} RED
 <div class="diff" id="diff">vs last snapshot: (computing…)</div></div>
<div class="grid">{''.join(cards)}</div>
<p class="sub">read-only · regenerate: <code>python3 scripts/harness_health.py</code> · design: docs/references/migration/20260531_harness_audit_dashboard_design_01.md</p>
<script>
 const SNAP = {data};
 try {{
   const prev = JSON.parse(localStorage.getItem('ism_harness_health') || 'null');
   const el = document.getElementById('diff');
   if (prev && prev.ts !== SNAP.ts) {{
     const pm = {{}}; (prev.lights||[]).forEach(l=>pm[l.name]=l.status);
     const ch = (SNAP.lights||[]).filter(l=>pm[l.name] && pm[l.name]!==l.status)
                  .map(l=>`${{l.name}}: ${{pm[l.name]}}→${{l.status}}`);
     el.textContent = 'vs last ('+prev.ts+'): ' + (ch.length? ch.join(' · ') : 'no light changed');
   }} else {{ el.textContent = 'vs last snapshot: (first run — baseline stored)'; }}
   localStorage.setItem('ism_harness_health', JSON.stringify(SNAP));
 }} catch(e) {{}}
</script>
</body></html>"""


def main():
    os.makedirs(SNAP_DIR, exist_ok=True)
    cmd = _read(CLAUDE_MD)
    settings_raw = _read(SETTINGS)
    try:
        settings = json.loads(settings_raw)
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
    l3 = light_tiergate(settings_raw)
    l4 = light_state(active, focus)
    l5 = light_ledger(focus)
    l6 = light_queue(queue, active)

    lights = [
        {"name": "COUNT DRIFT", "status": l1[0], "rows": l1[1]},
        {"name": "HARD-GATE TRUTH", "status": l2[0], "rows": l2[1]},
        {"name": "TIER-GATE WIRED", "status": l3[0], "rows": l3[1]},
        {"name": "STATE<->DOC SYNC", "status": l4[0], "rows": l4[1]},
        {"name": "LEDGER FRESH", "status": l5[0], "rows": l5[1]},
        {"name": "QUEUE HYGIENE", "status": l6[0], "rows": l6[1]},
    ]
    mk = compile_marker_status()
    if mk.get("present"):
        lights[1]["rows"].append(f"stale compile marker: {mk['files']} files @ {mk['mtime']} — {mk['note']}")

    sts = {st: sum(1 for L in lights if L["status"] == st) for st in ("GREEN", "YELLOW", "RED")}
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    snap = {
        "ts": now,
        "disk": disk,
        "settings": {"model": settings.get("model", "?"), "effort": settings.get("effortLevel", "?")},
        "lights": lights,
        "overall": {"green": sts["GREEN"], "yellow": sts["YELLOW"], "red": sts["RED"]},
        "neutered_hardgates": l2[2],
        "compile_marker": mk,
    }
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M")
    json_path = os.path.join(SNAP_DIR, f"{stamp}.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(snap, f, ensure_ascii=False, indent=1)
    html_path = os.path.join(SNAP_DIR, "dashboard.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(render_html(snap))

    # terminal L0 table
    sym = {"GREEN": "●", "YELLOW": "●", "RED": "●"}
    print(f"# InterSubMod Harness Health — {now}")
    print(f"  disk: {disk['skills']} skills / {disk['agents']} agents / {disk['hooks']} hooks / {disk['workflows']} workflows")
    print(f"  model={snap['settings']['model']} effort={snap['settings']['effort']}")
    print()
    for i, L in enumerate(lights, 1):
        print(f"  {sym[L['status']]} {i} {L['name']:<18} [{L['status']}]")
        for r in L["rows"]:
            print(f"        {r}")
    print()
    print(f"  OVERALL: {sts['GREEN']} GREEN · {sts['YELLOW']} YELLOW · {sts['RED']} RED")
    print(f"  snapshot → {json_path}")
    print(f"  dashboard → {html_path}")


if __name__ == "__main__":
    main()
