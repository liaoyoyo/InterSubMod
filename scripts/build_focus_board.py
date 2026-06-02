#!/usr/bin/env python3
"""build_focus_board.py — InterSubMod 焦點板自動 derive 渲染器.

§13 Layer A by-construction: 所有研究數字 (verdict / p-value / tier / counts) 由
data 來源注入 — 本腳本不 hardcode 任何研究 metric. 缺資料 → 顯示 '?' 不捏造.

Reads (機械事實 SoT, 不改寫):
  state/active.json                         active cycles 索引
  state/cycles/{id}/state.json              per-cycle phase/tier/gate/links
  research/autoresearch/evidence_ledger.jsonl  結果譜系 (verdict/decision/parent_cycles/corrects)
  research/autoresearch/hypothesis_queue.json  假說佇列計數
  state/health_snapshots/*.json (newest)    harness 健康 strip
Reads (人工策展, 唯一手維護層):
  state/focus.json                          now / pinned+next_step / background / dead
Writes (DERIVED, gitignored, 非 SoT):
  state/focus_board.html

用法:  python3 scripts/build_focus_board.py   (從 repo root 跑)
"""
import json
import glob
import html
import os
from datetime import datetime, timezone

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _p(*parts):
    return os.path.join(ROOT, *parts)


def load_json(path, default=None):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return default


def load_jsonl(path):
    rows = []
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rows.append(json.loads(line))
                except Exception:
                    continue
    except Exception:
        pass
    return rows


def esc(x):
    return html.escape(str(x if x is not None else ""))


PHASES = ["P0", "P1", "P2", "P3", "P4", "P5", "P6"]
PHASE_FULL = ["P0_REGISTER", "P1_PLAN", "P2_PRECHECK", "P3_PILOT",
              "P4_GENERALIZE", "P5_EVALUATE", "P6_COMMIT"]


def verdict_class(v):
    """Sentiment colour from verdict text (no hardcoded research value)."""
    s = (v or "").lower()
    if any(k in s for k in ("negative", "dead", "retracted", "fail")):
        return "v-neg"
    if any(k in s for k in ("provisional", "caveat", "modest", "conditional", "pending", "superseded")):
        return "v-warn"
    if any(k in s for k in ("positive", "validated", "robust", "confirm", "strengthened")):
        return "v-good"
    return "v-mut"


def stepper(phase):
    idx = PHASE_FULL.index(phase) if phase in PHASE_FULL else -1
    out = []
    for i, p in enumerate(PHASES):
        cls = "done" if i < idx else ("cur" if i == idx else "todo")
        out.append(f'<span class="step {cls}" title="{esc(PHASE_FULL[i])}">{p}</span>')
    return '<div class="stepper">' + "".join(out) + "</div>"


def tier_label(tier):
    if tier in (None, "", "pending"):
        return "tier pending"
    if tier in ("NEGATIVE", "RETRACTED"):
        return tier
    return "⭐" * int(tier) if str(tier).isdigit() else esc(tier)


def main():
    focus = load_json(_p("state", "focus.json"), {}) or {}
    active = load_json(_p("state", "active.json"), {}) or {}
    ledger = load_jsonl(_p("research", "autoresearch", "evidence_ledger.jsonl"))
    queue = load_json(_p("research", "autoresearch", "hypothesis_queue.json"), []) or []

    # ledger indexed by cycle_id
    led = {}
    for e in ledger:
        cid = e.get("cycle_id")
        if cid and cid not in led:
            led[cid] = e

    # cycle state.json per active cycle
    cyc_state = {}
    for c in active.get("cycles", []):
        cid = c.get("cycle_id")
        st = load_json(_p("state", "cycles", cid, "state.json"))
        if st:
            cyc_state[cid] = st

    # health strip (newest snapshot)
    snaps = sorted(glob.glob(_p("state", "health_snapshots", "*.json")))
    health = load_json(snaps[-1], {}) if snaps else {}
    ov = (health or {}).get("overall", {}) if isinstance(health, dict) else {}
    hg, hy, hr = ov.get("green", "?"), ov.get("yellow", "?"), ov.get("red", "?")

    # queue counts
    LIVE = {"pending", "adopted", "adopted_annotation", "testing", "queued", "active"}
    q_live = sum(1 for h in queue if h.get("status") in LIVE)
    q_total = len(queue)
    q_closed = q_total - q_live

    now_id = focus.get("now")
    pinned = focus.get("pinned", [])
    nextstep = {p["cycle"]: p.get("next_step", "") for p in pinned if isinstance(p, dict)}
    short = {p["cycle"]: p.get("short", "") for p in pinned if isinstance(p, dict)}

    gen = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # ---- L0 NOW hero ----
    now_st = cyc_state.get(now_id, {})
    now_full = now_st.get("title", now_id or "（focus.json 未設 now）")
    now_title = short.get(now_id) or now_full
    now_phase = now_st.get("phase", "")
    now_tier = tier_label(now_st.get("tier"))
    now_next = nextstep.get(now_id, "")
    l0 = f"""
  <section class="l0">
    <div class="l0-label">▶ 現在做什麼</div>
    <div class="l0-title" title="{esc(now_full)}">{esc(now_title)}</div>
    <div class="l0-next">{esc(now_next)}</div>
    <div class="l0-meta"><span class="badge b-tier">{now_tier}</span> {stepper(now_phase)}</div>
  </section>"""

    # ---- L1 pinned cards ----
    cards = []
    for p in pinned:
        cid = p["cycle"] if isinstance(p, dict) else p
        st = cyc_state.get(cid, {})
        gates = st.get("open_gates", []) or []
        gate_html = ""
        for g in gates:
            if g.get("status") == "open":
                gate_html += f'<span class="badge b-gate" title="{esc(g.get("reason",""))}">⚠ {esc(g.get("gate",""))}</span>'
        cards.append(f"""
      <div class="card">
        <h3 title="{esc(st.get('title', cid))}">{esc(short.get(cid) or st.get('title', cid))}</h3>
        <div class="row">
          <span class="badge b-tier">{tier_label(st.get('tier'))}</span>
          <span class="badge b-phase">{esc(st.get('phase',''))}</span>
          {gate_html}
        </div>
        {stepper(st.get('phase',''))}
        <div class="cid">{esc(cid)}</div>
        <div class="next"><b>下一步：</b>{esc(nextstep.get(cid,''))}</div>
      </div>""")
    l1 = f'<section><h2>★ 兩條主軸 <span class="sub">（pinned · 機械 cycle）</span></h2><div class="grid2">{"".join(cards)}</div></section>'

    # ---- 6 goals overview (見林) ----
    goals = focus.get("goals", [])
    GST = {"ACTIVE-PRIMARY": ("g-primary", "主軸"), "ACTIVE": ("g-active", "進行中"),
           "DEFERRED": ("g-defer", "暫緩"), "BLOCKED": ("g-block", "受阻"), "DEAD": ("g-dead", "關閉")}
    grows = []
    for g in goals:
        cls, lab = GST.get(g.get("status"), ("g-active", g.get("status", "")))
        cyc = g.get("cycle")
        cyclink = f'<span class="g-cyc">▶ {esc(str(cyc)[:32])}</span>' if cyc else ""
        grows.append(
            f'<div class="goal {cls}"><span class="g-id">{esc(g.get("id"))}</span>'
            f'<span class="g-st">{lab}</span><span class="g-short">{esc(g.get("short",""))}</span>{cyclink}</div>')
    goals_sec = ('<section><h2>🎯 6 目標全景 <span class="sub">（見林 · 哪些活 / 暫緩 / 受阻 / 關閉）</span></h2>'
                 f'<div class="goals">{"".join(grows)}</div></section>') if goals else ""

    # ---- blockers strip ----
    blockers = focus.get("blockers", [])
    blk_html = "".join(
        f'<span class="blk blk-{esc(b.get("sev","high"))}">🚧 <b>{esc(b.get("title",""))}</b>：{esc(b.get("detail",""))}</span>'
        for b in blockers)
    blk_sec = f'<div class="blockers">{blk_html}</div>' if blockers else ""

    # ---- L2 result lineage ----
    lin_blocks = []
    for p in pinned:
        cid = p["cycle"] if isinstance(p, dict) else p
        st = cyc_state.get(cid, {})
        refs = (st.get("links", {}) or {}).get("ledger_entries", []) or []
        entries = [led[r] for r in refs if r in led]
        entries.sort(key=lambda e: str(e.get("timestamp", "")))
        nodes = []
        for i, e in enumerate(entries):
            v = e.get("verdict", "?")
            dec = (e.get("decision") or "")[:90]
            corrects = e.get("corrects")
            tag = f'<span class="anno-corr">corrects {esc(str(corrects)[:18])}…</span>' if corrects else ""
            if i:
                nodes.append('<span class="arr">→</span>')
            nodes.append(
                f'<div class="node"><span class="nc">{esc(str(e.get("cycle_id",""))[:34])}</span>'
                f'<span class="nv {verdict_class(v)}">{esc(v[:46])}</span>'
                f'<span class="ndec">{esc(dec)}</span>{tag}</div>')
        rp = ""
        if entries:
            rps = [e.get("research_potential") for e in entries if e.get("research_potential")]
            if rps:
                rp = f'<div class="anno">→ <span class="spawn">research_potential</span>：{esc(str(rps[-1])[:160])}</div>'
        body = "".join(nodes) if nodes else '<span class="v-mut">（此 cycle state.json 未連 ledger_entries）</span>'
        lin_blocks.append(
            f'<details open><summary>{esc(st.get("title", cid))}</summary>'
            f'<div class="chain">{body}</div>{rp}</details>')
    l2 = '<section><h2>🧬 結果怎麼累積的 <span class="sub">（result lineage · 來自 evidence_ledger）</span></h2>' + "".join(lin_blocks) + "</section>"

    # ---- L3 detail row: background / queue / recent ----
    bg = focus.get("background", [])
    bg_html = "".join(
        f'<div class="mini"><span class="t">{esc(b.get("title",""))}</span><br>{esc(b.get("detail",""))}</div>'
        for b in bg) or '<div class="mini v-mut">（無）</div>'

    recent = sorted(ledger, key=lambda e: str(e.get("timestamp", "")), reverse=True)[:5]
    rec_html = "".join(
        f'<li><span class="nc">{esc(str(e.get("timestamp",""))[:10])}</span> '
        f'<span class="{verdict_class(e.get("verdict"))}">{esc(str(e.get("verdict",""))[:40])}</span> '
        f'<span class="v-mut">{esc(str(e.get("cycle_id",""))[:30])}</span></li>'
        for e in recent) or "<li class='v-mut'>（無）</li>"

    l3 = f"""
  <section><h2>▸ 背景 · 佇列 · 近期 <span class="sub">（L3 detail）</span></h2>
    {blk_sec}
    <div class="grid3">
      <div class="panel"><div class="ph">▸ 背景孵化（≤2）</div>{bg_html}</div>
      <div class="panel"><div class="ph">📋 假說佇列</div>
        <div class="big">{q_live}<span class="unit"> live</span></div>
        <div class="v-mut">{q_closed} closed / {q_total} total</div></div>
      <div class="panel"><div class="ph">🕐 近期結果（ledger tail）</div><ul class="rec">{rec_html}</ul></div>
    </div>
  </section>"""

    # ---- L3 dead ----
    dead = focus.get("dead", [])
    dead_html = "".join(
        f'<div><s>{esc(d.get("title",""))}</s> — {esc(d.get("why",""))}</div>' for d in dead)
    l3dead = f'<section><details><summary>❌ 已死路 / Concluded（透明保留 · 勿再開）</summary><div class="dead">{dead_html}</div></details></section>'

    # ---- assemble ----
    doc = f"""<!DOCTYPE html>
<!-- DERIVED by scripts/build_focus_board.py at {gen}. NON-SoT. 機械 SoT = active.json/cycles/evidence_ledger; 人工策展 = state/focus.json. 重生: python3 scripts/build_focus_board.py -->
<html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>InterSubMod 焦點板</title>
<style>
:root{{--bg:#f6f8fa;--card:#fff;--ink:#1f2328;--mut:#656d76;--line:#d0d7de;
--accent:#0969da;--accbg:#ddf4ff;--good:#1a7f37;--warn:#9a6700;--neg:#cf222e;--dead:#6e7781}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--ink);font-size:15px;line-height:1.5;padding-bottom:48px;
font-family:-apple-system,"Segoe UI","Droid Sans Fallback","Noto Sans CJK TC",sans-serif}}
.wrap{{max-width:1000px;margin:0 auto;padding:0 16px}}
.topbar{{background:var(--card);border-bottom:1px solid var(--line);padding:12px 16px}}
.topbar .wrap{{display:flex;flex-wrap:wrap;gap:8px 16px;align-items:baseline;justify-content:space-between}}
h1{{font-size:18px;margin:0}}
.gen{{color:var(--mut);font-size:12px}}
.hstrip{{font-size:12px;font-weight:600}}
.hg{{color:var(--good)}}.hy{{color:var(--warn)}}.hr{{color:var(--neg)}}
/* L0 hero */
.l0{{background:var(--accbg);border:1px solid #aac7ff;border-left:5px solid var(--accent);
border-radius:10px;padding:16px 18px;margin:20px 0}}
.l0-label{{font-size:12px;font-weight:700;color:var(--accent);letter-spacing:.05em}}
.l0-title{{font-size:21px;font-weight:700;margin:4px 0 2px}}
.l0-next{{font-size:15px;color:#24292f}}
.l0-meta{{display:flex;align-items:center;gap:10px;margin-top:10px;flex-wrap:wrap}}
h2{{font-size:15px;margin:24px 0 10px;padding-bottom:5px;border-bottom:2px solid var(--line)}}
h2 .sub{{font-weight:400;color:var(--mut);font-size:12.5px}}
.grid2{{display:grid;grid-template-columns:1fr 1fr;gap:14px}}
.grid3{{display:grid;grid-template-columns:1fr 1fr 1fr;gap:12px}}
@media(max-width:760px){{.grid2,.grid3{{grid-template-columns:1fr}}}}
.card{{background:var(--card);border:1px solid var(--line);border-left:4px solid var(--accent);border-radius:8px;padding:13px 15px}}
.card h3{{font-size:15px;margin:0 0 7px}}
.row{{display:flex;flex-wrap:wrap;gap:6px;align-items:center;margin-bottom:8px}}
.badge{{font-size:11.5px;font-weight:600;padding:1px 7px;border-radius:10px;white-space:nowrap}}
.b-tier{{background:#ddf4ff;color:#0a3069}}.b-phase{{background:#eaeafc;color:#3a3aa0}}.b-gate{{background:#fff8c5;color:var(--warn)}}
.stepper{{display:flex;gap:3px;margin:6px 0}}
.step{{flex:1;text-align:center;font-size:10px;font-weight:600;padding:3px 0;border-radius:3px;background:#eaecf0;color:#a0a6ad}}
.step.done{{background:#cae7d0;color:var(--good)}}.step.cur{{background:var(--accent);color:#fff}}
.cid{{font-family:ui-monospace,Menlo,monospace;font-size:10.5px;color:var(--mut);margin-top:6px}}
.next{{font-size:13.5px;margin-top:7px;padding-top:7px;border-top:1px dashed var(--line)}}.next b{{color:var(--accent)}}
details{{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:8px 14px;margin-top:8px}}
summary{{cursor:pointer;font-weight:600;font-size:13.5px}}
.chain{{display:flex;flex-wrap:wrap;align-items:stretch;margin:10px 0 2px}}
.node{{flex:0 1 auto;min-width:155px;max-width:250px;border:1px solid var(--line);border-radius:6px;padding:7px 9px;font-size:12px;background:#fbfcfd}}
.node .nc{{font-family:ui-monospace,monospace;font-size:10px;color:var(--mut);display:block}}
.node .nv{{font-weight:600;display:block;margin:2px 0}}
.node .ndec{{font-size:11px;color:var(--mut);display:block}}
.anno-corr{{font-size:10.5px;color:var(--warn)}}
.arr{{align-self:center;color:var(--mut);padding:0 5px}}
.anno{{font-size:11.5px;color:var(--mut);margin-top:4px}}.spawn{{color:var(--accent);font-weight:600}}
.panel{{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:11px 13px;font-size:13px}}
.panel .ph{{font-weight:600;font-size:12.5px;margin-bottom:6px;color:#24292f}}
.mini{{margin-bottom:7px}}.mini .t{{font-weight:600}}
.big{{font-size:26px;font-weight:700;color:var(--accent)}}.unit{{font-size:13px;color:var(--mut);font-weight:400}}
.rec{{margin:0;padding-left:16px}}.rec li{{margin:2px 0;font-size:11.5px}}.rec .nc{{font-family:ui-monospace,monospace;color:var(--mut)}}
.dead{{font-size:13px;color:var(--dead);padding:4px 0}}.dead s{{color:var(--neg)}}.dead div{{margin:4px 0}}
.v-good{{color:var(--good)}}.v-warn{{color:var(--warn)}}.v-neg{{color:var(--neg)}}.v-mut{{color:var(--mut)}}
.goals{{display:flex;flex-direction:column;gap:6px}}
.goal{{display:flex;align-items:center;gap:8px;background:var(--card);border:1px solid var(--line);border-left:4px solid var(--mut);border-radius:6px;padding:7px 11px;font-size:13px}}
.goal .g-id{{font-weight:700;font-family:ui-monospace,monospace;min-width:30px}}
.goal .g-st{{font-size:11px;font-weight:600;padding:1px 7px;border-radius:9px;white-space:nowrap}}
.goal .g-short{{flex:1;color:#24292f;min-width:180px}}
.goal .g-cyc{{font-family:ui-monospace,monospace;font-size:10.5px;color:var(--accent)}}
.g-primary{{border-left-color:var(--accent);background:#f0f7ff}}.g-primary .g-st{{background:var(--accent);color:#fff}}
.g-active{{border-left-color:var(--good)}}.g-active .g-st{{background:#cae7d0;color:var(--good)}}
.g-defer{{border-left-color:#8b949e}}.g-defer .g-st{{background:#eaecf0;color:var(--mut)}}
.g-block{{border-left-color:var(--warn)}}.g-block .g-st{{background:#fff8c5;color:var(--warn)}}
.g-dead{{border-left-color:var(--neg);opacity:.72}}.g-dead .g-st{{background:#ffd7d5;color:var(--neg)}}.g-dead .g-short{{text-decoration:line-through}}
.blockers{{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:10px}}
.blk{{font-size:12px;padding:4px 9px;border-radius:6px;border:1px solid #eedfa0;background:#fff8c5;color:var(--warn)}}
.blk-high{{background:#ffebe9;color:var(--neg);border-color:#ffc1bc}}
footer{{margin-top:28px;padding-top:12px;border-top:1px solid var(--line);font-size:12px;color:var(--mut)}}
code{{background:#eff1f3;padding:1px 5px;border-radius:4px;font-size:11.5px}}
</style></head><body>
<div class="topbar"><div class="wrap">
  <h1>🎯 InterSubMod Mission Control</h1>
  <span class="gen">DERIVED {esc(gen)} · 機械 SoT: active.json + /cycle-state</span>
  <span class="hstrip">harness <span class="hg">{hg}G</span>·<span class="hy">{hy}Y</span>·<span class="hr">{hr}R</span></span>
</div></div>
<div class="wrap">
{l0}
{l1}
{goals_sec}
{l2}
{l3}
{l3dead}
  <footer>
    <b>分層</b>：L0 現在做什麼 · L1 兩條主軸 · L2 結果譜系 · L3 背景/佇列/近期 + 已死路。<br>
    <b>SoT</b>：機械態 = <code>evidence_ledger.jsonl</code>(衝突時為準) + <code>state/cycles/*/state.json</code> + <code>hypothesis_queue.json</code>；
    人工策展 = <code>state/focus.json</code>(now/pinned/背景/dead，可手改可 commit)；
    本檔 = <b>derived</b>(gitignored，勿手改)。<br>
    <b>重生</b>：<code>python3 scripts/build_focus_board.py</code>。harness 系統健康獨立看 <code>/harness-health</code>；專案全局 <code>/research-dashboard</code>。
  </footer>
</div></body></html>"""

    out = _p("state", "focus_board.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(doc)
    print(f"[build_focus_board] wrote {out} ({len(doc)} bytes)")
    print(f"  now={now_id} | pinned={len(pinned)} | active_cycles={len(active.get('cycles',[]))}")
    print(f"  queue {q_live} live / {q_closed} closed / {q_total} total | ledger {len(ledger)} entries")
    print(f"  health {hg}G/{hy}Y/{hr}R (snapshot {os.path.basename(snaps[-1]) if snaps else 'none'})")


if __name__ == "__main__":
    main()
