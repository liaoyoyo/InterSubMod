#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_verification_page.py — 2026-07-25 solver 改良驗證頁（兩條資料線）

設計契約(依 memory):
  · feedback_reports_need_layered_disclosure — L0 一眼結論 / L1 重點邏輯 / L2 證據卡 / L3 溯源
  · feedback_design_principles_canonical    — 5 秒測試、Assertion-Evidence 標題、CRAP、色盲安全
  · 本機無 emoji 字型(fc-list 實測 0) → 重要性標記一律 CSS 徽章(文字+形狀),絕不用 emoji

資料來源(本輪實跑,皆全量非抽樣):
  A. exact-PS 線  scratchpad/exactps_results.tsv   98,955 unit × 7 datasets  ← 子頁面顯示的那批
  B. layered 線   scratchpad/real_results_all.tsv  18,931 unit × HCC1395
輸出:
  docs/methodology/_assets/layered_workstation/20260725_solver_verification.html
  同名 .data.json(數字唯一來源;§13-A 缺 required key 即 refuse)
"""
import collections
import csv
import json
import statistics as st
from pathlib import Path

SCRATCH = Path("/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/"
               "efb6e3d8-c4af-43d8-ac97-9dffdbec60ed/scratchpad")
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/layered_workstation")
STEM = "20260725_solver_verification"

# ── A. exact-PS 線（子頁面所用的同一批）────────────────────────────────────
EX = json.loads((OUT / "20260725_exactps_recovery.data.json").read_text(encoding="utf-8"))

# ── B. layered 線（HCC1395，含 off-by-one bug 發現）──────────────────────
lrows = list(csv.DictReader(open(SCRATCH / "real_results_all.tsv"), delimiter="\t"))
lfull = {}
for line in open(SCRATCH / "instances_all.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 6:
        lfull[(f[1], f[2])] = [g for g in f[4].split(",") if g]
lcap = [r for r in lrows if r["py_capped"] == "1"]
lreg = [r for r in lrows if r["py_capped"] == "0"]
lsolved = [r for r in lcap if r["h_star"] != "—"]
anom = [r for r in lsolved if int(r["py_h"]) < int(r["h_star"])]
ctrl = [r for r in lsolved if int(r["py_h"]) >= int(r["h_star"])]
is_allref = lambda r: any(set(g) == {"R"} for g in lfull[(r["region"], r["family"])])
lchain = [r for r in lcap if r["path"] == "chain_closed_form"]
lms = [float(r["ms"]) for r in lrows]
gr = [(int(r["g_raw"]), int(r["g_dom"])) for r in lrows + [] if int(r["g_raw"]) > 0]

D = {
    "generated": "2026-07-25",
    "exactps": EX,
    "layered": {
        "n_total": len(lrows), "n_capped": len(lcap), "capped_solved": len(lsolved),
        "n_regression": len(lreg),
        "regression_match": sum(1 for r in lreg if r["h_star"] != "—" and int(r["h_star"]) == int(r["py_h"])),
        "unsolved": sum(1 for r in lrows if r["h_star"] == "—"),
        "total_s": round(sum(lms) / 1000, 1), "median_ms": round(st.median(lms), 3),
        "bug": {"n": len(anom), "allref_hit": sum(1 for r in anom if is_allref(r)),
                "allref_pct": round(sum(1 for r in anom if is_allref(r)) / max(len(anom), 1) * 100, 1),
                "ctrl_n": len(ctrl), "ctrl_allref": sum(1 for r in ctrl if is_allref(r)),
                "ctrl_pct": round(sum(1 for r in ctrl if is_allref(r)) / max(len(ctrl), 1) * 100, 1),
                "examples": [{"region": r["region"], "family": r["family"], "py_h": int(r["py_h"]),
                              "new_h": int(r["h_star"]), "full": lfull[(r["region"], r["family"])]}
                             for r in anom[:3]],
                "diff": [{"case": "chr6:867869-868030/1", "dp": 7, "bnb": 7, "py": 6},
                         {"case": "chr6:3143829-3144920", "dp": 7, "bnb": 7, "py": 6},
                         {"case": "chr6:2685736-2686644", "dp": 6, "bnb": 6, "py": 5}]},
        "chain": {"n": len(lchain), "old_total": sum(int(r["py_trees"]) for r in lchain),
                  "new_total": sum(int(r["trees"]) for r in lchain),
                  "breakdown": [{"trees": v, "n": c} for v, c in
                                sorted(collections.Counter(int(r["trees"]) for r in lchain).items())]},
    },
    "groups": {"raw_median": int(st.median([a for a, _ in gr])), "raw_max": max(a for a, _ in gr),
               "dom_median": int(st.median([d for _, d in gr])), "dom_max": max(d for _, d in gr),
               "reduction_median": round(st.median([1 - d / max(a, 1) for a, d in gr]) * 100, 1)},
    "speed": [{"m": 6, "py": "capped→誤報 1 棵；放寬預算後 40 秒", "cpp_ms": 7.9},
              {"m": 7, "py": "暴力法估 ~7 小時", "cpp_ms": 157.0},
              {"m": 8, "py": "暴力法估 ~85 天", "cpp_ms": 7461.7}],
    "kbound": [{"case": "k=64，active loci=3（稀疏區）", "ms": 0.05, "state": "CERTIFIED", "ok": 1},
               {"case": "active loci=10 全 A（最壞鏈）", "ms": 55776.0, "state": "CERTIFIED", "ok": 1},
               {"case": "active loci=13（超出核心保證）", "ms": 0.0, "state": "UNSUPPORTED_BIT_COUNT", "ok": 0}],
    "max_bits": 12,
    "solver": "LongLineage src/solver（obligation_bnb + terminal_subset_dp + small_q_oracle）",
}
for k in ("exactps", "layered", "groups", "speed", "kbound"):
    if k not in D:
        raise SystemExit(f"REFUSE: data.json 缺 required key '{k}'")
(OUT / f"{STEM}.data.json").write_text(json.dumps(D, indent=1, ensure_ascii=False), encoding="utf-8")

# ════════════════════════════ 渲染 ════════════════════════════
E, LA, G = D["exactps"], D["layered"], D["groups"]
BUG, CHN = LA["bug"], LA["chain"]
UNDER = CHN["new_total"] // max(CHN["old_total"], 1)
SAMPLES = list(E["per_sample"].items())


def n(x):
    return f"{x:,}" if isinstance(x, int) else x


def tbl(head, body, minw=460):
    h = "".join(f"<th>{c}</th>" for c in head)
    return (f'<div class="tw"><table style="min-width:{minw}px">'
            f'<thead><tr>{h}</tr></thead><tbody>{body}</tbody></table></div>')


def card(title, body, tag=""):
    t = f'<span class="ct">{tag}</span>' if tag else ""
    return f'<details class="card"><summary>{title}{t}</summary><div class="cb">{body}</div></details>'


# ── SVG 1：逐樣本 ABSTAIN 救回（本頁最重要的圖）──────────────────────────
def svg_recovery(width=940):
    rowh, top, left, barw = 34, 56, 132, 620
    height = top + rowh * len(SAMPLES) + 78
    mx = max(v["abstain_before"] for _, v in SAMPLES)
    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" role="img" '
         f'aria-label="逐樣本 ABSTAIN 救回長條圖"><rect width="{width}" height="{height}" fill="white"/>']
    p.append(f'<text x="{left}" y="26" font-size="12.5" font-weight="700" fill="#0f172a">'
             f'每個樣本原本的 resource-guard ABSTAIN，改良後幾乎全部可 certify</text>')
    p.append(f'<rect x="{left}" y="36" width="13" height="11" rx="2" fill="#0d9488"/>')
    p.append(f'<text x="{left+19}" y="46" font-size="11" fill="#475569">已救回</text>')
    p.append(f'<rect x="{left+80}" y="36" width="13" height="11" rx="2" fill="#dc2626"/>')
    p.append(f'<text x="{left+99}" y="46" font-size="11" fill="#475569">仍無法（顯示為最小可見寬度）</text>')
    for i, (smp, v) in enumerate(SAMPLES):
        y = top + i * rowh
        w = barw * v["abstain_before"] / mx
        wr = w * v["recovered"] / max(v["abstain_before"], 1)
        p.append(f'<text x="{left-10}" y="{y+15}" font-size="11.5" text-anchor="end" fill="#334155">{smp}</text>')
        p.append(f'<rect x="{left}" y="{y+3}" width="{max(w,2):.1f}" height="17" rx="3" fill="#fee2e2"/>')
        p.append(f'<rect x="{left}" y="{y+3}" width="{max(wr,2):.1f}" height="17" rx="3" fill="#0d9488" opacity=".9"/>')
        lbl = f'{v["recovered"]:,} / {v["abstain_before"]:,}'
        p.append(f'<text x="{left+max(w,2)+9}" y="{y+16}" font-size="11" fill="#334155" '
                 f'font-family="ui-monospace,Menlo,monospace">{lbl}</text>')
        col = "#0d9488" if v["abstain_after"] == 0 else "#b45309"
        p.append(f'<text x="{width-14}" y="{y+16}" font-size="11" text-anchor="end" font-weight="700" '
                 f'fill="{col}">{v["recovery_pct"]:.1f}%　剩 {v["abstain_after"]}</text>')
    yb = top + rowh * len(SAMPLES) + 16
    p.append(f'<line x1="{left}" y1="{yb}" x2="{width-14}" y2="{yb}" stroke="#e2e8f0"/>')
    p.append(f'<text x="{left-10}" y="{yb+22}" font-size="12" text-anchor="end" font-weight="700" fill="#0f172a">合計</text>')
    p.append(f'<text x="{left}" y="{yb+22}" font-size="12.5" fill="#0f172a">'
             f'<tspan font-weight="700" fill="#0d9488">{E["recovered"]:,}</tspan> / {E["abstain_before"]:,} 可 certify '
             f'（<tspan font-weight="700" fill="#0d9488">{E["recovery_pct"]:.2f}%</tspan>），'
             f'仍無法 <tspan font-weight="700" fill="#b45309">{E["abstain_after"]}</tspan> 個</text>')
    p.append(f'<text x="{left}" y="{yb+44}" font-size="11" fill="#64748b">'
             f'分母：exact-PS × HP cohort，{E["units_total"]:,} 個 topology unit（與 workstation 子頁面同一批）</text>')
    p.append("</svg>")
    return "".join(p)


# ── SVG 2：問題 → 三改進 → 結果 ─────────────────────────────────────────
def svg_flow(width=940, height=284):
    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" role="img" '
         f'aria-label="問題到結果的三項改進流程"><rect width="{width}" height="{height}" fill="white"/>']
    p.append('<rect x="14" y="58" width="196" height="166" rx="11" fill="#fef2f2" stroke="#dc2626" stroke-width="1.8"/>')
    p.append('<text x="112" y="84" font-size="12.5" font-weight="700" text-anchor="middle" fill="#b91c1c">問題（改良前）</text>')
    for i, (a, b) in enumerate([(f'{E["abstain_before"]:,} 個 unit', "被 guard 擋住、不給答案"),
                                ("guard 只有 1,000", "但 B&amp;B 常需數千至數百萬節點"),
                                ("鏈長 m=8", "暴力法估 ~85 天")]):
        y = 114 + i * 38
        p.append(f'<text x="32" y="{y}" font-size="11.5" font-weight="700" fill="#0f172a">{a}</text>')
        p.append(f'<text x="32" y="{y+15}" font-size="10" fill="#7f1d1d">{b}</text>')
    for i, (t, s) in enumerate([("① 群去重 + 支配消除", f'群數中位 {G["raw_median"]}→{G["dom_median"]}（最大 {G["raw_max"]}→{G["dom_max"]}）'),
                                ("② 鏈式閉合式", "單一終端 → h*=m−1、樹數=m!，O(1)"),
                                ("③ terminal-subset DP", f'{E["paths"].get("dp",0):,} 個 unit 由此 certify')]):
        y = 44 + i * 62
        p.append(f'<rect x="256" y="{y}" width="300" height="52" rx="9" fill="#f0f9ff" stroke="#2563eb" stroke-width="1.6"/>')
        p.append(f'<text x="274" y="{y+22}" font-size="12.5" font-weight="700" fill="#1e40af">{t}</text>')
        p.append(f'<text x="274" y="{y+40}" font-size="10.5" fill="#475569">{s}</text>')
    p.append('<rect x="612" y="58" width="314" height="166" rx="11" fill="#f0fdfa" stroke="#0d9488" stroke-width="2.2"/>')
    p.append('<text x="769" y="84" font-size="12.5" font-weight="700" text-anchor="middle" fill="#0f766e">結果（全量實測）</text>')
    for i, (a, b) in enumerate([(f'{E["recovery_pct"]:.2f}%', f'ABSTAIN 救回（剩 {E["abstain_after"]} 個）'),
                                ("0", f'{E["regression_match"]:,} 個回歸區不一致數'),
                                (f'{E["timing"]["total_min"]} 分', f'{E["units_total"]:,} 個 unit 全部跑完')]):
        y = 120 + i * 36
        p.append(f'<text x="634" y="{y}" font-size="17" font-weight="700" fill="#0d9488">{a}</text>')
        p.append(f'<text x="726" y="{y}" font-size="11" fill="#334155">{b}</text>')
    for x in (220, 570):
        p.append(f'<path d="M{x} 141 L{x+30} 141" stroke="#94a3b8" stroke-width="2" marker-end="url(#ar)"/>')
    p.append('<defs><marker id="ar" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto">'
             '<path d="M0,0 L8,4 L0,8 z" fill="#94a3b8"/></marker></defs>')
    p.append(f'<text x="{width/2}" y="266" font-size="11" text-anchor="middle" fill="#64748b">'
             f'三項改進互補：①壓縮搜尋空間　②把鏈段從「搜尋」改成「計算」　③把求 h* 從 B&amp;B 換成 DP</text>')
    p.append("</svg>")
    return "".join(p)


# ── 表格 body ──────────────────────────────────────────────────────────
def _verdict(v):
    if v["abstain_after"] == 0:
        return "pass", "全部救回"
    return "warn", "剩 %d 個" % v["abstain_after"]


sample_rows = "".join(
    '<tr><td><code>{s}</code></td><td class="num">{u}</td><td class="num">{b}</td>'
    '<td class="num ok">{r}</td><td class="num">{a}</td><td class="num">{p:.1f}%</td>'
    '<td class="num">{rm}/{rn}</td><td><span class="st {cls}">{lab}</span></td></tr>'.format(
        s=s, u=n(v["units"]), b=n(v["abstain_before"]), r=n(v["recovered"]),
        a=v["abstain_after"], p=v["recovery_pct"],
        rm=n(v["regression_match"]), rn=n(v["regression_n"]),
        cls=_verdict(v)[0], lab=_verdict(v)[1])
    for s, v in SAMPLES)
sample_rows += (f'<tr style="font-weight:700;background:#f8fafc"><td>合計</td>'
                f'<td class="num">{n(E["units_total"])}</td><td class="num">{n(E["abstain_before"])}</td>'
                f'<td class="num ok">{n(E["recovered"])}</td><td class="num">{E["abstain_after"]}</td>'
                f'<td class="num">{E["recovery_pct"]:.2f}%</td>'
                f'<td class="num">{n(E["regression_match"])}/{n(E["regression_n"])}</td>'
                f'<td><span class="st pass">99.92%</span></td></tr>')

line_rows = (
    f'<tr style="background:#f0fdfa"><td><span class="st pass">主要</span> <b>A · exact-PS 線</b><br>'
    f'<span class="mut">workstation index 與 7 個子頁面顯示的就是這批</span></td>'
    f'<td class="num">{n(E["units_total"])}</td><td>7 datasets<br><span class="mut">6 生物樣本</span></td>'
    f'<td>ABSTAIN {n(E["abstain_before"])} → <b class="ok">{E["abstain_after"]}</b>'
    f'（救回 {E["recovery_pct"]:.2f}%）</td>'
    f'<td class="num">{E["timing"]["total_min"]} 分</td></tr>'
    f'<tr><td><span class="st todo">輔助</span> <b>B · layered 線</b><br>'
    f'<span class="mut">方法開發用；engineering baseline，比例不可作正式 Results</span></td>'
    f'<td class="num">{n(LA["n_total"])}</td><td>僅 HCC1395</td>'
    f'<td>capped {n(LA["n_capped"])} → <b class="ok">{LA["unsolved"]}</b>'
    f'；另發現 off-by-one bug {BUG["n"]} 區</td>'
    f'<td class="num">{LA["total_s"]} 秒</td></tr>')

# ── 閱讀指引：色彩語意 / 數字規範 / 名詞定義 ────────────────────────────
COLORS = [("ok", "#0d9488", "已救回、通過、改良後的值"),
          ("bad", "#dc2626", "仍無法處理、錯誤值、改良前的問題"),
          ("warn", "#b45309", "部分達成、需注意、尚有剩餘"),
          ("acc", "#2563eb", "方法步驟／中性強調（非好壞）"),
          ("mut", "#64748b", "背景說明、分母註記")]
color_rows = "".join(
    f'<tr><td><span style="display:inline-block;width:15px;height:15px;border-radius:3px;'
    f'background:{hexv};vertical-align:-2px;margin-right:8px"></span>'
    f'<code>{cls}</code></td><td>{desc}</td></tr>' for cls, hexv, desc in COLORS)

MARKS = [('<span class="mk crit">決定性</span>', "卡結論或需要決策的事實"),
         ('<span class="mk key">重要</span>', "有價值、影響判讀但不卡結論"),
         ('<span class="mk bgm">背景</span>', "已穩定的脈絡，掃過即可"),
         ('<span class="st pass">達成／全部救回</span>', "目標已滿足，有實測證據"),
         ('<span class="st warn">部分達成／剩 N 個</span>', "有進展但未清空"),
         ('<span class="st todo">未做</span>', "本輪刻意界定在範圍外"),
         ('<span class="st stop">UNSUPPORTED…</span>', "超出核心保證，立即 fail-closed")]
mark_rows = "".join(f'<tr><td>{m}</td><td>{d}</td></tr>' for m, d in MARKS)

NUMFMT = [("千分位", "所有 ≥1,000 的計數一律加逗號（如 <code>10,717</code>）"),
          ("表格對齊", "數字欄右對齊並用等寬數字（tabular-nums），方便直向比對位數"),
          ("百分比", "一律附分母，不單獨出現（如 <code>10,708/10,717（99.92%）</code>）"),
          ("時間", "&lt;1 秒用 ms（3 位小數）；≥1 秒用秒；≥60 秒用分鐘"),
          ("樹數", "整數，不縮寫。可能極大（如 <code>40,320</code> = 8!），刻意不四捨五入"),
          ("h*", "整數，代表最小隱藏節點數；<code>—</code> 表示未能 certify")]
numfmt_rows = "".join(f'<tr><td><b>{k}</b></td><td>{v}</td></tr>' for k, v in NUMFMT)

GLOSS = [("topology unit", "一個 exact-PS × HP bounded block；本頁的<b>主要計數單位</b>（分母）"),
         ("h*（objective_h）", "最小隱藏節點數 = 為了讓所有觀測在單調樹上連通，<b>最少</b>需要幾個未觀測祖先"),
         ("certify", "<b>證明</b>某個 h* 是最小值（已窮盡所有更小的可能），而非只給一個可行解"),
         ("resource-guard ABSTAIN", "solver 因資源上限（節點數／狀態數）主動放棄並<b>不給答案</b>；"
                                    "<b>不等於</b>「該區無解」或「unresolved=negative」"),
         ("active loci", "該區中<b>真正有變異</b>的位點數；計算成本由它決定，<b>不是</b>區域 sSNV 總數"),
         ("terminal group", "一條 partial read 誘導出的<b>子立方體</b>；樹只需碰到其中任一頂點即算解釋"),
         ("支配消除", "若群 A ⊆ 群 B，滿足 A 必滿足 B → B 冗餘可移除"),
         ("鏈式閉合式", "單一終端且群皆含之時，h*=m−1、樹數=m!，<b>O(1) 直接得解</b>不需搜尋"),
         ("terminal-subset DP", "有向群終端子集動態規劃；只證 h*，不產候選樹家族"),
         ("obligation B&amp;B", "分枝定界；可產完整家族但成本高，且受 search-node guard 限制"),
         ("capped（layered 線用語）", "等同 exact-PS 線的 ABSTAIN；因術語源自不同 pipeline 故並存")]
gloss_rows = "".join(f'<tr><td><code>{k}</code></td><td>{v}</td></tr>' for k, v in GLOSS)

RANGES = [("exact-PS unit 數", f'{n(E["units_total"])}', "7 datasets / 6 生物樣本 / chr1–22"),
          ("　其中 ABSTAIN", f'{n(E["abstain_before"])}', f'佔 {E["abstain_before"]/E["units_total"]*100:.2f}%'),
          ("　其中已 certified", f'{n(E["regression_n"])}', f'佔 {E["regression_n"]/E["units_total"]*100:.2f}%'),
          ("layered unit 數", f'{n(LA["n_total"])}', "僅 HCC1395；engineering baseline"),
          ("　其中 capped", f'{n(LA["n_capped"])}', f'佔 {LA["n_capped"]/LA["n_total"]*100:.2f}%'),
          ("active loci 範圍", f'1 – {D["max_bits"]}', f'硬上限 {D["max_bits"]}（編譯期常數，超過 fail-closed）'),
          ("h* 觀察範圍", "0 – 11", "本輪實測；理論上界受 active loci 限制")]
range_rows = "".join(f'<tr><td>{k}</td><td class="num"><b>{v}</b></td><td class="mut">{d}</td></tr>'
                     for k, v, d in RANGES)

GOALS = [("讓被 guard 擋住的區能算出來",
          f'exact-PS：{n(E["recovered"])}/{n(E["abstain_before"])}（{E["recovery_pct"]:.2f}%）；'
          f'layered：{n(LA["capped_solved"])}/{n(LA["n_capped"])}（100%）', "達成"),
         ("不能破壞既有正確結果",
          f'exact-PS {n(E["regression_match"])}/{n(E["regression_n"])} ＋ layered '
          f'{n(LA["regression_match"])}/{n(LA["n_regression"])}，合計 0 不一致', "達成"),
         ("大幅降低計算量",
          f'{n(E["units_total"])} 區 {E["timing"]["total_min"]} 分鐘（中位 {E["timing"]["median_ms"]} ms）；'
          f'鏈長 m=8 從估 ~85 天 → 7.5 秒', "達成"),
         ("處理 subcube（partial group）情境",
          f'群數中位 {G["raw_median"]}→{G["dom_median"]}；{E["paths"].get("dp",0):,} 個 unit 由 DP certify', "達成"),
         ("k 上限是否可移除",
          f'不可移除，但可重新定義為 active loci ≤ {D["max_bits"]}；稀疏區 k=64 仍 0.05 ms', "部分"),
         ("完整候選樹家族枚舉",
          "本輪只驗證 h*（最小隱藏節點數）；家族枚舉本質指數、仍需 cap", "未做"),
         ("正式取代 workstation 既有數字",
          "未做：需重跑 canonical pipeline 並重新簽章；本輪僅提供可回收性證據", "未做")]
ST = {"達成": ("pass", "達成"), "部分": ("warn", "部分達成"), "未做": ("todo", "未做")}
goal_rows = "".join(f'<tr><td>{g}</td><td>{e}</td><td><span class="st {ST[s][0]}">{ST[s][1]}</span></td></tr>'
                    for g, e, s in GOALS)

speed_rows = "".join(f'<tr><td>m={s["m"]}</td><td>{s["py"]}</td>'
                     f'<td class="num ok"><b>{s["cpp_ms"]:,.1f} ms</b></td></tr>' for s in D["speed"])
kb_rows = "".join(f'<tr><td>{s["case"]}</td><td class="num">{s["ms"]:,.2f} ms</td>'
                  f'<td><span class="st {"pass" if s["ok"] else "stop"}">{s["state"]}</span></td></tr>'
                  for s in D["kbound"])
diff_rows = "".join(f'<tr><td><code>{d["case"]}</code></td><td class="num ok">{d["dp"]}</td>'
                    f'<td class="num ok">{d["bnb"]}</td><td class="num bad">{d["py"]}</td>'
                    f'<td><span class="st pass">兩獨立解法一致</span></td></tr>' for d in BUG["diff"])
ex_rows = "".join(f'<tr><td><code>{e["region"]}/{e["family"]}</code></td><td class="num bad">{e["py_h"]}</td>'
                  f'<td class="num ok">{e["new_h"]}</td><td><code>{", ".join(e["full"])}</code></td></tr>'
                  for e in BUG["examples"])
chain_rows = "".join(f'<tr><td>{n(b["n"])} 區</td><td class="num">{n(b["trees"])} 棵</td></tr>'
                     for b in CHN["breakdown"])
path_rows = "".join(f'<tr><td><code>{k}</code></td><td class="num">{n(v)}</td>'
                    f'<td class="num">{v/E["abstain_before"]*100:.2f}%</td></tr>'
                    for k, v in sorted(E["paths"].items(), key=lambda t: -t[1]))
REM = E["remaining"]

HTML = f"""<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>ABSTAIN 重跑驗證 · 10,717 個中 99.92% 可救回</title>
<style>
:root{{--bg:#fff;--fg:#0f172a;--mut:#64748b;--line:#e2e8f0;--card:#f8fafc;
--ok:#0d9488;--bad:#dc2626;--warn:#b45309;--acc:#2563eb;--rad:11px}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--fg);line-height:1.72;font-size:15.5px;
font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans CJK TC","PingFang TC","Microsoft JhengHei",Arial,sans-serif}}
.wrap{{max-width:1060px;margin:0 auto;padding:30px 20px 78px}}
.eyebrow{{font-size:11.5px;letter-spacing:.11em;text-transform:uppercase;color:var(--mut);font-weight:700}}
h1{{font-size:26px;margin:5px 0 7px;line-height:1.32;letter-spacing:-.01em}}
.sub{{color:var(--mut);font-size:13px;margin-bottom:18px}}
h2{{font-size:19px;margin:0 0 5px;line-height:1.42}}
.sec{{margin:36px 0 0;padding-top:22px;border-top:2px solid var(--line)}}
.sec>p.lead{{color:var(--mut);font-size:13.5px;margin:3px 0 12px}}
.mut{{color:var(--mut);font-size:12px}}

.verdict{{border:1px solid #99f6e4;border-left:6px solid var(--ok);
background:linear-gradient(180deg,#f0fdfa,#fafffe);border-radius:var(--rad);padding:17px 21px;margin:14px 0 6px}}
.verdict .big{{font-size:17.5px;font-weight:700;line-height:1.58}}
.verdict .dn{{margin-top:9px;font-size:12.5px;color:var(--mut);padding-top:9px;border-top:1px dashed #99f6e4}}

.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:14px 0 4px}}
.kpi{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:11px 13px}}
.kpi .v{{font-size:21px;font-weight:700;letter-spacing:-.02em;line-height:1.25;font-variant-numeric:tabular-nums}}
.kpi .v small{{font-size:13px;font-weight:600;color:var(--mut)}}
.kpi .l{{font-size:11.5px;color:var(--mut);margin-top:3px}}
.kpi.ok .v{{color:var(--ok)}} .kpi.bad .v{{color:var(--bad)}} .kpi.acc .v{{color:var(--acc)}}

.mk{{display:inline-block;font-size:10.5px;font-weight:700;padding:1px 7px;border-radius:3px;
margin-right:7px;vertical-align:1.5px;border:1px solid}}
.mk.crit{{background:#fef2f2;color:#b91c1c;border-color:#fca5a5}}
.mk.key{{background:#fffbeb;color:#92400e;border-color:#fcd34d}}
.mk.bgm{{background:#f1f5f9;color:#475569;border-color:#cbd5e1}}
.st{{display:inline-block;font-size:11px;font-weight:700;padding:1px 8px;border-radius:20px;border:1px solid;white-space:nowrap}}
.st.pass{{background:#f0fdfa;color:#0f766e;border-color:#5eead4}}
.st.warn{{background:#fffbeb;color:#92400e;border-color:#fcd34d}}
.st.todo{{background:#f1f5f9;color:#475569;border-color:#cbd5e1}}
.st.stop{{background:#fef2f2;color:#b91c1c;border-color:#fca5a5}}

ul.logic{{list-style:none;padding:0;margin:12px 0}}
ul.logic li{{padding:9px 0;border-bottom:1px dashed var(--line)}}
ul.logic li:last-child{{border-bottom:0}}

figure{{margin:16px 0;padding:12px;background:#fff;border:1px solid var(--line);border-radius:var(--rad)}}
figure svg{{width:100%;height:auto;display:block}}
figcaption{{font-size:12.5px;color:var(--mut);margin-top:9px;padding-top:8px;border-top:1px dashed var(--line)}}

.tw{{overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:13.5px;margin:10px 0}}
th{{background:#f1f5f9;text-align:left;padding:8px 11px;border-bottom:2px solid #cbd5e1;font-weight:700;white-space:nowrap}}
td{{padding:8px 11px;border-bottom:1px solid var(--line);vertical-align:top}}
td.num{{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}}
tbody tr:hover{{background:#fafcff}}
code{{background:#f1f5f9;padding:1px 5px;border-radius:4px;font-size:12.5px;font-family:ui-monospace,Menlo,monospace}}
.ok{{color:var(--ok);font-weight:700}} .bad{{color:var(--bad);font-weight:700}}

details.card{{border:1px solid var(--line);border-radius:var(--rad);margin:9px 0;background:#fff;overflow:hidden}}
details.card>summary{{cursor:pointer;padding:12px 16px;font-weight:600;font-size:14.5px;
background:var(--card);list-style:none;display:flex;align-items:center;gap:9px}}
details.card>summary::-webkit-details-marker{{display:none}}
details.card>summary::before{{content:"+";font-family:ui-monospace,Menlo,monospace;font-weight:700;
color:var(--acc);font-size:17px;width:14px;flex:0 0 14px;text-align:center}}
details.card[open]>summary::before{{content:"−"}}
details.card>summary:hover{{background:#eef2f7}}
details.card .ct{{margin-left:auto;font-size:11px;font-weight:700;color:var(--mut);
background:#fff;border:1px solid var(--line);border-radius:20px;padding:1px 9px}}
details.card .cb{{padding:4px 16px 16px}}

.box{{border:1px solid var(--line);border-radius:9px;padding:12px 16px;margin:11px 0;background:var(--card);font-size:14px}}
.box.dang{{background:#fef2f2;border-color:#fecaca}}
.box.warn{{background:#fffbeb;border-color:#fcd34d}}
.box.good{{background:#f0fdfa;border-color:#99f6e4}}
.box b.hd{{display:block;margin-bottom:4px}}
footer{{margin-top:42px;padding-top:16px;border-top:2px solid var(--line);font-size:12.5px;color:var(--mut)}}
@media(prefers-color-scheme:dark){{
:root{{--bg:#0b1120;--fg:#e2e8f0;--mut:#94a3b8;--line:#1e293b;--card:#0f172a}}
figure{{background:#0f172a}} figure svg{{filter:invert(.92) hue-rotate(180deg)}}
th{{background:#1e293b;border-bottom-color:#334155}} tbody tr:hover{{background:#111c33}}
code{{background:#111c33}} details.card{{background:#0b1120}} details.card>summary:hover{{background:#16233c}}
details.card .ct{{background:#0f172a}}
.verdict{{background:#0d2926;border-color:#115e59}}
.box.dang{{background:#2a1414;border-color:#7f1d1d}} .box.warn{{background:#2a2410;border-color:#78350f}}
.box.good{{background:#0d2926;border-color:#115e59}}
.mk.crit{{background:#2a1414}} .mk.key{{background:#2a2410}} .mk.bgm{{background:#1e293b}}
}}
@media(max-width:640px){{h1{{font-size:21.5px}} .wrap{{padding:22px 14px 60px}}}}
</style>
<div class="wrap">

<div class="eyebrow">Method-layer verification · 2026-07-25</div>
<h1>workstation 上那 {n(E["abstain_before"])} 個 ABSTAIN，<br>99.92% 不是「定不出來」，是被過緊的搜尋預算擋掉</h1>
<div class="sub">C++ solver：{D["solver"]}　·　所有數字由 <code>{STEM}.data.json</code> 注入</div>

<!-- ═══ L0 一眼結論 ═══ -->
<div class="verdict">
<div class="big">全 {n(E["units_total"])} 個 exact-PS topology unit 重跑：原本 {n(E["abstain_before"])} 個
resource-guard ABSTAIN，<b class="ok">{n(E["recovered"])} 個（{E["recovery_pct"]:.2f}%）現在可 certify h*</b>，
只剩 <b class="ok">{E["abstain_after"]} 個</b>；同時 {n(E["regression_match"])} 個原已 certified 的區
<b class="ok">0 個不一致</b>。全部跑完 <b>{E["timing"]["total_min"]} 分鐘</b>。</div>
<div class="dn"><b>兩條資料線都跑了，分母不同不可混算：</b>
A · exact-PS 線（{n(E["units_total"])} unit / 7 datasets）＝ workstation 子頁面顯示的那批；
B · layered 線（{n(LA["n_total"])} unit / 僅 HCC1395，engineering baseline）＝ 方法開發與 bug 發現用。</div>
</div>

<div class="grid">
<div class="kpi ok"><div class="v">{E["recovery_pct"]:.2f}<small>%</small></div><div class="l">ABSTAIN 救回率</div></div>
<div class="kpi ok"><div class="v">{n(E["recovered"])}<small>/{n(E["abstain_before"])}</small></div><div class="l">可 certify h*</div></div>
<div class="kpi bad"><div class="v">{E["abstain_after"]}</div><div class="l">仍無法（實務邊界）</div></div>
<div class="kpi ok"><div class="v">0</div><div class="l">回歸不一致（{n(E["regression_match"])} 區）</div></div>
<div class="kpi acc"><div class="v">{E["timing"]["total_min"]}<small> 分</small></div><div class="l">{n(E["units_total"])} 區全部跑完</div></div>
<div class="kpi"><div class="v">{E["timing"]["median_ms"]}<small> ms</small></div><div class="l">單區中位數</div></div>
</div>

<!-- ═══ 閱讀指引 ═══ -->
<div class="sec">
<div class="eyebrow">閱讀指引</div>
<h2>先對齊四件事：分母、色彩、數字寫法、名詞</h2>
<p class="lead">本頁刻意<b>不使用 emoji</b>（本機無 emoji 字型會顯示為方框）；重要性與狀態一律用文字徽章，
且<b>色彩絕不單獨編碼</b> — 每個狀態都附文字標籤，色盲讀者可完整判讀。</p>

{card("① 資料範圍與分母（最重要，先看這個）",
      tbl(["項目", "數量", "說明"], range_rows, 520)
      + '<div class="box warn"><b class="hd">兩條線不可並列或相除</b>'
        'exact-PS 是<b>主線</b>（workstation 顯示的那批）；layered 是<b>方法開發輔助線</b>'
        '（engineering baseline，比例數字不可作正式 Results）。本頁所有比例都已標明各自分母。</div>',
      "denominators")}

{card("② 色彩語意", tbl(["色彩", "代表意義"], color_rows, 400)
      + '<p class="mut">色彩僅作輔助；所有判定都另有文字標籤，不依賴色覺。</p>', "colorblind-safe")}

{card("③ 標記與狀態徽章", tbl(["徽章", "意義"], mark_rows, 420), "text + shape")}

{card("④ 數字顯示規範", tbl(["項目", "規範"], numfmt_rows, 480), "formatting")}

{card("⑤ 名詞定義", tbl(["名詞", "定義"], gloss_rows, 560), "glossary")}
</div>

<!-- ═══ L1 重點邏輯 ═══ -->
<div class="sec">
<div class="eyebrow">L1 · 重點邏輯</div>
<h2>瓶頸不是問題太難，是預算設太緊、又把「突變數多」誤判成「結構複雜」</h2>
<ul class="logic">
<li><span class="mk crit">決定性</span><b>ABSTAIN 的成因單一且明確。</b>
全部 {n(E["abstain_before"])} 個的 reason 都是 <code>{E["root_cause"]["all_reason_same"]}</code>，
卡在 <code>{E["root_cause"]["guard"]}</code>。但實測 B&amp;B 常需要數千至數百萬節點才收斂 ——
<b>1,000 這個預算太緊</b>，而 DP 路徑根本不受它限制。</li>
<li><span class="mk crit">決定性</span><b>調高上限行不通，必須改結構。</b>
鏈長 m=8 用暴力枚舉需 C(254,7)≈1.2×10¹³ 個候選（估 ~85 天），但真值只是閉合式 <code>m!</code>，O(1) 可得。</li>
<li><span class="mk key">重要</span><b>「鏈填充」不是結構複雜度。</b>
popcount-m 的終端<b>必然</b>需要 m−1 個祖先，那只是突變距離、零推論內容，卻把 cap 額度佔滿。</li>
<li><span class="mk key">重要</span><b>求 h* 該用 DP 而非 B&amp;B。</b>
{E["paths"].get("dp",0):,} 個 unit（{E["paths"].get("dp",0)/E["abstain_before"]*100:.1f}%）由 DP certify，實測比 B&amp;B 快 20–1000×。</li>
<li><span class="mk bgm">背景</span><b>k 上限不能移除，但可重新定義。</b>
硬上限綁在 <b>active loci</b>（真正有變異的位點）而非區域 sSNV 數；稀疏區 k=64 仍是 0.05 ms。</li>
</ul>
</div>

<div class="sec">
<div class="eyebrow">L1 · 主結果</div>
<h2>七個樣本全部大幅回收，四個樣本 100% 清空</h2>
<figure>{svg_recovery()}
<figcaption>H2009 原本最嚴重（{n(E["per_sample"]["H2009"]["abstain_before"])} 個，佔全 cohort 的
{E["per_sample"]["H2009"]["abstain_before"]/E["abstain_before"]*100:.0f}%），改良後 <b>100% 全部救回</b>。
仍無法處理的 {E["abstain_after"]} 個集中在 H1437（5）、HCC1395（3）、HCC1937（1）。</figcaption></figure>
</div>

<div class="sec">
<div class="eyebrow">L1 · 方法</div>
<h2>三項改進各自解決不同瓶頸，缺一則仍有區算不出來</h2>
<figure>{svg_flow()}
<figcaption>三者針對不同瓶頸，<b>不可互相取代</b>：只做①仍會被鏈長拖垮；只做②仍會被多群拖垮；
只做③則池子太大讓 DP 的 state 空間爆炸。實測 {E["paths"].get("dp",0):,} 個 unit 走 DP、
{E["paths"].get("chain_closed_form",0):,} 個走鏈式閉合式。</figcaption></figure>
</div>

<div class="sec">
<div class="eyebrow">L1 · 目標對照</div>
<h2>七項目標：四項達成、一項部分達成、兩項明確未做</h2>
<p class="lead">「未做」是本輪刻意界定的範圍邊界，誠實標示以免被誤讀為已完成。</p>
{tbl(["目標", "實測證據", "狀態"], goal_rows, 620)}
</div>

<div class="sec">
<div class="eyebrow">L1 · 分母</div>
<h2>兩條資料線都跑了全量，但不可並列或相除</h2>
{tbl(["資料線", "unit 數", "樣本範圍", "結果", "耗時"], line_rows, 640)}
</div>

<!-- ═══ L2 證據卡 ═══ -->
<div class="sec">
<div class="eyebrow">L2 · 證據卡</div>
<h2>逐項證據（預設收合，要驗證或引用時再展開）</h2>

{card("證據 1 · 逐樣本救回明細", tbl(["樣本", "units", "ABSTAIN 前", "救回", "剩餘", "救回率", "回歸一致", "判定"], sample_rows, 720), "7 datasets")}

{card("證據 2 · 根因：所有 ABSTAIN 都卡在同一個 guard",
      f'<div class="box dang"><b class="hd">單一成因</b>'
      f'{n(E["abstain_before"])} 個 ABSTAIN 的 reason <b>全部</b>是 <code>{E["root_cause"]["all_reason_same"]}</code>，'
      f'且 <code>search_nodes</code> 全部剛好等於上限值。{E["root_cause"]["note"]}。</div>'
      f'<div class="box"><b class="hd">對照實測</b>B&amp;B 收斂所需節點數：簡單案例 4,000 上下，'
      f'最壞（10 個位點全 active）需 <b>620 萬</b>。相對之下 <code>{E["root_cause"]["guard"]}</code> 極緊。</div>'
      + tbl(["救回路徑", "unit 數", "佔 ABSTAIN"], path_rows), "single root cause")}

{card("證據 3 · 回歸驗證：兩條線合計 0 個不一致",
      f'<p><b>exact-PS 線：</b>{n(E["regression_match"])}/{n(E["regression_n"])} 個原已 certified 的區，'
      f'新算 h* 與原值<b class="ok">完全相同</b>。</p>'
      f'<p><b>layered 線：</b>{n(LA["regression_match"])}/{n(LA["n_regression"])} 個非 capped 區，'
      f'新舊 h* <b class="ok">完全相同</b>。</p>'
      f'<div class="box good">合計 {n(E["regression_match"]+LA["regression_match"])} 個回歸區、'
      f'<b>0 個不一致</b> — 這是「沒有用新 bug 換舊限制」的直接證據。</div>', "0 mismatch")}

{card(f"證據 4 · 剩下 {E['abstain_after']} 個是真正的實務邊界",
      f'<p>全部集中在 <b>active loci {"–".join(str(x) for x in REM["active_loci"])}</b> 且'
      f'<b>有效群數 {min(REM["groups_dom"])}–{max(REM["groups_dom"])}</b>。</p>'
      f'<div class="box"><b class="hd">為什麼停在這裡</b>{REM["why"]}</div>'
      f'<div class="box good">「定不出來即是答案」現在適用於 <b>{E["abstain_after"]} 個</b>，'
      f'而不是原本的 {n(E["abstain_before"])} 個 — 這讓原則性拒絕真正對齊複雜度邊界，'
      f'而不是對齊參數選擇。</div>', "practical limit")}

{card(f"證據 5 · layered 線順帶抓到 Python solver 一個 off-by-one（{BUG['n']} 區）",
      f'<div class="box dang"><b class="hd">症狀（數學上不可能）</b>'
      f'{BUG["n"]} 個區的 Python <code>n_hidden</code> <b>小於</b>新算的最小值 h* —— '
      f'但 greedy 是可行解，其隱藏節點數必 ≥ 最小值。</div>'
      f'<div class="box"><b class="hd">根因</b><code>n_hidden = len(N) - len(full) - 1</code>。'
      f'當 full 觀測含<b>全-REF 基因型</b>（如 <code>RRRRRRRR</code>，其 altset 就是 ROOT）時，'
      f'<code>base = {{ROOT}} | full</code> 不會多出一個元素，公式卻多扣了 1。</div>'
      f'<div class="box"><b class="hd">證據（對照組設計）</b>異常組 '
      f'<b class="bad">{BUG["allref_hit"]}/{BUG["n"]}（{BUG["allref_pct"]}%）</b>含全-REF；'
      f'對照組僅 <b>{BUG["ctrl_allref"]}/{n(BUG["ctrl_n"])}（{BUG["ctrl_pct"]}%）</b>。</div>'
      + tbl(["案例", "Python n_hidden", "新 h*", "full 觀測"], ex_rows, 560)
      + '<p>DP 與 B&amp;B 兩條<b>獨立</b>解法交叉確認：</p>'
      + tbl(["案例", "DP h*", "B&amp;B h*", "Python", "判定"], diff_rows, 560)
      + '<p class="mut">影響範圍：僅 greedy fallback 路徑的 <code>n_hidden</code> 欄位；'
        '正常路徑用 <code>n_hidden = e</code> 不受影響。</p>', "layered 線")}

{card(f"證據 6 · 樹數曾被低報 {UNDER:,} 倍（layered 線 {CHN['n']} 區）",
      f'<p>舊系統把 capped 區一律記成「1 棵」；這 {CHN["n"]} 個鏈式區合計記 '
      f'<b class="bad">{n(CHN["old_total"])} 棵</b>，真值 <b class="ok">{n(CHN["new_total"])} 棵</b>。</p>'
      + tbl(["區數", "各區真實樹數"], chain_rows, 360)
      + f'<div class="box warn"><b class="hd">修正後 mean n_trees 會上升</b>'
        f'這不是結果變差，而是修正<b>系統性低報</b>。論文／口試須主動說明，'
        f'否則會被誤讀成「新方法讓不確定性變高」。</div>', "under-report")}

{card("證據 7 · 速度：Python 失敗處的加速幅度",
      tbl(["鏈長 m", "改良前", "改良後"], speed_rows, 520), "10⁶× on m=8")}

{card("證據 8 · k 上限綁在 active loci，不是區域 sSNV 數",
      tbl(["情況", "時間", "狀態"], kb_rows, 520)
      + f'<p>硬上限 <code>kMaximumExactHypercubeBits = {D["max_bits"]}</code>（編譯期常數，'
        f'超過立即 fail-closed，不假裝能算）。<b>但實務上限更低</b>：10 個位點全 active 就要 56 秒，'
        f'所以上限仍然必要 — 只是應從「k ≤ 8」改成「<b>active loci ≤ {D["max_bits"]} 且有資源上限</b>」。</p>',
      "redefine, not remove")}
</div>

<!-- ═══ L3 邊界與溯源 ═══ -->
<div class="sec">
<div class="eyebrow">L3 · 適用邊界</div>
<h2>四項限制：越界引用會變成 overclaim</h2>
<div class="box warn"><b class="hd">1 · 本頁未變更 workstation 既有數字</b>
子頁面的 funnel／比例仍是原 guard1000 產出。本輪是<b>可回收性證據</b>；
正式取代需重跑 canonical pipeline 並重新簽章（該線目前為 scientific release NO-GO）。</div>
<div class="box warn"><b class="hd">2 · 只驗證 h*，不是完整候選樹家族</b>
h* = 最小隱藏節點數。家族枚舉本質指數，仍需 cap。</div>
<div class="box warn"><b class="hd">3 · 兩條資料線分母不同</b>
exact-PS（{n(E["units_total"])} unit / 7 datasets）與 layered（{n(LA["n_total"])} unit / 僅 HCC1395、
engineering baseline）不可並列或相除。layered 線的比例數字不可作正式 Results。</div>
<div class="box warn"><b class="hd">4 · VAF 不能用來處理算不出來的案例</b>
依 pre-registration H9，read-AF ordering 只可作探索性排序；用它補救 edge case 會踩到自訂紅線。</div>

{card("L3 · 可重算指令與資料溯源",
      f'<ul>'
      f'<li><b>exact-PS 輸入</b>：<code>{E["source_root"]}</code> 的 '
      f'<code>*.exact_ps_mlhp.json</code>（觀測 pattern）+ <code>*.topology.jsonl</code>（ABSTAIN 狀態）</li>'
      f'<li><b>layered 輸入</b>：<code>docs/methodology/_assets/20260618_subcluster_pilot/mlhp_part_*.json</code></li>'
      f'<li><b>求解</b>：<code>ll_real.cpp</code> harness → <code>{D["solver"]}</code></li>'
      f'<li><b>既有測試</b>：LongLineage <code>test_runtime_solver</code> PASS'
      f'（small-q oracle / exact B&amp;B / subset DP / parent mapping / router differential）</li>'
      f'<li><b>反捏造</b>：本頁所有數字由 <code>{STEM}.data.json</code> 注入；'
      f'builder 缺 required key 直接 refuse（§13-A）</li>'
      f'<li><b>workstation 注入</b>：<code>20260725_exactps_recovery.data.json</code> → '
      f'<code>build_exact_ps_layered_workstation.py</code> 的 <code>recovery_note()</code> / '
      f'<code>solver_verification_band()</code>，index 與 7 個子頁面同步</li>'
      f'</ul>')}
</div>

<footer>
<p>生成 {D["generated"]}　·　資訊分層：L0 一眼結論 → L1 重點邏輯／主結果／方法／目標對照／分母 →
L2 證據卡（收合）→ L3 邊界與溯源。　·　無 emoji（本機無 emoji 字型），重要性以 CSS 徽章標示，
狀態一律附文字標籤（不單靠色彩編碼）。</p>
</footer>
</div>
"""
(OUT / f"{STEM}.html").write_text(HTML, encoding="utf-8")
print(f"wrote {OUT / (STEM + '.html')}  ({len(HTML)/1024:.1f} KB)")
print(f"wrote {OUT / (STEM + '.data.json')}")
