#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ⛔ DEPRECATED 2026-07-10 — 舊 is_somatic 骨幹。canonical driver 改 build_layered_per_sample.py
#    (ClairS v0.4.0 paired PASS + 分層樹枚舉 + 全面板)。此 driver 僅供歷史重生 topology_workstation/。
"""
每樣本工作站 driver（2026-07-06）— 把 66MB 單檔拆成「每樣本一檔(~9MB) + 主頁 index.html」。

背景:單檔 20260629_multisample_topology_workstation.standalone.html 已達 63MB(99.96% 是 7 樣本
inline JSON)→ 換到低 RAM 電腦瀏覽器 parse 不動/白畫面、傳輸易被 25-50MB 上限截斷。拆每樣本後各 ~9MB
任何機器可開。

做法:複用 build_topology_workstation.py 的 SM_SAMPLES(子集)+ SM_OUT(輸出路徑)env,對 7 樣本各跑一次
→ topology_workstation/{sample}.html;再產 topology_workstation/index.html(主頁:跨樣本摘要表 + 開啟連結)。
不改 build_topology_workstation.py。

用法:python3 build_per_sample_workstation.py     # 重生整個 topology_workstation/ 資料夾
"""
import os, sys, json, subprocess, gzip, shutil
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "build_topology_workstation.py")
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "topology_workstation"))
PY = sys.executable


def sample_dirs():
    pairs = [("HCC1395", DATA)]
    if os.path.isdir(MSROOT):
        for s in sorted(os.listdir(MSROOT)):
            if os.path.exists(os.path.join(MSROOT, s, "topology_per_region.json")):
                pairs.append((s, os.path.join(MSROOT, s)))
    return pairs


CONF = {"full_tree", "linear_nested", "sibling_only", "co_linked_lineage"}


def _canon(edges):
    """AHU 樹同構標準式(只看形狀不看標籤);有環回 '(X)'。"""
    from collections import defaultdict as dd
    ch = dd(list)
    for p, c in edges:
        ch[p].append(c)
    seen = set()

    def rec(nd):
        if nd in seen:
            return "(X)"
        seen.add(nd)
        return "(" + "".join(sorted(rec(x) for x in ch.get(nd, []))) + ")"
    return rec("ROOT")


def sample_summary(dr):
    """從 topology_per_region.json 算主頁多維度比較數字(fresh,不抄中間檔)。"""
    d = json.load(open(os.path.join(dr, "topology_per_region.json"), encoding="utf-8"))["detail"]
    n = len(d)
    inc = sum(1 for r in d if r.get("determinacy") == "incompatible")
    c = Counter(r["n_clusters"] for r in d)
    c2 = [r for r in d if r["n_clusters"] == 2]
    n2 = len(c2) or 1
    br = sum(1 for r in c2 if "branched" in (r.get("topology_type") or ""))
    lin = sum(1 for r in c2 if "linear" in (r.get("topology_type") or ""))
    cf = sum(1 for r in c2 if r.get("tree_shape") in CONF)
    # canonical 拓撲組成(單群/2平行/2直系)+ 不同形狀種類
    shapes = Counter()
    for r in d:
        e = r.get("edges") or []
        if e:
            shapes[_canon(e)] += 1
    return {"n": n, "inc_pct": round(100 * inc / n, 1), "inc_n": inc,
            "c1": c.get(1, 0), "c2": len(c2), "c3p": sum(v for k, v in c.items() if k >= 3),
            "br_pct": round(100 * br / n2, 1), "conf_pct": round(100 * cf / n2, 1),
            "br_n": br, "lin_n": lin,
            "single_pct": round(100 * shapes.get("(())", 0) / n, 1),
            "par2_pct": round(100 * shapes.get("(()())", 0) / n, 1),
            "lin2_pct": round(100 * shapes.get("((()))", 0) / n, 1),
            "nshapes": len(shapes)}


def build_one(name, dr):
    out = os.path.join(OUT, name + ".html")
    env = dict(os.environ, SM_SAMPLES="{}:{}".format(name, dr), SM_OUT=out)
    r = subprocess.run([PY, BUILD], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    ok = r.returncode == 0 and os.path.exists(out)
    size = os.path.getsize(out) if os.path.exists(out) else 0
    return {"ok": ok, "path": out, "size": size, "log": r.stdout.decode("utf-8", "replace")[-300:]}


def human(b):
    return "{:.1f} MB".format(b / 1048576) if b >= 1048576 else "{:.0f} KB".format(b / 1024)


TI = {"HCC1395": "SEQC2 truth·凍結主樣本", "COLO829": "⚠低coread artifact", "HCC1954": "⚠undetermined",
      "H2009": "⚠資料品質(incompatible 最高)", "HCC1937": "census 標籤弱", "H1437": "census 標籤弱",
      "HCC1395_DORADO": "同細胞株·DORADO basecaller"}


def cramerv(rows):
    """c=2 branched-vs-linear × 樣本 列聯表 → CramérV(跨樣本效應量)。"""
    obs = [[x["sum"]["br_n"], x["sum"]["lin_n"]] for x in rows if (x["sum"]["br_n"] + x["sum"]["lin_n"]) > 0]
    try:
        from scipy.stats import chi2_contingency
        import numpy as np
        M = np.array(obs)
        chi2, pv, dof, _ = chi2_contingency(M)
        return round((chi2 / (M.sum() * (min(M.shape) - 1))) ** 0.5, 3)
    except Exception:
        return None


def _cbar(s):  # canonical 組成 stacked mini-bar(單群灰/2平行橙/2直系藍)
    return ('<div style="display:flex;height:11px;width:100px;border-radius:2px;overflow:hidden;background:#f1f3f5" '
            'title="單群 {sg}% / 2平行 {p2}% / 2直系 {l2}%"><div style="width:{sg}%;background:#adb5bd"></div>'
            '<div style="width:{p2}%;background:#f59f00"></div><div style="width:{l2}%;background:#1c7ed6"></div></div>'
            '<div class="note" style="font-size:9.5px">{sg}/{p2}/{l2}%</div>').format(
        sg=s["single_pct"], p2=s["par2_pct"], l2=s["lin2_pct"])


def _bcbar(s):  # branched(幾何上界)橙 vs confirmed(read)綠 paired bar
    row = ('<div style="display:flex;align-items:center;gap:3px"><span style="width:34px;font-size:10px;color:{col}">{v}%</span>'
           '<div style="width:62px;background:#f1f3f5;border-radius:2px"><div style="width:{v}%;height:7px;background:{bg}"></div></div></div>')
    return (row.format(v=s["br_pct"], col="#f08c00", bg="#f59f00") +
            row.format(v=s["conf_pct"], col="#2b8a3e", bg="#37b24d"))


def index_html(rows, cv):
    tr = ""
    for x in rows:
        s = x["sum"]
        cp = lambda k: round(100 * s[k] / (s["n"] or 1))
        tr += ('<tr><td><b>{n}</b><div class="note">舊骨幹單樣本快照</div></td>'
               '<td>{nn}</td><td class="note">{p1}/{p2}/{p3}%</td><td>{cbar}</td><td>{nsh}</td>'
               '<td{ic}>{inc}%</td><td>{bcbar}</td><td class="note">{sz}</td>'
               '<td><a class="snapshot-link" href="{n}.html" aria-label="查看 {n} 舊骨幹歷史快照">查看歷史快照</a></td></tr>').format(
            n=x["name"], nn=s["n"], p1=cp("c1"), p2=cp("c2"), p3=cp("c3p"),
            cbar=_cbar(s), nsh=s["nshapes"],
            inc=s["inc_pct"], ic=' style="color:#c92a2a;font-weight:700"' if s["inc_pct"] >= 10 else "",
            bcbar=_bcbar(s), sz=human(x["size"]))
    return """<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>[ARCHIVED] 舊克隆樹拓樸快照索引</title><style>
:root{{--archive:#641e1e;--paper:#fffdf8;--ink:#24201d;--muted:#6e6964;--line:#ddd4c8}}
*{{box-sizing:border-box}}html,body{{max-width:100%;overflow-x:hidden}}body{{margin:0;font-family:-apple-system,"Segoe UI","Noto Sans TC","Microsoft JhengHei",sans-serif;color:var(--ink);background:#f2eee8;line-height:1.55}}
.wrap{{width:100%;max-width:1120px;margin:0 auto;padding:22px}}.archive-hero{{overflow:hidden;background:var(--archive);color:#fff;border-top:8px solid #f7c948;border-radius:12px;padding:24px 26px;box-shadow:0 16px 40px rgba(57,20,20,.18)}}
.serial{{margin:0 0 10px;color:#ffe8a3;font-size:11px;font-weight:800;letter-spacing:.13em;text-transform:uppercase}}h1{{font-size:26px;line-height:1.25;margin:0 0 8px}}.hero-copy{{max-width:760px;margin:0;color:#ffecec;font-size:13.5px}}
.hero-actions{{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-top:18px}}.canonical{{display:inline-flex;flex-direction:column;padding:10px 14px;background:#fff;color:var(--archive);border:2px solid #fff;border-radius:8px;text-decoration:none}}
.canonical:hover,.canonical:focus-visible{{background:#fff8de;border-color:#f7c948;outline:3px solid rgba(247,201,72,.35);outline-offset:2px}}.canonical small{{font-size:9px;font-weight:800;letter-spacing:.1em;text-transform:uppercase}}.canonical strong{{font-size:14px}}.stamp{{font-size:11px;color:#ffdede}}
.rules{{display:grid;grid-template-columns:repeat(3,1fr);gap:10px;margin:14px 0}}.rule{{background:var(--paper);border:1px solid var(--line);border-top:4px solid #8f3c3c;border-radius:8px;padding:12px 14px;font-size:12px}}.rule b{{display:block;margin-bottom:3px;color:#5c1d1d;font-size:12.5px}}.rule.current{{border-top-color:#2f855a}}.rule.current b{{color:#256b49}}
.archive-list,.docs{{background:var(--paper);border:1px solid var(--line);border-radius:10px;margin:12px 0;overflow:hidden}}.archive-list>summary,.docs>summary{{display:flex;justify-content:space-between;gap:14px;align-items:center;padding:14px 16px;cursor:pointer;font-weight:750;color:#5c1d1d;list-style-position:inside}}
.archive-list>summary:hover,.archive-list>summary:focus-visible,.docs>summary:hover,.docs>summary:focus-visible{{background:#fbf2ec;outline:none}}.summary-note{{color:var(--muted);font-size:11px;font-weight:400}}.table-scroll{{max-width:100%;overflow-x:auto;border-top:1px solid var(--line)}}
table{{border-collapse:collapse;width:100%;min-width:960px;font-size:12px;background:#fff}}th,td{{padding:7px 9px;text-align:left;border-bottom:1px solid #eee8e0}}th{{background:#f3eee8;font-size:10.5px}}td{{vertical-align:middle}}.note{{color:var(--muted);font-size:10.5px}}
.snapshot-link{{display:inline-block;padding:4px 9px;border:1px solid #8f7f70;border-radius:6px;color:#5f554c;text-decoration:none;font-size:11px;white-space:nowrap}}.snapshot-link:hover,.snapshot-link:focus-visible{{background:#f3eee8;border-color:#5f554c;outline:2px solid rgba(95,85,76,.18);outline-offset:2px}}
.definitions{{padding:12px 16px;border-top:1px solid var(--line);font-size:11.5px;color:#5f5a55;line-height:1.8}}.docs-body{{padding:12px 16px;border-top:1px solid var(--line);font-size:11.5px;line-height:1.9;color:#5f5a55}}.mono{{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;overflow-wrap:anywhere}}.footer{{margin:16px 2px;color:var(--muted);font-size:10.5px}}
@media(max-width:720px){{.wrap{{padding:0 0 16px}}.archive-hero{{border-radius:0;padding:19px 16px;border-left:0;border-right:0}}h1{{font-size:22px}}.rules{{grid-template-columns:1fr;padding:0 12px;margin-top:12px}}.archive-list,.docs{{margin:12px;border-radius:8px}}.archive-list>summary,.docs>summary{{align-items:flex-start;flex-direction:column;gap:3px}}.footer{{margin:14px}}}}
</style></head><body><main class="wrap">
<header class="archive-hero"><p class="serial">Research archive · deprecated 2026-07-10</p><h1>舊克隆樹拓樸工作站 · 歷史快照索引</h1>
<p class="hero-copy">此資料庫使用已停用的 <span class="mono">is_somatic</span> backbone（normal VAF&lt;5% 粗重檢，已知誤殺 429 個 SEQC2-TP）。保留內容只供方法沿革與介面追溯；所有頁面均不可作為目前 validation evidence 或可引用結果。</p>
<div class="hero-actions"><a class="canonical" href="../layered_workstation/index.html"><small>Current canonical</small><strong>前往分層樹枚舉工作站 →</strong></a><span class="stamp">Archive ID · topology-v1 · readonly</span></div></header>
<section class="rules" aria-label="封存使用規則"><div class="rule"><b>為何封存</b>舊 backbone 會混淆 allelic / clonal 關係，且 normal VAF 粗重檢已有已知真陽性損失。</div><div class="rule"><b>可以怎麼用</b>僅供重現歷史畫面、理解舊方法與比對介面演進；不可引用頁內統計或延續人工判讀。</div><div class="rule current"><b>目前去哪裡</b>研究解讀、證據匯出與新工作一律使用 layered workstation。</div></section>
<details class="archive-list"><summary><span>查看 {nsamp} 個歷史單樣本快照</span><span class="summary-note">次要資料 · 預設收合 · 每頁皆有封存警示</span></summary>
<div class="table-scroll"><table><caption class="note" style="text-align:left;padding:10px 12px">舊 backbone 動態重生值；只描述此封存資料，不代表目前研究口徑。</caption><tr><th>樣本</th><th>舊總區</th><th>群數 c=1/2/≥3%</th><th>舊拓樸組成<br>單群/2平行/2直系</th><th>形狀<br>種類</th><th>舊 incompatible%</th><th>舊 branched / confirmed</th><th>檔案</th><th>歷史入口</th></tr>{tr}</table></div>
<div class="definitions"><b>舊欄位定義：</b>群數 c＝含 ≥1 ALT 向量的群數；拓樸組成＝此舊模型的 canonical 形狀；incompatible＝四配子／perfect-phylogeny 違反率；branched 與 confirmed 僅是封存模型內的描述。這些欄位不應跨版本解讀。</div></details>
<details class="docs"><summary><span>歷史文件與重現資訊</span><span class="summary-note">供追溯，不代表現行規格</span></summary><div class="docs-body">三軸分析：<span class="mono">InterSubMod/docs/methodology/20260702_topology_cluster_shape_three_axis_analysis_01.md</span><br>對抗驗證：<span class="mono">InterSubMod/docs/methodology/20260702_topology_three_axis_verification_log_01.md</span><br>建置地圖：<span class="mono">InterSubMod/docs/methodology/20260701_multisample_topology_workstation_build_map_01.md</span><br>重現命令：<span class="mono">python3 build_per_sample_workstation.py</span>（只為 archive reproduction）</div></details>
<p class="footer">ARCHIVED / NON-CITABLE · index 只列動態產生的歷史快照，不再重述容易漂移的跨樣本 narrative。</p>
</main></body></html>""".format(nsamp=len(rows), tr=tr)


def main():
    index_only = "--index-only" in sys.argv
    os.makedirs(OUT, exist_ok=True)
    pairs = sample_dirs()
    print("=== {} {} 樣本 → {} ===".format("index-only 重算" if index_only else "拆", len(pairs), OUT))
    rows = []
    for name, dr in pairs:
        summ = sample_summary(dr)
        if index_only:
            out = os.path.join(OUT, name + ".html")
            r = {"ok": os.path.exists(out), "path": out, "size": os.path.getsize(out) if os.path.exists(out) else 0}
        else:
            r = build_one(name, dr)
            print("  [{}] {} → {} ({})".format("OK" if r["ok"] else "FAIL", name, os.path.basename(r["path"]), human(r["size"])))
            if not r["ok"]:
                print("    LOG:", r["log"])
        rows.append({"name": name, "sum": summ, "size": r["size"], "ok": r["ok"]})
    cv = cramerv(rows)
    idx = os.path.join(OUT, "index.html")
    with open(idx, "w", encoding="utf-8") as f:
        f.write(index_html(rows, cv))
    print("  [OK] index.html ({}) CramérV={}".format(human(os.path.getsize(idx)), cv))
    total = sum(x["size"] for x in rows) + os.path.getsize(idx)
    print("=== 完成:{} 檔,合計 {};全部 ok={} ===".format(len(rows) + 1, human(total), all(x["ok"] for x in rows)))


if __name__ == "__main__":
    main()
