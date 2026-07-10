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
        tr += ('<tr><td><a href="{n}.html"><b>{n}</b></a><div class="note">{ti}</div></td>'
               '<td>{nn}</td><td class="note">{p1}/{p2}/{p3}%</td><td>{cbar}</td><td>{nsh}</td>'
               '<td{ic}>{inc}%</td><td>{bcbar}</td><td class="note">{sz}</td>'
               '<td><a class="btn" href="{n}.html">開啟 ▶</a></td></tr>').format(
            n=x["name"], ti=TI.get(x["name"], ""), nn=s["n"], p1=cp("c1"), p2=cp("c2"), p3=cp("c3p"),
            cbar=_cbar(s), nsh=s["nshapes"],
            inc=s["inc_pct"], ic=' style="color:#c92a2a;font-weight:700"' if s["inc_pct"] >= 10 else "",
            bcbar=_bcbar(s), sz=human(x["size"]))
    tot_n = sum(x["sum"]["n"] for x in rows)
    tot_inc = sum(x["sum"]["inc_n"] for x in rows)
    # 跨樣本觀察 — 兩離群 + 效應量
    by_br = sorted(rows, key=lambda x: -x["sum"]["br_pct"])
    top2 = "、".join("{}({}%)".format(x["name"], x["sum"]["br_pct"]) for x in by_br[:2])
    return """<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>多樣本克隆樹拓樸工作站 — 主頁（跨樣本多維度比較）</title><style>
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Segoe UI","Noto Sans TC","Microsoft JhengHei",sans-serif;color:#212529;background:#f8f9fa;line-height:1.5}}
.wrap{{max-width:1120px;margin:0 auto;padding:22px}}h1{{font-size:22px;margin:0 0 4px}}h2{{font-size:15px;color:#1971c2;margin:20px 0 8px}}
.sub{{color:#868e96;font-size:13px;margin:0 0 14px}}table{{border-collapse:collapse;width:100%;font-size:12.5px;background:#fff;border-radius:8px;overflow:hidden}}
th,td{{padding:6px 9px;text-align:left;border-bottom:1px solid #f1f3f5}}th{{background:#f1f3f5;font-size:11px}}td{{vertical-align:middle}}
.note{{color:#868e96;font-size:10.5px}}a{{color:#1971c2;text-decoration:none}}a:hover{{text-decoration:underline}}
.btn{{background:#1971c2;color:#fff;padding:4px 11px;border-radius:6px;font-size:11.5px;white-space:nowrap}}.btn:hover{{text-decoration:none;background:#1864ab}}
.card{{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:11px 15px;margin:10px 0;font-size:12px}}
.ribbon{{background:#fffbf5;border:1px solid #ffe0a3;border-radius:6px;padding:8px 12px;margin:10px 0;font-size:12px;color:#a37200}}
.mono{{font-family:ui-monospace,monospace}}
</style></head><body><div class="wrap">
<h1>🧬 多樣本克隆樹拓樸工作站 — 主頁（跨樣本多維度比較）</h1>
<p class="sub">{nsamp} ONT 樣本 · somatic sSNV 單分子共現重建 · 每樣本獨立檔(任何電腦可開) · 生成 by build_per_sample_workstation.py</p>
<div class="ribbon">🏷️ <b>characterization(非 subclone 判別) · single-bulk · L3/⭐3 · 跨樣本差異 partly-artifact</b> · CN 6/7 樣本補自 SAVANA(cna-only) · 統計正確 ≠ 生物差異證據</div>

<h2>📊 跨樣本多維度比較表</h2>
<table><tr><th>樣本</th><th>總區</th><th>群數 c=1/2/≥3%</th><th>拓撲組成<br>單群/2平行/2直系</th><th>形狀<br>種類</th><th>有效性<br>incompatible%</th><th>branched%(橙·幾何上界)<br>confirmed%(綠·read)</th><th>檔案</th><th></th></tr>{tr}
<tr style="background:#f8f9fa"><td><b>合計</b></td><td><b>{tot_n}</b></td><td class="note">7 樣本同版本</td><td class="note">11 種 canonical(全樣本合計)</td><td>—</td><td><b>{tot_inc_pct}%</b>({tot_inc})</td><td class="note">跨樣本 CramérV={cv}</td><td>—</td><td></td></tr></table>

<div class="card" style="background:#f4f6ff;border-color:#c3d0f5"><b>🔭 跨樣本觀察（不同維度整合）</b>
<div style="line-height:1.85;margin-top:4px">
• <b>維度①群數 c</b>：c=1(單群)為主(36–67%)、c=2 次之、c≥3 稀少(≤6%)→ 拓撲空間淺(全樣本僅 <b>11 種</b> canonical 樹形)<br>
• <b>維度②拓撲組成</b>：多數樣本 單群+2平行+2直系 三分;<b>COLO829 單群最高(67%)=結構最簡</b>、HCC1937/H2009 最複雜<br>
• <b>維度③有效性</b>：incompatible% <b>差 0.4↔19.1%</b>(COLO829 最乾淨、<b>H2009 19.1% 資料品質最差</b>)<br>
• <b>維度④read 驗證</b>：🔴 <b>branched%(幾何上界)常 > confirmed%(read 驗證)</b>—偏 branched 前二 <b>{top2}</b>,但 COLO829 confirmed 僅 39%(缺口最大=低 coread artifact)、HCC1954=undetermined<br>
• <b>跨樣本效應</b>：c=2 branched-vs-linear <b>CramérV={cv}</b>(small-med);<b>不報 chi2 p 值</b>(pseudoreplication·真 n=7);6/7 樣本 CN 未控 → 裁決 <b>partly-artifact / L3</b>(統計正確≠生物差異)
</div></div>

<div class="card"><b>📖 欄位意義</b><div class="note" style="line-height:1.85;margin-top:3px">
• <b>群數 c</b>=含≥1 ALT 向量的細胞群數(germline 不計);上界 c≤k<br>
• <b>拓撲組成</b>=canonical 樹形(單群 <span style="color:#868e96">■</span>/2平行 <span style="color:#f59f00">■</span>/2直系 <span style="color:#1c7ed6">■</span>);形狀種類=此樣本不同 canonical 樹形數(全 7 樣本合計 11 種、0 未分類)<br>
• <b>incompatible%</b>=四配子/perfect-phylogeny 違反率(🔴 誠實有效性訊號)<br>
• <b>branched%</b>=兩 ALT 向量非嵌套(幾何上界,不需 read) ｜ <b>confirmed%</b>=有 read 驗證結構的比例。缺口=幾何斷言但 read 未驗
</div></div>

<h2>📄 相關文件</h2><div class="card note" style="line-height:1.95">
• 三軸分析 + 定義附錄: <span class="mono">InterSubMod/docs/methodology/20260702_topology_cluster_shape_three_axis_analysis_01.md</span><br>
• 對抗驗證記錄(R1-R4 DRY): <span class="mono">InterSubMod/docs/methodology/20260702_topology_three_axis_verification_log_01.md</span><br>
• 建置地圖 / 本資料夾說明: <span class="mono">20260701_multisample_topology_workstation_build_map_01.md</span> / <span class="mono">README.md</span>(同資料夾)
</div>
<p class="note" style="margin-top:16px">重生:<span class="mono">python3 build_per_sample_workstation.py</span>(全部)或 <span class="mono">--index-only</span>(只主頁)。</p>
</div></body></html>""".format(nsamp=len(rows), tr=tr, tot_n=tot_n, tot_inc=tot_inc,
                                tot_inc_pct=round(100 * tot_inc / (tot_n or 1), 1), cv=cv if cv is not None else "—", top2=top2)


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
