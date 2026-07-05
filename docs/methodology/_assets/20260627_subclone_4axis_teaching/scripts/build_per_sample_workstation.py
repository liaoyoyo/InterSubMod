#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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


def sample_summary(dr):
    """從 topology_per_region.json 算主頁要顯示的關鍵數字(fresh,不抄中間檔)。"""
    d = json.load(open(os.path.join(dr, "topology_per_region.json"), encoding="utf-8"))["detail"]
    n = len(d)
    inc = sum(1 for r in d if r.get("determinacy") == "incompatible")
    c = Counter(r["n_clusters"] for r in d)
    c2 = [r for r in d if r["n_clusters"] == 2]
    n2 = len(c2) or 1
    br = sum(1 for r in c2 if "branched" in (r.get("topology_type") or ""))
    cf = sum(1 for r in c2 if r.get("tree_shape") in CONF)
    return {"n": n, "inc_pct": round(100 * inc / n, 1), "c1": c.get(1, 0), "c2": len(c2),
            "c3p": sum(v for k, v in c.items() if k >= 3),
            "br_pct": round(100 * br / n2, 1), "conf_pct": round(100 * cf / n2, 1)}


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


def index_html(rows):
    tr = ""
    for x in rows:
        s = x["sum"]
        tr += ('<tr><td><a href="{n}.html"><b>{n}</b></a></td><td class="note">{ti}</td>'
               '<td>{nn}</td><td>{c1}/{c2}/{c3p}</td><td{ic}>{inc}%</td><td>{br}%</td><td>{cf}%</td>'
               '<td>{sz}</td><td><a class="btn" href="{n}.html">開啟 ▶</a></td></tr>').format(
            n=x["name"], ti=TI.get(x["name"], ""), nn=s["n"], c1=s["c1"], c2=s["c2"], c3p=s["c3p"],
            inc=s["inc_pct"], ic=' style="color:#c92a2a;font-weight:700"' if s["inc_pct"] >= 10 else "",
            br=s["br_pct"], cf=s["conf_pct"], sz=human(x["size"]))
    tot_n = sum(x["sum"]["n"] for x in rows)
    tot_inc = sum(round(x["sum"]["inc_pct"] * x["sum"]["n"] / 100) for x in rows)
    return """<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>多樣本克隆樹拓樸工作站 — 主頁（每樣本檢視）</title><style>
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Segoe UI","Noto Sans TC","Microsoft JhengHei",sans-serif;color:#212529;background:#f8f9fa;line-height:1.55}}
.wrap{{max-width:1080px;margin:0 auto;padding:22px}}h1{{font-size:22px;margin:0 0 4px}}h2{{font-size:15px;color:#1971c2;margin:22px 0 8px}}
.sub{{color:#868e96;font-size:13px;margin:0 0 16px}}table{{border-collapse:collapse;width:100%;font-size:13px;background:#fff;border-radius:8px;overflow:hidden}}
th,td{{padding:7px 10px;text-align:left;border-bottom:1px solid #f1f3f5}}th{{background:#f1f3f5;font-size:11.5px}}td{{vertical-align:middle}}
.note{{color:#868e96;font-size:11px}}a{{color:#1971c2;text-decoration:none}}a:hover{{text-decoration:underline}}
.btn{{background:#1971c2;color:#fff;padding:4px 12px;border-radius:6px;font-size:12px;white-space:nowrap}}.btn:hover{{text-decoration:none;background:#1864ab}}
.card{{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:12px 16px;margin:10px 0}}
.ribbon{{background:#fffbf5;border:1px solid #ffe0a3;border-radius:6px;padding:8px 12px;margin:10px 0;font-size:12px;color:#a37200}}
</style></head><body><div class="wrap">
<h1>🧬 多樣本克隆樹拓樸工作站 — 主頁</h1>
<p class="sub">{nsamp} ONT 樣本 · somatic sSNV 單分子共現重建 · 每樣本獨立檔案(~9MB,任何電腦可開) · 生成 by build_per_sample_workstation.py</p>
<div class="ribbon">🏷️ characterization(非 subclone 判別) · single-bulk · L3/⭐3 · 跨樣本差異 partly-artifact(見驗證記錄) · CN 6/7 樣本補自 SAVANA(cna-only)</div>
<h2>逐樣本檢視（點「開啟」看該樣本完整工作站）</h2>
<table><tr><th>樣本</th><th>備註</th><th>總區</th><th>c=1/2/≥3</th><th>incompatible%</th><th>c=2 branched%</th><th>c=2 confirmed%</th><th>檔案</th><th></th></tr>{tr}
<tr style="background:#f8f9fa"><td><b>合計</b></td><td class="note">7 樣本同版本</td><td><b>{tot_n}</b></td><td>—</td><td><b>{tot_inc_pct}%</b>({tot_inc})</td><td>—</td><td>—</td><td>—</td><td></td></tr></table>
<div class="card"><b>📖 欄位意義</b>
<div class="note" style="line-height:1.9;margin-top:4px">
• <b>總區</b>=n_sSNV≥2 的最大單分子連鎖區數 ｜ <b>c=1/2/≥3</b>=群數分布(c=含≥1 ALT 向量數,germline 不計)<br>
• <b>incompatible%</b>=四配子/perfect-phylogeny 違反率(🔴 誠實有效性訊號;逐樣本 0.4%–19.1%,H2009 最差)<br>
• <b>c=2 branched%</b>=兩 ALT 向量非嵌套(=幾何上界,<b>非</b> read 驗證) ｜ <b>c=2 confirmed%</b>=有 read 驗證結構的比例(對照上界)<br>
• 🔴 branched%>confirmed% 的缺口=幾何斷言但 read 未驗(COLO829 87%→僅 39% 驗證)
</div></div>
<h2>📄 相關文件（分析 + 驗證）</h2>
<div class="card note" style="line-height:2">
• 三軸分析報告(c×拓撲×樣本 + 定義附錄): <span class="mono">InterSubMod/docs/methodology/20260702_topology_cluster_shape_three_axis_analysis_01.md</span><br>
• 迭代對抗驗證記錄(R1-R4 DRY 收斂): <span class="mono">InterSubMod/docs/methodology/20260702_topology_three_axis_verification_log_01.md</span><br>
• 建置地圖(產生程式/資料流): <span class="mono">InterSubMod/docs/methodology/20260701_multisample_topology_workstation_build_map_01.md</span><br>
• 本資料夾說明 + 重生指令: <span class="mono">README.md</span>(同資料夾)
</div>
<p class="note" style="margin-top:18px">重生:<span class="mono">cd .../20260627_subclone_4axis_teaching/scripts &amp;&amp; python3 build_per_sample_workstation.py</span> → 覆寫本資料夾全部檔案。</p>
</div></body></html>""".format(nsamp=len(rows), tr=tr, tot_n=tot_n, tot_inc=tot_inc,
                                tot_inc_pct=round(100 * tot_inc / (tot_n or 1), 1))


def main():
    os.makedirs(OUT, exist_ok=True)
    pairs = sample_dirs()
    print("=== 拆 {} 樣本 → {} ===".format(len(pairs), OUT))
    rows = []
    for name, dr in pairs:
        summ = sample_summary(dr)
        r = build_one(name, dr)
        status = "OK" if r["ok"] else "FAIL"
        print("  [{}] {} → {} ({})".format(status, name, os.path.basename(r["path"]), human(r["size"])))
        if not r["ok"]:
            print("    LOG:", r["log"])
        rows.append({"name": name, "sum": summ, "size": r["size"], "ok": r["ok"]})
    idx = os.path.join(OUT, "index.html")
    with open(idx, "w", encoding="utf-8") as f:
        f.write(index_html(rows))
    print("  [OK] index.html ({})".format(human(os.path.getsize(idx))))
    total = sum(x["size"] for x in rows) + os.path.getsize(idx)
    print("=== 完成:{} 檔,合計 {} (原單檔 63.6MB);全部 ok={} ===".format(
        len(rows) + 1, human(total), all(x["ok"] for x in rows)))


if __name__ == "__main__":
    main()
