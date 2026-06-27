#!/usr/bin/env python3
"""[全位點儀表板] 從 C++ 全基因組輸出建觀察儀表板 — 統計 + 分佈圖 + 全位點 compact 表 + 代表多群圖.
資料: phylo_cpp_wg_records.json(全位點) + summary.json + figs_cpp_wg/(多群圖). 自含 HTML(代表圖內嵌)."""
import json, os, base64, html as _h
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import sys; sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H  # noqa
A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
rec = json.load(open(f"{A}/phylo_cpp_wg_records.json"))
summ = json.load(open(f"{A}/phylo_cpp_wg_summary.json"))
tp, fp = summ["TP"], summ["FP"]
TP = [r for r in rec if r["set"] == "TP"]; FP = [r for r in rec if r["set"] == "FP"]

# ---- charts ----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
# (1) structure TP/FP
cats = ["aligned\ncis-ASM", "unaligned\nsubclone候選", "no_structure"]
x = np.arange(3); w = 0.38
ax[0].bar(x - w / 2, [tp["aligned_pct"], tp["unaligned_pct"], tp["no_structure_pct"]], w, label=f"TP({tp['n']})", color="#2563eb")
ax[0].bar(x + w / 2, [fp["aligned_pct"], fp["unaligned_pct"], fp["no_structure_pct"]], w, label=f"FP({fp['n']})", color="#dc2626")
ax[0].set_xticks(x); ax[0].set_xticklabels(cats, fontsize=9); ax[0].legend(fontsize=8); ax[0].set_ylabel("% loci")
ax[0].set_title("C++ 全基因組結構分類 TP vs FP", fontsize=10)
# (2) ngroups dist
gd = tp["ngroups_dist"]; ks = sorted(int(k) for k in gd)
ax[1].bar([str(k) for k in ks], [gd[str(k)] for k in ks], color="#16a34a")
ax[1].set_xlabel("coarse 群數"); ax[1].set_ylabel("位點數"); ax[1].set_title("TP 群數分佈", fontsize=10)
for i, k in enumerate(ks): ax[1].text(i, gd[str(k)], str(gd[str(k)]), ha="center", va="bottom", fontsize=7)
# (3) new layers
lyr = ["fine_multi", "other", "unstable", "hidden_het"]
ax[2].bar(np.arange(4) - w / 2, [tp["fine_multi_pct"], tp["other_pct"], tp["unstable_pct"], tp["hidden_het_pct"]], w, label="TP", color="#2563eb")
ax[2].bar(np.arange(4) + w / 2, [fp["fine_multi_pct"], fp["other_pct"], fp["unstable_pct"], fp["hidden_het_pct"]], w, label="FP", color="#dc2626")
ax[2].set_xticks(np.arange(4)); ax[2].set_xticklabels(lyr, fontsize=8.5, rotation=12); ax[2].legend(fontsize=8); ax[2].set_ylabel("% loci")
ax[2].set_title("v4.1 新層 (other/hidden_het FP偏多)", fontsize=10)
fig.suptitle(f"phylo-v4.1 C++ 原生全基因組 (HCC1395 tumor-only, TP {tp['n']}/FP {fp['n']}, {summ.get('n_figs',0)} 多群圖)", fontsize=11, y=1.02)
fig.tight_layout(); fig.savefig(f"{A}/figs_cpp_wg/_wg_charts.png", dpi=140, bbox_inches="tight"); plt.close(fig)
chart_b64 = "data:image/png;base64," + base64.b64encode(open(f"{A}/figs_cpp_wg/_wg_charts.png", "rb").read()).decode()

# ---- 代表多群圖 (內嵌 ~12: 各群數+aligned/unaligned 代表) ----
def has_fig(r):
    p = f"{A}/figs_cpp_wg/cpp_chr{r['chrom'][3:]}_{int(r['pos'])}_{int(r['pos'])+10000}.png"
    # region fig naming = cpp_{chrN_start_end}.png ; start=pos
    cand = [x for x in os.listdir(f"{A}/figs_cpp_wg") if x.startswith(f"cpp_{r['chrom']}_{r['pos']}_")]
    return cand[0] if cand else None
samples = []
multi = [r for r in TP if r["coarse_ng"] >= 2]
seen_ng = {}
for r in sorted(multi, key=lambda z: (-z["coarse_ng"], -z["n"])):
    f = has_fig(r)
    if not f: continue
    key = (r["coarse_ng"], r["aligned"])
    if seen_ng.get(key, 0) >= 2: continue
    seen_ng[key] = seen_ng.get(key, 0) + 1
    samples.append((r, f))
    if len(samples) >= 12: break
samp_html = ""
for r, f in samples:
    b = "data:image/png;base64," + base64.b64encode(open(f"{A}/figs_cpp_wg/{f}", "rb").read()).decode()
    tag = "aligned cis-ASM" if r["aligned"] else "🔴unaligned"
    samp_html += f"<div style='margin:10px 0'><div class='sub'><b>{r['chrom']}:{int(r['pos'])+5000}</b> n={r['n']} coarse {r['coarse_ng']}群 {tag} (V_allele {r['V_allele']}){' ⚠unstable' if r['unstable'] else ''}{' other'+str(r['n_other']) if r['n_other'] else ''}</div><img src='{b}' style='width:100%;border:1px solid #2a3344;border-radius:6px'></div>"

# ---- 全位點 compact 表 (JS 篩選) ----
rows_js = json.dumps([{"c": r["chrom"], "p": int(r["pos"]) + 5000, "s": r["set"], "n": r["n"], "g": r["coarse_ng"],
                       "f": r["fine_ng"], "o": r["n_other"], "u": 1 if r["unstable"] else 0,
                       "a": 1 if r["aligned"] else 0, "h": 1 if r["hidden_het"] else 0} for r in rec])

HTML = f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>phylo-v4.1 C++ 全位點觀察儀表板</title><style>
:root{{--bg:#0a0e16;--card:#111826;--fg:#e7edf5;--mut:#8b9bb4;--line:#2a3344;--ac:#D97757}}
body{{font-family:system-ui,"PingFang TC","Noto Sans CJK TC",sans-serif;background:var(--bg);color:var(--fg);margin:0;line-height:1.55}}
.wrap{{max-width:1280px;margin:0 auto;padding:22px 18px 80px}} h1{{font-size:20px}} h2{{font-size:15px;border-left:4px solid var(--ac);padding-left:9px;margin-top:24px}}
.sub{{color:var(--mut);font-size:12.5px}} .box{{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:13px 16px;margin:10px 0}}
table{{border-collapse:collapse;width:100%;font-size:12px}} th,td{{border:1px solid var(--line);padding:3px 7px;text-align:center}} th{{background:#0b1222;cursor:pointer;position:sticky;top:0}}
td:first-child{{text-align:left}} .tw{{max-height:560px;overflow:auto;border:1px solid var(--line);border-radius:8px}}
input,select{{background:#0b1222;color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:5px 8px;font-size:12.5px}}
.k{{font-weight:700;color:var(--ac)}} .ok{{color:#4ade80}} .no{{color:#f87171}}</style></head><body><div class="wrap">
<h1>phylo-v4.1 C++ 原生全基因組觀察儀表板</h1>
<div class="sub">HCC1395 tumor-only · binary 原生 phylo_groups.tsv(K=10 modal/fine/other) · Python 只讀+畫 · TP {tp['n']}/FP {fp['n']} · {summ.get('n_figs',0)} 多群圖 · build {summ.get('source','')}</div>

<div class="box"><b>🔴 一句話</b>：C++ 原生確認 v3.1/v4.1 結論 —— 有結構 TP <span class="k">{tp['structure_pct']}%</span>，其中 aligned cis-ASM <span class="k">{tp['aligned_pct']}%</span>(FP {fp['aligned_pct']}%, 偏TP) / unaligned subclone候選 <span class="k">{tp['unaligned_pct']}%</span>(FP <span class="no">{fp['unaligned_pct']}%</span>, 偏FP=非subclone)。</div>

<h2>① 全基因組統計圖</h2><img src="{chart_b64}" style="width:100%;border-radius:8px">

<h2>② 全位點表（{len(rec)} 位點，可篩選/排序）</h2>
<div class="box"><span class="sub">篩選:</span>
 set <select id="fset"><option value="">全</option><option>TP</option><option>FP</option></select>
 群數≥ <input id="fg" type="number" value="0" min="0" max="9" style="width:55px">
 <label><input type="checkbox" id="fa"> 僅aligned</label>
 <label><input type="checkbox" id="fu"> 僅unstable</label>
 <label><input type="checkbox" id="fh"> 僅hidden_het</label>
 <span class="sub" id="cnt"></span></div>
<div class="tw"><table id="tb"><thead><tr>
<th data-k="c">chrom</th><th data-k="p">pos</th><th data-k="s">set</th><th data-k="n">n</th><th data-k="g">coarse群</th><th data-k="f">fine群</th><th data-k="o">other</th><th data-k="u">unstable</th><th data-k="a">aligned</th><th data-k="h">hidden_het</th>
</tr></thead><tbody id="tbody"></tbody></table></div>

<h2>③ 代表多群圖（內嵌；完整 {summ.get('n_figs',0)} 圖在 figs_cpp_wg/）</h2>
<div class="box sub">點放大看 C++ coarse 標籤(side bar)是否對齊 HP/ALT。完整多群圖在本機 <code>docs/methodology/_assets/20260618_subcluster_pilot/figs_cpp_wg/</code>。</div>
{samp_html}

<div class="box sub">資料源 phylo_cpp_wg_records.json + summary.json(C++ native)。單樣本 ⭐2-3 characterization。</div>
</div>
<script>
const D={rows_js};
const tb=document.getElementById('tbody'); let sortK='g',sortDir=-1;
function draw(){{
 const fs=document.getElementById('fset').value, fg=+document.getElementById('fg').value;
 const fa=document.getElementById('fa').checked, fu=document.getElementById('fu').checked, fh=document.getElementById('fh').checked;
 let r=D.filter(x=>(!fs||x.s===fs)&&x.g>=fg&&(!fa||x.a)&&(!fu||x.u)&&(!fh||x.h));
 r.sort((a,b)=>(a[sortK]>b[sortK]?1:a[sortK]<b[sortK]?-1:0)*sortDir);
 document.getElementById('cnt').textContent=' 顯示 '+r.length+' / '+D.length;
 tb.innerHTML=r.slice(0,3000).map(x=>`<tr><td>${{x.c}}</td><td>${{x.p}}</td><td>${{x.s}}</td><td>${{x.n}}</td><td>${{x.g>=2?'<b>'+x.g+'</b>':x.g}}</td><td>${{x.f}}</td><td>${{x.o||''}}</td><td>${{x.u?'<span class=no>✓</span>':''}}</td><td>${{x.a?'<span class=ok>✓</span>':''}}</td><td>${{x.h?'<span class=no>✓</span>':''}}</td></tr>`).join('')
  +(r.length>3000?`<tr><td colspan=10 class=sub>(顯示前3000, 共${{r.length}}; 用篩選縮小)</td></tr>`:'');
}}
document.querySelectorAll('th[data-k]').forEach(t=>t.onclick=()=>{{const k=t.dataset.k; sortDir=(sortK===k)?-sortDir:-1; sortK=k; draw();}});
['fset','fg','fa','fu','fh'].forEach(id=>document.getElementById(id).oninput=draw);
draw();
</script></body></html>"""
out = f"{A}/20260623_phylo_cpp_all_loci_dashboard.standalone.html"
open(out, "w").write(HTML)
print(f"WROTE {out} ({round(len(HTML)/1024)} KB, {len(rec)} loci, {len(samples)} 代表圖)")
