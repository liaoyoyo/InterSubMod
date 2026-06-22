#!/usr/bin/env python3
"""切群重設計 觀察 HTML（§13-A：數字全從 cluster_redesign.json 注入；圖 base64）。
三閘(balance群≥3 + null-excess明顯 + alignment歸因) + coarse/fine 雙輸出 + A/B/C 三類 + per-read。"""
import json, base64, os, sys, html
from collections import Counter
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
d=json.load(open(f"{A}/cluster_redesign.json")); loci=d["loci"]; P=d["params"]
try: figidx=json.load(open(f"{A}/cluster_redesign_figindex.json"))["figs"]
except Exception: figidx=[]
try: dendro=json.load(open(f"{A}/dendro_figindex.json"))["gaps"]
except Exception: dendro={}
def b64(rel):
    p=os.path.join(A,rel)
    return "data:image/png;base64,"+base64.b64encode(open(p,"rb").read()).decode() if os.path.exists(p) else None
def cls_of(r):
    cf=r["coarse_confidence"]; ff=r["fine_confidence"]; nov=r["n_novel"]
    if ff=="NO_CLEAR_SPLIT": return "D"
    if cf=="CONFIRMED" and ff=="CONFIRMED" and nov==0: return "A"
    if cf=="CONFIRMED": return "B"                       # 骨幹 + novel/diffuse 候選
    if cf=="NONE" and ff in("REAL_NOVEL","REAL_DIFFUSE"): return "C"
    return "—"
CC={"A":("#0891b2","germline-cis(非subclone)"),"B":("#7c3aed","骨幹+subclone候選"),
    "C":("#db2777","純novel/diffuse(候選/需cis-control)"),"D":("#6b7280","真無法切(全無真實)"),"—":("#888","其他")}
clscnt=Counter(cls_of(r) for r in loci)
n=d["n_loci"]; novn=d.get("has_REAL_NOVEL(subclone候選)"); edn=d.get("has_edge_group"); mrn=d.get("FM3_multiresolution_confirmed(≥2)")

# 代表圖卡
cards=""
for fk in ["chr8_64387206","chr8_138619384","chr4_190112507","chr8_132099309","chr17_79567788","chr1_218115208","chr10_55938697","chr2_44798355"]:
    img=b64(f"figs_redesign/redesign_{fk}.png"); r=next((x for x in loci if x["key"]==fk),None)
    if not img or not r: continue
    cl=cls_of(r); cc=CC[cl][0]
    cards+=f"""<div class="card"><div class="ch"><span class="cb" style="background:{cc}">類{cl}</span>
      <b>{html.escape(fk)}</b> <span class="mut">coarse k={r['coarse_k']}({r['coarse_confidence']}) | fine k={r['fine_k']}({r['fine_confidence']}) · confirmed={r['confirmed_ks']} novel={r['novel_ks']}</span></div>
      <img src="{img}"></div>"""

# 全 29 表
rows=""
for r in sorted(loci,key=lambda r:(cls_of(r),r["key"])):
    cl=cls_of(r); cc=CC[cl][0]
    nfine=Counter(p["type"] for p in (r["perread_fine"] or []))
    rows+=(f"<tr><td><span class='cb' style='background:{cc};font-size:10px'>{cl}</span></td>"
           f"<td style='text-align:left'>{html.escape(r['group'])}</td><td>{html.escape(r['key'])}</td><td>{r['n']}</td>"
           f"<td>{r['old']['k_groups']}</td><td>{r['coarse_k']}<br><span class='mut'>{r['coarse_confidence']}</span></td>"
           f"<td>{r['fine_k']}<br><span class='mut'>{r['fine_confidence']}</span></td>"
           f"<td>{r['confirmed_ks']}</td><td>{r['novel_ks']}</td>"
           f"<td>{nfine.get('core_confirmed',0)+nfine.get('core_novel',0)}核/{nfine.get('edge',0)}邊/{nfine.get('outlier',0)}離</td></tr>")

# dendrogram 區（gap 階梯表 + 12 樹圖）
GAP_TH=8.0
grows=""
for key,g in sorted(dendro.items(),key=lambda kv:kv[1].get("group","")):
    gr=g.get("gap_ratio",{})
    cells="".join(f"<td style='{'background:#d1fae5;font-weight:700' if (gr.get(str(k)) or gr.get(k) or 0)>=GAP_TH else ''}'>{gr.get(str(k),gr.get(k,'-'))}</td>" for k in range(2,7))
    grows+=(f"<tr><td style='text-align:left'>{html.escape(g.get('group',''))}</td><td>{html.escape(key)}</td><td>{g['n']}</td>"
            f"<td>{g.get('coarse_k')}</td><td>{g.get('fine_k')}</td>{cells}</tr>")
dimgs=""
for key,g in dendro.items():
    img=b64(g["png"])
    if img: dimgs+=f"""<div class="card"><div class="ch"><b>{html.escape(key)}</b> <span class="mut">{html.escape(g.get('group',''))} n={g['n']} · coarse k={g.get('coarse_k')} / fine k={g.get('fine_k')} · gap 大跳 k={[k for k in range(2,7) if (g.get('gap_ratio',{}).get(str(k)) or 0)>=GAP_TH]}</span></div><img src="{img}"></div>"""
dendro_section=f"""<h2>4. UPGMA 距離樹圖 + gap 階梯（精細度天花板候選標準）</h2>
<div class="box key"><b>問題</b>:excess(扣 null)在高覆蓋位點隨 k 自動上升、定不出天花板(修 stab_excess bug 後 REAL_NOVEL 衝 27/29)。
<b>解</b>:樹的<b>「分支跳躍 gap」</b>才是客觀訊號 — 大跳=真分界、≈1×=非真分界。下表 gap-ratio(×中位 gap),<b>≥{GAP_TH}×</b> 標綠=結構支持的切點;搭配 excess(真實)+ alignment(歸因)。
例 chr4:k2={dendro.get('chr4_190112507',{}).get('gap_ratio',{}).get('2','?')}× + k5={dendro.get('chr4_190112507',{}).get('gap_ratio',{}).get('5','?')}× 兩大跳 → ALT 子群在 k=5,該切;k3/4/6 中等。chr2_44798355(NO_CLEAR)全 ≤5× → 無清楚切。⚠ loose 但對齊者 gap 小仍可能有意義 → 標準 = <b>gap 大 OR 可靠對齊</b>。</div>
<table><tr><th>group</th><th>locus</th><th>n</th><th>coarse k</th><th>fine k</th><th>gap k2</th><th>k3</th><th>k4</th><th>k5</th><th>k6</th></tr>{grows}</table>
<div class="sub">每張:左=UPGMA 樹(分支色=fine k cut,橙虛線=k2/3/4 cut 高度)+ allele/hp/carrier 葉色帶 + 右距離熱圖</div>
{dimgs}
"""

HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>切群重設計觀察 — 三閘 + coarse/fine</title><style>
:root{{--ac:#D97757;--tx:#141413;--bg:#FAF9F5;--bd:#E3DACC}}
*{{box-sizing:border-box}} body{{font-family:system-ui,-apple-system,"PingFang TC","Microsoft JhengHei","Noto Sans CJK TC","Droid Sans Fallback",sans-serif;color:var(--tx);background:var(--bg);margin:0;line-height:1.6}}
.wrap{{max-width:1200px;margin:0 auto;padding:26px 20px}}
h1{{font-size:22px;margin:0 0 4px}} h2{{font-size:17px;border-left:4px solid var(--ac);padding-left:10px;margin-top:28px}}
.sub{{color:#6b6b6b;font-size:13px;margin-bottom:12px}}
.chip{{display:inline-block;color:#fff;border-radius:12px;padding:3px 10px;margin:3px 4px 0 0;font-size:12px;font-weight:600}}
table{{border-collapse:collapse;width:100%;font-size:12px;margin:10px 0}} th,td{{border:1px solid var(--bd);padding:4px 6px;text-align:center}} th{{background:#f0ebe2}}
.card{{border:1px solid var(--bd);border-radius:10px;margin:14px 0;background:#fff;overflow:hidden}}
.ch{{padding:8px 12px;background:#faf6ef;border-bottom:1px solid var(--bd);font-size:13px}} .card img{{max-width:100%;display:block}}
.cb{{display:inline-block;color:#fff;border-radius:5px;padding:2px 7px;font-size:11px;font-weight:600}}
.mut{{color:#888;font-size:11.5px;font-weight:400}}
.box{{border:1px solid var(--bd);border-radius:8px;padding:12px 16px;background:#fff;margin:12px 0}}
.box.key{{border-color:var(--ac);background:#fdf6f1}} code{{background:#f0ebe2;padding:1px 5px;border-radius:3px;font-size:12px}}
.foot{{margin-top:30px;border-top:1px solid var(--bd);padding-top:12px;color:#888;font-size:12px}}
</style></head><body><div class="wrap">
<h1>切群重設計觀察 — 三閘 + coarse/fine 雙輸出</h1>
<div class="sub">HCC1395 tumor-only · n={n} 代表位點 · 單樣本 ⭐2-3 characterization · 純分析(未改 C++)</div>

<div class="box key"><b>四閘「切更細」標準</b>:
<b>① balance</b> 群≥3(業界 caller 最低) · <b>② null-excess</b> 扣 within-1-group null ≥{P['exc_min']}(「明顯」=真實) · <b>③ gap-significance</b> 分支跳躍 ≥max(8×中位, 0.4×最大gap)(scale-invariant,定精細度天花板) · <b>④ alignment</b> 可靠對齊 germline。
每群標 <b>CONFIRMED</b>(真實+對齊)/ <b>REAL_NOVEL</b>(真實+大跳不對齊=<b>subclone候選</b>)/ <b>REAL_DIFFUSE</b>(真實但無大跳+無對齊,低信心候選)/ <b>NO_CLEAR</b>(全無真實=真無法切)/ edge(≥3離群)/ outlier。
<b>雙輸出</b>:coarse=最粗 CONFIRMED 骨幹 · fine=最細 supported。</div>
{f'<div class="card"><div class="ch"><b>方法流程圖(workflow 確認)</b></div><img src="{b64("figs_dendro/method_flowchart.png")}"></div>' if b64("figs_dendro/method_flowchart.png") else ""}

<h2>0. 三類位點(對 subclone 目標)</h2>
<div>{"".join(f'<span class="chip" style="background:{CC[k][0]}">類{k} {CC[k][1]}: {clscnt.get(k,0)}</span>' for k in ["A","B","C","—"])}</div>
<div class="box"><b>🔑 重點</b>:<b>B+{clscnt.get('B',0)} + C {clscnt.get('C',0)} = {clscnt.get('B',0)+clscnt.get('C',0)} 個有 REAL_NOVEL</b>(真實但不跟 germline 走)= subclone 藏身處(真 subclone 甲基不該被 germline haplotype 決定)。
類 A {clscnt.get('A',0)} 個各 k 都對齊 germline = cis-ASM,非 subclone。has_REAL_NOVEL={novn} · 邊緣群={edn} · 多解析度 confirmed={mrn}。</div>

<h2>1. 精細度標準(完整+快速確認)</h2>
<div class="box">每位點掃 k=2..{P['kmax']} 算 excess(扣 null)曲線 → <b>精細度天花板 = excess 掉到 {P['exc_min']} 以下、或開始出現 <3 單離群那一階</b>。完整(掃所有 k)、快速(僅 null 成本秒級)、客觀(門檻明寫)。「小群沒切」= 看卡哪閘:excess<{P['exc_min']}→不夠明顯;群<3→歸邊緣;不對齊→REAL_NOVEL 記錄。</div>
{f'<img src="{b64("figs_redesign/granularity_curve.png")}" style="max-width:100%;border:1px solid var(--bd);border-radius:6px">' if b64("figs_redesign/granularity_curve.png") else ""}

<h2>2. 代表位點 old vs new(coarse/fine + per-read)</h2>
<div class="sub">側欄 OLD maxclust | NEW primary(=fine) | hp 藍/橙 · carrier 綠/紅 · allele 青/桃;標題含 excess(扣 null)+ 邊緣群</div>
{cards}

<h2>3. 全 {n} 位點表</h2>
<table><tr><th>類</th><th>group</th><th>locus</th><th>n</th><th>OLD k</th><th>coarse</th><th>fine</th><th>confirmed_ks</th><th>novel_ks</th><th>per-read(核/邊/離)</th></tr>{rows}</table>

{dendro_section}
<h2>5. 限制</h2>
<ul>
<li>代表位點皆已 splittable(無純無訊號);<b>FM1 的 CANT_SPLIT→可切 rescue 要全基因組才看得到</b>。</li>
<li>單樣本 HCC1395 ⭐2-3:REAL_NOVEL=subclone <b>候選</b>,需 normal cis-control 排 cis-ASM 才能確認為 subclone。</li>
<li>門檻(群≥3 / excess≥{P['exc_min']} / e≥5 / gap×{P['gap_mult']})為約定;per-k 全數字在 json 可 re-threshold。</li>
</ul>
<div class="foot">資料源 <code>cluster_redesign.json</code>(seed={P['seed']})。重生:<code>cluster_redesign.py → plot_cluster_redesign.py → granularity_diagnostic.py → build_redesign_observation.py</code>。三閘 = balance+null-excess+alignment。單樣本 characterization 非 subclone 確認。</div>
</div></body></html>"""
out=f"{A}/20260622_cluster_redesign_observation_01.standalone.html"
open(out,"w").write(HTML)
print(f"WROTE {out} ({len(HTML)//1024} KB)  A={clscnt.get('A',0)} B={clscnt.get('B',0)} C={clscnt.get('C',0)}")
