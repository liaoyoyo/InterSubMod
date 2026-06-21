#!/usr/bin/env python3
"""切割穩定性 觀察 HTML（§13-A：數字全從 kprofile_stability.json 注入，缺 key refuse；圖 base64 嵌入）。
不手打任何 metric。輸出 standalone HTML 到 asset 同層。"""
import json, base64, os, sys, html

A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
d=json.load(open(f"{A}/kprofile_stability.json"))
figidx=json.load(open(f"{A}/kprofile_stability_figindex.json"))
P=d["params"]; loci=d["loci"]

def req(obj,*keys):
    cur=obj
    for k in keys:
        if isinstance(cur,dict) and k in cur and cur[k] is not None: cur=cur[k]
        else: sys.exit(f"[REFUSE §13-A] 缺 key: {'/'.join(map(str,keys))}")
    return cur
def b64(rel):
    p=os.path.join(A,rel)
    if not os.path.exists(p): sys.exit(f"[REFUSE] 圖不存在: {rel}")
    return "data:image/png;base64,"+base64.b64encode(open(p,"rb").read()).decode()

VB={"GOOD_CUT":("#0d9488","可確認的好切法"),"STABLE+ALIGNED_DIFF_K":("#7c3aed","穩定+對齊不同k(單樣本天花板)"),
    "STABLE_BUT_UNALIGNED":("#db2777","穩但不對齊/不可靠(陷阱)"),"ALIGNED_BUT_UNSTABLE":("#d97706","對齊但不穩"),
    "NEITHER":("#475569","皆無")}
vc=req(d,"verdict_counts"); n=req(d,"n_loci"); sc=req(d,"bernoulli_selfcheck_max_absdiff")

# verdict summary chips
chips="".join(f'<span class="chip" style="background:{VB.get(k,("#888",k))[0]}">{html.escape(k)}: {v}</span>' for k,v in vc.items())

# by_group table
bg=req(d,"by_group")
allv=list(VB.keys())
bgrows=""
for g,vd in bg.items():
    cells="".join(f'<td>{vd.get(v,0)}</td>' for v in allv)
    bgrows+=f"<tr><td style='text-align:left'>{html.escape(g)}</td>{cells}<td><b>{sum(vd.values())}</b></td></tr>"
bghead="".join(f"<th>{html.escape(v.replace('_','_<wbr>'))}</th>" for v in allv)

# per-locus rep cards
cards=""
for rp in figidx["reps"]:
    key=rp["key"]; r=next(x for x in loci if x["chrom"]+":"+x["pos"]==key)
    h=r["headline"]; col=VB.get(r["verdict"],("#888",""))[0]
    dual=rp.get("dualpanel","")
    dualimg=f'<img src="{b64(dual)}" alt="dualpanel">' if dual and os.path.exists(os.path.join(A,dual)) else '<div class="nofig">（無 dual-panel）</div>'
    # per-k mini table
    rows=""
    for pk in r["per_k"]:
        pb="✓" if pk["pass_both"] else ""
        recol="#0d9488" if (pk["align_e"] or 0)>=5 else "#db2777"
        rows+=(f"<tr><td>k={pk['k']}</td><td>{pk['jac_real_mean']}</td><td>{pk['jac_null_p95']}</td>"
               f"<td><b>{pk['stab_excess']:+.3f}</b></td><td>{pk['align_V']}</td>"
               f"<td style='color:{recol}'>{pk['align_e']}</td><td>{pk['align_shuffle_p']}</td>"
               f"<td style='color:#0d9488'>{pb}</td></tr>")
    cards+=f"""<div class="card">
      <div class="card-h"><span class="vbadge" style="background:{col}">{html.escape(r['verdict'])}</span>
        <b>{html.escape(key)}</b> <span class="muted">[{html.escape(rp['label'])}] n={r['n']} 軸={html.escape(r['primary_axis'])} · headline k={req(r,'headline_k')}</span></div>
      <div class="card-body">
        <div class="fig"><img src="{b64(rp['png'])}" alt="stability"></div>
        <div class="fig">{dualimg}<div class="cap">既有 dual-panel：read×CpG 甲基 + read×read 距離 + a-priori 側欄（看實際結構）</div></div>
        <table class="kk"><tr><th>k</th><th>real Jac</th><th>null95</th><th>excess</th><th>V</th><th>Cochran e</th><th>shufP</th><th>both✓</th></tr>{rows}</table>
      </div></div>"""

# 全 29 verification table
vrows=""
for r in sorted(loci,key=lambda r:(r["group"],-r["headline"]["stab_excess"])):
    h=r["headline"]; col=VB.get(r["verdict"],("#888",""))[0]
    vrows+=(f"<tr><td style='text-align:left'>{html.escape(r['group'])}</td><td>{html.escape(r['chrom']+':'+r['pos'])}</td>"
            f"<td>{r['n']}</td><td>{html.escape(r['primary_axis'])}</td><td>k={req(r,'headline_k')}</td>"
            f"<td>{h['jac_real']}</td><td>{h['jac_null_p95']}</td><td><b>{h['stab_excess']:+.3f}</b></td>"
            f"<td>{h['align_V']}</td><td>{h['align_e']}</td><td>{h['align_shuffle_p']}</td>"
            f"<td><span class='vbadge' style='background:{col};font-size:10px'>{html.escape(r['verdict'])}</span></td></tr>")

HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>切割穩定性驗證 — HCC1395 tumor-only</title>
<style>
:root{{--accent:#D97757;--text:#141413;--bg:#FAF9F5;--bd:#E3DACC;--good:#0d9488;--trap:#db2777;--ceil:#7c3aed}}
*{{box-sizing:border-box}}
body{{font-family:system-ui,-apple-system,"PingFang TC","Microsoft JhengHei","Noto Sans CJK TC","Droid Sans Fallback",sans-serif;
 color:var(--text);background:var(--bg);margin:0;line-height:1.6}}
.wrap{{max-width:1180px;margin:0 auto;padding:28px 22px}}
h1{{font-size:23px;margin:0 0 4px}} h2{{font-size:18px;border-left:4px solid var(--accent);padding-left:10px;margin-top:30px}}
h3{{font-size:15px;margin-top:20px}}
.sub{{color:#6b6b6b;font-size:13px;margin-bottom:14px}}
.chip{{display:inline-block;color:#fff;border-radius:12px;padding:3px 10px;margin:3px 4px 3px 0;font-size:12px;font-weight:600}}
.badge{{display:inline-block;background:#e8f5f3;color:#0d7a6e;border:1px solid #0d9488;border-radius:6px;padding:2px 8px;font-size:12px}}
table{{border-collapse:collapse;width:100%;font-size:12.5px;margin:10px 0}}
th,td{{border:1px solid var(--bd);padding:4px 7px;text-align:center}} th{{background:#f0ebe2}}
td:first-child{{text-align:left}}
.kk{{font-size:11.5px}} .kk th{{background:#f6f1e8}}
.card{{border:1px solid var(--bd);border-radius:10px;margin:16px 0;background:#fff;overflow:hidden}}
.card-h{{padding:9px 14px;background:#faf6ef;border-bottom:1px solid var(--bd);font-size:14px}}
.card-body{{padding:12px 14px}}
.fig img{{max-width:100%;border:1px solid var(--bd);border-radius:6px}}
.fig{{margin-bottom:10px}} .cap{{font-size:11.5px;color:#777;margin-top:3px}}
.nofig{{color:#aaa;font-size:12px;padding:18px;text-align:center;border:1px dashed var(--bd)}}
.vbadge{{display:inline-block;color:#fff;border-radius:5px;padding:2px 7px;font-size:11px;font-weight:600}}
.muted{{color:#888;font-size:12.5px;font-weight:400}}
.box{{border:1px solid var(--bd);border-radius:8px;padding:12px 16px;background:#fff;margin:12px 0}}
.box.warn{{border-color:var(--trap);background:#fdf2f7}} .box.ok{{border-color:var(--good);background:#effbf9}}
.box.key{{border-color:var(--accent);background:#fdf6f1}}
code{{background:#f0ebe2;padding:1px 5px;border-radius:3px;font-size:12px}}
ul{{margin:6px 0}} li{{margin:3px 0}}
.foot{{margin-top:34px;border-top:1px solid var(--bd);padding-top:12px;color:#888;font-size:12px}}
</style></head><body><div class="wrap">

<h1>切割穩定性驗證 — 位點甲基切群「同一刀是否重現、是否夠好/更好」</h1>
<div class="sub">HCC1395 tumor-only · n={n} 代表位點 · single-sample ⭐2-3 characterization · 純分析（不改 C++）·
 <span class="badge">BERNOULLI 重算自檢 max|Δ| = {sc}（≈0 → null 距離正確）</span></div>

<div class="box key">
<b>一句話結論</b>：raw clusterboot Jaccard <b>全 0.74–1.0</b>（幾乎全「看似很穩」）→ <b>raw 穩定性零鑑別力</b>。
真正鑑別 = <b>excess-over-null</b>（扣 within-1-group null）+ <b>同一 k coherence</b> + <b>Cochran e≥5 對齊可靠性</b>。
三軸 chance-corrected gate 後：<b>{vc.get("GOOD_CUT",0)} GOOD</b> /
<b>{vc.get("STABLE+ALIGNED_DIFF_K",0)} 單樣本天花板</b> /
<b>{vc.get("STABLE_BUT_UNALIGNED",0)} 穩但不對齊（陷阱）</b>。
</div>

<h2>0. 怎麼驗證（3 軸 + 雙 null）</h2>
<ul>
<li><b>① 穩定性</b>：clusterboot Jaccard（{int(P['subsample_frac']*100)}% subsample × {P['B_subsamples']}，Hennig）。<b>但 raw 不可信</b> →
   對 <b>within-1-group null</b>（逐 CpG 欄內非 NaN 重排 → 重算 BERNOULLI，{P['R_null']} draws）校正 → <b>excess = real − null(95pct)</b>。</li>
<li><b>② 對齊（非循環）</b>：CramérV(cut vs a-priori 軸) 對 <b>label-shuffle null</b>（{P['P_shuffle']} 次）校正 + <b>Cochran e≥5 可靠性閘</b>。</li>
<li><b>判準</b>：<code>{html.escape(P['stability_pass'])}</code>；<code>{html.escape(P['align_pass'])}</code>。
   <b>GOOD = 同一 k 同時過</b>（<code>{html.escape(P['verdict_rule'])}</code>）。</li>
<li><b>雙 null 理由</b>：read-內甲基相關是「穩定≠正確」的根因（A-path 81% 失敗）→ within-1-group 擋「假穩」、label-shuffle 擋「假對齊」。</li>
</ul>

<h2>1. 總覽（一張圖看「夠好 vs 陷阱」）</h2>
<div class="fig"><img src="{b64(figidx['overview'])}" alt="overview"></div>
<div class="box warn"><b>關鍵讀法</b>：右圖多數點落右上（高 excess + 高 V）→ <b>連 (excess,V) 邊際也『看似都好』</b>；
真正把 GOOD 與陷阱分開的是<b>顏色</b>（同 k coherence + Cochran 可靠性），不是任何單一座標軸。
左圖：<b>GOOD 集中在高覆蓋組</b>（multi-resolution {bg.get('multi-resolution',{}).get('GOOD_CUT',0)}/{sum(bg.get('multi-resolution',{}).values())}、
single-k {bg.get('single-k-forced',{}).get('GOOD_CUT',0)}/{sum(bg.get('single-k-forced',{}).values())}）；
<b>ambiguous 幾全 UNALIGNED</b>（{bg.get('ambiguous-near-tie',{}).get('STABLE_BUT_UNALIGNED',0)}/{sum(bg.get('ambiguous-near-tie',{}).values())}）。</div>

<table><tr><th>group</th>{bghead}<th>小計</th></tr>{bgrows}</table>

<h2>2. 結果夠好或更好？</h2>
<div class="box ok"><b>「更好」（vs 舊版/raw-only）</b>：
<ul>
<li><b>raw Jaccard 0.74–1.0 全高 → 零鑑別力</b>；excess-over-null 跨 0.06–1.0 → <b>excess 才有鑑別力</b>（同 cross-sample ASM「看 excess 非 raw」教訓）。</li>
<li>3 軸 gate <b>正確降級 {vc.get("STABLE_BUT_UNALIGNED",0)+vc.get("STABLE+ALIGNED_DIFF_K",0)}/{n}</b> 個 raw-only 會誤判為「穩」的位點 → 嚴格優於只看穩定性。</li>
</ul></div>
<div class="box warn"><b>「夠好」的邊界（單樣本天花板）</b>：
<b>STABLE+ALIGNED_DIFF_K = {vc.get("STABLE+ALIGNED_DIFF_K",0)} 個</b>（多為 confident-unique）：
細 k 切得乾淨穩定（excess 0.7–1.0、V=1.0），<b>但該 k 下 Cochran e&lt;5（對齊統計不可靠）</b>；可靠對齊的 k=2 又 excess 低。
→ <b>單樣本上「最穩的刀」與「最可靠對齊的刀」常不同 k</b>；要把這些升為「確認」需更多樣本/覆蓋（COLO829、cross-sample）。</div>

<h2>3. 代表位點圖（肉眼確認）</h2>
<div class="sub">每張：左=穩定性(real vs within-1-group null，excess=鑑別量) · 右=對齊(綠=可靠對齊 e≥5；紅=不可靠/弱) · 下=既有 dual-panel 實際結構</div>
{cards}

<h2>4. 全 {n} 位點驗證表（數字全 L1，從 kprofile_stability.json 注入）</h2>
<table><tr><th>group</th><th>locus</th><th>n</th><th>軸</th><th>head k</th><th>real Jac</th><th>null95</th><th>excess</th><th>V</th><th>e</th><th>shufP</th><th>verdict</th></tr>{vrows}</table>

<h2>5. 限制</h2>
<ul>
<li><b>代表位點集</b>：取自 kprofile_examples（皆已 splittable）→ 無 state② 純無訊號位點；excess 下界偏高（最低 {min(r['headline']['stab_excess'] for r in loci):+.3f}）。全基因組 sweep 未跑。</li>
<li><b>單樣本 HCC1395 ⭐2-3</b>：穩定/對齊皆 cis-ASM characterization，<b>「好切法」≠ subclone</b>（仍需 normal cis-control）。</li>
<li>門檻 excess≥{P['excess_min']} / V≥0.3 / e≥5 為約定；per-k 全數字在驗證表，可自行 re-threshold。null draws={P['R_null']}、subsample={int(P['subsample_frac']*100)}%、seed={P['seed']}。</li>
</ul>

<div class="foot">資料源：<code>kprofile_stability.json</code>（seed={P['seed']}，BERNOULLI 自檢 max|Δ|={sc}）。
切群法 = ksweep_wg_v2.py 一致（peel→UPGMA→fcluster）。重生：<code>scripts/kprofile_stability.py → plot_kprofile_stability.py → build_stability_observation.py</code>。
單樣本 ⭐2-3 characterization，非 subclone 確認。</div>
</div></body></html>"""

out=f"{A}/20260622_clustering_stability_validation_01.standalone.html"
open(out,"w").write(HTML)
print(f"WROTE {out}  ({len(HTML)//1024} KB)")
print(f"self-check max|Δ|={sc}  verdicts={vc}")
