#!/usr/bin/env python3
"""
82 - Interactive capstone explorer: charts (caller-SEQC2 F1 3-stage, condition→FP OR,
CramersV×Δβ 2×2 scatter) + LIVE-adjustable filter sliders (CramersV/|Δβ|/tests/LOH) with
real-time TP/FP pass count + filterable/sortable table + click-locus→methylation heatmap.
Vanilla JS standalone. §13 layer-A.

Output: docs/experiments/in_progress/2026/06/20260608_capstone_interactive_explorer_01.standalone.html
"""
import os, csv, json, html, io, base64
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
       "20260608_capstone_interactive_explorer_01.standalone.html")
CFP = json.load(open(f"{CS}/condition_fp_consensus.json"))
FIGS = json.load(open(f"{CS}/capstone_loci_figs.json")) if os.path.exists(f"{CS}/capstone_loci_figs.json") else {}
TESTS = ["ClusterPermanovaP", "LabelHPPermanovaP", "LabelAllelePermanovaP", "HPFineP"]

# caller-SEQC2 3-stage benchmark (from benchmark_comparison.tsv)
BENCH = [("ClairS 原始", 0.9543, 0.7571, 0.8443, 1430),
         ("LongPhase-S", 0.9794, 0.7543, 0.8522, 627),
         ("InterSubMod", 0.9820, 0.7542, 0.8532, 544)]


def num(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def load(cls):
    rows = list(csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")))
    out = []
    for r in rows:
        out.append(dict(
            locus=f"{r.get('Chr')}:{r.get('Pos')}", cls=cls.upper(),
            cv=round(num(r, "CramersV") or 0.0, 3),
            db=round(num(r, "HPMergedDelta") or 0.0, 3),
            nt=sum(1 for t in TESTS if (num(r, t) is not None and num(r, t) <= 0.05)),
            loh=1 if (tb(r, "Potential_LOH") or tb(r, "LOH_Bed_Overlap")) else 0,
            sig=1 if tb(r, "Significant") else 0,
            reads=int(num(r, "NumReads") or 0),
        ))
    return out


ALLROWS = load("tp") + load("fp")
DATA = json.dumps(ALLROWS)
FIGJSON = json.dumps({k: v["fig"] for k, v in FIGS.items()})


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=125, bbox_inches="tight")
    plt.close(fig); return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def fig_bench():
    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    x = np.arange(3); names = [b[0] for b in BENCH]
    ax.bar(x - 0.27, [b[1] for b in BENCH], 0.25, label="Precision", color="#2563eb")
    ax.bar(x, [b[2] for b in BENCH], 0.25, label="Recall", color="#f59e0b")
    ax.bar(x + 0.27, [b[3] for b in BENCH], 0.25, label="F1", color="#15803d")
    for i, b in enumerate(BENCH):
        ax.text(i + 0.27, b[3] + 0.005, f"{b[3]:.3f}", ha="center", fontsize=7, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=8); ax.set_ylim(0.7, 1.0)
    ax.set_title("caller vs SEQC2 — 3 階段 (FP 過濾 1430→627→544)")
    ax.legend(fontsize=7, ncol=3); ax.spines[["top", "right"]].set_visible(False)
    return b64(fig)


def fig_or():
    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    conds = [c for c in CFP["A_condition_fp"] if c["fp_vs_tp_OR"] is not None]
    labs = [c["condition"].split(" (")[0].replace("_", " ") for c in conds]
    ors = [c["fp_vs_tp_OR"] for c in conds]
    cols = ["#dc2626" if o >= 2 else "#9ca3af" for o in ors]
    ax.barh(range(len(labs)), ors, color=cols)
    ax.axvline(1, color="#374151", ls="--", lw=1)
    for i, o in enumerate(ors):
        ax.text(o + 0.1, i, f"{o}", va="center", fontsize=8, fontweight="bold")
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=7.5)
    ax.set_xlabel("FP/TP OR (高Δβ 位點中)"); ax.set_title("條件→FP：無分群/LOH/小子群 聚集 FP，非低覆蓋")
    ax.spines[["top", "right"]].set_visible(False)
    return b64(fig)


def fig_scatter():
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    tp = [r for r in ALLROWS if r["cls"] == "TP"]
    fp = [r for r in ALLROWS if r["cls"] == "FP"]
    ax.scatter([abs(r["db"]) for r in tp], [r["cv"] for r in tp], s=4, alpha=0.18, c="#15803d", label="TP")
    ax.scatter([abs(r["db"]) for r in fp], [r["cv"] for r in fp], s=14, alpha=0.7, c="#dc2626", label="FP")
    ax.axhline(0.1, color="#1E3A8A", ls="--", lw=1)
    ax.axvline(0.1, color="#b45309", ls="--", lw=1)
    ax.set_xlabel("|Δβ| (平均偏移)"); ax.set_ylabel("CramersV (分群強度)")
    ax.set_title("2×2: CramersV × Δβ（右下=Δβ-only FP-prone 區，FP 紅點聚集）")
    ax.legend(fontsize=8); ax.set_xlim(-0.01, 0.6); ax.set_ylim(-0.02, 1.02)
    ax.spines[["top", "right"]].set_visible(False)
    return b64(fig)


F_BENCH, F_OR, F_SCAT = fig_bench(), fig_or(), fig_scatter()

HTML_DOC = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Capstone 互動 explorer — 圖表 + 可調篩選 + 位點圖</title>
<style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC",sans-serif;color:var(--ink);background:var(--bg);line-height:1.5}}
.wrap{{max-width:1200px;margin:0 auto;padding:20px 14px 80px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:12px;padding:18px 22px;margin-bottom:12px}}
header.top h1{{margin:0 0 4px;font-size:1.24rem}}.meta{{font-size:.8rem;opacity:.88}}
h2{{color:var(--accent);font-size:1.05rem;border-bottom:2px solid var(--line);padding-bottom:4px;margin-top:22px}}
.charts{{display:grid;grid-template-columns:1fr 1fr;gap:12px}}@media(max-width:820px){{.charts{{grid-template-columns:1fr}}}}
.chart{{background:#fff;border:1px solid var(--line);border-radius:10px;padding:8px}}.chart img{{width:100%}}
.panel{{position:sticky;top:0;z-index:5;background:#fff;border:1px solid var(--line);border-radius:10px;padding:12px 14px;margin:12px 0;box-shadow:0 2px 8px rgba(0,0,0,.05)}}
.sliders{{display:flex;gap:18px;flex-wrap:wrap;align-items:center;font-size:.82rem}}
.sliders label{{display:flex;flex-direction:column;gap:2px;font-size:.74rem;color:var(--mut)}}
.sliders input[type=range]{{width:130px}}
.live{{display:flex;gap:14px;flex-wrap:wrap;margin-top:8px;font-size:.85rem}}
.kpi{{background:var(--bg);border-radius:8px;padding:6px 12px}}.kpi b{{color:var(--accent);font-size:1.1rem}}
table{{width:100%;border-collapse:collapse;font-size:.8rem;background:#fff}}
th,td{{padding:4px 7px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}}th{{background:#f1f5f9;font-size:.72rem;cursor:pointer;position:sticky;top:0}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}.mono{{font-family:ui-monospace,monospace;font-size:.74rem}}
.tag{{font-size:.68rem;padding:1px 6px;border-radius:5px;color:#fff}}.tp{{background:#15803d}}.fp{{background:#dc2626}}
.pass{{color:#15803d;font-weight:700}}.fail{{color:#9ca3af}}
.haspic{{color:#1E3A8A;cursor:zoom-in;text-decoration:underline}}
.tablewrap{{max-height:62vh;overflow:auto;border:1px solid var(--line);border-radius:10px}}
#modal{{display:none;position:fixed;inset:0;background:rgba(0,0,0,.82);z-index:50;align-items:center;justify-content:center;cursor:zoom-out}}
#modal img{{max-width:92vw;max-height:88vh;border:3px solid #fff;border-radius:6px}}
.note{{font-size:.76rem;color:var(--mut);margin:4px 0}}
</style></head><body><div class="wrap">

<header class="top"><h1>Capstone 互動 explorer — 圖表 + 可調篩選 + 點位點看甲基圖</h1>
<div class="meta">HCC1395 全 TP+FP（{len(ALLROWS):,} 位點）· 拖滑桿即時看通過數/TP-FP 比 · 點有圖位點放大看 read×CpG 甲基熱圖 · 2026-06-08 · §13 layer-A</div></header>

<h2>① 數據圖表</h2>
<div class="charts">
  <div class="chart"><img src="{F_BENCH}"></div>
  <div class="chart"><img src="{F_OR}"></div>
</div>
<div class="chart" style="margin-top:12px"><img src="{F_SCAT}"></div>

<h2>② 可調整篩選 + 即時觀察（拖滑桿）</h2>
<div class="panel">
  <div class="sliders">
    <label>CramersV ≥ <span id="vcv">0.10</span><input id="cv" type="range" min="0" max="1" step="0.05" value="0.1"></label>
    <label>|Δβ| ≥ <span id="vdb">0.10</span><input id="db" type="range" min="0" max="0.5" step="0.02" value="0.1"></label>
    <label>≥ N/4 檢定 <span id="vnt">2</span><input id="nt" type="range" min="0" max="4" step="1" value="2"></label>
    <label>聯集模式<select id="mode"><option value="or">CramersV OR Δβ</option><option value="and">兩者皆需</option><option value="cv">只 CramersV</option><option value="db">只 Δβ</option></select></label>
    <label>排除 LOH<input id="exloh" type="checkbox"></label>
  </div>
  <div class="live" id="live"></div>
  <div class="note">通過 = (依模式 CramersV/Δβ 門檻) AND ≥N/4 檢定 AND (排除 LOH 勾選時非 LOH)。比 = TP通過率 / FP通過率（&gt;1 才有判別力）。</div>
</div>

<div class="tablewrap"><table id="tbl">
<thead><tr><th data-k="locus">locus</th><th data-k="cls">類</th><th data-k="cv" class="num">CramersV</th><th data-k="db" class="num">Δβ</th>
  <th data-k="nt" class="num">檢定</th><th data-k="reads" class="num">reads</th><th data-k="loh">LOH</th><th data-k="sig">ISM-Sig</th><th data-k="pass">通過</th><th>圖</th></tr></thead>
<tbody id="tb"></tbody></table></div>
<div class="note">表只列「通過」的位點（依當前滑桿），最多 300 列依 |Δβ|+CramersV 排序。藍色「看圖」可點開甲基熱圖（{len(FIGS)} 個 curated 位點有圖）。</div>

<div id="modal"><img id="mimg"></div>

<script>
const DATA={DATA};
const FIGS={FIGJSON};
const $=id=>document.getElementById(id);
let sortK="db",sortDir=-1;
function passes(r,cv,db,nt,mode,exloh){{
  if(r.nt<nt)return false;
  if(exloh&&r.loh)return false;
  const c=r.cv>=cv, d=Math.abs(r.db)>=db;
  if(mode==="or")return c||d;
  if(mode==="and")return c&&d;
  if(mode==="cv")return c;
  return d;
}}
function update(){{
  const cv=+$("cv").value, db=+$("db").value, nt=+$("nt").value, mode=$("mode").value, exloh=$("exloh").checked;
  $("vcv").textContent=cv.toFixed(2);$("vdb").textContent=db.toFixed(2);$("vnt").textContent=nt;
  let tpN=0,tpP=0,fpN=0,fpP=0;
  const sel=[];
  for(const r of DATA){{
    const p=passes(r,cv,db,nt,mode,exloh);
    if(r.cls==="TP"){{tpN++;if(p)tpP++;}}else{{fpN++;if(p)fpP++;}}
    if(p)sel.push(r);
  }}
  const tpr=tpP/tpN, fpr=fpP/fpN, ratio=fpr>0?(tpr/fpr):Infinity;
  $("live").innerHTML=
    `<div class="kpi">TP 通過 <b>${{tpP}}</b>/${{tpN}} (${{(tpr*100).toFixed(2)}}%)</div>`+
    `<div class="kpi">FP 通過 <b>${{fpP}}</b>/${{fpN}} (${{(fpr*100).toFixed(2)}}%)</div>`+
    `<div class="kpi">TP/FP 判別比 <b>${{isFinite(ratio)?ratio.toFixed(2):'∞'}}</b> ${{ratio>=2?'✅':(ratio<1.3?'⚠ 弱':'')}}</div>`+
    `<div class="kpi">總通過 <b>${{sel.length}}</b></div>`;
  sel.sort((a,b)=>sortDir*((sortK==="db"?Math.abs(b.db)+b.cv:b[sortK])>(sortK==="db"?Math.abs(a.db)+a.cv:a[sortK])?1:-1));
  const tb=$("tb");tb.innerHTML="";
  const fr=document.createDocumentFragment();
  for(const r of sel.slice(0,300)){{
    const tr=document.createElement("tr");
    const pic=FIGS[r.locus]?`<span class="haspic" data-l="${{r.locus}}">看圖🔬</span>`:"—";
    tr.innerHTML=`<td class="mono">${{r.locus}}</td><td><span class="tag ${{r.cls.toLowerCase()}}">${{r.cls}}</span></td>`+
      `<td class="num">${{r.cv.toFixed(2)}}</td><td class="num">${{(r.db>=0?'+':'')+r.db.toFixed(2)}}</td>`+
      `<td class="num">${{r.nt}}</td><td class="num">${{r.reads}}</td><td>${{r.loh?'LOH':'—'}}</td>`+
      `<td class="${{r.sig?'pass':'fail'}}">${{r.sig?'✓':'·'}}</td><td class="pass">✓</td><td>${{pic}}</td>`;
    fr.appendChild(tr);
  }}
  tb.appendChild(fr);
}}
["cv","db","nt","mode","exloh"].forEach(id=>$(id).oninput=update);
document.querySelectorAll("#tbl th").forEach(th=>th.onclick=()=>{{const k=th.dataset.k;if(k==="pass"||!k)return;if(sortK===k)sortDir*=-1;else{{sortK=k;sortDir=-1;}}update();}});
document.addEventListener("click",e=>{{if(e.target.classList.contains("haspic")){{$("mimg").src=FIGS[e.target.dataset.l];$("modal").style.display="flex";}}}});
$("modal").onclick=()=>$("modal").style.display="none";
update();
</script>
</div></body></html>'''

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write(HTML_DOC)
print(f"[82] wrote {OUT} ({len(HTML_DOC)//1024} KB, {len(ALLROWS)} loci, {len(FIGS)} figs)")
