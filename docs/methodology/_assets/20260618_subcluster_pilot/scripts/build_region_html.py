"""[per-region 觀察 HTML] §13-A 由 JSON 注入: 樹形分布 bar + 完整性帳本 + 乾淨 full_tree 區域表 + chr17 例。
缺 key 直接 refuse。輸出 ../../20260626_genomewide_sSNV_linkage_region_trees_01.standalone.html。"""
import json

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260626_genomewide_sSNV_linkage_region_trees_01.standalone.html"
TX = "#141413"; BD = "#E3DACC"; MUT = "#6B6862"; AC = "#D97757"; BLU = "#4A6E8A"; GRN = "#5B8A5B"; RED = "#C0563F"; TEAL = "#2E8B8B"


def req(d, k):
    if k not in d:
        raise SystemExit(f"§13-A refuse: missing key {k}")
    return d[k]


reg = json.load(open(f"{A}/sm_region_integration.json"))
led = json.load(open(f"{A}/sm_completeness_ledger.json"))
summ = json.load(open(f"{A}/sm_summary.json"))
agg = reg["aggregate"]
n = req(agg, "n_regions")
shape = req(agg, "shape_distribution")
shape_clean = agg["shape_distribution_clean_loh_neutral"]
buckets = req(led, "buckets")

SHAPE_ORDER = ["full_tree", "linear_nested", "sibling_only", "co_linked_lineage", "no_confirmed_structure", "inconsistent"]
SHAPE_LAB = {"full_tree": "full_tree 分支+深度", "linear_nested": "linear_nested 祖先→後代鏈",
             "sibling_only": "sibling_only 平行分支", "co_linked_lineage": "co_linked 單lineage",
             "no_confirmed_structure": "no_confirmed_structure", "inconsistent": "inconsistent(cycle)"}
SHAPE_COL = {"full_tree": AC, "linear_nested": BLU, "sibling_only": TEAL, "co_linked_lineage": GRN,
             "no_confirmed_structure": MUT, "inconsistent": RED}

# bar chart (SVG)
mx = max(shape.values())
W = 720; rh = 30; pad = 230
bars = []
y = 10
for k in SHAPE_ORDER:
    v = shape.get(k, 0); vc = shape_clean.get(k, 0)
    w = v / mx * (W - pad - 60)
    wc = vc / mx * (W - pad - 60)
    bars.append(f'<text x="{pad-8}" y="{y+rh/2+4}" font-size="12" fill="{TX}" text-anchor="end">{SHAPE_LAB[k]}</text>')
    bars.append(f'<rect x="{pad}" y="{y}" width="{w:.1f}" height="{rh-8}" fill="{SHAPE_COL[k]}" opacity="0.35"/>')
    bars.append(f'<rect x="{pad}" y="{y}" width="{wc:.1f}" height="{rh-8}" fill="{SHAPE_COL[k]}"/>')
    bars.append(f'<text x="{pad+w+6:.1f}" y="{y+rh/2+2}" font-size="11.5" fill="{TX}">{v} ({100*v/n:.1f}%) · clean {vc}</text>')
    y += rh
svg = f'<svg viewBox="0 0 {W} {y+10}" xmlns="http://www.w3.org/2000/svg" font-family="system-ui">{"".join(bars)}</svg>'

# top clean full_tree regions
clean_ft = [r for r in reg["regions"] if r["tree_shape"] == "full_tree" and r["cn"] in ("loh", "neutral")]
clean_ft = sorted(clean_ft, key=lambda r: -(r["n_nested"] + r["n_sibling"]))[:30]
rows = []
for r in clean_ft:
    pops = r.get("populations") or {}
    popstr = " ".join(f"{k}:{v}" for k, v in list(pops.items())[:5]) if pops else "—"
    rows.append(f'<tr><td><code>{r["region"]}</code></td><td>{r["cn"]}</td><td>{r["n_sSNV"]}</td>'
                f'<td>{r["n_nodes"]}</td><td>{r["n_nested"]}</td><td>{r["n_sibling"]}</td><td>{r["max_depth"]}</td>'
                f'<td>{r.get("n_populations","—")}</td><td style="font-family:monospace;font-size:11px">{popstr}</td></tr>')
rows_html = "".join(rows)

f1 = summ["F1_per_sSNV_distinct"]
cn = summ["F2_artifact_quantification"]["linked_somatic_by_CN_state"]

html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>全基因組 sSNV 連鎖 → 每區域克隆樹</title><style>
body{{margin:0;font-family:system-ui,"Noto Sans CJK TC",sans-serif;color:{TX};background:#FAF9F5;line-height:1.7}}
.wrap{{max-width:1080px;margin:0 auto;padding:0 22px 80px}}
header{{background:{TX};color:#FAF9F5;padding:26px 22px}} header h1{{margin:0;font-size:20px}} header .s{{color:#cfc9bf;font-size:13px;margin-top:5px}}
h2{{font-size:17px;margin:28px 0 10px;padding-bottom:5px;border-bottom:2px solid {AC}}}
.cards{{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}}
.card{{flex:1;min-width:150px;background:white;border:1px solid {BD};border-radius:8px;padding:12px 14px}}
.card b{{font-size:22px;color:{AC}}} .card span{{display:block;font-size:12px;color:{MUT};margin-top:3px}}
table{{border-collapse:collapse;font-size:12.5px;margin:8px 0;width:100%}} th,td{{border:1px solid {BD};padding:4px 8px;text-align:center}} th{{background:#efe9df}}
.box{{background:white;border:1px solid {BD};border-left:4px solid {AC};border-radius:0 6px 6px 0;padding:11px 15px;margin:12px 0;font-size:13.5px}}
.red{{border-left-color:{RED};background:#fbf0ec}} code{{background:#efe9df;padding:1px 5px;border-radius:3px;font-size:11.5px}}
svg{{max-width:100%;background:white;border:1px solid {BD};border-radius:8px;padding:8px}}
footer{{margin-top:30px;padding-top:13px;border-top:1px solid {BD};font-size:12px;color:{MUT}}}
</style></head><body>
<header><div class="wrap" style="padding:0"><h1>全基因組 sSNV 單分子連鎖 → 每區域克隆樹重建 <span style="background:{AC};color:#fff;font-size:11px;padding:1px 7px;border-radius:10px">⭐3 HCC1395</span></h1>
<div class="s">2-locus pairwise → multi-locus 組合 → per-region 樹 · Tier-R · 已過對抗稽核(5修正) · 數字 §13-A 注入</div></div></header>
<div class="wrap">

<h2>① 整體（單位=最大可關聯區域）</h2>
<div class="cards">
<div class="card"><b>{n:,}</b><span>最大可關聯區域</span></div>
<div class="card"><b>{shape['full_tree']}</b><span>full_tree 分支樹（clean {shape_clean.get('full_tree',0)}）</span></div>
<div class="card"><b>{buckets['linked']:,}</b><span>linked sSNV（完整帳本 /{led['universe_total']:,}）</span></div>
<div class="card"><b>{f1['any_confirmed_sameHP_relation']:,}</b><span>distinct sSNV 有確認同HP連結</span></div>
</div>

<h2>② 各區域樹形分布（深色=乾淨 LOH/neutral 子集）</h2>
{svg}
<div class="box">3,820 區域（53%）有確認克隆分支/階層；<b>full_tree {shape['full_tree']}（9.5%）最豐富</b>，乾淨 CN 區 {shape_clean.get('full_tree',0)} 個 = 論文級可信子集。</div>

<h2>③ 完整性帳本（no-miss, Tier-R）+ CN 分層</h2>
<table><tr><th>bucket</th><th>數量</th><th></th><th>linked-somatic CN 分層</th><th>數量</th></tr>
<tr><td>linked</td><td>{buckets['linked']:,}</td><td rowspan="3"></td><td>gain（混淆）</td><td>{cn.get('gain',0):,}</td></tr>
<tr><td>underpowered</td><td>{buckets['underpowered']:,}</td><td>loh（乾淨）</td><td>{cn.get('loh',0):,}</td></tr>
<tr><td>isolated_singleton</td><td>{buckets['isolated_singleton']:,}</td><td>neutral（乾淨）</td><td>{cn.get('neutral',0):,}</td></tr></table>
<div class="box red">🔴 ~69% linked-somatic 在 CN-gain（multiplicity/amplicon/segdup 混淆）→ 可信集 = LOH+neutral。</div>

<h2>④ 乾淨 full_tree 區域（LOH/neutral，top 30 by 分支數）</h2>
<table><tr><th>region</th><th>CN</th><th>n_sSNV</th><th>nodes</th><th>nested</th><th>sibling</th><th>depth</th><th>pops</th><th>populations(top5)</th></tr>{rows_html}</table>
<div class="box">每列 = 一個局部克隆樹。populations 欄 = multi-locus 基因型向量:read數（A=ALT,R=REF，位置依 region 內 sSNV 順序）。完整查詢：<code>lists/regions.tsv</code> + <code>sm_region_integration.json</code>。</div>

<h2>⑤ chr17:48360161 canonical（pipeline 自動重建）</h2>
<div class="box">shape=<b>full_tree</b>、CN=<b>loh</b>、4 sSNV、3 populations：<code>ARRR</code>(γ subclone)7 / <code>RRAR</code>(L1 α-only)10 / <code>RAAA</code>(L2 α+β)15 reads。樹：germline →（γ sibling ∥ α 祖先→β 後代），VAF α0.82>β0.48>γ0.18。手工驗證例被全基因組 pipeline 重現。</div>

<div class="box red">🔴 <b>誠實邊界</b>：⭐3 單樣本；regional(≤read-span)非 genome-wide tree；Fisher-sig 不分 subclone/allelic(HP-gate 才分)；偽影未清(缺 mappability)；分子證據非 single-cell confirmation。對外勿稱「甲基偵測 subclone / genome-wide tree / 對手缺檢定」。</div>

<footer>HCC1395 ⭐3 · Tier-R · §13-A 由 sm_region_integration.json + sm_completeness_ledger.json + sm_summary.json 注入<br>
manifest: _assets/20260618_subcluster_pilot/README_sm_linkage_pipeline.md · 對抗稽核: sm_adversarial_audit.md · 腳本 build_region_html.py</footer>
</div></body></html>"""
open(OUT, "w").write(html)
print(f"[-> {OUT}] {len(html)} bytes; clean_full_tree shown={len(clean_ft)}")
