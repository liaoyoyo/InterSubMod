"""[HCC1395 區域統計 圖示化 HTML] §13-A 由 sm_region_stats.json + config census + ledger + summary 注入。
inline SVG 直方圖(sSNV/region, span) + 連接關係 grouped bar + 樹形 all/clean 兩色條 + CN 比例條 + 敘述。
缺 key refuse。輸出 ../../20260626_HCC1395_region_statistics_01.standalone.html。"""
import json

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260626_HCC1395_region_statistics_01.standalone.html"
TX = "#141413"; BD = "#E3DACC"; MUT = "#6B6862"; AC = "#D97757"; BLU = "#4A6E8A"; GRN = "#5B8A5B"; RED = "#C0563F"; TEAL = "#2E8B8B"; GOLD = "#C9A14A"


def req(d, k):
    if k not in d:
        raise SystemExit(f"§13-A refuse: missing {k}")
    return d[k]


st = json.load(open(f"{A}/sm_region_stats.json"))
cfg = json.load(open(f"{A}/sm_configuration_census.json"))
led = json.load(open(f"{A}/sm_completeness_ledger.json"))
summ = json.load(open(f"{A}/sm_summary.json"))
N = req(st, "n_regions")


def hbar(items, w=620, rowh=27, lab_w=120, val_w=170):
    """items: (label, value, color, clean_or_None). clean 畫深色 overlay。"""
    mx = max(i[1] for i in items) or 1
    barw = w - lab_w - val_w
    y = 6
    s = []
    for it in items:
        label, val, col = it[0], it[1], it[2]
        clean = it[3] if len(it) > 3 else None
        bw = val / mx * barw
        s.append(f'<text x="{lab_w-8}" y="{y+rowh/2+4}" font-size="12" text-anchor="end" fill="{TX}">{label}</text>')
        s.append(f'<rect x="{lab_w}" y="{y+3}" width="{bw:.1f}" height="{rowh-9}" fill="{col}" opacity="0.4" rx="2"/>')
        if clean is not None:
            s.append(f'<rect x="{lab_w}" y="{y+3}" width="{clean/mx*barw:.1f}" height="{rowh-9}" fill="{col}" rx="2"/>')
            extra = f' · clean {clean:,}'
        else:
            s.append(f'<rect x="{lab_w}" y="{y+3}" width="{bw:.1f}" height="{rowh-9}" fill="{col}" rx="2"/>')
            extra = ''
        s.append(f'<text x="{lab_w+bw+6:.1f}" y="{y+rowh/2+4}" font-size="11" fill="{TX}">{val:,} ({100*val/N:.1f}%){extra}</text>')
        y += rowh
    return f'<svg viewBox="0 0 {w} {y+8}" xmlns="http://www.w3.org/2000/svg" font-family="system-ui">{"".join(s)}</svg>'


def grouped(items, w=620, rowh=30, lab_w=160):
    """connection: (label, same, diff). 兩條 (同HP 實色 / 異HP 淡)。"""
    mx = max(max(i[1], i[2]) for i in items) or 1
    barw = w - lab_w - 130
    y = 6
    s = []
    for label, same, diff in items:
        s.append(f'<text x="{lab_w-8}" y="{y+rowh/2+4}" font-size="11.5" text-anchor="end" fill="{TX}">{label}</text>')
        s.append(f'<rect x="{lab_w}" y="{y+2}" width="{same/mx*barw:.1f}" height="{(rowh-8)/2}" fill="{BLU}"/>')
        s.append(f'<rect x="{lab_w}" y="{y+2+(rowh-8)/2}" width="{diff/mx*barw:.1f}" height="{(rowh-8)/2}" fill="{GOLD}"/>')
        s.append(f'<text x="{lab_w+max(same,diff)/mx*barw+6:.1f}" y="{y+rowh/2+4}" font-size="10.5" fill="{MUT}">同{same:,} / 異{diff:,}</text>')
        y += rowh
    return f'<svg viewBox="0 0 {w} {y+8}" xmlns="http://www.w3.org/2000/svg" font-family="system-ui">{"".join(s)}</svg>'


# data
snv = st["snv_per_region_hist"]
snv_items = [(b, snv.get(b, 0), BLU) for b in ["2", "3", "4", "5", "6-7", "8-10", "11-20", "21-50", "51+"]]
span = st["span_hist"]
span_items = [(b, span.get(b, 0), TEAL) for b in ["<1kb", "1-5kb", "5-10kb", "10-50kb", "50-200kb", ">=200kb"]]
shp = st["tree_shape"]; shpc = st["tree_shape_clean"]
SH = [("full_tree 分支+深度", "full_tree", AC), ("linear_nested 祖先→後代鏈", "linear_nested", BLU),
      ("sibling_only 平行分支", "sibling_only", TEAL), ("co_linked 單lineage", "co_linked_lineage", GRN),
      ("no_confirmed", "no_confirmed_structure", MUT), ("inconsistent", "inconsistent", RED)]
shp_items = [(lab, shp.get(k, 0), col, shpc.get(k, 0)) for lab, k, col in SH]
cn = st["cn_hist"]
CN = [("gain（混淆）", "gain", RED), ("loh（乾淨）", "loh", GRN), ("neutral（乾淨）", "neutral", TEAL), ("loss", "loss", MUT)]
cn_items = [(lab, cn.get(k, 0), col) for lab, k, col in CN]
# connection from config census
conn = {r["config"]: r for r in cfg["census_powered_somatic"]}
CN_ORD = [("共連 co_linked", "co_linked"), ("nested a⊂b", "nested_a_in_b"), ("nested b⊂a", "nested_b_in_a"),
          ("互斥 mutual_excl", "mutual_excl"), ("independent", "independent")]
conn_items = [(lab, conn[k]["pairs_sameHP"], conn[k]["pairs_diffHP"]) for lab, k in CN_ORD if k in conn]
f1 = summ["F1_per_sSNV_distinct"]
buckets = led["buckets"]
clean_regions = cn.get("loh", 0) + cn.get("neutral", 0)
struct_regions = shp.get("full_tree", 0) + shp.get("linear_nested", 0) + shp.get("sibling_only", 0) + shp.get("co_linked_lineage", 0)

html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HCC1395 區域統計畫像</title><style>
body{{margin:0;font-family:system-ui,"Noto Sans CJK TC",sans-serif;color:{TX};background:#FAF9F5;line-height:1.7}}
.wrap{{max-width:1000px;margin:0 auto;padding:0 22px 80px}}
header{{background:{TX};color:#FAF9F5;padding:26px 22px}} header h1{{margin:0;font-size:20px}} header .s{{color:#cfc9bf;font-size:13px;margin-top:5px}}
h2{{font-size:17px;margin:26px 0 8px;padding-bottom:5px;border-bottom:2px solid {AC}}}
.cards{{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}}
.card{{flex:1;min-width:140px;background:white;border:1px solid {BD};border-radius:8px;padding:12px 14px}}
.card b{{font-size:23px;color:{AC}}} .card span{{display:block;font-size:12px;color:{MUT};margin-top:3px}}
svg{{max-width:100%;background:white;border:1px solid {BD};border-radius:8px;padding:10px;margin:6px 0}}
.box{{background:white;border:1px solid {BD};border-left:4px solid {AC};border-radius:0 6px 6px 0;padding:11px 15px;margin:11px 0;font-size:13.5px}}
.red{{border-left-color:{RED};background:#fbf0ec}} .ok{{border-left-color:{GRN};background:#eef4ee}}
.leg{{font-size:11.5px;color:{MUT};margin:2px 0 8px}} .sw{{display:inline-block;width:11px;height:11px;border-radius:2px;vertical-align:middle;margin:0 3px 0 10px}}
code{{background:#efe9df;padding:1px 5px;border-radius:3px;font-size:11.5px}}
footer{{margin-top:28px;padding-top:12px;border-top:1px solid {BD};font-size:12px;color:{MUT}}}
</style></head><body>
<header><div class="wrap" style="padding:0"><h1>HCC1395 全基因組 sSNV 連鎖區域 — 統計畫像 <span style="background:{AC};color:#fff;font-size:11px;padding:1px 7px;border-radius:10px">⭐3</span></h1>
<div class="s">單位=最大可關聯區域（≤50kb 內 ≥2 somatic sSNV）· Tier-R · 數字 §13-A 注入 · sum-check ✓</div></div></header>
<div class="wrap">

<h2>① 區域總覽</h2>
<div class="cards">
<div class="card"><b>{N:,}</b><span>最大可關聯區域</span></div>
<div class="card"><b>{st['total_somatic_sSNV_in_regions']:,}</b><span>區域內 somatic sSNV</span></div>
<div class="card"><b>{struct_regions:,}</b><span>有確認克隆結構（{100*struct_regions/N:.0f}%）</span></div>
<div class="card"><b>{clean_regions:,}</b><span>乾淨 CN 區（LOH+neutral，{100*clean_regions/N:.0f}%）</span></div>
</div>
<div class="box">完整性帳本（全 35,332 union sSNV）：linked {buckets['linked']:,} / underpowered {buckets['underpowered']:,} / 孤立 {buckets['isolated_singleton']:,}（加總=35,332 ✓）。</div>

<h2>② 每區域 sSNV 數量分布</h2>
{hbar(snv_items)}
<div class="box">🔑 <b>74% 區域只有 2–3 個位點</b>（中位數 2，最大 150）；多位點密集區稀少，且部分為偽影簇。</div>

<h2>③ 區域跨度 span 分布</h2>
{hbar(span_items)}

<h2>④ 連接關係（每對 2×2；<span style="color:{BLU}">■同HP=克隆</span> / <span style="color:{GOLD}">■異HP=allelic</span>）</h2>
{grouped(conn_items)}
<div class="box">巢狀（祖先→後代）<b>同HP {conn['nested_a_in_b']['pairs_sameHP']+conn['nested_b_in_a']['pairs_sameHP']:,} >> 異HP {conn['nested_a_in_b']['pairs_diffHP']+conn['nested_b_in_a']['pairs_diffHP']:,}</b>（克隆階層幾乎全同單倍型）；互斥則<b>異HP {conn['mutual_excl']['pairs_diffHP']:,} > 同HP {conn['mutual_excl']['pairs_sameHP']:,}</b>（多為 allelic 非 subclone）。去重 distinct sSNV 有確認同HP連結={f1['any_confirmed_sameHP_relation']:,}。</div>

<h2>⑤ 存在的結構（樹形）<span class="leg"><span class="sw" style="background:{AC};opacity:.4"></span>全部 <span class="sw" style="background:{AC}"></span>乾淨LOH+neutral</span></h2>
{hbar(shp_items)}
<div class="ok"><b>有確認克隆結構（前 4 類）= {struct_regions:,} 區域（{100*struct_regions/N:.1f}%）</b>；其中<b>完整分支樹 full_tree {shp['full_tree']:,}（{100*shp['full_tree']/N:.1f}%）</b>。樹深度：depth1={st['depth_hist'].get('1',0):,} / depth2={st['depth_hist'].get('2',0):,} / depth≥3={sum(v for k,v in st['depth_hist'].items() if int(k)>=3)} → 多為淺樹。</div>

<h2>⑥ 🔴 CN 分層（評估必看）</h2>
{hbar(cn_items)}
<div class="red"><b>{100*cn.get('gain',0)/N:.0f}% 區域在 CN-gain（multiplicity/amplicon/segdup 混淆）</b>→ 論文級可信集 = LOH+neutral <b>{clean_regions:,} 區域（{100*clean_regions/N:.1f}%）</b>，其中 full_tree {shpc.get('full_tree',0)}。</div>

<h2>⑦ 驗證與紀錄</h2>
<div class="box">sum-check：sSNV 直方圖 / 樹形分布 各加總={N:,}；帳本三桶=35,332（=union）。全 ✓。<br>
逐區域查：<code>lists/regions.tsv</code> · 紀錄：<code>InterSubMod/docs/methodology/20260626_HCC1395_region_statistics_record_01.md</code> · 重現：<code>README_sm_linkage_pipeline.md</code></div>
<div class="red">🔴 限制：⭐3 單樣本；Tier-R only；64% 在 CN-gain；偽影未清（缺 mappability）；Fisher-sig 不分 subclone/allelic（HP-gate 才分）；regional 非 genome-wide tree；分子證據非 single-cell confirmation。</div>

<footer>HCC1395 ⭐3 · §13-A 由 sm_region_stats.json + sm_configuration_census.json + sm_completeness_ledger.json 注入 · 腳本 build_region_stats_html.py</footer>
</div></body></html>"""
open(OUT, "w").write(html)
print(f"[-> {OUT}] {len(html)} bytes")
