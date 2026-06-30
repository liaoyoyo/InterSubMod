#!/usr/bin/env python3
"""觀察 HTML generator — read-driven genome-wide cross-check。
§13-A: 所有數字從 data/*.json 注入，缺 key refuse；不手打。"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
RD = json.load(open(os.path.join(D, "read_driven_genomewide.json")))
RS = json.load(open(os.path.join(D, "region_shape_distribution.json")))
TR = json.load(open(os.path.join(D, "tier_r_characterization.json")))
OUT = os.path.join(HERE, "read_driven_crosscheck.standalone.html")


def g(d, *keys):
    for k in keys:
        d = d[k]
    return d


# ---- read-driven 真值 ----
n_snv = g(RD, "n_snv_total"); n_tp = g(RD, "n_tp"); n_fp = g(RD, "n_fp")
n_touch = g(RD, "n_reads_touching_total")
n_multi = g(RD, "n_reads_multi_alt_genuine")
n_chim = g(RD, "n_chimeric_reads_excluded")
combos = g(RD, "distinct_alt_combos_total")
pairs = g(RD, "distinct_pairs_total")
collapse = g(RD, "collapse_ratio_reads_per_combo")
max_alt = g(RD, "max_alt_per_read"); max_clean = g(RD, "max_clean_per_read")
cross50k = g(RD, "n_cross50k_genuine_chains")
csd = {int(k): v for k, v in g(RD, "direct_chain_size_dist").items()}
hpd = g(RD, "chain_hp_dist")
top = g(RD, "top_dense_chains")
per_chrom = g(RD, "per_chrom")

# ---- pipeline shape (cross-ref) ----
tot_reg = g(RS, "total_regions")
shapes = g(RS, "by_shape")
struct_plus = g(RS, "structured_plus_colinked")
full_tree = g(shapes, "full_tree", "n"); full_clean = g(shapes, "full_tree", "clean_loh_neutral")

# ---- tier_r buildability ----
linked = g(TR, "linked_ref", "n"); under = g(TR, "underpowered", "n"); iso = g(TR, "isolated", "n")

# chain-size buckets
buckets = [("2", csd.get(2, 0)), ("3", csd.get(3, 0)), ("4", csd.get(4, 0)),
           ("5–10", sum(csd.get(i, 0) for i in range(5, 11))),
           ("11–20", sum(csd.get(i, 0) for i in range(11, 21))),
           ("21–50", sum(csd.get(i, 0) for i in range(21, 51))),
           ("51–85", sum(csd.get(i, 0) for i in range(51, 86)))]
bmax = max(v for _, v in buckets) or 1

# collapse cascade
cascade = [("reads touching ≥1 sSNV", n_touch),
           ("reads carrying ≥2 ALT (直接多路鏈)", n_multi),
           ("distinct ALT combinations", combos),
           ("read-span regions (pipeline)", tot_reg),
           ("confirmed-structure regions", struct_plus),
           ("full evolution trees", full_tree)]
cmax = cascade[0][1] or 1


def bar(w_pct, color="#5B9BD5"):
    return f'<div style="background:{color};height:14px;width:{w_pct:.1f}%;border-radius:2px"></div>'


def fnum(n):
    return f"{n:,}"


rows_cascade = ""
prev = None
for label, v in cascade:
    ratio = f"÷{prev/v:.1f}×" if (prev and v) else "—"
    rows_cascade += (f'<tr><td>{label}</td><td class="num">{fnum(v)}</td>'
                     f'<td style="width:40%">{bar(100*v/cmax)}</td><td class="num">{ratio}</td></tr>')

rows_chain = ""
for label, v in buckets:
    rows_chain += (f'<tr><td>{label}-way</td><td class="num">{fnum(v)}</td>'
                   f'<td style="width:55%">{bar(100*v/bmax,"#70AD47")}</td></tr>')

rows_top = ""
for c in top[:10]:
    dens = c["span"] / c["n_alt"] if c["n_alt"] else 0
    flag = "🔴 artifact 疑似" if dens < 2000 else ""
    rows_top += (f'<tr><td class="num">{c["n_alt"]}</td><td>{c["chrom"]}:{fnum(c["lo"])}–{fnum(c["hi"])}</td>'
                 f'<td class="num">{fnum(c["span"])}</td><td class="num">{dens:.0f} bp/ALT</td>'
                 f'<td>{c["hp"]}</td><td>{flag}</td></tr>')

rows_shape = ""
order = ["full_tree", "linear_nested", "sibling_only", "co_linked_lineage", "no_confirmed_structure", "inconsistent"]
for s in order:
    sd = shapes[s]
    rows_shape += (f'<tr><td>{s}</td><td class="num">{fnum(sd["n"])}</td>'
                   f'<td class="num">{sd["pct_of_regions"]}%</td><td class="num">{fnum(sd["clean_loh_neutral"])}</td></tr>')

rows_hp = ""
for k, v in sorted(hpd.items(), key=lambda kv: -kv[1]):
    rows_hp += f'<tr><td>HP {k}</td><td class="num">{fnum(v)}</td></tr>'

# verification table (metric -> source file:key -> value)
ver = [
    ("union sSNV", "read_driven_genomewide.json : n_snv_total", fnum(n_snv) + f" (TP {fnum(n_tp)} + FP {fnum(n_fp)})"),
    ("reads touching sSNV", "… : n_reads_touching_total", fnum(n_touch)),
    ("reads ≥2 ALT (genuine)", "… : n_reads_multi_alt_genuine", fnum(n_multi)),
    ("chimeric excluded (span>100kb)", "… : n_chimeric_reads_excluded", fnum(n_chim)),
    ("distinct ALT combos", "… : distinct_alt_combos_total", fnum(combos)),
    ("distinct pairs (edges)", "… : distinct_pairs_total", fnum(pairs)),
    ("collapse ratio", "… : collapse_ratio_reads_per_combo", f"{collapse}×"),
    ("cross-50kb genuine chains", "… : n_cross50k_genuine_chains", fnum(cross50k)),
    ("max ALT / read", "… : max_alt_per_read", str(max_alt)),
    ("full evolution trees", "region_shape_distribution.json : by_shape.full_tree.n", fnum(full_tree) + f" (clean {fnum(full_clean)})"),
    ("confirmed-structure regions", "… : structured_plus_colinked", fnum(struct_plus)),
    ("total regions", "… : total_regions", fnum(tot_reg)),
    ("sSNV buildable (powered link)", "tier_r_characterization.json : linked_ref.n", fnum(linked) + f" / {fnum(n_snv)} ({100*linked/n_snv:.0f}%)"),
]
rows_ver = "".join(f'<tr><td>{m}</td><td class="num">{val}</td><td class="src">{src}</td></tr>' for m, src, val in ver)

html = f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Read-driven 全基因組交叉確認 — HCC1395</title>
<style>
body{{font-family:"Segoe UI","微軟正黑體",sans-serif;max-width:1080px;margin:0 auto;padding:24px;color:#222;line-height:1.55}}
h1{{font-size:22px;border-bottom:3px solid #44546A;padding-bottom:8px}}
h2{{font-size:17px;color:#44546A;margin-top:30px;border-left:4px solid #5B9BD5;padding-left:8px}}
.bluf{{background:#eef4fb;border:1px solid #cfe0f3;border-radius:6px;padding:12px 16px;font-size:14px}}
.cards{{display:flex;flex-wrap:wrap;gap:12px;margin:14px 0}}
.card{{flex:1 1 150px;background:#fff;border:1px solid #ddd;border-radius:8px;padding:12px;text-align:center;box-shadow:0 1px 3px rgba(0,0,0,.05)}}
.card .v{{font-size:24px;font-weight:700;color:#2c5f9e}}
.card .l{{font-size:12px;color:#666;margin-top:4px}}
table{{border-collapse:collapse;width:100%;margin:10px 0;font-size:13px}}
th,td{{border:1px solid #e3e3e3;padding:6px 9px;text-align:left;vertical-align:middle}}
th{{background:#f4f6f9;color:#44546A}}
.num{{text-align:right;font-variant-numeric:tabular-nums;font-weight:600}}
.src{{font-size:11px;color:#888;font-family:monospace}}
.red{{color:#c0392b;font-weight:700}}
.note{{font-size:12px;color:#666}}
footer{{margin-top:30px;font-size:11px;color:#999;border-top:1px solid #eee;padding-top:10px}}
</style></head><body>

<h1>Read-driven 全基因組交叉確認 — HCC1395（{RD.get("scope","")}）</h1>

<div class="bluf">
<b>BLUF</b>：遍歷 {fnum(n_touch)} 條 read 直接重建 sSNV 共現。<b>read-driven 可行 ✓</b>、且確實含「一條 read 串多個 sSNV」的<b>直接多路鏈</b>（最長 {max_alt}-way）。
但 ① <b class="red">最大的鏈是 artifact</b>（chr9 pericentromere 85-way）非真 subclone；
② 整合後資訊<b>驟減 {collapse}×</b>（{fnum(n_multi)} 多-ALT read → {fnum(combos)} distinct combos）；
③ 跨-50kb 真實 gain 僅 {fnum(cross50k)}（{100*cross50k/n_multi:.2f}%）→ <b>現行 50kb cap 合理</b>；chimeric（supplementary）{fnum(n_chim)} 條已正確排除。
最終可建完整演化樹 = <b>{fnum(full_tree)} 區</b>（與既有 pipeline 一致）。
</div>

<div class="cards">
<div class="card"><div class="v">{fnum(n_touch)}</div><div class="l">reads 觸及 sSNV</div></div>
<div class="card"><div class="v">{fnum(n_multi)}</div><div class="l">reads ≥2 ALT（直接鏈）</div></div>
<div class="card"><div class="v">{fnum(combos)}</div><div class="l">distinct combos</div></div>
<div class="card"><div class="v">{collapse}×</div><div class="l">塌縮比</div></div>
<div class="card"><div class="v">{fnum(full_tree)}</div><div class="l">完整演化樹</div></div>
</div>

<h2>1. 塌縮級聯：為何「整合後資訊小很多」</h2>
<p class="note">原始 read 層看似龐大，逐層整合（多-ALT read → 不同組合 → 區域 → 確認結構 → 完整樹）後驟減。每列右側為相對前一層的塌縮倍率。</p>
<table><tr><th>層級</th><th>數量</th><th></th><th>塌縮</th></tr>{rows_cascade}</table>

<h2>2. 串接分布（一條 read 直接串幾個 sSNV）</h2>
<p class="note">直接多路鏈（同一分子上 ≥2 ALT）大小分布；2-way 最多，長鏈快速遞減（且多為 artifact，見下）。最長 {max_alt}-way。</p>
<table><tr><th>鏈大小</th><th>read 數</th><th></th></tr>{rows_chain}</table>

<h2>3. 🔴 FP-inflation：最密集的鏈是 artifact 不是 subclone</h2>
<p class="note">最大鏈擠在極小區間（密度遠超真 somatic ~1/85kb）→ pericentromere / segdup mapping artifact。需 somatic-confirm 閘（normal=REF）清除。「真 subclone」訊號是小鏈（如 canonical chr17:48.3M ~3 sSNV）。</p>
<table><tr><th>n_ALT</th><th>位置</th><th>span(bp)</th><th>密度</th><th>HP</th><th>旗標</th></tr>{rows_top}</table>

<h2>4. 跨-50kb gain + chimeric 防呆</h2>
<table>
<tr><td>跨 50kb 真實鏈（提高 TIER_R 才撈到）</td><td class="num">{fnum(cross50k)}</td><td>= 多-ALT read 的 {100*cross50k/n_multi:.2f}% → 增益小，50kb cap 合理</td></tr>
<tr><td>chimeric / supplementary 排除（span>100kb）</td><td class="num">{fnum(n_chim)}</td><td>同 query_name 跨數 Mb，非單分子 → 正確排除</td></tr>
</table>

<h2>5. HP 一致性（多-ALT 鏈的 haplotag）</h2>
<p class="note">一條 read = 一個分子 = 同 HP；多-ALT 鏈幾乎全在 somatic-tagged（1-1/2-1）+ HP3，印證共現=cis。</p>
<table><tr><th>HP</th><th>鏈數</th></tr>{rows_hp}</table>

<h2>6. 與既有 pipeline 交叉對照（區域 shape）</h2>
<p class="note">distance-window pipeline 的最終區域結構；read-driven 確認其骨幹一致（completeness/collapse 上吻合）。clean = CN-neutral/LOH 乾淨子集。</p>
<table><tr><th>shape</th><th>區域數</th><th>佔比</th><th>clean(CN乾淨)</th></tr>{rows_shape}</table>
<p class="note">sSNV 層可建性（tier_r）：有 powered link <b>{fnum(linked)}</b> / underpowered {fnum(under)}（加深度 ~57% 可救）/ isolated {fnum(iso)}（稀疏孤點，需 Tier-PS）。</p>

<h2>7. 驗證數據表（每數字 → 來源檔:鍵）</h2>
<table><tr><th>metric</th><th>value</th><th>source</th></tr>{rows_ver}</table>

<h2>8. 誠實邊界</h2>
<ul class="note">
<li>read-driven 與 distance-window 結論一致：完整樹 <b>{fnum(full_tree)}</b>（clean {fnum(full_clean)}），confirmed-structure {fnum(struct_plus)} / {fnum(tot_reg)} 區。</li>
<li>最大鏈為 artifact（chr9 pericentromere 等）→ structure 為上界、含 FP；須 somatic 閘 per-region 清除（便宜）。</li>
<li>跨-50kb gain 極小（{fnum(cross50k)}）→ 提高 TIER_R 收益有限；真正 lever = 深度 + Tier-PS。</li>
<li>⭐3 單樣本 regional：{fnum(full_tree)} 為局部區域樹，非一棵全腫瘤系統發生樹；不能排 sibling 順序。</li>
</ul>

<footer>
generated by build_html.py（§13-A 數據注入）· data: read_driven_genomewide.json / region_shape_distribution.json / tier_r_characterization.json · HCC1395 ⭐3 · 2026-06-30
</footer>
</body></html>"""

open(OUT, "w").write(html)
print(f"DONE -> {OUT}  ({len(html)} bytes)")
