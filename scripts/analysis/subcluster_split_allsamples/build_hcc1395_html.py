#!/usr/bin/env python3
"""HCC1395 per-SNV-locus 全面整理 standalone HTML — 數字全注入自 hcc1395_full_accounting.json (§13-A)。
C++ 語義與 caveats 來自對抗驗證 workflow。缺 JSON 直接 refuse。
用法: build_hcc1395_html.py <ACCOUNTING_JSON> <COMMIT> <OUT_HTML>
所有動態片段先預算為變數，再組裝（避免 f-string 巢狀）。
"""
import sys, json, os, html as H

JP = sys.argv[1]; COMMIT = sys.argv[2] if len(sys.argv) > 2 else "2d5692f"
OUT = sys.argv[3] if len(sys.argv) > 3 else "HCC1395_per_locus.standalone.html"
if not os.path.exists(JP):
    sys.exit(f"REFUSE: 缺真值 JSON {JP}（§13-A 不可捏造）")
d = json.load(open(JP))
if "schema_contract" not in d:
    sys.exit("REFUSE: accounting JSON 缺 schema_contract；不可猜 current/legacy/region 語意")
schema = d["schema_contract"]
required_schema = {
    "verification_selection_field",
    "verification_schema_status",
    "verification_categories",
    "loh_selection_field",
    "region_status",
}
missing_schema = sorted(required_schema - set(schema))
if missing_schema:
    sys.exit(f"REFUSE: accounting JSON schema_contract 缺欄 {missing_schema}")
N = d["N_significance"]

def pc(n): return f"{round(100*n/N,1)}%"
def e(s): return H.escape(str(s))
def cm(n): return f"{n:,}"

split = d["dim1_split"]; venn = d["dim5_venn"]; tl = venn["three_lens"]
perm = d["dim2_permanova"]; ex = d["dim7_extra_axes"]; ctx = d["dim4_context"]; ct = d["dim6_crosstab"]
mode = d["mode"]; sigrec = d["dim3_classification"]["Significant_reconciliation"]
vc = d["dim3_classification"]["VerificationClass_Current"]
sg = d["dim3_classification"]["StrengthGrade_AtoE"]
cov = ctx["Coverage_Category"]; loh = ctx["LOH_Subtype_LegacyVC"]

vc_categories = list(schema["verification_categories"])
noise_classes = ["Noise_Uncorrelated", "Noise_Uniform", "Noise_Chaotic"]
def grp(keys): return sum(vc.get(k, 0) for k in keys)

VC_DEF = {
 "Strong_Bidirectional": "legacy label-first 與 cluster-first evidence 同時成立",
 "ClusterFirstOnly": "cluster-first evidence 成立、label-first 不成立；不是 cellular clone",
 "LOH-Structure": "潛在 LOH 位點 + HP-AUC 結構或顯著 Δβ",
 "MultiGroupNoLabel": "within-HP 乾淨多群子結構但無標籤對齊",
 "LabelShift": "Δβ mean-shift evidence",
 "PermanovaLocation": "label PERMANOVA 顯著且 PERMDISP 無離散度警告（乾淨位置位移）",
 "StructureNoLabel": "HP-AUC≥0.7 結構但無標籤/LOH/within-HP/PERMANOVA 匹配",
 "DispersionStructure": "PERMANOVA 顯著但離散度警告 — 有結構但歸因模糊（不算 Significant）",
 "Noise_Uncorrelated": "中間帶距離/epipoly，無乾淨群（預設 Noise）",
 "Noise_Uniform": "甲基均勻不可分（pairwise<0.15 或 epipoly<0.2）",
 "Noise_Chaotic": "高熵混亂（epipoly>0.5 或 pairwise>0.35）",
 "UnknownCurrentClass": "schema 已標記但 class enum 未知；保留且不得折疊為 Noise",
 "Strong": "UNVERSIONED v1 raw class；只可出現在明確 historical mode",
 "Subclone": "UNVERSIONED v1 raw class；不是 cellular clone",
 "Weak": "UNVERSIONED v1 raw class",
 "Noise": "UNVERSIONED v1 raw class",
}
GRADE_DEF = "A≥0.50 · B≥0.35 · C≥0.22 · D≥0.12 · E<0.12（StrengthScore live cutoffs）"

def bar(n, color="#8a5a2b", tot=None):
    tot = tot or N
    w = round(100*n/tot, 1)
    return f"<div class='barwrap'><div class='bar' style='width:{w}%;background:{color}'></div><span class='barlbl'>{cm(n)} ({w}%)</span></div>"

def table(headers, rows_html, cls=""):
    th = "".join(f"<th>{e(h)}</th>" for h in headers)
    return f"<table class='{cls}'><thead><tr>{th}</tr></thead><tbody>{rows_html}</tbody></table>"

# ---------- 預算所有動態片段 ----------
# Schema-explicit current taxonomy; Significant is read independently, not reconstructed here.
vc_rows = ""
for k in vc_categories:
    n = vc.get(k, 0)
    vc_rows += f"<tr><td>{e(k)}</td><td class='num'>{cm(n)}</td><td class='num'>{pc(n)}</td><td class='def'>{e(VC_DEF.get(k,''))}</td></tr>"

# split table
split_tbl = table(["狀態", "位點數", "佔比"],
  f"<tr><td>✅ cansplit（切得出 k≥2）</td><td>{bar(split['cansplit_total'],'#1f7a4d')}</td><td class='num'>{split['cansplit_pct']}%</td></tr>"
  f"<tr><td>✋ cant（reads≥6 但切不出平衡 k≥2）</td><td>{bar(split['cant_total'],'#b3541e')}</td><td class='num'>{split['cant_pct']}%</td></tr>"
  f"<tr><td>◌ insuff（tumor reads &lt;6，無法分群）</td><td>{bar(split['insuff'],'#999')}</td><td class='num'>{split['insuff_pct']}%</td></tr>")

# PERMANOVA axes
perm_rows = ""
for ax, v in [("LabelHP (HP 標籤軸)", perm["LabelHP"]), ("LabelAllele (REF/ALT 軸)", perm["LabelAllele"]), ("Cluster (無監督群軸)", perm["Cluster"])]:
    perm_rows += f"<tr><td>{e(ax)}</td><td class='num'>{cm(v['sig'])} ({v['sig_pct']}%)</td><td class='num'>{cm(v['clean'])} ({v['clean_pct']}%)</td><td class='num'>{cm(v['disp'])}</td></tr>"
perm_tbl = table(["軸", "顯著 sig (p&lt;.05)", "乾淨 clean(扣離散度)", "離散度 dispersion"], perm_rows)

# CramersV bands
cv_gated = perm["CramersV_exactly0_gated"]
cv_rows = ""
for k, v in d["dim2_permanova"]["CramersV_bands"].items():
    note = ""
    if k == "<0.3":
        note = f"其中 exactly 0（Cochran-gated/無關聯）= {cm(cv_gated['n'])}（非連續弱關聯，是 0 處尖峰）"
    cv_rows += f"<tr><td>{e(k)}</td><td class='num'>{cm(v)}</td><td class='def'>{e(note)}</td></tr>"
cv_tbl = table(["帶", "位點數", "說明"], cv_rows)

# StrengthGrade
sg_rows = ""
for k, v in sorted(sg.items()):
    color = "#8a5a2b" if k in ("A", "B") else ("#b89b6a" if k == "C" else "#cfc4ad")
    sg_rows += f"<tr><td>{e(k)}</td><td>{bar(v, color)}</td><td class='num'>{pc(v)}</td></tr>"
sg_tbl = table(["級", "位點數", "佔比"], sg_rows)

# gating table
pg = d['dim3_classification']['PassedGating'].get('true', 0)
sf = d['dim3_classification']['SuggestFilter'].get('true', 0)
gating_tbl = table(["欄位", "true", "佔比", "定義"],
  f"<tr><td>PassedGating</td><td class='num'>{cm(pg)}</td><td class='num'>{pc(pg)}</td><td class='def'>Fisher p≤0.1 且 CramérV≥0.1（alt 或 hp-family 任一）</td></tr>"
  f"<tr><td>SuggestFilter</td><td class='num'>{cm(sf)}</td><td class='num'>{pc(sf)}</td><td class='def'>label_delta&gt;0.3（legacy F1 啟發式）</td></tr>")

# 3-lens table
lens_tbl = table(["判準", "位點數", "佔比", "本質"],
  f"<tr><td>① cansplit（無監督 k≥2）</td><td class='num'>{cm(tl['cansplit'])}</td><td class='num'>{pc(tl['cansplit'])}</td><td class='def'>pooled reads 找到分群（可為 cis-ASM/germline-HP/CN）</td></tr>"
  f"<tr><td>② C clean PERMANOVA</td><td class='num'>{cm(tl['C_clean'])}</td><td class='num'>{pc(tl['C_clean'])}</td><td class='def'>label 軸乾淨位置位移</td></tr>"
  f"<tr><td>③ Significant verdict</td><td class='num'>{cm(tl['Significant'])}</td><td class='num'>{pc(tl['Significant'])}</td><td class='def'>直接讀 canonical Significant boolean</td></tr>")

# HP_AUC bands
auc_rows = "".join(f"<tr><td>{e(k)}</td><td class='num'>{cm(v)}</td><td class='num'>{pc(v)}</td></tr>" for k, v in ex['HP_AUC_Tumor_bands'].items())
auc_tbl = table(["帶", "位點數", "佔比"], auc_rows)
auc70 = ex['HP_AUC_Tumor_bands'].get('≥0.7', 0)

# TumorOnly K
tok = ex['TumorOnlyClusterK']
tok_rows = "".join(f"<tr><td>K={e(k)}</td><td class='num'>{cm(v)}</td></tr>" for k, v in sorted(tok.items(), key=lambda x: int(x[0]) if str(x[0]).isdigit() else 99))
tok_tbl = table(["K", "位點數"], tok_rows)
tok2_pct = round(100*tok.get('2', 0)/N)

# StrengthScore components
ssc = ex['StrengthScore_component_means']
ssc_tbl = table(["分量"] + list(ssc.keys()), "<tr><td>均值</td>" + "".join(f"<td class='num'>{v}</td>" for v in ssc.values()) + "</tr>")

# extra axes (Fisher / NHP)
extra_tbl = table(["軸", "非零位點", "佔比"],
  f"<tr><td>Fisher per-CpG ASM（≥1 顯著 CpG）</td><td class='num'>{cm(ex['Fisher_perCpG_ASM_nonzero']['n'])}</td><td class='num'>{ex['Fisher_perCpG_ASM_nonzero']['pct']}%</td></tr>"
  f"<tr><td>NHP_Somatic 1-1</td><td class='num'>{cm(ex['NHP_Somatic11_nonzero']['n'])}</td><td class='num'>{ex['NHP_Somatic11_nonzero']['pct']}%</td></tr>"
  f"<tr><td>NHP_Somatic 2-1</td><td class='num'>{cm(ex['NHP_Somatic21_nonzero']['n'])}</td><td class='num'>{ex['NHP_Somatic21_nonzero']['pct']}%</td></tr>"
  f"<tr><td>NHP_Somatic 3-3</td><td class='num'>{cm(ex['NHP_Somatic33_nonzero']['n'])}</td><td class='num'>{ex['NHP_Somatic33_nonzero']['pct']}%</td></tr>")

# CNV / LOH
cov_tbl = table(["類", "位點數", "佔比"], "".join(f"<tr><td>{e(k)}</td><td class='num'>{cm(v)}</td><td class='num'>{pc(v)}</td></tr>" for k, v in cov.items()))
loh_tbl = table(["類", "位點數", "佔比"], "".join(f"<tr><td>{e(k)}</td><td class='num'>{cm(v)}</td><td class='num'>{pc(v)}</td></tr>" for k, v in loh.items()))
ploh = ctx['Potential_LOH'].get('true', 0)
within_tbl = table(["欄位", "true", "佔比"],
  f"<tr><td>WithinHP_CleanMultigroup</td><td class='num'>{cm(ctx['WithinHP_CleanMultigroup'].get('true',0))}</td><td class='num'>{pc(ctx['WithinHP_CleanMultigroup'].get('true',0))}</td></tr>"
  f"<tr><td>WithinHP_LevelBimodal</td><td class='num'>{cm(ctx['WithinHP_LevelBimodal'].get('true',0))}</td><td class='num'>{pc(ctx['WithinHP_LevelBimodal'].get('true',0))}</td></tr>"
  f"<tr><td>WithinHP_SubcloneSig（a-priori germline-vs-carrier）</td><td class='num'>{cm(ctx['WithinHP_SubcloneSig'].get('true',0))}</td><td class='num'>{pc(ctx['WithinHP_SubcloneSig'].get('true',0))}</td></tr>")
legacy_subclone = d['dim3_classification']['VerificationClass_Legacy_4state'].get('Subclone', 0)

# split × VC crosstab
order = vc_categories
sx = ct["split_x_VerificationClass"]
cross_rows = ""
for state, label in [("cansplit", "✅ cansplit(切得出)"), ("cant", "✋ cant(切不出)"), ("insuff", "◌ insuff(reads<6)")]:
    c = sx.get(state, {}); tot = sum(c.values())
    cells = "".join(f"<td class='num'>{cm(c.get(k,0))}</td>" for k in order)
    cross_rows += f"<tr><td>{e(label)}</td><td class='num'><b>{cm(tot)}</b></td>{cells}</tr>"
cross_tbl = table(["切割狀態", "小計"] + order, cross_rows, "cross")

# empty audit
empty_rows = ""
for k, v in d['dim8_tumor_only_empty_audit']['columns_not_computable'].items():
    st = "全部 = 不可計算" if v == N else "部分"
    empty_rows += f"<tr><td>{e(k)}</td><td class='num'>{cm(v)}</td><td class='def'>{e(st)}</td></tr>"
empty_tbl = table(["欄位", "空/0/−1 位點數", "狀態"], empty_rows)

# KPI computed values
sig_flag = sigrec['significant_flag']; sigE = sigrec['of_which_StrengthGrade_E']
disp_struct = vc.get('DispersionStructure', 0); noise_tot = grp(noise_classes)
cluster_supported = vc.get('Strong_Bidirectional', 0) + vc.get('ClusterFirstOnly', 0)
ab = sg.get('A', 0) + sg.get('B', 0)
region_status = schema["region_status"]
region_note = (
    f"RegionStratifier status={region_status.get('status')} / reason={region_status.get('reason')}"
    if region_status.get("status") != "UNVERSIONED_REGION_SCHEMA"
    else "Historical unversioned input：RegionStratification schema unavailable；deprecated Subclone_ID 不作語意解讀"
)

doc = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>HCC1395 per-SNV-locus 切割與分類整理</title>
<style>
:root{{--ink:#1c1c1c;--mut:#6b6b6b;--line:#e4dfd4;--bg:#faf8f3;--acc:#8a5a2b;--ok:#1f7a4d;--warn:#b3541e;--card:#fff}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Helvetica Neue","PingFang TC","Microsoft JhengHei",sans-serif;color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:1240px;margin:0 auto;padding:28px 26px 90px}}
h1{{font-size:23px;margin:0 0 4px}}.sub{{color:var(--mut);font-size:12.5px;margin:2px 0}}
.modestamp{{background:#fdf0e4;border:1px solid var(--warn);border-left:5px solid var(--warn);border-radius:8px;padding:12px 16px;margin:14px 0;font-size:13px}}
.modestamp b{{color:var(--warn)}}
nav{{position:sticky;top:0;background:rgba(250,248,243,.96);backdrop-filter:blur(6px);border-bottom:1px solid var(--line);padding:9px 0;margin:16px 0 22px;font-size:12.5px;z-index:20}}
nav a{{color:var(--acc);text-decoration:none;margin-right:13px;white-space:nowrap}}nav a:hover{{text-decoration:underline}}
h2{{font-size:18px;border-left:4px solid var(--acc);padding-left:10px;margin:32px 0 12px}}
h3{{font-size:14px;margin:18px 0 8px;color:#4a4036}}
.head{{background:linear-gradient(180deg,#fff,#fcfaf5);border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:8px 0}}
.head .big{{font-size:15px;font-weight:600}}.head b{{color:var(--acc)}}
table{{width:100%;border-collapse:collapse;font-size:12.5px;background:var(--card);border:1px solid var(--line);border-radius:8px;overflow:hidden;margin:8px 0 6px}}
th,td{{padding:6px 9px;border-bottom:1px solid var(--line);text-align:left}}
td.num,th.num{{text-align:right;font-variant-numeric:tabular-nums}}
thead th{{background:#f1ece2;color:#4a4036;font-size:11.5px}}
tbody tr:hover{{background:#fcf8f0}}
tr.grphead td{{background:#f6f1e8;font-weight:700;font-size:12px;color:#5a4a36}}
td.def{{color:var(--mut);font-size:11px}}
.barwrap{{position:relative;background:#efe9dd;border-radius:5px;height:20px;min-width:210px}}
.bar{{height:100%;border-radius:5px}}.barlbl{{position:absolute;left:8px;top:1px;font-size:11.5px;color:#222;line-height:18px}}
.card{{background:var(--card);border:1px solid var(--line);border-left:4px solid var(--acc);border-radius:8px;padding:11px 15px;margin:10px 0;font-size:12.5px}}
.card.warn{{border-left-color:var(--warn);background:#fdf6ef}}.card.ok{{border-left-color:var(--ok)}}
.card b{{color:var(--acc)}}.card.warn b{{color:var(--warn)}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:14px}}
.note{{font-size:11.5px;color:var(--mut);margin:4px 0}}
.pill{{display:inline-block;background:#efe7d8;color:#6a4a22;border-radius:20px;padding:1px 9px;font-size:11px;margin-left:6px}}
.kpi{{display:flex;gap:10px;flex-wrap:wrap;margin:10px 0}}
.kpi .k{{flex:1;min-width:150px;background:var(--card);border:1px solid var(--line);border-radius:9px;padding:10px 13px}}
.kpi .v{{font-size:20px;font-weight:700;color:var(--acc)}}.kpi .l{{font-size:11px;color:var(--mut)}}
ul{{margin:6px 0;padding-left:19px}}li{{margin:3px 0;font-size:12.5px}}
code{{background:#f1ece2;padding:1px 5px;border-radius:4px;font-size:11px}}
footer{{margin-top:38px;padding-top:14px;border-top:1px solid var(--line);font-size:11px;color:var(--mut)}}
@media(max-width:780px){{.grid{{grid-template-columns:1fr}}}}
</style></head><body><div class="wrap">

<h1>HCC1395 — 以 SNV 位點為單位的切割 × 統計 × 方法學分類整理<span class="pill">N={cm(N)} loci</span></h1>
<div class="sub">canonical paired_full longphase-s tagged BAM · 全基因組 TP somatic SNV · binary <code>feat/summary-nreadsvalid {COMMIT}</code></div>
<div class="sub">ISM <code>w=5000 BERNOULLI C_min=3 SKIP</code> · 數字注入自 <code>hcc1395_full_accounting.json</code>（§13-A 反捏造）· 獨立重算 6/6 ALL_MATCH</div>

<div class="modestamp">
<b>⚠ 模式：TUMOR-ONLY run（Normal BAM: None）。</b> NNormalReads = 0 於全部 {cm(mode['NNormalReads_zero'])} 位點。
所有 <b>normal-anchored</b> 指標（SomaticResidualDbeta / GermlineAsmDbeta / SubcloneDbeta / SampleASM normal arm）為「<b>不可計算</b>」而非「陰性發現」——見 §8。{e(region_note)}。
</div>

<nav>
<a href="#tldr">總覽</a><a href="#split">①切割</a><a href="#perm">②統計</a><a href="#verdict">③分類</a>
<a href="#lens">④三透鏡</a><a href="#extra">⑤補軸</a><a href="#ctx">⑥CNV/LOH</a><a href="#cross">⑦交叉表</a><a href="#empty">⑧不可計算</a><a href="#caveat">口徑</a>
</nav>

<section id="tldr"><div class="head">
<div class="big">全基因組 {cm(N)} 個 TP somatic SNV — 三個獨立判準把同一批位點切得很不一樣：</div>
<div class="kpi">
<div class="k"><div class="v">{cm(split['cansplit_total'])}</div><div class="l">① 無監督切得出 cansplit ({split['cansplit_pct']}%)</div></div>
<div class="k"><div class="v">{cm(venn['C_clean'])}</div><div class="l">② 乾淨 PERMANOVA C ({venn['C_pct']}%)</div></div>
<div class="k"><div class="v">{cm(sig_flag)}</div><div class="l">③ verdict Significant ({pc(sig_flag)})</div></div>
<div class="k"><div class="v">{venn['jaccard']}</div><div class="l">Jaccard(①,②) 偏低</div></div>
</div>
<p style="margin:8px 0 0;font-size:12.5px">三者<b>部分重疊、非巢狀漏斗</b>（§4）。① cansplit 只代表「pooled reads 找到 k≥2 分群」，<b>可能是 cis-ASM / germline-HP / CN，不等於 tumor subclone</b>。</p>
</div></section>

<section id="split"><h2>① 切割維度（無監督 UPGMA + silhouette，tumor reads）</h2>
<p class="note">每位點對 tumor reads 算 read×read 距離矩陣 → UPGMA → silhouette 選 best_k。best_k≥2 = cansplit。</p>
{split_tbl}
<div class="card"><b>切得出 ≈ 1/5。</b>「切不出來」{split['cant_pct']}% 多數非失敗 = reads 相對同質。<b>切得出 ≠ subclone</b>（品質分帶與對齊軸見全樣本 master 報告）。</div>
</section>

<section id="perm"><h2>② 統計維度（PERMANOVA / GlobalTest，read×CpG 距離）</h2>
{perm_tbl}
<div class="card warn"><b>「顯著」多數是離散度而非位置位移。</b> LabelAllele sig {perm['LabelAllele']['sig_pct']}% 但乾淨僅 {perm['LabelAllele']['clean_pct']}%（{cm(perm['LabelAllele']['disp'])} 個 dispersion warning）；LabelHP sig {perm['LabelHP']['sig_pct']}% → 乾淨 {perm['LabelHP']['clean_pct']}%。PERMDISP 扣離散度後才是可信位置差異。<b>C = 任一 label 軸乾淨 = {cm(perm['C_clean_anyaxis']['n'])} ({perm['C_clean_anyaxis']['pct']}%)</b>。</div>
<h3>CramérV 關聯強度（cluster×label）</h3>
{cv_tbl}
<p class="note">⚠ &lt;0.3 帶的 {cv_gated['pct']}% 是 exactly 0（Cochran-gated），非平滑梯度。</p>
</section>

<section id="verdict"><h2>③ 方法學分類維度（C++ reclassify-v2 verdict）</h2>
<h3>VerificationClass（field={e(schema['verification_selection_field'])}；status={e(schema['verification_schema_status'])}）</h3>
{table(["VerificationClass","位點數","佔比","定義（源碼 RegionProcessor.cpp）"], vc_rows)}
<div class="card ok"><b>Significant = {cm(sig_flag)}（{pc(sig_flag)}），直接讀 canonical <code>Significant</code>。</b> 本頁不再由 class 名稱重建 truth；DispersionStructure={cm(disp_struct)}、Noise={cm(noise_tot)} 只作 current taxonomy 描述。</div>
<div class="card warn"><b>「Significant {pc(sig_flag)}」≠「{pc(sig_flag)} 有真結構」。</b> 其中 <b>{cm(sigE)}（{sigrec['E_pct_of_sig']}%）是最弱的 E 級</b>。Significant 是寬鬆位置-結構 OR-gate，非強度 verdict。</div>
<h3>StrengthGrade（強度分級）<span class="note"> — {GRADE_DEF}</span></h3>
{sg_tbl}
<p class="note">A+B 僅 {cm(ab)}（{pc(ab)}）。E 級 {pc(sg.get('E',0))} 偏高部分因 tumor-only：StrengthScore 5 分量中 StrengthSomatic/StrengthGermline 為 normal-anchored，此 run 恆 0（§5）。</p>
<h3>PassedGating / SuggestFilter</h3>
{gating_tbl}
</section>

<section id="lens"><h2>④ 三透鏡非巢狀（同批位點，不同判準）</h2>
{lens_tbl}
<div class="card warn"><b>三者部分重疊，不是漏斗。</b> ①∩③ = {cm(tl['cansplit_and_Sig'])}；只③不① = {cm(tl['Sig_only'])}；只①不③ = {cm(tl['cansplit_only'])}；①∩② = {cm(tl['cansplit_and_Cclean'])}（Jaccard(①,②)={venn['jaccard']}）。看到三個百分比<b>勿讀成子集鏈</b>。</div>
</section>

<section id="extra"><h2>⑤ 補充軸（read-level 訊號）</h2>
<div class="grid">
<div><h3>HP_AUC_Tumor（HP 可分性）</h3>{auc_tbl}<p class="note">≥0.7（HP 強可分）= {cm(auc70)}（{pc(auc70)}）</p></div>
<div><h3>TumorOnly 分群 K 分佈</h3>{tok_tbl}<p class="note">⚠ K=2 占 {tok2_pct}% — 近 chance 二分；TumorIntrinsic/TumorOnly PERMANOVA 為 double-dip（memory 已證 tumor-only 軸 NEGATIVE），勿當陽性。</p></div>
</div>
<h3>StrengthScore 5 分量均值（解釋 E 級偏高）</h3>
{ssc_tbl}
<p class="note">StrengthSomatic / StrengthGermline = 0（normal-anchored，tumor-only 恆 0）→ 5 分量中 2 個結構性為 0，壓低總分。</p>
<h3>per-CpG / 體細胞 haplotype 拓樸</h3>
{extra_tbl}
</section>

<section id="ctx"><h2>⑥ CNV / LOH / within-HP 脈絡</h2>
<div class="grid">
<div><h3>Coverage_Category（CNV 代理）</h3>{cov_tbl}</div>
<div><h3>{e(schema['loh_selection_field'])}（legacy-derived provenance）</h3>{loh_tbl}<p class="note">Potential_LOH = {cm(ploh)}（{pc(ploh)}）</p></div>
</div>
<h3>within-HP 子結構（subclone-direction，無監督）</h3>
{within_tbl}
<div class="card warn"><b>legacy「Subclone」類（{cm(legacy_subclone)}）名過其實</b> = within-HP multi-group candidate（無監督），<b>非</b> normal-anchored subclone 軸；確認 ≥2 within-HP 群 ≠ subclone（可為 cis-ASM，~85% 對齊 germline-carrier）。</div>
</section>

<section id="cross"><h2>⑦ 交叉表：切割狀態 × VerificationClass</h2>
{cross_tbl}
<div class="card"><b>切割與 verification evidence 是不同判準。</b> cluster-supported current classes 共 {cm(cluster_supported)}；cansplit 仍只代表無監督 k≥2，不能據此聲稱 cellular subclone，兩者互相不蘊含。</div>
</section>

<section id="empty"><h2>⑧ Tumor-only 模式：不可計算欄位（非陰性）</h2>
<p class="note">以下 normal-anchored 欄位在此 tumor-only run 全為空/0/−1，是「未計算」非「無發現」。解讀 subclone 時務必排除。</p>
{empty_tbl}
</section>

<section id="caveat"><h2>口徑（給 PI / 審稿）</h2>
<ul>
<li><b>單樣本 tumor-only</b>：normal-anchored subclone 軸不可計算（§8）；tier 上限 ⭐3 proof-of-concept；無跨樣本/無 paired 確認。</li>
<li><b>cansplit ≠ subclone</b>：無監督 k≥2 可為 cis-ASM / germline-HP / CN；確認 ≥2 within-HP 群亦 ≠ subclone。</li>
<li><b>Significant 是寬鬆 OR-gate</b>（{pc(sig_flag)}），57% 是 E 級；勿讀成「有真結構」。</li>
<li><b>「顯著」多為離散度</b>：PERMDISP 扣除後乾淨位置位移大幅縮減（LabelAllele {perm['LabelAllele']['sig_pct']}%→{perm['LabelAllele']['clean_pct']}%）。</li>
<li><b>TumorIntrinsic / TumorOnly PERMANOVA = double-dip</b>（K=2 近 chance 占多數），memory 已證 tumor-only 軸 NEGATIVE。</li>
<li><b>CramérV &lt;0.3 帶</b> {cv_gated['pct']}% 是 gated-to-0 尖峰，非弱關聯梯度。</li>
<li><b>denominator</b> = {cm(N)} matched（of 29,754 records；14 未匹配）。</li>
</ul>
</section>

<footer>
provenance：binary <code>inter_sub_mod @ feat/summary-nreadsvalid {COMMIT}</code> · verification field <code>{e(schema['verification_selection_field'])}</code> / schema <code>{e(schema['verification_schema_status'])}</code> · LOH field <code>{e(schema['loh_selection_field'])}</code> · region status <code>{e(region_status.get('status'))}</code> · 數字注入自 <code>hcc1395_full_accounting.json</code> · tumor-only 單樣本 tier⭐3<br>
報告目錄 <code>InterSubMod/docs/experiments/in_progress/2026/06/20260620_allsample_subcluster_split/</code>
</footer>
</div></body></html>"""

open(OUT, "w").write(doc)
print(f"[HCC1395 HTML] {OUT} ({len(doc)} chars)")
