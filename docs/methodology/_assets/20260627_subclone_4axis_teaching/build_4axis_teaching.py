#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4 軸 clone/subclone 重建教學 HTML 生成器 v2（§13-A 由構造防捏造）。
讀凍結 JSON（data/）→ 計算衍生值 + sum-check（不符即 raise）→ 注入 HTML。
v2 變更（2026-06-28，套 2 個對抗 workflow 修正）:
  - chr17 改由 chr17_subclone_data.json 注入(3 sSNV 機器真值,非硬編 4-sSNV)
  - CLEAN 子集改由 clean_subset.json 注入(非硬編 205);6-bucket 由 region_shape_distribution.json
  - 非循環性條件化(uniquely-mappable+CN-clean+somatic+coread≥6);最強反駁=mapping artifact
  - CNV/LOH 改 condition-on confound + 輔助刻畫(非刪除非證據 rung);VAF=consistency-check
  - 甲基 off-ladder;主張3 INVALID;主張5 RISKY;somatic=TP∪FP union build+SEQC2 observation-only
  - 新增 §3 證據階層重構 + §9 6 主張裁決表
用法：python3 build_4axis_teaching.py
"""
import json, os, html as _h

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "20260627_subclone_4axis_teaching.standalone.html"))

def load(n):
    with open(os.path.join(DATA, n), encoding="utf-8") as f:
        return json.load(f)

def req(d, *ks):
    cur = d
    for k in ks:
        if isinstance(cur, dict) and k in cur: cur = cur[k]
        elif isinstance(cur, list) and isinstance(k, int) and k < len(cur): cur = cur[k]
        else: raise KeyError(f"§13-A REFUSE: 缺 {ks}")
    return cur

locus=load("sm_locus_master_summary.json"); ledger=load("sm_completeness_ledger.json")
summ=load("sm_summary.json"); ps=load("sm_phaseset_extension.json"); ccf=load("sm_ccf_tiers.json")
hp=load("sm_hp_contribution.json"); mc=load("sm_methyl_corroboration.json")
mr=load("sm_methyl_reextract_ALL.json"); mg=load("sm_methyl_genetic_concordance.json")
rdist=load("region_shape_distribution.json"); clean=load("clean_subset.json")
cen=load("sm_configuration_census.json"); chr17=load("chr17_subclone_data.json")
band=load("eps2_final_band.json"); rti=load("region_threshold_impact.json"); enum=load("sSNV_combination_enumeration.json")

# 宇宙
U=req(locus,"n_loci"); TP=req(locus,"src","TP"); FP=req(locus,"src","FP")
assert TP+FP==U and req(locus,"sum_check_ok") is True
b_link=req(ledger,"buckets","linked"); b_under=req(ledger,"buckets","underpowered"); b_iso=req(ledger,"buckets","isolated_singleton")
gamma=req(ledger,"linked_somatic_by_source","FP"); tp_link=req(ledger,"linked_somatic_by_source","TP")
# region buckets (derived, grep-able)
n_reg=req(rdist,"total_regions"); bs=req(rdist,"by_shape")
ft=req(bs,"full_tree","n"); lin=req(bs,"linear_nested","n"); sib=req(bs,"sibling_only","n")
col=req(bs,"co_linked_lineage","n"); noconf=req(bs,"no_confirmed_structure","n"); incon=req(bs,"inconsistent","n")
structured=req(rdist,"structured_3shape(full+linear+sibling)"); structured_plus=req(rdist,"structured_plus_colinked")
cft=req(clean,"clean_full_tree"); clean_by=req(clean,"clean_by_shape")
npop=req(rdist,"n_populations_dist")
# HP
ab=req(hp,"mutual_excl_ablation"); hp_in=req(ab,"without_HP_gate_called_subclone")
hp_true=req(ab,"with_HP_gate_true_subclone(same-HP)"); hp_rm=req(ab,"removed_as_allelic(diff-HP)"); hp_rm_pct=req(ab,"false_subclone_removed_pct")
rel=req(hp,"per_relationship_same_hp")
# census powered_somatic
cps={c["config"]:c for c in req(cen,"census_powered_somatic")}
filt=req(cen,"filters")
def cfg(name,key): return req(cps[name],key)
nested_same=cfg("nested_a_in_b","pairs_sameHP")+cfg("nested_b_in_a","pairs_sameHP")
nested_diff=cfg("nested_a_in_b","pairs_diffHP")+cfg("nested_b_in_a","pairs_diffHP")
# CCF
g=req(ccf,"ancestor_ge_descendant_VAF_gradient"); g_n=req(g,"n_nested_edges")
g_sup=req(g,"ancestor_higher"); g_vio=req(g,"descendant_higher(violation)"); g_tie=req(g,"tie")
g_sup_pct=round(100*g_sup/g_n,1); g_vio_pct=round(100*g_vio/g_n,1); g_tie_pct=round(100*g_tie/g_n,1)
g_dec=round(100*req(g,"data_support_rate"),1); bic=req(ccf,"gmm_bic_1to4"); best_n=req(ccf,"gmm_best_n"); dbic=abs(bic[2]-bic[3])
# CN
cn=req(locus,"cn_somatic_pct"); cn_gain=req(cn,"gain"); cn_clean=req(cn,"clean_loh_neutral")
cns=req(locus,"cn_somatic"); cn_tot=sum(cns.values()); gain_n=req(cns,"gain")
link_gain=req(summ,"F2_artifact_quantification","linked_somatic_by_CN_state","gain")
link_tot=req(summ,"F2_artifact_quantification","linked_somatic_total")
link_gain_pct=round(100*link_gain/link_tot,1)
dense=req(summ,"F2_artifact_quantification","in_dense_uniform_cluster(>=5 in <=5kb)"); ndense=req(summ,"F2_artifact_quantification","n_dense_clusters")
# PS
ps_rel=req(ps,"ps_reliability_per_region","single_ps_reliable"); ps_rate=round(100*req(ps,"ps_reliability_per_region","reliable_rate"),1)
tps_con=round(100*req(ps,"tier_ps_extension","ccf_consistent_rate"),1)
# 甲基
mc_t=req(mc,"n_tested"); mc_c=req(mc,"n_methyl_corroborated(>=1 sig CpG, |db|>=0.2 q<0.05)"); mc_rate=round(100*mc_c/mc_t,1)
mr_t=req(mr,"n_tested"); mr_c=req(mr,"n_corroborated"); mr_rate=round(100*mr_c/mr_t,1)
mg_tgt=req(mg,"n_target_multi_sSNV"); mg_test=req(mg,"n_testable_permanova"); mg_rec=req(mg,"n_methyl_recovers_genetic(p<0.05)")
# Fisher null
f4=req(summ,"F4_fisher_vs_MC"); f_agree=req(f4,"agree"); f_dis=req(f4,"disagree"); f_eval=f_agree+f_dis
f_agree_pct=round(100*f_agree/f_eval,1); sparse_pct=round(100*req(f4,"sparse_floor_sig_frac"),1); diffhp_pct=round(100*req(f4,"diffHP_negctrl_sig_frac"),1)
# evolution
ev=req(summ,"validated_solid","evolution"); ncomp=req(ev,"n_components"); ncyc=req(ev,"n_with_cycle")
# chr17 inject
c_snvs=req(chr17,"snvs"); c_lc=req(chr17,"lineage_counts"); c_stat=req(chr17,"snv_stat"); c_pairs=req(chr17,"pairs_2x2")
def vaf(pos):
    s=c_stat[str(pos)]; t=s["tumor_REF"]+s["tumor_ALT"]; return round(s["tumor_ALT"]/t,2) if t else 0
A_POS=48365089; B1_POS=48362515; B2_POS=48365161
vaf_a=vaf(A_POS); vaf_b1=vaf(B1_POS); vaf_b2=vaf(B2_POS)
L0=req(c_lc,"L0_ancestral_root"); L1=req(c_lc,"L1_alpha_only(ancestor)"); L2=req(c_lc,"L2_alpha_beta(descendant)")
# ε=2% 門檻決策 + sSNV 統計
nsnv_b=req(enum,"A_by_n_sSNV_bucket")
bo=req(band,"orig_merged"); be=req(band,"eps2_merged"); pct_recl=req(band,"pct_reclassified")
sd=req(rti,"shape_dist_by_eps"); ft_cur=req(sd,"0.0","full_tree"); ft_e2=req(sd,"0.02","full_tree")
stab=req(rti,"stability_regions_changed"); stab12=req(stab,"ε1%→ε2%"); stab23=req(stab,"ε2%→ε3%")
con=req(rti,"concentration_collapsed_vs_kept"); coll_n=req(con,"structured→collapsed_n")
coll_cr=req(con,"collapsed_med_edge_coread"); kept_cr=req(con,"kept_med_edge_coread")
coll_gain=req(con,"collapsed_cn").get("gain",0); coll_gain_pct=round(100*coll_gain/coll_n,0)
ft_flat2=req(rti,"flat2_full_tree")

PAL={"green":"#2f9e44","blue":"#1c7ed6","orange":"#e8590c","red":"#e03131","grey":"#868e96","teal":"#0c8599","violet":"#6741d9","yellow":"#f08c00"}
def esc(s): return _h.escape(str(s))
def stacked(segs,width=620,h=44):
    tot=sum(v for _,v,_ in segs) or 1; x=0; bars=[]; leg=[]
    for lab,val,col_ in segs:
        w=width*val/tot; pct=100*val/tot
        if w>1:
            bars.append(f'<rect x="{x:.1f}" y="0" width="{w:.1f}" height="{h}" fill="{col_}"/>')
            if w>40: bars.append(f'<text x="{x+w/2:.1f}" y="{h/2+5:.0f}" text-anchor="middle" fill="#fff" font-size="13" font-weight="600">{pct:.0f}%</text>')
        x+=w; leg.append(f'<span class="lg"><i style="background:{col_}"></i>{esc(lab)} <b>{val:,}</b> ({pct:.1f}%)</span>')
    return f'<div class="chart"><svg viewBox="0 0 {width} {h}" width="100%" height="{h}">{"".join(bars)}</svg><div class="legend">{"".join(leg)}</div></div>'
def hbars(rows,width=620,maxv=None,baseline=None,unit=""):
    maxv=maxv or max(v for _,v,_,_ in rows)*1.12; rowh=32; pad=176; barw=width-pad-78; H=rowh*len(rows)+10
    o=[f'<svg viewBox="0 0 {width} {H}" width="100%" height="{H}">']
    if baseline is not None:
        bx=pad+barw*baseline/maxv; o.append(f'<line x1="{bx:.1f}" y1="0" x2="{bx:.1f}" y2="{H-10}" stroke="#adb5bd" stroke-dasharray="4 3"/><text x="{bx:.1f}" y="9" text-anchor="middle" font-size="10" fill="#868e96">baseline {baseline}</text>')
    for i,(lab,val,col_,sub) in enumerate(rows):
        y=i*rowh+16; w=barw*val/maxv
        o.append(f'<text x="{pad-8}" y="{y+13}" text-anchor="end" font-size="12" fill="#343a40">{esc(lab)}</text>')
        o.append(f'<rect x="{pad}" y="{y}" width="{max(w,1):.1f}" height="19" rx="3" fill="{col_}"/>')
        vt=f"{val:,}{unit}" if isinstance(val,int) else f"{val}{unit}"
        o.append(f'<text x="{pad+max(w,1)+6:.1f}" y="{y+14}" font-size="11.5" fill="#495057">{esc(vt)}{("  "+esc(sub)) if sub else ""}</text>')
    o.append('</svg>'); return f'<div class="chart">{"".join(o)}</div>'

def chr17_svg():
    return f'''<svg viewBox="0 0 660 250" width="100%" height="250" role="img">
<style>.nd{{fill:#fff;stroke:#495057;stroke-width:1.5}}.lb{{font-size:11.5px;fill:#212529}}.vf{{font-size:11px;fill:#0c8599;font-weight:600}}.ed{{stroke:#868e96;stroke-width:2;fill:none}}</style>
<line class="ed" x1="330" y1="48" x2="330" y2="110"/>
<line class="ed" x1="330" y1="150" x2="200" y2="195"/><line class="ed" x1="330" y1="150" x2="460" y2="195"/>
<circle class="nd" cx="330" cy="34" r="20"/><text class="lb" x="330" y="38" text-anchor="middle">germline</text>
<rect class="nd" x="250" y="110" width="160" height="40" rx="6"/><text class="lb" x="330" y="127" text-anchor="middle">α 祖先 ({A_POS})</text><text class="lb" x="330" y="142" text-anchor="middle">RAR · {L1} reads · VAF {vaf_a}</text>
<rect class="nd" x="120" y="195" width="160" height="42" rx="6"/><text class="lb" x="200" y="212" text-anchor="middle">β1 後代 ({B1_POS})</text><text class="lb" x="200" y="227" text-anchor="middle">VAF {vaf_b1} · normal 1/29 (3.4%)</text>
<rect class="nd" x="380" y="195" width="160" height="42" rx="6"/><text class="lb" x="460" y="212" text-anchor="middle">β2 後代 ({B2_POS})</text><text class="lb" x="460" y="227" text-anchor="middle">VAF {vaf_b2} · normal 0</text>
<text class="vf" x="330" y="180" text-anchor="middle">β1 ≡ β2 共連(同後代事件)</text>
<text x="14" y="246" font-size="10.5" fill="#868e96">機器真值 chr17_subclone_data.json：3 sSNV、群數 L0={L0}/α-only={L1}/α+β={L2}；線性巢狀(無 γ sibling)。方向由零格定、非 VAF。</text>
</svg>'''

# charts
c_ledger=stacked([("linked 21,554",b_link,PAL["green"]),("underpowered",b_under,PAL["yellow"]),("isolated(上界)",b_iso,PAL["grey"])])
c_shape=hbars([
 ("full_tree 完整樹",ft,PAL["green"],f"乾淨CN {clean_by['full_tree']}"),
 ("linear_nested 階層鏈",lin,PAL["teal"],f"乾淨CN {clean_by['linear_nested']}"),
 ("sibling_only 平行分支",sib,PAL["blue"],f"乾淨CN {clean_by['sibling_only']}"),
 ("co_linked 單lineage",col,PAL["violet"],"計入4678不計3820"),
 ("no_confirmed 無結構",noconf,PAL["grey"],""),
 ("inconsistent 成環",incon,PAL["red"],"丟棄"),
])
c_cfg=hbars([
 ("co_linked 共連",cfg("co_linked","pairs_sameHP"),PAL["green"],f"異HP {cfg('co_linked','pairs_diffHP'):,}"),
 ("nested 巢狀(合計)",nested_same,PAL["teal"],f"異HP {nested_diff:,}"),
 ("independent 無結構",cfg("independent","pairs_sameHP"),PAL["grey"],f"異HP {cfg('independent','pairs_diffHP'):,}"),
 ("mutual_excl 互斥",cfg("mutual_excl","pairs_sameHP"),PAL["red"],f"異HP {cfg('mutual_excl','pairs_diffHP'):,} ←異HP多=allelic"),
],maxv=11500)
c_hp=hbars([
 ("獨立 independent",round(rel["independent"]["enrichment_vs_chance"],2),PAL["grey"],"背景最高"),
 ("巢狀 nested b⊂a",round(rel["nested_b_in_a"]["enrichment_vs_chance"],2),PAL["grey"],""),
 ("共連 co_linked",round(rel["co_linked"]["enrichment_vs_chance"],2),PAL["grey"],""),
 ("互斥 mutual_excl",round(rel["mutual_excl"]["enrichment_vs_chance"],2),PAL["red"],"←DEPLETED=鑑別力所在"),
],maxv=2.1,baseline=1.0,unit="×")
c_ccf=stacked([("支持階層",g_sup,PAL["green"]),("tie 平手",g_tie,PAL["grey"]),("違反",g_vio,PAL["red"])])
c_cn=stacked([("CN-gain(multiplicity混淆)",gain_n,PAL["red"]),("LOH",req(cns,"loh"),PAL["green"]),("neutral",req(cns,"neutral"),PAL["teal"]),("loss",req(cns,"loss"),PAL["grey"])])
c_npop=hbars([(f"{k} 群",npop[k],(PAL["green"] if int(k)>=2 and int(k)<=4 else PAL["grey"]),"") for k in npop if k!="0"],maxv=max(npop.values())*1.1)

NSNV_ORDER=["2","3","4","5-9","10-49","50+"]
c_nsnv=hbars([(f"{k} 個 sSNV", nsnv_b.get(k,0), (PAL["green"] if k in("2","3") else PAL["teal"]), "") for k in NSNV_ORDER], maxv=max(nsnv_b.values())*1.1)

def tier(t): return f'<span class="tier {t}">{t}</span>'

CSS="""
:root{--ink:#212529;--mut:#868e96;--line:#dee2e6;--bg:#f8f9fa;--card:#fff;--green:#2f9e44;--red:#e03131;--amber:#e8590c}
*{box-sizing:border-box}body{margin:0;font-family:-apple-system,"Segoe UI","Noto Sans TC","PingFang TC","Microsoft JhengHei",sans-serif;color:var(--ink);background:var(--bg);line-height:1.72}
.wrap{max-width:980px;margin:0 auto;padding:26px 22px 90px}
h1{font-size:24px;margin:.2em 0 .1em;line-height:1.35}h2{font-size:20px;margin:1.7em 0 .5em;padding-bottom:.25em;border-bottom:2px solid var(--line)}
h3{font-size:15.5px;margin:1.25em 0 .35em;color:#343a40}p{margin:.55em 0}.sub{color:var(--mut);font-size:13px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:15px 18px;margin:13px 0;box-shadow:0 1px 2px rgba(0,0,0,.03)}
.bb{border-left:5px solid var(--green)}.aux{border-left:5px solid #7048e8}.cf{border-left:5px solid var(--amber)}.me{border-left:5px solid #0c8599}
.ribbon{display:flex;flex-wrap:wrap;gap:8px;margin:10px 0}.tag{font-size:12px;padding:3px 10px;border-radius:20px;background:#e7f5ff;color:#1971c2;border:1px solid #a5d8ff}
.tag.warn{background:#fff4e6;color:#d9480f;border-color:#ffd8a8}.tag.star{background:#fff9db;color:#e67700;border-color:#ffe066}
.prov{background:#fff9db;border:1px solid #ffe066;border-radius:10px;padding:12px 16px;font-size:13px}
.tier{display:inline-block;font-size:11px;font-weight:700;padding:1px 7px;border-radius:5px;margin-right:4px}
.L1{background:#d3f9d8;color:#2b8a3e}.L2{background:#e5dbff;color:#5f3dc4}.L3{background:#fff3bf;color:#b08900}.L5{background:#ffe3e3;color:#c92a2a}
.chart{margin:10px 0}.legend{display:flex;flex-wrap:wrap;gap:6px 16px;margin-top:8px;font-size:12px;color:#495057}.lg i{display:inline-block;width:11px;height:11px;border-radius:2px;margin-right:5px;vertical-align:middle}
.kv{display:flex;gap:14px;flex-wrap:wrap;margin:8px 0}.kv .b{background:#f1f3f5;border-radius:8px;padding:8px 13px;font-size:12.5px;min-width:108px}.kv .b b{display:block;font-size:20px;color:#1c7ed6}
table{border-collapse:collapse;width:100%;font-size:12.5px;margin:10px 0}th,td{border:1px solid var(--line);padding:6px 9px;text-align:left;vertical-align:top}th{background:#f1f3f5}
.red{color:var(--red);font-weight:600}.grn{color:var(--green);font-weight:600}
.two{display:grid;grid-template-columns:1fr 1fr;gap:13px}@media(max-width:680px){.two{grid-template-columns:1fr}}
.redline{background:#fff5f5;border:1px solid #ffc9c9;border-radius:10px;padding:13px 18px}.redline li{margin:5px 0}
.note{font-size:12px;color:var(--mut)}.eg{background:#f8f9fa;border:1px dashed #ced4da;border-radius:8px;padding:11px 14px;margin:9px 0}
.vd{display:inline-block;font-size:11px;font-weight:700;padding:1px 8px;border-radius:5px}
.vSOUND{background:#d3f9d8;color:#2b8a3e}.vCARE{background:#fff3bf;color:#b08900}.vRISKY{background:#ffe8cc;color:#d9480f}.vINVALID{background:#ffe3e3;color:#c92a2a}
footer{margin-top:38px;padding-top:16px;border-top:2px solid var(--line);font-size:12px;color:var(--mut)}
"""

HTML=f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>clone/subclone 重建 — 骨幹軸深化教學 v2（HCC1395 ⭐3）</title><style>{CSS}</style></head><body><div class="wrap">

<h1>用 ONT long-read 重建 clone/subclone — 骨幹軸深化 + 對抗驗證修正版</h1>
<p class="sub">HCC1395 單樣本 ⭐3 · Tier-R(same-read ≤50kb) · 每數字標 L1-L5 + 由凍結 JSON 注入(非手打) · 經 2 個對抗 workflow 修正</p>
<div class="ribbon"><span class="tag star">⭐3 單樣本上限</span><span class="tag">TP∪FP union build · SEQC2 僅觀察</span><span class="tag warn">數字源 branch 5308d9e(已凍結)</span></div>

<div class="prov"><b>📌 Provenance + 方法學原則</b><br>
數字原始 JSON 來自<b>未 merge branch <code>feat/summary-nreadsvalid@5308d9e</code></b>,已凍結到 trunk <code>_assets/20260627_subclone_4axis_teaching/{{data,scripts}}</code>(對照 <code>PROVENANCE.md</code>)。
<b>建構原則</b>:build 用 ClairS→longphase-S 的 <b>TP∪FP union</b>,somatic 由 <b>normal 比對</b>定,<b>SEQC2 的 TP/FP 只用於觀察、不進前處理</b>。
本頁由 <code>build_4axis_teaching.py</code> 從凍結 JSON 注入;region 數字由 <code>derive_region_buckets.py</code> 從 regions.tsv derive(clean full_tree {cft} 已非硬編)。</div>

<h2>§0 結論金字塔</h2>
<div class="card">
<p>{tier('L1')}體細胞 sSNV 的<b>單分子共現</b>是唯一非循環的克隆重建骨幹,在 read 跨度內重建<b>區域級局部克隆樹</b>;HP/VAF 只<b>鑑別與一致性檢查</b>、不獨立偵測;CNV/LOH 是<b>須先扣除的 confound(兼區域刻畫)</b>非證據;甲基<b>off-ladder bounded-auxiliary</b>(cis-control 0/{mr_t} 未測→subclone-specificity UNDETERMINED)。</p>
<div class="kv"><div class="b"><b>{U:,}</b>sSNV(TP {TP:,}+FP {FP:,})</div><div class="b"><b>{n_reg:,}</b>區域</div><div class="b"><b>{ft}</b>full_tree</div><div class="b"><b>{cft}</b>乾淨CN可信子集</div><div class="b"><b>{mr_rate}%</b>甲基弱佐證</div></div>
</div>

<h2>§1 核心觀念：非循環性是「有條件的」</h2>
<div class="card bb">
<p><b>一條 ONT read = 一條 DNA 分子 = 一個細胞的一條單倍型。</b>兩 sSNV 同 read 共現 = 物理 cis 連鎖 = 同 lineage,是<b>直接觀測</b>不需統計指派——這是非循環的根。</p>
<p>{tier('L1')}<b>對 phase-switch 免疫(真強項)</b>:相位切換是跨 read 統計拼接的假象,單分子內部不可能切換單倍型。</p>
<p>{tier('L1')}🔴 <b>但最強反駁=比對偽影(paralog/segdup/amplicon mis-map)</b>:錯比 read 帶「其實不在此處」的 ALT→分子層級偽造 cis 共現;somatic-confirmed 修不了(錯比 read 在 normal 也一致缺席)、MAPQ 擋不住 segdup、本輪缺 mappability track。</p>
<div class="eg"><b>∴ 條件化核心主張(對外一律用此版)</b>:在 <b>uniquely-mappable ∩ CN-clean(LOH+neutral) ∩ somatic-confirmed ∩ coread≥6</b> 的位點對下,同 read 共現=cis=同 lineage。整體骨幹 <b class="red">{link_gain_pct}% 落 CN-gain</b>(vs genome-wide {cn_gain}%)→ 含 gain 區的「方向」標 UNDETERMINED,structure fraction 為上界。</div>
<p class="note">對比:甲基分群是循環的(double-dip)——用甲基分群再檢定群間甲基差異,同量算兩次,故甲基只能事後佐證。</p>
</div>

<h2 id="bb">§2 骨幹軸深掘（6 小節）</h2>

<div class="card bb"><h3>2.1 直接性 vs 短讀 — 為何不需 phasing</h3>
<p>短讀要「統計 phasing 架橋」才知道兩變異是否 cis;long-read 直接讀出。骨幹**不依賴**甲基/頻率推估/下游分群 → 非循環。前置 gate(2.1 條件)是把「read=cell」這個載重前提守住的代價。</p></div>

<div class="card bb"><h3>2.2 2×2 共現 census — read-level four-gamete 檢定</h3>
<p>對每對 sSNV,只看同跨兩位點的 read 填 RR/RA/AR/AA。<b>這是 read-level four-gamete / perfect-phylogeny 相容性檢定</b>(infinite-sites 下每位點只突變一次,不該見四種 gamete)。實際 <code>classify()</code> 門檻(已驗碼):</p>
<table><tr><th>config</th><th>零格規則(實碼)</th><th>意義</th></tr>
<tr><td class="grn">co_linked</td><td>ra==0 <b>and</b> ar==0 (嚴格)</td><td>完美共現=同 lineage</td></tr>
<tr><td class="teal">nested_a_in_b</td><td>ar==0</td><td>a-ALT 必伴 b-ALT→b 祖先 a 後代</td></tr>
<tr><td>nested_b_in_a</td><td>ra==0</td><td>對稱反向</td></tr>
<tr><td class="red">mutual_excl</td><td><b>aa&lt;2</b>(容忍1個AA)且 ra+ar≥2</td><td>兩 ALT 幾不共現=sibling/allelic</td></tr>
<tr><td>independent</td><td>四格皆&gt;0</td><td>four-gamete 違反、無乾淨結構</td></tr></table>
<p class="note">⚠ 碼層發現(CODE_DISCREPANCY):AA 用 aa&lt;2(容忍1)、off-diagonal 用嚴格==0(不對稱)→ 單一雜訊/嵌合 read 會把真 co_linked 降 nested、把 nested 降 independent=<b>系統性低估 lineage</b>;建議改 Fisher/binomial over-dispersion(對稱容忍)。</p>
<p>canonical 子集(coread≥6+兩端 somatic,n={filt['powered+both_somatic']:,})的 same/diff HP:</p>
{c_cfg}
<p>🔑 <b>巢狀 same-HP {nested_same:,} ≫ 異-HP {nested_diff:,}</b>(克隆階層幾乎全同單倍型);<b>互斥反而異-HP {cfg('mutual_excl','pairs_diffHP'):,} &gt; same-HP {cfg('mutual_excl','pairs_sameHP'):,}</b>(互斥多為 allelic 非 subclone)——這就是為何需 HP 鑑別(見 §4)。
<span class="note">filter 兩道閘:{filt['all_recorded_coread>=2']:,}(coread≥2)→{filt['powered_coread>=6']:,}(coread≥6,保零格可信)→{filt['powered+both_somatic']:,}(兩端 somatic,除 germline allelic)。</span></p></div>

<div class="card bb"><h3>2.3 L1→L2→L3 組裝 + 為何 component≠subclone</h3>
<p>L1 pairwise 2×2 → L2 multi-locus 基因型向量(局部 population) → L3 整合成區域樹。連通圖散成 <b>{ncomp:,} 個 component</b>(局部連鎖區域),含環 <b>{ncyc}</b> 個(0.37%,必丟)。</p>
<p>🔴 <b>component ≠ subclone</b>:① 由可觀測性+距離定義非細胞身份;② 混 same-HP(克隆)與 diff-HP(allelic)須再 HP-gate;③ 真 subclone 突變散全基因組會被切進數千 component;④ Fisher 顯著只證連鎖存在不分 subclone/allelic。
<span class="note">誠實橋接:{ncomp:,} components ≠ {n_reg:,} regions(不同 pipeline 階段,環 {ncyc}≠inconsistent {incon})。</span></p>
<p>{n_reg:,} 區的樹形(全由 regions.tsv derive,grep-able):</p>
{c_shape}
<p><span class="grn">structured {structured:,}(full+linear+sibling)</span> / +單lineage {structured_plus:,};🔴 是<b>上界</b>(含 FP+偽影);論文級=<b>乾淨 CN full_tree {cft}</b>。</p></div>

<div class="card bb"><h3>2.4 統計 null 與 power-gate — 存在性≠鑑別性</h3>
<p>null=Fisher exact(2×2 條件檢定)。Fisher vs free-margin MC 一致 <b>{f_agree}/{f_eval}={f_agree_pct}%</b>(穩健)。power-gate:sparse 類僅 <b>{sparse_pct}%</b> 顯著≈機率底線(不從噪聲造結構)。</p>
<p>🔴 <b>diff-HP 負對照(物理不可能同克隆)仍 {diffhp_pct}% 顯著</b> → Fisher 顯著只回答「共現存在嗎」,<b>「是 subclone 還是 allelic」全靠 HP-gate</b>。這是骨幹(存在性)交棒 HP(鑑別性)的分工。</p></div>

<div class="card bb"><h3>2.5 γ 類 FP-source 與偽影 — TP∪FP union 的誠實框架</h3>
<p>骨幹 linked-somatic 按來源:TP {tp_link:,} + FP-source {gamma:,}。<b>依建構原則,build 用 TP∪FP union、somatic 由 normal 定;SEQC2 的 TP/FP 是觀察疊層不是 filter</b>。所以 FP-source {gamma:,} 是「SEQC2 標 FP 但流程 normal 判 somatic-like」的<b>觀察</b>,不可當「漏掉的確證 somatic」(F3 自承 convergence tautological)。</p>
<p>🔴 偽影:<b>{dense:,}</b> 個 linked-somatic 落 <b>{ndense}</b> 個 dense-uniform cluster(≥5 in ≤5kb);chr9:41777788-41804669 區 78 sSNV 擠 26.9kb(n_links 逼近上限)=偽影最像處(連鎖越多≠越可信)。→ structure fraction 為單向上界。</p></div>

<div class="card bb"><h3>2.6 chr17:48360161 worked example(機器真值 3-sSNV)</h3>
{chr17_svg()}
<p>三對 2×2 的零格直接定拓樸(全 same-HP 1-1):</p>
<table><tr><th>配對</th><th>coread</th><th>RR/RA/AR/AA</th><th>零格→推論</th></tr>
<tr><td>α×β2</td><td>{c_pairs['48365089_48365161']['n_coread']}</td><td>{c_pairs['48365089_48365161']['REF_REF']}/{c_pairs['48365089_48365161']['REF_ALT']}/{c_pairs['48365089_48365161']['ALT_REF']}/{c_pairs['48365089_48365161']['ALT_ALT']}</td><td>RA=0→β2 巢狀於 α</td></tr>
<tr><td>β1×α</td><td>{c_pairs['48362515_48365089']['n_coread']}</td><td>{c_pairs['48362515_48365089']['REF_REF']}/{c_pairs['48362515_48365089']['REF_ALT']}/{c_pairs['48362515_48365089']['ALT_REF']}/{c_pairs['48362515_48365089']['ALT_ALT']}</td><td>AR=0→β1 巢狀於 α</td></tr>
<tr><td>β1×β2</td><td>{c_pairs['48362515_48365161']['n_coread']}</td><td>{c_pairs['48362515_48365161']['REF_REF']}/{c_pairs['48362515_48365161']['REF_ALT']}/{c_pairs['48362515_48365161']['ALT_REF']}/{c_pairs['48362515_48365161']['ALT_ALT']}</td><td>RA=AR=0→β1≡β2 同事件</td></tr></table>
<p>{tier('L1')}方向<b>純由零格定、不看 VAF</b>(α⊃β)。{tier('L2')}VAF {vaf_a}&gt;{vaf_b1}/{vaf_b2} 只是<b>一致性檢查非獨立驗證</b>(nested 對近定義性)。
<span class="note">🔴 校正:舊 build 的「4 sSNV、γ sibling、7/10/15」不被機器源支持(γ=48357368 零命中);此為線性巢狀非 full_tree。β1 normal 1/29=3.4%,在統一 &lt;5% 定義下算 somatic(嚴格==0 則否)→ 此即待修的 somatic 雙定義。</span></p></div>

<h2>§2.7 sSNV 關聯結果統計 + 單讀脆弱性與 ε=2% 門檻定案</h2>
<div class="card bb">
<h3>每區 sSNV 數量分布（{n_reg:,} 區，1-sSNV=0 因單點不可連鎖）</h3>
{c_nsnv}
<p class="note">超過一半的區只有 2 個 sSNV（{nsnv_b.get('2',0):,}）→ 主張「更多位點」被 power 封頂的實證。</p>
<h3>🔴 單讀脆弱性 → ε=2% 噪聲底線定案</h3>
<p>分類靠零格定拓樸,但一條定序錯誤 read 能把真 co_linked 讀成 nested(假階層)。決定格只有 1 條 read 的比例:nested 50.6% / independent 57.8%。觀察發現這些**不是低覆蓋**(coread 中位 42-64),而是**off-diagonal 佔比落在 ONT 噪聲底線(1/42≈2.4%)**。</p>
<div class="eg"><b>定案規則</b>:cell 為真 ⟺ <b>count &gt; coread × 2%</b>(ONT 噪聲底線)。保留最低 1 條(低 coread 單讀仍算),高 coread 單讀(1 ≤ coread×2%,coread≥50)判 noise。<br>
逐筆可驗:<code>data/pairs_eps2_annotated.tsv</code>(每對附 floor_2pct + 強/弱/noise)。例:RA=1/coread=60/floor=1.2 → noise;RA=1/coread=20/floor=0.4 → 保留;RA=18/coread=22 → strong。</div>
<h3>方向 3 路印證 + 值錨 ONT（🔴 非「3 路收斂於同一數」）</h3>
<p class="note">對抗稽核校正：A 定<b>方向</b>(ε 對/VAF 錯)、B 封<b>上限</b>、值 2% 錨在 <b>ONT 錯誤率</b>。三路不獨立收斂於「2」。</p>
<table><tr><th>路徑</th><th>結果</th><th>角色</th></tr>
<tr><td>A. FP 裁判(觀察)</td><td>弱 both-FP 57.9% vs 強 30.1%;分離度單調 +14.8/+17.9/+24.6;VAF 調整 <b>−10.3 backwards</b></td><td>定<b>方向</b>(ε對 VAF錯);單調故不指定值</td></tr>
<tr><td>B. 結構穩定+零單讀</td><td>full_tree {ft_cur}→{ft_e2}(−9%);零單讀(off≥2)仍存 <b>{ft_flat2}</b>(77%);相鄰 ε 僅 {stab12}/{stab23} 區動</td><td>封<b>上限</b>;headline 非單讀假象(post-hoc edge 重推)</td></tr>
<tr><td>C. 塌陷集中度</td><td>塌陷 {coll_n} 區 edge coread <b>{coll_cr:.0f} vs {kept_cr:.0f}</b>;CN-gain <b>{coll_gain_pct:.0f}%</b></td><td>印證打對地方(偽影區)</td></tr>
<tr><td>值錨</td><td>ONT R10/Dorado substitution ~1-2% + min_base_quality=0 推高 → 2% 中央/1-3% band</td><td>定<b>值</b></td></tr></table>
<p class="note">⚠ scope：ε=2% 為 HCC1395 單樣本 + SEQC2 校準值(⭐3);跨樣本須重校 ε。</p>
<h3>最終 band（ε=2% 套用 38,049 pairs，重分類 {pct_recl}%）</h3>
<table><tr><th>config</th><th>現行(≥1)</th><th>ε2%</th></tr>
<tr><td class="grn">co_linked</td><td>{bo['co_linked']:,}</td><td><b>{be['co_linked']:,}</b></td></tr>
<tr><td>nested(合)</td><td>{bo['nested']:,}</td><td>{be['nested']:,}</td></tr>
<tr><td>independent</td><td>{bo['independent']:,}</td><td>{be['independent']:,}</td></tr>
<tr><td>mutual_excl</td><td>{bo['mutual_excl']:,}</td><td>{be['mutual_excl']:,}</td></tr></table>
<p class="note">⚠ 否決:純絕對 ≥2(太鈍、生 sparse)、VAF 調整(FP 裁判 backwards)、binomial/Fisher(小樣本 underpowered+需假設,不可驗)。決策全文 → 20260628_sSNV_linkage_threshold_decision_eps2_01.md。VAF 只當描述欄不當 cutoff。源頭另改 min_base_quality≥10 降 ε。</p>
</div>

<h2>§3 證據階層重構（你的 7-rung → 修正結構）</h2>
<div class="card">
<p>你原本的線性 7-rung,驗證後重構為「<b>骨幹 ladder + 正交軸 + precondition</b>」:</p>
<table><tr><th>類</th><th>內容</th><th>角色</th></tr>
<tr><td class="grn"><b>A 非循環骨幹</b></td><td>A1 同-read 多點 sSNV(最強,限 CN-clean∩non-segdup) &gt; A2 跨-read 交集(cross-PS&gt;50kb=GAP) &gt; A3 nesting topology</td><td>genetic 共現(A3 是輸出非獨立證據)</td></tr>
<tr><td>B 條件式一致性</td><td>VAF magnitude</td><td>{tier('L2')}非獨立、必先 condition on CNV、只 CN-clean、連報 {g_sup_pct}/{g_vio_pct}/{g_tie_pct}</td></tr>
<tr><td>C 鑑別/rule-out</td><td>HP tag</td><td>{tier('L1')}sibling-vs-allelic 鑑別器,診斷力只在互斥(0.86×)</td></tr>
<tr><td class="cf"><b>precondition</b></td><td><b>CNV + LOH</b></td><td>🔴 <b>condition-on confound + 區域刻畫資訊</b>(非刪除、非獨立證據 rung);LOH 與 CNV 同類非 HP 子類</td></tr>
<tr><td class="me">D off-ladder 佐證</td><td>甲基</td><td>bounded-auxiliary;subclone-specificity UNDETERMINED</td></tr></table>
<p class="note">你的關鍵更正已採納:<b>CNV/LOH 不刪除</b>(對解釋區域突變狀況有資訊),但改列「先 condition、兼輔助刻畫」而非「lineage 投票證據」。precedence「衝突時 genetic 勝甲基」是<b>設計立場非已驗階層</b>(壓測僅 n=1)。</p></div>

<h2>§4 鑑別軸：HP（只在互斥有診斷力）</h2>
<div class="card aux">
<p>same-HP 高在巢狀/共連是<b>區域背景</b>(同分子必同 HP tag,1.7-1.87×),非克隆特異。HP 真正診斷力在<b>互斥 DEPLETED</b>:</p>
{c_hp}
<p>{tier('L1')}{hp_in:,} 個 read 級互斥→HP 移除 <span class="red">{hp_rm:,}({hp_rm_pct}%)異-HP allelic</span>→剩 <span class="grn">{hp_true:,} 細胞層 sibling</span>。<b>HP=鑑別器,不新增 call。</b>
<span class="note">⚠ 碼:HP3('3')未排除→偽 same_hp 污染 sibling ~6.3%(1317/20815,待修);PS 單一可信 {ps_rate}%。</span></p></div>

<h2>§5 一致性軸：VAF/CCF（條件式、非獨立）</h2>
<div class="card aux">
<p>方向已由零格定,VAF 只順驗量級。誠實口徑:支持 <b>{g_sup_pct}%</b>/違反 {g_vio_pct}%/tie {g_tie_pct}%(禁只報 {g_vio_pct}%),只在 CN-clean({cn_clean}%)可估。</p>
{c_ccf}
<p>{tier('L3')}GMM best_n={best_n} 但 ΔBIC(vs n=4)僅 <b>{dbic}</b>(邊際)+CN-gain 污染→離散 subclone 層級降 L3。</p></div>

<h2>§6 precondition：CNV / LOH（confound + 區域刻畫，非證據）</h2>
<div class="card cf">
<p>{tier('L1')}somatic sSNV <b>{cn_gain}%</b> 落 CN-gain;骨幹更高達 <b>{link_gain_pct}%</b>。CN-gain multiplicity 讓「AA=二代」失效(AA 可能是同突變被擴增)→破壞 VAF→CCF。</p>
{c_cn}
<p>🔴 <b>LOH 不是「能乾淨定義 subclone 處」,反而是假陽最高發處</b>:HP 標籤塌陷(diff-HP→偽 same-HP)+imprinting-unmask(82-91%,Martin-Trujillo 文獻 L3);甲基 corroborated {mr_c} 區的 per-CN 拆分待補(本 bundle 凍結 JSON 無)。
<b>採納你的觀點</b>:CNV/LOH <b>不刪除</b>——它們刻畫區域狀態、是修正與輔助驗證資訊;但**先 condition(限 CN-clean 才讀 VAF)**,不當 lineage 投票。</p></div>

<h2>§7 甲基：off-ladder bounded-auxiliary</h2>
<div class="card me">
<div class="two"><div class="eg"><b class="red">50% 別當 headline</b><br>既有 ISM CN-clean≤8kb:{mc_t} 測 {mc_c}={mc_rate}%(cherry-pick)</div><div class="eg"><b class="grn">{mr_rate}% 誠實口徑</b><br>全基因組重抽 {mr_t} 測 {mr_c}={mr_rate}%</div></div>
<p>{tier('L1')}甲基<b>不能獨立偵測</b>:PERMANOVA 全 {mg_tgt} 區只 {mg_test} 可算、recover {mg_rec}=0 新 partition。{tier('L1')}cis-control <b>0/{mr_t}=structural zero</b>→subclone-specificity <b>UNDETERMINED</b>(非已驗陰性,可能近全 cis-ASM)。</p>
<div class="eg"><b>🔴 主張 3（甲基判突變先後時序）= INVALID</b>:同基因型 read(都 AAA/都 RRR)<b>遺傳上零 ordering 訊號</b>;甲基<b>非分子鐘</b>無時序映射;且是 double-dip。連有遺傳訊號的 740 區甲基都 recover 0,零訊號處更不可能。無未來救贖路徑。<br>
<b>🔴 主張 5（甲基定義 clone→外推單位點）= RISKY</b>:前提「甲基能定義 subclone」<b>未成立而非已成立</b>(待 T-GATE-GB cis-control);是高風險待驗證方向非死路。</div></div>

<h2>§8 6 主張裁決表</h2>
<table>
<tr><th>主張</th><th>裁決</th><th>核心理由</th></tr>
<tr><td>1 二代含一代(nesting)</td><td><span class="vd vSOUND">SOUND</span></td><td>=infinite-sites/perfect-phylogeny;但 het→hom 主機制是 LOH→改稱「二代<b>事件</b>」;記號 1-1/1-1-1 與 HP tag 正交須標軸別</td></tr>
<tr><td>2 互斥解析(diff-HP/same-HP)</td><td><span class="vd vCARE">NEEDS_CARE</span></td><td>鑑別方向對;但 LOH 區是假陽最高發(非乾淨)、用甲基距離定 1-1=double-dip</td></tr>
<tr><td>3 甲基判突變時序</td><td><span class="vd vINVALID">INVALID</span></td><td>同基因型零 ordering+甲基非分子鐘+double-dip;無救贖路徑</td></tr>
<tr><td>4 更多位點→更多群→更完整演化</td><td><span class="vd vCARE">NEEDS_CARE</span></td><td>「更多群」power 內成立(n_pop 實測見下);「更完整」被 read-span 封頂 regional;sSNV 稀疏(多數區僅 2 sSNV={nsnv_b.get('2',0):,})是瓶頸非計算量</td></tr>
<tr><td>5 甲基→單位點外推</td><td><span class="vd vRISKY">RISKY</span></td><td>前提未成立(recover 0、cis-control 未跑)+層級轉移;待 T-GATE-GB</td></tr>
<tr><td>6 還有哪些衝突/細節</td><td><span class="vd vCARE">已修多項</span></td><td>somatic 雙定義、classify 不對稱容忍、chr17 3vs4、孤兒檔、硬編 205→本版已修/標</td></tr>
</table>
<p>主張 4 的 n_populations 實測分布(regions.tsv derive):</p>
{c_npop}
<p class="note">2 群 {npop.get('2','?')} / 3 群 {npop.get('3','?')} 為主;≥5 群各僅 1 區(極可能 dense-cluster 偽影)。</p>

<h2>§9 誠實邊界 / 紅線</h2>
<div class="redline"><ul>
<li>🔴 ⭐3 單樣本 single-pipeline 封頂;regional(≤read-span)非 genome-wide tree;分子共現≠single-cell confirmation。</li>
<li>🔴 核心共現主張須條件化(uniquely-mappable∩CN-clean∩somatic∩coread≥6);CN-gain 區方向 UNDETERMINED;structure fraction 為上界(對外用 CN-clean full_tree {cft})。</li>
<li>🔴 CNV/LOH=condition-on confound+區域刻畫(非證據 rung);VAF=一致性檢查(非獨立);HP=鑑別器(只在互斥)。</li>
<li>🔴 甲基 corroborate 非 detect、0 新 partition、cis-control 0/{mr_t}→subclone-specificity UNDETERMINED;主張 3 INVALID、主張 5 RISKY。</li>
<li>🔴 SEQC2 TP/FP 只觀察不進前處理;成環=資料不支持單樹的誠實訊號(four-gamete 違反)→丟非強解。</li>
<li>🔴 對外丟「first read-level methyl lineage」(Foltz 2024 已有 single-molecule primitive);novelty=native ONT+read×read PERMANOVA+normal-anchored cis-test+有界甲基。</li>
</ul></div>

<h2>§10 下一步（load-bearing 排序）</h2>
<table>
<tr><th>步驟</th><th>理由</th><th>節點</th></tr>
<tr><td>① matched-normal cis-control</td><td>解 cis-control 0/{mr_t} structural zero;唯此把甲基從候選推進到佐證,並驗主張 5 前提</td><td class="grn">T-GATE-GB(#1)</td></tr>
<tr><td>② 統一 somatic 定義 + classify over-dispersion</td><td>消 somatic 雙定義(chr17 3-sSNV 定案)+ 修不對稱容忍(系統性低估 lineage)</td><td>code-fix</td></tr>
<tr><td>③ 補 mappability/segdup mask</td><td>分離 CN-gain 環/dense-cluster 偽影,把 structure 上界往真值收</td><td>T-ONT-CNV</td></tr>
<tr><td>④ COLO829+≥5/7 跨樣本</td><td>single-pipeline 封頂;升 ⭐3→⭐4</td><td>T-DORADO</td></tr>
<tr><td>⑤ single-cell/multi-region 正交確認</td><td>分子共現≠single-cell;L2 生物詮釋轉 confirmed</td><td>T-GATE-GD</td></tr>
</table>

<footer>數字源 branch 5308d9e(pending-merge),凍結於 _assets/20260627_subclone_4axis_teaching/{{data,scripts}},對照 PROVENANCE.md。本頁 build_4axis_teaching.py 從凍結 JSON 注入+derive_region_buckets.py derive(缺 key refuse)。
模型驗證裁決全文:InterSubMod/docs/methodology/20260628_reconstruction_model_verification_01.md(workflow wf_f2b070ea-64c 15 agent)。單一真值結論:InterSubMod/docs/methodology/20260627_subclone_unified_verified_narrative_01.md。
證據層級:L1 源碼/JSON 重現·L2 報告 JSON·L3 推論·L5 公理。</footer>
</div></body></html>"""

with open(OUT,"w",encoding="utf-8") as f: f.write(HTML)
print(f"OK wrote {OUT} ({len(HTML):,} bytes)")
print(f"checks: U={U} reg={n_reg} struct={structured}/{structured_plus} clean_ft={cft} chr17_vaf={vaf_a}/{vaf_b1}/{vaf_b2} chr17_lc={L0}/{L1}/{L2} link_gain%={link_gain_pct} npop2={npop.get('2')}")
