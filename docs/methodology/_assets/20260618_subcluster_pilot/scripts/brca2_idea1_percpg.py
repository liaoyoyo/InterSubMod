import json
import numpy as np
from scipy.stats import fisher_exact
d=json.load(open('research/tsg_promoter_asm_reviewer/output/step4_ism_results.json'))
brca2=d[0]; assert brca2['candidate']['gene']=='BRCA2'
cr=brca2['cpg_records']
def bh(ps):
    p=np.array(ps); m=len(p); o=np.argsort(p); r=p[o]*m/np.arange(1,m+1)
    q=np.minimum.accumulate(r[::-1])[::-1]; out=np.empty(m); out[o]=np.clip(q,0,1); return out
MIN=5  # 每組最少 read
recs=[]
for c in cr:
    # carrier 軸: germline = HP1 (+HP2), somatic = HP1-1 (+HP2-1)
    gm_m=c['HP1_meth']+c['HP2_meth']; gm_u=c['HP1_unmeth']+c['HP2_unmeth']
    so_m=c['HP1-1_meth']+c['HP2-1_meth']; so_u=c['HP1-1_unmeth']+c['HP2-1_unmeth']
    gm_n=gm_m+gm_u; so_n=so_m+so_u
    if gm_n<MIN or so_n<MIN: continue
    gb=gm_m/gm_n; sb=so_m/so_n; db=sb-gb  # som - germ; 負=somatic hypo (符 -0.122)
    try: _,p=fisher_exact([[so_m,so_u],[gm_m,gm_u]])
    except: continue
    recs.append({'pos':c['cpg_pos'],'dist':c['dist_to_var'],'germ_b':gb,'som_b':sb,'dbeta':db,'p':p,'gn':gm_n,'sn':so_n})
ps=[r['p'] for r in recs]; qs=bh(ps)
for r,q in zip(recs,qs): r['q']=q
sig=[r for r in recs if r['q']<0.05]
neg=[r for r in sig if r['dbeta']<0]  # 符合 -0.122 方向
print(f"=== BRCA2 chr13:32315128 想法1 per-CpG carrier 軸定位 (germline vs somatic-HP) ===")
print(f"truth (stats HP1_vs_HP1-1): region Δβ=-0.1219, wilcoxon p=6.1e-11, somatic HYPOMETHYLATED")
print(f"decisionflow: tumor ⑤ aligned=carrier (F=15.4); 此 pilot 定位驅動的 CpG")
print(f"\nCpG 測試池 (兩組各≥{MIN} read): {len(recs)} / 412 cpg_records")
print(f"顯著 CpG (BH-FDR q<0.05): {len(sig)}")
print(f"  其中方向符 -0.122 (somatic hypo, Δβ<0): {len(neg)}/{len(sig)} ({100*len(neg)/len(sig):.0f}%)")
import numpy as np
adb=np.array([r['dbeta'] for r in sig])
print(f"  顯著 CpG |Δβ| 中位={np.median(np.abs(adb)):.3f} 範圍 [{adb.min():.3f},{adb.max():.3f}]")
print(f"  pool 全體 Δβ 均值 (加權近似) = {np.mean([r['dbeta'] for r in recs]):.4f} (對照 region -0.122)")
dists=[r['dist'] for r in sig]
print(f"  顯著 CpG 距變異位點 dist_to_var: 中位={int(np.median(dists))} 範圍[{min(dists)},{max(dists)}]")
print(f"\n--- 前 12 顯著 CpG (依 |Δβ| 排序) ---")
for r in sorted(sig,key=lambda x:-abs(x['dbeta']))[:12]:
    print(f"  pos={r['pos']} dist={r['dist']:+5d} germ_β={r['germ_b']:.2f} som_β={r['som_b']:.2f} Δβ={r['dbeta']:+.3f} q={r['q']:.1e} (g{r['gn']}/s{r['sn']})")
json.dump({'n_tested':len(recs),'n_sig':len(sig),'n_neg':len(neg),'sig':sig},open('/tmp/_brca2_idea1_out.json','w'))
