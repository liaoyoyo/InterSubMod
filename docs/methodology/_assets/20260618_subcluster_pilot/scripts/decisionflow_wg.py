#!/usr/bin/env python3
"""全基因組判別流程 compute (dual-mode: tumor-only + merged tumor/normal 分開算)。
每位點: ① outlier-tolerant 切群(T=3/4/5) + a-priori 對齊(離散結構) ; ② per-axis PERMANOVA+PERMDISP
(均值差/散開度, mean-shift 軸). 兩軸合併供 analyze 做完整 5-態分類。disk-safe 逐 chr。
輸出 decisionflow_records_{tumor,merged}.json。

5 態(analyze 端組裝)：①不可驗證(n_complete<6) ②1群/無訊號(切不出 & PERMANOVA 不顯著)
③監督可分 mean-shift(切不出但 PERMANOVA-FDR 顯著; 內分 location-clean vs dispersion)
④可分但未對齊(離散切得出但不對 a-priori = epiallele?雜訊?) ⑤確認真結構(離散切得出且對齊生物軸)。

③切群: valid群≥T(≤T-1=outlier), ≥2 valid群, outlier≤20%, valid-read silhouette 最佳。
⑤對齊: valid-cut vs HP/carrier/allele CramerV≥0.3 & chi2 p<.05 & Cochran e≥5。
PERMANOVA: distance-based pseudo-F (Anderson2001/McArdle), permutation p (nperm=99, 種子=pos)。
PERMDISP: PCoA distance-to-centroid betadisper (vegan 式 real/imag 軸), permutation p。"""
import os, csv, glob, json, subprocess, shutil, time
for _v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[_v]="1"  # 單線程 BLAS → 由 multiprocessing Pool 跨位點平行(避免小矩陣 oversubscription)
import multiprocessing as _mp
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.metrics import silhouette_score
from collections import Counter
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
BIN=f"{WT}/build/bin/inter_sub_mod"
TUMOR="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"; NORMAL="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF="/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"; VCFDIR="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
ASSET=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
os.environ["TMPDIR"]="/big7_disk/liaoyoyo2001/tmp"
MAXOUT=0.20; NPERM=99; MIN_GRP=3; NPROC=24
LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}
def _f(x):
    x=x.strip(); return np.nan if x in("","NA","nan","NaN") else float(x)
def loadm(p):
    rows=open(p).read().strip().split("\n"); ids=[];M=[]
    for ln in rows[1:]:
        q=ln.split(","); ids.append(q[0]); M.append([_f(x) for x in q[1:]])
    return ids,np.array(M)
def peel(s):
    idx=list(range(s.shape[0]))
    while True:
        S=s[np.ix_(idx,idx)];b=(S<0)|np.isnan(S);np.fill_diagonal(b,False)
        if not b.any(): return idx
        idx.remove(idx[int(np.argmax(b.sum(1)))])
        if len(idx)<6: return idx
def cramerv(g,c):
    G=sorted(set(g));C=sorted(set(c))
    if len(G)<2 or len(C)<2: return None,None,None
    t=np.zeros((len(G),len(C)));gi={x:i for i,x in enumerate(G)};ci={x:i for i,x in enumerate(C)}
    for a,b in zip(g,c): t[gi[a],ci[b]]+=1
    t=t[t.sum(1)>0][:,t.sum(0)>0]
    if min(t.shape)<2: return None,None,None
    try: chi2,p,_,e=chi2_contingency(t)
    except: return None,None,None
    nn=t.sum(); return float(np.sqrt(chi2/(nn*(min(t.shape)-1)))),float(p),float(e.min())
# ---- mean-shift 軸: PERMANOVA + PERMDISP (distance-based) ----
def _ssw(D2,lab,a):
    s=0.0
    for k in range(a):
        idx=np.where(lab==k)[0]; ng=len(idx)
        if ng<2: continue
        sub=D2[np.ix_(idx,idx)]; s+=sub[np.triu_indices(ng,1)].sum()/ng
    return s
def permanova(D,lab,rng):
    """distance-based pseudo-F + permutation p。lab=int array(0..a-1)。2 群走向量化二次型。"""
    n=len(lab); a=int(lab.max())+1
    if n<6 or a<2: return None,None
    D2=D**2; SST=0.5*D2.sum()/n
    if a==2:
        m1=(lab==1).astype(float); m0=1.0-m1; n1=m1.sum(); n0=n-n1
        if n1<2 or n0<2: return None,None
        SSW=0.5*(m1@D2@m1)/n1+0.5*(m0@D2@m0)/n0
        if SSW<=0: return None,None
        F=(SST-SSW)/(SSW/(n-2))
        P=rng.permuted(np.broadcast_to(lab,(NPERM,n)).copy(),axis=1).astype(float)
        Q1=((P@D2)*P).sum(1); Pm0=1.0-P; Q0=((Pm0@D2)*Pm0).sum(1)
        nn1=P.sum(1); nn0=n-nn1; ok=(nn1>=2)&(nn0>=2)
        SSWp=0.5*Q1/np.where(nn1>0,nn1,1)+0.5*Q0/np.where(nn0>0,nn0,1)
        Fp=np.where(ok&(SSWp>0),(SST-SSWp)/(SSWp/(n-2)),-np.inf)
        return float(F),(1+int((Fp>=F).sum()))/(NPERM+1)
    SSW=_ssw(D2,lab,a); SSA=SST-SSW
    if SSW<=0: return None,None
    F=(SSA/(a-1))/(SSW/(n-a)); cnt=1
    for _ in range(NPERM):
        lp=rng.permutation(lab); SSWp=_ssw(D2,lp,a)
        if SSWp<=0: continue
        Fp=((SST-SSWp)/(a-1))/(SSWp/(n-a))
        if Fp>=F: cnt+=1
    return float(F),cnt/(NPERM+1)
def _anovaF(z,lab,a):
    gm=z.mean(); ssb=0.0; ssw=0.0
    for k in range(a):
        zi=z[lab==k]
        if len(zi)<1: continue
        ssb+=len(zi)*(zi.mean()-gm)**2; ssw+=((zi-zi.mean())**2).sum()
    n=len(z)
    if ssw<=0 or a<2 or n-a<1: return None
    return (ssb/(a-1))/(ssw/(n-a))
def permdisp(D,lab,rng):
    """PCoA betadisper: 各點到組心距 z 的 ANOVA-F + permutation p (vegan real/imag 軸)。"""
    n=len(lab); a=int(lab.max())+1
    if n<6 or a<2: return None,None
    A=-0.5*(D**2); J=np.eye(n)-np.ones((n,n))/n; G=J@A@J
    try: w,V=np.linalg.eigh(G)
    except: return None,None
    keep=np.abs(w)>1e-8; w=w[keep]; V=V[:,keep]
    if w.size==0: return None,None
    X=V*np.sqrt(np.abs(w)); sgn=np.sign(w)
    z=np.zeros(n)
    for k in range(a):
        m=lab==k
        if m.sum()<1: continue
        c=X[m].mean(0); diff=X[m]-c; ss=((diff**2)*sgn).sum(1); z[m]=np.sqrt(np.abs(ss))
    F=_anovaF(z,lab,a)
    if F is None: return None,None
    if a==2:  # 向量化: z 固定、permute 標籤的 ANOVA-F
        gm=z.mean(); tss=((z-gm)**2).sum(); zsum=z.sum()
        if tss<=0: return float(F),1.0
        P=rng.permuted(np.broadcast_to((lab==1),(NPERM,n)).copy(),axis=1)
        n1=P.sum(1); s1=(P*z).sum(1); n0=n-n1
        ok=(n1>=1)&(n0>=1)
        ssb=np.where(ok,n1*(s1/np.where(n1>0,n1,1)-gm)**2+n0*((zsum-s1)/np.where(n0>0,n0,1)-gm)**2,0.0)
        ssw=tss-ssb
        Fp=np.where(ok&(ssw>0),ssb/(ssw/(n-2)),-np.inf)
        return float(F),(1+int((Fp>=F).sum()))/(NPERM+1)
    cnt=1
    for _ in range(NPERM):
        lp=rng.permutation(lab); Fp=_anovaF(z,lp,a)
        if Fp is not None and Fp>=F: cnt+=1
    return float(F),cnt/(NPERM+1)
def axis_meanshift(sub,kids,reads,rng):
    """3 軸(hp/carrier/allele) PERMANOVA+PERMDISP。回 {axis:{n_grp,F,p,dispF,dispP}}。"""
    out={}
    hp=[LABMAP[reads[r]["hp"]] for r in kids]
    defs={"hp":["A" if l in("1","1-1") else "B" for l in hp],
          "carrier":["G" if l in("1","2") else "C" for l in hp],
          "allele":[reads[r]["alt_support"] for r in kids]}
    for ax,lab in defs.items():
        if ax=="allele": sel=[i for i,v in enumerate(lab) if v in("REF","ALT")]
        else: sel=list(range(len(lab)))
        gl=[lab[i] for i in sel]; cnt=Counter(gl)
        grp=[g for g,c in cnt.items() if c>=MIN_GRP]
        if len(grp)<2: out[ax]=None; continue
        sel=[sel[j] for j in range(len(sel)) if gl[j] in grp]
        if len(sel)<6: out[ax]=None; continue
        gmap={g:i for i,g in enumerate(sorted(grp))}
        li=np.array([gmap[lab[i]] for i in sel]); sm=sub[np.ix_(sel,sel)]
        F,p=permanova(sm,li,rng); dF,dP=permdisp(sm,li,rng)
        out[ax]=None if F is None else {"n_grp":len(grp),"F":round(F,3),"p":round(p,4),
            "dispF":None if dF is None else round(dF,3),"dispP":None if dP is None else round(dP,4)}
    return out
# ---- 離散結構軸: outlier-tolerant 切群 + 對齊 ----
def best_split(sub,Z,n,T):
    best=None
    for k in range(2,min(9,n)+1):
        cl=fcluster(Z,k,'maxclust');sz=Counter(cl)
        valid=[c for c,s in sz.items() if s>=T]
        if len(valid)<2: continue
        nout=sum(s for c,s in sz.items() if s<T)
        if nout/n>MAXOUT: continue
        mask=[i for i in range(n) if cl[i] in valid]
        if len(set(cl[i] for i in mask))<2: continue
        try: sil=silhouette_score(sub[np.ix_(mask,mask)],[cl[i] for i in mask],metric="precomputed")
        except: continue
        if best is None or sil>best[0]: best=(sil,len(valid),sorted([sz[c] for c in valid],reverse=True),nout,list(cl),mask)
    if best is None: return None
    return dict(sil=round(best[0],3),n_valid=best[1],sizes=best[2],n_out=best[3],cut=best[4],mask=best[5])
def aligned(labels,bs):
    cl=bs["cut"]; mask=bs["mask"]; clm=[cl[i] for i in mask]; out=None
    for ax in ("hp","carrier","allele"):
        gl=[labels[ax][i] for i in mask]
        idx=[j for j in range(len(mask)) if (ax!="allele" or gl[j] in("REF","ALT"))]
        if len(idx)<6: continue
        v,p,e=cramerv([gl[j] for j in idx],[clm[j] for j in idx])
        if v is not None and v>=0.3 and p is not None and p<0.05 and e is not None and e>=5:
            if out is None or v>out[1]: out=(ax,round(v,3))
    return out
def process(M,ids,reads,kidx,chrom,pos):
    kids=[ids[i] for i in kidx]
    labset=set(LABMAP[reads[r]["hp"]] for r in kids)
    rec={"chrom":chrom,"pos":pos,"n":len(kids),"lab_n":len(labset),"n_complete":None,"byT":{},"perm":{}}
    if len(kids)<6: return rec
    sub=M[np.ix_(kidx,kidx)]; kp=peel(sub); rec["n_complete"]=len(kp)
    if len(kp)<6: return rec
    kids=[kids[i] for i in kp]; sub=sub[np.ix_(kp,kp)].copy(); np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T)
    n=len(kids); Z=linkage(squareform(sub,checks=False),method="average")
    hp=[LABMAP[reads[r]["hp"]] for r in kids]
    labels={"hp":["A" if l in("1","1-1") else "B" for l in hp],"carrier":["G" if l in("1","2") else "C" for l in hp],"allele":[reads[r]["alt_support"] for r in kids]}
    rng=np.random.default_rng(int(pos)%(2**32) if pos and pos.isdigit() else 12345)
    rec["perm"]=axis_meanshift(sub,kids,reads,rng)
    for T in (3,4,5):
        bs=best_split(sub,Z,n,T)
        if bs:
            al=aligned(labels,bs)
            rec["byT"][T]={"split":True,"sil":bs["sil"],"n_valid":bs["n_valid"],"minority":min(bs["sizes"]),"n_out":bs["n_out"],"aligned":bool(al),"ax":al[0] if al else None}
        else: rec["byT"][T]={"split":False}
    return rec

def _istum(t): return str(t) in("1","true","True")
def _work(arg):
    """worker: 一個位點 matrix → (tumor_rec, merged_rec)。單線程 BLAS, Pool 跨位點平行。"""
    mp,chrom=arg
    rdir=os.path.dirname(os.path.dirname(os.path.dirname(mp))); rdt=os.path.join(rdir,"reads","reads.tsv")
    if not os.path.exists(rdt): return None
    pos=None
    for part in rdir.split("/"):
        if part.startswith(f"{chrom}_") and part.count("_")==1: pos=part.split("_")[1]; break
    try:
        reads={r["read_id"]:r for r in csv.DictReader(open(rdt),delimiter="\t")}
        ids,M=loadm(mp)
        tum=[i for i,rid in enumerate(ids) if rid in reads and _istum(reads[rid]["is_tumor"]) and reads[rid]["hp"] in LABMAP]
        mer=[i for i,rid in enumerate(ids) if rid in reads and reads[rid]["hp"] in LABMAP]
        return (process(M,ids,reads,tum,chrom,pos), process(M,ids,reads,mer,chrom,pos))
    except Exception as e:
        return {"_err":f"{chrom}:{pos}:{type(e).__name__}:{e}"}
def main():
  for d in glob.glob(f"{WT}/output/_df_tmp_*"): shutil.rmtree(d,ignore_errors=True)
  rt_=[]; rm_=[]; t0=time.time(); nerr=0
  pool=_mp.Pool(NPROC)
  for c in range(1,23):
    chrom=f"chr{c}"; vcf=f"{VCFDIR}/filtered_snv_tp_{chrom}.vcf.gz"
    if not os.path.exists(vcf): continue
    out=f"{WT}/output/_df_tmp_{chrom}"; shutil.rmtree(out,ignore_errors=True); os.makedirs(out)
    subprocess.run([BIN,"-t",TUMOR,"-n",NORMAL,"-r",REF,"-v",vcf,"-w","5000","-j","16",
        "--distance-metric","BERNOULLI","--nan-distance-strategy","SKIP","-o",out],
        stdout=open(f"{out}/run.log","w"),stderr=subprocess.STDOUT)
    n0=len(rt_)
    paths=[(p,chrom) for p in glob.glob(f"{out}/**/distance/BERNOULLI/matrix.csv",recursive=True)]
    for res in pool.imap_unordered(_work,paths,chunksize=8):
        if res is None: continue
        if isinstance(res,dict) and "_err" in res:
            nerr+=1
            if nerr<=5: print("  ERR",res["_err"],flush=True)
            continue
        rt_.append(res[0]); rm_.append(res[1])
    shutil.rmtree(out,ignore_errors=True)
    json.dump(rt_,open(f"{ASSET}/decisionflow_records_tumor.json","w")); json.dump(rm_,open(f"{ASSET}/decisionflow_records_merged.json","w"))
    print(f"[{chrom}] +{len(rt_)-n0} tumor={len(rt_)} merged={len(rm_)} err={nerr} elapsed={int(time.time()-t0)}s",flush=True)
  pool.close(); pool.join()
  json.dump(rt_,open(f"{ASSET}/decisionflow_records_tumor.json","w")); json.dump(rm_,open(f"{ASSET}/decisionflow_records_merged.json","w"))
  print(f"DONE tumor={len(rt_)} merged={len(rm_)} err={nerr} elapsed={int(time.time()-t0)}s")
if __name__=="__main__":
    main()
