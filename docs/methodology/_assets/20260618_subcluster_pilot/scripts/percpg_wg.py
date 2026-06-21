#!/usr/bin/env python3
"""全基因組 per-CpG-by-label (modkit-equiv marginal) 顯著 CpG 數 — 每位點三軸 (hp/carrier/allele)。
逐位點逐 CpG Fisher (meth/unmeth × axis-group) + BH-FDR within locus → n_sig per axis。
之後與 decisionflow_records (joint PERMANOVA) merge 做 Venn:
 ISM-only = joint-sig 但 per-CpG 貧乏 / modkit-only = per-CpG-sig 但 joint 不顯著 / both / neither。
disk-safe 逐 chr + multiprocessing。輸出 percpg_records.json。"""
import os, csv, glob, json, subprocess, shutil, time
for _v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[_v]="1"
import multiprocessing as _mp
import numpy as np
from scipy.stats import fisher_exact
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
BIN=f"{WT}/build/bin/inter_sub_mod"
TUMOR="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"; NORMAL="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF="/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"; VCFDIR="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
ASSET=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
os.environ["TMPDIR"]="/big7_disk/liaoyoyo2001/tmp"
NPROC=24; MIN_GRP=5; LABMAP={"1","2","1-1","2-1","HP1","HP2"}
def bh(ps):
    p=np.asarray(ps,float); m=len(p)
    if m==0: return p
    o=np.argsort(p); r=p[o]*m/np.arange(1,m+1); q=np.minimum.accumulate(r[::-1])[::-1]
    out=np.empty(m); out[o]=np.clip(q,0,1); return out
def axlab(r,axis):
    h=r['hp']; a=r['alt_support']
    if axis=='hp': return 'A' if h in('1','1-1','HP1') else 'B'
    if axis=='carrier': return 'g' if h in('1','2','HP1','HP2') else 's'
    if axis=='allele': return a if a in('REF','ALT') else None
def percpg_axis(M,ids,reads,tum_idx,cpgn):
    """回 {axis:{n_test,n_sig}} — 逐 CpG Fisher + BH within locus。M=read×CpG。"""
    out={}
    for axis in ('hp','carrier','allele'):
        ps=[]
        for ci in range(cpgn):
            g={}
            for i in tum_idx:
                l=axlab(reads[ids[i]],axis)
                if l is None: continue
                v=M[i,ci]
                if np.isnan(v): continue
                g.setdefault(l,[]).append(1 if v>=0.5 else 0)
            ks=sorted(g)
            if len(ks)<2 or any(len(g[k])<MIN_GRP for k in ks): continue
            a0,a1=g[ks[0]],g[ks[1]]; m0,m1=sum(a0),sum(a1)
            try: _,p=fisher_exact([[m1,len(a1)-m1],[m0,len(a0)-m0]])
            except: continue
            ps.append(p)
        if not ps: out[axis]={'n_test':0,'n_sig':0}; continue
        q=bh(ps); out[axis]={'n_test':len(ps),'n_sig':int((q<0.05).sum())}
    return out
def loadmeth(p):
    lines=open(p).read().strip().split('\n')
    cpgs=lines[0].split(',')[1:]; ids=[]; M=[]
    for ln in lines[1:]:
        q=ln.split(','); ids.append(q[0]); M.append([np.nan if v in('','NA','nan','NaN') else float(v) for v in q[1:]])
    return ids,np.array(M),len(cpgs)
def _work(arg):
    mp,chrom=arg
    rdir=mp.rsplit('/methylation/',1)[0]; rdt=os.path.join(rdir,'reads','reads.tsv')
    if not os.path.exists(rdt): return None
    pos=None
    for part in rdir.split('/'):
        if part.startswith(f"{chrom}_") and part.count('_')==1: pos=part.split('_')[1]; break
    try:
        reads={r['read_id']:r for r in csv.DictReader(open(rdt),delimiter='\t')}
        ids,M,cpgn=loadmeth(mp)
        def istum(t): return str(t) in('1','true','True')
        tum=[i for i,rid in enumerate(ids) if rid in reads and istum(reads[rid]['is_tumor']) and reads[rid]['hp'] in LABMAP]
        if len(tum)<6 or cpgn<1: return {'chrom':chrom,'pos':pos,'n_tum':len(tum),'n_cpg':cpgn,'axes':None}
        ax=percpg_axis(M,ids,reads,tum,cpgn)
        return {'chrom':chrom,'pos':pos,'n_tum':len(tum),'n_cpg':cpgn,'axes':ax}
    except Exception as e:
        return {'_err':f"{chrom}:{pos}:{type(e).__name__}:{e}"}
def main():
    for d in glob.glob(f"{WT}/output/_pcpg_tmp_*"): shutil.rmtree(d,ignore_errors=True)
    rec=[]; t0=time.time(); nerr=0; pool=_mp.Pool(NPROC)
    for c in range(1,23):
        chrom=f"chr{c}"; vcf=f"{VCFDIR}/filtered_snv_tp_{chrom}.vcf.gz"
        if not os.path.exists(vcf): continue
        out=f"{WT}/output/_pcpg_tmp_{chrom}"; shutil.rmtree(out,ignore_errors=True); os.makedirs(out)
        subprocess.run([BIN,"-t",TUMOR,"-n",NORMAL,"-r",REF,"-v",vcf,"-w","5000","-j","16",
            "--distance-metric","BERNOULLI","--nan-distance-strategy","SKIP","-o",out],
            stdout=open(f"{out}/run.log","w"),stderr=subprocess.STDOUT)
        n0=len(rec)
        paths=[(p,chrom) for p in glob.glob(f"{out}/**/methylation/methylation.csv",recursive=True)]
        for res in pool.imap_unordered(_work,paths,chunksize=8):
            if res is None: continue
            if '_err' in res:
                nerr+=1
                if nerr<=5: print("  ERR",res['_err'],flush=True)
                continue
            rec.append(res)
        shutil.rmtree(out,ignore_errors=True)
        json.dump(rec,open(f"{ASSET}/percpg_records.json","w"))
        print(f"[{chrom}] +{len(rec)-n0} total={len(rec)} err={nerr} elapsed={int(time.time()-t0)}s",flush=True)
    pool.close(); pool.join()
    json.dump(rec,open(f"{ASSET}/percpg_records.json","w"))
    print(f"DONE total={len(rec)} err={nerr} elapsed={int(time.time()-t0)}s")
if __name__=="__main__": main()
