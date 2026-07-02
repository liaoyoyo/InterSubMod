#!/usr/bin/env python3
"""選代表位點驗證「切不出有訊號」: 切不出·有訊號(HP p<=.001+|germΔβ|>.15) / 切不出·真null / 切得出對照。
→ kprofile... → cantsplit_examples.json + mini-VCF。供 dual-panel 熱圖 + 距離分佈。"""
import json, csv, subprocess
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
SIG=f"{WT}/output/_wg_bdcprime_verify/significance_summary.csv"
VCFDIR="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
recs={f"{r['chrom']}_{r['pos']}":r for r in json.load(open(f"{A}/records_wg2.json"))}
def lab_n(r): return len([l for l,c in r["all"]["labels"].items() if c>0])
def ff(v):
    try: return float(v)
    except: return None
def truthy(v): return str(v).strip().lower() in("true","1","yes")
sig={}
for row in csv.DictReader(open(SIG)):
    sig[f"{row['Chr']}_{row['Pos']}"]=row
cant_sig=[]; cant_null=[]; split=[]
for key,r in recs.items():
    bk=r["all"]["best_k"]; n=r["all"]["n"]; s=sig.get(key)
    if not s: continue
    p=ff(s.get("LabelHPPermanovaP")); valid=truthy(s.get("LabelHPPermanovaValid"))
    gm=abs(ff(s.get("GermlineAsmDbeta")) or 0)
    if bk is None and n>=20 and lab_n(r)>=2 and valid:
        if p is not None and p<=0.0011 and gm>0.15: cant_sig.append((key,r,gm,p))
        elif (p is None or p>0.3) and gm<0.05: cant_null.append((key,r,gm,p))
    elif bk is not None and bk>=2 and n>=20: split.append((key,r,gm,p))
cant_sig.sort(key=lambda x:-x[2]); cant_null.sort(key=lambda x:x[2]); split.sort(key=lambda x:-x[2])
sel=([dict(key=k,group="切不出·有訊號",germ=round(g,3),hpP=p) for k,r,g,p in cant_sig[:10]]
   + [dict(key=k,group="切不出·真null",germ=round(g,3),hpP=p) for k,r,g,p in cant_null[:5]]
   + [dict(key=k,group="切得出·對照",germ=round(g,3),hpP=p) for k,r,g,p in split[:5]])
for s in sel: s["chrom"],s["pos"]=s["key"].rsplit("_",1)
json.dump({"n":len(sel),"items":sel},open(f"{A}/cantsplit_examples.json","w"),indent=1)
from collections import Counter
print("選取:",dict(Counter(s["group"] for s in sel)),"共",len(sel))
# mini-VCF
bychr={}
for s in sel: bychr.setdefault(s["chrom"],set()).add(s["pos"])
hdr=None; body=set()
for chrom,poss in bychr.items():
    raw=subprocess.run(["zcat",f"{VCFDIR}/filtered_snv_tp_{chrom}.vcf.gz"],capture_output=True,text=True).stdout.splitlines()
    if hdr is None: hdr=[l for l in raw if l.startswith("#")]
    for l in raw:
        if not l.startswith("#") and l.split("\t")[1] in poss: body.add(l)
def sk(l): f=l.split("\t"); return (int(f[0].replace("chr","").replace("X","23").replace("Y","24")),int(f[1]))
mini=f"{A}/cantsplit_examples_mini.vcf"
open(mini,"w").write("\n".join(hdr)+"\n"+"\n".join(sorted(body,key=sk))+"\n")
subprocess.run(["bgzip","-f",mini]); subprocess.run(["tabix","-f","-p","vcf",mini+".gz"])
print(f"mini-VCF records={len(body)}")
