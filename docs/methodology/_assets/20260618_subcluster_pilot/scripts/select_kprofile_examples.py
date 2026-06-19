#!/usr/bin/env python3
"""Phase 3 選代表: 各 k-profile class 抽 nk>=3 代表 → kprofile_examples.json + mini-VCF。
confident-unique / multi-resolution(優先不同軸) / ambiguous-near-tie / single-k-forced(示意 k=2 only)。"""
import json, subprocess
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
VCFDIR="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
loci=json.load(open(f"{A}/kprofile_loci_tumor.json"))
by={}
for l in loci: by.setdefault(l["cls"],[]).append(l)

def pick(cls,n,key,need_nk3=True,prefer=None):
    c=[l for l in by.get(cls,[]) if (l["nk"]>=3 if need_nk3 else True)]
    if prefer: c=sorted(c,key=prefer)
    else: c=sorted(c,key=key)
    return c[:n]

sel=[]
# confident-unique: 最大 margin
sel+= [dict(l,group="confident-unique") for l in pick("confident-unique",8,
        key=lambda l:-(l["align_margin"] or l["sil_margin"] or 0))]
# multi-resolution: 優先 meaningful_ks 跨不同軸 (mk_axes 軸種類數多)
def naxes(l): return len({v[0] for v in l["mk_axes"].values()})
sel+= [dict(l,group="multi-resolution") for l in pick("multi-resolution",10,
        key=lambda l:(-naxes(l),-len(l["meaningful_ks"])))]
# ambiguous-near-tie: margin 最小 (近平手最典型)
sel+= [dict(l,group="ambiguous-near-tie") for l in pick("ambiguous-near-tie",8,
        key=lambda l:(l["sil_margin"] if l["sil_margin"] is not None else 9))]
# single-k-forced: 示意 (nk=1, k=2 only) — 取 aligned 的
sel+= [dict(l,group="single-k-forced") for l in pick("single-k-forced",4,
        key=lambda l:-l["n"],need_nk3=False) if l["sub"]=="aligned"][:4]

json.dump({"n":len(sel),"items":sel},open(f"{A}/kprofile_examples.json","w"),indent=1)
from collections import Counter
print("選取:",dict(Counter(s["group"] for s in sel)),"共",len(sel))

# build mini-VCF
bychr={}
for s in sel: bychr.setdefault(s["chrom"],set()).add(s["pos"])
hdr=None; body=set()
for chrom,poss in bychr.items():
    raw=subprocess.run(["zcat",f"{VCFDIR}/filtered_snv_tp_{chrom}.vcf.gz"],capture_output=True,text=True).stdout.splitlines()
    if hdr is None: hdr=[ln for ln in raw if ln.startswith("#")]
    for ln in raw:
        if not ln.startswith("#") and ln.split("\t")[1] in poss: body.add(ln)
def sk(ln): f=ln.split("\t"); return (int(f[0].replace("chr","").replace("X","23").replace("Y","24")),int(f[1]))
mini=f"{A}/kprofile_examples_mini.vcf"
open(mini,"w").write("\n".join(hdr)+"\n"+"\n".join(sorted(body,key=sk))+"\n")
subprocess.run(["bgzip","-f",mini]); subprocess.run(["tabix","-f","-p","vcf",mini+".gz"])
print(f"mini-VCF records={len(body)} → {mini}.gz")
