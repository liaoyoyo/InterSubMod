#!/usr/bin/env python3
"""抽出全部 FDR-顯著「silhouette-k 不是最對齊 k」位點(三軸 union) → ksweep_fdr_loci.json
+ 建 mini-VCF(extract 自 per-chr filtered_snv_tp) 供重跑 binary 拿甲基矩陣畫熱圖。
重現 ksweep_analyze.py 的 gate+Bonferroni+BH-FDR 邏輯, dump 被 reject 的位點全清單。"""
import json, subprocess, os
ASSET="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
VCFDIR="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
DELTA=0.10; FDRq=0.05
recs=json.load(open(f"{ASSET}/ksweep_wg_records.json"))

def bh_reject(pvals,q):
    m=len(pvals)
    if m==0: return [False]*0
    order=sorted(range(m),key=lambda i:pvals[i]); cut=-1
    for rank,i in enumerate(order,1):
        if pvals[i]<=q*rank/m: cut=rank
    rej=[False]*m
    for rank,i in enumerate(order,1):
        if rank<=cut: rej[i]=True
    return rej

def per_axis(axis):
    Vk,pk,ek=f"V_{axis}",f"p_{axis}",f"e_{axis}"
    rows=[]; pvals=[]
    for r in recs:
        valid=[d for d in r["per_k"] if d[Vk] is not None]
        if not valid: continue
        sil_opt=max(r["per_k"],key=lambda d:d["sil"]); V_sil=sil_opt[Vk]
        if V_sil is None: continue
        K=len(valid)
        imp=[d for d in valid if d["k"]!=sil_opt["k"] and (d[Vk]-V_sil)>DELTA and d[ek] is not None and d[ek]>=5]
        if imp:
            best=min(imp,key=lambda d:d[pk]); pb=min(1.0,best[pk]*K)
            rows.append({"chrom":r["chrom"],"pos":r["pos"],"n":r["n"],"axis":axis,
                "sil_k":sil_opt["k"],"V_sil":round(V_sil,3),"align_k":best["k"],
                "V_align":round(best[Vk],3),"dV":round(best[Vk]-V_sil,3),"p_bonf":pb})
        else:
            rows.append(None)  # placeholder, p=1
        pvals.append(rows[-1]["p_bonf"] if rows[-1] else 1.0)
    rej=bh_reject(pvals,FDRq)
    return [rows[i] for i in range(len(rows)) if rows[i] and rej[i]]

fdr={ax:per_axis(ax) for ax in ("hp","carrier","allele")}
# union by (chrom,pos), 記錄該位點在哪些軸 sig + 取 dV 最大的軸當 primary
uni={}
for ax,rows in fdr.items():
    for r in rows:
        key=(r["chrom"],r["pos"])
        if key not in uni: uni[key]={"chrom":r["chrom"],"pos":r["pos"],"n":r["n"],"axes":{},"primary":None}
        uni[key]["axes"][ax]={"sil_k":r["sil_k"],"V_sil":r["V_sil"],"align_k":r["align_k"],"V_align":r["V_align"],"dV":r["dV"],"p_bonf":r["p_bonf"]}
for v in uni.values():
    v["primary"]=max(v["axes"].items(),key=lambda kv:kv[1]["dV"])[0]
loci=sorted(uni.values(),key=lambda v:-max(a["dV"] for a in v["axes"].values()))
out={"n_hp":len(fdr["hp"]),"n_carrier":len(fdr["carrier"]),"n_allele":len(fdr["allele"]),
     "n_union":len(loci),"loci":loci}
json.dump(out,open(f"{ASSET}/ksweep_fdr_loci.json","w"),indent=1)
print(f"FDR-sig: HP={out['n_hp']} CARRIER={out['n_carrier']} ALLELE={out['n_allele']} union={out['n_union']}")

# --- build mini-VCF (union loci) ---
bychr={}
for v in loci: bychr.setdefault(v["chrom"],set()).add(v["pos"])
hdr=None; body=[]
for chrom,poss in bychr.items():
    vcf=f"{VCFDIR}/filtered_snv_tp_{chrom}.vcf.gz"
    raw=subprocess.run(["zcat",vcf],capture_output=True,text=True).stdout.splitlines()
    for ln in raw:
        if ln.startswith("#"):
            if hdr is None and ln.startswith("##") or ln.startswith("#CHROM"): pass
            continue
        f=ln.split("\t")
        if f[1] in poss: body.append(ln)
    if hdr is None:
        hdr=[ln for ln in raw if ln.startswith("#")]
mini=f"{ASSET}/fdr_loci_mini.vcf"
with open(mini,"w") as fo:
    fo.write("\n".join(hdr)+"\n")
    # sort body by chrom,pos
    def sk(ln): f=ln.split("\t"); return (int(f[0].replace("chr","").replace("X","23").replace("Y","24")),int(f[1]))
    for ln in sorted(set(body),key=sk): fo.write(ln+"\n")
subprocess.run(["bgzip","-f",mini]); subprocess.run(["tabix","-f","-p","vcf",mini+".gz"])
print(f"mini-VCF: {mini}.gz  records={len(set(body))}")
