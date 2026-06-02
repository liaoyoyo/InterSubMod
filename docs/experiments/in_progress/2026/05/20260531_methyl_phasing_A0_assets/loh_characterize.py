#!/usr/bin/env python3
"""LOH 純-LOH 異常定性：對每個已分析的 LOH 區量化判別指標，回答
『分得開的 LOH 為何分得開』。對齊用戶問：補上研究 + 可能性 + 驗證觀察。

對每個 v2 LOH region 量化：
  1. het_vaf_median: 該 ±2kb 區 germline het SNP 的 VAF 中位數
     （真純 LOH 該無 het 或 VAF 極偏；VAF~0.5 = 功能上雜合，SEQC2 LOH 標記與 read 層不符）
  2. n_het: het SNP 數（純 LOH 該少）
  3. loh_boundary_dist: 到最近 SEQC2 LOH segment 邊界的距離（小=邊界效應）
  4. hp_balance: min(HP1,HP2)/max(HP1,HP2)（純 LOH 該極不平衡）
  5. 已有: anchor_auc_raw, separable, gmm_bic

輸出：分得開 vs 分不開的這些指標對比 → 哪個指標最能解釋「為何分得開」。
唯讀。輸出 loh_characterize.json。
"""
import json, glob, bisect
import numpy as np
import pysam

GVCF="/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
LOHBED="/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
AD="/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets"

# LOH segments per chr
loh_seg={}
with open(LOHBED) as f:
    for line in f:
        p=line.rstrip().split("\t")
        if len(p)>=4 and p[3]=="loh":
            loh_seg.setdefault(p[0],[]).append((int(p[1]),int(p[2])))
for c in loh_seg: loh_seg[c].sort()

def boundary_dist(chrom,pos):
    segs=loh_seg.get(chrom,[])
    best=None
    for s,e in segs:
        if s<=pos<=e:
            d=min(pos-s,e-pos)
            best=d if best is None else min(best,d)
    return best

vcf=pysam.VariantFile(GVCF)
def het_stats(chrom,pos,win=2000):
    vafs=[]; n=0
    try:
        for rec in vcf.fetch(chrom,max(0,pos-win),pos+win):
            s0=rec.samples[0];gt=s0.get("GT");af=s0.get("AF")
            if gt and len(gt)==2 and gt[0]!=gt[1] and None not in gt:
                n+=1
                if af is not None:
                    vafs.append(af[0] if isinstance(af,tuple) else af)
    except Exception: pass
    return (round(float(np.median(vafs)),3) if vafs else None), n

# 載入 v2 LOH regions
rows=[]
for fp in glob.glob(AD+"/seqc2_cn_methyl_v2_chr*.json"):
    rows.extend(json.load(open(fp))["regions"])
loh=[r for r in rows if r["seqc2_status"]=="loh"]
# 去重疊
loh.sort(key=lambda r:(r["chrom"],r["pos"]))
ded=[];last=(None,-1e18)
for r in loh:
    if r["chrom"]!=last[0] or r["pos"]-last[1]>=4000:
        ded.append(r);last=(r["chrom"],r["pos"])

def sep(r): return r["anchor_auc_raw"] and r["shuffle_p95_raw"] and r["anchor_auc_raw"]-r["shuffle_p95_raw"]>0.05

enriched=[]
for r in ded:
    vaf,nhet=het_stats(r["chrom"],r["pos"])
    bd=boundary_dist(r["chrom"],r["pos"])
    hpb=round(min(r["n_hp1"],r["n_hp2"])/max(r["n_hp1"],r["n_hp2"]),3) if max(r["n_hp1"],r["n_hp2"])>0 else None
    enriched.append({"region":f"{r['chrom']}:{r['pos']}","separable":sep(r),
                     "auc":r["anchor_auc_raw"],"het_vaf_median":vaf,"n_het":nhet,
                     "loh_boundary_dist":bd,"hp_balance":hpb,"gmm_bic":r["gmm_bic_diff"],
                     "n_hp1":r["n_hp1"],"n_hp2":r["n_hp2"]})

def summ(grp,key):
    vals=[x[key] for x in grp if x[key] is not None]
    return round(float(np.median(vals)),3) if vals else None
sepg=[x for x in enriched if x["separable"]]
notg=[x for x in enriched if not x["separable"]]
out={
  "n_loh_dedup":len(enriched),"n_separable":len(sepg),"n_notsep":len(notg),
  "separable_vs_notsep":{
    "het_vaf_median":{"sep":summ(sepg,"het_vaf_median"),"notsep":summ(notg,"het_vaf_median")},
    "n_het":{"sep":summ(sepg,"n_het"),"notsep":summ(notg,"n_het")},
    "loh_boundary_dist":{"sep":summ(sepg,"loh_boundary_dist"),"notsep":summ(notg,"loh_boundary_dist")},
    "hp_balance":{"sep":summ(sepg,"hp_balance"),"notsep":summ(notg,"hp_balance")},
    "gmm_bic":{"sep":summ(sepg,"gmm_bic"),"notsep":summ(notg,"gmm_bic")},
  },
  "regions":enriched,
}
json.dump(out,open(AD+"/loh_characterize.json","w"),ensure_ascii=False,indent=2)
print(f"LOH 定性: {len(enriched)} 區 (sep={len(sepg)} notsep={len(notg)})")
print("指標對比 (sep vs notsep median):")
for k,v in out["separable_vs_notsep"].items():
    print(f"  {k}: sep={v['sep']} notsep={v['notsep']}")
