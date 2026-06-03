#!/usr/bin/env python3
"""全 LOH 區 tumor VAF 系統分布 + 甲基救援×LOH 型別交叉（V1+V4）。
回答: SEQC2 LOH 區裡, 多少 tumor 仍雜合(chr15型,無真cnLOH) vs 真偏移(chr8型);
      兩型的甲基救援/分離能力是否不同。

對每個 LOH het 位點: normal VAF (germline VCF) + tumor VAF (BAM pileup) → 分型:
  balanced (tumor |VAF-0.5|<0.15) = 無真 cnLOH / 雜合保留
  imbalanced (0.15<=|VAF-0.5|<0.35) = 部分 allelic imbalance
  homozygous (|VAF-0.5|>=0.35, tumor VAF→0/1) = 真 cnLOH/LOH
唯讀。輸出 loh_tumor_vaf_{chrom}.json。
"""
import sys, json, argparse, bisect
import numpy as np
import pysam

BAM="/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GVCF="/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
LOHBED="/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
AD="/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets"

loh_seg={}
with open(LOHBED) as f:
    for line in f:
        p=line.rstrip().split("\t")
        if len(p)>=4 and p[3]=="loh": loh_seg.setdefault(p[0],[]).append((int(p[1]),int(p[2])))
for c in loh_seg: loh_seg[c].sort()
def in_loh(c,pos):
    for s,e in loh_seg.get(c,[]):
        if s<=pos<=e: return True
    return False

def tumor_vaf(bam,chrom,pos,ref,alt):
    cnt={"r":0,"a":0}
    try:
        for pc in bam.pileup(chrom,pos-1,pos,truncate=True,min_base_quality=10,stepper="samtools"):
            for r in pc.pileups:
                if r.is_del or r.is_refskip or r.query_position is None: continue
                b=r.alignment.query_sequence[r.query_position].upper()
                if b==ref: cnt["r"]+=1
                elif b==alt: cnt["a"]+=1
    except Exception: return None,0
    tot=cnt["r"]+cnt["a"]
    return (cnt["a"]/tot if tot else None), tot

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--chrom",required=True)
    ap.add_argument("--max-sites",type=int,default=400,help="每染色體最多取樣 LOH het 位點")
    ap.add_argument("--min-depth",type=int,default=15)
    a=ap.parse_args()
    chrom=a.chrom
    rng=np.random.RandomState(20260603)
    vcf=pysam.VariantFile(GVCF)
    # 收 LOH 區 het 位點
    hets=[]
    for rec in vcf.fetch(chrom):
        s0=rec.samples[0];gt=s0.get("GT");af=s0.get("AF")
        if not(gt and len(gt)==2 and gt[0]!=gt[1] and None not in gt): continue
        if len(rec.ref)!=1 or len(rec.alts[0])!=1: continue
        if not in_loh(chrom,rec.pos): continue
        nvaf=af[0] if isinstance(af,tuple) else af
        hets.append((rec.pos,rec.ref,rec.alts[0],nvaf))
    if not hets:
        json.dump({"chrom":chrom,"n":0},open(f"{AD}/loh_tumor_vaf_{chrom}.json","w"))
        print(f"[lohvaf] {chrom}: no LOH het");return
    if len(hets)>a.max_sites:
        idx=rng.choice(len(hets),a.max_sites,replace=False); hets=[hets[i] for i in idx]
    bam=pysam.AlignmentFile(BAM,"rb")
    rows=[]
    for pos,ref,alt,nvaf in hets:
        tvaf,tdepth=tumor_vaf(bam,chrom,pos,ref,alt)
        if tvaf is None or tdepth<a.min_depth: continue
        dev=abs(tvaf-0.5)
        typ="balanced" if dev<0.15 else ("imbalanced" if dev<0.35 else "homozygous")
        rows.append({"pos":pos,"normal_vaf":round(nvaf,3) if nvaf else None,
                     "tumor_vaf":round(tvaf,3),"tumor_depth":tdepth,"dev":round(dev,3),"type":typ})
    bam.close()
    from collections import Counter
    tc=Counter(r["type"] for r in rows)
    tvafs=[r["tumor_vaf"] for r in rows]
    out={"chrom":chrom,"n_loh_het":len(rows),
         "tumor_vaf_median":round(float(np.median(tvafs)),3) if tvafs else None,
         "type_counts":dict(tc),
         "frac_balanced":round(tc.get("balanced",0)/len(rows),3) if rows else None,
         "frac_imbalanced":round(tc.get("imbalanced",0)/len(rows),3) if rows else None,
         "frac_homozygous":round(tc.get("homozygous",0)/len(rows),3) if rows else None,
         "sites":rows}
    json.dump(out,open(f"{AD}/loh_tumor_vaf_{chrom}.json","w"),ensure_ascii=False,indent=2)
    print(f"[lohvaf] {chrom}: {len(rows)} LOH het, tumor_vaf_med={out['tumor_vaf_median']}, "
          f"bal/imb/hom={out['frac_balanced']}/{out['frac_imbalanced']}/{out['frac_homozygous']}")

if __name__=="__main__":main()
