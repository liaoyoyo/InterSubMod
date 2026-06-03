#!/usr/bin/env python3
"""本資料直接算 aDMR（allele-specific DMR）+ LOH 富集，對照文獻 79%。

aDMR 定義（per ±win 窗，需 HP1/HP2 anchor read 充足）:
  對窗內每個 CpG，比 HP1 vs HP2 read 的 mean 甲基 β → |Δβ|。
  窗級 aDMR = 該窗有 >=min_sig_cpg 個 CpG 達 |Δβ|>=delta_thr (Mann-Whitney p<0.05)。
  → 即「該窗 HP1/HP2 甲基顯著不同」= allele-specific DMR。

富集對照: aDMR 窗落在 SEQC2 LOH 區的比例 vs 全窗落 LOH 區的比例（背景率）。
  文獻: 79% aDMR 落 CNV/LOH 區 + aDMR 數與 LOH 比例正相關。
  本資料報: frac(aDMR in LOH) + odds ratio(aDMR vs non-aDMR 落 LOH) + 背景 LOH 覆蓋率。

⚠ confound: read 數 → 易顯著。報告同時記 n_read，並可分層。
唯讀。輸出 admr_{chrom}.json。
"""
import sys, json, argparse, bisect
import numpy as np
import pysam
from scipy.stats import mannwhitneyu

BAM="/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GVCF="/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
LOHBED="/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
AD="/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets"

# SEQC2 所有 CNV/LOH 區（gain/loss/loh 都算"CNV/LOH 區"，對齊文獻 "CNV/LOH"）
cnv_loh={}
with open(LOHBED) as f:
    for line in f:
        p=line.rstrip().split("\t")
        if len(p)>=4 and p[3] in ("gain","loss","loh"):
            cnv_loh.setdefault(p[0],[]).append((int(p[1]),int(p[2])))
for c in cnv_loh: cnv_loh[c].sort()
# 也單獨 LOH
loh_only={}
with open(LOHBED) as f:
    for line in f:
        p=line.rstrip().split("\t")
        if len(p)>=4 and p[3]=="loh": loh_only.setdefault(p[0],[]).append((int(p[1]),int(p[2])))
for c in loh_only: loh_only[c].sort()
def in_seg(segs,pos):
    for s,e in segs:
        if s<=pos<=e: return True
    return False

def get_hp(r):
    try: return str(r.get_tag("HP"))
    except KeyError: return "."

def window_admr(bam,chrom,ws,we,min_anchor=8,delta_thr=0.25,min_sig_cpg=2):
    """回傳 (is_admr, n_sig_cpg, max_delta, n_hp1, n_hp2) 或 None。"""
    hp1=[];hp2=[];allpos=set();seen=set()
    for r in bam.fetch(chrom,ws,we):
        if r.is_secondary or r.is_supplementary or r.is_unmapped:continue
        if r.query_name in seen:continue
        hp=get_hp(r)
        if hp not in ("1","2"):continue
        mods=r.modified_bases or {}
        mc=None
        for k,c in mods.items():
            if k[0] in ("C",b"C") and k[2] in ("m",27551):mc=c;break
        if not mc:continue
        q2r={q:rr for q,rr in r.get_aligned_pairs(matches_only=True)}
        meth={}
        for qp,qual in mc:
            rp=q2r.get(qp)
            if rp is None or rp<ws or rp>we:continue
            meth[rp]=qual/255.0
        if len(meth)<3:continue
        seen.add(r.query_name);allpos.update(meth.keys())
        (hp1 if hp=="1" else hp2).append(meth)
    if len(hp1)<min_anchor or len(hp2)<min_anchor:return None
    positions=sorted(allpos)
    n_sig=0;maxd=0
    for p in positions:
        v1=[m[p] for m in hp1 if p in m];v2=[m[p] for m in hp2 if p in m]
        if len(v1)<4 or len(v2)<4:continue
        d=abs(np.mean(v1)-np.mean(v2))
        if d>maxd:maxd=d
        if d>=delta_thr:
            try:
                _,pp=mannwhitneyu(v1,v2,alternative="two-sided")
                if pp<0.05:n_sig+=1
            except Exception:pass
    return (n_sig>=min_sig_cpg, n_sig, round(float(maxd),3), len(hp1), len(hp2))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--chrom",required=True)
    ap.add_argument("--n-windows",type=int,default=200)
    ap.add_argument("--win",type=int,default=2000)
    a=ap.parse_args()
    chrom=a.chrom
    rng=np.random.RandomState(20260603)
    vcf=pysam.VariantFile(GVCF)
    hets=[]
    for rec in vcf.fetch(chrom):
        s0=rec.samples[0];gt=s0.get("GT")
        if gt and len(gt)==2 and gt[0]!=gt[1] and None not in gt:hets.append(rec.pos)
    if len(hets)<20:
        json.dump({"chrom":chrom,"n":0},open(f"{AD}/admr_{chrom}.json","w"));print(f"[admr] {chrom}: few hets");return
    # 隨機抽窗（het 為中心，確保有 anchor）
    centers=[hets[i] for i in rng.choice(len(hets),min(a.n_windows*3,len(hets)),replace=False)]
    bam=pysam.AlignmentFile(BAM,"rb")
    wins=[]
    for c in centers:
        if len(wins)>=a.n_windows:break
        r=window_admr(bam,chrom,c-a.win,c+a.win)
        if r is None:continue
        is_admr,nsig,maxd,n1,n2=r
        wins.append({"center":int(c),"is_admr":is_admr,"n_sig_cpg":nsig,"max_delta":maxd,
                     "n_hp1":n1,"n_hp2":n2,
                     "in_loh":in_seg(loh_only.get(chrom,[]),c),
                     "in_cnvloh":in_seg(cnv_loh.get(chrom,[]),c)})
    bam.close()
    admr=[w for w in wins if w["is_admr"]]
    nonadmr=[w for w in wins if not w["is_admr"]]
    def frac_loh(grp,key): return round(sum(1 for w in grp if w[key])/len(grp),3) if grp else None
    out={"chrom":chrom,"n_windows":len(wins),"n_admr":len(admr),
         "frac_admr":round(len(admr)/len(wins),3) if wins else None,
         "admr_frac_in_loh":frac_loh(admr,"in_loh"),
         "admr_frac_in_cnvloh":frac_loh(admr,"in_cnvloh"),
         "nonadmr_frac_in_loh":frac_loh(nonadmr,"in_loh"),
         "nonadmr_frac_in_cnvloh":frac_loh(nonadmr,"in_cnvloh"),
         "bg_frac_in_cnvloh":frac_loh(wins,"in_cnvloh"),
         "windows":wins}
    json.dump(out,open(f"{AD}/admr_{chrom}.json","w"),ensure_ascii=False,indent=2)
    print(f"[admr] {chrom}: {len(wins)}窗 aDMR={len(admr)}({out['frac_admr']}) "
          f"aDMR落CNV/LOH={out['admr_frac_in_cnvloh']} vs 背景={out['bg_frac_in_cnvloh']} "
          f"(非aDMR落={out['nonadmr_frac_in_cnvloh']})")

if __name__=="__main__":main()
