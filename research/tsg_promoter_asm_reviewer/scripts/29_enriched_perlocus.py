#!/usr/bin/env python3
"""
Enriched per-locus dataset for Q1/Q2/Q3 (paired HCC1395, has HP tags).
Per locus compute:
  - blind k=2 cluster (observed-only distance)
  - ARI(blind, HP-tag)  vs  ARI(blind, ALT/REF)   [Q1: which labeling does clustering align with?]
  - ARI(HP-tag, ALT/REF)                            [how much HP and allele agree]
  - confounds: LOH status, coverage(n), VAF, mean MAPQ, mean NM/len (error proxy) per allele group
Output: enriched JSON for workflow agents (Q2 attribution, Q3 map).
"""
import numpy as np, random, subprocess, json
import pysam
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score
from collections import Counter
random.seed(42); np.random.seed(42)
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
PILE = "/big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup"
GERMLINE = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
LOH_BED = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/output/seqc2_loh_only.bed"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/enriched_perlocus.json"
WINDOW = 600; ML_HIGH, ML_LOW = 200, 50; MIN_CPG = 3; MIN_SHARED = 3; MIN_GROUP = 8; N_PER = 50
bam = pysam.AlignmentFile(BAM, "rb")

# load LOH
loh = {}
for line in open(LOH_BED):
    p = line.split("\t")
    if len(p) >= 3:
        loh.setdefault(p[0], []).append((int(p[1]), int(p[2])))
def in_loh(c, pos):
    return any(s <= pos <= e for s, e in loh.get(c, []))

def hp(r):
    for t, v in r.tags:
        if t == "HP": return str(v)
    return None
def base_at(r, ref0):
    for rd, rf in r.get_aligned_pairs(matches_only=True):
        if rf == ref0: return r.query_sequence[rd] if r.query_sequence else None
    return None
def mcalls(r):
    o = {}
    try: mod = r.modified_bases
    except: return o
    if not mod: return o
    r2 = {a: b for a, b in r.get_aligned_pairs(matches_only=False) if a is not None and b is not None}
    for k, calls in mod.items():
        if k[2] != 'm': continue
        for rp, ml in calls:
            rf = r2.get(rp)
            if rf is not None: o[rf] = 1 if ml >= ML_HIGH else (0 if ml <= ML_LOW else -1)
    return o

def collect(chrom, pos, ref, alt):
    var0 = pos - 1; s, e = var0 - WINDOW, var0 + WINDOW
    reads = []
    for r in bam.fetch(chrom, max(0, s), e):
        if r.flag & 0x900 or r.flag & 0x400: continue
        m = {p: st for p, st in mcalls(r).items() if s <= p <= e and st >= 0}
        if len(m) < MIN_CPG: continue
        b = base_at(r, var0)
        allele = 'ALT' if b == alt else ('REF' if b == ref else None)
        nm = dict(r.tags).get('NM', 0)
        reads.append(dict(m=m, hp=hp(r), allele=allele, mapq=r.mapping_quality,
                          nmrate=nm / max(r.query_length or 1, 1)))
    return reads

def observed_D(mlist):
    n = len(mlist); cov = Counter()
    for r in mlist:
        for c in r: cov[c] += 1
    core = {c for c, k in cov.items() if k >= max(3, 0.25 * n)}
    mlist = [{c: r[c] for c in r if c in core} for r in mlist]
    idx = [i for i, r in enumerate(mlist) if len(r) >= MIN_CPG]
    mlist = [mlist[i] for i in idx]; n = len(mlist)
    if n < 2 * MIN_GROUP: return None, idx
    while n > 2 * MIN_GROUP:
        bad = np.array([sum(1 for j in range(n) if i != j and len(set(mlist[i]) & set(mlist[j])) < MIN_SHARED) / (n - 1) for i in range(n)])
        if bad.max() <= 0.45: break
        d = int(bad.argmax()); mlist.pop(d); idx.pop(d); n -= 1
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            sh = set(mlist[i]) & set(mlist[j])
            D[i, j] = D[j, i] = np.mean([mlist[i][c] != mlist[j][c] for c in sh]) if len(sh) >= MIN_SHARED else np.nan
    if np.isnan(D).any(): D[np.isnan(D)] = np.nanmax(D) if not np.isnan(D).all() else 0.5
    return D, idx

def analyze(chrom, pos, ref, alt, cls):
    reads = collect(chrom, pos, ref, alt)
    # need allele labels for both groups
    al = [r for r in reads if r['allele']]
    if len([r for r in al if r['allele'] == 'REF']) < MIN_GROUP or len([r for r in al if r['allele'] == 'ALT']) < MIN_GROUP:
        return None
    D, idx = observed_D([r['m'] for r in al])
    if D is None: return None
    sub = [al[i] for i in idx]
    allele_lab = np.array([0 if r['allele'] == 'REF' else 1 for r in sub])
    hp_lab = np.array([{'1': 0, '1-1': 0, '2': 1, '2-1': 1}.get(r['hp'], -1) for r in sub])
    if allele_lab.sum() < MIN_GROUP or (len(allele_lab) - allele_lab.sum()) < MIN_GROUP:
        return None
    # blind clusters
    blind = AgglomerativeClustering(n_clusters=2, metric='precomputed', linkage='average').fit_predict(D)
    ari_allele = adjusted_rand_score(allele_lab, blind)
    hp_ok = hp_lab >= 0
    ari_hp = adjusted_rand_score(hp_lab[hp_ok], blind[hp_ok]) if hp_ok.sum() >= 2 * MIN_GROUP and len(set(hp_lab[hp_ok])) == 2 else None
    ari_hp_allele = adjusted_rand_score(hp_lab[hp_ok], allele_lab[hp_ok]) if hp_ok.sum() >= 2 * MIN_GROUP and len(set(hp_lab[hp_ok])) == 2 else None
    def gstat(g):
        sel = [r for r in sub if r['allele'] == g]
        return dict(n=len(sel), mapq=float(np.mean([r['mapq'] for r in sel])), nmrate=float(np.mean([r['nmrate'] for r in sel])))
    return dict(cls=cls, locus=f"{chrom}:{pos}", loh=bool(in_loh(chrom, pos)),
                n=len(sub), vaf=round(float(allele_lab.mean()), 3),
                ari_blind_allele=round(float(ari_allele), 3),
                ari_blind_hp=round(float(ari_hp), 3) if ari_hp is not None else None,
                ari_hp_allele=round(float(ari_hp_allele), 3) if ari_hp_allele is not None else None,
                ref_stat=gstat('REF'), alt_stat=gstat('ALT'))

def load_v(vcf, n):
    out = subprocess.run(f"bcftools view {vcf} 2>/dev/null | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' 2>/dev/null | awk 'length($3)==1&&length($4)==1' | shuf | head -{n*3}", shell=True, capture_output=True, text=True)
    return [(p.split('\t')[0], int(p.split('\t')[1]), p.split('\t')[2], p.split('\t')[3]) for p in out.stdout.splitlines()]

CLASSES = {'TP': f"{PILE}/filtered_snv_tp.vcf", 'FP': f"{PILE}/filtered_snv_fp.vcf", 'FN': f"{PILE}/filtered_snv_fn.vcf"}
rows = []
for cls, vcf in CLASSES.items():
    done = 0
    for chrom, pos, ref, alt in load_v(vcf, N_PER):
        if done >= N_PER: break
        try: r = analyze(chrom, pos, ref, alt, cls)
        except Exception: continue
        if r: rows.append(r); done += 1
    print(f"[{cls}] {done} loci")
# het null (HP1 vs HP2 = both germline)
out = subprocess.run(f"bcftools view -r chr1,chr2,chr3,chr5,chr11 {GERMLINE} 2>/dev/null | bcftools view -i 'GT=\"het\"' 2>/dev/null | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' 2>/dev/null | awk 'length($3)==1&&length($4)==1' | shuf | head -150", shell=True, capture_output=True, text=True)
done = 0
for line in out.stdout.splitlines():
    if done >= N_PER: break
    p = line.split('\t')
    try: r = analyze(p[0], int(p[1]), p[2], p[3], 'het-NULL')
    except Exception: continue
    if r: rows.append(r); done += 1
print(f"[het-NULL] {done} loci")
json.dump(rows, open(OUT, "w"), indent=1)
print(f"saved {len(rows)} -> enriched_perlocus.json")
# quick Q1 summary
for cls in ['TP', 'FP', 'FN', 'het-NULL']:
    cr = [r for r in rows if r['cls'] == cls]
    if not cr: continue
    hp_aris = [r['ari_blind_hp'] for r in cr if r['ari_blind_hp'] is not None]
    al_aris = [r['ari_blind_allele'] for r in cr]
    print(f"  {cls}: median ARI(blind,HP)={np.median(hp_aris):.3f}(n={len(hp_aris)})  ARI(blind,ALT/REF)={np.median(al_aris):.3f}(n={len(al_aris)})")
