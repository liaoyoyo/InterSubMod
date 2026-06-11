#!/usr/bin/env python3
"""
B2 broad scan with VALIDATED ruler: blind-ARI + observed-only + placebo across TP/FP/FN/het-null.
Split: ALT-allele reads vs REF-allele reads at the variant (universal across classes).
Goal: (1) push BRCA2-like loci; (2) which variant class (TP/FP/FN) shows B2 clustering.
M0 relaxed (core-CpG threshold lowered) so fewer loci over-drop.
Output: per-locus TSV + per-class summary.
"""
import numpy as np, random, subprocess, json
import pysam
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score
from collections import Counter
random.seed(42); np.random.seed(42)

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
PILE = "/big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup"
GERMLINE = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
OUTDIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2"
WINDOW = 600; ML_HIGH, ML_LOW = 200, 50
MIN_CORE_CPG = 3; MIN_SHARED = 3; MIN_GROUP = 8; N_PER_CLASS = 45
bam = pysam.AlignmentFile(BAM, "rb")

def base_at(read, ref0):
    for rd, rf in read.get_aligned_pairs(matches_only=True):
        if rf == ref0:
            return read.query_sequence[rd] if read.query_sequence else None
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

def collect_by_allele(chrom, pos, ref, alt):
    var0 = pos - 1; s, e = var0 - WINDOW, var0 + WINDOW
    g = {'REF': [], 'ALT': []}; lens = {'REF': [], 'ALT': []}
    for r in bam.fetch(chrom, max(0, s), e):
        if r.flag & 0x900 or r.flag & 0x400: continue
        b = base_at(r, var0)
        if b == ref: grp = 'REF'
        elif b == alt: grp = 'ALT'
        else: continue
        m = {p: st for p, st in mcalls(r).items() if s <= p <= e and st >= 0}
        if len(m) >= MIN_CORE_CPG:
            g[grp].append(m); lens[grp].append(r.query_length or (r.reference_end - r.reference_start))
    return g, lens

def observed_only_D(reads):
    n = len(reads)
    cov = Counter()
    for r in reads:
        for c in r: cov[c] += 1
    core = {c for c, k in cov.items() if k >= max(3, 0.25 * n)}
    reads = [{c: r[c] for c in r if c in core} for r in reads]
    keep = [i for i, r in enumerate(reads) if len(r) >= MIN_CORE_CPG]
    reads = [reads[i] for i in keep]; n = len(reads)
    if n < 2 * MIN_GROUP: return None, keep
    for _ in range(n):
        if n <= 2 * MIN_GROUP: break
        bad = np.zeros(n)
        for i in range(n):
            inc = sum(1 for j in range(n) if i != j and len(set(reads[i]) & set(reads[j])) < MIN_SHARED)
            bad[i] = inc / (n - 1)
        if bad.max() <= 0.45: break
        d = int(bad.argmax()); reads.pop(d); keep.pop(d); n -= 1
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            sh = set(reads[i]) & set(reads[j])
            D[i, j] = D[j, i] = np.mean([reads[i][c] != reads[j][c] for c in sh]) if len(sh) >= MIN_SHARED else np.nan
    if np.isnan(D).any(): D[np.isnan(D)] = np.nanmax(D) if not np.isnan(D).all() else 0.5
    return D, keep

def blind_ari(D, lab):
    preds = []
    for link in ('average', 'complete'):
        preds.append(AgglomerativeClustering(n_clusters=2, metric='precomputed', linkage=link).fit_predict(D))
    try:
        preds.append(SpectralClustering(n_clusters=2, affinity='precomputed', random_state=42).fit_predict(np.exp(-D / (D.std() + 1e-9))))
    except Exception: pass
    return max(adjusted_rand_score(lab, p) for p in preds), max(adjusted_mutual_info_score(lab, p) for p in preds)

def eval_split(g0, g1):
    reads = g0 + g1; lab0 = [0] * len(g0) + [1] * len(g1)
    D, keep = observed_only_D(reads)
    if D is None: return None
    lab = np.array([lab0[i] for i in keep])
    if lab.sum() < MIN_GROUP or (len(lab) - lab.sum()) < MIN_GROUP: return None
    ari, ami = blind_ari(D, lab)
    sil = silhouette_score(D, lab, metric='precomputed')
    return dict(n0=int((lab == 0).sum()), n1=int((lab == 1).sum()), ari=ari, ami=ami, sil=sil)

def placebo(g_ref, lens_ref):
    """split REF reads by length (mimic ALT-spanning selection) -> blind-ARI should be low."""
    if len(g_ref) < 2 * MIN_GROUP: return np.nan
    med = np.median(lens_ref)
    a = [g_ref[i] for i in range(len(g_ref)) if lens_ref[i] >= med]
    b = [g_ref[i] for i in range(len(g_ref)) if lens_ref[i] < med]
    if len(a) < MIN_GROUP or len(b) < MIN_GROUP: return np.nan
    r = eval_split(b, a)
    return r['ari'] if r else np.nan

def load_variants(vcf, n):
    out = subprocess.run(f"bcftools view -f PASS {vcf} 2>/dev/null | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' 2>/dev/null | "
                         f"awk 'length($3)==1&&length($4)==1' | shuf | head -{n*3}", shell=True, capture_output=True, text=True)
    v = []
    for line in out.stdout.splitlines():
        p = line.split("\t"); v.append((p[0], int(p[1]), p[2], p[3]))
        if len(v) >= n * 3: break
    return v

CLASSES = {
    'TP': f"{PILE}/filtered_snv_tp.vcf",
    'FP': f"{PILE}/filtered_snv_fp.vcf",
    'FN': f"{PILE}/filtered_snv_fn.vcf",
    'het-NULL': GERMLINE,
}
rows = []
print(f"{'class':<10}{'locus':<18}{'n0/n1':>9}{'blind_ARI':>10}{'sil':>7}{'placebo':>9}{'cluster?':>10}")
for cls, vcf in CLASSES.items():
    if cls == 'het-NULL':
        out = subprocess.run(f"bcftools view -r chr1,chr2,chr5,chr11,chr17 {vcf} 2>/dev/null | bcftools view -i 'GT=\"het\"' 2>/dev/null | "
                             f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' 2>/dev/null | awk 'length($3)==1&&length($4)==1' | shuf | head -150",
                             shell=True, capture_output=True, text=True)
        variants = [(p.split('\t')[0], int(p.split('\t')[1]), p.split('\t')[2], p.split('\t')[3]) for p in out.stdout.splitlines()]
    else:
        variants = load_variants(vcf, N_PER_CLASS)
    done = 0
    for chrom, pos, ref, alt in variants:
        if done >= N_PER_CLASS: break
        try:
            g, lens = collect_by_allele(chrom, pos, ref, alt)
        except Exception: continue
        if len(g['REF']) < MIN_GROUP or len(g['ALT']) < MIN_GROUP: continue
        r = eval_split(g['REF'], g['ALT'])
        if not r: continue
        pl = placebo(g['REF'], lens['REF'])
        is_cluster = r['ari'] >= 0.30 and (np.isnan(pl) or pl < 0.15)
        rows.append(dict(cls=cls, locus=f"{chrom}:{pos}", ari=r['ari'], ami=r['ami'], sil=r['sil'],
                         placebo=float(pl) if not np.isnan(pl) else None, n0=r['n0'], n1=r['n1'], is_cluster=bool(is_cluster)))
        done += 1
        if r['ari'] >= 0.30:
            print(f"{cls:<10}{chrom+':'+str(pos):<18}{str(r['n0'])+'/'+str(r['n1']):>9}{r['ari']:>10.3f}{r['sil']:>7.3f}{(pl if not np.isnan(pl) else -1):>9.3f}{'YES' if is_cluster else 'placebo?':>10}")
    print(f"  [{cls}] evaluated {done} loci")

# per-class summary
print("\n=== PER-CLASS B2-cluster RATE (blind-ARI>=0.30 + placebo<0.15) ===")
import numpy as np
for cls in CLASSES:
    cr = [r for r in rows if r['cls'] == cls]
    if not cr: print(f"  {cls}: 0 evaluated"); continue
    rate = np.mean([r['is_cluster'] for r in cr]); med_ari = np.median([r['ari'] for r in cr])
    n_clu = sum(r['is_cluster'] for r in cr)
    print(f"  {cls:<10} n={len(cr):3d}  median_ARI={med_ari:.3f}  B2-cluster {n_clu}/{len(cr)} = {100*rate:.1f}%")

json.dump(rows, open(f"{OUTDIR}/b2_broad_scan_results.json", "w"), indent=1)
print(f"\nsaved {len(rows)} loci -> b2_broad_scan_results.json")
# key comparison
tp = [r['ari'] for r in rows if r['cls'] == 'TP']; fp = [r['ari'] for r in rows if r['cls'] == 'FP']
fn = [r['ari'] for r in rows if r['cls'] == 'FN']; het = [r['ari'] for r in rows if r['cls'] == 'het-NULL']
from scipy.stats import mannwhitneyu
if tp and het:
    u, p = mannwhitneyu(tp, het, alternative='greater')
    print(f"TP-ARI median={np.median(tp):.3f} vs het-NULL={np.median(het):.3f}  MW(TP>het) p={p:.3f}")
if tp and fp:
    print(f"FP-ARI median={np.median(fp):.3f}  FN-ARI median={np.median(fn) if fn else float('nan'):.3f}")
