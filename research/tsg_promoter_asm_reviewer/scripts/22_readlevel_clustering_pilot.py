#!/usr/bin/env python3
"""
B2 嚴格 pilot — read-level methylation CLUSTERING (Layer B2, NOT per-CpG mean).
Tests whether somatic-subclone reads (HP1-1/HP2-1) form a DISTINCT, cohesive
methylation cluster vs germline-haplotype reads (HP1/HP2), with 3 controls:
  (1) n_reads control: subsample germline group to somatic group size
  (2) LOH/CNV annotation
  (3) germline-het null: same analysis on het SNPs (HP1 vs HP2 = baseline allelic)

Stats: skbio PERMANOVA (pseudo-F, R2, p via permutation) + dispersion ratio + silhouette.
Distance: pairwise Hamming over shared CpGs (read-read methylation distance matrix).

Key question: do somatic loci show HIGHER PERMANOVA R2 (n_reads-matched) than germline-het null?
"""
import sys, gzip, random
import numpy as np
import pysam
from scipy import stats
from skbio.stats.distance import DistanceMatrix, permanova
from sklearn.metrics import silhouette_score

random.seed(42); np.random.seed(42)
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GERMLINE = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
LOH_BED = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/output/seqc2_loh_only.bed"
WINDOW = 600
ML_HIGH, ML_LOW = 200, 50
MIN_CPG_PER_READ = 4      # read 需覆蓋 >=4 個 window CpG
MIN_READS_PER_GROUP = 8   # 每組 >=8 reads
MIN_SHARED_CPG = 3        # 兩 read 算距離需 >=3 共同 CpG
N_PERM = 499

bam = pysam.AlignmentFile(BAM, "rb")

# real-case somatic loci (top robust from phase2)
SOMATIC = [
    ("chr13", 32315128, "BRCA2/ZAR1L"),
    ("chr9", 42376881, "biggest HP effect"),
    ("chr22", 34684238, "strong hypo"),
    ("chr12", 31601630, "strong hyper"),
    ("chr16", 17774746, "hypo"),
    ("chr2", 94874876, "strong hypo n=86"),
    ("chr1", 200029019, "biggest |Δβ|"),
    ("chr20", 50267392, "n=74 strong"),
]

def extract_hp(read):
    for t, v in read.tags:
        if t == "HP":
            return str(v)
    return None

def meth_calls(read):
    out = {}
    try:
        mod = read.modified_bases
    except Exception:
        return out
    if not mod:
        return out
    r2r = {rd: rf for rd, rf in read.get_aligned_pairs(matches_only=False) if rd is not None and rf is not None}
    for key, calls in mod.items():
        if key[2] != 'm':
            continue
        for rp, ml in calls:
            rf = r2r.get(rp)
            if rf is not None:
                out[rf] = 1 if ml >= ML_HIGH else (0 if ml <= ML_LOW else -1)  # -1 ambiguous
    return out

def collect(chrom, pos, groups_wanted):
    """Return {hp: [ {cpg:state}, ... ]} for reads with HP in groups_wanted."""
    var0 = pos - 1
    start, end = var0 - WINDOW, var0 + WINDOW
    g = {k: [] for k in groups_wanted}
    for read in bam.fetch(chrom, max(0, start), end):
        if read.flag & 0x900 or read.flag & 0x400:
            continue
        hp = extract_hp(read)
        if hp not in groups_wanted:
            continue
        m = {p: s for p, s in meth_calls(read).items() if start <= p <= end and s >= 0}
        if len(m) >= MIN_CPG_PER_READ:
            g[hp].append(m)
    return g

def build_distance(reads):
    """reads = list of {cpg:state}; pairwise Hamming over shared CpG."""
    n = len(reads)
    D = np.zeros((n, n))
    valid_pairs = 0
    dlist = []
    for i in range(n):
        for j in range(i + 1, n):
            shared = set(reads[i]) & set(reads[j])
            if len(shared) >= MIN_SHARED_CPG:
                dis = np.mean([reads[i][c] != reads[j][c] for c in shared])
                D[i, j] = D[j, i] = dis
                dlist.append(dis); valid_pairs += 1
            else:
                D[i, j] = D[j, i] = np.nan
    # impute nan with median
    if dlist:
        med = np.median(dlist)
        D[np.isnan(D)] = med
    return D, valid_pairs

def analyze(chrom, pos, g_germ, g_som, label):
    g = collect(chrom, pos, [g_germ, g_som])
    n_germ, n_som = len(g[g_germ]), len(g[g_som])
    if n_germ < MIN_READS_PER_GROUP or n_som < MIN_READS_PER_GROUP:
        return None
    def run(germ_reads, som_reads):
        reads = germ_reads + som_reads
        labels = [0] * len(germ_reads) + [1] * len(som_reads)
        D, vp = build_distance(reads)
        ids = [str(i) for i in range(len(reads))]
        dm = DistanceMatrix(D, ids)
        res = permanova(dm, labels, permutations=N_PERM)
        # R2 = SS_between/SS_total
        F = res['test statistic']; N = len(reads); k = 2
        R2 = (F * (k - 1)) / (F * (k - 1) + (N - k))
        # dispersion: within-group mean dist
        def wdisp(idx):
            sub = D[np.ix_(idx, idx)]
            return np.mean(sub[np.triu_indices(len(idx), 1)]) if len(idx) > 1 else np.nan
        gi = [i for i, l in enumerate(labels) if l == 0]; si = [i for i, l in enumerate(labels) if l == 1]
        disp_g, disp_s = wdisp(gi), wdisp(si)
        try:
            sil = silhouette_score(D, labels, metric='precomputed')
        except Exception:
            sil = np.nan
        return dict(F=F, R2=R2, p=res['p-value'], disp_germ=disp_g, disp_som=disp_s, sil=sil)
    full = run(g[g_germ], g[g_som])
    # n_reads control: subsample BOTH groups to min size, median over 5 draws
    m = min(n_germ, n_som)
    matched = []
    for _ in range(3):
        gsub = random.sample(g[g_germ], m)
        ssub = random.sample(g[g_som], m)
        matched.append(run(gsub, ssub))
    R2m = np.median([m['R2'] for m in matched]); pm = np.median([m['p'] for m in matched])
    return dict(locus=f"{chrom}:{pos}", label=label, n_germ=n_germ, n_som=n_som,
                R2_full=full['R2'], p_full=full['p'], R2_matched=R2m, p_matched=pm,
                disp_germ=full['disp_germ'], disp_som=full['disp_som'],
                disp_ratio=full['disp_som'] / full['disp_germ'] if full['disp_germ'] else np.nan,
                sil=full['sil'])

# germline-het null sites (HP1 vs HP2 baseline)
def sample_het(n=20):
    sites = []
    import subprocess
    out = subprocess.run(
        f"bcftools view -r chr1,chr2,chr3 {GERMLINE} 2>/dev/null | bcftools view -i 'GT=\"het\"' 2>/dev/null | "
        f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' 2>/dev/null | awk 'length($3)==1&&length($4)==1' | shuf | head -200",
        shell=True, capture_output=True, text=True)
    for line in out.stdout.splitlines():
        p = line.split("\t")
        sites.append((p[0], int(p[1])))
        if len(sites) >= n * 4:
            break
    return sites

print("=== SOMATIC loci: HP1/HP2(germline) vs HP1-1/HP2-1(somatic subclone) ===")
print(f"{'locus':<18}{'label':<18}{'n_g':>5}{'n_s':>5}{'R2_full':>9}{'p_full':>8}{'R2_match':>10}{'p_match':>9}{'disp_ratio':>11}{'sil':>7}")
som_results = []
for chrom, pos, lab in SOMATIC:
    for gg, gs in [("1", "1-1"), ("2", "2-1")]:
        r = analyze(chrom, pos, gg, gs, lab)
        if r:
            som_results.append(r)
            print(f"{r['locus']:<18}{(lab+'/'+gg):<18}{r['n_germ']:>5}{r['n_som']:>5}{r['R2_full']:>9.3f}{r['p_full']:>8.3f}{r['R2_matched']:>10.3f}{r['p_matched']:>9.3f}{r['disp_ratio']:>11.3f}{r['sil']:>7.3f}")

print("\n=== GERMLINE-HET NULL: HP1 vs HP2 (baseline allelic clustering) ===")
het_results = []
for chrom, pos in sample_het(15):
    r = analyze(chrom, pos, "1", "2", "het-null")
    if r:
        het_results.append(r)
if het_results:
    print(f"  n={len(het_results)} het sites | R2_matched median={np.median([r['R2_matched'] for r in het_results]):.3f} | p_matched<0.05 frac={np.mean([r['p_matched']<0.05 for r in het_results]):.2f}")

print("\n=== 關鍵比較: somatic vs het-null (n_reads-matched PERMANOVA R2) ===")
if som_results and het_results:
    sr = [r['R2_matched'] for r in som_results]; hr = [r['R2_matched'] for r in het_results]
    sp = np.mean([r['p_matched'] < 0.05 for r in som_results]); hp = np.mean([r['p_matched'] < 0.05 for r in het_results])
    u, pu = stats.mannwhitneyu(sr, hr, alternative='greater')
    print(f"  somatic R2_matched median={np.median(sr):.3f} (p<.05 frac={sp:.2f}, n={len(sr)})")
    print(f"  het-null R2_matched median={np.median(hr):.3f} (p<.05 frac={hp:.2f}, n={len(hr)})")
    print(f"  Mann-Whitney somatic>het: U={u:.0f} p={pu:.3f}")
    print(f"  VERDICT: {'somatic subclone 結構 > baseline (B2 POSITIVE)' if pu<0.05 and np.median(sr)>np.median(hr) else 'somatic ≈ baseline (B2 NEGATIVE, 確認 O11)'}")
