#!/usr/bin/env python3
"""
B2 clustering 人眼核圖: HP1(germline) vs HP1-1(somatic subclone) read-level methylation.
Left: MDS 2D scatter of read-read methylation distance (看兩群是否分開).
Right: read×CpG methylation heatmap, reads grouped by HP (看甲基化模式差).
變異附近 CpG 標星號 (CpG-context artifact 提示).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam
from sklearn.manifold import MDS
from skbio.stats.distance import DistanceMatrix, permanova
from sklearn.metrics import silhouette_score

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/figures/b2_clusters"
import os; os.makedirs(OUT, exist_ok=True)
WINDOW = 600; ML_HIGH, ML_LOW = 200, 50; MIN_CPG = 4; MIN_SHARED = 3
bam = pysam.AlignmentFile(BAM, "rb")
LOCI = [("chr13", 32315128, "BRCA2/ZAR1L", "1", "1-1"),
        ("chr20", 50267392, "chr20 (best silhouette)", "1", "1-1"),
        ("chr12", 31601630, "chr12 strong hyper", "2", "2-1")]

def hp(r):
    for t, v in r.tags:
        if t == "HP": return str(v)
    return None
def mc(r):
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

def viz(chrom, pos, label, gg, gs):
    var0 = pos - 1; s, e = var0 - WINDOW, var0 + WINDOW
    reads, grp = [], []
    for r in bam.fetch(chrom, max(0, s), e):
        if r.flag & 0x900 or r.flag & 0x400: continue
        h = hp(r)
        if h not in (gg, gs): continue
        m = {p: st for p, st in mc(r).items() if s <= p <= e and st >= 0}
        if len(m) >= MIN_CPG:
            reads.append(m); grp.append(0 if h == gg else 1)
    if grp.count(0) < 8 or grp.count(1) < 8:
        print(f"{chrom}:{pos} insufficient"); return
    n = len(reads); D = np.zeros((n, n)); dl = []
    for i in range(n):
        for j in range(i + 1, n):
            sh = set(reads[i]) & set(reads[j])
            if len(sh) >= MIN_SHARED:
                d = np.mean([reads[i][c] != reads[j][c] for c in sh]); D[i, j] = D[j, i] = d; dl.append(d)
            else: D[i, j] = D[j, i] = np.nan
    D[np.isnan(D)] = np.median(dl)
    grp = np.array(grp)
    dm = DistanceMatrix(D, [str(i) for i in range(n)])
    res = permanova(dm, grp.tolist(), permutations=499)
    F = res['test statistic']; R2 = (F) / (F + (n - 2)); sil = silhouette_score(D, grp, metric='precomputed')
    # MDS
    coords = MDS(n_components=2, dissimilarity='precomputed', random_state=42, normalized_stress='auto').fit_transform(D)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [1, 1.3]})
    # MDS scatter
    for gi, col, nm in [(0, '#2563eb', f'HP{gg} germline'), (1, '#dc2626', f'HP{gs} somatic-subclone')]:
        ax1.scatter(coords[grp == gi, 0], coords[grp == gi, 1], c=col, s=40, alpha=0.7, label=f'{nm} (n={int((grp==gi).sum())})', edgecolors='k', linewidths=0.3)
    ax1.set_title(f"{label}  {chrom}:{pos}\nMDS of read-read methylation distance\nPERMANOVA R²={R2:.3f} p={res['p-value']:.3f}  silhouette={sil:.3f}", fontsize=11)
    ax1.set_xlabel("MDS-1"); ax1.set_ylabel("MDS-2"); ax1.legend(fontsize=10)
    # heatmap: reads × CpG, sorted by group
    allcpg = sorted(set().union(*[set(r) for r in reads]))
    order = list(np.argsort(grp, kind='stable'))
    M = np.full((n, len(allcpg)), 0.5)
    ci = {c: k for k, c in enumerate(allcpg)}
    for ri, r in enumerate(reads):
        for c, st in r.items(): M[ri, ci[c]] = st
    M = M[order]
    ax2.imshow(M, aspect='auto', cmap='Greys', interpolation='none', vmin=0, vmax=1)
    nsplit = int((grp[order] == 0).sum())
    ax2.axhline(nsplit - 0.5, color='orange', lw=2)
    # variant + nearby CpG marker
    vx = min(range(len(allcpg)), key=lambda k: abs(allcpg[k] - var0))
    ax2.axvline(vx, color='red', lw=1.5, ls='--')
    ax2.set_title(f"read×CpG methylation (black=meth white=unmeth grey=missing)\n上={gg} germline / 下={gs} somatic; 紅虛線=somatic SNV; 橘線=分組界", fontsize=10)
    ax2.set_xlabel(f"CpG sites ({len(allcpg)})"); ax2.set_ylabel("reads (grouped)")
    plt.tight_layout()
    fn = f"{OUT}/{chrom}_{pos}_b2cluster.png"; plt.savefig(fn, dpi=120, bbox_inches='tight'); plt.close()
    print(f"{chrom}:{pos} -> {os.path.basename(fn)}  R2={R2:.3f} p={res['p-value']:.3f} sil={sil:.3f} n=({grp.tolist().count(0)},{grp.tolist().count(1)})")
    return fn

for L in LOCI: viz(*L)
print("done")
