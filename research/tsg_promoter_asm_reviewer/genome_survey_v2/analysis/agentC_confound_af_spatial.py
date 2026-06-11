#!/usr/bin/env python3
"""
Agent C — adversarial confound test (focus: AF / spatial autocorrelation).
Goal: try to REFUTE Agent B's strong-ASM / clear-consistent loci as confound artifacts.
read-only on master TSVs; writes only into analysis/.

Confound battery:
 C1. AF availability check  — is caller_af actually present? (refutes "CNV/LOH allelic-imbalance" mechanism if AF is a constant stub)
 C2. Spatial autocorrelation — are strong-ASM loci genomically CLUSTERED (non-independent pseudo-replication)?
 C3. LOH-block enrichment    — strong-ASM over-represented in LOH (where Agent B itself calls ALLELE Δβ an artifact)?
 C4. chr8 LOH hotspot overlap — HCC1395 known chr8 LOH/HPSig 7.4x FP enrichment hotspot.
 C5. Baseline-beta ceiling   — does extreme germline beta (near 0 or 1) inflate the Δβ ceiling (regression-to-extreme / AF-proxy confound)?
 C6. Coverage (n_paired_cpg) — do strong-ASM loci sit at low CpG count (small-n unstable Δβ)?
"""
import csv, math
import numpy as np
from collections import defaultdict

TP = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/asm_dualaxis_tp.tsv"

# ---- load ----
rows = []
with open(TP) as f:
    r = csv.DictReader(f, delimiter="\t")
    for d in r:
        try:
            d["somatic_pos"] = int(d["somatic_pos"])
            d["n_paired_cpg"] = int(d["n_paired_cpg"])
            d["mean_germ_beta"] = float(d["mean_germ_beta"])
            d["mean_som_beta"] = float(d["mean_som_beta"])
            d["mean_delta"] = float(d["mean_delta"])
            d["max_abs_delta"] = float(d["max_abs_delta"])
            d["wilcoxon_p"] = float(d["wilcoxon_p"]) if d["wilcoxon_p"] not in ("", "nan", "NA") else float("nan")
        except (ValueError, KeyError):
            continue
        rows.append(d)

N = len(rows)
valid_p = [r for r in rows if not math.isnan(r["wilcoxon_p"])]
BONF = 0.05 / len(valid_p)          # match Agent B alpha
print(f"[load] N={N} valid_p={len(valid_p)} Bonferroni_alpha={BONF:.3e}")

def is_strong(r):
    return (not math.isnan(r["wilcoxon_p"])) and r["wilcoxon_p"] < BONF and abs(r["mean_delta"]) >= 0.1

strong = [r for r in rows if is_strong(r)]
nonstrong = [r for r in rows if not is_strong(r)]
print(f"[strong] strong-ASM records={len(strong)}  nonstrong={len(nonstrong)}")

# unique strong loci (chrom,pos) — multiple axes per locus
strong_loci = set((r["chrom"], r["somatic_pos"]) for r in strong)
all_loci = set((r["chrom"], r["somatic_pos"]) for r in rows)
print(f"[strong] unique strong loci={len(strong_loci)}  all unique loci={len(all_loci)}")

print("\n" + "="*70)
print("C2. SPATIAL AUTOCORRELATION — are strong-ASM loci genomically clustered?")
print("="*70)
# nearest-neighbour distance among strong loci, per chrom, vs null = same count of random loci from all_loci
def nn_distances(loci):
    bychrom = defaultdict(list)
    for c, p in loci:
        bychrom[c].append(p)
    ds = []
    for c, ps in bychrom.items():
        if len(ps) < 2:
            continue
        ps = sorted(ps)
        for i, p in enumerate(ps):
            cand = []
            if i > 0: cand.append(p - ps[i-1])
            if i < len(ps)-1: cand.append(ps[i+1] - p)
            ds.append(min(cand))
    return ds

strong_nn = nn_distances(strong_loci)
all_loci_list = list(all_loci)
rng = np.random.default_rng(42)
# null: draw same number of unique loci at random from the full TP locus pool, repeat
null_med = []
null_close = []   # fraction within 500 bp
CLOSE = 500
for _ in range(2000):
    samp = set()
    idx = rng.choice(len(all_loci_list), size=len(strong_loci), replace=False)
    for i in idx:
        samp.add(all_loci_list[i])
    d = nn_distances(samp)
    if d:
        null_med.append(np.median(d))
        null_close.append(np.mean([x <= CLOSE for x in d]))

obs_med = np.median(strong_nn) if strong_nn else float("nan")
obs_close = np.mean([x <= CLOSE for x in strong_nn]) if strong_nn else float("nan")
null_med = np.array(null_med); null_close = np.array(null_close)
# p = P(null median <= observed median)  (clustering => smaller NN distance)
p_med = (np.sum(null_med <= obs_med) + 1) / (len(null_med) + 1)
p_close = (np.sum(null_close >= obs_close) + 1) / (len(null_close) + 1)
print(f"strong-ASM unique loci NN: median={obs_med:.0f} bp, frac<= {CLOSE}bp = {obs_close:.3f}")
print(f"null  (random same-N loci): median={np.median(null_med):.0f} bp [p2.5={np.percentile(null_med,2.5):.0f}, p97.5={np.percentile(null_med,97.5):.0f}]")
print(f"null  frac<= {CLOSE}bp: mean={np.mean(null_close):.3f} [p97.5={np.percentile(null_close,97.5):.3f}]")
print(f"--> clustering test: p(median<=obs)={p_med:.4f}   p(closeFrac>=obs)={p_close:.4f}")

# explicit adjacent-run detection (CpG within 500bp form a cluster)
def clusters(loci, window=500):
    bychrom = defaultdict(list)
    for c, p in loci:
        bychrom[c].append(p)
    runs = []
    for c, ps in bychrom.items():
        ps = sorted(ps)
        cur = [ps[0]]
        for p in ps[1:]:
            if p - cur[-1] <= window:
                cur.append(p)
            else:
                if len(cur) >= 2:
                    runs.append((c, cur[0], cur[-1], len(cur)))
                cur = [p]
        if len(cur) >= 2:
            runs.append((c, cur[0], cur[-1], len(cur)))
    return runs

runs = sorted(clusters(strong_loci), key=lambda x: -x[3])
in_cluster = sum(r[3] for r in runs)
print(f"\nstrong loci in >=2-member clusters (<=500bp): {in_cluster}/{len(strong_loci)} = {in_cluster/len(strong_loci):.1%}")
print("top clusters (chrom, start, end, n_loci):")
for c, s, e, n in runs[:8]:
    print(f"   {c}:{s}-{e}  n={n}  span={e-s}bp")

print("\n" + "="*70)
print("C3. LOH-BLOCK ENRICHMENT — strong-ASM over-represented in LOH?")
print("="*70)
# at RECORD level using ALLELE-axis only (valid for both LOH/nonLOH), match Agent B common baseline
allele = [r for r in rows if r["axis_type"] == "ALLELE"]
a_strong = [r for r in allele if is_strong(r)]
def loh_table(recs):
    loh = sum(1 for r in recs if r["loh_status"] == "LOH")
    non = len(recs) - loh
    return loh, non
sL, sN = loh_table(a_strong)
bL, bN = loh_table([r for r in allele if not is_strong(r)])
print(f"ALLELE-axis strong:    LOH={sL} nonLOH={sN}  -> LOH frac={sL/(sL+sN):.3f}")
print(f"ALLELE-axis nonstrong: LOH={bL} nonLOH={bN}  -> LOH frac={bL/(bL+bN):.3f}")
# fisher
try:
    from scipy.stats import fisher_exact, mannwhitneyu, binomtest, spearmanr
    OR, pf = fisher_exact([[sL, sN], [bL, bN]])
    print(f"Fisher LOH-enrichment in strong vs nonstrong: OR={OR:.3f} p={pf:.3e}")
except Exception as e:
    print("scipy missing:", e)

print("\n" + "="*70)
print("C4. chr8 LOH HOTSPOT — HCC1395 known chr8 LOH/HPSig FP-enrichment hotspot")
print("="*70)
chr8_strong = [l for l in strong_loci if l[0] == "chr8"]
chr8_all = [l for l in all_loci if l[0] == "chr8"]
print(f"chr8 strong loci={len(chr8_strong)}  chr8 all loci={len(chr8_all)}")
print(f"strong-on-chr8 fraction={len(chr8_strong)/len(strong_loci):.3f}  vs all-loci-on-chr8={len(chr8_all)/len(all_loci):.3f}")
try:
    a, b = len(chr8_strong), len(strong_loci)-len(chr8_strong)
    c8, d8 = len(chr8_all)-len(chr8_strong), (len(all_loci)-len(chr8_all))-(len(strong_loci)-len(chr8_strong))
    OR8, p8 = fisher_exact([[a, b],[c8, d8]])
    print(f"Fisher chr8 enrichment in strong: OR={OR8:.3f} p={p8:.3e}")
except Exception as e:
    print(e)
print("chr8 strong loci positions:", sorted(p for c,p in chr8_strong))

print("\n" + "="*70)
print("C5. BASELINE-BETA CEILING — does extreme germline beta inflate |Δβ|?")
print("="*70)
# regression-to-extreme / AF proxy: if germ beta is near 0 or 1, |Δβ| is mechanically bounded;
# but if strong loci cluster at extreme germ beta, the 'effect' may be a beta-ceiling artifact.
gb_strong = np.array([r["mean_germ_beta"] for r in strong])
gb_all = np.array([r["mean_germ_beta"] for r in rows])
def extremity(x):   # distance to nearest bound, small = extreme
    return np.minimum(x, 1-x)
ext_strong = extremity(gb_strong)
ext_all = extremity(gb_all)
print(f"germ beta extremity (min(b,1-b)): strong median={np.median(ext_strong):.3f}  all median={np.median(ext_all):.3f}")
frac_ext_strong = np.mean(ext_strong < 0.1)
frac_ext_all = np.mean(ext_all < 0.1)
print(f"frac with germ beta in extreme tail (<0.1 from 0/1): strong={frac_ext_strong:.3f}  all={frac_ext_all:.3f}")
try:
    U, pmw = mannwhitneyu(ext_strong, ext_all, alternative="less")
    print(f"Mann-Whitney (strong more extreme?): U={U:.3e} p={pmw:.3e}")
    # correlation: across ALL records, does |delta| track germ-beta extremity?
    rho, prho = spearmanr(np.abs([r["mean_delta"] for r in rows]), ext_all)
    print(f"Spearman |Δβ| vs germ-extremity (all records): rho={rho:.3f} p={prho:.3e}  (negative rho => extreme beta inflates |Δβ|)")
except Exception as e:
    print(e)

print("\n" + "="*70)
print("C6. COVERAGE (n_paired_cpg) — strong-ASM at low CpG count (unstable)?")
print("="*70)
np_strong = np.array([r["n_paired_cpg"] for r in strong])
np_all = np.array([r["n_paired_cpg"] for r in rows])
print(f"n_paired_cpg: strong median={np.median(np_strong):.0f} mean={np.mean(np_strong):.1f}  | all median={np.median(np_all):.0f} mean={np.mean(np_all):.1f}")
print(f"strong with n_paired_cpg < 5: {np.mean(np_strong<5):.3f}  | all: {np.mean(np_all<5):.3f}")
try:
    U, pc = mannwhitneyu(np_strong, np_all, alternative="less")
    print(f"Mann-Whitney (strong lower coverage?): p={pc:.3e}")
except Exception as e:
    print(e)

# ---- focused check on Agent B clear-consistent 32 loci: how many survive C2(cluster)+C4(chr8)+C5(extreme)? ----
print("\n" + "="*70)
print("VERDICT TABLE — Agent B 32 clear-consistent loci vs confound flags")
print("="*70)
cc = []
with open("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/analysis/clear_consistent_loci.csv") as f:
    for d in csv.DictReader(f):
        cc.append((d["chrom"], int(d["somatic_pos"]), d["al_loh"], float(d["combined_abs_delta"])))

cluster_members = set()
for c, s, e, n in clusters(strong_loci):
    bychrom = defaultdict(list)
    for cc2, p in strong_loci:
        bychrom[cc2].append(p)
    for p in sorted(bychrom[c]):
        if s <= p <= e:
            cluster_members.add((c, p))

# beta extremity per cc locus (use the strong record at that locus)
rec_by_locus = defaultdict(list)
for r in strong:
    rec_by_locus[(r["chrom"], r["somatic_pos"])].append(r)

n_loh = n_cluster = n_chr8 = n_extreme = 0
flag_clean = 0
for c, p, loh, cab in cc:
    f_loh = (loh == "LOH")
    f_cluster = (c, p) in cluster_members
    f_chr8 = (c == "chr8")
    recs = rec_by_locus.get((c, p), [])
    gb = np.mean([rr["mean_germ_beta"] for rr in recs]) if recs else float("nan")
    f_extreme = (not math.isnan(gb)) and (min(gb, 1-gb) < 0.1)
    n_loh += f_loh; n_cluster += f_cluster; n_chr8 += f_chr8; n_extreme += f_extreme
    flags = []
    if f_loh: flags.append("LOH")
    if f_cluster: flags.append("CLUSTER")
    if f_chr8: flags.append("chr8")
    if f_extreme: flags.append("BETA-EXTREME")
    if not flags:
        flag_clean += 1
    print(f"  {c}:{p}  |Δ|={cab:.3f}  germβ={gb:.2f}  flags={flags if flags else 'CLEAN'}")

print(f"\nof 32 clear-consistent loci: LOH={n_loh}  in-cluster={n_cluster}  chr8={n_chr8}  beta-extreme={n_extreme}  fully-CLEAN(no confound flag)={flag_clean}")
