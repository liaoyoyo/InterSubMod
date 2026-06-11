#!/usr/bin/env python3
"""
M2 - Methylation MEASUREMENT scrutiny + subclone proxy (HCC1395, single sample, A pilot).

Question (user): is the current beta measurement biased, and does it hide subclone?
Current beta (18_dual_axis_pivot.py) = frac reads with P>=0.5 AFTER MAX-collapsing
5mC('+m') and 5hmC('+h') per-read ML, then binarizing at ML>=128.

Reuses BAM-read machinery from scripts/24_cluster_eval_core.py (modified_bases, get_aligned_pairs)
but fetches 5mC and 5hmC ML SEPARATELY (0-255) instead of max-collapse + ternary.

ANALYSES:
  (A) BINARIZATION sensitivity: per-locus somatic-group beta 3 ways:
        (i)  binarized P>=0.5  (current: ML>=128 on MAX-collapsed)
        (ii) CONTINUOUS mean(ML/255) on max-collapsed
        (iii) binarized at 0.3 (ML>=76.5) and 0.7 (ML>=178.5)
      Correlate (i) vs (ii); flag loci that diverge most (many intermediate-ML reads).
  (B) 5mC vs 5hmC: per-locus ASM |db| (germ vs somatic HP group) with 5mC-ONLY,
      5hmC-ONLY, MAX-collapsed. Correlation across loci. Frac of CpG-read calls
      where 5hmC ML > 5mC ML.
  (C) DEAD-ZONE: frac of read-CpG calls (max-collapsed ML) in [50,200] overall + per locus.
  (D) SUBCLONE PROXY: per-read mean methylation (continuous, max-collapse ML/255) for the
      somatic HP group. Bimodality via GMM 1-vs-2 BIC (diptest unavailable). Frac of loci
      bimodal. Cross-tab bimodal-frac by is_tp and by cn_class.

Output:
  json:  genome_survey_v2/cn_confound/m2_meth_measurement_audit.json
  fig:   figures/cn_confound/m2_meth_measurement.png
"""
import os, sys, json, random
os.environ.setdefault("TMPDIR", "/big7_disk")
import numpy as np
import pandas as pd
import pysam
from scipy.stats import pearsonr, spearmanr
from sklearn.mixture import GaussianMixture
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

random.seed(42); np.random.seed(42)

CJK = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
FP = FontProperties(fname=CJK) if os.path.exists(CJK) else None

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
MASTER = f"{ROOT}/genome_survey_v2/cn_confound/master_o1_cn.tsv"
OUT_JSON = f"{ROOT}/genome_survey_v2/cn_confound/m2_meth_measurement_audit.json"
OUT_FIG = f"{ROOT}/figures/cn_confound/m2_meth_measurement.png"
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"

WINDOW = 600
MIN_CPG_PER_READ = 4         # min CpG calls for a read to count (matches 24_core MIN_CORE_CPG)
MIN_GROUP = 8                # min reads per group
ML_BIN = 128                 # P>=0.5 binarization (current beta)
ML_LO, ML_HI = 50, 200       # dead-zone bounds (matches 24_core ML_LOW/ML_HIGH)
ML_30, ML_70 = 0.3 * 255, 0.7 * 255

bam = pysam.AlignmentFile(BAM, "rb")


def hp_tag(r):
    for t, v in r.tags:
        if t == "HP":
            return str(v)
    return None


def read_ml_separate(r):
    """Return dict refpos -> (ml_5mC, ml_5hmC) (0-255, or None if that mod absent at that pos)."""
    try:
        mod = r.modified_bases
    except Exception:
        return {}
    if not mod:
        return {}
    r2 = {a: b for a, b in r.get_aligned_pairs(matches_only=False)
          if a is not None and b is not None}
    out = {}
    for k, calls in mod.items():
        code = k[2]  # 'm' (5mC) or 'h' (5hmC)
        if code not in ('m', 'h'):
            continue
        for rp, ml in calls:
            rf = r2.get(rp)
            if rf is None:
                continue
            if rf not in out:
                out[rf] = [None, None]
            if code == 'm':
                out[rf][0] = ml
            else:
                out[rf][1] = ml
    return {p: tuple(v) for p, v in out.items()}


def collect_loci(chrom, pos, groups):
    """For a locus, fetch reads in +/-WINDOW, return per-group list of per-read CpG dicts.
    Each entry: dict refpos -> (ml_5mC, ml_5hmC). Filters by MIN_CPG_PER_READ (any-mod present)."""
    var0 = pos - 1
    s, e = var0 - WINDOW, var0 + WINDOW
    g = {k: [] for k in groups}
    for r in bam.fetch(chrom, max(0, s), e):
        if r.flag & 0x900 or r.flag & 0x400:
            continue
        h = hp_tag(r)
        if h not in groups:
            continue
        calls = {p: v for p, v in read_ml_separate(r).items()
                 if s <= p <= e and (v[0] is not None or v[1] is not None)}
        if len(calls) >= MIN_CPG_PER_READ:
            g[h].append(calls)
    return g


def collapse_max(v):
    """max-collapse a (ml_5mC, ml_5hmC) tuple -> single ML (0-255). Skip None."""
    a, b = v
    vals = [x for x in (a, b) if x is not None]
    return max(vals) if vals else None


def group_calls_flat(reads, which):
    """Flatten all CpG calls of a group into ML list.
    which in {'max','mC','hC'}. Returns list of ML (0-255), skipping None for that channel."""
    out = []
    for rd in reads:
        for v in rd.values():
            if which == 'max':
                m = collapse_max(v)
            elif which == 'mC':
                m = v[0]
            else:
                m = v[1]
            if m is not None:
                out.append(m)
    return out


def beta_binarized(mls, thr):
    """frac of calls with ML >= thr."""
    if not mls:
        return np.nan
    return float(np.mean([1 if m >= thr else 0 for m in mls]))


def beta_continuous(mls):
    """mean ML / 255."""
    if not mls:
        return np.nan
    return float(np.mean(mls) / 255.0)


def per_read_mean_meth(reads, which='max'):
    """per-read mean methylation (continuous ML/255). Returns array of read-level betas."""
    out = []
    for rd in reads:
        vals = []
        for v in rd.values():
            if which == 'max':
                m = collapse_max(v)
            elif which == 'mC':
                m = v[0]
            else:
                m = v[1]
            if m is not None:
                vals.append(m / 255.0)
        if vals:
            out.append(float(np.mean(vals)))
    return np.array(out)


def gmm_bimodality(x):
    """1-vs-2 component GMM BIC test. Returns dict with is_bimodal, bic1, bic2, delta_bic,
    component means/weights. Bimodal if BIC2 < BIC1 by >=10 (strong) AND both comps weight>0.10
    AND mean separation > 0.10 (on [0,1] scale)."""
    n = len(x)
    res = dict(n=int(n), is_bimodal=False, bic1=None, bic2=None, delta_bic=None,
               comp_means=None, comp_weights=None, mean_sep=None)
    if n < 12:
        res['reason'] = 'n<12'
        return res
    X = x.reshape(-1, 1)
    try:
        g1 = GaussianMixture(n_components=1, covariance_type='full',
                             random_state=42, n_init=2).fit(X)
        g2 = GaussianMixture(n_components=2, covariance_type='full',
                             random_state=42, n_init=5).fit(X)
    except Exception as ex:
        res['reason'] = f'gmm_fail:{ex}'
        return res
    bic1, bic2 = g1.bic(X), g2.bic(X)
    means = sorted(g2.means_.ravel().tolist())
    w = g2.weights_.ravel()
    minw = float(min(w))
    sep = float(abs(means[1] - means[0]))
    res.update(bic1=float(bic1), bic2=float(bic2), delta_bic=float(bic1 - bic2),
               comp_means=[round(m, 4) for m in means],
               comp_weights=[round(float(x), 4) for x in sorted(w.tolist())],
               mean_sep=round(sep, 4), min_weight=round(minw, 4))
    res['is_bimodal'] = bool((bic1 - bic2) >= 10.0 and minw >= 0.10 and sep >= 0.10)
    return res


# ---------- locus sampling ----------
df = pd.read_csv(MASTER, sep="\t", dtype={'chrom': str})
hp = df[df.axis_type == "HP"].copy()
hp['wilcoxon_p'] = pd.to_numeric(hp['wilcoxon_p'], errors='coerce')

CANON = [("chr13", 32315128, "BRCA2"), ("chr20", 50267392, "chr20_50267392"),
         ("chr12", 31601630, "chr12_31601630"), ("chr9", 42376881, "chr9_42376881"),
         ("chr16", 17774746, "chr16_17774746")]
canon_keys = {(c, p) for c, p, _ in CANON}


def axis_to_groups(axis):
    # 'HP1_vs_HP1-1' -> ('1','1-1'); 'HP2_vs_HP2-1' -> ('2','2-1')
    left, right = axis.split("_vs_")
    return left.replace("HP", ""), right.replace("HP", "")


# Build candidate locus rows: prefer the strongest axis per canonical, then random sig / non-sig
def row_to_spec(row, label=None):
    gg, gs = axis_to_groups(row['axis'])
    return dict(chrom=str(row['chrom']), pos=int(row['somatic_pos']), axis=row['axis'],
                gg=gg, gs=gs, is_tp=int(row['is_tp']),
                cn_class=str(row['cn_class']), median_cn=float(row['median_cn']),
                loh_status=str(row['loh_status']), wilcoxon_p=float(row['wilcoxon_p'])
                if pd.notna(row['wilcoxon_p']) else None,
                n_paired_cpg=int(row['n_paired_cpg']), abs_delta=float(row['abs_delta']),
                label=label, is_canonical=label is not None)


specs = []
seen = set()
# canonical: pick the row with smallest wilcoxon_p among that position's HP rows
for c, p, lab in CANON:
    sub = hp[(hp.chrom == c) & (hp.somatic_pos == p)]
    if len(sub) == 0:
        continue
    best = sub.sort_values('wilcoxon_p').iloc[0]
    specs.append(row_to_spec(best, label=lab))
    seen.add((c, p))

# remaining HP rows, exclude canonical positions; one row per (chrom,pos) keeping strongest axis
hp_rest = hp[~hp.apply(lambda r: (r['chrom'], r['somatic_pos']) in seen, axis=1)].copy()
hp_rest = hp_rest.sort_values('wilcoxon_p').drop_duplicates(['chrom', 'somatic_pos'], keep='first')

sig = hp_rest[(hp_rest.wilcoxon_p.notna()) & (hp_rest.wilcoxon_p < 0.05)]
nonsig = hp_rest[(hp_rest.wilcoxon_p.notna()) & (hp_rest.wilcoxon_p >= 0.05)]

sig_samp = sig.sample(n=min(150, len(sig)), random_state=42)
nonsig_samp = nonsig.sample(n=min(150, len(nonsig)), random_state=42)
for _, row in sig_samp.iterrows():
    specs.append(row_to_spec(row))
for _, row in nonsig_samp.iterrows():
    specs.append(row_to_spec(row))

print(f"[sampling] {len(specs)} loci: {len([s for s in specs if s['is_canonical']])} canonical "
      f"+ {len(sig_samp)} sig + {len(nonsig_samp)} non-sig")

# ---------- per-locus computation ----------
per_locus = []
all_deadzone_calls = 0
all_total_calls = 0
all_5hmC_gt_5mC = 0
all_pairwise_calls = 0   # calls where BOTH 5mC and 5hmC present (for fair 5hmC>5mC comparison)

for i, sp in enumerate(specs):
    g = collect_loci(sp['chrom'], sp['pos'], [sp['gg'], sp['gs']])
    n_germ, n_som = len(g[sp['gg']]), len(g[sp['gs']])
    rec = dict(sp)
    rec.update(n_germ_reads=n_germ, n_som_reads=n_som)
    if n_som < MIN_GROUP or n_germ < MIN_GROUP:
        rec['skipped'] = f'insufficient reads germ={n_germ} som={n_som}'
        per_locus.append(rec)
        continue

    som_reads = g[sp['gs']]
    germ_reads = g[sp['gg']]

    # (A) binarization sensitivity on somatic group (max-collapsed)
    som_max = group_calls_flat(som_reads, 'max')
    beta_i = beta_binarized(som_max, ML_BIN)        # current
    beta_ii = beta_continuous(som_max)              # continuous
    beta_03 = beta_binarized(som_max, ML_30)
    beta_07 = beta_binarized(som_max, ML_70)
    frac_interm = float(np.mean([1 if ML_LO < m < ML_HI else 0 for m in som_max])) if som_max else np.nan

    # (B) 5mC vs 5hmC ASM |db| (germ vs som)
    def asm_db(which):
        gm = group_calls_flat(germ_reads, which)
        sm = group_calls_flat(som_reads, which)
        if not gm or not sm:
            return np.nan
        return abs(beta_continuous(gm) - beta_continuous(sm))
    db_mC = asm_db('mC')
    db_hC = asm_db('hC')
    db_max = asm_db('max')

    # 5hmC>5mC frac among calls where both present (locus-level), both groups
    loc_pair = 0; loc_hgt = 0
    for rd in (germ_reads + som_reads):
        for v in rd.values():
            if v[0] is not None and v[1] is not None:
                loc_pair += 1
                if v[1] > v[0]:
                    loc_hgt += 1
    all_pairwise_calls += loc_pair
    all_5hmC_gt_5mC += loc_hgt

    # (C) dead-zone (max-collapsed, both groups)
    both_max = group_calls_flat(germ_reads, 'max') + som_max
    loc_dz = sum(1 for m in both_max if ML_LO <= m <= ML_HI)
    all_deadzone_calls += loc_dz
    all_total_calls += len(both_max)
    deadzone_locus = loc_dz / len(both_max) if both_max else np.nan

    # (D) subclone proxy: per-read mean meth (somatic group), GMM bimodality
    pr_beta = per_read_mean_meth(som_reads, 'max')
    bim = gmm_bimodality(pr_beta)

    rec.update(
        beta_binarized=round(beta_i, 4), beta_continuous=round(beta_ii, 4),
        beta_03=round(beta_03, 4), beta_07=round(beta_07, 4),
        frac_intermediate_reads=round(frac_interm, 4),
        db_5mC=None if np.isnan(db_mC) else round(db_mC, 4),
        db_5hmC=None if np.isnan(db_hC) else round(db_hC, 4),
        db_max=None if np.isnan(db_max) else round(db_max, 4),
        loc_5hmC_gt_5mC_frac=round(loc_hgt / loc_pair, 4) if loc_pair else None,
        loc_pairwise_calls=loc_pair,
        deadzone_frac_locus=round(deadzone_locus, 4) if both_max else None,
        n_calls_max=len(both_max),
        per_read_beta_mean=round(float(pr_beta.mean()), 4) if len(pr_beta) else None,
        per_read_beta_std=round(float(pr_beta.std()), 4) if len(pr_beta) else None,
        bimodal=bim,
    )
    per_locus.append(rec)
    if (i + 1) % 50 == 0:
        print(f"  [{i+1}/{len(specs)}] processed")

# ---------- aggregate ----------
valid = [r for r in per_locus if 'skipped' not in r]
print(f"[done] {len(valid)}/{len(per_locus)} loci passed read-count filter")

# (A) binarized vs continuous corr
bi = np.array([r['beta_binarized'] for r in valid])
co = np.array([r['beta_continuous'] for r in valid])
mask = ~(np.isnan(bi) | np.isnan(co))
if mask.sum() >= 3:
    r_pear, p_pear = pearsonr(bi[mask], co[mask])
    r_spear, _ = spearmanr(bi[mask], co[mask])
else:
    r_pear = p_pear = r_spear = np.nan
# divergence: |binarized - continuous|, where high frac_intermediate expected
for r in valid:
    r['_binconti_div'] = abs(r['beta_binarized'] - r['beta_continuous']) \
        if not (np.isnan(r['beta_binarized']) or np.isnan(r['beta_continuous'])) else np.nan
top_div = sorted([r for r in valid if not np.isnan(r['_binconti_div'])],
                 key=lambda r: -r['_binconti_div'])[:10]

# (B) 5mC vs 5hmC corr (db_max vs db_5mC, db_max vs db_5hmC, db_5mC vs db_5hmC)
def arr(key):
    return np.array([r[key] if r[key] is not None else np.nan for r in valid])
dmC, dhC, dmax = arr('db_5mC'), arr('db_5hmC'), arr('db_max')
def safecorr(a, b):
    m = ~(np.isnan(a) | np.isnan(b))
    if m.sum() < 3:
        return None
    return round(float(pearsonr(a[m], b[m])[0]), 4)
corr_max_mC = safecorr(dmax, dmC)
corr_max_hC = safecorr(dmax, dhC)
corr_mC_hC = safecorr(dmC, dhC)
frac_5hmC_gt_5mC = all_5hmC_gt_5mC / all_pairwise_calls if all_pairwise_calls else np.nan

# (C) dead-zone overall
deadzone_frac = all_deadzone_calls / all_total_calls if all_total_calls else np.nan

# (D) bimodality
bim_loci = [r for r in valid if r['bimodal']['is_bimodal']]
testable = [r for r in valid if r['bimodal'].get('n', 0) >= 12 and r['bimodal'].get('bic1') is not None]
n_bim = len(bim_loci)
frac_bim = n_bim / len(testable) if testable else np.nan

def xtab(keyfn):
    d = {}
    for r in testable:
        k = keyfn(r)
        d.setdefault(k, {'n': 0, 'bimodal': 0})
        d[k]['n'] += 1
        if r['bimodal']['is_bimodal']:
            d[k]['bimodal'] += 1
    return {str(k): {'n': v['n'], 'bimodal': v['bimodal'],
                     'frac': round(v['bimodal'] / v['n'], 4) if v['n'] else None}
            for k, v in d.items()}

bimodal_by_tp = xtab(lambda r: f"tp={r['is_tp']}")
bimodal_by_cn = xtab(lambda r: r['cn_class'])

# correlation between binarized somatic beta and continuous, restated as the key return metric
binarized_vs_continuous_corr = round(float(r_pear), 4) if not np.isnan(r_pear) else None

# ---------- recommendation ----------
recs = []
if binarized_vs_continuous_corr is not None and binarized_vs_continuous_corr >= 0.97:
    recs.append(f"Binarized P>=0.5 tracks continuous mean-ML tightly (r={binarized_vs_continuous_corr}); "
                "binarization is NOT introducing material bias at the per-locus beta level.")
elif binarized_vs_continuous_corr is not None:
    recs.append(f"Binarized vs continuous corr r={binarized_vs_continuous_corr} is only moderate; "
                "report continuous mean-ML alongside binarized beta.")
if not np.isnan(deadzone_frac):
    recs.append(f"Dead-zone (ML in [50,200]) = {deadzone_frac:.1%} of read-CpG calls; "
                + ("LOW -> calls are confidently meth/unmeth, binarization safe."
                   if deadzone_frac < 0.15 else
                   "NON-trivial -> intermediate calls exist; a 3-state (meth/unmeth/ambiguous) or continuous "
                   "model is preferable for borderline loci."))
if not np.isnan(frac_5hmC_gt_5mC):
    recs.append(f"5hmC>5mC in {frac_5hmC_gt_5mC:.1%} of dual-called CpGs; "
                + (f"MAX-collapse is dominated by 5mC (corr db_max~db_5mC={corr_max_mC}); "
                   "5hmC contributes little signal, mostly low-ML noise."
                   if (corr_max_mC or 0) >= 0.9 else
                   "5hmC materially shifts max-collapse; report 5mC-only as primary, 5hmC as secondary."))
if not np.isnan(frac_bim):
    recs.append(f"Per-read methylation bimodal in {n_bim}/{len(testable)} testable loci ({frac_bim:.1%}); "
                + ("subclone/allelic mixture structure is COMMON -> per-read modality should be reported, "
                   "single-beta summary hides it." if frac_bim >= 0.15 else
                   "most loci are unimodal -> single-beta summary is adequate for the majority; "
                   "flag the bimodal minority for per-read inspection."))
measurement_recommendation = " ".join(recs)

summary = dict(
    n_loci_sampled=len(specs),
    n_loci_valid=len(valid),
    n_canonical=len([s for s in specs if s['is_canonical']]),
    binarized_vs_continuous_corr=binarized_vs_continuous_corr,
    binarized_vs_continuous_spearman=round(float(r_spear), 4) if not np.isnan(r_spear) else None,
    corr_db_max_vs_5mC=corr_max_mC,
    corr_db_max_vs_5hmC=corr_max_hC,
    corr_db_5mC_vs_5hmC=corr_mC_hC,
    frac_5hmC_gt_5mC=round(float(frac_5hmC_gt_5mC), 4) if not np.isnan(frac_5hmC_gt_5mC) else None,
    total_dual_called_cpg=all_pairwise_calls,
    deadzone_frac=round(float(deadzone_frac), 4) if not np.isnan(deadzone_frac) else None,
    deadzone_calls=all_deadzone_calls, total_max_calls=all_total_calls,
    n_loci_testable_bimodal=len(testable),
    n_loci_bimodal=n_bim,
    frac_loci_bimodal=round(float(frac_bim), 4) if not np.isnan(frac_bim) else None,
    bimodal_by_tp=bimodal_by_tp,
    bimodal_by_cn=bimodal_by_cn,
    top_divergent_loci=[dict(chrom=r['chrom'], pos=r['pos'], label=r['label'],
                             beta_bin=r['beta_binarized'], beta_cont=r['beta_continuous'],
                             div=round(r['_binconti_div'], 4),
                             frac_intermediate=r['frac_intermediate_reads'])
                        for r in top_div],
    measurement_recommendation=measurement_recommendation,
)

out = dict(meta=dict(script="49_meth_measurement_audit.py", sample="HCC1395",
                     task_type="A pilot", bam=BAM, window=WINDOW,
                     ml_bin=ML_BIN, deadzone=[ML_LO, ML_HI],
                     bimodal_test="GMM 1-vs-2 BIC (diptest unavailable); "
                                  "bimodal if dBIC>=10 & min_weight>=0.10 & mean_sep>=0.10",
                     note="HP-axis = somatic-controlled (germ HP-group vs somatic HP-subgroup)"),
           summary=summary, per_locus=per_locus)

os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
with open(OUT_JSON, "w") as f:
    json.dump(out, f, indent=2, default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
print(f"[json] {OUT_JSON}")

# ---------- figure ----------
os.makedirs(os.path.dirname(OUT_FIG), exist_ok=True)
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# (1) binarized vs continuous scatter
ax = axes[0, 0]
ax.scatter(co[mask], bi[mask], s=18, alpha=0.5, c='#2c7fb8')
ax.plot([0, 1], [0, 1], 'r--', lw=1)
ax.set_xlabel("continuous beta = mean(ML/255)")
ax.set_ylabel("binarized beta P>=0.5 (ML>=128)")
ax.set_title(f"(A) Binarized vs continuous  r={binarized_vs_continuous_corr}")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

# (2) 5mC vs 5hmC |db| scatter
ax = axes[0, 1]
m2 = ~(np.isnan(dmC) | np.isnan(dhC))
ax.scatter(dmC[m2], dhC[m2], s=18, alpha=0.5, c='#d95f0e')
mx = max(0.05, np.nanmax(np.concatenate([dmC[m2], dhC[m2]])) if m2.sum() else 0.05)
ax.plot([0, mx], [0, mx], 'k--', lw=1)
ax.set_xlabel("ASM |db| 5mC-only")
ax.set_ylabel("ASM |db| 5hmC-only")
ax.set_title(f"(B) 5mC vs 5hmC ASM  corr={corr_mC_hC}\nfrac 5hmC>5mC={frac_5hmC_gt_5mC:.1%}")

# (3) dead-zone histogram of ALL max-collapsed ML across all loci
ax = axes[0, 2]
all_ml = []
for r in valid:
    g = collect_loci(r['chrom'], r['pos'], [r['gg'], r['gs']])
    for grp in (r['gg'], r['gs']):
        all_ml += group_calls_flat(g[grp], 'max')
all_ml = np.array(all_ml)
ax.hist(all_ml, bins=51, color='#756bb1', alpha=0.85)
ax.axvspan(ML_LO, ML_HI, color='red', alpha=0.12)
ax.axvline(ML_BIN, color='k', ls='--', lw=1)
ax.set_xlabel("max-collapsed ML (0-255)")
ax.set_ylabel("read-CpG call count")
ax.set_title(f"(C) Dead-zone [50,200] = {deadzone_frac:.1%}\n(red band = ambiguous; dashed = P=0.5)")

# (4-6) per-read beta distributions for 3 example loci (incl BRCA2 + a bimodal + a unimodal)
example_specs = []
brca2 = next((r for r in valid if r.get('label') == 'BRCA2'), None)
if brca2:
    example_specs.append(brca2)
# add the strongest bimodal (largest dBIC) and a clear unimodal
bim_sorted = sorted([r for r in testable if r['bimodal']['is_bimodal']],
                    key=lambda r: -(r['bimodal']['delta_bic'] or 0))
if bim_sorted:
    example_specs.append(bim_sorted[0])
uni_sorted = sorted([r for r in testable if not r['bimodal']['is_bimodal']
                     and r['bimodal'].get('n', 0) >= 25],
                    key=lambda r: r['bimodal'].get('delta_bic', 0))
if uni_sorted:
    example_specs.append(uni_sorted[0])
example_specs = example_specs[:3]
while len(example_specs) < 3 and len(valid) > len(example_specs):
    for r in valid:
        if r not in example_specs:
            example_specs.append(r); break

panel_ax = [axes[1, 0], axes[1, 1], axes[1, 2]]
for k, r in enumerate(example_specs):
    ax = panel_ax[k]
    g = collect_loci(r['chrom'], r['pos'], [r['gg'], r['gs']])
    pr = per_read_mean_meth(g[r['gs']], 'max')
    ax.hist(pr, bins=20, range=(0, 1), color='#31a354', alpha=0.8)
    bim = r['bimodal']
    tag = "BIMODAL" if bim['is_bimodal'] else "unimodal"
    lab = r.get('label') or f"{r['chrom']}:{r['pos']}"
    ax.set_title(f"(D) {lab}  {tag}\nn_reads={bim.get('n')} dBIC={bim.get('delta_bic')}"
                 f" sep={bim.get('mean_sep')}", fontsize=9)
    ax.set_xlabel("per-read mean methylation (somatic HP grp)")
    ax.set_ylabel("read count")
    if bim['is_bimodal'] and bim.get('comp_means'):
        for cm in bim['comp_means']:
            ax.axvline(cm, color='red', ls='--', lw=1)

fig.suptitle("M2 - Methylation MEASUREMENT audit + subclone proxy (HCC1395, HP-axis, single sample)",
             fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT_FIG, dpi=130)
print(f"[fig] {OUT_FIG}")

# ---------- console summary ----------
print("\n" + "=" * 70)
print("M2 SUMMARY")
print("=" * 70)
print(f"binarized_vs_continuous_corr : {binarized_vs_continuous_corr}")
print(f"frac_5hmC_gt_5mC             : {summary['frac_5hmC_gt_5mC']}  (corr db_max~5mC={corr_max_mC})")
print(f"deadzone_frac                : {summary['deadzone_frac']}")
print(f"n_loci_bimodal / testable    : {n_bim} / {len(testable)}  (frac={summary['frac_loci_bimodal']})")
print(f"bimodal_by_tp                : {bimodal_by_tp}")
print(f"bimodal_by_cn                : {bimodal_by_cn}")
print(f"RECOMMENDATION: {measurement_recommendation}")
