#!/usr/bin/env python3
"""
30 — Regime-A credible-subset probe (G1 framing chosen by user 2026-06-02).

Question (decisive): do the CLEANEST candidate ASM loci — obs23 regime-A
(n_paired_cpg>=100 AND |germ_beta-0.5|<=0.3 AND nonLOH, HP-axis, sig p<0.05),
i.e. exactly the conditions that AVOID regression-to-extreme — show a REAL
somatic methylation cluster (HP1 vs HP1-1) that survives the collider gate?

Pipeline per locus (HP-axis = somatic-controlled, the right axis):
  M0  observed-only Hamming distance over shared CpG (NO median-impute;
      iterative drop poorly-connected reads; residual NaN -> conservative nanmax).
      [reuses proven extraction from 24_cluster_eval_core / 25_b2_broad_scan]
  M1  blind-ARI/AMI: k=2 Agglomerative(average,complete)+Spectral, max ARI vs true split (PRIMARY).
  M8  placebo-anchor (collider gate): split the germline-tag (HP1) reads by read-length
      median (mimics the spanning/length selection that DEFINES HP1-1) -> blind-ARI.
      If placebo ARI is high, the "cluster" is a length/spanning collider, not methylation.
  aux R2 (skbio permanova), silhouette (supervised; reference only).
Aggregate:
  M5  coverage-effect signature: Spearman(ARI, log n_paired_cpg) + (ARI, log total_reads).
      Monotone-positive = artifact signature.
  pass rule (Tier-A site candidate): ARI>=0.30 AND placebo_ARI<0.10.

Single-sample HCC1395 paired_full. Tier-A ceiling (M9 cross-sample not available).
Output: genome_survey_v2/regimeA_credible_probe.json (+ printed summary).
"""
import json, random
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
from collections import Counter
from scipy.stats import spearmanr, wilcoxon
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score
from skbio.stats.distance import DistanceMatrix, permanova

random.seed(42); np.random.seed(42)

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
TPP = WS / "genome_survey_v2/asm_dualaxis_tp.tsv"
OUT = WS / "genome_survey_v2/regimeA_credible_probe.json"
HP_AXES = {"HP1_vs_HP1-1", "HP2_vs_HP2-1"}
WINDOW = 600; ML_HIGH, ML_LOW = 200, 50
MIN_CORE_CPG = 4; MIN_SHARED = 3; MIN_GROUP = 8; N_PERM = 499
ARI_PASS = 0.30; PLACEBO_MAX = 0.10

bam = pysam.AlignmentFile(BAM, "rb")


def hp(r):
    for t, v in r.tags:
        if t == "HP":
            return str(v)
    return None


def mcalls(r):
    o = {}
    try:
        mod = r.modified_bases
    except Exception:
        return o
    if not mod:
        return o
    r2 = {a: b for a, b in r.get_aligned_pairs(matches_only=False) if a is not None and b is not None}
    for k, calls in mod.items():
        if k[2] != 'm':
            continue
        for rp, ml in calls:
            rf = r2.get(rp)
            if rf is not None:
                o[rf] = 1 if ml >= ML_HIGH else (0 if ml <= ML_LOW else -1)
    return o


def collect(chrom, pos, groups):
    var0 = pos - 1; s, e = var0 - WINDOW, var0 + WINDOW
    g = {k: [] for k in groups}; meta = {k: [] for k in groups}
    for r in bam.fetch(chrom, max(0, s), e):
        if r.flag & 0x900 or r.flag & 0x400:
            continue
        h = hp(r)
        if h not in groups:
            continue
        m = {p: st for p, st in mcalls(r).items() if s <= p <= e and st >= 0}
        if len(m) >= MIN_CORE_CPG:
            g[h].append(m)
            meta[h].append(r.query_length or (r.reference_end - r.reference_start))
    return g, meta


def observed_only_distance(reads):
    """M0: Hamming over shared CpG, NO median-impute. Drop poorly-connected reads iteratively."""
    n = len(reads)
    cov = Counter()
    for r in reads:
        for c in r:
            cov[c] += 1
    core = {c for c, k in cov.items() if k >= max(3, 0.3 * n)}
    reads = [{c: r[c] for c in r if c in core} for r in reads]
    keep = [i for i, r in enumerate(reads) if len(r) >= MIN_CORE_CPG]
    reads = [reads[i] for i in keep]; n = len(reads)
    if n < 2 * MIN_GROUP:
        return None, keep
    while n > 2 * MIN_GROUP:
        bad = np.zeros(n)
        for i in range(n):
            inc = sum(1 for j in range(n) if i != j and len(set(reads[i]) & set(reads[j])) < MIN_SHARED)
            bad[i] = inc / (n - 1)
        if bad.max() <= 0.4:
            break
        d = int(bad.argmax()); reads.pop(d); keep.pop(d); n -= 1
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            sh = set(reads[i]) & set(reads[j])
            D[i, j] = D[j, i] = (np.mean([reads[i][c] != reads[j][c] for c in sh]) if len(sh) >= MIN_SHARED else np.nan)
    if np.isnan(D).any():
        D[np.isnan(D)] = np.nanmax(D) if not np.isnan(D).all() else 0.5
    return D, keep


def blind_ari(D, lab):
    preds = []
    for link in ('average', 'complete'):
        preds.append(AgglomerativeClustering(n_clusters=2, metric='precomputed', linkage=link).fit_predict(D))
    try:
        preds.append(SpectralClustering(n_clusters=2, affinity='precomputed', random_state=42).fit_predict(np.exp(-D / (D.std() + 1e-9))))
    except Exception:
        pass
    return (max(adjusted_rand_score(lab, p) for p in preds),
            max(adjusted_mutual_info_score(lab, p) for p in preds))


def eval_split(g0, g1):
    reads = g0 + g1; lab0 = [0] * len(g0) + [1] * len(g1)
    D, keep = observed_only_distance(reads)
    if D is None:
        return None
    lab = np.array([lab0[i] for i in keep])
    if lab.sum() < MIN_GROUP or (len(lab) - lab.sum()) < MIN_GROUP:
        return None
    ari, ami = blind_ari(D, lab)
    try:
        dm = DistanceMatrix(D, [str(i) for i in range(len(lab))])
        res = permanova(dm, lab.tolist(), permutations=N_PERM)
        F = res['test statistic']; R2 = F / (F + (len(lab) - 2)); pp = res['p-value']
    except Exception:
        R2, pp = float('nan'), float('nan')
    try:
        sil = silhouette_score(D, lab, metric='precomputed')
    except Exception:
        sil = float('nan')
    return dict(n0=int((lab == 0).sum()), n1=int((lab == 1).sum()), ari=ari, ami=ami, R2=R2, p=pp, sil=sil)


def placebo(g_germ, lens_germ):
    """M8 collider gate: split germline reads by length median (mimics ALT-spanning selection)."""
    if len(g_germ) < 2 * MIN_GROUP:
        return None
    med = np.median(lens_germ)
    a = [g_germ[i] for i in range(len(g_germ)) if lens_germ[i] >= med]
    b = [g_germ[i] for i in range(len(g_germ)) if lens_germ[i] < med]
    if len(a) < MIN_GROUP or len(b) < MIN_GROUP:
        return None
    return eval_split(b, a)


def main():
    d = pd.read_csv(TPP, sep="\t").dropna(subset=["wilcoxon_p", "mean_delta"])
    d = d[d.axis.isin(HP_AXES)].copy()
    d["extremity"] = (d.mean_germ_beta - 0.5).abs()
    regimeA = d[(d.wilcoxon_p < 0.05) & (d.n_paired_cpg >= 100) &
                (d.extremity <= 0.3) & (d.loh_status != "LOH")].copy()
    print(f"regime-A credible loci (recomputed): n={len(regimeA)}")

    rows = []
    for _, r in regimeA.iterrows():
        chrom, pos, axis = r.chrom, int(r.somatic_pos), r.axis
        gg, gs = ("1", "1-1") if axis == "HP1_vs_HP1-1" else ("2", "2-1")
        try:
            g, meta = collect(chrom, pos, [gg, gs])
        except Exception as ex:
            continue
        if len(g[gg]) < MIN_GROUP or len(g[gs]) < MIN_GROUP:
            continue
        res = eval_split(g[gg], g[gs])
        if not res:
            continue
        pl = placebo(g[gg], meta[gg])
        placebo_ari = pl['ari'] if pl else float('nan')
        is_pass = (res['ari'] >= ARI_PASS) and (np.isnan(placebo_ari) or placebo_ari < PLACEBO_MAX)
        collider = (not np.isnan(placebo_ari)) and (placebo_ari >= ARI_PASS)
        rows.append(dict(locus=f"{chrom}:{pos}", axis=axis,
                         n_cpg_survey=int(r.n_paired_cpg), germ_beta=float(r.mean_germ_beta),
                         delta=float(r.mean_delta), wilcoxon_p_survey=float(r.wilcoxon_p),
                         n_germ=res['n0'], n_som=res['n1'], ari=res['ari'], ami=res['ami'],
                         R2=res['R2'], perm_p=res['p'], sil=res['sil'],
                         placebo_ari=(None if np.isnan(placebo_ari) else placebo_ari),
                         pass_tierA=bool(is_pass), collider_flag=bool(collider)))

    n_eval = len(rows)
    n_pass = sum(r['pass_tierA'] for r in rows)
    n_collider = sum(r['collider_flag'] for r in rows)
    aris = [r['ari'] for r in rows]
    placebos = [r['placebo_ari'] for r in rows if r['placebo_ari'] is not None]
    # M5 coverage-effect signature within regime-A
    cov = [r['n_cpg_survey'] for r in rows]
    rho_cov, p_cov = (spearmanr(aris, np.log(cov)) if n_eval >= 5 else (float('nan'), float('nan')))
    nreads = [r['n_germ'] + r['n_som'] for r in rows]
    rho_nr, p_nr = (spearmanr(aris, np.log(nreads)) if n_eval >= 5 else (float('nan'), float('nan')))
    # paired real-ari vs placebo-ari (same loci, where placebo available)
    paired = [(r['ari'], r['placebo_ari']) for r in rows if r['placebo_ari'] is not None]
    if len(paired) >= 6:
        ra = [x[0] for x in paired]; pa = [x[1] for x in paired]
        try:
            w_stat, w_p = wilcoxon(ra, pa, alternative='greater')
        except Exception:
            w_stat, w_p = float('nan'), float('nan')
    else:
        w_stat, w_p = float('nan'), float('nan')

    summary = dict(
        regimeA_n_total=int(len(regimeA)), n_evaluated=n_eval,
        n_pass_tierA=n_pass, pass_rate=(n_pass / n_eval if n_eval else None),
        n_collider_flag=n_collider,
        median_ari=(float(np.median(aris)) if aris else None),
        max_ari=(float(np.max(aris)) if aris else None),
        median_placebo_ari=(float(np.median(placebos)) if placebos else None),
        M5_spearman_ari_vs_logcov=dict(rho=float(rho_cov), p=float(p_cov)),
        M5_spearman_ari_vs_lognreads=dict(rho=float(rho_nr), p=float(p_nr)),
        real_vs_placebo_wilcoxon_greater=dict(stat=float(w_stat), p=float(w_p), n_pairs=len(paired)),
        axis="HP-axis (HP1 vs HP1-1 / HP2 vs HP2-1, somatic-controlled)",
        sample="HCC1395 paired_full (single-sample; Tier-A ceiling, M9 cross-sample N/A)",
        pass_rule=f"ARI>={ARI_PASS} AND placebo_ARI<{PLACEBO_MAX}",
    )
    json.dump(dict(summary=summary, loci=rows), open(OUT, "w"), indent=1)

    print("\n=== REGIME-A CREDIBLE PROBE SUMMARY ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print(f"\n=== per-locus (sorted by ARI desc), top 20 ===")
    print(f"{'locus':<18}{'axis':<16}{'ng/ns':>9}{'ARI':>7}{'placebo':>9}{'R2':>7}{'pp':>7}{'sil':>7}  verdict")
    for r in sorted(rows, key=lambda x: -x['ari'])[:20]:
        v = "Tier-A候選" if r['pass_tierA'] else ("COLLIDER" if r['collider_flag'] else "弱ARI")
        pa = f"{r['placebo_ari']:.3f}" if r['placebo_ari'] is not None else "NA"
        print(f"{r['locus']:<18}{r['axis']:<16}{str(r['n_germ'])+'/'+str(r['n_som']):>9}"
              f"{r['ari']:>7.3f}{pa:>9}{r['R2']:>7.3f}{r['perm_p']:>7.3f}{r['sil']:>7.3f}  {v}")
    print(f"\nsaved -> {OUT}")


if __name__ == "__main__":
    main()
