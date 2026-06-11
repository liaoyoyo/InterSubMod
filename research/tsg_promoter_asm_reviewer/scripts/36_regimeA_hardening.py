#!/usr/bin/env python3
"""
36 — Regime-A hardening: M2 matched germline-het null + M3 rarefied silhouette.

Context: 30_regimeA_credible_probe.py found regime-A (cleanest subset, HP-axis)
shows 23/62 loci pass Tier-A cut (ARI>=0.30 & placebo<0.10), real>>placebo p=5e-8,
no coverage signature (M5 rho=-0.18 NS). The DECISIVE remaining gate is M2:
is the somatic-HP-axis clustering ABOVE the germline-het (baseline-allelic) null?
This directly addresses the 5/31 concern "ASM = baseline allelic methylation".

M2: build germline-het null (HP1 vs HP2 at het SNPs, = baseline allelic, NO somatic),
    caliper-matched to regime-A coverage (total reads), n>=50; same observed-only +
    blind-ARI pipeline. Compare regime-A ARI distribution vs het-null:
    Mann-Whitney(greater) + Cliff's delta + bootstrap CI of median diff.
M3: rarefied silhouette — downsample both groups to min size, B=200, for regime-A
    pass loci vs het-null; confirms clustering not a read-count artifact (O11 killer).

Loads regime-A ARIs from regimeA_credible_probe.json (no recompute).
Single-sample HCC1395. Tier-A ceiling.
Output: genome_survey_v2/regimeA_hardening.json
"""
import json, random, subprocess
from pathlib import Path
import numpy as np
import pysam
from collections import Counter
from scipy.stats import mannwhitneyu
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.metrics import adjusted_rand_score, silhouette_score

random.seed(42); np.random.seed(42)

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GERMLINE = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
PROBE = WS / "genome_survey_v2/regimeA_credible_probe.json"
OUT = WS / "genome_survey_v2/regimeA_hardening.json"
WINDOW = 600; ML_HIGH, ML_LOW = 200, 50
MIN_CORE_CPG = 4; MIN_SHARED = 3; MIN_GROUP = 8
N_HET_TARGET = 60; RAREFY_B = 200

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
    g = {k: [] for k in groups}
    for r in bam.fetch(chrom, max(0, s), e):
        if r.flag & 0x900 or r.flag & 0x400:
            continue
        h = hp(r)
        if h not in groups:
            continue
        m = {p: st for p, st in mcalls(r).items() if s <= p <= e and st >= 0}
        if len(m) >= MIN_CORE_CPG:
            g[h].append(m)
    return g


def observed_only_distance(reads):
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
        return None, keep, None
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
    return D, keep, reads


def blind_ari(D, lab):
    preds = []
    for link in ('average', 'complete'):
        preds.append(AgglomerativeClustering(n_clusters=2, metric='precomputed', linkage=link).fit_predict(D))
    try:
        preds.append(SpectralClustering(n_clusters=2, affinity='precomputed', random_state=42).fit_predict(np.exp(-D / (D.std() + 1e-9))))
    except Exception:
        pass
    return max(adjusted_rand_score(lab, p) for p in preds)


def eval_split(g0, g1, rarefy=False):
    reads = g0 + g1; lab0 = [0] * len(g0) + [1] * len(g1)
    D, keep, kreads = observed_only_distance(reads)
    if D is None:
        return None
    lab = np.array([lab0[i] for i in keep])
    if lab.sum() < MIN_GROUP or (len(lab) - lab.sum()) < MIN_GROUP:
        return None
    ari = blind_ari(D, lab)
    try:
        sil = silhouette_score(D, lab, metric='precomputed')
    except Exception:
        sil = float('nan')
    out = dict(ari=ari, sil=sil, n0=int((lab == 0).sum()), n1=int((lab == 1).sum()))
    if rarefy:
        idx0 = np.where(lab == 0)[0]; idx1 = np.where(lab == 1)[0]
        m = min(len(idx0), len(idx1))
        sils = []
        for _ in range(RAREFY_B):
            s0 = np.random.choice(idx0, m, replace=False); s1 = np.random.choice(idx1, m, replace=False)
            sub = np.concatenate([s0, s1]); subD = D[np.ix_(sub, sub)]
            subl = np.array([0] * m + [1] * m)
            try:
                sils.append(silhouette_score(subD, subl, metric='precomputed'))
            except Exception:
                pass
        out['rarefied_sil_median'] = float(np.median(sils)) if sils else float('nan')
    return out


def cliffs_delta(a, b):
    a = np.asarray(a); b = np.asarray(b)
    gt = sum((x > b).sum() for x in a); lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (len(a) * len(b))


def sample_het_sites(n_target):
    """Sample het SNP sites genome-wide (baseline-allelic null), 1bp REF/ALT."""
    out = subprocess.run(
        f"bcftools view -r chr1,chr2,chr3,chr5,chr7,chr11,chr15,chr17 {GERMLINE} 2>/dev/null | "
        f"bcftools view -i 'GT=\"het\"' 2>/dev/null | "
        f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' 2>/dev/null | "
        f"awk 'length($3)==1&&length($4)==1' | shuf --random-source=<(yes 42) | head -{n_target * 8}",
        shell=True, capture_output=True, text=True, executable="/bin/bash")
    sites = []
    for line in out.stdout.splitlines():
        p = line.split("\t")
        sites.append((p[0], int(p[1])))
    return sites


def main():
    probe = json.load(open(PROBE))
    regimeA_loci = probe['loci']
    regimeA_ari = [r['ari'] for r in regimeA_loci]
    # coverage floor for caliper-match: 25th pct of regime-A total reads
    regimeA_totreads = [r['n_germ'] + r['n_som'] for r in regimeA_loci]
    cov_floor = int(np.percentile(regimeA_totreads, 25))
    print(f"regime-A: n={len(regimeA_loci)}, median ARI={np.median(regimeA_ari):.3f}, "
          f"total-reads 25pct(cov_floor)={cov_floor}")

    # ---- M2: matched germline-het null (HP1 vs HP2, baseline allelic) ----
    het_sites = sample_het_sites(N_HET_TARGET)
    print(f"sampled {len(het_sites)} candidate het sites; evaluating (cov-matched, target n>={N_HET_TARGET})...")
    null_rows = []
    for chrom, pos in het_sites:
        if len(null_rows) >= N_HET_TARGET:
            break
        try:
            g = collect(chrom, pos, ["1", "2"])
        except Exception:
            continue
        if len(g["1"]) < MIN_GROUP or len(g["2"]) < MIN_GROUP:
            continue
        if (len(g["1"]) + len(g["2"])) < cov_floor:   # caliper-match coverage to regime-A
            continue
        r = eval_split(g["1"], g["2"], rarefy=True)
        if not r:
            continue
        null_rows.append(dict(locus=f"{chrom}:{pos}", ari=r['ari'], sil=r['sil'],
                              rarefied_sil=r.get('rarefied_sil_median'), n0=r['n0'], n1=r['n1']))
    null_ari = [r['ari'] for r in null_rows]
    print(f"het-null evaluated: n={len(null_rows)}, median ARI={np.median(null_ari):.3f}" if null_ari else "het-null: 0")

    # ---- M3 rarefied for regime-A pass loci (recompute rarefied sil) ----
    # (regime-A pass loci rarefied silhouette — recompute on the fly for pass loci)
    pass_loci = [r for r in regimeA_loci if r['pass_tierA']]
    print(f"regime-A pass loci (Tier-A cut): n={len(pass_loci)}; computing rarefied silhouette...")
    rare_rows = []
    for r in pass_loci:
        chrom, pos = r['locus'].split(":"); pos = int(pos)
        gg, gs = ("1", "1-1") if r['axis'] == "HP1_vs_HP1-1" else ("2", "2-1")
        try:
            g = collect(chrom, pos, [gg, gs])
        except Exception:
            continue
        if len(g[gg]) < MIN_GROUP or len(g[gs]) < MIN_GROUP:
            continue
        rr = eval_split(g[gg], g[gs], rarefy=True)
        if rr:
            rare_rows.append(dict(locus=r['locus'], rarefied_sil=rr.get('rarefied_sil_median'), ari=rr['ari']))

    # ---- M2 stats: regime-A vs het-null ----
    stats = {}
    if null_ari and len(regimeA_ari) >= 6:
        u, p = mannwhitneyu(regimeA_ari, null_ari, alternative='greater')
        delta = cliffs_delta(regimeA_ari, null_ari)
        # bootstrap CI of median diff
        diffs = []
        for _ in range(2000):
            sa = np.random.choice(regimeA_ari, len(regimeA_ari), replace=True)
            sn = np.random.choice(null_ari, len(null_ari), replace=True)
            diffs.append(np.median(sa) - np.median(sn))
        ci = (float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5)))
        stats = dict(mw_u=float(u), mw_p_greater=float(p), cliffs_delta=float(delta),
                     median_diff=float(np.median(regimeA_ari) - np.median(null_ari)),
                     median_diff_ci95=ci,
                     regimeA_median_ari=float(np.median(regimeA_ari)),
                     hetnull_median_ari=float(np.median(null_ari)),
                     n_regimeA=len(regimeA_ari), n_hetnull=len(null_ari))

    # also compare pass-rate at ARI>=0.30 between regime-A and het-null
    passrate_regimeA = float(np.mean([a >= 0.30 for a in regimeA_ari]))
    passrate_hetnull = float(np.mean([a >= 0.30 for a in null_ari])) if null_ari else None
    rare_med_pass = float(np.median([x['rarefied_sil'] for x in rare_rows if x['rarefied_sil'] is not None])) if rare_rows else None
    rare_med_null = float(np.median([x['rarefied_sil'] for x in null_rows if x['rarefied_sil'] is not None])) if null_rows else None

    summary = dict(
        M2_regimeA_vs_hetnull=stats,
        passrate_ARI30=dict(regimeA=passrate_regimeA, hetnull=passrate_hetnull),
        M3_rarefied_silhouette=dict(regimeA_pass_median=rare_med_pass, hetnull_median=rare_med_null,
                                    note="rarefied to min group size, B=200; >hetnull = not read-count artifact"),
        cov_floor_total_reads=cov_floor,
        axis="HP-axis somatic-controlled (regimeA) vs germline-het HP1-vs-HP2 baseline-allelic (null)",
        sample="HCC1395 paired_full single-sample; Tier-A ceiling (M9 cross-sample N/A)",
    )
    json.dump(dict(summary=summary, het_null=null_rows, regimeA_rarefied=rare_rows), open(OUT, "w"), indent=1)

    print("\n=== HARDENING SUMMARY ===")
    print(json.dumps(summary, indent=1, ensure_ascii=False))
    print(f"\nsaved -> {OUT}")


if __name__ == "__main__":
    main()
