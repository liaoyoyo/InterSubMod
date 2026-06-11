#!/usr/bin/env python3
"""
38 — G2 LOH diagnostic classifier (per-locus). User-chosen 診斷分類 (NOT forced reassign).

Premise: longphase-S does NOT validate LOH (unlike longphase-TO V6). In LOH regions
(physically one haplotype), longphase-S still emits HP1 vs HP1-1 splits driven by the
somatic variant -> "apparent two-haplotype". This classifies each LOH-region somatic
locus into the CAUSE of its apparent HP1/HP1-1 methylation split:

  (a) self_phasing_artifact : the tag split exists but methylation does NOT form a real
       cluster (blind-ARI < 0.30) OR the split is a length/spanning collider (placebo high).
       = same physical haplotype split by the variant; no real subclone methylation.
  (b) CN_regression         : real-ish |Δβ| but driven by low coverage + extreme baseline
       (regression-to-extreme), modulated by CN class (cnLOH/gainLOH/lossLOH). NOT a clean
       somatic subclone signal.
  (c) candidate_subclone    : blind-ARI>=0.30 AND placebo<0.10 AND mid-baseline AND
       decent coverage -> a real methylation cluster on the somatic axis even inside LOH.

Cross with CN class (SEQC2 BED set-ops, obs25 logic) + baseline extremity + coverage.
Stratified sample of sig LOH-regime HP-axis loci (cap per CN class) for tractable compute.
Picks one exemplar per class for the IGV triptych (A6 fig 3).

Single-sample HCC1395. Output: genome_survey_v2/loh_diagnostic_classifier.json
"""
import json, random, csv
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
import pysam
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.metrics import adjusted_rand_score

random.seed(42); np.random.seed(42)
WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
TPP = WS / "genome_survey_v2/asm_dualaxis_tp.tsv"
CNV_BED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
OUT = WS / "genome_survey_v2/loh_diagnostic_classifier.json"
HP_AXES = {"HP1_vs_HP1-1", "HP2_vs_HP2-1"}
WINDOW = 600; ML_HIGH, ML_LOW = 200, 50
MIN_CORE_CPG = 4; MIN_SHARED = 3; MIN_GROUP = 8
CAP_PER_CN = 55          # cap loci per CN class for tractable compute
ARI_PASS = 0.30; PLACEBO_MAX = 0.10; HI_COV = 100; EXTREME = 0.3
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


def collect(chrom, var0, groups):
    s, e = var0 - WINDOW, var0 + WINDOW
    g = {k: [] for k in groups}; lens = {k: [] for k in groups}
    for r in bam.fetch(chrom, max(0, s), e):
        if r.flag & 0x900 or r.flag & 0x400:
            continue
        h = hp(r)
        if h not in groups:
            continue
        m = {p: st for p, st in mcalls(r).items() if s <= p <= e and st >= 0}
        if len(m) >= MIN_CORE_CPG:
            g[h].append(m); lens[h].append(r.query_length or (r.reference_end - r.reference_start))
    return g, lens


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


def blind_ari(g0, g1):
    reads = g0 + g1; lab0 = [0] * len(g0) + [1] * len(g1)
    D, keep = observed_only_distance(reads)
    if D is None:
        return None, None
    lab = np.array([lab0[i] for i in keep])
    if lab.sum() < MIN_GROUP or (len(lab) - lab.sum()) < MIN_GROUP:
        return None, None
    preds = []
    for link in ('average', 'complete'):
        preds.append(AgglomerativeClustering(n_clusters=2, metric='precomputed', linkage=link).fit_predict(D))
    try:
        preds.append(SpectralClustering(n_clusters=2, affinity='precomputed', random_state=42).fit_predict(np.exp(-D / (D.std() + 1e-9))))
    except Exception:
        pass
    return max(adjusted_rand_score(lab, p) for p in preds), (int((lab == 0).sum()), int((lab == 1).sum()))


def placebo(g_germ, lens_germ):
    if len(g_germ) < 2 * MIN_GROUP:
        return None
    med = np.median(lens_germ)
    a = [g_germ[i] for i in range(len(g_germ)) if lens_germ[i] >= med]
    b = [g_germ[i] for i in range(len(g_germ)) if lens_germ[i] < med]
    if len(a) < MIN_GROUP or len(b) < MIN_GROUP:
        return None
    ari, _ = blind_ari(b, a)
    return ari


def load_bed_types(path):
    by_type = {"gain": defaultdict(list), "loss": defaultdict(list), "loh": defaultdict(list)}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[3] in by_type:
                by_type[p[3]][p[0]].append((int(p[1]), int(p[2])))
    return by_type


def in_set(regs, chrom, pos):
    for s, e in regs.get(chrom, []):
        if s <= pos <= e:
            return True
    return False


def cn_class(bt, chrom, pos):
    g = in_set(bt["gain"], chrom, pos); l = in_set(bt["loss"], chrom, pos)
    if g:
        return "gainLOH"
    if l:
        return "lossLOH"
    return "cnLOH"   # loh AND not gain AND not loss


def classify(ari, plc, n_cpg, extremity):
    if ari is None:
        return "uneval"
    collider = (plc is not None and plc >= ARI_PASS)
    if ari < ARI_PASS or collider:
        return "self_phasing_artifact"
    if extremity > EXTREME and n_cpg < HI_COV:
        return "CN_regression"
    if ari >= ARI_PASS and (plc is None or plc < PLACEBO_MAX):
        return "candidate_subclone"
    return "intermediate"


def main():
    bt = load_bed_types(CNV_BED)
    # load sig LOH-regime HP-axis loci
    loh_loci = []
    with open(TPP) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            if row['axis'] not in HP_AXES:
                continue
            if row['loh_status'] != "LOH":
                continue
            try:
                wp = float(row['wilcoxon_p'])
            except Exception:
                continue
            if wp >= 0.05:
                continue
            loh_loci.append(dict(chrom=row['chrom'], pos=int(row['somatic_pos']), axis=row['axis'],
                                 n_cpg=int(row['n_paired_cpg']), germ_beta=float(row['mean_germ_beta']),
                                 delta=float(row['mean_delta']), wp=wp))
    # assign CN class, stratified cap
    for r in loh_loci:
        r['cn_class'] = cn_class(bt, r['chrom'], r['pos'])
        r['extremity'] = abs(r['germ_beta'] - 0.5)
    by_cn = defaultdict(list)
    for r in loh_loci:
        by_cn[r['cn_class']].append(r)
    sample = []
    for cn, lst in by_cn.items():
        random.shuffle(lst)
        sample.extend(lst[:CAP_PER_CN])
    print(f"sig LOH HP-axis loci total={len(loh_loci)} (cnLOH={len(by_cn['cnLOH'])}, "
          f"gainLOH={len(by_cn['gainLOH'])}, lossLOH={len(by_cn['lossLOH'])}); evaluating sample n={len(sample)}")

    rows = []
    for r in sample:
        var0 = r['pos'] - 1
        gg, gs = ("1", "1-1") if r['axis'] == "HP1_vs_HP1-1" else ("2", "2-1")
        try:
            g, lens = collect(r['chrom'], var0, [gg, gs])
        except Exception:
            continue
        if len(g[gg]) < MIN_GROUP or len(g[gs]) < MIN_GROUP:
            r2 = dict(r); r2.update(ari=None, placebo=None, klass="uneval", n_germ=len(g[gg]), n_som=len(g[gs]))
            rows.append(r2); continue
        ari, ns = blind_ari(g[gg], g[gs])
        plc = placebo(g[gg], lens[gg])
        klass = classify(ari, plc, r['n_cpg'], r['extremity'])
        r2 = dict(r)
        r2.update(ari=(None if ari is None else round(ari, 3)),
                  placebo=(None if plc is None else round(plc, 3)),
                  klass=klass, n_germ=(ns[0] if ns else len(g[gg])), n_som=(ns[1] if ns else len(g[gs])))
        rows.append(r2)

    klass_ct = Counter(x['klass'] for x in rows)
    # class x CN cross-tab
    cross = defaultdict(lambda: Counter())
    for x in rows:
        cross[x['klass']][x['cn_class']] += 1
    # exemplars for IGV triptych: best candidate_subclone (max ari, low placebo),
    # clearest self_phasing_artifact (low ari), clearest CN_regression (extreme baseline + lowcov)
    def pick(klass, key, reverse=True):
        cand = [x for x in rows if x['klass'] == klass and x['ari'] is not None]
        if not cand:
            return None
        return sorted(cand, key=key, reverse=reverse)[0]
    exemplars = {
        "candidate_subclone": pick("candidate_subclone", lambda x: (x['ari'], -(x['placebo'] or 0))),
        "self_phasing_artifact": pick("self_phasing_artifact", lambda x: -x['ari']),  # lowest ari
        "CN_regression": pick("CN_regression", lambda x: (x['extremity'], -x['n_cpg'])),
    }
    summary = dict(
        n_sig_loh_total=len(loh_loci),
        cn_class_totals={k: len(v) for k, v in by_cn.items()},
        n_evaluated=len([x for x in rows if x['klass'] != 'uneval']),
        class_counts=dict(klass_ct),
        class_by_cn={k: dict(v) for k, v in cross.items()},
        exemplars={k: (None if v is None else dict(locus=f"{v['chrom']}:{v['pos']}", axis=v['axis'],
                                                   ari=v['ari'], placebo=v['placebo'], cn_class=v['cn_class'],
                                                   n_cpg=v['n_cpg'], germ_beta=v['germ_beta'], delta=v['delta']))
                   for k, v in exemplars.items()},
        rule=f"self_phasing: ari<{ARI_PASS} OR placebo>={ARI_PASS}; CN_regression: extremity>{EXTREME} & n_cpg<{HI_COV}; candidate_subclone: ari>={ARI_PASS} & placebo<{PLACEBO_MAX} & mid-baseline & hi-cov",
        sample="HCC1395 single-sample; longphase-S HP-axis; LOH-regime sig loci (stratified sample)",
    )
    json.dump(dict(summary=summary, loci=rows), open(OUT, "w"), indent=1)
    print("\n=== LOH DIAGNOSTIC CLASSIFIER SUMMARY ===")
    print(json.dumps(summary, indent=1, ensure_ascii=False))
    print(f"saved -> {OUT}")


if __name__ == "__main__":
    main()
