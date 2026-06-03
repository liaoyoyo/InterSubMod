#!/usr/bin/env python3
"""
62 - Per-sample CREDIBLE ASM loci discovery (replicates the HCC1395 method per sample).

GOAL (user 2026-06-03): apply the SAME pipeline that discovered HCC1395's
credible_loci_annotation.tsv to EACH of the other samples, de novo, to find each
sample's OWN credible ASM loci — and confirm the method transfers cleanly.

FAITHFUL REPLICATION of the HCC1395 5-stage pipeline (Explore-mapped):
  Stage 1-2 (survey, normally MSA 19 + pivot 18): here via VALIDATED pysam
            (crossval Pearson=1.0) — per-CpG paired HP-axis beta + Wilcoxon.
  Stage 3 (regimeA credible probe, 30): blind read-clustering ARI (Hamming +
           agglomerative/spectral k=2) + placebo collider gate. SAME constants.
  Stage 4 (LOH/CN class, 38): HCC1395-ONLY (SEQC2 CNV BED) -> SKIPPED for non-HCC1395
           (cn_class left blank), exactly as Explore recommended. Not a blocker.
  Stage 5 (annotation, 39): RefSeq TSS + UCSC CpG islands (universal reference files).

TRACTABILITY: credible loci require n_paired_cpg>=100 (very CpG-dense regions) so they
  concentrate in/near CpG islands. We PRE-FILTER each sample's TP somatic SNVs to those
  within +/-CPG_PROX of a CpG island, then run the full credible pipeline EXHAUSTIVELY on
  that biologically-relevant universe (random subsetting would miss the ~0.06% credible
  rate). Reported as a funnel so the method's transfer + difficulties are explicit.

SAME thresholds as HCC1395 regimeA_tierA:
  survey: METH_THR=0.5, MIN_N=3/group/CpG, MIN_PAIRED=5
  regimeA filter: wilcoxon_p<0.05, n_paired_cpg>=100 (strict) | >=30 (relaxed tier),
                  extremity=|germ_beta-0.5|<=0.3
  ARI gate: ARI_PASS=0.30, PLACEBO_MAX=0.10 (pass_tierA = ari>=0.30 AND placebo<0.10)

USAGE: python3 62_per_sample_credible_discovery.py <SAMPLE> [SURVEY_WINDOW=1000]
  -> genome_survey_v2/cn_confound/cross_sample/<SAMPLE>_credible_discovery.json
"""
import os, sys, json
os.environ.setdefault("TMPDIR", "/big7_disk")
import gzip
import numpy as np
import pysam
from collections import Counter, defaultdict
from scipy import stats
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
GTF = f"{ROOT}/data/hg38.ncbiRefSeq.gtf.gz"
CPG = f"{ROOT}/data/cpgIslandExt_hg38.txt.gz"

RUNDIR = {
    "HCC1395": "20260314_HCC1395_paired_full_full_complete_matrix",
    "COLO829": "20260315_COLO829_paired_full_full_complete_matrix",
    "H1437":   "20260315_H1437_paired_full_full_complete_matrix",
    "H2009":   "20260315_H2009_paired_full_full_complete_matrix",
    "HCC1937": "20260315_HCC1937_paired_full_full_complete_matrix",
    "HCC1954": "20260315_HCC1954_paired_full_full_complete_matrix",
}
CANCER = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
          "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}

# survey
METH_THR = 0.5
MIN_N = 3
MIN_PAIRED = 5
# regimeA filter
WP_MAX = 0.05
NCPG_STRICT = 100
NCPG_RELAX = 30
EXTREMITY_MAX = 0.3
# ARI (exact from script 30)
WINDOW = 600
ML_HIGH, ML_LOW = 200, 50
MIN_CORE_CPG = 4
MIN_SHARED = 3
MIN_GROUP = 8
ARI_PASS = 0.30
PLACEBO_MAX = 0.10
# pre-filter
CPG_PROX = 2000
SHORE = 2000


def paths_for(s):
    base = (f"/big7_disk/liaoyoyo2001/big7_disk_output/canonical/{s}/"
            f"paired_full/{RUNDIR[s]}/longphase_s")
    return f"{base}/{s}_tagged.bam", f"{base}/filtered_snv_tp.vcf.gz"


# ---------- annotation (stage 5, script 39) ----------
def load_cpg():
    isl = defaultdict(list)
    with gzip.open(CPG, "rt") as op:
        for line in op:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[1].startswith("chr"):
                try:
                    isl[f[1]].append((int(f[2]), int(f[3])))
                except ValueError:
                    continue
            elif len(f) >= 3 and f[0].startswith("chr"):
                try:
                    isl[f[0]].append((int(f[1]), int(f[2])))
                except ValueError:
                    continue
    for c in isl:
        isl[c].sort()
    return isl


def load_tss():
    tss = defaultdict(list)
    seen = set()
    with gzip.open(GTF, "rt") as op:
        for line in op:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            chrom, start, end, strand, attr = f[0], int(f[3]), int(f[4]), f[6], f[8]
            gene = None
            for kv in attr.split(";"):
                kv = kv.strip()
                if kv.startswith("gene_id"):
                    gene = kv.split('"')[1] if '"' in kv else kv.split()[-1]
                    break
            if gene is None:
                continue
            t = start if strand == "+" else end
            key = (chrom, t, gene)
            if key in seen:
                continue
            seen.add(key)
            tss[chrom].append((t, strand, gene))
    for c in tss:
        tss[c].sort()
    return tss


def nearest_tss(tss, chrom, pos):
    arr = tss.get(chrom, [])
    if not arr:
        return None, None, None
    best = None
    for t, strand, gene in arr:
        d = pos - t if strand == "+" else t - pos
        if best is None or abs(d) < abs(best[0]):
            best = (d, gene, strand)
    return best[1], best[0], best[2]


def cpg_context(isl, chrom, pos):
    arr = isl.get(chrom, [])
    mind = None
    for s, e in arr:
        if s <= pos <= e:
            return "island", 0
        d = min(abs(pos - s), abs(pos - e))
        if mind is None or d < mind:
            mind = d
    if mind is None:
        return "open", None
    return ("shore" if mind <= SHORE else "open"), mind


def near_cpg_island(isl, chrom, pos, prox=CPG_PROX):
    arr = isl.get(chrom, [])
    for s, e in arr:
        if s - prox <= pos <= e + prox:
            return True
    return False


# ---------- ARI core (stage 3, script 30 verbatim logic) ----------
def hp(r):
    for t, v in r.tags:
        if t == "HP":
            return str(v)
    return None


def mcalls_binary(r):
    o = {}
    try:
        mod = r.modified_bases
    except Exception:
        return o
    if not mod:
        return o
    r2 = {a: b for a, b in r.get_aligned_pairs(matches_only=False)
          if a is not None and b is not None}
    for k, calls in mod.items():
        if k[2] != 'm':
            continue
        for rp, ml in calls:
            rf = r2.get(rp)
            if rf is not None:
                o[rf] = 1 if ml >= ML_HIGH else (0 if ml <= ML_LOW else -1)
    return o


def mcalls_prob(r):
    """per ref-pos 5mC prob (ML/255) for survey beta."""
    o = {}
    try:
        mod = r.modified_bases
    except Exception:
        return o
    if not mod:
        return o
    r2 = {a: b for a, b in r.get_aligned_pairs(matches_only=False)
          if a is not None and b is not None}
    for k, calls in mod.items():
        if k[2] != 'm':
            continue
        for rp, ml in calls:
            rf = r2.get(rp)
            if rf is not None:
                o[rf] = ml / 255.0
    return o


def observed_only_distance(reads):
    n = len(reads)
    cov = Counter()
    for r in reads:
        for c in r:
            cov[c] += 1
    core = {c for c, k in cov.items() if k >= max(3, 0.3 * n)}
    reads = [{c: r[c] for c in r if c in core} for r in reads]
    keep = [i for i, r in enumerate(reads) if len(r) >= MIN_CORE_CPG]
    reads = [reads[i] for i in keep]
    n = len(reads)
    if n < 2 * MIN_GROUP:
        return None, keep
    while n > 2 * MIN_GROUP:
        bad = np.zeros(n)
        for i in range(n):
            inc = sum(1 for j in range(n) if i != j
                      and len(set(reads[i]) & set(reads[j])) < MIN_SHARED)
            bad[i] = inc / (n - 1)
        if bad.max() <= 0.4:
            break
        d = int(bad.argmax()); reads.pop(d); keep.pop(d); n -= 1
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            sh = set(reads[i]) & set(reads[j])
            D[i, j] = D[j, i] = (np.mean([reads[i][c] != reads[j][c] for c in sh])
                                 if len(sh) >= MIN_SHARED else np.nan)
    if np.isnan(D).any():
        D[np.isnan(D)] = np.nanmax(D) if not np.isnan(D).all() else 0.5
    return D, keep


def blind_ari(D, lab):
    preds = []
    for link in ('average', 'complete'):
        preds.append(AgglomerativeClustering(n_clusters=2, metric='precomputed',
                                             linkage=link).fit_predict(D))
    try:
        preds.append(SpectralClustering(n_clusters=2, affinity='precomputed',
                                        random_state=42).fit_predict(np.exp(-D / (D.std() + 1e-9))))
    except Exception:
        pass
    return (max(adjusted_rand_score(lab, p) for p in preds),
            max(adjusted_mutual_info_score(lab, p) for p in preds))


def eval_split(g0, g1):
    reads = g0 + g1
    lab0 = [0] * len(g0) + [1] * len(g1)
    D, keep = observed_only_distance(reads)
    if D is None:
        return None
    lab = np.array([lab0[i] for i in keep])
    if lab.sum() < MIN_GROUP or (len(lab) - lab.sum()) < MIN_GROUP:
        return None
    ari, ami = blind_ari(D, lab)
    try:
        sil = silhouette_score(D, lab, metric='precomputed')
    except Exception:
        sil = float('nan')
    return dict(n0=int((lab == 0).sum()), n1=int((lab == 1).sum()),
                ari=float(ari), ami=float(ami), sil=float(sil))


def placebo(g_germ, lens_germ):
    if len(g_germ) < 2 * MIN_GROUP:
        return None
    med = np.median(lens_germ)
    a = [g_germ[i] for i in range(len(g_germ)) if lens_germ[i] >= med]
    b = [g_germ[i] for i in range(len(g_germ)) if lens_germ[i] < med]
    if len(a) < MIN_GROUP or len(b) < MIN_GROUP:
        return None
    return eval_split(b, a)


def collect_ari(bam, chrom, pos, groups):
    var0 = pos - 1
    s, e = var0 - WINDOW, var0 + WINDOW
    g = {k: [] for k in groups}
    meta = {k: [] for k in groups}
    for r in bam.fetch(chrom, max(0, s), e):
        if r.flag & 0x900 or r.flag & 0x400:
            continue
        h = hp(r)
        if h not in groups:
            continue
        m = {p: st for p, st in mcalls_binary(r).items() if s <= p <= e and st >= 0}
        if len(m) >= MIN_CORE_CPG:
            g[h].append(m)
            meta[h].append(r.query_length or (r.reference_end - r.reference_start))
    return g, meta


# ---------- survey (stage 1-2, per-CpG paired beta + Wilcoxon) ----------
def survey_locus(bam, chrom, pos, win):
    """Return best HP-axis survey dict (n_paired_cpg, mean_delta, mean_germ_beta,
    wilcoxon_p, axis) over the two possible axes, or None."""
    var0 = pos - 1
    s, e = max(0, var0 - win), var0 + win
    # per CpG -> per group -> {read_id: prob}
    per = defaultdict(lambda: defaultdict(dict))
    for r in bam.fetch(chrom, s, e):
        if r.flag & 0x900 or r.flag & 0x400:
            continue
        h = hp(r)
        if h is None:
            continue
        rid = r.query_name
        for rf, prob in mcalls_prob(r).items():
            if s <= rf <= e:
                per[rf][h][rid] = prob

    def beta_of(cell):
        n = len(cell)
        if n == 0:
            return None, 0
        return sum(1 for v in cell.values() if v >= METH_THR) / n, n

    best = None
    for gg, gs, axis in (("1", "1-1", "HP1_vs_HP1-1"), ("2", "2-1", "HP2_vs_HP2-1")):
        pg, ps = [], []
        for rf, gd in per.items():
            bg, ng = beta_of(gd.get(gg, {}))
            bs, ns = beta_of(gd.get(gs, {}))
            if bg is not None and bs is not None and ng >= MIN_N and ns >= MIN_N:
                pg.append(bg); ps.append(bs)
        n = len(pg)
        if n < MIN_PAIRED:
            continue
        deltas = [b - a for a, b in zip(pg, ps)]
        try:
            wp = float(stats.wilcoxon(pg, ps).pvalue)
        except Exception:
            wp = float("nan")
        cand = dict(axis=axis, n_paired_cpg=n,
                    mean_germ_beta=round(float(np.mean(pg)), 4),
                    mean_delta=round(float(np.mean(deltas)), 4),
                    max_abs_delta=round(float(np.max(np.abs(deltas))), 4),
                    wilcoxon_p=wp)
        if best is None or cand["n_paired_cpg"] > best["n_paired_cpg"]:
            best = cand
    return best


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in RUNDIR:
        print(f"usage: {sys.argv[0]} <{'|'.join(RUNDIR)}> [SURVEY_WINDOW=1000]")
        sys.exit(2)
    sample = sys.argv[1]
    survey_win = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    os.makedirs(CS, exist_ok=True)
    bam_path, vcf_path = paths_for(sample)
    print(f"[62] sample={sample} cancer={CANCER[sample]} survey_window=+/-{survey_win}")

    print("[62] loading CpG islands + RefSeq TSS ...")
    isl = load_cpg()
    tss = load_tss()

    # load TP somatic SNVs + pre-filter to CpG-island-proximal
    vf = pysam.VariantFile(vcf_path)
    n_tp = 0
    cand_loci = []
    for rec in vf:
        if not (rec.alts and len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts)):
            continue
        n_tp += 1
        if near_cpg_island(isl, rec.chrom, rec.pos):
            cand_loci.append((rec.chrom, rec.pos))
    print(f"[62] TP somatic SNV total={n_tp}, CpG-island-proximal candidates={len(cand_loci)}")

    bam = pysam.AlignmentFile(bam_path, "rb")
    # funnel counters
    n_surveyed = 0
    n_hp_axis = 0
    n_regimeA_strict = 0
    n_regimeA_relax = 0
    n_ari_eval = 0
    credible_strict = []
    credible_relax = []

    for i, (chrom, pos) in enumerate(cand_loci):
        sv = survey_locus(bam, chrom, pos, survey_win)
        n_surveyed += 1
        if sv is None:
            continue
        n_hp_axis += 1
        extremity = abs(sv["mean_germ_beta"] - 0.5)
        passes_wp = (not np.isnan(sv["wilcoxon_p"])) and sv["wilcoxon_p"] < WP_MAX
        passes_ext = extremity <= EXTREMITY_MAX
        is_strict = passes_wp and passes_ext and sv["n_paired_cpg"] >= NCPG_STRICT
        is_relax = passes_wp and passes_ext and sv["n_paired_cpg"] >= NCPG_RELAX
        if not is_relax:
            continue
        n_regimeA_relax += 1
        if is_strict:
            n_regimeA_strict += 1
        # ARI gate
        gg, gs = ("1", "1-1") if sv["axis"] == "HP1_vs_HP1-1" else ("2", "2-1")
        g, meta = collect_ari(bam, chrom, pos, [gg, gs])
        if len(g[gg]) < MIN_GROUP or len(g[gs]) < MIN_GROUP:
            continue
        res = eval_split(g[gg], g[gs])
        if not res:
            continue
        n_ari_eval += 1
        pl = placebo(g[gg], meta[gg])
        placebo_ari = pl["ari"] if pl else float("nan")
        pass_tierA = (res["ari"] >= ARI_PASS) and (np.isnan(placebo_ari) or placebo_ari < PLACEBO_MAX)
        gene, tssd, strand = nearest_tss(tss, chrom, pos)
        ctx, cdist = cpg_context(isl, chrom, pos)
        rec = dict(
            locus=f"{chrom}:{pos}", axis=sv["axis"],
            n_paired_cpg=sv["n_paired_cpg"], mean_germ_beta=sv["mean_germ_beta"],
            delta=sv["mean_delta"], max_abs_delta=sv["max_abs_delta"],
            wilcoxon_p=round(sv["wilcoxon_p"], 6),
            ari=round(res["ari"], 4), ami=round(res["ami"], 4), sil=round(res["sil"], 4),
            n_germ=res["n0"], n_som=res["n1"],
            placebo_ari=(None if np.isnan(placebo_ari) else round(placebo_ari, 4)),
            pass_tierA=bool(pass_tierA),
            nearest_gene=gene, tss_dist=tssd, strand=strand,
            cpg_context=ctx, cpg_dist=cdist,
            tier=("strict_ge100" if is_strict else "relaxed_ge30"),
        )
        credible_relax.append(rec)
        if is_strict:
            credible_strict.append(rec)
        if (i + 1) % 200 == 0:
            print(f"     ...{i+1}/{len(cand_loci)} surveyed; regimeA_relax={n_regimeA_relax} "
                  f"ari_eval={n_ari_eval} credible(pass_tierA)={sum(r['pass_tierA'] for r in credible_relax)}")

    bam.close()
    n_pass_strict = sum(1 for r in credible_strict if r["pass_tierA"])
    n_pass_relax = sum(1 for r in credible_relax if r["pass_tierA"])
    credible_relax.sort(key=lambda r: (-int(r["pass_tierA"]), -abs(r["delta"]), -r["ari"]))

    out = dict(
        meta=dict(
            script="62_per_sample_credible_discovery.py", sample=sample,
            cancer_type=CANCER[sample], task_type="A pilot (per-sample credible discovery)",
            bam=bam_path, tp_vcf=vcf_path, survey_window=survey_win,
            method=("faithful HCC1395 replication: pysam survey (per-CpG paired HP-axis "
                    "beta+Wilcoxon) -> regimeA filter (wp<0.05, n_paired_cpg>=100|30, "
                    "extremity<=0.3) -> blind-ARI gate (Hamming+agglo/spectral k=2, "
                    "ARI>=0.30 & placebo<0.10) -> RefSeq/CpG-island annotation; "
                    "CpG-island pre-filter for tractability; cn_class SKIPPED "
                    "(SEQC2 CN is HCC1395-only)"),
            thresholds=dict(METH_THR=METH_THR, MIN_N=MIN_N, MIN_PAIRED=MIN_PAIRED,
                            wp_max=WP_MAX, ncpg_strict=NCPG_STRICT, ncpg_relax=NCPG_RELAX,
                            extremity_max=EXTREMITY_MAX, ari_pass=ARI_PASS,
                            placebo_max=PLACEBO_MAX, window=WINDOW, cpg_prox=CPG_PROX),
        ),
        funnel=dict(
            n_tp_somatic=n_tp, n_cpg_island_proximal=len(cand_loci),
            n_surveyed=n_surveyed, n_hp_axis_survey=n_hp_axis,
            n_regimeA_relax=n_regimeA_relax, n_regimeA_strict=n_regimeA_strict,
            n_ari_evaluable=n_ari_eval,
            n_credible_pass_tierA_strict=n_pass_strict,
            n_credible_pass_tierA_relax=n_pass_relax,
        ),
        credible_loci=credible_relax,
    )
    outp = f"{CS}/{sample}_credible_discovery.json"
    with open(outp, "w") as f:
        json.dump(out, f, indent=2,
                  default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
    print(f"[62] wrote {outp}")
    print(f"[62] FUNNEL {sample}: TP={n_tp} -> island-prox={len(cand_loci)} -> "
          f"HP-axis survey={n_hp_axis} -> regimeA(relax/strict)={n_regimeA_relax}/{n_regimeA_strict} "
          f"-> ARI-eval={n_ari_eval} -> CREDIBLE pass_tierA(relax/strict)={n_pass_relax}/{n_pass_strict}")


if __name__ == "__main__":
    main()
