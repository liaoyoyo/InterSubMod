#!/usr/bin/env python3
"""
Stage B - BAM single-pass: O2 read-level clustering ARI + Q5 error features.

Adapts the EXISTING template 24_cluster_eval_core.py (observed-only Hamming distance,
blind k=2 ARI, placebo-by-length, PERMANOVA, silhouette, PC1/NC1 ruler-validation,
MAPQ at meta). HP-axis ONLY (somatic-controlled clean axis).

Per locus: fetch BAM WINDOW=600 around somatic_pos, group reads by HP tag (gg vs gs),
compute clustering metrics + per-read error features (MAPQ, NM, supplementary frac).

Inputs : genome_survey_v2/cn_confound/master_o1_cn.tsv (Stage A; HP-axis rows used)
Outputs: genome_survey_v2/cn_confound/master_o2_error.tsv
         genome_survey_v2/cn_confound/ruler_validation.json
Return : compact summary (stdout JSON block)
"""
import os, sys, json, random
import numpy as np
import pandas as pd
import pysam
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score
from skbio.stats.distance import DistanceMatrix, permanova

random.seed(42); np.random.seed(42)

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
OUTDIR = f"{ROOT}/genome_survey_v2/cn_confound"
MASTER_O1 = f"{OUTDIR}/master_o1_cn.tsv"
MASTER_O2 = f"{OUTDIR}/master_o2_error.tsv"
RULER_JSON = f"{OUTDIR}/ruler_validation.json"
os.makedirs(OUTDIR, exist_ok=True)

BAM = ("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
       "20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam")
WINDOW = 600
ML_HIGH, ML_LOW = 200, 50
MIN_CORE_CPG = 4; MIN_SHARED = 3; MIN_GROUP = 8

bam = pysam.AlignmentFile(BAM, "rb")

# ===================== template-reused read helpers =====================
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

def collect(chrom, pos, groups):
    """Fetch primary reads in window, group by HP, return meth dicts + meta.
    meta tuple = (ref_start, ref_end, mapq, qlen, nm). Returns g, meta."""
    var0 = pos - 1
    s, e = var0 - WINDOW, var0 + WINDOW
    g = {k: [] for k in groups}; meta = {k: [] for k in groups}
    for r in bam.fetch(chrom, max(0, s), e):
        if r.flag & 0x900 or r.flag & 0x400:  # skip supp/secondary/dup here
            continue
        h = hp(r)
        if h not in groups:
            continue
        m = {p: st for p, st in mcalls(r).items() if s <= p <= e and st >= 0}
        if len(m) >= MIN_CORE_CPG:
            g[h].append(m)
            nm = r.get_tag('NM') if r.has_tag('NM') else None
            qlen = r.query_length or (r.reference_end - r.reference_start)
            meta[h].append((r.reference_start, r.reference_end,
                            r.mapping_quality, qlen, nm))
    return g, meta

def error_features(chrom, pos, groups):
    """Q5 per-site error aggregation over ALL reads (both HP groups) in window.
    Separate fetch counts supplementary/secondary (flag & 0x900) for frac_supp."""
    var0 = pos - 1
    s, e = var0 - WINDOW, var0 + WINDOW
    mapqs = []; nms = []; nm_per_kb = []
    n_total = 0; n_supp = 0
    for r in bam.fetch(chrom, max(0, s), e):
        if r.flag & 0x400:  # skip PCR/optical duplicates from all counts
            continue
        h = hp(r)
        if h not in groups:
            continue
        n_total += 1
        if r.flag & 0x900:  # supplementary or secondary
            n_supp += 1
            continue  # primary-only stats for mapq/nm below
        mapqs.append(r.mapping_quality)
        if r.has_tag('NM'):
            nm = r.get_tag('NM')
            nms.append(nm)
            qlen = r.query_length or (r.reference_end - r.reference_start)
            if qlen and qlen > 0:
                nm_per_kb.append(1000.0 * nm / qlen)
    median_mapq = float(np.median(mapqs)) if mapqs else np.nan
    frac_mapq_lt20 = float(np.mean(np.array(mapqs) < 20)) if mapqs else np.nan
    median_nm = float(np.median(nms)) if nms else np.nan
    mean_nm_per_kb = float(np.mean(nm_per_kb)) if nm_per_kb else np.nan
    frac_supp = float(n_supp / n_total) if n_total > 0 else np.nan
    return dict(median_mapq=median_mapq, frac_mapq_lt20=frac_mapq_lt20,
                median_nm=median_nm, mean_nm_per_kb=mean_nm_per_kb,
                frac_supp=frac_supp)

# ===================== template-reused distance / cluster =====================
def observed_only_distance(reads):
    """M0: Hamming over shared CpG, NO impute. Drop poorly-connected reads."""
    from collections import Counter
    n = len(reads)
    cov = Counter()
    for r in reads:
        for c in r:
            cov[c] += 1
    core = [c for c, k in cov.items() if k >= max(3, 0.3 * n)]
    reads = [{c: r[c] for c in r if c in core} for r in reads]
    keep = [i for i, r in enumerate(reads) if len(r) >= MIN_CORE_CPG]
    reads = [reads[i] for i in keep]
    n = len(reads)
    if n < 2 * MIN_GROUP:
        return None, None, keep
    while True:
        bad = np.zeros(n)
        for i in range(n):
            inc = 0
            for j in range(n):
                if i == j:
                    continue
                if len(set(reads[i]) & set(reads[j])) < MIN_SHARED:
                    inc += 1
            bad[i] = inc / (n - 1)
        if bad.max() <= 0.4 or n <= 2 * MIN_GROUP:
            break
        drop = int(bad.argmax())
        reads.pop(drop); keep.pop(drop); n -= 1
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            sh = set(reads[i]) & set(reads[j])
            if len(sh) < MIN_SHARED:
                D[i, j] = D[j, i] = np.nan
            else:
                d = np.mean([reads[i][c] != reads[j][c] for c in sh])
                D[i, j] = D[j, i] = d
    if np.isnan(D).any():
        D[np.isnan(D)] = np.nanmax(D)
    return D, reads, keep

def blind_ari(D, true_labels):
    """M1: blind k=2 cluster, max ARI vs true labels."""
    preds = []
    for link in ('average', 'complete'):
        preds.append(AgglomerativeClustering(
            n_clusters=2, metric='precomputed', linkage=link).fit_predict(D))
    try:
        aff = np.exp(-D / (D.std() + 1e-9))
        preds.append(SpectralClustering(
            n_clusters=2, affinity='precomputed',
            random_state=42).fit_predict(aff))
    except Exception:
        pass
    aris = [adjusted_rand_score(true_labels, p) for p in preds]
    amis = [adjusted_mutual_info_score(true_labels, p) for p in preds]
    return max(aris), max(amis)

def permR2(D, labels):
    dm = DistanceMatrix(D, [str(i) for i in range(len(labels))])
    res = permanova(dm, labels, permutations=499)
    F = res['test statistic']; R2 = F / (F + (len(labels) - 2))
    return R2, res['p-value']

def evaluate(reads_g0, reads_g1, name):
    reads = reads_g0 + reads_g1
    labels0 = [0] * len(reads_g0) + [1] * len(reads_g1)
    D, kreads, keep = observed_only_distance(reads)
    if D is None:
        return None
    labels = np.array([labels0[i] for i in keep])
    if labels.sum() < MIN_GROUP or (len(labels) - labels.sum()) < MIN_GROUP:
        return None
    ari, ami = blind_ari(D, labels)
    R2, p = permR2(D, labels.tolist())
    sil = silhouette_score(D, labels, metric='precomputed')
    return dict(name=name, n0=int((labels == 0).sum()), n1=int((labels == 1).sum()),
                ari=ari, ami=ami, R2=R2, p=p, sil=sil, n_drop=len(reads) - len(keep))

# ===================== META-VALIDATION (驗尺 — run once) =====================
def sim_reads(n, p_meth, ncpg=20, miss=0.4):
    out = []
    for _ in range(n):
        r = {}
        for c in range(ncpg):
            if random.random() > miss:
                r[c] = 1 if random.random() < p_meth else 0
        if len(r) >= MIN_CORE_CPG:
            out.append(r)
    return out

print("=" * 70); print("META-VALIDATION (先驗尺再用尺)"); print("=" * 70)
pc1 = evaluate(sim_reads(40, 0.65), sim_reads(40, 0.35), "PC1 simulated delta-beta=0.3")
one = sim_reads(80, 0.5)
nc1 = evaluate(one[:40], one[40:], "NC1 random-split (no real cluster)")
for r in [pc1, nc1]:
    if r:
        print(f"  {r['name']:<34} ARI={r['ari']:.3f} AMI={r['ami']:.3f} "
              f"R2={r['R2']:.3f} p={r['p']:.3f} sil={r['sil']:.3f}")
ruler_pc1_ari = float(pc1['ari']) if pc1 else None
ruler_nc1_ari = float(nc1['ari']) if nc1 else None
ruler_ok = bool(pc1 and nc1 and pc1['ari'] > 0.5 and nc1['ari'] < 0.15)
print(f"  >>> 尺驗證: PC1 ARI>0.5={pc1['ari']>0.5 if pc1 else '?'} AND "
      f"NC1 ARI<0.15={nc1['ari']<0.15 if nc1 else '?'} -> "
      f"{'VALIDATED' if ruler_ok else 'RULER BROKEN'}")
with open(RULER_JSON, "w") as fh:
    json.dump(dict(pc1_ari=ruler_pc1_ari, nc1_ari=ruler_nc1_ari,
                   pc1_pass=(ruler_pc1_ari is not None and ruler_pc1_ari > 0.5),
                   nc1_pass=(ruler_nc1_ari is not None and ruler_nc1_ari < 0.15),
                   ruler_ok=ruler_ok,
                   thresholds=dict(pc1_min=0.5, nc1_max=0.15)), fh, indent=2)

# ===================== LOCUS SELECTION =====================
AXIS_GG_GS = {
    "HP1_vs_HP1-1": ("1", "1-1"),
    "HP2_vs_HP2-1": ("2", "2-1"),
}
df = pd.read_csv(MASTER_O1, sep="\t")
hpdf = df[df.axis_type == "HP"].copy()
hpdf = hpdf[hpdf.axis.isin(AXIS_GG_GS)]

# (i) sig-ASM: wilcoxon_p<0.05 ; if >800 keep top 800 by abs_delta
sig = hpdf[hpdf.wilcoxon_p < 0.05].copy()
if len(sig) > 800:
    sig = sig.sort_values("abs_delta", ascending=False).head(800)
sel_sig = sig[["chrom", "somatic_pos", "axis"]].copy()
sel_sig["sel_src"] = "sig"

# (ii) CN-stratified random sample (seed=42): up to 150 each from {loss,neutral,gain}
strat_parts = []
for cn in ["loss", "neutral", "gain"]:
    pool = hpdf[hpdf.cn_class == cn]
    take = min(150, len(pool))
    if take > 0:
        s = pool.sample(n=take, random_state=42)[["chrom", "somatic_pos", "axis"]].copy()
        s["sel_src"] = f"strat_{cn}"
        strat_parts.append(s)
sel_strat = pd.concat(strat_parts, ignore_index=True) if strat_parts else pd.DataFrame()

# (iii) canonical trees (always include)
CANON = [("chr13", 32315128, "HP1_vs_HP1-1"),
         ("chr20", 50267392, "HP1_vs_HP1-1"),
         ("chr12", 31601630, "HP2_vs_HP2-1"),
         ("chr9",  42376881, "HP2_vs_HP2-1"),
         ("chr16", 17774746, "HP2_vs_HP2-1")]
sel_canon = pd.DataFrame(CANON, columns=["chrom", "somatic_pos", "axis"])
sel_canon["sel_src"] = "canon"

sel = pd.concat([sel_sig, sel_strat, sel_canon], ignore_index=True)
# dedupe by (chrom,somatic_pos,axis), prefer canon > sig > strat priority order
src_rank = {"canon": 0, "sig": 1, "strat_loss": 2, "strat_neutral": 2, "strat_gain": 2}
sel["rank"] = sel["sel_src"].map(src_rank)
sel = sel.sort_values("rank").drop_duplicates(
    subset=["chrom", "somatic_pos", "axis"], keep="first").drop(columns="rank")

# attach Stage A columns
meta_cols = ["chrom", "somatic_pos", "axis", "is_tp", "cn_class", "median_cn",
             "loh_status", "n_paired_cpg", "abs_delta", "mean_delta"]
sel = sel.merge(hpdf[meta_cols].drop_duplicates(subset=["chrom", "somatic_pos", "axis"]),
                on=["chrom", "somatic_pos", "axis"], how="left")
n_loci = len(sel)
print(f"\nLocus selection: sig={len(sel_sig)} strat={len(sel_strat)} "
      f"canon={len(sel_canon)} -> deduped total = {n_loci}")
print("  cn_class breakdown:\n", sel.cn_class.value_counts().to_string())

# ===================== APPLY per-locus =====================
print("\n" + "=" * 70)
print("APPLY: per-locus blind-ARI + placebo + error features (HP-axis)")
print("=" * 70)
rows = []
n_processed = 0; n_skipped = 0
for idx, row in sel.reset_index(drop=True).iterrows():
    if idx % 200 == 0:
        print(f"  [{idx}/{n_loci}] processed={n_processed} skipped={n_skipped} ...",
              flush=True)
    chrom = row["chrom"]; pos = int(row["somatic_pos"]); axis = row["axis"]
    gg, gs = AXIS_GG_GS[axis]
    try:
        g, meta = collect(chrom, pos, [gg, gs])
        if len(g[gg]) < MIN_GROUP or len(g[gs]) < MIN_GROUP:
            n_skipped += 1
            continue
        r = evaluate(g[gg], g[gs], f"{chrom}:{pos}")
        if not r:
            n_skipped += 1
            continue
        # M8 placebo: split germline-gg reads by read-length median
        gm = meta[gg]
        lengths = np.array([m[3] for m in gm])
        med = np.median(lengths)
        pseudo = [i for i in range(len(g[gg])) if lengths[i] >= med]
        rest = [i for i in range(len(g[gg])) if lengths[i] < med]
        placebo_ari = np.nan
        if len(pseudo) >= MIN_GROUP and len(rest) >= MIN_GROUP:
            pr = evaluate([g[gg][i] for i in rest],
                          [g[gg][i] for i in pseudo], "placebo")
            if pr:
                placebo_ari = pr["ari"]
        # Q5 error features over the same window (both groups)
        ef = error_features(chrom, pos, [gg, gs])
        rows.append(dict(
            chrom=chrom, somatic_pos=pos, axis=axis,
            is_tp=int(row["is_tp"]) if pd.notna(row["is_tp"]) else "",
            cn_class=row["cn_class"], median_cn=row["median_cn"],
            loh_status=row["loh_status"], n_paired_cpg=row["n_paired_cpg"],
            abs_delta=row["abs_delta"], mean_delta=row["mean_delta"],
            blind_ari=round(float(r["ari"]), 5),
            placebo_ari=(round(float(placebo_ari), 5)
                         if not np.isnan(placebo_ari) else ""),
            permanova_r2=round(float(r["R2"]), 5),
            silhouette=round(float(r["sil"]), 5),
            n0=int(r["n0"]), n1=int(r["n1"]),
            median_mapq=ef["median_mapq"],
            frac_mapq_lt20=(round(ef["frac_mapq_lt20"], 5)
                            if not np.isnan(ef["frac_mapq_lt20"]) else ""),
            median_nm=(ef["median_nm"] if not np.isnan(ef["median_nm"]) else ""),
            mean_nm_per_kb=(round(ef["mean_nm_per_kb"], 4)
                            if not np.isnan(ef["mean_nm_per_kb"]) else ""),
            frac_supp=(round(ef["frac_supp"], 5)
                       if not np.isnan(ef["frac_supp"]) else ""),
        ))
        n_processed += 1
    except Exception as ex:
        n_skipped += 1
        if idx < 20 or idx % 500 == 0:
            print(f"    WARN locus {chrom}:{pos} {axis} failed: {ex}", flush=True)
        continue

print(f"  [DONE] processed={n_processed} skipped={n_skipped}", flush=True)

out_cols = ["chrom", "somatic_pos", "axis", "is_tp", "cn_class", "median_cn",
            "loh_status", "n_paired_cpg", "abs_delta", "mean_delta",
            "blind_ari", "placebo_ari", "permanova_r2", "silhouette",
            "n0", "n1", "median_mapq", "frac_mapq_lt20", "median_nm",
            "mean_nm_per_kb", "frac_supp"]
odf = pd.DataFrame(rows, columns=out_cols)
odf.to_csv(MASTER_O2, sep="\t", index=False)

# ===================== compact summary =====================
ba = odf["blind_ari"].astype(float)
mq = pd.to_numeric(odf["median_mapq"], errors="coerce")
summary = dict(
    n_loci_processed=int(n_processed),
    n_skipped=int(n_skipped),
    ruler_pc1_ari=ruler_pc1_ari,
    ruler_nc1_ari=ruler_nc1_ari,
    ruler_ok=ruler_ok,
    blind_ari_median=(round(float(ba.median()), 4) if len(ba) else None),
    blind_ari_q25=(round(float(ba.quantile(0.25)), 4) if len(ba) else None),
    blind_ari_q75=(round(float(ba.quantile(0.75)), 4) if len(ba) else None),
    median_mapq_overall=(round(float(mq.median()), 2) if mq.notna().any() else None),
    n_selected=int(n_loci),
    output_master_o2=MASTER_O2,
    output_ruler_json=RULER_JSON,
)
print("\n" + "=" * 70)
print("COMPACT_SUMMARY_JSON")
print(json.dumps(summary, indent=2))
