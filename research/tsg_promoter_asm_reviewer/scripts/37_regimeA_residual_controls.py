#!/usr/bin/env python3
"""
37 — Regime-A residual controls on the 23 Tier-A pass loci: M4c CpG-context-removed
+ stronger M8 random-non-somatic-anchor placebo. Locks the Tier-A claim beyond the
same-locus length-split placebo already passed (real>>placebo p=5e-8).

Already passed for regime-A credible subset (30+36):
  M1 blind-ARI median 0.229; M2 vs het-null p=2.3e-4, Cliff δ=0.37, CI excl 0;
  M3 rarefied sil 0.308 vs null 0.080; M5 coverage ρ=-0.18 NS; M8(length-split) real>>placebo.

This script adds:
  M4c CpG-context-removed: drop CpGs within +/-20bp of the variant (ALT could rewrite
       local CpG motif) -> recompute blind-ARI; cluster must survive (|ΔARI| small, ARI>=0.30).
  M8-strong random-anchor: for each locus pick K random non-somatic, non-LOH anchors
       within +/-20kb; at each, split HP1(germline) reads by length median -> ARI.
       max anchor-ARI = genomic-context collider estimate; pass if < 0.10.

Input: regimeA_credible_probe.json (pass loci). Single-sample HCC1395. Tier-A.
Output: genome_survey_v2/regimeA_residual_controls.json
"""
import json, random
from pathlib import Path
import numpy as np
import pysam
from collections import defaultdict, Counter
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.metrics import adjusted_rand_score

random.seed(42); np.random.seed(42)
WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
PROBE = WS / "genome_survey_v2/regimeA_credible_probe.json"
LOH_BED = WS / "output/seqc2_loh_only.bed"
TPP = WS / "genome_survey_v2/asm_dualaxis_tp.tsv"
OUT = WS / "genome_survey_v2/regimeA_residual_controls.json"
WINDOW = 600; ML_HIGH, ML_LOW = 200, 50
MIN_CORE_CPG = 4; MIN_SHARED = 3; MIN_GROUP = 8
CTX_BP = 20; N_ANCHOR = 4; ANCHOR_SPAN = 20000
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


def collect(chrom, center0, groups):
    s, e = center0 - WINDOW, center0 + WINDOW
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


def blind_ari_from_groups(g0, g1):
    reads = g0 + g1; lab0 = [0] * len(g0) + [1] * len(g1)
    D, keep = observed_only_distance(reads)
    if D is None:
        return None
    lab = np.array([lab0[i] for i in keep])
    if lab.sum() < MIN_GROUP or (len(lab) - lab.sum()) < MIN_GROUP:
        return None
    preds = []
    for link in ('average', 'complete'):
        preds.append(AgglomerativeClustering(n_clusters=2, metric='precomputed', linkage=link).fit_predict(D))
    try:
        preds.append(SpectralClustering(n_clusters=2, affinity='precomputed', random_state=42).fit_predict(np.exp(-D / (D.std() + 1e-9))))
    except Exception:
        pass
    return max(adjusted_rand_score(lab, p) for p in preds)


def load_loh():
    regs = defaultdict(list)
    with open(LOH_BED) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                regs[p[0]].append((int(p[1]), int(p[2])))
    return regs


def in_loh(regs, chrom, pos):
    for s, e in regs.get(chrom, []):
        if s <= pos <= e:
            return True
    return False


def main():
    probe = json.load(open(PROBE))
    pass_loci = [r for r in probe['loci'] if r['pass_tierA']]
    loh = load_loh()
    # somatic positions (avoid anchoring near them)
    som_pos = defaultdict(set)
    import csv
    with open(TPP) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            som_pos[row['chrom']].add(int(row['somatic_pos']))

    print(f"residual controls on {len(pass_loci)} Tier-A pass loci")
    rows = []
    for r in pass_loci:
        chrom, pos = r['locus'].split(":"); pos = int(pos); var0 = pos - 1
        gg, gs = ("1", "1-1") if r['axis'] == "HP1_vs_HP1-1" else ("2", "2-1")
        g, lens = collect(chrom, var0, [gg, gs])
        if len(g[gg]) < MIN_GROUP or len(g[gs]) < MIN_GROUP:
            continue
        ari_full = blind_ari_from_groups(g[gg], g[gs])
        # M4c: drop CpGs within +/-CTX_BP of variant
        def drop_ctx(reads):
            return [{c: v for c, v in rd.items() if abs(c - var0) > CTX_BP} for rd in reads]
        g0c = [x for x in drop_ctx(g[gg]) if len(x) >= MIN_CORE_CPG]
        g1c = [x for x in drop_ctx(g[gs]) if len(x) >= MIN_CORE_CPG]
        ari_noctx = (blind_ari_from_groups(g0c, g1c) if len(g0c) >= MIN_GROUP and len(g1c) >= MIN_GROUP else None)
        # M8-strong: random non-somatic non-LOH anchors +/-ANCHOR_SPAN, split germline reads by length
        anchor_aris = []
        tries = 0
        while len(anchor_aris) < N_ANCHOR and tries < 40:
            tries += 1
            off = random.randint(-ANCHOR_SPAN, ANCHOR_SPAN)
            apos = pos + off
            if abs(off) < 2 * WINDOW:
                continue
            if in_loh(loh, chrom, apos):
                continue
            if any(abs(apos - sp) < 1000 for sp in som_pos.get(chrom, ())):
                continue
            ag, alens = collect(chrom, apos - 1, [gg])
            greads = ag[gg]; gl = alens[gg]
            if len(greads) < 2 * MIN_GROUP:
                continue
            med = np.median(gl)
            a = [greads[i] for i in range(len(greads)) if gl[i] >= med]
            b = [greads[i] for i in range(len(greads)) if gl[i] < med]
            if len(a) < MIN_GROUP or len(b) < MIN_GROUP:
                continue
            aa = blind_ari_from_groups(b, a)
            if aa is not None:
                anchor_aris.append(aa)
        max_anchor = max(anchor_aris) if anchor_aris else None
        m4c_survive = (ari_noctx is not None and ari_noctx >= 0.30 and
                       (ari_full is None or abs(ari_noctx - ari_full) / max(ari_full, 1e-9) < 0.5))
        m8s_pass = (max_anchor is None) or (max_anchor < 0.10)
        rows.append(dict(locus=r['locus'], axis=r['axis'], ari_full=ari_full,
                         ari_noctx=ari_noctx, m4c_survive=bool(m4c_survive),
                         anchor_aris=[round(a, 3) for a in anchor_aris],
                         max_anchor_ari=(None if max_anchor is None else round(max_anchor, 3)),
                         m8strong_pass=bool(m8s_pass),
                         all_pass=bool(m4c_survive and m8s_pass)))
        print(f"  {r['locus']:<18} ari_full={ari_full:.3f} ari_noCpGctx={ari_noctx if ari_noctx is None else round(ari_noctx,3)} "
              f"max_anchor={max_anchor if max_anchor is None else round(max_anchor,3)} "
              f"M4c={'OK' if m4c_survive else 'FAIL'} M8strong={'OK' if m8s_pass else 'COLLIDER'}")

    n = len(rows)
    n_m4c = sum(x['m4c_survive'] for x in rows)
    n_m8s = sum(x['m8strong_pass'] for x in rows)
    n_all = sum(x['all_pass'] for x in rows)
    anchors_all = [a for x in rows for a in x['anchor_aris']]
    summary = dict(
        n_pass_loci_tested=n,
        M4c_CpGcontext_survive=f"{n_m4c}/{n} ({100*n_m4c/n:.0f}%)" if n else "0",
        M8strong_random_anchor_pass=f"{n_m8s}/{n} ({100*n_m8s/n:.0f}%)" if n else "0",
        median_anchor_ari=(float(np.median(anchors_all)) if anchors_all else None),
        n_survive_BOTH=f"{n_all}/{n}" if n else "0",
        note="M4c: cluster survives dropping CpGs +/-20bp of variant. M8-strong: random non-somatic non-LOH anchor length-split ARI<0.10 = no genomic-context collider.",
        sample="HCC1395 single-sample Tier-A ceiling",
    )
    json.dump(dict(summary=summary, loci=rows), open(OUT, "w"), indent=1)
    print("\n=== RESIDUAL CONTROLS SUMMARY ===")
    print(json.dumps(summary, indent=1, ensure_ascii=False))
    print(f"saved -> {OUT}")


if __name__ == "__main__":
    main()
