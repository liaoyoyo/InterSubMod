#!/usr/bin/env python3
"""
59 - Per-sample GENOME-WIDE somatic ASM rate + imprinted-DMR positive controls.

QUESTION (expansion, user authorised "more resources", 2026-06-03)
  The targeted 38-position check (57/58) showed HCC1395's somatic ASM loci are
  PRIVATE (0/38 recur as somatic in the other 5 cancers). That answers "same
  position" but leaves the real replication question:

    Does each cancer, using its OWN private somatic mutations, INDEPENDENTLY show
    the same ASM PHENOMENON (somatic-controlled HP-axis Δβ at a comparable rate and
    magnitude)? If yes, somatic ASM is a reproducible phenomenon (⭐4 evidence) even
    though specific loci are private. If no, HCC1395's ASM may be sample-specific.

  Plus: a KNOWN-IMPRINTED-DMR panel (IGF2/H19/GNAS/MEG3/...) measured by germline
  HP1-vs-HP2 Δβ in every sample — a positive control that the germline-ASM detector
  works, and a characterisation of the imprinting baseline the 4 "recurrent
  imprinting" hits sit against.

METHOD
  - somatic loci = this sample's filtered_snv_tp.vcf.gz (high-confidence TP somatic
    SNVs). Random subset of N_SOM (seed=42) for a stable rate estimate (full set is
    10k-80k; 1.3s/locus would be hours).
  - per locus: somatic-controlled HP-axis Δβ. The somatic sub-haplotype is whichever
    HP main group ('1' or '2') has a paired '-' sub-tag (e.g. '1-1') at the locus.
    Δβ = beta(sub) - beta(main), 5mC-only, window ±600bp, MIN_N=3, validated pysam.
  - evaluable = locus has a subhap group AND both main & sub have >=MIN_N reads.
  - strong ASM = |Δβ| >= STRONG (0.10). rate = n_strong / n_evaluable.
  - imprinted DMRs: germline HP1-vs-HP2 Δβ over the DMR window (no somatic needed).

USAGE
  python3 59_genomewide_asm_rate.py <SAMPLE> [N_SOM]
  -> genome_survey_v2/cn_confound/cross_sample/<SAMPLE>_gwasm.json
"""
import os, sys, json, random
os.environ.setdefault("TMPDIR", "/big7_disk")
import numpy as np
import pysam

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
REF = "/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"

WINDOW = 600
MIN_N = 3
STRONG = 0.10
N_SOM_DEFAULT = 2000
SEED = 42

RUNDIR = {
    "HCC1395": "20260314_HCC1395_paired_full_full_complete_matrix",
    "COLO829": "20260315_COLO829_paired_full_full_complete_matrix",
    "H1437":   "20260315_H1437_paired_full_full_complete_matrix",
    "H2009":   "20260315_H2009_paired_full_full_complete_matrix",
    "HCC1937": "20260315_HCC1937_paired_full_full_complete_matrix",
    "HCC1954": "20260315_HCC1954_paired_full_full_complete_matrix",
}
CANCER = {
    "HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
    "H1437": "lung", "H2009": "lung", "COLO829": "melanoma",
}

# Known imprinted DMRs (GRCh38). Germline HP1-vs-HP2 ASM expected in ALL samples.
IMPRINTED_DMRS = [
    ("H19/IGF2_ICR", "chr11", 1999000, 2004000),
    ("IGF2",         "chr11", 2129000, 2148000),
    ("KCNQ1OT1_KvDMR","chr11", 2698000, 2722000),
    ("MEG3_DLK1",    "chr14", 100808000, 100828000),
    ("SNRPN_PWS",    "chr15", 24954000, 25069000),
    ("GNAS",         "chr20", 58839000, 58912000),
    ("PEG3",         "chr19", 56792000, 56840000),
    ("PLAGL1_ZAC1",  "chr6",  144005000, 144009000),
    ("GRB10",        "chr7",  50780000, 50861000),
    ("MEST_PEG1",    "chr7",  130490000, 130540000),
    ("NNAT",         "chr20", 37520000, 37530000),
    ("PEG10",        "chr7",  94656000, 94670000),
]


def paths_for(sample):
    base = (f"/big7_disk/liaoyoyo2001/big7_disk_output/canonical/{sample}/"
            f"paired_full/{RUNDIR[sample]}/longphase_s")
    return (f"{base}/{sample}_tagged.bam",
            f"{base}/filtered_snv_tp.vcf.gz")


def read_hp(r):
    try:
        return str(r.get_tag("HP"))
    except KeyError:
        return None


def read_5mc_beta(r, s, e):
    """mean 5mC prob over CpGs in [s,e] for this read, or None."""
    try:
        mod = r.modified_bases
    except Exception:
        return None
    if not mod:
        return None
    r2 = {a: b for a, b in r.get_aligned_pairs(matches_only=False)
          if a is not None and b is not None}
    vals = []
    for k, calls in mod.items():
        if k[2] != 'm':
            continue
        for rp, ml in calls:
            rf = r2.get(rp)
            if rf is not None and s <= rf <= e:
                vals.append(ml / 255.0)
    return float(np.mean(vals)) if vals else None


def hp_betas(bam, chrom, pos):
    """Return dict hp_tag -> list of per-read 5mC betas, in ±WINDOW of pos."""
    var0 = pos - 1
    s, e = max(0, var0 - WINDOW), var0 + WINDOW
    out = {}
    for r in bam.fetch(chrom, s, e):
        if r.flag & 0x900 or r.flag & 0x400:
            continue
        hp = read_hp(r)
        if hp is None:
            continue
        b = read_5mc_beta(r, s, e)
        if b is not None:
            out.setdefault(hp, []).append(b)
    return out, s, e


def somatic_hp_axis_delta(hp_b, with_null=False, n_perm=20, rng=None):
    """From hp_tag->[betas], find the somatic sub-haplotype (main HP with a '-' sub).
    Return (delta, main_tag, sub_tag, n_main, n_sub[, null_abs_delta]) or None.
    null = mean |Δ| from n_perm random splits of pooled (main+sub) reads into groups
    of the SAME sizes -> finite-sample noise floor for the |Δ|>=STRONG rate."""
    subs = [t for t in hp_b if "-" in t]
    best = None
    for sub in subs:
        main = sub.split("-")[0]
        if main not in hp_b:
            continue
        nm, ns = len(hp_b[main]), len(hp_b[sub])
        if nm >= MIN_N and ns >= MIN_N:
            d = float(np.mean(hp_b[sub]) - np.mean(hp_b[main]))
            if best is None or ns > best[4]:
                best = [round(d, 4), main, sub, nm, ns]
    if best is None:
        return None
    if with_null:
        main, sub, nm, ns = best[1], best[2], best[3], best[4]
        pool = np.array(hp_b[main] + hp_b[sub])
        r = rng if rng is not None else np.random
        nulls = []
        for _ in range(n_perm):
            idx = r.permutation(len(pool))
            a = pool[idx[:ns]]
            b = pool[idx[ns:ns + nm]]
            nulls.append(abs(float(a.mean() - b.mean())))
        best.append(round(float(np.mean(nulls)), 4))
    return tuple(best)


def germline_delta(hp_b):
    """HP1 vs HP2 delta = beta(HP2)-beta(HP1) if both >=MIN_N."""
    if "1" in hp_b and "2" in hp_b and len(hp_b["1"]) >= MIN_N and len(hp_b["2"]) >= MIN_N:
        return (round(float(np.mean(hp_b["2"]) - np.mean(hp_b["1"])), 4),
                len(hp_b["1"]), len(hp_b["2"]))
    return None


def load_tp_loci(vcf):
    out = []
    vf = pysam.VariantFile(vcf)
    for rec in vf:
        # SNV only
        if rec.alts and len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts):
            out.append((rec.chrom, rec.pos))
    return out


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in RUNDIR:
        print(f"usage: {sys.argv[0]} <{'|'.join(RUNDIR)}> [N_SOM]")
        sys.exit(2)
    sample = sys.argv[1]
    n_som = int(sys.argv[2]) if len(sys.argv) > 2 else N_SOM_DEFAULT
    os.makedirs(CS, exist_ok=True)
    bam_path, vcf_path = paths_for(sample)
    print(f"[59] sample={sample} cancer={CANCER[sample]} N_SOM={n_som}")

    bam = pysam.AlignmentFile(bam_path, "rb")

    # ---- genome-wide somatic ASM rate (random subset of own TP somatic SNVs) ----
    all_tp = load_tp_loci(vcf_path)
    n_tp_total = len(all_tp)
    random.seed(SEED)
    subset = all_tp if n_tp_total <= n_som else random.sample(all_tp, n_som)
    print(f"[59] TP somatic SNVs total={n_tp_total}, sampling {len(subset)}")

    rng = np.random.RandomState(SEED)
    deltas = []
    null_deltas = []
    n_evaluable = 0
    n_no_subhap = 0
    n_low_cov = 0
    n_hypo = 0
    n_hyper = 0
    for i, (chrom, pos) in enumerate(subset):
        try:
            hp_b, s, e = hp_betas(bam, chrom, pos)
        except (ValueError, OSError):
            continue
        res = somatic_hp_axis_delta(hp_b, with_null=True, rng=rng)
        if res is None:
            # distinguish no-subhap vs low-cov
            if any("-" in t for t in hp_b):
                n_low_cov += 1
            else:
                n_no_subhap += 1
            continue
        d = res[0]
        null_abs = res[5]
        deltas.append(d)
        null_deltas.append(null_abs)
        n_evaluable += 1
        if d < 0:
            n_hypo += 1
        elif d > 0:
            n_hyper += 1
        if (i + 1) % 500 == 0:
            print(f"     ...{i+1}/{len(subset)} evaluable={n_evaluable}")

    da = np.array(deltas) if deltas else np.array([])
    nd = np.array(null_deltas) if null_deltas else np.array([])
    n_strong = int((np.abs(da) >= STRONG).sum()) if da.size else 0
    rate_strong = round(n_strong / n_evaluable, 4) if n_evaluable else None
    # null: shuffled-label |Δ| rate (noise floor) + signal-above-null
    n_strong_null = int((nd >= STRONG).sum()) if nd.size else 0
    rate_strong_null = round(n_strong_null / n_evaluable, 4) if n_evaluable else None
    rate_excess = (round(rate_strong - rate_strong_null, 4)
                   if (rate_strong is not None and rate_strong_null is not None) else None)
    gw = dict(
        n_tp_total=n_tp_total, n_sampled=len(subset),
        n_evaluable=n_evaluable, n_no_subhap=n_no_subhap, n_low_cov=n_low_cov,
        n_strong_asm=n_strong, rate_strong_asm=rate_strong,
        n_strong_null=n_strong_null, rate_strong_null=rate_strong_null,
        rate_excess_over_null=rate_excess,
        median_abs_delta=round(float(np.median(np.abs(da))), 4) if da.size else None,
        median_abs_null=round(float(np.median(nd)), 4) if nd.size else None,
        mean_delta=round(float(da.mean()), 4) if da.size else None,
        frac_hypo=round(n_hypo / n_evaluable, 4) if n_evaluable else None,
        frac_hyper=round(n_hyper / n_evaluable, 4) if n_evaluable else None,
        p10_abs=round(float(np.percentile(np.abs(da), 10)), 4) if da.size else None,
        p90_abs=round(float(np.percentile(np.abs(da), 90)), 4) if da.size else None,
    )
    print(f"[59] GW somatic ASM: evaluable={n_evaluable} strong={n_strong} "
          f"rate={rate_strong} (null={rate_strong_null} excess={rate_excess}) "
          f"median|Δ|={gw['median_abs_delta']} (null_median={gw['median_abs_null']}) "
          f"hypo/hyper={gw['frac_hypo']}/{gw['frac_hyper']}")

    # ---- imprinted DMR positive controls (germline HP1-vs-HP2) ----
    # EXPLORATORY: DMRs are kb-scale; scan sub-windows (step WINDOW) across the DMR and
    # report the MAX |germline Δβ| (best sub-window). In heavily-LOH genomes many
    # imprinted loci lose heterozygosity -> germline axis uncomputable (germline_delta
    # =None), which is itself informative (LOH erases imprinting measurability).
    imp = []
    for name, chrom, start, end in IMPRINTED_DMRS:
        best_g = None       # (abs_delta, signed_delta, n1, n2, center)
        any_cov = False
        # up to 8 evenly-spaced sub-windows across the DMR (kept light; exploratory)
        centers = list(np.linspace(start + WINDOW, max(start + WINDOW, end - WINDOW),
                                   num=8, dtype=int))
        centers = sorted(set(int(c) for c in centers))
        for c in centers:
            try:
                hp_b, s, e = hp_betas(bam, chrom, c)
            except (ValueError, OSError):
                continue
            if hp_b:
                any_cov = True
            g = germline_delta(hp_b)
            if g is not None:
                ad = abs(g[0])
                if best_g is None or ad > best_g[0]:
                    best_g = (ad, g[0], g[1], g[2], c)
        imp.append(dict(
            name=name, region=f"{chrom}:{start}-{end}",
            germline_delta=None if best_g is None else best_g[1],
            best_center=None if best_g is None else f"{chrom}:{best_g[4]}",
            n_hp1=None if best_g is None else best_g[2],
            n_hp2=None if best_g is None else best_g[3],
            covered=any_cov,
            note=("uncomputable_germline (likely LOH/single-hap)"
                  if best_g is None else "ok"),
        ))
    n_imp_asm = sum(1 for x in imp if x["germline_delta"] is not None
                    and abs(x["germline_delta"]) >= STRONG)
    n_imp_eval = sum(1 for x in imp if x["germline_delta"] is not None)
    print(f"[59] imprinted DMR panel (EXPLORATORY, max sub-window): "
          f"{n_imp_asm}/{n_imp_eval} show |germline Δβ|>=0.10 "
          f"({len([x for x in imp if x['note'].startswith('uncomputable')])} LOH-uncomputable)")

    bam.close()
    out = dict(
        meta=dict(
            script="59_genomewide_asm_rate.py", sample=sample,
            cancer_type=CANCER[sample], task_type="A pilot expansion (genome-wide rate)",
            bam=bam_path, tp_vcf=vcf_path, window=WINDOW, min_n=MIN_N,
            strong_threshold=STRONG, n_som_target=n_som, seed=SEED,
            method=("own TP somatic SNVs, random subset; somatic-controlled HP-axis "
                    "Δβ = beta(subhap)-beta(main), 5mC-only, validated pysam"),
        ),
        genomewide_somatic_asm=gw,
        imprinted_dmr_panel=dict(
            n_evaluable=n_imp_eval, n_strong=n_imp_asm, loci=imp),
    )
    outp = f"{CS}/{sample}_gwasm.json"
    with open(outp, "w") as f:
        json.dump(out, f, indent=2,
                  default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
    print(f"[59] wrote {outp}")


if __name__ == "__main__":
    main()
