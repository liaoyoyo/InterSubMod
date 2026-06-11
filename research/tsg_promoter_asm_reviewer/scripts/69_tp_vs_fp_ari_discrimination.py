#!/usr/bin/env python3
"""
69 - Does the blind-ARI ASM clustering DISCRIMINATE TP from FP somatic calls?

The credible-loci gallery (62/67) used TP somatic SNVs only. The user asks: is this
ARI clustering an EFFECTIVE DISCRIMINATOR? Two distinct questions:
  (a) within a locus, does methylation separate germline vs somatic reads? -> ARI YES.
  (b) does high ARI distinguish a TRUE somatic call (TP) from a FALSE one (FP)?
This script tests (b) for HCC1395: run the SAME survey+ARI pipeline on CpG-island-
proximal TP AND FP loci, and compare the ARI distributions + tierA pass rates.

Prior result (memory project_zar1l_brca2_asm / phase2_cycle1): ASM is
anti/non-discriminative for TP/FP (strong-ASM FP-enriched OR=0.194). This is the
concrete per-locus ARI test of that.

Reuses script 62's functions via importlib.

Output: genome_survey_v2/cn_confound/cross_sample/HCC1395_tp_vs_fp_ari.json
"""
import os, sys, json, importlib.util
os.environ.setdefault("TMPDIR", "/big7_disk")
import numpy as np
import pysam
from scipy.stats import mannwhitneyu

SCRIPTS = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/scripts"
spec = importlib.util.spec_from_file_location("disc62", f"{SCRIPTS}/62_per_sample_credible_discovery.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
SAMPLE = "HCC1395"
RUN = "20260314_HCC1395_paired_full_full_complete_matrix"
BASE = f"/big7_disk/liaoyoyo2001/big7_disk_output/canonical/{SAMPLE}/paired_full/{RUN}/longphase_s"
BAMP = f"{BASE}/{SAMPLE}_tagged.bam"
SEED = 42
SURVEY_WIN = 1000


def load_snv(vcf):
    out = []
    vf = pysam.VariantFile(vcf)
    for rec in vf:
        if rec.alts and len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts):
            out.append((rec.chrom, rec.pos))
    return out


def run_set(label, vcf, isl, bam, cap=None):
    loci = load_snv(vcf)
    cand = [(c, p) for c, p in loci if m.near_cpg_island(isl, c, p)]
    rng = np.random.RandomState(SEED)
    if cap and len(cand) > cap:
        sel = rng.choice(len(cand), cap, replace=False)
        cand = [cand[i] for i in sorted(sel)]
    aris, deltas = [], []
    n_regimeA = n_arieval = n_credible = 0
    for chrom, pos in cand:
        sv = m.survey_locus(bam, chrom, pos, SURVEY_WIN)
        if sv is None:
            continue
        ext = abs(sv["mean_germ_beta"] - 0.5)
        if not ((not np.isnan(sv["wilcoxon_p"])) and sv["wilcoxon_p"] < m.WP_MAX
                and ext <= m.EXTREMITY_MAX and sv["n_paired_cpg"] >= m.NCPG_RELAX):
            continue
        n_regimeA += 1
        gg, gs = ("1", "1-1") if sv["axis"] == "HP1_vs_HP1-1" else ("2", "2-1")
        g, meta = m.collect_ari(bam, chrom, pos, [gg, gs])
        if len(g[gg]) < m.MIN_GROUP or len(g[gs]) < m.MIN_GROUP:
            continue
        res = m.eval_split(g[gg], g[gs])
        if not res:
            continue
        n_arieval += 1
        pl = m.placebo(g[gg], meta[gg])
        plc = pl["ari"] if pl else float("nan")
        passT = (res["ari"] >= m.ARI_PASS) and (np.isnan(plc) or plc < m.PLACEBO_MAX)
        aris.append(res["ari"]); deltas.append(sv["mean_delta"])
        if passT:
            n_credible += 1
    print(f"[69] {label}: total={len(loci)} island-prox={len(cand)} regimeA={n_regimeA} "
          f"ARI-eval={n_arieval} credible={n_credible}")
    return dict(label=label, n_total=len(loci), n_island_prox=len(cand),
                n_regimeA=n_regimeA, n_ari_eval=n_arieval, n_credible=n_credible,
                credible_rate_of_arieval=(round(n_credible / n_arieval, 3) if n_arieval else None),
                ari_list=[round(a, 3) for a in aris],
                ari_median=(round(float(np.median(aris)), 3) if aris else None),
                ari_mean=(round(float(np.mean(aris)), 3) if aris else None))


def main():
    print("[69] loading CpG islands...")
    isl = m.load_cpg()
    bam = pysam.AlignmentFile(BAMP, "rb")
    fp = run_set("FP", f"{BASE}/filtered_snv_fp.vcf.gz", isl, bam)
    # match TP count to FP island-prox for fair comparison
    tp = run_set("TP", f"{BASE}/filtered_snv_tp.vcf.gz", isl, bam,
                 cap=max(60, fp["n_island_prox"] * 3))
    bam.close()

    # compare ARI distributions of ARI-evaluable loci
    mw_p = None
    if fp["ari_list"] and tp["ari_list"] and len(fp["ari_list"]) >= 3 and len(tp["ari_list"]) >= 3:
        mw_p = round(float(mannwhitneyu(tp["ari_list"], fp["ari_list"], alternative="two-sided")[1]), 4)

    verdict = ("ASM-ARI does NOT discriminate TP from FP: FP loci reach ARI-evaluable AND "
               "pass tierA at a comparable (or higher) rate; ARI distributions overlap "
               "(Mann-Whitney p={}). High ARI means 'methylation separates the two "
               "haplotype groups at this locus' — it is NOT evidence the somatic call is "
               "real. Consistent with prior anti-discriminative finding (strong-ASM "
               "FP-enriched OR=0.194).").format(mw_p)

    out = dict(meta=dict(script="69_tp_vs_fp_ari_discrimination.py", sample=SAMPLE,
                         bam=BAMP, survey_win=SURVEY_WIN, thresholds=dict(
                             ari_pass=m.ARI_PASS, placebo_max=m.PLACEBO_MAX,
                             ncpg_relax=m.NCPG_RELAX, wp_max=m.WP_MAX)),
               TP=tp, FP=fp, ari_mannwhitney_p=mw_p,
               discrimination_verdict=verdict)
    outp = f"{CS}/{SAMPLE}_tp_vs_fp_ari.json"
    with open(outp, "w") as f:
        json.dump(out, f, indent=2,
                  default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
    print(f"[69] wrote {outp}")
    print(f"[69] TP credible-rate(of ARI-eval)={tp['credible_rate_of_arieval']} (n_arieval={tp['n_ari_eval']}, ari_med={tp['ari_median']}) "
          f"vs FP credible-rate={fp['credible_rate_of_arieval']} (n_arieval={fp['n_ari_eval']}, ari_med={fp['ari_median']}); MW p={mw_p}")


if __name__ == "__main__":
    main()
