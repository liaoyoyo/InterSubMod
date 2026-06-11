#!/usr/bin/env python3
"""
52_subclone_assessment.py  (S - subclone assessment synthesis)

Synthesizes:
  - M2 (measurement audit): per-read methylation bimodality + binarized-vs-continuous
    divergence + 5mC/5hmC separation
  - C2 (clustering context): blind_ari (real read clustering) vs placebo (length-matched)
  - master_o2_error.tsv: per-locus blind_ari / placebo_ari / permanova / silhouette

Answers (HCC1395, single sample, A pilot):
  (1) Do read-level meth clustering (blind_ari) + per-read bimodality jointly indicate a
      SUBCLONE (minor methylation-distinct read subpopulation), or only an HP/germ-vs-som split?
  (2) Is binarized MAX-collapse HIDING subclone signal that continuous / separated 5mC|5hmC
      would reveal? (cross-ref M2 binarized-vs-continuous divergence + 5hmC orthogonality)
  (3) Top subclone-candidate loci (high blind_ari + bimodal continuous beta + TP); BRCA2?

KEY METHODOLOGICAL POINT (subclone vs HP-split discriminator)
  The M2 per-read bimodal GMM is fit on ALL reads at a locus (germline-group + somatic-subgroup
  POOLED). A bimodal per-read beta is therefore NOT automatically a subclone -- it can be the
  trivial germ-vs-som HP split. We discriminate:
    * HP-split-explained : GMM component weights ~ match the germ/som read fraction
      (the two modes ARE the two HP groups) -> NOT a novel subclone.
    * subclone-candidate : modes do NOT align with germ/som fraction (a minor read
      subpopulation exists WITHIN/ACROSS HP groups) -> genuine sub-structure.
  We cross-corroborate with blind_ari >> placebo_ari (C2/H-E: real read clustering, not a
  length/coverage artifact, p=1.8e-43 globally).

NOTE: HP-axis = somatic-controlled axis (germ HP-group vs its somatic HP-subgroup). FP power is
  very low (43 FP in 596 BAM-pass loci) -> per-cell counts reported, no classifier claimed.
"""
import json
import os
import numpy as np
import pandas as pd
import pysam
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CN = f"{BASE}/genome_survey_v2/cn_confound"
FIGDIR = f"{BASE}/figures/cn_confound"
os.makedirs(CN, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

M2_PATH = f"{CN}/m2_meth_measurement_audit.json"
C2_PATH = f"{CN}/c2_context_tpfp_clustering.json"
O2_PATH = f"{CN}/master_o2_error.tsv"
BAM = ("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
       "20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam")
OUT_JSON = f"{CN}/s_subclone_assessment.json"
OUT_FIG = f"{FIGDIR}/s_subclone_assessment.png"

CJK = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
fp = FontProperties(fname=CJK) if os.path.exists(CJK) else None


def main():
    m2 = json.load(open(M2_PATH))
    c2 = json.load(open(C2_PATH))
    o2 = pd.read_csv(O2_PATH, sep="\t")
    o2["key"] = o2["chrom"] + ":" + o2["somatic_pos"].astype(str) + ":" + o2["axis"]

    msum = m2["summary"]

    # ---- per-locus M2 -> DataFrame (skip insufficient-read loci) -------------
    rows = []
    for p in m2["per_locus"]:
        if "bimodal" not in p:  # skipped (insufficient germ reads)
            continue
        bm = p["bimodal"]
        # HP-split vs subclone discriminator:
        # expected HP-split weight = som / (germ + som); the somatic subgroup is the
        # smaller / displaced mode. We test whether either GMM weight ~ matches the
        # germ-or-som fraction (mode = HP group) vs deviates (extra structure).
        ng, ns = p.get("n_germ_reads", np.nan), p.get("n_som_reads", np.nan)
        frac_som = ns / (ng + ns) if (ng + ns) else np.nan
        weights = bm.get("comp_weights") or [np.nan, np.nan]
        wmin = min(weights) if len(weights) == 2 else np.nan
        # distance of the minor mode weight from BOTH the germ-frac and som-frac.
        # if minor-mode weight matches frac_som (or frac_germ) -> HP split explains it.
        d_to_som = abs(wmin - frac_som) if frac_som == frac_som else np.nan
        d_to_germ = abs(wmin - (1 - frac_som)) if frac_som == frac_som else np.nan
        hp_split_alignment = min(d_to_som, d_to_germ) if d_to_som == d_to_som else np.nan
        rows.append({
            "chrom": p["chrom"], "pos": p["pos"], "axis": p.get("axis"),
            "key": f"{p['chrom']}:{p['pos']}:{p.get('axis')}",
            "is_tp": p["is_tp"], "loh": p.get("loh_status"), "cn_class": p.get("cn_class"),
            "median_cn": p.get("median_cn"), "label": p.get("label"),
            "n_paired_cpg": p.get("n_paired_cpg"),
            "n_germ": ng, "n_som": ns, "frac_som": frac_som,
            "abs_delta": p.get("abs_delta"),
            "beta_bin": p.get("beta_binarized"), "beta_cont": p.get("beta_continuous"),
            "frac_intermediate": p.get("frac_intermediate_reads"),
            "db_5mC": p.get("db_5mC"), "db_5hmC": p.get("db_5hmC"), "db_max": p.get("db_max"),
            "per_read_beta_mean": p.get("per_read_beta_mean"),
            "per_read_beta_std": p.get("per_read_beta_std"),
            "binconti_div": p.get("_binconti_div"),
            "is_bimodal": bm.get("is_bimodal"),
            "bm_dbic": bm.get("delta_bic"), "bm_sep": bm.get("mean_sep"),
            "bm_minw": bm.get("min_weight"),
            "bm_means": bm.get("comp_means"), "bm_weights": weights,
            "minor_mode_w": wmin,
            "hp_split_alignment": hp_split_alignment,  # small => HP-split explains the modes
        })
    df = pd.DataFrame(rows)

    # join blind_ari / placebo_ari from o2 (BAM-pass clustering subset)
    df = df.merge(
        o2[["key", "blind_ari", "placebo_ari", "permanova_r2", "silhouette", "n0", "n1"]],
        on="key", how="left")

    # STRUCTURAL CAVEAT: M2 (per-read bimodality) and C2 (blind_ari clustering) sampled
    # largely DISJOINT locus sets. Only ~15 loci overlap; of the 43 bimodal loci only 4 have
    # blind_ari. We therefore CANNOT gate candidates on (bimodal AND high blind_ari) -- that
    # intersection is near-empty by sampling design, not by biology. We instead rank candidates
    # by the DIRECT per-read subpopulation test (M2 bimodality, subclone-like) and ANNOTATE
    # blind_ari where present, and report the overlap separately.
    n_m2_testable = int(df.shape[0])
    n_with_ari = int(df["blind_ari"].notna().sum())
    n_bimodal_with_ari = int(((df["is_bimodal"] == True) & df["blind_ari"].notna()).sum())

    # ---------------------------------------------------------------------
    # Q1: subclone vs HP-split.
    #  - bimodal AND blind_ari high AND minor mode NOT aligned to germ/som frac
    #    => structure beyond the 2 HP groups => subclone candidate.
    #  - HP_SPLIT_TOL: minor-mode weight within this of germ-or-som frac => HP-split-explained.
    # ---------------------------------------------------------------------
    HP_SPLIT_TOL = 0.12
    ARI_HI = c2["thresholds"]["blind_ari_median"]  # 0.2584
    bim = df[df["is_bimodal"] == True].copy()
    bim_with_ari = bim[bim["blind_ari"].notna()].copy()

    bim["hp_split_explained"] = bim["hp_split_alignment"] <= HP_SPLIT_TOL
    n_bim = len(bim)
    n_bim_hp_explained = int(bim["hp_split_explained"].sum())
    n_bim_subclone_like = int((~bim["hp_split_explained"]).sum())

    # blind_ari vs placebo separation among bimodal loci that have clustering data
    bim_ari = bim[bim["blind_ari"].notna()]
    blind_med = float(bim_ari["blind_ari"].median()) if len(bim_ari) else None
    plac_med = (float(bim_ari["placebo_ari"].median())
                if bim_ari["placebo_ari"].notna().any() else None)

    # ---------------------------------------------------------------------
    # Q2: does binarized MAX-collapse HIDE subclone signal?
    #  evidence axes:
    #   (a) binarized vs continuous per-locus beta corr (global, M2)  -> 0.9979 (tight)
    #   (b) BUT top_divergent loci show div up to 0.13 (continuous shifts beta) when
    #       intermediate-methylation reads (frac_intermediate) are present -> these are exactly
    #       the loci where a minor partially-methylated subpopulation gets crushed to 0/1.
    #   (c) 5hmC orthogonality: corr(db_max, db_5hmC)=0.033 => MAX-collapse discards 5hmC axis;
    #       loci with high loc_5hmC_gt_5mC_frac carry an orthogonal mod-channel that binarized
    #       max hides entirely.
    # ---------------------------------------------------------------------
    bin_cont_corr = msum["binarized_vs_continuous_corr"]
    corr_max_5hmc = msum["corr_db_max_vs_5hmC"]
    frac_5hmc_gt_5mc = msum["frac_5hmC_gt_5mC"]
    deadzone_frac = msum["deadzone_frac"]
    # loci where continuous would reveal hidden intermediate structure that binarized crushes:
    #   bimodal in per-read beta AND high frac_intermediate (reads stuck mid-methylation).
    df["intermediate_hi"] = df["frac_intermediate"] > 0.15
    hidden = df[(df["is_bimodal"] == True) & (df["intermediate_hi"])].copy()
    n_hidden = len(hidden)

    # cross-tab the top_divergent loci with bimodality
    topdiv = pd.DataFrame(msum["top_divergent_loci"])
    topdiv["key2"] = topdiv["chrom"] + ":" + topdiv["pos"].astype(str)
    df["key2"] = df["chrom"] + ":" + df["pos"].astype(str)
    topdiv_bim = topdiv.merge(df[["key2", "is_bimodal", "blind_ari", "is_tp"]],
                              on="key2", how="left")

    measurement_hides = bool((corr_max_5hmc < 0.1 and frac_5hmc_gt_5mc > 0.05)
                             or (n_hidden >= 3))

    # ---------------------------------------------------------------------
    # Q3: top subclone-candidate loci.
    #  Direct subpopulation test = M2 per-read bimodality. Gate = TP & bimodal & subclone-like
    #  (minor mode NOT explained by germ/som HP split). blind_ari is ANNOTATED where present
    #  (only ~4 bimodal loci have it -- M2/C2 sampled disjoint sets) and BOOSTS the score when
    #  it confirms real clustering, but is NOT a hard requirement (would force an empty set).
    #  Score prioritizes a real, sizeable, well-separated minor subpopulation.
    # ---------------------------------------------------------------------
    cand = df[(df["is_bimodal"] == True) & (df["is_tp"] == 1)].copy()
    cand["hp_split_explained"] = cand["hp_split_alignment"] <= HP_SPLIT_TOL
    # blind_ari confirmation bonus: 1.5x if blind_ari>=median (real clustering corroborated),
    # 1.0x if no clustering data, 0.6x if blind_ari present but low (clustering refutes).
    def _ari_factor(a):
        if pd.isna(a):
            return 1.0
        return 1.5 if a >= ARI_HI else 0.6
    cand["ari_factor"] = cand["blind_ari"].apply(_ari_factor)
    # base subpopulation strength = separation * minor-mode mass; downweight HP-split-explained.
    cand["subclone_score"] = (cand["bm_sep"].fillna(0)
                              * cand["minor_mode_w"].fillna(0)
                              * cand["ari_factor"]
                              * np.where(cand["hp_split_explained"], 0.4, 1.0))
    cand = cand.sort_values("subclone_score", ascending=False)
    # subclone CANDIDATES (strict) = TP & bimodal & subclone-like (not HP-split-explained)
    cand["is_candidate"] = ~cand["hp_split_explained"]

    top_candidates = []
    for _, r in cand[cand["is_candidate"]].head(12).iterrows():
        top_candidates.append({
            "locus": f"{r['chrom']}:{int(r['pos'])}", "axis": r["axis"],
            "label": r["label"], "loh": r["loh"], "cn_class": r["cn_class"],
            "median_cn": r["median_cn"], "is_tp": int(r["is_tp"]),
            "blind_ari": (round(float(r["blind_ari"]), 4)
                          if pd.notna(r["blind_ari"]) else None),
            "blind_ari_available": bool(pd.notna(r["blind_ari"])),
            "placebo_ari": (round(float(r["placebo_ari"]), 4)
                            if pd.notna(r["placebo_ari"]) else None),
            "beta_cont": round(float(r["beta_cont"]), 4),
            "bm_dbic": round(float(r["bm_dbic"]), 1) if pd.notna(r["bm_dbic"]) else None,
            "bm_sep": round(float(r["bm_sep"]), 4) if pd.notna(r["bm_sep"]) else None,
            "bm_means": r["bm_means"], "bm_weights": r["bm_weights"],
            "minor_mode_w": round(float(r["minor_mode_w"]), 4)
            if pd.notna(r["minor_mode_w"]) else None,
            "frac_som": round(float(r["frac_som"]), 3) if pd.notna(r["frac_som"]) else None,
            "n_germ": int(r["n_germ"]) if pd.notna(r["n_germ"]) else None,
            "n_som": int(r["n_som"]) if pd.notna(r["n_som"]) else None,
            "hp_split_explained": bool(r["hp_split_explained"]),
            "frac_intermediate": round(float(r["frac_intermediate"]), 4)
            if pd.notna(r["frac_intermediate"]) else None,
            "subclone_score": round(float(r["subclone_score"]), 5),
        })

    # ---- BRCA2 verdict ------------------------------------------------------
    brca2_m2 = next((p for p in m2["per_locus"]
                     if p["chrom"] == "chr13" and p["pos"] == 32315128), None)
    brca2_c2 = c2.get("brca2_canonical", {})
    brca2_bim = brca2_m2.get("bimodal", {}) if brca2_m2 else {}
    brca2 = {
        "locus": "chr13:32315128", "label": "BRCA2",
        "is_tp": brca2_c2.get("is_tp"), "loh": brca2_c2.get("loh"),
        "cn_class": brca2_c2.get("cn_class"), "median_cn": brca2_c2.get("median_cn"),
        "blind_ari": brca2_c2.get("blind_ari"),
        "abs_delta": brca2_c2.get("abs_delta"),
        "beta_bin": brca2_m2.get("beta_binarized") if brca2_m2 else None,
        "beta_cont": brca2_m2.get("beta_continuous") if brca2_m2 else None,
        "binconti_div": brca2_m2.get("_binconti_div") if brca2_m2 else None,
        "frac_intermediate": brca2_m2.get("frac_intermediate_reads") if brca2_m2 else None,
        "is_bimodal": brca2_bim.get("is_bimodal"),
        "bm_delta_bic": brca2_bim.get("delta_bic"),
        "bm_means": brca2_bim.get("comp_means"),
        "bm_weights": brca2_bim.get("comp_weights"),
        "bm_sep": brca2_bim.get("mean_sep"),
    }
    # BRCA2 qualifies as subclone-candidate?
    brca2_qualifies = bool(brca2["is_tp"] == 1
                           and brca2["blind_ari"] is not None
                           and brca2["blind_ari"] >= ARI_HI
                           and brca2["is_bimodal"] is True)
    brca2["qualifies_subclone_candidate"] = brca2_qualifies
    brca2["why"] = (
        "BRCA2 is TP, nonLOH, CN-gain (4n), high blind_ari=0.790 (real read clustering: "
        "blind 0.79 vs placebo 0.13, permanova R2=0.70, silhouette 0.61) AND continuous beta "
        "(0.162) exceeds binarized beta (0.072) -> a partially-methylated/intermediate read "
        "subpopulation (frac_intermediate=0.147) is being crushed by MAX-collapse binarization. "
        "BUT per-read GMM is NOT bimodal (dBIC=-6.6 favours 1 component): the intermediate reads "
        "form a graded shoulder, not a cleanly separated minor mode. So BRCA2 shows the "
        "MEASUREMENT-HIDING signature (continuous reveals what binarized hides) but does NOT meet "
        "the strict per-read-bimodality bar for a discrete subclone."
    ) if not brca2_qualifies else (
        "BRCA2 meets bimodal + high blind_ari + TP subclone bar."
    )

    # =====================================================================
    # FIGURE: (A) blind_ari vs per-read bimodality strength, colored by is_tp,
    #             marker = HP-split-explained vs subclone-like
    #         (B) per-read continuous-beta distribution (GMM modes) for top candidates,
    #             fetched live from BAM (3-5 loci)
    # =====================================================================
    fig = plt.figure(figsize=(15, 6.2))
    gsA = fig.add_axes([0.06, 0.12, 0.40, 0.78])
    # bimodality strength proxy = delta_bic (only meaningful where testable)
    plot_df = df[df["blind_ari"].notna() & df["bm_dbic"].notna()].copy()
    plot_df["bic_clip"] = plot_df["bm_dbic"].clip(-20, 120)
    for tp, col, lab in [(1, "#2166ac", "TP"), (0, "#b2182b", "FP")]:
        sub = plot_df[plot_df["is_tp"] == tp]
        sub_sub = sub[(sub["is_bimodal"] == True) & (sub["hp_split_alignment"] > HP_SPLIT_TOL)]
        sub_hp = sub[~((sub["is_bimodal"] == True)
                       & (sub["hp_split_alignment"] > HP_SPLIT_TOL))]
        gsA.scatter(sub_hp["blind_ari"], sub_hp["bic_clip"], s=26, c=col, alpha=0.45,
                    edgecolors="none", label=f"{lab} (uni/HP-split)")
        gsA.scatter(sub_sub["blind_ari"], sub_sub["bic_clip"], s=85, marker="*", c=col,
                    edgecolors="k", linewidths=0.5, label=f"{lab} subclone-like")
    gsA.axhline(10, color="grey", ls="--", lw=0.8)
    gsA.axvline(ARI_HI, color="grey", ls=":", lw=0.8)
    gsA.text(ARI_HI + 0.01, -18, f"blind_ari median={ARI_HI:.2f}", fontsize=8, color="grey")
    gsA.text(0.62, 12, "bimodal threshold dBIC=10", fontsize=8, color="grey")
    gsA.set_xlabel("Read-level blind methylation clustering (blind_ari)")
    gsA.set_ylabel("Per-read bimodality strength (GMM delta-BIC, clipped)")
    gsA.set_title("(A) Clustering vs per-read bimodality\n(star = subclone-like: modes NOT "
                  "aligned to germ/som HP split)", fontsize=10)
    gsA.legend(fontsize=7, loc="upper left", framealpha=0.9)

    # (B) live per-read beta histograms for top candidates + BRCA2
    pick = []
    for c in top_candidates:
        ch, ps = c["locus"].split(":")
        pick.append((ch, int(ps), c["axis"], c["label"] or c["locus"], c["bm_means"]))
        if len(pick) >= 4:
            break
    # always include BRCA2 as reference
    pick = pick[:3] + [("chr13", 32315128, "HP1_vs_HP1-1", "BRCA2 (ref)",
                        brca2["bm_means"])]

    gsB = fig.add_axes([0.54, 0.12, 0.43, 0.78])
    offsets = np.arange(len(pick)) * 1.15
    colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(pick)))
    try:
        bam = pysam.AlignmentFile(BAM, "rb")
        for i, (ch, ps, axis, lab, means) in enumerate(pick):
            betas = fetch_perread_beta(bam, ch, ps, win=600)
            if len(betas) < 5:
                continue
            hist, edges = np.histogram(betas, bins=np.linspace(0, 1, 26), density=True)
            centers = (edges[:-1] + edges[1:]) / 2
            gsB.fill_between(centers, offsets[i], offsets[i] + hist / (hist.max() + 1e-9) * 0.95,
                             color=colors[i], alpha=0.7, step="mid")
            gsB.plot(centers, offsets[i] + hist / (hist.max() + 1e-9) * 0.95,
                     color="k", lw=0.5)
            if means:
                for m in means:
                    gsB.plot([m, m], [offsets[i], offsets[i] + 0.95], color="red",
                             ls="--", lw=0.9)
            gsB.text(0.01, offsets[i] + 0.55, f"{lab}\n{ch}:{ps} (n={len(betas)})",
                     fontsize=7.5, va="center")
        bam.close()
    except Exception as e:
        gsB.text(0.5, 0.5, f"BAM fetch failed:\n{e}", ha="center", transform=gsB.transAxes)
    gsB.set_yticks([])
    gsB.set_xlabel("Per-read mean methylation beta (continuous, 5mC|5hmC max)")
    gsB.set_title("(B) Per-read beta distributions (live BAM)\n"
                  "red dashes = GMM component means", fontsize=10)
    gsB.set_xlim(0, 1)

    fig.suptitle("S - Subclone assessment: read-level clustering + per-read bimodality "
                 "(HCC1395, A pilot)", fontsize=12)
    fig.savefig(OUT_FIG, dpi=140)
    plt.close(fig)

    # ---------------------------------------------------------------------
    # VERDICT synthesis
    # ---------------------------------------------------------------------
    n_subclone_candidates = int(cand["is_candidate"].sum())  # TP & bimodal & subclone-like
    n_cand_ari_confirmed = int((cand["is_candidate"]
                                & (cand["blind_ari"] >= ARI_HI)).sum())
    frac_subclone_of_bim = (n_bim_subclone_like / n_bim) if n_bim else 0

    # Strong MIXED_PARTIAL requires the bimodal subclone-like loci to be CORROBORATED by the
    # independent blind_ari clustering signal. Here corroboration n=0 (disjoint sampling), so the
    # subclone-like label rests only on the GMM-weight heuristic -> stay conservative.
    if n_cand_ari_confirmed >= 3 and n_subclone_candidates >= 3 and frac_subclone_of_bim >= 0.4:
        verdict = ("MIXED_PARTIAL: per-read bimodality + read clustering JOINTLY flag a minority "
                   "of loci (n=%d TP, subclone-like, %d blind_ari-corroborated) with sub-structure "
                   "NOT explained by the germ-vs-som HP split, but most bimodality co-localizes "
                   "with the HP split and FP power is too low to claim a TP/FP-discriminating "
                   "subclone signal." % (n_subclone_candidates, n_cand_ari_confirmed))
    elif n_subclone_candidates >= 1:
        verdict = ("PREDOMINANTLY_HP_SPLIT: per-read bimodality is largely the germline-group vs "
                   "somatic-subgroup HP split (%d of %d bimodal loci); a minority (n=%d TP) show "
                   "extra subclone-like sub-structure where the minor mode does not match the "
                   "germ/som read fraction -- low-confidence (single-sample, only %d of these have "
                   "blind_ari clustering corroboration)."
                   % (n_bim_hp_explained, n_bim, n_subclone_candidates, n_cand_ari_confirmed))
    else:
        verdict = ("HP_SPLIT_ONLY: bimodality tracks the germ-vs-som HP split; no robust evidence "
                   "of a methylation-distinct minor subclone read subpopulation.")

    out = {
        "meta": {
            "script": "52_subclone_assessment.py", "sample": "HCC1395",
            "task_type": "A pilot", "axis_note": "HP-axis = somatic-controlled",
            "inputs": [M2_PATH, C2_PATH, O2_PATH], "bam": BAM,
            "hp_split_tol": HP_SPLIT_TOL, "ari_hi_threshold": ARI_HI,
            "subclone_discriminator": (
                "bimodal per-read beta where GMM minor-mode weight does NOT match germ-or-som "
                "read fraction (within HP_SPLIT_TOL) => structure beyond the 2 HP groups."),
        },
        "q1_subclone_vs_hpsplit": {
            "n_loci_bimodal_testable": int(msum["n_loci_testable_bimodal"]),
            "n_bimodal": n_bim,
            "n_bimodal_hp_split_explained": n_bim_hp_explained,
            "n_bimodal_subclone_like": n_bim_subclone_like,
            "frac_bimodal_subclone_like": round(frac_subclone_of_bim, 3),
            "bimodal_blind_ari_median": blind_med,
            "bimodal_placebo_ari_median": plac_med,
            "blind_vs_placebo_global": "blind_ari>>placebo (C2 H-E p=1.8e-43): clustering is real",
            "bimodal_tp_vs_fp": {
                "tp_frac": msum["bimodal_by_tp"]["tp=1"]["frac"],
                "fp_frac": msum["bimodal_by_tp"]["tp=0"]["frac"],
                "fp_n_bimodal": msum["bimodal_by_tp"]["tp=0"]["bimodal"],
                "note": "FP n=16 testable -> bimodality cannot discriminate TP/FP (underpowered)",
            },
            "interpretation": (
                "Most per-read bimodality (n=%d of %d) aligns with the germ-vs-som HP split; "
                "%d loci show a minor mode whose weight deviates from the germ/som fraction "
                "(subclone-like sub-structure). blind_ari>>placebo confirms the clustering is "
                "real read structure (not length/coverage artifact), but the dominant structure "
                "IS the somatic HP subgroup, not a novel minor clone."
                % (n_bim_hp_explained, n_bim, n_bim_subclone_like)),
        },
        "q2_measurement_hides_subclone": {
            "binarized_vs_continuous_corr": bin_cont_corr,
            "corr_db_max_vs_5hmC": corr_max_5hmc,
            "frac_5hmC_gt_5mC": frac_5hmc_gt_5mc,
            "deadzone_frac": deadzone_frac,
            "n_bimodal_with_high_intermediate": n_hidden,
            "top_divergent_loci_bimodal_overlap": topdiv_bim.to_dict("records"),
            "measurement_hides_subclone": measurement_hides,
            "why": (
                "Two distinct hiding mechanisms: "
                "(1) BINARIZATION: global binarized~continuous corr is tight (0.998), so for the "
                "AVERAGE locus max-collapse is safe; BUT %d bimodal loci carry a sizeable "
                "intermediate-methylation read fraction (>0.15) that binarization crushes to 0/1, "
                "and the top-divergent loci (div up to 0.13, incl BRCA2) are exactly where a "
                "partially-methylated subpopulation is flattened -> continuous beta DOES recover "
                "structure binarized hides at these specific loci. "
                "(2) MAX-COLLAPSE of 5mC|5hmC: corr(db_max, db_5hmC)=0.033 (orthogonal) with "
                "5hmC>5mC in %.1f%% of CpGs -> the max-collapse measurement is 5mC-dominated and "
                "DISCARDS the 5hmC channel entirely; any 5hmC-defined subpopulation is invisible "
                "to the current binarized-max measurement."
                % (n_hidden, frac_5hmc_gt_5mc * 100)),
        },
        "q3_top_candidates": {
            "n_candidates": n_subclone_candidates,
            "n_candidates_blind_ari_confirmed": n_cand_ari_confirmed,
            "m2_c2_locus_overlap": {
                "m2_testable": n_m2_testable, "c2_bam_pass": int(o2.shape[0]),
                "n_with_blind_ari": n_with_ari, "n_bimodal_with_blind_ari": n_bimodal_with_ari,
                "note": ("M2 (bimodality) and C2 (blind_ari) sampled largely DISJOINT loci; "
                         "only ~15 overlap, 4 of 43 bimodal loci have blind_ari. Candidates are "
                         "ranked by the DIRECT per-read bimodality test and annotated with "
                         "blind_ari where available -- NOT gated on the near-empty intersection."),
            },
            "ranking_criteria": "TP & is_bimodal & subclone-like(minor mode NOT germ/som HP-split); "
                                "score = bm_sep*minor_mode_w*ari_factor "
                                "(1.5 if blind_ari>=median, 1.0 if absent, 0.6 if low)",
            "top_candidate_loci": top_candidates,
            "brca2": brca2,
        },
        "verdict": {
            "subclone_evidence_verdict": verdict,
            "n_subclone_candidates": n_subclone_candidates,
            "n_subclone_candidates_blind_ari_confirmed": n_cand_ari_confirmed,
            "measurement_hides_subclone": measurement_hides,
            "measurement_hides_why": (
                "Binarized-max is safe on average (r=0.998) but (a) crushes intermediate-meth "
                "subpopulations at specific bimodal/high-frac_intermediate loci incl BRCA2, and "
                "(b) discards the orthogonal 5hmC channel (corr 0.033). Separated continuous "
                "5mC + 5hmC would expose structure the current measurement hides at these loci."),
            "brca2_qualifies": brca2_qualifies,
            "recommendation": (
                "Do NOT claim a discrete methylation subclone from this single-sample HCC1395 "
                "pilot: the dominant read structure is the somatic HP subgroup (germ-vs-som "
                "split), bimodality does not discriminate TP/FP (FP underpowered, n=16 testable), "
                "and subclone-like loci are a minority. HOWEVER the measurement IS lossy at "
                "specific loci -> (1) switch from binarized-max to SEPARATED continuous 5mC and "
                "5hmC beta for the %d top-divergent / high-frac_intermediate loci (incl BRCA2); "
                "(2) treat the %d subclone-like loci as hypotheses to confirm in a second sample "
                "(COLO829) with per-read GMM on continuous beta SEPARATED by mod-type, before any "
                "subclone claim." % (len(msum["top_divergent_loci"]), n_subclone_candidates)),
        },
        "power_caveats": {
            "single_sample": "HCC1395 only; subclone claims need >=2 samples",
            "fp_power": "43 FP in 596 BAM-pass loci; 16 testable for bimodality -> no TP/FP claim",
            "pseudoreplication": "per-position; HP1/HP2 records of same locus correlated",
        },
        "paths": {"json": OUT_JSON, "figure": OUT_FIG},
    }
    json.dump(out, open(OUT_JSON, "w"), indent=2, default=_json_default)

    # console summary
    print("=== S subclone assessment ===")
    print("verdict:", verdict)
    print(f"n_bimodal={n_bim} | HP-split-explained={n_bim_hp_explained} | "
          f"subclone-like={n_bim_subclone_like} | n_subclone_candidates={n_subclone_candidates}")
    print("measurement_hides_subclone:", measurement_hides)
    print("BRCA2 qualifies subclone-candidate:", brca2_qualifies,
          "| is_bimodal:", brca2["is_bimodal"], "| blind_ari:", brca2["blind_ari"])
    print(f"n_subclone_candidates_blind_ari_confirmed={n_cand_ari_confirmed}")
    print("top candidate loci (TP, bimodal, subclone-like):")
    for c in top_candidates[:8]:
        ari = "%.3f" % c["blind_ari"] if c["blind_ari"] is not None else "NA"
        print("  ", c["locus"], c["axis"], "tp=%s" % c["is_tp"], c["loh"], c["cn_class"],
              "blind_ari=%s" % ari, "sep=%s" % c["bm_sep"], "minorW=%s" % c["minor_mode_w"],
              "score=%s" % c["subclone_score"])
    print("JSON:", OUT_JSON)
    print("FIG :", OUT_FIG)
    return out


def fetch_perread_beta(bam, chrom, pos, win=600):
    """Per-read mean methylation beta over [pos-win, pos+win], 5mC|5hmC MAX-collapse, P>=0.5."""
    betas = []
    start, end = max(0, pos - win), pos + win
    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        mods = read.modified_bases
        if not mods:
            continue
        # collect per-position max ML across 5mC(+m)/5hmC(+h)
        best = {}
        for (base, strand, modname), data in mods.items():
            if modname not in ("m", "h", "C+m", "C+h", 27551):
                # accept '+m' / '+h' tag forms
                if str(modname) not in ("m", "h"):
                    continue
            for qpos, ml in data:
                best[qpos] = max(best.get(qpos, 0), ml)
        if not best:
            continue
        calls = [1 if ml >= 128 else 0 for ml in best.values()]  # P>=0.5 (ML 0-255)
        if calls:
            betas.append(np.mean(calls))
    return np.array(betas)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if pd.isna(o):
        return None
    return str(o)


if __name__ == "__main__":
    main()
