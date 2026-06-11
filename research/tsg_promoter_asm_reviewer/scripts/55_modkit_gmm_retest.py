#!/usr/bin/env python3
"""
55_modkit_gmm_retest.py  (P2 GMM re-test on modkit SEPARATED CONTINUOUS beta)
=============================================================================

HCC1395 single sample, A pilot. Gate already passed (54_modkit_extract_crossval.py
verdict=AGREE; modkit<->pysam per-read 5mC Pearson=1.0).

GOAL
----
M2 (49_meth_measurement_audit.py) tested per-read methylation bimodality using
pysam-binarized **MAX-collapse** of (5mC, 5hmC) ML -> single channel. This script
re-tests the SAME candidate loci on modkit's per-read CONTINUOUS `mod_qual`, with
5mC and 5hmC SEPARATED, asking:
  (b) does separating the two mod channels (and using continuous prob instead of a
      MAX-collapse) GAIN / LOSE bimodality vs M2?
  (c) does a 5hmC-defined subpopulation exist that MAX-collapse hides?
  (d) for candidates with blind_ari, does GMM-bimodal agree with the independent
      blind clustering ARI?

MEASUREMENT (modkit continuous, separated)
------------------------------------------
- Input: modkit `extract full` table (modkit_extract.tsv); cols read_id / ref_position
  / chrom / mod_qual (continuous prob in [0,1]) / mod_code (m=5mC | h=5hmC, SEPARATE rows).
- modkit output does NOT carry the HP tag -> join read_id -> HP via pysam over the BAM
  (same `hp_tag` convention as 49_meth_measurement_audit.py: read the "HP" aux tag,
  string values '1' / '1-1' / '2' / '2-1' ...).
- For each candidate (axis = e.g. HP1_vs_HP1-1), the SOMATIC HP SUBGROUP is the RIGHT
  side of the axis (gs, e.g. '1-1'). Per-read mean 5mC beta = mean of mod_qual over
  mod_code=='m' rows of that read within the locus +/-WINDOW; per-read mean 5hmC beta
  = same over mod_code=='h'. (mod_qual already [0,1]; no /255 rescale needed.)

GMM (identical thresholds to M2 / 49_meth_measurement_audit.gmm_bimodality)
---------------------------------------------------------------------------
1-vs-2 component GaussianMixture BIC. bimodal iff:
  delta_BIC = bic1 - bic2 >= 10  AND  min component weight >= 0.10  AND  mean_sep >= 0.10
Run on (i) per-read 5mC-continuous, (ii) per-read 5hmC-continuous,
(iii) joint 2D [5mC, 5hmC].

COMPARISON
----------
M2 verdict = m2_meth_measurement_audit.json -> per_locus[...].bimodal.is_bimodal
(pysam-binarized MAX-collapse). "modkit_continuous bimodal" = ANY of (i)/(ii)/(iii)
bimodal under the better (continuous, separated) measurement.
n_5hmC_new_structure = loci where 5hmC-channel (or joint via 5hmC axis) is bimodal but
M2 MAX-collapse was NOT -> a 5hmC-defined subpopulation invisible to MAX-collapse.

Output: gmm_retest.json + modkit_gmm_retest.png. Single sample -> A pilot tier, no
cross-sample claim.
"""
import os, json, csv
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from sklearn.mixture import GaussianMixture

os.environ.setdefault("TMPDIR", "/big7_disk")

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
TRIAL = f"{ROOT}/genome_survey_v2/cn_confound/modkit_trial"
EXTRACT = f"{TRIAL}/modkit_extract.tsv"
M2_JSON = f"{ROOT}/genome_survey_v2/cn_confound/m2_meth_measurement_audit.json"
SUBCLONE = f"{ROOT}/genome_survey_v2/cn_confound/s_subclone_assessment.json"
CROSSVAL = f"{TRIAL}/crossval.json"
MASTER_O2 = f"{ROOT}/genome_survey_v2/cn_confound/master_o2_error.tsv"
BAM = ("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
       "20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam")
OUT_JSON = f"{TRIAL}/gmm_retest.json"
FIG = f"{ROOT}/figures/cn_confound/modkit_gmm_retest.png"

WINDOW = 600              # +/- bp, matches 49/54
MIN_CPG_PER_READ = 1      # any-mod CpG present (modkit emits both m/h per CpG; count CpGs)
GMM_DBIC = 10.0
GMM_MINW = 0.10
GMM_SEP = 0.10
MIN_READS_GMM = 12        # matches 49 n<12 guard

# ---------------------------------------------------------------- candidates
subc = json.load(open(SUBCLONE))
cand_loci = subc["q3_top_candidates"]["top_candidate_loci"]
# +BRCA2
BRCA2 = dict(locus="chr13:32315128", axis="HP1_vs_HP1-1", label="BRCA2",
             is_tp=1, beta_cont=None, blind_ari=None)

def axis_to_groups(axis):
    # 'HP1_vs_HP1-1' -> ('1','1-1'); 'HP2_vs_HP2-1' -> ('2','2-1')
    left, right = axis.split("_vs_")
    return left.replace("HP", ""), right.replace("HP", "")

cands = []
seen = set()
for c in cand_loci:
    ch, pos = c["locus"].split(":")
    pos = int(pos)
    gg, gs = axis_to_groups(c["axis"])
    cands.append(dict(chrom=ch, pos=pos, axis=c["axis"], gg=gg, gs=gs,
                      label=c.get("label"), is_tp=int(c["is_tp"]),
                      beta_cont=c.get("beta_cont"), blind_ari=c.get("blind_ari")))
    seen.add((ch, pos))
# add BRCA2 if not already present
bc, bp = BRCA2["locus"].split(":"); bp = int(bp)
if (bc, bp) not in seen:
    gg, gs = axis_to_groups(BRCA2["axis"])
    cands.append(dict(chrom=bc, pos=bp, axis=BRCA2["axis"], gg=gg, gs=gs,
                      label="BRCA2", is_tp=1, beta_cont=None, blind_ari=None))

# ---- attach blind_ari from master_o2_error.tsv (somatic_pos + axis match) -----
ari_idx = {}
with open(MASTER_O2) as f:
    rd = csv.DictReader(f, delimiter="\t")
    for row in rd:
        try:
            key = (row["chrom"], int(row["somatic_pos"]), row["axis"])
        except Exception:
            continue
        ba = row.get("blind_ari", "")
        try:
            ari_idx[key] = float(ba) if ba not in ("", "None", "nan") else None
        except Exception:
            ari_idx[key] = None
for c in cands:
    c["blind_ari"] = ari_idx.get((c["chrom"], c["pos"], c["axis"]), c.get("blind_ari"))

# ------------------------------------------------ M2 pysam-binarized verdicts
m2 = json.load(open(M2_JSON))
m2_idx = {(r["chrom"], r["pos"]): r for r in m2["per_locus"]}

# ------------------------------------------------ read_id -> HP via pysam
print("[hp-join] building read_id -> HP map over candidate windows ...")
bam = pysam.AlignmentFile(BAM, "rb")

def hp_tag(r):
    for t, v in r.tags:
        if t == "HP":
            return str(v)
    return None

read_hp = {}   # read_id -> HP string
for c in cands:
    var0 = c["pos"] - 1
    s, e = var0 - WINDOW, var0 + WINDOW
    for r in bam.fetch(c["chrom"], max(0, s), e):
        if r.flag & 0x900 or r.flag & 0x400:   # secondary/supplementary/dup
            continue
        h = hp_tag(r)
        if h is None:
            continue
        read_hp[r.query_name] = h
print(f"[hp-join] {len(read_hp)} primary reads tagged with HP across candidate windows")

# ------------------------------------------------ load modkit extract table
# group rows by (chrom, read_id) -> {'m':[(refpos,q)], 'h':[(refpos,q)]}
print("[modkit] loading extract table ...")
per_read = {}   # (chrom, read_id) -> {'m':[...], 'h':[...]}
with open(EXTRACT) as f:
    rd = csv.DictReader(f, delimiter="\t")
    for row in rd:
        ch = row["chrom"]; rid = row["read_id"]
        code = row["mod_code"]
        if code not in ("m", "h"):
            continue
        try:
            refpos = int(row["ref_position"]); q = float(row["mod_qual"])
        except Exception:
            continue
        d = per_read.setdefault((ch, rid), {"m": [], "h": []})
        d[code].append((refpos, q))
print(f"[modkit] {len(per_read)} (chrom,read) groups")

# ------------------------------------------------ GMM helper (M2-identical)
def gmm_1d(x):
    n = len(x)
    res = dict(n=int(n), is_bimodal=False, bic1=None, bic2=None, delta_bic=None,
               comp_means=None, comp_weights=None, mean_sep=None, min_weight=None)
    if n < MIN_READS_GMM:
        res["reason"] = "n<%d" % MIN_READS_GMM
        return res
    X = np.asarray(x, float).reshape(-1, 1)
    if np.nanstd(X) < 1e-6:
        res["reason"] = "degenerate"
        return res
    try:
        g1 = GaussianMixture(1, covariance_type="full", random_state=42, n_init=2).fit(X)
        g2 = GaussianMixture(2, covariance_type="full", random_state=42, n_init=5).fit(X)
    except Exception as ex:
        res["reason"] = "gmm_fail:%s" % ex
        return res
    b1, b2 = g1.bic(X), g2.bic(X)
    means = sorted(g2.means_.ravel().tolist())
    w = g2.weights_.ravel()
    minw = float(min(w)); sep = float(abs(means[1] - means[0]))
    res.update(bic1=float(b1), bic2=float(b2), delta_bic=float(b1 - b2),
               comp_means=[round(m, 4) for m in means],
               comp_weights=[round(float(z), 4) for z in sorted(w.tolist())],
               mean_sep=round(sep, 4), min_weight=round(minw, 4))
    res["is_bimodal"] = bool((b1 - b2) >= GMM_DBIC and minw >= GMM_MINW and sep >= GMM_SEP)
    return res

def gmm_2d(xy):
    n = len(xy)
    res = dict(n=int(n), is_bimodal=False, bic1=None, bic2=None, delta_bic=None,
               comp_means=None, comp_weights=None, mean_sep=None, min_weight=None)
    if n < MIN_READS_GMM:
        res["reason"] = "n<%d" % MIN_READS_GMM
        return res
    X = np.asarray(xy, float)
    if np.nanstd(X[:, 0]) < 1e-6 and np.nanstd(X[:, 1]) < 1e-6:
        res["reason"] = "degenerate"
        return res
    try:
        g1 = GaussianMixture(1, covariance_type="full", random_state=42, n_init=2).fit(X)
        g2 = GaussianMixture(2, covariance_type="full", random_state=42, n_init=5).fit(X)
    except Exception as ex:
        res["reason"] = "gmm_fail:%s" % ex
        return res
    b1, b2 = g1.bic(X), g2.bic(X)
    mu = g2.means_              # 2x2
    w = g2.weights_.ravel()
    minw = float(min(w))
    sep = float(np.linalg.norm(mu[0] - mu[1]))   # euclidean separation in [5mC,5hmC]
    res.update(bic1=float(b1), bic2=float(b2), delta_bic=float(b1 - b2),
               comp_means=[[round(float(mu[i, 0]), 4), round(float(mu[i, 1]), 4)]
                           for i in range(2)],
               comp_weights=[round(float(z), 4) for z in sorted(w.tolist())],
               mean_sep=round(sep, 4), min_weight=round(minw, 4))
    res["is_bimodal"] = bool((b1 - b2) >= GMM_DBIC and minw >= GMM_MINW and sep >= GMM_SEP)
    return res

# ------------------------------------------------ per-candidate computation
results = []
scatter_store = {}   # locus -> (mC_arr, hC_arr, gmm verdicts) for figure

for c in cands:
    ch, pos, gs = c["chrom"], c["pos"], c["gs"]
    var0 = pos - 1
    s, e = var0 - WINDOW, var0 + WINDOW
    # collect per-read mean 5mC / 5hmC for reads in the SOMATIC HP subgroup (gs)
    mC, hC, both = [], [], []
    n_reads_som = 0
    for (rch, rid), d in per_read.items():
        if rch != ch:
            continue
        if read_hp.get(rid) != gs:           # somatic HP subgroup only
            continue
        m_in = [q for (rp, q) in d["m"] if s <= rp <= e]
        h_in = [q for (rp, q) in d["h"] if s <= rp <= e]
        if len(m_in) < MIN_CPG_PER_READ and len(h_in) < MIN_CPG_PER_READ:
            continue
        n_reads_som += 1
        mm = float(np.mean(m_in)) if m_in else np.nan
        hh = float(np.mean(h_in)) if h_in else np.nan
        mC.append(mm); hC.append(hh)
        if not np.isnan(mm) and not np.isnan(hh):
            both.append((mm, hh))
    mC = np.array(mC); hC = np.array(hC)
    mC_valid = mC[~np.isnan(mC)]
    hC_valid = hC[~np.isnan(hC)]
    both_arr = np.array(both) if both else np.empty((0, 2))

    g_mC = gmm_1d(mC_valid)
    g_hC = gmm_1d(hC_valid)
    g_2d = gmm_2d(both_arr)

    modkit_bimodal = bool(g_mC["is_bimodal"] or g_hC["is_bimodal"] or g_2d["is_bimodal"])

    # M2 pysam-binarized-MAX verdict
    m2r = m2_idx.get((ch, pos))
    m2_bimodal = bool(m2r["bimodal"]["is_bimodal"]) if m2r else None
    m2_dbic = m2r["bimodal"].get("delta_bic") if m2r else None

    changed = (m2_bimodal is not None) and (modkit_bimodal != m2_bimodal)
    gained = (m2_bimodal is False) and modkit_bimodal
    lost = (m2_bimodal is True) and (not modkit_bimodal)
    # 5hmC-defined NEW structure: 5hmC channel (or joint) bimodal AND MAX-collapse (M2) was not,
    # AND it is NOT merely re-discovering the 5mC bimodality
    hmc_new = bool((g_hC["is_bimodal"] or g_2d["is_bimodal"])
                   and (m2_bimodal is False)
                   and (not g_mC["is_bimodal"]))

    rec = dict(
        chrom=ch, pos=pos, locus=f"{ch}:{pos}", label=c.get("label"),
        axis=c["axis"], somatic_hp=gs, is_tp=c["is_tp"], beta_cont=c.get("beta_cont"),
        blind_ari=c.get("blind_ari"),
        n_reads_somatic=n_reads_som,
        n_with_5mC=int(len(mC_valid)), n_with_5hmC=int(len(hC_valid)),
        n_with_both=int(len(both_arr)),
        mean_5mC=round(float(np.nanmean(mC)), 4) if len(mC_valid) else None,
        mean_5hmC=round(float(np.nanmean(hC)), 4) if len(hC_valid) else None,
        gmm_5mC=g_mC, gmm_5hmC=g_hC, gmm_joint2d=g_2d,
        modkit_continuous_bimodal=modkit_bimodal,
        m2_pysam_binarized_bimodal=m2_bimodal, m2_delta_bic=m2_dbic,
        verdict_changed_vs_pysam=changed, gained_bimodal=gained, lost_bimodal=lost,
        hmc_new_structure=hmc_new,
    )
    results.append(rec)
    scatter_store[f"{ch}:{pos}"] = (mC, hC, g_mC, g_hC, g_2d, rec)
    print(f"  {ch}:{pos:>10} {c.get('label') or '':8} hp={gs:4} "
          f"n_som={n_reads_som:3} | mod-bimodal={modkit_bimodal} "
          f"(5mC={g_mC['is_bimodal']},5hmC={g_hC['is_bimodal']},2D={g_2d['is_bimodal']}) "
          f"| M2={m2_bimodal} | changed={changed} 5hmC_new={hmc_new}")

# ------------------------------------------------ summary
n_tested = len(results)
n_mod_bimodal = sum(r["modkit_continuous_bimodal"] for r in results)
n_changed = sum(r["verdict_changed_vs_pysam"] for r in results)
n_gained = sum(r["gained_bimodal"] for r in results)
n_lost = sum(r["lost_bimodal"] for r in results)
n_hmc_new = sum(r["hmc_new_structure"] for r in results)
brca2 = next((r for r in results if r["label"] == "BRCA2"), None)

# top confirmed = bimodal under BOTH modkit-continuous AND M2 pysam (robust to measurement)
top_confirmed = [r["locus"] + (f"({r['label']})" if r["label"] else "")
                 for r in results
                 if r["modkit_continuous_bimodal"] and r["m2_pysam_binarized_bimodal"]]

# blind_ari cross-check (d)
ari_rows = [r for r in results if r["blind_ari"] is not None]
ari_check = [dict(locus=r["locus"], label=r["label"], blind_ari=r["blind_ari"],
                  modkit_bimodal=r["modkit_continuous_bimodal"],
                  m2_bimodal=r["m2_pysam_binarized_bimodal"])
             for r in ari_rows]

brca2_verdict = None
if brca2:
    parts = []
    parts.append("5mC=%s" % brca2["gmm_5mC"]["is_bimodal"])
    parts.append("5hmC=%s" % brca2["gmm_5hmC"]["is_bimodal"])
    parts.append("2D=%s" % brca2["gmm_joint2d"]["is_bimodal"])
    brca2_verdict = ("BIMODAL" if brca2["modkit_continuous_bimodal"] else "UNIMODAL") + \
        f" (modkit-continuous; {','.join(parts)}; M2-pysam-MAX={brca2['m2_pysam_binarized_bimodal']})"

summary = dict(
    n_candidates_tested=n_tested,
    n_bimodal_modkit_continuous=n_mod_bimodal,
    n_bimodal_changed_vs_pysam=n_changed,
    n_gained_bimodal=n_gained, n_lost_bimodal=n_lost,
    n_5hmC_new_structure=n_hmc_new,
    brca2_modkit_verdict=brca2_verdict,
    top_confirmed_candidates=top_confirmed,
    n_with_blind_ari=len(ari_rows),
    blind_ari_crosscheck=ari_check,
)

out = dict(
    meta=dict(
        script="55_modkit_gmm_retest.py", sample="HCC1395", task_type="A pilot",
        bam=BAM, modkit_extract=EXTRACT, window=WINDOW,
        measurement="modkit extract full per-read CONTINUOUS mod_qual [0,1]; "
                    "5mC(m)/5hmC(h) SEPARATE; per-read mean over CpG within +/-WINDOW; "
                    "somatic HP subgroup = right side of axis (gs)",
        hp_join="read_id -> HP via pysam 'HP' aux tag over candidate windows "
                "(same convention as 49_meth_measurement_audit.py)",
        gmm=f"1-vs-2 GMM BIC; bimodal iff dBIC>={GMM_DBIC} & min_weight>={GMM_MINW} & "
            f"sep>={GMM_SEP} (1D sep=|mean diff|, 2D sep=euclidean); identical thresholds to M2",
        m2_baseline=M2_JSON + " (pysam-binarized MAX-collapse, ML>=128)",
        gate="54_modkit_extract_crossval verdict=AGREE (modkit<->pysam 5mC Pearson=1.0)",
        note="Single sample HCC1395 -> A pilot; NO cross-sample generalization claim",
    ),
    summary=summary,
    per_locus=results,
)
with open(OUT_JSON, "w") as f:
    json.dump(out, f, indent=2)
print("\n[write]", OUT_JSON)

# ------------------------------------------------ FIGURE
# Layout: top = per-read 5mC vs 5hmC scatter w/ GMM for top candidates incl BRCA2
#         bottom = verdict comparison table (modkit-continuous vs pysam-binarized)
order = sorted(results, key=lambda r: (r["label"] != "BRCA2",
                                       -(r["n_reads_somatic"])))
# pick BRCA2 + up to 7 others with most somatic reads
panel_loci = []
for r in order:
    panel_loci.append(r["locus"])
    if len(panel_loci) >= 8:
        break

ncol = 4
nrow_scatter = 2
fig = plt.figure(figsize=(18, 13))
gs_fig = fig.add_gridspec(nrow_scatter + 2, ncol, height_ratios=[1, 1, 0.12, 1.4])

def draw_ellipse(ax, mean, cov, color):
    try:
        v, w = np.linalg.eigh(cov)
        v = 2.0 * np.sqrt(np.maximum(v, 1e-9))
        ang = np.degrees(np.arctan2(w[1, 0], w[0, 0]))
        e = Ellipse(mean, v[0], v[1], angle=ang, fc="none", ec=color, lw=2, ls="--")
        ax.add_patch(e)
    except Exception:
        pass

for i, loc in enumerate(panel_loci):
    mC, hC, g_mC, g_hC, g_2d, rec = scatter_store[loc]
    rr = i // ncol; cc = i % ncol
    ax = fig.add_subplot(gs_fig[rr, cc])
    m = ~np.isnan(mC) & ~np.isnan(hC)
    ax.scatter(mC[m], hC[m], s=22, alpha=0.6,
               c=("#d7191c" if rec["label"] == "BRCA2" else "#2c7fb8"),
               edgecolors="none")
    # 2D GMM ellipses if bimodal
    if g_2d.get("comp_means") and len(rec) and g_2d["n"] >= MIN_READS_GMM and np.sum(m) >= MIN_READS_GMM:
        try:
            X = np.column_stack([mC[m], hC[m]])
            gg = GaussianMixture(2, covariance_type="full", random_state=42, n_init=5).fit(X)
            for k in range(2):
                draw_ellipse(ax, gg.means_[k], gg.covariances_[k], "#222222")
        except Exception:
            pass
    title = (rec["label"] or loc)
    tag = []
    if g_mC["is_bimodal"]: tag.append("5mC*")
    if g_hC["is_bimodal"]: tag.append("5hmC*")
    if g_2d["is_bimodal"]: tag.append("2D*")
    ax.set_title(f"{title}\n{loc} hp={rec['somatic_hp']} n={rec['n_reads_somatic']} "
                 f"[{','.join(tag) if tag else 'unimodal'}]", fontsize=9)
    ax.set_xlabel("per-read mean 5mC (modkit prob)", fontsize=8)
    ax.set_ylabel("per-read mean 5hmC (modkit prob)", fontsize=8)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, max(0.3, np.nanmax(hC) * 1.15 if len(hC) else 0.3))
    ax.tick_params(labelsize=7)

# verdict table spanning bottom
axt = fig.add_subplot(gs_fig[3, :])
axt.axis("off")
cols = ["locus", "label", "hp", "n_som", "modkit\n5mC", "modkit\n5hmC", "modkit\n2D",
        "modkit\nANY", "M2 pysam\nMAX", "changed", "5hmC\nNEW", "blind_ari"]
tbl = []
for r in sorted(results, key=lambda x: (x["label"] != "BRCA2", x["chrom"], x["pos"])):
    tbl.append([
        r["locus"], r["label"] or "-", r["somatic_hp"], r["n_reads_somatic"],
        "Y" if r["gmm_5mC"]["is_bimodal"] else "·",
        "Y" if r["gmm_5hmC"]["is_bimodal"] else "·",
        "Y" if r["gmm_joint2d"]["is_bimodal"] else "·",
        "BIMOD" if r["modkit_continuous_bimodal"] else "uni",
        "BIMOD" if r["m2_pysam_binarized_bimodal"] else "uni",
        "CHG" if r["verdict_changed_vs_pysam"] else "",
        "NEW" if r["hmc_new_structure"] else "",
        f"{r['blind_ari']:.3f}" if r["blind_ari"] is not None else "-",
    ])
table = axt.table(cellText=tbl, colLabels=cols, loc="center", cellLoc="center")
table.auto_set_font_size(False); table.set_fontsize(8); table.scale(1, 1.5)
for j in range(len(cols)):
    table[0, j].set_facecolor("#dddddd"); table[0, j].set_text_props(weight="bold")
for ri, r in enumerate(sorted(results, key=lambda x: (x["label"] != "BRCA2", x["chrom"], x["pos"])), start=1):
    if r["verdict_changed_vs_pysam"]:
        for j in range(len(cols)):
            table[ri, j].set_facecolor("#fff2cc")
    if r["label"] == "BRCA2":
        for j in range(len(cols)):
            table[ri, j].set_facecolor("#fde0dc")

fig.suptitle("M55 - modkit SEPARATED CONTINUOUS GMM re-test vs M2 pysam-binarized MAX-collapse "
             f"(HCC1395, A pilot)\n"
             f"tested={n_tested}  modkit-bimodal={n_mod_bimodal}  changed-vs-pysam={n_changed} "
             f"(gain={n_gained}/lost={n_lost})  5hmC-NEW={n_hmc_new}",
             fontsize=13, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.97])
os.makedirs(os.path.dirname(FIG), exist_ok=True)
fig.savefig(FIG, dpi=130, bbox_inches="tight")
print("[write]", FIG)

print("\n=== SUMMARY ===")
print(json.dumps(summary, indent=2))
