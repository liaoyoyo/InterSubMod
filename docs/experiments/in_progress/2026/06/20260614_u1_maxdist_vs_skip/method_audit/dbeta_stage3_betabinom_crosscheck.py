#!/usr/bin/env python3
"""Δβ module stage-3: beta-binomial / GLM cross-check vs permutation (DSS interaction alignment).

For each region, independently re-test the 3 Δβ contrasts using two methods that do NOT
use the module's permutation, then measure concordance with the module's perm p/sig columns:

  Contrast            module column (perm)        cross-check
  ------------------  --------------------------  ------------------------------------------
  somatic residual    SomaticResidualDbeta_*      quasi-binomial GLM cbind(k,n-k)~HP*sample,
                                                  test HP:sample interaction (= DSS interaction)
  germline ASM        GermlineAsmDbeta_*          GLM ~HP on normal reads (+ MWU on mean β)
  subclone HP1/HP2    SubcloneDbeta_HP{1,2}_*     GLM ~group on tumor germline-vs-carrier (+ MWU)

quasi-binomial GLM (family=Binomial, scale='X2') is the over-dispersion-aware standard that
DSS/methylKit implement as full beta-binomial; it is also exactly the over-dispersion correction
that per-CpG Fisher (#4 gap) lacks. Reads per region are the unit (each read = one molecule).

Inputs: per-region reads.tsv (hp,is_tumor) + methylation.csv (read x CpG β), matched by SNV pos
to significance_summary.csv. Output: concordance JSON (no report writing).
"""
import csv
import glob
import json
import math
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
import statsmodels.api as sm  # noqa: E402
from scipy.stats import mannwhitneyu  # noqa: E402

try:
    from statsmodels.othermod.betareg import BetaModel

    HAS_BETA = True
except Exception:
    HAS_BETA = False


def _squeeze(y):
    """Smithson-Verkuilen squeeze of proportions into the open interval (0,1) for beta regression."""
    y = np.asarray(y, float)
    nn = y.size
    return (y * (nn - 1) + 0.5) / nn


def beta_reg_p(yvals, xvals):
    """Beta regression mean β ~ x (single 0/1 covariate), read-level (N=#reads). Return coef p."""
    if not HAS_BETA:
        return None
    y = np.asarray(yvals, float)
    x = np.asarray(xvals, float)
    if y.size < 4 or len(set(x.tolist())) < 2:
        return None
    try:
        exog = sm.add_constant(x)
        res = BetaModel(_squeeze(y), exog).fit(disp=0)
        return float(res.pvalues[1])
    except Exception:
        return None


def beta_reg_interaction_p(yvals, hp, samp):
    """Beta regression mean β ~ hp + samp + hp:samp (read-level), interaction p (DSS analog, calibrated N=#reads)."""
    if not HAS_BETA:
        return None
    y = np.asarray(yvals, float)
    hp = np.asarray(hp, float)
    samp = np.asarray(samp, float)
    if y.size < 6 or len(set(hp.tolist())) < 2 or len(set(samp.tolist())) < 2:
        return None
    inter = hp * samp
    if len(set(inter.tolist())) < 2:
        return None
    try:
        exog = np.column_stack([np.ones_like(y), hp, samp, inter])
        res = BetaModel(_squeeze(y), exog).fit(disp=0)
        return float(res.pvalues[3])
    except Exception:
        return None

SUMMARY = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dbeta_chr1_s2/significance_summary.csv"
REGION_ROOT = sys.argv[2] if len(sys.argv) > 2 else "/tmp/dbeta_chr1_s2/filtered_snv_tp_chr1/chr1"
MAX_REGIONS = int(sys.argv[3]) if len(sys.argv) > 3 else 0  # 0 = all
SIG_A = 0.05


def fnum(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (ValueError, TypeError):
        return None


def hp_merged(hp):
    if hp in ("1", "HP1", "1-1"):
        return 0
    if hp in ("2", "HP2", "2-1"):
        return 1
    return -1


def hp_fine(hp):
    return {"1": 0, "HP1": 0, "1-1": 1, "2": 2, "HP2": 2, "2-1": 3}.get(hp, -1)


def read_region(rdir):
    """Return list of dict(read_id, hp, is_tumor, mean_beta, k, n) for one region inner dir."""
    rt = os.path.join(rdir, "reads", "reads.tsv")
    mc = os.path.join(rdir, "methylation", "methylation.csv")
    if not (os.path.exists(rt) and os.path.exists(mc)):
        return None
    meta = {}
    with open(rt) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            meta[row["read_id"]] = (row["hp"], row["is_tumor"] == "1")
    reads = []
    with open(mc) as f:
        header = f.readline()  # read_id,<pos>,<pos>,...
        for line in f:
            parts = line.rstrip("\n").split(",")
            rid = parts[0]
            if rid not in meta:
                continue
            betas = [float(v) for v in parts[1:] if v != "NA" and v != ""]
            if not betas:
                continue
            hp, is_t = meta[rid]
            arr = np.array(betas)
            reads.append(
                {
                    "hp": hp,
                    "is_tumor": is_t,
                    "mean_beta": float(arr.mean()),
                    "k": int((arr > 0.5).sum()),  # binarize at 0.5 for count-based beta-binomial/GLM
                    "n": int(arr.size),
                }
            )
    return reads


def glm_quasi_coef_p(k, n, x):
    """Quasi-binomial GLM cbind(k, n-k) ~ x (single 0/1 covariate). Return (coef_sign_effect, p)."""
    k = np.asarray(k, float)
    n = np.asarray(n, float)
    x = np.asarray(x, float)
    if len(set(x.tolist())) < 2 or (n <= 0).any():
        return None, None
    endog = np.column_stack([k, n - k])
    exog = sm.add_constant(x)
    try:
        m = sm.GLM(endog, exog, family=sm.families.Binomial())
        res = m.fit(scale="X2")  # quasi-binomial: over-dispersion via Pearson chi2
        return float(res.params[1]), float(res.pvalues[1])
    except Exception:
        return None, None


def glm_interaction_p(k, n, hp, samp):
    """Quasi-binomial GLM cbind(k,n-k) ~ hp + samp + hp:samp. Return interaction p (DSS analog)."""
    k = np.asarray(k, float)
    n = np.asarray(n, float)
    hp = np.asarray(hp, float)
    samp = np.asarray(samp, float)
    # need both HP and both sample present with the 4-cell design non-degenerate
    if len(set(hp.tolist())) < 2 or len(set(samp.tolist())) < 2:
        return None
    inter = hp * samp
    if len(set(inter.tolist())) < 2:
        return None
    exog = np.column_stack([np.ones_like(k), hp, samp, inter])
    try:
        res = sm.GLM(np.column_stack([k, n - k]), exog, family=sm.families.Binomial()).fit(scale="X2")
        return float(res.pvalues[3])
    except Exception:
        return None


def mwu_p(a, b):
    if len(a) < 1 or len(b) < 1:
        return None
    if len(set(a + b)) < 2:
        return None
    try:
        return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except Exception:
        return None


# ---- module perm results keyed by SNV pos ----
mod = {}
with open(SUMMARY) as f:
    for row in csv.DictReader(f):
        pos = row.get("Pos") or row.get("pos")
        if not pos:
            continue
        mod[pos] = {
            "som_p": fnum(row.get("SomaticResidualDbeta_P")),
            "som_sig": row.get("SomaticResidualDbeta_Sig") == "true",
            "germ_p": fnum(row.get("GermlineAsmDbeta_P")),
            "germ_sig": row.get("GermlineAsmDbeta_Sig") == "true",
            "s1_p": fnum(row.get("SubcloneDbeta_HP1_P")),
            "s1_sig": row.get("SubcloneDbeta_HP1_Sig") == "true",
            "s2_p": fnum(row.get("SubcloneDbeta_HP2_P")),
            "s2_sig": row.get("SubcloneDbeta_HP2_Sig") == "true",
        }

region_dirs = sorted(glob.glob(os.path.join(REGION_ROOT, "chr*")))
if MAX_REGIONS:
    region_dirs = region_dirs[:MAX_REGIONS]

# concordance accumulators per contrast: list of (perm_sig, cc_sig, perm_p, cc_p)
acc = {"som": [], "germ": [], "sub": []}
n_processed = 0

for rd in region_dirs:
    # pos = trailing number of chr1_<pos>
    base = os.path.basename(rd)
    pos = base.split("_")[-1]
    if pos not in mod:
        continue
    inner = glob.glob(os.path.join(rd, "chr*"))
    if not inner:
        continue
    reads = read_region(inner[0])
    if not reads:
        continue
    n_processed += 1

    # somatic residual: interaction HP(merged) x sample, all HP-labeled reads
    hp_l, samp_l, k_l, n_l, beta_l = [], [], [], [], []
    for r in reads:
        hm = hp_merged(r["hp"])
        if hm < 0:
            continue
        hp_l.append(hm)
        samp_l.append(1 if r["is_tumor"] else 0)
        k_l.append(r["k"])
        n_l.append(r["n"])
        beta_l.append(r["mean_beta"])
    p_int = glm_interaction_p(k_l, n_l, hp_l, samp_l) if hp_l else None
    p_breg = beta_reg_interaction_p(beta_l, hp_l, samp_l) if hp_l else None
    if p_int is not None and mod[pos]["som_p"] is not None:
        acc["som"].append((mod[pos]["som_sig"], p_int <= SIG_A, mod[pos]["som_p"], p_int, p_breg))

    # germline ASM: normal reads, HP merged
    nk, nn, nx, nbeta0, nbeta1 = [], [], [], [], []
    for r in reads:
        if r["is_tumor"]:
            continue
        hm = hp_merged(r["hp"])
        if hm < 0:
            continue
        nk.append(r["k"])
        nn.append(r["n"])
        nx.append(hm)
        (nbeta0 if hm == 0 else nbeta1).append(r["mean_beta"])
    _, gp = glm_quasi_coef_p(nk, nn, nx) if nx else (None, None)
    gmwu = mwu_p(nbeta0, nbeta1)
    gbreg = beta_reg_p(nbeta0 + nbeta1, [0] * len(nbeta0) + [1] * len(nbeta1))
    if gp is not None and mod[pos]["germ_p"] is not None:
        acc["germ"].append((mod[pos]["germ_sig"], gp <= SIG_A, mod[pos]["germ_p"], gp, gmwu, gbreg))

    # subclone HP1 & HP2: tumor germline vs carrier
    for fine_g, fine_c, key in [(0, 1, "s1"), (2, 3, "s2")]:
        sk, sn, sx, sb0, sb1 = [], [], [], [], []
        for r in reads:
            if not r["is_tumor"]:
                continue
            hf = hp_fine(r["hp"])
            if hf == fine_g:
                grp = 0
            elif hf == fine_c:
                grp = 1
            else:
                continue
            sk.append(r["k"])
            sn.append(r["n"])
            sx.append(grp)
            (sb0 if grp == 0 else sb1).append(r["mean_beta"])
        _, sp = glm_quasi_coef_p(sk, sn, sx) if sx else (None, None)
        smwu = mwu_p(sb0, sb1)
        sbreg = beta_reg_p(sb0 + sb1, [0] * len(sb0) + [1] * len(sb1))
        mp, ms = mod[pos][key + "_p"], mod[pos][key + "_sig"]
        if sp is not None and mp is not None:
            acc["sub"].append((ms, sp <= SIG_A, mp, sp, smwu, sbreg))


def _agree_rate(rows, perm_sig, idx):
    """Sig-agreement (perm vs method at p-index `idx`) on the subset where that p exists."""
    sub = [(ps, r[idx] <= SIG_A) for ps, r in zip(perm_sig, rows) if r[idx] is not None]
    if not sub:
        return None, None
    sig_rate = round(float(np.mean([s[1] for s in sub])), 4)
    agree = round(float(np.mean([s[0] == s[1] for s in sub])), 4)
    return sig_rate, agree


def concord(rows, mwu_idx=None, breg_idx=None):
    if not rows:
        return {"n": 0}
    perm_sig = np.array([r[0] for r in rows])
    cc_sig = np.array([r[1] for r in rows])
    perm_p = np.array([r[2] for r in rows])
    cc_p = np.array([r[3] for r in rows])
    n = len(rows)
    agree = int((perm_sig == cc_sig).sum())
    a = int((perm_sig & cc_sig).sum())
    b = int((perm_sig & ~cc_sig).sum())
    c = int((~perm_sig & cc_sig).sum())
    d = int((~perm_sig & ~cc_sig).sum())
    # Cohen's kappa
    po = (a + d) / n
    pe = ((a + b) * (a + c) + (c + d) * (b + d)) / (n * n)
    kappa = (po - pe) / (1 - pe) if (1 - pe) > 1e-9 else None
    # Spearman of -log10 p (robust to floor ties); compute manually
    lp1 = -np.log10(np.clip(perm_p, 1e-6, 1))
    lp2 = -np.log10(np.clip(cc_p, 1e-300, 1))
    rp1 = np.argsort(np.argsort(lp1))
    rp2 = np.argsort(np.argsort(lp2))
    sp_corr = float(np.corrcoef(rp1, rp2)[0, 1]) if n > 2 else None
    out = {
        "n": n,
        "perm_sig_rate": round(float(perm_sig.mean()), 4),
        "glm_quasibinom_sig_rate": round(float(cc_sig.mean()), 4),
        "glm_quasibinom_vs_perm_agreement": round(agree / n, 4),
        "glm_quasibinom_kappa": round(kappa, 4) if kappa is not None else None,
        "confusion_perm_x_glm": {"both_sig": a, "perm_only": b, "glm_only": c, "both_ns": d},
        "glm_quasibinom_spearman_neglog10p": round(sp_corr, 4) if sp_corr is not None else None,
    }
    if mwu_idx is not None:
        ms, ma = _agree_rate(rows, perm_sig, mwu_idx)
        out["mwu_sig_rate"] = ms
        out["mwu_vs_perm_agreement"] = ma
    if breg_idx is not None:
        bs, ba = _agree_rate(rows, perm_sig, breg_idx)
        out["betareg_sig_rate"] = bs
        out["betareg_vs_perm_agreement"] = ba
    return out


result = {
    "summary": SUMMARY,
    "region_root": REGION_ROOT,
    "n_regions_processed": n_processed,
    "has_betareg": HAS_BETA,
    "method": "3 read-level calibrated tests should agree (perm[module] / MWU[non-param] / beta-regression[beta-family, N=#reads]); count-based quasi-binomial GLM (binarize β>0.5, scale=X2) shown as naive per-CpG-count DSS-style which over-calls via within-read CpG correlation",
    "somatic_residual_interaction": concord(acc["som"], breg_idx=4),
    "germline_asm": concord(acc["germ"], mwu_idx=4, breg_idx=5),
    "subclone": concord(acc["sub"], mwu_idx=4, breg_idx=5),
}
print(json.dumps(result, indent=2, ensure_ascii=False))
