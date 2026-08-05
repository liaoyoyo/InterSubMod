#!/usr/bin/env python3
# DEMO / PILOT (subset, simulated — validates the STATISTICAL PRINCIPLE of the two V6 FIXes,
# NOT real BAM data). 2026-06-21.
# FIX-1: beta-binomial vs Fisher Type-I under read-level over-dispersion (epiallele correlation).
# FIX-2: correlation-preserving null vs naive null for over-clustering 1 clone.
# Honest caveat: per memory project_subcluster_cluster_count_determination, the simple
# correlation-preserving null OVER-REJECTED on real chr1 (81%); real data needs model-based.
import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar
from scipy.special import betaln, gammaln
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import json

rng = np.random.default_rng(42)
OUT = {}

# ---------- helpers ----------
def bb_logpmf(k, n, mu, rho):
    if rho <= 1e-9:
        return stats.binom.logpmf(k, n, mu)
    a = mu*(1-rho)/rho; b = (1-mu)*(1-rho)/rho
    return (gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1) + betaln(k+a, n-k+b) - betaln(a, b))

def bb_rvs(n, mu, rho, rng):
    if rho <= 1e-9:
        return rng.binomial(n, mu)
    a = mu*(1-rho)/rho; b = (1-mu)*(1-rho)/rho
    return rng.binomial(n, rng.beta(a, b))

# ================= TEST 1: beta-binomial vs Fisher (Type-I under over-dispersion) =================
def test1(N=3000, nA=30, nB=30, mu=0.5, rho_true=0.15):
    kA = np.array([bb_rvs(nA, mu, rho_true, rng) for _ in range(N)])
    kB = np.array([bb_rvs(nB, mu, rho_true, rng) for _ in range(N)])
    # Fisher per locus (assumes binomial / independent reads)
    fp_fisher = np.mean([stats.fisher_exact([[kA[i], nA-kA[i]], [kB[i], nB-kB[i]]])[1] < 0.05
                         for i in range(N)])
    # estimate shared dispersion rho across loci under H0 (shared mu per locus)
    def nll(logr):
        rho = 1/(1+np.exp(-logr)); ll = 0.0
        for i in range(N):
            mui = min(max((kA[i]+kB[i])/(nA+nB), 1e-4), 1-1e-4)
            ll += bb_logpmf(kA[i], nA, mui, rho) + bb_logpmf(kB[i], nB, mui, rho)
        return -ll
    res = minimize_scalar(nll, bounds=(-8, 4), method='bounded')
    rho_hat = 1/(1+np.exp(-res.x))
    # per-locus beta-binomial LRT (mu_A=mu_B vs differ, shared rho_hat), df=1
    def bb_p(k1, n1, k2, n2):
        mu0 = min(max((k1+k2)/(n1+n2), 1e-4), 1-1e-4)
        ll0 = bb_logpmf(k1, n1, mu0, rho_hat) + bb_logpmf(k2, n2, mu0, rho_hat)
        m1 = min(max(k1/n1, 1e-4), 1-1e-4); m2 = min(max(k2/n2, 1e-4), 1-1e-4)
        ll1 = bb_logpmf(k1, n1, m1, rho_hat) + bb_logpmf(k2, n2, m2, rho_hat)
        return stats.chi2.sf(max(2*(ll1-ll0), 0), 1)
    fp_bb = np.mean([bb_p(kA[i], nA, kB[i], nB) < 0.05 for i in range(N)])
    # POWER check under a TRUE difference (mu_A=0.5 vs mu_B=0.68)
    muB2 = 0.68
    kA2 = np.array([bb_rvs(nA, mu,  rho_true, rng) for _ in range(N)])
    kB2 = np.array([bb_rvs(nB, muB2, rho_true, rng) for _ in range(N)])
    pow_fisher = np.mean([stats.fisher_exact([[kA2[i], nA-kA2[i]], [kB2[i], nB-kB2[i]]])[1] < 0.05
                          for i in range(N)])
    pow_bb = np.mean([bb_p(kA2[i], nA, kB2[i], nB) < 0.05 for i in range(N)])
    return dict(N=N, rho_true=rho_true, rho_hat=round(float(rho_hat), 4),
                typeI_fisher=round(float(fp_fisher), 4), typeI_betabinom=round(float(fp_bb), 4),
                power_fisher=round(float(pow_fisher), 4), power_betabinom=round(float(pow_bb), 4))

# ================= TEST 2: correlation-preserving vs naive null (over-clustering 1 clone) =================
def gen_one_clone(n_reads=50, n_cpg=30, sigma_read=1.2, rng=rng):
    # ONE clone, but reads carry a latent global offset z_i (epiallele correlation)
    base = rng.normal(0, 0.8, n_cpg)            # per-CpG baseline logit
    z = rng.normal(0, sigma_read, n_reads)      # per-READ offset -> within-clone correlation
    logit = base[None, :] + z[:, None]
    p = 1/(1+np.exp(-logit))
    return rng.binomial(1, p).astype(float)     # n_reads x n_cpg binary methylation

def split2_stat(M):
    # pseudo-F-like: 2-cluster hierarchical split, between/within distance ratio
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        D[i] = np.mean(np.abs(M - M[i]), axis=1)  # mean L1 (Hamming-rate) distance
    Z = linkage(squareform(D, checks=False), method='average')
    lab = fcluster(Z, 2, criterion='maxclust')
    if len(np.unique(lab)) < 2:
        return 0.0
    g1, g2 = M[lab == 1], M[lab == 2]
    within = (np.sum([np.mean(np.abs(g1 - r)) for r in g1]) + np.sum([np.mean(np.abs(g2 - r)) for r in g2]))
    tot = np.sum([np.mean(np.abs(M - r)) for r in M])
    between = tot - within
    return between / (within + 1e-9)

def test2(M_DATASETS=150, B=149, n_reads=50, n_cpg=30, sigma_read=1.2):
    naive_hits = 0; corr_hits = 0
    for _ in range(M_DATASETS):
        X = gen_one_clone(n_reads, n_cpg, sigma_read, rng)   # ONE clone (null truth: k=1)
        obs = split2_stat(X)
        # naive null: permute each CpG column independently (destroys read-level correlation)
        naive_null = []
        for _ in range(B):
            Xp = X.copy()
            for j in range(n_cpg):
                Xp[:, j] = rng.permutation(Xp[:, j])
            naive_null.append(split2_stat(Xp))
        p_naive = (1 + np.sum(np.array(naive_null) >= obs)) / (B + 1)
        # correlation-preserving null: regenerate 1-clone-with-offset (parametric bootstrap)
        corr_null = [split2_stat(gen_one_clone(n_reads, n_cpg, sigma_read, rng)) for _ in range(B)]
        p_corr = (1 + np.sum(np.array(corr_null) >= obs)) / (B + 1)
        naive_hits += (p_naive < 0.05); corr_hits += (p_corr < 0.05)
    return dict(M_datasets=M_DATASETS, B=B, truth="k=1 (one clone, epiallele-correlated reads)",
                false_kge2_naive_null=round(naive_hits / M_DATASETS, 4),
                false_kge2_corr_preserving_null=round(corr_hits / M_DATASETS, 4))

OUT["FIX1_betabinom_vs_fisher"] = test1()
OUT["FIX2_corrnull_vs_naive"] = test2()
OUT["caveat"] = ("Simulated principle-validation only. Per memory project_subcluster_cluster_count_determination: "
                 "the SIMPLE correlation-preserving null OVER-REJECTED on real chr1 (81%; +AR1 -> 86%) because real "
                 "ONT methylation is discrete + doubly-correlated -> real data needs full model-based (Path B). "
                 "beta-binomial (FIX-1) is the cleaner real-data win.")
print(json.dumps(OUT, indent=2, ensure_ascii=False))
with open("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/20260621_v6fix_pilot/result.json", "w") as f:
    json.dump(OUT, f, indent=2, ensure_ascii=False)
