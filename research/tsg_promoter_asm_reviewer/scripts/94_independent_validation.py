#!/usr/bin/env python3
"""
94 - TP/FP-INDEPENDENT validation that a locus has a REAL methylation difference
that is genuinely associated with the (SNV-derived) HP tag — no SEQC2 TP/FP used.

Non-circularity: HP labels come from LongPhase-S germline-SNP phasing
(germline_phased_merged.vcf.gz), NOT from methylation. We test on BASE HP
(HP1-family vs HP2-family); HP-fine sub-haplotypes (which ISM derives partly from
methylation) are NOT used as the label, to avoid circularity.

Per curated locus (read x CpG = methylation.csv; base HP + strand = reads.tsv):
  L1 significance  : HP-permutation null on |Δβ| (N perms) + Mann-Whitney U on
                     per-read mean methylation (HP1 vs HP2).  -> is it beyond chance?
  L2 replication   : CpG split-half (even vs odd CpGs) sign-concordant
                     + strand consistency (fwd-only vs rev-only Δβ same sign).
                     -> not driven by one CpG / one strand artifact.
  L3 stability     : read bootstrap, fraction of resamples with same Δβ sign.
  validated = perm_p<=0.05 & MWU_p<=0.05 & splithalf_concordant
              & strand_concordant(where measurable) & bootstrap_stable

Then: BH-FDR across loci; pass-rate per Tier (A/B/none) and for Δβ-only-no-HP-struct.
Output: display_v2/validation.json
"""
import os, csv, json, glob
os.environ.setdefault("TMPDIR", "/big7_disk/liaoyoyo2001/tmp")
from multiprocessing import Pool
import numpy as np
from scipy.stats import mannwhitneyu

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
DV = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample/display_v2"
SCAN = "/big7_disk/liaoyoyo2001/ism_display_scan"
N_PERM, N_BOOT = 2000, 1000
MIN_PER_GROUP = 4          # min reads per HP family
MIN_HALF = 0.03            # min |Δβ| in each split-half to count as concordant


def region_dir(cls, chrom, pos):
    for d in glob.glob(f"{SCAN}/HCC1395_{cls}/curated_{cls}/{chrom}/{chrom}_{pos}/*"):
        if os.path.exists(f"{d}/methylation/methylation.csv"):
            return d
    return None


def load(rd):
    rows = list(csv.reader(open(f"{rd}/methylation/methylation.csv")))
    rids = [int(r[0]) for r in rows[1:]]
    M = np.array([[np.nan if x == "NA" else float(x) for x in r[1:]] for r in rows[1:]])
    base, strand = {}, {}
    for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t"):
        h = r.get("hp", "")
        b = 1 if h in ("1", "1-1") else (2 if h in ("2", "2-1") else 0)
        base[int(r["read_id"])] = b
        strand[int(r["read_id"])] = r.get("strand", "")
    return rids, M, base, strand


def dbeta(means, lab):
    a = means[lab == 1]; b = means[lab == 2]
    if len(a) < 1 or len(b) < 1:
        return np.nan
    return np.nanmean(a) - np.nanmean(b)


def work(rec):
    cls, chrom, pos = rec["cls"], rec["chr"], rec["pos"]
    key = f"{chrom}_{pos}"
    rd = region_dir(cls, chrom, pos)
    if rd is None:
        return key, None
    try:
        rids, M, base, strand = load(rd)
        lab = np.array([base.get(i, 0) for i in rids])
        str_ = np.array([strand.get(i, "") for i in rids])
        keep = lab > 0
        M, lab, str_ = M[keep], lab[keep], str_[keep]
        n1, n2 = int((lab == 1).sum()), int((lab == 2).sum())
        if n1 < MIN_PER_GROUP or n2 < MIN_PER_GROUP:
            return key, dict(ok=False, reason="too few per HP family", n1=n1, n2=n2)
        rmean = np.nanmean(M, axis=1)
        obs = dbeta(rmean, lab)
        if not np.isfinite(obs):
            return key, dict(ok=False, reason="nan dbeta")
        # L1 permutation null on |Δβ|
        rng = np.random.RandomState(abs(hash(key)) % (2**31))
        cnt = 0
        for _ in range(N_PERM):
            if abs(dbeta(rmean, rng.permutation(lab))) >= abs(obs):
                cnt += 1
        perm_p = (cnt + 1) / (N_PERM + 1)
        # Mann-Whitney on per-read means
        a, b = rmean[lab == 1], rmean[lab == 2]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        try:
            mwu_p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        except ValueError:
            mwu_p = 1.0
        # L2a split-half by CpG index
        ncol = M.shape[1]
        ev, od = np.arange(ncol) % 2 == 0, np.arange(ncol) % 2 == 1
        db_ev = dbeta(np.nanmean(M[:, ev], axis=1), lab) if ev.sum() else np.nan
        db_od = dbeta(np.nanmean(M[:, od], axis=1), lab) if od.sum() else np.nan
        sh_conc = bool(np.isfinite(db_ev) and np.isfinite(db_od)
                       and np.sign(db_ev) == np.sign(db_od)
                       and abs(db_ev) >= MIN_HALF and abs(db_od) >= MIN_HALF)
        # L2b strand consistency
        def strand_db(s):
            m = str_ == s
            if (lab[m] == 1).sum() >= 3 and (lab[m] == 2).sum() >= 3:
                return dbeta(rmean[m], lab[m])
            return np.nan
        df, dr = strand_db("+"), strand_db("-")
        strand_measurable = np.isfinite(df) and np.isfinite(dr)
        strand_conc = bool(strand_measurable and np.sign(df) == np.sign(dr))
        # L3 bootstrap sign stability
        idx = np.arange(len(rmean)); sb = 0
        for _ in range(N_BOOT):
            s = rng.choice(idx, len(idx), replace=True)
            d = dbeta(rmean[s], lab[s])
            if np.isfinite(d) and np.sign(d) == np.sign(obs):
                sb += 1
        boot_stab = sb / N_BOOT
        return key, dict(ok=True, n1=n1, n2=n2, obs_dbeta=round(float(obs), 4),
                         perm_p=round(perm_p, 4), mwu_p=round(float(mwu_p), 4),
                         splithalf_concordant=sh_conc, db_even=round(float(db_ev), 3), db_odd=round(float(db_od), 3),
                         strand_measurable=bool(strand_measurable), strand_concordant=strand_conc,
                         db_fwd=None if not np.isfinite(df) else round(float(df), 3),
                         db_rev=None if not np.isfinite(dr) else round(float(dr), 3),
                         boot_stability=round(boot_stab, 3))
    except Exception as e:
        return key, dict(ok=False, reason=f"{type(e).__name__}: {e}")


def bh_fdr(pvals):
    n = len(pvals); order = np.argsort(pvals)
    q = np.empty(n); prev = 1.0
    for rank, i in enumerate(reversed(order)):
        k = n - rank
        prev = min(prev, pvals[i] * n / k)
        q[i] = prev
    return q


def main():
    curated = json.load(open(f"{DV}/curated_loci.json"))
    tier = json.load(open(f"{DV}/tier_assignment.json"))
    man = {m["chr"] + "_" + str(m["pos"]): m for m in json.load(open(f"{DV}/manifest.json"))}
    print(f"[94] validating {len(curated)} loci (perm={N_PERM} boot={N_BOOT}) ...")
    res = {}
    with Pool(12) as p:
        for i, (key, r) in enumerate(p.imap_unordered(work, curated, chunksize=8)):
            if r:
                res[key] = r
            if (i + 1) % 300 == 0:
                print(f"   ...{i+1}/{len(curated)}")

    # BH-FDR on perm_p for the testable loci
    ok = {k: v for k, v in res.items() if v.get("ok")}
    keys = list(ok); pv = np.array([ok[k]["perm_p"] for k in keys])
    q = bh_fdr(pv)
    for k, qq in zip(keys, q):
        ok[k]["perm_q"] = round(float(qq), 4)
        v = ok[k]
        v["validated"] = bool(v["perm_q"] <= 0.05 and v["mwu_p"] <= 0.05
                              and v["splithalf_concordant"]
                              and (v["strand_concordant"] or not v["strand_measurable"])
                              and v["boot_stability"] >= 0.90)

    # aggregate by tier + by Δβ-only-no-structure
    def tier_of(k):
        return tier.get(k, "")
    agg = {}
    for t in ("A", "B", ""):
        ks = [k for k in ok if tier_of(k) == t]
        agg[t or "none"] = dict(n=len(ks), validated=sum(1 for k in ks if ok[k]["validated"]),
                                strand_meas=sum(1 for k in ks if ok[k]["strand_measurable"]),
                                strand_conc=sum(1 for k in ks if ok[k]["strand_concordant"]))
    # Δβ-only loci: |Δβ|>=0.2 but NOT HP-structured (latent=False & tier not A/B)
    dbonly = [k for k in ok if abs(man.get(k, {}).get("db", 0)) >= 0.2 and tier_of(k) not in ("A", "B")]
    dbonly_val = sum(1 for k in dbonly if ok[k]["validated"])

    summary = dict(n_total=len(curated), n_testable=len(ok),
                   tier_agg=agg, dbeta_only_n=len(dbonly), dbeta_only_validated=dbonly_val,
                   params=dict(N_PERM=N_PERM, N_BOOT=N_BOOT, MIN_PER_GROUP=MIN_PER_GROUP,
                               criteria="perm_q<=0.05 & MWU<=0.05 & splithalf & strand(if measurable) & boot>=0.90"))
    out = dict(summary=summary, loci=ok)
    with open(f"{DV}/validation.json", "w") as f:
        json.dump(out, f)

    print(f"\n[94] testable={len(ok)}/{len(curated)}")
    for t in ("A", "B", "none"):
        a = agg[t]
        print(f"   Tier {t:4s}: n={a['n']:4d}  validated={a['validated']:4d} "
              f"({100*a['validated']/max(a['n'],1):.0f}%)  strand-concordant={a['strand_conc']}/{a['strand_meas']}")
    print(f"   Δβ-only (|Δβ|≥0.2, no HP structure): n={len(dbonly)} validated={dbonly_val} "
          f"({100*dbonly_val/max(len(dbonly),1):.0f}%)  <- should be LOW if strand/splithalf catch artifacts")
    print(f"[94] wrote {DV}/validation.json")


if __name__ == "__main__":
    main()
