#!/usr/bin/env python3
"""Audit B' (clean over-dispersion test on REAL reads): per-CpG Fisher exact (ISM-style)
   under a read-level HP-shuffle null that PRESERVES each read's full CpG profile
   (= preserves within-read/clonal correlation). If the null mean 'fraction of CpGs
   with Fisher p<0.05' is >> nominal 5%, Fisher is anti-conservative at the AGGREGATE
   ASM-existence-rate level (audit FISHER-1). Tests it on real data, not simulation."""
import gzip, csv, glob, json
from collections import defaultdict
import numpy as np
from scipy.stats import fisher_exact
rng = np.random.default_rng(20260610)
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod"
B = f"{ROOT}/docs/method_comparison/20260609_ism_vs_external_methylation_tools/benchmark"
THR = 0.5; NPERM = 500; MIN_N = 5
LOCI = {
    "BRCA2":   glob.glob(f"{ROOT}/research/tsg_promoter_asm_reviewer/pipeline/cache/level1/chr13_32315128_*.tsv.gz")[0],
    "TBC1D16": glob.glob(f"{ROOT}/research/tsg_promoter_asm_reviewer/pipeline/cache/level1/chr17_79991120_*.tsv.gz")[0],
}
res = {}
for locus, level1 in LOCI.items():
    reads = {}
    with gzip.open(level1, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["bam_source"] != "tumor" or row["mod_code"] != "m":
                continue
            hp = row["haplotype_tag"]
            if hp not in ("1", "1-1"):
                continue
            rid = row["read_id"]
            reads.setdefault(rid, {"hp": hp, "cpg": {}})
            reads[rid]["cpg"][int(row["methyl_pos"])] = 1 if float(row["meth_call"]) >= THR else 0
    rids = list(reads)
    labels = np.array([reads[r]["hp"] for r in rids])

    def frac_sig(lab):
        bycpg = defaultdict(lambda: {"1": [0, 0], "1-1": [0, 0]})  # [meth, unmeth]
        for r, l in zip(rids, lab):
            for pos, b in reads[r]["cpg"].items():
                bycpg[pos][l][0 if b == 1 else 1] += 1
        ps = []
        for pos, g in bycpg.items():
            n1 = g["1"][0] + g["1"][1]; n2 = g["1-1"][0] + g["1-1"][1]
            if n1 >= MIN_N and n2 >= MIN_N:
                _, p = fisher_exact([[g["1"][0], g["1"][1]], [g["1-1"][0], g["1-1"][1]]])
                ps.append(p)
        if not ps:
            return None, 0
        return float(np.mean(np.array(ps) < 0.05)), len(ps)

    obs_frac, n_cpg = frac_sig(labels)
    null = []
    for _ in range(NPERM):
        sh = labels.copy(); rng.shuffle(sh)
        f, _ = frac_sig(sh)
        if f is not None:
            null.append(f)
    null = np.array(null)
    res[locus] = {
        "n_reads_HP1+HP1-1": len(rids), "n_cpg_tested": n_cpg,
        "observed_frac_fisher_sig": round(obs_frac, 4),
        "null_mean_frac_sig": round(float(null.mean()), 4),
        "null_median": round(float(np.median(null)), 4),
        "null_95pct": round(float(np.percentile(null, 95)), 4),
        "nominal_expected_if_valid": 0.05,
        "inflation_ratio_null/nominal": round(float(null.mean() / 0.05), 2),
        "obs_perm_p_vs_null": round(float((np.sum(null >= obs_frac) + 1) / (len(null) + 1)), 4),
        "verdict": ("Fisher ANTI-CONSERVATIVE at aggregate (null mean >> 5%) -> ASM-existence-rate inflated"
                    if null.mean() > 0.08 else
                    "Fisher aggregate FP ~nominal (null mean ~5%) -> per-CpG Fisher OK"),
    }
json.dump(res, open(f"{B}/runs/audit_B2.json", "w"), indent=2)
print(json.dumps(res, indent=2))
