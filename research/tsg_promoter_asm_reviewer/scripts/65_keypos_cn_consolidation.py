#!/usr/bin/env python3
"""
65 - Consolidate the 38 key positions across 6 samples: existence (coverage) +
copy-number (CN) consistency.

QUESTION (user 2026-06-03): for these key positions, (a) do they all CLEARLY EXIST
(are they covered/detectable in every sample?) and (b) is there a CONSISTENT
phenomenon/problem with the COPY-NUMBER situation across samples?

DATA (read-back, §13 layer-A):
  - 6x <sample>_keypos.json (script 57): per (position, sample) n_reads_window,
    hp_groups (per-HP-tag read counts), somatic_in_sample, has_somatic_subhap.
  - master_o1_cn.tsv (cn_confound): HCC1395 SEQC2 CN ground-truth per position
    (median_cn, cn_class, cnloh_flag, loh_status). HCC1395-ONLY (SEQC2).

CN handling:
  - HCC1395: actual SEQC2 median_cn + cn_class.
  - ALL samples: relative coverage (n_reads / sample-median over the 38 positions) as a
    crude CN/depth proxy; + local zygosity from HP structure (single-HP=LOH-like,
    both-HP=het, subhap=somatic). NO orthogonal CN for non-HCC1395 (disclosed).

ANALYSES:
  1. EXISTENCE: coverage matrix; are all 38 positions covered >=COV_MIN in all samples?
  2. HCC1395 CN VALIDATION: does HCC1395 relative coverage track SEQC2 median_cn
     (Spearman)? does HP-structure (LOH) match cn_class (cnLOH/loss)?
  3. CROSS-SAMPLE CN CONSISTENCY: per position, is the zygosity (LOH vs het) the SAME
     across samples (consistent CN) or cancer-specific (mixed)? fraction-LOH per position.

Output:
  genome_survey_v2/cn_confound/cross_sample/keypos_cn_consolidation.json
"""
import os, json
import numpy as np
from collections import Counter
from scipy.stats import spearmanr

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
MASTER_CN = f"{ROOT}/genome_survey_v2/cn_confound/master_o1_cn.tsv"
OUT = f"{CS}/keypos_cn_consolidation.json"

ALL = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
CANCER = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
          "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}
COV_MIN = 10
LOH_DOM = 10     # dominant hap >= 10 reads
LOH_MINOR = 3    # minor hap < 3 reads => LOH-like


def zygosity(hp):
    """from hp_groups dict -> ('LOH'|'het'|'subhap'|'low') + (n1,n2,nsub)."""
    n1 = hp.get("1", {}).get("n", 0)
    n2 = hp.get("2", {}).get("n", 0)
    nsub = sum(v.get("n", 0) for k, v in hp.items() if "-" in k)
    total = sum(v.get("n", 0) for v in hp.values())
    if total < COV_MIN:
        return "low", n1, n2, nsub
    if nsub >= LOH_DOM:
        return "subhap", n1, n2, nsub
    hi, lo = max(n1, n2), min(n1, n2)
    if hi >= LOH_DOM and lo < LOH_MINOR:
        return "LOH", n1, n2, nsub
    if n1 >= LOH_DOM and n2 >= LOH_DOM:
        return "het", n1, n2, nsub
    return "het", n1, n2, nsub   # default partial -> het-ish


# ---- load HCC1395 SEQC2 CN ----
cn395 = {}
with open(MASTER_CN) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(hdr)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(ci.values()):
            continue
        key = f"{p[ci['chrom']]}:{p[ci['somatic_pos']]}"
        try:
            mcn = float(p[ci["median_cn"]]) if p[ci["median_cn"]] not in ("", "nan") else None
        except ValueError:
            mcn = None
        cn395[key] = dict(median_cn=mcn, cn_class=p[ci["cn_class"]],
                          cnloh=p[ci["cnloh_flag"]], loh_status=p[ci["loh_status"]])

# ---- load 6 keypos JSONs ----
data = {s: json.load(open(f"{CS}/{s}_keypos.json")) for s in ALL}
positions = [(r["chrom"], r["pos"], r["gene"]) for r in data["HCC1395"]["per_position"]]
idx = {s: {f"{r['chrom']}:{r['pos']}": r for r in data[s]["per_position"]} for s in ALL}

# sample-median coverage over the 38 positions (for relative-CN proxy)
sample_med_cov = {}
for s in ALL:
    covs = [r["n_reads_window"] for r in data[s]["per_position"]]
    sample_med_cov[s] = float(np.median(covs)) if covs else 1.0

# ---- per-position consolidation ----
per_pos = []
for chrom, pos, gene in positions:
    key = f"{chrom}:{pos}"
    cells = {}
    for s in ALL:
        r = idx[s].get(key, {})
        cov = r.get("n_reads_window", 0)
        zyg, n1, n2, nsub = zygosity(r.get("hp_groups", {}))
        cells[s] = dict(cancer=CANCER[s], cov=cov,
                        rel_cov=round(cov / sample_med_cov[s], 2) if sample_med_cov[s] else None,
                        zygosity=zyg, n_hp1=n1, n_hp2=n2, n_subhap=nsub,
                        covered=cov >= COV_MIN)
    # cross-sample consistency of zygosity
    zygs = [cells[s]["zygosity"] for s in ALL]
    n_loh = sum(1 for z in zygs if z == "LOH")
    n_het = sum(1 for z in zygs if z == "het")
    n_sub = sum(1 for z in zygs if z == "subhap")
    n_covered = sum(1 for s in ALL if cells[s]["covered"])
    zyg_mode, zyg_modecount = Counter(zygs).most_common(1)[0]
    per_pos.append(dict(
        pos=key, gene=gene,
        hcc1395_cn=cn395.get(key, {}),
        n_covered=n_covered, all_covered=(n_covered == len(ALL)),
        min_cov=min(cells[s]["cov"] for s in ALL),
        zyg_counts=dict(LOH=n_loh, het=n_het, subhap=n_sub,
                        low=sum(1 for z in zygs if z == "low")),
        zyg_consistent=(zyg_modecount == len(ALL)),
        zyg_mode=zyg_mode, zyg_modecount=zyg_modecount,
        per_sample=cells,
    ))

# ---- analysis 1: existence ----
n_all_covered = sum(1 for p in per_pos if p["all_covered"])
min_cov_overall = min(p["min_cov"] for p in per_pos)

# ---- analysis 2: HCC1395 CN validation ----
xcn, yrelcov = [], []
loh_match = 0
loh_total = 0
for p in per_pos:
    cn = p["hcc1395_cn"]
    mcn = cn.get("median_cn")
    relcov = p["per_sample"]["HCC1395"]["rel_cov"]
    if mcn is not None and relcov is not None:
        xcn.append(mcn); yrelcov.append(relcov)
    # LOH structure vs SEQC2 loh_status
    if cn.get("loh_status"):
        loh_total += 1
        is_loh_struct = p["per_sample"]["HCC1395"]["zygosity"] in ("LOH", "subhap")
        is_loh_seqc2 = cn["loh_status"] == "LOH"
        if is_loh_struct == is_loh_seqc2:
            loh_match += 1
cn_cov_rho = (round(float(spearmanr(xcn, yrelcov)[0]), 3) if len(xcn) >= 5 else None)
cn_cov_p = (round(float(spearmanr(xcn, yrelcov)[1]), 4) if len(xcn) >= 5 else None)

# ---- analysis 3: cross-sample CN/zygosity consistency ----
n_zyg_consistent = sum(1 for p in per_pos if p["zyg_consistent"])
# per-position fraction LOH across samples
frac_loh_dist = [p["zyg_counts"]["LOH"] / len(ALL) for p in per_pos]

out = dict(
    meta=dict(script="65_keypos_cn_consolidation.py", samples=ALL, cancer_types=CANCER,
              n_positions=len(per_pos), cov_min=COV_MIN,
              sample_median_coverage=sample_med_cov,
              note=("Existence = coverage>=10. CN: HCC1395 uses SEQC2 median_cn/cn_class "
                    "(ground-truth); other samples have NO orthogonal CN -> rel_cov "
                    "(n_reads/sample-median) is a crude depth/CN proxy + HP-structure "
                    "zygosity (LOH/het/subhap). Cross-sample 'CN consistency' = whether the "
                    "same position has the same zygosity across cancers.")),
    existence=dict(n_positions=len(per_pos), n_all_covered=n_all_covered,
                   frac_all_covered=round(n_all_covered / len(per_pos), 3),
                   min_cov_overall=min_cov_overall),
    hcc1395_cn_validation=dict(
        n_cn_pairs=len(xcn), cn_relcov_spearman=cn_cov_rho, cn_relcov_p=cn_cov_p,
        loh_structure_vs_seqc2_match=loh_match, loh_total=loh_total,
        loh_match_frac=(round(loh_match / loh_total, 3) if loh_total else None)),
    cross_sample_cn_consistency=dict(
        n_zyg_consistent=n_zyg_consistent,
        frac_zyg_consistent=round(n_zyg_consistent / len(per_pos), 3),
        mean_frac_loh_per_pos=round(float(np.mean(frac_loh_dist)), 3),
        interpretation=("zygosity_consistent positions have the SAME LOH/het structure "
                        "across all 6 cancers; low consistency => CN/zygosity is "
                        "cancer-specific (each cancer its own CN profile).")),
    per_position=per_pos,
)
with open(OUT, "w") as f:
    json.dump(out, f, indent=2,
              default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
print(f"[65] wrote {OUT}")
print(f"[65] EXISTENCE: {n_all_covered}/{len(per_pos)} positions covered>=10 in ALL 6 samples "
      f"(min_cov_overall={min_cov_overall})")
print(f"[65] HCC1395 CN validation: rel_cov~median_cn Spearman={cn_cov_rho} (p={cn_cov_p}, n={len(xcn)}); "
      f"LOH-structure vs SEQC2 match={loh_match}/{loh_total}")
print(f"[65] CROSS-SAMPLE zygosity consistency: {n_zyg_consistent}/{len(per_pos)} positions "
      f"same zygosity across all 6 (mean frac-LOH per pos={np.mean(frac_loh_dist):.3f})")
print(f"[65] sample median coverage: " + ", ".join(f"{s}={sample_med_cov[s]:.0f}" for s in ALL))
