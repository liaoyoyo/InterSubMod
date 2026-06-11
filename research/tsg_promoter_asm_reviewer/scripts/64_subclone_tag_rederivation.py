#!/usr/bin/env python3
"""
64_subclone_tag_rederivation.py — corrects the copy-partition framing after the
HP1-1 tag-definition audit (longphase-S HaplotagStrategy.cpp:505-516: HP1-1 =
germline-H1 + somatic-ALT = a SOMATIC SUBCLONE tag, not a copy).

Part A (validity of d_within as subclone-control): for BRCA2, confirm HP1-1/ref reads
are genuine subclone reads (REF at focal, carry somatic ALT elsewhere, methyl clusters
with subclone). [Numbers hard-recorded from the BAM/VCF check; see ledger 20260607_hp11_tag_definition_audit.]

Part B (clean re-derivation, all key loci): four measures per locus
  d_HP        = HP1 vs HP1-1                      (tag-based; subclone-ALT vs germline-REF)
  d_focal_CLEAN = germline-H1 {1,1-1} split by ACTUAL focal allele (ref vs alt)  <- bypasses subclone tag, computable everywhere
  d_within    = HP1-1: alt vs ref                (subclone-internal focal allele; VALID subclone control)
  d_copy      = HP1/ref vs HP1-1/ref             (germline vs subclone, same REF; = the 'subclone/copy' component)
Relabel: 'copy-partition' -> 'subclone/copy partition'. Quantification (d_within/d_HP) retained (d_within validated).
"""
import sys, glob, json
import numpy as np
sys.path.insert(0, "pipeline")
from lib.cis_asm_core import load_level1_plus, collapse_modtype, _beta_cpg, MIN_N
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"

# 10 Bonferroni-clean T3 + neutral note; group from the decision table
LOCI = [
    ("chr17", "79991120", "CIS-CAND clean (TBC1D16)"),
    ("chr13", "32315128", "BRCA2/ZAR1L"),
    ("chr5", "6201328", "copy-art"),
    ("chr19", "14434617", "mechanical"),
    ("chr4", "133868344", "untestable-strong"),
    ("chr7", "79941963", "untestable"),
    ("chr20", "59439285", "untestable"),
    ("chr20", "61415690", "untestable"),
    ("chr20", "61564264", "untestable"),
    ("chr11", "64557316", "untestable+highNM"),
    ("chr18", "11741161", "CpG-creation"),
]


def paired(A, B, n=MIN_N):
    sh = [c for c in set(A) & set(B) if A[c][1] >= n and B[c][1] >= n]
    return (round(float(np.mean([B[c][0] - A[c][0] for c in sh])), 3), len(sh)) if len(sh) >= 4 else (None, len(sh))


rows = []
for c, p, tag in LOCI:
    cache = glob.glob(f"{ROOT}/pipeline/cache/level1/{c}_{p}_*.tsv.gz")
    if not cache:
        rows.append({"locus": f"{c}:{p}", "tag": tag, "error": "no cache"}); continue
    Dp = load_level1_plus(cache[0]); D = collapse_modtype(Dp)
    b = lambda hps, al: _beta_cpg(D, p, "tumor", set(hps), set(al) if al else None)
    d_hp, n_hp = paired(b({"1"}, None), b({"1-1"}, None))
    d_clean, n_clean = paired(b({"1", "1-1"}, {"ref"}), b({"1", "1-1"}, {"alt"}))
    d_within, n_w = paired(b({"1-1"}, {"ref"}), b({"1-1"}, {"alt"}))
    d_copy, n_c = paired(b({"1"}, {"ref"}), b({"1-1"}, {"ref"}))
    # subclone fraction = how much of |d_HP| is the subclone/copy component (1 - |d_within|/|d_HP|)
    subclone_frac = (round(1 - abs(d_within) / abs(d_hp), 2)
                     if d_within is not None and d_hp not in (None, 0) else None)
    rows.append({"locus": f"{c}:{p}", "tag": tag,
                 "d_HP": d_hp, "d_focal_CLEAN": d_clean, "d_within": d_within, "d_copy": d_copy,
                 "n_clean": n_clean, "n_within": n_w,
                 "within_testable": d_within is not None,
                 "subclone_copy_frac": subclone_frac})

# Part A validity (recorded from the BAM/VCF audit this session)
validity = {
    "locus": "chr13:32315128 (BRCA2)",
    "HP11_ref_reads": 19, "HP11_alt_reads": 19,
    "focal_base_all_REF": "G x20 (REF), median_BQ=26.5 (confident REF, not ALT mis-call)",
    "carry_other_somatic_ALT": "19/19 (avg 1.53 other somatic ALT at chr13:32317522/32324831/32339132)",
    "methyl_clusters_with": "subclone (d_copy=-0.110 large vs HP1/ref; d_within=-0.023 small vs HP1-1/alt)",
    "verdict": "HP1-1/ref = genuine somatic-subclone reads, REF at focal -> d_within is a VALID subclone-controlled focal-allele estimate",
    "source": "longphase-s/src/haplotag/HaplotagStrategy.cpp:505-516 (HP1-1 = germline-H1 AND somatic-H3/ALT)",
}

OUT = {"validity_d_within": validity, "rederivation": rows}
json.dump(OUT, open(f"{ROOT}/genome_survey_v2/subclone_tag_rederivation.json", "w"), indent=1)

print(f"{'locus':18s}{'tag':26s}{'d_HP':>7}{'d_clean':>9}{'d_within':>9}{'d_copy':>8}{'subcl_frac':>11}")
for r in rows:
    if "error" in r:
        print(f"{r['locus']:18s}{r['tag']:26s}  {r['error']}"); continue
    f = lambda x: f"{x:.3f}" if isinstance(x, float) else str(x)
    print(f"{r['locus']:18s}{r['tag']:26s}{f(r['d_HP']):>7}{f(r['d_focal_CLEAN']):>9}{f(r['d_within']):>9}{f(r['d_copy']):>8}{f(r['subclone_copy_frac']):>11}")
nt = sum(1 for r in rows if r.get("within_testable"))
print(f"\nwithin-testable (subclone control possible): {nt}/{len(rows)}")
print(f"d_focal_CLEAN computable: {sum(1 for r in rows if r.get('d_focal_CLEAN') is not None)}/{len(rows)} (raw focal contrast works everywhere)")
print("DONE")
