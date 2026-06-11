#!/usr/bin/env python3
"""
65_dwithin_validity_reproducible.py — REPRODUCIBLE d_within validity check.

Fixes the provenance gap flagged by the capstone adversary: script 64 hard-coded
'19/19 carry other somatic ALT' as a string literal. This re-derives it from the
tumor BAM + somatic VCF so the number is grep-able / re-runnable.

For BRCA2 (chr13:32315128): from the level1 cache get HP1-1/ref read_ids, then from
the BAM check (a) focal base + BQ, (b) how many OTHER somatic SNVs (from somatic_pass.vcf.gz,
±25kb) each read carries ALT at. HP1-1 reads MUST carry somatic ALT somewhere
(longphase-S HaplotagStrategy.cpp:505,627-628,658-660: hpCount[3]++ only at base==TumorAltBase),
so this confirms HP1-1/ref reads are genuine subclone reads -> d_within is subclone-controlled.
"""
import pysam, glob, gzip, json
import numpy as np
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
DIR = ("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
       "20260314_HCC1395_paired_full_full_complete_matrix/longphase_s")
BAM = f"{DIR}/HCC1395_tagged.bam"
VCF = f"{DIR}/somatic_pass.vcf.gz"
CHROM, SPOS = "chr13", 32315128

# somatic SNVs near BRCA2 (+-25kb) from somatic_pass.vcf.gz (no tabix -> linear scan)
sompos = {}
with gzip.open(VCF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if p[0] == CHROM and len(p) >= 5 and len(p[3]) == 1 and len(p[4]) == 1:
            pos = int(p[1])
            if SPOS - 25000 <= pos <= SPOS + 25000:
                sompos[pos] = (p[3], p[4])

# HP1-1/ref read_ids from level1 cache
cache = glob.glob(f"{ROOT}/pipeline/cache/level1/{CHROM}_{SPOS}_*.tsv.gz")[0]
hp11ref = set()
with gzip.open(cache, "rt") as f:
    h = f.readline().rstrip().split("\t"); ix = {c: i for i, c in enumerate(h)}
    for line in f:
        p = line.rstrip().split("\t")
        if (p[ix["bam_source"]] == "tumor" and p[ix["haplotype_tag"]] == "1-1"
                and p[ix["somatic_allele_type"]] == "ref"):
            hp11ref.add(p[ix["read_id"]])

bam = pysam.AlignmentFile(BAM)
focal_ref, focal_alt = sompos.get(SPOS, ("G", "A"))
bqs, focal_bases, carry = [], {}, []
seen = set()
for r in bam.fetch(CHROM, SPOS - 25000, SPOS + 25000):
    if r.query_name not in hp11ref or r.is_secondary or r.is_supplementary or r.query_name in seen:
        continue
    seen.add(r.query_name)
    r2q = {rp: q for q, rp in r.get_aligned_pairs(matches_only=True)}
    # focal base + BQ
    if (SPOS - 1) in r2q:
        q = r2q[SPOS - 1]
        focal_bases[r.query_sequence[q].upper()] = focal_bases.get(r.query_sequence[q].upper(), 0) + 1
        bqs.append(r.query_qualities[q])
    # other somatic ALT count
    nalt = ncov = 0
    for pos, (ref, alt) in sompos.items():
        if pos == SPOS:
            continue
        if (pos - 1) in r2q:
            ncov += 1
            if r.query_sequence[r2q[pos - 1]].upper() == alt.upper():
                nalt += 1
    carry.append({"read": r.query_name, "n_other_somatic_alt": nalt, "n_other_somatic_cov": ncov})

n = len(carry)
nz = sum(1 for c in carry if c["n_other_somatic_alt"] > 0)
OUT = {
    "locus": f"{CHROM}:{SPOS}", "ref": focal_ref, "alt": focal_alt,
    "n_other_somatic_snv_near": len(sompos) - 1, "other_somatic_positions": sorted(p for p in sompos if p != SPOS),
    "n_HP11_ref_reads_checked": n,
    "focal_base_distribution": focal_bases,
    "focal_BQ_median": round(float(np.median(bqs)), 1) if bqs else None,
    "focal_BQ_min": int(min(bqs)) if bqs else None,
    "carry_other_somatic_ALT": f"{nz}/{n}",
    "mean_other_somatic_alt_per_read": round(float(np.mean([c["n_other_somatic_alt"] for c in carry])), 2) if carry else None,
    "verdict": ("VALID — HP1-1/ref reads are genuine somatic-subclone reads (focal REF + carry other somatic ALT) "
                "-> d_within is a subclone-controlled focal-allele estimate (QUALITATIVE; % quantification NOT robust)"
                if carry and nz >= n * 0.5 else "CHECK — many HP1-1/ref reads lack other somatic ALT"),
}
json.dump(OUT, open(f"{ROOT}/genome_survey_v2/dwithin_validity.json", "w"), indent=1)
print(json.dumps(OUT, ensure_ascii=False, indent=1))
print("DONE")
