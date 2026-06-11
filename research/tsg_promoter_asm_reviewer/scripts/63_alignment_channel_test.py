#!/usr/bin/env python3
"""
63_alignment_channel_test.py — Test ③: does alignment difficulty explain the 6
UNTESTABLE strongest signals (pure-ALT tag, CGI-desert)? Closes the 06-02 MAPQ blind
spot (MAPQ saturated at 60) using soft-clip rate / secondary+supplementary / NM-per-kb.

Per locus window (±1kb): primary-read soft-clip frac, NM/kb, MAPQ, n; secondary/supp
record counts; per-tag (HP1 vs HP1-1 from cache read_id->tag map) alignment delta.
Compare groups: untestable (6) vs control (chr17 clean / BRCA2,chr5 copy) vs neutral baseline (15).

Runs in background (region-restricted fetch on the 278G tumor BAM; minutes).
"""
import json, glob, gzip
import numpy as np
import pysam
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
TUMOR = ("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
         "20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam")
WIN = 1000
manifest = json.load(open(f"{ROOT}/genome_survey_v2/alignment_test_manifest.json"))


def read_tag_map(chrom, pos):
    """read_id -> haplotype_tag (1/1-1/2) from the cache (MSA-assigned somatic sub-tag)."""
    c = glob.glob(f"{ROOT}/pipeline/cache/level1/{chrom}_{pos}_*.tsv.gz")
    if not c:
        return {}
    m = {}
    with gzip.open(c[0], "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {x: i for i, x in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[ix["bam_source"]] == "tumor":
                m[p[ix["read_id"]]] = p[ix["haplotype_tag"]]
    return m


def softclip_frac(r):
    if not r.cigartuples:
        return None
    s = sum(l for op, l in r.cigartuples if op == 4)        # soft
    aln = r.query_alignment_length or 0
    tot = s + aln
    return s / tot if tot else None


def nm_per_kb(r):
    try:
        nm = r.get_tag("NM")
    except KeyError:
        return None
    aln = r.query_alignment_length or 0
    return (nm / (aln / 1000.0)) if aln >= 200 else None


def analyze(locus):
    chrom, pos = locus["chrom"], locus["pos"]
    tagmap = read_tag_map(chrom, pos)
    bam = pysam.AlignmentFile(TUMOR)
    prim, sec, supp = 0, 0, 0
    sc, nm, mapq = [], [], []
    by_tag = {"1": {"sc": [], "nm": []}, "1-1": {"sc": [], "nm": []}}
    for r in bam.fetch(chrom, max(0, pos - WIN), pos + WIN):
        if r.is_unmapped:
            continue
        if r.is_secondary:
            sec += 1; continue
        if r.is_supplementary:
            supp += 1; continue
        prim += 1
        s = softclip_frac(r); n = nm_per_kb(r)
        if s is not None: sc.append(s)
        if n is not None: nm.append(n)
        mapq.append(r.mapping_quality)
        t = tagmap.get(r.query_name)
        if t in by_tag:
            if s is not None: by_tag[t]["sc"].append(s)
            if n is not None: by_tag[t]["nm"].append(n)
    bam.close()
    def m(x): return round(float(np.mean(x)), 4) if x else None
    out = {"locus": locus["locus"], "group": locus["group"], "d_cis": locus.get("d_cis"),
           "total_cn": locus.get("total_cn"), "n_primary": prim, "n_secondary": sec, "n_supplementary": supp,
           "frac_supp_records": round(supp / max(prim, 1), 4), "frac_sec_records": round(sec / max(prim, 1), 4),
           "mean_softclip_frac": m(sc), "mean_NM_per_kb": m(nm), "mean_MAPQ": m(mapq),
           "HP1_softclip": m(by_tag["1"]["sc"]), "HP11_softclip": m(by_tag["1-1"]["sc"]),
           "HP1_NMkb": m(by_tag["1"]["nm"]), "HP11_NMkb": m(by_tag["1-1"]["nm"])}
    # HP1-1 vs HP1 alignment delta (positive = somatic tag worse aligned = alignment-driven cluster)
    if out["HP1_softclip"] is not None and out["HP11_softclip"] is not None:
        out["d_softclip_HP11_minus_HP1"] = round(out["HP11_softclip"] - out["HP1_softclip"], 4)
    if out["HP1_NMkb"] is not None and out["HP11_NMkb"] is not None:
        out["d_NMkb_HP11_minus_HP1"] = round(out["HP11_NMkb"] - out["HP1_NMkb"], 4)
    return out


results = []
for i, loc in enumerate(manifest):
    try:
        results.append(analyze(loc))
        print(f"[{i+1}/{len(manifest)}] {loc['locus']} ({loc['group']}) done", flush=True)
    except Exception as e:
        print(f"[{i+1}/{len(manifest)}] {loc['locus']} ERROR {str(e)[:100]}", flush=True)

# group aggregates
def gstat(grp, key):
    v = [r[key] for r in results if r["group"] == grp and r.get(key) is not None]
    return {"n": len(v), "mean": round(float(np.mean(v)), 4) if v else None,
            "median": round(float(np.median(v)), 4) if v else None}
agg = {}
for grp in ("untestable", "control", "neutral_baseline"):
    agg[grp] = {k: gstat(grp, k) for k in ("mean_softclip_frac", "mean_NM_per_kb", "frac_supp_records", "frac_sec_records",
                                            "d_softclip_HP11_minus_HP1", "d_NMkb_HP11_minus_HP1")}
OUT = {"per_locus": results, "group_aggregate": agg}
json.dump(OUT, open(f"{ROOT}/genome_survey_v2/alignment_channel_test.json", "w"), indent=1)
print("\n=== group aggregate (mean of per-locus means) ===")
for grp in ("untestable", "control", "neutral_baseline"):
    a = agg[grp]
    print(f"  {grp:18s} softclip={a['mean_softclip_frac']['mean']} NM/kb={a['mean_NM_per_kb']['mean']} "
          f"supp_rate={a['frac_supp_records']['mean']} sec_rate={a['frac_sec_records']['mean']}")
print("DONE_ALIGNMENT_TEST")
