#!/usr/bin/env python3
"""
Extract per-read base x SNV x methylation matrix at chr2:18.06M-18.11M for HCC1395
subclone reconstruction verification.

For each read overlapping the region, record:
  - HP / PS tag, strand, ref span
  - base call at each of the 6 candidate sSNV positions (DEL / None if not covered)
  - 5mC / 5hmC modification probability at user-specified CpG sites AND region-wide CpGs

Usage:
  extract_locus_matrix.py <tumor_bam> <normal_bam> <out_json> [label]

Output is JSON; NO interpretation here (pure extraction, per data-integrity rule).
"""
import sys, json
import pysam

# ---- 1-based reference coordinates (from user manual observation) ----
SNV_1BASED = {
    "1": (18068480, "C", "G"),   # TP HighConf
    "2": (18072546, "G", "C"),   # TP HighConf
    "3": (18086020, "G", "?"),   # FP? (HC gap)
    "4": (18096269, "C", "G/T/DEL"),  # TP MedConf, multi-allele
    "5": (18099697, "G", "C"),   # TP (low VAF)
    "6": (18108828, "C", "G"),   # TP HighConf
}
CPG_1BASED = {
    "2.1": 18073987,
    "2.2": 18076409,
    "3.1": 18090947,
    "3.2": 18091584,
    "3.3": 18093414,
    "3.4": 18094114,
    "3.5": 18096033,
    "4.1": 18096041,
    "5.1": 18100544,
    "5.2": 18100683,
}
CHROM = "chr2"
REGION_START0 = 18068480 - 2000   # pad
REGION_END0   = 18108828 + 2000


def extract(bam_path, snv0, cpg0, region_wide_cpg=True):
    bam = pysam.AlignmentFile(bam_path, "rb")
    records = []
    # gather region-wide methylation ref positions too (for independent association analysis)
    for read in bam.fetch(CHROM, REGION_START0, REGION_END0):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        seq = read.query_sequence
        # build ref->query and query->ref maps
        ref2query = {}
        query2ref = {}
        for qp, rp in read.get_aligned_pairs(matches_only=False):
            if rp is not None:
                ref2query[rp] = qp      # qp may be None -> deletion at this ref pos
            if qp is not None and rp is not None:
                query2ref[qp] = rp
        rec = {
            "name": read.query_name,
            "HP": (read.get_tag("HP") if read.has_tag("HP") else None),
            "PS": (read.get_tag("PS") if read.has_tag("PS") else None),
            "strand": ("-" if read.is_reverse else "+"),
            "ref_start": read.reference_start,
            "ref_end": read.reference_end,
            "mapq": read.mapping_quality,
            "bases": {},
            "meth": {},
            "meth_all": {},   # ref_pos -> {m:prob,h:prob}  region-wide
        }
        # SNV bases
        for label, p0 in snv0.items():
            if p0 in ref2query:
                qp = ref2query[p0]
                rec["bases"][label] = ("DEL" if qp is None else seq[qp])
            else:
                rec["bases"][label] = None
        # methylation: build refpos -> mod prob from modified_bases
        modbases = read.modified_bases or {}
        refpos_mod = {}
        for (canon, mstrand, mod), calls in modbases.items():
            for qp, qual in calls:
                rp = query2ref.get(qp)
                if rp is None:
                    continue
                refpos_mod.setdefault(rp, {})[mod] = qual
        # user CpG sites (check site, +1, -1 to catch reverse-strand CpG dinuc)
        for label, p0 in cpg0.items():
            found = None; matched_rp = None
            for cand in (p0, p0 + 1, p0 - 1):
                if cand in refpos_mod:
                    found = refpos_mod[cand]; matched_rp = cand; break
            if found is not None:
                m = found.get("m"); h = found.get("h")
                rec["meth"][label] = {
                    "ref_pos1": matched_rp + 1,
                    "m_prob": (round(m / 255.0, 3) if m is not None else None),
                    "h_prob": (round(h / 255.0, 3) if h is not None else None),
                }
            else:
                rec["meth"][label] = None
        # region-wide CpG methylation (for independent associated-site discovery)
        if region_wide_cpg:
            for rp, mods in refpos_mod.items():
                m = mods.get("m"); h = mods.get("h")
                rec["meth_all"][str(rp + 1)] = {
                    "m": (round(m / 255.0, 3) if m is not None else None),
                    "h": (round(h / 255.0, 3) if h is not None else None),
                }
        records.append(rec)
    bam.close()
    return records


def main():
    tumor_bam, normal_bam, out_json = sys.argv[1], sys.argv[2], sys.argv[3]
    label = sys.argv[4] if len(sys.argv) > 4 else "sample"
    snv0 = {k: (v[0] - 1) for k, v in SNV_1BASED.items()}
    cpg0 = {k: (v - 1) for k, v in CPG_1BASED.items()}

    out = {
        "label": label,
        "region": f"{CHROM}:{REGION_START0+1}-{REGION_END0}",
        "snv_1based": {k: {"pos": v[0], "ref": v[1], "alt": v[2]} for k, v in SNV_1BASED.items()},
        "cpg_1based": CPG_1BASED,
        "tumor_bam": tumor_bam,
        "normal_bam": normal_bam,
    }
    out["tumor_reads"] = extract(tumor_bam, snv0, cpg0, region_wide_cpg=True)
    # normal: bases only (no need region-wide meth for normal, but keep user CpGs)
    out["normal_reads"] = extract(normal_bam, snv0, cpg0, region_wide_cpg=False)

    with open(out_json, "w") as f:
        json.dump(out, f)
    print(f"[OK] {label}: tumor_reads={len(out['tumor_reads'])} normal_reads={len(out['normal_reads'])} -> {out_json}")


if __name__ == "__main__":
    main()
