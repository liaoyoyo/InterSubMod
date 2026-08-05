#!/usr/bin/env python3
"""V1: HP tag x ALT-support crosstab for HCC1395 tagged tumor BAM.
Method faithful to ReadParser::determine_alt_support_with_reason:
 walk CIGAR to base at somatic pos; base==ALT->ALT, ==REF->REF, else UNKNOWN;
 BQ<20 or in-deletion/refskip -> UNKNOWN. MAPQ>=20, no secondary/supp/unmapped.
Per-read class over all somatic sites it covers: ALT-touching (>=1 ALT),
 REF-only (>=1 REF, 0 ALT), UNKNOWN-only (covered but no confident REF/ALT).
Parallel across chromosomes. argv: chroms(comma|ALL) out_prefix [nproc]
"""
import sys, json, bisect
from collections import defaultdict
from multiprocessing import Pool
import pysam

TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
VCF_DIR = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
TP_VCF = f"{VCF_DIR}/filtered_snv_tp.vcf.gz"
FP_VCF = f"{VCF_DIR}/filtered_snv_fp.vcf.gz"
MAPQ_MIN = 20
BQ_MIN = 20

# CIGAR consumes: M=0 I=1 D=2 N=3 S=4 H=5 P=6 ==7 X=8
CONSUME_REF   = {0,2,3,7,8}
CONSUME_QUERY = {0,1,4,7,8}

def canon_hp(aln):
    if not aln.has_tag("HP"):
        return "0"
    s = str(aln.get_tag("HP"))
    return {"1":"1","2":"2","11":"1-1","21":"2-1","33":"3",
            "1-1":"1-1","2-1":"2-1","3":"3"}.get(s, s)

def hp_family(hp):
    return {"0":"unphase(0)","1":"germline HP1","2":"germline HP2",
            "1-1":"HP1-1","2-1":"HP2-1","3":"HP3"}.get(hp, f"OTHER({hp})")

def base_at(a, ref_pos0):
    """Return ('ALT'|'REF'|'UNK'|None-code, ) faithful to ReadParser.
    Returns tuple (base_char_or_None, bq_or_None, code) where code in
    {'ok','del','nocov','lowq'}. base only meaningful when code=='ok'."""
    read_start = a.reference_start
    ref_pos = read_start
    seq_pos = 0
    read_offset = -1
    for (op, length) in a.cigartuples:
        if op in (0,7,8):  # M/=/X : consume both
            if ref_pos <= ref_pos0 < ref_pos + length:
                read_offset = seq_pos + (ref_pos0 - ref_pos)
                break
            ref_pos += length; seq_pos += length
        elif op in (1,4):  # I/S : consume query only
            seq_pos += length
        elif op in (2,3):  # D/N : consume ref only
            if ref_pos <= ref_pos0 < ref_pos + length:
                return (None, None, 'del')
            ref_pos += length
        # H(5)/P(6): consume nothing
    if read_offset == -1:
        return (None, None, 'nocov')
    q = a.query_qualities
    if q is None or q[read_offset] < BQ_MIN:
        return (None, None, 'lowq')
    base = a.query_sequence[read_offset].upper()
    return (base, q[read_offset], 'ok')

def load_sites(vcf, label, chroms):
    sites = defaultdict(dict)
    tb = pysam.TabixFile(vcf)
    want = tb.contigs if chroms == ["ALL"] else [c for c in chroms if c in tb.contigs]
    n = 0
    for chrom in want:
        for row in tb.fetch(chrom):
            f = row.split("\t")
            pos = int(f[1]); ref = f[3].upper(); alt = f[4].upper()
            if len(ref) == 1 and len(alt) == 1:
                sites[chrom][pos-1] = (ref, alt, label)  # 0-based
                n += 1
    tb.close()
    return sites, n

# globals set per worker
_SITES = None

def process_chrom(chrom):
    posmap = _SITES[chrom]                 # {pos0: (ref,alt,label)}
    positions = sorted(posmap.keys())
    ct = defaultdict(lambda: defaultdict(int))       # hp -> cat -> n (union)
    ct_tp = defaultdict(lambda: defaultdict(int))    # hp -> cat -> n (TP only)
    covered_positions = set()
    bam = pysam.AlignmentFile(TBAM, "rb")
    nreads = 0; nreads_hit = 0
    for a in bam.fetch(chrom):
        nreads += 1
        if a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        if a.mapping_quality < MAPQ_MIN:
            continue
        rs = a.reference_start; re = a.reference_end
        if re is None:
            continue
        lo = bisect.bisect_left(positions, rs)
        hi = bisect.bisect_left(positions, re)
        if lo >= hi:
            continue
        nreads_hit += 1
        hp = canon_hp(a)
        # per-read flags
        af=rf=uf=0; aftp=rftp=0
        for k in range(lo, hi):
            p0 = positions[k]
            ref, alt, lab = posmap[p0]
            base, bq, code = base_at(a, p0)
            if code != 'ok':
                uf = 1
                continue
            covered_positions.add(p0)
            if base == alt:
                af = 1
                if lab == "tp": aftp = 1
            elif base == ref:
                rf = 1
                if lab == "tp": rftp = 1
            else:
                uf = 1
        cat = "ALT-touching" if af else ("REF-only" if rf else "UNKNOWN-only")
        ct[hp][cat] += 1
        cat2 = "ALT-touching" if aftp else ("REF-only" if rftp else "UNKNOWN-only")
        ct_tp[hp][cat2] += 1
    bam.close()
    # to plain dict
    return (chrom, {h: dict(v) for h,v in ct.items()},
            {h: dict(v) for h,v in ct_tp.items()},
            nreads, nreads_hit, len(covered_positions), len(positions))

def _init(sites):
    global _SITES
    _SITES = sites

def main(chroms, out_prefix, nproc):
    tp_sites, n_tp = load_sites(TP_VCF, "tp", chroms)
    fp_sites, n_fp = load_sites(FP_VCF, "fp", chroms)
    sites = defaultdict(dict)
    for d in (fp_sites, tp_sites):
        for ch, m in d.items():
            sites[ch].update(m)
    total_sites = sum(len(m) for m in sites.values())
    chrom_list = sorted(sites.keys())
    print(f"[sites] TP={n_tp} FP={n_fp} merged_positions={total_sites} chroms={len(chrom_list)}", flush=True)

    ct = defaultdict(lambda: defaultdict(int))
    ct_tp = defaultdict(lambda: defaultdict(int))
    per_chrom = {}
    with Pool(min(nproc, len(chrom_list)), initializer=_init, initargs=(dict(sites),)) as pool:
        for (chrom, c1, c2, nr, nrh, covp, np_) in pool.imap_unordered(process_chrom, chrom_list):
            for h, v in c1.items():
                for cat, n in v.items(): ct[h][cat] += n
            for h, v in c2.items():
                for cat, n in v.items(): ct_tp[h][cat] += n
            per_chrom[chrom] = {"reads_total": nr, "reads_overlap_site": nrh,
                                "positions_covered": covp, "positions_total": np_}
            print(f"[done] {chrom} reads={nr} hit={nrh} pos_cov={covp}/{np_}", flush=True)

    out = {"bam": TBAM, "somatic_vcf": {"tp": TP_VCF, "fp": FP_VCF},
           "method": "CIGAR-walk faithful to ReadParser::determine_alt_support_with_reason",
           "n_sites": {"tp": n_tp, "fp": n_fp, "merged_positions": total_sites},
           "filters": {"mapq_min": MAPQ_MIN, "bq_min": BQ_MIN,
                       "exclude": "unmapped/secondary/supplementary",
                       "hp_mapping": "1->1,2->2,11->1-1,21->2-1,33->3; no HP tag->0(unphase)"},
           "chroms": chroms, "per_chrom": per_chrom,
           "crosstab_union_TPFP": {}, "crosstab_TPonly": {}}
    hp_order = ["0","1","2","1-1","2-1","3"]
    extra = sorted([h for h in set(list(ct)+list(ct_tp)) if h not in hp_order])
    tot_reads = 0
    for hp in hp_order + extra:
        row = ct.get(hp, {})
        alt = row.get("ALT-touching",0); rfo = row.get("REF-only",0); unk = row.get("UNKNOWN-only",0)
        cov = alt+rfo+unk; tot_reads += cov
        out["crosstab_union_TPFP"][hp] = {
            "family": hp_family(hp), "ALT_touching": alt, "REF_only": rfo,
            "UNKNOWN_only": unk, "covered_total": cov,
            "pct_ALT_of_covered": round(100*alt/cov,4) if cov else None,
            "pct_REFonly_of_covered": round(100*rfo/cov,4) if cov else None}
        rt = ct_tp.get(hp, {})
        talt=rt.get("ALT-touching",0); trfo=rt.get("REF-only",0); tunk=rt.get("UNKNOWN-only",0)
        tcov=talt+trfo+tunk
        out["crosstab_TPonly"][hp] = {"family": hp_family(hp), "ALT_touching": talt,
            "REF_only": trfo, "UNKNOWN_only": tunk, "covered_total": tcov,
            "pct_ALT_of_covered": round(100*talt/tcov,4) if tcov else None}
    out["n_covered_reads_union"] = tot_reads
    json.dump(out, open(out_prefix+".json","w"), ensure_ascii=False, indent=1)
    with open(out_prefix+"_crosstab.tsv","w") as f:
        f.write("view\thp_raw\tfamily\tALT_touching\tREF_only\tUNKNOWN_only\tcovered_total\tpct_ALT_of_covered\tpct_REFonly_of_covered\n")
        for hp in hp_order+extra:
            d=out["crosstab_union_TPFP"][hp]
            f.write(f"union_TPFP\t{hp}\t{d['family']}\t{d['ALT_touching']}\t{d['REF_only']}\t{d['UNKNOWN_only']}\t{d['covered_total']}\t{d['pct_ALT_of_covered']}\t{d['pct_REFonly_of_covered']}\n")
        for hp in hp_order+extra:
            d=out["crosstab_TPonly"][hp]
            f.write(f"TPonly\t{hp}\t{d['family']}\t{d['ALT_touching']}\t{d['REF_only']}\t{d['UNKNOWN_only']}\t{d['covered_total']}\t{d['pct_ALT_of_covered']}\t\n")
    print("WROTE", out_prefix+".json", flush=True)
    print("=== UNION TP+FP ===", flush=True)
    for hp in hp_order:
        d=out["crosstab_union_TPFP"][hp]
        print(f"  {hp:4s} {d['family']:14s} ALT={d['ALT_touching']:8d} REFonly={d['REF_only']:8d} UNK={d['UNKNOWN_only']:7d} cov={d['covered_total']:8d} pctALT={d['pct_ALT_of_covered']} pctREFonly={d['pct_REFonly_of_covered']}", flush=True)

if __name__ == "__main__":
    chroms = sys.argv[1].split(",")
    out_prefix = sys.argv[2]
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    main(chroms, out_prefix, nproc)
