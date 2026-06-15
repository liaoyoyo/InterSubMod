#!/usr/bin/env python3
"""
seqc2_subclone_concordance.py — compare our subclone-defining loci against the
SEQC2 HCC1395 external truth (somatic-SNV VAF + CNV/LOH benchmark).

WHAT IT DOES (honest, 3-layer concordance — NOT a 1:1 gold-standard match):
  L1  CN/LOH concordance : is each locus inside a SEQC2 gain/loss/LOH segment?
  L2  SNV clonal/subclonal: look up locus in SEQC2 high-confidence sSNV truth,
       read consensus TVAF (tumor VAF) -> approximate CCF = TVAF*CN/(purity*mult).
  HC  evaluability       : is the locus inside a SEQC2 High-Confidence region?
       (out-of-HC => truth-unevaluable, e.g. chr2:18,086,020).

CCF is an APPROXIMATION (multiplicity=1, CN from SEQC2 median; copy-neutral LOH
assumed CN=2). Cross-platform caveat: SEQC2 TVAF is short-read bulk; our reads are
ONT. Use as concordance/sanity check, not literal per-read truth. SEQC2's own
subclone tree is a single PhyloWGS pipeline (DREAM: 19-35% algorithm variance) —
it is a reference context, not a gold standard.

Reusable genome-wide: feed any loci TSV (id,chrom,pos1,ref,alt,our_label).

Usage:
  python3 seqc2_subclone_concordance.py --loci demo_loci_chr2_18M.tsv \
      --out-prefix ../data/chr2_18M_seqc2_concordance
"""
import argparse
import bisect
import os
import re
import sys

SEQC2 = "/big8_disk/data/HCC1395/SEQC2"
DEFAULTS = {
    "ssnv_vcf": f"{SEQC2}/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf",
    "hc_bed":   f"{SEQC2}/High-Confidence_Regions_v1.2.bed",
    "cnv_bed":  f"{SEQC2}/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed",
    "gain_cn":  f"{SEQC2}/CNV/Additional_file_3_cnv_gain_cn_median.txt",
    "loss_cn":  f"{SEQC2}/CNV/Additional_file_4_cnv_loss_cn_median.txt",
}
TVAF_RE = re.compile(r"(?:^|;)TVAF=([0-9.]+)")
NVAF_RE = re.compile(r"(?:^|;)NVAF=([0-9.]+)")
CONF_RE = re.compile(r"(HighConf|MedConf|LowConf)")


def load_loci(path):
    loci = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith("#"):
                continue
            f = ln.split("\t")
            loci.append({
                "id": f[idx["id"]], "chrom": f[idx["chrom"]],
                "pos1": int(f[idx["pos1"]]), "ref": f[idx["ref"]],
                "alt": f[idx["alt"]], "our_label": f[idx["our_label"]],
            })
    return loci


def scan_ssnv(vcf, want):
    """want = set of (chrom,pos1). Returns {(chrom,pos1): dict}."""
    hits = {}
    chroms = {c for c, _ in want}
    with open(vcf) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if f[0] not in chroms:
                continue
            try:
                pos = int(f[1])
            except (IndexError, ValueError):
                continue
            key = (f[0], pos)
            if key not in want:
                continue
            info = f[7] if len(f) > 7 else ""
            filt = f[6] if len(f) > 6 else ""
            tv = TVAF_RE.search(info)
            nv = NVAF_RE.search(info)
            cf = CONF_RE.search(filt) or CONF_RE.search(info)
            hits[key] = {
                "ref": f[3], "alt": f[4],
                "conf": cf.group(1) if cf else "PASS",
                "tvaf": float(tv.group(1)) if tv else None,
                "nvaf": float(nv.group(1)) if nv else None,
            }
    return hits


def load_intervals(path, chroms, val_col=None):
    """0-based half-open BED. Returns {chrom: (starts[], recs[])} sorted by start."""
    by = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(("#", "track", "browser")):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] not in chroms:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            val = f[val_col] if (val_col is not None and len(f) > val_col) else None
            by.setdefault(f[0], []).append((s, e, val))
    for c in by:
        by[c].sort()
    return by


def query(by, chrom, pos1):
    """Return list of overlapping records for 1-based pos1 (BED 0-based)."""
    p0 = pos1 - 1
    recs = by.get(chrom, [])
    if not recs:
        return []
    starts = [r[0] for r in recs]
    i = bisect.bisect_right(starts, p0)
    out = []
    for j in range(i - 1, -1, -1):
        s, e, v = recs[j]
        if e <= p0:
            # since sorted by start, earlier ones may still overlap if long; scan a small window
            if s < p0 - 30_000_000:
                break
            continue
        if s <= p0 < e:
            out.append((s, e, v))
    return out


def classify_ccf(ccf):
    if ccf is None:
        return "n/a"
    if ccf >= 0.85:
        return "clonal"
    if ccf >= 0.25:
        return "subclonal-major"
    if ccf > 0.0:
        return "subclonal-minor"
    return "n/a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--purity", type=float, default=0.99)
    ap.add_argument("--mult", type=int, default=1, help="assumed mutation multiplicity")
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k.replace('_', '-')}", default=v)
    a = ap.parse_args()

    for k in DEFAULTS:
        p = getattr(a, k)
        if not os.path.exists(p):
            sys.exit(f"[ERR] missing truth file ({k}): {p}")

    loci = load_loci(a.loci)
    if not loci:
        sys.exit("[ERR] no loci loaded")
    chroms = {l["chrom"] for l in loci}

    want = {(l["chrom"], l["pos1"]) for l in loci}
    ssnv = scan_ssnv(a.ssnv_vcf, want)
    hc = load_intervals(a.hc_bed, chroms)
    cnv = load_intervals(a.cnv_bed, chroms, val_col=3)
    gain = load_intervals(a.gain_cn, chroms, val_col=3)
    loss = load_intervals(a.loss_cn, chroms, val_col=3)

    rows = []
    for l in loci:
        key = (l["chrom"], l["pos1"])
        s = ssnv.get(key)
        in_hc = bool(query(hc, *key))
        cnv_hit = query(cnv, *key)
        cn_type = cnv_hit[0][2] if cnv_hit else "diploid?"
        # determine integer CN
        if cn_type == "gain":
            g = query(gain, *key)
            cn = int(g[0][2]) if g else 3
        elif cn_type == "loss":
            ls = query(loss, *key)
            cn = int(ls[0][2]) if ls else 1
        else:  # loh (copy-neutral) or no record
            cn = 2
        tvaf = s["tvaf"] if s else None
        ccf = None
        if tvaf is not None:
            ccf = round(tvaf * cn / (a.purity * a.mult), 3)
        # truth status
        if s:
            status = f"SEQC2-{s['conf']}"
        elif in_hc:
            status = "in-HC,not-in-truth(likely-FP/germline)"
        else:
            status = "out-of-HC(truth-unevaluable)"
        rows.append({
            "id": l["id"], "chrom": l["chrom"], "pos1": l["pos1"],
            "ref": l["ref"], "alt": l["alt"], "our_label": l["our_label"],
            "seqc2_status": status, "in_HC": "Y" if in_hc else "N",
            "cn_type": cn_type, "cn": cn,
            "tvaf": tvaf if tvaf is not None else ".",
            "nvaf": s["nvaf"] if s and s["nvaf"] is not None else ".",
            "ccf_est": ccf if ccf is not None else ".",
            "clonal_class": classify_ccf(ccf),
        })

    cols = ["id", "chrom", "pos1", "ref", "alt", "our_label", "seqc2_status",
            "in_HC", "cn_type", "cn", "tvaf", "nvaf", "ccf_est", "clonal_class"]
    tsv_path = a.out_prefix + ".tsv"
    with open(tsv_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # region-level LOH/CN concordance summary
    pmin = min(l["pos1"] for l in loci)
    pmax = max(l["pos1"] for l in loci)
    chrom0 = loci[0]["chrom"]
    region_cnv = query(cnv, chrom0, pmin) + query(cnv, chrom0, pmax)
    region_types = sorted({v for _, _, v in region_cnv})
    n_truth = sum(1 for r in rows if r["seqc2_status"].startswith("SEQC2-"))
    n_eval = sum(1 for r in rows if r["in_HC"] == "Y")
    n_sub = sum(1 for r in rows if r["clonal_class"].startswith("subclonal"))

    md = a.out_prefix + ".md"
    with open(md, "w") as fh:
        fh.write(f"# SEQC2 concordance — {a.loci}\n\n")
        fh.write(f"- loci: {len(loci)} | region {chrom0}:{pmin:,}-{pmax:,}\n")
        fh.write(f"- region SEQC2 CN/LOH type(s): {', '.join(region_types) or 'none'}\n")
        fh.write(f"- in SEQC2 sSNV truth: {n_truth}/{len(loci)} | "
                 f"in HC (evaluable): {n_eval}/{len(loci)} | "
                 f"subclonal (CCF<0.85): {n_sub}/{len(loci)}\n")
        fh.write(f"- purity={a.purity}, multiplicity={a.mult}; "
                 f"CCF≈TVAF*CN/(purity*mult); LOH assumed CN=2\n\n")
        fh.write("| id | locus | our_label | SEQC2 status | inHC | CN | TVAF | CCF≈ | class |\n")
        fh.write("|----|-------|-----------|--------------|:--:|:--:|:--:|:--:|------|\n")
        for r in rows:
            fh.write(f"| {r['id']} | {r['chrom']}:{r['pos1']:,} {r['ref']}>{r['alt']} | "
                     f"{r['our_label']} | {r['seqc2_status']} | {r['in_HC']} | "
                     f"{r['cn_type']}({r['cn']}) | {r['tvaf']} | {r['ccf_est']} | {r['clonal_class']} |\n")

    print(f"[done] {tsv_path}")
    print(f"[done] {md}")
    # echo table to stdout for immediate inspection
    print("\n--- per-locus ---")
    print("\t".join(cols))
    for r in rows:
        print("\t".join(str(r[c]) for c in cols))
    print(f"\n--- region summary ---")
    print(f"region={chrom0}:{pmin}-{pmax} CN/LOH={region_types} "
          f"in_truth={n_truth}/{len(loci)} evaluable={n_eval}/{len(loci)} "
          f"subclonal={n_sub}/{len(loci)}")


if __name__ == "__main__":
    main()
