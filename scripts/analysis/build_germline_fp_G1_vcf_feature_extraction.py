#!/usr/bin/env python3
"""G1: Extract under-utilized VCF features from ClairS-TO raw VCFs for germline FP identification.

Extracts strand counts (FAU/FCU/FGU/FTU/RAU/RCU/RGU/RTU), SB, GT, H flag,
individual PoN hits, trinucleotide context, and mutation type for all 7 samples.
Only PASS + truth-matched (TP/FP) variants are included.

Output:
  {WORKSPACE}/G1_vcf_features/{sample}_features.tsv.gz   (per sample)
  {WORKSPACE}/G1_vcf_features/all_samples_merged.tsv.gz  (combined)
"""

from __future__ import annotations

import gzip
import math
import sys
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pysam

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

sys.path.insert(0, str(Path(__file__).resolve().parent))

from big7_output_catalog import SAMPLE_SPECS, SampleSpec
from observation_common import OUTPUT_ROOT, SAMPLE_ORDER, ensure_dir

WORKSPACE = ensure_dir(OUTPUT_ROOT / "20260401_germline_fp_identification")
G1_DIR = ensure_dir(WORKSPACE / "G1_vcf_features")

REF_FASTA = Path("/big8_disk/zhenyu112/GRCh38_no_alt_analysis_set.fasta")

N_WORKERS = min(7, max(1, cpu_count() - 2))

# ---------------------------------------------------------------------------
# Column schema
# ---------------------------------------------------------------------------

G1_COLUMNS = [
    "sample", "variant_key", "chrom", "pos", "ref", "alt", "truth_label",
    # caller basics
    "qual", "gq", "dp", "af", "ad_ref", "ad_alt",
    # GT
    "gt_field", "is_het", "is_hom_alt",
    # flags
    "h_flag", "pon_1", "pon_2", "pon_3", "pon_4", "pon_count",
    "verdict_germline", "verdict_somatic", "verdict_subclonal",
    # FILTER
    "filter_raw", "filter_is_pass",
    # strand counts
    "fau", "fcu", "fgu", "ftu", "rau", "rcu", "rgu", "rtu",
    "sb_clairs",
    # derived strand
    "forward_ref", "reverse_ref", "forward_alt", "reverse_alt",
    "sb_ratio", "sb_fisher_p",
    # trinucleotide context
    "trinuc_context", "mutation_type_6class", "is_transition",
    "is_cpg_site", "is_cpg_ct",
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _to_int(val, default: int = 0) -> int:
    if val is None:
        return default
    try:
        return int(val)
    except (ValueError, TypeError):
        return default


def _to_float(val, default: float = float("nan")) -> float:
    if val is None:
        return default
    try:
        v = float(val)
        return v if math.isfinite(v) else default
    except (ValueError, TypeError):
        return default


def variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}:{alt}"


# ---------------------------------------------------------------------------
# Truth set loading
# ---------------------------------------------------------------------------


def load_truth_keys_from_vcf(vcf_path: str, bed_path: str = "") -> Set[str]:
    """Load truth SNV keys (chr:pos:ref:alt) from truth VCF, optionally restricted to BED."""
    keys: Set[str] = set()
    vf = pysam.VariantFile(vcf_path)

    if bed_path and Path(bed_path).exists():
        bed_intervals = _load_bed(Path(bed_path))
        for chrom in sorted(bed_intervals):
            if chrom not in vf.header.contigs:
                continue
            for start, end in bed_intervals[chrom]:
                for rec in vf.fetch(chrom, start, end):
                    if len(rec.ref) != 1:
                        continue
                    if not rec.alts or len(rec.alts) != 1 or len(rec.alts[0]) != 1:
                        continue
                    keys.add(variant_key(rec.chrom, rec.pos, rec.ref, rec.alts[0]))
    else:
        for rec in vf.fetch():
            if len(rec.ref) != 1:
                continue
            if not rec.alts or len(rec.alts) != 1 or len(rec.alts[0]) != 1:
                continue
            keys.add(variant_key(rec.chrom, rec.pos, rec.ref, rec.alts[0]))

    vf.close()
    print(f"  [truth] Loaded {len(keys):,} truth SNV keys from {Path(vcf_path).name}")
    return keys


def _load_bed(path: Path) -> Dict[str, List[Tuple[int, int]]]:
    by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    with path.open("r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            by_chrom.setdefault(chrom, []).append((int(parts[1]), int(parts[2])))
    return by_chrom


# ---------------------------------------------------------------------------
# Trinucleotide context
# ---------------------------------------------------------------------------

COMPLEMENT = str.maketrans("ACGT", "TGCA")

# Canonical pyrimidine-context 6-class mapping
MUTATION_6CLASS = {
    ("C", "A"): "C>A", ("G", "T"): "C>A",
    ("C", "G"): "C>G", ("G", "C"): "C>G",
    ("C", "T"): "C>T", ("G", "A"): "C>T",
    ("T", "A"): "T>A", ("A", "T"): "T>A",
    ("T", "C"): "T>C", ("A", "G"): "T>C",
    ("T", "G"): "T>G", ("A", "C"): "T>G",
}

TRANSITIONS = {"C>T", "T>C"}


def get_trinuc_context(fasta: pysam.FastaFile, chrom: str, pos_1based: int) -> str:
    """Return 3-mer context (e.g., ACG) centered on the variant position."""
    start = pos_1based - 2  # 0-based, one base before
    end = pos_1based + 1    # 0-based, one base after (exclusive for fetch)
    try:
        seq = fasta.fetch(chrom, start, end).upper()
    except (ValueError, KeyError):
        return "NNN"
    if len(seq) != 3:
        return "NNN"
    return seq


def is_cpg_site_check(trinuc: str, ref: str) -> bool:
    """Check if variant is at a CpG site (ref=C with next=G, or ref=G with prev=C)."""
    if len(trinuc) != 3:
        return False
    if ref == "C" and trinuc[2] == "G":
        return True
    if ref == "G" and trinuc[0] == "C":
        return True
    return False


def is_cpg_ct_check(ref: str, alt: str, trinuc: str) -> bool:
    """Check if this is a CpG>TpG transition (C>T at CpG, or G>A at CpG on reverse)."""
    if ref == "C" and alt == "T" and len(trinuc) == 3 and trinuc[2] == "G":
        return True
    if ref == "G" and alt == "A" and len(trinuc) == 3 and trinuc[0] == "C":
        return True
    return False


# ---------------------------------------------------------------------------
# Strand bias helpers
# ---------------------------------------------------------------------------


def compute_sb_ratio(forward_alt: int, reverse_alt: int) -> float:
    total = forward_alt + reverse_alt
    if total == 0:
        return float("nan")
    return forward_alt / total


def compute_sb_fisher(forward_ref: int, reverse_ref: int,
                      forward_alt: int, reverse_alt: int) -> float:
    """Fisher exact test p-value for strand bias (2x2 table)."""
    from scipy.stats import fisher_exact
    table = [[forward_ref, reverse_ref], [forward_alt, reverse_alt]]
    if sum(sum(row) for row in table) == 0:
        return float("nan")
    try:
        _, p = fisher_exact(table)
        return p
    except Exception:
        return float("nan")


# ---------------------------------------------------------------------------
# Per-record feature extraction
# ---------------------------------------------------------------------------


def _print(msg: str) -> None:
    """Print with flush for unbuffered output."""
    print(msg, flush=True)


def _get_strand_count(info, key: str) -> int:
    val = info.get(key)
    if val is None:
        return 0
    if isinstance(val, (list, tuple)):
        return int(val[0]) if val else 0
    return int(val)


def _get_ref_alt_strand_counts(
    ref: str, alt: str,
    fau: int, fcu: int, fgu: int, ftu: int,
    rau: int, rcu: int, rgu: int, rtu: int,
) -> Tuple[int, int, int, int]:
    """Map strand counts to forward_ref, reverse_ref, forward_alt, reverse_alt."""
    base_to_forward = {"A": fau, "C": fcu, "G": fgu, "T": ftu}
    base_to_reverse = {"A": rau, "C": rcu, "G": rgu, "T": rtu}
    forward_ref = base_to_forward.get(ref, 0)
    reverse_ref = base_to_reverse.get(ref, 0)
    forward_alt = base_to_forward.get(alt, 0)
    reverse_alt = base_to_reverse.get(alt, 0)
    return forward_ref, reverse_ref, forward_alt, reverse_alt


def extract_features_from_record(
    rec: pysam.VariantRecord,
    fasta: pysam.FastaFile,
) -> Dict:
    """Extract all G1 features from a single VCF record."""
    chrom = rec.chrom
    pos = rec.pos
    ref_allele = rec.ref
    alt_allele = rec.alts[0] if rec.alts else "."

    # FORMAT fields
    sample_data = list(rec.samples.values())[0] if rec.samples else {}
    gt_tuple = sample_data.get("GT", (None, None))
    if gt_tuple and gt_tuple != (None, None):
        gt_str = "/".join(str(a) if a is not None else "." for a in gt_tuple)
        is_het = (gt_tuple[0] != gt_tuple[1]) if len(gt_tuple) >= 2 else False
        is_hom_alt = (gt_tuple[0] == gt_tuple[1] and gt_tuple[0] is not None and gt_tuple[0] > 0)
    else:
        gt_str = "./."
        is_het = False
        is_hom_alt = False

    # AD field
    ad = sample_data.get("AD")
    if ad is not None and isinstance(ad, (list, tuple)) and len(ad) >= 2:
        ad_ref, ad_alt = int(ad[0]), int(ad[1])
    else:
        ad_ref, ad_alt = 0, 0

    # INFO flags
    info = rec.info
    h_flag = "H" in info
    pon_1 = "PoN_1" in info
    pon_2 = "PoN_2" in info
    pon_3 = "PoN_3" in info
    pon_4 = "PoN_4" in info
    pon_count = int(pon_1) + int(pon_2) + int(pon_3) + int(pon_4)

    verdict_germline = "Verdict_Germline" in info
    verdict_somatic = "Verdict_Somatic" in info
    verdict_subclonal = "Verdict_SubclonalSomatic" in info

    # FILTER
    filter_terms = list(rec.filter.keys()) or ["PASS"]
    filter_raw = ";".join(filter_terms)
    filter_is_pass = (filter_raw == "PASS")

    # Strand counts
    fau = _get_strand_count(info, "FAU")
    fcu = _get_strand_count(info, "FCU")
    fgu = _get_strand_count(info, "FGU")
    ftu = _get_strand_count(info, "FTU")
    rau = _get_strand_count(info, "RAU")
    rcu = _get_strand_count(info, "RCU")
    rgu = _get_strand_count(info, "RGU")
    rtu = _get_strand_count(info, "RTU")
    sb_clairs = _to_float(info.get("SB"))

    forward_ref, reverse_ref, forward_alt, reverse_alt = _get_ref_alt_strand_counts(
        ref_allele, alt_allele, fau, fcu, fgu, ftu, rau, rcu, rgu, rtu,
    )
    sb_ratio = compute_sb_ratio(forward_alt, reverse_alt)
    sb_fisher_p = compute_sb_fisher(forward_ref, reverse_ref, forward_alt, reverse_alt)

    # Trinucleotide context
    trinuc = get_trinuc_context(fasta, chrom, pos)
    mut_6class = MUTATION_6CLASS.get((ref_allele, alt_allele), "unknown")
    is_ti = mut_6class in TRANSITIONS
    cpg_site = is_cpg_site_check(trinuc, ref_allele)
    cpg_ct = is_cpg_ct_check(ref_allele, alt_allele, trinuc)

    return {
        "chrom": chrom,
        "pos": pos,
        "ref": ref_allele,
        "alt": alt_allele,
        "variant_key": variant_key(chrom, pos, ref_allele, alt_allele),
        "qual": _to_float(rec.qual, 0.0),
        "gq": _to_int(sample_data.get("GQ")),
        "dp": _to_int(sample_data.get("DP")),
        "af": _to_float(sample_data.get("AF")),
        "ad_ref": ad_ref,
        "ad_alt": ad_alt,
        "gt_field": gt_str,
        "is_het": is_het,
        "is_hom_alt": is_hom_alt,
        "h_flag": h_flag,
        "pon_1": pon_1,
        "pon_2": pon_2,
        "pon_3": pon_3,
        "pon_4": pon_4,
        "pon_count": pon_count,
        "verdict_germline": verdict_germline,
        "verdict_somatic": verdict_somatic,
        "verdict_subclonal": verdict_subclonal,
        "filter_raw": filter_raw,
        "filter_is_pass": filter_is_pass,
        "fau": fau, "fcu": fcu, "fgu": fgu, "ftu": ftu,
        "rau": rau, "rcu": rcu, "rgu": rgu, "rtu": rtu,
        "sb_clairs": sb_clairs,
        "forward_ref": forward_ref, "reverse_ref": reverse_ref,
        "forward_alt": forward_alt, "reverse_alt": reverse_alt,
        "sb_ratio": sb_ratio,
        "sb_fisher_p": sb_fisher_p,
        "trinuc_context": trinuc,
        "mutation_type_6class": mut_6class,
        "is_transition": is_ti,
        "is_cpg_site": cpg_site,
        "is_cpg_ct": cpg_ct,
    }


# ---------------------------------------------------------------------------
# Per-sample worker
# ---------------------------------------------------------------------------


def extract_features_lightweight(rec: pysam.VariantRecord) -> Dict:
    """Lightweight feature extraction for Tier-2 records (skip fisher + trinuc FASTA)."""
    chrom = rec.chrom
    pos = rec.pos
    ref_allele = rec.ref
    alt_allele = rec.alts[0] if rec.alts else "."

    sample_data = list(rec.samples.values())[0] if rec.samples else {}
    gt_tuple = sample_data.get("GT", (None, None))
    if gt_tuple and gt_tuple != (None, None):
        gt_str = "/".join(str(a) if a is not None else "." for a in gt_tuple)
        is_het = (gt_tuple[0] != gt_tuple[1]) if len(gt_tuple) >= 2 else False
        is_hom_alt = (gt_tuple[0] == gt_tuple[1] and gt_tuple[0] is not None and gt_tuple[0] > 0)
    else:
        gt_str = "./."
        is_het = False
        is_hom_alt = False

    ad = sample_data.get("AD")
    if ad is not None and isinstance(ad, (list, tuple)) and len(ad) >= 2:
        ad_ref, ad_alt = int(ad[0]), int(ad[1])
    else:
        ad_ref, ad_alt = 0, 0

    info = rec.info
    filter_terms = list(rec.filter.keys()) or ["PASS"]
    filter_raw = ";".join(filter_terms)

    fau = _get_strand_count(info, "FAU")
    fcu = _get_strand_count(info, "FCU")
    fgu = _get_strand_count(info, "FGU")
    ftu = _get_strand_count(info, "FTU")
    rau = _get_strand_count(info, "RAU")
    rcu = _get_strand_count(info, "RCU")
    rgu = _get_strand_count(info, "RGU")
    rtu = _get_strand_count(info, "RTU")

    forward_ref, reverse_ref, forward_alt, reverse_alt = _get_ref_alt_strand_counts(
        ref_allele, alt_allele, fau, fcu, fgu, ftu, rau, rcu, rgu, rtu,
    )

    mut_6class = MUTATION_6CLASS.get((ref_allele, alt_allele), "unknown")

    return {
        "chrom": chrom, "pos": pos, "ref": ref_allele, "alt": alt_allele,
        "variant_key": variant_key(chrom, pos, ref_allele, alt_allele),
        "qual": _to_float(rec.qual, 0.0),
        "gq": _to_int(sample_data.get("GQ")),
        "dp": _to_int(sample_data.get("DP")),
        "af": _to_float(sample_data.get("AF")),
        "ad_ref": ad_ref, "ad_alt": ad_alt,
        "gt_field": gt_str, "is_het": is_het, "is_hom_alt": is_hom_alt,
        "h_flag": "H" in info,
        "pon_1": "PoN_1" in info, "pon_2": "PoN_2" in info,
        "pon_3": "PoN_3" in info, "pon_4": "PoN_4" in info,
        "pon_count": sum(1 for i in range(1, 5) if f"PoN_{i}" in info),
        "verdict_germline": "Verdict_Germline" in info,
        "verdict_somatic": "Verdict_Somatic" in info,
        "verdict_subclonal": "Verdict_SubclonalSomatic" in info,
        "filter_raw": filter_raw,
        "filter_is_pass": (filter_raw == "PASS"),
        "fau": fau, "fcu": fcu, "fgu": fgu, "ftu": ftu,
        "rau": rau, "rcu": rcu, "rgu": rgu, "rtu": rtu,
        "sb_clairs": _to_float(info.get("SB")),
        "forward_ref": forward_ref, "reverse_ref": reverse_ref,
        "forward_alt": forward_alt, "reverse_alt": reverse_alt,
        "sb_ratio": compute_sb_ratio(forward_alt, reverse_alt),
        "sb_fisher_p": float("nan"),  # Skip fisher for Tier-2
        "trinuc_context": "", "mutation_type_6class": mut_6class,
        "is_transition": mut_6class in TRANSITIONS,
        "is_cpg_site": False, "is_cpg_ct": False,  # Skip trinuc for Tier-2
    }


def _fast_scan_vcf_keys(vcf_path: str) -> Tuple[Dict[str, dict], int, int]:
    """Fast text-based scan: build dict of {variant_key -> lightweight info} for all SNV records.

    Returns (key_info_dict, n_snv_total, n_pass).
    key_info has: filter_terms, has_pon, pon_1..4.
    """
    import gzip as _gzip

    opener = _gzip.open if vcf_path.endswith(".gz") else open
    key_info: Dict[str, dict] = {}
    n_snv = 0
    n_pass = 0

    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t", 10)
            if len(fields) < 8:
                continue
            ref = fields[3]
            alt = fields[4]
            # SNV only
            if len(ref) != 1 or len(alt) != 1 or alt == ".":
                continue
            n_snv += 1
            filt = fields[6]
            info_str = fields[7]
            is_pass = (filt == "PASS" or filt == ".")
            if is_pass:
                n_pass += 1

            vkey = f"{fields[0]}:{fields[1]}:{ref}:{alt}"
            pon1 = "PoN_1" in info_str
            pon2 = "PoN_2" in info_str
            pon3 = "PoN_3" in info_str
            pon4 = "PoN_4" in info_str

            # Extract AF from FORMAT field for PoN summary
            af_val = float("nan")
            if len(fields) >= 10:
                fmt_keys = fields[8].split(":")
                fmt_vals = fields[9].split(":")
                if "AF" in fmt_keys:
                    idx = fmt_keys.index("AF")
                    if idx < len(fmt_vals):
                        try:
                            af_val = float(fmt_vals[idx])
                        except (ValueError, TypeError):
                            pass

            # Extract GT from FORMAT
            gt_str = "./."
            if len(fields) >= 10:
                fmt_keys = fields[8].split(":")
                fmt_vals = fields[9].split(":")
                if "GT" in fmt_keys:
                    idx = fmt_keys.index("GT")
                    if idx < len(fmt_vals):
                        gt_str = fmt_vals[idx]

            # H flag
            h_flag = ";H;" in f";{info_str};" or info_str.startswith("H;") or info_str.endswith(";H") or info_str == "H"

            key_info[vkey] = {
                "is_pass": is_pass,
                "pon_1": pon1, "pon_2": pon2, "pon_3": pon3, "pon_4": pon4,
                "pon_count": int(pon1) + int(pon2) + int(pon3) + int(pon4),
                "filter_raw": filt,
                "af": af_val,
                "gt": gt_str,
                "h_flag": h_flag,
            }

    return key_info, n_snv, n_pass


def process_sample(sample_name: str) -> Tuple[str, int, int, int]:
    """Process one sample: load truth, scan TO VCF, extract features, write output.

    Strategy:
      1. Fast text scan to build key->filter lookup for ALL SNVs
      2. Match truth keys to identify TP
      3. Select Tier-1 records: PASS TP + PASS FP (full feature extraction)
      4. Select Tier-2 records: non-PASS TP + PoN-filtered FP (lightweight: no fisher, no trinuc FASTA)
      5. Use pysam fetch-by-position only for Tier-1 + Tier-2 selected records
    """
    spec = SAMPLE_SPECS[sample_name]
    to_vcf_path = spec.to_somatic_vcf_pileup
    truth_vcf_path = spec.truth_vcf
    truth_bed_path = spec.truth_bed

    _print(f"\n{'='*60}")
    _print(f"[G1] Processing {sample_name}")
    _print(f"  TO VCF:    {to_vcf_path}")

    # Step 1: Load truth keys
    tp_keys = load_truth_keys_from_vcf(truth_vcf_path, truth_bed_path)
    _print(f"  [truth] {len(tp_keys):,} TP keys loaded")

    # Step 2: Fast text scan for all SNV keys
    _print(f"  [scan] Fast text scan of {Path(to_vcf_path).name} ...")
    all_keys, n_snv_total, n_pass = _fast_scan_vcf_keys(to_vcf_path)
    _print(f"  [scan] {n_snv_total:,} SNVs total, {n_pass:,} PASS, {len(all_keys):,} unique keys")

    # Step 3: Classify — only PASS variants need full pysam extraction
    # PoN-filtered FP (~millions) handled by text-scan summary only
    pass_keys: Set[str] = set()
    for vkey, info in all_keys.items():
        if info["is_pass"]:
            pass_keys.add(vkey)
        elif vkey in tp_keys:
            pass_keys.add(vkey)  # non-PASS TP also needed for completeness

    _print(f"  [classify] PASS + non-PASS-TP keys for pysam: {len(pass_keys):,}")

    # Step 4: pysam extraction (full features) — only for selected keys
    _print(f"  [extract] Opening pysam for full feature extraction ...")

    fasta = pysam.FastaFile(str(REF_FASTA))
    vf = pysam.VariantFile(to_vcf_path)
    rows: List[Dict] = []
    n_extracted = 0

    for rec in vf.fetch():
        if len(rec.ref) != 1:
            continue
        if not rec.alts or len(rec.alts) != 1 or len(rec.alts[0]) != 1:
            continue

        vkey = variant_key(rec.chrom, rec.pos, rec.ref, rec.alts[0])
        if vkey not in pass_keys:
            continue

        is_tp = vkey in tp_keys
        features = extract_features_from_record(rec, fasta)
        features["sample"] = sample_name
        features["truth_label"] = "TP" if is_tp else "FP"
        rows.append(features)
        n_extracted += 1

        if n_extracted % 10000 == 0:
            _print(f"    ... extracted {n_extracted:,} / {len(pass_keys):,}")

    vf.close()
    fasta.close()

    # Step 5: Write output
    out_path = G1_DIR / f"{sample_name}_features.tsv.gz"
    _write_tsv_gz(rows, out_path)

    n_tp = sum(1 for r in rows if r["truth_label"] == "TP")
    n_fp = sum(1 for r in rows if r["truth_label"] == "FP")
    n_pass_tp = sum(1 for r in rows if r["truth_label"] == "TP" and r.get("filter_is_pass"))
    n_pass_fp = sum(1 for r in rows if r["truth_label"] == "FP" and r.get("filter_is_pass"))
    _print(f"  [G1] {sample_name}: {n_snv_total:,} SNVs -> {len(rows):,} kept")
    _print(f"        TP={n_tp:,} (PASS={n_pass_tp:,}), FP={n_fp:,} (PASS={n_pass_fp:,})")

    # Write PoN summary for G2
    _write_pon_summary(sample_name, all_keys, tp_keys)

    return sample_name, len(rows), n_tp, n_fp


def _write_pon_summary(sample_name: str, all_keys: Dict, tp_keys: Set[str]) -> None:
    """Write PoN breakdown + AF/GT summary for all variants (from text scan, no pysam)."""
    import numpy as np

    summary_path = G1_DIR / f"{sample_name}_pon_summary.tsv"
    detail_path = G1_DIR / f"{sample_name}_all_variants_lite.tsv.gz"

    # Bucket counts
    counts = {"TP": {}, "FP": {}}
    for vkey, info in all_keys.items():
        label = "TP" if vkey in tp_keys else "FP"
        pc = info["pon_count"]
        is_pass = info["is_pass"]
        bucket = f"pon{pc}_{'pass' if is_pass else 'filtered'}"
        counts[label][bucket] = counts[label].get(bucket, 0) + 1

    with open(summary_path, "w") as f:
        f.write("sample\ttruth_label\tbucket\tcount\n")
        for label in ["TP", "FP"]:
            for bucket, count in sorted(counts[label].items()):
                f.write(f"{sample_name}\t{label}\t{bucket}\t{count}\n")

    # Lightweight detail: variant_key, truth_label, is_pass, pon_count, af, gt, h_flag
    with gzip.open(detail_path, "wt") as f:
        f.write("variant_key\ttruth_label\tis_pass\tpon_count\taf\tgt\th_flag\n")
        for vkey, info in all_keys.items():
            label = "TP" if vkey in tp_keys else "FP"
            af_str = f"{info['af']:.4f}" if not math.isnan(info['af']) else ""
            f.write(f"{vkey}\t{label}\t{1 if info['is_pass'] else 0}\t{info['pon_count']}\t{af_str}\t{info['gt']}\t{1 if info['h_flag'] else 0}\n")

    _print(f"  [PoN] Written PoN summary to {summary_path.name}")
    _print(f"  [PoN] Written lightweight detail to {detail_path.name} ({detail_path.stat().st_size / 1024 / 1024:.1f} MB)")


def _write_tsv_gz(rows: List[Dict], path: Path) -> None:
    if not rows:
        return
    ensure_dir(path.parent)
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write("\t".join(G1_COLUMNS) + "\n")
        for row in rows:
            vals = []
            for col in G1_COLUMNS:
                v = row.get(col, "")
                if isinstance(v, bool):
                    vals.append("1" if v else "0")
                elif isinstance(v, float):
                    if math.isnan(v):
                        vals.append("")
                    else:
                        vals.append(f"{v:.6g}")
                else:
                    vals.append(str(v))
            f.write("\t".join(vals) + "\n")
    print(f"  [G1] Written {len(rows):,} rows to {path.name} ({path.stat().st_size / 1024 / 1024:.1f} MB)")


# ---------------------------------------------------------------------------
# Merge all samples
# ---------------------------------------------------------------------------


def merge_all_samples() -> int:
    """Merge per-sample TSVs into one combined file."""
    merged_path = G1_DIR / "all_samples_merged.tsv.gz"
    total = 0

    with gzip.open(merged_path, "wt", encoding="utf-8") as out:
        header_written = False
        for sample in SAMPLE_ORDER:
            sample_path = G1_DIR / f"{sample}_features.tsv.gz"
            if not sample_path.exists():
                print(f"  [merge] WARNING: {sample_path.name} not found, skipping")
                continue
            with gzip.open(sample_path, "rt", encoding="utf-8") as inp:
                for i, line in enumerate(inp):
                    if i == 0:
                        if not header_written:
                            out.write(line)
                            header_written = True
                        continue
                    out.write(line)
                    total += 1

    print(f"\n[G1] Merged {total:,} rows -> {merged_path.name} ({merged_path.stat().st_size / 1024 / 1024:.1f} MB)")
    return total


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    _print("=" * 70)
    _print("G1: VCF Feature Extraction for Germline FP Identification")
    _print(f"  Samples: {', '.join(SAMPLE_ORDER)}")
    _print(f"  Reference: {REF_FASTA}")
    _print(f"  Output: {G1_DIR}")
    _print("=" * 70)

    results = []
    for sample in SAMPLE_ORDER:
        if sample not in SAMPLE_SPECS:
            _print(f"  [SKIP] {sample} not in SAMPLE_SPECS")
            continue
        if not SAMPLE_SPECS[sample].to_somatic_vcf_pileup:
            _print(f"  [SKIP] {sample} has no TO VCF")
            continue
        result = process_sample(sample)
        results.append(result)

    # Summary
    _print("\n" + "=" * 70)
    _print("G1 Per-Sample Summary:")
    _print(f"{'Sample':<20} {'Total':>8} {'TP':>8} {'FP':>8}")
    _print("-" * 50)
    for name, total, tp, fp in results:
        _print(f"{name:<20} {total:>8,} {tp:>8,} {fp:>8,}")

    # Merge
    total_merged = merge_all_samples()

    _print(f"\n[G1] DONE. Total merged rows: {total_merged:,}")
    _print(f"[G1] Output: {G1_DIR}")


if __name__ == "__main__":
    main()
