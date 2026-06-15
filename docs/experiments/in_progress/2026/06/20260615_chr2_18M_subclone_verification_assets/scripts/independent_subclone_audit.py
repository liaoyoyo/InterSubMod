#!/usr/bin/env python3
"""
Independent audit of the HCC1395 chr2:18.07-18.11 Mb regional subclone claim.

The audit intentionally separates:
  1. observed read-level variant linkage,
  2. methylation association with a local allele/lineage anchor,
  3. matched-normal haplotype methylation,
  4. inferred lineage topology and evolutionary order.

It uses exact GRCh38 CpG cytosine coordinates. In particular, the figure's
"4.1" label is chr2:18,096,341 (G coordinate), whose CpG C coordinate is
chr2:18,096,340. The previously transcribed chr2:18,096,041 is not a CpG.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam
from scipy.stats import fisher_exact, mannwhitneyu


CHROM = "chr2"
REGION_START0 = 18_066_480
REGION_END0 = 18_110_828
REFERENCE = "/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
SEQC2_TRUTH = (
    "/big8_disk/data/HCC1395/SEQC2/"
    "high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
)
SEQC2_HC_BED = "/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"
SEQC2_LOH_BED = (
    "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
)

SNVS = {
    "1": {"pos": 18_068_480, "ref": "C", "alt": "G"},
    "2": {"pos": 18_072_546, "ref": "G", "alt": "C"},
    "3": {"pos": 18_086_020, "ref": "G", "alt": "A"},
    "4": {"pos": 18_096_269, "ref": "C", "alt": "G/T/DEL"},
    "5": {"pos": 18_099_697, "ref": "G", "alt": "C"},
    "6": {"pos": 18_108_828, "ref": "C", "alt": "G"},
}

# Exact 1-based CpG cytosine coordinates after checking GRCh38.
CPGS = {
    "2.1": 18_073_987,
    "2.2": 18_076_408,  # user's 18,076,409 is the CpG G coordinate
    "3.1": 18_090_947,
    "3.2": 18_091_583,  # user's 18,091,584 is the CpG G coordinate
    "3.3": 18_093_413,  # user's 18,093,414 is the CpG G coordinate
    "3.4": 18_094_113,  # user's 18,094,114 is the CpG G coordinate
    "3.5": 18_096_032,  # user's 18,096,033 is the CpG G coordinate
    "4.1": 18_096_340,  # figure: 18,096,341 G coordinate; text had a typo
    "5.1": 18_100_543,  # user's 18,100,544 is the CpG G coordinate
    "5.2": 18_100_682,  # user's 18,100,683 is the CpG G coordinate
}

SAMPLES = {
    "HKU_T": {
        "role": "tumor",
        "bam": "/big8_disk/liaoyoyo2001/data/bam/"
        "HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam",
    },
    "HKU_N": {
        "role": "normal",
        "bam": "/big8_disk/liaoyoyo2001/data/bam/"
        "HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam",
    },
    "DORADO_TAGGED_T": {
        "role": "tumor",
        "bam": "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395_DORADO/"
        "paired_full/20260315_HCC1395_DORADO_paired_full_full_complete_matrix/"
        "longphase_s/HCC1395_DORADO_tagged.bam",
    },
    "DORADO_RAW_T": {
        "role": "tumor",
        "bam": "/big8_disk/data/HCC1395/ONT_Dorado/HCC1395.bam",
    },
    "DORADO_N": {
        "role": "normal",
        "bam": "/big8_disk/data/HCC1395/ONT_Dorado/HCC1395BL.bam",
    },
}


def bed_intervals(path: str, chrom: str) -> list[dict]:
    intervals = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[0] != chrom:
                continue
            intervals.append(
                {
                    "chrom": fields[0],
                    "start0": int(fields[1]),
                    "end0": int(fields[2]),
                    "extra": fields[3:],
                }
            )
    return intervals


def benchmark_context() -> dict:
    hc_intervals = bed_intervals(SEQC2_HC_BED, CHROM)
    loh_intervals = bed_intervals(SEQC2_LOH_BED, CHROM)
    truth_by_pos = defaultdict(list)
    truth_vcf = pysam.VariantFile(SEQC2_TRUTH)
    for record in truth_vcf.fetch(CHROM, REGION_START0, REGION_END0):
        truth_by_pos[record.pos].append(
            {
                "ref": record.ref,
                "alts": list(record.alts or []),
                "filters": list(record.filter.keys()),
            }
        )
    truth_vcf.close()

    site_status = {}
    for label, site in SNVS.items():
        pos0 = site["pos"] - 1
        containing_hc = [
            interval
            for interval in hc_intervals
            if interval["start0"] <= pos0 < interval["end0"]
        ]
        truth_records = truth_by_pos.get(site["pos"], [])
        if truth_records:
            status = "truth_confirmed"
        elif not containing_hc:
            status = "out_of_hc_truth_unevaluable"
        else:
            status = "in_hc_not_in_truth"
        site_status[label] = {
            "position": site["pos"],
            "status": status,
            "in_hc": bool(containing_hc),
            "hc_intervals": containing_hc,
            "truth_records": truth_records,
        }

    overlapping_loh = [
        interval
        for interval in loh_intervals
        if interval["start0"] < REGION_END0 and interval["end0"] > REGION_START0
    ]
    return {
        "truth_vcf": SEQC2_TRUTH,
        "hc_bed": SEQC2_HC_BED,
        "loh_bed": SEQC2_LOH_BED,
        "site_status": site_status,
        "overlapping_loh_intervals": overlapping_loh,
    }


def bh_fdr(pvalues: list[float | None]) -> list[float | None]:
    valid = [(i, p) for i, p in enumerate(pvalues) if p is not None and math.isfinite(p)]
    valid.sort(key=lambda x: x[1])
    out: list[float | None] = [None] * len(pvalues)
    running = 1.0
    m = len(valid)
    for rank_from_end, (idx, p) in enumerate(reversed(valid), start=1):
        rank = m - rank_from_end + 1
        running = min(running, p * m / rank)
        out[idx] = min(1.0, running)
    return out


def hp_parent(hp) -> str | None:
    if hp is None:
        return None
    value = str(hp)
    if value.startswith("1"):
        return "1"
    if value.startswith("2"):
        return "2"
    return value


def deduplicate(records: list[dict]) -> tuple[list[dict], dict]:
    by_name: dict[str, list[dict]] = defaultdict(list)
    for record in records:
        by_name[record["name"]].append(record)

    kept = []
    inconsistent = 0
    for name, group in by_name.items():
        signatures = {
            (
                r["ref_start"],
                r["ref_end"],
                r["strand"],
                json.dumps(r["bases"], sort_keys=True),
                str(r["HP"]),
                str(r["PS"]),
            )
            for r in group
        }
        if len(signatures) > 1:
            inconsistent += 1
        group.sort(
            key=lambda r: (
                r["mapq"],
                sum(v["base"] is not None for v in r["bases"].values()),
                r["ref_end"] - r["ref_start"],
            ),
            reverse=True,
        )
        kept.append(group[0])

    kept.sort(key=lambda r: (r["ref_start"], r["name"]))
    stats = {
        "records": len(records),
        "unique_names": len(kept),
        "duplicate_extra_records": len(records) - len(kept),
        "duplicate_names_with_inconsistent_records": inconsistent,
    }
    return kept, stats


def extract_records(bam_path: str, fasta: pysam.FastaFile) -> tuple[list[dict], dict]:
    records = []
    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch(CHROM, REGION_START0, REGION_END0):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        ref_to_query = {}
        query_to_ref = {}
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
            if ref_pos is not None:
                ref_to_query[ref_pos] = query_pos
            if query_pos is not None and ref_pos is not None:
                query_to_ref[query_pos] = ref_pos

        bases = {}
        for label, site in SNVS.items():
            pos0 = site["pos"] - 1
            query_pos = ref_to_query.get(pos0, "NOT_COVERED")
            if query_pos == "NOT_COVERED":
                base = None
                baseq = None
            elif query_pos is None:
                base = "DEL"
                baseq = None
            else:
                base = read.query_sequence[query_pos].upper()
                baseq = int(read.query_qualities[query_pos]) if read.query_qualities is not None else None
            bases[label] = {"base": base, "baseq": baseq}

        methylation = {}
        normalized_mods: dict[int, dict[str, float]] = defaultdict(dict)
        for (_canon, _mod_strand, mod), calls in (read.modified_bases or {}).items():
            mod_name = chr(mod) if isinstance(mod, int) else str(mod)
            for query_pos, qual in calls:
                ref_pos0 = query_to_ref.get(query_pos)
                if ref_pos0 is None:
                    continue
                base = fasta.fetch(CHROM, ref_pos0, ref_pos0 + 1).upper()
                prev_base = fasta.fetch(CHROM, ref_pos0 - 1, ref_pos0).upper()
                next_base = fasta.fetch(CHROM, ref_pos0 + 1, ref_pos0 + 2).upper()
                if base == "C" and next_base == "G":
                    c_pos1 = ref_pos0 + 1
                elif base == "G" and prev_base == "C":
                    c_pos1 = ref_pos0
                else:
                    continue
                normalized_mods[c_pos1][mod_name] = qual / 255.0

        for label, c_pos1 in CPGS.items():
            mods = normalized_mods.get(c_pos1)
            methylation[label] = None if not mods else {
                "m": mods.get("m"),
                "h": mods.get("h"),
            }

        records.append(
            {
                "name": read.query_name,
                "mapq": int(read.mapping_quality),
                "strand": "-" if read.is_reverse else "+",
                "ref_start": int(read.reference_start),
                "ref_end": int(read.reference_end),
                "HP": read.get_tag("HP") if read.has_tag("HP") else None,
                "PS": read.get_tag("PS") if read.has_tag("PS") else None,
                "bases": bases,
                "meth": methylation,
            }
        )
    bam.close()
    return deduplicate(records)


def callable_base(record: dict, site: str, min_mapq=20, min_baseq=10) -> str | None:
    if record["mapq"] < min_mapq:
        return None
    call = record["bases"][site]
    base = call["base"]
    if base is None:
        return None
    if base == "DEL":
        return base
    if call["baseq"] is None or call["baseq"] < min_baseq:
        return None
    return base


def is_event(record: dict, event: str, min_baseq=10) -> bool | None:
    site = event[1]
    base = callable_base(record, site, min_baseq=min_baseq)
    if base is None:
        return None
    if event == "E1":
        return base == "G"
    if event == "E2":
        return base == "C"
    if event == "E3":
        return base == "A"
    if event == "E4G":
        return base == "G"
    if event == "E4T":
        return base == "T"
    if event == "E4D":
        return base == "DEL"
    if event == "E4ANY":
        return base != SNVS["4"]["ref"]
    if event == "E5":
        return base == "C"
    if event == "E6":
        return base == "G"
    raise ValueError(event)


def event_table(records: list[dict], event_a: str, event_b: str, min_baseq=10) -> dict:
    table = Counter()
    for record in records:
        a = is_event(record, event_a, min_baseq)
        b = is_event(record, event_b, min_baseq)
        if a is None or b is None:
            continue
        table[f"{int(a)}{int(b)}"] += 1
    return {
        "n": sum(table.values()),
        "00": table["00"],
        "10": table["10"],
        "01": table["01"],
        "11": table["11"],
    }


def allele_counts(records: list[dict], min_baseq: int) -> dict:
    result = {}
    for site in SNVS:
        counts = Counter()
        strands = defaultdict(Counter)
        hp_by_base = defaultdict(Counter)
        for record in records:
            base = callable_base(record, site, min_baseq=min_baseq)
            if base is not None:
                counts[base] += 1
                strands[base][record["strand"]] += 1
                hp_by_base[base][str(record["HP"])] += 1
        result[site] = {
            "counts": dict(counts),
            "strands": {base: dict(value) for base, value in strands.items()},
            "hp_by_base": {base: dict(value) for base, value in hp_by_base.items()},
        }
    return result


def group_alpha(record: dict) -> bool:
    return is_event(record, "E3") is True


def group_pos3_ref(record: dict) -> bool:
    return callable_base(record, "3") == "G"


def group_beta_strict(record: dict) -> bool:
    return any(is_event(record, event) is True for event in ("E1", "E2", "E6"))


def group_normal_hp1(record: dict) -> bool:
    return hp_parent(record["HP"]) == "1"


def group_normal_hp2(record: dict) -> bool:
    return hp_parent(record["HP"]) == "2"


def variant_group_functions(site: str):
    def alt_group(record: dict) -> bool:
        base = callable_base(record, site)
        if base is None:
            return False
        if site == "4":
            return base != SNVS[site]["ref"]
        return base == SNVS[site]["alt"]

    def ref_group(record: dict) -> bool:
        return callable_base(record, site) == SNVS[site]["ref"]

    return alt_group, ref_group


def group_methyl_test(records: list[dict], group_a, group_b) -> dict:
    rows = []
    for label in CPGS:
        a = [
            record["meth"][label]["m"]
            for record in records
            if group_a(record)
            and record["meth"][label] is not None
            and record["meth"][label].get("m") is not None
        ]
        b = [
            record["meth"][label]["m"]
            for record in records
            if group_b(record)
            and record["meth"][label] is not None
            and record["meth"][label].get("m") is not None
        ]
        if len(a) >= 2 and len(b) >= 2:
            mw_p = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
            ah = sum(value >= 0.5 for value in a)
            bh = sum(value >= 0.5 for value in b)
            fisher_p = float(
                fisher_exact([[ah, len(a) - ah], [bh, len(b) - bh]], alternative="two-sided").pvalue
            )
        else:
            mw_p = None
            fisher_p = None
        rows.append(
            {
                "site": label,
                "c_pos1": CPGS[label],
                "a_n": len(a),
                "a_mean": float(np.mean(a)) if a else None,
                "a_high_fraction": sum(value >= 0.5 for value in a) / len(a) if a else None,
                "b_n": len(b),
                "b_mean": float(np.mean(b)) if b else None,
                "b_high_fraction": sum(value >= 0.5 for value in b) / len(b) if b else None,
                "mean_delta_a_minus_b": float(np.mean(a) - np.mean(b)) if a and b else None,
                "mann_whitney_p": mw_p,
                "fisher_binary_p": fisher_p,
            }
        )
    mw_fdr = bh_fdr([row["mann_whitney_p"] for row in rows])
    fisher_fdr = bh_fdr([row["fisher_binary_p"] for row in rows])
    for row, mw_q, fisher_q in zip(rows, mw_fdr, fisher_fdr):
        row["mann_whitney_fdr"] = mw_q
        row["fisher_binary_fdr"] = fisher_q
    return {"rows": rows}


def two_set_methyl_test(records_a: list[dict], records_b: list[dict]) -> dict:
    rows = []
    for label in CPGS:
        values = []
        for records in (records_a, records_b):
            values.append(
                [
                    record["meth"][label]["m"]
                    for record in records
                    if record["meth"][label] is not None
                    and record["meth"][label].get("m") is not None
                ]
            )
        a, b = values
        if len(a) >= 2 and len(b) >= 2:
            mw_p = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
            ah = sum(value >= 0.5 for value in a)
            bh = sum(value >= 0.5 for value in b)
            fisher_p = float(
                fisher_exact([[ah, len(a) - ah], [bh, len(b) - bh]], alternative="two-sided").pvalue
            )
        else:
            mw_p = None
            fisher_p = None
        rows.append(
            {
                "site": label,
                "c_pos1": CPGS[label],
                "a_n": len(a),
                "a_mean": float(np.mean(a)) if a else None,
                "a_high_fraction": sum(value >= 0.5 for value in a) / len(a) if a else None,
                "b_n": len(b),
                "b_mean": float(np.mean(b)) if b else None,
                "b_high_fraction": sum(value >= 0.5 for value in b) / len(b) if b else None,
                "mean_delta_a_minus_b": float(np.mean(a) - np.mean(b)) if a and b else None,
                "mann_whitney_p": mw_p,
                "fisher_binary_p": fisher_p,
            }
        )
    mw_fdr = bh_fdr([row["mann_whitney_p"] for row in rows])
    fisher_fdr = bh_fdr([row["fisher_binary_p"] for row in rows])
    for row, mw_q, fisher_q in zip(rows, mw_fdr, fisher_fdr):
        row["mann_whitney_fdr"] = mw_q
        row["fisher_binary_fdr"] = fisher_q
    return {"rows": rows}


def all_ref_coverage_check(records: list[dict]) -> dict:
    result = {}
    for minimum_coverage in (2, 3, 4, 5, 6):
        checked = 0
        all_ref = 0
        for record in records:
            calls = {
                site: callable_base(record, site)
                for site in SNVS
            }
            covered = {site: base for site, base in calls.items() if base is not None}
            if len(covered) < minimum_coverage:
                continue
            checked += 1
            if all(base == SNVS[site]["ref"] for site, base in covered.items()):
                all_ref += 1
        result[str(minimum_coverage)] = {"checked": checked, "all_ref": all_ref}
    return result


def pos4_subtype_summary(records: list[dict]) -> dict:
    result = {}
    for subtype in ("G", "T", "DEL"):
        selected = [record for record in records if callable_base(record, "4") == subtype]
        methyl = {}
        for cpg in CPGS:
            values = [
                record["meth"][cpg]["m"]
                for record in selected
                if record["meth"][cpg] is not None and record["meth"][cpg].get("m") is not None
            ]
            methyl[cpg] = {
                "n": len(values),
                "mean": float(np.mean(values)) if values else None,
                "high_fraction": sum(value >= 0.5 for value in values) / len(values) if values else None,
            }
        result[subtype] = {
            "n": len(selected),
            "strand_counts": dict(Counter(record["strand"] for record in selected)),
            "with_E6_G": sum(is_event(record, "E6") is True for record in selected),
            "methyl": methyl,
        }
    return result


def sample_summary(records: list[dict], role: str, dedup_stats: dict) -> dict:
    hp = Counter(str(record["HP"]) for record in records)
    hp_parent_counts = Counter(str(hp_parent(record["HP"])) for record in records)
    summary = {
        "role": role,
        "dedup": dedup_stats,
        "hp_counts": dict(hp),
        "hp_parent_counts": dict(hp_parent_counts),
        "allele_counts_bq10": allele_counts(records, 10),
        "allele_counts_bq20": allele_counts(records, 20),
    }
    if role == "tumor":
        pairs = [
            ("E1", "E2"),
            ("E1", "E3"),
            ("E2", "E3"),
            ("E1", "E6"),
            ("E2", "E6"),
            ("E3", "E5"),
            ("E3", "E4ANY"),
            ("E3", "E6"),
            ("E4ANY", "E5"),
            ("E4ANY", "E6"),
            ("E4G", "E6"),
            ("E4T", "E6"),
            ("E4D", "E6"),
            ("E5", "E6"),
        ]
        summary["linkage_bq10"] = {
            f"{a}_vs_{b}": event_table(records, a, b, 10) for a, b in pairs
        }
        summary["lineage_anchor_counts"] = {
            "alpha_pos3_A": sum(group_alpha(record) for record in records),
            "pos3_ref_G": sum(group_pos3_ref(record) for record in records),
            "beta_strict_E1_or_E2_or_E6": sum(group_beta_strict(record) for record in records),
            "alpha_and_beta_strict_violation": sum(
                group_alpha(record) and group_beta_strict(record) for record in records
            ),
        }
        summary["all_ref_coverage_check"] = all_ref_coverage_check(records)
        summary["pos4_subtype_summary"] = pos4_subtype_summary(records)
        summary["methyl_alpha_vs_pos3_ref"] = group_methyl_test(
            records, group_alpha, group_pos3_ref
        )
        summary["methyl_alpha_vs_beta_strict"] = group_methyl_test(
            records, group_alpha, group_beta_strict
        )
        summary["methyl_by_variant_anchor"] = {}
        for site in SNVS:
            alt_group, ref_group = variant_group_functions(site)
            summary["methyl_by_variant_anchor"][site] = group_methyl_test(
                records, alt_group, ref_group
            )
    else:
        summary["methyl_normal_hp1_vs_hp2"] = group_methyl_test(
            records, group_normal_hp1, group_normal_hp2
        )
    return summary


def fmt_float(value, digits=3):
    return "NA" if value is None else f"{value:.{digits}g}"


def write_markdown(result: dict, output_path: Path) -> None:
    lines = [
        "# Independent chr2:18M Subclone Audit",
        "",
        "## Coordinate correction",
        "",
        "- Figure `(4.1)` is `chr2:18,096,341` (CpG G coordinate; C coordinate "
        "`chr2:18,096,340`).",
        "- Transcribed `chr2:18,096,041` is an `A`, not a CpG.",
        "",
        "## SEQC2 truth, HC, and LOH status",
        "",
        f"- Truth VCF: `{result['benchmark_context']['truth_vcf']}`",
        f"- HC BED: `{result['benchmark_context']['hc_bed']}`",
        f"- LOH BED: `{result['benchmark_context']['loh_bed']}`",
        f"- Overlapping LOH intervals: "
        f"`{result['benchmark_context']['overlapping_loh_intervals']}`",
        "",
        "| Site | Position | Status | In HC | Truth records |",
        "|---|---:|---|---|---|",
    ]
    for site, status in result["benchmark_context"]["site_status"].items():
        lines.append(
            f"| {site} | {status['position']} | `{status['status']}` | "
            f"`{status['in_hc']}` | `{status['truth_records']}` |"
        )
    lines.append("")

    for sample, summary in result["samples"].items():
        lines.extend(
            [
                f"## {sample}",
                "",
                f"- Role: `{summary['role']}`",
                f"- Records / unique / duplicate extras: "
                f"`{summary['dedup']['records']}` / `{summary['dedup']['unique_names']}` / "
                f"`{summary['dedup']['duplicate_extra_records']}`",
                f"- HP parent counts: `{summary['hp_parent_counts']}`",
                "",
                "### BQ20 allele counts",
                "",
                "| Site | Position | REF | Counts |",
                "|---|---:|---|---|",
            ]
        )
        for site, info in SNVS.items():
            counts = summary["allele_counts_bq20"][site]["counts"]
            lines.append(f"| {site} | {info['pos']} | {info['ref']} | `{counts}` |")

        if summary["role"] == "tumor":
            lines.extend(
                [
                    "",
                    "### Key linkage tables",
                    "",
                    "`00/10/01/11` indicate absence/presence of the two events among reads "
                    "that cover both sites.",
                    "",
                    "| Pair | n | 00 | 10 | 01 | 11 |",
                    "|---|---:|---:|---:|---:|---:|",
                ]
            )
            for pair, table in summary["linkage_bq10"].items():
                lines.append(
                    f"| {pair} | {table['n']} | {table['00']} | {table['10']} | "
                    f"{table['01']} | {table['11']} |"
                )
            lines.extend(
                [
                    "",
                    f"- Lineage anchors: `{summary['lineage_anchor_counts']}`",
                    f"- All-REF coverage check: `{summary['all_ref_coverage_check']}`",
                    "- Position 4 subtype counts: `"
                    + str(
                        {
                            subtype: {
                                "n": value["n"],
                                "strand_counts": value["strand_counts"],
                                "with_E6_G": value["with_E6_G"],
                            }
                            for subtype, value in summary["pos4_subtype_summary"].items()
                        }
                    )
                    + "`",
                    "",
                    "### Methylation: pos3 A vs pos3 reference G",
                    "",
                    "| CpG | A n/mean | G n/mean | delta | MW FDR | Fisher FDR |",
                    "|---|---|---|---:|---:|---:|",
                ]
            )
            rows = summary["methyl_alpha_vs_pos3_ref"]["rows"]
            for row in rows:
                lines.append(
                    f"| {row['site']} | {row['a_n']}/{fmt_float(row['a_mean'])} | "
                    f"{row['b_n']}/{fmt_float(row['b_mean'])} | "
                    f"{fmt_float(row['mean_delta_a_minus_b'])} | "
                    f"{fmt_float(row['mann_whitney_fdr'])} | "
                    f"{fmt_float(row['fisher_binary_fdr'])} |"
                )
        else:
            lines.extend(
                [
                    "",
                    "### Matched-normal methylation: HP1 vs HP2",
                    "",
                    "| CpG | HP1 n/mean | HP2 n/mean | delta | MW FDR | Fisher FDR |",
                    "|---|---|---|---:|---:|---:|",
                ]
            )
            rows = summary["methyl_normal_hp1_vs_hp2"]["rows"]
            for row in rows:
                lines.append(
                    f"| {row['site']} | {row['a_n']}/{fmt_float(row['a_mean'])} | "
                    f"{row['b_n']}/{fmt_float(row['b_mean'])} | "
                    f"{fmt_float(row['mean_delta_a_minus_b'])} | "
                    f"{fmt_float(row['mann_whitney_fdr'])} | "
                    f"{fmt_float(row['fisher_binary_fdr'])} |"
                )
        lines.append("")

    lines.extend(
        [
            "## Cross-sample methylation: tumor versus matched normal",
            "",
            "These are all-read tumor/normal comparisons. They establish sample-level "
            "differential methylation, but do not by themselves separate haplotype ASM "
            "from tumor-acquired methylation.",
            "",
        ]
    )
    for comparison, summary in result["cross_sample_methyl"].items():
        lines.extend(
            [
                f"### {comparison}",
                "",
                "| CpG | Tumor n/mean | Normal n/mean | delta | MW FDR | Fisher FDR |",
                "|---|---|---|---:|---:|---:|",
            ]
        )
        for row in summary["rows"]:
            lines.append(
                f"| {row['site']} | {row['a_n']}/{fmt_float(row['a_mean'])} | "
                f"{row['b_n']}/{fmt_float(row['b_mean'])} | "
                f"{fmt_float(row['mean_delta_a_minus_b'])} | "
                f"{fmt_float(row['mann_whitney_fdr'])} | "
                f"{fmt_float(row['fisher_binary_fdr'])} |"
            )
        lines.append("")

    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", required=True)
    parser.add_argument("--md-out", required=True)
    args = parser.parse_args()

    fasta = pysam.FastaFile(REFERENCE)
    reference_checks = {}
    for label, c_pos1 in CPGS.items():
        reference_checks[label] = {
            "c_pos1": c_pos1,
            "dinucleotide": fasta.fetch(CHROM, c_pos1 - 1, c_pos1 + 1).upper(),
        }

    result = {
        "reference": REFERENCE,
        "region": f"{CHROM}:{REGION_START0 + 1}-{REGION_END0}",
        "snvs": SNVS,
        "cpg_c_positions": CPGS,
        "cpg_reference_checks": reference_checks,
        "benchmark_context": benchmark_context(),
        "samples": {},
        "cross_sample_methyl": {},
    }
    record_store = {}
    for sample, config in SAMPLES.items():
        print(f"[extract] {sample}: {config['bam']}", flush=True)
        records, dedup_stats = extract_records(config["bam"], fasta)
        record_store[sample] = records
        result["samples"][sample] = sample_summary(records, config["role"], dedup_stats)

    result["cross_sample_methyl"]["HKU_T_vs_HKU_N_all_reads"] = two_set_methyl_test(
        record_store["HKU_T"], record_store["HKU_N"]
    )
    result["cross_sample_methyl"]["DORADO_TAGGED_T_vs_DORADO_N_all_reads"] = two_set_methyl_test(
        record_store["DORADO_TAGGED_T"], record_store["DORADO_N"]
    )

    fasta.close()
    json_path = Path(args.json_out)
    md_path = Path(args.md_out)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    write_markdown(result, md_path)
    print(f"[done] {json_path}")
    print(f"[done] {md_path}")


if __name__ == "__main__":
    main()
