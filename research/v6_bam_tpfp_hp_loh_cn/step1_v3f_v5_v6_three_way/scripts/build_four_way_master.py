#!/usr/bin/env python3
"""Step 1 — baseline/V3F/V5/V6 FOUR-way ISM integration (HCC1395 paired-pileup).

Extends build_three_way_master.py to include baseline (LongPhase-TO upstream, pre-V3F)
as the 4th BAM. Output: step1_master_four_way.tsv

Reuses identical infrastructure (per-region reads.tsv aggregation, caller_af lookup,
LOH/CN master join). Only diff is BAMS list (4 BAMs instead of 3) and output filename.
"""
from __future__ import annotations

import csv
import gzip
import os
import subprocess
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
PHASEC = REPO / "research/paired_priority_bug_audit/phaseC_genome_three_way"
STEP1 = REPO / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way"
INTERMEDIATE = STEP1 / "intermediate"

VCF_TP = Path("/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz")
VCF_FP = Path("/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz")
MASTER = REPO / "research/tpfp_loh_af_kde_discrimination/data/master.tsv.gz"

BAMS = ["baseline", "V3F", "V5", "V6"]  # 4-way: baseline added
FLAGS = ["off", "on"]
LABELS = ["tp", "fp"]
RUNS = [(f, l) for f in FLAGS for l in LABELS]

NG_VALID = {"1", "2", "1-1", "2-1", "3", "11", "21", "33"}
HP_FAMILY_KEYS = ["0", "1", "2", "1-1", "2-1", "3", "11", "21", "33", "other"]


def get_region_hp_counts(reads_tsv: str) -> tuple[str, dict[str, int], int]:
    region_dir = os.path.dirname(os.path.dirname(reads_tsv))
    region_id = os.path.basename(region_dir)
    counter: Counter[str] = Counter()
    n_reads = 0
    try:
        with open(reads_tsv, "rb") as fh:
            header = fh.readline().decode().rstrip().split("\t")
            try:
                hp_idx = header.index("hp")
            except ValueError:
                return (region_id, dict(counter), 0)
            for line in fh:
                cols = line.decode().rstrip().split("\t")
                if len(cols) > hp_idx:
                    counter[cols[hp_idx]] += 1
                    n_reads += 1
    except Exception:
        pass
    return (region_id, dict(counter), n_reads)


def find_reads_tsv(run_dir: Path) -> list[str]:
    out = []
    for root, _, files in os.walk(run_dir):
        if "reads.tsv" in files:
            out.append(os.path.join(root, "reads.tsv"))
    return out


def aggregate_run(run_dir: Path, n_workers: int = 16) -> dict[str, dict]:
    files = find_reads_tsv(run_dir)
    out: dict[str, dict] = {}
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for region_id, counter, n_reads in ex.map(get_region_hp_counts, files, chunksize=200):
            out[region_id] = {"counter": counter, "n_reads": n_reads}
    return out


def compute_ng(counter: dict[str, int]) -> int:
    return sum(1 for k, v in counter.items() if k in NG_VALID and v > 0)


def bucket_hp(counter: dict[str, int]) -> dict[str, int]:
    out = {k: 0 for k in HP_FAMILY_KEYS}
    for k, v in counter.items():
        if k in out:
            out[k] += v
        else:
            out["other"] += v
    return out


def extract_caller_af(vcf_path: Path) -> dict[tuple[str, int], float]:
    out = {}
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t[%AF]\n", str(vcf_path)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return _extract_caller_af_fallback(vcf_path)
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        try:
            out[(parts[0], int(parts[1]))] = float(parts[2])
        except ValueError:
            continue
    return out


def _extract_caller_af_fallback(vcf_path: Path) -> dict[tuple[str, int], float]:
    out = {}
    with subprocess.Popen(["zcat", str(vcf_path)], stdout=subprocess.PIPE, text=True) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            fmt = cols[8].split(":")
            sample = cols[9].split(":")
            if "AF" in fmt:
                idx = fmt.index("AF")
                if idx < len(sample):
                    try:
                        out[(cols[0], int(cols[1]))] = float(sample[idx])
                    except ValueError:
                        pass
    return out


def load_master_paired_full() -> dict[tuple[str, int], dict]:
    out = {}
    with gzip.open(MASTER, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("sample") != "HCC1395" or row.get("mode") != "paired_full":
                continue
            try:
                key = (row["Chr"], int(row["Pos"]))
            except (KeyError, ValueError):
                continue
            out[key] = {
                "LOH_Bed_Overlap": row.get("LOH_Bed_Overlap", ""),
                "LOH_Bed_Annotation": row.get("LOH_Bed_Annotation", ""),
                "LOH_Subtype": row.get("LOH_Subtype", ""),
                "Coverage_Multiple": row.get("Coverage_Multiple", ""),
                "Coverage_Category": row.get("Coverage_Category", ""),
                "Diploid_Coverage_Used": row.get("Diploid_Coverage_Used", ""),
                "AF_master": row.get("AF", ""),
                "loh_side": row.get("loh_side", ""),
                "NumReads_master": row.get("NumReads", ""),
                "HPFineNGroups_master": row.get("HPFineNGroups", ""),
            }
    return out


def parse_region_id(rid: str) -> tuple[str, int]:
    parts = rid.split("_")
    return parts[0], (int(parts[1]) + int(parts[2])) // 2


def main() -> int:
    INTERMEDIATE.mkdir(parents=True, exist_ok=True)

    print("[Step 1 four-way] Stage 1: extracting per-region HP counts from 16 ISM runs ...", flush=True)
    data: dict[str, dict[str, dict[str, dict]]] = {b: {f: {} for f in FLAGS} for b in BAMS}
    for bam in BAMS:
        for flag, label in RUNS:
            run_dir = PHASEC / f"{bam}_{flag}_{label}"
            if not run_dir.exists():
                print(f"  [SKIP] {bam}_{flag}_{label}: missing")
                continue
            agg = aggregate_run(run_dir, n_workers=16)
            data[bam].setdefault(flag, {})[label] = agg
            print(f"  {bam}_{flag}_{label}: {len(agg)} regions, "
                  f"total reads={sum(d['n_reads'] for d in agg.values())}", flush=True)

    print("[Step 1 four-way] Stage 2: extracting caller_af ...", flush=True)
    af_tp = extract_caller_af(VCF_TP)
    af_fp = extract_caller_af(VCF_FP)
    print(f"  TP VCF: {len(af_tp)} | FP VCF: {len(af_fp)}", flush=True)

    af_by_region: dict[tuple[str, str], float] = {}
    all_regions: set[tuple[str, str]] = set()
    for bam in BAMS:
        for flag in FLAGS:
            for label in LABELS:
                for rid in data[bam].get(flag, {}).get(label, {}):
                    all_regions.add((label.upper(), rid))
    for label_upper, rid in all_regions:
        chrom, pos = parse_region_id(rid)
        af = (af_tp if label_upper == "TP" else af_fp).get((chrom, pos))
        af_by_region[(label_upper, rid)] = af if af is not None else float("nan")

    print("[Step 1 four-way] Stage 3: loading master.tsv covariates ...", flush=True)
    master = load_master_paired_full()
    print(f"  master.tsv paired_full rows: {len(master)}", flush=True)

    print("[Step 1 four-way] Stage 4: assembling wide-format master_four_way TSV ...", flush=True)
    regions_by_label: dict[str, set[str]] = {"TP": set(), "FP": set()}
    for bam in BAMS:
        for flag in FLAGS:
            for label in LABELS:
                regions_by_label[label.upper()].update(
                    data[bam].get(flag, {}).get(label, {}).keys()
                )

    feature_keys = HP_FAMILY_KEYS + ["NG", "n_reads"]
    wide_cols = ["region_id", "chr", "pos", "label"]
    for bam in BAMS:
        for flag in FLAGS:
            for k in feature_keys:
                wide_cols.append(f"{bam}_{flag}_{k}")
    wide_cols += [
        "caller_af", "master_join_ok",
        "LOH_Bed_Overlap", "LOH_Bed_Annotation", "LOH_Subtype",
        "Coverage_Multiple", "Coverage_Category", "Diploid_Coverage_Used",
        "loh_side", "AF_master", "HPFineNGroups_master", "NumReads_master",
    ]

    master_path = STEP1 / "step1_master_four_way.tsv"
    n_joined = 0
    with master_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=wide_cols, delimiter="\t")
        w.writeheader()
        for label in ["TP", "FP"]:
            for rid in sorted(regions_by_label[label]):
                chrom, pos = parse_region_id(rid)
                row: dict[str, object] = {"region_id": rid, "chr": chrom, "pos": pos, "label": label}
                for bam in BAMS:
                    for flag in FLAGS:
                        rec = data[bam].get(flag, {}).get(label.lower(), {}).get(rid)
                        if rec is None:
                            for k in feature_keys:
                                row[f"{bam}_{flag}_{k}"] = "NA"
                        else:
                            bucket = bucket_hp(rec["counter"])
                            for k in HP_FAMILY_KEYS:
                                row[f"{bam}_{flag}_{k}"] = bucket[k]
                            row[f"{bam}_{flag}_NG"] = compute_ng(rec["counter"])
                            row[f"{bam}_{flag}_n_reads"] = rec["n_reads"]
                af = af_by_region.get((label, rid))
                row["caller_af"] = f"{af:.6f}" if (af is not None and af == af) else "NA"

                m = master.get((chrom, pos))
                if m is None:
                    row["master_join_ok"] = 0
                    for c in ["LOH_Bed_Overlap", "LOH_Bed_Annotation", "LOH_Subtype",
                              "Coverage_Multiple", "Coverage_Category",
                              "Diploid_Coverage_Used", "loh_side",
                              "AF_master", "HPFineNGroups_master", "NumReads_master"]:
                        row[c] = "NA"
                else:
                    row["master_join_ok"] = 1
                    n_joined += 1
                    for c in ["LOH_Bed_Overlap", "LOH_Bed_Annotation", "LOH_Subtype",
                              "Coverage_Multiple", "Coverage_Category",
                              "Diploid_Coverage_Used", "loh_side",
                              "AF_master", "HPFineNGroups_master", "NumReads_master"]:
                        row[c] = m.get(c, "")
                w.writerow(row)

    total_rows = sum(len(s) for s in regions_by_label.values())
    print(f"  wrote {master_path} ({total_rows} rows, master-joined={n_joined} = {n_joined/total_rows:.1%})", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
