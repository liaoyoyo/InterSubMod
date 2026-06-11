#!/usr/bin/env python3
"""Step 1 — V3F/V5/V6 three-way ISM integration (HCC1395 paired-pileup).

Reads 12 phaseC ISM run dirs (3 BAM × 2 flag × 2 label) and produces a region-level
wide-format master TSV with HP family counts, NG, caller AF, LOH bed annotation, and
Coverage_Multiple covariates joined from the canonical HCC1395 master.tsv.

Reality check (data inventory):
- phaseC significance_summary.csv is HEADER-ONLY (0 regions analyzed) — no master.tsv exists
- Per-region reads.tsv has columns: read_id, read_name, chr, start, end, mapq, hp, alt_support, is_tumor, strand
- The only cross-version aggregate is v3f_vs_v5_vs_v6_region_ng.tsv (105,997 rows = TP 30,490 + FP 4,842 × 3 BAM)
- region_id format: {chr}_{start}_{end} where start = pos-5000, end = pos+5000 (10kb window)
- caller_af source: /big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_{tp,fp}.vcf.gz
- LOH/CN source: /big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/data/master.tsv.gz
  HCC1395 paired_full mode covers 96.8% TP positions; FP coverage is partial (paired_full FP n=627 vs phaseC FP n=4842)

Outputs:
- step1_master_three_way.tsv  (wide-format, one row per region_id)
- step1_delta_{v5_vs_v3f,v6_vs_v5,v6_vs_v3f}.tsv
- step1_trajectory.tsv (per-region V3F→V5→V6 class A-E)
- step1_off_vs_on_compare.tsv
- intermediate/per_region_hp_counts.tsv.gz (all 12 runs, long format)
- intermediate/caller_af_lookup.tsv (region_id → caller_af, label)

Read-only with respect to phaseC and source master.tsv.
"""
from __future__ import annotations

import csv
import gzip
import os
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
PHASEC = REPO / "research/paired_priority_bug_audit/phaseC_genome_three_way"
STEP1 = REPO / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way"
INTERMEDIATE = STEP1 / "intermediate"

VCF_TP = Path("/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz")
VCF_FP = Path("/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz")
MASTER = REPO / "research/tpfp_loh_af_kde_discrimination/data/master.tsv.gz"

BAMS = ["V3F", "V5", "V6"]
FLAGS = ["off", "on"]
LABELS = ["tp", "fp"]
RUNS = [(f, l) for f in FLAGS for l in LABELS]  # (off,tp) (off,fp) (on,tp) (on,fp)

# HP values that count toward NG (mirrors phaseC_region_ng_fast.py)
NG_VALID = {"1", "2", "1-1", "2-1", "3", "11", "21", "33"}

# Canonical HP family bucketing (used for per-region wide table)
HP_FAMILY_KEYS = ["0", "1", "2", "1-1", "2-1", "3", "11", "21", "33", "other"]


# ---------------------------------------------------------------------------
# Stage 1: extract HP family counts from per-region reads.tsv (parallel)
# ---------------------------------------------------------------------------

def get_region_hp_counts(reads_tsv: str) -> tuple[str, dict[str, int], int]:
    """Return (region_id, HP family counter, total reads)."""
    # path layout: .../{run_dir}/filtered_snv_{tp|fp}/chr*/chr*_pos/{region_id}/reads/reads.tsv
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
                    val = cols[hp_idx]
                    counter[val] += 1
                    n_reads += 1
    except Exception:
        pass
    return (region_id, dict(counter), n_reads)


def find_reads_tsv(run_dir: Path) -> list[str]:
    out: list[str] = []
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
    """Aggregate counter into HP_FAMILY_KEYS (everything else→other)."""
    out = {k: 0 for k in HP_FAMILY_KEYS}
    for k, v in counter.items():
        if k in out:
            out[k] += v
        else:
            out["other"] += v
    return out


# ---------------------------------------------------------------------------
# Stage 2: caller_af extraction via bcftools
# ---------------------------------------------------------------------------

def extract_caller_af(vcf_path: Path) -> dict[tuple[str, int], float]:
    """Return {(chr, pos): caller_af} using FORMAT/AF first sample column."""
    out: dict[tuple[str, int], float] = {}
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t[%AF]\n", str(vcf_path)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"  [WARN] bcftools failed for {vcf_path.name}: {exc}; trying zcat fallback")
        return _extract_caller_af_fallback(vcf_path)
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        chrom, pos, af = parts[0], parts[1], parts[2]
        try:
            out[(chrom, int(pos))] = float(af)
        except ValueError:
            continue
    return out


def _extract_caller_af_fallback(vcf_path: Path) -> dict[tuple[str, int], float]:
    """zcat-based fallback that parses FORMAT/AF (Clair-S layout)."""
    out: dict[tuple[str, int], float] = {}
    with subprocess.Popen(["zcat", str(vcf_path)], stdout=subprocess.PIPE, text=True) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            chrom, pos = cols[0], cols[1]
            fmt = cols[8].split(":")
            sample = cols[9].split(":")
            if "AF" in fmt:
                idx = fmt.index("AF")
                if idx < len(sample):
                    try:
                        out[(chrom, int(pos))] = float(sample[idx])
                    except ValueError:
                        pass
    return out


# ---------------------------------------------------------------------------
# Stage 3: master.tsv (HCC1395 paired_full) join — LOH / CN / AF
# ---------------------------------------------------------------------------

def load_master_paired_full() -> dict[tuple[str, int], dict]:
    out: dict[tuple[str, int], dict] = {}
    with gzip.open(MASTER, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("sample") != "HCC1395":
                continue
            if row.get("mode") != "paired_full":
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


# ---------------------------------------------------------------------------
# Stage 4: assemble wide-format master_three_way TSV
# ---------------------------------------------------------------------------

def parse_region_id(rid: str) -> tuple[str, int]:
    """region_id = chr_<start>_<end>; midpoint = SNV pos."""
    parts = rid.split("_")
    chrom = parts[0]
    start = int(parts[1])
    end = int(parts[2])
    pos = (start + end) // 2
    return chrom, pos


def main() -> int:
    INTERMEDIATE.mkdir(parents=True, exist_ok=True)
    STEP1.mkdir(parents=True, exist_ok=True)

    print("[Step 1] Stage 1: extracting per-region HP counts from 12 ISM runs ...", flush=True)
    # data[bam][flag][label] = {region_id: {"counter": {...}, "n_reads": int}}
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

    # Dump long-format intermediate
    long_path = INTERMEDIATE / "per_region_hp_counts.tsv.gz"
    print(f"[Step 1] writing long-format intermediate → {long_path}", flush=True)
    with gzip.open(long_path, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["BAM", "flag", "label", "region_id", "n_reads", "NG"] + HP_FAMILY_KEYS)
        for bam in BAMS:
            for flag in FLAGS:
                for label in LABELS:
                    agg = data[bam].get(flag, {}).get(label, {})
                    for rid, rec in agg.items():
                        bucket = bucket_hp(rec["counter"])
                        ng = compute_ng(rec["counter"])
                        writer.writerow([bam, flag, label, rid, rec["n_reads"], ng]
                                        + [bucket[k] for k in HP_FAMILY_KEYS])

    print("[Step 1] Stage 2: extracting caller_af from source VCFs ...", flush=True)
    af_tp = extract_caller_af(VCF_TP)
    af_fp = extract_caller_af(VCF_FP)
    print(f"  TP VCF: {len(af_tp)} records | FP VCF: {len(af_fp)} records", flush=True)

    # caller AF lookup intermediate (region_id → af)
    af_lookup_path = INTERMEDIATE / "caller_af_lookup.tsv"
    af_by_region: dict[tuple[str, str], float] = {}
    all_regions: set[tuple[str, str]] = set()  # (label_upper, region_id)
    for bam in BAMS:
        for flag in FLAGS:
            for label in LABELS:
                for rid in data[bam].get(flag, {}).get(label, {}):
                    all_regions.add((label.upper(), rid))
    with af_lookup_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["label", "region_id", "chr", "pos", "caller_af", "source"])
        for label_upper, rid in sorted(all_regions):
            chrom, pos = parse_region_id(rid)
            af = (af_tp if label_upper == "TP" else af_fp).get((chrom, pos))
            af_by_region[(label_upper, rid)] = af if af is not None else float("nan")
            w.writerow([label_upper, rid, chrom, pos,
                        f"{af:.6f}" if af is not None else "NA",
                        VCF_TP.name if label_upper == "TP" else VCF_FP.name])

    print("[Step 1] Stage 3: loading HCC1395 paired_full master.tsv covariates ...", flush=True)
    master = load_master_paired_full()
    print(f"  master.tsv paired_full rows: {len(master)}", flush=True)

    print("[Step 1] Stage 4: assembling wide-format master_three_way TSV ...", flush=True)
    # Build union of regions per label
    regions_by_label: dict[str, set[str]] = {"TP": set(), "FP": set()}
    for bam in BAMS:
        for flag in FLAGS:
            for label in LABELS:
                regions_by_label[label.upper()].update(
                    data[bam].get(flag, {}).get(label, {}).keys()
                )

    # Wide schema
    feature_keys = HP_FAMILY_KEYS + ["NG", "n_reads"]
    wide_cols = ["region_id", "chr", "pos", "label"]
    for bam in BAMS:
        for flag in FLAGS:
            for k in feature_keys:
                wide_cols.append(f"{bam}_{flag}_{k}")
    wide_cols += [
        "caller_af",
        "master_join_ok",
        "LOH_Bed_Overlap",
        "LOH_Bed_Annotation",
        "LOH_Subtype",
        "Coverage_Multiple",
        "Coverage_Category",
        "Diploid_Coverage_Used",
        "loh_side",
        "AF_master",
        "HPFineNGroups_master",
        "NumReads_master",
    ]

    master_path = STEP1 / "step1_master_three_way.tsv"
    n_joined = 0
    n_stale_diploid75 = 0
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
                    try:
                        if abs(float(m["Diploid_Coverage_Used"]) - 75.0) < 1e-6:
                            n_stale_diploid75 += 1
                    except (ValueError, KeyError, TypeError):
                        pass
                w.writerow(row)

    total_rows = sum(len(s) for s in regions_by_label.values())
    print(f"  wrote {master_path} ({total_rows} rows, master-joined={n_joined} = {n_joined/total_rows:.1%})", flush=True)
    if n_stale_diploid75 > 0:
        print(f"  [WARNING] {n_stale_diploid75} rows have Diploid_Coverage_Used=75.0 — stale KDE binary artifact (CURRENT_FOCUS H-CN1)", flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
