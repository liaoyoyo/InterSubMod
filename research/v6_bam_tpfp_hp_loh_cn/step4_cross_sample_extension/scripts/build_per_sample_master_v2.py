#!/usr/bin/env python3
"""Step 4 stage 1 v2 — Build per-sample V6-only master TSV using pre-generated path lists.

Reads from intermediate/file_lists/{sample}_{flag}_{label}.txt instead of running find.
This avoids the slow Python subprocess find overhead.

Output (per sample):
  step4_cross_sample_extension/per_sample_master/{sample}_v6_master.tsv
"""
from __future__ import annotations

import csv
import gzip
import os
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
PHASED = REPO / "research/paired_priority_bug_audit/phaseD_v6_5sample"
STEP4 = REPO / "research/v6_bam_tpfp_hp_loh_cn/step4_cross_sample_extension"
PER_SAMPLE = STEP4 / "per_sample_master"
FILE_LISTS = STEP4 / "intermediate" / "file_lists"
MASTER_TSV = REPO / "research/tpfp_loh_af_kde_discrimination/data/master.tsv.gz"

SAMPLES = ["H1437", "H2009", "HCC1954", "HCC1937"]
FLAGS = ["off", "on"]
LABELS = ["tp", "fp"]

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


def aggregate_from_list(file_list: Path, n_workers: int = 32) -> dict[str, dict]:
    if not file_list.exists():
        print(f"    [WARN] missing list: {file_list}", flush=True)
        return {}
    with file_list.open() as fh:
        files = [line.strip() for line in fh if line.strip()]
    print(f"    {len(files):,} files in list; parsing with {n_workers} workers ...", flush=True)
    out: dict[str, dict] = {}
    if not files:
        return out
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for region_id, counter, n_reads in ex.map(
            get_region_hp_counts, files, chunksize=500
        ):
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


def _vcf_path(sample: str, label: str) -> Path:
    return Path(
        f"/big8_disk/liaoyoyo2001/Somatic_methylation_analysis_output/{sample}/VCF/"
        f"filtered_snv_{label}.vcf.gz"
    )


def extract_caller_af(vcf_path: Path) -> dict[tuple[str, int], float]:
    out: dict[tuple[str, int], float] = {}
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
    out: dict[tuple[str, int], float] = {}
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


def load_master_for_sample(sample: str) -> dict[tuple[str, int], dict]:
    out: dict[tuple[str, int], dict] = {}
    with gzip.open(MASTER_TSV, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("sample") != sample:
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
                "kde_status": row.get("kde_status", ""),
                "AF_master": row.get("AF", ""),
                "loh_side": row.get("loh_side", ""),
                "NumReads_master": row.get("NumReads", ""),
                "HPFineNGroups_master": row.get("HPFineNGroups", ""),
            }
    return out


def parse_region_id(rid: str) -> tuple[str, int]:
    parts = rid.split("_")
    chrom = parts[0]
    start = int(parts[1])
    end = int(parts[2])
    pos = (start + end) // 2
    return chrom, pos


def build_sample(sample: str) -> dict:
    t0 = time.time()
    print(f"\n[Step 4 build v2] {sample} START", flush=True)
    data: dict[str, dict[str, dict]] = {f: {} for f in FLAGS}
    for flag in FLAGS:
        for label in LABELS:
            file_list = FILE_LISTS / f"{sample}_{flag}_{label}.txt"
            print(f"  {sample}/{flag}_{label}: reading {file_list.name}", flush=True)
            agg = aggregate_from_list(file_list, n_workers=32)
            data[flag][label] = agg
            print(
                f"    {sample}/{flag}_{label}: {len(agg):,} regions, "
                f"reads={sum(d['n_reads'] for d in agg.values()):,}",
                flush=True,
            )

    af_by_region: dict[tuple[str, str], float] = {}
    for label in LABELS:
        vcf = _vcf_path(sample, label)
        if not vcf.exists():
            print(f"  [WARN] VCF missing for {sample}/{label}: {vcf}", flush=True)
            continue
        af_map = extract_caller_af(vcf)
        for flag in FLAGS:
            for rid in data[flag].get(label, {}):
                chrom, pos = parse_region_id(rid)
                af_by_region[(label.upper(), rid)] = af_map.get(
                    (chrom, pos), float("nan")
                )
        print(f"  {sample}/{label}: caller_af records={len(af_map):,}", flush=True)

    master = load_master_for_sample(sample)
    print(f"  {sample}: master.tsv paired_full rows={len(master):,}", flush=True)

    regions_by_label: dict[str, set[str]] = {"TP": set(), "FP": set()}
    for flag in FLAGS:
        for label in LABELS:
            regions_by_label[label.upper()].update(
                data[flag].get(label, {}).keys()
            )

    feature_keys = HP_FAMILY_KEYS + ["NG", "n_reads"]
    wide_cols = ["region_id", "chr", "pos", "label"]
    for flag in FLAGS:
        for k in feature_keys:
            wide_cols.append(f"V6_{flag}_{k}")
    wide_cols += [
        "caller_af",
        "master_join_ok",
        "LOH_Bed_Overlap",
        "LOH_Bed_Annotation",
        "LOH_Subtype",
        "Coverage_Multiple",
        "Coverage_Category",
        "Diploid_Coverage_Used",
        "kde_status",
        "loh_side",
        "AF_master",
        "HPFineNGroups_master",
        "NumReads_master",
    ]

    PER_SAMPLE.mkdir(parents=True, exist_ok=True)
    out_path = PER_SAMPLE / f"{sample}_v6_master.tsv"
    n_total = 0
    n_joined = 0
    n_stale_diploid75 = 0
    n_af_missing = 0
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=wide_cols, delimiter="\t")
        w.writeheader()
        for label in ["TP", "FP"]:
            for rid in sorted(regions_by_label[label]):
                chrom, pos = parse_region_id(rid)
                row: dict[str, object] = {
                    "region_id": rid,
                    "chr": chrom,
                    "pos": pos,
                    "label": label,
                }
                for flag in FLAGS:
                    rec = data[flag].get(label.lower(), {}).get(rid)
                    if rec is None:
                        for k in feature_keys:
                            row[f"V6_{flag}_{k}"] = "NA"
                    else:
                        bucket = bucket_hp(rec["counter"])
                        for k in HP_FAMILY_KEYS:
                            row[f"V6_{flag}_{k}"] = bucket[k]
                        row[f"V6_{flag}_NG"] = compute_ng(rec["counter"])
                        row[f"V6_{flag}_n_reads"] = rec["n_reads"]
                af = af_by_region.get((label, rid))
                if af is None or af != af:
                    row["caller_af"] = "NA"
                    n_af_missing += 1
                else:
                    row["caller_af"] = f"{af:.6f}"

                m = master.get((chrom, pos))
                if m is None:
                    row["master_join_ok"] = 0
                    for c in [
                        "LOH_Bed_Overlap",
                        "LOH_Bed_Annotation",
                        "LOH_Subtype",
                        "Coverage_Multiple",
                        "Coverage_Category",
                        "Diploid_Coverage_Used",
                        "kde_status",
                        "loh_side",
                        "AF_master",
                        "HPFineNGroups_master",
                        "NumReads_master",
                    ]:
                        row[c] = "NA"
                else:
                    row["master_join_ok"] = 1
                    n_joined += 1
                    for c in [
                        "LOH_Bed_Overlap",
                        "LOH_Bed_Annotation",
                        "LOH_Subtype",
                        "Coverage_Multiple",
                        "Coverage_Category",
                        "Diploid_Coverage_Used",
                        "kde_status",
                        "loh_side",
                        "AF_master",
                        "HPFineNGroups_master",
                        "NumReads_master",
                    ]:
                        row[c] = m.get(c, "")
                    try:
                        if abs(float(m["Diploid_Coverage_Used"]) - 75.0) < 1e-6:
                            n_stale_diploid75 += 1
                    except (ValueError, KeyError, TypeError):
                        pass
                w.writerow(row)
                n_total += 1

    summary = {
        "sample": sample,
        "n_rows": n_total,
        "n_TP": len(regions_by_label["TP"]),
        "n_FP": len(regions_by_label["FP"]),
        "master_joined_n": n_joined,
        "master_join_rate": n_joined / max(n_total, 1),
        "diploid75_stale_n": n_stale_diploid75,
        "af_missing_n": n_af_missing,
        "output": str(out_path),
        "elapsed_sec": round(time.time() - t0, 1),
    }
    print(
        f"  wrote {out_path.name}: rows={n_total:,}, joined={n_joined:,} "
        f"({summary['master_join_rate']:.1%}), AF_missing={n_af_missing:,}, "
        f"stale_diploid75={n_stale_diploid75}, elapsed={summary['elapsed_sec']}s",
        flush=True,
    )
    return summary


def main() -> int:
    PER_SAMPLE.mkdir(parents=True, exist_ok=True)
    (STEP4 / "intermediate").mkdir(parents=True, exist_ok=True)
    summaries = []
    only = os.environ.get("ONLY_SAMPLE", "").strip()
    samples = [only] if only else SAMPLES
    for s in samples:
        # Skip if already done
        out_path = PER_SAMPLE / f"{s}_v6_master.tsv"
        if out_path.exists():
            print(f"[Step 4 build v2] {s} already done ({out_path.stat().st_size:,} bytes), skipping")
            continue
        try:
            summaries.append(build_sample(s))
        except Exception as e:
            print(f"[ERROR] {s} failed: {e}", flush=True)
            raise

    out = STEP4 / "intermediate" / "step4_per_sample_build_summary.tsv"
    keys = [
        "sample",
        "n_rows",
        "n_TP",
        "n_FP",
        "master_joined_n",
        "master_join_rate",
        "diploid75_stale_n",
        "af_missing_n",
        "elapsed_sec",
        "output",
    ]
    existing = []
    if out.exists():
        with out.open() as fh:
            r = csv.DictReader(fh, delimiter="\t")
            existing = [row for row in r]
    keep = {row["sample"]: row for row in existing}
    for s in summaries:
        keep[s["sample"]] = {k: s.get(k, "") for k in keys}
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        for s_name, row in sorted(keep.items()):
            w.writerow(row)
    print(f"\n[Step 4 build v2] wrote summary → {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
